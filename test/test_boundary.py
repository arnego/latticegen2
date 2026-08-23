"""Boundary trimming and the symmetric interface rule (docs/algorithm.md §7).

A cap is a hole one junction punches for its neighbour to fill. Punching it is
only sound if the neighbour punched the matching one, and the two sides are
decided by two *independent* booleans — so the decision cannot be made locally
from classification. These tests pin the rule that replaced that: a cap becomes
an interface only when both sides present material there and the two regions
agree, and anything else stays as exterior surface rather than becoming a hole.
"""

import os

import numpy as np
import pytest
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut
from OCP.TopTools import TopTools_ListOfShape

from OCP.BRepTools import BRepTools

from latticegen2 import occ, weld
from latticegen2.boundary import (
    BoundaryPiece,
    CAP_AREA_REL_TOL,
    PINHOLE_WIRE_TOL,
    _open_shell,
    _piece_from,
    finalize_pieces,
    fuse_disagreeing_pairs,
    resolve_interfaces,
    trim_boundary,
    trim_junction,
)
from latticegen2.connect import build_components
from latticegen2.errors import ProcessingError
from latticegen2.junction import build_template, is_cap_plane_face
from latticegen2.parallel import WorkerPool
from latticegen2.lattice import (
    HALF_STRUTS,
    OPPOSITE_HALF,
    half_strut_offset,
    lattice_params,
    neighbor_step,
    node as lattice_node,
    nodes as lattice_nodes,
    profile_vertices,
)

CC, T = 10.0, 1.5


@pytest.fixture(scope="module")
def lp():
    return lattice_params(CC, T)


def across(node, h):
    step = neighbor_step(h)
    return (node[0] + int(step[0]), node[1] + int(step[1]), node[2] + int(step[2]))


def cap_face(lp, node, h, scale=1.0):
    """A face lying in node ``h``'s cap plane, ``scale`` of the full quad's size.

    ``scale=1`` reproduces the whole ``t x t`` cap; a smaller value stands in for
    a cap a trim cut down, which has the same plane and a smaller area.
    """
    k, _ = HALF_STRUTS[h]
    centre = lattice_node(lp, node) + half_strut_offset(lp, h)
    corners = profile_vertices(lp, centre, k)
    return occ.polygon_face(centre + scale * (corners - centre))


def piece(lp, node, caps, volume=1.0, scale=1.0):
    """A boundary piece presenting ``caps``, with no other faces."""
    p = BoundaryPiece(node=node, volume=volume)
    for h in caps:
        p.cap_faces[(node, h)] = [cap_face(lp, node, h, scale)]
    return p


def nodes(*rows):
    return np.array(rows, dtype=np.int64) if rows else np.empty((0, 3), dtype=np.int64)


# --- the trim tags caps rather than dropping them ---------------------------


def test_trim_tags_every_cap_plane_face_and_keeps_it(lp):
    """No face is discarded in the worker — the master decides, later."""
    tpl = build_template(lp)
    box = occ.prism(
        occ.polygon_face(np.array([[-9.0, -9.0, -9.0], [9.0, -9.0, -9.0],
                                   [9.0, 9.0, -9.0], [-9.0, 9.0, -9.0]])),
        np.array([0.0, 0.0, 18.0]),
    )
    trim = trim_junction(lp, tpl, np.zeros(3), box)
    trimmed = trim.pieces
    assert len(trimmed) == 1
    assert trim.n_pinholes_removed == 0, (
        "a junction wholly inside a box has no grazing pinholes"
    )
    faces, tags, volume = trimmed[0]
    assert len(faces) == len(tags)
    assert volume == pytest.approx(tpl.volume, rel=1e-9), "the box contains the junction whole"
    # A junction wholly inside the body comes back with all six caps present.
    assert sorted(h for h in tags if h >= 0) == list(range(6))
    assert all(
        (is_cap_plane_face(lp, f, np.zeros(3)) is None) == (h < 0)
        for f, h in zip(faces, tags)
    )


# --- resolve_interfaces -----------------------------------------------------


def test_two_whole_caps_facing_each_other_become_an_interface(lp):
    a = piece(lp, (0, 0, 0), [0])
    b = piece(lp, across((0, 0, 0), 0), [OPPOSITE_HALF[0]])
    got = resolve_interfaces(lp, nodes(), [a, b])
    assert got.n_pairs == 1
    assert got.unpaired == [] and got.mismatched == []


def test_the_interface_set_always_holds_both_sides(lp):
    """Either side can test membership directly, so the two can never disagree."""
    a = piece(lp, (0, 0, 0), [0, 1])
    b = piece(lp, across((0, 0, 0), 0), [OPPOSITE_HALF[0]])
    got = resolve_interfaces(lp, nodes(across((0, 0, 0), 1)), [a, b])
    assert got.n_pairs == 2
    for node, h in got.interfaces:
        assert (across(node, h), OPPOSITE_HALF[h]) in got.interfaces


def test_an_interior_node_presents_all_six_caps_whole(lp):
    """docs/algorithm.md §5.3(b): an INTERIOR node's cap planes are strictly inside."""
    got = resolve_interfaces(lp, nodes((0, 0, 0), (1, 0, 0)), [])
    assert got.n_pairs == 1
    assert ((0, 0, 0), 0) in got.interfaces


def test_a_cap_facing_a_junction_that_produced_none_is_not_an_interface(lp):
    """The failure the whole rework exists for.

    The neighbour produced geometry, just nothing at this cap — the two booleans
    disagreed about a face they share. Opening the interface anyway punches a
    hole with nothing behind it.
    """
    a = piece(lp, (0, 0, 0), [0])
    b = piece(lp, across((0, 0, 0), 0), [1])  # geometry, but not at the shared cap
    got = resolve_interfaces(lp, nodes(), [a, b])
    assert got.n_pairs == 0
    assert ((0, 0, 0), 0) in got.unpaired


def test_a_cap_facing_empty_space_is_not_reported_as_an_anomaly(lp):
    """An exterior cap is ordinary. Only a *disagreement* is worth reporting."""
    a = piece(lp, (0, 0, 0), [0])
    got = resolve_interfaces(lp, nodes(), [a])
    assert got.n_pairs == 0 and got.unpaired == []


def test_caps_that_disagree_in_area_are_not_stitched_across(lp):
    a = piece(lp, (0, 0, 0), [0], scale=1.0)
    b = piece(lp, across((0, 0, 0), 0), [OPPOSITE_HALF[0]], scale=0.5)
    got = resolve_interfaces(lp, nodes(), [a, b])
    assert got.n_pairs == 0
    assert len(got.mismatched) == 1
    node, h, area_a, area_b = got.mismatched[0]
    assert area_a == pytest.approx(T * T, rel=1e-9)
    assert area_b == pytest.approx(T * T / 4.0, rel=1e-9)


def test_equally_partial_caps_are_a_perfectly_good_interface(lp):
    """Partiality is not the problem; asymmetry is. A cut cap still stitches."""
    a = piece(lp, (0, 0, 0), [0], scale=0.5)
    b = piece(lp, across((0, 0, 0), 0), [OPPOSITE_HALF[0]], scale=0.5)
    got = resolve_interfaces(lp, nodes(), [a, b])
    assert got.n_pairs == 1
    assert got.mismatched == []


def test_fragments_of_one_cap_split_across_pieces_are_summed(lp):
    """A trim can leave part of one cap on each of two pieces of a junction.

    What faces the neighbour is the union, so the comparison has to be against
    the total. Judging each fragment on its own would reject a sound interface.
    """
    node, other = (0, 0, 0), across((0, 0, 0), 0)
    half = 1.0 / np.sqrt(2.0)  # area t^2 / 2 each
    halves = [piece(lp, node, [0], scale=half), piece(lp, node, [0], scale=half)]
    whole = piece(lp, other, [OPPOSITE_HALF[0]])
    got = resolve_interfaces(lp, nodes(), halves + [whole])
    assert got.n_pairs == 1
    assert got.mismatched == []
    finalize_pieces(halves + [whole], got.interfaces)
    assert all(p.caps == frozenset({(node, 0)}) for p in halves)


# --- finalize_pieces --------------------------------------------------------


def test_an_interface_cap_is_dropped_and_any_other_cap_is_kept(lp):
    a = piece(lp, (0, 0, 0), [0, 1])
    b = piece(lp, across((0, 0, 0), 0), [OPPOSITE_HALF[0]])
    got = resolve_interfaces(lp, nodes(), [a, b])
    finalize_pieces([a, b], got.interfaces)

    assert a.caps == frozenset({((0, 0, 0), 0)})
    # Cap 0 is stitched across so it goes; cap 1 faces nothing, so it stays and
    # closes the piece there rather than leaving a hole.
    assert len(a.faces) == 1
    assert is_cap_plane_face(lp, a.faces[0], lattice_node(lp, a.node)) == 1
    assert set(a.cap_faces) == {((0, 0, 0), 0)}, \
        "only interface caps are held back, for their rings"


def test_a_piece_whose_caps_all_stay_gives_up_nothing(lp):
    a = piece(lp, (0, 0, 0), [0])
    got = resolve_interfaces(lp, nodes(), [a])
    finalize_pieces([a], got.interfaces)
    assert a.caps == frozenset()
    assert len(a.faces) == 1
    assert a.cap_faces == {}


# --- the two halves agree ---------------------------------------------------


def test_resolved_interfaces_never_trip_the_connectivity_invariant(lp):
    """What resolve_interfaces produces, build_components must always accept.

    The mix here is deliberately awkward — an interior neighbour, a disagreeing
    cap, a cap facing nothing — because those are exactly the cases the old
    one-sided rule got wrong.
    """
    interior = nodes((0, 0, 0))
    a = piece(lp, across((0, 0, 0), 0), [OPPOSITE_HALF[0], 1])
    b = piece(lp, across(across((0, 0, 0), 0), 1), [OPPOSITE_HALF[1]], scale=0.25)
    got = resolve_interfaces(lp, interior, [a, b])
    finalize_pieces([a, b], got.interfaces)

    comps = build_components(interior, 10.0, [a, b], got.interfaces)
    assert len(comps.volumes) == 2, "the disagreeing cap leaves b on its own"
    assert len(got.mismatched) == 1


def test_an_unresolved_cap_still_fails_loudly(lp):
    """The invariant is a real check, not a formality that resolution defeats."""
    a = piece(lp, (0, 0, 0), [0])
    a.caps = frozenset({((0, 0, 0), 0)})
    with pytest.raises(ProcessingError, match="hole with no matching hole"):
        build_components(nodes(), 0.0, [a], interfaces={((0, 0, 0), 0)})


@pytest.mark.parametrize(
    "excess, stitched", [(0.5 * CAP_AREA_REL_TOL, True), (2.0 * CAP_AREA_REL_TOL, False)]
)
def test_the_agreement_bar_sits_where_the_constant_says(lp, excess, stitched):
    """Quadrature noise passes; a real difference in the region does not."""
    a = piece(lp, (0, 0, 0), [0])
    b = piece(lp, across((0, 0, 0), 0), [OPPOSITE_HALF[0]], scale=np.sqrt(1.0 - excess))
    got = resolve_interfaces(lp, nodes(), [a, b])
    assert (got.n_pairs == 1) is stitched
    assert (got.mismatched == []) is stitched


# --- fuse_disagreeing_pairs (docs/algorithm.md §7.1) -------------------------


def _notched_junction(lp, tpl, node_pos, oh):
    """A whole junction with a real corner clipped off cap ``oh``.

    A genuine boolean cut, not the synthetic uniform ``scale`` the other tests
    in this file use — the regression this repair exists for needs a real
    *mismatched partial* cap, not two proportionally-shrunk copies of the same
    region, which would still agree in shape.
    """
    k, _ = HALF_STRUTS[oh]
    cap_c = node_pos + half_strut_offset(lp, oh)
    u, v, e = lp.u[k], lp.v[k], lp.e[k]
    r = lp.r
    lo, hi, depth = 0.2 * r, 3.0 * r, 1.5 * r
    corners = np.array([
        cap_c - depth * e + lo * u + lo * v,
        cap_c - depth * e + hi * u + lo * v,
        cap_c - depth * e + hi * u + hi * v,
        cap_c - depth * e + lo * u + hi * v,
    ])
    notch = occ.prism(occ.polygon_face(corners), 2 * depth * e)

    instance = tpl.solid.Moved(occ.translation(node_pos))
    cut = BRepAlgoAPI_Cut()
    args = TopTools_ListOfShape()
    args.Append(instance)
    tools = TopTools_ListOfShape()
    tools.Append(notch)
    cut.SetArguments(args)
    cut.SetTools(tools)
    cut.Build()
    assert cut.IsDone()
    solids = occ.solids(cut.Shape())
    assert len(solids) == 1, "the notch must not sever the junction"
    return solids[0]


def _tag(lp, node_pos, solid):
    faces = occ.faces(solid)
    tags = [is_cap_plane_face(lp, f, node_pos) for f in faces]
    return faces, [-1 if h is None else h for h in tags]


@pytest.fixture(scope="module")
def disagreeing_pieces(lp):
    """Two real trimmed junctions whose shared cap genuinely disagrees.

    ``a`` is a whole, untouched junction (full ``t x t`` cap 0). ``b`` has a
    real corner cut off its matching cap 3, so the two sides present different
    *partial* regions of the same nominal cap — not the "whole vs sliver"
    scale mismatch the synthetic tests use, but the kind two independent
    ``BRepAlgoAPI_Common`` calls actually produce when they disagree
    (docs/algorithm.md §7.1).
    """
    tpl = build_template(lp)
    a, b = (0, 0, 0), across((0, 0, 0), 0)
    pos_a, pos_b = lattice_nodes(lp, np.array([a, b], dtype=np.int64))

    a_solid = tpl.solid.Moved(occ.translation(pos_a))
    faces_a, tags_a = _tag(lp, pos_a, a_solid)
    piece_a = _piece_from(a, faces_a, tags_a, occ.volume(a_solid))

    b_solid = _notched_junction(lp, tpl, pos_b, OPPOSITE_HALF[0])
    faces_b, tags_b = _tag(lp, pos_b, b_solid)
    piece_b = _piece_from(b, faces_b, tags_b, occ.volume(b_solid))

    return a, b, piece_a, piece_b


def test_a_notched_cap_is_reported_as_a_genuine_mismatch(lp, disagreeing_pieces):
    a, b, piece_a, piece_b = disagreeing_pieces
    got = resolve_interfaces(lp, nodes(), [piece_a, piece_b])
    assert got.mismatched == [(a, 0, pytest.approx(T * T), pytest.approx(T * T, abs=0.3))]
    assert got.mismatched[0][2] != pytest.approx(got.mismatched[0][3])
    assert got.n_pairs == 0


def test_fuse_disagreeing_pairs_produces_one_agreed_solid(lp, disagreeing_pieces):
    """The regression this repair exists for: two pieces whose shared cap
    disagrees must assemble into a closed orientable shell, edge-use tally
    clean, rather than declining and leaving non-manifold overlap plus a hole.
    """
    a, b, piece_a, piece_b = disagreeing_pieces
    mismatched = resolve_interfaces(lp, nodes(), [piece_a, piece_b]).mismatched

    fused, n_groups = fuse_disagreeing_pairs(lp, [piece_a, piece_b], mismatched)
    assert n_groups == 1
    assert len(fused) == 1
    merged = fused[0]
    assert merged.volume == pytest.approx(piece_a.volume + piece_b.volume, rel=1e-9)

    # Resolved a second time, the disagreement is simply gone: the two nodes
    # are already one solid, so neither presents that cap as a boundary face
    # any more, and there is nothing left to decline.
    iface2 = resolve_interfaces(lp, nodes(), fused)
    assert iface2.mismatched == []
    assert (a, 0) not in iface2.interfaces
    assert (b, OPPOSITE_HALF[0]) not in iface2.interfaces

    # Every other cap of both nodes survived the fuse correctly attributed —
    # the proximity disambiguation in `_owning_cap` earning its keep, since a
    # plain one-axis `is_cap_plane_face` test alone cannot tell node a's caps
    # from node b's for any half-strut whose axis is orthogonal to the one
    # separating them.
    remaining = set(merged.cap_faces)
    assert remaining == {(a, h) for h in range(6)} | {
        (b, h) for h in range(6) if h != OPPOSITE_HALF[0]
    }

    finalize_pieces(fused, iface2.interfaces)
    shell = _open_shell(merged.faces)
    open_edges, misoriented = weld.shell_defects(shell)[:2]
    assert (open_edges, misoriented) == (0, 0)
    shell.Closed(True)
    assert occ.is_valid(occ.make_solid(shell))


def test_fuse_disagreeing_pairs_is_a_no_op_without_mismatches(lp):
    a = piece(lp, (0, 0, 0), [0])
    b = piece(lp, across((0, 0, 0), 0), [OPPOSITE_HALF[0]])
    fused, n_groups = fuse_disagreeing_pairs(lp, [a, b], [])
    assert n_groups == 0
    assert fused == [a, b]


# --- zero-area pinhole wires (docs/algorithm.md §7) -------------------------
#
# A near-tangential trim leaves a face carrying an extra inner wire of one
# few-micron edge whose endpoints do not meet — a pinhole bounding no area. The
# edge is used once, so the every-edge-twice proof rejects the shell for it.
#
# These tests use the **real** part and the real node, not a synthetic
# reproduction. That is deliberate: an earlier attempt at this defect built a
# synthetic near-tangential graze which reproduced the symptom's scale to four
# significant figures (a non-degenerate ~3e-06 mm edge) and was still the wrong
# object — the synthetic edges were ordinary two-owner small edges, which OCCT's
# ShapeFix removes happily, while the real ones are one-owner pinhole wires it
# cannot touch at all. See tools/prototypes/RESULTS.md G10.

PINHOLE_NODE = (591, -46, -70)
"""The junction on `TD_HX_rehearsal_test` at `cc=5, t=1` whose trim leaves the two
pinholes docs/specification.md §10 documents, at 3.171690e-06 and 5.808982e-06
mm."""


@pytest.fixture(scope="module")
def pinhole_case():
    """`(lp, template, body, node position)` for the real failing junction."""
    lp5 = lattice_params(5.0, 1.0)
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    body = occ.read_step(os.path.join(here, "test", "TD_HX_rehearsal_test.step"))
    pos = lattice_node(lp5, np.array(PINHOLE_NODE))
    return lp5, build_template(lp5), body, pos


def raw_trim(tpl, node_pos, body):
    """The boolean alone, with no repair — what `trim_junction` starts from."""
    algo = BRepAlgoAPI_Common()
    args = TopTools_ListOfShape()
    args.Append(tpl.solid.Moved(occ.translation(node_pos)))
    tools = TopTools_ListOfShape()
    tools.Append(body)
    algo.SetArguments(args)
    algo.SetTools(tools)
    algo.Build()
    assert algo.IsDone()
    solids = occ.solids(algo.Shape())
    assert len(solids) == 1
    return solids[0]


def test_the_raw_trim_really_does_leave_an_unpaired_pinhole(pinhole_case):
    """Guard on the guard: if the boolean stops producing them, the repair test
    below would pass while testing nothing."""
    _, tpl, body, pos = pinhole_case
    raw = raw_trim(tpl, pos, body)
    # The defect is visible on the piece alone, before any sewing — an edge used
    # by exactly one face inside a solid OCCT itself calls valid.
    assert occ.is_valid(raw)
    assert weld.shell_defects(_open_shell(occ.faces(raw)))[:2] == (2, 0)


def test_trim_junction_removes_the_pinholes_and_closes_the_piece(pinhole_case):
    lp5, tpl, body, pos = pinhole_case
    raw = raw_trim(tpl, pos, body)

    trim = trim_junction(lp5, tpl, pos, body)
    pieces, n_pinholes = trim.pieces, trim.n_pinholes_removed

    assert n_pinholes == 2
    assert len(pieces) == 1
    faces, tags, _ = pieces[0]
    assert len(faces) == len(tags), "tags stay aligned with the repaired face list"
    assert weld.shell_defects(_open_shell(faces))[:2] == (0, 0)
    # Exactly the two pinhole edges are gone; no face is lost with them.
    assert occ.count_subshapes(_open_shell(faces))[1] == occ.count_subshapes(raw)[1] - 2
    assert len(faces) == len(occ.faces(raw))


def test_removing_a_pinhole_moves_no_geometry(pinhole_case):
    """A wire bounding no area cannot change the area of the face carrying it,
    so this is checked exactly rather than within a comfortable margin."""
    lp5, tpl, body, pos = pinhole_case
    raw = raw_trim(tpl, pos, body)
    faces, tags, volume = trim_junction(lp5, tpl, pos, body).pieces[0]

    assert sum(occ.area(f) for f in faces) == sum(occ.area(f) for f in occ.faces(raw))
    assert occ.is_valid(occ.make_solid(_closed(_open_shell(faces))))

    # And the caps resolve_interfaces compares are untouched, so a repaired
    # piece still agrees with an unrepaired neighbour across a shared cap.
    before = _cap_area_map(lp5, pos, occ.faces(raw))
    after = {h: a for h, a in _cap_area_map(lp5, pos, faces).items()}
    assert set(before) == set(after)
    for h in before:
        assert after[h] == pytest.approx(before[h], abs=0.0, rel=0.0)


def _closed(shell):
    shell.Closed(True)
    return shell


def _cap_area_map(lp, node_pos, faces):
    per: dict[int, float] = {}
    for f in faces:
        h = is_cap_plane_face(lp, f, node_pos)
        if h is not None:
            per[h] = per.get(h, 0.0) + occ.area(f)
    return per


REFUSED_CASE = (12.0, 2.5, (262, -27, -37))
"""`(cc, t, node)` — the junction a relative-volume bar of 1e-9 refused.

The only pinhole-bearing junction of `TD_HX_rehearsal_test` at these
parameters, out of 2,907 boundary junctions. Its wire is 1.010641e-06 mm, on a
6.01 mm^2 planar face of a 77.4 mm^3 piece, and removing it moves the volume
OCCT reports by 1.235e-09 relative — over a bar the same repair on the same
part at `cc=5, t=1` clears by six orders. Both parameter sets are inside the
CLI's ranges, so this was a valid input being refused (docs/algorithm.md §11).
"""


@pytest.fixture(scope="module")
def refused_case():
    """`(lp, template, body, node position)` for the junction above."""
    cc, t, ijk = REFUSED_CASE
    lp12 = lattice_params(cc, t)
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    body = occ.read_step(os.path.join(here, "test", "TD_HX_rehearsal_test.step"))
    return lp12, build_template(lp12), body, lattice_node(lp12, np.array(ijk))


def test_a_pinhole_that_shifts_the_reported_volume_is_still_repaired(refused_case):
    """The regression: this junction must trim, not raise.

    What moves is the volume OCCT *reports*, not the piece. The unrepaired piece
    has a free boundary — that is what a pinhole is — and
    `BRepGProp::VolumeProperties` requires a shape "exempt of any free boundary",
    so the pre-repair figure carries a bias that disappears with the wire. Gate
    G19 measures that bias reproduced within 1 % by adding a synthetic open wire
    to a clean face, at a magnitude independent of the wire's length over three
    decades and set by which face carries it.
    """
    lp12, tpl, body, pos = refused_case
    raw = raw_trim(tpl, pos, body)

    result = trim_junction(lp12, tpl, pos, body)

    assert result.n_pinholes_removed == 1
    assert len(result.pieces) == 1
    faces, _, volume = result.pieces[0]
    assert weld.shell_defects(_open_shell(faces))[:2] == (0, 0)
    # Area is the quantity that can be trusted here, and it is bit-identical.
    assert sum(occ.area(f) for f in faces) == sum(occ.area(f) for f in occ.faces(raw))
    # The volume genuinely differs, well past what any relative bar could absorb
    # without also absorbing a real change. It is the defect's footprint on the
    # integral, and the repaired figure is the one the pipeline goes on to use.
    assert abs(volume - occ.volume(raw)) / occ.volume(raw) > 1e-9


def test_the_repair_is_proved_structurally_rather_than_by_an_integral(refused_case):
    """`only_inner_wires_dropped` is what replaced the volume bar."""
    lp12, tpl, body, pos = refused_case
    raw = raw_trim(tpl, pos, body)
    cleaned, n_removed = occ.remove_pinhole_wires(raw, PINHOLE_WIRE_TOL)
    assert n_removed == 1

    assert occ.only_inner_wires_dropped(raw, cleaned, n_removed) is None
    # Only the repaired face is rebuilt; every other face comes back as the
    # same object, which is what makes this a walk over objects, not geometry.
    same = sum(1 for a, b in zip(occ.faces(raw), occ.faces(cleaned)) if a.IsSame(b))
    assert same == len(occ.faces(raw)) - 1
    # And it does not simply return None: an unaccounted-for wire is rejected.
    assert occ.only_inner_wires_dropped(raw, cleaned, n_removed + 1) is not None
    assert occ.only_inner_wires_dropped(cleaned, raw, n_removed) is not None


def test_the_structural_check_rejects_a_face_that_lost_a_real_wire(lp):
    """A wire that bounded something changes an area, and is caught by name."""
    tpl = build_template(lp)
    faces = occ.faces(tpl.solid)
    # Dropping a whole face is the coarsest possible violation of "only inner
    # wires went"; the check must not need geometry to see it.
    reason = occ.only_inner_wires_dropped(tpl.solid, _open_shell(faces[:-1]), 0)
    assert reason is not None and "face count" in reason


def test_a_paired_wire_is_never_removed_however_small_the_bar_says(lp):
    """The safety property that makes the length threshold not load-bearing.

    Run at a tolerance larger than *every* edge in the template, so length alone
    would condemn all of them. Nothing is removed, because every edge of a
    closed solid is used twice and the repair only ever touches already-unpaired
    edges on a non-outer wire.
    """
    tpl = build_template(lp)
    huge = 10.0 * lp.a
    fixed, n_removed = occ.remove_pinhole_wires(tpl.solid, huge)
    assert n_removed == 0
    assert fixed is tpl.solid, "a no-op must not rebuild the shape either"


def test_a_clean_trim_has_no_pinholes(lp):
    """Negative control: a junction wholly inside a box is left alone."""
    tpl = build_template(lp)
    box = occ.prism(
        occ.polygon_face(np.array([[-9.0, -9.0, -9.0], [9.0, -9.0, -9.0],
                                   [9.0, 9.0, -9.0], [-9.0, 9.0, -9.0]])),
        np.array([0.0, 0.0, 18.0]),
    )
    assert trim_junction(lp, tpl, np.zeros(3), box).n_pinholes_removed == 0


def test_a_repair_that_moved_geometry_is_a_named_failure(lp, monkeypatch):
    """Silence here is the failure the guard exists to prevent: a face quietly
    reshaped in the worker surfaces much later as a wrong volume."""
    tpl = build_template(lp)
    box = occ.prism(
        occ.polygon_face(np.array([[-9.0, -9.0, -9.0], [9.0, -9.0, -9.0],
                                   [9.0, 9.0, -9.0], [-9.0, 9.0, -9.0]])),
        np.array([0.0, 0.0, 18.0]),
    )
    # Stand in for a repair that claims to have removed a pinhole but returns a
    # different solid. No synthetic pinhole reproduces the real pathology (G10),
    # so the guard is exercised directly rather than through a fake defect.
    smaller = occ.prism(
        occ.polygon_face(np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                                   [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])),
        np.array([0.0, 0.0, 1.0]),
    )
    monkeypatch.setattr(occ, "remove_pinhole_wires", lambda shape, tol: (smaller, 1))
    with pytest.raises(ProcessingError, match="changed its surface area"):
        trim_junction(lp, tpl, np.zeros(3), box)


# --- the trim reports progress while it runs, not after ---------------------


def _trim_a_grid(lp, tmp_path, workers, progress):
    """Trim nine junctions against a containing box, with or without a pool."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    tpl = build_template(lp)
    box = occ.prism(
        occ.polygon_face(np.array([[-9.0, -9.0, -9.0], [9.0, -9.0, -9.0],
                                   [9.0, 9.0, -9.0], [-9.0, 9.0, -9.0]])),
        np.array([0.0, 0.0, 18.0]),
    )
    body_path = str(tmp_path / "input.brep")
    BRepTools.Write_s(box, body_path)
    grid = nodes(*[(i, j, 0) for i in range(-1, 2) for j in range(-1, 2)])
    args = (lp, tpl, grid, box, body_path, str(tmp_path))
    if workers > 1:
        with WorkerPool(workers) as pool:
            result = trim_boundary(*args, workers=workers, progress=progress, pool=pool)
    else:
        result = trim_boundary(*args, workers=1, progress=progress)
    return result, len(grid)


def test_progress_counts_junctions_on_the_sequential_path(lp, tmp_path):
    ticks = []
    _result, total = _trim_a_grid(
        lp, tmp_path / "serial", 1, lambda done, tot: ticks.append((done, tot))
    )
    assert ticks == [(i, total) for i in range(1, total + 1)]


def test_progress_counts_junctions_on_the_parallel_path_too(lp, tmp_path):
    ticks = []
    result, total = _trim_a_grid(
        lp, tmp_path / "parallel", 2, lambda done, tot: ticks.append((done, tot))
    )
    assert ticks[-1] == (total, total)
    assert [d for d, _ in ticks] == sorted(d for d, _ in ticks)
    assert all(1 <= d <= total and tot == total for d, tot in ticks)
    assert len(result.pieces) == total, "every junction is inside the box"


def test_the_parallel_path_reports_from_the_dispatch_loop(lp, tmp_path, monkeypatch):
    """The regression this hook exists for.

    Progress used to be reported from ``_consume``, which runs only *after*
    ``pool.run`` has completely returned — so a 13-minute stage emitted its whole
    counter in the final millisecond, which is no use to a bar. Pinned by the
    mechanism rather than by timing: the ticks must come from ``on_result``,
    which :func:`test_parallel.test_on_result_fires_once_per_job_with_a_rising_count`
    separately pins as firing inside the dispatch loop.
    """
    seen = {}
    real_run = WorkerPool.run

    def spy(self, worker_fn, jobs, **kwargs):
        seen["on_result"] = kwargs.get("on_result")
        return real_run(self, worker_fn, jobs, **kwargs)

    monkeypatch.setattr(WorkerPool, "run", spy)
    _trim_a_grid(lp, tmp_path / "spy", 2, lambda done, tot: None)
    assert seen["on_result"] is not None


def test_the_pieces_are_the_same_whether_or_not_anyone_is_watching(lp, tmp_path):
    """Progress is an observer: it may not change what the stage produces."""
    watched, _ = _trim_a_grid(lp, tmp_path / "watched", 2, lambda d, t: None)
    unwatched, _ = _trim_a_grid(lp, tmp_path / "unwatched", 2, None)
    assert [p.node for p in watched.pieces] == [p.node for p in unwatched.pieces]
    assert [p.volume for p in watched.pieces] == [p.volume for p in unwatched.pieces]
    assert watched.n_empty == unwatched.n_empty
