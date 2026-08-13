"""Boundary trimming and the symmetric interface rule (docs/algorithm.md §7).

A cap is a hole one junction punches for its neighbour to fill. Punching it is
only sound if the neighbour punched the matching one, and the two sides are
decided by two *independent* booleans — so the decision cannot be made locally
from classification. These tests pin the rule that replaced that: a cap becomes
an interface only when both sides present material there and the two regions
agree, and anything else stays as exterior surface rather than becoming a hole.
"""

import numpy as np
import pytest
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCP.TopTools import TopTools_ListOfShape

from latticegen2 import occ, weld
from latticegen2.boundary import (
    BoundaryPiece,
    CAP_AREA_REL_TOL,
    _open_shell,
    _piece_from,
    finalize_pieces,
    fuse_disagreeing_pairs,
    resolve_interfaces,
    trim_junction,
)
from latticegen2.connect import build_components
from latticegen2.errors import ProcessingError
from latticegen2.junction import build_template, is_cap_plane_face
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
    trimmed = trim_junction(lp, tpl, np.zeros(3), box)
    assert len(trimmed) == 1
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


# --- fuse_disagreeing_pairs (specification.md §10 "Fuse junction pairs whose
# two booleans disagree", docs/algorithm.md §7.1) ----------------------------


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
