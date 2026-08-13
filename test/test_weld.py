"""Assembly by adoption (docs/algorithm.md §8).

The boundary layer is sewn to itself, and the interior is then built *onto* the
topology that came out — adopting its vertices and edges rather than making its
own and reconciling afterwards. These tests pin the two properties that took the
longest to get right, both of which fail silently when they are wrong:

* the assembled surface must be closed **and orientable**, because a shell whose
  faces are joined back-to-front still gives every edge exactly two users;
* an adopted edge's own direction is not ours to choose, so the index has to
  observe it rather than assume it.
"""

import numpy as np
import pytest

from OCP.BRep import BRep_Builder
from OCP.TopAbs import TopAbs_ShapeEnum
from OCP.TopoDS import TopoDS, TopoDS_Shell

from latticegen2 import occ, weld
from latticegen2.boundary import BoundaryPiece, trim_junction
from latticegen2.connect import lattice_interfaces
from latticegen2.interior import build_interior_shell, extract_template_mesh
from latticegen2.junction import build_template, is_cap_plane_face
from latticegen2.lattice import OPPOSITE_HALF, lattice_params, neighbor_step, nodes

CC, T = 10.0, 1.5


@pytest.fixture(scope="module")
def template():
    lp = lattice_params(CC, T)
    tpl = build_template(lp)
    return lp, tpl, extract_template_mesh(lp, tpl)


def big_box():
    return occ.prism(
        occ.polygon_face(np.array([[-50.0, -50.0, -50.0], [50.0, -50.0, -50.0],
                                   [50.0, 50.0, -50.0], [-50.0, 50.0, -50.0]])),
        np.array([0.0, 0.0, 100.0]),
    )


# --- ring matching ----------------------------------------------------------


def test_a_trimmed_cap_matches_the_ring_predicted_from_lattice_maths(template):
    """The interior's side of an interface is known before any geometry exists.

    That is what lets an interface that will not join be rejected while both
    sides still have their cap face to fall back on.
    """
    lp, tpl, tmesh = template
    node = (0, 0, 0)
    faces, tags, _ = trim_junction(lp, tpl, np.zeros(3), big_box())[0]
    cap = next(f for f, t in zip(faces, tags) if t == 0)
    predicted = weld.template_cap_ring(lp, tmesh, node, 0)
    assert weld.match_rings(predicted, weld.ring_of_face(cap)) is not None


def test_rings_of_different_size_never_match(template):
    lp, _, tmesh = template
    four = weld.template_cap_ring(lp, tmesh, (0, 0, 0), 0)
    three = weld.Ring(edges=four.edges[:3], verts=four.verts[:3], points=four.points[:3])
    assert weld.match_rings(four, three) is None


def test_matching_is_indifferent_to_which_way_each_edge_runs(template):
    """The two sides are free to have parametrised the same segment either way."""
    lp, _, tmesh = template
    ring = weld.template_cap_ring(lp, tmesh, (0, 0, 0), 0)
    flipped = weld.Ring(
        edges=list(ring.edges),
        verts=list(ring.verts),
        points=[p[::-1] for p in ring.points],
    )
    assert weld.match_rings(ring, flipped) is not None


# --- the watertightness proof ----------------------------------------------


def shell_of(faces):
    shell = TopoDS_Shell()
    builder = BRep_Builder()
    builder.MakeShell(shell)
    for f in faces:
        builder.Add(shell, f)
    return shell


def test_a_closed_shell_has_no_defects(template):
    lp, tpl, _ = template
    faces, _, _ = trim_junction(lp, tpl, np.zeros(3), big_box())[0]
    assert weld.shell_defects(shell_of(faces))[:2] == (0, 0)


def test_a_face_joined_back_to_front_is_caught(template):
    """Counting uses alone would pass this: every edge still has two faces.

    Measured on the 80 mm ball while this was being built — 0 open edges, 954
    edges traversed the same way twice, and a volume of 29,111 mm³ against a true
    51,393 mm³, with nothing else looking wrong.
    """
    lp, tpl, _ = template
    faces, _, _ = trim_junction(lp, tpl, np.zeros(3), big_box())[0]
    flipped = [faces[0].Reversed()] + list(faces[1:])
    open_edges, misoriented, _, _ = weld.shell_defects(shell_of(flipped))
    assert open_edges == 0, "the flip does not open the shell — that is the point"
    assert misoriented > 0


def test_a_missing_face_is_caught(template):
    lp, tpl, _ = template
    faces, _, _ = trim_junction(lp, tpl, np.zeros(3), big_box())[0]
    open_edges, _, samples, by_use = weld.shell_defects(shell_of(faces[1:]))
    assert set(by_use) == {1}, "a missing face leaves edges on one face, not three"
    assert open_edges > 0
    assert samples, "an open edge must be located, not just counted"


# --- adoption ---------------------------------------------------------------


def test_the_interior_joins_a_trimmed_junction_by_adopting_its_topology(template):
    """The whole mechanism, at two nodes.

    Node (0,0,0) is instanced by the index; node (1,0,0) comes out of a boolean,
    exactly as a real boundary piece does. The index adopts the boundary's ring,
    so the two share four edges rather than owning two coincident copies — and
    the assembled solid's volume is the exact sum, which is only true if nothing
    was double-counted or lost.
    """
    lp, tpl, tmesh = template
    a, b = (0, 0, 0), (1, 0, 0)
    h = 0
    interfaces = {(a, h), (b, OPPOSITE_HALF[h])}

    pos_b = nodes(lp, np.array([b], dtype=np.int64))[0]
    faces, tags, vol_b = trim_junction(lp, tpl, pos_b, big_box())[0]
    assert vol_b == pytest.approx(tpl.volume, rel=1e-12), "the box contains it whole"
    boundary_faces = [f for f, t in zip(faces, tags) if t != OPPOSITE_HALF[h]]

    rings = weld.interface_rings(lp, tmesh, {0: boundary_faces}, {(a, h): 0})
    assert set(rings) == {(a, h)}

    built = build_interior_shell(
        lp, tpl, tmesh, np.array([a], dtype=np.int64), interfaces,
        groups={a: 0}, adopted=rings,
    )
    shells = weld.assemble(built.shells, {0: boundary_faces})
    solids, _ = weld.close_shells(shells)

    assert len(solids) == 1
    assert occ.is_valid(solids[0])
    assert occ.volume(solids[0]) == pytest.approx(2 * tpl.volume, rel=1e-9)


def test_adoption_shares_the_edges_rather_than_duplicating_them(template):
    """Four edges are shared, so the assembly has four fewer than the two parts."""
    lp, tpl, tmesh = template
    a, b = (0, 0, 0), (1, 0, 0)
    interfaces = {(a, 0), (b, OPPOSITE_HALF[0])}
    pos_b = nodes(lp, np.array([b], dtype=np.int64))[0]
    faces, tags, _ = trim_junction(lp, tpl, pos_b, big_box())[0]
    boundary_faces = [f for f, t in zip(faces, tags) if t != OPPOSITE_HALF[0]]

    rings = weld.interface_rings(lp, tmesh, {0: boundary_faces}, {(a, 0): 0})
    built = build_interior_shell(
        lp, tpl, tmesh, np.array([a], dtype=np.int64), interfaces,
        groups={a: 0}, adopted=rings,
    )
    assembled = weld.assemble(built.shells, {0: boundary_faces})[0]

    separate = len(occ.faces(shell_of(boundary_faces)))
    assert len(occ.faces(assembled)) == separate + built.stats["interior_faces"]
    # Every edge of the assembly is used twice: nothing is coincident-but-distinct.
    assert weld.shell_defects(assembled)[:2] == (0, 0)


def test_a_ring_the_boundary_does_not_present_is_a_named_failure(template):
    """Silence here would become an unclosed shell much later, with no location."""
    lp, tpl, tmesh = template
    faces, tags, _ = trim_junction(lp, tpl, np.zeros(3), big_box())[0]
    with pytest.raises(Exception, match="no unique free edge"):
        weld.interface_rings(lp, tmesh, {0: faces}, {((7, 7, 7), 0): 0})


# --- tiling the boundary sew (docs/specification.md §10) --------------------


def _long_box():
    """A box wide enough along ``+e0`` to hold a chain of ~20 nodes untrimmed.

    ``big_box()`` is only 100 mm on a side, which a chain of even a dozen nodes
    outgrows at this module's ``cc=10`` — trimming would then genuinely cut into
    the last few junctions, which is not what these tests are isolating.
    """
    return occ.prism(
        occ.polygon_face(np.array([[-150.0, -150.0, -150.0], [150.0, -150.0, -150.0],
                                   [150.0, 150.0, -150.0], [-150.0, 150.0, -150.0]])),
        np.array([0.0, 0.0, 300.0]),
    )


def _line_pieces(lp, tpl, n: int) -> list[BoundaryPiece]:
    """``n`` trimmed junctions in a row along ``+e0``, interfaces already dropped.

    Each node is fully inside :func:`_long_box`, so trimming leaves the whole
    template; the shared cap between consecutive nodes (``h=0`` outgoing,
    ``h=3`` incoming — the same pair the adoption tests above use) is removed
    from both sides, exactly what :func:`latticegen2.boundary.finalize_pieces`
    leaves behind for a real interface. Sewing must rediscover that pairing from
    the resulting holes: every interior hole has a matching one to close against,
    and the two ends keep their outward cap intact, so the whole chain sews into
    one fully **closed** shell — a tube capped at both ends, not an open pipe.
    """
    positions = nodes(lp, np.array([(i, 0, 0) for i in range(n)], dtype=np.int64))
    pieces = []
    for i in range(n):
        faces, tags, vol = trim_junction(lp, tpl, positions[i], _long_box())[0]
        keep = [
            f for f, t in zip(faces, tags)
            if not (t == 3 and i > 0) and not (t == 0 and i < n - 1)
        ]
        pieces.append(BoundaryPiece(node=(i, 0, 0), volume=vol, faces=keep))
    return pieces


def test_tile_pieces_declines_to_split_a_small_or_compact_group(template):
    lp, tpl, _ = template
    pieces = _line_pieces(lp, tpl, 4)
    assert weld._tile_pieces(pieces, target=1, min_to_tile=100) is None, \
        "below the piece-count threshold, tiling must not engage at all"
    assert weld._tile_pieces(pieces, target=1000, min_to_tile=1) is None, \
        "a target far bigger than the group's own footprint is one tile"


def test_tile_pieces_splits_a_large_spread_group_and_keeps_every_piece(template):
    lp, tpl, _ = template
    pieces = _line_pieces(lp, tpl, 12)
    tiles = weld._tile_pieces(pieces, target=2, min_to_tile=1)
    assert tiles is not None
    assert len(tiles) > 1, "12 nodes spread over a target of 2 must split"
    assert sorted(p.node for tile in tiles for p in tile) == sorted(p.node for p in pieces), \
        "tiling only regroups pieces, it must never drop or duplicate one"


def test_tiled_and_untiled_sew_produce_the_same_watertight_result(template):
    """The point of tiling: identical geometry, reached by a different route.

    Forced into several tiles by a tiny target on one side, forced into a single
    untiled sew by a huge ``min_to_tile`` on the other, both starting from the
    same pieces built the same way.
    """
    lp, tpl, _ = template
    pieces = _line_pieces(lp, tpl, 10)
    groups = [0] * len(pieces)

    tiled, tiled_stats = weld.sew_boundary(pieces, groups, tile_target=2, min_to_tile=1)
    plain, plain_stats = weld.sew_boundary(pieces, groups, tile_target=10**9, min_to_tile=10**9)

    assert tiled_stats.tiles > 1 and tiled_stats.tiled_components == 1
    assert plain_stats.tiles == 0 and plain_stats.tiled_components == 0

    assert len(tiled[0]) == len(plain[0]), "sewing merges topology, never faces"
    for faces, stats in ((tiled[0], tiled_stats), (plain[0], plain_stats)):
        # A chain capped at both ends is fully closed: every interior hole has a
        # match, and the two end caps were never dropped in the first place.
        open_edges, misoriented, *_ = weld.shell_defects(shell_of(faces))
        assert (open_edges, misoriented) == (0, 0), f"stats={stats}"

    solid_tiled = occ.make_solid(shell_of(tiled[0]))
    solid_plain = occ.make_solid(shell_of(plain[0]))
    assert occ.volume(solid_tiled) == pytest.approx(10 * tpl.volume, rel=1e-9)
    assert occ.volume(solid_plain) == pytest.approx(occ.volume(solid_tiled), rel=1e-12)


def test_tiled_sew_across_worker_processes_matches_the_sequential_path(template, tmp_path):
    """The one code path that actually crosses a process boundary.

    Every other tiling test takes the ``workers=1``/``tmpdir=None`` branch of
    :func:`latticegen2.weld._sew_all_tiles`, same as :func:`_tile_pieces`'s own
    unit tests — none of them exercise the ``.brep`` round-trip through a real
    ``multiprocessing`` pool. This does, the same way
    :func:`latticegen2.boundary.trim_boundary`'s worker path is exercised for
    real: not by a unit test (spawning processes per-test is slow and, per this
    project's own findings, fragile under some launch methods), but at least
    once, here, so a regression in the IPC plumbing itself — a bad path, a
    mismatched job tuple, a worker exception swallowed instead of propagated —
    fails fast in CI rather than only showing up at rehearsal scale.
    """
    lp, tpl, _ = template
    pieces = _line_pieces(lp, tpl, 8)
    groups = [0] * len(pieces)

    parallel, stats = weld.sew_boundary(
        pieces, groups, workers=4, tmpdir=str(tmp_path), background=False,
        tile_target=2, min_to_tile=1,
    )
    assert stats.tiles > 1
    assert stats.max_worker_rss > 0, "a worker actually ran and reported its RSS"

    open_edges, misoriented, *_ = weld.shell_defects(shell_of(parallel[0]))
    assert (open_edges, misoriented) == (0, 0)
    solid = occ.make_solid(shell_of(parallel[0]))
    assert occ.volume(solid) == pytest.approx(8 * tpl.volume, rel=1e-9)
