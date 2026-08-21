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
from latticegen2.errors import ProcessingError
from latticegen2.interior import build_interior_shell, extract_template_mesh
from latticegen2.junction import build_template, is_cap_plane_face
from latticegen2.lattice import OPPOSITE_HALF, lattice_params, neighbor_step, nodes
from latticegen2.parallel import WorkerPool

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
    faces, tags, _ = trim_junction(lp, tpl, np.zeros(3), big_box()).pieces[0]
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


def test_free_edges_does_not_count_a_degenerate_edge_as_a_hole():
    """A degenerate edge has no extent, and the one face that owns it uses it
    once by construction — so counting it as a free edge counts a hole that is
    not there. `shell_defects` has always skipped them; `free_edges` did not,
    and that difference is what made even a correct unsplit sew miss the
    round-2 check by exactly the 10 degenerate edges `TD_HX_rehearsal_test`
    leaves at cc=5, t=1 (docs/specification.md §10).

    A sphere is the cheapest real carrier: the kernel gives it two degenerate
    pole edges and a seam its single face uses twice, so a correct count finds
    no free edge at all.
    """
    from OCP.BRep import BRep_Tool
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeSphere
    from OCP.TopoDS import TopoDS

    face = occ.faces(BRepPrimAPI_MakeSphere(10.0).Shape())[0]
    degenerate = [
        e for e in occ._explore(face, TopAbs_ShapeEnum.TopAbs_EDGE)
        if BRep_Tool.Degenerated_s(TopoDS.Edge_s(e))
    ]
    assert degenerate, "the fixture must actually carry degenerate edges"
    assert weld.free_edges([face]) == []


def test_a_closed_shell_has_no_defects(template):
    lp, tpl, _ = template
    faces, _, _ = trim_junction(lp, tpl, np.zeros(3), big_box()).pieces[0]
    assert weld.shell_defects(shell_of(faces))[:2] == (0, 0)


def test_a_face_joined_back_to_front_is_caught(template):
    """Counting uses alone would pass this: every edge still has two faces.

    Measured on the 80 mm ball while this was being built — 0 open edges, 954
    edges traversed the same way twice, and a volume of 29,111 mm³ against a true
    51,393 mm³, with nothing else looking wrong.
    """
    lp, tpl, _ = template
    faces, _, _ = trim_junction(lp, tpl, np.zeros(3), big_box()).pieces[0]
    flipped = [faces[0].Reversed()] + list(faces[1:])
    open_edges, misoriented, _, _ = weld.shell_defects(shell_of(flipped))
    assert open_edges == 0, "the flip does not open the shell — that is the point"
    assert misoriented > 0


def test_a_missing_face_is_caught(template):
    lp, tpl, _ = template
    faces, _, _ = trim_junction(lp, tpl, np.zeros(3), big_box()).pieces[0]
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
    faces, tags, vol_b = trim_junction(lp, tpl, pos_b, big_box()).pieces[0]
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
    faces, tags, _ = trim_junction(lp, tpl, pos_b, big_box()).pieces[0]
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
    faces, tags, _ = trim_junction(lp, tpl, np.zeros(3), big_box()).pieces[0]
    with pytest.raises(Exception, match="no unique free edge"):
        weld.interface_rings(lp, tmesh, {0: faces}, {((7, 7, 7), 0): 0})


# --- tiling the boundary sew (docs/algorithm.md §8) -------------------------


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
        faces, tags, vol = trim_junction(lp, tpl, positions[i], _long_box()).pieces[0]
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

    # A closed chain adopts no interior interfaces of its own, so the correctly
    # sewn result has zero free edges either way — an empty ``want_rings`` states
    # that expectation and lets the round-2 verification run without spuriously
    # triggering the repair path on a sound split.
    tiled, tiled_stats = weld.sew_boundary(
        pieces, groups, tile_target=2, min_to_tile=1, want_rings={}
    )
    plain, plain_stats = weld.sew_boundary(
        pieces, groups, tile_target=10**9, min_to_tile=10**9, want_rings={}
    )

    assert tiled_stats.tiles > 1 and tiled_stats.tiled_components == 1
    assert plain_stats.tiles == 0 and plain_stats.tiled_components == 0
    assert tiled_stats.repaired_components == 0, "a correct split must not be redone"
    assert plain_stats.repaired_components == 0, "an untiled component is never checked"

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


# --- seam-only round 2 (G8, docs/algorithm.md §8) ---------------------------


def test_split_seam_interior_separates_tile_seams_from_already_closed_faces(template):
    """Only faces bearing a free edge of their *own tile's* sew end up in
    ``seam``; a tile with an inter-tile boundary has some (the pieces facing
    the neighbouring tile), and a bulk of already-closed faces goes straight
    to ``interior`` untouched.
    """
    lp, tpl, _ = template
    pieces = _line_pieces(lp, tpl, 12)
    tiles = weld._tile_pieces(pieces, target=3, min_to_tile=1)
    assert tiles is not None and len(tiles) > 1

    tile_results = [weld._sew_faces([p.faces for p in tile], weld.SEW_TOLERANCE) for tile in tiles]
    seam_lists, interior = weld._split_seam_interior(tile_results)

    total_in = sum(len(tr) for tr in tile_results)
    seam_ids = [id(f) for group in seam_lists for f in group]
    interior_ids = [id(f) for f in interior]
    # Every face is accounted for exactly once, split between the two groups.
    assert len(seam_ids) + len(interior_ids) == total_in
    assert set(seam_ids).isdisjoint(interior_ids)
    all_ids = {id(f) for tr in tile_results for f in tr}
    assert set(seam_ids) | set(interior_ids) == all_ids
    # A chain of several tiles genuinely has inter-tile seams to find.
    assert len(seam_ids) > 0
    assert len(interior_ids) > 0


def test_split_seam_interior_finds_nothing_to_split_on_a_fully_closed_tile(template):
    """A tile with zero free edges (already fully closed) puts everything in
    ``interior`` and returns no seam list for it at all — nothing for round 2
    to even consider."""
    lp, tpl, tmesh = template
    faces, _, _ = trim_junction(lp, tpl, np.zeros(3), big_box()).pieces[0]
    assert weld.shell_defects(shell_of(faces))[:2] == (0, 0)  # sanity: closed

    seam_lists, interior = weld._split_seam_interior([faces])
    assert seam_lists == []
    assert sorted(id(f) for f in interior) == sorted(id(f) for f in faces)


def test_seam_only_round_two_matches_a_full_round_two(template):
    """The identity check G8 (`tools/prototypes/RESULTS.md`) ran ad hoc, pinned
    as a permanent regression: sewing only the seam subset and carrying the
    rest through unsewn must reproduce a full round-2 sew exactly — same face
    count, still fully closed, same volume.
    """
    lp, tpl, _ = template
    pieces = _line_pieces(lp, tpl, 16)
    plan = {0: weld._tile_pieces(pieces, target=3, min_to_tile=1)}
    assert plan[0] is not None and len(plan[0]) > 1

    tile_results = {0: [weld._sew_faces([p.faces for p in tile], weld.SEW_TOLERANCE) for tile in plan[0]]}
    by_group = {0: pieces}

    baseline = weld._sew_faces(tile_results[0], weld.SEW_TOLERANCE)
    seam_only, _max_rss, _repaired = weld._sew_round_two(
        by_group, plan, tile_results, weld.SEW_TOLERANCE, workers=1, tmpdir=None,
        pool=None,
    )

    assert len(seam_only[0]) == len(baseline)
    assert weld.free_edges(baseline) == []
    assert weld.free_edges(seam_only[0]) == []
    baseline_volume = occ.volume(occ.make_solid(shell_of(baseline)))
    seam_only_volume = occ.volume(occ.make_solid(shell_of(seam_only[0])))
    assert seam_only_volume == pytest.approx(baseline_volume, rel=1e-9)
    assert seam_only_volume == pytest.approx(16 * tpl.volume, rel=1e-9)


def test_round_two_repairs_a_component_the_seam_split_got_wrong(template, monkeypatch):
    """The safety net a production regression showed the split needed.

    `_split_seam_interior`'s argument — that sewing the seam-only subset in
    isolation reproduces exactly what a full round 2 would have done — held at
    every prototype scale tried (G8), all of them on lightly trimmed junctions,
    but not on the real, heavily trimmed geometry of the `cc=5, t=1` production
    rehearsal (docs/specification.md §10): there, `BRepBuilderAPI_Sewing` could
    rebuild a "straddling" edge shared with a carried-through face onto a new
    `TopoDS_Edge`, leaving 118,760 open edges at `assemble` where 10 were
    expected. A setup this small cannot reproduce that OCCT behaviour directly
    (no prototype scale has) — so the split is monkeypatched to drop one seam
    face, standing in for "the split silently produced a wrong result", and the
    fix is verified by its symptom: a wrong free-edge count is caught against
    ``want_rings`` and that component is redone on the unsplit tile results,
    ending up fully closed regardless.
    """
    lp, tpl, _ = template
    pieces = _line_pieces(lp, tpl, 12)
    groups = [0] * len(pieces)

    real_split = weld._split_seam_interior
    called = []

    def broken_split(face_lists):
        seam, interior = real_split(face_lists)
        called.append(1)
        if len(called) == 1:
            # Drop one seam face from its tile's list -- an incomplete split,
            # standing in for the real defect's edge substitution. The dropped
            # face is neither sewn nor carried, so the seam-only result alone
            # would come back with extra free edges where it used to border.
            for i, group_faces in enumerate(seam):
                if group_faces:
                    seam = list(seam)
                    seam[i] = group_faces[1:]
                    break
        return seam, interior

    monkeypatch.setattr(weld, "_split_seam_interior", broken_split)

    # A closed chain (both ends capped) adopts no interior interfaces, so the
    # correct free-edge count for this component is zero.
    out, stats = weld.sew_boundary(
        pieces, groups, tile_target=3, min_to_tile=1, want_rings={}
    )

    assert stats.repaired_components == 1, "the broken split must be caught and redone"
    open_edges, misoriented, *_ = weld.shell_defects(shell_of(out[0]))
    assert (open_edges, misoriented) == (0, 0), "the repaired component must still close"
    solid = occ.make_solid(shell_of(out[0]))
    assert occ.volume(solid) == pytest.approx(12 * tpl.volume, rel=1e-9)


def test_a_hole_the_unsplit_sew_cannot_close_fails_in_stitch_not_in_assemble(template):
    """A free-edge count still wrong after the full unsplit sew is a hard failure.

    The unsplit sew is what a component was sewn with before
    :func:`weld._split_seam_interior` existed, so a count still wrong there is
    not the split's doing and no further re-sew can help. It is a hole in the
    boundary layer itself, and the interior cannot close it: it adopts the rings
    it was told about and nothing else. Left to carry on, the run spends a whole
    `instance` stage building an interior for a layer that can never close, then
    fails in `assemble` naming assembly rather than the sew — which is exactly
    what a v3.0.0 report of this showed, at 17 and 14 edges over two runs, in
    both cases precisely the excess `_sew_round_two` had already measured.

    The hole is made by deleting one ordinary face from one piece, which is the
    production symptom rather than a wrong expectation: the geometry really is
    short of a face and every route to sewing it produces the same free edges.
    """
    lp, tpl, _ = template
    pieces = _line_pieces(lp, tpl, 12)
    # Not an end piece and not a cap: an interior piece's lateral face, whose
    # removal leaves a hole with no partner anywhere in the chain.
    del pieces[5].faces[0]
    groups = [0] * len(pieces)

    with pytest.raises(ProcessingError) as excinfo:
        weld.sew_boundary(pieces, groups, tile_target=3, min_to_tile=1, want_rings={})

    message = str(excinfo.value)
    assert "free edge(s)" in message
    assert "unsplit sew" in message,         "the message must say the split is not the cause, or it misdirects the reader"
    assert "position" in message and "[[" in message,         "positions are the whole point: counts alone cannot locate a hole in a part"


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
        pieces, groups, workers=4, tmpdir=str(tmp_path),
        tile_target=2, min_to_tile=1,
    )
    assert stats.tiles > 1
    assert stats.max_worker_rss > 0, "a worker actually ran and reported its RSS"

    open_edges, misoriented, *_ = weld.shell_defects(shell_of(parallel[0]))
    assert (open_edges, misoriented) == (0, 0)
    solid = occ.make_solid(shell_of(parallel[0]))
    assert occ.volume(solid) == pytest.approx(8 * tpl.volume, rel=1e-9)


def test_tiled_sew_falls_back_to_sequential_when_the_shared_pool_is_inert(template, tmp_path):
    """``--cores 1`` on a part big enough to tile: the production shape of it.

    ``pipeline._run`` builds ``WorkerPool(args.workers)`` unconditionally, and
    ``WorkerPool(1)`` is inert by design — it creates no ``mp.Pool``, so
    ``.active`` is ``False`` and ``.run()`` refuses. So ``_sew_all_tiles`` is
    handed a pool object that exists and cannot run anything, which is a case
    neither ``pool is None`` nor ``pool.active`` alone describes.

    Testing only ``pool is None`` let that fall through to the transient-pool
    branch, which built a *second* inert pool and raised ``ProcessingError`` out
    of ``run()`` — exit 4 on input the CLI accepts (``--cores`` is 1-128 per
    specification.md §3), which docs/algorithm.md §11 rules out. It needed a
    component past ``MIN_PIECES_TO_TILE`` to reach, so no committed scenario
    does: ``dense-lattice`` tiles nothing, and every other test here passes
    ``pool=None``. ``min_to_tile=1`` buys that condition at eight pieces.
    """
    lp, tpl, _ = template
    pieces = _line_pieces(lp, tpl, 8)
    groups = [0] * len(pieces)

    with WorkerPool(1) as inert:
        assert not inert.active, "WorkerPool(1) is inert — the premise of this test"
        out, stats = weld.sew_boundary(
            pieces, groups, workers=1, tmpdir=str(tmp_path),
            tile_target=2, min_to_tile=1, pool=inert,
        )

    assert stats.tiles > 1, "the tiling path was actually taken"
    assert stats.max_worker_rss == 0, "no worker ran — this is the sequential route"

    # And the route still produces the shell, not merely an absence of crash.
    open_edges, misoriented, *_ = weld.shell_defects(shell_of(out[0]))
    assert (open_edges, misoriented) == (0, 0)
    solid = occ.make_solid(shell_of(out[0]))
    assert occ.volume(solid) == pytest.approx(8 * tpl.volume, rel=1e-9)


# --- vertex tolerances the sew leaves wrong (docs/algorithm.md §8) ----------
#
# Sewing can leave an edge whose vertex is recorded as sitting off the edge's
# own 3D curve, with that vertex's tolerance inflated to *exactly* the
# distance, so BRepCheck's test sits on the knife edge and rejects both faces
# sharing the edge. On the cc=5, t=1 rehearsal that is 17 edges and 34 faces.
#
# The two fixtures are the real faces, lifted out of that run's assembled
# solid: one on a cylinder (elliptical trim curve, deviation 2.474044e-05 mm)
# and one on a B-spline surface (deviation 3.316370e-04 mm). Both are kept
# because they discriminate between the candidate repairs — BRepLib's
# UpdateTolerances fixes the first and not the second, which is why the fix
# uses ShapeFix_Edge (tools/prototypes/RESULTS.md G11).

BAD_FACES = ("invalid-vertex-tolerance-ellipse.brep",
             "invalid-vertex-tolerance-bspline.brep")


def load_face(name):
    import os

    from OCP.TopoDS import TopoDS as _TopoDS

    from latticegen2.parallel import read_brep

    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return _TopoDS.Face_s(read_brep(os.path.join(here, "test", name)))


@pytest.mark.parametrize("name", BAD_FACES)
def test_the_fixture_really_is_invalid(name):
    """Guard on the guard: if these ever load as valid, the repair test below
    would pass while proving nothing."""
    assert not occ.is_valid(load_face(name))


@pytest.mark.parametrize("name", BAD_FACES)
def test_correcting_the_vertex_tolerance_makes_the_face_valid(name):
    face = load_face(name)
    before = occ.area(face)

    repaired, still_invalid = occ.fix_vertex_tolerances([face])

    assert (repaired, still_invalid) == (1, 0)
    assert occ.is_valid(face)
    # A tolerance is metadata, not geometry: the area must be bit-identical,
    # not merely close.
    assert occ.area(face) == before


def test_a_sound_face_is_left_alone(template):
    """Negative control — the repair only looks at faces the analyzer rejects."""
    lp, tpl, _ = template
    faces, _tags, _vol = trim_junction(lp, tpl, np.zeros(3), big_box()).pieces[0]
    assert all(occ.is_valid(f) for f in faces), "the trim itself is clean"
    assert occ.fix_vertex_tolerances(faces) == (0, 0)


def test_the_batch_scan_is_the_same_predicate_as_the_per_face_one(template):
    """`fix_vertex_tolerances` finds its work with one parallel analyzer per
    window rather than one per face (docs/algorithm.md §8, G22), and then
    confirms every candidate with the standalone check. What is pinned here is
    the whole two-stage result: it must contain exactly the faces the
    standalone predicate rejects — no more, no less.

    The two stages have never been observed to disagree — the confirmation is
    insurance on the case no corpus of *loose* fault fixtures can reach, where a
    compound analyzer shares one subshape's status between neighbouring faces.
    What this test pins is the contract, so that stays true. Pinned against the
    real committed faults as well as clean faces, because a scan that only ever
    agrees on sound geometry proves nothing (G10), and with a chunk small enough
    to put a window boundary between faces that share an edge.
    """
    from OCP.BRepCheck import BRepCheck_Analyzer

    lp, tpl, _ = template
    clean, _tags, _vol = trim_junction(lp, tpl, np.zeros(3), big_box()).pieces[0]
    corpus = list(clean) + [load_face(n) for n in BAD_FACES + SELF_INTERSECTING_FACES]

    per_face = {
        id(f) for f in corpus if not BRepCheck_Analyzer(f).IsValid()
    }
    assert len(per_face) == len(BAD_FACES) + len(SELF_INTERSECTING_FACES), (
        "the fixtures must load invalid, or this test proves nothing"
    )
    for chunk in (3, len(corpus)):
        assert {id(f) for f in occ.invalid_faces(corpus, chunk=chunk)} == per_face


def test_a_repair_that_moved_geometry_is_a_named_failure(monkeypatch):
    """The bound this repair is held to: it adjusts recorded tolerances and
    must move nothing. Silence here would let a real geometry change through
    on the one code path that is allowed to touch a proven-watertight shell."""
    from latticegen2.errors import ProcessingError

    face = load_face(BAD_FACES[0])
    areas = iter([1.0, 2.0])
    monkeypatch.setattr(occ, "area", lambda _shape: next(areas))
    with pytest.raises(ProcessingError, match="must move no geometry"):
        occ.fix_vertex_tolerances([face])


# --- falsely self-intersecting wires (docs/algorithm.md §8) -----------------
#
# What the repair above leaves behind: a face whose edges and vertices are all
# valid *standalone* and which passes every named face check, yet the analyzer
# rejects. The fault is BRepCheck_SelfIntersectingWire, reported for two edges
# adjacent in the wire — a tight-tolerance trim edge against the fat-tolerance
# B-spline the boolean fitted to the strut/input-surface intersection. It is
# not a real self-intersection: the two pcurves cross once, at the shared
# vertex, inside its tolerance. The vertex tolerance is simply recorded a
# little too tight to swallow the crossing.
#
# Both fixtures are the real faces from the cc=5, t=1 rehearsal, and both are
# kept because they discriminate between candidate repairs: BRepLib's
# SameParameter fixes the B-spline one and *not* the cylinder one, which is
# why the fix widens the shared vertex instead (tools/prototypes/RESULTS.md
# G12).

SELF_INTERSECTING_FACES = ("self-intersecting-wire-cylinder.brep",
                           "self-intersecting-wire-bspline.brep")


@pytest.mark.parametrize("name", SELF_INTERSECTING_FACES)
def test_the_self_intersecting_fixture_really_is_invalid(name):
    """Guard on the guard, and on the *mechanism* rather than just the symptom:
    a fixture that were invalid for some other reason would let the repair test
    below pass while proving nothing about self-intersection."""
    face = load_face(name)
    assert not occ.is_valid(face)
    # Invalid *only* contextually: nothing is wrong with any part on its own.
    assert all(
        occ.is_valid(e) for e in occ._explore(face, TopAbs_ShapeEnum.TopAbs_EDGE)
    )
    assert all(
        occ.is_valid(v) for v in occ._explore(face, TopAbs_ShapeEnum.TopAbs_VERTEX)
    )
    pair = occ._self_intersecting_pair(face)
    assert pair is not None, "the fault is a self-intersecting wire"
    assert occ._shared_vertices(*pair), "the reported pair is adjacent in the wire"


@pytest.mark.parametrize("name", SELF_INTERSECTING_FACES)
def test_widening_the_shared_vertex_makes_the_face_valid(name):
    face = load_face(name)
    before = occ.area(face)

    repaired, still_invalid = occ.fix_vertex_tolerances([face])

    assert (repaired, still_invalid) == (1, 0)
    assert occ.is_valid(face)
    # Same bound as the repair above: a tolerance is metadata, so the area must
    # be bit-identical rather than merely close.
    assert occ.area(face) == before


@pytest.mark.parametrize("name", SELF_INTERSECTING_FACES)
def test_the_repair_replaces_no_topology(name):
    """The property that makes this safe on an already-proven-watertight shell.
    ShapeFix_Shape also fixes these faces, and is rejected precisely because it
    mints new edge objects — the mechanism behind the seam-split regression."""
    face = load_face(name)
    before = {e.TShape() for e in occ._explore(face, TopAbs_ShapeEnum.TopAbs_EDGE)}
    verts = {v.TShape() for v in occ._explore(face, TopAbs_ShapeEnum.TopAbs_VERTEX)}

    occ.fix_vertex_tolerances([face])

    assert {e.TShape() for e in occ._explore(face, TopAbs_ShapeEnum.TopAbs_EDGE)} == before
    assert {
        v.TShape() for v in occ._explore(face, TopAbs_ShapeEnum.TopAbs_VERTEX)
    } == verts


@pytest.mark.parametrize("name", SELF_INTERSECTING_FACES)
def test_the_widening_is_bounded(name):
    """Its failure mode must be "leave the face for validate to report", never
    an unbounded tolerance (docs/algorithm.md §11)."""
    from OCP.BRep import BRep_Tool

    face = load_face(name)
    a, b = occ._self_intersecting_pair(face)
    shared = occ._shared_vertices(a, b)
    start = {v.TShape(): BRep_Tool.Tolerance_s(v) for v in shared}

    occ.fix_vertex_tolerances([face])

    for v in shared:
        grown = BRep_Tool.Tolerance_s(v)
        assert grown >= start[v.TShape()], "UpdateVertex never lowers a tolerance"
        assert grown <= start[v.TShape()] * occ.SELF_INTERSECT_TOL_GROWTH
        assert grown <= occ.SELF_INTERSECT_MAX_VERTEX_TOL


@pytest.mark.parametrize("name", SELF_INTERSECTING_FACES)
def test_widening_a_vertex_is_monotonically_permissive(name):
    """Why widening is safe on a shell where the vertex is shared with other
    faces: every check that reads a vertex tolerance is a "within tolerance"
    test, so a neighbour can only become more valid. Measured well past the
    repair's own cap rather than argued."""
    from OCP.BRep import BRep_Builder

    face = load_face(name)
    occ.fix_vertex_tolerances([face])
    assert occ.is_valid(face)
    before = occ.area(face)

    builder = BRep_Builder()
    for v in occ._explore(face, TopAbs_ShapeEnum.TopAbs_VERTEX):
        builder.UpdateVertex(TopoDS.Vertex_s(v), 25 * occ.SELF_INTERSECT_MAX_VERTEX_TOL)

    assert occ.is_valid(face)
    assert occ.area(face) == before


# --- ...and the same rung on a part whose trims are fat ----------------------
#
# `SpiralTest.step` at cc=5, t=1 carries the same false self-intersection, but
# its swept B-spline surface makes the boolean record tolerances two orders
# above the rehearsal's: the shared vertex below starts at 6.573e-02 mm, where
# the rehearsal's four are 8.7e-04 to 1.5e-03 mm (G12). That is *sixteen times*
# the fixed absolute cap the rung used to carry -- so the first candidate step
# was refused, nothing grew, and the run reported the face as residual without
# ever having tried it. A bound that can sit below the value it is bounding is
# not a bound, and this is the part that showed it.

FAT_VERTEX_FACE = "self-intersecting-wire-fat-vertex.brep"
SPIRAL_MAX_VERTEX_TOL = occ.SELF_INTERSECT_MAX_VERTEX_TOL_FRACTION * 1.0  # t = 1 mm


def test_the_fat_vertex_fixture_is_the_same_fault_at_a_different_scale():
    """Guard on the guard, exactly as for the two rehearsal fixtures: the same
    mechanism, so the repair below is not passing for some unrelated reason."""
    from OCP.BRep import BRep_Tool

    face = load_face(FAT_VERTEX_FACE)
    assert not occ.is_valid(face)
    assert all(occ.is_valid(e) for e in occ._explore(face, TopAbs_ShapeEnum.TopAbs_EDGE))
    assert all(occ.is_valid(v) for v in occ._explore(face, TopAbs_ShapeEnum.TopAbs_VERTEX))
    pair = occ._self_intersecting_pair(face)
    assert pair is not None, "the fault is a self-intersecting wire"
    shared = occ._shared_vertices(*pair)
    assert shared, "the reported pair is adjacent in the wire"
    assert BRep_Tool.Tolerance_s(shared[0]) > occ.SELF_INTERSECT_MAX_VERTEX_TOL, (
        "the point of this fixture: the kernel's own recorded tolerance is "
        "already above the fixed cap, so a fixed cap can never widen it"
    )


def test_the_fixed_cap_silently_disables_the_rung_on_a_fat_vertex():
    """The defect itself, pinned so it cannot come back as a default.

    Not a raise and not a wrong repair -- a repair that reports "still invalid"
    having done nothing, which is how it reached a released version.
    """
    face = load_face(FAT_VERTEX_FACE)
    before = occ.area(face)

    assert occ.fix_vertex_tolerances([face]) == (0, 1)

    assert not occ.is_valid(face)
    assert occ.area(face) == before


def test_a_cap_scaled_to_the_run_repairs_the_fat_vertex():
    face = load_face(FAT_VERTEX_FACE)
    before = occ.area(face)

    repaired, still_invalid = occ.fix_vertex_tolerances([face], SPIRAL_MAX_VERTEX_TOL)

    assert (repaired, still_invalid) == (1, 0)
    assert occ.is_valid(face)
    assert occ.area(face) == before


def test_the_fat_vertex_repair_stays_inside_both_bounds():
    """One 1.25x step is all it needs, so the *relative* bound is what governs
    here too -- the scaled absolute one only has to stop forbidding it."""
    from OCP.BRep import BRep_Tool

    face = load_face(FAT_VERTEX_FACE)
    a, b = occ._self_intersecting_pair(face)
    shared = occ._shared_vertices(a, b)
    start = {v.TShape(): BRep_Tool.Tolerance_s(v) for v in shared}

    occ.fix_vertex_tolerances([face], SPIRAL_MAX_VERTEX_TOL)

    for v in shared:
        grown = BRep_Tool.Tolerance_s(v)
        assert grown >= start[v.TShape()], "UpdateVertex never lowers a tolerance"
        assert grown <= start[v.TShape()] * occ.SELF_INTERSECT_TOL_GROWTH
        assert grown <= SPIRAL_MAX_VERTEX_TOL


def test_the_fat_vertex_repair_replaces_no_topology():
    face = load_face(FAT_VERTEX_FACE)
    before = {e.TShape() for e in occ._explore(face, TopAbs_ShapeEnum.TopAbs_EDGE)}
    verts = {v.TShape() for v in occ._explore(face, TopAbs_ShapeEnum.TopAbs_VERTEX)}

    occ.fix_vertex_tolerances([face], SPIRAL_MAX_VERTEX_TOL)

    assert {e.TShape() for e in occ._explore(face, TopAbs_ShapeEnum.TopAbs_EDGE)} == before
    assert {v.TShape() for v in occ._explore(face, TopAbs_ShapeEnum.TopAbs_VERTEX)} == verts


def test_a_capped_widening_leaves_the_face_for_validate(monkeypatch):
    """With no headroom the repair must give up quietly and report the face as
    residual, not raise and not loop forever."""
    monkeypatch.setattr(occ, "SELF_INTERSECT_TOL_GROWTH", 1.0)

    face = load_face(SELF_INTERSECTING_FACES[0])
    before = occ.area(face)

    assert occ.fix_vertex_tolerances([face]) == (0, 1)
    assert not occ.is_valid(face)
    assert occ.area(face) == before
