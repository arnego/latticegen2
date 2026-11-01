"""Frustum supports on the lattice grid (docs/algorithm.md §14.3).

The grid rule, both manifold repairs and the fine grid are checked against
their definitions by brute force; whole supports against properties that hold
by construction -- every face on the body or on a lattice plane, nothing
outside the body or inside a strut, no edge used by more than two faces --
asked of the solids, never of the plan that built them.
"""

from __future__ import annotations

import itertools
import os

import numpy as np
import pytest

from latticegen2.lattice import lattice_params

NO_NODES = np.empty((0, 3), dtype=np.int64)


def _box(lo, size):
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt

    from latticegen2 import occ

    return occ.solids(BRepPrimAPI_MakeBox(gp_Pnt(*lo), *size).Shape())[0]


def _plane(face):
    """``(unit normal, offset)`` of a planar face's plane, ``n . x = offset``."""
    from OCP.BRepAdaptor import BRepAdaptor_Surface

    pl = BRepAdaptor_Surface(face).Plane()
    d, o = pl.Axis().Direction(), pl.Location()
    n = np.array([d.X(), d.Y(), d.Z()])
    return n, float(n @ np.array([o.X(), o.Y(), o.Z()]))


# --- the grid ----------------------------------------------------------------------


def test_the_grid_is_the_finest_whose_rod_gap_holds_the_minimum_thickness():
    from latticegen2.support.frustum import grid_divisions

    lp = lattice_params(10.0, 1.5)
    m = grid_divisions(lp, 1.5)
    assert m == 2
    assert lp.a / m - lp.t >= 1.5 > lp.a / (m + 1) - lp.t
    assert grid_divisions(lattice_params(20.0, 4.0), 4.0) == 1
    # Never below one, even where no gap could hold the thickness.
    assert grid_divisions(lattice_params(10.0, 6.0), 6.0) == 1


# --- manifold repair -----------------------------------------------------------------


def _remaining(S, L):
    """Edge and corner contacts of support nodes ``S`` that a fill could still resolve:
    brute force over every 2x2 and 2x2x2 window, against the definition."""
    S = np.pad(S, 1)
    L = np.pad(L, 1)
    n = 0
    for ax in range(3):
        o = [i for i in range(3) if i != ax]
        for c in itertools.product(*[range(s - 1) for s in S.shape]):
            def node(di, dj):
                p = list(c)
                p[o[0]] += di
                p[o[1]] += dj
                return tuple(p)
            A, B, C, D = node(0, 0), node(1, 0), node(0, 1), node(1, 1)
            if S[A] and S[D] and not S[B] and not S[C] and not (L[B] and L[C]):
                n += 1
            if S[B] and S[C] and not S[A] and not S[D] and not (L[A] and L[D]):
                n += 1
    for c in itertools.product(*[range(s - 1) for s in S.shape]):
        cube = {d: tuple(np.add(c, d)) for d in itertools.product((0, 1), repeat=3)}
        for d0 in cube:
            d1 = tuple(1 - x for x in d0)
            if d0 > d1 or not (S[cube[d0]] and S[cube[d1]]):
                continue
            if any(S[cube[d]] for d in cube if d not in (d0, d1)):
                continue
            steps = [tuple((d1[k] - d0[k]) if k == i else 0 for k in range(3)) for i in range(3)]
            for a1, a2, _ in itertools.permutations(steps):
                p1 = tuple(np.add(d0, a1))
                p2 = tuple(np.add(p1, a2))
                if not L[cube[p1]] and not L[cube[p2]]:
                    n += 1
                    break
    return n


def test_a_diagonal_pair_is_joined_through_the_nodes_above_the_lower_one():
    from latticegen2.support.frustum import make_manifold

    S = np.zeros((2, 2, 1), dtype=bool)
    S[0, 0, 0] = S[1, 1, 0] = True
    got, filled, stuck = make_manifold(S)
    assert got[1, 0, 0] and got[0, 1, 0] and (filled, stuck) == (2, 0)


def test_beside_a_rod_the_node_below_is_filled_instead():
    """Support on two sides of a rod's corner line, the rod on the third and
    free space on the fourth: the cylinder's top support had 69 such edges."""
    from latticegen2.support.frustum import make_manifold

    S = np.zeros((2, 2, 1), dtype=bool)
    L = np.zeros_like(S)
    S[1, 0, 0] = S[0, 1, 0] = True
    L[1, 1, 0] = True
    got, filled, stuck = make_manifold(S, L)
    assert got[0, 0, 0] and not got[1, 1, 0] and (filled, stuck) == (1, 0)


def test_a_corner_contact_is_joined_by_a_path_round_the_rods():
    from latticegen2.support.frustum import make_manifold

    S = np.zeros((2, 2, 2), dtype=bool)
    L = np.zeros_like(S)
    S[0, 0, 0] = S[1, 1, 1] = True
    L[1, 0, 0] = True
    got, _filled, stuck = make_manifold(S, L)
    assert stuck == 0 and not (got & L).any()
    assert _remaining(got, L) == 0


@pytest.mark.parametrize("seed", range(6))
def test_the_repair_leaves_no_contact_it_could_fill(seed):
    from latticegen2.support.frustum import make_manifold

    rng = np.random.default_rng(seed)
    S = rng.random((7, 6, 5)) < 0.3
    L = (rng.random(S.shape) < 0.12) & ~S
    got, _filled, _stuck = make_manifold(S, L)
    assert (got >= S).all() and not (got & L).any()
    assert _remaining(got, L) == 0


# --- the fine grid ------------------------------------------------------------------


def test_rods_are_whole_nodes_of_the_fine_grid():
    """Every rod face and cap plane is a plane of the grid: a fine support
    node's centre and the points just inside each of its corners agree about
    the lattice."""
    from latticegen2.support.assemble import _LocalLattice
    from latticegen2.support.frustum import _Grid, _planes

    lp = lattice_params(10.0, 1.5, (35.0, 20.0, -110.0))
    nodes = np.array(list(itertools.product(range(-1, 2), repeat=3)))
    local = _LocalLattice(lp, nodes)
    g = _Grid(lp, 2)
    coords = [_planes(g, -6, 12) for _ in range(3)]
    centres = [(c[:-1] + c[1:]) / 2.0 for c in coords]
    C = np.stack(np.meshgrid(*centres, indexing="ij"), -1).reshape(-1, 3)
    half = [np.diff(c) / 2.0 for c in coords]
    H = np.stack(np.meshgrid(*half, indexing="ij"), -1).reshape(-1, 3)
    at_centre = local.contains(g.world(C))
    assert at_centre.any() and not at_centre.all()
    for corner in itertools.product((-1, 1), repeat=3):
        P = C + (1.0 - 1e-6) * H * np.array(corner)
        assert np.array_equal(local.contains(g.world(P)), at_centre)


def test_a_staircase_has_only_lattice_plane_faces_and_its_nodes_volume():
    from latticegen2 import occ
    from latticegen2.support.frustum import _Grid, _planes, _staircase

    lp = lattice_params(10.0, 1.5)
    g = _Grid(lp, 2)
    coords = [_planes(g, 0, 4) for _ in range(3)]
    shape = tuple(len(c) - 1 for c in coords)
    S = np.zeros(shape, dtype=bool)
    S[1:4, 1:3, 0:2] = True
    S[2, 2:4, 1:4] = True
    solids = _staircase(S, coords, g)
    assert len(solids) == 1 and occ.is_valid(solids[0])
    widths = [np.diff(c) for c in coords]
    nodes = sum(widths[0][i] * widths[1][j] * widths[2][k] for i, j, k in np.argwhere(S))
    assert occ.volume(solids[0]) == pytest.approx(nodes, rel=1e-12)
    E = np.asarray(lp.e)
    for f in occ.faces(solids[0]):
        n, _off = _plane(f)
        assert np.isclose(np.abs(E @ n).max(), 1.0, atol=1e-9)



def test_a_staircase_splits_on_a_grid_plane_into_halves_that_add_up():
    from latticegen2 import occ
    from latticegen2.support.frustum import _Grid, _planes, _split, _staircase

    lp = lattice_params(10.0, 1.5)
    g = _Grid(lp, 2)
    coords = [_planes(g, 0, 6) for _ in range(3)]
    shape = tuple(len(c) - 1 for c in coords)
    S = np.zeros(shape, dtype=bool)
    S[0:shape[0], 1:3, 1:3] = True
    whole = _staircase(S, coords, g)[0]
    halves = _split(whole, g)
    assert halves is not None and len(halves) == 2
    assert sum(occ.volume(h) for h in halves) == pytest.approx(occ.volume(whole), rel=1e-9)
    E = np.asarray(lp.e)
    for h in halves:
        assert occ.is_valid(h)
        for f in occ.faces(h):
            n, _off = _plane(f)
            assert np.isclose(np.abs(E @ n).max(), 1.0, atol=1e-9)


def test_a_body_no_trim_passes_is_trimmed_in_halves_and_joined_again(
        chamber, monkeypatch, tmp_path):
    """The whole body's intersection is refused, as SpiralTest's was (0.1 mm^3
    where sampling found 1,627.5): its halves are trimmed instead, and the
    result is one body again, of the same volume."""
    from latticegen2.support import frustum
    from latticegen2.support.generate import supports_for

    body, mesh, lp, kept, occ = chamber
    real = frustum._trim
    refused = []

    def first_refused(piece, *args, **kw):
        if not refused:
            refused.append(args[-1])
            return None, ["refused"], 0.0
        return real(piece, *args, **kw)

    monkeypatch.setattr(frustum, "_trim", first_refused)
    solids, stats, lines = supports_for(body, mesh, [], kept, lp, 1.5, tmpdir=str(tmp_path))
    assert refused and any("split in 2" in line for line in lines)
    assert len(solids) == stats["support_bodies"]
    assert all(occ.is_valid(s) and frustum.non_manifold_edges(s) == 0 for s in solids)
    monkeypatch.setattr(frustum, "_trim", real)
    (tmp_path / "plain").mkdir()
    plain, _st, _l = supports_for(body, mesh, [], kept, lp, 1.5, tmpdir=str(tmp_path / "plain"))
    assert len(solids) == len(plain)
    assert sum(occ.volume(s) for s in solids) == pytest.approx(
        sum(occ.volume(s) for s in plain), rel=1e-6)

# --- whole supports ----------------------------------------------------------------------


@pytest.fixture(scope="module")
def chamber():
    """A box chamber 30 x 30 x 20 mm with its lattice: every node inside it. Its
    junctions reach out of the box, as boundary junctions do before their trim;
    a support meets them only inside the body."""
    from latticegen2 import lattice as L
    from latticegen2 import occ
    from latticegen2.classify import tessellate_surface

    lp = lattice_params(10.0, 1.5)
    body = _box((0.0, 0.0, 0.0), (30.0, 30.0, 20.0))
    mesh = tessellate_surface(body, lp)
    lo, hi = np.array([0.0, 0.0, 0.0]), np.array([30.0, 30.0, 20.0])
    span = np.arange(-12, 13)
    idx = np.array(list(itertools.product(span, span, span)))
    pos = L.nodes(lp, idx)
    kept = idx[np.all((pos > lo) & (pos < hi), axis=1)]
    assert len(kept) > 20
    return body, mesh, lp, kept, occ


@pytest.fixture(scope="module")
def chamber_supports(chamber, tmp_path_factory):
    from latticegen2.support.generate import supports_for

    body, mesh, lp, kept, _occ = chamber
    return supports_for(body, mesh, [], kept, lp, 1.5,
                        tmpdir=str(tmp_path_factory.mktemp("frustum")))


def test_a_chamber_lid_is_carried_down_to_the_rods(chamber, chamber_supports):
    body, mesh, lp, kept, occ = chamber
    solids, stats, _lines = chamber_supports
    assert solids and stats["support_grid_divisions"] == 2
    # The rods are near: few if any ceiling points need the wall.
    assert stats["support_wall_anchors"] < stats["support_ceiling_points"] / 2
    assert stats["support_ceiling_points_uncovered"] == 0
    for s in solids:
        assert occ.is_valid(s)


def test_chamber_supports_are_manifold_and_overlap_no_strut(chamber, chamber_supports):
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Common

    from latticegen2.support.assemble import _LocalLattice
    from latticegen2.support.build import _boolean
    from latticegen2.support.frustum import non_manifold_edges

    body, mesh, lp, kept, occ = chamber
    solids, _stats, _lines = chamber_supports
    local = _LocalLattice(lp, kept)
    for s in solids:
        assert non_manifold_edges(s) == 0
        near = local.near(occ.bounding_box(s))
        got = _boolean(BRepAlgoAPI_Common, s, near)
        assert sum(occ.volume(x) for x in occ.solids(got)) < 1e-6


def test_chamber_supports_pass_their_gates_with_every_face_on_a_lattice_plane(
        chamber, chamber_supports):
    """Self-support by construction: a face not on the body lies on a plane of
    the lattice frame, so it overhangs exactly the rod angle or less."""
    from latticegen2.support.build import Gates

    body, mesh, lp, kept, occ = chamber
    solids, _stats, _lines = chamber_supports
    gates = Gates(mesh, lp.b)
    E = np.asarray(lp.e)
    walls = [(np.eye(3)[k] * s, w) for k, top in enumerate((30.0, 30.0, 20.0))
             for s, w in ((-1.0, 0.0), (1.0, top))]
    for s in solids:
        fails, _gs = gates.check(s)
        assert fails == []
        for f in occ.faces(s):
            n, off = _plane(f)
            on_body = any(abs(abs(n @ wn) - 1.0) < 1e-9 and abs(abs(off) - w) < 1e-6
                          for wn, w in walls)
            assert on_body or np.isclose(np.abs(E @ n).max(), 1.0, atol=1e-9)


def test_a_piece_no_trim_produces_is_dropped_and_reported(chamber, monkeypatch, tmp_path):
    """Where the kernel trims a piece at no rung and no size -- three pieces of
    SpiralTest at -90 45 0 -- the piece is dropped, named, and the ceiling
    points it held are counted as left unbacked; the run goes on (§14.3 S6)."""
    from latticegen2.support import frustum
    from latticegen2.support.generate import supports_for

    body, mesh, lp, kept, occ = chamber
    real = frustum._trim

    def refuse_first_body(piece, *args, **kw):
        if args[-1].startswith("body-0"):
            return None, ["refused"], 123.0
        return real(piece, *args, **kw)

    monkeypatch.setattr(frustum, "_trim", refuse_first_body)
    monkeypatch.setattr(frustum, "SPLIT_DEPTH", 0)
    solids, stats, lines = supports_for(body, mesh, [], kept, lp, 1.5, tmpdir=str(tmp_path))
    assert stats["support_pieces_dropped"] == 1
    assert stats["support_dropped_volume_mm3"] == pytest.approx(123.0)
    warning = [line for line in lines if line.startswith(frustum.DROPPED_WARNING)]
    assert len(warning) == 1 and "123.0 mm^3" in warning[0]
    assert any(line.startswith("  body-0 near [") for line in lines)
    assert stats["support_ceiling_points_uncovered"] > 0

def test_supports_are_assembled_beside_the_lattice(chamber, chamber_supports):
    from latticegen2.support import assemble

    body, mesh, lp, kept, occ = chamber
    solids, _stats, _lines = chamber_supports
    got = assemble.assemble([], solids, lp, kept)
    assert got.support_volume == pytest.approx(sum(occ.volume(s) for s in solids))


def _box(lo, hi):
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt

    return BRepPrimAPI_MakeBox(gp_Pnt(*lo), gp_Pnt(*hi)).Solid()


def test_a_lost_junction_piece_is_not_lattice():
    """A strut tip the floating-body rule dropped from a kept node: neither the
    point test nor the junction instances may claim it is there."""
    from latticegen2 import occ
    from latticegen2.support.assemble import _LocalLattice

    lp = lattice_params(10.0, 1.5, (0.0, 0.0, 0.0))
    nodes = np.zeros((1, 3), dtype=np.int64)
    tip = 3.0 * lp.e[0]
    lost = _box(tip - 1.5, tip + 1.5)
    whole = _LocalLattice(lp, nodes)
    local = _LocalLattice(lp, nodes, [lost])

    assert whole.contains(tip[None])[0] and whole.in_template(tip[None])[0]
    assert not local.contains(tip[None])[0] and local.in_template(tip[None])[0]
    hub = 0.1 * lp.e[1]
    assert local.contains(hub[None])[0]

    box = (tip - 5.0, tip + 5.0)
    full = sum(occ.volume(s) for s in whole.near(box))
    less = sum(occ.volume(s) for s in local.near(box))
    assert 0.0 < less < full


def test_a_fine_node_is_freed_only_if_no_kept_lattice_is_in_it():
    """A wall can cross a rod's fine node with a lost tip on one side and kept
    lattice on the other: support there would overlap the lattice, so that node
    stays lattice, while one wholly inside the lost piece is freed."""
    from latticegen2.support.assemble import _LocalLattice
    from latticegen2.support.frustum import _free_of_lattice, _Grid, _planes, _staircase

    lp = lattice_params(10.0, 1.5, (0.0, 0.0, 0.0))
    g = _Grid(lp, 2)
    coords = [_planes(g, -2, 4) for _ in range(3)]
    centres = [(c[:-1] + c[1:]) / 2.0 for c in coords]
    C = np.stack(np.meshgrid(*centres, indexing="ij"), -1)
    # The rod along e_0: fine nodes on its axis, past the hub, in order.
    on_axis = (np.abs(C[..., 1]) < 1e-9) & (np.abs(C[..., 2]) < 1e-9) & (C[..., 0] > lp.t)
    inner, outer = [tuple(i) for i in np.argwhere(on_axis)[np.argsort(C[on_axis][:, 0])][:2]]
    # All of the rod from the middle of the first of them onward is lost.
    start = (coords[0][inner[0]] + coords[0][inner[0] + 1]) / 2.0
    cube = [np.array([start, coords[0][-1] + 1.0]),
            np.array([coords[1][0] - 1.0, coords[1][-1] + 1.0]),
            np.array([coords[2][0] - 1.0, coords[2][-1] + 1.0])]
    lost = _staircase(np.ones((1, 1, 1), dtype=bool), cube, g)[0]
    local = _LocalLattice(lp, np.zeros((1, 3), dtype=np.int64), [lost])
    L0 = np.ones(tuple(len(c) for c in centres), dtype=bool)

    L = _free_of_lattice(g, coords, local, lambda P: np.ones(len(P), dtype=bool), L0)

    assert L[inner], "kept lattice remains in the node the wall crosses"
    assert not L[outer], "nothing of the lattice is left in the node beyond it"


def test_a_keel_is_supported_though_no_face_of_it_is_ceiling():
    """The flanks of test/support-gaps.step's keel are 45 deg, under the bar,
    so its concave edge alone asks for support, and gets it."""
    from latticegen2 import occ
    from latticegen2.classify import tessellate_surface
    from latticegen2.support import overhang
    from latticegen2.support.frustum import exact_inside
    from latticegen2.support.generate import supports_for

    lp = lattice_params(10.0, 1.5, (0.0, 0.0, 0.0))
    body = occ.read_step(os.path.join(os.path.dirname(__file__), "support-gaps.step"))
    mesh = tessellate_surface(body, lp)
    samples = overhang.detect(body, lp.b)
    keel = samples.points[samples.kind == overhang.EDGE]
    assert len(keel) > 10
    solids, _stats, _lines = supports_for(body, mesh, [], np.empty((0, 3), dtype=np.int64),
                                          lp, 1.5, samples=samples)
    backed = np.zeros(len(keel), dtype=bool)
    for s in solids:
        backed |= exact_inside(s, keel - 2e-3 * lp.b)
    assert backed.all(), f"{int((~backed).sum())} of {len(keel)} keel samples unbacked"
