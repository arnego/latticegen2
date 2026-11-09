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


def test_the_grid_is_the_finest_whose_thinnest_member_holds_the_thickness():
    from latticegen2.cli import support_grid
    from latticegen2.errors import ProcessingError
    from latticegen2.support.frustum import grid_divisions

    lp = lattice_params(10.0, 1.5)
    m = grid_divisions(lp, 1.5)
    assert m == 3 == support_grid(lp.a, lp.t, 1.5)[0]
    assert lp.a / m - lp.t / 2 >= 1.5 > lp.a / (m + 1) - lp.t / 2
    assert grid_divisions(lattice_params(20.0, 4.0), 4.0) == 2
    # m = 1 holds a - t; beyond that nothing does, and the generator refuses
    # rather than shipping thinner support than was asked for.
    lp = lattice_params(10.0, 6.0)
    assert grid_divisions(lp, lp.a - lp.t) == 1
    with pytest.raises(ProcessingError):
        grid_divisions(lp, 6.0)


@pytest.mark.parametrize("m", range(1, 7))
def test_a_support_node_carries_a_rod_on_at_most_one_side_per_axis(m):
    """The lemma the thinnest-member figure rests on (cli.support_grid). Rods
    lie on every m-th grid plane, so at m >= 2 no support node has rod planes
    on both of its sides along any axis; at m = 1 every node has both."""
    for i in range(-2 * m, 2 * m):
        sides = (i % m == 0) + ((i + 1) % m == 0)
        assert sides == 2 if m == 1 else sides <= 1


@pytest.mark.parametrize("cc, t, m", [(10.0, 1.5, 3), (7.0, 1.4, 2), (10.0, 1.5, 1)])
def test_the_thinnest_member_beside_the_rods_is_what_the_rule_says(cc, t, m):
    """A support node at a lattice node, less every rod through it, built and
    measured rather than argued: its volume is the node less a t/2 quadrant
    per rod, and the chord across its thinnest arm is the figure
    ``cli.support_grid`` reports -- q - t/2 at m >= 2, a - t at m = 1."""
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.TopAbs import TopAbs_IN
    from OCP.TopTools import TopTools_ListOfShape
    from OCP.gp import gp_Pnt

    from latticegen2 import occ
    from latticegen2.cli import support_grid

    a = cc / np.sqrt(2.0)
    q, h = a / m, t / 2
    planes = [c for c in (0.0, a) if c <= q + 1e-9]      # rod axes on the node
    objs, tools = TopTools_ListOfShape(), TopTools_ListOfShape()
    objs.Append(BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), q, q, q).Shape())
    n_rods = 0
    for ax in range(3):
        o = [k for k in range(3) if k != ax]
        for c0, c1 in itertools.product(planes, planes):
            lo, size = [0.0] * 3, [t] * 3
            lo[ax], size[ax] = -1.0, q + 2.0
            lo[o[0]], lo[o[1]] = c0 - h, c1 - h
            tools.Append(BRepPrimAPI_MakeBox(gp_Pnt(*lo), *size).Shape())
            n_rods += 1
    cut = BRepAlgoAPI_Cut()
    cut.SetArguments(objs)
    cut.SetTools(tools)
    cut.Build()
    support = cut.Shape()
    assert n_rods == (12 if m == 1 else 3)
    # Inclusion-exclusion at each corner where three rods meet: -3 h^3 + h^3.
    corners = 8 if m == 1 else 1
    assert occ.volume(support) == pytest.approx(
        q ** 3 - n_rods * h * h * q + corners * 2 * h ** 3, rel=1e-9)

    # The arm beside the rod along e_2, at mid-height: a chord along y at
    # x = h/2, from the rod's face to the node's far face (or the next rod).
    ys = np.linspace(0.0, q, 4001)
    inside = [BRepClass3d_SolidClassifier(support, gp_Pnt(h / 2, float(y), q / 2), 1e-9)
              .State() == TopAbs_IN for y in ys]
    span = ys[np.nonzero(inside)[0]]
    chord = span.max() - span.min()
    want = a - t if m == 1 else q - h
    assert chord == pytest.approx(want, abs=2 * q / 4000)
    # And the rule, asked for exactly that thickness, picks this m and says so.
    got_m, thinnest = support_grid(a, t, want)
    assert got_m == m and thinnest == pytest.approx(want)


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


@pytest.mark.parametrize("origin, fills", [(0.0, 0), (-1.0, 1)])
def test_boxes_passing_a_rod_on_its_lower_sides_meet_inside_it(monkeypatch, origin, fills):
    """Two support nodes passing a rod on its two lower sides, as boxes from
    two neighbouring rods do. On a grid through the rods' axes they meet along
    the axis, inside the rod, and with the rod left out they do not meet at
    all: nothing to repair. On a grid on the rods' -t/2 faces they meet along
    the rod's lowest corner line, outside it, and the repair fills the node
    below them -- a node with nothing under it (specification.md §6.2's
    final-print check found 33 on the ball and 110 on the cylinder)."""
    from latticegen2.support import frustum

    monkeypatch.setattr(frustum, "GRID_ORIGIN", origin)
    lp = lattice_params(10.0, 1.5)
    g = frustum._Grid(lp, 2)
    coords = [frustum._planes(g, -1, 2), frustum._planes(g, -1, 2), frustum._planes(g, 0, 1)]
    centres = [(c[:-1] + c[1:]) / 2.0 for c in coords]
    C = np.stack(np.meshgrid(*centres, indexing="ij"), -1)
    coarse = np.floor((C - g.o) / g.q).astype(np.int64)
    S = (np.all(coarse == (0, -1, 0), axis=-1) | np.all(coarse == (-1, 0, 0), axis=-1))
    L = (np.abs(C[..., 0]) < g.h) & (np.abs(C[..., 1]) < g.h)   # the rod along e_2
    S &= ~L
    assert (len(frustum._contacts(S)) > 0) == (fills > 0)
    _got, filled, stuck = frustum.make_manifold(S, L)
    assert stuck == 0
    assert filled >= fills and (filled == 0) == (fills == 0)


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


@pytest.mark.parametrize("seed", range(6))
def test_contacts_are_what_the_repair_leaves(seed):
    from latticegen2.support.frustum import _contacts, make_manifold

    rng = np.random.default_rng(seed)
    S = rng.random((7, 6, 5)) < 0.3
    L = (rng.random(S.shape) < 0.3) & ~S
    got, _filled, stuck = make_manifold(S, L)
    assert (stuck > 0) == bool(_contacts(got))
    for nodes, upper in _contacts(got):
        assert got[tuple(upper)] and any(tuple(n) == tuple(upper) for n in nodes.tolist())


def _stuck_pair():
    """A diagonal pair with rod nodes on both of the other two sides."""
    S = np.zeros((2, 2, 1), dtype=bool)
    L = np.zeros_like(S)
    S[0, 0, 0] = S[1, 1, 0] = True
    L[1, 0, 0] = L[0, 1, 0] = True
    return S, L


def test_a_contact_beside_a_rod_node_the_body_does_not_reach_is_filled(monkeypatch):
    """Three such contacts made the rehearsal part's 58,000 mm3 support body
    invalid, so that no trim and no split could take it."""
    from latticegen2.support import frustum

    monkeypatch.setattr(frustum, "_holds_kept", lambda *a: False)
    S, L = _stuck_pair()
    got, filled, removed = frustum._resolve_contacts(S, L, None, None, None, None, None)
    assert removed == 0 and filled >= 1 and got[0, 0, 0] and got[1, 1, 0]
    assert not frustum._contacts(got)


def test_a_contact_between_real_rods_loses_its_upper_node(monkeypatch):
    from latticegen2.support import frustum

    monkeypatch.setattr(frustum, "_holds_kept", lambda *a: True)
    S, L = _stuck_pair()
    got, _filled, removed = frustum._resolve_contacts(S, L, None, None, None, None, None)
    assert removed == 1 and got[0, 0, 0] and not got[1, 1, 0] and not (got & L).any()
    assert not frustum._contacts(got)


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
            return None, ["refused"], (0.0, 0.0)
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
    assert solids and stats["support_grid_divisions"] == 3
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
            return None, ["refused"], (123.0, 1.0)
        return real(piece, *args, **kw)

    monkeypatch.setattr(frustum, "_trim", refuse_first_body)
    monkeypatch.setattr(frustum, "SPLIT_DEPTH", 0)
    solids, stats, lines = supports_for(body, mesh, [], kept, lp, 1.5, tmpdir=str(tmp_path))
    assert stats["support_pieces_dropped"] == 1
    assert stats["support_dropped_volume_mm3"] == pytest.approx(123.0)
    warning = [line for line in lines
               if line.startswith(f"WARNING: {frustum.DROPPED_WARNING}")]
    assert len(warning) == 1 and "123.0 mm^3" in warning[0]
    assert "\n  body-0 near [" in warning[0]
    assert stats["support_ceiling_points_uncovered"] > 0


def _sampled_empty_first_body(frustum, monkeypatch):
    real = frustum._trim

    def empty_first_body(piece, *args, **kw):
        if args[-1].startswith("body-0"):
            return [], ["sampled empty"], (0.0, 1.0)
        return real(piece, *args, **kw)

    monkeypatch.setattr(frustum, "_trim", empty_first_body)
    monkeypatch.setattr(frustum, "SPLIT_DEPTH", 0)


def test_a_piece_sampling_finds_empty_that_holds_ceiling_is_dropped_and_reported(
        chamber, monkeypatch, tmp_path):
    """Sampling can miss a sliver of support, and the coverage count reads the
    grid the piece was built from: going quietly, its ceiling points would be
    counted as backed with no material shipped under them."""
    from latticegen2.support import frustum
    from latticegen2.support.generate import supports_for

    body, mesh, lp, kept, occ = chamber
    _sampled_empty_first_body(frustum, monkeypatch)
    solids, stats, lines = supports_for(body, mesh, [], kept, lp, 1.5, tmpdir=str(tmp_path))
    assert stats["support_pieces_dropped"] == 1
    warning = [line for line in lines
               if line.startswith(f"WARNING: {frustum.DROPPED_WARNING}")]
    assert len(warning) == 1 and "sampling found no material" in warning[0]
    assert stats["support_ceiling_points_uncovered"] > 0


def test_a_piece_sampling_finds_empty_that_backs_nothing_is_discarded_quietly(
        chamber, monkeypatch, tmp_path):
    from latticegen2.support import frustum
    from latticegen2.support.generate import supports_for

    body, mesh, lp, kept, occ = chamber
    _sampled_empty_first_body(frustum, monkeypatch)
    real_inside = frustum._tessellated_inside
    first = []

    def nothing_in_first(piece):
        if not first:  # body-0 is the first piece asked
            first.append(piece)
            return lambda pts: np.zeros(len(pts), dtype=bool)
        return real_inside(piece)

    monkeypatch.setattr(frustum, "_tessellated_inside", nothing_in_first)
    solids, stats, lines = supports_for(body, mesh, [], kept, lp, 1.5, tmpdir=str(tmp_path))
    assert stats["support_pieces_dropped"] == 0
    assert not any(line.startswith(f"WARNING: {frustum.DROPPED_WARNING}") for line in lines)


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
    # The rod along e_0: its fine nodes in one quadrant round its axis, which
    # the support grid's planes run through, past the hub, in order.
    q4 = lp.t / 4.0
    on_axis = (np.abs(C[..., 1] - q4) < 1e-9) & (np.abs(C[..., 2] - q4) < 1e-9) & (C[..., 0] > lp.t)
    node = tuple(np.argwhere(on_axis)[np.argsort(C[on_axis][:, 0])][0])
    L0 = np.ones(tuple(len(c) for c in centres), dtype=bool)
    body = _box((-50.0, -50.0, -50.0), (50.0, 50.0, 50.0))

    def stays_lattice(start):
        """Whether the node stays lattice with all of the rod from ``start`` on lost."""
        cube = [np.array([start, coords[0][-1] + 1.0]),
                np.array([coords[1][0] - 1.0, coords[1][-1] + 1.0]),
                np.array([coords[2][0] - 1.0, coords[2][-1] + 1.0])]
        lost = _staircase(np.ones((1, 1, 1), dtype=bool), cube, g)[0]
        local = _LocalLattice(lp, np.zeros((1, 3), dtype=np.int64), [lost])
        return _free_of_lattice(g, coords, local, lambda P: np.ones(len(P), dtype=bool),
                                body, L0)[node]

    lo, hi = coords[0][node[0]], coords[0][node[0] + 1]
    assert stays_lattice((lo + hi) / 2.0), "kept lattice remains in the node the wall crosses"
    assert not stays_lattice(lo - 0.01), "nothing of the lattice is left in a node wholly lost"


@pytest.mark.parametrize("wedge", [True, False])
def test_a_wedge_of_kept_lattice_between_the_samples_keeps_its_node(wedge):
    """A node the samples call free can still hold kept lattice, across a wall from
    the lost tip: on the rehearsal part a 0.003 mm3 wedge became an overlap."""
    from latticegen2 import occ
    from latticegen2.support.assemble import _LocalLattice
    from latticegen2.support.frustum import FREE_SAMPLES, _Grid, _holds_kept, _planes, _staircase

    lp = lattice_params(10.0, 1.5, (0.0, 0.0, 0.0))
    g = _Grid(lp, 2)
    coords = [_planes(g, -2, 4) for _ in range(3)]
    centres = [(c[:-1] + c[1:]) / 2.0 for c in coords]
    C = np.stack(np.meshgrid(*centres, indexing="ij"), -1)
    q4 = lp.t / 4.0
    on_axis = (np.abs(C[..., 1] - q4) < 1e-9) & (np.abs(C[..., 2] - q4) < 1e-9) & (C[..., 0] > lp.t)
    node = tuple(np.argwhere(on_axis)[np.argsort(C[on_axis][:, 0])][0])
    lo = np.array([coords[a][node[a]] for a in range(3)])
    hi = np.array([coords[a][node[a] + 1] for a in range(3)])
    side = 0.4 * (hi - lo).min() / FREE_SAMPLES     # short of the first sample
    wall = lo[0] + 0.3 * (hi - lo)[0]

    def box(a, b):
        cube = [np.array([a[k], b[k]]) for k in range(3)]
        return _staircase(np.ones((1, 1, 1), dtype=bool), cube, g)[0]

    tip = box([wall, lo[1] - 1.0, lo[2] - 1.0], [hi[0] + 1.0, hi[1] + 1.0, hi[2] + 1.0])
    chip = box(lo, lo + side)
    body = occ.compound([tip, chip]) if wedge else tip
    local = _LocalLattice(lp, np.zeros((1, 3), dtype=np.int64), [tip])

    def in_body(P):
        return g.cube(P)[:, 0] >= wall

    assert _holds_kept(g, coords, local, in_body, body, node) is wedge


def test_an_overlap_that_is_a_lost_piece_is_not_one():
    from latticegen2.support.assemble import _LocalLattice
    from latticegen2.support.frustum import _Grid, _staircase

    lp = lattice_params(10.0, 1.5, (0.0, 0.0, 0.0))
    g = _Grid(lp, 2)

    def box(a, b):
        return _staircase(np.ones((1, 1, 1), dtype=bool),
                          [np.array([a, b]) for _ in range(3)], g)[0]

    local = _LocalLattice(lp, np.zeros((1, 3), dtype=np.int64), [box(0.0, 1.0)])
    assert local.is_lost(box(0.25, 0.75))
    assert not local.is_lost(box(2.0, 2.5))


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


def test_a_ceiling_clipping_a_node_corner_gets_a_ceiling_point():
    """A triangle with no vertex in a node, clipping its corner: that node needs
    a ceiling point. One under the bar, or on a strut, does not."""
    from latticegen2.support.frustum import _clipped_ceiling, _Grid

    lp = lattice_params(10.0, 1.5, (0.0, 0.0, 0.0))
    g = _Grid(lp, 2)
    # In cube coordinates, across the corner of node (1, 1, 1) nearest the origin.
    corner = np.array([1, 1, 1]) * g.q + g.o
    tri = corner + np.array([[0.2, -0.1, -0.05], [-0.1, 0.2, -0.05], [-0.05, -0.1, 0.2]])
    verts = g.world(tri)
    tris = np.array([[0, 1, 2]])
    b = np.asarray(lp.b, dtype=float)
    up = np.broadcast_to(b, (3, 3))
    nowhere = lambda P: np.zeros(len(P), dtype=bool)                    # noqa: E731
    held = {tuple(x) for x in g.node_of(verts).tolist()}
    assert (1, 1, 1) not in held

    p, n = _clipped_ceiling(g, verts, up, tris, held, nowhere, b, 0.866)
    assert (1, 1, 1) in {tuple(x) for x in g.node_of(p).tolist()}
    assert np.allclose(n, b)

    assert not len(_clipped_ceiling(g, verts, up, tris, held, nowhere, b, 1.1)[0])
    on_strut = lambda P: np.ones(len(P), dtype=bool)                    # noqa: E731
    assert not len(_clipped_ceiling(g, verts, up, tris, held, on_strut, b, 0.866)[0])


def test_a_crease_of_faces_under_the_bar_gets_no_clipped_ceiling_point():
    """The vertex normals are averaged across a crease, and two faces each under
    the bar can average to a flat ceiling: the rehearsal part's roof apex of 13
    and 30 deg faces averaged to 8.5 deg from flat, and a ceiling point there was
    reported unbacked. The face's own normal decides."""
    from latticegen2.support.frustum import _clipped_ceiling, _Grid

    lp = lattice_params(10.0, 1.5, (0.0, 0.0, 0.0))
    g = _Grid(lp, 2)
    corner = np.array([1, 1, 1]) * g.q + g.o
    tri = corner + np.array([[0.2, -0.1, -0.05], [-0.1, 0.2, -0.05], [-0.05, -0.1, 0.2]])
    verts = g.world(tri)
    tris = np.array([[0, 1, 2]])
    b = np.asarray(lp.b, dtype=float)
    averaged = np.broadcast_to(b, (3, 3))
    nowhere = lambda P: np.zeros(len(P), dtype=bool)                    # noqa: E731
    held = {tuple(x) for x in g.node_of(verts).tolist()}
    side = np.cross(b, [1.0, 0.0, 0.0])
    side /= np.linalg.norm(side)
    # 30 deg from vertical: n . b = sin 30.
    roof = np.sin(np.radians(30.0)) * b + np.cos(np.radians(30.0)) * side

    def face(n):
        return lambda t, w: (None, np.broadcast_to(n, (len(w), 3)))

    assert not len(_clipped_ceiling(g, verts, averaged, tris, held, nowhere, b, 0.866,
                                    surface=face(roof))[0])
    p, n = _clipped_ceiling(g, verts, averaged, tris, held, nowhere, b, 0.866,
                            surface=face(b))
    assert len(p) and np.allclose(n, b)


def test_a_core_is_the_staircase_less_what_its_trim_did_not_keep():
    """The core the pocket search falls back on: the staircase less the debris
    and dropped pieces inside it; empty where the cut does not take out exactly
    those, rather than a core that is not the support."""
    from latticegen2 import occ
    from latticegen2.support.frustum import _core

    piece = _box((0.0, 0.0, 0.0), (10.0, 10.0, 10.0))
    chip = _box((1.0, 1.0, 1.0), (2.0, 2.0, 2.0))
    assert _core(piece, []) == [piece]
    got = _core(piece, [chip])
    assert sum(occ.volume(s) for s in got) == pytest.approx(999.0, rel=1e-9)
    elsewhere = _box((20.0, 0.0, 0.0), (21.0, 1.0, 1.0))
    assert _core(piece, [elsewhere]) == []


def test_a_point_resting_on_a_strut_face_is_backed_by_it():
    """A ceiling point in a strut a hair below it is backed though the strut's
    face is nearer than PROBE_DEPTH: on the rehearsal part a keel lay in a strut
    to 0.08 mm, and was reported unbacked. One just outside the strut is not."""
    from latticegen2.support.assemble import _LocalLattice
    from latticegen2.support.frustum import PROBE_DEPTH, lattice_backs

    lp = lattice_params(10.0, 1.5, (0.0, 0.0, 0.0))
    local = _LocalLattice(lp, np.zeros((1, 3), dtype=np.int64))
    e0, e1 = (np.asarray(lp.e[k], dtype=float) for k in (0, 1))
    # On the +e0 half-strut, 0.05 mm inside its +e1 face; outward normal -e1, so
    # the probe moves towards the face and leaves the strut at PROBE_DEPTH.
    inside = 0.3 * lp.a * e0 + (lp.t / 2.0 - 0.05) * e1
    outside = 0.3 * lp.a * e0 + (lp.t / 2.0 + 0.05) * e1
    N = np.array([-e1, -e1])
    P = np.array([inside, outside])
    assert not local.contains(P - PROBE_DEPTH * N).any()
    assert lattice_backs(local.contains, P, N).tolist() == [True, False]


def test_a_lost_piece_holding_a_bare_point_joins_the_support_it_touches():
    from latticegen2 import occ
    from latticegen2.support.frustum import exact_inside
    from latticegen2.support.generate import absorb_lost

    support = _box((0.0, 0.0, 0.0), (2.0, 2.0, 2.0))
    tip = _box((1.0, 0.0, 0.0), (3.0, 2.0, 2.0))
    elsewhere = _box((10.0, 0.0, 0.0), (11.0, 1.0, 1.0))
    bare = np.array([[2.5, 1.0, 1.0]])

    solids, n = absorb_lost([support], [elsewhere, tip], bare)

    assert n == 1 and len(solids) == 1
    assert occ.volume(solids[0]) == pytest.approx(12.0)     # overlap counted once
    assert exact_inside(solids[0], bare).all()
