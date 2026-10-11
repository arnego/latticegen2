"""Print support generation (specification.md §4.6, docs/algorithm.md §14).

Three layers, cheapest first: the convex primitives and their clip (pure
NumPy), overhang detection against shapes whose answer is known in closed form,
and one whole run on a small box that goes through the real kernel and the real
pipeline.
"""

from __future__ import annotations

import os
import subprocess
import sys

import numpy as np
import pytest

from latticegen2.lattice import lattice_params, nodes
from latticegen2.support import polyhedra as ph

FRAME = (np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0]))
B = FRAME[2]


def _angle_from_vertical(v: np.ndarray) -> float:
    return float(np.degrees(np.arccos(abs(float(v @ B)) / np.linalg.norm(v))))


# --- primitives -------------------------------------------------------------


@pytest.mark.parametrize("overhang", [30.0, 45.0, 60.0, 75.0])
def test_capital_ridges_sit_at_the_limit_and_faces_inside_it(overhang):
    """Neither a face nor an edge of a support may exceed the limit."""
    tau = np.tan(np.radians(overhang))
    cap = ph.capital(np.array([1.0, 2.0, 3.0]), 5.0, tau, 0.37, FRAME, "node")
    apex = cap.verts[0]
    for corner in cap.verts[1:]:
        assert _angle_from_vertical(corner - apex) == pytest.approx(overhang, abs=1e-9)
    normals, _ = ph.face_planes(cap.verts, cap.faces)
    side = normals[normals @ B < 0.0]          # the eight down-and-out faces
    assert len(side) == ph.SIDES
    # A face's tilt from vertical is the angle its normal makes with horizontal.
    tilt = np.degrees(np.arcsin(-(side @ B)))
    expected = np.degrees(np.arctan(tau * ph.INSCRIBED))
    assert np.allclose(tilt, expected, atol=1e-9)
    assert expected < overhang


def test_post_is_never_thinner_than_the_minimum_thickness():
    thickness = 1.3
    p = ph.post(np.zeros(3), 9.0, thickness, np.tan(np.radians(60.0)), 0.2, FRAME, "node")
    normals, offsets = ph.face_planes(p.verts, p.faces)
    walls = np.abs(normals @ B) < 1e-9             # the eight vertical faces
    assert int(walls.sum()) == ph.SIDES
    # Opposite walls are parallel planes `thickness` apart.
    assert np.allclose(offsets[walls], thickness / 2.0, atol=1e-9)
    # The tip is a cone at the limit, so the post starts from a point.
    for corner in p.verts[1:1 + ph.SIDES]:
        assert _angle_from_vertical(corner - p.verts[0]) == pytest.approx(60.0, abs=1e-9)


@pytest.mark.parametrize("angles", [(0.0, 0.0, 0.0), (20.0, 30.0, 40.0)])
def test_clipping_to_lattice_cells_partitions_the_primitive_exactly(angles):
    """The cells' pieces add up to the primitive: nothing lost, nothing doubled."""
    lp = lattice_params(10.0, 1.5, angles)
    cap = ph.capital(np.array([0.3, 0.2, -1.0]), 9.0, np.tan(np.radians(60.0)),
                     0.3, FRAME, "node")
    whole = ph.volume(cap.verts, cap.faces)
    total = 0.0
    for i in range(-6, 7):
        for j in range(-6, 7):
            for k in range(-6, 7):
                got = ph.clip_to_cell(cap.verts, cap.faces, lp,
                                      nodes(lp, np.array([[i, j, k]]))[0])
                if got is not None:
                    total += ph.volume(*got)
    assert total == pytest.approx(whole, rel=1e-10)


def test_a_piece_presents_no_cross_section_where_its_neighbour_keeps_nothing():
    """A capital poking a hair through a cell face: the far cell drops its
    share, so the near cell must not present a cross-section in that plane."""
    lp = lattice_params(10.0, 1.5)
    pos = nodes(lp, np.array([[0, 0, 0]]))[0]
    half, e = lp.a / 2.0, lp.e[0]
    hx = np.cross(B, e); hx /= np.linalg.norm(hx)
    frame = (hx, np.cross(B, hx), B)
    # Apex inside the cell, the top ring's furthest corner 2e-4 mm past the face.
    tau, height = np.tan(np.radians(60.0)), 1.0
    prim = ph.capital(pos, height, tau, 0.0, frame, "test")
    reach = float((prim.verts @ e).max() - pos @ e)
    prim = ph.capital(pos + (half + 2e-4 - reach) * e, height, tau, 0.0, frame, "test")
    plain = ph.clip_to_cell(prim.verts, prim.faces, lp, pos, 1e-9)
    beyond = ph.clip_to_cell(prim.verts, prim.faces, lp, pos + lp.a * e, 1e-9)
    assert beyond is None                                     # too small to build
    on_face = np.abs(plain[0] @ e - (pos @ e + half)) <= 1e-9
    assert on_face.sum() >= 3                                 # yet a real cross-section

    matched = ph.clip_to_cell_matched(prim.verts, prim.faces, lp, pos, 1e-9)
    assert (matched[0] @ e).max() == pytest.approx(pos @ e + half - ph.RETREAT, abs=1e-12)
    assert ph.volume(*matched) == pytest.approx(ph.volume(*plain), rel=1e-6)

    # Where the neighbour keeps its share, nothing changes at all.
    deep = ph.capital(pos + (half - 0.5 * reach) * e, height, tau, 0.0, frame, "test")
    a = ph.clip_to_cell(deep.verts, deep.faces, lp, pos, 1e-9)
    b = ph.clip_to_cell_matched(deep.verts, deep.faces, lp, pos, 1e-9)
    assert np.array_equal(a[0], b[0]) and a[1] == b[1]


def test_a_piece_inside_another_is_dropped_and_the_rest_kept():
    tau = np.tan(np.radians(60.0))
    big = ph.capital(np.zeros(3), 10.0, tau, 0.1, FRAME, "node")
    small = ph.capital(np.array([0.0, 0.0, 6.0]), 2.0, tau, 0.4, FRAME, "node")
    apart = ph.capital(np.array([40.0, 0.0, 0.0]), 3.0, tau, 0.2, FRAME, "node")
    kept = ph.drop_contained([(p.verts, p.faces) for p in (small, big, apart)])
    assert len(kept) == 2
    assert sorted(len(v) for v, _f in kept) == [9, 9]
    assert not any(np.array_equal(v, small.verts) for v, _f in kept)


# --- placement (kernel-free given a mesh) ------------------------------------


def _box_mesh(lo, hi):
    from latticegen2.classify import TriMesh

    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    c = np.array([[x, y, z] for x in (lo[0], hi[0]) for y in (lo[1], hi[1])
                  for z in (lo[2], hi[2])])
    quads = [(0, 1, 3, 2), (4, 6, 7, 5), (0, 4, 5, 1), (2, 3, 7, 6),
             (0, 2, 6, 4), (1, 5, 7, 3)]
    tris = [t for q in quads for t in ((q[0], q[1], q[2]), (q[0], q[2], q[3]))]
    return TriMesh(verts=c, tris=np.array(tris), deviation=0.0)


def test_ray_caster_finds_the_nearest_wall_and_its_outward_normal():
    from latticegen2.support.plan import RayCaster

    mesh = _box_mesh((0, 0, 0), (10, 10, 10))
    t, n = RayCaster(mesh, np.array([0.0, 0.0, 1.0])).first_hit(
        np.array([[5.0, 5.0, 2.0], [3.0, 7.0, 9.5], [50.0, 5.0, 5.0]]))
    assert t[0] == pytest.approx(8.0)
    assert t[1] == pytest.approx(0.5)
    assert np.isinf(t[2])
    assert np.allclose(n[:2], [0.0, 0.0, 1.0])


def _flat_ceiling_plan(angles=(0.0, 0.0, 0.0)):
    """Anchors on a full grid of interior nodes under a flat lid."""
    from latticegen2.support import plan

    lp = lattice_params(10.0, 1.5, angles)
    mesh = _box_mesh((-40, -40, -40), (40, 40, 12.3))
    params = plan.SupportParams(overhang=60.0, thickness=1.0, pitch=0.5, deviation=0.0)
    xs = np.arange(-14.0, 14.01, 0.5)
    gx, gy = np.meshgrid(xs, xs, indexing="ij")
    points = np.column_stack([gx.ravel(), gy.ravel(), np.full(gx.size, 12.3)])
    idx = np.array([[i, j, k] for i in range(-9, 10) for j in range(-9, 10)
                    for k in range(-9, 10)])
    pos = nodes(lp, idx)
    inside = np.all(np.abs(pos[:, :2]) < 34.0, axis=1) & (pos[:, 2] < 6.0) & (pos[:, 2] > -39.5)
    interior = idx[inside]
    anchors = plan.find_anchors(lp, mesh, params, interior, [], set())
    return lp, params, points, plan.plan_supports(lp, mesh, params, points, anchors)


def test_every_sample_is_covered_with_margin_by_a_planned_capital():
    """Sample coverage with the margin is what makes it continuum coverage."""
    lp, params, points, placed = _flat_ceiling_plan()
    assert placed.stats["support_drop_columns"] == 0
    covered = np.zeros(len(points), dtype=bool)
    ring = [np.array([np.cos(a), np.sin(a), 0.0]) * params.margin
            for a in np.linspace(0, 2 * np.pi, 12, endpoint=False)]
    for prim in placed.primitives:
        if prim.kind != "capital":
            continue
        normals, offsets = ph.face_planes(prim.verts, prim.faces)
        ok = np.ones(len(points), dtype=bool)
        for shift in ring:
            # Just beneath the surface, a margin's radius in every direction.
            q = points + shift - np.array([0.0, 0.0, 1e-6])
            ok &= ((q @ normals.T) <= offsets[None, :] + 1e-9).all(axis=1)
        covered |= ok
    assert covered.all()


def test_planned_supports_respect_thickness_and_are_deterministic():
    lp, params, _points, placed = _flat_ceiling_plan()
    _lp, _pa, _pt, again = _flat_ceiling_plan()
    assert len(placed.primitives) == len(again.primitives) > 0
    for a, b in zip(placed.primitives, again.primitives):
        assert a.kind == b.kind and np.array_equal(a.verts, b.verts)
    for prim in placed.primitives:
        if prim.kind == "capital":
            top = prim.verts[1:]
            across = 2.0 * np.linalg.norm(top[0] - top.mean(axis=0)) * ph.INSCRIBED
            assert across >= params.thickness - 1e-9
        # No vertex sits on a lattice cell plane (the generic-position rule).
        cell = np.linalg.solve(lp.B, prim.verts.T).T
        assert np.abs((cell - 0.5) - np.round(cell - 0.5)).min() > 0.0


def test_no_two_support_faces_share_an_azimuth_with_the_lattice():
    from latticegen2.support.plan import LATTICE_AZIMUTH_CLEARANCE, _phase

    sector = 2.0 * np.pi / ph.SIDES
    phases = [_phase(i) for i in range(200)]
    assert len({round(p, 9) for p in phases}) > 150       # not one shared phase
    for phase in phases:
        normals = phase + sector / 2.0 + np.arange(ph.SIDES) * sector
        off = np.abs((normals + np.pi / 6.0) % (np.pi / 3.0) - np.pi / 6.0)
        assert off.min() >= LATTICE_AZIMUTH_CLEARANCE - 1e-12


# --- overhang detection (kernel) ---------------------------------------------


def _meshed(shape, cc=10.0, t=1.5):
    from latticegen2.classify import tessellate_surface

    lp = lattice_params(cc, t)
    tessellate_surface(shape, lp)
    return lp


def test_box_overhang_is_its_lid_and_nothing_else():
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox

    from latticegen2.support import overhang

    box = BRepPrimAPI_MakeBox(30.0, 20.0, 10.0).Shape()
    lp = _meshed(box)
    got = overhang.detect(box, lp.b, 60.0, 0.5)
    assert got.area == pytest.approx(600.0, rel=1e-9)
    assert np.allclose(got.points[:, 2], 10.0)            # never the underside
    assert np.allclose(got.inward, [0.0, 0.0, -1.0])
    assert got.counts()["edge"] == got.counts()["vertex"] == 0


@pytest.mark.parametrize("overhang_deg", [45.0, 60.0])
def test_sphere_overhang_matches_the_closed_form_cap(overhang_deg):
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeSphere

    from latticegen2.support import overhang

    radius = 30.0
    ball = BRepPrimAPI_MakeSphere(radius).Shape()
    lp = _meshed(ball)
    got = overhang.detect(ball, lp.b, overhang_deg, 0.5, margin_deg=0.0)
    # n.b > sin(phi) is the cap within (90 - phi) of the pole.
    exact = 2.0 * np.pi * radius ** 2 * (1.0 - np.sin(np.radians(overhang_deg)))
    assert got.area == pytest.approx(exact, rel=0.03)
    # Every sample really is on the surface and over the bar: exact normals.
    assert np.allclose(np.linalg.norm(got.points, axis=1), radius, atol=1e-9)
    assert (got.points[:, 2] / radius >= np.sin(np.radians(overhang_deg)) - 1e-9).all()


def test_a_keel_edge_is_an_overhang_even_though_its_faces_are_not():
    """A V-notch in the lid: both faces steep, the ridge between them hanging."""
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox

    from latticegen2 import occ
    from latticegen2.support import overhang

    profile = np.array([[0, 0, 0], [20, 0, 0], [20, 0, 10], [12, 0, 10],
                        [10, 0, 6], [8, 0, 10], [0, 0, 10]], dtype=float)
    notched = occ.prism(occ.polygon_face(profile), np.array([0.0, 15.0, 0.0]))
    assert occ.volume(notched) == pytest.approx((200.0 - 8.0) * 15.0)
    lp = _meshed(notched)
    got = overhang.detect(notched, lp.b, 60.0, 0.5)
    keel = got.points[got.kind == overhang.EDGE]
    assert len(keel) > 10
    assert np.allclose(keel[:, 0], 10.0) and np.allclose(keel[:, 2], 6.0)
    assert keel[:, 1].min() == pytest.approx(0.0) and keel[:, 1].max() == pytest.approx(15.0)

    plain = BRepPrimAPI_MakeBox(20.0, 15.0, 10.0).Shape()
    lp = _meshed(plain)
    assert overhang.detect(plain, lp.b, 60.0, 0.5).counts()["edge"] == 0


def test_a_hanging_tip_is_an_overhang_vertex():
    """A steep cone of shell material hanging into the volume: only its tip overhangs."""
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCone
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

    from latticegen2 import occ
    from latticegen2.support import overhang

    box = BRepPrimAPI_MakeBox(20.0, 20.0, 12.0).Shape()
    cone = BRepPrimAPI_MakeCone(
        gp_Ax2(gp_Pnt(10.0, 10.0, 12.0), gp_Dir(0.0, 0.0, -1.0)), 3.0, 0.0, 6.0
    ).Shape()
    body = occ.solids(BRepAlgoAPI_Cut(box, cone).Shape())[0]
    lp = _meshed(body)
    got = overhang.detect(body, lp.b, 60.0, 0.5)
    tips = got.points[got.kind == overhang.VERTEX]
    assert len(tips) == 1
    assert np.allclose(tips[0], [10.0, 10.0, 6.0], atol=1e-6)


def test_a_surface_at_the_limit_is_sampled_but_not_required():
    """A 120 degree drill point in the lid is tilted exactly 60 degrees.

    It is inside the margin, so supports are placed with it in mind, and it is
    not over the limit, so nothing has to be proven about it.
    """
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCone
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

    from latticegen2 import occ
    from latticegen2.support import overhang

    box = BRepPrimAPI_MakeBox(20.0, 20.0, 10.0).Shape()
    radius = 4.0
    cone = BRepPrimAPI_MakeCone(
        gp_Ax2(gp_Pnt(10.0, 10.0, 10.0), gp_Dir(0.0, 0.0, 1.0)),
        radius, 0.0, radius / np.tan(np.radians(60.0))).Shape()
    body = occ.solids(BRepAlgoAPI_Fuse(box, cone).Shape())[0]
    lp = _meshed(body)
    got = overhang.detect(body, lp.b, 60.0, 0.5)
    face = got.kind == overhang.FACE
    on_cone = face & (got.points[:, 2] > 10.0 + 1e-6)
    on_lid = face & (np.abs(got.points[:, 2] - 10.0) <= 1e-6)
    assert on_cone.sum() > 20 and on_lid.sum() > 100
    assert not got.required[on_cone].any()
    assert got.required[on_lid].all()
    # One degree stricter and the same cone is over the limit.
    assert overhang.detect(body, lp.b, 59.0, 0.5).required[on_cone].all()


# --- the trim that checks itself ----------------------------------------------


def _cube_poly(lo, hi):
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    verts = np.array([[x, y, z] for x in (lo[0], hi[0]) for y in (lo[1], hi[1])
                      for z in (lo[2], hi[2])])
    faces = [[0, 1, 3, 2], [4, 6, 7, 5], [0, 4, 5, 1], [2, 3, 7, 6],
             [0, 2, 6, 4], [1, 5, 7, 3]]
    return verts, ph.orient_outward(verts, faces)


def test_a_trim_whose_two_answers_add_up_is_consistent():
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox

    from latticegen2 import occ
    from latticegen2.support import build

    body = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
    for lo, hi, want in (((2, 2, 2), (4, 4, 4), 8.0),          # inside
                         ((8, 2, 2), (12, 4, 4), 8.0),         # straddling
                         ((11, 2, 2), (13, 4, 4), 0.0)):       # outside
        clock = {}
        solid = build.polyhedron_solid(*_cube_poly(lo, hi))
        inside, ok = build.checked_trim(solid, body, "test", clock)
        assert ok
        assert sum(occ.volume(x) for x in inside) == pytest.approx(want, abs=1e-9)
        assert clock["trim"] > 0.0


def test_an_inconsistent_trim_is_bisected_and_loses_nothing(monkeypatch):
    """The kernel cannot be made to fail on demand, so the first answer is
    declared inconsistent and the halves have to reproduce it."""
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox

    from latticegen2 import occ
    from latticegen2.support import build

    body = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
    poly = _cube_poly((7, 1, 1), (13, 5, 4))
    real, calls = build.checked_trim, []

    def lying(shape, body_, where, clock):
        inside, ok = real(shape, body_, where, clock)
        calls.append(occ.volume(shape))
        return ([], False) if len(calls) == 1 else (inside, ok)

    monkeypatch.setattr(build, "checked_trim", lying)
    # The other route to the same region is taken out of the way here; it has
    # its own test below.
    monkeypatch.setattr(build, "trim_by_planes", lambda *a: ([], False))
    clock = {}
    got = build.trim_convex(poly, body, "test", clock)
    assert clock["bisections"] == 1 and len(calls) == 3
    assert calls[1] + calls[2] == pytest.approx(calls[0], rel=1e-9)
    assert sum(occ.volume(x) for x in got) == pytest.approx(3.0 * 4.0 * 3.0, rel=1e-9)

    # And a leaf that never agrees stops at the depth bound instead of forever.
    monkeypatch.setattr(build, "checked_trim", lambda *a: ([], False))
    clock = {}
    assert build.trim_convex(poly, body, "test", clock) == []
    assert clock["bisections"] == 2 ** build.BISECT_DEPTH - 1


@pytest.mark.parametrize("up", [(0.0, 0.0, 1.0), (0.0, 0.0, -1.0), (0.3, -0.5, 0.81)])
def test_a_bisection_plane_is_vertical_whatever_the_shape(up):
    """An exposed cut has to be self-supporting, so it contains the build direction."""
    from latticegen2.support import build

    up = np.array(up) / np.linalg.norm(up)
    for lo, hi in (((0, 0, 0), (1, 1, 9)), ((0, 0, 0), (9, 1, 1)), ((0, 0, 0), (2, 5, 3))):
        poly = _cube_poly(lo, hi)
        halves = build._halves(*poly, up=up)
        assert len(halves) == 2
        assert sum(ph.volume(*h) for h in halves) == pytest.approx(ph.volume(*poly))
        normals, offsets = ph.face_planes(*halves[0])
        old, _ = ph.face_planes(*poly)
        cut = [n for n in normals if not any(np.allclose(n, m) for m in old)]
        assert len(cut) == 1 and abs(float(cut[0] @ up)) < 1e-12


def test_plane_by_plane_gives_the_same_region_as_the_direct_trim():
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder

    from latticegen2 import occ
    from latticegen2.support import build

    body = BRepPrimAPI_MakeCylinder(10.0, 10.0).Shape()
    # A capital leaning out through the curved wall.
    prim = ph.capital(np.array([7.0, 0.5, 2.0]), 5.0, np.tan(np.radians(60.0)), 0.1,
                      FRAME, "test")
    poly = (prim.verts, prim.faces)
    direct, ok = build.checked_trim(build.polyhedron_solid(*poly), body, "test", {})
    assert ok
    planed, ok = build.trim_by_planes(poly, body, "test", {})
    assert ok
    want = sum(occ.volume(x) for x in direct)
    assert 0.0 < want < ph.volume(*poly)
    assert sum(occ.volume(x) for x in planed) == pytest.approx(want, rel=1e-6)


def test_the_judge_sees_what_the_kernels_own_cross_check_cannot():
    """An operand returned whole, and one returned empty: both add up."""
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox

    from latticegen2.classify import OutsideProbe, tessellate_surface
    from latticegen2.support import build

    body = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
    mesh = tessellate_surface(body, lattice_params(10.0, 1.5))
    judge = build.Judge(OutsideProbe(mesh, 0.05))
    poly = _cube_poly((7, 3, 3), (13, 6, 6))                 # half in, half out
    whole = [build.polyhedron_solid(*poly)]
    right, ok = build.checked_trim(whole[0], body, "test", {})
    assert ok
    assert judge.outside(right) == 0.0 and judge.missing([poly], right) == 0
    assert not judge.refutes([poly], right)
    assert judge.outside(whole) == pytest.approx(3.0)        # untrimmed
    assert judge.missing([poly], []) > 0                     # lost
    assert judge.refutes([poly], whole) and judge.refutes([poly], [])

    # And trim_convex acts on it: a first answer of "nothing" is not accepted.
    real, calls = build.checked_trim, []

    def lying(shape, body_, where, clock):
        calls.append(1)
        return ([], True) if len(calls) == 1 else real(shape, body_, where, clock)

    import latticegen2.support.build as b
    saved, b.checked_trim = b.checked_trim, lying
    try:
        clock = {}
        got = build.trim_convex(poly, body, "test", clock, judge=judge)
    finally:
        b.checked_trim = saved
    from latticegen2 import occ
    assert clock["judge_refusals"] == 1 and clock["by_planes"] == 1
    assert sum(occ.volume(x) for x in got) == pytest.approx(27.0, rel=1e-9)


def test_a_repair_fuse_that_loses_volume_is_folded_instead(monkeypatch):
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt

    from latticegen2 import boundary, occ

    boxes = [BRepPrimAPI_MakeBox(gp_Pnt(float(i), 0.0, 0.0), 1.0, 1.0, 1.0).Shape()
             for i in range(4)]
    got = boundary._fuse_solids(boxes, "test")
    assert len(got) == 1 and occ.volume(got[0]) == pytest.approx(4.0)

    # The one-shot call returns one of its operands; the fold must take over.
    real, calls = boundary.occ.solids, []

    def lying(shape):
        calls.append(1)
        return [boxes[0]] if len(calls) == 1 else real(shape)

    monkeypatch.setattr(boundary.occ, "solids", lying)
    got = boundary._fuse_solids(boxes, "test")
    assert len(calls) == 4                                   # one shot + three steps
    assert len(got) == 1 and occ.volume(got[0]) == pytest.approx(4.0)


def test_a_fragment_beside_the_main_body_is_kept_as_its_own_shell():
    """Two closed surfaces in one component: the junction graph over-connects
    supported cells, so the second is a body and not a contradiction."""
    from OCP.BRep import BRep_Builder
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.TopoDS import TopoDS_Shell
    from OCP.gp import gp_Pnt

    from latticegen2 import occ
    from latticegen2.support.build import fill_voids

    big = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
    crumb = BRepPrimAPI_MakeBox(gp_Pnt(20.0, 0.0, 0.0), 0.5, 0.5, 0.5).Shape()
    shell = TopoDS_Shell()
    builder = BRep_Builder()
    builder.MakeShell(shell)
    for face in occ.faces(big) + occ.faces(crumb):
        builder.Add(shell, face)
    out, n_filled, filled, extra = fill_voids({0: shell, 1: occ.shells(big)[0]})
    assert n_filled == 0 and filled == 0.0
    assert extra == [pytest.approx(0.125)]
    assert sorted(out) == [0, 1, 2]                       # sortable integer keys
    assert len(occ.faces(out[0])) == 6 and len(occ.faces(out[2])) == 6
    assert out[1].IsSame(occ.shells(big)[0])              # an ordinary shell is untouched


# --- the whole pipeline -------------------------------------------------------


@pytest.fixture(scope="module")
def box_step(tmp_path_factory):
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt

    from latticegen2 import occ

    path = tmp_path_factory.mktemp("support") / "box.step"
    occ.write_step(BRepPrimAPI_MakeBox(gp_Pnt(1.3, 0.7, 0.4), 26.0, 23.0, 19.0).Shape(),
                   str(path), "box")
    return str(path)


def _run(args):
    main = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "src", "main.py")
    return subprocess.run([sys.executable, main, *args], capture_output=True, text=True)


def test_a_supported_box_is_valid_covered_and_only_heavier(box_step, tmp_path):
    """End to end through the real kernel: the lid of a box, fully supported."""
    from latticegen2 import occ

    plain = str(tmp_path / "plain.step")
    held = str(tmp_path / "held.step")
    a = _run(["-i", box_step, "-cc", "10", "-t", "1.5", "-o", plain, "--cores", "1"])
    b = _run(["-i", box_step, "-cc", "10", "-t", "1.5", "-o", held, "--cores", "1",
              "--support"])
    assert a.returncode == 0, a.stderr
    assert b.returncode == 0, b.stderr

    log = open(held[:-5] + ".log", encoding="utf-8").read()
    assert "stage overhang" in log and "stage plan" in log and "stage support" in log
    assert "stage overhang" not in open(plain[:-5] + ".log", encoding="utf-8").read()
    samples = int(log.split("support_samples: ")[1].split()[0])
    assert samples > 1000
    assert f"support_samples_verified: {samples}" in log

    lattice = occ.solids(occ.read_step(plain))
    supported = occ.solids(occ.read_step(held))
    assert all(occ.is_valid(s) for s in supported)
    v_lattice = sum(occ.volume(s) for s in lattice)
    v_supported = sum(occ.volume(s) for s in supported)
    assert v_supported > v_lattice * 1.01
    # Nothing outside the body: the lid is the highest the output may reach.
    lo, hi = occ.bounding_box(occ.compound(supported))
    assert hi[2] <= 19.4 + 1e-6 and lo[2] >= 0.4 - 1e-6


def test_support_code_is_not_loaded_without_the_flag(box_step, tmp_path):
    """specification.md §4.6: a run without `--support` never imports it."""
    src = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src")
    out = str(tmp_path / "o.step")
    code = (
        "import sys; sys.path.insert(0, r'%s'); "
        "from latticegen2.__main__ import main; "
        "rc = main(['-i', r'%s', '-cc', '10', '-t', '1.5', '-o', r'%s', "
        "'--cores', '1', '--orient', '10', '0', '0']); "
        "print('LOADED' if any(m.startswith('latticegen2.support') for m in sys.modules) "
        "else 'ABSENT', rc)" % (src, box_step, out)
    )
    proc = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert proc.stdout.strip().endswith("ABSENT 0"), proc.stdout + proc.stderr
