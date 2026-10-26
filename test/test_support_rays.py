"""Overhang regions and first-contact rays (docs/algorithm.md §14.3 S1, S3a).

Everything here is checked against an answer known in closed form, or against
an independent method -- a dense fan of plain ray casts for the exact cone
query -- never against the code's own earlier output.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from latticegen2.lattice import lattice_params

TESTDIR = os.path.dirname(os.path.abspath(__file__))
SQRT2 = np.sqrt(2.0)
THETA = np.arcsin(np.sqrt(2.0 / 3.0))


def _meshed(shape, orient=(0.0, 0.0, 0.0), cc=10.0, t=1.5):
    from latticegen2.classify import tessellate_surface

    lp = lattice_params(cc, t, orient)
    tessellate_surface(shape, lp)
    return lp


def _box(lo, size):
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt

    from latticegen2 import occ

    return occ.solids(BRepPrimAPI_MakeBox(gp_Pnt(*lo), *size).Shape())[0]


def _obstacles(body):
    from latticegen2.support import obstacles

    return obstacles.from_shapes(body, [], deflection=0.05)


# --- regions --------------------------------------------------------------------


def test_a_box_lid_is_one_region_refined_to_the_step_with_one_closed_loop():
    from latticegen2.support import regions

    box = _box((0.0, 0.0, 0.0), (30.0, 20.0, 10.0))
    lp = _meshed(box)
    mesh = regions.build(box, lp.b)
    regs = regions.regions(mesh)
    assert len(regs) == 1
    assert regs[0].area == pytest.approx(600.0, rel=1e-12)
    assert np.allclose(mesh.verts[:, 2], 10.0)
    e = np.concatenate([mesh.tris[:, [0, 1]], mesh.tris[:, [1, 2]], mesh.tris[:, [2, 0]]])
    assert np.linalg.norm(mesh.verts[e[:, 0]] - mesh.verts[e[:, 1]], axis=1).max() <= 0.5
    assert len(regs[0].loops) == 1
    loop = mesh.verts[regs[0].loops[0]]
    # The loop is the lid's border: every point on one of its four sides.
    on_side = (np.isclose(loop[:, 0], 0) | np.isclose(loop[:, 0], 30)
               | np.isclose(loop[:, 1], 0) | np.isclose(loop[:, 1], 20))
    assert on_side.all()


def test_the_refined_mesh_is_well_shaped_not_bisected_slivers():
    """BRepMesh triangulates a large plane into slivers; Delaunay flips between
    bisection passes are what keep the refined mesh's triangles round."""
    from latticegen2 import occ
    from latticegen2.support import regions

    body = occ.read_step(os.path.join(TESTDIR, "Hole-block-basic.STEP"))
    _meshed(body, (90.0, 0.0, 0.0))
    roof = max(occ.faces(body), key=occ.area)
    P, _uv, T = regions._refine_face(roof, 0.5)
    a, b, c = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
    L2 = (np.linalg.norm(b - a, axis=1) ** 2 + np.linalg.norm(c - b, axis=1) ** 2
          + np.linalg.norm(a - c, axis=1) ** 2)
    quality = 4.0 * np.sqrt(3.0) * 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1) / L2
    assert np.median(quality) > 0.7
    assert np.sqrt(L2 / 3.0).max() <= 0.5


def test_a_sphere_cap_region_ends_on_the_exact_sixty_degree_isocline():
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeSphere

    from latticegen2.support import regions

    radius = 20.0
    ball = BRepPrimAPI_MakeSphere(radius).Shape()
    lp = _meshed(ball)
    mesh = regions.build(ball, lp.b)
    regs = regions.regions(mesh)
    assert len(regs) == 1
    cap = 2.0 * np.pi * radius ** 2 * (1.0 - np.sin(np.radians(60.0)))
    # The region is whole triangles inside the bar, so it falls short of the
    # cap by at most a band one step wide along the isocline.
    rim = 2.0 * np.pi * radius * np.cos(np.radians(60.0))
    assert cap - rim * 0.5 <= regs[0].area <= cap
    iso = mesh.isocline
    assert len(iso) > 100
    assert np.allclose(np.linalg.norm(iso, axis=1), radius, atol=1e-9)
    assert np.abs(iso[:, 2] - radius * np.sin(np.radians(60.0))).max() < 1e-4


def test_a_tube_underside_is_supported_within_thirty_degrees_of_its_bottom():
    """A horizontal hole: its lower surface is a ceiling of the body below it."""
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

    from latticegen2 import occ
    from latticegen2.support import regions

    r = 4.0
    tube = BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(-1.0, 10.0, 10.0), gp_Dir(1.0, 0.0, 0.0)), r, 32.0).Shape()
    body = occ.solids(BRepAlgoAPI_Cut(_box((0.0, 0.0, 0.0), (30.0, 20.0, 20.0)), tube).Shape())[0]
    lp = _meshed(body)
    mesh = regions.build(body, lp.b)
    regs = [g for g in regions.regions(mesh) if g.area < 500.0]   # not the lid
    assert len(regs) == 1
    # The arc within 30 deg either side of the bottom, along 30 mm of tube.
    exact = 30.0 * r * 2.0 * np.radians(30.0)
    assert exact - 2 * 30.0 * 0.5 <= regs[0].area <= exact
    iso = mesh.isocline[np.abs(mesh.isocline[:, 2] - 10.0) < r + 1e-6]
    assert len(iso) > 50
    # On the circle, at exactly 30 deg from the bottom.
    rel = iso[:, 1:] - [10.0, 10.0]
    assert np.allclose(np.linalg.norm(rel, axis=1), r, atol=1e-9)
    assert np.abs(rel[:, 1] + r * np.cos(np.radians(30.0))).max() < 1e-4


def test_hole_block_chamfered_supports_its_flat_and_fifteen_degree_roof_only():
    """The 30 deg facet overhangs exactly 60 deg from vertical, and the bar is strict."""
    from latticegen2 import occ
    from latticegen2.support import regions

    body = occ.read_step(os.path.join(TESTDIR, "Hole-block-chamfered.STEP"))
    lp = _meshed(body, (90.0, 0.0, 0.0))
    mesh = regions.build(body, lp.b)
    regs = regions.regions(mesh)
    faces = occ.faces(body)
    tilt = {}
    for i, f in enumerate(faces):
        if occ.is_planar(f):
            from latticegen2.support.overhang import _SurfaceNormals

            uv = occ.face_uv_triangulation(f)[1][0]
            _p, n = _SurfaceNormals(f).at(*uv)
            tilt[i] = float(np.degrees(np.arccos(np.clip(n @ lp.b, -1, 1))))
    flat = [i for i, a in tilt.items() if a < 1e-6]
    fifteen = [i for i, a in tilt.items() if abs(a - 15.0) < 1e-6]
    thirty = [i for i, a in tilt.items() if abs(a - 30.0) < 1e-6]
    assert len(flat) == len(fifteen) == len(thirty) == 1
    supported = set(np.unique(mesh.tri_face[mesh.required]).tolist())
    assert supported == {flat[0], fifteen[0]}
    total = sum(g.area for g in regs)
    assert total == pytest.approx(occ.area(faces[flat[0]]) + occ.area(faces[fifteen[0]]),
                                  rel=2e-4)


# --- first-contact rays -----------------------------------------------------------


def test_a_wide_flat_ceiling_over_a_floor_reaches_straight_down():
    """Every point of the floor in the cone ties on drop; the tie goes to the
    nearest, which is straight below."""
    from latticegen2.support import rays

    body = _box((0.0, 0.0, 0.0), (100.0, 100.0, 5.0))
    lp = _meshed(body)
    ob = _obstacles(body)
    apex = np.array([[50.0, 50.0, 5.0 - 1e-3], [40.0, 61.3, 5.0 - 1e-3]])
    got = rays.first_contact(apex, lp.b, ob)
    assert np.allclose(got.drop, 5.0 - 1e-3, atol=1e-12)
    assert np.allclose(got.target, np.column_stack([apex[:, :2], [0.0, 0.0]]), atol=1e-9)
    assert (got.anchor == 0).all()


def test_beside_a_vertical_wall_the_ray_runs_on_the_cone_at_exactly_the_rod_angle():
    from latticegen2.support import rays

    body = _box((0.0, 0.0, 0.0), (100.0, 100.0, 50.0))
    lp = _meshed(body)
    ob = _obstacles(body)
    rho = np.array([0.7, 2.0, 5.5, 9.25])
    apex = np.column_stack([rho, np.full(4, 50.0), np.full(4, 50.0 - 1e-3)])
    got = rays.first_contact(apex, lp.b, ob)
    assert np.allclose(got.drop, rho / SQRT2, atol=1e-12)
    d = got.target - apex
    angle = np.arccos(-d[:, 2] / np.linalg.norm(d, axis=1))
    assert np.allclose(angle, THETA, atol=1e-12)
    assert np.allclose(got.target[:, 0], 0.0, atol=1e-12)


def test_over_a_sloped_floor_the_ray_lands_on_the_high_side():
    """First contact by drop, not by distance: the cone's rim reaches the rising
    floor before its axis does."""
    from OCP.BRepPrimAPI import BRepPrimAPI_MakePrism
    from OCP.gp import gp_Vec

    from latticegen2 import occ
    from latticegen2.support import rays

    s = np.tan(np.radians(20.0))                 # the floor rises 20 deg toward +x
    prof = np.array([[0, 0, 0], [80, 0, 80 * s], [80, 0, 40], [0, 0, 40]], dtype=float)
    body = occ.solids(BRepPrimAPI_MakePrism(
        occ.polygon_face(prof), gp_Vec(0.0, 60.0, 0.0)).Shape())[0]
    lp = _meshed(body)
    ob = _obstacles(body)
    x0, H = 40.0, 40.0 - 1e-3
    got = rays.first_contact(np.array([[x0, 30.0, H]]), lp.b, ob)
    # Along the generator toward +x: t sin(th) * s = (H - x0 s) - t cos(th).
    t = (H - x0 * s) / (np.cos(THETA) + s * np.sin(THETA))
    assert got.drop[0] == pytest.approx(t * np.cos(THETA), abs=1e-9)
    assert got.target[0, 0] > x0
    assert got.drop[0] < H - x0 * s                 # shallower than straight down


def test_the_exact_query_is_never_beaten_by_a_dense_fan_of_rays():
    """The fan samples the cone, so it can only find a contact as deep or
    deeper; it approaches the exact one as it grows."""
    from latticegen2 import occ
    from latticegen2.support import rays

    body = occ.read_step(os.path.join(TESTDIR, "Hole-block-basic.STEP"))
    lp = _meshed(body, (90.0, 0.0, 0.0))
    ob = _obstacles(body)
    rng = np.random.default_rng(7)
    pts = np.column_stack([rng.uniform(-24, 24, 3), np.full(3, 30.0 - 1e-3),
                           rng.uniform(-24, 24, 3)])
    got = rays.first_contact(pts, lp.b, ob)
    for i, x in enumerate(pts):
        fan, _p = rays.brute_force(x, lp.b, ob, polar=120, azimuth=360)
        assert got.drop[i] <= fan + 1e-9
        assert fan - got.drop[i] < 0.02


def test_a_surface_within_the_limit_is_its_own_first_contact():
    """Tilted 45 deg off vertical, a lid is steep enough to print on its own:
    the cone from just inside it meets the lid itself, at once."""
    from latticegen2.support import rays

    body = _box((0.0, 0.0, 0.0), (40.0, 40.0, 40.0))
    lp = _meshed(body, (45.0, 0.0, 0.0))
    ob = _obstacles(body)
    apex = np.array([[20.0, 20.0, 40.0 - 1e-3], [11.0, 27.0, 40.0 - 1e-3]])
    got = rays.first_contact(apex, lp.b, ob)
    assert (got.length < 0.01).all()


def test_the_hole_block_roof_reaches_the_tubes_and_walls_as_the_closed_form_says():
    """Vertical tubes and walls: first contact is the nearest one horizontally,
    at drop rho / tan(theta) = rho / sqrt(2)."""
    from latticegen2 import occ
    from latticegen2.support import rays

    body = occ.read_step(os.path.join(TESTDIR, "Hole-block-basic.STEP"))
    lp = _meshed(body, (90.0, 0.0, 0.0))
    ob = _obstacles(body)
    centres = np.array([[x, z] for x in (-15, -5, 5, 15) for z in (-15, -5, 5, 15)], float)
    g = np.arange(-24.5, 25.0, 1.0)
    pts = np.array([[x, 30.0 - 1e-3, z] for x in g for z in g])
    tube = np.linalg.norm(pts[:, None, [0, 2]] - centres[None], axis=2).min(axis=1) - 1.5
    pts, tube = pts[tube > 0.2], tube[tube > 0.2]
    wall = np.minimum(25 - np.abs(pts[:, 0]), 25 - np.abs(pts[:, 2]))
    got = rays.first_contact(pts, lp.b, ob)
    expect = np.minimum(tube, wall) / SQRT2
    # The tubes are meshed at 0.05 mm deflection, polygons inside the circle,
    # so a contact on one is reached at most that much later.
    assert (got.drop - expect >= -1e-9).all()
    assert (got.drop - expect <= 0.05 / SQRT2 + 1e-9).all()
    # The deepest point is between four tubes: (5 sqrt 2 - 1.5) / sqrt 2.
    assert got.drop.max() <= (5 * SQRT2 - 1.5) / SQRT2 + 0.05
