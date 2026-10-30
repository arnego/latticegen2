"""Ceiling points: the refined surface the support generator backs (docs/algorithm.md §14.3).

Everything here is checked against an answer known in closed form, never
against the code's own earlier output.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from latticegen2.lattice import lattice_params

TESTDIR = os.path.dirname(os.path.abspath(__file__))


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


def _required_area(mesh, faces=None) -> float:
    T = mesh.tris[mesh.required]
    if faces is not None:
        T = T[np.isin(mesh.tri_face[mesh.required], list(faces))]
    a, b, c = mesh.verts[T[:, 0]], mesh.verts[T[:, 1]], mesh.verts[T[:, 2]]
    return float(0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1).sum())


def test_a_box_lid_is_its_ceiling_refined_to_the_step():
    from latticegen2.support import regions

    box = _box((0.0, 0.0, 0.0), (30.0, 20.0, 10.0))
    lp = _meshed(box)
    mesh = regions.build(box, lp.b)
    assert _required_area(mesh) == pytest.approx(600.0, rel=1e-12)
    assert np.allclose(mesh.verts[mesh.vert_required][:, 2], 10.0)
    assert np.allclose(mesh.normals[mesh.vert_required], [0.0, 0.0, 1.0])
    e = np.concatenate([mesh.tris[:, [0, 1]], mesh.tris[:, [1, 2]], mesh.tris[:, [2, 0]]])
    assert np.linalg.norm(mesh.verts[e[:, 0]] - mesh.verts[e[:, 1]], axis=1).max() <= 0.5


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


def test_a_sphere_cap_ends_on_the_exact_sixty_degree_isocline():
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeSphere

    from latticegen2.support import regions

    radius = 20.0
    ball = BRepPrimAPI_MakeSphere(radius).Shape()
    lp = _meshed(ball)
    mesh = regions.build(ball, lp.b)
    cap = 2.0 * np.pi * radius ** 2 * (1.0 - np.sin(np.radians(60.0)))
    # Whole triangles inside the bar, so it falls short of the cap by at most a
    # band one step wide along the isocline.
    rim = 2.0 * np.pi * radius * np.cos(np.radians(60.0))
    assert cap - rim * 0.5 <= _required_area(mesh) <= cap
    iso = mesh.isocline
    assert len(iso) > 100
    assert np.allclose(np.linalg.norm(iso, axis=1), radius, atol=1e-9)
    assert np.abs(iso[:, 2] - radius * np.sin(np.radians(60.0))).max() < 1e-4


def test_a_tube_underside_is_ceiling_within_thirty_degrees_of_its_bottom():
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
    faces = occ.faces(body)
    tube_faces = [i for i, f in enumerate(faces) if not occ.is_planar(f)]
    # The arc within 30 deg either side of the bottom, along 30 mm of tube.
    exact = 30.0 * r * 2.0 * np.radians(30.0)
    assert exact - 2 * 30.0 * 0.5 <= _required_area(mesh, tube_faces) <= exact
    iso = mesh.isocline[np.abs(mesh.isocline[:, 2] - 10.0) < r + 1e-6]
    assert len(iso) > 50
    # On the circle, at exactly 30 deg from the bottom.
    rel = iso[:, 1:] - [10.0, 10.0]
    assert np.allclose(np.linalg.norm(rel, axis=1), r, atol=1e-9)
    assert np.abs(rel[:, 1] + r * np.cos(np.radians(30.0))).max() < 1e-4


def test_hole_block_chamfered_ceiling_is_its_flat_and_fifteen_degree_roof_only():
    """The 30 deg facet overhangs exactly 60 deg from vertical, and the bar is strict."""
    from latticegen2 import occ
    from latticegen2.support import regions
    from latticegen2.support.overhang import _SurfaceNormals

    body = occ.read_step(os.path.join(TESTDIR, "Hole-block-chamfered.STEP"))
    lp = _meshed(body, (90.0, 0.0, 0.0))
    mesh = regions.build(body, lp.b)
    faces = occ.faces(body)
    tilt = {}
    for i, f in enumerate(faces):
        if occ.is_planar(f):
            uv = occ.face_uv_triangulation(f)[1][0]
            _p, n = _SurfaceNormals(f).at(*uv)
            tilt[i] = float(np.degrees(np.arccos(np.clip(n @ lp.b, -1, 1))))
    flat = [i for i, a in tilt.items() if a < 1e-6]
    fifteen = [i for i, a in tilt.items() if abs(a - 15.0) < 1e-6]
    thirty = [i for i, a in tilt.items() if abs(a - 30.0) < 1e-6]
    assert len(flat) == len(fifteen) == len(thirty) == 1
    supported = set(np.unique(mesh.tri_face[mesh.required]).tolist())
    assert supported == {flat[0], fifteen[0]}
    assert _required_area(mesh) == pytest.approx(
        occ.area(faces[flat[0]]) + occ.area(faces[fifteen[0]]), rel=2e-4)
