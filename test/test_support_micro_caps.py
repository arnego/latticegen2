"""Sub-micron caps on support bodies (docs/algorithm.md §14.3 S6).

A support node's corner within a fraction of a micron of the input's wall
leaves the trimmed support carrying the corner's tip: three faces of ~1e-8
mm^2 closing a 3-edge inner wire of the wall face. The body is valid and
manifold, and no mesher closes the wall around a hole that small, so §9's
STEP round-trip check refuses the run. The repair removes the cap and its wire,
which is exact by construction; these tests ask that of the solids.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from latticegen2 import occ
from latticegen2.lattice import strut_direction
from latticegen2.support import frustum, generate

HERE = os.path.dirname(__file__)
FIXTURE = os.path.join(HERE, "support-micro-cap.brep")


def load_fixture():
    """The support body that failed the STEP round-trip check.

    `test-cylinder.STEP` at ``cc=10, t=1.5, --orient 140.5 -76.7 103.6
    --support``, 2026-11-07: the third of three support bodies, 1,047 faces,
    6,442 mm^3, its cap at [171.544, 34.737, 100.905]. Committed whole rather
    than cut down: the tessellation fails only on the whole wall face -- the
    body cut to an 80 mm box around the cap still fails, cut to 50 mm it does
    not -- so a smaller fixture would test a body that never failed.
    """
    from OCP.BRep import BRep_Builder
    from OCP.BRepTools import BRepTools
    from OCP.TopoDS import TopoDS_Shape

    shape = TopoDS_Shape()
    BRepTools.Read_s(shape, FIXTURE, BRep_Builder())
    (solid,) = occ.solids(shape)
    return solid


def _dimpled_cylinder(depth: float, seam_on_crest: bool = False):
    """A cylinder with a support node's corner pressed ``depth`` into its crest.

    The corner is a cube on the lattice's own edges, whose lowest corner points
    straight down -- the trihedral a support node presents to a ceiling."""
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

    R = 20.0
    ax = (gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0)) if seam_on_crest else
          gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0), gp_Dir(0, 1, 0)))
    cyl = BRepPrimAPI_MakeCylinder(ax, R, 40.0).Shape()
    e = [np.asarray(strut_direction(k)) * 3.0 for k in range(3)]
    c = np.array([20.3, 0.0, R - depth])
    base = occ.polygon_face(np.array([c, c + e[0], c + e[0] + e[1], c + e[1]]))
    (solid,) = occ.solids(BRepAlgoAPI_Cut(cyl, occ.prism(base, e[2])).Shape())
    return solid


@pytest.mark.parametrize("depth", [5e-5, 2e-4])
def test_a_micro_dimple_is_filled_exactly(depth):
    s = _dimpled_cylinder(depth)
    assert len(occ.faces(s)) == 6
    (out,), removed, refused = generate.remove_micro_caps([s])
    assert refused == []
    assert len(removed) == 1
    assert np.allclose(removed[0], (20.3, 0.0, 20.0), atol=1e-3)
    # The cylinder's own three faces, and only those: the wall continues.
    assert len(occ.faces(out)) == 3
    assert occ.is_valid(out)
    assert frustum.edges_not_used_twice(out) == 0
    # The cap's own volume is under depth^3; what moves beyond that is the
    # host face re-integrated, which is quadrature noise.
    assert abs(occ.volume(out) - occ.volume(s)) < 1e-6


def test_a_repair_that_fails_its_check_keeps_the_body_as_built(monkeypatch):
    """The repair is checked, not trusted: a result that is not one solid is
    refused, and the body goes on exactly as it came, its cap reported."""
    s = _dimpled_cylinder(1e-4)
    real = occ.remove_micro_caps

    def broken(shape, *a, **k):
        _fixed, where = real(shape, *a, **k)
        return occ.compound([shape, shape]), where

    monkeypatch.setattr(occ, "remove_micro_caps", broken)
    (out,), removed, refused = generate.remove_micro_caps([s])
    assert removed == []
    assert len(refused) == 1
    assert np.allclose(refused[0], (20.3, 0.0, 20.0), atol=1e-3)
    assert out is s


def test_a_dimple_that_is_material_is_left_alone():
    s = _dimpled_cylinder(0.5)
    (out,), removed, refused = generate.remove_micro_caps([s])
    assert removed == [] and refused == []
    assert out is s


def test_a_cap_on_an_outer_wire_is_collapsed_instead():
    """On the seam the hole is part of the wall face's outer wire: removing that
    wire would remove the face, so it is not removed as a cap. It is collapsed
    as a corner facet instead (``test_lattice_micro_corners.py``): the cap's
    faces and edges go, their vertices become one, and the outer wire stays."""
    s = _dimpled_cylinder(1e-4, seam_on_crest=True)
    assert len(occ.faces(s)) == 6
    (out,), removed, refused = generate.remove_micro_caps([s])
    assert refused == []
    assert np.allclose(removed, [(20.3, 0.0, 20.0)], atol=1e-3)
    assert len(occ.faces(out)) == 3
    assert occ.is_valid(out)
    assert frustum.edges_not_used_twice(out) == 0
    assert abs(occ.volume(out) - occ.volume(s)) / occ.area(s) < 1e-4


def test_the_real_body_fails_export_truth_and_passes_once_repaired(tmp_path):
    s = load_fixture()
    assert occ.is_valid(s)
    before = occ.exported_mesh_defects(s, str(tmp_path / "before.step"))
    assert before.bad == 3
    assert before.refinement is not None and not before.refinement.resolved

    (out,), removed, refused = generate.remove_micro_caps([s])
    assert refused == []
    assert np.allclose(removed, [(171.544, 34.737, 100.905)], atol=1e-3)
    assert len(occ.faces(out)) == len(occ.faces(s)) - 3
    assert occ.is_valid(out)
    assert frustum.edges_not_used_twice(out) == 0
    assert abs(occ.volume(out) - occ.volume(s)) / occ.area(s) < occ.MICRO_CAP_TOL

    after = occ.exported_mesh_defects(out, str(tmp_path / "after.step"))
    assert after.bad == 0


LENS = os.path.join(HERE, "support-lens.brep")
LENS_AT = (1865.836, -45.051, 919.752)


def load_lens():
    """A lens on a support body: two descriptions of one edge enclosing a sliver.

    Support body 32 of `TDHX_error_report_test.stp` (a user's part) at ``cc=7,
    t=1.4, --orient 180 0 0 --support``, 2026-11-25, cut to a 2 mm cube around
    the lens: a cylinder face of area -5.2e-06 mm^2 bounded by the input's own
    circle and the boolean's B-spline for the same intersection, 0.025 mm long
    and 3.1e-04 mm apart. The whole body, 8.3 MB, read 4 edges used once at
    every rung there and nowhere else; the cut reads the same 4.
    """
    from OCP.BRep import BRep_Builder
    from OCP.BRepTools import BRepTools
    from OCP.TopoDS import TopoDS_Shape

    shape = TopoDS_Shape()
    BRepTools.Read_s(shape, LENS, BRep_Builder())
    (solid,) = occ.solids(shape)
    return solid


def test_the_lens_fails_the_round_trip_check_at_every_rung(tmp_path):
    s = load_lens()
    assert occ.is_valid(s)
    assert frustum.edges_not_used_twice(s) == 0
    md = occ.exported_mesh_defects(s, str(tmp_path / "rt.step"))
    assert md.bad == 4
    assert [n for _d, n in md.refinement.counts] == [4] * len(occ.MESH_REFINEMENT_LADDER)


def test_the_lens_is_collapsed_onto_one_edge(tmp_path):
    s = load_lens()
    hosts: list = []
    corners: list = []
    fixed, where = occ.remove_micro_caps(s, hosts=hosts, corners=corners)
    assert hosts == [] and len(corners) == 1
    touched, n = corners[0]
    assert n == 1 and len(touched) == 2
    assert np.allclose(where, [LENS_AT], atol=1e-3)
    assert occ.micro_cap_repair_fault(s, fixed, hosts, corners=corners) is None

    (out,), removed, refused = generate.remove_micro_caps([s])
    assert refused == [] and np.allclose(removed, [LENS_AT], atol=1e-3)
    assert len(occ.faces(out)) == len(occ.faces(s)) - 1
    assert occ.is_valid(out)
    assert frustum.edges_not_used_twice(out) == 0
    assert abs(occ.volume(out) - occ.volume(s)) / occ.area(s) < 1e-6
    md = occ.exported_mesh_defects(out, str(tmp_path / "rt.step"))
    assert md.bad == 0 or md.refinement.resolved


def test_two_edged_faces_that_are_not_slivers_are_left_alone():
    """A half-disc is two edges between two vertices too; its edges are far
    apart, so it is a face and not a lens."""
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

    half = BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 5.0, 2.0, np.pi).Shape()
    (solid,) = occ.solids(half)
    fixed, where = occ.remove_micro_caps(solid)
    assert where == [] and fixed is solid
