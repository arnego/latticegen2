"""Sub-micron caps on support bodies (docs/algorithm.md §14.3 S6).

A support node's corner within a fraction of a micron of the input's wall
leaves the trimmed support carrying the corner's tip: three faces of ~1e-8
mm^2 closing a 3-edge inner wire of the wall face. The body is valid and
manifold, and no mesher closes the wall around a hole that small, so §9's
export-truth gate refuses the run. The repair removes the cap and its wire,
which is exact by construction; these tests ask that of the solids.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from latticegen2 import occ
from latticegen2.lattice import strut_direction
from latticegen2.support import generate

HERE = os.path.dirname(__file__)
FIXTURE = os.path.join(HERE, "support-micro-cap.brep")


def load_fixture():
    """The support body that failed export truth.

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
    assert generate._edges_not_used_twice(out) == 0
    # The cap's own volume is under depth^3; what moves beyond that is the
    # host face re-integrated, which is quadrature noise.
    assert abs(occ.volume(out) - occ.volume(s)) < 1e-6


def test_a_dimple_that_is_material_is_left_alone():
    s = _dimpled_cylinder(0.5)
    (out,), removed, refused = generate.remove_micro_caps([s])
    assert removed == [] and refused == []
    assert out is s


def test_a_cap_on_an_outer_wire_is_left_alone():
    """On the seam the hole is part of the wall face's outer wire: removing that
    wire would remove the face, so the repair does not apply and the body is
    returned as it came, for export truth to judge."""
    s = _dimpled_cylinder(1e-4, seam_on_crest=True)
    assert len(occ.faces(s)) == 6
    (out,), removed, refused = generate.remove_micro_caps([s])
    assert removed == [] and refused == []
    assert out is s


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
    assert generate._edges_not_used_twice(out) == 0
    assert abs(occ.volume(out) - occ.volume(s)) / occ.area(s) < occ.MICRO_CAP_TOL

    after = occ.exported_mesh_defects(out, str(tmp_path / "after.step"))
    assert after.bad == 0
