"""Sub-micron corner facets on lattice bodies (docs/algorithm.md §9).

The "curve on curve" error. A strut's corner poking a fraction of a micron
*out* through a curved wall of the input is shaved by the trim, which leaves a
facet of the wall -- a face of ~1e-10 mm^2 lying in the outer wires of the strut
faces that met at the corner. The body is valid and manifold, and ``BRepMesh``
leaves the facet's rim used once at every deflection, so §9's STEP round-trip
check refused the run after the whole pipeline had been paid for.

It is the sub-micron cap of ``test_lattice_micro_caps.py`` turned inside out:
the cap closes an *inner* wire of one wall face; the facet sits between several
faces. `simplify` now collapses it to one vertex.

The fixture is the real failing geometry, per docs/testing.md: junction
(197, -258, -279) of `TDHX_error_report_test.stp` (a user's cut-out, received
2026-11-23 as `TDHX_(5).stp`) at ``cc=7, t=1.4, --orient 180 0 0``, trimmed by
:func:`latticegen2.boundary.trim_junction` and closed into a solid. Its facet
is a 3.568e-10 mm^2 cylinder triangle with sides of 12 to 69 nm at
[1880.289, -74.490, 972.343]; the whole 61,582-face lattice body of that run
read the same 7 edges, and nothing else, at every rung.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from latticegen2 import occ, pipeline
from latticegen2.errors import ProcessingError
from latticegen2.lattice import strut_direction
from latticegen2.pipeline import _simplify_one
from latticegen2.support import frustum, generate

HERE = os.path.dirname(__file__)
FIXTURE = os.path.join(HERE, "corner-facet-junction.brep")
FACET_AT = (1880.289, -74.490, 972.343)


def load_fixture():
    from OCP.BRep import BRep_Builder
    from OCP.BRepTools import BRepTools
    from OCP.TopoDS import TopoDS_Shape

    shape = TopoDS_Shape()
    BRepTools.Read_s(shape, FIXTURE, BRep_Builder())
    (solid,) = occ.solids(shape)
    return solid


@pytest.fixture(scope="module")
def junction():
    return load_fixture()


def _shaved_cylinder(depth: float):
    """A cylinder's crest with a cube's corner poking ``depth`` out of it,
    intersected: the cube comes back with its tip shaved off by the wall.

    The cube runs along the lattice's own edges, so the trihedral is the one a
    junction presents to a wall -- the cap fixture's cube, turned round."""
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Common
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

    R = 20.0
    cyl = BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0), gp_Dir(0, 1, 0)), R, 40.0).Shape()
    e = [np.asarray(strut_direction(k)) * 3.0 for k in range(3)]
    c = np.array([20.3, 0.0, R + depth]) - e[0] - e[1] - e[2]
    base = occ.polygon_face(np.array([c, c + e[0], c + e[0] + e[1], c + e[1]]))
    (solid,) = occ.solids(BRepAlgoAPI_Common(occ.prism(base, e[2]), cyl).Shape())
    return solid


def _displacement(before, after) -> float:
    """How far the boundary moved, in mm, read as ``|dV| / area``.

    A volume difference alone is not a bar here: measured on the shaved cube at
    5e-5 mm, collapsing a facet of ~1e-14 mm^3 shifts the integrated volume by
    1.3e-4 mm^3, which is 2.4e-6 mm spread over the body's surface."""
    return abs(occ.volume(after) - occ.volume(before)) / occ.area(before)


def test_the_junction_is_sound_and_the_round_trip_check_refuses_it(junction, tmp_path):
    """The defect as the user met it: a valid, manifold body that no
    tessellation closes, at any rung of the ladder."""
    assert occ.is_valid(junction)
    assert frustum.edges_not_used_twice(junction) == 0
    facets = [f for f in occ.faces(junction) if occ.area(f) < 1e-8]
    assert len(facets) == 1
    lo, hi = occ.bounding_box(facets[0])
    assert np.linalg.norm(hi - lo) < 1e-4
    assert np.allclose((lo + hi) / 2, FACET_AT, atol=1e-3)

    md = occ.exported_mesh_defects(junction, str(tmp_path / "rt.step"))
    assert md.bad == 7 and md.by_use == {1: 7}
    assert all(np.allclose(p, FACET_AT, atol=1e-3) for p in md.where)
    assert md.refinement is not None and not md.refinement.resolved
    assert [n for _d, n in md.refinement.counts] == [7] * len(occ.MESH_REFINEMENT_LADDER)


def test_the_cap_repair_alone_did_not_see_it(junction):
    """A cap's rim is one inner wire of one host. This facet's rim is in the
    outer wires of three faces, so it is collapsed as a corner, and the record
    says so."""
    hosts: list = []
    corners: list = []
    _fixed, where = occ.remove_micro_caps(junction, hosts=hosts, corners=corners)
    assert hosts == []
    assert len(corners) == 1
    touched, n = corners[0]
    assert n == 1 and len(touched) == 3
    assert len(where) == 1 and np.allclose(where[0], FACET_AT, atol=1e-3)


def test_collapsing_the_facet_lets_the_body_be_written(junction, tmp_path):
    hosts: list = []
    corners: list = []
    fixed, where = occ.remove_micro_caps(junction, hosts=hosts, corners=corners)
    assert occ.micro_cap_repair_fault(junction, fixed, hosts, corners=corners) is None
    (out,) = occ.solids(fixed)
    assert len(occ.faces(out)) == len(occ.faces(junction)) - 1
    assert occ.is_valid(out)
    assert frustum.edges_not_used_twice(out) == 0
    # The facet's own volume is ~1e-14 mm^3. What moves beyond it is the
    # rebuilt faces re-integrated around wires that now close within the
    # merged vertex's tolerance rather than exactly, so it is read as §9 reads
    # unification: as a displacement of the boundary, ``|dV| / area``.
    assert _displacement(junction, out) < 1e-5
    md = occ.exported_mesh_defects(out, str(tmp_path / "rt.step"))
    assert md.bad == 0 and md.refinement is None


def test_simplify_collapses_it_and_reports_where(junction):
    out, ran, degraded, removed, refused = _simplify_one(junction)
    assert not degraded and refused == []
    assert len(removed) == 1 and np.allclose(removed[0], FACET_AT, atol=1e-3)
    assert occ.is_valid(out)
    assert frustum.edges_not_used_twice(out) == 0


def test_support_bodies_get_the_same_repair(junction):
    """Support bodies are trimmed against the same walls and run the same
    repair, checked whole (docs/algorithm.md §14.3 S6)."""
    (out,), removed, refused = generate.remove_micro_caps([junction])
    assert refused == []
    assert len(removed) == 1 and np.allclose(removed[0], FACET_AT, atol=1e-3)
    assert frustum.edges_not_used_twice(out) == 0


def test_the_local_check_names_a_wrong_repair(junction):
    hosts: list = []
    corners: list = []
    fixed, _where = occ.remove_micro_caps(junction, hosts=hosts, corners=corners)
    touched, n = corners[0]
    # One touched face fewer than were rebuilt: the counts no longer agree.
    assert "rebuilt" in occ.micro_cap_repair_fault(
        junction, fixed, hosts, corners=[(touched[:-1], n)])
    # The input itself is not the input less a facet.
    assert occ.micro_cap_repair_fault(junction, junction, hosts, corners=corners) is not None


def test_a_repair_failing_the_local_check_keeps_the_facet(junction, monkeypatch):
    monkeypatch.setattr(occ, "micro_cap_repair_fault", lambda *a, **k: "refused")
    out, ran, degraded, removed, refused = _simplify_one(junction)
    assert removed == []
    assert len(refused) == 1 and np.allclose(refused[0], FACET_AT, atol=1e-3)
    assert len(occ.faces(out)) == len(occ.faces(junction))


@pytest.mark.parametrize("depth", [5e-5, 2e-4])
def test_a_shaved_corner_is_collapsed_exactly(depth):
    """The construction in general, on a corner whose facet is known: three
    faces of the cube and the shaved tip of the wall between them."""
    s = _shaved_cylinder(depth)
    assert len(occ.faces(s)) == 7
    hosts: list = []
    corners: list = []
    fixed, where = occ.remove_micro_caps(s, hosts=hosts, corners=corners)
    assert len(where) == 1 and np.allclose(where[0], (20.3, 0.0, 20.0), atol=1e-3)
    assert occ.micro_cap_repair_fault(s, fixed, hosts, corners=corners) is None
    (out,) = occ.solids(fixed)
    assert len(occ.faces(out)) == 6
    assert occ.is_valid(out)
    assert frustum.edges_not_used_twice(out) == 0
    assert _displacement(s, out) < 1e-5


def test_a_corner_further_out_than_the_bar_is_left_alone():
    """A facet larger than ``MICRO_CAP_TOL`` is real geometry, not noise."""
    s = _shaved_cylinder(0.05)
    fixed, where = occ.remove_micro_caps(s)
    assert where == [] and fixed is s


def test_the_round_trip_failure_names_where_the_body_breaks(junction, tmp_path):
    """The run that found this said where the *input's* seam gap was and not
    where the output broke, and the user went to look at the gap."""
    from latticegen2.runlog import RunLog

    log = RunLog(str(tmp_path / "x.log")).open()
    far_gap = occ.SurfaceSeamGap(4.8745e-3, (1909.791, -60.231, 937.482), 1395, 31)
    with pytest.raises(ProcessingError) as err:
        pipeline._check_export_truth(log, [junction], str(tmp_path), {}, far_gap,
                                     1.4, 7 / np.sqrt(2))
    msg = str(err.value)
    assert "[1880.289, -74.490, 972.343] mm" in msg
    assert "not what this body failed on" in msg
    assert "the input is what has to change" not in msg
