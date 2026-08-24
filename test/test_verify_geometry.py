"""The E2E harness's own containment check (docs/testing.md, specification.md §6.2).

`tools/` is normally tested by running it — `tools/e2e.py` *is* the gate for the
harness. This file exists for the one part of it whose failure mode is a
plausible number rather than an exception, and which therefore cannot be caught
by reading the output of a passing run.

`material_outside` reported **354,733 mm³** of material outside the input body —
the entire lattice — on the `cc=12, t=2.5` output of `TD_HX_rehearsal_test`,
whose every junction was *built* by intersecting with that body. The boolean did
not decline: it reported `IsDone`, `HasModified` and `HasGenerated`, and returned
43,672 faces where 43,530 went in. A lattice trimmed from a body has a large
share of its faces lying exactly *on* that body's surface, which is the classic
ill-conditioned case for a cut.

What these tests pin is the shape of the repair rather than that one part: a
trustworthy answer is *exactly* zero, so any remainder is contradicted against a
boolean-free containment check, and a disagreement is reported as unmeasured
rather than as a pass. The real part is too large to carry in a unit test — it is
covered by `tools/e2e.py` and by the measurements in specification.md §11 — so
these use boxes, where "contained" and "sticking out" are known by construction.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tools"))

import verify_geometry as vg  # noqa: E402

from latticegen2 import occ  # noqa: E402


def box(dx: float, dy: float, dz: float, at=(0.0, 0.0, 0.0)):
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt

    return BRepPrimAPI_MakeBox(gp_Pnt(*at), dx, dy, dz).Solid()


@pytest.fixture(scope="module")
def outer():
    return box(100.0, 100.0, 100.0)


def test_a_contained_solid_is_measured_exactly_and_at_zero(outer, tmp_path):
    """The ordinary case must stay exact — the repair must not have turned every
    measurement into the weaker fallback."""
    inner = box(10.0, 10.0, 10.0, at=(20.0, 20.0, 20.0))
    cand, body = tmp_path / "cand.step", tmp_path / "body.step"
    occ.write_step(inner, str(cand), "cand")
    occ.write_step(outer, str(body), "body")

    result = vg.material_outside(str(cand), str(body))

    assert result["exact"] is True
    assert result["unmeasured"] == []
    assert result["volume_mm3"] == 0.0
    assert all(s["trusted"] for s in result["per_solid"])


def test_a_solid_sticking_out_is_reported_not_swallowed(outer, tmp_path):
    """The contradiction path must never be able to hide real material outside.

    A boolean-free containment check sees a protrusion too, so whichever route
    the result takes, the scenario fails rather than passes.
    """
    poking = box(10.0, 10.0, 10.0, at=(95.0, 20.0, 20.0))   # 5 mm of 10 mm is out
    cand, body = tmp_path / "cand.step", tmp_path / "body.step"
    occ.write_step(poking, str(cand), "cand")
    occ.write_step(outer, str(body), "body")

    result = vg.material_outside(str(cand), str(body))

    if result["exact"]:
        assert result["volume_mm3"] == pytest.approx(500.0, rel=1e-6)
    else:
        # Taken as unmeasured, the containment evidence must still condemn it.
        assert len(result["unmeasured"]) == 1
        assert result["unmeasured"][0]["containment"]["contained"] is False


def test_containment_is_boolean_free_and_agrees_with_construction(outer):
    """`surface_points_outside` is the arbiter, so it has to be right on its own."""
    inside = vg.surface_points_outside(box(10.0, 10.0, 10.0, at=(20.0, 20.0, 20.0)), outer)
    assert inside["sampled"] > 0
    assert inside["outside"] == 0
    assert inside["contained"] is True
    assert inside["cleared_at_mm"] == vg.CONTAINMENT_SWEEP[0], (
        "a solid well clear of the boundary must clear at the tightest tolerance"
    )

    out = vg.surface_points_outside(box(10.0, 10.0, 10.0, at=(95.0, 20.0, 20.0)), outer)
    assert out["outside"] > 0
    assert out["contained"] is False
    assert out["cleared_at_mm"] is None


def test_a_face_lying_on_the_boundary_is_a_tie_not_a_protrusion(outer):
    """The case the whole tolerance sweep exists for.

    A trimmed lattice's faces lie exactly *on* the input surface, so its points
    there are ties. They must not read as material outside — otherwise the
    fallback would condemn every part it was built to rescue.
    """
    flush = box(10.0, 10.0, 10.0, at=(90.0, 20.0, 20.0))    # right face coincident
    ev = vg.surface_points_outside(flush, outer)

    assert ev["contained"] is True, ev
    assert ev["cleared_at_mm"] <= vg.CONTAINMENT_TOL


def test_a_real_protrusion_is_not_swallowed_by_the_distance_contradiction(outer):
    """The load-bearing half: the contradiction must not become an off switch.

    `surface_points_outside` now contradicts a point the tolerance sweep cannot
    clear by measuring its exact distance to the body, because the sweep tops
    out at 2e-03 mm and a 1e-09 mm tie reads there exactly like a 1 mm
    protrusion. A contradiction that only ever exonerates would be worse than
    no check at all (docs/testing.md, G10), so this pins the other direction:
    a box genuinely sticking out measures a real distance and is still
    condemned.
    """
    out = vg.surface_points_outside(box(10.0, 10.0, 10.0, at=(95.0, 20.0, 20.0)), outer)
    assert out["contained"] is False
    assert out["worst_distance_mm"] > vg.PROTRUSION_MIN_MM
    assert out["worst_distance_at"] is not None, "and it says where"


def test_the_distance_contradiction_exonerates_only_points_actually_on_the_body(outer):
    """And the direction it was built for, stated as the quantity rather than a count.

    A face flush with the boundary has its points *on* the body, so their exact
    distance to it is zero however the classifier votes about which side they
    are on. That is what makes the bar tight enough to be meaningless as a
    calibration: the ties measured on `SpiralTest.step` sit at 0.000000 mm
    (median 9.95e-10) against protrusions of 0.239 to 0.771 mm.
    """
    flush = box(10.0, 10.0, 10.0, at=(90.0, 20.0, 20.0))
    ev = vg.surface_points_outside(flush, outer)
    assert ev["contained"] is True
    assert ev["worst_distance_mm"] <= vg.PROTRUSION_MIN_MM
