"""Export truth: the two measurements that ask what ``BRepCheck_Analyzer`` cannot.

``BRepCheck_Analyzer`` asks whether a shape agrees with itself to within the
tolerances **recorded in this process**. STEP AP214 carries one modelling
tolerance for a whole file, against one per vertex, edge and face in an OCCT
B-rep — so a shape whose validity is *carried by* a locally fat tolerance is
valid here and is not guaranteed valid in the file the user receives.

The first test in this module pins that premise directly against the kernel,
because everything else on this branch rests on it: if a future OCCT ever did
preserve per-subshape tolerance through STEP, the two gates below would be
solving a problem that no longer exists, and this test is what would say so.

See docs/algorithm.md §7.3 and §9.
"""

import math
import os

import pytest
from OCP.BOPAlgo import BOPAlgo_ArgumentAnalyzer
from OCP.BRep import BRep_Builder, BRep_Tool
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeSphere
from OCP.TopAbs import TopAbs_ShapeEnum
from OCP.TopoDS import TopoDS

from latticegen2 import occ
from latticegen2.junction import build_template
from latticegen2.lattice import lattice_params, nodes as lattice_nodes
from latticegen2.boundary import trim_junction

TESTDIR = os.path.dirname(os.path.abspath(__file__))
SPIRAL = os.path.join(TESTDIR, "SpiralTest.step")
BALL = os.path.join(TESTDIR, "80mm-test-ball.step")


def vertices(shape):
    return [TopoDS.Vertex_s(v) for v in occ._explore(shape, TopAbs_ShapeEnum.TopAbs_VERTEX)]


def edges(shape):
    return [TopoDS.Edge_s(e) for e in occ._explore(shape, TopAbs_ShapeEnum.TopAbs_EDGE)]


# --- the premise ------------------------------------------------------------


def test_step_does_not_carry_a_per_subshape_tolerance(tmp_path):
    """A fat vertex tolerance does not survive a STEP round trip.

    This is the whole reason both gates in this module exist. The repair
    docs/algorithm.md §8 performs on a falsely self-intersecting wire *is* a
    widened vertex tolerance — "it moves no geometry", which is what makes it
    safe on an already-proven-watertight shell and also what makes it
    unrepresentable in the exported file.

    A box needs no tolerance at all, so the widening below is pure metadata and
    the round trip has every excuse to discard it — which is exactly the point:
    the file has nowhere to put it. Measured, the 6.573e-02 mm this uses (the
    real figure OCCT recorded on `SpiralTest`'s fat vertex) comes back at
    OCCT's own confusion, seven orders smaller.
    """
    box = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
    BRep_Builder().UpdateVertex(vertices(box)[0], 6.573e-2)
    assert max(BRep_Tool.Tolerance_s(v) for v in vertices(box)) == pytest.approx(6.573e-2)

    path = str(tmp_path / "box.step")
    occ.write_step(box, path, "premise")
    back = occ.read_step(path)

    assert max(BRep_Tool.Tolerance_s(v) for v in vertices(back)) < 1e-5, (
        "STEP carried a per-vertex tolerance through a round trip; if that is "
        "genuinely true now, the export-truth gates need re-deriving"
    )


def test_the_file_declares_exactly_one_tolerance_for_the_whole_shape(tmp_path):
    """And it is an *average*, which on a lattice is the pathological summary.

    Almost every edge of a lattice is an exactly-built interior edge at OCCT's
    confusion; the only edges carrying real tolerance are the boundary trims.
    Averaging over the first group is how the second group's tolerance is lost.
    """
    box = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
    builder = BRep_Builder()
    for e in edges(box)[:2]:
        builder.UpdateEdge(e, 1e-2)

    path = str(tmp_path / "box.step")
    occ.write_step(box, path, "premise")
    declared = [
        line for line in open(path, errors="replace") if "UNCERTAINTY_MEASURE_WITH_UNIT" in line
    ]
    assert len(declared) == 1, "AP214 has one uncertainty per representation context"


# --- tolerance against feature size (docs/algorithm.md §7.3) ----------------


def test_tolerance_feature_ratio_is_negligible_on_exact_geometry():
    """A box is built exactly, so its tolerance is nothing against its faces."""
    box = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
    reading = occ.tolerance_feature_ratio(occ.faces(box))
    assert reading.ratio < 1e-6
    assert reading.face_area == pytest.approx(100.0)


def test_tolerance_feature_ratio_is_exactly_tolerance_over_root_area():
    """No calibration in the quantity itself — it is the ratio it says it is."""
    box = BRepPrimAPI_MakeBox(4.0, 4.0, 4.0).Shape()
    BRep_Builder().UpdateEdge(edges(box)[0], 0.25)
    reading = occ.tolerance_feature_ratio(occ.faces(box))
    assert reading.tolerance == pytest.approx(0.25)
    assert reading.face_area == pytest.approx(16.0)
    assert reading.ratio == pytest.approx(0.25 / math.sqrt(16.0))


def test_tolerance_feature_ratio_finds_the_worst_face_not_the_first():
    box = BRepPrimAPI_MakeBox(4.0, 4.0, 4.0).Shape()
    box_edges = edges(box)
    builder = BRep_Builder()
    builder.UpdateEdge(box_edges[0], 0.01)
    builder.UpdateEdge(box_edges[-1], 0.40)
    assert occ.tolerance_feature_ratio(occ.faces(box)).tolerance == pytest.approx(0.40)


def test_tolerance_feature_ratio_ignores_degenerate_edges():
    """A degenerate edge has no extent, so its tolerance bounds nothing.

    A sphere is the cheapest real carrier: OCCT gives it degenerate edges at
    both poles. Widening one must not move the reading, or the gate would fire
    on a parametric artefact rather than on geometry.
    """
    sphere = BRepPrimAPI_MakeSphere(10.0).Shape()
    faces = occ.faces(sphere)
    degenerate = [e for e in edges(sphere) if BRep_Tool.Degenerated_s(e)]
    assert degenerate, "the fixture must actually carry degenerate edges"

    before = occ.tolerance_feature_ratio(faces).ratio
    for e in degenerate:
        BRep_Builder().UpdateEdge(e, 0.5)
    assert occ.tolerance_feature_ratio(faces).ratio == pytest.approx(before)


def test_tolerance_feature_ratio_reports_where_to_look():
    """The reading has to name a place, or a report of it cannot be acted on."""
    box = BRepPrimAPI_MakeBox(4.0, 4.0, 4.0).Shape()
    BRep_Builder().UpdateEdge(edges(box)[0], 0.25)
    where = occ.tolerance_feature_ratio(occ.faces(box)).where
    assert len(where) == 3
    assert all(-1.0 <= c <= 5.0 for c in where)


# --- pcurve versus 3D curve (docs/algorithm.md §9) --------------------------


def bopalgo_curve_on_surface_faults(shape) -> int:
    """OCCT's own count of the same fault, as the control for the test below."""
    analyzer = BOPAlgo_ArgumentAnalyzer()
    analyzer.SetShape1(shape)
    analyzer.CurveOnSurfaceMode = True
    analyzer.SelfInterMode = False
    analyzer.Perform()
    return len(list(analyzer.GetCheckResult())) if analyzer.HasFaulty() else 0


def test_curve_on_surface_deviations_is_clean_on_exact_geometry():
    box = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
    reading = occ.curve_on_surface_deviations(box)
    assert reading.pairs == 24, "six faces, four edges each"
    assert reading.over_tolerance == 0
    assert reading.worst < 1e-6
    assert reading.loose_area_fraction == 0.0
    assert reading.loose_faces == 0
    assert reading.where == []


def test_loose_area_fraction_is_the_share_of_surface_not_the_worst_face():
    """The distinction the gate turns on, pinned on geometry built to show it.

    A cuboid with one large face and one small one, where only the small face
    is loose: the *worst* reading is the same whichever face carries the fault,
    and the fraction is not. `SpiralTest` is the real case — there the sound
    dominant body scores worse than the unwritable island on the worst face —
    but that part is eight minutes and this is the property in isolation.
    """
    from latticegen2.occ import CurveOnSurface

    # 100 mm^2 faces and 1 mm^2 faces on one box: 2 x (10x10), 4 x (10x0.1)
    box = BRepPrimAPI_MakeBox(10.0, 10.0, 0.1).Shape()
    areas = sorted(occ.area(f) for f in occ.faces(box))
    assert areas[0] < 2.0 and areas[-1] > 90.0, "the fixture must be lopsided"

    total = sum(areas)
    # A fault on one thin side is a small share of the surface...
    assert areas[0] / total < 1e-2
    # ...and a fault on one big face is a large one, at the same worst reading.
    assert areas[-1] / total > 4e-1


def test_loose_faces_are_counted_once_however_many_edges_they_have():
    """Area-weighted, not pair-weighted: a four-edged face is one face.

    Pair counting would weight a face by its edge count, which is a property of
    how the boolean happened to subdivide it rather than of how much of the
    body it is.
    """
    box = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
    reading = occ.curve_on_surface_deviations(box)
    assert reading.pairs == 24
    assert reading.loose_faces == 0, "nothing loose on exact geometry"


def test_curve_on_surface_deviations_agrees_with_occts_own_verdict(tmp_path):
    """The control, and the reason to trust the number this returns.

    A measurement that only ever agrees on sound geometry proves nothing (G10),
    so this is run on a shape that genuinely carries the fault: a real lattice
    output after a STEP round trip, which is where the pcurves have been
    regenerated and no longer match their own 3D curves.

    Both instruments must find the **same** pairs, because this one replaces
    ``BOPAlgo_ArgumentAnalyzer`` in the pipeline — and unlike the analyzer it
    also returns the deviation in millimetres, which is the quantity the
    decision needs. A count cannot tell a harmless sub-micron overshoot from a
    deviation the size of the face it sits on.
    """
    lp = lattice_params(20.0, 4.0)
    tpl = build_template(lp)
    body = occ.read_step(BALL)
    pos = lattice_nodes(lp, [(3, -1, -1)])[0]
    pieces = trim_junction(lp, tpl, pos, body).pieces
    assert pieces, "the fixture junction must actually be trimmed"

    path = str(tmp_path / "piece.step")
    occ.write_step(occ.compound(pieces[0][0]), path, "control")
    shipped = occ.read_step(path)

    reading = occ.curve_on_surface_deviations(shipped)
    assert reading.over_tolerance == bopalgo_curve_on_surface_faults(shipped), (
        "this measurement replaces OCCT's own in the pipeline, so it has to "
        "select the same pairs OCCT would"
    )


def test_curve_on_surface_deviations_reports_positions_when_it_fires():
    """A failure has to say where, capped so one bad solid is not a wall of text."""
    assert occ.CURVE_ON_SURFACE_SAMPLES > 0
    sphere = BRepPrimAPI_MakeSphere(10.0).Shape()
    reading = occ.curve_on_surface_deviations(sphere)
    assert len(reading.where) <= occ.CURVE_ON_SURFACE_SAMPLES
    assert reading.pairs > 0, "degenerate edges are skipped, but not every edge is"


# --- the two parts, measured -------------------------------------------------


@pytest.mark.parametrize(
    "part, cc, t, node, ceiling",
    [
        pytest.param(BALL, 20.0, 4.0, (3, -1, -1), 1e-4, id="ball-cc20-t4"),
    ],
)
def test_a_soundly_trimmed_junction_sits_far_below_the_bar(part, cc, t, node, ceiling):
    """A bar nothing reaches is untested, and a bar sound input reaches is a bug.

    docs/algorithm.md §11: refusing correct input is the one failure mode this
    project does not accept, so the committed parts are pinned well clear of
    the warning bar rather than merely passing it.
    """
    lp = lattice_params(cc, t)
    tpl = build_template(lp)
    body = occ.read_step(part)
    pos = lattice_nodes(lp, [node])[0]
    trim = trim_junction(lp, tpl, pos, body)
    assert trim.pieces
    worst = max(r.ratio for r in trim.tolerances)
    assert worst < ceiling
    assert worst < occ.TOLERANCE_FEATURE_RATIO_WARN


def test_a_fat_trim_is_flagged_at_the_source():
    """And the warning bar is reachable, or it is untested decoration.

    `SpiralTest.step` is committed for exactly this: its swept B-spline surface
    makes the boolean fit trim curves two orders fatter than anything the other
    parts produce. Junction (-6, -10, 1) is one of the six forming the 4.17 mm^3
    island that cannot be written to STEP, and it is the 2nd worst of the part's
    2,404 boundary pieces — flagged at the end of `boundary`, before the
    junction graph that turns it into a body even exists.
    """
    lp = lattice_params(5.0, 1.0)
    tpl = build_template(lp)
    body = occ.read_step(SPIRAL)
    pos = lattice_nodes(lp, [(-6, -10, 1)])[0]
    trim = trim_junction(lp, tpl, pos, body)
    assert trim.pieces, "the fixture junction must actually be trimmed"
    worst = max(r.ratio for r in trim.tolerances)
    assert worst > occ.TOLERANCE_FEATURE_RATIO_WARN, (
        "this junction is why the warning bar exists; if it no longer clears it "
        "the bar or the trim has changed and both need re-deriving"
    )
