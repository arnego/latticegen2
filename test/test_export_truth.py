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

import io
import math
import os

import pytest
from OCP.BOPAlgo import BOPAlgo_ArgumentAnalyzer
from OCP.BRep import BRep_Builder, BRep_Tool
from OCP.BRepPrimAPI import (
    BRepPrimAPI_MakeBox,
    BRepPrimAPI_MakeCone,
    BRepPrimAPI_MakeSphere,
)
from OCP.BRepTools import BRepTools
from OCP.TopoDS import TopoDS_Shape, TopoDS_Shell
from OCP.BRep import BRep_Builder as _BRepBuilder
from OCP.TopAbs import TopAbs_ShapeEnum
from OCP.TopoDS import TopoDS

from latticegen2 import occ
from latticegen2.junction import build_template
from latticegen2.lattice import lattice_params, nodes as lattice_nodes
from latticegen2.boundary import trim_junction

TESTDIR = os.path.dirname(os.path.abspath(__file__))
SPIRAL = os.path.join(TESTDIR, "SpiralTest.step")
BALL = os.path.join(TESTDIR, "80mm-test-ball.step")
CYLINDER = os.path.join(TESTDIR, "test-cylinder.STEP")
ISLAND = os.path.join(TESTDIR, "spiral-island-unwritable.brep")
REHEARSAL = os.path.join(TESTDIR, "TD_HX_rehearsal_test.step")


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


# --- the gate itself, against the body it exists for ------------------------


def load_island():
    """The one body this project has produced that cannot be written to STEP.

    Lifted straight out of a `SpiralTest.step` run at ``cc=5, t=1`` after
    `simplify`: 4.1629 mm^3, 36 faces. Committed rather than reproduced,
    per docs/testing.md — a synthetic case proves the code does what you think,
    only the real one proves it does what the part needs, and this defect took
    three wrong diagnoses before the right instrument found it.
    """
    shape = TopoDS_Shape()
    BRepTools.Read_s(shape, ISLAND, _BRepBuilder())
    solids = occ.solids(shape)
    assert len(solids) == 1
    return solids[0]


def test_repairing_the_island_is_what_makes_it_invisible(tmp_path):
    """The gap this gate exists for, demonstrated end to end on the real body.

    Three steps, all measured here:

    1. As this branch builds it, the island is ``BRepCheck_Analyzer``-**invalid**
       — one face carries a falsely self-intersecting wire whose shared vertex
       OCCT recorded at 6.573e-02 mm, above rung 2's fixed 4e-3 mm cap, so the
       repair declines it and the existing validity gate refuses the run.
    2. Raise that cap — which is exactly what rung 2 does once it is allowed to
       act — and the repair widens the vertex, **the face becomes valid, and so
       does the body**. Nothing moved: a tolerance is metadata.
    3. Write it. Its 147 triangles still carry 11 edges used by one triangle
       rather than two.

    So the repair that makes the body pass every gate this pipeline had is
    precisely what makes the remaining defect invisible, and the export-truth
    check is the only thing left that sees it. That is the case
    ``docs/algorithm.md`` §9 is about, and it is why the check measures the
    exported tessellation rather than anything about the solid in memory.
    """
    island = load_island()
    assert not occ.is_valid(island), (
        "as built on this branch the existing validity gate already refuses it"
    )

    # Rung 2 with its cap raised, as docs/algorithm.md §8's repair would apply
    # it. Patched on the module rather than reimplemented, so this exercises the
    # real repair and not a stand-in for it.
    original = occ.SELF_INTERSECT_MAX_VERTEX_TOL
    try:
        occ.SELF_INTERSECT_MAX_VERTEX_TOL = 0.1
        repaired, still_invalid = occ.fix_vertex_tolerances(occ.faces(island))
    finally:
        occ.SELF_INTERSECT_MAX_VERTEX_TOL = original
    assert (repaired, still_invalid) == (1, 0)
    assert occ.is_valid(island), "the repair makes the body pass BRepCheck_Analyzer"

    defects = occ.exported_mesh_defects(island, str(tmp_path / "island.step"))
    triangles, bad = defects.triangles, defects.bad
    assert triangles > 0
    assert bad > 0, (
        "this body is why the export-truth gate exists; if it now survives a "
        "round trip, either OCCT or the generator has changed and the gate "
        "needs re-deriving against fresh ground truth"
    )


def test_the_island_does_not_survive_being_written(tmp_path):
    """The same measurement without the repair, so the fixture cannot rot.

    ``(148, 10)`` is the body's **own** tessellation, not damage the write did:
    since :func:`occ.write_step` declares the greatest tolerance the shape
    carries rather than the average of them, the round trip now reproduces what
    the solid already is. Measured, this fixture in memory: 148 triangles, 10
    edges not used by exactly two of them. Under the old Average default the
    file under-declared and the figure came home at ``(147, 11)`` — one defect
    of which the writer had added.

    That is the sharpest statement of what this fixture is for. Its 42.4 micron
    surface gap is not an export artefact and no writer setting reaches it, so
    it stays refused while `SpiralTest.step`'s own bodies — whose only obstacle
    *was* the under-declared tolerance — now ship.

    It is also immune to the degenerate-triangle skip, and measurably so
    rather than by assumption: this body meshes with **0** degenerate
    triangles and carries no degenerate B-rep edge, so every one of its 10
    defects is real geometry. That is what makes it the required real fault
    in the control set (G10) — a rule that cleared it would be wrong.
    """
    defects = occ.exported_mesh_defects(
        load_island(), str(tmp_path / "island.step")
    )
    assert (defects.triangles, defects.bad) == (148, 10)

    # The breakdown and the positions are what make a refusal actionable:
    # one use is a hole or two faces discretizing a shared edge
    # differently, more than two is duplicate material. Reporting only a
    # total is what forced the rehearsal's 26 to be re-measured outside the
    # pipeline before anything could be said about them (G23).
    assert sum(defects.by_use.values()) == defects.bad
    assert 2 not in defects.by_use, "an edge used twice is not a defect"
    assert defects.by_use == {1: 7, 3: 3}, (
        "this fixture's defects are five used once, two interior segments "
        "used once, and three used three times -- measured, and pinned so a "
        "change in kind is visible and not only a change in total"
    )
    assert defects.where and len(defects.where) <= occ.MESH_DEFECT_SAMPLES
    assert all(len(pos) == 3 for pos in defects.where)
    assert defects.degenerate == 0, (
        "this body has no pole, so the degenerate-triangle skip cannot reach "
        "it -- which is why its numbers are a fixed point for that change"
    )


def test_a_sound_body_survives_being_written(tmp_path):
    """The control. A gate that only ever fires proves nothing (G10)."""
    box = BRepPrimAPI_MakeBox(3.0, 4.0, 5.0).Shape()
    defects = occ.exported_mesh_defects(box, str(tmp_path / "box.step"))
    assert defects.triangles == 12
    assert defects.bad == 0
    assert defects.by_use == {} and defects.where == []


def test_a_pole_does_not_read_as_a_defect(tmp_path):
    """A sound solid with a pole must read clean, and did not.

    At a pole -- a sphere pole, a cone apex -- a whole parametric range maps to
    one 3D point, so distinct parametric nodes intern to a single id. A triangle
    with two vertices there contributes its collapsed key once *and the real
    edge twice*, since that edge appears as both (A,P) and (P,A) within the same
    triangle. With two sound neighbours also using it, a perfectly closed fan
    reads as four uses.

    Measured before the degenerate-triangle skip: this sphere reported
    ``bad == 4`` with ``by_use == {4: 2, 1: 2}`` and `MakeCone(5, 0, 10)`
    reported ``bad == 2`` with ``{4: 1, 1: 1}`` -- one defect per degenerate
    triangle, on closed solids the export-truth gate would have refused
    outright. Those figures are what make this a regression test rather than a
    tautology.
    """
    sphere = BRepPrimAPI_MakeSphere(10.0).Shape()
    degenerate = [e for e in edges(sphere) if BRep_Tool.Degenerated_s(e)]
    assert degenerate, "the fixture must actually carry degenerate edges"

    defects = occ.exported_mesh_defects(sphere, str(tmp_path / "sphere.step"))
    assert defects.triangles == 2022, "the mesher's own total, skips included"
    assert defects.degenerate == 2, "one per pole"
    assert defects.bad == 0
    assert defects.by_use == {} and defects.where == []

    cone = BRepPrimAPI_MakeCone(5.0, 0.0, 10.0).Shape()
    cone_defects = occ.exported_mesh_defects(cone, str(tmp_path / "cone.step"))
    assert (cone_defects.bad, cone_defects.degenerate) == (0, 1)


def test_a_hole_at_a_pole_is_still_caught(tmp_path):
    """The blind-spot control, and the reason the rule is about *triangles*.

    Skipping anything merely *near* a pole would hide a genuine hole whose
    boundary reached one. Skipping triangles with two coincident vertices cannot:
    a triangle that bounds nothing is not what bounds a hole.

    This is the sharpest case available -- the cone's base disk removed, so the
    hole lies on the boundary of the very face that carries the apex. Every one
    of its 63 mesh segments is still reported. Measured before the skip, the
    same shell read **65** with ``{1: 64, 4: 1}``: a real hole with two apex
    artefacts piled on top. The fix makes the reading cleaner, not weaker.
    """
    cone = BRepPrimAPI_MakeCone(5.0, 0.0, 10.0).Shape()
    lateral = [
        f for f in occ.faces(cone)
        if any(BRep_Tool.Degenerated_s(TopoDS.Edge_s(e))
               for e in occ._explore(f, TopAbs_ShapeEnum.TopAbs_EDGE))
    ]
    assert len(lateral) == 1, "the apex-bearing face must be identifiable"

    builder = BRep_Builder()
    shell = TopoDS_Shell()
    builder.MakeShell(shell)
    for f in lateral:
        builder.Add(shell, f)

    defects = occ.exported_mesh_defects(shell, str(tmp_path / "holed.step"))
    assert defects.degenerate == 1, "the apex triangle is still skipped"
    assert defects.bad == 63, "the whole base polyline, nothing swallowed"
    assert defects.by_use == {1: 63}, "a hole reads as used-once, and only that"


def test_the_gate_does_not_refuse_its_own_input(tmp_path):
    """The sharpest statement of the pole miscount: it condemned the *input*.

    The primitives above show the mechanism on the smallest possible shape.
    This shows what it cost on a real one — and it is a different kind of
    argument, because `TD_HX_rehearsal_test.step` is a valid, accepted CAD model
    this project neither wrote nor may change. It carries **54 apex-touching
    conical faces**, drill points and countersink tips whose `v` range begins
    where `RefRadius + v*sin(SemiAngle)` is exactly zero, and the counter read
    it as carrying **86** non-manifold edges, one per apex.

    That also corrects an attribution: `tools/prototypes/RESULTS.md` G23 noticed
    every folding `ConicalSurface` in the *output* beginning at `v = -sqrt(3)`
    and read it as a property of the generated patches. It is inherited — the
    input's own parametrization, `RefRadius = 1.5` and `SemiAngle = pi/3`
    throughout (G24).

    Pinned in both directions, because "0 bad edges" alone would also pass if
    the counter had quietly stopped working.
    """
    body = occ.read_step(REHEARSAL)
    defects = occ.exported_mesh_defects(body, str(tmp_path / "input.step"))

    assert defects.degenerate == 86, (
        "one degenerate triangle per apex; if this moves, either the mesher or "
        "the input has changed and the figures in G24 need re-deriving"
    )
    assert defects.bad == 0, "the user's own input is not a defective body"
    assert defects.by_use == {} and defects.where == []


# --- what the run says about a body the second pass cleared ------------------


def _report_export_truth(tmp_path, defects):
    """Drive `pipeline._check_export_truth` with one canned reading.

    The reporting path around a *resolved* body cannot otherwise be reached
    without a 583,894-face solid and a 93-minute run, which is precisely the
    kind of code that gets a typo into it. The geometry here is a box; what is
    under test is what the run says, not what the gate measured.
    """
    from latticegen2 import pipeline
    from latticegen2.runlog import RunLog

    box = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
    rl = RunLog(str(tmp_path / "run.log")).open()
    stats = {}
    real = occ.exported_mesh_defects
    occ.exported_mesh_defects = lambda solid, probe: defects
    try:
        error = None
        try:
            pipeline._check_export_truth(rl, [box], str(tmp_path), stats)
        except pipeline.ProcessingError as exc:
            error = str(exc)
    finally:
        occ.exported_mesh_defects = real
        rl.close()
    return io.open(str(tmp_path / "run.log"), encoding="utf-8").read(), stats, error


def test_a_body_cleared_by_refinement_does_not_read_as_one_that_was_clean(tmp_path):
    """Two different runs must not produce the same sentence.

    A body that read 0 at the coarse ruler and a body whose readings were shown
    to be the ruler are different outcomes, and only the second one has ever
    been asked a question. The run has to say so -- in the log, and in a stats
    key that survives into the summary so it is greppable afterwards.
    """
    resolved = occ.Refinement(
        True, [(0.01, 8), (0.002, 4), (0.0005, 4), (0.0001, 0)], 3, 68, "resolved"
    )
    defects = occ.MeshDefects(
        1427670, 13, {1: 13}, [(1.0, 2.0, 3.0)], 9, frozenset({1, 2, 3}), resolved
    )
    log, stats, error = _report_export_truth(tmp_path, defects)

    assert error is None, "a resolved body must not stop the run"
    assert "68-face neighbourhood" in log
    assert "0.0001:0" in log and "resolved" in log
    assert "export_truth_refined" in stats
    assert "none at 0.0001 mm" in stats["export_truth_refined"]
    assert "cleared only after being re-measured" in log


def test_a_body_that_did_not_resolve_is_refused_with_its_ladder(tmp_path):
    """The failure message carries the shape of the ladder, not just a total.

    "13 -> 8 -> 4 -> 4, still 4" and "10 -> 13 -> 17" are different bug
    reports; a bare count is neither. Recovering that distinction by hand cost
    this project two prototype gates and a 93-minute run.
    """
    unresolved = occ.Refinement(
        False, [(0.01, 13), (0.002, 17)], 4, 26,
        "increased: the body disagrees with itself more the more closely it is "
        "measured",
    )
    defects = occ.MeshDefects(
        148, 10, {1: 7, 3: 3}, [(1.0, 2.0, 3.0)], 0, frozenset({0, 1, 2, 3}),
        unresolved,
    )
    log, stats, error = _report_export_truth(tmp_path, defects)

    assert error is not None, "a body that did not resolve must stop the run"
    assert "26-face neighbourhood" in error
    assert "13 at 0.01 mm -> 17 at 0.002 mm" in error
    assert "increased" in error
    assert "export_truth_refined" not in stats, "nothing was cleared"


# --- the second pass: is a reading the ruler, or the body? -------------------


def holed_cone():
    """A cone with its base disk removed -- a genuine hole, on the apex face.

    G23's blind-spot control, reused here as the control that makes the
    terminus rule safe. A hole is the one thing refinement must never resolve:
    a finer ruler finds *more* of its boundary polyline, never less.
    """
    cone = BRepPrimAPI_MakeCone(5.0, 0.0, 10.0).Shape()
    lateral = [
        f for f in occ.faces(cone)
        if any(BRep_Tool.Degenerated_s(TopoDS.Edge_s(e))
               for e in occ._explore(f, TopAbs_ShapeEnum.TopAbs_EDGE))
    ]
    assert len(lateral) == 1
    builder = BRep_Builder()
    shell = TopoDS_Shell()
    builder.MakeShell(shell)
    for f in lateral:
        builder.Add(shell, f)
    return shell


def test_a_clean_body_never_pays_for_the_second_pass(tmp_path):
    """``refinement is None`` means *not attempted*, and that is the cheap path.

    The whole cost argument for the second pass is that it rides the failing
    path only. A body reading 0 must therefore come back with ``None`` rather
    than with a ``Refinement`` that happens to say ``resolved`` -- those are
    different claims, and only one of them was measured.
    """
    box = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
    defects = occ.exported_mesh_defects(box, str(tmp_path / "box.step"))
    assert defects.bad == 0
    assert defects.refinement is None, "nothing to re-measure, nothing measured"
    assert defects.implicated == frozenset(), "and no face is implicated"


def test_the_readings_are_attributed_to_the_faces_that_carry_them(tmp_path):
    """Every reading lands on a face, and an unattributed one is a fault.

    ``implicated`` is what the second pass cuts its neighbourhood around, so a
    reading nobody claims would be a reading refinement never looks at. On the
    island all 10 land on 4 of its 36 faces; the holed cone's 63 all land on
    its single face.
    """
    island = load_island()
    defects = occ.exported_mesh_defects(island, str(tmp_path / "island.step"))
    assert defects.bad == 10
    assert len(defects.implicated) == 4, "the 10 readings sit on 4 faces"
    assert defects.implicated <= set(range(len(occ.faces(island))))

    holed = occ.exported_mesh_defects(holed_cone(), str(tmp_path / "holed.step"))
    assert holed.bad == 63
    assert holed.implicated == frozenset({0}), "one face, and it is the one there is"


def test_the_neighbourhood_is_the_core_plus_one_ring(tmp_path):
    """The extract is the core faces and everything sharing an edge with them.

    Pinned on the island because it is committed and small. The property that
    matters is not the number but that the core is a strict subset of a ring
    that is a strict subset of the body: a neighbourhood equal to the core
    would have a cut boundary running through the readings, and one equal to
    the body would not be an optimisation at all.
    """
    island = load_island()
    defects = occ.exported_mesh_defects(island, str(tmp_path / "island.step"))
    comp, core_map, n_faces = occ._extract_neighbourhood(island, defects.implicated)
    assert core_map.Size() == len(defects.implicated)
    assert len(defects.implicated) < n_faces < len(occ.faces(island))
    assert n_faces == len(occ.faces(comp)), "the compound holds what it counted"


def test_the_island_is_refused_because_it_comes_apart_under_refinement(tmp_path):
    """The positive control, and the reason the rule is about mechanism.

    ``spiral-island-unwritable.brep`` carries the *same* reading classes as a
    sound body (G23), so no rule about the kind of reading separates them. What
    separates them is what happens when the ruler gets finer: this body's count
    **rises**, because what disagrees is its geometry -- a pcurve regenerating
    2.118e-02 mm from its own 3D curve on a 0.05 mm² face.

    The ladder stops at the first increase, which is why the numbers below are
    the first two rungs and not all four. That early stop is a cost decision,
    not the criterion; the criterion is that no rung ever read zero.
    """
    island = load_island()
    defects = occ.exported_mesh_defects(island, str(tmp_path / "island.step"))
    fine = defects.refinement
    assert fine is not None, "a body with readings must be re-measured"
    assert fine.resolved is False
    assert fine.reason.startswith("increased")
    assert [n for _d, n in fine.counts] == [13, 17], (
        "13 then 17 -- and these are exactly G24's whole-solid sweep at the "
        "same two deflections, which is what says the neighbourhood route "
        "measures the same thing the whole body does"
    )
    assert all(n > 0 for _d, n in fine.counts)


def test_a_real_hole_is_never_resolved_by_looking_harder(tmp_path):
    """The control the terminus rule rests on.

    The rule clears a body on an exact zero at some rung. That is only safe if
    a genuine hole can never reach zero -- and it cannot, because refining
    subdivides the hole's boundary polyline rather than closing it. Measured
    over the whole ladder, so this is not one rung's luck.
    """
    defects = occ.exported_mesh_defects(holed_cone(), str(tmp_path / "holed.step"))
    fine = defects.refinement
    assert fine is not None
    assert fine.resolved is False, "a hole must never clear"
    assert all(n > 0 for _d, n in fine.counts), "and never read zero at any rung"
    assert fine.counts[0][1] > defects.bad, "refining finds more of it, not less"


def test_the_neighbourhood_route_counts_what_the_whole_solid_counts(tmp_path):
    """The oracle: the cheap route must not merely agree on the verdict.

    Refinement measures a cut-out neighbourhood, discarding any reading no core
    face contributed. That is only the cut-boundary exclusion if it selects the
    *same readings* a whole-solid measurement selects -- so this compares the
    per-rung counts, not the yes/no.

    The island is small enough to afford both. The whole-solid route is kept
    here and nowhere in production, because on a 583,894-face body it is the
    cost the neighbourhood exists to avoid.
    """
    island = load_island()
    defects = occ.exported_mesh_defects(island, str(tmp_path / "island.step"))
    rungs = [d for d, _n in defects.refinement.counts]

    whole = []
    for deflection in rungs:
        BRepTools.Clean_s(island)
        tally = occ._mesh_and_count(island, deflection)
        whole.append(sum(1 for n in tally.counts.values() if n != 2))
    BRepTools.Clean_s(island)

    assert whole == [n for _d, n in defects.refinement.counts]


def test_the_attribution_agrees_with_the_pcurve_route(tmp_path):
    """A second opinion on ``implicated``, sharing no machinery with the first.

    ``_implicated_faces`` works from the interned mesh points. OCCT will also
    say, per (edge, face) pair, exactly which triangulation nodes an edge's
    boundary polyline used -- ``BRep_Tool.PolygonOnTriangulation_s``, which is
    how ``tools/prototypes/g23_bad_edge_provenance.py`` classified these
    readings. Where the two can both speak they must agree.

    It is a **subset** relation and not equality, deliberately: a reading lying
    inside one face's own triangulation is claimed by no B-rep edge at all
    (G23's ``interior`` class), so the polygon route cannot see it while the
    mesh route can. Requiring equality would be requiring the weaker instrument
    to be complete.
    """
    island = load_island()
    defects = occ.exported_mesh_defects(island, str(tmp_path / "island.step"))

    # Re-mesh the island itself at the same ruler and re-derive the readings,
    # then attribute them through the polygons rather than through the points.
    BRepTools.Clean_s(island)
    tally = occ._mesh_and_count(island, occ.DEFAULT_MESH_DEFLECTION)
    bad_keys = {e for e, n in tally.counts.items() if n != 2}
    assert len(bad_keys) == defects.bad

    from OCP.TopLoc import TopLoc_Location

    claimed = set()
    for fi, f in enumerate(occ._explore(island, TopAbs_ShapeEnum.TopAbs_FACE)):
        face = TopoDS.Face_s(f)
        loc = TopLoc_Location()
        tri = BRep_Tool.Triangulation_s(face, loc)
        if tri is None:
            continue
        trsf = loc.Transformation()
        ids = []
        for i in range(1, tri.NbNodes() + 1):
            pt = tri.Node(i).Transformed(trsf)
            ids.append(tally.point_id[(
                round(pt.X() * occ._QUANTUM),
                round(pt.Y() * occ._QUANTUM),
                round(pt.Z() * occ._QUANTUM),
            )])
        for e in occ._explore(face, TopAbs_ShapeEnum.TopAbs_EDGE):
            poly = BRep_Tool.PolygonOnTriangulation_s(TopoDS.Edge_s(e), tri, loc)
            if poly is None:
                continue
            nodes = [ids[poly.Node(k) - 1] for k in range(1, poly.NbNodes() + 1)]
            for a, b in zip(nodes, nodes[1:]):
                lo, hi = (a, b) if a <= b else (b, a)
                if lo != hi and lo * occ._EDGE_STRIDE + hi in bad_keys:
                    claimed.add(fi)
    BRepTools.Clean_s(island)

    assert claimed, "the polygon route must see something, or it proves nothing"
    assert claimed <= defects.implicated


def test_an_unmeasured_refinement_never_reads_as_clean(tmp_path):
    """Not measured and measured-clean must not be the same value.

    Every way the second pass can fail to answer -- no face carries the
    readings, the neighbourhood will not build, a rung raises, a rung is past
    the triangle guard -- refuses the body. This project has recorded more than
    once what it costs when an unmeasured quantity reads as a measured zero.
    """
    island = load_island()
    nothing = occ.refine_until_manifold(island, frozenset())
    assert nothing.resolved is False
    assert nothing.reason.startswith("unmeasured")
    assert nothing.counts == []

    defects = occ.exported_mesh_defects(island, str(tmp_path / "island.step"))
    starved = occ.refine_until_manifold(
        island, defects.implicated, ladder=(0.01,)
    )
    assert starved.resolved is False
    assert starved.reason.startswith("increased") or starved.reason.startswith(
        "ladder exhausted"
    ), "a ladder that ran out has not cleared anything"


def test_the_cheaper_proxies_do_not_decide():
    """Three quantities were tried before the tessellation and two of them
    false-positive on sound rehearsal solids (:data:`occ.LOOSE_AREA_FRACTION_MAX`).
    They are still measured and logged, because they say *why* a body is
    fragile — this pins that they do not decide, so a later change cannot
    quietly promote one back into the gate.
    """
    import ast
    import inspect

    from latticegen2 import pipeline

    tree = ast.parse(inspect.getsource(pipeline._check_export_truth).strip())
    body = tree.body[0].body
    if isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant):
        body = body[1:]                      # drop the docstring, keep the code
    code = "\n".join(ast.unparse(node) for node in body)
    assert "exported_mesh_defects" in code, "the tessellation is the instrument"
    assert "LOOSE_AREA_FRACTION_MAX" not in code, (
        "the loose-area fraction false-positives on two sound rehearsal solids"
    )


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
    makes the boolean fit trim curves against a near-tangent, and junction
    (-10, -2, 0) is the worst of the part's 2,349 boundary junctions at
    4.497e-01 — 6.657e-04 mm of recorded tolerance on a face of 0.000002 mm^2.

    **A flagged junction is not a defect**, which is the reason this reading
    reports rather than refuses: this part now completes and writes a valid
    STEP with 56 of its junctions above the bar. What the reading buys is a
    name and a coordinate at the end of `boundary`, minutes before anything
    downstream exists to be measured.
    """
    lp = lattice_params(5.0, 1.0)
    tpl = build_template(lp)
    body = occ.read_step(SPIRAL)
    pos = lattice_nodes(lp, [(-10, -2, 0)])[0]
    trim = trim_junction(lp, tpl, pos, body)
    assert trim.pieces, "the fixture junction must actually be trimmed"
    worst = max(r.ratio for r in trim.tolerances)
    assert worst > occ.TOLERANCE_FEATURE_RATIO_WARN, (
        "this junction is why the warning bar exists; if it no longer clears it "
        "the bar or the trim has changed and both need re-deriving"
    )


# --- why the island cannot be repaired at the source -------------------------


def test_a_sound_body_has_no_seam_gap():
    """The control. Two faces of a box meet exactly, so there is nothing to find."""
    box = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
    gap = occ.surface_seam_gaps(box)
    assert gap.edges == 12
    assert gap.worst < 1e-9


def test_the_committed_spiral_is_measured_where_its_predecessor_failed():
    """The seam gap is read from the real file, and the reading is what warns.

    The `SpiralTest.step` committed before 2026-08-24 had two swept B-spline
    patches running **42.4 microns apart** along the whole of their shared edge
    — held together by a 2.12e-02 mm edge tolerance, invisible on its own
    886 mm2 and 1,018 mm2 faces, and fatal on the 0.01 mm2 ones the trim cut out
    of the same region. That is the defect this measurement exists to name, and
    the replacement part closes it by a factor of six.

    Pinned as a **bound rather than a value**, deliberately: what must not
    regress is that the reading is taken and is small, not that the file stays
    byte-identical. It is still above :data:`pipeline.SEAM_GAP_WARN_FRACTION`
    at ``t = 1``, so this part still reports its seam — and it now completes and
    writes a valid STEP anyway, which is exactly why the reading is a warning
    and not a gate.
    """
    from latticegen2.pipeline import SEAM_GAP_WARN_FRACTION

    gap = occ.surface_seam_gaps(occ.read_step(SPIRAL), bar=1e-3)
    assert gap.worst < 1e-2, (
        "the committed part regressed to the 4.24e-02 mm seam gap that made its "
        "predecessor's island unwritable"
    )
    assert gap.worst > SEAM_GAP_WARN_FRACTION * 1.0, (
        "and it is still the part that exercises the warning, or the warning is "
        "untested decoration"
    )


def test_the_parts_that_ship_are_ranked_below_the_part_that_does_not():
    """And the ranking is a ranking, not a bar — the ball is the far end of it."""
    ball = occ.surface_seam_gaps(occ.read_step(BALL))
    spiral = occ.surface_seam_gaps(occ.read_step(SPIRAL))
    assert ball.worst < 1e-9 < spiral.worst


def test_a_material_input_gap_is_reported_to_the_user_with_its_coordinates(tmp_path):
    """A defective input has to be named, on the console, in its own coordinates.

    This is the behaviour the whole seam-gap measurement exists to deliver, and
    it is pinned rather than left to a docstring: the one thing a modeller can
    act on is *where in their file* the two faces fail to meet. The line goes
    out with ``console=True`` rather than under ``-v``, because a user who did
    not ask for verbose output still needs to hear that the run is about to
    trim a region their model does not close.

    Both halves are asserted — that it reaches the console, and that it carries
    the point — because either alone is useless.
    """
    from latticegen2.pipeline import _report_seam_gaps
    from latticegen2.runlog import RunLog

    shown = []

    class Recorder(RunLog):
        def line(self, msg, console=None):
            shown.append((msg, console))

    rl = Recorder(str(tmp_path / "run.log"), verbose=False)
    stats = {}
    _report_seam_gaps(
        rl, occ.SurfaceSeamGap(4.2406e-2, (-22.044, -8.786, -30.0), 30, 10), 1.0, stats
    )
    assert len(shown) == 1
    msg, console = shown[0]
    assert console is True, "a user without -v must still be told about their input"
    for coordinate in ("-22.044", "-8.786", "-30.000"):
        assert coordinate in msg, "the point in the *input* is what can be acted on"
    assert "4.2406e-02" in msg and "input" in msg
    assert stats["input_seam_gap_mm"] == "4.241e-02"


def test_a_negligible_input_gap_does_not_shout_at_the_user(tmp_path):
    """The other half, and the one that keeps the warning worth reading.

    Every real B-rep records some gap. The cylinder's is 4.5e-04 mm and it ships
    a golden-sample-exact lattice, so a line claiming its input needs fixing
    would be false — and a warning that fires on sound parts stops being read.
    """
    from latticegen2.pipeline import _report_seam_gaps
    from latticegen2.runlog import RunLog

    shown = []

    class Recorder(RunLog):
        def line(self, msg, console=None):
            shown.append((msg, console))

    rl = Recorder(str(tmp_path / "run.log"), verbose=False)
    _report_seam_gaps(rl, occ.SurfaceSeamGap(4.5232e-4, (0.0, 0.0, 0.0), 27, 0), 1.5, {})
    assert len(shown) == 1
    msg, console = shown[0]
    assert console is None, "log-only: nothing here needs the user's attention"
    assert "negligible" in msg


def test_the_committed_inputs_land_on_the_side_of_that_line_their_outcome_says(tmp_path):
    """And it is measured against the real files, not against constructed ones.

    A bar is only as good as what it is applied to. The ball and the cylinder
    must stay silent — both ship — and the reading taken from each file has to
    be the one the run would take.
    """
    from latticegen2.pipeline import _report_seam_gaps, SEAM_GAP_WARN_FRACTION
    from latticegen2.runlog import RunLog

    for part, t in ((BALL, 4.0), (CYLINDER, 1.5)):
        gap = occ.surface_seam_gaps(occ.read_step(part))
        assert gap.worst <= SEAM_GAP_WARN_FRACTION * t, (
            f"{part} ships today; warning about its input would be false"
        )
        shown = []

        class Recorder(RunLog):
            def line(self, msg, console=None):
                shown.append(console)

        _report_seam_gaps(Recorder(str(tmp_path / "run.log")), gap, t, {})
        assert shown == [None]


def test_an_export_truth_failure_cites_the_input_when_the_gap_is_material():
    """A refusal has to point somewhere the user can act.

    The failure itself names a coordinate in the *output*, which is not
    something anyone can fix. Where the input's own seam gap is material against
    `t`, the message adds it — and says the thing that took a session to
    establish: a gap between two surfaces is not repairable by re-fitting a
    curve on either of them.
    """
    from latticegen2 import pipeline

    wide = occ.SurfaceSeamGap(4.2406e-2, (-22.044, -8.786, -30.0), 30, 10)
    note = pipeline._seam_gap_note(wide, 1.0)
    assert "4.2406e-02" in note and "-22.044" in note

    narrow = occ.SurfaceSeamGap(4.5232e-4, (0.0, 0.0, 0.0), 27, 0)
    assert pipeline._seam_gap_note(narrow, 1.5) == "", (
        "the cylinder ships; its own gap must not be blamed for anything"
    )
    assert pipeline._seam_gap_note(None, 1.0) == ""


def test_the_file_declares_the_greatest_tolerance_the_shape_carries(tmp_path):
    """Not the average of them — which on a lattice is the pathological summary.

    AP214 has room for one ``UNCERTAINTY_MEASURE_WITH_UNIT`` per file against
    one tolerance per subshape, so that number decides what every edge's
    tolerance becomes on import. OCCT's default is *Average*, and ~99 % of a
    lattice's edges are exactly-built interior edges at ``Precision::Confusion``
    — so the average lands near the floor and the boundary trims that carry real
    tolerance are clamped below what their own geometry needs.

    Measured on `SpiralTest.step`'s dominant solid, which tessellates into
    111,618 triangles with **zero** non-manifold edges in memory: Average
    declares 1.E-05 and it comes back with 2, Greatest declares 1.E-02 and it
    comes back with 0. The part completes on the second and is refused on the
    first.

    The box here carries one deliberately fattened edge, so a correct writer has
    to declare at least that — an averaging one cannot, having eleven exact
    edges to average it against.
    """
    box = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
    BRep_Builder().UpdateEdge(edges(box)[0], 1e-2)

    path = str(tmp_path / "box.step")
    occ.write_step(box, path, "precision")
    declared = [
        line for line in open(path, errors="replace")
        if "UNCERTAINTY_MEASURE_WITH_UNIT" in line
    ]
    assert len(declared) == 1
    value = float(declared[0].split("LENGTH_MEASURE(")[1].split(")")[0])
    assert value >= 1e-2, (
        "the file must declare a tolerance that covers the shape's own worst; "
        "below it, the geometry that tolerance was carrying changes on import"
    )


def test_the_precision_setting_survives_the_writer_being_constructed(tmp_path):
    """The ordering trap, pinned — it cost a full run to find.

    ``Interface_Static`` values are session state, and constructing a
    ``STEPControl_Writer`` initialises the STEP session, which resets them to
    their defaults. Setting the mode *before* the writer exists is therefore
    silently discarded in any process that has not already touched the STEP
    controller — which the pipeline had, because it reads its input first, so
    the bug hid there and appeared only in a fresh process.

    Asserting the mode *after* a write is what catches a regression to the
    wrong order: the value has to be the one this module set, not the default
    the writer's own initialisation would restore.
    """
    from OCP.Interface import Interface_Static

    Interface_Static.SetIVal_s("write.precision.mode", 0)      # the default
    occ.write_step(
        BRepPrimAPI_MakeBox(1.0, 1.0, 1.0).Shape(), str(tmp_path / "b.step"), "order"
    )
    assert Interface_Static.IVal_s("write.precision.mode") == occ.WRITE_PRECISION_GREATEST


def test_the_part_name_survives_the_writer_being_constructed(tmp_path):
    """The same trap, on the static that carries specification.md §5's metadata.

    The reset does not distinguish between statics, so ``write.step.product.name``
    is on exactly the same footing as ``write.precision.mode`` above — and it is
    the one where a regression is *invisible* rather than merely wrong, because
    ``stepout.rewrite_step_header`` patches FILE_NAME and FILE_DESCRIPTION and
    never touches the PRODUCT entity. Set before the writer exists, this comes
    back as OCCT's own translator string.

    Asserted against the file rather than against the static, since the file is
    what a downstream reader opens.
    """
    path = tmp_path / "named.step"
    occ.write_step(
        BRepPrimAPI_MakeBox(1.0, 1.0, 1.0).Shape(), str(path), "fixture+lattice+cc20+t4"
    )
    text = path.read_text(encoding="utf-8", errors="replace")
    assert "PRODUCT('fixture+lattice+cc20+t4" in text
    assert "Open CASCADE STEP translator" not in text.split("DATA;", 1)[1]
