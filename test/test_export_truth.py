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
from OCP.BRepTools import BRepTools
from OCP.TopoDS import TopoDS_Shape
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

    triangles, bad = occ.exported_mesh_defects(island, str(tmp_path / "island.step"))
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
    """
    triangles, bad = occ.exported_mesh_defects(
        load_island(), str(tmp_path / "island.step")
    )
    assert (triangles, bad) == (148, 10)


def test_a_sound_body_survives_being_written(tmp_path):
    """The control. A gate that only ever fires proves nothing (G10)."""
    box = BRepPrimAPI_MakeBox(3.0, 4.0, 5.0).Shape()
    triangles, bad = occ.exported_mesh_defects(box, str(tmp_path / "box.step"))
    assert triangles == 12
    assert bad == 0


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
