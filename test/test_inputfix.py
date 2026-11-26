"""A vertex the input file places at the wrong point (docs/algorithm.md §7.5).

A user's file names, for an edge whose curve starts at one point, a vertex
1.334 mm away on a straight edge running on past it. The body the reader builds
uses every edge twice and, on a planar face, even reads valid; but one face's
boundary is open in its own parameters, two faces claim the same strip of the
straight edge, and support there could not be trimmed. ``test/input-defect.step``
(tools/make_input_defect_part.py) carries the same contradiction on a block's
roof, 13 mm wide; these tests ask the repair of it, and that it leaves sound
inputs alone.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from OCP.BRep import BRep_Tool
from OCP.BRepTools import BRepTools

from latticegen2 import inputfix, occ

HERE = os.path.dirname(__file__)
DEFECT = os.path.join(HERE, "input-defect.step")
BALL = os.path.join(HERE, "80mm-test-ball.step")
CYLINDER = os.path.join(HERE, "test-cylinder.STEP")
GAPS = os.path.join(HERE, "support-gaps.step")
REHEARSAL = os.path.join(HERE, "TD_HX_rehearsal_test.step")
HOLE_BLOCKS = [os.path.join(HERE, n) for n in ("Hole-block-basic.STEP",
                                                "Hole-block-chamfered.STEP")]

#: The block of tools/make_input_defect_part.py.
BLOCK_VOLUME = 40.0 * 30.0 * 20.0


def _near_roof_face(shape):
    """The roof face that lost its piece of the front edge: x 0 to 13."""
    return next(f for f in occ.faces(shape)
                if abs(occ.centroid(f)[2] - 20.0) < 1e-6 and occ.centroid(f)[0] < 13.0)


def test_the_fixture_carries_the_defect():
    """So that the tests below cannot pass by no longer testing it: one edge
    names a vertex 13 mm from where its curve ends, that vertex's tolerance
    covers the gap, and the near roof face meshes to half its area."""
    body = occ.read_step(DEFECT)
    found = inputfix.misplaced_vertices(body)
    assert len(found) == 1
    _, x, q, gap = found[0]
    assert gap == pytest.approx(13.0, abs=1e-9)
    assert np.allclose(q, (13.0, 0.0, 20.0))
    assert np.allclose(inputfix._vertex_xyz(x), (0.0, 0.0, 20.0))
    assert BRep_Tool.Tolerance_s(x) >= 13.0
    roof = _near_roof_face(body)
    assert inputfix._mesh_area(roof) == pytest.approx(0.5 * occ.area(roof), rel=1e-6)


def test_the_repair_splits_the_edge_where_the_curve_ends():
    body = occ.read_step(DEFECT)
    rep = inputfix.repair_input(body)
    assert [d.repaired for d in rep.defects] == [True]
    assert rep.defects[0].where == (13.0, 0.0, 20.0)
    assert rep.defects[0].span == pytest.approx(13.0)
    fixed = rep.shape
    assert inputfix.misplaced_vertices(fixed) == []
    assert rep.max_vertex_tol <= inputfix.REPAIR_MAX_TOL
    assert len(occ.solids(fixed)) == 1 and len(occ.faces(fixed)) == len(occ.faces(body))
    assert occ.is_valid(fixed)
    assert inputfix._edge_use_fault(fixed) is None
    assert occ.volume(fixed) == pytest.approx(BLOCK_VOLUME, rel=1e-12)
    roof = _near_roof_face(fixed)
    assert occ.area(roof) == pytest.approx(13.0 * 30.0, rel=1e-12)
    assert inputfix._mesh_area(roof) == pytest.approx(13.0 * 30.0, rel=1e-9)


def test_the_repaired_body_meshes_closed():
    """The whole body, at the ruler the classification mesh is built with and
    at the STEP round-trip check's: no mesh edge used other than twice."""
    rep = inputfix.repair_input(occ.read_step(DEFECT))
    for deflection in (0.15, occ.DEFAULT_MESH_DEFLECTION):
        BRepTools.Clean_s(rep.shape)
        tally = occ._mesh_and_count(rep.shape, deflection)
        assert tally.triangles > 0
        assert [k for k, n in tally.counts.items() if n != 2] == [], deflection


def test_a_refused_repair_leaves_the_body_as_built(monkeypatch):
    """A repair its check refuses must not leave the body half-changed: the
    misplaced vertex keeps the tolerance that covers the gap, since the
    defective body is consistent only with it, and the defect is reported."""
    body = occ.read_step(DEFECT)
    monkeypatch.setattr(inputfix, "_check", lambda *a: "refused for the test")
    rep = inputfix.repair_input(body)
    assert rep.shape.IsSame(body)
    assert [(d.repaired, d.reason) for d in rep.defects] == [(False, "refused for the test")]
    _, x, _, _ = inputfix.misplaced_vertices(rep.shape)[0]
    assert BRep_Tool.Tolerance_s(x) >= 13.0
    assert rep.max_vertex_tol >= 13.0


def test_no_split_is_planned_where_no_edge_runs_through_the_curve_end():
    """The repair acts only where the description it rests on holds: a curve
    end no edge at the vertex runs through is not this defect."""
    body = occ.read_step(DEFECT)
    e, x, q, _ = inputfix.misplaced_vertices(body)[0]
    plan, why = inputfix._plan(inputfix._Body(body), e, x, q + np.array([0.0, 0.5, 0.0]))
    assert plan is None and "not one" in why


@pytest.mark.parametrize("path", [BALL, CYLINDER, GAPS, *HOLE_BLOCKS])
def test_sound_inputs_are_left_alone(path):
    """Every run goes through the repair (no bypass), so on every committed
    part it must find nothing and return the body itself."""
    body = occ.read_step(path)
    rep = inputfix.repair_input(body)
    assert rep.defects == []
    assert rep.shape.IsSame(body)


def test_the_rehearsal_part_is_left_alone():
    """The one committed part with a vertex gap past :data:`GAP_TOL`: 0.016 mm,
    where no edge runs on to be split. It is declined, not repaired."""
    body = occ.read_step(REHEARSAL)
    rep = inputfix.repair_input(body)
    assert rep.shape.IsSame(body)
    assert len(rep.defects) == 1 and not rep.defects[0].repaired
    assert rep.defects[0].span == pytest.approx(0.0162, abs=1e-4)


def _recorder(tmp_path):
    from latticegen2.runlog import RunLog

    shown = []

    class Recorder(RunLog):
        def line(self, msg, console=None):
            shown.append((msg, console))

    return Recorder(str(tmp_path / "run.log"), verbose=False), shown


def test_a_repaired_input_is_a_note_the_user_sees(tmp_path):
    from latticegen2.pipeline import _report_input_repair

    rl, shown = _recorder(tmp_path)
    stats = {}
    rep = inputfix.repair_input(occ.read_step(DEFECT))
    _report_input_repair(rl, rep, 1.5, stats)
    notes = [(m, c) for m, c in shown if m.startswith("note: the input file")]
    assert len(notes) == 1 and notes[0][1] is True, "about the user's file: always shown"
    assert "[13.000, 0.000, 20.000] mm" in notes[0][0]
    assert not any(m.startswith("WARNING") for m, _ in shown)
    assert stats["input_defects_repaired"] == 1 and stats["input_defects_left"] == 0


def test_an_unrepaired_input_is_one_warning_with_its_position(tmp_path, monkeypatch):
    """Wider than a tenth of ``t``, the contradiction left as built is a
    warning naming both points. The vertex it fattens is the same place, so it
    gets no second warning."""
    from latticegen2.pipeline import _report_input_repair

    monkeypatch.setattr(inputfix, "_check", lambda *a: "refused for the test")
    rl, shown = _recorder(tmp_path)
    stats = {}
    _report_input_repair(rl, inputfix.repair_input(occ.read_step(DEFECT)), 1.5, stats)
    warnings = [m for m, c in shown if m.startswith("WARNING") and c is True]
    assert len(warnings) == 1
    assert "refused for the test" in warnings[0]
    assert "curve ends at [13.000, 0.000, 20.000] mm" in warnings[0]
    assert "vertex at [0.000, 0.000, 20.000] mm" in warnings[0]
    assert stats["input_defects_left"] == 1


def test_a_narrow_contradiction_left_as_built_is_logged_only(tmp_path):
    """`TD_HX_rehearsal_test.step` has an edge ending 0.016 mm from its vertex,
    where no edge runs on to split: left as built, and not worth a warning on a
    part that has been inspected and accepted."""
    from latticegen2.pipeline import _report_input_repair

    defect = inputfix.InputDefect((1.0, 2.0, 3.0), 0.0162, False, "not ours",
                                  (1.0, 2.0, 3.016))
    repair = inputfix.InputRepair(None, [defect], 0.0162, (1.0, 2.0, 3.016))
    rl, shown = _recorder(tmp_path)
    _report_input_repair(rl, repair, 1.4, {})
    assert not any(m.startswith("WARNING") for m, _ in shown)
    assert all(c is None for _, c in shown)
    assert any("left as built (not ours)" in m for m, _ in shown)


def test_a_sound_input_logs_its_worst_vertex_only(tmp_path):
    from latticegen2.pipeline import _report_input_repair

    rl, shown = _recorder(tmp_path)
    _report_input_repair(rl, inputfix.repair_input(occ.read_step(BALL)), 4.0, {})
    assert len(shown) == 1 and shown[0][1] is None
    assert shown[0][0].startswith("input vertex tolerances")
