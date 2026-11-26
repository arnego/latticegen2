"""A face of the input that crosses its own boundary at a cusp (docs/algorithm.md §7.5).

Where two holes of a part touch tangentially, a planar face between them is
bounded by two circular edges leaving their shared vertex in the same direction.
A user's ``TDHX_error_report_test.stp`` writes those circles to six decimals,
so on its faces 166, 168 and 266 they overlap by up to 4.2e-9 mm and cross
again 8.3e-5 mm from the vertex, inside both edges: ``BRepCheck`` rightly calls
each such face invalid. A boolean against the body widens the vertex's
tolerance to cover the crossing, so the lattice body of the user's run at
``-cc 7 -t 1.4 --orient 180 0 0`` was valid in memory and failed ``BRepCheck``
once read back from its STEP file, on one 0.083 mm² face at that vertex.

The fixture is the real geometry, per docs/testing.md: the 19 x 11 x 15 mm block
``[1879, 1898] x [84, 95] x [978, 993]`` of the user's part, as
``BRepAlgoAPI_Common`` of that box with the part as read, holding both cusps of
face 266. The cut keeps the defect: the block is valid in memory, with the
vertex widened to 8.5e-5 mm, and invalid after a round trip.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from OCP.BRep import BRep_Builder, BRep_Tool
from OCP.BRepBuilderAPI import (
    BRepBuilderAPI_MakeEdge,
    BRepBuilderAPI_MakeFace,
    BRepBuilderAPI_MakeVertex,
    BRepBuilderAPI_MakeWire,
)
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakePrism
from OCP.BRepTools import BRepTools
from OCP.Geom import Geom_Circle
from OCP.TopoDS import TopoDS, TopoDS_Shape, TopoDS_Wire
from OCP.gp import gp_Ax2, gp_Circ, gp_Dir, gp_Pnt, gp_Vec

from latticegen2 import inputfix, occ
from latticegen2.boundary import trim_junction
from latticegen2.junction import build_template
from latticegen2.lattice import lattice_params

HERE = os.path.dirname(__file__)
FIXTURE = os.path.join(HERE, "cusp-block.brep")

#: The two cusps of face 266 in the block, at the R4.5 bore's two ends.
CUSPS = [(1893.261, 89.682, 983.6), (1884.411, 89.682, 983.6)]

#: The junction whose trim carried the invalid face in the user's run.
NODE = (196, -283, -257)


def _block():
    shape = TopoDS_Shape()
    BRepTools.Read_s(shape, FIXTURE, BRep_Builder())
    (block,) = occ.solids(shape)
    return block


def _round_trip(shape, tmp_path, name) -> list[bool]:
    path = str(tmp_path / f"{name}.step")
    occ.write_step(shape, path, name)
    return [occ.is_valid(s) for s in occ.solids(occ.read_step(path))]


def _trimmed(block):
    lp = lattice_params(7.0, 1.4, (180.0, 0.0, 0.0))
    tpl = build_template(lp)
    trim = trim_junction(lp, tpl, lp.B @ np.array(NODE, dtype=float), block, None)
    (faces, _tags, volume), = trim.pieces
    (shell,) = occ.shells(occ.sew(faces, 1e-6))
    return occ.make_solid(shell), volume


@pytest.fixture(scope="module")
def repaired():
    block = _block()
    return block, inputfix.repair_input(block)


def test_the_fixture_carries_the_defect(tmp_path):
    """Both cusps cross again 8.3e-5 mm from their vertex over a sliver 4.2e-9
    mm wide; the block is valid here and not once written. Without this the
    tests below prove nothing."""
    block = _block()
    found = inputfix.crossed_cusps(block)
    assert sorted(inputfix._round(c.q) for c in found) == sorted(CUSPS)
    for c in found:
        assert c.span == pytest.approx(8.326e-5, rel=1e-3)
        assert c.width == pytest.approx(4.24e-9, rel=1e-2)
    assert occ.is_valid(block)
    assert _round_trip(block, tmp_path, "block") == [False]


def test_the_trimmed_junction_is_valid_here_and_not_once_written(tmp_path):
    """The user's failure, on one junction: what the in-memory validity gate
    passes, the file's reader refuses."""
    solid, _ = _trimmed(_block())
    assert occ.is_valid(solid)
    assert _round_trip(solid, tmp_path, "junction") == [False]


def test_the_repair_moves_each_vertex_to_where_its_arcs_cross(repaired, tmp_path):
    block, rep = repaired
    assert [(c.vertex, c.repaired) for c in rep.crossings] == [(p, True) for p in CUSPS] \
        or sorted((c.vertex, c.repaired) for c in rep.crossings) == [(p, True) for p in sorted(CUSPS)]
    assert inputfix.crossed_cusps(rep.shape) == []
    assert occ.is_valid(rep.shape)
    assert _round_trip(rep.shape, tmp_path, "repaired") == [True]
    # The material moved is the sliver's own.
    assert occ.volume(rep.shape) == pytest.approx(occ.volume(block), abs=1e-6)
    # The vertex the boolean widened to cover the crossing is gone: the one at
    # the crossing needs no more than the kernel's confusion.
    for p in CUSPS:
        assert BRep_Tool.Tolerance_s(_vertex_near(block, p)) > 8e-5
        assert BRep_Tool.Tolerance_s(_vertex_near(rep.shape, p)) < 1e-6


def test_the_junction_trimmed_against_the_repaired_body_survives_being_written(repaired, tmp_path):
    block, rep = repaired
    solid, volume = _trimmed(rep.shape)
    assert volume == pytest.approx(_trimmed(block)[1], abs=1e-6)
    assert occ.is_valid(solid)
    assert _round_trip(solid, tmp_path, "junction") == [True]


def test_a_refused_repair_leaves_the_body_as_built(monkeypatch):
    """Every cusp of the face is reported, and nothing it touched keeps a
    changed tolerance."""
    block = _block()
    before = BRep_Tool.Tolerance_s(_vertex_near(block, CUSPS[1]))
    monkeypatch.setattr(inputfix, "_check_cusps", lambda *a: "refused for the test")
    rep = inputfix.repair_input(block)
    assert rep.shape.IsSame(block)
    assert sorted(c.vertex for c in rep.crossings) == sorted(CUSPS)
    assert all(not c.repaired and c.reason == "refused for the test" for c in rep.crossings)
    assert BRep_Tool.Tolerance_s(_vertex_near(rep.shape, CUSPS[1])) == before


def _vertex_near(shape, p):
    from OCP.TopAbs import TopAbs_ShapeEnum
    from OCP.TopExp import TopExp_Explorer

    ex = TopExp_Explorer(shape, TopAbs_ShapeEnum.TopAbs_VERTEX)
    while ex.More():
        v = TopoDS.Vertex_s(ex.Current())
        if np.linalg.norm(inputfix._vertex_xyz(v) - np.array(p)) < 1e-3:
            return v
        ex.Next()
    raise AssertionError(p)


# -- Synthetic: the user's numbers, the overlap swept --------------------------

ORIGIN = np.array([1883.427673, 89.5, 983.6])


def _cusp_prism(overlap: float, at_crossing: bool = False):
    """Two holes of radius 1 and 4.5 whose walls touch, ``overlap`` mm into
    each other, and the floor between them: the shape of face 266's cusps. The
    vertex is at the nominal tangent point, or with ``at_crossing`` where the
    two circles really cross."""
    r1, r2 = 1.0, 4.5
    d = r1 + r2 - overlap

    def pnt(x, y):
        return gp_Pnt(*(ORIGIN + np.array([x, y, 0.0])))

    def vtx(x, y, tol=1e-7):
        v = BRepBuilderAPI_MakeVertex(pnt(x, y)).Vertex()
        BRep_Builder().UpdateVertex(v, tol)
        return v

    if at_crossing:
        x = (d * d + r1 * r1 - r2 * r2) / (2 * d)
        y = -np.sqrt(r1 * r1 - x * x)
    else:
        x, y = (r1 + d - r2) / 2, 0.0
    vc = vtx(x, y, max(1e-7, abs(overlap)))
    small = Geom_Circle(gp_Circ(gp_Ax2(pnt(0, 0), gp_Dir(0, 0, 1)), r1))
    big = Geom_Circle(gp_Circ(gp_Ax2(pnt(d, 0), gp_Dir(0, 0, 1)), r2))
    vb1, vb2 = vtx(0.0, -r1), vtx(d, -r2)
    a_big = np.arctan2(y, x - d) % (2 * np.pi)
    edges = [
        BRepBuilderAPI_MakeEdge(small, vb1, vc, -np.pi / 2, np.arctan2(y, x)).Edge(),
        BRepBuilderAPI_MakeEdge(big, vc, vb2, a_big, 1.5 * np.pi).Edge(),
    ]
    vl, vk = vtx(d, -6.0), vtx(0.0, -6.0)
    edges += [BRepBuilderAPI_MakeEdge(vb2, vl).Edge(), BRepBuilderAPI_MakeEdge(vl, vk).Edge(),
              BRepBuilderAPI_MakeEdge(vk, vb1).Edge()]
    wire = BRepBuilderAPI_MakeWire()
    for e in edges:
        wire.Add(e)
    face = BRepBuilderAPI_MakeFace(wire.Wire(), True).Face()
    return BRepPrimAPI_MakePrism(face, gp_Vec(0, 0, 7.5)).Shape()


def test_the_synthetic_cusp_reproduces_the_user_part(tmp_path):
    """The user's overlap is invalid at the nominal vertex and valid at the
    crossing, in memory and once written: the repair's whole premise."""
    built = _cusp_prism(4.2e-9)
    assert not occ.is_valid(built)
    moved = _cusp_prism(4.2e-9, at_crossing=True)
    assert occ.is_valid(moved)
    assert _round_trip(moved, tmp_path, "moved") == [True]


@pytest.mark.parametrize("overlap", [4.2e-9, 1e-7])
def test_a_synthetic_cusp_is_repaired(overlap, tmp_path):
    body = _cusp_prism(overlap)
    rep = inputfix.repair_input(body)
    # The cusp is on the prism's floor and roof: the line is moved twice.
    assert [c.repaired for c in rep.crossings] == [True, True], rep.crossings
    assert occ.is_valid(rep.shape)
    assert _round_trip(rep.shape, tmp_path, "rep") == [True]
    assert occ.volume(rep.shape) == pytest.approx(occ.volume(body), abs=1e-6)


@pytest.mark.parametrize("overlap", [0.0, -1e-9, -5e-8])
def test_circles_that_touch_or_part_are_left_alone(overlap):
    """Exactly tangent or apart, the arcs never cross: nothing to repair, and
    the body comes back as the same object."""
    body = _cusp_prism(overlap)
    assert occ.is_valid(body)
    rep = inputfix.repair_input(body)
    assert rep.crossings == () and rep.shape.IsSame(body)


def test_a_sliver_wider_than_a_rounding_is_reported_not_moved():
    """Two holes that really overlap are the model, not its rounding: the
    crossing is reported and the body left as built."""
    body = _cusp_prism(2e-6)
    rep = inputfix.repair_input(body)
    assert rep.shape.IsSame(body)
    assert [c.repaired for c in rep.crossings] == [False, False]
    assert "sliver" in rep.crossings[0].reason


def test_a_face_whose_wire_the_file_stores_reversed_is_rebuilt_as_it_was():
    """The user's faces 166 and 168 are reversed and store their wire reversed.
    A rebuild once turned such a wire round twice, so every repair there came
    out invalid. Built here the same way from a box face."""
    box = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
    face = occ.faces(box)[0]
    builder = BRep_Builder()
    (wire, edges), = inputfix._wire_edges(face)
    stored = TopoDS_Wire()
    builder.MakeWire(stored)
    for e in reversed(edges):
        builder.Add(stored, e.Reversed())
    stored.Closed(True)
    flipped = TopoDS.Face_s(TopoDS.Face_s(face.Oriented(inputfix.FORWARD)).EmptyCopied())
    builder.Add(flipped, stored.Reversed())
    flipped = TopoDS.Face_s(flipped.Oriented(face.Orientation()))
    assert BRepCheck_Analyzer(flipped).IsValid()
    (w, _), = inputfix._wire_edges(flipped)
    assert w.Orientation() == inputfix.REVERSED
    rebuilt = inputfix._rebuild_face(flipped, [])
    assert BRepCheck_Analyzer(rebuilt).IsValid()
    assert occ.area(rebuilt) == pytest.approx(occ.area(face), rel=1e-12)


def _recorder(tmp_path):
    from latticegen2.runlog import RunLog

    shown = []

    class Recorder(RunLog):
        def line(self, msg, console=None):
            shown.append((msg, console))

    return Recorder(str(tmp_path / "run.log"), verbose=False), shown


def test_a_repaired_crossing_is_a_note_the_user_sees(repaired, tmp_path):
    from latticegen2.pipeline import _report_input_repair

    _block_, rep = repaired
    rl, shown = _recorder(tmp_path)
    stats = {}
    _report_input_repair(rl, rep, 1.4, stats)
    notes = [(m, c) for m, c in shown if m.startswith("note: a face of the input")]
    assert len(notes) == 1 and notes[0][1] is True
    assert "[1884.411, 89.682, 983.600] mm" in notes[0][0]
    assert not any(m.startswith("WARNING") for m, _ in shown)
    assert stats["input_crossings_repaired"] == 2 and stats["input_crossings_left"] == 0


def test_an_unrepaired_crossing_is_one_warning_with_its_positions(tmp_path, monkeypatch):
    from latticegen2.pipeline import _report_input_repair

    monkeypatch.setattr(inputfix, "_check_cusps", lambda *a: "refused for the test")
    rl, shown = _recorder(tmp_path)
    stats = {}
    _report_input_repair(rl, inputfix.repair_input(_block()), 1.4, stats)
    warnings = [m for m, c in shown if m.startswith("WARNING") and c is True]
    assert len(warnings) == 1
    assert "crosses itself" in warnings[0] and "refused for the test" in warnings[0]
    assert "vertex at [1893.261, 89.682, 983.600] mm" in warnings[0]
    assert stats["input_crossings_left"] == 2
