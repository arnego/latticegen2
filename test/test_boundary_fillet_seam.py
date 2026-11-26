"""Junctions the exact trim keeps whole across a seam of the input (§7.2).

Three fillets around a bore meet one another 3.5e-05 mm apart, a gap only their
edges' tolerance covers, at the corner where they meet a plane. Where a strut
crosses that corner, ``BRepAlgoAPI_Common`` at fuzzy 0 leaves the strut face
unsplit and keeps it whole, so the trimmed junction reaches up to 2.4 mm outside
the body. §7.2's per-half-strut and local-block repairs ask the same kernel at
the same tolerance and fail the same way, and the junction used to be dropped:
a gap in the lattice, reported as "lattice junctions left out".

The fixture is the real geometry, per docs/testing.md: a 18 x 18 x 12 mm block
of `TDHX_error_report_test.stp` (a user's part, received 2026-11-23) cut around
the three junctions that run dropped, at ``cc=7, t=1.4, --orient 180 0 0``.
Cut from the body, it still reproduces all three exactly.
"""

from __future__ import annotations

import os

import numpy as np
import pytest
from OCP.BRep import BRep_Builder
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.BRepTools import BRepTools
from OCP.TopAbs import TopAbs_State
from OCP.TopTools import TopTools_ListOfShape
from OCP.TopoDS import TopoDS_Shape
from OCP.gp import gp_Pnt

from latticegen2 import boundary, occ
from latticegen2.boundary import _piece_vertices, trim_boundary, trim_junction
from latticegen2.classify import OutsideProbe, stage_mesh, tessellate_surface
from latticegen2.junction import build_template
from latticegen2.lattice import lattice_params
from latticegen2.parallel import WorkerPool

HERE = os.path.dirname(__file__)
FIXTURE = os.path.join(HERE, "fillet-seam-block.brep")

# The three junctions, and what each trims to at fuzzy 1e-4. Point
# classification of the block agrees with each of these trims at every one of
# 20,000 points (measured on the whole part), so they are the truth here.
JUNCTIONS = {
    (196, -266, -274): 9.6398,
    (197, -268, -273): 15.0748,
    (197, -267, -274): 10.1194,
}


def _read_block():
    shape = TopoDS_Shape()
    BRepTools.Read_s(shape, FIXTURE, BRep_Builder())
    (block,) = occ.solids(shape)
    return block


@pytest.fixture(scope="module")
def setup():
    # The probe's mesh is made on a copy: the trims are asked of a block with no
    # triangulation on it, as the workers' staged `input.brep` has none.
    block = _read_block()
    lp = lattice_params(7.0, 1.4, (180.0, 0.0, 0.0))
    tpl = build_template(lp)
    mesh = tessellate_surface(_read_block(), lp)
    probe = OutsideProbe(mesh, 2.0 * mesh.deviation)
    return lp, tpl, block, mesh, probe


def _at(lp, node):
    return lp.B @ np.asarray(node, dtype=float)


def _exact_common(tpl, pos, block):
    algo = BRepAlgoAPI_Common()
    args = TopTools_ListOfShape()
    args.Append(tpl.solid.Moved(occ.translation(pos)))
    tools = TopTools_ListOfShape()
    tools.Append(block)
    algo.SetArguments(args)
    algo.SetTools(tools)
    algo.Build()
    return algo.Shape()


@pytest.mark.parametrize("node", list(JUNCTIONS))
def test_the_fixture_carries_the_defect(setup, node):
    """The exact trim keeps (nearly) the whole junction and leaves material
    outside, which the probe sees: without this the tests below prove nothing."""
    lp, tpl, block, _mesh, probe = setup
    got = _exact_common(tpl, _at(lp, node), block)
    vol = sum(occ.volume(s) for s in occ.solids(got))
    assert vol > JUNCTIONS[node] + 5.0
    pieces = [(occ.faces(s), None, None) for s in occ.solids(got)]
    assert probe.worst_outside(_piece_vertices(pieces)) > 1.0


@pytest.mark.parametrize("node", list(JUNCTIONS))
def test_a_junction_across_the_seam_is_trimmed_not_dropped(setup, node):
    lp, tpl, block, _mesh, probe = setup
    trim = trim_junction(lp, tpl, _at(lp, node), block, probe)
    assert not trim.dropped
    assert trim.outside_mm == 0.0
    assert trim.fuzzy in boundary.JUNCTION_FUZZY_LADDER
    assert sum(v for _f, _t, v in trim.pieces) == pytest.approx(JUNCTIONS[node], abs=1e-3)


def test_the_kept_trim_is_the_junction_inside_the_block(setup):
    """Asked of points, not of a boolean: no point of the trim lies outside the
    block, and no point of the junction inside the block is missing from it."""
    lp, tpl, block, _mesh, probe = setup
    node = (197, -267, -274)
    pos = _at(lp, node)
    trim = trim_junction(lp, tpl, pos, block, probe)
    (faces, _tags, _vol), = trim.pieces
    (shell,) = occ.shells(occ.sew(faces, 1e-6))
    solid = occ.make_solid(shell)
    assert occ.volume(solid) == pytest.approx(JUNCTIONS[node], abs=1e-3)
    junction = tpl.solid.Moved(occ.translation(pos))
    lo, hi = (np.asarray(x) for x in occ.bounding_box(junction))
    pts = lo + np.random.default_rng(5).random((1500, 3)) * (hi - lo)
    inside = {}
    for name, shape in (("junction", junction), ("block", block), ("trim", solid)):
        c = BRepClass3d_SolidClassifier(shape)
        got = []
        for p in pts:
            c.Perform(gp_Pnt(*p), 1e-7)
            got.append(c.State() == TopAbs_State.TopAbs_IN)
        inside[name] = np.array(got)
    want = inside["junction"] & inside["block"]
    assert want.sum() > 100 and (inside["junction"] & ~inside["block"]).sum() > 100
    assert np.array_equal(inside["trim"], want)


def test_a_junction_the_exact_trim_handles_never_climbs_the_ladder(setup, monkeypatch):
    """The ladder runs only where the probe has already found material outside:
    a sound junction pays nothing and keeps the exact trim."""
    lp, tpl, block, _mesh, probe = setup
    calls = []
    real = boundary._trim_fuzzy
    monkeypatch.setattr(boundary, "_trim_fuzzy",
                        lambda *a: calls.append(a) or real(*a))
    # A node of the block well inside it, trimmed by the box faces only.
    trim = trim_junction(lp, tpl, _at(lp, (196, -267, -272)), block, probe)
    assert trim.pieces and not trim.dropped and trim.fuzzy == 0.0
    assert calls == []


def test_a_rung_the_kernel_refuses_leaves_the_junction_dropped(setup, monkeypatch):
    """No rung that succeeds means the old outcome, the drop and its report --
    never a junction kept with material outside."""
    lp, tpl, block, _mesh, probe = setup
    monkeypatch.setattr(boundary, "_trim_fuzzy", lambda *a: None)
    trim = trim_junction(lp, tpl, _at(lp, (197, -267, -274)), block, probe)
    assert trim.dropped and trim.pieces == [] and trim.outside_mm > 1.0


def test_a_kernel_throw_is_a_refused_rung(setup, monkeypatch):
    lp, tpl, block, _mesh, _probe = setup

    class Boom:
        def __init__(self, *a):
            raise RuntimeError("Standard_Failure")

    monkeypatch.setattr(boundary, "BRepAlgoAPI_Common", Boom)
    assert boundary._trim_fuzzy(tpl, _at(lp, (197, -267, -274)), block, 1e-4) is None


@pytest.mark.parametrize("workers", [1, 2])
def test_the_count_and_the_nodes_reach_the_result(setup, tmp_path, workers):
    """On the sequential path and through the workers alike."""
    lp, tpl, block, mesh, _probe = setup
    body_path = str(tmp_path / "input.brep")
    BRepTools.Write_s(block, body_path)
    mesh_path = stage_mesh(mesh, str(tmp_path))
    nodes = np.array(list(JUNCTIONS), dtype=np.int64)
    args = (lp, tpl, nodes, block, body_path, str(tmp_path))
    kw = dict(mesh_path=mesh_path, outside_margin=2.0 * mesh.deviation)
    if workers > 1:
        with WorkerPool(workers) as pool:
            result = trim_boundary(*args, workers=workers, pool=pool, **kw)
    else:
        result = trim_boundary(*args, workers=1, **kw)
    assert result.n_dropped_junctions == 0
    assert result.n_fuzzy_junctions == 3
    assert sorted(n for n, _fz in result.fuzzy_nodes) == sorted(JUNCTIONS)
    assert all(fz in boundary.JUNCTION_FUZZY_LADDER for _n, fz in result.fuzzy_nodes)
