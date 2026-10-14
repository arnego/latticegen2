"""Print support (specification.md §4.6, docs/algorithm.md §14).

Three layers, cheapest first: overhang detection against shapes whose answer
is known in closed form, the support stage's contract and its merge checked on
their own, and
whole runs on a small box that go through the real kernel and the real
pipeline -- including one with a stand-in generator, which pins what the
pipeline hands the generator and that what it returns is what gets validated
and written.
"""

from __future__ import annotations

import os
import subprocess
import sys

import numpy as np
import pytest

from latticegen2.lattice import lattice_params


def _meshed(shape, cc=10.0, t=1.5):
    from latticegen2.classify import tessellate_surface

    lp = lattice_params(cc, t)
    tessellate_surface(shape, lp)
    return lp


def test_box_overhang_is_its_lid_and_nothing_else():
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox

    from latticegen2.support import overhang

    box = BRepPrimAPI_MakeBox(30.0, 20.0, 10.0).Shape()
    lp = _meshed(box)
    got = overhang.detect(box, lp.b, 60.0, 0.5)
    assert got.area == pytest.approx(600.0, rel=1e-9)
    assert np.allclose(got.points[:, 2], 10.0)            # never the underside
    assert np.allclose(got.inward, [0.0, 0.0, -1.0])
    assert got.counts()["edge"] == got.counts()["vertex"] == 0


@pytest.mark.parametrize("overhang_deg", [45.0, 60.0])
def test_sphere_overhang_matches_the_closed_form_cap(overhang_deg):
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeSphere

    from latticegen2.support import overhang

    radius = 30.0
    ball = BRepPrimAPI_MakeSphere(radius).Shape()
    lp = _meshed(ball)
    got = overhang.detect(ball, lp.b, overhang_deg, 0.5, margin_deg=0.0)
    # n.b > sin(phi) is the cap within (90 - phi) of the pole.
    exact = 2.0 * np.pi * radius ** 2 * (1.0 - np.sin(np.radians(overhang_deg)))
    assert got.area == pytest.approx(exact, rel=0.03)
    # Every sample really is on the surface and over the bar: exact normals.
    assert np.allclose(np.linalg.norm(got.points, axis=1), radius, atol=1e-9)
    assert (got.points[:, 2] / radius >= np.sin(np.radians(overhang_deg)) - 1e-9).all()


def test_a_keel_edge_is_an_overhang_even_though_its_faces_are_not():
    """A V-notch in the lid: both faces steep, the ridge between them hanging."""
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox

    from latticegen2 import occ
    from latticegen2.support import overhang

    profile = np.array([[0, 0, 0], [20, 0, 0], [20, 0, 10], [12, 0, 10],
                        [10, 0, 6], [8, 0, 10], [0, 0, 10]], dtype=float)
    notched = occ.prism(occ.polygon_face(profile), np.array([0.0, 15.0, 0.0]))
    assert occ.volume(notched) == pytest.approx((200.0 - 8.0) * 15.0)
    lp = _meshed(notched)
    got = overhang.detect(notched, lp.b, 60.0, 0.5)
    keel = got.points[got.kind == overhang.EDGE]
    assert len(keel) > 10
    assert np.allclose(keel[:, 0], 10.0) and np.allclose(keel[:, 2], 6.0)
    assert keel[:, 1].min() == pytest.approx(0.0) and keel[:, 1].max() == pytest.approx(15.0)

    plain = BRepPrimAPI_MakeBox(20.0, 15.0, 10.0).Shape()
    lp = _meshed(plain)
    assert overhang.detect(plain, lp.b, 60.0, 0.5).counts()["edge"] == 0


def test_a_hanging_tip_is_an_overhang_vertex():
    """A steep cone of shell material hanging into the volume: only its tip overhangs."""
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCone
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

    from latticegen2 import occ
    from latticegen2.support import overhang

    box = BRepPrimAPI_MakeBox(20.0, 20.0, 12.0).Shape()
    cone = BRepPrimAPI_MakeCone(
        gp_Ax2(gp_Pnt(10.0, 10.0, 12.0), gp_Dir(0.0, 0.0, -1.0)), 3.0, 0.0, 6.0
    ).Shape()
    body = occ.solids(BRepAlgoAPI_Cut(box, cone).Shape())[0]
    lp = _meshed(body)
    got = overhang.detect(body, lp.b, 60.0, 0.5)
    tips = got.points[got.kind == overhang.VERTEX]
    assert len(tips) == 1
    assert np.allclose(tips[0], [10.0, 10.0, 6.0], atol=1e-6)


def test_a_surface_at_the_limit_is_sampled_but_not_required():
    """A 120 degree drill point in the lid is tilted exactly 60 degrees.

    It is inside the margin, so supports are placed with it in mind, and it is
    not over the limit, so nothing has to be proven about it.
    """
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCone
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

    from latticegen2 import occ
    from latticegen2.support import overhang

    box = BRepPrimAPI_MakeBox(20.0, 20.0, 10.0).Shape()
    radius = 4.0
    cone = BRepPrimAPI_MakeCone(
        gp_Ax2(gp_Pnt(10.0, 10.0, 10.0), gp_Dir(0.0, 0.0, 1.0)),
        radius, 0.0, radius / np.tan(np.radians(60.0))).Shape()
    body = occ.solids(BRepAlgoAPI_Fuse(box, cone).Shape())[0]
    lp = _meshed(body)
    got = overhang.detect(body, lp.b, 60.0, 0.5)
    face = got.kind == overhang.FACE
    on_cone = face & (got.points[:, 2] > 10.0 + 1e-6)
    on_lid = face & (np.abs(got.points[:, 2] - 10.0) <= 1e-6)
    assert on_cone.sum() > 20 and on_lid.sum() > 100
    assert not got.required[on_cone].any()
    assert got.required[on_lid].all()
    # One degree stricter and the same cone is over the limit.
    assert overhang.detect(body, lp.b, 59.0, 0.5).required[on_cone].all()


# --- the support stage's contract ---------------------------------------------


def _box(lo=(0.0, 0.0, 0.0), size=(5.0, 5.0, 5.0)):
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt

    from latticegen2 import occ

    return occ.solids(BRepPrimAPI_MakeBox(gp_Pnt(*lo), *size).Shape())[0]


def test_a_result_must_be_support_solids_and_never_the_lattice():
    from latticegen2 import occ
    from latticegen2.errors import ProcessingError
    from latticegen2.support import SupportResult, check_result

    lattice = [_box(), _box((10.0, 0.0, 0.0))]
    check_result(SupportResult(supports=[]), lattice)
    check_result(SupportResult(supports=[_box((20.0, 0.0, 0.0))]), lattice)
    with pytest.raises(ProcessingError, match="other than a solid"):
        check_result(SupportResult(supports=[occ.faces(lattice[0])[0]]), lattice)
    with pytest.raises(ProcessingError, match="lattice solid as support"):
        check_result(SupportResult(supports=[lattice[1]]), lattice)


def _no_unify(solids):
    return list(solids), {"volume_drift": 0.0}


def test_supports_are_fused_into_the_lattice_they_touch_and_nothing_else():
    from latticegen2 import occ
    from latticegen2.support import merge

    touched, alone = _box(), _box((10.0, 0.0, 0.0))
    reaching = _box((4.0, 1.0, 1.0), (2.0, 2.0, 2.0))     # 1 mm into `touched`
    wall_attached = _box((20.0, 0.0, 0.0), (1.0, 1.0, 1.0))
    got = merge.merge_with_lattice([touched, alone], [reaching, wall_attached],
                                   unify=_no_unify)
    assert any(s.IsSame(alone) for s in got.solids)       # untouched: same object
    assert not any(s.IsSame(touched) for s in got.solids)
    assert len(got.solids) == 3                            # merged, alone, wall-attached
    assert got.lattice_volume == pytest.approx(250.0)
    assert got.support_volume == pytest.approx(9.0)
    assert got.output_volume == pytest.approx(250.0 + 4.0 + 1.0)
    assert sorted(round(occ.volume(s), 6) for s in got.solids) == [1.0, 125.0, 129.0]


def test_a_merge_with_no_supports_passes_the_lattice_through_unmeasured():
    from latticegen2.support import merge

    lattice = [_box(), _box((10.0, 0.0, 0.0))]
    got = merge.merge_with_lattice(lattice, [], unify=_no_unify)
    assert all(a.IsSame(b) for a, b in zip(got.solids, lattice))
    assert got.lattice_volume is None and got.output_volume is None


def test_a_merge_that_loses_lattice_fails_rather_than_shipping_less(monkeypatch):
    from latticegen2.errors import ProcessingError
    from latticegen2.support import merge

    lattice = [_box()]
    support = _box((4.0, 1.0, 1.0), (2.0, 2.0, 2.0))
    # A fuse that reports success and returns only the support -- the shape of
    # the kernel failure docs/algorithm.md §7.2 records.
    monkeypatch.setattr(merge, "_fuse", lambda members, where: [members[-1]])
    with pytest.raises(ProcessingError, match="lost lattice"):
        merge.merge_with_lattice(lattice, [support], unify=_no_unify)


def test_volume_bounds_refuse_both_loss_and_invented_material():
    from latticegen2.errors import ProcessingError
    from latticegen2.support.merge import check_volumes

    check_volumes(100.0, 10.0, 100.0)
    check_volumes(100.0, 10.0, 110.0)
    with pytest.raises(ProcessingError, match="lost lattice"):
        check_volumes(100.0, 10.0, 99.0)
    with pytest.raises(ProcessingError, match="gained material"):
        check_volumes(100.0, 10.0, 111.0)


# --- the whole pipeline -------------------------------------------------------


@pytest.fixture(scope="module")
def box_step(tmp_path_factory):
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt

    from latticegen2 import occ

    path = tmp_path_factory.mktemp("support") / "box.step"
    occ.write_step(BRepPrimAPI_MakeBox(gp_Pnt(1.3, 0.7, 0.4), 26.0, 23.0, 19.0).Shape(),
                   str(path), "box")
    return str(path)


def _run(args):
    main = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "src", "main.py")
    return subprocess.run([sys.executable, main, *args], capture_output=True, text=True)


def test_a_supported_box_runs_through_the_stub_and_writes_the_lattice(box_step, tmp_path):
    """Until the generator exists: overhangs found, a warning, the lattice unchanged."""
    from latticegen2 import occ
    from latticegen2.runlog import stages_for
    from latticegen2.support import NOT_IMPLEMENTED_WARNING

    plain = str(tmp_path / "plain.step")
    held = str(tmp_path / "held.step")
    a = _run(["-i", box_step, "-cc", "10", "-t", "1.5", "-o", plain, "--cores", "1"])
    b = _run(["-i", box_step, "-cc", "10", "-t", "1.5", "-o", held, "--cores", "1",
              "--support"])
    assert a.returncode == 0, a.stderr
    assert b.returncode == 0, b.stderr
    assert NOT_IMPLEMENTED_WARNING in b.stdout

    log = open(held[:-5] + ".log", encoding="utf-8").read()
    order = [log.index(f"stage {name}") for name in
             ("simplify", "overhang", "support", "validate", "export")]
    assert order == sorted(order)
    assert stages_for(True).index("overhang") == stages_for(True).index("simplify") + 1
    assert "stage overhang" not in open(plain[:-5] + ".log", encoding="utf-8").read()
    assert int(log.split("overhang_samples: ")[1].split()[0]) > 1000
    assert NOT_IMPLEMENTED_WARNING in log
    assert "support_volume_mm3: 0.0" in log
    assert "lattice_volume_mm3" in log and "output_volume_mm3" in log

    lattice = occ.solids(occ.read_step(plain))
    supported = occ.solids(occ.read_step(held))
    assert len(supported) == len(lattice)
    assert sum(occ.volume(s) for s in supported) == pytest.approx(
        sum(occ.volume(s) for s in lattice), rel=1e-9)


def test_the_generator_gets_the_completed_lattice_and_its_result_is_what_ships(
        box_step, tmp_path, monkeypatch):
    """In process, with a stand-in generator: what goes in, and what comes out.

    The stand-in returns one cube straddling a lattice node. That it reaches
    the written file merged into the lattice, with the volumes in the summary
    adding up, is what proves the support stage merges the generator's result
    and that the product -- not the lattice it was given -- is gated and
    exported.
    """
    from latticegen2 import occ
    from latticegen2.cli import parse_args
    from latticegen2.pipeline import run_pipeline
    from latticegen2.runlog import RunLog, stages_for
    from latticegen2.support import SupportResult
    from latticegen2.support import generate as gen

    seen = {}

    def stand_in(inp):
        from latticegen2.lattice import nodes

        seen["inp"] = inp
        # Asked now: the temp folder is removed once the run succeeds.
        seen["staged"] = os.path.isfile(inp.body_path) and os.path.isfile(inp.mesh_path)
        # A 3 mm cube on the kept node nearest the body's middle: partly inside
        # lattice material and partly in the void around it, so the merge has
        # something to fuse and something to add.
        pos = nodes(inp.lp, inp.kept_nodes)
        centre = pos[np.argmin(np.linalg.norm(pos - [14.3, 12.2, 9.9], axis=1))]
        return SupportResult(supports=[_box(tuple(centre - 1.5), (3.0, 3.0, 3.0))],
                             stats={"support_stand_in": 1},
                             notes=["note: stand-in generator ran"])

    monkeypatch.setattr(gen, "generate", stand_in)
    out = str(tmp_path / "held.step")
    args = parse_args(["-i", box_step, "-cc", "10", "-t", "1.5", "-o", out,
                       "--cores", "1", "--support", "--orient", "180", "0", "0",
                       "--overhang", "50", "--support-thickness", "1.2"])
    rl = RunLog(args.log_path).open()
    rl.stage_names = stages_for(True)
    try:
        stats = run_pipeline(args, rl)
    finally:
        rl.close()

    inp = seen["inp"]
    assert len(inp.lattice) >= 1
    assert all(occ.is_valid(s) for s in inp.lattice)
    assert inp.overhang_deg == 50.0 and inp.thickness_mm == 1.2
    assert inp.lp.orient == (180.0, 0.0, 0.0)
    assert np.allclose(inp.lp.b, [0.0, 0.0, -1.0], atol=1e-12)
    # Turned over, the box's bottom face is the ceiling: z = 0.4 in part terms.
    assert len(inp.samples) > 1000
    assert np.allclose(inp.samples.points[inp.samples.kind == 0][:, 2], 0.4)
    # Every kept junction once: the interior nodes and the boundary pieces' nodes.
    assert inp.kept_nodes.ndim == 2 and inp.kept_nodes.shape[1] == 3
    assert len({tuple(r) for r in inp.kept_nodes.tolist()}) == len(inp.kept_nodes)
    assert stats["nodes_interior"] < len(inp.kept_nodes) <= (
        stats["nodes_interior"] + stats["nodes_boundary"])
    assert seen["staged"]

    assert stats["support_solids_generated"] == 1 and stats["support_stand_in"] == 1
    assert stats["support_merge_groups"] == 1
    # The summary splits what ships by origin, and the parts add up.
    assert 0.0 < stats["support_volume_mm3"] < 27.0
    assert stats["lattice_volume_mm3"] + stats["support_volume_mm3"] == pytest.approx(
        stats["output_volume_mm3"], abs=1e-3)
    written = occ.solids(occ.read_step(out))
    assert len(written) == stats["solids_written"]
    assert sum(occ.volume(w) for w in written) == pytest.approx(
        stats["output_volume_mm3"], rel=1e-6)
    # Face and edge counts are retaken on the product, not left at the lattice's.
    assert stats["output_faces"] == sum(occ.count_subshapes(w)[0] for w in written)
    log = open(args.log_path, encoding="utf-8").read()
    assert "note: stand-in generator ran" in log
    assert "no support algorithm is implemented" not in log


def test_support_code_is_not_loaded_without_the_flag(box_step, tmp_path):
    """specification.md §4.6: a run without `--support` never imports it."""
    src = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src")
    out = str(tmp_path / "o.step")
    code = (
        "import sys; sys.path.insert(0, r'%s'); "
        "from latticegen2.__main__ import main; "
        "rc = main(['-i', r'%s', '-cc', '10', '-t', '1.5', '-o', r'%s', "
        "'--cores', '1', '--orient', '10', '0', '0']); "
        "print('LOADED' if any(m.startswith('latticegen2.support') for m in sys.modules) "
        "else 'ABSENT', rc)" % (src, box_step, out)
    )
    proc = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert proc.stdout.strip().endswith("ABSENT 0"), proc.stdout + proc.stderr
