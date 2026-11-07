"""Print support (specification.md §4.6, docs/algorithm.md §14).

Three layers, cheapest first: overhang detection against shapes whose answer
is known in closed form, the support stage's contract and its assembly checked
on their own, and whole runs on a small box that go through the real kernel and
the real pipeline -- including one with a stand-in generator, which pins what
the pipeline hands the generator and that what it returns is what gets
validated and written, as separate bodies beside the lattice.
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
    got = overhang.detect(box, lp.b, 0.5)
    assert got.area == pytest.approx(600.0, rel=1e-9)
    assert np.allclose(got.points[:, 2], 10.0)            # never the underside
    assert np.allclose(got.inward, [0.0, 0.0, -1.0])
    assert got.counts()["edge"] == got.counts()["vertex"] == 0


def test_sphere_overhang_matches_the_closed_form_cap():
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeSphere

    from latticegen2.support import overhang

    radius = 30.0
    ball = BRepPrimAPI_MakeSphere(radius).Shape()
    lp = _meshed(ball)
    got = overhang.detect(ball, lp.b, 0.5)
    # n.b > sin 60 deg is the cap within 30 deg of the pole.
    exact = 2.0 * np.pi * radius ** 2 * (1.0 - np.sin(np.radians(60.0)))
    assert got.area == pytest.approx(exact, rel=0.03)
    # Every sample really is on the surface and over the bar: exact normals.
    assert np.allclose(np.linalg.norm(got.points, axis=1), radius, atol=1e-9)
    assert (got.points[:, 2] / radius > np.sin(np.radians(60.0))).all()


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
    got = overhang.detect(notched, lp.b, 0.5)
    keel = got.points[got.kind == overhang.EDGE]
    assert len(keel) > 10
    assert np.allclose(keel[:, 0], 10.0) and np.allclose(keel[:, 2], 6.0)
    assert keel[:, 1].min() == pytest.approx(0.0) and keel[:, 1].max() == pytest.approx(15.0)

    plain = BRepPrimAPI_MakeBox(20.0, 15.0, 10.0).Shape()
    lp = _meshed(plain)
    assert overhang.detect(plain, lp.b, 0.5).counts()["edge"] == 0


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
    got = overhang.detect(body, lp.b, 0.5)
    tips = got.points[got.kind == overhang.VERTEX]
    assert len(tips) == 1
    assert np.allclose(tips[0], [10.0, 10.0, 6.0], atol=1e-6)


def test_a_surface_at_sixty_degrees_is_not_supported():
    """A 120 degree drill point in the lid is tilted exactly 60 degrees.

    The bar is strict (docs/algorithm.md §14.1), so it is not detected at all,
    while the flat lid around it is. A point a degree flatter would be.
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
    got = overhang.detect(body, lp.b, 0.5)
    face = got.kind == overhang.FACE
    on_cone = face & (got.points[:, 2] > 10.0 + 1e-6)
    on_lid = face & (np.abs(got.points[:, 2] - 10.0) <= 1e-6)
    assert on_cone.sum() == 0 and on_lid.sum() > 100
    # The lid's area is the box's lid less the cone's footprint.
    assert got.area == pytest.approx(400.0 - np.pi * radius ** 2, rel=0.02)


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


_NO_NODES = np.empty((0, 3), dtype=np.int64)


def test_supports_ship_beside_the_lattice_as_separate_bodies():
    """Nothing is fused: the lattice comes through as the same objects, and
    supports touching it face to face -- or touching each other -- are bodies
    of their own, with additive volumes."""
    from latticegen2.support import assemble

    lp = lattice_params(10.0, 1.5)
    lattice = [_box(), _box((10.0, 0.0, 0.0))]
    on_top = _box((0.0, 0.0, 5.0), (5.0, 5.0, 2.0))       # shares the lid face
    beside = _box((5.0, 0.0, 5.0), (2.0, 5.0, 2.0))       # shares a face with on_top
    got = assemble.assemble(lattice, [on_top, beside], lp, _NO_NODES)
    assert len(got.solids) == 4
    assert all(a.IsSame(b) for a, b in zip(got.solids[:2], lattice))
    assert got.solids[2].IsSame(on_top) and got.solids[3].IsSame(beside)
    assert got.lattice_volume == pytest.approx(250.0)
    assert got.support_volume == pytest.approx(50.0 + 20.0)
    assert got.stats["support_bodies"] == 2 and got.stats["support_pairs_checked"] == 1


def test_no_supports_pass_the_lattice_through_unmeasured():
    from latticegen2.support import assemble

    lattice = [_box(), _box((10.0, 0.0, 0.0))]
    got = assemble.assemble(lattice, [], lattice_params(10.0, 1.5), _NO_NODES)
    assert all(a.IsSame(b) for a, b in zip(got.solids, lattice))
    assert got.lattice_volume is None


def test_a_support_overlapping_another_fails_rather_than_shipping():
    from latticegen2.errors import ProcessingError
    from latticegen2.support import assemble

    a = _box((0.0, 0.0, 5.0), (5.0, 5.0, 2.0))
    b = _box((4.0, 0.0, 5.0), (2.0, 5.0, 2.0))            # 1 mm into `a`
    with pytest.raises(ProcessingError, match="overlap by"):
        assemble.assemble([_box()], [a, b], lattice_params(10.0, 1.5), _NO_NODES)


def test_a_support_overlapping_the_lattice_fails_rather_than_shipping():
    """Checked against the junctions at the kept nodes, and then confirmed
    against the lattice solids themselves before it is believed."""
    from latticegen2.errors import ProcessingError
    from latticegen2.junction import build_template
    from latticegen2.lattice import nodes
    from latticegen2 import occ
    from latticegen2.support import assemble

    lp = lattice_params(10.0, 1.5)
    kept = np.array([[0, 0, 0]], dtype=np.int64)
    junction = build_template(lp).solid.Moved(occ.translation(nodes(lp, kept)[0]))
    clear = _box((20.0, 20.0, 20.0), (1.0, 1.0, 1.0))
    into = _box((-0.5, -0.5, -0.5), (1.0, 1.0, 1.0))      # the node region itself
    assemble.assemble([junction], [clear], lp, kept)
    with pytest.raises(ProcessingError, match="overlaps the lattice"):
        assemble.assemble([junction], [into], lp, kept)


def test_overlap_with_the_lattice_sums_over_junction_groups():
    """The junctions share no volume, so a support's overlap with all of them
    is the sum over any grouping of them: the check can be split into jobs."""
    from latticegen2.lattice import nodes
    from latticegen2.support.assemble import _LocalLattice, _Overlaps

    lp = lattice_params(10.0, 1.5)
    kept = np.array([[i, j, k] for i in range(3) for j in range(3) for k in range(3)],
                    dtype=np.int64)
    pos = nodes(lp, kept)
    lo = pos.min(axis=0) - 1.0
    support = _box(tuple(lo), tuple((pos.mean(axis=0) - lo).tolist()))
    local = _LocalLattice(lp, kept)
    idx = local.near_solid(support)
    groups = local.groups(idx, 2)
    chk = _Overlaps([support], local)
    whole = chk.run(("lattice", 0, idx))[0]
    parts = sum(chk.run(("lattice", 0, g))[0] for g in groups)
    assert len(groups) > 1 and sum(len(g) for g in groups) == len(idx)
    assert whole > 1.0
    assert parts == pytest.approx(whole, rel=1e-9)


def _fused(a, b):
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse

    from latticegen2 import occ

    return occ.solids(BRepAlgoAPI_Fuse(a, b).Shape())[0]


def test_a_junction_only_in_a_supports_box_is_not_asked():
    """An L-shaped support's box holds a junction no face of it comes near: it
    lies outside the support, and only the support's own faces decide that."""
    from latticegen2.support.assemble import _LocalLattice

    lp = lattice_params(10.0, 1.5)
    kept = np.zeros((1, 3), dtype=np.int64)          # one junction, at the origin
    ell = _fused(_box((-20.0, -20.0, -2.0), (10.0, 40.0, 4.0)),
                 _box((-20.0, 10.0, -2.0), (40.0, 10.0, 4.0)))
    local = _LocalLattice(lp, kept)
    assert len(local.near((np.array([-20.0, -20.0, -2.0]), np.array([20.0, 20.0, 2.0])))) == 1
    assert len(local.near_solid(ell)) == 0


def test_a_junction_wholly_inside_a_support_is_still_caught():
    """No face of the support comes near it, and it is asked anyway: its node
    lies inside the support."""
    from latticegen2.errors import ProcessingError
    from latticegen2.junction import build_template
    from latticegen2.support import assemble
    from latticegen2.support.assemble import _LocalLattice

    lp = lattice_params(10.0, 1.5)
    kept = np.zeros((1, 3), dtype=np.int64)
    around = _box((-10.0, -10.0, -10.0), (20.0, 20.0, 20.0))
    assert _LocalLattice(lp, kept).near_solid(around).tolist() == [0]
    with pytest.raises(ProcessingError, match="overlaps the lattice"):
        assemble.assemble([build_template(lp).solid], [around], lp, kept)


def test_the_overlap_check_reads_alike_across_the_worker_pool(tmp_path):
    """The same jobs in workers as in this process: the same refusal, word for
    word, for the lattice and for two supports."""
    from latticegen2 import occ
    from latticegen2.errors import ProcessingError
    from latticegen2.junction import build_template
    from latticegen2.lattice import nodes
    from latticegen2.parallel import WorkerPool
    from latticegen2.support import assemble

    lp = lattice_params(10.0, 1.5)
    kept = np.array([[i, j, 0] for i in range(3) for j in range(3)], dtype=np.int64)
    template = build_template(lp).solid
    lattice = [template.Moved(occ.translation(p)) for p in nodes(lp, kept)]
    into = _box((-0.5, -0.5, -0.5), (1.0, 1.0, 1.0))
    a = _box((40.0, 0.0, 5.0), (5.0, 5.0, 2.0))
    b = _box((44.0, 0.0, 5.0), (2.0, 5.0, 2.0))      # 1 mm into `a`

    def refusal(supports, pool=None):
        with pytest.raises(ProcessingError) as got:
            assemble.assemble(lattice, supports, lp, kept, pool=pool, tmpdir=str(tmp_path))
        return str(got.value)

    here = [refusal([into, a]), refusal([a, b])]
    with WorkerPool(2) as pool:
        there = [refusal([into, a], pool), refusal([a, b], pool)]
    assert "overlaps the lattice" in here[0] and "overlap by" in here[1]
    assert there == here


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
    # The run writes UTF-8 whatever the console's code page; decode it as such,
    # or the warning's section sign arrives as two characters.
    return subprocess.run([sys.executable, main, *args], capture_output=True, text=True,
                          encoding="utf-8", errors="replace")


def test_a_supported_box_ships_supports_beside_the_unchanged_lattice(box_step, tmp_path):
    """The real generator, end to end: overhangs found, support bodies written
    after the lattice, which is the same as a run without supports."""
    from latticegen2 import occ
    from latticegen2.runlog import stages_for
    from latticegen2.support import LATTICE_BACKS_ALL_NOTE

    plain = str(tmp_path / "plain.step")
    held = str(tmp_path / "held.step")
    a = _run(["-i", box_step, "-cc", "10", "-t", "1.5", "-o", plain, "--cores", "1"])
    b = _run(["-i", box_step, "-cc", "10", "-t", "1.5", "-o", held, "--cores", "1",
              "--support"])
    assert a.returncode == 0, a.stderr
    assert b.returncode == 0, b.stderr
    assert LATTICE_BACKS_ALL_NOTE not in b.stdout

    log = open(held[:-5] + ".log", encoding="utf-8").read()
    order = [log.index(f"stage {name}") for name in
             ("simplify", "overhang", "support", "pockets", "validate", "export")]
    assert order == sorted(order)
    assert stages_for(True).index("overhang") == stages_for(True).index("simplify") + 1
    assert "stage overhang" not in open(plain[:-5] + ".log", encoding="utf-8").read()
    assert int(log.split("overhang_samples: ")[1].split()[0]) > 1000
    n = int(log.split("supports: ")[1].split()[0])
    assert n >= 1
    support = float(log.split("support_volume_mm3: ")[1].split()[0])
    lattice_v = float(log.split("lattice_volume_mm3: ")[1].split()[0])
    output_v = float(log.split("output_volume_mm3: ")[1].split()[0])
    assert support > 0.0
    assert output_v == pytest.approx(lattice_v + support, rel=1e-6)

    lattice = occ.solids(occ.read_step(plain))
    supported = occ.solids(occ.read_step(held))
    assert len(supported) == len(lattice) + n
    assert sum(occ.volume(s) for s in supported[:len(lattice)]) == pytest.approx(
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
        seen["inp"] = inp
        # Asked now: the temp folder is removed once the run succeeds.
        seen["staged"] = os.path.isfile(inp.body_path) and os.path.isfile(inp.mesh_path)
        # A thin slab against the ceiling, which turned over is the box's
        # bottom face (z = 0.4), placed between the struts that reach it: inside
        # the body and overlapping no lattice, as a real support must be. It is
        # found by asking the lattice where it is not.
        from latticegen2 import occ as _occ

        slab = None
        for x in np.arange(2.0, 24.0, 0.25):
            for y in np.arange(1.5, 21.0, 0.25):
                trial = _box((x, y, 0.4), (0.5, 0.5, 0.3))
                if all(_occ.volume(_common(trial, l)) == 0.0 for l in inp.lattice):
                    slab = trial
                    break
            if slab is not None:
                break
        assert slab is not None
        seen["slab"] = slab
        return SupportResult(supports=[slab],
                             stats={"support_stand_in": 1},
                             notes=["note: stand-in generator ran"])

    monkeypatch.setattr(gen, "generate", stand_in)
    out = str(tmp_path / "held.step")
    args = parse_args(["-i", box_step, "-cc", "10", "-t", "1.5", "-o", out,
                       "--cores", "1", "--support", "--orient", "180", "0", "0",
                       "--support-thickness", "1.2"])
    rl = RunLog(args.log_path).open()
    rl.stage_names = stages_for(True)
    try:
        stats = run_pipeline(args, rl)
    finally:
        rl.close()

    inp = seen["inp"]
    assert len(inp.lattice) >= 1
    assert all(occ.is_valid(s) for s in inp.lattice)
    assert inp.thickness_mm == 1.2
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
    assert stats["support_bodies"] == 1
    # The summary splits what ships by origin, and the parts add up: the slab
    # is a body of its own, its whole volume added.
    assert stats["support_volume_mm3"] == pytest.approx(0.075, abs=1e-6)
    assert stats["solids_written"] == len(inp.lattice) + 1
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
    assert "returned no supports" not in log


def _common(a, b):
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Common

    return BRepAlgoAPI_Common(a, b).Shape()


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
