"""Pockets the supports seal off, filled (docs/algorithm.md §14.3 S5).

The fixture is a cup: a ring wall standing on the floor of a round chamber,
its top 3 mm below the roof. A lid from the roof down onto the ring seals the
cup's inside, which is a cylinder of known volume.
"""

from __future__ import annotations

import numpy as np
import pytest

from latticegen2.lattice import lattice_params

NO_NODES = np.empty((0, 3), dtype=np.int64)
R1, R2, R, H, TOP = 2.0, 6.0, 10.0, 8.0, 5.0
"""Cup radius, ring's outer radius, chamber radius, roof height, ring top."""


def _cyl(r, h, z0=0.0):
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

    return BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, z0), gp_Dir(0, 0, 1)), r, h).Shape()


@pytest.fixture(scope="module")
def cup():
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut

    from latticegen2 import occ
    from latticegen2.classify import tessellate_surface

    ring = occ.solids(BRepAlgoAPI_Cut(_cyl(R2, TOP), _cyl(R1, TOP)).Shape())[0]
    body = occ.solids(BRepAlgoAPI_Cut(_cyl(R, H), ring).Shape())[0]
    lp = lattice_params(10.0, 1.5, (0.0, 0.0, 0.0))
    return body, tessellate_surface(body, lp), lp


def test_a_lid_on_the_cup_fills_exactly_the_cup(cup):
    """The envelope starts 1.5 mm below the lid, where the cup's inside and the
    chamber around the ring both reach its faces; the cup's inside, the smaller,
    reaches only the lower one, so the envelope is widened past it once."""
    from latticegen2 import occ
    from latticegen2.support import voids

    body, mesh, lp = cup
    lid = _cyl(R2, H - TOP, TOP)
    got = voids.find_pockets(body, mesh, [lid], lambda box: [], 1.5)
    assert got.stats["support_pocket_widenings"] == 1
    assert len(got.fills) == 1
    fill = got.fills[0]
    assert occ.volume(fill) == pytest.approx(np.pi * R1 ** 2 * TOP, rel=1e-6)
    assert occ.is_valid(fill)
    # Every face of the fill lies on the ring, the floor or the lid: none is free.
    from latticegen2.support.build import _mesh_of

    for _f, pts, _u, _t in _mesh_of(fill, 0.01):
        on_ring = np.all(np.abs(np.hypot(pts[:, 0], pts[:, 1]) - R1) < 1e-6)
        assert on_ring or np.all(np.abs(pts[:, 2]) < 1e-6) or np.all(np.abs(pts[:, 2] - TOP) < 1e-6)



def test_the_cup_found_in_slabs_is_the_cup(cup, monkeypatch):
    """Four slabs put a plane through the middle of the cup: the pieces either
    side of it are joined by their faces on it, so the pocket is found whole
    and filled exactly, in two pieces that add up to it."""
    from latticegen2 import occ
    from latticegen2.support import voids

    body, mesh, lp = cup
    monkeypatch.setattr(voids, "_slab_count", lambda *args: 4)
    lid = _cyl(R2, H - TOP, TOP)
    got = voids.find_pockets(body, mesh, [lid], lambda box: [], 1.5)
    assert got.stats["support_pocket_slabs"] == 4
    assert len(got.fills) == 2
    assert sum(occ.volume(f) for f in got.fills) == pytest.approx(np.pi * R1 ** 2 * TOP, rel=1e-6)



def test_a_refused_combined_cut_is_redone_in_two_steps(cup, monkeypatch):
    """The region cut by lattice and supports in one call can be refused --
    on the rehearsal part, a support beside junctions it touches -- and the
    lattice-free space is then cut by the supports alone, to the same pocket."""
    from latticegen2 import occ
    from latticegen2.support import voids

    body, mesh, lp = cup
    real = voids._free
    refused = []

    def refuse_combined(region, tools, removes, what, one_by_one=True):
        if what == "the lattice and the supports":
            refused.append(what)
            raise voids.CutRefused("refused")
        return real(region, tools, removes, what, one_by_one)

    monkeypatch.setattr(voids, "_free", refuse_combined)
    got = voids.find_pockets(body, mesh, [_cyl(R2, H - TOP, TOP)], lambda box: [], 1.5)
    assert refused
    assert sum(occ.volume(f) for f in got.fills) == pytest.approx(np.pi * R1 ** 2 * TOP, rel=1e-6)


def _refusing(voids, monkeypatch, refuse):
    """``voids._slab_free`` refusing every box ``refuse(lo, hi)`` names; the
    x extents of the boxes asked for are recorded."""
    real = voids._slab_free
    asked = []

    def slab_free(body, probe, lo, hi, *args):
        asked.append((float(lo[0]), float(hi[0])))
        if refuse(lo, hi):
            raise voids.CutRefused("refused")
        return real(body, probe, lo, hi, *args)

    monkeypatch.setattr(voids, "_slab_free", slab_free)
    return asked


def test_a_refused_slab_is_split_and_searched_again(cup, monkeypatch):
    """Four slabs, each refused while wider than 3 mm: every one is split in
    two, the halves are searched, and the pocket is found whole."""
    from latticegen2 import occ
    from latticegen2.support import voids

    body, mesh, lp = cup
    monkeypatch.setattr(voids, "_slab_count", lambda *args: 4)
    asked = _refusing(voids, monkeypatch, lambda lo, hi: hi[0] - lo[0] > 3.0)
    got = voids.find_pockets(body, mesh, [_cyl(R2, H - TOP, TOP)], lambda box: [], 1.5)
    # Once before the envelope is widened past the cup and once after.
    assert len(asked) == 2 * (4 + 8)
    assert got.stats["support_pocket_slabs"] == 8
    assert got.unsearched == []
    assert sum(occ.volume(f) for f in got.fills) == pytest.approx(np.pi * R1 ** 2 * TOP, rel=1e-6)


def test_a_slab_refused_at_every_split_is_left_open(cup, monkeypatch):
    """The slabs through the cup's middle refused however finely they are
    split: the free space beside them counts as open, so nothing is filled,
    and the slabs left are named."""
    from latticegen2.support import voids

    body, mesh, lp = cup
    monkeypatch.setattr(voids, "_slab_count", lambda *args: 4)
    _refusing(voids, monkeypatch, lambda lo, hi: lo[0] < 0.5 and hi[0] > -0.5)
    got = voids.find_pockets(body, mesh, [_cyl(R2, H - TOP, TOP)], lambda box: [], 1.5)
    assert got.fills == []
    assert got.unsearched and got.stats["support_pocket_slabs_unsearched"] == len(got.unsearched)
    # Whole coordinates, so the region can be found in the part; the kernel's
    # reason is for the log.
    for line in got.unsearched:
        assert all(f"{axis} " in line for axis in "xyz") and "refused" not in line
    assert len(got.refusals) == len(got.unsearched) and "refused" in got.refusals[0]


def test_the_slabs_run_across_a_worker_pool(cup, monkeypatch, tmp_path):
    """The same four slabs, each in a worker reading the staged body and mesh,
    give the same fills as in this process."""
    from latticegen2 import occ
    from latticegen2.classify import save_mesh
    from latticegen2.parallel import WorkerPool, write_brep
    from latticegen2.support import voids

    body, mesh, lp = cup
    monkeypatch.setattr(voids, "_slab_count", lambda *args: 4)
    write_brep(body, str(tmp_path / "body.brep"))
    save_mesh(mesh, str(tmp_path / "mesh.npz"))
    staged = {"body_path": str(tmp_path / "body.brep"), "mesh_path": str(tmp_path / "mesh.npz"),
              "lp_args": (lp.cc, lp.t, tuple(lp.orient)), "kept_nodes": NO_NODES,
              "tmpdir": str(tmp_path)}
    lid = _cyl(R2, H - TOP, TOP)
    with WorkerPool(2) as pool:
        got = voids.find_pockets(body, mesh, [lid], lambda box: [], 1.5,
                                 pool=pool, staged=staged)
    assert len(got.fills) == 2
    assert sum(occ.volume(f) for f in got.fills) == pytest.approx(np.pi * R1 ** 2 * TOP, rel=1e-6)

def test_a_pocket_the_lattice_already_sealed_is_left_alone(cup):
    from latticegen2.support import voids

    body, mesh, lp = cup
    plug = _cyl(R1, 0.5, TOP - 0.5)        # lattice in the cup's mouth
    lid = _cyl(R2, H - TOP, TOP)
    got = voids.find_pockets(body, mesh, [lid], lambda box: [plug], 1.5)
    assert got.fills == []


def test_a_cell_that_seals_nothing_fills_nothing(cup):
    from latticegen2 import occ
    from latticegen2.support import voids

    body, mesh, lp = cup
    post = _cyl(1.0, H - TOP, TOP).Moved(occ.translation(np.array([8.0, 0.0, 0.0])))
    got = voids.find_pockets(body, mesh, [post], lambda box: [], 1.5)
    assert got.fills == []
    assert got.stats["support_pocket_widenings"] == 0


def test_another_groups_support_in_the_envelope_is_not_free_space(cup):
    """A support of another group standing in the cup is support, not powder:
    the fill is the cup less that support, and overlaps it nowhere."""
    from latticegen2 import occ
    from latticegen2.support import voids

    body, mesh, lp = cup
    other = _cyl(1.0, 2.0, 0.0)            # on the cup's floor
    lid = _cyl(R2, H - TOP, TOP)
    got = voids.find_pockets(body, mesh, [lid], lambda box: [], 1.5, others=[other])
    assert len(got.fills) == 1
    assert occ.volume(got.fills[0]) == pytest.approx(
        np.pi * (R1 ** 2 * TOP - 1.0 ** 2 * 2.0), rel=1e-6)
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Common

    common = occ.solids(BRepAlgoAPI_Common(got.fills[0], other).Shape())
    assert sum(occ.volume(s) for s in common) < 1e-6


def test_a_witness_the_cut_disagrees_with_fails_loudly(cup, monkeypatch):
    from latticegen2.errors import ProcessingError
    from latticegen2.support import voids

    body, mesh, lp = cup
    monkeypatch.setattr(voids, "witness", lambda piece: np.array([0.0, 0.0, -5.0]))
    with pytest.raises(ProcessingError, match="disagree"):
        voids.find_pockets(body, mesh, [_cyl(R2, H - TOP, TOP)], lambda box: [], 1.5)


def test_the_generator_leaves_no_pocket_sealed_on_the_cup(cup, tmp_path):
    """Whatever the frustum supports seal off in the cup is filled and joined to
    them: asked again of the result, the pocket finder finds nothing."""
    from latticegen2 import occ
    from latticegen2.support import assemble, join, voids
    from latticegen2.support.generate import supports_for

    body, mesh, lp = cup
    solids, stats, _lines = supports_for(body, mesh, [], NO_NODES, lp, 1.5,
                                         tmpdir=str(tmp_path))
    assert solids and len(solids) == stats["support_bodies"]
    assert all(occ.is_valid(s) for s in solids)
    for g in join.groups(solids):
        again = voids.find_pockets(body, mesh, [solids[i] for i in g],
                                   lambda box: [], 1.5)
        assert again.fills == []
    assemble.assemble([], solids, lp, NO_NODES)
