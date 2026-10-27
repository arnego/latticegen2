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
    from latticegen2 import occ
    from latticegen2.support import voids

    body, mesh, lp = cup
    lid = _cyl(R2, H - TOP, TOP)
    got = voids.find_pockets(body, mesh, [lid], lambda box: [], 1.5, lp.b)
    assert len(got.fills) == 1
    fill = got.fills[0]
    assert occ.volume(fill) == pytest.approx(np.pi * R1 ** 2 * TOP, rel=1e-6)
    assert occ.is_valid(fill)
    # Every face of the fill lies on the ring, the floor or the lid: none is free.
    from latticegen2.support.build import _mesh_of

    for _f, pts, _u, _t in _mesh_of(fill, 0.01):
        on_ring = np.all(np.abs(np.hypot(pts[:, 0], pts[:, 1]) - R1) < 1e-6)
        assert on_ring or np.all(np.abs(pts[:, 2]) < 1e-6) or np.all(np.abs(pts[:, 2] - TOP) < 1e-6)


def test_a_pocket_the_lattice_already_sealed_is_left_alone(cup):
    from latticegen2.support import voids

    body, mesh, lp = cup
    plug = _cyl(R1, 0.5, TOP - 0.5)        # lattice in the cup's mouth
    lid = _cyl(R2, H - TOP, TOP)
    got = voids.find_pockets(body, mesh, [lid], lambda box: [plug], 1.5, lp.b)
    assert got.fills == []


def test_a_cell_that_seals_nothing_fills_nothing(cup):
    from latticegen2 import occ
    from latticegen2.support import voids

    body, mesh, lp = cup
    post = _cyl(1.0, H - TOP, TOP).Moved(occ.translation(np.array([8.0, 0.0, 0.0])))
    got = voids.find_pockets(body, mesh, [post], lambda box: [], 1.5, lp.b)
    assert got.fills == []


def test_a_witness_the_cut_disagrees_with_fails_loudly(cup, monkeypatch):
    from latticegen2.errors import ProcessingError
    from latticegen2.support import voids

    body, mesh, lp = cup
    monkeypatch.setattr(voids, "witness", lambda piece: np.array([0.0, 0.0, -5.0]))
    with pytest.raises(ProcessingError, match="disagree"):
        voids.find_pockets(body, mesh, [_cyl(R2, H - TOP, TOP)], lambda box: [], 1.5, lp.b)


def test_the_generator_seals_the_cup_and_fills_it(cup):
    """Rays from above the cup's mouth converge on its rim, so that patch is a
    converging cell of its own, closed by the cone of its rays (§14.3 S3b, S4);
    it seals the cup, which S5 fills."""
    from latticegen2 import occ
    from latticegen2.support import assemble
    from latticegen2.support.generate import supports_for

    body, mesh, lp = cup
    solids, stats, lines = supports_for(body, mesh, [], NO_NODES, lp, 1.5)
    assert "'converging_cells': 1" in lines[0]
    assert stats["support_pockets"] >= 1
    assert len(solids) == stats["support_cells"] + stats["support_pockets"]
    fills = solids[stats["support_cells"]:]
    cup_fill = [f for f in fills if np.hypot(*occ.centroid(f)[:2]) < 0.1
                and occ.centroid(f)[2] < TOP]
    assert len(cup_fill) == 1
    # The cup, up to the lid's cone, which meets the ring a little below its top.
    v = occ.volume(cup_fill[0])
    assert np.pi * R1 ** 2 * (TOP - 1.0) <= v <= np.pi * R1 ** 2 * H
    # Nothing else is sealed: the wall cells reach the roof's crease, so no
    # ring is left there for S5 to fill (§14.3 S3d).
    assert stats["support_pockets"] == 1
    assert all(occ.is_valid(s) for s in solids)
    assemble.assemble([], solids, lp, NO_NODES)


def test_supports_on_the_ledge_lean_in_and_keep_a_footing(cup):
    """The ledge's cells stand on a level face: their outer sides lean in from
    the parting toward the cup, and stop where the footing on the ledge would
    narrow below t_s (§14.3 S3d)."""
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.TopAbs import TopAbs_State
    from OCP.gp import gp_Pnt

    from latticegen2 import occ
    from latticegen2.support.generate import supports_for
    from latticegen2.support.params import SUPPORT_OVERHANG_DEG

    body, mesh, lp = cup
    solids, stats, lines = supports_for(body, mesh, [], NO_NODES, lp, 1.5)
    region = lines[0]
    lean = float(region.split("'lean_max_deg': ")[1].split("}")[0].split(",")[0])
    assert 30.0 < lean <= SUPPORT_OVERHANG_DEG + 1e-6
    cells = solids[:stats["support_cells"]]
    tests = [BRepClass3d_SolidClassifier(s) for s in cells]

    def held(p):
        for t in tests:
            t.Perform(gp_Pnt(*map(float, p)), 1e-7)
            if t.State() == TopAbs_State.TopAbs_IN:
                return True
        return False

    z = TOP + 0.02
    rs = np.arange(R1 + 0.01, R2, 0.02)
    outer = []
    for a in np.radians(np.arange(0.0, 360.0, 15.0)):
        on = [r for r in rs if held((r * np.cos(a), r * np.sin(a), z))]
        assert on, f"nothing stands on the ledge at {np.degrees(a):.0f} deg"
        # The footing: the material on the ledge, from the cup's rim outward.
        assert on[-1] - on[0] >= 1.5 - 0.05
        outer.append(on[-1])
    # Leaning in: the footing ends well inside the parting at r = 5.76, where
    # it stood before -- at the cuts between the ledge's cells too, which lean
    # once for both cells (§14.3 S3d). They stood vertical, and the footing
    # reached out to the parting in four lobes.
    assert max(outer) < 4.5


def test_cuts_lean_once_for_both_cells(cup):
    """A cut is an inner face of one support split in two: both cells carry the
    same ray at every point of it, and its ends lean as far as the parting they
    end on (§14.3 S3d). They kept a vertical ray, the lean beside them tapered
    up from it, and the footing reached the parting in four lobes."""
    from latticegen2.support import build, cells, generate, lean, obstacles, rays, regions

    body, mesh, lp = cup
    rm = regions.build(body, lp.b)
    reg = regions.regions(rm)[0]
    ob = obstacles.from_shapes(body, [], deflection=generate.OBSTACLE_DEFLECTION)
    delta = generate.overshoot_for(1.5, mesh.deviation)
    got = rays.first_contact(rm.verts[reg.verts] - rays.APEX_NUDGE * rm.normals[reg.verts],
                             lp.b, ob)
    part = cells.partition(rm, reg, got, ob, lp.b, 1.5,
                           inside_body=build.Gates(mesh, lp.b).inside_body, overshoot=delta)
    before = {c.id: c.loops[0].targets.copy() for c in part.cells}
    lean.lean_cells(part.cells, ob, lp.b, 1.5, delta)
    shared = {}
    for c in part.cells:
        L = c.loops[0]
        for k, (pid, kd) in enumerate(zip(L.keys, L.kind)):
            if kd in (cells.CUT, cells.TRIPLE):
                shared.setdefault((pid, c.key), []).append(L.targets[k])
    cut_pairs = [v for (pid, _key), v in shared.items() if len(v) == 2]
    assert len(cut_pairs) > 20
    for a, b_ in cut_pairs:
        assert np.array_equal(a, b_)
    # Where a cut meets the ledge's parting, its end leans with the parting.
    n_ends = 0
    for c in part.cells:
        L = c.loops[0]
        n = len(L.kind)
        s = np.linalg.norm(L.targets - before[c.id], axis=1)
        if s.max() < 1.0:
            continue                       # the wall annulus stands on walls
        for k in np.nonzero(L.kind == cells.TRIPLE)[0]:
            prev, nxt = L.kind[(k - 1) % n], L.kind[(k + 1) % n]
            if {int(prev), int(nxt)} == {cells.CUT, cells.PARTING}:
                beside = (k - 1) % n if prev == cells.PARTING else (k + 1) % n
                assert s[beside] > 1.0
                assert s[k] == pytest.approx(s[beside], rel=0.1)
                n_ends += 1
    assert n_ends == 8
