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
    assert all(occ.is_valid(s) for s in solids)
    assemble.assemble([], solids, lp, NO_NODES)
