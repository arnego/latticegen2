"""Support cell solids (docs/algorithm.md §14.3 S4), without a lattice.

The pillar in a round chamber has supports whose volume is known in closed
form; the hole-blocks are the parts this design was first inspected on.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from latticegen2.lattice import lattice_params

TESTDIR = os.path.dirname(os.path.abspath(__file__))
NO_NODES = np.empty((0, 3), dtype=np.int64)


def _supports(body, orient=(0.0, 0.0, 0.0), t=1.5):
    from latticegen2.classify import tessellate_surface
    from latticegen2.support.generate import supports_for

    lp = lattice_params(10.0, t, orient)
    mesh = tessellate_surface(body, lp)
    solids, stats, _lines = supports_for(body, mesh, [], NO_NODES, lp, t)
    return lp, solids, stats


def test_a_pillar_in_a_chamber_gets_the_closed_form_crowns():
    """Rays at the rod angle part midway (R = 5) between pillar (r = 2) and
    wall (W = 8); under the roof each support is the solid of revolution
    between the roof and the cone from the parting circle to its anchor."""
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder

    from latticegen2 import occ

    body = occ.solids(BRepAlgoAPI_Cut(BRepPrimAPI_MakeCylinder(8.0, 6.0).Shape(),
                                      BRepPrimAPI_MakeCylinder(2.0, 6.0).Shape()).Shape())[0]
    _lp, solids, stats = _supports(body)
    r, R, W = 2.0, 5.0, 8.0
    k = 2.0 * np.pi / np.sqrt(2.0)
    pillar = k * (R ** 3 / 6 - R * r ** 2 / 2 + r ** 3 / 3)
    wall = k * (W ** 3 / 3 - R * W ** 2 / 2 + R ** 3 / 6)
    # The wall's annulus is split around the pillar's cell, and again wherever
    # a piece's loft would close through the body: the count is not fixed,
    # the volumes are. Pocket fills come last; the one here is the ring the
    # wall cells close off in the roof's corner, inside the exempt band.
    assert stats["support_cells"] + stats["support_pockets"] == len(solids)
    assert stats["support_cells"] >= 3
    solids = solids[:stats["support_cells"]]
    # Measured +5.4 % on the crown and +1.8 % on the wall: the rays land on the
    # obstacle mesh, a polygon inscribed in each cylinder, so they reach a
    # little further than to the true surface, and they are smoothed. Both err
    # toward more support, never less.
    centres = [np.hypot(*occ.centroid(s)[:2]) for s in solids]
    crown = sum(v for v, c in zip([occ.volume(s) for s in solids], centres) if c < 4.0)
    walls = sum(v for v, c in zip([occ.volume(s) for s in solids], centres) if c >= 4.0)
    assert pillar <= crown <= 1.06 * pillar
    assert wall <= walls <= 1.06 * wall


@pytest.mark.parametrize("name", ["Hole-block-basic", "Hole-block-chamfered"])
def test_hole_block_supports_pass_their_gates_and_overlap_nothing(name):
    from latticegen2 import occ
    from latticegen2.support import assemble

    body = occ.read_step(os.path.join(TESTDIR, f"{name}.STEP"))
    lp, solids, stats = _supports(body, (90.0, 0.0, 0.0))
    assert len(solids) == stats["support_cells"] >= 15
    assert all(occ.is_valid(s) for s in solids)
    # Separate bodies, sharing faces, never volume.
    assemble.assemble([], solids, lp, NO_NODES)
    # Bounded faces: a handful per cell, never a mesh.
    assert max(len(occ.faces(s)) for s in solids) < 60
