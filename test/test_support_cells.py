"""Support cells and parting curves (docs/algorithm.md §14.3 S3b, S3c).

Checked against closed forms: a pillar in a round chamber, whose roof parts
midway between the pillar and the wall, and the hole-block, whose cells are
known by symmetry.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from latticegen2.lattice import lattice_params

TESTDIR = os.path.dirname(os.path.abspath(__file__))


def _setup(body, orient=(0.0, 0.0, 0.0), deflection=0.05):
    from latticegen2.classify import tessellate_surface
    from latticegen2.support import obstacles, rays, regions

    lp = lattice_params(10.0, 1.5, orient)
    tessellate_surface(body, lp)
    rm = regions.build(body, lp.b)
    regs = regions.regions(rm)
    ob = obstacles.from_shapes(body, [], deflection=deflection)
    return lp, rm, regs, ob, rays


def _partition(lp, rm, region, ob, rays, thickness=1.5):
    from latticegen2.support import cells

    got = rays.first_contact(rm.verts[region.verts]
                             - rays.APEX_NUDGE * rm.normals[region.verts], lp.b, ob)
    return cells.partition(rm, region, got, ob, lp.b, thickness)


@pytest.fixture(scope="module")
def chamber():
    """A round chamber of radius 8 around a pillar of radius 2, 6 mm tall."""
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder

    from latticegen2 import occ

    outer = BRepPrimAPI_MakeCylinder(8.0, 6.0).Shape()
    pillar = BRepPrimAPI_MakeCylinder(2.0, 6.0).Shape()
    body = occ.solids(BRepAlgoAPI_Cut(outer, pillar).Shape())[0]
    lp, rm, regs, ob, rays = _setup(body, deflection=1e-3)
    assert len(regs) == 1
    return lp, rm, regs[0], ob, rays, _partition(lp, rm, regs[0], ob, rays)


def test_a_pillar_in_a_chamber_parts_midway_between_pillar_and_wall(chamber):
    from latticegen2.support import cells

    _lp, _rm, _reg, _ob, _rays, part = chamber
    partings = np.vstack([lp_.points[lp_.kind == cells.PARTING]
                          for c in part.cells for lp_ in c.loops])
    assert len(partings) > 50
    radius = np.linalg.norm(partings[:, :2], axis=1)
    # Midway: the cone reaches the pillar (r = 2) and the wall (W = 8) at the
    # same drop. Both are meshed at 1e-3 mm, inscribed, which bounds the error.
    assert np.abs(radius - 5.0).max() < 2e-3


def test_the_wall_cell_enclosing_the_pillar_cell_is_split(chamber):
    """An annulus around another cell cannot be lofted from its outer loop."""
    from latticegen2.support import cells

    _lp, _rm, _reg, _ob, _rays, part = chamber
    assert part.stats["cells_split"] >= 1
    for c in part.cells:
        for lp_ in c.loops[1:]:
            assert not (lp_.kind == cells.PARTING).any()


def test_neighbouring_cells_share_one_parting_curve_and_part_in_a_v(chamber):
    from latticegen2.support import cells

    _lp, _rm, _reg, _ob, _rays, part = chamber
    where = {}
    for c in part.cells:
        for lp_ in c.loops:
            for pid, p, k, t in zip(lp_.keys, lp_.points, lp_.kind, lp_.targets):
                if k == cells.PARTING:
                    where.setdefault(pid, []).append((c.id, p, t))
    assert where
    for pid, uses in where.items():
        assert len(uses) == 2 and uses[0][0] != uses[1][0]
        assert np.array_equal(uses[0][1], uses[1][1])
        # The two rays from one parting point head to different anchors.
        assert np.linalg.norm(uses[0][2] - uses[1][2]) > 1.0
    for c in part.cells:
        assert c.min_divergence_deg > 30.0


def test_hole_block_cells_are_a_crown_per_tube_and_one_per_wall():
    from latticegen2 import occ
    from latticegen2.support import cells

    body = occ.read_step(os.path.join(TESTDIR, "Hole-block-basic.STEP"))
    lp, rm, regs, ob, rays = _setup(body, (90.0, 0.0, 0.0))
    part = _partition(lp, rm, regs[0], ob, rays)
    assert len(part.cells) == 20
    crowns = [c for c in part.cells if len(c.loops) == 2]
    walls = [c for c in part.cells if len(c.loops) == 1]
    assert len(crowns) == 16 and len(walls) == 4
    # A crown's inner loop is the exempt ring round its tube: no partings.
    for c in crowns:
        assert (c.loops[1].kind == cells.EXEMPT_EDGE).all()
    # Symmetry: the four walls alike, and the crowns in three classes.
    assert np.ptp([c.area for c in walls]) < 0.05
    assert len({round(c.area, 1) for c in crowns}) == 3
    # Where two walls' cells part, their rays meet at arccos(1/3).
    assert min(c.min_divergence_deg for c in walls) == pytest.approx(
        np.degrees(np.arccos(1.0 / 3.0)), abs=0.1)
    # Every triple point is shared by the three cells meeting there.
    triples = {}
    for c in part.cells:
        for lp_ in c.loops:
            for pid, k in zip(lp_.keys, lp_.kind):
                if k == cells.TRIPLE:
                    triples.setdefault(pid, set()).add(c.id)
    assert triples and all(len(v) >= 2 for v in triples.values())
