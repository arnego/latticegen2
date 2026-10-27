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
    # a piece's loft would close through the body, but its pieces share the
    # faces along their cuts and are one body again (§14.3 S6). The crown meets
    # it only along the parting circle, and stays apart.
    assert stats["support_cells"] >= 3
    assert stats["support_pockets"] == 0
    assert len(solids) == stats["support_bodies"] == 2
    # Measured +5.4 % on the crown and +1.8 % on the wall: the rays land on the
    # obstacle mesh, a polygon inscribed in each cylinder, so they reach a
    # little further than to the true surface, and they are smoothed. Both err
    # toward more support, never less.
    # The wall's ring is centred on the axis like the crown: tell them by reach.
    reach = [float(np.abs(np.vstack(occ.bounding_box(s))[:, :2]).max()) for s in solids]
    crown = sum(occ.volume(s) for s, r_ in zip(solids, reach) if r_ < 6.0)
    walls = sum(occ.volume(s) for s, r_ in zip(solids, reach) if r_ >= 6.0)
    assert pillar <= crown <= 1.06 * pillar
    assert wall <= walls <= 1.06 * wall


@pytest.mark.parametrize("name", ["Hole-block-basic", "Hole-block-chamfered"])
def test_hole_block_supports_pass_their_gates_and_overlap_nothing(name):
    from latticegen2 import occ
    from latticegen2.support import assemble

    body = occ.read_step(os.path.join(TESTDIR, f"{name}.STEP"))
    lp, solids, stats = _supports(body, (90.0, 0.0, 0.0))
    assert len(solids) == stats["support_bodies"]
    assert all(occ.is_valid(s) for s in solids)
    # Cells meet across partings along an edge, and stay separate bodies; on
    # the chamfered block cells of the two regions either side of the bare
    # facet overlap where their rays cross, and are fused (§14.3 S6).
    want = {"Hole-block-basic": (20, 0), "Hole-block-chamfered": (13, 6)}[name]
    assert (len(solids), stats["support_overlapping_pairs"]) == want
    # Never sharing volume.
    assemble.assemble([], solids, lp, NO_NODES)
    # Bounded faces: a handful per cell, never a mesh.
    assert max(len(occ.faces(s)) for s in solids) < 60
    # The supports run into every crease where the roof meets a wall or a
    # tube: no bare band along it (§14.3 S3d). It was a 0.408 mm bevel, 0.236 mm
    # from the crease, while ceiling within 0.5 mm reach was left out.
    gaps = _crease_gaps(body, solids, lp.b)
    assert len(gaps) >= 20
    assert max(gaps) < 1e-6


def _crease_gaps(body, solids, b):
    """Distance from points along every roof/wall crease to the nearest support."""
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
    from OCP.BRepExtrema import BRepExtrema_DistShapeShape
    from OCP.BRepTools import BRepTools
    from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE
    from OCP.TopExp import TopExp
    from OCP.TopTools import TopTools_IndexedDataMapOfShapeListOfShape
    from OCP.TopoDS import TopoDS

    from latticegen2 import occ
    from latticegen2.support.overhang import _SurfaceNormals

    supports = occ.compound(list(solids))
    b = np.asarray(b, dtype=float)

    def normal(face):
        u0, u1, v0, v1 = BRepTools.UVBounds_s(face)
        return _SurfaceNormals(face).at(0.5 * (u0 + u1), 0.5 * (v0 + v1))[1]

    owners = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(occ.solids(body)[0], TopAbs_EDGE, TopAbs_FACE, owners)
    out = []
    for i in range(1, owners.Extent() + 1):
        faces = [TopoDS.Face_s(f) for f in owners.FindFromIndex(i)]
        if len(faces) != 2:
            continue
        nb = [float(normal(f) @ b) for f in faces]
        if not (max(nb) > np.sin(np.radians(60.0)) and abs(min(nb, key=abs)) < 0.2):
            continue
        curve = BRepAdaptor_Curve(TopoDS.Edge_s(owners.FindKey(i)))
        for t in np.linspace(curve.FirstParameter(), curve.LastParameter(), 5)[1:-1]:
            v = BRepBuilderAPI_MakeVertex(curve.Value(t)).Vertex()
            out.append(BRepExtrema_DistShapeShape(v, supports).Value())
    return out
