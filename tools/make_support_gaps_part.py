"""Write ``test/support-gaps.step``: the two coverage failures of the supported
rehearsal (docs/specification.md §6.1), in a part that runs in under a minute.

For ``-cc 10 -t 1.5`` at the default orientation. Its profile, in the XZ plane,
extruded along Y:

* **The keel**, left. The roof comes down in two planes at 45 deg and meets at a
  concave edge along Y: a ridge of shell hanging into the volume. Neither plane
  overhangs more than 60 deg, so no face of it is ceiling, and only the edge
  needs support (docs/algorithm.md §14.1).
* **The rib**, right. A flat roof at ``z = ROOF`` carries a short rib hanging
  from it, placed to cut the strut leaving lattice node (0, 0, 0) along ``e_0``
  clean through. What is left of that strut beyond the rib is a tip under the
  roof, smaller than ``t^3``, so the floating-body rule drops it while the node
  stays kept through its other arms. The ceiling above the tip is backed by
  nothing.

Run from the repository root; the file it writes is committed.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut                        # noqa: E402
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakePolygon  # noqa: E402
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakePrism  # noqa: E402
from OCP.gp import gp_Pnt, gp_Vec                                  # noqa: E402

from latticegen2 import occ                                        # noqa: E402

ROOF = 1.4
"""mm. Height of the flat roof: just above node (0, 0, 0)."""
FLOOR = -14.0
Y_HALF = 8.0
KEEL = (-10.0, -8.0)
"""``(x, z)`` of the keel. Its flanks rise at 45 deg to the roof's height."""
RIB_X = (1.2, 2.0)
RIB_BOTTOM = 0.0
RIB_Y_HALF = 3.0
OUT = os.path.join(os.path.dirname(__file__), "..", "test", "support-gaps.step")


def part():
    rise = ROOF - KEEL[1]
    profile = [(-20.0, FLOOR), (12.0, FLOOR), (12.0, ROOF), (KEEL[0] + rise, ROOF),
               KEEL, (KEEL[0] - rise, ROOF), (-20.0, ROOF)]
    wire = BRepBuilderAPI_MakePolygon()
    for x, z in profile:
        wire.Add(gp_Pnt(x, -Y_HALF, z))
    wire.Close()
    face = BRepBuilderAPI_MakeFace(wire.Wire(), True).Face()
    block = BRepPrimAPI_MakePrism(face, gp_Vec(0.0, 2.0 * Y_HALF, 0.0)).Shape()
    rib = BRepPrimAPI_MakeBox(
        gp_Pnt(RIB_X[0], -RIB_Y_HALF, RIB_BOTTOM),
        RIB_X[1] - RIB_X[0], 2.0 * RIB_Y_HALF, ROOF + 1.0 - RIB_BOTTOM,
    ).Shape()
    return occ.solids(BRepAlgoAPI_Cut(block, rib).Shape())[0]


if __name__ == "__main__":
    solid = part()
    assert occ.is_valid(solid), "the part is not a valid solid"
    occ.write_step(solid, OUT, "support-gaps")
    print(f"{os.path.normpath(OUT)}: {occ.volume(solid):.1f} mm^3, "
          f"{occ.count_subshapes(solid)[0]} faces")
