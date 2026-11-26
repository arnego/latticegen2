"""Write ``test/input-defect.step``: a block whose roof carries the defect of a
user's file (docs/algorithm.md §7.5), for ``-cc 10 -t 1.5 --support``.

The user's ``TDHX_error_report_test.stp`` has an edge whose curve starts
1.334 mm from the vertex the file names for it, that vertex lying on a straight
edge that runs on past the curve's real start. The reader then gives the face
beside it a degenerated edge sweeping the gap, two faces claim the same strip of
the straight edge, and support near it could not be trimmed.

This builds the same contradiction on a block, where the roof -- a ceiling of
the volume, so supported -- is split by an edge across it at ``x = SPLIT_X``:

1. the block is written to STEP cleanly, its roof split into two faces;
2. the split edge, which starts at ``(SPLIT_X, 0, H)`` on the roof's front
   edge, is made to name the front edge's corner vertex ``(0, 0, H)`` instead;
3. the two pieces of the front edge are made one edge from that corner to the
   far end, used by the front face and the far roof face, and the near roof
   face loses its piece.

Every loop stays closed in the file's topology, as in the user's; only the
geometry disagrees, by ``SPLIT_X``.

Run from the repository root; the file it writes is committed.
"""

from __future__ import annotations

import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge             # noqa: E402
from OCP.BRepFeat import BRepFeat_SplitShape                        # noqa: E402
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox                     # noqa: E402
from OCP.TopoDS import TopoDS                                       # noqa: E402
from OCP.gp import gp_Pnt                                           # noqa: E402

from latticegen2 import occ                                         # noqa: E402

LENGTH, DEPTH, H = 40.0, 30.0, 20.0
"""mm. The block, from the origin."""
SPLIT_X = 13.0
"""mm. Where the split edge really starts on the front edge: the gap."""
OUT = os.path.join(os.path.dirname(__file__), "..", "test", "input-defect.step")


def write_intended(path: str) -> None:
    """The block as the defective file means it, its roof split, written cleanly."""
    box = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), LENGTH, DEPTH, H).Shape()
    roof = next(f for f in occ.faces(box)
                if abs(occ.centroid(f)[2] - H) < 1e-9)
    split = BRepFeat_SplitShape(box)
    edge = BRepBuilderAPI_MakeEdge(gp_Pnt(SPLIT_X, 0, H), gp_Pnt(SPLIT_X, DEPTH, H)).Edge()
    split.Add(edge, TopoDS.Face_s(roof))
    split.Build()
    occ.write_step(split.Shape(), path, "input-defect")


def _entities(text: str) -> dict[int, str]:
    return {int(m.group(1)): m.group(2)
            for m in re.finditer(r"^#(\d+)\s*=\s*(.*?);\s*$", text, re.M | re.S)}


def _refs(body: str) -> list[int]:
    return [int(r) for r in re.findall(r"#(\d+)", body)]


def inject(text: str) -> str:
    ents = _entities(text)
    points = {}
    for i, body in ents.items():
        m = re.match(r"CARTESIAN_POINT\('[^']*',\(([^)]*)\)\)", body)
        if m:
            points[i] = tuple(round(float(v), 6) for v in m.group(1).split(","))
    vertex = {}
    for i, body in ents.items():
        if body.startswith("VERTEX_POINT"):
            vertex[points[_refs(body)[0]]] = i
    corner, q, far = (vertex[(0.0, 0.0, H)], vertex[(SPLIT_X, 0.0, H)],
                      vertex[(LENGTH, 0.0, H)])
    back = vertex[(SPLIT_X, DEPTH, H)]
    edges = {i: _refs(b)[:2] for i, b in ents.items() if b.startswith("EDGE_CURVE")}

    def edge(a, b):
        return next(i for i, (u, v) in edges.items() if {u, v} == {a, b})

    near, rest, cross = edge(corner, q), edge(q, far), edge(q, back)

    def swap_vertex(i, old, new):
        ents[i] = re.sub(rf"#{old}\b", f"#{new}", ents[i], count=1)

    swap_vertex(cross, q, corner)      # the split edge names the wrong vertex
    swap_vertex(rest, q, corner)       # the front edge is one edge again
    oriented = [i for i, b in ents.items()
                if b.startswith("ORIENTED_EDGE") and _refs(b)[-1] == near]
    for i, b in list(ents.items()):
        if b.startswith("EDGE_LOOP"):
            kept = [r for r in _refs(b) if r not in oriented]
            if len(kept) != len(_refs(b)):
                ents[i] = "EDGE_LOOP('',(" + ",".join(f"#{r}" for r in kept) + "))"
    for i in oriented + [near, q]:
        del ents[i]
    head, _, tail = text.partition("DATA;")
    _, _, end = tail.partition("ENDSEC;")
    data = "\n".join(f"#{i}={b};" for i, b in sorted(ents.items()))
    return f"{head}DATA;\n{data}\nENDSEC;{end}"


def main() -> None:
    tmp = OUT + ".clean.step"
    write_intended(tmp)
    with open(tmp, encoding="utf-8") as fh:
        text = fh.read()
    os.remove(tmp)
    with open(OUT, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(inject(text))
    print(f"wrote {os.path.normpath(OUT)}")


if __name__ == "__main__":
    main()
