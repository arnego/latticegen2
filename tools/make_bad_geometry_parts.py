"""Write the two parts the bad-geometry viewer's e2e scenarios run
(docs/specification.md §3.1, §6.1), for ``-cc 10 -t 1.5``:

* ``test/bad-geometry-warning.step`` -- a block one of whose floor corners the
  file places ``CORNER_SHIFT`` mm from where its three edges end. No edge runs
  on through the true corner, so the input repair (docs/algorithm.md §7.5)
  declines it, and the run *succeeds* with a warning naming both points: a
  yellow segment in the viewer.

* ``test/bad-geometry-failure.step`` -- ``input-defect.step``'s block, its roof
  split at ``x = SPLIT_X``, the split edge naming the front edge's corner vertex
  as there -- but with the front edge left in its two pieces. No single edge at
  the corner runs on through the split edge's true start, so the repair declines
  this one too, with a warning; and the near roof face's boundary is open, so
  tessellation covers part of it and the run *fails* at exit 3 naming that face:
  a yellow segment and a red face in one run.

Both are measured, not assumed: ``tools/e2e.py`` checks that each says what is
expected of it, so neither can pass by no longer testing the viewer.

Run from the repository root; the files it writes are committed.
"""

from __future__ import annotations

import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.dirname(__file__))

from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox                     # noqa: E402
from OCP.gp import gp_Pnt                                           # noqa: E402

import make_input_defect_part as defect                             # noqa: E402
from latticegen2 import occ                                         # noqa: E402

HERE = os.path.dirname(__file__)
WARNING_OUT = os.path.join(HERE, "..", "test", "bad-geometry-warning.step")
FAILURE_OUT = os.path.join(HERE, "..", "test", "bad-geometry-failure.step")

BLOCK = (40.0, 30.0, 20.0)
"""mm. The warning part's block, from the origin."""
CORNER_SHIFT = 0.5
"""mm. How far the file moves the block's corner vertex at the origin, along
x: past a tenth of ``t = 1.5`` mm, so the run warns rather than only logs."""


def warning_text(text: str) -> str:
    """Move the vertex at the origin, and only the vertex: its three edges'
    lines still end at the true corner."""
    ents = defect._entities(text)
    points = {}
    for i, body in ents.items():
        m = re.match(r"CARTESIAN_POINT\('[^']*',\(([^)]*)\)\)", body)
        if m:
            points[i] = tuple(round(float(v), 6) for v in m.group(1).split(","))
    vertex_point = next(defect._refs(body)[0] for body in ents.values()
                        if body.startswith("VERTEX_POINT")
                        and points[defect._refs(body)[0]] == (0.0, 0.0, 0.0))
    ents[vertex_point] = f"CARTESIAN_POINT('',({CORNER_SHIFT:g},0.,0.))"
    return _rebuild(text, ents)


def failure_text(text: str) -> str:
    """``input-defect``'s contradiction, the front edge left in two pieces."""
    ents = defect._entities(text)
    points = {}
    for i, body in ents.items():
        m = re.match(r"CARTESIAN_POINT\('[^']*',\(([^)]*)\)\)", body)
        if m:
            points[i] = tuple(round(float(v), 6) for v in m.group(1).split(","))
    vertex = {points[defect._refs(b)[0]]: i for i, b in ents.items()
              if b.startswith("VERTEX_POINT")}
    h, sx, depth = defect.H, defect.SPLIT_X, defect.DEPTH
    corner, q, back = vertex[(0.0, 0.0, h)], vertex[(sx, 0.0, h)], vertex[(sx, depth, h)]
    edges = {i: defect._refs(b)[:2] for i, b in ents.items() if b.startswith("EDGE_CURVE")}
    cross = next(i for i, (u, v) in edges.items() if {u, v} == {q, back})
    ents[cross] = re.sub(rf"#{q}\b", f"#{corner}", ents[cross], count=1)
    return _rebuild(text, ents)


def _rebuild(text: str, ents: dict) -> str:
    head, _, tail = text.partition("DATA;")
    _, _, end = tail.partition("ENDSEC;")
    data = "\n".join(f"#{i}={b};" for i, b in sorted(ents.items()))
    return f"{head}DATA;\n{data}\nENDSEC;{end}"


def _write(path: str, text: str) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)
    print(f"wrote {os.path.normpath(path)}")


def _clean(build, path: str) -> str:
    tmp = path + ".clean.step"
    build(tmp)
    with open(tmp, encoding="utf-8") as fh:
        text = fh.read()
    os.remove(tmp)
    return text


def main() -> None:
    box = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), *BLOCK).Shape()
    _write(WARNING_OUT, warning_text(_clean(
        lambda p: occ.write_step(box, p, "bad-geometry-warning"), WARNING_OUT)))
    _write(FAILURE_OUT, failure_text(_clean(defect.write_intended, FAILURE_OUT)))


if __name__ == "__main__":
    main()
