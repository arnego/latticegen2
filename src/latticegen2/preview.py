"""A coarse triangle mesh of the input, for the window's orientation view.

The window imports no geometry kernel (:mod:`latticegen2.gui`), so it cannot
read a STEP file itself. It runs this instead, as a short-lived child:

    python src/main.py --preview-mesh <input.step> <out.npz>

which reads the file, tessellates it coarsely, and writes plain arrays the
window can draw. Like ``--progress-stream`` it is transport rather than a
parameter, so it is absent from the usage text (specification.md §3.1).

The mesh is for looking at, not for measuring: it is capped at a few thousand
triangles so a canvas can redraw it, and nothing in a run reads it.
"""

from __future__ import annotations

import numpy as np

MAX_TRIANGLES = 4000
"""What a ``tkinter.Canvas`` redraws comfortably."""


def preview_mesh(step_path: str):
    """``(verts, tris)`` with every triangle wound outward, capped in size."""
    from OCP.BRepTools import BRepTools
    from OCP.TopAbs import TopAbs_Orientation

    from . import occ

    occ.quiet_kernel()
    shape = occ.read_step(step_path)
    lo, hi = occ.bounding_box(shape)
    deflection = float(np.linalg.norm(hi - lo)) / 150.0
    verts = tris = None
    for _ in range(6):
        BRepTools.Clean_s(shape)
        occ.mesh_shape(shape, deflection, 0.5)
        vs, ts, offset = [], [], 0
        for face in occ.faces(shape):
            got = occ.face_triangulation(face)
            if got is None:
                continue
            v, f = got
            if face.Orientation() == TopAbs_Orientation.TopAbs_REVERSED:
                f = f[:, ::-1]
            vs.append(v)
            ts.append(f + offset)
            offset += len(v)
        if not ts:
            raise ValueError("the input produced no triangles")
        verts, tris = np.vstack(vs), np.vstack(ts)
        if len(tris) <= MAX_TRIANGLES:
            break
        deflection *= 2.5
    if len(tris) > MAX_TRIANGLES:
        # Curvature the mesher will not coarsen past. Thin evenly: the view
        # shows gaps rather than refusing to show the part at all.
        tris = tris[np.linspace(0, len(tris) - 1, MAX_TRIANGLES).astype(np.int64)]
    return verts, tris


def main(argv: list[str]) -> int:
    """``--preview-mesh <input.step> <out.npz>``. Exit 0 and a file, or 3."""
    import sys

    if len(argv) != 2:
        print("FAILED: --preview-mesh needs <input.step> <out.npz>", file=sys.stderr)
        return 2
    try:
        verts, tris = preview_mesh(argv[0])
        np.savez(argv[1], verts=verts, tris=tris)
    except Exception as exc:                                   # noqa: BLE001
        print(f"FAILED: could not build a preview of {argv[0]}: {exc}",
              file=sys.stderr)
        return 3
    return 0
