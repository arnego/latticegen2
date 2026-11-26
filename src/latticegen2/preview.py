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
"""What a ``tkinter.Canvas`` redraws comfortably, in the window's small view."""

VIEWER_TRIANGLES = 12000
"""The bad-geometry viewer's cap (specification.md §3.1). Finer, because a
mark on a face is only as findable as that face is drawn, and that view is
redrawn only when the camera moves -- with a thinned copy while it does."""


def preview_mesh(step_path: str, max_triangles: int = MAX_TRIANGLES):
    """``(verts, tris, face, face_lo, face_hi)``: every triangle wound outward
    and capped in number, the input face each came from -- by its order in
    :func:`latticegen2.occ.faces` -- and each face's bounding box.

    The boxes are what the bad-geometry viewer matches a failure's face
    against, rather than the index: the pipeline repairs the input before it
    names a face (docs/algorithm.md §7.5), and a face's place in the traversal
    order is not something either side promises.
    """
    from OCP.BRepTools import BRepTools
    from OCP.TopAbs import TopAbs_Orientation

    from . import occ

    occ.quiet_kernel()
    shape = occ.read_step(step_path)
    lo, hi = occ.bounding_box(shape)
    # A curved face's triangle count goes as one over the deflection, so a
    # larger cap starts proportionally finer; the small view's is unchanged.
    deflection = float(np.linalg.norm(hi - lo)) / 150.0 * (MAX_TRIANGLES / max_triangles)
    faces = list(occ.faces(shape))
    boxes = [occ.bounding_box(f) for f in faces]
    face_lo = np.array([b[0] for b in boxes], dtype=float).reshape(-1, 3)
    face_hi = np.array([b[1] for b in boxes], dtype=float).reshape(-1, 3)
    verts = tris = owner = None
    for _ in range(6):
        BRepTools.Clean_s(shape)
        occ.mesh_shape(shape, deflection, 0.5)
        vs, ts, fs, offset = [], [], [], 0
        for i, face in enumerate(faces):
            got = occ.face_triangulation(face)
            if got is None:
                continue
            v, f = got
            if face.Orientation() == TopAbs_Orientation.TopAbs_REVERSED:
                f = f[:, ::-1]
            vs.append(v)
            ts.append(f + offset)
            fs.append(np.full(len(f), i, dtype=np.int64))
            offset += len(v)
        if not ts:
            raise ValueError("the input produced no triangles")
        verts, tris, owner = np.vstack(vs), np.vstack(ts), np.concatenate(fs)
        if len(tris) <= max_triangles:
            break
        deflection *= 2.5
    if len(tris) > max_triangles:
        # Curvature the mesher will not coarsen past. Thin evenly: the view
        # shows gaps rather than refusing to show the part at all.
        keep = thin_keeping_faces(owner, max_triangles)
        tris, owner = tris[keep], owner[keep]
    return verts, tris, owner, face_lo, face_hi


def thin_keeping_faces(owner: np.ndarray, limit: int) -> np.ndarray:
    """Indices of at most ``limit`` triangles, spread evenly, keeping at least
    one of every face -- so no face a mark names can vanish from the view.
    Where there are more faces than ``limit``, one triangle of each is kept."""
    owner = np.asarray(owner)
    n = len(owner)
    if n <= limit:
        return np.arange(n)
    _, first = np.unique(owner, return_index=True)
    rest = np.setdiff1d(np.arange(n), first, assume_unique=False)
    room = max(limit - len(first), 0)
    extra = rest[np.linspace(0, len(rest) - 1, room).astype(np.int64)] if room and len(rest) else []
    return np.unique(np.concatenate([first, np.asarray(extra, dtype=np.int64)]))


def main(argv: list[str]) -> int:
    """``--preview-mesh <input.step> <out.npz> [--max-triangles N]``.
    Exit 0 and a file, 2 for bad arguments, or 3."""
    import sys

    cap = MAX_TRIANGLES
    if len(argv) == 4 and argv[2] == "--max-triangles":
        try:
            cap = int(argv[3])
        except ValueError:
            cap = 0
        argv = argv[:2]
    if len(argv) != 2 or cap < 1:
        print("FAILED: --preview-mesh needs <input.step> <out.npz> "
              "[--max-triangles N]", file=sys.stderr)
        return 2
    try:
        verts, tris, face, face_lo, face_hi = preview_mesh(argv[0], cap)
        np.savez(argv[1], verts=verts, tris=tris, face=face,
                 face_lo=face_lo, face_hi=face_hi)
    except Exception as exc:                                   # noqa: BLE001
        print(f"FAILED: could not build a preview of {argv[0]}: {exc}",
              file=sys.stderr)
        return 3
    return 0
