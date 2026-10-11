"""The print-orientation view (specification.md §3.1).

A small canvas showing the print bed, the part's coordinate axes and — once the
preview mesh has arrived — the input geometry standing on the bed the way the
three angles say it will. Dragging turns the part; the boxes beside it follow,
and editing the boxes turns the part. With support enabled, every triangle that
overhangs at the current orientation and limit is tinted, and the tint is
recomputed on each change so the user sees what will be supported *before*
starting an hour-long run.

A ``tkinter.Canvas`` and a painter's sort, on purpose: no dependency, and the
same on Windows and Linux. All the arithmetic is in
:mod:`latticegen2.gui.view3d`.
"""

from __future__ import annotations

import tkinter as tk

import numpy as np

from .. import orient as _orient
from . import view3d

WIDTH, HEIGHT = 210, 190
DRAG_TRIANGLES = 700
"""Triangles drawn while the mouse is down. A canvas redraws a few thousand
polygons in tens of milliseconds, which is fine for a click and not for a drag;
the full mesh comes back on release."""

AXIS_COLOURS = ("#d33", "#2a2", "#36c")
OVERHANG_COLOUR = (0.86, 0.33, 0.24)
BODY_COLOUR = (0.62, 0.68, 0.76)


def _hex(rgb, k: float = 1.0) -> str:
    return "#%02x%02x%02x" % tuple(int(max(0.0, min(1.0, c * k)) * 255) for c in rgb)


class OrientationView(tk.Canvas):
    def __init__(self, master, on_drag, **kw):
        super().__init__(master, width=WIDTH, height=HEIGHT, highlightthickness=1,
                         highlightbackground="#bdbdbd", background="#fbfbfb", **kw)
        self._on_drag = on_drag
        self._verts = None
        self._tris = None
        self._normals = None
        self._orient = _orient.IDENTITY
        self._overhang: float | None = None
        self._enabled = True
        self._last = None
        self._R_drag = None
        self.bind("<ButtonPress-1>", self._press)
        self.bind("<B1-Motion>", self._motion)
        self.bind("<ButtonRelease-1>", self._release)
        self.redraw()

    # -- state -------------------------------------------------------------

    def set_mesh(self, verts, tris) -> None:
        """Give the view the input geometry, or ``None`` to show axes only."""
        if verts is None or tris is None or len(tris) == 0:
            self._verts = self._tris = self._normals = None
        else:
            self._verts = np.asarray(verts, dtype=float)
            self._tris = np.asarray(tris, dtype=np.int64)
            self._normals = view3d.triangle_normals(self._verts, self._tris)
        self.redraw()

    def set_state(self, orient, overhang: float | None) -> None:
        """The orientation to show, and the overhang limit or ``None`` for no tint."""
        orient = tuple(float(a) for a in orient)
        if orient == self._orient and overhang == self._overhang:
            return
        self._orient = orient
        self._overhang = overhang
        self.redraw()

    def set_enabled(self, enabled: bool) -> None:
        self._enabled = bool(enabled)

    # -- dragging ----------------------------------------------------------

    def _press(self, event) -> None:
        if not self._enabled:
            return
        self._last = (event.x, event.y)
        self._R_drag = _orient.rotation_matrix(self._orient)

    def _motion(self, event) -> None:
        if self._last is None:
            return
        dx, dy = event.x - self._last[0], event.y - self._last[1]
        self._last = (event.x, event.y)
        # Accumulated as a matrix and only *reported* as angles: converting to
        # the three boxes and back on every pixel would walk the part toward
        # the nearest tenth of a degree and stick at the gimbal singularity.
        self._R_drag = view3d.drag_rotation(self._R_drag, dx, dy)
        self._orient = view3d.snapped_angles(self._R_drag)
        self.redraw(dragging=True)
        self._on_drag(self._orient)

    def _release(self, _event) -> None:
        if self._last is None:
            return
        self._last = None
        self._R_drag = None
        self.redraw()

    # -- drawing -----------------------------------------------------------

    def redraw(self, dragging: bool = False) -> None:
        self.delete("all")
        R = (self._R_drag if self._R_drag is not None
             else _orient.rotation_matrix(self._orient))
        if self._verts is not None:
            verts, _shift = view3d.place_on_bed(self._verts, R)
            extent = float(np.abs(verts).max())
        else:
            verts, extent = None, 1.0
        half = max(extent, 1e-6) * 1.25
        scale = min(WIDTH, HEIGHT) * 0.5 / (half * 1.45)
        cx, cy = WIDTH / 2.0, HEIGHT * 0.56

        def to_screen(points):
            xy, depth = view3d.project(points)
            return np.column_stack([cx + xy[:, 0] * scale, cy - xy[:, 1] * scale]), depth

        # The bed: a square plate with a grid, drawn first — everything stands on it.
        corners = np.array([[-half, -half, 0], [half, -half, 0],
                            [half, half, 0], [-half, half, 0]], dtype=float)
        s, _ = to_screen(corners)
        self.create_polygon(*s.ravel(), fill="#e9edf2", outline="#9aa6b2")
        for i in range(1, 6):
            f = -half + 2.0 * half * i / 6.0
            for a, b in (([f, -half, 0], [f, half, 0]), ([-half, f, 0], [half, f, 0])):
                s, _ = to_screen(np.array([a, b], dtype=float))
                self.create_line(*s.ravel(), fill="#cfd6de")

        if verts is not None:
            keep = (view3d.thin(self._tris, DRAG_TRIANGLES) if dragging
                    else np.arange(len(self._tris)))
            tris = self._tris[keep]
            normals_bed = self._normals[keep] @ R.T
            light = view3d.shade(normals_bed)
            if self._overhang is not None:
                tint = view3d.overhang_mask(self._normals[keep],
                                            _orient.euler_from_matrix(R)
                                            if self._R_drag is not None else self._orient,
                                            self._overhang)
            else:
                tint = np.zeros(len(tris), dtype=bool)
            screen, depth = to_screen(verts)
            order = np.argsort(depth[tris].mean(axis=1))
            # Back faces are hidden by the faces in front of them; drawing only
            # the ones turned toward the eye halves the polygons for a closed body.
            front = normals_bed @ view3d.EYE > 0.0
            for i in order:
                if not front[i]:
                    continue
                pts = screen[tris[i]].ravel()
                colour = _hex(OVERHANG_COLOUR if tint[i] else BODY_COLOUR, light[i])
                self.create_polygon(*pts, fill=colour, outline=colour)

        # The part's own axes, always on top. Drawn at a corner of the bed
        # rather than at the part's origin: what they show is which way the
        # part is turned, and a CAD origin can lie far outside the part.
        origin = np.array([-half * 0.8, -half * 0.8, 0.0])
        length = half * 0.45
        for k, name in enumerate("XYZ"):
            tip = origin + R[:, k] * length
            s, _ = to_screen(np.array([origin, tip]))
            self.create_line(*s.ravel(), fill=AXIS_COLOURS[k], width=2, arrow="last")
            self.create_text(s[1, 0] + 6, s[1, 1] - 6, text=name,
                             fill=AXIS_COLOURS[k], font=("TkDefaultFont", 8, "bold"))
        self.create_text(6, HEIGHT - 6, anchor="sw", fill="#667",
                         font=("TkDefaultFont", 7),
                         text="print bed" if verts is not None
                         else "print bed — choose an input to see it here")
