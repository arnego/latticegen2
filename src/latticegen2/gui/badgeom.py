"""The bad-geometry viewer (specification.md §3.1).

A large, resizable window over the main one, showing the input at the run's
print orientation and marking on it every place the run's warnings (yellow) and
failure (red) point at. Click a mark for what was said about it.

It draws what the orientation view draws -- the same shading, through the same
renderer (:func:`latticegen2.gui.orientview.draw_mesh`), the same blue for what
will be supported -- without the print bed, which adds nothing when looking for
a defect. Where the small view turns the *part*, this one turns the *view*: a
left drag here is the small view's drag (:func:`view3d.drag_rotation`), applied
to a view rotation in front of the fixed camera, so the run's orientation is
never touched.

The marks are drawn after the mesh, on top of it, so a mark behind the part is
still seen -- and can be found through the issue list, which centres the view on
whatever is picked there.
"""

from __future__ import annotations

import tkinter as tk
from tkinter import ttk

import numpy as np

from .. import orient as _orient
from . import marks as _marks
from . import view3d
from .orientview import AXIS_COLOURS, BODY_COLOUR, OVERHANG_COLOUR, _hex, draw_mesh

BACKGROUND = "#fbfbfb"
#: The discreet line of help in the top-left corner.
CONTROLS_TEXT = ("Rotate: left drag  ·  Pan: middle drag  ·  Zoom: wheel  ·  "
                 "F: fit all  ·  Click a mark for details")
#: How far a left press may move and still be a click on a mark, in pixels.
#: Past it the press is a rotation, so a click never turns the part.
CLICK_SLOP = 3
#: Triangles drawn while a button is held; the full mesh comes back on release.
DRAG_TRIANGLES = 1500
#: The wheel's zoom per notch.
WHEEL_FACTOR = 1.15
MIN_SIZE = (480, 360)
POINT_RADIUS = 5
REDRAW_DEBOUNCE_MS = 30


class BadGeometryViewer(tk.Toplevel):
    """The viewer window. ``mesh`` is ``(verts, tris, face, face_lo, face_hi)``
    with ``face`` and the boxes possibly ``None``, or ``None`` for no mesh."""

    def __init__(self, master, orient=_orient.IDENTITY, mesh=None, marks=(),
                 on_close=None):
        super().__init__(master)
        self.title("latticegen2 — bad geometry")
        # Kept over the window that opened it -- but only a shown one: on
        # Windows a window transient to a withdrawn master is never mapped.
        if master.winfo_viewable():
            self.transient(master)
        sw, sh = self.winfo_screenwidth(), self.winfo_screenheight()
        w, h = max(MIN_SIZE[0], int(sw * 0.7)), max(MIN_SIZE[1], int(sh * 0.7))
        self.geometry(f"{w}x{h}+{max(0, (sw - w) // 2)}+{max(0, (sh - h) // 3)}")
        self.minsize(*MIN_SIZE)
        self._on_close = on_close

        self._R = _orient.rotation_matrix(tuple(float(a) for a in orient))
        self._orient = tuple(float(a) for a in orient)
        self._view = np.eye(3)
        self._scale = 1.0
        self._pan = (0.0, 0.0)
        self._fitted = False
        self._press = None
        self._mode = None
        self._last = None
        self._detail = None            # (mark index, place index) of the open box
        self._redraw_pending = None
        self._marks: list = []
        self._face_match: dict = {}    # (mark, place) -> preview face index
        self._verts = self._tris = self._normals = self._owner = None
        self._face_lo = self._face_hi = None
        self._bed = None
        self._shift = np.zeros(3)
        self._pivot = np.zeros(3)
        self._extent = 1.0

        self.canvas = tk.Canvas(self, background=BACKGROUND, highlightthickness=0)
        self.canvas.pack(fill="both", expand=True)

        self.overhangs_var = tk.BooleanVar(value=False)
        self._build_overlays()

        c = self.canvas
        c.bind("<Configure>", lambda _e: self._schedule_redraw())
        c.bind("<ButtonPress-1>", self._left_press)
        c.bind("<B1-Motion>", self._left_motion)
        c.bind("<ButtonRelease-1>", self._left_release)
        c.bind("<ButtonPress-2>", self._middle_press)
        c.bind("<B2-Motion>", self._middle_motion)
        c.bind("<ButtonRelease-2>", self._middle_release)
        c.bind("<MouseWheel>", self._wheel)
        c.bind("<Button-4>", lambda e: self._zoom(WHEEL_FACTOR, (e.x, e.y)))
        c.bind("<Button-5>", lambda e: self._zoom(1.0 / WHEEL_FACTOR, (e.x, e.y)))
        for key in ("<KeyPress-f>", "<KeyPress-F>", "<Home>"):
            self.bind(key, lambda _e: self.fit())
        self.bind("<Escape>", self._escape)
        self.protocol("WM_DELETE_WINDOW", self.close)

        self.set_mesh(mesh)
        self.set_marks(marks)
        c.focus_set()

    # -- overlays -----------------------------------------------------------

    def _build_overlays(self) -> None:
        """Two boxes over the viewport's top-right corner, one above the other:
        the overhang tick box, always shown, and the issue list, which folds to
        its title bar."""
        holder = tk.Frame(self.canvas, background=BACKGROUND)
        holder.place(relx=1.0, x=-10, y=10, anchor="ne")
        self.overlay = holder

        over = tk.Frame(holder, background="white", highlightthickness=1,
                        highlightbackground="#bdbdbd")
        over.pack(side="top", anchor="e", fill="x")
        self.overhangs_check = ttk.Checkbutton(
            over, text="Show overhangs", variable=self.overhangs_var,
            command=self._schedule_redraw)
        self.overhangs_check.pack(side="left", padx=6, pady=3)
        self.overhangs_box = over

        issues = tk.Frame(holder, background="white", highlightthickness=1,
                          highlightbackground="#bdbdbd")
        issues.pack(side="top", anchor="e", fill="x", pady=(6, 0))
        self.issues_box = issues
        bar = tk.Frame(issues, background="#eef0f3")
        bar.pack(side="top", fill="x")
        self.issues_title = tk.Label(bar, text="Issues", background="#eef0f3",
                                     font=("TkDefaultFont", 9, "bold"))
        self.issues_title.pack(side="left", padx=6, pady=2)
        self.collapse_button = ttk.Button(bar, text="▾", width=3,
                                          command=self.toggle_issues)
        self.collapse_button.pack(side="right", padx=2, pady=1)

        body = tk.Frame(issues, background="white")
        body.pack(side="top", fill="both")
        self.issues_body = body
        self.tree = ttk.Treeview(body, show="tree", height=8, selectmode="browse")
        self.tree.column("#0", width=300, stretch=True)
        self.tree.tag_configure("error", background="#ffd9d9")
        self.tree.tag_configure("warning", background="#fff4c2")
        scroll = ttk.Scrollbar(body, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=scroll.set)
        self.tree.pack(side="left", fill="both")
        scroll.pack(side="right", fill="y")
        self.tree.bind("<<TreeviewSelect>>", self._picked)
        self.issues_expanded = True

    def toggle_issues(self) -> None:
        """Fold the issue list to its title bar, or open it again."""
        if self.issues_expanded:
            self.issues_body.pack_forget()
            self.collapse_button.configure(text="▸")
        else:
            self.issues_body.pack(side="top", fill="both")
            self.collapse_button.configure(text="▾")
        self.issues_expanded = not self.issues_expanded

    def _fill_tree(self) -> None:
        self.tree.delete(*self.tree.get_children())
        order = sorted(range(len(self._marks)), key=lambda i: not self._marks[i].is_error)
        for i in order:
            mark = self._marks[i]
            tag = "error" if mark.is_error else "warning"
            parent = self.tree.insert("", "end", iid=f"m{i}", text=mark.headline(),
                                      open=len(mark.places) <= 3, tags=(tag,))
            for j, place in enumerate(mark.places):
                self.tree.insert(parent, "end", iid=f"m{i}_{j}",
                                 text=place.text or place.kind, tags=(tag,))
        n_err = sum(m.is_error for m in self._marks)
        n_warn = len(self._marks) - n_err
        self.issues_title.configure(
            text=f"Issues: {n_err} error(s), {n_warn} warning(s)")

    def _picked(self, _event=None) -> None:
        sel = self.tree.selection()
        if not sel:
            return
        iid = sel[0]
        i, _, j = iid[1:].partition("_")
        self.focus_on(int(i), int(j) if j else 0)

    # -- state --------------------------------------------------------------

    def set_mesh(self, mesh) -> None:
        """Show this mesh of the input, or none (marks only)."""
        if mesh is None or mesh[0] is None or len(mesh[1]) == 0:
            self._verts = self._tris = self._normals = self._owner = None
            self._face_lo = self._face_hi = None
        else:
            verts, tris, owner, face_lo, face_hi = mesh
            self._verts = np.asarray(verts, dtype=float)
            self._tris = np.asarray(tris, dtype=np.int64)
            self._normals = view3d.triangle_normals(self._verts, self._tris)
            self._owner = None if owner is None else np.asarray(owner, dtype=np.int64)
            self._face_lo = None if face_lo is None else np.asarray(face_lo, dtype=float)
            self._face_hi = None if face_hi is None else np.asarray(face_hi, dtype=float)
        self._place()
        self._match_faces()
        self._schedule_redraw()

    def set_marks(self, marks) -> None:
        """Show these marks -- all of them, replacing what was shown."""
        self._marks = list(marks)
        if self._detail is not None and self._detail[0] >= len(self._marks):
            self._detail = None
        if self._verts is None:
            self._place()
        self._match_faces()
        self._fill_tree()
        self._schedule_redraw()

    def _place(self) -> None:
        """Take the mesh to the bed frame, as the small view does, and find the
        point the view turns about and the sphere it fits."""
        if self._verts is not None:
            self._bed, self._shift = view3d.place_on_bed(self._verts, self._R)
            pts = self._bed
        else:
            self._bed, self._shift = None, np.zeros(3)
            pts = self._mark_points()
        if len(pts) == 0:
            self._pivot, self._extent = np.zeros(3), 1.0
            return
        self._pivot = (pts.min(axis=0) + pts.max(axis=0)) / 2.0
        self._extent = max(float(np.linalg.norm(pts - self._pivot, axis=1).max()), 1e-3)

    def _mark_points(self) -> np.ndarray:
        pts = [self.to_bed(p.points) for m in self._marks for p in m.places]
        return np.vstack(pts) if pts else np.zeros((0, 3))

    def _match_faces(self) -> None:
        self._face_match = {}
        if self._owner is None or self._face_lo is None:
            return
        for i, mark in enumerate(self._marks):
            for j, place in enumerate(mark.places):
                if place.kind != "face":
                    continue
                k = _marks.match_face(place.points[1], place.points[2],
                                      self._face_lo, self._face_hi)
                if k is not None:
                    self._face_match[(i, j)] = k

    def to_bed(self, points) -> np.ndarray:
        """Part-space points in the frame the mesh is drawn in."""
        return np.atleast_2d(np.asarray(points, dtype=float)) @ self._R.T + self._shift

    # -- the camera -----------------------------------------------------------

    def _size(self) -> tuple[int, int]:
        return max(self.canvas.winfo_width(), 2), max(self.canvas.winfo_height(), 2)

    def fit(self) -> None:
        """Fit the whole part in the view, keeping its rotation."""
        w, h = self._size()
        self._scale = view3d.fit_scale(self._extent, w, h)
        self._pan = (0.0, 0.0)
        self._fitted = True
        self.redraw()

    def _zoom(self, factor: float, cursor) -> None:
        w, h = self._size()
        self._scale, self._pan = view3d.zoom_about(
            self._scale, self._pan, factor, cursor, (w / 2.0, h / 2.0))
        self.redraw()

    def screen_of(self, bed_points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Screen positions and depths of bed-frame points, as drawn now."""
        w, h = self._size()
        turned = view3d.view_points(bed_points, self._view, self._pivot)
        xy, depth = view3d.project(turned)
        c, _ = view3d.project(self._pivot[None])
        sx = w / 2.0 + self._pan[0] + (xy[:, 0] - c[0, 0]) * self._scale
        sy = h / 2.0 + self._pan[1] - (xy[:, 1] - c[0, 1]) * self._scale
        return np.column_stack([sx, sy]), depth

    def focus_on(self, i: int, j: int = 0) -> None:
        """Centre the view on place ``j`` of mark ``i`` and open its detail box."""
        if not (0 <= i < len(self._marks)) or not self._marks[i].places:
            return
        j = min(max(j, 0), len(self._marks[i].places) - 1)
        if not self._fitted:
            self.fit()
        w, h = self._size()
        at, _ = self.screen_of(self.to_bed(self._marks[i].places[j].anchor()))
        self._pan = (self._pan[0] + w / 2.0 - at[0, 0], self._pan[1] + h / 2.0 - at[0, 1])
        self._detail = (i, j)
        self.redraw()

    # -- the mouse ------------------------------------------------------------

    def _left_press(self, event) -> None:
        self.canvas.focus_set()
        self._press = (event.x, event.y)
        self._last = (event.x, event.y)
        self._mode = None

    def _left_motion(self, event) -> None:
        if self._press is None:
            return
        if self._mode is None:
            if max(abs(event.x - self._press[0]), abs(event.y - self._press[1])) <= CLICK_SLOP:
                return
            self._mode = "rotate"
        dx, dy = event.x - self._last[0], event.y - self._last[1]
        self._last = (event.x, event.y)
        self._view = view3d.drag_rotation(self._view, dx, dy)
        self.redraw(dragging=True)

    def _left_release(self, event) -> None:
        if self._press is None:
            return
        mode, self._press, self._mode = self._mode, None, None
        if mode == "rotate":
            self.redraw()
            return
        self.click(event.x, event.y)

    def click(self, x: float, y: float) -> None:
        """A left click: open the detail box of the mark under it, or close it."""
        hit = None
        items = self.canvas.find_overlapping(x - 3, y - 3, x + 3, y + 3)
        for item in reversed(items):
            for tag in self.canvas.gettags(item):
                if tag.startswith("place:"):
                    i, j = tag[6:].split(":")
                    hit = (int(i), int(j))
                    break
            if hit is not None:
                break
        self._detail = hit
        self.redraw()

    def _middle_press(self, event) -> None:
        self._last = (event.x, event.y)

    def _middle_motion(self, event) -> None:
        if self._last is None:
            return
        dx, dy = event.x - self._last[0], event.y - self._last[1]
        self._last = (event.x, event.y)
        self._pan = (self._pan[0] + dx, self._pan[1] + dy)
        self.redraw(dragging=True)

    def _middle_release(self, _event) -> None:
        self._last = None
        self.redraw()

    def _wheel(self, event) -> None:
        notches = event.delta / 120.0 if event.delta else 0.0
        if notches:
            self._zoom(WHEEL_FACTOR ** notches, (event.x, event.y))

    def _escape(self, _event=None) -> None:
        if self._detail is not None:
            self._detail = None
            self.redraw()
        else:
            self.close()

    def close(self) -> None:
        if self._redraw_pending is not None:
            try:
                self.after_cancel(self._redraw_pending)
            except tk.TclError:
                pass
            self._redraw_pending = None
        callback, self._on_close = self._on_close, None
        self.destroy()
        if callback is not None:
            callback()

    # -- drawing --------------------------------------------------------------

    def _schedule_redraw(self) -> None:
        if self._redraw_pending is None:
            self._redraw_pending = self.after(REDRAW_DEBOUNCE_MS, self._debounced)

    def _debounced(self) -> None:
        self._redraw_pending = None
        if not self._fitted and self.canvas.winfo_width() > 2:
            self.fit()
        else:
            self.redraw()

    def redraw(self, dragging: bool = False) -> None:
        c = self.canvas
        c.delete("all")
        w, h = self._size()
        if self._bed is not None:
            self._draw_mesh(dragging)
        anchors = self._draw_marks()
        self._draw_axes(h)
        c.create_text(10, 10, anchor="nw", text=CONTROLS_TEXT, fill="#8a8f96",
                      font=("TkDefaultFont", 8), tags=("controls",))
        if self._detail is not None and self._detail in anchors:
            self._draw_detail(*self._detail, anchors[self._detail], w, h)

    def _draw_mesh(self, dragging: bool) -> None:
        keep = (view3d.thin(self._tris, DRAG_TRIANGLES) if dragging
                else np.arange(len(self._tris)))
        tris = self._tris[keep]
        colours = np.tile(np.asarray(BODY_COLOUR, dtype=float), (len(tris), 1))
        if self.overhangs_var.get():
            tint = view3d.overhang_mask(self._normals[keep], self._orient)
            colours[tint] = OVERHANG_COLOUR
        tags = None
        if self._face_match and self._owner is not None:
            owner = self._owner[keep]
            tags = [()] * len(tris)
            # Warnings first, failures over them: red wins where both name a face.
            for (i, j), k in sorted(self._face_match.items(),
                                    key=lambda kv: self._marks[kv[0][0]].is_error):
                hit = np.nonzero(owner == k)[0]
                colours[hit] = self._marks[i].colour
                for t in hit:
                    tags[t] = ("mark", f"place:{i}:{j}")
        screen, depth = self.screen_of(self._bed)
        normals_view = self._normals[keep] @ self._R.T @ self._view.T
        draw_mesh(self.canvas, screen, depth, tris, normals_view, colours, tags=tags)

    def _draw_marks(self) -> dict:
        """Every mark on top of the mesh; returns each place's screen anchor."""
        c = self.canvas
        anchors = {}
        order = sorted(range(len(self._marks)), key=lambda i: self._marks[i].is_error)
        for i in order:
            mark = self._marks[i]
            fill = _hex(mark.colour)
            for j, place in enumerate(mark.places):
                tag = ("mark", f"place:{i}:{j}")
                pts, _ = self.screen_of(self.to_bed(place.points))
                box = place.box()
                if place.kind == "box" or (place.kind == "face"
                                           and (i, j) not in self._face_match):
                    self._draw_box(box[0], box[1], fill, tag)
                if place.kind == "segment":
                    c.create_line(*pts[0], *pts[1], fill=fill, width=3, tags=tag)
                    c.create_line(*pts[0], *pts[1], fill="#333", width=1, tags=tag)
                    for p in pts:
                        self._dot(p, fill, tag)
                at, _ = self.screen_of(self.to_bed(place.anchor()))
                if place.kind != "segment":
                    self._dot(at[0], fill, tag, ring=place.kind != "point")
                anchors[(i, j)] = at[0]
        return anchors

    def _dot(self, p, fill: str, tag, ring: bool = False) -> None:
        r = POINT_RADIUS + (1 if ring else 0)
        self.canvas.create_oval(p[0] - r, p[1] - r, p[0] + r, p[1] + r, fill=fill,
                                outline="#ffffff" if ring else "#222222",
                                width=2 if ring else 1, tags=tag)

    def _draw_box(self, lo, hi, fill: str, tag) -> None:
        lo, hi = np.asarray(lo, dtype=float), np.asarray(hi, dtype=float)
        corners = np.array([[(lo, hi)[(k >> a) & 1][a] for a in range(3)] for k in range(8)])
        s, _ = self.screen_of(self.to_bed(corners))
        for a in range(8):
            for axis in range(3):
                b = a | (1 << axis)
                if b != a:
                    self.canvas.create_line(*s[a], *s[b], fill=fill, width=2, tags=tag)

    def _draw_axes(self, h: int) -> None:
        """The part's X, Y and Z in the bottom-left corner, turned as it is."""
        origin = np.array([46.0, h - 46.0])
        for k, name in enumerate("XYZ"):
            d = self._view @ self._R[:, k]
            xy, _ = view3d.project(d[None])
            tip = origin + np.array([xy[0, 0], -xy[0, 1]]) * 32.0
            self.canvas.create_line(*origin, *tip, fill=AXIS_COLOURS[k], width=2,
                                    arrow="last", tags=("axes",))
            self.canvas.create_text(tip[0] + 6, tip[1] - 6, text=name,
                                    fill=AXIS_COLOURS[k], tags=("axes",),
                                    font=("TkDefaultFont", 8, "bold"))

    def _draw_detail(self, i: int, j: int, at, w: int, h: int) -> None:
        mark = self._marks[i]
        c = self.canvas
        text = c.create_text(0, 0, anchor="nw", width=340, fill="#222",
                             text=_marks.detail_text(mark, mark.places[j]),
                             font=("TkDefaultFont", 9), tags=("detail",))
        x0, y0, x1, y1 = c.bbox(text)
        bw, bh = x1 - x0 + 20, y1 - y0 + 16
        x = at[0] + 14 if at[0] + 14 + bw < w else at[0] - 14 - bw
        y = min(max(at[1] - bh / 2.0, 6), max(h - bh - 6, 6))
        c.coords(text, x + 10, y + 8)
        box = c.create_rectangle(x, y, x + bw, y + bh, fill="white",
                                 outline=_hex(mark.colour), width=2, tags=("detail",))
        c.tag_lower(box, text)

    # -- for tests and the window -----------------------------------------

    @property
    def detail(self):
        """``(mark, place)`` whose detail box is open, or ``None``."""
        return self._detail

    @property
    def view_rotation(self) -> np.ndarray:
        return self._view.copy()

    @property
    def camera(self) -> tuple:
        return self._scale, self._pan
