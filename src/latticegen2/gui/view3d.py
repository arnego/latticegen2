"""The orientation view's mathematics. No Tk in this module.

Split from :mod:`latticegen2.gui.orientview` for the reason the front-end is
split everywhere: a projection, a depth sort and an overhang mask can be tested
without a display, and a canvas cannot.

Everything is in **bed coordinates** — the bed is the ``z = 0`` plane and its
normal is ``+z`` — because that is the frame the user is looking at. The part is
taken there with the same rotation the pipeline uses
(:func:`latticegen2.orient.rotation_matrix`), so what the view shows and what a
run does cannot be two different conventions.
"""

from __future__ import annotations

import numpy as np

from .. import orient as _orient
from ..cli import OVERHANG_MARGIN

#: Direction from the scene toward the eye: a fixed three-quarter view from
#: above. Fixed on purpose — the thing being manipulated is the part, and a
#: camera that also moved would make the same drag mean two things.
EYE = np.array([1.0, -1.35, 0.95]) / np.linalg.norm([1.0, -1.35, 0.95])
RIGHT = np.cross([0.0, 0.0, 1.0], EYE)
RIGHT = RIGHT / np.linalg.norm(RIGHT)
UP = np.cross(EYE, RIGHT)

#: Where the light comes from, for flat shading.
LIGHT = np.array([0.35, -0.5, 0.8]) / np.linalg.norm([0.35, -0.5, 0.8])

DRAG_DEGREES_PER_PIXEL = 0.6


def project(points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Orthographic ``(xy, depth)`` of bed-space points; larger depth is nearer."""
    points = np.asarray(points, dtype=float)
    return np.column_stack([points @ RIGHT, points @ UP]), points @ EYE


def place_on_bed(verts: np.ndarray, R: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotate part-space vertices onto the bed.

    Returns ``(bed_verts, shift)``: the part is turned by ``R``, then moved so
    it rests on the bed (lowest point at ``z = 0``) centred over the origin.
    ``shift`` is that translation, so the part's own origin can be drawn too.
    """
    turned = np.asarray(verts, dtype=float) @ R.T
    lo, hi = turned.min(axis=0), turned.max(axis=0)
    shift = np.array([-(lo[0] + hi[0]) / 2.0, -(lo[1] + hi[1]) / 2.0, -lo[2]])
    return turned + shift, shift


def triangle_normals(verts: np.ndarray, tris: np.ndarray) -> np.ndarray:
    """Unit normals of triangles, zero for a degenerate one."""
    a, b, c = verts[tris[:, 0]], verts[tris[:, 1]], verts[tris[:, 2]]
    n = np.cross(b - a, c - a)
    mag = np.linalg.norm(n, axis=1, keepdims=True)
    return np.where(mag > 0.0, n / np.where(mag > 0.0, mag, 1.0), 0.0)


def overhang_mask(normals_part: np.ndarray, orient, overhang: float) -> np.ndarray:
    """Which preview triangles are inner overhangs at this orientation.

    The same test a run applies (docs/algorithm.md §14.1): the outward normal
    leans toward the build direction by more than the limit allows, with the
    run's own margin. The preview's triangles are coarse, so this is an
    indication of what will be supported, not the measurement the run makes on
    the true surface.
    """
    b = _orient.build_direction(orient)
    bar = float(np.sin(np.radians(overhang - OVERHANG_MARGIN)))
    return np.asarray(normals_part) @ b >= bar


def shade(normals_bed: np.ndarray) -> np.ndarray:
    """Flat-shading intensity in ``[0.35, 1]`` per triangle."""
    return 0.35 + 0.65 * np.clip(np.abs(np.asarray(normals_bed) @ LIGHT), 0.0, 1.0)


def drag_rotation(R: np.ndarray, dx: float, dy: float) -> np.ndarray:
    """The part's rotation after the mouse moved ``(dx, dy)`` pixels.

    Horizontal movement turns the part about the bed's normal, vertical
    movement tips it about the screen's horizontal axis — both fixed in the bed
    frame, so a drag does what it looks like it does whichever way the part is
    already turned.
    """
    yaw = np.radians(dx * DRAG_DEGREES_PER_PIXEL)
    tip = np.radians(dy * DRAG_DEGREES_PER_PIXEL)
    return _axis_rotation(RIGHT, tip) @ _axis_rotation(np.array([0.0, 0.0, 1.0]), yaw) @ R


def _axis_rotation(axis: np.ndarray, angle: float) -> np.ndarray:
    x, y, z = axis / np.linalg.norm(axis)
    c, s = np.cos(angle), np.sin(angle)
    k = np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
    return np.eye(3) + s * k + (1.0 - c) * (k @ k)


def snapped_angles(R: np.ndarray, step: float = 0.1) -> tuple[float, float, float]:
    """Euler angles for ``R``, rounded to what the three boxes display."""
    out = []
    for a in _orient.euler_from_matrix(R):
        a = round(a / step) * step
        out.append(0.0 if abs(a) < step / 2.0 else float(round(a, 6)))
    return tuple(out)


def thin(tris: np.ndarray, limit: int) -> np.ndarray:
    """At most ``limit`` triangles, evenly taken — the copy drawn during a drag."""
    if len(tris) <= limit:
        return np.arange(len(tris))
    return np.linspace(0, len(tris) - 1, limit).astype(np.int64)
