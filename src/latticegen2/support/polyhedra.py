"""Convex support primitives and their analytic clip (docs/algorithm.md §14.3, §14.5).

Pure NumPy. A primitive is a convex polyhedron held as a vertex array plus one
index loop per face, wound so its normal points outward. Two kinds exist:

* a **capital** — the upward 8-sided pyramid ``K(s)`` truncated at a height;
* a **post** — ``K(q)`` intersected with a vertical octagonal prism, i.e. a
  pencil standing on its tip.

Both have their **ridge edges at exactly the maximum overhang angle** from the
build direction, so their faces are slightly steeper and neither a face nor an
edge of a support can exceed the limit (specification.md §4.6).

Clipping a primitive to one lattice cell is a convex polyhedron against six
planes, done here rather than in the kernel: it is exact to floating point and
costs microseconds, and it is what keeps every boolean in the support stage
local to one junction (§14.5).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

SIDES = 8
"""Sides of every support cross-section (specification.md §4.6: faceted)."""

INSCRIBED = float(np.cos(np.pi / SIDES))
"""Inradius over circumradius of the regular octagon, ``cos(pi/8)``.

A sample is only counted as covered when it lies inside the *inscribed* circle
of a capital, so coverage holds whatever azimuth phase the capital was given."""

POST_AREA_FACTOR = float(SIDES * np.tan(np.pi / SIDES) / 4.0)
"""Cross-section area of an octagon over the square of its across-flats width."""

CLIP_EPS = 1e-10
"""Millimetres. A vertex this close to a clip plane is treated as lying on it."""


@dataclass
class Primitive:
    """One convex support solid, before it is clipped to lattice cells."""

    kind: str
    """``"capital"`` or ``"post"``."""
    verts: np.ndarray
    faces: list
    anchor: str
    """What the apex is rooted in: ``"strut"``, ``"node"`` or ``"wall"``."""


def _ring(centre: np.ndarray, radius: float, phase: float, hx: np.ndarray,
          hy: np.ndarray) -> np.ndarray:
    ang = phase + np.arange(SIDES) * (2.0 * np.pi / SIDES)
    return centre + radius * (np.cos(ang)[:, None] * hx + np.sin(ang)[:, None] * hy)


def orient_outward(verts: np.ndarray, faces: list) -> list:
    """Wind every face of a convex polyhedron so its normal points outward."""
    centre = verts.mean(axis=0)
    out = []
    for loop in faces:
        pts = verts[loop]
        normal = np.zeros(3)
        for i in range(len(pts)):
            normal += np.cross(pts[i], pts[(i + 1) % len(pts)])
        if float(np.dot(normal, pts.mean(axis=0) - centre)) < 0.0:
            loop = loop[::-1]
        out.append(list(loop))
    return out


def capital(apex: np.ndarray, height: float, tau: float, phase: float,
            frame, anchor: str) -> Primitive:
    """``K(apex)`` truncated ``height`` above its apex.

    ``frame`` is ``(hx, hy, b)``: two horizontal axes and the build direction,
    all in part coordinates. The ring's *circumradius* is ``height * tau``, so
    the eight ridge edges recline exactly ``atan(tau)`` from ``b``.
    """
    hx, hy, b = frame
    apex = np.asarray(apex, dtype=float)
    ring = _ring(apex + height * b, height * tau, phase, hx, hy)
    verts = np.vstack([apex[None, :], ring])
    faces = [[0, 1 + k, 1 + (k + 1) % SIDES] for k in range(SIDES)]
    faces.append(list(range(1, 1 + SIDES)))
    return Primitive("capital", verts, orient_outward(verts, faces), anchor)


def post(base: np.ndarray, length: float, across_flats: float, tau: float,
         phase: float, frame, anchor: str) -> Primitive:
    """A vertical pencil: ``K(base)`` inside an octagonal prism.

    ``across_flats`` is the minimum support thickness, so no section of the
    prism is thinner than it. The tip is a cone at the overhang limit, which is
    what lets the post start from a point buried in its anchor instead of from
    a flat face that would itself overhang.
    """
    hx, hy, b = frame
    base = np.asarray(base, dtype=float)
    radius = across_flats / (2.0 * INSCRIBED)
    shoulder = radius / tau
    length = max(float(length), shoulder * 1.05)
    lower = _ring(base + shoulder * b, radius, phase, hx, hy)
    upper = _ring(base + length * b, radius, phase, hx, hy)
    verts = np.vstack([base[None, :], lower, upper])
    faces = [[0, 1 + k, 1 + (k + 1) % SIDES] for k in range(SIDES)]
    for k in range(SIDES):
        k2 = (k + 1) % SIDES
        faces.append([1 + k, 1 + k2, 1 + SIDES + k2, 1 + SIDES + k])
    faces.append(list(range(1 + SIDES, 1 + 2 * SIDES)))
    return Primitive("post", verts, orient_outward(verts, faces), anchor)


def volume(verts: np.ndarray, faces: list) -> float:
    """Volume of a closed, outward-wound polyhedron (divergence theorem)."""
    total = 0.0
    for loop in faces:
        p0 = verts[loop[0]]
        for i in range(1, len(loop) - 1):
            total += float(np.dot(p0, np.cross(verts[loop[i]], verts[loop[i + 1]])))
    return total / 6.0


def clip_halfspace(verts: np.ndarray, faces: list, normal: np.ndarray,
                   offset: float):
    """The part of a convex polyhedron on the inner side of one plane.

    Keeps ``normal . x <= offset`` and returns ``(verts, faces)``, or ``None``
    when nothing with volume is left. A vertex within :data:`CLIP_EPS` of the
    plane is *on* it and is never cut again, which is what stops a plane grazing
    a vertex from minting a second, nearly coincident one.
    """
    d = verts @ normal - offset
    inside = d < -CLIP_EPS
    outside = d > CLIP_EPS
    if not outside.any():
        return verts, faces
    if not inside.any():
        return None

    new_verts = [v for v in verts]
    cuts: dict[tuple[int, int], int] = {}

    def cut(i: int, j: int) -> int:
        # Keyed and evaluated on the ordered pair, so the two faces sharing an
        # edge get the same point rather than two near-equal ones.
        a, b_ = (i, j) if i < j else (j, i)
        got = cuts.get((a, b_))
        if got is None:
            t = d[a] / (d[a] - d[b_])
            new_verts.append(verts[a] + t * (verts[b_] - verts[a]))
            got = len(new_verts) - 1
            cuts[(a, b_)] = got
        return got

    new_faces = []
    on_plane: set[int] = set()
    for loop in faces:
        out = []
        n = len(loop)
        for k in range(n):
            i, j = loop[k], loop[(k + 1) % n]
            if not outside[i]:
                out.append(i)
                if not inside[i]:
                    on_plane.add(i)
            if (inside[i] and outside[j]) or (outside[i] and inside[j]):
                c = cut(i, j)
                out.append(c)
                on_plane.add(c)
        if len(out) >= 3:
            new_faces.append(out)

    used = {i for loop in new_faces for i in loop}
    cap = sorted(on_plane & used)
    all_verts = np.array(new_verts)
    if len(cap) >= 3:
        pts = all_verts[cap]
        centre = pts.mean(axis=0)
        ref = (np.array([1.0, 0.0, 0.0]) if abs(normal[0]) < 0.9
               else np.array([0.0, 1.0, 0.0]))
        ax = np.cross(normal, ref)
        ax /= np.linalg.norm(ax)
        ay = np.cross(normal, ax)
        rel = pts - centre
        order = np.argsort(np.arctan2(rel @ ay, rel @ ax))
        new_faces.append([cap[i] for i in order])

    # Compact: drop vertices nothing references any more.
    keep = sorted({i for f in new_faces for i in f})
    remap = {old: new for new, old in enumerate(keep)}
    out_verts = all_verts[keep]
    out_faces = [[remap[i] for i in f] for f in new_faces]
    return out_verts, out_faces


def clip_to_cell(verts: np.ndarray, faces: list, lp, node_pos: np.ndarray,
                 min_volume: float = 0.0):
    """Clip a convex polyhedron to the lattice cell of the node at ``node_pos``.

    The cell is the cube ``|e_k . (x - node)| <= a/2`` whose faces are the
    junction's six cap planes (docs/algorithm.md §14.5). Returns
    ``(verts, faces)`` or ``None`` when the primitive does not reach the cell
    with more than ``min_volume``.
    """
    half = lp.a / 2.0
    for k in range(3):
        e = lp.e[k]
        centre = float(np.dot(e, node_pos))
        for sign in (1.0, -1.0):
            got = clip_halfspace(verts, faces, sign * e, sign * centre + half)
            if got is None:
                return None
            verts, faces = got
    if len(faces) < 4:
        return None
    faces = orient_outward(verts, faces)
    if volume(verts, faces) <= min_volume:
        return None
    return verts, faces


def face_planes(verts: np.ndarray, faces: list) -> tuple[np.ndarray, np.ndarray]:
    """``(normals, offsets)`` of an outward-wound convex polyhedron's faces:
    a point ``x`` is inside iff ``normals @ x <= offsets`` for every face."""
    normals = np.empty((len(faces), 3))
    offsets = np.empty(len(faces))
    for i, loop in enumerate(faces):
        pts = verts[loop]
        n = np.zeros(3)
        for k in range(len(pts)):
            n += np.cross(pts[k], pts[(k + 1) % len(pts)])
        n /= np.linalg.norm(n)
        normals[i] = n
        offsets[i] = float(n @ pts.mean(axis=0))
    return normals, offsets


CONTAIN_TOL = 1e-9
"""Millimetres of slack when asking whether one convex piece contains another."""


def drop_contained(pieces: list) -> list:
    """Remove every convex piece that lies wholly inside another one.

    A convex polyhedron is inside another exactly when all its vertices are, so
    this is exact, and it removes nothing from the union. It matters because
    capitals nest: a deep one largely contains the shallower ones around it,
    and within one lattice cell most of the clipped pieces are then redundant
    operands of a fuse whose cost grows much faster than its operand count
    (docs/algorithm.md §12).
    """
    order = sorted(range(len(pieces)),
                   key=lambda i: volume(*pieces[i]), reverse=True)
    kept: list = []
    planes: list = []
    for i in order:
        verts = pieces[i][0]
        inside = False
        for normals, offsets in planes:
            if bool(((verts @ normals.T) <= offsets[None, :] + CONTAIN_TOL).all()):
                inside = True
                break
        if not inside:
            kept.append(i)
            planes.append(face_planes(*pieces[i]))
    return [pieces[i] for i in sorted(kept)]
