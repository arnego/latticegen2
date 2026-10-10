"""Overhang detection on the input body (docs/algorithm.md §14.1).

A point of the body's surface is an inner overhang when its outward normal
``n`` satisfies ``n . b > sin(phi)``: the surface acts as a ceiling of the
volume, and the shell beyond it overhangs the lattice core more flatly than the
limit ``phi`` (measured from the build direction ``b``). A flat ceiling has
``n = b`` and always qualifies; a vertical wall never does. The underside of
the body is an *outside* overhang and is never examined
(specification.md §4.6).

Three kinds of sample come out of this module, and only the first is what a
face-normal test alone would find:

* **faces** — sampled on the *true* surface with its exact normal, not on the
  classification mesh's triangles, whose angular deflection would smear the
  boundary of a region by more than ten degrees;
* **edges** — a down-pointing ridge of shell material (a concave edge of the
  body) whose two faces are each steep enough but which itself runs flatter
  than the limit: a keel hanging in mid-air;
* **vertices** — a down-pointing tip, the lowest point of a stalactite.

Everything flatter than ``phi - margin`` is reported. The margin is
conservative by design (more support, never less) and it keeps a support face
from meeting the wall tangentially at the border of an overhang region, which
is the grazing-trim regime docs/algorithm.md §7 documents.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.GCPnts import GCPnts_AbscissaPoint
from OCP.GeomAbs import GeomAbs_SurfaceType
from OCP.TopAbs import TopAbs_Orientation, TopAbs_ShapeEnum, TopAbs_State
from OCP.TopExp import TopExp
from OCP.TopoDS import TopoDS
from OCP.TopTools import TopTools_IndexedDataMapOfShapeListOfShape
from OCP.gp import gp_Pnt, gp_Vec

from .. import occ
from ..cli import OVERHANG_MARGIN as MARGIN_DEG

FACE, EDGE, VERTEX = 0, 1, 2

PREFILTER_SLACK = 0.25
"""How far below the threshold a mesh node's ``n . b`` may sit while its
triangle is still refined. The classification mesh allows 0.2 rad of angular
deflection inside one triangle, so a triangle whose corners all read below the
bar can still cross it in the middle."""


@dataclass
class OverhangSamples:
    """What has to be covered by support."""

    points: np.ndarray
    """``(N, 3)`` points on the body's surface."""
    inward: np.ndarray
    """``(N, 3)`` unit direction from each point into the body — ``-n`` for a
    face sample, ``-b`` for an edge or vertex sample. Moving a point a hair
    along it gives a location that must end up inside support."""
    kind: np.ndarray
    """``(N,)`` :data:`FACE`, :data:`EDGE` or :data:`VERTEX`."""
    area: float
    """Surface area the face samples stand for, in mm²."""
    pitch: float

    def __len__(self) -> int:
        return len(self.points)

    def counts(self) -> dict:
        return {
            "face": int(np.count_nonzero(self.kind == FACE)),
            "edge": int(np.count_nonzero(self.kind == EDGE)),
            "vertex": int(np.count_nonzero(self.kind == VERTEX)),
        }


def sample_pitch(t: float, thickness: float) -> float:
    """Spacing of the surface samples, in mm (docs/algorithm.md §14.1)."""
    return float(min(max(min(t, thickness) / 2.0, 0.25), 1.0))


def threshold(overhang_deg: float, margin_deg: float = MARGIN_DEG) -> float:
    """``sin(phi - margin)``: the bar ``n . b`` is compared against."""
    return float(np.sin(np.radians(overhang_deg - margin_deg)))


class _SurfaceNormals:
    """Exact outward normals of one face of the body."""

    def __init__(self, face):
        self.surf = BRepAdaptor_Surface(face)
        self.flip = -1.0 if face.Orientation() == TopAbs_Orientation.TopAbs_REVERSED else 1.0
        self._p = gp_Pnt()
        self._du = gp_Vec()
        self._dv = gp_Vec()

    def at(self, u: float, v: float):
        """``(point, unit outward normal)``, the normal ``None`` at a pole."""
        self.surf.D1(float(u), float(v), self._p, self._du, self._dv)
        n = self._du.Crossed(self._dv)
        mag = n.Magnitude()
        p = np.array([self._p.X(), self._p.Y(), self._p.Z()])
        if mag < 1e-12:
            return p, None
        return p, self.flip * np.array([n.X(), n.Y(), n.Z()]) / mag


def _subdivide(n_long: int, n_high: int):
    """Barycentric weights filling a triangle with ``n_long * n_high`` samples.

    Rows run parallel to the triangle's longest edge (between its first two
    corners) and are spaced evenly toward the third. Sampling the two
    directions separately is what keeps a long thin triangle — and the mesher
    makes many of them on a cylinder or a large flat face — from being sampled
    at the square of its long side.
    """
    u = (np.arange(n_long) + 0.5) / n_long
    v = (np.arange(n_high) + 0.5) / n_high
    uu, vv = np.meshgrid(u, v, indexing="ij")
    s = (uu * (1.0 - vv)).ravel()
    t = vv.ravel()
    # Rows shorten toward the third corner while holding the same number of
    # samples, so a sample there stands for less area: its share is the row's
    # width, (1 - v).
    share = (1.0 - vv).ravel()
    return np.column_stack([1.0 - s - t, s, t]), share / share.sum()


def _face_samples(face, b: np.ndarray, bar: float, pitch: float):
    tri = occ.face_uv_triangulation(face)
    if tri is None:
        return None
    pts, uvs, tris, _evaluate = tri
    normals = _SurfaceNormals(face)
    planar = normals.surf.GetType() == GeomAbs_SurfaceType.GeomAbs_Plane

    if planar:
        _p, n = normals.at(*uvs[0])
        if n is None or float(n @ b) < bar:
            return None
        node_dot = None
    else:
        node_dot = np.full(len(uvs), -2.0)
        for i, (u, v) in enumerate(uvs):
            _p, nn = normals.at(u, v)
            # A pole has no normal of its own; let it pass the prefilter and be
            # decided by the samples inside its triangles.
            node_dot[i] = 2.0 if nn is None else float(nn @ b)

    out_p, out_n = [], []
    area = 0.0
    A, B, C = pts[tris[:, 0]], pts[tris[:, 1]], pts[tris[:, 2]]
    tri_area = 0.5 * np.linalg.norm(np.cross(B - A, C - A), axis=1)
    sides = np.sqrt(np.stack([
        np.einsum("ij,ij->i", B - A, B - A),
        np.einsum("ij,ij->i", C - B, C - B),
        np.einsum("ij,ij->i", A - C, A - C),
    ], axis=1))
    # Rotate each triangle so its longest side joins its first two corners.
    roll = np.argmax(sides, axis=1)
    for ti in range(len(tris)):
        if tri_area[ti] <= 0.0:
            continue
        if node_dot is not None and float(node_dot[tris[ti]].max()) < bar - PREFILTER_SLACK:
            continue
        corners = np.roll(tris[ti], -int(roll[ti]))
        base = float(sides[ti, roll[ti]])
        w, share = _subdivide(
            max(1, int(np.ceil(base / pitch))),
            max(1, int(np.ceil(2.0 * tri_area[ti] / base / pitch))),
        )
        if planar:
            out_p.append(w @ pts[corners])
            out_n.append(np.broadcast_to(n, (len(w), 3)))
            area += tri_area[ti]
            continue
        for (u, v), part in zip(w @ uvs[corners], share):
            p, nn = normals.at(u, v)
            if nn is not None and float(nn @ b) >= bar:
                out_p.append(p[None, :])
                out_n.append(nn[None, :])
                area += tri_area[ti] * part
    if not out_p:
        return None
    return np.vstack(out_p), np.vstack(out_n), area


def _edge_orientation_in(face, edge):
    for e in occ._explore(face, TopAbs_ShapeEnum.TopAbs_EDGE):
        if e.IsSame(edge):
            return e.Orientation()
    return None


def _edge_samples(body, b: np.ndarray, bar: float, pitch: float) -> list:
    """Points on keel edges: concave edges of the body flatter than the limit.

    The outward normals between the two faces of an edge sweep an arc. Where
    the shell protrudes into the volume along that edge (the body is concave
    there) the shell's surface really does face every direction on the arc, so
    the edge overhangs iff the arc contains a direction with ``n . b`` over the
    bar while neither face reaches it on its own — the faces are then covered
    as faces. Equivalently, the edge runs flatter than the limit and faces
    down.
    """
    owners = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(
        body, TopAbs_ShapeEnum.TopAbs_EDGE, TopAbs_ShapeEnum.TopAbs_FACE, owners
    )
    out = []
    for i in range(1, owners.Extent() + 1):
        edge = TopoDS.Edge_s(owners.FindKey(i))
        if BRep_Tool.Degenerated_s(edge):
            continue
        faces = [TopoDS.Face_s(f) for f in owners.FindFromIndex(i)]
        if len(faces) != 2 or faces[0].IsSame(faces[1]):
            continue
        ori = _edge_orientation_in(faces[0], edge)
        if ori is None:
            continue
        sign = -1.0 if ori == TopAbs_Orientation.TopAbs_REVERSED else 1.0
        curve = BRepAdaptor_Curve(edge)
        first, last = curve.FirstParameter(), curve.LastParameter()
        try:
            length = GCPnts_AbscissaPoint.Length_s(curve)
        except Exception:                                      # noqa: BLE001
            continue
        if not np.isfinite(length) or length <= 0.0:
            continue
        pcurves = []
        for f in faces:
            c2d = BRep_Tool.CurveOnSurface_s(edge, f, 0.0, 0.0)
            pcurves.append(c2d)
        if any(c is None for c in pcurves):
            continue
        normals = [_SurfaceNormals(f) for f in faces]
        steps = max(2, int(np.ceil(length / pitch)) + 1)
        p3, d3 = gp_Pnt(), gp_Vec()
        for s in np.linspace(first, last, steps):
            uv1 = pcurves[0].Value(float(s))
            uv2 = pcurves[1].Value(float(s))
            _p1, n1 = normals[0].at(uv1.X(), uv1.Y())
            _p2, n2 = normals[1].at(uv2.X(), uv2.Y())
            if n1 is None or n2 is None:
                continue
            if float(n1 @ b) >= bar or float(n2 @ b) >= bar:
                continue                      # covered as a face already
            curve.D1(float(s), p3, d3)
            if d3.Magnitude() < 1e-12:
                continue
            tangent = sign * np.array([d3.X(), d3.Y(), d3.Z()]) / d3.Magnitude()
            # `into` points across face 0 away from the edge. The body is
            # concave along the edge when face 1's outward normal leans toward
            # that direction rather than away from it.
            into = np.cross(n1, tangent)
            if float(into @ n2) <= 1e-6:
                continue
            axis = np.cross(n1, n2)
            mag = float(np.linalg.norm(axis))
            if mag < 1e-9:
                continue
            axis /= mag
            proj = b - float(b @ axis) * axis
            reach = float(np.linalg.norm(proj))
            if reach < bar:
                continue
            proj /= reach
            if float(np.cross(n1, proj) @ axis) < 0.0 or float(np.cross(proj, n2) @ axis) < 0.0:
                continue                      # steepest direction is off the arc
            out.append((p3.X(), p3.Y(), p3.Z()))
    return out


def _vertex_samples(body, b: np.ndarray) -> list:
    """Down-pointing tips: vertices lower than everything attached to them, with
    the body directly beneath."""
    owners = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(
        body, TopAbs_ShapeEnum.TopAbs_VERTEX, TopAbs_ShapeEnum.TopAbs_EDGE, owners
    )
    lo, hi = occ.bounding_box(body)
    eps = 1e-4 * float(np.linalg.norm(hi - lo))
    out = []
    classifier = None
    for i in range(1, owners.Extent() + 1):
        vertex = TopoDS.Vertex_s(owners.FindKey(i))
        pv = BRep_Tool.Pnt_s(vertex)
        v = np.array([pv.X(), pv.Y(), pv.Z()])
        rises = 0
        tip = True
        for e in owners.FindFromIndex(i):
            edge = TopoDS.Edge_s(e)
            if BRep_Tool.Degenerated_s(edge):
                continue
            curve = BRepAdaptor_Curve(edge)
            first, last = curve.FirstParameter(), curve.LastParameter()
            ends = [curve.Value(first), curve.Value(last)]
            near = min((0, 1), key=lambda k: ends[k].Distance(pv))
            s = first + 0.05 * (last - first) if near == 0 else last - 0.05 * (last - first)
            q = curve.Value(float(s))
            if float((np.array([q.X(), q.Y(), q.Z()]) - v) @ b) <= eps * 1e-3:
                tip = False
                break
            rises += 1
        if not tip or rises == 0:
            continue
        below = v - eps * b
        if classifier is None:
            classifier = BRepClass3d_SolidClassifier(body)
        classifier.Perform(gp_Pnt(*below), 1e-7)
        if classifier.State() == TopAbs_State.TopAbs_IN:
            out.append(tuple(v))
    return out


def detect(body, b: np.ndarray, overhang_deg: float, pitch: float,
           margin_deg: float = MARGIN_DEG, report=None) -> OverhangSamples:
    """Every inner overhang of ``body`` for build direction ``b``.

    ``body`` must already carry the classification triangulation
    (:func:`latticegen2.classify.tessellate_surface`), whose UV nodes are what
    the face samples are refined from.
    """
    b = np.asarray(b, dtype=float)
    bar = threshold(overhang_deg, margin_deg)
    points, inward, kinds = [], [], []
    area = 0.0
    faces = occ.faces(body)
    for fi, face in enumerate(faces):
        got = _face_samples(face, b, bar, pitch)
        if report is not None:
            report("sampling overhangs", fi + 1, len(faces))
        if got is None:
            continue
        p, n, a = got
        points.append(p)
        inward.append(-n)
        kinds.append(np.full(len(p), FACE))
        area += a

    for kind, pts in ((EDGE, _edge_samples(body, b, bar, pitch)),
                      (VERTEX, _vertex_samples(body, b))):
        if pts:
            arr = np.array(pts, dtype=float)
            points.append(arr)
            inward.append(np.broadcast_to(-b, arr.shape).copy())
            kinds.append(np.full(len(arr), kind))

    if not points:
        return OverhangSamples(np.empty((0, 3)), np.empty((0, 3)),
                               np.empty(0, dtype=np.int64), 0.0, pitch)
    return OverhangSamples(np.vstack(points), np.vstack(inward),
                           np.concatenate(kinds).astype(np.int64), area, pitch)
