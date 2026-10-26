"""What a first-contact ray can land on (docs/algorithm.md §14.3 S3a).

The input body's faces, tessellated finely enough that a target found on the
mesh is within a fraction of the support step of the true surface, and the
lattice's faces, which are planar and so tessellate exactly. Both are read off
the shapes themselves: the lattice is never rebuilt from the analytic struts,
which would include material the floating-body rule dropped.

**Anchor groups.** Every obstacle triangle also carries the *anchor group* of
its B-rep face: faces of one shape joined across every edge that is smooth or
convex **as seen from the anchor's own material** -- the outside of the input
body, the inside of a lattice solid. A support cell lands on one group only
(§14.3 S3b), and that is what lets its loft close inside the anchor: a surface
spanning a loop on a convex patch stays inside the convex material behind it.
Two walls meeting at a box's inside corner are two groups, since outside the
body that corner is reflex; a tube's two half-cylinders, meeting smoothly, are
one; so are the four faces of a strut, convex from inside the strut.
"""

from __future__ import annotations

import numpy as np
from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Curve
from OCP.BRepBuilderAPI import BRepBuilderAPI_Copy
from OCP.TopAbs import TopAbs_Orientation, TopAbs_ShapeEnum
from OCP.TopExp import TopExp
from OCP.TopoDS import TopoDS
from OCP.TopTools import TopTools_IndexedDataMapOfShapeListOfShape, TopTools_IndexedMapOfShape
from OCP.gp import gp_Pnt, gp_Vec

from .. import occ
from ..classify import TriMesh
from .overhang import _SurfaceNormals, _edge_orientation_in
from .rays import Obstacles

BODY = 0
"""Shape id of the input body; lattice solid ``k`` is ``k + 1``."""

SMOOTH_TOL = 1e-3
"""Radians. Faces meeting at less than this are treated as smooth."""


def _face_groups(shape, material_inside: bool) -> np.ndarray:
    """Anchor group per face of ``shape`` (in ``occ.faces`` order)."""
    faces = occ.faces(shape)
    fmap = TopTools_IndexedMapOfShape()
    for f in faces:
        fmap.Add(f)
    parent = list(range(len(faces)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    owners = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(shape, TopAbs_ShapeEnum.TopAbs_EDGE,
                                   TopAbs_ShapeEnum.TopAbs_FACE, owners)
    for i in range(1, owners.Extent() + 1):
        edge = TopoDS.Edge_s(owners.FindKey(i))
        if BRep_Tool.Degenerated_s(edge):
            continue
        adj = [TopoDS.Face_s(f) for f in owners.FindFromIndex(i)]
        if len(adj) != 2:
            continue
        ia, ib = fmap.FindIndex(adj[0]) - 1, fmap.FindIndex(adj[1]) - 1
        if ia < 0 or ib < 0 or find(ia) == find(ib):
            continue
        if _joins(edge, adj[0], adj[1], material_inside):
            parent[max(find(ia), find(ib))] = min(find(ia), find(ib))
    return np.array([find(i) for i in range(len(faces))], dtype=np.int64)


def _joins(edge, f0, f1, material_inside: bool) -> bool:
    """Smooth, or convex seen from the material, at the edge's midpoint."""
    got = _edge_shape(edge, f0, f1, material_inside)
    return got is None or got[0] < SMOOTH_TOL or got[1]


def _edge_shape(edge, f0, f1, material_inside: bool):
    """``(angle, convex)`` at the edge's midpoint: the angle between the two
    faces' normals, and whether the edge is convex seen from the material;
    ``None`` for a seam or where the faces cannot be evaluated (then the caller
    treats it as smooth)."""
    if f0.IsSame(f1):
        return None                                # a seam
    curve = BRepAdaptor_Curve(edge)
    t = 0.5 * (curve.FirstParameter() + curve.LastParameter())
    p, d = gp_Pnt(), gp_Vec()
    curve.D1(t, p, d)
    if d.Magnitude() < 1e-12:
        return None
    normals = []
    for f in (f0, f1):
        c2d = BRep_Tool.CurveOnSurface_s(edge, f, 0.0, 0.0)
        if c2d is None:
            return (np.pi, False)
        uv = c2d.Value(t)
        _p, n = _SurfaceNormals(f).at(uv.X(), uv.Y())
        if n is None:
            return (np.pi, False)
        normals.append(n)
    n0, n1 = normals
    angle = float(np.arccos(np.clip(float(n0 @ n1), -1.0, 1.0)))
    if angle < SMOOTH_TOL:
        return (angle, True)
    ori = _edge_orientation_in(f0, edge)
    if ori is None:
        return (angle, False)
    sign = -1.0 if ori == TopAbs_Orientation.TopAbs_REVERSED else 1.0
    tangent = sign * np.array([d.X(), d.Y(), d.Z()]) / d.Magnitude()
    # `into` runs across face 0, away from the edge. The solid is concave along
    # the edge -- reflex inside, convex outside -- when face 1's outward normal
    # leans toward it (the test `overhang._edge_samples` makes for keels).
    into = np.cross(n0, tangent)
    concave_inside = float(into @ n1) > 0.0
    # Material inside the solid (a lattice strut): join where it is convex
    # inside. Material outside (the input body's surroundings): join where the
    # body is concave inside, which is convex outside.
    return (angle, concave_inside != material_inside)


def convex_edge_points(shape, material_inside: bool, min_angle_deg: float,
                       spacing: float) -> np.ndarray:
    """Points every ``spacing`` along the edges of ``shape`` that are convex seen
    from the material and sharper than ``min_angle_deg``: a ledge's rim, which
    first-contact rays converge on (docs/algorithm.md §14.3 S3b)."""
    from OCP.GCPnts import GCPnts_AbscissaPoint

    owners = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(shape, TopAbs_ShapeEnum.TopAbs_EDGE,
                                   TopAbs_ShapeEnum.TopAbs_FACE, owners)
    out = []
    for i in range(1, owners.Extent() + 1):
        edge = TopoDS.Edge_s(owners.FindKey(i))
        if BRep_Tool.Degenerated_s(edge):
            continue
        adj = [TopoDS.Face_s(f) for f in owners.FindFromIndex(i)]
        if len(adj) != 2:
            continue
        got = _edge_shape(edge, adj[0], adj[1], material_inside)
        if got is None or not got[1] or got[0] < np.radians(min_angle_deg):
            continue
        curve = BRepAdaptor_Curve(edge)
        a, b = curve.FirstParameter(), curve.LastParameter()
        n = max(2, int(np.ceil(GCPnts_AbscissaPoint.Length_s(curve) / spacing)) + 1)
        for t in np.linspace(a, b, n):
            p = curve.Value(float(t))
            out.append((p.X(), p.Y(), p.Z()))
    return np.array(out, dtype=float).reshape(-1, 3)


def _triangles(shape, deflection: float):
    """``(verts, tris, face_index)`` of ``shape`` meshed on a copy, wound outward."""
    copy = BRepBuilderAPI_Copy(shape).Shape()
    occ.mesh_shape(copy, deflection)
    verts, tris, owner, offset = [], [], [], 0
    for fi, face in enumerate(occ.faces(copy)):
        got = occ.face_triangulation(face)
        if got is None:
            continue
        v, t = got
        if face.Orientation() == TopAbs_Orientation.TopAbs_REVERSED:
            t = t[:, ::-1]
        verts.append(v)
        tris.append(t + offset)
        owner.append(np.full(len(t), fi, dtype=np.int64))
        offset += len(v)
    if not verts:
        return (np.empty((0, 3)), np.empty((0, 3), dtype=np.int64),
                np.empty(0, dtype=np.int64))
    return np.vstack(verts), np.vstack(tris), np.concatenate(owner)


def from_shapes(body, lattice: list, deflection: float) -> Obstacles:
    """Obstacles for ``body`` (shape 0) and each lattice solid (shape ``k + 1``).

    The body is meshed on a copy so the classification triangulation it carries
    -- which the overhang regions are refined from -- is left as it was. Each
    triangle's anchor key numbers its (shape, anchor group) pair.
    """
    verts, tris, shape_id, key = [], [], [], []
    offset, key_base = 0, 0
    for k, shape in enumerate([body, *lattice]):
        v, t, face = _triangles(shape, deflection)
        groups = _face_groups(shape, material_inside=(k != BODY))
        verts.append(v)
        tris.append(t + offset)
        shape_id.append(np.full(len(t), k, dtype=np.int64))
        key.append(key_base + groups[face])
        offset += len(v)
        key_base += len(groups)
    obs = Obstacles(TriMesh(np.vstack(verts), np.vstack(tris)),
                    np.concatenate(shape_id), np.concatenate(key))
    from .cells import FOLD_EDGE_DEG, FOLD_SAMPLE

    obs.feature_points = convex_edge_points(body, False, FOLD_EDGE_DEG, FOLD_SAMPLE)
    return obs
