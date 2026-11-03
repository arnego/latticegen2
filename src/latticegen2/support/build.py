"""Booleans and gates for support solids (docs/algorithm.md §14.3).

A support body is trimmed to the input body with an intersection, and the
kernel can return a wrong answer while reporting success (§7.2). So the result
is never taken on the boolean's word: :class:`Gates` asks it questions the
boolean did not answer -- validity, material outside the body by the
classification mesh, and the overhang of every free face at its exact normals.
"""

from __future__ import annotations

import numpy as np
from OCP.TopTools import TopTools_ListOfShape

from .. import occ
from .params import SUPPORT_OVERHANG_DEG

FUZZY_LADDER = (0.0, 1e-5, 1e-4, 1e-3)
"""mm. Fuzzy values a trim is tried with, in order, until one passes the gates.
A fuzzy value moves nothing by more than itself, and the gates check the result
independently of the boolean that made it."""


def _boolean(op, a, tools: list, fuzzy: float = 0.0):
    algo = op()
    args = TopTools_ListOfShape()
    args.Append(a)
    tl = TopTools_ListOfShape()
    for s in tools:
        tl.Append(s)
    algo.SetArguments(args)
    algo.SetTools(tl)
    if fuzzy > 0.0:
        algo.SetFuzzyValue(fuzzy)
    algo.Build()
    if not algo.IsDone():
        return None
    return algo.Shape()


GATE_DEFLECTION = 0.02
"""mm. Tessellation the gates read a support solid through."""

ON_BODY_TOL = 0.02
"""mm, beyond the body mesh's own deviation: a face this close to the body's
surface lies on it."""

SELF_SUPPORT_TOL_DEG = 0.05
"""Degrees past the rod angle a free face may read before it fails: a planar
face on a lattice plane reads exactly the rod angle, to rounding."""


def _mesh_of(shape, deflection: float):
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Copy
    from OCP.TopAbs import TopAbs_Orientation

    copy = BRepBuilderAPI_Copy(shape).Shape()
    occ.mesh_shape(copy, deflection)
    per_face = []
    for face in occ.faces(copy):
        got = occ.face_uv_triangulation(face)
        if got is None:
            continue
        pts, uvs, tris, _ev = got
        if face.Orientation() == TopAbs_Orientation.TopAbs_REVERSED:
            tris = tris[:, ::-1]
        per_face.append((face, pts, uvs, tris))
    return per_face


class Gates:
    """The checks a support solid must pass, none of them the boolean that made it."""

    def __init__(self, body_mesh, b: np.ndarray):
        from ..classify import OutsideProbe, PointInside, SpatialHash

        self.body_mesh = body_mesh
        self.probe = OutsideProbe(body_mesh, 2.0 * body_mesh.deviation)
        self.hash = SpatialHash(body_mesh, min_cell=1.0)
        self.b = np.asarray(b, dtype=float)
        self.bar = float(np.sin(np.radians(SUPPORT_OVERHANG_DEG + SELF_SUPPORT_TOL_DEG)))
        self.on_tol = body_mesh.deviation + ON_BODY_TOL
        self.inside_body = PointInside(body_mesh)

    def _dist_to_body(self, p: np.ndarray) -> float:
        from ..classify import _point_triangle_dist

        r = self.on_tol + 1e-6
        cand = self.hash.query(p - r, p + r)
        if not len(cand):
            return np.inf
        a, b, c = self.body_mesh.triangle_points
        return float(_point_triangle_dist(p, a[cand], b[cand], c[cand]).min())

    def outside(self, solid) -> float:
        """How far ``solid``'s surface reaches outside the body, by the probe;
        ``0.0`` within its margin."""
        pts = np.vstack([p for _f, p, _u, _t in _mesh_of(solid, GATE_DEFLECTION)])
        return self.probe.worst_outside(pts)

    def check(self, solid) -> tuple[list, dict]:
        """``(failures, stats)`` of ``solid``: validity, material outside the
        body, and free faces overhanging past the rod angle."""
        from .overhang import _SurfaceNormals

        fails = []
        if not occ.is_valid(solid, parallel=False):
            fails.append("not BRepCheck-valid")
        per_face = _mesh_of(solid, GATE_DEFLECTION)
        allp = np.vstack([p for _f, p, _u, _t in per_face])
        worst = self.probe.worst_outside(allp)
        if worst > 0.0:
            fails.append(f"material {worst:.4f} mm outside the body")
        worst_deg = 0.0
        free_faces = 0
        for face, pts, uvs, _tris in per_face:
            # On the body's surface, the face is shared with a wall or ceiling.
            probe_pts = pts[:: max(1, len(pts) // 12)]
            if all(self._dist_to_body(p) <= self.on_tol for p in probe_pts):
                continue
            free_faces += 1
            sn = _SurfaceNormals(face)
            for u, v in uvs:
                _p, n = sn.at(u, v)
                if n is None:
                    continue
                down = -float(n @ self.b)
                if down > self.bar:
                    worst_deg = max(worst_deg, float(np.degrees(np.arcsin(min(1.0, down)))))
        if worst_deg > 0.0:
            fails.append(f"a free face overhangs {worst_deg:.2f} deg from vertical")
        return fails, {"free_faces": free_faces, "worst_free_face_deg": round(worst_deg, 3)}
