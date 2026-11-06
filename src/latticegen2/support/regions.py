"""Ceiling points: the surface the support generator works on (docs/algorithm.md §14.3).

The detector (:mod:`.overhang`) says *where* the body needs support, as points.
The generator needs them denser and on a mesh: a ceiling point every ``h`` of
that surface, each with its exact normal, so that a support backs all of it.
This module builds that mesh.

Every face that could reach the 60 deg bar is taken from the body's
classification triangulation and refined until no edge is longer than ``h``:

* **on the true surface** -- an interior split point is evaluated at the UV
  midpoint, through the face's own surface;
* **exactly on the edge curve along a face's boundary** -- a split point there
  is evaluated at the edge-parameter midpoint, read from the polygon BRepMesh
  stored for that edge. Both faces sharing the edge compute the same point from
  the same curve, so the refined faces weld into one mesh without a crack;
* **conformingly, by longest-edge bisection** -- a triangle is only ever cut
  across its longest edge, and an edge is split for both triangles using it, so
  the mesh never acquires a T-junction;
* **into well-shaped triangles** -- interior edges are flipped to the Delaunay
  criterion before and between bisection passes. Bisection alone keeps the
  shape class it starts from, and BRepMesh triangulates a large plane into
  slivers: `Hole-block-basic`'s roof starts at a median quality of 0.03 and
  bisected without flips into 172,196 triangles where 40,000 do.

Each vertex carries its exact outward normal, and the 60 deg **isocline** -- where
``n . b`` crosses the bar inside a curved face -- is located by bisection on the
true surface rather than read off the mesh staircase.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
from OCP.GeomAbs import GeomAbs_SurfaceType
from OCP.ShapeAnalysis import ShapeAnalysis_Surface
from OCP.gp import gp_Pnt, gp_Vec
from OCP.TopAbs import TopAbs_ShapeEnum
from OCP.TopLoc import TopLoc_Location
from OCP.TopoDS import TopoDS

from .. import occ
from .overhang import PREFILTER_SLACK, _SurfaceNormals
from .params import REQUIRED_SLACK, SIN_START, STEP

WELD_DECIMALS = 8
"""Decimal places vertices are rounded to when welding faces together. Shared
boundary points are computed from one edge curve by both faces, so they agree
to the last bits; this only has to absorb that."""

ISOCLINE_TOL = 1e-4
"""mm. How exactly the 60 deg isocline is located on the true surface."""


COVER_STEP = STEP / 4.0
"""mm. Edge length the region mesh is refined to across the outline of a strut
end meeting the ceiling: within the crease push of 0.2 mm, so a border vertex
pushed into the strut lands in it."""


@dataclass
class RegionMesh:
    """The refined surface that may need support, welded across faces."""

    verts: np.ndarray
    """``(V, 3)`` points on the true surface."""
    normals: np.ndarray
    """``(V, 3)`` unit outward normals (averaged where a vertex lies on a crease
    between two faces)."""
    tris: np.ndarray
    """``(T, 3)`` vertex indices, wound with the outward normal."""
    tri_face: np.ndarray
    """``(T,)`` index of the body face each triangle lies on."""
    required: np.ndarray
    """``(T,)`` bool: every corner of the triangle overhangs more than 60 deg,
    measured with its own face's normal there."""
    vert_required: np.ndarray
    """``(V,)`` bool: the vertex is a corner of a required triangle."""
    isocline: np.ndarray = field(default_factory=lambda: np.empty((0, 3)))
    """``(K, 3)`` points of the 60 deg isocline inside curved faces, on the true
    surface to :data:`ISOCLINE_TOL`."""
    isocline_edges: np.ndarray = field(
        default_factory=lambda: np.empty((0, 2), dtype=np.int64))
    """``(K, 2)`` the mesh edge each isocline point was found on."""
    step: float = STEP
    tri_uv: np.ndarray = field(default_factory=lambda: np.empty((0, 3, 2)))
    """``(T, 3, 2)`` each triangle corner's parameters on the triangle's own face."""
    faces: list = field(default_factory=list)
    """The body's faces, indexed by :attr:`tri_face`."""

    def __len__(self) -> int:
        return len(self.verts)

    def surface_points(self, tri: np.ndarray, bary: np.ndarray):
        """``(points, normals)`` on the true surface at barycentric ``bary`` of ``tri``.

        Interpolated in the triangle's own face parameters and evaluated there,
        so a point between two mesh vertices lies on the surface, not on the chord.
        """
        tri = np.asarray(tri, dtype=np.int64)
        bary = np.atleast_2d(bary)
        uv = np.einsum("ik,ikj->ij", bary, self.tri_uv[tri])
        pts = np.empty((len(tri), 3))
        nrm = np.empty((len(tri), 3))
        cache = {}
        for k, (t, (u, v)) in enumerate(zip(tri, uv)):
            fi = int(self.tri_face[t])
            sn = cache.get(fi)
            if sn is None:
                sn = cache[fi] = _SurfaceNormals(self.faces[fi])
            p, n = sn.at(u, v)
            pts[k] = p
            nrm[k] = n if n is not None else self.normals[self.tris[t]].mean(axis=0)
        return pts, nrm


# --- per face ----------------------------------------------------------------


def _could_reach_bar(face, b: np.ndarray, uvs: np.ndarray, normals: _SurfaceNormals,
                     bar: float) -> bool:
    if normals.surf.GetType() == GeomAbs_SurfaceType.GeomAbs_Plane:
        _p, n = normals.at(*uvs[0])
        return n is not None and float(n @ b) > bar
    best = -2.0
    for u, v in uvs:
        _p, n = normals.at(u, v)
        best = max(best, 2.0 if n is None else float(n @ b))
    return best >= bar - PREFILTER_SLACK


def _boundary_polygons(face, tri, loc):
    """``{(i, j): (curve, t_i, t_j, pcurve)}`` for every mesh edge on a B-rep edge.

    Seam edges are left out: their pcurves are isolines, so the UV midpoint the
    interior rule uses is already on the seam, and the two copies of the seam
    in one face would otherwise collide in this map.
    """
    out = {}
    for e in occ._explore(face, TopAbs_ShapeEnum.TopAbs_EDGE):
        edge = TopoDS.Edge_s(e)
        if BRep_Tool.Degenerated_s(edge) or BRep_Tool.IsClosed_s(edge, face):
            continue
        poly = BRep_Tool.PolygonOnTriangulation_s(edge, tri, loc)
        if poly is None or not poly.HasParameters():
            continue
        pcurve = BRep_Tool.CurveOnSurface_s(edge, face, 0.0, 0.0)
        if pcurve is None:
            continue
        curve = BRepAdaptor_Curve(edge)
        nodes = [poly.Node(k) - 1 for k in range(1, poly.NbNodes() + 1)]
        params = [poly.Parameter(k) for k in range(1, poly.NbNodes() + 1)]
        for k in range(len(nodes) - 1):
            i, j = nodes[k], nodes[k + 1]
            ti, tj = params[k], params[k + 1]
            if i > j:
                i, j, ti, tj = j, i, tj, ti
            out[(i, j)] = (curve, ti, tj, pcurve)
    return out


FLIP_PASSES_MAX = 400
"""Rounds of simultaneous independent flips before giving up on a face. A round
flips every improving edge whose two triangles no better candidate claims, so
convergence takes tens of rounds, not hundreds."""

STRAY_FRACTION = 0.25
"""A split point further than this fraction of its edge from the edge's own
midpoint is taken from a projection instead of the UV midpoint."""

REFINE_PASSES_MAX = 64
"""Bisection passes per face. Each halves the longest edges, so 64 would refine
a kilometre to a micron; the cap is a backstop against a degenerate face, not
a limit any real face reaches."""


def _min_angle(A, B, C):
    def angle(o, p, q):
        u, v = p - o, q - o
        nu, nv = np.linalg.norm(u, axis=1), np.linalg.norm(v, axis=1)
        cosv = np.einsum("ij,ij->i", u, v) / np.maximum(nu * nv, 1e-300)
        return np.where((nu > 0) & (nv > 0), np.arccos(np.clip(cosv, -1.0, 1.0)), 0.0)

    return np.minimum(np.minimum(angle(A, B, C), angle(B, C, A)), angle(C, A, B))


def _uv_turn(uv: np.ndarray, a, b, c) -> np.ndarray:
    """Signed doubled area of UV triangles ``(a, b, c)``."""
    return ((uv[b, 0] - uv[a, 0]) * (uv[c, 1] - uv[a, 1])
            - (uv[b, 1] - uv[a, 1]) * (uv[c, 0] - uv[a, 0]))


FLIP_FIDELITY = 1.5
"""Most a flipped triangle's area may exceed the area of surface it stands for
(its parameter triangle times the surface's area scale there). A faithful
triangle on a curved face is close to 1; a chord bridging a boundary curve that
is straight in the parameters stands for almost nothing."""


def _area_scale(surf, uv: np.ndarray) -> np.ndarray:
    """``|S_u x S_v|`` of ``surf`` at each parameter pair: surface area per unit
    parameter area."""
    p, du, dv = gp_Pnt(), gp_Vec(), gp_Vec()
    out = np.empty(len(uv))
    for k, (u, v) in enumerate(uv):
        surf.D1(float(u), float(v), p, du, dv)
        out[k] = du.Crossed(dv).Magnitude()
    return out


def _flip_to_delaunay(P: np.ndarray, tris: np.ndarray, uv=None, scale=None) -> np.ndarray:
    """Flip interior edges while a flip raises the pair's smallest angle.

    An edge ``(a, b)`` between triangles ``(a, b, c)`` and ``(b, a, d)`` is
    flipped to ``(c, d)`` when the two new triangles' smallest angle beats the
    old pair's -- the max-min-angle criterion, which on a plane is the Delaunay
    one -- and only when both new triangles keep the winding of the old, so a
    flip can never fold the mesh. Edges used once (the face's boundary, seams)
    are never touched. Angles are measured in 3D on vertices that lie on the
    true surface.

    The criterion is max-min rather than the opposite-angle sum because the
    latter is ill-posed on a triangle with no extent: at a sphere's pole a face's
    triangulation has corners that coincide in 3D, and the angle-sum test flipped
    edges there back and forth against the bisection for as long as it ran. A
    degenerate triangle's smallest angle is zero, and a flip that cannot better
    zero is not made.

    **With ``uv``, a flip must also keep both new triangles wound as the old
    pair in the face's parameters**, which is the test that cannot be fooled:
    the parameter plane is flat. The 3D test compares against the sum of the
    pair's normals, and on long thin triangles each spanning a wide arc of a
    curved face that sum points nowhere useful. Measured on a half-cylinder of
    radius 7.5 mm and 79 mm length (test-cylinder.STEP's hole): the first round
    of flips folded 32 of the classification mesh's 96 triangles, and
    refinement never ended -- 680,000 triangles on 1,850 mm^2, and growing.

    **With ``scale`` too, a flip must keep each new triangle faithful to the
    surface** (:data:`FLIP_FIDELITY`). Measured on the same face: its end
    curves are arcs in 3D and almost straight in the parameters, and the 3D
    criterion flipped 80 of 96 triangles into fans along them -- well shaped in
    3D, slivers in the parameters, standing for no surface. Splitting such a
    triangle's chord at its parameter midpoint lands on the boundary point
    already there, and the duplicate breeds flat triangles that never refine
    away: 135,000 triangles, half of them flat. On a plane the two areas are
    equal and the test never refuses; at a pole it can only refuse a flip.
    """
    tris = tris.copy()
    if uv is not None:
        uv = np.asarray(uv, dtype=float)
    for _ in range(FLIP_PASSES_MAX):
        m = len(tris)
        a = np.concatenate([tris[:, 0], tris[:, 1], tris[:, 2]])
        b = np.concatenate([tris[:, 1], tris[:, 2], tris[:, 0]])
        c = np.concatenate([tris[:, 2], tris[:, 0], tris[:, 1]])
        owner = np.tile(np.arange(m), 3)
        key = np.minimum(a, b) * (len(P) + 1) + np.maximum(a, b)
        order = np.argsort(key, kind="stable")
        ks = key[order]
        pair = np.nonzero(ks[1:] == ks[:-1])[0]
        if len(pair) == 0:
            return tris
        i1, i2 = order[pair], order[pair + 1]
        A, B, C, D = P[a[i1]], P[b[i1]], P[c[i1]], P[c[i2]]

        before = np.minimum(_min_angle(A, B, C), _min_angle(B, A, D))
        after = np.minimum(_min_angle(C, A, D), _min_angle(C, D, B))
        excess = after - before
        nrm = np.cross(B - A, C - A) + np.cross(A - B, D - B)
        new1 = np.cross(A - C, D - C)
        new2 = np.cross(D - C, B - C)
        ok = ((excess > 1e-9)
              & (np.einsum("ij,ij->i", new1, nrm) > 0.0)
              & (np.einsum("ij,ij->i", new2, nrm) > 0.0))
        if uv is not None:
            side = np.sign(_uv_turn(uv, a[i1], b[i1], c[i1]))
            turn1 = _uv_turn(uv, c[i1], a[i1], c[i2])
            turn2 = _uv_turn(uv, c[i1], c[i2], b[i1])
            ok &= (np.sign(turn1) == side) & (np.sign(turn2) == side) & (side != 0)
            if scale is not None:
                s1 = np.abs(turn1) * (scale[c[i1]] + scale[a[i1]] + scale[c[i2]]) / 3.0
                s2 = np.abs(turn2) * (scale[c[i1]] + scale[c[i2]] + scale[b[i1]]) / 3.0
                ok &= ((np.linalg.norm(new1, axis=1) <= FLIP_FIDELITY * s1)
                       & (np.linalg.norm(new2, axis=1) <= FLIP_FIDELITY * s2))
        cand = np.nonzero(ok)[0]
        if len(cand) == 0:
            return tris
        # Independent set: each triangle goes to the candidate violating most.
        cand = cand[np.argsort(-excess[cand], kind="stable")]
        rank = np.empty(len(cand), dtype=np.int64)
        rank[:] = np.arange(len(cand))
        claim = np.full(m, len(cand), dtype=np.int64)
        np.minimum.at(claim, owner[i1[cand]], rank)
        np.minimum.at(claim, owner[i2[cand]], rank)
        win = (claim[owner[i1[cand]]] == rank) & (claim[owner[i2[cand]]] == rank)
        cand = cand[win]
        t1, t2 = owner[i1[cand]], owner[i2[cand]]
        av, bv, cv, dv = a[i1[cand]], b[i1[cand]], c[i1[cand]], c[i2[cand]]
        tris[t1] = np.column_stack([cv, av, dv])
        tris[t2] = np.column_stack([cv, dv, bv])
    return tris


def _mark_longest(edge_of: np.ndarray, longest: np.ndarray, marked: np.ndarray) -> None:
    """Close ``marked`` under "a triangle with a marked edge marks its longest".

    That rule is what makes every split a bisection of the longest edge, and it
    terminates: each round only adds marks, and there are finitely many edges.
    """
    m = len(edge_of)
    rows = np.arange(m)
    while True:
        touched = marked[edge_of].any(axis=1)
        need = touched & ~marked[edge_of[rows, longest]]
        if not need.any():
            return
        marked[edge_of[rows[need], longest[need]]] = True


def _split_patterns(tris: np.ndarray, mid: np.ndarray, longest: np.ndarray) -> np.ndarray:
    """New triangles after bisecting, longest edge first.

    ``mid[t, k]`` is the midpoint vertex of triangle ``t``'s edge ``k`` (from
    corner ``k`` to corner ``k + 1``), or -1; whenever any edge of ``t`` is split
    its longest edge, ``longest[t]``, is too. The triangle is rotated so that edge
    is edge 0 and cut from its midpoint to the opposite corner; a split of either
    remaining edge then cuts the half it belongs to from that midpoint. Winding
    is kept.
    """
    count = (mid >= 0).sum(axis=1)
    out = [tris[count == 0]]
    for k in range(3):
        sel = (count > 0) & (longest == k)
        if not sel.any():
            continue
        t = np.roll(tris[sel], -k, axis=1)
        mm = np.roll(mid[sel], -k, axis=1)
        v0, v1, v2 = t[:, 0], t[:, 1], t[:, 2]
        m0, m1, m2 = mm[:, 0], mm[:, 1], mm[:, 2]
        # The half (m0, v1, v2) holds edge 1; the half (v0, m0, v2) holds edge 2.
        a = m1 < 0
        out.append(np.column_stack([m0[a], v1[a], v2[a]]))
        out.append(np.column_stack([m0[~a], v1[~a], m1[~a]]))
        out.append(np.column_stack([m0[~a], m1[~a], v2[~a]]))
        c = m2 < 0
        out.append(np.column_stack([v0[c], m0[c], v2[c]]))
        out.append(np.column_stack([v0[~c], m0[~c], m2[~c]]))
        out.append(np.column_stack([m0[~c], v2[~c], m2[~c]]))
    return np.vstack(out)


def _refine_face(face, step: float, active=None, covered=None, fine_step=None):
    """``(points, uvs, tris)`` of ``face``'s triangulation, refined to ``step``.

    ``active(uv)`` -- ``(n, 2) -> (n,) bool`` -- says which nodes could be over
    the bar; only triangles with an active corner are refined (their neighbours
    only as far as conformity needs), so a large face with a small overhang is
    not refined all over. ``None`` refines everything.

    ``covered(points)`` -- ``(n, 3) -> (n,) bool`` -- says which nodes lie under
    lattice material; a triangle whose corners disagree straddles the outline
    of a strut end and is refined further, to ``fine_step``.
    """
    loc = TopLoc_Location()
    tri = BRep_Tool.Triangulation_s(face, loc)
    got = occ.face_uv_triangulation(face)
    if tri is None or got is None:
        return None
    pts, uvs, tris, evaluate = got
    surf = BRepAdaptor_Surface(face, False)
    scale = list(_area_scale(surf, uvs))
    tris = _flip_to_delaunay(pts, tris, uvs, np.asarray(scale))
    pts, uvs = list(pts), list(uvs)
    bmap = _boundary_polygons(face, tri, loc)
    trsf = loc.Transformation()
    projector = None
    node_active = (np.ones(len(pts), dtype=bool) if active is None
                   else np.asarray(active(np.asarray(uvs)), dtype=bool))
    node_covered = (None if covered is None
                    else np.asarray(covered(np.asarray(pts)), dtype=bool))
    for _pass in range(REFINE_PASSES_MAX):
        P = np.asarray(pts)
        m = len(tris)
        e = np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]])
        e.sort(axis=1)
        uniq, inverse = np.unique(e, axis=0, return_inverse=True)
        edge_of = inverse.reshape(-1).reshape(3, m).T
        length = np.linalg.norm(P[uniq[:, 0]] - P[uniq[:, 1]], axis=1)
        longest = np.argmax(length[edge_of], axis=1)
        split = np.zeros(len(uniq), dtype=bool)
        worst = length[edge_of[np.arange(m), longest]]
        too_long = (worst > step) & node_active[tris].any(axis=1)
        if node_covered is not None:
            straddle = node_covered[tris].any(axis=1) & ~node_covered[tris].all(axis=1)
            too_long |= (worst > fine_step) & straddle & node_active[tris].any(axis=1)
        if not too_long.any():
            break
        split[edge_of[np.nonzero(too_long)[0], longest[too_long]]] = True
        _mark_longest(edge_of, longest, split)
        mid_of = np.full(len(uniq), -1, dtype=np.int64)
        interior_uv = []
        interior_slot = []
        interior_edge = []
        for u_idx in np.nonzero(split)[0]:
            i, j = int(uniq[u_idx, 0]), int(uniq[u_idx, 1])
            hit = bmap.get((i, j))
            new = len(pts)
            mid_of[u_idx] = new
            if hit is None:
                pts.append(None)
                uvs.append(0.5 * (np.asarray(uvs[i]) + np.asarray(uvs[j])))
                interior_uv.append(uvs[-1])
                interior_slot.append(new)
                interior_edge.append(u_idx)
                continue
            curve, ti, tj, pcurve = hit
            tm = 0.5 * (ti + tj)
            p = curve.Value(tm).Transformed(trsf)
            uv = pcurve.Value(tm)
            pts.append(np.array([p.X(), p.Y(), p.Z()]))
            uvs.append(np.array([uv.X(), uv.Y()]))
            bmap[(i, new)] = (curve, ti, tm, pcurve) if i < new else (curve, tm, ti, pcurve)
            bmap[(j, new)] = (curve, tj, tm, pcurve) if j < new else (curve, tm, tj, pcurve)
        if interior_uv:
            values = evaluate(np.asarray(interior_uv))
            for slot, value, u_idx in zip(interior_slot, values, interior_edge):
                i, j = int(uniq[u_idx, 0]), int(uniq[u_idx, 1])
                chord = 0.5 * (P[i] + P[j])
                if np.linalg.norm(value - chord) > STRAY_FRACTION * length[u_idx]:
                    # A UV midpoint that lands away from its edge: at a pole,
                    # where one node's parameters are arbitrary along the pole
                    # line, the UV midpoint of an edge to it lies at an
                    # unrelated azimuth, the split never shortens that edge,
                    # and refinement does not end. Project the 3D midpoint.
                    if projector is None:
                        projector = ShapeAnalysis_Surface(BRep_Tool.Surface_s(face))
                    loc_p = gp_Pnt(*chord).Transformed(trsf.Inverted())
                    uv2 = projector.ValueOfUV(loc_p, 1e-9)
                    uvs[slot] = np.array([uv2.X(), uv2.Y()])
                    value = evaluate(uvs[slot][None])[0]
                pts[slot] = value
        added = len(pts) - len(node_active)
        if added:
            fresh = (np.ones(added, dtype=bool) if active is None
                     else np.asarray(active(np.asarray(uvs[-added:])), dtype=bool))
            node_active = np.concatenate([node_active, fresh])
            if node_covered is not None:
                node_covered = np.concatenate([node_covered, np.asarray(
                    covered(np.asarray(pts[-added:])), dtype=bool)])
        tris = _split_patterns(tris, mid_of[edge_of], longest)
        scale += list(_area_scale(surf, np.asarray(uvs[len(scale):])))
        tris = _flip_to_delaunay(np.asarray(pts), tris, np.asarray(uvs), np.asarray(scale))
    return np.asarray(pts), np.asarray(uvs), tris


def _bisect_isocline(normals: _SurfaceNormals, b, uv_a, uv_b, f_a, bar, length):
    """UV point where ``n . b`` crosses ``bar`` between two parameters."""
    lo, hi = np.asarray(uv_a, float), np.asarray(uv_b, float)
    # Halve until the 3D segment is shorter than the tolerance.
    for _ in range(max(1, int(np.ceil(np.log2(max(length, ISOCLINE_TOL) / ISOCLINE_TOL))) + 1)):
        mid = 0.5 * (lo + hi)
        _p, n = normals.at(*mid)
        f = -1.0 if n is None else float(n @ b) - bar
        if (f > 0) == (f_a > 0):
            lo = mid
        else:
            hi = mid
    mid = 0.5 * (lo + hi)
    p, _n = normals.at(*mid)
    return p


# --- the whole body ------------------------------------------------------------


def build(body, b: np.ndarray, step: float = STEP, report=None,
          covered=None) -> RegionMesh:
    """The refined, welded mesh of every face of ``body`` that may need support.

    ``body`` must carry the classification triangulation
    (:func:`latticegen2.classify.tessellate_surface`).

    ``covered(points)`` says which surface points lie under lattice material:
    where a strut meets the ceiling, the lattice backs it already
    (specification.md §4.6), so those points are not required and the strut's
    end is a hole in its region. Its outline is not an edge of the body, which
    is how a wall's crease is found exactly, so it is found by refining: the
    triangles straddling it go down to :data:`COVER_STEP`, putting the hole's
    border within that of the outline, where the crease push (§14.3 S3d) carries
    the support's edge on into the strut.
    """
    b = np.asarray(b, dtype=float)
    bar = SIN_START + REQUIRED_SLACK
    all_p, all_n, all_tris, all_face, all_corner_f, all_uv = [], [], [], [], [], []
    all_point_face = []
    iso_p, iso_e = [], []
    offset = 0
    faces = occ.faces(body)
    for fi, face in enumerate(faces):
        if report is not None:
            report("refining overhang faces", fi, len(faces))
        got = occ.face_uv_triangulation(face)
        if got is None:
            continue
        normals = _SurfaceNormals(face)
        if not _could_reach_bar(face, b, got[1], normals, bar):
            continue
        if normals.surf.GetType() == GeomAbs_SurfaceType.GeomAbs_Plane:
            active = None
        else:
            def active(uv, normals=normals):
                out = np.empty(len(uv), dtype=bool)
                for k, (u, v) in enumerate(uv):
                    _p, n = normals.at(u, v)
                    out[k] = n is None or float(n @ b) >= bar - PREFILTER_SLACK
                return out
        refined = _refine_face(face, step, active, covered, COVER_STEP)
        if refined is None:
            continue
        pts, uvs, tris = refined
        nrm = np.empty_like(pts)
        f = np.empty(len(pts))
        for k, (u, v) in enumerate(uvs):
            _p, n = normals.at(u, v)
            if n is None:
                # A pole: take the normal a hair inside the face's parameter
                # range from the nearest triangle's centroid.
                owner = tris[np.nonzero((tris == k).any(axis=1))[0][0]]
                cu, cv = uvs[owner].mean(axis=0)
                _p, n = normals.at(cu, cv)
            nrm[k] = n
            f[k] = float(n @ b) - bar
        if covered is not None:
            # Under a strut end: not required, whatever its slope.
            f = np.where(np.asarray(covered(pts), dtype=bool), -1.0, f)
        # Isocline points on edges whose ends straddle the bar.
        e = np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]])
        e = np.unique(np.sort(e, axis=1), axis=0)
        cross = (f[e[:, 0]] > 0) != (f[e[:, 1]] > 0)
        for i, j in e[cross]:
            length = float(np.linalg.norm(pts[i] - pts[j]))
            iso_p.append(_bisect_isocline(normals, b, uvs[i], uvs[j], f[i], bar, length))
            iso_e.append((offset + i, offset + j))
        all_p.append(pts)
        all_n.append(nrm)
        all_tris.append(tris + offset)
        all_face.append(np.full(len(tris), fi, dtype=np.int64))
        all_corner_f.append(f[tris])
        all_uv.append(uvs[tris])
        all_point_face.append(np.full(len(pts), fi, dtype=np.int64))
        offset += len(pts)
    if report is not None:
        report("refining overhang faces", len(faces), len(faces))
    if not all_tris:
        return RegionMesh(np.empty((0, 3)), np.empty((0, 3)),
                          np.empty((0, 3), dtype=np.int64), np.empty(0, dtype=np.int64),
                          np.empty(0, dtype=bool), np.empty(0, dtype=bool), step=step)

    P = np.vstack(all_p)
    N = np.vstack(all_n)
    tris = np.vstack(all_tris)
    tri_face = np.concatenate(all_face)
    corner_f = np.vstack(all_corner_f)
    tri_uv = np.vstack(all_uv)
    # Welded within a smooth face group only. Across a crease -- the edge
    # between Hole-block-chamfered's flat and 15 deg roof facets -- the mesh
    # stays split, and the copies of a crease point share one averaged normal,
    # so the point probed into the body from it lies between the two faces.
    group = _smooth_groups(body)
    point_face = np.concatenate(all_point_face)
    pos = np.round(P, WELD_DECIMALS)
    _up, pos_inv = np.unique(pos, axis=0, return_inverse=True)
    pos_inv = pos_inv.reshape(-1)
    pos_normal = np.zeros((pos_inv.max() + 1, 3))
    np.add.at(pos_normal, pos_inv, N)
    pos_normal /= np.linalg.norm(pos_normal, axis=1, keepdims=True)
    keys = np.column_stack([pos, group[point_face].astype(float)])
    _u, first, inverse = np.unique(keys, axis=0, return_index=True, return_inverse=True)
    inverse = inverse.reshape(-1)
    verts = P[first]
    normals = pos_normal[pos_inv[first]]
    tris = inverse[tris]
    keep = (tris[:, 0] != tris[:, 1]) & (tris[:, 1] != tris[:, 2]) & (tris[:, 2] != tris[:, 0])
    tris, tri_face, corner_f, tri_uv = tris[keep], tri_face[keep], corner_f[keep], tri_uv[keep]
    required = (corner_f > 0).all(axis=1)
    vert_required = np.zeros(len(verts), dtype=bool)
    vert_required[tris[required].ravel()] = True
    iso_edges = (inverse[np.asarray(iso_e, dtype=np.int64)] if iso_e
                 else np.empty((0, 2), dtype=np.int64))
    return RegionMesh(verts, normals, tris, tri_face, required, vert_required,
                      np.asarray(iso_p) if iso_p else np.empty((0, 3)),
                      iso_edges, step, tri_uv, faces)


SMOOTH_DEG = 1.0
"""Degrees. Faces meeting at less than this are one smooth surface to a region."""


def _smooth_groups(body) -> np.ndarray:
    """Group per face of ``body``: faces joined across edges where they meet smoothly."""
    from OCP.TopExp import TopExp
    from OCP.TopTools import (TopTools_IndexedDataMapOfShapeListOfShape,
                              TopTools_IndexedMapOfShape)

    faces = occ.faces(body)
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
    TopExp.MapShapesAndAncestors_s(body, TopAbs_ShapeEnum.TopAbs_EDGE,
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
        curve = BRepAdaptor_Curve(edge)
        t = 0.5 * (curve.FirstParameter() + curve.LastParameter())
        normals = []
        for f in adj:
            c2d = BRep_Tool.CurveOnSurface_s(edge, f, 0.0, 0.0)
            if c2d is None:
                break
            uv = c2d.Value(t)
            _p, n = _SurfaceNormals(f).at(uv.X(), uv.Y())
            if n is None:
                break
            normals.append(n)
        if len(normals) == 2 and np.degrees(np.arccos(
                np.clip(float(normals[0] @ normals[1]), -1.0, 1.0))) < SMOOTH_DEG:
            parent[max(find(ia), find(ib))] = min(find(ia), find(ib))
    return np.array([find(i) for i in range(len(faces))], dtype=np.int64)
