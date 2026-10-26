"""Overhang regions: the surface the support generator works on (docs/algorithm.md §14.3 S1).

The detector (:mod:`.overhang`) says *where* the body needs support, as points.
The generator needs more: a mesh of that surface fine enough to send a ray from
every ``h`` of it, with neighbours known, so that cells and the partings between
them can be read off it. This module builds that mesh.

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
from OCP.BRepAdaptor import BRepAdaptor_Curve
from OCP.GeomAbs import GeomAbs_SurfaceType
from OCP.ShapeAnalysis import ShapeAnalysis_Surface
from OCP.gp import gp_Pnt
from OCP.TopAbs import TopAbs_ShapeEnum
from OCP.TopLoc import TopLoc_Location
from OCP.TopoDS import TopoDS

from .. import occ
from ..errors import ProcessingError
from .overhang import PREFILTER_SLACK, _SurfaceNormals
from .params import REQUIRED_SLACK, SIN_START, STEP

WELD_DECIMALS = 8
"""Decimal places vertices are rounded to when welding faces together. Shared
boundary points are computed from one edge curve by both faces, so they agree
to the last bits; this only has to absorb that."""

ISOCLINE_TOL = 1e-4
"""mm. How exactly the 60 deg isocline is located on the true surface."""


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


@dataclass
class Region:
    """One enclosed overhang area: a connected set of required triangles."""

    tris: np.ndarray
    """Indices into :attr:`RegionMesh.tris`."""
    verts: np.ndarray
    """Indices into :attr:`RegionMesh.verts` of the triangles' corners."""
    area: float
    loops: list
    """Boundary loops, each an ``(n,)`` array of vertex indices, closed (the
    first vertex is not repeated)."""


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


def _flip_to_delaunay(P: np.ndarray, tris: np.ndarray) -> np.ndarray:
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
    """
    tris = tris.copy()
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


def _refine_face(face, step: float, active=None):
    """``(points, uvs, tris)`` of ``face``'s triangulation, refined to ``step``.

    ``active(uv)`` -- ``(n, 2) -> (n,) bool`` -- says which nodes could be over
    the bar; only triangles with an active corner are refined (their neighbours
    only as far as conformity needs), so a large face with a small overhang is
    not refined all over. ``None`` refines everything.
    """
    loc = TopLoc_Location()
    tri = BRep_Tool.Triangulation_s(face, loc)
    got = occ.face_uv_triangulation(face)
    if tri is None or got is None:
        return None
    pts, uvs, tris, evaluate = got
    tris = _flip_to_delaunay(pts, tris)
    pts, uvs = list(pts), list(uvs)
    bmap = _boundary_polygons(face, tri, loc)
    trsf = loc.Transformation()
    projector = None
    node_active = (np.ones(len(pts), dtype=bool) if active is None
                   else np.asarray(active(np.asarray(uvs)), dtype=bool))
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
        too_long = (length[edge_of[np.arange(m), longest]] > step) & node_active[tris].any(axis=1)
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
        tris = _split_patterns(tris, mid_of[edge_of], longest)
        tris = _flip_to_delaunay(np.asarray(pts), tris)
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


def build(body, b: np.ndarray, step: float = STEP, report=None) -> RegionMesh:
    """The refined, welded mesh of every face of ``body`` that may need support.

    ``body`` must carry the classification triangulation
    (:func:`latticegen2.classify.tessellate_surface`).
    """
    b = np.asarray(b, dtype=float)
    bar = SIN_START + REQUIRED_SLACK
    all_p, all_n, all_tris, all_face, all_corner_f, all_uv = [], [], [], [], [], []
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
        refined = _refine_face(face, step, active)
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
    keys = np.round(P, WELD_DECIMALS)
    _u, first, inverse = np.unique(keys, axis=0, return_index=True, return_inverse=True)
    inverse = inverse.reshape(-1)
    verts = P[first]
    normals = np.zeros_like(verts)
    np.add.at(normals, inverse, N)
    normals /= np.linalg.norm(normals, axis=1, keepdims=True)
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


def _edge_table(tris: np.ndarray):
    """``(edges, owner_tris)``: each triangle edge, sorted, with its triangle."""
    e = np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]])
    owner = np.tile(np.arange(len(tris)), 3)
    return np.sort(e, axis=1), owner


def _components(n: int, pairs: np.ndarray) -> np.ndarray:
    """Union-find labels for ``n`` items joined by ``pairs``."""
    parent = np.arange(n)

    def find(i):
        root = i
        while parent[root] != root:
            root = parent[root]
        while parent[i] != root:
            parent[i], i = root, parent[i]
        return root

    for a, c in pairs:
        ra, rc = find(int(a)), find(int(c))
        if ra != rc:
            parent[max(ra, rc)] = min(ra, rc)
    return np.array([find(i) for i in range(n)])


def _loops(boundary: np.ndarray) -> list:
    """Chain boundary edges into closed loops, or ``None`` if one does not close.

    A vertex where the boundary pinches (degree 4) is left by the edge that turns
    least, which is enough to keep loops apart for the closed-loop check; an odd
    degree anywhere means an open boundary.
    """
    from collections import defaultdict

    adj = defaultdict(list)
    for a, c in boundary:
        adj[int(a)].append(int(c))
        adj[int(c)].append(int(a))
    if any(len(v) % 2 for v in adj.values()):
        return None
    unused = {(min(a, c), max(a, c)): 0 for a, c in boundary}
    for a, c in boundary:
        unused[(min(a, c), max(a, c))] += 1
    loops = []
    for start in list(adj):
        while any(unused[(min(start, o), max(start, o))] for o in adj[start]):
            loop = [start]
            prev, cur = None, start
            while True:
                nxt = None
                for o in adj[cur]:
                    key = (min(cur, o), max(cur, o))
                    if unused[key] and o != prev:
                        nxt = o
                        break
                if nxt is None:
                    for o in adj[cur]:
                        key = (min(cur, o), max(cur, o))
                        if unused[key]:
                            nxt = o
                            break
                if nxt is None:
                    return None
                unused[(min(cur, nxt), max(cur, nxt))] -= 1
                if nxt == start:
                    break
                loop.append(nxt)
                prev, cur = cur, nxt
            loops.append(np.asarray(loop, dtype=np.int64))
    return loops


def regions(mesh: RegionMesh) -> list[Region]:
    """The enclosed overhang areas of ``mesh``, largest first (gate G1 inside).

    Fails (exit 4) naming the region if a boundary does not close into loops.
    """
    req = np.nonzero(mesh.required)[0]
    if len(req) == 0:
        return []
    sub = mesh.tris[req]
    edges, owner = _edge_table(sub)
    order = np.lexsort((edges[:, 1], edges[:, 0]))
    edges, owner = edges[order], owner[order]
    same = np.all(edges[1:] == edges[:-1], axis=1)
    pairs = np.column_stack([owner[:-1][same], owner[1:][same]])
    labels = _components(len(sub), pairs)

    a, c, d = (mesh.verts[sub[:, k]] for k in range(3))
    tri_area = 0.5 * np.linalg.norm(np.cross(c - a, d - a), axis=1)
    out = []
    for lab in np.unique(labels):
        local = np.nonzero(labels == lab)[0]
        e, _o = _edge_table(sub[local])
        uniq, counts = np.unique(e, axis=0, return_counts=True)
        boundary = uniq[counts == 1]
        loops = _loops(boundary)
        if loops is None:
            centre = mesh.verts[np.unique(sub[local])].mean(axis=0)
            raise ProcessingError(
                f"An overhang region near {np.round(centre, 3).tolist()} has a "
                f"boundary that does not close into loops; it cannot be divided "
                f"into support cells (docs/algorithm.md §14.3 S1)."
            )
        out.append(Region(tris=req[local], verts=np.unique(sub[local]),
                          area=float(tri_area[local].sum()), loops=loops))
    out.sort(key=lambda r: -r.area)
    return out
