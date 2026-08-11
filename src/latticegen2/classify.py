"""Boundary classification — deciding what needs a boolean and what does not.

This is a port of docs/algorithm.md §5, with one structural change the fuse-free
architecture requires: the unit of classification is the **half-strut**, and node
classes are derived from the six half-struts a node owns.

The margin analysis is unchanged and is what keeps the whole design safe. A
half-strut is only promoted to INTERIOR/OUTSIDE when its axis segment is further
than ``r + d`` from the tessellated surface — ``r`` because the solid extends
that far from the axis, ``d`` because the mesh is only a chordal approximation
of the true B-rep. Anything closer degrades to BOUNDARY, which costs one small
intersection and can never produce a wrong answer (docs/algorithm.md §10).

Two consequences of that margin are load-bearing downstream and worth stating
explicitly, because the rest of the pipeline is built on them:

1. **One inside/outside test per node suffices.** If a half-strut's whole axis
   segment stays further than ``r + d`` from the surface, the segment cannot
   cross the surface, so its midpoint and its node are on the same side. The
   expensive ray-parity test therefore runs once per node, not once per
   half-strut.
2. **An INTERIOR node's neighbours are always kept.** If node ``A``'s half-strut
   toward ``B`` is INTERIOR, then ``B``'s half-strut back toward ``A`` covers the
   other half of a strut that never approaches the surface and whose far end is
   inside, so it cannot be OUTSIDE — meaning ``B`` is never classified OUTSIDE.
   That is why every cap of an INTERIOR node is an interface to real geometry
   (:mod:`latticegen2.interior`), and why an intact cap always has material on
   both sides.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum

import numpy as np
from OCP.TopoDS import TopoDS_Shape

from . import occ
from .errors import InputGeometryError
from .lattice import LatticeParams, half_strut_offset, nodes


class Klass(IntEnum):
    """Classification of a half-strut or node relative to the input solid."""

    INTERIOR = 0
    BOUNDARY = 1
    OUTSIDE = 2


@dataclass
class TriMesh:
    """A plain-data triangle surface mesh with welded (de-duplicated) vertices."""

    verts: np.ndarray
    """``(V, 3)`` float vertex positions."""
    tris: np.ndarray
    """``(T, 3)`` int vertex indices."""
    deviation: float = 0.0
    """Measured worst chordal deviation of this mesh from the true surface.

    This is the ``d`` of docs/algorithm.md §5, established by measurement rather
    than by trusting the mesher to have honoured the deflection it was asked
    for. See :func:`tessellate_surface`."""

    @property
    def triangle_points(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        v = self.verts
        return v[self.tris[:, 0]], v[self.tris[:, 1]], v[self.tris[:, 2]]

    def bounds(self) -> tuple[np.ndarray, np.ndarray]:
        return self.verts.min(axis=0), self.verts.max(axis=0)


def chordal_target(lp: LatticeParams) -> float:
    """The chordal deviation tolerance ``d`` (docs/algorithm.md §5.1)."""
    return min(lp.t, lp.a) / 10.0


# --- Tessellation and its completeness gate ---------------------------------


def _weld(verts: np.ndarray, tris: np.ndarray, decimals: int = 9):
    """Merge coincident vertices so shared edges get shared indices.

    OCCT triangulates each face independently, so a shared edge's nodes appear
    once per adjacent face at identical coordinates. Welding them is what makes
    the mesh usable for the manifold check in ``tools/verify_geometry.py``, and
    it shrinks the spatial hash.
    """
    keys = np.round(verts, decimals)
    _, first, inverse = np.unique(keys, axis=0, return_index=True, return_inverse=True)
    welded = verts[first]
    remapped = inverse[tris]
    # Drop triangles that welding collapsed to a degenerate sliver.
    keep = (
        (remapped[:, 0] != remapped[:, 1])
        & (remapped[:, 1] != remapped[:, 2])
        & (remapped[:, 2] != remapped[:, 0])
    )
    return welded, remapped[keep]


def check_surface_mesh_coverage(
    shape: TopoDS_Shape,
    lp: LatticeParams,
    area_rel_tol: float = 0.25,
    min_deficit_factor: float = 1.0,
    bbox_slack: float = 1.0,
) -> None:
    """Defensive completeness gate on the freshly generated surface mesh.

    Every classification decision depends on this mesh faithfully representing
    the input solid's boundary; a silently incomplete mesh would misclassify
    struts beyond the missing region as OUTSIDE, which is exactly the "wrong
    answer" failure mode the rest of the design rules out. This converts it into
    a loud one (:class:`InputGeometryError`, exit 3).

    Two per-face tests, each using a quantity that is sound *in the direction it
    is used*:

    1. **Coverage against exact trimmed area.** ``BRepGProp`` gives OCCT's exact
       Gauss-quadrature area of the trimmed face; the mesh must not fall short of
       it by both more than ``area_rel_tol`` of the face **and** more than one
       coarsest element's worth of area. Both bars are load-bearing: a relative
       bar alone false-positives on the thousands of sliver faces a *lattice*
       has (this function also runs on generated output during verification),
       and an absolute bar alone misses a large fractional loss on a large face.
    2. **Containment against the CAD bounding box, as an upper bound only.**
       Every mesh node must lie inside the face's box inflated by ``bbox_slack``.

    The bounding box is deliberately never used to test *coverage*. OCCT
    over-estimates it — control-point hull for B-spline edges, untrimmed UV
    rectangle for planes — and an earlier version of this gate that required the
    mesh to *reach* it rejected a perfectly valid input file outright
    (docs/algorithm.md §5.1). A conservative over-estimate is sound to test
    containment against and meaningless to test coverage against.
    """
    min_deficit = min_deficit_factor * min(lp.t, lp.a) ** 2
    for idx, face in enumerate(occ.faces(shape)):
        exact = occ.area(face)
        tri = occ.face_triangulation(face)
        if tri is None or len(tri[1]) == 0:
            raise InputGeometryError(
                f"Surface tessellation produced zero elements for face {idx} "
                f"(exact CAD area {exact:.4f} mm^2). The input geometry could not be "
                f"faithfully meshed for boundary classification."
            )
        verts, tris = tri
        a0, b0, c0 = verts[tris[:, 0]], verts[tris[:, 1]], verts[tris[:, 2]]
        meshed = float(
            0.5 * np.linalg.norm(np.cross(b0 - a0, c0 - a0), axis=1).sum()
        )
        if exact > 0 and meshed < (1 - area_rel_tol) * exact and (exact - meshed) > min_deficit:
            raise InputGeometryError(
                f"Surface tessellation for face {idx} covers only "
                f"{100 * meshed / exact:.2f}% of that face's exact area (meshed "
                f"{meshed:.4f} mm^2 vs. exact {exact:.4f} mm^2, shortfall "
                f"{exact - meshed:.4f} mm^2; tolerances {100 * area_rel_tol:.1f}% and "
                f"{min_deficit:.4f} mm^2). Classification against this mesh would "
                f"misclassify struts near the missing region as OUTSIDE rather than "
                f"fail. The input face likely needs repair in the originating CAD tool."
            )

        lo, hi = occ.bounding_box(face)
        if np.any(verts < lo - bbox_slack) or np.any(verts > hi + bbox_slack):
            raise InputGeometryError(
                f"Surface tessellation for face {idx} placed a mesh node outside that "
                f"face's own CAD bounding box lo={tuple(lo)} hi={tuple(hi)} (slack "
                f"{bbox_slack} mm). The triangulation is mis-parametrized and "
                f"classification against it would be unreliable."
            )


def measure_face_deviation(face) -> float:
    """How far the current triangulation of ``face`` departs from the surface.

    Sampled where a flat approximation is worst — the parametric midpoint of
    every triangle edge and the parametric centroid of every triangle — and
    measured as the distance from the true surface point to the *mesh element*
    (the chord, or the triangle), not to the mesh point carrying the same
    parameters.

    That distinction matters: comparing same-parameter points would measure how
    skewed the surface's parametrization is as much as how flat the mesh is.
    Measuring to the element removes most of that tangential component.

    The result is a deliberate **over-estimate**, for two reasons. It compares
    each sample only against the triangle it parametrically belongs to, so on a
    strongly non-uniform parametrization some lateral skew still leaks in; and
    it takes the worst sample on the face rather than an average. Both err the
    same way, and that way is the safe one — the value becomes the ``d`` of the
    classification margin, where too large merely sends more junctions down the
    (correct, slower) boundary path, while too small could let a strut that
    genuinely touches the surface be treated as clear of it.

    Measuring rather than assuming is not academic. Asked for 0.4 mm on
    ``80mm-test-ball.step``, OCCT's mesher leaves samples up to ~2.8 mm off the
    sphere; asked for 0.15 mm on ``TD_HX_Indre_Volum.step``, up to ~0.5 mm. A
    margin built on the requested figure would have been too small in both
    cases.

    Planar faces return 0 with no work — a plane's triangulation is exact.
    """
    if occ.is_planar(face):
        return 0.0
    data = occ.face_uv_triangulation(face)
    if data is None:
        return 0.0
    pts, uvs, tris, evaluate = data
    if len(tris) == 0:
        return 0.0

    edges = np.vstack([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]])
    edges = np.unique(np.sort(edges, axis=1), axis=0)
    worst = 0.0

    if len(edges):
        mids = evaluate(0.5 * (uvs[edges[:, 0]] + uvs[edges[:, 1]]))
        A, B = pts[edges[:, 0]], pts[edges[:, 1]]
        ab = B - A
        denom = np.einsum("ij,ij->i", ab, ab)
        safe = np.where(denom > 1e-300, denom, 1.0)
        s = np.clip(np.einsum("ij,ij->i", mids - A, ab) / safe, 0.0, 1.0)
        worst = float(np.linalg.norm(mids - (A + s[:, None] * ab), axis=1).max())

    # `_point_triangle_dist` broadcasts row-wise, so one point per triangle is a
    # single vectorized call.
    centroids = evaluate(uvs[tris].mean(axis=1))
    A, B, C = pts[tris[:, 0]], pts[tris[:, 1]], pts[tris[:, 2]]
    return max(worst, float(_point_triangle_dist(centroids, A, B, C).max()))


MESH_ANGLE = 0.2
"""Angular deflection in radians (~11.5°) passed to OCCT alongside the linear
one. It bounds how coarse a gently-curved face may get where linear deflection
alone would allow very long chords."""


def tessellate_surface(shape: TopoDS_Shape, lp: LatticeParams) -> TriMesh:
    """Mesh the input solid's boundary and return it as plain triangle data.

    ``min(t, a)/10`` is what OCCT is *asked* for as a linear chordal deviation.
    What the classification margin then uses is the deviation the mesher
    actually delivered, measured per face by :func:`occ.face_chordal_deviation`
    and returned on the mesh.

    That distinction is not pedantic. Measured on ``test-cylinder.STEP`` at
    ``cc=10, t=1.5``, asking for 0.15 mm produced a mesh whose real worst-case
    sagitta was around 0.4 mm — comfortably enough to make a margin of ``r +
    0.15`` too small, and enough to misclassify a node that a ten-times finer
    mesh puts on the boundary. Taking the larger of the two keeps the margin an
    upper bound on the mesh's own error, which is the property the whole
    "ambiguity always degrades to BOUNDARY" guarantee rests on.

    No maximum element size is needed: a planar face meshes to a couple of exact
    triangles however large it is, and only curvature drives refinement.
    """
    requested = chordal_target(lp)
    occ.mesh_shape(shape, requested, MESH_ANGLE)
    check_surface_mesh_coverage(shape, lp)

    all_verts: list[np.ndarray] = []
    all_tris: list[np.ndarray] = []
    offset = 0
    measured = 0.0
    for face in occ.faces(shape):
        measured = max(measured, measure_face_deviation(face))
        tri = occ.face_triangulation(face)
        if tri is None:
            continue
        verts, tris = tri
        all_verts.append(verts)
        all_tris.append(tris + offset)
        offset += len(verts)
    if not all_tris:
        raise InputGeometryError("Surface tessellation produced zero triangles.")
    verts, tris = _weld(np.vstack(all_verts), np.vstack(all_tris))
    if len(tris) == 0:
        raise InputGeometryError("Surface tessellation produced zero usable triangles.")

    deviation = max(requested, measured)
    if deviation > lp.a / 4.0:
        raise InputGeometryError(
            f"The surface tessellation deviates from the true geometry by up to "
            f"{deviation:.4f} mm, which is more than a quarter of the lattice cell "
            f"edge ({lp.a:.4f} mm). Classification against this mesh would be "
            f"dominated by meshing error rather than by the geometry."
        )
    return TriMesh(verts=verts, tris=tris, deviation=deviation)


# --- Distance primitives (vectorized over triangles) ------------------------


def _point_segment_dist(p: np.ndarray, A: np.ndarray, B: np.ndarray) -> np.ndarray:
    ab = B - A
    denom = np.einsum("ij,ij->i", ab, ab)
    safe = np.where(denom > 1e-300, denom, 1.0)
    t = np.clip(np.einsum("ij,ij->i", p - A, ab) / safe, 0.0, 1.0)
    t = np.where(denom > 1e-300, t, 0.0)
    return np.linalg.norm(p - (A + t[:, None] * ab), axis=1)


def _point_triangle_dist(p: np.ndarray, A: np.ndarray, B: np.ndarray, C: np.ndarray) -> np.ndarray:
    """Distance from one point to each of many triangles.

    The closest point is either the perpendicular foot (when its barycentric
    coordinates are interior) or the closest point on one of the three edges —
    the two cases are computed unconditionally and selected between, which
    vectorizes cleanly and avoids the branch-heavy region decomposition.
    """
    ab = B - A
    ac = C - A
    n = np.cross(ab, ac)
    nn = np.einsum("ij,ij->i", n, n)
    safe_nn = np.where(nn > 1e-300, nn, 1.0)
    ap = p - A
    # Barycentric coordinates of the perpendicular projection.
    d_ab_ab = np.einsum("ij,ij->i", ab, ab)
    d_ab_ac = np.einsum("ij,ij->i", ab, ac)
    d_ac_ac = np.einsum("ij,ij->i", ac, ac)
    d_ap_ab = np.einsum("ij,ij->i", ap, ab)
    d_ap_ac = np.einsum("ij,ij->i", ap, ac)
    det = d_ab_ab * d_ac_ac - d_ab_ac * d_ab_ac
    safe_det = np.where(np.abs(det) > 1e-300, det, 1.0)
    u = (d_ac_ac * d_ap_ab - d_ab_ac * d_ap_ac) / safe_det
    v = (d_ab_ab * d_ap_ac - d_ab_ac * d_ap_ab) / safe_det
    interior = (np.abs(det) > 1e-300) & (u >= 0) & (v >= 0) & (u + v <= 1) & (nn > 1e-300)
    perp = np.abs(np.einsum("ij,ij->i", ap, n)) / np.sqrt(safe_nn)
    edge = np.minimum(
        np.minimum(_point_segment_dist(p, A, B), _point_segment_dist(p, B, C)),
        _point_segment_dist(p, C, A),
    )
    return np.where(interior, perp, edge)


def _segment_segment_dist(p1: np.ndarray, q1: np.ndarray, P2: np.ndarray, Q2: np.ndarray) -> np.ndarray:
    """Distance from one segment to each of many segments (Ericson, RTCD §5.1.9)."""
    d1 = q1 - p1
    d2 = Q2 - P2
    r = p1 - P2
    a = float(np.dot(d1, d1))
    e = np.einsum("ij,ij->i", d2, d2)
    f = np.einsum("ij,ij->i", d2, r)
    c = r @ d1
    b = d2 @ d1
    EPS = 1e-12
    if a <= EPS:  # our segments are half-struts; never degenerate
        return np.linalg.norm(p1 - (P2 + np.clip(f / np.where(e > EPS, e, 1.0), 0, 1)[:, None] * d2), axis=1)
    denom = a * e - b * b
    s = np.where(denom > EPS, np.clip((b * f - c * e) / np.where(denom > EPS, denom, 1.0), 0.0, 1.0), 0.0)
    t = np.where(e > EPS, (b * s + f) / np.where(e > EPS, e, 1.0), 0.0)
    s = np.where(t < 0.0, np.clip(-c / a, 0.0, 1.0), s)
    s = np.where(t > 1.0, np.clip((b - c) / a, 0.0, 1.0), s)
    t = np.clip(t, 0.0, 1.0)
    degenerate = e <= EPS
    s = np.where(degenerate, np.clip(-c / a, 0.0, 1.0), s)
    t = np.where(degenerate, 0.0, t)
    c1 = p1 + s[:, None] * d1
    c2 = P2 + t[:, None] * d2
    return np.linalg.norm(c1 - c2, axis=1)


def _segment_triangle_intersects(s0, s1, A, B, C) -> np.ndarray:
    """Möller-Trumbore, bounded to ``t in [0,1]`` and orientation-agnostic."""
    EPS = 1e-9
    d = s1 - s0
    e1 = B - A
    e2 = C - A
    h = np.cross(d, e2)
    det = np.einsum("ij,ij->i", e1, h)
    ok = np.abs(det) >= EPS
    inv = np.where(ok, 1.0 / np.where(ok, det, 1.0), 0.0)
    s = s0 - A
    u = inv * np.einsum("ij,ij->i", s, h)
    q = np.cross(s, e1)
    v = inv * (q @ d)
    t = inv * np.einsum("ij,ij->i", e2, q)
    return ok & (u >= -EPS) & (u <= 1 + EPS) & (v >= -EPS) & (u + v <= 1 + EPS) & (t >= -EPS) & (t <= 1 + EPS)


def segment_triangle_dist(s0, s1, A, B, C) -> np.ndarray:
    """Exact distance from segment ``[s0,s1]`` to each triangle.

    A piercing segment has distance 0 and is checked first: the vertex/edge
    feature decomposition below is only valid for the disjoint case.
    """
    d = np.minimum(_point_triangle_dist(s0, A, B, C), _point_triangle_dist(s1, A, B, C))
    d = np.minimum(d, _segment_segment_dist(s0, s1, A, B))
    d = np.minimum(d, _segment_segment_dist(s0, s1, B, C))
    d = np.minimum(d, _segment_segment_dist(s0, s1, C, A))
    return np.where(_segment_triangle_intersects(s0, s1, A, B, C), 0.0, d)


# --- Spatial hash -----------------------------------------------------------


_KEY_OFFSET = 1 << 20
_KEY_BASE = 1 << 21


def _pack_cells(cells: np.ndarray) -> np.ndarray:
    """Pack signed integer cell coordinates into one int64 key each."""
    shifted = cells + _KEY_OFFSET
    if np.any(shifted < 0) or np.any(shifted >= _KEY_BASE):
        raise InputGeometryError(
            "Geometry extends beyond the spatial index's addressable range; the "
            "input is implausibly large relative to the requested cell size."
        )
    return (shifted[:, 0] * _KEY_BASE + shifted[:, 1]) * _KEY_BASE + shifted[:, 2]


def _cell_assignments(lo: np.ndarray, hi: np.ndarray, cell: float):
    """``(item_ids, cells)`` for every grid cell each AABB overlaps.

    Vectorized ragged expansion — the alternative, a Python loop over items with
    a nested loop over cells, costs tens of seconds once the mesh reaches the
    hundreds of thousands of triangles a *lattice* has (this index is also built
    over generated output during verification, not just over input CAD).
    """
    klo = np.floor(lo / cell).astype(np.int64)
    khi = np.floor(hi / cell).astype(np.int64)
    extent = khi - klo + 1
    counts = extent[:, 0] * extent[:, 1] * extent[:, 2]
    total = int(counts.sum())
    item_ids = np.repeat(np.arange(len(counts), dtype=np.int64), counts)
    starts = np.cumsum(counts) - counts
    within = np.arange(total, dtype=np.int64) - np.repeat(starts, counts)
    ny = extent[item_ids, 1]
    nz = extent[item_ids, 2]
    ox = within // (ny * nz)
    rem = within % (ny * nz)
    oy = rem // nz
    oz = rem % nz
    cells = klo[item_ids] + np.stack([ox, oy, oz], axis=1)
    return item_ids, cells


class SpatialHash:
    """Uniform grid over triangle AABBs, for exact near-surface queries.

    Cell size follows the triangles (2x the median edge), but ``min_cell`` lets a
    caller stop it going finer than its own queries need. That floor matters:
    query cost is proportional to the number of cells the query box covers, and
    the box is a half-strut inflated by the margin — roughly the cell edge ``a``.
    On a part where the mesh is far finer than the lattice (large ``cc`` with
    small ``t``, where ``d = t/10`` drives a fine tessellation of a coarse
    lattice), an edge-derived cell size alone would put hundreds of cells along
    each axis of every query and turn an index lookup into a sweep of the whole
    grid. A coarser cell only means more candidate triangles per query, which is
    a graceful, bounded cost.
    """

    def __init__(self, mesh: TriMesh, min_cell: float = 0.0):
        a, b, c = mesh.triangle_points
        edges = np.concatenate(
            [
                np.linalg.norm(b - a, axis=1),
                np.linalg.norm(c - b, axis=1),
                np.linalg.norm(a - c, axis=1),
            ]
        )
        self.cell = max(2.0 * float(np.median(edges)), float(min_cell), 1e-9)
        self.mesh = mesh
        lo = np.minimum(np.minimum(a, b), c)
        hi = np.maximum(np.maximum(a, b), c)
        item_ids, cells = _cell_assignments(lo, hi, self.cell)
        keys = _pack_cells(cells)
        order = np.argsort(keys, kind="stable")
        self._keys = keys[order]
        self._items = item_ids[order]
        self._unique, self._starts = np.unique(self._keys, return_index=True)

    def cell_groups(self):
        """Yield the triangle indices sharing each occupied cell.

        Two triangles whose AABBs overlap always share at least one cell, so
        sweeping cells and testing within each is equivalent to querying per
        triangle — and far cheaper when *every* triangle has to be checked, as
        in the self-intersection check.
        """
        for pos in range(len(self._unique)):
            start = self._starts[pos]
            end = self._starts[pos + 1] if pos + 1 < len(self._starts) else len(self._items)
            if end - start > 1:
                yield self._items[start:end]

    def query(self, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
        """Triangle indices whose cells overlap the AABB ``[lo, hi]``."""
        klo = np.floor(np.asarray(lo) / self.cell).astype(np.int64)
        khi = np.floor(np.asarray(hi) / self.cell).astype(np.int64)
        ranges = [np.arange(klo[i], khi[i] + 1) for i in range(3)]
        grid = np.stack([g.ravel() for g in np.meshgrid(*ranges, indexing="ij")], axis=1)
        wanted = _pack_cells(grid)
        idx = np.searchsorted(self._unique, wanted)
        keep = idx < len(self._unique)
        idx, wanted = idx[keep], wanted[keep]
        keep = self._unique[idx] == wanted
        idx = idx[keep]
        if len(idx) == 0:
            return np.empty(0, dtype=np.int64)
        found = []
        for pos in idx:
            start = self._starts[pos]
            end = self._starts[pos + 1] if pos + 1 < len(self._starts) else len(self._items)
            found.append(self._items[start:end])
        return np.unique(np.concatenate(found))


class Occupancy:
    """Coarse "could this point be near the surface?" index.

    Deliberately a separate grid from :class:`SpatialHash`, sized so that one
    cell already covers the query radius. That keeps every query a fixed
    27-cell neighbourhood lookup instead of a radius-dependent sweep, and the
    whole sweep vectorizes over all points at once.

    The answer is conservative by construction: a ``True`` may still turn out to
    be far away, which only costs the exact test that follows. A false ``False``
    would skip that test, and cannot happen — the cell size is at least the
    radius, so anything within the radius is in the 27-neighbourhood.
    """

    def __init__(self, mesh: TriMesh, radius: float):
        a, b, c = mesh.triangle_points
        self.cell = max(float(radius), 1e-9)
        lo = np.minimum(np.minimum(a, b), c)
        hi = np.maximum(np.maximum(a, b), c)
        _, cells = _cell_assignments(lo, hi, self.cell)
        self._occupied = np.unique(_pack_cells(cells))

    def near(self, points: np.ndarray) -> np.ndarray:
        out = np.zeros(len(points), dtype=bool)
        if len(self._occupied) == 0:
            return out
        base = np.floor(points / self.cell).astype(np.int64)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    keys = _pack_cells(base + np.array([dx, dy, dz]))
                    pos = np.searchsorted(self._occupied, keys)
                    pos = np.clip(pos, 0, len(self._occupied) - 1)
                    out |= self._occupied[pos] == keys
        return out


# --- Point-in-solid by ray parity -------------------------------------------

RAY_DIRS = np.array(
    [
        [0.8574, 0.2196, 0.4671],
        [-0.3157, 0.8842, 0.3427],
        [0.1732, -0.4899, 0.8547],
    ]
)
"""Three fixed, mutually non-parallel ray directions (docs/algorithm.md §5.2).

Fixed rather than sampled so a run is exactly reproducible; chosen off-axis and
off-diagonal so they avoid grazing the axis-aligned and 45° features typical of
CAD geometry. Three rays with a majority vote defeat the single-ray degeneracies
(a ray through a shared edge or vertex is counted by both incident triangles).
"""


class _DirectionalCaster:
    """Ray-parity counter for one fixed direction.

    Triangles are projected onto the plane perpendicular to the ray and bucketed
    into a 2D grid, so a query only tests triangles whose projection could
    contain the query point. Points are then processed bucket by bucket, which
    turns the whole sweep into a handful of vectorized ``(points x triangles)``
    Möller-Trumbore evaluations instead of a per-point Python loop.
    """

    def __init__(self, mesh: TriMesh, direction: np.ndarray, target_per_cell: float = 6.0):
        self.d = direction / np.linalg.norm(direction)
        ref = np.array([0.0, 0.0, 1.0]) if abs(self.d[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
        self.p = np.cross(self.d, ref)
        self.p /= np.linalg.norm(self.p)
        self.q = np.cross(self.d, self.p)
        self.A, self.B, self.C = mesh.triangle_points
        basis = np.stack([self.p, self.q])
        pa, pb, pc = self.A @ basis.T, self.B @ basis.T, self.C @ basis.T
        lo2 = np.minimum(np.minimum(pa, pb), pc)
        hi2 = np.maximum(np.maximum(pa, pb), pc)
        diag = np.linalg.norm(hi2 - lo2, axis=1)
        self.cell = max(float(np.median(diag)) * float(np.sqrt(target_per_cell)), 1e-9)
        self.basis = basis
        klo = np.floor(lo2 / self.cell).astype(np.int64)
        khi = np.floor(hi2 / self.cell).astype(np.int64)
        buckets: dict[tuple[int, int], list[int]] = {}
        for ti in range(len(pa)):
            for ix in range(klo[ti, 0], khi[ti, 0] + 1):
                for iy in range(klo[ti, 1], khi[ti, 1] + 1):
                    buckets.setdefault((ix, iy), []).append(ti)
        self.buckets = {k: np.array(v, dtype=np.int64) for k, v in buckets.items()}

    def hit_counts(self, points: np.ndarray) -> np.ndarray:
        """Number of forward triangle crossings for each point."""
        keys = np.floor((points @ self.basis.T) / self.cell).astype(np.int64)
        counts = np.zeros(len(points), dtype=np.int64)
        order = np.lexsort((keys[:, 1], keys[:, 0]))
        start = 0
        while start < len(order):
            end = start + 1
            k0 = keys[order[start]]
            while end < len(order) and keys[order[end], 0] == k0[0] and keys[order[end], 1] == k0[1]:
                end += 1
            idx = order[start:end]
            tris = self.buckets.get((int(k0[0]), int(k0[1])))
            if tris is not None and len(tris):
                counts[idx] = self._count(points[idx], tris)
            start = end
        return counts

    def _count(self, origins: np.ndarray, tris: np.ndarray) -> np.ndarray:
        EPS = 1e-9
        A = self.A[tris]
        e1 = self.B[tris] - A
        e2 = self.C[tris] - A
        h = np.cross(self.d, e2)
        det = np.einsum("ij,ij->i", e1, h)
        ok = np.abs(det) >= EPS
        inv = np.where(ok, 1.0 / np.where(ok, det, 1.0), 0.0)
        s = origins[:, None, :] - A[None, :, :]
        u = inv[None, :] * np.einsum("pij,ij->pi", s, h)
        qv = np.cross(s, e1[None, :, :])
        v = inv[None, :] * (qv @ self.d)
        t = inv[None, :] * np.einsum("pij,ij->pi", qv, e2)
        hit = (
            ok[None, :]
            & (u >= -EPS)
            & (u <= 1 + EPS)
            & (v >= -EPS)
            & (u + v <= 1 + EPS)
            & (t > EPS)
        )
        return hit.sum(axis=1)


class PointInside:
    """Majority-vote point-in-solid test over the three fixed ray directions."""

    def __init__(self, mesh: TriMesh):
        self.casters = [_DirectionalCaster(mesh, d) for d in RAY_DIRS]

    def __call__(self, points: np.ndarray) -> np.ndarray:
        votes = np.zeros(len(points), dtype=np.int64)
        for caster in self.casters:
            votes += (caster.hit_counts(points) % 2).astype(np.int64)
        return votes >= 2


# --- Top-level classification -----------------------------------------------


@dataclass
class Classification:
    """Per-node classes plus the counts the run log reports."""

    node_index: np.ndarray
    """``(N, 3)`` integer node indices, in the same order as ``node_class``."""
    node_class: np.ndarray
    """``(N,)`` :class:`Klass` values."""
    n_triangles: int

    def of(self, k: Klass) -> np.ndarray:
        return self.node_index[self.node_class == int(k)]

    def counts(self) -> dict[str, int]:
        return {
            "candidates": int(len(self.node_class)),
            "interior": int(np.count_nonzero(self.node_class == Klass.INTERIOR)),
            "boundary": int(np.count_nonzero(self.node_class == Klass.BOUNDARY)),
            "outside": int(np.count_nonzero(self.node_class == Klass.OUTSIDE)),
        }


def classify_nodes(lp: LatticeParams, mesh: TriMesh, candidates: np.ndarray) -> Classification:
    """Classify every candidate node as INTERIOR, BOUNDARY or OUTSIDE.

    The sweep is staged cheapest-first:

    1. A grid-resolution occupancy query rules out every node whose junction
       cannot come near the surface at all (the overwhelming majority — the
       interior and the padding around the part).
    2. Those far nodes are decided purely by the inside/outside test, since none
       of their six half-struts can touch the surface.
    3. Only the remaining near-surface nodes pay for exact segment-triangle
       distances, six segments each.
    """
    margin = lp.r + mesh.deviation
    # Keep the grid no finer than a quarter of the largest query box, so each
    # half-strut query touches a handful of cells regardless of mesh density.
    sh = SpatialHash(mesh, min_cell=(lp.a / 2.0 + 2.0 * margin) / 4.0)
    positions = nodes(lp, candidates)

    inside = PointInside(mesh)(positions)

    # A junction reaches at most a/2 from its node, so anything further than
    # a/2 + margin from the surface has all six half-struts provably clear.
    near = Occupancy(mesh, lp.a / 2.0 + margin).near(positions)

    node_class = np.where(inside, int(Klass.INTERIOR), int(Klass.OUTSIDE)).astype(np.int64)

    offsets = np.array([half_strut_offset(lp, h) for h in range(6)])
    A, B, C = mesh.triangle_points
    for i in np.nonzero(near)[0]:
        p0 = positions[i]
        touched = False
        for h in range(6):
            p1 = p0 + offsets[h]
            lo = np.minimum(p0, p1) - margin
            hi = np.maximum(p0, p1) + margin
            cand = sh.query(lo, hi)
            if len(cand) == 0:
                continue
            dist = segment_triangle_dist(p0, p1, A[cand], B[cand], C[cand])
            if float(dist.min()) <= margin:
                touched = True
                break
        if touched:
            node_class[i] = int(Klass.BOUNDARY)

    return Classification(node_index=candidates, node_class=node_class,
                          n_triangles=int(len(mesh.tris)))
