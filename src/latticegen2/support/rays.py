"""First-contact rays (docs/algorithm.md §14.3 S3a). Pure NumPy.

From an apex ``x`` inside the body, the closed cone of half-angle exactly the
lattice rod angle about ``w = -b`` grows downward. Its **first contact** is the
obstacle point ``y`` inside the cone with the smallest drop ``z = (y - x) . w``;
ties go to the smallest distance ``|y - x|``, then to the point straight below.
That point is the ray's target: how deep the support under ``x`` must reach.

The query is exact. Drop is linear, so over (triangle ∩ cone) -- a compact
convex set -- its minimum lies at an extreme point, and every extreme point is
one of:

1. a triangle corner inside the cone;
2. a triangle edge crossing the cone's surface (a quadratic in the edge
   parameter);
3. the point where the triangle's plane meets the cone's *first-contact
   generator*: the generator whose horizontal part runs along the horizontal
   part of the plane normal (oriented away from the apex), when it falls
   inside the triangle;
4. for a horizontal plane, where drop is constant over the whole disc, the
   point straight below -- which is also what the tie-break picks;
5. the apex itself, when it lies on the triangle: drop zero. That is not an
   exotic case. A ceiling vertex on the crease where the roof meets a vertical
   wall has its apex on the wall's plane, where the cone meets the plane in two
   generators through the apex and the quadratic of case 2 degenerates to a
   double root -- whose discriminant rounds either side of zero. Measured on
   `Hole-block-basic`'s outer border before this case existed: the wall was
   missed and the corner rays reached 8.94 mm down to the nearest tube.

Candidates are pruned with a :class:`latticegen2.classify.SpatialHash` over the
obstacles, inside a drop bound that doubles until something is found: anything
within drop ``Z`` of the apex lies in the truncated cone's bounding box, so the
first candidate found inside a bound is the global first contact.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..classify import SpatialHash, TriMesh
from ..errors import ProcessingError
from .params import SUPPORT_OVERHANG_DEG

COS_T = float(np.cos(np.radians(SUPPORT_OVERHANG_DEG)))
SIN_T = float(np.sin(np.radians(SUPPORT_OVERHANG_DEG)))
TAN_T = SIN_T / COS_T

APEX_NUDGE = 1e-3
"""mm. How far into the body a ray's apex sits below its ceiling point, so the
ceiling it starts on is not its own first contact."""

TIE_TOL = 1e-9
"""mm. Drops this close count as a tie, decided by distance."""

INSIDE_TOL = 1e-12
"""Barycentric slack when asking whether a plane point falls in its triangle."""

ON_PLANE_TOL = 1e-9
"""mm. An apex this close to a triangle's plane, inside it, is on it."""

DISC_REL_TOL = 1e-12
"""Relative slack on case 2's discriminant, for the double root of an edge
running along a generator."""


@dataclass
class Obstacles:
    """Everything a ray can land on: triangles, each labelled with its shape."""

    mesh: TriMesh
    anchor: np.ndarray
    """``(T,)`` int: 0 for the input body, ``k + 1`` for lattice solid ``k``."""
    key: np.ndarray = None  # type: ignore[assignment]
    """``(T,)`` int: the anchor group (:mod:`.obstacles`) -- the unit a support
    cell may land on. Defaults to the shape id."""

    def __post_init__(self):
        if self.key is None:
            self.key = np.asarray(self.anchor).copy()
        self.hash = SpatialHash(self.mesh)
        a, b, c = self.mesh.triangle_points
        self._A, self._B, self._C = a, b, c
        n = np.cross(b - a, c - a)
        self._n = n / np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-300)
        self._lo = np.minimum(np.minimum(a, b), c)
        self._hi = np.maximum(np.maximum(a, b), c)
        lo, hi = self.mesh.bounds()
        self.diagonal = float(np.linalg.norm(hi - lo))

    def subset(self, mask: np.ndarray) -> "Obstacles":
        """The obstacles where ``mask`` holds: a cast restricted to one anchor."""
        idx = np.nonzero(mask)[0]
        tris = self.mesh.tris[idx]
        used, inv = np.unique(tris.ravel(), return_inverse=True)
        sub = Obstacles(TriMesh(self.mesh.verts[used], inv.reshape(-1, 3)),
                        self.anchor[idx], self.key[idx])
        sub.parent_index = idx
        sub.diagonal = self.diagonal
        return sub


@dataclass
class Rays:
    """One first-contact ray per apex."""

    apex: np.ndarray
    """``(N, 3)``."""
    target: np.ndarray
    """``(N, 3)`` the first-contact point."""
    drop: np.ndarray
    """``(N,)`` ``(target - apex) . w``."""
    length: np.ndarray
    """``(N,)`` ``|target - apex|``."""
    anchor: np.ndarray
    """``(N,)`` the anchor shape the ray lands on (:attr:`Obstacles.anchor`)."""
    tri: np.ndarray
    """``(N,)`` index of the obstacle triangle it lands on."""
    key: np.ndarray = None  # type: ignore[assignment]
    """``(N,)`` the anchor group it lands on (:attr:`Obstacles.key`)."""

    def __len__(self) -> int:
        return len(self.apex)


def _inside(P, A, B, C, n) -> np.ndarray:
    """Is ``P`` (on the plane of ``A B C``) inside the triangle?"""
    c0 = np.einsum("ij,ij->i", np.cross(B - A, P - A), n)
    c1 = np.einsum("ij,ij->i", np.cross(C - B, P - B), n)
    c2 = np.einsum("ij,ij->i", np.cross(A - C, P - C), n)
    scale = np.einsum("ij,ij->i", np.cross(B - A, C - A), n)
    tol = INSIDE_TOL * np.abs(scale)
    return (c0 >= -tol) & (c1 >= -tol) & (c2 >= -tol)


def _candidates(X: np.ndarray, w: np.ndarray, A, B, C, n):
    """Per (apex, triangle) pair: ``(drop, point, distance)`` of its first contact.

    Row ``i`` pairs apex ``X[i]`` with triangle ``A[i] B[i] C[i]`` (normal
    ``n[i]``); ``drop`` is ``inf`` where the triangle has no point in the cone.
    """
    m = len(A)
    best_z = np.full(m, np.inf)
    best_p = np.zeros((m, 3))
    best_d = np.full(m, np.inf)

    def offer(mask, z, p):
        if not mask.any():
            return
        idx = np.nonzero(mask)[0]
        z, p = z[mask], p[mask]
        d = np.linalg.norm(p - X[idx], axis=1)
        better = (z < best_z[idx] - TIE_TOL) | (
            (np.abs(z - best_z[idx]) <= TIE_TOL) & (d < best_d[idx]))
        idx, z, p, d = idx[better], z[better], p[better], d[better]
        best_z[idx], best_p[idx], best_d[idx] = z, p, d

    def in_cone(P):
        q = P - X
        z = q @ w
        h = np.linalg.norm(q - np.outer(z, w), axis=1)
        return z, (z >= 0.0) & (h <= z * TAN_T * (1.0 + 1e-12) + 1e-12)

    # 1. corners inside the cone.
    for V in (A, B, C):
        z, ok = in_cone(V)
        offer(ok, z, V)

    # 5. the apex on the triangle.
    c_ = np.einsum("ij,ij->i", n, A - X)
    on = np.abs(c_) <= ON_PLANE_TOL
    if on.any():
        on &= _inside(X - c_[:, None] * n, A, B, C, n)
        offer(on, np.zeros(m), X)

    # 2. edges crossing the cone surface: |H(s)|^2 = tan^2 * z(s)^2.
    for P0, P1 in ((A, B), (B, C), (C, A)):
        q, d = P0 - X, P1 - P0
        qz, dz = q @ w, d @ w
        qh, dh = q - np.outer(qz, w), d - np.outer(dz, w)
        a2 = np.einsum("ij,ij->i", dh, dh) - TAN_T ** 2 * dz ** 2
        a1 = 2.0 * (np.einsum("ij,ij->i", qh, dh) - TAN_T ** 2 * qz * dz)
        a0 = np.einsum("ij,ij->i", qh, qh) - TAN_T ** 2 * qz ** 2
        disc = a1 ** 2 - 4.0 * a2 * a0
        # A double root -- an edge through the apex's generator -- rounds
        # either side of zero; within that rounding it is a root.
        real = (disc >= -DISC_REL_TOL * (a1 ** 2 + np.abs(4.0 * a2 * a0))) & (
            np.abs(a2) > 1e-300)
        root = np.sqrt(np.where(real, np.maximum(disc, 0.0), 0.0))
        for sgn in (-1.0, 1.0):
            s = np.where(real, (-a1 + sgn * root) / np.where(real, 2.0 * a2, 1.0), -1.0)
            ok = real & (s >= 0.0) & (s <= 1.0)
            z = qz + s * dz
            ok &= z >= 0.0
            offer(ok, z, P0 + s[:, None] * d)
        # A degenerate quadratic (edge parallel to a generator) is linear.
        lin = (np.abs(a2) <= 1e-300) & (np.abs(a1) > 1e-300)
        if lin.any():
            s = np.where(lin, -a0 / np.where(lin, a1, 1.0), -1.0)
            ok = lin & (s >= 0.0) & (s <= 1.0)
            z = qz + s * dz
            offer(ok & (z >= 0.0), z, P0 + s[:, None] * d)

    # 3. the plane's first-contact generator, and 4. straight below.
    c = c_
    m_ = np.where((c < 0.0)[:, None], -n, n)
    c = np.abs(c)
    mw = m_ @ w
    mh = m_ - np.outer(mw, w)
    mhn = np.linalg.norm(mh, axis=1)
    reach = COS_T * mw + SIN_T * mhn
    ok = (c > 1e-15) & (reach > 1e-12) & (mhn > 1e-12)
    g = COS_T * w + SIN_T * mh / np.maximum(mhn, 1e-300)[:, None]
    lam = c / np.where(ok, reach, 1.0)
    P = X + lam[:, None] * g
    ok &= _inside(P, A, B, C, n)
    offer(ok, lam * COS_T, P)
    ok = (c > 1e-15) & (mw > 1e-12)
    s = c / np.where(ok, mw, 1.0)
    P = X + s[:, None] * w
    ok &= _inside(P, A, B, C, n)
    offer(ok, s, P)
    return best_z, best_p, best_d


def _cone_box(x: np.ndarray, w: np.ndarray, Z):
    """World AABB of the cone truncated at drop ``Z``, for one apex or many."""
    x = np.atleast_2d(x)
    Z = np.broadcast_to(np.asarray(Z, dtype=float), (len(x),))[:, None]
    centre = x + Z * w
    ext = Z * TAN_T * np.sqrt(np.clip(1.0 - w ** 2, 0.0, 1.0))
    return np.minimum(x, centre - ext), np.maximum(x, centre + ext)


BUCKET = 2.0
"""mm. Apexes are cast together in cubes of this side: one index query and one
vectorised candidate pass serve every apex in the cube."""

PAIRS_MAX = 400_000
"""Apex x triangle pairs per vectorised pass, to bound the pass's memory."""


def first_contact(apex: np.ndarray, b: np.ndarray, obstacles: Obstacles,
                  start: float = 0.5, report=None, allow_miss: bool = False) -> Rays:
    """The first-contact ray from every row of ``apex`` (docs/algorithm.md §14.3 S3a).

    Apexes are bucketed in cubes of :data:`BUCKET`; a bucket shares one index
    query over the union of its members' cone boxes, each member at its own
    drop bound. An apex is settled once a candidate lies within its own bound:
    the guarantee is per apex, so a larger shared query only adds candidates
    and never removes the one that decides.
    """
    apex = np.asarray(apex, dtype=float)
    w = -np.asarray(b, dtype=float)
    w = w / np.linalg.norm(w)
    n = len(apex)
    target = np.zeros((n, 3))
    drop = np.full(n, np.inf)
    length = np.full(n, np.inf)
    anchor = np.full(n, -1, dtype=np.int64)
    key_ = np.full(n, -1, dtype=np.int64)
    tri = np.full(n, -1, dtype=np.int64)
    bound = np.full(n, float(start))
    limit = 2.0 * obstacles.diagonal
    A, B, C, N = obstacles._A, obstacles._B, obstacles._C, obstacles._n
    TLO, THI = obstacles._lo, obstacles._hi
    pending = np.arange(n)
    while len(pending):
        if report is not None:
            report("casting first-contact rays", n - len(pending), n)
        keys = np.floor(apex[pending] / BUCKET).astype(np.int64)
        _u, group = np.unique(keys, axis=0, return_inverse=True)
        group = group.reshape(-1)
        order = np.argsort(group, kind="stable")
        starts = np.r_[0, np.nonzero(np.diff(group[order]))[0] + 1, len(order)]
        for g0, g1 in zip(starts[:-1], starts[1:]):
            members = pending[order[g0:g1]]
            lo, hi = _cone_box(apex[members], w, bound[members])
            cand = obstacles.hash.query(lo.min(axis=0), hi.max(axis=0))
            if not len(cand):
                continue
            chunk = max(1, PAIRS_MAX // len(cand))
            for c0 in range(0, len(members), chunk):
                mem = members[c0:c0 + chunk]
                ia = np.repeat(mem, len(cand))
                it = np.tile(cand, len(mem))
                # Cheap and exact pruning before the candidate pass: the
                # triangle's box must meet this apex's own cone box, and its
                # least corner drop -- a lower bound on any drop over it, drop
                # being linear -- must be within this apex's bound.
                alo, ahi = lo[c0:c0 + chunk], hi[c0:c0 + chunk]
                k = np.repeat(np.arange(len(mem)), len(cand))
                keep = np.all((TLO[it] <= ahi[k]) & (THI[it] >= alo[k]), axis=1)
                X = apex[ia]
                lb = np.minimum(np.minimum((A[it] - X) @ w, (B[it] - X) @ w),
                                (C[it] - X) @ w)
                keep &= lb <= bound[ia]
                ia, it = ia[keep], it[keep]
                if not len(ia):
                    continue
                z, p, d = _candidates(apex[ia], w, A[it], B[it], C[it], N[it])
                z = np.where(z <= bound[ia], z, np.inf)
                # Per apex: least drop, ties to least distance.
                zq = np.round(z / TIE_TOL)
                key = np.lexsort((d, zq, ia))
                first = np.r_[True, ia[key][1:] != ia[key][:-1]]
                win = key[first]
                win = win[np.isfinite(z[win])]
                a_ = ia[win]
                target[a_], drop[a_], length[a_] = p[win], z[win], d[win]
                tri[a_] = it[win]
                anchor[a_] = obstacles.anchor[it[win]]
                key_[a_] = obstacles.key[it[win]]
        settled = np.isfinite(drop[pending])
        rest = pending[~settled]
        if allow_miss and len(rest):
            # A cast restricted to one anchor may legitimately find nothing:
            # from beyond a parting, the other cell's anchor need not be in the
            # cone at all. Those apexes keep an infinite drop and no target.
            gone = rest[bound[rest] > limit]
            target[gone] = np.nan
            rest = rest[bound[rest] <= limit]
        if len(rest) and (bound[rest] > limit).any():
            x = apex[rest[np.argmax(bound[rest])]]
            raise ProcessingError(
                f"A first-contact ray from {np.round(x, 3).tolist()} found "
                f"nothing below it: the body's surface mesh does not enclose "
                f"that point (docs/algorithm.md §14.3 S3a)."
            )
        bound[rest] *= 2.0
        pending = rest
    if report is not None:
        report("casting first-contact rays", n, n)
    parent = getattr(obstacles, "parent_index", None)
    if parent is not None:
        tri = np.where(tri >= 0, parent[np.maximum(tri, 0)], -1)
    return Rays(apex, target, drop, length, anchor, tri, key_)


def brute_force(x: np.ndarray, b: np.ndarray, obstacles: Obstacles,
                polar: int = 200, azimuth: int = 720):
    """Reference answer for tests: a dense fan of rays, nearest drop wins.

    Approximate by construction -- it samples the cone -- so it is only ever an
    upper bound on the exact drop, approached as the fan grows.
    """
    w = -np.asarray(b, dtype=float)
    w /= np.linalg.norm(w)
    e1 = np.cross(w, [1.0, 0.0, 0.0])
    if np.linalg.norm(e1) < 1e-6:
        e1 = np.cross(w, [0.0, 1.0, 0.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(w, e1)
    A, B, C = obstacles._A, obstacles._B, obstacles._C
    best = (np.inf, None)
    for alpha in np.linspace(0.0, np.radians(SUPPORT_OVERHANG_DEG), polar):
        for phi in np.linspace(0.0, 2 * np.pi, azimuth, endpoint=False):
            d = np.cos(alpha) * w + np.sin(alpha) * (np.cos(phi) * e1 + np.sin(phi) * e2)
            # Moller-Trumbore against every triangle.
            e_1, e_2 = B - A, C - A
            pv = np.cross(d, e_2)
            det = np.einsum("ij,ij->i", e_1, pv)
            ok = np.abs(det) > 1e-14
            inv = np.where(ok, 1.0 / np.where(ok, det, 1.0), 0.0)
            tv = x - A
            u = np.einsum("ij,ij->i", tv, pv) * inv
            qv = np.cross(tv, e_1)
            v = (qv @ d) * inv
            t = np.einsum("ij,ij->i", e_2, qv) * inv
            hit = ok & (u >= 0) & (v >= 0) & (u + v <= 1) & (t > 0)
            if hit.any():
                tt = t[hit].min()
                z = tt * float(d @ w)
                if z < best[0]:
                    best = (z, x + tt * d)
    return best
