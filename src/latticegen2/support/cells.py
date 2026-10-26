"""Support cells and the partings between them (docs/algorithm.md §14.3 S3b, S3c).

A **cell** is a patch of one overhang region whose rays all land on one anchor
group and vary continuously: one support body will carry it. Cells meet on
**parting curves** -- where two anchors are reached at the same drop -- computed
exactly rather than read off the mesh, so that two neighbouring bodies share
one boundary curve on the ceiling instead of two copies of a staircase.

What comes out per cell is its perimeter: closed loops of points on the true
ceiling surface, each with **this cell's** ray from it. That is all a cell solid
needs (§14.3 S4).
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np

from .params import REACH, STEP
from .rays import APEX_NUDGE, first_contact

KAPPA = 4.0
"""Continuity: neighbours stay in one cell while their targets are no further
apart than ``KAPPA`` times their own spacing, plus ``2 h``."""

PARTING_TOL = 1e-4
"""mm. How exactly a parting point is located along its mesh edge."""

SPLIT_ROUNDS_MAX = 16
"""Rounds of splitting cells that enclose another cell before giving up."""

FOLD_EDGE_DEG = 20.0
"""Degrees. Two anchor facets meeting at a sharper angle than this make a
feature edge -- a ledge's rim -- that rays can converge on."""

FOLD_TOL = 0.5 * STEP
"""mm. A ray landing this close to a convex edge of its anchor is a converging
one; the boundary between converging and following rays is where a target is
exactly this far from the edge."""

FOLD_SAMPLE = 0.25 * FOLD_TOL
"""mm. Spacing of the points feature edges are sampled at."""

EXEMPT = -1
OUTSIDE = -2

VERTEX, PARTING, EXEMPT_EDGE, TRIPLE, CUT, BORDER, FOLD = 0, 1, 2, 3, 4, 5, 6


@dataclass
class Loop:
    """One closed perimeter loop of a cell, wound with the ceiling's outward normal."""

    points: np.ndarray
    """``(n, 3)`` on the true ceiling surface."""
    normals: np.ndarray
    """``(n, 3)`` outward normals there."""
    targets: np.ndarray
    """``(n, 3)`` where this cell's ray from each point lands."""
    kind: np.ndarray
    """``(n,)`` :data:`VERTEX`, :data:`PARTING`, :data:`EXEMPT_EDGE`,
    :data:`TRIPLE` or :data:`CUT`."""
    keys: list = field(default_factory=list)
    """Point identity per entry, shared with the neighbouring cell's loop at a
    parting, so the two bodies can be checked to share their boundary."""


@dataclass
class Cell:
    id: int
    key: int
    """The anchor group every ray of this cell lands on."""
    verts: np.ndarray
    """Region-mesh vertex indices in the cell."""
    area: float
    loops: list = field(default_factory=list)
    """:class:`Loop` s, the outer one first."""
    seed: np.ndarray = None  # type: ignore[assignment]
    """A point a hair inside the body under the cell: the cell's solid is the
    piece of its trimmed loft that holds this point."""
    min_divergence_deg: float = 180.0
    """Smallest angle between this cell's ray and a neighbour's at a parting."""


@dataclass
class Partition:
    labels: np.ndarray
    """``(V,)`` cell id per region-mesh vertex, :data:`EXEMPT` or :data:`OUTSIDE`."""
    cells: list
    stats: dict = field(default_factory=dict)


# --- helpers ----------------------------------------------------------------


def _union_find(n: int, pairs: np.ndarray) -> np.ndarray:
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


class _Caster:
    """First-contact casts restricted to one anchor group, cached per group."""

    def __init__(self, obstacles, b):
        self.obstacles = obstacles
        self.b = np.asarray(b, dtype=float)
        self._subs = {}

    def __call__(self, points, normals, key: int):
        sub = self._subs.get(key)
        if sub is None:
            sub = self._subs[key] = self.obstacles.subset(self.obstacles.key == key)
        if not len(points):
            return np.empty((0, 3)), np.empty(0), np.empty(0)
        got = first_contact(np.asarray(points) - APEX_NUDGE * np.asarray(normals),
                            self.b, sub, allow_miss=True)
        return got.target, got.drop, got.length

    def full(self, points, normals):
        got = first_contact(np.asarray(points) - APEX_NUDGE * np.asarray(normals),
                            self.b, self.obstacles)
        return got.target, got.drop, got.length


def _edge_triangles(tris: np.ndarray) -> dict:
    """``{(i, j) sorted: (tri, k_i, k_j)}`` -- one triangle using each edge."""
    out = {}
    for t, (a, b, c) in enumerate(tris):
        for (p, q), (kp, kq) in (((a, b), (0, 1)), ((b, c), (1, 2)), ((c, a), (2, 0))):
            key = (p, q) if p < q else (q, p)
            if key not in out:
                out[key] = (t, kp, kq) if p < q else (t, kq, kp)
    return out


def _bary(k_i: int, k_j: int, s: np.ndarray) -> np.ndarray:
    out = np.zeros((len(s), 3))
    out[:, k_i] = 1.0 - s
    out[:, k_j] = s
    return out


# --- the partition ------------------------------------------------------------


def _split(rmesh, labels, l, b):
    """Cut cell ``l`` in two by a vertical plane across its long axis."""
    mine = np.nonzero(labels == l)[0]
    P = rmesh.verts[mine]
    c = P.mean(axis=0)
    Q = P - c
    Q -= np.outer(Q @ b, b)
    _s, _v, vt = np.linalg.svd(Q, full_matrices=False)
    m = vt[0] - (vt[0] @ b) * b
    m /= np.linalg.norm(m)
    side = (P - c) @ m < 0.0
    new = int(labels.max()) + 1
    labels[mine[side]] = new
    return new, (c, m)


def _record_split(split_planes: dict, l: int, new: int, plane) -> None:
    """Note that cell ``new`` was cut from ``l`` along ``plane``.

    A cut cell's earlier cuts stay where they were, so each of ``l``'s existing
    planes now also bounds ``new``: without that, a boundary between two pieces
    of a twice-split cell found no plane, fell back to the edge midpoints and
    zig-zagged by 0.4 mm on the chamber fixture.
    """
    for (a, c), pl in list(split_planes.items()):
        if a == l or c == l:
            other = c if a == l else a
            split_planes[(min(other, new), max(other, new))] = pl
    split_planes[(min(l, new), max(l, new))] = plane


def _fan_enters_body(P, T, overshoot, inside_body) -> bool:
    """Would a cell landing on targets ``T`` from ``P`` close its loft through the body?

    The prediction of :func:`.build.loft`'s own check, from the cell's vertex
    rays rather than its perimeter: bottom points ``overshoot`` past the
    targets, an apex past their mean, and a few points along each segment.
    """
    D = T - P
    u = D / np.maximum(np.linalg.norm(D, axis=1, keepdims=True), 1e-12)
    B = T + overshoot * u
    if len(B) > 400:
        B = B[np.linspace(0, len(B) - 1, 400).astype(int)]
    mean = u.mean(axis=0)
    apex = B.mean(axis=0) + overshoot * mean / max(np.linalg.norm(mean), 1e-12)
    probe = np.vstack([B + f * (apex - B) for f in (0.25, 0.5, 0.75)])
    return bool(inside_body(probe).any())


class _Features:
    """Points along the convex edges of the body's anchors, and distance to them.

    Rays converge on such an edge: from everywhere above a cup's mouth, the
    nearest point of the ledge round it at the first drop that reaches it is the
    rim. Under the mouth what the support needs is the cone of those rays, which
    no loft from a cell's perimeter is; so a cell never mixes rays landing on a
    rim with rays that follow their sources (docs/algorithm.md §14.3 S3b), and
    the first kind are closed by that cone (§14.3 S4).

    The edges are read off the B-rep (``obstacles.convex_edge_points``), convex
    seen from outside the body. Lattice anchors are not looked at yet.
    """

    def __init__(self, points: np.ndarray):
        pts = np.asarray(points, dtype=float).reshape(-1, 3)
        self.points = pts
        self.cell = FOLD_TOL * 2.0
        self.grid = defaultdict(list)
        for k, g in enumerate(np.floor(pts / self.cell).astype(np.int64)):
            self.grid[tuple(g)].append(k)

    def distance(self, X: np.ndarray, cap: float = 4.0 * FOLD_TOL) -> np.ndarray:
        """Distance from each point to the nearest feature sample, capped."""
        out = np.full(len(X), cap)
        if not len(self.points):
            return out
        reach = int(np.ceil(cap / self.cell))
        for i, x in enumerate(np.asarray(X)):
            if not np.isfinite(x).all():
                continue
            g = np.floor(x / self.cell).astype(np.int64)
            idx = []
            for dx in range(-reach, reach + 1):
                for dy in range(-reach, reach + 1):
                    for dz in range(-reach, reach + 1):
                        idx += self.grid.get((g[0] + dx, g[1] + dy, g[2] + dz), [])
            if idx:
                out[i] = min(cap, float(np.linalg.norm(self.points[idx] - x, axis=1).min()))
        return out


def _features(obstacles):
    got = getattr(obstacles, "_features", None)
    if got is None:
        got = obstacles._features = _Features(getattr(obstacles, "feature_points",
                                                      np.empty((0, 3))))
    return got


def partition(rmesh, region, rays, obstacles, b, thickness: float,
              step: float = STEP, report=None, inside_body=None,
              overshoot: float = 0.0) -> Partition:
    """Cells of one overhang ``region``, from ``rays`` cast at ``region.verts``.

    With ``inside_body`` given, a cell whose loft would close through the body
    -- one landing on a wall curving away from it, like a round chamber's wall
    seen from inside -- is split across its long axis until it would not.
    """
    V = len(rmesh.verts)
    b = np.asarray(b, dtype=float)
    caster = _Caster(obstacles, b)
    tgt = np.full((V, 3), np.nan)
    drop = np.full(V, np.nan)
    length = np.full(V, np.nan)
    key = np.full(V, -1, dtype=np.int64)
    rv = region.verts
    tgt[rv], drop[rv], length[rv], key[rv] = rays.target, rays.drop, rays.length, rays.key
    exempt = np.zeros(V, dtype=bool)
    exempt[rv] = rays.length <= REACH

    tris = rmesh.tris[region.tris]
    e = np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]])
    e = np.unique(np.sort(e, axis=1), axis=0)
    i, j = e[:, 0], e[:, 1]
    live = ~exempt[i] & ~exempt[j]
    dx = np.linalg.norm(rmesh.verts[i] - rmesh.verts[j], axis=1)
    dy = np.linalg.norm(tgt[i] - tgt[j], axis=1)
    join = live & (key[i] == key[j]) & (dy <= KAPPA * dx + 2.0 * step)
    # Rays landing on a convex edge of their anchor and rays following their
    # sources are never one cell (_Features).
    conv = np.zeros(V, dtype=bool)
    feats = None
    if inside_body is not None:
        feats = _features(obstacles)
        conv[rv] = feats.distance(tgt[rv]) < FOLD_TOL
        join &= conv[i] == conv[j]

    labels = np.full(V, OUTSIDE, dtype=np.int64)
    labels[rv] = EXEMPT
    cell_v = rv[~exempt[rv]]
    if len(cell_v) == 0:
        return Partition(labels, [], {"cells": 0})
    local = {int(v): k for k, v in enumerate(cell_v)}
    pairs = np.array([(local[int(a)], local[int(c)]) for a, c in e[join]],
                     dtype=np.int64).reshape(-1, 2)
    comp = _union_find(len(cell_v), pairs)
    _u, lab = np.unique(comp, return_inverse=True)
    labels[cell_v] = lab
    cell_key = {int(l): int(key[cell_v[np.nonzero(lab == l)[0][0]]]) for l in np.unique(lab)}

    a_, b_, c_ = (rmesh.verts[tris[:, k]] for k in range(3))
    tri_area = 0.5 * np.linalg.norm(np.cross(b_ - a_, c_ - a_), axis=1)

    def cell_areas():
        out = defaultdict(float)
        for k in range(3):
            lk = labels[tris[:, k]]
            for l, aa in zip(lk, tri_area / 3.0):
                if l >= 0:
                    out[int(l)] += aa
        return out

    # Fold cells too small to stand on their own into a neighbour.
    folded = 0
    areas = cell_areas()
    for l, aa in sorted(areas.items(), key=lambda kv: kv[1]):
        if aa >= thickness ** 2:
            continue
        mine = np.nonzero(labels == l)[0]
        if not len(mine):
            continue
        nb = defaultdict(int)
        for a, c in e:
            la, lc = labels[a], labels[c]
            if la == l and lc >= 0 and lc != l:
                nb[int(lc)] += 1
            elif lc == l and la >= 0 and la != l:
                nb[int(la)] += 1
        if not nb:
            continue
        other = max(nb, key=nb.get)
        t2, d2, l2 = caster(rmesh.verts[mine], rmesh.normals[mine], cell_key[other])
        if not np.isfinite(d2).all():
            continue                       # the neighbour's anchor is out of reach
        # Refused where the re-cast ray runs through something before its
        # target: then the neighbour's anchor is not what that ray reaches.
        tf, df, lf = caster.full(rmesh.verts[mine], rmesh.normals[mine])
        if np.any(df < d2 - 1e-6) and not _clear(rmesh.verts[mine], t2, caster):
            continue
        labels[mine] = other
        tgt[mine], drop[mine], length[mine] = t2, d2, l2
        key[mine] = cell_key[other]
        folded += 1

    # Split cells that enclose another cell, until none does.
    split_planes = {}
    splits = 0
    for _round in range(SPLIT_ROUNDS_MAX):
        nested = _nested_cells(rmesh, tris, labels)
        if not nested:
            break
        for l in nested:
            new, plane = _split(rmesh, labels, l, b)
            cell_key[new] = cell_key[l]
            _record_split(split_planes, l, new, plane)
            splits += 1
    else:
        nested = _nested_cells(rmesh, tris, labels)
        if nested:
            centre = rmesh.verts[labels == nested[0]].mean(axis=0)
            from ..errors import ProcessingError

            raise ProcessingError(
                f"A support cell near {np.round(centre, 3).tolist()} still "
                f"encloses another after {SPLIT_ROUNDS_MAX} splits "
                f"(docs/algorithm.md §14.3 S3b)."
            )

    if inside_body is not None:
        for _round in range(SPLIT_ROUNDS_MAX):
            bad = []
            for l in np.unique(labels[labels >= 0]):
                mine = np.nonzero(labels == l)[0]
                # Converging cells close upward, by the cone of their rays (§14.3
                # S4), so a fan beyond their anchor says nothing about them.
                if conv[mine].mean() > 0.5:
                    continue
                if len(mine) >= 8 and _fan_enters_body(
                        rmesh.verts[mine], tgt[mine], overshoot, inside_body):
                    bad.append(int(l))
            if not bad:
                break
            for l in bad:
                new, plane = _split(rmesh, labels, l, b)
                cell_key[new] = cell_key[l]
                _record_split(split_planes, l, new, plane)
                splits += 1

    cell_conv = {}
    for l in np.unique(labels[labels >= 0]):
        cell_conv[int(l)] = bool(conv[labels == l].mean() > 0.5)
    cells = _perimeters(rmesh, tris, labels, cell_key, tgt, drop, length, exempt,
                        caster, split_planes, b, region.tris, feats, cell_conv)
    stats = {"cells": len(cells), "cells_folded": folded, "cells_split": splits,
             "exempt_vertices": int(exempt[rv].sum()),
             "converging_cells": int(sum(cell_conv.values()))}
    return Partition(labels, cells, stats)


def _clear(apex, target, caster) -> bool:
    """Do the segments from ``apex`` to ``target`` meet no obstacle on the way?"""
    obs = caster.obstacles
    A, B, C = obs._A, obs._B, obs._C
    for x, y in zip(apex, target):
        d = y - x
        L = float(np.linalg.norm(d))
        if L < 1e-9:
            continue
        d /= L
        lo, hi = np.minimum(x, y), np.maximum(x, y)
        cand = obs.hash.query(lo, hi)
        if not len(cand):
            continue
        e1, e2 = B[cand] - A[cand], C[cand] - A[cand]
        pv = np.cross(d, e2)
        det = np.einsum("ij,ij->i", e1, pv)
        ok = np.abs(det) > 1e-14
        inv = np.where(ok, 1.0 / np.where(ok, det, 1.0), 0.0)
        tv = x - A[cand]
        u = np.einsum("ij,ij->i", tv, pv) * inv
        qv = np.cross(tv, e1)
        v = (qv @ d) * inv
        t = np.einsum("ij,ij->i", e2, qv) * inv
        if np.any(ok & (u >= 0) & (v >= 0) & (u + v <= 1) & (t > 1e-6) & (t < L - 1e-6)):
            return False
    return True


def _segments(tris, labels, cross_id, centre_id):
    """Directed boundary segments per cell, as pairs of point ids.

    Point ids are ``("v", i)`` for a mesh vertex, ``("x", i, j)`` (sorted) for
    the crossing on edge ``(i, j)``, ``("c", t)`` for triangle ``t``'s triple
    point. Segments are oriented with the cell on their left, seen from the
    outward normal.
    """
    from collections import Counter

    edge_use = Counter()
    for a, b_, c in tris:
        for p, q in ((a, b_), (b_, c), (c, a)):
            edge_use[(min(p, q), max(p, q))] += 1
    segs = defaultdict(list)
    for t, (a, b_, c) in enumerate(tris):
        v = (a, b_, c)
        lab = [int(labels[x]) for x in v]

        def x_(p, q):
            return ("x", min(p, q), max(p, q))

        distinct = set(lab)
        if len(distinct) == 2:
            # The odd corner out.
            for k in range(3):
                if lab[k] != lab[(k + 1) % 3] and lab[k] != lab[(k + 2) % 3]:
                    odd = k
                    break
            o, n1, n2 = v[odd], v[(odd + 1) % 3], v[(odd + 2) % 3]
            lo, lm = lab[odd], lab[(odd + 1) % 3]
            # Corners in winding order: o, n1, n2. The majority's part is
            # n1, n2, x(n2,o), x(o,n1); its interior boundary runs
            # x(n2,o) -> x(o,n1). The odd corner's runs the other way.
            if lm >= 0:
                segs[lm].append((x_(n2, o), x_(o, n1)))
            if lo >= 0:
                segs[lo].append((x_(o, n1), x_(n2, o)))
        elif len(distinct) == 3:
            cid = ("c", t)
            for k in range(3):
                p, q, r = v[k], v[(k + 1) % 3], v[(k + 2) % 3]
                lk = lab[k]
                if lk >= 0:
                    # Corner p's part: p, x(p,q), centre, x(r,p).
                    segs[lk].append((x_(p, q), cid))
                    segs[lk].append((cid, x_(r, p)))
        # Edges on the region's own border.
        for k in range(3):
            p, q = v[k], v[(k + 1) % 3]
            if edge_use[(min(p, q), max(p, q))] != 1:
                continue
            lp, lq = lab[k], lab[(k + 1) % 3]
            if lp == lq:
                if lp >= 0:
                    segs[lp].append((("v", p), ("v", q)))
            else:
                if lp >= 0:
                    segs[lp].append((("v", p), x_(p, q)))
                if lq >= 0:
                    segs[lq].append((x_(p, q), ("v", q)))
    return segs


def _chain(segments):
    """Directed segments into closed loops of point ids; ``None`` if they do not close."""
    out_of = defaultdict(list)
    for s, t in segments:
        out_of[s].append(t)
    loops = []
    used = defaultdict(int)
    for start in list(out_of):
        while used[start] < len(out_of[start]):
            loop = [start]
            cur = start
            while True:
                if used[cur] >= len(out_of[cur]):
                    return None
                nxt = out_of[cur][used[cur]]
                used[cur] += 1
                if nxt == start:
                    break
                loop.append(nxt)
                cur = nxt
                if len(loop) > 10 * len(segments) + 10:
                    return None
            loops.append(loop)
    return loops


def _nested_cells(rmesh, tris, labels) -> list:
    """Cells with an inner loop touching another cell (crossings at edge midpoints)."""
    segs = _segments(tris, labels, None, None)
    pos = {}

    def point(pid):
        if pid[0] == "v":
            return rmesh.verts[pid[1]]
        if pid[0] == "x":
            return 0.5 * (rmesh.verts[pid[1]] + rmesh.verts[pid[2]])
        a, b_, c = tris[pid[1]]
        return rmesh.verts[[a, b_, c]].mean(axis=0)

    out = []
    for l, s in segs.items():
        loops = _chain(s)
        if loops is None or len(loops) < 2:
            continue
        areas = []
        for loop in loops:
            P = np.array([pos.setdefault(p, point(p)) for p in loop])
            n = rmesh.normals[labels == l].mean(axis=0)
            areas.append(0.5 * float(np.cross(P, np.roll(P, -1, axis=0)).sum(axis=0) @ n))
        outer = int(np.argmax(np.abs(areas)))
        for k, loop in enumerate(loops):
            if k == outer:
                continue
            for p in loop:
                if p[0] == "x" and labels[p[1]] >= 0 and labels[p[2]] >= 0 \
                        and labels[p[1]] != labels[p[2]]:
                    out.append(l)
                    break
            if out and out[-1] == l:
                break
    return out


def _perimeters(rmesh, tris, labels, cell_key, tgt, drop, length, exempt, caster,
                split_planes, b, tri_global, feats=None, cell_conv=None):
    """Exact crossings, triple points and rays; then each cell's loops."""
    etri = _edge_triangles(tris)
    # Crossing edges, grouped by kind.
    xs = {}
    parting_edges, exempt_edges, cut_edges, fold_edges = [], [], [], []
    for (p, q), (t, kp, kq) in etri.items():
        lp, lq = int(labels[p]), int(labels[q])
        if lp == lq or (lp < 0 and lq < 0):
            continue
        if lp >= 0 and lq >= 0:
            pair = (min(lp, lq), max(lp, lq))
            if (cell_conv is not None and pair not in split_planes
                    and cell_key[lp] == cell_key[lq] and cell_conv[lp] != cell_conv[lq]):
                fold_edges.append((p, q, t, kp, kq))
            elif pair in split_planes or cell_key[lp] == cell_key[lq]:
                cut_edges.append((p, q, t, kp, kq))
            else:
                parting_edges.append((p, q, t, kp, kq))
        else:
            exempt_edges.append((p, q, t, kp, kq))

    def surface(edges, s):
        # Triangle indices here are the region's own; the mesh is the whole
        # body's. They agreed while a part had one region, and stopped agreeing
        # on Hole-block-chamfered, whose crease splits it in two: parting points
        # were evaluated on the wrong facet, 40 mm away.
        tt = tri_global[np.array([x[2] for x in edges], dtype=np.int64)]
        bary = np.stack([_bary(x[3], x[4], np.array([si]))[0] for x, si in zip(edges, s)])
        return rmesh.surface_points(tt, bary)

    def cast_groups(points, normals, keys):
        out_t = np.zeros((len(points), 3))
        out_d = np.zeros(len(points))
        out_l = np.zeros(len(points))
        for k in np.unique(keys):
            sel = keys == k
            t_, d_, l_ = caster(points[sel], normals[sel], int(k))
            out_t[sel], out_d[sel], out_l[sel] = t_, d_, l_
        return out_t, out_d, out_l

    # Partings: the zero of drop_A - drop_B along the edge, by bisection.
    if parting_edges:
        P_ = np.array([x[0] for x in parting_edges])
        Q_ = np.array([x[1] for x in parting_edges])
        kA = np.array([cell_key[int(labels[p])] for p in P_])
        kB = np.array([cell_key[int(labels[q])] for q in Q_])
        lo = np.zeros(len(P_))
        hi = np.ones(len(P_))
        L = np.linalg.norm(rmesh.verts[P_] - rmesh.verts[Q_], axis=1)
        iters = int(np.ceil(np.log2(max(L.max(), PARTING_TOL) / PARTING_TOL))) + 1
        for _ in range(iters):
            mid = 0.5 * (lo + hi)
            pts, nrm = surface(parting_edges, mid)
            _ta, da, _la = cast_groups(pts, nrm, kA)
            _tb, db, _lb = cast_groups(pts, nrm, kB)
            # At p (A's own vertex) A is reached first: f = dA - dB <= 0 there.
            # Where neither anchor is in reach both drops are infinite; that is
            # A's side, as at p itself, rather than a NaN steering the search.
            left = ~(da > db)
            lo = np.where(left, mid, lo)
            hi = np.where(left, hi, mid)
        s = 0.5 * (lo + hi)
        pts, nrm = surface(parting_edges, s)
        ta, _da, _ = cast_groups(pts, nrm, kA)
        tb, _db, _ = cast_groups(pts, nrm, kB)
        for k, x in enumerate(parting_edges):
            p, q = x[0], x[1]
            lp, lq = int(labels[p]), int(labels[q])
            xs[(min(p, q), max(p, q))] = (pts[k], nrm[k], {lp: ta[k], lq: tb[k]}, PARTING)

    # Exempt borders: where the unrestricted ray length crosses the reach.
    if exempt_edges:
        s = []
        keys = []
        for p, q, *_r in exempt_edges:
            lp = length[p] if labels[p] >= 0 else length[p]
            Lp, Lq = length[p], length[q]
            # Measured from p: p in a cell (long), q exempt (short), or reverse.
            frac = (Lp - REACH) / (Lp - Lq) if abs(Lp - Lq) > 1e-12 else 0.5
            s.append(float(np.clip(frac, 0.0, 1.0)))
            cl = int(labels[p]) if labels[p] >= 0 else int(labels[q])
            keys.append(cell_key[cl])
        pts, nrm = surface(exempt_edges, np.array(s))
        t_, _d, _l = cast_groups(pts, nrm, np.array(keys))
        for k, x in enumerate(exempt_edges):
            p, q = x[0], x[1]
            cl = int(labels[p]) if labels[p] >= 0 else int(labels[q])
            xs[(min(p, q), max(p, q))] = (pts[k], nrm[k], {cl: t_[k]}, EXEMPT_EDGE)

    # Cuts inside one anchor group: where a split plane crosses, or midway.
    if cut_edges:
        s = []
        keys = []
        for p, q, *_r in cut_edges:
            lp, lq = int(labels[p]), int(labels[q])
            plane = split_planes.get((min(lp, lq), max(lp, lq)))
            if plane is None:
                s.append(0.5)
            else:
                c, m = plane
                fp = float((rmesh.verts[p] - c) @ m)
                fq = float((rmesh.verts[q] - c) @ m)
                s.append(float(np.clip(fp / (fp - fq), 0.0, 1.0)) if fp != fq else 0.5)
            keys.append(cell_key[lp])
        pts, nrm = surface(cut_edges, np.array(s))
        t_, _d, _l = cast_groups(pts, nrm, np.array(keys))
        for k, x in enumerate(cut_edges):
            p, q = x[0], x[1]
            lp, lq = int(labels[p]), int(labels[q])
            xs[(min(p, q), max(p, q))] = (pts[k], nrm[k], {lp: t_[k], lq: t_[k]}, CUT)

    # Folds: where a target is FOLD_TOL from a convex edge, by bisection.
    if fold_edges:
        keys = np.array([cell_key[int(labels[x[0]])] for x in fold_edges])
        conv_p = np.array([cell_conv[int(labels[x[0]])] for x in fold_edges])
        lo = np.zeros(len(fold_edges))
        hi = np.ones(len(fold_edges))
        L = max(float(np.linalg.norm(rmesh.verts[x[0]] - rmesh.verts[x[1]]))
                for x in fold_edges)
        for _ in range(int(np.ceil(np.log2(max(L, PARTING_TOL) / PARTING_TOL))) + 1):
            mid = 0.5 * (lo + hi)
            pts, nrm = surface(fold_edges, mid)
            t_, _d, _l = cast_groups(pts, nrm, keys)
            near = feats.distance(t_) < FOLD_TOL
            left = near == conv_p
            lo = np.where(left, mid, lo)
            hi = np.where(left, hi, mid)
        pts, nrm = surface(fold_edges, 0.5 * (lo + hi))
        t_, _d, _l = cast_groups(pts, nrm, keys)
        for k, x in enumerate(fold_edges):
            p, q = x[0], x[1]
            lp, lq = int(labels[p]), int(labels[q])
            xs[(min(p, q), max(p, q))] = (pts[k], nrm[k], {lp: t_[k], lq: t_[k]}, FOLD)

    # Triple points: the mean of the triangle's three crossings, cast per cell.
    centres = {}
    for t, (a, b_, c) in enumerate(tris):
        lab = {int(labels[a]), int(labels[b_]), int(labels[c])}
        if len(lab) != 3:
            continue
        ps = [xs[(min(p, q), max(p, q))] for p, q in ((a, b_), (b_, c), (c, a))]
        P = np.mean([x[0] for x in ps], axis=0)
        bary_guess = np.full((1, 3), 1.0 / 3.0)
        _pp, nn = rmesh.surface_points(np.array([tri_global[t]]), bary_guess)
        targets = {}
        for l in lab:
            if l >= 0:
                tt, _d, _l = caster(P[None], nn, cell_key[l])
                targets[l] = tt[0]
        centres[t] = (P, nn[0], targets)

    def own_ray(l, P):
        """The ray of cell ``l``'s vertex nearest ``P``: for a perimeter point
        from which a cast restricted to the cell's anchor found nothing."""
        mine = np.nonzero(labels == l)[0]
        k = mine[np.argmin(np.linalg.norm(rmesh.verts[mine] - P, axis=1))]
        return P + (tgt[k] - rmesh.verts[k])

    segs = _segments(tris, labels, None, None)
    cells = []
    for l in sorted(segs):
        loops_ids = _chain(segs[l])
        mine = np.nonzero(labels == l)[0]
        if loops_ids is None or not len(mine):
            continue
        loops = []
        divergence = 180.0
        for ids in loops_ids:
            pts, nrm, tg, kind, keys = [], [], [], [], []
            for pid in ids:
                if pid[0] == "v":
                    v = pid[1]
                    pts.append(rmesh.verts[v])
                    nrm.append(rmesh.normals[v])
                    tg.append(tgt[v])
                    kind.append(VERTEX)
                elif pid[0] == "x":
                    P, N, tmap, kd = xs[(pid[1], pid[2])]
                    pts.append(P)
                    nrm.append(N)
                    tg.append(tmap[l] if np.isfinite(tmap[l]).all() else own_ray(l, P))
                    kind.append(kd)
                    if kd == PARTING and all(np.isfinite(v_).all() for v_ in tmap.values()):
                        other = [v_ for k_, v_ in tmap.items() if k_ != l][0]
                        u1, u2 = tmap[l] - P, other - P
                        ang = np.degrees(np.arccos(np.clip(
                            u1 @ u2 / max(np.linalg.norm(u1) * np.linalg.norm(u2), 1e-300),
                            -1.0, 1.0)))
                        divergence = min(divergence, float(ang))
                else:
                    P, N, tmap = centres[pid[1]]
                    pts.append(P)
                    nrm.append(N)
                    tg.append(tmap[l] if np.isfinite(tmap[l]).all() else own_ray(l, P))
                    kind.append(TRIPLE)
                keys.append(pid)
            loops.append(Loop(np.array(pts), np.array(nrm), np.array(tg),
                              np.array(kind, dtype=np.int64), keys))
        n_avg = rmesh.normals[mine].mean(axis=0)
        n_avg /= np.linalg.norm(n_avg)
        signed = [0.5 * float(np.cross(lp_.points, np.roll(lp_.points, -1, axis=0)).sum(axis=0)
                              @ n_avg) for lp_ in loops]
        order = np.argsort(-np.abs(signed))
        loops = [loops[k] for k in order]
        # Seed: the cell vertex furthest from its own perimeter, a hair inside.
        perim = np.vstack([lp_.points for lp_ in loops])
        far = mine[np.argmax(np.min(np.linalg.norm(
            rmesh.verts[mine][:, None, :] - perim[None, :, :], axis=2), axis=1))] \
            if len(perim) * len(mine) < 4e7 else mine[0]
        seed = rmesh.verts[far] - 2.0 * APEX_NUDGE * rmesh.normals[far]
        area = float(abs(signed[order[0]]) - sum(abs(signed[k]) for k in order[1:]))
        cells.append(Cell(int(l), int(cell_key[l]), mine, area, loops, seed, divergence))
    return cells
