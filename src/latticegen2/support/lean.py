"""Perimeter rays that lean into their cell (docs/algorithm.md §14.3 S3d). Pure NumPy.

A ray landing on a level face of its anchor is one of a whole disc of targets
at the same drop, and the first-contact tie-break picks the one straight
below. Only the rays on a cell's perimeter shape its solid, and there that
makes the side face vertical: on the cup fixture, a cylinder from the ceiling
down to the ledge, carrying a block of support the ceiling does not need.

Instead the target moves along the perimeter's horizontal inward normal ``m``,
``y = y0 + s m``, as far as the cone, the level face, a clear line of sight and
the cell's **footing** allow: the footprint left on the anchor keeps at least
``t_s`` across. A cut is an inner face of one support split in two, so its
rays lean once for both cells, interpolated between its two ends; fold lines
keep the one ray their two cells share, and crease points
(:data:`.cells.CREASE`) have been moved into the wall already.
"""

from __future__ import annotations

from collections import defaultdict
from types import SimpleNamespace

import numpy as np

from . import cells as cells_
from .rays import APEX_NUDGE, TAN_T

LEVEL_TOL = 1e-9
"""A facet whose normal is within this of ``+-b`` (in ``1 - |n . b|``) is level."""

HEIGHT_TOL = 1e-6
"""mm. A level facet this close to a target's height is at the same height."""

HALVINGS = 6
"""Times a cell's leans are halved together before a folding footprint gives up
leaning altogether."""

CLEAR_HALVINGS = 3
"""Times one point's lean is halved to clear an obstacle before it stays put."""

TAPER = 1.0
"""How fast a lean may change along the perimeter, in mm per mm. Away from a
point that does not lean -- a triple point keeps its ray while the parting
beside it leans -- without a taper the band twisted from one to the other
within a span; and between two leaning points, footings measured to different
far sides moved neighbours past each other."""

TANGENT_ARC = 1.0
"""mm. The inward normal is taken across this much perimeter either side. A
parting's points zig-zag by tens of microns, so a neighbour-to-neighbour
tangent swung ``m`` by 10 deg or more, and a 2 mm lean turned that into
neighbours swapping places: the leaned footprint crossed itself on the cup."""

LEANING = (cells_.PARTING, cells_.VERTEX)
"""Point kinds whose rays lean. A triple point leans when both its neighbours do."""


def lean_cells(cells: list, obstacles, b, thickness: float, overshoot: float) -> dict:
    """Lean every cell's outer-loop rays in place; return summary stats.

    Cells are leaned together, not one at a time: a cut is an inner face of
    what would have been one support, so the cells on either side of it carry
    the same rays along it -- leaned once for both, from its two ends
    (:func:`_cut_leans`) -- and a lean one of them must give up, the other
    gives up too (:func:`_sync`).
    """
    b = np.asarray(b, dtype=float)
    b = b / np.linalg.norm(b)
    e1, e2 = _frame(b)
    to2 = np.column_stack([e1, e2])
    ws = [_prepare(c.loops[0], c.key, obstacles, b, to2, thickness, overshoot) for c in cells]
    shared = defaultdict(list)
    for i, (c, w) in enumerate(zip(cells, ws)):
        loop = c.loops[0]
        for k in np.nonzero(np.isin(loop.kind, (cells_.CUT, cells_.TRIPLE)))[0]:
            shared[(loop.keys[k], c.key)].append((i, int(k)))
    shared = [v for v in shared.values() if len(v) > 1]
    # Shared points lean along one direction.
    for group in shared:
        m = sum(ws[i].m[k] for i, k in group)
        L = float(np.linalg.norm(m))
        for i, k in group:
            ws[i].m[k] = m / L if L > 1e-12 else 0.0
    # Families: the cells one support was split into.
    family = list(range(len(ws)))

    def root(i):
        while family[i] != i:
            i = family[i]
        return i

    for group in shared:
        for i, _k in group[1:]:
            family[root(i)] = root(group[0][0])
    members = defaultdict(list)
    for i, w in enumerate(ws):
        if not w.frozen:
            members[root(i)].append(i)
    for fam in members.values():
        if len(fam) > 1:
            _union_normals([ws[i] for i in fam], b)
        _footings([ws[i] for i in fam], thickness, to2)
    _sync(ws, shared)
    for w in ws:
        _cut_leans(w, obstacles, b, to2)
    _sync(ws, shared)
    for _ in range(4 * HALVINGS):
        changed = False
        for w in ws:
            new = _taper(w.P, w.s)
            changed |= bool(np.any(new < w.s - 1e-12))
            w.s = new
        changed |= _sync(ws, shared)
        if not changed:
            break
    caster = SimpleNamespace(obstacles=obstacles)
    for w in ws:
        for k in np.nonzero(w.s)[0]:
            for _ in range(CLEAR_HALVINGS + 1):
                if cells_._clear(w.apex[k][None], (w.T[k] + w.s[k] * w.m[k])[None], caster):
                    break
                w.s[k] *= 0.5
            else:
                w.s[k] = 0.0
    _sync(ws, shared)
    # A footprint that folds halves the lean of every cell in its family.
    for fam in members.values():
        for _ in range(HALVINGS + 1):
            if not any(_crosses((ws[i].T + ws[i].s[:, None] * ws[i].m) @ to2) for i in fam):
                break
            for i in fam:
                ws[i].s *= 0.5
        else:
            for i in fam:
                ws[i].s[:] = 0.0
    for c, w in zip(cells, ws):
        c.loops[0].targets = w.T + w.s[:, None] * w.m
    # One value per shared point, bit for bit: the two cells interpolate a cut
    # from opposite ends and agree only to the last ulp otherwise.
    for group in shared:
        i0, k0 = group[0]
        for i, k in group[1:]:
            cells[i].loops[0].targets[k] = cells[i0].loops[0].targets[k0]
    n_points = 0
    worst = 0.0
    for c, w in zip(cells, ws):
        loop = c.loops[0]
        leaned = w.s > 1e-9
        if not leaned.any():
            continue
        d = loop.targets[leaned] - w.apex[leaned]
        ang = np.degrees(np.arctan2(np.linalg.norm(_horizontal(d, b), axis=1), -(d @ b)))
        n_points += int(leaned.sum())
        worst = max(worst, float(ang.max()))
    return {"leaned_points": n_points, "lean_max_deg": round(worst, 2)}


def _sync(ws: list, shared: list) -> bool:
    """Every shared point to the smallest of its leans; whether any changed."""
    changed = False
    for group in shared:
        lo = min(ws[i].s[k] for i, k in group)
        for i, k in group:
            if ws[i].s[k] > lo + 1e-12:
                ws[i].s[k] = lo
                changed = True
    return changed


def _horizontal(v: np.ndarray, b: np.ndarray) -> np.ndarray:
    return v - np.outer(v @ b, b) if v.ndim == 2 else v - (v @ b) * b


def _frame(b: np.ndarray):
    e1 = np.cross(b, [1.0, 0.0, 0.0])
    if np.linalg.norm(e1) < 1e-6:
        e1 = np.cross(b, [0.0, 1.0, 0.0])
    e1 /= np.linalg.norm(e1)
    return e1, np.cross(b, e1)


def _which(loop) -> np.ndarray:
    """Which of the loop's points lean."""
    kind, pushed = loop.kind, loop.pushed
    ok = np.isin(kind, LEANING) & ~pushed
    n = len(kind)
    for k in np.nonzero(kind == cells_.TRIPLE)[0]:
        a, c = (k - 1) % n, (k + 1) % n
        cut_a, cut_c = kind[a] == cells_.CUT, kind[c] == cells_.CUT
        if cut_a != cut_c:
            # A cut's end leans with the parting it ends on (_cut_leans).
            ok[k] = bool(ok[a] if cut_c else ok[c]) and not pushed[k]
        else:
            ok[k] = bool(ok[a] and ok[c])
    return ok


def _inward(P: np.ndarray, b: np.ndarray, arc: float = TANGENT_ARC,
            stop: np.ndarray | None = None) -> np.ndarray:
    """Unit horizontal inward normal at every loop point, over ``arc`` either side.

    The window never runs onto a ``stop`` point: a cut's corner is the cut's,
    not the support's, so beside a cut the normal is the parting's alone. Taken
    across the corner it swung up to 45 deg toward the cut, the lean then
    carried the point across the cut line, and the cell's footprint folded.
    """
    area = 0.5 * float(np.cross(P, np.roll(P, -1, axis=0)).sum(axis=0) @ b)
    n = len(P)
    seg = np.linalg.norm(np.roll(P, -1, axis=0) - P, axis=1)
    total = float(seg.sum())
    s = np.r_[0.0, np.cumsum(seg)[:-1]]
    ss = np.r_[s - total, s, s + total]          # three turns, for the wrap
    idx = np.r_[np.arange(n), np.arange(n), np.arange(n)]
    w = min(arc, total / 4.0)
    ahead = idx[np.clip(np.searchsorted(ss, s + w), 0, 3 * n - 1)]
    behind = idx[np.clip(np.searchsorted(ss, s - w, side="right") - 1, 0, 3 * n - 1)]
    if stop is not None and stop.any():
        for k in np.nonzero(~stop)[0]:
            for step, out in ((1, ahead), (-1, behind)):
                j, run = k, 0.0
                while run < w:
                    jn = (j + step) % n
                    if stop[jn] or jn == k:
                        break
                    run += seg[j] if step == 1 else seg[jn]
                    j = jn
                out[k] = j
    T = P[ahead] - P[behind]
    m = np.cross(b, T) if area > 0.0 else np.cross(T, b)
    m = _horizontal(m, b)
    L = np.linalg.norm(m, axis=1, keepdims=True)
    return np.where(L > 1e-12, m / np.maximum(L, 1e-300), 0.0)


def _level_facets(obstacles, b: np.ndarray, key: int):
    """Indices of the obstacle triangles that are level, upward and of ``key``."""
    cache = getattr(obstacles, "_level", None)
    if cache is None:
        cache = obstacles._level = {}
    got = cache.get(key)
    if got is None:
        n = obstacles._n
        level = (1.0 - np.abs(n @ b) <= LEVEL_TOL) & (obstacles.key == key)
        got = cache[key] = level
    return got


def _on_level(y0: np.ndarray, m: np.ndarray, s_max: float, obstacles, b, level,
              e1, e2) -> float:
    """How far ``y0 + s m`` stays on level facets at ``y0``'s height, from ``s = 0``.

    The union of every facet's interval along the line, merged from 0; 0 when
    ``y0`` itself is on none.
    """
    y1 = y0 + s_max * m
    lo = np.minimum(y0, y1) - 1e-6
    hi = np.maximum(y0, y1) + 1e-6
    cand = obstacles.hash.query(lo, hi)
    if not len(cand):
        return 0.0
    cand = cand[level[cand]]
    A, B, C = obstacles._A[cand], obstacles._B[cand], obstacles._C[cand]
    same = np.abs((A - y0) @ b) <= HEIGHT_TOL
    A, B, C = A[same], B[same], C[same]
    if not len(A):
        return 0.0
    to2 = np.column_stack([e1, e2])

    def q(X):
        return X @ to2

    a2, b2, c2 = q(A), q(B), q(C)
    p0, d = q(y0[None])[0], q(m[None])[0]
    t0 = np.zeros(len(A))
    t1 = np.full(len(A), s_max)
    for (U, V, W) in ((a2, b2, c2), (b2, c2, a2), (c2, a2, b2)):
        e = V - U
        side = np.sign(e[:, 0] * (W - U)[:, 1] - e[:, 1] * (W - U)[:, 0])
        # f(t) = side * cross(e, p0 + t d - U) >= 0, linear in t.
        f0 = side * (e[:, 0] * (p0 - U)[:, 1] - e[:, 1] * (p0 - U)[:, 0])
        f1 = side * (e[:, 0] * d[1] - e[:, 1] * d[0])
        eps = 1e-12
        grow = f1 > eps
        shrink = f1 < -eps
        flat = ~grow & ~shrink
        t0 = np.where(grow, np.maximum(t0, -f0 / np.where(grow, f1, 1.0)), t0)
        t1 = np.where(shrink, np.minimum(t1, -f0 / np.where(shrink, f1, 1.0)), t1)
        t1 = np.where(flat & (f0 < -eps), -1.0, t1)
    ok = t1 >= t0 - 1e-12
    if not ok.any():
        return 0.0
    iv = sorted(zip(t0[ok], t1[ok]))
    end = None
    for a, c in iv:
        if end is None:
            if a > 1e-9:
                return 0.0
            end = c
        elif a <= end + 1e-9:
            end = max(end, c)
        else:
            break
    return max(0.0, float(end))


def _cone(x: np.ndarray, y0: np.ndarray, m: np.ndarray, b: np.ndarray) -> float:
    """Largest ``s`` with ``y0 + s m`` inside the cone from ``x`` at ``y0``'s drop."""
    drop = float((x - y0) @ b)
    if drop <= 0.0:
        return 0.0
    h0 = _horizontal(y0 - x, b)
    R = drop * TAN_T
    hm = float(h0 @ m)
    disc = hm * hm - float(h0 @ h0) + R * R
    if disc < 0.0:
        return 0.0
    return max(0.0, -hm + float(np.sqrt(disc)))


def _ray_hits(p, d, Q, R_):
    """``t`` along ray ``p + t d`` where it crosses each segment ``Q -> R_`` (2D), or inf."""
    e = R_ - Q
    den = d[0] * e[:, 1] - d[1] * e[:, 0]
    w = Q - p
    ok = np.abs(den) > 1e-15
    safe = np.where(ok, den, 1.0)
    t = (w[:, 0] * e[:, 1] - w[:, 1] * e[:, 0]) / safe
    u = (w[:, 0] * d[1] - w[:, 1] * d[0]) / safe
    hit = ok & (u >= -1e-12) & (u <= 1.0 + 1e-12) & (t > 1e-9)
    return np.where(hit, t, np.inf)


def _crosses(Q: np.ndarray) -> bool:
    """Does the closed 2D polygon ``Q`` cross itself?"""
    n = len(Q)
    A, B = Q, np.roll(Q, -1, axis=0)

    def orient(p, q, r):
        return (q[..., 0] - p[..., 0]) * (r[..., 1] - p[..., 1]) - \
            (q[..., 1] - p[..., 1]) * (r[..., 0] - p[..., 0])

    for i in range(n):
        j = np.arange(i + 2, n)
        if i == 0:
            j = j[j != n - 1]
        if not len(j):
            continue
        d1 = orient(A[i], B[i], A[j])
        d2 = orient(A[i], B[i], B[j])
        d3 = orient(A[j], B[j], A[i][None])
        d4 = orient(A[j], B[j], B[i][None])
        if np.any((d1 * d2 < 0) & (d3 * d4 < 0)):
            return True
    return False


def _taper(P: np.ndarray, s: np.ndarray, slope: float = TAPER) -> np.ndarray:
    """The largest ``s' <= s`` changing by at most ``slope`` per mm along the loop."""
    n = len(P)
    seg = np.linalg.norm(np.roll(P, -1, axis=0) - P, axis=1)
    out = s.copy()
    for _ in range(2):                           # twice round covers the wrap
        for k in range(n):                       # forward
            j = (k + 1) % n
            out[j] = min(out[j], out[k] + slope * seg[k])
        for k in range(n - 1, -1, -1):           # backward
            j = (k - 1) % n
            out[j] = min(out[j], out[k] + slope * seg[j])
    return out


def _prepare(loop, key: int, obstacles, b, to2, thickness: float, overshoot: float):
    """A loop's leaning state: every point's direction ``m`` and its lean ``s``
    as far as the cone, the level face and the footing allow.

    ``frozen`` loops lean nowhere -- too few points, a ray that found nothing, a
    footprint folded before any lean -- but still take part in sharing, so a
    cell across a cut from one does not lean alone.
    """
    P, T = loop.points, np.array(loop.targets, dtype=float)
    n = len(P)
    w = SimpleNamespace(P=P, T=T, m=np.zeros((n, 3)), s=np.zeros(n),
                        apex=P - APEX_NUDGE * loop.normals, kind=loop.kind,
                        keys=loop.keys, frozen=True, key=key)
    if n < 3 or not np.isfinite(T).all():
        return w
    F = T @ to2                                  # the footprint, before leaning
    if _crosses(F):
        return w
    w.frozen = False
    ok = _which(loop)
    kind = loop.kind
    m = _inward(P, b, stop=kind == cells_.CUT)
    w.m, w.ok, w.F = m, ok, F
    level = _level_facets(obstacles, b, key)
    e1, e2 = to2[:, 0], to2[:, 1]
    w.cap = np.zeros(n)
    for k in np.nonzero(ok)[0]:
        if not np.any(m[k]):
            continue
        s_cone = _cone(w.apex[k], T[k], m[k], b)
        if s_cone <= 0.0:
            continue
        s_face = _on_level(T[k], m[k], s_cone, obstacles, b, level, e1, e2) - overshoot
        w.cap[k] = max(0.0, min(s_cone, s_face))
    w.level, w.overshoot = level, overshoot
    return w


def _union_normals(fam: list, b) -> None:
    """Inward normals of a family's leaning points, off the support's own outline.

    The cells' loops less their cuts, and less every edge two of them run in
    opposite directions, chain into the perimeter of the support they were
    split from. Taken off each cell's loop instead, a cut's end got the mean of
    two one-sided tangents and the point beside it one of them, 5 deg apart on
    the cup, and a 2 mm lean put them in each other's place.
    """
    pos, where, edges = {}, defaultdict(list), set()
    for i, w in enumerate(fam):
        n = len(w.P)
        for k in range(n):
            pos[w.keys[k]] = w.P[k]
            where[w.keys[k]].append((i, k))
            k1 = (k + 1) % n
            if w.kind[k] != cells_.CUT and w.kind[k1] != cells_.CUT:
                edges.add((w.keys[k], w.keys[k1]))
    succ = {}
    for a, c in edges:
        if (c, a) in edges:
            continue
        succ[a] = None if a in succ else c       # None: more than one way on
    done = set()
    for start in list(succ):
        if start in done or succ[start] is None:
            continue
        loop, cur = [], start
        while cur is not None and cur not in done:
            done.add(cur)
            loop.append(cur)
            cur = succ.get(cur)
        if cur != start or len(loop) < 3:
            continue
        m = _inward(np.array([pos[p] for p in loop]), b)
        for p, mp in zip(loop, m):
            for i, k in where[p]:
                fam[i].m[k] = mp


def _footings(fam: list, thickness: float, to2) -> None:
    """Each leaning point's lean, held to its footing across its whole family.

    The footing is the distance along ``m`` to the far side of the footprint,
    less ``t_s``; half of it where the far side leans too. A cut is no far
    side -- the support carries on across it into the other cell -- so the far
    side is sought among every cell of the family and none of their cuts: at a
    cut's end ``m`` runs along the cut itself, and the far side lies in the
    neighbouring cell.
    """
    Q, R_, opp, owner = [], [], [], []
    for i, w in enumerate(fam):
        kind = w.kind
        keep = ~((kind == cells_.CUT) | (np.roll(kind, -1) == cells_.CUT))
        idx = np.nonzero(keep)[0]
        Q.append(w.F[idx])
        R_.append(np.roll(w.F, -1, axis=0)[idx])
        opp.append((w.ok & np.roll(w.ok, -1))[idx])
        owner.append(np.c_[np.full(len(idx), i), idx])
    Q, R_, opp, owner = (np.concatenate(x) for x in (Q, R_, opp, owner))
    for i, w in enumerate(fam):
        n = len(w.P)
        for k in np.nonzero(w.cap > 0.0)[0]:
            t = _ray_hits(w.F[k], w.m[k] @ to2, Q, R_)
            own = owner[:, 0] == i
            t[own & np.isin(owner[:, 1], [(k - 1) % n, k])] = np.inf   # the edges at k
            j = int(np.argmin(t))
            if not np.isfinite(t[j]):
                continue
            share = 2.0 if opp[j] else 1.0
            w.s[k] = max(0.0, min(w.cap[k], (float(t[j]) - thickness) / share))


def _cut_leans(w, obstacles, b, to2) -> None:
    """Leans along every cut of ``w``'s loop, interpolated from its two ends.

    A cut is an inner face of one support split in two, so its rays are not a
    perimeter's: each is the lean of the cut's ends, interpolated along the
    cut by arc length, then held to the cone and the level face as any lean
    is. The two cells sharing the cut compute the same values, from the same
    points and ends.
    """
    if w.frozen:
        return
    kind = w.kind
    n = len(kind)
    is_cut = kind == cells_.CUT
    if not is_cut.any() or is_cut.all():
        return
    e1, e2 = to2[:, 0], to2[:, 1]
    start = int(np.nonzero(~is_cut)[0][0])
    k = start
    for _ in range(n):
        k1 = (k + 1) % n
        if is_cut[k1] and not is_cut[k]:
            run = []
            j = k1
            while is_cut[j]:
                run.append(j)
                j = (j + 1) % n
            a, z = k, j
            chain = [a] + run + [z]
            d = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(w.P[chain], axis=0), axis=1))]
            u = d / d[-1] if d[-1] > 0 else np.linspace(0.0, 1.0, len(chain))
            Da, Dz = w.s[a] * w.m[a], w.s[z] * w.m[z]
            for jj, uu in zip(run, u[1:-1]):
                D = (1.0 - uu) * Da + uu * Dz
                L = float(np.linalg.norm(D))
                if L <= 1e-12:
                    continue
                m = D / L
                s_cone = _cone(w.apex[jj], w.T[jj], m, b)
                s_face = _on_level(w.T[jj], m, max(s_cone, 0.0), obstacles, b, w.level,
                                   e1, e2) - w.overshoot
                w.m[jj] = m
                w.s[jj] = max(0.0, min(L, s_cone, s_face))
        k = k1
        if k == start:
            break
