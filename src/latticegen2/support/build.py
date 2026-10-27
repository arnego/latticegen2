"""Cell solids (docs/algorithm.md §14.3 S4).

A cell's perimeter loop is cut into **runs** at its kinks -- triple points, and
wherever the kind of boundary changes -- and every run becomes one ruled
B-spline face,

    S(s, v) = X(s) + v E(s),        v in [v_top, 1]

with ``X`` the run's curve on the ceiling and ``E`` its rays, lengthened by the
overshoot, in **the same basis**. Every control vector of ``E`` is clamped into
the closed cone of the rod angle, and every ruling ``E(s)`` is a positive
combination of control vectors, so it lies in that convex cone too: **each side
face is self-supporting by construction**, at the rod angle or steeper, and a
cell has a handful of faces however long its perimeter is.

A parting run is fitted once, in a canonical direction, so the two cells on
either side build on one and the same ``X``: their faces leave the ceiling from
the same curve and diverge, sharing an edge rather than two copies of one.

The ruled band is carried on at both ends -- straight up into the shell above
the ceiling, along the cell's mean ray into the anchor -- and closed there by
planar faces (:func:`loft`). The intersection with the body removes the first,
and the second lies in the anchor's own material, outside the body or inside
the lattice. What remains is the support, cut to the ceiling, the walls and the
struts along their own surfaces.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut
from OCP.BRepBuilderAPI import (BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace,
                                BRepBuilderAPI_MakeSolid, BRepBuilderAPI_MakeWire,
                                BRepBuilderAPI_Sewing)
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.Geom import Geom_BSplineCurve, Geom_BSplineSurface
from OCP.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
from OCP.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCP.TopAbs import TopAbs_State
from OCP.TopTools import TopTools_ListOfShape
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Pnt

from .. import occ
from ..errors import ProcessingError
from . import cells as cells_
from .params import STEP, SUPPORT_OVERHANG_DEG

LIFT = 0.2
"""mm. How far up its own rays a cell's band reaches above the ceiling, at its
longest ray: enough for the intersection with the body to cut the band cleanly,
never grazing."""

FIT_TOL = 2e-3
"""mm. How closely a run's ceiling curve follows its perimeter points."""

RAY_FIT_TOL = 0.3
"""mm. How closely a run's ray field follows each cell's rays along it --
loosely, by design. Targets are first contacts on the faceted obstacle mesh and
snap from facet to facet: round a 2 mm pillar meshed at 0.05 mm they jump in
26 deg steps, and a ray field fitted tightly to them zig-zagged the loft's lower
cap over itself. Smoothing them costs nothing the overshoot does not absorb, and
smooth bands are what the intersection with the body handles reliably."""

SEW_TOL = 1e-6
MAX_DEGREE = 3
COS_T = float(np.cos(np.radians(SUPPORT_OVERHANG_DEG)))
SIN_T = float(np.sin(np.radians(SUPPORT_OVERHANG_DEG)))


# --- B-spline fitting -----------------------------------------------------------


def _knots(u: np.ndarray, n_ctrl: int, degree: int) -> np.ndarray:
    """Clamped knot vector on [0, 1] placed by the data parameters ``u``.

    Piegl and Tiller's averaging (The NURBS Book, eq. 9.69): every knot span
    holds data, which is the Schoenberg-Whitney condition that keeps the fit
    well-posed. Uniform knots do not, and on a cut run whose points crowd at one
    end the interpolation they gave was singular -- control points 1e5 mm away.
    """
    n = len(u)
    inner = []
    d = n / (n_ctrl - degree)
    for j in range(1, n_ctrl - degree):
        i = int(j * d)
        a = j * d - i
        inner.append((1.0 - a) * u[max(i - 1, 0)] + a * u[min(i, n - 1)])
    return np.r_[np.zeros(degree + 1), inner, np.ones(degree + 1)]


def _basis(u: np.ndarray, knots: np.ndarray, degree: int, n_ctrl: int) -> np.ndarray:
    """``(len(u), n_ctrl)`` B-spline basis values (Cox-de Boor)."""
    u = np.clip(np.asarray(u, dtype=float), 0.0, 1.0)
    N = np.zeros((len(u), len(knots) - 1))
    for i in range(len(knots) - 1):
        lo, hi = knots[i], knots[i + 1]
        N[:, i] = (u >= lo) & ((u < hi) | ((hi == 1.0) & (u == 1.0) & (lo < hi)))
    for p in range(1, degree + 1):
        nxt = np.zeros((len(u), len(knots) - 1 - p))
        for i in range(len(knots) - 1 - p):
            d1 = knots[i + p] - knots[i]
            d2 = knots[i + p + 1] - knots[i + 1]
            a = (u - knots[i]) / d1 * N[:, i] if d1 > 0 else 0.0
            b = (knots[i + p + 1] - u) / d2 * N[:, i + 1] if d2 > 0 else 0.0
            nxt[:, i] = a + b
        N = nxt
    return N[:, :n_ctrl]


@dataclass
class Fit:
    degree: int
    knots: np.ndarray
    ctrl: np.ndarray
    """``(m, k)`` control points of every fitted quantity stacked side by side."""


MIN_SPACING = 0.02
"""mm. Fitting data closer than this to the previous point kept is dropped."""

SMOOTHING = 1e-6
"""Weight, per data point, of the poles' second differences in a fit."""


def _thin(points: np.ndarray) -> np.ndarray:
    """Indices of ``points`` to fit: no two consecutive closer than :data:`MIN_SPACING`.

    Parting points crowd together where three cells meet, a few microns apart,
    and knots placed by their parameters crowd with them. The fit there hit its
    tolerance at the data with poles 1e7 mm away, cancelling between the points
    and nowhere else. Both ends are always kept.
    """
    keep = [0]
    for k in range(1, len(points) - 1):
        if np.linalg.norm(points[k] - points[keep[-1]]) >= MIN_SPACING:
            keep.append(k)
    if len(points) > 1:
        if len(keep) > 1 and np.linalg.norm(points[-1] - points[keep[-1]]) < MIN_SPACING:
            keep[-1] = len(points) - 1
        else:
            keep.append(len(points) - 1)
    return np.asarray(keep)


def _fit(u: np.ndarray, data: np.ndarray, tol_cols: list) -> Fit:
    """Least-squares clamped B-spline through ``data`` (``(n, k)``) at ``u``.

    The two ends are interpolated exactly, so consecutive runs meet at their
    shared point with its own ray. Control points are added until every column
    group ``(cols, tol)`` in ``tol_cols`` is within its tolerance, or the curve
    interpolates every point. A small penalty on the poles' second differences
    keeps the fit from buying accuracy at the data with poles that run wild
    between it.
    """
    n = len(u)
    degree = min(MAX_DEGREE, n - 1)
    m = degree + 1
    while True:
        knots = _knots(u, m, degree)
        N = _basis(u, knots, degree, m)
        ctrl = np.empty((m, data.shape[1]))
        ctrl[0], ctrl[-1] = data[0], data[-1]
        if m > 2:
            rhs = data - np.outer(N[:, 0], data[0]) - np.outer(N[:, -1], data[-1])
            A = N[:, 1:-1]
            if m > 3:
                D2 = np.zeros((m - 2, m))
                for i in range(m - 2):
                    D2[i, i:i + 3] = (1.0, -2.0, 1.0)
                lam = np.sqrt(SMOOTHING * n)
                A = np.vstack([A, lam * D2[:, 1:-1]])
                rhs = np.vstack([rhs, -lam * (np.outer(D2[:, 0], data[0])
                                              + np.outer(D2[:, -1], data[-1]))])
            sol, *_r = np.linalg.lstsq(A, rhs, rcond=None)
            ctrl[1:-1] = sol
        err = N @ ctrl - data
        ok = all(np.linalg.norm(err[:, cols], axis=1).max() <= tol for cols, tol in tol_cols)
        if ok or m >= n:
            return Fit(degree, knots, ctrl)
        m = min(n, m + max(1, m // 2))


def _smooth_in_basis(u: np.ndarray, data: np.ndarray, knots: np.ndarray, degree: int,
                     m: int, tol: float) -> np.ndarray:
    """The smoothest fit of ``data`` in a given basis that stays within ``tol``.

    Ends interpolated; the second-difference penalty is taken as large as the
    tolerance allows, searched downward by decades. A ray field needs no more
    than that: its targets snap from facet to facet of the obstacle mesh, and
    following them closely zig-zagged a converging bundle's lower end -- a 1.4 mm
    loop inside a pillar -- over itself.
    """
    N = _basis(u, knots, degree, m)
    n = len(u)
    ctrl = np.empty((m, data.shape[1]))
    ctrl[0], ctrl[-1] = data[0], data[-1]
    if m <= 2:
        return ctrl
    rhs0 = data - np.outer(N[:, 0], data[0]) - np.outer(N[:, -1], data[-1])
    if m == 3:
        sol, *_r = np.linalg.lstsq(N[:, 1:-1], rhs0, rcond=None)
        ctrl[1:-1] = sol
        return ctrl
    D2 = np.zeros((m - 2, m))
    for i in range(m - 2):
        D2[i, i:i + 3] = (1.0, -2.0, 1.0)
    best = None
    for lam in (1e2, 1e1, 1.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-6):
        s = np.sqrt(lam * n)
        A = np.vstack([N[:, 1:-1], s * D2[:, 1:-1]])
        rhs = np.vstack([rhs0, -s * (np.outer(D2[:, 0], data[0])
                                     + np.outer(D2[:, -1], data[-1]))])
        sol, *_r = np.linalg.lstsq(A, rhs, rcond=None)
        ctrl[1:-1] = sol
        best = ctrl.copy()
        if np.linalg.norm(N @ ctrl - data, axis=1).max() <= tol:
            return best
    return best


# --- runs ------------------------------------------------------------------------


_OWN = (cells_.VERTEX, cells_.CREASE)
"""Point kinds on the region's own border: one cell's, never shared."""


def _runs(loop) -> list:
    """Split a closed loop into runs: index lists, each sharing its ends."""
    n = len(loop.points)
    kind = loop.kind
    keys = loop.keys
    cut = np.zeros(n, dtype=bool)
    cut |= kind == cells_.TRIPLE
    for k in range(n):
        # A change of kind starts a run -- except right after a triple point,
        # which already ended the last one. Between a crossing, shared with the
        # neighbouring cell, and a border vertex, this cell's own, the run
        # breaks at the crossing whichever way the loop runs: so both cells end
        # their shared run on the same point and fit it as one curve. Breaking
        # at the change itself ended it on the vertex in one cell and on the
        # crossing in the other, and the two fits left slivers between them.
        if kind[k] != kind[k - 1] and kind[k - 1] != cells_.TRIPLE \
                and kind[k] != cells_.TRIPLE:
            if kind[k] in _OWN and kind[k - 1] not in _OWN:
                cut[k - 1] = True
            else:
                cut[k] = True
    starts = np.nonzero(cut)[0]
    if len(starts) == 0:
        starts = np.array([0, n // 2]) if n >= 4 else np.array([0])
    elif len(starts) == 1:
        starts = np.array([starts[0], (starts[0] + n // 2) % n])
    starts = np.sort(starts)
    runs = []
    for a, b in zip(starts, np.r_[starts[1:], starts[0] + n]):
        idx = np.arange(a, b + 1) % n
        runs.append(idx)
    return runs


def _chord(points: np.ndarray) -> np.ndarray:
    d = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
    return d / d[-1] if d[-1] > 0 else np.linspace(0.0, 1.0, len(points))


def _clamp_cone(vec: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Each row turned toward ``w`` until it is within the rod angle, length kept."""
    out = np.array(vec, dtype=float)
    L = np.linalg.norm(out, axis=1)
    z = out @ w
    h = out - np.outer(z, w)
    hn = np.linalg.norm(h, axis=1)
    over = (hn > z * np.tan(np.radians(SUPPORT_OVERHANG_DEG))) | (z <= 0)
    if over.any():
        hd = np.where(hn[:, None] > 0, h / np.maximum(hn, 1e-300)[:, None], 0.0)
        out[over] = L[over, None] * (COS_T * w + SIN_T * hd[over])
    return out


@dataclass
class Band:
    """One run's ruled face, as control rows."""

    top: np.ndarray
    bottom: np.ndarray
    degree: int
    knots: np.ndarray


def _canonical(loop, idx):
    """``(key, forward)`` of a run: its point ids in their canonical order."""
    ids = [loop.keys[k] for k in idx]
    forward = ids[0] <= ids[-1]
    return (tuple(ids) if forward else tuple(reversed(ids))), forward


SPAN_FIT_TOL = 0.2
"""mm. How far one span's ceiling curve may stray from its run before the run is
split -- inside the cell gate's :data:`COVER_SLACK`. Loose on purpose, and
measured: on `Hole-block-basic`, one cubic per run built all 20 cells, splitting
runs to 0.01 mm (302 spans of 57 runs) left 13 trims failing and to 0.08 mm
still some, and giving one span more control points instead (up to 6) failed
as often. The intersection with the body is reliable on simple bands and on
nothing else measured. A cubic follows a straight parting exactly and a 45 deg
arc to 1e-5 of its radius; only a genuinely wavy run is split, and the worst
fit is reported (``support_span_fit_mm``)."""

SPAN_MAX_POLES = 4
"""Control points per span: one cubic Bezier piece. Measured on
`Hole-block-basic`: bands fitted with one B-spline per run, 13 to 28 poles each,
left 16 of 20 cells' trims empty or invalid -- their faces' boundary wires did
not close after the intersection with the body -- where one cubic per span left
2. A run a single cubic cannot follow is split into more spans instead, so a
complex border costs faces, never wilder ones."""


SPAN_POLES_MAX = 4

_WORST_SPAN_FIT = 0.0
"""The worst ceiling-curve fit any span has had in this process, in mm."""


def worst_span_fit(reset: bool = False) -> float:
    """The worst span fit so far, optionally resetting it."""
    global _WORST_SPAN_FIT
    value = _WORST_SPAN_FIT
    if reset:
        _WORST_SPAN_FIT = 0.0
    return value
"""The most control points one span may carry before its run is split. A few
more poles in one face are what the trims tolerate; more faces are not."""


def _span(u: np.ndarray, data: np.ndarray, tol: float = np.inf):
    """One cubic B-spline through ``data`` at ``u``, ends interpolated.

    Starts as a single Bezier piece and takes control points, up to
    :data:`SPAN_POLES_MAX`, until its first three columns are within ``tol``.
    """
    n = len(u)
    degree = min(SPAN_MAX_POLES - 1, n - 1)
    m = degree + 1
    while True:
        knots = _knots(u, m, degree)
        N = _basis(u, knots, degree, m)
        ctrl = np.empty((m, data.shape[1]))
        ctrl[0], ctrl[-1] = data[0], data[-1]
        if m > 2:
            rhs = data - np.outer(N[:, 0], data[0]) - np.outer(N[:, -1], data[-1])
            sol, *_r = np.linalg.lstsq(N[:, 1:-1], rhs, rcond=None)
            ctrl[1:-1] = sol
        err = N @ ctrl - data
        if (np.linalg.norm(err[:, :3], axis=1).max() <= tol or m >= min(n, SPAN_POLES_MAX)):
            return Fit(degree, knots, ctrl), err
        m += 1


def _split_spans(P: np.ndarray, Ds: list, tol: float) -> list:
    """``[(a, b, fit)]``: index ranges of ``P`` each one span within tolerance.

    Ranges share their end points. A range is halved until its curve is within
    ``tol`` and every ray field in ``Ds`` within :data:`RAY_FIT_TOL`, or until it
    holds only two points.
    """
    out = []
    stack = [(0, len(P) - 1)]
    global _WORST_SPAN_FIT
    while stack:
        a, b = stack.pop()
        sl = slice(a, b + 1)
        u = _chord(P[sl])
        fit, err = _span(u, np.hstack([P[sl]] + [D[sl] for D in Ds]), tol)
        # Split on the ceiling curve alone: the rays are smoothed in the span's
        # basis afterwards, and splitting to follow their facet-to-facet noise
        # only made the bands busier and the trims fail again.
        ok = np.linalg.norm(err[:, :3], axis=1).max() <= tol
        if ok or b - a <= SPAN_MAX_POLES - 1:
            _WORST_SPAN_FIT = max(_WORST_SPAN_FIT, float(np.linalg.norm(err[:, :3], axis=1).max()))
            out.append((a, b, Fit(fit.degree, fit.knots, fit.ctrl[:, :3])))
        else:
            mid = (a + b) // 2
            stack.append((mid, b))
            stack.append((a, mid))
    out.sort(key=lambda r: r[0])
    return out


def prepare_fits(cells: list, cache: dict) -> None:
    """Fit every run's ceiling curve once, against every cell that uses it.

    A parting run belongs to two cells, each with its own rays along it. The
    curve is cut into spans and fitted here, with all of their rays in view, so
    both cells build on the same curve: within :data:`SPAN_FIT_TOL` of the
    run, and with every ray field within :data:`RAY_FIT_TOL` in its basis.
    """
    runs = {}
    for c in cells:
        loop = c.loops[0]
        for idx in _runs(loop):
            key, forward = _canonical(loop, idx)
            P = loop.points[idx]
            D = loop.targets[idx] - P
            if not forward:
                P, D = P[::-1], D[::-1]
            entry = runs.setdefault(key, {"P": P, "D": [], "kind": loop.kind[idx]})
            entry["D"].append(D)
    for key, entry in runs.items():
        if key not in cache:
            cache[key] = _fit_run(entry["P"], entry["D"], entry["kind"])


def _fit_run(P: np.ndarray, Ds: list, kind: np.ndarray) -> list:
    """Spans of one run in canonical order, indices into its unthinned points."""
    keep = _thin(P)
    spans = _split_spans(P[keep], [D[keep] for D in Ds], SPAN_FIT_TOL)
    return [(int(keep[a]), int(keep[b]), fit) for a, b, fit in spans]


def _bands(loop, w: np.ndarray, overshoot: float, cache: dict, extend=None) -> list:
    """The ruled bands of one loop, one per span, sharing curves through ``cache``.

    Each control ray is carried ``overshoot`` past its target into the anchor:
    along itself, or, given ``extend``, along that one direction for the whole
    loop (see :func:`loft`).
    """
    out = []
    for idx in _runs(loop):
        P = loop.points[idx]
        D = loop.targets[idx] - P
        key, forward = _canonical(loop, idx)
        spans = cache.get(key)
        if spans is None:
            Pc, Dc = (P, D) if forward else (P[::-1], D[::-1])
            spans = cache[key] = _fit_run(Pc, [Dc], loop.kind[idx])
        n = len(P)
        if not forward:
            spans = [(n - 1 - b, n - 1 - a, Fit(f.degree, 1.0 - f.knots[::-1], f.ctrl[::-1]))
                     for a, b, f in reversed(spans)]
        for a, b, fit in spans:
            sl = slice(a, b + 1)
            # Chord length is symmetric, so the span's own parameters are also
            # the canonical curve's, reversed when the run is.
            u = _chord(P[sl])
            ctrl_d = _smooth_in_basis(u, D[sl], fit.knots, fit.degree, len(fit.ctrl),
                                      RAY_FIT_TOL)
            ctrl_d = _clamp_cone(ctrl_d, w)
            if extend is None:
                L = np.linalg.norm(ctrl_d, axis=1)
                e = ctrl_d * (1.0 + overshoot / np.maximum(L, 1e-9))[:, None]
            else:
                e = ctrl_d + overshoot * extend
            out.append((fit.ctrl, e, fit.degree, fit.knots))
    return out


def _mean_direction(loop, w: np.ndarray) -> np.ndarray:
    """The loop's mean ray direction, inside the rod-angle cone about ``w``."""
    D = loop.targets - loop.points
    L = np.linalg.norm(D, axis=1)
    m = (D[L > 1e-9] / L[L > 1e-9][:, None]).mean(axis=0) if (L > 1e-9).any() else w
    m = _clamp_cone(m[None], w)[0]
    n = np.linalg.norm(m)
    return m / n if n > 1e-9 else w


PINCH_TOL = 1e-3
"""mm. Two neighbouring control rulings closer than this meet."""


def _rulings_cross(bands) -> bool:
    """Do neighbouring control rulings of some band cross below the ceiling?

    Rays converging on one target -- every ray from a radial cut line landing
    on the same point of a ledge's rim -- carried on along themselves fan out
    again below it: the band becomes a bowtie pinched at the target, valid
    as a B-rep and not trimmable. Seen as control rulings: two neighbours meet.
    Asking only whether the bottom row runs against the top one fires on
    ordinary curved runs too, and sent nine of Hole-block-basic's twenty cells
    down the wrong extension.
    """
    from ..classify import _segment_segment_dist

    for x, e, _d, _k in bands:
        b_ = x + e
        for i in range(len(x) - 1):
            gap = float(_segment_segment_dist(x[i], b_[i], x[i + 1:i + 2], b_[i + 1:i + 2])[0])
            if gap < PINCH_TOL:
                return True
    return False


# --- OCCT faces ------------------------------------------------------------------


def _surface(row0: np.ndarray, row1: np.ndarray, degree: int, knots: np.ndarray):
    m = len(row0)
    poles = TColgp_Array2OfPnt(1, m, 1, 2)
    for i in range(m):
        poles.SetValue(i + 1, 1, gp_Pnt(*map(float, row0[i])))
        poles.SetValue(i + 1, 2, gp_Pnt(*map(float, row1[i])))
    uniq, mult = np.unique(np.round(knots, 14), return_counts=True)
    uk = TColStd_Array1OfReal(1, len(uniq))
    um = TColStd_Array1OfInteger(1, len(uniq))
    for i, (k, c) in enumerate(zip(uniq, mult)):
        uk.SetValue(i + 1, float(k))
        um.SetValue(i + 1, int(c))
    vk = TColStd_Array1OfReal(1, 2)
    vm = TColStd_Array1OfInteger(1, 2)
    vk.SetValue(1, 0.0)
    vk.SetValue(2, 1.0)
    vm.SetValue(1, 2)
    vm.SetValue(2, 2)
    return Geom_BSplineSurface(poles, uk, vk, um, vm, degree, 1)


def _face(row0, row1, degree, knots):
    mk = BRepBuilderAPI_MakeFace(_surface(row0, row1, degree, knots), 1e-7)
    if not mk.IsDone():
        raise ProcessingError("A support cell face could not be built (docs/algorithm.md §14.3 S4).")
    return mk.Face()


@dataclass
class CellSolid:
    cell: object
    solid: object
    """The trimmed support body, or ``None``."""
    loft: object = None
    stats: dict = field(default_factory=dict)
    failure: str = ""


APEX_SEARCH_POLES = 24
"""Bottom poles below which a lower fan's apex is also looked for."""

FAN_SAMPLES = (0.25, 0.5, 0.75)
"""Fractions along the loft's extension beyond the anchor asked whether they are
inside the body."""


def _curve(poles: np.ndarray, degree: int, knots: np.ndarray):
    m = len(poles)
    arr = TColgp_Array1OfPnt(1, m)
    for i in range(m):
        arr.SetValue(i + 1, gp_Pnt(*map(float, poles[i])))
    uniq, mult = np.unique(np.round(knots, 14), return_counts=True)
    uk = TColStd_Array1OfReal(1, len(uniq))
    um = TColStd_Array1OfInteger(1, len(uniq))
    for i, (k, c) in enumerate(zip(uniq, mult)):
        uk.SetValue(i + 1, float(k))
        um.SetValue(i + 1, int(c))
    return Geom_BSplineCurve(arr, uk, um, degree)


def _folds(rows: list, degree_knots: list, axis: np.ndarray, per_run: int = 24) -> bool:
    """Does the closed chain of curves, seen along ``axis``, cross itself?"""
    pts = []
    for poles, (deg, knots) in zip(rows, degree_knots):
        u = np.linspace(0.0, 1.0, per_run, endpoint=False)
        pts.append(_basis(u, knots, deg, len(poles)) @ poles)
    P = np.vstack(pts)
    a = np.cross(axis, [1.0, 0.0, 0.0])
    if np.linalg.norm(a) < 1e-6:
        a = np.cross(axis, [0.0, 1.0, 0.0])
    a /= np.linalg.norm(a)
    c = np.cross(axis, a)
    Q = np.column_stack([P @ a, P @ c])
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


def _cap(rows: list, degree_knots: list):
    """A planar face bounded by the closed chain of curves through ``rows``."""
    wire = BRepBuilderAPI_MakeWire()
    for poles, (deg, knots) in zip(rows, degree_knots):
        wire.Add(BRepBuilderAPI_MakeEdge(_curve(poles, deg, knots)).Edge())
    if not wire.IsDone():
        raise ProcessingError("A support cell cap could not be closed (docs/algorithm.md §14.3 S4).")
    mk = BRepBuilderAPI_MakeFace(wire.Wire(), True)
    if not mk.IsDone():
        raise ProcessingError("A support cell cap is not planar (docs/algorithm.md §14.3 S4).")
    return mk.Face()


def _fan_faces(rows: list, bands: list, apex: np.ndarray) -> list:
    return [_face(r, np.repeat(apex[None], len(r), axis=0), d, k)
            for r, (_x, _e, d, k) in zip(rows, bands)]


def loft(cell, b: np.ndarray, overshoot: float, cache: dict, inside_body=None,
         extend: str = "auto"):
    """The closed, untrimmed loft of a cell (its outer loop), and how it was closed.

    The ruled band runs from just above the ceiling to ``overshoot`` past the
    anchor. Above, it is closed by a fan to an apex high in the shell; the body
    never reaches that fan. Below, it is closed one of two ways, chosen by
    ``inside_body(points)``:

    * ``"fan"`` -- to one apex just past the anchor, when that fan stays out of
      the body: behind a planar wall, inside a tube, inside a strut;
    * ``"plane"`` -- each control ray carried on, along itself, to a plane beyond
      the anchor, and capped there. Used where the material behind the anchor
      is not convex -- a round chamber's wall seen from inside -- and the fan's
      apex would land in the body. Rays to such a wall spread apart, so the
      carried-on loop keeps its order and the band its rulings.

    Projecting the whole bottom loop along one direction instead was tried and
    folds it: round a tube, where rays converge, and along a curved wall alike.

    How far past the anchor the rays reach is ``extend``: ``"ray"`` carries
    each control ray on along itself, which is what enters a wall the ray meets
    at an angle -- the tube crowns need it, since their mean ray is straight
    down a vertical wall; ``"mean"`` carries all of them along the loop's mean
    ray, which is what keeps rays converging on one target from crossing there
    (:func:`_rulings_cross`); ``"auto"`` is the first unless its rulings cross.
    Both stay self-supporting: a control ray plus a vector inside the cone is
    inside the cone.
    """
    w = -np.asarray(b, dtype=float)
    w /= np.linalg.norm(w)
    up = -w
    loop = cell.loops[0]
    direction = {"auto": None, "ray": None, "down": w}.get(extend)
    if extend == "mean":
        direction = _mean_direction(loop, w)
    bands = _bands(loop, w, overshoot, cache, direction)
    if extend == "auto" and _rulings_cross(bands):
        direction = _mean_direction(loop, w)
        bands = _bands(loop, w, overshoot, cache, direction)
    longest = max(float(np.linalg.norm(e, axis=1).max()) for _x, e, _d, _k in bands)
    v_top = -LIFT / longest
    tops = [x + v_top * e for x, e, _d, _k in bands]
    bottoms = [x + e for x, e, _d, _k in bands]
    dk = [(d, k) for _x, _e, d, k in bands]

    T = np.vstack(tops)
    c = T.mean(axis=0)
    rise = float(((T - c) @ up).max())
    radius = float(np.linalg.norm((T - c) - np.outer((T - c) @ up, up), axis=1).max())
    apex_top = c + (rise + radius * np.tan(np.radians(31.0)) + LIFT) * up

    E = np.vstack([e for _x, e, _d, _k in bands])
    mean_ray = E.mean(axis=0)
    mean_ray /= np.linalg.norm(mean_ray)
    Bm = np.vstack(bottoms)
    mode = "fan"
    # Candidate apexes for the lower fan, the first that keeps the fan out of
    # the body winning. Along the mean ray is right for a tube under a level
    # roof; under a sloping one the mean ray tilts with it and that apex can
    # leave the tube, where the bottom points' own centre still sits inside.
    if direction is None:
        T_ = np.vstack([x + e / (1.0 + overshoot / np.maximum(np.linalg.norm(e, axis=1), 1e-9))[:, None]
                        for x, e, _d, _k in bands])
    else:
        T_ = np.vstack([x + e - overshoot * direction for x, e, _d, _k in bands])
    candidates = [Bm.mean(axis=0) + overshoot * mean_ray, Bm.mean(axis=0),
                  T_.mean(axis=0) + overshoot * mean_ray, Bm.mean(axis=0) + overshoot * w]
    # Then below each bottom pole, deeper and deeper, nearest the centre first.
    # Rays converging on a ledge's rim pull the centre onto the rim itself, and
    # a fan from there reaches back over the ledge's edge into the body; an apex
    # further into the ledge sees the whole bottom loop.
    pick = Bm[np.linspace(0, len(Bm) - 1, min(len(Bm), APEX_SEARCH_POLES)).astype(int)]
    deeper = np.vstack([pick + k * overshoot * mean_ray for k in (1.0, 2.0, 4.0)])
    order = np.argsort(np.linalg.norm(deeper - Bm.mean(axis=0), axis=1))
    candidates += list(deeper[order])
    apex_bot = candidates[0]
    if inside_body is not None:
        for apex_bot in candidates:
            probe = np.vstack([Bm + f * (apex_bot - Bm) for f in FAN_SAMPLES])
            if not inside_body(probe).any():
                break
        else:
            mode = "plane"
    if mode == "plane":
        # Rays converging on a ledge's rim: what lies under the cell's middle is
        # the cone of rays from its middle down to the rim, which no closure
        # beyond the anchor can be. A fan from the bottom loop up to a point just
        # over the ceiling is that cone -- when every ruling of it is inside the
        # rod-angle cone, so that it can be a support's lower face.
        apex_up = c + (rise + LIFT / 2.0) * up
        rul = apex_up - Bm
        if np.all(rul @ up >= np.cos(np.radians(SUPPORT_OVERHANG_DEG)) * np.linalg.norm(rul, axis=1)):
            mode = "up"
            apex_bot = apex_up
    sew = BRepBuilderAPI_Sewing(SEW_TOL)
    for (x, e, deg, knots), t1, b1 in zip(bands, tops, bottoms):
        sew.Add(_face(t1, b1, deg, knots))
    for f in _fan_faces(tops, bands, apex_top):
        sew.Add(f)
    if mode in ("fan", "up"):
        for f in _fan_faces(bottoms, bands, apex_bot):
            sew.Add(f)
    else:
        d0 = float((Bm @ mean_ray).max()) + overshoot
        bottoms2 = []
        for b1, (_x, e, _d, _k) in zip(bottoms, bands):
            u = e / np.linalg.norm(e, axis=1, keepdims=True)
            along = u @ mean_ray
            if (along < 0.2).any():
                return None, "a ray runs too flat to the plane beyond its anchor", mode
            bottoms2.append(b1 + ((d0 - b1 @ mean_ray) / along)[:, None] * u)
        if _folds(bottoms2, dk, mean_ray):
            return None, "the loft's cap beyond the anchor folds over itself", mode
        if inside_body(np.vstack([Bm + f * (np.vstack(bottoms2) - Bm)
                                  for f in FAN_SAMPLES])).any():
            return None, "the loft's extension beyond the anchor reaches into the body", mode
        for (x, e, deg, knots), b1, b2 in zip(bands, bottoms, bottoms2):
            sew.Add(_face(b1, b2, deg, knots))
        sew.Add(_cap(bottoms2, dk))
    sew.Perform()
    shells = occ.shells(sew.SewedShape())
    if len(shells) != 1:
        return None, f"sewing the loft gave {len(shells)} shells", mode
    mk = BRepBuilderAPI_MakeSolid(shells[0])
    solid = mk.Solid()
    if not occ.is_valid(solid, parallel=False):
        return None, "the loft is not a valid solid", mode
    if occ.volume(solid) < 0.0:
        # Wound the other way round: BRepLib's orientation fix does not reach
        # the caller's handle through the bindings, so the solid is reversed
        # here, where the sign of its volume says it must be.
        solid = TopoDS.Solid_s(solid.Reversed())
    return solid, "", mode


FUZZY_LADDER = (0.0, 1e-5, 1e-4, 1e-3)
"""mm. Fuzzy values a cell's trim is tried with, in order, until one gives a valid
piece holding the seed. A fuzzy value moves nothing by more than itself, and
the cell gates check the result independently of the boolean that made it."""


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


FACE_TRIM_SEW_TOLS = (1e-5, 1e-4, 1e-3, 1e-2)
"""mm. Sewing tolerances the face-level trim is tried with, tightest first.

The loosest is what Hole-block-chamfered needs: a band meeting a tube wall where
the two separate face booleans place the intersection curve up to ~0.01 mm
apart, a 0.09 mm lens. It is a hundredth of the thinnest member this run can
ask for and far below the 0.5 mm support step, and the result still has to
pass every gate."""


def _face_trim(solid, body, seed: np.ndarray, lattice_near: list, tol: float):
    """The trim done face by face, for a cell the solid boolean fails on.

    Measured on the hole-blocks: intersecting a loft with the body can fail to
    rebuild a solid -- empty, or a face whose wire does not close -- while every
    face of it intersects the body correctly on its own. Each loft face is cut to
    the body and outside the lattice, each nearby body face to the loft, and the
    pieces are sewn; the shell holding the seed is the support.
    """
    lo, hi = occ.bounding_box(solid)
    pieces = []
    for f in occ.faces(solid):
        part = _boolean(BRepAlgoAPI_Common, f, [body])
        if part is None:
            return None
        if lattice_near:
            part = _boolean(BRepAlgoAPI_Cut, part, lattice_near)
            if part is None:
                return None
        pieces += occ.faces(part)
    for others in ([body], lattice_near):
        for shape in others:
            for g in occ.faces(shape):
                glo, ghi = occ.bounding_box(g)
                if np.any(glo > hi) or np.any(ghi < lo):
                    continue
                part = _boolean(BRepAlgoAPI_Common, g, [solid])
                if part is not None:
                    pieces += occ.faces(part)
    sew = BRepBuilderAPI_Sewing(tol)
    for p in pieces:
        sew.Add(p)
    sew.Perform()
    for sh in occ.shells(sew.SewedShape()):
        s = BRepBuilderAPI_MakeSolid(sh).Solid()
        if occ.volume(s) < 0.0:
            s = TopoDS.Solid_s(s.Reversed())
        held = _holding(s, seed)
        if held is not None:
            return held
    return None


def _holding(shape, point: np.ndarray):
    """The solid of ``shape`` containing ``point``, or ``None``."""
    for s in occ.solids(shape):
        cls = BRepClass3d_SolidClassifier(s)
        cls.Perform(gp_Pnt(*map(float, point)), 1e-7)
        if cls.State() == TopAbs_State.TopAbs_IN:
            return s
    return None


EXTEND_LADDER = ("auto", "mean", "ray", "down")
"""How a cell's rays are carried past their anchor (:func:`loft`), tried in turn
until one trims: ``auto`` first; then each fixed choice."""


def build_cell(cell, body, lattice_near: list, b: np.ndarray, overshoot: float,
               cache: dict, inside_body=None) -> CellSolid:
    """One support body: the cell's loft, within the body, outside the lattice."""
    whys = []
    solid = None
    for extend in EXTEND_LADDER:
        solid, why, mode = loft(cell, b, overshoot, cache, inside_body, extend)
        if solid is None:
            whys.append(f"{extend}: {why}")
            continue
        got = _trim(cell, solid, body, lattice_near)
        if got is not None:
            piece, stats = got
            stats.update({"closure": mode, "extend": extend})
            return CellSolid(cell, piece, solid, stats)
        whys.append(f"{extend}: no fuzzy value or face-by-face trim gave a valid "
                    f"piece holding the cell's seed point")
    return CellSolid(cell, None, solid, failure="; ".join(whys))


def _trim(cell, solid, body, lattice_near: list):
    """``(piece, stats)``: the loft within the body, outside the lattice, or ``None``."""
    for fuzzy in FUZZY_LADDER:
        inside = _boolean(BRepAlgoAPI_Common, solid, [body], fuzzy)
        if inside is None:
            continue
        if lattice_near:
            inside = _boolean(BRepAlgoAPI_Cut, inside, lattice_near, fuzzy)
            if inside is None:
                continue
        piece = _holding(inside, cell.seed)
        if piece is None or not occ.is_valid(piece, parallel=False):
            continue
        return piece, {"pieces": len(occ.solids(inside)), "fuzzy": fuzzy, "trim": "solid"}
    for tol in FACE_TRIM_SEW_TOLS:
        piece = _face_trim(solid, body, cell.seed, lattice_near, tol)
        if piece is not None and occ.is_valid(piece, parallel=False):
            return piece, {"trim": "faces", "sew_tol": tol}
    return None


# --- gates -------------------------------------------------------------------------

GATE_DEFLECTION = 0.02
"""mm. Tessellation the gates read a cell solid through."""

ON_BODY_TOL = 0.02
"""mm, beyond the body mesh's own deviation: a face this close to the body's
surface lies on it."""

SELF_SUPPORT_TOL_DEG = 0.05

COVER_SLACK = 0.25
"""mm. A ceiling point the cell's solid does not contain may still lie this close
to it. The solid's border is the fitted span, which follows the parting to
:data:`SPAN_FIT_TOL`, so the points along a parting fall either side of it,
and the neighbour across the parting holds those that fall outside. Whether
every overhang point is backed by material is asked of the written file, by
the verifier (specification.md §6.2)."""


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
    """The checks a cell solid must pass (docs/algorithm.md §14.3 S4), none of
    them the boolean that made it."""

    def __init__(self, body_mesh, b: np.ndarray):
        from ..classify import OutsideProbe, SpatialHash

        self.body_mesh = body_mesh
        self.probe = OutsideProbe(body_mesh, 2.0 * body_mesh.deviation)
        self.hash = SpatialHash(body_mesh, min_cell=1.0)
        self.b = np.asarray(b, dtype=float)
        self.bar = float(np.sin(np.radians(SUPPORT_OVERHANG_DEG + SELF_SUPPORT_TOL_DEG)))
        self.on_tol = body_mesh.deviation + ON_BODY_TOL
        from ..classify import PointInside

        self.inside_body = PointInside(body_mesh)

    def _dist_to_body(self, p: np.ndarray) -> float:
        from ..classify import _point_triangle_dist

        r = self.on_tol + 1e-6
        cand = self.hash.query(p - r, p + r)
        if not len(cand):
            return np.inf
        a, b, c = self.body_mesh.triangle_points
        return float(_point_triangle_dist(p, a[cand], b[cand], c[cand]).min())

    def check(self, solid, inside_points: np.ndarray) -> tuple[list, dict]:
        from ..classify import PointInside, TriMesh, _weld
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
        for face, pts, uvs, tris in per_face:
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
        verts, ts, offset = [], [], 0
        for _f, p, _u, t in per_face:
            verts.append(p)
            ts.append(t + offset)
            offset += len(p)
        v, t = _weld(np.vstack(verts), np.vstack(ts))
        covered = PointInside(TriMesh(v, t))(inside_points) if len(inside_points) else np.ones(0, bool)
        if (~covered).any():
            from ..classify import _point_triangle_dist

            a, b_, c = v[t[:, 0]], v[t[:, 1]], v[t[:, 2]]
            for k in np.nonzero(~covered)[0]:
                if float(_point_triangle_dist(inside_points[k], a, b_, c).min()) <= COVER_SLACK:
                    covered[k] = True
        missing = int((~covered).sum())
        if missing:
            fails.append(f"{missing} of {len(inside_points)} ceiling points not covered")
        return fails, {"free_faces": free_faces, "uncovered": missing,
                       "worst_free_face_deg": round(worst_deg, 3)}
