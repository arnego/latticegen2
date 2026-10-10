"""Support placement: where the apexes go (docs/algorithm.md §14.2-§14.4).

Kernel-free. Everything here is NumPy over the classification mesh, the lattice
expressions and the overhang samples; the geometry kernel only sees the convex
primitives this module returns.

**The validity principle.** Let ``K(s)`` be the upward 8-sided pyramid with apex
``s`` and ridge edges at the maximum overhang angle. ``K(s)`` intersected with
the body is self-supporting whenever ``s`` lies in grounded material or outside
the body: every point of it reaches ``s`` along a segment inside ``K(s)``, and
that segment either ends in material or leaves the body through a wall whose
own tilt it forces below the limit. So the whole support is a union of such
pyramids, each rooted in material that exists, and placement reduces to choosing
apexes — no visibility test, and no support surface that can exceed the limit.

**The lattice is the trunk system.** Apexes are taken, cheapest first, from

1. *strut stations* — points on the axes of surviving, grounded struts;
2. *king posts* — a thin vertical post from a grounded node up its cell's body
   diagonal, which is vertical and free of struts, carrying a capital;
3. *drop columns* — a post straight down to the wall. Going down from a ceiling
   point always leaves the body through an up-facing surface, so this fallback
   always exists.

Each sample chooses the candidate with the lowest estimated volume per covered
area. A capital's volume per unit of the area it covers grows linearly with its
depth below the ceiling, so shallowest-first is the right greedy.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .. import orient as _orient
from ..classify import (
    PointInside,
    SpatialHash,
    TriMesh,
    _DirectionalCaster,
    _point_triangle_dist,
    segment_triangle_dist,
)
from ..lattice import HALF_STRUTS, LatticeParams, half_strut_offset, nodes
from . import polyhedra as ph

NodeKey = tuple[int, int, int]

STATION_FRACTIONS = (0.15, 0.3, 0.42, 0.58, 0.72, 0.86)
"""Where along a strut, as a fraction of ``a`` from its node, an apex may sit.

Never at 0.5: that is the cap plane, a lattice cell face, and an apex on a cell
face is exactly the degenerate position the per-cell clip must never see.

The fractions past 0.5 lie on the *neighbouring* node's half of the strut. They
matter more than the near ones: a strut that runs up to an overhanging ceiling
ends at a node above it, so its upper half is where the lattice comes closest to
the surface, and an apex there makes a capital a fraction of the size of one
rooted a whole half-strut lower. They are only used where the strut is whole
all the way from the node (:func:`find_anchors`)."""

POST_LEVELS = 6
"""Heights tried for a king post's capital, spread over one cell's diagonal."""

PHASE_STEP = 0.6180339887498949
"""Golden-ratio conjugate: each primitive's azimuth phase is its index times
this, folded into one octagon sector, so no two support faces are coplanar by
construction of a regular array."""

LATTICE_AZIMUTH_CLEARANCE = np.radians(2.0)
"""Minimum azimuth between a support face normal and any lattice face or cap
normal, all of which sit at multiples of 60 degrees about the build direction."""

PREFERRED_DEPTH = 0.3
"""Depth below a sample, as a fraction of the cell edge ``a``, below which a
shallower apex is no longer treated as cheaper.

Volume alone would always prefer more and smaller capitals; this is the bound
on how small. It trades a little material for a number of support faces per
lattice cell the kernel can fuse in bounded time (docs/algorithm.md §14.4)."""

REACH = 1.25
"""Deepest an apex may sit below the sample it covers, in cell diagonals."""

ABSORB = 5.0
"""How far past the preferred depth a capital may grow to take over samples its
neighbours proposed for, as a multiple of that depth (docs/algorithm.md §14.4)."""

CELL_CLEARANCE = 2e-3
"""How close, as a fraction of the cell edge, a primitive vertex may come to a
lattice cell plane before its apex is nudged down. A vertex on a cell plane
makes one cell's clipped piece a sliver, and slivers are where a boolean and a
sew stop agreeing."""


@dataclass(frozen=True)
class SupportParams:
    overhang: float
    """Maximum overhang angle from the build direction, degrees."""
    thickness: float
    """Minimum support thickness, mm."""
    pitch: float
    """Overhang sample spacing, mm."""
    deviation: float
    """Measured deviation of the classification mesh, mm."""

    @property
    def tau(self) -> float:
        return float(np.tan(np.radians(self.overhang)))

    @property
    def tau_in(self) -> float:
        """Slope of a capital's *inscribed* cone — what coverage is judged by."""
        return self.tau * ph.INSCRIBED

    @property
    def margin(self) -> float:
        """Radius around a sample that must lie inside the capital covering it.

        Sample spacing plus mesh deviation: every surface point is within that
        of some sample, so covering each sample's neighbourhood covers the
        surface (docs/algorithm.md §14.2)."""
        return self.pitch + self.deviation

    @property
    def post_radius(self) -> float:
        return self.thickness / (2.0 * ph.INSCRIBED)

    @property
    def min_capital(self) -> float:
        """Smallest capital height: at least ``thickness`` across its top, and
        tall enough to swallow the top of the post it may stand on."""
        return max(self.thickness / (2.0 * self.tau_in),
                   1.5 * self.post_radius / self.tau)


def bed_frame(orient) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(hx, hy, b)``: the bed's two horizontal axes and its normal, in part
    coordinates."""
    R = _orient.rotation_matrix(orient)
    return R.T @ np.array([1.0, 0.0, 0.0]), R.T @ np.array([0.0, 1.0, 0.0]), \
        R.T @ np.array([0.0, 0.0, 1.0])


class RayCaster:
    """First hit of parallel rays against the classification mesh.

    Built on :class:`latticegen2.classify._DirectionalCaster`'s bucketing, so a
    ray only meets the triangles whose projection could contain it.
    """

    def __init__(self, mesh: TriMesh, direction: np.ndarray):
        self._c = _DirectionalCaster(mesh, np.asarray(direction, dtype=float))
        self._normals = np.cross(self._c.B - self._c.A, self._c.C - self._c.A)

    def first_hit(self, origins: np.ndarray):
        """``(t, outward_dot)`` per ray: distance to the nearest forward hit
        (``inf`` if none) and the unit triangle normal there, oriented along the
        ray. A ray leaving the body meets its outward normal that way."""
        c = self._c
        origins = np.asarray(origins, dtype=float)
        t_out = np.full(len(origins), np.inf)
        n_out = np.zeros((len(origins), 3))
        if len(origins) == 0:
            return t_out, n_out
        keys = np.floor((origins @ c.basis.T) / c.cell).astype(np.int64)
        order = np.lexsort((keys[:, 1], keys[:, 0]))
        start = 0
        EPS = 1e-9
        while start < len(order):
            end = start + 1
            k0 = keys[order[start]]
            while (end < len(order) and keys[order[end], 0] == k0[0]
                   and keys[order[end], 1] == k0[1]):
                end += 1
            idx = order[start:end]
            start = end
            tris = c.buckets.get((int(k0[0]), int(k0[1])))
            if tris is None or len(tris) == 0:
                continue
            A = c.A[tris]
            e1 = c.B[tris] - A
            e2 = c.C[tris] - A
            h = np.cross(c.d, e2)
            det = np.einsum("ij,ij->i", e1, h)
            ok = np.abs(det) >= EPS
            inv = np.where(ok, 1.0 / np.where(ok, det, 1.0), 0.0)
            s = origins[idx][:, None, :] - A[None, :, :]
            u = inv[None, :] * np.einsum("pij,ij->pi", s, h)
            qv = np.cross(s, e1[None, :, :])
            v = inv[None, :] * (qv @ c.d)
            t = inv[None, :] * np.einsum("pij,ij->pi", qv, e2)
            hit = (ok[None, :] & (u >= -EPS) & (u <= 1 + EPS) & (v >= -EPS)
                   & (u + v <= 1 + EPS) & (t > EPS))
            t = np.where(hit, t, np.inf)
            best = np.argmin(t, axis=1)
            tt = t[np.arange(len(idx)), best]
            t_out[idx] = tt
            nn = self._normals[tris[best]]
            nn = nn / np.maximum(np.linalg.norm(nn, axis=1, keepdims=True), 1e-300)
            flip = np.sign(nn @ c.d)
            flip[flip == 0.0] = 1.0
            n_out[idx] = nn * flip[:, None]
        return t_out, n_out


def near_samples(lp: LatticeParams, points: np.ndarray):
    """A predicate: is a position close enough to any sample to matter?

    An apex is never placed deeper below a sample than ``REACH`` cell
    diagonals, so a node further than that from every sample cannot root a
    support, and the per-station clearance checks — the expensive part of
    :func:`find_anchors` — are skipped for it. Decided on a coarse grid, which
    errs toward keeping a node.
    """
    cell = REACH * lp.a * np.sqrt(3.0)
    occupied = set()
    for c in np.unique(np.floor(points / cell).astype(np.int64), axis=0):
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    occupied.add((int(c[0]) + dx, int(c[1]) + dy, int(c[2]) + dz))

    def near(pos: np.ndarray) -> np.ndarray:
        keys = np.floor(pos / cell).astype(np.int64)
        return np.array([tuple(k) in occupied for k in keys.tolist()], dtype=bool)

    return near


@dataclass
class Anchors:
    """The lattice a support may root in (docs/algorithm.md §14.4)."""

    nodes: np.ndarray
    """``(M, 3)`` integer indices of eligible, grounded nodes."""
    stations: np.ndarray
    """``(M, 6, len(STATION_FRACTIONS))`` bool: is that strut station usable."""
    n_lattice: int = 0
    n_grounded: int = 0


def find_anchors(
    lp: LatticeParams,
    mesh: TriMesh,
    params: SupportParams,
    interior_nodes: np.ndarray,
    piece_nodes: list,
    interfaces: set,
    near=None,
    pool=None,
    mesh_path: str | None = None,
) -> Anchors:
    """Which lattice nodes and strut stations a support may root in.

    ``piece_nodes`` is the node of every surviving boundary piece (one entry
    per piece), ``interfaces`` the accepted cap interfaces, both as the lattice
    actually came out of the boundary trim — anchoring is evaluated against the
    lattice that exists, not the one that was planned
    (specification.md §4.6).

    **Grounded** means there is a way down to the wall through lattice that
    exists: from a node, down one of its three lower struts, either to a lower
    node that is itself grounded and joined across an accepted interface, or to
    where that strut emerges from a part of the wall that is not itself an
    overhang. Lattice hanging only from an overhanging ceiling is not a place
    to stand a support on.

    ``near``, if given, is a predicate over node positions that limits the
    (more expensive) station checks to nodes that could matter.
    """
    interior = {(int(r[0]), int(r[1]), int(r[2])) for r in interior_nodes}
    count: dict[NodeKey, int] = {}
    for n in piece_nodes:
        count[tuple(n)] = count.get(tuple(n), 0) + 1
    present = sorted(interior | set(count))
    if not present:
        return Anchors(np.empty((0, 3), dtype=np.int64),
                       np.empty((0, 6, len(STATION_FRACTIONS)), dtype=bool))
    idx = np.array(present, dtype=np.int64)
    pos = nodes(lp, idx)
    b = lp.b
    bar = float(np.sin(np.radians(params.overhang)))

    # Rays down each of the three lower struts, from every node at once.
    root = np.zeros((len(idx), 3), dtype=bool)
    open_ = np.zeros((len(idx), 3), dtype=bool)
    for k in range(3):
        t_hit, normal = RayCaster(mesh, -lp.e[k]).first_hit(pos)
        walled = t_hit <= lp.a
        root[:, k] = walled & (normal @ b <= bar)
        open_[:, k] = ~walled

    where = {n: i for i, n in enumerate(present)}
    grounded = np.zeros(len(idx), dtype=bool)
    for i in np.argsort(pos @ b, kind="stable"):
        n = present[i]
        if root[i].any():
            grounded[i] = True
            continue
        for k in range(3):
            if not open_[i, k]:
                continue
            lower = (n[0] - (k == 0), n[1] - (k == 1), n[2] - (k == 2))
            j = where.get(lower)
            if j is None or not grounded[j]:
                continue
            # Interior-to-interior caps share one index and are always joined;
            # anything touching a boundary piece must have been accepted.
            if (n in interior and lower in interior) or (n, k + 3) in interfaces:
                grounded[i] = True
                break

    clearance = lp.r + params.deviation
    inside = PointInside(mesh)(pos)
    eligible = grounded & inside
    if near is not None:
        eligible &= near(pos)
    # What the clearance checks need to know about each candidate node, decided
    # here from the junction graph so that the checks themselves depend on the
    # mesh alone and can be handed to worker processes as plain arrays.
    rows, whole_flags, joined_flags = [], [], []
    for i in np.nonzero(eligible)[0]:
        n = present[i]
        whole = n in interior
        # A boundary node is only a place to stand if its trim left exactly one
        # piece; whether that piece holds the node's centre is a clearance
        # question and is asked with the rest.
        if not whole and count.get(n, 0) != 1:
            continue
        joined = []
        for h in range(6):
            k, sgn = HALF_STRUTS[h]
            far = (n[0] + sgn * (k == 0), n[1] + sgn * (k == 1), n[2] + sgn * (k == 2))
            # Past the cap the material belongs to the neighbour's junction, so
            # it has to exist and to have been joined to this one.
            joined.append(far in where and (
                (whole and far in interior) or (n, h) in interfaces))
        rows.append(i)
        whole_flags.append(whole)
        joined_flags.append(joined)
    rows = np.array(rows, dtype=np.int64)
    whole_flags = np.array(whole_flags, dtype=bool)
    joined_flags = np.array(joined_flags, dtype=bool).reshape(-1, 6)

    if pool is not None and pool.active and mesh_path and len(rows) >= 4 * pool.workers:
        bounds = np.linspace(0, len(rows), pool.workers * 4 + 1).astype(int)
        jobs = [
            (mesh_path, lp.cc, lp.t, lp.orient, clearance,
             pos[rows[bounds[j]:bounds[j + 1]]],
             whole_flags[bounds[j]:bounds[j + 1]],
             joined_flags[bounds[j]:bounds[j + 1]])
            for j in range(len(bounds) - 1) if bounds[j + 1] > bounds[j]
        ]
        results, _rss = pool.run(_worker_stations, jobs)
        ok_node = np.concatenate([r[0] for r in results])
        ok_station = np.concatenate([r[1] for r in results])
    else:
        ok_node, ok_station = station_checks(
            lp, mesh, clearance, pos[rows], whole_flags, joined_flags)
    keep = rows[ok_node]
    stations = ok_station[ok_node]

    return Anchors(
        nodes=idx[keep],
        stations=stations,
        n_lattice=len(idx),
        n_grounded=int(grounded.sum()),
    )


_STATION_INDEX: dict = {}


def station_checks(lp, mesh, clearance, positions, whole, joined, hash_=None):
    """Which candidate nodes and strut stations are clear of the surface.

    Per node, independent of every other node, and a function of the mesh and
    plain arrays alone — the property that lets :func:`find_anchors` spread it
    across the worker pool the way classification is (docs/algorithm.md §5.4).
    Returns ``(node_ok, station_ok)``.
    """
    A, B, C = mesh.triangle_points
    if hash_ is None:
        hash_ = SpatialHash(mesh, min_cell=(lp.a / 2.0 + 2.0 * clearance) / 4.0)
    offsets = np.array([half_strut_offset(lp, h) for h in range(6)]) / (lp.a / 2.0)
    node_ok = np.ones(len(positions), dtype=bool)
    station_ok = np.ones((len(positions), 6, len(STATION_FRACTIONS)), dtype=bool)
    for i in range(len(positions)):
        p0 = positions[i]
        if not whole[i]:
            cand = hash_.query(p0 - clearance, p0 + clearance)
            if len(cand) and float(
                _point_triangle_dist(p0, A[cand], B[cand], C[cand]).min()
            ) <= clearance:
                node_ok[i] = False
                continue
        for h in range(6):
            for fi, frac in enumerate(STATION_FRACTIONS):
                if frac > 0.5 and not joined[i, h]:
                    station_ok[i, h, fi:] = False
                    break
                if whole[i] and frac < 0.5:
                    continue      # an interior node's own half-struts are clear
                p1 = p0 + offsets[h] * (frac * lp.a)
                lo = np.minimum(p0, p1) - clearance
                hi = np.maximum(p0, p1) + clearance
                cand = hash_.query(lo, hi)
                if len(cand) and float(segment_triangle_dist(
                    p0, p1, A[cand], B[cand], C[cand]
                ).min()) <= clearance:
                    station_ok[i, h, fi:] = False
                    break
    return node_ok, station_ok


def _worker_stations(job):
    """One chunk of :func:`station_checks` in a worker process.

    No geometry crosses the process boundary — the mesh is the staged ``.npz``
    and everything else is plain arrays — so, as for classification, neither
    the GIL finding nor the topology-identity finding of docs/algorithm.md §12
    has anything to attach to.
    """
    from ..classify import load_mesh
    from ..lattice import lattice_params
    from ..runlog import peak_rss_bytes

    mesh_path, cc, t, orient, clearance, positions, whole, joined = job
    lp = lattice_params(cc, t, orient)
    key = (mesh_path, cc, t, orient, clearance)
    got = _STATION_INDEX.get(key)
    if got is None:
        mesh = load_mesh(mesh_path)
        got = (mesh, SpatialHash(mesh, min_cell=(lp.a / 2.0 + 2.0 * clearance) / 4.0))
        _STATION_INDEX.clear()
        _STATION_INDEX[key] = got
    mesh, hash_ = got
    node_ok, station_ok = station_checks(
        lp, mesh, clearance, positions, whole, joined, hash_=hash_)
    return node_ok, station_ok, peak_rss_bytes()


@dataclass
class Plan:
    primitives: list = field(default_factory=list)
    stats: dict = field(default_factory=dict)


def _phase(index: int) -> float:
    """A deterministic azimuth phase clear of every lattice azimuth."""
    sector = 2.0 * np.pi / ph.SIDES
    phase = (index * PHASE_STEP) % 1.0 * sector
    for _ in range(ph.SIDES * 4):
        # Face normals sit half a sector off the ring corners.
        normals = phase + sector / 2.0 + np.arange(ph.SIDES) * sector
        off = np.abs((normals + np.pi / 6.0) % (np.pi / 3.0) - np.pi / 6.0)
        if float(off.min()) >= LATTICE_AZIMUTH_CLEARANCE:
            break
        phase = (phase + 1.25 * LATTICE_AZIMUTH_CLEARANCE) % sector
    return float(phase)


def _clear_of_cells(lp: LatticeParams, verts: np.ndarray) -> bool:
    idx = np.linalg.solve(lp.B, verts.T).T
    off = np.abs((idx - 0.5) - np.round(idx - 0.5))
    return bool(off.min() >= CELL_CLEARANCE)


def _settle(lp: LatticeParams, build, b: np.ndarray):
    """Build a primitive, nudging its apex down until it is in generic position.

    Lowering an apex is always safe: ``K`` only grows, so whatever it covered it
    still covers, and an apex buried in its anchor is still buried.
    """
    prim = build(np.zeros(3))
    for step in range(1, 12):
        if _clear_of_cells(lp, prim.verts):
            return prim
        prim = build(-step * 1.37 * CELL_CLEARANCE * lp.a * b)
    return prim


def _group(cells: np.ndarray) -> dict:
    out: dict[tuple, list] = {}
    for i, c in enumerate(map(tuple, cells.tolist())):
        out.setdefault(c, []).append(i)
    return {k: np.array(v, dtype=np.int64) for k, v in out.items()}


def plan_supports(
    lp: LatticeParams,
    mesh: TriMesh,
    params: SupportParams,
    points: np.ndarray,
    anchors: Anchors,
) -> Plan:
    """Choose apexes for every overhang sample and return the primitives."""
    plan = Plan()
    hx, hy, b = frame = bed_frame(lp.orient)
    tau, tau_in, margin = params.tau, params.tau_in, params.margin
    n_samples = len(points)
    if n_samples == 0:
        return plan

    reach = lp.a * np.sqrt(3.0) * REACH        # deepest apex considered
    grid = np.array([reach * tau_in, reach * tau_in, reach])
    to_bed = np.stack([hx, hy, b], axis=1)     # world -> bed coordinates

    P = points @ to_bed

    # ---- candidates -------------------------------------------------------
    cand_pos, cand_pen, cand_meta = [], [], []
    if len(anchors.nodes):
        npos = nodes(lp, anchors.nodes)
        unit = np.array([half_strut_offset(lp, h) for h in range(6)]) / (lp.a / 2.0)
        ceiling, _n = RayCaster(mesh, b).first_hit(npos)
        diag = lp.a * np.sqrt(3.0)
        post_pen = 3.0 * ph.POST_AREA_FACTOR * params.thickness ** 2 / (
            2.0 * np.sqrt(2.0) * tau * tau)
        for i in range(len(npos)):
            cand_pos.append(npos[i])
            cand_pen.append(0.0)
            cand_meta.append(("node", i, 0.0))
            for h in range(6):
                for fi, frac in enumerate(STATION_FRACTIONS):
                    if anchors.stations[i, h, fi]:
                        cand_pos.append(npos[i] + unit[h] * (frac * lp.a))
                        cand_pen.append(0.0)
                        cand_meta.append(("strut", i, 0.0))
            top = min(float(ceiling[i]) - 0.05 * lp.a, diag)
            for lvl in range(1, POST_LEVELS + 1):
                g = diag * lvl / (POST_LEVELS + 1)
                if g >= top:
                    break
                cand_pos.append(npos[i] + g * b)
                cand_pen.append(post_pen * g)
                cand_meta.append(("post", i, g))
    cand_pos = np.array(cand_pos).reshape(-1, 3)
    cand_pen = np.array(cand_pen, dtype=float)
    Q = cand_pos @ to_bed if len(cand_pos) else np.empty((0, 3))

    sample_cells = _group(np.floor(P / grid).astype(np.int64))
    cand_cells = _group(np.floor(Q / grid).astype(np.int64)) if len(Q) else {}

    def neighbours(cell):
        got = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0):
                    arr = cand_cells.get((cell[0] + dx, cell[1] + dy, cell[2] + dz))
                    if arr is not None:
                        got.append(arr)
        return np.concatenate(got) if got else np.empty(0, dtype=np.int64)

    def geometry(si, ci):
        dz = P[si, 2][:, None] - Q[ci, 2][None, :]
        rho = np.hypot(P[si, 0][:, None] - Q[ci, 0][None, :],
                       P[si, 1][:, None] - Q[ci, 1][None, :])
        feasible = (dz > 0.0) & (dz <= reach) & (rho + margin <= dz * tau_in)
        return dz, feasible

    # ---- round 1: every sample proposes an apex ----------------------------
    # Its cheapest feasible one, where "cheap" is depth below the sample: a
    # capital's volume per unit of covered area grows linearly with depth. Depth
    # is floored at `preferred`, though, because the true minimum-volume answer
    # is an unbounded number of vanishing capitals -- every one of them more
    # faces for the kernel to fuse. Below the floor a shallower apex is not
    # cheaper, so samples agree on apexes near it instead of each finding its
    # own.
    preferred = PREFERRED_DEPTH * lp.a
    absorb = ABSORB
    choice = np.full(n_samples, -1, dtype=np.int64)
    depth = np.zeros(n_samples)
    for cell, si in sample_cells.items():
        ci = neighbours(cell)
        if len(ci) == 0:
            continue
        dz, feasible = geometry(si, ci)
        base = np.where(dz >= preferred, dz, preferred + 0.5 * (preferred - dz))
        with np.errstate(divide="ignore", invalid="ignore"):
            price = np.where(feasible, base + cand_pen[ci][None, :] / (dz * dz), np.inf)
        best = np.argmin(price, axis=1)
        ok = np.isfinite(price[np.arange(len(si)), best])
        choice[si[ok]] = ci[best[ok]]
        depth[si[ok]] = dz[np.arange(len(si)), best][ok]

    # ---- round 2: greedy set cover among the proposed apexes ---------------
    # Each proposed apex gets the height its proposers need; then apexes are
    # taken largest-cover-first until every sample is covered, and the rest are
    # dropped. Nearby supports combine here: one capital already tall enough to
    # reach its neighbours' samples takes them over.
    top_margin = margin
    height = np.zeros(len(cand_pos))
    np.maximum.at(height, choice[choice >= 0], depth[choice >= 0])
    # A proposed capital may grow to `ABSORB` times the preferred depth to take
    # over its neighbours' samples; what it actually ends up as tall as is
    # decided by the samples it is finally given, below.
    height = np.where(
        height > 0.0,
        np.maximum.reduce([height + top_margin,
                           np.full(len(height), params.min_capital),
                           np.full(len(height), absorb * preferred)]),
        0.0)
    proposed = height > 0.0
    cover_s: list = [[] for _ in range(len(cand_pos))]
    cover_d: list = [[] for _ in range(len(cand_pos))]
    for cell, si in sample_cells.items():
        ci = neighbours(cell)
        ci = ci[proposed[ci]] if len(ci) else ci
        if len(ci) == 0:
            continue
        dz, feasible = geometry(si, ci)
        feasible &= dz + top_margin <= height[ci][None, :]
        for col in np.nonzero(feasible.any(axis=0))[0]:
            rows = np.nonzero(feasible[:, col])[0]
            cover_s[ci[col]].append(si[rows])
            cover_d[ci[col]].append(dz[rows, col])
    import heapq

    # Samples covered per unit of material is what is maximised, not samples
    # covered: by count alone the deepest capital always wins, because it has
    # the widest footprint, and the result is a ceiling backed almost solid
    # (measured on `test-cylinder.STEP`: 3.2 mm of solid-equivalent support
    # under capitals whose own cones average a third of that). Volume alone has
    # the opposite failure, so every capital also carries the fixed cost of one
    # at the preferred depth -- the same floor round 1 uses.
    cone = (2.0 * np.sqrt(2.0) / 3.0) * tau * tau
    fixed = cone * preferred ** 3
    covered = choice < 0          # samples nothing reaches are the fallback's
    weight = np.zeros(len(cand_pos))
    heap = []
    for c in np.nonzero(proposed)[0]:
        cover_s[c] = np.concatenate(cover_s[c]) if cover_s[c] else np.empty(0, np.int64)
        cover_d[c] = np.concatenate(cover_d[c]) if cover_d[c] else np.empty(0)
        if len(cover_s[c]):
            deep = float(cover_d[c].mean())
            weight[c] = 1.0 / (cone * deep ** 3 + cand_pen[c] * cone / 3.0 * 2.0
                               * np.sqrt(2.0) + fixed)
            heap.append((-len(cover_s[c]) * weight[c], int(c)))
    heapq.heapify(heap)
    height = np.zeros(len(cand_pos))
    while heap:
        neg, c = heapq.heappop(heap)
        fresh = ~covered[cover_s[c]]
        gain = float(fresh.sum()) * weight[c]
        if gain == 0.0:
            continue
        if heap and gain < -heap[0][0]:
            heapq.heappush(heap, (-gain, c))     # stale count: re-queue
            continue
        take = cover_s[c][fresh]
        covered[take] = True
        choice[take] = c
        depth[take] = cover_d[c][fresh]
        height[c] = float(cover_d[c][fresh].max())

    counts = {"node": 0, "strut": 0, "post": 0, "column": 0}
    serial = 0

    def add_capital(apex, h, anchor):
        nonlocal serial
        phase = _phase(serial)
        serial += 1
        plan.primitives.append(_settle(
            lp, lambda shift: ph.capital(apex + shift, h + float(-(shift @ b)),
                                         tau, phase, frame, anchor), b))

    def add_post(base, length, anchor):
        nonlocal serial
        phase = _phase(serial)
        serial += 1
        plan.primitives.append(_settle(
            lp, lambda shift: ph.post(base + shift, length + float(-(shift @ b)),
                                      params.thickness, tau, phase, frame, anchor), b))

    for c in np.nonzero(height > 0.0)[0]:
        kind, node_i, g = cand_meta[c]
        h = max(float(height[c]) + top_margin, params.min_capital)
        add_capital(cand_pos[c], h, "node" if kind == "post" else kind)
        counts[kind] += 1
        if kind == "post":
            # From the node's centre up into the capital, ending above the
            # height at which the capital has grown wider than the post.
            base = nodes(lp, anchors.nodes[node_i:node_i + 1])[0]
            add_post(base, g + 1.2 * params.post_radius / tau, "node")

    # ---- fallback: drop columns for whatever nothing reached ---------------
    left = np.nonzero(choice < 0)[0]
    if len(left):
        pitch = max(4.0 * params.thickness, 6.0 * params.pitch)
        down = RayCaster(mesh, -b)
        inside = PointInside(mesh)
        cells = _group(np.floor(P[left, :2] / pitch).astype(np.int64))
        for cell, li in cells.items():
            si = left[li]
            centre = (np.array(cell) + 0.5) * pitch
            rho = np.hypot(P[si, 0] - centre[0], P[si, 1] - centre[1])
            z_apex = float((P[si, 2] - (rho + margin) / tau_in).min())
            h = max(float(P[si, 2].max()) - z_apex + top_margin, params.min_capital)
            apex = centre[0] * hx + centre[1] * hy + z_apex * b
            add_capital(apex, h, "wall")
            counts["column"] += 1
            if bool(inside(apex[None, :])[0]):
                t_hit, _n = down.first_hit(apex[None, :])
                if np.isfinite(t_hit[0]):
                    drop = float(t_hit[0]) + 2.0 * params.deviation + 0.05
                    add_post(apex - drop * b,
                             drop + 1.2 * params.post_radius / tau, "wall")

    plan.stats = {
        "support_samples": n_samples,
        "support_capitals_on_nodes": counts["node"],
        "support_capitals_on_struts": counts["strut"],
        "support_king_posts": counts["post"],
        "support_drop_columns": counts["column"],
        "support_primitives": len(plan.primitives),
        "support_volume_untrimmed_mm3": round(
            sum(ph.volume(p.verts, p.faces) for p in plan.primitives), 3),
    }
    return plan
