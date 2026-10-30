"""Frustum supports on the lattice grid (docs/algorithm.md §14.3).

Every ceiling point -- a point of the input's inner surface overhanging more
than 60 deg -- is carried down to an **anchor**: a support node on a lattice
rod, or, where no rod is reached, a support node just outside the body, in the
enveloping wall. The support for that point is the **box of support nodes**
between its anchor and the support node holding the point. The support grid is
the lattice's own cube frame divided ``m`` times (:func:`grid_divisions`), so
every face of every box lies on a plane of the lattice frame: downward faces
are rod faces (35.26 deg from horizontal), every edge leans at most the rod
angle, and the lowest corner of every box is its anchor. Self-support is a
property of the construction, not a check.

The union of the boxes is made **manifold** twice: once on the support grid,
and once on a fine grid refined with every rod face and cap plane, on which
every rod is made of whole fine support nodes -- because the support grid's
planes are rod face planes, a rod's corner line can otherwise lie exactly where
two support nodes meet diagonally, and once the rod is left out the support
touches itself along that line (an edge four faces share, which no CAD kernel
downstream accepts as a solid). On the fine grid the support leaves the rod's
nodes out by construction, so it shares faces with the lattice and needs no cut
against it.

The support bodies are then trimmed to the input body, one intersection per body,
each result held to a point-sampled estimate and to the material-outside probe,
because the kernel can return nothing, or material outside, while reporting
success (§7.2).
"""

from __future__ import annotations

import itertools
import math
import os
import time
from dataclasses import dataclass, field

import numpy as np
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common
from OCP.TopoDS import TopoDS

from .. import occ
from ..errors import ProcessingError

ANCHOR_DROP_PITCHES = 3
"""How many grid pitches *per division* a ceiling point's anchor is first
searched below its support node, along each rod direction. Every ceiling point
of the committed parts anchors within it; a point that does not is searched
further (:data:`ANCHOR_DROP_MAX_PITCHES`)."""

ANCHOR_DROP_MAX_PITCHES = 4096
"""The furthest an anchor is searched, in grid pitches: past any body this tool
can mesh. A ceiling point with no anchor within it fails the run."""

ANCHOR_EPS = 1e-3
"""mm. A support node is tested this far inside its lowest corner, along its
diagonal, so that a corner exactly on a rod face reads as inside the rod."""

PROBE_DEPTH = 0.1
"""mm. How far into the body a ceiling point is moved to ask what backs it:
past the sag of a 0.02 mm tessellation of a curved ceiling face."""

ESTIMATE_SAMPLES = 200_000
"""Points sampled in a support body's box to estimate its intersection with the
body, which every trim is held to."""

SPLIT_DEPTH = 4
"""How many times a support body no trim passes is split in two along a grid
plane, each half trimmed on its own, before the run fails."""

ESTIMATE_REL_TOL = 0.03
"""A trim's volume must agree with the sampled estimate to this fraction (or
four standard errors of the estimate, whichever is larger)."""


def grid_divisions(lp, thickness: float) -> int:
    """``m``: the largest division of the cube edge ``a`` whose gap between
    parallel rod faces, ``a/m - t``, is still at least ``thickness``; 1 at the
    least. ``a = 7.07``, ``t = thickness = 1.5`` gives 2."""
    return max(1, int(math.floor(lp.a / (lp.t + thickness) + 1e-9)))


@dataclass
class Frustums:
    """What :func:`build` hands back."""

    solids: list
    """The trimmed support bodies, disjoint, each backing at least one ceiling
    point."""
    covered: np.ndarray
    """Per ceiling point: backed by one of :attr:`solids`."""
    stats: dict = field(default_factory=dict)
    lines: list = field(default_factory=list)


class _Grid:
    """The support grid: support node ``(i, j, k)`` spans
    ``[i q - t/2, (i+1) q - t/2]`` along each cube axis, so its planes are the
    rods' ``-t/2`` face planes."""

    def __init__(self, lp, m: int):
        self.E = np.asarray(lp.e, dtype=float)      # rows e_k: cube coords C = x @ E.T
        self.m = m
        self.a = lp.a
        self.q = lp.a / m
        self.h = lp.t / 2.0

    def cube(self, x):
        return np.asarray(x, dtype=float) @ self.E.T

    def world(self, c):
        return np.asarray(c, dtype=float) @ self.E

    def node_of(self, x):
        return np.floor((self.cube(x) + self.h) / self.q).astype(np.int64)

    def corner(self, idx):
        """World position of support node ``idx``'s lowest corner, nudged inside it."""
        return self.world(np.asarray(idx) * self.q - self.h + ANCHOR_EPS)


def _anchors(g: _Grid, v: np.ndarray, in_body, in_lattice):
    """``(anchor, on_wall)``: per support node ``v`` holding a ceiling point,
    its anchor support node.

    Rod anchors lie on keel lines -- both transverse indices multiples of ``m``
    -- and are searched down the rod direction; wall anchors are any support
    node outside the body, searched down each rod direction. Of every candidate the
    highest (largest index sum) wins: the nearest anchor by drop.
    """
    n = len(v)
    best = np.full(n, np.iinfo(np.int64).min)
    anchor = np.zeros((n, 3), dtype=np.int64)
    wall = np.zeros(n, dtype=bool)
    reach = ANCHOR_DROP_PITCHES * g.m

    def take(gk, is_wall):
        todo = gk.sum(axis=1) > best
        if not todo.any():
            return
        P = g.corner(gk[todo])
        if is_wall:
            ok = ~in_body(P)
        else:
            ok = in_lattice(P)
            if ok.any():
                ok[ok] = in_body(P[ok])
        idx = np.nonzero(todo)[0][ok]
        best[idx] = gk[idx].sum(axis=1)
        anchor[idx] = gk[idx]
        wall[idx] = is_wall

    while True:
        # Rods first, so that a rod and a wall at the same height go to the rod.
        for k in range(3):
            keel = v.copy()
            for i in (i for i in range(3) if i != k):
                keel[:, i] = np.floor_divide(v[:, i], g.m) * g.m
            for drop in range(reach + 1):
                gk = keel.copy()
                gk[:, k] = v[:, k] - drop
                take(gk, False)
        for drop in range(reach + 1):
            for k in range(3):
                gk = v.copy()
                gk[:, k] -= drop
                take(gk, True)
        missing = best == np.iinfo(np.int64).min
        if not missing.any():
            return anchor, wall
        if reach >= ANCHOR_DROP_MAX_PITCHES:
            where = np.round(g.world(v[missing][0] * g.q), 3).tolist()
            raise ProcessingError(
                f"{int(missing.sum())} ceiling point(s) found no anchor on a rod or a wall "
                f"within {ANCHOR_DROP_MAX_PITCHES} grid pitches, the first near {where} "
                f"(docs/algorithm.md §14.3)."
            )
        reach *= 4


def _roll(X, d):
    return np.roll(X, [-x for x in d], axis=(0, 1, 2))


def _put(mask, d):
    return np.roll(mask, d, axis=(0, 1, 2))


def make_manifold(S: np.ndarray, L: np.ndarray | None = None):
    """``(S, filled, stuck)``: support nodes added to ``S`` until no two of its
    nodes meet only along an edge or at a corner, never filling a node of ``L``.

    Of a diagonal pair ``A``, ``D`` the two nodes between are filled, which rest
    on ``A``; of an anti-diagonal pair ``B``, ``C`` the node ``D`` above both,
    which rests on them -- or, where ``D`` is a rod, ``A`` below them. Two nodes
    touching at a corner are joined through a monotone path of two nodes.
    ``stuck`` counts contacts no fill could resolve, every candidate being
    lattice.
    """
    L = np.zeros_like(S) if L is None else L
    filled = 0
    while True:
        P = np.pad(S, 1)
        Lp = np.pad(L, 1)
        fill = np.zeros_like(P)
        stuck = np.zeros_like(P)
        for i, j in ((0, 1), (0, 2), (1, 2)):
            ei = [0, 0, 0]
            ei[i] = 1
            ej = [0, 0, 0]
            ej[j] = 1
            eij = [a + b for a, b in zip(ei, ej)]
            A, B, C, D = P, _roll(P, ei), _roll(P, ej), _roll(P, eij)
            LA, LB, LC, LD = Lp, _roll(Lp, ei), _roll(Lp, ej), _roll(Lp, eij)
            bad = A & D & ~B & ~C
            fill |= _put(bad & ~LB, ei) | _put(bad & ~LC, ej)
            stuck |= bad & LB & LC
            bad = B & C & ~A & ~D
            fill |= _put(bad & ~LD, eij) | (bad & LD & ~LA)
            stuck |= bad & LD & LA
        for s1 in (1, -1):
            for s2 in (1, -1):
                d = (1, s1, s2)
                others = [(x, y, z) for x in (0, 1) for y in (0, s1) for z in (0, s2)
                          if (x, y, z) not in ((0, 0, 0), d)]
                bad = P & _roll(P, d)
                for o in others:
                    bad &= ~_roll(P, o)
                todo = bad.copy()
                for a1, a2, _a3 in itertools.permutations([(1, 0, 0), (0, s1, 0), (0, 0, s2)]):
                    c2 = tuple(x + y for x, y in zip(a1, a2))
                    ok = todo & ~_roll(Lp, a1) & ~_roll(Lp, c2)
                    fill |= _put(ok, a1) | _put(ok, c2)
                    todo &= ~ok
                stuck |= todo
        fill &= ~P & ~Lp
        k = int(fill[1:-1, 1:-1, 1:-1].sum())
        if not k:
            return S, filled, int(stuck[1:-1, 1:-1, 1:-1].sum())
        S = S | fill[1:-1, 1:-1, 1:-1]
        filled += k


def _planes(g: _Grid, lo: int, n: int) -> np.ndarray:
    """Cube coordinates of the fine grid's planes along one axis: the support
    grid's planes of support nodes ``lo .. lo + n - 1`` plus every rod face
    (``±t/2``) and cap (``±a/2``) plane between them."""
    c0 = lo * g.q - g.h
    c1 = (lo + n) * g.q - g.h
    out = {round(float(x), 9) for x in (lo + np.arange(n + 1)) * g.q - g.h}
    for lattice_node in range(int(np.floor(c0 / g.a)) - 1, int(np.ceil(c1 / g.a)) + 2):
        for off in (-g.h, g.h, -g.a / 2.0, g.a / 2.0):
            x = lattice_node * g.a + off
            if c0 < x < c1:
                out.add(round(float(x), 9))
    return np.array(sorted(out))


def _staircase(S: np.ndarray, coords: list, g: _Grid) -> list:
    """One solid per connected body of support nodes ``S`` on the grid ``coords``."""
    quads = []
    P = np.pad(S, 1)
    for ax in range(3):
        o = [i for i in range(3) if i != ax]
        for sgn in (1, -1):
            face = P & ~np.roll(P, -sgn, axis=ax)
            for c in np.argwhere(face[1:-1, 1:-1, 1:-1]):
                base = c.copy()
                if sgn > 0:
                    base[ax] += 1
                corners = []
                for d0, d1 in ((0, 0), (1, 0), (1, 1), (0, 1)):
                    p = base.copy()
                    p[o[0]] += d0
                    p[o[1]] += d1
                    corners.append([coords[0][p[0]], coords[1][p[1]], coords[2][p[2]]])
                quads.append(g.world(np.array(corners)))
    if not quads:
        return []
    out = []
    for sh in occ.shells(occ.sew([occ.polygon_face(q) for q in quads], 1e-6)):
        s = occ.make_solid(sh)
        if occ.volume(s) < 0.0:
            s = TopoDS.Solid_s(s.Reversed())
        out.append(occ.unify_same_domain(s))
    return out


def _tessellated_inside(solid, deflection: float = 0.02):
    from ..classify import PointInside, TriMesh, _weld
    from .build import _mesh_of

    per = _mesh_of(solid, deflection)
    V = np.vstack([p for _f, p, _u, _t in per])
    off = np.cumsum([0] + [len(p) for _f, p, _u, _t in per])
    T = np.vstack([t + off[i] for i, (_f, _p, _u, t) in enumerate(per)])
    pin = PointInside(TriMesh(*_weld(V, T)))
    return lambda P: (np.concatenate([pin(P[i:i + 4000]) for i in range(0, len(P), 4000)])
                      if len(P) else np.zeros(0, dtype=bool))


def exact_inside(solid, P: np.ndarray) -> np.ndarray:
    """Which ``P`` lie inside ``solid`` or on its boundary, by OCCT's classifier
    on the exact B-rep. A freshly trimmed solid's welded tessellation can carry
    cracks, and ray parity through one misread 2,026 of 19,739 points as
    outside. A point on the boundary counts: beside a wall the support's face
    lies on the wall, and a ceiling point at the rim probes onto it."""
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.TopAbs import TopAbs_State
    from OCP.gp import gp_Pnt

    P = np.asarray(P, dtype=float)
    lo, hi = (np.asarray(x, dtype=float) for x in occ.bounding_box(solid))
    out = np.zeros(len(P), dtype=bool)
    cls = BRepClass3d_SolidClassifier(solid)
    for i in np.nonzero(np.all((P >= lo - 1e-6) & (P <= hi + 1e-6), axis=1))[0]:
        cls.Perform(gp_Pnt(*map(float, P[i])), 1e-7)
        out[i] = cls.State() in (TopAbs_State.TopAbs_IN, TopAbs_State.TopAbs_ON)
    return out


def _estimate(piece, in_body, seed: int = 1):
    """``(volume, standard error)`` of ``piece`` within the body, by sampling."""
    lo, hi = (np.asarray(x, dtype=float) for x in occ.bounding_box(piece))
    box = float(np.prod(hi - lo))
    P = lo + np.random.default_rng(seed).random((ESTIMATE_SAMPLES, 3)) * (hi - lo)
    ok = _tessellated_inside(piece)(P)
    if ok.any():
        ok[ok] = in_body(P[ok])
    return ok.mean() * box, math.sqrt(max(int(ok.sum()), 1)) / ESTIMATE_SAMPLES * box


def non_manifold_edges(solid) -> int:
    """Edges of ``solid`` used by more than two of its faces."""
    from OCP.TopAbs import TopAbs_ShapeEnum
    from OCP.TopExp import TopExp
    from OCP.TopTools import TopTools_IndexedDataMapOfShapeListOfShape

    m = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(solid, TopAbs_ShapeEnum.TopAbs_EDGE,
                                   TopAbs_ShapeEnum.TopAbs_FACE, m)
    return sum(1 for k in range(1, m.Extent() + 1) if m.FindFromIndex(k).Size() > 2)


def _unified(s, gates):
    """``s`` with its coplanar faces merged -- unless unification moved anything.
    Measured on the cylinder: a valid unified solid with material 2.25 mm
    outside the body."""
    from .build import GATE_DEFLECTION, _mesh_of

    u = occ.unify_same_domain(s)
    if not occ.is_valid(u, parallel=False):
        return s
    if abs(occ.volume(u) - occ.volume(s)) > 1e-6 * occ.volume(s):
        return s
    pts = np.vstack([p for _f, p, _u, _t in _mesh_of(u, GATE_DEFLECTION)])
    return u if gates.probe.worst_outside(pts) == 0.0 else s


def _trim(piece, body, gates, tmpdir: str, tag: str):
    """``(solids, lines)``: ``piece`` intersected with the body, checked;
    ``solids`` is ``None`` when no rung passed.

    The ladder runs the fuzzy values on the piece as built, then again on the
    piece as a STEP round trip gives it back: measured on the cylinder at m=1,
    the in-memory piece trims to material 2.25 mm outside the body at every
    fuzzy value and the round-tripped one trims cleanly at fuzzy 0. Every rung
    is held to the sampled estimate and to the gates.
    """
    from .build import FUZZY_LADDER, _boolean

    lines = []
    est, sig = _estimate(piece, gates.inside_body)
    if est < 3.0 * sig:
        return [], [f"  {tag}: no material inside the body by sampling; dropped"]
    path = os.path.join(tmpdir, f"support-{tag}.step")
    occ.write_step(piece, path, "support-piece")
    piece_rt = occ.solids(occ.read_step(path))[0]
    for operand, fuzzy in ([(piece, f) for f in FUZZY_LADDER]
                           + [(piece_rt, f) for f in FUZZY_LADDER]):
        got = _boolean(BRepAlgoAPI_Common, operand, [body], fuzzy)
        vin = 0.0 if got is None else sum(occ.volume(x) for x in occ.solids(got))
        rung = f"fuzzy {fuzzy:g}{' after a STEP round trip' if operand is piece_rt else ''}"
        if abs(vin - est) > max(ESTIMATE_REL_TOL * est, 4.0 * sig):
            lines.append(f"  {tag}, {rung}: {vin:.1f} mm^3 against {est:.1f} sampled; rejected")
            continue
        sol = [_unified(s, gates) for s in occ.solids(got)]
        why = []
        for s in sol:
            fails = list(gates.check(s)[0])
            nm = non_manifold_edges(s)
            if nm:
                fails.append(f"{nm} edge(s) used by more than two faces")
            if fails:
                why.append(f"{occ.volume(s):.1f} mm^3: {'; '.join(fails)}")
        if why:
            lines.append(f"  {tag}, {rung}: rejected, {' | '.join(why)}")
            continue
        return sol, lines
    return None, lines


def _split(piece, g: _Grid):
    """The two halves of ``piece`` either side of a support grid plane near its
    middle, across its longest extent in the lattice frame that a grid plane
    crosses; ``None`` where none does. A staircase cut along one of its own
    planes is a well-conditioned boolean, and the halves are checked to add up.
    Extents are read off the piece's vertices: its world box, turned into the
    lattice frame, is loose enough that its middle can miss the piece."""
    from OCP.BRep import BRep_Tool
    from OCP.TopAbs import TopAbs_ShapeEnum
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS as _T

    from .build import _boolean

    pts = []
    ex = TopExp_Explorer(piece, TopAbs_ShapeEnum.TopAbs_VERTEX)
    while ex.More():
        q = BRep_Tool.Pnt_s(_T.Vertex_s(ex.Current()))
        pts.append((q.X(), q.Y(), q.Z()))
        ex.Next()
    C = g.cube(np.array(pts))
    clo, chi = C.min(axis=0), C.max(axis=0)
    v = occ.volume(piece)
    for k in np.argsort(clo - chi):                 # longest extent first
        c = round(((clo[k] + chi[k]) / 2.0 + g.h) / g.q) * g.q - g.h
        if not clo[k] + 1e-6 < c < chi[k] - 1e-6:
            continue
        halves = []
        for a, b in ((clo[k] - 1.0, c), (c, chi[k] + 1.0)):
            coords = [np.array([clo[i] - 1.0, chi[i] + 1.0]) for i in range(3)]
            coords[k] = np.array([a, b])
            box = _staircase(np.ones((1, 1, 1), dtype=bool), coords, g)[0]
            got = _boolean(BRepAlgoAPI_Common, piece, [box])
            halves += [] if got is None else occ.solids(got)
        if len(halves) >= 2 and abs(sum(occ.volume(h) for h in halves) - v) <= 1e-9 * v:
            return halves
    return None


def build(body, mesh, rmesh, local, lp, thickness: float, gates, tmpdir: str,
          report=None) -> Frustums:
    """The support bodies for the ceiling points of ``rmesh``, before pockets."""
    need = rmesh.verts[rmesh.vert_required]
    need_n = rmesh.normals[rmesh.vert_required]
    m = grid_divisions(lp, thickness)
    g = _Grid(lp, m)
    stats = {"support_ceiling_points": int(len(need)), "support_grid_divisions": m}
    lines = []
    if not len(need):
        return Frustums([], np.zeros(0, dtype=bool), stats, lines)

    def in_lattice(P):
        return local.contains(P) if len(local.positions) else np.zeros(len(P), dtype=bool)

    if report is not None:
        report("anchoring ceiling points", 0, None)
    t0 = time.perf_counter()
    v = g.node_of(need)
    anchor, wall = _anchors(g, v, gates.inside_body, in_lattice)
    stats["support_wall_anchors"] = int(wall.sum())

    # Boxes [anchor, ceiling point's node], then growth along the ceiling: a
    # support node holding a ceiling point, or the point PROBE_DEPTH into the
    # body, that rests on support directly below it along a rod direction.
    pairs = np.unique(np.hstack([anchor, v]), axis=0)
    lo = pairs[:, :3].min(axis=0) - 1
    hi = pairs[:, 3:].max(axis=0) + 2
    shape = tuple((hi - lo).tolist())
    U = np.zeros(shape, dtype=bool)
    for r in pairs:
        a_, c_ = r[:3] - lo, r[3:] - lo
        U[a_[0]:c_[0] + 1, a_[1]:c_[1] + 1, a_[2]:c_[2] + 1] = True
    ceiling = np.zeros(shape, dtype=bool)
    for P in (need, need - PROBE_DEPTH * need_n):
        vi = g.node_of(P) - lo
        ok = np.all((vi >= 0) & (vi < np.array(shape)), axis=1)
        ceiling[vi[ok, 0], vi[ok, 1], vi[ok, 2]] = True
    grown = 0
    while True:
        below = np.zeros(shape, dtype=bool)
        Pd = np.pad(U, 1)
        for k in range(3):
            below |= np.roll(Pd, 1, axis=k)[1:-1, 1:-1, 1:-1]
        add = below & ceiling & ~U
        if not add.any():
            break
        grown += int(add.sum())
        U |= add
    U, filled, _ = make_manifold(U)

    # The fine grid: every rod is whole fine support nodes, and the support
    # leaves them out.
    coords = [_planes(g, int(lo[ax]), shape[ax]) for ax in range(3)]
    centres = [(c[:-1] + c[1:]) / 2.0 for c in coords]
    rshape = tuple(len(c) for c in centres)
    CC = np.stack(np.meshgrid(*centres, indexing="ij"), -1).reshape(-1, 3)
    own = np.floor((CC + g.h) / g.q).astype(np.int64) - lo
    S = U[own[:, 0], own[:, 1], own[:, 2]].reshape(rshape)
    L = np.concatenate([in_lattice(g.world(CC[i:i + 50_000]))
                        for i in range(0, len(CC), 50_000)]).reshape(rshape)
    S &= ~L
    S, rod_filled, stuck = make_manifold(S, L)
    stats.update({"support_nodes": int(U.sum()), "support_nodes_grown": grown,
                  "support_nodes_joined": filled, "support_nodes_beside_rods": rod_filled})
    lines.append(
        f"support grid: m = {m}, pitch {g.q:.3f} mm; {len(need)} ceiling point(s), "
        f"{int(wall.sum())} anchored on the wall; {int(U.sum())} support node(s) "
        f"({grown} grown along the ceiling, {filled} joining an edge or corner contact), "
        f"{rod_filled} fine support node(s) filled beside the rods"
    )
    if stuck:
        lines.append(f"support grid: {stuck} contact(s) beside the rods could not be filled")

    t1 = time.perf_counter()
    if report is not None:
        report("building support bodies", 0, None)
    pieces = _staircase(S, coords, g)
    t2 = time.perf_counter()
    trim_s = backed_s = 0.0
    solids, debris = [], []
    covered = np.zeros(len(need), dtype=bool)
    probes = need - PROBE_DEPTH * need_n
    def trim_or_split(piece, tag: str, depth: int = 0):
        got, trim_lines = _trim(piece, body, gates, tmpdir, tag)
        if got is not None:
            return got, trim_lines
        halves = _split(piece, g) if depth < SPLIT_DEPTH else None
        if halves is None:
            where = np.round(np.asarray(occ.centroid(piece), dtype=float), 3).tolist()
            raise ProcessingError(
                f"The support body near {where} could not be trimmed to the input body at "
                f"any fuzzy value, nor in halves "
                f"({' '.join(x.strip() for x in trim_lines)}) (docs/algorithm.md §14.3)."
            )
        from .join import join

        trim_lines.append(f"  {tag}: no trim passed; split in {len(halves)} along a grid plane")
        out = []
        for j, half in enumerate(halves):
            part, more = trim_or_split(half, f"{tag}.{j}", depth + 1)
            out += part
            trim_lines += more
        # The halves meet on the split plane: one body again.
        joined, _js = join(out, overlaps=False)
        return joined, trim_lines

    for i, piece in enumerate(pieces):
        if report is not None:
            report("trimming support bodies", i, len(pieces))
        ta = time.perf_counter()
        got, trim_lines = trim_or_split(piece, f"body-{i}")
        tb = time.perf_counter()
        trim_s += tb - ta
        lines += trim_lines
        for s in got:
            hold = exact_inside(s, probes)
            if hold.any():
                solids.append(s)
                covered |= hold
            else:
                debris.append(s)
        backed_s += time.perf_counter() - tb
    if debris:
        lines.append(f"support debris discarded (backs no ceiling point): {len(debris)} "
                     f"solid(s), {sum(occ.volume(s) for s in debris):.3f} mm^3")
    covered |= in_lattice(probes)
    # A probe outside the body is a ceiling point on the ceiling's own rim,
    # where a wall continues it: nothing below it can back it, and the wall does.
    covered |= ~gates.inside_body(probes)
    stats.update({"support_debris_solids": len(debris),
                  "support_ceiling_points_uncovered": int((~covered).sum()),
                  "support_grid_s": round(t1 - t0, 2), "support_staircase_s": round(t2 - t1, 2),
                  "support_trim_s": round(trim_s, 2), "support_backed_s": round(backed_s, 2)})
    return Frustums(solids, covered, stats, lines)
