"""Frustum supports on the lattice grid (docs/algorithm.md §14.3).

Every ceiling point -- a point of the input's inner surface overhanging more
than 60 deg -- is carried down to an **anchor**: a support node on a lattice
rod, or, where no rod is reached, a support node just outside the body, in the
enveloping wall. The support for that point is the **box of support nodes**
between its anchor and the support node holding the point. The support grid is
the lattice's own cube frame divided ``m`` times (:func:`grid_divisions`), so
every face of every box lies on a plane of the lattice frame: downward faces
lean 35.26 deg from vertical, every edge leans at most the rod angle, and the
lowest corner of every box is its anchor. The grid's planes run through the
rods' centre axes, so a box anchored on a rod has its lowest edge on the rod's
axis: it sits on the rod. Self-support is a property of the construction, not
a check, and so is the thinnest free-standing member: a support node beside a
rod loses only a ``t/2`` square corner of its section, so no member is thinner
than ``a/m - t/2`` (``a - t`` at ``m = 1``), which ``m`` is chosen to keep at
least ``--support-thickness``.

The union of the boxes is made **manifold** twice: once on the support grid,
and once on a fine grid refined with every rod face and cap plane, on which
every rod is made of whole fine support nodes, so the support leaves the rods'
nodes out by construction, shares faces with the lattice and needs no cut
against it. Two boxes passing a rod on its two lower sides meet along the
rod's axis, inside the rod, so leaving the rod out leaves them apart rather
than touching along an edge.

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
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt

from .. import occ
from ..errors import ProcessingError
from ..cli import support_grid
from ..runlog import point_text, warning_text
from .params import REQUIRED_SLACK, SIN_START

ANCHOR_DROP_PITCHES = 3
"""How many grid pitches *per division* a ceiling point's anchor is first
searched below its support node, along each rod direction. Every ceiling point
of the committed parts anchors within it; a point that does not is searched
further (:data:`ANCHOR_DROP_MAX_PITCHES`)."""

ANCHOR_DROP_MAX_PITCHES = 4096
"""The furthest an anchor is searched, in grid pitches: past any body this tool
can mesh. A ceiling point with no anchor within it fails the run."""

GRID_ORIGIN = 0.0
"""Where the support grid's planes sit, in half strut widths off the rods'
centre planes. ``0`` puts them through the rod axes; ``-1``, on the rods'
``-t/2`` faces, is the grid this generator used until 2026-11-05, kept
selectable for the comparison docs/algorithm.md §14.3 S4 records: there two
boxes passing a rod on its lower sides meet along its corner line, the
manifold repair fills the node below them, and that node rests on nothing --
33 of them on the ball, 110 on the cylinder."""

ANCHOR_EPS = 1e-3
"""mm. A support node is tested this far inside its lowest corner, along its
diagonal, so that a corner exactly on a rod face reads as inside the rod."""

FREE_SAMPLES = 6
"""Points per axis a fine support node near a lost junction piece is sampled at,
to ask whether any kept lattice is in it. A point in kept lattice settles it;
none settles nothing, and the node is then cut against the body."""

KEPT_MIN_VOLUME = 1e-8
"""mm³. A solid of a fine support node inside the body that is not lost counts as
kept lattice from this size, well under ``assemble.SUPPORT_OVERLAP_TOL``, so that
several in one body cannot add up to it."""

PROBE_DEPTH = 0.1
"""mm. How far into the body a ceiling point is moved to ask what backs it:
past the sag of a 0.02 mm tessellation of a curved ceiling face."""

SHALLOW_DEPTH = 2e-3
"""mm. The depth the harness's exact check asks at (`SUPPORT_EXACT_STEP`)."""

def lattice_backs(contains, P: np.ndarray, N: np.ndarray) -> np.ndarray:
    """Which ceiling points ``P``, outward normals ``N``, a strut backs: asked
    at both depths a support is asked at, :data:`PROBE_DEPTH` and
    :data:`SHALLOW_DEPTH`. A ceiling point resting on a strut's face lies in the
    strut a hair below it and can be past it at ``PROBE_DEPTH``: on the rehearsal
    part a keel lay in a strut to 0.08 mm, and was reported unbacked."""
    P = np.asarray(P, dtype=float)
    N = np.asarray(N, dtype=float)
    if not len(P):
        return np.zeros(0, dtype=bool)
    return (np.asarray(contains(P - PROBE_DEPTH * N), dtype=bool)
            | np.asarray(contains(P - SHALLOW_DEPTH * N), dtype=bool))


ESTIMATE_SAMPLES = 200_000
"""Points sampled in a support body's box to estimate its intersection with the
body, which every trim is held to."""

SPLIT_DEPTH = 4
"""How many times a support body no trim passes is split in two along a grid
plane, each half trimmed on its own, before the piece is dropped."""

DROPPED_WARNING = "support pieces dropped"
"""Opens the note a run prints when it drops support pieces (docs/algorithm.md
§14.3 S6), regardless of ``-v``."""

ESTIMATE_REL_TOL = 0.03
"""A trim's volume must agree with the sampled estimate to this fraction (or
four standard errors of the estimate, whichever is larger)."""


def grid_divisions(lp, thickness: float) -> int:
    """``m``: the largest division of the cube edge ``a`` whose thinnest
    free-standing support member is still at least ``thickness`` --
    ``a/m - t/2`` at ``m >= 2``, ``a - t`` at ``m = 1``
    (:func:`latticegen2.cli.support_grid`, the one copy of the rule).
    ``a = 7.07``, ``t = thickness = 1.5`` gives 3."""
    got = support_grid(lp.a, lp.t, thickness)
    if got is None:
        # `cli.parse_args` refuses such a value; reaching here is a caller that
        # bypassed it, and falling back to m = 1 would ship thinner support
        # than was asked for.
        raise ProcessingError(
            f"support thickness {thickness} mm cannot be met at a = {lp.a:.6f}, "
            f"t = {lp.t} mm: at most a - t = {lp.a - lp.t:.6f} mm")
    return got[0]


@dataclass
class Frustums:
    """What :func:`build` hands back."""

    solids: list
    """The trimmed support bodies, disjoint, each backing at least one ceiling
    point."""
    covered: np.ndarray
    """Per ceiling point: backed by one of :attr:`solids`, by the grid they were
    built from, or by a strut. A first answer: :func:`.generate.supports_for`
    asks the final bodies about the rest."""
    points: np.ndarray
    """The ceiling points."""
    normals: np.ndarray
    """Their outward normals."""
    stats: dict = field(default_factory=dict)
    lines: list = field(default_factory=list)
    notes: list = field(default_factory=list)
    """Lines for the console regardless of ``-v``: the pieces dropped."""
    cores: list = field(default_factory=list)
    """The staircase solids the trims were cut from, less the debris and the
    pieces dropped: inside the body exactly :attr:`solids`, with no face on the
    body's surface. What the pocket search cuts with where a cut by the
    trimmed solids is refused (§14.3 S7)."""


class _Grid:
    """The support grid: support node ``(i, j, k)`` spans
    ``[i q + o, (i+1) q + o]`` along each cube axis, ``o`` from
    :data:`GRID_ORIGIN` -- ``0``, so every ``m``-th plane runs through a rod's
    centre axis."""

    def __init__(self, lp, m: int):
        self.E = np.asarray(lp.e, dtype=float)      # rows e_k: cube coords C = x @ E.T
        self.m = m
        self.a = lp.a
        self.q = lp.a / m
        self.h = lp.t / 2.0
        self.o = GRID_ORIGIN * self.h

    def cube(self, x):
        return np.asarray(x, dtype=float) @ self.E.T

    def world(self, c):
        return np.asarray(c, dtype=float) @ self.E

    def node_of(self, x):
        return np.floor((self.cube(x) - self.o) / self.q).astype(np.int64)

    def corner(self, idx):
        """World position of support node ``idx``'s lowest corner, nudged inside it."""
        return self.world(np.asarray(idx) * self.q + self.o + ANCHOR_EPS)


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
    c0 = lo * g.q + g.o
    c1 = (lo + n) * g.q + g.o
    out = {round(float(x), 9) for x in (lo + np.arange(n + 1)) * g.q + g.o}
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


def backs_any(solid, P: np.ndarray) -> bool:
    """Whether ``solid`` holds any of ``P``, by the exact classifier, stopping at
    the first. A real support holds most of the ceiling points in its box, so
    this asks a handful; a piece of debris is small, and few lie in its box."""
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.TopAbs import TopAbs_State
    from OCP.gp import gp_Pnt

    P = np.asarray(P, dtype=float)
    lo, hi = (np.asarray(x, dtype=float) for x in occ.bounding_box(solid))
    cls = BRepClass3d_SolidClassifier(solid)
    for i in np.nonzero(np.all((P >= lo - 1e-6) & (P <= hi + 1e-6), axis=1))[0]:
        cls.Perform(gp_Pnt(*map(float, P[i])), 1e-7)
        if cls.State() in (TopAbs_State.TopAbs_IN, TopAbs_State.TopAbs_ON):
            return True
    return False


def _on_grid(g: "_Grid", coords: list, S: np.ndarray, P: np.ndarray) -> np.ndarray:
    """Which ``P`` lie in a support node of ``S`` on the fine grid ``coords``."""
    C = g.cube(P)
    idx = np.stack([np.searchsorted(coords[k], C[:, k]) - 1 for k in range(3)], axis=1)
    ok = np.all((idx >= 0) & (idx < np.array(S.shape)), axis=1)
    out = np.zeros(len(P), dtype=bool)
    out[ok] = S[idx[ok, 0], idx[ok, 1], idx[ok, 2]]
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


def _edge_uses(solid) -> list:
    """``(edge, number of face uses)`` for every edge of ``solid``."""
    from OCP.TopAbs import TopAbs_ShapeEnum
    from OCP.TopExp import TopExp
    from OCP.TopTools import TopTools_IndexedDataMapOfShapeListOfShape

    m = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(solid, TopAbs_ShapeEnum.TopAbs_EDGE,
                                   TopAbs_ShapeEnum.TopAbs_FACE, m)
    return [(m.FindKey(k), m.FindFromIndex(k).Size()) for k in range(1, m.Extent() + 1)]


def non_manifold_edges(solid) -> int:
    """Edges of ``solid`` used by more than two of its faces."""
    return sum(1 for _e, n in _edge_uses(solid) if n > 2)


def edges_not_used_twice(solid) -> int:
    """Non-degenerate edges of ``solid`` not used by exactly two of its faces."""
    from OCP.BRep import BRep_Tool
    from OCP.TopoDS import TopoDS

    return sum(1 for e, n in _edge_uses(solid)
               if n != 2 and not BRep_Tool.Degenerated_s(TopoDS.Edge_s(e)))


def _unified(s, gates):
    """``s`` with its coplanar faces merged -- unless unification moved anything.
    Measured on the cylinder: a valid unified solid with material 2.25 mm
    outside the body."""
    u = occ.unify_same_domain(s)
    if not occ.is_valid(u, parallel=False):
        return s
    if abs(occ.volume(u) - occ.volume(s)) > 1e-6 * occ.volume(s):
        return s
    return u if gates.outside(u) == 0.0 else s


def _trim(piece, body, gates, tmpdir: str, tag: str, full: bool = True, refused=None):
    """``(solids, lines, estimate)``: ``piece`` intersected with the body,
    checked; ``solids`` is ``None`` when no rung passed. ``estimate`` is
    ``(volume, standard error)``, the sampled volume of ``piece`` inside the
    body, mm^3.

    ``refused`` is the ``estimate`` of a ``full=False`` call already refused on
    this same piece: it is reused, and the fuzzy-0 rung that call ran is not
    run again.

    The full ladder runs the fuzzy values on the piece as built, then again on
    the piece as a STEP round trip gives it back: measured on the cylinder at
    m=1, the in-memory piece trims to material 2.25 mm outside the body at
    every fuzzy value and the round-tripped one trims cleanly at fuzzy 0. Every
    rung is held to the sampled estimate and to the gates.

    ``full=False`` tries fuzzy 0 alone, for a piece that can still be split:
    a refused rung on a large piece costs as much as the trim itself, and
    eight of them before every split ran a SpiralTest support stage past five
    hours, where its halves cost a fraction of that.
    """
    from .build import FUZZY_LADDER, _boolean

    lines = []
    est, sig = refused if refused is not None else _estimate(piece, gates.inside_body)
    if est < 3.0 * sig:
        return [], [f"  {tag}: no material inside the body by sampling; dropped"], (est, sig)
    rungs = [(piece, f) for f in (FUZZY_LADDER if full else FUZZY_LADDER[:1])]
    if refused is not None:
        rungs = rungs[1:]
    if full:
        path = os.path.join(tmpdir, f"support-{tag}.step")
        occ.write_step(piece, path, "support-piece")
        piece_rt = occ.solids(occ.read_step(path))[0]
        rungs += [(piece_rt, f) for f in FUZZY_LADDER]
    else:
        piece_rt = None
    for operand, fuzzy in rungs:
        got = _boolean(BRepAlgoAPI_Common, operand, [body], fuzzy)
        rung = f"fuzzy {fuzzy:g}{' after a STEP round trip' if operand is piece_rt else ''}"
        if got is None:
            lines.append(f"  {tag}, {rung}: the intersection did not complete; rejected")
            continue
        vin = sum(occ.volume(x) for x in occ.solids(got))
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
        return sol, lines, (est, sig)
    return None, lines, (est, sig)


def _core(piece, left: list) -> list:
    """``piece`` less ``left``, what its trim did not keep: inside the body, the
    support the trim produced. Empty where the cut does not take out exactly
    ``left``'s volume, which leaves the pocket search without this route there
    rather than with a wrong one."""
    if not left:
        return [piece]
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut

    from .build import _boolean

    got = _boolean(BRepAlgoAPI_Cut, piece, left)
    if got is None:
        return []
    out = occ.solids(got)
    gone = occ.volume(piece) - sum(occ.volume(s) for s in out)
    want = sum(occ.volume(s) for s in left)
    return out if abs(gone - want) <= max(1e-6 * occ.volume(piece), 1e-4) else []


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
        c = round(((clo[k] + chi[k]) / 2.0 - g.o) / g.q) * g.q + g.o
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


def _free_of_lattice(g: _Grid, coords: list, local, in_body, body, L: np.ndarray) -> np.ndarray:
    """``L`` with the fine support nodes a lost junction piece leaves empty taken
    out of it.

    Asked of the whole node, not its centre: a wall can cross a rod's fine node
    with a lost tip on one side and kept lattice on the other, and support there
    would overlap the lattice. Such a node stays lattice, and the tip's share of
    it stays bare."""
    if not local.lost:
        return L
    from .build import GATE_DEFLECTION, _mesh_of

    n = len(coords[0]) - 1, len(coords[1]) - 1, len(coords[2]) - 1
    L = L.copy()
    seen = set()
    for piece in local.lost:
        # Read off the piece's own surface: its world box, turned into the
        # lattice frame, is several times the piece. A node this misses stays
        # lattice, which leaves it empty and never overlaps anything.
        C = g.cube(np.vstack([p for _f, p, _u, _t in _mesh_of(piece, GATE_DEFLECTION)]))
        rng = [(int(np.searchsorted(coords[k], C[:, k].min() - GATE_DEFLECTION) - 1),
                int(np.searchsorted(coords[k], C[:, k].max() + GATE_DEFLECTION)))
               for k in range(3)]
        for i in range(max(rng[0][0], 0), min(rng[0][1], n[0])):
            for j in range(max(rng[1][0], 0), min(rng[1][1], n[1])):
                for k in range(max(rng[2][0], 0), min(rng[2][1], n[2])):
                    if (i, j, k) in seen:
                        continue
                    seen.add((i, j, k))
                    got = _holds_kept(g, coords, local, in_body, body, (i, j, k))
                    if got is not None:
                        L[i, j, k] = got
    return L


def _holds_kept(g: _Grid, coords: list, local, in_body, body, node) -> "bool | None":
    """Whether fine support node ``node`` holds kept lattice inside the body;
    ``None`` when none of it lies in a junction.

    Sampling can only say yes: a wedge of kept lattice slips between the samples,
    and support there overlaps the lattice. Measured on the rehearsal part: 23 of
    the 364 nodes the samples called free, wedges of 7e-05 to 2.4e-02 mm³. So a
    no is settled exactly, by cutting the node with the body. Each solid that
    leaves lies in one connected piece of the junction's trim, since all of a rod's
    fine node is junction, and a piece is lost or kept whole, so one point inside
    it decides. Where the kernel will not answer the node stays lattice: that
    leaves it empty, and overlaps nothing."""
    w = (np.arange(FREE_SAMPLES) + 0.5) / FREE_SAMPLES
    W = np.stack(np.meshgrid(w, w, w, indexing="ij"), -1).reshape(-1, 3)
    lo = np.array([coords[a][node[a]] for a in range(3)])
    hi = np.array([coords[a][node[a] + 1] for a in range(3)])
    Q = g.world(lo + W * (hi - lo))
    Q = Q[local.in_template(Q)]
    if not len(Q):
        return None
    Q = Q[in_body(Q)]
    if len(Q) and local.any_kept(Q):
        return True
    cut = BRepAlgoAPI_Common(_node_box(g, lo, hi), body)
    if not cut.IsDone():
        return True
    solids = occ.solids(cut.Shape())
    if len(Q) and not solids:
        return True
    return any(occ.volume(s) >= KEPT_MIN_VOLUME and not local.is_lost(s) for s in solids)


def _node_box(g: _Grid, lo: np.ndarray, hi: np.ndarray):
    """The box ``[lo, hi]`` in cube coordinates, as a solid. The frame is
    right-handed (`e0 x e1 = e2`), which is what the box's axes need."""
    axes = gp_Ax2(gp_Pnt(*map(float, g.world(lo))), gp_Dir(*map(float, g.E[2])),
                  gp_Dir(*map(float, g.E[0])))
    d = hi - lo
    return BRepPrimAPI_MakeBox(axes, float(d[0]), float(d[1]), float(d[2])).Shape()


def _contacts(S: np.ndarray) -> list:
    """``(nodes, upper)`` of every edge or corner contact left in ``S``: the
    window of nodes it lies in, and the one of the pair nearest the ceiling."""
    P = np.pad(S, 1)
    out = []
    for i, j in ((0, 1), (0, 2), (1, 2)):
        ei = np.eye(3, dtype=np.int64)[i]
        ej = np.eye(3, dtype=np.int64)[j]
        A, B, C, D = P, _roll(P, ei), _roll(P, ej), _roll(P, ei + ej)
        window = [np.zeros(3, dtype=np.int64), ei, ej, ei + ej]
        out += [(c - 1, window, ei + ej) for c in np.argwhere(A & D & ~B & ~C)]
        out += [(c - 1, window, ei) for c in np.argwhere(B & C & ~A & ~D)]
    for s1 in (1, -1):
        for s2 in (1, -1):
            d = np.array((1, s1, s2))
            window = [np.array((x, y, z)) for x in (0, 1) for y in (0, s1) for z in (0, s2)]
            bad = P & _roll(P, d)
            for o in window[1:-1]:
                bad &= ~_roll(P, o)
            upper = d if d.sum() > 0 else np.zeros(3, dtype=np.int64)
            out += [(c - 1, window, upper) for c in np.argwhere(bad)]
    return [(c + np.array(window), c + upper) for c, window, upper in out]


def _resolve_contacts(S: np.ndarray, L: np.ndarray, coords: list, g: _Grid, local, in_body, body):
    """``(S, filled, removed)``: ``S`` manifold on the fine grid beside the rods.

    :func:`make_manifold` will not fill a node of ``L``, and not all of ``L`` is
    material: a rod's node the body does not reach reads as lattice too, unless a
    lost piece happened to bring :func:`_free_of_lattice` to it. Measured on the
    rehearsal part: three contacts between such a node and one left standing,
    which made the 58,000 mm³ body invalid, so that no trim and no split could
    take it and it was dropped whole. So the nodes around a stuck contact are
    asked. One with kept lattice on both sides is real, and cannot be filled: the
    node of the pair nearest the ceiling goes, for the coverage count to report
    anything that leaves bare."""
    L = L.copy()
    S, filled, stuck = make_manifold(S.copy(), L)
    removed = 0
    asked = set()
    while stuck:
        freed = False
        found = _contacts(S)
        for nodes, _upper in found:
            for node in nodes:
                t = tuple(int(x) for x in node)
                if t in asked or min(t) < 0 or any(t[a] >= S.shape[a] for a in range(3)):
                    continue
                asked.add(t)
                if L[t] and not _holds_kept(g, coords, local, in_body, body, t):
                    L[t] = False
                    freed = True
        if not freed:
            for _nodes, upper in found:
                t = tuple(int(x) for x in upper)
                if S[t]:
                    S[t] = False
                    L[t] = True
                    asked.add(t)
                    removed += 1
        S, more, stuck = make_manifold(S, L)
        filled += more
    return S, filled, removed


def _clip_to_node(tri: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """The part of triangle ``tri`` (cube coordinates) inside the box ``[lo, hi]``,
    as a polygon; empty when none is. Sutherland-Hodgman, one plane at a time."""
    poly = tri
    for k in range(3):
        for sign, c in ((1.0, lo[k]), (-1.0, hi[k])):
            if not len(poly):
                return poly
            d = sign * (poly[:, k] - c)
            keep = []
            for i in range(len(poly)):
                j = (i + 1) % len(poly)
                if d[i] >= 0.0:
                    keep.append(poly[i])
                if (d[i] >= 0.0) != (d[j] >= 0.0):
                    keep.append(poly[i] + (poly[j] - poly[i]) * d[i] / (d[i] - d[j]))
            poly = np.array(keep) if keep else np.empty((0, 3))
    return poly


def _clipped_ceiling(g: _Grid, verts: np.ndarray, normals: np.ndarray, tris: np.ndarray,
                     held: set, backed, b: np.ndarray, bar: float, surface=None):
    """``(points, normals)``: a ceiling point in every support node the ceiling
    passes through and that holds none yet.

    The ceiling's refined vertices are ``h`` apart and a node is several times
    that, so a node nearly always holds one. Not where the surface only clips a
    corner, which takes a grazing contact: on the rehearsal part, cones and
    cylinders just past the bar, leaving ceiling bare 0.02 to 0.07 mm from the
    support beside it. So the triangles are walked across the grid, each clipped
    to the nodes it spans that hold no ceiling point; the point taken is the
    polygon's corner or centre that overhangs most, and it must itself be over
    the bar and not on a strut.

    ``surface(t, bary)`` -- ``(points, normals)`` on the true surface of triangle
    ``t``'s own face -- gives the normal there; without it the vertex normals
    are interpolated. Those are averaged where a vertex lies on a crease, and a
    crease between two faces each well under the bar can average to a ceiling:
    on the rehearsal part a roof apex of 13 and 30 deg faces read as 8.5 deg
    from flat, and its ceiling point was reported unbacked, there being nothing
    to support. A crease that does overhang is a keel, which §14.1 samples."""
    C = g.cube(verts)
    node = np.floor((C - g.o) / g.q).astype(np.int64)
    tn = node[tris]
    spans = np.nonzero((tn.min(axis=1) != tn.max(axis=1)).any(axis=1))[0]
    pts, nrm = [], []
    held = set(held)
    for t in spans:
        nlo, nhi = tn[t].min(axis=0), tn[t].max(axis=0)
        for i, j, k in itertools.product(*(range(nlo[a], nhi[a] + 1) for a in range(3))):
            if (i, j, k) in held:
                continue
            idx = np.array([i, j, k])
            lo = idx * g.q + g.o
            tri = C[tris[t]]
            poly = _clip_to_node(tri, lo, lo + g.q)
            if len(poly) < 3:
                continue
            cands = np.vstack([poly, poly.mean(axis=0)])
            # Barycentric weights of each candidate, for the normal there.
            A = np.vstack([tri.T, np.ones(3)])
            w = np.linalg.lstsq(A, np.vstack([cands.T, np.ones(len(cands))]), rcond=None)[0].T
            if surface is not None:
                n = np.array(surface(t, np.clip(w, 0.0, 1.0))[1], dtype=float)
            else:
                n = w @ normals[tris[t]]
            n /= np.linalg.norm(n, axis=1, keepdims=True)
            dots = n @ b
            pick = int(np.argmax(dots))
            if dots[pick] <= bar:
                continue
            # Nudged off the node's own faces, so that it reads as in this node.
            p = g.world(cands[pick] + 1e-9 * (lo + g.q / 2.0 - cands[pick]))
            if backed(p[None])[0]:
                continue
            held.add((i, j, k))
            pts.append(p)
            nrm.append(n[pick])
    if not pts:
        return np.empty((0, 3)), np.empty((0, 3))
    return np.array(pts), np.array(nrm)


def build(body, mesh, rmesh, local, lp, thickness: float, gates, tmpdir: str,
          report=None, extra=None) -> Frustums:
    """The support bodies for the ceiling points of ``rmesh``, before pockets.

    ``extra`` is ``(points, normals)``: further ceiling points, unless a strut
    already backs them."""
    need = rmesh.verts[rmesh.vert_required]
    need_n = rmesh.normals[rmesh.vert_required]
    if extra is not None and len(extra[0]):
        unbacked = ~lattice_backs(local.contains, extra[0], extra[1])
        need = np.vstack([need, extra[0][unbacked]])
        need_n = np.vstack([need_n, extra[1][unbacked]])
    m = grid_divisions(lp, thickness)
    g = _Grid(lp, m)
    bar = SIN_START + REQUIRED_SLACK
    over = rmesh.normals @ np.asarray(lp.b, dtype=float) > bar
    near_bar = over[rmesh.tris].any(axis=1) if len(rmesh.tris) else np.zeros(0, dtype=bool)
    near_ids = np.nonzero(near_bar)[0]
    exact = None
    if len(rmesh.tri_uv) == len(rmesh.tris) and rmesh.faces:
        def exact(t, w):
            return rmesh.surface_points(np.full(len(w), near_ids[t]), w)
    clipped, clipped_n = _clipped_ceiling(
        g, rmesh.verts, rmesh.normals, rmesh.tris[near_bar],
        {tuple(x) for x in g.node_of(need).tolist()}, local.contains,
        np.asarray(lp.b, dtype=float), bar, surface=exact)
    need = np.vstack([need, clipped])
    need_n = np.vstack([need_n, clipped_n])
    stats = {"support_ceiling_points": int(len(need)), "support_grid_divisions": m,
             "support_min_member_mm": round(support_grid(lp.a, lp.t, thickness)[1], 4),
             "support_ceiling_points_clipped": int(len(clipped))}
    lines = []
    if not len(need):
        return Frustums([], np.zeros(0, dtype=bool), need, need_n, stats, lines)

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
    own = np.floor((CC - g.o) / g.q).astype(np.int64) - lo
    S = U[own[:, 0], own[:, 1], own[:, 2]].reshape(rshape)
    L = np.concatenate([in_lattice(g.world(CC[i:i + 50_000]))
                        for i in range(0, len(CC), 50_000)]).reshape(rshape)
    L = _free_of_lattice(g, coords, local, gates.inside_body, body, L)
    S &= ~L
    S, rod_filled, rod_removed = _resolve_contacts(S, L, coords, g, local, gates.inside_body, body)
    stats.update({"support_nodes": int(U.sum()), "support_nodes_grown": grown,
                  "support_nodes_joined": filled, "support_nodes_beside_rods": rod_filled,
                  "support_nodes_removed_beside_rods": rod_removed})
    lines.append(
        f"support grid: m = {m}, pitch {g.q:.3f} mm; {len(need)} ceiling point(s), "
        f"{int(wall.sum())} anchored on the wall; {int(U.sum())} support node(s) "
        f"({grown} grown along the ceiling, {filled} joining an edge or corner contact), "
        f"{rod_filled} fine support node(s) filled beside the rods"
    )
    if rod_removed:
        lines.append(f"support grid: {rod_removed} fine support node(s) removed at contacts "
                     f"beside the rods that could not be filled")

    t1 = time.perf_counter()
    if report is not None:
        report("building support bodies", 0, None)
    pieces = _staircase(S, coords, g)
    t2 = time.perf_counter()
    trim_s = backed_s = 0.0
    solids, debris = [], []
    probes = need - PROBE_DEPTH * need_n
    # A sliver of support under a ceiling that only clips a node is thinner than
    # a probe is deep, so debris is also asked of the point just inside the ceiling.
    shallow = need - SHALLOW_DEPTH * need_n
    dropped = []      # (tag, centroid, mm^3 inside the body, ceiling points it held)
    dropped_pieces = []

    def trim_or_split(piece, tag: str, depth: int = 0):
        # One rung while the piece can still be split; the full ladder on a
        # piece that cannot.
        halves = None
        if depth < SPLIT_DEPTH:
            got, trim_lines, est = _trim(piece, body, gates, tmpdir, tag, full=False)
            if got is not None:
                return got, trim_lines
            halves = _split(piece, g)
            if halves is None:
                got, more, est = _trim(piece, body, gates, tmpdir, tag, refused=est)
                trim_lines += more
                if got is not None:
                    return got, trim_lines
        else:
            got, trim_lines, est = _trim(piece, body, gates, tmpdir, tag)
            if got is not None:
                return got, trim_lines
        if halves is None:
            # The kernel cannot trim it at any rung or size: dropped and
            # reported, as a junction no trim produces is (§7.2).
            where = np.round(np.asarray(occ.centroid(piece), dtype=float), 3).tolist()
            # Through its tessellation, not the exact classifier: a figure for a
            # report, and on a 10,000-face piece the classifier would take hours.
            held = _tessellated_inside(piece)(probes) & gates.inside_body(probes)
            dropped.append((tag, where, est[0], held))
            dropped_pieces.append(piece)
            trim_lines.append(f"  {tag}: no trim passed and no split is left; dropped")
            return [], trim_lines
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

    cores = []
    for i, piece in enumerate(pieces):
        if report is not None:
            report("trimming support bodies", i, len(pieces))
        ta = time.perf_counter()
        n_debris, n_dropped = len(debris), len(dropped_pieces)
        got, trim_lines = trim_or_split(piece, f"body-{i}")
        tb = time.perf_counter()
        trim_s += tb - ta
        lines += trim_lines
        for s in got:
            (solids if backs_any(s, probes) or backs_any(s, shallow) else debris).append(s)
        if got:
            cores += _core(piece, debris[n_debris:] + dropped_pieces[n_dropped:])
        backed_s += time.perf_counter() - tb
    if debris:
        lines.append(f"support debris discarded (backs no ceiling point): {len(debris)} "
                     f"solid(s), {sum(occ.volume(s) for s in debris):.3f} mm^3")
    # Backed by the support grid the bodies were built from -- every trim is
    # held to that grid's volume by its sampled estimate and gates -- or by a
    # strut. Classifying every ceiling point against the exact solids instead
    # took over half an hour on the rehearsal part: 6 ms a point on a
    # 5,000-face body, and 125,000 points.
    covered = _on_grid(g, coords, S, probes) | lattice_backs(in_lattice, need, need_n)
    # A probe outside the body is a ceiling point on the ceiling's own rim,
    # where a wall continues it: nothing below it can back it, and the wall does.
    covered |= ~gates.inside_body(probes)
    for _tag, _where, _est, held in dropped:
        covered &= ~held
    notes = []
    if dropped:
        # Every point a dropped piece held is uncovered by now, just above.
        lost = np.zeros(len(need), dtype=bool)
        for _tag, _where, _est, held in dropped:
            lost |= held
        notes.append(warning_text(
            DROPPED_WARNING,
            f"the geometry kernel could not trim {len(dropped)} piece(s) of support to "
            f"the input at any fuzzy value, after a STEP round trip, or in halves "
            f"(docs/algorithm.md §14.3 S6)",
            f"{sum(d[2] for d in dropped):.1f} mm^3 of support is missing, and "
            f"{int(lost.sum())} ceiling point(s) may print unsupported",
            [f"{tag} near {point_text(where)}: {est:.1f} mm^3, "
             f"{int(held.sum())} ceiling point(s)"
             for tag, where, est, held in dropped],
        ))
    stats.update({"support_debris_solids": len(debris),
                  "support_pieces_dropped": len(dropped),
                  "support_dropped_volume_mm3": round(sum(d[2] for d in dropped), 3),
                  "support_ceiling_points_uncovered": int((~covered).sum()),
                  "support_grid_s": round(t1 - t0, 2), "support_staircase_s": round(t2 - t1, 2),
                  "support_trim_s": round(trim_s, 2), "support_backed_s": round(backed_s, 2)})
    return Frustums(solids, covered, need, need_n, stats, lines, notes, cores)
