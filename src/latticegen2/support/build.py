"""Building supports as B-rep, one lattice cell at a time (docs/algorithm.md §14.5).

Fusing supports into a finished lattice is not an option: that is one boolean
against hundreds of thousands of faces, the operation the whole architecture
exists to avoid. Supports therefore enter **per junction, before assembly**, and
the lattice already provides the partition. Node ``n`` owns the cube

    Q_n = { x : |e_k . (x - n)| <= a/2 for k = 0, 1, 2 }

which contains its junction ``J_n`` and whose six faces *are* the junction's cap
planes. So for every cell a support primitive reaches:

1. each primitive is clipped to ``Q_n`` analytically
   (:func:`latticegen2.support.polyhedra.clip_to_cell`);
2. the clipped pieces are fused with ``J_n`` — a handful of small convex
   operands, the same class of call as the junction template's own fuse;
3. the result is intersected with the input body, one object operand at a time
   (docs/algorithm.md §7), unless nothing of the body's surface is near;
4. the faces lying in the cube's planes are tagged exactly as a trimmed
   junction's caps are.

What comes out is an ordinary :class:`latticegen2.boundary.BoundaryPiece` that
replaces the node's lattice-only one. Everything downstream — interface
resolution, the sew, the interior build, unification, validation — treats a
support cross-section in a cap plane the way it has always treated a trimmed
cap: both sides must present agreeing material, and then it is sewn across.

**Coverage is proven here, against the real geometry.** Every overhang sample
belonging to the cell is moved a hair into the body and classified against the
piece the kernel actually produced. A miss fails the run by name; it is never
left for someone to find in the part.
"""

from __future__ import annotations

import os
import pickle
import time

import numpy as np
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCP.BRepClass3d import BRepClass3d, BRepClass3d_SolidClassifier
from OCP.BRepLib import BRepLib
from OCP.BRepTools import BRepTools
from OCP.TopAbs import TopAbs_ShapeEnum, TopAbs_State
from OCP.TopExp import TopExp
from OCP.TopoDS import TopoDS, TopoDS_Shape
from OCP.TopTools import TopTools_IndexedMapOfShape, TopTools_ListOfShape
from OCP.gp import gp_Pnt

from .. import occ
from ..boundary import _piece_from, _piece_vertices, _remove_pinholes, _worker_probe
from ..classify import _point_triangle_dist
from ..errors import ProcessingError
from ..junction import build_template, is_cap_plane_face
from ..lattice import (
    HALF_STRUTS,
    LatticeParams,
    half_strut_offset,
    lattice_params,
    nodes,
    profile_vertices,
)
from ..parallel import WorkerPool
from ..parallel import compound_children as _compound_children
from ..parallel import read_brep as _read_brep
from . import polyhedra as ph

NodeKey = tuple[int, int, int]

MIN_PIECE_VOLUME = 1e-6
"""mm³. A clipped piece smaller than this is not built, and the cell on the
other side of the face it pokes through is cut back so that it presents no
cross-section there either (`polyhedra.clip_to_cell_matched`).

It was 1e-9, on the assumption that both cells sharing such a piece drop it
alike. They do not: the near cell keeps a real cross-section. And 1e-9 is below
what the kernel holds on to -- a 3.6e-09 mm³ share of a capital in the cell of
node (206, -283, -264) of `TD_HX_rehearsal_test` was nominally kept and then
lost in the fuse, leaving its neighbour presenting 5e-6 mm² to nothing. 1e-6 is
a cube a hundredth of a millimetre on a side."""

INWARD_STEP = 2e-3
"""Millimetres a sample is moved into the body before it is classified."""

SOLID_SEW_TOL = 1e-7


def polyhedron_solid(verts: np.ndarray, faces: list) -> TopoDS_Shape:
    """A convex polyhedron as an OCCT solid with outward-facing faces."""
    shell = occ.sew([occ.polygon_face(verts[loop]) for loop in faces],
                    SOLID_SEW_TOL, cutting=False)
    shells = occ.shells(shell)
    if len(shells) != 1:
        raise ProcessingError(
            f"A clipped support piece did not close into one shell "
            f"({len(shells)} found)."
        )
    solid = occ.make_solid(shells[0])
    BRepLib.OrientClosedSolid_s(TopoDS.Solid_s(solid))
    return solid


def _boolean(algo, objects: list, tools: list, what: str) -> TopoDS_Shape:
    args = TopTools_ListOfShape()
    for s in objects:
        args.Append(s)
    tl = TopTools_ListOfShape()
    for s in tools:
        tl.Append(s)
    algo.SetArguments(args)
    algo.SetTools(tl)
    algo.Build()
    if not algo.IsDone():
        raise ProcessingError(f"{what} did not complete.")
    return algo


TRIM_VOLUME_REL_TOL = 1e-3
TRIM_VOLUME_ABS_TOL = 1e-7
"""How far ``|common| + |cut|`` may differ from the operand's own volume, as a
fraction of it plus an absolute floor in mm³, before the two are treated as
contradicting each other.

Calibrated on both sides, on `TD_HX_rehearsal_test` at ``cc=7, t=1.4`` turned
over. Sound trims against that part's B-spline wall disagree by 1.8e-05 to
2.9e-05 of the operand -- the two booleans each fit their own intersection
curve -- and by under 1e-8 against its analytic faces. The contradictions this
exists for are whole regions going missing: 2.3e-01 to 3.9e-01. The bar was
1e-5 for one run, which put every cell on the spline wall on the wrong side of
it: one cell was halved 124 times chasing noise no halving removes, and the
stage went from 6 minutes to 31."""

BISECT_DEPTH = 6
"""How many times a convex operand is halved while its trim contradicts itself."""


def checked_trim(shape: TopoDS_Shape, body: TopoDS_Shape, where, clock: dict):
    """``(solids, consistent)``: ``shape`` inside ``body``, and whether the
    kernel's two answers about it add up.

    **An intersection can come back wrong while reporting success**
    (docs/algorithm.md §7.2), and for a support that is a hole in the coverage
    rather than a strut reaching outside. Measured on `TD_HX_rehearsal_test` at
    ``cc=7, t=1.4`` turned over, in the cell of node (227, -273, -298): a valid
    65.9 mm³ fused cell intersects the body to **0.021 mm³**, and one 12.18 mm³
    clipped capital to **nothing at all**, with every one of swapped operands,
    fuzzy values from 1e-7 to 1e-3, ``CheckInverted`` off, a local block of the
    body and ``BOPAlgo_CellsBuilder`` returning the same. Point classification
    puts about 4.9 mm³ of that capital inside the body.

    What does expose it is asking the complementary question. ``Cut`` of the
    same capital returns 7.43 mm³, so 4.74 mm³ of a 12.18 mm³ operand is in
    neither answer. Agreement between the two is not proof — §7.2 records a
    case where both were wrong the same way, which is what the coverage proof
    and the outside probe are still for — but a contradiction is conclusive,
    and it costs one more boolean.
    """
    t1 = time.perf_counter()
    try:
        inside = occ.solids(_boolean(
            BRepAlgoAPI_Common(), [shape], [body],
            f"Intersecting the cell at {where} with the input body").Shape())
        outside = _boolean(
            BRepAlgoAPI_Cut(), [shape], [body],
            f"Cutting the cell at {where} against the input body").Shape()
    finally:
        clock["trim"] = clock.get("trim", 0.0) + time.perf_counter() - t1
    whole = occ.volume(shape)
    got = sum(occ.volume(x) for x in inside) + occ.volume(outside)
    ok = abs(got - whole) <= TRIM_VOLUME_REL_TOL * abs(whole) + TRIM_VOLUME_ABS_TOL
    return inside, ok


def _halves(verts: np.ndarray, faces: list, up: np.ndarray | None = None) -> list:
    """A convex polyhedron cut in two across its longest extent.

    The cutting plane is tilted slightly off the axis so it is not parallel to
    anything the lattice or the supports are built from.

    **With ``up`` — the build direction — the plane is vertical**, across the
    longer of the polyhedron's two horizontal extents. A cutting plane is not
    always buried: where the kernel returns less of one half than of the other,
    part of it is left as a face of the result, and that face has to be
    self-supporting like every other support surface (docs/algorithm.md §14.2).
    A vertical one always is. The first version cut across the longest extent
    whatever its direction, and on `TD_HX_rehearsal_test` with supports that
    left 0.5 mm² of a near-horizontal cut exposed in the cell of node
    (228, -273, -298) — found by the written file's own check, at 89.13° from
    the build direction, which is that plane's tilt to five figures.
    """
    if up is None:
        extent = verts.max(axis=0) - verts.min(axis=0)
        normal = np.array([0.0137, 0.0071, 0.0113])
        normal[int(np.argmax(extent))] += 1.0
    else:
        up = np.asarray(up, dtype=float)
        seed = np.eye(3)[int(np.argmin(np.abs(up)))]
        hx = np.cross(up, seed)
        hx /= np.linalg.norm(hx)
        hy = np.cross(up, hx)
        span = [float(np.ptp(verts @ h)) for h in (hx, hy)]
        major, minor = (hx, hy) if span[0] >= span[1] else (hy, hx)
        normal = major + 0.0137 * minor          # still exactly perpendicular to `up`
    normal /= np.linalg.norm(normal)
    offset = float(normal @ verts.mean(axis=0))
    out = []
    for sign in (1.0, -1.0):
        got = ph.clip_halfspace(verts, faces, sign * normal, sign * offset)
        if got is None or len(got[1]) < 4:
            continue
        got = (got[0], ph.orient_outward(*got))
        if ph.volume(*got) > MIN_PIECE_VOLUME:
            out.append(got)
    return out


BLOCK_PAD = 0.137
"""How far, as a fraction of the polyhedron's extent, the local block reaches
beyond its bounding box on every side.

The block is axis-aligned on purpose. A block turned to a frame aligned with
nothing was tried first, on the argument that it could not be coplanar with a
machined body's faces, and in the cell of node (209, -266, -284) it intersected
the body to nothing where axis-aligned boxes of four sizes round the same
polyhedron all intersect correctly."""


OUTSIDE_TOL = 2e-3
"""Millimetres within which a vertex counts as on the body's surface."""

JUDGE_SAMPLES = 24
"""Sample points per operand the judge asks about."""


class Judge:
    """What is known about a trim without asking the kernel.

    Three intersections of one cell of `TD_HX_rehearsal_test` (node
    (209, -266, -284), ``cc=7, t=1.4``, turned over) show why the kernel's own
    cross-check is not enough. One 61.9 mm³ capital comes back **whole**, its
    cut empty, reaching 1.46 mm outside the body; another of 47.1 mm³ comes
    back **empty**, its cut whole, where point sampling puts 39 mm³ of it
    inside. Both pairs add up exactly. No overhang sample happens to lie in
    either, so the coverage proof says nothing.

    This asks the classification mesh instead — ray parity and a distance, the
    two halves of :class:`latticegen2.classify.OutsideProbe` and for the same
    reasons:

    * **outside**: no vertex of the result may lie outside the body. The mesh
      proposes and the body's own point classifier disposes: the probe's margin
      is twice the mesh's deviation, 0.68 mm on that part, and a capital in the
      cell of node (228, -273, -298) came back reaching **0.62 mm** outside —
      under the margin, so through the probe, into a written file whose
      verification then found it. So every vertex ray parity puts outside, or
      within the margin of the mesh, is asked of
      ``BRepClass3d_SolidClassifier`` at :data:`OUTSIDE_TOL`. That is a point
      query, not a boolean, and it is what every Monte Carlo figure in this
      module was measured with;
    * **missing**: a point inside an operand and inside the body must be inside
      the result. Parity within the mesh's own deviation of the surface means
      nothing, so a missed point counts only when it is further than the margin
      from the mesh — which is computed for the misses alone, never for all.
    """

    def __init__(self, probe, body: TopoDS_Shape | None = None):
        self.probe = probe
        self._classifier = BRepClass3d_SolidClassifier(body) if body is not None else None

    def outside(self, solids: list) -> float:
        """How far the result reaches outside the body, in mm; 0.0 if it does not."""
        if not solids:
            return 0.0
        points = np.unique(np.round(
            _piece_vertices([(occ.faces(s), None, 0) for s in solids]), 9), axis=0)
        if self._classifier is None:
            return self.probe.worst_outside(points)
        a, b, c = self.probe._tris
        suspect = points[~self.probe._inside(points)]
        worst = 0.0
        for p in suspect:
            distance = float(_point_triangle_dist(p, a, b, c).min())
            if distance <= OUTSIDE_TOL:
                continue                     # on the surface, to the mesh itself
            self._classifier.Perform(gp_Pnt(*p), OUTSIDE_TOL)
            if self._classifier.State() == TopAbs_State.TopAbs_OUT:
                worst = max(worst, distance)
        return worst

    def missing(self, polys: list, solids: list) -> int:
        rng = np.random.default_rng(20261010)
        found = []
        for verts, faces in polys:
            lo, hi = verts.min(axis=0), verts.max(axis=0)
            if float(np.prod(hi - lo)) <= MIN_PIECE_VOLUME:
                continue
            normals, offsets = ph.face_planes(verts, faces)
            q = lo + rng.random((8 * JUDGE_SAMPLES, 3)) * (hi - lo)
            q = q[((q @ normals.T) <= offsets - WITNESS_DEPTH).all(axis=1)]
            found.append(q[:JUDGE_SAMPLES])
        if not found:
            return 0
        q = np.vstack(found)
        q = q[self.probe._inside(q)]
        if not len(q):
            return 0
        lost = q[~_covered(solids, q)] if solids else q
        if not len(lost):
            return 0
        a, b, c = self.probe._tris
        return int(sum(float(_point_triangle_dist(p, a, b, c).min()) > self.probe.margin
                       for p in lost))

    def refutes(self, polys: list, solids: list) -> bool:
        return self.outside(solids) > 0.0 or self.missing(polys, solids) > 0


def _halfspace(normal: np.ndarray, offset: float, inside: bool):
    """The half-space ``normal . x <= offset`` as a solid, or its complement."""
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeHalfSpace
    from OCP.gp import gp_Dir, gp_Pln

    on = normal * offset
    face = BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(*on), gp_Dir(*normal))).Face()
    ref = on - normal if inside else on + normal
    return BRepPrimAPI_MakeHalfSpace(face, gp_Pnt(*ref)).Solid()


def trim_by_planes(poly, body: TopoDS_Shape, where, clock: dict):
    """``(solids, consistent)``: a convex polyhedron inside ``body``, reached
    without ever intersecting the polyhedron with the body.

    The polyhedron is the intersection of its face half-spaces, so the same
    region is a block of the body cut by one plane after another. Every
    operation is then either a box against the body -- well conditioned at
    exactly the places a slanted operand is not (docs/algorithm.md §7.2) -- or
    a single plane against what is left.

    Measured on the capital :func:`trim_convex` describes, the one whose
    intersection and cut were both wrong and agreed: 10.085 mm³ in 0.17 s,
    against 9.83 by point sampling, covering all 140 of the samples inside it,
    where bisecting it to :data:`BISECT_DEPTH` still left two uncovered.

    Each cut is checked the way :func:`checked_trim` checks: the two sides of
    the plane must add up to what was cut.
    """
    verts, faces = poly
    t1 = time.perf_counter()
    try:
        lo, hi = verts.min(axis=0), verts.max(axis=0)
        pad = BLOCK_PAD * float((hi - lo).max())
        lo, hi = lo - pad, hi + pad
        corners = np.array([[x, y, z] for x in (lo[0], hi[0]) for y in (lo[1], hi[1])
                            for z in (lo[2], hi[2])])
        cube = [[0, 1, 3, 2], [4, 6, 7, 5], [0, 4, 5, 1], [2, 3, 7, 6],
                [0, 2, 6, 4], [1, 5, 7, 3]]
        block = polyhedron_solid(corners, ph.orient_outward(corners, cube))
    finally:
        clock["trim"] = clock.get("trim", 0.0) + time.perf_counter() - t1
    current, ok = checked_trim(block, body, where, clock)
    normals, offsets = ph.face_planes(verts, faces)
    t1 = time.perf_counter()
    try:
        for normal, offset in zip(normals, offsets):
            if not current:
                break
            keep_side = _halfspace(normal, float(offset), True)
            drop_side = _halfspace(normal, float(offset), False)
            kept = []
            for solid in current:
                mine = occ.solids(_boolean(
                    BRepAlgoAPI_Common(), [solid], [keep_side],
                    f"Cutting the local block at {where} by a support plane").Shape())
                rest = _boolean(
                    BRepAlgoAPI_Common(), [solid], [drop_side],
                    f"Cutting the local block at {where} by a support plane").Shape()
                whole = occ.volume(solid)
                got = sum(occ.volume(x) for x in mine) + occ.volume(rest)
                ok = ok and (abs(got - whole)
                             <= TRIM_VOLUME_REL_TOL * abs(whole) + TRIM_VOLUME_ABS_TOL)
                kept.extend(x for x in mine if occ.volume(x) > MIN_PIECE_VOLUME)
            current = kept
    finally:
        clock["trim"] = clock.get("trim", 0.0) + time.perf_counter() - t1
    return current, ok


WITNESS_DEPTH = 1e-6
"""How far inside a polyhedron, in mm, a point must lie to be a witness of it."""


def _witnesses_in(poly, witnesses) -> np.ndarray:
    """Those of ``witnesses`` lying strictly inside the convex ``poly``."""
    if witnesses is None or not len(witnesses):
        return np.empty((0, 3))
    normals, offsets = ph.face_planes(*poly)
    return witnesses[((witnesses @ normals.T) <= offsets - WITNESS_DEPTH).all(axis=1)]


def trim_convex(poly, body: TopoDS_Shape, where, clock: dict, depth: int = 0,
                witnesses: np.ndarray | None = None, judge: "Judge | None" = None,
                up: np.ndarray | None = None) -> list:
    """A convex polyhedron inside ``body``, bisected while the kernel
    contradicts itself about it.

    The same idea as §7.2's repairs — give the kernel a smaller, simpler
    operand — applied recursively and decided by :func:`checked_trim` rather
    than by which operand happened to work. On the capital above: the whole is
    inconsistent, one half is consistent at 1.71 mm³, the other is halved
    again, and three levels down everything adds up, for 4.71 mm³ against a
    point-sampled 4.87. A leaf still inconsistent at :data:`BISECT_DEPTH` is
    kept as the kernel returned it; whether that left a gap is the coverage
    proof's question, and it fails the run by name if so.

    **Two answers that add up can both be wrong**, which §7.2 already records
    for a junction and which happens here too. In the cell of node
    (240, -277, -280) of the same run a 32.95 mm³ clipped capital intersects
    the body to nothing *and* cuts against it to all of itself, perfectly
    consistent, with 140 overhang samples inside it. ``witnesses`` is the
    evidence that needs no kernel: points known to be inside the body — the
    cell's own samples, moved a hair inward — so any of them inside this
    polyhedron and outside its trim is a contradiction as conclusive as a
    volume that does not add up, and is treated the same way. ``judge`` is the
    same kind of evidence for the operands no sample lies in (:class:`Judge`).
    """
    solid = polyhedron_solid(*poly)
    inside, ok = checked_trim(solid, body, where, clock)
    mine = _witnesses_in(poly, witnesses)

    def refuted(solids):
        if len(mine) and not _covered(solids, mine).all():
            clock["witness_refusals"] = clock.get("witness_refusals", 0) + 1
            return True
        if judge is not None and judge.refutes([poly], solids):
            clock["judge_refusals"] = clock.get("judge_refusals", 0) + 1
            return True
        return False

    if ok and not refuted(inside):
        return inside
    # The same region by another route before the operand is made smaller:
    # plane by plane, which never intersects the polyhedron with the body.
    planed, ok = trim_by_planes(poly, body, where, clock)
    if ok and not refuted(planed):
        clock["by_planes"] = clock.get("by_planes", 0) + 1
        return planed
    if depth >= BISECT_DEPTH:
        # Nothing has produced a trim that stands up, and what is left is a
        # sixty-fourth of the operand. Whether it is *missing* is the coverage
        # proof's question; whether it reaches *outside* is not negotiable
        # (specification.md §1), so of the two answers the one that stays
        # inside is kept, and if neither does the leaf is left out -- the same
        # decision §7.2 takes for a junction no intersection keeps inside.
        if judge is None:
            return inside
        for candidate in (inside, planed):
            if judge.outside(candidate) == 0.0:
                return candidate
        clock["leaves_dropped"] = clock.get("leaves_dropped", 0) + 1
        return []
    clock["bisections"] = clock.get("bisections", 0) + 1
    out = []
    for half in _halves(*poly, up=up):
        out.extend(trim_convex(half, body, where, clock, depth + 1, mine, judge, up))
    return out


def _half_strut_polys(lp: LatticeParams, node_pos: np.ndarray) -> list:
    """The junction's six half-struts as convex polyhedra (docs/algorithm.md §3.2)."""
    out = []
    for h in range(6):
        k = HALF_STRUTS[h][0]
        base = profile_vertices(lp, node_pos, k)
        verts = np.vstack([base, base + half_strut_offset(lp, h)])
        faces = [[0, 1, 2, 3], [4, 5, 6, 7]] + [
            [i, (i + 1) % 4, 4 + (i + 1) % 4, 4 + i] for i in range(4)]
        out.append((verts, ph.orient_outward(verts, faces)))
    return out


def build_cell(
    lp: LatticeParams,
    tpl,
    node_pos: np.ndarray,
    clipped: list,
    body: TopoDS_Shape | None,
    with_junction: bool,
    per_operand: bool = False,
    clock: dict | None = None,
    witnesses: np.ndarray | None = None,
    judge: "Judge | None" = None,
) -> list:
    """The solids of one lattice cell: junction plus support, inside the body.

    ``clipped`` is ``(verts, faces)`` per primitive already cut to the cell.
    ``body`` is ``None`` when nothing of its surface is near the cell, and the
    intersection is then skipped.

    The usual route is one fuse and one intersection of the result, verified by
    :func:`checked_trim`. Where the kernel contradicts itself about that
    intersection — or when the caller asks for it with ``per_operand``, after
    the usual route left a sample uncovered — every operand is trimmed on its
    own with :func:`trim_convex` and the trimmed pieces are fused instead.

    ``witnesses`` are points known to lie inside the body, and ``judge`` the
    kernel-free check of :class:`Judge`; both are handed to :func:`trim_convex`
    on the operand-by-operand route.

    ``clock``, if given, accumulates seconds under ``"solids"``, ``"fuse"`` and
    ``"trim"``, and counts under ``"fallback_cells"`` and ``"bisections"``.
    """
    clock = clock if clock is not None else {}
    where = tuple(np.round(node_pos, 4))

    def fuse(shapes):
        if len(shapes) == 1:
            return occ.solids(shapes[0])
        t1 = time.perf_counter()
        try:
            return _fuse(shapes)
        finally:
            clock["fuse"] = clock.get("fuse", 0.0) + time.perf_counter() - t1

    def _fuse(shapes):
        algo = _boolean(BRepAlgoAPI_Fuse(), shapes[:1], shapes[1:],
                        f"Fusing the support pieces of the cell at {where}")
        # Merge the coplanar faces the operands leave side by side. Several
        # clipped pieces present overlapping polygons in the same cell plane,
        # and how the kernel happens to partition that overlap is not the same
        # from the two cells sharing the plane -- while the *outline* of the
        # merged region is. The interface check compares rings edge for edge
        # (latticegen2.weld.unweldable), so it must be handed outlines.
        try:
            algo.SimplifyResult()
        except Exception:                                      # noqa: BLE001
            pass      # a size optimisation that refuses is not a failure (§11)
        return occ.solids(algo.Shape())

    if not clipped and not with_junction:
        return []

    if not per_operand:
        t0 = time.perf_counter()
        operands = []
        if with_junction:
            operands.append(tpl.solid.Moved(occ.translation(node_pos)))
        operands.extend(polyhedron_solid(v, f) for v, f in clipped)
        clock["solids"] = clock.get("solids", 0.0) + time.perf_counter() - t0
        fused = fuse(operands)
        if body is None:
            return fused
        out, consistent = [], True
        for solid in fused:
            inside, ok = checked_trim(solid, body, where, clock)
            out.extend(inside)
            consistent = consistent and ok
        if consistent:
            return out

    # Every operand on its own. The junction first as the one solid it is; only
    # if the kernel contradicts itself about that too, as its six convex
    # half-struts, which is §7.2's own repair.
    clock["fallback_cells"] = clock.get("fallback_cells", 0) + 1
    if body is None:
        return fuse([tpl.solid.Moved(occ.translation(node_pos))] * with_junction
                    + [polyhedron_solid(v, f) for v, f in clipped])
    trimmed = []
    if with_junction:
        struts = _half_strut_polys(lp, node_pos)
        inside, ok = checked_trim(
            tpl.solid.Moved(occ.translation(node_pos)), body, where, clock)
        if ok and not (judge is not None and judge.refutes(struts, inside)):
            trimmed.extend(inside)
        else:
            for poly in struts:
                trimmed.extend(trim_convex(poly, body, where, clock, judge=judge,
                                           up=lp.b))
    for poly in clipped:
        trimmed.extend(trim_convex(poly, body, where, clock,
                                   witnesses=witnesses, judge=judge, up=lp.b))
    return fuse(trimmed) if trimmed else []


def _covered(solids: list, points: np.ndarray) -> np.ndarray:
    """Which ``points`` lie in or on one of ``solids``."""
    done = np.zeros(len(points), dtype=bool)
    for solid in solids:
        if done.all():
            break
        classifier = BRepClass3d_SolidClassifier(solid)
        for i in np.nonzero(~done)[0]:
            classifier.Perform(gp_Pnt(*points[i]), 1e-6)
            if classifier.State() in (TopAbs_State.TopAbs_IN, TopAbs_State.TopAbs_ON):
                done[i] = True
    return done


_CTX_CACHE: dict = {}


def _context(ctx_path: str) -> dict:
    """The staged support context, read and prepared once per worker process."""
    got = _CTX_CACHE.get(ctx_path)
    if got is None:
        with open(ctx_path, "rb") as fh:
            got = pickle.load(fh)
        got["lp"] = lattice_params(got["cc"], got["t"], got["orient"])
        got["tpl"] = build_template(got["lp"])
        got["body"] = _read_brep(got["body_path"])
        _CTX_CACHE.clear()
        _CTX_CACHE[ctx_path] = got
    return got


def _worker_support(job):
    """Build one batch of supported cells in a worker process.

    Same small-IPC discipline as :func:`latticegen2.boundary._worker_trim`: the
    geometry comes back as a ``.brep`` of per-piece face compounds beside
    parallel metadata, and only paths and plain data cross the process
    boundary.
    """
    ctx_path, batch, out_path = job
    ctx = _context(ctx_path)
    lp, tpl, body = ctx["lp"], ctx["tpl"], ctx["body"]
    prims = ctx["prims"]
    probe = _worker_probe(ctx["mesh_path"], ctx["margin"])
    judge = Judge(probe, body) if probe is not None else None

    bundles: list[TopoDS_Shape] = []
    meta = []
    n_retried = 0
    n_points = 0
    n_voids = 0
    n_clipped = n_operands = 0
    worst_outside = 0.0
    clock: dict = {}
    untouched: list = []
    for node in batch:
        node = tuple(int(x) for x in node)
        info = ctx["cells"][node]
        pos = nodes(lp, np.array([node], dtype=np.int64))[0]
        clipped = []
        t0 = time.perf_counter()
        for pid in info["prims"]:
            got = ph.clip_to_cell_matched(prims[pid][0], prims[pid][1], lp, pos,
                                          MIN_PIECE_VOLUME)
            if got is not None:
                clipped.append(got)
        n_clipped += len(clipped)
        clipped = ph.drop_contained(clipped)
        n_operands += len(clipped)
        clock["clip"] = clock.get("clip", 0.0) + time.perf_counter() - t0
        points = info["points"]
        if not clipped:
            # The primitive's bounding box reached this cell and the primitive
            # itself did not. The cell keeps whatever the boundary trim gave it
            # -- rebuilding a junction no support touches would only replace a
            # piece by its equal, re-described.
            if len(points):
                p = points[0]
                raise ProcessingError(
                    f"{len(points)} overhang sample(s) lie in the lattice cell "
                    f"of node {node}, the first at [{p[0]:.3f}, {p[1]:.3f}, "
                    f"{p[2]:.3f}], and no support reaches that cell."
                )
            untouched.append(node)
            continue
        local = body if info["trim"] else None
        solids = build_cell(lp, tpl, pos, clipped, local, info["junction"],
                            clock=clock)
        t0 = time.perf_counter()
        ok = (not len(points)) or bool(_covered(solids, points).all())
        clock["cover"] = clock.get("cover", 0.0) + time.perf_counter() - t0
        operands = clipped
        if judge is not None and local is not None:
            # What no sample speaks for: an operand returned untrimmed, or lost.
            t0 = time.perf_counter()
            if info["junction"]:
                operands = clipped + _half_strut_polys(lp, pos)
            if ok and judge.refutes(operands, solids):
                ok = False
                clock["judged_cells"] = clock.get("judged_cells", 0) + 1
            clock["judge"] = clock.get("judge", 0.0) + time.perf_counter() - t0
        if not ok:
            n_retried += 1
            solids = build_cell(lp, tpl, pos, clipped, local, info["junction"],
                                per_operand=True, clock=clock,
                                witnesses=points, judge=judge)
            if judge is not None and local is not None:
                clock["cells_missing_material"] = (
                    clock.get("cells_missing_material", 0)
                    + int(judge.missing(operands, solids) > 0))
            miss = ~_covered(solids, points)
            if miss.any():
                p = points[np.nonzero(miss)[0][0]]
                raise ProcessingError(
                    f"Support generation left {int(miss.sum())} overhang "
                    f"sample(s) uncovered in the lattice cell of node {node}, "
                    f"the first at [{p[0]:.3f}, {p[1]:.3f}, {p[2]:.3f}], with "
                    f"every operand trimmed on its own as well. The output would have an unsupported "
                    f"overhang there, so the run stops instead of writing it."
                )
        n_points += len(points)

        t0 = time.perf_counter()
        pieces = []
        for solid in solids:
            cleaned, _removed = _remove_pinholes(pos, solid)

            # Outer shell only. Overlapping capitals can seal off a pocket of
            # empty space between them; inside one cell that is an inner shell
            # of this solid, and keeping it would ship an enclosed void no
            # powder can leave. Dropping the shell fills it. A pocket spanning
            # several cells is filled at assembly instead
            # (latticegen2.support.build.fill_voids).
            shells = occ.shells(cleaned)
            if len(shells) > 1:
                outer = BRepClass3d.OuterShell_s(TopoDS.Solid_s(cleaned))
                n_voids += len(shells) - 1
                cleaned = occ.make_solid(outer)
            faces = occ.faces(cleaned)
            tags = [is_cap_plane_face(lp, f, pos) for f in faces]
            tags = [-1 if h is None else h for h in tags]
            pieces.append((faces, tags, occ.volume(cleaned)))
        if judge is not None and local is not None and solids:
            worst_outside = max(worst_outside, judge.outside(solids))
        for faces, tags, vol in pieces:
            bundles.append(occ.compound(faces))
            meta.append((node, tags, vol, tuple(occ.tolerance_feature_ratio(faces))))
        clock["finish"] = clock.get("finish", 0.0) + time.perf_counter() - t0

    from ..runlog import peak_rss_bytes

    extra = (n_retried, n_points, worst_outside, n_voids, clock,
             n_clipped, n_operands, untouched)
    if bundles:
        BRepTools.Write_s(occ.compound(bundles), out_path)
        return out_path, meta, extra, peak_rss_bytes()
    return None, meta, extra, peak_rss_bytes()


def cells_of(lp: LatticeParams, verts: np.ndarray) -> list:
    """Every lattice cell a primitive's bounding box reaches, in index space.

    A superset of the cells it actually meets — the exact clip decides that —
    but never a subset: a cell left out here would keep its lattice-only piece
    while its neighbour presented a support cross-section to it.
    """
    idx = np.linalg.solve(lp.B, verts.T).T
    lo = np.floor(idx.min(axis=0) + 0.5 - 1e-6).astype(np.int64)
    hi = np.floor(idx.max(axis=0) + 0.5 + 1e-6).astype(np.int64)
    return [(i, j, k)
            for i in range(lo[0], hi[0] + 1)
            for j in range(lo[1], hi[1] + 1)
            for k in range(lo[2], hi[2] + 1)]


def build_supports(
    lp: LatticeParams,
    primitives: list,
    samples,
    mesh,
    lattice_nodes: set,
    interior: set,
    body_path: str,
    mesh_path: str,
    outside_margin: float,
    tmpdir: str,
    pool: WorkerPool | None,
    workers: int,
    report=None,
):
    """Build every supported cell. Returns ``(pieces, supported_nodes, stats)``.

    ``lattice_nodes`` are the nodes that hold lattice after the boundary trim
    (interior nodes and nodes with a surviving piece): only those get their
    junction rebuilt alongside the support. A node the trim left empty stays
    empty of lattice, so a junction dropped for reaching outside the body
    (docs/algorithm.md §7.2) is not resurrected by a support passing through
    its cell.
    """
    from ..classify import PointInside, SpatialHash

    cells: dict[NodeKey, dict] = {}
    for pid, prim in enumerate(primitives):
        for node in cells_of(lp, prim.verts):
            cells.setdefault(node, {"prims": []})["prims"].append(pid)
    if not cells:
        return [], set(), {}

    keys = sorted(cells)
    idx = np.array(keys, dtype=np.int64)
    pos = nodes(lp, idx)
    # Half-extent of the cube's world bounding box, plus the mesh's own error.
    half = (lp.a / 2.0) * np.abs(lp.e).sum(axis=0) + mesh.deviation + INWARD_STEP
    hash_ = SpatialHash(mesh, min_cell=lp.a / 2.0)
    centre_inside = PointInside(mesh)(pos)

    # Each sample is checked in the cell its inward-moved copy falls in.
    # Only what exceeds the limit is proven covered; the margin band is there
    # to place supports generously, not to be proven (`OverhangSamples.required`).
    moved = (samples.points + INWARD_STEP * samples.inward)[samples.required]
    sources = samples.points[samples.required]
    owner = np.round(np.linalg.solve(lp.B, moved.T).T).astype(np.int64)
    by_cell: dict[NodeKey, list] = {}
    for i, n in enumerate(map(tuple, owner.tolist())):
        by_cell.setdefault(n, []).append(i)

    live = []
    for i, node in enumerate(keys):
        near = len(hash_.query(pos[i] - half, pos[i] + half)) > 0
        if not near and not centre_inside[i]:
            continue                      # the whole cube is outside the body
        info = cells[node]
        info["trim"] = bool(near)
        info["junction"] = node in lattice_nodes
        info["points"] = moved[by_cell.get(node, [])]
        live.append(node)

    orphans = [n for n in by_cell if n not in cells or "points" not in cells[n]]
    if orphans:
        n = orphans[0]
        p = sources[by_cell[n][0]]
        raise ProcessingError(
            f"{sum(len(by_cell[o]) for o in orphans)} overhang sample(s) fall in "
            f"lattice cells no support reaches, the first at "
            f"[{p[0]:.3f}, {p[1]:.3f}, {p[2]:.3f}] in the cell of node {n}."
        )

    ctx_path = os.path.join(tmpdir, "support_ctx.pkl")
    with open(ctx_path, "wb") as fh:
        pickle.dump({
            "cc": lp.cc, "t": lp.t, "orient": lp.orient,
            "body_path": body_path, "mesh_path": mesh_path,
            "margin": outside_margin,
            "prims": [(p.verts, p.faces) for p in primitives],
            "cells": {n: cells[n] for n in live},
        }, fh, protocol=pickle.HIGHEST_PROTOCOL)

    parallel = pool is not None and pool.active and workers > 1
    n_batches = max(1, min(len(live), workers * 4)) if parallel else 1
    bounds = np.linspace(0, len(live), n_batches + 1).astype(int)
    batches = [live[bounds[i]:bounds[i + 1]] for i in range(n_batches)
               if bounds[i + 1] > bounds[i]]
    jobs = [(ctx_path, batch, os.path.join(tmpdir, f"support_{bi}.brep"))
            for bi, batch in enumerate(batches)]

    def tick(done: int, total: int) -> None:
        if report is not None:
            report("building supported cells", done, total)

    if parallel:
        results, max_rss = pool.run(_worker_support, jobs, on_result=tick)
    else:
        results, max_rss = [], 0
        for i, job in enumerate(jobs):
            results.append(_worker_support(job))
            tick(i + 1, len(jobs))

    pieces = []
    n_retried = n_points = n_voids = n_clipped = n_operands = 0
    worst_outside = 0.0
    phases: dict = {}
    untouched: set = set()
    for path, meta, extra, _rss in results:
        n_retried += extra[0]
        n_points += extra[1]
        worst_outside = max(worst_outside, extra[2])
        n_voids += extra[3]
        n_clipped += extra[5]
        n_operands += extra[6]
        untouched.update(tuple(n) for n in extra[7])
        for key, seconds in extra[4].items():
            phases[key] = phases.get(key, 0) + seconds
        if path is None:
            continue
        children = _compound_children(_read_brep(path))
        if len(children) != len(meta):
            raise ProcessingError(
                f"Support worker result mismatch in {path}: {len(children)} "
                f"pieces vs {len(meta)} metadata records."
            )
        for bundle, (node, tags, vol, tf) in zip(children, meta):
            piece = _piece_from(tuple(node), _compound_children(bundle), tags, vol,
                                tuple(tf))
            piece.support = True
            pieces.append(piece)
    if worst_outside > 0.0:
        raise ProcessingError(
            f"A supported lattice cell reaches {worst_outside:.4f} mm outside "
            f"the input body after trimming. specification.md §1 requires the "
            f"output to fit within it, so the run stops rather than writing it."
        )
    if n_points != len(moved):
        raise ProcessingError(
            f"Only {n_points} of {len(moved)} overhang samples were checked "
            f"for coverage; the rest belong to cells that produced nothing."
        )
    live = [n for n in live if n not in untouched]
    stats = {
        "support_cells": len(live),
        "support_pieces": len(pieces),
        "support_cells_retried": n_retried,
        "support_voids_filled_in_cells": n_voids,
        "support_clipped_pieces": n_clipped,
        "support_fuse_operands": n_operands,
        # Cells whose one-shot intersection contradicted itself and were
        # rebuilt operand by operand, and how many operands had to be halved.
        "support_cells_per_operand": int(phases.get("fallback_cells", 0)),
        "support_operand_bisections": int(phases.get("bisections", 0)),
        "support_witness_refusals": int(phases.get("witness_refusals", 0)),
        "support_operands_by_planes": int(phases.get("by_planes", 0)),
        # Cells sent down that route by the kernel-free judge alone, operands it
        # refused, and cells where it still finds material missing afterwards.
        "support_cells_judged": int(phases.get("judged_cells", 0)),
        "support_judge_refusals": int(phases.get("judge_refusals", 0)),
        "support_cells_missing_material": int(phases.get("cells_missing_material", 0)),
        "support_leaves_dropped": int(phases.get("leaves_dropped", 0)),
        # Summed over workers, so core-seconds rather than wall clock: what
        # says which part of a cell's build a proposal would actually help.
        "support_core_seconds": " ".join(
            f"{k} {phases[k]:.1f}" for k in
            ("clip", "solids", "fuse", "trim", "cover", "finish") if k in phases),
        "support_samples_verified": n_points,
        "support_cells_interior": sum(1 for n in live if n in interior),
        "peak_support_worker_rss": max_rss,
    }
    return pieces, set(live), stats


def fill_voids(shells: dict) -> tuple[dict, int, float, list]:
    """Drop enclosed voids from assembled shells. Returns ``(shells, n, volume)``.

    Overlapping capitals can seal off a pocket of empty space. A component's
    faces then form more than one closed surface: the outside of the material,
    and the wall of each pocket. One shell holding both is not a connected
    shell, which `BRepCheck_Analyzer` rejects, and the pocket itself is a void
    powder could never leave.

    Each component's faces are grouped by shared edges. The group enclosing the
    most volume is the outside; every other group must enclose *negative*
    volume as oriented -- its normals point out of the material, into the
    pocket -- and is dropped, which fills the pocket with material. A second
    group of positive volume is a separate body, and it is kept as one.

    **The junction graph does not prove such a body cannot exist**, though this
    function once failed the run on that reading. The graph joins every piece
    presenting a cap to every piece presenting its partner
    (`connect.build_components`), which over-connects by design: it can only
    keep a body, never delete one. A bare junction has one piece per cap, so
    there the graph is exact. A supported cell can hold several, and where the
    body's wall cuts a sliver of support off from the rest, that sliver shares a
    cap *plane* with material it does not touch. Measured on
    `TD_HX_rehearsal_test` at ``cc=7, t=1.4`` with supports, turned over: one
    0.1455 mm³ fragment beside a 373,000 mm³ lattice. It is closed, it lies
    against the wall it was cut by, and specification.md §5 says a
    wall-attached support is not a floating body — so it ships, as its own
    solid, and is reported. Its closedness does not rest on the graph either
    way: every edge of every shell here has already been proven used twice.

    Returns ``(shells, n_filled, filled_volume, extra_volumes)``.
    """
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Shell

    from ..connect import UnionFind

    builder = BRep_Builder()
    out: dict = {}
    n_filled = 0
    filled_volume = 0.0
    extra: list = []
    next_key = max(shells, default=-1) + 1      # keys stay sortable integers
    for group, shell in shells.items():
        faces = occ.faces(shell)
        # Identity through OCCT's own indexed map, as `weld.shell_defects` does:
        # a Python-side hash of the underlying shape is not an identity.
        edges = TopTools_IndexedMapOfShape()
        TopExp.MapShapes_s(shell, TopAbs_ShapeEnum.TopAbs_EDGE, edges)
        owner = [-1] * (edges.Extent() + 1)
        uf = UnionFind(len(faces))
        for fi, face in enumerate(faces):
            for edge in occ._explore(face, TopAbs_ShapeEnum.TopAbs_EDGE):
                key = edges.FindIndex(edge)
                if owner[key] < 0:
                    owner[key] = fi
                elif owner[key] != fi:
                    uf.union(owner[key], fi)
        parts: dict = {}
        for fi in range(len(faces)):
            parts.setdefault(uf.find(fi), []).append(fi)
        if len(parts) == 1:
            out[group] = shell
            continue
        measured = []
        for members in parts.values():
            sub = TopoDS_Shell()
            builder.MakeShell(sub)
            for fi in members:
                builder.Add(sub, faces[fi])
            sub.Closed(True)
            measured.append((occ.volume(occ.make_solid(sub)), sub))
        measured.sort(key=lambda m: m[0], reverse=True)
        for vol, sub in measured[1:]:
            if vol > 0.0:
                out[next_key + len(extra)] = sub
                extra.append(vol)
                continue
            n_filled += 1
            filled_volume += -vol
        out[group] = measured[0][1]
    return out, n_filled, filled_volume, extra
