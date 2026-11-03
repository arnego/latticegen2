"""Geometry-validity checks for the E2E harness (specification.md §6.2).

These are QA tooling, not part of the generation pipeline, and they reason only
about the finished STEP file: they re-read it, re-tessellate it, and check the
resulting triangles. That keeps them independent of *how* the geometry was
built, which is the property that makes them worth running at all.

The strongest of them is :func:`brepcheck`, which asks OCCT to validate the
B-rep exactly rather than inferring validity from a tessellation.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import numpy as np
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCP.TopTools import TopTools_ListOfShape

from latticegen2 import occ
from latticegen2.classify import (
    SpatialHash,
    TriMesh,
    _segment_triangle_intersects,
    _weld,
    chordal_target,
)
from latticegen2.lattice import lattice_params


def mesh_of(path: str, cc: float, t: float, refine: float = 1.0) -> TriMesh:
    """Re-read a STEP file and tessellate it for the mesh-based checks.

    ``refine`` divides the chordal deflection. It exists for one caller: a
    crossing found at the default deflection is re-asked at a finer one before
    it is believed (:func:`tools.e2e.geometry_checks`).

    **Each solid is welded on its own.** Support bodies share faces with the
    lattice and with each other (specification.md §4.6), and welded together
    the edges of a shared face would be used four times — two closed bodies
    read as one non-manifold one. Kept apart, every body has to close by
    itself, which is the property the manifold check asks for.
    """
    shape = occ.read_step(path)
    lp = lattice_params(cc, t)
    occ.mesh_shape(shape, chordal_target(lp) / refine)
    verts, tris, owner, offset = [], [], [], 0
    for solid in occ.solids(shape):
        sv, sf, soff = [], [], 0
        for face in occ.faces(solid):
            got = occ.face_triangulation(face)
            if got is None:
                continue
            v, f = got
            sv.append(v)
            sf.append(f + soff)
            soff += len(v)
        if not sf:
            continue
        v, f = _weld(np.vstack(sv), np.vstack(sf))
        verts.append(v)
        tris.append(f + offset)
        owner.append(np.full(len(f), len(owner), dtype=np.int64))
        offset += len(v)
    if not tris:
        raise RuntimeError(f"No triangles produced for {path}")
    mesh = TriMesh(verts=np.vstack(verts), tris=np.vstack(tris))
    mesh.solid_of = np.concatenate(owner)
    return mesh


def manifold_check(mesh: TriMesh):
    """A closed manifold triangle mesh uses every edge exactly twice."""
    edges = np.vstack([mesh.tris[:, [0, 1]], mesh.tris[:, [1, 2]], mesh.tris[:, [2, 0]]])
    edges = np.sort(edges, axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    bad = int(np.count_nonzero(counts != 2))
    return bad == 0, bad


def _straddles(points: np.ndarray, plane_pt: np.ndarray, plane_n: np.ndarray, eps: float) -> bool:
    d = (points - plane_pt) @ plane_n
    return bool(np.any(d > eps) and np.any(d < -eps))


def triangles_properly_cross(t1: np.ndarray, t2: np.ndarray, eps: float = 1e-7) -> bool:
    """Do two triangles genuinely, transversally cross — not merely touch?

    Both triangles must straddle each other's plane before the edge-piercing
    test is attempted. Without that pre-check, two separate solids that touch
    along a coincident face report hundreds of spurious "intersections" purely
    because their independent triangulations are not vertex-aligned. Measured on
    two boxes sharing one exact face with zero volume overlap: 344 false
    positives without the pre-check, 0 with it.
    """
    n1 = np.cross(t1[1] - t1[0], t1[2] - t1[0])
    if not _straddles(t2, t1[0], n1, eps):
        return False
    n2 = np.cross(t2[1] - t2[0], t2[2] - t2[0])
    if not _straddles(t1, t2[0], n2, eps):
        return False
    for a, b in ((0, 1), (1, 2), (2, 0)):
        if _segment_triangle_intersects(t1[a], t1[b], t2[0:1], t2[1:2], t2[2:3])[0]:
            return True
        if _segment_triangle_intersects(t2[a], t2[b], t1[0:1], t1[1:2], t1[2:3])[0]:
            return True
    return False


def self_intersection_check(mesh: TriMesh, limit: int = 50):
    """Count genuinely crossing triangle pairs that share no vertex.

    Sweeps the spatial hash cell by cell rather than querying per triangle: any
    two triangles whose bounding boxes overlap share a cell, so this covers the
    same pairs while touching each cell's contents once. Adjacent triangles are
    *expected* to touch along shared edges, so pairs sharing any vertex are
    skipped.

    **Within one solid only, where the mesh says which solid a triangle is of**
    (``mesh.solid_of``, :func:`mesh_of`). A self-intersection is a property of a
    solid. Support bodies share faces with the lattice by design, and where the
    two meet on a curved ceiling their independent tessellations cross --
    measured on the supported ball: 316 pairs, every one between the lattice
    and the support, none within either. Whether two bodies overlap is asked
    exactly, of the B-rep, by :func:`supports_disjoint`.
    """
    solid_of = getattr(mesh, "solid_of", None)
    sh = SpatialHash(mesh)
    A, B, C = mesh.triangle_points
    bad: list[tuple[int, int]] = []
    seen: set[tuple[int, int]] = set()
    for group in sh.cell_groups():
        for a_pos in range(len(group)):
            ti = int(group[a_pos])
            for b_pos in range(a_pos + 1, len(group)):
                tj = int(group[b_pos])
                pair = (ti, tj) if ti < tj else (tj, ti)
                if pair in seen:
                    continue
                seen.add(pair)
                if solid_of is not None and solid_of[pair[0]] != solid_of[pair[1]]:
                    continue
                if set(mesh.tris[pair[0]]) & set(mesh.tris[pair[1]]):
                    continue
                t1 = np.array([A[pair[0]], B[pair[0]], C[pair[0]]])
                t2 = np.array([A[pair[1]], B[pair[1]], C[pair[1]]])
                if triangles_properly_cross(t1, t2):
                    bad.append(pair)
                    if len(bad) >= limit:
                        return False, bad
    return not bad, bad


def brepcheck(path: str):
    """OCCT's exact B-rep validity check on every solid in a STEP file."""
    shape = occ.read_step(path)
    sols = occ.solids(shape)
    invalid = [i for i, s in enumerate(sols) if not occ.is_valid(s)]
    return not invalid, len(sols), invalid


LOOSE_AREA_FRACTION_MAX = occ.LOOSE_AREA_FRACTION_MAX
"""The pipeline's own bar, re-used rather than restated.

Unlike ``CUT_MAX_FACES`` and ``SELF_INTERSECT_MAX_FACES`` beside it - limits on
what this harness can afford to measure, which are properly its own - this is a
bar on the *product*, and the pipeline enforces it too. Two copies of a bar are
two bars (docs/testing.md)."""


def curve_on_surface(path: str) -> dict:
    """Is the *shipped file* still describable, once its tolerances are gone?

    The strongest form of the export-truth check, because it is asked of the
    artefact rather than of the process that wrote it: the file is read back
    with OCCT's own reader, tolerances and all, and every edge/face pair is
    measured exactly (``latticegen2.occ.curve_on_surface_deviations``).

    The pipeline asks the same question of the solids in memory, which is where
    it can still name a junction. This asks it of what a reader reconstructs,
    which is what Solidworks and Catia will see - and the two are not the same
    shape, because STEP AP214 declares one tolerance for a whole file where the
    B-rep carries one per subshape (docs/algorithm.md S9).

    **Everything this returns is a diagnostic, not a verdict.** Three cheaper
    quantities were tried as gates and two false-positive on sound rehearsal
    solids (``latticegen2.occ.LOOSE_AREA_FRACTION_MAX``); what decides is
    whether the exported body still tessellates, which ``manifold_check`` above
    already asks of the whole file. These numbers say *why* a body is fragile
    and rank the candidates; they do not pass or fail anything.
    """
    shape = occ.read_step(path)
    worst, ratio, over, pairs = 0.0, 0.0, 0, 0
    fraction, loose_faces = 0.0, 0
    where, worst_area = [], 0.0
    for solid in occ.solids(shape):
        reading = occ.curve_on_surface_deviations(solid)
        worst = max(worst, reading.worst)
        over += reading.over_tolerance
        pairs += reading.pairs
        loose_faces += reading.loose_faces
        if reading.ratio > ratio:
            ratio, worst_area = reading.ratio, reading.face_area
        if reading.loose_area_fraction > fraction:
            fraction, where = reading.loose_area_fraction, list(reading.where)
    return {"worst_mm": worst, "ratio": ratio, "fraction": fraction,
            "loose_faces": loose_faces, "over": over, "pairs": pairs,
            "where": where, "face_area": worst_area}


def _cut_volume(a, b) -> float:
    alg = BRepAlgoAPI_Cut()
    la, lb = TopTools_ListOfShape(), TopTools_ListOfShape()
    la.Append(a)
    lb.Append(b)
    alg.SetArguments(la)
    alg.SetTools(lb)
    alg.Build()
    if not alg.IsDone():
        return float("nan")
    return occ.volume(alg.Shape())


def golden_sample_volume_diff(candidate_path: str, golden_path: str) -> float:
    """Larger of the two one-way difference volumes, in mm³.

    Near zero means the two files occupy essentially the same volume
    (specification.md §6.2). This is the exact test, and it is what should be
    used whenever it is affordable — but it is a general boolean between two
    complete lattices, so its cost grows with the square of the lattice's
    complexity. See :func:`golden_sample_agreement` for what to do when it is
    not affordable.
    """
    cand = occ.read_step(candidate_path)
    gold = occ.read_step(golden_path)
    return max(_cut_volume(cand, gold), _cut_volume(gold, cand))


def _diff_worker(queue, candidate_path: str, golden_path: str) -> None:
    try:
        queue.put(("ok", golden_sample_volume_diff(candidate_path, golden_path)))
    except Exception as exc:  # report rather than hang the parent
        queue.put(("error", f"{type(exc).__name__}: {exc}"))


def golden_sample_volume_diff_bounded(candidate_path: str, golden_path: str,
                                      timeout_s: float = 300.0):
    """The exact difference, or ``None`` if it does not finish in ``timeout_s``.

    The boolean runs in a child process because OCCT cannot be interrupted
    mid-call; the parent gives up on the child rather than on the operation.
    A ``None`` result means "not measured", never "measured as fine" — the
    caller must fall back to :func:`golden_sample_agreement` and say so.
    """
    import multiprocessing as mp

    ctx = mp.get_context("spawn")
    queue = ctx.Queue()
    proc = ctx.Process(target=_diff_worker, args=(queue, candidate_path, golden_path))
    proc.start()
    proc.join(timeout_s)
    if proc.is_alive():
        proc.terminate()
        proc.join(10)
        return None
    if queue.empty():
        return None
    kind, value = queue.get()
    if kind == "error":
        raise RuntimeError(value)
    return value


BBOX_TOL = 0.01
"""Millimetres. How far two representations of the same solid may differ in
extent before it means something. Two independently-tessellated, independently
-booleaned copies of one shape routinely differ by microns at an extreme; the
smallest feature this tool produces is `t >= 0.4` mm, so 0.01 mm is far below
anything geometrically meaningful and far above float noise."""

SURFACE_TIE_TOL = 0.05
"""Millimetres. A sampled point closer than this to *both* surfaces is on the
shared boundary, where an inside/outside test is at a tie and the two models can
legitimately disagree. Only disagreements further out than this from both
surfaces indicate an actual difference in shape."""


def golden_sample_agreement(candidate_path: str, golden_path: str, cc: float, t: float,
                            samples: int = 200_000, seed: int = 20260810) -> dict:
    """Sampled equivalence between two solids, for when the exact cut is too slow.

    Measured: the exact boolean difference between two ~100 k-face lattices had
    not finished after 20 minutes, where the same test on the ~10 k-face smoke
    output takes seconds. Rather than leave the larger scenario unchecked, this
    compares the two solids by:

    * total volume and bounding box (exact, cheap);
    * point membership at ``samples`` deterministic pseudo-random points drawn
      across the union bounding box, classified against each solid's own
      tessellation by the same three-ray parity test classification uses.

    This is **weaker than the boolean** and must be reported as such: it bounds
    the symmetric difference statistically rather than computing it. It is not a
    substitute for the exact test where the exact test can run — it is a way to
    say something quantified instead of nothing.

    Returns a dict of findings; ``disagreements`` of 0 with a matching volume is
    strong evidence the two solids are the same.
    """
    from latticegen2.classify import PointInside

    cand_mesh = mesh_of(candidate_path, cc, t)
    gold_mesh = mesh_of(golden_path, cc, t)
    cand_shape = occ.read_step(candidate_path)
    gold_shape = occ.read_step(golden_path)

    lo_c, hi_c = occ.bounding_box(cand_shape)
    lo_g, hi_g = occ.bounding_box(gold_shape)
    lo = np.minimum(lo_c, lo_g)
    hi = np.maximum(hi_c, hi_g)

    rng = np.random.default_rng(seed)
    pts = lo + rng.random((samples, 3)) * (hi - lo)
    in_cand = PointInside(cand_mesh)(pts)
    in_gold = PointInside(gold_mesh)(pts)
    disputed = np.nonzero(in_cand != in_gold)[0]

    # A disagreement only means something if the point is clear of both
    # surfaces. Points sitting on the shared boundary are ties, and with
    # hundreds of thousands of samples against a lattice's enormous surface
    # area, some always will be.
    real = 0
    if len(disputed):
        # Built once, outside the loop: each of these indexes every triangle in
        # the mesh, so building them per disputed point would make the cost of
        # this check scale with the number of disagreements times the mesh size.
        cand_hash = SpatialHash(cand_mesh, min_cell=_SEARCH_RADIUS)
        gold_hash = SpatialHash(gold_mesh, min_cell=_SEARCH_RADIUS)
        for i in disputed:
            p = pts[i]
            near_cand = _distance_to_mesh(cand_mesh, cand_hash, p)
            near_gold = _distance_to_mesh(gold_mesh, gold_hash, p)
            if min(near_cand, near_gold) > SURFACE_TIE_TOL:
                real += 1

    box_volume = float(np.prod(hi - lo))
    bbox_delta = float(max(np.abs(lo_c - lo_g).max(), np.abs(hi_c - hi_g).max()))
    return {
        "candidate_volume": occ.volume(cand_shape),
        "golden_volume": occ.volume(gold_shape),
        "volume_diff": abs(occ.volume(cand_shape) - occ.volume(gold_shape)),
        "bbox_delta": bbox_delta,
        "bbox_match": bbox_delta <= BBOX_TOL,
        "samples": samples,
        "disagreements": int(len(disputed)),
        "real_disagreements": real,
        "implied_symmetric_difference_mm3": box_volume * real / samples,
        "resolution_mm3": box_volume / samples,
    }


_SEARCH_RADIUS = 1.0
"""Millimetres. How far ``_distance_to_mesh`` looks before giving up and calling
the point far from the surface. Well above ``SURFACE_TIE_TOL``, which is the only
threshold the answer is ever compared against."""


def _distance_to_mesh(mesh: TriMesh, sh: SpatialHash, point: np.ndarray) -> float:
    """Distance from a point to the nearest triangle within ``_SEARCH_RADIUS`` mm.

    Takes a prebuilt index: the caller is expected to reuse one across points.
    """
    from latticegen2.classify import segment_triangle_dist

    candidates = sh.query(point - _SEARCH_RADIUS, point + _SEARCH_RADIUS)
    if len(candidates) == 0:
        return float("inf")
    A, B, C = mesh.triangle_points
    return float(segment_triangle_dist(point, point, A[candidates], B[candidates],
                                       C[candidates]).min())


CONTAINMENT_TOL = 2e-3
"""Millimetres. How far outside the input body a sampled point may sit before it
counts as material outside rather than as a point *on* the shared boundary.

A lattice is trimmed against the input, so a large part of its surface lies
exactly on the input's surface and every point there is a tie. What decides the
tie is how precisely "on" is recorded, and OCCT's own answer to that is the
tolerance it puts on those faces: 8.7e-04 to 1.5e-03 mm on the trimmed B-splines
of `TD_HX_rehearsal_test` (docs/algorithm.md §8, G12). This sits just above the
largest of them, so a point it calls outside is outside by more than the kernel's
own uncertainty about where the surface is — and still 200x below the smallest
strut the CLI accepts (`t >= 0.4` mm), so nothing geometrically meaningful can
hide under it."""

CONTAINMENT_DEFLECTION = 0.05
"""Millimetres. Chordal deviation for the tessellation whose vertices are
sampled by :func:`surface_points_outside`. Fine enough that a protrusion of
`CONTAINMENT_TOL` cannot slip between samples on a strut of side `t >= 0.4` mm."""

PROTRUSION_MIN_MM = 1e-6
"""Exact point-to-body distance below which a point outside is a tie, not material.

A point *on* the input's boundary cannot be material outside it, whatever the
classifier says about which side it is on. This is the bar on the exact
distance, and it is deliberately tight rather than tuned: the ties measured on
`SpiralTest.step` sit at 0.000000 mm (median 9.95e-10), and the genuine
protrusions docs/algorithm.md S7.2 was built to remove measured 0.239, 0.771 and
0.380 mm. Eight orders separate them, so nothing here rests on where in that gap
the bar sits -- which is the property docs/specification.md S11 keeps asking for
and rarely gets.

It is also the *only* bar that can be tight. The classifier sweep cannot be,
because it is asking a tolerance question about a point lying on a surface, and
its answer there is a coin flip an unclosed input seam can bias."""

#: Points measured exactly when the sweep cannot clear them. `BRepExtrema` is
#: ~milliseconds per point against a 13-face body, and a real containment
#: failure shows up in the first few -- the residuals docs/algorithm.md S7.2
#: found were whole junction faces, hundreds of points each, not singletons.
DISTANCE_PROBE_MAX = 200

CONTAINMENT_SWEEP = (1e-6, 1e-5, 1e-4, 1e-3, 2e-3)
"""Tolerances tried in order, so the result reports *how far* out the worst point
is rather than only whether it cleared one bar."""


def surface_points_outside(solid, body, deflection: float = CONTAINMENT_DEFLECTION) -> dict:
    """Boolean-free containment: are any of ``solid``'s surface points outside ``body``?

    Every point of ``solid``'s tessellation is classified against ``body`` with
    OCCT's exact solid classifier — no boolean, so none of a boolean's
    conditioning problems apply. The sweep over tolerances reports the smallest at
    which nothing is outside, which is an upper bound on how far the worst point
    protrudes.

    A point the sweep cannot clear is then **contradicted rather than believed**:
    its exact distance to ``body`` is measured with ``BRepExtrema_DistShapeShape``
    — no tolerance, no boolean — and reported as ``worst_distance_mm``. That is
    the quantity that separates a tie from a protrusion, and the classifier's own
    answer cannot, being a tolerance question asked about a point lying on a
    surface.

    **A candidate the distance probe did not answer for is not contained.**
    ``contained`` is only ever True on that arm when every point the sweep left
    over was actually measured -- ``probed == unswept`` -- and every measurement
    came back a tie. Two ways it can be less: a ``BRepExtrema_DistShapeShape``
    that does not report ``IsDone``, and a candidate set longer than
    :data:`DISTANCE_PROBE_MAX`, of which only the first slice is probed. Both
    used to leave ``worst_distance_mm`` at its initial 0.0 and read as a pass,
    which is the "unmeasured must never read as measured" rule ``material_outside``
    beside this was written to obey (docs/testing.md), broken on the check that
    stands in for it. ``unprobed`` reports the shortfall so a caller can say
    which of the two it is.

    **This is weaker than the exact cut and must be reported as such.** It samples
    the boundary rather than integrating the volume, and it is sound *here*
    because a lattice's material is everywhere within `t/2` of its own surface, so
    material outside the body implies surface points outside it. That argument is
    specific to this geometry, not general.
    """
    from OCP.BRep import BRep_Tool
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.BRepExtrema import BRepExtrema_DistShapeShape
    from OCP.BRepMesh import BRepMesh_IncrementalMesh
    from OCP.TopAbs import TopAbs_ShapeEnum, TopAbs_State
    from OCP.TopLoc import TopLoc_Location
    from OCP.TopoDS import TopoDS
    from OCP.gp import gp_Pnt

    BRepMesh_IncrementalMesh(solid, deflection, False, 0.2, True)
    pts = []
    for f in occ._explore(solid, TopAbs_ShapeEnum.TopAbs_FACE):
        loc = TopLoc_Location()
        tri = BRep_Tool.Triangulation_s(TopoDS.Face_s(f), loc)
        if tri is None:
            continue
        trsf = loc.Transformation()
        for i in range(1, tri.NbNodes() + 1):
            p = tri.Node(i).Transformed(trsf)
            pts.append((p.X(), p.Y(), p.Z()))
    if not pts:
        return {"sampled": 0, "outside": None, "cleared_at_mm": None,
                "unswept": 0, "worst_distance_mm": 0.0, "worst_distance_at": None,
                "probed": 0, "contained": False}
    pts = np.unique(np.array(pts), axis=0)

    # Each tolerance re-tests only what the previous one called OUT, which is
    # exact rather than an approximation: `Perform(p, tol)` calls a point ON
    # rather than OUT once it is within `tol` of the boundary, so raising the
    # tolerance can only ever turn OUT into ON, never the reverse. The OUT set
    # therefore shrinks monotonically along CONTAINMENT_SWEEP, and a point the
    # tightest tolerance already accepts cannot be reported by a looser one.
    #
    # Testing every point at every tolerance instead costs up to 5x for
    # nothing, and that is not academic at lattice scale: this is the check
    # that runs precisely where the exact cut was skipped, i.e. on the largest
    # solids, and one pass over them is already minutes.
    classifier = BRepClass3d_SolidClassifier(body)
    counts: dict[float, int] = {}
    cleared = None
    candidates = pts
    for tol in CONTAINMENT_SWEEP:
        still = []
        for p in candidates:
            classifier.Perform(gp_Pnt(*p), tol)
            if classifier.State() == TopAbs_State.TopAbs_OUT:
                still.append(p)
        counts[tol] = len(still)
        candidates = still
        if not still:
            cleared = tol
            break

    # **A point the sweep cannot clear is contradicted, not believed.** The
    # sweep tops out at CONTAINMENT_SWEEP[-1], so past it a 1e-09 mm tie and a
    # 1 mm protrusion report identically -- and the classifier does return OUT
    # for points that are exactly on the boundary, which is the failure mode
    # docs/testing.md already records ("nothing was actually wrong with the
    # geometry; the check simply could not say so").
    #
    # `BRepExtrema_DistShapeShape` answers the question the classifier is being
    # asked to approximate, without a tolerance and without a boolean: how far
    # is this point from the body? Material outside has somewhere to be, so a
    # real protrusion measures a real distance. Measured on `SpiralTest.step`
    # at cc=5, t=1, the 15 points its sweep could not clear are 0.000000 mm from
    # the input, median 9.95e-10 -- against the 0.239 to 0.771 mm the genuine
    # protrusions of docs/algorithm.md S7.2 measured. Eight orders apart.
    # **A probe that did not answer is not a probe that said "fine".** Two ways
    # this arm can measure less than it was asked to, and both have to reach the
    # verdict rather than the initial value of `worst_mm`: a
    # `BRepExtrema_DistShapeShape` that does not report `IsDone`, and a candidate
    # set longer than DISTANCE_PROBE_MAX, of which only the first slice is
    # measured at all. Counting successful probes is what separates "all of them
    # are ties" from "none of them was measured" -- the same distinction
    # `material_outside` draws when it reports a solid unmeasured rather than
    # passing it (docs/testing.md), on the check that stands in for it.
    worst_mm = 0.0
    worst_at = None
    probed = 0
    for p in candidates[:DISTANCE_PROBE_MAX]:
        vertex = BRepBuilderAPI_MakeVertex(gp_Pnt(*p)).Vertex()
        probe = BRepExtrema_DistShapeShape(vertex, body)
        probe.Perform()
        if not probe.IsDone():
            continue
        probed += 1
        if probe.Value() > worst_mm:
            worst_mm = probe.Value()
            worst_at = (round(float(p[0]), 4), round(float(p[1]), 4),
                        round(float(p[2]), 4))

    return {
        "sampled": int(len(pts)),
        "outside": counts[CONTAINMENT_SWEEP[0]],
        "counts": counts,
        "cleared_at_mm": cleared,
        "unswept": len(candidates),
        "worst_distance_mm": worst_mm,
        "worst_distance_at": worst_at,
        "probed": probed,
        "unprobed": len(candidates) - probed,
        "contained": (
            (cleared is not None and cleared <= CONTAINMENT_TOL)
            # Every candidate measured, and every measurement a tie.
            or (probed == len(candidates) and worst_mm <= PROTRUSION_MIN_MM)
        ),
    }


CUT_MAX_FACES = 25_000
"""Above this many faces, a solid's containment is not put to a boolean at all.

The per-solid cut below is the exact test and is preferred wherever it works,
but "wherever it works" has a measured ceiling and this is it. Three data
points, all on real output:

* `dense-lattice`, **15,966** faces — measures exactly 0 mm^3, in seconds.
* `TD_HX_rehearsal_test` at `cc=12, t=2.5`, **43,530** faces — ran to
  completion and reported the *entire* solid outside (354,733 mm^3), which is
  what the contradiction test below exists to catch.
* `SpiralTest` at `cc=5, t=1`, **45,897** faces — did not finish: 17 minutes of
  CPU and a 14.5 GB working set before it was stopped.

So past this scale the cut is not merely untrustworthy, it is not affordable
either, and attempting it buys nothing: the answer would be discarded by the
contradiction test and replaced by :func:`surface_points_outside` regardless.
Skipping straight to that check reaches the same verdict in seconds. The bar
sits between the scale that works and the scale that does not, and every
committed scenario stays below it; the only output that exceeds it is
`SpiralTest.step`'s, which is a fixture rather than a scenario.

**This does not weaken the guarantee, it relabels it.** A skipped solid is
reported as *unmeasured* with its containment evidence attached, exactly like
one whose cut was contradicted — never as a pass. What would weaken the
guarantee is a number nobody can afford to check.
"""


def material_outside(candidate_path: str, input_path: str) -> dict:
    """Generated material lying outside the input body (specification.md §1, §6.2).

    §1 requires the lattice to fit *exactly* within the input geometry, so this
    should be zero. Returns a dict rather than a float because on a large part the
    exact answer is not always available, and "not measured" must never be able to
    read as "measured as fine":

    * ``volume_mm3`` — the exact figure, summed over the output's solids, or
      ``None`` if any solid's cut could not be trusted;
    * ``unmeasured`` — one entry per such solid, with the containment evidence
      gathered for it instead;
    * ``exact`` — whether every solid was measured by boolean.

    **The cut is taken per solid, not over the whole output**, and that is the
    whole point of the shape of this function. Cutting a 43,530-face lattice
    against the body it was trimmed from returns a wrong answer: on the
    `cc=12, t=2.5` output of `TD_HX_rehearsal_test` it reported **354,733 mm³**
    outside — the entire solid — where the same cut, solid by solid, returns
    **exactly 0** for eight of the nine and only the largest misbehaves.

    That result is not a kernel no-op, which is worth stating because it is the
    obvious reading and it is wrong: the call reports ``IsDone``, ``HasModified``
    and ``HasGenerated``, and returns 43,672 faces where 43,530 went in, having
    removed 3.6e-03 mm³. It genuinely ran and genuinely mis-classified, which is
    unsurprising on this input — a lattice trimmed from a body has a large share
    of its faces lying *exactly on* that body's surface, the classic
    ill-conditioned case for a boolean. So a detector keyed on "the volume came
    back unchanged" would not fire: the relative difference is 1.02e-08.

    What does fire is a contradiction test. A trustworthy zero is exactly zero, so
    any non-zero remainder is checked against :func:`surface_points_outside`,
    which uses no boolean at all. If that finds the solid contained, the two
    disagree, the cut is not believed, and the solid is reported as unmeasured
    with the containment figures attached — labelled weaker, the way
    :func:`golden_sample_agreement` already is. A real defect survives this,
    because real material outside the body shows up in both tests.

    **A solid past :data:`CUT_MAX_FACES` skips the cut entirely** and goes
    straight to the containment evidence, reported as unmeasured with
    ``reason`` of ``"not attempted"`` rather than ``"contradicted"``. See that
    constant for the measurements behind the bar.
    """
    body = occ.read_step(input_path)
    solids = occ.solids(occ.read_step(candidate_path))
    total = 0.0
    unmeasured = []
    per_solid = []
    for i, solid in enumerate(solids):
        own = occ.volume(solid)
        faces = occ.count_subshapes(solid)[0]
        if faces > CUT_MAX_FACES:
            cut = float("nan")
            trustworthy = False
            reason = "not attempted"
        else:
            cut = _cut_volume(solid, body)
            # Exactly zero, in effect -- either sign: a cut can come back
            # negative (-1.5e-03 mm^3 on a 151 mm^3 support body of the
            # rehearsal part, 0 of 273 surface points outside), and a negative
            # remainder is no more a measurement than a positive one.
            trustworthy = cut == cut and abs(cut) <= own * 1e-9
            reason = "contradicted"
        if not trustworthy:
            evidence = surface_points_outside(solid, body)
            unmeasured.append({
                "solid": i,
                "faces": faces,
                "volume_mm3": own,
                "cut_claimed_mm3": cut,
                "reason": reason,
                "containment": evidence,
            })
        else:
            total += cut
        per_solid.append({"solid": i, "cut_mm3": cut, "trusted": trustworthy})
    return {
        "volume_mm3": total if not unmeasured else None,
        "per_solid": per_solid,
        "unmeasured": unmeasured,
        "exact": not unmeasured,
        "solids": len(solids),
    }


def bounding_box_within(candidate_path: str, input_path: str, tol: float) -> bool:
    lo_c, hi_c = occ.bounding_box(occ.read_step(candidate_path))
    lo_i, hi_i = occ.bounding_box(occ.read_step(input_path))
    return bool(np.all(lo_c >= lo_i - tol) and np.all(hi_c <= hi_i + tol))


# --- print supports (specification.md §4.6, §6.2) ----------------------------

SUPPORT_MESH_DEFLECTION = 0.02
"""Millimetres. Chordal deviation of the meshes the two support checks use.

Much finer than the classification mesh, because both checks ask about material
within a few hundredths of a millimetre of a curved surface."""

SUPPORT_INWARD_STEP = 0.06
"""Millimetres an overhang sample is moved into the body before it is asked
whether the output has material there: three times the mesh deflection, so a
point that is really inside is not put outside by the mesh's own error."""


SUPPORT_EXACT_STEP = 2e-3
"""Millimetres. Probe depth for the exact re-check of a mesh miss — the depth
the pipeline's own coverage proof uses."""

COVERAGE_SLIVER_TOL = 0.05
"""Millimetres. A bare sample this near an output solid is reported as a sliver
and does not fail coverage. Where a ceiling only grazes a support node, trimming
the node can leave a chip of about 1e-3 mm^3 on the far side of a grid plane,
which the generator discards rather than ship a floating micro-solid: well under
a powder grain, and no support could print it."""

EXACT_RECHECK_MAX = 2000
"""Most mesh misses the exact classifier is asked about. Past this the misses
are reported as they stand: that many is a real gap, not a thin wedge."""


def _oriented_mesh(shape, deflection: float):
    """``(verts, tris, normals)`` with every triangle's normal pointing out of
    the solid — `face_triangulation` reports a face's triangles in the surface's
    own sense, which a reversed face inverts."""
    from OCP.TopAbs import TopAbs_Orientation

    occ.mesh_shape(shape, deflection)
    verts, tris, offset = [], [], 0
    for face in occ.faces(shape):
        got = occ.face_triangulation(face)
        if got is None:
            continue
        v, f = got
        if face.Orientation() == TopAbs_Orientation.TopAbs_REVERSED:
            f = f[:, ::-1]
        verts.append(v)
        tris.append(f + offset)
        offset += len(v)
    v = np.vstack(verts)
    f = np.vstack(tris)
    n = np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]])
    area = 0.5 * np.linalg.norm(n, axis=1)
    keep = area > 0.0
    return v, f[keep], n[keep] / (2.0 * area[keep])[:, None], area[keep]


EXACT_NORMAL_PREFILTER_DEG = 5.0
"""Degrees. How much flatter than the limit a facet may read and still be
re-read at the exact normal: far more than a chord at
:data:`SUPPORT_MESH_DEFLECTION` tilts on any face a support carries."""


def _oriented_mesh_exact(shape, deflection: float):
    """:func:`_oriented_mesh` plus ``exact(i)``: the true outward normal of
    triangle ``i``'s face at the triangle's centroid in the face's parameters,
    ``None`` at a pole or where the face has no parameters. ``exact(i,
    point=True)`` is the point there instead, on the face itself."""
    from OCP.TopAbs import TopAbs_Orientation

    from latticegen2.support.overhang import _SurfaceNormals

    occ.mesh_shape(shape, deflection)
    verts, tris, owner, uvc, faces, offset = [], [], [], [], [], 0
    for face in occ.faces(shape):
        got = occ.face_uv_triangulation(face)
        if got is None:
            continue
        p, uv, f, _evaluate = got
        if face.Orientation() == TopAbs_Orientation.TopAbs_REVERSED:
            f = f[:, ::-1]
        verts.append(p)
        tris.append(f + offset)
        uvc.append(uv[f].mean(axis=1))
        owner.append(np.full(len(f), len(faces), dtype=np.int64))
        faces.append(face)
        offset += len(p)
    v = np.vstack(verts)
    f = np.vstack(tris)
    owner = np.concatenate(owner)
    uvc = np.vstack(uvc)
    nrm = np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]])
    area = 0.5 * np.linalg.norm(nrm, axis=1)
    keep = area > 0.0
    f, owner, uvc = f[keep], owner[keep], uvc[keep]
    normals = {}

    def exact(i: int, point: bool = False):
        k = int(owner[i])
        if k not in normals:
            normals[k] = _SurfaceNormals(faces[k])
        return normals[k].at(*uvc[i])[0 if point else 1]

    return v, f, nrm[keep] / (2.0 * area[keep])[:, None], area[keep], exact


def overhang_coverage(candidate_path: str, input_path: str, cc: float, t: float,
                      orient, thickness: float, pitch: float | None = None) -> dict:
    """Is every inner overhang of the input backed by material in the output?

    Asked of the **written file**, and independently of how the supports were
    placed or built: the input's overhangs are sampled afresh, each sample is
    moved a hair into the body, and the output's own tessellation is asked
    whether that point is inside it. Nothing here reads the plan. ``pitch`` is
    the sample spacing, the support step by default: a finer one finds ceiling
    bare over less than that, which clipped corners of support nodes are.
    """
    from latticegen2.classify import PointInside, tessellate_surface
    from latticegen2.support import overhang as oh

    lp = lattice_params(cc, t, orient)
    body = occ.read_step(input_path)
    tessellate_surface(body, lp)
    samples = oh.detect(body, lp.b) if pitch is None else oh.detect(body, lp.b, pitch)
    out = occ.read_step(candidate_path)
    v, f, _n, _a = _oriented_mesh(out, SUPPORT_MESH_DEFLECTION)
    v, f = _weld(v, f)
    inside = PointInside(TriMesh(verts=v, tris=f))
    probe = samples.points + SUPPORT_INWARD_STEP * samples.inward
    covered = inside(probe)
    # **A miss on the mesh is contradicted against the exact solid before it is
    # reported.** The mesh needs a probe well clear of its own deflection, and
    # near the border of an overhang region the support under the wall is a
    # wedge thinner than that: the wall's tilt is within a few degrees of the
    # support's own face there, so the two diverge slowly. Measured on
    # `test-cylinder.STEP`: 9 of 51,156 probes fall out of the bottom of a
    # wedge 0.01-0.06 mm thick, every one of them on a sample the solid touches
    # (distance exactly 0) and contains 2e-3 mm below. The exact classifier is
    # asked at that depth; it costs tens of milliseconds a point, which is why
    # it is the second question and not the first. (Measured on v4.0.0's first
    # support algorithm, since removed; the mechanism is a property of any
    # support meeting a wall near the limit, not of that algorithm.)
    doubtful = np.nonzero(~covered)[0]
    thin = 0
    if 0 < len(doubtful) <= EXACT_RECHECK_MAX:
        from OCP.BRepClass3d import BRepClass3d_SolidClassifier
        from OCP.TopAbs import TopAbs_State
        from OCP.gp import gp_Pnt

        classifiers = [BRepClass3d_SolidClassifier(s) for s in occ.solids(out)]
        for i in doubtful:
            q = samples.points[i] + SUPPORT_EXACT_STEP * samples.inward[i]
            for classifier in classifiers:
                classifier.Perform(gp_Pnt(*q), 1e-7)
                if classifier.State() in (TopAbs_State.TopAbs_IN, TopAbs_State.TopAbs_ON):
                    covered[i] = True
                    thin += 1
                    break
    # A sample within a hair of a body edge can probe out of the body: its rim
    # is continued by a wall, and nothing inside the body can lie behind it.
    # Bare samples are then measured against the output, which is exact.
    rim = slivers = 0
    if 0 < int((~covered).sum()) <= EXACT_RECHECK_MAX:
        from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
        from OCP.BRepClass3d import BRepClass3d_SolidClassifier
        from OCP.BRepExtrema import BRepExtrema_DistShapeShape
        from OCP.TopAbs import TopAbs_State
        from OCP.gp import gp_Pnt

        in_body = BRepClass3d_SolidClassifier(body)
        solids = occ.solids(out)
        boxes = [occ.bounding_box(s) for s in solids]
        for i in np.nonzero(~covered)[0]:
            q = samples.points[i] + SUPPORT_EXACT_STEP * samples.inward[i]
            in_body.Perform(gp_Pnt(*map(float, q)), 1e-7)
            if in_body.State() == TopAbs_State.TopAbs_OUT:
                covered[i] = True
                rim += 1
                continue
            p = samples.points[i]
            vertex = BRepBuilderAPI_MakeVertex(gp_Pnt(*map(float, p))).Vertex()
            for s, (lo, hi) in zip(solids, boxes):
                if (np.any(p < lo - COVERAGE_SLIVER_TOL) or np.any(p > hi + COVERAGE_SLIVER_TOL)):
                    continue
                d = BRepExtrema_DistShapeShape(vertex, s)
                if d.IsDone() and d.Value() <= COVERAGE_SLIVER_TOL:
                    covered[i] = True
                    slivers += 1
                    break
    missing = samples.points[~covered]
    return {
        "samples": int(len(probe)),
        "uncovered": int(len(missing)),
        "thin": thin,
        "rim": rim,
        "slivers": slivers,
        "area_mm2": float(samples.area),
        "where": [np.round(p, 3).tolist() for p in missing[:5]],
    }


def self_supporting(candidate_path: str, input_path: str, cc: float, t: float,
                    orient, tol_deg: float = 0.05) -> dict:
    """Does any free surface of the output face down more flatly than the limit?

    A facet is exempt only for one of two reasons, both geometric: it lies on
    the input's own surface, where the output is closed against the wall rather
    than hanging free; or it is a lattice strut face, which is recognised by its
    normal and which the limit does not apply to (specification.md §4.6 — the
    lattice is not analysed). Everything else is support, and a single such
    facet over the limit is a failure.
    """
    from latticegen2.cli import SUPPORT_OVERHANG_DEG

    lp = lattice_params(cc, t, orient)
    b = lp.b
    out = occ.read_step(candidate_path)
    v, f, n, area, exact = _oriented_mesh_exact(out, SUPPORT_MESH_DEFLECTION)

    # **Read at exact normals** (specification.md §6.2). A support face built at
    # the limit is curved, and its chords tilt either side of the true surface,
    # so a facet can read a hair past the bar where the face is not. The facet
    # normals only pick candidates, against a bar loose enough to hold every
    # chord the mesh can make; each candidate is then read on the face itself.
    bar = float(np.sin(np.radians(SUPPORT_OVERHANG_DEG + tol_deg)))
    loose = float(np.sin(np.radians(SUPPORT_OVERHANG_DEG - EXACT_NORMAL_PREFILTER_DEG)))
    suspect = np.nonzero(-(n @ b) > loose)[0]
    n = n.copy()
    for i in suspect:
        got = exact(int(i))
        if got is not None:
            n[i] = got
    down = -(n @ b)
    suspect = suspect[down[suspect] > bar]

    strut = np.array([s * (su * lp.u[k] + sv * lp.v[k]) / np.sqrt(2.0)
                      for k in range(3) for su in (1, -1) for sv in (1, -1)
                      for s in (1,)])
    is_strut = (np.abs(n[suspect] @ strut.T) > 1.0 - 1e-9).any(axis=1)
    suspect = suspect[~is_strut]

    body = occ.read_step(input_path)
    occ.mesh_shape(body, SUPPORT_MESH_DEFLECTION)
    bv, bf = [], []
    offset = 0
    for face in occ.faces(body):
        got = occ.face_triangulation(face)
        if got is None:
            continue
        bv.append(got[0])
        bf.append(got[1] + offset)
        offset += len(got[0])
    bv, bf = _weld(np.vstack(bv), np.vstack(bf))
    wall = TriMesh(verts=bv, tris=bf)
    sh = SpatialHash(wall, min_cell=1.0)
    bad_area = 0.0
    where = []
    worst = 0.0
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
    from OCP.BRepExtrema import BRepExtrema_DistShapeShape
    from OCP.gp import gp_Pnt

    def on_body(i) -> bool:
        # The facet's own point on its face, measured against the input exactly:
        # a facet centre sits off a curved face by the chord's sag, and on
        # `TD_HX_rehearsal_test` OCCT's real sag is a multiple of the deflection
        # asked for (docs/algorithm.md §5.1) -- a support face lying on a
        # countersink's 60 deg cone read 0.097 mm off the input's mesh.
        q = exact(int(i), point=True)
        d = BRepExtrema_DistShapeShape(BRepBuilderAPI_MakeVertex(gp_Pnt(*map(float, q))).Vertex(),
                                       body)
        return bool(d.IsDone() and d.Value() <= CONTAINMENT_TOL)

    for i in suspect:
        centre = v[f[i]].mean(axis=0)
        if _distance_to_mesh(wall, sh, centre) <= 3.0 * SUPPORT_MESH_DEFLECTION:
            continue                      # closed against the wall
        if on_body(i):
            continue
        bad_area += float(area[i])
        worst = max(worst, float(np.degrees(np.arcsin(min(1.0, down[i])))))
        if len(where) < 5:
            where.append(np.round(centre, 3).tolist())
    return {
        "triangles": int(len(f)),
        "bad_area_mm2": bad_area,
        "worst_deg": worst,
        "where": where,
    }


SUPPORT_OVERLAP_TOL = 1e-6
"""mm³. Volume a support may share with the lattice or another support: the
stage's own bar (docs/algorithm.md §14.2)."""

CONTACT_PLANE_GAP = 1e-2
"""mm. A planar support face parallel to a strut face and nearer to its plane
than this is taken to rest on it, and must then lie on it."""

CONTACT_PLANE_TOL = 1e-6
"""mm. How far off a strut's plane a support face resting on it may lie."""


def _split_supports(candidate_path: str, n_supports: int):
    """``(lattice, supports)`` solids of a written file: the generator writes the
    lattice first and the supports after it (docs/algorithm.md §14.2)."""
    solids = occ.solids(occ.read_step(candidate_path))
    if not 0 < n_supports <= len(solids):
        raise RuntimeError(f"{candidate_path} holds {len(solids)} solid(s), "
                           f"not the {n_supports} support(s) its log reports")
    return solids[:-n_supports], solids[-n_supports:]


def _common(a, b) -> float:
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Common

    alg = BRepAlgoAPI_Common(a, b)
    if not alg.IsDone():
        return float("nan")
    return occ.volume(alg.Shape())


def supports_disjoint(candidate_path: str, n_supports: int) -> dict:
    """Do the supports overlap the lattice or each other? Exact, one boolean
    per pair whose boxes touch."""
    lattice, supports = _split_supports(candidate_path, n_supports)
    gap = 1e-3

    def box(s):
        lo, hi = occ.bounding_box(s)
        return np.asarray(lo) - gap, np.asarray(hi) + gap

    def touch(p, q):
        return bool(np.all(p[0] <= q[1]) and np.all(q[0] <= p[1]))

    lb = [box(s) for s in lattice]
    sb = [box(s) for s in supports]
    worst, pairs, where = 0.0, 0, None
    for i, s in enumerate(supports):
        others = [(lattice[j], lb[j]) for j in range(len(lattice))]
        others += [(supports[j], sb[j]) for j in range(i + 1, len(supports))]
        for o, ob in others:
            if not touch(sb[i], ob):
                continue
            pairs += 1
            v = _common(s, o)
            # NaN counts as the worst, and stays the worst: a failed boolean is
            # an unverified pair, which no later reading can vouch for.
            if worst == worst and not v <= worst:
                worst = v
                where = np.round((sb[i][0] + sb[i][1]) / 2.0, 3).tolist()
    return {"pairs": pairs, "worst_mm3": worst, "where": where}


def supports_on_strut_planes(candidate_path: str, n_supports: int) -> dict:
    """Does every support face resting on a strut lie on that strut's plane?

    Planar faces only, since a strut face is planar: a support face parallel to
    a lattice face, overlapping it and nearer than :data:`CONTACT_PLANE_GAP`,
    rests on it, and must lie within :data:`CONTACT_PLANE_TOL` of its plane, so
    the downstream merge finds one shared face rather than a sliver.
    """
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.GeomAbs import GeomAbs_SurfaceType
    from OCP.TopAbs import TopAbs_Orientation

    def planes(solids):
        out = []
        for s in solids:
            for face in occ.faces(s):
                surf = BRepAdaptor_Surface(face)
                if surf.GetType() != GeomAbs_SurfaceType.GeomAbs_Plane:
                    continue
                pl = surf.Plane()
                d, o = pl.Axis().Direction(), pl.Location()
                n = np.array([d.X(), d.Y(), d.Z()])
                if face.Orientation() == TopAbs_Orientation.TopAbs_REVERSED:
                    n = -n
                lo, hi = occ.bounding_box(face)
                out.append((n, float(n @ np.array([o.X(), o.Y(), o.Z()])),
                            np.asarray(lo), np.asarray(hi)))
        return out

    lattice, supports = _split_supports(candidate_path, n_supports)
    lat = planes(lattice)
    if not lat:
        return {"resting": 0, "worst_mm": 0.0, "where": None}
    ln = np.array([p[0] for p in lat])
    lo_ = np.array([p[1] for p in lat])
    llo = np.array([p[2] for p in lat])
    lhi = np.array([p[3] for p in lat])
    resting, worst, where = 0, 0.0, None
    for n, off, lo, hi in planes(supports):
        # Resting means facing the strut: outward normals opposite.
        near = ((ln @ n) < -1.0 + 1e-9) \
            & np.all(llo <= hi + CONTACT_PLANE_GAP, axis=1) \
            & np.all(lo <= lhi + CONTACT_PLANE_GAP, axis=1)
        if not near.any():
            continue
        gaps = np.abs(off + lo_[near])             # opposite normals: offsets cancel
        hit = gaps < CONTACT_PLANE_GAP
        if not hit.any():
            continue
        resting += 1
        g = float(gaps[hit].max())
        if g > worst:
            worst = g
            where = np.round((lo + hi) / 2.0, 3).tolist()
    return {"resting": resting, "worst_mm": worst, "where": where}


def contained_in(inner_path: str, outer_path: str) -> float:
    """Volume of ``inner`` lying outside ``outer``, in mm³ — exact, one way."""
    inner = occ.read_step(inner_path)
    outer = occ.read_step(outer_path)
    return _cut_volume(inner, outer)
