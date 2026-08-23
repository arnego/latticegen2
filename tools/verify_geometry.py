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


def mesh_of(path: str, cc: float, t: float) -> TriMesh:
    """Re-read a STEP file and tessellate it for the mesh-based checks."""
    shape = occ.read_step(path)
    lp = lattice_params(cc, t)
    occ.mesh_shape(shape, chordal_target(lp))
    verts, tris, offset = [], [], 0
    for face in occ.faces(shape):
        got = occ.face_triangulation(face)
        if got is None:
            continue
        v, f = got
        verts.append(v)
        tris.append(f + offset)
        offset += len(v)
    if not tris:
        raise RuntimeError(f"No triangles produced for {path}")
    v, f = _weld(np.vstack(verts), np.vstack(tris))
    return TriMesh(verts=v, tris=f)


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
    """
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

    ``fraction`` is the quantity to judge on: how much of a body's surface is
    described more loosely than its own feature size. ``over`` is OCCT's own
    per-edge predicate, reported rather than gated on - the 80 mm ball's
    accepted output has 73 of 5,760 pairs past their recorded tolerance, every
    one by ~4 % at 3.8e-07 mm, which is a knife edge and not a defect. Neither
    is ``worst``: measured on `SpiralTest`, the sound dominant body scores
    *worse* on it than the island that cannot be written.
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

    **This is weaker than the exact cut and must be reported as such.** It samples
    the boundary rather than integrating the volume, and it is sound *here*
    because a lattice's material is everywhere within `t/2` of its own surface, so
    material outside the body implies surface points outside it. That argument is
    specific to this geometry, not general.
    """
    from OCP.BRep import BRep_Tool
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
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
        return {"sampled": 0, "outside": None, "cleared_at_mm": None}
    pts = np.unique(np.array(pts), axis=0)

    classifier = BRepClass3d_SolidClassifier(body)
    counts: dict[float, int] = {}
    cleared = None
    for tol in CONTAINMENT_SWEEP:
        n = 0
        for p in pts:
            classifier.Perform(gp_Pnt(*p), tol)
            if classifier.State() == TopAbs_State.TopAbs_OUT:
                n += 1
        counts[tol] = n
        if n == 0:
            cleared = tol
            break
    return {
        "sampled": int(len(pts)),
        "outside": counts[CONTAINMENT_SWEEP[0]],
        "counts": counts,
        "cleared_at_mm": cleared,
        "contained": cleared is not None and cleared <= CONTAINMENT_TOL,
    }


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
    """
    body = occ.read_step(input_path)
    solids = occ.solids(occ.read_step(candidate_path))
    total = 0.0
    unmeasured = []
    per_solid = []
    for i, solid in enumerate(solids):
        own = occ.volume(solid)
        cut = _cut_volume(solid, body)
        trustworthy = cut == cut and cut <= own * 1e-9      # exactly zero, in effect
        if not trustworthy:
            evidence = surface_points_outside(solid, body)
            unmeasured.append({
                "solid": i,
                "faces": occ.count_subshapes(solid)[0],
                "volume_mm3": own,
                "cut_claimed_mm3": cut,
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
