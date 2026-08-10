"""Geometry-validity checks for the E2E harness (specification.md §6.2).

Ported from the Julia implementation's ``tools/verify_geometry.jl``. These are QA
tooling, not part of the generation pipeline, and they reason only about the
finished STEP file: they re-read it, re-tessellate it, and check the resulting
triangles. That keeps them independent of *how* the geometry was built, which is
the property that makes them worth running at all.

The Julia originals ran on a different language runtime as well. That
independence went away with the Julia toolchain (keeping it would have meant
keeping Julia and gmsh for one script), and is more than repaid by
:func:`brepcheck` — OCCT's own exact B-rep validity test, which the gmsh-based
implementation could not reach at all and which left docs/algorithm.md §11.1's
self-intersection question answerable only by indirect evidence.
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
    because their independent triangulations are not vertex-aligned — the
    false-positive mode diagnosed in docs/algorithm.md §11.1.
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


def golden_sample_agreement(candidate_path: str, golden_path: str, samples: int = 200_000,
                            seed: int = 20260810) -> dict:
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

    cand_mesh = mesh_of(candidate_path, 10.0, 1.5)
    gold_mesh = mesh_of(golden_path, 10.0, 1.5)
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
    for i in disputed:
        p = pts[i]
        if min(_distance_to_mesh(cand_mesh, p), _distance_to_mesh(gold_mesh, p)) > SURFACE_TIE_TOL:
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


def _distance_to_mesh(mesh: TriMesh, point: np.ndarray, search: float = 1.0) -> float:
    """Distance from a point to the nearest triangle within ``search`` mm."""
    from latticegen2.classify import segment_triangle_dist

    sh = SpatialHash(mesh, min_cell=search)
    candidates = sh.query(point - search, point + search)
    if len(candidates) == 0:
        return float("inf")
    A, B, C = mesh.triangle_points
    return float(segment_triangle_dist(point, point, A[candidates], B[candidates],
                                       C[candidates]).min())


def material_outside(candidate_path: str, input_path: str) -> float:
    """Volume of generated material lying outside the input body, in mm³.

    specification.md §1 requires the lattice to fit *exactly* within the input
    geometry, so this should be zero. It is a direct check of that requirement,
    independent of any golden sample.
    """
    return _cut_volume(occ.read_step(candidate_path), occ.read_step(input_path))


def bounding_box_within(candidate_path: str, input_path: str, tol: float) -> bool:
    lo_c, hi_c = occ.bounding_box(occ.read_step(candidate_path))
    lo_i, hi_i = occ.bounding_box(occ.read_step(input_path))
    return bool(np.all(lo_c >= lo_i - tol) and np.all(hi_c <= hi_i + tol))
