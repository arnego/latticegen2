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
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Fuse
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
from ..errors import ProcessingError
from ..junction import build_template, is_cap_plane_face
from ..lattice import LatticeParams, lattice_params, nodes
from ..parallel import WorkerPool
from ..parallel import compound_children as _compound_children
from ..parallel import read_brep as _read_brep
from . import polyhedra as ph

NodeKey = tuple[int, int, int]

MIN_PIECE_VOLUME = 1e-9
"""mm³. A clipped piece smaller than this is not built. It is far below any
feature the sew can resolve, so both cells sharing it drop it alike."""

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


def build_cell(
    lp: LatticeParams,
    tpl,
    node_pos: np.ndarray,
    clipped: list,
    body: TopoDS_Shape | None,
    with_junction: bool,
    trim_first: bool = False,
    clock: dict | None = None,
) -> list:
    """The solids of one lattice cell: junction plus support, inside the body.

    ``clipped`` is ``(verts, faces)`` per primitive already cut to the cell.
    ``body`` is ``None`` when nothing of its surface is near the cell, and the
    intersection is then skipped. ``trim_first`` swaps the operand order — each
    operand is intersected with the body and the results fused — and is the
    retry used when the usual order leaves a sample uncovered.

    ``clock``, if given, accumulates seconds under ``"solids"``, ``"fuse"`` and
    ``"trim"``.
    """
    clock = clock if clock is not None else {}
    t0 = time.perf_counter()
    operands = []
    if with_junction:
        operands.append(tpl.solid.Moved(occ.translation(node_pos)))
    operands.extend(polyhedron_solid(v, f) for v, f in clipped)
    clock["solids"] = clock.get("solids", 0.0) + time.perf_counter() - t0
    if not operands:
        return []
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

    def common(shape):
        t1 = time.perf_counter()
        algo = _boolean(BRepAlgoAPI_Common(), [shape], [body],
                        f"Intersecting the cell at {where} with the input body")
        clock["trim"] = clock.get("trim", 0.0) + time.perf_counter() - t1
        return occ.solids(algo.Shape())

    if body is None:
        return fuse(operands)
    if trim_first:
        trimmed = [s for op in operands for s in common(op)]
        return fuse(trimmed) if trimmed else []
    return [s for solid in fuse(operands) for s in common(solid)]


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
            got = ph.clip_to_cell(prims[pid][0], prims[pid][1], lp, pos,
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
        if not ok:
            n_retried += 1
            solids = build_cell(lp, tpl, pos, clipped, local, info["junction"],
                                trim_first=True, clock=clock)
            miss = ~_covered(solids, points)
            if miss.any():
                p = points[np.nonzero(miss)[0][0]]
                raise ProcessingError(
                    f"Support generation left {int(miss.sum())} overhang "
                    f"sample(s) uncovered in the lattice cell of node {node}, "
                    f"the first at [{p[0]:.3f}, {p[1]:.3f}, {p[2]:.3f}], with "
                    f"both operand orders. The output would have an unsupported "
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
        if probe is not None and pieces:
            worst_outside = max(worst_outside,
                                probe.worst_outside(_piece_vertices(pieces)))
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
    moved = samples.points + INWARD_STEP * samples.inward
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
        p = samples.points[by_cell[n][0]]
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
            phases[key] = phases.get(key, 0.0) + seconds
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
    if n_points != len(samples):
        raise ProcessingError(
            f"Only {n_points} of {len(samples)} overhang samples were checked "
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


def fill_voids(shells: dict) -> tuple[dict, int, float]:
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
    group of positive volume would be a separate body the junction graph says
    cannot exist, and is a hard failure rather than something to discard.
    """
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Shell

    from ..connect import UnionFind

    builder = BRep_Builder()
    out: dict = {}
    n_filled = 0
    filled_volume = 0.0
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
        for vol, _sub in measured[1:]:
            if vol > 0.0:
                raise ProcessingError(
                    f"Component {group} assembled into more than one body "
                    f"({vol:.4f} mm^3 beside the main one) where the junction "
                    f"graph proves it is one."
                )
            n_filled += 1
            filled_volume += -vol
        out[group] = measured[0][1]
    return out, n_filled, filled_volume
