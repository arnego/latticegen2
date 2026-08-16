"""Pipeline orchestration.

    parse -> template -> import -> tessellate -> classify
          -> trim boundary -> connect -> instance interior -> weld
          -> simplify -> validate -> write STEP -> summary

The shape of this is deliberately flat: no hierarchical assembly, no distributed
merge rounds, no fuse-based cleanup, because there is nothing left for them to
do — the interior is instanced rather than fused, connectivity is a graph
property rather than a boolean experiment, and assembly is an index lookup
rather than a geometric search. The one exception is inside "weld": the
boundary-layer sew tiles a large component's pieces before sewing them
(docs/algorithm.md §8), because sewing itself is the one remaining
operation that genuinely scales worse than linearly (docs/algorithm.md §8). That
stays an internal detail of one stage, not a pipeline-level tiling stage of its
own — the flow above is still accurate at the granularity it describes.
"""

from __future__ import annotations

import datetime as _dt
import os
import shutil
import sys
from dataclasses import dataclass

import numpy as np
from OCP.BRepTools import BRepTools
from OCP.Standard import Standard_Failure
from OCP.TopoDS import TopoDS_Shape

from . import occ, weld
from .boundary import (
    CAP_AREA_REL_TOL,
    finalize_pieces,
    fuse_disagreeing_pairs,
    resolve_interfaces,
    trim_boundary,
)
from .classify import Klass, classify_nodes, tessellate_surface
from .cli import Args
from .connect import build_components
from .errors import InputGeometryError, OutputError, ProcessingError
from .interior import build_interior_shell, extract_template_mesh
from .junction import build_template
from .lattice import candidate_nodes, lattice_params, neighbor_step, part_name
from .lattice import node as lattice_node
from .parallel import WorkerPool
from .parallel import read_brep as _read_brep
from .parallel import write_brep as _write_brep
from .runlog import RunLog, Timer, format_bytes
from .stepout import generation_params_text, rewrite_step_header

def _make_tmpdir(output_path: str) -> str:
    stamp = _dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    tmpdir = os.path.join(os.path.dirname(output_path) or ".", "temp", stamp)
    os.makedirs(tmpdir, exist_ok=True)
    return tmpdir


def run_pipeline(args: Args, rl: RunLog) -> dict:
    """Generate the lattice and write it, returning the stats for the summary.

    On any failure the temporary directory is left in place for post-mortem
    analysis and its path is reported on the console, not just in the log
    (specification.md §4.4, §7).
    """
    tmpdir = _make_tmpdir(args.output)
    try:
        return _run(args, rl, tmpdir)
    except BaseException:
        rl.always(f"Intermediate files kept for analysis in: {tmpdir}")
        raise


def _run(args: Args, rl: RunLog, tmpdir: str) -> dict:
    occ.quiet_kernel()
    lp = lattice_params(args.cc, args.t)
    stats: dict[str, object] = {}
    rl.line(f"temp directory: {tmpdir}")

    # The template also validates the parameter pair geometrically: it fails
    # fast (exit 2) if the mid-strut cap interfaces do not survive the fuse.
    with Timer(rl, "template"):
        tpl = build_template(lp)
        tmesh = extract_template_mesh(lp, tpl)
    rl.line(
        f"junction template: {tpl.n_faces} faces, {len(tmesh.verts)} vertices, "
        f"volume {tpl.volume:.6f} mm^3, worst cap area error {tpl.cap_area_error:.2e}"
    )

    with Timer(rl, "import"):
        body = occ.read_step(args.input)
        if not occ.solids(body):
            raise InputGeometryError(
                f"Input STEP contains no solid bodies: {args.input}. The lattice "
                f"bounds must be defined by a solid, not by surfaces alone."
            )
        lo, hi = occ.bounding_box(body)
    rl.line(f"input bounding box: lo={np.round(lo, 4).tolist()} hi={np.round(hi, 4).tolist()}")

    body_path = os.path.join(tmpdir, "input.brep")
    BRepTools.Write_s(body, body_path)

    with Timer(rl, "tessellate"):
        mesh = tessellate_surface(body, lp)
    rl.line(
        f"surface mesh: {len(mesh.tris)} triangles, {len(mesh.verts)} vertices, "
        f"measured chordal deviation d={mesh.deviation:.5f} mm "
        f"(classification margin r+d={lp.r + mesh.deviation:.5f} mm)"
    )

    with Timer(rl, "classify"):
        candidates = candidate_nodes(lp, lo, hi)
        classification = classify_nodes(lp, mesh, candidates)
    counts = classification.counts()
    rl.line(
        f"classification: {counts['candidates']} candidate nodes -> "
        f"interior {counts['interior']}, boundary {counts['boundary']}, "
        f"outside {counts['outside']}"
    )
    stats.update({f"nodes_{k}": v for k, v in counts.items()})

    interior_nodes = classification.of(Klass.INTERIOR)
    boundary_nodes = classification.of(Klass.BOUNDARY)
    if len(interior_nodes) == 0 and len(boundary_nodes) == 0:
        raise ProcessingError(
            "No lattice nodes fall inside the input geometry. The volume is likely "
            "smaller than one lattice cell at these parameters — try a smaller -cc."
        )

    # One worker pool for the whole rest of the run (docs/algorithm.md §8, §12),
    # shared by every parallel stage from here on — boundary trim, boundary-sew,
    # same-domain unification and validation — instead of each stage building and
    # tearing down its own. Built even when `args.workers <= 1`: `WorkerPool` is
    # then inert (`.active` is False) and every stage below takes its own
    # sequential path exactly as if no pool existed, so this costs nothing when
    # there is nothing to parallelise. Entered and exited explicitly rather than
    # as a `with` block wrapping the rest of this function, purely to avoid
    # re-indenting the whole boundary-through-validate span; the effect is
    # identical, and a Ctrl+C here is still caught at the top level as
    # `CancelledError` after the pool tears its workers down (docs/algorithm.md
    # §10) — `Pool.__exit__` terminates unconditionally, so it does not matter
    # that no real exception info is threaded through this `finally`.
    pool = WorkerPool(args.workers)
    pool.__enter__()
    try:
        return _run_with_pool(
            args, rl, tmpdir, lp, tpl, tmesh, body, body_path, stats,
            interior_nodes, boundary_nodes, pool,
        )
    finally:
        pool.__exit__(*sys.exc_info())


def _run_with_pool(
    args: Args,
    rl: RunLog,
    tmpdir: str,
    lp,
    tpl,
    tmesh,
    body,
    body_path: str,
    stats: dict,
    interior_nodes,
    boundary_nodes,
    pool: WorkerPool,
) -> dict:
    """The boundary-through-export span of :func:`_run`, run under one shared pool.

    Split out only so that span can be wrapped in a plain function call inside
    :func:`_run`'s ``try``/``finally`` rather than a deeply re-indented ``with``
    block — the pool's lifetime is what matters, not where this split falls.
    """
    with Timer(rl, "boundary"):
        # Reported per decile crossed rather than on an exact modulo: the parallel
        # path advances in whole batches, so a `done % step == 0` test lands on a
        # multiple only by luck and most of the run goes unreported.
        reported = [0]

        def progress(done: int, total: int) -> None:
            step = max(1, total // 10)
            if done >= total or done - reported[0] >= step:
                reported[0] = done
                rl.line(f"  boundary trim: {done}/{total}")

        boundary = trim_boundary(
            lp, tpl, boundary_nodes, body, body_path, tmpdir,
            workers=args.workers, progress=progress,
            pool=pool,
        )
    rl.line(
        f"boundary trim: {len(boundary.pieces)} pieces from {len(boundary_nodes)} "
        f"junctions ({boundary.n_empty} produced no geometry)"
    )
    stats["boundary_pieces"] = len(boundary.pieces)
    stats["boundary_empty"] = boundary.n_empty
    stats["workers"] = args.workers
    if boundary.max_worker_rss:
        # Folded into the run's high-water mark as well as reported separately:
        # the summary's "Peak memory" must be the maximum this run used anywhere,
        # not just in the master (specification.md §3).
        rl.note_worker_rss(boundary.max_worker_rss)
        rl.line(f"peak worker memory: {format_bytes(boundary.max_worker_rss)}")
        stats["peak_worker_memory"] = format_bytes(boundary.max_worker_rss)

    with Timer(rl, "connect"):
        interior_set = {(int(a), int(b), int(c)) for a, b, c in interior_nodes} \
            if len(interior_nodes) else set()
        iface = resolve_interfaces(lp, interior_nodes, boundary.pieces)
        n_fused = 0
        if iface.mismatched:
            # Declining a mismatched cap is not a safe degradation on its own:
            # where the two sides present *different* partial regions, keeping
            # both keeps the overlap as non-manifold material and the remainder
            # as an unfilled hole (docs/algorithm.md §7.1). The repair falls back
            # to the kernel's own general boolean and fuses the two disagreeing
            # junctions into one solid, then interfaces are resolved again —
            # the merged piece presents one agreed region, so nothing is
            # declined there the second time.
            boundary.pieces, n_fused = fuse_disagreeing_pairs(
                lp, boundary.pieces, iface.mismatched
            )
            iface = resolve_interfaces(lp, interior_nodes, boundary.pieces)
            if iface.mismatched:
                raise ProcessingError(
                    f"{len(iface.mismatched)} cap(s) still disagree after fusing "
                    f"the junctions on either side; the local repair could not "
                    f"reconcile them."
                )
        # Weldability is settled here, while a rejection still only costs an
        # extra solid: once a cap face has been given up, putting it back is an
        # undo, and the caps are given up in `finalize_pieces` on the next line.
        for node, h in weld.unweldable(
            lp, tmesh, interior_set, boundary.pieces, iface.interfaces
        ):
            iface.decline(node, h)
        finalize_pieces(boundary.pieces, iface.interfaces)
        comps = build_components(
            interior_nodes, tpl.volume, boundary.pieces, iface.interfaces
        )
        threshold = lp.t ** 3
        dropped = comps.dropped(threshold)
        keep_labels = set(comps.volumes) - dropped
    if n_fused:
        rl.always(
            f"note: {n_fused} disagreeing cap cluster(s) repaired with a local "
            f"boolean fuse (docs/algorithm.md §7.1); the merged junctions present "
            f"one agreed region rather than an extra solid."
        )
    _log_interfaces(rl, lp, iface)
    stats["interfaces"] = iface.n_pairs
    stats["caps_unpaired"] = len(iface.unpaired)
    stats["caps_mismatched"] = len(iface.mismatched)
    stats["caps_unweldable"] = len(iface.unweldable)
    stats["fused_cap_clusters"] = n_fused
    _log_dropped(rl, comps, dropped, threshold)
    stats["components"] = len(comps.volumes)
    stats["components_dropped"] = len(dropped)

    if not keep_labels:
        raise ProcessingError(
            f"Every connected component fell below the floating-body threshold of "
            f"{threshold:g} mm^3 (specification.md §5), leaving nothing to write."
        )

    kept = _partition_kept(comps, keep_labels, boundary.pieces)
    want_rings = _rings_needed(kept.interior_groups, iface.interfaces)

    # The boundary layer is sewn to itself first, because it is the only part
    # whose pairing is genuinely unknown, and because the interior is then built
    # onto whatever topology it ends up with (docs/algorithm.md §8).
    with Timer(rl, "stitch"):
        boundary_faces, sew_stats = weld.sew_boundary(
            kept.pieces, kept.piece_groups,
            workers=args.workers, tmpdir=tmpdir,
            pool=pool, want_rings=want_rings,
        )
        rings = weld.interface_rings(lp, tmesh, boundary_faces, want_rings)
    if sew_stats.max_worker_rss:
        # Same fold-plus-report pattern as the boundary-trim stage above: folded
        # into the run's high-water mark, and also reported per-stage so a peak
        # can be attributed to trimming or to tiled sewing rather than only to
        # "somewhere in this run" (specification.md §3).
        rl.note_worker_rss(sew_stats.max_worker_rss)
        rl.line(f"peak stitch worker memory: {format_bytes(sew_stats.max_worker_rss)}")
        stats["peak_stitch_worker_memory"] = format_bytes(sew_stats.max_worker_rss)
    rl.line(
        f"boundary layer: {len(kept.pieces)} piece(s) sewn into "
        f"{sum(len(f) for f in boundary_faces.values())} faces across "
        f"{len(boundary_faces)} component(s)"
        + (
            f", tiled into {sew_stats.tiles} tile(s) across "
            f"{sew_stats.tiled_components} component(s)"
            if sew_stats.tiles else ""
        )
        + f"; {len(rings)} interface ring(s) located"
    )
    if sew_stats.repaired_components:
        rl.always(
            f"note: {sew_stats.repaired_components} tiled component(s)' boundary-sew "
            f"round 2 (docs/algorithm.md §8) left a free-edge count other than its "
            f"interior interfaces after the seam-only split and were redone with a "
            f"full unsplit sew; the output is unaffected, only slower for those "
            f"components."
        )
    stats["interior_interfaces"] = len(rings)
    stats["stitch_tiles"] = sew_stats.tiles
    stats["stitch_tiled_components"] = sew_stats.tiled_components
    stats["stitch_repaired_components"] = sew_stats.repaired_components

    with Timer(rl, "instance"):
        interior = build_interior_shell(
            lp, tpl, tmesh, kept.interior_nodes, iface.interfaces,
            groups=kept.interior_groups, adopted=rings,
        )
    istats = interior.stats
    stats.update(istats)
    rl.line(
        f"interior shell: {istats['interior_faces']} faces, "
        f"{istats['interior_vertices']} shared vertices, {istats['interior_edges']} shared edges, "
        f"{istats['interior_open_edges']} open edges"
    )

    with Timer(rl, "assemble"):
        shells = weld.assemble(interior.shells, boundary_faces)
        result_solids, _ = weld.close_shells(shells)
    rl.line(f"assembled {len(result_solids)} watertight solid(s)")

    if len(result_solids) != len(keep_labels):
        raise ProcessingError(
            f"Assembly produced {len(result_solids)} solids where the junction "
            f"graph proves there are {len(keep_labels)} connected components. The "
            f"faces were grouped wrongly, so the output cannot be trusted."
        )

    with Timer(rl, "simplify"):
        result_solids, simplify_stats = _unify(result_solids, pool=pool, tmpdir=tmpdir)
    rl.line(
        f"same-domain unification: {simplify_stats['faces_before']} -> "
        f"{simplify_stats['output_faces']} faces, {simplify_stats['edges_before']} -> "
        f"{simplify_stats['output_edges']} edges "
        f"({100 * (1 - simplify_stats['output_faces'] / max(simplify_stats['faces_before'], 1)):.0f}% fewer); "
        f"volume drift {simplify_stats['volume_drift']:.2e} "
        f"(tolerance {UNIFY_VOLUME_TOL:g})"
    )
    if simplify_stats["unmerged_solids"]:
        rl.always(
            f"note: the geometry kernel declined to unify "
            f"{simplify_stats['unmerged_solids']} of {len(result_solids)} solid(s); "
            f"they are exported as built. The output is larger, not different."
        )
    if simplify_stats["max_worker_rss"]:
        rl.note_worker_rss(simplify_stats["max_worker_rss"])
        rl.line(f"peak simplify worker memory: {format_bytes(simplify_stats['max_worker_rss'])}")
        stats["peak_simplify_worker_memory"] = format_bytes(simplify_stats["max_worker_rss"])
    stats["output_faces"] = simplify_stats["output_faces"]
    stats["output_edges"] = simplify_stats["output_edges"]
    stats["unmerged_solids"] = simplify_stats["unmerged_solids"]

    with Timer(rl, "validate"):
        invalid, total_volume, validate_rss = _validate(result_solids, pool=pool, tmpdir=tmpdir)
    if validate_rss:
        rl.note_worker_rss(validate_rss)
        rl.line(f"peak validate worker memory: {format_bytes(validate_rss)}")
        stats["peak_validate_worker_memory"] = format_bytes(validate_rss)
    if invalid:
        raise ProcessingError(
            f"{len(invalid)} of {len(result_solids)} output solids failed OCCT's "
            f"BRepCheck_Analyzer validity check."
        )
    rl.line(f"validity: all {len(result_solids)} solid(s) pass BRepCheck_Analyzer")
    stats["solids_written"] = len(result_solids)
    stats["lattice_volume_mm3"] = round(total_volume, 4)

    name = part_name(args.input, args.cc, args.t)
    with Timer(rl, "export"):
        shape = result_solids[0] if len(result_solids) == 1 else occ.compound(result_solids)
        occ.write_step(shape, args.output, name)
        rewrite_step_header(
            args.output, name, generation_params_text(args.input, args.cc, args.t)
        )
    if not os.path.isfile(args.output) or os.path.getsize(args.output) == 0:
        raise OutputError(f"Output STEP file was not written or is empty: {args.output}")
    stats["output_size"] = format_bytes(os.path.getsize(args.output))

    # No round-trip re-import here (removed: docs/algorithm.md §9). It cost
    # 22 m 29 s of CPU re-parsing 2.00 GB of STEP to full B-rep purely to count
    # solids on the docs/specification.md §10 rehearsal — the single most
    # expensive stage in the run — to check something `tools/e2e.py` already
    # checks on every committed scenario in dev/CI (`vg.brepcheck`, a real
    # `STEPControl_Reader` round trip). Paying that cost again on every
    # production run of a part this size bought back only the reassurance that
    # export and the filesystem did not silently corrupt the file, which the
    # "written and non-empty" check just above already covers for the write
    # failure this tool can actually cause and correct for (exit 6).

    stats["output"] = args.output
    shutil.rmtree(tmpdir, ignore_errors=True)
    return stats


UNIFY_VOLUME_TOL = 1e-5
"""Relative tolerance on the volume that same-domain unification must preserve.

Unification re-describes the boundary without moving it, so the only expected
drift is quadrature noise: merged faces are larger, more complex trimmed regions,
and OCCT integrates them slightly differently. Calibrated against measurement —
on purely planar geometry, where the volume is known analytically, the drift is
**1.9e-15**, i.e. exact; the drift only appears on boundary solids carrying
curved trimmed faces, at up to ~2e-7 relative on ``dense-lattice``. The bar sits
roughly fifty times above that, and still well below ``t³``, the smallest volume
this tool considers meaningful at all.

The real guards against a bad merge are the solid-count check beside this one and
the ``BRepCheck_Analyzer`` gate immediately after: merging two faces that are not
the same surface distorts the boundary, which shows up as an invalid solid long
before it shows up as a volume this close to unchanged. Every run logs the
observed drift, so the margin this bar actually has is visible rather than
assumed.
"""


def _unify_one(solid: TopoDS_Shape) -> tuple[TopoDS_Shape, bool]:
    """Unify one solid, giving up on the parts of the job that will not run.

    Returns the solid and whether the kernel ran at all, so the caller reports
    what happened rather than inferring it from the face count. Those are not
    the same thing: an already-minimal solid unifies *successfully* and changes
    nothing — docs/algorithm.md §9 records the junction template doing exactly
    that, 30 faces to 30 — so "nothing changed" must never be read as "the
    kernel refused".

    Unification is a size optimization, not a correctness step, so a kernel that
    refuses to perform it must not end the run — spec §11's principle that "do
    more work" is an acceptable failure mode while refusing sound input is not.
    The output is simply described more verbosely, and every downstream gate
    still applies to it.

    Two rungs, because the failure is not all-or-nothing. Merging *edges*
    (concatenating the collinear pairs left inside a merged wire) is the part
    that throws, and it is nearly worthless here: measured on the 80 mm ball at
    ``cc=10, t=1``, dropping it still merges 20,268 faces down to 10,554 and
    81,816 edges down to 62,376, while OCCT's edge pass on its own removes 4
    edges out of 81,816 — and, with face merging enabled alongside it, raises
    ``Standard_Failure: Courbes non jointives`` instead.
    """
    try:
        return occ.unify_same_domain(solid), True
    except Standard_Failure:
        pass
    try:
        return occ.unify_same_domain(solid, unify_edges=False), True
    except Standard_Failure:
        return solid, False


def _check_unify_result(pre_volume: float, post_volume: float, n_produced: int) -> float:
    """The two guards every unified solid must clear, wherever it ran.

    Both guards are load-bearing regardless of which stage produced the
    numbers (docs/algorithm.md §9): the solid-count check also protects the
    junction-graph cross-check that precedes assembly, and the volume-drift
    bar is what actually catches a merge across faces that were not genuinely
    the same surface. Shared by the serial and parallel paths so the bars —
    and their error text — cannot drift apart between them.
    """
    if n_produced != 1:
        raise ProcessingError(
            f"Same-domain unification turned one solid into {n_produced}. "
            f"It must only re-describe the boundary, never re-partition the body."
        )
    drift = abs(post_volume - pre_volume) / max(abs(pre_volume), 1.0)
    if drift > UNIFY_VOLUME_TOL:
        raise ProcessingError(
            f"Same-domain unification changed a solid's volume from "
            f"{pre_volume:.6f} to {post_volume:.6f} mm^3 ({drift:.2e} relative, "
            f"tolerance {UNIFY_VOLUME_TOL:g}). Faces that are not genuinely the "
            f"same surface were merged, so the boundary moved."
        )
    return drift


def _unify_serial(solids: list[TopoDS_Shape]):
    """Compact each solid's B-rep on the master, one at a time.

    Each solid is unified on its own rather than as one compound: it keeps a
    1:1 mapping so :func:`_check_unify_result`'s count guard is exact.
    """
    faces_before = edges_before = 0
    faces_after = edges_after = 0
    worst_drift = 0.0
    skipped = 0
    out: list[TopoDS_Shape] = []
    for solid in solids:
        f, e = occ.count_subshapes(solid)
        faces_before += f
        edges_before += e
        pre_volume = occ.volume(solid)

        merged, ran = _unify_one(solid)
        if not ran:
            skipped += 1
        produced = occ.solids(merged)
        post_volume = occ.volume(produced[0]) if len(produced) == 1 else float("nan")
        worst_drift = max(worst_drift, _check_unify_result(pre_volume, post_volume, len(produced)))

        f, e = occ.count_subshapes(produced[0])
        faces_after += f
        edges_after += e
        out.append(produced[0])

    return out, {
        "faces_before": faces_before,
        "edges_before": edges_before,
        "output_faces": faces_after,
        "output_edges": edges_after,
        "volume_drift": worst_drift,
        "unmerged_solids": skipped,
        "max_worker_rss": 0,
    }


def _worker_unify(job):
    """Unify one solid in a worker process — the same ladder as :func:`_unify_one`.

    Only paths and small plain data cross the process boundary: the produced
    solid is written back to ``out_path`` and everything :func:`_check_unify_result`
    needs is returned as scalars, so the master never has to guess at what
    happened inside the worker.
    """
    in_path, out_path = job
    solid = _read_brep(in_path)
    faces_before, edges_before = occ.count_subshapes(solid)
    pre_volume = occ.volume(solid)

    merged, ran = _unify_one(solid)
    produced = occ.solids(merged)
    n_produced = len(produced)
    if n_produced == 1:
        post_volume = occ.volume(produced[0])
        faces_after, edges_after = occ.count_subshapes(produced[0])
        _write_brep(produced[0], out_path)
    else:
        # The master raises on this before touching out_path, so nothing else
        # here needs to be meaningful — just present, so the tuple shape holds.
        post_volume = float("nan")
        faces_after = edges_after = 0

    from .runlog import peak_rss_bytes

    return (
        out_path, faces_before, edges_before, faces_after, edges_after,
        pre_volume, post_volume, n_produced, ran, peak_rss_bytes(),
    )


def _unify_parallel(solids: list[TopoDS_Shape], pool: WorkerPool, tmpdir: str):
    """Compact each solid's B-rep across the shared worker pool.

    Dispatch, guards and the log-facing stats are identical to
    :func:`_unify_serial`'s — only *where* :func:`_unify_one` runs differs — so
    the two paths cannot silently diverge in what they consider a pass.
    """
    jobs = []
    sizes: dict[str, int] = {}
    for i, solid in enumerate(solids):
        in_path = os.path.join(tmpdir, f"unify_{i}.brep")
        _write_brep(solid, in_path)
        jobs.append((in_path, os.path.join(tmpdir, f"unify_{i}_out.brep")))
        # Cheap: this shape is already in memory. Used only to dispatch the
        # biggest job first — the rehearsal's 14 solids are very unequal
        # (specification.md §10), and `imap` assigns queued jobs to workers
        # as they idle, so a large job left late in submission order can only
        # start once an unrelated small job frees a worker (WorkerPool.run).
        sizes[in_path] = occ.count_subshapes(solid)[0]

    results, max_rss = pool.run(_worker_unify, jobs, sort_by=lambda job: sizes[job[0]])

    faces_before = edges_before = 0
    faces_after = edges_after = 0
    worst_drift = 0.0
    skipped = 0
    out: list[TopoDS_Shape] = []
    for (out_path, fb, eb, fa, ea, pre_volume, post_volume, n_produced, ran, _rss) in results:
        faces_before += fb
        edges_before += eb
        if not ran:
            skipped += 1
        worst_drift = max(worst_drift, _check_unify_result(pre_volume, post_volume, n_produced))
        faces_after += fa
        edges_after += ea
        out.append(_read_brep(out_path))

    return out, {
        "faces_before": faces_before,
        "edges_before": edges_before,
        "output_faces": faces_after,
        "output_edges": edges_after,
        "volume_drift": worst_drift,
        "unmerged_solids": skipped,
        "max_worker_rss": max_rss,
    }


def _unify(solids: list[TopoDS_Shape], *, pool: WorkerPool | None = None, tmpdir: str | None = None):
    """Compact every result solid's B-rep, verifying the geometry is unchanged.

    Dispatched across ``pool`` when one is active and there is more than one
    solid to spread across it; a single solid, or no pool, runs on the master
    exactly as before this existed (docs/specification.md §10 path 1). G7
    (`tools/prototypes/RESULTS.md`) measured that OCP holds the GIL around
    ``ShapeUpgrade_UnifySameDomain``, so this is a process pool with a `.brep`
    round-trip, the same mechanism boundary trimming and the boundary sew use
    — not threads, which showed no real speedup.
    """
    if pool is not None and pool.active and tmpdir is not None and len(solids) >= 2:
        return _unify_parallel(solids, pool, tmpdir)
    return _unify_serial(solids)


def _worker_validate(path: str):
    """Validate one solid in a worker process."""
    solid = _read_brep(path)
    valid = occ.is_valid(solid)
    volume = occ.volume(solid)

    from .runlog import peak_rss_bytes

    return valid, volume, peak_rss_bytes()


def _validate(
    solids: list[TopoDS_Shape], *, pool: WorkerPool | None = None, tmpdir: str | None = None
) -> tuple[list[int], float, int]:
    """Run ``BRepCheck_Analyzer`` and sum the volume of every result solid.

    Same dispatch rule as :func:`_unify`: a shared, active pool with more than
    one solid to spread across it and somewhere to stage ``.brep`` files, or
    the master, unchanged from before this existed. G7 found no thread
    speedup here either (`tools/prototypes/RESULTS.md`), so this reuses the
    same process-pool mechanism rather than introducing a second one.

    Returns the indices of any invalid solids, the total volume, and the
    highest peak worker RSS observed (0 on the serial path, folded into the
    run's high-water mark by the caller exactly like every other stage's
    worker RSS).
    """
    if pool is not None and pool.active and tmpdir is not None and len(solids) >= 2:
        paths = []
        sizes: dict[str, int] = {}
        for i, solid in enumerate(solids):
            path = os.path.join(tmpdir, f"validate_{i}.brep")
            _write_brep(solid, path)
            paths.append(path)
            sizes[path] = occ.count_subshapes(solid)[0]  # dispatch order only; see _unify_parallel
        results, max_rss = pool.run(_worker_validate, paths, sort_by=lambda p: sizes[p])
        invalid = [i for i, (valid, _v, _rss) in enumerate(results) if not valid]
        total_volume = sum(v for _valid, v, _rss in results)
        return invalid, total_volume, max_rss
    invalid = [i for i, s in enumerate(solids) if not occ.is_valid(s)]
    total_volume = sum(occ.volume(s) for s in solids)
    return invalid, total_volume, 0


@dataclass
class _Kept:
    """The surviving graph vertices, split by kind and tagged with a component."""

    interior_nodes: np.ndarray
    interior_groups: dict[tuple[int, int, int], int]
    pieces: list
    piece_groups: list[int]


def _partition_kept(comps, keep_labels: set[int], boundary_pieces) -> _Kept:
    """Split surviving graph vertices back into interior nodes and boundary pieces.

    The component label travels with each one, because that is what the assembly
    groups faces by: two components share no interface, so grouping while
    building is exact and free, where rediscovering the split afterwards would be
    a geometric search over every face in the output.
    """
    keep_interior: list[tuple[int, int, int]] = []
    interior_groups: dict[tuple[int, int, int], int] = {}
    keep_pieces = []
    piece_groups: list[int] = []
    for i, v in enumerate(comps.vertices):
        label = int(comps.labels[i])
        if label not in keep_labels:
            continue
        if v.is_interior:
            keep_interior.append(v.node)
            interior_groups[v.node] = label
        else:
            keep_pieces.append(boundary_pieces[v.piece_index])
            piece_groups.append(label)
    nodes_arr = (
        np.array(keep_interior, dtype=np.int64)
        if keep_interior
        else np.empty((0, 3), dtype=np.int64)
    )
    return _Kept(nodes_arr, interior_groups, keep_pieces, piece_groups)


def _rings_needed(interior_groups, interfaces) -> dict:
    """``(node, half-strut) -> component`` for each cap the interior adopts.

    Interior-to-interior caps are excluded: both sides come out of one index and
    already share their topology, so there is nothing to adopt. At rehearsal
    scale that is the difference between locating 18 thousand rings and 97
    thousand.
    """
    steps = [tuple(int(x) for x in neighbor_step(h)) for h in range(6)]
    out = {}
    for node, group in interior_groups.items():
        for h in range(6):
            if (node, h) not in interfaces:
                continue
            step = steps[h]
            if (node[0] + step[0], node[1] + step[1], node[2] + step[2]) not in interior_groups:
                out[(node, h)] = group
    return out


def _log_interfaces(rl: RunLog, lp, iface) -> None:
    """Report the interface tally, and anything the two sides disagreed about.

    The two anomaly counts are the diagnostic this pipeline previously lacked.
    The old one-sided rule dropped both of these categories regardless — leaving,
    respectively, a hole with nothing behind it and two holes that do not match —
    and either way the consequence surfaced only much later, as an unclosed shell
    out of the stitcher with no indication of where it was. Sample world positions
    are logged so the region can be found in the part.
    """
    rl.line(f"interfaces: {iface.n_pairs} cap(s) stitched across")
    if iface.unpaired:
        sample = ", ".join(
            f"{node}h{h}@{np.round(lattice_node(lp, node), 3).tolist()}"
            for node, h in iface.unpaired[:10]
        )
        rl.always(
            f"note: {len(iface.unpaired)} cap(s) face a junction that produced no "
            f"matching cap, so they are kept as exterior surface rather than opened "
            f"as an interface. This is a boolean disagreeing with itself across a "
            f"shared face; the output stays watertight. Samples: {sample}"
            + (" ..." if len(iface.unpaired) > 10 else "")
        )
    if iface.mismatched:
        sample = ", ".join(
            f"{node}h{h} {a:.6g} vs {b:.6g} mm^2" for node, h, a, b in iface.mismatched[:10]
        )
        rl.always(
            f"note: {len(iface.mismatched)} cap(s) were presented by both sides but "
            f"with regions that disagree beyond {CAP_AREA_REL_TOL:g} relative, so they "
            f"are not stitched across. Each leaves an extra solid in the output rather "
            f"than a hole. Samples: {sample}"
            + (" ..." if len(iface.mismatched) > 10 else "")
        )
    if iface.unweldable:
        sample = ", ".join(
            f"{node}h{h}@{np.round(lattice_node(lp, node), 3).tolist()}"
            for node, h in iface.unweldable[:10]
        )
        rl.always(
            f"note: {len(iface.unweldable)} cap(s) have holes on the two sides that "
            f"do not correspond edge for edge, so they are not welded. Each leaves an "
            f"extra solid in the output rather than a hole. Samples: {sample}"
            + (" ..." if len(iface.unweldable) > 10 else "")
        )


def _log_dropped(rl: RunLog, comps, dropped: set[int], threshold: float) -> None:
    """Report floating-body removals as one aggregate line, never one per body.

    A line per removed solid turns a run that drops many of them into a
    multi-thousand-line tail, which buries the rest of the log and slows the run
    down for no gain. One aggregate line carries the same information.
    """
    if not dropped:
        rl.line(f"floating bodies: none below the {threshold:g} mm^3 threshold")
        return
    vols = sorted(comps.volumes[c] for c in dropped)
    members = sum(len(comps.members[c]) for c in dropped)
    sample = ", ".join(f"{v:.4g}" for v in vols[:20])
    rl.always(
        f"Dropped {len(dropped)} floating body/bodies ({members} junction pieces) "
        f"below the {threshold:g} mm^3 threshold: total {sum(vols):.4f} mm^3, "
        f"min {vols[0]:.4g}, max {vols[-1]:.4g}; volumes: {sample}"
        + (" ..." if len(vols) > 20 else "")
    )

