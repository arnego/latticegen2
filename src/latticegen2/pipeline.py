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
from .classify import Klass, classify_nodes, stage_mesh, tessellate_surface
from .cli import Args
from .connect import build_components
from .errors import InputGeometryError, OutputError, ProcessingError
from .interior import build_interior_shell, extract_template_mesh
from .junction import build_template
from .lattice import candidate_nodes, lattice_params, neighbor_step, part_name
from .lattice import node as lattice_node
from .parallel import WorkerPool, set_thread_budget
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
    rl.tmpdir = tmpdir
    try:
        stats = _run(args, rl, tmpdir)
    except BaseException:
        rl.always(f"Intermediate files kept for analysis in: {tmpdir}")
        raise
    # Cleared only on the success path, and only after `_run` has actually
    # removed it, so the attribute means "there is a folder there" rather than
    # "there was one".
    rl.tmpdir = None
    return stats


def _run(args: Args, rl: RunLog, tmpdir: str) -> dict:
    occ.quiet_kernel()
    # OCCT's own native thread pool sizes itself to the machine, not to what
    # this run was told it may use, so `--cores 2` on a six-core box would
    # otherwise let `validate` launch six threads (docs/algorithm.md §9,
    # specification.md §3). Set once here, on the master, which is the only
    # process that asks OCCT for threads.
    set_thread_budget(args.workers)
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
        seam_gap = occ.surface_seam_gaps(body, bar=SEAM_GAP_WARN_FRACTION * args.t)
    rl.line(f"input bounding box: lo={np.round(lo, 4).tolist()} hi={np.round(hi, 4).tolist()}")
    _report_seam_gaps(rl, seam_gap, args.t, stats)

    body_path = os.path.join(tmpdir, "input.brep")
    BRepTools.Write_s(body, body_path)

    with Timer(rl, "tessellate"):
        mesh = tessellate_surface(body, lp)
    # Staged here rather than by `classify`, which needs the same file but only
    # on its parallel path: `boundary` needs it either way, to ask whether a
    # trim left material outside the body without asking a boolean
    # (docs/algorithm.md §7.2). `classify_nodes` is handed the path so the two
    # stages share one staging write rather than making the same one twice.
    mesh_path = stage_mesh(mesh, tmpdir)
    rl.line(
        f"surface mesh: {len(mesh.tris)} triangles, {len(mesh.verts)} vertices, "
        f"measured chordal deviation d={mesh.deviation:.5f} mm "
        f"(classification margin r+d={lp.r + mesh.deviation:.5f} mm)"
    )

    # One worker pool for the whole rest of the run (docs/algorithm.md §8, §12),
    # shared by every parallel stage from here on — classification, boundary
    # trim, boundary-sew, same-domain unification and validation — instead of
    # each stage building and tearing down its own. Built even when
    # `args.workers <= 1`: `WorkerPool` is then inert (`.active` is False) and
    # every stage below takes its own sequential path exactly as if no pool
    # existed, so this costs nothing when there is nothing to parallelise.
    # Entered and exited explicitly rather than as a `with` block wrapping the
    # rest of this function, purely to avoid re-indenting the whole
    # classify-through-validate span; the effect is identical, and a Ctrl+C here
    # is still caught at the top level as `CancelledError` after the pool tears
    # its workers down (docs/algorithm.md §10) — `Pool.__exit__` terminates
    # unconditionally, so it does not matter that no real exception info is
    # threaded through this `finally`.
    #
    # This is built *before* `classify` rather than after it because
    # `classify_nodes` is itself dispatched across it: the sweep decides every
    # node independently of every other, and it was measured single-threaded at
    # 126 s of the `cc=5, t=1` rehearsal (specification.md §10). Nothing OCCT
    # crosses this stage's process boundary — the mesh and the node indices are
    # plain arrays — so it is the one parallel stage here that needs neither
    # G7's GIL finding nor G15's identity argument.
    pool = WorkerPool(args.workers)
    pool.__enter__()
    try:
        with Timer(rl, "classify"):
            candidates = candidate_nodes(lp, lo, hi)
            classification = classify_nodes(
                lp, mesh, candidates, tmpdir=tmpdir, mesh_path=mesh_path,
                pool=pool, report=rl.substage,
            )
        counts = classification.counts()
        rl.line(
            f"classification: {counts['candidates']} candidate nodes -> "
            f"interior {counts['interior']}, boundary {counts['boundary']}, "
            f"outside {counts['outside']}"
        )
        stats.update({f"nodes_{k}": v for k, v in counts.items()})
        if classification.max_worker_rss:
            rl.note_worker_rss(classification.max_worker_rss)
            rl.line(
                f"peak classify worker memory: "
                f"{format_bytes(classification.max_worker_rss)}"
            )
            stats["peak_classify_worker_memory"] = format_bytes(
                classification.max_worker_rss
            )

        interior_nodes = classification.of(Klass.INTERIOR)
        boundary_nodes = classification.of(Klass.BOUNDARY)
        if len(interior_nodes) == 0 and len(boundary_nodes) == 0:
            raise ProcessingError(
                "No lattice nodes fall inside the input geometry. The volume is likely "
                "smaller than one lattice cell at these parameters — try a smaller -cc."
            )

        return _run_with_pool(
            args, rl, tmpdir, lp, tpl, tmesh, body, body_path, stats,
            interior_nodes, boundary_nodes, pool, mesh_path, mesh.deviation,
            seam_gap,
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
    mesh_path: str,
    mesh_deviation: float,
    seam_gap,
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
            # Two consumers with different appetites. The front-end wants every
            # tick it can get and is rate-limited by `RunLog.substage` itself;
            # the log wants ten lines for the whole stage, whatever its size.
            rl.substage("boundary trim", done, total)
            step = max(1, total // 10)
            if done >= total or done - reported[0] >= step:
                reported[0] = done
                rl.line(f"  boundary trim: {done}/{total}")

        boundary = trim_boundary(
            lp, tpl, boundary_nodes, body, body_path, tmpdir,
            workers=args.workers, progress=progress,
            pool=pool, mesh_path=mesh_path,
            # Twice the mesh's own measured deviation: below `d` the mesh
            # cannot tell inside from outside at all (docs/algorithm.md §5.1),
            # and doubling it leaves the tie-vs-protrusion separation this
            # relies on with an order of magnitude to spare — measured
            # 0.004-0.007 mm on sound junctions against 0.47-0.69 mm on
            # protruding ones (§7.2).
            outside_margin=2.0 * mesh_deviation,
        )
    rl.line(
        f"boundary trim: {len(boundary.pieces)} pieces from {len(boundary_nodes)} "
        f"junctions ({boundary.n_empty} produced no geometry)"
    )
    stats["boundary_pieces"] = len(boundary.pieces)
    stats["boundary_empty"] = boundary.n_empty
    if boundary.n_pinhole_wires:
        # One aggregate line, never one per junction (docs/algorithm.md §7).
        rl.line(
            f"pinhole wires removed: {boundary.n_pinhole_wires} zero-area wire(s) "
            f"from {boundary.n_pinhole_junctions} junction(s) grazing the input surface"
        )
        stats["pinhole_wires_removed"] = boundary.n_pinhole_wires
        stats["pinhole_junctions_repaired"] = boundary.n_pinhole_junctions
    _report_tolerance_ratio(rl, boundary, stats)
    if boundary.n_retrimmed_junctions:
        # `always`, not `line`: this one records the kernel contradicting
        # itself about a junction it was asked to trim, and a run where that
        # happens should say so without -v (docs/algorithm.md §7.2).
        rl.always(
            f"note: {boundary.n_retrimmed_junctions} junction(s)' intersection "
            f"with the input body came back untrimmed and were redone per "
            f"half-strut (docs/algorithm.md §7.2). Without that repair each "
            f"would have left a whole junction's material outside the body."
        )
        stats["retrimmed_junctions"] = boundary.n_retrimmed_junctions
    if boundary.n_localized_junctions:
        rl.always(
            f"note: {boundary.n_localized_junctions} junction(s) left material "
            f"outside the input body and were re-trimmed against a local block "
            f"of it (docs/algorithm.md §7.2)."
        )
        stats["localized_junctions"] = boundary.n_localized_junctions
    if boundary.n_dropped_junctions:
        # Loud rather than quiet: this is the one place the tool knowingly
        # leaves lattice out, and specification.md §1's containment is why.
        rl.always(
            f"note: {boundary.n_dropped_junctions} junction(s) were discarded "
            f"because no intersection available to this stage kept them inside "
            f"the input body; the worst reached {boundary.worst_outside_mm:.4f} mm "
            f"outside at {boundary.worst_outside_node} (docs/algorithm.md §7.2). "
            f"The output is that much smaller and stays within the boundary, "
            f"which specification.md §1 requires and material outside would not."
        )
        stats["dropped_junctions"] = boundary.n_dropped_junctions
        stats["worst_outside_mm"] = round(boundary.worst_outside_mm, 6)
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
        n_weld_fused = 0
        n_invented = 0
        if iface.mismatched:
            # Declining a mismatched cap is not a safe degradation on its own:
            # where the two sides present *different* partial regions, keeping
            # both keeps the overlap as non-manifold material and the remainder
            # as an unfilled hole (docs/algorithm.md §7.1). The repair falls back
            # to the kernel's own general boolean and fuses the two disagreeing
            # junctions into one solid, then interfaces are resolved again —
            # the merged piece presents one agreed region, so nothing is
            # declined there the second time.
            boundary.pieces, n_fused, invented = fuse_disagreeing_pairs(
                lp, boundary.pieces, iface.mismatched
            )
            n_invented += len(invented)
            iface = resolve_interfaces(lp, interior_nodes, boundary.pieces)
            if iface.mismatched:
                raise ProcessingError(
                    f"{len(iface.mismatched)} cap(s) still disagree after fusing "
                    f"the junctions on either side; the local repair could not "
                    f"reconcile them."
                )
        # Weldability is settled here, while a rejection still costs only an
        # extra solid or a local fuse: once a cap face has been given up,
        # putting it back is an undo, and the caps are given up in
        # `finalize_pieces` below.
        rejected = weld.unweldable(
            lp, tmesh, interior_set, boundary.pieces, iface.interfaces
        )
        # **Declining is not free where the two sides are coincident**, which is
        # the case this repair exists for. `decline` withdraws both keys and
        # `finalize_pieces` then keeps *both* sides' cap faces — two sheets at
        # the same nominal quad. That costs "an extra solid" only when the
        # decline actually separates the bodies; where the two pieces stay
        # connected through their other interfaces the duplicate ends up inside
        # one solid, as material counted twice. Measured on
        # `TD_HX_rehearsal_test` at `cc=5, t=1`: one declined cap, both sides
        # 7.368093e-01 mm² and planar, and the two mesh edges between them are
        # the only defects in that body this generator is responsible for
        # (tools/prototypes/RESULTS.md G23).
        #
        # It is the same hazard `fuse_disagreeing_pairs` already exists to avoid
        # for a mismatched cap, so it gets the same repair — non-strict, since
        # unlike a mismatched cap this one has a safe fallback: the decline that
        # stood before, which is what happens to any group the fuse cannot merge.
        if rejected:
            boundary.pieces, n_weld_fused, invented = fuse_disagreeing_pairs(
                lp, boundary.pieces, rejected, strict=False
            )
            n_invented += len(invented)
            if n_weld_fused:
                iface = resolve_interfaces(lp, interior_nodes, boundary.pieces)
                if iface.mismatched:
                    boundary.pieces, extra, invented = fuse_disagreeing_pairs(
                        lp, boundary.pieces, iface.mismatched
                    )
                    n_fused += extra
                    n_invented += len(invented)
                    iface = resolve_interfaces(lp, interior_nodes, boundary.pieces)
                rejected = weld.unweldable(
                    lp, tmesh, interior_set, boundary.pieces, iface.interfaces
                )
        for node, h in rejected:
            iface.decline(node, h)
        finalize_pieces(boundary.pieces, iface.interfaces)
        comps = build_components(
            interior_nodes, tpl.volume, boundary.pieces, iface.interfaces
        )
        threshold = lp.t ** 3
        dropped = comps.dropped(threshold)
        keep_labels = set(comps.volumes) - dropped
        _check_component_tolerance(rl, comps, keep_labels, boundary.pieces, stats)
    if n_fused:
        rl.always(
            f"note: {n_fused} disagreeing cap cluster(s) repaired with a local "
            f"boolean fuse (docs/algorithm.md §7.1); the merged junctions present "
            f"one agreed region rather than an extra solid."
        )
    if n_weld_fused:
        rl.always(
            f"note: {n_weld_fused} cap cluster(s) whose two sides do not "
            f"correspond edge for edge were repaired with the same local fuse "
            f"(docs/algorithm.md §8) rather than declined. Declining leaves both "
            f"sides' cap faces in place, which is duplicate material wherever the "
            f"two pieces stay connected by their other interfaces."
        )
    if n_invented:
        rl.always(
            f"note: a local fuse re-tagged {n_invented} cap key(s) that neither "
            f"operand presented (docs/specification.md §10). The output is "
            f"unaffected unless a boundary layer later fails to close, in which "
            f"case this is the first thing to look at."
        )
    _log_interfaces(rl, lp, iface)
    stats["interfaces"] = iface.n_pairs
    stats["caps_unpaired"] = len(iface.unpaired)
    stats["caps_mismatched"] = len(iface.mismatched)
    stats["caps_unweldable"] = len(iface.unweldable)
    stats["fused_cap_clusters"] = n_fused
    stats["fused_unweldable_clusters"] = n_weld_fused
    stats["invented_cap_keys"] = n_invented
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
            pool=pool, want_rings=want_rings, lp=lp,
            max_vertex_tol=occ.SELF_INTERSECT_MAX_VERTEX_TOL_FRACTION * lp.t,
            report=rl.substage,
        )
        rl.substage("locating interface rings", 0, None)
        t_rings = _dt.datetime.now()
        rings = weld.interface_rings(lp, tmesh, boundary_faces, want_rings)
        t_rings = (_dt.datetime.now() - t_rings).total_seconds()
    # `stitch` is one Timer covering five phases with very different characters
    # — round 1 is worker-parallel, everything after it is serial on the master
    # — which is why the profile reports this stage at 1.09 mean cores against a
    # 5.42 peak (docs/profiling-reports.md). Reported per phase so that gap is
    # attributable rather than merely visible, and in particular so the cost of
    # a discarded seam-only attempt (`round2` when `repair` is nonzero) can be
    # read off directly.
    rl.line(
        f"  stitch phases: round1 {sew_stats.t_round1:.1f}s, split "
        f"{sew_stats.t_split:.1f}s, round2 {sew_stats.t_round2:.1f}s, repair "
        f"{sew_stats.t_repair:.1f}s, retolerance {sew_stats.t_retolerance:.1f}s, "
        f"rings {t_rings:.1f}s"
    )
    stats["stitch_round1_s"] = round(sew_stats.t_round1, 2)
    stats["stitch_split_s"] = round(sew_stats.t_split, 2)
    stats["stitch_round2_s"] = round(sew_stats.t_round2, 2)
    stats["stitch_repair_s"] = round(sew_stats.t_repair, 2)
    stats["stitch_retolerance_s"] = round(sew_stats.t_retolerance, 2)
    stats["stitch_rings_s"] = round(t_rings, 2)
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
        # What the repair discarded, and where — see `SewStats.repair_evidence`.
        # The counts say the split was wrong; only the positions say where, and
        # without them a report of this can be acted on no further than "the
        # split failed on a part I do not have".
        for group, want, got_split, got_unsplit, where in sew_stats.repair_evidence:
            rl.always(
                f"  component {group}: expected {want} free edge(s), seam-only "
                f"split gave {got_split}, full unsplit sew gives {got_unsplit}"
            )
            if where:
                rl.always(
                    f"  component {group}: sample positions of the "
                    f"{got_split} edge(s) the split left free: "
                    f"{[p.tolist() for p in where]}"
                    + (" ..." if got_split > len(where) else "")
                )
    if sew_stats.retoleranced_faces or sew_stats.still_invalid_faces:
        rl.always(
            f"vertex tolerances corrected on {sew_stats.retoleranced_faces} sewn "
            f"boundary face(s) (docs/algorithm.md §8); no geometry moved"
            + (
                f". {sew_stats.still_invalid_faces} face(s) remain invalid for "
                f"some other reason and will be reported by the validity gate"
                if sew_stats.still_invalid_faces else ""
            )
        )
    stats["interior_interfaces"] = len(rings)
    stats["stitch_tiles"] = sew_stats.tiles
    stats["stitch_tiled_components"] = sew_stats.tiled_components
    stats["stitch_repaired_components"] = sew_stats.repaired_components
    if sew_stats.retoleranced_faces:
        stats["retoleranced_faces"] = sew_stats.retoleranced_faces
    if sew_stats.still_invalid_faces:
        stats["still_invalid_faces"] = sew_stats.still_invalid_faces

    with Timer(rl, "instance"):
        rl.substage("building the interior shell", 0, None)
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
        rl.substage("assembling shells", 0, None)
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
        result_solids, simplify_stats = _unify(
            result_solids, pool=pool, tmpdir=tmpdir, report=rl.substage
        )
    rl.line(
        f"same-domain unification: {simplify_stats['faces_before']} -> "
        f"{simplify_stats['output_faces']} faces, {simplify_stats['edges_before']} -> "
        f"{simplify_stats['output_edges']} edges "
        f"({100 * (1 - simplify_stats['output_faces'] / max(simplify_stats['faces_before'], 1)):.0f}% fewer); "
        f"volume drift {simplify_stats['volume_drift']:.2e}"
        + (
            f" (past the {UNIFY_VOLUME_TOL:g} pre-filter; the boundary "
            f"displacement it implies was measured and is within "
            f"{UNIFY_MAX_DISPLACEMENT:g} mm)"
            if simplify_stats["volume_drift"] > UNIFY_VOLUME_TOL
            else f" (within the {UNIFY_VOLUME_TOL:g} pre-filter)"
        )
    )
    if simplify_stats["unmerged_solids"]:
        rl.always(
            f"note: the geometry kernel declined to unify "
            f"{simplify_stats['unmerged_solids']} of {len(result_solids)} solid(s); "
            f"they are exported as built. The output is larger, not different."
        )
    if simplify_stats["invalid_merges"]:
        # Unconditional, like the "declined to unify" note above it: a run whose
        # output is larger than it should be has to say why, and this one says
        # the kernel produced something the validity gate would have refused.
        rl.always(
            f"note: unifying {simplify_stats['invalid_merges']} of "
            f"{len(result_solids)} solid(s) produced an invalid result, so the "
            f"un-unified solid was kept instead (docs/algorithm.md S9). The "
            f"output is larger, not different -- and this is the one guard on "
            f"that step the volume bars cannot stand in for."
        )
    if simplify_stats["max_worker_rss"]:
        rl.note_worker_rss(simplify_stats["max_worker_rss"])
        rl.line(f"peak simplify worker memory: {format_bytes(simplify_stats['max_worker_rss'])}")
        stats["peak_simplify_worker_memory"] = format_bytes(simplify_stats["max_worker_rss"])
    stats["output_faces"] = simplify_stats["output_faces"]
    stats["output_edges"] = simplify_stats["output_edges"]
    stats["unmerged_solids"] = simplify_stats["unmerged_solids"]
    stats["invalid_merges"] = simplify_stats["invalid_merges"]

    with Timer(rl, "validate"):
        invalid, total_volume = _validate(result_solids, report=rl.substage)
        if invalid:
            # The established gate speaks first. A body that is already invalid
            # here needs no argument about what the file would do with it, and
            # reporting the export check instead would name the second-order
            # symptom of a fault BRepCheck has already found.
            raise ProcessingError(
                f"{len(invalid)} of {len(result_solids)} output solids failed OCCT's "
                f"BRepCheck_Analyzer validity check."
            )
        _check_export_truth(rl, result_solids, tmpdir, stats, seam_gap, args.t)
    rl.line(f"validity: all {len(result_solids)} solid(s) pass BRepCheck_Analyzer")
    stats["solids_written"] = len(result_solids)
    stats["lattice_volume_mm3"] = round(total_volume, 4)

    name = part_name(args.input, args.cc, args.t)
    with Timer(rl, "export"):
        # One opaque `STEPControl_Writer` call over the whole result, so there is
        # nothing to count — 10.7 % of the rehearsal's clock with no fraction to
        # offer. Saying so beats inventing one (docs/algorithm.md §10).
        rl.substage("writing STEP", 0, None)
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


def _report_tolerance_ratio(rl: RunLog, boundary, stats: dict) -> None:
    """Report which trims lean hardest on their own recorded tolerance.

    Reported, never failed on — :data:`latticegen2.occ.TOLERANCE_FEATURE_RATIO_WARN`
    says why, with the measurement. The quantity is
    :func:`latticegen2.occ.tolerance_feature_ratio` and docs/algorithm.md §7.3
    explains why it is taken here, in the worker, where the junction still has
    a name.
    """
    worst_pieces = boundary.worst_tolerance_pieces()
    if not worst_pieces:
        return
    worst = worst_pieces[0]
    tol, face_area, where = worst.tolerance_evidence
    flagged = sum(
        1 for p in boundary.pieces
        if p.tolerance_ratio > occ.TOLERANCE_FEATURE_RATIO_WARN
    )
    stats["worst_tolerance_ratio"] = f"{worst.tolerance_ratio:.3e}"
    stats["tolerance_ratio_flagged"] = flagged
    rl.line(
        f"tolerance vs feature size: worst {worst.tolerance_ratio:.3e} at junction "
        f"{worst.node} ({tol:.4e} mm recorded on a {face_area:.6f} mm^2 face at "
        f"[{where[0]:.3f}, {where[1]:.3f}, {where[2]:.3f}]); {flagged} of "
        f"{len(boundary.pieces)} piece(s) above {occ.TOLERANCE_FEATURE_RATIO_WARN:g}"
    )
    for piece in worst_pieces[1:]:
        t, a, w = piece.tolerance_evidence
        rl.line(
            f"  next: {piece.tolerance_ratio:.3e} at {piece.node} "
            f"({t:.4e} mm on {a:.6f} mm^2 at [{w[0]:.3f}, {w[1]:.3f}, {w[2]:.3f}])"
        )


def _check_component_tolerance(rl: RunLog, comps, keep_labels, pieces, stats: dict) -> None:
    """Name the bodies whose *whole* description leans on tolerance.

    The sharp predicate, and the reason both halves of §7.3 exist. A grazing
    trim with a fat recorded tolerance is common and usually harmless: welded
    into a body of tens of thousands of mm³ it is one small region described
    loosely, and the exported solid is sound. The same trim alone in a small
    floating component is a body whose *entire* description is that loose — and
    that is the body §9's export cannot write faithfully.

    **What separates them is the fraction, not the worst reading, and the first
    version of this got that wrong.** On `SpiralTest.step` at ``cc=5, t=1``, 79
    of 2,404 pieces clear the warning bar and *both* surviving components
    contain one — the dominant body holds the single worst junction in the whole
    part (4.041e-01 at ``(-4, -8, 2)``). With 2,348 boundary junctions in it, a
    body that large is almost certain to contain a bad one, so a maximum says
    nothing about the body. How much *of* the body is described that way does:

    ======================  =========  ==========  ==========
    component               volume     flagged     fraction
    ======================  =========  ==========  ==========
    0, the lattice proper   27,864 mm³  82 / 2,348  **3.5 %**
    14, the island          4.17 mm³    2 / 6       **33.3 %**
    ======================  =========  ==========  ==========

    **Both figures are reported and neither is a bar.** That is one part, and a
    two-junction component reaches 100 % trivially — every such component here
    is already below ``t^3`` and dropped by §8 before this runs. A tenfold
    separation on a single measurement is a ranking, not a calibration, and
    dressing it up as a threshold would be the mistake docs/specification.md §11
    records four times over.

    So this predicts and :func:`_check_export_truth` decides. What it buys is
    the report arriving at `connect`, about 40 s into a part that takes eight
    minutes, naming junctions of a body the run has not yet built.
    """
    worst_at_node: dict = {}
    for p in pieces:
        node = tuple(p.node)
        if p.tolerance_ratio > worst_at_node.get(node, (0.0,))[0]:
            worst_at_node[node] = (p.tolerance_ratio, p.tolerance_evidence)

    flagged = []
    for cid in sorted(keep_labels):
        readings = []
        for vid in comps.members[cid]:
            vertex = comps.vertices[vid]
            if vertex.is_interior:
                continue        # built by index, exactly; no boolean, no slack
            hit = worst_at_node.get(tuple(vertex.node))
            if hit is not None:
                readings.append((hit[0], tuple(vertex.node), hit[1]))
        if not readings:
            continue
        ratio, node, evidence = max(readings)
        above = sum(1 for r, _n, _e in readings
                    if r > occ.TOLERANCE_FEATURE_RATIO_WARN)
        if above:
            flagged.append((cid, comps.volumes[cid], above, len(readings),
                            ratio, node, evidence))

    stats["tolerance_flagged_components"] = len(flagged)
    if not flagged:
        return
    rl.always(
        f"note: {len(flagged)} output body/bodies contain junctions whose trim "
        f"needed a tolerance comparable to the feature it bounds "
        f"(docs/algorithm.md §7.3). Read the percentage, not the worst reading: "
        f"a large body will contain one, a body largely made of them is where an "
        f"unwritable output comes from. The export-truth gate at `validate` "
        f"decides:"
    )
    for cid, volume, above, total, ratio, node, (tol, face_area, where) in flagged:
        rl.always(
            f"  body {cid}: {volume:.4f} mm^3, {above} of {total} boundary "
            f"junction(s) above the bar ({100.0 * above / max(total, 1):.1f}% of "
            f"the body); worst {ratio:.3e} at {node} ({tol:.4e} mm on a "
            f"{face_area:.6f} mm^2 face at [{where[0]:.3f}, {where[1]:.3f}, "
            f"{where[2]:.3f}])"
        )


SEAM_GAP_WARN_FRACTION = 1e-3
"""When an input body's own seam gap is worth saying out loud, against ``t``.

A thousandth of the strut side: 4e-04 mm at ``t = 0.4``, 1e-03 mm at ``t = 1``.
Nothing is refused on it — see :func:`latticegen2.occ.surface_seam_gaps` for
why the committed parts do not support a bar — but a gap this size relative to
the features the trim is about to cut is the thing to know about *before* an
eight-minute run, rather than from a coordinate at the end of one.

Measured, gap divided by the ``t`` each part is normally run at: the ball
2e-14, the cylinder 3.0e-04, `TD_HX_rehearsal_test` 8.1e-03, `SpiralTest`
4.2e-02. The first two are silent, the last two speak up, and only the last is
actually refused later — which is the honest state of knowledge here.
"""


def _report_seam_gaps(rl: RunLog, gap, t: float, stats: dict) -> None:
    """Say whether the *input*'s own faces meet, before anything is built on it.

    docs/specification.md S11: the one body this project has produced that
    cannot be written to STEP inherited a 42.4 micron gap between two of its
    input's swept B-spline patches. That gap is invisible on the input's own
    1,000 mm2 faces and fatal on the 0.01 mm2 ones the trim cuts out of the same
    region, and it is the reason no re-fitting of a pcurve repairs that body:
    there is no curve on either surface that closes a gap between the surfaces.

    **Reported to the console, and never gated on.** The warning goes out with
    ``console=True`` rather than under ``-v``, because it is the one line in the
    run that is about the user's *file* rather than about this tool's progress,
    and because the failure it predicts arrives minutes later carrying an output
    coordinate that nobody can act on. Both places name a point in the input's
    own coordinate system, which is what a modeller needs to find the seam.

    It is not a gate. The separation between the part that fails and the worst
    part that ships is a factor of five on two data points, and refusing on that
    would refuse `TD_HX_rehearsal_test`, which has been inspected and accepted —
    docs/algorithm.md §11's one unacceptable failure mode. The run continues
    exactly as before.
    """
    if gap.edges == 0:
        return
    stats["input_seam_gap_mm"] = f"{gap.worst:.3e}"
    if gap.worst <= SEAM_GAP_WARN_FRACTION * t:
        rl.line(
            f"input seam gaps: worst {gap.worst:.3e} mm over {gap.edges} shared "
            f"edge(s) -- negligible against t={t:g} mm"
        )
        return
    rl.line(
        f"input seam gaps: the two faces meeting at one of this body's own edges "
        f"are {gap.worst:.4e} mm apart at [{gap.where[0]:.3f}, {gap.where[1]:.3f}, "
        f"{gap.where[2]:.3f}] -- {gap.worst / t:.2e} of t={t:g} mm, on {gap.over} "
        f"of {gap.edges} shared edge(s). That is a property of the input file, "
        f"not of this run: its own tolerance covers it, and the trim is about to "
        f"cut that region into features smaller than the gap. If this run fails "
        f"at `export truth`, that coordinate is the geometry to fix",
        console=True,
    )


def _check_export_truth(rl: RunLog, solids: list[TopoDS_Shape], tmpdir: str,
                        stats: dict, seam_gap=None, t: float = 0.0) -> None:
    """Can each output body survive being written to STEP? Measured, per body.

    `BRepCheck_Analyzer` asks whether a shape agrees with itself to within the
    tolerances **recorded in this process**. STEP AP214 cannot carry those: one
    ``UNCERTAINTY_MEASURE_WITH_UNIT`` per file, against one per vertex, edge and
    face in the B-rep. So a solid can pass every gate this pipeline has and
    still not describe, in the file, what it describes here — and the pipeline's
    own repairs manufacture exactly that geometry, docs/algorithm.md S8's second
    rung being a widened vertex tolerance whose safety rests on moving nothing.

    **The instrument is a tessellation of the exported body**, not a proxy for
    one (:func:`latticegen2.occ.exported_mesh_defects`), and three cheaper
    proxies were tried and disproved against ground truth on 16 real bodies
    before settling here: a fault count after the round trip false-positives on
    two sound rehearsal solids, and so does the share of a body's surface
    described more loosely than its own feature size. Both would refuse, or
    delete, geometry from a part that has been inspected and accepted. Only
    "does it still tessellate" matches what actually breaks downstream, and it
    is the symptom that found the one genuinely unwritable body this project has
    produced.

    **Every solid is measured, with no size bound, and that is the expensive
    decision on this branch.** An earlier revision skipped solids above a face
    count and reported them *unmeasured* — honest, but it put the dominant body
    of every production part outside the only detector with a clean record,
    which is precisely where an unwritable description would do the most damage.
    Per the user's decision the cost is accepted instead. docs/algorithm.md S9
    removed a whole-output re-import for costing 22 minutes, and this is that
    cost returning knowingly, for a correctness check rather than for a count of
    solids. The rehearsal-scale figure is **not yet measured**.

    The pcurve reading below is kept and logged because it is cheap, exact and
    informative about *why* a body is fragile - but it decides nothing, and
    :data:`latticegen2.occ.LOOSE_AREA_FRACTION_MAX` records the measurement that
    took it out of the deciding seat.
    """
    started = _dt.datetime.now()
    broken: list[tuple[int, int, int, occ.Refinement]] = []
    refined: list[tuple[int, int, float]] = []
    worst = 0.0
    fraction = 0.0
    loose_faces = 0
    pairs = 0
    probe = os.path.join(tmpdir, "export_truth_probe.step")
    for i, solid in enumerate(solids):
        reading = occ.curve_on_surface_deviations(solid)
        pairs += reading.pairs
        loose_faces += reading.loose_faces
        worst = max(worst, reading.worst)
        fraction = max(fraction, reading.loose_area_fraction)
        n_faces, _ = occ.count_subshapes(solid)
        rl.substage("export truth", i, len(solids))
        defects = occ.exported_mesh_defects(solid, probe)
        triangles, bad = defects.triangles, defects.bad
        rl.line(
            f"  solid {i}: {n_faces} faces -> {triangles} triangle(s) after a "
            f"STEP round trip, {bad} non-manifold edge(s)"
            # Unconditional, not folded into the `if bad:` branch below: a solid
            # reading 0 *because* triangles were skipped is precisely what a
            # reader must be able to see.
            + (f" ({defects.degenerate} degenerate triangle(s) skipped)"
               if defects.degenerate else "")
            + f"; {reading.loose_area_fraction:.4e} of its surface loosely "
            f"described, worst deviation {reading.worst:.4e} mm"
        )
        if bad:
            # The breakdown by use count and the positions, because a total
            # alone cannot be acted on: one use is a hole or two faces
            # discretizing a shared edge differently, more than two is
            # duplicate material, and they have different causes
            # (docs/specification.md §10).
            rl.line(
                f"    by use count {dict(sorted(defects.by_use.items()))}; "
                f"positions {[list(p) for p in defects.where]}"
            )
            # `exported_mesh_defects` always re-measures a body with
            # readings, so `None` here would mean the instrument did not run.
            # That is refused rather than crashed on: an unmeasured body must
            # never reach the *clean* branch, and it must not take the run down
            # with it either.
            fine = defects.refinement or occ.Refinement(
                False, [], len(defects.implicated), 0,
                "unmeasured: the body was never re-measured")
            ladder = " ".join(f"{d}:{n}" for d, n in fine.counts) or "none"
            rl.line(
                f"    on {fine.core_faces} face(s); re-measured over a "
                f"{fine.extract_faces}-face neighbourhood at {ladder} "
                f"-- {fine.reason}"
            )
            if fine.resolved:
                # Deliberately *not* the sentence a solid reading 0 at the
                # coarse ruler prints. This body did carry readings; what the
                # ladder establishes is that they were the ruler and not the
                # geometry, and a reader has to be able to tell the two runs
                # apart afterwards.
                refined.append((i, bad, fine.counts[-1][0]))
            else:
                broken.append((i, bad, triangles, fine))
    rl.substage("export truth", len(solids), len(solids))
    elapsed = (_dt.datetime.now() - started).total_seconds()
    stats["export_truth_s"] = round(elapsed, 2)
    stats["worst_pcurve_deviation_mm"] = f"{worst:.3e}"
    stats["loose_area_fraction"] = f"{fraction:.3e}"
    stats["loose_faces"] = loose_faces
    if refined:
        # Greppable afterwards, because "cleared by refinement" and "was clean"
        # are different runs and a summary that conflates them would hide the
        # only bodies this second pass has ever been asked about.
        stats["export_truth_refined"] = "; ".join(
            f"solid {i}: {bad} reading(s) at "
            f"{occ.DEFAULT_MESH_DEFLECTION} mm, none at {d} mm"
            for i, bad, d in refined
        )
    rl.line(
        f"export truth: all {len(solids)} solid(s) round-tripped and tessellated; "
        f"worst pcurve-vs-3D deviation {worst:.4e} mm over {pairs} edge/face "
        f"pair(s) [{elapsed:.1f}s]"
        + (f"; {len(refined)} solid(s) cleared only after being re-measured at "
           f"a finer deflection" if refined else "")
    )
    if broken:
        first, bad, triangles, fine = broken[0]
        ladder = " -> ".join(f"{n} at {d} mm" for d, n in fine.counts) or "none"
        raise ProcessingError(
            f"Solid {first} does not survive being written to STEP: after a "
            f"round trip its {triangles} triangles carry {bad} edge(s) not used "
            f"by exactly two of them, so nothing downstream can tessellate it "
            f"consistently"
            # The ladder, not just the count. "13 -> 8 -> 4 -> 4, still 4 at
            # 0.0001 mm" says the readings are geometry; "10 -> 13 -> 17" says
            # the body comes apart the more closely it is measured. A bare
            # total says neither, and this project spent two gates and a
            # 93-minute run recovering that distinction by hand.
            + f". Re-measured over a {fine.extract_faces}-face neighbourhood: "
            f"{ladder} ({fine.reason})"
            + (f" ({len(broken)} solid(s) affected)" if len(broken) > 1 else "")
            + f". STEP AP214 declares one tolerance for a whole file where this "
            f"B-rep carries one per subshape, so a body whose validity rests on "
            f"a locally fat tolerance is valid here and not in the file "
            f"(docs/algorithm.md S9). The run stops rather than writing it: "
            f"nothing has been discarded, the temporary folder is kept, and the "
            f"`connect` stage above named the junctions this body was built from."
            + _seam_gap_note(seam_gap, t)
        )


def _seam_gap_note(gap, t: float) -> str:
    """Name the input's own seam gap in an export-truth failure, when it is big.

    The failure above gives a coordinate in the *output*; this gives the reason,
    in the *input*, and it is the difference between a user re-running with
    different parameters and a user fixing their model. `SpiralTest.step` at
    ``cc=5, t=1`` is the whole case: a 4.2406e-02 mm gap between two of its own
    swept patches, cut into 0.01 mm2 faces by a 1 mm strut.
    """
    if gap is None or t <= 0.0 or gap.worst <= SEAM_GAP_WARN_FRACTION * t:
        return ""
    return (
        f" Note also the `import` stage's reading: two faces of the *input* body "
        f"are {gap.worst:.4e} mm apart at their shared edge "
        f"[{gap.where[0]:.3f}, {gap.where[1]:.3f}, {gap.where[2]:.3f}], "
        f"{gap.worst / t:.2e} of t={t:g} mm. A gap between two surfaces cannot be "
        f"repaired by re-fitting a curve on either of them, so where that region "
        f"is what was trimmed, the input is what has to change."
    )


UNIFY_VOLUME_TOL = 1e-4
"""Relative volume drift above which unification is measured properly.

**A pre-filter, not the bar.** It decides only whether it is worth paying for
the surface area needed to express the drift as a boundary displacement, which
is what :data:`UNIFY_MAX_DISPLACEMENT` then judges. Everything below is the
calibration history of the figure, and it is kept because it is what
established that the drift is a re-description rather than movement -- but
the reason this stayed at 1e-4 rather than being loosened again the next time a
part exceeded it is that loosening it was always the wrong response to a
quantity biased by solid size (:func:`_check_unify_result`).

Relative tolerance on the volume that same-domain unification must preserve.

Unification re-describes the boundary rather than moving it, and the drift is
what that re-description costs: on purely planar geometry, where the volume is
known analytically, it is **1.9e-15** — exact — and it only appears at all on
boundary solids carrying curved trimmed faces, at ~2e-7 relative on
``dense-lattice``.

**The bar was 1e-5 and refused a valid run at 1.381e-05**, on a 181 mm³
floating island of `TD_HX_rehearsal_test` at ``cc=12, t=2.5``. Nothing had
moved: the exact symmetric difference between the two solids, cut both ways, is
**0.000000000 mm³**, and both are ``BRepCheck_Analyzer``-valid. What the number
is measuring is a genuine but sub-tolerance re-description — surface area shifts
too (3.16e-06), and adaptive Gauss-Kronrod integration to a requested 1e-11
does *not* bring the two figures together, so this is not the integrator's
truncation error and cannot be integrated away.

**Read as a displacement it is nothing at all**, and that is the reading this
bar is now set by. ``|ΔV| / surface area`` is the boundary movement the drift
implies: **6.96e-06 mm** on the solid that failed, against the **8.7e-04 to
1.5e-03 mm** tolerances OCCT itself records on the trimmed B-spline faces being
merged (docs/algorithm.md §8, G12). The two descriptions are the same surface
to well within the kernel's own idea of the same surface. At 1e-4, and across
the surface-to-volume ratios these parts produce (1.5–3.7 per mm, measured over
the nine solids of that run), the bar admits at most ~3e-05 to 7e-05 mm of
boundary movement — still more than an order of magnitude inside those face
tolerances, so a merge across genuinely different surfaces cannot hide under it.

That the failing island's mirror twin — same volume and area to 0.1 % — drifts
**29x less** is why no tighter bar is defensible: the magnitude is a property of
which merge the kernel happened to perform, not of the geometry, so it cannot be
predicted from the part.

The real guards against a bad merge remain the solid-count check beside this one
and the ``BRepCheck_Analyzer`` gate immediately after: merging two faces that are
not the same surface distorts the boundary, which shows up as an invalid solid
long before it shows up as a volume this close to unchanged. Every run logs the
observed drift, so the margin this bar actually has is visible rather than
assumed.
"""


UNIFY_MAX_DISPLACEMENT = 1e-3
"""Millimetres. How far same-domain unification may move a solid's boundary.

The quantity :data:`UNIFY_VOLUME_TOL` was always a proxy for, now measured
directly as ``|dV| / surface area`` — see :func:`_check_unify_result` for why
the proxy is biased by solid size and for the pair of solids that showed it,
90x apart in relative drift and both sound.

Calibrated against what the kernel itself records. Worst displacement measured
anywhere: **2.494e-04 mm**, on `SpiralTest.step`'s 4.17 mm^3 island at
``cc=5, t=1``, whose merged edges carry recorded tolerances of 2.1e-02 mm — so
the movement is 85x inside OCCT's own idea of where that surface is. The
rehearsal's 181 mm^3 island at ``cc=12, t=2.5`` measures 6.96e-06 mm, and
`SpiralTest`'s own dominant body 2.752e-06 mm. This bar clears the worst of
them by 4x, sits at or below the 8.7e-04 to 1.5e-03 mm tolerances OCCT records
on the rehearsal's merged faces (docs/algorithm.md §8, G12), and is 400x below
the smallest legal strut (``t = 0.4`` mm, specification.md §3).

The stronger guards are unchanged and are still the ones that would catch a
merge across genuinely different surfaces: the solid-count check beside this
one, and ``BRepCheck_Analyzer`` at `validate`.
"""


def _unify_one(solid: TopoDS_Shape) -> tuple[TopoDS_Shape, bool, bool]:
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

    **The two passes are run separately, and what that now buys is degradation,
    not speed.** ``ShapeUpgrade_UnifySameDomain`` merges coplanar faces and then
    concatenates the collinear edge pairs left inside each merged wire.
    Splitting them was built to let the face merge run with edge merging *off*,
    which the sub-body tiling of docs/specification.md §10 required — and that
    tiling is now disproved (G15: tiles reassemble by shared topology only while
    they stay in one process). The split is kept anyway, for the reason at the
    end of this docstring, and because measurement says it costs nothing:
    `simplify` 13.87 s and 13.94 s split against 13.21 s and 16.18 s combined on
    `dense-lattice`, with `test_pipeline.py` pinning the B-rep it produces as
    identical to the combined call's, faces *and* edges.

    **Dropping the edge pass outright was tried and rejected**, so it is not
    retried. Edge merging is not the near-no-op docs/algorithm.md §9 once
    recorded from the 80 mm ball (4 edges of 81,816): at lattice scale G13
    measured it taking 307,200 edges down to 215,040. Skipping it took
    `simplify` to 9.45 s on `dense-lattice` and handed every second of that back
    to `validate` (6.24 -> 8.21 s) and `export` (6.25 -> 10.50 s), which scale
    with edge count too, for a 35 % larger file (52.80 -> 71.29 MB) and no net
    run-time change (57.28 -> 57.57 s). Face count was identical throughout,
    which is what identifies edges as the whole of the difference.

    What the split is actually worth is here. Edge merging is
    the part that raises ``Standard_Failure: Courbes non jointives`` on geometry
    this tool legitimately produces; now that it runs last and alone, a refusal
    costs only the edge concatenation, where the old ladder threw away a
    completed face merge and paid for a second one.

    **A kernel that does not throw can still hand back an invalid solid, and
    that is degraded from rather than shipped.** docs/algorithm.md §9 already
    names the validity gate as the *stronger* of the guards on this step —
    "merging faces that are not the same surface moves the boundary, which shows
    up as an invalid solid long before it shows up as a changed volume" — but
    until this the run simply *failed* there, which is the one response §11
    forbids for a step whose only job is to make the output smaller. Measured on
    `TD_HX_rehearsal_test` at ``cc=7, t=1.4``: the 279,358-face dominant body is
    `BRepCheck_Analyzer`-valid and tessellates with **0** non-manifold edges,
    unification takes it to 193,721 faces with exactly **one** invalid face of
    3.426439 mm², and the run died at `validate` having discarded the sound
    solid it started from. The volume guards cannot see it — the drift is
    6.20e-06 relative, well inside the 1e-4 pre-filter.

    So the result is checked and the *input* is kept when it does not hold up.
    The failure mode becomes a larger file, exactly as it already is for a
    kernel that throws. The check costs one whole-solid ``is_valid`` per solid,
    in the worker where the solid already lives and parallel across solids —
    23.3 s on that 193,721-face body. It is deliberately **not** compared
    against the input's validity: an input that is itself invalid gains nothing
    from being kept, and `validate` refuses it exactly as it did before, so
    paying 34.7 s to find that out would buy nothing.
    """
    ran = True
    try:
        merged = occ.unify_same_domain(solid, unify_edges=False)
    except Standard_Failure:
        return solid, False, False
    try:
        merged = occ.unify_same_domain(merged, unify_edges=True, unify_faces=False)
    except Standard_Failure:
        pass
    # `parallel=False` because this runs *inside* a worker: the pool around it
    # already occupies the whole `--cores` budget, and OCCT's own thread pool is
    # deliberately left unbounded there (`parallel.set_thread_budget` is not
    # called in the worker initializer), so asking for threads here would be W
    # processes x W threads on W cores -- the over-subscription
    # docs/algorithm.md S9 keeps the *gate* on the master to avoid.
    if not occ.is_valid(merged, parallel=False):
        return solid, ran, True
    return merged, ran, False


def _check_unify_result(
    pre_volume: float, post_volume: float, n_produced: int, area_of=None
) -> float:
    """The two guards every unified solid must clear, wherever it ran.

    Both guards are load-bearing regardless of which stage produced the
    numbers (docs/algorithm.md §9): the solid-count check also protects the
    junction-graph cross-check that precedes assembly, and the volume bar is
    what actually catches a merge across faces that were not genuinely the same
    surface. Shared by the serial and parallel paths so the bars — and their
    error text — cannot drift apart between them.

    **The volume bar is a two-stage test, and only the second stage decides.**
    Relative drift is cheap — both volumes have already been measured — but it
    is not the quantity this guard is about: a merge that moves the boundary by
    a given distance changes a small solid's volume by a much larger *fraction*
    than a big one's, so a relative bar is biased by size rather than by the
    thing it is trying to detect. Measured on one `SpiralTest.step` run, whose
    two solids were unified by the same code in the same call:

    ======================  ===============  ==============
    solid                   relative drift   displacement
    ======================  ===============  ==============
    27,864 mm^3, 64k faces  1.013e-05        2.752e-06 mm
    4.17 mm^3, 53 faces     1.837e-03        2.494e-04 mm
    ======================  ===============  ==============

    Ninety times apart in the first column, and both far inside the 2.1e-02 mm
    tolerances OCCT itself records on the edges it was merging. The small one
    failed a 1e-4 relative bar; its exact symmetric difference against the
    un-unified solid, cut both ways, is **0 mm^3**, and their intersection
    measures the un-unified volume to twelve decimals — the same region, read
    by a quadrature that re-describing the boundary changed.

    So drift is used only as a **pre-filter**, and a solid that trips it is then
    judged on ``|dV| / area``, a displacement in millimetres, which is what
    docs/algorithm.md §9 already argued was the honest reading. ``area_of`` is a
    callable rather than a number because measuring it is expensive — 7.7 s on
    the 64k-face solid above, which at rehearsal scale would put ~70 s on every
    run for a guard that almost never fires. Called only on the path that was
    about to raise, it costs nothing the rest of the time. Without it the
    pre-filter alone decides, which is the old behaviour.
    """
    if n_produced != 1:
        raise ProcessingError(
            f"Same-domain unification turned one solid into {n_produced}. "
            f"It must only re-describe the boundary, never re-partition the body."
        )
    drift = abs(post_volume - pre_volume) / max(abs(pre_volume), 1.0)
    if drift > UNIFY_VOLUME_TOL:
        area = area_of() if area_of is not None else 0.0
        displacement = abs(post_volume - pre_volume) / area if area > 0.0 else float("inf")
        if displacement > UNIFY_MAX_DISPLACEMENT:
            raise ProcessingError(
                f"Same-domain unification changed a solid's volume from "
                f"{pre_volume:.6f} to {post_volume:.6f} mm^3 ({drift:.2e} relative), "
                f"which over its {area:.4f} mm^2 of surface is a boundary "
                f"displacement of {displacement:.3e} mm against a tolerance of "
                f"{UNIFY_MAX_DISPLACEMENT:g} mm. Faces that are not genuinely the "
                f"same surface were merged, so the boundary moved."
            )
    return drift


def _unify_serial(solids: list[TopoDS_Shape], report=None):
    """Compact each solid's B-rep on the master, one at a time.

    Each solid is unified on its own rather than as one compound: it keeps a
    1:1 mapping so :func:`_check_unify_result`'s count guard is exact.
    """
    faces_before = edges_before = 0
    faces_after = edges_after = 0
    worst_drift = 0.0
    skipped = 0
    invalid_merges = 0
    out: list[TopoDS_Shape] = []
    for done, solid in enumerate(solids):
        if report is not None:
            report("unifying solids", done, len(solids))
        f, e = occ.count_subshapes(solid)
        faces_before += f
        edges_before += e
        pre_volume = occ.volume(solid)

        merged, ran, degraded = _unify_one(solid)
        if not ran:
            skipped += 1
        if degraded:
            invalid_merges += 1
        produced = occ.solids(merged)
        post_volume = occ.volume(produced[0]) if len(produced) == 1 else float("nan")
        worst_drift = max(worst_drift, _check_unify_result(
            pre_volume, post_volume, len(produced),
            area_of=lambda: occ.area(produced[0]) if len(produced) == 1 else 0.0,
        ))

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
        "invalid_merges": invalid_merges,
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

    merged, ran, degraded = _unify_one(solid)
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
        # `degraded` goes before the RSS reading, not after: `WorkerPool.run`
        # takes the peak from `rss_index=-1`, so anything appended past it is
        # silently read as a memory figure.
        out_path, faces_before, edges_before, faces_after, edges_after,
        pre_volume, post_volume, n_produced, ran, degraded, peak_rss_bytes(),
    )


def _unify_parallel(solids: list[TopoDS_Shape], pool: WorkerPool, tmpdir: str, report=None):
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

    # Deliberately reported as a *label* with no fraction. Jobs are dispatched
    # largest-first precisely because the solids are wildly unequal — the
    # rehearsal's 14 are one dominant body plus 13 scraps — so "3 of 14 done"
    # after fifteen minutes would be a number that actively misleads.
    if report is not None:
        report(f"unifying {len(jobs)} solid(s), largest first", 0, None)
    results, max_rss = pool.run(_worker_unify, jobs, sort_by=lambda job: sizes[job[0]])

    faces_before = edges_before = 0
    faces_after = edges_after = 0
    worst_drift = 0.0
    skipped = 0
    invalid_merges = 0
    out: list[TopoDS_Shape] = []
    for (out_path, fb, eb, fa, ea, pre_volume, post_volume, n_produced, ran,
         degraded, _rss) in results:
        faces_before += fb
        edges_before += eb
        if not ran:
            skipped += 1
        if degraded:
            invalid_merges += 1
        # Read back only if the pre-filter trips: the shape is loaded a few
        # lines below anyway on the path that does not, so this costs a second
        # read solely on the path that is already about to fail.
        worst_drift = max(worst_drift, _check_unify_result(
            pre_volume, post_volume, n_produced,
            area_of=lambda path=out_path: occ.area(_read_brep(path)),
        ))
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
        "invalid_merges": invalid_merges,
        "max_worker_rss": max_rss,
    }


def _unify(solids: list[TopoDS_Shape], *, pool: WorkerPool | None = None,
           tmpdir: str | None = None, report=None):
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
        return _unify_parallel(solids, pool, tmpdir, report=report)
    return _unify_serial(solids, report=report)


def _validate(solids: list[TopoDS_Shape], report=None) -> tuple[list[int], float]:
    """Run ``BRepCheck_Analyzer`` and sum the volume of every result solid.

    **On the master, in one process, with OCCT's own threads** — not dispatched
    per solid across :class:`WorkerPool`, which is what this did before gate
    G18. That dispatch was docs/specification.md §10's "path 4", kept despite
    measuring *slower* (2 m 59 s -> 3 m 29.6 s) on the reasoning that a part
    with evenly-sized components would benefit even though the rehearsal's 14
    solids — one dominant body plus 13 scraps — could not. G18 changed the
    arithmetic behind that trade in two ways:

    * ``BRepCheck_Analyzer`` has a parallel flag of its own
      (:func:`latticegen2.occ.is_valid`), so the dominant solid is no longer
      stuck on one core just because it is one job. That is the case per-solid
      dispatch could never reach, and it is the case that sets this stage's
      floor on any real part.
    * The two cannot simply be combined. ``--cores`` is "honoured exactly"
      (specification.md §3), and W worker processes each launching W OCCT
      threads is W² threads on W cores. Keeping the dispatch would mean
      bounding each worker to one thread, which gives up exactly the win above.

    Running here also deletes a whole ``.brep`` round trip that existed only to
    reach the workers: the master wrote every solid out and each worker read it
    back, 464 MB each way on the rehearsal, to compute two scalars per solid.
    The solids are already in memory here.

    What this gives up is real and worth stating: on a part whose components
    *are* evenly sized, per-solid dispatch across W processes would beat 1.6x.
    No such part has been measured, and the failure mode is "slower than it
    could have been", never a wrong verdict — docs/algorithm.md §11's rule.

    Returns the indices of any invalid solids and the total volume.
    """
    invalid = []
    for i, solid in enumerate(solids):
        if report is not None:
            report("checking solids", i, len(solids))
        if not occ.is_valid(solid):
            invalid.append(i)
    if report is not None:
        report("checking solids", len(solids), len(solids))
    total_volume = sum(occ.volume(s) for s in solids)
    return invalid, total_volume


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

