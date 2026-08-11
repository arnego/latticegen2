"""Pipeline orchestration.

    parse -> template -> import -> tessellate -> classify
          -> instance interior -> trim boundary -> connect -> sew
          -> validate -> write STEP -> round-trip -> summary

The shape of this is deliberately flat. In the fuse-free architecture there is
no tiling stage, no hierarchical assembly, no distributed merge rounds and no
fuse-based cleanup, because there is nothing left for them to do: the interior is
instanced rather than fused, and connectivity is a graph property rather than a
boolean experiment.
"""

from __future__ import annotations

import datetime as _dt
import os
import shutil

import numpy as np
from OCP.BRepTools import BRepTools
from OCP.TopoDS import TopoDS_Shape

from . import occ
from .boundary import trim_boundary
from .classify import Klass, classify_nodes, tessellate_surface
from .cli import Args
from .connect import build_components
from .errors import InputGeometryError, OutputError, ProcessingError
from .interior import build_interior_shell, extract_template_mesh
from .junction import build_template
from .lattice import candidate_nodes, lattice_params, part_name
from .runlog import RunLog, Timer, format_bytes
from .stepout import generation_params_text, rewrite_step_header, round_trip_check

SEW_TOLERANCE = 1e-6
"""Millimetres. Interfaces are geometrically identical by construction; this
only has to absorb the last-ulp difference between computing a shared cap from
one side of a strut or the other, which is ~1e-14 mm at realistic coordinates.
It stays far below the smallest real feature (``t >= 0.4`` mm), so it can never
weld two genuinely distinct vertices together."""


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

    kept = {(int(a), int(b), int(c)) for a, b, c in np.vstack(
        [x for x in (interior_nodes, boundary_nodes) if len(x)]
    )}

    with Timer(rl, "boundary"):
        def progress(done: int, total: int) -> None:
            if done == total or done % max(1, total // 10) == 0:
                rl.line(f"  boundary trim: {done}/{total}")

        boundary = trim_boundary(
            lp, tpl, boundary_nodes, kept, body, body_path, tmpdir,
            workers=args.workers, background=args.background, progress=progress,
        )
    rl.line(
        f"boundary trim: {len(boundary.pieces)} pieces from {len(boundary_nodes)} "
        f"junctions ({boundary.n_empty} produced no geometry)"
    )
    stats["boundary_pieces"] = len(boundary.pieces)
    stats["boundary_empty"] = boundary.n_empty
    stats["workers"] = args.workers
    if boundary.max_worker_rss:
        rl.line(f"peak worker memory: {format_bytes(boundary.max_worker_rss)}")
        stats["peak_worker_memory"] = format_bytes(boundary.max_worker_rss)

    with Timer(rl, "connect"):
        comps = build_components(interior_nodes, tpl.volume, boundary.pieces)
        threshold = lp.t ** 3
        dropped = comps.dropped(threshold)
        keep_labels = set(comps.volumes) - dropped
    _log_dropped(rl, comps, dropped, threshold)
    stats["components"] = len(comps.volumes)
    stats["components_dropped"] = len(dropped)

    if not keep_labels:
        raise ProcessingError(
            f"Every connected component fell below the floating-body threshold of "
            f"{threshold:g} mm^3 (specification.md §5), leaving nothing to write."
        )

    keep_interior, keep_pieces = _partition_kept(comps, keep_labels, boundary.pieces)
    final_kept = {v.node for i, v in enumerate(comps.vertices) if comps.labels[i] in keep_labels}

    with Timer(rl, "instance"):
        shell, istats = build_interior_shell(lp, tpl, tmesh, keep_interior, final_kept)
        pieces = ([shell] if shell is not None else []) + [p.shell for p in keep_pieces]
    stats.update(istats)
    rl.line(
        f"interior shell: {istats['interior_faces']} faces, "
        f"{istats['interior_vertices']} shared vertices, {istats['interior_edges']} shared edges"
    )
    rl.line(f"assembly input: {len(pieces)} shell(s) to stitch")

    with Timer(rl, "sew"):
        sewn = occ.sew(pieces, SEW_TOLERANCE)
        result_solids = _close_solids(sewn)
    rl.line(f"sewn into {len(result_solids)} solid(s)")

    if len(result_solids) > len(keep_labels):
        raise ProcessingError(
            f"Stitching produced {len(result_solids)} solids where the junction "
            f"graph proves there are {len(keep_labels)} connected components. Some "
            f"interface did not close, which would mean a non-watertight body."
        )

    with Timer(rl, "simplify"):
        result_solids, simplify_stats = _unify(result_solids)
    rl.line(
        f"same-domain unification: {simplify_stats['faces_before']} -> "
        f"{simplify_stats['output_faces']} faces, {simplify_stats['edges_before']} -> "
        f"{simplify_stats['output_edges']} edges "
        f"({100 * (1 - simplify_stats['output_faces'] / max(simplify_stats['faces_before'], 1)):.0f}% fewer); "
        f"volume drift {simplify_stats['volume_drift']:.2e} "
        f"(tolerance {UNIFY_VOLUME_TOL:g})"
    )
    stats["output_faces"] = simplify_stats["output_faces"]
    stats["output_edges"] = simplify_stats["output_edges"]

    with Timer(rl, "validate"):
        invalid = [i for i, s in enumerate(result_solids) if not occ.is_valid(s)]
        total_volume = sum(occ.volume(s) for s in result_solids)
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

    with Timer(rl, "verify"):
        n_solids = round_trip_check(args.output)
    rl.line(f"round-trip re-import found {n_solids} solid(s)")
    if n_solids != len(result_solids):
        raise OutputError(
            f"Round-trip check disagrees with what this run believes it wrote: "
            f"{n_solids} solids in the file vs {len(result_solids)} expected."
        )

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


def _unify(solids: list[TopoDS_Shape]):
    """Compact each solid's B-rep, verifying the geometry is unchanged.

    Each solid is unified on its own rather than as one compound: it keeps a 1:1
    mapping so the count guard below is exact, and leaves the step
    straightforward to parallelise if it ever becomes the bottleneck at scale.
    """
    faces_before = edges_before = 0
    faces_after = edges_after = 0
    worst_drift = 0.0
    out: list[TopoDS_Shape] = []
    for solid in solids:
        f, e = occ.count_subshapes(solid)
        faces_before += f
        edges_before += e
        pre_volume = occ.volume(solid)

        merged = occ.unify_same_domain(solid)
        produced = occ.solids(merged)
        if len(produced) != 1:
            raise ProcessingError(
                f"Same-domain unification turned one solid into {len(produced)}. "
                f"It must only re-describe the boundary, never re-partition the body."
            )
        post_volume = occ.volume(produced[0])
        drift = abs(post_volume - pre_volume) / max(abs(pre_volume), 1.0)
        worst_drift = max(worst_drift, drift)
        if drift > UNIFY_VOLUME_TOL:
            raise ProcessingError(
                f"Same-domain unification changed a solid's volume from "
                f"{pre_volume:.6f} to {post_volume:.6f} mm^3 ({drift:.2e} relative, "
                f"tolerance {UNIFY_VOLUME_TOL:g}). Faces that are not genuinely the "
                f"same surface were merged, so the boundary moved."
            )

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
    }


def _partition_kept(comps, keep_labels: set[int], boundary_pieces):
    """Split surviving graph vertices back into interior nodes and boundary pieces."""
    keep_interior: list[tuple[int, int, int]] = []
    keep_pieces = []
    for i, v in enumerate(comps.vertices):
        if int(comps.labels[i]) not in keep_labels:
            continue
        if v.is_interior:
            keep_interior.append(v.node)
        else:
            keep_pieces.append(boundary_pieces[v.piece_index])
    nodes_arr = (
        np.array(keep_interior, dtype=np.int64)
        if keep_interior
        else np.empty((0, 3), dtype=np.int64)
    )
    return nodes_arr, keep_pieces


def _log_dropped(rl: RunLog, comps, dropped: set[int], threshold: float) -> None:
    """Report floating-body removals as one aggregate line, never one per body.

    The Julia implementation logged a line per removed solid, which turned a
    pathological run into a multi-hour, multi-thousand-line tail
    (docs/algorithm.md §11.2). One line carries the same information.
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


def _close_solids(sewn: TopoDS_Shape) -> list[TopoDS_Shape]:
    """Turn the sewing result into closed solids, rejecting any open shell."""
    out: list[TopoDS_Shape] = []
    found = occ.shells(sewn)
    if not found:
        raise ProcessingError(
            "Stitching produced no shell at all — every interface failed to close."
        )
    open_shells = 0
    for shell in found:
        if not shell.Closed():
            open_shells += 1
            continue
        out.append(occ.make_solid(shell))
    if open_shells:
        raise ProcessingError(
            f"{open_shells} of {len(found)} stitched shells are not closed, so the "
            f"output would not be watertight. This means two junction interfaces "
            f"that should have matched did not."
        )
    return out
