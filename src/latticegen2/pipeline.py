"""Pipeline orchestration.

    parse -> template -> import -> tessellate -> classify
          -> trim boundary -> connect -> instance interior -> weld
          -> validate -> write STEP -> round-trip -> summary

The shape of this is deliberately flat. In the fuse-free architecture there is
no tiling stage, no hierarchical assembly, no distributed merge rounds and no
fuse-based cleanup, because there is nothing left for them to do: the interior is
instanced rather than fused, connectivity is a graph property rather than a
boolean experiment, and assembly is an index lookup rather than a geometric
search.
"""

from __future__ import annotations

import datetime as _dt
import os
import shutil
from dataclasses import dataclass

import numpy as np
from OCP.BRepTools import BRepTools
from OCP.Standard import Standard_Failure
from OCP.TopoDS import TopoDS_Shape

from . import occ, weld
from .boundary import CAP_AREA_REL_TOL, finalize_pieces, resolve_interfaces, trim_boundary
from .classify import Klass, classify_nodes, tessellate_surface
from .cli import Args
from .connect import build_components
from .errors import InputGeometryError, OutputError, ProcessingError
from .interior import build_interior_shell, extract_template_mesh
from .junction import build_template
from .lattice import candidate_nodes, lattice_params, neighbor_step, part_name
from .lattice import node as lattice_node
from .runlog import RunLog, Timer, format_bytes
from .stepout import generation_params_text, rewrite_step_header, round_trip_check

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
    _log_interfaces(rl, lp, iface)
    stats["interfaces"] = iface.n_pairs
    stats["caps_unpaired"] = len(iface.unpaired)
    stats["caps_mismatched"] = len(iface.mismatched)
    stats["caps_unweldable"] = len(iface.unweldable)
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
        boundary_faces = weld.sew_boundary(kept.pieces, kept.piece_groups)
        rings = weld.interface_rings(lp, tmesh, boundary_faces, want_rings)
    rl.line(
        f"boundary layer: {len(kept.pieces)} piece(s) sewn into "
        f"{sum(len(f) for f in boundary_faces.values())} faces across "
        f"{len(boundary_faces)} component(s); {len(rings)} interface ring(s) located"
    )
    stats["interior_interfaces"] = len(rings)

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
        result_solids, simplify_stats = _unify(result_solids)
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
    stats["output_faces"] = simplify_stats["output_faces"]
    stats["output_edges"] = simplify_stats["output_edges"]
    stats["unmerged_solids"] = simplify_stats["unmerged_solids"]

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


def _unify(solids: list[TopoDS_Shape]):
    """Compact each solid's B-rep, verifying the geometry is unchanged.

    Each solid is unified on its own rather than as one compound: it keeps a 1:1
    mapping so the count guard below is exact, and leaves the step
    straightforward to parallelise if it ever becomes the bottleneck at scale.
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
        "unmerged_solids": skipped,
    }


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

