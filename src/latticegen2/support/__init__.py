"""Print support generation (specification.md §4.6, docs/algorithm.md §14).

Loaded only by a run with ``--support``: :mod:`latticegen2.pipeline` imports
this package inside the one branch that needs it, so a run without supports
never executes — or even imports — any of it.

Three stages, between ``boundary`` and ``connect``:

``overhang``
    Find where the input body's inner surface overhangs
    (:mod:`.overhang`).
``plan``
    Decide where support apexes go, against the lattice that actually survived
    the boundary trim (:mod:`.plan`).
``support``
    Build every lattice cell a support reaches as B-rep and hand it to the
    rest of the pipeline in place of the lattice-only piece (:mod:`.build`).
"""

from __future__ import annotations

import numpy as np

from ..boundary import resolve_interfaces
from ..runlog import Timer, format_bytes
from . import build, overhang, plan


def run_support(args, rl, lp, mesh, body, body_path, mesh_path, tmpdir, pool,
                interior_nodes, boundary, outside_margin, stats):
    """Run the three support stages.

    Returns the piece list and the interior node array the rest of the pipeline
    should continue with: supported cells replace their lattice-only pieces, and
    a supported node that was interior leaves the instanced set, since its
    junction now comes out of a boolean like any boundary junction's.
    """
    params = plan.SupportParams(
        overhang=args.overhang,
        thickness=args.support_thickness,
        pitch=overhang.sample_pitch(lp.t, args.support_thickness),
        deviation=mesh.deviation,
    )

    with Timer(rl, "overhang"):
        samples = overhang.detect(body, lp.b, params.overhang, params.pitch,
                                  report=rl.substage)
    kinds = samples.counts()
    rl.line(
        f"overhangs: {len(samples)} sample(s) at {params.pitch:g} mm pitch over "
        f"{samples.area:.1f} mm^2 of inner surface flatter than "
        f"{params.overhang - overhang.MARGIN_DEG:g} deg from the build direction "
        f"({kinds['edge']} on edges, {kinds['vertex']} on vertices)"
    )
    stats["overhang_area_mm2"] = round(samples.area, 2)
    stats["overhang_samples"] = len(samples)
    if len(samples) == 0:
        # Not a failure: a body with no inner overhang at this orientation
        # needs no support, and saying so is the whole result.
        rl.always("note: no inner overhang exceeds the limit at this "
                  "orientation; no support was generated.")
        for name in ("plan", "support"):
            with Timer(rl, name):
                pass
        return boundary.pieces, interior_nodes

    with Timer(rl, "plan"):
        rl.substage("finding anchors", 0, None)
        iface = resolve_interfaces(lp, interior_nodes, boundary.pieces)
        anchors = plan.find_anchors(
            lp, mesh, params, interior_nodes,
            [p.node for p in boundary.pieces], iface.interfaces,
            near=plan.near_samples(lp, samples.points),
            pool=pool, mesh_path=mesh_path,
        )
        rl.substage("placing supports", 0, None)
        placed = plan.plan_supports(lp, mesh, params, samples.points, anchors)
    rl.line(
        f"support plan: {len(anchors.nodes)} anchor node(s) of "
        f"{anchors.n_lattice} lattice node(s) ({anchors.n_grounded} grounded); "
        f"{placed.stats.get('support_capitals_on_struts', 0)} capital(s) on "
        f"struts, {placed.stats.get('support_capitals_on_nodes', 0)} on nodes, "
        f"{placed.stats.get('support_king_posts', 0)} king post(s), "
        f"{placed.stats.get('support_drop_columns', 0)} drop column(s)"
    )
    stats.update(placed.stats)

    interior_set = {(int(r[0]), int(r[1]), int(r[2])) for r in interior_nodes}
    lattice_nodes = interior_set | {tuple(p.node) for p in boundary.pieces}
    with Timer(rl, "support"):
        pieces, supported, bstats = build.build_supports(
            lp, placed.primitives, samples, mesh, lattice_nodes, interior_set,
            body_path, mesh_path, outside_margin, tmpdir, pool, args.workers,
            report=rl.substage,
        )
    rss = bstats.pop("peak_support_worker_rss", 0)
    if rss:
        rl.note_worker_rss(rss)
        stats["peak_support_worker_memory"] = format_bytes(rss)
    stats.update(bstats)
    rl.line(
        f"supports built: {bstats.get('support_cells', 0)} lattice cell(s) "
        f"({bstats.get('support_cells_interior', 0)} of them interior) -> "
        f"{len(pieces)} piece(s); coverage verified at "
        f"{bstats.get('support_samples_verified', 0)} sample(s)"
        + (f", {bstats['support_cells_retried']} cell(s) needed the alternative "
           f"operand order" if bstats.get("support_cells_retried") else "")
    )

    kept = [p for p in boundary.pieces if tuple(p.node) not in supported]
    keep_interior = np.array(
        [r for r in interior_nodes
         if (int(r[0]), int(r[1]), int(r[2])) not in supported],
        dtype=np.int64,
    ).reshape(-1, 3)
    return kept + pieces, keep_interior
