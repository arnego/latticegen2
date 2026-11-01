"""The support generator (docs/algorithm.md §14.3).

Ceiling points (:mod:`.regions`), frustum supports on the lattice grid
(:mod:`.frustum`), the pockets they seal (:mod:`.voids`) and the joining of
bodies that share a face (:mod:`.join`), in that order.

:func:`generate` receives the completed lattice and the input body
(:class:`.contract.SupportInput`) and returns the support solids alone
(:class:`.contract.SupportResult`).

What happens to them is not the generator's job: the support stage checks
they overlap neither the lattice nor each other and appends them beside the
lattice as separate bodies (:mod:`.assemble`); then ``validate`` ->
export-truth -> ``export`` run on the final product. What an implementation
owes that chain:

* **Inside the body.** No material outside the input solid (specification.md
  §1). A boolean against the body can be wrong while reporting success
  (docs/algorithm.md §7.2): every trim is held to a sampled estimate and to
  ``classify.OutsideProbe`` over ``inp.mesh``.
* **Coverage.** Every ceiling point inside support or lattice material --
  ``tools/verify_geometry.overhang_coverage`` asks this of the written file in
  e2e.
* **Self-supporting.** No support surface or edge overhangs more than the rod
  angle :data:`.params.SUPPORT_OVERHANG_DEG` from ``inp.lp.b`` -- by
  construction, every face lying on a plane of the lattice frame.
* **Separate, manifold bodies.** Separate from the lattice, sharing faces with
  it and with the input's surfaces and overlapping nothing; supports sharing a
  face or volume with each other are one body (:mod:`.join`). No edge of a
  body is used by more than two of its faces. Support bodies are exempt from
  the floating-body rule (specification.md §5).
* **Valid, writable solids.** Closed and ``BRepCheck_Analyzer``-valid, and
  surviving a STEP round trip, which ``validate`` checks.
* **Supports only.** Never a copy of the lattice.
"""

from __future__ import annotations

import os
import tempfile
import time

import numpy as np

from .contract import SupportInput, SupportResult


def supports_for(body, mesh, lattice: list, kept_nodes, lp, thickness: float,
                 report=None, tmpdir: str | None = None, pool=None, staged=None,
                 pockets_begin=None, samples=None, lost=()):
    """``(solids, stats, lines)``: every support body for ``body``.

    The generator proper, shared by :func:`generate` and by
    ``tools/support_dev.py``, which calls it with no lattice: every anchor is
    then on a wall. ``pool`` and ``staged`` -- the body's ``.brep`` and the
    mesh's ``.npz`` -- let the pocket search run across the run's workers;
    without them it runs in this process. ``samples`` is what the detector
    found, detected here if not given; ``lost`` the junction pieces the
    floating-body rule took from kept nodes.
    """
    from .. import occ
    from . import build, frustum, join, overhang, regions, voids
    from .assemble import _LocalLattice

    t0 = time.perf_counter()
    local = _LocalLattice(lp, np.asarray(kept_nodes), lost)
    if samples is None:
        samples = overhang.detect(body, lp.b)
    # Ceiling a strut ends on is backed by the lattice already.
    rmesh = regions.build(body, lp.b, report=report,
                          covered=local.contains if len(kept_nodes) else None)
    gates = build.Gates(mesh, lp.b)
    ceiling_s = time.perf_counter() - t0
    if tmpdir is None:
        tmpdir = tempfile.mkdtemp(prefix="latticegen2-support-")
    # Keels and tips have no face to refine: their samples are ceiling points
    # as they stand, asked what backs them straight below (§14.1).
    other = samples.kind != overhang.FACE
    got = frustum.build(body, mesh, rmesh, local, lp, thickness, gates, tmpdir,
                        report=report, extra=(samples.points[other], -samples.inward[other]))
    stats = dict(got.stats)
    stats["support_ceiling_s"] = round(ceiling_s, 2)
    lines = list(got.lines) + list(got.notes)
    if not got.solids:
        stats["support_bodies"] = 0
        return [], stats, lines

    # S7 per group of touching bodies: free space the supports seal off is
    # filled, and a fill joins the bodies whose faces it lies on.
    if pockets_begin is not None:
        pockets_begin()
    t3 = time.perf_counter()
    out, fills, unsearched = [], [], []
    gs = join.groups(got.solids)
    for k, g in enumerate(gs):
        if report is not None:
            report("looking for sealed pockets", k, len(gs))
        mine = [got.solids[i] for i in g]
        pockets = voids.find_pockets(body, mesh, mine, local.near, thickness, lp.b,
                                     pool=pool, staged=staged, report=report,
                                     reach=local.reach)
        fills += pockets.fills
        unsearched += pockets.unsearched
        if pockets.fills:
            mine, _js = join.join(mine + pockets.fills, overlaps=False)
        out += mine
    for s in out:
        nm = frustum.non_manifold_edges(s)
        if nm:
            from ..errors import ProcessingError

            where = np.round(np.asarray(occ.centroid(s), dtype=float), 3).tolist()
            raise ProcessingError(
                f"A support body near {where} has {nm} edge(s) used by more than two "
                f"faces: it is not a manifold solid (docs/algorithm.md §14.3)."
            )
    stats["support_pockets_s"] = round(time.perf_counter() - t3, 2)
    stats["support_pockets"] = len(fills)
    stats["support_pocket_volume_mm3"] = round(sum(occ.volume(f) for f in fills), 4)
    stats["support_bodies"] = len(out)
    stats["support_pocket_slabs_unsearched"] = len(unsearched)
    if unsearched:
        lines.append(f"{voids.UNSEARCHED_WARNING} {len(unsearched)} slab(s) of the "
                     f"support envelope at any fuzzy value; free space touching them "
                     f"is taken as open and no pocket there is filled:")
        lines += unsearched
    lines.append(
        f"support: {len(out)} bod(y/ies), {sum(occ.volume(s) for s in out):.3f} mm^3, "
        f"{len(fills)} sealed pocket(s) filled; "
        f"{int(got.covered.sum())} of {len(got.covered)} ceiling point(s) backed "
        f"by the support grid or a strut"
    )
    lines.append(
        "support timing: " + ", ".join(
            f"{k[len('support_'):-2]} {stats[k]:.1f} s" for k in
            ("support_ceiling_s", "support_grid_s", "support_staircase_s",
             "support_trim_s", "support_backed_s", "support_pockets_s"))
    )
    return out, stats, lines


def generate(inp: SupportInput) -> SupportResult:
    """Support solids backing the ceiling of ``inp.body``, anchored on ``inp.lattice``."""
    from .frustum import DROPPED_WARNING
    from .voids import UNSEARCHED_WARNING

    staged = {"body_path": inp.body_path, "mesh_path": inp.mesh_path,
              "lp_args": (inp.lp.cc, inp.lp.t, tuple(inp.lp.orient)),
              "kept_nodes": np.asarray(inp.kept_nodes), "tmpdir": inp.tmpdir}
    if inp.lost_lattice:
        from .. import occ
        from ..parallel import write_brep

        staged["lost_path"] = os.path.join(inp.tmpdir, "lost-lattice.brep")
        write_brep(occ.compound(list(inp.lost_lattice)), staged["lost_path"])
    solids, stats, lines = supports_for(inp.body, inp.mesh, list(inp.lattice),
                                        inp.kept_nodes, inp.lp, inp.thickness_mm,
                                        report=inp.report, tmpdir=inp.tmpdir,
                                        pool=inp.pool, staged=staged,
                                        pockets_begin=inp.pockets_begin,
                                        samples=inp.samples, lost=inp.lost_lattice)
    # The dropped pieces' and unsearched slabs' reports go to the console
    # whatever -v says.
    notes = []
    for head, count in ((DROPPED_WARNING, stats.get("support_pieces_dropped", 0)),
                        (UNSEARCHED_WARNING, stats.get("support_pocket_slabs_unsearched", 0))):
        k = next((i for i, line in enumerate(lines) if line.startswith(head)), None)
        if k is not None:
            notes += lines[k:k + 1 + count]
            lines = lines[:k] + lines[k + 1 + count:]
    return SupportResult(supports=solids, stats=stats, lines=lines, notes=notes)
