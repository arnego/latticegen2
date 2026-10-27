"""The support generator (docs/algorithm.md §14.3).

Regions (:mod:`.regions`), first-contact rays (:mod:`.rays`), cells and
parting curves (:mod:`.cells`) and cell solids (:mod:`.build`), in that order.

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
  (docs/algorithm.md §7.2): check a trim with something that is not the same
  kernel call, such as ``classify.OutsideProbe`` over ``inp.mesh``.
* **Coverage.** Every overhang sample inside support or lattice material,
  up to the overhang's own edge — ``tools/verify_geometry.overhang_coverage``
  asks this of the written file in e2e.
* **Self-supporting.** No support surface or edge overhangs more than the rod
  angle :data:`.params.SUPPORT_OVERHANG_DEG` from ``inp.lp.b`` —
  ``verify_geometry.self_supporting``.
* **Thickness.** No free-standing member thinner than ``inp.thickness_mm``.
* **Separate bodies.** Separate from the lattice, sharing faces with it and
  with the input's surfaces and overlapping nothing; supports sharing a face
  or volume with each other are one body (:mod:`.join`). Support bodies are exempt from
  the floating-body rule (specification.md §5).
* **Valid, writable solids.** Closed and ``BRepCheck_Analyzer``-valid, and
  surviving a STEP round trip, which ``validate`` checks.
* **Supports only.** Never a copy of the lattice.
"""

from __future__ import annotations

import numpy as np

from ..errors import ProcessingError
from .contract import SupportInput, SupportResult
from .params import STEP

OBSTACLE_DEFLECTION = 0.02
"""mm. Tessellation the rays land on. First contacts snap to its facets, so it
bounds how far a target sits from where the true surface would put it."""

OVERSHOOT_RETRIES = 2
"""Times a failing cell is rebuilt with twice the overshoot before the run fails."""


def overshoot_for(thickness: float, deviation: float) -> float:
    """``delta``, mm: how far a cell's loft reaches past its rays' targets."""
    return max(thickness / 2.0, 2.0 * deviation, 0.2)


def supports_for(body, mesh, lattice: list, kept_nodes, lp, thickness: float,
                 report=None):
    """``(solids, stats, lines)``: every support body for ``body``.

    The generator proper, shared by :func:`generate` and by
    ``tools/support_dev.py``, which calls it with no lattice.
    """
    from .. import occ
    from . import build, cells, join, lean, obstacles, rays, regions, voids
    from .assemble import _LocalLattice

    stats = {}
    lines = []
    rmesh = regions.build(body, lp.b, report=report)
    regs = regions.regions(rmesh)
    stats["support_regions"] = len(regs)
    if not regs:
        return [], stats, lines
    obs = obstacles.from_shapes(body, list(lattice), deflection=OBSTACLE_DEFLECTION)
    local = _LocalLattice(lp, np.asarray(kept_nodes))
    gates = build.Gates(mesh, lp.b)
    inside_body = gates.inside_body
    delta = overshoot_for(thickness, mesh.deviation)
    solids = []
    n_cells = 0
    retried = 0
    for k, reg in enumerate(regs):
        apex = rmesh.verts[reg.verts] - rays.APEX_NUDGE * rmesh.normals[reg.verts]
        got = rays.first_contact(apex, lp.b, obs, report=report)
        part = cells.partition(rmesh, reg, got, obs, lp.b, thickness, report=report,
                               inside_body=gates.inside_body, overshoot=delta)
        leaned = lean.lean_cells(part.cells, obs, lp.b, thickness, delta)
        part.stats.update(leaned)
        cache = {}
        build.prepare_fits(part.cells, cache)
        for c in part.cells:
            n_cells += 1
            if report is not None:
                report("building support cells", n_cells, None)
            inside = rmesh.verts[c.verts] - 2.0 * rays.APEX_NUDGE * rmesh.normals[c.verts]
            d = delta
            whys = []
            for attempt in range(OVERSHOOT_RETRIES + 1):
                lo = np.min([lp_.targets.min(axis=0) for lp_ in c.loops], axis=0) - d
                hi = np.max([lp_.points.max(axis=0) for lp_ in c.loops], axis=0) + d
                built = build.build_cell(c, body, local.near((lo, hi)), lp.b, d, cache,
                                         inside_body)
                if built.solid is not None:
                    fails, _gs = gates.check(built.solid, inside)
                    if not fails:
                        break
                    why = "; ".join(fails)
                else:
                    why = built.failure
                whys.append(f"overshoot {d:g} mm: {why}")
                d *= 2.0
                retried += 1
            else:
                where = np.round(c.loops[0].points.mean(axis=0), 3).tolist()
                raise ProcessingError(
                    f"The support cell near {where} could not be built ({' | '.join(whys)}) "
                    f"(docs/algorithm.md §14.3 S4)."
                )
            solids.append(built.solid)
        lines.append(f"support region {k}: {reg.area:.1f} mm^2, {len(part.cells)} "
                     f"cell(s), {part.stats}")
    stats["support_cells"] = n_cells
    stats["support_span_fit_mm"] = round(build.worst_span_fit(reset=True), 4)
    stats["support_cells_retried"] = retried
    # S6 before S5: pockets are found by cutting the supports out of the free
    # space, which only removes their volume exactly once they are disjoint.
    if report is not None:
        report("joining supports", 0, None)
    bodies, jstats = join.join(solids)
    stats.update(jstats)
    # S5 per group of touching bodies: a region's cells always touch, and a
    # body joined across two regions takes both into one envelope.
    out = []
    fills = []
    gs = join.groups(bodies)
    for k, g in enumerate(gs):
        if report is not None:
            report("looking for sealed pockets", k, len(gs))
        mine = [bodies[i] for i in g]
        pockets = voids.find_pockets(body, mesh, mine, local.near, thickness, lp.b)
        fills += pockets.fills
        if pockets.fills:
            # A fill's faces lie on the supports that seal it: one body with them.
            mine, fstats = join.join(mine + pockets.fills, overlaps=False)
            stats["support_joins_declined"] += fstats["support_joins_declined"]
        out += mine
    stats["support_pockets"] = len(fills)
    stats["support_pocket_volume_mm3"] = round(sum(occ.volume(f) for f in fills), 4)
    stats["support_bodies"] = len(out)
    return out, stats, lines


def generate(inp: SupportInput) -> SupportResult:
    """Support solids backing ``inp.samples``, anchored on ``inp.lattice``."""
    if len(inp.lattice):
        # Lattice as anchor and cutter is Phase 6 (specification.md §10). Until
        # then the product -- which always has a lattice -- ships none, and the
        # stage says so with NOT_IMPLEMENTED_WARNING; tools/support_dev.py
        # drives the generator without one.
        return SupportResult(supports=[], stats={}, lines=[])
    solids, stats, lines = supports_for(inp.body, inp.mesh, list(inp.lattice),
                                        inp.kept_nodes, inp.lp, inp.thickness_mm,
                                        report=inp.report)
    return SupportResult(supports=solids, stats=stats, lines=lines)
