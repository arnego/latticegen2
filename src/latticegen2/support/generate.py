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
* **Coverage.** Every overhang sample within 0.5 mm of support or lattice
  (:data:`.params.REACH`) — ``tools/verify_geometry.overhang_coverage`` asks
  this of the written file in e2e.
* **Self-supporting.** No support surface or edge overhangs more than the rod
  angle :data:`.params.SUPPORT_OVERHANG_DEG` from ``inp.lp.b`` —
  ``verify_geometry.self_supporting``.
* **Thickness.** No free-standing member thinner than ``inp.thickness_mm``.
* **Separate bodies.** One body per cell, sharing faces with the lattice and
  the input's surfaces and overlapping nothing. Support bodies are exempt from
  the floating-body rule (specification.md §5).
* **Valid, writable solids.** Closed and ``BRepCheck_Analyzer``-valid, and
  surviving a STEP round trip, which ``validate`` checks.
* **Supports only.** Never a copy of the lattice.
"""

from __future__ import annotations

import numpy as np

from ..errors import ProcessingError
from .contract import SupportInput, SupportResult
from .params import REACH, STEP

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
    from . import build, cells, obstacles, rays, regions
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
        exempt = np.zeros(len(rmesh.verts), dtype=bool)
        exempt[reg.verts] = got.length <= REACH
        cache = {}
        build.prepare_fits(part.cells, cache)
        for c in part.cells:
            n_cells += 1
            if report is not None:
                report("building support cells", n_cells, None)
            live = c.verts[~exempt[c.verts]]
            inside = rmesh.verts[live] - 2.0 * rays.APEX_NUDGE * rmesh.normals[live]
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
    solids, resolved = _resolve_overlaps(solids)
    stats["support_overlaps_resolved"] = resolved
    return solids, stats, lines


OVERLAP_VOLUME_REL_TOL = 1e-5
"""Relative slack, on the body cut, between the volume a cut removes and the
overlap measured before it."""


def _resolve_overlaps(solids: list):
    """S6: ``(solids, resolved)``, no two bodies sharing volume.

    Neighbouring cells meet on a shared run and diverge below it, so they
    normally share faces only. Two things leave them sharing some volume:
    first-contact rays that cross, and a cell's own edge -- along the exempt
    band, fitted through the mesh's staircase -- swinging across the shared
    run near the triple point where both end. Measured: slivers of 1e-5 to
    3e-4 mm^3, one to seven microns thick, on every wall cell of the test
    chamber and on most cells of Hole-block-chamfered, and 0.11 mm^3 twice
    there where rays cross.

    The later body of each pair is cut by the earlier, so the shared volume
    stays in exactly one of them: their union, and so the coverage, is
    unchanged, and the new face lies on the other body's side face, which is
    self-supporting. Fusing the two instead was tried first and fails on this
    input -- chains of 8 and 16 cells that each overlap the next fused into
    nothing valid. The cut is checked by volume and validity; a pair it fails
    on is fused, and a pair neither resolves fails the run
    (docs/algorithm.md §14.3 S6).
    """
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse

    from .. import occ
    from .assemble import SUPPORT_OVERLAP_TOL, _boxes, _common_volume, _touch, _where
    from .build import _boolean

    solids = list(solids)
    boxes = _boxes(solids)
    resolved = 0
    i = 0
    while i < len(solids):
        j = i + 1
        while j < len(solids):
            if solids[i] is None or solids[j] is None or not _touch(boxes[i], boxes[j]):
                j += 1
                continue
            overlap = _common_volume(solids[i], [solids[j]])
            if overlap <= SUPPORT_OVERLAP_TOL:
                j += 1
                continue
            vj = occ.volume(solids[j])
            slack = max(OVERLAP_VOLUME_REL_TOL * vj, SUPPORT_OVERLAP_TOL)
            got = _boolean(BRepAlgoAPI_Cut, solids[j], [solids[i]])
            pieces = sorted((s for s in (occ.solids(got) if got is not None else [])
                             if occ.volume(s) > SUPPORT_OVERLAP_TOL),
                            key=occ.volume, reverse=True)
            kept = sum(occ.volume(s) for s in pieces)
            if (pieces and all(occ.is_valid(s, parallel=False) for s in pieces)
                    and abs((vj - kept) - overlap) <= slack):
                solids[j] = pieces[0]
                boxes[j] = _boxes([pieces[0]])[0]
                for extra in pieces[1:]:
                    solids.append(extra)
                    boxes.append(_boxes([extra])[0])
                resolved += 1
                j += 1
                continue
            vi = occ.volume(solids[i])
            got = _boolean(BRepAlgoAPI_Fuse, solids[i], [solids[j]])
            fused = occ.solids(got) if got is not None else []
            total = vi + vj
            if (len(fused) != 1 or not occ.is_valid(fused[0], parallel=False)
                    or abs(occ.volume(fused[0]) - (total - overlap))
                    > max(OVERLAP_VOLUME_REL_TOL * total, SUPPORT_OVERLAP_TOL)):
                raise ProcessingError(
                    f"Two support bodies near {_where(boxes[i])} share {overlap:.3e} "
                    f"mm^3, and neither cutting one by the other nor fusing them "
                    f"gave valid bodies (docs/algorithm.md §14.3 S6)."
                )
            solids[i] = fused[0]
            boxes[i] = _boxes([fused[0]])[0]
            solids[j] = None
            resolved += 1
            # The fused body is new: ask every later body about it again.
            j = i + 1
        i += 1
    return [s for s in solids if s is not None], resolved


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
