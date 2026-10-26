"""Print support generation (specification.md §4.6, docs/algorithm.md §14).

Loaded only by a run with ``--support``: :mod:`latticegen2.pipeline` imports
this package inside the one branch that needs it, so a run without supports
never executes — or even imports — any of it.

Two stages, between ``simplify`` and ``validate``:

``overhang``
    Find where the input body's inner surface needs support (:mod:`.overhang`).
``support``
    Hand the completed, unified lattice, the body and the overhangs to the
    generator (:mod:`.generate`) through the contract in :mod:`.contract`;
    check that what comes back is support material that overlaps nothing, and
    append it to the output beside the lattice (:mod:`.assemble`).

Supports are **not** merged with the lattice: they ship as separate bodies that
share faces with it. ``validate`` then gates lattice and supports together.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..runlog import Timer
from . import assemble, generate, overhang
from .contract import SupportInput, SupportResult
from .params import STEP, SUPPORT_START_DEG

NOT_IMPLEMENTED_WARNING = (
    "WARNING: the support generator returned no supports although the input "
    "has inner overhangs; the output carries no supports (docs/algorithm.md "
    "§14.3)."
)
"""Printed by every run whose generator returns no supports while there were
overhangs to support. ``tools/e2e.py`` reads it to report the support checks
as skipped rather than failed."""

NO_OVERHANG_NOTE = (
    "note: no inner overhang exceeds 60 deg from vertical at this orientation; "
    "no support was generated."
)
"""Printed when the body has nothing to support. ``tools/e2e.py`` reads it too:
coverage of nothing is not a check that can pass or fail."""


@dataclass
class SupportOutcome:
    """What the support stage hands back to the pipeline."""

    solids: list
    """The final product, to be validated and written: the lattice solids as
    they came in, followed by the support bodies."""
    lattice_volume: "float | None"
    """The lattice's volume; ``None`` when no support was generated, the output
    then being the lattice."""


def run_support(args, rl, lp, body, body_path, mesh, mesh_path, lattice,
                kept_nodes, pool, tmpdir, stats) -> SupportOutcome:
    """Run the two support stages on the completed lattice."""
    with Timer(rl, "overhang"):
        samples = overhang.detect(body, lp.b, STEP, report=rl.substage)
    kinds = samples.counts()
    rl.line(
        f"overhangs: {len(samples)} sample(s) at {STEP:g} mm over "
        f"{samples.area:.1f} mm^2 of inner surface overhanging more than "
        f"{SUPPORT_START_DEG:g} deg from vertical "
        f"({kinds['edge']} on edges, {kinds['vertex']} on vertices)"
    )
    stats["overhang_area_mm2"] = round(samples.area, 2)
    stats["overhang_samples"] = len(samples)

    lattice = list(lattice)
    with Timer(rl, "support"):
        if len(samples) == 0:
            # Not a failure: a body with no inner overhang at this orientation
            # needs no support, and saying so is the whole result.
            result = SupportResult(supports=[])
        else:
            result = generate.generate(SupportInput(
                lattice=tuple(lattice), body=body, body_path=body_path,
                mesh=mesh, mesh_path=mesh_path, lp=lp,
                thickness_mm=args.support_thickness,
                samples=samples, kept_nodes=kept_nodes, pool=pool,
                workers=args.workers, tmpdir=tmpdir, report=rl.substage,
            ))
        assembled = assemble.assemble(lattice, list(result.supports), lp,
                                      kept_nodes, report=rl.substage)

    for line in result.lines:
        rl.line(line)
    for note in result.notes:
        rl.always(note)
    stats.update(result.stats)
    stats.update(assembled.stats)
    stats["support_solids_generated"] = len(result.supports)
    if len(samples) == 0:
        rl.always(NO_OVERHANG_NOTE)
    elif not result.supports:
        rl.always(NOT_IMPLEMENTED_WARNING)
    else:
        rl.line(
            f"supports: {len(result.supports)} separate bod(y/ies), "
            f"{assembled.support_volume:.4f} mm^3, beside "
            f"{assembled.lattice_volume:.4f} mm^3 of lattice"
        )
    return SupportOutcome(assembled.solids, assembled.lattice_volume)


# Kept importable from the package, where the stage's tests and the generator's
# own checks look for it.
check_result = assemble.check_result
