"""Print support generation (specification.md §4.6, docs/algorithm.md §14).

Loaded only by a run with ``--support``: :mod:`latticegen2.pipeline` imports
this package inside the one branch that needs it, so a run without supports
never executes — or even imports — any of it.

Three stages, between ``simplify`` and ``validate``:

``overhang``
    Find where the input body's inner surface needs support (:mod:`.overhang`).
``support``
    Hand the completed, unified lattice, the body and the overhangs to the
    generator (:mod:`.generate`) through the contract in :mod:`.contract`,
    which builds the support bodies and trims them to the body.
``pockets``
    The generator's pocket search and join, which it opens by calling
    ``SupportInput.pockets_begin``; then check that what comes back is support
    material that overlaps nothing, and append it to the output beside the
    lattice (:mod:`.assemble`).

Supports are **not** merged with the lattice: they ship as separate bodies that
share faces with it. ``validate`` then gates lattice and supports together.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..runlog import Timer
from . import assemble, generate, overhang
from .contract import SupportInput, SupportResult
from .params import STEP, SUPPORT_START_DEG

LATTICE_BACKS_ALL_NOTE = (
    "note: the lattice backs every inner overhang at this orientation; no "
    "support was generated."
)
"""Printed when there were overhangs and the generator returned no supports:
every ceiling point lies in a strut already (docs/algorithm.md §14.3).
``tools/e2e.py`` reads it: coverage is still checked, the body checks have no
body to ask."""

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
                kept_nodes, pool, tmpdir, stats, lost=()) -> SupportOutcome:
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
    # Two stages around one generator call: ``support`` until the generator
    # says its bodies are built, ``pockets`` from then on, and the overlap
    # check at the end of it.
    timer = [Timer(rl, "support").__enter__()]

    def pockets_begin():
        if timer[0].name == "support":
            timer[0].__exit__(None, None, None)
            timer[0] = Timer(rl, "pockets").__enter__()

    try:
        if len(samples) == 0:
            # Not a failure: a body with no inner overhang at this orientation
            # needs no support, and saying so is the whole result.
            result = SupportResult(supports=[])
        else:
            result = generate.generate(SupportInput(
                lattice=tuple(lattice), body=body, body_path=body_path,
                mesh=mesh, mesh_path=mesh_path, lp=lp,
                thickness_mm=args.support_thickness,
                samples=samples, kept_nodes=kept_nodes, lost_lattice=tuple(lost), pool=pool,
                workers=args.workers, tmpdir=tmpdir, report=rl.substage,
                pockets_begin=pockets_begin,
            ))
        pockets_begin()
        # Before the overlap check, which can fail the run: what the generator
        # did is what a failure there has to be read against.
        for line in result.lines:
            rl.line(line)
        for note in result.notes:
            rl.always(note)
        stats.update(result.stats)
        assembled = assemble.assemble(lattice, list(result.supports), lp,
                                      kept_nodes, report=rl.substage, lost=lost)
    finally:
        timer[0].__exit__(None, None, None)

    stats.update(assembled.stats)
    stats["support_solids_generated"] = len(result.supports)
    if len(samples) == 0:
        rl.always(NO_OVERHANG_NOTE)
    elif not result.supports:
        # Only when it is true: with every piece dropped, the generator's own
        # warnings say what is bare, and this note would contradict them.
        if not result.stats.get("support_ceiling_points_uncovered"):
            rl.always(LATTICE_BACKS_ALL_NOTE)
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
