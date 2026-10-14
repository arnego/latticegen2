"""Print support generation (specification.md §4.6, docs/algorithm.md §14).

Loaded only by a run with ``--support``: :mod:`latticegen2.pipeline` imports
this package inside the one branch that needs it, so a run without supports
never executes — or even imports — any of it.

Two stages, between ``simplify`` and ``validate``:

``overhang``
    Find where the input body's inner surface overhangs (:mod:`.overhang`).
``support``
    Hand the completed, unified lattice, the body and the overhangs to the
    generator (:mod:`.generate`) through the contract in :mod:`.contract`;
    then merge the supports it returns with the lattice into the final
    product, unify what changed and check that no lattice was lost
    (:mod:`.merge`).

``validate`` then gates the final product, in one pass.
"""

from __future__ import annotations

from dataclasses import dataclass

from OCP.TopAbs import TopAbs_ShapeEnum

from ..errors import ProcessingError
from ..runlog import Timer
from . import generate, merge, overhang
from .contract import SupportInput, SupportResult

NOT_IMPLEMENTED_WARNING = (
    "WARNING: no support algorithm is implemented yet; the output carries no "
    "supports (docs/algorithm.md §14.3)."
)
"""Printed by every run whose generator returns no supports while there were
overhangs to support. ``tools/e2e.py`` reads it to report the support checks
as skipped rather than failed."""

NO_OVERHANG_NOTE = (
    "note: no inner overhang exceeds the limit at this orientation; no support "
    "was generated."
)
"""Printed when the body has nothing to support. ``tools/e2e.py`` reads it too:
coverage of nothing is not a check that can pass or fail."""


@dataclass
class SupportOutcome:
    """What the support stage hands back to the pipeline."""

    solids: list
    """The final product, to be validated and written."""
    lattice_volume: "float | None"
    """The lattice's volume, measured before supports were merged in; ``None``
    when nothing was merged, the output then being the lattice."""


def run_support(args, rl, lp, body, body_path, mesh, mesh_path, lattice,
                kept_nodes, pool, tmpdir, stats, unify) -> SupportOutcome:
    """Run the two support stages on the completed lattice.

    ``unify(solids) -> (solids, stats)`` is the pipeline's same-domain
    unification, bound to the run's pool: the merged solids are compacted
    exactly as the lattice was.
    """
    pitch = overhang.sample_pitch(lp.t, args.support_thickness)
    with Timer(rl, "overhang"):
        samples = overhang.detect(body, lp.b, args.overhang, pitch,
                                  report=rl.substage)
    kinds = samples.counts()
    rl.line(
        f"overhangs: {len(samples)} sample(s) at {pitch:g} mm pitch over "
        f"{samples.area:.1f} mm^2 of inner surface flatter than "
        f"{args.overhang - overhang.MARGIN_DEG:g} deg from the build direction "
        f"({kinds['edge']} on edges, {kinds['vertex']} on vertices); "
        f"{int(samples.required.sum())} of them exceed the {args.overhang:g} deg "
        f"limit itself and must be proven covered"
    )
    stats["overhang_area_mm2"] = round(samples.area, 2)
    stats["overhang_samples"] = len(samples)
    stats["overhang_samples_required"] = int(samples.required.sum())

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
                overhang_deg=args.overhang, thickness_mm=args.support_thickness,
                samples=samples, kept_nodes=kept_nodes, pool=pool,
                workers=args.workers, tmpdir=tmpdir, report=rl.substage,
            ))
        check_result(result, lattice)
        merged = merge.merge_with_lattice(
            lattice, list(result.supports), unify=unify, report=rl.substage)

    for line in result.lines:
        rl.line(line)
    for note in result.notes:
        rl.always(note)
    stats.update(result.stats)
    stats.update(merged.stats)
    stats["support_solids_generated"] = len(result.supports)
    if len(samples) == 0:
        rl.always(NO_OVERHANG_NOTE)
    elif not result.supports:
        rl.always(NOT_IMPLEMENTED_WARNING)
    else:
        rl.line(
            f"supports: {len(result.supports)} solid(s), "
            f"{merged.support_volume:.4f} mm^3 as generated, merged into "
            f"{merged.stats['support_merge_groups']} group(s); output "
            f"{merged.output_volume:.4f} mm^3 against "
            f"{merged.lattice_volume:.4f} mm^3 of lattice"
        )
    return SupportOutcome(merged.solids, merged.lattice_volume)


def check_result(result: SupportResult, lattice: list) -> None:
    """Fail (exit 4) if ``result`` breaks :class:`SupportResult`'s contract.

    Only what the contract itself promises is checked here: support solids,
    none of them the lattice itself. Whether the merge kept every bit of
    lattice is :func:`.merge.check_volumes`'s question, and whether the
    product is sound is ``validate``'s.
    """
    supports = list(result.supports)
    for i, s in enumerate(supports):
        if s is None or s.IsNull() or s.ShapeType() != TopAbs_ShapeEnum.TopAbs_SOLID:
            raise ProcessingError(
                f"The support generator returned something other than a solid "
                f"at position {i} of {len(supports)} (docs/algorithm.md §14.2)."
            )
        if any(s.IsSame(l) for l in lattice):
            raise ProcessingError(
                f"The support generator returned a lattice solid as support "
                f"at position {i} of {len(supports)}; it must return support "
                f"material only (docs/algorithm.md §14.2)."
            )
