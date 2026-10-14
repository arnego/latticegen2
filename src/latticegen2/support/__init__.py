"""Print support generation (specification.md §4.6, docs/algorithm.md §14).

Loaded only by a run with ``--support``: :mod:`latticegen2.pipeline` imports
this package inside the one branch that needs it, so a run without supports
never executes — or even imports — any of it.

Two stages, between ``simplify`` and ``validate``:

``overhang``
    Find where the input body's inner surface overhangs (:mod:`.overhang`).
``support``
    Hand the completed lattice, the body and the overhangs to the generator
    (:mod:`.generate`) through the contract in :mod:`.contract`, and check
    that what comes back keeps it.

``validate`` then gates lattice and supports together, in one pass.
"""

from __future__ import annotations

from OCP.TopAbs import TopAbs_ShapeEnum

from ..errors import ProcessingError
from ..runlog import Timer
from . import generate, overhang
from .contract import SupportInput, SupportResult

NOT_IMPLEMENTED_WARNING = (
    "WARNING: no support algorithm is implemented yet; the output carries no "
    "supports (docs/algorithm.md §14.3)."
)
"""Printed by every run whose generator reports ``generated=False`` while
there were overhangs to support. ``tools/e2e.py`` reads it to report the
support checks as skipped rather than failed."""


def run_support(args, rl, lp, body, body_path, mesh, mesh_path, lattice,
                kept_nodes, pool, tmpdir, stats) -> list:
    """Run the two support stages on the completed lattice.

    Returns the solids the rest of the pipeline should validate and write.
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
            result = SupportResult(solids=lattice, generated=False)
        else:
            result = generate.generate(SupportInput(
                lattice=tuple(lattice), body=body, body_path=body_path,
                mesh=mesh, mesh_path=mesh_path, lp=lp,
                overhang_deg=args.overhang, thickness_mm=args.support_thickness,
                samples=samples, kept_nodes=kept_nodes, pool=pool,
                workers=args.workers, tmpdir=tmpdir, report=rl.substage,
            ))
        check_result(result, lattice)

    for line in result.lines:
        rl.line(line)
    for note in result.notes:
        rl.always(note)
    stats.update(result.stats)
    stats["support_generated"] = result.generated
    if len(samples) == 0:
        rl.always("note: no inner overhang exceeds the limit at this "
                  "orientation; no support was generated.")
    elif not result.generated:
        rl.always(NOT_IMPLEMENTED_WARNING)
    return list(result.solids)


def check_result(result: SupportResult, lattice: list) -> None:
    """Fail (exit 4) if ``result`` breaks :class:`SupportResult`'s contract.

    Only what the contract itself promises is checked here; whether the
    solids are sound is ``validate``'s question, asked next.
    """
    solids = list(result.solids)
    if not solids:
        raise ProcessingError(
            "The support generator returned no solids; it must return the "
            "complete output, lattice included (docs/algorithm.md §14.2)."
        )
    for i, s in enumerate(solids):
        if s is None or s.IsNull() or s.ShapeType() != TopAbs_ShapeEnum.TopAbs_SOLID:
            raise ProcessingError(
                f"The support generator returned something other than a solid "
                f"at position {i} of {len(solids)} (docs/algorithm.md §14.2)."
            )
    if not result.generated and (
        len(solids) != len(lattice)
        or not all(a.IsSame(b) for a, b in zip(solids, lattice))
    ):
        raise ProcessingError(
            "The support generator reported that it generated nothing but "
            "returned solids other than the lattice it was given "
            "(docs/algorithm.md §14.2)."
        )
