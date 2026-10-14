"""The support generator (docs/algorithm.md §14.3). **Not yet implemented.**

:func:`generate` receives the completed lattice and the input body
(:class:`.contract.SupportInput`) and returns every solid to write
(:class:`.contract.SupportResult`). Until the algorithm exists it returns the
lattice unchanged with ``generated=False``, and the run warns that its output
carries no supports.

What an implementation owes the rest of the pipeline, which runs
``validate`` -> export-truth -> ``export`` on whatever it returns:

* **Inside the body.** No material outside the input solid (specification.md
  §1). A boolean against the body can be wrong while reporting success
  (docs/algorithm.md §7.2): check a trim with something that is not the same
  kernel call, such as ``classify.OutsideProbe`` over ``inp.mesh``.
* **Coverage.** Every *required* overhang sample (``inp.samples.required``)
  backed by material — ``tools/verify_geometry.overhang_coverage`` asks this
  of the written file in e2e.
* **Self-supporting.** No support surface or edge exceeds
  ``inp.overhang_deg`` from ``inp.lp.b`` — ``verify_geometry.self_supporting``.
* **Thickness.** No free-standing member thinner than ``inp.thickness_mm``.
* **Floating bodies.** The lattice's own floating-body rule has already run.
  A support body attached to the wall is exempt from it (specification.md §5);
  any other new body under ``t³`` must not be produced.
* **Valid, writable solids.** Closed, ``BRepCheck_Analyzer``-valid, and
  surviving a STEP round trip — ``validate`` refuses the run otherwise.
* **Unified.** ``simplify`` has already run, so solids this stage changes
  should be unified again (``pipeline._unify`` takes a list and the pool).
"""

from __future__ import annotations

from .contract import SupportInput, SupportResult


def generate(inp: SupportInput) -> SupportResult:
    """Supports for ``inp.samples``, combined with ``inp.lattice``."""
    return SupportResult(solids=list(inp.lattice), generated=False)
