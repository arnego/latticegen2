"""The support generator (docs/algorithm.md §14.3). **Not yet implemented.**

Implementing it: follow docs/algorithm.md §14.4, step by step.

:func:`generate` receives the completed lattice and the input body
(:class:`.contract.SupportInput`) and returns the support solids alone
(:class:`.contract.SupportResult`). Until the
algorithm exists it returns none, and the run warns that its output carries no
supports.

What happens to them is not the generator's job: the support stage fuses them
with the lattice they touch, re-unifies the merged solids and checks that the
merge lost no lattice (:mod:`.merge`); then ``validate`` -> export-truth ->
``export`` run on the final product. What an implementation owes that chain:

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
* **Valid, writable solids.** Closed and ``BRepCheck_Analyzer``-valid, so the
  merge fuse is well posed, and — once merged — surviving a STEP round trip,
  which ``validate`` checks.
* **Supports only.** Never a copy of the lattice: a lattice solid returned as
  "support" would be fused with itself and counted twice in the volume
  bounds the merge is checked against.
"""

from __future__ import annotations

from .contract import SupportInput, SupportResult


def generate(inp: SupportInput) -> SupportResult:
    """Support solids backing ``inp.samples``, anchored on ``inp.lattice``."""
    return SupportResult(supports=[])
