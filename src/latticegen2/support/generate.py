"""The support generator (docs/algorithm.md §14.3). **Not yet implemented.**

Being implemented in phases (specification.md §10).

:func:`generate` receives the completed lattice and the input body
(:class:`.contract.SupportInput`) and returns the support solids alone
(:class:`.contract.SupportResult`). Until the
algorithm exists it returns none, and the run warns that its output carries no
supports.

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

from .contract import SupportInput, SupportResult


def generate(inp: SupportInput) -> SupportResult:
    """Support solids backing ``inp.samples``, anchored on ``inp.lattice``."""
    return SupportResult(supports=[])
