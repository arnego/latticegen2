"""The fixed numbers of the support algorithm (docs/algorithm.md §14).

None of these is a parameter (specification.md §4.6). They live in one place so
that the detector, the generator and the verifier in ``tools/`` cannot drift
apart; the two angles come from :mod:`latticegen2.cli`, which the window can
import without the geometry kernel.
"""

from __future__ import annotations

import math

from ..cli import SUPPORT_OVERHANG_DEG, SUPPORT_START_DEG

STEP = 0.5
"""``h``, mm: the support step. Spacing of overhang samples and of ray sources,
and the edge length the region mesh is refined to (§14.3 S1)."""

REQUIRED_SLACK = 1e-6
"""How far past ``sin(SUPPORT_START_DEG)`` a sample's ``n . b`` must be to count
as over the bar rather than at it — so a surface tilted at exactly 60 deg, like
`Hole-block-chamfered`'s 30 deg roof facet, is not over it by a rounding error."""

SIN_START = math.sin(math.radians(SUPPORT_START_DEG))
"""``sin(60 deg)``: a point needs support iff ``n . b > SIN_START``."""

TAN_CONE = math.tan(math.radians(SUPPORT_OVERHANG_DEG))
"""Tangent of the support cone's half-angle, the lattice rod angle: exactly
``sqrt(2)``, computed from the expression rather than typed."""

__all__ = ["STEP", "REQUIRED_SLACK", "SIN_START", "TAN_CONE",
           "SUPPORT_OVERHANG_DEG", "SUPPORT_START_DEG"]
