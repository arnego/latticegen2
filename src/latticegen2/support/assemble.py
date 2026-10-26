"""Putting the supports beside the lattice (docs/algorithm.md §14.2).

The generator returns support material only, and the stage does **not** merge
it with the lattice: the output carries the lattice solids exactly as
``simplify`` left them and, after them, the support bodies, built to share faces
with the lattice and the input's surfaces so the merge in the downstream CAD
tool succeeds (specification.md §4.6).

What a separate body can be held to without a second opinion on its geometry is
that it **overlaps nothing**: not the lattice, and not another support. Volumes
are then additive by construction, which is what the run summary reports.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common
from OCP.TopAbs import TopAbs_ShapeEnum
from OCP.TopTools import TopTools_ListOfShape

from .. import lattice as _lattice
from .. import occ
from ..errors import ProcessingError

SUPPORT_OVERLAP_TOL = 1e-6
"""mm³. The most volume a support may share with the lattice or with another
support. Faces may coincide; volume may not."""

CONTACT_GAP = 1e-3
"""mm. Bounding boxes closer than this are treated as touching when choosing
which pairs to intersect. Over-selecting costs a boolean, never a wrong answer."""


@dataclass
class AssembleResult:
    solids: list
    """The final product: the lattice solids, then the support bodies."""
    lattice_volume: "float | None"
    """``None`` when there were no supports: the volume is then the output's,
    which ``validate`` measures."""
    support_volume: float
    stats: dict = field(default_factory=dict)


def check_result(result, lattice: list) -> None:
    """Fail (exit 4) unless ``result.supports`` are solids and none is the lattice."""
    check_supports(list(result.supports), lattice)


def check_supports(supports: list, lattice: list) -> None:
    """Fail (exit 4) unless every support is a solid and none is the lattice."""
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


def _common_volume(a, tools: list) -> float:
    if not tools:
        return 0.0
    algo = BRepAlgoAPI_Common()
    args = TopTools_ListOfShape()
    args.Append(a)
    tl = TopTools_ListOfShape()
    for s in tools:
        tl.Append(s)
    algo.SetArguments(args)
    algo.SetTools(tl)
    algo.Build()
    if not algo.IsDone():
        raise ProcessingError(
            "Checking a support for overlap did not complete: the intersection "
            "failed (docs/algorithm.md §14.2)."
        )
    return float(sum(occ.volume(s) for s in occ.solids(algo.Shape())))


def _boxes(solids):
    out = []
    for s in solids:
        lo, hi = occ.bounding_box(s)
        out.append((np.asarray(lo) - CONTACT_GAP, np.asarray(hi) + CONTACT_GAP))
    return out


def _touch(box_a, box_b) -> bool:
    return bool(np.all(box_a[0] <= box_b[1]) and np.all(box_b[0] <= box_a[1]))


def _where(box) -> list:
    return np.round((box[0] + box[1]) / 2.0, 3).tolist()


class _LocalLattice:
    """Junction instances near a box: the lattice's material, by construction.

    Built from the junction template at the kept nodes rather than read off the
    lattice solids, whose dominant body has faces in the hundreds of thousands:
    an intersection against it per support would put that body on the support
    path once per cell. The template's faces are the lattice's own faces
    (docs/algorithm.md §3.2), so the two answer the overlap question alike
    inside the body, which is where supports are.
    """

    def __init__(self, lp, kept_nodes: np.ndarray):
        self.lp = lp
        self.positions = (_lattice.nodes(lp, kept_nodes) if len(kept_nodes)
                          else np.empty((0, 3)))
        self._template = None
        # A junction reaches a/2 along a strut, plus the profile's circumradius.
        self.reach = lp.a / 2.0 + lp.r

    def near(self, box) -> list:
        if not len(self.positions):
            return []
        lo, hi = box
        hit = np.all((self.positions >= lo - self.reach)
                     & (self.positions <= hi + self.reach), axis=1)
        if not hit.any():
            return []
        if self._template is None:
            from ..junction import build_template

            self._template = build_template(self.lp).solid
        return [self._template.Moved(occ.translation(p))
                for p in self.positions[hit]]


def assemble(lattice: list, supports: list, lp, kept_nodes: np.ndarray,
             report=None) -> AssembleResult:
    """Check ``supports`` overlap nothing, then return lattice + supports."""
    if not supports:
        # Nothing to add and nothing to measure: the lattice passes through as
        # the same objects, and `validate` measures its volume anyway.
        return AssembleResult(list(lattice), None, 0.0, {})
    check_supports(supports, lattice)

    local = _LocalLattice(lp, np.asarray(kept_nodes))
    boxes = _boxes(supports)
    n = len(supports)
    for i, s in enumerate(supports):
        if report is not None:
            report("checking supports against the lattice", i, n)
        overlap = _common_volume(s, local.near(boxes[i]))
        if overlap > SUPPORT_OVERLAP_TOL:
            # The junction instances include material the floating-body rule
            # may have dropped from a kept node, so a reading here is confirmed
            # against the lattice solids themselves before it is believed.
            real = _common_volume(
                s, [l for l in lattice if _touch(boxes[i], _boxes([l])[0])])
            if real > SUPPORT_OVERLAP_TOL:
                raise ProcessingError(
                    f"A support body overlaps the lattice by {real:.3e} mm^3 near "
                    f"{_where(boxes[i])}; supports must share faces with the "
                    f"lattice, never volume (docs/algorithm.md §14.2)."
                )
    pairs = 0
    for i in range(n):
        if report is not None:
            report("checking supports against each other", i, n)
        for j in range(i + 1, n):
            if not _touch(boxes[i], boxes[j]):
                continue
            pairs += 1
            overlap = _common_volume(supports[i], [supports[j]])
            if overlap > SUPPORT_OVERLAP_TOL:
                raise ProcessingError(
                    f"Two support bodies overlap by {overlap:.3e} mm^3 near "
                    f"{_where(boxes[i])}; neighbouring supports must share faces, "
                    f"never volume (docs/algorithm.md §14.2)."
                )
    if report is not None:
        report("checking supports against each other", n, n)

    lattice_volume = float(sum(occ.volume(s) for s in lattice))
    support_volume = float(sum(occ.volume(s) for s in supports))
    return AssembleResult(list(lattice) + list(supports), lattice_volume,
                          support_volume, {
                              "support_bodies": n,
                              "support_pairs_checked": pairs,
                          })
