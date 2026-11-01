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
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.TopAbs import TopAbs_ShapeEnum, TopAbs_State
from OCP.TopTools import TopTools_ListOfShape
from OCP.gp import gp_Pnt

from .. import lattice as _lattice
from .. import occ
from ..errors import ProcessingError

SUPPORT_OVERLAP_TOL = 1e-6
"""mm³. The most volume a support may share with the lattice or with another
support. Faces may coincide; volume may not."""

KEY_OFFSET = 1 << 20
KEY_SPAN = 1 << 21
"""Packing of a node index triple into one integer: each index shifted by
:data:`KEY_OFFSET` into ``[0, KEY_SPAN)``, a million cells either way of the
origin -- 0.4 mm cells reach 400 m."""

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

    ``lost`` is the exception: pieces of a kept node's junction that are **not**
    in the lattice, because the floating-body rule dropped them (specification.md
    §5). A trim can cut a strut's tip off from its node; the node stays kept
    through its other arms, and the tip is gone. Measured on
    `TD_HX_rehearsal_test`: 146 ceiling samples under such tips, read as backed
    by a strut that is not there. Neither :meth:`contains` nor :meth:`near`
    counts that material.
    """

    def __init__(self, lp, kept_nodes: np.ndarray, lost=()):
        self.lp = lp
        self._nodes = np.asarray(kept_nodes, dtype=np.int64).reshape(-1, 3)
        self.positions = (_lattice.nodes(lp, kept_nodes) if len(kept_nodes)
                          else np.empty((0, 3)))
        self._template = None
        self._keys = None
        # A junction reaches a/2 along a strut, plus the profile's circumradius.
        self.reach = lp.a / 2.0 + lp.r
        self.lost = list(lost)
        self.lost_boxes = [tuple(np.asarray(x, dtype=float) for x in occ.bounding_box(s))
                           for s in self.lost]
        self._instances = {}

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
        return [self._instance(int(i)) for i in np.nonzero(hit)[0]]

    def _instance(self, i: int):
        """Junction ``i``'s template, less any lost piece of it.

        A cut that does not take out what it should leaves the template: the
        overlap check then falls back on the lattice solids, and the pocket
        search on its refusal path, so neither is misled."""
        got = self._instances.get(i)
        if got is not None:
            return got
        p = self.positions[i]
        got = self._template.Moved(occ.translation(p))
        for s, (lo, hi) in zip(self.lost, self.lost_boxes):
            if np.any(lo > p + self.reach) or np.any(hi < p - self.reach):
                continue
            v = occ.volume(got)
            cut = BRepAlgoAPI_Cut(got, s)
            if not cut.IsDone():
                continue
            gone = v - occ.volume(cut.Shape())
            if -1e-9 <= gone <= occ.volume(s) * (1.0 + 1e-3):
                got = cut.Shape()
        self._instances[i] = got
        return got

    def contains(self, points) -> np.ndarray:
        """Which ``points`` lie in the lattice's material: the junction
        instances', less what was lost. :meth:`in_template` is the first half."""
        P = np.asarray(points, dtype=float).reshape(-1, 3)
        out = self.in_template(P)
        if not self.lost or not out.any():
            return out
        for s, (lo, hi) in zip(self.lost, self.lost_boxes):
            at = np.nonzero(out & np.all((P >= lo - 1e-6) & (P <= hi + 1e-6), axis=1))[0]
            if not len(at):
                continue
            cls = BRepClass3d_SolidClassifier(s)
            for i in at:
                cls.Perform(gp_Pnt(*map(float, P[i])), 1e-7)
                if cls.State() in (TopAbs_State.TopAbs_IN, TopAbs_State.TopAbs_ON):
                    out[i] = False
        return out

    def in_template(self, points) -> np.ndarray:
        """Which ``points`` lie in the junction instances' material, exactly.

        ``J`` is the union of six half-strut prisms (docs/algorithm.md §3.2), so
        a point is inside it iff, for some strut family ``k``, it lies within
        ``a/2`` of the node along ``e_k`` and inside the diamond profile across
        it: ``|x·u_k| + |x·v_k| <= r``. A point within a junction is less than
        one cell from its node in every index coordinate, so the nodes to ask
        are the corners of the index cell holding it.
        """
        P = np.asarray(points, dtype=float).reshape(-1, 3)
        out = np.zeros(len(P), dtype=bool)
        if not len(self.positions) or not len(P):
            return out
        if self._keys is None:
            self._keys = np.sort(self._key(self._nodes))
        lp = self.lp
        q = np.linalg.solve(lp.B, P.T).T
        base = np.floor(q).astype(np.int64)
        for corner in np.array(np.meshgrid([0, 1], [0, 1], [0, 1])).T.reshape(-1, 3):
            n = base + corner
            key = self._key(n)
            at = np.searchsorted(self._keys, key)
            kept = (at < len(self._keys)) & (self._keys[np.minimum(at, len(self._keys) - 1)] == key)
            if not kept.any():
                continue
            d = P[kept] - (lp.B @ n[kept].T).T
            hit = np.zeros(len(d), dtype=bool)
            for k in range(3):
                hit |= ((np.abs(d @ lp.e[k]) <= lp.a / 2.0)
                        & (np.abs(d @ lp.u[k]) + np.abs(d @ lp.v[k]) <= lp.r))
            out[np.nonzero(kept)[0][hit]] = True
        return out

    @staticmethod
    def _key(n: np.ndarray) -> np.ndarray:
        n = np.asarray(n, dtype=np.int64) + KEY_OFFSET
        return (n[:, 0] * KEY_SPAN + n[:, 1]) * KEY_SPAN + n[:, 2]


def assemble(lattice: list, supports: list, lp, kept_nodes: np.ndarray,
             report=None, lost=()) -> AssembleResult:
    """Check ``supports`` overlap nothing, then return lattice + supports."""
    if not supports:
        # Nothing to add and nothing to measure: the lattice passes through as
        # the same objects, and `validate` measures its volume anyway.
        return AssembleResult(list(lattice), None, 0.0, {})
    check_supports(supports, lattice)

    local = _LocalLattice(lp, np.asarray(kept_nodes), lost)
    boxes = _boxes(supports)
    n = len(supports)
    for i, s in enumerate(supports):
        if report is not None:
            report("checking supports against the lattice", i, n)
        overlap = _common_volume(s, local.near(boxes[i]))
        if overlap > SUPPORT_OVERLAP_TOL:
            # A junction instance is the template, untrimmed, so a reading here
            # is confirmed against the lattice solids themselves before it is
            # believed.
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
