"""Merging generated supports into the lattice (docs/algorithm.md §14.2).

The generator returns support material only. This module turns lattice plus
supports into the final product: every support is fused with the lattice
solids and other supports it touches, the solids that changed are unified
again — the lattice's own were unified by ``simplify`` before the stage began —
and the result is checked against the one thing a merge can be held to without
a second opinion on the geometry: **volume**. A union can neither lose
material it was given nor gain material nobody gave it, so

    lattice volume  <=  output volume  <=  lattice volume + support volume

each within :data:`MERGE_VOLUME_REL_TOL`. The lower bound is what stops a
merge — or a generator — from silently dropping lattice (specification.md §5:
no body is dropped to make an export succeed).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCP.TopTools import TopTools_ListOfShape

from .. import occ
from ..errors import ProcessingError

MERGE_VOLUME_REL_TOL = 1e-4
"""Relative slack on both volume bounds. The same figure as
``pipeline.UNIFY_VOLUME_TOL``: the merged solids are unified again, and that
re-description is the largest legitimate drift in play (docs/algorithm.md §9)."""

CONTACT_GAP = 1e-3
"""Millimetres. Bounding boxes closer than this are treated as touching when
grouping solids for the fuse. Over-grouping only costs a larger fuse; a pair
that does not actually touch comes back from it as separate solids."""


@dataclass
class MergeResult:
    solids: list
    """The final product: lattice and supports, merged and unified."""
    lattice_volume: "float | None"
    """``None`` when there was nothing to merge: the volume is then the
    output's, which ``validate`` measures."""
    support_volume: float
    """Total volume of the support solids as the generator returned them,
    overlap with the lattice included."""
    output_volume: "float | None"
    stats: dict = field(default_factory=dict)


def _boxes(solids):
    out = []
    for s in solids:
        lo, hi = occ.bounding_box(s)
        out.append((np.asarray(lo) - CONTACT_GAP, np.asarray(hi) + CONTACT_GAP))
    return out


def _groups(n_lattice: int, boxes) -> list[list[int]]:
    """Connected groups of solids whose boxes touch, each containing a support.

    Indices below ``n_lattice`` are lattice solids. A lattice solid no support
    reaches is in no group: it passes through untouched, as the same object.
    """
    parent = list(range(len(boxes)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(n_lattice, len(boxes)):
        lo_i, hi_i = boxes[i]
        for j in range(len(boxes)):
            if j == i:
                continue
            lo_j, hi_j = boxes[j]
            if np.all(lo_i <= hi_j) and np.all(lo_j <= hi_i):
                parent[find(i)] = find(j)
    by_root: dict[int, list[int]] = {}
    for i in range(len(boxes)):
        by_root.setdefault(find(i), []).append(i)
    return [sorted(g) for g in by_root.values() if any(i >= n_lattice for i in g)]


def _fuse(solids: list, where: str) -> list:
    if len(solids) == 1:
        return list(solids)
    algo = BRepAlgoAPI_Fuse()
    args = TopTools_ListOfShape()
    args.Append(solids[0])
    tools = TopTools_ListOfShape()
    for s in solids[1:]:
        tools.Append(s)
    algo.SetArguments(args)
    algo.SetTools(tools)
    algo.Build()
    if not algo.IsDone():
        raise ProcessingError(
            f"Merging the supports with the lattice {where} did not complete: "
            f"the fuse failed (docs/algorithm.md §14.2)."
        )
    algo.SimplifyResult()
    return occ.solids(algo.Shape())


def check_volumes(lattice_volume: float, support_volume: float,
                  output_volume: float, where: str = "") -> None:
    """Fail (exit 4) unless ``lattice <= output <= lattice + support``."""
    slack = MERGE_VOLUME_REL_TOL * max(lattice_volume + support_volume, 1e-12)
    if output_volume < lattice_volume - slack:
        raise ProcessingError(
            f"Adding supports{where} lost lattice: the output is "
            f"{output_volume:.4f} mm^3 where the lattice alone was "
            f"{lattice_volume:.4f} mm^3. Supports only ever add material, so "
            f"this is lattice dropped by the merge or the generator "
            f"(specification.md §5, docs/algorithm.md §14.2)."
        )
    if output_volume > lattice_volume + support_volume + slack:
        raise ProcessingError(
            f"Adding supports{where} gained material nobody generated: the "
            f"output is {output_volume:.4f} mm^3 against "
            f"{lattice_volume:.4f} mm^3 of lattice and {support_volume:.4f} "
            f"mm^3 of support (docs/algorithm.md §14.2)."
        )


def merge_with_lattice(lattice: list, supports: list, *, unify,
                       report=None) -> MergeResult:
    """Fuse ``supports`` into ``lattice`` and unify what changed.

    ``unify(solids) -> (solids, stats)`` is the pipeline's same-domain
    unification (``pipeline._unify`` bound to the run's pool), so a merged
    solid is compacted, guarded and degraded exactly as the lattice was.
    """
    if not supports:
        # Nothing to merge, and nothing to measure: the lattice passes through
        # as the same objects, and `validate` measures its volume anyway.
        return MergeResult(list(lattice), None, 0.0, None,
                           {"support_merge_groups": 0})
    lattice_volumes = [occ.volume(s) for s in lattice]
    support_volumes = [occ.volume(s) for s in supports]
    lattice_volume = float(sum(lattice_volumes))
    support_volume = float(sum(support_volumes))

    everything = list(lattice) + list(supports)
    volumes = lattice_volumes + support_volumes
    n = len(lattice)
    groups = _groups(n, _boxes(everything))
    grouped = {i for g in groups for i in g}
    passed = [s for i, s in enumerate(lattice) if i not in grouped]
    changed = []
    for gi, group in enumerate(groups):
        if report is not None:
            report("merging supports with the lattice", gi, len(groups))
        members = [everything[i] for i in group]
        lo, hi = _boxes(members[:1])[0]
        where = f" near {np.round((lo + hi) / 2.0, 3).tolist()}"
        fused = _fuse(members, where)
        check_volumes(
            sum(volumes[i] for i in group if i < n),
            sum(volumes[i] for i in group if i >= n),
            sum(occ.volume(s) for s in fused), where,
        )
        changed.extend(fused)
    if report is not None:
        report("merging supports with the lattice", len(groups), len(groups))
        report("unifying merged solids", 0, None)
    unified, ustats = unify(changed)
    solids = passed + list(unified)
    output_volume = float(sum(volumes[i] for i in range(n) if i not in grouped)
                          + sum(occ.volume(s) for s in unified))
    check_volumes(lattice_volume, support_volume, output_volume)
    return MergeResult(solids, lattice_volume, support_volume, output_volume, {
        "support_merge_groups": len(groups),
        "support_merged_solids": len(unified),
        "support_unify_volume_drift": f"{ustats.get('volume_drift', 0.0):.2e}",
    })
