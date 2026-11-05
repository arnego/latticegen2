"""Support solids that share a face become one body (docs/algorithm.md §14.3).

Three places in the generator hand over solids that belong to one body:

* S6, the trimmed halves of a support body that was split on a support grid
  plane, which meet on that plane;
* S6, a lost lattice piece absorbed into the supports beside it
  (``generate.absorb_lost``);
* S7, a pocket fill and the support bodies that seal it.

**The join is one boolean fuse per group of solids whose boxes touch.** The
kernel's fuse makes exactly the partition wanted: solids sharing a face or a
volume come back as one, and solids meeting only along an edge or at a vertex
stay apart, since one body joined there would not be a manifold solid. Each
fused body's coplanar faces are then merged, and the merge is kept only if it
is valid and leaves the volume unchanged.
"""

from __future__ import annotations

import numpy as np
from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse

from .. import occ
from ..errors import ProcessingError
from .assemble import SUPPORT_OVERLAP_TOL, _boxes, _common_volume, _touch, _where

JOIN_VOLUME_REL_TOL = 1e-4
"""Relative slack between a fused group's volume and its members' sum less
their overlaps. The members' own volumes carry the quadrature noise of their
curved faces: on the cup's eight cut cells the sum is 350.1275 mm^3 and both
the fuse and a stitch that moves nothing give 350.1328, 1.5e-5 apart. What the
check is for is a fuse that returns garbage -- one of its operands, as
:func:`latticegen2.boundary._fuse_solids` records -- which is orders larger."""


def groups(solids: list, grow: float = 0.0) -> list:
    """Indices of ``solids`` in groups whose bounding boxes, each inflated by
    ``grow``, touch, transitively."""
    boxes = [(lo - grow, hi + grow) for lo, hi in _boxes(solids)]
    parent = list(range(len(solids)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(len(solids)):
        for j in range(i + 1, len(solids)):
            if _touch(boxes[i], boxes[j]):
                parent[root(j)] = root(i)
    out = {}
    for i in range(len(solids)):
        out.setdefault(root(i), []).append(i)
    return list(out.values())


def _fuse(solids: list):
    algo = BRepAlgoAPI_Fuse()
    from OCP.TopTools import TopTools_ListOfShape

    args = TopTools_ListOfShape()
    args.Append(solids[0])
    tools = TopTools_ListOfShape()
    for s in solids[1:]:
        tools.Append(s)
    algo.SetArguments(args)
    algo.SetTools(tools)
    algo.Build()
    if not algo.IsDone():
        return None
    return occ.solids(algo.Shape())


def _adds_up(result, want: float) -> bool:
    if result is None or not all(occ.is_valid(s, parallel=False) for s in result):
        return False
    got = sum(occ.volume(s) for s in result)
    return abs(got - want) <= max(JOIN_VOLUME_REL_TOL * abs(want), SUPPORT_OVERLAP_TOL)


def join(solids: list, overlaps: bool = True) -> tuple[list, dict]:
    """``(bodies, stats)``: ``solids`` with every face-sharing or overlapping group fused.

    ``overlaps`` measures the volume each touching pair shares, which the
    fused volume is checked against; pass ``False`` where the members are
    disjoint by construction (pocket fills against the supports that seal
    them). The fuse is tried in one call, then folded one member at a time,
    each step checked the same way. A group whose members overlap and cannot be
    fused fails the run (exit 4): disjoint support bodies are a requirement
    (docs/algorithm.md §14.2). A group that only shares faces is left as it
    was, which is sound, and counted.
    """
    from ..pipeline import _unify_one

    stats = {"support_joined_groups": 0, "support_overlapping_pairs": 0,
             "support_joins_declined": 0}
    out = []
    boxes = _boxes(solids)
    for g in groups(solids):
        if len(g) == 1:
            out.append(solids[g[0]])
            continue
        members = [solids[i] for i in g]
        vols = [occ.volume(s) for s in members]
        shared = np.zeros((len(g), len(g)))
        if overlaps:
            for a in range(len(g)):
                for b in range(a + 1, len(g)):
                    if _touch(boxes[g[a]], boxes[g[b]]):
                        v = _common_volume(members[a], [members[b]])
                        if v > SUPPORT_OVERLAP_TOL:
                            shared[a, b] = shared[b, a] = v
                            stats["support_overlapping_pairs"] += 1
        fused = _fuse(members)
        if not _adds_up(fused, sum(vols) - shared.sum() / 2.0):
            fused, want = [members[0]], vols[0]
            for k in range(1, len(g)):
                want += vols[k] - shared[k, :k].sum()
                fused = _fuse(fused + [members[k]])
                if not _adds_up(fused, want):
                    fused = None
                    break
        if fused is None:
            if shared.any():
                raise ProcessingError(
                    f"Support bodies near {_where(boxes[g[0]])} share "
                    f"{shared.sum() / 2.0:.3e} mm^3 and could not be fused into "
                    f"valid bodies of the volume they cover "
                    f"(docs/algorithm.md §14.3 S6)."
                )
            stats["support_joins_declined"] += 1
            out += members
            continue
        if len(fused) < len(members):
            stats["support_joined_groups"] += 1
        for s in fused:
            unified, _ran, _edges = _unify_one(s)
            if occ.is_valid(unified, parallel=False) and abs(
                    occ.volume(unified) - occ.volume(s)) <= JOIN_VOLUME_REL_TOL * occ.volume(s):
                s = unified
            out.append(s)
    return out, stats
