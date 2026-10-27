"""S6: supports that share a face or volume become one body (docs/algorithm.md §14.3).

Cells that share rays along a cut or a fold share the face those rays span,
and on the cup fixture the two sides of every such face coincide to 1e-15 mm.
Cells of different regions can overlap where their rays cross under an
unsupported band between the regions -- 2.8e-4 to 0.155 mm^3 on
Hole-block-chamfered. A pocket fill shares faces with the supports that seal
it. All of them are joined here; cells that meet only along a parting, where
their side faces diverge in a V from one edge on the ceiling, stay separate
bodies, since one body joined along an edge would not be a manifold solid.

**The join is one boolean fuse per group of supports whose boxes touch.** The
kernel's fuse makes exactly that partition: solids sharing a face or volume
come back as one, solids meeting along an edge stay apart. Stitching without a
boolean -- dropping the coincident faces and sewing the rest -- was measured
first and holds only where the two sides' shared faces correspond one for one:
valid on the cup's ring of eight cut cells, 0.63 s against the fuse's 1.00 s,
and invalid, 0.85 mm^3 short, on the ledge and the lid it borders by folds,
whose shared faces the two cells subdivide differently; likewise on every pocket
fill and every overlapping pair. A path that needs the fuse behind it saves
less than its own cost in code.
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


def groups(solids: list) -> list:
    """Indices of ``solids`` in groups whose bounding boxes touch, transitively."""
    boxes = _boxes(solids)
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
