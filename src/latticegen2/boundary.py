"""Boundary junctions — the only place a boolean touches the input body.

Each BOUNDARY node contributes exactly one intersection: its instanced junction
solid against the input solid. That is deliberately **one object operand per
call**. OCCT's general boolean runs over all object operands together, so two
overlapping objects in one call are *partitioned* against the tool rather than
each trimmed independently. Measured directly: three struts sharing one lattice
node, intersected against a containing box in a single call, come back as **7**
fragment solids instead of the 1 solid produced by fusing them first. A single
already-fused junction cannot trigger that, by construction — which is why this
module needs no machinery to keep operands disjoint (docs/algorithm.md §7).

After trimming, every face lying in a cap plane is **tagged, not dropped**, and
the decision to drop it is taken later on the master by
:func:`resolve_interfaces`, once both sides of every cap are known.

That split is the whole point. A dropped cap is a hole this junction punches for
its neighbour to fill, so it is only sound if the neighbour punched the matching
one. The two sides are two *independent* ``BRepAlgoAPI_Common`` calls, in
different processes, against the same nominal cap quad — and OCCT gives no
guarantee that a shared face comes back the same from both. Deciding locally
from classification alone ("the node across it is kept, so drop") let a junction
punch a hole with nothing behind it wherever the two booleans disagreed — and
nothing between the trim and the stitcher looked, so the consequence surfaced
only at the very end, as an unclosed shell with no indication of where. That is
the leading explanation for the ``cc=5, t=1`` rehearsal ending in "1 of 14
stitched shells are not closed" after 5 hours. Resolving interfaces symmetrically
makes "every hole has a partner" true by construction, and any residual
inconsistency fails immediately and by name (docs/algorithm.md §5.3, §7.1, §8).

The work is embarrassingly parallel — constant-size, independent jobs — so it is
distributed over worker processes under a small-IPC discipline: the input body
goes to disk once as a ``.brep`` and workers read it directly; only file paths
and small metadata cross the process boundary, never geometry.
"""

from __future__ import annotations

import itertools
import os
from dataclasses import dataclass, field
from typing import NamedTuple

import numpy as np
from OCP.BRep import BRep_Builder
from OCP.BRep import BRep_Tool
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Fuse
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.BRepTools import BRepTools
from OCP.TopAbs import TopAbs_ShapeEnum
from OCP.TopoDS import TopoDS, TopoDS_Shape, TopoDS_Shell
from OCP.gp import gp_Pnt
from OCP.TopTools import TopTools_ListOfShape

from . import occ
from .connect import UnionFind
from .errors import ProcessingError
from .junction import JunctionTemplate, build_template, is_cap_plane_face
from .lattice import (
    OPPOSITE_HALF,
    LatticeParams,
    half_strut_offset,
    lattice_params,
    neighbor_step,
    nodes,
)
from .parallel import WorkerPool
from .parallel import compound_children as _compound_children
from .parallel import read_brep as _read_brep

NodeKey = tuple[int, int, int]

CAP_AREA_REL_TOL = 1e-6
"""Relative agreement two sides of a cap must show before it is stitched across.

The two sides are the same nominal region computed by two independent booleans,
so quadrature noise is the only difference expected and it lands far below this;
the same bar :func:`latticegen2.junction._identify_caps` uses on the template.
A real disagreement means one trim kept material the other did not, and stitching
across it would weld two boundaries that are not the same curve."""

PINHOLE_WIRE_TOL = 3e-5
"""Millimetres. An inner wire whose every edge is shorter than this, and whose
edges are already unpaired, bounds no area and is removed (docs/algorithm.md §7).

The two the `TD_HX_rehearsal_test` rehearsal leaves at `cc=5, t=1` are 3.171690e-06
and 5.808982e-06 mm, so this sits ~5x above the largest of them. It does not
need the wide margin a length threshold usually would, because length is not
what makes a wire removable here: :func:`latticegen2.occ.remove_pinhole_wires`
additionally requires that every edge of the wire is used exactly once, i.e. is
already a defect rather than a shared boundary. A real feature is paired and is
therefore out of reach of this repair at any threshold. See
``tools/prototypes/RESULTS.md`` G10."""

PINHOLE_AREA_TOL = 1e-12
"""Relative surface-area change allowed when pinhole wires are removed.

Deliberately near zero rather than a comfortable margin, and that is sound
*here* specifically: a wire bounding no area cannot change the area of the face
that carried it, so unlike same-domain unification (docs/algorithm.md §9) there
is no larger merged region to re-integrate and no quadrature noise to absorb.
Measured on the piece this was built for, the drift is 0.0 — bit-identical, not
merely small. This bar exists to catch a wire that turned out to bound
something, which is the only way this repair could go wrong."""

# There is deliberately no volume tolerance here, and the absence is
# load-bearing enough to be worth a note where the one that used to live here
# was. A relative-volume bar of 1e-9 stood here until it refused a valid run
# (`-cc 12 -t 2.5` on `TD_HX_rehearsal_test`, drift 1.235e-09 on a junction of
# 77.4 mm³) — not because anything moved, but because OCCT cannot measure the
# volume of the *unrepaired* piece: `BRepGProp::VolumeProperties` requires a
# shape "exempt of any free boundary", and a pinhole wire is one by definition.
# The drift is that measurement's bias, not the repair's effect, and G19
# measured it length-independent over three decades and set per face, so no bar
# expressed relative to volume, `t` or junction size can bound it.
# :func:`latticegen2.occ.only_inner_wires_dropped` replaces it and is exact.
# See docs/algorithm.md §7.


@dataclass
class BoundaryPiece:
    """One connected solid produced by trimming one boundary junction.

    ``node`` is the piece's own trimmed junction, or — after
    :func:`fuse_disagreeing_pairs` has merged two disagreeing pieces into one —
    an arbitrary representative of the nodes it now spans. Nothing downstream
    reads ``node`` to attribute a specific cap; ``caps`` and ``cap_faces`` carry
    the owning node explicitly for exactly that reason, so a fused piece can
    hold faces belonging to more than one node without ambiguity.
    """

    node: NodeKey
    volume: float
    tolerance_ratio: float = 0.0
    """Worst ``edge tolerance / sqrt(face area)`` over this piece's faces.

    Measured in the worker, on the piece as the boolean produced it, because
    this is the one moment at which the junction that produced it is still
    named. See :func:`latticegen2.occ.tolerance_feature_ratio` and
    docs/algorithm.md §7.3."""
    tolerance_evidence: tuple[float, float, tuple[float, float, float]] = (
        0.0, 0.0, (0.0, 0.0, 0.0)
    )
    """``(tolerance, face area, face centroid)`` behind ``tolerance_ratio``."""
    faces: list = field(default_factory=list)
    """Every face of the trimmed piece that does not lie in a cap plane."""
    cap_faces: dict[tuple[NodeKey, int], list] = field(default_factory=dict)
    """``(node, half-strut id)`` -> the cap-plane face(s) the trim left there.

    A trim can leave several disjoint fragments of one cap on the same piece, so
    this is a list per key, not a single face. Keyed by node as well as
    half-strut id because a piece :func:`fuse_disagreeing_pairs` has merged can
    hold caps belonging to either of the nodes it spans."""
    caps: frozenset[tuple[NodeKey, int]] = frozenset()
    """``(node, half-strut id)`` pairs whose cap face this piece gave up as an
    interface.

    Filled by :func:`finalize_pieces`. This is the piece's contribution to the
    junction graph (:mod:`latticegen2.connect`), and by construction every entry
    has a matching cap on the other side."""


@dataclass
class BoundaryResult:
    pieces: list[BoundaryPiece] = field(default_factory=list)
    n_empty: int = 0
    """Junctions whose intersection with the input body was empty."""
    n_pinhole_junctions: int = 0
    """Junctions that carried zero-area pinhole wires and were repaired."""
    n_pinhole_wires: int = 0
    """Pinhole wires removed in total (docs/algorithm.md §7).

    Logged as one aggregate line rather than one per junction, the same way the
    floating-body rule reports its removals: a part dense in grazing trims can
    produce many, and a line each would bury the rest of the log."""
    n_retrimmed_junctions: int = 0
    """Junctions whose fused intersection came back untrimmed and had to be
    redone per half-strut (docs/algorithm.md §7.2).

    Aggregated for the same reason the pinhole count is: a part that provokes
    it once tends to provoke it repeatedly, and one line each would bury the
    log. Zero on every committed scenario except `spiral-stress`."""
    n_localized_junctions: int = 0
    """Junctions re-trimmed against a locally cut block of the body because the
    whole-body intersection left material outside it (docs/algorithm.md §7.2)."""
    n_dropped_junctions: int = 0
    """Junctions discarded because no intersection kept them inside the body
    (docs/algorithm.md §7.2, :attr:`TrimResult.dropped`)."""
    worst_outside_mm: float = 0.0
    """Furthest outside the body any *discarded* junction reached. Nothing in
    the output reaches outside at all: specification.md §1 requires the lattice
    to fit exactly within the input, so what could not be trimmed was dropped."""
    worst_outside_node: tuple | None = None
    max_worker_rss: int = 0
    """Highest peak RSS reported by any worker process.

    Reported separately from the master's own peak because they are different
    processes: the run's true memory footprint is the master's peak plus what
    the workers held concurrently, and specification.md §3 asks for maximum
    memory usage.
    """
    diagnostics: list[str] = field(default_factory=list)

    def worst_tolerance_pieces(self, n: int = 5) -> list[BoundaryPiece]:
        """The ``n`` pieces whose description leans hardest on tolerance.

        Worst first. This is the source-side half of the export-truth question
        (docs/algorithm.md §7.3): a piece near the top of this list was trimmed
        with a slack approaching the size of the feature it bounds, and slack is
        exactly what §9's export cannot carry.
        """
        return sorted(self.pieces, key=lambda p: p.tolerance_ratio, reverse=True)[:n]


@dataclass
class InterfaceSet:
    """Which caps are stitched across, and what was rejected on the way there."""

    interfaces: set[tuple[NodeKey, int]] = field(default_factory=set)
    """Both sides of every stitched cap, so either side can test membership
    directly: ``(node, h)`` is present exactly when ``(node + step(h),
    OPPOSITE_HALF[h])`` is."""
    unpaired: list[tuple[NodeKey, int]] = field(default_factory=list)
    """Caps whose neighbouring junction produced geometry but no matching cap.

    Each of these is a hole the previous implementation would have punched with
    nothing behind it. They are kept as exterior surface instead."""
    mismatched: list[tuple[NodeKey, int, float, float]] = field(default_factory=list)
    """``(node, h, area, partner area)`` for caps both sides presented but with
    regions that disagree beyond :data:`CAP_AREA_REL_TOL`."""
    unweldable: list[tuple[NodeKey, int]] = field(default_factory=list)
    """Caps whose two holes do not correspond edge for edge
    (:func:`latticegen2.weld.unweldable`)."""

    @property
    def n_pairs(self) -> int:
        return len(self.interfaces) // 2

    def decline(self, node: NodeKey, h: int) -> None:
        """Withdraw both sides of an interface, before either gives up its cap.

        Called while a rejection is still free: once a cap face has been dropped,
        putting it back is the kind of undo that goes wrong quietly. Both sides
        keep their cap and the output carries one more solid (docs/algorithm.md
        §11).
        """
        step = neighbor_step(h)
        other = _neighbour(node, tuple(int(x) for x in step))
        self.interfaces.discard((node, h))
        self.interfaces.discard((other, OPPOSITE_HALF[h]))
        self.unweldable.append((node, h))


def _open_shell(faces) -> TopoDS_Shell:
    shell = TopoDS_Shell()
    builder = BRep_Builder()
    builder.MakeShell(shell)
    for f in faces:
        builder.Add(shell, f)
    return shell


class TrimResult(NamedTuple):
    """What one junction's trim produced."""

    pieces: list
    """``(faces, tags, volume)`` per connected solid the intersection left."""
    n_pinholes_removed: int
    """Zero-area pinhole wires dropped across all of them (docs/algorithm.md §7)."""
    tolerances: list
    """:class:`latticegen2.occ.ToleranceFeature` per piece, aligned with ``pieces``.

    A separate list rather than a fourth field of each piece tuple, because
    those tuples are what every caller and test in this codebase destructures
    and the reading is not part of the geometry they describe. Alignment is by
    construction: both lists are appended to once per solid, in one loop."""
    retrimmed: bool = False
    """Whether the fused junction's intersection had to be redone per half-strut
    (:func:`_retrim_per_half`, docs/algorithm.md §7.2)."""
    localized: bool = False
    """Whether the trim had to be redone against a locally cut block of the body
    (:func:`_retrim_against_local_block`, docs/algorithm.md §7.2)."""
    outside_mm: float = 0.0
    """How far the junction reached outside the body when it was given up, or
    zero. Non-zero implies :attr:`dropped`."""
    dropped: bool = False
    """Whether this junction was discarded because no intersection available to
    this stage would keep it inside the body (docs/algorithm.md §7.2).

    **Discarding is the safe direction and failing the run is not.**
    specification.md §1 requires the lattice to fit *exactly* within the input,
    so geometry that reaches outside it cannot be shipped; but a junction that
    contributes nothing costs a fraction of the lattice, while a hard failure
    costs the whole part. It is structurally the same as the junctions whose
    intersection is legitimately empty — of which this part already has 215 —
    so the neighbours' caps simply stay closed and the output stays watertight.
    Reported loudly, never silently: a run that has to do this should say so.
    """


LOCAL_BLOCK_CELLS = 1.0
"""Half-size of the block cut out of the input body before re-trimming, in cells.

The block has to be **tight**. Measured on `SpiralTest.step`, against Monte
Carlo ground truth: at one cell the re-trim lands within 1.2 sigma of the truth
at every junction tried, and from four cells outward it reproduces the wrong
answer the whole-body intersection gave. Two cells works at some junctions and
not others. There is no margin to widen here and no version of this that shares
one block across a worker's batch — that idea is what the measurement killed.

A junction reaches ``a/2``, so one cell either side of the node contains it with
room to spare.
"""


# How close to the untrimmed junction's own volume an intersection result has to
# be before it is treated as "the kernel returned its operand unchanged". The
# failing case comes back re-partitioned rather than byte-identical -- measured
# 5.7e-08 mm^3 from the template's own figure on `SpiralTest`, a relative
# 6.6e-09 -- so this has to sit above quadrature noise while staying far below
# any real trim. It is not load-bearing either way: anything inside it is
# re-derived per half-strut, and a junction that really is whole comes back
# whole and keeps its original result.
UNTRIMMED_VOLUME_REL_TOL = 1e-6


def _piece_vertices(pieces) -> np.ndarray:
    """Every distinct vertex position across a trim's pieces."""
    pts = []
    for faces, _tags, _vol in pieces:
        for face in faces:
            for v in occ._explore(face, TopAbs_ShapeEnum.TopAbs_VERTEX):
                p = BRep_Tool.Pnt_s(TopoDS.Vertex_s(v))
                pts.append((p.X(), p.Y(), p.Z()))
    if not pts:
        return np.empty((0, 3))
    return np.unique(np.array(pts), axis=0)


def _retrim_against_local_block(
    lp: LatticeParams,
    tpl: JunctionTemplate,
    node_pos: np.ndarray,
    body: TopoDS_Shape,
) -> TopoDS_Shape | None:
    """Re-trim one junction against a block cut out of the body around it.

    The residual case §7.2's per-half-strut repair cannot reach: the fused
    intersection and the six half-strut ones **agree with each other** and are
    both wrong, so there is no second construction left to compare against and
    no way to grade either with another boolean against the same body.

    What does work is giving the kernel a smaller tool. A box/body intersection
    is well conditioned at exactly the nodes where the junction/body one is not
    — measured — so the body is cut down to a :data:`LOCAL_BLOCK_CELLS` block
    around the node first and the junction is intersected with that. Checked
    against Monte Carlo ground truth built from point classification alone:

    ==================  ==============  ====================  ==============
    junction            whole body      truth (Monte Carlo)   local block
    ==================  ==============  ====================  ==============
    cc=5, (6,4,6)       8.3772 (5 s.d.) 8.010 +/- 0.074       **8.0662**
    cc=7, (2,4,2)       23.6165 (whole) 11.513 +/- 0.398      **11.0446**
    ==================  ==============  ====================  ==============

    Returns ``None`` if the block cannot be cut, leaving the caller to report
    the junction rather than silently keeping geometry that reaches outside.
    """
    half = LOCAL_BLOCK_CELLS * lp.a
    lo = node_pos - half
    box = BRepPrimAPI_MakeBox(gp_Pnt(*lo), 2 * half, 2 * half, 2 * half).Shape()

    cut = BRepAlgoAPI_Common()
    args = TopTools_ListOfShape()
    args.Append(box)
    tools = TopTools_ListOfShape()
    tools.Append(body)
    cut.SetArguments(args)
    cut.SetTools(tools)
    cut.Build()
    if not cut.IsDone() or not occ.solids(cut.Shape()):
        return None
    local = cut.Shape()

    algo = BRepAlgoAPI_Common()
    args = TopTools_ListOfShape()
    args.Append(tpl.solid.Moved(occ.translation(node_pos)))
    tools = TopTools_ListOfShape()
    tools.Append(local)
    algo.SetArguments(args)
    algo.SetTools(tools)
    algo.Build()
    return algo.Shape() if algo.IsDone() else None


def _retrim_per_half(
    lp: LatticeParams,
    tpl: JunctionTemplate,
    node_pos: np.ndarray,
    body: TopoDS_Shape,
) -> TopoDS_Shape | None:
    """Redo one junction's trim on its six half-struts, or ``None`` to keep the
    original.

    ``None`` means the check cleared the original result: every half-strut came
    back whole, so the junction really is wholly inside the body and an
    intersection returning it unchanged was right.

    Otherwise the surviving half-strut intersections are fused back into one
    solid and returned. That fuse is a general boolean, which this pipeline
    otherwise performs only twice (docs/algorithm.md §3.2's template fuse and
    §7.1's disagreeing-cap repair) — it is affordable here for the same reason
    as §7.1's: it runs only where the kernel has already contradicted itself,
    measured at 9.2 % of boundary junctions reaching the check and 0.4 %
    actually needing the repair.
    """
    loc = occ.translation(node_pos)
    kept = []
    whole = True
    for h in range(6):
        algo = BRepAlgoAPI_Common()
        args = TopTools_ListOfShape()
        args.Append(tpl.half_solids[h].Moved(loc))
        tools = TopTools_ListOfShape()
        tools.Append(body)
        algo.SetArguments(args)
        algo.SetTools(tools)
        algo.Build()
        if not algo.IsDone():
            raise ProcessingError(
                f"Re-trimming half-strut {h} of the junction at "
                f"{tuple(node_pos)} failed. The fused junction's own "
                f"intersection had already returned it untrimmed, so there is "
                f"no result left to fall back on."
            )
        solids = occ.solids(algo.Shape())
        got = sum(occ.volume(sd) for sd in solids)
        want = occ.volume(tpl.half_solids[h])
        if abs(got - want) > UNTRIMMED_VOLUME_REL_TOL * want:
            whole = False
        kept.extend(solids)
    if whole:
        return None  # genuinely wholly inside: the original result stands
    if not kept:
        return occ.compound([])
    if len(kept) == 1:
        return kept[0]
    fuse = BRepAlgoAPI_Fuse()
    args = TopTools_ListOfShape()
    args.Append(kept[0])
    tools = TopTools_ListOfShape()
    for sd in kept[1:]:
        tools.Append(sd)
    fuse.SetArguments(args)
    fuse.SetTools(tools)
    fuse.Build()
    if not fuse.IsDone():
        raise ProcessingError(
            f"Could not fuse the re-trimmed half-struts of the junction at "
            f"{tuple(node_pos)} back into one solid."
        )
    fuse.SimplifyResult()
    return fuse.Shape()


def trim_junction(
    lp: LatticeParams,
    tpl: JunctionTemplate,
    node_pos: np.ndarray,
    body: TopoDS_Shape,
    probe=None,
) -> TrimResult:
    """Intersect one instanced junction with the body and tag its cap faces.

    ``pieces`` holds a ``(faces, tags, volume)`` triple per connected solid the
    intersection produced. ``tags[i]`` is the half-strut id of the cap plane
    ``faces[i]`` lies in, or ``-1`` for every other face. A boundary junction can
    legitimately split into several pieces when the input surface cuts between
    its arms; each piece becomes its own vertex in the connectivity graph.

    Caps are tagged rather than dropped because whether a cap is an interface is
    not a local fact — it depends on what the junction on the other side of that
    cap produced, which is a different boolean in a different process. See the
    module docstring and :func:`resolve_interfaces`.

    Each piece has its zero-area pinhole wires removed *before* it is tagged, so
    ``faces`` and ``tags`` are built from one pass over the repaired face list
    and stay index-aligned by construction — which is what :func:`_piece_from`
    relies on.
    """
    instance = tpl.solid.Moved(occ.translation(node_pos))
    algo = BRepAlgoAPI_Common()
    args = TopTools_ListOfShape()
    args.Append(instance)
    tools = TopTools_ListOfShape()
    tools.Append(body)
    algo.SetArguments(args)
    algo.SetTools(tools)
    algo.Build()
    if not algo.IsDone():
        raise ProcessingError(
            f"Boundary intersection failed for junction at {tuple(node_pos)}."
        )

    result = algo.Shape()
    trimmed = occ.solids(result)

    # **The kernel can return this intersection's own operand, untrimmed.**
    # Measured on `SpiralTest.step` at cc=5, t=1: three junctions come back at
    # exactly the template volume with `IsDone` true, leaving lattice material
    # up to 1.29 mm outside the input body. It is specific to the *fused*
    # junction operand -- spheres and boxes of the same scale at the same node
    # trim correctly, the input solid is `BRepCheck_Analyzer`-valid, and
    # neither a fuzzy value, swapped operands nor baking the location away
    # changes it (docs/algorithm.md §7.2). Whatever the cause, a result equal
    # to the operand is the one case that can be checked cheaply and redone, so
    # it is.
    #
    # This does not decide whether the junction is inside the body -- that is
    # what the kernel just got wrong, so asking it again the same way would be
    # worthless. `_retrim_per_half` re-derives the answer from operands the
    # same call handles correctly, and clears the original result whenever the
    # junction really is wholly inside.
    retrimmed = False
    if len(trimmed) == 1 and abs(
        occ.volume(trimmed[0]) - tpl.volume
    ) <= UNTRIMMED_VOLUME_REL_TOL * tpl.volume:
        repaired = _retrim_per_half(lp, tpl, node_pos, body)
        if repaired is not None:
            trimmed = occ.solids(repaired)
            retrimmed = True

    def finish(solids):
        pieces = []
        tolerances = []
        n_pinholes = 0
        for solid in solids:
            cleaned, removed = _remove_pinholes(node_pos, solid)
            n_pinholes += removed
            faces = occ.faces(cleaned)
            tags = [is_cap_plane_face(lp, f, node_pos) for f in faces]
            tags = [-1 if h is None else h for h in tags]
            pieces.append((faces, tags, occ.volume(cleaned)))
            # Measured here rather than on the assembled output, over the face
            # list this pass already holds: it costs one area and one centroid
            # on the single worst face, it runs in the worker for free
            # alongside the trim that produced the geometry, and it is the last
            # point at which this piece can still be reported by the junction
            # it came from (docs/algorithm.md §7.3).
            tolerances.append(occ.tolerance_feature_ratio(faces))
        return pieces, tolerances, n_pinholes

    out, tolerances, n_pinholes = finish(trimmed)

    # **Whether the trim left material outside the body is asked directly, and
    # not of a boolean.** §7.2's per-half-strut repair covers the case where
    # the fused intersection returns its operand; it cannot cover the case
    # where the fused and per-half results *agree with each other* and are both
    # wrong, because it has nothing left to compare against. `probe` answers it
    # from the classification mesh instead — ray parity plus a distance bar,
    # sharing no machinery with the kernel that produced the error — at a few
    # milliseconds per junction, so every junction can be asked.
    localized = False
    outside_mm = 0.0
    if probe is not None and out:
        outside_mm = probe.worst_outside(_piece_vertices(out))
        if outside_mm > 0.0:
            repaired = _retrim_against_local_block(lp, tpl, node_pos, body)
            if repaired is not None:
                solids = occ.solids(repaired)
                candidate, cand_tol, n_ph = finish(solids)
                # Kept only if it is actually better: the block is a smaller
                # tool for the same kernel, not a different kind of answer, so
                # it is checked by the same probe rather than assumed to work.
                still = probe.worst_outside(_piece_vertices(candidate)) if candidate else 0.0
                if still < outside_mm:
                    out, tolerances, n_pinholes = candidate, cand_tol, n_ph
                    outside_mm = still
                    localized = True
        if outside_mm > 0.0:
            # Nothing this stage can do trims it correctly, so it does not go
            # into the output at all. See `TrimResult.dropped`. The pinhole
            # count is still reported: the repair really did run on this
            # junction's faces, and docs/algorithm.md §7's aggregate line is the
            # only evidence a run gives that it ran at all -- dropping the
            # junction's geometry is no reason to drop the measurement with it.
            return TrimResult(
                [], n_pinholes, [], retrimmed, localized, outside_mm, True
            )
    return TrimResult(out, n_pinholes, tolerances, retrimmed, localized, outside_mm)


def _remove_pinholes(node_pos: np.ndarray, solid):
    """Drop this piece's zero-area pinhole wires, proving nothing else moved.

    Returns ``(solid, n_removed)``. A pinhole is an inner wire of a few microns
    that bounds no area, left by a near-tangential trim; without this the shell
    is rejected by :func:`latticegen2.weld.shell_defects` for an edge used once
    (docs/algorithm.md §7). :func:`latticegen2.occ.remove_pinhole_wires` will
    only ever delete an edge that is *already* unpaired and bounds nothing, so
    the repair cannot open a hole — but "cannot" is checked rather than trusted,
    the way every other gate in this pipeline is (docs/algorithm.md §11).

    **Surface area is the quantity checked, and it is checked exactly.**
    Removing a wire that bounds no area cannot change the area of the face that
    carried it, so unlike a merge or a unification there is no quadrature noise
    to allow for: the measured drift on the piece this repair was built for is
    0.0, bit-identical, and anything else means a wire that did bound something
    was removed.

    **Volume is deliberately not the second signal it once was.** OCCT can only
    integrate the volume of a shape "exempt of any free boundary", and the
    pinhole is a free boundary — so the pre-repair figure carries a bias that
    disappears with the wire, and comparing the two measures the defect rather
    than the repair (see the note above the constants, and G19). What stands in
    its place is :func:`latticegen2.occ.only_inner_wires_dropped`, which proves
    the same thing structurally and exactly: same faces, same outer wires, same
    orientations, and exactly the accounted-for wires gone.

    This is a correctness repair, not an optimization, so it does not degrade
    silently: skipping it leaves the shell unclosed hours later with no
    indication of where, which is exactly the failure it exists to prevent.
    """
    try:
        cleaned, n_removed = occ.remove_pinhole_wires(solid, PINHOLE_WIRE_TOL)
    except Exception as exc:                                   # noqa: BLE001
        raise ProcessingError(
            f"Pinhole-wire removal failed for the junction at {tuple(node_pos)}: {exc}"
        ) from exc
    if n_removed <= 0:
        return solid, 0

    area_before = sum(occ.area(f) for f in occ.faces(solid))
    area_after = sum(occ.area(f) for f in occ.faces(cleaned))
    drift = abs(area_after - area_before) / max(area_before, 1e-30)
    if drift > PINHOLE_AREA_TOL:
        raise ProcessingError(
            f"Removing {n_removed} pinhole wire(s) from the junction at "
            f"{tuple(node_pos)} changed its surface area from {area_before:.9g} "
            f"to {area_after:.9g} mm^2 ({drift:.3e} relative, tolerance "
            f"{PINHOLE_AREA_TOL:g}). A wire that bounds no area cannot change "
            f"it, so one that bounded something was removed."
        )

    reason = occ.only_inner_wires_dropped(solid, cleaned, n_removed)
    if reason is not None:
        raise ProcessingError(
            f"Removing {n_removed} pinhole wire(s) from the junction at "
            f"{tuple(node_pos)} did not leave the piece's boundary intact: "
            f"{reason}."
        )
    return cleaned, n_removed


# --- Interface resolution (master side) -------------------------------------


def _neighbour(node: NodeKey, step) -> NodeKey:
    return (node[0] + step[0], node[1] + step[1], node[2] + step[2])


def resolve_interfaces(
    lp: LatticeParams,
    interior_nodes: np.ndarray,
    pieces: list[BoundaryPiece],
) -> InterfaceSet:
    """Decide from both sides at once which caps are stitched across.

    A cap becomes an interface — a hole both sides punch, to be closed by
    stitching — only when both sides actually present material there *and* the
    two regions agree. Anything else stays as exterior surface on whichever side
    has it, which leaves that side closed rather than holed.

    An INTERIOR node presents all six of its caps whole: its half-struts are all
    further than ``r + d`` from the surface, so its cap planes lie strictly
    inside the solid (docs/algorithm.md §5.3(b)). Only trimmed boundary pieces
    can present a partial cap, or none at all.

    The failure mode this exists to remove is asymmetry, not partiality. A
    partial cap matched by an equally partial cap on the other side is a perfectly
    good interface; a whole cap facing nothing is a hole in the output.
    """
    interior_set: set[NodeKey] = {
        (int(r[0]), int(r[1]), int(r[2])) for r in interior_nodes
    }
    cap_faces_at: dict[tuple[NodeKey, int], list] = {}
    for piece in pieces:
        for key, faces in piece.cap_faces.items():
            cap_faces_at.setdefault(key, []).extend(faces)

    # Nodes that produced any material at all. A neighbour that produced nothing
    # cannot be missing a cap, so it is not evidence of anything going wrong. A
    # piece's own `node` covers the ordinary case; the `cap_faces` keys cover a
    # fused piece's other node too (docs/algorithm.md §7.1 "Fuse junction pairs").
    occupied: set[NodeKey] = (
        interior_set
        | {p.node for p in pieces}
        | {key[0] for p in pieces for key in p.cap_faces}
    )

    whole_cap = lp.t * lp.t
    area_cache: dict[tuple[NodeKey, int], float] = {}

    def cap_area(key: tuple[NodeKey, int]) -> float:
        """Total cap area at ``key``, computed only for caps that have a partner."""
        if key[0] in interior_set:
            return whole_cap
        got = area_cache.get(key)
        if got is None:
            got = sum(occ.area(f) for f in cap_faces_at[key])
            area_cache[key] = got
        return got

    present: set[tuple[NodeKey, int]] = set(cap_faces_at)
    present.update((node, h) for node in interior_set for h in range(6))

    steps = [tuple(int(x) for x in neighbor_step(h)) for h in range(6)]
    out = InterfaceSet()
    for node, h in present:
        other = _neighbour(node, steps[h])
        ho = OPPOSITE_HALF[h]
        if (other, ho) not in present:
            if other in occupied:
                out.unpaired.append((node, h))
            continue
        if h >= 3:
            continue  # each pair is decided once, from its outgoing side
        area_a = cap_area((node, h))
        area_b = cap_area((other, ho))
        if abs(area_a - area_b) > CAP_AREA_REL_TOL * max(area_a, area_b, whole_cap):
            out.mismatched.append((node, h, area_a, area_b))
            continue
        out.interfaces.add((node, h))
        out.interfaces.add((other, ho))
    return out


def finalize_pieces(pieces: list[BoundaryPiece], interfaces: set[tuple[NodeKey, int]]) -> None:
    """Settle each piece's faces now that the interfaces are known.

    A cap face is given up exactly when its cap is an interface. Every other face
    moves into ``faces`` — including a cap the other side never matched, which is
    exterior surface and closes the piece there.

    What stays in ``cap_faces`` afterwards is therefore only the interface caps,
    which are *not* part of the output: they are kept because their wires are the
    hole each weld has to recognise (:mod:`latticegen2.weld`).
    """
    for piece in pieces:
        given_up: set[tuple[NodeKey, int]] = set()
        for key, cap in list(piece.cap_faces.items()):
            if key in interfaces:
                given_up.add(key)
            else:
                piece.faces.extend(cap)
                del piece.cap_faces[key]
        piece.caps = frozenset(given_up)


# --- Fusing pieces the two booleans disagreed about --------------------------
#
# resolve_interfaces declines a cap whose two sides present regions that
# disagree, and declining is not by itself a safe degradation: where the two
# regions are only *partially* the same, keeping both caps leaves the overlap
# as non-manifold material and the remainder as an unfilled hole
# (docs/algorithm.md §7.1). The repair is to fall back to the kernel's own general
# boolean: fuse the two disagreeing pieces into one solid, which is sound
# wherever instancing's exactness argument has broken down because the kernel
# contradicted itself about a face the two share.


def _rebuild_solid(piece: BoundaryPiece) -> TopoDS_Shape:
    """Reconstruct a piece's original trimmed solid from its faces.

    Only valid *before* :func:`finalize_pieces` has run: up to that point
    ``faces`` plus every entry of ``cap_faces`` together are exactly the faces
    :func:`trim_junction` produced, nothing yet dropped or reassigned, so their
    union is the same closed boundary the boolean returned — safe to mark
    ``Closed`` without re-proving it.
    """
    all_faces = list(piece.faces)
    for group in piece.cap_faces.values():
        all_faces.extend(group)
    shell = _open_shell(all_faces)
    shell.Closed(True)
    return occ.make_solid(shell)


def _owning_cap(
    lp: LatticeParams, face, node_positions: list[tuple[NodeKey, np.ndarray]]
) -> tuple[NodeKey, int] | None:
    """Which ``(node, half-strut)`` cap ``face`` belongs to, or ``None``.

    :func:`is_cap_plane_face` tests only the one axis its half-strut id names,
    so it can pass for more than one node in the group: two nodes that share a
    coordinate along an axis orthogonal to the one separating them look
    identical to that test — every node in the same lattice "row" along ``e1``
    does, when the group is separated along ``e0``. Proximity of the face's
    centroid to each matching candidate's own ideal cap centre disambiguates,
    and is decisive: candidate centres are separated by at least ``a`` along
    some axis, while a genuine (if trimmed) cap face's centroid never strays
    from its own ideal centre by more than the profile's own extent.
    """
    best = None
    best_d = None
    centre = occ.centroid(face)
    for node, pos in node_positions:
        h = is_cap_plane_face(lp, face, pos)
        if h is None:
            continue
        ideal = pos + half_strut_offset(lp, h)
        d = float(np.linalg.norm(centre - ideal))
        if best is None or d < best_d:
            best, best_d = (node, h), d
    return best


def _fuse_group(lp: LatticeParams, group: list[BoundaryPiece]) -> BoundaryPiece:
    """Fuse one cluster of mutually disagreeing pieces into a single solid.

    Rebuilds each piece as a solid, fuses them with one ``BRepAlgoAPI_Fuse``
    call, and re-tags the result's faces with :func:`_owning_cap` against every
    node in the group — so a cap belonging to a node's *other*, non-disagreeing
    neighbours survives the fuse correctly tagged, while the disagreeing cap
    itself simply stops existing as a boundary face: it is now interior
    material shared by construction, needing no interface at all.
    """
    solids = [_rebuild_solid(p) for p in group]
    fuse = BRepAlgoAPI_Fuse()
    args = TopTools_ListOfShape()
    args.Append(solids[0])
    tools = TopTools_ListOfShape()
    for s in solids[1:]:
        tools.Append(s)
    fuse.SetArguments(args)
    fuse.SetTools(tools)
    fuse.Build()
    involved = [p.node for p in group]
    if not fuse.IsDone():
        raise ProcessingError(
            f"Could not fuse the disagreeing boundary pieces at {involved}: the "
            f"local repair boolean did not complete."
        )
    result_solids = occ.solids(fuse.Shape())
    if len(result_solids) != 1:
        raise ProcessingError(
            f"Fusing the disagreeing boundary pieces at {involved} produced "
            f"{len(result_solids)} solid(s) instead of 1. The two pieces do not "
            f"even overlap consistently, which is beyond what this repair can fix."
        )
    solid = result_solids[0]

    group_nodes = np.array([p.node for p in group], dtype=np.int64)
    node_positions = list(zip((p.node for p in group), nodes(lp, group_nodes)))
    # Re-measured on the fused result rather than inherited from the group: the
    # fuse rebuilds faces along the seam it closes, so the operands' readings
    # describe geometry that no longer exists (docs/algorithm.md §7.3).
    fused_faces = occ.faces(solid)
    tf = occ.tolerance_feature_ratio(fused_faces)
    merged = BoundaryPiece(
        node=group[0].node,
        volume=occ.volume(solid),
        tolerance_ratio=tf.ratio,
        tolerance_evidence=(tf.tolerance, tf.face_area, tf.where),
    )
    for face in fused_faces:
        tag = _owning_cap(lp, face, node_positions)
        if tag is None:
            merged.faces.append(face)
        else:
            merged.cap_faces.setdefault(tag, []).append(face)
    return merged


def fuse_disagreeing_pairs(
    lp: LatticeParams,
    pieces: list[BoundaryPiece],
    mismatched: list[tuple[NodeKey, int, float, float]],
) -> tuple[list[BoundaryPiece], int]:
    """Fuse the pieces on either side of every cap the two booleans disagreed
    about, replacing them with their fused union.

    Must run *before* :func:`finalize_pieces`, while every piece's ``faces``
    plus ``cap_faces`` still form its complete closed boundary — once a cap
    face has been sorted into "given up" or "kept" there is nothing left to
    rebuild a solid from.

    A piece can be party to more than one disagreement (docs/algorithm.md
    §7.1's rehearsal example: one node's caps disagreed with two different
    neighbours at once), so pieces are grouped by shared membership across all
    of ``mismatched`` and each connected cluster is fused in a single
    ``BRepAlgoAPI_Fuse`` call, rather than fusing pairs one at a time and
    risking the same piece being consumed twice.

    Returns the new piece list (every fused cluster's operands replaced by
    their fused result; every other piece untouched) and the number of local
    fuses performed.
    """
    if not mismatched:
        return pieces, 0

    holders: dict[tuple[NodeKey, int], list[int]] = {}
    for i, p in enumerate(pieces):
        for key in p.cap_faces:
            holders.setdefault(key, []).append(i)

    steps = [tuple(int(x) for x in neighbor_step(h)) for h in range(6)]
    uf = UnionFind(len(pieces))
    touched: set[int] = set()
    for node, h, _, _ in mismatched:
        other = _neighbour(node, steps[h])
        idxs = holders.get((node, h), []) + holders.get((other, OPPOSITE_HALF[h]), [])
        if not idxs:
            raise ProcessingError(
                f"resolve_interfaces reported a mismatched cap at {node} h{h} with "
                f"no piece holding it on either side; the two disagree about "
                f"whether it is even internally consistent."
            )
        touched.update(idxs)
        for i in idxs[1:]:
            uf.union(idxs[0], i)

    groups: dict[int, list[int]] = {}
    for i in touched:
        groups.setdefault(uf.find(i), []).append(i)

    kept = [p for i, p in enumerate(pieces) if i not in touched]
    fused = [_fuse_group(lp, [pieces[i] for i in idxs]) for idxs in groups.values()]
    return kept + fused, len(groups)


# --- Worker-process plumbing ------------------------------------------------


_PROBE_CACHE: dict = {}


def _worker_probe(mesh_path: str | None, margin: float):
    """The staged classification mesh, as an :class:`~latticegen2.classify.OutsideProbe`.

    Memoised at module scope in the worker, exactly as `classify` memoises its
    own index and for the same reason: the pool hands a worker one batch at a
    time, so without this the mesh would be re-read and its ray-casting index
    rebuilt once per batch rather than once per process.
    """
    if mesh_path is None:
        return None
    got = _PROBE_CACHE.get(mesh_path)
    if got is None:
        from .classify import OutsideProbe, load_mesh

        got = OutsideProbe(load_mesh(mesh_path), margin)
        _PROBE_CACHE[mesh_path] = got
    return got


def _worker_trim(job):
    """Trim one batch of boundary junctions in a worker process.

    Returns the path of a ``.brep`` holding the batch's geometry plus parallel
    metadata lists, so nothing but paths and small plain data crosses the
    process boundary. The ``.brep`` is a compound of per-piece compounds, each
    holding that piece's faces in the same order as the piece's cap tags.

    Below-normal priority is not part of this job tuple: it is set once per
    worker process by :class:`latticegen2.parallel.WorkerPool`'s own initializer
    rather than once per job.
    """
    (body_path, cc, t, node_batch, out_path, mesh_path, margin) = job
    lp = lattice_params(cc, t)
    tpl = build_template(lp)
    body = _read_brep(body_path)
    probe = _worker_probe(mesh_path, margin)

    node_batch = np.asarray(node_batch, dtype=np.int64)
    positions = nodes(lp, node_batch)

    bundles: list[TopoDS_Shape] = []
    meta: list[tuple[NodeKey, list[int], float, tuple]] = []
    n_empty = 0
    n_pinhole_junctions = 0
    n_pinhole_wires = 0
    n_retrimmed = 0
    n_localized = 0
    n_dropped = 0
    worst_outside = 0.0
    worst_node = None
    for i in range(len(node_batch)):
        trim = trim_junction(lp, tpl, positions[i], body, probe)
        results, n_pinholes = trim.pieces, trim.n_pinholes_removed
        if trim.retrimmed:
            n_retrimmed += 1
        if trim.localized:
            n_localized += 1
        if trim.dropped:
            n_dropped += 1
            if trim.outside_mm > worst_outside:
                worst_outside = trim.outside_mm
                worst_node = tuple(int(x) for x in node_batch[i])
        if n_pinholes:
            n_pinhole_junctions += 1
            n_pinhole_wires += n_pinholes
        if not results:
            n_empty += 1
            continue
        node = (int(node_batch[i][0]), int(node_batch[i][1]), int(node_batch[i][2]))
        for (faces, tags, vol), tf in zip(results, trim.tolerances):
            bundles.append(occ.compound(faces))
            meta.append((node, tags, vol, tuple(tf)))

    from .runlog import peak_rss_bytes

    rss = peak_rss_bytes()
    pinholes = (n_pinhole_junctions, n_pinhole_wires, n_retrimmed,
                n_localized, worst_outside, worst_node, n_dropped)
    # `rss` stays last: WorkerPool.run reads it positionally at rss_index=-1, a
    # convention shared by every worker function in this codebase.
    if bundles:
        BRepTools.Write_s(occ.compound(bundles), out_path)
        return out_path, meta, n_empty, pinholes, rss
    return None, meta, n_empty, pinholes, rss


def _piece_from(
    node: NodeKey,
    faces: list,
    tags: list[int],
    volume: float,
    tolerance: tuple[float, float, float, tuple[float, float, float]] = (
        0.0, 0.0, 0.0, (0.0, 0.0, 0.0)
    ),
) -> BoundaryPiece:
    """Split one trimmed solid's faces into plain faces and tagged cap faces.

    ``tolerance`` is :func:`latticegen2.occ.tolerance_feature_ratio`'s reading
    for this piece, flattened to plain numbers so it crosses the worker boundary
    under the same small-IPC discipline as the cap tags beside it.
    """
    if len(faces) != len(tags):
        raise ProcessingError(
            f"Boundary piece at {node} has {len(faces)} faces but {len(tags)} cap "
            f"tags; the worker result and its metadata disagree."
        )
    ratio, tol, face_area, where = tolerance
    piece = BoundaryPiece(
        node=node,
        volume=volume,
        tolerance_ratio=float(ratio),
        tolerance_evidence=(float(tol), float(face_area), tuple(where)),
    )
    for face, h in zip(faces, tags):
        if h < 0:
            piece.faces.append(face)
        else:
            piece.cap_faces.setdefault((node, int(h)), []).append(face)
    return piece


def trim_boundary(
    lp: LatticeParams,
    tpl: JunctionTemplate,
    boundary_nodes: np.ndarray,
    body: TopoDS_Shape,
    body_path: str,
    tmpdir: str,
    workers: int,
    progress=None,
    pool: WorkerPool | None = None,
    mesh_path: str | None = None,
    outside_margin: float = 0.0,
) -> BoundaryResult:
    """Trim every boundary junction, sequentially or across worker processes.

    ``pool``, if given, is used instead of building a transient one — this is
    what lets a single :class:`latticegen2.parallel.WorkerPool` serve every
    parallel stage of one run rather than each stage paying `spawn`'s
    process-creation cost on its own. When ``pool`` is ``None`` a transient
    pool is still built exactly as before, so a caller with its own worker
    count and no run-wide pool to share (a unit test, for instance) needs no
    change.
    """
    result = BoundaryResult()
    if len(boundary_nodes) == 0:
        return result

    if workers <= 1:
        # Built only on the path that uses it. On the parallel path each worker
        # builds (and memoises) its own from the same staged `.npz`, and a copy
        # here would be a mesh plus a ray-casting index the master never queries
        # -- and `_PROBE_CACHE` would then hold both alive through `stitch`,
        # `simplify` and `export`, which are where this run's peak memory is.
        probe = _worker_probe(mesh_path, outside_margin)
        positions = nodes(lp, boundary_nodes)
        for i in range(len(boundary_nodes)):
            trim = trim_junction(lp, tpl, positions[i], body, probe)
            if trim.retrimmed:
                result.n_retrimmed_junctions += 1
            if trim.localized:
                result.n_localized_junctions += 1
            if trim.dropped:
                result.n_dropped_junctions += 1
                if trim.outside_mm > result.worst_outside_mm:
                    result.worst_outside_mm = trim.outside_mm
                    result.worst_outside_node = tuple(int(x) for x in boundary_nodes[i])
            if trim.n_pinholes_removed:
                result.n_pinhole_junctions += 1
                result.n_pinhole_wires += trim.n_pinholes_removed
            if not trim.pieces:
                result.n_empty += 1
            else:
                node = tuple(int(x) for x in boundary_nodes[i])
                for (faces, tags, vol), tf in zip(trim.pieces, trim.tolerances):
                    result.pieces.append(_piece_from(node, faces, tags, vol, tuple(tf)))
            # Reported unconditionally: a junction that produced nothing is still
            # a junction processed, and skipping it would make the counter jump.
            if progress is not None:
                progress(i + 1, len(boundary_nodes))
        return result

    batches = _split_batches(boundary_nodes, workers)
    jobs = [
        (body_path, lp.cc, lp.t, nb.tolist(),
         os.path.join(tmpdir, f"boundary_{bi}.brep"), mesh_path, outside_margin)
        for bi, nb in enumerate(batches)
    ]

    # Junctions finished once the first `k` batches are done. Reporting
    # junctions rather than batches is what makes the parallel path's counter
    # mean the same thing as the sequential path's, and the cumulative sum is
    # precomputed because the tick below runs inside the dispatch loop.
    cumulative = list(itertools.accumulate(len(b) for b in batches))

    def _tick(completed: int, _jobs: int) -> None:
        # Called by `WorkerPool.run` as each batch lands, rather than from
        # `_consume` below, which runs only after the entire pool has finished:
        # reported from there, a 13-minute stage emitted its whole counter in
        # the final millisecond, which is no use to a progress bar. The sequence
        # of `(done, total)` pairs is identical either way — this stage passes
        # no `sort_by`, so ordered `imap` yields batch `k` only once batches
        # `0..k` are all complete, which is the same arithmetic `_consume` did.
        if progress is not None:
            progress(cumulative[completed - 1], len(boundary_nodes))

    def _consume(results) -> None:
        # Ordered results, not completion order: `result.pieces` — and therefore
        # the shape list handed to sewing — is identical run to run. Sewing
        # resolves near-coincident vertices in the order it receives them, so an
        # arbitrary completion order would let two runs of the same command
        # produce byte-different output.
        for path, meta, n_empty, pinholes, rss in results:
            result.n_pinhole_junctions += pinholes[0]
            result.n_pinhole_wires += pinholes[1]
            result.n_retrimmed_junctions += pinholes[2]
            result.n_localized_junctions += pinholes[3]
            result.n_dropped_junctions += pinholes[6]
            if pinholes[4] > result.worst_outside_mm:
                result.worst_outside_mm = pinholes[4]
                result.worst_outside_node = pinholes[5]
            result.n_empty += n_empty
            result.max_worker_rss = max(result.max_worker_rss, rss)
            if path is None:
                continue
            shape = _read_brep(path)
            children = _compound_children(shape)
            if len(children) != len(meta):
                raise ProcessingError(
                    f"Boundary worker result mismatch in {path}: {len(children)} pieces "
                    f"vs {len(meta)} metadata records."
                )
            for bundle, (node, tags, vol, tf) in zip(children, meta):
                result.pieces.append(
                    _piece_from(
                        tuple(node), _compound_children(bundle), tags, vol, tuple(tf)
                    )
                )

    if pool is not None:
        results, _ = pool.run(_worker_trim, jobs, on_result=_tick)
        _consume(results)
    else:
        with WorkerPool(workers) as owned:
            results, _ = owned.run(_worker_trim, jobs, on_result=_tick)
        _consume(results)
    return result


def _split_batches(node_index: np.ndarray, workers: int):
    """Split jobs into a few batches per worker, keeping spatial locality.

    More batches than workers keeps the pool busy when junction costs vary;
    keeping each batch contiguous in the node ordering means a worker's
    junctions share input-body regions, which helps OCCT's own caching.
    """
    n = len(node_index)
    target = max(1, min(n, workers * 4))
    bounds = np.linspace(0, n, target + 1).astype(int)
    return [
        node_index[bounds[i]:bounds[i + 1]]
        for i in range(target)
        if bounds[i + 1] > bounds[i]
    ]
