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

import os
from dataclasses import dataclass, field

import numpy as np
from OCP.BRep import BRep_Builder
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Fuse
from OCP.BRepTools import BRepTools
from OCP.TopoDS import TopoDS_Shape, TopoDS_Shell
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
    max_worker_rss: int = 0
    """Highest peak RSS reported by any worker process.

    Reported separately from the master's own peak because they are different
    processes: the run's true memory footprint is the master's peak plus what
    the workers held concurrently, and specification.md §3 asks for maximum
    memory usage.
    """
    diagnostics: list[str] = field(default_factory=list)


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


def trim_junction(
    lp: LatticeParams,
    tpl: JunctionTemplate,
    node_pos: np.ndarray,
    body: TopoDS_Shape,
):
    """Intersect one instanced junction with the body and tag its cap faces.

    Returns a list of ``(faces, tags, volume)`` triples, one per connected solid
    the intersection produced. ``tags[i]`` is the half-strut id of the cap plane
    ``faces[i]`` lies in, or ``-1`` for every other face. A boundary junction can
    legitimately split into several pieces when the input surface cuts between
    its arms; each piece becomes its own vertex in the connectivity graph.

    Caps are tagged rather than dropped because whether a cap is an interface is
    not a local fact — it depends on what the junction on the other side of that
    cap produced, which is a different boolean in a different process. See the
    module docstring and :func:`resolve_interfaces`.
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

    out = []
    for solid in occ.solids(algo.Shape()):
        faces = occ.faces(solid)
        tags = [is_cap_plane_face(lp, f, node_pos) for f in faces]
        tags = [-1 if h is None else h for h in tags]
        out.append((faces, tags, occ.volume(solid)))
    return out


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
    merged = BoundaryPiece(node=group[0].node, volume=occ.volume(solid))
    for face in occ.faces(solid):
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


def _worker_trim(job):
    """Trim one batch of boundary junctions in a worker process.

    Returns the path of a ``.brep`` holding the batch's geometry plus parallel
    metadata lists, so nothing but paths and small plain data crosses the
    process boundary. The ``.brep`` is a compound of per-piece compounds, each
    holding that piece's faces in the same order as the piece's cap tags.

    ``background`` priority is no longer part of this job tuple: it is set
    once per worker process by :class:`latticegen2.parallel.WorkerPool`'s own
    initializer rather than once per job.
    """
    (body_path, cc, t, node_batch, out_path) = job
    lp = lattice_params(cc, t)
    tpl = build_template(lp)
    body = _read_brep(body_path)

    node_batch = np.asarray(node_batch, dtype=np.int64)
    positions = nodes(lp, node_batch)

    bundles: list[TopoDS_Shape] = []
    meta: list[tuple[NodeKey, list[int], float]] = []
    n_empty = 0
    for i in range(len(node_batch)):
        results = trim_junction(lp, tpl, positions[i], body)
        if not results:
            n_empty += 1
            continue
        node = (int(node_batch[i][0]), int(node_batch[i][1]), int(node_batch[i][2]))
        for faces, tags, vol in results:
            bundles.append(occ.compound(faces))
            meta.append((node, tags, vol))

    from .runlog import peak_rss_bytes

    rss = peak_rss_bytes()
    if bundles:
        BRepTools.Write_s(occ.compound(bundles), out_path)
        return out_path, meta, n_empty, rss
    return None, meta, n_empty, rss


def _piece_from(node: NodeKey, faces: list, tags: list[int], volume: float) -> BoundaryPiece:
    """Split one trimmed solid's faces into plain faces and tagged cap faces."""
    if len(faces) != len(tags):
        raise ProcessingError(
            f"Boundary piece at {node} has {len(faces)} faces but {len(tags)} cap "
            f"tags; the worker result and its metadata disagree."
        )
    piece = BoundaryPiece(node=node, volume=volume)
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
    background: bool = False,
    progress=None,
    pool: WorkerPool | None = None,
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
        positions = nodes(lp, boundary_nodes)
        for i in range(len(boundary_nodes)):
            trimmed = trim_junction(lp, tpl, positions[i], body)
            if not trimmed:
                result.n_empty += 1
            else:
                node = tuple(int(x) for x in boundary_nodes[i])
                for faces, tags, vol in trimmed:
                    result.pieces.append(_piece_from(node, faces, tags, vol))
            # Reported unconditionally: a junction that produced nothing is still
            # a junction processed, and skipping it would make the counter jump.
            if progress is not None:
                progress(i + 1, len(boundary_nodes))
        return result

    batches = _split_batches(boundary_nodes, workers)
    jobs = [
        (body_path, lp.cc, lp.t, nb.tolist(), os.path.join(tmpdir, f"boundary_{bi}.brep"))
        for bi, nb in enumerate(batches)
    ]

    def _consume(results) -> None:
        done = 0
        # Ordered results, not completion order: `result.pieces` — and therefore
        # the shape list handed to sewing — is identical run to run. Sewing
        # resolves near-coincident vertices in the order it receives them, so an
        # arbitrary completion order would let two runs of the same command
        # produce byte-different output.
        for batch_index, (path, meta, n_empty, rss) in enumerate(results):
            result.n_empty += n_empty
            result.max_worker_rss = max(result.max_worker_rss, rss)
            # Report junctions, the same unit the sequential path reports, rather
            # than batches — results are ordered, so batch i is the i'th job.
            done += len(batches[batch_index])
            if progress is not None:
                progress(done, len(boundary_nodes))
            if path is None:
                continue
            shape = _read_brep(path)
            children = _compound_children(shape)
            if len(children) != len(meta):
                raise ProcessingError(
                    f"Boundary worker result mismatch in {path}: {len(children)} pieces "
                    f"vs {len(meta)} metadata records."
                )
            for bundle, (node, tags, vol) in zip(children, meta):
                result.pieces.append(
                    _piece_from(tuple(node), _compound_children(bundle), tags, vol)
                )

    if pool is not None:
        results, _ = pool.run(_worker_trim, jobs)
        _consume(results)
    else:
        with WorkerPool(workers, background) as owned:
            results, _ = owned.run(_worker_trim, jobs)
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
