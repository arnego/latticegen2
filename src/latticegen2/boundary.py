"""Boundary junctions — the only place a boolean touches the input body.

Each BOUNDARY node contributes exactly one intersection: its instanced junction
solid against the input solid. That is deliberately **one object operand per
call**. OCCT's general boolean runs over all object operands together, so two
overlapping objects in one call are *partitioned* against the tool rather than
each trimmed independently — the fragmentation that cost the Julia pipeline
dearly (docs/algorithm.md §6.3). A single already-fused junction cannot trigger
it, by construction.

After trimming, every face lying in an interface cap plane is dropped, exactly
as the interior path drops caps, so the trimmed junction presents the same
square hole to its neighbour and the two stitch together. Which caps count as
interfaces is decided by classification, not by inspecting geometry: a cap is an
interface iff the node across it is itself kept.

The work is embarrassingly parallel — constant-size, independent jobs — so it is
distributed over worker processes with the same small-IPC discipline the Julia
implementation used: the input body goes to disk once as a ``.brep`` and workers
read it directly; only file paths and small metadata cross the process boundary.
"""

from __future__ import annotations

import multiprocessing as mp
import os
from dataclasses import dataclass, field

import numpy as np
from OCP.BRep import BRep_Builder
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common
from OCP.BRepTools import BRepTools
from OCP.TopoDS import TopoDS_Shape, TopoDS_Shell
from OCP.TopTools import TopTools_ListOfShape

from . import occ
from .errors import ProcessingError
from .junction import JunctionTemplate, build_template, is_cap_plane_face
from .lattice import LatticeParams, lattice_params, neighbor_step, nodes


@dataclass
class BoundaryPiece:
    """One connected solid produced by trimming one boundary junction."""

    node: tuple[int, int, int]
    shell: TopoDS_Shape
    """The piece's boundary with its interface caps removed, ready to sew."""
    caps: frozenset[int]
    """Half-strut ids whose cap face this piece carried and gave up as an interface."""
    volume: float


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
    interface_mask: int,
):
    """Intersect one instanced junction with the body and drop its interface caps.

    ``interface_mask`` has bit ``h`` set when the node across half-strut ``h`` is
    itself kept, i.e. when that cap is an interface rather than exterior surface.

    Returns a list of ``(shell, caps)`` pairs, one per connected solid the
    intersection produced. A boundary junction can legitimately split into
    several pieces when the input surface cuts between its arms; each piece
    becomes its own vertex in the connectivity graph.
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
        keep_faces = []
        dropped: set[int] = set()
        for face in occ.faces(solid):
            h = is_cap_plane_face(lp, face, node_pos)
            if h is not None and (interface_mask >> h) & 1:
                dropped.add(h)
                continue
            keep_faces.append(face)
        out.append((_open_shell(keep_faces), frozenset(dropped), occ.volume(solid)))
    return out


def interface_masks(node_index: np.ndarray, kept: set[tuple[int, int, int]]) -> np.ndarray:
    """Per-node 6-bit masks marking which caps face another kept node."""
    steps = [neighbor_step(h) for h in range(6)]
    masks = np.zeros(len(node_index), dtype=np.int64)
    for i in range(len(node_index)):
        n = node_index[i]
        m = 0
        for h in range(6):
            nb = (int(n[0] + steps[h][0]), int(n[1] + steps[h][1]), int(n[2] + steps[h][2]))
            if nb in kept:
                m |= 1 << h
        masks[i] = m
    return masks


# --- Worker-process plumbing ------------------------------------------------


def _worker_trim(job):
    """Trim one batch of boundary junctions in a worker process.

    Returns the path of a ``.brep`` holding the batch's shells plus parallel
    metadata lists, so nothing but paths and small plain data crosses the
    process boundary.
    """
    (body_path, cc, t, node_batch, mask_batch, out_path, background) = job
    if background:
        _set_background_priority()
    lp = lattice_params(cc, t)
    tpl = build_template(lp)
    body = _read_brep(body_path)

    node_batch = np.asarray(node_batch, dtype=np.int64)
    positions = nodes(lp, node_batch)

    shells: list[TopoDS_Shape] = []
    meta: list[tuple[tuple[int, int, int], list[int], float]] = []
    n_empty = 0
    for i in range(len(node_batch)):
        results = trim_junction(lp, tpl, positions[i], body, int(mask_batch[i]))
        if not results:
            n_empty += 1
            continue
        node = (int(node_batch[i][0]), int(node_batch[i][1]), int(node_batch[i][2]))
        for shell, caps, vol in results:
            shells.append(shell)
            meta.append((node, sorted(caps), vol))

    from .runlog import peak_rss_bytes

    rss = peak_rss_bytes()
    if shells:
        BRepTools.Write_s(occ.compound(shells), out_path)
        return out_path, meta, n_empty, rss
    return None, meta, n_empty, rss


def _read_brep(path: str) -> TopoDS_Shape:
    shape = TopoDS_Shape()
    builder = BRep_Builder()
    if not BRepTools.Read_s(shape, path, builder):
        raise ProcessingError(f"Could not read intermediate BREP file: {path}")
    return shape


def _set_background_priority() -> None:
    """Drop this process to below-normal scheduling priority (``-bg``)."""
    try:
        if os.name == "nt":
            import ctypes

            BELOW_NORMAL_PRIORITY_CLASS = 0x00004000
            ctypes.windll.kernel32.SetPriorityClass(
                ctypes.windll.kernel32.GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS
            )
        else:
            os.nice(5)
    except Exception:
        pass  # priority is a courtesy, never a reason to fail a run


def trim_boundary(
    lp: LatticeParams,
    tpl: JunctionTemplate,
    boundary_nodes: np.ndarray,
    kept: set[tuple[int, int, int]],
    body: TopoDS_Shape,
    body_path: str,
    tmpdir: str,
    workers: int,
    background: bool = False,
    progress=None,
) -> BoundaryResult:
    """Trim every boundary junction, sequentially or across worker processes."""
    result = BoundaryResult()
    if len(boundary_nodes) == 0:
        return result

    masks = interface_masks(boundary_nodes, kept)

    if workers <= 1:
        positions = nodes(lp, boundary_nodes)
        for i in range(len(boundary_nodes)):
            pieces = trim_junction(lp, tpl, positions[i], body, int(masks[i]))
            if not pieces:
                result.n_empty += 1
                continue
            node = tuple(int(x) for x in boundary_nodes[i])
            for shell, caps, vol in pieces:
                result.pieces.append(BoundaryPiece(node=node, shell=shell, caps=caps, volume=vol))
            if progress is not None:
                progress(i + 1, len(boundary_nodes))
        return result

    batches = _split_batches(boundary_nodes, masks, workers)
    jobs = [
        (
            body_path,
            lp.cc,
            lp.t,
            nb.tolist(),
            mb.tolist(),
            os.path.join(tmpdir, f"boundary_{bi}.brep"),
            background,
        )
        for bi, (nb, mb) in enumerate(batches)
    ]

    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=workers) as pool:
        done = 0
        for path, meta, n_empty, rss in pool.imap_unordered(_worker_trim, jobs):
            result.n_empty += n_empty
            result.max_worker_rss = max(result.max_worker_rss, rss)
            done += 1
            if progress is not None:
                progress(done, len(jobs))
            if path is None:
                continue
            shape = _read_brep(path)
            children = _compound_children(shape)
            if len(children) != len(meta):
                raise ProcessingError(
                    f"Boundary worker result mismatch in {path}: {len(children)} shells "
                    f"vs {len(meta)} metadata records."
                )
            for shell, (node, caps, vol) in zip(children, meta):
                result.pieces.append(
                    BoundaryPiece(node=tuple(node), shell=shell, caps=frozenset(caps), volume=vol)
                )
    return result


def _compound_children(shape: TopoDS_Shape) -> list[TopoDS_Shape]:
    """Direct children of a compound, in the order they were added."""
    from OCP.TopoDS import TopoDS_Iterator

    out = []
    it = TopoDS_Iterator(shape)
    while it.More():
        out.append(it.Value())
        it.Next()
    return out


def _split_batches(node_index: np.ndarray, masks: np.ndarray, workers: int):
    """Split jobs into a few batches per worker, keeping spatial locality.

    More batches than workers keeps the pool busy when junction costs vary;
    keeping each batch contiguous in the node ordering means a worker's
    junctions share input-body regions, which helps OCCT's own caching.
    """
    n = len(node_index)
    target = max(1, min(n, workers * 4))
    bounds = np.linspace(0, n, target + 1).astype(int)
    return [
        (node_index[bounds[i]:bounds[i + 1]], masks[bounds[i]:bounds[i + 1]])
        for i in range(target)
        if bounds[i + 1] > bounds[i]
    ]
