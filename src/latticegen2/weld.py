"""Assembly: sew the boundary layer, then build the interior onto it.

The problem this solves is measured in `tools/prototypes/RESULTS.md` G5.
`BRepBuilderAPI_Sewing` costs about `n^1.8` in piece count, no combination of its
optional phases changes that by more than 2 %, and — decisively — adding one
**closed** 194,400-face shell with zero free edges to a 4,000-piece sew took it
from 76.5 s to 716.6 s. Total face count dominates. The interior shell's face
count scales with the *volume* of the part, so handing it to sewing put a
superlinear volume term on the pipeline: 4 h 45 m of a 5 h 04 m run.

So the interior shell never goes near sewing. The assembly runs the other way
round:

1. **Sew the boundary pieces to each other**, per junction-graph component.
   That is the surface-area-scaling part, and it is the only place a pairing is
   genuinely unknown: both sides of a boundary-to-boundary cap come out of
   independent booleans with no shared topology to exploit.
2. **Read the interface rings off the sewn result.** Its free edges are exactly
   the holes facing the interior, and each one is the whole template cap quad
   (docs/algorithm.md §5.3(b)), so they can be identified by their corners
   without any search beyond a dictionary lookup.
3. **Build the interior shell on those rings.** The instancing index adopts the
   boundary's vertices and edges at every interface instead of making its own,
   so the two are the same objects from the start.
4. **Assemble with `BRep_Builder.Add`** and prove watertightness: every edge used
   exactly twice, once each way.

**Why the interior is the side that adapts.** The obvious alternative is to
rewrite the boundary pieces to use the interior's topology with
`BRepTools_ReShape`. That does not work, and the way it fails is quiet.
`ReShape` will swap an edge inside a face happily, but replacing that edge's
**vertices** leaves the neighbouring edges still pointing at the old ones and the
wire comes apart: measured on the 80 mm ball, `BRepCheck_NotConnected` wires, the
solid invalid, and the volume wrong — 29,111 mm³ against a true 51,393 mm³ —
while every edge still had exactly two faces and the shell still "closed".
Keeping the vertices makes the same swap perfectly clean. Since a cap's two sides
cannot each keep their own vertices, the side that gives way has to be the one
whose faces we build ourselves.

**Step 1 is itself tiled**, per component, once a component is large enough for
it to matter (docs/algorithm.md §8). G5a
(`tools/prototypes/RESULTS.md`) measured sewing at about `n^1.8` in piece count
with no configuration that changes that by more than 2 % — the cost is a search,
and the pairing it searches for is already known exactly from the junction graph
(§7.1). Splitting a component's pieces into spatial tiles by lattice-index
block, sewing each tile on its own — across worker processes when there are
enough of them — and then sewing the tile results applies that superlinear term
to a much smaller `n` in round 1. Round 2 is a real but smaller saving, not a
free one: G6 measured it staying close to a monolithic sew's own cost rather
than shrinking with tile count, because it still sews shells whose *combined*
face count equals the untiled input's — 1.3–1.45× overall, not the
order-of-magnitude round 1 alone would suggest. Either way, tiling only changes
*how* the component is sewn, never *what* it sews: every piece still ends up in
exactly one final sew call's input, so the result is the same shell regardless
of how it was partitioned to get there — a wall-clock lever, not a geometry
decision, matching docs/algorithm.md §11's rule that an optimization's failure
mode must be "do more work", never "produce a different result".
"""

from __future__ import annotations

import os
import time
from dataclasses import dataclass, field

import numpy as np
from OCP.BRep import BRep_Builder, BRep_Tool
from OCP.TopAbs import TopAbs_Orientation, TopAbs_ShapeEnum
from OCP.TopExp import TopExp
from OCP.TopoDS import TopoDS, TopoDS_Shell
from OCP.TopTools import TopTools_IndexedDataMapOfShapeListOfShape, TopTools_IndexedMapOfShape

from . import occ
from .errors import ProcessingError
from .lattice import OPPOSITE_HALF, LatticeParams, neighbor_step
from .parallel import WorkerPool
from .parallel import compound_children as _compound_children
from .parallel import read_brep as _read_brep
from .parallel import write_brep as _write_brep

NodeKey = tuple[int, int, int]

WELD_TOL = 1e-6
"""Millimetres. How far apart two sides of an interface may place the same corner
and still be recognised as the same corner.

A *recognition* bar, not a repair: what follows is adoption, so a corner that
matches becomes the very same object rather than two objects a tolerance apart.
Measured, the two sides land on top of each other — an untouched cap comes back
from `BRepAlgoAPI_Common` at the template's own coordinates, to 0.0 mm. The bar
stays far below the smallest real feature (``t >= 0.4`` mm), so it cannot fuse two
genuinely distinct corners."""

SEW_TOLERANCE = 1e-6
"""Millimetres, for the boundary-to-boundary sew. Both sides of such a cap are
the same nominal region computed by two independent booleans, so this only has to
absorb their disagreement, which is far below it.

Raising this does **not** fix the micron-scale slivers a grazing trim leaves
(specification.md §10). Tried at 1e-5 on the `TD_HX_rehearsal_test` rehearsal, which
carries two such edges of 3.171690e-06 mm and 5.808982e-06 mm: both survived
unchanged. The tolerance governs whether two *different* faces' free edges are
paired up, so it cannot remove an edge that has no partner to be paired with —
which is what a sliver on one side of a near-tangential trim is."""


# --- rings ------------------------------------------------------------------


@dataclass
class Ring:
    """The boundary of one cap hole: its edges, as endpoint-position pairs.

    Deliberately unordered and wire-agnostic. Two independently computed versions
    of the same trimmed region can traverse it in either direction and can split
    it across several wires, so matching them by an ordered outer wire would
    reject sound geometry. Matching edge by edge on endpoints does not.
    """

    edges: list = field(default_factory=list)
    verts: list = field(default_factory=list)
    """Per edge, its two ``TopoDS_Vertex`` objects."""
    points: list = field(default_factory=list)
    """Per edge, its two endpoint positions as a ``(2, 3)`` array."""

    def __len__(self) -> int:
        return len(self.points)


def _pnt(vertex) -> np.ndarray:
    p = BRep_Tool.Pnt_s(TopoDS.Vertex_s(vertex))
    return np.array([p.X(), p.Y(), p.Z()])


def ring_of_edges(edges) -> Ring:
    """A ring from edges, each reported in **its own** natural direction.

    Every edge is normalised to FORWARD before its endpoints are read, so that
    "which way round" means the same thing on both sides of a comparison: an edge
    explored out of a face otherwise carries that face's orientation composed
    into it.
    """
    ring = Ring()
    for edge in edges:
        fwd = TopoDS.Edge_s(edge).Oriented(TopAbs_Orientation.TopAbs_FORWARD)
        vs = [TopoDS.Vertex_s(v) for v in occ._explore(fwd, TopAbs_ShapeEnum.TopAbs_VERTEX)]
        if len(vs) != 2:
            # A closed or degenerate edge has no two distinct corners to match on.
            return Ring(edges=[None], verts=[None], points=[None])
        ring.edges.append(fwd)
        ring.verts.append(vs)
        ring.points.append(np.array([_pnt(vs[0]), _pnt(vs[1])]))
    return ring


def ring_of_face(face) -> Ring:
    """Every edge bounding ``face``, with its endpoints."""
    return ring_of_edges(list(occ._explore(face, TopAbs_ShapeEnum.TopAbs_EDGE)))


def ring_of_faces(faces) -> Ring:
    """The union of several faces' edges — one cap can arrive as fragments."""
    out = Ring()
    for face in faces:
        part = ring_of_face(face)
        out.edges.extend(part.edges)
        out.verts.extend(part.verts)
        out.points.extend(part.points)
    return out


def template_cap_corners(lp: LatticeParams, tmesh, node: NodeKey, h: int) -> np.ndarray:
    """The cap quad's corners at ``(node, h)``, in the template's loop order."""
    fi = next(i for i, c in enumerate(tmesh.face_cap) if c == h)
    base = lp.B @ np.array(node, dtype=float)
    return np.array([base + tmesh.verts[local] for local in tmesh.loops[fi]])


def template_cap_ring(lp: LatticeParams, tmesh, node: NodeKey, h: int) -> Ring:
    """The whole cap quad at ``(node, h)``, from the lattice expressions alone.

    Available before any geometry is assembled, which is what lets an interface
    whose two sides do not correspond be rejected *before* either gives up its
    cap face.
    """
    corners = template_cap_corners(lp, tmesh, node, h)
    ring = Ring()
    for i in range(len(corners)):
        ring.edges.append(None)
        ring.verts.append(None)
        ring.points.append(np.array([corners[i], corners[(i + 1) % len(corners)]]))
    return ring


def _same_edge(a: np.ndarray, b: np.ndarray, tol: float) -> bool:
    """Do two endpoint pairs describe the same segment, in either direction?"""
    forward = max(np.linalg.norm(a[0] - b[0]), np.linalg.norm(a[1] - b[1]))
    backward = max(np.linalg.norm(a[0] - b[1]), np.linalg.norm(a[1] - b[0]))
    return min(forward, backward) <= tol


def match_rings(a: Ring, b: Ring, tol: float = WELD_TOL) -> list[int] | None:
    """``b``'s edges in terms of ``a``'s, or ``None`` if they do not correspond.

    Requires a bijection. Ambiguity — two candidates within tolerance — is
    rejected rather than resolved, because picking the wrong one would be silent.
    """
    if len(a) != len(b) or len(a) == 0:
        return None
    if any(p is None for p in a.points) or any(p is None for p in b.points):
        return None
    taken = [False] * len(a)
    out = []
    for pb in b.points:
        found = [i for i, pa in enumerate(a.points) if not taken[i] and _same_edge(pa, pb, tol)]
        if len(found) != 1:
            return None
        taken[found[0]] = True
        out.append(found[0])
    return out


# --- planning ---------------------------------------------------------------


def _neighbour(node: NodeKey, step) -> NodeKey:
    return (node[0] + step[0], node[1] + step[1], node[2] + step[2])


def unweldable(
    lp: LatticeParams,
    tmesh,
    interior_set: set[NodeKey],
    pieces,
    interfaces: set[tuple[NodeKey, int]],
) -> list[tuple[NodeKey, int]]:
    """Interior interfaces whose boundary side is not the whole cap quad.

    Only interior-to-boundary caps are checked. Those are the ones the interior
    will be built onto, so their rings have to be exactly the template quad —
    §5.3(b) says they always are, and this is where that stops being an
    assumption. Boundary-to-boundary caps are sewn rather than adopted, so they
    need no correspondence of their own; interior-to-interior caps share one
    index already.

    Run *before* the caps are dropped, so a rejection costs an extra solid rather
    than a hole.
    """
    caps_at: dict[tuple[NodeKey, int], list] = {}
    for piece in pieces:
        for key, faces in piece.cap_faces.items():
            caps_at.setdefault(key, []).extend(faces)

    steps = [tuple(int(x) for x in neighbor_step(h)) for h in range(6)]
    rejected = []
    for node, h in interfaces:
        if node not in interior_set:
            continue
        other = (_neighbour(node, steps[h]), OPPOSITE_HALF[h])
        if other[0] in interior_set:
            continue
        expected = template_cap_ring(lp, tmesh, node, h)
        if match_rings(expected, ring_of_faces(caps_at.get(other, []))) is None:
            rejected.append((node, h))
    return rejected


# --- step 1: sew the boundary layer -----------------------------------------


TILE_TARGET_PIECES = 500
"""Pieces aimed for per boundary-sew tile, by choosing the tile's index-space
edge length.

Tied to measurement rather than picked blind (`tools/prototypes/RESULTS.md`):
G5a clocks 500 pieces at 1.46 s, well before the `n^1.8` term costs anything,
and G6 — which actually tiles and re-sews — finds the total cost bottoms out
around a few hundred to ~1,000 pieces per tile at both 4,000 and 8,000 pieces
tried, because round 2 (sewing the tiles' results together) does not shrink
with finer tiling the way round 1 does. 500 sits inside that measured plateau
rather than at either edge of it."""

MIN_PIECES_TO_TILE = 3 * TILE_TARGET_PIECES
"""Below this a component sews faster in one call than tiling could ever save.

Tiling adds a second sewing round, plus a worker round-trip when parallel, so it
only pays for itself once a component is big enough to split into several tiles
in the first place — one tile is definitionally not worth tiling
(:func:`_tile_pieces` returns ``None`` for it regardless). Set at three tiles'
worth so tiling never engages for a component that would only ever produce one
or two, where the second sewing round is pure overhead."""


def _tile_edge_length(node_index: np.ndarray, target: int) -> int:
    """Index-space edge length whose cube holds about ``target`` pieces on average.

    A density estimate, not a promise: the boundary layer is a thin shell in
    index space, so occupied cells are far sparser than the bounding box this
    divides implies, and any one real tile will hold more or fewer than
    ``target``. That only ever costs or saves wall time — every tiling, balanced
    or not, sews to the same result (module docstring), so a rough estimate is
    all this needs to be.
    """
    lo = node_index.min(axis=0)
    hi = node_index.max(axis=0)
    extent = (hi - lo + 1).astype(np.float64)
    n = len(node_index)
    v_idx = float(np.prod(extent))
    if n <= target or v_idx <= 0:
        return int(max(extent.max(), 1))  # one tile, covering everything
    edge = (target * v_idx / n) ** (1.0 / 3.0)
    return max(1, int(round(edge)))


def _tile_pieces(pieces: list, target: int, min_to_tile: int) -> list[list] | None:
    """Bucket ``pieces`` into spatial tiles by lattice-index block.

    Returns ``None`` when tiling would not be worth it — too few pieces, or a
    footprint small enough that every piece lands in one block anyway — so the
    caller falls back to a single sew, unchanged from before this existed.

    Bucketing is by each piece's own node (or, for a piece
    :func:`latticegen2.boundary.fuse_disagreeing_pairs` has merged, its
    representative node — close enough for a spatial *tile assignment*, where
    :func:`latticegen2.boundary._owning_cap`'s exactness is not needed). Iteration
    order is ``pieces``' own order throughout, so which tile gets which pieces —
    and therefore the order tiles are sewn in — is identical run to run.
    """
    if len(pieces) < min_to_tile:
        return None
    node_index = np.array([p.node for p in pieces], dtype=np.int64)
    lo = node_index.min(axis=0)
    edge = _tile_edge_length(node_index, target)
    buckets: dict[tuple, list] = {}
    for piece, node in zip(pieces, node_index):
        key = tuple(int((node[d] - lo[d]) // edge) for d in range(3))
        buckets.setdefault(key, []).append(piece)
    if len(buckets) <= 1:
        return None
    return list(buckets.values())


def _sew_faces(face_lists, tolerance: float) -> list:
    """Sew several already-computed face lists together into one."""
    shells = [occ.faces_shell(faces) for faces in face_lists]
    return occ.faces(occ.sew(shells, tolerance, cutting=False))


def _worker_sew_tile(job):
    """Sew one bundle of already-computed face lists in a worker process.

    Shared by both sewing rounds: round 1's tiles and round 2's per-component
    merges are both "sew a handful of face-list bundles together", so one
    worker body does both. Mirrors
    :func:`latticegen2.boundary._worker_trim`'s small-IPC discipline: only a
    file path and a peak-RSS number cross the process boundary, never
    geometry. The input ``.brep`` is a compound of per-bundle compounds.

    Below-normal priority is not part of this job tuple: it is set once per
    worker process by :class:`latticegen2.parallel.WorkerPool`'s own initializer
    rather than once per job.
    """
    in_path, out_path, tolerance = job
    shape = _read_brep(in_path)
    bundle_faces = [occ.faces(pc) for pc in _compound_children(shape)]
    result_faces = _sew_faces(bundle_faces, tolerance)
    _write_brep(occ.compound(result_faces), out_path)

    from .runlog import peak_rss_bytes

    return out_path, peak_rss_bytes()


def _sew_all_tiles(
    plan: dict[int, list[list] | None],
    tolerance: float,
    workers: int,
    tmpdir: str | None,
    pool: WorkerPool | None = None,
) -> tuple[dict[int, list[list]], int]:
    """Round 1 for every tiled component in ``plan``, in **one** worker pool.

    ``plan`` maps each component to its tiles (from :func:`_tile_pieces`) or to
    ``None`` for a component that is not being tiled — those are skipped here,
    ``sew_boundary`` sews them directly. One pool for the whole call, not one per
    component, is the literal reading of docs/algorithm.md §8 ("in parallel
    across the run's worker processes") and matters in practice: a part with
    several large tiled components — the `TD_HX_rehearsal_test` rehearsal has 14 — would
    otherwise pay `spawn`'s process-creation cost once per component instead of
    once for the whole stitch stage.

    ``pool``, if given, is the run's shared :class:`latticegen2.parallel.WorkerPool`
    (docs/algorithm.md §8) rather than one built and torn down for this call
    alone. When ``pool`` is ``None`` a transient one is still built, so a caller
    with its own worker count and no run-wide pool to share needs no change.

    Sequential, in-process when there is no worker pool to use, or too few jobs
    to justify one — still a real saving over one monolithic sew per component,
    since it is smaller sews that are cheap, not parallel ones.

    Jobs are consumed in job order, not completion order: which order round 2
    receives each component's tile results in is identical run to run
    (:func:`latticegen2.boundary.trim_boundary` relies on the same property for
    the same reason).

    Filenames carry the component id as well as the tile index — reused index-only
    names across components would never collide at runtime (`sew_boundary`
    finishes writing and reading back one component's jobs before the next
    component's are built), only in the temp folder a failed run leaves behind,
    which defeats the post-mortem analysis it exists for (specification.md §4.4).
    """
    jobs_meta = [
        (group, i, tile)
        for group, tiles in plan.items() if tiles is not None
        for i, tile in enumerate(tiles)
    ]
    results: dict[int, list[list]] = {
        group: [None] * len(tiles) for group, tiles in plan.items() if tiles is not None
    }
    if not jobs_meta:
        return results, 0

    if (pool is None and workers <= 1) or tmpdir is None or len(jobs_meta) < 2:
        for group, i, tile in jobs_meta:
            results[group][i] = _sew_faces([p.faces for p in tile], tolerance)
        return results, 0

    jobs = []
    for group, i, tile in jobs_meta:
        in_path = os.path.join(tmpdir, f"sew_tile_{group}_{i}.brep")
        out_path = os.path.join(tmpdir, f"sew_tile_{group}_{i}_out.brep")
        _write_brep(occ.compound(occ.compound(p.faces) for p in tile), in_path)
        jobs.append((in_path, out_path, tolerance))

    def _run(p: WorkerPool):
        raw, max_rss = p.run(_worker_sew_tile, jobs)
        for (group, i, _tile), (out_path, _rss) in zip(jobs_meta, raw):
            results[group][i] = occ.faces(_read_brep(out_path))
        return results, max_rss

    if pool is not None and pool.active:
        return _run(pool)
    with WorkerPool(min(workers, len(jobs))) as owned:
        return _run(owned)


def _split_seam_interior(face_lists: list[list]) -> tuple[list[list], list]:
    """Split each tile's round-1 result into faces round 2 can still affect.

    After round 1, every tile's own result is itself a sewn shell: each of its
    edges is either used **twice** already (joined to a neighbour within the
    same tile — nothing left for round 2 to do with the face that owns it) or
    used **once** — a genuine free edge, which is either a tile-to-tile seam
    round 2 exists to close, or an interface hole meant to stay open for the
    interior shell, and round 2 cannot (and does not need to) tell those two
    apart in advance any more than a full round 2 call does. So the only faces
    round 2 can possibly affect are the ones bearing at least one free edge;
    everything else — the ``interior`` return value — is carried into the
    final shell unchanged, by direct reference, with no sewing call ever
    touching it.

    G8 (`tools/prototypes/RESULTS.md`) confirmed this by identity, not just by
    argument: sewing only the free-edge-bearing subset and concatenating the
    rest by reference reproduces a full round 2 exactly (same face count, same
    free-edge count, same volume to machine precision) at two scales tried,
    and the seam fraction measured there — 13–14 % of a tile's faces — is what
    makes this worth doing: round 2 stops paying its flat per-face cost
    (docs/algorithm.md §8) for the 86–87 % of faces it could never touch anyway.

    The split itself must stay `O(faces)`, not `O(faces²)`: an earlier version
    called :func:`free_edges` for a plain Python list and then tested every
    face's every edge against it with `.IsSame()` — each test `O(len(fe))` and
    ``TopoDS_Shape`` not being cheaply hashable in plain Python — which cost
    51 min on the `cc=5, t=1` rehearsal's dominant component, against 8 m 57 s
    *before* this optimization existed. Rejected on that measurement.
    `TopTools_IndexedMapOfShape` is OCCT's own shape-identity map (same
    underlying TShape and location, ignoring orientation — what "the same
    edge" or "the same face" means here), so every membership test below is
    the map's own near-`O(1)` lookup, not a Python-level scan.
    """
    seam_lists: list[list] = []
    interior: list = []
    for faces in face_lists:
        shell = occ.faces_shell(faces)
        edge_faces = TopTools_IndexedDataMapOfShapeListOfShape()
        TopExp.MapShapesAndAncestors_s(
            shell, TopAbs_ShapeEnum.TopAbs_EDGE, TopAbs_ShapeEnum.TopAbs_FACE, edge_faces
        )
        # The single owning face of every free (used-once) edge in this tile —
        # exactly the faces round 2 can still affect. Collected once, up front,
        # from the map's own ancestor lists rather than re-derived per face.
        seam_faces = TopTools_IndexedMapOfShape()
        for i in range(1, edge_faces.Extent() + 1):
            owners = edge_faces.FindFromIndex(i)
            if owners.Extent() == 1:
                seam_faces.Add(owners.First())
        if seam_faces.Extent() == 0:
            interior.extend(faces)
            continue
        seam, rest = [], []
        for f in faces:
            (seam if seam_faces.Contains(f) else rest).append(f)
        seam_lists.append(seam)
        interior.extend(rest)
    return seam_lists, interior


def _sew_round_two(
    by_group: dict[int, list],
    plan: dict[int, list[list] | None],
    tile_results: dict[int, list[list]],
    tolerance: float,
    workers: int,
    tmpdir: str | None,
    pool: WorkerPool | None,
    expected_rings: dict[int, int] | None = None,
    stats: "SewStats | None" = None,
) -> tuple[dict[int, list], int, int]:
    """Sew each component's round-1 result (or its own pieces, if untiled)
    into that component's final boundary-layer shell.

    Dispatched per component across the shared worker pool rather than run
    serially on the master, the same lever round 1 already uses — components
    share no interface, so this is exact, embarrassingly parallel, and free
    generality. It is generality rather than *the* win on a part like the
    docs/specification.md §10 rehearsal, whose 21,955 pieces sit almost
    entirely in one dominant component: with one job there is nothing to
    parallelise across, and the gain is close to zero there specifically. The
    real win for that shape of part is :func:`_split_seam_interior`, applied
    below regardless of whether there is a pool to dispatch across — it is a
    genuine reduction in the work round 2 does, not merely a parallelism lever.

    **A hierarchical tree reduction was considered for round 2 itself — sewing
    tile results together pairwise across levels instead of in one call — and
    rejected without being built.** G6 (`tools/prototypes/RESULTS.md`) found
    round 2's cost tracks total face count almost flatly in shape count: going
    from 8 tiles to 8,000 pieces' worth of unified shells barely moves it,
    which is the signature of a cost dominated by a flat `a·F` term rather than
    G5a's `n^1.8` shape-count term (`n^1.8` at `n` in the tens is negligible
    next to `F` in the hundred-thousands). A tree cannot beat a floor every
    level pays in full: `L = ceil(log2(T))` levels of pairwise merges each pass
    all `F` faces through a sew, so a tree costs `L * a*F` against one call's
    `a*F` — and the final, root-level merge alone, over both halves of the
    tree, already costs the entirety of what one call costs today, before any
    of the other `L-1` levels are counted. Parallelism across workers hides the
    early levels but never the root. So a tree is strictly worse than what is
    built here, at every tile count that has been measured.

    **``expected_rings``, if given, is the safety net `_split_seam_interior`'s
    argument turned out to need in production.** It maps each component to how
    many interior-to-boundary interfaces it must present as free edges once
    correctly sewn — four apiece, since every such interface is the whole
    template cap quad (docs/algorithm.md §5.3(b)); every other cap a boundary
    piece carries (declined as unpaired, mismatched or unweldable, §7.1) keeps
    its own face, so it contributes no free edge, making this count exact
    rather than a lower bound. `_split_seam_interior`'s own argument — that
    sewing the seam-only subset in isolation reproduces exactly what a full
    round 2 would have done to those faces — held at every prototype scale
    tried (G8), all of them on lightly trimmed junctions. It does not hold on
    the real, heavily trimmed geometry a production part is made of: a seam
    face there can share an edge with a face the split carried through
    unchanged (a "straddling" edge, free within the seam-only subset even
    though the full tile uses it twice), and sewing that subset without its
    carried neighbour present lets `BRepBuilderAPI_Sewing` rebuild the edge
    onto a new `TopoDS_Edge` while the carried face keeps the original — one
    shared edge becomes two, each used once. Measured on the `cc=5, t=1`
    production rehearsal: 118,760 open edges at `assemble`, effectively all of
    them this (see the "Micron-scale debris edges" entry's superseding note,
    docs/specification.md §10). Verifying free-edge count against
    ``expected_rings`` and, on a mismatch, redoing that one component's round 2
    on the unsplit tile results (the behaviour before the split existed) makes
    that failure mode unreachable, at the cost of the split's saving only for
    the components where it was actually wrong.
    """
    groups = list(by_group.keys())
    face_lists = {
        group: (tile_results[group] if plan[group] is not None else [p.faces for p in by_group[group]])
        for group in groups
    }

    # Seam-only round 2 (G8): only where round 1 actually tiled the component —
    # an untiled component's "tile result" is just its own pieces' raw, never-
    # yet-sewn faces, which have no round-1 free-edge structure to exploit, so
    # it is sewn exactly as before this existed.
    t0 = time.perf_counter()
    seam_face_lists: dict[int, list[list]] = {}
    interior_faces: dict[int, list] = {}
    for group in groups:
        if plan[group] is not None:
            seam_face_lists[group], interior_faces[group] = _split_seam_interior(face_lists[group])
        else:
            seam_face_lists[group] = face_lists[group]
            interior_faces[group] = []
    if stats is not None:
        stats.t_split = time.perf_counter() - t0
    t0 = time.perf_counter()

    def _finish(sewn: dict[int, list]) -> dict[int, list]:
        return {g: sewn[g] + interior_faces[g] for g in groups}

    def _serial() -> tuple[dict[int, list], int]:
        return _finish({g: _sew_faces(seam_face_lists[g], tolerance) for g in groups}), 0

    if tmpdir is None or len(groups) < 2:
        out, max_rss = _serial()
    else:
        jobs = []
        for group in groups:
            in_path = os.path.join(tmpdir, f"sew_round2_{group}.brep")
            out_path = os.path.join(tmpdir, f"sew_round2_{group}_out.brep")
            _write_brep(occ.compound(occ.compound(fl) for fl in seam_face_lists[group]), in_path)
            jobs.append((in_path, out_path, tolerance))

        def _run(p: WorkerPool) -> tuple[dict[int, list], int]:
            raw, rss = p.run(_worker_sew_tile, jobs)
            sewn = {group: occ.faces(_read_brep(out_path)) for group, (out_path, _rss) in zip(groups, raw)}
            return _finish(sewn), rss

        if pool is not None and pool.active:
            out, max_rss = _run(pool)
        elif workers <= 1:
            out, max_rss = _serial()
        else:
            with WorkerPool(min(workers, len(jobs))) as owned:
                out, max_rss = _run(owned)

    if stats is not None:
        stats.t_round2 = time.perf_counter() - t0

    t0 = time.perf_counter()
    repaired = 0
    if expected_rings is not None:
        for group in groups:
            if plan[group] is None:
                continue  # never split, so there is nothing the split could have broken
            want = 4 * expected_rings.get(group, 0)
            got = len(free_edges(out[group]))
            if got != want:
                out[group] = _sew_faces(face_lists[group], tolerance)
                repaired += 1
                # The unsplit sew has just run, so its own free-edge count is
                # free to take — and it is the one number that says *which* of
                # the two things this check can catch actually happened. If the
                # unsplit result meets `want` and the split one did not, the
                # split was wrong (G9's straddling-edge mechanism). If neither
                # meets it, the split was not the problem and the ~11 min this
                # repair costs at rehearsal scale bought nothing
                # (docs/specification.md §10). Without these three numbers the
                # log could only say a count was wrong, never which.
                if stats is not None:
                    stats.repair_evidence.append(
                        (group, want, got, len(free_edges(out[group])))
                    )
    if stats is not None:
        stats.t_repair = time.perf_counter() - t0

    return out, max_rss, repaired


@dataclass
class SewStats:
    """What :func:`sew_boundary` did, for the run log (specification.md §1, §3)."""

    tiles: int = 0
    tiled_components: int = 0
    max_worker_rss: int = 0
    repaired_components: int = 0
    """How many tiled components' seam-only round 2 (:func:`_split_seam_interior`)
    left a free-edge count other than ``4 * its interior interfaces`` and were
    redone with a full, unsplit sew (:func:`_sew_round_two`'s ``expected_rings``
    check). Zero on every committed scenario; the check exists for real, heavily
    trimmed geometry where it has measured nonzero (docs/specification.md §10)."""
    repair_evidence: list = field(default_factory=list)
    """``(component, want, got_split, got_unsplit)`` for every repaired component.

    Recorded because :attr:`repaired_components` alone cannot distinguish the
    two things the check catches, and they call for opposite responses. If
    ``got_unsplit == want != got_split`` the seam-only split really did produce
    a different shell and the repair earned its cost. If ``got_unsplit ==
    got_split != want`` the split reproduced the unsplit sew exactly and the
    check fired on an expectation neither route can meet — the repair then
    costs a full sew (651 s of the `cc=5, t=1` rehearsal's `stitch`) and
    changes nothing. Both numbers are free: the unsplit sew has to run either
    way before this can be known."""
    retoleranced_faces: int = 0
    """Faces made valid again by correcting a vertex recorded off its edge's
    curve (:func:`latticegen2.occ.fix_vertex_tolerances`, docs/algorithm.md §8)."""
    t_round1: float = 0.0
    """Seconds in round 1 — the tiled, worker-parallel sew."""
    t_split: float = 0.0
    """Seconds in :func:`_split_seam_interior`, on the master."""
    t_round2: float = 0.0
    """Seconds in round 2's *attempt* — the seam-only sew, worker or serial."""
    t_repair: float = 0.0
    """Seconds redoing round 2 as a full unsplit sew, on the master.

    Zero unless ``repaired_components`` is nonzero. Split out from
    :attr:`t_round2` because the two answer different questions: `t_round2` is
    what the seam-only optimization costs, and on a part where the check fails
    that cost is **discarded** — the repair recomputes the component from
    scratch. Until these were measured separately, the run log could say
    `stitch_repaired_components: 1` without saying what that had cost, and
    profiling could only report the stage's 1.09 mean cores against a 5.42 peak
    without attributing the gap."""
    t_retolerance: float = 0.0
    """Seconds in :func:`latticegen2.occ.fix_vertex_tolerances`, on the master.

    One ``BRepCheck_Analyzer`` per boundary face, measured at 0.215 ms on real
    trimmed faces — ~1.1 min across the rehearsal's 301,505 (specification.md
    §10)."""
    still_invalid_faces: int = 0
    """Faces the analyzer rejects that this repair did **not** account for.

    Not a hard failure here: ``validate`` is already the gate for output
    validity (docs/algorithm.md §9), and failing twice for one cause would only
    obscure which check found it. Reported so a nonzero count is visible in the
    log *before* validate spends minutes rediscovering it."""


def sew_boundary(
    pieces,
    groups: list[int],
    *,
    workers: int = 1,
    tmpdir: str | None = None,
    tile_target: int = TILE_TARGET_PIECES,
    min_to_tile: int = MIN_PIECES_TO_TILE,
    pool: WorkerPool | None = None,
    want_rings: dict[tuple[NodeKey, int], int] | None = None,
) -> tuple[dict[int, list], SewStats]:
    """Sew each component's boundary pieces to each other, returning their faces.

    Per component rather than all at once: components share no interface, so this
    is exact, and it keeps every face attributable to the component it belongs to
    without a second search.

    A component whose piece count clears ``min_to_tile`` is split into spatial
    tiles by lattice-index block first (docs/algorithm.md §8): each tile is sewn
    on its own — every tiled component's round 1
    shares one worker pool for the whole call, not one pool per component — then
    each component's tile results are sewn together, itself dispatched per
    component across the same pool (:func:`_sew_round_two`). Round 1 shrinks
    roughly as G5a's measured `n^1.8` scaling predicts — a component of size `N`
    split into `T` tiles of about `N/T` each costs roughly `T * (N/T)^1.8` there,
    against `N^1.8` for one call. Round 2 is a real but smaller saving, not a
    free one: it sews shells whose *combined* face count still equals `N`, and G6
    (`tools/prototypes/RESULTS.md`) measured that cost staying close to a
    monolithic sew's rather than shrinking with `T` — 1.3–1.45× overall, not the
    order-of-magnitude round 1 alone would suggest. A component below the
    threshold, or whose whole footprint lands in one tile anyway, is sewn exactly
    as before this existed.

    ``pool``, if given, is the run's shared :class:`latticegen2.parallel.WorkerPool`,
    used for both rounds rather than each stage building and tearing down its own
    (docs/algorithm.md §8, §12).

    ``want_rings`` is ``(node, half-strut) -> component``, the same dict the
    caller passes to :func:`interface_rings` right after this — every interior
    interface the finished boundary shell must present as a free-edge ring. Given,
    it is turned into a per-component ring count and handed to
    :func:`_sew_round_two` as ``expected_rings``, so a seam-only split that
    silently produced the wrong result for a component is caught and repaired
    here rather than surfacing later as an unclosed shell out of ``assemble``
    with no indication of which component or why (docs/specification.md §10).
    Omitted, no verification runs — the pre-fix behaviour, kept for callers (and
    tests) that have no rings to check against.

    ``Cutting`` is switched off throughout — splitting free edges so they match is
    wasted work when they already match by construction, and G5a measures the
    whole optional-phase group at under 2 % either way.
    """
    by_group: dict[int, list] = {}
    for piece, group in zip(pieces, groups):
        by_group.setdefault(group, []).append(piece)

    plan = {
        group: _tile_pieces(group_pieces, tile_target, min_to_tile)
        for group, group_pieces in by_group.items()
    }
    stats = SewStats()
    for tiles in plan.values():
        if tiles is not None:
            stats.tiled_components += 1
            stats.tiles += len(tiles)

    ring_counts: dict[int, int] | None = None
    if want_rings is not None:
        ring_counts = {}
        for group in want_rings.values():
            ring_counts[group] = ring_counts.get(group, 0) + 1

    t0 = time.perf_counter()
    tile_results, max_rss1 = _sew_all_tiles(plan, SEW_TOLERANCE, workers, tmpdir, pool=pool)
    stats.t_round1 = time.perf_counter() - t0

    out, max_rss2, repaired = _sew_round_two(
        by_group, plan, tile_results, SEW_TOLERANCE, workers, tmpdir, pool,
        expected_rings=ring_counts, stats=stats,
    )
    stats.max_worker_rss = max(max_rss1, max_rss2)
    stats.repaired_components = repaired

    # Sewing can leave an edge whose vertex is recorded as sitting off the
    # edge's own curve, which fails BRepCheck on both faces sharing it. Repair
    # it here, on the sewn result, because that is where it is created — the
    # trimmed pieces going in are clean (docs/algorithm.md §8, G11). Doing it
    # before the interface rings are read means the interior adopts the
    # corrected vertices rather than a copy that would have to be fixed twice.
    t0 = time.perf_counter()
    for group_faces in out.values():
        fixed, residual = occ.fix_vertex_tolerances(group_faces)
        stats.retoleranced_faces += fixed
        stats.still_invalid_faces += residual
    stats.t_retolerance = time.perf_counter() - t0
    return out, stats


# --- step 2: read the interface rings off the sewn boundary ------------------


def free_edges(faces) -> list:
    """Edges used by exactly one of ``faces`` — the holes facing the interior.

    **Degenerate edges are not holes and are excluded**, for the same reason
    and by the same test :func:`shell_defects` uses: an edge with no extent is
    a parametric artefact whose owning face uses it once by construction, so
    counting it as a free edge counts something that is not there.

    That is not cosmetic. This count is what `_sew_round_two` compares against
    ``4 x interior interfaces`` to decide whether the seam-only split produced a
    correct shell, and on `TD_HX_rehearsal_test` at ``cc=5, t=1`` the trim
    against a grazing surface leaves exactly **10** degenerate edges (3.0e-9 to
    8.3e-8 mm, the same ones `shell_defects` records skipping). Counting them
    made even a *correct* unsplit sew read 73,994 against an expected 73,984 —
    a test that fires whatever the split does cannot report which of the two
    happened, and would force a full re-sew on a part whose split was fine
    (docs/specification.md §10).
    """
    shell = occ.faces_shell(faces)
    edge_faces = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(
        shell, TopAbs_ShapeEnum.TopAbs_EDGE, TopAbs_ShapeEnum.TopAbs_FACE, edge_faces
    )
    return [
        edge_faces.FindKey(i)
        for i in range(1, edge_faces.Extent() + 1)
        if not BRep_Tool.Degenerated_s(TopoDS.Edge_s(edge_faces.FindKey(i)))
        if edge_faces.FindFromIndex(i).Extent() == 1
    ]


def _corner(p: np.ndarray) -> tuple:
    """A hashable key for a corner position, quantised well inside WELD_TOL."""
    return (round(float(p[0]), 5), round(float(p[1]), 5), round(float(p[2]), 5))


def interface_rings(
    lp: LatticeParams,
    tmesh,
    boundary_faces: dict[int, list],
    wanted: dict[tuple[NodeKey, int], int],
) -> dict[tuple[NodeKey, int], tuple[list, list]]:
    """Locate each interior interface's hole on the sewn boundary shell.

    ``wanted`` maps ``(node, half-strut)`` to the component it belongs to.
    Returns, per interface, the ``(vertices, edges)`` in the template cap's loop
    order, ready for the instancing index to adopt.

    The lookup is by corner position rather than by search: an interior interface
    is the whole template quad, so its four corners are known exactly from the
    lattice expressions, and one dictionary keyed on those positions resolves
    every interface in the run.
    """
    index: dict[int, dict[tuple, list]] = {}
    for group, faces in boundary_faces.items():
        table: dict[tuple, list] = {}
        for edge in free_edges(faces):
            fwd = TopoDS.Edge_s(edge).Oriented(TopAbs_Orientation.TopAbs_FORWARD)
            vs = [TopoDS.Vertex_s(v)
                  for v in occ._explore(fwd, TopAbs_ShapeEnum.TopAbs_VERTEX)]
            if len(vs) != 2:
                continue
            key = frozenset((_corner(_pnt(vs[0])), _corner(_pnt(vs[1]))))
            table.setdefault(key, []).append((fwd, vs))
        index[group] = table

    out: dict[tuple[NodeKey, int], tuple[list, list]] = {}
    for (node, h), group in wanted.items():
        corners = template_cap_corners(lp, tmesh, node, h)
        table = index.get(group, {})
        verts, edges = [], []
        for i in range(len(corners)):
            a, b = corners[i], corners[(i + 1) % len(corners)]
            found = table.get(frozenset((_corner(a), _corner(b))))
            if not found or len(found) != 1:
                raise ProcessingError(
                    f"The boundary shell has no unique free edge for interface "
                    f"{node} cap {h} between {np.round(a, 4).tolist()} and "
                    f"{np.round(b, 4).tolist()}. The interior cannot be joined to "
                    f"it, so the output would not be watertight."
                )
            edge, vs = found[0]
            edges.append(edge)
            # In loop order the ring's i'th vertex is corner i.
            verts.append(vs[0] if np.linalg.norm(_pnt(vs[0]) - a) <= WELD_TOL else vs[1])
        out[(node, h)] = (verts, edges)
    return out


# --- step 4: assembly -------------------------------------------------------


def assemble(
    interior_shells: dict[int, TopoDS_Shell],
    boundary_faces: dict[int, list],
) -> dict[int, TopoDS_Shell]:
    """One shell per component, from the interior shells plus the sewn boundary."""
    builder = BRep_Builder()
    out: dict[int, TopoDS_Shell] = dict(interior_shells)
    for group, faces in boundary_faces.items():
        shell = out.get(group)
        if shell is None:
            shell = TopoDS_Shell()
            builder.MakeShell(shell)
            out[group] = shell
        for face in faces:
            builder.Add(shell, face)
    return out


def shell_defects(shell: TopoDS_Shell) -> tuple[int, int, list]:
    """``(open_edges, misoriented_edges, sample positions)`` for one shell.

    The same proof the instancing index applies to itself, and it needs **both**
    halves: a closed orientable surface uses every edge exactly twice *and*
    traverses it once each way. Counting uses alone is not enough — a face joined
    in back-to-front still gives every edge two users, so the shell "closes" and
    the volume silently comes out wrong. Measured while this module was being
    built: 0 open edges, 954 traversed the same way twice, volume 29,111 mm³
    against a true 51,393 mm³.
    """
    edges = TopTools_IndexedMapOfShape()
    TopExp.MapShapes_s(shell, TopAbs_ShapeEnum.TopAbs_EDGE, edges)
    n = edges.Extent()
    uses = [0] * (n + 1)
    net = [0] * (n + 1)
    for face in occ._explore(shell, TopAbs_ShapeEnum.TopAbs_FACE):
        # The explorer composes the face's own orientation into each edge, which
        # is exactly the traversal direction this test needs.
        for edge in occ._explore(face, TopAbs_ShapeEnum.TopAbs_EDGE):
            i = edges.FindIndex(edge)
            uses[i] += 1
            net[i] += 1 if edge.Orientation() == TopAbs_Orientation.TopAbs_FORWARD else -1

    open_edges = miso = 0
    samples: list[np.ndarray] = []
    by_use: dict[int, int] = {}
    for i in range(1, n + 1):
        if uses[i] == 2 and net[i] == 0:
            continue
        # A degenerate edge is a parametric artefact, not geometry: it has no
        # length, and the face that owns it uses it exactly once by
        # construction. Requiring two uses of it reports sound geometry as
        # broken. Measured on `TD_HX_rehearsal_test` at cc=5, t=1, where the trim
        # against a grazing surface leaves them: 10 of the 12 edges this test
        # rejected were degenerate, 3.0e-9 to 8.3e-8 mm long, and the shell was
        # closed everywhere they appeared. Skipping them is not a relaxation of
        # the proof — an edge with no extent cannot be a hole — and the check is
        # only as trustworthy as the quantity it compares (docs/algorithm.md
        # §11), which is exactly what issue #6 was about.
        if BRep_Tool.Degenerated_s(TopoDS.Edge_s(edges.FindKey(i))):
            continue
        if uses[i] != 2:
            open_edges += 1
            by_use[uses[i]] = by_use.get(uses[i], 0) + 1
        else:
            miso += 1
        if len(samples) < 10:
            pts = [_pnt(v) for v in occ._explore(edges.FindKey(i),
                                                 TopAbs_ShapeEnum.TopAbs_VERTEX)]
            if pts:
                samples.append(np.round(sum(pts) / len(pts), 3))
    # How many faces an edge has says *which* failure it is, and they need
    # opposite fixes: one face is a hole nothing filled, three or more is
    # material meeting where it should have been joined into one boundary.
    return open_edges, miso, samples, by_use


def close_shells(shells: dict[int, TopoDS_Shell]) -> tuple[list, dict]:
    """Verify every shell is watertight and coherent, then present them as solids."""
    solids = []
    open_edges = miso = 0
    samples: list[np.ndarray] = []
    by_use: dict[int, int] = {}
    for group in sorted(shells):
        shell = shells[group]
        bad_open, bad_orient, where, uses = shell_defects(shell)
        if bad_open or bad_orient:
            open_edges += bad_open
            miso += bad_orient
            for k, v in uses.items():
                by_use[k] = by_use.get(k, 0) + v
            samples.extend(where[: max(0, 10 - len(samples))])
            continue
        shell.Closed(True)
        solids.append(occ.make_solid(shell))

    if open_edges or miso:
        breakdown = ", ".join(f"{v} edge(s) on {k} face(s)" for k, v in sorted(by_use.items()))
        raise ProcessingError(
            f"The assembled lattice is not a closed orientable surface: "
            f"{open_edges} edge(s) are not used by exactly two faces ({breakdown}) "
            f"and {miso} are used twice the same way round. An edge on one face is "
            f"a hole nothing filled; on three or more, material met where it should "
            f"have been joined. Sample positions: {[p.tolist() for p in samples]}"
        )
    return solids, {"assembled_shells": len(shells)}
