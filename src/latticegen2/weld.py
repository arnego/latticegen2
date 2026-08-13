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
"""

from __future__ import annotations

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
absorb their disagreement, which is far below it."""


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
        for h, faces in piece.cap_faces.items():
            caps_at.setdefault((piece.node, h), []).extend(faces)

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


def sew_boundary(pieces, groups: list[int]) -> dict[int, list]:
    """Sew each component's boundary pieces to each other, returning their faces.

    Per component rather than all at once: components share no interface, so this
    is exact, it keeps every face attributable to the component it belongs to
    without a second search, and it is the natural unit to parallelise over if
    this ever needs it.

    ``Cutting`` is switched off — splitting free edges so they match is wasted
    work when they already match by construction, and G5a measures the whole
    optional-phase group at under 2 % either way.
    """
    by_group: dict[int, list] = {}
    for piece, group in zip(pieces, groups):
        by_group.setdefault(group, []).append(occ.faces_shell(piece.faces))
    out: dict[int, list] = {}
    for group, shells in by_group.items():
        sewn = occ.sew(shells, SEW_TOLERANCE, cutting=False)
        out[group] = occ.faces(sewn)
    return out


# --- step 2: read the interface rings off the sewn boundary ------------------


def free_edges(faces) -> list:
    """Edges used by exactly one of ``faces`` — the holes facing the interior."""
    shell = occ.faces_shell(faces)
    edge_faces = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(
        shell, TopAbs_ShapeEnum.TopAbs_EDGE, TopAbs_ShapeEnum.TopAbs_FACE, edge_faces
    )
    return [
        edge_faces.FindKey(i)
        for i in range(1, edge_faces.Extent() + 1)
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
