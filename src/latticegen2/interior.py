"""Interior lattice construction — one shell, built by index, with no booleans.

Every INTERIOR node contributes a translated copy of the junction template, and
adjacent copies meet exactly on their shared mid-strut cap quads. Both sides drop
that quad, so the two junctions' remaining faces bound a single continuous
volume.

The construction is *indexed*, not geometric. A global vertex is identified by
``(owning node, local template vertex)``, where a vertex sitting on an incoming
cap is re-attributed to the neighbour that owns the outgoing side of the same
cap. Two adjacent junctions therefore reference the *same* ``TopoDS_Vertex`` and
``TopoDS_Edge`` objects, and the shell is watertight by construction — there is
no tolerance, no search, and no possibility of a near-miss.

Why not just sew the instances together: measured, ``BRepBuilderAPI_Sewing``
takes 14.9 s for 1,000 junctions and grows clearly superlinearly (over 250 s of
CPU at 8,000), because it has to *discover* by geometric search the pairing that
is already known here exactly. Sewing is still the right tool where the pairing
genuinely isn't known — stitching trimmed boundary junctions on — but it must not
be on the path whose size scales with the volume of the part.

Why all six caps of an INTERIOR node can always be dropped: a node is INTERIOR
only when all six of its half-struts are, which puts every one of its cap planes
strictly inside the input solid; the node across such a cap therefore has a
half-strut that is not OUTSIDE and is never classified OUTSIDE itself
(:mod:`latticegen2.classify`). The "neighbour absent" branch below is unreachable
for a real INTERIOR node — it is kept so that a synthetic grid still closes and
so a classification bug degrades into a closed, oversized solid rather than a
hole in the output.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from OCP.BRep import BRep_Builder, BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeVertex
from OCP.BRepTools import BRepTools, BRepTools_WireExplorer
from OCP.GeomAbs import GeomAbs_SurfaceType
from OCP.TopAbs import TopAbs_Orientation
from OCP.TopoDS import TopoDS_Shell, TopoDS_Wire
from OCP.gp import gp_Dir, gp_Pln

from . import occ
from .errors import ProcessingError
from .junction import JunctionTemplate
from .lattice import HALF_STRUTS, LatticeParams, half_strut_offset, neighbor_step, nodes, profile_vertices

_WELD = 9
"""Decimal places used to identify coincident template vertices. The template is
built from exact expressions at millimetre scale, so coincident vertices agree
far beyond this; it only guards against a last-ulp difference."""


@dataclass
class TemplateMesh:
    """The junction template as plain indexed polyhedral data."""

    verts: np.ndarray
    """``(V, 3)`` local vertex coordinates."""
    loops: list[list[int]] = field(default_factory=list)
    """Per face, its outer wire as vertex indices wound so the right-hand
    normal points *out* of the solid."""
    face_cap: list[int] = field(default_factory=list)
    """Per face, the half-strut id of the cap it is, or ``-1``."""
    face_normal: list[np.ndarray] = field(default_factory=list)
    """Per face, its outward unit normal in template-local coordinates.

    Carried explicitly because the face has to be *rebuilt* on a plane of this
    exact orientation. Letting OCCT infer the plane from the wire picks an
    arbitrary normal direction, which silently produces a shell whose faces
    point every which way — closed, but enclosing zero volume."""
    vertex_cap: list[int] = field(default_factory=list)
    """Per vertex, the half-strut id of the cap it lies on, or ``-1``."""
    cap_partner: dict[int, int] = field(default_factory=dict)
    """Vertex on an incoming cap -> the matching vertex on the outgoing cap."""


def _wire_loop(face) -> list:
    """Ordered vertex points of a face's outer wire."""
    points = []
    exp = BRepTools_WireExplorer(BRepTools.OuterWire_s(face))
    while exp.More():
        p = BRep_Tool.Pnt_s(exp.CurrentVertex())
        points.append((p.X(), p.Y(), p.Z()))
        exp.Next()
    return points


def _outward_normal(face) -> np.ndarray:
    surf = BRepAdaptor_Surface(face)
    if surf.GetType() != GeomAbs_SurfaceType.GeomAbs_Plane:
        raise ProcessingError(
            "The junction template contains a non-planar face. Struts are prisms "
            "over planar profiles, so every face of the template must be a plane."
        )
    d = surf.Plane().Axis().Direction()
    n = np.array([d.X(), d.Y(), d.Z()], dtype=float)
    if face.Orientation() == TopAbs_Orientation.TopAbs_REVERSED:
        n = -n
    return n


def _newell(points: np.ndarray) -> np.ndarray:
    """Polygon normal by Newell's method — robust for non-convex loops."""
    rolled = np.roll(points, -1, axis=0)
    n = np.array(
        [
            np.sum((points[:, 1] - rolled[:, 1]) * (points[:, 2] + rolled[:, 2])),
            np.sum((points[:, 2] - rolled[:, 2]) * (points[:, 0] + rolled[:, 0])),
            np.sum((points[:, 0] - rolled[:, 0]) * (points[:, 1] + rolled[:, 1])),
        ]
    )
    return n


def extract_template_mesh(lp: LatticeParams, tpl: JunctionTemplate) -> TemplateMesh:
    """Turn the fused junction into indexed data the instancing loop can use."""
    index: dict[tuple, int] = {}
    coords: list[tuple[float, float, float]] = []

    def vid(p) -> int:
        key = (round(p[0], _WELD), round(p[1], _WELD), round(p[2], _WELD))
        if key not in index:
            index[key] = len(coords)
            coords.append((float(p[0]), float(p[1]), float(p[2])))
        return index[key]

    mesh = TemplateMesh(verts=np.empty((0, 3)))
    for face in occ.faces(tpl.solid):
        points = _wire_loop(face)
        if len(points) < 3:
            raise ProcessingError("Junction template has a degenerate face wire.")
        loop = [vid(p) for p in points]
        arr = np.array(points, dtype=float)
        normal = _outward_normal(face)
        if float(np.dot(_newell(arr), normal)) < 0:
            loop.reverse()
        mesh.loops.append(loop)
        mesh.face_normal.append(normal)
        mesh.face_cap.append(_which_cap(lp, face, tpl))
    mesh.verts = np.array(coords, dtype=float)

    if sum(1 for c in mesh.face_cap if c >= 0) != 6:
        raise ProcessingError(
            f"Expected exactly 6 cap faces in the junction template, found "
            f"{sum(1 for c in mesh.face_cap if c >= 0)}."
        )

    mesh.vertex_cap = _classify_vertices(lp, mesh.verts)
    mesh.cap_partner = _pair_caps(lp, mesh)
    return mesh


def _which_cap(lp: LatticeParams, face, tpl: JunctionTemplate) -> int:
    for h, cap in enumerate(tpl.cap_faces):
        if face.IsSame(cap):
            return h
    return -1


def _classify_vertices(lp: LatticeParams, verts: np.ndarray) -> list[int]:
    """Mark each template vertex with the cap it lies on, if any."""
    out = [-1] * len(verts)
    for h, (k, _) in enumerate(HALF_STRUTS):
        corners = profile_vertices(lp, half_strut_offset(lp, h), k)
        for corner in corners:
            d = np.linalg.norm(verts - corner, axis=1)
            i = int(np.argmin(d))
            if d[i] > 1e-7:
                raise ProcessingError(
                    f"Junction template is missing a corner of cap {h}; the cap "
                    f"quad did not survive the fuse intact."
                )
            out[i] = h
    return out


def _pair_caps(lp: LatticeParams, mesh: TemplateMesh) -> dict[int, int]:
    """Map each incoming-cap vertex to the outgoing-cap vertex it coincides with.

    Node ``n``'s cap ``h >= 3`` occupies the same square as node
    ``n + neighbor_step(h)``'s cap ``h - 3``; in template-local terms the two
    differ by exactly ``a * e_k``. The pairing is computed once here so the
    instancing loop performs integer lookups only, never coordinate matching.
    """
    partner: dict[int, int] = {}
    for h in range(3, 6):
        k, _ = HALF_STRUTS[h]
        shift = lp.a * lp.e[k]
        for i, cap in enumerate(mesh.vertex_cap):
            if cap != h:
                continue
            target = mesh.verts[i] + shift
            d = np.linalg.norm(mesh.verts - target, axis=1)
            j = int(np.argmin(d))
            if d[j] > 1e-7 or mesh.vertex_cap[j] != h - 3:
                raise ProcessingError(
                    f"Could not pair cap {h} with cap {h - 3} in the junction "
                    f"template; the two mid-strut faces are not congruent."
                )
            partner[i] = j
    return partner


class _ShellBuilder:
    """Accumulates shared vertices, edges and faces into one OCCT shell.

    Edge usage is tallied as faces are added, which gives an exact closure test
    for free: a closed orientable surface uses every edge exactly twice, once in
    each direction. Edges used once are the shell's genuine open boundary — the
    square holes where boundary junctions attach — and anything else (an edge
    used three times, or twice the same way) means the index is wrong and the
    shell must not be presented as a solid.
    """

    def __init__(self):
        self.builder = BRep_Builder()
        self.shell = TopoDS_Shell()
        self.builder.MakeShell(self.shell)
        self._vertices: dict[tuple, object] = {}
        self._edges: dict[tuple[int, int], object] = {}
        self._vid: dict[tuple, int] = {}
        self._use: dict[tuple[int, int], int] = {}
        self._net: dict[tuple[int, int], int] = {}
        self._pos: dict[int, np.ndarray] = {}
        self.n_faces = 0

    def vertex(self, key: tuple, position) -> int:
        """Intern a global vertex, creating its ``TopoDS_Vertex`` on first use."""
        got = self._vid.get(key)
        if got is None:
            got = len(self._vid)
            self._vid[key] = got
            pos = np.asarray(position, dtype=float)
            self._pos[got] = pos
            self._vertices[got] = BRepBuilderAPI_MakeVertex(occ._pnt(pos)).Vertex()
        return got

    def _edge(self, a: int, b: int):
        key = (a, b) if a < b else (b, a)
        edge = self._edges.get(key)
        if edge is None:
            edge = BRepBuilderAPI_MakeEdge(self._vertices[key[0]], self._vertices[key[1]]).Edge()
            self._edges[key] = edge
        return edge, (a < b)

    def face(self, loop: list[int], normal: np.ndarray) -> None:
        wire = TopoDS_Wire()
        self.builder.MakeWire(wire)
        for i in range(len(loop)):
            a, b = loop[i], loop[(i + 1) % len(loop)]
            if a == b:
                raise ProcessingError("Degenerate edge in an instanced junction face.")
            edge, forward = self._edge(a, b)
            key = (a, b) if forward else (b, a)
            self._use[key] = self._use.get(key, 0) + 1
            self._net[key] = self._net.get(key, 0) + (1 if forward else -1)
            self.builder.Add(
                wire,
                edge if forward else edge.Oriented(TopAbs_Orientation.TopAbs_REVERSED),
            )
        # The plane is stated, never inferred: its normal is the template face's
        # outward normal, so the resulting face is FORWARD with that normal and
        # its wire traverses by the right-hand rule about it — which is exactly
        # the winding `loop` already has. Every shared edge is therefore
        # traversed oppositely by its two faces, as a valid shell requires.
        plane = gp_Pln(occ._pnt(self._pos[loop[0]]),
                       gp_Dir(float(normal[0]), float(normal[1]), float(normal[2])))
        maker = BRepBuilderAPI_MakeFace(plane, wire)
        if not maker.IsDone():
            raise ProcessingError("Could not build a planar face from an instanced wire.")
        self.builder.Add(self.shell, maker.Face())
        self.n_faces += 1

    def closure(self) -> tuple[int, int]:
        """``(free_edges, malformed_edges)`` from the edge-use tally.

        A free edge is used exactly once — the expected open boundary. A
        malformed edge is used more than twice, or twice in the same direction,
        either of which means the surface is not a valid orientable manifold.
        """
        free = 0
        bad = 0
        for key, count in self._use.items():
            if count == 1:
                free += 1
            elif count == 2 and self._net[key] == 0:
                continue
            else:
                bad += 1
        return free, bad

    def result(self) -> tuple[TopoDS_Shell, int]:
        """``(shell, free_edge_count)``, with OCCT's ``Closed`` flag set correctly.

        ``BRep_Builder`` leaves the flag false on a hand-built shell regardless
        of the geometry, and downstream code (``BRepBuilderAPI_MakeSolid``,
        validity checking) reads the flag rather than recomputing it.

        The free-edge count is returned rather than left for the caller to ask
        for again: ``closure()`` walks every edge, and at scale that map holds
        millions of them.
        """
        free, bad = self.closure()
        if bad:
            raise ProcessingError(
                f"{bad} edge(s) of the interior shell are used more than twice or "
                f"twice in the same direction. The instanced lattice is not a valid "
                f"orientable surface, which would make the output non-manifold."
            )
        self.shell.Closed(free == 0)
        return self.shell, free

    @property
    def n_vertices(self) -> int:
        return len(self._vid)

    @property
    def n_edges(self) -> int:
        return len(self._edges)


def build_interior_shell(
    lp: LatticeParams,
    tpl: JunctionTemplate,
    tmesh: TemplateMesh,
    interior_nodes: np.ndarray,
    kept: set[tuple[int, int, int]],
) -> tuple[TopoDS_Shell, dict]:
    """Build the whole interior lattice as one shared-topology shell.

    Returns the shell (open exactly where it meets a boundary junction) and a
    small stats dict for the run log.
    """
    if len(interior_nodes) == 0:
        return None, {
            "interior_faces": 0,
            "interior_vertices": 0,
            "interior_edges": 0,
            "interior_open_edges": 0,
        }

    positions = nodes(lp, interior_nodes)
    steps = [tuple(int(x) for x in neighbor_step(h)) for h in range(6)]
    vertex_cap = tmesh.vertex_cap
    partner = tmesh.cap_partner
    verts = tmesh.verts
    builder = _ShellBuilder()

    # Owning-node coordinates are cached so a shared vertex is positioned from
    # one expression only, whichever junction reaches it first.
    node_pos: dict[tuple[int, int, int], np.ndarray] = {}
    for idx in range(len(interior_nodes)):
        node_pos[tuple(int(x) for x in interior_nodes[idx])] = positions[idx]

    def owner(node: tuple[int, int, int], local: int):
        cap = vertex_cap[local]
        if cap >= 3:
            step = steps[cap]
            neighbour = (node[0] + step[0], node[1] + step[1], node[2] + step[2])
            return neighbour, partner[local]
        return node, local

    def position(node: tuple[int, int, int], local: int) -> np.ndarray:
        base = node_pos.get(node)
        if base is None:
            base = lp.B @ np.array(node, dtype=float)
            node_pos[node] = base
        return base + verts[local]

    for idx in range(len(interior_nodes)):
        node = tuple(int(x) for x in interior_nodes[idx])
        gid: dict[int, int] = {}
        for fi, loop in enumerate(tmesh.loops):
            cap = tmesh.face_cap[fi]
            if cap >= 0:
                step = steps[cap]
                neighbour = (node[0] + step[0], node[1] + step[1], node[2] + step[2])
                if neighbour in kept:
                    continue  # interface: dropped from both sides
            out_loop = []
            for local in loop:
                got = gid.get(local)
                if got is None:
                    key = owner(node, local)
                    got = builder.vertex(key, position(*key))
                    gid[local] = got
                out_loop.append(got)
            builder.face(out_loop, tmesh.face_normal[fi])

    shell, free = builder.result()
    stats = {
        "interior_faces": builder.n_faces,
        "interior_vertices": builder.n_vertices,
        "interior_edges": builder.n_edges,
        "interior_open_edges": free,
    }
    return shell, stats
