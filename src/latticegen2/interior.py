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

Which caps get dropped is **not** decided here. It is handed in as the interface
set resolved from both sides at once
(:func:`latticegen2.boundary.resolve_interfaces`), because a cap may only be
dropped when the junction across it drops the matching one. An INTERIOR node's
six caps normally all qualify — a node is INTERIOR only when all six of its
half-struts are, which puts every one of its cap planes strictly inside the input
solid, so the trimmed junction across such a cap keeps that cap whole
(docs/algorithm.md §5.3(b)). Where a boolean has nonetheless left nothing on the
other side, the cap is simply kept, and the shell closes over it rather than
opening a hole nothing will fill.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from OCP.BRep import BRep_Builder, BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeVertex
from OCP.BRepTools import BRepTools, BRepTools_WireExplorer
from OCP.GeomAbs import GeomAbs_SurfaceType
from OCP.TopAbs import TopAbs_Orientation, TopAbs_ShapeEnum
from OCP.TopoDS import TopoDS, TopoDS_Shell, TopoDS_Wire
from OCP.gp import gp_Dir, gp_Pln

from . import occ
from .errors import ProcessingError
from .junction import JunctionTemplate
from .lattice import (
    HALF_STRUTS,
    OPPOSITE_HALF,
    LatticeParams,
    half_strut_offset,
    neighbor_step,
    nodes,
    profile_vertices,
)

_WELD = 9
"""Decimal places used to identify coincident template vertices. The template is
built from exact expressions at millimetre scale, so coincident vertices agree
far beyond this; it only guards against a last-ulp difference."""

_ADOPT_TOL = 1e-6
"""Millimetres. Whether an adopted edge starts at one corner or the other — a
question about which of two corners a whole strut-width apart it is nearer, so
the bar only has to be far below that."""


def _first_vertex_position(edge) -> np.ndarray:
    """Where the edge's own parametrisation starts, ignoring how it is used."""
    fwd = TopoDS.Edge_s(edge).Oriented(TopAbs_Orientation.TopAbs_FORWARD)
    first = next(occ._explore(fwd, TopAbs_ShapeEnum.TopAbs_VERTEX))
    p = BRep_Tool.Pnt_s(TopoDS.Vertex_s(first))
    return np.array([p.X(), p.Y(), p.Z()])


@dataclass
class MergedLateral:
    """One full-strut lateral face, spliced from the two half-strut faces.

    ``loop`` holds ``(side, local vertex)`` pairs, where side 0 is the node that
    owns the outgoing half-strut and side 1 is its neighbour across the shared
    cap. Resolving a side to a node is the caller's job, so one of these serves
    every strut in the lattice — the same instancing argument as the template
    itself (docs/algorithm.md §3.2).
    """

    loop: list[tuple[int, int]] = field(default_factory=list)
    normal: np.ndarray = None


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
    face_half: list[int] = field(default_factory=list)
    """Per face, the half-strut whose cap edge this lateral face carries, or
    ``-1`` for a cap face (or a lateral face that could not be attributed)."""
    merged: dict[int, MergedLateral] = field(default_factory=dict)
    """Outgoing-side lateral face index -> its merged full-strut face."""
    merged_halves: set[int] = field(default_factory=set)
    """Outgoing half-struts (``0..2``) whose four merged faces all validated.

    Membership is per strut family rather than global, so a splice that fails
    validation costs that family its merge and nothing else — the failure mode
    docs/algorithm.md §11 requires of an optimization."""
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
    mesh.face_half = _classify_lateral_faces(mesh)
    mesh.merged, mesh.merged_halves = _merge_lateral_pairs(lp, mesh)
    return mesh


# --- Full-strut lateral faces (docs/algorithm.md §6) -------------------------
#
# Instancing emits four lateral faces per half-strut, so every strut whose two
# ends are both instanced carries eight where four suffice, and same-domain
# unification spends the `simplify` stage merging them back
# (docs/algorithm.md §9). Emitting the merged face directly is strictly cheaper:
# the pair is known by construction — one per surviving mid-strut interface —
# where the kernel has to rediscover it by search over the whole solid.
#
# The splice is computed once here, on the template, as a fixed pattern of
# `(side, local vertex)` indices, exactly the shape of precomputation
# `_pair_caps` already does. At run time it is integer lookups only.

_COLLINEAR_SIN_TOL = 1e-9
"""Bar on ``|u x v| / (|u| |v|)`` — the sine of the turn at a vertex — for
calling three consecutive loop points collinear.

The points being tested are exact expressions at millimetre scale and the two
edges meeting at a dropped cap corner are both *along the strut axis*, so the
true sine is at ulp level; a genuine corner of this geometry turns by tens of
degrees. Nothing real sits between the two, which is why a single loose-looking
bar is safe here."""

_PLANAR_TOL = 1e-9
"""Millimetres a merged loop's points may deviate from the face plane."""

_AREA_REL_TOL = 1e-12
"""Relative bar on ``area(merged) == area(half) + area(half)``.

Both sides are computed from the same exact vertex coordinates by the same
Newell sum, so only floating-point summation order separates them."""


def _polygon_area(points: np.ndarray) -> float:
    return 0.5 * float(np.linalg.norm(_newell(points)))


def _classify_lateral_faces(mesh: TemplateMesh) -> list[int]:
    """Per face, the half-strut whose cap quad supplies its far edge.

    A lateral face of half-strut ``h`` carries exactly two vertices of cap
    ``h`` — one edge of that cap square — and a cap face carries four. Anything
    else is left unattributed (``-1``) and simply never merges, which is why
    this needs no tolerance: it is a count, not a measurement.
    """
    out: list[int] = []
    for fi, loop in enumerate(mesh.loops):
        if mesh.face_cap[fi] >= 0:
            out.append(-1)
            continue
        counts: dict[int, int] = {}
        for v in loop:
            c = mesh.vertex_cap[v]
            if c >= 0:
                counts[c] = counts.get(c, 0) + 1
        owners = [h for h, n in counts.items() if n == 2]
        out.append(owners[0] if len(owners) == 1 else -1)
    return out


def _cap_edge_index(mesh: TemplateMesh, loop: list[int], h: int) -> int | None:
    """Index ``i`` where ``loop[i] -> loop[i+1]`` is the cap-``h`` edge."""
    n = len(loop)
    hits = [
        i
        for i in range(n)
        if mesh.vertex_cap[loop[i]] == h and mesh.vertex_cap[loop[(i + 1) % n]] == h
    ]
    return hits[0] if len(hits) == 1 else None


def _drop_collinear(loop: list[tuple[int, int]], points) -> list[tuple[int, int]]:
    """Remove loop entries whose two neighbours are collinear with them.

    The two cap corners are the point of this: in the merged face the edges
    meeting there both run along the strut axis, so the corner stops being one.
    Dropping them is what makes the merged wire *minimal* and delivers the edge
    reduction the kernel's own edge pass would otherwise have to be paid for
    (docs/algorithm.md §9). Written as a general collinearity test rather than
    "delete the two cap corners" so the loop is minimal whatever the template
    turns out to contain.
    """
    keep = list(loop)
    pts = [np.asarray(p, dtype=float) for p in points]
    i = 0
    while len(keep) > 3 and i < len(keep):
        p, c, q = pts[i - 1], pts[i], pts[(i + 1) % len(keep)]
        u, v = c - p, q - c
        nu, nv = float(np.linalg.norm(u)), float(np.linalg.norm(v))
        if nu > 0 and nv > 0 and float(np.linalg.norm(np.cross(u, v))) <= _COLLINEAR_SIN_TOL * nu * nv:
            del keep[i]
            del pts[i]
            i = 0  # a removal can make its neighbours collinear in turn
            continue
        i += 1
    return keep


def _merge_lateral_pairs(
    lp: LatticeParams, mesh: TemplateMesh
) -> tuple[dict[int, MergedLateral], set[int]]:
    """Splice each outgoing half-strut's lateral faces with the neighbour's.

    Returns the merged faces keyed by the *outgoing* side's template face index,
    plus the set of half-struts for which all four splices validated. A family
    that fails validation is simply absent, and the build falls back to the two
    half-faces for it.
    """
    merged: dict[int, MergedLateral] = {}
    halves: set[int] = set()
    for h in range(3):
        k, _ = HALF_STRUTS[h]
        shift = lp.a * lp.e[k]
        opp = OPPOSITE_HALF[h]
        mine = [fi for fi, hh in enumerate(mesh.face_half) if hh == h]
        theirs = [fi for fi, hh in enumerate(mesh.face_half) if hh == opp]
        if len(mine) != len(theirs) or not mine:
            continue
        spliced = {}
        for fi in mine:
            got = _splice_lateral(mesh, fi, theirs, h, opp, shift)
            if got is None:
                spliced = None
                break
            spliced[fi] = got
        if spliced:
            merged.update(spliced)
            halves.add(h)
    return merged, halves


def _splice_lateral(mesh, fi, candidates, h, opp, shift) -> MergedLateral | None:
    """Join lateral face ``fi`` to the neighbour face sharing its cap edge.

    The two faces share that edge already — it is why the instanced shell is
    watertight — and a closed orientable surface traverses a shared edge
    oppositely from its two sides, so the neighbour presents ``b -> a`` where
    this face presents ``a -> b``. The union's boundary is therefore each
    face's boundary minus that edge, walked one after the other.
    """
    f1 = mesh.loops[fi]
    i = _cap_edge_index(mesh, f1, h)
    if i is None:
        return None
    n = len(f1)
    a, b = f1[i], f1[(i + 1) % n]

    for fj in candidates:
        f2 = mesh.loops[fj]
        j = _cap_edge_index(mesh, f2, opp)
        if j is None:
            continue
        m = len(f2)
        x, y = f2[j], f2[(j + 1) % m]
        # Mapped through the cap correspondence, the neighbour's cap edge must
        # be this face's, reversed. Anything else is a different lateral side.
        if mesh.cap_partner.get(x) != b or mesh.cap_partner.get(y) != a:
            continue

        normal = mesh.face_normal[fi]
        if float(np.dot(normal, mesh.face_normal[fj])) < 1 - 1e-12:
            return None  # not the same lateral plane after all

        # f1 walked from b round to a (every f1 edge but a->b), then f2 walked
        # from just past a round to just before b (every f2 edge but b->a).
        loop = [(0, f1[(i + 1 + s) % n]) for s in range(n)]
        loop += [(1, f2[(j + 2 + s) % m]) for s in range(m - 2)]

        def pos(entry):
            side, v = entry
            return mesh.verts[v] + (shift if side else 0.0)

        loop = _drop_collinear(loop, [pos(e) for e in loop])
        pts = np.array([pos(e) for e in loop])
        if len(pts) < 3:
            return None
        if float(np.abs((pts - pts[0]) @ normal).max()) > _PLANAR_TOL:
            return None
        nv = _newell(pts)
        area = 0.5 * float(np.linalg.norm(nv))
        if area <= 0 or float(np.dot(nv / np.linalg.norm(nv), normal)) < 1 - 1e-9:
            return None
        want = _polygon_area(mesh.verts[f1]) + _polygon_area(mesh.verts[f2])
        if abs(area - want) > _AREA_REL_TOL * want:
            return None
        return MergedLateral(loop=loop, normal=normal)
    return None


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
        self.shells: dict[int, TopoDS_Shell] = {}
        self._vertices: dict[tuple, object] = {}
        self._edges: dict[tuple[int, int], object] = {}
        self._vid: dict[tuple, int] = {}
        self._use: dict[tuple[int, int], int] = {}
        self._net: dict[tuple[int, int], int] = {}
        self._pos: dict[int, np.ndarray] = {}
        self.n_faces = 0

    def shell_for(self, group: int) -> TopoDS_Shell:
        """The shell collecting ``group``'s faces, created on first use.

        Faces are separated by junction-graph component as they are built, which
        is exact and free: two components share no interface, so no edge and no
        vertex can be common to both. The alternative — assemble everything and
        rediscover the split geometrically — is the kind of search this module
        exists to avoid.
        """
        got = self.shells.get(group)
        if got is None:
            got = TopoDS_Shell()
            self.builder.MakeShell(got)
            self.shells[group] = got
        return got

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

    def adopt(self, keys: list[tuple], positions, verts, edges) -> None:
        """Take an existing ring's vertices and edges as this index's own.

        Used where the interior meets a boundary junction: rather than build our
        own topology there and try to reconcile the two afterwards, the interior
        is built *on* the boundary's, so the two are the same objects from the
        start and nothing needs reconciling.

        The direction is not arbitrary. `BRepTools_ReShape` will happily swap an
        edge inside a face, but replacing that edge's **vertices** leaves the
        neighbouring edges still referring to the old ones, and the wire comes
        apart — `BRepCheck_NotConnected`, measured, with the volume quietly wrong
        and every edge still used exactly twice. So the side that adapts has to
        be the side whose faces we build ourselves, which is this one.

        ``keys``, ``positions``, ``verts`` and ``edges`` are in loop order;
        ``edges[i]`` joins ``keys[i]`` to ``keys[i + 1]``.
        """
        ids = []
        for key, pos, vert in zip(keys, positions, verts):
            got = self._vid.get(key)
            if got is None:
                got = len(self._vid)
                self._vid[key] = got
            self._pos[got] = np.asarray(pos, dtype=float)
            # A vertex explored out of an edge carries FORWARD/REVERSED to say
            # whether it is that edge's start or end. Passed on to
            # `BRepBuilderAPI_MakeEdge` that flag swaps the new edge's endpoints,
            # so edges we build on an adopted corner come out running the wrong
            # way — and it is the *interior's own* edges that then disagree.
            self._vertices[got] = TopoDS.Vertex_s(
                vert.Oriented(TopAbs_Orientation.TopAbs_FORWARD))
            ids.append(got)
        for i, edge in enumerate(edges):
            a, b = ids[i], ids[(i + 1) % len(ids)]
            # `_edge` keys on the sorted id pair and `face` decides which way to
            # traverse from that order alone, so every stored edge must run from
            # the lower id to the higher one. An edge we made ourselves does by
            # construction; an adopted one runs whichever way its own maker chose,
            # and storing it unexamined makes the interior traverse it the same
            # way round as the boundary does — the shell still closes, and the
            # surface is no longer orientable.
            natural = _first_vertex_position(edge)
            forward = float(np.linalg.norm(natural - self._pos[a])) <= _ADOPT_TOL
            if forward != (a < b):
                edge = edge.Reversed()
            self._edges[(a, b) if a < b else (b, a)] = edge

    def _edge(self, a: int, b: int):
        key = (a, b) if a < b else (b, a)
        edge = self._edges.get(key)
        if edge is None:
            edge = BRepBuilderAPI_MakeEdge(self._vertices[key[0]], self._vertices[key[1]]).Edge()
            self._edges[key] = edge
        return edge, (a < b)

    def ring(self, loop: list[int]) -> tuple[list, list]:
        """``(vertices, edges)`` of a cap loop that was dropped, by index only.

        The four edges already exist: the cap face was skipped, but the four
        lateral faces around it each created one of them. So the topology a
        boundary junction has to be welded onto is four dictionary lookups, with
        no geometric search and no tolerance anywhere in it (docs/algorithm.md
        §8).
        """
        verts = [self._vertices[v] for v in loop]
        edges = []
        for i in range(len(loop)):
            edge, _ = self._edge(loop[i], loop[(i + 1) % len(loop)])
            edges.append(edge)
        return verts, edges

    def face(self, loop: list[int], normal: np.ndarray, group: int = 0) -> None:
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
            # `Reversed()` flips relative to the edge's current flag; `Oriented`
            # would set it absolutely and throw away an adopted edge's own
            # direction, which is not ours to choose.
            self.builder.Add(wire, edge if forward else edge.Reversed())
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
        self.builder.Add(self.shell_for(group), maker.Face())
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

    def result(self) -> tuple[dict[int, TopoDS_Shell], int]:
        """``(shells by group, free_edge_count)``.

        The free-edge count is returned rather than left for the caller to ask
        for again: ``closure()`` walks every edge, and at scale that map holds
        millions of them. The ``Closed`` flag is deliberately *not* set here —
        these shells are open at every cap a boundary junction will be welded
        into, and the flag is set once the assembly is complete
        (:mod:`latticegen2.weld`).
        """
        free, bad = self.closure()
        if bad:
            raise ProcessingError(
                f"{bad} edge(s) of the interior shell are used more than twice or "
                f"twice in the same direction. The instanced lattice is not a valid "
                f"orientable surface, which would make the output non-manifold."
            )
        return self.shells, free

    @property
    def n_vertices(self) -> int:
        return len(self._vid)

    @property
    def n_edges(self) -> int:
        return len(self._edges)


@dataclass
class InteriorShells:
    """The instanced interior, split by component, plus its weldable rings."""

    shells: dict[int, TopoDS_Shell] = field(default_factory=dict)
    """Component id -> that component's shell, open at every boundary interface."""
    stats: dict = field(default_factory=dict)


def build_interior_shell(
    lp: LatticeParams,
    tpl: JunctionTemplate,
    tmesh: TemplateMesh,
    interior_nodes: np.ndarray,
    interfaces: set[tuple[tuple[int, int, int], int]],
    groups: dict[tuple[int, int, int], int] | None = None,
    adopted: dict[tuple[tuple[int, int, int], int], tuple[list, list]] | None = None,
) -> InteriorShells:
    """Build the interior lattice as shared-topology shells, one per component.

    ``interfaces`` holds ``(node, half-strut)`` for every cap that is stitched
    across; a cap face is dropped exactly when its key is in it, and kept
    otherwise. ``groups`` maps a node to its junction-graph component, so the
    faces are separated as they are built rather than rediscovered afterwards.

    ``adopted`` supplies, per interface cap, the ``(vertices, edges)`` the
    boundary side already has there, in the template cap's loop order. The
    interior is then built directly on them, which is what joins the two without
    any search, tolerance or repair — see :meth:`_ShellBuilder.adopt`.
    """
    if len(interior_nodes) == 0:
        return InteriorShells(stats={
            "interior_faces": 0,
            "interior_vertices": 0,
            "interior_edges": 0,
            "interior_open_edges": 0,
        })

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

    # Taken from `interior_nodes` and never from `node_pos`, which `position`
    # below grows with any neighbour it is asked about — including boundary
    # nodes reached through the cap correspondence, which are emphatically not
    # instanced here. Reading membership out of that cache instead let an
    # interior junction merge its lateral faces with a boundary neighbour that
    # was never built, leaving the interface ring with nothing behind it: 366
    # unmatched edges on the 80 mm ball, caught by weld's every-edge-twice proof.
    interior_set = frozenset(node_pos)

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

    def globals_of(node, loop, gid) -> list[int]:
        out_loop = []
        for local in loop:
            got = gid.get(local)
            if got is None:
                key = owner(node, local)
                got = builder.vertex(key, position(*key))
                gid[local] = got
            out_loop.append(got)
        return out_loop

    cap_face_of = {c: i for i, c in enumerate(tmesh.face_cap) if c >= 0}
    for (node, h), (ring_verts, ring_edges) in (adopted or {}).items():
        loop = tmesh.loops[cap_face_of[h]]
        keys = [owner(node, local) for local in loop]
        builder.adopt(keys, [position(*k) for k in keys], ring_verts, ring_edges)

    def merge_across(node: tuple[int, int, int], h: int):
        """The neighbour this half-strut's lateral faces merge with, or ``None``.

        Merging needs the cap to be an interface — otherwise there is no second
        half to merge with — *and* the node across it to be interior. At a
        boundary interface the other side's faces come out of a boolean, so
        there is no template loop to splice and the half-faces stand as built.
        """
        if (h if h < 3 else h - 3) not in tmesh.merged_halves:
            return None
        if (node, h) not in interfaces:
            return None
        step = steps[h]
        nb = (node[0] + step[0], node[1] + step[1], node[2] + step[2])
        return nb if nb in interior_set else None

    for idx in range(len(interior_nodes)):
        node = tuple(int(x) for x in interior_nodes[idx])
        group = 0 if groups is None else groups[node]
        gid: dict[int, int] = {}
        for fi, loop in enumerate(tmesh.loops):
            cap = tmesh.face_cap[fi]
            if cap >= 0:
                if (node, cap) in interfaces:
                    continue  # interface: dropped from both sides
                builder.face(globals_of(node, loop, gid), tmesh.face_normal[fi], group)
                continue

            half = tmesh.face_half[fi]
            across = None if half < 0 else merge_across(node, half)
            if across is not None:
                if half >= 3:
                    continue  # built once, from the outgoing side's node
                mf = tmesh.merged[fi]
                keys = [owner(node if side == 0 else across, v) for side, v in mf.loop]
                builder.face(
                    [builder.vertex(key, position(*key)) for key in keys], mf.normal, group
                )
                continue

            builder.face(globals_of(node, loop, gid), tmesh.face_normal[fi], group)

    shells, free = builder.result()
    return InteriorShells(
        shells=shells,
        stats={
            "interior_faces": builder.n_faces,
            "interior_vertices": builder.n_vertices,
            "interior_edges": builder.n_edges,
            "interior_open_edges": free,
        },
    )
