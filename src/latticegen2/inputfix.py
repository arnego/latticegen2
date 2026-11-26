"""Repair of one defect class in the *input* body, found at ``import``.

The geometry kernel builds what a STEP file says. Where the file contradicts
itself, the body it builds carries the contradiction into every boolean this
tool runs against it, and the support and lattice near it come out wrong or not
at all. docs/algorithm.md §7.5 describes the defect, the repair and its checks;
this module implements them.

**The defect: a vertex the file places at the wrong point.** A user's file
(``TDHX_error_report_test.stp``) has an edge ``E`` whose B-spline curve starts
at x = 1919.736 while the vertex ``X`` it names sits at x = 1918.402, 1.334 mm
away, on a straight edge ``G`` that runs on past the curve's real start ``q``.
The reader builds it faithfully: ``X``'s tolerance grows to 1.334 mm to cover
both points. The face ``F`` beside ``E`` that should have run along ``G`` from
``X`` to ``q`` has no edge there, and the face ``H`` on ``E``'s other side runs
out along ``G`` to ``X`` and back, a spike of no width. Every edge is still used
by two faces, so no closure check sees it; on a planar ``F`` the body even reads
valid (``tools/make_input_defect_part.py``). But ``F``'s boundary is open in its
own parameters, two faces claim the same 1.3 mm of ``G``, and on the user's part
two support pieces there could not be trimmed: 62 mm³ of support missing, 439
ceiling points unbacked. Where ``F`` is curved the reader also patches its gap
with a *degenerated* edge -- one OCCT treats as a single point -- whose
parameter curve nevertheless runs the 1.334 mm across ``F``.

**The repair splits ``G`` where ``E``'s curve really ends.** With ``K`` the
other face using ``G``, a vertex ``V`` is made at ``q`` and ``G`` is split there
into ``A`` (from ``X``) and ``B``. ``F`` takes ``A`` beside ``E``, in place of
the degenerated edge if there is one; ``H`` takes ``B`` only, dropping its
spike; ``K`` takes ``A`` and ``B``; and ``E`` ends at ``V``. Every edge used by
two faces before is used by two after, and no surface or curve moves: the faces
are re-bounded along the curves they already had.

**It acts only where every part of that description holds**, and is checked
afterwards: the rebuilt faces must each be valid, every edge of the body must
be used by exactly two faces, the vertex tolerances it leaves must be within
:data:`REPAIR_MAX_TOL`, each rebuilt face's exact area must agree with the
area of its own mesh (:data:`AREA_MESH_REL_TOL`), and the neighbourhood
must mesh closed at some rung of the ladder the STEP round-trip check uses
(:func:`latticegen2.occ.refine_until_manifold`). A defect that fails any of
it is left as built and reported as a warning with its position, so a run on
such a file says where its output may be missing.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np

from OCP.BRep import BRep_Builder, BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Curve
from OCP.BRepBuilderAPI import BRepBuilderAPI_Copy, BRepBuilderAPI_MakeEdge
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepTools import BRepTools, BRepTools_ReShape, BRepTools_WireExplorer
from OCP.GeomAPI import GeomAPI_ProjectPointOnCurve
from OCP.Precision import Precision
from OCP.ShapeBuild import ShapeBuild_Edge
from OCP.ShapeFix import ShapeFix_Edge, ShapeFix_ShapeTolerance
from OCP.TopAbs import TopAbs_Orientation, TopAbs_ShapeEnum
from OCP.TopExp import TopExp
from OCP.TopTools import (
    TopTools_IndexedDataMapOfShapeListOfShape,
    TopTools_IndexedMapOfShape,
)
from OCP.TopoDS import (
    TopoDS,
    TopoDS_Edge,
    TopoDS_Face,
    TopoDS_Iterator,
    TopoDS_Shape,
    TopoDS_Vertex,
    TopoDS_Wire,
)
from OCP.gp import gp_Pnt

from . import occ

FORWARD = TopAbs_Orientation.TopAbs_FORWARD
REVERSED = TopAbs_Orientation.TopAbs_REVERSED

#: An edge whose 3D curve ends further than this from the vertex it names is
#: a contradiction in the file, mm. A vertex's tolerance covers the gap, and
#: ordinary CAD exports leave gaps up to ~1e-3 mm (9.6e-4 is the largest on the
#: user's part besides the defect, at 1.334 mm). Ten microns is ten times that
#: and forty times under the smallest legal strut.
GAP_TOL = 1e-2

#: How closely the pieces of the description must agree for the repair to
#: act, mm: the curve's real end against the edge it splits, that edge's own
#: end against the vertex, and a degenerated edge's path against the edge it
#: runs along. On the user's part they agree to 1e-9 mm.
MATCH_TOL = 1e-4

#: The largest vertex or edge tolerance the repair may leave, mm.
REPAIR_MAX_TOL = occ.MICRO_CAP_TOL

#: A rebuilt face's exact area against the area of its own mesh at
#: :data:`AREA_MESH_DEFLECTION`: relative and absolute parts, mm². This stands
#: where a volume bar would, because the defective body has none to compare
#: with. Its open parameter-space wire makes OCCT's boundary integral count a
#: strip that is not there: face 423 of the user's part reads 3.8605 mm², its
#: 1 mm fillet radius times 0.2839 rad times 13.598 mm, out to the misplaced
#: vertex, while its mesh and the repaired face read 3.4817, the same times the
#: 12.264 mm it really spans -- docs/algorithm.md §7's lesson again (G19): a
#: quantity OCCT does not define on the shape cannot be a bar. The mesh is the
#: independent reading, and it agrees with the repaired face to 2e-7.
AREA_MESH_REL_TOL = 1e-3
AREA_MESH_ABS_TOL = 1e-6
AREA_MESH_DEFLECTION = 1e-4

#: Samples along a degenerated edge's path, to compare it with the edge it
#: runs along.
PATH_SAMPLES = 9


class InputDefect(NamedTuple):
    """One place where the input file contradicts itself."""

    where: tuple[float, float, float]   #: where the edge's curve really ends
    span: float          #: mm the curve ends from the vertex it names
    repaired: bool
    reason: str          #: what was done, or why it was not
    vertex: tuple[float, float, float] = (0.0, 0.0, 0.0)  #: the vertex it names


class InputRepair(NamedTuple):
    shape: TopoDS_Shape
    defects: list[InputDefect]
    max_vertex_tol: float               #: after the repair
    max_vertex_tol_at: tuple[float, float, float]


def _xyz(p: gp_Pnt) -> np.ndarray:
    return np.array([p.X(), p.Y(), p.Z()])


def _round(p: np.ndarray) -> tuple[float, float, float]:
    return (round(float(p[0]), 3), round(float(p[1]), 3), round(float(p[2]), 3))


def _vertex_xyz(v: TopoDS_Vertex) -> np.ndarray:
    return _xyz(BRep_Tool.Pnt_s(v))


def _ends(edge: TopoDS_Edge) -> tuple[np.ndarray, np.ndarray]:
    """The 3D curve's start and end, in ``edge``'s own orientation."""
    c = BRepAdaptor_Curve(edge)
    a, b = _xyz(c.Value(c.FirstParameter())), _xyz(c.Value(c.LastParameter()))
    return (a, b) if edge.Orientation() != REVERSED else (b, a)


def _path(edge: TopoDS_Edge, face: TopoDS_Face) -> np.ndarray:
    """The 3D points ``edge``'s parameter curve sweeps on ``face``, in the
    edge's own orientation."""
    c = BRepAdaptor_Curve(edge, face)
    f, l = c.FirstParameter(), c.LastParameter()
    pts = np.array([_xyz(c.Value(f + (l - f) * k / (PATH_SAMPLES - 1)))
                    for k in range(PATH_SAMPLES)])
    return pts if edge.Orientation() != REVERSED else pts[::-1]


def _wire_edges(face: TopoDS_Face) -> list[tuple[TopoDS_Wire, list[TopoDS_Edge]]]:
    """Each wire of ``face`` with its edges in connection order, oriented as
    the face traverses them."""
    fwd = TopoDS.Face_s(face.Oriented(FORWARD))
    out = []
    it = TopoDS_Iterator(fwd)
    while it.More():
        w = TopoDS.Wire_s(it.Value())
        ex = BRepTools_WireExplorer(w, fwd)
        edges = []
        while ex.More():
            edges.append(TopoDS.Edge_s(ex.Current()))
            ex.Next()
        out.append((w, edges))
        it.Next()
    return out


def _rebuild_face(face: TopoDS_Face, subs) -> TopoDS_Face:
    """``face`` with each edge ``subs`` names replaced.

    ``subs`` maps an edge to a function from the edge as this face traverses
    it to the oriented edges that take its place, in traversal order.
    """
    builder = BRep_Builder()
    fwd = TopoDS.Face_s(face.Oriented(FORWARD))
    new = TopoDS.Face_s(fwd.EmptyCopied())
    for wire, edges in _wire_edges(face):
        if len(edges) != wire.NbChildren():
            raise RuntimeError("a wire of a face to rebuild could not be walked in order")
        nw = TopoDS_Wire()
        builder.MakeWire(nw)
        for e in edges:
            repl = next((fn for key, fn in subs if key.IsSame(e)), None)
            for r in (repl(e) if repl is not None else [e]):
                builder.Add(nw, r)
        nw.Closed(True)
        builder.Add(new, nw.Oriented(wire.Orientation()))
    return TopoDS.Face_s(new.Oriented(face.Orientation()))


def _sub_edge(g: TopoDS_Edge, v1: TopoDS_Vertex, t1: float,
              v2: TopoDS_Vertex, t2: float, faces: list[TopoDS_Face]) -> TopoDS_Edge:
    """The piece of ``g``'s curve from ``t1`` to ``t2`` (``t1 < t2``), with
    ``g``'s parameter curves on ``faces``."""
    curve = BRep_Tool.Curve_s(g, 0.0, 0.0)
    e = BRepBuilderAPI_MakeEdge(curve, v1, v2, t1, t2).Edge()
    builder = BRep_Builder()
    tol = BRep_Tool.Tolerance_s(g)
    for f in faces:
        c2d = BRep_Tool.CurveOnSurface_s(g, f, 0.0, 0.0)
        builder.UpdateEdge(e, c2d, f, tol)
        builder.Range(e, f, t1, t2)
    builder.SameParameter(e, True)
    builder.SameRange(e, True)
    return e


def _orient_from(edge: TopoDS_Edge, start: np.ndarray) -> TopoDS_Edge:
    """``edge`` oriented to start at the end nearer ``start``."""
    a = _vertex_xyz(TopExp.FirstVertex_s(edge))
    b = _vertex_xyz(TopExp.LastVertex_s(edge))
    o = FORWARD if np.linalg.norm(a - start) <= np.linalg.norm(b - start) else REVERSED
    return TopoDS.Edge_s(edge.Oriented(o))


class _Body:
    """Adjacency of one shape, built once."""

    def __init__(self, shape: TopoDS_Shape):
        self.faces = occ.faces(shape)
        self.face_index = TopTools_IndexedMapOfShape()
        for f in self.faces:
            self.face_index.Add(f)
        self.edge_faces = TopTools_IndexedDataMapOfShapeListOfShape()
        TopExp.MapShapesAndAncestors_s(shape, TopAbs_ShapeEnum.TopAbs_EDGE,
                                       TopAbs_ShapeEnum.TopAbs_FACE, self.edge_faces)
        self.vertex_edges = TopTools_IndexedDataMapOfShapeListOfShape()
        TopExp.MapShapesAndAncestors_s(shape, TopAbs_ShapeEnum.TopAbs_VERTEX,
                                       TopAbs_ShapeEnum.TopAbs_EDGE, self.vertex_edges)

    def faces_of(self, edge: TopoDS_Shape) -> list[TopoDS_Face]:
        out: list[TopoDS_Face] = []
        for f in self.edge_faces.FindFromKey(edge):
            f = TopoDS.Face_s(f)
            if not any(g.IsSame(f) for g in out):
                out.append(self.faces[self.face_index.FindIndex(f) - 1])
        return out

    def edges_at(self, vertex: TopoDS_Vertex) -> list[TopoDS_Edge]:
        out: list[TopoDS_Edge] = []
        for e in self.vertex_edges.FindFromKey(vertex):
            e = TopoDS.Edge_s(e)
            if not any(g.IsSame(e) for g in out):
                out.append(e)
        return out


def misplaced_vertices(shape: TopoDS_Shape):
    """``(edge, vertex, end, gap)`` for every edge whose 3D curve ends more
    than :data:`GAP_TOL` from the vertex the edge names there."""
    body = _Body(shape)
    out = []
    for i in range(1, body.edge_faces.Extent() + 1):
        e = TopoDS.Edge_s(body.edge_faces.FindKey(i).Oriented(FORWARD))
        if BRep_Tool.Degenerated_s(e) or BRep_Tool.Curve_s(e, 0.0, 0.0) is None:
            continue
        ends = _ends(e)
        for v, p in ((TopExp.FirstVertex_s(e), ends[0]), (TopExp.LastVertex_s(e), ends[1])):
            gap = float(np.linalg.norm(p - _vertex_xyz(v)))
            if gap > GAP_TOL:
                out.append((e, v, p, gap))
    return out


def _wire_with(face: TopoDS_Face, edge: TopoDS_Shape) -> list[TopoDS_Edge] | None:
    hits = [edges for _, edges in _wire_edges(face) if any(e.IsSame(edge) for e in edges)]
    return hits[0] if len(hits) == 1 else None


def _neighbour_at(edges: list[TopoDS_Edge], edge: TopoDS_Shape, x: TopoDS_Vertex):
    """The edge before or after ``edge`` in the wire ``edges``, on the side
    where the wire passes through ``x``."""
    i = next(k for k, e in enumerate(edges) if e.IsSame(edge))
    here = edges[i]
    if TopExp.FirstVertex_s(here, True).IsSame(x):
        return edges[i - 1]
    if TopExp.LastVertex_s(here, True).IsSame(x):
        return edges[(i + 1) % len(edges)]
    return None


def _plan(body: _Body, e: TopoDS_Edge, x: TopoDS_Vertex, q: np.ndarray):
    """The split that repairs ``e``'s misplaced vertex ``x``, or the reason
    there is none."""
    xp = _vertex_xyz(x)
    e_faces = body.faces_of(e)
    if len(e_faces) != 2:
        return None, "the edge is not used by two faces"
    # G: the edge at X whose curve runs on through q.
    found = []
    for g in body.edges_at(x):
        g = TopoDS.Edge_s(g.Oriented(FORWARD))
        if g.IsSame(e) or BRep_Tool.Degenerated_s(g):
            continue
        curve = BRep_Tool.Curve_s(g, 0.0, 0.0)
        if curve is None:
            continue
        x_first = TopExp.FirstVertex_s(g).IsSame(x)
        if x_first == TopExp.LastVertex_s(g).IsSame(x):
            continue                                  # closed at X, or not at X
        g_ends = _ends(g)
        if np.linalg.norm((g_ends[0] if x_first else g_ends[1]) - xp) > MATCH_TOL:
            continue                                  # G itself does not end at X
        adaptor = BRepAdaptor_Curve(g)
        f0, l0 = adaptor.FirstParameter(), adaptor.LastParameter()
        proj = GeomAPI_ProjectPointOnCurve(gp_Pnt(*q), curve)
        if proj.NbPoints() < 1 or proj.LowerDistance() > MATCH_TOL:
            continue
        tq = proj.LowerDistanceParameter()
        eps = 1e-9 * max(1.0, abs(l0 - f0))
        if not f0 + eps < tq < l0 - eps:
            continue
        found.append((g, curve, tq, f0, l0, x_first))
    if len(found) != 1:
        return None, (f"{len(found)} edges at the vertex run on through the curve's "
                      f"end, not one")
    g, curve, tq, f0, l0, x_first = found[0]
    if not BRep_Tool.SameParameter_s(g) or not BRep_Tool.SameRange_s(g):
        return None, "the edge to split is not same-parameter"
    g_faces = body.faces_of(g)
    if len(g_faces) != 2:
        return None, "the edge to split is not used by two faces"
    shared = [h for h in e_faces if any(k.IsSame(h) for k in g_faces)]
    if len(shared) != 1:
        return None, "not exactly one face uses both the edge and the edge to split"
    h = shared[0]
    f = next(k for k in e_faces if not k.IsSame(h))
    k = next(k for k in g_faces if not k.IsSame(h))
    if f.IsSame(k):
        return None, "one face lies on both sides"
    for face in (f, h, k):
        for edge in (g, e):
            if BRep_Tool.IsClosed_s(edge, face):
                return None, "a seam edge is involved"
    h_wire = _wire_with(h, e)
    beside = None if h_wire is None else _neighbour_at(h_wire, e, x)
    if beside is None or not beside.IsSame(g):
        return None, "the edge and the edge to split do not meet at the vertex"
    f_wire = _wire_with(f, e)
    d = None if f_wire is None else _neighbour_at(f_wire, e, x)
    if d is None:
        return None, "the edge's other face does not pass through the vertex"
    if BRep_Tool.Degenerated_s(d):
        # The reader's patch over the gap: it must sweep X -> q along G.
        path = _path(d, f)
        ends_match = min(np.linalg.norm(path[0] - xp) + np.linalg.norm(path[-1] - q),
                         np.linalg.norm(path[-1] - xp) + np.linalg.norm(path[0] - q))
        if ends_match > 2 * MATCH_TOL:
            return None, "a degenerated edge at the vertex sweeps somewhere else"
        if any(GeomAPI_ProjectPointOnCurve(gp_Pnt(*p), curve).LowerDistance() > MATCH_TOL
               for p in path):
            return None, ("a degenerated edge at the vertex does not run along the "
                          "edge to split")
    else:
        d = None
    return dict(x=x, q=q, e=e, g=g, tq=tq, f0=f0, l0=l0, x_first=x_first,
                d=d, f=f, h=h, k=k), None


def _apply(shape: TopoDS_Shape, plan) -> tuple[TopoDS_Shape, list[TopoDS_Face]]:
    x, q, e, g, tq = plan["x"], plan["q"], plan["e"], plan["g"], plan["tq"]
    d, f, h, k = plan["d"], plan["f"], plan["h"], plan["k"]
    curve = BRep_Tool.Curve_s(g, 0.0, 0.0)
    builder = BRep_Builder()
    vq = TopoDS_Vertex()
    builder.MakeVertex(vq, curve.Value(tq), Precision.Confusion_s())
    other = TopExp.LastVertex_s(g) if plan["x_first"] else TopExp.FirstVertex_s(g)
    if plan["x_first"]:
        a = _sub_edge(g, x, plan["f0"], vq, tq, [k])
        b = _sub_edge(g, vq, tq, other, plan["l0"], [k, h])
        pieces = [a, b]                    # in g's own sense
    else:
        b = _sub_edge(g, other, plan["f0"], vq, tq, [k, h])
        a = _sub_edge(g, vq, tq, x, plan["l0"], [k])
        pieces = [b, a]
    sfe = ShapeFix_Edge()
    sfe.FixAddPCurve(a, f, False, MATCH_TOL)
    if BRep_Tool.CurveOnSurface_s(a, f, 0.0, 0.0) is None:
        raise RuntimeError("no parameter curve for the split piece on its new face")
    sfe.FixSameParameter(a)
    # E ends at the new vertex instead of X.
    first_is_x = TopExp.FirstVertex_s(e).IsSame(x)
    e2 = ShapeBuild_Edge().CopyReplaceVertices(
        e, vq if first_is_x else TopoDS_Vertex(), TopoDS_Vertex() if first_is_x else vq)
    xp = _vertex_xyz(x)

    def same_sense(old):
        return [TopoDS.Edge_s(e2.Oriented(old.Orientation()))]

    def g_pieces(only_b: bool):
        def fn(old):
            ps = [b] if only_b else pieces
            if old.Orientation() == REVERSED:
                ps = ps[::-1]
            return [TopoDS.Edge_s(p.Oriented(old.Orientation())) for p in ps]
        return fn

    def e_with_a(old):
        # F's wire passes through X beside E: A now runs between X and q there.
        new = TopoDS.Edge_s(e2.Oriented(old.Orientation()))
        if TopExp.FirstVertex_s(old, True).IsSame(x):
            return [_orient_from(a, xp), new]
        return [new, _orient_from(a, q)]

    f_subs = [(e, e_with_a)] + ([(d, lambda old: [])] if d is not None else [])
    new_f = _rebuild_face(f, f_subs)
    new_h = _rebuild_face(h, [(g, g_pieces(True)), (e, same_sense)])
    new_k = _rebuild_face(k, [(g, g_pieces(False))])
    # X kept the tolerance that covered q; it needs only what is left at it.
    ShapeFix_ShapeTolerance().SetTolerance(x, Precision.Confusion_s(),
                                           TopAbs_ShapeEnum.TopAbs_VERTEX)
    for face in (new_f, new_h, new_k):
        for _, edges in _wire_edges(face):
            for edge in edges:
                ends = (TopExp.FirstVertex_s(edge), TopExp.LastVertex_s(edge))
                if any(w.IsSame(x) or w.IsSame(vq) for w in ends):
                    sfe.FixVertexTolerance(edge, face)
    reshape = BRepTools_ReShape()
    for old, new in ((f, new_f), (h, new_h), (k, new_k)):
        reshape.Replace(old, new)
    return reshape.Apply(shape), [new_f, new_h, new_k]


def _edge_use_fault(shape: TopoDS_Shape) -> str | None:
    emap = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(shape, TopAbs_ShapeEnum.TopAbs_EDGE,
                                   TopAbs_ShapeEnum.TopAbs_FACE, emap)
    for i in range(1, emap.Extent() + 1):
        e = TopoDS.Edge_s(emap.FindKey(i))
        if BRep_Tool.Degenerated_s(e):
            continue
        faces: list = []
        for f in emap.FindFromIndex(i):
            if not any(g.IsSame(f) for g in faces):
                faces.append(f)
        uses = 0
        for f in faces:
            uses += 2 if BRep_Tool.IsClosed_s(e, TopoDS.Face_s(f)) else 1
        if uses != 2:
            p = _round(_ends(e)[0])
            return f"an edge at {p} is used {uses} time(s)"
    return None


def _mesh_area(face: TopoDS_Face) -> float:
    """The area of ``face``'s own mesh, on a copy so nothing is left on it."""
    copy = TopoDS.Face_s(BRepBuilderAPI_Copy(face, True, False).Shape())
    occ.mesh_shape(copy, AREA_MESH_DEFLECTION, 0.1)
    tri = occ.face_triangulation(copy)
    if tri is None:
        return 0.0
    v, t = tri
    return float(0.5 * np.linalg.norm(
        np.cross(v[t[:, 1]] - v[t[:, 0]], v[t[:, 2]] - v[t[:, 0]]), axis=1).sum())


class _Snapshot:
    """What the three faces a repair rebuilds were before it: each vertex's
    and edge's tolerance, and whether each face was valid on its own.

    The checks hold the repair to what it changes, not to the file: an input
    whose faces already carry a 4.9e-3 mm vertex (the user's part has them) or
    fail ``BRepCheck`` on their own (its faces 166, 168 and 266 do, where arcs
    meet tangentially) must not have a sound repair refused for that. Taken
    before :func:`_apply`, which changes tolerances in place, and used to put
    them back if the repair is refused."""

    def __init__(self, faces: list[TopoDS_Face]):
        self.tolerance: list[tuple[TopoDS_Shape, float]] = []
        self.index = TopTools_IndexedMapOfShape()
        for kind in (TopAbs_ShapeEnum.TopAbs_VERTEX, TopAbs_ShapeEnum.TopAbs_EDGE):
            for f in faces:
                TopExp.MapShapes_s(f, kind, self.index)
        for i in range(1, self.index.Extent() + 1):
            sub = self.index.FindKey(i)
            tol = (BRep_Tool.Tolerance_s(TopoDS.Vertex_s(sub))
                   if sub.ShapeType() == TopAbs_ShapeEnum.TopAbs_VERTEX
                   else BRep_Tool.Tolerance_s(TopoDS.Edge_s(sub)))
            self.tolerance.append((sub, tol))
        self.valid = all(BRepCheck_Analyzer(f).IsValid() for f in faces)

    def allowed(self, sub: TopoDS_Shape, fresh: TopoDS_Vertex) -> float:
        """The largest tolerance ``sub`` may keep: :data:`REPAIR_MAX_TOL`, or
        what it already had if more -- except ``fresh``, the vertex whose
        inflated tolerance is the defect."""
        i = self.index.FindIndex(sub)
        if i == 0 or sub.IsSame(fresh):
            return REPAIR_MAX_TOL
        return max(REPAIR_MAX_TOL, self.tolerance[i - 1][1])

    def restore(self) -> None:
        fix = ShapeFix_ShapeTolerance()
        for sub, tol in self.tolerance:
            kind = sub.ShapeType()
            fix.SetTolerance(sub, tol, kind)


def _check(before: TopoDS_Shape, after: TopoDS_Shape, rebuilt: list[TopoDS_Face],
           snap: _Snapshot, x: TopoDS_Vertex) -> str | None:
    """Why ``after`` is not a sound repair of ``before``, or ``None``."""
    if len(occ.solids(after)) != len(occ.solids(before)):
        return "the number of solids changed"
    if snap.valid and not all(BRepCheck_Analyzer(f).IsValid() for f in rebuilt):
        return "a rebuilt face is not valid on its own"
    for face in rebuilt:
        for _, edges in _wire_edges(face):
            for edge in edges:
                if BRep_Tool.Tolerance_s(edge) > snap.allowed(edge, x):
                    return "a rebuilt face has an edge tolerance past the bar"
                for w in (TopExp.FirstVertex_s(edge), TopExp.LastVertex_s(edge)):
                    if BRep_Tool.Tolerance_s(w) > snap.allowed(w, x):
                        return f"a vertex at {_round(_vertex_xyz(w))} keeps a tolerance past the bar"
    fault = _edge_use_fault(after)
    if fault:
        return fault
    for face in rebuilt:
        exact, meshed = occ.area(face), _mesh_area(face)
        if abs(exact - meshed) > AREA_MESH_REL_TOL * exact + AREA_MESH_ABS_TOL:
            return (f"a rebuilt face's area {exact:.6f} mm^2 disagrees with its "
                    f"mesh's {meshed:.6f} mm^2")
    index = TopTools_IndexedMapOfShape()
    for fc in occ.faces(after):
        index.Add(fc)
    core = frozenset(index.FindIndex(fc) - 1 for fc in rebuilt)
    if -1 in core:
        return "a rebuilt face is not in the repaired body"
    ref = occ.refine_until_manifold(after, core)
    if not ref.resolved:
        return f"the repaired faces do not mesh closed ({ref.reason})"
    return None


def _scan(work: TopoDS_Shape):
    """:func:`misplaced_vertices`, or ``None`` where the kernel will not answer:
    this runs on every input, so a failure here must not end a run that would
    have succeeded without it."""
    try:
        return misplaced_vertices(work)
    except Exception:                     # pragma: no cover - kernel refusal
        return None


def repair_input(shape: TopoDS_Shape) -> InputRepair:
    """``shape`` with each misplaced vertex repaired where the repair holds,
    and every one found reported."""
    defects: list[InputDefect] = []
    work = shape
    tried: set[tuple] = set()
    checked = False
    # A repair rebuilds faces, so the body is scanned again after each one; a
    # declined candidate changes nothing, so the scan it came from still holds.
    while True:
        repaired = False
        for e, x, q, gap in _scan(work) or []:
            where = _round(q)
            at = _round(_vertex_xyz(x))
            if (where, at) in tried:
                continue
            tried.add((where, at))
            try:
                plan, why = _plan(_Body(work), e, x, q)
            except Exception as exc:      # kernel refusal: leave it as built
                plan, why = None, f"the kernel refused to examine it ({exc})"
            if plan is None:
                defects.append(InputDefect(where, gap, False, why, at))
                continue
            snap = _Snapshot([plan["f"], plan["h"], plan["k"]])
            checked = True
            try:
                after, rebuilt = _apply(work, plan)
                why = _check(work, after, rebuilt, snap, x)
            except Exception as exc:      # kernel refusal: leave it as built
                why = f"the kernel refused the repair ({exc})"
            if why is not None:
                snap.restore()
                defects.append(InputDefect(where, gap, False, why, at))
                continue
            defects.append(InputDefect(where, gap, True, "split where the curve ends", at))
            work, repaired = after, True
            break
        if not repaired:
            break
    if checked:
        # The checks mesh faces of the body in place, down to 1e-4 mm, and
        # `tessellate` would keep a mesh finer than it asks for.
        BRepTools.Clean_s(work)
    tol, at = max_vertex_tolerance(work)
    return InputRepair(work, defects, tol, at)


def max_vertex_tolerance(shape: TopoDS_Shape) -> tuple[float, tuple[float, float, float]]:
    vm = TopTools_IndexedMapOfShape()
    TopExp.MapShapes_s(shape, TopAbs_ShapeEnum.TopAbs_VERTEX, vm)
    best, at = 0.0, (0.0, 0.0, 0.0)
    for i in range(1, vm.Extent() + 1):
        v = TopoDS.Vertex_s(vm.FindKey(i))
        t = BRep_Tool.Tolerance_s(v)
        if t > best:
            best, at = t, _round(_vertex_xyz(v))
    return best, at
