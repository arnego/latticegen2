"""Repair of two defect classes in the *input* body, found at ``import``.

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

**The second defect: a face whose boundary crosses itself at a cusp.** Where
two holes of a part touch tangentially, a planar face between them is bounded
by two circular edges that leave their shared vertex ``V`` in the same
direction. The user's file writes those circles to six decimals, so on its
faces 166, 168 and 266 they overlap by 4.2e-9 mm and cross again at ``Q``,
8.3e-5 mm from ``V``, inside both edges: the wire crosses itself, enclosing an
inverted sliver of no material, and ``BRepCheck`` calls the face invalid. A
boolean against the body widens ``V``'s tolerance to cover ``Q``, so a body
trimmed there is valid in this process -- and not once written, since STEP
cannot carry that tolerance: the lattice body of the user's run failed
``BRepCheck`` after a round trip there and nowhere else.

**The repair moves ``V`` to ``Q``**, where the two circles really meet. Both
arcs are re-bounded there on their own circles, and the third edge at ``V`` --
the line where the two hole walls touch -- is replaced by the segment from its
other, unmoved end to ``Q``. ``Q`` lies along the tangent from ``V``, so that
segment departs from both walls by the square of a micron's tenth: 3.5e-9 mm
here. No surface moves, and the material moved is the sliver's own. It acts only
on a planar face, two circles in its plane, a vertex of exactly those three
edges and three faces, and a sliver narrower than :data:`SLIVER_MAX`; and it is
checked as the other repair is, plus the one property it exists for: each
rebuilt face, written to STEP and read back, is valid on its own
(docs/algorithm.md §7.5).
"""

from __future__ import annotations

import os
import tempfile
from typing import NamedTuple

import numpy as np

from OCP.BRep import BRep_Builder, BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
from OCP.BRepBuilderAPI import BRepBuilderAPI_Copy, BRepBuilderAPI_MakeEdge
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepTools import BRepTools, BRepTools_ReShape, BRepTools_WireExplorer
from OCP.ElCLib import ElCLib
from OCP.Geom2d import Geom2d_BSplineCurve
from OCP.GeomAbs import GeomAbs_CurveType, GeomAbs_SurfaceType
from OCP.GeomAPI import GeomAPI_ProjectPointOnCurve, GeomAPI_ProjectPointOnSurf
from OCP.Precision import Precision
from OCP.ShapeBuild import ShapeBuild_Edge
from OCP.ShapeFix import ShapeFix_Edge, ShapeFix_ShapeTolerance
from OCP.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
from OCP.TColgp import TColgp_Array1OfPnt2d
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
from OCP.gp import gp_Pnt, gp_Pnt2d

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


class WireCrossing(NamedTuple):
    """One place where a face of the input file crosses its own boundary."""

    where: tuple[float, float, float]   #: where the two arcs really cross
    vertex: tuple[float, float, float]  #: the vertex the file gives them
    span: float          #: mm between the two
    width: float         #: mm, the widest the sliver they enclose gets
    repaired: bool
    reason: str          #: what was done, or why it was not


class InputRepair(NamedTuple):
    shape: TopoDS_Shape
    defects: list[InputDefect]
    max_vertex_tol: float               #: after the repair
    max_vertex_tol_at: tuple[float, float, float]
    crossings: tuple[WireCrossing, ...] = ()


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
        # The explorer already gave the edges in the sense the face runs them,
        # its wire's orientation included, so the new wire is added forward:
        # orienting it as the old one would turn a reversed wire round twice.
        builder.Add(new, nw)
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
    fail ``BRepCheck`` on their own (its faces 166, 168 and 266 do, crossing
    their own boundary at a cusp until :func:`_repair_cusps` runs) must not have a sound repair refused for that. Taken
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


# -- The second defect: a face boundary crossing itself at a cusp ------------

#: Two edges leaving their shared vertex closer than this to one direction
#: form a cusp, rad. The user's part's leave theirs 2e-8 apart.
CUSP_ANGLE = 1e-3

#: The furthest from its vertex a cusp's crossing may lie for the repair to act,
#: mm. The user's part's lie 8.3e-5 mm away, the square root of a 4.2e-9 mm
#: overlap; it grows only as that root, so this is far past any rounding.
CROSSING_MAX = GAP_TOL

#: The widest sliver between the two arcs the repair takes for a rounding of
#: the file, mm: ten times the kernel's confusion, and 240 times the 4.2e-9 mm
#: of the user's part. Wider, the holes overlap in the model itself, which the
#: CAD model must settle; the crossing is reported instead.
SLIVER_MAX = 1e-6

#: How close a computed crossing must lie to both circles, mm.
ON_CIRCLE_TOL = 1e-9

#: Samples along an arc from the vertex to the crossing, to measure the sliver
#: the two arcs enclose, and along a replaced line, to measure it against its
#: two faces.
SLIVER_SAMPLES = 9


class _Cusp(NamedTuple):
    face: TopoDS_Face
    v: TopoDS_Vertex
    a: TopoDS_Edge       #: the edge before ``v`` in the face's wire
    b: TopoDS_Edge       #: the edge after it
    q: np.ndarray        #: where their circles really cross
    ta: float            #: ``q``'s parameter on ``a``'s forward curve
    tb: float
    span: float          #: mm from ``v`` to ``q``
    width: float         #: mm, the widest the sliver between them gets


def _leaving(edge: TopoDS_Edge, v: TopoDS_Vertex) -> tuple[float, np.ndarray]:
    """``v``'s parameter on ``edge``'s forward curve, and the unit direction
    the edge leaves ``v`` in."""
    ef = TopoDS.Edge_s(edge.Oriented(FORWARD))
    c = BRepAdaptor_Curve(ef)
    at_first = TopExp.FirstVertex_s(ef).IsSame(v)
    t = c.FirstParameter() if at_first else c.LastParameter()
    d = _xyz(gp_Pnt(*c.DN(t, 1).Coord()))
    d = d / np.linalg.norm(d)
    return t, (d if at_first else -d)


def _circle_in(edge: TopoDS_Edge, pln):
    """``edge``'s circle, if it lies in the plane ``pln``, else ``None``."""
    c = BRepAdaptor_Curve(edge)
    if c.GetType() != GeomAbs_CurveType.GeomAbs_Circle:
        return None
    circ = c.Circle()
    if not circ.Axis().Direction().IsParallel(pln.Axis().Direction(), 1e-9):
        return None
    if pln.Distance(circ.Location()) > Precision.Confusion_s():
        return None
    return circ


def _in_range(edge: TopoDS_Edge, circ, p: np.ndarray) -> float | None:
    """``p``'s parameter on ``edge``'s forward circle, if strictly inside the
    edge."""
    f, l = BRep_Tool.Range_s(TopoDS.Edge_s(edge.Oriented(FORWARD)))
    t = ElCLib.InPeriod_s(ElCLib.Parameter_s(circ, gp_Pnt(*p)), f, f + 2.0 * np.pi)
    eps = Precision.PConfusion_s()
    return t if f + eps < t < l - eps else None


def _crossings(c1, c2, normal: np.ndarray) -> list[np.ndarray]:
    """Where two coplanar circles cross. The half-chord is formed as a product
    of differences so that it keeps its digits however nearly they touch."""
    p1, p2 = _xyz(c1.Location()), _xyz(c2.Location())
    r1, r2 = c1.Radius(), c2.Radius()
    d = float(np.linalg.norm(p2 - p1))
    if d == 0.0:
        return []
    h2 = (r1 + r2 - d) * (r1 + r2 + d) * (d - r1 + r2) * (d + r1 - r2) / (4.0 * d * d)
    if h2 <= 0.0:
        return []
    e = (p2 - p1) / d
    perp = np.cross(normal, e)
    along = (d * d + r1 * r1 - r2 * r2) / (2.0 * d)
    h = np.sqrt(h2)
    return [p1 + along * e + h * perp, p1 + along * e - h * perp]


def _cusp_at(face: TopoDS_Face, pln, v: TopoDS_Vertex,
             a: TopoDS_Edge, b: TopoDS_Edge) -> _Cusp | None:
    ca, cb = _circle_in(a, pln), _circle_in(b, pln)
    if ca is None or cb is None:
        return None
    va, da = _leaving(a, v)
    _vb, db = _leaving(b, v)
    if float(np.dot(da, db)) < np.cos(CUSP_ANGLE):
        return None
    vp = _vertex_xyz(v)
    best = None
    for q in _crossings(ca, cb, _xyz(gp_Pnt(*pln.Axis().Direction().Coord()))):
        if max(abs(np.linalg.norm(q - _xyz(c.Location())) - c.Radius())
               for c in (ca, cb)) > ON_CIRCLE_TOL:
            continue
        span = float(np.linalg.norm(q - vp))
        if not Precision.Confusion_s() < span <= CROSSING_MAX:
            continue
        ta, tb = _in_range(a, ca, q), _in_range(b, cb, q)
        if ta is None or tb is None:
            continue
        if best is None or span < best[1]:
            best = (q, span, ta, tb)
    if best is None:
        return None
    q, span, ta, tb = best
    # The sliver: how far a's arc from v to q strays from b's circle.
    cpb, rb = _xyz(cb.Location()), cb.Radius()
    width = 0.0
    for k in range(SLIVER_SAMPLES):
        p = _xyz(ElCLib.Value_s(va + (ta - va) * k / (SLIVER_SAMPLES - 1), ca))
        width = max(width, abs(float(np.linalg.norm(p - cpb)) - rb))
    return _Cusp(face, v, a, b, q, ta, tb, span, width)


def crossed_cusps(shape: TopoDS_Shape) -> list[_Cusp]:
    """Every vertex of a planar face where two circular edges leave in the same
    direction and cross again inside both, further than the kernel's confusion
    from it: a wire crossing itself."""
    out: list[_Cusp] = []
    for face in occ.faces(shape):
        surf = BRepAdaptor_Surface(face)
        if surf.GetType() != GeomAbs_SurfaceType.GeomAbs_Plane:
            continue
        pln = surf.Plane()
        for _w, edges in _wire_edges(face):
            n = len(edges)
            for i in range(n):
                a, b = edges[i], edges[(i + 1) % n]
                if n < 2 or a.IsSame(b):
                    continue
                v = TopExp.LastVertex_s(a, True)
                if not v.IsSame(TopExp.FirstVertex_s(b, True)):
                    continue
                if any(TopExp.FirstVertex_s(e).IsSame(TopExp.LastVertex_s(e)) for e in (a, b)):
                    continue
                cusp = _cusp_at(face, pln, v, a, b)
                if cusp is not None:
                    out.append(cusp)
    return out


def _plan_cusp(body: _Body, c: _Cusp):
    """What the repair of one cusp rebuilds, or the reason there is none."""
    if c.width > SLIVER_MAX:
        return None, f"the two arcs enclose a sliver {c.width:.1e} mm wide, not a rounding"
    edges = body.edges_at(c.v)
    if len(edges) != 3:
        return None, f"{len(edges)} edges meet at the vertex, not three"
    rest = [e for e in edges if not e.IsSame(c.a) and not e.IsSame(c.b)]
    if len(rest) != 1:
        return None, "the vertex's third edge is one of its arcs"
    line = TopoDS.Edge_s(rest[0].Oriented(FORWARD))
    if (BRep_Tool.Degenerated_s(line)
            or BRepAdaptor_Curve(line).GetType() != GeomAbs_CurveType.GeomAbs_Line):
        return None, "the vertex's third edge is not a straight line"
    w = TopExp.LastVertex_s(line) if TopExp.FirstVertex_s(line).IsSame(c.v) \
        else TopExp.FirstVertex_s(line)
    if w.IsSame(c.v):
        return None, "the line is closed"
    p = body.faces[body.face_index.FindIndex(c.face) - 1]
    fa, fb, fl = body.faces_of(c.a), body.faces_of(c.b), body.faces_of(line)
    if not (len(fa) == len(fb) == len(fl) == 2):
        return None, "an edge at the vertex is not used by two faces"
    c1 = next((f for f in fa if not f.IsSame(p)), None)
    c2 = next((f for f in fb if not f.IsSame(p)), None)
    if (c1 is None or c2 is None or c1.IsSame(c2)
            or not any(f.IsSame(c1) for f in fl) or not any(f.IsSame(c2) for f in fl)):
        return None, "the three faces at the vertex do not meet as two walls and a floor"
    for face in (p, c1, c2):
        for edge in (c.a, c.b, line):
            if BRep_Tool.IsClosed_s(edge, face):
                return None, "a seam edge is involved"
    for edge in (c.a, c.b, line):
        if not BRep_Tool.SameParameter_s(edge) or not BRep_Tool.SameRange_s(edge):
            return None, "an edge at the vertex is not same-parameter"
    # The new line, from its unmoved end to q, against both walls.
    wp = _vertex_xyz(w)
    bar = max(BRep_Tool.Tolerance_s(line), Precision.Confusion_s())
    for face in (c1, c2):
        surf = BRep_Tool.Surface_s(face)
        for k in range(SLIVER_SAMPLES):
            pt = wp + (c.q - wp) * k / (SLIVER_SAMPLES - 1)
            proj = GeomAPI_ProjectPointOnSurf(gp_Pnt(*pt), surf)
            if proj.NbPoints() < 1 or proj.LowerDistance() > bar:
                return None, (f"the line moved to the crossing would leave a wall by "
                              f"more than {bar:.1e} mm")
    return dict(c=c, p=p, c1=c1, c2=c2, line=line, w=w), None


def _rebound(edge: TopoDS_Edge, moves, faces: list[TopoDS_Face]) -> TopoDS_Edge:
    """``edge`` on its own curve and parameter curves, in the forward edge's
    sense, each end ``moves`` names -- ``(v, vq, t)`` -- moved from ``v`` to
    ``vq`` at parameter ``t``. One arc can end in a cusp at both ends."""
    ef = TopoDS.Edge_s(edge.Oriented(FORWARD))
    f, l = BRep_Tool.Range_s(ef)
    first, last = TopExp.FirstVertex_s(ef), TopExp.LastVertex_s(ef)
    for v, vq, t in moves:
        if first.IsSame(v):
            first, f = vq, t
        elif last.IsSame(v):
            last, l = vq, t
    return _sub_edge(ef, first, f, last, l, faces)


def _uv(edge: TopoDS_Edge, face: TopoDS_Face, t: float) -> gp_Pnt2d:
    c2d = BRep_Tool.CurveOnSurface_s(TopoDS.Edge_s(edge.Oriented(FORWARD)), face, 0.0, 0.0)
    return c2d.Value(t)


def _moved_line(plan, vq: TopoDS_Vertex) -> TopoDS_Edge:
    """The line from its unmoved end to ``q``, in the old line's sense, with
    parameter curves on both walls that run straight from where the old one
    left ``w`` to where each wall's arc now ends: so they join the walls'
    other edges exactly in each wall's own parameters."""
    c, line, w = plan["c"], plan["line"], plan["w"]
    from_v = TopExp.FirstVertex_s(line).IsSame(c.v)
    v1, v2 = (vq, w) if from_v else (w, vq)
    edge = BRepBuilderAPI_MakeEdge(v1, v2).Edge()
    t0, t1 = BRep_Tool.Range_s(edge)
    tw = BRep_Tool.Parameter_s(w, line)
    builder = BRep_Builder()
    for face, arc, t in ((plan["c1"], c.a, c.ta), (plan["c2"], c.b, c.tb)):
        uw, uq = _uv(line, face, tw), _uv(arc, face, t)
        poles = TColgp_Array1OfPnt2d(1, 2)
        poles.SetValue(1, uq if from_v else uw)
        poles.SetValue(2, uw if from_v else uq)
        knots = TColStd_Array1OfReal(1, 2)
        knots.SetValue(1, t0)
        knots.SetValue(2, t1)
        mults = TColStd_Array1OfInteger(1, 2)
        mults.SetValue(1, 2)
        mults.SetValue(2, 2)
        builder.UpdateEdge(edge, Geom2d_BSplineCurve(poles, knots, mults, 1), face,
                           Precision.Confusion_s())
        builder.Range(edge, face, t0, t1)
    builder.SameParameter(edge, True)
    builder.SameRange(edge, True)
    ShapeFix_Edge().FixSameParameter(edge)
    return edge


def _touched(plans) -> list[TopoDS_Face]:
    """The faces a repair of ``plans`` rebuilds: the floor first, then each
    wall once."""
    out: list[TopoDS_Face] = []
    for f in [plans[0]["p"]] + [f for pl in plans for f in (pl["c1"], pl["c2"])]:
        if not any(f.IsSame(g) for g in out):
            out.append(f)
    return out


def _apply_cusps(shape: TopoDS_Shape, plans,
                 faces: list[TopoDS_Face]) -> tuple[TopoDS_Shape, list[TopoDS_Face]]:
    """``shape`` with ``plans`` applied, and ``faces`` rebuilt, in their order."""
    builder = BRep_Builder()
    subs: list[list] = [[] for _ in faces]

    def sub(face, old, new):
        i = next(k for k, f in enumerate(faces) if f.IsSame(face))
        subs[i].append((old, lambda o, n=new: [TopoDS.Edge_s(n.Oriented(o.Orientation()))]))

    # Every move of one arc first, so an arc with a cusp at each end is
    # rebuilt once with both ends moved.
    arcs: list[tuple[TopoDS_Edge, list, list[TopoDS_Face]]] = []

    def arc_move(edge, move, arc_faces):
        hit = next((a for a in arcs if a[0].IsSame(edge)), None)
        if hit is None:
            arcs.append((edge, [move], arc_faces))
        else:
            hit[1].append(move)

    lines = []
    for plan in plans:
        c = plan["c"]
        vq = TopoDS_Vertex()
        builder.MakeVertex(vq, gp_Pnt(*c.q), Precision.Confusion_s())
        arc_move(c.a, (c.v, vq, c.ta), [plan["p"], plan["c1"]])
        arc_move(c.b, (c.v, vq, c.tb), [plan["p"], plan["c2"]])
        lines.append((plan, _moved_line(plan, vq)))
    moved = []
    for edge, moves, arc_faces in arcs:
        new = _rebound(edge, moves, arc_faces)
        for face in arc_faces:
            sub(face, edge, new)
        moved.append(new)
    for plan, line2 in lines:
        sub(plan["c1"], plan["line"], line2)
        sub(plan["c2"], plan["line"], line2)
        moved.append(line2)
    rebuilt = [_rebuild_face(f, s) for f, s in zip(faces, subs)]
    sfe = ShapeFix_Edge()
    for face in rebuilt:
        for _, edges in _wire_edges(face):
            for edge in edges:
                if any(edge.IsSame(m) for m in moved):
                    sfe.FixVertexTolerance(edge, face)
    reshape = BRepTools_ReShape()
    for old, new in zip(faces, rebuilt):
        reshape.Replace(old, new)
    return reshape.Apply(shape), rebuilt


def _survives_step(faces: list[TopoDS_Face]) -> bool:
    """Each face, written to STEP alone and read back, still valid: the file
    keeps no vertex tolerance a boolean widened to cover the crossing."""
    with tempfile.TemporaryDirectory(prefix="latticegen2-inputfix-") as d:
        for i, face in enumerate(faces):
            path = os.path.join(d, f"face{i}.step")
            occ.write_step(face, path, "inputfix-probe")
            back = occ.faces(occ.read_step(path))
            if len(back) != 1 or not BRepCheck_Analyzer(back[0]).IsValid():
                return False
    return True


def _check_cusps(before: TopoDS_Shape, after: TopoDS_Shape, rebuilt: list[TopoDS_Face],
                 valid_before: list[bool], snap: _Snapshot) -> str | None:
    """Why ``after`` is not a sound repair of the cusps, or ``None``. The face
    the cusps lie in must come out valid; each wall, valid if it was."""
    if len(occ.solids(after)) != len(occ.solids(before)):
        return "the number of solids changed"
    for i, (face, was) in enumerate(zip(rebuilt, valid_before)):
        if (was or i == 0) and not BRepCheck_Analyzer(face).IsValid():
            return "a rebuilt face is not valid on its own"
    for face in rebuilt:
        for _, edges in _wire_edges(face):
            for edge in edges:
                if BRep_Tool.Tolerance_s(edge) > snap.allowed(edge, TopoDS_Vertex()):
                    return "a rebuilt face has an edge tolerance past the bar"
                for w in (TopExp.FirstVertex_s(edge), TopExp.LastVertex_s(edge)):
                    if BRep_Tool.Tolerance_s(w) > snap.allowed(w, TopoDS_Vertex()):
                        return f"a vertex at {_round(_vertex_xyz(w))} keeps a tolerance past the bar"
    fault = _edge_use_fault(after)
    if fault:
        return fault
    for face in rebuilt:
        exact, meshed = occ.area(face), _mesh_area(face)
        if abs(exact - meshed) > AREA_MESH_REL_TOL * exact + AREA_MESH_ABS_TOL:
            return (f"a rebuilt face's area {exact:.6f} mm^2 disagrees with its "
                    f"mesh's {meshed:.6f} mm^2")
    if not _survives_step(rebuilt):
        return "a rebuilt face is not valid once written to STEP"
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


def _crossing(c: _Cusp, repaired: bool, reason: str) -> WireCrossing:
    return WireCrossing(_round(c.q), _round(_vertex_xyz(c.v)), c.span, c.width,
                        repaired, reason)


def _repair_cusps(work: TopoDS_Shape) -> tuple[TopoDS_Shape, list[WireCrossing], bool]:
    """``work`` with the cusps of each face moved to where their arcs cross,
    one face at a time, where the repair holds. All of a face's cusps are
    repaired together: the face is valid only once every one is."""
    found: list[WireCrossing] = []
    tried: set[tuple] = set()
    checked = False
    while True:
        try:
            cusps = crossed_cusps(work)
        except Exception:                 # pragma: no cover - kernel refusal
            return work, found, checked
        groups: list[list[_Cusp]] = []
        for c in cusps:
            g = next((g for g in groups if g[0].face.IsSame(c.face)), None)
            if g is None:
                groups.append([c])
            else:
                g.append(c)
        repaired = False
        for group in groups:
            key = tuple(sorted(_round(c.q) for c in group))
            if key in tried:
                continue
            tried.add(key)
            body = _Body(work)
            plans, why = [], None
            for c in group:
                try:
                    plan, why = _plan_cusp(body, c)
                except Exception as exc:  # kernel refusal: leave it as built
                    plan, why = None, f"the kernel refused to examine it ({exc})"
                if plan is None:
                    break
                plans.append(plan)
            if why is None:
                faces = _touched(plans)
                valid_before = [BRepCheck_Analyzer(f).IsValid() for f in faces]
                snap = _Snapshot(faces)
                checked = True
                try:
                    after, rebuilt = _apply_cusps(work, plans, faces)
                    why = _check_cusps(work, after, rebuilt, valid_before, snap)
                except Exception as exc:  # kernel refusal: leave it as built
                    why = f"the kernel refused the repair ({exc})"
                if why is None:
                    found += [_crossing(c, True, "vertex moved to where the arcs cross")
                              for c in group]
                    work, repaired = after, True
                    break
                snap.restore()
            found += [_crossing(c, False, why) for c in group]
        if not repaired:
            return work, found, checked


def _scan(work: TopoDS_Shape):
    """:func:`misplaced_vertices`, or ``None`` where the kernel will not answer:
    this runs on every input, so a failure here must not end a run that would
    have succeeded without it."""
    try:
        return misplaced_vertices(work)
    except Exception:                     # pragma: no cover - kernel refusal
        return None


def repair_input(shape: TopoDS_Shape) -> InputRepair:
    """``shape`` with each misplaced vertex and each crossed cusp repaired
    where the repair holds, and every one found reported."""
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
    work, crossings, cusp_checked = _repair_cusps(work)
    if checked or cusp_checked:
        # The checks mesh faces of the body in place, down to 1e-4 mm, and
        # `tessellate` would keep a mesh finer than it asks for.
        BRepTools.Clean_s(work)
    tol, at = max_vertex_tolerance(work)
    return InputRepair(work, defects, tol, at, tuple(crossings))


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
