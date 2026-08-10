"""Thin, typed helpers over the OCCT kernel (via the OCP bindings).

Everything in this module is a small wrapper whose only job is to keep OCCT's
C++-shaped API (static ``*_s`` methods, out-parameters, 1-based arrays) out of
the rest of the codebase. No lattice logic lives here.

Why OCP rather than gmsh: the fuse-free architecture needs sewing, shared
topology construction, located instances and exact validity checking, none of
which gmsh's scripting API exposes — see
docs/research/perf-rearchitecture-proposal.md §2 "Insight 4".
"""

from __future__ import annotations

import os
from typing import Iterable, Iterator

import numpy as np

from OCP.BRep import BRep_Builder, BRep_Tool
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepBuilderAPI import (
    BRepBuilderAPI_MakeFace,
    BRepBuilderAPI_MakePolygon,
    BRepBuilderAPI_MakeSolid,
    BRepBuilderAPI_Sewing,
)
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.BRepPrimAPI import BRepPrimAPI_MakePrism
from OCP.Bnd import Bnd_Box
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_ReturnStatus
from OCP.Interface import Interface_Static
from OCP.STEPControl import STEPControl_Reader, STEPControl_StepModelType, STEPControl_Writer
from OCP.TopAbs import TopAbs_ShapeEnum
from OCP.TopExp import TopExp_Explorer
from OCP.TopLoc import TopLoc_Location
from OCP.TopoDS import TopoDS, TopoDS_Compound, TopoDS_Face, TopoDS_Shape, TopoDS_Shell
from OCP.gp import gp_Pnt, gp_Trsf, gp_Vec

from .errors import InputGeometryError, OutputError

# --- Shape traversal --------------------------------------------------------


def _explore(shape: TopoDS_Shape, kind: TopAbs_ShapeEnum) -> Iterator[TopoDS_Shape]:
    exp = TopExp_Explorer(shape, kind)
    while exp.More():
        yield exp.Current()
        exp.Next()


def solids(shape: TopoDS_Shape) -> list[TopoDS_Shape]:
    """Every solid in ``shape`` (the shape itself if it is already one)."""
    if shape.ShapeType() == TopAbs_ShapeEnum.TopAbs_SOLID:
        return [shape]
    return [TopoDS.Solid_s(s) for s in _explore(shape, TopAbs_ShapeEnum.TopAbs_SOLID)]


def faces(shape: TopoDS_Shape) -> list[TopoDS_Face]:
    """Every face in ``shape``."""
    return [TopoDS.Face_s(f) for f in _explore(shape, TopAbs_ShapeEnum.TopAbs_FACE)]


def compound(shapes: Iterable[TopoDS_Shape]) -> TopoDS_Compound:
    """Pack shapes into a single ``TopoDS_Compound`` (no geometric work)."""
    comp = TopoDS_Compound()
    builder = BRep_Builder()
    builder.MakeCompound(comp)
    for s in shapes:
        builder.Add(comp, s)
    return comp


# --- Measurement ------------------------------------------------------------


def volume(shape: TopoDS_Shape) -> float:
    """Exact volume in mm³ (OCCT Gauss quadrature, unit density)."""
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, props)
    return props.Mass()


def area(shape: TopoDS_Shape) -> float:
    """Exact surface area in mm² of a face or shell (OCCT Gauss quadrature).

    On a trimmed face this is the *trimmed* area — the tight quantity the
    mesh-coverage gate needs, as opposed to the deliberately conservative
    bounding box (docs/algorithm.md §11.3).
    """
    props = GProp_GProps()
    BRepGProp.SurfaceProperties_s(shape, props)
    return props.Mass()


def bounding_box(shape: TopoDS_Shape, tol: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """Axis-aligned ``(lo, hi)`` bounds.

    OCCT's box is a deliberate **over-estimate** (control-point hull for
    B-splines, untrimmed UV rectangle for planes), so it is only ever sound to
    use as an upper bound — never to test whether geometry *reaches* something
    (docs/algorithm.md §11.3).
    """
    box = Bnd_Box()
    box.SetGap(tol)
    # useTriangulation=False: compute from the exact geometry. A triangulation,
    # if one happens to be present, gives a box that can sit slightly *inside*
    # the true surface, and this box is what bounds candidate node enumeration.
    BRepBndLib.Add_s(shape, box, False)
    if box.IsVoid():
        raise InputGeometryError("Shape has an empty bounding box (no geometry).")
    xmin, ymin, zmin, xmax, ymax, zmax = box.Get()
    return np.array([xmin, ymin, zmin]), np.array([xmax, ymax, zmax])


def is_valid(shape: TopoDS_Shape) -> bool:
    """OCCT's exact B-rep validity check.

    This is the gate the Julia/gmsh implementation could not express at all —
    docs/algorithm.md §11.1 had to settle for indirect, mesh-based evidence
    because ``BRepCheck_Analyzer`` is not reachable through gmsh's scripting
    API.
    """
    return BRepCheck_Analyzer(shape).IsValid()


# --- Construction -----------------------------------------------------------


def _pnt(p: np.ndarray) -> gp_Pnt:
    return gp_Pnt(float(p[0]), float(p[1]), float(p[2]))


def polygon_face(verts: np.ndarray) -> TopoDS_Face:
    """Planar face bounded by the closed polygon through ``verts`` (``(N,3)``)."""
    poly = BRepBuilderAPI_MakePolygon()
    for p in verts:
        poly.Add(_pnt(p))
    poly.Close()
    return BRepBuilderAPI_MakeFace(poly.Wire()).Face()


def prism(base: TopoDS_Face, direction: np.ndarray) -> TopoDS_Shape:
    """Extrude ``base`` along ``direction`` (length = ``|direction|``)."""
    vec = gp_Vec(float(direction[0]), float(direction[1]), float(direction[2]))
    return BRepPrimAPI_MakePrism(base, vec).Shape()


def translation(offset: np.ndarray) -> TopLoc_Location:
    """A ``TopLoc_Location`` pure translation, for O(1) instancing.

    ``shape.Moved(loc)`` re-uses the underlying geometry and only records the
    transform, so instancing a template costs no geometric work.
    """
    trsf = gp_Trsf()
    trsf.SetTranslation(gp_Vec(float(offset[0]), float(offset[1]), float(offset[2])))
    return TopLoc_Location(trsf)


def sew(shapes: Iterable[TopoDS_Shape], tolerance: float) -> TopoDS_Shape:
    """Stitch shapes into shells along coincident free edges.

    Feeding already-built shells (rather than loose faces) matters: sewing only
    works on *free* edges, so a shell whose interior edges are already shared
    contributes only its holes to the workload. That is what keeps the cost
    proportional to the boundary region instead of to the whole lattice.
    """
    sewing = BRepBuilderAPI_Sewing(tolerance, True, True, True, False)
    for s in shapes:
        sewing.Add(s)
    sewing.Perform()
    return sewing.SewedShape()


def make_solid(shell: TopoDS_Shell) -> TopoDS_Shape:
    """Wrap a closed shell as a solid."""
    return BRepBuilderAPI_MakeSolid(shell).Solid()


def shells(shape: TopoDS_Shape) -> list[TopoDS_Shell]:
    """Every shell in ``shape``."""
    return [TopoDS.Shell_s(s) for s in _explore(shape, TopAbs_ShapeEnum.TopAbs_SHELL)]


# --- Meshing ----------------------------------------------------------------


def mesh_shape(shape: TopoDS_Shape, deflection: float, angle: float = 0.35) -> None:
    """Tessellate ``shape`` in place with a chordal-deviation bound.

    ``deflection`` is OCCT's linear (chordal) deviation target, so it *is* the
    ``d`` of docs/algorithm.md §5 rather than a proxy for it — a planar face
    meshes to a couple of exact triangles regardless of size, while curvature
    drives refinement automatically. That is why this needs no separate maximum
    element size the way gmsh's curvature-based sizing did.
    """
    BRepMesh_IncrementalMesh(shape, deflection, False, angle, True)


def is_planar(face: TopoDS_Face) -> bool:
    """Is this face a plane? Planes are tessellated exactly, with zero deviation."""
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.GeomAbs import GeomAbs_SurfaceType

    return BRepAdaptor_Surface(face).GetType() == GeomAbs_SurfaceType.GeomAbs_Plane


def face_uv_triangulation(face: TopoDS_Face):
    """``(points, uv, tris, evaluate)`` for a meshed face, or ``None``.

    ``points`` are world-space triangulation nodes, ``uv`` their surface
    parameters, ``tris`` the 0-based index triples, and ``evaluate(uv_array)``
    maps parameters to world-space points on the *true* surface. Together these
    are what :func:`latticegen2.classify.measure_face_deviation` needs to
    compare the mesh against the geometry it approximates.
    """
    from OCP.BRepAdaptor import BRepAdaptor_Surface

    loc = TopLoc_Location()
    tri = BRep_Tool.Triangulation_s(face, loc)
    if tri is None or not tri.HasUVNodes():
        return None
    trsf = loc.Transformation()
    n = tri.NbNodes()
    pts = np.empty((n, 3))
    uvs = np.empty((n, 2))
    for i in range(1, n + 1):
        p = tri.Node(i).Transformed(trsf)
        pts[i - 1] = (p.X(), p.Y(), p.Z())
        uv = tri.UVNode(i)
        uvs[i - 1] = (uv.X(), uv.Y())
    m = tri.NbTriangles()
    tris = np.empty((m, 3), dtype=np.int64)
    for i in range(1, m + 1):
        t = tri.Triangle(i)
        tris[i - 1] = (t.Value(1) - 1, t.Value(2) - 1, t.Value(3) - 1)

    surf = BRepAdaptor_Surface(face)

    def evaluate(uv_array: np.ndarray) -> np.ndarray:
        out = np.empty((len(uv_array), 3))
        for i, (u, v) in enumerate(uv_array):
            p = surf.Value(float(u), float(v)).Transformed(trsf)
            out[i] = (p.X(), p.Y(), p.Z())
        return out

    return pts, uvs, tris, evaluate


def face_triangulation(face: TopoDS_Face) -> tuple[np.ndarray, np.ndarray] | None:
    """``(verts, tris)`` of a meshed face in world coordinates, or ``None``.

    ``verts`` is ``(N, 3)`` float, ``tris`` is ``(M, 3)`` int with 0-based
    indices into ``verts``. The face's ``TopLoc_Location`` is applied here so
    callers always see world coordinates.
    """
    loc = TopLoc_Location()
    tri = BRep_Tool.Triangulation_s(face, loc)
    if tri is None:
        return None
    trsf = loc.Transformation()
    n = tri.NbNodes()
    verts = np.empty((n, 3), dtype=float)
    for i in range(1, n + 1):
        p = tri.Node(i).Transformed(trsf)
        verts[i - 1] = (p.X(), p.Y(), p.Z())
    m = tri.NbTriangles()
    tris = np.empty((m, 3), dtype=np.int64)
    for i in range(1, m + 1):
        t = tri.Triangle(i)
        tris[i - 1] = (t.Value(1) - 1, t.Value(2) - 1, t.Value(3) - 1)
    return verts, tris


# --- STEP I/O ---------------------------------------------------------------


def quiet_kernel() -> None:
    """Stop OCCT's own informational chatter from polluting the console.

    OCCT's STEP writer reports transfer statistics through the default
    messenger at Info level. The run's own log is the place for that, so the
    console printer is raised to warnings-and-worse; genuine problems still get
    through.
    """
    try:
        from OCP.Message import Message, Message_Gravity

        printers = Message.DefaultMessenger_s().Printers()
        for i in range(printers.Length()):
            printers.Value(i + 1).SetTraceLevel(Message_Gravity.Message_Warning)
    except Exception:
        pass  # cosmetic only


def _configure_step_units() -> None:
    """Pin STEP I/O to millimetres (specification.md §5).

    Set defensively even though mm is the default: a mismatched unit would
    silently corrupt every downstream dimension rather than fail.
    """
    Interface_Static.SetCVal_s("xstep.cascade.unit", "MM")
    Interface_Static.SetCVal_s("write.step.unit", "MM")


def read_step(path: str) -> TopoDS_Shape:
    """Read a STEP file into a single shape, raising :class:`InputGeometryError`."""
    if not os.path.isfile(path):
        raise InputGeometryError(f"Input STEP file not found: {path}")
    _configure_step_units()
    reader = STEPControl_Reader()
    status = reader.ReadFile(path)
    if status != IFSelect_ReturnStatus.IFSelect_RetDone:
        raise InputGeometryError(
            f"Could not read STEP file (OCCT status {status}): {path}"
        )
    reader.TransferRoots()
    if reader.NbShapes() < 1:
        raise InputGeometryError(f"STEP file contains no transferable shapes: {path}")
    shape = reader.OneShape()
    if shape.IsNull():
        raise InputGeometryError(f"STEP file produced a null shape: {path}")
    return shape


def write_step(shape: TopoDS_Shape, path: str, part_name: str) -> None:
    """Write ``shape`` as an AP214 STEP file in millimetres.

    ``part_name`` is applied to the file's header by
    :func:`latticegen2.stepout.rewrite_step_header`, which runs after this.
    """
    _configure_step_units()
    Interface_Static.SetCVal_s("write.step.schema", "AP214IS")
    Interface_Static.SetCVal_s("write.step.product.name", part_name)
    writer = STEPControl_Writer()
    status = writer.Transfer(shape, STEPControl_StepModelType.STEPControl_AsIs)
    if status != IFSelect_ReturnStatus.IFSelect_RetDone:
        raise OutputError(f"STEP transfer failed (OCCT status {status}): {path}")
    status = writer.Write(path)
    if status != IFSelect_ReturnStatus.IFSelect_RetDone:
        raise OutputError(f"STEP write failed (OCCT status {status}): {path}")
