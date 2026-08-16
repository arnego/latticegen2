"""Thin, typed helpers over the OCCT kernel (via the OCP bindings).

Everything in this module is a small wrapper whose only job is to keep OCCT's
C++-shaped API (static ``*_s`` methods, out-parameters, 1-based arrays) out of
the rest of the codebase. No lattice logic lives here.

OCP is used rather than a higher-level wrapper because this design needs the
whole OCCT surface: sewing, shared-topology construction, located instances and
exact validity checking.
"""

from __future__ import annotations

import os
from typing import Iterable, Iterator

import numpy as np

from OCP.BRep import BRep_Builder, BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Surface
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
from OCP.BRepTools import BRepTools, BRepTools_ReShape
from OCP.Bnd import Bnd_Box
from OCP.GProp import GProp_GProps
from OCP.GeomAbs import GeomAbs_SurfaceType
from OCP.IFSelect import IFSelect_ReturnStatus
from OCP.Interface import Interface_Static
from OCP.STEPControl import STEPControl_Reader, STEPControl_StepModelType, STEPControl_Writer
from OCP.ShapeFix import ShapeFix_Edge
from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain
from OCP.TopAbs import TopAbs_ShapeEnum
from OCP.TopExp import TopExp, TopExp_Explorer
from OCP.TopTools import TopTools_IndexedDataMapOfShapeListOfShape
from OCP.TopLoc import TopLoc_Location
from OCP.TopoDS import TopoDS, TopoDS_Compound, TopoDS_Face, TopoDS_Shape, TopoDS_Shell
from OCP.gp import gp_Pnt, gp_Trsf, gp_Vec

from .errors import InputGeometryError, OutputError, ProcessingError

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
    bounding box (docs/algorithm.md §5.1).
    """
    props = GProp_GProps()
    BRepGProp.SurfaceProperties_s(shape, props)
    return props.Mass()


def centroid(shape: TopoDS_Shape) -> np.ndarray:
    """Area-weighted centroid of a face (or any shape ``GProp`` can measure).

    Used to disambiguate which of several candidate nodes a re-tagged face
    belongs to after :func:`latticegen2.boundary.fuse_disagreeing_pairs` fuses
    two junctions: :func:`latticegen2.junction.is_cap_plane_face` tests only the
    one axis a cap's half-strut id names, so it can pass for more than one node
    in a group that share a coordinate along an orthogonal axis. Proximity to
    each candidate's own ideal cap centre resolves it.
    """
    props = GProp_GProps()
    BRepGProp.SurfaceProperties_s(shape, props)
    p = props.CentreOfMass()
    return np.array([p.X(), p.Y(), p.Z()])


def bounding_box(shape: TopoDS_Shape, tol: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """Axis-aligned ``(lo, hi)`` bounds.

    OCCT's box is a deliberate **over-estimate** (control-point hull for
    B-splines, untrimmed UV rectangle for planes), so it is only ever sound to
    use as an upper bound — never to test whether geometry *reaches* something
    (docs/algorithm.md §5.1).
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

    This is an *exact* check on the B-rep itself, not a mesh-based approximation
    of one, which is what makes it worth running on every output solid before a
    run reports success (docs/algorithm.md §9).
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


def faces_shell(faces: Iterable[TopoDS_Face]) -> TopoDS_Shell:
    """Bag ``faces`` into a shell, with no geometric work and no closure claim."""
    shell = TopoDS_Shell()
    builder = BRep_Builder()
    builder.MakeShell(shell)
    for f in faces:
        builder.Add(shell, f)
    return shell


def sew(shapes: Iterable[TopoDS_Shape], tolerance: float,
        cutting: bool = True) -> TopoDS_Shape:
    """Stitch shapes into shells along coincident free edges.

    Feeding already-built shells rather than loose faces is still worth doing —
    an edge already shared inside a shell is not a free edge and is never a
    candidate for merging.

    **But it does not make an already-sewn shell free to pass in, and the
    architecture was designed as though it did.** Measured (G5 in
    ``tools/prototypes/RESULTS.md``): adding one *closed* instanced shell — zero
    free edges, nothing for sewing to pair up — to a 4,000-piece sew multiplies
    its cost by more than an order of magnitude, purely on face count. The whole
    call also scales at about ``n^1.8`` in piece count, and OCCT's optional
    phases (``Cutting``, ``Analysis``, ``SameParameter``) account for under 2 %
    of it, so there is no configuration that rescues this.

    The consequence for callers: **never hand this the interior shell.** Its size
    scales with the volume of the part, which puts a superlinear term on the one
    path docs/algorithm.md §6 exists to keep linear. At ``cc=5, t=1`` on
    ``TD_HX_rehearsal_test`` that cost 4 h 45 m of a 5 h 04 m run. Sewing is now
    confined to the boundary layer (:mod:`latticegen2.weld`).

    ``cutting`` is OCCT's phase that splits free edges so they match. Where the
    two sides already match by construction it is pure waste, and G5a measures
    the whole optional-phase group at under 2 % either way, so switching it off
    costs nothing and removes a chance to re-partition an edge we intend to
    identify afterwards.
    """
    sewing = BRepBuilderAPI_Sewing(tolerance, True, True, cutting, False)
    for s in shapes:
        sewing.Add(s)
    sewing.Perform()
    return sewing.SewedShape()


def make_solid(shell: TopoDS_Shell) -> TopoDS_Shape:
    """Wrap a closed shell as a solid."""
    return BRepBuilderAPI_MakeSolid(shell).Solid()


def unify_same_domain(shape: TopoDS_Shape, unify_edges: bool = True) -> TopoDS_Shape:
    """Merge adjacent faces (and edges) that lie on the same underlying surface.

    A pure *representation* change: the point set is identical, only its
    partition into faces differs. Instancing produces a lot of this redundancy
    by construction — across every shared mid-strut interface the two junctions'
    lateral faces are coplanar and share an edge, but nothing ever merges them
    (docs/algorithm.md §9).

    ``unify_edges`` also concatenates the collinear edge pairs left inside a
    merged face's wire. It is the part of the algorithm that can throw on
    otherwise sound geometry, and it is worth far less than the face merging —
    see :func:`latticegen2.pipeline._unify`, which retries without it.

    OCCT's default tolerances are used deliberately. The faces this is expected
    to merge are exactly coplanar by construction, so the defaults suffice;
    loosening them would let genuinely distinct faces merge, which would change
    the geometry rather than just how it is described.
    """
    upgrade = ShapeUpgrade_UnifySameDomain(shape, unify_edges, True, False)
    upgrade.Build()
    return upgrade.Shape()


def fix_vertex_tolerances(faces: Iterable[TopoDS_Face]) -> tuple[int, int]:
    """Correct vertices recorded as sitting off their edge's curve.

    Returns ``(repaired, still_invalid)`` counted in faces.

    Sewing can leave an edge whose vertex lies off the edge's own 3D curve, with
    that vertex's tolerance inflated to *exactly* the distance — so the validity
    test sits on the knife edge and falls the wrong way, and
    ``BRepCheck_Analyzer`` rejects both faces sharing the edge. Measured on
    `TD_HX_rehearsal_test` at ``cc=5, t=1``: 17 such edges, 34 faces, deviations
    of 2.474044e-05 and 3.316370e-04 mm. See docs/algorithm.md §8 and
    ``tools/prototypes/RESULTS.md`` G11.

    ``ShapeFix_Edge.FixVertexTolerance`` is OCCT's own repair for it and adjusts
    the *recorded tolerance* in place — it moves no geometry, which is why this
    is safe to run on a shell that is already proven watertight. That is checked
    rather than assumed: a repaired face's surface area must come back
    bit-identical, and anything else is a hard failure.

    Only faces the analyzer already rejects are examined, so a sound boundary
    layer pays one validity check per face and nothing more.
    ``BRepLib.UpdateTolerances`` was measured as the obvious alternative and
    rejected: it repairs the planar and cylindrical cases but not the B-spline
    one (G11).
    """
    fixer = ShapeFix_Edge()
    repaired = still_invalid = 0
    for face in faces:
        if BRepCheck_Analyzer(face).IsValid():
            continue
        before = area(face)
        touched = False
        for e in _explore(face, TopAbs_ShapeEnum.TopAbs_EDGE):
            edge = TopoDS.Edge_s(e)
            if BRepCheck_Analyzer(edge).IsValid():
                continue
            fixer.FixVertexTolerance(edge, face)
            touched = True
        if not touched:
            still_invalid += 1
            continue
        after = area(face)
        if after != before:
            raise ProcessingError(
                f"Correcting a vertex tolerance changed a boundary face's area "
                f"from {before:.12g} to {after:.12g} mm^2. This repair adjusts "
                f"only recorded tolerances and must move no geometry at all."
            )
        if BRepCheck_Analyzer(face).IsValid():
            repaired += 1
        else:
            still_invalid += 1
    return repaired, still_invalid


def remove_pinhole_wires(shape: TopoDS_Shape, tol: float) -> tuple[TopoDS_Shape, int]:
    """Drop inner wires that bound no area. ``(shape, n_removed)``.

    A strut grazing the input surface almost tangentially leaves a face carrying
    an extra **inner wire** made of a single edge a few microns long whose two
    endpoints do not even meet — a pinhole. It encloses nothing, but it is a
    wire, so the edge is used once and :func:`latticegen2.weld.shell_defects`
    rejects the shell for it. Measured on `TD_HX_rehearsal_test` at ``cc=5, t=1``:
    two of them, 3.171690e-06 and 5.808982e-06 mm, on planar faces of 1.19 and
    1.25 mm², in a solid ``BRepCheck_Analyzer`` calls valid. See
    docs/algorithm.md §7 and ``tools/prototypes/RESULTS.md`` G10.

    **OCCT's own repairs do not touch these, which is why this is hand-rolled.**
    ``ShapeFix_Wireframe`` targets small *edges* between two faces and reports
    no candidates here at any precision from 1e-5 to 1e-2;
    ``ShapeFix_Face.FixSmallAreaWire`` targets small *wires* and returns false
    without removing anything. Both expect a well-formed closed wire, and a
    single non-closing edge is not one.

    A wire is removed only when **all** of these hold, which is what makes the
    repair unable to open a hole:

    * it is not the face's outer wire;
    * every one of its edges is shorter than ``tol``, so the region it could
      bound is smaller than ``tol²`` — nothing, against a strut of side ``t``;
    * every one of its edges is used exactly **once** in ``shape``, i.e. is
      already unpaired. An inner wire properly shared with a neighbouring face
      is left alone, however small.

    So this only ever deletes edges that are already defects, and only when they
    bound nothing. Measured on the piece above: surface area and cap areas
    unchanged **exactly**, volume drift 2.7e-15 (machine precision, the same
    order docs/algorithm.md §9 records for planar geometry), the solid still
    valid, and ``shell_defects`` going from ``(2, 0)`` to ``(0, 0)``.
    """
    edge_faces = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(
        shape, TopAbs_ShapeEnum.TopAbs_EDGE, TopAbs_ShapeEnum.TopAbs_FACE, edge_faces
    )

    def is_pinhole(wire) -> bool:
        n = 0
        for e in _explore(wire, TopAbs_ShapeEnum.TopAbs_EDGE):
            n += 1
            edge = TopoDS.Edge_s(e)
            if BRep_Tool.Degenerated_s(edge):
                continue
            props = GProp_GProps()
            BRepGProp.LinearProperties_s(edge, props)
            if props.Mass() >= tol:
                return False
            i = edge_faces.FindIndex(edge)
            if i == 0 or edge_faces.FindFromIndex(i).Extent() != 1:
                return False        # shared with another face: not ours to remove
        return n > 0

    reshape = BRepTools_ReShape()
    n_removed = 0
    for f in _explore(shape, TopAbs_ShapeEnum.TopAbs_FACE):
        face = TopoDS.Face_s(f)
        outer = BRepTools.OuterWire_s(face)
        for wire in _explore(face, TopAbs_ShapeEnum.TopAbs_WIRE):
            if wire.IsSame(outer):
                continue
            if is_pinhole(wire):
                reshape.Remove(wire)
                n_removed += 1
    if not n_removed:
        return shape, 0
    fixed = reshape.Apply(shape)
    if fixed is None or fixed.IsNull():
        return shape, 0
    return fixed, n_removed


def count_subshapes(shape: TopoDS_Shape) -> tuple[int, int]:
    """``(faces, edges)`` in ``shape`` — the compactness of its B-rep."""
    return (
        sum(1 for _ in _explore(shape, TopAbs_ShapeEnum.TopAbs_FACE)),
        sum(1 for _ in _explore(shape, TopAbs_ShapeEnum.TopAbs_EDGE)),
    )


def shells(shape: TopoDS_Shape) -> list[TopoDS_Shell]:
    """Every shell in ``shape``."""
    return [TopoDS.Shell_s(s) for s in _explore(shape, TopAbs_ShapeEnum.TopAbs_SHELL)]


# --- Meshing ----------------------------------------------------------------


def mesh_shape(shape: TopoDS_Shape, deflection: float, angle: float = 0.35) -> None:
    """Tessellate ``shape`` in place with a chordal-deviation bound.

    ``deflection`` is OCCT's linear (chordal) deviation target, so it *is* the
    ``d`` of docs/algorithm.md §5 rather than a proxy for it — a planar face
    meshes to a couple of exact triangles regardless of size, while curvature
    drives refinement automatically, so no separate maximum element size is
    needed.
    """
    BRepMesh_IncrementalMesh(shape, deflection, False, angle, True)


def is_planar(face: TopoDS_Face) -> bool:
    """Is this face a plane? Planes are tessellated exactly, with zero deviation."""
    return BRepAdaptor_Surface(face).GetType() == GeomAbs_SurfaceType.GeomAbs_Plane


def face_uv_triangulation(face: TopoDS_Face):
    """``(points, uv, tris, evaluate)`` for a meshed face, or ``None``.

    ``points`` are world-space triangulation nodes, ``uv`` their surface
    parameters, ``tris`` the 0-based index triples, and ``evaluate(uv_array)``
    maps parameters to world-space points on the *true* surface. Together these
    are what :func:`latticegen2.classify.measure_face_deviation` needs to
    compare the mesh against the geometry it approximates.
    """
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
