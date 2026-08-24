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
from typing import Iterable, Iterator, NamedTuple

import numpy as np

from OCP.Adaptor3d import Adaptor3d_CurveOnSurface
from OCP.BRep import BRep_Builder, BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepBuilderAPI import (
    BRepBuilderAPI_MakeFace,
    BRepBuilderAPI_MakePolygon,
    BRepBuilderAPI_MakeSolid,
    BRepBuilderAPI_Sewing,
)
from OCP.BRepCheck import BRepCheck_Analyzer, BRepCheck_Status, BRepCheck_Wire
from OCP.BRepGProp import BRepGProp
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.BRepPrimAPI import BRepPrimAPI_MakePrism
from OCP.BRepTools import BRepTools, BRepTools_ReShape
from OCP.Bnd import Bnd_Box
from OCP.GProp import GProp_GProps
from OCP.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCP.GeomAbs import GeomAbs_SurfaceType
from OCP.GeomAdaptor import GeomAdaptor_Curve, GeomAdaptor_Surface
from OCP.GeomLib import GeomLib_CheckCurveOnSurface
from OCP.IFSelect import IFSelect_ReturnStatus
from OCP.Interface import Interface_Static
from OCP.STEPControl import STEPControl_Reader, STEPControl_StepModelType, STEPControl_Writer
from OCP.ShapeFix import ShapeFix_Edge
from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain
from OCP.TopAbs import TopAbs_ShapeEnum
from OCP.TopExp import TopExp, TopExp_Explorer
from OCP.TopTools import (
    TopTools_IndexedDataMapOfShapeListOfShape,
    TopTools_IndexedMapOfShape,
)
from OCP.TopLoc import TopLoc_Location
from OCP.TopoDS import (
    TopoDS,
    TopoDS_Compound,
    TopoDS_Edge,
    TopoDS_Face,
    TopoDS_Shape,
    TopoDS_Shell,
    TopoDS_Vertex,
)
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

    **Run with OCCT's own parallel flag on** (``theIsParallel``, the third
    constructor argument). This is the one heavy call in the pipeline that has
    such a flag: docs/specification.md §11 records that
    ``ShapeUpgrade_UnifySameDomain`` has none — "no ``SetRunParallel``, no
    thread-pool hook, so there is nothing to switch on either" — and that is
    half of why sub-body parallelism is closed for `simplify` and open here.
    The threads are OCCT's own native ones, so G17's GIL result does not bind
    them, and the verdict stays OCCT's rather than a conjunction this project
    assembles itself. Gate G18 measured **1.60x** at 3.43 core-equivalents on a
    real trimmed lattice solid, with the verdict unchanged on valid solids *and*
    on all four of the real invalid faces committed from the `cc=5, t=1`
    rehearsal (``test/*.brep``, docs/algorithm.md §8) — a control that matters
    more than the speedup, since this is the exit-4 gate and a check that only
    ever agrees on good geometry would prove nothing (G10).

    Thread count is bounded by :func:`latticegen2.parallel.set_thread_budget`
    so ``--cores`` keeps meaning what specification.md §3 says it means. That
    bound is why this call is made on the master rather than in a worker
    (docs/algorithm.md §9): W worker processes each launching W OCCT threads
    would be W² threads on W cores.
    """
    return BRepCheck_Analyzer(shape, True, True).IsValid()


FACE_SCAN_CHUNK = 20_000
"""Faces per window in :func:`invalid_faces`.

Sized from memory, not from speed: G22 measured the analyzer holding ~14 kB per
face while it is alive, so scanning the rehearsal's 301,505 boundary faces in
one call would hold ~4.2 GB at a point in `stitch` already carrying ~3 GB. A
20,000-face window is ~0.3 GB and measures within 1 % of the unchunked scan —
the speedup is flat from 1,000 faces upward, so there is nothing to trade away
by keeping the window small.
"""


def invalid_faces(
    faces: Iterable[TopoDS_Face], chunk: int = FACE_SCAN_CHUNK
) -> list[TopoDS_Face]:
    """The faces ``BRepCheck_Analyzer`` rejects standalone, in input order.

    The predicate is ``BRepCheck_Analyzer(face).IsValid()``, one face at a
    time, unchanged — what changes is how the faces that need it are found. A
    boundary layer is hundreds of thousands of faces of which a handful are
    ever rejected (19 of 301,505 on the `cc=5, t=1` rehearsal), and asking one
    analyzer per face cost **44.6 s** of that run's `stitch`, serial on the
    master (docs/algorithm.md §8). So a batch scan proposes candidates and the
    standalone check disposes of them:

    * **Propose.** One ``BRepCheck_Analyzer`` per window of ``chunk`` faces,
      packed into a ``TopoDS_Compound``, with ``theIsParallel`` on — OCCT's own
      native threads, the flag G18 measured for the validity gate.
    * **Dispose.** Every candidate is then re-checked **alone**, which is the
      predicate everything downstream is calibrated against. Candidates are
      tens of faces, so this costs nothing measurable.

    **The confirmation stage is insurance rather than a measured necessity, and
    it is worth being precise about which.** G22 compared the two face for face
    over a 17,308-face corpus carrying the four real committed faults and found
    zero disagreements, so the batch scan has never been *observed* to differ.
    What that corpus cannot exercise is the case where it plausibly could: its
    faults are loose faces, while a sewn layer's faults have neighbours, and a
    compound analyzer holds one result per subshape shared between every face
    using it. Confirming tens of candidates costs nothing next to scanning
    hundreds of thousands, so the check that decides is the one this repair was
    calibrated against, and the batch scan only ever proposes.

    Whatever the two stages disagree about, the disagreement can only run one
    way: the batch scan records at least the statuses the standalone one does,
    so it can over-report and cannot hide a fault the standalone check would
    find.

    **This is not what produced the rehearsal's 15 phantom "still invalid"
    faces** — that was an ordering effect in the caller, and
    :func:`fix_vertex_tolerances` says where.
    """
    faces = list(faces)
    out: list[TopoDS_Face] = []
    for i in range(0, len(faces), chunk):
        window = faces[i:i + chunk]
        analyzer = BRepCheck_Analyzer(compound(window), True, True)
        out.extend(
            f for f in window
            if not analyzer.IsValid(f) and not BRepCheck_Analyzer(f).IsValid()
        )
    return out


# --- Export truth: is this shape describable outside OCCT's memory model? ---
#
# STEP AP214 carries exactly **one** modelling tolerance for a whole file — the
# `UNCERTAINTY_MEASURE_WITH_UNIT` of its `geometric_representation_context` —
# where an OCCT B-rep carries one per vertex, per edge and per face. Export
# therefore collapses N tolerances into 1 and import re-derives all N from that
# single number. Any face whose validity is *carried by* a locally fat tolerance
# is valid in this process and is not guaranteed valid in the file, and nothing
# in `BRepCheck_Analyzer` can see the difference: it asks whether the geometry
# agrees with itself to within the tolerances **recorded here**.
#
# Both measurements below exist to ask the question the analyzer cannot, and
# both are cheap, exact and local. See docs/algorithm.md §7.3 and §9.


TOLERANCE_FEATURE_RATIO_WARN = 1e-2
"""When a trim's recorded tolerance is worth reporting against its feature.

A **warning bar, not a gate**, and the distinction is measured rather than
cautious. On `SpiralTest.step` at ``cc=5, t=1`` this flags ~1 % of the 2,404
boundary pieces — a family of grazing trims against the part's swept B-spline
surface, all carrying the same 2.120e-02 mm the boolean fits there. Most of
those junctions are then welded into the 27,864 mm³ dominant body, where a fat
local description is absorbed by everything around it and the exported solid is
fine. Failing on this alone would refuse a part whose output is sound, which
docs/algorithm.md §11 forbids more strongly than it asks for any gate.

What it is *for* is naming junctions. Two of the six junctions forming that
part's 4.17 mm³ floating island rank **2nd and 4th of 2,404** here, at 2.907e-01
and 2.093e-01 — so by the end of the boundary stage, five minutes into the run
and before the junction graph exists, the geometry that cannot be exported has
already been pointed at. :func:`latticegen2.pipeline._check_component_tolerance`
combines that with the component the junction lands in, which is the sharp
predicate; :func:`latticegen2.pipeline._check_export_truth` is what actually
decides, on the finished solid, where the question is well posed.

Worst piece measured, for scale: 80 mm ball at ``cc=20 t=4`` **8.80e-06**,
`SpiralTest` at ``cc=5 t=1`` median **2.01e-05**, p90 **1.09e-03**, p99
**3.71e-02**, max **4.04e-01**."""


class ToleranceFeature(NamedTuple):
    """The worst "tolerance as a fraction of the feature" on one shape."""

    ratio: float
    """``max edge tolerance on a face / sqrt(that face's area)``, worst face."""
    tolerance: float
    """The recorded edge tolerance behind ``ratio``, in mm."""
    face_area: float
    """That face's exact trimmed area, in mm²."""
    where: tuple[float, float, float]
    """That face's centroid, so the junction can be found and inspected."""


_NO_TOLERANCE_FEATURE = ToleranceFeature(0.0, 0.0, 0.0, (0.0, 0.0, 0.0))


def tolerance_feature_ratio(faces: Iterable[TopoDS_Face]) -> ToleranceFeature:
    """How large the kernel's own recorded tolerance is against the feature.

    For every face, the largest tolerance recorded on any of its edges divided
    by ``sqrt`` of the face's exact trimmed area — a dimensionless reading of
    "the kernel needed *this* much slack to describe a feature *that* big".
    The worst face wins.

    **Why this quantity, and why at the source.** A boolean trimming a strut
    almost tangentially to a fat B-spline surface fits an intersection curve it
    can only place to within some distance, and records that distance as the
    edge's tolerance. While the tolerance stays far below the feature it bounds,
    the description is sound and the file is too. Once it approaches the feature
    size, the B-rep is only "valid" in the sense that everything is within a
    slack the size of the thing itself — and that is precisely the geometry that
    cannot survive §9's export, because the slack is the part STEP throws away.

    Measured, worst face per solid:

    ==========================================  ==========
    body                                        ratio
    ==========================================  ==========
    80 mm ball, ``cc=20 t=4``, 1,338 faces      2.11e-07
    test cylinder, ``cc=10 t=1.5``, 14,790 f    1.18e-03
    ...its second solid, 1,176 faces            8.01e-07
    ==========================================  ==========

    ``sqrt(area)`` rather than a bounding-box diagonal deliberately: a long thin
    sliver has a small area and a large diagonal, so the diagonal would flatter
    exactly the shape whose thin direction is the one the tolerance swallows.
    Erring high sends a junction down a reporting path; erring low ships it.

    Degenerate edges are skipped — they have no extent, so their tolerance
    bounds nothing — and a face of zero area contributes nothing rather than an
    infinity, since it has no feature to compare against and §8's own machinery
    is what deals with it.
    """
    worst = _NO_TOLERANCE_FEATURE
    for face in faces:
        tol = 0.0
        for e in _explore(face, TopAbs_ShapeEnum.TopAbs_EDGE):
            edge = TopoDS.Edge_s(e)
            if BRep_Tool.Degenerated_s(edge):
                continue
            tol = max(tol, BRep_Tool.Tolerance_s(edge))
        if tol <= 0.0:
            continue
        a = area(face)
        if a <= 0.0:
            continue
        ratio = tol / (a ** 0.5)
        if ratio > worst.ratio:
            c = centroid(face)
            worst = ToleranceFeature(
                ratio, tol, a, (float(c[0]), float(c[1]), float(c[2]))
            )
    return worst


class CurveOnSurface(NamedTuple):
    """How far a shape's pcurves stray from their own 3D curves."""

    worst: float
    """Largest max-deviation over all edge/face pairs, in mm. Reported only."""
    ratio: float
    """That deviation against the feature carrying it, worst pair. Reported
    only — :data:`LOOSE_AREA_FRACTION_MAX` records the measurement showing that a
    bar on this column refuses a sound body and passes a broken one."""
    loose_area_fraction: float
    """**The quantity the gate reads.** Share of the shape's surface area on
    faces that are loose at :data:`LOOSE_PCURVE_RATIO`."""
    loose_faces: int
    """How many faces that is."""
    over_tolerance: int
    """Pairs whose deviation exceeds the edge's own recorded tolerance.

    OCCT's own predicate — the one ``BOPAlgo_ArgumentAnalyzer``'s
    ``CurveOnSurfaceMode`` reports as ``InvalidCurveOnSurface`` — reproduced
    here so the two can be compared face for face (``test_export_truth.py``).
    Reported, never gated on: see ``ratio``."""
    pairs: int
    """Edge/face pairs carrying a pcurve, i.e. the pairs actually measured."""
    where: list[tuple[float, float, float]]
    """Sample midpoints of the offending edges, for the failure message."""
    where_worst: tuple[float, float, float]
    """Midpoint of the edge behind ``ratio``, so a failure names one place."""
    face_area: float
    """Area of the face behind ``ratio``, in mm²."""


#: Offending edge positions reported when the gate fires. Enough to localise the
#: region without turning one failure into a wall of coordinates — the same
#: reasoning behind §8's sampled free-edge positions.
CURVE_ON_SURFACE_SAMPLES = 10

LOOSE_PCURVE_RATIO = 1e-3
"""When one face counts as *loosely described*, against its own feature size.

Dimensionless: a face is loose when some edge's pcurve strays further than this
times ``sqrt`` of the face's own area from that edge's 3D curve. An absolute
figure in millimetres cannot express this — 2e-03 mm is a defect on a
0.05 mm² face and nothing at all on a 100 mm² one, and this tool's legal
parameters span both.

This is not a bar on its own. It is the predicate
:data:`LOOSE_AREA_FRACTION_MAX` then measures the *extent* of."""


LOOSE_AREA_FRACTION_MAX = 1e-2
"""How much of a body's surface may be loosely described. **Disproved as a gate.**

Kept as the level at which the reported figure is worth reading, and kept
documented because what it cost to learn is the useful part.

The idea was sound and the separation looked decisive. `SpiralTest.step` at
``cc=5, t=1`` produces two solids — the 27,864 mm3 lattice, which exports and
tessellates cleanly, and a 4.17 mm3 island whose tessellation carries 11
non-manifold edges. By worst face the **sound** body scores worse (3.07e-02
against 2.45e-02), and by fault count the **broken** body scores zero; only the
share of surface described more loosely than its own feature size ranked them
correctly, at 5.6178e-04 against 3.9685e-01, a 706x separation.

**Then it was measured against ground truth on a second part and it
false-positives.** Every one of the `TD_HX_rehearsal_test` rehearsal's fourteen
solids was round-tripped and tessellated:

======================  ======  ==========  ===========  ==========
body                    faces   loose area  faults ship  bad edges
======================  ======  ==========  ===========  ==========
`SpiralTest` island     36      3.97e-01    29           **11**
rehearsal unify 3       12      1.76e-01    0            0
rehearsal unify 5       12      1.76e-01    0            0
rehearsal unify 10      29      0           8            0
rehearsal unify 13      7       0           1            0
rehearsal, nine others  --      0           0            0
======================  ======  ==========  ===========  ==========

Exactly one body is broken. This quantity refuses two that are not, and a
**fault count after the round trip** — the other obvious candidate — refuses two
different ones that are not. Both would have refused, or in the shape of rule
this replaces deleted, geometry from a part that has been inspected and
accepted.

So the gate is :func:`exported_mesh_defects`, which asks the question directly,
and this stays as a diagnostic: it is cheap, it is exact, and it says *why* a
body is fragile. It decides nothing, and ``test_export_truth.py`` pins that so a
later change cannot quietly promote it back.

**The lesson is the one this project keeps re-learning** (specification.md §11):
a quantity that separates a known-good body from a known-bad one on the part
that motivated it is a ranking, not a calibration, and the second part is where
you find out."""


def curve_on_surface_deviations(
    shape: TopoDS_Shape, samples: int = CURVE_ON_SURFACE_SAMPLES
) -> CurveOnSurface:
    """Exact pcurve-versus-3D-curve deviation over every edge/face pair.

    The second half of the export-truth question, and the one that is about the
    *file* rather than about the kernel's slack. STEP stores both
    representations of an edge — a 3D curve and a pcurve per face — inside one
    ``surface_curve``, and OCCT's writer marks ``master_representation`` as
    ``.PCURVE_S1.``, so where the two disagree the file instructs the reader to
    believe the pcurve. The reader then re-derives every tolerance from the
    file's single uncertainty. An edge whose two representations disagree by
    more than the tolerance recorded *here* is therefore an edge whose geometry
    changes on import, by an amount nothing in the file bounds.

    ``GeomLib_CheckCurveOnSurface`` gives the true maximum by optimisation
    rather than by sampling, and it is the same predicate
    ``BOPAlgo_ArgumentAnalyzer``'s ``CurveOnSurfaceMode`` applies: measured on
    the 80 mm ball's output it reproduces that analyzer's count **exactly**
    (73 of 5,760 pairs) — but it also returns the *magnitude*, which the fault
    count does not, and which is the whole difference between a harmless
    9.08e-07 mm and a deviation the size of the face it sits on.

    **This needs no round trip, no boolean and no volume bar.** It is
    ``O(edge/face pairs)``, measured at 46–57 µs per pair, and it is affordable
    on every solid of every run — where a round-trip check is affordable only on
    small ones, and ``BOPAlgo_ArgumentAnalyzer``'s ``SelfInterMode`` (the half
    that makes that call superlinear: 0.67 s against ``CurveOnSurfaceMode``'s
    0.12 s on the same 1,176-face solid) is already covered by §9's validity
    gate.
    """
    worst = 0.0
    ratio = 0.0
    where_worst = (0.0, 0.0, 0.0)
    worst_area = 0.0
    over = 0
    pairs = 0
    total_area = 0.0
    loose_area = 0.0
    loose_faces = 0
    where: list[tuple[float, float, float]] = []
    for f in _explore(shape, TopAbs_ShapeEnum.TopAbs_FACE):
        face = TopoDS.Face_s(f)
        surface = BRep_Tool.Surface_s(face)
        if surface is None:
            continue
        face_area = area(face)
        if face_area <= 0.0:
            continue
        total_area += face_area
        root = face_area ** 0.5
        loose = False
        for e in _explore(face, TopAbs_ShapeEnum.TopAbs_EDGE):
            edge = TopoDS.Edge_s(e)
            if BRep_Tool.Degenerated_s(edge):
                continue
            deviation = _pcurve_deviation(edge, face, surface)
            if deviation is None:
                continue
            pairs += 1
            if deviation > worst:
                worst = deviation
            if deviation > BRep_Tool.Tolerance_s(edge):
                over += 1
                if len(where) < samples:
                    where.append(_edge_midpoint(edge))
            r = deviation / root
            if r > ratio:
                ratio = r
                where_worst = _edge_midpoint(edge)
                worst_area = face_area
            if r > LOOSE_PCURVE_RATIO and not loose:
                loose = True
                if len(where) < samples:
                    where.append(_edge_midpoint(edge))
        if loose:
            loose_faces += 1
            loose_area += face_area
    fraction = loose_area / total_area if total_area > 0.0 else 0.0
    return CurveOnSurface(worst, ratio, fraction, loose_faces, over, pairs,
                          where, where_worst, worst_area)



class SurfaceSeamGap(NamedTuple):
    """How far apart two faces' surfaces are along an edge they share."""

    worst: float
    """Largest separation found, in mm. Zero on a body whose faces meet."""
    where: tuple[float, float, float]
    """A point on the worse side of the widest gap, so it can be inspected."""
    edges: int
    """Shared edges measured, i.e. edges with exactly two faces and two pcurves."""
    over: int
    """How many of them exceed ``bar``."""


_NO_SEAM_GAP = SurfaceSeamGap(0.0, (0.0, 0.0, 0.0), 0, 0)

#: Samples per edge. The gap this exists to find is a modelling error rather
#: than a wobble, so it is essentially constant along the edge that carries it —
#: measured 4.2406e-02 mm at all 17 samples of `SpiralTest`'s worst seam, to six
#: figures. Seventeen is therefore generous rather than tuned.
SEAM_GAP_SAMPLES = 17


def surface_seam_gaps(
    shape: TopoDS_Shape, bar: float = 0.0, samples: int = SEAM_GAP_SAMPLES
) -> SurfaceSeamGap:
    """Do the two faces meeting at each edge actually *meet*?

    For every edge with exactly two faces, both faces' pcurves are evaluated at
    the same edge parameter and the two surface points compared. A B-rep is
    allowed to record a gap here and cover it with the edge's tolerance; that is
    legal, ordinary, and what almost every CAD file does at the micron scale.

    **It stops being harmless when a later operation cuts the region carrying it
    into features smaller than the gap**, which is exactly what this tool does
    to its input. Measured on the committed parts, with the run each is normally
    given:

    =============================  ============  ==============  ==========
    input                          worst gap     max edge tol    ships?
    =============================  ============  ==============  ==========
    `80mm-test-ball.step`          7.7e-14       1.0e-05         yes
    `test-cylinder.STEP`           4.5e-04       2.5e-04         yes
    `TD_HX_rehearsal_test.step`    8.1e-03       8.1e-03         yes
    `SpiralTest.step`              **4.2e-02**   2.1e-02         **no**
    =============================  ============  ==============  ==========

    The last row is a real defect in that file: two swept B-spline patches that
    run **42.4 microns apart along the whole of their shared edge** — constant to
    six figures at every sample, and confirmed independently by projecting each
    face's boundary onto the other's surface, which returns the same 4.2406e-02
    both ways. The two surfaces never meet; a 2.12e-02 mm edge tolerance is what
    holds the model together there. At ``t = 1`` mm the lattice cuts that seam
    into faces of 0.01 to 0.23 mm2, where the same gap is a hole, and the body
    built from them is the one docs/specification.md S11's export-truth gate
    refuses.

    **This is why re-fitting a pcurve cannot repair that body**, and the
    measurement is here so that conclusion is reproducible rather than
    remembered. No curve on either surface can close a gap between the surfaces
    themselves.

    ``bar`` counts how many edges exceed it and is otherwise not used: nothing
    here refuses anything. The separation between the part that fails and the
    part that ships is a factor of five on two data points, which is a ranking
    and not a calibration — docs/specification.md S11 records what happens when
    those are confused.

    Evaluating both pcurves at the **edge**'s own parameter is what makes this
    orientation-safe: a pcurve is stored same-parameter with the edge, so the
    two are compared at corresponding points whether or not the faces traverse
    the edge the same way.
    """
    ancestors = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(
        shape, TopAbs_ShapeEnum.TopAbs_EDGE, TopAbs_ShapeEnum.TopAbs_FACE, ancestors
    )
    worst = 0.0
    where = (0.0, 0.0, 0.0)
    measured = 0
    over = 0
    for i in range(1, ancestors.Extent() + 1):
        edge = TopoDS.Edge_s(ancestors.FindKey(i))
        if BRep_Tool.Degenerated_s(edge):
            continue
        faces = list(ancestors.FindFromIndex(i))
        if len(faces) != 2:
            continue
        fa, fb = TopoDS.Face_s(faces[0]), TopoDS.Face_s(faces[1])
        sa, sb = BRep_Tool.Surface_s(fa), BRep_Tool.Surface_s(fb)
        if sa is None or sb is None:
            continue
        ca = BRep_Tool.CurveOnSurface_s(edge, fa, 0.0, 0.0)
        cb = BRep_Tool.CurveOnSurface_s(edge, fb, 0.0, 0.0)
        if ca is None or cb is None:
            continue
        first, last = BRep_Tool.Range_s(edge)
        measured += 1
        widest = 0.0
        widest_at = None
        for k in range(samples):
            t = first + (last - first) * k / (samples - 1)
            ua, ub = ca.Value(t), cb.Value(t)
            pa = sa.Value(ua.X(), ua.Y())
            gap = pa.Distance(sb.Value(ub.X(), ub.Y()))
            if gap > widest:
                widest = gap
                widest_at = pa
        if bar > 0.0 and widest > bar:
            over += 1
        if widest > worst:
            worst = widest
            where = (round(widest_at.X(), 4), round(widest_at.Y(), 4),
                     round(widest_at.Z(), 4))
    if measured == 0:
        return _NO_SEAM_GAP
    return SurfaceSeamGap(worst, where, measured, over)

#: No size bound, deliberately. A solid this cannot afford to measure used to be
#: reported *unmeasured*, which put the dominant body of every production part —
#: the one where an unwritable description matters most — outside the only
#: detector with a clean record. Per the user's decision the cost is accepted
#: instead: docs/algorithm.md §9 removed a whole-output re-import for costing
#: 22 minutes, and this is that cost coming back, knowingly, for a correctness
#: check rather than for a count of solids.

DEFAULT_MESH_DEFLECTION = 0.05
"""Chordal deviation used when tessellating an exported body to check it.

Fine enough that a face of ~0.05 mm2 — the scale of the geometry this catches —
gets more than one triangle, which is what makes a non-manifold edge visible
at all."""


def exported_mesh_defects(
    shape: TopoDS_Shape, path: str, deflection: float = DEFAULT_MESH_DEFLECTION
) -> tuple[int, int]:
    """Tessellate ``shape`` **as the file writes it** and count broken edges.

    Returns ``(triangles, edges not used by exactly two triangles)``. The shape
    is written to STEP and read back first, because the question is about the
    artefact and not about the solid in memory: STEP AP214 declares one
    tolerance for a whole file where the B-rep carries one per subshape, so the
    body a reader reconstructs is not the body written.

    **This is the one instrument with a clean record on both parts, and three
    cheaper proxies were disproved against it.** Ground truth, measured on the
    14 solids of the `cc=5, t=1` rehearsal and the 2 of `SpiralTest.step`:

    ===========================  ======  ==========  ===========  ==========
    body                         faces   loose area  faults ship  bad edges
    ===========================  ======  ==========  ===========  ==========
    `SpiralTest` island          36      3.97e-01    29           **11**
    rehearsal unify 3            12      1.76e-01    0            0
    rehearsal unify 5            12      1.76e-01    0            0
    rehearsal unify 10           29      0           8            0
    rehearsal unify 13           7       0           1            0
    rehearsal, nine others       --      0           0            0
    ===========================  ======  ==========  ===========  ==========

    Exactly one of those bodies is broken, and only the last column says so.
    A **fault count after the round trip** fires on unify 10 and 13, which
    tessellate perfectly. The **share of surface loosely described**
    (:data:`LOOSE_AREA_FRACTION_MAX`) fires on unify 3 and 5, which also
    tessellate perfectly. Both would refuse — or, in the shape of the rule this
    replaces, delete — sound geometry from a part that has been accepted.

    Non-manifold tessellation is also the symptom that found the island in the
    first place: 11 bad edges in 147 triangles, against 105,628 triangles of
    lattice beside it carrying none. It is the property that makes a body
    useless downstream, rather than a proxy for it.

    The mesh is compared on **rounded coordinates** rather than on topology,
    deliberately: what breaks in a reader is that two triangles which should
    share an edge do not, and that is a question about points, not about the
    B-rep's opinion of them.

    **Run on every solid, with no size bound.** That is the expensive decision
    on this branch and it was taken deliberately: bounding it by face count is
    what puts a production part's dominant body outside the check, and a body
    nobody looks at is exactly where an unwritable description would do the most
    damage. Points are interned to integers and edges counted by integer key
    rather than by tuples of floats, which is the one concession made to the
    size of shape this now has to survive; the cost at rehearsal scale is
    **not yet measured** and is expected to be tens of minutes.
    """
    write_step(shape, path, "export-truth-probe")
    shipped = read_step(path)
    BRepMesh_IncrementalMesh(shipped, deflection, False, 0.2, True)
    # A point is interned to a small integer and an edge to one integer built
    # from its two endpoints, rather than keeping tuples of floats as keys. On a
    # lattice this is the difference between a dict of a few hundred megabytes
    # and one of several gigabytes, and the check now runs on solids of ~600,000
    # faces where it used to be skipped.
    point_id: dict[tuple[int, int, int], int] = {}
    counts: dict[int, int] = {}
    triangles = 0
    quantum = 1e7                       # 1e-7 mm, OCCT's own confusion
    for f in _explore(shipped, TopAbs_ShapeEnum.TopAbs_FACE):
        face = TopoDS.Face_s(f)
        loc = TopLoc_Location()
        tri = BRep_Tool.Triangulation_s(face, loc)
        if tri is None:
            continue
        trsf = loc.Transformation()
        ids = []
        for i in range(1, tri.NbNodes() + 1):
            p = tri.Node(i).Transformed(trsf)
            key = (round(p.X() * quantum), round(p.Y() * quantum),
                   round(p.Z() * quantum))
            got = point_id.get(key)
            if got is None:
                got = point_id[key] = len(point_id)
            ids.append(got)
        for i in range(1, tri.NbTriangles() + 1):
            a, b, c = tri.Triangle(i).Get()
            triangles += 1
            for u, v in ((a, b), (b, c), (c, a)):
                lo, hi = ids[u - 1], ids[v - 1]
                if lo > hi:
                    lo, hi = hi, lo
                edge = lo * 4_294_967_296 + hi
                counts[edge] = counts.get(edge, 0) + 1
    return triangles, sum(1 for n in counts.values() if n != 2)


def _pcurve_deviation(edge: TopoDS_Edge, face: TopoDS_Face, surface) -> float | None:
    """Max distance between one edge's 3D curve and its pcurve on ``face``.

    ``None`` where the pair carries no pcurve or the optimiser does not
    converge — both of which mean "nothing measured", never "measured as fine",
    so the caller counts neither as a pass.
    """
    curve3d = BRep_Tool.Curve_s(edge, 0.0, 0.0)
    if curve3d is None:
        return None
    curve2d = BRep_Tool.CurveOnSurface_s(edge, face, 0.0, 0.0)
    if curve2d is None:
        return None
    first3d, last3d = BRep_Tool.Range_s(edge)
    first2d, last2d = BRep_Tool.Range_s(edge, face)
    check = GeomLib_CheckCurveOnSurface(GeomAdaptor_Curve(curve3d, first3d, last3d))
    check.Perform(
        Adaptor3d_CurveOnSurface(
            Geom2dAdaptor_Curve(curve2d, first2d, last2d), GeomAdaptor_Surface(surface)
        )
    )
    return check.MaxDistance() if check.IsDone() else None


def _edge_midpoint(edge: TopoDS_Edge) -> tuple[float, float, float]:
    """A point on the edge, for a failure message to point at."""
    curve = BRep_Tool.Curve_s(edge, 0.0, 0.0)
    first, last = BRep_Tool.Range_s(edge)
    p = curve.Value(0.5 * (first + last))
    return (round(p.X(), 4), round(p.Y(), 4), round(p.Z(), 4))


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


def unify_same_domain(
    shape: TopoDS_Shape, unify_edges: bool = True, unify_faces: bool = True
) -> TopoDS_Shape:
    """Merge adjacent faces (and edges) that lie on the same underlying surface.

    A pure *representation* change: the point set is identical, only its
    partition into faces differs. Instancing produces a lot of this redundancy
    by construction — across every shared mid-strut interface the two junctions'
    lateral faces are coplanar and share an edge, but nothing ever merges them
    (docs/algorithm.md §9).

    The two passes are separable, and :func:`latticegen2.pipeline._unify_one`
    runs them as two calls. ``unify_faces`` merges the coplanar faces — the part
    that shrinks the output. ``unify_edges`` then concatenates the collinear edge
    pairs left inside a merged face's wire. Splitting them costs nothing and
    describes the solid identically; it exists so the face merge can run with
    edge merging off, which is what a tiled unification needs (see
    ``_unify_one``, and G13 in `tools/prototypes/RESULTS.md`).

    OCCT's default tolerances are used deliberately. The faces this is expected
    to merge are exactly coplanar by construction, so the defaults suffice;
    loosening them would let genuinely distinct faces merge, which would change
    the geometry rather than just how it is described.
    """
    upgrade = ShapeUpgrade_UnifySameDomain(shape, unify_edges, unify_faces, False)
    upgrade.Build()
    return upgrade.Shape()


# A vertex tolerance this repair widened may grow by at most this factor over
# the one the kernel itself recorded. The repair is allowed to nudge a figure
# OCCT chose, never to invent a new scale for it. Measured need on the four
# `cc=5, t=1` rehearsal faces: 1.05x, 1.05x, 1.25x, 1.10x (G12).
SELF_INTERSECT_TOL_GROWTH = 4.0

# ...and never above this in absolute terms, whatever the kernel recorded. The
# CLI's smallest legal strut is ``t = 0.4`` mm (specification.md §3), so this
# sits a hundredfold below the smallest feature any legal run can produce.
SELF_INTERSECT_MAX_VERTEX_TOL = 4e-3

# ...except that a *fixed* absolute cap is only right for the part it was
# calibrated on, and this one was calibrated on `TD_HX_rehearsal_test`, whose
# vertices carry tolerances of 8.7e-04 to 1.5e-03 mm (G12). A vertex whose
# recorded tolerance already **exceeds** the cap cannot be widened at all: the
# first candidate step is refused, nothing grows, and the repair reports the
# face as residual without having tried. That is not a bound on growth, it is
# an off switch, and it fired the first time a part with fatter trims arrived.
# On `SpiralTest.step` at ``cc=5, t=1`` the falsely self-intersecting wire's
# shared vertex carries **6.573e-02 mm** — sixteen times the cap — and needs
# exactly one 1.25x step to clear, which the relative bound above allows
# comfortably and this one forbade outright.
#
# So the absolute cap scales with the run's own smallest real feature instead
# of with one part's measurements. A tenth of ``t`` reproduces the same order
# as the fixed figure at small ``t`` while admitting what the kernel itself
# records at larger ``t``: 0.04 mm at ``t = 0.4``, 0.1 mm at ``t = 1`` (the
# 8.216e-02 mm the spiral needs, with 22 % to spare). The relative bound stays
# the one that actually governs — every case measured needs 1.05x to 1.25x,
# far inside 4x — and this remains the backstop against a kernel-recorded
# tolerance that is itself absurd.
SELF_INTERSECT_MAX_VERTEX_TOL_FRACTION = 0.1

# Geometric step used to search for the smallest widening that satisfies OCCT's
# own predicate. A search rather than a formula because the predicate is
# OCCT's, not ours: `BRepCheck_Wire::SelfIntersect` is asked directly whether it
# is satisfied, instead of this code re-deriving the rule it applies.
SELF_INTERSECT_TOL_STEP = 1.25


def _shared_vertices(a: TopoDS_Shape, b: TopoDS_Shape) -> list[TopoDS_Vertex]:
    """Vertices belonging to both edges, by topological identity."""
    first = {v.TShape(): v for v in _explore(a, TopAbs_ShapeEnum.TopAbs_VERTEX)}
    return [
        TopoDS.Vertex_s(first[v.TShape()])
        for v in _explore(b, TopAbs_ShapeEnum.TopAbs_VERTEX)
        if v.TShape() in first
    ]


def _self_intersecting_pair(
    face: TopoDS_Face,
) -> tuple[TopoDS_Edge, TopoDS_Edge] | None:
    """The two edges ``BRepCheck_Wire`` reports as self-intersecting, if any."""
    for w in _explore(face, TopAbs_ShapeEnum.TopAbs_WIRE):
        a, b = TopoDS_Edge(), TopoDS_Edge()
        status = BRepCheck_Wire(TopoDS.Wire_s(w)).SelfIntersect(face, a, b, False)
        if status != BRepCheck_Status.BRepCheck_NoError:
            return a, b
    return None


def _widen_self_intersection_vertices(
    face: TopoDS_Face, max_vertex_tol: float | None = None
) -> bool:
    """Widen the shared vertex of a falsely self-intersecting adjacent pair.

    Returns whether anything was widened. Moves no geometry: it only raises a
    recorded tolerance, and ``BRep_Builder.UpdateVertex`` never lowers one.

    ``max_vertex_tol`` is the absolute ceiling, which the caller derives from
    the run's own ``t`` (:data:`SELF_INTERSECT_MAX_VERTEX_TOL_FRACTION`) rather
    than taking a fixed figure calibrated on one part — see that constant for
    the case where a fixed one silently disabled this rung.
    """
    builder = BRep_Builder()
    as_built: dict = {}
    touched = False
    while True:
        pair = _self_intersecting_pair(face)
        if pair is None:
            return touched
        grew = False
        for vertex in _shared_vertices(*pair):
            key = vertex.TShape()
            start = as_built.setdefault(key, BRep_Tool.Tolerance_s(vertex))
            want = BRep_Tool.Tolerance_s(vertex) * SELF_INTERSECT_TOL_STEP
            # Resolved here rather than bound as a default argument, so that
            # patching the module constant still reaches this — which is how
            # `test_export_truth.py` exercises the real repair with a raised
            # cap rather than a stand-in for it.
            ceiling = (
                SELF_INTERSECT_MAX_VERTEX_TOL if max_vertex_tol is None
                else max_vertex_tol
            )
            if want > min(start * SELF_INTERSECT_TOL_GROWTH, ceiling):
                continue  # capped: leave it for the caller to count as residual
            builder.UpdateVertex(vertex, want)
            grew = touched = True
        if not grew:
            return touched


def fix_vertex_tolerances(
    faces: Iterable[TopoDS_Face],
    max_vertex_tol: float | None = None,
) -> tuple[int, int]:
    """Correct vertex tolerances the boundary sew leaves wrong.

    Returns ``(repaired, still_invalid)`` counted in faces.

    ``max_vertex_tol`` bounds rung 2's widening absolutely. Callers inside the
    pipeline pass :data:`SELF_INTERSECT_MAX_VERTEX_TOL_FRACTION` times the
    run's ``t``. ``None`` falls back to :data:`SELF_INTERSECT_MAX_VERTEX_TOL`,
    read at call time so the constant can still be patched — right only for
    parts whose recorded tolerances are of the rehearsal's order, and kept so a
    caller with no lattice in hand still gets a bound.

    Two faults, both of them a *recorded tolerance* being wrong rather than any
    geometry being wrong, and both repaired the same way: by correcting the
    number, never by moving a point or rebuilding a face. That is what makes
    this safe to run on a shell ``assemble`` has already proven watertight — no
    ``TopoDS_Edge`` or ``TopoDS_Vertex`` object is ever replaced, so the
    every-edge-twice proof cannot be disturbed. It is checked rather than
    assumed: a repaired face's surface area must come back **bit-identical**,
    and anything else is a hard failure. That bar is exact rather than
    approximate because a tolerance is metadata; there is no quadrature noise
    to allow for.

    **Rung 1 — a vertex sitting off its edge's 3D curve.** Sewing can leave an
    edge whose vertex lies off the edge's own 3D curve, with that vertex's
    tolerance inflated to *exactly* the distance — so the validity
    test sits on the knife edge and falls the wrong way, and
    ``BRepCheck_Analyzer`` rejects both faces sharing the edge. Measured on
    `TD_HX_rehearsal_test` at ``cc=5, t=1``: 17 such edges, 34 faces, deviations
    of 2.474044e-05 and 3.316370e-04 mm. See docs/algorithm.md §8 and
    ``tools/prototypes/RESULTS.md`` G11.

    ``ShapeFix_Edge.FixVertexTolerance`` is OCCT's own repair for it and adjusts
    the recorded tolerance in place. ``BRepLib.UpdateTolerances`` was measured
    as the obvious alternative and rejected: it repairs the planar and
    cylindrical cases but not the B-spline one (G11).

    **Rung 2 — a false self-intersection at a shared vertex.** What rung 1
    leaves behind is a face every one of whose edges and vertices is valid
    *standalone*, and which passes ``IntersectWires``, ``ClassifyWires`` and
    ``OrientationOfWires`` — yet the analyzer rejects it. The fault is
    ``BRepCheck_SelfIntersectingWire``, reported for two edges that are
    **adjacent in the wire**: a short, tight-tolerance trim edge against the
    fat-tolerance B-spline the boolean fitted to the strut/input-surface
    intersection. Measured on the same rehearsal: 4 faces, one pair each, fat
    tolerances 8.741e-04 to 1.540e-03 mm (G12).

    It is not a real self-intersection. The two pcurves cross at exactly one
    point, and that point lies **at the shared vertex, inside its tolerance**
    (3.5e-04 mm against 8.7e-04 mm on the worst of the four) — the edges meet
    where they are supposed to. What OCCT keys on is the shared vertex's
    recorded tolerance, which the sew left a little too tight to swallow the
    crossing: widening it clears the check on all four, while widening the fat
    edge instead does nothing at any factor up to 5x (G12).

    So this rung widens that vertex, and asks OCCT's own predicate whether the
    result is enough rather than re-deriving the rule OCCT applies. The search
    is bounded twice — at :data:`SELF_INTERSECT_TOL_GROWTH` times the tolerance
    the kernel itself recorded, and at ``max_vertex_tol`` absolutely — so its
    failure mode is a face left in ``still_invalid`` for ``validate`` to
    report, never an unbounded tolerance. Widening is also monotonically *permissive*: every
    check that reads a vertex tolerance (the contextual vertex check, this
    self-intersection check) is a "within tolerance" test, so a neighbouring
    face sharing the vertex can only become more valid, never less.

    ``ShapeFix_Wire.FixSelfIntersection`` is the obvious tool and was measured:
    it is a **no-op** on all four, at either of its two modes.
    ``ShapeFix_Shape.Perform`` does fix all four, and is rejected — it moves
    geometry (up to 6.4e-04 relative area) and rebuilds faces, replacing 2 to 3
    edge objects per face. Minting new edges on an already-proven shell is the
    exact mechanism behind the seam-split regression that produced 118,760 open
    edges (G9), so it is not an acceptable price for a repair whose whole
    justification is that nothing needs to move. ``BRepLib.SameParameter``
    fixes three of the four at ~1e-09 drift, by re-fitting pcurves until the
    pair happens to separate; it is not used either, since it treats the
    symptom on a subset while the widening treats the mechanism on all four at
    zero drift.

    Only faces the analyzer already rejects are examined, so a sound boundary
    layer pays nothing beyond finding that out — and finding it out is the
    stage's own cost, 44.6 s of the `cc=5, t=1` rehearsal's `stitch` over
    301,505 faces (docs/specification.md §10). :func:`invalid_faces` finds them
    with a parallel batch scan per window, so the set below is unchanged and
    only the search for it got cheaper: measured as a controlled pair on that
    part, **44.1 s -> 22.6 s**, with the same 19 faces repaired and the output
    byte-identical.
    """
    fixer = ShapeFix_Edge()
    repaired = still_invalid = 0
    for face in invalid_faces(faces):
        # Re-checked here, not only when the candidate list was built, because
        # the two are not the same question. Rungs 1 and 2 widen tolerances on
        # vertices and edges that neighbouring faces share, and widening is
        # monotonically permissive — so repairing one face can make the next
        # one valid before this loop reaches it. Asking once up front counted
        # those neighbours as unrepaired: the rehearsal reported 19 corrected
        # and 15 "still invalid" where the per-face loop had always reported
        # none, on a run whose validity gate then passed all 14 solids.
        # Nothing here can move the other way: no repair invalidates a face.
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
        if not BRepCheck_Analyzer(face).IsValid():
            touched |= _widen_self_intersection_vertices(face, max_vertex_tol)
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


def _wire_closes(wire, tol: float = 1e-7) -> bool:
    """Whether ``wire``'s edges form a closed loop.

    Asked of the vertices rather than of OCCT's ``Closed`` flag, which
    ``BRep_Builder`` leaves unset on wires a boolean produced — the same reason
    :mod:`latticegen2.weld` computes closure itself rather than reading the
    flag. A closed loop uses every vertex an even number of times; a wire with
    loose ends does not.
    """
    counts: dict = {}
    for e in _explore(wire, TopAbs_ShapeEnum.TopAbs_EDGE):
        edge = TopoDS.Edge_s(e)
        if BRep_Tool.Degenerated_s(edge):
            continue
        for v in _explore(edge, TopAbs_ShapeEnum.TopAbs_VERTEX):
            key = TopoDS.Vertex_s(v).TShape()
            counts[key] = counts.get(key, 0) + 1
    return bool(counts) and all(n % 2 == 0 for n in counts.values())


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
    * every one of its edges is used exactly **once** in ``shape``, i.e. is
      already unpaired. An inner wire properly shared with a neighbouring face
      is left alone, however small. **This is the condition that carries the
      argument**: the repair can only ever delete edges that are already
      defects;
    * and it bounds nothing — either because every edge is shorter than ``tol``,
      so the region it could enclose is under ``tol²``, **or** because the wire
      **does not close**, which is a stronger statement and needs no threshold
      at all. An open wire encloses nothing at any length.

    So this only ever deletes edges that are already defects, and only when they
    bound nothing. Measured on the piece above: surface area and cap areas
    unchanged **exactly**, the solid still valid, and ``shell_defects`` going
    from ``(2, 0)`` to ``(0, 0)``.

    **The "does not close" arm exists because the length bar is not what makes
    this safe, and a part eventually proved it.** `test-cylinder.STEP` at
    ``cc=4, t=0.8`` leaves junction (48,16,-27) carrying exactly this defect at
    **0.3267 mm** — a single-edge inner wire whose endpoints are that far
    apart, four orders above ``tol``, on a face of 1.29 mm². The trim's own
    output is an open shell for it, and `stitch` reported it as one free edge
    with nothing to weld to. Removing it leaves the face count unchanged, the
    surface area **bit-identical**, and the shell closed — the same evidence the
    micron-scale cases give, because the mechanism is the same and only the
    scale differs. Raising ``tol`` to cover it instead would have made a
    threshold load-bearing that the paragraph above is careful to say is not.

    **What the caller must not check is volume**, and the reason is a property
    of OCCT rather than of this repair: ``BRepGProp::VolumeProperties`` documents
    that its shape "must be exempt of any free boundary", and a pinhole wire *is*
    a free boundary — that is the definition of the defect. So the volume read
    off the piece before this runs carries a spurious term, which vanishes when
    the wire goes. It is length-independent and per-face, so it cannot be bounded
    by anything the repair controls (G19). :func:`only_inner_wires_dropped` is
    the check that replaces it, and it is exact.
    """
    edge_faces = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(
        shape, TopAbs_ShapeEnum.TopAbs_EDGE, TopAbs_ShapeEnum.TopAbs_FACE, edge_faces
    )

    def is_pinhole(wire) -> bool:
        n = 0
        tiny = True
        for e in _explore(wire, TopAbs_ShapeEnum.TopAbs_EDGE):
            n += 1
            edge = TopoDS.Edge_s(e)
            if BRep_Tool.Degenerated_s(edge):
                continue
            props = GProp_GProps()
            BRepGProp.LinearProperties_s(edge, props)
            if props.Mass() >= tol:
                tiny = False
            i = edge_faces.FindIndex(edge)
            if i == 0 or edge_faces.FindFromIndex(i).Extent() != 1:
                return False        # shared with another face: not ours to remove
        # Short enough to bound nothing, *or* structurally incapable of bounding
        # anything because it does not close. See the docstring.
        return n > 0 and (tiny or not _wire_closes(wire))

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


def only_inner_wires_dropped(
    before: TopoDS_Shape, after: TopoDS_Shape, n_removed: int
) -> str | None:
    """``None`` when ``after`` is ``before`` minus exactly ``n_removed`` inner wires.

    Returns a human-readable reason otherwise. This is the guard on
    :func:`remove_pinhole_wires`, and it is deliberately **structural and exact**
    rather than a comparison of two measured quantities: it asks whether the
    region the piece encloses is the same object-for-object, which is a question
    with a bit-exact answer, where volume is a question OCCT can only answer with
    a bias while the defect is still present (see that function, and G19).

    Four things are required of every face, and together they pin the enclosed
    region completely:

    * the face count is unchanged — nothing was merged away or split;
    * each face's **outer wire is the same object** as some face's in ``before``,
      one for one, so no face's outer boundary was rebuilt (which is the
      mechanism behind the seam-split regression, docs/algorithm.md §8);
    * the paired faces have the same orientation, so no face was flipped —
      a shell can stay closed while enclosing the wrong volume, and
      docs/algorithm.md §8's every-edge-twice proof exists because that has
      happened here before;
    * every wire ``after`` keeps was already a wire of its counterpart, and the
      total lost across the piece is exactly ``n_removed``. Nothing was gained,
      and nothing was lost that this repair did not account for.

    In practice ``BRepTools_ReShape`` rebuilds only the faces it touches — 27 of
    28 faces of the junction G19 was measured on come back ``IsSame`` — so this
    walk is over objects, not geometry, and costs one area evaluation per face.
    """
    fb, fa = faces(before), faces(after)
    if len(fb) != len(fa):
        return f"the face count changed, {len(fb)} -> {len(fa)}"

    outers = TopTools_IndexedMapOfShape()
    for f in fb:
        outers.Add(BRepTools.OuterWire_s(f))
    if outers.Extent() != len(fb):
        return "two faces share one outer wire, so they cannot be paired up"

    lost = 0
    matched: set[int] = set()
    for f in fa:
        i = outers.FindIndex(BRepTools.OuterWire_s(f))
        if i == 0:
            return "a face's outer wire was replaced, so its boundary was rebuilt"
        if i in matched:
            return "two faces came back with the same outer wire"
        matched.add(i)
        src = fb[i - 1]
        if f.Orientation() != src.Orientation():
            return "a face changed orientation"
        a_src, a_out = area(src), area(f)
        if a_out != a_src:
            return (
                f"a face's area changed, {a_src:.17g} -> {a_out:.17g} mm^2; a wire "
                f"that bounds no area cannot change it, so one that bounded "
                f"something was removed"
            )
        kept = TopTools_IndexedMapOfShape()
        for w in _explore(src, TopAbs_ShapeEnum.TopAbs_WIRE):
            kept.Add(w)
        n_after = 0
        for w in _explore(f, TopAbs_ShapeEnum.TopAbs_WIRE):
            n_after += 1
            if kept.FindIndex(w) == 0:
                return "a face gained a wire it did not have before"
        lost += kept.Extent() - n_after

    if lost != n_removed:
        return f"{lost} wire(s) went missing where {n_removed} were removed"
    return None


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


#: OCCT's ``write.precision.mode`` value for *Greatest* -- declare the file's
#: single uncertainty as the largest tolerance actually recorded anywhere in the
#: shape, rather than as the average of them all.
WRITE_PRECISION_GREATEST = 1


def write_step(shape: TopoDS_Shape, path: str, part_name: str) -> None:
    """Write ``shape`` as an AP214 STEP file in millimetres.

    ``part_name`` is applied to the file's header by
    :func:`latticegen2.stepout.rewrite_step_header`, which runs after this.

    **The file's one declared tolerance is the greatest the shape carries, not
    the average.** AP214 has room for exactly one
    ``UNCERTAINTY_MEASURE_WITH_UNIT`` per representation context, against one
    tolerance per vertex, per edge and per face in an OCCT B-rep, and the reader
    re-derives all of them from that single number. OCCT's default is *Average*,
    which on a lattice is pathological by construction: ~99 % of the edges are
    exactly-built interior edges at ``Precision::Confusion``, so the average
    lands near 1e-05 mm and the 1 % of boundary trims carrying real tolerance
    are silently clamped. The geometry those trims describe then changes on
    import, by an amount nothing in the file bounds.

    **Measured on the body that made this change**, `SpiralTest.step`'s dominant
    solid at ``cc=5, t=1`` -- 45,801 faces, and a solid that tessellates into
    111,618 triangles with **zero** non-manifold edges *in this process*:

    ==============  ======================  ==================
    mode            declared uncertainty    bad edges after RT
    ==============  ======================  ==================
    Average         1.E-05                  **2**
    Session 1e-4    1.E-04                  **2**
    Session 5e-4    5.E-04                  0
    **Greatest**    1.E-02                  **0**
    ==============  ======================  ==================

    The two defects sat on one 0.037 mm2 sliver whose edge records 1.512e-04 mm
    and whose pcurve strays 3.53e-05 mm -- comfortably consistent while the
    tolerance survives, and inconsistent once it is clamped to a figure below
    it. Both modes that declare *more* than 1.512e-04 fix it and both that
    declare less do not, which identifies the mechanism rather than merely
    correlating with it.

    **Why Greatest rather than a tuned constant.** It is the only setting that
    is a property of the shape instead of a guess about it: it declares what
    this body actually needs, so a body needing nothing declares nothing. The
    figure it produces is bounded by the pipeline's own repairs, which are
    capped at a tenth of ``t`` (docs/algorithm.md S8), so it can never exceed
    that -- 1e-02 mm here against a 1 mm strut, forty times below the smallest
    feature the CLI accepts.

    **docs/specification.md S11 recorded this lever as tried and useless, and
    that entry was measured against the wrong quantity** -- the count of
    edge/face pairs past their own recorded tolerance, which is one of the three
    proxies the same section later disproved as a gate. Against the gate that
    actually decides, a tessellation of the written body, it is the difference
    between a part that ships and a part that is refused.
    """
    _configure_step_units()
    Interface_Static.SetCVal_s("write.step.schema", "AP214IS")
    Interface_Static.SetCVal_s("write.step.product.name", part_name)
    writer = STEPControl_Writer()
    # **After the writer exists, deliberately.** Constructing one initialises the
    # STEP session, and that initialisation resets these statics to their
    # defaults -- so a value set beforehand is silently discarded in any process
    # that has not already touched the STEP controller. It cost a run to find:
    # the pipeline reads its input first, so the session was already up and the
    # setting held, while the same call in a fresh process wrote an Average
    # file and reported `write.precision.mode` still reading 0.
    Interface_Static.SetIVal_s("write.precision.mode", WRITE_PRECISION_GREATEST)
    status = writer.Transfer(shape, STEPControl_StepModelType.STEPControl_AsIs)
    if status != IFSelect_ReturnStatus.IFSelect_RetDone:
        raise OutputError(f"STEP transfer failed (OCCT status {status}): {path}")
    status = writer.Write(path)
    if status != IFSelect_ReturnStatus.IFSelect_RetDone:
        raise OutputError(f"STEP write failed (OCCT status {status}): {path}")
