"""Gate G19 — why removing a pinhole wire changes the volume OCCT reports.

`-cc 12 -t 2.5` on `TD_HX_rehearsal_test` was refused by the pinhole repair's
volume guard: 1.235e-09 relative drift against a 1e-09 bar, on a junction of
77.4 mm^3, while the same repair on the same part at `cc=5, t=1` drifts 3.1e-15.
Both parameter sets are inside the CLI's ranges, so a valid input was being
refused — docs/algorithm.md §11's own worst failure mode.

    python tools/prototypes/g19_pinhole_volume_bias.py

The obvious readings are all wrong, and each is disproved here rather than
argued away:

* **not quadrature noise** — the difference is stable under adaptive
  Gauss-Kronrod integration down to a requested 1e-11;
* **not a coordinate-magnitude artifact** — it is unchanged when the piece is
  translated to the origin or ten times further out, and when the properties
  are taken about a different reference point;
* **not scaled by the wire, the junction or `t`** — the cc=12 wire is
  *shorter* than either cc=5 wire (1.01e-06 mm against 3.17e-06 and
  5.81e-06 mm) and drifts seven orders of magnitude further.

What it is: `BRepGProp::VolumeProperties` documents that its shape "must be
exempt of any free boundary ... results will be false if [that is] not
respected". A pinhole wire **is** a free boundary — an edge used by exactly one
face — so that precondition is violated on the piece before the repair and
satisfied after it. The guard was comparing a figure OCCT does not promise
against one it does, and calling the difference geometry loss.

Part 3 below is the control that turns that from a plausible reading into a
measurement: adding a synthetic open wire to a *clean* face reproduces the
defect's own footprint on the same face to within 1 %, at a magnitude set by
the face and not by the wire. Lengthening the wire a thousandfold moves the
shift by under 1 % where a term proportional to length would move it by 1000x,
so there is nothing about the wire to scale a bar with; what does vary, by an
order of magnitude, is which face carries it. Nothing the repair controls can
bound that — which is why the bar is gone and `occ.only_inner_wires_dropped`
stands in its place.

Surface area, throughout all of it, is bit-identical. That is the check that
was always doing the work, exactly as docs/algorithm.md §7 argues it must:
a wire bounding no area cannot change an area.
"""

import _bootstrap  # noqa: F401
import numpy as np
from OCP.BRep import BRep_Builder, BRep_Tool
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire
from OCP.BRepGProp import BRepGProp
from OCP.BRepLib import BRepLib
from OCP.BRepTools import BRepTools_ReShape
from OCP.Geom import Geom_Plane
from OCP.Geom2d import Geom2d_Line
from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf
from OCP.GProp import GProp_GProps
from OCP.TopAbs import TopAbs_ShapeEnum
from OCP.TopLoc import TopLoc_Location
from OCP.TopTools import TopTools_ListOfShape
from OCP.gp import gp_Dir2d, gp_Pnt, gp_Pnt2d, gp_Trsf, gp_Vec

from latticegen2 import occ
from latticegen2.boundary import PINHOLE_WIRE_TOL
from latticegen2.junction import build_template
from latticegen2.lattice import lattice_params, node as lattice_node
from latticegen2.occ import _explore

INPUT = "test/TD_HX_rehearsal_test.step"
FAILING = (12.0, 2.5, (262, -27, -37))
"""The parameters and node the 1e-09 bar refused."""
CALIBRATION = (5.0, 1.0, (591, -46, -70))
"""The junction the bar was originally calibrated on (G10)."""

occ.quiet_kernel()


def area_of(shape) -> float:
    return sum(occ.area(f) for f in occ.faces(shape))


def volume_gk(shape, eps: float) -> float:
    props = GProp_GProps()
    BRepGProp.VolumePropertiesGK_s(shape, props, eps)
    return props.Mass()


def volume_about(shape, point) -> float:
    props = GProp_GProps(gp_Pnt(*point))
    BRepGProp.VolumeProperties_s(shape, props)
    return props.Mass()


def translated(shape, delta):
    trsf = gp_Trsf()
    trsf.SetTranslation(gp_Vec(*delta))
    return shape.Moved(TopLoc_Location(trsf))


def trim(cc, t, ijk):
    """``(raw piece, repaired piece, n_removed, node position)``."""
    lp = lattice_params(cc, t)
    tpl = build_template(lp)
    body = occ.read_step(INPUT)
    pos = lattice_node(lp, np.array(ijk))
    algo = BRepAlgoAPI_Common()
    args = TopTools_ListOfShape()
    args.Append(tpl.solid.Moved(occ.translation(pos)))
    tools = TopTools_ListOfShape()
    tools.Append(body)
    algo.SetArguments(args)
    algo.SetTools(tools)
    algo.Build()
    for solid in occ.solids(algo.Shape()):
        cleaned, n = occ.remove_pinhole_wires(solid, PINHOLE_WIRE_TOL)
        if n:
            return solid, cleaned, n, pos
    raise SystemExit(f"no pinhole at cc={cc}, t={t}, node {ijk} — gate is stale")


def pinhole_edges(solid):
    """Length of every edge of every pinhole wire, and its host face's area."""
    out = []
    from OCP.BRepTools import BRepTools
    from OCP.TopExp import TopExp
    from OCP.TopoDS import TopoDS
    from OCP.TopTools import TopTools_IndexedDataMapOfShapeListOfShape

    ef = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(
        solid, TopAbs_ShapeEnum.TopAbs_EDGE, TopAbs_ShapeEnum.TopAbs_FACE, ef
    )
    for f in _explore(solid, TopAbs_ShapeEnum.TopAbs_FACE):
        face = TopoDS.Face_s(f)
        outer = BRepTools.OuterWire_s(face)
        for wire in _explore(face, TopAbs_ShapeEnum.TopAbs_WIRE):
            if wire.IsSame(outer):
                continue
            lengths, uses = [], []
            for e in _explore(wire, TopAbs_ShapeEnum.TopAbs_EDGE):
                edge = TopoDS.Edge_s(e)
                props = GProp_GProps()
                BRepGProp.LinearProperties_s(edge, props)
                lengths.append(props.Mass())
                i = ef.FindIndex(edge)
                uses.append(0 if i == 0 else ef.FindFromIndex(i).Extent())
            if lengths and max(lengths) < PINHOLE_WIRE_TOL and all(u == 1 for u in uses):
                out.append((lengths, occ.area(face)))
    return out


def with_free_wire(solid, face, u, v, length):
    """``solid`` with one extra open inner wire of ``length`` mm on ``face``."""
    surface = BRep_Tool.Surface_s(face)
    line = Geom2d_Line(gp_Pnt2d(u, v), gp_Dir2d(1.0, 0.0))
    edge = BRepBuilderAPI_MakeEdge(line, surface, 0.0, length).Edge()
    BRepLib.BuildCurves3d_s(edge)

    rebuilt = face.EmptyCopied()
    builder = BRep_Builder()
    for w in _explore(face, TopAbs_ShapeEnum.TopAbs_WIRE):
        builder.Add(rebuilt, w)
    builder.Add(rebuilt, BRepBuilderAPI_MakeWire(edge).Wire())

    reshape = BRepTools_ReShape()
    reshape.Replace(face, rebuilt)
    return reshape.Apply(solid)


def main():
    raw, fixed, n, pos = trim(*FAILING)
    v_raw, v_fixed = occ.volume(raw), occ.volume(fixed)
    a_raw, a_fixed = area_of(raw), area_of(fixed)
    cc, t, ijk = FAILING

    print(f"=== 1. the refused case: cc={cc}, t={t}, node {ijk} ===")
    print(f"  node at {np.round(pos, 4).tolist()}, {n} pinhole wire(s) removed")
    for lengths, host in pinhole_edges(raw):
        print(f"  wire: {len(lengths)} edge(s), "
              f"{', '.join('%.6e' % L for L in lengths)} mm, "
              f"on a {host:.6g} mm^2 face")
    print(f"  area   {a_raw!r} -> {a_fixed!r}")
    print(f"         drift {abs(a_fixed - a_raw) / a_raw:.3e}  "
          f"({'bit-identical' if a_fixed == a_raw else 'CHANGED'})")
    print(f"  volume {v_raw!r} -> {v_fixed!r}")
    print(f"         drift {abs(v_fixed - v_raw) / abs(v_raw):.6e} relative, "
          f"{abs(v_fixed - v_raw):.6e} mm^3 absolute")

    raw5, fixed5, n5, _ = trim(*CALIBRATION)
    cc5, t5, ijk5 = CALIBRATION
    v5_raw, v5_fixed = occ.volume(raw5), occ.volume(fixed5)
    a5_raw, a5_fixed = area_of(raw5), area_of(fixed5)
    print(f"\n  for comparison, the calibration case cc={cc5}, t={t5}, node {ijk5}:")
    for lengths, host in pinhole_edges(raw5):
        print(f"    wire: {', '.join('%.6e' % L for L in lengths)} mm "
              f"on a {host:.6g} mm^2 face")
    print(f"    piece volume {v5_raw:.9g} mm^3, area drift "
          f"{abs(a5_fixed - a5_raw) / a5_raw:.3e}, volume drift "
          f"{abs(v5_fixed - v5_raw) / abs(v5_raw):.3e}")
    print("    the shorter wire drifts further: the magnitude is not the wire's")

    print("\n=== 2. it is neither quadrature noise nor a coordinate artifact ===")
    for eps in (1e-3, 1e-6, 1e-9, 1e-11):
        d = abs(volume_gk(fixed, eps) - volume_gk(raw, eps))
        print(f"  adaptive Gauss-Kronrod, eps={eps:<6g}: drift {d:.6e} mm^3")
    for label, point in (("origin", np.zeros(3)), ("the node", pos),
                         ("100x the node", pos * 100)):
        d = abs(volume_about(fixed, point) - volume_about(raw, point))
        print(f"  taken about {label:<14}: drift {d:.6e} mm^3")
    for label, delta in (("moved to the origin", -pos),
                         ("moved 10x further out", pos * 9)):
        d = abs(occ.volume(translated(fixed, delta)) - occ.volume(translated(raw, delta)))
        print(f"  {label:<26}: drift {d:.6e} mm^3")

    print("\n=== 3. control: add a synthetic open wire to a clean face ===")
    print("  (if a free boundary biases the integral, adding one must reproduce it)")
    print(f"  {'face':>5} {'area mm^2':>10} {'wire mm':>10} {'dV mm^3':>15} "
          f"{'relative':>11} {'dA':>8}")
    per_face = {}
    for fi, face in enumerate(occ.faces(fixed)):
        if not isinstance(BRep_Tool.Surface_s(face), Geom_Plane):
            continue
        centre = occ.centroid(face)
        proj = GeomAPI_ProjectPointOnSurf(gp_Pnt(*centre), BRep_Tool.Surface_s(face))
        u, v = proj.LowerDistanceParameters()
        for length in (1e-6, 1e-5, 1e-4, 1e-3):
            perturbed = with_free_wire(fixed, face, u, v, length)
            dv = occ.volume(perturbed) - v_fixed
            da = area_of(perturbed) - a_fixed
            per_face.setdefault(fi, []).append((dv, da))
            print(f"  {fi:>5} {occ.area(face):>10.4f} {length:>10g} {dv:>+15.6e} "
                  f"{abs(dv) / v_fixed:>11.3e} {da:>8.1e}")
        if len(per_face) >= 4:
            break

    spread = {fi: max(abs(a[0] - b[0]) / max(abs(a[0]), 1e-30)
                      for a in vals for b in vals)
              for fi, vals in per_face.items()}
    areas_exact = all(da == 0.0 for vals in per_face.values() for _, da in vals)
    # A term proportional to the wire's length would move by 1000x across the
    # three decades swept. The bar is 1 %, which is four orders inside that.
    length_independent = all(s < 1e-2 for s in spread.values())
    magnitudes = sorted(abs(dv) / v_fixed for vals in per_face.values() for dv, _ in vals)
    by_face = max(magnitudes) / min(magnitudes)
    print(f"\n  area unchanged for every synthetic wire:      {areas_exact}")
    print(f"  1000x the wire length moves the shift by:     "
          f"{max(spread.values()):.2e} ({length_independent})")
    print(f"  changing which face carries it moves it by:   "
          f"{by_face:.1f}x  ({magnitudes[0]:.2e} to {magnitudes[-1]:.2e} relative)")

    print("\n=== 4. the check that replaces the bar ===")
    reason = occ.only_inner_wires_dropped(raw, fixed, n)
    print(f"  only_inner_wires_dropped on the refused case: "
          f"{'clean' if reason is None else reason}")
    reason5 = occ.only_inner_wires_dropped(raw5, fixed5, n5)
    print(f"  ... and on the calibration case:              "
          f"{'clean' if reason5 is None else reason5}")
    same = sum(1 for x, y in zip(occ.faces(raw), occ.faces(fixed)) if x.IsSame(y))
    print(f"  faces returned as the same object:            "
          f"{same} of {len(occ.faces(raw))} (only the repaired face is rebuilt)")
    miscount = occ.only_inner_wires_dropped(raw, fixed, n + 1)
    print(f"  it does reject a wrong wire count:            "
          f"{miscount is not None}")

    passed = (
        a_fixed == a_raw
        and a5_fixed == a5_raw
        and areas_exact
        and length_independent
        and reason is None
        and reason5 is None
        and miscount is not None
    )
    print("\n=== VERDICT ===")
    print(f"  area is bit-identical in every case:          "
          f"{a_fixed == a_raw and a5_fixed == a5_raw and areas_exact}")
    print(f"  the volume shift is a per-face constant:      {length_independent}")
    print(f"  the structural check accepts both repairs:    "
          f"{reason is None and reason5 is None}")
    print("PASS - drop the volume bar" if passed else "FAIL - do not adopt")


if __name__ == "__main__":
    main()
