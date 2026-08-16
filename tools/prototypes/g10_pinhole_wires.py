"""Gate G10 — what the rehearsal's last 2 open edges actually are, and what removes them.

docs/specification.md §10 tracked these as "micron-scale debris edges" and
proposed `ShapeFix` small-edge removal as the last untried candidate. Both the
diagnosis and the proposal turn out to be wrong, and this gate records why,
because the wrong version is very convincing: the edges really are
non-degenerate and really are ~3e-06 mm long, so a synthetic reproduction
matching those two properties passes every test you can think of while
repairing a different defect.

    python tools/prototypes/g10_pinhole_wires.py

Run against the real part at the real node — no synthetic stand-in, on purpose.

What they actually are: each is an **inner wire** of a planar face, made of a
single edge whose two endpoints do not meet. It bounds no area at all. It is a
pinhole, not a sliver. That is why it is used by exactly one face, and why
`weld.shell_defects` rejects the shell for it while `BRepCheck_Analyzer` calls
the same solid valid.

Why OCCT's own repairs decline: `ShapeFix_Wireframe` looks for small *edges*
between two faces, and `ShapeFix_Face.FixSmallAreaWire` for small *wires* with
an area to measure. Both want a well-formed closed wire; a single non-closing
edge is not one.

The fix is therefore explicit and narrow (`occ.remove_pinhole_wires`): drop a
non-outer wire when every edge in it is shorter than the bar **and** every edge
is already used exactly once. The second condition is what makes the repair
unable to open a hole — it only ever deletes edges that are already unpaired —
and it is what this gate checks last, by running at a bar larger than every
edge in a closed solid and confirming nothing is touched.
"""

import _bootstrap  # noqa: F401
import numpy as np
from OCP.BRep import BRep_Tool
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common
from OCP.BRepGProp import BRepGProp
from OCP.BRepTools import BRepTools
from OCP.GProp import GProp_GProps
from OCP.ShapeFix import ShapeFix_Face, ShapeFix_Wireframe
from OCP.TopAbs import TopAbs_ShapeEnum
from OCP.TopoDS import TopoDS
from OCP.TopTools import (
    TopTools_DataMapOfShapeListOfShape,
    TopTools_ListOfShape,
    TopTools_MapOfShape,
)

from latticegen2 import occ, weld
from latticegen2.boundary import PINHOLE_WIRE_TOL, _open_shell
from latticegen2.junction import build_template, is_cap_plane_face
from latticegen2.lattice import lattice_params, node as lattice_node

CC, T = 5.0, 1.0
NODE = (591, -46, -70)
INPUT = "test/TD_HX_rehearsal_test.step"

occ.quiet_kernel()
lp = lattice_params(CC, T)
tpl = build_template(lp)


def edge_length(edge) -> float:
    props = GProp_GProps()
    BRepGProp.LinearProperties_s(edge, props)
    return props.Mass()


def real_edges(shape):
    return [TopoDS.Edge_s(e) for e in occ._explore(shape, TopAbs_ShapeEnum.TopAbs_EDGE)
            if not BRep_Tool.Degenerated_s(TopoDS.Edge_s(e))]


def cap_areas(solid, pos):
    per = {}
    for f in occ.faces(solid):
        h = is_cap_plane_face(lp, f, pos)
        if h is not None:
            per[h] = per.get(h, 0.0) + occ.area(f)
    return per


def main() -> None:
    pos = lattice_node(lp, np.array(NODE))
    body = occ.read_step(INPUT)
    algo = BRepAlgoAPI_Common()
    args = TopTools_ListOfShape()
    args.Append(tpl.solid.Moved(occ.translation(pos)))
    tools = TopTools_ListOfShape()
    tools.Append(body)
    algo.SetArguments(args)
    algo.SetTools(tools)
    algo.Build()
    solid = occ.solids(algo.Shape())[0]

    faces0, edges0 = occ.count_subshapes(solid)
    vol0 = occ.volume(solid)
    area0 = sum(occ.area(f) for f in occ.faces(solid))
    caps0 = cap_areas(solid, pos)
    print(f"cc={CC} t={T}, node {NODE}")
    print(f"raw trim: {faces0} faces, {edges0} edges, volume {vol0:.9f}, "
          f"area {area0:.9f}")
    print(f"  BRepCheck_Analyzer says valid: {occ.is_valid(solid)}")
    print(f"  shell_defects says: {weld.shell_defects(_open_shell(occ.faces(solid)))[:2]}"
          "   <-- (open, misoriented); the defect, on the piece alone")

    # --- 1. what they are -------------------------------------------------
    print("\n=== 1. the offending edges ===")
    for edge in real_edges(solid):
        L = edge_length(edge)
        if L >= PINHOLE_WIRE_TOL:
            continue
        pts = [BRep_Tool.Pnt_s(TopoDS.Vertex_s(v))
               for v in occ._explore(edge, TopAbs_ShapeEnum.TopAbs_VERTEX)]
        gap = 0.0
        if len(pts) == 2:
            gap = np.linalg.norm(np.array([pts[0].X() - pts[1].X(),
                                           pts[0].Y() - pts[1].Y(),
                                           pts[0].Z() - pts[1].Z()]))
        print(f"  length {L:.6e} mm, endpoints {gap:.3e} mm apart "
              f"-> {'closes' if gap < 1e-12 else 'does NOT close'}")
        for fi, f in enumerate(occ.faces(solid)):
            face = TopoDS.Face_s(f)
            outer = BRepTools.OuterWire_s(face)
            for w in occ._explore(face, TopAbs_ShapeEnum.TopAbs_WIRE):
                if any(e.IsSame(edge) for e in occ._explore(w, TopAbs_ShapeEnum.TopAbs_EDGE)):
                    n = sum(1 for _ in occ._explore(w, TopAbs_ShapeEnum.TopAbs_EDGE))
                    print(f"    on face #{fi} (area {occ.area(face):.6f} mm^2), "
                          f"in a {n}-edge "
                          f"{'OUTER' if w.IsSame(outer) else 'INNER'} wire")

    # --- 2. OCCT's own repairs decline ------------------------------------
    print("\n=== 2. what OCCT's repair tools do about it ===")
    for prec in (1e-5, 3e-5, 1e-4, 1e-3, 1e-2):
        sfw = ShapeFix_Wireframe(solid)
        sfw.SetPrecision(prec)
        small = TopTools_MapOfShape()
        e2f = TopTools_DataMapOfShapeListOfShape()
        fws = TopTools_DataMapOfShapeListOfShape()
        multi = TopTools_MapOfShape()
        sfw.CheckSmallEdges(small, e2f, fws, multi)
        sfw.SetLimitAngle(-1.0)
        sfw.ModeDropSmallEdges = True
        sfw.FixSmallEdges()
        after = occ.count_subshapes(sfw.Shape())[1]
        print(f"  ShapeFix_Wireframe(prec={prec:<7g}): CheckSmallEdges finds "
              f"{small.Extent()} candidate(s), edges {edges0}->{after}")
    for prec in (1e-4, 1e-2):
        removed = 0
        for f in occ.faces(solid):
            face = TopoDS.Face_s(f)
            before = sum(1 for _ in occ._explore(face, TopAbs_ShapeEnum.TopAbs_WIRE))
            sff = ShapeFix_Face(face)
            sff.SetPrecision(prec)
            sff.FixSmallAreaWireMode = True
            sff.Perform()
            removed += before - sum(1 for _ in occ._explore(sff.Face(),
                                                            TopAbs_ShapeEnum.TopAbs_WIRE))
        print(f"  ShapeFix_Face.FixSmallAreaWire(prec={prec:<7g}): "
              f"{removed} wire(s) removed")

    # --- 3. the fix --------------------------------------------------------
    print("\n=== 3. occ.remove_pinhole_wires ===")
    fixed, n_removed = occ.remove_pinhole_wires(solid, PINHOLE_WIRE_TOL)
    faces1, edges1 = occ.count_subshapes(fixed)
    area1 = sum(occ.area(f) for f in occ.faces(fixed))
    vol1 = occ.volume(fixed)
    caps1 = cap_areas(fixed, pos)
    cap_drift = max((abs(caps1.get(h, 0.0) - caps0.get(h, 0.0)) / max(caps0.get(h, 0.0), T * T)
                     for h in set(caps0) | set(caps1)), default=0.0)
    print(f"  removed {n_removed} pinhole wire(s) at bar {PINHOLE_WIRE_TOL:g} mm")
    print(f"  faces {faces0}->{faces1}, edges {edges0}->{edges1}")
    print(f"  surface area drift: {abs(area1 - area0) / area0:.3e}")
    print(f"  volume drift:       {abs(vol1 - vol0) / vol0:.3e}")
    print(f"  cap set unchanged:  {set(caps0) == set(caps1)}, worst cap drift {cap_drift:.3e}")
    print(f"  still valid:        {occ.is_valid(fixed)}")
    print(f"  shell_defects:      {weld.shell_defects(_open_shell(occ.faces(fixed)))[:2]}")

    # --- 4. the safety property -------------------------------------------
    print("\n=== 4. a paired wire is never removed, at any bar ===")
    huge = 10.0 * lp.a
    _, n_bad = occ.remove_pinhole_wires(tpl.solid, huge)
    print(f"  template solid at bar {huge:g} mm (larger than every edge in it): "
          f"{n_bad} wire(s) removed")

    defects = weld.shell_defects(_open_shell(occ.faces(fixed)))[:2]
    passed = (
        n_removed == 2
        and defects == (0, 0)
        and faces1 == faces0
        and abs(area1 - area0) == 0.0
        and cap_drift == 0.0
        and occ.is_valid(fixed)
        and n_bad == 0
    )
    print("\n=== VERDICT ===")
    print(f"  both pinholes removed, piece closes:  {n_removed == 2 and defects == (0, 0)}")
    print(f"  no face lost, area exactly unchanged: {faces1 == faces0 and abs(area1 - area0) == 0.0}")
    print(f"  caps untouched:                       {cap_drift == 0.0}")
    print(f"  paired wires out of reach:            {n_bad == 0}")
    print("PASS - adopt" if passed else "FAIL - do not adopt")


if __name__ == "__main__":
    main()
