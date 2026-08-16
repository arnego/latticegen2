"""G12: the residual invalid boundary faces are falsely self-intersecting wires.

Establishes the mechanism behind the 4 faces PR #17's vertex-tolerance repair
left invalid on the ``cc=5, t=1`` rehearsal, and measures every candidate
repair against it. Results and reasoning: ``tools/prototypes/RESULTS.md`` G12.

Run against the faces lifted out of that run:

    python tools/prototypes/g12_self_intersecting_wire.py \
        --faces diagnostics/residual-invalid-faces

Falls back to the two committed fixtures in ``test/`` when ``--faces`` is not
given, so the gate is reproducible from a clean checkout.
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from OCP.BRep import BRep_Builder, BRep_Tool  # noqa: E402
from OCP.BRepCheck import BRepCheck_Analyzer, BRepCheck_Face, BRepCheck_Wire  # noqa: E402
from OCP.BRepLib import BRepLib  # noqa: E402
from OCP.ShapeFix import ShapeFix_Shape, ShapeFix_Wire  # noqa: E402
from OCP.TopAbs import TopAbs_ShapeEnum  # noqa: E402
from OCP.TopExp import TopExp_Explorer  # noqa: E402
from OCP.TopoDS import TopoDS, TopoDS_Edge  # noqa: E402

from latticegen2 import occ  # noqa: E402
from latticegen2.parallel import read_brep  # noqa: E402

FALLBACK = [
    ("self-intersecting-wire-cylinder", "test/self-intersecting-wire-cylinder.brep"),
    ("self-intersecting-wire-bspline", "test/self-intersecting-wire-bspline.brep"),
]


def explore(shape, kind):
    exp = TopExp_Explorer(shape, kind)
    out = []
    while exp.More():
        out.append(exp.Current())
        exp.Next()
    return out


def face_of(path):
    return TopoDS.Face_s(explore(read_brep(path), TopAbs_ShapeEnum.TopAbs_FACE)[0])


def self_intersect(face):
    """``(status_name, edge_a, edge_b)`` from ``BRepCheck_Wire``."""
    for w in explore(face, TopAbs_ShapeEnum.TopAbs_WIRE):
        a, b = TopoDS_Edge(), TopoDS_Edge()
        st = BRepCheck_Wire(TopoDS.Wire_s(w)).SelfIntersect(face, a, b, False)
        return str(st).split(".")[-1], a, b
    return "no wire", None, None


def topo(face):
    return (
        [e for e in explore(face, TopAbs_ShapeEnum.TopAbs_EDGE)],
        [v for v in explore(face, TopAbs_ShapeEnum.TopAbs_VERTEX)],
    )


def kept(before, after):
    be, bv = before
    ae, av = after
    e = sum(1 for x in ae if any(x.IsSame(y) for y in be))
    v = sum(1 for x in av if any(x.IsSame(y) for y in bv))
    return f"edges {e}/{len(be)}, verts {v}/{len(bv)}"


def describe(name, path):
    face = face_of(path)
    st, a, b = self_intersect(face)
    chk = BRepCheck_Face(face)
    print(f"\n### {name}")
    print(f"  face valid                 : {BRepCheck_Analyzer(face).IsValid()}")
    print(f"  area                       : {occ.area(face):.12f} mm^2")
    print(f"  wires / edges              : "
          f"{len(explore(face, TopAbs_ShapeEnum.TopAbs_WIRE))} / "
          f"{len(explore(face, TopAbs_ShapeEnum.TopAbs_EDGE))}")
    for check in ("IntersectWires", "ClassifyWires", "OrientationOfWires"):
        print(f"  {check:27s}: {str(getattr(chk, check)()).split('.')[-1]}")
    bad_e = [
        i for i, e in enumerate(explore(face, TopAbs_ShapeEnum.TopAbs_EDGE))
        if not BRepCheck_Analyzer(e).IsValid()
    ]
    bad_v = [
        i for i, v in enumerate(explore(face, TopAbs_ShapeEnum.TopAbs_VERTEX))
        if not BRepCheck_Analyzer(v).IsValid()
    ]
    print(f"  standalone-invalid sub     : {len(bad_e)} edges, {len(bad_v)} vertices")
    print(f"  BRepCheck_Wire.SelfIntersect: {st}")
    if a is None or st.endswith("NoError"):
        return
    edges = [TopoDS.Edge_s(e) for e in explore(face, TopAbs_ShapeEnum.TopAbs_EDGE)]

    def idx(x):
        return next((i for i, e in enumerate(edges) if e.IsSame(x)), None)

    print(f"    reported pair            : e{idx(a)} <-> e{idx(b)}")
    for tag, x in (("A", a), ("B", b)):
        print(f"    {tag} tol={BRep_Tool.Tolerance_s(x):.6e} "
              f"3D curve={type(BRep_Tool.Curve_s(x, 0.0, 0.0)).__name__}")
    for v in occ._shared_vertices(a, b):
        p = BRep_Tool.Pnt_s(v)
        print(f"    shared vertex            : tol={BRep_Tool.Tolerance_s(v):.6e} "
              f"at [{p.X():.5f}, {p.Y():.5f}, {p.Z():.5f}]")


def sweep(name, path):
    """Which tolerance does the check key on? Widening moves no geometry."""
    print(f"\n### {name}: tolerance sweep (c=check clean, S=still self-intersecting)")
    for target in ("shared vertex", "fat edge"):
        cells = []
        for mult in (1.05, 1.1, 1.25, 1.5, 2.0, 5.0):
            face = face_of(path)
            st, a, b = self_intersect(face)
            if st.endswith("NoError"):
                cells.append("n/a")
                continue
            a0 = occ.area(face)
            builder = BRep_Builder()
            if target == "shared vertex":
                for v in occ._shared_vertices(a, b):
                    builder.UpdateVertex(v, BRep_Tool.Tolerance_s(v) * mult)
            else:
                fat = max((a, b), key=BRep_Tool.Tolerance_s)
                builder.UpdateEdge(fat, BRep_Tool.Tolerance_s(fat) * mult)
            st2, _, _ = self_intersect(face)
            ok = BRepCheck_Analyzer(face).IsValid()
            moved = occ.area(face) != a0
            cells.append(
                f"x{mult:g}:{'c' if st2.endswith('NoError') else 'S'}"
                f"{'V' if ok else '.'}{'~' if moved else ''}"
            )
        print(f"  {target:16s} " + "  ".join(f"{c:>9s}" for c in cells))


def repairs(name, path):
    print(f"\n### {name}: candidate repairs")

    def r_widen(f):
        occ._widen_self_intersection_vertices(f)
        return f

    def r_sameparam(f):
        BRepLib.SameParameter_s(f, 1e-7, True)
        return f

    def r_wirefix(f):
        wire = TopoDS.Wire_s(explore(f, TopAbs_ShapeEnum.TopAbs_WIRE)[0])
        fx = ShapeFix_Wire(wire, f, 1e-7)
        fx.FixSelfIntersectionMode = True
        fx.FixSelfIntersection()
        return f

    def r_shapefix(f):
        fx = ShapeFix_Shape(f)
        fx.Perform()
        return TopoDS.Face_s(fx.Shape())

    for label, fn in (
        ("occ._widen_self_intersection_vertices", r_widen),
        ("occ.fix_vertex_tolerances", lambda f: (occ.fix_vertex_tolerances([f]), f)[1]),
        ("BRepLib.SameParameter", r_sameparam),
        ("ShapeFix_Wire.FixSelfIntersection", r_wirefix),
        ("ShapeFix_Shape.Perform", r_shapefix),
    ):
        face = face_of(path)
        a0, t0 = occ.area(face), topo(face)
        try:
            out = fn(face)
        except Exception as exc:  # noqa: BLE001
            print(f"  {label:38s} RAISED {type(exc).__name__}")
            continue
        drift = abs(occ.area(out) - a0) / a0
        print(
            f"  {label:38s} valid={str(BRepCheck_Analyzer(out).IsValid()):5s} "
            f"area drift={drift:.3e}  kept {kept(t0, topo(out))}"
        )


def control(path):
    """Negative control: the repair must be a no-op on a face that is valid."""
    print("\n### negative control (a valid face)")
    face = face_of(path)
    before = occ.area(face)
    occ.fix_vertex_tolerances([face])
    print(f"  valid before/after         : True / {BRepCheck_Analyzer(face).IsValid()}")
    print(f"  area drift                 : {abs(occ.area(face) - before):.3e}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--faces", help="directory of residual_*.brep")
    ap.add_argument("--control", help="a .brep holding a valid face")
    args = ap.parse_args()

    if args.faces:
        cases = [
            (f"residual_{i}", os.path.join(args.faces, f"residual_{i}.brep"))
            for i in range(4)
        ]
    else:
        cases = [(n, p) for n, p in FALLBACK if os.path.exists(p)]
        if not cases:
            print("No faces to measure: pass --faces, or run from the repo root.")
            return 1

    print("G12 -- residual invalid boundary faces: falsely self-intersecting wires")
    for name, path in cases:
        describe(name, path)
    for name, path in cases:
        sweep(name, path)
    for name, path in cases:
        repairs(name, path)
    if args.control:
        control(args.control)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
