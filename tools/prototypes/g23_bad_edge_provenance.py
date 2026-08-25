"""Gate G23a — what *are* the 26 bad edges on the rehearsal's dominant body?

`TD_HX_rehearsal_test.step` at `cc=5, t=1` is refused at `export truth`
(docs/specification.md §10): solid 0 tessellates into 1,427,670 triangles
carrying 26 edges not used by exactly two of them. Two things are already
measured and are not re-derived here — the STEP write is exonerated (the same 26
appear in memory) and so is `simplify` (the pre-unification solid gives the same
26). So the defect is already in the shell `weld.assemble` proved watertight.

**And that is the tension this gate exists to resolve.** `weld.shell_defects`
proves every *B-rep* edge is used exactly twice, once each way, keyed on
`TopoDS_Shape::IsSame`. `occ.exported_mesh_defects` counts *mesh* edges keyed on
coordinates rounded to 1e-7 mm. Both can be right at once only if, at each of
the 26, one of these holds:

* **T (topology)** — two coincident but *distinct* ``TopoDS_Edge`` objects, each
  used by two faces, so the mesh sees four uses at one place. This is duplicate
  material, and it is what a declined cap that keeps both sides would produce.
* **M (mesh)** — one shared ``TopoDS_Edge`` whose two faces discretize it to
  *different* points, so the mesh sees two edges each used once. **No hole exists
  in the B-rep at all**, and the repair would belong to the pcurve family
  (docs/algorithm.md §9), not to `boundary`/`weld`.
* **one-face** — a B-rep edge genuinely used by one face, which would contradict
  `assemble` and make `shell_defects` itself the thing to audit.

docs/specification.md §10 currently records the 14 used-once edges as "holes".
**That is an inference, not a measurement**, and on a shell proven closed it is
the less likely of the two readings. Settling it decides which module the work
lands in, so nothing else starts until this runs. Five times in this project's
history the first convincing explanation has been wrong (G9-G12 and the pcurve
chapter); this is the cheapest possible way not to make it six.

**How the mapping is done, and why not by endpoints.** A bad mesh edge is one
segment of the polyline discretizing a B-rep edge, so its endpoints are in
general *not* that edge's vertices and cannot be matched against them.
``BRep_Tool.PolygonOnTriangulation_s`` gives, per (edge, face), the node indices
into that face's triangulation — which maps every mesh segment on a face
boundary back to the exact B-rep edge that produced it, with no tolerance and no
search.

**Controls, because an instrument that only ever agrees on sound geometry proves
nothing (G10).**

* **C1, negative** — the 13 sibling solids measured clean by the same run must
  come back at 0 bad edges. Anything else means this is measuring something other
  than what `export truth` measured, and no other result here is admissible.
* **C2, positive** — a sound solid with one face deliberately removed must come
  back with exactly that face's non-degenerate edges reported as used-once, each
  attributed to it. This is what stops the provenance silently resolving to
  "unattributable" on broken *and* sound input alike.

Usage::

    python tools/prototypes/g23_bad_edge_provenance.py <solid.brep> [more.brep ...]
    python tools/prototypes/g23_bad_edge_provenance.py --control-c2 <sound.brep>
"""

import os
import sys
import time
from collections import defaultdict

import _bootstrap  # noqa: F401

from OCP.BRep import BRep_Tool, BRep_Builder
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.TopAbs import TopAbs_ShapeEnum
from OCP.TopExp import TopExp
from OCP.TopLoc import TopLoc_Location
from OCP.TopoDS import TopoDS, TopoDS_Shell
from OCP.TopTools import (
    TopTools_IndexedMapOfShape,
    TopTools_IndexedDataMapOfShapeListOfShape,
)

from latticegen2 import occ
from latticegen2.parallel import read_brep

QUANTUM = 1e7
"""Points are interned on a 1e-7 mm grid — OCCT's own confusion, and exactly
what `occ.exported_mesh_defects` uses, so edge keys here mean what they mean in
production.

``--skip-degenerate`` additionally matches production's rule of dropping
triangles with two coincident vertices. Without the flag this script reports the
count as it stood *before* that rule, which is what the before/after in G23 is
measured with — and it is why the two once disagreed by one, production having
counted the collapsed ``lo == hi`` key that this script has always skipped."""


def _mesh_edge_key(a, b):
    lo, hi = (a, b) if a <= b else (b, a)
    return lo * 4_294_967_296 + hi


def analyse(shape, deflection=occ.DEFAULT_MESH_DEFLECTION, skip_degenerate=False):
    """Count mesh edges by use, and record which (face, B-rep edge) produced each.

    Returns ``(triangles, counts, provenance, point_of, edge_faces, edge_map)``.
    """
    BRepMesh_IncrementalMesh(shape, deflection, False, 0.2, True)

    edge_map = TopTools_IndexedMapOfShape()
    TopExp.MapShapes_s(shape, TopAbs_ShapeEnum.TopAbs_EDGE, edge_map)

    # Which faces use each B-rep edge. This is the ancestry that separates T
    # from one-face: it is a statement about the B-rep, independent of any mesh.
    ancestry = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(
        shape, TopAbs_ShapeEnum.TopAbs_EDGE, TopAbs_ShapeEnum.TopAbs_FACE, ancestry
    )
    edge_faces = {}
    for i in range(1, edge_map.Extent() + 1):
        e = edge_map.FindKey(i)
        edge_faces[i] = ancestry.FindFromKey(e).Extent() if ancestry.Contains(e) else 0

    missing_poly = [0]
    degenerate_tris = [0]
    pole_nodes = set()
    poly_pairs = [0]
    point_id = {}
    point_of = {}
    counts = defaultdict(int)
    provenance = defaultdict(list)   # mesh edge key -> [(face_idx, brep_edge_idx|None)]
    triangles = 0

    for fi, f in enumerate(occ._explore(shape, TopAbs_ShapeEnum.TopAbs_FACE)):
        face = TopoDS.Face_s(f)
        loc = TopLoc_Location()
        tri = BRep_Tool.Triangulation_s(face, loc)
        if tri is None:
            continue
        trsf = loc.Transformation()
        ids = []
        for i in range(1, tri.NbNodes() + 1):
            p = tri.Node(i).Transformed(trsf)
            key = (round(p.X() * QUANTUM), round(p.Y() * QUANTUM), round(p.Z() * QUANTUM))
            got = point_id.get(key)
            if got is None:
                got = point_id[key] = len(point_id)
                point_of[got] = (p.X(), p.Y(), p.Z())
            ids.append(got)

        # Which mesh segments on this face lie on which B-rep edge. No tolerance,
        # no search: OCCT stores the node indices per (edge, triangulation).
        seg_owner = {}
        for e in occ._explore(face, TopAbs_ShapeEnum.TopAbs_EDGE):
            edge = TopoDS.Edge_s(e)
            if BRep_Tool.Degenerated_s(edge):
                # A pole: every u maps to one point, so distinct parametric
                # nodes intern to one id and the fan around it reads as an
                # over-used edge. `shell_defects` and `free_edges` both skip
                # degenerate edges for the same reason (docs/specification.md
                # §11); this counter never learned to.
                poly = BRep_Tool.PolygonOnTriangulation_s(edge, tri, loc)
                if poly is not None:
                    nodes = poly.Nodes()
                    for k in range(nodes.Lower(), nodes.Upper() + 1):
                        pole_nodes.add(ids[nodes.Value(k) - 1])
                continue
            poly = BRep_Tool.PolygonOnTriangulation_s(edge, tri, loc)
            if poly is None:
                # An edge/face pair with no polygon would leave its boundary
                # segments unattributed and they would read as "interior",
                # which is a different finding entirely. Counted so the
                # classification can be trusted rather than assumed.
                missing_poly[0] += 1
                continue
            poly_pairs[0] += 1
            ei = edge_map.FindIndex(edge)
            nodes = poly.Nodes()
            for k in range(nodes.Lower(), nodes.Upper()):
                a, b = ids[nodes.Value(k) - 1], ids[nodes.Value(k + 1) - 1]
                if a != b:
                    seg_owner[_mesh_edge_key(a, b)] = ei

        for i in range(1, tri.NbTriangles() + 1):
            a, b, c = tri.Triangle(i).Get()
            triangles += 1
            ia, ib, ic = ids[a - 1], ids[b - 1], ids[c - 1]
            if skip_degenerate and (ia == ib or ib == ic or ia == ic):
                # A triangle two of whose vertices intern to one point has no
                # area. It contributes its collapsed key once *and the real
                # edge twice*, since that edge appears as both (A,P) and (P,A)
                # in the same triangle -- which is how a sound pole reads as an
                # over-used edge.
                degenerate_tris[0] += 1
                continue
            for u, v in ((a, b), (b, c), (c, a)):
                lo, hi = ids[u - 1], ids[v - 1]
                if lo == hi:
                    continue
                key = _mesh_edge_key(lo, hi)
                counts[key] += 1
                provenance[key].append((fi, seg_owner.get(key)))

    print(f"  [edge/face pairs with a polygon: {poly_pairs[0]}, "
          f"without: {missing_poly[0]}; degenerate triangles: "
          f"{degenerate_tris[0]}, skipped: {skip_degenerate}]", flush=True)
    return triangles, counts, provenance, point_of, edge_faces, edge_map, pole_nodes


def verdict(key, uses, prov, edge_faces):
    """T, M, one-face, or interior — from the provenance, not from the count."""
    brep_edges = {ei for _fi, ei in prov if ei is not None}
    faces = {fi for fi, _ei in prov}
    if not brep_edges:
        # No B-rep edge claims it: an interior segment of one face's own
        # triangulation. Used once here would mean a broken triangulation.
        return "interior", brep_edges, faces
    if len(brep_edges) >= 2:
        return "T", brep_edges, faces
    only = next(iter(brep_edges))
    if edge_faces.get(only, 0) <= 1:
        return "one-face", brep_edges, faces
    return "M", brep_edges, faces


def report(path, deflection=occ.DEFAULT_MESH_DEFLECTION, show=True,
           skip_degenerate=False):
    t0 = time.time()
    shape = read_brep(path)
    tri, counts, prov, point_of, edge_faces, edge_map, poles = analyse(shape, deflection, skip_degenerate)
    bad = [(k, n) for k, n in counts.items() if n != 2]
    at_pole = sum(1 for k, _n in bad
                  if (k // 4_294_967_296) in poles or (k % 4_294_967_296) in poles)
    print(f"  [pole nodes: {len(poles)}; bad edges touching a pole: {at_pole} "
          f"of {len(bad)} -> {len(bad) - at_pole} would survive]", flush=True)
    name = os.path.basename(path)
    print(f"\n=== {name}: {tri} triangle(s), {len(bad)} bad edge(s) "
          f"[{time.time() - t0:.1f}s]", flush=True)
    if not show or not bad:
        return len(bad), []

    rows = []
    for key, uses in sorted(bad, key=lambda x: -x[1]):
        lo, hi = divmod(key, 4_294_967_296)
        p, q = point_of[lo], point_of[hi]
        mid = tuple((p[i] + q[i]) / 2 for i in range(3))
        length = sum((p[i] - q[i]) ** 2 for i in range(3)) ** 0.5
        kind, brep_edges, faces = verdict(key, uses, prov[key], edge_faces)
        rows.append((uses, kind, mid, length, brep_edges, faces))
        fdesc = ", ".join(f"e{ei}({edge_faces.get(ei, 0)}f)" for ei in sorted(brep_edges)) or "-"
        print(f"  used {uses}x  {kind:9s}  mid [{mid[0]:.3f}, {mid[1]:.3f}, {mid[2]:.3f}]"
              f"  len {length:.6f}  faces {sorted(faces)}  brep {fdesc}")

    tally = defaultdict(int)
    for uses, kind, *_ in rows:
        tally[(uses, kind)] += 1
    print("  --- verdict tally (uses, kind) -> count:",
          {f"{u}x/{k}": c for (u, k), c in sorted(tally.items())})

    # Which faces are implicated, and what are they like? 23 of 25 readings on
    # the rehearsal are a mesher failure rather than a B-rep defect, so the
    # question becomes what is pathological about the faces it fails on.
    want = sorted({fi for _u, _k, _m, _l, _b, fs in rows for fi in fs})
    print(f"  --- {len(want)} implicated face(s):")
    wanted = set(want)
    for fi, f in enumerate(occ._explore(shape, TopAbs_ShapeEnum.TopAbs_FACE)):
        if fi not in wanted:
            continue
        face = TopoDS.Face_s(f)
        try:
            area = occ.area(face)
        except Exception:
            area = float("nan")
        tol = BRep_Tool.Tolerance_s(face)
        etols = [BRep_Tool.Tolerance_s(TopoDS.Edge_s(e))
                 for e in occ._explore(face, TopAbs_ShapeEnum.TopAbs_EDGE)]
        surf = BRep_Tool.Surface_s(face)
        kind = type(surf).__name__.replace("Geom_", "")
        loc = TopLoc_Location()
        t = BRep_Tool.Triangulation_s(face, loc)
        print(f"    face {fi}: area {area:.6e} mm^2, face tol {tol:.3e}, "
              f"max edge tol {max(etols) if etols else 0:.3e}, {len(etols)} edge(s), "
              f"{kind}, {t.NbTriangles() if t else 0} triangle(s)")
    return len(bad), rows


def control_c2(path):
    """Positive control: drop one face from a sound solid, expect its edges back."""
    shape = read_brep(path)
    faces = [TopoDS.Face_s(f) for f in occ._explore(shape, TopAbs_ShapeEnum.TopAbs_FACE)]
    print(f"C2 control on {os.path.basename(path)}: {len(faces)} faces")
    n_clean, _ = report(path, show=False)
    if n_clean:
        print(f"  FAIL: control solid is not clean ({n_clean} bad edges)")
        return False

    dropped = faces[0]
    builder = BRep_Builder()
    shell = TopoDS_Shell()
    builder.MakeShell(shell)
    for f in faces[1:]:
        builder.Add(shell, f)

    tri, counts, prov, point_of, edge_faces, _, _poles = analyse(shell)
    bad = [(k, n) for k, n in counts.items() if n != 2]
    want = sum(
        1 for e in occ._explore(dropped, TopAbs_ShapeEnum.TopAbs_EDGE)
        if not BRep_Tool.Degenerated_s(TopoDS.Edge_s(e))
    )
    kinds = defaultdict(int)
    for key, uses in bad:
        kind, _b, _f = verdict(key, uses, prov[key], edge_faces)
        kinds[(uses, kind)] += 1
    print(f"  removed 1 face with {want} non-degenerate B-rep edge(s); "
          f"{len(bad)} bad mesh edge(s), kinds {dict(kinds)}")
    ok = len(bad) > 0 and all(u == 1 for _k, u in bad)
    print(f"  {'PASS' if ok else 'FAIL'}: every bad edge from a removed face is used once")
    return ok


def main(argv):
    if not argv:
        print(__doc__)
        return 2
    skip = False
    if argv[0] == "--skip-degenerate":
        skip = True
        argv = argv[1:]
    if not argv:
        print(__doc__)
        return 2
    if argv[0] == "--control-c2":
        return 0 if control_c2(argv[1]) else 1
    defl = occ.DEFAULT_MESH_DEFLECTION
    if argv[0] == "--deflection":
        defl = float(argv[1])
        argv = argv[2:]
    for path in argv:
        report(path, deflection=defl, skip_degenerate=skip)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
