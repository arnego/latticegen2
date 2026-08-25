"""Gate G24 - what are the 13 readings that survive the pole fix?

`TD_HX_rehearsal_test.step` at `cc=5, t=1` is refused at `export truth`. G23
classified solid 0's 26 readings and the degenerate-triangle skip (#36) plus the
fused cap (#34) took them to **13**, every survivor used once: 9 `M` and 4
`interior`, on faces `BRepCheck_Analyzer` calls valid with 1e-07 tolerances.
docs/specification.md 10 closes there, noting that
`spiral-island-unwritable.brep` shows the same two classes, so no narrowing of
the gate **by defect kind** is defensible.

This gate answers what those 13 are, and supplies a discriminator that is not
about defect kind at all.

**`M` is not two faces discretizing a shared edge differently.** Dumping both
faces' `PolygonOnTriangulation` for every implicated edge gives bit-identical
node coordinates, interning to identical keys. Nothing disagrees. What actually
happens is that the triangulation **does not cover the face** -- exact trimmed
area against summed triangle area, at the gate's own 0.05 mm deflection:

    face 317973  Plane     1.601941e+00 mm^2   covers 51.27 %
    face  77748  Cylinder  4.226842e-07 mm^2   covers 68.75 %
    face 500747  Cylinder  8.372943e-07 mm^2   covers 87.50 %

One per residual cluster, all three valid. The boundary segments in the
un-meshed part are what read as used-once. It is purely a function of
deflection, identical meshed alone or in context to five decimals, and each
face recovers once the deflection drops below its own feature size.

**The discriminator: refine, and a sound body converges while a broken one
diverges.** Across the whole residual the count falls 13 -> 8 -> 4 -> 4 -> **0**
at 0.05, 0.01, 0.002, 0.0005, 0.0001. `spiral-island-unwritable.brep` over the
same sweep goes **10 -> 13 -> 17 -> 21 -> 39**. A body whose *description* is
broken disagrees with itself more the harder you look; a body merely finer than
the ruler agrees as soon as the ruler is fine enough. Recorded rather than
applied -- it changes what the gate meshes at, and it wants a third part before
it becomes a rule.

**One attribution corrected on the way.** G23 noticed every folding
`ConicalSurface` beginning at `v = -sqrt(3)` and read it as a property of the
generated patches. It is the *input file's* own parametrization: 54
apex-touching conical faces, `RefRadius = 1.5`, `SemiAngle = pi/3`, radius
exactly zero there. Nothing about the lattice put an apex there. That is also
why the counter, before #36, read the accepted input body itself as carrying 86
non-manifold edges.

**Two cheap explanations ruled out**, recorded so they are not retried. The gate
meshes with `isInParallel=True` and an `M` is exactly the fault a per-face
parallel mesher could produce -- measured on solid 0 both ways, 25 bad edges,
identical `{4x: 12, 1x: 13}`. And face 317973 has a vertex 4.013e-04 mm from a
non-adjacent edge of its own wire whose tolerance is 3.346e-04 mm, the same
order and a convincing story for a face the mesher abandons half of -- but it is
*outside* that tolerance by 1.2x, and the deflection sweep is what explains it.

Usage::

    python tools/prototypes/g24_residual_coverage.py <shape.brep|shape.step> ...
    python tools/prototypes/g24_residual_coverage.py --deflection <d> <shape> ...
    python tools/prototypes/g24_residual_coverage.py --sweep <shape> ...
    python tools/prototypes/g24_residual_coverage.py --coverage <shape> [face,face]
    python tools/prototypes/g24_residual_coverage.py --extract <big.brep> <out.brep> <faces>

``--extract`` is the inner loop docs/specification.md 10 points at: reading and
meshing the rehearsal's 486 MB dominant body costs ~3 minutes an iteration,
while the residual clusters plus their edge-neighbours are 68 faces that
reproduce every remaining reading in 0.0 s.
"""

import math
import os
import sys
import time
from collections import defaultdict

import _bootstrap  # noqa: F401

from OCP.BRep import BRep_Tool, BRep_Builder
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepGProp import BRepGProp
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.GProp import GProp_GProps
from OCP.TopAbs import TopAbs_ShapeEnum
from OCP.TopLoc import TopLoc_Location
from OCP.TopoDS import TopoDS, TopoDS_Compound
from OCP.TopTools import TopTools_IndexedMapOfShape

from latticegen2 import occ
from latticegen2.parallel import read_brep, write_brep

QUANTUM = 1e7
"""1e-7 mm, OCCT's own confusion - the same grid `occ.exported_mesh_defects`
interns on, so a key here means what it means in production."""

SWEEP = (0.05, 0.01, 0.002, 0.0005, 0.0001)


def _key(a, b):
    lo, hi = (a, b) if a <= b else (b, a)
    return lo * 4_294_967_296 + hi


def _load(path):
    return read_brep(path) if path.endswith(".brep") else occ.read_step(path)


def analyse(shape, deflection=occ.DEFAULT_MESH_DEFLECTION):
    """Count mesh edges by use, skipping degenerate triangles as production does.

    Mirrors `occ.exported_mesh_defects` deliberately: a prototype that counts a
    different quantity from the gate cannot be used to reason about the gate.

    Returns ``(triangles, counts, point_of, degenerate)``.
    """
    BRepMesh_IncrementalMesh(shape, deflection, False, 0.2, True)

    point_id = {}
    point_of = {}
    counts = defaultdict(int)
    triangles = 0
    degenerate = 0

    for f in occ._explore(shape, TopAbs_ShapeEnum.TopAbs_FACE):
        face = TopoDS.Face_s(f)
        loc = TopLoc_Location()
        tri = BRep_Tool.Triangulation_s(face, loc)
        if tri is None:
            continue
        trsf = loc.Transformation()
        ids = []
        for i in range(1, tri.NbNodes() + 1):
            p = tri.Node(i).Transformed(trsf)
            key = (round(p.X() * QUANTUM), round(p.Y() * QUANTUM),
                   round(p.Z() * QUANTUM))
            got = point_id.get(key)
            if got is None:
                got = point_id[key] = len(point_id)
                point_of[got] = (p.X(), p.Y(), p.Z())
            ids.append(got)

        for i in range(1, tri.NbTriangles() + 1):
            a, b, c = tri.Triangle(i).Get()
            triangles += 1
            ia, ib, ic = ids[a - 1], ids[b - 1], ids[c - 1]
            if ia == ib or ib == ic or ia == ic:
                degenerate += 1
                continue
            for lo, hi in ((ia, ib), (ib, ic), (ic, ia)):
                counts[_key(lo, hi)] += 1

    return triangles, counts, point_of, degenerate


def report(path, deflection=occ.DEFAULT_MESH_DEFLECTION, show=14):
    t0 = time.time()
    tri, counts, point_of, degen = analyse(_load(path), deflection)
    bad = [(k, n) for k, n in counts.items() if n != 2]
    by_use = defaultdict(int)
    for _k, n in bad:
        by_use[f"{n}x"] += 1
    print(f"\n=== {os.path.basename(path)}: {tri} triangle(s), {len(bad)} bad "
          f"edge(s) {dict(by_use)}, {degen} degenerate skipped "
          f"[{time.time() - t0:.1f}s]")
    for k, n in sorted(bad, key=lambda x: -x[1])[:show]:
        lo, hi = divmod(k, 4_294_967_296)
        p, q = point_of[lo], point_of[hi]
        mid = tuple((p[i] + q[i]) / 2 for i in range(3))
        print(f"    used {n}x  mid [{mid[0]:.3f}, {mid[1]:.3f}, {mid[2]:.3f}]  "
              f"len {math.dist(p, q):.3e}")
    return len(bad)


def sweep(path):
    """Convergence or divergence under refinement - the discriminator.

    **Run this on a whole solid, never on an ``--extract``.** A cut-out patch is
    open, so most of its readings are its own cut boundary -- 150 of 166 on the
    rehearsal's extract, and they climb with refinement simply because a finer
    mesh puts more segments along the same open edge. On an extract, classify
    with ``g23_bad_edge_provenance.py`` and read the ``M``/``interior``/``T``
    tally, which ignores ``one-face``.
    """
    print(f"\n=== sweep {os.path.basename(path)}")
    for d in SWEEP:
        tri, counts, _pt, degen = analyse(_load(path), d)
        bad = sum(1 for n in counts.values() if n != 2)
        print(f"  deflection {d:<8}: {tri:>9} triangles, {bad:>4} bad, "
              f"{degen:>4} degenerate skipped", flush=True)


def coverage(path, only=None, deflection=occ.DEFAULT_MESH_DEFLECTION):
    """Exact trimmed area against summed triangle area, per face.

    This is what showed the residual readings to be an *incomplete*
    triangulation rather than an inconsistent one. Without ``only``, reports
    just the faces the mesher does not cover.
    """
    shape = _load(path)
    BRepMesh_IncrementalMesh(shape, deflection, False, 0.2, True)
    print(f"{'face':>7} {'exact mm^2':>13} {'ratio':>9} {'tri':>5}  surface")
    for fi, f in enumerate(occ._explore(shape, TopAbs_ShapeEnum.TopAbs_FACE)):
        if only is not None and fi not in only:
            continue
        face = TopoDS.Face_s(f)
        props = GProp_GProps()
        BRepGProp.SurfaceProperties_s(face, props)
        exact = props.Mass()
        loc = TopLoc_Location()
        tri = BRep_Tool.Triangulation_s(face, loc)
        if tri is None or not exact:
            continue
        trsf = loc.Transformation()
        total = 0.0
        for i in range(1, tri.NbTriangles() + 1):
            a, b, c = tri.Triangle(i).Get()
            pa = tri.Node(a).Transformed(trsf)
            pb = tri.Node(b).Transformed(trsf)
            pc = tri.Node(c).Transformed(trsf)
            u = (pb.X() - pa.X(), pb.Y() - pa.Y(), pb.Z() - pa.Z())
            v = (pc.X() - pa.X(), pc.Y() - pa.Y(), pc.Z() - pa.Z())
            cr = (u[1] * v[2] - u[2] * v[1],
                  u[2] * v[0] - u[0] * v[2],
                  u[0] * v[1] - u[1] * v[0])
            total += 0.5 * sum(x * x for x in cr) ** 0.5
        ratio = total / exact
        if only is None and ratio >= 0.999:
            continue
        kind = str(BRepAdaptor_Surface(face).GetType()).split("GeomAbs_")[-1]
        print(f"{fi:>7} {exact:13.6e} {ratio:9.5f} {tri.NbTriangles():>5}  {kind}")


def extract(src, out, want):
    """Cut ``want`` and everything sharing an edge with it into one compound.

    Adjacency is built with an edge map rather than
    `TopExp::MapShapesAndAncestors`, whose `TopTools_ListOfShape` OCP exposes
    with no iterator.
    """
    shape = read_brep(src)
    faces = [TopoDS.Face_s(f)
             for f in occ._explore(shape, TopAbs_ShapeEnum.TopAbs_FACE)]
    print(f"{len(faces)} faces", flush=True)

    edge_index = TopTools_IndexedMapOfShape()
    edges_of = [[edge_index.Add(e)
                 for e in occ._explore(face, TopAbs_ShapeEnum.TopAbs_EDGE)]
                for face in faces]
    by_edge = defaultdict(list)
    for fi, ids in enumerate(edges_of):
        for ei in ids:
            by_edge[ei].append(fi)

    keep = set(want)
    for fi in want:
        for ei in edges_of[fi]:
            keep.update(by_edge[ei])

    builder = BRep_Builder()
    comp = TopoDS_Compound()
    builder.MakeCompound(comp)
    for j in sorted(keep):
        builder.Add(comp, faces[j])
    write_brep(comp, out)
    print(f"kept {len(keep)} face(s) -> {out}", flush=True)


def main(argv):
    if not argv:
        print(__doc__)
        return 2
    if argv[0] == "--sweep":
        if len(argv) < 2:
            print("--sweep needs at least one shape")
            return 2
        for path in argv[1:]:
            sweep(path)
        return 0
    if argv[0] == "--extract":
        if len(argv) < 4:
            print("--extract needs <big.brep> <out.brep> <face,face,...>")
            return 2
        extract(argv[1], argv[2], [int(x) for x in argv[3].split(",")])
        return 0
    if argv[0] == "--coverage":
        if len(argv) < 2:
            print("--coverage needs a shape")
            return 2
        only = {int(x) for x in argv[2].split(",")} if len(argv) > 2 else None
        coverage(argv[1], only)
        return 0
    defl = occ.DEFAULT_MESH_DEFLECTION
    if argv[0] == "--deflection":
        if len(argv) < 3:
            print("--deflection needs a value and at least one shape")
            return 2
        defl = float(argv[1])
        argv = argv[2:]
    for path in argv:
        report(path, deflection=defl)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
