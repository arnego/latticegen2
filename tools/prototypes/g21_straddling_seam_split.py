"""Gate G21 — can the seam-only round 2 be made correct on real trimmed parts?

docs/specification.md §10's one remaining live item. Boundary-sew round 2 sews
only the free-edge-bearing subset of each tile's round-1 result
(`weld._split_seam_interior`, docs/algorithm.md §8), which G8 proved identical
to a full round 2 at every prototype scale tried — and G9 found it wrong on the
`cc=5, t=1` rehearsal's real, heavily trimmed pieces, where it left 118,760
open edges. The `expected_rings` check catches that and redoes the component
with a full unsplit sew, which on that part costs **651.2 s of `stitch`'s
769 s**. §10's proposed fix, recorded there as unproven:

    carry a seam face's straddling neighbours into the sewn subset, so
    BRepBuilderAPI_Sewing cannot rebuild a shared edge onto a new TopoDS_Edge
    while the carried face keeps the original.

**The doubt this gate exists to settle.** A "straddling" edge is one used twice
inside the tile but only once inside the seam subset — i.e. an edge of the cut
between the subset and the carried remainder. A subset with no such edge is
exactly a union of connected components of the tile's face-adjacency graph, and
a tile's round-1 result is normally one connected shell. If that is the whole
story, the fixpoint of "add my straddling neighbours" is the entire tile — the
unsplit sew the fallback already runs, at no saving — and one hop merely moves
the frontier outward. Whether it does depends on what actually triggers the
rebuild inside `BRepBuilderAPI_Sewing`, which is a question about OCCT, not
about this codebase. So it is measured.

Four parts, in the order that decides it:

* **A — reproduce.** Round 2 computed both ways on the same round-1 results.
  A gate that cannot reproduce the divergence proves nothing about a fix for
  it, and G9 records that no synthetic setup has reproduced it — which is why
  the geometry here is the real part, trimmed for real, and only the *extent*
  is reduced.
* **B — mechanism.** Which edges come back as new objects, and why. The
  decisive probe is the interface holes: those are free in the input and have
  no partner anywhere, so if they are rebuilt too, the trigger is "free in my
  input" rather than "merged with something", and every hop of any closure
  creates a fresh frontier of exactly that kind.
* **C — closure size.** |S0| (today's subset), |S1| (one hop), and the
  fixpoint, as a fraction of tile faces, with the hop count.
* **D — cost.** Full round 2 against seam-only and against each closure, so a
  correct variant can still be rejected for costing what it saves.

Pass bar, set in advance: some k-hop closure must reproduce the full round-2
result **exactly** — same face count, same free-edge count, same total area to
machine precision — **and** cost materially less than the full unsplit sew.
Both, or this is a disproof and the fallback's cost is structural.

    python tools/prototypes/g21_straddling_seam_split.py --cores 6
    python tools/prototypes/g21_straddling_seam_split.py --tiles <dir>

Building the geometry is the slow part (a real classify and a real trim of a
block of the rehearsal part), so round 1's results are cached as one `.brep`
per tile — which is also exactly the form production round 1 hands round 2,
since its tiles come back from worker processes through `.brep` files. A
directory of such files, from the cache or lifted out of a failed run's temp
folder (`sew_tile_<group>_<i>_out.brep`), can be replayed with ``--tiles``.
"""

import argparse
import os
import shutil
import time

import _bootstrap  # noqa: F401

from OCP.TopAbs import TopAbs_ShapeEnum
from OCP.TopExp import TopExp
from OCP.TopTools import (
    TopTools_IndexedDataMapOfShapeListOfShape,
    TopTools_IndexedMapOfShape,
)

from latticegen2 import occ, weld
from latticegen2.boundary import resolve_interfaces, trim_boundary
from latticegen2.classify import Klass, classify_nodes, tessellate_surface
from latticegen2.junction import build_template
from latticegen2.lattice import candidate_nodes, lattice_params
from latticegen2.parallel import WorkerPool, read_brep, write_brep

ROOT = _bootstrap.ROOT
DEFAULT_INPUT = os.path.join(ROOT, "test", "TD_HX_rehearsal_test.step")


# --- the geometry -----------------------------------------------------------


def _sub_box(lo, hi, frac):
    """A corner block of the part's bounding box, ``frac`` of it on each axis.

    A corner rather than a centred block: the part is a heat exchanger whose
    interior is mostly empty, and a corner reliably holds real surface.
    """
    return lo, lo + (hi - lo) * frac


def build_tiles(args, cache: str | None) -> list:
    """Round 1's per-tile results for a block of the real part."""
    tmpdir = os.path.join(_bootstrap.ROOT, "temp", "g21")
    os.makedirs(tmpdir, exist_ok=True)
    print(f"  reading {os.path.basename(args.input)} at cc={args.cc}, t={args.t}")
    body = occ.read_step(args.input)
    lp = lattice_params(args.cc, args.t)
    tpl = build_template(lp)
    mesh = tessellate_surface(body, lp)
    print(f"  surface mesh: {len(mesh.tris)} triangles, d={mesh.deviation:.5f} mm")

    lo, hi = occ.bounding_box(body)
    blo, bhi = _sub_box(lo, hi, args.block)
    cands = candidate_nodes(lp, blo, bhi)
    print(f"  block {args.block:.2f} of the bounding box: {len(cands)} candidate nodes")

    body_path = os.path.join(tmpdir, "input.brep")
    write_brep(body, body_path)
    with WorkerPool(args.cores) as pool:
        t0 = time.perf_counter()
        cls = classify_nodes(lp, mesh, cands, tmpdir=tmpdir, pool=pool)
        counts = cls.counts()
        print(f"  classify {time.perf_counter() - t0:.1f}s: {counts}")
        boundary_nodes = cls.of(Klass.BOUNDARY)
        t0 = time.perf_counter()
        trimmed = trim_boundary(
            lp, tpl, boundary_nodes, body, body_path, tmpdir,
            workers=args.cores, pool=pool,
        )
        print(f"  trim {time.perf_counter() - t0:.1f}s: {len(trimmed.pieces)} pieces "
              f"from {len(boundary_nodes)} junctions ({trimmed.n_empty} empty)")

    resolve_interfaces(lp, cls.of(Klass.INTERIOR), trimmed.pieces)
    plan = weld._tile_pieces(trimmed.pieces, args.tile_target, min_to_tile=1)
    if plan is None:
        raise SystemExit(
            "this block did not tile - raise --block, or lower --tile-target"
        )
    print(f"  tiled into {len(plan)} tile(s) of "
          f"{min(len(t) for t in plan)}-{max(len(t) for t in plan)} pieces")

    t0 = time.perf_counter()
    results, _rss = weld._sew_all_tiles(
        {0: plan}, weld.SEW_TOLERANCE, workers=args.cores, tmpdir=tmpdir, pool=None
    )
    print(f"  round 1 {time.perf_counter() - t0:.1f}s")
    tiles = results[0]
    if cache:
        os.makedirs(cache, exist_ok=True)
        for i, faces in enumerate(tiles):
            write_brep(occ.compound(faces), os.path.join(cache, f"sew_tile_0_{i}_out.brep"))
        print(f"  cached {len(tiles)} tile(s) to {cache}")
    shutil.rmtree(tmpdir, ignore_errors=True)
    return tiles


def load_tiles(path: str) -> list:
    """Replay round 1 from `.brep` files — production's own transport.

    Each file is read separately, exactly as `_sew_all_tiles` reads its
    workers' output, so sharing is preserved within a tile and absent between
    tiles. That is the topology round 2 really sees (G15).
    """
    names = sorted(n for n in os.listdir(path) if n.endswith(".brep"))
    tiles = [occ.faces(read_brep(os.path.join(path, n))) for n in names]
    print(f"  {len(tiles)} tile(s) from {path}: "
          f"{sum(len(t) for t in tiles)} faces total")
    return tiles


# --- the split, parameterised by how far it closes --------------------------


def _edge_faces(faces):
    m = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(
        occ.faces_shell(faces), TopAbs_ShapeEnum.TopAbs_EDGE, TopAbs_ShapeEnum.TopAbs_FACE, m
    )
    return m


def _seam_faces(edge_faces) -> TopTools_IndexedMapOfShape:
    """Faces bearing a free edge — `_split_seam_interior`'s own S0."""
    seam = TopTools_IndexedMapOfShape()
    for i in range(1, edge_faces.Extent() + 1):
        owners = edge_faces.FindFromIndex(i)
        if owners.Extent() == 1:
            seam.Add(owners.First())
    return seam


def _straddling(edge_faces, seam) -> list:
    """Edges used twice in the tile with exactly one owner in ``seam``."""
    out = []
    for i in range(1, edge_faces.Extent() + 1):
        owners = edge_faces.FindFromIndex(i)
        if owners.Extent() != 2:
            continue
        inside = sum(1 for f in (owners.First(), owners.Last()) if seam.Contains(f))
        if inside == 1:
            out.append((edge_faces.FindKey(i), owners.First(), owners.Last()))
    return out


def _grow(edge_faces, seam, hops: int):
    """Add the other owner of every straddling edge, ``hops`` times.

    Returns ``(faces, hops_taken)``. ``hops=0`` is today's behaviour; a
    negative value means "to the fixpoint", and ``hops_taken`` then says how
    many rounds that took.
    """
    grown = TopTools_IndexedMapOfShape()
    for i in range(1, seam.Extent() + 1):
        grown.Add(seam.FindKey(i))
    n = 0
    while hops < 0 or n < hops:
        fresh = [
            b if grown.Contains(a) else a
            for _e, a, b in _straddling(edge_faces, grown)
        ]
        fresh = [f for f in fresh if not grown.Contains(f)]
        if not fresh:
            break
        for f in fresh:
            grown.Add(f)
        n += 1
    return grown, n


def split(face_lists, hops: int):
    """`_split_seam_interior`, with the seam set closed ``hops`` times."""
    seam_lists, interior = [], []
    for faces in face_lists:
        edge_faces = _edge_faces(faces)
        keep, _n = _grow(edge_faces, _seam_faces(edge_faces), hops)
        if keep.Extent() == 0:
            interior.extend(faces)
            continue
        seam, rest = [], []
        for f in faces:
            (seam if keep.Contains(f) else rest).append(f)
        seam_lists.append(seam)
        interior.extend(rest)
    return seam_lists, interior


def measure(faces) -> tuple:
    return len(faces), len(weld.free_edges(faces)), sum(occ.area(f) for f in faces)


def timed(fn):
    t0 = time.perf_counter()
    out = fn()
    return out, time.perf_counter() - t0


# --- the parts --------------------------------------------------------------


def part_a(tiles):
    print()
    print("=" * 72)
    print("A - reproduce: does the seam-only split still agree with a full round 2?")
    print("=" * 72)
    full, t_full = timed(lambda: weld._sew_faces(tiles, weld.SEW_TOLERANCE))
    seam_lists, interior = split(tiles, hops=0)
    seam, t_seam = timed(lambda: weld._sew_faces(seam_lists, weld.SEW_TOLERANCE))
    seam = seam + interior

    print(f"  {'round 2':<28} {'faces':>8} {'free edges':>11} {'area mm^2':>16} {'time':>8}")
    for label, faces, dt in (("full, unsplit", full, t_full), ("seam-only (today)", seam, t_seam)):
        n, fe, area = measure(faces)
        print(f"  {label:<28} {n:>8} {fe:>11} {area:>16.6f} {dt:>7.2f}s")
    agree = measure(full)[:2] == measure(seam)[:2]
    print(f"\n  seam-only {'agrees with' if agree else 'DIVERGES from'} the full sew"
          f" on this input.")
    if agree:
        print(
            "  The divergence G9 measured is not reproduced on this face set. Parts C\n"
            "  and D below are unaffected: what a straddle-free subset costs is a\n"
            "  property of the tiles, not of whether this particular sew diverged."
        )
    return full, (t_full, t_seam), agree


def part_b(tiles, full):
    print()
    print("=" * 72)
    print("B - mechanism: which edges come back as new objects")
    print("=" * 72)
    # Over every tile at once, because that is what round 2 sews: a single
    # tile's seam subset has no partner to merge with, so sewing it alone
    # rebuilds nothing and would make any mechanism look absent.
    strad_edges, both_in, free_in = [], [], []
    seam_lists = []
    for tile in tiles:
        edge_faces = _edge_faces(tile)
        seam = _seam_faces(edge_faces)
        seam_lists.append([f for f in tile if seam.Contains(f)])
        strad_map = TopTools_IndexedMapOfShape()
        for e, _a, _b in _straddling(edge_faces, seam):
            strad_map.Add(e)
            strad_edges.append(e)
        # Only edges the subset actually carries can be judged: an edge whose
        # two owners were both carried through is absent from the output
        # because it was never handed to the sew, which is not the same thing
        # as rebuilt.
        for i in range(1, edge_faces.Extent() + 1):
            edge = edge_faces.FindKey(i)
            owners = edge_faces.FindFromIndex(i)
            if owners.Extent() == 2:
                if (not strad_map.Contains(edge) and seam.Contains(owners.First())
                        and seam.Contains(owners.Last())):
                    both_in.append(edge)
            elif owners.Extent() == 1 and seam.Contains(owners.First()):
                free_in.append(edge)
    print(f"  {len(tiles)} tile(s), {sum(len(t) for t in tiles)} faces: "
          f"{sum(len(s) for s in seam_lists)} in the seam subset, "
          f"{len(strad_edges)} straddling edge(s)")

    sewn = weld._sew_faces(seam_lists, weld.SEW_TOLERANCE)
    out_edges = TopTools_IndexedMapOfShape()
    for f in sewn:
        for e in occ._explore(f, TopAbs_ShapeEnum.TopAbs_EDGE):
            out_edges.Add(e)

    def survived(edges):
        kept = sum(1 for e in edges if out_edges.Contains(e))
        return kept, len(edges) - kept
    print(f"\n  {'edge class (only edges the subset carries)':<44} {'survived':>9} {'rebuilt':>8}")
    for label, edges in (
        ("free in the tile (holes and seams)", free_in),
        ("straddling (2 owners, 1 in subset)", strad_edges),
        ("used twice, both owners in the subset", both_in),
    ):
        kept, gone = survived(edges)
        print(f"  {label:<44} {kept:>9} {gone:>8}")
    print(
        "\n  A free edge with no partner anywhere - an interface hole - is the decisive\n"
        "  row: if those are rebuilt too, the trigger is 'free in my input', and every\n"
        "  hop of any closure hands sewing a fresh frontier of exactly that kind."
    )


def part_c(tiles):
    print()
    print("=" * 72)
    print("C - closure size: how far 'add my straddling neighbours' spreads")
    print("=" * 72)
    print(f"  {'tile':>6} {'faces':>8} {'|S0|':>8} {'|S1|':>8} {'|S2|':>8} "
          f"{'fixpoint':>9} {'hops':>5}")
    for i, faces in enumerate(sorted(tiles, key=len, reverse=True)[:5]):
        edge_faces = _edge_faces(faces)
        s0 = _seam_faces(edge_faces)
        sizes = [_grow(edge_faces, s0, h)[0].Extent() for h in (0, 1, 2)]
        fix, hops = _grow(edge_faces, s0, -1)
        print(f"  {i:>6} {len(faces):>8} {sizes[0]:>8} {sizes[1]:>8} {sizes[2]:>8} "
              f"{fix.Extent():>9} {hops:>5}")
    print(
        "\n  A fixpoint equal to the tile's face count means the only straddle-free\n"
        "  subset is the whole tile - which is the unsplit sew the fallback already\n"
        "  runs, at no saving."
    )


def part_d(tiles, full, times):
    print()
    print("=" * 72)
    print("D - cost, and whether any correct variant is cheaper")
    print("=" * 72)
    t_full, t_seam = times
    want = measure(full)
    print(f"  {'variant':<28} {'faces':>8} {'free edges':>11} {'area mm^2':>16} "
          f"{'time':>8} {'exact':>6}")
    n, fe, area = want
    print(f"  {'full, unsplit':<28} {n:>8} {fe:>11} {area:>16.6f} {t_full:>7.2f}s "
          f"{'-':>6}")
    ok = False
    for hops in (0, 1, 2, -1):
        seam_lists, interior = split(tiles, hops)
        (sewn, dt) = timed(lambda: weld._sew_faces(seam_lists, weld.SEW_TOLERANCE))
        got = measure(sewn + interior)
        exact = got[:2] == want[:2] and abs(got[2] - want[2]) <= 1e-9 * max(1.0, want[2])
        label = {0: "seam-only (today)", -1: "closure to fixpoint"}.get(hops, f"closure, {hops} hop(s)")
        print(f"  {label:<28} {got[0]:>8} {got[1]:>11} {got[2]:>16.6f} {dt:>7.2f}s "
              f"{str(exact):>6}")
        # Only the fixpoint counts toward the verdict. The intermediate rows are
        # there to show the cost curve, not as candidates: a partial closure
        # still leaves straddling edges, so its coming out exact on one input
        # says nothing about the next one - which is precisely the trap G8 fell
        # into. Today's seam-only split is the hops=0 row and is the thing under
        # test, not a way to pass it.
        if hops == -1:
            ok = exact and dt < 0.9 * t_full
    print(
        "\n  Only the fixpoint row can pass: it is the one variant with no straddling\n"
        "  edge left, so it is the only one whose exactness is structural rather than\n"
        "  incidental. Exact-but-not-cheaper is the fallback with extra steps;\n"
        "  cheaper-but-not-exact is what G9 already measured at 118,760 open edges."
    )
    return ok


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-i", "--input", default=DEFAULT_INPUT)
    ap.add_argument("-cc", type=float, default=5.0)
    ap.add_argument("-t", type=float, default=1.0)
    ap.add_argument("--block", type=float, default=0.22,
                    help="fraction of the bounding box, per axis, to build")
    ap.add_argument("--tile-target", type=int, default=weld.TILE_TARGET_PIECES)
    ap.add_argument("--cores", type=int, default=6)
    ap.add_argument("--cache", default=os.path.join(ROOT, "temp", "g21-tiles"))
    ap.add_argument("--tiles", default=None, help="replay round 1 from this directory")
    args = ap.parse_args()
    occ.quiet_kernel()

    print("Round 1 input")
    print("-" * 72)
    if args.tiles:
        tiles = load_tiles(args.tiles)
    elif args.cache and os.path.isdir(args.cache) and os.listdir(args.cache):
        print(f"  reusing the cache at {args.cache} (delete it to rebuild)")
        tiles = load_tiles(args.cache)
    else:
        tiles = build_tiles(args, args.cache)

    full, times, agree = part_a(tiles)
    part_b(tiles, full)
    part_c(tiles)
    cheaper = part_d(tiles, full, times)

    print()
    print("=" * 72)
    if cheaper:
        print("PASS - a straddle-free subset reproduces the full sew and costs less")
    else:
        print("FAIL - a straddle-free subset is the whole tile, at the unsplit sew's cost")
    if agree:
        print(
            "NOTE - this input does not exhibit the divergence itself, which does not\n"
            "       touch the verdict above: that is about what a *correct* subset\n"
            "       costs. Do not read it as the split being sound either. The\n"
            "       pipeline sews the face set `boundary.finalize_pieces` leaves, and\n"
            "       a gate stopping short of that call sews a different one;\n"
            "       production's own evidence line is the reference\n"
            "       (docs/specification.md §11)."
        )
    print("=" * 72)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
