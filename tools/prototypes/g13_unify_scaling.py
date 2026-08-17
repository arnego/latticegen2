"""Gate G13 — can same-domain unification be tiled *below* the solid?

docs/specification.md §10 records optimization path 1 (parallelising `simplify`
across the shared `WorkerPool`) as implemented, correct, and **not a wall-clock
win** on the `cc=5, t=1` rehearsal: 17 m 17 s -> 18 m 39 s, still 0.99 cores,
because that part's 14 solids are one dominant body plus 13 scraps and "the
largest single solid is the floor, not the sum". Dispatching body-for-body
cannot lower that floor. Tiling *within* a solid is the only lever that can.

Same-domain unification is a **local** operation — a merge group is a connected
set of coplanar adjacent faces — so partitioning a solid's faces and unifying
each subset can only ever *miss* merges that straddle a tile seam, never merge
something it should not. That failure mode is "the output is slightly larger",
which is what docs/algorithm.md §11 requires of an optimization and what §9
already accepts when the kernel declines to unify at all.

Two questions decide whether it is worth building, and neither is answered
anywhere in this project's measurements today:

**A. Is the cost linear or superlinear in face count?** If linear, tiling is
worth at most `W` and the case rests entirely on parallelism. If superlinear,
tiling wins even serially and wins more than `W` in parallel. The two rates
that exist point in opposite directions and cannot both be right: G7 unified 8
solids of 702-8,526 faces in 5.138 s (~0.17 ms/face), while the rehearsal did
1,006,505 faces in 18 m 39 s (~1.11 ms/face) — 6.5x worse per face, which is
either superlinearity or just trimmed curved boundary faces costing more.
docs/algorithm.md §9's quoted "~0.24 ms/face" reconciles with neither. This
gate measures the exponent directly on one geometry family, so the comparison
is like-for-like.

**B. Do unified tiles reassemble without sewing?** This is the real risk. The
whole architecture depends on never putting a volume-scaling face set through
`BRepBuilderAPI_Sewing` (docs/algorithm.md §6, §8). Reassembling tiles is free
*iff* OCCT preserves the identity (`TShape` + location) of the edges and
vertices it did not touch, so the two sides of a tile seam remain literally the
same edge and `BRep_Builder.Add` alone yields a closed shell. If identity is
not preserved, every seam edge appears twice as two distinct objects, the shell
has free edges, and the idea collapses back into the cost it exists to avoid.

The test for B is binary and needs no tolerance: the input is an all-interior
grid, which is **closed** (zero free edges). Unify its faces in tiles,
concatenate the results, and count edges used by exactly one face. Zero means
identity held. Anything else means it did not.

`unify_edges` is measured both ways because it is the obvious mitigation for B:
edge concatenation is the one pass that rewrites edges *inside* a merged wire,
so switching it off is what would keep a seam edge untouched — and
docs/algorithm.md §9 already measured OCCT's edge pass as worth 4 edges out of
81,816 here, so it would cost almost nothing to give up.

    python tools/prototypes/g13_unify_scaling.py
    python tools/prototypes/g13_unify_scaling.py --max-m 12   # quicker

Bars set in advance:

* **A**: exponent >= 1.15 over the measured range => superlinear, tiling wins
  serially too. Below that, treat as linear and judge on parallelism alone.
* **B**: zero free edges after reassembly, and the reassembled solid passes
  `BRepCheck_Analyzer`, with a merge loss under 5 %. Anything else is a fail
  and the fallback (seam-only sew, or a merge-group partition that cannot
  straddle) has to be designed instead.
"""

import argparse
import time

import _bootstrap  # noqa: F401
import numpy as np

from latticegen2 import occ
from latticegen2.connect import lattice_interfaces
from latticegen2.interior import build_interior_shell, extract_template_mesh
from latticegen2.junction import build_template
from latticegen2.lattice import lattice_params
from latticegen2.weld import free_edges

CC, T = 10.0, 1.5
SIZES = [4, 6, 8, 10, 12, 14, 16]
SUPERLINEAR_BAR = 1.15
MERGE_LOSS_BAR = 0.05


def grid(m: int) -> np.ndarray:
    a = np.arange(m, dtype=np.int64)
    g = np.meshgrid(a, a, a, indexing="ij")
    return np.stack([x.ravel() for x in g], axis=1)


def build_solid(lp, tpl, tmesh, m: int):
    """One closed, all-interior m x m x m grid — the same family G7 used.

    Closed matters for part B: a shell with zero free edges by construction is
    what turns "did the seam survive" into a count rather than a judgement.
    """
    ns = grid(m)
    kept = {(int(a), int(b), int(c)) for a, b, c in ns}
    built = build_interior_shell(lp, tpl, tmesh, ns, lattice_interfaces(kept))
    shell = list(built.shells.values())[0]
    assert built.stats["interior_open_edges"] == 0, "grid should be closed"
    shell.Closed(True)
    return occ.make_solid(shell)


# --- part A: scaling ---------------------------------------------------------


def scaling(lp, tpl, tmesh, sizes, keep: int = 2):
    """Time a full-solid unification at each scale, keeping the largest solids.

    ``keep`` solids are retained for part B, so the merge loss tiling costs can
    be read against *scale* as well as against tile count — the seam fraction
    of a tiled partition falls as the part grows, so a loss measured on a small
    grid is an over-estimate of the production one and must not be quoted as if
    it were not.
    """
    print("--- A: unify_same_domain cost against face count ---")
    print(f"{'m':>4} {'junctions':>10} {'faces':>9} {'build s':>9} "
          f"{'unify s':>9} {'ms/face':>9} {'faces out':>10}")
    rows = []
    kept: dict[int, object] = {}
    for m in sizes:
        t0 = time.perf_counter()
        solid = build_solid(lp, tpl, tmesh, m)
        t_build = time.perf_counter() - t0

        n_faces, _ = occ.count_subshapes(solid)
        t0 = time.perf_counter()
        merged = occ.unify_same_domain(solid)
        t_unify = time.perf_counter() - t0
        faces_out, _ = occ.count_subshapes(merged)

        rows.append((m, n_faces, t_unify))
        print(f"{m:>4} {m**3:>10} {n_faces:>9} {t_build:>9.2f} "
              f"{t_unify:>9.3f} {1000 * t_unify / n_faces:>9.3f} {faces_out:>10}")
        kept[m] = solid
        for stale in sorted(kept)[:-keep]:
            del kept[stale]
    return rows, kept


def exponents(rows):
    """Local slopes plus one least-squares fit over the whole range.

    Local slopes are reported because a single fit hides a cost that is linear
    until some internal structure blows up and superlinear after it — which is
    exactly the shape that would make a 1 M-face solid behave unlike an 8 k one.
    """
    print("\nlog-log slope d(log t)/d(log faces):")
    for (m0, f0, t0), (m1, f1, t1) in zip(rows, rows[1:]):
        k = np.log(t1 / t0) / np.log(f1 / f0)
        print(f"  m {m0:>2} -> {m1:<2}  ({f0:>7} -> {f1:<7} faces)   {k:.3f}")
    lf = np.log([r[1] for r in rows])
    lt = np.log([r[2] for r in rows])
    k, _ = np.polyfit(lf, lt, 1)
    print(f"  overall least-squares fit:                     {k:.3f}")
    return float(k)


# --- part B: tiling ----------------------------------------------------------


def face_centroids(faces) -> np.ndarray:
    out = np.empty((len(faces), 3))
    for i, f in enumerate(faces):
        lo, hi = occ.bounding_box(f)
        out[i] = 0.5 * (lo + hi)
    return out


def tile_faces(faces, n_per_axis: int):
    """Bucket faces into an ``n^3`` spatial grid by centroid.

    Deliberately the *generic* partition rather than a lattice-aware one: it
    splits merge groups at every tile seam, so whatever merge loss it reports is
    an upper bound on what a strut-aware partition (which can be made not to
    straddle at all) would cost.
    """
    c = face_centroids(faces)
    lo, hi = c.min(axis=0), c.max(axis=0)
    span = np.maximum(hi - lo, 1e-12)
    idx = np.clip(((c - lo) / span * n_per_axis).astype(int), 0, n_per_axis - 1)
    buckets: dict[tuple, list] = {}
    for f, key in zip(faces, map(tuple, idx)):
        buckets.setdefault(key, []).append(f)
    return list(buckets.values())


def tiling_trial(all_faces, whole_faces_out, t_whole, n_per_axis, unify_edges):
    tiles = tile_faces(all_faces, n_per_axis)
    out_faces, times = [], []
    for tf in tiles:
        shell = occ.faces_shell(tf)
        t0 = time.perf_counter()
        merged = occ.unify_same_domain(shell, unify_edges=unify_edges)
        times.append(time.perf_counter() - t0)
        out_faces.extend(occ.faces(merged))

    fe = len(free_edges(out_faces))
    loss = (len(out_faces) - whole_faces_out) / max(whole_faces_out, 1)
    volume = valid = None
    if fe == 0:
        shell = occ.faces_shell(out_faces)
        shell.Closed(True)
        rebuilt = occ.make_solid(shell)
        volume = occ.volume(rebuilt)
        valid = occ.is_valid(rebuilt)

    print(f"  {len(tiles):>4} tiles  serial {sum(times):>8.3f}s  slowest "
          f"{max(times):>7.3f}s  ceiling {t_whole / max(times):>6.2f}x  "
          f"faces {len(out_faces):>7} ({100 * loss:+.2f}%)  free edges {fe:>6}"
          + (f"  valid {valid}" if valid is not None else ""))
    return fe, loss, volume, valid, max(times)


def tiling(solid, tile_counts, label):
    print(f"\n--- B[{label}]: does a tiled unification reassemble without sewing? ---")
    all_faces = occ.faces(solid)
    v0 = occ.volume(solid)
    print(f"input: {len(all_faces)} faces, volume {v0:.6f} mm^3, "
          f"{len(free_edges(all_faces))} free edges")

    results = {}
    for unify_edges in (True, False):
        t0 = time.perf_counter()
        whole = occ.unify_same_domain(solid, unify_edges=unify_edges)
        t_whole = time.perf_counter() - t0
        whole_faces = occ.faces(whole)
        f_out, e_out = occ.count_subshapes(whole)
        print(f"\nunify_edges={unify_edges}: whole solid {t_whole:.3f}s -> "
              f"{f_out} faces, {e_out} edges, {len(free_edges(whole_faces))} free, "
              f"volume drift {abs(occ.volume(whole) - v0) / v0:.2e}")
        for n in tile_counts:
            results[(unify_edges, n)] = tiling_trial(
                all_faces, len(whole_faces), t_whole, n, unify_edges
            )
    return results, v0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-m", type=int, default=max(SIZES),
                    help="largest grid edge in cells (cost grows as m^3)")
    ap.add_argument("--tiles", type=int, nargs="+", default=[2, 3],
                    help="tiles per axis to try in part B (n^3 tiles each)")
    args = ap.parse_args()
    sizes = [m for m in SIZES if m <= args.max_m]

    print(f"Building interior grids (cc={CC}, t={T}) ...\n")
    lp = lattice_params(CC, T)
    tpl = build_template(lp)
    tmesh = extract_template_mesh(lp, tpl)

    rows, kept = scaling(lp, tpl, tmesh, sizes)
    k = exponents(rows)
    trials = {m: tiling(solid, args.tiles, f"m={m}") for m, solid in sorted(kept.items())}

    print("\n=== VERDICT ===")
    superlinear = k >= SUPERLINEAR_BAR
    reading = ("SUPERLINEAR: tiling wins serially as well as in parallel"
               if superlinear else
               "effectively LINEAR: tiling is worth at most W, no more")
    print(f"A. exponent {k:.3f} -> {reading}")

    passed = False
    for m, (results, v0) in sorted(trials.items()):
        for (unify_edges, n), (fe, loss, volume, valid, _t) in sorted(results.items()):
            good = fe == 0 and valid and loss <= MERGE_LOSS_BAR
            passed = passed or good
            drift = abs(volume - v0) / v0 if volume is not None else float("nan")
            print(f"B. m={m:<3} unify_edges={str(unify_edges):<5} {n**3:>4} tiles: "
                  f"{'PASS' if good else 'FAIL'}  free edges {fe:>6}, merge loss "
                  f"{100 * loss:+.2f}%, volume drift {drift:.2e}")
    if passed:
        print("\nTiles reassemble by shared topology — no sewing needed on the "
              "volume-scaling face set.")
    else:
        print("\nNo configuration reassembles cleanly within the bars. Sub-body "
              "tiling needs a seam repair, and its cost has to be measured "
              "before anything is built on it.")


if __name__ == "__main__":
    main()
