"""G26: does round 2 of the boundary sew fail only because it crosses a process?

Controlled measurement for docs/algorithm.md §8. The run is a real pipeline
run with the tiling threshold lowered, so a committed part tiles. At round 2
the real round-1 tile results of every tiled component are sewn three ways,
from the same inputs:

  serial   the seam-only subset sewn in this process, interior faces carried
  pooled   the same subset sewn in a pool worker and read back, as the pooled
           branch of `_sew_round_two` does
  unsplit  the whole tile results sewn in this process (the repair)

Each route's free-edge count is compared with the count a correct sew leaves,
4 x the component's interior interfaces, and timed. Then the run goes on
through the real `_sew_round_two`, so its own repair count is reported too.

    python tools/prototypes/g26_round2_pool_identity.py <input> <cc> <t> \
        [--tile-target 100] [--min-to-tile 300] [--cores 6]
"""

import argparse
import os
import sys
import tempfile
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("cc")
    ap.add_argument("t")
    ap.add_argument("--tile-target", type=int, default=100)
    ap.add_argument("--min-to-tile", type=int, default=300)
    ap.add_argument("--cores", default="6")
    a = ap.parse_args()

    from latticegen2 import occ, weld
    from latticegen2.cli import parse_args
    from latticegen2.pipeline import run_pipeline
    from latticegen2.runlog import RunLog, stages_for

    real_sew_boundary = weld.sew_boundary
    real_round_two = weld._sew_round_two
    results = []

    def sew_boundary(*args, **kw):
        kw["tile_target"], kw["min_to_tile"] = a.tile_target, a.min_to_tile
        got = real_sew_boundary(*args, **kw)
        results.append(("run", got[1].repaired_components, got[1].tiled_components,
                        got[1].tiles, got[1].t_round2, got[1].t_repair))
        return got

    def round_two(by_group, plan, tile_results, tolerance, workers, tmpdir, pool,
                  expected_rings=None, **kw):
        print(f"G26: {len(by_group)} component(s), "
              f"{sum(1 for p in plan.values() if p is not None)} tiled", flush=True)
        for group, tiles in plan.items():
            if tiles is None:
                continue
            want = 4 * (expected_rings or {}).get(group, 0)
            seam, interior = weld._split_seam_interior(tile_results[group])
            n_seam = sum(len(s) for s in seam)
            row = {"group": group, "tiles": len(tiles), "want": want,
                   "seam_faces": n_seam, "interior_faces": len(interior)}

            t0 = time.perf_counter()
            sewn = weld._sew_faces(seam, tolerance)
            row["serial_s"] = time.perf_counter() - t0
            row["serial_free"] = len(weld.free_edges(sewn + interior))

            if pool is not None and pool.active:
                in_path = os.path.join(tmpdir, f"g26_{group}.brep")
                out_path = os.path.join(tmpdir, f"g26_{group}_out.brep")
                t0 = time.perf_counter()
                weld._write_brep(occ.compound(occ.compound(fl) for fl in seam), in_path)
                raw, _rss = pool.run(weld._worker_sew_tile, [(in_path, out_path, tolerance)])
                back = occ.faces(weld._read_brep(raw[0][0]))
                row["pooled_s"] = time.perf_counter() - t0
                row["pooled_free"] = len(weld.free_edges(back + interior))

            t0 = time.perf_counter()
            full = weld._sew_faces(tile_results[group], tolerance)
            row["unsplit_s"] = time.perf_counter() - t0
            row["unsplit_free"] = len(weld.free_edges(full))
            results.append(("route", row))
            print(f"G26: {row}", flush=True)
        return real_round_two(by_group, plan, tile_results, tolerance, workers, tmpdir,
                              pool, expected_rings=expected_rings, **kw)

    weld.sew_boundary = sew_boundary
    weld._sew_round_two = round_two

    out_dir = tempfile.mkdtemp(prefix="g26-")
    argv = ["-i", a.input, "-cc", a.cc, "-t", a.t, "--cores", a.cores,
            "-o", os.path.join(out_dir, "g26.step")]
    args = parse_args(argv)
    rl = RunLog(args.log_path, verbose=False).open()
    rl.stage_names = stages_for(False)
    try:
        run_pipeline(args, rl)
    finally:
        rl.close()
    for r in results:
        print("G26 result:", r)
    print("G26 log:", args.log_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
