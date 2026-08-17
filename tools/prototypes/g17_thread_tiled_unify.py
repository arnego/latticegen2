"""Gate G17 — can tiled unification be threaded instead of processed?

G15 killed docs/specification.md §11's Phase 3 on the *transport*: unified tiles
reassemble by shared topology only while they stay in one process, and
dispatching them to workers as `.brep` files gives every seam edge two copies.

Threads would sidestep that completely. One process, one heap, tiles pointing at
literally the same `TShape` objects, no serialization anywhere — the reassembly
G13 and G14 measured would just work, because it would be the same reassembly.

So the obvious question is whether the transport was ever the real obstacle.
Three things have to hold, and this gate measures the first two directly rather
than inheriting them:

**1. Does OCP release the GIL?** G7 recorded that it does not, for
``unify_same_domain`` and ``is_valid``, on the strength of a 0.91–1.01x speedup
over 6 threads. That is a *symptom* measured once, on a different call mix, on
an older build. Probe A tests the mechanism head-on: spin a Python counter in a
background thread while one OCCT call runs on the main thread. If the counter
advances, the GIL is released during the call and threading is on the table. If
it stalls, no arrangement of threads can help, because only one of them can be
inside Python — or inside OCP — at a time.

**2. Does threading actually go faster?** Probe B tiles a real solid and unifies
the tiles across a `ThreadPoolExecutor`, against the same work done serially. It
also counts free edges after reassembly, because the whole appeal of threads is
that this should be zero without any repair.

**3. Is it safe?** This gate cannot answer that, and says so rather than
implying otherwise. Tiles that share a seam edge share a `TShape`, so two
threads unifying neighbouring tiles touch the same reference count and the same
object concurrently — precisely the case OCCT's thread-safety notes do not
cover. A green run here is *not* evidence of safety: a data race that does not
fire is indistinguishable from no race at all. That question only becomes worth
asking if probes A and B both pass.

    python tools/prototypes/g17_thread_tiled_unify.py
    python tools/prototypes/g17_thread_tiled_unify.py --m 12 --threads 6
"""

import argparse
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import _bootstrap  # noqa: F401

import g13_unify_scaling as g13
from latticegen2 import occ
from latticegen2.interior import extract_template_mesh
from latticegen2.junction import build_template
from latticegen2.lattice import lattice_params
from latticegen2.weld import free_edges

CC, T = 10.0, 1.5


def probe_gil(solid) -> float:
    """Fraction of a free-running counter's rate that survives during one OCCT call.

    Calibrated against the same counter with nothing else running, so the figure
    is "how much Python still gets done", not a raw count that would depend on
    the machine. ~1.0 means the GIL is released for the duration of the call;
    ~0.0 means it is held and threads cannot overlap with it at all.
    """
    stop = threading.Event()
    counts = [0]

    def spin():
        while not stop.is_set():
            counts[0] += 1

    # Baseline: the counter alone, no OCCT work competing with it.
    t = threading.Thread(target=spin, daemon=True)
    t.start()
    time.sleep(0.5)
    stop.set()
    t.join()
    baseline_rate = counts[0] / 0.5

    stop = threading.Event()
    counts[0] = 0
    t = threading.Thread(target=spin, daemon=True)
    t.start()
    t0 = time.perf_counter()
    occ.unify_same_domain(solid, unify_edges=False)
    elapsed = time.perf_counter() - t0
    stop.set()
    t.join()

    during_rate = counts[0] / elapsed
    print(f"  OCCT call took {elapsed:.2f}s")
    print(f"  counter alone      : {baseline_rate:>12,.0f} iterations/s")
    print(f"  counter during call: {during_rate:>12,.0f} iterations/s")
    return during_rate / baseline_rate if baseline_rate else 0.0


def unify_tiles(tiles, workers: int):
    """Unify each tile's faces, serially or across a thread pool."""

    def one(tile):
        return occ.faces(occ.unify_same_domain(occ.faces_shell(tile), unify_edges=False))

    t0 = time.perf_counter()
    if workers <= 1:
        results = [one(tile) for tile in tiles]
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            results = list(pool.map(one, tiles))
    return time.perf_counter() - t0, [f for r in results for f in r]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--m", type=int, default=12, help="grid edge in cells")
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--tiles", type=int, default=3, help="tiles per axis (n^3 total)")
    args = ap.parse_args()

    print("G17 — threads instead of processes for a tiled unification\n")
    lp = lattice_params(CC, T)
    tpl = build_template(lp)
    tmesh = extract_template_mesh(lp, tpl)
    solid = g13.build_solid(lp, tpl, tmesh, args.m)
    all_faces = occ.faces(solid)
    print(f"input: {len(all_faces)} faces, {len(free_edges(all_faces))} free edges\n")

    print("--- probe A: does OCP release the GIL during unify_same_domain? ---")
    share = probe_gil(solid)
    print(f"  Python throughput retained during the call: {100 * share:.1f}%")
    released = share > 0.5
    print(f"  => GIL is {'RELEASED' if released else 'HELD'}\n")

    print(f"--- probe B: {args.tiles ** 3} tiles, serial against {args.threads} threads ---")
    tiles = g13.tile_faces(all_faces, args.tiles)
    t_serial, faces_serial = unify_tiles(tiles, 1)
    t_threaded, faces_threaded = unify_tiles(tiles, args.threads)
    fe = len(free_edges(faces_threaded))
    speedup = t_serial / t_threaded if t_threaded else float("nan")
    print(f"  {len(tiles)} tiles, serial   : {t_serial:>8.3f}s  {len(faces_serial)} faces")
    print(f"  {len(tiles)} tiles, threaded : {t_threaded:>8.3f}s  {len(faces_threaded)} faces")
    print(f"  speedup on {args.threads} threads: {speedup:.2f}x")
    print(f"  free edges after reassembly: {fe}")

    print("\n=== VERDICT ===")
    if not released and speedup < 1.5:
        print("Threads cannot help: OCP holds the GIL for the whole call, so the")
        print("tiles run one at a time however many threads dispatch them. This")
        print("confirms G7's finding by mechanism rather than by symptom, and it")
        print("closes the last route to sub-body parallelism in this stage —")
        print("processes break tile identity (G15), threads do not run.")
    elif released and speedup >= 1.5:
        print(f"Threads DO overlap ({speedup:.2f}x) and reassembly left {fe} free")
        print("edge(s). If that is 0, the transport was the only obstacle and")
        print("docs/specification.md §11's Phase 3 should be reopened — but only")
        print("after OCCT thread-safety on SHARED seam TShapes is established.")
        print("This gate does not test that, and a run that happens not to crash")
        print("is not evidence for it.")
    else:
        print(f"Mixed: GIL {'released' if released else 'held'}, speedup "
              f"{speedup:.2f}x, {fe} free edges. Re-run before drawing anything")
        print("from it; a single timing on a small solid is not a scaling result.")


if __name__ == "__main__":
    main()
