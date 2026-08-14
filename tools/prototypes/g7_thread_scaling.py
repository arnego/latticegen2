"""Gate G7 — does OCP release the GIL around `unify_same_domain` / `is_valid`?

docs/specification.md §10 ranks parallelising `simplify` (same-domain
unification) and `validate` (`BRepCheck_Analyzer`) as optimization paths 1 and
4 — together ~16.5 min of the `cc=5, t=1` rehearsal's 73.1. Every other
parallel stage in this codebase (boundary trim, boundary sew) uses a
`multiprocessing` spawn-context pool with a `.brep` file round-trip, because
that is the only mechanism *known* to work — but for `simplify`/`validate` the
payload is the run's entire output, up to 2 GB, so a process-pool round-trip
adds a serial master-side read-back that did not exist before.

If OCP releases the GIL around these two calls — plausible, since they are
long-running calls into native OCCT code, exactly the shape of call pybind11
extensions typically release the GIL for — a plain `ThreadPoolExecutor` gets
the same wall-clock win with zero IPC, zero serialization, and zero extra
memory. This gate measures it directly rather than assuming either way, per
docs/algorithm.md §11.

Bar: threads win outright if threaded wall time is under ~60% of serial (i.e.
close to what >=2 real cores running is meant to look like on this 6-core
machine); "GIL held" is the presumption otherwise, and the process+.brep path
is used instead.

    python tools/prototypes/g7_thread_scaling.py
"""

import time
from concurrent.futures import ThreadPoolExecutor

import _bootstrap  # noqa: F401
import numpy as np

from latticegen2 import occ
from latticegen2.connect import lattice_interfaces
from latticegen2.interior import build_interior_shell, extract_template_mesh
from latticegen2.junction import build_template
from latticegen2.lattice import lattice_params, nodes

CC, T = 10.0, 1.5
WORKERS = 6


def grid(m: int) -> np.ndarray:
    a = np.arange(m, dtype=np.int64)
    g = np.meshgrid(a, a, a, indexing="ij")
    return np.stack([x.ravel() for x in g], axis=1)


def build_solid(lp, tpl, tmesh, m: int):
    """One valid, reasonably complex solid: an m x m x m all-interior grid."""
    ns = grid(m)
    kept = {(int(a), int(b), int(c)) for a, b, c in ns}
    built = build_interior_shell(lp, tpl, tmesh, ns, lattice_interfaces(kept))
    shell = list(built.shells.values())[0]
    shell.Closed(built.stats["interior_open_edges"] == 0)
    return occ.make_solid(shell)


def make_solids() -> list:
    """A handful of unequal solids, echoing the rehearsal's 14 very-unequal ones."""
    lp = lattice_params(CC, T)
    tpl = build_template(lp)
    tmesh = extract_template_mesh(lp, tpl)
    sizes = [3, 4, 4, 5, 5, 6, 6, 7]  # unequal, biggest last is not assumed
    return [build_solid(lp, tpl, tmesh, m) for m in sizes]


def bench(label: str, fn, solids: list) -> float:
    t0 = time.perf_counter()
    fn(solids)
    return time.perf_counter() - t0


def run_unify_serial(solids: list) -> None:
    for s in solids:
        occ.unify_same_domain(s)


def run_unify_threaded(solids: list) -> None:
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        list(ex.map(occ.unify_same_domain, solids))


def run_validate_serial(solids: list) -> None:
    for s in solids:
        occ.is_valid(s)


def run_validate_threaded(solids: list) -> None:
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        list(ex.map(occ.is_valid, solids))


def main() -> None:
    print(f"Building test solids (cc={CC}, t={T}) ...")
    solids = make_solids()
    faces = [occ.count_subshapes(s)[0] for s in solids]
    print(f"  {len(solids)} solids, face counts: {faces}")

    print("\n--- unify_same_domain ---")
    t_serial = bench("serial", run_unify_serial, solids)
    solids2 = make_solids()  # fresh copies: unify mutates via replacement, not in place, but be safe
    t_threaded = bench("threaded", run_unify_threaded, solids2)
    ratio = t_threaded / t_serial if t_serial else float("nan")
    print(f"  serial:   {t_serial:.3f}s")
    print(f"  threaded: {t_threaded:.3f}s ({WORKERS} workers)")
    print(f"  ratio (threaded/serial): {ratio:.2f}  [<0.6 => GIL released, threads win]")
    unify_threads_win = ratio < 0.6

    print("\n--- is_valid (BRepCheck_Analyzer) ---")
    t_serial_v = bench("serial", run_validate_serial, solids)
    t_threaded_v = bench("threaded", run_validate_threaded, solids)
    ratio_v = t_threaded_v / t_serial_v if t_serial_v else float("nan")
    print(f"  serial:   {t_serial_v:.3f}s")
    print(f"  threaded: {t_threaded_v:.3f}s ({WORKERS} workers)")
    print(f"  ratio (threaded/serial): {ratio_v:.2f}  [<0.6 => GIL released, threads win]")
    validate_threads_win = ratio_v < 0.6

    print("\n=== VERDICT ===")
    print(f"unify_same_domain:  {'THREADS' if unify_threads_win else 'PROCESSES'}")
    print(f"is_valid:           {'THREADS' if validate_threads_win else 'PROCESSES'}")


if __name__ == "__main__":
    main()
