"""Gates G3 and G4 — per-junction intersection latency, and STEP writer throughput.

G3: instance the junction template at nodes straddling a real input surface and
time one single-operand ``BRepAlgoAPI_Common`` each.
    Pass bar: median < 50 ms, p95 < 250 ms.

G4: write an instanced interior shell through ``STEPControl_Writer`` and re-read
it.
    Pass bar: >= 3,000 faces/s, and the round trip preserves solid count and
    total volume.
"""

import os
import statistics
import sys
import tempfile
import time

import _bootstrap  # noqa: F401
import numpy as np

from latticegen2 import occ
from latticegen2.boundary import trim_junction
from latticegen2.classify import Klass, classify_nodes, tessellate_surface
from latticegen2.interior import build_interior_shell, extract_template_mesh
from latticegen2.junction import build_template
from latticegen2.lattice import candidate_nodes, lattice_params, nodes

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
CYLINDER = os.path.join(ROOT, "test", "test-cylinder.STEP")


def gate_g3() -> bool:
    lp = lattice_params(10.0, 1.5)
    tpl = build_template(lp)
    body = occ.read_step(CYLINDER)
    lo, hi = occ.bounding_box(body)
    mesh = tessellate_surface(body, lp)
    cls = classify_nodes(lp, mesh, candidate_nodes(lp, lo, hi))
    boundary = cls.of(Klass.BOUNDARY)
    sample = boundary[:: max(1, len(boundary) // 200)][:200]
    positions = nodes(lp, sample)

    timings = []
    for i in range(len(sample)):
        t0 = time.perf_counter()
        trim_junction(lp, tpl, positions[i], body, 0b111111)
        timings.append((time.perf_counter() - t0) * 1e3)

    timings.sort()
    median = statistics.median(timings)
    p95 = timings[int(0.95 * (len(timings) - 1))]
    ok = median < 50.0 and p95 < 250.0
    print(f"G3: {len(timings)} boundary junctions of {len(boundary)}; "
          f"median {median:.1f} ms, p95 {p95:.1f} ms, max {timings[-1]:.1f} ms "
          f"-> {'PASS' if ok else 'FAIL'} (bars: median < 50, p95 < 250)")
    return ok


def gate_g4() -> bool:
    lp = lattice_params(10.0, 1.5)
    tpl = build_template(lp)
    tmesh = extract_template_mesh(lp, tpl)
    m = 12
    a = np.arange(m, dtype=np.int64)
    g = np.meshgrid(a, a, a, indexing="ij")
    ns = np.stack([x.ravel() for x in g], axis=1)
    kept = {(int(x), int(y), int(z)) for x, y, z in ns}
    shell, stats = build_interior_shell(lp, tpl, tmesh, ns, kept)
    solid = occ.make_solid(shell)
    n_faces = stats["interior_faces"]
    volume = occ.volume(solid)

    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "g4.step")
        t0 = time.perf_counter()
        occ.write_step(solid, path, "g4+cc10+t1.5")
        write_s = time.perf_counter() - t0
        size = os.path.getsize(path)
        t0 = time.perf_counter()
        back = occ.read_step(path)
        read_s = time.perf_counter() - t0
        sols = occ.solids(back)
        round_volume = sum(occ.volume(s) for s in sols)

    rate = n_faces / write_s
    vol_ok = abs(round_volume - volume) / volume < 1e-9
    ok = rate >= 3000 and len(sols) == 1 and vol_ok
    print(f"G4: {n_faces} faces, wrote {size / 1e6:.1f} MB in {write_s:.2f}s "
          f"({rate:.0f} faces/s), re-read in {read_s:.2f}s -> {len(sols)} solid(s), "
          f"volume rel.err {abs(round_volume - volume) / volume:.2e} "
          f"-> {'PASS' if ok else 'FAIL'} (bar: >= 3000 faces/s, exact round trip)")
    return ok


def main() -> int:
    occ.quiet_kernel()
    results = [gate_g3(), gate_g4()]
    print()
    print("G3/G4: PASS" if all(results) else "G3/G4: FAIL")
    return 0 if all(results) else 1


if __name__ == "__main__":
    sys.exit(main())
