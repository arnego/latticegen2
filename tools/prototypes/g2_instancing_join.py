"""Gate G2 — how the interior lattice gets joined, and how that scales.

Three candidate mechanisms, over an ``m x m x m`` grid of junctions:

  (a) ``indexed``  — build one shell directly, sharing ``TopoDS_Vertex`` and
      ``TopoDS_Edge`` objects between neighbours via the template's precomputed
      cap correspondence (``latticegen2.interior``);
  (b) ``sew``      — instance the template's open shell and let
      ``BRepBuilderAPI_Sewing`` discover the pairing;
  (c) ``glue``     — ``BOPAlgo_Builder`` with ``SetGlue(BOPAlgo_GlueFull)``, OCCT's
      dedicated mode for operands meeting only on coincident faces.

Pass bar: some mechanism joins 6.4e4 junctions in under 5 minutes and under
16 GB RSS, producing one valid solid whose volume is exactly
``n_nodes * volume(J)`` — adjacent junctions have zero-volume contact, so the
union's volume is the plain sum, which makes this an exact, independent check
rather than an approximate one.
"""

import sys
import time

import _bootstrap  # noqa: F401
import numpy as np
from OCP.BOPAlgo import BOPAlgo_Builder, BOPAlgo_GlueEnum
from OCP.TopTools import TopTools_ListOfShape

from latticegen2 import occ
from latticegen2.connect import lattice_interfaces
from latticegen2.interior import build_interior_shell, extract_template_mesh

def _one_shell(built):
    """The single shell a synthetic all-interior grid produces."""
    shells = list(built.shells.values())
    assert len(shells) == 1, shells
    shells[0].Closed(built.stats["interior_open_edges"] == 0)
    return shells[0]

from latticegen2.junction import build_template
from latticegen2.lattice import lattice_params, neighbor_step, nodes

SEW_TOL = 1e-6
GLUE_MAX = 1000
"""Glue mode is only measured at the smallest size: it does not merge at all
(measured below), so larger runs would only burn time confirming that."""
SEW_MAX = 1000
"""Sewing is only measured at the smallest size once its scaling is known:
14.9 s at 1,000 junctions, still running after 250 s of CPU at 8,000."""


def rss_gb() -> float:
    from latticegen2.runlog import peak_rss_bytes

    return peak_rss_bytes() / 1e9


def grid(m: int) -> np.ndarray:
    a = np.arange(m, dtype=np.int64)
    g = np.meshgrid(a, a, a, indexing="ij")
    return np.stack([x.ravel() for x in g], axis=1)


def run_indexed(lp, tpl, tmesh, ns):
    kept = {(int(a), int(b), int(c)) for a, b, c in ns}
    t0 = time.perf_counter()
    built = build_interior_shell(lp, tpl, tmesh, ns, lattice_interfaces(kept))
    shell, stats = _one_shell(built), built.stats
    t_build = time.perf_counter() - t0
    if not shell.Closed():
        # A grid built this way is closed only if every outward cap was capped
        # off; report it rather than silently making an invalid solid.
        return None, t_build, 0.0, stats["interior_faces"]
    t0 = time.perf_counter()
    solid = occ.make_solid(shell)
    t_close = time.perf_counter() - t0
    return [solid], t_build, t_close, stats["interior_faces"]


def _open_shell_pieces(lp, tpl, ns):
    """Template open shells plus the outward caps a finite grid needs."""
    positions = nodes(lp, ns)
    kept = {(int(a), int(b), int(c)) for a, b, c in ns}
    steps = [neighbor_step(h) for h in range(6)]
    pieces = []
    for i in range(len(ns)):
        loc = occ.translation(positions[i])
        pieces.append(tpl.open_shell.Moved(loc))
        n = ns[i]
        for h in range(6):
            nb = (int(n[0] + steps[h][0]), int(n[1] + steps[h][1]), int(n[2] + steps[h][2]))
            if nb not in kept:
                pieces.append(tpl.cap_faces[h].Moved(loc))
    return pieces


def run_sew(lp, tpl, tmesh, ns):
    t0 = time.perf_counter()
    pieces = _open_shell_pieces(lp, tpl, ns)
    t_build = time.perf_counter() - t0
    t0 = time.perf_counter()
    sewn = occ.sew(pieces, SEW_TOL)
    t_join = time.perf_counter() - t0
    sols = [occ.make_solid(s) for s in occ.shells(sewn) if s.Closed()]
    return sols, t_build, t_join, len(occ.faces(sewn))


def run_glue(lp, tpl, tmesh, ns):
    positions = nodes(lp, ns)
    t0 = time.perf_counter()
    args = TopTools_ListOfShape()
    for p in positions:
        args.Append(tpl.solid.Moved(occ.translation(p)))
    t_build = time.perf_counter() - t0
    t0 = time.perf_counter()
    builder = BOPAlgo_Builder()
    builder.SetArguments(args)
    builder.SetGlue(BOPAlgo_GlueEnum.BOPAlgo_GlueFull)
    builder.SetRunParallel(True)
    builder.Perform()
    t_join = time.perf_counter() - t0
    shape = builder.Shape()
    return occ.solids(shape), t_build, t_join, len(occ.faces(shape))


def main() -> int:
    lp = lattice_params(10.0, 1.5)
    tpl = build_template(lp)
    tmesh = extract_template_mesh(lp, tpl)
    print(f"template: {tpl.n_faces} faces, {len(tmesh.verts)} vertices, "
          f"volume {tpl.volume:.6f} mm^3")
    print(f"{'nodes':>8} {'method':>8} {'build s':>9} {'join s':>9} {'solids':>7} "
          f"{'faces':>9} {'vol err':>10} {'valid':>6} {'RSS GB':>7}")

    ok = True
    for m in (10, 20, 40):
        ns = grid(m)
        expected = len(ns) * tpl.volume
        for name, fn, limit in (
            ("indexed", run_indexed, None),
            ("sew", run_sew, SEW_MAX),
            ("glue", run_glue, GLUE_MAX),
        ):
            if limit is not None and len(ns) > limit:
                print(f"{len(ns):8d} {name:>8}   skipped above {limit} nodes (see module docstring)")
                continue
            try:
                sols, t_build, t_join, nfaces = fn(lp, tpl, tmesh, ns)
            except Exception as exc:
                print(f"{len(ns):8d} {name:>8} ERROR {type(exc).__name__}: {exc}")
                ok = False
                continue
            if not sols:
                print(f"{len(ns):8d} {name:>8} {t_build:9.2f} {t_join:9.2f} "
                      f"{0:7d} {nfaces:9d}   no closed solid produced")
                ok = False
                continue
            vol = sum(occ.volume(s) for s in sols)
            rel = abs(vol - expected) / expected
            valid = all(occ.is_valid(s) for s in sols)
            print(f"{len(ns):8d} {name:>8} {t_build:9.2f} {t_join:9.2f} {len(sols):7d} "
                  f"{nfaces:9d} {rel:10.2e} {str(valid):>6} {rss_gb():7.2f}")
            if name == "indexed" and m == 40:
                if t_build + t_join > 300 or rel > 1e-9 or len(sols) != 1 or not valid:
                    ok = False
    print()
    print("G2: PASS" if ok else "G2: FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
