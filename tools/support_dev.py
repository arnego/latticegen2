"""Run the support generator's stages on a part, with or without a lattice.

A development harness, not a product feature: the product always generates the
lattice (specification.md §4.6). Without ``--lattice`` this drives the generator
on parts whose only anchors are their own surfaces -- `test/Hole-block-*.STEP`
at `--orient 90 0 0`. With it, the real pipeline is run once to the support
hook and the completed lattice and its kept nodes are cached beside the output
(``<out>-lattice.brep``, ``<out>-lattice-nodes.npz``), so the generator can be
rerun against them in seconds. Either way it writes what each stage produced as
STEP files for inspection in CAD:

* ``<out>-supports.step`` the support bodies, as the generator returns them
* ``<out>-input.step``   the input body, for context beside them
* ``<out>-summary.log``  the generator's own log lines and statistics

Usage::

    python tools/support_dev.py -i test/Hole-block-basic.STEP --orient 90 0 0 \\
        -o out/hole-block-basic [-t 1.5] [--support-thickness 1.5]

Without ``--lattice``, ``-t`` sets only the classification mesh's deflection
(``min(t, a)/10``), the grid and the default thickness; there are no struts,
and every anchor is on a wall.
"""

from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "src"))

from latticegen2 import occ  # noqa: E402
from latticegen2.classify import tessellate_surface  # noqa: E402
from latticegen2.lattice import lattice_params  # noqa: E402
from latticegen2.support import overhang  # noqa: E402
from latticegen2.support.params import STEP  # noqa: E402


def _lattice(args):
    """``(solids, kept_nodes)`` of the lattice the product builds for these
    arguments, from the cache beside the output or from one pipeline run.

    The run is the product's own, stopped short at nothing: the generator is
    swapped for one that records what it is handed and returns no supports, and
    the run then finishes as an unsupported one would.
    """
    from latticegen2.parallel import read_brep, write_brep

    brep = args.output + "-lattice.brep"
    npz = args.output + "-lattice-nodes.npz"
    if not (os.path.isfile(brep) and os.path.isfile(npz)):
        from latticegen2 import __main__ as entry
        from latticegen2.support import generate as gen
        from latticegen2.support.contract import SupportResult

        def record(inp):
            write_brep(occ.compound(list(inp.lattice)), brep)
            np.savez(npz, kept_nodes=np.asarray(inp.kept_nodes))
            return SupportResult(supports=[])

        real, gen.generate = gen.generate, record
        try:
            entry.main(["-i", args.input, "-cc", str(args.cc), "-t", str(args.t),
                        "--orient", *[str(a) for a in args.orient], "--support",
                        "--cores", str(args.cores), "-o", args.output + "-lattice-run.step"])
        finally:
            gen.generate = real
        if not os.path.isfile(brep):
            raise SystemExit("the pipeline did not reach the support stage")
    return occ.solids(read_brep(brep)), np.load(npz)["kept_nodes"]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", required=True, help="output path stem")
    ap.add_argument("--orient", nargs=3, type=float, default=(0.0, 0.0, 0.0))
    ap.add_argument("-t", type=float, default=1.5)
    ap.add_argument("-cc", type=float, default=10.0,
                    help="the lattice frame; and the lattice, with --lattice")
    ap.add_argument("--lattice", action="store_true",
                    help="build the product's lattice (cached) and anchor on it")
    ap.add_argument("--cores", type=int, default=6)
    ap.add_argument("--support-thickness", type=float, default=None)
    args = ap.parse_args(argv)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    lines = []

    def say(s):
        print(s, flush=True)
        lines.append(s)

    t0 = time.time()
    body = occ.read_step(args.input)
    lp = lattice_params(args.cc, args.t, tuple(args.orient))
    mesh = tessellate_surface(body, lp)
    say(f"input {args.input}, orient {tuple(args.orient)}, b = {np.round(lp.b, 9).tolist()}")
    samples = overhang.detect(body, lp.b, STEP)
    say(f"overhang: {len(samples)} sample(s), {samples.area:.2f} mm^2 "
        f"({samples.counts()})")

    lattice, kept = [], np.empty((0, 3), dtype=np.int64)
    if args.lattice:
        t2 = time.time()
        lattice, kept = _lattice(args)
        say(f"lattice: {len(lattice)} solid(s), {len(kept)} kept node(s) "
            f"({time.time() - t2:.1f} s)")
    thickness = args.support_thickness or args.t
    t5 = time.time()
    from latticegen2.support.generate import supports_for

    solids, gstats, glines = supports_for(body, mesh, lattice, kept, lp, thickness,
                                          tmpdir=os.path.dirname(os.path.abspath(args.output)))
    for line in glines:
        say(line)
    vols = [occ.volume(s) for s in solids] or [0.0]
    say(f"supports: {len(solids)} bodies, {sum(vols):.2f} mm^3 "
        f"(min {min(vols):.2f}, max {max(vols):.2f}), faces per body "
        f"{[len(occ.faces(s)) for s in solids]}, {gstats} ({time.time() - t5:.1f} s)")
    occ.write_step(occ.compound(solids), args.output + "-supports.step", "supports")
    occ.write_step(body, args.output + "-input.step", "input")
    say(f"wrote {args.output}-supports.step and -input.step ({time.time() - t0:.1f} s in all)")
    with open(args.output + "-summary.log", "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
