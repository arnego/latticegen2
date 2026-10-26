"""Run the support generator's stages on a part **without a lattice**.

A development harness, not a product feature: the product always generates the
lattice (specification.md §4.6). This drives the generator on parts whose only
anchors are their own surfaces -- `test/Hole-block-*.STEP` at `--orient 90 0 0`
-- and writes what each stage produced as STEP files for inspection in CAD:

* ``<out>-rays.step``    first-contact rays, apex to target, grouped by region;
                         the region boundary loops as polylines
* ``<out>-cells.step``   one body per cell: its perimeter loops as polylines and
                         the cell's own rays from them
* ``<out>-supports.step`` the support bodies, one per cell
* ``<out>-input.step``   the input body, for context beside them
* ``<out>-summary.log``  counts, areas, ray drops and lengths

Usage::

    python tools/support_dev.py -i test/Hole-block-basic.STEP --orient 90 0 0 \\
        -o inspection/p2-hole-block-basic [-t 1.5] [--support-thickness 1.5]

``-t`` sets only the classification mesh's deflection (``min(t, a)/10``) and the
default thickness; there are no struts.
"""

from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "src"))

from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakePolygon  # noqa: E402
from OCP.gp import gp_Pnt  # noqa: E402

from latticegen2 import occ  # noqa: E402
from latticegen2.classify import tessellate_surface  # noqa: E402
from latticegen2.lattice import lattice_params  # noqa: E402
from latticegen2.support import cells as cells_  # noqa: E402
from latticegen2.support import obstacles as obst  # noqa: E402
from latticegen2.support import overhang, rays, regions  # noqa: E402
from latticegen2.support.params import REACH, STEP  # noqa: E402


def _pnt(p):
    return gp_Pnt(float(p[0]), float(p[1]), float(p[2]))


def _segments(a: np.ndarray, b: np.ndarray):
    out = []
    for p, q in zip(a, b):
        if np.linalg.norm(q - p) > 1e-9:
            out.append(BRepBuilderAPI_MakeEdge(_pnt(p), _pnt(q)).Edge())
    return occ.compound(out)


def _polyline(points: np.ndarray):
    poly = BRepBuilderAPI_MakePolygon()
    for p in points:
        poly.Add(_pnt(p))
    poly.Close()
    return poly.Wire()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", required=True, help="output path stem")
    ap.add_argument("--orient", nargs=3, type=float, default=(0.0, 0.0, 0.0))
    ap.add_argument("-t", type=float, default=1.5)
    ap.add_argument("-cc", type=float, default=10.0,
                    help="only shapes the lattice frame; no lattice is built")
    ap.add_argument("--support-thickness", type=float, default=None)
    ap.add_argument("--every", type=int, default=1,
                    help="write every n-th ray only, to keep the file light")
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

    t1 = time.time()
    rmesh = regions.build(body, lp.b)
    regs = regions.regions(rmesh)
    say(f"region mesh: {len(rmesh.verts)} vertices, {len(rmesh.tris)} triangles, "
        f"{int(rmesh.required.sum())} required, {len(rmesh.isocline)} isocline points "
        f"({time.time() - t1:.1f} s)")
    for k, r in enumerate(regs):
        say(f"  region {k}: {r.area:.2f} mm^2, {len(r.verts)} vertices, "
            f"{len(r.loops)} boundary loop(s) of {[len(l) for l in r.loops]}")

    t2 = time.time()
    obstacles = obst.from_shapes(body, [], deflection=STEP / 10.0)
    say(f"obstacles: {len(obstacles.mesh.tris)} triangles ({time.time() - t2:.1f} s)")

    ray_edges = []
    loop_wires = []
    cell_shapes = []
    thickness = args.support_thickness or args.t
    for k, r in enumerate(regs):
        t3 = time.time()
        idx = r.verts
        apex = rmesh.verts[idx] - rays.APEX_NUDGE * rmesh.normals[idx]
        got = rays.first_contact(apex, lp.b, obstacles)
        exempt = got.length <= REACH
        say(f"  region {k}: {len(got)} rays in {time.time() - t3:.1f} s; "
            f"{int(exempt.sum())} exempt (<= {REACH} mm); drop "
            f"min {got.drop.min():.3f} median {np.median(got.drop):.3f} "
            f"max {got.drop.max():.3f} mm; length max {got.length.max():.3f} mm; "
            f"anchors {sorted(set(got.anchor.tolist()))}")
        keep = np.nonzero(~exempt)[0][::max(1, args.every)]
        ray_edges.append(_segments(got.apex[keep], got.target[keep]))
        for loop in r.loops:
            loop_wires.append(_polyline(rmesh.verts[loop]))
        t4 = time.time()
        part = cells_.partition(rmesh, r, got, obstacles, lp.b, thickness)
        say(f"  region {k}: {part.stats} ({time.time() - t4:.1f} s)")
        for c in sorted(part.cells, key=lambda c: -c.area):
            say(f"    cell {c.id}: anchor group {c.key}, {c.area:.2f} mm^2, loops "
                f"{[len(lp_.points) for lp_ in c.loops]}, narrowest parting "
                f"{c.min_divergence_deg:.1f} deg")
            parts = []
            for lp_ in c.loops:
                parts.append(_polyline(lp_.points))
                parts.append(_segments(lp_.points, lp_.targets))
            cell_shapes.append(occ.compound(parts))
    occ.write_step(occ.compound(ray_edges + loop_wires),
                   args.output + "-rays.step", "support-dev-rays")
    occ.write_step(occ.compound(cell_shapes), args.output + "-cells.step",
                   "support-dev-cells")
    say(f"wrote {args.output}-rays.step and -cells.step ({time.time() - t0:.1f} s total)")

    t5 = time.time()
    from latticegen2.support.generate import supports_for

    solids, gstats, glines = supports_for(body, mesh, [], np.empty((0, 3), dtype=np.int64),
                                          lp, thickness)
    for line in glines:
        say(line)
    vols = [occ.volume(s) for s in solids]
    say(f"supports: {len(solids)} bodies, {sum(vols):.2f} mm^3 "
        f"(min {min(vols):.2f}, max {max(vols):.2f}), faces per body "
        f"max {max(len(occ.faces(s)) for s in solids)}, {gstats} "
        f"({time.time() - t5:.1f} s)")
    occ.write_step(occ.compound(solids), args.output + "-supports.step", "supports")
    occ.write_step(body, args.output + "-input.step", "input")
    say(f"wrote {args.output}-supports.step and -input.step")
    with open(args.output + "-summary.log", "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
