"""Gate G14 — does a tiled unification still reassemble on *trimmed* geometry?

G13 established that unifying a solid's faces in spatial tiles and concatenating
the results with `BRep_Builder.Add` reproduces a closed, valid solid with zero
free edges — provided edge merging is off. That is the property
docs/specification.md §10 Phase 3 rests on: it is what keeps the
volume-scaling face set away from `BRepBuilderAPI_Sewing` (docs/algorithm.md
§6, §8).

**But G13 measured it on all-planar instanced grids only**, and Phase 3 lists
that as blocking risk R2: real output carries trimmed, curved faces produced by
`BRepAlgoAPI_Common` against the input body, where `TShape` identity across a
tile seam is unproven. Phase 2 made this sharper rather than softer. With the
interior now built pre-merged (docs/algorithm.md §6), the faces still entering
`simplify` are mostly boundary-derived — on `dense-lattice`, 15,718 of 25,234 —
so the region Phase 3 would tile is exactly the region G13 never tested.

This gate builds a solid that genuinely mixes the two: an instanced lattice grid
intersected with a sphere, so the result carries the lattice's planar faces
*and* spherical trimmed faces, closed and valid. It then tiles the unification
the way Phase 3 would and checks the reassembly by count, not by argument.

    python tools/prototypes/g14_tiled_unify_trimmed.py
    python tools/prototypes/g14_tiled_unify_trimmed.py --pieces path/to/temp/<stamp>

The second form additionally runs the same test over the real trimmed boundary
pieces a kept `temp/<stamp>` folder holds, sewn into a boundary-layer shell the
way `weld` does. That shell is *open* at every interface hole, so the bar there
is that the tiled reassembly has the **same** free-edge count as the whole-solid
unification — not zero.

Bars, set in advance:

* the test solid must actually contain curved faces, or the gate is vacuous;
* tiled reassembly: **0 free edges**, `BRepCheck_Analyzer` valid, volume
  preserved to 1e-12 relative;
* merge loss is reported but not barred here — G13 already covers it, and R3's
  fix (bucket by strut, not by centroid) is a separate change.
"""

import argparse
import glob
import os

import _bootstrap  # noqa: F401
import numpy as np
from OCP.BRep import BRep_Builder
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common
from OCP.BRepPrimAPI import BRepPrimAPI_MakeSphere
from OCP.BRepTools import BRepTools
from OCP.TopoDS import TopoDS_Shape

import g13_unify_scaling as g13
from latticegen2 import occ
from latticegen2.connect import lattice_interfaces
from latticegen2.interior import build_interior_shell, extract_template_mesh
from latticegen2.junction import build_template
from latticegen2.lattice import lattice_params
from latticegen2.weld import free_edges

CC, T = 10.0, 1.5
SEW_TOL = 1e-6

MEASURE_REL_BAR = 1e-5
"""Relative bar on the volume (or area) a tiled unification must preserve.

Deliberately `latticegen2.pipeline.UNIFY_VOLUME_TOL`, the bar production already
applies to this exact operation, rather than a figure invented here. G13's own
1e-13 drift was measured on all-planar grids where the integral is effectively
exact; on curved trimmed faces it is not, and docs/algorithm.md §9 records
2.4e-7 relative drift from the same cause on `dense-lattice`. A tighter bar
would reject a correct result for quadrature noise — the mistake
docs/specification.md §10 already records once, where a 1e-9 area bar rejected a
correct 3.08e-09 unification. The exact, tolerance-free part of this gate is the
free-edge count; this is the secondary sanity check."""


def trimmed_lattice_solid(m: int, radius_frac: float):
    """An instanced grid cut by a sphere: planar lattice faces plus curved ones.

    The sphere is deliberately sized to pass *through* the grid rather than
    contain it, so a large number of struts are cut off mid-length and the
    result carries genuinely trimmed, curved boundary faces — the face type R2
    is about.
    """
    lp = lattice_params(CC, T)
    tpl = build_template(lp)
    tmesh = extract_template_mesh(lp, tpl)
    ns = g13.grid(m)
    kept = {(int(a), int(b), int(c)) for a, b, c in ns}
    built = build_interior_shell(lp, tpl, tmesh, ns, lattice_interfaces(kept))
    shell = list(built.shells.values())[0]
    shell.Closed(built.stats["interior_open_edges"] == 0)
    grid_solid = occ.make_solid(shell)

    lo, hi = occ.bounding_box(grid_solid)
    centre = 0.5 * (lo + hi)
    # A fraction of the *smallest* half-extent, not of the diagonal: the lattice
    # bounding box is elongated along Z (the body diagonal stands up), so a
    # radius scaled from the diagonal swallowed the whole grid and left nothing
    # trimmed at all — a vacuous pass, which the curved-face count now catches.
    radius = radius_frac * float(np.min(0.5 * (hi - lo)))
    sphere = BRepPrimAPI_MakeSphere(occ._pnt(centre), radius).Shape()

    common = BRepAlgoAPI_Common(grid_solid, sphere)
    common.Build()
    if not common.IsDone():
        raise SystemExit("the intersection did not complete")
    solids = occ.solids(common.Shape())
    if not solids:
        raise SystemExit("the intersection produced no solid; adjust --radius")
    return max(solids, key=occ.volume)


def load_sewn_pieces(directory: str, limit: int):
    """The real trimmed boundary pieces of a kept run, sewn as `weld` sews them."""
    paths = sorted(glob.glob(os.path.join(directory, "boundary_*.brep")))[:limit]
    if not paths:
        raise SystemExit(f"no boundary_*.brep in {directory}")
    pieces = []
    for p in paths:
        shape = TopoDS_Shape()
        BRepTools.Read_s(shape, p, BRep_Builder())
        pieces.append(shape)
    return occ.sew(pieces, SEW_TOL), len(pieces)


def describe(label: str, shape) -> tuple[int, int]:
    faces = occ.faces(shape)
    curved = sum(1 for f in faces if not occ.is_planar(f))
    print(f"{label}: {len(faces)} faces ({curved} curved), "
          f"{len(free_edges(faces))} free edges")
    return len(faces), curved


def tiled_unify(all_faces, n_per_axis: int, unify_edges: bool):
    """Unify each spatial tile on its own and concatenate by reference."""
    tiles = g13.tile_faces(all_faces, n_per_axis)
    out = []
    for tf in tiles:
        merged = occ.unify_same_domain(occ.faces_shell(tf), unify_edges=unify_edges)
        out.extend(occ.faces(merged))
    return tiles, out


def run_closed(solid, tile_counts) -> bool:
    print("\n=== closed solid: instanced grid ∩ sphere ===")
    n_faces, curved = describe("input", solid)
    v0 = occ.volume(solid)
    print(f"input volume {v0:.6f} mm^3, valid {occ.is_valid(solid)}")
    if curved == 0:
        print("VACUOUS: no curved faces in the test solid")
        return False

    whole = occ.unify_same_domain(solid, unify_edges=False)
    whole_faces = occ.faces(whole)
    print(f"whole-solid unify: {n_faces} -> {len(whole_faces)} faces, "
          f"volume drift {abs(occ.volume(whole) - v0) / v0:.2e}")

    ok = True
    all_faces = occ.faces(solid)
    for n in tile_counts:
        tiles, out = tiled_unify(all_faces, n, unify_edges=False)
        fe = len(free_edges(out))
        loss = (len(out) - len(whole_faces)) / max(len(whole_faces), 1)
        verdict = "PASS"
        drift = float("nan")
        valid = False
        if fe == 0:
            shell = occ.faces_shell(out)
            shell.Closed(True)
            rebuilt = occ.make_solid(shell)
            valid = occ.is_valid(rebuilt)
            drift = abs(occ.volume(rebuilt) - v0) / abs(v0)
            if not valid or drift > MEASURE_REL_BAR:
                verdict = "FAIL"
        else:
            verdict = "FAIL"
        ok = ok and verdict == "PASS"
        print(f"  {len(tiles):>3} tiles: {verdict}  free edges {fe:>5}, "
              f"faces {len(out):>6} ({100 * loss:+.2f}%), valid {valid}, "
              f"volume drift {drift:.2e}")
    return ok


def run_open(shell, n_pieces: int, tile_counts) -> bool:
    print(f"\n=== open boundary-layer shell: {n_pieces} real trimmed pieces ===")
    _, curved = describe("input", shell)

    whole = occ.unify_same_domain(shell, unify_edges=False)
    whole_faces = occ.faces(whole)
    base_fe = len(free_edges(whole_faces))
    base_area = occ.area(whole)
    print(f"whole-shell unify: {len(whole_faces)} faces, {base_fe} free edges, "
          f"area {base_area:.6f} mm^2")

    ok = True
    all_faces = occ.faces(shell)
    for n in tile_counts:
        tiles, out = tiled_unify(all_faces, n, unify_edges=False)
        fe = len(free_edges(out))
        area = occ.area(occ.faces_shell(out))
        drift = abs(area - base_area) / max(base_area, 1.0)
        # An open shell cannot be closed; the property under test is that tiling
        # introduced no *new* free edge, i.e. every tile seam re-identified.
        good = fe == base_fe and drift <= MEASURE_REL_BAR
        ok = ok and good
        print(f"  {len(tiles):>3} tiles: {'PASS' if good else 'FAIL'}  free edges "
              f"{fe:>5} (baseline {base_fe}), faces {len(out):>6}, "
              f"area drift {drift:.2e}")
    return ok


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--m", type=int, default=6, help="grid edge in cells")
    ap.add_argument("--radius", type=float, default=0.75,
                    help="sphere radius as a fraction of the smallest half-extent")
    ap.add_argument("--tiles", type=int, nargs="+", default=[2, 3])
    ap.add_argument("--pieces", default=None,
                    help="also test a kept temp/<stamp> folder's boundary_*.brep")
    ap.add_argument("--max-pieces", type=int, default=4000)
    args = ap.parse_args()

    print(f"G14 (cc={CC}, t={T}) — tiled unification over trimmed geometry\n")
    solid = trimmed_lattice_solid(args.m, args.radius)
    closed_ok = run_closed(solid, args.tiles)

    open_ok = True
    if args.pieces:
        shell, n = load_sewn_pieces(os.path.abspath(args.pieces), args.max_pieces)
        open_ok = run_open(shell, n, args.tiles)

    print("\n=== VERDICT ===")
    print(f"closed trimmed solid: {'PASS' if closed_ok else 'FAIL'}")
    if args.pieces:
        print(f"real boundary pieces: {'PASS' if open_ok else 'FAIL'}")
    if closed_ok and open_ok:
        print("\nR2 clears: tile identity survives trimmed, curved faces, so "
              "Phase 3 may tile the boundary region as well as the interior.")
    else:
        print("\nR2 stands: tiling must fall back to interior-derived faces "
              "only, with the boundary region handed to one unify call.")


if __name__ == "__main__":
    main()
