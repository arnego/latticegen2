"""Gate G15 — does tile identity survive the `.brep` round trip Phase 3 needs?

G13 and G14 established the property docs/specification.md §10 Phase 3 rests on:
unify a solid's faces in spatial tiles with ``unify_edges=False``, concatenate
the results with ``BRep_Builder.Add``, and you get the whole-solid shell back —
0 free edges, valid, volume preserved. That works because
``ShapeUpgrade_UnifySameDomain`` leaves the edges it did not merge *as the same
objects*, so the two sides of a tile seam remain literally one ``TShape`` and
nothing has to be re-identified.

**Both gates measured that in one process, where identity is pointer identity.**
Phase 3 does not run in one process. Its whole purpose is to spread tiles across
:class:`latticegen2.parallel.WorkerPool`, and every stage that does so moves
geometry as ``.brep`` files — small-IPC discipline, docs/algorithm.md §7, §8.
A ``.brep`` file is a serialization: sharing is preserved *within* one file and
cannot possibly be preserved *between* two, because each file writes its own
copy of every edge it references. Two tiles that share a seam edge are written
by two different workers into two different files.

So the question this gate asks is the one neither G13 nor G14 could: **after the
round trip, is the seam still one edge, or two?** If it is two, tiles no longer
reassemble by ``BRep_Builder.Add`` alone, and Phase 3 needs a seam-repair step
whose cost has to be measured before anything is built on it — the same shape of
finding G8 recorded when a property proven at prototype scale failed on real
production geometry.

The test is binary and needs no tolerance, exactly as G13's was. The input is a
closed solid, so the whole-solid unification has zero free edges by
construction; count edges used by exactly one face after each reassembly route.

Three routes are measured, and the first two are controls rather than
candidates:

* **in-process** — reproduces G14. Must be 0, or the harness is wrong and
  nothing else here can be read.
* **one file, all tiles** — every tile's result written to a *single* ``.brep``
  as one compound and read back. This is what a single worker handling every
  tile would produce, and it isolates "does serialization preserve sharing at
  all" from "does it preserve sharing across files".
* **one file per tile** — what Phase 3 actually does.

    python tools/prototypes/g15_tiled_unify_ipc.py
    python tools/prototypes/g15_tiled_unify_ipc.py --tiles 2 3 --m 6
"""

import argparse
import os
import tempfile

import _bootstrap  # noqa: F401

import g13_unify_scaling as g13
import g14_tiled_unify_trimmed as g14
from latticegen2 import occ
from latticegen2.parallel import compound_children, read_brep, write_brep
from latticegen2.weld import free_edges


def merged_tiles(all_faces, n_per_axis: int) -> list[list]:
    """Unify each spatial tile on its own, returning each tile's faces separately.

    Kept per tile rather than flattened because the whole point of this gate is
    what happens *between* tiles; flattening first would throw away the boundary
    the round trip is being tested across.
    """
    out = []
    for tf in g13.tile_faces(all_faces, n_per_axis):
        merged = occ.unify_same_domain(occ.faces_shell(tf), unify_edges=False)
        out.append(occ.faces(merged))
    return out


def route_in_process(tiles: list[list]) -> list:
    """G14's route: concatenate the live face objects."""
    return [f for tile in tiles for f in tile]


def route_one_file(tiles: list[list], tmpdir: str) -> list:
    """Every tile into one ``.brep`` as a compound of compounds, then back."""
    path = os.path.join(tmpdir, "all_tiles.brep")
    write_brep(occ.compound([occ.compound(tile) for tile in tiles]), path)
    return [f for child in compound_children(read_brep(path)) for f in occ.faces(child)]


def route_file_per_tile(tiles: list[list], tmpdir: str) -> list:
    """One ``.brep`` per tile — what dispatching tiles across workers produces."""
    out = []
    for i, tile in enumerate(tiles):
        path = os.path.join(tmpdir, f"tile_{i}.brep")
        write_brep(occ.compound(tile), path)
        out.extend(occ.faces(read_brep(path)))
    return out


ROUTES = [
    ("in-process (G14's route)", route_in_process, False),
    ("one file, all tiles", route_one_file, True),
    ("one file per tile", route_file_per_tile, True),
]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--m", type=int, default=6, help="grid edge in cells")
    ap.add_argument("--radius", type=float, default=0.75,
                    help="sphere radius as a fraction of the smallest half-extent")
    ap.add_argument("--tiles", type=int, nargs="+", default=[2, 3])
    args = ap.parse_args()

    print("G15 — does a tiled unification survive the .brep round trip?\n")
    solid = g14.trimmed_lattice_solid(args.m, args.radius)
    all_faces = occ.faces(solid)
    curved = sum(1 for f in all_faces if not occ.is_planar(f))
    v0 = occ.volume(solid)
    print(f"input: {len(all_faces)} faces ({curved} curved), "
          f"{len(free_edges(all_faces))} free edges, volume {v0:.6f} mm^3")
    if curved == 0:
        raise SystemExit("VACUOUS: no curved faces in the test solid")

    whole = occ.unify_same_domain(solid, unify_edges=False)
    n_whole = len(occ.faces(whole))
    print(f"whole-solid unify: {len(all_faces)} -> {n_whole} faces\n")

    verdicts: dict[str, bool] = {}
    with tempfile.TemporaryDirectory() as tmpdir:
        for n in args.tiles:
            tiles = merged_tiles(all_faces, n)
            n_faces = sum(len(t) for t in tiles)
            print(f"{len(tiles)} tiles, {n_faces} merged faces "
                  f"({100 * (n_faces - n_whole) / max(n_whole, 1):+.2f}% against whole):")
            for label, route, uses_files in ROUTES:
                sub = os.path.join(tmpdir, f"{len(tiles)}_{label.split()[0]}")
                os.makedirs(sub, exist_ok=True)
                faces = route(tiles, sub) if uses_files else route(tiles)
                fe = len(free_edges(faces))
                ok = fe == 0
                verdicts[label] = verdicts.get(label, True) and ok
                print(f"    {label:<26} free edges {fe:>7}  {'OK' if ok else 'BROKEN'}")
            print()

    print("=== VERDICT ===")
    for label, _route, _f in ROUTES:
        print(f"  {label:<26} {'PASS' if verdicts[label] else 'FAIL'}")
    if not verdicts["in-process (G14's route)"]:
        print("\nThe control failed: this harness does not reproduce G14 and none "
              "of the other rows can be read.")
    elif verdicts["one file per tile"]:
        print("\nSerialization preserves the seam. Phase 3 may dispatch tiles "
              "across the pool and reassemble with BRep_Builder.Add alone.")
    else:
        print("\nThe seam does not survive one file per tile. Phase 3 as planned "
              "cannot reassemble by shared topology across processes, and needs a "
              "seam repair whose cost must be measured before it is built.")


if __name__ == "__main__":
    main()
