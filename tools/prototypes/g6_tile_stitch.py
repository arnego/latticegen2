"""Gate G6 — does tiling the boundary sew actually save wall time?

docs/specification.md §10 "Tile the boundary sew" proposes splitting a
component's pieces into spatial tiles, sewing each tile on its own, then sewing
the tile results together — on the theory that G5a's measured `n^1.8` scaling
applied twice to a smaller `n` beats applying it once to the full piece count.

G5b is a reason to be suspicious of that theory rather than to assume it: it
found that ``BRepBuilderAPI_Sewing`` pays for total face count even when a shape
contributes **zero** free edges — adding one closed, unmatchable 194,400-face
shell to a 4,000-piece sew cost +640 s for nothing to merge. Round 2 of tiling
feeds sewing several already-mostly-closed tile shells whose *total* face count
equals the untiled input's — only the tile-boundary edges are genuinely free —
so if G5b's finding generalises, round 2 could cost close to what one monolithic
sew costs, and tiling would only add round 1's cost on top for nothing.

This gate measures it directly, on the same real trimmed boundary pieces G5 used
(the `cc=5, t=1` rehearsal's leftovers), rather than assuming either story::

    python tools/prototypes/g6_tile_stitch.py --pieces path/to/temp/<stamp>

Tiling here is by contiguous chunks of the loaded piece order, not by lattice
index — the saved ``.brep`` files carry geometry only, not the node each piece
belongs to, so :func:`latticegen2.weld._tile_pieces`'s actual bucketing cannot be
reproduced from this leftover data. `boundary._split_batches`'s own docstring
records that a worker's batch already has spatial locality ("a worker's
junctions share input-body regions"), and this loads pieces in that same batch
order, so contiguous chunks are a reasonable stand-in for real spatial tiles —
close enough to answer the question this gate asks, which is about the sewing
call's cost model, not about the exact partition.
"""

import argparse
import os
import sys
import time

import _bootstrap  # noqa: F401
import g5_stitch_scaling as g5

from latticegen2 import occ

SEW_TOL = g5.SEW_TOL


def sew_all(shapes, tolerance: float) -> tuple[list, float]:
    t0 = time.perf_counter()
    sewn = occ.sew(shapes, tolerance, cutting=False)
    elapsed = time.perf_counter() - t0
    return occ.faces(sewn), elapsed


def tiled_sew(shapes, n_tiles: int, tolerance: float):
    """Round 1 over ``n_tiles`` contiguous chunks, then round 2 over their results."""
    chunk = max(1, -(-len(shapes) // n_tiles))  # ceil division
    chunks = [shapes[i:i + chunk] for i in range(0, len(shapes), chunk)]
    round1 = 0.0
    tile_shells = []
    for c in chunks:
        faces, elapsed = sew_all(c, tolerance)
        round1 += elapsed
        tile_shells.append(occ.faces_shell(faces))
    _, round2 = sew_all(tile_shells, tolerance)
    return len(chunks), round1, round2


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pieces", default=g5.DEFAULT_PIECES,
                     help="kept temp/<stamp> folder holding boundary_*.brep")
    ap.add_argument("--n", type=int, default=4000, help="pieces to load and sew")
    ap.add_argument("--tiles", default="4,8,16", help="comma-separated tile counts to try")
    args = ap.parse_args()

    directory = os.path.abspath(args.pieces)
    print(f"G6: reading trimmed boundary pieces from {directory}")
    pieces = g5.load_pieces(directory, args.n)
    print(f"G6: {len(pieces)} pieces loaded\n")

    _, baseline = sew_all(pieces, SEW_TOL)
    print(f"baseline, one call, {len(pieces)} pieces: {baseline:8.2f} s\n")

    for spec in args.tiles.split(","):
        n_tiles = int(spec)
        got_tiles, round1, round2 = tiled_sew(pieces, n_tiles, SEW_TOL)
        total = round1 + round2
        verdict = "faster" if total < baseline else "SLOWER"
        print(
            f"tiles={got_tiles:>3}  round1={round1:8.2f}s  round2={round2:8.2f}s  "
            f"total={total:8.2f}s  vs baseline {baseline:8.2f}s -> {verdict} "
            f"({baseline / total:.2f}x)"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
