"""Gate G8 — can boundary-sew round 2 skip faces that are already fully sewn?

docs/specification.md §10 ranks a cheaper round 2 (optimization path 3) after
rejecting a hierarchical tree reduction on paper (docs/algorithm.md §8): round
2's cost tracks total face count almost flatly, not shape count, so a tree
would pay that flat cost once per level for nothing. The lever that survives
that argument is different — stop handing round 2 faces it cannot do anything
with.

After round 1, every tile's own result is itself a sewn shell: each of its
edges is either used **twice** already (fully joined to a neighbour within the
same tile — nothing left for round 2 to do with the face that owns it) or used
**once** (a genuine free edge — either a tile-to-tile seam round 2 exists to
close, or an interface hole meant to stay open for the interior shell, and
round 2 cannot tell those two apart in advance any more than a full round 2
call can). So the only faces round 2 can possibly affect are the ones bearing
at least one free edge; everything else can be carried into the final shell
unchanged, by direct reference, with no sewing call touching it at all.

That argument only holds if this reflects reality correctly, so it is
prototyped and checked identity-for-identity against a real (small) tiled sew,
not just argued:

    python tools/prototypes/g8_seam_only_sew.py

Pass bar: seam-only round 2 (sew only the free-edge-bearing faces from every
tile; concatenate the rest straight in) produces **exactly** what today's
round 2 (sew every tile's full result together) produces — same total face
count, same free-edge count (zero, for a chain closed at both ends), and the
same volume to machine precision. Anything short of that is a disproof, and
this is not adopted into :mod:`latticegen2.weld` at all if it comes out that
way — see the module docstring's honesty rule (docs/algorithm.md §11):
optimizations are only allowed to fail by doing more work, never by producing
a different result.
"""

import time

import _bootstrap  # noqa: F401
import numpy as np
from OCP.TopAbs import TopAbs_ShapeEnum

from latticegen2 import occ, weld
from latticegen2.boundary import BoundaryPiece, trim_junction
from latticegen2.junction import build_template
from latticegen2.lattice import lattice_params, nodes

CC, T = 10.0, 1.5
N_NODES = 40
TILE_TARGET = 5


def _long_box():
    """Wide enough along +e0 to hold the whole chain untrimmed (mirrors test_weld.py)."""
    return occ.prism(
        occ.polygon_face(np.array([[-400.0, -400.0, -400.0], [400.0, -400.0, -400.0],
                                   [400.0, 400.0, -400.0], [-400.0, 400.0, -400.0]])),
        np.array([0.0, 0.0, 800.0]),
    )


def _line_pieces(lp, tpl, n: int) -> list[BoundaryPiece]:
    """``n`` trimmed junctions in a row along +e0, closed at both ends.

    Identical construction to test_weld.py's own ``_line_pieces`` — a tube
    capped at both ends, so a correct round 2 (of any kind) leaves the whole
    thing fully closed and its volume is exactly ``n * volume(J)``.
    """
    box = _long_box()
    positions = nodes(lp, np.array([(i, 0, 0) for i in range(n)], dtype=np.int64))
    pieces = []
    for i in range(n):
        faces, tags, vol = trim_junction(lp, tpl, positions[i], box)[0]
        keep = [
            f for f, t in zip(faces, tags)
            if not (t == 3 and i > 0) and not (t == 0 and i < n - 1)
        ]
        pieces.append(BoundaryPiece(node=(i, 0, 0), volume=vol, faces=keep))
    return pieces


def split_by_free_edges(faces: list):
    """``(seam, interior)``: faces touching >=1 free edge of this face set, and the rest.

    O(faces x edges x free_edges) shape comparisons — fine for this gate's
    scale, not written for production use.
    """
    fe = weld.free_edges(faces)
    seam, interior = [], []
    for f in faces:
        edges = list(occ._explore(f, TopAbs_ShapeEnum.TopAbs_EDGE))
        if any(e.IsSame(x) for e in edges for x in fe):
            seam.append(f)
        else:
            interior.append(f)
    return seam, interior


def main() -> None:
    lp = lattice_params(CC, T)
    tpl = build_template(lp)
    pieces = _line_pieces(lp, tpl, N_NODES)
    print(f"{N_NODES} pieces built.")

    plan = weld._tile_pieces(pieces, target=TILE_TARGET, min_to_tile=1)
    assert plan is not None and len(plan) > 1, "must actually tile for this gate to mean anything"
    print(f"tiled into {len(plan)} tiles")

    t0 = time.perf_counter()
    tile_results = [weld._sew_faces([p.faces for p in tile], weld.SEW_TOLERANCE) for tile in plan]
    t_round1 = time.perf_counter() - t0
    total_faces_in = sum(len(t) for t in tile_results)
    print(f"round 1: {t_round1:.3f}s, {total_faces_in} faces across {len(tile_results)} tile results")

    # --- baseline: today's round 2, sew every tile's full result together ---
    t0 = time.perf_counter()
    baseline_faces = weld._sew_faces(tile_results, weld.SEW_TOLERANCE)
    t_round2_full = time.perf_counter() - t0
    baseline_free = weld.free_edges(baseline_faces)
    baseline_solid = occ.make_solid(occ.faces_shell(baseline_faces))
    baseline_volume = occ.volume(baseline_solid)

    # --- candidate: seam-only round 2 ---
    seam_lists, interior_lists = [], []
    for tr in tile_results:
        seam, interior = split_by_free_edges(tr)
        seam_lists.append(seam)
        interior_lists.append(interior)
    n_seam = sum(len(s) for s in seam_lists)
    n_interior = sum(len(i) for i in interior_lists)
    print(f"seam faces: {n_seam} ({100 * n_seam / total_faces_in:.1f}%), "
          f"interior (carried through unsewn): {n_interior}")

    t0 = time.perf_counter()
    sewn_seam = weld._sew_faces(seam_lists, weld.SEW_TOLERANCE)
    seam_only_faces = sewn_seam + [f for lst in interior_lists for f in lst]
    t_round2_seam = time.perf_counter() - t0
    seam_only_free = weld.free_edges(seam_only_faces)
    seam_only_solid = occ.make_solid(occ.faces_shell(seam_only_faces))
    seam_only_volume = occ.volume(seam_only_solid)

    print(f"\nround 2, full:       {t_round2_full:.3f}s, {len(baseline_faces)} faces, "
          f"{len(baseline_free)} free edges, volume {baseline_volume:.9f} mm^3")
    print(f"round 2, seam-only:  {t_round2_seam:.3f}s, {len(seam_only_faces)} faces, "
          f"{len(seam_only_free)} free edges, volume {seam_only_volume:.9f} mm^3")

    same_face_count = len(baseline_faces) == len(seam_only_faces)
    same_free_edges = len(baseline_free) == len(seam_only_free) == 0
    same_volume = abs(baseline_volume - seam_only_volume) <= 1e-9 * abs(baseline_volume)
    passed = same_face_count and same_free_edges and same_volume

    print("\n=== VERDICT ===")
    print(f"face count matches:  {same_face_count}")
    print(f"both fully closed:   {same_free_edges}")
    print(f"volume matches:      {same_volume} (diff {abs(baseline_volume - seam_only_volume):.3e})")
    print("PASS - seam-only round 2 is identity-safe, adopt it" if passed else
          "FAIL - seam-only round 2 diverges from a full round 2, do not adopt it")


if __name__ == "__main__":
    main()
