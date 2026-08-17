"""Gate G16 — is same-domain unification priced by the *size* of its input?

docs/specification.md §11 records two attempts to make `simplify` cheaper by
giving `ShapeUpgrade_UnifySameDomain` less to look at, and both under-delivered
badly against projections built on face counts. Phase 2 cut its input ~31 % and
the stage fell 12.4 %. The restricted face merge cut it **46 %** and the stage
did not fall at all, while producing a byte-identical output.

The natural conclusion — "the call is not priced by its input" — is **wrong**,
and this gate is what disproves it. On a generic subset the cost is very nearly
proportional to the face count:

    elasticity = (relative time saved) / (relative input removed)

and this measures **~1.0**. Remove a representative 40 % of the faces and you
save very nearly 40 % of the time.

    python tools/prototypes/g16_unify_elasticity.py
    python tools/prototypes/g16_unify_elasticity.py --m 8 --fractions 0.9 0.75 0.5

**That is what makes the pipeline's own measurement interesting rather than
merely disappointing.** On a real `dense-lattice` solid, removing 20 % of the
faces saved only 6 % of the time — an elasticity near 0.3, against the ~1.0
here. The faces removed were not a representative sample: a *correct*
restriction skips exactly the faces unification would return unchanged, which
are exactly the cheap ones, and keeps exactly the faces that merge, which are
the expensive ones. So the restriction is self-defeating by construction, and
it is the **selection**, not the kernel's pricing, that destroys the saving.

Two things this gate deliberately does not do. It does not measure the real
elasticity of a correct restriction — that needs geometry which actually merges
heavily, and the trimmed test solid here merges barely at all now that Phase 2
builds the interior pre-merged (docs/algorithm.md §6); the number to cite for
that is the pipeline measurement above, taken on `dense-lattice`. And it does
not time the edge pass, which no restriction may touch anyway (§9).

**Run this before believing any proposal to feed this stage less** — but read
its ~1.0 as the *ceiling* a perfectly representative removal would reach, never
as what a correct restriction will get.
"""

import argparse
import time

import _bootstrap  # noqa: F401
import numpy as np

import g13_unify_scaling as g13
import g14_tiled_unify_trimmed as g14
from latticegen2 import occ

FRACTIONS = [0.9, 0.8, 0.6, 0.4]


def timed_unify(faces) -> tuple[float, int]:
    """Face merge only — the pass a restriction could ever attack (§9)."""
    shell = occ.faces_shell(faces)
    t0 = time.perf_counter()
    merged = occ.unify_same_domain(shell, unify_edges=False)
    return time.perf_counter() - t0, len(occ.faces(merged))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--m", type=int, default=8, help="grid edge in cells")
    ap.add_argument("--radius", type=float, default=0.75,
                    help="sphere radius as a fraction of the smallest half-extent")
    ap.add_argument("--fractions", type=float, nargs="+", default=FRACTIONS,
                    help="input fractions to time, as a share of all faces")
    ap.add_argument("--repeats", type=int, default=3,
                    help="timings per point; the fastest is kept, since a slow "
                         "run is noise and a fast one cannot be")
    args = ap.parse_args()

    print("G16 — unification's cost against the size of its input\n")
    solid = g14.trimmed_lattice_solid(args.m, args.radius)
    all_faces = occ.faces(solid)
    curved = sum(1 for f in all_faces if not occ.is_planar(f))
    print(f"input: {len(all_faces)} faces ({curved} curved)")

    # Centroid order, so a fraction is a spatially coherent slab rather than a
    # scatter of isolated faces with no neighbours left to merge with.
    order = np.argsort(g13.face_centroids(all_faces)[:, 0])
    ordered = [all_faces[i] for i in order]

    def best(faces):
        runs = [timed_unify(faces) for _ in range(args.repeats)]
        return min(t for t, _n in runs), runs[0][1]

    t_full, n_full = best(ordered)
    print(f"whole-solid face merge: {len(ordered)} -> {n_full} faces "
          f"({len(ordered) - n_full} merged away)\n")

    print(f"{'input':>8} {'share':>7} {'time s':>9} {'saved':>8} {'faces out':>10} "
          f"{'elasticity':>11}")
    print(f"{len(ordered):>8} {1.0:>7.2f} {t_full:>9.3f} {'—':>8} {n_full:>10} {'—':>11}")

    elasticities = []
    for frac in sorted(args.fractions, reverse=True):
        n = max(1, int(round(frac * len(ordered))))
        t, out = best(ordered[:n])
        removed = 1.0 - n / len(ordered)
        saved = 1.0 - t / t_full
        e = saved / removed
        elasticities.append(e)
        print(f"{n:>8} {n / len(ordered):>7.2f} {t:>9.3f} {100 * saved:>7.1f}% "
              f"{out:>10} {e:>11.2f}")

    mean = float(np.mean(elasticities))
    print(f"\nmean elasticity on a generic subset: {mean:.2f}")
    if mean >= 0.8:
        print("\nThe kernel prices its input close to linearly. So a restriction "
              "that\nsaves far less than it removes — the pipeline measured 6 % of "
              "the time\nfor 20 % of the faces — is not being let down by the "
              "kernel: it is\nremoving the cheap faces and keeping the expensive "
              "ones, which is what\nany *correct* restriction must do "
              "(docs/specification.md §11).")
    else:
        print("\nThe kernel does not price its input linearly here, so a "
              "restriction cannot\nrecover its own bookkeeping even on a "
              "representative subset.")


if __name__ == "__main__":
    main()
