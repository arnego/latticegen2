"""Gate G18 — can `validate` use more than one core, and at what risk?

`validate` (docs/algorithm.md §9) is 225 s of the `cc=5, t=1` rehearsal at 0.98
cores. docs/testing.md lumps it with `simplify` — "parallelising `simplify` and
`validate` is correct but was not a wall-clock win on this part" — which is true
of dispatching *per body*, where one dominant solid sets the floor. It says
nothing about going below the body, and the two stages differ there in a way
that matters:

    `simplify` must hand back geometry; `validate` hands back a scalar.

That is exactly the property G15 killed Phase 3 over. Unified tiles have to
reassemble by shared topology, and a `.brep` per tile duplicates every seam
edge. `_worker_validate` returns ``(bool, float)``. Nothing reassembles, so
duplicated edges across chunk files cannot hurt — G15 has nothing to attach to.
G17 (no thread speedup, OCP holds the GIL) still applies, so any Python-level
split has to be processes.

**But there is a cheaper option than any split, and this gate checks it first.**
`ShapeUpgrade_UnifySameDomain` "exposes no internal parallel mode — no
SetRunParallel, no thread-pool hook — so there is nothing to switch on either"
(docs/specification.md §11). `BRepCheck_Analyzer` *does*::

    BRepCheck_Analyzer(S, GeomControls=True, theIsParallel=False, theIsExact=False)

OCCT's own native threads are not Python threads, so the GIL result does not
bind them. If the flag works, `validate` needs no decomposition at all and the
verdict stays OCCT's own — no chunking, no reassembly, no new correctness
surface.

Three parts, in the order that decides the question:

* **A — verdict safety.** The parallel flag must not change a single verdict.
  Checked against the four *real* invalid faces committed from the rehearsal
  (``test/*.brep``, docs/algorithm.md §8 rungs 1 and 2) as well as valid solids.
  This is the negative control, and it is not optional: G10's lesson in this
  project is that a blind scanner reporting "no faults" costs real time, and
  `validate` is the exit-4 gate before a 2 GB file ships. A check that only ever
  passes on good geometry proves nothing.
* **B — what the flag buys**, in wall clock and core-equivalents.
* **C — whether a manual chunked split could buy more**, by timing the
  whole-solid check against the sum of its per-face parts and the structural
  (shell/solid) remainder. Only if C shows per-face work dominating *and* B
  falls well short of `W` is the riskier split worth designing.

Note on scale: the largest committed-scenario solid runs ~0.30 ms/face against
the rehearsal's ~0.39 ms/face (225 s / 584,028), so per-face cost here is
representative to within a small factor even though the solid is 40x smaller.
What does *not* transfer is any claim about the whole stage, which is why this
gate reports a split rather than a projected run time.

    python tools/prototypes/g18_validate_decomposition.py <output.step> [...]
"""

import argparse
import os
import time

import _bootstrap  # noqa: F401

from OCP.BRepCheck import BRepCheck_Analyzer, BRepCheck_NoError, BRepCheck_Shell
from OCP.TopoDS import TopoDS

from latticegen2 import occ
from latticegen2.parallel import read_brep

ROOT = _bootstrap.ROOT

# The real faces the rehearsal's two repair rungs were built against
# (docs/algorithm.md §8, G11 and G12). Every one of these is rejected by
# BRepCheck_Analyzer, which is what makes them usable as a negative control.
BAD_FACES = (
    "invalid-vertex-tolerance-ellipse.brep",
    "invalid-vertex-tolerance-bspline.brep",
    "self-intersecting-wire-cylinder.brep",
    "self-intersecting-wire-bspline.brep",
)


def timed(fn):
    t0 = time.perf_counter()
    t0c = time.process_time()
    out = fn()
    return out, time.perf_counter() - t0, time.process_time() - t0c


def part_a_verdict_safety() -> bool:
    """The parallel flag must agree with the serial one — including on faults."""
    print("=" * 72)
    print("A - verdict safety: does theIsParallel change any answer?")
    print("=" * 72)
    ok = True

    print(f"\n  {'fixture':<40} {'serial':>8} {'parallel':>9}  {'agree':>6}")
    for name in BAD_FACES:
        face = TopoDS.Face_s(read_brep(os.path.join(ROOT, "test", name)))
        serial = BRepCheck_Analyzer(face, True, False).IsValid()
        parallel = BRepCheck_Analyzer(face, True, True).IsValid()
        agree = serial == parallel
        ok &= agree and not serial
        print(f"  {name:<40} {str(serial):>8} {str(parallel):>9}  {str(agree):>6}")

    print(
        "\n  Every fixture must read False on BOTH sides. A True anywhere means the\n"
        "  control is broken, not that the geometry improved — these faces are the\n"
        "  ones the rehearsal's repair rungs exist for."
    )
    return ok


def structural_only(solid) -> bool:
    """The part of the check that is genuinely global, and cannot be chunked.

    Face/wire/edge/vertex checks are local to a face and would divide across
    workers. Shell closure and orientation consistency are properties of the
    whole shell — they are what a chunked scan would still have to run in one
    piece, so their cost sets the floor any split could reach.
    """
    for shell in occ.shells(solid):
        checker = BRepCheck_Shell(TopoDS.Shell_s(shell))
        if checker.Closed() != BRepCheck_NoError:
            return False
        if checker.Orientation() != BRepCheck_NoError:
            return False
    return True


def part_bc(path: str) -> None:
    print("\n" + "=" * 72)
    print(f"B/C - cost, on {os.path.basename(path)}")
    print("=" * 72)

    shape = occ.read_step(path) if path.lower().endswith((".step", ".stp")) else read_brep(path)
    solids = occ.solids(shape)

    for i, solid in enumerate(solids):
        n_faces, n_edges = occ.count_subshapes(solid)
        if n_faces < 500:
            continue  # a scrap; its timings are noise either way
        faces = occ.faces(solid)

        (v_ser, w_ser, c_ser) = timed(lambda: BRepCheck_Analyzer(solid, True, False).IsValid())
        (v_par, w_par, c_par) = timed(lambda: BRepCheck_Analyzer(solid, True, True).IsValid())
        (v_face, w_face, _) = timed(
            lambda: all(BRepCheck_Analyzer(f, True, False).IsValid() for f in faces)
        )
        (v_str, w_str, _) = timed(lambda: structural_only(solid))
        (_, w_vol, _) = timed(lambda: occ.volume(solid))

        print(f"\n  solid {i}: {n_faces} faces, {n_edges} edges")
        print(f"    {'whole solid, serial':<28} {w_ser:>8.3f} s  {c_ser / max(w_ser, 1e-9):>5.2f} cores  -> {v_ser}")
        print(f"    {'whole solid, theIsParallel':<28} {w_par:>8.3f} s  {c_par / max(w_par, 1e-9):>5.2f} cores  -> {v_par}")
        print(f"    {'sum of per-face checks':<28} {w_face:>8.3f} s  {'':>11}  -> {v_face}")
        print(f"    {'structural only (shell)':<28} {w_str:>8.3f} s  {'':>11}  -> {v_str}")
        print(f"    {'occ.volume':<28} {w_vol:>8.3f} s")
        print(f"    speedup from the flag       {w_ser / max(w_par, 1e-9):>8.2f}x")
        print(f"    per-face share of serial    {100 * w_face / max(w_ser, 1e-9):>8.1f} %"
              "   <- what a chunked split could divide")
        print(f"    per-face cost               {1000 * w_ser / n_faces:>8.3f} ms/face")

        if v_ser != v_par:
            print("    *** VERDICT DISAGREEMENT - the flag is not safe ***")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("solids", nargs="*", help="generated .step/.brep files to measure")
    args = ap.parse_args()

    safe = part_a_verdict_safety()
    for path in args.solids:
        part_bc(path)

    print("\n" + "=" * 72)
    print(f"A verdict safety: {'PASS' if safe else 'FAIL'}")
    print("=" * 72)


if __name__ == "__main__":
    main()
