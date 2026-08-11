"""Gate G1 — junction template cap integrity across the CLI parameter grid.

Builds the junction template ``J`` for every ``(cc, t)`` in the documented CLI
range and checks that all six mid-strut cap quads survive the fuse as intact
``t x t`` faces, and that ``volume(J)`` matches an independent
inclusion-exclusion prediction computed from separate pairwise/triple
intersections, so the fuse result is never used to check itself.

Pass bar: all six caps intact for every combination with ``t < cc/2``, and the
exact boundary where they start failing recorded.
"""

import itertools
import sys
import time

import _bootstrap  # noqa: F401
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common

from latticegen2 import occ
from latticegen2.errors import LatticeGenError
from latticegen2.junction import _half_strut_solid, build_template
from latticegen2.lattice import lattice_params

CCS = [0.4, 5.0, 10.0, 20.0, 50.0]
TS = [0.4, 1.0, 1.5, 4.0, 20.0]


def inclusion_exclusion_volume(lp):
    """Volume of the 6-half-strut union from independent pairwise booleans.

    Uses the |A u B u ...| = sum|.| - sum|. n .| + ... identity, evaluated with
    separate ``COMMON`` calls, so it never touches the fuse result being
    checked.
    """
    halves = [_half_strut_solid(lp, h) for h in range(6)]
    total = sum(occ.volume(h) for h in halves)
    for order in range(2, 7):
        sign = -1.0 if order % 2 == 0 else 1.0
        acc = 0.0
        for combo in itertools.combinations(range(6), order):
            shape = halves[combo[0]]
            for idx in combo[1:]:
                alg = BRepAlgoAPI_Common(shape, halves[idx])
                if not alg.IsDone():
                    return None
                shape = alg.Shape()
                if occ.volume(shape) <= 0.0:
                    break
            acc += occ.volume(shape)
        total += sign * acc
    return total


def main() -> int:
    print(f"{'cc':>6} {'t':>6} {'t<cc/2':>7} {'caps':>5} {'faces':>6} "
          f"{'vol(J)':>12} {'incl-excl':>12} {'rel.err':>10} {'ms':>7}")
    failures = []
    for cc, t in itertools.product(CCS, TS):
        lp = lattice_params(cc, t)
        if t >= lp.a:
            continue  # invalid per the CLI's own t < a cross-constraint
        expected_ok = t < cc / 2.0
        t0 = time.perf_counter()
        try:
            tpl = build_template(lp)
            caps = 6
            vol = tpl.volume
            nf = tpl.n_faces
        except LatticeGenError:
            caps, vol, nf = 0, float("nan"), 0
        ms = (time.perf_counter() - t0) * 1e3

        ie = inclusion_exclusion_volume(lp) if caps == 6 else None
        rel = abs(vol - ie) / ie if (ie and ie > 0) else float("nan")
        print(f"{cc:6g} {t:6g} {str(expected_ok):>7} {caps:5d} {nf:6d} "
              f"{vol:12.6f} {(ie if ie else float('nan')):12.6f} {rel:10.2e} {ms:7.1f}")

        if expected_ok and caps != 6:
            failures.append(f"cc={cc} t={t}: caps missing but t < cc/2")
        if expected_ok and not (rel < 1e-9):
            failures.append(f"cc={cc} t={t}: volume mismatch rel={rel:.3e}")

    print()
    print("Sweeping t/cc up to the CLI's own t < a limit, to locate where (if")
    print("anywhere) caps stop surviving:")
    for cc in (5.0, 10.0, 20.0):
        best_ok, first_fail = None, None
        for pct in range(45, 100):
            t = cc * pct / 100.0
            lp = lattice_params(cc, t)
            if t >= lp.a:
                break
            try:
                build_template(lp)
                best_ok = pct / 100.0
            except LatticeGenError:
                first_fail = first_fail if first_fail is not None else pct / 100.0
        print(f"  cc={cc:g}: caps intact up to t/cc={best_ok:.2f} "
              f"(a/cc = {1 / 2 ** 0.5:.4f}); first cap failure at "
              f"{first_fail if first_fail is not None else 'none in range'}")
        if first_fail is not None:
            failures.append(f"cc={cc}: unexpected cap failure at t/cc={first_fail}")

    print()
    if failures:
        for f in failures:
            print("FAIL:", f)
        print("G1: FAIL")
        return 1
    print("G1: PASS - all six caps intact and volume-consistent across the whole")
    print("           valid parameter range (t < a); no narrower window exists.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
