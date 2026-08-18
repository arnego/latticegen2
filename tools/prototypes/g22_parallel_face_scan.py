"""Gate G22 — can the sew's per-face validity scan be one parallel call?

`occ.fix_vertex_tolerances` (docs/algorithm.md §8) opens by asking, once per
boundary face, whether ``BRepCheck_Analyzer`` rejects it — and only the rejects
are examined further. On the `cc=5, t=1` rehearsal that scan is **44.6 s** over
301,505 faces, 5.8 % of `stitch` and ~1.2 % of the run, all of it serial on the
master (docs/specification.md §10).

docs/specification.md §10 rejects the two obvious ways to speed it up, and both
rejections are worth restating because this gate has to get past them:

* **Dispatching it to workers is not free I/O.** No file already on disk holds
  exactly the face list that needs scanning in the order it needs it, so
  dispatch means writing ~380 MB first, against a ~37 s saving.
* **Querying faces inside the assembled solid changes the predicate.** G18
  found ``BRepCheck_Analyzer`` has a ``theIsParallel`` flag, so one analyzer
  over everything would parallelise the scan with no I/O at all — but
  ``IsValid(subshape)`` on a solid is the **in-context** overload, and this code
  checks each face **standalone**. Those are different predicates; that
  difference is exactly what G12 used to diagnose repair rung 2, and this
  part's result (19 corrected, 0 residual) is the only evidence the repair is
  complete.

**What this gate tests is a third thing, and the distinction is the whole
point.** A ``TopoDS_Compound`` of loose faces carries no shell and no solid, so
a face inside it has no context to be checked in: its parent in the exploration
is the compound itself. If that is right, one analyzer over a compound of every
face computes exactly today's standalone predicate, face for face, while
``theIsParallel`` divides it across OCCT's own native threads — no dispatch, no
extra I/O, no predicate change.

"If that is right" is not something to assume. So:

* **A — predicate identity**, face for face, over a corpus that contains real
  faults. Per G10, a scanner that only ever agrees on sound geometry proves
  nothing, so the four real invalid faces committed from the rehearsal
  (``test/*.brep``, rungs 1 and 2 of docs/algorithm.md §8) are mixed into the
  corpus and must read invalid on **both** sides.
* **B — identity after repair.** ``fix_vertex_tolerances`` re-checks each face
  it touched to decide ``repaired`` against ``still_invalid``. A repaired face
  must read the same way in a compound as it does alone, or the counts the run
  log reports would change meaning.
* **C — what it buys**: wall clock and core-equivalents against the serial
  loop, with the flag on and off, so the speedup is attributed to the threads
  rather than to having made one call instead of N.
* **D — what it costs in memory.** The serial loop holds one analyzer at a
  time; a compound scan holds one result per subshape for the whole corpus at
  once, and this stage runs where the rehearsal is already carrying ~3 GB. So
  the scan is chunked in production, and this part measures the bytes per face
  that says how large a chunk may be.

Pass bar, set in advance: **zero** disagreements in A and B — this is a
correctness gate first — and a speedup in C worth the change. Anything else is
a disproof, and the scan stays as it is.

    python tools/prototypes/g22_parallel_face_scan.py [extra.step ...]

The default corpus is the committed 80 mm ball golden sample, which is a real
lattice output whose faces share edges exactly as the sewn boundary layer's do,
so the gate is reproducible from a clean checkout. Extra STEP files are added
to the corpus when given — a real rehearsal output is the corpus this is
ultimately about.
"""

import argparse
import gc
import os
import time

import _bootstrap  # noqa: F401

from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.TopoDS import TopoDS

from latticegen2 import occ
from latticegen2.parallel import read_brep

ROOT = _bootstrap.ROOT

# The real faults, committed from the `cc=5, t=1` rehearsal: two vertices
# recorded off their edge's own curve (rung 1, G11) and two falsely
# self-intersecting wires (rung 2, G12).
BAD_FACES = (
    "invalid-vertex-tolerance-ellipse.brep",
    "invalid-vertex-tolerance-bspline.brep",
    "self-intersecting-wire-cylinder.brep",
    "self-intersecting-wire-bspline.brep",
)

GOLDEN = os.path.join(ROOT, "test", "80mm-test-ball-cc20t4-golden-sample.step")


def load_bad() -> list:
    return [
        TopoDS.Face_s(read_brep(os.path.join(ROOT, "test", name)))
        for name in BAD_FACES
    ]


def serial_scan(faces) -> list:
    """Today's predicate, exactly as `fix_vertex_tolerances` computes it."""
    return [BRepCheck_Analyzer(f).IsValid() for f in faces]


def compound_scan(faces, parallel: bool) -> list:
    """One analyzer over a compound of loose faces, queried face by face."""
    analyzer = BRepCheck_Analyzer(occ.compound(faces), True, parallel)
    return [analyzer.IsValid(f) for f in faces]


def chunked_scan(faces, chunk: int, parallel: bool = True) -> list:
    """What production does: the same scan, over bounded slices.

    Chunking is what keeps part D's per-face memory bounded, so it has to be
    timed rather than assumed neutral — a chunk small enough to starve OCCT's
    threads would give the memory back and the speedup with it.
    """
    out = []
    for i in range(0, len(faces), chunk):
        window = faces[i:i + chunk]
        analyzer = BRepCheck_Analyzer(occ.compound(window), True, parallel)
        out.extend(analyzer.IsValid(f) for f in window)
    return out


def timed(fn):
    t0, t0c = time.perf_counter(), time.process_time()
    out = fn()
    return out, time.perf_counter() - t0, time.process_time() - t0c


def report_disagreements(label, faces, a, b, bad_index) -> bool:
    wrong = [i for i in range(len(faces)) if a[i] != b[i]]
    print(f"  {label:<44} {len(wrong)} disagreement(s) over {len(faces)} face(s)")
    for i in wrong[:10]:
        origin = "committed fault" if i in bad_index else "corpus"
        print(f"      face {i} ({origin}): standalone={a[i]} compound={b[i]}")
    return not wrong


def part_a(faces, bad_index) -> bool:
    print("=" * 72)
    print("A - predicate identity: standalone face vs face inside a compound")
    print("=" * 72)
    standalone = serial_scan(faces)
    ok = True
    for parallel in (False, True):
        got = compound_scan(faces, parallel)
        ok &= report_disagreements(
            f"compound scan, theIsParallel={parallel}", faces, standalone, got, bad_index
        )
    # Chunking is what production actually runs (part D), and a chunk boundary
    # falling between two faces that share an edge is exactly where a context
    # difference would show up if there were one.
    ok &= report_disagreements(
        "chunked scan at 1000", faces, standalone, chunked_scan(faces, 1000), bad_index
    )
    n_bad = sum(1 for i in bad_index if not standalone[i])
    print(f"\n  {n_bad} of {len(bad_index)} committed fault(s) read invalid standalone.")
    if n_bad != len(bad_index):
        print("  CONTROL BROKEN: a fault fixture loaded as valid. Nothing below is evidence.")
        return False
    par = compound_scan(faces, True)
    seen = sum(1 for i in bad_index if not par[i])
    print(f"  {seen} of {len(bad_index)} committed fault(s) read invalid in the compound.")
    print(
        "  Both counts must be 4. A scan that misses a real fault is worse than a\n"
        "  slow one: these are the faces the repair exists for (G10's lesson)."
    )
    return ok and seen == len(bad_index)


def part_b(faces) -> bool:
    print()
    print("=" * 72)
    print("B - identity after repair: the counts the run log reports")
    print("=" * 72)
    repaired, still_invalid = occ.fix_vertex_tolerances(faces)
    standalone = serial_scan(faces)
    compound = compound_scan(faces, True)
    print(f"  fix_vertex_tolerances -> repaired={repaired} still_invalid={still_invalid}")
    ok = report_disagreements("after repair", faces, standalone, compound, set())
    print(
        "  The repair decides `repaired` vs `still_invalid` by re-checking each face\n"
        "  it touched, so a post-repair disagreement would silently change what the\n"
        "  run log means, not merely how long the scan takes."
    )
    return ok


def part_c(faces) -> None:
    print()
    print("=" * 72)
    print("C - what the compound scan buys")
    print("=" * 72)
    print(f"  corpus: {len(faces)} faces\n")
    print(f"  {'scan':<40} {'wall':>9} {'cpu':>9} {'cores':>7} {'speedup':>8}")
    _, w0, c0 = timed(lambda: serial_scan(faces))
    print(f"  {'serial, one analyzer per face':<40} {w0:>8.3f}s {c0:>8.3f}s "
          f"{c0 / w0:>7.2f} {1.0:>8.2f}")
    for parallel in (False, True):
        _, w, c = timed(lambda: compound_scan(faces, parallel))
        label = f"compound, theIsParallel={parallel}"
        print(f"  {label:<40} {w:>8.3f}s {c:>8.3f}s {c / w:>7.2f} {w0 / w:>8.2f}")
    for chunk in (1000, 5000, 20000):
        _, w, c = timed(lambda: chunked_scan(faces, chunk))
        label = f"chunked at {chunk}, theIsParallel=True"
        print(f"  {label:<40} {w:>8.3f}s {c:>8.3f}s {c / w:>7.2f} {w0 / w:>8.2f}")
    print(
        "\n  Two rows for the compound so the speedup is attributed: the "
        "theIsParallel=False\n  row is what one call instead of N buys on its own, "
        "and the difference between\n  the two rows is what OCCT's threads buy."
    )


def part_d(faces) -> None:
    print()
    print("=" * 72)
    print("D - what the compound scan costs in memory")
    print("=" * 72)
    try:
        import psutil
    except ImportError:
        print("  psutil not installed; skipped (python -m pip install psutil)")
        return
    proc = psutil.Process()
    gc.collect()
    before = proc.memory_info().rss
    analyzer = BRepCheck_Analyzer(occ.compound(faces), True, True)
    held = [analyzer.IsValid(f) for f in faces]
    peak = proc.memory_info().rss
    del analyzer, held
    gc.collect()
    after = proc.memory_info().rss
    grew = peak - before
    print(f"  RSS before {before / 1e6:.0f} MB, holding the analyzer {peak / 1e6:.0f} MB, "
          f"released {after / 1e6:.0f} MB")
    print(f"  {grew / max(len(faces), 1):.0f} bytes per face while the analyzer is alive")
    print(
        f"\n  At the rehearsal's 301,505 boundary faces that extrapolates to "
        f"~{grew / max(len(faces), 1) * 301505 / 1e9:.2f} GB\n"
        "  held at once, against a serial loop's one analyzer. That is the reason the\n"
        "  production scan is chunked rather than one call over everything: `stitch`\n"
        "  already peaks near 3 GB, and this stage must not be where a run runs out."
    )


def corpus(paths) -> tuple:
    faces = []
    for path in paths:
        got = occ.faces(occ.read_step(path))
        print(f"  {os.path.basename(path)}: {len(got)} faces")
        faces.extend(got)
    bad = load_bad()
    bad_index = set(range(len(faces), len(faces) + len(bad)))
    faces.extend(bad)
    print(f"  {len(bad)} committed fault fixture(s); corpus is {len(faces)} faces\n")
    return faces, bad_index


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("steps", nargs="*", help="extra STEP files to add to the corpus")
    args = ap.parse_args()
    occ.quiet_kernel()

    print("Corpus")
    print("-" * 72)
    faces, bad_index = corpus([GOLDEN, *args.steps])

    ok = part_a(faces, bad_index)
    ok &= part_b(load_bad())
    part_c(faces)
    part_d(faces)

    print()
    print("=" * 72)
    print("PASS - the compound scan is the same predicate" if ok else
          "FAIL - the compound scan is NOT the same predicate; leave the scan alone")
    print("=" * 72)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
