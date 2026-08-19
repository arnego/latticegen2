"""End-to-end verification harness (specification.md §6, docs/testing.md).

Runs every scenario in specification.md §6.1 as a subprocess and applies every
applicable check from §6.2. Each scenario is independent: one failing scenario
never prevents a later one from running, and the process exit code reflects the
run as a whole.

Usage:
    python tools/e2e.py [--only smoke-fast,dense-lattice]
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, HERE)

import verify_geometry as vg  # noqa: E402

from latticegen2 import occ, progress  # noqa: E402
from latticegen2.runlog import STAGES  # noqa: E402

PYTHON = os.environ.get("LATTICEGEN2_PYTHON", sys.executable)
MAIN = os.path.join(ROOT, "src", "main.py")
TESTDIR = os.path.join(ROOT, "test")
BALL = os.path.join(TESTDIR, "80mm-test-ball.step")
CYLINDER = os.path.join(TESTDIR, "test-cylinder.STEP")
BALL_GOLDEN = os.path.join(TESTDIR, "80mm-test-ball-cc20t4-golden-sample.step")
CYLINDER_GOLDEN = os.path.join(TESTDIR, "test-cylinder-cc10t1.5-golden-sample.step")


class Report:
    def __init__(self, name: str):
        self.name = name
        self.checks: list[tuple[str, bool, str]] = []

    def check(self, label: str, ok: bool, detail: str = "") -> bool:
        self.checks.append((label, bool(ok), detail))
        mark = "PASS" if ok else "FAIL"
        print(f"  [{mark}] {label}" + (f" — {detail}" if detail else ""), flush=True)
        return bool(ok)

    @property
    def ok(self) -> bool:
        return all(c[1] for c in self.checks)


def run_generator(args: list[str]):
    t0 = time.time()
    proc = subprocess.run(
        [PYTHON, MAIN, *args], capture_output=True, text=True, cwd=ROOT
    )
    return proc, time.time() - t0


def geometry_checks(rep: Report, output: str, input_path: str, cc: float, t: float) -> None:
    """The §6.2 checks that apply to any successfully generated STEP file."""
    rep.check("STEP file written and non-empty",
              os.path.isfile(output) and os.path.getsize(output) > 0,
              f"{os.path.getsize(output) / 1e6:.1f} MB" if os.path.isfile(output) else "missing")
    if not (os.path.isfile(output) and os.path.getsize(output) > 0):
        return

    valid, n_solids, invalid = vg.brepcheck(output)
    rep.check("round-trip read succeeds and yields solids", n_solids > 0, f"{n_solids} solid(s)")
    rep.check("exact B-rep validity (BRepCheck_Analyzer)", valid,
              "all valid" if valid else f"invalid solids: {invalid}")

    outside_check(rep, output, input_path, t)
    rep.check("bounding box within input + (cc+t)",
              vg.bounding_box_within(output, input_path, cc + t))

    mesh = vg.mesh_of(output, cc, t)
    manifold_ok, bad_edges = vg.manifold_check(mesh)
    rep.check("mesh is a closed manifold", manifold_ok,
              f"{len(mesh.tris)} triangles" + ("" if manifold_ok else f", {bad_edges} bad edges"))
    no_cross, pairs = vg.self_intersection_check(mesh)
    rep.check("no self-intersections", no_cross,
              "0 crossing pairs" if no_cross else f"{len(pairs)}+ crossing pairs")


def outside_check(rep: Report, output: str, input_path: str, t: float) -> None:
    """§6.2's "no generated material lies outside the input body", exactly where
    that is available and by a labelled weaker check where it is not.

    The exact figure is a boolean cut per output solid. On a large part the cut
    can mis-classify — it returned the whole solid as outside on a 43,530-face
    lattice, having genuinely run rather than declined — so `vg.material_outside`
    contradicts any non-zero remainder against a boolean-free containment check
    and reports the solid as unmeasured when the two disagree. Reporting that as
    a pass would be exactly the failure this check exists to catch, so it is
    reported as its own named check with the containment evidence attached.
    """
    result = vg.material_outside(output, input_path)
    if result["exact"]:
        rep.check("no material outside the input body",
                  abs(result["volume_mm3"]) < (t ** 3) * 1e-6,
                  f"{result['volume_mm3']:.6g} mm^3 over {result['solids']} solid(s)")
        return

    measured = sum(s["cut_mm3"] for s in result["per_solid"] if s["trusted"])
    n_ok = sum(1 for s in result["per_solid"] if s["trusted"])
    rep.check("no material outside the input body (solids the cut could measure)",
              abs(measured) < (t ** 3) * 1e-6,
              f"{measured:.6g} mm^3 over {n_ok} of {result['solids']} solid(s)")
    for item in result["unmeasured"]:
        ev = item["containment"]
        print(f"  [NOTE] solid {item['solid']} ({item['faces']} faces): the cut claimed "
              f"{item['cut_claimed_mm3']:.6g} mm^3 of its own {item['volume_mm3']:.6g} mm^3 "
              f"lies outside, which a boolean-free containment check contradicts; "
              f"falling back to that check")
        rep.check(f"solid {item['solid']} is contained in the input "
                  f"(weaker than the exact cut)",
                  bool(ev.get("contained")),
                  f"{ev['outside']} of {ev['sampled']} surface points outside at "
                  f"{vg.CONTAINMENT_SWEEP[0]:g} mm, none at "
                  f"{ev['cleared_at_mm']:g} mm (bar {vg.CONTAINMENT_TOL:g} mm)"
                  if ev.get("cleared_at_mm") is not None else
                  f"{ev['outside']} of {ev['sampled']} surface points still outside at "
                  f"{vg.CONTAINMENT_TOL:g} mm")


GOLDEN_EXACT_BUDGET_S = 1800.0
"""How long the exact boolean comparison gets before falling back.

Not a correctness knob, and deliberately generous. The exact test is a general
boolean between two complete lattices, so its cost grows far faster than the
lattice does: seconds on the ~10 k-face smoke output, **306 s** on the
~100 k-face dense one. An earlier 300 s budget missed that by six seconds and
sent a scenario that passes exactly (0 mm³) down the weaker sampled path — the
budget exists to stop a runaway, not to decide the common case. The fallback
never converts an unmeasured result into a pass.
"""


def golden_check(rep: Report, output: str, golden: str, cc: float, t: float) -> None:
    if not os.path.isfile(golden):
        print(f"  [SKIP] golden-sample comparison — {os.path.basename(golden)} not present")
        return
    tol = t ** 3
    diff = vg.golden_sample_volume_diff_bounded(output, golden, GOLDEN_EXACT_BUDGET_S)
    if diff is not None:
        rep.check("matches golden sample (exact symmetric-difference volume)", diff < tol,
                  f"{diff:.6g} mm^3 vs tolerance {tol:g} mm^3")
        return

    # The exact test did not finish. Report a quantified weaker check rather than
    # nothing — and label it, so a reader never mistakes it for the exact one.
    print(f"  [NOTE] exact golden comparison exceeded {GOLDEN_EXACT_BUDGET_S:.0f}s; "
          f"falling back to a sampled equivalence check")
    agree = vg.golden_sample_agreement(output, golden, cc, t)
    rep.check("golden sample: total volume matches", agree["volume_diff"] < tol,
              f"candidate {agree['candidate_volume']:.4f} vs golden "
              f"{agree['golden_volume']:.4f} mm^3 (diff {agree['volume_diff']:.6g})")
    rep.check("golden sample: bounding boxes match", agree["bbox_match"],
              f"worst extent difference {agree['bbox_delta']:.3g} mm "
              f"(tolerance {vg.BBOX_TOL} mm)")
    rep.check("golden sample: sampled membership agrees (weaker than the exact test)",
              agree["real_disagreements"] == 0,
              f"{agree['real_disagreements']} substantive of {agree['disagreements']} raw "
              f"disagreements over {agree['samples']} points; one sample = "
              f"{agree['resolution_mm3']:.4g} mm^3")


def scenario_smoke_fast(outdir: str) -> Report:
    rep = Report("smoke-fast")
    print(f"=== {rep.name} ===", flush=True)
    out = os.path.join(outdir, "smoke-fast.step")
    proc, elapsed = run_generator(
        ["-i", BALL, "-cc", "20", "-t", "4", "-o", out, "--cores", "4"]
    )
    rep.check("process exits 0", proc.returncode == 0, proc.stderr.strip()[-300:])
    rep.check("runtime under 10 minutes", elapsed < 600, f"{elapsed:.1f}s")
    if proc.returncode == 0:
        geometry_checks(rep, out, BALL, 20.0, 4.0)
    return rep


def scenario_smoke_verified(outdir: str) -> Report:
    rep = Report("smoke-verified")
    print(f"=== {rep.name} ===", flush=True)
    out = os.path.join(outdir, "smoke-verified.step")
    proc, elapsed = run_generator(
        ["-i", BALL, "-cc", "20", "-t", "4", "-o", out, "--cores", "4"]
    )
    rep.check("process exits 0", proc.returncode == 0, proc.stderr.strip()[-300:])
    rep.check("runtime under 20 minutes", elapsed < 1200, f"{elapsed:.1f}s")
    if proc.returncode == 0:
        geometry_checks(rep, out, BALL, 20.0, 4.0)
        golden_check(rep, out, BALL_GOLDEN, 20.0, 4.0)
    return rep


def scenario_dense_lattice(outdir: str) -> Report:
    rep = Report("dense-lattice")
    print(f"=== {rep.name} ===", flush=True)
    out = os.path.join(outdir, "dense-lattice.step")
    proc, elapsed = run_generator(
        ["-i", CYLINDER, "-cc", "10", "-t", "1.5", "-o", out, "--cores", "6"]
    )
    rep.check("process exits 0", proc.returncode == 0, proc.stderr.strip()[-300:])
    rep.check("runtime under 10 minutes", elapsed < 600, f"{elapsed:.1f}s")
    if proc.returncode == 0:
        geometry_checks(rep, out, CYLINDER, 10.0, 1.5)
        golden_check(rep, out, CYLINDER_GOLDEN, 10.0, 1.5)
    return rep


def _mask_volatile(text: str) -> str:
    """Blank the two timestamps a STEP header carries, and nothing else."""
    text = re.sub(r"'\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d'", "'<ts>'", text)
    return re.sub(r"generated=\d{4}-\d\d-\d\d \d\d:\d\d:\d\d", "generated=<ts>", text)


def scenario_progress_stream(outdir: str) -> Report:
    """The gate protecting the whole progress-event design (docs/algorithm.md §10).

    Watching a run must not change it. The event stream reaches into
    :class:`~latticegen2.runlog.RunLog` and into
    :meth:`~latticegen2.parallel.WorkerPool.run`'s dispatch loop, which is the
    loop every parallel stage goes through — so "additive" is a claim about the
    hot path of the entire pipeline, and the only way to hold it is to run the
    same scenario both ways and compare the bytes.

    Cheap enough to run on every invocation: ``smoke-fast`` is ~7 s, so this
    doubles one of the four scenarios rather than adding a new class of cost.
    """
    rep = Report("progress-stream")
    print(f"=== {rep.name} ===", flush=True)
    # Same basename in two directories, so the only thing that differs about
    # the two invocations is the flag under test.
    quiet_dir = os.path.join(outdir, "progress-quiet")
    loud_dir = os.path.join(outdir, "progress-loud")
    os.makedirs(quiet_dir, exist_ok=True)
    os.makedirs(loud_dir, exist_ok=True)
    plain = os.path.join(quiet_dir, "progress.step")
    streamed = os.path.join(loud_dir, "progress.step")
    base = ["-i", BALL, "-cc", "20", "-t", "4", "--cores", "4"]

    quiet, _ = run_generator([*base, "-o", plain])
    loud, _ = run_generator([*base, "-o", streamed, "--progress-stream"])

    rep.check("both runs exit 0", quiet.returncode == 0 and loud.returncode == 0,
              f"{quiet.returncode}/{loud.returncode} {loud.stderr.strip()[-200:]}")
    if quiet.returncode != 0 or loud.returncode != 0:
        return rep

    # 1. The output is the same file.
    left = _mask_volatile(open(plain, encoding="utf-8", errors="replace").read())
    right = _mask_volatile(open(streamed, encoding="utf-8", errors="replace").read())
    rep.check("STEP identical outside the header timestamp", left == right,
              f"{os.path.getsize(plain)} vs {os.path.getsize(streamed)} bytes")

    # 2. The log is the same log. Only the front-end's channel may differ.
    def log_body(step_path: str, own_dir: str) -> str:
        """The log with everything that legitimately varies between two runs
        blanked: clock times, elapsed times, measured memory, and the run's own
        directory. What survives is every count, every geometric figure and the
        order they appear in — which is precisely what an observer may not
        change.
        """
        with open(os.path.splitext(step_path)[0] + ".log", encoding="utf-8") as fh:
            body = fh.read()
        body = re.sub(r"^\[\d\d:\d\d:\d\d\] ", "", body, flags=re.M)
        body = re.sub(r"\d{4}-\d\d-\d\d[ T]\d\d:\d\d:\d\d", "<ts>", body)
        body = re.sub(r"\d{8}-\d{6}", "<stamp>", body)
        body = re.sub(r"\d+(\.\d+)?s\b", "<dur>", body)
        body = re.sub(r"\d+(\.\d+)? [KMGT]?B\b", "<size>", body)
        # `stitch` reports its six phase timings as bare seconds under keys
        # ending `_s`, which the duration pattern above cannot reach.
        body = re.sub(r"(_s: )[\d.]+", r"\1<dur>", body)
        return body.replace(own_dir, "<outdir>")

    rep.check(".log identical outside clock, memory and path",
              log_body(plain, quiet_dir) == log_body(streamed, loud_dir))

    # 3. And the stream it produced is actually usable — otherwise the two
    #    checks above would also pass if the flag did nothing at all.
    events, junk = [], []
    for line in loud.stdout.splitlines():
        if not line.strip():
            continue
        parsed = progress.parse_line(line)
        (events if parsed is not None else junk).append(parsed or line)

    kinds = {e["ev"] for e in events}
    rep.check("stdout is pure NDJSON", not junk, f"{len(junk)} unparsed: {junk[:2]}")
    rep.check("every stage is announced and closed",
              sum(e["ev"] == "stage_begin" for e in events) == len(STAGES)
              and sum(e["ev"] == "stage_end" for e in events) == len(STAGES),
              f"{sum(e['ev'] == 'stage_begin' for e in events)} begun")
    rep.check("the run reports itself finished",
              [e for e in events if e["ev"] == "result"][-1:] != []
              and events[-1]["ev"] == "result" and events[-1]["status"] == "ok")
    rep.check("sub-stage progress arrives", "progress" in kinds,
              f"{sum(e['ev'] == 'progress' for e in events)} tick(s)")
    rep.check("the quiet run printed no events", "@@" not in quiet.stdout
              and not any(progress.parse_line(l) for l in quiet.stdout.splitlines()))
    return rep


def scenario_invalid_input(outdir: str) -> Report:
    rep = Report("invalid-input")
    print(f"=== {rep.name} ===", flush=True)
    out = os.path.join(outdir, "invalid-input.step")
    log = os.path.join(outdir, "invalid-input.log")
    # t = 4 against cc = 5 gives a = 3.54 mm, so the strut cannot fit the cell.
    proc, _ = run_generator(["-i", BALL, "-cc", "5", "-t", "4", "-o", out])
    rep.check("exits with the parameter-error code", proc.returncode == 2, f"exit {proc.returncode}")
    rep.check("prints one human-readable reason", "FAILED:" in proc.stderr,
              proc.stderr.strip().splitlines()[0][:160] if proc.stderr.strip() else "(no stderr)")
    rep.check("no .step written", not os.path.isfile(out))
    rep.check("no .log written", not os.path.isfile(log))
    return rep


SCENARIOS = {
    "smoke-fast": scenario_smoke_fast,
    "invalid-input": scenario_invalid_input,
    "progress-stream": scenario_progress_stream,
    "smoke-verified": scenario_smoke_verified,
    "dense-lattice": scenario_dense_lattice,
}


def main(argv: list[str]) -> int:
    only = None
    if "--only" in argv:
        only = set(argv[argv.index("--only") + 1].split(","))
    occ.quiet_kernel()

    reports = []
    with tempfile.TemporaryDirectory(prefix="latticegen2-e2e-") as outdir:
        for name, fn in SCENARIOS.items():
            if only and name not in only:
                continue
            reports.append(fn(outdir))
            print(flush=True)

    print("=== e2e summary ===")
    for rep in reports:
        failed = [c[0] for c in rep.checks if not c[1]]
        status = "PASSED" if rep.ok else f"FAILED ({len(failed)}): " + "; ".join(failed)
        print(f"  {rep.name}: {status}")
    all_ok = all(r.ok for r in reports)
    print()
    print("=== e2e: ALL SCENARIOS PASSED ===" if all_ok
          else "=== e2e: ONE OR MORE SCENARIOS FAILED ===")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
