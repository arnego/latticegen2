"""End-to-end verification harness (specification.md §6, docs/testing.md).

Runs every scenario in specification.md §6.1 as a subprocess and applies every
applicable check from §6.2. Each scenario is independent: one failing scenario
never prevents a later one from running, and the process exit code reflects the
run as a whole.

Usage:
    python tools/e2e.py [--only smoke-fast,dense-lattice]

`support-cylinder` is much the slowest. Use `--only` while iterating.
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
GAPS = os.path.join(TESTDIR, "support-gaps.step")
CLIPPED = {name: os.path.join(TESTDIR, f"support-clipped-{name}.step")
           for name in ("cone", "cylinder")}
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
    # Reported, never gated: the manifold check above is what decides whether
    # the shipped file survived being written, and the three cheaper quantities
    # below were each tried as a bar and two of them false-positive on sound
    # geometry (latticegen2.occ.LOOSE_AREA_FRACTION_MAX). They say why a body is
    # fragile and rank the candidates; they do not pass or fail anything.
    cos = vg.curve_on_surface(output)
    print(f"  [NOTE] export fidelity: {100 * cos['fraction']:.4f}% of the worst "
          f"solid's surface is described more loosely than its own feature size "
          f"({cos['loose_faces']} face(s)); worst pcurve deviation "
          f"{cos['worst_mm']:.3e} mm over {cos['pairs']} edge/face pair(s); "
          f"{cos['over']} pair(s) past their own recorded tolerance")

    rep.check("bounding box within input + (cc+t)",
              vg.bounding_box_within(output, input_path, cc + t))

    mesh = vg.mesh_of(output, cc, t)
    manifold_ok, bad_edges = vg.manifold_check(mesh)
    rep.check("mesh is a closed manifold", manifold_ok,
              f"{len(mesh.tris)} triangles" + ("" if manifold_ok else f", {bad_edges} bad edges"))
    no_cross, pairs = vg.self_intersection_check(mesh)
    detail = "0 crossing pairs"
    if not no_cross:
        # A crossing between two *triangles* is only a crossing between two
        # *faces* if it survives a finer ruler. Where a strut ends on a curved
        # wall, the wall patch is meshed by chords that sag below the true
        # surface while the strut's lateral triangles end on it, so at a coarse
        # deflection the two can cross with no vertex in common. Measured on
        # the 80 mm ball at `cc=20, t=4` with the lattice rotated: 3 pairs at
        # the default 0.4 mm, 0 at 0.1, 0.02 and 0.005 mm, on a solid
        # `BRepCheck_Analyzer` calls valid. Geometry that really intersects
        # itself does so at every deflection — the same rule docs/algorithm.md
        # §9 uses to turn a non-manifold *reading* into a verdict.
        coarse = len(pairs)
        no_cross, pairs = vg.self_intersection_check(vg.mesh_of(output, cc, t, 8.0))
        detail = (f"{coarse} crossing pair(s) at the default deflection, 0 when "
                  f"re-measured at an eighth of it" if no_cross
                  else f"{len(pairs)}+ crossing pairs at an eighth of the default "
                       f"deflection ({coarse}+ at the default)")
    rep.check("no self-intersections", no_cross, detail)


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
        if item.get("reason") == "not attempted":
            print(f"  [NOTE] solid {item['solid']} ({item['faces']} faces): past "
                  f"{vg.CUT_MAX_FACES} faces the exact cut is neither affordable nor "
                  f"trustworthy on this input, so it was not attempted; using the "
                  f"boolean-free containment check instead")
        else:
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
                  f"{vg.CONTAINMENT_SWEEP[-1]:g} mm; the {ev.get('probed', 0)} of "
                  f"{ev.get('unswept', 0)} measured exactly are "
                  f"{ev.get('worst_distance_mm', 0.0):.3e} mm from the input at worst"
                  + (f" (worst at {ev['worst_distance_at']})"
                     if ev.get("worst_distance_at") else "")
                  + f", against a {vg.PROTRUSION_MIN_MM:g} mm bar"
                  # An unprobed candidate is not a cleared one, and the verdict
                  # says so; the detail says how many, so a failure here reads
                  # as "not measured" rather than as "measured and outside".
                  + (f"; {ev['unprobed']} could not be measured"
                     if ev.get("unprobed") else ""))


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


CYLINDER_SUPPORT_GOLDEN = os.path.join(
    TESTDIR, "test-cylinder-cc10t1.5-rotx0y-90z0-sup1-golden-sample.step")


def support_checks(rep: Report, output: str, input_path: str, cc: float, t: float,
                   orient, thickness: float, final_print: bool = False) -> None:
    """specification.md §6.2's support checks, asked of the written file.

    A run with nothing to support at its orientation says so in its log, and
    the checks are reported as `SKIP` rather than failed: coverage of no samples
    can neither pass nor fail. A run whose lattice backs every overhang already
    says that too; coverage is checked, and the checks of the support bodies
    have no body to ask.
    """
    from latticegen2.support import LATTICE_BACKS_ALL_NOTE, NO_OVERHANG_NOTE

    log = open(os.path.splitext(output)[0] + ".log", encoding="utf-8").read()
    order = [log.find(f"stage {name}") for name in
             ("simplify", "overhang", "support", "pockets", "validate")]
    rep.check("the support stages ran between simplify and validate",
              -1 not in order and order == sorted(order))
    # A run these checks have nothing to say about, recognised by the line the
    # run itself prints for it: nothing to support at this orientation.
    if NO_OVERHANG_NOTE in log:
        for label in ("every inner overhang is backed by material",
                      "no support surface exceeds the overhang limit"):
            print(f"  [SKIP] {label} — no inner overhang at this orientation", flush=True)
        return
    cov = vg.overhang_coverage(output, input_path, cc, t, orient, thickness)
    rep.check("every inner overhang is backed by material",
              cov["samples"] > 0 and cov["uncovered"] == 0,
              f"{cov['samples']} sample(s) over {cov['area_mm2']:.0f} mm^2, "
              f"{cov['uncovered']} uncovered"
              + (f" ({cov['thin']} over material thinner than the mesh probe, "
                 f"confirmed by exact classification)" if cov["thin"] else "")
              + (f", first at {cov['where'][0]}" if cov["where"] else ""))
    ss = vg.self_supporting(output, input_path, cc, t, orient)
    rep.check("no support surface exceeds the overhang limit",
              ss["bad_area_mm2"] == 0.0,
              f"{ss['triangles']} triangles"
              + (f", {ss['bad_area_mm2']:.3f} mm^2 over the limit, worst "
                 f"{ss['worst_deg']:.2f} deg at {ss['where'][0]}"
                 if ss["bad_area_mm2"] else ""))
    # Which bodies are support is what the run says it wrote; whether they
    # overlap anything is asked of the file alone.
    if LATTICE_BACKS_ALL_NOTE in log:
        print("  [SKIP] support bodies — the lattice backs every overhang", flush=True)
        return
    found = re.search(r"^\[[0-9:]+\] supports: (\d+) separate", log, re.MULTILINE)
    rep.check("the run reports its support bodies", found is not None)
    if found is None:
        return
    n = int(found.group(1))
    dis = vg.supports_disjoint(output, n)
    rep.check("supports overlap neither the lattice nor each other",
              dis["worst_mm3"] <= vg.SUPPORT_OVERLAP_TOL,
              f"{n} support bod(y/ies), {dis['pairs']} pair(s) intersected, worst "
              f"{dis['worst_mm3']:.3e} mm^3"
              + (f" near {dis['where']}" if dis["worst_mm3"] > vg.SUPPORT_OVERLAP_TOL
                 else ""))
    sp = vg.supports_on_strut_planes(output, n)
    rep.check("support faces resting on a strut lie on its plane",
              sp["worst_mm"] <= vg.CONTACT_PLANE_TOL,
              f"{sp['resting']} face(s) resting on a strut, worst "
              f"{sp['worst_mm']:.3e} mm off its plane"
              + (f" near {sp['where']}" if sp["worst_mm"] > vg.CONTACT_PLANE_TOL
                 else ""))
    if final_print:
        final_print_check(rep, output, input_path, cc, t, orient, n, cov)


def final_print_check(rep: Report, output: str, input_path: str, cc: float, t: float,
                      orient, n_supports: int, cov: dict) -> None:
    """No overhang, keel or tip of the finished print hangs over powder: the
    input's inner surface past 60 deg, and the support bodies past the rod
    angle, each unless it rests on the wall, the lattice or another support."""
    fp = vg.final_print_overhangs(output, input_path, cc, t, orient, n_supports, cov)
    inp, sup, cand = fp["input_free"], fp["free"], fp["candidates"]
    bare = sum(inp.values()) + sum(sup.values())
    first = next((f"support {w[0][0]} at {w[0][1]}" for w in fp["where"].values() if w),
                 None) or (f"input at {fp['input_where'][0]}" if sum(inp.values()) else None)
    rep.check("the final print has no unsupported overhang, keel or tip",
              bare == 0,
              f"input bare: {inp.get('face', 0)} face / {inp.get('edge', 0)} keel / "
              f"{inp.get('vertex', 0)} tip sample(s); supports free: "
              f"{sup['face']} of {cand['face']} face / {sup['edge']} of {cand['edge']} keel / "
              f"{sup['vertex']} of {cand['vertex']} tip sample(s)"
              + (f", {sum(fp['near'].values())} within {vg.NEAR_TOL} mm of material"
                 if sum(fp["near"].values()) else "")
              + (f", first {first}" if first else ""))


def scenario_oriented_lattice(outdir: str) -> Report:
    """Print orientation on its own (specification.md §4.5), without support.

    The lattice frame is rotated to the bed and nothing else changes, so every
    §6.2 check applies as it stands. The run must also load no support code and
    announce no support stage.
    """
    rep = Report("oriented-lattice")
    print(f"=== {rep.name} ===", flush=True)
    out = os.path.join(outdir, "oriented-lattice.step")
    proc, elapsed = run_generator(
        ["-i", BALL, "-cc", "20", "-t", "4", "-o", out, "--cores", "4",
         "--orient", "35", "20", "-110"]
    )
    rep.check("process exits 0", proc.returncode == 0, proc.stderr.strip()[-300:])
    rep.check("runtime under 10 minutes", elapsed < 600, f"{elapsed:.1f}s")
    if proc.returncode == 0:
        geometry_checks(rep, out, BALL, 20.0, 4.0)
        log = open(os.path.splitext(out)[0] + ".log", encoding="utf-8").read()
        rep.check("the orientation is recorded and no support stage ran",
                  "orient: 35 20 -110 deg" in log and "stage overhang" not in log)
    return rep


def _unbacked(log: str) -> str:
    got = re.search(r"support_ceiling_points_uncovered: (\d+)", log)
    return f"{got.group(1)} reported" if got else "not reported"


def scenario_support_clipped(outdir: str) -> Report:
    """Ceiling just past the bar that clips a support node's corner, leaving
    ceiling bare beside a support: two small cuts of the rehearsal part
    (tools/make_support_gaps_part.py). Coverage is asked at 0.1 mm, since the
    0.5 mm the harness otherwise samples at can step over such a sliver."""
    rep = Report("support-clipped")
    print(f"=== {rep.name} ===", flush=True)
    orient = (180.0, 0.0, 0.0)
    for name, part in CLIPPED.items():
        out = os.path.join(outdir, f"support-clipped-{name}.step")
        proc, elapsed = run_generator(
            ["-i", part, "-cc", "7", "-t", "1.4", "-o", out, "--cores", "4",
             "--orient", "180", "0", "0", "--support", "--support-thickness", "1"]
        )
        rep.check(f"{name}: process exits 0", proc.returncode == 0, proc.stderr.strip()[-300:])
        rep.check(f"{name}: runtime under 2 minutes", elapsed < 120, f"{elapsed:.1f}s")
        if proc.returncode != 0:
            continue
        with open(os.path.splitext(out)[0] + ".log", encoding="utf-8") as fh:
            log = fh.read()
        got = re.search(r"support_ceiling_points_clipped: (\d+)", log)
        rep.check(f"{name}: the ceiling clips a support node", got and int(got.group(1)) > 0,
                  f"{got.group(1) if got else 0} point(s)")
        cov = vg.overhang_coverage(out, part, 7.0, 1.4, orient, 1.0, pitch=0.1)
        rep.check(f"{name}: backed at 0.1 mm sampling", cov["uncovered"] == 0,
                  f"{cov['uncovered']} of {cov['samples']} bare, {cov['rim']} on a rim, "
                  f"{cov['slivers']} sliver(s) within {vg.COVERAGE_SLIVER_TOL:g} mm" +
                  (f", first at {cov['where'][0]}" if cov["where"] else ""))
        geometry_checks(rep, out, part, 7.0, 1.4)
        support_checks(rep, out, part, 7.0, 1.4, orient, 1.0, final_print=True)
    return rep


def scenario_support_gaps(outdir: str) -> Report:
    """The two ways the supported rehearsal left ceiling bare, in a part that
    takes seconds (tools/make_support_gaps_part.py): a keel no face of which is
    ceiling, and a strut tip the floating-body rule dropped from a kept node.
    Both are checked to be present, so that this cannot pass by no longer
    testing them."""
    from latticegen2.classify import tessellate_surface
    from latticegen2.lattice import lattice_params
    from latticegen2.support import overhang

    rep = Report("support-gaps")
    print(f"=== {rep.name} ===", flush=True)
    out = os.path.join(outdir, "support-gaps.step")
    proc, elapsed = run_generator(
        ["-i", GAPS, "-cc", "10", "-t", "1.5", "-o", out, "--cores", "4", "--support"]
    )
    rep.check("process exits 0", proc.returncode == 0, proc.stderr.strip()[-300:])
    rep.check("runtime under 2 minutes", elapsed < 120, f"{elapsed:.1f}s")
    if proc.returncode == 0:
        with open(os.path.splitext(out)[0] + ".log", encoding="utf-8") as fh:
            log = fh.read()
        lost = re.search(r"lattice_pieces_lost: (\d+)", log)
        rep.check("a strut tip was dropped from a kept node", lost is not None,
                  f"{lost.group(1)} piece(s)" if lost else "none: the part no longer tests this")
        rep.check("the generator's own count agrees: no ceiling point bare",
                  "support_ceiling_points_uncovered: 0" in log, _unbacked(log))
        lp = lattice_params(10.0, 1.5, (0.0, 0.0, 0.0))
        body = occ.read_step(GAPS)
        tessellate_surface(body, lp)
        keel = overhang.detect(body, lp.b).counts()["edge"]
        rep.check("the part has a keel to support", keel > 0, f"{keel} edge sample(s)")
        geometry_checks(rep, out, GAPS, 10.0, 1.5)
        support_checks(rep, out, GAPS, 10.0, 1.5, (0.0, 0.0, 0.0), 1.0, final_print=True)
    return rep


def scenario_support_tilted(outdir: str) -> Report:
    """Orientation and support together, on the ball: a curved overhang, tilted."""
    rep = Report("support-tilted")
    print(f"=== {rep.name} ===", flush=True)
    out = os.path.join(outdir, "support-tilted.step")
    orient = (35.0, 20.0, -110.0)
    proc, elapsed = run_generator(
        ["-i", BALL, "-cc", "10", "-t", "1.5", "-o", out, "--cores", "4",
         "--orient", "35", "20", "-110", "--support"]
    )
    rep.check("process exits 0", proc.returncode == 0, proc.stderr.strip()[-300:])
    rep.check("runtime under 10 minutes", elapsed < 600, f"{elapsed:.1f}s")
    if proc.returncode == 0:
        geometry_checks(rep, out, BALL, 10.0, 1.5)
        support_checks(rep, out, BALL, 10.0, 1.5, orient, 1.0, final_print=True)
    return rep


def scenario_support_cylinder(outdir: str) -> Report:
    """Print support on the test cylinder, prepared for a golden sample.

    Default thickness (1 mm), turned `--orient 0 -90 0`. Two things are
    checked beyond the usual: the support checks of §6.2, and that **the same
    run without `--support` is contained in the output** — supports are separate
    bodies beside an unchanged lattice, so this is what says adding them lost
    none of it. The committed lattice golden sample cannot say that here: it is
    the unrotated lattice.

    **No golden sample exists yet, deliberately** (docs/testing.md): one is only
    worth having once a human has inspected the output it pins. `golden_check`
    reports `SKIP` until `test-cylinder-cc10t1.5-rotx0y-90z0-sup1-golden-sample.step`
    is committed, and compares against it from then on.
    """
    rep = Report("support-cylinder")
    print(f"=== {rep.name} ===", flush=True)
    out = os.path.join(outdir, "support-cylinder.step")
    orient = (0.0, -90.0, 0.0)
    plain = os.path.join(outdir, "support-cylinder-lattice-only.step")
    proc, elapsed = run_generator(
        ["-i", CYLINDER, "-cc", "10", "-t", "1.5", "-o", out, "--cores", "6",
         "--orient", "0", "-90", "0", "--support"]
    )
    rep.check("process exits 0", proc.returncode == 0, proc.stderr.strip()[-300:])
    rep.check("runtime under 20 minutes", elapsed < 1200, f"{elapsed:.1f}s")
    if proc.returncode == 0:
        geometry_checks(rep, out, CYLINDER, 10.0, 1.5)
        support_checks(rep, out, CYLINDER, 10.0, 1.5, orient, 1.0, final_print=True)
        ref, _ = run_generator(
            ["-i", CYLINDER, "-cc", "10", "-t", "1.5", "-o", plain,
             "--cores", "6", "--orient", "0", "-90", "0"]
        )
        rep.check("the reference run without --support exits 0",
                  ref.returncode == 0, ref.stderr.strip()[-300:])
        if ref.returncode == 0:
            lost = vg.contained_in(plain, out)
            rep.check("the same run without --support is contained in the output",
                      lost < 1.5 ** 3, f"{lost:.6g} mm^3 of it outside")
        golden_check(rep, out, CYLINDER_SUPPORT_GOLDEN, 10.0, 1.5)
    return rep


REHEARSAL = os.path.join(TESTDIR, "TD_HX_rehearsal_test.step")
REHEARSAL_SUPPORT_GOLDEN = os.path.join(
    TESTDIR, "TD_HX_rehearsal_test-cc7t1.4-rotx180y0z0-sup1-golden-sample.step")


def scenario_rehearsal_support(outdir: str) -> Report:
    """Print support at production scale: the rehearsal part, turned over.

    `TD_HX_rehearsal_test.step` at `cc=7, t=1.4` is the cheapest real part
    docs/testing.md names for anything touching `validate`, and with
    `--orient 180 0 0` what was its floor becomes its ceiling. This is the only
    scenario that puts the support stages through a part of that size.

    Much the slowest scenario — hours for the run and more for the checks — so
    it is **opt-in** (:data:`OPT_IN`): it runs only when named with `--only`,
    and only for final verification. Pass `--keep <dir>` with it to write into a folder that outlives
    the harness, since a run this long is not one to repeat for its log.

    **No golden sample exists yet** (docs/testing.md): `golden_check` reports
    `SKIP` until an inspected output is committed.
    """
    rep = Report("rehearsal-support")
    print(f"=== {rep.name} ===", flush=True)
    out = os.path.join(outdir, "rehearsal-support.step")
    orient = (180.0, 0.0, 0.0)
    proc, elapsed = run_generator(
        ["-i", REHEARSAL, "-cc", "7", "-t", "1.4", "-o", out, "--cores", "6",
         "--support", "--orient", "180", "0", "0", "--support-thickness", "1"]
    )
    rep.check("process exits 0", proc.returncode == 0, proc.stderr.strip()[-600:])
    rep.check("runtime under 12 hours", elapsed < 43200, f"{elapsed:.1f}s")
    if proc.returncode == 0:
        geometry_checks(rep, out, REHEARSAL, 7.0, 1.4)
        support_checks(rep, out, REHEARSAL, 7.0, 1.4, orient, 1.0)
        golden_check(rep, out, REHEARSAL_SUPPORT_GOLDEN, 7.0, 1.4)
    return rep


SCENARIOS = {
    "smoke-fast": scenario_smoke_fast,
    "invalid-input": scenario_invalid_input,
    "progress-stream": scenario_progress_stream,
    "smoke-verified": scenario_smoke_verified,
    "dense-lattice": scenario_dense_lattice,
    "oriented-lattice": scenario_oriented_lattice,
    "support-gaps": scenario_support_gaps,
    "support-clipped": scenario_support_clipped,
    "support-tilted": scenario_support_tilted,
    "support-cylinder": scenario_support_cylinder,
    "rehearsal-support": scenario_rehearsal_support,
}


#: Scenarios that run only when asked for by name with `--only`.
OPT_IN = {"rehearsal-support"}


def main(argv: list[str]) -> int:
    only = None
    if "--only" in argv:
        only = set(argv[argv.index("--only") + 1].split(","))
    keep = argv[argv.index("--keep") + 1] if "--keep" in argv else None
    occ.quiet_kernel()

    reports = []
    with tempfile.TemporaryDirectory(prefix="latticegen2-e2e-") as outdir:
        if keep:
            os.makedirs(keep, exist_ok=True)
            outdir = os.path.abspath(keep)
        for name, fn in SCENARIOS.items():
            if only and name not in only:
                continue
            if name in OPT_IN and not (only and name in only):
                continue
            try:
                reports.append(fn(outdir))
            except Exception as exc:  # a check that raises is a failed scenario
                # Scenarios are independent: one helper raising (an unreadable
                # STEP, a support count the file does not hold) is recorded as
                # that scenario's failure, and the rest still run and report.
                rep = Report(name)
                rep.check("scenario completed without an exception", False,
                          f"{type(exc).__name__}: {exc}")
                reports.append(rep)
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
