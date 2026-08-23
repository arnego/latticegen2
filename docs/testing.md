# Testing Reference Guide

Required verification and testing procedures for this project. Test code and
assets live in `test/`; the harnesses that drive whole runs live in `tools/`.

All commands assume a Python 3.11+ interpreter with the dependencies from
[../README.md](../README.md) installed. Set `LATTICEGEN2_PYTHON` if the
interpreter you want is not the default `python`.

## Unit testing

Parameter validation, lattice mathematics, the junction template, the interior
shell build, classification, the mesh gates, the connectivity rule, and the STEP
header rewrite.

- Run everything:
  ```
  python -m pytest test -q
  ```
- Run the fast, kernel-free subset (pure maths — CLI, lattice, connectivity,
  STEP header) during iterative development:
  ```
  python -m pytest test/test_cli.py test/test_lattice.py test/test_connect.py test/test_stepmeta.py -q
  ```
- Run one file:
  ```
  python -m pytest test/test_junction.py -q
  ```

| File | Covers | Needs the geometry kernel |
|---|---|---|
| `test_lattice.py` | docs/algorithm.md §2.3's verified identities, including the mutual orthogonality the architecture rests on; index ranges; naming | no |
| `test_cli.py` | Flags, ranges, the `t < a` cross-constraint, path resolution, preflight checks | no |
| `test_connect.py` | The junction graph and the floating-body rule's three outcomes | no |
| `test_stepmeta.py` | Quote-aware STEP header editing, including never overwriting a populated `FILE_SCHEMA` | no |
| `test_junction.py` | Cap integrity across the parameter range, the inradius argument behind it, and the exact `N x volume(J)` identity for instanced grids | yes |
| `test_weld.py` | Ring matching, adoption of boundary topology by the instancing index, the every-edge-twice-and-once-each-way proof, and that tiling the boundary sew (docs/specification.md §10) produces the same watertight result as sewing in one call. Also both rungs of the sew's vertex-tolerance repair (docs/algorithm.md §8), against the **real** rehearsal faces rather than synthetic stand-ins — including that neither rung replaces a topology object, which is what makes it safe on an already-proven-watertight shell, and `SpiralTest`'s own fat-tolerance face, where a *fixed* absolute cap disabled rung 2 outright. Also that a boundary layer still short of its expected free edges after the unsplit sew fails in `stitch` rather than in `assemble`, and the two properties the round-2 check rests on: that the batch validity scan behind the repair returns exactly what the per-face predicate does, and that `free_edges` does not count a degenerate edge as a hole | yes |
| `test_boundary.py` | The per-half-strut re-trim behind docs/algorithm.md §7.2, against the four real `SpiralTest.step` junctions whose intersection returns them untrimmed — including that a correctly trimmed junction is left bit-for-bit alone and that a junction wholly inside is checked and kept. Also the symmetric interface rule (docs/algorithm.md §7.1): caps are tagged not dropped, an interface needs both sides to present agreeing material, and what `resolve_interfaces` produces never trips `connect`'s invariant. Also pinhole-wire removal (§7) and the two guards on it, tested against the **real** failing junctions in `TD_HX_rehearsal_test.step` rather than synthetic stand-ins — both the `cc=5, t=1` junction the repair was built for and the `cc=12, t=2.5` one whose repair a relative-volume bar wrongly refused (G19) — see the note below | yes |
| `test_classify.py` | Distance primitives, spatial indices, ray parity, node classes, and both mesh gates — including the pole-degeneracy regression from issue #6. Also that the strided parallel sweep (docs/algorithm.md §5.4) returns *identical* classes to the serial one, across a real process boundary — identical rather than equivalent, because stride arithmetic invites off-by-ones that a tolerance would hide | yes |
| `test_main.py` | Exit codes and the "exactly one reason line" rule, before and after the log file opens | yes |
| `test_pipeline.py` | Same-domain unification's fallback ladder — a kernel that refuses to merge must yield a larger file, never a failed run | yes |
| `test_progress.py` | The NDJSON event schema (docs/algorithm.md §10): every event carries exactly the fields it declares, the reader returns `None` on kernel chatter and truncated lines rather than raising, and a dead consumer cannot abandon a run mid-pipeline | no |
| `test_runlog_events.py` | **The guarantee that watching a run does not change it.** The same sequence of writes with and without an emitter produces an identical `.log`; `stage_begin` emits without logging; `substage` rate-limits without ever dropping a stage's final count | no |
| `test_gui.py` | The front-end's testable half (specification.md §3.1): the stage weights, the reduction from an event stream to what the window shows — fed deliberately hostile input — the argv handed to the child, the cancel sentinel, and the **verbose** tick box: that unticked the pane holds only what the child wrote outside the event stream and is hidden when that is empty, that ticking it reveals the run's whole log, and that it stays enabled while a run is in flight. Several tests build widgets in one shared withdrawn Tk root, guarded by `importorskip` | no |
| `test_verify_geometry.py` | The one part of the harness whose failure mode is a plausible number rather than an exception: `material_outside`'s per-solid cut, the contradiction against a boolean-free containment check, and that a face lying *on* the input surface is a tie rather than a protrusion. The only test file that reaches into `tools/` — see the note below | yes |

## E2E verification

```
python tools/e2e.py
```

Runs every scenario in [specification.md](specification.md) §6.1 as a subprocess
and applies every applicable §6.2 check. Scenarios are independent — one failing
scenario does not prevent the others from running — and the exit code reflects the
run as a whole. `--only smoke-fast,invalid-input` restricts the set.

Each generated run writes its own `<output-stem>.log` containing the required run
data (parameters, start time, duration, run characteristics, peak memory). Save
it for analysis and attach it to the pull request.

To drive a single scenario by hand:

```
python src/main.py -i test/80mm-test-ball.step -cc 20 -t 4 -o /tmp/smoke.step -v
```

### What the harness checks

| Check | How |
|---|---|
| Exit code and one human-readable reason line | subprocess return code and stderr |
| STEP written and non-empty | filesystem |
| Parses back successfully | `STEPControl_Reader` round trip |
| **Exact B-rep validity** | `BRepCheck_Analyzer` on every solid |
| Closed manifold | every mesh edge used by exactly 2 triangles; a failure is contradicted against **exact B-rep closure** (`weld.shell_defects`), which is what the mesh test approximates |
| No self-intersections | `triangles_properly_cross` — plane-straddle pre-check, then edge piercing; a failure is contradicted against OCCT's exact `SelfInterMode` where the solid is under `SELF_INTERSECT_MAX_FACES`, and stands where it is not |
| No material outside the input body | boolean cut of output against input **per solid**, volume ≈ 0; a non-zero remainder is contradicted against a boolean-free containment check and reported as unmeasured if the two disagree (see below) |
| Bounding box within input + (cc+t) | direct comparison |
| Runtime budget | wall clock: 10 min for `smoke-fast` and `dense-lattice`, 20 min for `smoke-verified` |
| Golden-sample match | symmetric-difference volume both ways, tolerance `t³` |
| **Watching a run does not change it** | the `progress-stream` scenario runs `smoke-fast` twice, with and without `--progress-stream`, and requires an identical `.step` and an identical `.log` |

The `progress-stream` scenario is the gate the whole front-end rests on, and it
is cheap enough to run every time because it doubles a seven-second scenario
rather than adding a new class of cost. It earns its place because the event
emission is not confined to a corner: it reaches into `RunLog.line` and into
`WorkerPool.run`'s dispatch loop, which every parallel stage goes through. Two
runs and a byte comparison is the only thing that can hold "additive" as a claim
about that. Its log comparison masks the four quantities that legitimately vary
between any two runs — clock times, elapsed times, measured memory and the run's
own directory — which leaves every count and every geometric figure under test.

The self-intersection check's plane-straddle pre-check is load-bearing, not
decoration: without it, two separate solids merely *touching* along a coincident
face report hundreds of false crossings purely because their independent
triangulations are not vertex-aligned. Measured on two boxes sharing one exact
face with zero volume overlap: 344 false positives without the pre-check, 0 with
it.

The material-outside check has a similar story and a sharper edge to it. A
lattice trimmed from a body has a large share of its faces lying *exactly on*
that body's surface, which is the classic ill-conditioned input for a boolean —
and at 43,530 faces the cut duly mis-classified, reporting the **entire** solid
as outside. It did not decline: `IsDone`, `HasModified` and `HasGenerated` were
all true and it returned 43,672 faces where 43,530 went in, so a test for "the
result came back unchanged" would not have caught it. What catches it is that a
trustworthy answer here is *exactly* zero, so any remainder is contradicted
against `surface_points_outside` — an exact solid classification of the solid's
own tessellation vertices, no boolean involved. Disagreement means unmeasured,
never a pass. The same output measures exactly 0 mm³ on eight of its nine solids
and, on the ninth, 1 of 55,513 surface points outside at 1e-06 mm and none at
1e-05 mm — so nothing was actually wrong with the geometry; the check simply
could not say so.

**`spiral-stress` is the slow one, and deliberately in the default set.** Its
generation is ~7 minutes, and the containment check that replaces the exact cut
on its 45,897-face solid is a further ~31 minutes: 43,935 surface points at
~43 ms each against a swept B-spline body. `--only` exists for iterating on the
others. It is kept in the set because it is the only committed part whose
kernel-recorded tolerances reach 1e-02 mm, the only one whose dominant
component tiles, and the only one that has ever exercised docs/algorithm.md
§7.2 — four generator defects came out of it in one session
(docs/specification.md §10).

Two things follow for anyone extending this file. A boolean-based check needs a
scale at which it is known to work, and `dense-lattice`'s 15,966 faces is not
evidence about 43,530. And a check whose failure mode is a *plausible number*
rather than an exception is worse than one that raises: the harness reported a
354,733 mm³ violation with total confidence.

That second point is why `test/test_verify_geometry.py` exists at all, and it is
the only test file that imports from `tools/`. The rest of the harness is tested
by running it — a broken `manifold_check` shows up as a failed scenario — but a
`material_outside` that quietly reports the wrong number cannot be caught by
reading a passing run's output, so the contradiction logic is pinned directly.
It uses boxes rather than the real part: "contained", "sticking out" and "flush
with the boundary" are known there by construction, and the real part is covered
by this harness and by specification.md §11.

### Testing against real geometry, not a reproduction of it

`test_boundary.py`'s pinhole-wire tests load `test/TD_HX_rehearsal_test.step`
and trim two named junctions from it, rather than constructing a small synthetic
case the way the rest of the suite does. That is deliberate, and it is worth
knowing why before "simplifying" it.

The defect those tests pin was misdiagnosed for two days as a *small edge*, and
the fix built for that diagnosis was validated by a synthetic reproduction which
matched the real symptom's scale to four significant figures — a genuine,
non-degenerate ~3e-06 mm edge produced by a real boolean — swept a tolerance,
measured drift across 25 configurations, and passed. It was repairing ordinary
two-owner small edges. The real defect is a one-owner *pinhole wire* bounding no
area, which OCCT's small-edge machinery cannot see at all
(`tools/prototypes/RESULTS.md` G10).

The second junction, at `cc=12, t=2.5`, is there for the opposite reason: not a
defect the repair missed, but a valid repair a *guard* refused. Its wire is
shorter than either of the first junction's and it shifts the volume OCCT
reports seven orders further (G19), which is exactly the kind of thing no
synthetic case would have suggested was possible.

`test_weld.py`'s `self-intersecting-wire-fat-vertex.brep` is the third instance
of the same policy, and the clearest. It is one 0.0053 mm² face lifted straight
out of a failed `SpiralTest.step` run, and the only thing distinguishing it
from the two rehearsal faces beside it is a number: the shared vertex OCCT
recorded at **6.573e-02 mm**, sixteen times the fixed cap rung 2 used to carry,
where the rehearsal's are 8.7e-04 to 1.5e-03 mm. The repair logic was correct
and unchanged; only the bound was wrong, and no amount of testing against the
rehearsal's faces could have shown that, because a bound is only exercised by
input that reaches it.

So where the geometry that actually fails is committed to the repo, test against
it. A synthetic case proves the code does what you think; only the real one
proves it does what the part needs. The cost here is about 2 s.

### Golden samples

`test/80mm-test-ball-cc20t4-golden-sample.step` and
`test/test-cylinder-cc10t1.5-golden-sample.step` are reviewed baselines,
committed after manual verification of the geometry they contain. Both are
reproduced by the current implementation with a symmetric-difference volume of
0 mm³.

Regenerating one is a review decision, not a routine step: a golden sample only
has value while somebody has actually looked at it. If you do replace one,
confirm the volume-based comparison against the *existing* sample first — that
is what distinguishes a change in how the geometry is described from a change in
the geometry itself.

Both are compared with `golden_sample_volume_diff`, which cuts candidate and
golden against each other both ways and takes the larger remainder. That
comparison is a general boolean between two complete lattices, so it is slow at
lattice scale: seconds for the ball, **306 s** for the cylinder. The harness
gives it a generous budget (`GOLDEN_EXACT_BUDGET_S`) and only falls back to a
sampled equivalence check if it is exceeded — the fallback is labelled as
weaker wherever it is reported, and never turns an unmeasured result into a
pass.

A mismatch is a stop-the-line result: **never adjust a tolerance to make it
pass**, and do not assume the generator is at fault either. Investigate, and
escalate to the user for manual verification.

To compare two files directly:

```python
import sys; sys.path.insert(0, "tools")
import verify_geometry as vg
vg.golden_sample_volume_diff("candidate.step", "test/80mm-test-ball-cc20t4-golden-sample.step")
```

## Performance work

Every run writes `<output-stem>.log` with a header, one line per pipeline stage
with its wall-clock duration, template and mesh statistics, classification counts,
boundary-trim progress, the aggregate floating-body line, and the end-of-run
summary including peak memory. To iterate:

1. Run with `-v`, or just read the `.log` afterwards — it is always written in
   full regardless of `-v`.
2. Compare per-stage timings across runs to find the current bottleneck.
3. Cross-reference [algorithm.md](algorithm.md) §12 for the lever that affects
   that stage.
4. Re-run and compare the same stage.

For reference, the two committed scenarios on a 6-core / 32 GB workstation:

| Scenario | Total | Dominant stages |
|---|---|---|
| 80 mm ball, `cc=20 t=4` | ~6 s | boundary trim, export |
| `SpiralTest`, `cc=5 t=1` | ~7 min generation | boundary trim (6 min) |
| test cylinder, `cc=10 t=1.5` | ~40 s | simplify, boundary trim, stitch |
| `TD_HX_rehearsal_test`, `cc=5 t=1` | 51.7 min, 19.3 GB peak, 2.00 GB output | simplify, boundary trim, stitch, validate |

Both scenario rows and the rehearsal are post-Phase-2 (specification.md §10):
building the interior's full-strut lateral faces already merged took the
cylinder from ~56 s to ~45 s and the rehearsal from 55.3 min to 51.7 min, with
an identical output in both cases. The rehearsal figure comes from a
**controlled pair run back to back on the same machine** rather than from
comparing two sessions — its five untouched stages agree to within 1 %, where
the 2026-08-14/15 pair swung 25-36 % on machine load alone. Prefer that method
for any future performance claim here.

The cylinder's ~45 s → ~40 s since comes from two changes, measured together as
a controlled pair on `dense-lattice` (**50.30 s → 40.50 s, −19.5 %**, output
byte-identical outside the header timestamp):

| Stage | before | after | |
|---|---|---|---|
| **classify** | 10.70 s | **4.08 s** | −61.9 % — the strided parallel sweep, §5.4 |
| **validate** | 5.49 s | **2.11 s** | −61.6 % — OCCT's own parallel flag, on the master, §9 |
| boundary | 8.10 s | 7.08 s | −12.6 % — a side effect: the pool is now built before `classify`, so `boundary` no longer pays worker spawn |
| everything else | | | +1.8 % to +7.1 % |

The untouched stages drifting up ~4–7 % rather than the ~1 % an ideal pair
shows is worth naming rather than glossing: it is thermal, the two runs being
seconds apart. It does not threaten the reading, since the two changed stages
moved an order of magnitude further and in the opposite direction — but a
smaller claim than −62 % would not survive that much drift, which is the reason
to run the pair on an otherwise idle machine.

**The rehearsal row is still 51.7 min deliberately, even though the part has
been re-run since**, and the reason is a useful worked example of when *not* to
update a number. The 2026-08-17 re-run with both changes in
([profiling-reports.md](profiling-reports.md)) measured 56.8 min — slower — and
neither half of that is usable as a whole-run figure:

* another process started three seconds before `stitch` ended, so `instance`
  onward ran under competition;
* `boundary` measured 889 s at **4.28 cores**, and the 2026-08-15 profile
  recorded 14 m 51 s at **4.22 cores** on *untouched* code, labelled "no —
  variance" in specification.md §10's table. The run reproduces that slow-band
  point to within 1 %.

The per-stage results from it *are* usable where the window was clean —
`classify` 126 s → 47.0 s at 4.02 cores, `validate` 225 s → 114.0 s with its
I/O falling to exactly zero — and they are quoted in
[profiling-reports.md](profiling-reports.md) rather than here. **A stage
measurement in a clean window and a whole-run total are not the same kind of
claim**, and this run yields the first and not the second.

The third row is the scale rehearsal, first run end to end on 2026-08-14,
re-profiled on 2026-08-15 after implementing specification.md §10's paths 1–4
(its full per-stage table, both dates side by side, and the honest result for
each path are in [specification.md](specification.md) §10), and re-measured on
2026-08-17 — the figures above — on the first run of this part to pass
`validate` and write its STEP.

**`stitch` grew from 1 m 13 s to 11 m 18 s between those two profiles, and
almost all of it is the price of a correct result rather than new overhead.**
The 08-15 figure was measured on a run whose round-2 seam-only split was
silently producing a broken shell (118,760 open edges — docs/specification.md
§10). With the free-edge check that catches that in place, this part's one
tiled component fails it and is redone with a **full unsplit sew**, which the
run reports as `stitch_repaired_components: 1` — exactly the documented
fallback cost in docs/algorithm.md §8, "at the cost of the saving only for the
repaired components". So the honest comparison is against the untiled round 2,
not against 1 m 13 s. The per-face `BRepCheck_Analyzer` scan that
`occ.fix_vertex_tolerances` (§8) added is the small remainder: measured at
0.215 ms on real trimmed boundary faces, or **~1.1 min** across this part's
301,505 of them.

**`stitch` now reports its own per-phase timings**, so this is measured rather
than apportioned. On the 2026-08-17 re-run:

    round1 49.1s   split 3.0s   round2 15.1s   repair 651.2s
    retolerance 44.6s   rings 6.7s

The repair is **85 %** of the stage. Two things follow, and the second is the
one worth carrying:

* The retolerance scan is 44.6 s — the easy half, worth ~1.2 % of the run.
  **Done, 2026-08-18: 44.1 s → 22.6 s in a controlled pair** (docs/algorithm.md
  §8). What bit was not the predicate change §10 warned about but *when* the
  predicate is evaluated — repairs widen shared tolerances, so scanning every
  face before repairing anything counts the neighbours a repair fixes for free
  as unrepaired. See specification.md §11: the first explanation fitted the
  symptom and was wrong, and shipping it is what disproved it.
* **Running the unsplit sew speculatively alongside the seam-only one, so the
  discarded attempt leaves the critical path, recovers at most 15.1 s.** That
  proposal was worth building only while the discarded attempt was assumed
  expensive. It is not. The 651 s means making the seam-only split correct on
  heavily trimmed geometry — and **that is now disproved rather than merely
  hard** (G21): the sewn subset has to grow until no edge straddles it, which is
  the whole tile, at the unsplit sew's own cost. `stitch`'s repair is the price
  of a correct shell on this part, not an unfixed inefficiency.

The general lesson is one this file keeps re-learning: a stage timer that lumps
six phases together tells you the stage is slow and nothing about which proposal
would help. These timers cost nothing and retired a proposal on their first run.

What is worth knowing before doing performance work on this project at all, now
that the paths 1–4 chapter is closed:

* **The `stitch` stage is no longer a top-3 cost.** Round 2's seam-only split
  (only sew the faces round 1 left with a free edge, carry the rest through
  unchanged) plus dispatching it across the shared pool took it from 8 m 57 s
  to **1 m 13.5 s** — this single change is 95 % of the run's total
  improvement (73.1 → 47.1 min).
* **`simplify` is the largest stage and both obvious levers are now
  disproved.** It cannot be spread *below* the body: unified tiles reassemble
  by shared topology only while they stay in one process, so dispatching them to
  worker processes duplicates every seam edge (G15), while threads keep identity
  perfectly and deliver 1.04x on six of them because OCP holds the GIL for the
  whole call (G17, which retains 3.7 % of Python throughput during it). The two
  transports fix and break exactly opposite things, and there is no third. It
  also cannot
  usefully be fed *less*: restricting its face merge to the region that can
  still merge was exact (byte-identical output, 46 % less input at rehearsal
  scale) and no faster, because a correct restriction skips exactly the faces
  the kernel would have returned unchanged — the cheap ones — and keeps the ones
  that merge. Elasticity ~0.3, against 0.98 for a generic subset of the same
  size (G16), so the selection is what defeats it and no implementation can do
  better. Both are written up in docs/specification.md §11. **Do not propose a
  new way to parallelise or restrict this call without reading them first.**
* **Parallelising `simplify` across the pool is correct but was not a
  wall-clock win on this part.** It still measures at 0.99 cores in
  `profile_report.py` — this part's 14 solids are one dominant body plus 13
  small scraps, so there is nothing to spread across workers, and the added
  `.brep` round trip cost a few percent (+8 %) rather than saving anything.
  Correct behaviour on a part with more evenly sized components; not this one.
  See specification.md §10 for the full account of why this was kept anyway.
* **`validate` no longer uses the worker pool at all, and this is the one
  distinction worth carrying away from that chapter.** Per-*body* dispatch
  failed for `simplify` and `validate` alike, for the same reason. Going
  *below* the body is where they part company: `simplify` must hand back
  geometry, so its tiles must reassemble by shared topology and G15's file
  boundary destroys that; `validate` hands back a boolean, so nothing
  reassembles. G18 then found `BRepCheck_Analyzer` has a `theIsParallel` flag
  of its own — OCCT's native threads, which G17's GIL result does not bind —
  worth **1.60×** with the verdict unchanged on all four committed invalid
  faces. It runs on the master because `--cores` cannot absorb `W` processes ×
  `W` threads, which also deletes a 464 MB round trip. Controlled pair on
  `dense-lattice`: **5.49 s → 2.11 s**. The rule this leaves behind: before
  reaching for the process pool, ask whether the stage returns geometry or a
  number, and whether the kernel call already has a parallel flag.
* **`boundary` was for a long time the only stage that used more than one
  core** — parallel by construction, not by that chapter's changes. `classify`
  now joins it (docs/algorithm.md §5.4): it is pure NumPy, every node is decided
  independently, and it measured 10.33 s → 3.39 s on `dense-lattice` (3.05× on
  six cores) for a bit-identical classification. It is worth knowing *why* that
  one was easy where the rest were not — nothing OCCT crosses its process
  boundary, so neither the GIL result (G7, G17) nor the tile-identity result
  (G15) has anything to attach to. Treat that as the test to apply to any new
  candidate stage before reaching for a pool.
* **The round-trip re-import stage is gone.** It cost 22 m 29 s on
  2026-08-14 — the single most expensive stage in that run — to re-establish
  in-process what `tools/e2e.py` already checks in dev/CI on every committed
  scenario (see the check table above). Removed by deliberate decision, not
  cheapened; docs/algorithm.md §9 has the reasoning in full.

Note that the `simplify` stage (same-domain unification, docs/algorithm.md §9)
has a *negative* net cost even setting parallelisation aside: it takes ~8 s at
`dense-lattice` scale but halves the face count, which takes more than that
back out of `export`. Removing it would make the run slower as well as the
output twice as large.

### Profiling a run

Wall-clock per stage comes from the `.log`. For CPU, memory and I/O over time,
wrap the run:

```
python tools/profile_run.py --out run-profile.csv -- python src/main.py -i part.step -cc 5 -t 1 -o out.step -v
```

It samples the whole process tree — master plus every boundary and sew worker —
every 2 s and propagates the child's exit code, so it is a transparent prefix.
Then join the samples to the stage boundaries:

```
python tools/profile_report.py run-profile.csv out.log --cores 6
```

That prints a per-stage table of duration, mean and peak CPU, **core-equivalents**
(1.00 is one core fully busy), RSS, process count and I/O, and ranks the stages
where parallelism would recover the most wall time. Core-equivalents is the
number to look at: it is what shows that `boundary` uses the machine and nothing
else does.

Reports kept from past runs, with a note between consecutive ones on what
changed in between, are in [profiling-reports.md](profiling-reports.md). Append
there rather than rewriting: the series is the value, and a stage delta only
means something when the untouched stages beside it agree.

`profile_run.py` needs `psutil` for **sampling the process tree**. `psutil` is
not a runtime dependency of latticegen2 itself — the CLI's former `--ram`
budget was its only use in `src/`, and that budget was removed (specification.md
§11) — so it is not installed by any of README.md's routes and does not
ship in a release bundle. It is purely a development-only dependency of
`tools/`, which is itself `export-ignore`d and never reaches a bundle. The
script says so plainly if `psutil` is missing:
`python -m pip install psutil`.

The measurements behind the architecture's load-bearing choices (junction cap
integrity, join-mechanism throughput, per-junction intersection latency, STEP
writer throughput, boundary-stitching scaling, the seam-split repair and the
pinhole-wire repair) are in
[`../tools/prototypes/RESULTS.md`](../tools/prototypes/RESULTS.md), with the
scripts that produced them alongside. Re-run one with, e.g.:

```
python tools/prototypes/g2_instancing_join.py
```

## Release bundles

The offline bundles are built by `tools/build_release.py` and verified by
`tools/smoke_bundle.py`. Both run locally, not only in CI:

```
python tools/build_release.py --flavour both
```

```
python tools/smoke_bundle.py --platform win64
```

The smoke script is the gate that stands between a green build and a published
release, and it deliberately tests the **archive** rather than the staging tree.
Extraction is where the executable bit is lost, where a CRLF shebang bites, and
where a relocated interpreter stops finding its own libraries — none of which
appear in a build log. It extracts outside the repository, scrubs
`LATTICEGEN2_PYTHON`/`PYTHONPATH`/`PYTHONHOME` from the environment so a bundle
cannot borrow the development machine's interpreter, runs one generation through
the bundle's own launcher, and checks the result with `tools/verify_geometry.py`
(exact B-rep validity, no material outside the input, closed manifold).

Note that the builder packages a **git ref**, not your working tree, so
uncommitted changes will not appear in a bundle. It must also run under the
Python version the bundle targets when downloading wheels, since wheels are
version-specific; it says so plainly if you get that wrong.

The full release procedure is [release.md](release.md).

## Testing the Linux bundles from Windows, via WSL

The Windows dev machine can build and smoke-test the **linux-x86_64** bundles
directly, using WSL, without a separate Linux box. This needs Python 3.11 on
the WSL side (`python3` there defaults to whatever Ubuntu ships — 3.12 on
24.04 "noble" — which is too new for the pinned wheels).

**Install Python 3.11 in WSL once, ahead of time:**

```bash
sudo add-apt-repository -y ppa:deadsnakes/ppa && sudo apt update && sudo apt install -y python3.11 python3.11-venv
```

Ubuntu 24.04's own repos have no `python3.11` package at all (`apt-cache
policy python3.11` returns nothing) — deadsnakes is the standard source for
it. This needs a password prompt for `sudo`, so run it interactively yourself
rather than through a non-interactive tool call.

**Do the work from WSL's native filesystem, not the Windows-mounted repo.**
A worktree's `.git` file points at an absolute Windows path
(`E:\Git\...\.git\worktrees\...`), which WSL's `git` cannot resolve — `git
archive`, which `tools/build_release.py` depends on, fails immediately with
`fatal: not a git repository`. Clone into WSL's own filesystem instead (also
faster than building on a DrvFs-mounted path):

```bash
git clone /mnt/e/Git/latticegen2 ~/build/latticegen2 && cd ~/build/latticegen2 && git checkout <commit-or-branch>
```

**Build and smoke-test:**

```bash
python3.11 tools/build_release.py --flavour both
python3.11 -m venv ~/build/testenv
~/build/testenv/bin/python -m pip install --no-index --find-links dist/wheels -r requirements-bundle.txt
PATH="/path/to/python3.11/bin:$PATH" ~/build/testenv/bin/python tools/smoke_bundle.py --platform linux-x86_64
```

Two things that are easy to trip on:

* `tools/smoke_bundle.py` imports `tools/verify_geometry.py`, which needs
  `numpy`/`OCP` — the same "development interpreter" requirement `release.md`
  describes for the host build. Running it with a bare `python3.11` that has
  no packages installed fails at that import, not at anything the gate is
  actually checking. Use the venv above.
* The **wheels** flavour's own launcher and `install.sh` require a `python3`
  on `PATH` that is exactly 3.11 (by design — see `docs/release.md`'s
  reproducibility notes). If WSL's default `python3` is newer, the bundle
  correctly refuses to install/run, and the smoke gate reports that as a
  failure — which is the gate doing its job, not a bug. Put the 3.11
  interpreter first on `PATH` for the smoke run so it matches what CI's
  runner provides.

**Unit tests and e2e** run the same way as anywhere else, just under the venv
interpreter:

```bash
~/build/testenv/bin/python -m pytest test -q
~/build/testenv/bin/python tools/e2e.py
```

## Verification checklist

1. `python -m pytest test -q` passes. (Python has no separate lint step here;
   test failures are the signal.)
2. Edge cases for any modified boundary logic, in particular:
   - Parameter bounds (`-cc`/`-t` at 0.4/50/20) and the `t < a` cross-constraint.
   - Classification margin edges (a strut exactly at `r + d`).
   - Cap integrity at high `t/cc` ratios (docs/algorithm.md §3.3).
   - The interior shell's edge-use tally: every edge twice, opposite directions.
   - The floating-body rule's three outcomes: dropped, kept-because-connected,
     and kept-because-large.
3. `python tools/e2e.py` passes.
4. If packaging, launchers or dependencies changed: `python tools/build_release.py`
   then `python tools/smoke_bundle.py --platform <host>` passes.
5. Run `/code-review` (per CLAUDE.md) and address findings before pushing.
