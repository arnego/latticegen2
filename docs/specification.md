# latticegen2 Specification

What we currently know about the latticegen2 project. Where we don't know yet we write [TODO: needs decision] rather
than leaving it blank, so gaps are visible instead of silently assumed.
When a feature or characteristic of this project has been proposed by claude, it must clearly state so by using [TODO: proposed] rather than leaving it blank, so the user can control what enters the specification.
Do not implement something that is tagged with [TODO: needs decision] or [TODO: proposed]. 
When these change to [TODO: implement] then go ahead with implementation and clean up the TODO tag.

---

## 1. Purpose & Scope

**Goal:**
The script must generate and output a parameterized lattice geometry based on user input that fits exactly within the users boundry geometry. 
The script must use a highly optimized and parameterized generation algorithm that can taking the running hardware into account in order to ensure minimum duration to output, and good stability of the run-time system.

**Primary output:** 
A single watertight STEP representing a lattice core, filling volume defined by the solid body of the input STEP geometry, with boundry against the surfaces of the input STEP geometry, placed within the same coordinate system as the input STEP geometry. The STEP file may contain multiple bodies if the input geometry cuts rods off in such a way that some rods become floating islands disconnected from the rest.

**Secondary output:** 
Run data from the script including:
  - The runs input parameters 
  - Date and time of run start
  - Duration from start to completion
  - Run characteristics (number of tiles, number of parallel threads per stage of the generation procedure, etc)
  - Maximum memory usage

---

## 2. Deployment Target & Constraints

- **Runtime environment:** Windows 11 offline workstation and Linux command line
- **Language/runtime:** **Python 3.11+**.
- **Offline requirement:** Package must run with **zero network access**. Satisfied two ways: a published release bundle needs no network at any point (the portable flavour needs no install either), and from a checkout the dependencies are ordinary wheels installable with `pip install --no-index --find-links` — see README.md "Installation". Nothing contacts the network at run time: no package manager, license check or telemetry. The CI smoke gate proves this rather than assuming it, by installing with `--no-index` and running the extracted bundle.
- **Packaging form:** invoked as `python src/main.py <args>`, with a thin `latticegen2.bat` (Windows) / `latticegen2.sh` (Linux) wrapper provided for convenience (implemented: [latticegen2.bat](../latticegen2.bat), [latticegen2.sh](../latticegen2.sh)). No install step is needed; `pip install .` optionally provides a `latticegen2` console script.
- **Distribution form:** per-platform **offline bundles**, published as GitHub release assets by [`.github/workflows/release.yml`](../.github/workflows/release.yml) on a `v*` tag, in two flavours for each of Windows and Linux x86-64:
  - *portable* — carries a relocatable CPython with every dependency installed. Extract and run: no Python on the target, no install step, no admin rights, no network.
  - *wheels* — source plus the dependency wheels and an `install` script, for a target that already has Python 3.11.

  Each release also publishes `SHA256SUMS.txt` for verification after transfer. Every asset is extracted and run end-to-end by a CI smoke gate before publication. Procedure: [release.md](release.md). Bundle contents are `git archive`-derived, so they contain committed files only, filtered by `.gitattributes`.
- **No single-file standalone executable is produced.** PyInstaller was evaluated and rejected: `boundary.py` uses `multiprocessing` with the `spawn` start method and the codebase has no `freeze_support()` call, which on Windows makes a frozen build re-launch its own launcher; OCP/OCCT is awkward to freeze (hidden imports, DLL discovery); and freezing dissolves the LGPL-2.1 relinking argument in [licenses/LICENSES.md](../licenses/LICENSES.md), which depends on OCCT remaining a stock, replaceable shared library. The portable bundle delivers the same "extract and run" property without those costs.
- **Target machine specs / limits:** Main development system: 32 GB RAM, 6 core CPU, Nvidia RTX 3080 GPU, disk space for intermediate files.
RAM and CPU cores may optionally be provided as input parameters, as *budgets* rather than a mandatory pair. Without `--cores` the worker count is the machine's logical core count, since boundary-junction jobs are constant-size and independent; without `--ram` the budget is the memory free at startup. See §3.
- **Allowed third-party libraries:** Must be compatible with the target OS/arch. License text must be obtained and put into /licenses folder, and @/licenses/LICENSES.md must be updated with the cross reference between the library used and the corresponding license text file valid for that library.
- **License constraints:** TBD

---

## 3. Command-Line Interface

Exact invocation the human will type. This is the user-facing surface.

For each parameter, specify: **name, type, units, valid range, default, required?**

| Flag | Type | Required | Units | Range | Default | Description |
|------|------|----------|-------|-------|---------|-------------|
| -i --input | path | required | NA | NA | NA | Path to STEP file defining the lattice bounds |
| -o --output | path | optional | NA | NA | `<input_stem>-cc<cc>t<t>.step` | Path and name of the output .step file. Must name a **file**, not a directory: `-o .\` and friends are rejected (exit 2) rather than turned into `.\.step`. `.step` is appended if absent. |
| -cc | float | required | mm | 0.4 - 50 | NA | Distance between the bottom nodes of two adjacent cells |
| -t | float | required  | mm | 0.4 - 20 | NA | Side length of the diamond rod profile. Must be smaller than the cell edge `a = cc/√2`; that is the only cross-constraint. |
| -v --verbose | flag | optional | NA | NA | disabled | Enable verbose console diagnostics while always writing a full `.log` file. |
| --cores | int | optional | count | 1 - 128 | logical cores on the machine | Maximum CPU cores this run may use. One worker process per core, honoured exactly — the master needs none reserved for it, being blocked waiting on results for effectively the whole boundary stage. Since workers always run at below-normal priority, this exists to further protect the response time of the system for other tasks. |
| --ram | float | optional | GB | 1 - total physical RAM detected | free RAM at startup | Maximum memory this run may use. May be set above or below what is currently free, but never above the machine's total physical RAM. Advisory: recorded in the run log next to the measured peak. |

Both are optional **budgets** rather than a mandatory pair, and both resolve to a
concrete figure either way: an explicit value is honoured exactly, and an omitted
one is taken from the machine — logical core count for `--cores`, free memory at
startup for `--ram`. Detection lives in
[`src/latticegen2/sysinfo.py`](../src/latticegen2/sysinfo.py).

**Process priority is not a parameter.** Every run — master and every worker —
executes at below-normal priority so the machine stays usable for other work.
This was the opt-in `-bg` flag through v2.x; it is unconditional now, since a
choice whose only alternative is "make the desktop unusable" is not worth
offering. Implemented by `latticegen2.parallel.set_background_priority`, called
once in `__main__` and once per worker from the pool initializer.

**Exit:** 

Upon success the script shall produce an end of run summary report in the .log file and to console independent of the -verbose flag. This shall include: 
 - The runs input parameters 
 - Date and time of run start
 - Duration from start to completion
 - Run characteristics (node classification counts, boundary pieces, worker count, connected components, face/vertex/edge counts of the assembled shell, etc)
 - Maximum memory usage
 - Path to output .step file 

Upon failure, the script shall output a human readable reason for the failure, e.g.: parameter bounds exceeded, issues with input geometry, issues with resarouces from the run-time system, write or read access issues, etc.

Upon cancellation by the user (Ctrl+C), the script shall shut down gracefully rather than terminate abruptly: worker processes are stopped in an orderly way (force-stopped only if they do not respond within a short grace period), a single human readable `CANCELLED` line is written to console and `.log` file, the temporary folder is left in place for analysis, and the exit code is 130. See docs/algorithm.md §10.
  
**Logging:**

A log file should be produced every run with the same name as the output file which is generated from the input file or provided by the -o flag. The log file should end with `.log` and should not include .step (that is only the last name for the geometry file.  

---

## 4. Geometry Domain Specification

### 4.1 Lattice unit cell type
The base geometry for the lattice is a strut-based uniform grid forming cube-like cells standing on its tip. The struts form the boundries of the cells along each edge. The struts have the profile of a square standing on one corner, like a diamond. In other words; for each strut the square profile is oriented so one diagonal lies in the vertical plane containing the strut axis and the Z-axis, and the other diagonal is horizontal. See docs/algorithm.md §3.1 for the exact frame construction. Verified: [`profile_vertices`](../src/latticegen2/lattice.py) builds each profile from `u_k = normalize(cross((0,0,1), e_k))` (horizontal) and `v_k = cross(e_k, u_k)` (vertical-plane) exactly per docs/algorithm.md §3.1, and `test/test_lattice.py` asserts both the frame orientation and that the profile is a square of side `t`.

The dimentions of the square profile is defined by input parameter `t`. Upon inspection of the end result, the struts are reclined from the normal axis (Z-axis) in degrees by the following calculation in numpy: np.degrees(np.arcsin(np.sqrt(2/3)))
Make sure to use the exact expression rather than a decimal literal. It should be close to 55 degrees (but not exactly).
The distance between base points of each cell on the XY plane is defined by input parameter `cc`.
Upon inspection the rods protruding up from the xy plane from an intersecting node are separated by an angle of 120 degrees around the Z-axis.

### 4.2 Parametrization
- The sides of the diamond shaped square rod is defined as `t` in millimaters
- `cc` is the XY-plane distance between the bottom nodes of two adjacent cells (consistent with §4.1). The cube edge length is therefore `a = cc / √2`.
- The bounds of the generated lattice and its placement in the xyz coordinate system is defined by the input step file. 



### 4.3 Boundary / shell requirements
- The lattice shall not have an outer solid shell generated around the build volume. It will be merged with the outer shell upon import into the enveloping part. However the truts must be closed against the geometry of provided input STEP file.
- No fillets/chamfers at strut junctions or bounding geometry

### 4.4 Performance & Optimization

See [algorithm.md](algorithm.md) for the full normative algorithm specification,
including the exact lattice math, the pipeline and classification diagrams, and
the detailed optimization strategy this section summarizes.

Since this involves computational geometry:
- Profile geometry generation routines to identify bottlenecks
- Consider vectorization or parallelization
- Cache expensive calculations — the dominant instance of this is the junction
  template, computed once per run and instanced at every node (algorithm.md §3.2)
- If caching to disk is used, put the files in a temporary folder `temp/<date><time>` where the output file is generated to. Clean up after a sucessful run. Leave for error analysis if the run fails.

---

## 5. STEP Output Requirements
- **Clean up** Never produce a floating body with volume < t³ mm³ (i.e. a cube of side `t`).
  This rule targets **floating (disconnected) bodies only** — a solid is only ever
  discarded once it is verified to have no geometric connection to the rest of the
  output, never merely because its own volume happens to fall below t³. A
  sub-threshold solid that is still connected to other geometry is kept.
  Connectivity is **proven by construction** rather than resolved
  experimentally: two junctions are joined exactly when they share a surviving
  mid-strut interface, so the rule is a connected-components query over a graph
  (docs/algorithm.md §8). There is consequently no "cannot determine
  connectivity" case. Note that the distinction this rule draws is not academic:
  a boolean intersection can leave sub-threshold junction wedges that are
  genuinely *connected* material, and reading the rule as an unconditional
  "volume < t³ → delete" would punch holes in the output.
  
- **STEP schema/AP:** AP214

- **Geometry representation in the file:** exact B-rep solid
  
- **Units:** mm
 
- **Metadata to embed:** Part name as concetenated <input_file_name>+cc<cc>+t<t> and generation parameters as STEP header

- **Downstream tool(s) that will open this file:** Soldiworks and Catia

---

## 6. Autonomous End-to-End Verification


### 6.1 Test scenarios
List concrete parameter sets that must be run automatically (at minimum: one small
case, one large/dense case, one edge case at parameter boundaries, one expected-failure
case for invalid input).

All four are implemented in [`tools/e2e.py`](../tools/e2e.py) and all four pass.

| Scenario | Parameters | Expected result |
|----------|-----------|------------------|
| smoke-fast | -i test/80mm-test-ball.step -cc 20 -t 4 --cores 4 | generation < 10 minutes. **Measured: 6.4 s.** |
| smoke-verified | -i test/80mm-test-ball.step -cc 20 -t 4 --cores 4 | valid STEP, generation < 20 minutes, matching golden sample test/80mm-test-ball-cc20t4-golden-sample.step. **Measured: 6.3 s, symmetric-difference volume 0.0000 mm³.** |
| dense-lattice | -i test/test-cylinder.STEP -cc 10 -t 1.5 --cores 6 --ram 20 | valid STEP, no self-intersections, matching golden sample test/test-cylinder-cc10t1.5-golden-sample.step, generation < 10 minutes. **Measured: 47.5 s, symmetric-difference volume 0 mm³.** |
| invalid-input | -i test/80mm-test-ball.step -cc 5 -t 4 (strut size `t` >= cell edge `a=cc/√2`) | exits 2, no `.step` or `.log` file written, one human-readable reason line. **Passes.** |

### 6.2 Automated pass/fail checks
For every scenario the harness must verify, without human intervention:
- Process exits with expected console output.
- STEP file is written and non-empty.
- STEP file parses back successfully (round-trip read).
- Geometry is a valid closed manifold solid (no open edges / non-manifold edges).
- **Geometry passes OCCT's exact B-rep validity check** (`BRepCheck_Analyzer`) —
  an exact test on the B-rep itself, not an inference from a tessellation.
- No self-intersections.
- **No generated material lies outside the input body** (boolean cut of output
  against input leaves ~zero volume) — a direct check of §1's "fits exactly
  within the user's boundary geometry", independent of any golden sample.
- Bounding box of output matches requested `--input` within tolerance.
- Runtime stays under an agreed performance budget: `smoke-fast` and
  `dense-lattice` < 10 minutes, `smoke-verified` < 20 minutes.
- If a golden sample is defined, check similarity of geometries by subtracting
  candidate and golden both ways; the larger remainder must be near zero.


### 6.3 How verification runs offline
- Verification runs only in the dev/CI environment.
- Test runner: `pytest` for unit tests in `test/` (alongside the STEP assets the
  scenarios reference); whole-run harnesses and geometry checks are in `tools/`
  (`e2e.py`, `verify_geometry.py`). Both run offline — the only extra dependency
  over the tool itself is `pytest`.
- Results are reported as console summary for analysis and addition to the pull-request.

---

## 7. Error Handling & Edge Cases

- Invalid/out-of-range parameters should be rejected before any computation starts.
- Read and write failures should be reported and result in a hard fail. Existing files can be overwritten.

---

## 8. Non-Functional Requirements

TBD

---

## 9. Open Questions / Decisions Needed

*Anything you're unsure about — list it here explicitly so it doesn't get silently
assumed by default. Delete each line once resolved.*

---

## 10. Roadmap features or bugs to fix in later sessions

*Concrete, actionable work items discovered but deliberately not fixed in the session
that found them. Each item should carry enough context (what's broken, where, why, and
how to verify the fix) that a later session can act on it without re-deriving the
diagnosis. Remove an item once it's fixed and verified.*

### `--ram` is accepted and validated but never enforced

**Found 2026-08-17**, while reading the rehearsal's memory profile
([profiling-reports.md](profiling-reports.md)). Not a regression — it has always
been this way — but the rehearsal came close enough to the number the user
supplied that the gap is worth closing or documenting honestly.

**What's wrong.** `--ram` is parsed, range-checked against the machine's
physical RAM, resolved to `Args.ram_budget_gb` (`cli.py`), and then **only
printed**. Nothing in `pipeline.py`, `parallel.py`, `weld.py`, `boundary.py` or
`runlog.py` ever reads it. §3's table is accurate — "Advisory: recorded in the
run log next to the measured peak" — but `cli.py`'s own module docstring opens
by calling both budgets "ceilings on what a run may use, not hints it may
exceed". That is true of `--cores`, which resolves to the worker count and is
genuinely enforced, and false of `--ram`. **The two statements contradict each
other and one of them has to change.**

**Why it matters now.** The 2026-08-17 rehearsal was run with `--ram 20` and
peaked at **19,291 MB** of whole-tree RSS — inside 3 % of the stated budget,
with nothing watching. Had it crossed, the run would have continued exactly as
it did: no warning, no throttle, no failure, and a summary reporting a peak
above the budget without remarking on it.

**What happens on real exhaustion, since it is not what the exit codes suggest.**
Windows grows the pagefile first, so the run degrades rather than fails — worst
of all in `export`, which is both the memory peak and CPU-bound serialisation.
If commit charge is then exhausted, a Python-level `MemoryError` is caught by
`__main__`'s `except Exception`, reported as one `FAILED: unexpected ...` line
and re-raised, so it surfaces as **exit 1 with a traceback**; an allocation
failure inside OCCT may instead abort natively with no reason line at all.
Exit 5 ("Resource limits") is unreachable (docs/algorithm.md §10), so a genuine
out-of-memory reads as a bug rather than as the resource limit it is.

**A second, independent defect: the reported peak under-counts.**
`RunLog.max_rss` is the master's own `PeakWorkingSetSize`, plus worker peaks
only where a stage explicitly folds them in with `note_worker_rss` — `boundary`,
`stitch`, `simplify` and `validate` do; no other stage does. On the rehearsal
the tool reported **18.12 GB** where external whole-tree sampling measured
**19,291 MB**. Any future budget check written against `max_rss` inherits that
~1.2 GB under-count, so this needs fixing first either way.

**Two ways forward. [TODO: needs decision]**

*(a) Keep it advisory, and say so consistently.* Correct `cli.py`'s docstring so
only `--cores` is described as a ceiling, and have the end-of-run summary state
plainly when the measured peak exceeded the budget. Cheap, no behaviour change,
and it makes the number honest rather than decorative.

*(b) Actually enforce it.* Harder than it looks, and worth being clear about why
before anyone starts: **the run's peak is not in the workers.** It is the master
holding the whole 2 GB result while `export` serialises it, so capping worker
count or tile concurrency — the obvious levers, and the ones Phase 3's risk R5
proposes — cannot lower the number that actually sets the peak. Real enforcement
would mean changing how the result is held and written (streaming the STEP
write, or exporting per solid), which is a substantially larger change than the
flag suggests. Anything less would be enforcement in name only.

**How to verify a fix.** For (a): run with `--ram` set deliberately below a known
peak and confirm the summary says so, with the run still succeeding. For either:
fix the under-count first — fold every stage's worker peak into `max_rss`, or
have the master sample its own process tree — and check the reported peak
against `tools/profile_run.py` on the same run; today they differ by ~1.2 GB and
should agree. `python -m pytest test -q` and `python tools/e2e.py` must stay
green; `test_cli.py` already pins the budget's parsing and defaulting.


### `stitch` pays the full round-2 cost on heavily trimmed parts

**Found 2026-08-17**, on the first rehearsal of `TD_HX_rehearsal_test.step` at
`cc=5, t=1` to complete. `stitch` took **11 m 18 s**, against 1 m 13 s in the
2026-08-15 profile. Not a regression, and worth recording precisely so it is not
re-diagnosed as one.

**Why.** Boundary-sew round 2 sews only the free-edge-bearing subset of each
tile's faces (algorithm.md §8, gate G8). That identity holds at every prototype
scale but **not** on this part's real, heavily trimmed pieces, where straddling
edges let sewing split one shared edge into two. The per-component free-edge
check against `want_rings` catches it and redoes round 2 on the unsplit tile
results — this run reports `stitch_repaired_components: 1`, so the dominant
component pays close to the untiled round-2 cost. That is the documented
fallback behaving exactly as designed; the saving is simply unavailable here.

**What would recover it.** Make the seam-only split correct in the presence of
straddling edges — i.e. carry a seam face's straddling neighbours into the sewn
subset so `BRepBuilderAPI_Sewing` cannot rebuild a shared edge onto a new
`TopoDS_Edge` while the carried face keeps the original (that mechanism is
written up in §8 and in `tools/prototypes/RESULTS.md` G8). 144–720 such edges
were measured in every prototype block tried once the check for them existed, so
they are easy to enumerate; what is unproven is whether including them keeps
round 2's cost below a full sew on a part this size.

**A separate, much smaller item in the same stage.** `occ.fix_vertex_tolerances`
(§8) runs one `BRepCheck_Analyzer` per boundary face, serially on the master:
0.215 ms measured on real trimmed faces, ~1.1 min across this part's 301,505.
It is embarrassingly parallel and could dispatch across the shared `WorkerPool`
like `simplify` and `validate` do — though note specification.md §10's own
finding that doing so bought nothing on this part's very unequal components.

**How to verify either fix.** Re-run the rehearsal; `stitch_repaired_components`
must be 0 (for the first) and `stitch` must drop, while the run still writes 14
valid solids and `python -m pytest test -q` plus `python tools/e2e.py` stay
green. `test_weld.py` already pins that tiling produces the same watertight
result as sewing in one call, which is the property any change here must keep.

### Scale rehearsal, chapter closed: paths 1–4 implemented and re-measured

**First run 2026-08-14**, **re-profiled 2026-08-15** after implementing paths
1–4 below, both on `TD_HX_rehearsal_test.step` at `cc=5, t=1`,
`--cores 6 --ram 20 -bg` on the 6-core / 32 GB development workstation. Both
runs used a temporary, uncommitted bypass of the `assemble`-stage watertightness
gate ("Micron-scale debris edges" below, still open and still out of scope for
this chapter) so every later stage could be measured; neither run's output is a
shippable file, but the stage costs are representative because the same work is
done either way, and both runs agree on everything the defect doesn't touch —
330,354 mm³, 14 solids, 705,000 interior faces, 301,505 boundary faces,
1,006,505 total faces pre-unification — which is itself evidence the bypass
didn't quietly change what work got measured.

Boundary-sew tiling (§8) and the disagreeing-cap fuse (§7.1) were verified at
this scale and retired into [algorithm.md](algorithm.md) on 2026-08-14. This
entry closes the chapter's remaining three items — parallelise `simplify`
(path 1), a cheaper round-trip gate (path 2), parallelise `stitch` round 2
(path 3), parallelise `validate` (path 4) — all four implemented, and all four
mechanisms documented normatively in [algorithm.md](algorithm.md) (§8 for
stitch round 2 and the shared pool, §9 for unification/validation and the
removed round-trip check, §12 for the updated cost model). What follows is the
measured result, including one negative one.

**A production-scale bug was caught before merge and is worth recording
alongside the win.** The first implementation of round 2's seam-only split
(§8) computed free edges as a plain Python list and tested every face against
it with `.IsSame()` — fine at prototype scale (G8, hundreds to low thousands
of faces) but `O(faces²)` in effect, and at this rehearsal's scale it took
`stitch` from 8 m 57 s to **51 m 07 s**, a 5.7× regression. Rewritten to use
`TopTools_IndexedMapOfShape` (OCCT's own shape-identity map) for near-`O(1)`
membership tests, `stitch` dropped to 1 m 13.5 s — see below, and the full
account in `tools/prototypes/RESULTS.md` G8. The lesson: a correctness gate
proven at a few hundred faces is not a performance gate, and a design meant
for hundred-thousand-face parts needs at least one measurement taken there.

#### Per-stage cost and resource profile, 2026-08-15

Wall time from the run's `.log`; CPU, memory and I/O from
[`tools/profile_run.py`](../tools/profile_run.py), joined to the stage boundaries by
[`tools/profile_report.py`](../tools/profile_report.py). "Cores" is mean CPU over the
stage, where 1.00 is one core fully busy and 6.00 is the machine. The `2026-08-14`
column is the pre-optimization baseline; stages this chapter did not touch are
included for the total but flagged, since some genuinely shifted between the two
sessions from ordinary machine-load variance rather than any code change.

| Stage | 2026-08-14 | 2026-08-15 | Cores (08-15) | RSS peak (08-15) | Touched? |
|---|---|---|---|---|---|
| template | 0.04 s | 0.06 s | — | — | no |
| import | 0.40 s | 0.46 s | — | — | no |
| tessellate | 2.4 s | 3.1 s | 0.41 | 262 MB | no |
| classify | 1 m 54 s | 2 m 27 s | 0.89 | 309 MB | no — variance |
| boundary | 10 m 57 s | 14 m 51 s | 4.22 | 2,000 MB | no — variance |
| connect | 11.0 s | 11.0 s | 0.99 | 2,236 MB | no |
| stitch | 8 m 57 s | **1 m 13.5 s** | 1.89 | 3,190 MB | **yes — path 3** |
| instance | 1 m 07 s | 1 m 07.2 s | 0.99 | 5,011 MB | no |
| assemble | 29.1 s | 30.0 s | 1.00 | 4,597 MB | no |
| simplify | 17 m 17 s | 18 m 39 s | 0.99 | 9,886 MB | **yes — path 1** |
| validate | 2 m 59 s | 3 m 29.6 s | 0.99 | 16,416 MB | **yes — path 4** |
| export | 6 m 42 s | 4 m 21.5 s | 0.96 | 21,226 MB | no — variance |
| verify | 22 m 29 s | *(removed)* | — | — | **yes — path 2** |
| **total** | **73.1 min** | **47.1 min** | | **19.61 GB** | **35.6 % shorter** |

`boundary` and `classify` are untouched by this chapter — same code, same
input — yet both measure 25–36 % slower on 08-15 than 08-14. That is
environmental (the two sessions ran on different days under different machine
load), not a regression, and it is the reason the table reports both dates
side by side rather than a single "before/after" pair: a stage-by-stage
comparison is only meaningful once it is clear which deltas are code and which
are noise. `export`'s 35 % *improvement* despite being equally untouched cuts
the other way and reinforces the same point — take the touched-stage deltas
below at face value, not the untouched ones.

#### Path 3 — `stitch`: the headline win, confirmed at scale

**8 m 57 s → 1 m 13.5 s, a 7.3× improvement over the already-tiled baseline**
(16.7× against the never-tiled 20 m 27 s from 2026-08-14's own control), on the
identical 21,955-piece / 35-tile component. This is what closes 95 % of the
25.9-minute total reduction (73.1 → 47.1 min) — every other stage's net change
roughly cancels (simplify and validate got slightly slower; export got faster;
classify and boundary's apparent slowdowns are the variance noted above).

Two independent levers, both in [algorithm.md](algorithm.md) §8: round 2 now
dispatches per component across the run's shared `WorkerPool` (generality —
this part's 14 components are 1 dominant tiled one plus 13 small untiled ones,
so there is only ever one real job to parallelise, and the gain from this
lever alone is close to zero here); and round 2 now sews only the
free-edge-bearing subset of each tile's result, carrying the rest through
unchanged (G8, `tools/prototypes/RESULTS.md`) — this is the lever that
actually moved the number, since it cuts what round 2 pays its flat per-face
cost on rather than just parallelising an unchanged cost. A hierarchical tree
reduction was considered and rejected on paper before either lever was built
(algorithm.md §8): round 2's cost tracks total face count almost flatly in
shape count, so a tree pays that flat cost once per level for nothing.

#### Paths 1 and 4 — `simplify` and `validate`: correct, but no wall-clock win on this part

**This is the honest negative result of the chapter.** Both stages dispatch
across the shared pool exactly as designed — G7 (`tools/prototypes/RESULTS.md`)
measured OCP holding the GIL around both `unify_same_domain` and `is_valid`, so
this is the same process-pool-plus-`.brep` mechanism as everywhere else, not
threads — and `profile_report.py` confirms the dispatch is real: both still
measure at 0.99 cores, i.e. *no* parallel speedup materialised. `simplify` went
from 17 m 17 s to 18 m 39 s (**8 % slower**); `validate` from 2 m 59 s to
3 m 29.6 s (**17 % slower**).

The cause is exactly what was flagged as a risk before this was built: **the
largest single solid is the floor, not the sum**, and on this part that floor
*is* essentially the whole workload. Of the 14 solids, 13 are small
floating-body-scale remnants that unify and validate in a fraction of a
second; one dominates almost completely. Dispatching 14 jobs across 6 workers
does not help when 13 of them are nearly free and the 14th has to run alone
regardless — and the process-pool path adds a `.brep` write/read round trip
per solid that a single in-process loop never paid, which is where the slowdown
comes from. This was measured, not assumed away, exactly as the risk note
said it should be.

This is not a reason to revert the change. Per [algorithm.md](algorithm.md)
§11, an optimization's failure mode must be "do more work", never "produce a
wrong result" — and a few percent of added `.brep` I/O on a part shaped like
this one is exactly that failure mode, not a correctness problem. A part whose
components are more evenly sized (several separate floating islands of
comparable scale, rather than one dominant body plus scraps) would see the
intended benefit; `TD_HX_rehearsal_test` at `cc=5, t=1` simply is not that
shape. The mechanism is sound and stays; the benefit is part-shape-dependent,
and that dependency is now documented rather than assumed.

#### Path 2 — the round-trip gate: removed, not cheapened

Per the user's decision, `round_trip_check` was deleted outright rather than
made cheaper (the original path 2 proposal). It cost **22 m 29 s** — the single
most expensive stage in the 2026-08-14 run — to re-parse the 2.00 GB output to
full B-rep purely to count solids, for a guarantee `tools/e2e.py` already
establishes independently, on every committed scenario, in dev/CI
(`vg.brepcheck`, a real `STEPControl_Reader` round trip). See
[algorithm.md](algorithm.md) §9 for the removal rationale in full.

#### Output size

2.11 GB, 611,651 faces, 14 solids (2026-08-14: 2.00 GB, 584,028 faces). The
face-count difference between the two runs is same-domain unification's own
representation choice, not a geometry difference — both runs report
`unmerged_solids: 0` and the same volume-drift figure (1.60e-07), so every
solid fully unified both times; which faces end up coincident enough to merge
can differ slightly depending on the exact face objects a run's boundary sew
happened to produce (round 2's seam-only split changed *which* face objects
those are, though not the shell they describe — see path 3 above), and
unification is a size optimization, not a correctness one, so this difference
is expected and harmless (algorithm.md §9, §11).

---

## 11. Closed — kept for the reasoning, not as work

### Invalid boundary faces from grazing trims — FIXED 2026-08-17 (34 → 0)

Tracked from 2026-08-16, when the seam-split (#13) and pinhole-wire (#15) fixes
first got the `cc=5, t=1` rehearsal of `TD_HX_rehearsal_test.step` to assemble
**14 watertight solids** — and it then failed `validate`, with 1 of those 14
solids carrying **34 individually invalid faces** (exit 4, no output written).

Closed in two passes, because it was **two unrelated faults wearing one
symptom**. Both are repaired by `occ.fix_vertex_tolerances`, as two rungs, and
both turn out to be a *recorded tolerance* being wrong rather than any geometry
being wrong (docs/algorithm.md §8; `tools/prototypes/RESULTS.md` G11 and G12).

* **Rung 1 (#17) — 34 → 4.** Not 34 bad faces but **17 bad edges**: each
  invalid face had exactly one invalid edge, and faces came in pairs because
  the two faces of a pair *shared* it. On each, a vertex sat off the edge's own
  3D curve (2.474044e-05 mm on an ellipse, 3.316370e-04 mm on a B-spline) with
  that vertex's tolerance inflated to *exactly* that distance, so the check sat
  on the knife edge. Repaired with `ShapeFix_Edge.FixVertexTolerance`.
* **Rung 2 (this change) — 4 → 0.** The remaining 4 had **no** standalone-invalid
  edge or vertex and passed every named face check, which G11 read as a
  *contextual pcurve-versus-3D-curve* fault. **That reading was wrong.** The
  fault is `BRepCheck_SelfIntersectingWire`, for two edges *adjacent in the
  wire* — and it is not a real self-intersection: their pcurves cross at exactly
  one point, at the shared vertex, *inside* its tolerance. The shared vertex's
  recorded tolerance was simply left too tight to swallow the crossing.
  Repaired by widening that vertex, bounded at 4× the tolerance the kernel
  itself recorded and at 4e-3 mm absolutely.

**What the wrong diagnosis cost, and how it was caught.** The pcurve deviation
G11 named is real — on each of the four faces one edge deviates by 98–100 % of
its own tolerance, which reads as a root cause. Two measurements disproved it:
widening that edge's tolerance fixes nothing at any factor up to 5×, and
`BRepLib.SameParameter` improves the deviation on the one face it cannot fix.
The mechanism was then found with `BRepCheck_Analyzer.IsValid(subshape)` — the
*in-context* overload — run as a **controlled** probe: it fires on all four and
clears on the three `SameParameter` does fix, so a negative on the fourth could
not have been the probe being blind. That control is the lesson from G10, where
a blind scanner reporting "no faults" cost real time.

This makes **four** defects in this family, each first diagnosed wrongly in a
way that matched the symptom quantitatively and convincingly (G9, G10, G11,
G12). The pattern worth keeping: OCCT's named repair for the named symptom did
not touch the actual defect in any of them — `ShapeFix_Wireframe` for the
"small edges" that were pinhole wires, `ShapeFix_Wire.FixSelfIntersection` for
this self-intersection (a complete no-op). Establish the mechanism on the real
failing geometry before building the repair.

**Result, measured on the full rehearsal 2026-08-17** (`-cc 5 -t 1 --cores 6
--ram 20 -v`, 58m 18s, 18.57 GB peak):

* `vertex tolerances corrected on 19 sewn boundary face(s); no geometry moved`,
  with **no residual** — against 15 corrected and 4 remaining before this change.
* `validity: all 14 solid(s) pass BRepCheck_Analyzer` — the gate passing for the
  first time on this part.
* The run **wrote its STEP**: 2.01 GB, 14 solids, 584,028 faces, lattice volume
  330,354.002 mm³, AP214, part name `TD_HX_rehearsal_test+cc5+t1`. Every figure
  the defect does not touch is unchanged from the 2026-08-15 profile (21,955
  boundary pieces, 122,180 interfaces, 2 pinhole wires, 1,006,505 faces before
  unification, same volume), which is evidence the repair changed only what it
  was aimed at.

`python -m pytest test -q` and `python tools/e2e.py` stay green (0 mm³ against
both golden samples); the committed scenarios contain no grazing trims this
severe, so neither rung fires on them, exactly as intended.

### Micron-scale debris edges from near-tangential trims — FIXED 2026-08-16

Tracked from 2026-08-14 as "2 edges used by exactly one face" at
`[1874.836, 60.370, 970.121]` and `[1874.836, 59.912, 969.775]`, with three
candidate repairs of which two were disproved and the third — `ShapeFix`
small-edge removal in the worker — was recommended as the remaining one.

**The diagnosis was wrong, and convincingly so.** They are not small *edges*.
Each is an **inner wire of a single edge whose endpoints do not meet**, bounding
no area — a pinhole. That is why each is used by exactly one face inside a solid
`BRepCheck_Analyzer` calls valid, and why OCCT's small-edge machinery cannot see
them: `ShapeFix_Wireframe.CheckSmallEdges` reports **0** candidates at every
precision from 1e-5 to 1e-2, and `ShapeFix_Face.FixSmallAreaWire` removes
nothing. **Do not retry `ShapeFix` here.**

Fixed by `occ.remove_pinhole_wires`, called from `boundary.trim_junction`
before cap tagging: it drops a non-outer wire when every edge in it is under
`PINHOLE_WIRE_TOL` **and** every edge is already used exactly once, so it can
only ever delete edges that are already unpaired. Surface area is preserved
bit-identically, volume to 2.7e-15, and `shell_defects` on the affected piece
goes from `(2, 0)` to `(0, 0)`. On the full rehearsal it removes 2 wires from 1
junction of 19,552 and `assemble` reports 14 watertight solids. Full account,
including why a synthetic reproduction of this defect passed a complete
measurement gate while repairing something else entirely, in
`tools/prototypes/RESULTS.md` G10 and docs/algorithm.md §7.

### `simplify` beyond body-for-body: three ranked optimizations

**Proposed by Claude 2026-08-15, approved for implementation 2026-08-16.**
**Phases 1 and 2 are done and measured at rehearsal scale; Phase 3 is not
built.** Phase 1's headline projection was wrong and Phase 2's was too
optimistic — both are recorded below rather than quietly dropped, because what
they cost to learn is the useful part. The measurement the proposals rest on is
gate G13
([`tools/prototypes/RESULTS.md`](../tools/prototypes/RESULTS.md),
`tools/prototypes/g13_unify_scaling.py`).

**Why this exists.** Path 1 above — parallelising `simplify` across the shared
pool, body for body — is implemented, correct, and measured at **no wall-clock
win** (17 m 17 s -> 18 m 39 s, still 0.99 cores), because this part's 14 solids
are one dominant body plus 13 scraps and "the largest single solid is the floor,
not the sum". Dispatching by body cannot lower that floor. G13 asked whether
tiling *within* a solid can, and answered two questions that were open.

#### What G13 measured

Closed all-interior `m x m x m` grids at `cc=10, t=1.5`, 1,632 -> 99,840 faces.

* **Cost is mildly superlinear and worsening.** Overall log-log slope **1.135**
  (under the 1.15 bar this gate set in advance), but local slopes climb with
  scale — 1.044 at 12 k->25 k faces, 1.193 at 42 k->67 k, **1.273** at
  67 k->100 k — and ms/face rises monotonically 0.128 -> 0.230. Serially, 8
  tiles cost 6.494 s against the whole solid's 7.054 s, an 8 % saving. **Plan on
  tiling being worth about `W` and no more**; the win is parallelism, not a
  smaller exponent. This trend is also why G7's 0.17 ms/face (<=8.5 k faces) and
  the rehearsal's 1.11 ms/face (1 M faces) cannot be reconciled by assuming
  linearity.
* **Tiles reassemble with no sewing — but only with `unify_edges=False`.** With
  edge unification on, 3,632–11,232 seam edges come back as two distinct objects
  and the reassembled shell is full of holes: the edge pass concatenates
  collinear pairs *on* the tile boundary, so the two sides stop being the same
  `TShape`. With it off, identity is exact at every tile count and both scales
  tried — **0 free edges**, `BRepCheck_Analyzer` valid, volume preserved to
  ~1e-13 — so `BRep_Builder.Add` alone suffices and nothing sews the
  volume-scaling face set (docs/algorithm.md §6, §8).
* **`unify_edges=False` is 3.1–3.3x faster and merges exactly the same faces.**
  23.007 s -> 7.054 s at m=16, both producing 53,760 faces; 14.183 s -> 4.512 s
  at m=14, both 36,456.
* **Correction to docs/algorithm.md §9.** Its "edge merging is worth almost
  nothing" (4 edges of 81,816) was measured on the 80 mm ball and **does not
  hold at lattice scale**: at m=16 the edge pass takes 307,200 edges down to
  215,040, a 30 % reduction. Dropping it is therefore a *trade*, not a free win.
* **Tiling costs merged faces at tile seams** — +2.1 % to +11.4 % here — but that
  is an upper bound from a deliberately *generic* partition (by face centroid).
  The merge pairs are known by construction, one per surviving mid-strut
  interface, so bucketing by **strut** makes the loss zero.

#### Phase 1 — drop edge unification: **R1 FAILED, landed as enabling work only**

Projected `simplify` 18 m 39 s -> ~5 m 45 s, ~27 % off the run. **Delivered
nothing of the sort.** Both candidate forms were measured on `dense-lattice`
and `smoke-verified` on 2026-08-16, against R1's pre-set bars — net run time
must improve, file size must not grow more than ~20 %:

| `simplify` on `dense-lattice` | time | output |
|---|---|---|
| combined, one call (baseline) | 13.21 s, 16.18 s | 67,898 edges, 52.80 MB |
| edge pass dropped (the proposal) | 9.45 s | 94,476 edges, 71.29 MB |
| split: faces, then edges alone (the fallback) | 13.87 s, 13.94 s | 67,898 edges, 52.80 MB |

**The proposal fails both bars.** Dropping the edge pass saves 3.8 s in
`simplify` and hands all of it to `validate` (6.24 -> 8.21 s) and `export`
(6.25 -> 10.50 s) — both scale with edge count, which the projection assumed
away — for a **35 % larger file** and no net run-time change (57.28 ->
57.57 s). The 80 mm ball behaves the same way (+43 % file).

**The fallback is neutral, not cheap.** Its premise — that the edge pass would
be far cheaper over an already-face-merged solid — is disproved: it costs
~4.4 s either way.

Two things went wrong in the projection, both worth remembering. First, the
3.1–3.3x speedup came from G13's synthetic all-planar grid and is only ~1.4x on
a real trimmed solid. Second, `simplify`'s "baseline" of 16.18 s was a
cold-cache first run; the true baseline is ~13.2 s, and the repeat measurement
that caught this is the only reason the fallback was not recorded as a win.

**What landed, and why.** The split (`_unify_one` calls face merging with
`unify_edges=False`, then edge merging alone) is kept — not as a speed win, but
because it is a precondition for Phase 3: the edge pass rewrites edges on a tile
boundary, so tiles only reassemble by shared topology if the face merge runs
with it off. It also degrades better — a throwing edge pass no longer discards a
completed face merge. `test/test_pipeline.py` pins its B-rep as identical to the
combined call's, faces and edges alike, and `tools/e2e.py` passes with both
golden samples at 0 mm³ and byte-identical file sizes.

**Do not retry dropping the edge pass.** The measurement above is the disproof.

#### Phase 2 — build the interior pre-merged: **DONE, and it delivered**

Implemented 2026-08-16 in `interior.py` (`_merge_lateral_pairs`,
`_splice_lateral`, `MergedLateral`) and documented normatively in
[algorithm.md](algorithm.md) §6, §9 and §12. Measured on the two committed
scenarios, with the output unchanged throughout — same face and edge counts,
both golden samples at 0 mm³, `tools/e2e.py` green:

| `dense-lattice` | before | after |
|---|---|---|
| interior faces | 14,256 | **9,516** (−33 %) |
| interior edges | 30,900 | **21,420** (−31 %) |
| `instance` | 1.60 s | 1.06 s |
| `assemble` | 0.97 s | 0.74 s |
| `simplify` | 13.21 s | **9.93 s** (−25 %) |
| `validate` | 6.24 s | **5.44 s** (−13 %) |
| `export` | 6.25 s | **5.32 s** (−15 %) |
| total | 55.6 s | **44.8 s** (−19 %) |

The reduction is 33 % rather than the projected 50 % because only
interior↔interior struts merge, and `dense-lattice` has 594 interior nodes
against 968 boundary ones. The rehearsal, with 29,375 interior nodes against
19,552 boundary, does land close to the ceiling: **44.8 %**, measured below.

**One bug worth recording, because the class of it will recur.** The merge
condition "is the node across this cap also interior?" was first read from the
node-position cache, which `position()` grows with any neighbour it is asked
about — including boundary nodes reached through the cap correspondence. Interior
junctions then merged onto neighbours that were never built, leaving 366
unmatched edges on the 80 mm ball. Every unit test still passed: they use
synthetic all-interior grids, where the polluted set and the correct one are the
same. It was `weld.assemble`'s every-edge-twice proof that caught it, in `e2e`.
The lesson is the same one G8 left about scale — a set that is *derived* rather
than *stated* will eventually be derived from the wrong thing, and only a test
with a real boundary can tell.

*Blocking risk R4 — merged-loop correctness — **CLEARED**.* The spliced wire has
to be planar, simple and correctly wound for every `(cc, t)`.
`test_junction.py::test_merged_lateral_faces_are_sound_across_the_parameter_range`
asserts exactly that over G1's own sweep — all 14 valid pairs, every merged loop
coplanar, Newell normal matching the template's outward normal, and area equal to
the sum of the two halves to a relative 1e-12 — alongside watertightness, the
`N x volume(J)` identity, and both golden samples at 0 mm³. The fallback is built
as designed and is per strut family, not global: `_splice_lateral` returns
``None`` on any of those checks and that family keeps its two half-faces, so a
pathological case costs the face count this exists to reduce and never
correctness (docs/algorithm.md §11).

##### The original proposal



`interior.py` emits 4 lateral faces per half-strut and `simplify` merges them
back across every mid-strut interface. Emit the merged full-strut lateral face
directly instead: at template-build time, splice junction A's loop with the
neighbour's matching loop and drop the shared cap edge — a fixed pattern of
`(node-local, neighbour-local)` vertex indices computed once, the same shape of
precomputation `interior._pair_caps` already does — with the collinear pairs
already dropped so the wire is minimal. Merge only where the cap is in
`interfaces` **and** both nodes are interior; at a boundary interface the
neighbour's face comes from a boolean.

This is the strongest item on the list, and Phase 1's failure does not touch it:
it removes work rather than reordering or parallelising it. It halves interior
faces *before they are built* (G13: unification achieves exactly 99,840 ->
53,760), so it shrinks `instance`, `assemble`, `simplify`, `validate`, `export`,
file size and peak memory together — and because the wire is ours to choose, it
delivers the edge reduction by construction, at no kernel cost.

##### Measured at rehearsal scale: **-6.5 %**, not the -19 % `dense-lattice` showed

The projection this once carried (47.1 -> ~37 min) is withdrawn, not restated:
it was anchored on the 2026-08-15 profile, which ran with the `assemble`
watertightness gate bypassed and before the G9-G12 defect family was fixed.
Phase 2 was instead measured directly, as a **controlled pair run back to back
on the same machine on 2026-08-17** — `82adbb1` against `928cc57`, same input,
same `--cores 6 --ram 20`.

| Stage | before (`82adbb1`) | after (Phase 2) | delta |
|---|---|---|---|
| classify | 2 m 06.3 s | 2 m 05.9 s | -0.3 % |
| boundary | 11 m 54.6 s | 12 m 01.3 s | +0.9 % |
| connect | 11.3 s | 11.3 s | -0.6 % |
| stitch | 10 m 29.1 s | 10 m 23.0 s | -1.0 % |
| **instance** | 1 m 14.2 s | **41.7 s** | **-43.8 %** |
| **assemble** | 30.9 s | **21.0 s** | **-31.9 %** |
| **simplify** | 20 m 41.3 s | **18 m 07.9 s** | **-12.4 %** |
| **validate** | 4 m 02.8 s | **3 m 44.9 s** | **-7.4 %** |
| export | 3 m 49.3 s | 3 m 52.5 s | +1.4 % |
| **total** | **55 m 17.9 s** | **51 m 43.3 s** | **-6.5 %** |
| peak tree RSS | 19,827 MB | 19,291 MB | -2.7 % |
| `simplify` peak RSS | 9,956 MB | **7,339 MB** | **-26.3 %** |

**The five untouched stages agree to within 1 %**, which is what makes the
touched-stage deltas readable at all — and is the control the 2026-08-14/15
pair lacked, where untouched stages swung 25-36 % from machine load alone.

**The output is identical**: 584,028 faces, 2,517,881 edges, 14 solids,
330,354.002 mm³, 2.00 GB. Interior faces fell 705,000 -> **389,492 (-44.8 %)**,
close to the 50 % ceiling this optimization can reach, and far above
`dense-lattice`'s 33 % — exactly as predicted from the interior/boundary node
ratio.

**Two findings worth carrying forward.**

*The gain does not transfer between parts.* `dense-lattice` measured -19 %; the
rehearsal measures -6.5 %, on a **larger** interior reduction. The reason is
that this part's run is dominated by stages Phase 2 does not touch — `boundary`
(12 m) and `stitch` (10 m 23 s) are 43 % of it between them — so halving the
interior cannot move most of the clock.

*`simplify` scales with its output, not its input.* Unification's input fell
~31 % (1,006,505 -> ~691,000 faces) and the stage fell only 12.4 %, while its
output face count was unchanged at 584,028 by construction. Its memory and I/O
fell much more (-26 % peak RSS, -21 % bytes). **What this changes is the
arithmetic of any face-count extrapolation, not Phase 3's case**, and an earlier
revision of this section got that wrong. Removing input work — Phase 2's lever —
buys less than face counts suggest. Tiling is a *parallelism* lever and divides
whatever work remains, output included, so it is unaffected by which of the two
dominates. If anything it leaves Phase 3 worth **more** in absolute terms than
projected, because the base it attacks measured 18 m 07.9 s rather than the
~11 m 45 s that was extrapolated.

#### Phase 3 — sub-body tiling of `simplify`

Partition the solid's faces into tiles, unify each across the shared
`WorkerPool`, reassemble with `BRep_Builder.Add`, then apply the existing
`pipeline._check_unify_result` guards to the reassembled solid *plus* the
every-edge-twice proof from `weld.assemble` — a strictly stronger gate than what
runs today. Tile **by strut, not by centroid**: tag each face with its owning
node during the interior build and bucket by lattice-index block, reusing
`weld._tile_edge_length`; boundary faces fall back to centroid bucketing.

**Its gate condition is met.** "Do this only if `simplify` is still a top-3 cost
after Phase 2" — after Phase 2 it is the **largest single stage**, 18 m 07.9 s
of a 51 m 43.3 s run (35 %), and one of only three stages left that are pinned
at 0.98 cores while five of six sit idle (profiling-reports.md).

Only the *face* merge tiles; the edge pass has to stay global, since
concatenating collinear pairs is exactly what breaks tile identity. Measured on
`dense-lattice` the two split roughly 68 % / 32 % (9.45 s face merge, ~4.4 s
edges), so tiling can only ever attack about two thirds of the stage and the
edge pass becomes the new floor. Carrying that split onto the rehearsal's
measured 18 m 07.9 s gives ~12 m 20 s of tileable face merge against ~5 m 48 s
of global edge pass; at the ~5x the shared pool achieves elsewhere that is
`simplify` **~18 m 08 s -> ~8 m 20 s, worth ~9-10 min**, or ~19 % of the run.
Treat it as an upper bound and quote any total against the measured 51 m 43 s,
not against the withdrawn 47.1-minute profile.

*Blocking risk R2 — **CLEARED** by gate G14, 2026-08-16.* G13 tested only planar
instanced faces, leaving `TShape` identity across a tile seam unproven for the
trimmed, curved faces a boolean produces — which Phase 2 made the *dominant*
part of what `simplify` still sees (on `dense-lattice`, 15,718 boundary-derived
faces of 25,234). `tools/prototypes/g14_tiled_unify_trimmed.py` tests both a
closed instanced-grid ∩ sphere solid (1,804 faces, 76 curved) and a shell sewn
from a kept run's real trimmed boundary pieces (2,342 faces, 176 curved):
**0 free edges at 8, 27 and 64 tiles on both**, valid, volume/area preserved.
So tiling need not fall back to interior-only. Two caveats are recorded in
`RESULTS.md` G14: the gate's inputs merge only lightly (Phase 2 having already
merged the interior), so "identity survives *heavy* merging" still rests on
G13's planar case; and its first volume bar (1e-12, carried over from planar
measurements) wrongly failed a correct tiling at 5.31e-10 drift before being
reset to `pipeline.UNIFY_VOLUME_TOL`.

*Risk R3 — merge loss.* Strut-aware bucketing must measure ~0 % loss against a
whole-solid unification at m=14 and m=16. *Fallback:* accept the centroid
partition and a few percent more faces (a §11-acceptable larger output), or skip
tiling.

*Risk R5 — memory and IPC.* Verify with `tools/profile_run.py` and
`tools/profile_report.py` on the rehearsal. **Bars, against the measured
post-Phase-2 profile rather than the withdrawn one:** whole-tree peak RSS at or
below 19,291 MB, and `simplify` core-equivalents above 3.0 against today's 0.98.
Note that `simplify`'s own peak is 7,339 MB, well under the run's 19,291 MB peak
— which `export` sets — so there is real headroom to spend here before the run's
high-water mark moves at all. *Fallback:* cap the number of concurrent tiles.

*Risk R6 — the G8 trap.* Any new per-face bookkeeping must use
`TopTools_IndexedMapOfShape`, never Python lists tested with `.IsSame()`. Not
detectable at gate scale; only the rehearsal shows it.

#### Common gate for every phase

`python -m pytest test -q`, `python tools/e2e.py`, both golden samples at 0 mm³,
then one full rehearsal under `profile_run.py`. On landing, update
docs/algorithm.md §9 and §12, this section, and docs/testing.md's performance
notes — including the §9 correction above.
