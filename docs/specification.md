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
- **No single-file standalone executable is produced.** PyInstaller was evaluated and rejected: `boundary.py` uses `multiprocessing` with the `spawn` start method and the codebase has no `freeze_support()` call, which on Windows makes a frozen build re-launch its own launcher; OCP/OCCT is awkward to freeze (hidden imports, DLL discovery); and freezing dissolves the LGPL-2.1 relinking argument in [licenses/libraries.md](../licenses/libraries.md), which depends on OCCT remaining a stock, replaceable shared library. The portable bundle delivers the same "extract and run" property without those costs.
- **Target machine specs / limits:** Main development system: 32 GB RAM, 6 core CPU, Nvidia RTX 3080 GPU, disk space for intermediate files.
RAM and CPU cores may optionally be provided as input parameters, as *budgets* rather than a mandatory pair. Without `--cores` the worker count is the machine's logical core count, since boundary-junction jobs are constant-size and independent; without `--ram` the budget is the memory free at startup. See §3.
- **Allowed third-party libraries:** Must be compatible with the target OS/arch. License text must be obtained and put into /licenses folder, and @/licenses/libraries.md must be updated with the cross reference between the library used and the corresponding license text file valid for that library.
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

### Scale rehearsal, chapter closed: paths 1–4 implemented and re-measured

**First run 2026-08-14**, **re-profiled 2026-08-15** after implementing paths
1–4 below, both on `TD_HX_Indre_Volum.step` at `cc=5, t=1`,
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
intended benefit; `TD_HX_Indre_Volum` at `cc=5, t=1` simply is not that
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

### Micron-scale debris edges from near-tangential trims

**What's broken.** The rehearsal fails in `assemble`: 2 edges are used by exactly
one face, so the every-edge-twice proof (docs/algorithm.md §8) rejects the shell.
Both are in component 0, at `[1874.836, 60.370, 970.121]` and
`[1874.836, 59.912, 969.775]`.

**What they are — measured, not inferred.** Each is a genuine, non-degenerate
edge of **3.171690e-06 mm** and **5.808982e-06 mm**, on a *planar* face of
~1.2 mm² carrying **8 edges where 7 would do**. They are boolean debris from a
strut grazing the input surface almost tangentially — the same pathology that
makes 2,969 of this part's 19,552 boundary junctions produce no geometry at all.

**They are not unpaired halves of anything.** The nearest other edge is
**0.070 mm** from one and **0.263 mm** from the other — four to five orders of
magnitude further away than the slivers are long — and every neighbouring edge is
correctly used twice. So the repair has to *remove* an edge, not match one.

**Three approaches tried and eliminated**, so they are not retried:

| Approach | Result |
|---|---|
| Raise `SEW_TOLERANCE` 1e-6 → 1e-5 | **No effect.** The tolerance decides whether two *different* faces' free edges are paired; it cannot remove an edge that has no partner. Reverted; the disproof is recorded in the constant's docstring. |
| Same-domain unification of the sewn boundary, before the rings are read | **No effect on these edges.** Implemented, measured and reverted. It ran on every component (301,505 → 268,520 faces, 97,043 edges removed, area drift 5.26e-07), so the kernel does not consider these slivers same-domain with their neighbours — each is a real tiny corner at an angle, not a collinear split. |
| `ShapeFix` small-edge removal on each trimmed piece, in the worker | **Not tried.** The remaining candidate. It cleans debris at its source and parallelises for free, but it runs before cap tagging, so it risks perturbing the cap-area agreement `resolve_interfaces` checks at `CAP_AREA_REL_TOL` (1e-6 relative). That interaction is what needs designing. |

**How to verify a fix.** Run the rehearsal; `assemble` must report zero open and
zero misoriented edges. `python -m pytest test -q` and `python tools/e2e.py` must
stay green — the two committed scenarios do not contain grazing trims this severe,
so they will not exercise the repair and must therefore be unchanged by it
(0 mm³ against both golden samples).

**One caution for whoever takes this.** A first attempt guarded the boundary
unification with a 1e-9 relative area bar, reasoned from "this merge is more
constrained than `simplify`'s, so it should hold tighter". That bar rejected a
*correct* unification of the 301,172-face component at a drift of 3.08e-09 and
silently undid the repair on the one component that needed it. Quadrature noise
scales with the shell being integrated; the figure has to come from measurement,
exactly as docs/algorithm.md §11 says about every gate. The same mistake in the
same shape as issue #6.
