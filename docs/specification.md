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
- **Julia version:** 1.10
- **Offline requirement:** Package must run with **zero network access**. List anything that currently requires network (package managers, license checks, telemetry) so it can be eliminated or vendored.
- **Packaging form:** Julia project invoked as `julia --project=. src/main.jl <args>`, with a thin `latticegen2.bat` (Windows) / `latticegen2.sh` (Linux) wrapper script provided for convenience (implemented: [latticegen2.bat](../latticegen2.bat), [latticegen2.sh](../latticegen2.sh)).
Using PackageCompiler.jl a standalone-executable distribution for windows must be created. Build tooling implemented: [`tools/build_app.jl`](../tools/build_app.jl) (isolated build-time dependency in [`tools/build/Project.toml`](../tools/build/Project.toml), entry point [`src/app_entry.jl`](../src/app_entry.jl)) — see README.md "Building a standalone executable". The actual compiled `.exe` has not yet been built or tested; running `tools/build_app.jl` takes ~15-45+ minutes and was deferred by user decision (2026-08-08) to keep this session's footprint small — do so before relying on this for distribution.
- **Target machine specs / limits:** Main development system: 32 GB RAM, 6 cure CPU, Nvidia RTX 3080 GPU , disk space for intermediate meshes
RAM and CPU cores should optinally be provided as input parameters, and optimization parameters should be determined automatically thereafter. If RAM and CPU cores are not provided, the optimization parameters must be provided explicitly instead.
- **Allowed third-party libraries:** Must be compatible with the target OS/arch. License text must be obtained and put into /licenses folder, and @/licenses/libraries.md must be updated with the cross reference between the library used and the corresponding license text file valid for that library.
- **License constraints:** TBD

---

## 3. Command-Line Interface

Exact invocation the human will type. This is the user-facing surface.

For each parameter, specify: **name, type, units, valid range, default, required?**

| Flag | Type | Required | Units | Range | Default | Description |
|------|------|----------|-------|-------|---------|-------------|
| -i --input | path | required | NA | NA | NA | Path to STEP file defining the lattice bounds |
| -o --output | path | optional | NA | NA | NA | Path and name of the output .step file, otherwise generated (e.g. input_file_name-cc5t1|
| -cc | float | required | mm | 0.4 - 50 | NA | Distance between the bottom nodes of two adjecent cells |
| -t | float | required  | mm | 0.4 - 20 | NA | Side length of the diamond rod profile |
| -bg --background | flag | optional | NA | NA | disabled | Run worker processes at below-normal priority to reduce desktop impact |
| -v --verbose | flag | optional | NA | NA | disabled | Enable verbose console diagnostics while always writing a full `.log` file. |
| --cores | int | optional | count | 1 - 128 | NA | Physical CPU cores available for auto-tuning optimization parameters (§2). Mutually exclusive with `--workers`/`--tile-cells`. |
| --ram | int | optional | GB | 1 - 1024 | NA | RAM budget available for auto-tuning optimization parameters (§2). Mutually exclusive with `--workers`/`--tile-cells`. |
| --workers | int | required if `--cores`/`--ram` not given, forbidden otherwise | count | 1 - 128 | NA | Explicit number of parallel worker processes for the tiling stage (§2). |
| --tile-cells | int | required if `--cores`/`--ram` not given, forbidden otherwise | count | 2 - 64 | NA | Explicit lattice-cell edge length (cells per tile axis) for the tiling stage (§2). |

**Exit:** 

Upon success the script shall produce an end of run summary report in the .log file and to console independent of the -verbose flag. This shall include: 
 - The runs input parameters 
 - Date and time of run start
 - Duration from start to completion
 - Run characteristics (number of tiles, number of parallel threads per stage of the generation procedure, etc)
 - Maximum memory usage
 - Path to output .step file 

Upon failure, the script shall output a human readable reason for the failure, e.g.: parameter bounds exceeded, issues with input geometry, issues with resarouces from the run-time system, write or read access issues, etc.

Upon cancellation by the user (Ctrl+C), the script shall shut down gracefully rather than terminate abruptly: worker processes are stopped in an orderly way (force-stopped only if they do not respond within a short grace period), a single human readable `CANCELLED` line is written to console and `.log` file, the temporary folder is left in place for analysis, and the exit code is 130. See docs/algorithm.md §9.1.
  
**Logging:**

A log file should be produced every run with the same name as the output file which is generated from the input file or provided by the -o flag. The log file should end with `.log` and should not include .step (that is only the last name for the geometry file.  

---

## 4. Geometry Domain Specification

### 4.1 Lattice unit cell type
The base geometry for the lattice is a strut-based uniform grid forming cube-like cells standing on its tip. The struts form the boundries of the cells along each edge. The struts have the profile of a square standing on one corner, like a diamond. In other words; for each strut the square profile is oriented so one diagonal lies in the vertical plane containing the strut axis and the Z-axis, and the other diagonal is horizontal. See docs/algorithm.md §2.2 for the exact frame construction. Verified: [`profile_vertices`](../src/lattice.jl) builds each profile from `u_k = normalize(cross((0,0,1), e_k))` (horizontal) and `v_k = cross(e_k, u_k)` (vertical-plane) exactly per docs/algorithm.md §3.1.

The dimentions of the square profile is defined by input parameter `t`. Upon inspection of the end result, the struts are reclined from the normal axis (Z-axis) in degrees by the following calculation in numpy: np.degrees(np.arcsin(np.sqrt(2/3)))
Make sure to sure the native exact definition in native Jula language. It should be close to 55 degrees (but not exactly).
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
including the exact lattice math, pipeline/classification/tiling diagrams, and the
detailed optimization strategy this section summarizes.

Since this involves computational geometry:
- Profile geometry generation routines to identify bottlenecks
- Consider vectorization or parallelization for lattice tiling
- Cache expensive calculations (e.g., basis functions for triply-periodic surfaces)
- If caching to disk is used, put the files in a temporary folder `temp/<date><time>` where the output file is generated to. Clean up after a sucessful run. Leave for error analysis if the fun failes.

---

## 5. STEP Output Requirements
- **Clean up** Never produce a floating body with volume < t³ mm³ (i.e. a cube of side `t`).
  This rule targets **floating (disconnected) bodies only** — a solid is only ever
  discarded once it is verified to have no geometric connection to the rest of the
  output, never merely because its own volume happens to fall below t³. A
  sub-threshold solid that is still connected to other geometry (e.g. a junction
  fragment produced by an upstream boolean operation that hasn't fully converged) is
  kept, and if that connectivity cannot be resolved automatically the run fails rather
  than guessing (exit 4). See docs/algorithm.md §8 for the exact resolve-then-classify
  algorithm and docs/algorithm.md §11.2 for the investigation that found the earlier,
  unconditional "volume < t³ → delete" reading of this rule was deleting connected
  junction material, not floating debris.
  
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

| Scenario | Parameters | Expected result |
|----------|-----------|------------------|
| smoke-fast | -i test/80mm-test-ball.step -cc 20 -t 4 -bg | generation < 10 minutes, quick test not applicable for output geometry verification |
| smoke-verified | -i test/80mm-test-ball.step -cc 20 -t 4 -bg | valid STEP, generation < 20 minutes, matching golden sample test/80mm-test-ball-cc20t4-golden-sample.step — implemented as `tools/e2e.jl`'s `smoke-verified e2e` testset. (Params corrected from an earlier -cc 10 -t 2, which didn't match the golden file's own cc20t4 name; test/80mm-test-ball-cc20t4-golden-sample.step was generated with cc=20/t=4, so the scenario's params were changed to match it. Golden file renamed from test/80mm-test-ball-cc20t4.step to test/80mm-test-ball-cc20t4-golden-sample.step 2026-08-09.) |
| dense-lattice | -i test/test-cylinder.STEP -cc 10 -t 1.5 --cores 6 --ram 20 -bg | valid STEP, no self-intersections, matching golden sample test/test-cylinder-cc10t1.5-golden-sample.step, generation < 60 minutes — harness implemented as `tools/e2e.jl`'s `dense-lattice e2e` testset, but it self-skips because the golden sample does not exist yet. Originally specified at -cc 5 -t 1, the one attempt to generate that denser golden sample ran for hours and was manually terminated (not a crash — the run's own log recorded exactly what was happening throughout; see docs/algorithm.md §11.2 for the full investigation): an auto-tuned tile size well past the fuse-time performance knee left assembly with 11k+ unfused solids, and the (now-fixed) unconditional sub-threshold cleanup rule was deleting connected junction material one solid at a time for over an hour. Root cause diagnosed and fixed (docs/algorithm.md §7.1, §6.3, §8); the scenario's params were changed 2026-08-09 to -cc 10 -t 1.5 (a less dense lattice capable of finishing within a reasonable time) with a 60-minute budget. Regenerating the golden sample at these params is a follow-up step, not automated here (committing a golden sample is a decision for whoever reviews the regenerated run's output). |
| invalid-input | -i test/80mm-test-ball.step -cc 5 -t 4 --workers 1 --tile-cells 4 (strut size `t` >= cell edge `a=cc/√2`) | exits nonzero (exit 2), no `.step` or `.log` file written — implemented as `tools/e2e.jl`'s `invalid-input e2e` testset |

### 6.2 Automated pass/fail checks
For every scenario the harness must verify, without human intervention:
- [ ] Process exits with expected console output
- [ ] STEP file is written and non-empty
- [ ] STEP file parses back successfully (round-trip read with the same or an
      independent library)
- [ ] Geometry is a valid closed manifold solid (no open edges / non-manifold edges)
- [ ] No self-intersections
- [ ] Bounding box of output matches requested `--input` within tolerance
- [ ] Runtime stays under an agreed performance budget for each scenario size: `smoke-fast` < 10 minutes and `smoke-verified` < 20 minutes (both measured against the 80mm test ball; enforced in `tools/e2e.jl`'s `BUDGET_SECONDS`/`VERIFIED_BUDGET_SECONDS`). `dense-lattice`'s budget remains TBD until a golden sample exists to establish a baseline against (§6.1).
- [ ] If golden sample is defined, check similarity of geometriers by running a geometry check script that inspects if the output geometry is inhibiting the same volume as the golden sample (e.g. subtraction either way should leave near zero volume). 


### 6.3 How verification runs offline
- Verification runs only in the dev/CI environment.
- Test runner: Test scripts and assets are separated into separate `tools/` folder.
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

### `filter_floating!` has no time budget (progress logging added 2026-08-09)

**What's broken:** `filter_floating!` (`src/pipeline.jl`, the export-stage
floating-body-only cleanup gate — docs/algorithm.md §8) has **no `max_seconds` circuit
breaker**, unlike `balanced_fuse!`, which has both a budget and progress logging
(docs/algorithm.md §6.5). Its per-component `overlap_components` resolution loop calls
`fuse_all` on every ambiguous multi-member component with nothing bounding how long
that can take in aggregate. **Progress logging (item 3 below) was added on 2026-08-09**
— one line up front (solid/component/ambiguous-component counts), one line every 10s
of wall time while resolving, one summary line when done, all via `log_line(rl, ...)`
so they land in the `.log` file always and on console when `-v`/`rl.verbose` — but the
loop still has no upper bound on total time (items 1, 2 below remain open).

**How it was found:** re-verifying the `test-cylinder-cc5t1` scenario
(`-i test/test-cylinder.STEP -cc 5 -t 1 --cores 6 --ram 20 -bg`) after the fixes in
docs/algorithm.md §11.2. `tile_stage` completed in 11 minutes (vs. 1h 38m before,
confirming those fixes work). Assembly then left roughly 239–325 solids unresolved per
merge group — curved-surface boundary-trim fragments that don't share exactly
coincident faces after trimming, a known OCCT robustness limit (docs/algorithm.md §6.5
"Fuse-failure fallback"), not a new bug. `filter_floating!` then ran **silently for 56+
minutes with zero log output** before being manually killed at ~2.5h total elapsed.
This reproduces, in new code, the exact class of problem (an invisible, unbounded
long-running cleanup step) the whole session's fix work was aimed at eliminating.

**How to fix it:**
1. Add a `max_seconds` parameter to `filter_floating!`, checked before each
   ambiguous component's `fuse_all` attempt — mirror `balanced_fuse!`'s pattern
   exactly (check-then-attempt, not mid-call preemption, since OCCT calls can't be
   interrupted once started).
2. On budget exhaustion, any sub-threshold member of a component that hasn't yet been
   resolved must be treated as **unresolved**, exactly like a fuse-failure case already
   is — hard-fail (`ProcessingError`, exit 4) per the existing policy of never
   guessing whether an unresolved connected fragment is safe to delete. Do not
   silently keep or drop it just because the clock ran out.
3. ~~Add per-component progress logging~~ **Done (2026-08-09):** `filter_floating!`
   now takes a `progress_seconds::Real=10.0` keyword and logs an up-front summary line,
   a progress line every `progress_seconds` of wall time while resolving ambiguous
   (multi-member) components, and a completion summary line — see the updated
   docstring in `src/pipeline.jl`. Deliberately scoped to logging only, per explicit
   user request; items 1, 2, 4, 5 below remain open follow-up work, not implied by
   this change.
4. Add a regression test exercising a slow/many-component scenario (an injectable
   `fuse_fn` — already supported for the hard-fail test in `test/test_cleanup.jl` —
   can simulate this without needing genuinely slow geometry).
5. Re-run the full `test-cylinder-cc5t1` scenario end-to-end afterward and confirm the
   export stage now completes (or hard-fails) within a bounded, logged time, then
   compare the whole run against the pass-criteria table used for this fix (tile_cells
   ≤ 8, full-interior > 0, no un-budgeted silent stage, tile_stage < 20 min).

**Where things stand:** progress logging (item 3) implemented and covered by a minimal
unit test in `test/test_cleanup.jl` (`filter_floating! logs progress for ambiguous
components`, using the existing injectable-`fuse_fn` pattern with a synthetic 2-member
component — no genuinely slow geometry needed). The `max_seconds` circuit breaker
(items 1, 2) and the full `test-cylinder-cc5t1` end-to-end re-run (item 5) are still
not started. The verification run that originally found this issue was killed by
explicit user choice rather than left running indefinitely or fixed on the spot. ~944 MB
of diagnostic temp files from that killed run may still be on disk in whatever
scratchpad directory that session used, under `cylinder_verify/temp/20260809-002344`,
if useful for follow-up analysis — treat that path as ephemeral and verify it still
exists before relying on it.
