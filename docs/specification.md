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
RAM and CPU cores may optionally be provided as input parameters. They are *hints* rather than a mandatory pair: without them the worker count is derived from the machine, since boundary-junction jobs are constant-size and independent.
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
| -bg --background | flag | optional | NA | NA | disabled | Run at below-normal priority to reduce desktop impact |
| -v --verbose | flag | optional | NA | NA | disabled | Enable verbose console diagnostics while always writing a full `.log` file. |
| --cores | int | optional | count | 1 - 128 | detected | Physical CPU cores available; `--workers` is derived from it as `min(cores, 8)` — one worker per core, since the master is blocked waiting on results for effectively the whole boundary stage. |
| --ram | float | optional | GB | 1 - 1024 | NA | Memory budget. Advisory: recorded in the run log. |
| --workers | int | optional | count | 1 - 128 | from `--cores`, else from the machine | Parallel worker processes for the boundary-junction stage. Overrides `--cores`. |

The optimization parameters are hints, not a mandatory pair: an explicit
`--workers` overrides everything, `--cores` derives it, and with neither given it
follows from the machine.

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
| smoke-fast | -i test/80mm-test-ball.step -cc 20 -t 4 -bg | generation < 10 minutes. **Measured: 8.6 s.** |
| smoke-verified | -i test/80mm-test-ball.step -cc 20 -t 4 -bg | valid STEP, generation < 20 minutes, matching golden sample test/80mm-test-ball-cc20t4-golden-sample.step. **Measured: 8.3 s, symmetric-difference volume 0.0000 mm³.** |
| dense-lattice | -i test/test-cylinder.STEP -cc 10 -t 1.5 --cores 6 --ram 20 -bg | valid STEP, no self-intersections, matching golden sample test/test-cylinder-cc10t1.5-golden-sample.step, generation < 10 minutes. **Measured: 61 s, symmetric-difference volume 0 mm³.** |
| invalid-input | -i test/80mm-test-ball.step -cc 5 -t 4 (strut size `t` >= cell edge `a=cc/√2`) | exits 2, no `.step` or `.log` file written, one human-readable reason line. **Passes.** |

### 6.2 Automated pass/fail checks
For every scenario the harness must verify, without human intervention:
- [x] Process exits with expected console output
- [x] STEP file is written and non-empty
- [x] STEP file parses back successfully (round-trip read)
- [x] Geometry is a valid closed manifold solid (no open edges / non-manifold edges)
- [x] **Geometry passes OCCT's exact B-rep validity check** (`BRepCheck_Analyzer`) —
      an exact test on the B-rep itself, not an inference from a tessellation
- [x] No self-intersections
- [x] **No generated material lies outside the input body** (boolean cut of output
      against input leaves ~zero volume) — a direct check of §1's "fits exactly
      within the user's boundary geometry", independent of any golden sample
- [x] Bounding box of output matches requested `--input` within tolerance
- [x] Runtime stays under an agreed performance budget: `smoke-fast` and
      `dense-lattice` < 10 minutes, `smoke-verified` < 20 minutes
- [x] If a golden sample is defined, check similarity of geometries by subtracting
      candidate and golden both ways; the larger remainder must be near zero


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

### Scale rehearsal: run end to end, profiled — one defect outstanding

**Run on 2026-08-14**, `TD_HX_Indre_Volum.step` at `cc=5, t=1`, `--cores 6 --ram 20 -bg`,
on the 6-core / 32 GB development workstation. This is the first time the part has
reached the end of the pipeline: the previous attempt (2026-08-12) died in `sew`
after 5 h 04 m, and the 2026-08-13 re-run was stopped at assembly input.

Two of this chapter's three items are now **verified at the scale that motivated
them** and have been retired into [algorithm.md](algorithm.md): the boundary-sew
tiling (§8) and the disagreeing-cap fuse (§7.1). What replaces them is the
measured result of the run and the profile taken during it.

**The run is not yet a clean pass.** It fails in `assemble` on 2 edges out of
~1.4 M — see "Micron-scale debris edges" below. The timings in this section come
from a diagnostic run that was allowed past that gate on purpose, so every later
stage could be measured; the geometry it produced is therefore *not* a shippable
output, but the stage costs are representative because the same work is done
either way.

#### Per-stage cost and resource profile

Wall time from the run's `.log`; CPU, memory and I/O from
[`tools/profile_run.py`](../tools/profile_run.py), joined to the stage boundaries by
[`tools/profile_report.py`](../tools/profile_report.py). "Cores" is mean CPU over the
stage, where 1.00 is one core fully busy and 6.00 is the machine.

| Stage | Duration | Cores | RSS peak | Notes |
|---|---|---|---|---|
| template | 0.04 s | — | — | |
| import | 0.40 s | — | — | |
| tessellate | 2.4 s | 0.50 | 280 MB | 28,654 triangles, deviation 0.409 mm |
| classify | 1 m 54 s | 0.99 | 308 MB | 527,425 candidates → 29,375 interior, 19,552 boundary, 478,498 outside |
| boundary | 10 m 57 s | **5.55** | 2,128 MB | 21,955 pieces from 19,552 junctions; 2,969 produced no geometry |
| connect | 11.0 s | 0.98 | 880 MB | 122,180 interfaces; 1 cap cluster fused; 3 declined; 183 floating bodies dropped |
| stitch | 8 m 57 s | 1.16 | 2,713 MB | 21,694 pieces → 301,505 faces, 14 components, 35 tiles, 18,496 rings |
| instance | 1 m 07 s | 1.00 | 3,537 MB | 705,000 faces, 624,492 shared vertices |
| assemble | 29.1 s | 1.00 | 3,039 MB | |
| simplify | 17 m 17 s | 0.99 | 6,810 MB | 1,006,505 → 584,028 faces (42 % fewer), drift 1.6e-07 |
| validate | 2 m 59 s | 0.99 | 13,260 MB | |
| export | 6 m 42 s | 0.99 | **18,925 MB** | 2.00 GB written |
| verify | 22 m 29 s | 0.99 | 16,539 MB | 2.00 GB re-read |
| **total** | **73.1 min** | | **18.48 GB** | 14 solids, 330,354 mm³ |

Three things this profile establishes that no earlier measurement could, because
no run had ever reached these stages:

1. **The cost centre has moved past assembly.** `simplify + validate + export +
   verify` is **49.4 min of the 73.1**, i.e. 68 % of the run. Everything before
   `assemble` totals 23.7 min. The old bottleneck (`sew`, 94 % of a 5 h run) is
   gone and something else is now dominant.
2. **Only one stage uses the machine.** `boundary` runs at 5.55 of 6 cores.
   *Every other stage is single-threaded*, `stitch` included once its serial
   round 2 is averaged in.
3. **Memory is now the binding constraint, and it never was before.** RSS climbs
   monotonically from 3.0 GB at `assemble` to **18.9 GB** during `export`, and
   the system was left with only 3.0 GB available. This is priority #2 in
   CLAUDE.md, and at a moderately larger part this run would swap. The `--ram 20`
   advisory was very nearly consumed.

#### Boundary-sew tiling, measured against its own control

The tiling was re-run with tiling disabled (same input, same parameters, via a
wrapper that only raises `min_to_tile`), to separate its effect from everything
else:

| | Tiled | Untiled |
|---|---|---|
| `stitch` wall time | **8 m 57 s** | 20 m 27 s |
| `stitch` mean cores | 1.16 (peak 5.96) | 0.99 (peak 1.01) |
| Total to `assemble` | 25.1 min | 35.4 min |
| Pieces / faces / components / rings | 21,694 / 301,505 / 14 / 18,496 | **identical** |

So **2.25×** at 21,955 pieces, against G6's 1.43–1.45× at 4,000–8,000 — better
than the prototype measured, because production round 1 runs across worker
processes and G6's serial sum could not credit that. Both of the design's claims
hold at real scale: the saving is real, and the result is byte-for-byte the same
shell whichever route produced it. The remaining `stitch` cost is now round 2,
which is strictly serial — that is what holds the stage's mean at 1.16 cores
despite round 1 peaking at 5.96.

#### Optimization paths, ranked

Recoverable time assumes perfect parallelism across 6 cores where that is the
lever, which is an upper bound, not a promise. Ranked by (recoverable time) ×
(confidence).

1. **Parallelise `simplify` across solids — up to ~14 min, high confidence.**
   1,036 s at 0.99 cores over 14 independent solids.
   [algorithm.md](algorithm.md) §9 already unifies each solid separately and
   states it "parallelises across solids if it ever becomes the bottleneck at
   scale". It now is: 24 % of the run. The pattern is the same `.brep`
   round-trip `boundary.py` and `weld.py` already use, and the per-solid volume
   guard stays exact because the 1:1 mapping is unchanged. The 14 solids are
   very unequal, so expect well under 6× — but the largest single solid is the
   floor, not the sum.
2. **Make the round-trip gate cheaper — up to ~22 min, medium confidence.**
   `verify` is the single most expensive stage at 1,349 s, and
   [`round_trip_check`](../src/latticegen2/stepout.py) spends all of it
   re-parsing 2.00 GB of STEP to full B-rep **purely to count solids**. It is
   also CPU-bound (99 % CPU), not I/O-bound. The gate itself is mandatory
   (algorithm.md §9) and must not be dropped — a run that has not read back what
   it wrote has not established that it wrote it. But full B-rep reconstruction
   is far more than counting requires: a text scan for `MANIFOLD_SOLID_BREP`
   entities answers the same question at I/O speed. That is a *weaker* check,
   so this is a design decision about what the gate is for, not a free win, and
   it should be taken deliberately rather than silently.
3. **Parallelise `stitch` round 2 — up to ~7 min, medium confidence.** Round 1
   already parallelises; round 2 merges its outputs in one serial call. A tree
   reduction (pairwise merges across workers) would apply the same lever again.
   G6 warns the ceiling is low — round 2 does not shrink with tile count because
   its input face count is fixed — so measure before building.
4. **Parallelise `validate` across solids — up to ~2.5 min, high confidence.**
   180 s at 0.99 cores; `BRepCheck_Analyzer` per solid is independent. Small, but
   it is the same change as (1) and would ride along with it.
5. **Reduce peak memory — no wall-time win, but this is priority #2.** 18.9 GB
   peak with 3.0 GB headroom left. The growth is monotonic across
   `simplify → validate → export`, which suggests intermediate geometry is being
   retained after it is needed rather than any single stage being inherently
   huge. Worth *measuring* where (the profile CSV localises it to a stage; it
   does not say which object), because it is what decides whether a part
   moderately larger than this one runs at all.

**Recorded as at their floor**, so they are not re-investigated:

* **`export` (6 m 42 s)** — 99 % CPU writing 2.00 GB, so it is serialization
  cost, not disk. Note this *contradicts* algorithm.md §12's description of
  export as irreducible `O(faces)` I/O: it is irreducible in the sense that the
  file must be produced, but it is CPU-bound, so it is not I/O-bound and a
  faster disk would not help.
* **`boundary` (10 m 57 s)** — already at 5.55 of 6 cores. Nothing to recover
  without more cores.
* **`classify` (1 m 54 s)** and **`instance` (1 m 07 s)** — single-threaded, but
  together only 4 % of the run. `instance` in particular builds one shared-topology
  index, which is exactly the structure that does not parallelise cleanly.

#### Output size

2.00 GB, 584,028 faces, 14 solids. Nothing about this is wasteful: same-domain
unification runs successfully on every solid (`unmerged_solids: 0`) and removes
42 % of the faces before export. The size is what an exact B-rep lattice of this
density costs. Whether SolidWorks/Catia can usefully import a 2 GB STEP is a
question about the output contract, and a user-level decision rather than
something to change here.

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
