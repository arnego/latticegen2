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
| --cores | int | optional | count | 1 - 128 | detected | Physical CPU cores available; `--workers` is derived from it as `min(cores-1, 8)`. |
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

### Scale rehearsal at the stated future size has not been run

**What's missing:** the stated future workload is parts ~8× larger in volume at
`cc=5, t=1` — roughly 64× more lattice cells than `dense-lattice`. Nothing at that
size has actually been generated. The architecture was chosen partly on projections
for it, and the interior path was measured standalone at that scale (64,000
junctions, 1.55 M faces, 199 s — `tools/prototypes/g2_instancing_join.py`), but a
full end-to-end run with classification, boundary trimming and STEP export at that
size has not happened.

**Why it matters:** STEP *writing* is expected to become the dominant cost at
that scale, producing a multi-GB file with ~3 M faces. If that
holds, the bottleneck moves outside this tool — into whether SolidWorks/Catia can
usefully import such a file — and that is a user-level decision about the output
contract, not something to change unilaterally.

**How to do it:** scale `test/test-cylinder.STEP` 2× linearly, run at `cc=5, t=1`,
and record the per-stage table, peak RSS, output file size and write time. Compare
against the proposal's projections. Report the file-size finding to the user rather
than acting on it.

**Watch the `simplify` stage there specifically.** Same-domain unification
(docs/algorithm.md §9) runs at roughly 0.24 ms/face, so ~3 M faces implies on the
order of 12 minutes. At `dense-lattice` scale it more than pays for itself by
halving the file that export and the round-trip check then have to handle, but
that trade has only been measured at ~30 k faces. It unifies each solid
independently, so it parallelises across solids if it needs to.

### Boundary stitching is the remaining superlinear step

**What it is:** the interior is built by indexed instancing and is linear, but
trimmed boundary junctions are joined with `BRepBuilderAPI_Sewing`, which was
measured at 14.9 s per 1,000 shells and growing worse than linearly
(docs/algorithm.md §6). Boundary count grows with surface area rather than volume,
so this is not currently a bottleneck — `dense-lattice` spends 12.5 s there — but it
is the term that will bite first at the 64× size above.

**How to fix it if needed:** the pairing is not arbitrary. Each trimmed junction
knows which cap interfaces it kept, so the stitching could be done per-interface
against known partners instead of by global geometric search, or the interior
shell's own index could be extended to absorb boundary pieces whose caps came
through the trim intact (which is the common case — see docs/algorithm.md §5.3).
Measure before building either: this only matters once the boundary count is in the
tens of thousands.

