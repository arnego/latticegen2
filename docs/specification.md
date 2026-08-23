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
- **Packaging form:** invoked as `python src/main.py <args>`, with a thin `latticegen2.bat` (Windows) / `latticegen2.sh` (Linux) wrapper provided for convenience (implemented: [latticegen2.bat](../latticegen2.bat), [latticegen2.sh](../latticegen2.sh)). No install step is needed; `pip install .` optionally provides a `latticegen2` console script. **Started with no arguments at all, and where a display exists, the launcher opens the graphical front-end instead** (§3.1); it is built on the standard library's `tkinter`, so it adds no packaged dependency, though portable bundles consequently redistribute Tcl/Tk (see [licenses/LICENSES.md](../licenses/LICENSES.md)).
- **Distribution form:** per-platform **offline bundles**, published as GitHub release assets by [`.github/workflows/release.yml`](../.github/workflows/release.yml) on a `v*` tag, in two flavours for each of Windows and Linux x86-64:
  - *portable* — carries a relocatable CPython with every dependency installed. Extract and run: no Python on the target, no install step, no admin rights, no network.
  - *wheels* — source plus the dependency wheels and an `install` script, for a target that already has Python 3.11.

  Each release also publishes `SHA256SUMS.txt` for verification after transfer. Every asset is extracted and run end-to-end by a CI smoke gate before publication. Procedure: [release.md](release.md). Bundle contents are `git archive`-derived, so they contain committed files only, filtered by `.gitattributes`.
- **No single-file standalone executable is produced.** PyInstaller was evaluated and rejected: `boundary.py` uses `multiprocessing` with the `spawn` start method and the codebase has no `freeze_support()` call, which on Windows makes a frozen build re-launch its own launcher (the graphical front-end does not weaken this argument — it runs the pipeline as a *child* of `src/main.py`, precisely so `spawn`'s re-import contract is untouched); OCP/OCCT is awkward to freeze (hidden imports, DLL discovery); and freezing dissolves the LGPL-2.1 relinking argument in [licenses/LICENSES.md](../licenses/LICENSES.md), which depends on OCCT remaining a stock, replaceable shared library. The portable bundle delivers the same "extract and run" property without those costs.
- **Target machine specs / limits:** Main development system: 32 GB RAM, 6 core CPU, Nvidia RTX 3080 GPU, disk space for intermediate files.
CPU cores may optionally be provided as an input parameter, as a *budget*. Without `--cores` the worker count is the machine's logical core count, since boundary-junction jobs are constant-size and independent. See §3. (A companion `--ram` budget existed through v2.x; it was removed as accepted-but-unenforced dead weight — see §11.)
- **Allowed third-party libraries:** Must be compatible with the target OS/arch. License text must be obtained and put into /licenses folder, and @/licenses/LICENSES.md must be updated with the cross reference between the library used and the corresponding license text file valid for that library.
- **License constraints:** TBD

---

## 3. Command-Line Interface

Exact invocation the human will type. This is the user-facing surface.

For each parameter, specify: **name, type, units, valid range, default, required?**

| Flag | Type | Required | Units | Range | Default | Description |
|------|------|----------|-------|-------|---------|-------------|
| -i --input | path | required | NA | NA | NA | Path to STEP file defining the lattice bounds |
| -o --output | path | optional | NA | NA | `<input_stem>-lattice-cc<cc>t<t>.step` | Path and name of the output .step file. Must name a **file**, not a directory: `-o .\` and friends are rejected (exit 2) rather than turned into `.\.step`. `.step` is appended if absent. |
| -cc | float | required | mm | 0.4 - 50 | NA | Distance between the bottom nodes of two adjacent cells |
| -t | float | required  | mm | 0.4 - 20 | NA | Side length of the diamond rod profile. Must be smaller than the cell edge `a = cc/√2`; that is the only cross-constraint. |
| -v --verbose | flag | optional | NA | NA | disabled | Enable verbose console diagnostics while always writing a full `.log` file. |
| --gui | flag | optional | NA | NA | NA | Open the graphical front-end (§3.1) instead of running. Implied when the launcher is started with **no arguments at all** and a display is available. |
| --cores | int | optional | count | 1 - 128 | logical cores on the machine | Maximum CPU cores this run may use. One worker process per core, honoured exactly — the master needs none reserved for it, being blocked waiting on results for effectively the whole boundary stage. Since workers always run at below-normal priority, this exists to further protect the response time of the system for other tasks. |

`--cores` is an optional **budget** and resolves to a concrete figure either
way: an explicit value is honoured exactly, and an omitted one is taken from
the machine — its logical core count. Detection lives in
[`src/latticegen2/sysinfo.py`](../src/latticegen2/sysinfo.py).

**There used to be a second budget, `--ram`, removed 2026-08-17.** It was
accepted, range-checked and recorded in the run log next to the measured peak,
but nothing in the pipeline ever read it to change what a run did — see §11
for the full account of why, and what was done instead.

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

**Logging:**

A log file should be produced every run with the same name as the output file which is generated from the input file or provided by the -o flag. The log file should end with `.log` and should not include .step (that is only the last name for the geometry file.

---

## 3.1 Graphical front-end

**What it is.** Clicking `latticegen2.bat` / `latticegen2.sh` with no arguments
opens a small window that collects the same parameters §3 lists, runs the same
`src/main.py`, and shows the run's progress while it goes. **It adds no
capability the command line lacks**, and any argument on the command line keeps
today's behaviour exactly — so nothing scripted changes.

Implemented in [`src/latticegen2/gui/`](../src/latticegen2/gui/); the event
stream it reads is [`src/latticegen2/progress.py`](../src/latticegen2/progress.py)
and docs/algorithm.md §10.

**Why it exists.** A production part takes the better part of an hour
(docs/profiling-reports.md) and the tool could not say where it was. Peak memory
and the per-stage timings were only readable after the fact, from the `.log`.

**The window.** Fixed width, and it grows only while a run is in flight, so it
can sit in a corner of the screen.

* **Input** — a file, chosen with a browse dialog.
* **Output** — a **folder**. The filename is derived from the input stem and the
  parameters exactly as §3's default rule says, shown read-only beneath the
  field, and recomputed as the input, `cc` or `t` change.
* **cc**, **t**, **cores** — spin controls bounded by the same constants
  `cli.parse_args` enforces, defaulting `cores` to the machine's logical count.
* **verbose** — a tick box beside `cores`, controlling the log pane described
  below. It is the window's equivalent of `-v`, and like `-v` it changes only
  what is *shown*: the `.log` is written in full either way, and so is the
  command line's own console output, which this does not touch.
* **?** — the parameter reference, which is `cli.USAGE` verbatim plus a
  description of what the tool produces. There is one description of the
  parameters, so the window and `--help` cannot disagree.
* **Start!** / **Exit**.

Every field is validated by calling `cli.parse_args` and `cli.preflight_checks`
on the command line the window is about to run — not by a second copy of those
rules. `Start!` is disabled while anything is invalid and the parser's own
message is shown, so `t < cc/√2` is a greyed-out button rather than a failed run.

**While a run is in flight** the fields go inactive but stay readable, `Start!`
becomes `Stop!`, and two bars appear:

* the **top** bar names the running stage and fills according to a fixed share
  per stage, taken from the 2026-08-18 controlled pair in
  docs/profiling-reports.md and declared in
  [`gui/weights.py`](../src/latticegen2/gui/weights.py);
* the **second** bar shows work within the stage — `boundary trim: 9,776 /
  19,552` — with the text drawn over it.

Beneath them, one line of resource data: peak memory and elapsed time. Peak
memory is genuinely all that is measured continuously (§3's summary list), and
the line says only that rather than padding itself out.

Beneath *that*, with **verbose** ticked, the run's own output in a scrolling
pane — every line of it, which is what the `.log` gets.

**Unticked, the pane holds only what the child wrote outside the event stream:**
kernel chatter on stdout, and anything at all on stderr, which is where a
failure's one reason line lands (§7). It is hidden entirely when that is empty,
so an ordinary run does not carry a blank box around for the hour it takes —
and the box appearing is itself the signal that something spoke up.

**Nothing else is shown unticked, deliberately.** The window is not a terminal,
and everything a clean run logs is either already on screen in a widget — the
parameters are in the input boxes, the stage in the top bar, duration and peak
memory in the resource line, the outcome in the result banner — or a statistic
that is neither a warning nor an error. Repeating it as text would be the
window competing with itself.

**This is a window-side rule, and it is deliberately not the `console` flag the
events carry.** That flag says whether a *command-line* run would have printed
the line, and §3 requires the whole end-of-run summary to print there
regardless of `-v` — which is right for a terminal and wrong here, the window
having drawn most of that summary already. The command line is unaffected by
any of this.

**The tick box is the one control that stays live while a run is in flight**,
and that is the point of it. The child is never given `-v`: every line crosses
the event stream either way (docs/algorithm.md §10), so verbosity here is a
filter over output already received rather than a flag that had to be decided
before the run started. Tick it forty minutes in — or after a failure, to read
what led to it — and the pane fills in what was hidden.

**Where a stage has no countable work — `export`'s single writer call, or
`simplify` while its one dominant solid is unified — the second bar sweeps and
says so, and the top bar holds.** Neither invents a fraction. The alternative,
filling on a timer, produces a number that measures nothing and sits at 100 % of
a stage that has not finished.

**The weights are one part's shape on one machine**, and a small part is
boundary- and export-dominated, so the bar will visibly jump there. That is
accepted: a fixed weighting is monotone and never claims a stage finished before
it did. Read the bar as how far through the work, never as an estimate of time
remaining.

**Stop is exactly Ctrl+C.** It causes the same graceful shutdown §3 already
specifies — workers stopped in order, one `CANCELLED` line, exit 130, and the
temp folder **left in place**, whose location the window then shows. It is not
instantaneous: the interrupt is delivered between bytecodes, so inside a long
kernel call it lands when that call returns, and the button says `Cancelling…`
rather than pretending otherwise. If the run has not stopped within a short
grace period the window force-stops the whole process tree, which on Windows
matters — killing the master alone would orphan every worker.

**When a run ends** the window reports success or failure, offers **Open export
folder**, re-enables the inputs and returns `Stop!` to `Start!`, ready for
another run with the same or different parameters. A failure shows the same
single reason line the command line prints (§7).

**Zero arguments open a window only where a window can exist.** On a machine
with no display a bare invocation still exits 2 with usage, exactly as it always
has, so scripts and CI are unaffected. `--gui` given explicitly always tries and
reports one line if it cannot.

**Two things deliberately absent from §3's table**, because they are transport
rather than parameters:

* `--progress-stream` makes a run report itself as machine-readable events on
  stdout. The window sets it on the child it launches; it changes nothing about
  the geometry, and `tools/e2e.py`'s `progress-stream` scenario proves that by
  running the same case with and without it and comparing the bytes.
* A **cancel sentinel file**, `<output-stem>.cancel`, is how Stop reaches the
  run. Neither of the obvious channels works: a windowed process has no console
  and so cannot send Ctrl+Break, and giving the child a pipe for stdin was
  measured stalling the worker pool outright (docs/algorithm.md §10).

**Toolkit.** `tkinter`, from the standard library — no new packaged dependency.
Portable bundles consequently redistribute Tcl/Tk, which
`tools/build_release.py` proves is present by importing it at build time and
`tools/smoke_bundle.py` re-checks in the extracted archive. See
[licenses/LICENSES.md](../licenses/LICENSES.md). The **wheels** flavour and a
source checkout use the operator's own interpreter, where on Debian and Ubuntu
`tkinter` is the separate `python3-tk` package; without it the command line is
unaffected and only the front-end is unavailable.

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
  
- **No body is ever dropped to make an export succeed, and the run fails
  instead.** A body the generator cannot write faithfully is a hard failure
  (exit 4) naming the face and its position, with the temporary folder kept —
  not a body quietly removed from the output. §1 asks for a lattice filling the
  user's volume; silently shipping less of one is a wrong answer rather than a
  degraded one, and the size of the piece does not change that.

  **What makes a body unwritable is a property of STEP, not of this generator.**
  AP214 carries exactly one modelling tolerance for a whole file — the
  `UNCERTAINTY_MEASURE_WITH_UNIT` of its representation context — where an OCCT
  B-rep carries one per vertex, per edge and per face. Export collapses them and
  import re-derives them all from that single number, so a body whose validity
  rests on a locally fat tolerance is valid in the generator and is not
  guaranteed valid in the file. Measured: a 6.573e-02 mm vertex tolerance comes
  back from a round trip at 1e-07, and `dense-lattice`'s dominant solid loses
  its 5.151e-04 mm edge tolerances to a declared 2.E-07.

  The run therefore measures the quantity that decides — how far each pcurve
  strays from its own 3D curve, against the size of the face carrying it — on
  **every** output solid, large and small alike, and refuses past a bar of 1e-2.
  It is also measured at the source, per trimmed junction, so a failure can name
  the junctions responsible rather than only a coordinate. See
  docs/algorithm.md §7.3 and §9.

- **STEP schema/AP:** AP214

- **Geometry representation in the file:** exact B-rep solid
  
- **Units:** mm
 
- **Metadata to embed:** Part name as concetenated <input_file_name>+lattice+cc<cc>+t<t> and generation parameters as STEP header.
  The part name carries the same four components as the default output file
  name (§3), in the same order, so a body opened in Solidworks or Catia is
  recognisable as the file it came from. The two differ only in punctuation:
  `+` between every component here, against `-` there with `cc` and `t` run
  together (`ball-lattice-cc20t4.step` carries `ball+lattice+cc20+t4`).

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
| dense-lattice | -i test/test-cylinder.STEP -cc 10 -t 1.5 --cores 6 | valid STEP, no self-intersections, matching golden sample test/test-cylinder-cc10t1.5-golden-sample.step, generation < 10 minutes. **Measured: 47.5 s, symmetric-difference volume 0 mm³.** |
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
- **The shipped file's pcurves still agree with their own 3D curves.** Read the
  output back with OCCT's own reader and measure, per edge/face pair, the exact
  distance between the two representations against the area of the face carrying
  it. This is asked of the *artefact* rather than of the process that wrote it,
  which is the only version the downstream tools see.
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
- The graphical front-end (§3.1) reports a failure with the same single reason
  line, in the window. A front-end that cannot *start* — no display, or an
  interpreter without `tkinter` — reports itself on stderr, or in a dialog box
  when there is no stderr to print to, which under `pythonw.exe` there is not.
  Without that, a broken interpreter would be a double-click that does nothing
  at all.

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

### Re-fit the pcurve at the source, so the body is kept rather than refused

**Found 2026-08-23**, alongside §11's export-truth gate, and deliberately not
attempted there. The gate is the right behaviour today and the wrong end state:
it tells the user a body cannot be written and stops, where the fault is
repairable and the body could be kept.

**What is wrong.** A boundary trim against a fat curved surface can leave an
edge whose pcurve does not match its own 3D curve — measured up to
**2.118e-02 mm** apart on a face of 0.05 mm², on `SpiralTest.step` at
`cc=5, t=1`. In this process that is legal, because the edge records a tolerance
large enough to cover it. In the exported file it is not, because STEP AP214
declares one tolerance for the whole file (§11), so the reader re-derives a
tight one and the two representations no longer agree. The body then tessellates
into a shell with holes: 11 edges of its 147 triangles used by one triangle
rather than two.

**Where.** The trim is `boundary.trim_junction` (docs/algorithm.md §7). That is
also the only place the junction is still identifiable — §7.3's
`tolerance_feature_ratio` already ranks the offenders there, and on this part it
puts two of the six junctions forming the unwritable island **2nd and 4th of
2,404 boundary pieces**, at the end of `boundary`, 5 m 55 s into an eight-minute
run. So the information needed to aim a repair is already measured and already
in hand before anything downstream exists.

**What to try.** Re-fit the offending pcurve so it agrees with its 3D curve to
within a tolerance the file *can* carry — `ShapeFix_Edge::FixAddPCurve` or
`ShapeConstruct_ProjectCurveOnSurface`, asked for a tight tolerance rather than
the one the boolean settled for. Failing that, re-derive the trim for that
junction the way docs/algorithm.md §7.1's disagreeing-cap repair does: give the
kernel operands it handles better and redo the intersection locally.

**Two things already ruled out, so they are not retried.** `ShapeFix_Shape` and
`BRepLib::SameParameter` were both measured against this defect and neither
helps — they were being asked to repair damage that does not exist until export
(§11). And no export-side setting reaches it: `write.precision.mode` in every
mode, `write.surfacecurve.mode = Off`, and `SameParameter` before writing were
each measured and each leave it (§11).

**How to verify.** `test/spiral-island-unwritable.brep` is the body, committed
for exactly this. Today it fails `occ.exported_mesh_defects` with **11**
non-manifold edges in 147 triangles
(`test_export_truth.py::test_the_island_does_not_survive_being_written`, which
pins `(147, 11)`). A working repair takes that to `(n, 0)` while preserving the
solid's volume, and the whole-part check is `SpiralTest` at `cc=5, t=1`
completing and writing its STEP. Note that two *other* defects stop that part
first on this branch — see the item below and §11.

### Re-run the rehearsal with the unbounded export-truth check

**Deferred by decision 2026-08-23.** §11's export-truth gate was measured on
`TD_HX_rehearsal_test` at `cc=5, t=1` *while it still skipped solids above a
face count*; the bound was then removed so the dominant body is checked too, and
the rehearsal has not been re-run since.

**What is unknown, precisely.** Two things, and they are different questions.
*Does it pass* — whether the 583,806-face dominant body's tessellation survives
its own round trip. Its thirteen small siblings do; the 80 mm ball's dominant
body ships 73 `InvalidCurveOnSurface` faults and `dense-lattice`'s picks up 4,
both harmless and both tessellating cleanly, so this is genuinely open rather
than rhetorical. *And what it costs* — docs/algorithm.md §9 removed a
whole-output re-import for costing **22 minutes**, and this adds a write and a
tessellation on top of that re-read, on one 2 GB solid.

**How to verify.** `python src/main.py -i test/TD_HX_rehearsal_test.step -cc 5
-t 1 --cores 6 -v`, and read `export_truth_s` in the run summary. Watch peak
memory as well as the clock: mesh points are interned to integers and edges
counted by integer key specifically so a solid this size can be attempted, and
that is an argument rather than a measurement until this run exists. A refusal
here is a real finding about the part, not a miscalibration — the instrument has
a clean record on all sixteen bodies measured so far.

---

## 11. Closed — kept for the reasoning, not as work

### A body can be valid here and not describable in the file — the export-truth gate

**Found and closed 2026-08-23.** `BRepCheck_Analyzer` passes every output solid,
`assemble` proves every one watertight, containment holds — and the file the
user receives can still be inconsistent, because **STEP AP214 has no
representation for per-subshape tolerance**. A file carries exactly one
`UNCERTAINTY_MEASURE_WITH_UNIT`, in its `geometric_representation_context`,
against one tolerance per vertex, per edge and per face in an OCCT B-rep. Export
collapses N into 1; import re-derives all N from that one number.

**This is not incidental to the pipeline, it is aimed at its own repairs.**
docs/algorithm.md §8's two-rung repair fixes a falsely self-intersecting wire by
*widening a recorded vertex tolerance*, and the property that makes it safe —
"it moves no geometry" — is exactly the property that makes it unexportable. The
same is true of every fat tolerance a boolean records when it trims a strut
almost tangentially to a curved input surface.

**Measured, and pinned by `test/test_export_truth.py`:**

* A vertex tolerance of 6.573e-02 mm — the figure OCCT itself records on
  `SpiralTest`'s fat vertex — written to STEP and read back comes home at
  **1e-07**.
* `dense-lattice`'s dominant body, the same solid before and after its own round
  trip: worst pcurve↔3D deviation 5.1514e-04 mm with a max edge tolerance of
  **5.151e-04** covering it exactly and **0** of 62,792 pairs over tolerance;
  afterwards the tolerances are clamped to **1.525e-04** and **4** pairs are
  over. The file had declared `2.E-07`, because OCCT's `write.precision.mode`
  defaults to *Average* and a lattice averages ~99 % exactly-built interior
  edges at `Precision::Confusion` against the 1 % of boundary trims that carry
  real tolerance.

**Three export-side levers were tried and none fixes it**, recorded so they are
not retried: `write.precision.mode` at Greatest, Least or an explicit session
value leaves the ball's over-tolerance count at 60–95 in every mode;
`write.surfacecurve.mode = Off` makes the worst deviation *worse*
(1.5e-06 → 6.3e-06 mm); `BRepLib::SameParameter` before writing changes nothing.
Coordinate precision is ruled out with a number — 14 significant digits, ~1e-11
mm at 2,000 mm coordinates, six orders below the tightest tolerance in play.
What OCCT does do is mark every `surface_curve`'s `master_representation` as
`.PCURVE_S1.`, so where the two representations disagree the file tells the
reader to believe the pcurve.

#### What was built: measure at the source, gate on the output, never drop

Three parts, and the middle one took two wrong turns before it measured right.

**1. At the source (docs/algorithm.md §7.3).** Every trimmed piece is measured
in the worker, on the face list the trim already produced: worst
`edge tolerance / sqrt(face area)` over its faces. One area and one centroid on
the single worst face — and the last moment at which the junction still has a
name. On `SpiralTest` at `cc=5, t=1`, **two of the six junctions forming the
4.17 mm³ island rank 2nd and 4th of 2,404 boundary pieces**, at 2.907e-01 and
2.093e-01. That is available at the end of `boundary`, 5 m 55 s into an
eight-minute run, before the junction graph that turns them into a body exists.

**It reports and does not refuse, and that is measured rather than cautious.**
79 of 2,404 pieces clear the 1e-2 warning bar, most of them welded into the
27,864 mm³ dominant body where a locally loose description is absorbed and the
exported solid is sound. Failing on it would refuse a part whose output is fine.

Combined with connectivity it is much sharper, but **not in the way the first
implementation assumed**: both surviving components contain a flagged junction —
the dominant body holds the single worst one in the whole part — so the maximum
says nothing. The fraction does: 82 of 2,348 boundary junctions (3.5 %) for the
lattice against 2 of 6 (33.3 %) for the island. Reported per body at `connect`,
about 40 s in.

**2. On the output (docs/algorithm.md §9), and the instrument took four tries.**
The question is whether a body survives being written, so the check asks exactly
that: write one solid to STEP, read it back, tessellate it, count edges not used
by exactly two triangles. Cheap on the bodies it exists for — the rehearsal's
thirteen small solids cost well under a minute between them.

Three cheaper quantities were tried first and are recorded because what they
cost to learn is the useful part. On `SpiralTest` alone, two look decisive: a
fault count is blind (the broken island has **0**, the accepted ball's own
output has 73), the worst face is backwards (the sound lattice scores 3.07e-02
against the island's 2.45e-02), and the share of surface described more loosely
than its own feature size separates them by **706×**.

**Then the rehearsal was used as ground truth and the third one
false-positives.** All fourteen of its solids were round-tripped and
tessellated:

| body | faces | loose area | faults after RT | **bad mesh edges** |
|---|---|---|---|---|
| `SpiralTest` island | 36 | 3.97e-01 | 29 | **11** |
| rehearsal unify 3 | 12 | 1.76e-01 | 0 | 0 |
| rehearsal unify 5 | 12 | 1.76e-01 | 0 | 0 |
| rehearsal unify 10 | 29 | 0 | 8 | 0 |
| rehearsal unify 13 | 7 | 0 | 1 | 0 |
| rehearsal, nine others | — | 0 | 0 | 0 |

Exactly one of the sixteen bodies is broken. The loose-area fraction refuses
unify 3 and 5; a fault count refuses unify 10 and 13; only the tessellation is
right about all of them. **Both rejected proxies would have refused — or, in the
shape of rule this replaces, deleted — geometry from a part that has been
inspected and accepted.** The pcurve figures are still measured and logged
because they say why a body is fragile; they decide nothing, and a unit test
pins that.

**Every solid is measured, with no size bound — the user's decision, and the
expensive one.** An earlier revision skipped solids above a face count and
reported them *unmeasured*, which was honest and wrong in the way that matters:
it put the dominant body of every production part outside the only detector with
a clean record. docs/algorithm.md §9 removed a whole-output re-import for
costing 22 minutes, and this is that cost returning knowingly, for a correctness
check rather than for a count of solids. The rehearsal-scale figure is **not yet
measured** and is expected to be tens of minutes; `export_truth_s` in the run
summary is what to read.

**Why the gate is needed at all, demonstrated on the real body.** As this
pipeline builds it, the island is `BRepCheck_Analyzer`-**invalid** — rung 2
declines its fat vertex and the validity gate refuses the run. Let rung 2 act
(widen the vertex; nothing moves) and **the body becomes valid**. Written, its
147 triangles still carry 11 broken edges. The repair that makes a body pass
every gate the pipeline had is precisely what hides the remaining defect.

**3. Past the bar the run fails (exit 4), and nothing is ever discarded.**
Deleting material so an export can succeed is "produce a wrong result", which
docs/algorithm.md §11 forbids; §1 asks for a lattice filling the user's volume,
and silently shipping less of one is a wrong answer rather than a degraded one,
whatever the size of the piece. The temp folder is kept, the message names the
face and its position, and `connect` has already named the junctions.

#### What this leaves open

**`SpiralTest` at `cc=5, t=1` does not complete on this branch**, and it stops
*before* the new gate: at `simplify`, on the same 4.170538 mm³ island, at the
unification volume guard — a *separate* defect this branch deliberately does not
touch. Relaxing that pre-filter in the measurement harness only, the run then
reaches `validate` and the island is refused there by the existing
`BRepCheck_Analyzer` gate, because rung 2's fixed cap declines its fat vertex.
So on this branch the part is refused twice over before the export-truth check
is what stops it — which is why the demonstration above is made on the committed
island fixture, where rung 2 can be allowed to act.

**The rehearsal's thirteen small solids are refused by nothing**, measured on
the solids the run actually produced and replayed from its kept temporary
folder. Its **583,806-face dominant body has not been put through the check
yet** — the size bound that used to exclude it was removed after that run, on
the user's decision. Re-running it is deferred rather than skipped, and is §10's
second item.

**What the dominant bodies carry is worth recording either way.** The 80 mm
ball's own dominant body ships **73** `InvalidCurveOnSurface` faults and
`dense-lattice`'s picks up 4 from its own round trip — harmless at their
magnitudes, both tessellating cleanly, and nothing in the pipeline was watching
before this.

**Refusing is the correct behaviour and not the end state.** It is strictly
better than removing the body — the user learns the part has a problem, learns
where it is, and loses nothing — but the fault is repairable, and repairing it
at the source would keep the body instead of stopping the run. That is real
work rather than closed reasoning, so it is written up as an item in §10 with
what to try, what is already ruled out, and how to verify it.


### `stitch`'s round-2 repair — chapter closed: the fix disproved, the check repaired, the scan parallelised

**Opened 2026-08-17, closed 2026-08-18.** §10's last live item, in two parts:
the ~545–651 s full unsplit sew that boundary-sew round 2 falls back to on this
part (85 % of `stitch`), and the 44.1 s serial validity scan beside it. Both
are resolved, neither the way the item expected, and the second half of each is
the part worth keeping.

#### The proposed fix is structurally the fallback — G21

§10 recorded one lever "worth more than ~1 %" and marked it unproven: carry a
seam face's **straddling** neighbours — edges used twice inside a tile but only
once inside the seam-only subset — into the sewn subset, so
`BRepBuilderAPI_Sewing` cannot rebuild one while the carried face keeps the
original. [`g21_straddling_seam_split.py`](../tools/prototypes/g21_straddling_seam_split.py)
measured it on the real part, at full extent.

**A subset with no straddling edge is a union of connected components of the
tile's face-adjacency graph**, and a tile's round-1 result is one connected
shell but for a few strays. So the closure is the tile:

| Tile | Faces | \|S0\| (today) | \|S1\| | \|S2\| | Fixpoint | Hops |
|---|---|---|---|---|---|---|
| 0 | 23,546 | 7,555 | 12,491 | 19,954 | **23,523** | 8 |
| 1 | 22,872 | 6,514 | 10,966 | 17,769 | **22,837** | 8 |
| 3 | 16,831 | 4,116 | 6,983 | 11,635 | **16,723** | 18 |

and closing to it costs what the fallback costs — **523.20 s against the full
unsplit sew's 520.89 s**, on the same 302,576 faces where the seam-only split
takes 9.34 s. One hop is cheap (43.42 s) and is not a fix: it moves the
frontier, leaving a fresh set of straddling edges of exactly the kind it exists
to remove. **Do not build this**; the table is the disproof.

That closes the stage. `stitch`'s repair is the price of a correct shell on
heavily trimmed geometry, not an unfixed inefficiency, and the seam-only split
stays because it is nearly free where it works and caught where it does not.

**Nor is it worth predicting the failure to skip the attempt.** On a part where
the split always fails, everything spent before the repair is waste — and the
run log now prices it: `split 2.8s, round2 14.5s`, **17.3 s of a 53-minute
run**, 0.5 %. A heuristic that guessed "this component is too heavily trimmed,
sew it unsplit" would save that at the risk of skipping the split on a part
where it works, which is the wrong side of docs/algorithm.md §11's rule. This is
the same arithmetic that closed the speculative-sew proposal in §10, and it
closes this one too.

#### The check could not have told success from failure — and that is now fixed

`_sew_round_two` now records `(component, want, got_split, got_unsplit)`
whenever it repairs, which costs nothing because the unsplit sew has to run
before either number exists. The rehearsal:

    component 0: expected 73984 free edge(s), seam-only split gave 192692,
                 full unsplit sew gives 73994

The first two numbers confirm G9's mechanism is live at production scale — the
split is wrong by a factor of 2.6, not by a rounding error. **The third is the
finding.** A correct sew missed the expectation too, by exactly 10, because
`free_edges` counted **degenerate** edges — parametric artefacts with no
extent, whose one owning face uses them once by construction. `weld.shell_defects`
has skipped them since the day it was written and records finding exactly 10 on
this part (3.0e-9 to 8.3e-8 mm); the two counts should always have agreed, and
now do.

Until this, the check fired *whatever* round 2 produced. It happened to fire
for the right reason here, but on a part whose split was sound it would have
forced a full re-sew anyway, and `stitch_repaired_components` could not have
meant what it says. This is docs/algorithm.md §11's rule about the quantity a
gate compares, in the one place where both quantities were already in the file.

#### The scan beside it: 44.1 s → 22.6 s, and two lessons

The other half of §10's item. `occ.fix_vertex_tolerances` asked one
`BRepCheck_Analyzer` per boundary face to find the 19 it repairs — 44.6 s over
301,505 faces when §10 recorded it, 44.1 s in the controlled pair below. §10 had rejected the two obvious speed-ups, and both rejections
stand: dispatch has no free-I/O route, and one analyzer over the assembled
*solid* is the **in-context** overload, a different predicate.

A `TopoDS_Compound` of loose faces is neither, and G22 measured **1.66x** at
3.57 core-equivalents — **44.1 s → 22.6 s** in the rehearsal pair — with zero
disagreements over a 17,308-face corpus
carrying all four committed faults — chunked at 20,000 faces, since the
analyzer holds ~14 kB per face and one call over the rehearsal's would hold
~4.2 GB.

**The change was nonetheless wrong on first landing, and the rehearsal said so
within the hour**: 19 repaired and **15 "still invalid"**, where the serial scan
had always reported none — on a run whose validity gate then passed all 14
solids, which is what marked them as phantoms rather than news.

**The first diagnosis of those 15 was wrong too, and the pair disproved it.**
The obvious reading is a predicate difference: G22's faults are all *loose*
faces, a sewn layer's have neighbours, and a compound analyzer holds one result
per subshape shared between the faces using it — so a face could be rejected in
batch for the fault beside it. A confirmation stage was added on that reading
and the count stayed at 15, which settles it: all 34 candidates really are
invalid standalone.

**The cause is *when* the predicate is evaluated, not which one.** Both repair
rungs widen tolerances on vertices and edges that neighbouring faces share, and
widening is monotonically permissive — so repairing one face can make the next
one valid before the loop reaches it. The serial loop asked at the moment it
arrived at each face and never saw them; a scan-then-repair pass asks before any
repair has happened. Re-checking each candidate as the repair reaches it
restores the original set exactly, and cannot fail the other way, since no
repair here invalidates a face. The confirmation stage stays as cheap insurance
on the case G22's corpus cannot reach, now labelled as insurance rather than as
a measured necessity.

**Two lessons, and the second is the one that keeps recurring.** A control whose
faults are isolated cannot test a predicate whose difference is about
neighbours — G10's lesson in a new place. And a mechanism that *fits* the
symptom is not the mechanism: this is the fifth time in this project's history
(G9, G10, G11, G12) that the first convincing explanation was the wrong one, and
what settled it here was the cheapest possible test — ship the fix the
explanation implies, and see whether the number moves.

It sits beside the other gate lesson this session produced. G21's part A sewed
the same tiles both ways and found them *identical*, because the gate rebuilt
the pipeline's front half and stopped one call short of `finalize_pieces`, which
is what opens the interface holes. A prototype that rebuilds part of the
pipeline has to rebuild all of it, or it measures a shape the pipeline never
sees.


### `--ram` removed — an accepted-but-unenforced budget, taken out rather than fixed

**Found 2026-08-17**, moved here and closed the same day. Originally logged in
§10 as "`--ram` is accepted and validated but never enforced": the flag was
parsed, range-checked against the machine's physical RAM, resolved to
`Args.ram_budget_gb`, and then **only printed**. Nothing in `pipeline.py`,
`parallel.py`, `weld.py`, `boundary.py` or `runlog.py` ever read it, which
directly contradicted `cli.py`'s own module docstring calling both budgets
"ceilings on what a run may use, not hints it may exceed" — true of `--cores`,
false of `--ram`. The write-up also found a second, independent defect: the
peak `--ram` would have been checked against already under-counted by ~1.2 GB
(`RunLog.max_rss` folds in worker RSS only where a stage explicitly calls
`note_worker_rss`, which not every stage does), and it recorded why real
enforcement is harder than the flag suggests — the run's peak sits in the
master holding the finished result while `export` serialises it, a place no
worker-count or tile-size lever can reach, so genuine enforcement would mean
streaming the STEP write or exporting per solid, a substantially larger change
than the flag ever implied.

**The decision.** Given a choice between (a) keeping the budget advisory and
making the docstring and summary honest about that, or (b) building real
enforcement, the user chose neither: **remove `--ram` outright.** A budget nothing
enforces and that costs a materially larger change to enforce for real is not
worth carrying as CLI surface, log fields, and documentation weight for what it
actually does today — print a number back at the user.

**What was done.** `--ram` is gone from the CLI: `cli.py` no longer parses the
flag, validates its range, or carries `ram`/`ram_budget_gb` on `Args`; the
`--cores` budget is unaffected and still the only one. `sysinfo.py`'s
`total_ram_gb`/`free_ram_gb` are deleted along with it, since the `--ram`
ceiling and default were their only callers. That in turn left `psutil` with no
remaining use anywhere in `src/` — core-count detection has always come from
the standard library — so `psutil` was dropped as a **packaged runtime
dependency**: out of `pyproject.toml`'s `dependencies`, `requirements-bundle.txt`,
`licenses/LICENSES.md`'s table (and `licenses/psutil-LICENSE.txt` deleted with
it, since nothing shipped needs the text), README.md's dependency table and
install command, and the `latticegen2.bat`/`latticegen2.sh` launchers' dependency
probe. `psutil` remains exactly what it already partly was: a development-only
dependency of `tools/profile_run.py`, which was never bundled either way — the
same footing as `pytest`.

Every other reference to `--ram` or the budget it fed was removed from CLI help
text, `docs/algorithm.md`, `docs/testing.md`, `docs/specification.md` §2/§3/§6.1,
`README.md` and `CLAUDE.md`, and from the committed test suite
(`test_cli.py`'s three budget tests and `test_sysinfo.py`'s memory-detection
test) and `tools/e2e.py`'s `dense-lattice` invocation. Historical measurement
records that quote real past commands run with `--ram` — `docs/profiling-reports.md`,
`tools/prototypes/RESULTS.md`, and the already-closed narrative entries later in
this chapter — were left untouched: they are verbatim accounts of what was
actually run at the time, not current usage instructions, and rewriting them
would misrepresent history rather than clarify it.

**Verification.** `python -m pytest test -q` and `python tools/e2e.py` (unit
tests plus the four committed end-to-end scenarios, `smoke-fast`,
`smoke-verified`, `dense-lattice`, `invalid-input` — golden-sample comparison
included, since e2e.py always runs it, but no separate rehearsal-scale run) both
green after the change. `--ram` now behaves like any other removed flag
(`-bg`, `--background`, `--workers`): rejected as `Unknown argument`, pinned by
`test_cli.py::test_removed_flags_are_rejected`.

### Pipeline parallelism between `classify` and `assemble` — chapter closed, two stages won, two proposals disproved

**Opened 2026-08-17** by the question: the pipeline runs at ~1.9 of 6 cores, so
is there work between classification and assembly that could run concurrently
rather than in sequence? Specifically — since interior nodes are handled
separately from boundary nodes before the interior is built onto the interface
rings, can the interior build be dispatched to the worker pool alongside the
boundary work?

**The originating hypothesis was structurally right and worth ~42 s.** The
dependency really is narrower than it looks: of `build_interior_shell`'s seven
arguments, six are available the moment `connect` finishes, and only `adopted`
comes from `stitch`. `_rings_needed` emits an entry solely where an interior
node's neighbour across a cap is *not* interior, so the adopted set is the
shell's outer **skin** — 18 k rings against 97 k interfaces. Every node whose
six neighbours are all interior depends on `stitch` for nothing.

**But it cannot use the worker pool, for the same reason Phase 3 died.** The
shell's watertightness *is* pointer identity, and `_ShellBuilder.adopt` must
take the boundary sew's actual `TopoDS_Edge` objects, which live in the master's
heap. And `instance` measures 46 s of a ~3,100 s run, so the ceiling is ~1 %
even done perfectly. It was not built.

#### What the search found instead

Asking "which stages are single-core and *why*" turned out to be the productive
form of the question, and the answer separates cleanly into two groups.

**Two stages shipped.**

| | `dense-lattice` (controlled pair) | rehearsal (clean window) |
|---|---|---|
| `classify` | 10.70 s → **4.08 s** | 126 s → **47.0 s**, 0.98 → **4.02 cores** |
| `validate` | 5.49 s → **2.11 s** | 225 s → **114.0 s**, I/O 464/464 MB → **0/0 MB** |

`dense-lattice` went 50.30 s → 40.50 s (−19.5 %) with the output
byte-identical outside the header timestamp; the rehearsal wrote its 14 valid
solids, 584,028 faces, 330,354.002 mm³ unchanged.

**The distinction that made both possible, and that generalises past them:**
every stage the earlier chapters failed to parallelise failed for one of two
reasons — OCP holds the GIL (G7, G17), or tiles must reassemble by shared
topology and a file boundary destroys it (G15). Neither applies to a stage that
moves no geometry:

* **`classify` is pure NumPy.** Nothing OCCT crosses its process boundary — the
  mesh and node indices are plain arrays over an `.npz` — so both findings are
  simply irrelevant, and the sweep divides exactly because every node is decided
  independently of every other. Slices are *strided* rather than contiguous:
  `candidate_nodes` ravels a meshgrid, so contiguous slices are slabs, and the
  only expensive step runs solely for near-surface nodes, which lie on a shell.
* **`validate` returns a scalar, not geometry.** That is precisely what
  `simplify` lacks, and it is why sub-body parallelism is closed there and open
  here. It also turned out to need no decomposition at all: G18 found
  `BRepCheck_Analyzer` has its own `theIsParallel` flag — OCCT's *native*
  threads, which the GIL result does not bind. 1.60×, verdict identical on all
  four real invalid faces committed from this part.

**Two proposals were disproved by measurement, and both had been ranked above
the two that shipped.**

* **Speculative round-2 sew — 15.1 s, not minutes.** On this part `stitch`
  computes a seam-only round 2, fails the `expected_rings` check, and re-sews
  the dominant component from scratch on the master while workers idle. The
  proposal was to run both jobs concurrently and keep whichever passes. Adding
  per-phase timers to `stitch` priced it: `round1 49.1s, split 3.0s, round2
  15.1s, repair 651.2s, retolerance 44.6s, rings 6.7s`. **The discarded attempt
  costs 15.1 s** — 0.5 % of the run — because the seam subset is small. The
  repair is 85 % of the stage and no scheduling change reaches it (§10).
* **A chunked per-face split of `validate` — measured worth ~5×, rejected on
  correctness.** G18 measured per-face checks at 94.4 % of the serial cost
  against a 4.8 % structural floor, so it would beat the flag substantially.
  Not built: `BRepCheck_Analyzer(solid)` checks subshapes **in context**, and a
  standalone per-face check is a different predicate — the very difference G12
  used to diagnose rung 2. Replacing the exit-4 gate with a hand-assembled
  conjunction would make its failure mode "produce a wrong result", which
  algorithm.md §11 forbids. The 94.4 % is recorded so the option stays visible
  if a control for in-context faults is ever built.

#### Three things worth keeping beyond this chapter

**A stage timer that lumps phases together prices the stage, not the proposal.**
`stitch` had sat at "0.96 mean cores against a 4.38 peak" across several
profiles — visible but unattributable. Six timers, costing nothing, retired a
proposal on their first run. Had they existed earlier, the speculative sew would
never have been ranked first.

**`--cores` and OCCT's own thread pool are two different budgets, and using both
at once double-counts.** Once `BRepCheck_Analyzer` runs threaded, `W` worker
processes each launching `W` threads is `W²` on `W` cores — so the flag and
path 4's per-solid dispatch cannot coexist. Choosing the flag means running
`validate` on the master, which also deleted a 464 MB round trip that existed
only to reach the workers. `parallel.set_thread_budget` caps the OCCT side.
Anything added later that asks OCCT for threads inherits this constraint.

**Cross-session stage comparison remains the trap this project keeps falling
into.** `boundary` in the verification run measured 889 s at 4.28 cores against
the previous session's 721 s at 5.20 — and was briefly written up here as a
possible regression, before the 2026-08-15 profile was found recording 14 m 51 s
at **4.22 cores** on untouched code, labelled "no — variance" in §10's own
table. The re-run reproduces that slow-band point to within 1 %. A stage
measured in a clean window and a whole-run total are not the same kind of claim,
and only the controlled pair supports the second.

#### What this chapter leaves open

Nothing large. `stitch`'s 651 s repair is §10's existing item and is
algorithmic, not a scheduling change. The retolerance scan (44.6 s, ~1.2 %) is
rescoped there with two findings that make the obvious implementation worse than
it looked. The interior-overlap hypothesis above is ~1 % and needs a background
submission path in `parallel.py` plus a two-phase shell build.

One loose thread, recorded rather than resolved: the verification run wrote
**2,517,853 edges against the 2,517,881** of the Phase 2 pair, a difference of
28 (0.001 %). Nothing in this chapter can plausibly cause it — the
classification is identical, `validate` is read-only, and `set_thread_budget(6)`
is a no-op where OCCT's default is already 6 — and every other figure matches
exactly with all 14 solids valid. It sits inside what algorithm.md §9 calls
unification's own representation choice. A future rehearsal should confirm it is
variance rather than drift.

### `simplify` beyond body-for-body — chapter closed, two levers disproved

**Proposed by Claude 2026-08-15, approved for implementation 2026-08-16,
closed 2026-08-17.** Of four things tried on this stage, **one paid**: Phase 2,
building the interior pre-merged, −6.5 % on the rehearsal. Phase 1's projection
was wrong, Phase 3 is impossible as designed, and the input-side alternative
built in its place is exact and worth nothing. All four are recorded here rather
than quietly dropped, because what they cost to learn is the useful part — and
because two of them close off directions that would otherwise look obvious to
the next person.

*Read this as one chapter: the later sections' bars are set by the earlier
sections' measurements, so splitting them apart would leave each half
unreadable.*

The measurement the proposals rest on is gate G13
([`tools/prototypes/RESULTS.md`](../tools/prototypes/RESULTS.md),
`tools/prototypes/g13_unify_scaling.py`).

**Why this exists.** Path 1 of the scale-rehearsal chapter (§11) —
parallelising `simplify` across the shared pool, body for body — is
implemented, correct, and measured at **no wall-clock win** (17 m 17 s -> 18 m 39 s, still 0.99 cores), because this part's 14 solids
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

#### Phase 3 — sub-body tiling of `simplify`: **BLOCKED, not built**

The plan was to partition a solid's faces into tiles, unify each across the
shared `WorkerPool`, and reassemble with `BRep_Builder.Add` — resting on G13's
and G14's finding that unification leaves the edges it did not merge as the same
objects, so a tile seam stays one `TShape` and nothing has to be sewn.

**Its gate condition was met and its risk list was wrong.** After Phase 2
`simplify` is the largest single stage — 18 m 07.9 s of a 51 m 43.3 s run
(35 %), pinned at 0.98 cores while five of six sit idle — so it was worth
attacking. R2 (curved faces) had been cleared by G14; R3, R5 and R6 were open
but manageable. **None of the four risks asked whether the mechanism survives
the process boundary the plan's own first sentence puts it across** — which is
the one that decides it, and which two gates had silently assumed away.

**G15 answers it: tile identity is pointer identity, and Phase 3 puts a process
boundary through the middle of it.** G13 and G14 both measured reassembly *in
one process*. Phase 3's entire purpose is to spread tiles across workers, and
every stage that does so moves geometry as `.brep` files (docs/algorithm.md §7,
§8). A `.brep` preserves sharing *within* a file and cannot preserve it
*between* two, because each file writes its own copy of every edge it
references. On G14's own trimmed test solid
(`tools/prototypes/g15_tiled_unify_ipc.py`):

| reassembly route | 8 tiles | 27 tiles | |
|---|---|---|---|
| in-process (G14's route) | **0** free edges | **0** | control: reproduces G14 |
| one file, all tiles as one compound | **0** | **0** | control: serialization itself is fine |
| **one file per tile** | **864** | **1,760** | what Phase 3 actually does |

The middle row is what makes this readable rather than arguable: serialization
does not destroy sharing, the *file boundary* does — and that is precisely the
boundary Phase 3 puts between workers. So this is not a `.brep` defect to work
around; it is what one-process-per-tile means.

**The three repairs are all closed off, each by a measurement this project
already has.** Re-identifying the duplicates on the master is a cheap exact
lookup, but *merging* them needs `BRepTools_ReShape` to replace edges **and**
their vertices, which docs/algorithm.md §8 measured coming apart
(`BRepCheck_NotConnected`, invalid solid, wrong volume, while the shell still
"closes"). Sewing only the seam-bearing faces is G8's split, whose production
failure mode is a full unsplit sew of a volume-scaling face set — the one thing
docs/algorithm.md §6 and §8 exist to prevent, at 4 h 45 m of a 5 h 04 m run.
Tiling inside a *single* worker is sound (the middle row proves it) but serial,
and G13 measured serial tiling at an 8 % saving having already concluded "plan
on tiling being worth about `W` and no more; the win is parallelism".

**Threads were then checked, because they would sidestep all of this** — one
heap, tiles pointing at the same objects, nothing serialized. Gate G17 measures
both halves and they cancel exactly:

| | parallelism | tile identity |
|---|---|---|
| processes (`.brep`) | real — `boundary` reaches 5.2 of 6 cores | **destroyed**: 864 free edges |
| threads (shared heap) | **none**: 1.04x on 6 threads | perfect: **0** free edges |

Threads fix precisely what processes break and break precisely what processes
fix. G17's probe A shows why, by mechanism rather than by symptom: a Python
counter running alongside one `unify_same_domain` call retains **3.7 %** of its
solo throughput, so OCP holds the GIL for essentially the whole call and the
tiles run one at a time however many threads dispatch them. (This is G7's
finding, which had been inferred from a 0.91–1.01x speedup, now observed
directly.) `ShapeUpgrade_UnifySameDomain` also exposes no internal parallel
mode — no `SetRunParallel`, no thread-pool hook — so there is nothing to switch
on either.

**So sub-body parallelism here is closed, not merely unbuilt.** There is no
third transport in this architecture. The only thing that would change the
answer is free-threaded CPython (PEP 703), and the project pins Python 3.11
with pinned wheels against which no free-threaded `cadquery-ocp` build exists;
a C extension releasing the GIL is ruled out by specification.md §2's
no-compiler, ordinary-wheels packaging. Shared memory is not an alternative
either — it addresses I/O, which was never the cost (940 MB against an
18-minute stage), and a `TopoDS_Shape` is a graph of pointer-linked C++ objects
that another process cannot use at a different base address, which is why the
identity loss is structural rather than a consequence of choosing files.

**Do not retry this without first disproving G15 and G17**, which take one
command each and about a minute together.

#### The input-side alternative: exact, and worth nothing

With tiling closed off, the other way to attack the same stage is to give the
kernel less to do. Phase 2 makes that *provable* rather than heuristic: since
the interior is built pre-merged, an interior face has no interior partner left
— within a junction its coplanar neighbours are not adjacent, and across a
mid-strut interface the two half-faces are already one — so the only merge it
can still take part in is with a boundary-derived face it touches. And which
faces are boundary-derived is known by construction: they are the objects
`weld.assemble` added. So the face merge can be restricted to the boundary layer
plus one hop, with the rest carried into the result by reference, exactly and
with no search.

**It was built, measured, and reverted.** The correctness half worked perfectly:

| | `dense-lattice` | rehearsal |
|---|---|---|
| faces reaching the kernel | 20,494 of 25,234 (−19 %) | **375,489 of 690,997 (−46 %)** |
| output | 15,966 faces, 67,898 edges, 52.67 MB | 584,028 faces, 2,517,881 edges, 2.00 GB |
| against the unrestricted run | **byte-identical** | **identical in every figure**, volume drift 1.60e-07 both |
| reassembly fell back | 0 | 0 |

So the restriction loses **no merge at all** — the proof about which faces can
merge is sound, and the reassembly by shared topology held on the real part,
checked with §8's every-edge-twice proof rather than trusted.

**And it is not faster.** `simplify` measured 19 m 06.6 s against the
post-Phase-2 18 m 07.9 s. That comparison is *not* a controlled pair — the
untouched stages moved between +0.1 % (`classify`) and +17.6 % (`export`) — so
the honest reading is "no measurable win", not "a 5 % loss". But the mechanism
was then measured directly, on one solid in one run, and it is unambiguous:

    cutting the face merge's input 20 % (23,236 -> 18,660 faces)
    cuts its time 6 %   (1.515 s -> 1.420 s)

**An elasticity near 0.3 — and the reason is not the one it looks like.** The
obvious reading is that `ShapeUpgrade_UnifySameDomain` is priced by what it
emits rather than what it consumes, so a smaller input buys little. **That
reading is wrong, and gate G16 disproves it**
(`tools/prototypes/g16_unify_elasticity.py`): on a *generic* subset the kernel
prices its input almost exactly linearly, at a mean elasticity of **0.98** over
removals of 10–60 %. Remove a representative 40 % of the faces and you save
very nearly 40 % of the time.

So the kernel is not letting the restriction down. **The selection is.** A
*correct* restriction skips exactly the faces unification would have returned
unchanged — which are exactly the cheap ones — and keeps exactly the faces that
merge, which are the expensive ones. It is therefore self-defeating by
construction: the better it is at removing only inert faces, the less of the
cost it removes with them. 0.98 generic against ~0.3 targeted is the size of
that effect, and no implementation can improve on it, because the property that
makes the restriction *correct* is the same property that makes it worthless.

Against that ~0.3, the restriction costs a *linear* ~0.045 ms/face of
bookkeeping — the adjacency index plus the closure proof, ~45 s at rehearsal
scale — which is the same order as the ~1 min the elasticity predicts it could
save. There is no version of this that wins.

This closes the input-reduction direction rather than leaving it open, and it
generalises past this one attempt: **any** proposal to feed this stage less must
skip only faces that would not have changed, so it inherits the same 0.3. Phase 2
above already recorded that "`simplify` scales with its output, not its input"
and that removing input buys less than face counts suggest; the mechanism here
is the precise version of that observation, and it is a sharper claim — the
earlier one blamed the kernel, and G16 shows the kernel is linear. The
implementation itself is not in the tree; it is recorded here because what it
cost to learn is the useful part, and the same is true of Phase 1.

#### Where `simplify` stands now

Still the largest single stage, still at ~0.98 cores, and now with both obvious
levers disproved: it cannot be spread below the body — not by processes, which
break tile identity (G15), and not by threads, which do not run in parallel at
all (G17) — and it cannot be *selectively* fed less, because the faces a correct
restriction may skip are the cheap ones (G16).

What is left is the direction Phase 2 already took: **have fewer faces exist at
all**. G16's 0.98 says the cost really is close to linear in the faces this
stage handles, so anything that reduces that count upstream converts almost
one-for-one — which is exactly why Phase 2, which removed 44.8 % of the interior
before it was ever built, is the one thing on this list that paid. The
distinction that matters, and that took two failed attempts to see clearly, is
between *not building* a face (works) and *not looking at* one that exists
(does not).

**There is no open plan here, and no obvious next candidate** — Phase 2 already
took the interior close to its 50 % ceiling, and the remainder is the boundary
layer, whose faces come out of booleans and are not ours to choose. Anyone
picking this up should start from that distinction and from G16, not by
proposing a new way to parallelise or restrict the call.

#### The gate every phase was held to

`python -m pytest test -q`, `python tools/e2e.py` with both golden samples at
0 mm³, then one full rehearsal under `profile_run.py`. Phases 1 and 2 landed
through it; Phase 3 never reached it, and the input-side alternative passed
every correctness check in it and was reverted on the performance measurement
alone. Worth keeping as the shape of gate this stage needs: **the rehearsal is
the only place any of these could be judged**, since `dense-lattice` is
boundary-dominated and a fortieth of the size.


### `material_outside` reported the whole lattice as outside — FIXED 2026-08-18

**Found 2026-08-17** while verifying the `cc=12, t=2.5` output by hand after
fixing the two guards below, and fixed the next day. A harness defect, not a
generator one — but it meant §6.2's "no generated material lies outside the input
body" could not be run at production scale, which is where it is worth most.

**The symptom.** One `BRepAlgoAPI_Cut` of the whole output against the input body
reported **354,733 mm³** outside — the entire lattice, on a part every junction
of which was *built* by intersecting with that body.

**The mechanism is not the obvious one, and the obvious one would have produced a
detector that never fires.** This entry first recorded the cut as having "quietly
did nothing" and returned "its own argument untouched". It did not. The call
reports `IsDone`, `HasModified` **and** `HasGenerated`, returns **43,672** faces
where 43,530 went in, and removes 3.6e-03 mm³. It ran, re-partitioned the solid,
and mis-classified nearly all of it — so a test for "the volume came back
unchanged" would have missed it, the relative difference being **1.02e-08**
rather than zero. The input is the classic ill-conditioned case for a boolean: a
lattice trimmed from a body has a large share of its faces lying *exactly on*
that body's surface.

**It is a scale limit, and it is bracketed.** Cutting the same output *solid by
solid* gives **exactly 0 mm³** for eight of the nine, at 18 to 68 faces each;
only the 43,530-face body misbehaves. The committed scenarios top out at 15,966
faces (`dense-lattice`), which is why this was never seen.

**The output is sound**, established without any boolean: meshing that body at
0.05 mm and classifying all 55,513 distinct surface points against the input with
`BRepClass3d_SolidClassifier` puts **1** point outside at 1e-06 mm and **0** at
1e-05 mm — that point lies *on* the input surface, two orders inside the 8.7e-04
to 1.5e-03 mm tolerances the part's own faces carry (docs/algorithm.md §8).

**The fix, in three parts.** `material_outside` now cuts **per solid** and sums;
it **contradicts** rather than trusts, checking any non-zero remainder against
`surface_points_outside` — which uses no boolean — and reporting the solid as
*unmeasured* when the two disagree; and it returns a dict, so "not measured" has
somewhere to live that cannot read as "measured as fine". `tools/e2e.py` reports
the exact sum over the solids it could measure plus a separate, explicitly
weaker containment check for the rest, mirroring how `golden_check` already
labels its own fallback. `tools/smoke_bundle.py` instead *fails* on an
unmeasured solid, because at the ball's scale that would mean something is wrong
with the bundle rather than with the check.

A real defect survives all of this: material genuinely outside the body shows up
in both tests, so the contradiction path cannot swallow one.

### Two guards that refused valid input — FIXED 2026-08-17

`-cc 12 -t 2.5` on `TD_HX_rehearsal_test` failed, twice in succession, on two
unrelated guards in two different stages. Both parameters are inside §3's ranges
and satisfy `t < a`, so both failures are the one thing docs/algorithm.md §11
says is never acceptable: refusing correct input. Neither was a regression —
both bars had been there since the checks were written, and this part is simply
the first to reach them.

**They are the same mistake at one remove from each other**, which is why they
are recorded together: a scalar proxy standing in for "did the geometry change",
with its bar expressed in units that do not match what the number actually
varies with. The fix differs because the answer to "is there a cheap exact test
instead?" differs.

**First, the pinhole repair's volume bar** (`PINHOLE_VOLUME_TOL`, 1e-9;
docs/algorithm.md §7, `tools/prototypes/RESULTS.md` G19). Removing a 1.010641e-06
mm pinhole wire "changed" a 77.4 mm³ junction's volume by 1.235e-09 relative.
It changed nothing. `BRepGProp::VolumeProperties` requires a shape "exempt of
any free boundary" and a pinhole wire *is* one, so the pre-repair figure was
never OCCT's to promise. Adding a synthetic open wire to a clean face reproduces
the same footprint to within 1 %, unchanged when the wire is lengthened a
thousandfold and varying tenfold with which face carries it — so there is
nothing the repair controls to size a bar with, and the 2.7e-15 the bar was
calibrated on was two wires that happened to sit on faces contributing nothing.
Three plausible readings — quadrature noise, a coordinate-magnitude artifact, a
bar mis-scaled for junction size — are each disproved by measurement in G19
rather than argued away; the third is the one worth remembering, since the
failing junction is 9x the calibration junction's volume with a *shorter* wire
and drifts seven orders further.

*Fixed by replacing the bar with an exact structural proof*,
`occ.only_inner_wires_dropped`: same face count, every face's outer wire the
same object, same orientations, bit-identical areas, and exactly the accounted-for
wires gone. That pins the enclosed region object-for-object, where volume could
only ever be a biased proxy for it. Area — the check §7 always argued was the
sound one — stays exactly as it was, and is bit-identical throughout.

**Second, same-domain unification's volume bar** (`UNIFY_VOLUME_TOL`, 1e-5;
docs/algorithm.md §9, G20), reached for the first time once the pinhole guard
was cleared. A 181 mm³ floating island drifted 1.381e-05. Here something
genuinely does change — surface area shifts too, and adaptive Gauss-Kronrod to
1e-11 does not converge the two figures, so it is a real re-description and not
the integrator — but the exact symmetric difference, cut both ways, is
**0.000000000 mm³**, and read as a displacement (`|ΔV| / area`) it is
**6.96e-06 mm**, two orders inside the 8.7e-04 to 1.5e-03 mm tolerances OCCT
itself records on the faces being merged. The island's mirror twin, identical to
0.1 % in volume and area, drifts 29x less, so the magnitude belongs to the merge
rather than to the part and cannot be predicted from it.

*Fixed by loosening the bar to 1e-4*, per the user's decision, which across the
1.5–3.7 per mm surface-to-volume ratios that run spans admits at most ~3e-05 to
7e-05 mm of boundary movement — still over an order of magnitude inside those
face tolerances. No exact structural test was available to replace it with: the
symmetric-difference boolean that settles the question takes seconds on a
99-face island and does not finish quickly on the 69,305-face body beside it.
The guards §9 always named as the stronger pair — the solid-count check and
`BRepCheck_Analyzer` — are untouched.

**Result.** `python src/main.py -i test/TD_HX_rehearsal_test.step -cc 12 -t 2.5
--cores 6` completes and writes a valid STEP. `python -m pytest test -q` and
`python tools/e2e.py` stay green, both golden samples at 0 mm³; the committed
scenarios never reach either bar, so neither change alters their output.

**The lesson, added to the four G9–G12 already carry.** Those were repairs aimed
at the wrong object. These are *checks* aimed at the wrong quantity —
docs/algorithm.md §5.1 already warns that a gate is only as trustworthy as the
tightness of what it compares, and these add that it is only as trustworthy as
whether that quantity is defined, and dimensionally right, on the shape it is
handed. Where a check can be made structural, it should be.

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

### Scale rehearsal, chapter closed: paths 1–4 implemented and re-measured

**First run 2026-08-14**, **re-profiled 2026-08-15** after implementing paths
1–4 below, both on `TD_HX_rehearsal_test.step` at `cc=5, t=1`,
`--cores 6 --ram 20 -bg` on the 6-core / 32 GB development workstation. Both
runs used a temporary, uncommitted bypass of the `assemble`-stage watertightness
gate (the defect closed above as "Micron-scale debris edges" — it was open and
out of scope at the time) so every later stage could be measured; neither run's output is a
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

> **Superseded as a baseline, kept as a record.** Both columns below were
> measured with the `assemble` watertightness gate bypassed and before the
> G9–G12 defect family was fixed, so neither is a figure current work should be
> compared against — the 47.1-minute total in particular has been withdrawn
> wherever it was used as an anchor. The current profile, measured as a
> controlled pair with no bypass, is in
> [profiling-reports.md](profiling-reports.md): **55 m 17.9 s** before Phase 2
> and **51 m 43.3 s** after. What this table is still good for is the *shape* of
> the chapter it documents — which of paths 1–4 paid and which did not.

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
