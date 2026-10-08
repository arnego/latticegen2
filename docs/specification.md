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

All are implemented in [`tools/e2e.py`](../tools/e2e.py).

**`SpiralTest.step` is the hardest committed part.** At 2,073 boundary
junctions it is a fraction of `dense-lattice`, but its swept B-spline surface
makes the boolean fit trim curves carrying recorded tolerances an order above
anything the other committed parts produce, and it is the only committed part
whose dominant component tiles. It is also the only *scenario* that exercises
§7.2 — the intersection that returns its own operand (8 junctions), the
local-block re-trim (10), and a junction dropped for containment (1) — which is
why the unit suite uses it heavily (`test_boundary.py`, `test_export_truth.py`)
as well.

**§7.2 is not peculiar to it, though, and reading it that way understates
where that repair matters.** `TD_HX_rehearsal_test.step` at `cc=5, t=1` returns
**35** junctions untrimmed and re-trims every one of them per half-strut. That
part is not an e2e scenario — it is far too slow — so the committed scenarios
are the only automated coverage §7.2 has, and they carry one part's worth of it.

Two properties of the toolchain are what let this part complete, and both are
load-bearing rather than incidental to it. Its worst **input seam gap** — the
distance between the two faces meeting at one of its own edges — is
6.8602e-03 mm, low enough that trimming does not cut the region into features
smaller than the gap. And the STEP writer declares the *greatest* tolerance a
shape carries rather than the average of them (§11), without which its dominant
body loses 2 edges to the file under-declaring what the geometry needs.

| Scenario | Parameters | Expected result |
|----------|-----------|------------------|
| smoke-fast | -i test/80mm-test-ball.step -cc 20 -t 4 --cores 4 | generation < 10 minutes. **Measured: 6.4 s.** |
| smoke-verified | -i test/80mm-test-ball.step -cc 20 -t 4 --cores 4 | valid STEP, generation < 20 minutes, matching golden sample test/80mm-test-ball-cc20t4-golden-sample.step. **Measured: 6.3 s, symmetric-difference volume 0.0000 mm³.** |
| dense-lattice | -i test/test-cylinder.STEP -cc 10 -t 1.5 --cores 6 | valid STEP, no self-intersections, matching golden sample test/test-cylinder-cc10t1.5-golden-sample.step, generation < 10 minutes. **Measured: 47.5 s, symmetric-difference volume 0 mm³.** |
| spiral-stress | -i test/SpiralTest.step -cc 5 -t 1 --cores 6 | valid STEP, no golden sample, generation < 20 minutes. **Measured: 11 m 01 s**, both solids `BRepCheck_Analyzer`-valid with 0 non-manifold edges after a round trip. Containment falls back to point sampling — the dominant body is past `verify_geometry.CUT_MAX_FACES` — and is reported as the weaker check, never as an unmeasured pass. |
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

**Nothing is open here.** The last item — `TD_HX_rehearsal_test` at `cc=5, t=1`
being refused at `export truth` — was closed 2026-08-26 and is the first chapter
of §11, together with the two guards found while closing it. Every committed
part now writes its output at every parameter set this project has run.
`docs/testing.md` names the cheapest inner loop for anything touching
`validate`, which is that part at `cc=7, t=1.4` (~20 minutes) rather than at
`cc=5, t=1` (~91).

## 11. Closed — kept for the reasoning, not as work

The closed chapters are in [specification-closed.md](specification-closed.md).
That file is deliberately not loaded with this one, so it is read on demand:
**before changing a stage one of the chapters below covers, read that chapter
first.** Each records what was tried, measured and disproved, and several end
in "do not retry this". A reference elsewhere to "specification.md §11" means
that file.

- `TD_HX_rehearsal_test` at `cc=5, t=1` — the readings were the ruler, not the body
- Two guards that refused sound geometry before `export truth` could speak
- Re-fitting the pcurve cannot keep the body — the gap is in the input file
- A body can be valid here and not describable in the file — the export-truth gate
- `stitch`'s round-2 repair — chapter closed: the fix disproved, the check repaired, the scan parallelised
- `--ram` removed — an accepted-but-unenforced budget, taken out rather than fixed
- Pipeline parallelism between `classify` and `assemble` — chapter closed, two stages won, two proposals disproved
- `simplify` beyond body-for-body — chapter closed, two levers disproved
- `material_outside` reported the whole lattice as outside — FIXED 2026-08-18
- Two guards that refused valid input — FIXED 2026-08-17
- Invalid boundary faces from grazing trims — FIXED 2026-08-17 (34 → 0)
- Micron-scale debris edges from near-tangential trims — FIXED 2026-08-16
- Scale rehearsal, chapter closed: paths 1–4 implemented and re-measured
