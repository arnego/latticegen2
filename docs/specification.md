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
| --orient | 3 × float | optional | degrees | -360 - 360 each | 0 0 0 | Print orientation `<rx> <ry> <rz>` of the input geometry relative to the print bed (§4.5). Orients the lattice to the bed. Usable with or without `--support`. |
| --support | flag | optional | NA | NA | disabled | Generate print supports for inner overhangs (§4.6). |
| --support-thickness | float | optional | mm | 0.2 - 20 | `t` | Minimum support thickness. Requires `--support`. |

**The overhang limit is not a parameter.** Support faces never overhang more
than the lattice rod angle `arcsin(sqrt(2/3))` ≈ 54.74° from vertical, and the
input's surface is supported where it overhangs more than 60° (§4.6). Both are
fixed; `--overhang` existed through the support skeleton and was removed with
the support algorithm, and passing it is refused (exit 2).

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
 - Volume of the end product, split by origin: the lattice's (`lattice_volume_mm3`), what supports add to it with `--support` (`support_volume_mm3`, the sum of the support bodies, which never overlap the lattice), and the total written (`output_volume_mm3`)
 - Maximum memory usage
 - Path to output .step file 

Upon failure, the script shall output a human readable reason for the failure, e.g.: parameter bounds exceeded, issues with input geometry, issues with resarouces from the run-time system, write or read access issues, etc.

A run that succeeds with an output that may be wrong, incomplete or degraded says so in a warning, on the console and in the `.log` regardless of `-v`: one message starting `WARNING: `, stating what went wrong, what it may cause, and where to inspect in the part's coordinates. A repair that left the output as intended is a `note:`.

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
* **Print orientation** — always shown, because the orientation applies with or
  without support (§4.5): three spin controls, `rot X`, `rot Y`, `rot Z`, in
  degrees with 5° steps, beside a small 3D view of the print bed. The view draws
  the bed as a gridded plate, the part's X/Y/Z axes, and the **input geometry
  itself standing on the bed** the way the three angles say it will. Dragging
  the view turns the part and the three boxes follow; editing the boxes turns
  the part.
* **Generate support** — a tick box, with **Min thickness** (0.1 mm steps)
  beside it. The field follows `t` until it is edited, which is the command
  line's default made visible; it is always shown and is **greyed out while
  the box is unticked**; unticked, no support argument reaches the run, which
  then loads no support code (§4.6). With the box ticked, every part of the
  previewed surface that overhangs more than 60° at the current orientation —
  what will be supported — is **tinted**, and the tint is recomputed on every
  change of orientation — boxes, steppers and mouse drag alike — so it can be
  checked before an hour-long run is started.
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

**How the input gets drawn.** The window loads no geometry kernel, so it cannot
read a STEP file. It starts a short-lived helper — `src/main.py --preview-mesh
<input> <npz>`, transport like `--progress-stream` and equally absent from §3's
table — which tessellates the input coarsely (at most 4,000 triangles) and
writes plain arrays; the view shows the axes alone until that arrives. The
preview is a `tkinter.Canvas` and a painter's sort, so it adds no dependency.
While the mouse is down a reduced copy is drawn so the drag stays smooth. The
tint is an indication from that coarse mesh; a run measures overhang on the
true surface (docs/algorithm.md §14.1).

**While a run is in flight** the fields go inactive but stay readable, `Start!`
becomes `Stop!`, and two bars appear:

* the **top** bar names the running stage and fills according to a fixed share
  per stage, declared in [`gui/weights.py`](../src/latticegen2/gui/weights.py).
  **There are two layouts, chosen by whether the run generates supports.** A
  supported run goes through three more stages — `overhang`, `support` and
  `pockets`, between `simplify` and `validate` — and each gets a band of its
  own. The plain layout comes from the plain run of the 2026-10-11 controlled
  pair in docs/profiling-reports.md, the supported layout from the 2026-11-06
  `rehearsal-support` run, every stage of it; there `pockets` alone is 53 % of
  the bar. Both are computed once when the window starts. A run announces its stages when it
  begins, and that list is what selects the layout;
* the **second** bar shows work within the stage — `boundary trim: 9,776 /
  19,552` — with the text drawn over it.

Beneath them, one line of resource data: peak memory and elapsed time. Peak
memory is genuinely all that is measured continuously (§3's summary list), and
the line says only that rather than padding itself out.

Beneath *that*, with **verbose** ticked, the run's own output in a scrolling
pane — every line of it, which is what the `.log` gets.

**Unticked, the pane holds only what the user has to know about:** the run's
warnings (below), kernel chatter on stdout, and anything at all on stderr, which
is where a failure's one reason line lands (§7). It is hidden entirely when that
is empty, so an ordinary run does not carry a blank box around for the hour it
takes — and the box appearing is itself the signal that something spoke up.

**A warning is any line starting `WARNING: `, and it always gets through.** It
means the output may be wrong, incomplete or degraded, and no widget says so. A
warning opens the pane whichever way the box is set, and the backlog the window
keeps never drops one, so an early warning is still there at the end of an hour.
Every warning takes one form (docs/algorithm.md §10): what went wrong, what it
may cause, and where to inspect in the part's coordinates.

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

### 4.5 Print orientation

The orientation of the struts and of the junction template is
determined by the orientation of the **print bed**, not by the coordinate system
of the input geometry. Everything §4.1 says about the Z-axis holds about the
bed's normal (the *build direction*) instead.

- **Input:** three angles in degrees, `rx ry rz`, rotating the input geometry
  relative to the bed about the bed's fixed X, then Y, then Z axis:
  `R = Rz·Ry·Rx`, `p_bed = R·p_part`, bed normal `+Z`. The build direction in
  the part's own coordinates is therefore `b = Rᵀẑ`.
- **Independent of support generation.** It can be set with or without
  `--support`.
- **No identity bypass.** The rotation is applied on every run, `0 0 0`
  included. An unrotated run must still reproduce the golden samples.
- The input body is never transformed, and the output stays in the input's
  coordinate system (§1).
- **File and part name:** `-rotx<rx>y<ry>z<rz>` (part name `+rotx…`) is added
  after `t<t>` when any angle is non-zero.

### 4.6 Print support generation

An opt-in stage that creates support structures, as used for
3D printing, for the inner overhangs of the input geometry.

- **Bypassed entirely unless enabled.** No support code is loaded by a run
  without `--support`.
- **Process.** Metal powder-bed fusion.
- **Detection.** The input geometry is analysed where a surface **or edge**
  overhangs more than **60° from vertical** *on the inside of the input body* —
  the regions of its surface acting as a ceiling of the volume, where the
  enveloping shell overhangs the lattice core. Outside overhangs are neither
  detected nor supported. The lattice itself is not analysed.
- **Two angles, both fixed.** Support surfaces never overhang more than the
  lattice rod angle, `arcsin(sqrt(2/3))` ≈ 54.74°. The input's walls and roofs
  are supported only past 60°, so surface between 54.74° and 60° is left bare
  by decision: a support ending on a ceiling at the rod angle would meet it
  tangentially, and ending at 60° keeps it at least 5.26° off the ceiling. A
  planar face tilted 55–60° is left unsupported entirely, and a curved one
  keeps a bare band (3.7 mm wide on the 80 mm ball). Neither angle is a
  parameter.
- **Full coverage.** Every point overhanging more than 60° is backed by support
  or lattice material, all the way to the edge of each overhang: every ceiling
  point lies in a support or a strut. **The one exception is a support piece
  the geometry kernel cannot trim to the input** at any fuzzy value, after a
  STEP round trip, or in halves down to a single support node: it is dropped,
  as a lattice junction no trim produces is (docs/algorithm.md §7.2), and the
  run reports it regardless of `-v` — by name and position, the support
  volume lost, and the ceiling points it leaves unbacked
  (docs/algorithm.md §14.3). A run also warns, with coordinates, of any ceiling
  point nothing in its output backs. Chips of about 1e-3 mm³ that trimming leaves
  where a ceiling grazes a support node are discarded rather than shipped as
  floating micro-solids; the harness reports the bare samples beside them as
  slivers, within 0.05 mm of a support (§6.2).
- **Anchoring and design depth.** A support runs from the ceiling down to the
  nearest lattice rod below it, or, where no rod is reached, to the enveloping
  wall — the nearest anchor by drop, searched down the three rod directions —
  so the lattice is generated first.
- **Pipeline position.** Supports are generated against the **completed
  lattice**: after `simplify` and before `validate`. The generator receives the
  lattice solids, the input geometry, `--orient`, `--support-thickness` and the
  inner overhangs already detected, and returns support solids only. One
  validity pass — including the export-truth gate — then judges the lattice and
  the supports together (docs/algorithm.md §14.2).
- **Separate from the lattice, manifold bodies.** Supports are **not** merged
  with the lattice. They share faces with the lattice and with the input's
  walls, floors and roofs, so that the merge into the enveloping part in the
  downstream CAD tool succeeds. Supports that would touch are **one body**, and
  a pocket fill is one body with the supports that seal it. Every support body
  is a manifold solid: no edge is shared by more than two of its faces, beside
  a strut as anywhere else, since a CAD kernel reads such a body as no solid.
  Support bodies never overlap the lattice or each other, and a run in which
  one would, or in which a body is not manifold, fails rather than ships.
- **Enclosed voids are filled.** Powder must not be trapped: a pocket of free
  space that supports seal off is filled with support material. Where the
  geometry kernel cannot cut a part of the support envelope even after splitting
  it, that part is left unsearched rather than the run failing, and the run
  warns regardless of `-v` that enclosed pockets may remain there, giving each
  region's x, y and z range (docs/algorithm.md §14.3 S7).
- **Self-supporting.** No support surface or edge exceeds the rod angle.
- **Thickness.** No free-standing support member is thinner than the minimum
  support thickness, which defaults to the strut size `t`. The support grid divides the lattice cube as finely as still leaves at least that thickness
  between parallel rod faces.
- **Volume.** Each support is the box of support nodes between its anchor and
  the ceiling point it backs: it reaches only as far as its nearest anchor, and no
  higher than the ceiling it serves.
- **Integral.** Once merged downstream, supports are a permanent part of the
  produced body, with solid interfaces and no feature designed to ease removal.
- **Shape.** Supports are staircases of support nodes on a grid of the lattice's
  own frame: every face lies on a plane of that frame — the grid's planes run
  through the rods' centre axes, so a support standing on a rod has its lowest
  edge on the rod's axis, and a face meeting a rod lies on the rod's own face —
  except where it lies on the input's own surface, and coplanar faces are
  merged.
- **Rooted.** Every support stands on a rod, a wall or another support: no
  corner of it hangs over powder.
- **File and part name:** `-sup` (part name `+sup`) is added last when support
  is enabled.
- **Parallelisation** is used wherever a run-time reduction above 5 % is likely.

The normative algorithm is docs/algorithm.md §14.

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

  **A support body is not a floating body.** Supports are separate bodies
  by design (§4.6), each attached to the ceiling it supports and to its
  anchor, so every support body is exempt from this rule, whether it is
  anchored on a wall or on the lattice. Supports are generated after the
  lattice's own floating-body rule has run. A fragment the support generator
  cuts off while building a body is construction debris, discarded inside the
  generator; it is never an output body, so the rule below is unaffected.
  
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

**`SpiralTest.step` was removed on 2026-10-30, with its `spiral-stress`
scenario and the unit tests that read it.** It was the only committed part whose
dominant component tiled and the only scenario exercising §7.2 — the
intersection that returns its own operand, the local-block re-trim, and a
junction dropped for containment — so **§7.2 now has no automated coverage**
beyond the box-only unit test that a junction wholly inside is checked and kept.
`TD_HX_rehearsal_test.step` at `cc=5, t=1` re-trims 35 junctions per
half-strut, but it is far too slow to be a scenario. The measurements made on
`SpiralTest` stay in docs/algorithm.md as the record they are.
`test/spiral-island-unwritable.brep`, cut from an older version of that part,
stays: it is the export-truth gate's one body that cannot be written, and reads
no `SpiralTest.step`.

| Scenario | Parameters | Expected result |
|----------|-----------|------------------|
| smoke-fast | -i test/80mm-test-ball.step -cc 20 -t 4 --cores 4 | generation < 10 minutes. **Measured: 6.4 s.** |
| smoke-verified | -i test/80mm-test-ball.step -cc 20 -t 4 --cores 4 | valid STEP, generation < 20 minutes, matching golden sample test/80mm-test-ball-cc20t4-golden-sample.step. **Measured: 6.3 s, symmetric-difference volume 0.0000 mm³.** |
| dense-lattice | -i test/test-cylinder.STEP -cc 10 -t 1.5 --cores 6 | valid STEP, no self-intersections, matching golden sample test/test-cylinder-cc10t1.5-golden-sample.step, generation < 10 minutes. **Measured: 47.5 s, symmetric-difference volume 0 mm³.** |
| oriented-lattice | -i test/80mm-test-ball.step -cc 20 -t 4 --cores 4 --orient 35 20 -110 | valid STEP through every §6.2 check with the lattice turned to the bed; the orientation recorded in the log and no support stage run. No golden sample. |
| support-gaps | -i test/support-gaps.step -cc 10 -t 1.5 --cores 4 --support | valid STEP through every §6.2 check, on a small part (`tools/make_support_gaps_part.py`) holding the two ways `rehearsal-support` left ceiling bare: a keel no face of which overhangs past 60°, and a strut tip dropped as a floating body from a node still kept. Both are checked to be present, so it cannot pass by no longer testing them. Before the fix: 16 of 1,836 face samples and 32 of 33 keel samples unbacked. **Measured: 21 s, none unbacked**, and none at 0.1 mm sampling either (42,690). |
| support-clipped | -i test/support-clipped-{cone,cylinder}.step -cc 7 -t 1.4 --cores 4 --orient 180 0 0 --support --support-thickness 1 | valid STEP through every §6.2 check on two small cuts of `TD_HX_rehearsal_test.step`, kept in its coordinates, where a ceiling just past the bar clips support nodes no refined mesh vertex lands in. Coverage is also asked at 0.1 mm, since 0.5 mm can step over such a sliver. Before the fix: 16 and 481 samples bare at 0.1 mm. **Measured: 13 s and 21 s; none bare, 1 sliver within 0.05 mm** beside the discarded chips, 11 samples on a rim. |
| support-tilted | -i test/80mm-test-ball.step -cc 10 -t 1.5 --cores 4 --orient 35 20 -110 --support | valid STEP through every §6.2 check, the support stages running between `simplify` and `validate`, and the support checks of §6.2: a curved overhang, tilted. No golden sample. **Measured: 107 s**, one support body of 5,296 mm³, 12,524 of 12,524 overhang samples backed, no support corner hanging free. |
| support-cylinder | -i test/test-cylinder.STEP -cc 10 -t 1.5 --cores 6 --orient 0 -90 0 --support | valid STEP, the support checks, the same run without `--support` contained in the output, generation < 20 minutes. **Measured: 235 s**, three support bodies of 15,206 mm³, 33,065 of 33,065 overhang samples backed, no support corner hanging free, 0 mm³ of the unsupported run outside the output. **Prepared for a golden sample, which does not exist yet:** the comparison against `test/test-cylinder-cc10t1.5-rotx0y-90z0-sup-golden-sample.step` reports `SKIP` until an inspected output is committed under that name. |
| rehearsal-support | -i test/TD_HX_rehearsal_test.step -cc 7 -t 1.4 --cores 6 --support --orient 180 0 0 --support-thickness 1 | **Opt-in: runs only when named with `--only rehearsal-support`**, and only for final verification. Print support at production scale, the part turned over: valid STEP through every §6.2 check including the support checks, generation < 12 hours on 6 cores. **Measured 2026-10-31: 1 h 40 m, and coverage failed** — 630 of 125,777 overhang samples unbacked. The causes are fixed and verified on `support-gaps`. **Re-run 2026-11-02: failed at exit 4 after 8 h 22 m**, a support body overlapping the lattice by 1.475e-02 mm³; no file written. **Run again 2026-11-02 after those fixes: passes every §6.2 check. Measured: 2 h 08 m**, 4 h 56 m with the harness's checks; 27 support bodies of 40,136 mm³, 125,777 of 125,777 overhang samples backed, no overlap, 3 pocket regions left unsearched and reported. **Run 2026-11-05 with the local pocket search: passes every §6.2 check. Measured: 1 h 42 m**, 4 h 31 m with the harness's checks; the same 27 support bodies and 40,136 mm³, all overhang samples backed, no overlap, 2 pocket regions left unsearched and reported. **Run 2026-11-06 with the support grid through the rod axes: passes every §6.2 check. Measured: 1 h 46 m**, 4 h 45 m with the harness's checks; 52 support bodies of 32,018 mm³, all 125,777 overhang samples backed, no overlap; the generator warns of 2 of 225,637 ceiling points unbacked and 4 pocket regions unsearched. **Run again 2026-11-06 with those warnings traced to detection and to the kernel, and fixed: passes every §6.2 check with no support warning. Measured: 1 h 43 m**, 4 h 39 m with the harness's checks; 48 support bodies of 31,618 mm³, all 225,547 ceiling points and all 125,777 overhang samples backed, no overlap, no pocket region left unsearched. **Prepared for a golden sample, which does not exist yet:** the comparison against `test/TD_HX_rehearsal_test-cc7t1.4-rotx180y0z0-sup-golden-sample.step` reports `SKIP` until an inspected output is committed under that name. |

| invalid-input | -i test/80mm-test-ball.step -cc 5 -t 4 (strut size `t` >= cell edge `a=cc/√2`) | exits 2, no `.step` or `.log` file written, one human-readable reason line. **Passes.** |

**The hole-block parts are not e2e scenarios.** `test/Hole-block-basic.STEP` and
`test/Hole-block-chamfered.STEP`, at `--orient 90 0 0`, are the support
generator's first test parts, run **without a lattice** — their only anchors are
their own walls — through `tools/support_dev.py`; the `pytest` suite uses them
for ceiling-point detection. The product cannot run them that way, since it
always generates the lattice.

### 6.2 Automated pass/fail checks
For every scenario the harness must verify, without human intervention:
- Process exits with expected console output.
- STEP file is written and non-empty.
- STEP file parses back successfully (round-trip read).
- Geometry is a valid closed manifold solid (no open edges / non-manifold edges).
- **Geometry passes OCCT's exact B-rep validity check** (`BRepCheck_Analyzer`) —
  an exact test on the B-rep itself, not an inference from a tessellation.
- No self-intersections. A crossing found between two triangles at the default
  deflection is re-asked at an eighth of it before it is believed: where a strut
  ends on a curved wall, the wall's chords sag below the true surface and can
  cross the strut's own triangles, which a finer mesh resolves and real
  self-intersection does not.
- **No generated material lies outside the input body** (boolean cut of output
  against input leaves ~zero volume) — a direct check of §1's "fits exactly
  within the user's boundary geometry", independent of any golden sample.
- **The shipped file's pcurves still agree with their own 3D curves.** Read the
  output back with OCCT's own reader and measure, per edge/face pair, the exact
  distance between the two representations against the area of the face carrying
  it. This is asked of the *artefact* rather than of the process that wrote it,
  which is the only version the downstream tools see.
- **With `--support`: every inner overhang of the input is backed by
  material in the output.** The input's surface overhanging more than 60° is
  sampled afresh at 0.5 mm, and every sample, moved a hair into the body, must
  lie inside the *written file's* support or lattice material — independently
  of how the supports were placed or built. A sample whose probe leaves the body
  is on a rim a wall continues, and one within 0.05 mm of output material is a
  sliver (§4.6); both are reported, and neither fails.
- **With `--support`: supports overlap neither the lattice nor each other,**
  and a support face resting on a strut lies on that strut's plane, so the
  downstream merge finds shared faces.
- **With `--support`: no support surface faces down more flatly than the rod
  angle.** Every face of the written file that is neither on the input's own
  surface nor a lattice strut face must be inside 54.74°, read at exact surface
  normals, since a tessellated curved face built at the limit can read a hair
  past it.
- **With `--support`: nothing in the finished print hangs over powder.** The input's inner surface past 60° and every
  support body are searched for faces facing down past their limit, keels and
  down-pointing tips; one touching the input's wall, the lattice or another
  support is rooted there and passes. The lattice is not analysed.
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

* **Finish the print support generator** (§4.6, docs/algorithm.md §14.3–§14.4).
  The frustum algorithm is implemented on branch
  `custom-designed-support-feature`; every support scenario passes, and
  `rehearsal-support` passes at production scale since 2026-11-02, and since
  the second run of 2026-11-06 with no support warning (1 h 43 m;
  docs/algorithm.md §14.4), its stage times in `gui/weights.py`. What remains:
  the overlap check of docs/algorithm.md §14.2 dispatched across the worker
  pool — ~35 min at one core on that run, its largest serial piece — and the
  trims likewise if a measurement shows more than 5 %; and the free-standing
  thickness gate. `rehearsal-support` does not run the final-print
  check (§6.2's "hangs over powder"), so the axis grid's freedom from free tips
  is verified on the four small support scenarios only.
  Golden samples for the three support scenarios wait for an inspected output.
* **The window has not been looked at on a screen in this session.** Its widget
  states and the orientation view's drawing were checked programmatically (the
  canvas was rasterised to an image), not by eye in a running window.

**Otherwise nothing is open here.** The last item — `TD_HX_rehearsal_test` at `cc=5, t=1`
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
