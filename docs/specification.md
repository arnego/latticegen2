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

### Scale rehearsal: attempted, failed in `sew`, needs re-running

**What happened.** `TD_HX_Indre_Volum.stp` at `cc=5, t=1` was run end to end on
2026-08-12. It reached the stitcher and failed there after **5 h 04 m**:

| Stage | Time | |
|---|---|---|
| tessellate | 2.4 s | 28,654 triangles, measured deviation 0.409 mm |
| classify | 2 m 00 s | 527,425 candidates → 29,375 interior, 19,552 boundary, 478,498 outside |
| boundary | 12 m 36 s | 5 workers; 21,955 pieces, **2,969 junctions produced no geometry** |
| connect | 0.5 s | 183 floating bodies dropped, down to 5.08e-06 mm³ |
| instance | 1 m 07 s | 705,000 faces, 624,492 vertices |
| **sew** | **4 h 45 m** | 21,697 shells in; `1 of 14 stitched shells are not closed` |

So two findings, and the second hid behind the first for five hours:

* **`sew` is 94 % of the run** and single-threaded (~16 % of a 6-core CPU, 2–3 GB
  RSS — algorithmic, not memory-bound). See the stitching item below.
* **The cap-drop decision was one-sided**, so a junction could punch an interface
  hole with nothing behind it. Fixed by resolving interfaces symmetrically
  (docs/algorithm.md §7.1) with a hard check in `connect` (§8), which moves this
  class of failure from hour 5 to minute 15 and names the junction.

**Re-run 2026-08-13 with the interface fix, confirming the second finding.** It
reached `assembly input` in **16 m 27 s** and was stopped there, since `sew` is
unchanged. `connect` reported 122,180 interfaces and exactly three caps the two
sides disagreed about, all in one cluster at nodes `(633,-97,-61)` and
`(633,-97,-62)` — the three unmatched holes the old rule would have opened, which
is `1 of 14 stitched shells are not closed`. Detail in docs/algorithm.md §7.1.

Per-stage cost of the fix, against the first run:

| Stage | Before | After | |
|---|---|---|---|
| classify | 2 m 00.2 s | 2 m 03.6 s | noise |
| boundary | 12 m 35.6 s | 13 m 08.9 s | +4.4 %: every cap-plane face is now tagged and shipped, not only the masked ones |
| connect | 0.51 s | 3.94 s | the cap-area quadrature that decides agreement |
| instance | 1 m 07.4 s | 1 m 08.3 s | noise |

3.4 s to make an entire class of watertightness failure unreachable, and to name
it in minute 15 instead of hour 5.

**Still to do:** re-run end to end once stitching is fixed, and record the full
per-stage table, peak RSS, output file size and write time. Everything before
`sew` now totals 16 m 27 s, so that is the floor the fixed run should approach.

**Why the file size still matters:** STEP *writing* is expected to become the
dominant cost once `sew` no longer is, producing a multi-GB file with ~3 M faces.
If that holds, the bottleneck moves outside this tool — into whether
SolidWorks/Catia can usefully import such a file — and that is a user-level
decision about the output contract, not something to change unilaterally.

**Watch the `simplify` stage there specifically.** Same-domain unification
(docs/algorithm.md §9) runs at roughly 0.24 ms/face, so ~3 M faces implies on the
order of 12 minutes. At `dense-lattice` scale it more than pays for itself by
halving the file that export and the round-trip check then have to handle, but
that trade has only been measured at ~30 k faces. It unifies each solid
independently, so it parallelises across solids if it needs to.

### Fuse junction pairs whose two booleans disagree

**Status: built and unit-tested; not yet re-verified at the scale that found
it.** Implemented in [`boundary.fuse_disagreeing_pairs`](../src/latticegen2/boundary.py)
and wired into the pipeline's `connect` stage, between `resolve_interfaces` and
`finalize_pieces`, exactly as designed below. `BoundaryPiece.caps` and
`.cap_faces` are now keyed by `(node, half-strut)` throughout boundary, connect,
pipeline and weld, as the structural change required.

`test/test_boundary.py` reproduces the failure with a real boolean (a genuine
corner cut off one side's cap via `BRepAlgoAPI_Cut`, not a synthetic uniform
scale — the earlier tests in that file already covered the "declined, kept as
exterior" path, so the new ones target specifically the fuse repair):
`resolve_interfaces` reports the notched cap as mismatched; `fuse_disagreeing_pairs`
fuses the pair into one solid whose volume is the exact sum, re-tagged so every
*other* cap of both nodes survives correctly attributed to its own node (the
`_owning_cap` proximity check earns its keep here — a plane-only
`is_cap_plane_face` test alone cannot tell node A's caps from node B's along an
axis orthogonal to the one separating them); resolved a second time the
disagreement is simply gone, with neither side presenting that cap as a
boundary face any more; and the merged piece assembles into a closed,
orientable shell — `weld.shell_defects` reports zero open and zero misoriented
edges, and `BRepCheck_Analyzer` passes it. `python -m pytest test -q` (160
passed) and `python tools/e2e.py` (all four scenarios, matching golden samples
exactly) both stay green — expected, since none of the committed scenarios
contain a disagreeing pair, so this path is exercised only by the new tests.

**Not yet done:** re-running the `TD_HX_Indre_Volum` rehearsal that found this
bug in the first place, to confirm the three real disagreeing caps at
`(633,-97,-61)`/`(633,-97,-62)` get repaired rather than declined. That re-run
is still blocked on the boundary-sew tiling item below — the rehearsal reaches
`sew` at minute ~35 today and `sew` itself is unchanged — so it cannot happen
until that lands too.

Original diagnosis, kept for reference:

**What happens.** `resolve_interfaces` (docs/algorithm.md §7.1) declines a cap
when the two sides present regions that disagree. On `TD_HX_Indre_Volum` at
`cc=5, t=1` that is 3 caps out of 122,180, all in one cluster around
`[2055.4, -90.0, 969.6]`:

| Cap | Side A | Side B |
|---|---|---|
| `(633,-97,-61)` h3 | present | absent |
| `(633,-97,-61)` h0 | 1.000000 mm² | 0.014613 mm² |
| `(633,-97,-62)` h2 | 0.736809 mm² | 1.000000 mm² |

Declining means both sides keep their cap face. **That degradation is unsound,
and the claim that it "leaves an extra solid rather than a hole" was wrong.**
Where the two caps are the same region it is harmless, but these are *mismatched
partial* caps, so keeping both leaves the overlap as non-manifold material and
the remainder as an unfilled hole. `assemble` reports exactly that: **12 edges on
1 face and 12 edges on 3 faces** — 3 caps × 4 ring edges, twice over.

**The fix: fuse the disagreeing pair with a local boolean.** Where instancing's
exactness argument has broken down because the kernel contradicted itself, fall
back to the kernel's own general operation. It is sound, it produces correct
geometry rather than trading a hole for a sliver, and at three occurrences its
cost is irrelevant — this is nowhere near the volume-scaling path §12 keeps
booleans off.

The alternatives were considered and rejected: failing the run refuses sound
input over a kernel defect, which is the one failure mode docs/algorithm.md §11
rules out; dropping the offending pieces removes material, which §5 reserves for
*floating* sub-`t³` bodies only.

**What it needs.** The work is small except for one structural change:

* `BoundaryPiece` assumes **one node per piece**. Its caps are keyed by
  half-strut id alone and `connect.py` reads them as `(piece.node, h)`. A fused
  pair spans two nodes, so `caps` and `cap_faces` must be keyed by `(node, h)`
  throughout — boundary, connect, pipeline, weld and their tests. Nothing else
  about the piece changes.
* The fuse itself runs *before* `finalize_pieces`, while `faces` plus every entry
  of `cap_faces` still form the piece's complete closed boundary, so each side
  can be rebuilt as a solid, fused with `BRepAlgoAPI_Fuse`, and re-tagged with
  `is_cap_plane_face` against **both** nodes.
* Interfaces are then resolved as usual. The merged piece presents one agreed
  region at the cap that disagreed, so nothing is declined there and the rest of
  the pipeline is unchanged.
* A fuse that returns more than one solid, or that throws, is a hard failure
  naming the junctions — it means the two pieces did not even overlap
  consistently, which is beyond what this repair can honestly fix.

Regression to keep: two pieces whose shared cap disagrees must assemble into a
closed orientable shell, with the edge-use tally clean. The synthetic case that
does *not* reproduce it is two pristine template instances — their caps are
identical, so nothing disagrees. The reproduction needs genuinely mismatched
partial caps.

### Tile the boundary sew

**Status: built and unit-tested; not yet re-verified at the scale that found
it.** Implemented in [`weld.sew_boundary`](../src/latticegen2/weld.py) and its
helpers (`_tile_pieces`, `_tile_edge_length`, `_sew_tiles`,
`_worker_sew_tile`), called from the pipeline's `stitch` stage exactly as
designed below.

`BRepBuilderAPI_Sewing` used to receive the whole interior shell along with the
trimmed boundary pieces, on the assumption that an already-shared shell costs it
nothing. `tools/prototypes/RESULTS.md` G5 disproved that — face count dominates,
and adding one closed 194,400-face shell with zero free edges took a 4,000-piece
sew from 76.5 s to 716.6 s — so the assembly was inverted: the boundary layer is
sewn to itself first, and the interior is then *built onto* the topology that
comes out (docs/algorithm.md §8). The volume-scaling shell no longer enters a
geometric search at all.

What remained was the boundary sew itself, which G5a measured at about `n^1.8`
in piece count: 195.8 s at 8,000 pieces, extrapolating to roughly **20 minutes**
at the rehearsal's 21,955. The sew is already per-component and every component
is independent, so each component whose piece count clears a threshold (1,500
pieces — three tiles' worth) is now split into spatial tiles by lattice-index
block (sized to average 500 pieces each, the largest size G5a measured cheap);
each tile is sewn on its own, in parallel across the run's worker processes via
the same `.brep` round-trip `boundary.py` already uses for the trim stage; and
only then are the tiles' results sewn together. A component below the threshold,
or too compact to produce more than one tile, sews exactly as it did before this
existed — confirmed by the two committed e2e scenarios, whose piece counts (a
few hundred to ~1,100) never clear it, so their output is byte-for-byte
unaffected and both still match their golden samples exactly (0 mm³ symmetric
difference).

**Measured, not assumed (`tools/prototypes/RESULTS.md` G6), on the same real
trimmed pieces G5 used, at two scales.** The saving is real but bounded: round 1
shrinks with tile count roughly as the `n^1.8` model predicts, but round 2 sews
shells whose *combined* face count equals the untiled input's, and G5b already
found that sewing pays a face-count cost even where there is nothing to merge —
so round 2 does not shrink to match round 1, and in fact grows slightly as
tiles get smaller. Best measured: **1.45×** at 4,000 pieces / 8 tiles of 500
(round 1 11.0 s + round 2 51.8 s = 62.8 s against a 91.3 s baseline) and
**1.43×** at 8,000 pieces / 8 tiles of 1,000 (27.8 s + 131.8 s = 159.6 s against
228.3 s) — real at both scales, but round 2 alone is more than half the
baseline either way, so there is a shallow optimum around a few hundred to
~1,000 pieces per tile rather than a runaway win from finer tiling, which is
why the tile target is pinned inside that plateau rather than pushed as small
as possible. Production round 1 also runs in parallel across workers, which
this serial measurement does not credit, so the real saving should exceed what
is measured here — by how much is exactly what the rehearsal re-run below would
confirm.

**Not yet done:** re-running the `TD_HX_Indre_Volum` rehearsal, both to see the
tiled boundary sew's wall time at the scale that motivated it (tens of thousands
of pieces per component, where G6 only measured up to 8,000) and to confirm the
three real disagreeing caps at `(633,-97,-61)`/`(633,-97,-62)` get repaired by
`fuse_disagreeing_pairs` rather than declined — see the item above, which was
blocked on this one landing first.

Note that the obvious-looking alternative — welding boundary pieces to each other
by index, the way the interior is joined — **does not work**, and the reason is
worth keeping. `BRepTools_ReShape` will swap an edge inside a face, but replacing
that edge's *vertices* leaves the neighbouring edges pointing at the old ones and
the wire comes apart (`BRepCheck_NotConnected`, volume wrong, every edge still
used exactly twice). Keeping the vertices makes the same swap exact. Two boolean
pieces cannot each keep their own vertices, so one side must be geometry the
program builds itself — which is true of the interior and false of the boundary.


