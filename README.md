# latticegen2

Parametric generation of PCM heat exchanger lattice geometry. `latticegen2` fills the
solid volume of an input STEP part with a diamond-strut lattice, trimmed exactly to the
input geometry's boundary, and writes the result as a single watertight AP214 STEP file
(possibly multiple bodies if the input geometry cuts the lattice into disconnected
islands). See [docs/specification.md](docs/specification.md) for the full requirements
this tool implements, and [docs/algorithm.md](docs/algorithm.md) for the exact algorithm
and optimization strategy.

## Contents

- [Dependencies](#dependencies)
- [Installation](#installation)
- [Offline deployment](#offline-deployment)
- [Usage](#usage)
- [Parameter reference](#parameter-reference)
- [Input format](#input-format)
- [Output format](#output-format)
- [Example workflows](#example-workflows)
- [Memory usage](#memory-usage)
- [Algorithm overview](#algorithm-overview)
- [Testing](#testing)

## Dependencies

| Dependency | Role | License |
|---|---|---|
| [Julia](https://julialang.org/) 1.10 | Language runtime | MIT |
| [Gmsh.jl](https://github.com/JuliaFEM/Gmsh.jl) 0.3 | Only third-party Julia dependency: binds the Gmsh SDK (which bundles the OCCT B-rep kernel) for STEP/BREP I/O, meshing, and boolean geometry operations | MIT (binding); GPL-2.0-or-later w/ exception (gmsh SDK); LGPL-2.1 (OCCT) |

Everything else used by `src/` is Julia standard library (`Dates`, `Distributed`,
`LinearAlgebra`, `Printf`, `Test`) — no other third-party packages. License texts for
all of the above are in [licenses/](licenses/), cross-referenced in
[licenses/libraries.md](licenses/libraries.md).

## Installation

Requires Julia 1.10, managed via [juliaup](https://github.com/JuliaLang/juliaup):

**Windows:**
```bash
winget install --id Julialang.Juliaup -e
juliaup add 1.10
juliaup default 1.10
```

**Linux:**
```bash
curl -fsSL https://install.julialang.org | sh
juliaup add 1.10
juliaup default 1.10
```

Then, from the repository root, resolve and precompile the project's one dependency
(this is the **only step in the whole workflow that needs network access** —
`specification.md` §2's offline requirement applies to every subsequent run):

```bash
julia --project=. -e "using Pkg; Pkg.instantiate(); Pkg.precompile()"
```

Verify the install:
```bash
julia --project=. -e "using Pkg; Pkg.test()"
```

## Offline deployment

`latticegen2` must run with zero network access once installed (specification.md §2).
`Manifest.toml` is committed to this repository specifically so that `Pkg.instantiate()`
always resolves to the exact same dependency versions — this is what makes the following
offline transfer procedure reproducible:

1. On a machine with network access, run the Installation steps above (this populates
   `~/.julia` — or `%USERPROFILE%\.julia` on Windows — with every package/artifact
   `latticegen2` needs, including the Gmsh SDK binaries).
2. Copy the entire Julia install directory (from `juliaup`) and the entire `~/.julia`
   depot directory to the offline machine, along with this repository.
3. On the offline machine, either restore `juliaup`'s directory layout, or point Julia
   at the copied depot via the `JULIA_DEPOT_PATH` environment variable.
4. Run `julia --project=. -e "using Pkg; Pkg.instantiate()"` once more **offline** — with
   the depot already populated, this only verifies/links the local cache and requires no
   network access.

## Building a standalone executable

As an alternative to the depot-copy procedure above, [PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl)
can produce a self-contained Windows distribution — a `latticegen2.exe` with its own
bundled Julia runtime, sysimage, and artifacts (the same `gmsh_jll`/`OCCT_jll` binaries
the normal invocation uses) — that runs on a machine with no Julia install at all
(specification.md §2). This is a **build-time-only** tool: `PackageCompiler` is declared
in its own environment, [tools/build/Project.toml](tools/build/Project.toml), never in
the root `Project.toml`, so it is never resolved by an ordinary
`julia --project=. src/main.jl` run.

```bash
# Once, on a machine with network access (fetches PackageCompiler.jl itself):
julia --project=tools/build -e "using Pkg; Pkg.instantiate()"

# Build (takes on the order of 15-45+ minutes; output is several hundred MB):
julia --project=tools/build tools/build_app.jl
```

This produces `build/latticegen2-app/bin/latticegen2.exe`, which takes the same
arguments as `latticegen2.bat`:

```bash
build\latticegen2-app\bin\latticegen2.exe -i <input.step> -cc <mm> -t <mm> [options]
```

`build/` is gitignored — it's a local build artifact, not something checked in. The
entry point PackageCompiler calls, [`julia_main()`](src/app_entry.jl), just `include`s
`src/main.jl`, so the built executable and the `julia --project=. src/main.jl`/wrapper-
script invocation run identical code (see `src/app_entry.jl`'s docstring for why it's
structured that way rather than as an ordinary function).

## Usage

```bash
julia --project=. src/main.jl -i <input.step> -cc <mm> -t <mm> [options]
```

or, using the convenience wrapper (equivalent, same arguments):

```bash
./latticegen2.sh -i <input.step> -cc <mm> -t <mm> [options]      # Linux
latticegen2.bat -i <input.step> -cc <mm> -t <mm> [options]       # Windows
```

Exactly one of the following two parameter pairs is required, in addition to `-i`,
`-cc`, and `-t`:

- `--cores <n> --ram <GB>` — auto-tune the number of parallel worker processes and tile
  size from the given hardware budget (docs/algorithm.md §7.1).
- `--workers <n> --tile-cells <n>` — use these optimization parameters directly, no
  auto-tuning/calibration probe.

Run `julia --project=. src/main.jl --help` for the full flag summary.

The wrapper scripts pass `-t auto` to Julia, enabling Julia's own thread pool (distinct
from the `--cores`/`--workers` `Distributed` *process* pool above) for the classification
stage's threaded per-strut loop. Invoking `src/main.jl` directly without `-t auto` (or
`JULIA_NUM_THREADS` set) still works, just single-threaded for that stage.

### Cancelling a run

Press **Ctrl+C** to stop a run in progress. The tool shuts down in order rather than
dying on the spot: in-flight tiles are allowed to unwind, the worker processes are
stopped (force-stopped only if one is stuck inside a long geometry-kernel operation
and does not exit within ~2 seconds), a single `CANCELLED: ...` line goes to the
console and the `.log` file, and the exit code is `130`. The run's `temp/<timestamp>/`
directory is kept, so completed tiles remain available for inspection. See
docs/algorithm.md §9.1.

## Parameter reference

| Flag | Type | Required | Units | Range | Default | Description |
|------|------|----------|-------|-------|---------|-------------|
| `-i`, `--input` | path | **yes** | — | — | — | Input STEP file defining the lattice bounds |
| `-o`, `--output` | path | no | — | — | `<input_stem>-cc<cc>t<t>.step`, next to the input | Output `.step` path |
| `-cc` | float | **yes** | mm | 0.4–50 | 5 (see note below) | XY-plane distance between the bottom nodes of two adjacent cells (cube edge `a = cc/√2`) |
| `-t` | float | **yes** | mm | 0.4–20, and `t < cc/√2` | 1 (see note below) | Side length of the diamond strut profile |
| `--cores` | int | one of these two pairs required | count | 1–128 | — | Physical cores available, for auto-tuning |
| `--ram` | float | " | GB | 1–1024 | — | RAM budget available, for auto-tuning |
| `--workers` | int | " | count | 1–128 | — | Explicit parallel worker process count |
| `--tile-cells` | int | " | count | 2–64 | — | Explicit tile edge length, in lattice cells. **Values above 8 are counter-productive**: per-tile fuse cost grows well past quadratic with strut count (~N^2.5, measured), so a tile past that point can take *much* longer to fuse than the memory savings are worth. Auto-tuning (`--cores`/`--ram`) never picks above 8 for this reason (docs/algorithm.md §7.1); prefer it over a large explicit `--tile-cells` unless you have a specific reason to pick your own. |
| `-bg`, `--background` | flag | no | — | — | disabled | Run at below-normal OS scheduling priority |
| `-v`, `--verbose` | flag | no | — | — | disabled | Verbose console output (the `.log` file is always written in full regardless) |
| `-h`, `--help` | flag | no | — | — | — | Show usage and exit |


Every parameter is validated — including the `t < cc/√2` cross-constraint (a strut
thicker than the cube edge cannot fit in one lattice cell) — **before any computation
starts** (specification.md §7). Invalid parameters exit with code `2` and a
human-readable reason.

## Input format

Any STEP file (AP203/AP214/AP242 — all import successfully via the OCCT kernel)
containing at least one solid volume, in millimeters, in whatever coordinate system the
output should be placed in (the lattice is generated directly in the input's own
coordinate system — no re-centering or scaling). Both files under
[test/](test/) are real example inputs (an 80 mm test sphere and a larger real
heat-exchanger cavity volume).

## Output format

- **Schema:** STEP AP214 (`AUTOMOTIVE_DESIGN`), exact B-rep solid geometry, millimeters
  (specification.md §5 — the spec's original `AP203` requirement was changed to AP214 by
  user decision, since the OCCT kernel used here writes AP214 natively).
- **Bodies:** normally one solid; the file may contain **multiple** disconnected solid
  bodies if the input geometry's boundary cuts some struts into islands disconnected
  from the main lattice (specification.md §1) — this is expected, not an error. No
  floating body smaller than `t³` mm³ is ever emitted (specification.md §5).
- **No outer shell:** the lattice is not wrapped in a bounding solid shell — struts are
  trimmed flush against the input boundary and are meant to be merged with an
  enveloping part on import (specification.md §4.3).
- **Metadata:** the STEP model/part name is `<input_stem>+cc<cc>+t<t>`, and the file's
  `FILE_DESCRIPTION` header entry additionally records the full generation parameters
  (input path, `cc`, `t`, generation timestamp) — see docs/algorithm.md §8.
- **Log file:** every run also writes `<output_stem>.log` (never `<output>.step.log`),
  containing the full run header, per-stage timings, and end-of-run summary
  (specification.md §3) — always written in full regardless of `-v`.

## Example workflows

Explicit worker/tile parameters (no auto-tuning, fully reproducible regardless of host):
```bash
julia --project=. src/main.jl -i test/80mm-test-ball.step -cc 10 -t 2 \
    --workers 4 --tile-cells 6 -v
```

Auto-tuned from a 6-core / 32 GB workstation, running at below-normal priority so it
doesn't disturb interactive desktop use:
```bash
julia --project=. src/main.jl -i test/TD_HX_Indre_Volum.step -cc 5 -t 1 \
    --cores 6 --ram 24 -bg -v
```

Custom output path:
```bash
julia --project=. src/main.jl -i part.step -cc 8 -t 1.5 -o out/part_lattice.step \
    --workers 8 --tile-cells 8
```

## Memory usage

Priority #2 (after correctness): the tool must never destabilize the host by
over-committing RAM (specification.md "Key Considerations"). docs/algorithm.md §7
describes the full model; in summary:

- With `--cores`/`--ram`, a small calibration probe (an isolated, synthetic 4×4×4-cell
  reference block — independent of the actual input geometry) measures both
  memory-per-strut and fuse time before any real work starts. Tile size is bounded by
  **both**: memory (`workers × struts_per_tile × mem_per_strut ≤ 0.6 × RAM budget` —
  the 0.6 factor reserves headroom for the master process, OS, and the geometry
  kernel's own working memory beyond steady-state) *and* fuse time (extrapolated from
  the probe, conservatively assuming worst-case quadratic growth), then hard-capped at
  8 cells regardless — per-tile fuse cost was measured growing well past quadratic
  (~N^2.5) with strut count, so memory headroom alone is not a safe sizing signal past
  that point (docs/algorithm.md §7.1, §11.2).
- A runtime watchdog pauses dispatch of new tiles (without killing in-flight work) if
  observed memory trends toward 0.8× the RAM budget, bounded to at most 120 s of
  continuous pause per dispatching task before resuming anyway (memory usage is a
  monotonic high-water mark within a run, so an unbounded wait here would be a
  guaranteed hang, not a slow pause).
- Large intermediate results — the imported input body, each tile's fused geometry, and
  each distributed assembly merge round's output — are staged to disk as `.brep` files
  under `temp/<yyyymmdd-HHMMSS>/`, next to the output file (specification.md §4.4) —
  kept for post-mortem analysis on failure, deleted automatically on success. The final
  hand-off from assembly to export is the one exception: it stays in one gmsh session
  rather than round-tripping the completed lattice through disk a last time.
- With `--workers`/`--tile-cells` given explicitly instead, no calibration runs; memory
  behavior is directly determined by the tile size you choose — smaller tiles use less
  memory per worker at some cost to speed (fewer struts fused per boolean operation
  amortizes worse — see the optimization table below).

## Algorithm overview

Full detail, exact math, and mermaid diagrams: **[docs/algorithm.md](docs/algorithm.md)**.
In brief, the pipeline is:

```
parse args -> import STEP -> tessellate surface -> classify every candidate
strut (interior / boundary / outside, via mesh-distance + ray-cast) -> tile
struts into n×n×n blocks -> process boundary tiles + one reference interior
tile in parallel worker processes (trimming boundary struts via the
operand-disjointness-safe trim_disjoint, dropping tile-local floating
islands) -> reproduce every other full-interior tile via cheap
copy+translate (periodicity shortcut) -> distributed hierarchical assembly
(merge rounds across worker processes, final fuse on the master) ->
floating-body-only cleanup (filter_floating!) -> export STEP + rewrite
header -> round-trip verify
```

Optimization strategy (docs/algorithm.md §11/§11.2 has the full complexity analysis
and investigation history):

| Lever | Effect |
|---|---|
| Curvature-adaptive surface tessellation | Avoids over-refining large gently-curved regions (~100x fewer triangles than a uniform fine target on the 80mm test ball) |
| Classify-before-boolean | Expensive OCCT booleans only run on the O(surface-area) boundary struts, not the whole O(volume) candidate set |
| Spatial-hash-accelerated ray casting (3D DDA) | Point-in-solid tests only visit mesh triangles in cells the ray actually crosses, not the whole mesh |
| Threaded classification (`Threads.@threads`, needs `-t auto`/`JULIA_NUM_THREADS` — the wrapper scripts set this) | Parallel per-strut classification within each process, layered under the process-level parallelism below |
| `trim_disjoint` / COMMON operand-disjointness invariant | Prevents boundary-strut trimming from fragmenting mutually-overlapping struts into extra pieces — a correctness fix, not just a speed one (docs/algorithm.md §6.3, §11.2) |
| Multi-operand fuse per tile | One n-way OCCT GeneralFuse beats incremental pairwise fusing |
| Unit-tile copy+translate (periodicity shortcut) | The dominant interior-volume fuse cost collapses to one fuse total, plus cheap copies |
| Worker-side floating-island removal | Drops provably-floating sub-threshold solids at the tile that produced them, before they ever reach assembly or export |
| Process-parallel tiles *and* assembly (`Distributed`) | Scales across cores despite the geometry kernel not being thread-safe; assembly merges tiles in rounds across the same worker pool instead of one giant single-session master fuse, bounding peak master memory to the last round's survivors |
| Longest-processing-time-first tile dispatch | Large tiles start first, so the tail of the tile stage isn't dominated by a handful of big tiles starting last while other workers sit idle |
| Balanced hierarchical fuse, spatial-locality batched, AABB pre-filtered (one shared sync per round) | Bounded per-call operand complexity; batches of spatially-adjacent solids; groups with no bounding-box overlap (the only unconditionally-safe fallback case) skip the OCCT call entirely — with the box-lookup cost paid once per round, not once per solid (a per-solid sync was measured 202× slower and made this filter a net loss for a long time, docs/algorithm.md §11.2); a generous wall-clock circuit breaker guards against a runaway fuse without cutting off legitimate work, since doing so risks leaving self-intersecting geometry un-fused |
| `filter_floating!`: resolve-then-classify cleanup | Only ever deletes solids proven disconnected (never merely "small"); an unresolvable connected sub-threshold fragment fails the run rather than being silently deleted or kept (docs/algorithm.md §8) |
| Calibration probe (memory *and* time) + dual-bounded tile sizing + RSS watchdog | Memory stability on the target workstation, and tile sizes that stay below the measured fuse-time performance knee |

## Testing

See [docs/testing.md](docs/testing.md) for the full procedure. Quick reference:

```bash
julia --project=. -e "using Pkg; Pkg.test()"    # unit tests
julia --project=. tools/e2e.jl                  # end-to-end smoke-fast scenario
```
