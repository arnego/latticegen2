# Testing Reference Guide

This file contains the required verification and testing procedures for this project.
All test files and assets shall reside within the test/ folder.

## Unit testing:

Run unit tests for parameter validation and intermediate calculations

- Run all tests:
  ```
  julia --project=. -e "using Pkg; Pkg.test()"
  ```
  (equivalently: `julia --project=. test/runtests.jl`)
- Run only the fast, gmsh-independent subset (pure math: cli.jl, runlog.jl, lattice.jl,
  tiling.jl, stepmeta.jl) during iterative development:
  ```
  julia --project=. -e "using Test; include(\"src/latticegen2.jl\"); using .latticegen2; @testset \"pure-math\" begin; include(\"test/test_runlog.jl\"); include(\"test/test_cli.jl\"); include(\"test/test_lattice.jl\"); include(\"test/test_tiling.jl\"); include(\"test/test_stepmeta.jl\"); end"
  ```
- Run a single test file directly (each is a self-contained `@testset`, included by `test/runtests.jl`):
  ```
  julia --project=. -e "using Test; include(\"src/latticegen2.jl\"); using .latticegen2; include(\"test/test_lattice.jl\")"
  ```
  Substitute the file for any of: `test_runlog.jl`, `test_cli.jl`, `test_lattice.jl`,
  `test_tiling.jl`, `test_stepmeta.jl` (fast, no gmsh), or `test_geomkernel.jl`,
  `test_classify.jl`, `test_pipeline.jl`, `test_cleanup.jl` (slower, require the
  gmsh/OCCT kernel). `test_cleanup.jl` covers `filter_floating!`'s floating-body-only
  cleanup gate specifically (docs/algorithm.md §8) — kept separate from
  `test_pipeline.jl`'s broader orchestration tests since it's a priority-#1-critical
  correctness gate on its own.

## E2E verification:

1. Run the script with the project's verification input geometry and parameters. The
   run's log file (`<output-stem>.log`, specification.md §3) already contains the
   required run data — input parameters, start date/time, duration, run characteristics
   (tile counts, worker counts, per-stage timings), and maximum memory usage — save this
   file for analysis and attach it to the pull request.

   ```
   julia --project=. src/main.jl -i test/80mm-test-ball.step -cc 20 -t 4 -o /tmp/smoke-fast.step --workers 4 --tile-cells 6 -v
   ```

   Or run the full automated harness, which invokes the scenario above as a subprocess
   and performs every check in §6.2 below automatically:

   ```
   julia --project=. tools/e2e.jl
   ```

2. Verification: `tools/e2e.jl` also runs `smoke-verified` (specification.md §6.1),
   which compares its output against the committed golden sample
   `test/80mm-test-ball-cc20t4-golden-sample.step` via `golden_sample_volume_diff`.
   `dense-lattice` is implemented the same way but self-skips: its golden sample
   (`test/test-cylinder-cc10t1.5.step`) does not exist yet. Originally specified at
   -cc 5 -t 1, the one attempt to generate that denser golden sample ran for hours
   before being manually terminated — root-caused and fixed (not a crash: an
   auto-tuned tile size past the fuse-time performance knee, plus an unconditional
   sub-threshold cleanup rule that was deleting connected junction material; see
   docs/algorithm.md §11.2 for the full investigation). The scenario's params were
   changed 2026-08-09 to -cc 10 -t 1.5 (a less dense lattice capable of finishing
   within a reasonable time) with a 60-minute runtime budget. Regenerating the
   golden sample and committing it is a follow-up step, tracked as an open item in
   specification.md §9, not automated here. Until that sample exists, `tools/e2e.jl`
   provides the `smoke-fast`/`smoke-verified` output `.step` files (and console/log
   summaries) for manual user verification, and automatically ensures the generated
   geometry is manifold, non-self-intersecting, and that sub-threshold-solid removals
   stayed a small fraction of total solids (a regression check on its own — this is
   exactly the signal that would have caught the `dense-lattice` blow-up before its
   multi-hour cleanup tail, docs/algorithm.md §11.2):

   ```
   julia --project=. tools/e2e.jl
   ```

   When a golden sample *is* configured for a future scenario, similarity is checked by
   subtracting the candidate and golden geometries both ways (near-zero remainder means
   equivalent volume) via `tools/verify_geometry.jl`'s `golden_sample_volume_diff`:

   ```julia
   include("tools/verify_geometry.jl")   # after include("src/latticegen2.jl"); using .latticegen2
   golden_sample_volume_diff("path/to/candidate.step", "path/to/golden.step")
   ```

## Goal oriented performance optimization:

The log file written by every run (`<output-stem>.log`) records, per docs/algorithm.md
§9: a full run header (all input parameters, start timestamp), one line per pipeline
stage (import, tessellate, classify, tiling, tile_stage, assembly, export, verify) with
its wall-clock duration, per-tile stats as each tile completes — strut counts, total
elapsed time, peak RSS, **and** the per-stage breakdown within the tile
(`t_interior`/`t_boundary`/`t_final`, the three `balanced_fuse!` calls a tile makes)
and its `dropped_islands` count (docs/algorithm.md §6.4a) — one line per distributed
assembly merge-round group as it completes (docs/algorithm.md §6.5), calibration-probe
results when `--cores`/`--ram` auto-tuning ran (now `mem_per_strut`, the probe's
`(struts, elapsed)` pair, and which of `n_mem`/`n_time` bound the chosen tile size —
docs/algorithm.md §7.1), and the mandatory end-of-run summary. Any warning a tile's
`balanced_fuse!`/boundary-trim step produced is also logged, prefixed with the tile's
key — previously invisible whenever a tile happened to run on a worker process, since
workers have no log file of their own to write to. To iterate on performance:

1. Run with `-v` for full console visibility, or just inspect the `.log` file afterward
   (it is always written in full regardless of `-v` — specification.md §3).
2. Compare the per-stage timings across runs/parameter sets to identify the current
   bottleneck stage.
3. Cross-reference docs/algorithm.md §11 (Complexity analysis and optimization strategy
   summary) for the specific lever expected to affect that stage.
4. Re-run and compare the same stage's timing before/after a change.

## Verification Checklist

1. Ensure all linting passes (`julia --project=. -e "using Pkg; Pkg.test()"` — Julia has
   no separate lint step for this project; test failures are the linting signal).
2. Verify edge cases for any modified boundary logic — in particular:
   - Parameter range boundaries (`-cc`/`-t` at 0.4/50/20, the `t < cc/√2` cross-constraint).
   - Classification margin edge cases (a strut exactly at the `r + d` boundary).
   - Tile-partition edges (negative lattice indices, `is_full_interior`'s fringe-tile guard).
   - AABB-overlap-graph edge cases (`overlap_components`): a chain that must merge into
     one component vs. two genuinely separate clusters (docs/algorithm.md §6.3, §8).
   - `filter_floating!`'s three-way outcome (removed / kept / hard-fail) at each
     boundary: singleton-and-small, connected-and-resolves, connected-and-unresolvable
     (docs/algorithm.md §8 — `test/test_cleanup.jl`).
   - `derive_tile_size`'s dual memory/time bound and its `n <= 8` hard cap
     (docs/algorithm.md §7.1).
3. Run `tools/e2e.jl` and confirm all checks pass before opening a pull request.
4. Run `/code-review` (per CLAUDE.md) and address any findings before pushing.
