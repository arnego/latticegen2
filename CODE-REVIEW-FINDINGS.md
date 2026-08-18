# Code review: v2.1.0 → main (01b4206)

Findings from `/code-review ultra`, and their disposition.

**Review target:** a throwaway branch `review/v2.1.0-to-main-code` — a squashed
tree of `main` parented on `v2.1.0` (`29f240b`), so the branch-vs-merge-base
diff was exactly `v2.1.0..main`, minus two measurement-record paths held back to
fit ultrareview's 8,000-line budget (the full 12,555-line diff was refused):

- `tools/prototypes/` — one-off gate scripts, `export-ignore`d, never shipped.
- `docs/profiling-reports.md` — an append-only log of profile runs.

Everything shipped is included in full, at its exact `v2.1.0..main` line count:
`src/` 2,546, `test/` 2,053, `tools/` (e2e, verify_geometry, build_release,
smoke_bundle, profile_*) 689, `docs/` 2,404 — both normative documents,
`algorithm.md` and `specification.md`, complete — and packaging 77. **7,769 lines
over 43 files.** A sibling branch `review/v2.1.0-to-main` carried the untrimmed
12,555-line target. Both were review scaffolding — they held `main`'s own tree
and nothing unique — and both were deleted once the review returned.

**Fixes land on:** `claude/code-review-v2-1-0-main-3ed6b6`, branched from `main`,
which merges back normally.

**Scope:** 14 commits (#10–#23).

| Area | Files |
|---|---|
| Source | `src/latticegen2/` (12 — `occ.py`, `pipeline.py`, `weld.py`, `boundary.py`, `classify.py`, `interior.py`, `parallel.py`, `cli.py`, `sysinfo.py`, `stepout.py`, …) |
| Tests | `test/` (13) |
| Harnesses | `tools/` (6) |
| Docs | `docs/` (4), `README.md`, `CLAUDE.md`, licenses, packaging |

**Not covered by this pass:** the 15 prototype scripts and the profiling log. They
are measurement records rather than shipped code, but if a finding turns on a gate
whose evidence lives in `tools/prototypes/RESULTS.md`, read it directly — it is in
the tree, just not in the diff.

---

## Findings

Three reports, two distinct findings — the reviewer filed the `classify` one
twice (`bug_003`, `bug_003_1`) with identical substance. Both were verified
against the code before anything was changed; both are real.

| # | Severity | Site | Finding | Disposition |
|---|---|---|---|---|
| 1 | normal | [weld.py:434](src/latticegen2/weld.py#L434) | `_sew_all_tiles` crashes under `--cores 1` on any component large enough to tile | **fixed** |
| 2 | nit | [classify.py:944](src/latticegen2/classify.py#L944) | `_worker_classify` rebuilt `_ClassifyIndex` per *slice* while its docstring and algorithm.md §5.4 both said per *worker* | **fixed** |

### 1 — `_sew_all_tiles` crashes under `--cores 1`

Confirmed exactly as reported. `pipeline._run` builds `WorkerPool(args.workers)`
unconditionally, and `WorkerPool(1)` is inert by design — `__enter__` creates no
`mp.Pool`, so `.active` is `False` and `.run()` raises. `_sew_all_tiles`'s guard
tested `pool is None and workers <= 1`, which is `False` for a pool that exists
and cannot run, so the call fell through to the transient-pool branch, built a
*second* inert pool, and raised `ProcessingError` — **exit 4, no output**.

`--cores` is 1–128 per specification.md §3, so this is refusing input the CLI
accepts, which docs/algorithm.md §11 rules out as a failure mode.

**Why nothing caught it.** It needs a component past `MIN_PIECES_TO_TILE`
(1,500). No committed scenario tiles at all — `dense-lattice` is well under —
and every tiling test in `test_weld.py` passes `pool=None`, which takes the
correct branch. The rehearsal part has 21,955 pieces in one component and would
have hit it immediately.

**Fix:** the guard now reads `(pool is None or not pool.active) and workers <= 1`.
I audited every pool-dispatch site in `src/` rather than trusting the report's
claim that this was the only one: `classify.py:1026`, `pipeline.py:738` and
`weld.py:619` all already gate on `pool.active`, and `boundary.py:757` is
unreachable with an inert pool because `boundary.py:701` returns sequentially at
`workers <= 1` first. This was the only hole.

I took the guard fix over the reviewer's suggested alternative — a third
`elif workers <= 1` branch mirroring `_sew_round_two` — because the early guard
also skips writing the per-tile `.brep` files the sequential path never reads.
Same behaviour, less wasted I/O, no duplicated loop.

### 2 — `_ClassifyIndex` rebuilt per slice, not per worker

Confirmed. `_classify_parallel` dispatches `workers * 4` slices, `WorkerPool.run`
uses `imap(chunksize=1)`, and `_worker_classify` called `classify_slice` with no
`index`, so `classify_slice` built a fresh `_ClassifyIndex` — plus a fresh
`load_mesh` — on every call. Roughly four rebuilds per worker where the
docstring and algorithm.md §5.4 both priced one.

**Fixed the code rather than the doc.** CLAUDE.md makes algorithm.md normative
("source code must match it exactly"), and §5.4's claim is the useful one — so
the worker now memoises `(mesh, _ClassifyIndex)` at module scope keyed on
`(mesh_path, cc, t)`. Reuse is sound because `_ClassifyIndex` is read-only once
constructed and a function of `lp` and the mesh alone, which is the same
property that makes the sweep divisible.

Wall clock only — the classification is bit-identical either way, which the
existing serial/parallel identity test already pins.

**One thing the review did not raise, checked because memory is priority #2.**
Caching means each worker now *retains* the mesh and its indices after
`classify` ends, for the rest of the run. Measured rather than assumed: a
`_ClassifyIndex` over 81,920 triangles retains ~19.5 MB, so the rehearsal's
28,654-triangle mesh is roughly 5–10 MB per worker, ~30–60 MB across six
against a ~19 GB peak. Peak *during* `classify` is unchanged — the index was
already live for the duration of every job. Not worth an eviction mechanism. The reviewer's ~1.1 s
estimate at rehearsal scale is its own arithmetic, not a measurement, and I have
not re-measured it; the fix is justified by the doc/code mismatch, not by a
projected saving.

### Verification

Both fixes carry a regression test, and each was **confirmed to fail on the
pre-fix code** — reverted, run, restored — rather than merely passing after:

- `test_weld.py::test_tiled_sew_falls_back_to_sequential_when_the_shared_pool_is_inert`
  raises the exact reported `ProcessingError` on the old guard.
- `test_classify.py::test_worker_reuses_one_classify_index_across_the_slices_it_is_given`
  fails `0 == 1` on the old worker body.

Gate, both green after the change:

- `python -m pytest test -q` — **244 passed, 1 skipped**.
- `python tools/e2e.py` — **all four scenarios passed**: `smoke-fast`,
  `smoke-verified`, `dense-lattice` (37.5 s), `invalid-input`; both golden
  samples at **0 mm³**, all solids `BRepCheck_Analyzer`-valid, 0 mm³ outside the
  input body, 0 crossing pairs.

Neither fix changes what a run produces — the guard only picks a different route
to the same sewn shell, and the cache only avoids rebuilding a read-only index —
so an unchanged golden-sample result is the expected outcome here, not evidence
that the fixes are exercised. What shows they are exercised is that both
regression tests fail on the pre-fix code.

---

## Disposition key

- **fixed** — change landed on the fix branch, with the commit named.
- **no change needed** — verified as a non-issue; record *why*, per the project's
  habit of keeping the reasoning rather than the conclusion.
- **deferred** — real but out of scope; add to `docs/specification.md` §10 with
  enough context to act on without re-deriving the diagnosis.

Anything touching the pipeline is gated by `python -m pytest test -q` and
`python tools/e2e.py` (both golden samples at 0 mm³) before it is called fixed —
see `docs/testing.md`.
