# Testing Reference Guide

Required verification and testing procedures for this project. Test code and
assets live in `test/`; the harnesses that drive whole runs live in `tools/`.

All commands assume a Python 3.11+ interpreter with the two dependencies from
[../README.md](../README.md) installed. Set `LATTICEGEN2_PYTHON` if the
interpreter you want is not the default `python`.

## Unit testing

Parameter validation, lattice mathematics, the junction template, the interior
shell build, classification, the mesh gates, the connectivity rule, and the STEP
header rewrite.

- Run everything:
  ```
  python -m pytest test -q
  ```
- Run the fast, kernel-free subset (pure maths — CLI, lattice, connectivity,
  STEP header) during iterative development:
  ```
  python -m pytest test/test_cli.py test/test_lattice.py test/test_connect.py test/test_stepmeta.py -q
  ```
- Run one file:
  ```
  python -m pytest test/test_junction.py -q
  ```

| File | Covers | Needs the geometry kernel |
|---|---|---|
| `test_lattice.py` | docs/algorithm.md §2.3's verified identities, including the mutual orthogonality the architecture rests on; index ranges; naming | no |
| `test_cli.py` | Flags, ranges, the `t < a` cross-constraint, path resolution, preflight checks | no |
| `test_connect.py` | The junction graph and the floating-body rule's three outcomes | no |
| `test_stepmeta.py` | Quote-aware STEP header editing, including never overwriting a populated `FILE_SCHEMA` | no |
| `test_junction.py` | Cap integrity across the parameter range, the inradius argument behind it, and the exact `N x volume(J)` identity for instanced grids | yes |
| `test_classify.py` | Distance primitives, spatial indices, ray parity, node classes, and both mesh gates — including the pole-degeneracy regression from issue #6 | yes |
| `test_main.py` | Exit codes and the "exactly one reason line" rule, before and after the log file opens | yes |
| `test_pipeline.py` | Same-domain unification's fallback ladder — a kernel that refuses to merge must yield a larger file, never a failed run | yes |

## E2E verification

```
python tools/e2e.py
```

Runs every scenario in [specification.md](specification.md) §6.1 as a subprocess
and applies every applicable §6.2 check. Scenarios are independent — one failing
scenario does not prevent the others from running — and the exit code reflects the
run as a whole. `--only smoke-fast,invalid-input` restricts the set.

Each generated run writes its own `<output-stem>.log` containing the required run
data (parameters, start time, duration, run characteristics, peak memory). Save
it for analysis and attach it to the pull request.

To drive a single scenario by hand:

```
python src/main.py -i test/80mm-test-ball.step -cc 20 -t 4 -o /tmp/smoke.step -v
```

### What the harness checks

| Check | How |
|---|---|
| Exit code and one human-readable reason line | subprocess return code and stderr |
| STEP written and non-empty | filesystem |
| Parses back successfully | `STEPControl_Reader` round trip |
| **Exact B-rep validity** | `BRepCheck_Analyzer` on every solid |
| Closed manifold | every mesh edge used by exactly 2 triangles |
| No self-intersections | `triangles_properly_cross` — plane-straddle pre-check, then edge piercing |
| No material outside the input body | boolean cut of output against input, volume ≈ 0 |
| Bounding box within input + (cc+t) | direct comparison |
| Runtime budget | wall clock: 10 min for `smoke-fast` and `dense-lattice`, 20 min for `smoke-verified` |
| Golden-sample match | symmetric-difference volume both ways, tolerance `t³` |

The self-intersection check's plane-straddle pre-check is load-bearing, not
decoration: without it, two separate solids merely *touching* along a coincident
face report hundreds of false crossings purely because their independent
triangulations are not vertex-aligned. Measured on two boxes sharing one exact
face with zero volume overlap: 344 false positives without the pre-check, 0 with
it.

### Golden samples

`test/80mm-test-ball-cc20t4-golden-sample.step` and
`test/test-cylinder-cc10t1.5-golden-sample.step` are reviewed baselines,
committed after manual verification of the geometry they contain. Both are
reproduced by the current implementation with a symmetric-difference volume of
0 mm³.

Regenerating one is a review decision, not a routine step: a golden sample only
has value while somebody has actually looked at it. If you do replace one,
confirm the volume-based comparison against the *existing* sample first — that
is what distinguishes a change in how the geometry is described from a change in
the geometry itself.

Both are compared with `golden_sample_volume_diff`, which cuts candidate and
golden against each other both ways and takes the larger remainder. That
comparison is a general boolean between two complete lattices, so it is slow at
lattice scale: seconds for the ball, **306 s** for the cylinder. The harness
gives it a generous budget (`GOLDEN_EXACT_BUDGET_S`) and only falls back to a
sampled equivalence check if it is exceeded — the fallback is labelled as
weaker wherever it is reported, and never turns an unmeasured result into a
pass.

A mismatch is a stop-the-line result: **never adjust a tolerance to make it
pass**, and do not assume the generator is at fault either. Investigate, and
escalate to the user for manual verification.

To compare two files directly:

```python
import sys; sys.path.insert(0, "tools")
import verify_geometry as vg
vg.golden_sample_volume_diff("candidate.step", "test/80mm-test-ball-cc20t4-golden-sample.step")
```

## Performance work

Every run writes `<output-stem>.log` with a header, one line per pipeline stage
with its wall-clock duration, template and mesh statistics, classification counts,
boundary-trim progress, the aggregate floating-body line, and the end-of-run
summary including peak memory. To iterate:

1. Run with `-v`, or just read the `.log` afterwards — it is always written in
   full regardless of `-v`.
2. Compare per-stage timings across runs to find the current bottleneck.
3. Cross-reference [algorithm.md](algorithm.md) §12 for the lever that affects
   that stage.
4. Re-run and compare the same stage.

For reference, the two committed scenarios on a 6-core / 32 GB workstation:

| Scenario | Total | Dominant stages |
|---|---|---|
| 80 mm ball, `cc=20 t=4` | ~7 s | boundary trim, export |
| test cylinder, `cc=10 t=1.5` | ~61 s | sewing, classify, verify, simplify |

Note that the `simplify` stage (same-domain unification, docs/algorithm.md §9)
has a *negative* net cost: it takes ~8 s but halves the face count, which takes
more than that back out of export and the round-trip check. Removing it would
make the run slower as well as the output twice as large.

The measurements behind the architecture's load-bearing choices (junction cap
integrity, join-mechanism throughput, per-junction intersection latency, STEP
writer throughput) are in
[`../tools/prototypes/RESULTS.md`](../tools/prototypes/RESULTS.md), with the
scripts that produced them alongside. Re-run one with, e.g.:

```
python tools/prototypes/g2_instancing_join.py
```

## Release bundles

The offline bundles are built by `tools/build_release.py` and verified by
`tools/smoke_bundle.py`. Both run locally, not only in CI:

```
python tools/build_release.py --flavour both
```

```
python tools/smoke_bundle.py --platform win64
```

The smoke script is the gate that stands between a green build and a published
release, and it deliberately tests the **archive** rather than the staging tree.
Extraction is where the executable bit is lost, where a CRLF shebang bites, and
where a relocated interpreter stops finding its own libraries — none of which
appear in a build log. It extracts outside the repository, scrubs
`LATTICEGEN2_PYTHON`/`PYTHONPATH`/`PYTHONHOME` from the environment so a bundle
cannot borrow the development machine's interpreter, runs one generation through
the bundle's own launcher, and checks the result with `tools/verify_geometry.py`
(exact B-rep validity, no material outside the input, closed manifold).

Note that the builder packages a **git ref**, not your working tree, so
uncommitted changes will not appear in a bundle. It must also run under the
Python version the bundle targets when downloading wheels, since wheels are
version-specific; it says so plainly if you get that wrong.

The full release procedure is [release.md](release.md).

## Testing the Linux bundles from Windows, via WSL

The Windows dev machine can build and smoke-test the **linux-x86_64** bundles
directly, using WSL, without a separate Linux box. This needs Python 3.11 on
the WSL side (`python3` there defaults to whatever Ubuntu ships — 3.12 on
24.04 "noble" — which is too new for the pinned wheels).

**Install Python 3.11 in WSL once, ahead of time:**

```bash
sudo add-apt-repository -y ppa:deadsnakes/ppa && sudo apt update && sudo apt install -y python3.11 python3.11-venv
```

Ubuntu 24.04's own repos have no `python3.11` package at all (`apt-cache
policy python3.11` returns nothing) — deadsnakes is the standard source for
it. This needs a password prompt for `sudo`, so run it interactively yourself
rather than through a non-interactive tool call.

**Do the work from WSL's native filesystem, not the Windows-mounted repo.**
A worktree's `.git` file points at an absolute Windows path
(`E:\Git\...\.git\worktrees\...`), which WSL's `git` cannot resolve — `git
archive`, which `tools/build_release.py` depends on, fails immediately with
`fatal: not a git repository`. Clone into WSL's own filesystem instead (also
faster than building on a DrvFs-mounted path):

```bash
git clone /mnt/e/Git/latticegen2 ~/build/latticegen2 && cd ~/build/latticegen2 && git checkout <commit-or-branch>
```

**Build and smoke-test:**

```bash
python3.11 tools/build_release.py --flavour both
python3.11 -m venv ~/build/testenv
~/build/testenv/bin/python -m pip install --no-index --find-links dist/wheels -r requirements-bundle.txt
PATH="/path/to/python3.11/bin:$PATH" ~/build/testenv/bin/python tools/smoke_bundle.py --platform linux-x86_64
```

Two things that are easy to trip on:

* `tools/smoke_bundle.py` imports `tools/verify_geometry.py`, which needs
  `numpy`/`OCP` — the same "development interpreter" requirement `release.md`
  describes for the host build. Running it with a bare `python3.11` that has
  no packages installed fails at that import, not at anything the gate is
  actually checking. Use the venv above.
* The **wheels** flavour's own launcher and `install.sh` require a `python3`
  on `PATH` that is exactly 3.11 (by design — see `docs/release.md`'s
  reproducibility notes). If WSL's default `python3` is newer, the bundle
  correctly refuses to install/run, and the smoke gate reports that as a
  failure — which is the gate doing its job, not a bug. Put the 3.11
  interpreter first on `PATH` for the smoke run so it matches what CI's
  runner provides.

**Unit tests and e2e** run the same way as anywhere else, just under the venv
interpreter:

```bash
~/build/testenv/bin/python -m pytest test -q
~/build/testenv/bin/python tools/e2e.py
```

## Verification checklist

1. `python -m pytest test -q` passes. (Python has no separate lint step here;
   test failures are the signal.)
2. Edge cases for any modified boundary logic, in particular:
   - Parameter bounds (`-cc`/`-t` at 0.4/50/20) and the `t < a` cross-constraint.
   - Classification margin edges (a strut exactly at `r + d`).
   - Cap integrity at high `t/cc` ratios (docs/algorithm.md §3.3).
   - The interior shell's edge-use tally: every edge twice, opposite directions.
   - The floating-body rule's three outcomes: dropped, kept-because-connected,
     and kept-because-large.
3. `python tools/e2e.py` passes.
4. If packaging, launchers or dependencies changed: `python tools/build_release.py`
   then `python tools/smoke_bundle.py --platform <host>` passes.
5. Run `/code-review` (per CLAUDE.md) and address findings before pushing.
