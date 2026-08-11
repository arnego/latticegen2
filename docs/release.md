# Release Manual

How to get from a working `main` to published, downloadable offline bundles.
This is the developer's procedure; the operator-facing instructions live inside
each bundle as `README-OFFLINE.md`.

This file is excluded from the bundles themselves via `.gitattributes`.

---

## What a release contains

Four assets plus a checksum file, built and verified by
[`.github/workflows/release.yml`](../.github/workflows/release.yml):

| Asset | Needs Python on the target? |
|---|---|
| `latticegen2-<ver>-win64-portable.zip` | no — carries its own interpreter |
| `latticegen2-<ver>-win64-wheels.zip` | yes, exactly Python 3.11 x86-64 |
| `latticegen2-<ver>-linux-x86_64-portable.tar.gz` | no |
| `latticegen2-<ver>-linux-x86_64-wheels.tar.gz` | yes, exactly Python 3.11 x86-64 |
| `SHA256SUMS.txt` | — |

Linux uses `.tar.gz` rather than `.zip` on purpose: ZIP loses the executable
bit in most extractors, and a portable Linux bundle whose `runtime/bin/python3`
is not executable is inert.

Measured on Windows: **182 MB** portable, **158 MB** wheels. Linux is larger —
its `vtk` wheel is 146 MB against Windows' 81 MB — but has not been built yet at
the time of writing; the first CI run will tell you.

About 60% of a bundle is VTK, which latticegen2 never calls but cannot drop. It
also drags in matplotlib and its dependency tree, taking the wheel set from 3
packages to 14. See [requirements-bundle.txt](../requirements-bundle.txt) for
the measurements behind that, including what happens if you try to remove it.

Builds are reproducible: archive entries are sorted and timestamps fixed, so
rebuilding the same commit gives a byte-identical archive (verified). That makes
the published SHA-256 a build fingerprint rather than just a transfer checksum.

---

## 1. Before you tag

```bash
python -m pytest test -q
```

```bash
python tools/e2e.py
```

Both must pass. Keep the per-run `.log` files `tools/e2e.py` produces and attach
them to the pull request, per [testing.md](testing.md).

Then check the two things a release can get wrong that testing will not catch:

- **Pins versus licences.** [`requirements-bundle.txt`](../requirements-bundle.txt)
  and [`licenses/libraries.md`](../licenses/libraries.md) name the same
  versions. If a dependency moved, both change together and the licence texts in
  `licenses/` are re-checked. A mismatch is a release blocker, not a tidiness
  issue: the licence texts are only meaningful if they describe the bytes
  actually shipped.
- **Docs describe what ships.** [`README.md`](../README.md) and
  [`specification.md`](specification.md) reflect reality.

## 2. Bump the version

Edit **one** file:

```
src/latticegen2/__init__.py     __version__ = "2.1.0"
```

`pyproject.toml` derives its version from that attribute, so there is no second
place to update. This is the step most likely to be done the old way out of
habit — there used to be a literal in `pyproject.toml`, and there is not any
more.

Commit on `main`:

```bash
git commit -am "Release 2.1.0"
```

## 3. Dry run (recommended, and required for the first release)

GitHub → **Actions** → **release** → **Run workflow** on `main`.

This runs the whole pipeline — tests, both bundles on both platforms, and the
smoke gate — and stops short of publishing. Download the run's artifacts and
spot-check one if you want a look before it becomes public.

## 4. Tag and push

```bash
git tag -a v2.1.0 -m "latticegen2 2.1.0"
```

```bash
git push origin main --follow-tags
```

The workflow refuses to build if the tag disagrees with `__version__`, so a
mistyped tag costs you a re-tag rather than a mislabelled release.

## 5. Reading a failure

Every step is a gate. None of them should be worked around.

| Failing step | What it means |
|---|---|
| *Check tag matches `__version__`* | The tag and `src/latticegen2/__init__.py` disagree. Delete the tag, fix one of them, re-tag. |
| *Install and run unit tests* | A real regression. Fix it; do not release. |
| *Build bundles* | Usually a wheel or interpreter that could not be downloaded, or a pin that no longer exists on PyPI. Check `requirements-bundle.txt` and `PBS_RELEASE` in `tools/build_release.py`. |
| *Smoke-test the built bundles* | **The bundle is broken.** This is the gate doing its job — it extracts the archive and runs a real generation. Never publish past it. |
| *publish* | Token permissions, or the release/tag already exists. |

## 6. After the release appears

1. Confirm all four assets plus `SHA256SUMS.txt` are attached.
2. Download **one portable bundle**, verify its hash, extract it on a machine
   that has never had the development environment, and run one generation. CI
   proves the archive works; this proves the *download and transfer* works, and
   it is the only check that does.
3. Edit the generated release notes if they need it.

## 7. Building bundles by hand

For an air-gapped prep machine, or when CI is unavailable. The script builds for
the host platform only, so a Linux bundle needs a Linux machine.

Online:

```bash
python tools/build_release.py --flavour both
```

Offline, with the inputs carried over on removable media:

```bash
python tools/build_release.py --flavour both --wheels-dir /media/usb/wheels --runtime /media/usb/cpython-3.11.15+20260807-x86_64-unknown-linux-gnu-install_only_stripped.tar.gz
```

To collect those inputs on a connected machine, run the online form once — it
leaves the wheels in `dist/wheels/` and the interpreter in `.runtime-cache/`.

Output lands in `dist/`. Then verify what you built:

```bash
python tools/smoke_bundle.py --platform win64
```

Two things to know about the builder:

- **It packages a git ref, not your working tree.** Uncommitted changes will not
  appear in the bundle. Commit first, or pass `--ref <sha>`.
- **It must run under the Python version the bundle targets** (3.11) when it
  downloads wheels, because wheels are version-specific. It says so plainly if
  you get this wrong. `--wheels-dir` bypasses the requirement.

## 8. Fixing a bad release

Delete the release and its tag, fix the problem, and re-tag with a **new patch
version**. Do not reuse a version number: assets someone may already have
downloaded must not change meaning underneath them.

```bash
gh release delete v2.1.0 --yes
```

```bash
git push --delete origin v2.1.0 && git tag -d v2.1.0
```

## 9. Withdrawing a release

Mark it as a pre-release, or delete it, from the GitHub releases page. Note that
bundles already downloaded are entirely self-contained and keep working — there
is no phone-home and nothing to revoke. Withdrawal stops new downloads; it does
not recall old ones. If a release is withdrawn for a correctness reason, say so
in the notes of the replacement.

---

## Changing what goes into a bundle

| To change | Edit |
|---|---|
| Which files ship | [`.gitattributes`](../.gitattributes) `export-ignore` list |
| Dependency versions | [`requirements-bundle.txt`](../requirements-bundle.txt) **and** [`licenses/libraries.md`](../licenses/libraries.md) |
| Bundled interpreter | `PBS_RELEASE` / `DEFAULT_TARGET_PYTHON` in [`tools/build_release.py`](../tools/build_release.py) |
| Operator instructions | `write_offline_readme` in the same file |
| What the gate checks | [`tools/smoke_bundle.py`](../tools/smoke_bundle.py) |
