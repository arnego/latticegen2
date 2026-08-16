# Third-party libraries and their licenses

Cross-reference between every third-party dependency latticegen2 uses and the
license text held in this directory (specification.md §2).

The versions below are the ones the release bundles are built against, and they
must stay in sync with [requirements-bundle.txt](../requirements-bundle.txt). A
mismatch between the two is a release blocker: the license texts here are only
meaningful if they describe the bytes actually shipped. See
[docs/release.md](../docs/release.md).

| Library | Version used | Role | License | Text |
|---|---|---|---|---|
| **Open CASCADE Technology (OCCT)** | 7.9.3 | The geometry kernel. STEP read/write, boolean intersection, sewing, meshing, B-rep validity checking. Reached through OCP; not vendored separately — the shared libraries ship inside the `cadquery-ocp` wheel. | LGPL-2.1 (with the Open CASCADE exception) | [occt-LICENSE_LGPL_21.txt](occt-LICENSE_LGPL_21.txt) |
| **OCP** (`cadquery-ocp`) | 7.9.3.1.1 | Python bindings for OCCT. The binding layer itself; it also redistributes the OCCT binaries above. | Apache-2.0 | [ocp-LICENSE_Apache_20.txt](ocp-LICENSE_Apache_20.txt) |
| **OCP proxy** (`cadquery-ocp-proxy`) | 7.9.3.1.1 | A 3 kB shim package `cadquery-ocp` requires; carries no binaries of its own. | Apache-2.0 | [ocp-LICENSE_Apache_20.txt](ocp-LICENSE_Apache_20.txt) |
| **VTK** | 9.6.2 | **Not used by latticegen2.** Present only because `cadquery-ocp` requires it and the OCP extension module links against its shared libraries — see the note below. | BSD-3-Clause | [vtk-LICENSE.txt](vtk-LICENSE.txt) |
| **NumPy** | 2.4.6 | Vectorised array maths for classification, spatial indexing and the lattice basis. | BSD-3-Clause | [numpy-LICENSE.txt](numpy-LICENSE.txt) |
| **psutil** | 7.1.0 | Machine resource detection for the CLI's `--cores`/`--ram` budgets (`src/latticegen2/sysinfo.py`). Core count comes from the stdlib; physical memory does not, and `--ram` needs it as both its upper bound and its default. Not on the geometry path. | BSD-3-Clause | [psutil-LICENSE.txt](psutil-LICENSE.txt) |

## Why VTK is redistributed despite being unused

latticegen2 never imports VTK; nothing in `src/` references it. It is shipped
because it cannot be removed, which was established by measurement rather than
assumed:

* `cadquery-ocp` 7.9.3.1.1 declares `vtk==9.6.2` as a hard requirement.
* Installing without it fails at import — `OCP/__init__.py` calls
  `os.add_dll_directory(<site-packages>/vtk.libs)` with no `isdir()` guard.
* Creating an empty `vtk.libs/` directory to satisfy that call gets one line
  further and then fails with `ImportError: DLL load failed while importing
  OCP`. The OCP extension module genuinely links against the VTK libraries.
* Downgrading is not an escape: `cadquery-ocp` 7.8.1.1 requires `vtk==9.3.1`.

It is therefore a redistribution obligation, and its BSD-3-Clause text is
included here accordingly. It is also roughly 60% of a release bundle's size
(328 MB installed).

## Bundled Python runtime (portable bundles only)

The **portable** release bundles carry a relocatable CPython interpreter under
`runtime/`, taken unmodified from
[python-build-standalone](https://github.com/astral-sh/python-build-standalone).
For those bundles Python is a redistributed component, not merely a
prerequisite, and the distribution's own license documents — covering CPython
under the PSF License Agreement plus the third-party components it embeds, such
as OpenSSL, zlib and libffi — are copied verbatim into `licenses/runtime/` by
`tools/build_release.py` at build time.

The **wheels** bundles and a plain source checkout bundle no interpreter. There,
Python remains a runtime prerequisite covered by the PSF License Agreement that
ships with whatever interpreter the operator installed.

## LGPL compliance note

OCCT is LGPL-2.1. latticegen2 uses it as a shared library through OCP's bindings
and does not modify or statically link it, so the LGPL's relinking obligation is
satisfied by the fact that the OCCT binaries are the stock ones distributed in
the `cadquery-ocp` wheel and can be replaced by installing a different build of
that wheel.

This argument has to remain true of the *release bundles*, not just of a
developer checkout, and it does:

* The **wheels** bundle ships the unmodified `cadquery-ocp` wheel in `wheels/`,
  and its `install` script is ordinary `pip`.
* The **portable** bundle ships the same stock binaries installed under
  `runtime/`, together with a working `pip`, so the recipient can install a
  different build of the wheel over them. The wheel files themselves are not
  duplicated inside this bundle — they would add 140–230 MB of content already
  present in installed form — but the wheels bundle for the same platform and
  version is published alongside it in the same GitHub release, so the exact
  wheel remains available to anyone who wants it.

Nothing in either bundle is a modified OCCT build, and nothing prevents
replacing it.

No dependency carries a GPL copyleft obligation: the set is LGPL-2.1,
Apache-2.0 and BSD-3-Clause.

Note that psutil is *also* used by `tools/profile_run.py`, which is a
development-only tool and is `export-ignore`d out of every bundle. That is a
separate use from the runtime one above and does not affect this table: psutil
is listed here because `src/latticegen2/sysinfo.py` imports it at run time, so
it ships whether or not `tools/` does.
