# Third-party libraries and their licenses

Cross-reference between every third-party dependency latticegen2 uses and the
license text held in this directory (specification.md §2).

The versions below are the ones the release bundles are built against, and they
must stay in sync with [requirements-bundle.txt](../requirements-bundle.txt). A
mismatch between the two is a release blocker: the license texts here are only
meaningful if they describe the bytes actually shipped. See
[docs/release.md](../docs/release.md).

**One row cannot be checked that way, and it is worth knowing which.** Tcl/Tk is
not a wheel: it arrives inside the python-build-standalone interpreter, so it is
outside `requirements-bundle.txt` and outside every check the dependency pins
get. What stands in for that is the build-time import gate in
`tools/build_release.py` — it loads `tkinter` in the freshly populated runtime
and prints the Tk version into `BUNDLE-INFO.txt`, so the version recorded in a
bundle is one that was actually observed rather than one that was assumed.

| Library | Version used | Role | License | Text |
|---|---|---|---|---|
| **Open CASCADE Technology (OCCT)** | 7.9.3 | The geometry kernel. STEP read/write, boolean intersection, sewing, meshing, B-rep validity checking. Reached through OCP; not vendored separately — the shared libraries ship inside the `cadquery-ocp` wheel. | LGPL-2.1 (with the Open CASCADE exception) | [occt-LICENSE_LGPL_21.txt](occt-LICENSE_LGPL_21.txt) |
| **OCP** (`cadquery-ocp`) | 7.9.3.1.1 | Python bindings for OCCT. The binding layer itself; it also redistributes the OCCT binaries above. | Apache-2.0 | [ocp-LICENSE_Apache_20.txt](ocp-LICENSE_Apache_20.txt) |
| **OCP proxy** (`cadquery-ocp-proxy`) | 7.9.3.1.1 | A 3 kB shim package `cadquery-ocp` requires; carries no binaries of its own. | Apache-2.0 | [ocp-LICENSE_Apache_20.txt](ocp-LICENSE_Apache_20.txt) |
| **VTK** | 9.6.2 | **Not used by latticegen2.** Present only because `cadquery-ocp` requires it and the OCP extension module links against its shared libraries — see the note below. | BSD-3-Clause | [vtk-LICENSE.txt](vtk-LICENSE.txt) |
| **NumPy** | 2.4.6 | Vectorised array maths for classification, spatial indexing and the lattice basis. | BSD-3-Clause | [numpy-LICENSE.txt](numpy-LICENSE.txt) |
| **Tcl/Tk** | 8.6 (as built into the runtime) | The graphical front-end's toolkit, reached through the standard library's `tkinter`. **Portable bundles only** — see the note below. | BSD-style (Tcl/Tk licence terms) | [tcl-tk-LICENSE.txt](tcl-tk-LICENSE.txt) |

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

That distribution also carries **Tcl/Tk**, which the graphical front-end
reaches through the standard library's `tkinter`. It is a separate project from
CPython and is *not* covered by the PSF agreement or by CPython's incorporated-
software appendix, so its own terms are listed in the table above and copied
into `licenses/runtime/TCL-TK-LICENSE.txt` alongside the interpreter's.

The **wheels** bundles and a plain source checkout bundle no interpreter. There,
Python remains a runtime prerequisite covered by the PSF License Agreement that
ships with whatever interpreter the operator installed — and `tkinter` likewise
comes from that interpreter. On Debian and Ubuntu it is the separate
`python3-tk` package; without it the command line is unaffected and only the
front-end is unavailable.

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

## psutil is no longer a runtime dependency

Earlier revisions listed `psutil` in the table above, for the CLI's `--ram`
budget (`src/latticegen2/sysinfo.py` used it to read total and free physical
memory). That budget was removed outright — validated and recorded but never
enforced (specification.md §11) — and with it `sysinfo.py`'s only use of
`psutil`. Core count, the one budget that remains, has always come from the
stdlib.

`psutil` is still used by `tools/profile_run.py`, a development-only tool that
is `export-ignore`d out of every bundle, the same as `pytest`. Neither carries
a license text in this directory, for the same reason: nothing in a release
bundle imports either of them.
