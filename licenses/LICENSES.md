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
| **Matplotlib** (`matplotlib`) | 3.10.6 | **Not used.** Required by `vtk`. Its wheel carries the first three texts. The rest cover fonts, colour maps and scripts it ships, from its source tree's `LICENSE/` at `v3.10.6`, and the code it compiles in is listed in the next table. | Matplotlib License (PSF-based); bundled fonts and assets under their own terms | [matplotlib-LICENSE.txt](matplotlib-LICENSE.txt), [matplotlib-LICENSE_DEJAVU.txt](matplotlib-LICENSE_DEJAVU.txt), [matplotlib-LICENSE_STIX.txt](matplotlib-LICENSE_STIX.txt), [matplotlib-LICENSE_AMSFONTS.txt](matplotlib-LICENSE_AMSFONTS.txt), [matplotlib-LICENSE_BAKOMA.txt](matplotlib-LICENSE_BAKOMA.txt), [matplotlib-LICENSE_CARLOGO.txt](matplotlib-LICENSE_CARLOGO.txt), [matplotlib-LICENSE_COLORBREWER.txt](matplotlib-LICENSE_COLORBREWER.txt), [matplotlib-LICENSE_COURIERTEN.txt](matplotlib-LICENSE_COURIERTEN.txt), [matplotlib-LICENSE_JSXTOOLS_RESIZE_OBSERVER.txt](matplotlib-LICENSE_JSXTOOLS_RESIZE_OBSERVER.txt), [matplotlib-LICENSE_QT4_EDITOR.txt](matplotlib-LICENSE_QT4_EDITOR.txt), [matplotlib-LICENSE_SOLARIZED.txt](matplotlib-LICENSE_SOLARIZED.txt), [matplotlib-LICENSE_YORICK.txt](matplotlib-LICENSE_YORICK.txt) |
| **ContourPy** (`contourpy`) | 1.3.3 | **Not used.** Required by `matplotlib`. | BSD-3-Clause | [contourpy-LICENSE.txt](contourpy-LICENSE.txt) |
| **Cycler** (`cycler`) | 0.12.1 | **Not used.** Required by `matplotlib`. | BSD-3-Clause | [cycler-LICENSE.txt](cycler-LICENSE.txt) |
| **fontTools** (`fonttools`) | 4.60.1 | **Not used.** Required by `matplotlib`. Its external notices cover the font data and code it embeds. | MIT; embedded: SIL OFL 1.1, Apache-2.0, MIT | [fonttools-LICENSE.txt](fonttools-LICENSE.txt), [fonttools-LICENSE.external.txt](fonttools-LICENSE.external.txt) |
| **Kiwi** (`kiwisolver`) | 1.4.9 | **Not used.** Required by `matplotlib`. | BSD-3-Clause | [kiwisolver-LICENSE.txt](kiwisolver-LICENSE.txt) |
| **packaging** (`packaging`) | 25.0 | **Not used.** Required by `matplotlib`. | Apache-2.0 or BSD-2-Clause, at the recipient's choice | [packaging-LICENSE.txt](packaging-LICENSE.txt), [packaging-LICENSE.APACHE.txt](packaging-LICENSE.APACHE.txt), [packaging-LICENSE.BSD.txt](packaging-LICENSE.BSD.txt) |
| **Pillow** (`pillow`) | 11.3.0 | **Not used.** Required by `matplotlib`. Each platform's wheel appends the licences of the libraries it embeds to its own licence file, so the two differ; both are kept as the wheels ship them. The Windows file lacks two codecs compiled into libavif, listed in the next table. | MIT-CMU (HPND); embedded libraries under their own terms | [pillow-LICENSE-windows.txt](pillow-LICENSE-windows.txt), [pillow-LICENSE-linux.txt](pillow-LICENSE-linux.txt) |
| **pyparsing** (`pyparsing`) | 3.2.5 | **Not used.** Required by `matplotlib`. | MIT | [pyparsing-LICENSE.txt](pyparsing-LICENSE.txt) |
| **dateutil** (`python-dateutil`) | 2.9.0.post0 | **Not used.** Required by `matplotlib`. | Apache-2.0 and BSD-3-Clause (dual, by contribution date) | [python-dateutil-LICENSE.txt](python-dateutil-LICENSE.txt) |
| **six** (`six`) | 1.17.0 | **Not used.** Required by `python-dateutil`. | MIT | [six-LICENSE.txt](six-LICENSE.txt) |
| **NumPy** | 2.4.6 | Vectorised array maths for classification, spatial indexing and the lattice basis. Each platform's wheel appends the licences of its embedded OpenBLAS and, on Linux, the GCC Fortran runtime; the code it compiles in keeps its own texts in `numpy-components/`, as the wheel does. | BSD-3-Clause; OpenBLAS BSD-3-Clause; libgfortran GPL-3 with the GCC Runtime Library Exception; libquadmath LGPL-2.1 (Linux) | [numpy-LICENSE-windows.txt](numpy-LICENSE-windows.txt), [numpy-LICENSE-linux.txt](numpy-LICENSE-linux.txt), [numpy-components/](numpy-components/) |
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

**What VTK requires is shipped too, for the same reason.** `vtk` 9.6.2 declares
`matplotlib`, and `matplotlib` nine more: the ten rows marked *Required by* in
the table. latticegen2 imports none of them. Until they were pinned in
`requirements-bundle.txt`, `pip download` took whatever version PyPI served on
the day, so the bundle held wheels whose licences this directory did not
record. `tools/build_release.py` now refuses a download holding any wheel not
pinned there, and `test/test_release_pins.py` holds every pin to a row of the
table above at the same version, with its licence files present.

## Libraries compiled into the wheels

A wheel with compiled code can carry other projects' code inside it, as a
shared library or linked in statically. Their licences ask for their notices to
travel with the binary, so each is listed here. Where a wheel's own licence file
already carries a library's text, that file is the one linked above. These are
the ones a wheel leaves out, taken from upstream at the version compiled in.
[`tools/collect_licences.py`](../tools/collect_licences.py) refetches all of
them when a pin moves.

| Library | Version | Compiled into | License | Text |
|---|---|---|---|---|
| **FreeType** | 2.6.1 | `matplotlib` (`ft2font`), both platforms. Matplotlib pins this version for reproducible text rendering (`subprojects/freetype-2.6.1.wrap`). | FreeType License (FTL), chosen over its GPL-2.0 alternative; asks for a credit in the documentation, given in README.md | [matplotlib-freetype-2.6.1-FTL.txt](matplotlib-freetype-2.6.1-FTL.txt) |
| **Qhull** | 8.0.2 (2020.2) | `matplotlib` (`_qhull`), both platforms (`subprojects/qhull.wrap`) | Qhull License | [matplotlib-qhull-8.0.2-COPYING.txt](matplotlib-qhull-8.0.2-COPYING.txt) |
| **Anti-Grain Geometry** | 2.4 | `matplotlib` (`_backend_agg`, `_image`), both platforms (`extern/agg24-svn`) | Agg permissive notice, carried in each source file's header | [matplotlib-agg-2.4-LICENSE.txt](matplotlib-agg-2.4-LICENSE.txt) |
| **pybind11** | 2.13.2 or later (not recorded in the wheels) | `matplotlib`, `contourpy`, both platforms; header-only, so compiled in | BSD-3-Clause, the same text throughout 2.x | [pybind11-LICENSE.txt](pybind11-LICENSE.txt) |
| **cppy** | 1.3.0 or later (not recorded in the wheel) | `kiwisolver`, both platforms; header-only | BSD-3-Clause | [cppy-LICENSE.txt](cppy-LICENSE.txt) |
| **AOM** (libaom) | 3.12.1 | `pillow`, inside libavif 1.3.0 (`AVIF_CODEC_AOM=LOCAL`). On Windows only this text is missing; the Linux wheel's own file carries it. | BSD-2-Clause, plus the Alliance for Open Media Patent License 1.0 | [pillow-libavif-aom-3.12.1-LICENSE.txt](pillow-libavif-aom-3.12.1-LICENSE.txt), [pillow-libavif-aom-3.12.1-PATENTS.txt](pillow-libavif-aom-3.12.1-PATENTS.txt) |
| **dav1d** | 1.5.1 | `pillow`, inside libavif 1.3.0 (`AVIF_CODEC_DAV1D=LOCAL`). As for AOM, missing from the Windows file only. | BSD-2-Clause | [pillow-libavif-dav1d-1.5.1-COPYING.txt](pillow-libavif-dav1d-1.5.1-COPYING.txt) |

**What the Pillow wheels carry themselves**, for reference:
- **Both platforms:** Brotli, FreeType 2.13.3, HarfBuzz, LittleCMS, libavif (with libyuv), libjpeg-turbo, libpng, libwebp (with sharpyuv), OpenJPEG, libtiff and xz.
- **Windows:** zlib-ng.
- **Linux:** zlib, bzip2, libxcb and libXau, and the AOM and dav1d texts above.

The Windows wheel's FreeType section also prints the GPL-2.0 text. That is FreeType's alternative licence, and the FTL in the same section is the one taken.

**Not covered yet: `cadquery-ocp` and `vtk`.** They embed far more than the
packages above. `cadquery_ocp.libs` holds FreeImage, OpenEXR and Imath, LibRaw,
Lerc, libdeflate, OpenJPH, zstd and its own copies of the image libraries, and
`vtk.libs` holds about 200 libraries. Neither wheel carries texts for these.
Two of OCP's need a decision rather than a copy:
- FreeImage is offered under GPL-2.0, GPL-3.0 or the FreeImage Public License.
- LibRaw is offered under LGPL-2.1 or CDDL-1.0.

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

No dependency covered above carries a GPL copyleft obligation:
- FreeType is taken under the FTL.
- NumPy's libgfortran is GPL-3 with the GCC Runtime Library Exception, which exists so that programs linking it are not bound by the GPL.
- NumPy's libquadmath, on Linux, is LGPL-2.1. It ships as a separate shared library in `numpy.libs/`, replaceable on the same terms as OCCT above.
- The rest are permissive: Apache-2.0, BSD-2/3-Clause, MIT, MIT-CMU, the PSF-based Matplotlib License, the Qhull and Agg notices, the AOM patent licence, and font licences (SIL Open Font License, Bitstream Vera and others).

`cadquery-ocp` and `vtk`'s embedded libraries are not yet assessed (above).

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
