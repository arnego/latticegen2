# Third-party libraries and their licenses

Cross-reference between every third-party dependency latticegen2 uses and the
license text held in this directory (specification.md §2).

| Library | Version used | Role | License | Text |
|---|---|---|---|---|
| **Open CASCADE Technology (OCCT)** | 7.9.3 | The geometry kernel. STEP read/write, boolean intersection, sewing, meshing, B-rep validity checking. Reached through OCP; not vendored separately — the shared libraries ship inside the `cadquery-ocp` wheel. | LGPL-2.1 (with the Open CASCADE exception) | [occt-LICENSE_LGPL_21.txt](occt-LICENSE_LGPL_21.txt) |
| **OCP** (`cadquery-ocp`) | 7.9.3.1.1 | Python bindings for OCCT. The binding layer itself; it also redistributes the OCCT binaries above. | Apache-2.0 | [ocp-LICENSE_Apache_20.txt](ocp-LICENSE_Apache_20.txt) |
| **NumPy** | 2.4.6 | Vectorised array maths for classification, spatial indexing and the lattice basis. | BSD-3-Clause | [numpy-LICENSE.txt](numpy-LICENSE.txt) |

Python itself is a runtime prerequisite rather than a bundled dependency, and is
covered by the PSF License Agreement that ships with the interpreter.

## LGPL compliance note

OCCT is LGPL-2.1. latticegen2 uses it as a shared library through OCP's bindings
and does not modify or statically link it, so the LGPL's relinking obligation is
satisfied by the fact that the OCCT binaries are the stock ones distributed in
the `cadquery-ocp` wheel and can be replaced by installing a different build of
that wheel.

No dependency carries a GPL copyleft obligation: the set is LGPL-2.1,
Apache-2.0 and BSD-3-Clause.
