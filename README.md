# latticegen2

Generates a parameterised diamond-strut lattice that fills the solid volume of an
input STEP file, and writes it back out as an exact B-rep AP214 STEP body in the
same coordinate system.

```bash
latticegen2.bat -i part.step -cc 10 -t 1.5
```

The lattice is a simple cubic grid of nodes rotated so its body diagonal points
along +Z — cube cells standing on a tip. Struts run along the cube edges, each
with a square profile turned 45° so one diagonal is horizontal (a "diamond"
cross-section). `-cc` sets the XY-plane distance between the bottom nodes of
adjacent cells and `-t` the strut profile's side length.

---

## Dependencies

| Dependency | Version | Why | License |
|---|---|---|---|
| Python | 3.11 or newer | Runtime | PSF |
| [OCP](https://github.com/CadQuery/OCP) (`cadquery-ocp`) | 7.7+ | Python bindings for the Open CASCADE (OCCT) geometry kernel — STEP I/O, booleans, sewing, meshing, exact validity checking | Apache-2.0 (OCP), LGPL-2.1 (OCCT) |
| NumPy | 1.24+ | Vectorised classification and indexing | BSD-3-Clause |

Nothing else. License texts are in [`licenses/`](licenses/), cross-referenced in
[`licenses/libraries.md`](licenses/libraries.md).

**Nothing contacts the network at runtime.** Installation does, once — unless
you use a release bundle, which never does.

## Releases (offline / air-gapped)

Each tagged release publishes ready-made offline bundles on the
[Releases page](../../releases), in two flavours per platform:

| Asset | Python needed on the target? | Use when |
|---|---|---|
| `latticegen2-<ver>-win64-portable.zip` | **no** | The machine has no Python, or you cannot install anything |
| `latticegen2-<ver>-win64-wheels.zip` | yes, exactly 3.11 x86-64 | The machine already has Python 3.11 |
| `latticegen2-<ver>-linux-x86_64-portable.tar.gz` | **no** | as above |
| `latticegen2-<ver>-linux-x86_64-wheels.tar.gz` | yes, exactly 3.11 x86-64 | as above |

The portable bundles carry their own interpreter with every dependency
installed: extract, run, no admin rights and no network at any point. Linux uses
`.tar.gz` because ZIP loses the executable bit that `runtime/bin/python3` needs.

After transferring the file, check it against the release's `SHA256SUMS.txt`:

```bash
sha256sum -c SHA256SUMS.txt
```

```bash
certutil -hashfile latticegen2-2.0.0-win64-portable.zip SHA256
```

Then extract and run — each bundle contains a `README-OFFLINE.md` with the
specifics for its flavour:

```bash
latticegen2.bat -i part.step -cc 10 -t 1.5
```

Every asset is extracted and run end to end by CI before it is published. The
release procedure is documented in [docs/release.md](docs/release.md).

## Installation

Not needed if you use a release bundle above. From a source checkout:

Online, on a machine that can reach PyPI:

```bash
python -m pip install numpy cadquery-ocp
```

Offline (the deployment target — specification.md §2), download the wheels on a
connected machine and carry them over. Use `requirements-bundle.txt` so you get
the exact versions the releases are built and tested against:

```bash
python -m pip download -r requirements-bundle.txt -d wheels
```

then, on the offline workstation:

```bash
python -m pip install --no-index --find-links wheels -r requirements-bundle.txt
```

Wheels are specific to a Python minor version, so download them with the same
version you will install into.

> `cadquery-ocp` requires `vtk`, which pulls in matplotlib and its dependencies.
> latticegen2 never imports any of it, but the OCP extension module links
> against the VTK libraries and will not load without them, so it cannot be
> omitted — see [requirements-bundle.txt](requirements-bundle.txt) for the
> measurements. It is most of the download size.

No install step is needed for latticegen2 itself — it runs straight from a
checkout via `src/main.py` or the wrapper scripts. To install it as a command
instead, `python -m pip install .` provides a `latticegen2` entry point.

## Usage

```bash
latticegen2.bat -i part.step -cc 10 -t 1.5 -v          # Windows
./latticegen2.sh -i part.step -cc 10 -t 1.5 -v         # Linux
python src/main.py -i part.step -cc 10 -t 1.5 -v       # directly
```

Set `LATTICEGEN2_PYTHON` if the interpreter you want is not the default `python`.

### Parameters

| Flag | Type | Required | Units | Range | Default | Description |
|---|---|---|---|---|---|---|
| `-i`, `--input` | path | yes | — | — | — | STEP file whose solid defines the lattice bounds |
| `-cc` | float | yes | mm | 0.4 – 50 | — | XY distance between the bottom nodes of adjacent cells |
| `-t` | float | yes | mm | 0.4 – 20 | — | Side length of the diamond strut profile |
| `-o`, `--output` | path | no | — | — | `<input_stem>-cc<cc>t<t>.step` | Output STEP path |
| `--workers` | int | no | count | 1 – 128 | from `--cores`, else from the machine | Worker processes for the boundary stage |
| `--cores` | int | no | count | 1 – 128 | detected | Physical cores available; `--workers` is derived as `min(cores-1, 8)` |
| `--ram` | float | no | GB | 1 – 1024 | — | Memory budget; recorded in the log |
| `-bg`, `--background` | flag | no | — | — | off | Run at below-normal process priority |
| `-v`, `--verbose` | flag | no | — | — | off | Verbose console output (a full `.log` is always written) |
| `-h`, `--help` | flag | no | — | — | — | Usage |

`-t` must be smaller than the cell edge `a = cc/√2`; a thicker strut cannot fit
inside one cell. That is the only cross-constraint.

### Output

* **`<output>.step`** — the lattice, AP214, millimetres, exact B-rep. Usually one
  solid; more if the input geometry trims struts into genuinely disconnected
  islands (specification.md §1). The STEP header carries the part name
  `<input_stem>+cc<cc>+t<t>` and the generation parameters.
* **`<output>.log`** — always written, with the full run header, per-stage
  timings, and the end-of-run summary. `-v` only raises *console* verbosity; it
  never changes what the log contains.

### Exit codes

| Code | Meaning |
|---|---|
| 0 | Success |
| 2 | Parameter validation failure (before any computation) |
| 3 | Input geometry could not be read, parsed, or faithfully meshed |
| 4 | Geometry processing failure |
| 5 | Resource limit |
| 6 | Output write failure or failed round-trip check |
| 130 | Cancelled by the user with Ctrl+C |

Every non-zero exit prints one human-readable reason line.

---

## How it works

Full detail is in [docs/algorithm.md](docs/algorithm.md); the measurements behind
the design choices are in
[tools/prototypes/RESULTS.md](tools/prototypes/RESULTS.md).

```
import → tessellate → classify → instance interior → trim boundary
       → connect → stitch → simplify → validate → write STEP → round-trip verify
```

The load-bearing idea is that **the lattice needs almost no boolean operations at
all.** The three strut directions are mutually orthogonal, so if every strut is
cut at its midpoint and each half assigned to its nearer node, every node owns
six half-struts whose union — the *junction solid* — is congruent at every node in
the lattice. All strut-strut overlap lives inside that solid; two adjacent
junctions never overlap in volume, they meet exactly on a shared square face.

So the whole interior is one solid, instanced. Neighbouring instances drop their
shared face and reference the *same* vertices and edges through a precomputed
index correspondence, which makes the result watertight by construction rather
than by tolerance. Only junctions that straddle the input surface need a boolean,
and each gets exactly one single-operand intersection, distributed across worker
processes.

Instancing does have one cost, and it is paid at the end. Because nothing is ever
merged, two neighbouring junctions leave their coplanar lateral faces sitting
side by side unmerged at every strut interface — twice as many faces as the
geometry needs. A same-domain unification pass before export merges them back,
halving the face count and the file size.

| Lever | Effect |
|---|---|
| Classify before intersecting | Booleans run only for the O(surface area) boundary junctions, not the O(volume) interior |
| One junction template, instanced | The only general fuse in the program is six operands, once per run (~40 ms) |
| Indexed shared-topology shell | Interior assembly is O(nodes) with exact watertightness — no sewing, no tolerance |
| One object operand per intersection | Makes OCCT's operand-fragmentation failure mode structurally unreachable |
| Connectivity by graph, not by boolean | The floating-body rule is a BFS over surviving interfaces |
| Sewing confined to the boundary | Stitching cost scales with surface area, not volume |
| Process-parallel boundary junctions | Constant-size independent jobs |
| Same-domain unification before export | Instancing merges nothing, so coplanar faces meet unmerged at every strut interface; unifying them halves the face count and file size — and makes the run *faster*, since export and verification then handle half as much |
| Measured, not assumed, mesh deviation | The classification margin is an upper bound on the mesher's real error |

### Memory

Peak memory is dominated by the finished B-rep held in the master process before
export. Measured: **270 MB** for the 80 mm test ball at `cc=20, t=4`, **762 MB**
for the test cylinder at `cc=10, t=1.5`. It scales roughly linearly with face
count, which scales with the number of lattice nodes. Worker processes each hold
one junction and the input body — 214 MB at peak on the cylinder — and the run
log reports their high-water mark alongside the master's.

### Performance

Measured on a 6-core / 32 GB Windows workstation:

| Scenario | Nodes (interior / boundary) | Faces | Output | Time |
|---|---|---|---|---|
| 80 mm ball, `cc=20 t=4` | 27 / 176 | 1,338 | 4.7 MB | **7 s** |
| test cylinder, `cc=10 t=1.5` | 594 / 968 | 15,966 | 52.6 MB | **1 m 01 s** |

Both match their golden samples with a symmetric-difference volume of 0 mm³, put
0 mm³ of material outside the input body, and pass `BRepCheck_Analyzer` with zero
self-intersections.

---

## Repository layout

| Path | Contents |
|---|---|
| `src/latticegen2/` | The package |
| `src/main.py` | Runnable entry point for a bare checkout |
| `test/` | pytest suite and the STEP test assets |
| `tools/` | `e2e.py`, `verify_geometry.py`, and `prototypes/` (standalone benchmarks for the design's load-bearing assumptions) |
| `docs/` | Specification, algorithm, and testing guide |
| `licenses/` | Dependency license texts |

## Testing

```bash
python -m pytest test -q      # unit tests
python tools/e2e.py           # end-to-end scenarios
```

See [docs/testing.md](docs/testing.md).
