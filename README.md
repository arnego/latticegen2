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

**Nothing contacts the network at runtime.** Installation does, once.

## Installation

Online, on a machine that can reach PyPI:

```bash
python -m pip install numpy cadquery-ocp
```

Offline (the deployment target — specification.md §2), download the wheels on a
connected machine and carry them over:

```bash
python -m pip download numpy cadquery-ocp -d wheels
```

then, on the offline workstation:

```bash
python -m pip install --no-index --find-links wheels numpy cadquery-ocp
```

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

Full detail is in [docs/algorithm.md](docs/algorithm.md); the reasoning behind
the design is in
[docs/research/perf-rearchitecture-proposal.md](docs/research/perf-rearchitecture-proposal.md).

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
which is what the removed boolean used to do for free.

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
export. Measured: **305 MB** for the 80 mm test ball at `cc=20, t=4`, **1.06 GB**
for the test cylinder at `cc=10, t=1.5` (~14 k interior faces plus 1.1 k trimmed
boundary pieces). It scales roughly linearly with face count, which scales with
the number of lattice nodes. Worker processes each hold one junction and the
input body, and are negligible by comparison.

### Performance

Measured on a 6-core / 32 GB Windows workstation, against the previous
Julia/gmsh implementation (now in [`old-julia/`](old-julia/)):

| Scenario | Julia/gmsh | This implementation |
|---|---|---|
| 80 mm ball, `cc=20 t=4` | 13 m 52 s | **7.0 s** |
| test cylinder, `cc=10 t=1.5` | 25 m 55 s | **1 m 01 s** |

Both produce the same geometry: against the golden samples the
symmetric-difference volume is 0 mm³.

Output compactness is also at parity — 15,966 faces / 52.6 MB for the test
cylinder, against the old pipeline's 15,969 / 53.1 MB. That does not come for
free from instancing, which merges nothing; it comes from a same-domain
unification pass before export, which recovers what the removed boolean used to
do as a side effect. See "How it works".

---

## Changes from the Julia implementation

The command-line surface, log format, exit codes and STEP metadata were treated
as a starting point to improve on rather than a contract to reproduce (a
deliberate decision — see ground rule 6 of
[the implementation guide](docs/research/perf-rearchitecture-implementation-guide.md)).
Everything that changed:

1. **`--tile-cells` is removed.** It sized the tiles of a tiling-and-fusion stage
   that no longer exists.
2. **`--workers`, `--cores` and `--ram` are optional hints, not a mandatory
   either-or pair.** The old CLI required either `--cores` *and* `--ram`, or
   `--workers` *and* `--tile-cells`, because tile sizing genuinely could not be
   guessed safely. Boundary-junction jobs are constant-size and independent, so a
   worker count follows from the machine. `--workers` still overrides everything.
3. **`--ram` is advisory.** It is recorded in the run log. There is no tile-sizing
   calculation left for it to feed, and no memory watchdog: the stage that used to
   need backpressure (a distributed assembly holding every tile at once) is gone.
4. **The `t < a` cross-constraint is the only one.** An intermediate revision of
   this implementation also rejected `t >= cc/2`, expecting the junction's
   mid-strut interfaces to be swallowed by the node overlap. Measurement showed
   they are not, at any valid parameters, so that restriction was removed rather
   than shipped — see [docs/algorithm.md](docs/algorithm.md) §3.
5. **Log lines are reformatted** (timestamped lines, one block per stage, a single
   summary table). Every field specification.md §3 requires is still there:
   parameters, start time, duration, run characteristics, peak memory, output
   path.
6. **Exit code 5 (resource limits) is retained but currently unreachable**, since
   the memory watchdog it belonged to is gone. It is kept rather than renumbered
   so existing codes keep their meanings.
7. **`filter_floating!`'s exit-4 "unresolvable fragment" case cannot occur.**
   Connectivity is proven by construction now, so there is no case where the tool
   has to refuse to guess. The exit code still covers other processing failures.
8. **New: exact B-rep validation.** Every output solid is checked with OCCT's
   `BRepCheck_Analyzer` before the run reports success. The gmsh-based
   implementation could not express this check at all.

The verification tooling moved with it: `tools/verify_geometry.jl` and
`tools/e2e.jl` are now `tools/verify_geometry.py` and `tools/e2e.py`.

## Repository layout

| Path | Contents |
|---|---|
| `src/latticegen2/` | The package |
| `src/main.py` | Runnable entry point for a bare checkout |
| `test/` | pytest suite and the STEP test assets |
| `tools/` | `e2e.py`, `verify_geometry.py`, and the Phase-0 `prototypes/` |
| `docs/` | Specification, algorithm, testing guide, and the re-architecture research |
| `licenses/` | Dependency license texts |
| `old-julia/` | The superseded Julia/gmsh implementation, kept for reference |

## Testing

```bash
python -m pytest test -q      # unit tests
python tools/e2e.py           # end-to-end scenarios
```

See [docs/testing.md](docs/testing.md).
