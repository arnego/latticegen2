# old-julia — the superseded Julia/gmsh implementation

This directory holds the original latticegen2 implementation, moved aside
**unmodified** when the tool was re-architected around fuse-free lattice
synthesis (see [../docs/research/perf-rearchitecture-proposal.md](../docs/research/perf-rearchitecture-proposal.md)).
It is here for reference during the transition, not as a maintained parallel
product.

## Status

**Not maintained. Do not fix bugs or add features here.** The new
implementation at the repository root is the product; this code exists so the
old behaviour can be consulted while the new one settles.

Its end of life is the **`v1.0` tag on `main`**, which is the archive of record.
If this implementation is ever needed again, branch from that tag — this
directory is a convenience copy and will be deleted, not preserved.

## What is here

| Path | What it was |
|---|---|
| `src/*.jl` | The pipeline: CLI, lattice math, gmsh/OCCT kernel wrapper, classification, tiling, distributed assembly, STEP metadata, logging |
| `test/*.jl` | The Julia unit tests |
| `tools/e2e.jl`, `tools/verify_geometry.jl` | E2E harness and geometry checks |
| `tools/build_app.jl`, `tools/build/` | PackageCompiler standalone-executable build (never actually run to completion) |
| `Project.toml`, `Manifest.toml` | The Julia project and its pinned dependency set |
| `latticegen2.bat`, `latticegen2.sh` | Wrapper scripts that invoked `julia --project=. src/main.jl` |

## What is missing, deliberately

* **Dependency license texts.** `licenses/` now describes the dependencies the
  shipped tool actually has (OCCT via OCP, and NumPy). The Julia, Gmsh,
  Gmsh.jl-binding and PackageCompiler license texts were removed along with the
  dependencies themselves; they are on the `v1.0` tag with the code that needed
  them.
* **Working paths.** The tests and tools here refer to `src/…` and `test/…`
  relative to what used to be the repository root. Running them from this
  directory needs those paths adjusted, or a checkout of the `v1.0` tag.

## Why it was replaced

The Julia implementation spent about 91% of its wall time inside OCCT boolean
fusion, which grows at roughly N^2.5 in operand count, and it could not reach the
scale the project needs. The replacement removes fusion from the hot path
entirely by instancing one pre-fused junction solid per lattice node and
stitching instances along their shared mid-strut faces. On the `dense-lattice`
scenario the same output takes **1m 17s** instead of **25m 55s**.

The detailed analysis, including the parts of this implementation that were
*right* and were carried over intact (the classification margin analysis, the
mesh-coverage gate, the STEP header rewrite, the floating-body rule), is in
[../docs/research/perf-rearchitecture-proposal.md](../docs/research/perf-rearchitecture-proposal.md)
and the history sections of [../docs/algorithm.md](../docs/algorithm.md).
