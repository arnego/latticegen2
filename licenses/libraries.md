# Third-Party Library License Cross-Reference

Per [specification.md](../docs/specification.md) §2 ("License text must be obtained and
put into /licenses folder, and /licenses/libraries.md must be updated with the cross
reference between the library used and the corresponding license text file valid for
that library"). License texts here were copied verbatim from the actual installed
package/artifact files on the development machine (paths noted below), not
retyped/reconstructed, so they are guaranteed to match exactly what is vendored into an
offline deployment's `~/.julia` depot.

| Library | Role in this project | Version (as pinned in `Manifest.toml`) | License | License file |
|---|---|---|---|---|
| [Gmsh.jl](https://github.com/JuliaFEM/Gmsh.jl) | Julia binding used to call the Gmsh SDK API | 0.3.1 | MIT | [gmsh-jl-binding-LICENSE.txt](gmsh-jl-binding-LICENSE.txt) |
| [gmsh (SDK)](https://gmsh.info/) | Geometry kernel driver: STEP/BREP I/O, meshing, model management. Bundled as a binary artifact (`gmsh_jll`) pulled in transitively by Gmsh.jl. | 4.15.2+0 | GPL-2.0-or-later, with an explicit exception permitting combination with OpenCASCADE (and others) under their own licenses | [gmsh-LICENSE.txt](gmsh-LICENSE.txt) |
| [OCCT (Open CASCADE Technology)](https://www.opencascade.com/) | B-rep geometry kernel: primitive construction, boolean operations (fuse/common/cut), STEP AP214 export. Bundled transitively as `OCCT_jll`, used internally by the gmsh SDK. | 7.9.3+0 | LGPL-2.1 (with the OCCT-specific exception permitting static linking without triggering LGPL's relinking requirement, per OCCT's own license text) | [occt-LICENSE_LGPL_21.txt](occt-LICENSE_LGPL_21.txt) |
| [Julia](https://julialang.org/) | Language runtime | 1.10.11 | MIT | [julia-LICENSE.md](julia-LICENSE.md) |
| [PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl) | Build-time only: produces the standalone-executable distribution (specification.md §2). Not a runtime dependency — declared in `tools/build/Project.toml`, a separate environment from the root `Project.toml`, so it is never resolved by an ordinary `julia --project=. src/main.jl` run. | 2.4.0 | MIT | [packagecompiler-LICENSE.txt](packagecompiler-LICENSE.txt) |

## Notes

- **Only Gmsh.jl is a direct runtime dependency** of this project (see `Project.toml`);
  `gmsh_jll` and `OCCT_jll` are transitive dependencies pulled in automatically by
  Gmsh.jl / `Pkg.instantiate()`, and everything else `latticegen2` uses (`Dates`,
  `Distributed`, `LinearAlgebra`, `Printf`, `Test`) is Julia standard library, covered by
  Julia's own MIT license above — no separate license file is needed for stdlib modules.
  PackageCompiler.jl (above) is a build-only dependency, isolated in
  `tools/build/Project.toml`; see [README.md](../README.md) "Building a standalone
  executable".
- **Offline deployment implication:** because gmsh is GPL-licensed (even with the
  combination exception) and OCCT is LGPL-licensed, `latticegen2` as a whole — being
  dynamically linked against, not statically embedding, these libraries via the
  standard `gmsh_jll`/`OCCT_jll` binary artifacts — remains compatible with keeping
  `src/` itself under a separate license if desired; this project does not currently
  declare its own license (see [specification.md](../docs/specification.md) §2,
  "License constraints: TBD"). Redistributing the compiled offline deployment bundle
  (the `~/.julia` depot copy described in [README.md](../README.md)) must include these
  license files, which is exactly what this `licenses/` folder is for.
- License files were copied on 2026-08-08 from the local Julia depot, from the exact
  artifact hashes pinned in `Manifest.toml`:
  - gmsh: `~/.julia/artifacts/688edbeef6af47e7819cd3283e18669870f6f6af/share/licenses/gmsh/LICENSE.txt` (Windows x86_64 build)
  - OCCT: `~/.julia/artifacts/78e4aecdcc8bf5308e4b38ace088f4fc64580daa/share/licenses/OCCT/LICENSE_LGPL_21.txt` (Windows x86_64 build)
  - Gmsh.jl: `~/.julia/packages/Gmsh/<hash>/LICENSE`
  - Julia: `<juliaup installation>/LICENSE.md`
  - PackageCompiler.jl: `~/.julia/packages/PackageCompiler/<hash>/LICENSE` (copied
    2026-08-08, version 2.4.0, resolved via `tools/build/Project.toml`)
