# tools/build_app.jl
#
# Builds the standalone-executable distribution of latticegen2 using
# PackageCompiler.jl (specification.md §2: "Using PackageCompiler.jl a
# standalone-executable distribution for windows must be created").
#
# This is a build-time tool, not part of the offline runtime: it must be run
# from a machine with network access (to fetch PackageCompiler.jl itself, once)
# and takes on the order of 15-45+ minutes, producing a multi-hundred-MB output
# tree under build/latticegen2-app/. Neither of those properties apply to the
# resulting app — the compiled executable itself runs fully offline, per
# specification.md §2's "zero network access" runtime requirement.
#
# Usage (from the repository root):
#   julia --project=tools/build -e "using Pkg; Pkg.instantiate()"   # once
#   julia --project=tools/build tools/build_app.jl
#
# Output: build/latticegen2-app/bin/latticegen2.exe (Windows) plus its bundled
# Julia runtime, sysimage, and artifacts (gmsh_jll/OCCT_jll — the same
# binaries the normal `julia --project=. src/main.jl` invocation uses, just
# bundled instead of resolved from ~/.julia). build/ is gitignored: this
# script's output is a local build artifact, not a checked-in binary.

using PackageCompiler

const ROOT = normpath(joinpath(@__DIR__, ".."))
const APP_DIR = joinpath(ROOT, "build", "latticegen2-app")

println("latticegen2: building standalone app")
println("  package dir: ", ROOT)
println("  app dir:     ", APP_DIR)
println()
println("This mirrors the entry point `julia_main()` in src/app_entry.jl, which")
println("delegates to src/main.jl (docs/algorithm.md §12) — the built executable")
println("takes the same CLI arguments as `latticegen2.bat`/`latticegen2.sh`.")
println()

isdir(APP_DIR) && (println("Removing previous build at ", APP_DIR); rm(APP_DIR; recursive=true, force=true))

create_app(
    ROOT, APP_DIR;
    executables=["latticegen2" => "julia_main"],
    force=true,
    incremental=false,        # fresh sysimage: reproducible, no stale-state risk (priority #1, precision)
    filter_stdlibs=false,     # keep the full stdlib set; latticegen2 uses Distributed, Dates, Printf, LinearAlgebra
)

println()
println("Build complete: ", joinpath(APP_DIR, "bin", Sys.iswindows() ? "latticegen2.exe" : "latticegen2"))
println("Run it exactly like the wrapper scripts, e.g.:")
println("  ", joinpath(APP_DIR, "bin", "latticegen2.exe"), " -i <input.step> -cc <mm> -t <mm> [options]")
