# src/app_entry.jl
#
# PackageCompiler.jl standalone-executable entry point (specification.md §2,
# "Using PackageCompiler.jl a standalone-executable distribution for windows
# must be created"). Built by tools/build_app.jl / README.md "Building a
# standalone executable" — not used by the normal
# `julia --project=. src/main.jl <args>` invocation path.
#
# `PackageCompiler.create_app` requires an exported `julia_main()::Cint`
# function in the target package's main module as its executable entry
# point; it does not run arbitrary top-level script code the way
# `julia script.jl` does. `src/main.jl` is deliberately written as a
# top-level script rather than a function — its own header comment explains
# why: `Distributed.@everywhere` must run at true top level, not from inside
# a function body. Rather than restructure that (tested, working) driver,
# `julia_main` delegates to it via `include`, which evaluates its contents
# at the top level of the `Main` module — the same semantics a direct
# `julia src/main.jl <args>` invocation has (Julia's own script runner does
# the same `include(Main, path)` internally), so `main.jl` sees no
# difference between the two invocation paths.

"""
    julia_main() -> Cint

Standalone-app entry point. Delegates to `src/main.jl`, which reads `ARGS`
(populated from the compiled executable's own command-line arguments by
Julia's startup machinery, identically to a plain `julia script.jl <args>`
run) and calls `exit(code)` itself on every path per specification.md §3 /
docs/algorithm.md §9 — so this function's own `return` only matters if
`include` fails before `main.jl` reaches its own exit call.
"""
function julia_main()::Cint
    try
        Base.include(Main, joinpath(@__DIR__, "main.jl"))
        return 0
    catch e
        Base.showerror(stderr, e)
        println(stderr)
        return 1
    end
end
