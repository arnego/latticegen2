# Unit test entry point. Run with:
#   julia --project=. -e "using Pkg; Pkg.test()"
# or directly:
#   julia --project=. test/runtests.jl
#
# Pure-math tests (no gmsh, fast) are kept in a separate testset from
# gmsh-dependent tests (slower, require the OCCT kernel) per
# docs/algorithm.md §12 / the project plan, so the fast subset can be run
# on its own during iterative development. Included files are added here as
# each corresponding src/ module is implemented.

using Test

include(joinpath(@__DIR__, "..", "src", "latticegen2.jl"))
using .latticegen2

@testset "latticegen2" begin
    @testset "pure-math (no gmsh)" begin
        include("test_runlog.jl")
        include("test_cli.jl")
        include("test_lattice.jl")
        include("test_tiling.jl")
        include("test_stepmeta.jl")
    end

    @testset "gmsh-dependent" begin
        include("test_geomkernel.jl")
        include("test_classify.jl")
        include("test_pipeline.jl")
        include("test_cleanup.jl")
        include("test_verify_geometry.jl")
    end
end
