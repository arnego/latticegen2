# Unit/regression tests for src/pipeline.jl.

import Gmsh: gmsh

@testset "pipeline.jl" begin
    @testset "determine_worker_count" begin
        a1 = parse_args(["-i", "x.step", "-cc", "5", "-t", "1", "--workers", "7", "--tile-cells", "8"])
        @test determine_worker_count(a1) == 7

        a2 = parse_args(["-i", "x.step", "-cc", "5", "-t", "1", "--cores", "6", "--ram", "24"])
        @test determine_worker_count(a2) == clamp(6 - 1, 1, 8)

        a3 = parse_args(["-i", "x.step", "-cc", "5", "-t", "1", "--cores", "1", "--ram", "8"])
        @test determine_worker_count(a3) == 1   # clamp floor

        a4 = parse_args(["-i", "x.step", "-cc", "5", "-t", "1", "--cores", "100", "--ram", "8"])
        @test determine_worker_count(a4) == 8   # clamp ceiling
    end

    @testset "derive_tile_size" begin
        # More RAM / fewer workers -> larger (or equal) memory-bound tile,
        # with a generous, non-binding time probe so the memory term is
        # what's being exercised.
        fast_probe = (1.0, 1000)   # probe_seconds, probe_struts: cheap/fast
        n_small_ram = derive_tile_size(1000.0, fast_probe..., 4, 1.0)
        n_large_ram = derive_tile_size(1000.0, fast_probe..., 4, 100.0)
        @test n_large_ram.n >= n_small_ram.n
        @test 2 <= n_small_ram.n <= 8
        @test 2 <= n_large_ram.n <= 8

        # Degenerate zero mem_per_strut / zero probe_seconds must not throw
        # or divide by zero.
        @test derive_tile_size(0.0, 1.0, 100, 4, 8.0).n >= 2
        @test derive_tile_size(1000.0, 0.0, 100, 4, 8.0).n >= 2

        # The hard cap: even with unlimited RAM and an instant probe, n
        # never exceeds 8 (docs/algorithm.md §7.1/§11.2 fuse-time knee).
        huge = derive_tile_size(1.0, 0.001, 1000, 1, 1024.0)
        @test huge.n <= 8

        # A slow probe (expensive fuse per strut) should push n_time, and
        # therefore the chosen n, down relative to a fast one at the same
        # memory budget.
        slow = derive_tile_size(1000.0, 60.0, 192, 4, 100.0)   # probe took as long as the whole target
        fast = derive_tile_size(1000.0, 0.1, 192, 4, 100.0)
        @test slow.n_time <= fast.n_time
    end

    @testset "shutdown_workers! is a safe no-op without a worker pool" begin
        # `workers()` returns `[1]` (the master) when no workers were ever
        # added, and `rmprocs` throws on the master — teardown must filter it
        # out and do nothing, since it runs on every exit path including a
        # plain successful sequential run.
        @test shutdown_workers!() === nothing
        @test shutdown_workers!() === nothing   # idempotent
    end

    @testset "_group_might_touch AABB pre-filter" begin
        with_gmsh() do
            new_model("group-might-touch")
            near1 = Int32(gmsh.model.occ.addBox(0.0, 0.0, 0.0, 1.0, 1.0, 1.0))
            near2 = Int32(gmsh.model.occ.addBox(0.5, 0.0, 0.0, 1.0, 1.0, 1.0))   # overlaps near1
            far = Int32(gmsh.model.occ.addBox(1000.0, 0.0, 0.0, 1.0, 1.0, 1.0))  # nowhere near either
            gmsh.model.occ.synchronize()
            boxes = bounding_boxes(Int32[near1, near2, far])

            @test latticegen2._group_might_touch(Int32[near1, near2], boxes)
            @test !latticegen2._group_might_touch(Int32[near1, far], boxes)
            @test !latticegen2._group_might_touch(Int32[far], boxes)   # single element: nothing to touch
        end
    end

    @testset "balanced_fuse! terminates and preserves disjoint solids" begin
        # Regression test for a real bug found during development: when a
        # group of solids consistently fails to fuse (or, as here, simply
        # never reduces in count because they're mutually disjoint),
        # re-grouping the *same* un-reduced elements into the *same*
        # batches every round looped forever. balanced_fuse! must detect a
        # no-progress round and stop rather than hang. The AABB pre-filter
        # should also mean this never even attempts an OCCT fuse call here
        # (all 12 boxes are mutually far apart), so it should be fast, not
        # just "eventually terminates".
        with_gmsh() do
            new_model("balanced-fuse-disjoint")
            tags = Int32[]
            for i in 0:11
                push!(tags, Int32(gmsh.model.occ.addBox(100.0 * i, 0.0, 0.0, 1.0, 1.0, 1.0)))
            end
            gmsh.model.occ.synchronize()

            elapsed = @elapsed result = balanced_fuse!(tags)
            @test elapsed < 5.0                  # AABB pre-filter: no OCCT calls needed at all
            @test length(result) == length(tags) # nothing to merge -> nothing merged
        end
    end

    @testset "balanced_fuse! actually fuses touching solids" begin
        with_gmsh() do
            new_model("balanced-fuse-touching")
            tags = Int32[]
            for i in 0:5
                push!(tags, Int32(gmsh.model.occ.addBox(0.5 * i, 0.0, 0.0, 1.0, 1.0, 1.0)))
            end
            gmsh.model.occ.synchronize()

            result = balanced_fuse!(tags)
            @test length(result) < length(tags)   # overlapping boxes must merge down
        end
    end

    @testset "balanced_fuse! fully converges a real interior block (docs/algorithm.md §11.2)" begin
        # A genuine 3x3x3-cell, all-INTERIOR block (81 struts, real pipeline
        # strut order and geometry) must fuse down to exactly 1 connected
        # solid -- every strut in the block shares a node with a neighbor,
        # so nothing should be left un-merged.
        lp = LatticeParams(5.0, 1.0)
        with_gmsh() do
            new_model("balanced-fuse-3x3x3")
            protos = build_prototypes(lp)
            n = 3
            tags = Int32[instantiate_strut(lp, protos, StrutRef(i,j,k,d))
                          for i in 0:n-1 for j in 0:n-1 for k in 0:n-1 for d in 1:3]
            remove_entities(collect(protos))
            @test length(tags) == 3*n^3

            result = balanced_fuse!(tags)
            @test length(result) == 1
        end
    end

    @testset "process_tile preserves volume against an untrimmed reference (docs/algorithm.md §6.3)" begin
        # Regression test for the COMMON-fragmentation bug (docs/algorithm.md
        # §6.3, §11.2): when the input body fully contains a tile's boundary
        # struts, trimming against it must not change the total volume
        # produced, compared to simply fusing the same struts with no
        # trimming at all.
        lp = LatticeParams(5.0, 1.0)
        td = TileDescriptor(TileKey(0,0,0), StrutRef[],
                             [StrutRef(0,0,0,1), StrutRef(0,0,0,2), StrutRef(0,0,0,3)])

        mktempdir() do dir
            body_brep = joinpath(dir, "input.brep")
            with_gmsh() do
                new_model("body")
                gmsh.model.occ.addBox(-20,-20,-20, 40,40,40)   # fully contains the tile
                write_model(body_brep)
            end

            out_path = joinpath(dir, "tile_0_0_0.brep")
            stat = process_tile(lp, body_brep, td, out_path, 2, (0,0,0))
            @test stat.solids >= 1
            @test stat.wrote

            expected_vol = with_gmsh() do
                new_model("expected")
                protos = build_prototypes(lp)
                tags = Int32[instantiate_strut(lp, protos, StrutRef(0,0,0,d)) for d in 1:3]
                f = fuse_all(tags)
                sum(volume_of(t) for t in f)
            end
            got_vol = with_gmsh() do
                new_model("got")
                vols = import_shapes(out_path)
                sum(volume_of(t) for (_,t) in vols)
            end
            @test got_vol ≈ expected_vol rtol=1e-6
        end
    end
    # filter_floating! decision-logic tests live in test_cleanup.jl, kept
    # separate since it's its own priority-#1-critical cleanup gate
    # (docs/algorithm.md §8, §11.2).

    @testset "assert_no_stray_solids (docs/algorithm.md §6.5, §8, §11.2)" begin
        # Regression coverage for the safeguard that replaced the implicit
        # guarantee assembly's old write-assembled.brep-then-reimport design
        # gave for free: with assembly and export now sharing one gmsh
        # session, nothing else forces the model to contain *only* the
        # tags the caller thinks it does.
        @testset "clean model: no strays, does not throw" begin
            with_gmsh() do
                new_model("stray-check-clean")
                a = Int32(gmsh.model.occ.addBox(0.0, 0.0, 0.0, 1.0, 1.0, 1.0))
                b = Int32(gmsh.model.occ.addBox(10.0, 0.0, 0.0, 1.0, 1.0, 1.0))
                gmsh.model.occ.synchronize()
                @test assert_no_stray_solids(Int32[a, b]) === nothing
            end
        end

        @testset "leaked solid not in final_tags: throws ProcessingError naming it" begin
            with_gmsh() do
                new_model("stray-check-leak")
                kept = Int32(gmsh.model.occ.addBox(0.0, 0.0, 0.0, 1.0, 1.0, 1.0))
                leaked = Int32(gmsh.model.occ.addBox(100.0, 0.0, 0.0, 1.0, 1.0, 1.0))
                gmsh.model.occ.synchronize()
                err = nothing
                try
                    assert_no_stray_solids(Int32[kept])
                catch e
                    err = e
                end
                @test err isa ProcessingError
                @test occursin(string(leaked), sprint(showerror, err))
            end
        end

        @testset "model_solids reflects exactly the model's dim-3 entities" begin
            with_gmsh() do
                new_model("model-solids")
                a = Int32(gmsh.model.occ.addBox(0.0, 0.0, 0.0, 1.0, 1.0, 1.0))
                b = Int32(gmsh.model.occ.addBox(10.0, 0.0, 0.0, 1.0, 1.0, 1.0))
                gmsh.model.occ.synchronize()
                @test sort(model_solids()) == sort(Int32[a, b])
            end
        end
    end
end
