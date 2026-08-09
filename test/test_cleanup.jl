# Unit tests for filter_floating! (src/pipeline.jl, docs/algorithm.md §8,
# §11.2) — the floating-body-only cleanup gate that replaced the earlier
# unconditional "volume < threshold -> delete" rule after that rule was
# found to delete connected junction material (docs/algorithm.md §6.3,
# §11.2's investigation of the test-cylinder-cc5t1 run). Kept as its own
# file, separate from test_pipeline.jl's broader orchestration tests, since
# this cleanup-correctness gate is priority #1 (precision) critical on its
# own (specification.md "Key Considerations").

import Gmsh: gmsh

@testset "cleanup (filter_floating!)" begin
    @testset "small disjoint solid is removed as a genuine floating body" begin
        with_gmsh() do
            new_model("filter-disjoint")
            big = Int32(gmsh.model.occ.addBox(0,0,0, 10,10,10))          # vol=1000
            small = Int32(gmsh.model.occ.addBox(100,0,0, 0.5,0.5,0.5))   # vol=0.125, far away
            gmsh.model.occ.synchronize()
            kept, removed, vols = filter_floating!(Int32[big, small], 1.0)
            @test removed == 1
            @test length(kept) == 1 && big in kept
            @test vols ≈ [0.125] atol=1e-9
        end
    end

    @testset "small solid touching a larger one is kept, not removed" begin
        with_gmsh() do
            new_model("filter-touching")
            big = Int32(gmsh.model.occ.addBox(0,0,0, 10,10,10))
            small = Int32(gmsh.model.occ.addBox(10,0,0, 0.5,0.5,0.5))    # shares the x=10 face
            gmsh.model.occ.synchronize()
            kept, removed, vols = filter_floating!(Int32[big, small], 1.0)
            @test removed == 0
            @test length(kept) == 1   # merged into one solid by the resolving fuse
            @test isempty(vols)
        end
    end

    @testset "a small solid sharing an overlap component with a kept large one is never counted as removed" begin
        # Multiple small solids, all touching one large solid transitively
        # (a chain), none of them individually large enough to survive the
        # threshold on their own -- but every one of them is genuinely
        # connected, so none may be deleted.
        with_gmsh() do
            new_model("filter-chain")
            big = Int32(gmsh.model.occ.addBox(0,0,0, 10,10,10))
            mid = Int32(gmsh.model.occ.addBox(10,0,0, 0.5,0.5,0.5))       # touches big
            tip = Int32(gmsh.model.occ.addBox(10.5,0,0, 0.5,0.5,0.5))     # touches mid, not big directly
            gmsh.model.occ.synchronize()
            kept, removed, vols = filter_floating!(Int32[big, mid, tip], 1.0)
            @test removed == 0
            @test isempty(vols)
            @test !isempty(kept)
        end
    end

    @testset "unresolvable connected sub-threshold component hard-fails" begin
        # Inject a fuse_fn that always fails, simulating a genuine OCCT
        # boolean robustness failure on an AABB-ambiguous (overlapping)
        # component -- must not silently keep or delete the ambiguous
        # sub-threshold solid, must raise ProcessingError instead
        # (specification.md priority #1, docs/algorithm.md §8).
        with_gmsh() do
            new_model("filter-unresolvable")
            big = Int32(gmsh.model.occ.addBox(0,0,0, 10,10,10))
            small = Int32(gmsh.model.occ.addBox(10,0,0, 0.5,0.5,0.5))    # AABB-touches big
            gmsh.model.occ.synchronize()
            failing_fuse(_) = error("simulated OCCT boolean robustness failure")
            @test_throws ProcessingError filter_floating!(Int32[big, small], 1.0; fuse_fn=failing_fuse)
        end
    end

    @testset "a resolvable-but-still-small component's fused result is re-evaluated, not assumed floating" begin
        # Two small solids that DO genuinely touch (no larger solid
        # involved) fuse into one solid still below threshold -- it must be
        # removed (a genuine small floating body), not kept just because it
        # came from a >1-member component.
        with_gmsh() do
            new_model("filter-small-fused-still-floating")
            a = Int32(gmsh.model.occ.addBox(100,0,0, 0.3,0.3,0.3))         # vol=0.027
            b = Int32(gmsh.model.occ.addBox(100.3,0,0, 0.3,0.3,0.3))       # touches a, vol=0.027
            gmsh.model.occ.synchronize()
            kept, removed, vols = filter_floating!(Int32[a, b], 1.0)
            @test removed == 1
            @test isempty(kept)
            @test length(vols) == 1
        end
    end

    @testset "empty input" begin
        kept, removed, vols = filter_floating!(Int32[], 1.0)
        @test isempty(kept) && removed == 0 && isempty(vols)
    end

    @testset "everything sub-threshold and floating -> kept is empty (caller decides how to react)" begin
        with_gmsh() do
            new_model("filter-all-floating")
            small = Int32(gmsh.model.occ.addBox(0,0,0, 0.5,0.5,0.5))
            gmsh.model.occ.synchronize()
            kept, removed, vols = filter_floating!(Int32[small], 1.0)
            @test isempty(kept)
            @test removed == 1
        end
    end
end
