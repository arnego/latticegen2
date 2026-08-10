# Unit/smoke tests for src/geomkernel.jl — requires the gmsh/OCCT kernel.

import Gmsh: gmsh

# @__DIR__ is this file's own directory (test/), independent of the current
# working directory the test runner happens to be invoked from — notably,
# `Pkg.test()` does NOT run with pwd() at the project root, unlike
# `julia --project=. test/runtests.jl`.
const _BALL_STEP = joinpath(@__DIR__, "80mm-test-ball.step")
const _HX_STEP = joinpath(@__DIR__, "TD_HX_Indre_Volum.step")

@testset "geomkernel.jl" begin
    @testset "capture_output (pipe-based, docs/algorithm.md §6.2)" begin
        # Regression coverage for the deadlock hazard documented on
        # capture_output itself: an OS pipe's kernel buffer is small (as
        # low as a few KB), so anything that writes more than that and is
        # only drained *after* it finishes would hang forever. This writes
        # well past any plausible pipe buffer size and must still return
        # promptly with the full text intact.
        @testset "large output does not deadlock and is captured intact" begin
            n = 20_000
            elapsed = @elapsed begin
                result, text = capture_output() do
                    for i in 1:n
                        println("line $i " * "x"^40)
                    end
                    return :done
                end
                @test result == :done
                @test occursin("line 1 ", text)
                @test occursin("line $n ", text)
                @test count(c -> c == '\n', text) == n
            end
            @test elapsed < 30.0   # generous; a real deadlock never returns at all
        end

        @testset "stdout is restored even when f() throws" begin
            old = stdout
            @test_throws ErrorException capture_output() do
                println("before throw")
                error("boom")
            end
            @test stdout === old
            # A normal call afterwards must still work — confirms the pipe
            # redirect machinery itself wasn't left in a broken state.
            result, text = capture_output(() -> 7)
            @test result == 7
        end
    end

    @testset "import test fixtures" begin
        with_gmsh() do
            new_model("ball")
            vols = import_shapes(_BALL_STEP)
            @test length(vols) == 1
            lo, hi = bounding_box()
            @test lo.x ≈ -40.0 atol=0.01
            @test hi.x ≈ 40.0 atol=0.01
        end
        with_gmsh() do
            new_model("hx")
            vols = import_shapes(_HX_STEP)
            @test length(vols) == 1
        end
    end

    @testset "import missing file throws InputGeometryError" begin
        with_gmsh() do
            new_model("missing")
            @test_throws InputGeometryError import_shapes("test/does_not_exist.step")
        end
    end

    @testset "prototype construction and volume" begin
        lp = LatticeParams(5.0, 1.0)
        with_gmsh() do
            new_model("proto")
            protos = build_prototypes(lp)
            @test length(protos) == 3
            @test length(unique(protos)) == 3
            expected_vol = lp.t^2 * lp.a   # square profile area t^2 * extrusion length a
            for p in protos
                @test volume_of(p) ≈ expected_vol rtol = 1e-6
            end
        end
    end

    @testset "instantiate + fuse struts sharing a node" begin
        lp = LatticeParams(5.0, 1.0)
        with_gmsh() do
            new_model("fuse")
            protos = build_prototypes(lp)
            s1 = StrutRef(0, 0, 0, 1)
            s2 = StrutRef(0, 0, 0, 2)
            t1 = instantiate_strut(lp, protos, s1)
            t2 = instantiate_strut(lp, protos, s2)
            v1 = volume_of(t1)
            v2 = volume_of(t2)
            fused = fuse_all([t1, t2])
            @test length(fused) == 1
            vfused = volume_of(fused[1])
            @test vfused < v1 + v2          # solids overlap near the shared node
            @test vfused > max(v1, v2)      # union is bigger than either alone
        end
    end

    @testset "common_with (boolean intersection)" begin
        lp = LatticeParams(5.0, 1.0)
        with_gmsh() do
            new_model("common")
            protos = build_prototypes(lp)
            s = StrutRef(0, 0, 0, 1)
            strut_tag = instantiate_strut(lp, protos, s)
            v_before = volume_of(strut_tag)   # capture before the boolean op consumes the tag
            box = gmsh.model.occ.addBox(-1, -1, -1, 2, 2, 2)  # small box around the node
            gmsh.model.occ.synchronize()
            trimmed = common_with([strut_tag], Int32(box))
            @test length(trimmed) >= 1
            for tt in trimmed
                @test volume_of(tt) < v_before
            end
        end
    end

    @testset "write/re-import round trip" begin
        lp = LatticeParams(5.0, 1.0)
        mktempdir() do dir
            brep = joinpath(dir, "proto.brep")
            with_gmsh() do
                new_model("export")
                protos = build_prototypes(lp)
                fuse_all(collect(protos))
                write_model(brep)
            end
            @test isfile(brep)
            with_gmsh() do
                new_model("reimport")
                vols = import_shapes(brep)
                @test length(vols) >= 1
            end
        end
    end

    @testset "copy_translate" begin
        lp = LatticeParams(5.0, 1.0)
        with_gmsh() do
            new_model("copytrans")
            protos = build_prototypes(lp)
            v0 = volume_of(protos[1])
            moved = copy_translate(protos[1], Vec3(100.0, 0.0, 0.0))
            @test volume_of(moved) ≈ v0 rtol = 1e-9
            lo, hi = bounding_box(3, moved)
            @test lo.x > 50.0
        end
    end

    @testset "common_with operand-disjointness hazard (docs/algorithm.md §6.3)" begin
        lp = LatticeParams(5.0, 1.0)

        @testset "overlapping operands fragment instead of union" begin
            # The exact scenario documented in common_with's docstring: 3
            # struts sharing one lattice node, intersected against a
            # containing box in a single call, fragment into more than 1
            # solid instead of producing the union's single trimmed result.
            with_gmsh() do
                new_model("common-overlap-hazard")
                protos = build_prototypes(lp)
                tags = Int32[instantiate_strut(lp, protos, StrutRef(0,0,0,d)) for d in 1:3]
                box = Int32(gmsh.model.occ.addBox(-20,-20,-20, 40,40,40))
                gmsh.model.occ.synchronize()
                out = common_with(tags, box)
                @test length(out) > 1   # documents the fragmentation hazard
            end
        end

        @testset "disjoint operands trim independently -> exactly one result per operand" begin
            with_gmsh() do
                new_model("common-disjoint")
                protos = build_prototypes(lp)
                # Two struts far apart (different, non-adjacent nodes) -> disjoint.
                t1 = instantiate_strut(lp, protos, StrutRef(0,0,0,1))
                t2 = instantiate_strut(lp, protos, StrutRef(100,0,0,1))
                box1 = Int32(gmsh.model.occ.addBox(-20,-20,-20, 40,40,40))
                gmsh.model.occ.synchronize()
                @test !aabbs_overlap(bounding_box_nosync(3, t1), bounding_box_nosync(3, t2))
                out = common_with(Int32[t1], box1; remove_tool=false)
                @test length(out) == 1
                remove_entities(Int32[box1])
            end
        end

        @testset "remove_tool=false keeps the tool alive for a second call" begin
            with_gmsh() do
                new_model("common-remove-tool-false")
                protos = build_prototypes(lp)
                t1 = instantiate_strut(lp, protos, StrutRef(0,0,0,1))
                t2 = instantiate_strut(lp, protos, StrutRef(0,0,0,2))
                box = Int32(gmsh.model.occ.addBox(-20,-20,-20, 40,40,40))
                gmsh.model.occ.synchronize()
                r1 = common_with(Int32[t1], box; remove_tool=false)
                r2 = common_with(Int32[t2], box; remove_tool=false)   # would throw if box was consumed
                @test length(r1) == 1 && length(r2) == 1
            end
        end
    end

    @testset "bounding_boxes matches bounding_box, with one shared sync" begin
        lp = LatticeParams(5.0, 1.0)
        with_gmsh() do
            new_model("bounding-boxes")
            protos = build_prototypes(lp)
            tags = Int32[instantiate_strut(lp, protos, StrutRef(i,0,0,1)) for i in 0:4]
            gmsh.model.occ.synchronize()
            individual = Dict(t => bounding_box(3, t) for t in tags)
            batched = bounding_boxes(tags)
            @test sort(collect(keys(batched))) == sort(tags)
            for t in tags
                lo, hi = individual[t]
                b = batched[t]
                @test [lo.x, lo.y, lo.z, hi.x, hi.y, hi.z] ≈ [b[1], b[2], b[3], b[4], b[5], b[6]]
            end
        end
    end

    @testset "remove_model_entities + write_model(sync=false) round trip (docs/algorithm.md §8)" begin
        lp = LatticeParams(5.0, 1.0)
        mktempdir() do dir
            path = joinpath(dir, "removed.brep")
            with_gmsh() do
                new_model("remove-model-entities")
                protos = build_prototypes(lp)
                keep = Int32[instantiate_strut(lp, protos, StrutRef(i,0,0,1)) for i in 0:2]
                drop = Int32[instantiate_strut(lp, protos, StrutRef(i,5,0,1)) for i in 0:1]
                remove_entities(collect(protos))   # avoid the prototype-leak confound in this test
                gmsh.model.occ.synchronize()
                remove_model_entities(drop)
                write_model(path; sync=false)
            end
            @test isfile(path)
            with_gmsh() do
                new_model("remove-model-entities-check")
                vols = import_shapes(path)
                @test length(vols) == 3   # exactly the 3 kept struts, none of the 2 dropped
                # Each strut's B-rep is 6 faces (a square-profile extrusion:
                # 4 sides + 2 caps) -> 3 * 6 = 18 surface entities.
                ents2 = gmsh.model.getEntities(2)
                @test length(ents2) == 18
            end
        end
    end
end
