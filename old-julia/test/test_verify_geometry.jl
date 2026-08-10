# Regression tests for tools/verify_geometry.jl's self_intersection_check.
#
# Not a src/ module (tools/ holds QA tooling per specification.md §6.3), but
# its correctness is load-bearing for every geometry-validity claim the E2E
# harness makes, so it gets its own gmsh-dependent testset here rather than
# being trusted un-tested. See docs/algorithm.md §11.1 for the false-positive
# bug this guards against: two independently-tessellated, genuinely SEPARATE
# solids that merely touch along a shared coincident face were originally
# reported as ~344 "self-intersecting" triangle pairs by a pure edge-piercing
# test, purely because their independent triangulations don't align
# vertex-for-vertex at the shared boundary. The fix (a plane-straddle
# pre-check, `triangles_properly_cross`) must reject that case while still
# catching genuine transversal crossings.

include(joinpath(@__DIR__, "..", "tools", "verify_geometry.jl"))
import Gmsh: gmsh

function _box_pair_mesh(dx::Real, dy::Real, dz::Real)
    with_gmsh() do
        new_model("verify-geometry-test")
        gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
        gmsh.model.occ.addBox(dx, dy, dz, 1, 1, 1)
        gmsh.model.occ.synchronize()
        gmsh.option.setNumber("Mesh.MeshSizeMax", 0.3)
        gmsh.model.mesh.generate(2)
        nodeTags, coords, _ = gmsh.model.mesh.getNodes()
        idx_of = Dict{UInt64,Int}()
        verts = Vector{Vec3}(undef, length(nodeTags))
        for i in eachindex(nodeTags)
            idx_of[nodeTags[i]] = i
            verts[i] = Vec3(coords[3i-2], coords[3i-1], coords[3i])
        end
        elemTypes, _, elemNodeTags = gmsh.model.mesh.getElements(2)
        tris = NTuple{3,Int}[]
        for ti in eachindex(elemTypes)
            elemTypes[ti] == 2 || continue
            tags = elemNodeTags[ti]
            for e in 1:(length(tags)÷3)
                push!(tris, (idx_of[tags[3*e-2]], idx_of[tags[3*e-1]], idx_of[tags[3*e]]))
            end
        end
        TriMesh(verts, tris)
    end
end

@testset "self_intersection_check" begin
    @testset "disjoint boxes: no false positive" begin
        mesh = _box_pair_mesh(100, 0, 0)
        ok, bad = self_intersection_check(mesh)
        @test ok
        @test isempty(bad)
    end

    @testset "boxes touching at a shared coincident face: no false positive" begin
        # Regression case for docs/algorithm.md §11.1: shares the x=1 face
        # exactly, zero volume overlap. Previously reported ~344 pairs.
        mesh = _box_pair_mesh(1, 0, 0)
        ok, bad = self_intersection_check(mesh)
        @test ok
        @test isempty(bad)
    end

    @testset "boxes offset on all three axes, genuinely overlapping: detected" begin
        # No face-plane is shared or parallel-coincident between the two
        # boxes, so the overlap boundary requires true transversal crossings.
        mesh = _box_pair_mesh(0.5, 0.5, 0.5)
        ok, bad = self_intersection_check(mesh)
        @test !ok
        @test !isempty(bad)
    end

    @testset "hand-built X-crossing triangle pair: detected" begin
        verts = [Vec3(0, 0, 0), Vec3(2, 0, 0), Vec3(1, 0, 2),
            Vec3(1, -1, 1), Vec3(1, 1, 1), Vec3(1, 0, -1)]
        tris = [(1, 2, 3), (4, 5, 6)]
        mesh = TriMesh(verts, tris)
        ok, bad = self_intersection_check(mesh)
        @test !ok
        @test bad == [(1, 2)]
    end
end
