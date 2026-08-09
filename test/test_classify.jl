# Unit tests for src/classify.jl.
#
# Most of this module is pure Julia acting on plain TriMesh data, so it is
# tested here against a synthetic, analytically-known icosphere mesh rather
# than gmsh output — this is the fast, gmsh-independent path. A small
# gmsh-dependent smoke test for `tessellate_surface` itself lives at the
# bottom, guarded so it only runs when this file is included from the
# gmsh-dependent testset.

import Gmsh: gmsh

# --- Synthetic icosphere mesh (test-only helper, not part of src/) --------

function _icosahedron()
    φ = (1 + sqrt(5)) / 2
    raw = [
        (-1, φ, 0), (1, φ, 0), (-1, -φ, 0), (1, -φ, 0),
        (0, -1, φ), (0, 1, φ), (0, -1, -φ), (0, 1, -φ),
        (φ, 0, -1), (φ, 0, 1), (-φ, 0, -1), (-φ, 0, 1),
    ]
    verts = [normalize3(Vec3(x, y, z)) for (x, y, z) in raw]
    faces = NTuple{3,Int}[
        (1,12,6),(1,6,2),(1,2,8),(1,8,11),(1,11,12),
        (2,6,10),(6,12,5),(12,11,3),(11,8,7),(8,2,9),
        (4,10,5),(4,5,3),(4,3,7),(4,7,9),(4,9,10),
        (5,10,6),(3,5,12),(7,3,11),(9,7,8),(10,9,2),
    ]
    return verts, faces
end

function _subdivide(verts::Vector{Vec3}, faces::Vector{NTuple{3,Int}})
    midcache = Dict{Tuple{Int,Int},Int}()
    newverts = copy(verts)
    function mididx(i, j)
        key = i < j ? (i, j) : (j, i)
        haskey(midcache, key) && return midcache[key]
        m = normalize3(0.5 * (newverts[i] + newverts[j]))
        push!(newverts, m)
        idx = length(newverts)
        midcache[key] = idx
        return idx
    end
    newfaces = NTuple{3,Int}[]
    for (a, b, c) in faces
        ab, bc, ca = mididx(a, b), mididx(b, c), mididx(c, a)
        push!(newfaces, (a, ab, ca))
        push!(newfaces, (b, bc, ab))
        push!(newfaces, (c, ca, bc))
        push!(newfaces, (ab, bc, ca))
    end
    return newverts, newfaces
end

function icosphere(radius::Real, center::Vec3; subdiv::Int=3)
    verts, faces = _icosahedron()
    for _ in 1:subdiv
        verts, faces = _subdivide(verts, faces)
    end
    return TriMesh([center + radius * v for v in verts], faces)
end

function _brute_force_min_distance(mesh::TriMesh, s0::Vec3, s1::Vec3)
    best = Inf
    for (a, b, c) in mesh.tris
        best = min(best, dist_segment_triangle(s0, s1, mesh.verts[a], mesh.verts[b], mesh.verts[c]))
    end
    return best
end

# --- Tests ------------------------------------------------------------

@testset "classify.jl" begin
    @testset "point/segment/triangle distance primitives" begin
        a, b, c = Vec3(0,0,0), Vec3(1,0,0), Vec3(0,1,0)

        @testset "closest_point_triangle" begin
            @test dist_point_triangle(Vec3(0.2,0.2,1.0), a, b, c) ≈ 1.0 atol=1e-10  # directly above interior
            @test dist_point_triangle(Vec3(-1,0,0), a, b, c) ≈ 1.0 atol=1e-10       # beyond vertex a
            @test dist_point_triangle(Vec3(0,0,0), a, b, c) ≈ 0.0 atol=1e-10        # on a vertex
        end

        @testset "dist_point_segment" begin
            @test dist_point_segment(Vec3(0.5,1.0,0.0), Vec3(0,0,0), Vec3(1,0,0)) ≈ 1.0 atol=1e-10
            @test dist_point_segment(Vec3(-1,0,0), Vec3(0,0,0), Vec3(1,0,0)) ≈ 1.0 atol=1e-10
        end

        @testset "dist_segment_segment" begin
            # two perpendicular unit segments offset by 1 in z, crossing in xy
            d = dist_segment_segment(Vec3(-1,0,1), Vec3(1,0,1), Vec3(0,-1,0), Vec3(0,1,0))
            @test d ≈ 1.0 atol=1e-10
            # coincident collinear overlapping segments -> 0
            @test dist_segment_segment(Vec3(0,0,0), Vec3(2,0,0), Vec3(1,0,0), Vec3(3,0,0)) ≈ 0.0 atol=1e-10
        end

        @testset "dist_segment_triangle: piercing segment is ~0" begin
            s0, s1 = Vec3(0.2,0.2,-1.0), Vec3(0.2,0.2,1.0)
            @test dist_segment_triangle(s0, s1, a, b, c) ≈ 0.0 atol=1e-10
        end

        @testset "dist_segment_triangle: far segment matches brute force" begin
            s0, s1 = Vec3(5,5,5), Vec3(6,5,5)
            d = dist_segment_triangle(s0, s1, a, b, c)
            @test d > 0
        end
    end

    @testset "mesh_bounds" begin
        mesh = icosphere(40.0, Vec3(5.0, -3.0, 2.0); subdiv=2)
        lo, hi = mesh_bounds(mesh)
        @test lo.x < 5.0 - 39.0 && hi.x > 5.0 + 39.0
        @test all(v.x >= lo.x && v.x <= hi.x for v in mesh.verts)
    end

    @testset "ray_triangle" begin
        a, b, c = Vec3(0,0,0), Vec3(1,0,0), Vec3(0,1,0)
        @test ray_triangle(Vec3(0.2,0.2,-1.0), Vec3(0,0,1.0), a, b, c) !== nothing
        @test ray_triangle(Vec3(5,5,-1.0), Vec3(0,0,1.0), a, b, c) === nothing
        @test ray_triangle(Vec3(0.2,0.2,1.0), Vec3(0,0,1.0), a, b, c) === nothing  # points away
    end

    @testset "icosphere point_inside classification" begin
        R = 40.0
        mesh = icosphere(R, Vec3(0.0,0.0,0.0); subdiv=3)
        sh = build_spatial_hash(mesh)
        max_t = 500.0   # generous, well beyond the R=40 sphere's extent
        @test point_inside(sh, mesh, Vec3(0,0,0), max_t) == true
        @test point_inside(sh, mesh, Vec3(30,0,0), max_t) == true     # well inside
        @test point_inside(sh, mesh, Vec3(0,0,-30), max_t) == true
        @test point_inside(sh, mesh, Vec3(60,0,0), max_t) == false    # well outside
        @test point_inside(sh, mesh, Vec3(0,0,100), max_t) == false
    end

    @testset "raycast_parity (DDA) matches brute-force ray-triangle scan" begin
        mesh = icosphere(40.0, Vec3(0.0,0.0,0.0); subdiv=3)
        sh = build_spatial_hash(mesh)
        max_t = 500.0
        function brute_force_parity(mesh, origin, dir)
            hits = 0
            for (a,b,c) in mesh.tris
                t = ray_triangle(origin, dir, mesh.verts[a], mesh.verts[b], mesh.verts[c])
                t === nothing || (hits += 1)
            end
            return isodd(hits)
        end
        test_dirs = (normalize3(Vec3(1.0,0.0,0.0)), normalize3(Vec3(0.3,0.7,0.2)), normalize3(Vec3(-0.4,0.1,0.9)))
        for origin in (Vec3(0,0,0), Vec3(30,0,0), Vec3(60,0,0), Vec3(0,0,-30))
            for dir in test_dirs
                @test raycast_parity(sh, mesh, origin, dir, max_t) == brute_force_parity(mesh, origin, dir)
            end
        end
    end

    @testset "spatial hash matches brute-force when early exit can't fire" begin
        # segment_mesh_min_distance is documented to return only *some* value
        # ≤ margin once it confirms the true minimum is within margin (an
        # intentional early-exit optimization) — it does not promise the
        # exact global minimum in that regime. To validate the accelerated
        # search's exactness, pick a margin below the true distance so the
        # early exit can never fire, forcing a full, comparable scan.
        mesh = icosphere(40.0, Vec3(0.0,0.0,0.0); subdiv=2)
        sh = build_spatial_hash(mesh)
        @test !isempty(sh.cells)

        for (s0, s1) in ((Vec3(0,0,35), Vec3(0,0,45)),      # crosses surface near pole
                         (Vec3(50,0,0), Vec3(60,0,0)),       # far outside
                         (Vec3(0,0,0), Vec3(5,0,0)))         # deep inside
            expected = _brute_force_min_distance(mesh, s0, s1)
            got = segment_mesh_min_distance(sh, mesh, s0, s1, max(expected - 1e-6, 0.0))
            @test got ≈ expected atol=1e-6
        end
    end

    @testset "segment_mesh_min_distance margin contract" begin
        # True distance from this segment to the R=40 sphere is exactly 10.0
        # (closest point of [50,60] along x to the surface is x=50).
        mesh = icosphere(40.0, Vec3(0.0,0.0,0.0); subdiv=2)
        sh = build_spatial_hash(mesh)
        s0, s1 = Vec3(50,0,0), Vec3(60,0,0)
        @test segment_mesh_min_distance(sh, mesh, s0, s1, 5.0) > 5.0    # margin below true distance
        @test segment_mesh_min_distance(sh, mesh, s0, s1, 20.0) <= 20.0 # margin above true distance
    end

    @testset "classify_strut end-to-end against icosphere" begin
        lp = LatticeParams(5.0, 1.0)   # a ≈ 3.5355, r ≈ 0.7071
        mesh = icosphere(40.0, Vec3(0.0,0.0,0.0); subdiv=3)
        sh = build_spatial_hash(mesh)
        d = mesh_chordal_target(lp)
        max_t = 500.0

        @test classify_strut(lp, sh, mesh, d, StrutRef(0,0,0,1); max_t=max_t) == INTERIOR
        @test classify_strut(lp, sh, mesh, d, StrutRef(30,0,0,1); max_t=max_t) == OUTSIDE
        # node(11,0,0) has |.| ≈ 38.9 (inside R=40), node(12,0,0) ≈ 42.4 (outside):
        # this strut's segment must cross the sphere surface.
        @test classify_strut(lp, sh, mesh, d, StrutRef(11,0,0,1); max_t=max_t) == BOUNDARY
    end
end

# --- gmsh-dependent smoke test for tessellate_surface ----------------------
# Only meaningful when gmsh is available; included from the gmsh-dependent
# testset in runtests.jl, after test_geomkernel.jl.

@testset "classify.jl (gmsh-dependent)" begin
    @testset "tessellate_surface on 80mm test ball" begin
        lp = LatticeParams(10.0, 2.0)
        with_gmsh() do
            new_model("tess")
            import_shapes(joinpath(@__DIR__, "80mm-test-ball.step"))
            mesh = tessellate_surface(lp)
            @test length(mesh.tris) > 0
            @test length(mesh.verts) > 0
            # Every triangle vertex should lie approximately on the sphere of
            # radius 40 centered at the origin (the known fixture geometry).
            for (a, b, c) in mesh.tris[1:min(20, end)]
                for vi in (a, b, c)
                    @test norm3(mesh.verts[vi]) ≈ 40.0 atol = 0.5
                end
            end
        end
    end

    @testset "check_surface_mesh_coverage" begin
        @testset "no false positive on a well-formed, fully-meshed box" begin
            lp = LatticeParams(10.0, 2.0)
            with_gmsh() do
                new_model("coverage-ok")
                gmsh.model.occ.addBox(0, 0, 0, 20, 20, 20)
                gmsh.model.occ.synchronize()
                # Must not throw: tessellate_surface calls this internally,
                # and a healthy mesh of a simple box must pass cleanly.
                mesh = tessellate_surface(lp)
                @test length(mesh.tris) > 0
            end
        end

        @testset "raises InputGeometryError when a face has no mesh at all" begin
            lp = LatticeParams(10.0, 2.0)
            with_gmsh() do
                new_model("coverage-bad")
                gmsh.model.occ.addBox(0, 0, 0, 20, 20, 20)
                gmsh.model.occ.synchronize()
                gmsh.option.setNumber("Mesh.MeshSizeMax", 2.0)
                gmsh.option.setNumber("Mesh.MeshSizeMin", 0.2)
                gmsh.model.mesh.generate(2)
                # Delete every mesh element belonging to one face: that face's
                # meshed area is then trivially zero against a nonzero exact
                # area -- the "zero elements" branch of the gate.
                gmsh.model.mesh.clear([(2, 1)])
                @test_throws InputGeometryError check_surface_mesh_coverage(lp)
            end
        end

        @testset "raises InputGeometryError when a face's meshed area falls short of its exact area" begin
            # The coverage branch, exercised against real geometry rather than
            # a doctored mesh: a coarsely-meshed sphere's secant triangles
            # genuinely under-cover its exact (curved) area, so a tolerance
            # tightened below that deficit must reject it. This is the same
            # comparison a truly truncated parametric domain would trip, just
            # produced by a mechanism that can be created on demand.
            lp = LatticeParams(10.0, 2.0)
            with_gmsh() do
                new_model("coverage-area-deficit")
                gmsh.model.occ.addSphere(0, 0, 0, 10)
                gmsh.model.occ.synchronize()
                gmsh.option.setNumber("Mesh.MeshSizeMax", 8.0)
                gmsh.option.setNumber("Mesh.MeshSizeMin", 8.0)
                gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
                gmsh.model.mesh.generate(2)
                exact = gmsh.model.occ.getMass(2, 1)
                meshed = meshed_face_areas()[1]
                deficit = 1 - meshed / exact
                # Sanity-check the fixture actually under-covers, so the
                # assertion below is testing the gate and not a no-op.
                @test deficit > 0.05
                @test_throws InputGeometryError check_surface_mesh_coverage(lp; area_rel_tol=0.01)
                # ...and passes once the tolerance covers the real deficit.
                @test check_surface_mesh_coverage(lp; area_rel_tol=0.5) === nothing
                # The absolute bar alone must also be able to spare it: the
                # same under-covering mesh passes when its shortfall is
                # declared smaller than one mesh element's worth of area.
                @test check_surface_mesh_coverage(lp; area_rel_tol=0.01,
                                                  min_deficit_factor=1e6) === nothing
            end
        end

        @testset "a sliver face's meshing noise is not mistaken for a lost region" begin
            # Regression for a false positive found by re-tessellating generated
            # lattice output (which tools/e2e.jl does for its manifold and
            # self-intersection checks). A lattice has thousands of tiny trimmed
            # sliver faces at strut junctions; on one real 0.0403 mm^2 face the
            # meshed area was 0.0302 mm^2 = 75.0%, which trips a pure
            # relative-tolerance test even though the 0.01 mm^2 shortfall is
            # physically irrelevant. The absolute `min_deficit` bar is what
            # keeps the gate meaningful on such faces (docs/algorithm.md §11.3).
            lp = LatticeParams(10.0, 1.5)
            with_gmsh() do
                new_model("coverage-sliver")
                # A long, extremely thin sliver: area is tiny, so even a
                # proportionally large meshing shortfall is negligible.
                gmsh.model.occ.addRectangle(0, 0, 0, 20.0, 0.002)
                gmsh.model.occ.synchronize()
                gmsh.option.setNumber("Mesh.MeshSizeMax", 1.5)
                gmsh.option.setNumber("Mesh.MeshSizeMin", 0.15)
                gmsh.model.mesh.generate(2)
                exact = gmsh.model.occ.getMass(2, 1)
                # Whatever the ratio comes out at, the absolute shortfall is far
                # below one coarsest element (min(t,a)^2 = 2.25 mm^2) and must
                # not be reported as an incomplete tessellation.
                @test exact < 1.0
                @test check_surface_mesh_coverage(lp) === nothing
            end
        end

        @testset "raises InputGeometryError when a mesh node lands outside its face's CAD box" begin
            lp = LatticeParams(10.0, 2.0)
            with_gmsh() do
                new_model("coverage-containment")
                gmsh.model.occ.addBox(0, 0, 0, 20, 20, 20)
                gmsh.model.occ.synchronize()
                gmsh.option.setNumber("Mesh.MeshSizeMax", 4.0)
                gmsh.option.setNumber("Mesh.MeshSizeMin", 4.0)
                gmsh.model.mesh.generate(2)
                # Displace one of face 1's nodes far outside the face's own CAD
                # bounding box: the containment half of the gate must catch a
                # mesh that sits in the wrong place even when its area is right.
                tags, _, _ = gmsh.model.mesh.getNodes(2, 1, true)
                gmsh.model.mesh.setNode(tags[1], [500.0, 500.0, 500.0], Float64[])
                @test_throws InputGeometryError check_surface_mesh_coverage(lp)
            end
        end

        @testset "regression: test-cylinder.STEP at cc=10/t=1.5 tessellates faithfully" begin
            # docs/algorithm.md §11.3. This exact input/parameter combination
            # was previously rejected outright (InputGeometryError, exit 3) by
            # a bounding-box-based version of this gate, blocking the
            # dense-lattice scenario (specification.md §6.1) entirely. The
            # mesh was correct all along: OCC's getBoundingBox is a deliberate
            # over-estimate (B-spline control-point hull for edges, untrimmed
            # UV rectangle for planar faces), so requiring the mesh to *reach*
            # it was never a sound test. Guards against reintroducing that
            # false positive, and pins the exact-area evidence that the mesh is
            # in fact complete.
            lp = LatticeParams(10.0, 1.5)
            with_gmsh() do
                new_model("coverage-test-cylinder")
                import_shapes(joinpath(@__DIR__, "test-cylinder.STEP"))
                mesh = tessellate_surface(lp)
                @test length(mesh.tris) > 0

                meshed = meshed_face_areas()
                for (_, tag) in gmsh.model.getEntities(2)
                    exact = gmsh.model.occ.getMass(2, tag)
                    @test isapprox(meshed[Int(tag)], exact; rtol=0.01)
                end
            end
        end

        @testset "OCC face bounding boxes over-estimate, so coverage cannot be tested against them" begin
            # Pins the root cause of the false positive above (docs/algorithm.md
            # §11.3): face 9 of test-cylinder.STEP is a Plane bounded by B-spline
            # curves. Its reported CAD box reaches x ~= 190.25, but sampling its
            # boundary B-spline (edge 9) directly shows the real surface stops at
            # x ~= 171.58. Any future "the mesh must reach the CAD box" check
            # would fail on this part again.
            with_gmsh() do
                new_model("bbox-overestimate")
                import_shapes(joinpath(@__DIR__, "test-cylinder.STEP"))
                sync_model()
                _, _, _, xmax, _, _ = gmsh.model.getBoundingBox(2, 9)
                lo, hi = gmsh.model.getParametrizationBounds(1, 9)
                sampled_xmax = maximum(gmsh.model.getValue(1, 9, [lo[1] + (hi[1] - lo[1]) * i / 500])[1]
                                       for i in 0:500)
                @test xmax > sampled_xmax + 10.0
                @test isapprox(sampled_xmax, 171.58; atol=0.1)
            end
        end
    end
end
