# Unit tests for src/lattice.jl — verifies every identity required by
# docs/algorithm.md §2.3. Pure Julia, no gmsh required.

@testset "lattice.jl" begin
    @testset "strut_direction identities (§2.3)" begin
        es = [strut_direction(k) for k in 0:2]

        @testset "unit length" begin
            for e in es
                @test norm3(e) ≈ 1.0 atol=1e-12
            end
        end

        @testset "recline angle == THETA" begin
            zhat = Vec3(0.0, 0.0, 1.0)
            for e in es
                @test acos(dot3(e, zhat)) ≈ THETA atol=1e-12
            end
            @test THETA ≈ deg2rad(54.7356) atol=1e-4
        end

        @testset "azimuthal separation is 120°" begin
            azimuth(e) = atan(e.y, e.x)
            φs = sort([mod(azimuth(e), 2π) for e in es])
            diffs = sort([mod(φs[2]-φs[1], 2π), mod(φs[3]-φs[2], 2π), mod(φs[1]-φs[3]+2π, 2π)])
            for d in diffs
                @test d ≈ 2π/3 atol=1e-10
            end
        end
    end

    @testset "node lattice identities (§2.3)" begin
        lp = LatticeParams(5.0, 1.0)  # cc=5, t=1

        @testset "in-plane neighbor spacing == cc, zero Δz" begin
            p10 = node(lp, 1, 0, 0)
            p01 = node(lp, 0, 1, 0)
            d = p10 - p01
            @test norm3(d) ≈ lp.cc atol=1e-10
            @test d.z ≈ 0.0 atol=1e-10
        end

        @testset "vertical space diagonal" begin
            p111 = node(lp, 1, 1, 1)
            p000 = node(lp, 0, 0, 0)
            diag = p111 - p000
            @test diag.x ≈ 0.0 atol=1e-10
            @test diag.y ≈ 0.0 atol=1e-10
            @test diag.z ≈ lp.a * sqrt(3) atol=1e-10
        end

        @testset "cube edge a = cc/√2" begin
            @test lp.a ≈ lp.cc / sqrt(2) atol=1e-12
        end

        @testset "node(i,j,k) == i*B1+j*B2+k*B3" begin
            @test node(lp, 2, -1, 3) ≈ (2*lp.B[1] + (-1)*lp.B[2] + 3*lp.B[3])
        end
    end

    @testset "profile frame orthonormality (§3.1)" begin
        lp = LatticeParams(5.0, 1.0)
        for d in 1:3
            e = lp.directions[d]
            u = lp.u[d]
            v = lp.v[d]
            @test norm3(u) ≈ 1.0 atol=1e-12
            @test norm3(v) ≈ 1.0 atol=1e-12
            @test dot3(u, e) ≈ 0.0 atol=1e-12
            @test dot3(v, e) ≈ 0.0 atol=1e-12
            @test dot3(u, v) ≈ 0.0 atol=1e-12
            @test u.z ≈ 0.0 atol=1e-12   # u is horizontal
        end
    end

    @testset "profile_vertices geometry" begin
        lp = LatticeParams(5.0, 2.0)  # t=2 -> diagonal = 2*sqrt2
        center = node(lp, 0, 0, 0)
        verts = profile_vertices(lp, center, 1)
        @test length(verts) == 4
        for vtx in verts
            @test norm3(vtx - center) ≈ lp.r atol=1e-12
        end
        # opposite corners: full diagonal length t*sqrt(2)
        @test norm3(verts[1] - verts[3]) ≈ lp.t * sqrt(2) atol=1e-10
        @test norm3(verts[2] - verts[4]) ≈ lp.t * sqrt(2) atol=1e-10
    end

    @testset "index_range / candidate_struts" begin
        lp = LatticeParams(5.0, 1.0)
        lo = Vec3(-10.0, -10.0, -10.0)
        hi = Vec3(10.0, 10.0, 10.0)
        lo_idx, hi_idx = index_range(lp, lo, hi)
        @test all(lo_idx .<= hi_idx)

        struts = candidate_struts(lp, lo, hi)
        @test length(struts) > 0
        @test all(s -> 1 <= s.dir <= 3, struts)
        @test all(s -> lo_idx[1] <= s.i <= hi_idx[1], struts)
        @test all(s -> lo_idx[2] <= s.j <= hi_idx[2], struts)
        @test all(s -> lo_idx[3] <= s.k <= hi_idx[3], struts)

        # A strut's endpoints must be exactly a*e_dir apart.
        s = struts[1]
        p0, p1 = endpoints(lp, s)
        @test norm3(p1 - p0) ≈ lp.a atol=1e-10
        @test midpoint(lp, s) ≈ 0.5*(p0+p1)

        # Every candidate node in a small, centered AABB must be covered:
        # in particular the origin node's three struts must be present.
        @test StrutRef(0,0,0,1) in struts
        @test StrutRef(0,0,0,2) in struts
        @test StrutRef(0,0,0,3) in struts
    end

    @testset "aabbs_overlap (docs/algorithm.md §6.3, §6.5)" begin
        unit = (0.0, 0.0, 0.0, 1.0, 1.0, 1.0)

        @testset "disjoint" begin
            far = (10.0, 10.0, 10.0, 11.0, 11.0, 11.0)
            @test !aabbs_overlap(unit, far)
            @test !aabbs_overlap(far, unit)   # symmetric
        end

        @testset "face-touching (exact shared face counts as overlap)" begin
            adjacent = (1.0, 0.0, 0.0, 2.0, 1.0, 1.0)   # shares the x=1 face exactly
            @test aabbs_overlap(unit, adjacent)
        end

        @testset "nested" begin
            inner = (0.25, 0.25, 0.25, 0.75, 0.75, 0.75)
            @test aabbs_overlap(unit, inner)
        end

        @testset "eps-separated: just inside eps counts as overlap, just outside does not" begin
            eps = 1e-6
            just_inside = (1.0 + eps/2, 0.0, 0.0, 2.0, 1.0, 1.0)
            just_outside = (1.0 + 10*eps, 0.0, 0.0, 2.0, 1.0, 1.0)
            @test aabbs_overlap(unit, just_inside; eps=eps)
            @test !aabbs_overlap(unit, just_outside; eps=eps)
        end
    end

    @testset "overlap_components (docs/algorithm.md §6.3, §6.5)" begin
        # Box helper: axis-aligned unit cube at integer offset (ox,oy,oz).
        boxat(ox, oy, oz) = (Float64(ox), Float64(oy), Float64(oz),
                              Float64(ox)+1.0, Float64(oy)+1.0, Float64(oz)+1.0)

        @testset "all singletons (mutually far apart)" begin
            tags = Int32[1, 2, 3]
            boxes = Dict{Int32,NTuple{6,Float64}}(1 => boxat(0,0,0), 2 => boxat(100,0,0), 3 => boxat(-100,0,0))
            comps = overlap_components(tags, boxes)
            @test length(comps) == 3
            @test all(c -> length(c) == 1, comps)
        end

        @testset "a chain (1-2 touch, 2-3 touch, 1-3 do not) merges into one component" begin
            tags = Int32[1, 2, 3]
            boxes = Dict{Int32,NTuple{6,Float64}}(1 => boxat(0,0,0), 2 => boxat(1,0,0), 3 => boxat(2,0,0))
            @test !aabbs_overlap(boxes[1], boxes[3])   # not directly touching
            comps = overlap_components(tags, boxes)
            @test length(comps) == 1
            @test sort(comps[1]) == Int32[1, 2, 3]
        end

        @testset "two separate clusters" begin
            tags = Int32[1, 2, 10, 11]
            boxes = Dict{Int32,NTuple{6,Float64}}(
                1 => boxat(0,0,0), 2 => boxat(1,0,0),        # cluster A: touching
                10 => boxat(1000,0,0), 11 => boxat(1001,0,0), # cluster B: touching, far from A
            )
            comps = overlap_components(tags, boxes)
            @test length(comps) == 2
            sizes = sort(length.(comps))
            @test sizes == [2, 2]
            clusterA = only(filter(c -> 1 in c, comps))
            @test sort(clusterA) == Int32[1, 2]
            clusterB = only(filter(c -> 10 in c, comps))
            @test sort(clusterB) == Int32[10, 11]
        end
    end
end
