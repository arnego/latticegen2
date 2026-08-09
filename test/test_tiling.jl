# Unit tests for src/tiling.jl — pure Julia, no gmsh required.

@testset "tiling.jl" begin
    @testset "tile_key uses floor division (handles negative indices)" begin
        n = 4
        @test tile_key(StrutRef(0, 0, 0, 1), n) == TileKey(0, 0, 0)
        @test tile_key(StrutRef(3, 3, 3, 1), n) == TileKey(0, 0, 0)
        @test tile_key(StrutRef(4, 0, 0, 1), n) == TileKey(1, 0, 0)
        @test tile_key(StrutRef(-1, 0, 0, 1), n) == TileKey(-1, 0, 0)   # not 0: fld(-1,4) == -1
        @test tile_key(StrutRef(-4, 0, 0, 1), n) == TileKey(-1, 0, 0)
        @test tile_key(StrutRef(-5, 0, 0, 1), n) == TileKey(-2, 0, 0)
    end

    @testset "tile_key origin offset avoids a spurious split at zero" begin
        # Regression test for a real bug found during development: bucketing
        # from absolute (0,0,0) always splits a range straddling zero into
        # >= 2 buckets per axis no matter how large n is (fld(-1,n) == -1
        # and fld(0,n) == 0 for every n >= 1) — for origin-centered input
        # geometry (the common case) this silently forced a minimum of 8
        # tiles regardless of --tile-cells, which turned out to be the root
        # cause of a persistent residual self-intersection in assembly
        # (docs/algorithm.md §11.1). With origin=(0,0,0) (the default), a
        # range from -6 to 6 needs 2 buckets even at n=100:
        n = 100
        @test tile_key(StrutRef(-6, 0, 0, 1), n) != tile_key(StrutRef(6, 0, 0, 1), n)
        # With origin set to the range's own minimum corner, one large-enough
        # tile covers the whole thing:
        origin = (-6, 0, 0)
        @test tile_key(StrutRef(-6, 0, 0, 1), n, origin) == tile_key(StrutRef(6, 0, 0, 1), n, origin)
    end

    @testset "partition_tiles basic bucketing" begin
        n = 2
        struts = [
            StrutRef(0,0,0,1), StrutRef(0,0,0,2), StrutRef(0,0,0,3),   # tile (0,0,0), interior
            StrutRef(1,1,1,1),                                         # tile (0,0,0), boundary
            StrutRef(2,0,0,1),                                         # tile (1,0,0), interior
            StrutRef(5,5,5,1),                                         # tile (2,2,2), outside -> dropped
        ]
        classes = [INTERIOR, INTERIOR, INTERIOR, BOUNDARY, INTERIOR, OUTSIDE]
        tiles = partition_tiles(struts, classes, n)

        @test length(tiles) == 2
        @test haskey(tiles, TileKey(0,0,0))
        @test haskey(tiles, TileKey(1,0,0))
        @test !haskey(tiles, TileKey(2,2,2))

        t00 = tiles[TileKey(0,0,0)]
        @test length(t00.interior) == 3
        @test length(t00.boundary) == 1
        @test StrutRef(1,1,1,1) in t00.boundary

        t10 = tiles[TileKey(1,0,0)]
        @test length(t10.interior) == 1
        @test isempty(t10.boundary)
    end

    @testset "partition_tiles rejects mismatched lengths / bad n" begin
        @test_throws ArgumentError partition_tiles([StrutRef(0,0,0,1)], StrutClass[], 2)
        @test_throws ArgumentError partition_tiles(StrutRef[], StrutClass[], 0)
    end

    @testset "is_full_interior" begin
        n = 2
        # A genuinely full n=2 tile has 3*2^3 = 24 interior struts spanning
        # cells i,j,k in {0,1}^3 and all 3 directions, no boundary struts.
        full_struts = StrutRef[]
        for i in 0:1, j in 0:1, k in 0:1, d in 1:3
            push!(full_struts, StrutRef(i,j,k,d))
        end
        full_classes = fill(INTERIOR, length(full_struts))
        tiles = partition_tiles(full_struts, full_classes, n)
        td = tiles[TileKey(0,0,0)]
        @test length(td.interior) == 3*n^3
        @test is_full_interior(td, n)

        @testset "one boundary strut disqualifies it" begin
            classes2 = copy(full_classes)
            classes2[1] = BOUNDARY
            tiles2 = partition_tiles(full_struts, classes2, n)
            td2 = tiles2[TileKey(0,0,0)]
            @test !is_full_interior(td2, n)
        end

        @testset "incomplete enumeration (fringe tile) disqualifies it" begin
            partial_struts = full_struts[1:end-1]   # missing one strut entirely
            partial_classes = fill(INTERIOR, length(partial_struts))
            tiles3 = partition_tiles(partial_struts, partial_classes, n)
            td3 = tiles3[TileKey(0,0,0)]
            @test isempty(td3.boundary)             # would wrongly look "full interior" by that check alone
            @test !is_full_interior(td3, n)          # but the count check catches it
        end
    end

    @testset "tile_translation matches node-position differences" begin
        lp = LatticeParams(5.0, 1.0)
        n = 3
        from = TileKey(0, 0, 0)
        to = TileKey(2, -1, 1)
        offset = tile_translation(lp, from, to, n)
        expected = node(lp, 2*n, -1*n, 1*n) - node(lp, 0, 0, 0)
        @test offset ≈ expected

        @test tile_translation(lp, from, from, n) ≈ Vec3(0.0, 0.0, 0.0)
    end

    @testset "interface_nodes (docs/algorithm.md §6.4a)" begin
        lp = LatticeParams(5.0, 1.0)

        @testset "count == shell of the (n+2)^3 index box" begin
            n = 4
            ifn = interface_nodes(lp, TileKey(0,0,0), n, (0,0,0))
            expected = (n+2)^3 - n^3   # full box minus its strictly-interior sub-box
            @test length(ifn) == expected
        end

        @testset "a strictly-interior node is absent" begin
            n = 4
            key = TileKey(0,0,0)
            origin = (0,0,0)
            ifn = interface_nodes(lp, key, n, origin)
            # Well inside [i0-1, i0+n] = [-1,4] on every axis.
            interior_pt = node(lp, 2, 2, 2)
            @test !any(p -> p ≈ interior_pt, ifn)
        end

        @testset "correct world positions via node(), and origin/key offset accounted for" begin
            n = 2
            key = TileKey(1, -1, 0)
            origin = (3, -2, 5)
            ifn = interface_nodes(lp, key, n, origin)
            i0, j0, k0 = origin[1] + key.bi*n, origin[2] + key.bj*n, origin[3] + key.bk*n
            # One of the 8 corners of the shell box must be present, at the
            # exact world position node() computes for it.
            corner = node(lp, i0-1, j0-1, k0-1)
            @test any(p -> p ≈ corner, ifn)
            # And it must equal the shell-box corner count formula too.
            @test length(ifn) == (n+2)^3 - n^3
        end

        @testset "shell count matches a face-by-face manual count for n=1" begin
            n = 1
            ifn = interface_nodes(lp, TileKey(0,0,0), n, (0,0,0))
            # Box is [-1,1]^3 (3x3x3 = 27 points), interior is [0,0]^3 (1 point).
            @test length(ifn) == 27 - 1
        end
    end
end
