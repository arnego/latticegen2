# src/tiling.jl
#
# Tile partitioning of classified struts into n×n×n-cell blocks in
# lattice-index space, and full-interior-tile detection for the periodicity
# shortcut. Implements docs/algorithm.md §6.1, §6.4.

"""
    TileKey

Identifies a tile by its block coordinates in lattice-index space; each
block spans `n` lattice cells along each axis (docs/algorithm.md §6.1).
"""
struct TileKey
    bi::Int
    bj::Int
    bk::Int
end

"""
    tile_key(s::StrutRef, n::Integer, origin::NTuple{3,Int}=(0,0,0)) -> TileKey

The tile block that strut `s`'s origin node belongs to, for tile edge `n`
cells, bucketed relative to `origin` rather than absolute index `(0,0,0)`.

This matters more than it looks: `fld(i, n)` always places a boundary at
`i = 0` for *any* `n ≥ 1` (`fld(-1,n) = -1`, `fld(0,n) = 0`, unconditionally
— there is no tile edge large enough to put both sides of that boundary in
one bucket). Real input parts are routinely centered on their own STEP
origin, so their candidate index range straddles zero in all three axes —
bucketing from absolute zero then *always* creates a minimum of 8 tiles (a
2×2×2 split through the geometry's own center), no matter how large `n` is.
This was discovered via a live diagnostic during development: a
`--tile-cells` deliberately set larger than the entire candidate range still
produced 8 tiles, not 1, and is suspected to be the root cause of a
persistent residual self-intersection in multi-tile assembly
(docs/algorithm.md §11.1). Bucketing relative to the candidate range's own
minimum corner (passed as `origin`, computed once in `run_pipeline`) removes
this spurious origin-centered split — a large enough `n` then genuinely
produces a single tile.
"""
tile_key(s::StrutRef, n::Integer, origin::NTuple{3,Int}=(0, 0, 0)) =
    TileKey(fld(s.i - origin[1], n), fld(s.j - origin[2], n), fld(s.k - origin[3], n))

"""
    TileDescriptor

Per-tile working set: struts split by classification (docs/algorithm.md
§6.1, §6.3). `OUTSIDE` struts are never stored here — they need no further
geometry work at all.
"""
struct TileDescriptor
    key::TileKey
    interior::Vector{StrutRef}
    boundary::Vector{StrutRef}
end

TileDescriptor(key::TileKey) = TileDescriptor(key, StrutRef[], StrutRef[])

"""
    partition_tiles(struts, classes, n; origin=(0,0,0)) -> Dict{TileKey,TileDescriptor}

Partition classified struts into `n×n×n`-cell tiles in lattice-index space
(docs/algorithm.md §6.1). `struts[i]`'s classification is `classes[i]`;
`OUTSIDE` struts are dropped. `origin` should be the candidate index range's
own minimum corner (not left at the default `(0,0,0)`) — see `tile_key` for
why this matters for origin-centered input geometry.
"""
function partition_tiles(struts::Vector{StrutRef}, classes::Vector{StrutClass}, n::Integer;
                          origin::NTuple{3,Int}=(0, 0, 0))
    length(struts) == length(classes) ||
        throw(ArgumentError("struts and classes must have the same length ($(length(struts)) vs $(length(classes)))"))
    n >= 1 || throw(ArgumentError("tile edge n must be >= 1 (got $n)"))

    tiles = Dict{TileKey,TileDescriptor}()
    for idx in eachindex(struts)
        c = classes[idx]
        c == OUTSIDE && continue
        s = struts[idx]
        key = tile_key(s, n, origin)
        td = get!(() -> TileDescriptor(key), tiles, key)
        if c == INTERIOR
            push!(td.interior, s)
        else
            push!(td.boundary, s)
        end
    end
    return tiles
end

"""
    is_full_interior(td::TileDescriptor, n::Integer) -> Bool

True when tile `td` contains the complete set of `3n³` interior struts (`n`
cells along each axis, 3 strut directions per cell) and no boundary struts
— i.e. it is geometrically congruent, via translation only, to any other
tile satisfying this condition, making it eligible for the periodicity
shortcut (docs/algorithm.md §6.4).

A tile with zero boundary struts but an *incomplete* interior strut count
(e.g. at the fringe of the candidate-enumeration index range, where some of
its struts were never enumerated at all) does **not** qualify: reusing a
cached unit-tile solid for it would silently add or omit struts that don't
actually belong there. This check is what makes the shortcut safe.
"""
is_full_interior(td::TileDescriptor, n::Integer) = isempty(td.boundary) && length(td.interior) == 3 * n^3

"""
    interface_nodes(lp, key, n, origin) -> Vector{Vec3}

World-space positions of every lattice node on tile `key`'s own shell —
every node index position a *neighbouring* tile's own struts could also
reach (docs/algorithm.md §6.4a). A solid produced entirely within one
tile whose AABB comes nowhere near any of these positions cannot possibly
merge with anything from outside its own tile, and — if also sub-threshold
and a singleton in the tile's own overlap graph — is therefore a genuine
floating body, safe to drop before the tile's `.brep` is even written (see
`drop_tile_local_islands` in pipeline.jl).

Covers the index box `[i0-1, i0+n] × [j0-1, j0+n] × [k0-1, k0+n]` (`i0 =
origin[1] + key.bi*n`, similarly for j0/k0) — one cell of margin beyond the
tile's own `n`-cell reach on every side, so a node just outside the tile's
own strut enumeration is still checked rather than silently excluded —
keeping only the nodes on that box's own outer shell (at least one
coordinate at an extreme). The interior of that box is never a node any
other tile could reach, so it is skipped.
"""
function interface_nodes(lp::LatticeParams, key::TileKey, n::Integer, origin::NTuple{3,Int})
    i0 = origin[1] + key.bi * n
    j0 = origin[2] + key.bj * n
    k0 = origin[3] + key.bk * n
    ilo, ihi = i0 - 1, i0 + n
    jlo, jhi = j0 - 1, j0 + n
    klo, khi = k0 - 1, k0 + n
    nodes = Vec3[]
    for i in ilo:ihi, j in jlo:jhi, k in klo:khi
        if i == ilo || i == ihi || j == jlo || j == jhi || k == klo || k == khi
            push!(nodes, node(lp, i, j, k))
        end
    end
    return nodes
end

"""
    tile_translation(lp, from, to, n) -> Vec3

World-space translation vector from tile `from` to tile `to`, given tile
edge `n` cells. For two full-interior tiles, applying this to a cached
unit-tile solid (via `copy_translate`) reproduces exactly what re-fusing
`to`'s own struts from scratch would have produced (docs/algorithm.md §6.4).
"""
function tile_translation(lp::LatticeParams, from::TileKey, to::TileKey, n::Integer)
    di = (to.bi - from.bi) * n
    dj = (to.bj - from.bj) * n
    dk = (to.bk - from.bk) * n
    return di * lp.B[1] + dj * lp.B[2] + dk * lp.B[3]
end
