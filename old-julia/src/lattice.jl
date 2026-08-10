# src/lattice.jl
#
# Lattice mathematics: strut directions, node/basis geometry, diamond profile
# frames, and candidate-strut enumeration from a world-space AABB. Implements
# docs/algorithm.md §2 and §3.1 exactly — do not approximate these formulas.

using LinearAlgebra

# --- Lightweight 3D vector type --------------------------------------------
#
# A small immutable, isbits struct (not a type-pirated NTuple/Base operator
# overload) used throughout lattice.jl, classify.jl, tiling.jl, and
# geomkernel.jl for cheap, allocation-free 3D vector arithmetic.

struct Vec3
    x::Float64
    y::Float64
    z::Float64
end

Vec3(t::NTuple{3,<:Real}) = Vec3(Float64(t[1]), Float64(t[2]), Float64(t[3]))

Base.:+(a::Vec3, b::Vec3) = Vec3(a.x + b.x, a.y + b.y, a.z + b.z)
Base.:-(a::Vec3, b::Vec3) = Vec3(a.x - b.x, a.y - b.y, a.z - b.z)
Base.:-(a::Vec3) = Vec3(-a.x, -a.y, -a.z)
Base.:*(s::Real, a::Vec3) = Vec3(s * a.x, s * a.y, s * a.z)
Base.:*(a::Vec3, s::Real) = s * a
Base.:/(a::Vec3, s::Real) = Vec3(a.x / s, a.y / s, a.z / s)
Base.isapprox(a::Vec3, b::Vec3; kwargs...) =
    isapprox(a.x, b.x; kwargs...) && isapprox(a.y, b.y; kwargs...) && isapprox(a.z, b.z; kwargs...)
Base.Tuple(a::Vec3) = (a.x, a.y, a.z)

dot3(a::Vec3, b::Vec3) = a.x * b.x + a.y * b.y + a.z * b.z
cross3(a::Vec3, b::Vec3) = Vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x)
norm3(a::Vec3) = sqrt(dot3(a, a))
normalize3(a::Vec3) = a / norm3(a)

# --- Strut recline angle (docs/algorithm.md §2.1) --------------------------

"""Exact strut recline angle from +Z, computed natively — not a decimal literal."""
const THETA = asin(sqrt(2 / 3))
const COS_THETA = inv(sqrt(3))
const SIN_THETA = sqrt(2 / 3)

"""
    strut_direction(k::Integer) -> Vec3

Unit direction vector for strut family `k ∈ {0,1,2}`, at azimuth `2πk/3`,
reclined from +Z by `THETA` (docs/algorithm.md §2.1).
"""
function strut_direction(k::Integer)
    φ = 2π * k / 3
    return Vec3(SIN_THETA * cos(φ), SIN_THETA * sin(φ), COS_THETA)
end

"""Cube edge length `a = cc / √2` (docs/algorithm.md §2.2)."""
cube_edge(cc::Real) = cc / sqrt(2)

# --- Precomputed lattice geometry ------------------------------------------

"""
    LatticeParams(cc, t)

All lattice geometry derived from the CLI parameters `cc` and `t`
(docs/algorithm.md §2, §3.1). Direction/profile-axis tuples are indexed
`1:3`, corresponding to strut families `e0, e1, e2` (`k = 0, 1, 2`) in the
math notation.
"""
struct LatticeParams
    cc::Float64
    t::Float64
    a::Float64                    # cube edge length
    directions::NTuple{3,Vec3}    # e0, e1, e2 (unit vectors)
    u::NTuple{3,Vec3}             # horizontal profile axis per direction
    v::NTuple{3,Vec3}             # vertical-plane profile axis per direction
    r::Float64                    # strut circumradius = t/√2
    B::NTuple{3,Vec3}             # basis columns: node(i,j,k) = i*B[1]+j*B[2]+k*B[3]
end

function LatticeParams(cc::Real, t::Real)
    a = cube_edge(cc)
    dirs = ntuple(m -> strut_direction(m - 1), 3)
    zhat = Vec3(0.0, 0.0, 1.0)
    us = ntuple(m -> normalize3(cross3(zhat, dirs[m])), 3)
    vs = ntuple(m -> cross3(dirs[m], us[m]), 3)
    r = t / sqrt(2)
    B = ntuple(m -> a * dirs[m], 3)
    return LatticeParams(Float64(cc), Float64(t), a, dirs, us, vs, r, B)
end

"""
    node(lp, i, j, k) -> Vec3

World-space position of lattice node `(i,j,k)` (docs/algorithm.md §2.2).
"""
node(lp::LatticeParams, i::Integer, j::Integer, k::Integer) =
    i * lp.B[1] + j * lp.B[2] + k * lp.B[3]

"""
    profile_vertices(lp, center, dir) -> NTuple{4,Vec3}

The four corner points of the diamond cross-section profile centered at
`center`, for strut direction index `dir ∈ 1:3` (docs/algorithm.md §3.1).
"""
function profile_vertices(lp::LatticeParams, center::Vec3, dir::Integer)
    hu = lp.r * lp.u[dir]
    hv = lp.r * lp.v[dir]
    return (center + hu, center + hv, center - hu, center - hv)
end

# --- Strut references and enumeration --------------------------------------

"""A strut identified by its origin node index and direction (`1:3`)."""
struct StrutRef
    i::Int
    j::Int
    k::Int
    dir::Int
end

"""World-space endpoints `(p0, p1)` of a strut, `p1 - p0 == lp.B[s.dir]`."""
function endpoints(lp::LatticeParams, s::StrutRef)
    p0 = node(lp, s.i, s.j, s.k)
    p1 = p0 + lp.B[s.dir]
    return p0, p1
end

"""World-space midpoint of a strut segment."""
function midpoint(lp::LatticeParams, s::StrutRef)
    p0, p1 = endpoints(lp, s)
    return 0.5 * (p0 + p1)
end

"""
    index_range(lp, lo, hi) -> (lo_idx, hi_idx)

Inclusive integer index range, as `NTuple{3,Int}` bounds, covering every
lattice node whose emanating struts could intersect the world-space AABB
`[lo, hi]`, padded by the strut circumradius (docs/algorithm.md §2.4).
"""
function index_range(lp::LatticeParams, lo::Vec3, hi::Vec3)
    Bmat = hcat(collect(Tuple(lp.B[1])), collect(Tuple(lp.B[2])), collect(Tuple(lp.B[3])))
    corners = (Vec3(x, y, z) for x in (lo.x, hi.x), y in (lo.y, hi.y), z in (lo.z, hi.z))
    idx_corners = [Bmat \ collect(Tuple(c)) for c in corners]
    mins = reduce((a, b) -> min.(a, b), idx_corners)
    maxs = reduce((a, b) -> max.(a, b), idx_corners)
    pad = ceil(Int, lp.r / lp.a) + 1
    lo_idx = floor.(Int, mins) .- pad
    hi_idx = ceil.(Int, maxs) .+ pad
    return (Tuple(lo_idx), Tuple(hi_idx))
end

# --- AABB overlap graph (docs/algorithm.md §6.3, §6.5) ---------------------
#
# Pure Julia, gmsh-free primitives shared by `balanced_fuse!`'s pre-filter
# (pipeline.jl), the boundary-COMMON operand-disjointness split
# (`trim_disjoint`, pipeline.jl), and the floating-body cleanup gate
# (`filter_floating!`, pipeline.jl). Kept here rather than in geomkernel.jl
# so they can be unit tested with no gmsh session at all.

"""
    aabbs_overlap(b1, b2; eps=1e-6) -> Bool

Do two axis-aligned bounding boxes, each `(xmin,ymin,zmin,xmax,ymax,zmax)`,
overlap (touching counts as overlapping, within `eps`)? No OCCT boolean math
— a cheap pre-filter over raw box coordinates."""
function aabbs_overlap(b1::NTuple{6,Float64}, b2::NTuple{6,Float64}; eps::Real=1e-6)
    return b1[1] <= b2[4] + eps && b1[4] >= b2[1] - eps &&
           b1[2] <= b2[5] + eps && b1[5] >= b2[2] - eps &&
           b1[3] <= b2[6] + eps && b1[6] >= b2[3] - eps
end

"""
    overlap_components(tags, boxes; eps=1e-6) -> Vector{Vector{Int32}}

Partition `tags` into connected components of the AABB-overlap graph (two
tags are in the same component iff there is a chain of pairwise
`aabbs_overlap` box overlaps connecting them), via union-find. `boxes` must
contain an entry for every tag in `tags` (see `bounding_boxes`).

This is the shared primitive behind the operand-disjointness invariant
(docs/algorithm.md §6.3): a singleton component is *provably* disjoint from
every other tag, so it is always safe to pass to a boolean call alone or
batched with other singletons; components with >1 member need to be
resolved (fused, or handled one at a time) before touching a boolean kernel
call, since OCCT's multi-operand booleans fragment (rather than union)
mutually-overlapping operands.

Naive O(n²) pairwise comparison, which is fine below a few thousand tags —
this is the common case (one tile or one assembly batch). Callers with
larger `tags` should pre-bucket by AABB center into a uniform grid sized at
the median box diagonal and only run this within/between neighboring
buckets, but no call site in this codebase currently needs that."""
function overlap_components(tags::Vector{Int32}, boxes::Dict{Int32,NTuple{6,Float64}}; eps::Real=1e-6)
    n = length(tags)
    parent = Dict{Int32,Int32}(t => t for t in tags)
    function find(x::Int32)
        while parent[x] != x
            parent[x] = parent[parent[x]]   # path halving
            x = parent[x]
        end
        return x
    end
    function union!(a::Int32, b::Int32)
        ra, rb = find(a), find(b)
        ra != rb && (parent[ra] = rb)
        return nothing
    end
    for i in 1:n, j in i+1:n
        ti, tj = tags[i], tags[j]
        aabbs_overlap(boxes[ti], boxes[tj]; eps=eps) && union!(ti, tj)
    end
    groups = Dict{Int32,Vector{Int32}}()
    for t in tags
        push!(get!(() -> Int32[], groups, find(t)), t)
    end
    return collect(values(groups))
end

"""
    candidate_struts(lp, lo, hi) -> Vector{StrutRef}

Every strut whose origin node index falls within the padded index range
covering the world-space AABB `[lo, hi]` (docs/algorithm.md §2.4). This may
include struts that are trivially outside the true volume — classification
(`classify.jl`) discards those cheaply.
"""
function candidate_struts(lp::LatticeParams, lo::Vec3, hi::Vec3)
    lo_idx, hi_idx = index_range(lp, lo, hi)
    n = (hi_idx[1] - lo_idx[1] + 1) * (hi_idx[2] - lo_idx[2] + 1) * (hi_idx[3] - lo_idx[3] + 1) * 3
    struts = Vector{StrutRef}(undef, 0)
    sizehint!(struts, max(n, 0))
    for i in lo_idx[1]:hi_idx[1], j in lo_idx[2]:hi_idx[2], k in lo_idx[3]:hi_idx[3]
        for d in 1:3
            push!(struts, StrutRef(i, j, k, d))
        end
    end
    return struts
end
