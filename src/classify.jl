# src/classify.jl
#
# Boundary classification: surface tessellation, spatial hash, exact
# segment-triangle and point-triangle distance, and ray-cast point-in-solid
# testing. Implements docs/algorithm.md §5 exactly.
#
# Design note: only `tessellate_surface` touches the gmsh/OCCT kernel
# directly (per docs/algorithm.md §12's module map, tessellation is
# classify.jl's responsibility). Every other function here is pure Julia
# operating on plain `TriMesh` data, which is what makes it possible to unit
# test the classification algorithm against a synthetic analytic mesh (an
# icosphere, in test_classify.jl) with no gmsh dependency at all.

import Gmsh: gmsh

@enum StrutClass INTERIOR BOUNDARY OUTSIDE

# --- Triangle mesh and spatial hash -----------------------------------

"""A plain-data triangle surface mesh: vertex positions plus 1-based vertex
index triples per triangle."""
struct TriMesh
    verts::Vector{Vec3}
    tris::Vector{NTuple{3,Int}}
end

"""Target element size used as an approximation of chordal deviation
tolerance `d` (docs/algorithm.md §5.1): `d = min(t, a) / 10`."""
mesh_chordal_target(lp::LatticeParams) = min(lp.t, lp.a) / 10

"""
    check_surface_mesh_coverage(lp::LatticeParams; tol_factor::Real=4.0,
                                 min_tol::Real=0.05, max_tol::Real=3.0)

Defensive completeness gate on the 2D surface mesh gmsh just produced,
called from `tessellate_surface` immediately after `gmsh.model.mesh.generate(2)`
(docs/algorithm.md §5.1, §11.3). For every face currently in the model,
compares the OCCT/CAD-reported bounding box of that face
(`gmsh.model.getBoundingBox(2, tag)`, computed from the exact B-rep, not the
mesh) against the bounding box of the mesh nodes gmsh actually generated for
it. A healthy, CAD-conforming mesh places its boundary nodes essentially
exactly on the true face boundary, so any face whose meshed extent falls
short of its true extent by more than `tol` indicates the mesher silently
produced an incomplete or mis-parametrized triangulation for that face.

`tol = clamp(tol_factor * mesh_chordal_target(lp), min_tol, max_tol)` — a few
of the *finest* elements' worth of slack (`mesh_chordal_target`, the
curvature-refinement floor `d`, not the coarser `min(t, a)` cap most of a
mesh is actually sized at), **clamped to an absolute `[min_tol, max_tol]`
mm range independent of `t`/`cc`**. This clamp matters: an earlier version of
this check scaled its tolerance directly off `min(t, a)` with no ceiling, so
at larger `-cc`/`-t` the tolerance itself grew past the size of the very
defect this check exists to catch — e.g. at `cc=10, t=5` the old formula gave
a 20 mm tolerance against the ~18.7 mm real-world truncation this check was
built to catch (§11.3), silently defeating the guard at exactly the
parameter range a larger lattice would use. A completeness gate whose
sensitivity depends on unrelated CLI parameters rather than on the actual
mesh/geometry is not a real guarantee; the absolute ceiling keeps `tol`
comfortably below that defect's scale (and any question about the CAD kernel
producing something similarly-sized on other geometry) regardless of how
large `t`/`cc` are chosen. The floor (`min_tol`) exists purely so an
extremely fine lattice doesn't get an unrealistically-tight, floating-point-
noise-level tolerance.

This is a real, observed gmsh/OCCT failure mode, not a hypothetical one:
investigating a "large volumes of missing lattice" report against
`test/test-cylinder.STEP` (docs/algorithm.md §11.3) found that gmsh's
Frontal-Delaunay mesher recovered a truncated parametric domain for one
large, curved-boundary planar face of that part, silently producing a
mesh that covers only part of the face's true extent — and, worse, the
resulting mesh is still a perfectly well-formed *closed 2-manifold* (a
`manifold_check`-style edge-sharing test does not catch it at all; the
mesh is geometrically wrong, not topologically broken). Every downstream
classification decision (`classify_strut`, docs/algorithm.md §5.2) depends
entirely on this mesh faithfully representing the input solid's boundary, so
this failure mode silently misclassified every strut beyond the
mis-tessellated region as `OUTSIDE`, with no error or warning anywhere —
exactly the opposite of the "worst case is more work, never a wrong result"
guarantee docs/algorithm.md §10 promises for every other optimization in this
pipeline. Multiple mitigation attempts (curvature-independent uniform sizing,
three different `Mesh.Algorithm` choices, OCCT's import-time `OCCFix*`/
`OCCSewFaces` healing options, and `gmsh.model.occ.healShapes`) all
reproduced the identical truncation, so no gmsh-level workaround is known;
this check exists to convert the failure from silent data loss into a loud,
diagnosable one (`InputGeometryError`, exit 3) until/unless a real fix is
found (docs/algorithm.md §11.3).
"""
function check_surface_mesh_coverage(lp::LatticeParams; tol_factor::Real=4.0,
                                      min_tol::Real=0.05, max_tol::Real=3.0)
    tol = clamp(tol_factor * mesh_chordal_target(lp), min_tol, max_tol)
    for (d, tag) in gmsh.model.getEntities(2)
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, tag)
        _, _, elemNodeTags = gmsh.model.mesh.getElements(2, tag)
        nodetags = Int[]
        for arr in elemNodeTags
            append!(nodetags, arr)
        end
        isempty(nodetags) && throw(InputGeometryError(
            "Surface tessellation produced zero elements for face $tag (CAD bounding box " *
            "lo=($xmin,$ymin,$zmin) hi=($xmax,$ymax,$zmax)); the input geometry could not be " *
            "faithfully meshed for boundary classification (docs/algorithm.md §11.3)."))
        mlo = Vec3(Inf, Inf, Inf)
        mhi = Vec3(-Inf, -Inf, -Inf)
        for t in unique(nodetags)
            coord, _, _, _ = gmsh.model.mesh.getNode(t)
            mlo = Vec3(min(mlo.x, coord[1]), min(mlo.y, coord[2]), min(mlo.z, coord[3]))
            mhi = Vec3(max(mhi.x, coord[1]), max(mhi.y, coord[2]), max(mhi.z, coord[3]))
        end
        if mlo.x > xmin + tol || mlo.y > ymin + tol || mlo.z > zmin + tol ||
           mhi.x < xmax - tol || mhi.y < ymax - tol || mhi.z < zmax - tol
            throw(InputGeometryError(
                "Surface tessellation for face $tag does not cover its own geometry " *
                "(docs/algorithm.md §11.3): CAD bounding box lo=($xmin,$ymin,$zmin) " *
                "hi=($xmax,$ymax,$zmax) vs. meshed bounding box lo=($(mlo.x),$(mlo.y),$(mlo.z)) " *
                "hi=($(mhi.x),$(mhi.y),$(mhi.z)) (tolerance $(round(tol, digits=4)) mm). This " *
                "indicates gmsh's mesher silently produced an incomplete or mis-parametrized " *
                "triangulation for this face — classification against this mesh would silently " *
                "misclassify struts near the missing region as OUTSIDE rather than fail. No " *
                "gmsh-level mitigation is currently known (docs/algorithm.md §11.3); the input " *
                "face likely needs repair in the originating CAD tool."))
        end
    end
    return nothing
end

"""
    tessellate_surface(lp::LatticeParams) -> TriMesh

Mesh the 2D boundary of every volume currently in the gmsh model, and
extract it as a plain-data `TriMesh` (docs/algorithm.md §5.1). Requires an
active gmsh session with the input solid already imported and synchronized.

Sizing is curvature-adaptive rather than a single uniform target size:
element size is capped at `min(t, a)` (never coarser than a lattice cell, so
a whole strut could never hide inside one flat facet) with
`Mesh.MeshSizeFromCurvature` refining down to a floor of
`mesh_chordal_target(lp)` on tightly-curved features (small fillets etc.). A
naive uniform target of `mesh_chordal_target(lp)` everywhere was tried first
and produces a chordal deviation far smaller than needed on large
gently-curved surfaces (e.g. ~1.2M triangles / 47s on the 80mm test ball,
against ~12K triangles / 0.4s with curvature-adaptive sizing) — see
docs/algorithm.md §11 for why this matters for the optimization strategy.

After meshing, `check_surface_mesh_coverage` verifies every face's mesh
actually covers that face's true CAD extent before this function returns —
see its docstring for the real, observed gmsh meshing-robustness failure
this guards against (docs/algorithm.md §11.3): a mesher-recovered
parametric domain that silently omits part of a face without any error, and
without even breaking the resulting mesh's own 2-manifold closure.
"""
function tessellate_surface(lp::LatticeParams)
    gmsh.model.occ.synchronize()
    cap = min(lp.t, lp.a)
    floor_ = mesh_chordal_target(lp)
    gmsh.option.setNumber("Mesh.MeshSizeMax", cap)
    gmsh.option.setNumber("Mesh.MeshSizeMin", floor_)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 20)
    gmsh.model.mesh.generate(2)
    check_surface_mesh_coverage(lp)

    nodeTags, coords, _ = gmsh.model.mesh.getNodes()
    n = length(nodeTags)
    index_of = Dict{UInt64,Int}()
    sizehint!(index_of, n)
    verts = Vector{Vec3}(undef, n)
    for idx in 1:n
        tag = nodeTags[idx]
        index_of[tag] = idx
        verts[idx] = Vec3(coords[3idx-2], coords[3idx-1], coords[3idx])
    end

    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(2)
    tris = NTuple{3,Int}[]
    for ti in eachindex(elemTypes)
        elemTypes[ti] == 2 || continue   # gmsh element type 2 == 3-node triangle
        tags = elemNodeTags[ti]
        m = length(tags) ÷ 3
        sizehint!(tris, length(tris) + m)
        for e in 1:m
            # NOTE: `3e-2`/`3e-1` are NOT valid here — Julia parses those as
            # the float literals 0.03/0.3 (scientific notation), not
            # `3*e - 2`. Multiplication must be explicit.
            a = index_of[tags[3*e-2]]
            b = index_of[tags[3*e-1]]
            c = index_of[tags[3*e]]
            push!(tris, (a, b, c))
        end
    end
    isempty(tris) && throw(InputGeometryError("Surface tessellation produced zero triangles"))
    return TriMesh(verts, tris)
end

"""Axis-aligned bounding box `(lo, hi)` of triangle `ti` in `mesh`."""
function tri_aabb(mesh::TriMesh, ti::Int)
    a, b, c = mesh.tris[ti]
    pa, pb, pc = mesh.verts[a], mesh.verts[b], mesh.verts[c]
    lo = Vec3(min(pa.x, pb.x, pc.x), min(pa.y, pb.y, pc.y), min(pa.z, pb.z, pc.z))
    hi = Vec3(max(pa.x, pb.x, pc.x), max(pa.y, pb.y, pc.y), max(pa.z, pb.z, pc.z))
    return lo, hi
end

"""Uniform spatial hash over a `TriMesh`'s triangles, cell size ≈ 2× the
median triangle edge length (docs/algorithm.md §5.1)."""
struct SpatialHash
    cellsize::Float64
    cells::Dict{NTuple{3,Int},Vector{Int}}
end

_cellkey(cellsize::Real, p::Vec3) =
    (floor(Int, p.x / cellsize), floor(Int, p.y / cellsize), floor(Int, p.z / cellsize))

function build_spatial_hash(mesh::TriMesh)
    ntri = length(mesh.tris)
    edge_lengths = Vector{Float64}(undef, 3ntri)
    @inbounds for ti in 1:ntri
        a, b, c = mesh.tris[ti]
        pa, pb, pc = mesh.verts[a], mesh.verts[b], mesh.verts[c]
        edge_lengths[3ti-2] = norm3(pa - pb)
        edge_lengths[3ti-1] = norm3(pb - pc)
        edge_lengths[3ti]   = norm3(pc - pa)
    end
    median_edge = isempty(edge_lengths) ? 1.0 : sort(edge_lengths)[cld(length(edge_lengths), 2)]
    cellsize = max(2 * median_edge, 1e-9)

    cells = Dict{NTuple{3,Int},Vector{Int}}()
    for ti in 1:ntri
        lo, hi = tri_aabb(mesh, ti)
        klo, khi = _cellkey(cellsize, lo), _cellkey(cellsize, hi)
        for kx in klo[1]:khi[1], ky in klo[2]:khi[2], kz in klo[3]:khi[3]
            push!(get!(() -> Int[], cells, (kx, ky, kz)), ti)
        end
    end
    return SpatialHash(cellsize, cells)
end

"""Indices of triangles whose cells overlap the AABB `[lo, hi]` (de-duplicated)."""
function query_cells(sh::SpatialHash, lo::Vec3, hi::Vec3)
    klo, khi = _cellkey(sh.cellsize, lo), _cellkey(sh.cellsize, hi)
    result = Int[]
    for kx in klo[1]:khi[1], ky in klo[2]:khi[2], kz in klo[3]:khi[3]
        tris = get(sh.cells, (kx, ky, kz), nothing)
        tris === nothing || append!(result, tris)
    end
    return unique!(result)
end

# --- Exact distance primitives -----------------------------------------

"""Closest point on segment `[a,b]` to point `p`."""
function closest_point_segment(p::Vec3, a::Vec3, b::Vec3)
    ab = b - a
    denom = dot3(ab, ab)
    denom < 1e-300 && return a
    t = clamp(dot3(p - a, ab) / denom, 0.0, 1.0)
    return a + t * ab
end
dist_point_segment(p::Vec3, a::Vec3, b::Vec3) = norm3(p - closest_point_segment(p, a, b))

"""Closest point on triangle `(a,b,c)` to point `p` (Ericson, *Real-Time
Collision Detection* §5.1.5 — the standard robust closed-form algorithm)."""
function closest_point_triangle(p::Vec3, a::Vec3, b::Vec3, c::Vec3)
    ab = b - a
    ac = c - a
    ap = p - a
    d1 = dot3(ab, ap)
    d2 = dot3(ac, ap)
    (d1 <= 0 && d2 <= 0) && return a

    bp = p - b
    d3 = dot3(ab, bp)
    d4 = dot3(ac, bp)
    (d3 >= 0 && d4 <= d3) && return b

    vc = d1 * d4 - d3 * d2
    if vc <= 0 && d1 >= 0 && d3 <= 0
        v = d1 / (d1 - d3)
        return a + v * ab
    end

    cp = p - c
    d5 = dot3(ab, cp)
    d6 = dot3(ac, cp)
    (d6 >= 0 && d5 <= d6) && return c

    vb = d5 * d2 - d1 * d6
    if vb <= 0 && d2 >= 0 && d6 <= 0
        w = d2 / (d2 - d6)
        return a + w * ac
    end

    va = d3 * d6 - d5 * d4
    if va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        return b + w * (c - b)
    end

    denom = 1.0 / (va + vb + vc)
    v = vb * denom
    w = vc * denom
    return a + v * ab + w * ac
end
dist_point_triangle(p::Vec3, a::Vec3, b::Vec3, c::Vec3) = norm3(p - closest_point_triangle(p, a, b, c))

"""Minimum distance between segments `[p1,q1]` and `[p2,q2]` (Ericson, RTCD
§5.1.9 — the standard robust closed-form algorithm)."""
function dist_segment_segment(p1::Vec3, q1::Vec3, p2::Vec3, q2::Vec3)
    d1 = q1 - p1
    d2 = q2 - p2
    r = p1 - p2
    a = dot3(d1, d1)
    e = dot3(d2, d2)
    f = dot3(d2, r)
    EPS = 1e-12

    local s::Float64, t::Float64
    if a <= EPS && e <= EPS
        return norm3(p1 - p2)
    elseif a <= EPS
        s = 0.0
        t = clamp(f / e, 0.0, 1.0)
    else
        c = dot3(d1, r)
        if e <= EPS
            t = 0.0
            s = clamp(-c / a, 0.0, 1.0)
        else
            b = dot3(d1, d2)
            denom = a * e - b * b
            s = denom > EPS ? clamp((b * f - c * e) / denom, 0.0, 1.0) : 0.0
            t = (b * s + f) / e
            if t < 0.0
                t = 0.0
                s = clamp(-c / a, 0.0, 1.0)
            elseif t > 1.0
                t = 1.0
                s = clamp((b - c) / a, 0.0, 1.0)
            end
        end
    end
    c1 = p1 + s * d1
    c2 = p2 + t * d2
    return norm3(c1 - c2)
end

"""Does segment `[s0,s1]` intersect triangle `(a,b,c)`? Möller–Trumbore
adapted to a bounded segment (`t ∈ [0,1]`) and either triangle-facing
orientation (unlike `ray_triangle`, which is direction-aware for the
ray-cast parity test and only reports forward hits)."""
function segment_triangle_intersect(s0::Vec3, s1::Vec3, a::Vec3, b::Vec3, c::Vec3)
    EPS = 1e-9
    dir = s1 - s0
    e1 = b - a
    e2 = c - a
    h = cross3(dir, e2)
    det = dot3(e1, h)
    abs(det) < EPS && return false   # segment parallel to the triangle's plane
    invdet = 1.0 / det
    s = s0 - a
    u = invdet * dot3(s, h)
    (u < -EPS || u > 1 + EPS) && return false
    q = cross3(s, e1)
    v = invdet * dot3(dir, q)
    (v < -EPS || u + v > 1 + EPS) && return false
    t = invdet * dot3(e2, q)
    return t >= -EPS && t <= 1 + EPS
end

"""
    dist_segment_triangle(s0, s1, a, b, c) -> Float64

Exact minimum distance between segment `[s0,s1]` and triangle `(a,b,c)`.

If the segment actually pierces the triangle, the minimum distance is 0 —
checked explicitly first, since a proper intersection is *not* generally
found at a vertex/edge feature pair (the vertex-face / edge-edge
decomposition below is only valid for the disjoint case; two intersecting
convex shapes can have their closest common point in either's interior).

Otherwise, decomposed into the standard convex-feature-pair set for a
segment (2 vertices) vs. a triangle (3 vertices, 3 edges, 1 face): the two
point-to-triangle-face distances plus the three segment-to-edge distances
together cover every possible closest-pair configuration when disjoint,
including vertex-vertex and vertex-edge cases (subsumed by the endpoint
clamping inside the edge-edge computation).
"""
function dist_segment_triangle(s0::Vec3, s1::Vec3, a::Vec3, b::Vec3, c::Vec3)
    segment_triangle_intersect(s0, s1, a, b, c) && return 0.0
    d = dist_point_triangle(s0, a, b, c)
    d = min(d, dist_point_triangle(s1, a, b, c))
    d = min(d, dist_segment_segment(s0, s1, a, b))
    d = min(d, dist_segment_segment(s0, s1, b, c))
    d = min(d, dist_segment_segment(s0, s1, c, a))
    return d
end

"""
    segment_mesh_min_distance(sh, mesh, s0, s1, margin) -> Float64

Minimum distance from segment `[s0,s1]` to `mesh`, short-circuiting the
moment any triangle within `margin` is found (docs/algorithm.md §5.2) — the
caller only needs to know whether the true minimum is `≤ margin`, so an
early exit is both correct (a lower bound is never needed) and the main
performance win of this function.
"""
function segment_mesh_min_distance(sh::SpatialHash, mesh::TriMesh, s0::Vec3, s1::Vec3, margin::Real)
    lo = Vec3(min(s0.x, s1.x) - margin, min(s0.y, s1.y) - margin, min(s0.z, s1.z) - margin)
    hi = Vec3(max(s0.x, s1.x) + margin, max(s0.y, s1.y) + margin, max(s0.z, s1.z) + margin)
    candidates = query_cells(sh, lo, hi)
    best = Inf
    for ti in candidates
        a, b, c = mesh.tris[ti]
        d = dist_segment_triangle(s0, s1, mesh.verts[a], mesh.verts[b], mesh.verts[c])
        if d < best
            best = d
            best <= margin && return best
        end
    end
    return best
end

# --- Point-in-solid ray casting -----------------------------------------

"""Möller–Trumbore ray-triangle intersection: returns the ray parameter `t`
(`> 0`) if `origin + t*dir` hits triangle `(a,b,c)`, else `nothing`."""
function ray_triangle(origin::Vec3, dir::Vec3, a::Vec3, b::Vec3, c::Vec3)
    EPS = 1e-9
    e1 = b - a
    e2 = c - a
    h = cross3(dir, e2)
    det = dot3(e1, h)
    abs(det) < EPS && return nothing
    invdet = 1.0 / det
    s = origin - a
    u = invdet * dot3(s, h)
    (u < -EPS || u > 1 + EPS) && return nothing
    q = cross3(s, e1)
    v = invdet * dot3(dir, q)
    (v < -EPS || u + v > 1 + EPS) && return nothing
    t = invdet * dot3(e2, q)
    return t > EPS ? t : nothing
end

"""
    raycast_parity(sh, mesh, origin, dir, max_t) -> Bool

Ray-cast parity test along one ray: odd number of triangle hits means
`origin` is inside the (closed) surface described by `mesh`.

Accelerated via 3D DDA (Amanatides & Woo) traversal of the spatial hash
grid: only triangles in cells the ray actually passes through (out to
`max_t`) are ever tested, instead of the whole mesh. This is the dominant
cost in classification (point-in-solid testing runs for every candidate
strut that isn't near the surface — i.e. most of them, for any reasonably
large lattice) — brute-forcing every triangle per ray was measured at
~1.4 billion ray-triangle tests / ~26s classifying ~39K candidate struts
against a ~12K-triangle mesh; DDA traversal only visits the handful of
cells actually crossed. A triangle spanning multiple cells is only tested
once (via a per-call `visited` set), since it would otherwise be found
again in every one of those cells.
"""
function raycast_parity(sh::SpatialHash, mesh::TriMesh, origin::Vec3, dir::Vec3, max_t::Real)
    cs = sh.cellsize
    ix, iy, iz = _cellkey(cs, origin)
    stepx = dir.x > 0 ? 1 : (dir.x < 0 ? -1 : 0)
    stepy = dir.y > 0 ? 1 : (dir.y < 0 ? -1 : 0)
    stepz = dir.z > 0 ? 1 : (dir.z < 0 ? -1 : 0)

    _next_boundary(i, s, cs) = s > 0 ? (i + 1) * cs : i * cs

    tmaxx = dir.x != 0 ? (_next_boundary(ix, stepx, cs) - origin.x) / dir.x : Inf
    tmaxy = dir.y != 0 ? (_next_boundary(iy, stepy, cs) - origin.y) / dir.y : Inf
    tmaxz = dir.z != 0 ? (_next_boundary(iz, stepz, cs) - origin.z) / dir.z : Inf
    tdeltax = dir.x != 0 ? cs / abs(dir.x) : Inf
    tdeltay = dir.y != 0 ? cs / abs(dir.y) : Inf
    tdeltaz = dir.z != 0 ? cs / abs(dir.z) : Inf

    visited = Set{Int}()
    hits = 0
    t = 0.0
    while t <= max_t
        tris = get(sh.cells, (ix, iy, iz), nothing)
        if tris !== nothing
            for ti in tris
                ti in visited && continue
                push!(visited, ti)
                a, b, c = mesh.tris[ti]
                th = ray_triangle(origin, dir, mesh.verts[a], mesh.verts[b], mesh.verts[c])
                th === nothing || (hits += 1)
            end
        end
        (stepx == 0 && stepy == 0 && stepz == 0) && break
        if tmaxx < tmaxy
            if tmaxx < tmaxz
                ix += stepx; t = tmaxx; tmaxx += tdeltax
            else
                iz += stepz; t = tmaxz; tmaxz += tdeltaz
            end
        else
            if tmaxy < tmaxz
                iy += stepy; t = tmaxy; tmaxy += tdeltay
            else
                iz += stepz; t = tmaxz; tmaxz += tdeltaz
            end
        end
    end
    return isodd(hits)
end

# Three fixed, mutually non-parallel, axis/diagonal-avoiding ray directions.
# Deterministic (not sampled at runtime) so classification — and therefore
# the whole generation run — is exactly reproducible given the same inputs,
# consistent with the precision/correctness priority (specification.md
# "Key Considerations"). Chosen to avoid grazing common degenerate cases
# (axis-aligned or 45°-diagonal edges/vertices in typical CAD geometry).
const _RAY_DIRS = (
    normalize3(Vec3(0.8574, 0.2196, 0.4671)),
    normalize3(Vec3(-0.3157, 0.8842, 0.3427)),
    normalize3(Vec3(0.1732, -0.4899, 0.8547)),
)

"""Axis-aligned bounding box `(lo, hi)` of every vertex in `mesh`."""
function mesh_bounds(mesh::TriMesh)
    lo = mesh.verts[1]
    hi = mesh.verts[1]
    for v in mesh.verts
        lo = Vec3(min(lo.x, v.x), min(lo.y, v.y), min(lo.z, v.z))
        hi = Vec3(max(hi.x, v.x), max(hi.y, v.y), max(hi.z, v.z))
    end
    return lo, hi
end

"""
    point_inside(sh, mesh, p, max_t) -> Bool

Is `p` inside the closed surface `mesh`? Majority vote over 3 fixed ray
directions defeats the classic degenerate single-ray failure modes (a ray
that grazes a triangle edge or passes through a vertex) (docs/algorithm.md
§5.2). `max_t` bounds how far each ray is traced — any value at least as
large as the mesh's bounding diagonal from `p` is safe, since a closed
surface has zero net crossings beyond its own extent (see `mesh_bounds`).
"""
function point_inside(sh::SpatialHash, mesh::TriMesh, p::Vec3, max_t::Real)
    votes = 0
    for d in _RAY_DIRS
        votes += raycast_parity(sh, mesh, p, d, max_t) ? 1 : 0
    end
    return votes >= 2
end

# --- Top-level classification --------------------------------------------

"""
    classify_strut(lp, sh, mesh, d, s; max_t) -> StrutClass

Classify strut `s` relative to the input solid described by `mesh`
(docs/algorithm.md §5.2). `d` is the mesh chordal deviation tolerance
(`mesh_chordal_target(lp)`); the margin `r + d` ensures mesh approximation
error can never mis-promote a strut that truly touches the surface into
`INTERIOR` — ambiguity always degrades to the safe, more expensive
`BOUNDARY` path. `max_t` bounds the ray-cast distance for `point_inside`
(see its docstring) — compute once per mesh with `mesh_bounds` and pass it
through rather than recomputing per strut.
"""
function classify_strut(lp::LatticeParams, sh::SpatialHash, mesh::TriMesh, d::Real, s::StrutRef;
                         max_t::Real)
    s0, s1 = endpoints(lp, s)
    margin = lp.r + d
    dist = segment_mesh_min_distance(sh, mesh, s0, s1, margin)
    if dist > margin
        return point_inside(sh, mesh, midpoint(lp, s), max_t) ? INTERIOR : OUTSIDE
    else
        return BOUNDARY
    end
end
