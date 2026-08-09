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

"""Fetch every mesh node in the model once, returning `(coords, index)` where
`coords` is gmsh's flat `[x1,y1,z1,x2,...]` array and `index` maps a node tag
to its 1-based slot in it. Bulk-fetching once and indexing is deliberate: the
previous bounding-box version of `check_surface_mesh_coverage` made one
`gmsh.model.mesh.getNode` API call per node, which on a fine mesh of a large
part is millions of individual round-trips."""
function _mesh_node_index()
    nodeTags, coords, _ = gmsh.model.mesh.getNodes()
    idx = Dict{UInt64,Int}()
    sizehint!(idx, length(nodeTags))
    for i in eachindex(nodeTags)
        idx[nodeTags[i]] = i
    end
    return coords, idx
end

"""
    meshed_face_areas([coords, idx]) -> Dict{Int,Float64}

Total triangle area of the 2D mesh gmsh currently holds for each face in the
model, keyed by face tag. Pass an existing `(coords, idx)` pair from
`_mesh_node_index()` to avoid re-fetching the model's whole node array.
"""
meshed_face_areas() = meshed_face_areas(_mesh_node_index()...)

function meshed_face_areas(coords::Vector{Float64}, idx::Dict{UInt64,Int})
    areas = Dict{Int,Float64}()
    for (_, tag) in gmsh.model.getEntities(2)
        elemTypes, _, elemNodeTags = gmsh.model.mesh.getElements(2, tag)
        total = 0.0
        for ti in eachindex(elemTypes)
            elemTypes[ti] == 2 || continue   # gmsh element type 2 == 3-node triangle
            tags = elemNodeTags[ti]
            for e in 1:(length(tags) ÷ 3)
                # NOTE: `3e-2`/`3e-1` are NOT valid here — Julia parses those
                # as float literals (0.03/0.3), not `3*e - 2`.
                i1, i2, i3 = idx[tags[3*e-2]], idx[tags[3*e-1]], idx[tags[3*e]]
                ux = coords[3i2-2] - coords[3i1-2]
                uy = coords[3i2-1] - coords[3i1-1]
                uz = coords[3i2]   - coords[3i1]
                vx = coords[3i3-2] - coords[3i1-2]
                vy = coords[3i3-1] - coords[3i1-1]
                vz = coords[3i3]   - coords[3i1]
                cx = uy * vz - uz * vy
                cy = uz * vx - ux * vz
                cz = ux * vy - uy * vx
                total += 0.5 * sqrt(cx * cx + cy * cy + cz * cz)
            end
        end
        areas[Int(tag)] = total
    end
    return areas
end

"""
    check_surface_mesh_coverage(lp::LatticeParams; area_rel_tol::Real=0.25,
                                 bbox_slack::Real=1.0)

Defensive completeness gate on the 2D surface mesh gmsh just produced, called
from `tessellate_surface` immediately after `gmsh.model.mesh.generate(2)`
(docs/algorithm.md §5.1, §11.3). Every classification decision
(`classify_strut`, docs/algorithm.md §5.2) depends entirely on this mesh
faithfully representing the input solid's boundary, so a silently incomplete
mesh would misclassify every strut beyond the missing region as `OUTSIDE` —
exactly the opposite of the "worst case is more work, never a wrong result"
guarantee docs/algorithm.md §10 makes for the rest of the pipeline. This
check converts that failure mode into a loud, diagnosable one
(`InputGeometryError`, exit 3).

Two independent per-face tests, each using a quantity that is *sound in the
direction it is used in*:

1. **Coverage — exact trimmed area.** `gmsh.model.occ.getMass(2, tag)` is
   OCCT's exact Gauss-quadrature area of the *trimmed* face; the meshed area
   is the sum of that face's triangle areas. A face is rejected only when its
   shortfall is **both** more than `area_rel_tol` of its own area **and** more
   than `min_deficit` mm² in absolute terms. A mesher that recovered a
   truncated parametric domain loses a large, contiguous chunk of the face,
   which clears both bars easily.
2. **Containment — CAD bounding box, upper bound only.** Every mesh node of a
   face must lie *inside* that face's `gmsh.model.getBoundingBox(2, tag)`,
   inflated by `bbox_slack` mm. This catches a mesh that covers roughly the
   right *amount* of area but in the wrong *place*, which test 1 alone cannot
   see.

**Why the bounding box is only ever used as an upper bound.** This check
previously worked the other way around — it required the *meshed* bbox to
reach the CAD bbox — and that was a false-positive bug that blocked
`test/test-cylinder.STEP` outright (docs/algorithm.md §11.3). OCC's
`getBoundingBox` is deliberately **conservative**: for a B-spline edge it
returns the hull of the control points, and for a planar face the rectangle
of the untrimmed UV parameter domain — neither is the true trimmed extent.
On `test-cylinder.STEP`, face 9 (a `Plane` with B-spline boundary curves) is
reported as reaching `x=190.25`, while the surface actually stops at
`x=171.58`; sampling its boundary B-spline directly confirms `x=171.58` is
the real maximum. The mesh was correct all along — per-face exact-area
comparison puts every face of that part, face 9 included, within 0.005% of
its true area. A conservative over-estimate is safe to test *containment*
against (nothing may lie outside it) but is meaningless to test *coverage*
against (the mesh is not required to reach it), which is the distinction the
two tests above respect.

**Why both a relative and an absolute bar.** `area_rel_tol` defaults to 25%,
calibrated against the largest legitimate deficit measured on *input* CAD —
`test-cylinder.STEP`, `80mm-test-ball.step` and `TD_HX_Indre_Volum.step` at
`(cc, t)` spanning the documented CLI range (specification.md §3) — which was
**4.45%**, on a 4 mm² face of the heat-exchanger part at the coarsest setting
(`cc=50, t=20`). That deficit is ordinary chordal-deviation noise: a curved
face's triangles are secants, so meshed area sits slightly below true area,
more so the coarser the mesh.

A ratio alone is **not** sufficient, because this function is also run on
generated lattice *output*, not just input CAD — `tools/e2e.jl` re-tessellates
the finished `.step` to run its manifold and self-intersection checks. A
lattice has thousands of tiny trimmed sliver faces at strut/boundary
junctions, and on a sliver the ratio is meaningless: one observed face of
**0.0403 mm²** meshed to 0.0302 mm², i.e. 75.0% — enough to trip a pure
25%-ratio test, while the shortfall it represents (0.01 mm²) is physically
irrelevant. Hence `min_deficit = min_deficit_factor * min(t, a)^2`, one
coarsest mesh element's worth of area (`min(t, a)` is the element-size cap
`tessellate_surface` sets): a shortfall smaller than a single element is
meshing noise by construction, never a lost region. A real truncation loses a
large contiguous chunk and clears both bars by orders of magnitude — the
~18.7 mm-scale defect this gate was originally built for would be ~1500 mm²
against a 2.25 mm² floor at `cc=10, t=1.5`.

A face with **zero** elements is always an error regardless of either
tolerance.

Note that `min_deficit` scales with `t`/`cc` while `area_rel_tol` does not,
and that asymmetry is deliberate. An earlier version of this gate scaled its
*whole* sensitivity off `min(t, a)` with no ceiling, so at `cc=10, t=5` it
computed a 20 mm tolerance against the ~18.7 mm defect it existed to catch.
The difference is that `min_deficit` is a floor on what counts as a *region*
(legitimately tied to the mesh's own element size, and bounded above by the
relative test), not the sole gate — a large fractional loss on a large face
still fails no matter how big `t` is.
"""
function check_surface_mesh_coverage(lp::LatticeParams; area_rel_tol::Real=0.25,
                                      min_deficit_factor::Real=1.0,
                                      bbox_slack::Real=1.0)
    min_deficit = min_deficit_factor * min(lp.t, lp.a)^2
    coords, node_index = _mesh_node_index()
    meshed = meshed_face_areas(coords, node_index)

    for (_, tag) in gmsh.model.getEntities(2)
        exact = gmsh.model.occ.getMass(2, tag)
        ma = get(meshed, Int(tag), 0.0)
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, tag)

        ma > 0 || throw(InputGeometryError(
            "Surface tessellation produced zero elements for face $tag (exact CAD area " *
            "$(round(exact, digits=4)) mm^2, CAD bounding box lo=($xmin,$ymin,$zmin) " *
            "hi=($xmax,$ymax,$zmax)); the input geometry could not be faithfully meshed " *
            "for boundary classification (docs/algorithm.md §11.3)."))

        # Reject only when the shortfall is BOTH a large fraction of the face
        # AND larger in absolute terms than one coarsest mesh element — see the
        # docstring for why a ratio alone false-positives on sliver faces.
        if exact > 0 && ma < (1 - area_rel_tol) * exact && (exact - ma) > min_deficit
            throw(InputGeometryError(
                "Surface tessellation for face $tag covers only " *
                "$(round(100 * ma / exact, digits=2))% of that face's exact area " *
                "(meshed $(round(ma, digits=4)) mm^2 vs. exact $(round(exact, digits=4)) mm^2, " *
                "shortfall $(round(exact - ma, digits=4)) mm^2; tolerances " *
                "$(round(100 * area_rel_tol, digits=1))% and $(round(min_deficit, digits=4)) mm^2 " *
                "— docs/algorithm.md §11.3). This indicates gmsh's mesher silently produced an " *
                "incomplete or mis-parametrized triangulation for this face — classification " *
                "against this mesh would misclassify struts near the missing region as OUTSIDE " *
                "rather than fail. The input face likely needs repair in the originating CAD tool."))
        end

        # Containment: the CAD bbox is a conservative over-estimate, so it is
        # only ever valid to require that the mesh stays *within* it.
        # `includeBoundary=true` so the face's own boundary nodes (owned by its
        # bounding edges/vertices, not by the face) are checked too.
        face_nodes, _, _ = gmsh.model.mesh.getNodes(2, tag, true)
        for nt in face_nodes
            i = node_index[nt]
            x, y, z = coords[3i-2], coords[3i-1], coords[3i]
            if x < xmin - bbox_slack || x > xmax + bbox_slack ||
               y < ymin - bbox_slack || y > ymax + bbox_slack ||
               z < zmin - bbox_slack || z > zmax + bbox_slack
                throw(InputGeometryError(
                    "Surface tessellation for face $tag placed a mesh node at " *
                    "($x, $y, $z), outside that face's own CAD bounding box " *
                    "lo=($xmin,$ymin,$zmin) hi=($xmax,$ymax,$zmax) (slack " *
                    "$(bbox_slack) mm — docs/algorithm.md §11.3). This indicates gmsh's " *
                    "mesher produced a mis-parametrized triangulation for this face; " *
                    "classification against this mesh would be unreliable."))
            end
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
actually covers that face's exact trimmed area, and stays within that face's
CAD bounding box, before this function returns — see its docstring for what
that guards against, and for why the coverage half of the test must use
exact area rather than the deliberately-conservative CAD bounding box
(docs/algorithm.md §11.3).
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
