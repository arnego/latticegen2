# tools/verify_geometry.jl
#
# Reusable geometry-validity checks for the E2E harness (specification.md
# §6.2, docs/testing.md): manifold/watertightness, self-intersection, and
# golden-sample volume-difference comparison. Not a `src/` module — these
# checks are QA tooling, not part of the generation pipeline itself, per
# specification.md §6.3 ("Test scripts and assets are separated into a
# separate `tools/` folder").
#
# Assumes `src/latticegen2.jl` has already been `include`d and `using
# .latticegen2` is in effect in the calling scope (see tools/e2e.jl).

"""
    manifold_check(mesh::TriMesh) -> (ok::Bool, bad_edges::Vector)

A closed, manifold triangle mesh has every edge shared by exactly 2
triangles. Returns `(true, [])` if so, else `(false, bad_edges)` where each
entry is `(edge_vertex_pair, incident_triangle_count)` for a failing edge
(specification.md §6.2 "valid closed manifold solid").
"""
function manifold_check(mesh::TriMesh)
    edge_count = Dict{Tuple{Int,Int},Int}()
    for (a, b, c) in mesh.tris
        for (u, v) in ((a, b), (b, c), (c, a))
            key = u < v ? (u, v) : (v, u)
            edge_count[key] = get(edge_count, key, 0) + 1
        end
    end
    bad = [(k, v) for (k, v) in edge_count if v != 2]
    return isempty(bad), bad
end

"""
    _straddles_plane(p1,p2,p3, plane_pt, plane_n; eps) -> Bool

True if the plane through `plane_pt` with normal `plane_n` genuinely
separates the three points — i.e. at least one has signed distance `> eps`
and at least one has signed distance `< -eps`. False both when all three
points are (within `eps`) on the plane (coplanar/touching) and when they're
all strictly on one side (no crossing possible).
"""
function _straddles_plane(p1::Vec3, p2::Vec3, p3::Vec3, plane_pt::Vec3, plane_n::Vec3; eps::Real)
    d1 = dot3(plane_n, p1 - plane_pt)
    d2 = dot3(plane_n, p2 - plane_pt)
    d3 = dot3(plane_n, p3 - plane_pt)
    has_pos = d1 > eps || d2 > eps || d3 > eps
    has_neg = d1 < -eps || d2 < -eps || d3 < -eps
    return has_pos && has_neg
end

"""
    triangles_properly_cross(a1,b1,c1, a2,b2,c2; eps=1e-7) -> Bool

True only if the two triangles genuinely, transversally cross in 3D — as
opposed to merely touching, being coplanar, or sharing a boundary. Requires
triangle 2's vertices to straddle triangle 1's plane *and* triangle 1's
vertices to straddle triangle 2's plane (`_straddles_plane`); only then is
the (already-validated) edge-piercing test `segment_triangle_intersect`
used to confirm and locate the crossing.

This two-stage check exists specifically to avoid a false-positive trap: two
independently-tessellated, genuinely separate solid bodies that merely
*touch* along a shared coincident face (perfectly valid geometry — e.g.
disconnected lattice islands, specification.md §1, or two adjacent tile
results before/without being fused into one topological solid) produce many
triangle pairs at that shared boundary whose edges cross in the mesh-only
sense (their independent triangulations aren't vertex-aligned) without any
volumetric overlap whatsoever. This was caught during development: the
single-triangle-pair edge-piercing test alone reported ~344 "intersections"
for two boxes sharing one exact face with zero volume overlap. The
plane-straddle pre-check correctly rejects that case (a touching/coplanar
configuration can never have both triangles straddle each other's plane).
"""
function triangles_properly_cross(a1::Vec3, b1::Vec3, c1::Vec3, a2::Vec3, b2::Vec3, c2::Vec3; eps::Real=1e-7)
    n1 = cross3(b1 - a1, c1 - a1)
    _straddles_plane(a2, b2, c2, a1, n1; eps=eps) || return false
    n2 = cross3(b2 - a2, c2 - a2)
    _straddles_plane(a1, b1, c1, a2, n2; eps=eps) || return false
    return segment_triangle_intersect(a1, b1, a2, b2, c2) ||
           segment_triangle_intersect(b1, c1, a2, b2, c2) ||
           segment_triangle_intersect(c1, a1, a2, b2, c2) ||
           segment_triangle_intersect(a2, b2, a1, b1, c1) ||
           segment_triangle_intersect(b2, c2, a1, b1, c1) ||
           segment_triangle_intersect(c2, a2, a1, b1, c1)
end

"""
    self_intersection_check(mesh::TriMesh) -> (ok::Bool, bad_pairs::Vector{Tuple{Int,Int}})

Detects genuinely crossing (not merely touching) triangle pairs: for each
pair of triangles that share no vertex (adjacent triangles are *expected*
to touch along a shared edge/vertex — that is not a self-intersection),
tests whether the two triangles properly cross in 3D (`triangles_properly_cross`).
Only candidate pairs whose AABBs are near each other (via the spatial hash)
are tested. Returns `(true, [])` if no crossings are found, else
`(false, pairs)` (specification.md §6.2 "no self-intersections").
"""
function self_intersection_check(mesh::TriMesh)
    sh = build_spatial_hash(mesh)
    bad = Tuple{Int,Int}[]
    for ti in eachindex(mesh.tris)
        a, b, c = mesh.tris[ti]
        lo, hi = tri_aabb(mesh, ti)
        for tj in query_cells(sh, lo, hi)
            tj <= ti && continue
            a2, b2, c2 = mesh.tris[tj]
            (a in (a2, b2, c2) || b in (a2, b2, c2) || c in (a2, b2, c2)) && continue
            pa, pb, pc = mesh.verts[a], mesh.verts[b], mesh.verts[c]
            pa2, pb2, pc2 = mesh.verts[a2], mesh.verts[b2], mesh.verts[c2]
            triangles_properly_cross(pa, pb, pc, pa2, pb2, pc2) && push!(bad, (ti, tj))
        end
    end
    return isempty(bad), bad
end

"""
    golden_sample_volume_diff(candidate_step_path, golden_step_path) -> Float64

Boolean-subtracts `candidate` and `golden` both ways (A-B and B-A) and
returns the larger of the two remainder volumes, in mm³ — near zero means
the two geometries occupy essentially the same volume
(specification.md §6.2). Only invoked when a golden sample is actually
configured (the `[TODO: for later]` scenarios in specification.md §6.1 are
not run automatically — see docs/testing.md).
"""
function golden_sample_volume_diff(candidate_step_path::AbstractString, golden_step_path::AbstractString)
    local diff_ab::Float64, diff_ba::Float64
    with_gmsh() do
        new_model("golden-diff-ab")
        cand = import_shapes(candidate_step_path)
        gold = import_shapes(golden_step_path)
        cand_tags = Int32[t for (_, t) in cand]
        gold_tags = Int32[t for (_, t) in gold]
        # Each cut consumes its operands, so subtract in one direction, then
        # reimport fresh copies of both files for the other direction rather
        # than trying to preserve/copy tags across two consuming operations.
        ab_out = cut_with(cand_tags, gold_tags)
        diff_ab = sum((volume_of(t) for t in ab_out); init=0.0)
    end
    with_gmsh() do
        new_model("golden-diff-ba")
        gold = import_shapes(golden_step_path)
        cand = import_shapes(candidate_step_path)
        gold_tags = Int32[t for (_, t) in gold]
        cand_tags = Int32[t for (_, t) in cand]
        ba_out = cut_with(gold_tags, cand_tags)
        diff_ba = sum((volume_of(t) for t in ba_out); init=0.0)
    end
    return max(diff_ab, diff_ba)
end
