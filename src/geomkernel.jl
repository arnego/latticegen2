# src/geomkernel.jl
#
# All interaction with the Gmsh/OCCT geometry kernel lives here. Implements
# docs/algorithm.md §3.2 (prototype construction) and the low-level
# primitives used by classify.jl (mesh access), tiling.jl/pipeline.jl (fuse,
# common, copy/translate), and stepmeta.jl (export). Nothing outside this
# file should call into the `gmsh` module directly.

import Gmsh: gmsh

# --- Session scoping ---------------------------------------------------

"""
    with_gmsh(f::Function; terminal=false)

Run `f()` inside a scoped gmsh session (`initialize` -> `f()` -> `finalize`,
the latter always running, even on error). Each worker process / tile owns
exactly one such session at a time — gmsh/OCCT state is process-global, so
sessions are never nested or shared across threads (docs/algorithm.md §6.2).
"""
function with_gmsh(f::Function; terminal::Bool=false)
    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Terminal", terminal ? 1 : 0)
        gmsh.option.setNumber("General.Verbosity", terminal ? 5 : 0)
        try
            # Defensive unit pin (docs/algorithm.md §8). Wrapped in try/catch
            # because option availability can vary across gmsh/OCCT builds;
            # mm is already OCCT's default, so failure here is non-fatal.
            gmsh.option.setString("Geometry.OCCTargetUnit", "MM")
        catch
        end
        return f()
    finally
        gmsh.finalize()
    end
end

"""
    capture_output(f::Function) -> (result, text::String)

Run `f()`, redirecting process-level stdout into a string instead of letting
it leak to the console. This is necessary because OCCT's own C++ trace
output (e.g. STEP writer "Statistics on Transfer") bypasses gmsh's own
`General.Terminal`/`General.Verbosity` options entirely. The caller routes
`text` into the run log at file-only verbosity (docs/algorithm.md §9) rather
than discarding it, so it remains available for diagnostics.
"""
function capture_output(f::Function)
    path, io = mktemp()
    local result
    try
        redirect_stdout(io) do
            result = f()
        end
    finally
        close(io)
    end
    text = read(path, String)
    rm(path; force=true)
    return result, text
end

# --- Import / export -----------------------------------------------------

"""
    import_shapes(path) -> Vector{Tuple{Int32,Int32}}

Import a STEP or BREP file into the current gmsh model and synchronize.
Returns the dim-tags of exactly the volumes (dim 3) *this call* just
imported, taken directly from `gmsh.model.occ.importShapes`'s own return
value — deliberately **not** `gmsh.model.getEntities(3)`, which returns
every entity in the whole model and would silently accumulate stale
entities from any earlier import in the same gmsh session (as happens when
assembling many tile files into one session — docs/algorithm.md §6.5).
Throws `InputGeometryError` on any failure (docs/algorithm.md §4 pipeline
stage C).
"""
function import_shapes(path::AbstractString)
    isfile(path) || throw(InputGeometryError("Geometry file not found: $path"))
    local out
    try
        out, _ = capture_output(() -> gmsh.model.occ.importShapes(String(path)))
    catch e
        # A Ctrl+C landing mid-import is a cancellation, not a bad input file:
        # let it through rather than mislabelling it as an exit-3 geometry error.
        is_interrupt(e) && rethrow()
        throw(InputGeometryError("Failed to import geometry from $path: $(sprint(showerror, e))"))
    end
    gmsh.model.occ.synchronize()
    vols = filter(dt -> dt[1] == 3, out)
    isempty(vols) && throw(InputGeometryError("No solid volumes found in $path"))
    return vols
end

"""Thin wrapper around `gmsh.model.occ.synchronize()`, so call sites that
need to batch several kernel queries behind a single sync can say so
explicitly (see `bounding_boxes` below)."""
sync_model() = gmsh.model.occ.synchronize()

"""Bounding box `(lo, hi)` as `Vec3` of the given dim-tag, or of the whole
model when `dim = -1, tag = -1` (the default). Self-synchronizing: OCC
kernel calls (`copy`, `translate`, `fuse`, ...) do not update gmsh's own
model-level bookkeeping until synchronized, so this always syncs first
rather than relying on every caller to remember to.

**Cost note:** `synchronize()` is O(whole model), not O(this query) — a
202× overhead was measured querying 200 boxes one at a time via this
function versus one `sync_model()` followed by 200 raw `getBoundingBox`
calls (docs/algorithm.md §6.5). Any caller that needs more than one box
from the *same*, already-synchronized model state should use
`bounding_box_nosync`/`bounding_boxes` instead."""
function bounding_box(dim::Integer=-1, tag::Integer=-1)
    gmsh.model.occ.synchronize()
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
    return Vec3(xmin, ymin, zmin), Vec3(xmax, ymax, zmax)
end

"""Bounding box of dim-tag `(dim, tag)` as a raw `NTuple{6,Float64}`
`(xmin,ymin,zmin,xmax,ymax,zmax)`, with **no** synchronize call. Only safe to
call when the caller already knows the model is synchronized (typically via
a preceding `sync_model()` — see `bounding_boxes`)."""
bounding_box_nosync(dim::Integer, tag::Integer) = gmsh.model.getBoundingBox(dim, tag)

"""
    bounding_boxes(tags::Vector{Int32}) -> Dict{Int32,NTuple{6,Float64}}

Bounding boxes of every dim-3 volume in `tags`, with exactly **one**
`sync_model()` call regardless of how many tags are given (docs/algorithm.md
§6.5) — the fix for the 202× `bounding_box`-per-call overhead. Use this (or
`bounding_box_nosync` after a manual `sync_model()`) instead of calling
`bounding_box` in a loop."""
function bounding_boxes(tags::Vector{Int32})
    sync_model()
    return Dict{Int32,NTuple{6,Float64}}(t => bounding_box_nosync(3, t) for t in tags)
end

"""
    write_model(path; sync=true) -> String

Export the current model to `path` (format inferred from extension — `.brep`
for fast exact intermediate staging per docs/algorithm.md §6.2, `.step` for
the final AP214 output per §8). Synchronizes first by default (see
`bounding_box`), since `gmsh.write` walks gmsh's model-level entity list, not
the raw OCC kernel state. Returns any captured OCCT trace text.

`sync=false` skips that synchronize call. This is **only** safe immediately
after `remove_model_entities`, whose whole point is to edit gmsh's
model-level entity list directly without touching the OCC kernel
(docs/algorithm.md §8) — a `synchronize()` in between would walk the *still
unmodified* OCC kernel state and resurrect the entities `remove_model_entities`
just removed from the model list. Do not pass `sync=false` anywhere else.
"""
function write_model(path::AbstractString; sync::Bool=true)
    sync && gmsh.model.occ.synchronize()
    _, text = capture_output(() -> gmsh.write(String(path)))
    return text
end

"""
    remove_model_entities(tags::Vector{Int32}) -> Nothing

Remove dim-3 entities `tags` from gmsh's **model-level** entity list via
`gmsh.model.removeEntities` (recursive: their bounding lower-dim entities go
too), *without* touching the OCC kernel's own data. Measured 16× faster than
`remove_entities`'s `gmsh.model.occ.remove` for a 300-solid removal from a
648-solid model (11.97 s -> 0.74 s, docs/algorithm.md §8) — the OCC-level
path recomputes kernel bookkeeping that model-level export doesn't need.

**Caveat, and the reason `remove_entities` still exists:** because this
never touches the OCC kernel, a subsequent `gmsh.model.occ.synchronize()`
would re-populate the model list from the (unmodified) kernel state and
silently undo the removal. Only ever follow this with `write_model(path;
sync=false)` — never with any function that synchronizes. If the model needs
further OCC-level boolean/copy/translate work after removal, use
`remove_entities` instead, which removes the entities from the kernel itself.
"""
function remove_model_entities(tags::Vector{Int32})
    isempty(tags) && return nothing
    gmsh.model.removeEntities([(3, t) for t in tags], true)
    return nothing
end

"""Start a fresh named model in the current gmsh session (used for the STEP
part name — docs/algorithm.md §8)."""
new_model(name::AbstractString) = gmsh.model.add(String(name))

# --- Strut prototype construction (docs/algorithm.md §3.2) ---------------

"""
    build_prototype(lp::LatticeParams, dir::Integer) -> Int32

Construct the direction-`dir` strut prototype solid at the origin: diamond
profile wire -> planar face -> extrude by `a * e_dir`. Returns the resulting
volume tag. `dir ∈ 1:3`.
"""
function build_prototype(lp::LatticeParams, dir::Integer)
    center = Vec3(0.0, 0.0, 0.0)
    verts = profile_vertices(lp, center, dir)
    pts = Int32[gmsh.model.occ.addPoint(v.x, v.y, v.z) for v in verts]
    lines = Int32[gmsh.model.occ.addLine(pts[m], pts[mod1(m + 1, 4)]) for m in 1:4]
    loop = gmsh.model.occ.addCurveLoop(lines)
    face = gmsh.model.occ.addPlaneSurface([loop])
    d = lp.a * lp.directions[dir]
    ext = gmsh.model.occ.extrude([(2, face)], d.x, d.y, d.z)
    vols = [t for (dim, t) in ext if dim == 3]
    length(vols) == 1 || throw(ProcessingError(
        "Strut prototype extrude for direction $dir produced $(length(vols)) volumes (expected exactly 1)"))
    return Int32(vols[1])
end

"""Build all three strut-direction prototypes once (docs/algorithm.md §3.2)."""
build_prototypes(lp::LatticeParams) = ntuple(dir -> build_prototype(lp, dir), 3)

"""
    instantiate_strut(lp, prototypes, s::StrutRef) -> Int32

Copy the direction-`s.dir` prototype and translate it to strut `s`'s origin
node. O(1): a copy + affine transform, never a re-tessellation of the
profile (docs/algorithm.md §3.2).
"""
function instantiate_strut(lp::LatticeParams, prototypes::NTuple{3,Int32}, s::StrutRef)
    proto = prototypes[s.dir]
    copies = gmsh.model.occ.copy([(3, proto)])
    p = node(lp, s.i, s.j, s.k)
    gmsh.model.occ.translate(copies, p.x, p.y, p.z)
    return Int32(copies[1][2])
end

# --- Boolean operations (docs/algorithm.md §6.3) --------------------------

"""
    fuse_all(tags::Vector{Int32}) -> Vector{Int32}

Multi-operand fuse of every volume in `tags` into as few solids as possible,
via OCCT's n-way GeneralFuse (one call handles the whole set — never
incremental pairwise fusing). Returns the resulting volume tags; there can be
more than one if the inputs are not all mutually connected.
"""
function fuse_all(tags::Vector{Int32})
    isempty(tags) && return Int32[]
    length(tags) == 1 && return Int32[tags[1]]
    objs = [(3, tags[1])]
    tools = [(3, t) for t in tags[2:end]]
    out, _ = gmsh.model.occ.fuse(objs, tools)
    return Int32[t for (dim, t) in out if dim == 3]
end

"""
    common_with(tags::Vector{Int32}, body_tag::Int32; remove_tool::Bool=true) -> Vector{Int32}

Boolean COMMON (intersection) of every volume in `tags` against `body_tag`.

**Operand-disjointness invariant (docs/algorithm.md §6.3):** `tags` must
never contain two volumes whose bounding boxes overlap. OCCT's `intersect`
runs its general boolean algorithm over *all* object operands, together, as
one step — when two objects overlap each other, the result is a *partition*
of both against `body_tag`, not a single union-then-intersect. Measured
directly: intersecting 3 mutually-overlapping struts sharing one lattice
node against a containing box returned **7** fragment solids (four
0.125 mm³ junction wedges plus three 3.16 mm³ pieces) instead of the 1 solid
(9.98 mm³) produced by fusing the struts first, then intersecting. Callers
with potentially-overlapping operands must partition them (e.g. via
`overlap_components`) and either pre-fuse each overlapping group or issue
one COMMON call per group/singleton — never hand the whole group to one
`common_with` call.

`remove_tool=false` keeps `body_tag` alive in the model after the call
(`gmsh.model.occ.intersect`'s `removeTool` option, which otherwise defaults
to consuming it) — required whenever more than one COMMON call reuses the
same imported body in one session, e.g. one call per overlap-component.
"""
function common_with(tags::Vector{Int32}, body_tag::Int32; remove_tool::Bool=true)
    isempty(tags) && return Int32[]
    # gmsh.model.occ.intersect(objectDimTags, toolDimTags, tag=-1, removeObject=true, removeTool=true)
    out, _ = gmsh.model.occ.intersect([(3, t) for t in tags], [(3, body_tag)], -1, true, remove_tool)
    return Int32[t for (dim, t) in out if dim == 3]
end

"""
    cut_with(tags::Vector{Int32}, tool_tags::Vector{Int32}) -> Vector{Int32}

Boolean CUT (subtraction) of every volume in `tool_tags` from every volume
in `tags`. Not used by the generation pipeline itself (which only ever
needs fuse/common) — provided for the golden-sample volume-difference
verification check (tools/verify_geometry.jl, specification.md §6.2).
"""
function cut_with(tags::Vector{Int32}, tool_tags::Vector{Int32})
    (isempty(tags) || isempty(tool_tags)) && return copy(tags)
    out, _ = gmsh.model.occ.cut([(3, t) for t in tags], [(3, t) for t in tool_tags])
    return Int32[t for (dim, t) in out if dim == 3]
end

"""Volume (mm³) of a solid, via OCCT mass properties at unit density."""
volume_of(tag::Integer) = gmsh.model.occ.getMass(3, Int32(tag))

"""Remove entities (and their bounding lower-dim entities) from the model."""
function remove_entities(tags::Vector{Int32})
    isempty(tags) && return nothing
    gmsh.model.occ.remove([(3, t) for t in tags], true)
    return nothing
end

"""Copy and translate an arbitrary volume by a `Vec3` offset (used for the
periodicity shortcut's unit-tile reuse — docs/algorithm.md §6.4)."""
function copy_translate(tag::Integer, offset::Vec3)
    copies = gmsh.model.occ.copy([(3, Int32(tag))])
    gmsh.model.occ.translate(copies, offset.x, offset.y, offset.z)
    return Int32(copies[1][2])
end
