# src/pipeline.jl
#
# Top-level orchestration: import -> classify -> tile -> parallel tile
# stage -> periodicity shortcut -> assembly -> cleanup -> export. Implements
# docs/algorithm.md §4 (pipeline), §6 (tiling/parallelism/fusion), §7
# (auto-tuning/memory model).

using Dates
using Distributed

# --- Temp directory & process priority -------------------------------

"""Create `temp/<yyyymmdd-HHMMSS>/` next to `output_path` (specification.md
§4.4) and return its path."""
function make_run_tempdir(output_path::AbstractString)
    dir = dirname(output_path)
    dir = isempty(dir) ? "." : dir
    ts = Dates.format(now(), "yyyymmdd-HHMMSS")
    tmp = joinpath(dir, "temp", ts)
    mkpath(tmp)
    return tmp
end

"""Request below-normal OS scheduling priority for the current process
(`-bg`/`--background`, docs/algorithm.md §7.3). Best-effort: failures are
swallowed since this is a niceness hint, not a correctness requirement."""
function set_below_normal_priority!()
    @static if Sys.iswindows()
        try
            proc = ccall(:GetCurrentProcess, Ptr{Cvoid}, ())
            ccall((:SetPriorityClass, "kernel32"), Cint, (Ptr{Cvoid}, Culong), proc, 0x00004000)
        catch
        end
    else
        try
            ccall(:nice, Cint, (Cint,), 5)
        catch
        end
    end
    return nothing
end

"""Default per-tile `balanced_fuse!` time budget (docs/algorithm.md §7.1,
§11.2): `max(120, 4 * target_tile_seconds)` using `derive_tile_size`'s own
default `target_tile_seconds=60.0`. `process_tile` makes up to three
separate `balanced_fuse!` calls (interior, boundary, final) and each is
independently budgeted at this value — bounding what was previously an
*unbudgeted* 600s-per-call worst case (three unlogged 10-minute stalls per
oversized tile, invisible in the log since `process_tile` never passed
`rl`) down to something that both fails fast and is now visible via the
tile's returned `diag` (docs/algorithm.md §9)."""
const DEFAULT_TILE_FUSE_BUDGET_SECONDS = max(120.0, 4 * 60.0)

"""Maximum continuous wait, per dispatching task, inside `dispatch_tiles`'
memory-watchdog backpressure loop (docs/algorithm.md §7.2) before dispatch
resumes regardless. `rl.max_rss` is a monotonic high-water mark (see
`update_rss!`/`observe_rss!` in runlog.jl) — it never falls back below the
0.8×budget trip threshold on its own, so without this bound the loop is a
genuine unconditional hang the first time RSS crosses that threshold, not
merely a slow pause."""
const MEMORY_WATCHDOG_MAX_PAUSE_SECONDS = 120.0

tile_brep_path(tmpdir::AbstractString, key::TileKey) =
    joinpath(tmpdir, "tile_$(key.bi)_$(key.bj)_$(key.bk).brep")

"""
    shutdown_workers!(; timeout=2.0) -> Nothing

Tear down the `Distributed` worker pool without ever throwing, so that
cleanup can run unconditionally on every exit path — success, failure, and
Ctrl+C alike — without an error *during cleanup* masking the run's real
outcome.

A plain `rmprocs(workers())` is not safe when a run is being cancelled:
`Distributed._rmprocs` holds the global worker lock while it waits for every
worker to reach a terminated state, and raises
`"rmprocs: pids [...] not terminated after N seconds"` if any does not — the
likely case here, since a worker interrupted mid-tile may be several seconds
deep inside a long-running OCCT boolean that cannot be preempted. So the
clean removal gets a short grace period (`timeout` seconds), and if that
fails the workers are force-removed without waiting (`waitfor=0` hands the
wait off to a background task instead of blocking the exit path), that being
itself wrapped against failure.

`myid()` is filtered out because `workers()` returns `[1]` (the master) when
no workers were ever added, and `rmprocs` refuses to remove the master.
"""
function shutdown_workers!(; timeout::Real=2.0)
    ws = filter(!=(myid()), workers())
    isempty(ws) && return nothing
    try
        rmprocs(ws; waitfor=timeout)
    catch
        # Force-remove whatever is still registered, without waiting. A
        # second Ctrl+C landing in here is swallowed too — at this point the
        # process is on its way out regardless.
        try
            rmprocs(filter(!=(myid()), workers()); waitfor=0)
        catch
        end
    end
    return nothing
end

"""
    determine_worker_count(args::CliArgs) -> Int

Worker count from explicit `--workers`, or derived from `--cores` as
`clamp(cores - 1, 1, 8)` (docs/algorithm.md §7.1) — one core reserved for
the master process and OS/desktop responsiveness. This does not depend on
the calibration probe (only tile size does), so it can be computed
immediately from parsed CLI arguments — notably, before any `Distributed`
worker pool needs to be set up (`@everywhere` can only run at true top
level, so pool setup lives in `main.jl`, not here — see `run_pipeline`).
"""
determine_worker_count(args::CliArgs) = args.cores !== nothing ? clamp(args.cores - 1, 1, 8) : args.workers

"""Part name for the STEP model/product name (specification.md §5):
`<input_stem>+cc<cc>+t<t>` — note this uses `+` separators, distinct from
the `-` convention used for the default *file* name (specification.md §3,
implemented in `resolve_output_paths`)."""
part_name(input::AbstractString, cc::Real, t::Real) =
    "$(splitext(basename(input))[1])+cc$(format_param(cc))+t$(format_param(t))"

# --- Calibration probe & tile sizing (docs/algorithm.md §7.1) -------------

"""
    calibrate(lp) -> (mem_per_strut_bytes, probe_seconds, probe_struts)

Build and fuse a small 4×4×4-cell, all-INTERIOR reference block —
independent of the actual input geometry — to estimate per-strut memory and
fuse-time cost, for auto-tuning tile size. Uses `Sys.maxrss()` (a monotonic
high-water mark, not an instantaneous reading) before/after as a proxy for
this operation's peak memory; run early in a process's lifetime so the
"before" baseline is representative.

Fuses via `balanced_fuse!` — the same code path `process_tile` actually
uses in production — not a single flat `fuse_all` call. Measured directly
(docs/algorithm.md §7.1, §11.2): `fuse_all`'s single-shot cost and
`balanced_fuse!`'s batched cost diverge sharply as strut count grows, so a
probe built on `fuse_all` systematically under-estimates real per-tile fuse
time and produces tile sizes well past the point real tiles blow their fuse
budget.
"""
function calibrate(lp::LatticeParams)
    n = 4
    struts = StrutRef[]
    for i in 0:n-1, j in 0:n-1, k in 0:n-1, d in 1:3
        push!(struts, StrutRef(i, j, k, d))
    end
    rss0 = Sys.maxrss()
    t0 = time()
    with_gmsh() do
        new_model("calibration")
        protos = build_prototypes(lp)
        tags = Int32[instantiate_strut(lp, protos, s) for s in struts]
        remove_entities(collect(protos))   # see process_tile's docstring note on the prototype leak
        balanced_fuse!(tags)
    end
    elapsed = time() - t0
    rss1 = Sys.maxrss()
    nstruts = length(struts)
    mem_per_strut = max(Float64(rss1) - Float64(rss0), 0.0) / nstruts
    return mem_per_strut, elapsed, nstruts
end

"""
    derive_tile_size(mem_per_strut_bytes, probe_seconds, probe_struts, W, ram_gb;
                      target_tile_seconds=60.0) -> (n, n_mem, n_time)

Tile edge (cells per axis), bounded by **both** memory and fuse time
(docs/algorithm.md §7.1, §11.2) — memory alone is not enough: a run that
comfortably fits in RAM can still pick a tile so large its interior/boundary
`balanced_fuse!` calls blow through their own time budget before converging
(measured directly: 192 struts fused in 12.4s, 648 struts in 84.7s — well
past quadratic, so tile size must be actively bounded by time, not just
assumed "more RAM budget = bigger tile is fine").

- `n_mem`: the original memory-only formula — `W * struts_per_tile *
  mem_per_strut ≤ 0.6 * ram_gb`.
- `n_time`: sized so a tile's own fuse is expected to finish within
  `target_tile_seconds`, extrapolating the calibration probe's
  `probe_seconds` for `probe_struts` assuming worst-case quadratic growth
  (`T ∝ S²`) — conservative given the measured ~N^2.5 trend, so this
  intentionally errs toward smaller tiles rather than larger ones.
- The chosen `n = clamp(min(n_mem, n_time), 2, 8)`. The **hard cap of 8**
  matters independently of both formulas: `8³ * 3 = 1536` struts is already
  past the measured superlinear-growth knee, so `n` is never allowed to
  exceed it even if both the memory and time estimates would otherwise
  allow a larger tile.

Returns all three so the caller can log which bound was binding
(docs/algorithm.md §7.1).
"""
function derive_tile_size(mem_per_strut::Real, probe_seconds::Real, probe_struts::Integer,
                           W::Integer, ram_gb::Real; target_tile_seconds::Real=60.0)
    ram_bytes = ram_gb * 1024.0^3
    budget = 0.6 * ram_bytes
    struts_mem = mem_per_strut > 0 ? floor(Int, budget / (W * mem_per_strut)) : typemax(Int) ÷ 4
    struts_mem = max(struts_mem, 3)
    n_mem = clamp(floor(Int, cbrt(struts_mem / 3)), 2, 32)

    struts_time = probe_seconds > 0 ?
        probe_struts * sqrt(target_tile_seconds / probe_seconds) : Float64(typemax(Int) ÷ 4)
    struts_time = max(struts_time, 3.0)
    n_time = clamp(floor(Int, cbrt(struts_time / 3)), 2, 32)

    n = clamp(min(n_mem, n_time), 2, 8)
    return (n=n, n_mem=n_mem, n_time=n_time)
end

"""
    trim_disjoint(tags, body_tag, diag) -> Vector{Int32}

Boolean COMMON of every volume in `tags` against `body_tag`, **never**
handing `common_with` two operands whose AABBs overlap (docs/algorithm.md
§6.3's operand-disjointness invariant — see `common_with`'s docstring for
the measured 7-fragment failure this avoids).

`tags` has typically already been through one `balanced_fuse!` pass (in
`process_tile`), so most of it is expected to already be reduced to mutually
disjoint solids — but `balanced_fuse!`'s no-progress guard, time budget, and
per-group fuse-failure fallback can all legitimately leave some of its
output still mutually overlapping (docs/algorithm.md §6.5). Rather than
assume the input is already clean:

1. Partition `tags` by `overlap_components` (built on the already-computed
   `bounding_boxes`, one sync total).
2. Every singleton component is provably disjoint from every other tag, so
   all singletons are trimmed together in **one** `common_with` call.
3. Every component with more than one member did not converge upstream; a
   diagnostic string is appended to `diag` (surfaced to the run log by the
   caller — docs/algorithm.md §9), and each of its members is trimmed with
   its **own** `common_with` call (one operand — cannot fragment, by
   construction).

`body_tag` is kept alive across every call (`remove_tool=false`) since it is
reused by however many calls step 3 needs, then removed once, explicitly, at
the end.
"""
function trim_disjoint(tags::Vector{Int32}, body_tag::Int32, diag::Vector{String})
    isempty(tags) && return Int32[]
    comps = overlap_components(tags, bounding_boxes(tags))
    out = Int32[]
    singles = Int32[c[1] for c in comps if length(c) == 1]
    isempty(singles) || append!(out, common_with(singles, body_tag; remove_tool=false))
    for c in comps
        length(c) == 1 && continue
        push!(diag, "boundary group of $(length(c)) solids did not converge before COMMON; " *
                     "trimming individually (one COMMON call per solid) to avoid fragmenting them")
        for t in c   # one operand per call -> cannot fragment (docs/algorithm.md §6.3)
            append!(out, common_with(Int32[t], body_tag; remove_tool=false))
        end
    end
    remove_entities(Int32[body_tag])
    return out
end

"""
    drop_tile_local_islands(lp, td, tags, n, origin, threshold) -> (kept, dropped_count)

Worker-side floating-island removal (docs/algorithm.md §6.4a): among a
tile's own final fused solids `tags`, drop any solid that is simultaneously

1. sub-threshold (`volume < threshold`);
2. a singleton in `tags`' own AABB-overlap graph — cannot merge with
   anything else *already in this tile*; and
3. nowhere near any of the tile's `interface_nodes` (padded by the strut
   circumradius `lp.r`) — cannot possibly merge with anything from a
   *neighbouring* tile either, since a solid can only ever share a node
   junction with another tile's geometry at one of those positions.

All three together prove the solid is a genuine floating body — exactly the
same standard `filter_floating!` applies at final assembly time (docs/
algorithm.md §8), checked here early so it never has to survive a tile
write + reimport + a second overlap resolution just to be discarded later.
A small solid that fails any of the three conditions is always kept — this
is a strictly-safe *subset* of `filter_floating!`'s check, not a separate,
weaker one; ambiguous cases are left for the master's fuse-and-resolve pass.
"""
function drop_tile_local_islands(lp::LatticeParams, td::TileDescriptor, tags::Vector{Int32},
                                  n::Integer, origin::NTuple{3,Int}, threshold::Real)
    isempty(tags) && return (tags, 0)
    boxes = bounding_boxes(tags)
    comps = overlap_components(tags, boxes)
    margin = lp.r + 1e-6
    ifn = interface_nodes(lp, td.key, n, origin)

    kept = Int32[]
    to_drop = Int32[]
    for c in comps
        if length(c) != 1
            append!(kept, c)
            continue
        end
        t = c[1]
        if volume_of(t) >= threshold
            push!(kept, t)
            continue
        end
        b = boxes[t]
        near_interface = any(ifn) do p
            b[1] - margin <= p.x <= b[4] + margin &&
            b[2] - margin <= p.y <= b[5] + margin &&
            b[3] - margin <= p.z <= b[6] + margin
        end
        if near_interface
            push!(kept, t)
        else
            push!(to_drop, t)
        end
    end
    # OCC-level removal (not the fast model-level path): counts here are
    # tiny (a handful of slivers per tile at most) and the tile's
    # write_model call right after this still synchronizes normally.
    isempty(to_drop) || remove_entities(to_drop)
    return (kept, length(to_drop))
end

# --- Per-tile worker function (docs/algorithm.md §6.2-6.3) ---------------

"""
    process_tile(lp, input_brep_path, td, out_path, n_tile, tile_origin;
                 max_seconds=DEFAULT_TILE_FUSE_BUDGET_SECONDS)
        -> NamedTuple or nothing

Self-contained per-tile unit of work, safe to call on any process (own gmsh
session lifecycle): instantiate+fuse the tile's INTERIOR struts (one
multi-operand fuse), instantiate+fuse+trim the tile's BOUNDARY struts (one
COMMON against the input body re-imported from `input_brep_path`), fuse both
parts, and write the result to `out_path` if non-empty. Each of the (up to
three) `balanced_fuse!` calls gets its own `max_seconds` budget, rather than
the previous unbudgeted default (docs/algorithm.md §7.1) — a genuinely stuck
tile now fails fast, visibly, per stage.

The returned `NamedTuple` carries per-stage timings (`t_interior`,
`t_boundary`, `t_final` — seconds spent in each `balanced_fuse!`/trim call)
and `diag::Vector{String}` (every warning `balanced_fuse!`/`trim_disjoint`
would otherwise only have been able to log if running on the master, which
holds the log file's `IOStream` — not something that crosses a `Distributed`
call, docs/algorithm.md §9) so the caller (`dispatch_tiles`) can surface
tile-level detail that was previously invisible in the log for any tile run
on a worker.

Returns `nothing` — never an exception — if the process running it is
interrupted (Ctrl+C). Catching the interrupt *inside* the unit of work means
the gmsh session is torn down through `with_gmsh`'s `finally` before this
returns, and the caller gets a plain sentinel value rather than an exception
it has to unwrap out of a `RemoteException` per in-flight tile.

This covers three cases: the sequential (no-worker) path, where the master
runs tiles itself and genuinely does receive Ctrl+C; an explicit
`Distributed.interrupt()`; and any pool whose workers share the console's
signal group. `addprocs` launches its workers *detached* (a new process
group on Windows, a new session on Unix), so a console Ctrl+C normally
reaches only the master — the master-side cancellation in `dispatch_tiles`
plus `shutdown_workers!` is what ends a distributed run, and this handler is
what keeps a worker that *does* see the signal from dying mid-session.
"""
function process_tile(lp::LatticeParams, input_brep_path::AbstractString, td::TileDescriptor,
                       out_path::AbstractString, n_tile::Integer, tile_origin::NTuple{3,Int};
                       max_seconds::Real=DEFAULT_TILE_FUSE_BUDGET_SECONDS)
    t0 = time()
    n_solids = 0
    dropped_islands = 0
    diag = String[]
    t_interior = 0.0
    t_boundary = 0.0
    t_final = 0.0
    try
        with_gmsh() do
            new_model("tile_$(td.key.bi)_$(td.key.bj)_$(td.key.bk)")
            protos = build_prototypes(lp)

            # Instantiate every strut this tile needs (both classes) *before*
            # touching the prototypes again, then drop the 3 prototype
            # solids from the OCC kernel immediately. instantiate_strut only
            # ever copies from a prototype, never consumes it, so the 3
            # originals would otherwise sit in the model, unused, for the
            # rest of this session — and `write_model` exports the *whole*
            # model, not just the tags a caller tracks. Verified directly:
            # left in place, every tile .brep silently gained 3 extra,
            # unclipped, full-length strut solids at world (0,0,0),
            # regardless of the tile's actual location — a real correctness
            # defect (phantom material outside the input volume), not just
            # wasted memory. See also `calibrate`, which has the same
            # pattern (no output file, so it only skews its own
            # mem-per-strut measurement rather than corrupting geometry).
            interior_tags = Int32[instantiate_strut(lp, protos, s) for s in td.interior]
            boundary_tags = Int32[instantiate_strut(lp, protos, s) for s in td.boundary]
            remove_entities(collect(protos))

            # balanced_fuse! (not a single flat fuse_all call) even within one
            # tile: a tile with hundreds+ of struts hitting OCCT's multi-operand
            # fuse in one shot reproduces the same "grows worse than linearly
            # with operand count" pitfall documented for assembly
            # (docs/algorithm.md §11) — observed directly when testing a
            # deliberately oversized single tile during development.
            ti0 = time()
            interior_fused = balanced_fuse!(interior_tags; max_seconds=max_seconds, diag=diag)
            t_interior = time() - ti0

            boundary_result = Int32[]
            if !isempty(td.boundary)
                tb0 = time()
                boundary_fused = balanced_fuse!(boundary_tags; max_seconds=max_seconds, diag=diag)
                body_vols = import_shapes(input_brep_path)
                body_tag = body_vols[1][2]
                # Never hand common_with a group with overlapping operands —
                # trim_disjoint partitions boundary_fused by AABB-overlap
                # component first (docs/algorithm.md §6.3).
                boundary_result = trim_disjoint(boundary_fused, body_tag, diag)
                t_boundary = time() - tb0
            end

            combined = vcat(interior_fused, boundary_result)
            tf0 = time()
            final = balanced_fuse!(combined; max_seconds=max_seconds, diag=diag)
            t_final = time() - tf0

            # Worker-side floating-island removal (docs/algorithm.md §6.4a,
            # §11.2): a sub-threshold solid that cannot touch *any*
            # neighbouring tile is provably a genuine floating body already
            # at this stage, cheaper to drop here than to carry it all the
            # way to the master's export-stage `filter_floating!` pass.
            kept, dropped_islands = drop_tile_local_islands(lp, td, final, n_tile, tile_origin, lp.t^3)

            n_solids = length(kept)
            n_solids > 0 && write_model(String(out_path))
        end
    catch e
        # Ctrl+C on this process: unwind quietly (the gmsh session is already
        # finalized by with_gmsh's `finally`) and report it as a sentinel, not
        # as an exception the master has to unwrap out of a RemoteException.
        is_interrupt(e) && return nothing
        rethrow()
    end
    elapsed = time() - t0
    return (key=td.key, interior_count=length(td.interior), boundary_count=length(td.boundary),
            solids=n_solids, elapsed=elapsed, rss=Int(Sys.maxrss()), wrote=n_solids > 0,
            t_interior=t_interior, t_boundary=t_boundary, t_final=t_final,
            dropped_islands=dropped_islands, diag=diag)
end

# --- Dispatch pool with memory backpressure (docs/algorithm.md §7.2) -----

"""
    dispatch_jobs(n, worker_call, sequential_call; rl, cancel_msg,
                  ram_budget_bytes=nothing, rss_of=res->res.rss,
                  on_result=(idx,res,wid)->nothing) -> Vector{Any}

Generic work-stealing dispatch of jobs `1:n` across `Distributed.workers()`,
or sequentially on the master when none are registered — the pool
management, memory-watchdog backpressure, and Ctrl+C handling extracted out
of what was originally `dispatch_tiles`'s only caller, so the distributed
assembly merge rounds (docs/algorithm.md §6.5) can share the *exact* same,
already-worked-out cancellation semantics instead of a hand-duplicated copy.

- `worker_call(wid, idx)` / `sequential_call(idx)` must each run job `idx`
  (remotely / locally) and return its result, or `nothing` if the job itself
  was interrupted (the sentinel convention `process_tile` established,
  docs/algorithm.md §6.2 — every dispatched unit-of-work function follows
  it).
- `rss_of(result)` extracts the peak RSS the job reported, fed into the
  shared memory-watchdog backpressure loop and the run's overall
  `max_rss` (docs/algorithm.md §7.2); defaults to `.rss`, which every job
  result in this codebase carries (`process_tile`, `merge_group`).
- `on_result(idx, result, wid)` runs once per successfully completed job
  (`wid` is `:sequential` on the no-worker path, otherwise the worker id) —
  the caller's place to log/record anything job-specific, in whatever order
  jobs actually finish (not necessarily index order).

**Cancellation (Ctrl+C):** an interrupt can arrive here in three different
shapes — raised in the master's own task while it waits in `@sync`, raised
inside one of the dispatching `@async` tasks (surfacing as a
`TaskFailedException`/`CompositeException`), or seen first by a worker
mid-call (`worker_call` returns `nothing`, or the call fails with a
`RemoteException` wrapping the interrupt). All three set a shared
cancellation flag that stops the remaining tasks from pulling further work,
and dispatch then ends with a single `CancelledError` rather than a burst of
unrelated-looking worker errors.
"""
function dispatch_jobs(n::Integer, worker_call::Function, sequential_call::Function;
                        rl::RunLog, cancel_msg::AbstractString,
                        ram_budget_bytes::Union{Nothing,Real}=nothing,
                        rss_of::Function=res -> res.rss,
                        on_result::Function=(idx, res, wid) -> nothing)
    results = Vector{Any}(undef, n)
    ws = workers()

    if isempty(ws) || ws == [1]
        # Sequential: no Distributed workers were added (single-worker run).
        for idx in 1:n
            res = try
                sequential_call(idx)
            catch e
                is_interrupt(e) || rethrow()
                nothing
            end
            res === nothing && throw(CancelledError(cancel_msg))
            results[idx] = res
            observe_rss!(rl, rss_of(res))
            on_result(idx, res, :sequential)
        end
        return results
    end

    next_idx = Ref(1)
    idx_lock = ReentrantLock()
    rss_lock = ReentrantLock()
    cancelled = Ref(false)

    try
        @sync for wid in ws
            @async begin
                while !cancelled[]
                    idx = lock(idx_lock) do
                        if next_idx[] > n
                            nothing
                        else
                            i = next_idx[]
                            next_idx[] += 1
                            i
                        end
                    end
                    idx === nothing && break

                    if ram_budget_bytes !== nothing
                        wait_start = time()
                        while !cancelled[] && lock(() -> rl.max_rss, rss_lock) > 0.8 * ram_budget_bytes
                            # rl.max_rss (docs/algorithm.md §7.2) is a
                            # monotonic high-water mark — it can only grow,
                            # never fall back below the threshold on its
                            # own — so an unbounded wait here is a genuine
                            # hang once RSS crosses 0.8x budget even once,
                            # not just a slow pause. Bound it the same way
                            # balanced_fuse!'s own time budget is a
                            # last-resort circuit breaker, not a routine
                            # speed control (docs/algorithm.md §6.5).
                            if time() - wait_start > MEMORY_WATCHDOG_MAX_PAUSE_SECONDS
                                log_line(rl, "Memory watchdog: pause budget " *
                                              "($(MEMORY_WATCHDOG_MAX_PAUSE_SECONDS)s) exhausted " *
                                              "(max_rss=$(format_bytes(rl.max_rss)) still over 0.8x budget); " *
                                              "continuing dispatch anyway."; console=true)
                                break
                            end
                            log_line(rl, "Memory watchdog: pausing new job dispatch (max_rss=$(format_bytes(rl.max_rss)))")
                            sleep(1.0)
                        end
                    end
                    cancelled[] && break

                    res = try
                        worker_call(wid, idx)
                    catch e
                        # An interrupt seen by the worker (or by this task
                        # while waiting on it) ends the whole dispatch. Once
                        # cancellation is already under way, a worker
                        # vanishing mid-call is expected too — the pool is
                        # being torn down. Anything else is a genuine job
                        # failure and must propagate.
                        (is_interrupt(e) || cancelled[]) || rethrow()
                        nothing
                    end
                    # `nothing` = interrupted, either remotely (the job
                    # function's sentinel) or locally (the catch above).
                    if res === nothing
                        cancelled[] = true
                        break
                    end

                    results[idx] = res
                    lock(rss_lock) do
                        observe_rss!(rl, rss_of(res))
                    end
                    on_result(idx, res, wid)
                end
            end
        end
    catch e
        # Ctrl+C delivered to the master's own task while it waited in @sync,
        # or re-raised out of an @async task. Flag it so any task still
        # running stops pulling new work, then report it as a cancellation.
        is_interrupt(e) || rethrow()
        cancelled[] = true
    end

    cancelled[] && throw(CancelledError(cancel_msg))
    return results
end

"""
    dispatch_tiles(lp, input_brep, tiles, tmpdir, rl, ram_budget_bytes, n_tile, tile_origin) -> Vector

Run `process_tile` over every `TileDescriptor` in `tiles` via `dispatch_jobs`
(docs/algorithm.md §6.2), logging each tile's collected `diag` warnings and a
per-tile summary line (including the per-stage timings and dropped-island
count `process_tile` now reports, docs/algorithm.md §7.1, §6.4a) as it
completes.
"""
function dispatch_tiles(lp::LatticeParams, input_brep::AbstractString, tiles::Vector{TileDescriptor},
                         tmpdir::AbstractString, rl::RunLog, ram_budget_bytes::Union{Nothing,Real},
                         n_tile::Integer, tile_origin::NTuple{3,Int})
    cancel_msg = "Run cancelled by user (Ctrl+C) during the tile stage; " *
                 "partial results are left in $tmpdir for analysis."

    sequential_call(idx) = process_tile(lp, input_brep, tiles[idx],
                                         tile_brep_path(tmpdir, tiles[idx].key), n_tile, tile_origin)
    worker_call(wid, idx) = remotecall_fetch(process_tile, wid, lp, input_brep, tiles[idx],
                                              tile_brep_path(tmpdir, tiles[idx].key), n_tile, tile_origin)

    on_result = (idx, stat, wid) -> begin
        # Every diagnostic string a tile's balanced_fuse!/trim_disjoint
        # calls collected (docs/algorithm.md §9) — previously invisible
        # whenever the tile ran on a worker, since process_tile had no
        # RunLog to write to directly.
        for msg in stat.diag
            log_line(rl, "tile $(stat.key): $msg")
        end
        where = wid === :sequential ? "(sequential)" : "on worker $wid"
        log_line(rl, "tile $(stat.key) done $where: interior=$(stat.interior_count) " *
                      "boundary=$(stat.boundary_count) solids=$(stat.solids) time=$(format_duration(stat.elapsed)) " *
                      "(t_interior=$(format_duration(stat.t_interior)) t_boundary=$(format_duration(stat.t_boundary)) " *
                      "t_final=$(format_duration(stat.t_final)) dropped_islands=$(stat.dropped_islands))")
    end

    return dispatch_jobs(length(tiles), worker_call, sequential_call;
                          rl=rl, cancel_msg=cancel_msg, ram_budget_bytes=ram_budget_bytes,
                          on_result=on_result)
end

# --- Periodicity shortcut + hierarchical assembly (docs/algorithm.md §6.4-6.5) --

"""Do any two of the given volumes' AABBs overlap (with a tiny numerical
epsilon for floating-point safety)? A cheap pre-filter — pure box-coordinate
comparison, no OCCT boolean math — used to skip fuse attempts on groups that
can be *proven* to share no geometry at all: struts sharing a lattice node
have bit-identical endpoint coordinates (`node()` is deterministic), so
genuinely touching solids' AABBs always overlap; there is no "almost
touching" case this could wrongly skip.

`boxes` must already contain every tag in `tags` — see `bounding_boxes`.
**Callers must not call `bounding_box`/`bounding_boxes` once per group**:
`gmsh.model.occ.synchronize()` (inside every `bounding_box` call) is
O(whole model), not O(this query) — measured at 202× the cost of one shared
sync followed by raw box lookups (docs/algorithm.md §6.5). `boxes` should be
computed exactly once per fuse *round* (`balanced_fuse!` below), not once
per group within a round."""
function _group_might_touch(tags::Vector{Int32}, boxes::Dict{Int32,NTuple{6,Float64}})
    n = length(tags)
    n <= 1 && return false
    for i in 1:n, j in i+1:n
        aabbs_overlap(boxes[tags[i]], boxes[tags[j]]) && return true
    end
    return false
end

"""
    balanced_fuse!(tags::Vector{Int32}; rl::Union{RunLog,Nothing}=nothing) -> Vector{Int32}

Fuse every volume in `tags` (already present in the active gmsh session)
using a balanced binary tree of pairwise-group fuses rather than one flat
`length(tags)`-operand fuse, bounding peak operand complexity per call
(docs/algorithm.md §6.5). Operates in batches of up to 8 to keep each
individual `fuse_all` call's operand count modest.

Before calling OCCT's boolean fuse on a group, `_group_might_touch` cheaply
rules out groups with no AABB overlap at all — nothing to merge, so the
group passes straight through unfused without ever invoking the (much more
expensive) boolean kernel.

If OCCT's boolean fuse throws on a particular group (a known robustness
failure mode — e.g. "Courbes non jointives" — that can occur on certain
degenerate/sliver inputs, most often fragments produced by trimming
boundary struts against a curved surface), that group is passed through
**un-fused** rather than aborting the whole run: specification.md §1
explicitly allows the output to contain multiple disconnected bodies. This
is unconditionally safe *only* for a group `_group_might_touch` already
ruled out (provably zero AABB overlap — nothing there could ever have been
self-intersecting). For a group that *does* overlap but whose fuse fails or
times out, the un-fused fallback can leave genuinely overlapping solids in
the output — struts sharing a lattice node are built to overlap near the
shared vertex until fused (docs/algorithm.md §3.2), so this is a real
self-intersection risk, not just "more bodies than ideal". This was caught
during development via `tools/verify_geometry.jl`'s self-intersection check:
an over-tight `max_seconds` cut assembly off before finishing genuinely
necessary fuses and produced ~65K self-intersecting triangle pairs. `catch`
here is a fallback for OCCT's own robustness limits — not a substitute for
giving the fuse process enough time to actually finish; see `max_seconds`.

**Termination guard:** if an entire round makes *no* progress at all (every
group either failed or fused to the same count it started with — e.g. a
persistently un-fusable sliver that keeps landing in a fresh group of
otherwise-unrelated solids each round), further rounds would regroup the
exact same unreduced elements into the exact same batches and repeat the
exact same failures forever. This was caught during development: a run with
several persistently un-fusable fragments hung indefinitely. The fix is to
detect a no-progress round and stop immediately, returning the current set
as separate bodies, rather than looping without bound.

**Wall-clock budget (`max_seconds`):** the no-progress guard bounds the
*number* of rounds, but not the total *time* — a round that makes some
progress (e.g. shrinking by one solid) still counts as progress and would be
retried, and a single slow, ultimately-failing OCCT attempt can itself take
several real seconds. So many "slow but nonzero progress" rounds can still
add up to a long run in aggregate. `max_seconds` (checked before each round
and, within a round, before each group) is a **last-resort circuit
breaker against a runaway assembly**, not a routine speed control — set
generously (default 10 minutes), because cutting fusing off early is not a
free "just produces more bodies" fallback here: it can leave genuinely
overlapping strut-junction solids un-fused, i.e. real self-intersections in
the output (see the note above). Prefer letting this run to completion;
only lower it if you have independently confirmed the remaining work truly
cannot make progress (in which case the no-progress guard would have
stopped it anyway).
"""
function balanced_fuse!(tags::Vector{Int32}; rl::Union{RunLog,Nothing}=nothing, max_seconds::Real=600.0,
                         diag::Union{Nothing,Vector{String}}=nothing)
    isempty(tags) && return Int32[]
    # Every warning goes to both sinks: `rl` when this call is happening on
    # the master (which owns the log file directly), and `diag` when it's
    # happening on a worker (an `IOStream` doesn't cross a `Distributed`
    # call — diag is a plain Vector{String} that does, surfaced by the
    # caller, e.g. process_tile/dispatch_tiles, docs/algorithm.md §9).
    _warn(msg) = begin
        rl !== nothing && log_line(rl, msg)
        diag !== nothing && push!(diag, msg)
    end
    current = tags
    batch = 8
    t_start = time()
    _budget_msg(n) = "Warning: balanced_fuse! exceeded its $(max_seconds)s time budget with $n " *
                      "solid(s) remaining; stopping further fuse attempts and keeping them as separate bodies."
    while length(current) > batch
        if time() - t_start > max_seconds
            _warn(_budget_msg(length(current)))
            return current
        end
        # One shared sync + box lookup per *round*, not per group (§6.5 —
        # see _group_might_touch's docstring for the measured 202× cost of
        # getting this wrong).
        boxes = bounding_boxes(current)
        next = Int32[]
        progressed = false
        for i in 1:batch:length(current)
            if time() - t_start > max_seconds
                _warn(_budget_msg(length(current) - i + 1))
                append!(next, current[i:end])
                return next
            end
            group = current[i:min(i + batch - 1, length(current))]
            if !_group_might_touch(group, boxes)
                # Provably nothing to merge (no AABB overlap at all) — skip
                # the expensive OCCT boolean call entirely rather than
                # letting it slowly discover the same thing (docs/algorithm.md §11).
                append!(next, group)
                continue
            end
            try
                fused = fuse_all(group)
                length(fused) < length(group) && (progressed = true)
                append!(next, fused)
            catch e
                # A Ctrl+C is not an OCCT robustness failure: swallowing it
                # here would silently downgrade a cancelled run into a
                # "kept as separate bodies" result, i.e. exactly the
                # un-fused, possibly self-intersecting output this fallback
                # is documented to be a last resort for.
                is_interrupt(e) && rethrow()
                _warn("Warning: fuse of a $(length(group))-solid group failed " *
                      "($(sprint(showerror, e))); keeping them as separate bodies.")
                append!(next, group)
            end
        end
        if !progressed
            _warn("Warning: a fuse round made no progress on $(length(current)) " *
                  "solid(s); stopping further fuse attempts and keeping them as separate bodies.")
            return next
        end
        current = next
    end
    if time() - t_start > max_seconds || !_group_might_touch(current, bounding_boxes(current))
        return current
    end
    try
        fused = fuse_all(current)
        return isempty(fused) ? current : fused
    catch e
        is_interrupt(e) && rethrow()   # cancellation, not a kernel failure — see above
        _warn("Warning: final fuse of a $(length(current))-solid group failed " *
              "($(sprint(showerror, e))); keeping them as separate bodies.")
        return current
    end
end

"""
    merge_group(specs, out_path, max_seconds) -> NamedTuple or nothing

Self-contained per-merge-round unit of work (docs/algorithm.md §6.5),
structurally mirroring `process_tile`: import every `.brep` path in `specs`
(a `Vector{Tuple{String,Vec3}}` of `(path, offset)` pairs), translate each
imported volume by its paired offset (`Vec3(0,0,0)` skips the
copy_translate — a "real" tile or a previous round's own merge result is
already at its correct world position; only a periodicity-shortcut
reference-tile duplicate carries a non-zero offset), `balanced_fuse!` the
combined set, and write `out_path` if the result is non-empty. Safe to call
on any process (own gmsh session lifecycle) — this is what lets assembly be
distributed across `Distributed` workers via `dispatch_jobs`, the same way
`process_tile` distributes the tile stage.

Returns `nothing` — never an exception — on Ctrl+C (docs/algorithm.md §6.2's
sentinel convention, so the caller doesn't have to unwrap a
`RemoteException` per in-flight merge).
"""
function merge_group(specs::Vector{Tuple{String,Vec3}}, out_path::AbstractString, max_seconds::Real)
    t0 = time()
    n_solids = 0
    diag = String[]
    try
        with_gmsh() do
            new_model("merge")
            all_tags = Int32[]
            for (path, offset) in specs
                vols = import_shapes(path)
                tags = Int32[t for (_, t) in vols]
                if offset.x == 0.0 && offset.y == 0.0 && offset.z == 0.0
                    append!(all_tags, tags)
                else
                    for t in tags
                        push!(all_tags, copy_translate(t, offset))
                    end
                    remove_entities(tags)   # drop the untranslated originals
                end
            end
            fused = balanced_fuse!(all_tags; max_seconds=max_seconds, diag=diag)
            n_solids = length(fused)
            n_solids > 0 && write_model(String(out_path))
        end
    catch e
        is_interrupt(e) && return nothing
        rethrow()
    end
    elapsed = time() - t0
    return (n_solids=n_solids, elapsed=elapsed, rss=Int(Sys.maxrss()), diag=diag)
end

"""
    do_merge_round(specs, round_idx, tmpdir, rl, max_seconds) -> Dict{TileKey,Tuple{String,Vec3}}

One round of the distributed hierarchical assembly (docs/algorithm.md
§6.5): group `specs` (a `TileKey => (path, offset)` map — the previous
round's outputs, or round 0's per-tile inputs) by 2×2×2 super-blocks
(`fld(bi,2), fld(bj,2), fld(bk,2)`), and for every super-block with more
than one member, dispatch one `merge_group` call (via `dispatch_jobs`,
spreading rounds across `Distributed` workers exactly like the tile stage)
that fuses them into a single `merge_r<round_idx>_<bi>_<bj>_<bk>.brep`. A
super-block with exactly one member has nothing to merge — its
`(path, offset)` is carried forward unchanged rather than spending a
worker call and a file write on a no-op.

Because every `TileKey` produced by `partition_tiles`/`tile_key` is bucketed
relative to the candidate range's own minimum corner (`tile_key`'s
docstring), every `bi/bj/bk` here is `>= 0`, so `fld(bi,2)` behaves as plain
integer division with no risk of the origin-centered spurious-split bug
`tile_key` itself works around.

Returns the next round's `TileKey => (path, offset)` map, keyed by
super-block coordinate. A super-block whose `merge_group` call produced zero
solids (all-empty input, not expected in practice but guarded regardless) is
simply dropped.
"""
function do_merge_round(specs::Dict{TileKey,Tuple{String,Vec3}}, round_idx::Integer,
                         tmpdir::AbstractString, rl::RunLog, max_seconds::Real)
    groups = Dict{TileKey,Vector{Tuple{TileKey,String,Vec3}}}()
    for (key, (path, offset)) in specs
        gkey = TileKey(fld(key.bi, 2), fld(key.bj, 2), fld(key.bk, 2))
        push!(get!(() -> Tuple{TileKey,String,Vec3}[], groups, gkey), (key, path, offset))
    end

    result = Dict{TileKey,Tuple{String,Vec3}}()
    to_merge = TileKey[]
    for gkey in sort(collect(keys(groups)); by=k -> (k.bi, k.bj, k.bk))
        members = groups[gkey]
        if length(members) == 1
            _, path, offset = members[1]
            result[gkey] = (path, offset)
        else
            push!(to_merge, gkey)
        end
    end
    isempty(to_merge) && return result

    _member_specs(gkey) = [(path, offset) for (_, path, offset) in
                            sort(groups[gkey]; by=m -> (m[1].bi, m[1].bj, m[1].bk))]
    out_paths = Dict{TileKey,String}(gkey => joinpath(tmpdir, "merge_r$(round_idx)_$(gkey.bi)_$(gkey.bj)_$(gkey.bk).brep")
                                      for gkey in to_merge)

    cancel_msg = "Run cancelled by user (Ctrl+C) during assembly round $round_idx; " *
                 "partial results are left in $tmpdir for analysis."
    sequential_call(i) = merge_group(_member_specs(to_merge[i]), out_paths[to_merge[i]], max_seconds)
    worker_call(wid, i) = remotecall_fetch(merge_group, wid, _member_specs(to_merge[i]),
                                            out_paths[to_merge[i]], max_seconds)
    on_result = (i, res, wid) -> begin
        gkey = to_merge[i]
        for msg in res.diag
            log_line(rl, "assembly round $round_idx group $gkey: $msg")
        end
        where = wid === :sequential ? "(sequential)" : "on worker $wid"
        log_line(rl, "assembly round $round_idx group $gkey done $where: " *
                      "solids=$(res.n_solids) time=$(format_duration(res.elapsed))")
    end

    merge_results = dispatch_jobs(length(to_merge), worker_call, sequential_call;
                                   rl=rl, cancel_msg=cancel_msg, on_result=on_result)
    for (i, gkey) in enumerate(to_merge)
        res = merge_results[i]
        res.n_solids > 0 && (result[gkey] = (out_paths[gkey], Vec3(0.0, 0.0, 0.0)))
    end
    return result
end

"""
    assemble!(lp, tile_results, ref_key, full_interior_keys, n, tmpdir, rl; max_seconds=600.0) -> Vector{Int32}

Distributed hierarchical assembly (docs/algorithm.md §6.5): builds the
round-0 `TileKey => (path, offset)` map (every tile that wrote a `.brep` at
`offset=Vec3(0,0,0)`, plus one entry per non-reference full-interior tile
pointing at the *reference* tile's `.brep` with the periodicity-shortcut
translation as its offset — docs/algorithm.md §6.4), then repeatedly calls
`do_merge_round` — each round dispatched across `Distributed` workers via
`dispatch_jobs`, exactly like the tile stage — until at most 8 files remain
or a round fails to reduce the file count (the same no-progress guard
`balanced_fuse!` itself uses; whatever remains is simply handed to the final
step below, which has its own complete safety net regardless of how many
files that is).

This keeps peak *master* memory bounded to the last round's survivors
(typically <= 8 already-substantially-fused files) instead of every tile
`.brep` at once — the original single-session design imported and held all
of them simultaneously.

**Requires an already-open gmsh session on the master**, with the model the
caller wants the result in already created via `new_model` — this function
does not open, close, or rename the session itself (previously it did, and
wrote its own result to a `tmpdir/assembled.brep` staging file that
`run_pipeline`'s export stage immediately re-imported into a *second* fresh
session, purely to get the right STEP part name; merging the two sessions
removes that full OCCT serialize+reparse of the entire finished lattice, at
the cost of the caller now being responsible for verifying the model
contains exactly what it expects — see `model_solids` — before writing it
out). The final step imports every surviving `(path, offset)` entry,
translates the non-zero-offset ones, and `balanced_fuse!`s everything
together into the caller's ambient model. Returns the resulting top-level
solid tags, live in that session. The volume-threshold cleanup gate is
**not** applied here — it runs exactly once, later, in the export stage's
`filter_floating!` (docs/algorithm.md §8), so there is a single place in the
pipeline that decides what gets discarded.
"""
function assemble!(lp::LatticeParams, tile_results::Dict{TileKey,Any}, ref_key::Union{TileKey,Nothing},
                    full_interior_keys::Vector{TileKey}, n::Integer, tmpdir::AbstractString, rl::RunLog;
                    max_seconds::Real=600.0)
    specs = Dict{TileKey,Tuple{String,Vec3}}()
    for (key, stat) in tile_results
        stat.wrote && (specs[key] = (tile_brep_path(tmpdir, key), Vec3(0.0, 0.0, 0.0)))
    end
    if ref_key !== nothing
        ref_path = tile_brep_path(tmpdir, ref_key)
        for key in full_interior_keys
            key == ref_key && continue
            specs[key] = (ref_path, tile_translation(lp, ref_key, key, n))
        end
    end

    log_line(rl, "Assembly: $(length(specs)) tile-level input(s) before distributed merge " *
                  "($(length(full_interior_keys)) full-interior tile(s) via periodicity shortcut)"; console=true)

    round_idx = 1
    while length(specs) > 8
        before = length(specs)
        next_specs = do_merge_round(specs, round_idx, tmpdir, rl, max_seconds)
        after = length(next_specs)
        log_line(rl, "Assembly round $round_idx: $before -> $after file(s)"; console=true)
        specs = next_specs
        if after >= before
            log_line(rl, "Assembly: distributed merge rounds plateaued at $after file(s) after " *
                          "$round_idx round(s) (no further 2x2x2 super-block grouping reduces the " *
                          "file count); finishing with a direct fuse on the master."; console=true)
            break
        end
        round_idx += 1
    end

    all_tags = Int32[]
    for (key, (path, offset)) in sort(collect(specs); by=kv -> (kv[1].bi, kv[1].bj, kv[1].bk))
        vols = import_shapes(path)
        tags = Int32[t for (_, t) in vols]
        if offset.x == 0.0 && offset.y == 0.0 && offset.z == 0.0
            append!(all_tags, tags)
        else
            for t in tags
                push!(all_tags, copy_translate(t, offset))
            end
            remove_entities(tags)
        end
    end
    return balanced_fuse!(all_tags; rl=rl, max_seconds=max_seconds)
end

"""
    assert_no_stray_solids(final_tags::Vector{Int32})

Defensive check, meant to be called immediately after `assemble!` returns
and before anything else touches the session: the model's dim-3 entity list
(`model_solids()`) must be *exactly* `final_tags`, no more and no fewer.
Throws `ProcessingError` naming the offending tags otherwise, refusing to
export ambiguous geometry.

Now that assembly and export share one gmsh session (docs/algorithm.md §6.5)
instead of assembly writing `assembled.brep` for export to re-import into a
fresh model, this replaces the implicit guarantee a virgin re-import used to
give for free — that the model contains exactly the returned tags and
nothing else. It exists to catch the same class of bug as the historical
`build_prototypes` leak that once silently exported un-consumed intermediate
solids into every tile `.brep`.

**Call only when safe to synchronize** — i.e. before any model-level-only
removal (`remove_model_entities`, `filter_floating!`) has happened yet in
this session; `model_solids`'s own `sync_model()` call would otherwise
resurrect exactly what such a removal just deleted (see `remove_model_entities`'s
docstring)."""
function assert_no_stray_solids(final_tags::Vector{Int32})
    stray = setdiff(model_solids(), final_tags)
    isempty(stray) || throw(ProcessingError(
        "Assembly left $(length(stray)) unexpected solid(s) in the model that are not " *
        "part of the fused lattice result (tags $stray); refusing to export ambiguous " *
        "geometry. This indicates a leaked intermediate solid upstream of assembly."))
    return nothing
end

# --- Floating-body cleanup (docs/algorithm.md §8, §11.2) -----------------

"""Classify one already-resolved solid `t` as `kept` (>= `threshold`) or
`floating` (< `threshold`, and its volume recorded), pushing into the
caller's accumulator vectors."""
function _classify_solid!(t::Int32, threshold::Real, kept::Vector{Int32},
                           floating::Vector{Int32}, removed_vols::Vector{Float64})
    v = volume_of(t)
    if v < threshold
        push!(floating, t)
        push!(removed_vols, v)
    else
        push!(kept, t)
    end
    return nothing
end

"""
    filter_floating!(tags, threshold; rl=nothing) -> (kept, n_removed, removed_vols)

Decide which of the model's top-level solids `tags` are safe to delete as
sub-threshold "floating" bodies (specification.md §5), and remove exactly
those from the model. This replaces the earlier unconditional "volume <
threshold -> delete" rule, which was deleting connected junction material
whenever upstream fusing hadn't fully converged (docs/algorithm.md §6.3,
§11.2) — the root cause of the original `test-cylinder-cc5t1` run's
multi-hour sub-threshold-removal tail.

Algorithm, mirroring the operand-disjointness discipline used everywhere
else booleans touch potentially-overlapping solids (docs/algorithm.md §6.3):

1. Partition `tags` by `overlap_components` (one `bounding_boxes` sync).
2. A singleton component is provably disjoint from everything else — safe
   to classify directly by volume.
3. A component with more than one member is *ambiguous*: AABB overlap does
   not prove genuine geometric contact. Resolve it by attempting
   `fuse_all` on the whole component in one call: OCCT's n-ary fuse fully
   unions every genuinely-touching subset within one call and leaves
   genuinely disjoint solids with their original tags/shapes untouched
   (verified directly — two truly disjoint boxes come back from `fuse_all`
   with identical tags; two boxes sharing an exact common face merge into
   one with the exact sum of their volumes). A *successful* fuse call is
   therefore always trustworthy either way: every resulting tag can be
   classified directly by volume, with no further ambiguity, because (a)
   OCCT already resolved all connectivity *within* the group, and (b) the
   group was already proven disjoint from every *other* component by the
   AABB partition in step 1.
4. If the fuse attempt throws (a known OCCT robustness limit — see
   `balanced_fuse!`'s docstring), the component's connectivity cannot be
   resolved. Any of its members that are sub-threshold are **unresolved**:
   this function refuses to guess, and raises `ProcessingError` (exit 4,
   docs/algorithm.md §9) rather than either silently deleting possibly-
   connected material or silently keeping it and calling the run a success.
   Priority #1 (precision) over completing the run (docs/algorithm.md's
   "Key Considerations").
5. Every remaining sub-threshold, fully-resolved solid is now a genuine
   singleton — safe to remove. Removed via `remove_model_entities` (the
   fast model-level path, docs/algorithm.md §8); the caller must follow up
   with `write_model(path; sync=false)`, never a syncing write.

`fuse_fn` defaults to `fuse_all` and exists so tests can inject a
deliberately-failing stub to exercise the hard-fail path (step 4) without
needing to manufacture a genuine OCCT boolean-robustness failure.

Returns `(kept, n_removed, removed_vols)`: `kept` is every surviving tag
(sub-threshold survivors of an ambiguous-but-still->=threshold component
included — a large solid can legitimately share an overlap component with a
tiny one, e.g. a main body plus one sliver that fuse legitimately left
separate), `n_removed` the count of deleted floating bodies, and
`removed_vols` their individual volumes (for the caller's aggregate log
line — never one log line per solid, docs/algorithm.md §9).
"""
function filter_floating!(tags::Vector{Int32}, threshold::Real; rl::Union{RunLog,Nothing}=nothing,
                           fuse_fn::Function=fuse_all)
    isempty(tags) && return (Int32[], 0, Float64[])

    boxes = bounding_boxes(tags)
    comps = overlap_components(tags, boxes)

    kept = Int32[]
    floating = Int32[]
    removed_vols = Float64[]
    unresolved_vols = Float64[]

    for c in comps
        if length(c) == 1
            _classify_solid!(c[1], threshold, kept, floating, removed_vols)
            continue
        end
        try
            fused = fuse_fn(c)
            for t in fused
                _classify_solid!(t, threshold, kept, floating, removed_vols)
            end
        catch e
            is_interrupt(e) && rethrow()
            rl !== nothing && log_line(rl, "Warning: fuse attempt on a $(length(c))-solid overlap-ambiguous " *
                                            "cleanup group failed ($(sprint(showerror, e))); its sub-threshold " *
                                            "member(s), if any, cannot be proven floating.")
            for t in c
                v = volume_of(t)
                v < threshold ? push!(unresolved_vols, v) : push!(kept, t)
            end
        end
    end

    if !isempty(unresolved_vols)
        sample = first(sort(unresolved_vols), min(5, length(unresolved_vols)))
        throw(ProcessingError(
            "$(length(unresolved_vols)) sub-threshold solid(s) remain connected to other geometry " *
            "and could not be resolved by fusing (sample volumes $(round.(sample, digits=4)) mm^3, " *
            "threshold=$(round(threshold, digits=4)) mm^3). Refusing to delete possibly-connected " *
            "material (docs/algorithm.md §8, §11.2) — this indicates an unfused boundary/junction " *
            "group upstream of assembly; inspect the kept temp directory for per-tile/assembly " *
            "diagnostics."))
    end

    isempty(floating) || remove_model_entities(floating)
    return (kept, length(floating), removed_vols)
end

# --- Top-level orchestration (docs/algorithm.md §4) -----------------------

"""
    run_pipeline(args::CliArgs, rl::RunLog) -> Dict{Symbol,Any}

Top-level orchestration implementing the full pipeline (docs/algorithm.md
§4): import -> tessellate -> classify -> tile -> parallel tile stage ->
periodicity shortcut + assembly -> cleanup -> export -> verify -> temp
cleanup. Returns the stats dict consumed by `log_summary`. On success, the
temp directory is removed; on any exception, it is left in place for
analysis (specification.md §4.4) and the exception propagates to the
caller (`main.jl`), which maps it to the appropriate exit code.

Any `Distributed` worker pool must already be set up by the caller before
this is called (`main.jl`, at true top level — `@everywhere` cannot run
inside a function body, so pool setup/teardown cannot live here). This
function only *consumes* whatever `Distributed.workers()` currently
reports; with none added, `dispatch_tiles` runs sequentially in-process.
"""
function run_pipeline(args::CliArgs, rl::RunLog)
    args.background && set_below_normal_priority!()

    lp = LatticeParams(args.cc, args.t)
    tmpdir = make_run_tempdir(args.output)
    log_line(rl, "Temp directory: $tmpdir"; console=true)
    input_brep = joinpath(tmpdir, "input.brep")

    stats = Dict{Symbol,Any}()

    lo, hi = @timed_stage rl "import" begin
        with_gmsh() do
            new_model("input")
            import_shapes(args.input)
            write_model(input_brep)
            bounding_box()
        end
    end
    log_line(rl, "Input bounding box: lo=($(lo.x), $(lo.y), $(lo.z)) hi=($(hi.x), $(hi.y), $(hi.z))")

    mesh, sh = @timed_stage rl "tessellate" begin
        with_gmsh() do
            new_model("tess")
            import_shapes(args.input)
            m = tessellate_surface(lp)
            (m, build_spatial_hash(m))
        end
    end
    dtol = mesh_chordal_target(lp)
    log_line(rl, "Surface mesh: $(length(mesh.tris)) triangles, $(length(mesh.verts)) vertices")

    struts, classes = @timed_stage rl "classify" begin
        ss = candidate_struts(lp, lo, hi)
        cs = Vector{StrutClass}(undef, length(ss))
        max_t = norm3(hi - lo) * 2.0   # safely exceeds any strut-midpoint-to-mesh ray distance
        # Pure Julia over read-only sh/mesh (classify.jl's design note: only
        # tessellate_surface touches gmsh; every classification primitive is
        # plain-data), each iteration writing only its own cs[i] — safe to
        # thread. Cost scales with candidate count (∝ enclosed volume), so
        # this matters increasingly for larger lattices even though it was a
        # minor stage (7s) on the reference run (docs/algorithm.md §11.2).
        Threads.@threads for i in eachindex(ss)
            cs[i] = classify_strut(lp, sh, mesh, dtol, ss[i]; max_t=max_t)
        end
        (ss, cs)
    end
    n_interior = count(==(INTERIOR), classes)
    n_boundary = count(==(BOUNDARY), classes)
    n_outside = count(==(OUTSIDE), classes)
    log_line(rl, "Classification: $(length(struts)) candidates -> interior=$n_interior " *
                  "boundary=$n_boundary outside=$n_outside"; console=true)

    W = determine_worker_count(args)
    local n_tile::Int
    if args.cores !== nothing
        mem_per_strut, probe_seconds, probe_struts = calibrate(lp)
        n_tile, n_mem, n_time = derive_tile_size(mem_per_strut, probe_seconds, probe_struts, W, args.ram)
        log_line(rl, "Calibration: mem_per_strut=$(round(mem_per_strut, digits=1)) bytes, " *
                      "probe=$(probe_struts) struts in $(format_duration(probe_seconds)) -> " *
                      "n_mem=$n_mem n_time=$n_time -> tile_cells=$n_tile (capped at 8)"; console=true)
    else
        n_tile = args.tile_cells
        n_tile > 8 && log_line(rl, "Warning: --tile-cells=$n_tile is above the empirically-measured fuse-time " *
                                    "knee (docs/algorithm.md §7.1, §11.2); consider --cores/--ram auto-tuning " *
                                    "instead, or a value <= 8."; console=true)
    end
    log_line(rl, "Workers: $W, tile edge: $n_tile cells"; console=true)
    if W > 1 && length(workers()) < W
        log_line(rl, "Note: $(length(workers())) Distributed worker(s) active (wanted $W); " *
                      "worker pool is set up by main.jl before run_pipeline is called."; console=true)
    end
    ram_budget_bytes = args.ram === nothing ? nothing : args.ram * 1024.0^3

    tiles = @timed_stage rl "tiling" begin
        # Bucket relative to the candidate range's own minimum corner, not
        # absolute index (0,0,0) — see tile_key's docstring: for geometry
        # centered on its own STEP origin (the common case), bucketing from
        # zero always creates a spurious 2×2×2 minimum tile split through
        # the middle of the part, regardless of how large n_tile is.
        tile_origin = (minimum(s.i for s in struts), minimum(s.j for s in struts), minimum(s.k for s in struts))
        partition_tiles(struts, classes, n_tile; origin=tile_origin)
    end
    fi_keys = TileKey[k for (k, td) in tiles if is_full_interior(td, n_tile)]
    fi_set = Set(fi_keys)
    ref_key = isempty(fi_keys) ? nothing : first(fi_keys)
    all_keys = collect(keys(tiles))
    dispatch_keys = ref_key === nothing ? all_keys : TileKey[k for k in all_keys if k == ref_key || !(k in fi_set)]
    dispatch_list = TileDescriptor[tiles[k] for k in dispatch_keys]
    # Longest-processing-time-first: `dispatch_tiles`' work-stealing pool
    # pulls tiles from `dispatch_list` in order via a shared index counter,
    # so sorting by descending strut count here directly makes it a classic
    # LPT list-scheduling policy. Without this, tile order follows
    # `Dict`-iteration (effectively arbitrary), which let a handful of the
    # largest tiles start late and finish long after every other worker had
    # gone idle (observed directly: the last tile to start in the
    # `test-cylinder-cc5t1` run began 20 minutes after the smallest tiles
    # had already finished, docs/algorithm.md §11.2).
    sort!(dispatch_list; by=td -> -(length(td.interior) + length(td.boundary)))
    log_line(rl, "Tiles: total=$(length(tiles)) full-interior=$(length(fi_keys)) dispatched=$(length(dispatch_list))"; console=true)

    tile_results_vec = @timed_stage rl "tile_stage" begin
        dispatch_tiles(lp, input_brep, dispatch_list, tmpdir, rl, ram_budget_bytes, n_tile, tile_origin)
    end
    tile_results = Dict{TileKey,Any}(r.key => r for r in tile_results_vec)

    n_written, n_removed = with_gmsh() do
        new_model(part_name(args.input, args.cc, args.t))

        final_tags = @timed_stage rl "assembly" begin
            assemble!(lp, tile_results, ref_key, fi_keys, n_tile, tmpdir, rl)
        end
        length(final_tags) == 0 && throw(ProcessingError(
            "Generation produced no geometry at all — the lattice may not intersect the input volume; " *
            "check that -cc/-t are appropriate for the input geometry's scale."))
        assert_no_stray_solids(final_tags)

        @timed_stage rl "export" begin
            threshold = args.t^3
            # filter_floating! only ever deletes solids proven to be
            # genuinely disconnected floating bodies (specification.md §5);
            # a sub-threshold solid that remains connected to other
            # geometry after an attempted fuse is a hard failure, not a
            # silent deletion (docs/algorithm.md §8, §11.2).
            kept, removed, removed_vols = filter_floating!(final_tags, threshold; rl=rl)
            kept_count = length(kept)
            kept_count == 0 && throw(ProcessingError("No solids remain after floating-body cleanup (all below t^3 = $threshold mm^3)"))
            if removed > 0
                total_vol = sum(removed_vols)
                sample = first(sort(removed_vols), min(20, length(removed_vols)))
                log_line(rl, "Removed $removed floating sub-threshold solid(s), total volume " *
                              "$(round(total_vol, digits=4)) mm^3 (min=$(round(minimum(removed_vols), digits=4)), " *
                              "max=$(round(maximum(removed_vols), digits=4)), threshold=$(round(threshold, digits=4)) mm^3); " *
                              "sample volumes: $(round.(sample, digits=4))"; console=true)
            end
            # Fast model-level removal (docs/algorithm.md §8) — must be
            # followed only by a sync-free write (a synchronizing write
            # would re-populate the model from the unmodified OCC kernel
            # state and resurrect the entities filter_floating! just
            # removed from the model-level list).
            write_model(args.output; sync=false)
            (kept_count, removed)
        end
    end
    rewrite_step_header!(args.output, part_name(args.input, args.cc, args.t),
                          generation_params_text(args.input, args.cc, args.t))
    stats[:solids_written] = n_written
    stats[:solids_removed] = n_removed
    stats[:output_path] = args.output
    stats[:tile_count] = length(tiles)
    stats[:interior_tile_count] = length(fi_keys)
    stats[:boundary_tile_count] = length(tiles) - length(fi_keys)
    stats[:worker_count] = W

    @timed_stage rl "verify" begin
        round_trip_check(args.output)
    end

    rm(tmpdir; recursive=true, force=true)
    log_line(rl, "Temp directory cleaned up: $tmpdir")

    return stats
end
