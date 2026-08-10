# src/runlog.jl
#
# Error taxonomy, exit codes, and run logging. Implements docs/algorithm.md §9
# (Logging and failure modes) and the "Exit" requirements of specification.md §3.

using Dates
using Printf

# --- Error taxonomy -> exit codes (docs/algorithm.md §9) --------------------

abstract type LatticeGenError <: Exception end

struct ParamError <: LatticeGenError        # exit 2: parameter validation
    msg::String
end
struct InputGeometryError <: LatticeGenError  # exit 3: input geometry read/parse
    msg::String
end
struct ProcessingError <: LatticeGenError     # exit 4: geometry processing
    msg::String
end
struct ResourceError <: LatticeGenError       # exit 5: resource limits
    msg::String
end
struct OutputError <: LatticeGenError         # exit 6: output write failure
    msg::String
end
struct CancelledError <: LatticeGenError      # exit 130: user cancellation (Ctrl+C)
    msg::String
end

Base.showerror(io::IO, e::LatticeGenError) = print(io, e.msg)

exit_code(::ParamError) = 2
exit_code(::InputGeometryError) = 3
exit_code(::ProcessingError) = 4
exit_code(::ResourceError) = 5
exit_code(::OutputError) = 6
# 128 + SIGINT(2), the POSIX convention for "terminated by Ctrl+C". Not a
# failure of the tool: the run did exactly what the user asked by stopping.
exit_code(::CancelledError) = 130
# Any unforeseen internal exception is bucketed as a geometry-processing
# failure (the pipeline stage where most runtime errors are expected to
# originate) rather than silently becoming a generic crash.
exit_code(::Exception) = 4

# --- Interrupt (Ctrl+C) recognition --------------------------------------

"""
    is_interrupt(e) -> Bool

Is `e` a user Ctrl+C, in any of the wrappers Julia can hand it back in?

A raw `InterruptException` is only what the *directly* interrupted task
sees. By the time one reaches an error handler in this codebase it is
typically wrapped one or more layers deep:

- `TaskFailedException` / `CompositeException` — an interrupt raised inside
  an `@async` task under `@sync` (`dispatch_tiles`' worker pool).
- `RemoteException` wrapping a `CapturedException` — an interrupt that hit a
  `Distributed` worker mid-`remotecall_fetch`.

`RemoteException` belongs to `Distributed`, which this module deliberately
does not depend on (error taxonomy sits below the parallel layer), so it is
matched structurally via its `captured` field rather than by type.

Note that `ProcessExitedException` is deliberately *not* treated as an
interrupt: a worker can also die from an out-of-memory kill or a crash, and
silently reporting those as "cancelled by the user" would hide a real
failure behind a benign exit code.
"""
function is_interrupt(e)
    e isa InterruptException && return true
    e isa CancelledError && return true
    e isa CapturedException && return is_interrupt(e.ex)
    e isa TaskFailedException && return is_interrupt(e.task.result)
    e isa CompositeException && return any(is_interrupt, e.exceptions)
    hasproperty(e, :captured) && return is_interrupt(getproperty(e, :captured))
    return false
end

# --- Formatting helpers -------------------------------------------------

"""
    format_param(x::Real) -> String

Format a numeric CLI parameter without trailing zeros, for use in generated
output file names and STEP metadata (docs/algorithm.md §8, §12), e.g.
`format_param(2.50) == "2.5"`, `format_param(10.0) == "10"`.
"""
function format_param(x::Real)
    s = @sprintf("%.6f", float(x))
    s = rstrip(s, '0')
    s = rstrip(s, '.')
    isempty(s) && return "0"
    return s
end

"""
    format_bytes(b) -> String

Human-readable byte count, switching from MB to GB at 1024 MB.
"""
function format_bytes(b::Integer)
    b < 1024^3 && return @sprintf("%.1f MB", b / 1024^2)
    return @sprintf("%.2f GB", b / 1024^3)
end

"""
    format_duration(seconds) -> String

Human-readable duration, scaling the unit breakdown by magnitude.
"""
function format_duration(seconds::Real)
    seconds < 60 && return @sprintf("%.2f s", seconds)
    m, s = divrem(seconds, 60)
    m < 60 && return @sprintf("%dm %.1fs", Int(m), s)
    h, m = divrem(m, 60)
    return @sprintf("%dh %dm %.1fs", Int(h), Int(m), s)
end

# --- Run log --------------------------------------------------------------

"""
    RunLog

Tracks the log file handle, run start time, per-stage timings, and the
running maximum RSS observed on the master process, for the duration of one
`latticegen2` invocation. `-v`/`--verbose` (the `verbose` field) only affects
what additionally prints to the console — everything is always written to
the log file in full (specification.md §3 "Logging").
"""
mutable struct RunLog
    io::IOStream
    path::String
    verbose::Bool
    start_time::DateTime
    t0::Float64
    stage_times::Vector{Tuple{String,Float64}}
    max_rss::UInt64
end

function open_runlog(path::AbstractString; verbose::Bool=false)
    io = open(path, "w")
    rl = RunLog(io, String(path), verbose, now(), time(), Tuple{String,Float64}[], Sys.maxrss())
    return rl
end

function close_runlog(rl::RunLog)
    flush(rl.io)
    close(rl.io)
end

"""Update the running maximum RSS estimate from the current process."""
function update_rss!(rl::RunLog)
    rl.max_rss = max(rl.max_rss, Sys.maxrss())
    return rl.max_rss
end

"""Record RSS reported by a worker process into the running maximum."""
function observe_rss!(rl::RunLog, worker_rss::Integer)
    rl.max_rss = max(rl.max_rss, UInt64(worker_rss))
    return rl.max_rss
end

_timestamp() = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")

"""
    log_line(rl, msg; console=false)

Always write `msg` to the log file. Also print to console when `console` is
true, or unconditionally when `rl.verbose` is true.
"""
function log_line(rl::RunLog, msg::AbstractString; console::Bool=false)
    line = "[$(_timestamp())] $msg"
    println(rl.io, line)
    flush(rl.io)
    (console || rl.verbose) && println(line)
    update_rss!(rl)
    return nothing
end

"""Log a message unconditionally to both the log file and the console
(used for the run header, stage summaries required by spec, and the final
success/failure report)."""
log_always(rl::RunLog, msg::AbstractString) = log_line(rl, msg; console=true)

"""Log a per-stage timing line and record it for the summary report."""
function log_stage(rl::RunLog, name::AbstractString, elapsed::Real)
    push!(rl.stage_times, (String(name), Float64(elapsed)))
    log_line(rl, "Stage '$name' completed in $(format_duration(elapsed))")
end

"""
    @timed_stage rl "stage name" begin ... end

Times the enclosed block, logs it via `log_stage`, and returns the block's
value.
"""
macro timed_stage(rl, name, block)
    quote
        local _t0 = time()
        local _result = $(esc(block))
        log_stage($(esc(rl)), $(esc(name)), time() - _t0)
        _result
    end
end

# --- Run header and end-of-run summary ------------------------------------

"""Print the run header: input parameters and start time (spec §3 "Exit")."""
function log_run_header(rl::RunLog, params::AbstractDict)
    log_always(rl, "=== latticegen2 run started ===")
    log_always(rl, "Start time: $(Dates.format(rl.start_time, "yyyy-mm-dd HH:MM:SS"))")
    for (k, v) in params
        log_always(rl, "Parameter $k = $v")
    end
end

"""
    log_summary(rl, params, stats)

Print the mandatory end-of-run summary (spec §3 "Exit"): input parameters,
start time, duration, run characteristics, maximum memory usage, and output
path. `stats` is a `Dict{Symbol,Any}` accumulated by the pipeline; missing
keys print as "n/a" so this never throws mid-report.
"""
function log_summary(rl::RunLog, params::AbstractDict, stats::AbstractDict)
    duration = time() - rl.t0
    g(k) = get(stats, k, "n/a")
    log_always(rl, "=== latticegen2 run summary ===")
    log_always(rl, "Start time: $(Dates.format(rl.start_time, "yyyy-mm-dd HH:MM:SS"))")
    log_always(rl, "Duration: $(format_duration(duration))")
    for (k, v) in params
        log_always(rl, "Parameter $k = $v")
    end
    log_always(rl, "Tiles: total=$(g(:tile_count)) interior=$(g(:interior_tile_count)) boundary=$(g(:boundary_tile_count))")
    log_always(rl, "Workers: $(g(:worker_count))")
    log_always(rl, "Stage timings:")
    for (name, elapsed) in rl.stage_times
        log_always(rl, "  - $name: $(format_duration(elapsed))")
    end
    log_always(rl, "Maximum memory usage: $(format_bytes(rl.max_rss))")
    log_always(rl, "Solids written: $(g(:solids_written)) (removed as sub-threshold: $(g(:solids_removed)))")
    log_always(rl, "Output STEP file: $(g(:output_path))")
    log_always(rl, "Log file: $(rl.path)")
end

"""Print the single human-readable failure line required by spec §3."""
function log_failure(rl::RunLog, err::Exception)
    msg = sprint(showerror, err)
    log_always(rl, "FAILED: $msg")
end

"""Print the single human-readable cancellation line for a user Ctrl+C
(exit 130). Distinct from `log_failure` because a cancelled run is not a
malfunction — the log should not read as though something went wrong."""
log_cancelled(rl::RunLog, msg::AbstractString) = log_always(rl, "CANCELLED: $msg")
