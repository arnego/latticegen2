# src/main.jl
#
# Entry point: `julia --project=. src/main.jl <args>`. Implements the exit
# code / logging / worker-pool-lifecycle surface required by
# specification.md §3 and docs/algorithm.md §9, §12.
#
# IMPORTANT: `Distributed`'s `@everywhere` can only run at true top level
# (never inside a function body — see docs/algorithm.md §7.1 /
# pipeline.jl's `run_pipeline` docstring), which is why the worker-pool
# setup below is written as plain top-level script code rather than
# wrapped in a function.

include(joinpath(@__DIR__, "latticegen2.jl"))
using .latticegen2
using Distributed

# A Julia script run without `-i` defaults to `exit_on_sigint(true)`: Ctrl+C
# terminates the process on the spot, with no `InterruptException` thrown, so
# no `finally` runs and the worker pool is orphaned — leaving half-torn-down
# gmsh sessions and the worker-lock errors that come with them. Turning it off
# routes Ctrl+C through the normal exception path instead, so the handlers
# below (and in `dispatch_tiles`/`process_tile`) can shut the run down in
# order. Worker processes are started as `julia --worker`, not as a program,
# so they already have this off by default and need no equivalent call.
Base.exit_on_sigint(false)

# --- 1. Parse + validate CLI arguments (exit 2). No log file exists yet:
#        specification.md §7 requires rejecting bad parameters "before any
#        computation starts", and we don't yet know a valid output/log path
#        to write one to if this step itself fails. -----------------------
try
    global parsed_args = parse_args(ARGS)
catch e
    if is_interrupt(e)
        println(stderr, "CANCELLED: run cancelled by user (Ctrl+C).")
        exit(130)
    end
    if e isa LatticeGenError
        println(stderr, "ERROR: ", sprint(showerror, e))
        exit(exit_code(e))
    end
    rethrow()
end

# --- 2. Open the run log. From here on every failure (preflight or
#        pipeline) is routed into it, per specification.md §3 "Logging". --
rl = open_runlog(parsed_args.log_path; verbose=parsed_args.verbose)

# --- 3. Worker pool setup (top level only — see module docstring above). -
W = determine_worker_count(parsed_args)
workers_added = false
if W > 1
    try
        project_dir = dirname(Base.active_project())
        addprocs(W; exeflags=["--project=$project_dir"])
        global workers_added = true
        # Target only the newly-added worker processes, NOT the master
        # (pid 1, the default @everywhere target): the master already
        # loaded and `using`'d latticegen2 above. Re-including it there too
        # creates a second, distinct `latticegen2` module object bound to
        # the same name, which then makes every exported identifier
        # ambiguous ("both latticegen2 and latticegen2 export ...") and
        # breaks the master's own already-working bindings.
        @everywhere workers() project_dir = $project_dir
        @everywhere workers() include(joinpath(Main.project_dir, "src", "latticegen2.jl"))
        @everywhere workers() using .latticegen2
        parsed_args.background && @everywhere workers() latticegen2.set_below_normal_priority!()
    catch e
        if is_interrupt(e)
            # Ctrl+C during pool setup: nothing has been computed yet, so
            # tear the partial pool down and stop rather than quietly
            # continuing into a full run the user just asked to abort.
            msg = "Run cancelled by user (Ctrl+C) during worker startup."
            log_cancelled(rl, msg)
            println(stderr, "CANCELLED: ", msg)
            shutdown_workers!()
            close_runlog(rl)
            exit(130)
        end
        log_line(rl, "Distributed worker startup failed ($(sprint(showerror, e))); falling back to sequential."; console=true)
        shutdown_workers!()
        global workers_added = false
    end
end

# --- 4. Run header, then preflight + pipeline, mapping exceptions to the
#        exit codes from docs/algorithm.md §9. -----------------------------
params = Dict(
    "input" => parsed_args.input, "output" => parsed_args.output, "cc" => parsed_args.cc,
    "t" => parsed_args.t, "background" => parsed_args.background, "verbose" => parsed_args.verbose,
    "cores" => parsed_args.cores, "ram" => parsed_args.ram,
    "workers" => parsed_args.workers, "tile_cells" => parsed_args.tile_cells,
)
log_run_header(rl, params)

code = 0
try
    preflight_checks(parsed_args)
    stats = run_pipeline(parsed_args, rl)
    log_summary(rl, params, stats)
    global code = 0
catch e
    if is_interrupt(e)
        # User cancellation (exit 130), not a malfunction: report it as such,
        # with no stacktrace. Whatever `temp/<timestamp>/` the run created is
        # left in place (specification.md §4.4) — its path was logged when the
        # pipeline started, and `CancelledError` from the tile stage repeats it.
        msg = e isa CancelledError ? sprint(showerror, e) :
              "Run cancelled by user (Ctrl+C); any temp directory from this run is left in place for analysis."
        log_cancelled(rl, msg)
        println(stderr, "CANCELLED: ", msg)
        global code = 130
    else
        log_failure(rl, e)
        println(stderr, "ERROR: ", sprint(showerror, e))
        if !(e isa LatticeGenError)
            log_line(rl, "Internal error stacktrace:\n" * sprint(Base.show_backtrace, catch_backtrace()))
        end
        global code = e isa LatticeGenError ? exit_code(e) : 4
    end
finally
    # Worker teardown happens here, before the log closes, so it is recorded
    # and so it runs on *every* path out of the pipeline — success, failure,
    # or Ctrl+C. `shutdown_workers!` never throws (see its docstring): a
    # cleanup error must not replace the run's real outcome.
    if workers_added
        log_line(rl, "Shutting down worker processes...")
        shutdown_workers!()
        log_line(rl, "Worker processes terminated.")
    end
    close_runlog(rl)
end

exit(code)
