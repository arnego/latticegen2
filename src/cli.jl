# src/cli.jl
#
# Hand-rolled CLI parser and validator for specification.md §3 (as amended by
# the additions documented there for --cores/--ram/--workers/--tile-cells).
# No third-party argument-parsing dependency is used, in keeping with
# CLAUDE.md's minimal-dependency guidance for this small, fixed surface.

"""
    CliArgs

Fully parsed and validated command-line arguments, with the output STEP path
and log path already resolved to their final values (specification.md §3).
"""
struct CliArgs
    input::String
    output::String          # resolved final .step path
    log_path::String        # resolved final .log path
    cc::Float64
    t::Float64
    background::Bool
    verbose::Bool
    cores::Union{Int,Nothing}
    ram::Union{Float64,Nothing}
    workers::Union{Int,Nothing}
    tile_cells::Union{Int,Nothing}
end

const USAGE = """
latticegen2 - generate a diamond-strut lattice filling an input STEP volume

Usage:
  julia --project=. src/main.jl -i <input.step> -cc <mm> -t <mm> [options]

Required:
  -i, --input <path>       Path to the STEP file defining the lattice bounds
  -cc <mm>                 Distance between bottom nodes of adjacent cells (0.4-50)
  -t <mm>                  Side length of the diamond rod profile (0.4-20)

  One of the following two optimization-parameter pairs is also required:
    --cores <n> --ram <GB>            Auto-tune workers/tile size from hardware
    --workers <n> --tile-cells <n>    Explicit optimization parameters

Optional:
  -o, --output <path>      Output .step path (default: <input_stem>-cc<cc>t<t>.step)
  -bg, --background        Run worker processes at below-normal priority
  -v, --verbose            Verbose console diagnostics (a full .log is always written)
  -h, --help                Show this message and exit
"""

_print_usage() = println(USAGE)

function _take_value(argv::Vector{String}, i::Int, flag::String)
    if i == length(argv)
        throw(ParamError("Missing value for $flag"))
    end
    return argv[i+1], i + 2
end

function _parse_float(flag::String, s::String)
    v = tryparse(Float64, s)
    v === nothing && throw(ParamError("Invalid numeric value for $flag: '$s'"))
    return v
end

function _parse_int(flag::String, s::String)
    v = tryparse(Int, s)
    v === nothing && throw(ParamError("Invalid integer value for $flag: '$s'"))
    return v
end

function _check_range(flag::String, v::Real, lo::Real, hi::Real)
    (v < lo || v > hi) && throw(ParamError("$flag = $v is out of the valid range [$lo, $hi]"))
    return v
end

"""
    parse_args(argv::Vector{String}) -> CliArgs

Parse and fully validate CLI arguments, throwing a `ParamError` (exit code 2)
on any malformed, missing, out-of-range, or conflicting argument. No
filesystem or geometry work happens here — see `preflight_checks` for the
checks that require touching disk (specification.md §7: "Invalid/out-of-range
parameters should be rejected before any computation starts").
"""
function parse_args(argv::Vector{String})::CliArgs
    input = nothing
    output = nothing
    cc = nothing
    t = nothing
    background = false
    verbose = false
    cores = nothing
    ram = nothing
    workers = nothing
    tile_cells = nothing

    i = 1
    n = length(argv)
    while i <= n
        a = argv[i]
        if a in ("-i", "--input")
            input, i = _take_value(argv, i, a)
        elseif a in ("-o", "--output")
            output, i = _take_value(argv, i, a)
        elseif a == "-cc"
            v, i = _take_value(argv, i, a)
            cc = _parse_float(a, v)
        elseif a == "-t"
            v, i = _take_value(argv, i, a)
            t = _parse_float(a, v)
        elseif a in ("-bg", "--background")
            background = true
            i += 1
        elseif a in ("-v", "--verbose")
            verbose = true
            i += 1
        elseif a == "--cores"
            v, i = _take_value(argv, i, a)
            cores = _parse_int(a, v)
        elseif a == "--ram"
            v, i = _take_value(argv, i, a)
            ram = _parse_float(a, v)
        elseif a == "--workers"
            v, i = _take_value(argv, i, a)
            workers = _parse_int(a, v)
        elseif a == "--tile-cells"
            v, i = _take_value(argv, i, a)
            tile_cells = _parse_int(a, v)
        elseif a in ("-h", "--help")
            _print_usage()
            exit(0)
        else
            throw(ParamError("Unknown argument: $a. Run with --help for usage."))
        end
    end

    input === nothing && throw(ParamError("-i/--input is required. " * USAGE))
    cc === nothing && throw(ParamError("-cc is required. " * USAGE))
    t === nothing && throw(ParamError("-t is required. " * USAGE))

    _check_range("-cc", cc, 0.4, 50.0)
    _check_range("-t", t, 0.4, 20.0)

    a_edge = cc / sqrt(2)
    if t >= a_edge
        throw(ParamError(
            "-t = $t mm must be smaller than the cube edge length a = cc/√2 = " *
            "$(round(a_edge; digits=6)) mm (derived from -cc = $cc); a strut this " *
            "thick would not fit within one lattice cell."
        ))
    end

    has_auto = (cores !== nothing) || (ram !== nothing)
    has_explicit = (workers !== nothing) || (tile_cells !== nothing)

    if has_auto && has_explicit
        throw(ParamError(
            "--cores/--ram and --workers/--tile-cells are mutually exclusive " *
            "(specification.md §3); choose exactly one optimization-parameter path."
        ))
    elseif has_auto
        (cores === nothing || ram === nothing) && throw(ParamError(
            "--cores and --ram must both be provided together for auto-tuning."
        ))
        _check_range("--cores", cores, 1, 128)
        _check_range("--ram", ram, 1.0, 1024.0)
    elseif has_explicit
        (workers === nothing || tile_cells === nothing) && throw(ParamError(
            "--workers and --tile-cells must both be provided together when not auto-tuning."
        ))
        _check_range("--workers", workers, 1, 128)
        _check_range("--tile-cells", tile_cells, 2, 64)
    else
        throw(ParamError(
            "Either both --cores and --ram, or both --workers and --tile-cells, " *
            "must be provided (specification.md §2)."
        ))
    end

    output_step, log_path = resolve_output_paths(input, output, cc, t)

    return CliArgs(input, output_step, log_path, cc, t, background, verbose,
                   cores, ram, workers, tile_cells)
end

"""
    resolve_output_paths(input, output_arg, cc, t) -> (step_path, log_path)

Default output naming per specification.md §3: `<input_stem>-cc<cc>t<t>.step`
in the input file's own directory when `-o` is not given. The log file always
shares the resolved output's stem with a `.log` extension, never
`<output>.step.log`.
"""
function resolve_output_paths(input::AbstractString, output_arg::Union{AbstractString,Nothing},
                               cc::Real, t::Real)
    local step_path::String
    if output_arg === nothing
        stem = splitext(basename(input))[1]
        dir = dirname(input)
        dir = isempty(dir) ? "." : dir
        step_path = joinpath(dir, "$(stem)-cc$(format_param(cc))t$(format_param(t)).step")
    else
        step_path = String(output_arg)
        if !endswith(lowercase(step_path), ".step")
            step_path = step_path * ".step"
        end
    end
    log_path = step_path[1:end-5] * ".log"   # ".step"/".STEP" are both 5 chars
    return step_path, log_path
end

"""
    preflight_checks(args::CliArgs)

Filesystem checks that must pass before any geometry computation begins, but
that are not pure parameter validation (so they raise `InputGeometryError` /
`OutputError`, not `ParamError`): input file existence/readability, and
output directory existence/writability (specification.md §7).
"""
function preflight_checks(args::CliArgs)
    isfile(args.input) || throw(InputGeometryError("Input STEP file not found: $(args.input)"))
    try
        open(args.input, "r") do io
            read(io, 4)
        end
    catch e
        throw(InputGeometryError("Input STEP file is not readable: $(args.input) ($(sprint(showerror, e)))"))
    end

    outdir = dirname(args.output)
    outdir = isempty(outdir) ? "." : outdir
    isdir(outdir) || throw(OutputError("Output directory does not exist: $outdir"))

    probe = joinpath(outdir, ".latticegen2_write_probe_$(getpid())")
    try
        open(probe, "w") do io
            write(io, 0x00)
        end
        rm(probe; force=true)
    catch e
        throw(OutputError("Output directory is not writable: $outdir ($(sprint(showerror, e)))"))
    end
    return nothing
end
