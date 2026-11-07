"""Command-line surface (specification.md §3).

``--cores`` is an optional *budget* — a ceiling on what a run may use, not a
hint it may exceed. It resolves to a concrete number either way: the worker
count for every parallel stage, defaulting to this machine's logical cores.
The detection behind that default lives in :mod:`latticegen2.sysinfo`.

There used to be a second budget, ``--ram``. It was removed (specification.md
§11): it was validated and recorded but never enforced, and the run's
actual peak lives in the master holding the finished result while ``export``
serialises it — a place no worker-count or tile-size lever can reach — so
enforcing it for real would have meant streaming the STEP write, a
substantially larger change than the flag ever suggested. A budget nothing
enforces is not a budget, so rather than leave a parameter that only ever
printed itself, it was taken out.

Every run executes at below-normal process priority, master and workers alike.
That used to be the opt-in ``-bg`` flag; it is unconditional now, because a run
that leaves the machine usable is what a user wants in every case that was ever
worth choosing between.

A hand-rolled parser is used rather than ``argparse`` so that ``-cc`` and ``-t``
can keep their single-dash spelling, which ``argparse`` would treat as clusters
of short flags.
"""

from __future__ import annotations

import os
from dataclasses import dataclass

import math

from . import orient as _orient
from . import sysinfo
from .errors import InputGeometryError, OutputError, ParamError
from .lattice import THETA, format_name_param, format_param, variant_suffix

#: The ranges specification.md §3's flag table declares, as data.
#:
#: Read by :func:`parse_args` below and by :mod:`latticegen2.gui`, which puts
#: them on its spin controls. Two copies of a bound is one copy too many: the
#: GUI would go on accepting a value the parser had started refusing, and the
#: user would meet the difference as a rejected run rather than as a greyed-out
#: field.
CC_RANGE = (0.4, 50.0)
T_RANGE = (0.4, 20.0)
CORES_RANGE = (1, 128)
ORIENT_RANGE = (-360.0, 360.0)
SUPPORT_THICKNESS_RANGE = (0.2, 20.0)
#: `--support-thickness` when it is not given, mm: the thinnest free-standing
#: support member a run is held to (specification.md §4.6). Fixed rather than
#: derived from `t`, so the window can show it before anything is entered.
SUPPORT_THICKNESS_DEFAULT = 1.0

#: The steepest a support face may overhang, from vertical, in degrees: the
#: lattice rod angle ``THETA`` itself (specification.md §4.6), an expression
#: rather than a literal. Not a parameter.
SUPPORT_OVERHANG_DEG = math.degrees(THETA)
#: The input's surface is supported where it overhangs more than this, from
#: vertical (docs/algorithm.md §14.1). Not a parameter either. Here rather than
#: with the detection code, which needs the geometry kernel, because the window
#: tints its preview with the same bar and must not import one.
SUPPORT_START_DEG = 60.0


def support_grid(a: float, t: float, thickness: float) -> tuple[int, float] | None:
    """``(m, thinnest)``: the support grid division and the thinnest
    free-standing support member it produces, in mm, for cube edge ``a``, strut
    size ``t`` and ``--support-thickness``; ``None`` when even ``m = 1`` is
    thinner than asked (docs/algorithm.md §14.3).

    The rod faces are perpendicular to the other two cube axes and the grid's
    planes run through the rod axes, so a support node beside a rod loses only a
    ``t/2 × t/2`` corner square. At ``m >= 2`` a support node carries a rod on
    at most one side per axis and its thinnest arm is ``q - t/2``, ``q = a/m``;
    at ``m = 1`` it is a whole cube cell with a rod on all twelve edges, and the
    thinnest member is ``a - t``. ``m`` is the largest division still at least
    ``thickness`` thick.

    Here rather than with the generator for the reason the two angles above
    are: the window shows the result before a run, and must not import the
    geometry kernel."""
    m = int(math.floor(a / (t / 2.0 + thickness) + 1e-9))
    if m >= 2:
        return m, a / m - t / 2.0
    if a - t >= thickness - 1e-9:
        return 1, a - t
    return None


def support_thickness_max(a: float, t: float, decimals: int | None = None) -> float:
    """The largest ``--support-thickness`` any division can deliver: ``a - t``,
    at ``m = 1``, which is never less than ``m = 2``'s ``(a - t)/2``.

    With ``decimals``, rounded *down*, so the figure a message names is itself
    accepted: ``a - t`` at ``cc = 5, t = 1`` is 2.5355, and "at most 2.536"
    would name a value the parser refuses."""
    top = a - t
    if decimals is None:
        return top
    scale = 10 ** decimals
    return math.floor(top * scale + 1e-9) / scale


USAGE = """\
latticegen2 - generate a diamond-strut lattice filling an input STEP volume

Usage:
  latticegen2 -i <input.step> -cc <mm> -t <mm> [options]

Required:
  -i, --input <path>    Path to the STEP file defining the lattice bounds
  -cc <mm>              Distance between bottom nodes of adjacent cells (0.4-50)
  -t <mm>               Side length of the diamond strut profile (0.4-20),
                        must be smaller than the cell edge a = cc/sqrt(2)

Optional:
  -o, --output <path>   Output .step path
                        (default: <input_stem>-lattice-cc<cc>t<t>.step)
  --cores <n>           Maximum cores to use (1-128), one worker process per
                        core. Default: this machine's logical core count.
  --orient <rx> <ry> <rz>
                        Print orientation of the input geometry relative to the
                        print bed, in degrees (-360-360 each, default 0 0 0):
                        rotations about the bed's X, then Y, then Z axis. The
                        lattice is oriented to the bed.
  --support             Generate print supports for inner overhangs: the
                        surface overhanging more than 60 deg from vertical,
                        with supports no steeper than the lattice rod angle
  --support-thickness <mm>
                        Thinnest free-standing support member (0.2-20, at
                        most a - t, default 1). Sets how finely the support
                        grid divides the lattice cube. Requires --support.
  -v, --verbose         Verbose console output (a full .log is always written)
  -h, --help            Show this message and exit

Every run executes at below-normal process priority so the machine stays
usable for other work.
"""


@dataclass(frozen=True)
class Args:
    """Parsed and fully validated arguments, with output paths resolved."""

    input: str
    output: str
    log_path: str
    cc: float
    t: float
    workers: int
    verbose: bool
    cores: int | None
    """Core budget as the user gave it, or ``None`` if they did not.

    Kept alongside the resolved :attr:`workers` rather than replaced by it, so
    the log can distinguish a run that was *told* to use six cores from one that
    happened to detect six.
    """

    progress_stream: bool = False
    """Report the run as machine-readable events on stdout (``--progress-stream``).

    An internal transport, set by :mod:`latticegen2.gui` on the child it
    launches, and deliberately **absent from** :data:`USAGE` and from
    specification.md §3's flag table: those describe the parameters a person
    chooses, and this one changes no parameter of the geometry. It is documented
    in specification.md §3.1 instead, with the rest of the front-end.
    """

    orient: tuple[float, float, float] = _orient.IDENTITY
    """Print orientation ``(rx, ry, rz)`` in degrees (specification.md §4.5)."""
    support: bool = False
    support_thickness: float | None = None
    """Minimum support thickness in mm; ``None`` means the strut size ``t``."""

    def __post_init__(self):
        if self.support_thickness is None:
            object.__setattr__(self, "support_thickness", self.t)

    def params_extra(self) -> str:
        """Orientation and support parameters for the STEP header, when in force."""
        out = ""
        if _orient.is_rotated(self.orient):
            out += ("orient=" + " ".join(format_param(a) for a in self.orient)
                    + " deg, ")
        if self.support:
            out += f"support thickness={format_param(self.support_thickness)} mm, "
        return out

    def as_dict(self) -> dict:
        """The parameter block the run header and summary report.

        ``progress_stream`` is deliberately not among them: it changes how the
        run is *watched*, not what it computes, and a log whose parameter block
        differed depending on whether a window was open would be misleading in
        exactly the place operators compare two runs.
        """
        out = {
            "input": self.input,
            "output": self.output,
            "cc": f"{format_param(self.cc)} mm",
            "t": f"{format_param(self.t)} mm",
            "workers": self.workers,
            "cores": self.cores if self.cores is not None else "auto",
        }
        # Added only when in force, so an unrotated, unsupported run's parameter
        # block is the one it has always had.
        if _orient.is_rotated(self.orient):
            out["orient"] = " ".join(format_param(a) for a in self.orient) + " deg"
        if self.support:
            out["support"] = "enabled"
            out["support_thickness"] = f"{format_param(self.support_thickness)} mm"
        return out


class HelpRequested(Exception):
    """Raised for ``-h``/``--help`` so the caller can exit 0 cleanly."""


def _value(argv: list[str], i: int, flag: str) -> tuple[str, int]:
    if i + 1 >= len(argv):
        raise ParamError(f"Missing value for {flag}")
    return argv[i + 1], i + 2


#: The ``ParamError.fix`` field naming each flag's parameter.
_FIX_FIELD = {
    "-cc": "cc", "-t": "t", "--cores": "cores", "--support-thickness": "thickness",
    "--orient": "orient", "--orient rx": "rx", "--orient ry": "ry", "--orient rz": "rz",
}
#: The unit a range in a fix is given in.
_FIX_UNIT = {"-cc": " mm", "-t": " mm", "--support-thickness": " mm",
             "--orient rx": "°", "--orient ry": "°", "--orient rz": "°"}


def _field(flag: str) -> str:
    return "{" + _FIX_FIELD[flag] + "}" if flag in _FIX_FIELD else flag


def _round_up(v: float, decimals: int = 2) -> float:
    """``v`` rounded up, so a bound named as "at least" is itself enough."""
    step = 10.0 ** decimals
    return math.ceil(v * step - 1e-9) / step


def _round_down(v: float, decimals: int = 2) -> float:
    """``v`` rounded down, so a bound named as "smaller than" is never passed."""
    step = 10.0 ** decimals
    return math.floor(v * step + 1e-9) / step


def _thickness_fix(cc: float, t: float, thickness: float, a: float) -> str:
    """What makes a support thickness above ``a - t`` acceptable: a thinner
    value, a larger cell, or a thinner strut -- each offered only where its
    own range allows it."""
    options = []
    most = support_thickness_max(a, t, decimals=2)
    if most >= SUPPORT_THICKNESS_RANGE[0]:
        options.append(f"lower {{thickness}} to at most {most:.2f} mm")
    cc_needed = _round_up((thickness + t) * 2.0 ** 0.5)
    if cc_needed <= CC_RANGE[1]:
        options.append(f"raise {{cc}} to at least {cc_needed:.2f} mm")
    t_most = _round_down(a - thickness)
    if not options and t_most >= T_RANGE[0]:
        options.append(f"lower {{t}} to at most {t_most:.2f} mm")
    if not options:
        return "lower {thickness} or {t}, or raise {cc}"
    return ", or ".join(options)


def _as_float(flag: str, s: str) -> float:
    try:
        return float(s)
    except ValueError:
        raise ParamError(f"Invalid numeric value for {flag}: '{s}'",
                         fix=f"enter a number for {_field(flag)}") from None


def _as_int(flag: str, s: str) -> int:
    try:
        return int(s)
    except ValueError:
        raise ParamError(f"Invalid integer value for {flag}: '{s}'",
                         fix=f"enter a whole number for {_field(flag)}") from None


def _in_range(flag: str, v: float, lo: float, hi: float) -> float:
    if v < lo or v > hi:
        raise ParamError(f"{flag} = {v} is out of the valid range [{lo}, {hi}]",
                         fix=f"set {_field(flag)} between {format_param(lo)} and "
                             f"{format_param(hi)}{_FIX_UNIT.get(flag, '')}")
    return v


def default_workers(cores: int | None) -> int:
    """Worker count from the core budget, else from this machine.

    **One worker per core, uncapped.** ``--cores`` is a budget the user set
    deliberately, so it is honoured exactly rather than second-guessed; with no
    ``--cores`` the budget is the machine's own logical core count
    (:func:`latticegen2.sysinfo.logical_core_count`).

    The master does not need a core reserved for it: for all but a percent or so
    of the boundary stage it is blocked in ``Pool.imap`` waiting on results, and
    the work it does between batches — deserialising each batch's ``.brep`` — is
    serial whether or not a core is held back for it. Measured on the scale
    rehearsal's own output: reading all 151 MB of its boundary ``.brep`` files
    takes 7.0 s against a 12 m 36 s boundary stage. Reserving one cost a fifth of
    this stage's throughput on a 6-core machine to protect an idle process.
    Desktop impact is handled instead by every process running at below-normal
    priority, which is now unconditional.

    Memory scales linearly with this number — 260 MB peak per worker at rehearsal
    scale — so it stays a small fraction of the master's own footprint.

    An earlier revision capped this at 8, on the empirical ground that boundary
    jobs are short enough for the pool's start-up and result-marshalling cost to
    start dominating past that. The cap is gone: it silently contradicted an
    explicit ``--cores``, and a budget the user has to guess is being ignored is
    worse than one that costs a little throughput on a very wide machine.
    """
    if cores is None:
        cores = sysinfo.logical_core_count()
    return max(1, cores)


def resolve_output_paths(input_path: str, output_arg: str | None, cc: float, t: float,
                         orient=_orient.IDENTITY, support: bool = False,
                         thickness: float | None = None):
    """``(step_path, log_path)`` per specification.md §3.

    ``orient`` and ``support`` add ``-rotx<rx>y<ry>z<rz>`` and ``-sup<thickness>``
    to the default name, each only when it applies
    (:func:`latticegen2.lattice.variant_suffix`).

    Default name is ``<input_stem>-lattice-cc<cc>t<t>.step`` beside the input.
    The ``-lattice-`` marks the file as this tool's output rather than another
    revision of the input, which matters because the default writes it into the
    input's own folder. The log always shares the output's stem with a ``.log``
    extension — never ``<output>.step.log``.

    ``-o`` names a *file*. A directory is rejected rather than quietly turned
    into one: ``-o .\\`` used to append the extension to nothing and write
    ``.\\.step`` beside a ``.\\.step.log``, because a basename that is all
    extension has no stem for the log to share. Silently producing two files
    under names the user never asked for is worse than refusing the argument
    (specification.md §7).
    """
    if output_arg is None:
        stem = os.path.splitext(os.path.basename(input_path))[0]
        directory = os.path.dirname(input_path) or "."
        step_path = os.path.join(
            directory,
            f"{stem}-lattice-cc{format_name_param(cc)}t{format_name_param(t)}"
            f"{variant_suffix(orient, support, thickness=thickness)}.step",
        )
    else:
        # The validated string is the one used. Trailing space is stripped
        # before both, or `-o "out.step "` would be approved as `out.step` and
        # then fail its own extension test, landing as `out.step .step`.
        step_path = _checked_output(output_arg)
        if not step_path.lower().endswith(".step"):
            step_path += ".step"
    return step_path, os.path.splitext(step_path)[0] + ".log"


def _checked_output(output_arg: str) -> str:
    """The ``-o`` argument, or :class:`ParamError` if it cannot name a file.

    Returns the string the caller should actually use, so the value that was
    validated and the value that gets built on are the same one.

    Two pure-syntax rules, so this stays free of the filesystem like the rest of
    parsing. A trailing separator (``.\\``, ``out/``) says "directory" on its
    face whatever is on disk — and *which characters those are is
    platform-dependent*: a backslash separates on Windows and is an ordinary
    filename character on POSIX, so the rule follows ``os.sep``/``os.altsep``
    rather than assuming. A basename that is empty or all dots (``.``, ``..``)
    has no stem to build a name from — appending the extension yields
    ``..step``, and the log then has no stem to share, so both files land under
    names the user never typed. Any real filename has at least one non-dot
    character in its basename, so neither rule can refuse a valid argument.

    An ``-o`` whose *resolved* path is an existing directory is caught by
    :func:`preflight_checks` instead, which is where filesystem checks belong.
    """
    trimmed = output_arg.rstrip()
    separators = (os.sep, os.altsep) if os.altsep else (os.sep,)
    if trimmed and not trimmed.endswith(separators):
        if os.path.basename(trimmed).strip("."):
            return trimmed
    # Quoted rather than repr()'d: on Windows the argument is full of
    # backslashes, and repr doubling them reads as if the tool mangled the
    # input. The quotes still make a whitespace-only argument visible.
    raise ParamError(
        f'-o/--output must name the output .step file, not a directory '
        f'(got "{output_arg}"). Give a filename, or omit -o to write '
        f'<input_stem>-lattice-cc<cc>t<t>.step beside the input.'
    )


def parse_args(argv: list[str]) -> Args:
    """Parse and fully validate arguments, raising :class:`ParamError` (exit 2).

    No filesystem or geometry work happens here — see :func:`preflight_checks`
    for the checks that need disk access. specification.md §7: invalid or
    out-of-range parameters are rejected before any computation starts.
    """
    input_path = output = None
    cc = t = cores = None
    verbose = False
    progress_stream = False
    orient = None
    support = False
    thickness = None

    i = 0
    while i < len(argv):
        a = argv[i]
        if a in ("-i", "--input"):
            input_path, i = _value(argv, i, a)
        elif a in ("-o", "--output"):
            output, i = _value(argv, i, a)
        elif a == "-cc":
            v, i = _value(argv, i, a)
            cc = _as_float(a, v)
        elif a == "-t":
            v, i = _value(argv, i, a)
            t = _as_float(a, v)
        elif a == "--cores":
            v, i = _value(argv, i, a)
            cores = _as_int(a, v)
        elif a == "--orient":
            if i + 3 >= len(argv):
                raise ParamError(
                    "--orient needs three values: <rx> <ry> <rz> in degrees"
                )
            orient = tuple(_as_float(a, v) for v in argv[i + 1:i + 4])
            i += 4
        elif a == "--support":
            support = True
            i += 1
        elif a == "--overhang":
            # Removed with the support algorithm (specification.md §3): named
            # rather than reported as unknown, so a script written against the
            # skeleton learns why.
            raise ParamError(
                "--overhang has been removed: support faces are held to the "
                "lattice rod angle and the input is supported past 60 deg from "
                "vertical; neither is a parameter (specification.md §4.6)."
            )
        elif a == "--support-thickness":
            v, i = _value(argv, i, a)
            thickness = _as_float(a, v)
        elif a in ("-v", "--verbose"):
            verbose = True
            i += 1
        elif a == "--progress-stream":
            progress_stream = True
            i += 1
        elif a in ("-h", "--help"):
            raise HelpRequested()
        else:
            raise ParamError(f"Unknown argument: {a}. Run with --help for usage.")

    if input_path is None:
        raise ParamError("-i/--input is required.\n\n" + USAGE, fix="choose an {input}")
    if cc is None:
        raise ParamError("-cc is required.\n\n" + USAGE, fix="enter {cc}")
    if t is None:
        raise ParamError("-t is required.\n\n" + USAGE, fix="enter {t}")

    _in_range("-cc", cc, *CC_RANGE)
    _in_range("-t", t, *T_RANGE)
    if cores is not None:
        _in_range("--cores", cores, *CORES_RANGE)

    a_edge = cc / (2.0 ** 0.5)
    if t >= a_edge:
        raise ParamError(
            f"-t = {t} mm must be smaller than the cube edge length "
            f"a = cc/sqrt(2) = {a_edge:.6f} mm (derived from -cc = {cc}); a strut "
            f"this thick would not fit within one lattice cell.",
            fix=f"make {{t}} smaller than {_round_down(a_edge):.2f} mm, or {{cc}} "
                f"larger than {_round_up(t * 2.0 ** 0.5):.2f} mm",
        )
    # No further restriction applies. The mid-strut cap faces this generator
    # stitches junctions along stay intact for the whole of `t < a` — see
    # `latticegen2.junction.build_template`, which verifies it geometrically at
    # the run's actual parameters anyway.
    if orient is None:
        orient = _orient.IDENTITY
    for axis, angle in zip("xyz", orient):
        _in_range(f"--orient r{axis}", angle, *ORIENT_RANGE)
    if not support and thickness is not None:
        # Refused rather than ignored: a support parameter that silently does
        # nothing would leave the user believing supports had been generated.
        raise ParamError("--support-thickness requires --support.")
    if thickness is None:
        thickness = SUPPORT_THICKNESS_DEFAULT
    _in_range("--support-thickness", thickness, *SUPPORT_THICKNESS_RANGE)
    if support and support_grid(a_edge, t, thickness) is None:
        # Asked only of a supported run: without `--support` the value is the
        # default and reaches nothing, so it must not refuse a lattice.
        raise ParamError(
            f"--support-thickness {format_param(thickness)} mm cannot be met at "
            f"cc={format_param(cc)}, t={format_param(t)}: the thinnest support "
            f"member is at most a - t = "
            f"{format_param(support_thickness_max(a_edge, t, decimals=3))} mm "
            f"(a = cc/sqrt(2)).",
            fix=_thickness_fix(cc, t, thickness, a_edge),
        )

    step_path, log_path = resolve_output_paths(
        input_path, output, cc, t, orient, support, thickness
    )
    return Args(
        orient=orient,
        support=support,
        support_thickness=thickness,
        input=input_path,
        output=step_path,
        log_path=log_path,
        cc=cc,
        t=t,
        workers=default_workers(cores),
        verbose=verbose,
        cores=cores,
        progress_stream=progress_stream,
    )


def preflight_checks(args: Args) -> None:
    """Filesystem checks that must pass before any geometry work begins."""
    if os.path.isdir(args.output):
        raise ParamError(
            f"-o/--output must name the output .step file, but {args.output} is "
            f"an existing directory. Give a filename inside it, or omit -o to "
            f"write <input_stem>-lattice-cc<cc>t<t>.step beside the input."
        )
    if not os.path.isfile(args.input):
        raise InputGeometryError(f"Input STEP file not found: {args.input}")
    try:
        with open(args.input, "rb") as fh:
            fh.read(4)
    except OSError as exc:
        raise InputGeometryError(f"Input STEP file is not readable: {args.input} ({exc})")

    outdir = os.path.dirname(args.output) or "."
    if not os.path.isdir(outdir):
        raise OutputError(f"Output directory does not exist: {outdir}")
    probe = os.path.join(outdir, f".latticegen2_write_probe_{os.getpid()}")
    try:
        with open(probe, "wb") as fh:
            fh.write(b"\0")
        os.remove(probe)
    except OSError as exc:
        raise OutputError(f"Output directory is not writable: {outdir} ({exc})")
