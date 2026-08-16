"""Command-line surface (specification.md §3).

``--cores`` and ``--ram`` are optional *budgets* — ceilings on what a run may
use, not hints it may exceed. Both resolve to a concrete number either way:
``--cores`` to the worker count for every parallel stage, defaulting to this
machine's logical cores; ``--ram`` to a memory budget, defaulting to what is
actually free at startup and capped at what the machine physically has. The
detection behind both defaults lives in :mod:`latticegen2.sysinfo`.

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

from . import sysinfo
from .errors import InputGeometryError, OutputError, ParamError
from .lattice import format_param

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
  -o, --output <path>   Output .step path (default: <input_stem>-cc<cc>t<t>.step)
  --cores <n>           Maximum cores to use (1-128), one worker process per
                        core. Default: this machine's logical core count.
  --ram <GB>            Memory budget, from 1 GB up to this machine's total
                        physical RAM. Default: RAM free at startup. Advisory:
                        recorded in the run log next to the measured peak.
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
    ram: float | None
    """Memory budget in GB as the user gave it, or ``None`` if they did not.

    Same distinction as :attr:`cores`, against the resolved
    :attr:`ram_budget_gb`.
    """
    ram_budget_gb: float
    """The memory budget actually in force, in GB — never ``None``.

    Advisory in the sense that nothing enforces it: there is no tile-sizing
    calculation left for it to feed and no memory watchdog, since the
    distributed assembly stage that needed backpressure does not exist in this
    architecture. It is recorded in the run log so a run's measured peak can be
    read against the budget it was given.
    """

    def as_dict(self) -> dict:
        """The parameter block the run header and summary report."""
        return {
            "input": self.input,
            "output": self.output,
            "cc": f"{format_param(self.cc)} mm",
            "t": f"{format_param(self.t)} mm",
            "workers": self.workers,
            "cores": self.cores if self.cores is not None else "auto",
            "ram": (
                f"{self.ram_budget_gb:.1f} GB"
                + ("" if self.ram is not None else " (free at startup)")
            ),
        }


class HelpRequested(Exception):
    """Raised for ``-h``/``--help`` so the caller can exit 0 cleanly."""


def _value(argv: list[str], i: int, flag: str) -> tuple[str, int]:
    if i + 1 >= len(argv):
        raise ParamError(f"Missing value for {flag}")
    return argv[i + 1], i + 2


def _as_float(flag: str, s: str) -> float:
    try:
        return float(s)
    except ValueError:
        raise ParamError(f"Invalid numeric value for {flag}: '{s}'") from None


def _as_int(flag: str, s: str) -> int:
    try:
        return int(s)
    except ValueError:
        raise ParamError(f"Invalid integer value for {flag}: '{s}'") from None


def _in_range(flag: str, v: float, lo: float, hi: float) -> float:
    if v < lo or v > hi:
        raise ParamError(f"{flag} = {v} is out of the valid range [{lo}, {hi}]")
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


def resolve_output_paths(input_path: str, output_arg: str | None, cc: float, t: float):
    """``(step_path, log_path)`` per specification.md §3.

    Default name is ``<input_stem>-cc<cc>t<t>.step`` beside the input. The log
    always shares the output's stem with a ``.log`` extension — never
    ``<output>.step.log``.

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
            directory, f"{stem}-cc{format_param(cc)}t{format_param(t)}.step"
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
        f'<input_stem>-cc<cc>t<t>.step beside the input.'
    )


def parse_args(argv: list[str]) -> Args:
    """Parse and fully validate arguments, raising :class:`ParamError` (exit 2).

    No filesystem or geometry work happens here — see :func:`preflight_checks`
    for the checks that need disk access. specification.md §7: invalid or
    out-of-range parameters are rejected before any computation starts.
    """
    input_path = output = None
    cc = t = ram = cores = None
    verbose = False

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
        elif a == "--ram":
            v, i = _value(argv, i, a)
            ram = _as_float(a, v)
        elif a in ("-v", "--verbose"):
            verbose = True
            i += 1
        elif a in ("-h", "--help"):
            raise HelpRequested()
        else:
            raise ParamError(f"Unknown argument: {a}. Run with --help for usage.")

    if input_path is None:
        raise ParamError("-i/--input is required.\n\n" + USAGE)
    if cc is None:
        raise ParamError("-cc is required.\n\n" + USAGE)
    if t is None:
        raise ParamError("-t is required.\n\n" + USAGE)

    _in_range("-cc", cc, 0.4, 50.0)
    _in_range("-t", t, 0.4, 20.0)
    if cores is not None:
        _in_range("--cores", cores, 1, 128)
    # The ceiling is the machine's own RAM rather than a static literal: a budget
    # above what physically exists is not a budget. The floor stays a plain
    # sanity check. Reported with its own message rather than through
    # `_in_range`, so the ceiling reads as the machine fact it is instead of as
    # a raw float with fifteen digits of detection noise.
    if ram is not None:
        total = sysinfo.total_ram_gb()
        if ram < 1.0 or ram > total:
            raise ParamError(
                f"--ram = {ram} GB is out of the valid range [1.0, {total:.1f}] "
                f"GB; the upper bound is this machine's total physical memory."
            )

    a_edge = cc / (2.0 ** 0.5)
    if t >= a_edge:
        raise ParamError(
            f"-t = {t} mm must be smaller than the cube edge length "
            f"a = cc/sqrt(2) = {a_edge:.6f} mm (derived from -cc = {cc}); a strut "
            f"this thick would not fit within one lattice cell."
        )
    # No further restriction applies. The mid-strut cap faces this generator
    # stitches junctions along stay intact for the whole of `t < a` — see
    # `latticegen2.junction.build_template`, which verifies it geometrically at
    # the run's actual parameters anyway.
    step_path, log_path = resolve_output_paths(input_path, output, cc, t)
    return Args(
        input=input_path,
        output=step_path,
        log_path=log_path,
        cc=cc,
        t=t,
        workers=default_workers(cores),
        verbose=verbose,
        cores=cores,
        ram=ram,
        ram_budget_gb=ram if ram is not None else sysinfo.free_ram_gb(),
    )


def preflight_checks(args: Args) -> None:
    """Filesystem checks that must pass before any geometry work begins."""
    if os.path.isdir(args.output):
        raise ParamError(
            f"-o/--output must name the output .step file, but {args.output} is "
            f"an existing directory. Give a filename inside it, or omit -o to "
            f"write <input_stem>-cc<cc>t<t>.step beside the input."
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
