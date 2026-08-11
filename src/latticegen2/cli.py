"""Command-line surface (specification.md §3).

``--workers``, ``--cores`` and ``--ram`` are optional hints. Boundary-junction
jobs are constant-size and independent, so a sensible worker count follows from
the machine; an explicit ``--workers`` overrides everything.

A hand-rolled parser is used rather than ``argparse`` so that ``-cc`` and ``-t``
can keep their single-dash spelling, which ``argparse`` would treat as clusters
of short flags.
"""

from __future__ import annotations

import os
from dataclasses import dataclass

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
  --workers <n>         Worker processes for the boundary stage (1-128).
                        Default: derived from --cores, else from this machine.
  --cores <n>           Physical cores available, used to derive --workers (1-128)
  --ram <GB>            Memory budget (1-1024). Advisory: recorded in the run
                        log next to the run's measured peak memory.
  -bg, --background     Run at below-normal process priority
  -v, --verbose         Verbose console output (a full .log is always written)
  -h, --help            Show this message and exit
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
    background: bool
    verbose: bool
    cores: int | None
    ram: float | None
    """Memory budget in GB, or ``None``.

    Advisory only. There is no tile-sizing calculation left for it to feed and no
    memory watchdog — the distributed assembly stage that needed backpressure does
    not exist in this architecture. It is recorded in the run log so a run's
    measured peak can be read against the budget it was given.
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
            "ram": f"{self.ram} GB" if self.ram is not None else "unspecified",
            "background": self.background,
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
    """Worker count from an explicit core count, else from this machine.

    One core is left to the master process and the desktop. The cap at 8 is
    empirical: boundary jobs are short, so past that the pool's own start-up and
    result-marshalling cost starts to dominate the work being distributed.
    """
    if cores is None:
        cores = os.cpu_count() or 2
    return max(1, min(cores - 1, 8))


def resolve_output_paths(input_path: str, output_arg: str | None, cc: float, t: float):
    """``(step_path, log_path)`` per specification.md §3.

    Default name is ``<input_stem>-cc<cc>t<t>.step`` beside the input. The log
    always shares the output's stem with a ``.log`` extension — never
    ``<output>.step.log``.
    """
    if output_arg is None:
        stem = os.path.splitext(os.path.basename(input_path))[0]
        directory = os.path.dirname(input_path) or "."
        step_path = os.path.join(
            directory, f"{stem}-cc{format_param(cc)}t{format_param(t)}.step"
        )
    else:
        step_path = output_arg
        if not step_path.lower().endswith(".step"):
            step_path += ".step"
    return step_path, os.path.splitext(step_path)[0] + ".log"


def parse_args(argv: list[str]) -> Args:
    """Parse and fully validate arguments, raising :class:`ParamError` (exit 2).

    No filesystem or geometry work happens here — see :func:`preflight_checks`
    for the checks that need disk access. specification.md §7: invalid or
    out-of-range parameters are rejected before any computation starts.
    """
    input_path = output = None
    cc = t = ram = None
    workers = cores = None
    background = verbose = False

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
        elif a == "--workers":
            v, i = _value(argv, i, a)
            workers = _as_int(a, v)
        elif a == "--cores":
            v, i = _value(argv, i, a)
            cores = _as_int(a, v)
        elif a == "--ram":
            v, i = _value(argv, i, a)
            ram = _as_float(a, v)
        elif a in ("-bg", "--background"):
            background = True
            i += 1
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
    if workers is not None:
        _in_range("--workers", workers, 1, 128)
    if cores is not None:
        _in_range("--cores", cores, 1, 128)
    if ram is not None:
        _in_range("--ram", ram, 1.0, 1024.0)

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
        workers=workers if workers is not None else default_workers(cores),
        background=background,
        verbose=verbose,
        cores=cores,
        ram=ram,
    )


def preflight_checks(args: Args) -> None:
    """Filesystem checks that must pass before any geometry work begins."""
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
