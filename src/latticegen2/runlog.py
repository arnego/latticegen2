"""Run logging: the always-on ``.log`` file, stage timings, and the summary.

A full log is written every run regardless of ``-v``, which only raises *console*
verbosity, and the end-of-run summary carries every field specification.md §3
requires: parameters, start time, duration, run characteristics, peak memory and
output path (docs/algorithm.md §10).

Optionally the same run also reports itself as machine-readable events, for
the graphical front-end (:mod:`latticegen2.progress`). This class is where
that happens because it is already the observer the pipeline injects -- ``rl``
is threaded through every stage -- so nothing new has to be passed around. The
emission is strictly additive: with no emitter attached every method behaves
exactly as it did, and the ``.log`` file receives the same bytes either way.
"""

from __future__ import annotations

import datetime as _dt
import os
import sys
import time
from dataclasses import dataclass, field

from . import progress

#: Every pipeline stage, in the order :func:`latticegen2.pipeline.run_pipeline`
#: runs them, which is also the order the summary reports them in.
#:
#: Declared here rather than left implicit in the twelve ``Timer(rl, name)``
#: call sites, because two things outside the pipeline now depend on it: the
#: ``hello`` event announces it, and :mod:`latticegen2.gui.weights` gives each
#: stage a share of the progress bar. :class:`Timer` checks its name against
#: this tuple, so a stage renamed in one place and not the other is an
#: immediate failure rather than a silently mis-weighted bar.
STAGES: tuple[str, ...] = (
    "template", "import", "tessellate", "classify", "boundary", "connect",
    "stitch", "instance", "assemble", "simplify", "validate", "export",
)

#: Shortest gap between two sub-stage events, in seconds.
#:
#: Not polish. ``boundary``'s sequential path calls its progress callback once
#: per junction (:func:`latticegen2.boundary.trim_boundary`) -- 19,552 times on
#: the docs/specification.md §10 rehearsal, which is more than a bar can use
#: and more than a pipe should carry. The terminal tick is always emitted
#: regardless, so a stage still ends on its true count.
SUBSTAGE_MIN_INTERVAL_S = 0.1


def peak_rss_bytes() -> int:
    """Peak resident set size of this process, or 0 when unavailable."""
    if os.name == "nt":
        try:
            import ctypes
            import ctypes.wintypes

            class _Counters(ctypes.Structure):
                _fields_ = [
                    ("cb", ctypes.wintypes.DWORD),
                    ("PageFaultCount", ctypes.wintypes.DWORD),
                    ("PeakWorkingSetSize", ctypes.c_size_t),
                    ("WorkingSetSize", ctypes.c_size_t),
                    ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                    ("PagefileUsage", ctypes.c_size_t),
                    ("PeakPagefileUsage", ctypes.c_size_t),
                ]

            kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
            psapi = ctypes.WinDLL("psapi", use_last_error=True)
            kernel32.GetCurrentProcess.restype = ctypes.c_void_p
            psapi.GetProcessMemoryInfo.argtypes = [
                ctypes.c_void_p,
                ctypes.POINTER(_Counters),
                ctypes.wintypes.DWORD,
            ]
            psapi.GetProcessMemoryInfo.restype = ctypes.wintypes.BOOL
            counters = _Counters()
            counters.cb = ctypes.sizeof(counters)
            ok = psapi.GetProcessMemoryInfo(
                kernel32.GetCurrentProcess(), ctypes.byref(counters), counters.cb
            )
            return int(counters.PeakWorkingSetSize) if ok else 0
        except Exception:
            return 0
    try:
        import resource

        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        # Linux reports kibibytes, macOS reports bytes.
        return int(peak) * (1 if sys.platform == "darwin" else 1024)
    except Exception:
        return 0


def format_bytes(n: int) -> str:
    if n <= 0:
        return "n/a"
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.2f} {unit}" if unit != "B" else f"{n} B"
        n /= 1024.0
    return f"{n:.2f} TB"


def format_duration(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.2f}s"
    minutes, secs = divmod(seconds, 60)
    if minutes < 60:
        return f"{int(minutes)}m {secs:04.1f}s"
    hours, minutes = divmod(int(minutes), 60)
    return f"{hours}h {minutes:02d}m {secs:04.1f}s"


@dataclass
class RunLog:
    """A run's log file plus the counters the end-of-run summary reports."""

    path: str
    verbose: bool = False
    start_time: _dt.datetime = field(default_factory=_dt.datetime.now)
    t0: float = field(default_factory=time.perf_counter)
    stages: list[tuple[str, float]] = field(default_factory=list)
    max_rss: int = 0
    _fh: object = None
    #: A :class:`latticegen2.progress.NdjsonEmitter`, or ``None`` for a plain
    #: CLI run. Declared last so no existing construction changes meaning.
    emit: object = None
    #: The run's temp folder while one exists, and ``None`` once it has been
    #: removed on success. Recorded rather than left to be scraped out of the
    #: log text, because "where are the intermediate files" is the first
    #: question after a failure and the front-end has to answer it in a window
    #: (specification.md §4.4).
    tmpdir: "str | None" = None
    _stage: "str | None" = field(default=None, repr=False)
    _last_substage: float = field(default=0.0, repr=False)

    def open(self) -> "RunLog":
        self._fh = open(self.path, "w", encoding="utf-8")
        return self

    def close(self) -> None:
        if self._fh is not None:
            self._fh.close()
            self._fh = None

    def __enter__(self):
        return self.open()

    def __exit__(self, *exc):
        self.close()
        return False

    # -- writing ----------------------------------------------------------

    def _emit(self, ev: str, **fields) -> None:
        if self.emit is not None:
            self.emit(ev, **fields)

    def line(self, msg: str, console: bool | None = None) -> None:
        """Write one line to the log, and to the console when appropriate.

        ``console=True`` forces it to the console regardless of ``-v``; the
        default shows it only in verbose mode.

        With an emitter attached the console *print* is replaced by a ``log``
        event carrying that same decision as a flag, rather than acted on here.
        Two consequences, both wanted. The front-end's stdout stays pure NDJSON,
        which matters because OCCT writes to file descriptor 1 below
        ``sys.stdout`` and a mixed stream could not be parsed reliably. And
        because every line is emitted rather than only the ones ``-v`` would
        have shown, verbosity becomes a filter the reader can change *during* a
        run instead of a flag fixed before it starts.

        The log-file write above is deliberately outside all of that: the same
        bytes either way (docs/algorithm.md §10).
        """
        stamp = _dt.datetime.now().strftime("%H:%M:%S")
        if self._fh is not None:
            self._fh.write(f"[{stamp}] {msg}\n")
            self._fh.flush()
        show = console if console is not None else self.verbose
        if self.emit is not None:
            self._emit(progress.LOG, msg=msg, console=show)
        elif show:
            print(msg, flush=True)

    def always(self, msg: str) -> None:
        self.line(msg, console=True)

    def stage_begin(self, name: str) -> None:
        """Announce a stage as it starts. Emits only; writes nothing to the log.

        The asymmetry with :meth:`stage` is the point. A reader has to be told
        that ``simplify`` is what has been running for the last nineteen
        minutes, but the ``.log``'s format is documented and is compared across
        sessions, so nothing may be added to it for the front-end's benefit.
        """
        self._stage = name
        self._last_substage = 0.0
        self._emit(progress.STAGE_BEGIN, name=name, i=STAGES.index(name), n=len(STAGES))

    def stage(self, name: str, elapsed: float) -> None:
        self.stages.append((name, elapsed))
        self.observe_rss()
        self.line(f"stage {name}: {format_duration(elapsed)}", console=self.verbose)
        self._emit(
            progress.STAGE_END, name=name, i=STAGES.index(name), n=len(STAGES),
            elapsed=elapsed, max_rss=self.max_rss,
        )
        self._stage = None

    def substage(self, label: str, done: int, total: "int | None") -> None:
        """Report progress within the running stage. Emits only, rate-limited.

        ``total=None`` means the stage has no countable unit of work -- one
        opaque kernel call, or one solid so much larger than the rest that a job
        count would mislead. Saying so is better than inventing a fraction.
        """
        if self.emit is None:
            return
        now = time.monotonic()
        if (total is not None and done < total
                and now - self._last_substage < SUBSTAGE_MIN_INTERVAL_S):
            return
        self._last_substage = now
        self._emit(
            progress.PROGRESS, stage=self._stage, label=label, done=done, total=total,
        )

    def observe_rss(self) -> None:
        self.max_rss = max(self.max_rss, peak_rss_bytes())

    def note_worker_rss(self, rss: int) -> None:
        """Fold a worker process's peak RSS into the run's high-water mark."""
        self.max_rss = max(self.max_rss, rss)

    # -- structured sections ----------------------------------------------

    def header(self, params: dict) -> None:
        self.line("=" * 72)
        self.line(f"latticegen2 run started {self.start_time:%Y-%m-%d %H:%M:%S}")
        self.line("=" * 72)
        for key, value in params.items():
            self.line(f"  {key}: {value}")

    def summary(self, params: dict, stats: dict) -> None:
        """The mandatory end-of-run report, printed to console and log."""
        total = time.perf_counter() - self.t0
        self.observe_rss()
        self.always("")
        self.always("=" * 72)
        self.always("RUN SUMMARY")
        self.always("=" * 72)
        self.always(f"  Started:        {self.start_time:%Y-%m-%d %H:%M:%S}")
        self.always(f"  Duration:       {format_duration(total)}")
        self.always("  Parameters:")
        for key, value in params.items():
            self.always(f"    {key}: {value}")
        self.always("  Run characteristics:")
        for key, value in stats.items():
            self.always(f"    {key}: {value}")
        self.always("  Stage timings:")
        for name, elapsed in self.stages:
            self.always(f"    {name}: {format_duration(elapsed)}")
        self.always(f"  Peak memory:    {format_bytes(self.max_rss)}")
        self.always("=" * 72)
        # Structured as well as formatted: the front-end shows the same figures
        # in a banner, and re-parsing the text above to get them would make this
        # method's layout a machine contract.
        self._emit(
            progress.SUMMARY, duration=total, peak_rss=self.max_rss,
            params={k: str(v) for k, v in params.items()},
            stats={k: str(v) for k, v in stats.items()},
            stages=[[name, elapsed] for name, elapsed in self.stages],
        )


class Timer:
    """Context manager that records one pipeline stage's wall time."""

    def __init__(self, rl: RunLog, name: str):
        if name not in STAGES:
            raise ValueError(
                f"unknown pipeline stage {name!r}. Add it to runlog.STAGES and "
                f"give it a share in latticegen2.gui.weights, or the progress "
                f"bar silently stops accounting for it."
            )
        self.rl = rl
        self.name = name

    def __enter__(self):
        self.t0 = time.perf_counter()
        self.rl.stage_begin(self.name)
        return self

    def __exit__(self, *exc):
        self.rl.stage(self.name, time.perf_counter() - self.t0)
        return False
