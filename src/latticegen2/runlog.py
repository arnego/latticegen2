"""Run logging: the always-on ``.log`` file, stage timings, and the summary.

A full log is written every run regardless of ``-v``, which only raises *console*
verbosity, and the end-of-run summary carries every field specification.md §3
requires: parameters, start time, duration, run characteristics, peak memory and
output path (docs/algorithm.md §10).
"""

from __future__ import annotations

import datetime as _dt
import os
import sys
import time
from dataclasses import dataclass, field


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

    def line(self, msg: str, console: bool | None = None) -> None:
        """Write one line to the log, and to the console when appropriate.

        ``console=True`` forces it to the console regardless of ``-v``; the
        default shows it only in verbose mode.
        """
        stamp = _dt.datetime.now().strftime("%H:%M:%S")
        if self._fh is not None:
            self._fh.write(f"[{stamp}] {msg}\n")
            self._fh.flush()
        if console if console is not None else self.verbose:
            print(msg, flush=True)

    def always(self, msg: str) -> None:
        self.line(msg, console=True)

    def stage(self, name: str, elapsed: float) -> None:
        self.stages.append((name, elapsed))
        self.observe_rss()
        self.line(f"stage {name}: {format_duration(elapsed)}", console=self.verbose)

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


class Timer:
    """Context manager that records one pipeline stage's wall time."""

    def __init__(self, rl: RunLog, name: str):
        self.rl = rl
        self.name = name

    def __enter__(self):
        self.t0 = time.perf_counter()
        return self

    def __exit__(self, *exc):
        self.rl.stage(self.name, time.perf_counter() - self.t0)
        return False
