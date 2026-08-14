"""Sample a latticegen2 run's resource use from outside the process.

Wraps a command, spawns it as a child, and polls the whole process tree
(master plus every boundary/sew worker) on a fixed interval, writing one CSV
row per sample. The child's exit code is propagated, so this is a transparent
prefix: anything that works on its own works under the profiler.

The point is to be *non-invasive*. It touches no pipeline code, and the run it
measures must stay comparable to one measured without it -- so sampling is a
single psutil poll per interval and nothing else.

Timestamps are written as ``HH:MM:SS`` to match ``runlog.RunLog.line``, which
is what lets a sample series be binned into the stage windows the ``.log``
records. Stage lines there are written when a stage *ends*, so a stage's window
runs from the previous stage's line to its own.

Usage::

    python tools/profile_run.py --out rehearsal/profile.csv -- \\
        python src/main.py -i part.step -cc 5 -t 1 -o out.step -v

``psutil`` is required, and is a **development-only** dependency: ``tools/`` is
``export-ignore``d in ``.gitattributes`` so it never reaches a release bundle,
and the shipped tool measures its own RSS through ctypes/psapi
(``runlog.peak_rss_bytes``) with no third-party package involved.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import os
import subprocess
import sys
import time

try:
    import psutil
except ImportError:  # pragma: no cover - environment-dependent
    sys.stderr.write(
        "profile_run.py needs psutil, which is a development-only dependency.\n"
        "Install it into the interpreter you run the tool with, e.g.:\n"
        "    python -m pip install psutil\n"
    )
    raise SystemExit(2)


FIELDS = [
    "stamp",          # HH:MM:SS, joins against the run's .log
    "elapsed_s",      # seconds since the child was spawned
    "procs",          # live processes in the tree, master included
    "rss_total_mb",   # summed RSS across the tree
    "rss_master_mb",
    "rss_max_worker_mb",
    "cpu_total_pct",  # summed; 100 == one core fully busy
    "cpu_master_pct",
    "threads_total",
    "sys_avail_mb",   # system-wide available RAM
    "read_mb",        # tree-cumulative, monotonic across worker exits
    "write_mb",
]

_MB = 1024.0 * 1024.0


class TreeSampler:
    """Polls a process tree, keeping cumulative I/O monotonic across exits.

    Worker processes come and go per stage. Their ``io_counters`` disappear
    with them, so retired processes' last-known totals are accumulated
    separately -- otherwise the tree's byte counts would drop every time a
    pool shut down.
    """

    def __init__(self, pid: int):
        self.master = psutil.Process(pid)
        self._seen: dict[int, psutil.Process] = {}
        self._io_last: dict[int, tuple[int, int]] = {}
        self._io_retired = [0, 0]

    def _tree(self) -> list[psutil.Process]:
        """Return the live tree as *cached* Process objects.

        Reusing the cached object per pid is load-bearing, not tidiness.
        ``cpu_percent(None)`` reports use since the last call *on that same
        object*, and ``children()`` hands back freshly constructed ones every
        poll -- so measuring through those reports 0.0 forever, while priming
        and reading one in the same pass reports a microsecond-wide window as
        a huge percentage. Both were observed before this was fixed. A newly
        appeared process is primed here and reads 0.0 for exactly one sample,
        which is honest.
        """
        found = [self.master]
        try:
            found.extend(self.master.children(recursive=True))
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass

        live: dict[int, psutil.Process] = {}
        for proc in found:
            cached = self._seen.get(proc.pid)
            if cached is None:
                self._seen[proc.pid] = cached = proc
                try:
                    cached.cpu_percent(None)  # prime; discard the 0.0
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    continue
            live[proc.pid] = cached

        for pid in list(self._io_last):
            if pid not in live:
                read, write = self._io_last.pop(pid)
                self._io_retired[0] += read
                self._io_retired[1] += write
        for pid in list(self._seen):
            if pid not in live and pid != self.master.pid:
                self._seen.pop(pid, None)
        return list(live.values())

    def sample(self) -> dict:
        procs = self._tree()
        rss_total = rss_master = rss_max_worker = 0
        cpu_total = cpu_master = 0.0
        threads = 0
        read_live = write_live = 0

        for proc in procs:
            is_master = proc.pid == self.master.pid
            try:
                with proc.oneshot():
                    rss = proc.memory_info().rss
                    cpu = proc.cpu_percent(None)
                    nthreads = proc.num_threads()
                    try:
                        io = proc.io_counters()
                        self._io_last[proc.pid] = (io.read_bytes, io.write_bytes)
                    except (psutil.AccessDenied, AttributeError):
                        pass
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue

            rss_total += rss
            cpu_total += cpu
            threads += nthreads
            if is_master:
                rss_master, cpu_master = rss, cpu
            else:
                rss_max_worker = max(rss_max_worker, rss)

        for read, write in self._io_last.values():
            read_live += read
            write_live += write

        return {
            "procs": len(procs),
            "rss_total_mb": round(rss_total / _MB, 1),
            "rss_master_mb": round(rss_master / _MB, 1),
            "rss_max_worker_mb": round(rss_max_worker / _MB, 1),
            "cpu_total_pct": round(cpu_total, 1),
            "cpu_master_pct": round(cpu_master, 1),
            "threads_total": threads,
            "sys_avail_mb": round(psutil.virtual_memory().available / _MB, 1),
            "read_mb": round((read_live + self._io_retired[0]) / _MB, 2),
            "write_mb": round((write_live + self._io_retired[1]) / _MB, 2),
        }


def _summarize(rows: list[dict], out_path: str, command: list[str],
               started: _dt.datetime, elapsed: float, returncode: int) -> str:
    """Write a plain-text companion to the CSV and return its path."""
    path = os.path.splitext(out_path)[0] + "-summary.txt"

    def series(name: str) -> list[float]:
        return [float(r[name]) for r in rows]

    with open(path, "w", encoding="utf-8") as fh:
        fh.write("latticegen2 run profile\n")
        fh.write("=" * 72 + "\n")
        fh.write(f"command:    {' '.join(command)}\n")
        fh.write(f"started:    {started:%Y-%m-%d %H:%M:%S}\n")
        fh.write(f"duration:   {elapsed:.1f} s\n")
        fh.write(f"exit code:  {returncode}\n")
        fh.write(f"samples:    {len(rows)}\n")
        fh.write(f"cores:      {psutil.cpu_count(logical=False)} physical, "
                 f"{psutil.cpu_count(logical=True)} logical\n")
        fh.write(f"system RAM: {psutil.virtual_memory().total / _MB / 1024:.1f} GB\n")
        if not rows:
            fh.write("\nNo samples taken -- the run finished inside one interval.\n")
            return path

        fh.write("\nPeaks and means\n")
        fh.write("-" * 72 + "\n")
        for name in ("rss_total_mb", "rss_master_mb", "rss_max_worker_mb",
                     "cpu_total_pct", "cpu_master_pct", "threads_total", "procs"):
            values = series(name)
            fh.write(f"  {name:<20} peak {max(values):>12,.1f}   "
                     f"mean {sum(values) / len(values):>12,.1f}\n")
        fh.write(f"  {'sys_avail_mb':<20} min  {min(series('sys_avail_mb')):>12,.1f}\n")
        fh.write(f"  {'read_mb':<20} final{series('read_mb')[-1]:>12,.2f}\n")
        fh.write(f"  {'write_mb':<20} final{series('write_mb')[-1]:>12,.2f}\n")

        # Windows where more than one process was alive: the parallel stages.
        fh.write("\nMulti-process windows (workers alive)\n")
        fh.write("-" * 72 + "\n")
        start = None
        for row in rows:
            multi = int(row["procs"]) > 1
            if multi and start is None:
                start = row
            elif not multi and start is not None:
                fh.write(f"  {start['stamp']} -> {row['stamp']}  "
                         f"({float(row['elapsed_s']) - float(start['elapsed_s']):.0f} s, "
                         f"up to {max(int(r['procs']) for r in rows):d} procs)\n")
                start = None
        if start is not None:
            fh.write(f"  {start['stamp']} -> end\n")
    return path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Sample a command's process tree and write a CSV profile.")
    parser.add_argument("--out", required=True,
                        help="CSV path; a -summary.txt companion is written alongside")
    parser.add_argument("--interval", type=float, default=2.0,
                        help="seconds between samples (default 2.0)")
    parser.add_argument("command", nargs=argparse.REMAINDER,
                        help="-- followed by the command to run")
    args = parser.parse_args(argv)

    command = args.command
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        parser.error("no command given; put it after a bare --")

    out_dir = os.path.dirname(os.path.abspath(args.out))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    started = _dt.datetime.now()
    t0 = time.perf_counter()
    child = subprocess.Popen(command)
    print(f"profiling pid {child.pid}: {' '.join(command)}", flush=True)

    rows: list[dict] = []
    sampler = TreeSampler(child.pid)
    with open(args.out, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS)
        writer.writeheader()
        try:
            while True:
                try:
                    returncode = child.wait(timeout=args.interval)
                    break
                except subprocess.TimeoutExpired:
                    pass
                row = sampler.sample()
                row["stamp"] = _dt.datetime.now().strftime("%H:%M:%S")
                row["elapsed_s"] = round(time.perf_counter() - t0, 1)
                rows.append(row)
                writer.writerow(row)
                fh.flush()
        except KeyboardInterrupt:
            # Let the child run its own graceful shutdown (spec §3, exit 130)
            # rather than killing it here; the partial profile is kept.
            print("profile_run: interrupted, waiting for the child to stop",
                  flush=True)
            returncode = child.wait()

    elapsed = time.perf_counter() - t0
    summary = _summarize(rows, args.out, command, started, elapsed, returncode)
    print(f"profile: {len(rows)} samples -> {args.out}", flush=True)
    print(f"summary: {summary}", flush=True)
    return returncode


if __name__ == "__main__":
    raise SystemExit(main())
