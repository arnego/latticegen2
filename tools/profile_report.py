"""Join a `profile_run.py` CSV against a run's `.log` to get per-stage resource use.

The `.log` records when each stage *ended*; the CSV records what the machine
was doing every couple of seconds. Together they answer the question neither
can alone: which stages are single-threaded, which are memory-hungry, and which
are waiting on disk rather than working.

Usage::

    python tools/profile_report.py rehearsal/profile.csv rehearsal/run.log

Prints a Markdown table (paste-ready for docs/) plus the derived findings:
core-equivalents per stage against the machine's physical core count, so a
stage pinned at 1.0 stands out as a parallelization candidate.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import os
import re
import sys

STAGE_RE = re.compile(r"^\[(\d{2}:\d{2}:\d{2})\] stage (\w+): (.+)$")
START_RE = re.compile(r"latticegen2 run started (\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})")


def _walk_clock(stamps: list[str], origin: _dt.datetime) -> list[_dt.datetime]:
    """Turn HH:MM:SS strings into datetimes, rolling over at midnight.

    Both the log and the CSV carry wall-clock times only. A multi-hour run
    started in the evening crosses midnight -- the 2026-08-12 rehearsal ran
    19:21 to 00:21 -- so a naive parse puts most of the run *before* its own
    start and every stage window comes out negative.
    """
    out: list[_dt.datetime] = []
    day = origin.date()
    previous = origin
    for stamp in stamps:
        clock = _dt.datetime.strptime(stamp, "%H:%M:%S").time()
        moment = _dt.datetime.combine(day, clock)
        if moment < previous - _dt.timedelta(seconds=1):
            day += _dt.timedelta(days=1)
            moment = _dt.datetime.combine(day, clock)
        out.append(moment)
        previous = moment
    return out


def parse_log(path: str) -> tuple[_dt.datetime, list[tuple[str, _dt.datetime]]]:
    """Return the run's start time and each stage's end time, in order."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read().splitlines()

    match = next((START_RE.search(line) for line in text if START_RE.search(line)), None)
    if match is None:
        raise SystemExit(f"{path}: no run-start line found; is this a latticegen2 log?")
    origin = _dt.datetime.strptime(match.group(1), "%Y-%m-%d %H:%M:%S")

    names, stamps = [], []
    for line in text:
        hit = STAGE_RE.match(line)
        if hit:
            names.append(hit.group(2))
            stamps.append(hit.group(1))
    if not names:
        raise SystemExit(f"{path}: no stage lines found.")
    return origin, list(zip(names, _walk_clock(stamps, origin)))


def parse_csv(path: str, origin: _dt.datetime) -> list[dict]:
    with open(path, encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise SystemExit(f"{path}: no samples.")
    moments = _walk_clock([r["stamp"] for r in rows], origin)
    for row, moment in zip(rows, moments):
        row["_at"] = moment
    return rows


def _stats(rows: list[dict], key: str) -> tuple[float, float]:
    if not rows:
        return 0.0, 0.0
    values = [float(r[key]) for r in rows]
    return sum(values) / len(values), max(values)


def report(csv_path: str, log_path: str, cores: int) -> int:
    origin, stages = parse_log(log_path)
    rows = parse_csv(csv_path, origin)

    print(f"run started : {origin:%Y-%m-%d %H:%M:%S}")
    print(f"samples     : {len(rows)}  ({rows[0]['stamp']} -> {rows[-1]['stamp']})")
    print(f"stages      : {len(stages)}")
    print(f"cores       : {cores} physical (100% CPU == 1 core)")
    print()

    header = ("| Stage | Duration | CPU mean | CPU peak | Cores used | RSS mean | "
              "RSS peak | Procs | Read | Written |")
    print(header)
    print("|" + "---|" * 10)

    findings: list[tuple[float, str]] = []
    window_start = origin
    for name, end in stages:
        window = [r for r in rows if window_start < r["_at"] <= end]
        seconds = (end - window_start).total_seconds()
        cpu_mean, cpu_peak = _stats(window, "cpu_total_pct")
        rss_mean, rss_peak = _stats(window, "rss_total_mb")
        _, procs = _stats(window, "procs")
        if window:
            read = float(window[-1]["read_mb"]) - float(window[0]["read_mb"])
            write = float(window[-1]["write_mb"]) - float(window[0]["write_mb"])
        else:
            read = write = 0.0

        cores_used = cpu_mean / 100.0
        print(f"| {name} | {seconds:,.1f}s | {cpu_mean:,.0f}% | {cpu_peak:,.0f}% | "
              f"{cores_used:.2f} | {rss_mean:,.0f} MB | {rss_peak:,.0f} MB | "
              f"{procs:.0f} | {read:,.0f} MB | {write:,.0f} MB |")

        # A long stage that never used more than about one core is where
        # parallelism would actually pay. Weight by the time it could recover.
        # `cores` directly, not `min(cores, 8)`: this used to mirror the worker
        # count's cap at 8, and that cap is gone — one worker per core, honoured
        # exactly (docs/specification.md §3).
        if seconds >= 60 and cores_used < 1.35:
            recoverable = seconds * (1 - 1 / cores)
            findings.append((recoverable, (
                f"{name}: {seconds:,.0f}s at {cores_used:.2f} cores "
                f"(single-threaded). Perfect {cores}-way parallelism would "
                f"recover up to ~{recoverable:,.0f}s.")))
        window_start = end

    total = (stages[-1][1] - origin).total_seconds()
    _, rss_peak_all = _stats(rows, "rss_total_mb")
    avail_min = min(float(r["sys_avail_mb"]) for r in rows)
    print()
    print(f"total to last stage : {total:,.1f}s ({total / 60:,.1f} min)")
    print(f"peak tree RSS       : {rss_peak_all:,.0f} MB")
    print(f"min system available: {avail_min:,.0f} MB")
    print(f"total written       : {float(rows[-1]['write_mb']):,.0f} MB")
    print(f"total read          : {float(rows[-1]['read_mb']):,.0f} MB")

    print()
    print("Parallelization candidates, by recoverable wall time")
    print("-" * 72)
    if not findings:
        print("  none: every stage over a minute already used more than ~1.35 cores.")
    for _, text in sorted(findings, reverse=True):
        print(f"  - {text}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("csv", help="profile CSV from tools/profile_run.py")
    parser.add_argument("log", help="the run's .log file")
    parser.add_argument("--cores", type=int, default=os.cpu_count() or 1,
                        help="physical cores, for the core-equivalent column")
    args = parser.parse_args(argv)
    for path in (args.csv, args.log):
        if not os.path.exists(path):
            sys.stderr.write(f"not found: {path}\n")
            return 2
    return report(args.csv, args.log, args.cores)


if __name__ == "__main__":
    raise SystemExit(main())
