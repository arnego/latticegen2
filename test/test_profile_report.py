"""tools/profile_report.py: joining wall-clock stamps to the run's start."""
import datetime as dt
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))

import profile_report as pr  # noqa: E402

ORIGIN = dt.datetime(2026, 11, 6, 8, 12, 37)


def test_a_sample_taken_just_before_the_start_line_stays_on_the_same_day():
    got = pr._walk_clock(["08:12:35", "08:12:38", "09:00:00"], ORIGIN)
    assert [m.date() for m in got] == [ORIGIN.date()] * 3


def test_a_run_crossing_midnight_rolls_over():
    origin = dt.datetime(2026, 8, 12, 19, 21, 0)
    got = pr._walk_clock(["19:21:02", "23:59:59", "00:00:01", "00:21:00"], origin)
    day = origin.date()
    assert [m.date() for m in got] == [day, day, day + dt.timedelta(days=1),
                                       day + dt.timedelta(days=1)]


def test_a_first_stamp_after_midnight_belongs_to_the_next_day():
    origin = dt.datetime(2026, 8, 12, 23, 59, 59)
    got = pr._walk_clock(["00:00:01"], origin)
    assert got[0] == dt.datetime(2026, 8, 13, 0, 0, 1)
