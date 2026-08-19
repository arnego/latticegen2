"""The guarantee that the progress stream changes nothing (docs/algorithm.md §10).

This file exists for one assertion: **a run with an emitter attached writes the
same ``.log`` as a run without one.** Everything the front-end needs is additive,
and the moment that stops being true the log — which is the only record a failed
production run leaves — starts depending on how the run was started.

The rest of the file pins the three places that property was easy to get wrong:
``stage_begin`` must emit without logging, ``line`` must report the console
decision rather than act on it, and ``substage`` must throttle without ever
dropping a stage's final count.
"""

import re

import pytest

from latticegen2 import progress
from latticegen2.runlog import STAGES, SUBSTAGE_MIN_INTERVAL_S, RunLog, Timer


class Recorder:
    """Stands in for :class:`~latticegen2.progress.NdjsonEmitter`, keeping events."""

    def __init__(self):
        self.events = []

    def __call__(self, ev, **fields):
        assert set(fields) == set(progress.FIELDS[ev]), ev
        self.events.append((ev, fields))

    def of(self, ev):
        return [fields for name, fields in self.events if name == ev]


def exercise(rl):
    """One run's worth of every kind of write RunLog offers."""
    rl.header({"input": "part.step", "cc": "10 mm"})
    with Timer(rl, "template"):
        rl.line("junction template: 30 faces")
    with Timer(rl, "boundary"):
        rl.line("  boundary trim: 1/2")
        rl.substage("boundary trim", 1, 2)
        rl.substage("boundary trim", 2, 2)
        rl.always("peak worker memory: 1.00 GB")
    rl.summary({"input": "part.step"}, {"solids_written": 3})


def log_body(path):
    """The log file with its two genuinely volatile quantities masked.

    Wall-clock stamps and elapsed times differ between any two runs, and the
    measured peak RSS differs between two runs *in the same process*. None of
    the three is something an emitter could influence, so masking them is what
    leaves the comparison testing the thing it claims to.
    """
    text = path.read_text(encoding="utf-8")
    text = re.sub(r"^\[\d\d:\d\d:\d\d\] ", "", text, flags=re.M)
    return re.sub(r"\d+(\.\d+)?s", "<dur>", text)


def test_the_log_file_is_identical_with_and_without_an_emitter(tmp_path, monkeypatch):
    monkeypatch.setattr("latticegen2.runlog.peak_rss_bytes", lambda: 4096)
    plain = tmp_path / "plain.log"
    streamed = tmp_path / "streamed.log"

    exercise(RunLog(str(plain)).open())
    exercise(RunLog(str(streamed), emit=Recorder()).open())

    assert log_body(plain) == log_body(streamed)
    assert log_body(plain).strip()  # and it is not empty, so the check has teeth


def test_stage_begin_emits_but_writes_nothing_to_the_log(tmp_path):
    path = tmp_path / "run.log"
    rec = Recorder()
    rl = RunLog(str(path), emit=rec).open()
    rl.stage_begin("simplify")
    rl.close()

    assert path.read_text(encoding="utf-8") == ""
    assert rec.of(progress.STAGE_BEGIN) == [
        {"name": "simplify", "i": STAGES.index("simplify"), "n": 12}
    ]


def test_each_line_becomes_one_event_carrying_the_console_decision(tmp_path, capsys):
    rec = Recorder()
    rl = RunLog(str(tmp_path / "run.log"), verbose=False, emit=rec).open()
    rl.line("quiet")                    # verbose-only: a CLI run would not show it
    rl.always("loud")                   # forced to the console
    rl.line("forced off", console=False)

    assert [(e["msg"], e["console"]) for e in rec.of(progress.LOG)] == [
        ("quiet", False), ("loud", True), ("forced off", False),
    ]
    # Nothing was printed: with an emitter the stream must stay pure NDJSON, so
    # the console decision travels as data rather than as a print.
    assert capsys.readouterr().out == ""


def test_without_an_emitter_the_console_still_works(tmp_path, capsys):
    rl = RunLog(str(tmp_path / "run.log"), verbose=False).open()
    rl.line("quiet")
    rl.always("loud")
    assert capsys.readouterr().out == "loud\n"


def test_substage_throttles_but_never_drops_the_final_count(tmp_path, monkeypatch):
    rec = Recorder()
    rl = RunLog(str(tmp_path / "run.log"), emit=rec).open()
    rl.stage_begin("boundary")

    now = [1000.0]
    monkeypatch.setattr("latticegen2.runlog.time.monotonic", lambda: now[0])
    for done in range(1, 100):
        rl.substage("boundary trim", done, 100)   # all within one interval
    now[0] += SUBSTAGE_MIN_INTERVAL_S * 2
    rl.substage("boundary trim", 100, 100)

    ticks = rec.of(progress.PROGRESS)
    assert len(ticks) < 99, "the sequential boundary path would flood the pipe"
    assert ticks[-1]["done"] == 100
    assert all(t["stage"] == "boundary" for t in ticks)


def test_a_terminal_tick_is_emitted_even_inside_the_throttle_window(tmp_path, monkeypatch):
    rec = Recorder()
    rl = RunLog(str(tmp_path / "run.log"), emit=rec).open()
    rl.stage_begin("classify")
    monkeypatch.setattr("latticegen2.runlog.time.monotonic", lambda: 1000.0)
    rl.substage("slices", 1, 24)
    rl.substage("slices", 2, 24)          # throttled away
    rl.substage("slices", 24, 24)         # done == total: never throttled
    assert [t["done"] for t in rec.of(progress.PROGRESS)] == [1, 24]


def test_an_indeterminate_substage_is_never_throttled(tmp_path, monkeypatch):
    # `total=None` is emitted by stages that have no countable work, so there is
    # no terminal tick to fall back on and a dropped event is a bar that never
    # updates its label.
    rec = Recorder()
    rl = RunLog(str(tmp_path / "run.log"), emit=rec).open()
    rl.stage_begin("export")
    monkeypatch.setattr("latticegen2.runlog.time.monotonic", lambda: 1000.0)
    rl.substage("writing STEP", 0, None)
    rl.substage("writing STEP", 0, None)
    assert len(rec.of(progress.PROGRESS)) == 2


def test_substage_costs_nothing_without_an_emitter(tmp_path):
    rl = RunLog(str(tmp_path / "run.log")).open()
    rl.stage_begin("boundary")
    rl.substage("boundary trim", 1, 2)     # must not raise, must not log
    rl.close()
    assert (tmp_path / "run.log").read_text(encoding="utf-8") == ""


def test_stage_end_reports_the_same_elapsed_time_the_log_does(tmp_path):
    rec = Recorder()
    rl = RunLog(str(tmp_path / "run.log"), emit=rec).open()
    with Timer(rl, "connect"):
        pass
    logged_name, logged_elapsed = rl.stages[-1]
    ended = rec.of(progress.STAGE_END)[-1]
    assert ended["name"] == logged_name == "connect"
    assert ended["elapsed"] == logged_elapsed


def test_summary_is_emitted_structured_as_well_as_formatted(tmp_path, capsys):
    rec = Recorder()
    rl = RunLog(str(tmp_path / "run.log"), emit=rec).open()
    with Timer(rl, "export"):
        pass
    rl.summary({"cc": "10 mm"}, {"solids_written": 14})
    summary = rec.of(progress.SUMMARY)[-1]
    assert summary["params"] == {"cc": "10 mm"}
    assert summary["stats"] == {"solids_written": "14"}
    assert [name for name, _ in summary["stages"]] == ["export"]
    assert summary["duration"] > 0


def test_an_unknown_stage_name_is_refused_at_the_timer(tmp_path):
    rl = RunLog(str(tmp_path / "run.log")).open()
    with pytest.raises(ValueError, match="unknown pipeline stage"):
        Timer(rl, "weld")
