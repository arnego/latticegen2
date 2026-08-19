"""The front-end's testable half (specification.md §3.1).

Everything that decides *what* the window shows lives outside
:mod:`latticegen2.gui.app` precisely so it can be tested without a display: the
stage weights, the reduction from an event stream to a state, and the command
line handed to the child. Only one test here needs Tk at all, and it withdraws
the window rather than showing it.

The stream these tests feed the model is deliberately hostile. Its producer is a
separate process that can be killed at any moment, so a missing ``stage_end``, a
half-written final line and kernel chatter on stdout are ordinary inputs, not
edge cases.
"""

import json

import pytest

from latticegen2 import progress
from latticegen2.cli import resolve_output_paths
from latticegen2.gui import model, runner, weights
from latticegen2.runlog import STAGES


def event(ev, **fields):
    """One serialized event, as it would arrive on the child's stdout."""
    payload = {"ev": ev, "t": fields.pop("t", 1.0)}
    payload.update({k: fields.get(k) for k in progress.FIELDS[ev]})
    payload.update(fields)
    return json.dumps(payload)


@pytest.fixture(scope="module")
def tk_root():
    """**One** withdrawn Tk root, shared by every test in this file.

    Module scope is the point, not an optimization. Creating and destroying a
    root per test made these tests skip *intermittently* — roughly one run in
    three — with a `TclError` reading "tk wasn't installed properly", which the
    old per-test guard turned into "no display available". There was a display
    and Tk was fine; repeatedly tearing a root down and standing another up in
    one interpreter simply is not reliable, and adding an `update()` after
    `destroy()` did not make it so. A test that quietly skips is a test that is
    not running, which is worse than one that fails.

    Widgets from earlier tests are left in place rather than cleared: the root
    is withdrawn so nothing is drawn, each test only inspects the objects it
    built itself, and destroying a widget with a pending `after` — the
    indeterminate bar always has one — just produces Tk callback noise on
    stderr.
    """
    tk = pytest.importorskip("tkinter")
    try:
        root = tk.Tk()
    except tk.TclError as exc:  # a genuinely headless machine
        pytest.skip(f"tkinter cannot open a window here: {exc}")
    root.withdraw()
    yield root
    root.destroy()


# --- the weight table --------------------------------------------------------


def test_every_stage_has_a_share_and_they_sum_to_one_thousand():
    """The guard against a stage being added and the bar silently ignoring it.

    :class:`latticegen2.runlog.Timer` already refuses an unknown stage name, so
    between the two a new stage cannot reach production without being weighted.
    """
    assert set(weights.STAGE_PERMILLE) == set(STAGES)
    assert sum(weights.STAGE_PERMILLE.values()) == 1000


def test_the_bands_tile_the_bar_in_pipeline_order():
    starts = [weights.BAND_START[name] for name in STAGES]
    assert starts == sorted(starts)
    assert starts[0] == 0
    last = STAGES[-1]
    assert weights.BAND_START[last] + weights.STAGE_PERMILLE[last] == 1000


def test_a_stage_with_no_countable_work_holds_at_its_band_start():
    """The user's decision, and the honest one: no invented fraction.

    ``export`` is 10.7 % of a rehearsal run behind one opaque kernel call. The
    top bar holds while the sub-bar sweeps, rather than creeping on a timer and
    reporting a number that measures nothing.
    """
    start = weights.BAND_START["export"]
    assert weights.overall_permille("export", None) == start
    assert weights.overall_permille("export", 0.5) > start
    assert weights.overall_permille("export", 1.0) == 1000


def test_a_fraction_outside_zero_to_one_cannot_escape_its_band():
    for stage in STAGES:
        band_end = weights.BAND_START[stage] + weights.STAGE_PERMILLE[stage]
        assert weights.overall_permille(stage, -5.0) == weights.BAND_START[stage]
        assert weights.overall_permille(stage, 99.0) == band_end


# --- the reduction -----------------------------------------------------------


def test_a_whole_run_reduces_to_a_finished_state():
    lines = [event(progress.HELLO, v=1, pid=1, stages=list(STAGES),
                   output="out.step", log="out.log", workers=6)]
    for i, name in enumerate(STAGES):
        lines.append(event(progress.STAGE_BEGIN, name=name, i=i, n=12, t=i))
        lines.append(event(progress.STAGE_END, name=name, i=i, n=12,
                           elapsed=1.0, max_rss=1000 * (i + 1), t=i + 0.5))
    lines.append(event(progress.SUMMARY, duration=42.0, peak_rss=99999,
                       params={}, stats={"solids_written": "14"},
                       stages=[], t=42.0))
    lines.append(event(progress.RESULT, status=progress.OK, exit=0, reason=None,
                       output="out.step", log="out.log", tmpdir=None, t=42.1))

    state = model.reduce_stream(lines)
    assert state.finished and state.status == progress.OK
    assert state.permille == 1000
    assert state.peak_rss == 99999
    assert "14 solid(s)" in state.result_text()


def test_the_bar_never_goes_backwards():
    """Not decoration. ``stage_end`` fills a band completely while a late tick
    from the same stage would compute a smaller value, and a bar that retreats
    reads as something having gone wrong — which here would be false.
    """
    lines = [
        event(progress.STAGE_BEGIN, name="boundary", i=4, n=12),
        event(progress.PROGRESS, stage="boundary", label="trim", done=9, total=10),
        event(progress.STAGE_END, name="boundary", i=4, n=12, elapsed=1.0, max_rss=0),
        event(progress.PROGRESS, stage="boundary", label="trim", done=1, total=10),
    ]
    seen = []
    state = model.RunState()
    for line in lines:
        state.apply(progress.parse_line(line))
        seen.append(state.permille)
    assert seen == sorted(seen)


def test_a_missing_stage_end_does_not_strand_the_bar():
    """A killed child leaves exactly this: a stage begun and never closed."""
    state = model.reduce_stream([
        event(progress.STAGE_BEGIN, name="classify", i=3, n=12),
        event(progress.STAGE_BEGIN, name="simplify", i=9, n=12),
    ])
    assert state.stage == "simplify"
    assert state.permille == weights.BAND_START["simplify"]


def test_junk_on_stdout_is_kept_as_text_and_never_raises():
    """OCCT writes to file descriptor 1 below ``sys.stdout``, and a process
    killed mid-write leaves a half-line. Neither may stop the reader."""
    state = model.reduce_stream([
        "Intel oneMKL FATAL ERROR",
        '{"ev":"stage_begin","name":"cla',
        event(progress.STAGE_BEGIN, name="connect", i=5, n=12),
        "",
    ])
    assert state.stage == "connect"
    assert [text for text, _console in state.lines] == [
        "Intel oneMKL FATAL ERROR", '{"ev":"stage_begin","name":"cla',
    ]


def test_the_log_backlog_is_bounded():
    state = model.reduce_stream(
        [event(progress.LOG, msg=f"line {i}", console=False)
         for i in range(model.LOG_BACKLOG * 2)]
    )
    assert len(state.lines) == model.LOG_BACKLOG
    assert state.lines[-1][0] == f"line {model.LOG_BACKLOG * 2 - 1}"


def test_a_console_flag_travels_with_every_line():
    """What makes verbosity a filter the reader can change during a run, rather
    than a flag that had to be decided before it started."""
    state = model.reduce_stream([
        event(progress.LOG, msg="quiet", console=False),
        event(progress.LOG, msg="loud", console=True),
    ])
    assert list(state.lines) == [("quiet", False), ("loud", True)]


def test_an_indeterminate_tick_switches_the_sub_bar_and_keeps_its_label():
    state = model.reduce_stream([
        event(progress.STAGE_BEGIN, name="export", i=11, n=12),
        event(progress.PROGRESS, stage="export", label="writing STEP",
              done=0, total=None),
    ])
    assert state.indeterminate and state.sub_fraction is None
    assert state.sub_text() == "writing STEP"


def test_a_countable_tick_reads_as_done_over_total():
    state = model.reduce_stream([
        event(progress.STAGE_BEGIN, name="boundary", i=4, n=12),
        event(progress.PROGRESS, stage="boundary", label="boundary trim",
              done=9776, total=19552),
    ])
    assert state.sub_fraction == pytest.approx(0.5)
    assert state.sub_text() == "boundary trim: 9,776 / 19,552"


def test_a_cancelled_run_reports_where_its_temp_folder_is():
    """specification.md §3 keeps it deliberately, so the window must say where."""
    state = model.reduce_stream([
        event(progress.RESULT, status=progress.CANCELLED, exit=130,
              reason="interrupted by user. Intermediate files kept.",
              output="out.step", log="out.log", tmpdir=r"C:\p\temp\20260819-004500",
              t=12.0),
    ])
    assert state.status == progress.CANCELLED
    assert state.tmpdir.endswith("20260819-004500")
    assert "Cancelled after" in state.result_text()


def test_a_failure_reports_the_single_reason_line():
    state = model.reduce_stream([
        event(progress.RESULT, status=progress.FAILED, exit=4,
              reason="1 of 14 output solids failed OCCT's check",
              output="out.step", log="out.log", tmpdir="/tmp/x", t=9.0),
    ])
    assert state.result_text() == "Failed: 1 of 14 output solids failed OCCT's check"


# --- the command line handed to the child ------------------------------------


def test_the_child_gets_exactly_the_expected_flags_and_never_dash_v():
    argv = runner.build_argv(input_path="in.step", output_path="out.step",
                             cc="10", t="1.5", cores=6)
    assert argv[1] == "-u"
    assert argv[2].replace("\\", "/").endswith("src/main.py")
    assert argv[3:] == ["-i", "in.step", "-cc", "10", "-t", "1.5",
                        "-o", "out.step", "--cores", "6", "--progress-stream"]
    # Verbosity is a client-side filter: every line arrives as an event carrying
    # whether a CLI run would have printed it.
    assert "-v" not in argv and "--verbose" not in argv


def test_the_child_argv_is_one_the_parser_accepts(tmp_path):
    """The window and the run must not disagree about what a valid command is."""
    from latticegen2.cli import parse_args

    step = tmp_path / "part.step"
    step.write_text("dummy")
    argv = runner.build_argv(input_path=str(step), output_path=str(tmp_path / "o.step"),
                             cc="12", t="2.5", cores=4)
    args = parse_args(argv[3:])
    assert args.cc == 12.0 and args.t == 2.5 and args.cores == 4
    assert args.progress_stream is True


@pytest.mark.parametrize("cc,t,expected", [
    (20.0, 4.0, "part-lattice-cc20t4.step"),
    (10.0, 1.5, "part-lattice-cc10t1.5.step"),
    (12.0, 2.5, "part-lattice-cc12t2.5.step"),
])
def test_the_derived_filename_is_the_cli_default(tmp_path, cc, t, expected):
    """The output name the window shows has to be the one the run will write.

    It comes from :func:`latticegen2.cli.resolve_output_paths` for that reason —
    ``format_param`` renders ``20.0`` as ``20`` and ``1.5`` as ``1.5``, and a
    second implementation of that rule in the front-end would drift from it.
    """
    step, _log = resolve_output_paths(str(tmp_path / "part.step"), None, cc, t)
    import os

    assert os.path.basename(step) == expected


# --- the window itself, minimally ---------------------------------------------


def test_the_window_builds_and_validates_with_the_parser(tmp_path, tk_root):
    root = tk_root
    from latticegen2.gui.app import App

    step = tmp_path / "part.step"
    step.write_text("dummy")
    app = App(root)
    app.input_var.set(str(step))
    app.cc_var.set("20")
    app.t_var.set("4")
    root.update()
    assert app.derived_var.get().endswith("part-lattice-cc20t4.step")
    assert "disabled" not in app.start.state()

    # The one cross-parameter rule, enforced by `parse_args` rather than by
    # a copy of it living here.
    app.t_var.set("18")
    root.update()
    assert "must be smaller than the cube edge" in app.message_var.get()
    assert "disabled" in app.start.state()


# --- the cancel channel (docs/algorithm.md §10) ------------------------------


def test_the_cancel_sentinel_sits_beside_the_log():
    """Derived rather than passed, so both sides compute it from one fact.

    A flag would be a third place the path could be spelled, and the two
    processes would have to agree about it on top of already agreeing about the
    log path.
    """
    assert progress.cancel_sentinel_path("/p/part-lattice-cc10t1.5.log") == \
        "/p/part-lattice-cc10t1.5.cancel"


def test_a_stale_sentinel_is_removed_rather_than_obeyed(tmp_path, monkeypatch):
    """Left over from a window that closed before it could tidy up.

    Honouring it would kill a fresh run instantly, for a reason nothing in the
    log or the window could explain.
    """
    from latticegen2.__main__ import _watch_for_cancel

    log = tmp_path / "part.log"
    sentinel = tmp_path / "part.cancel"
    sentinel.write_text("from a previous session")

    interrupted = []
    monkeypatch.setattr("_thread.interrupt_main", lambda: interrupted.append(True))
    _watch_for_cancel(str(log))

    import time

    time.sleep(0.05)
    assert not sentinel.exists()
    assert not interrupted


def test_writing_the_sentinel_interrupts_the_run(tmp_path, monkeypatch):
    from latticegen2.__main__ import CANCEL_POLL_S, _watch_for_cancel

    log = tmp_path / "part.log"
    interrupted = []
    monkeypatch.setattr("_thread.interrupt_main", lambda: interrupted.append(True))
    _watch_for_cancel(str(log))

    import time

    (tmp_path / "part.cancel").write_text("stop")
    deadline = time.time() + 5
    while not interrupted and time.time() < deadline:
        time.sleep(CANCEL_POLL_S / 2)
    assert interrupted, "Stop did not reach the run"
    # Consumed as it is acted on, so it cannot linger into the next run.
    assert not (tmp_path / "part.cancel").exists()


# --- what the code review found (all three anchored here) --------------------


def test_the_run_is_not_finished_until_both_streams_reach_eof(tmp_path, tk_root):
    """A successful run must never be reported as a failure.

    The child's exit code is available the instant it dies, while its pipes can
    still hold buffered lines and a reader thread can be between its read and
    its `put`. Finishing on "the queue is empty" alone therefore raced the
    terminal `result` event, and a run that wrote its STEP correctly showed a
    red banner reading "the run ended with exit code 0".
    """
    root = tk_root
    from latticegen2.gui.app import App

    app = App(root)

    class DeadChild:
        """Exited, but with its result still in flight behind one reader."""

        running = False
        cancel_path = str(tmp_path / "x.cancel")
        stderr_lines: list = []

        def __init__(self):
            import queue

            self.events = queue.Queue()

        def returncode(self):
            return 0

        def clear_cancel(self):
            pass

    app.run = DeadChild()
    app._closed_seen = 1          # only stdout has ended so far
    app._drain()
    assert app.run is not None, "finished while a reader was still live"

    app.run.events.put(("event", progress.parse_line(event(
        progress.RESULT, status=progress.OK, exit=0, reason=None,
        output="out.step", log="out.log", tmpdir=None))))
    app.run.events.put(("closed", None))
    app._drain()
    assert app.run is None
    assert app.state.status == progress.OK


def test_force_stop_is_scheduled_on_the_wall_clock(tmp_path, monkeypatch, tk_root):
    """The grace period has to expire even when the child says nothing.

    It was computed from `state.elapsed`, which only advances when an event
    arrives — so during the long GIL-holding kernel call that is the whole
    reason Stop is not instant, the deadline could never be reached and the
    window would sit on "Cancelling..." forever.
    """
    root = tk_root
    from latticegen2.gui.app import App
    from latticegen2.gui.runner import GRACE_SECONDS

    app = App(root)

    class SilentChild:
        running = True
        cancel_path = str(tmp_path / "x.cancel")
        stderr_lines: list = []
        killed = False

        def __init__(self):
            import queue

            self.events = queue.Queue()

        def cancel(self):
            pass

        def force_stop(self):
            self.killed = True

        def returncode(self):
            return None

    app.run = SilentChild()
    app.state.elapsed = 0.0        # nothing has been reported, and won't be
    app._on_stop()
    assert app._force_at is not None

    # Capture the real clock first: `latticegen2.gui.app.time` is the same
    # module object as `time`, so a lambda calling through it after the
    # patch would call itself.
    import time as _time

    real_monotonic = _time.monotonic
    monkeypatch.setattr(
        _time, "monotonic",
        lambda: real_monotonic() + GRACE_SECONDS + 1,
    )
    app._drain()
    assert app.run.killed, "the grace period never expired"


def test_a_missing_main_script_is_a_message_rather_than_a_silent_no_op(monkeypatch):
    """`pip install .` puts this package where there is no `src/` above it.

    The console script then reaches the window with no arguments, and `Popen`
    would raise `FileNotFoundError` into a Tk callback — where a windowed
    process has nowhere to print it, so pressing Start did nothing at all.
    """
    from latticegen2.gui import runner

    monkeypatch.setattr(runner.os.path, "isfile", lambda _p: False)
    with pytest.raises(runner.LaunchError, match="src/main.py"):
        runner.main_script()


def test_the_sub_bar_stops_and_empties_when_a_run_ends(tk_root):
    """Success, cancellation and failure alike.

    After the terminal event there is no sub-stage left, so ``sub_fraction`` is
    ``None`` — the same value that means "sweeping" during ``export``. Left to
    itself the bar therefore kept animating for as long as the window stayed
    open, which is both untrue and, being the only thing still moving, reads as
    a run that has not stopped.
    """
    root = tk_root
    from latticegen2.gui.app import App

    for status in (progress.OK, progress.CANCELLED, progress.FAILED):
        app = App(root)
        app.state.apply(progress.parse_line(event(
            progress.STAGE_BEGIN, name="export", i=11, n=12)))
        app.state.apply(progress.parse_line(event(
            progress.PROGRESS, stage="export", label="writing STEP",
            done=0, total=None)))
        app._refresh()
        assert app.sub_bar._sweep is not None, "the sweep never started"

        app.state.apply(progress.parse_line(event(
            progress.RESULT, status=status, exit=0, reason="x",
            output="o.step", log="o.log", tmpdir=None)))
        app._refresh()
        assert app.sub_bar._sweep is None, f"{status}: still sweeping"
        assert app.sub_bar._fraction == 0.0, f"{status}: bar not empty"
        assert app.sub_bar._text == "", f"{status}: label not cleared"


def test_going_indeterminate_twice_does_not_double_the_animation(tk_root):
    """A pending ``after`` cannot be cancelled by clearing the sweep state.

    It still fires; if the bar has gone indeterminate again by then it would see
    a live sweep, reschedule itself, and run alongside the new loop at twice the
    speed. The generation counter is what retires it.
    """
    root = tk_root
    from latticegen2.gui.app import Bar

    bar = Bar(root)
    bar.set(None, "first")
    stale = bar._sweep_id
    bar.set(0.5, "determinate")
    bar.set(None, "second")
    assert bar._sweep_id != stale

    before = bar._sweep
    bar._animate(stale)          # the callback left over from the first run
    assert bar._sweep == before, "a stale callback advanced the sweep"
