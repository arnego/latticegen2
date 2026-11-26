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
from latticegen2.cli import parse_args, resolve_output_paths
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


def test_the_supported_layout_has_bands_for_its_own_stages():
    """Two precomputed layouts, each summing to 1000 over its own stages."""
    from latticegen2.runlog import SUPPORT_STAGES, stages_for

    held = stages_for(True)
    assert set(weights.SUPPORT_PERMILLE) == set(held)
    assert sum(weights.SUPPORT_PERMILLE.values()) == 1000
    assert all(v >= 1 for v in weights.SUPPORT_PERMILLE.values())
    assert all(v >= 1 for v in weights.STAGE_PERMILLE.values())
    # Each support stage has a band to begin and end, and the measured
    # rehearsal puts the pocket search ahead of every other stage.
    assert all(weights.SUPPORT_PERMILLE[s] >= 1 for s in SUPPORT_STAGES)
    assert max(weights.SUPPORT_PERMILLE, key=weights.SUPPORT_PERMILLE.get) == "pockets"
    # The layout is selected by the announced stages, and by nothing else.
    assert weights.layout(held) == (weights.SUPPORT_PERMILLE, weights.SUPPORT_BAND_START)
    assert weights.layout(STAGES) == weights.layout() == (
        weights.STAGE_PERMILLE, weights.BAND_START)
    for stage in held:
        end = weights.SUPPORT_BAND_START[stage] + weights.SUPPORT_PERMILLE[stage]
        assert weights.overall_permille(stage, None, held) == weights.SUPPORT_BAND_START[stage]
        assert weights.overall_permille(stage, 1.0, held) == end
    assert weights.overall_permille("support", 0.0, STAGES) == 0     # not a plain stage


def test_shares_are_apportioned_exactly_and_never_to_nothing():
    stages = ("a", "b", "c", "d")
    share = weights._permille({"a": 0.0, "b": 1.0, "c": 1.0, "d": 1.0}, stages)
    assert sum(share.values()) == 1000 and share["a"] == 1
    assert sorted(share[s] for s in "bcd") == [333, 333, 333]
    # Nothing measured at all: an even split, still exact.
    assert sum(weights._permille({}, stages).values()) == 1000
    # A list this window has never seen still tiles the bar from 0 to 1000.
    odd = ("template", "something-new", "export")
    share, start = weights.layout(odd)
    assert sum(share.values()) == 1000 and share["something-new"] == 1
    assert start["template"] == 0 and start["export"] + share["export"] == 1000


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
    assert [text for text, _source in state.lines] == [
        "Intel oneMKL FATAL ERROR", '{"ev":"stage_begin","name":"cla',
    ]


def test_the_log_backlog_is_bounded():
    state = model.reduce_stream(
        [event(progress.LOG, msg=f"line {i}", console=False)
         for i in range(model.LOG_BACKLOG * 2)]
    )
    assert len(state.lines) == model.LOG_BACKLOG
    assert state.lines[-1][0] == f"line {model.LOG_BACKLOG * 2 - 1}"


def test_a_line_remembers_where_it_came_from():
    """What makes verbosity a filter the reader can change during a run, rather
    than a flag that had to be decided before it started.

    The run's own logging is one source and everything the child wrote outside
    the event stream is the other. The event's ``console`` flag is deliberately
    *not* the axis: a clean run marks its whole 66-line summary console-true,
    and every figure in it is either already drawn in a widget or a statistic
    that is neither a warning nor an error.
    """
    state = model.reduce_stream([
        event(progress.LOG, msg="stage boundary: 8.1s", console=False),
        event(progress.LOG, msg="  Peak memory:    282.73 MB", console=True),
        "Intel oneMKL FATAL ERROR",
    ])
    assert list(state.lines) == [
        ("stage boundary: 8.1s", model.LOG),
        ("  Peak memory:    282.73 MB", model.LOG),
        ("Intel oneMKL FATAL ERROR", model.OUT),
    ]


def test_only_what_the_run_wrote_outside_the_event_stream_shows_unticked():
    """The user's decision, and the reason for it.

    A clean run's every log line is either something the window already draws —
    the parameters sit in the input boxes, the stage in the top bar, duration
    and peak memory in the resource line — or a statistic that is neither a
    warning nor an error. What is left is what the child said unprompted, which
    is exactly the thing worth surfacing unasked.
    """
    state = model.reduce_stream([
        event(progress.LOG, msg="classification: 2925 nodes -> 33 interior",
              console=False),
        event(progress.LOG, msg="  Duration:       5.28s", console=True),
        "FAILED: 1 of 14 output solids failed OCCT's check",
    ])
    assert state.visible_lines(verbose=False) == [
        "FAILED: 1 of 14 output solids failed OCCT's check"]
    assert state.visible_lines(verbose=True) == [
        "classification: 2925 nodes -> 33 interior",
        "  Duration:       5.28s",
        "FAILED: 1 of 14 output solids failed OCCT's check",
    ]


def test_a_clean_run_leaves_the_unticked_view_completely_empty():
    """Nothing to read is the point: the result banner already carries the
    outcome, and the pane is hidden rather than shown empty."""
    state = model.reduce_stream(
        [event(progress.LOG, msg=f"  stat_{i}: {i}", console=True)
         for i in range(66)]
    )
    assert state.visible_lines(verbose=False) == []
    assert len(state.visible_lines(verbose=True)) == 66


WARNING = "WARNING: a title: what went wrong. What it may cause.\n  [1.000, 2.000, 3.000] mm"


def test_a_warning_shows_whichever_way_the_box_is_set():
    """It carries no flag of its own: the run's other console-true lines must
    stay hidden unticked, and a line merely mentioning the word is not one."""
    state = model.reduce_stream([
        event(progress.LOG, msg="stage boundary: 8.1s", console=False),
        event(progress.LOG, msg=WARNING, console=True),
        event(progress.LOG, msg="note: a repair, WARNING: not a warning", console=True),
    ])
    assert state.visible_lines(verbose=False) == [WARNING]
    assert len(state.visible_lines(verbose=True)) == 3


def test_a_warning_outlives_the_backlog():
    """The first thing a run says can be the most important thing an hour later."""
    state = model.reduce_stream(
        [event(progress.LOG, msg=WARNING, console=True)]
        + [event(progress.LOG, msg=f"line {i}", console=False)
           for i in range(model.LOG_BACKLOG * 2)]
    )
    assert state.visible_lines(verbose=False) == [WARNING]
    assert len(state.lines) == model.LOG_BACKLOG
    assert state.lines[0][0] == WARNING
    assert state.lines[-1][0] == f"line {model.LOG_BACKLOG * 2 - 1}"


def test_warning_text_says_what_what_it_may_cause_and_where():
    from latticegen2.runlog import WARNING_PLACES, box_text, point_text, warning_text

    text = warning_text("a title", "it broke", "It may leak",
                         [point_text((1, 2, 3)), box_text((0, 0, 0), (1, 2, 3))])
    assert text.split("\n") == [
        "WARNING: a title: it broke. It may leak. Inspect, in the part's coordinates:",
        "  [1.000, 2.000, 3.000] mm",
        "  x 0.000 to 1.000, y 0.000 to 2.000, z 0.000 to 3.000 mm",
    ]
    assert warning_text("a title", "it broke", "It may leak") == (
        "WARNING: a title: it broke. It may leak.")
    many = warning_text("t", "w", "r", [f"place {i}" for i in range(WARNING_PLACES + 3)])
    assert many.endswith("  ... and 3 more")


def test_lines_seen_keeps_counting_after_the_backlog_is_full():
    """What the window redraws its log pane on.

    ``len(lines)`` stops changing the moment the deque saturates, so keying the
    redraw on it would freeze the pane after 400 lines — at exactly the point in
    a long run where there is most to read.
    """
    state = model.reduce_stream(
        [event(progress.LOG, msg=f"line {i}", console=False)
         for i in range(model.LOG_BACKLOG + 25)]
    )
    assert len(state.lines) == model.LOG_BACKLOG
    assert state.lines_seen == model.LOG_BACKLOG + 25


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
    assert state.result_text() == "Cancelled after 0h 00m 12s"


@pytest.mark.parametrize("seconds, text", [
    (0.0, "0h 00m 00s"),
    (12.99, "0h 00m 12s"),
    (59.999, "0h 00m 59s"),
    (60.0, "0h 01m 00s"),
    (3725.4, "1h 02m 05s"),
    (36000.0, "10h 00m 00s"),
])
def test_the_window_counts_in_whole_seconds_as_hours_minutes_seconds(seconds, text):
    """Truncated, never rounded up: a counter shows no second that has not
    passed. The .log keeps its finer figure (runlog.format_duration)."""
    assert model.format_elapsed(seconds) == text
    state = model.RunState()
    state.elapsed = seconds
    assert state.resource_text().endswith(f"elapsed {text}")


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
    # Before a file is chosen the parameters are still checked, and the red
    # line asks for the one thing missing.
    app.cc_var.set("20")
    app.t_var.set("4")
    root.update()
    assert app.message_var.get() == "To start, choose an input file."
    assert "disabled" in app.start.state()
    app.t_var.set("18")
    root.update()
    assert app.message_var.get() == ("To start, make t smaller than 14.14 mm, "
                                     "or cc larger than 25.46 mm.")
    app.t_var.set("4")
    app.input_var.set(str(step))
    root.update()
    assert app.message_var.get() == ""
    assert app.derived_var.get().endswith("part-lattice-cc20t4.step")
    assert "disabled" not in app.start.state()

    # The one cross-parameter rule, enforced by `parse_args` rather than by
    # a copy of it living here.
    app.t_var.set("18")
    root.update()
    assert app.message_var.get().startswith("To start, make t smaller than")
    assert "disabled" in app.start.state()
    # The line is red: it is shown only while Start! is greyed out.
    assert str(app.message.cget("foreground")) == "#a11"


def test_start_is_greyed_out_while_the_output_folder_is_too_deep(tmp_path, tk_root):
    """The path limit is checked on the file the run will write, in the
    chosen folder, and named as the folder to change."""
    import os

    from latticegen2 import cli
    from latticegen2.gui.app import App

    step = tmp_path / "part.step"
    step.write_text("dummy")
    app = App(tk_root)
    app.cc_var.set("20")
    app.t_var.set("4")
    app.input_var.set(str(step))
    tk_root.update()
    assert "disabled" not in app.start.state()

    limit = cli.WINDOWS_PATH_MAX if os.name == "nt" else cli.POSIX_PATH_MAX
    deep = os.path.abspath(str(tmp_path))
    while len(deep) + 1 + cli.TEMP_PATH_RESERVE <= limit:
        deep = os.path.join(deep, "d" * 200)
    app.outdir_var.set(deep)
    tk_root.update()
    assert app.message_var.get().startswith(
        "To start, choose an output folder with a path ")
    assert "disabled" in app.start.state()

    app.outdir_var.set(str(tmp_path))
    tk_root.update()
    assert app.message_var.get() == ""
    assert "disabled" not in app.start.state()


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


# --- the verbose tick box (specification.md §3.1) ----------------------------


def test_ticking_the_box_reveals_the_run_s_own_log_and_unticking_hides_it(tk_root):
    """And the pane itself comes and goes with its contents.

    An ordinary run has nothing to put in it, so carrying an empty bordered box
    for the hour the run takes would read as something broken. The box
    appearing is itself the signal that the child said something.
    """
    from latticegen2.gui.app import App

    app = App(tk_root)
    app.panel.grid(row=7, column=0, columnspan=3, sticky="ew")
    app.state.apply(progress.parse_line(event(
        progress.LOG, msg="stage boundary: 8.1s", console=False)))
    app.state.apply(progress.parse_line(event(
        progress.LOG, msg="  Peak memory:    282.73 MB", console=True)))

    app._render_log()
    assert app.log_text.get("1.0", "end").strip() == ""
    assert not app.log_frame.grid_info()

    app.verbose_var.set(True)
    tk_root.update()
    assert app.log_text.get("1.0", "end").strip().splitlines() == [
        "stage boundary: 8.1s", "  Peak memory:    282.73 MB"]
    assert app.log_frame.grid_info()

    app.verbose_var.set(False)
    tk_root.update()
    assert app.log_text.get("1.0", "end").strip() == ""
    assert not app.log_frame.grid_info()


def test_the_tick_box_stays_live_while_a_run_is_in_flight(tmp_path, tk_root):
    """The one control that must *not* be frozen with the parameters.

    Every line crosses as an event whichever way the box is set, so verbosity is
    a filter over output already received (docs/algorithm.md §10). Disabling it
    with the rest of the fields would turn it back into a flag that has to be
    decided before the run starts — which is what the child never gets.
    """
    from latticegen2.gui.app import App

    app = App(tk_root)
    assert app.verbose_check not in app.fields
    for widget in app.fields:
        widget.state(["disabled"])
    assert "disabled" not in app.verbose_check.state()


def test_a_new_line_reaches_the_pane_after_the_backlog_saturates(tk_root):
    """The redraw is keyed on ``lines_seen``, not on the deque's length."""
    from latticegen2.gui.app import App

    app = App(tk_root)
    app.panel.grid(row=7, column=0, columnspan=3, sticky="ew")
    for i in range(model.LOG_BACKLOG):
        app.state.add_line(f"line {i}", model.OUT)
    app._render_log()
    app.state.add_line("the newest line", model.OUT)
    app._render_log()
    assert app.log_text.get("1.0", "end").strip().endswith("the newest line")


def test_a_warning_opens_the_pane_with_the_box_unticked(tk_root):
    from latticegen2.gui.app import App

    app = App(tk_root)
    app.panel.grid(row=7, column=0, columnspan=3, sticky="ew")
    app.state.apply(json.loads(event(progress.LOG, msg="a statistic", console=True)))
    app._render_log()
    assert not app.log_frame.winfo_manager()
    app.state.apply(json.loads(event(progress.LOG, msg=WARNING, console=True)))
    app._render_log()
    assert app.log_frame.winfo_manager() == "grid"
    assert app.log_text.get("1.0", "end").strip() == WARNING


def test_stderr_and_kernel_chatter_are_shown_whatever_the_box_says(tmp_path, tk_root):
    """They are not `-v` output; they are the run saying something went wrong.

    A failure's one reason line (specification.md §7) reaches the window this
    way, so it must not depend on the box being ticked.
    """
    from latticegen2.gui.app import App

    app = App(tk_root)
    app.panel.grid(row=7, column=0, columnspan=3, sticky="ew")

    class Child:
        running = False
        cancel_path = str(tmp_path / "x.cancel")
        stderr_lines = ["Intel oneMKL FATAL ERROR"]

        def __init__(self):
            import queue

            self.events = queue.Queue()

        def returncode(self):
            return 1

        def clear_cancel(self):
            pass

    app.run = Child()
    app.run.events.put(("stderr", "Intel oneMKL FATAL ERROR"))
    app._closed_seen = 2
    app._drain()
    assert not app.verbose_var.get()
    assert "Intel oneMKL FATAL ERROR" in app.log_text.get("1.0", "end")
    assert app.log_frame.grid_info()


# --- print orientation and support (specification.md §3.1) -----------------


def test_a_supported_run_lays_the_bar_out_over_the_stages_it_announced():
    """Two more stages than the plain pipeline, and the bar still ends at 1000."""
    from latticegen2.runlog import stages_for

    plain, held = stages_for(False), stages_for(True)
    assert plain == STAGES
    assert held[held.index("simplify") + 1:held.index("validate")] == (
        "overhang", "support", "pockets")
    marks = [weights.overall_permille(s, 0.0, held) for s in held]
    assert marks == sorted(marks) and marks[0] == 0
    assert weights.overall_permille(held[-1], 1.0, held) == 1000
    # The plain pipeline's table is untouched by any of this.
    assert weights.overall_permille("export", 1.0, plain) == 1000
    assert weights.overall_permille("boundary", 0.0, plain) == weights.BAND_START["boundary"]

    state = model.reduce_stream([
        event(progress.HELLO, v=1, pid=1, stages=list(held),
              output="o.step", log="o.log", workers=1),
        event(progress.STAGE_BEGIN, name="support", i=held.index("support"),
              n=len(held)),
    ])
    assert state.stage_text() == f"support  ({held.index('support') + 1} of {len(held)})"


def test_orientation_and_support_reach_the_child_and_the_parser_accepts_them(tmp_path):
    step = tmp_path / "in.step"
    step.write_text("x")
    extra = ["--orient", "30", "0", "-45", "--support", "--support-thickness", "1.2"]
    argv = runner.build_argv(input_path=str(step), output_path=str(tmp_path / "o.step"),
                             cc="10", t="1.5", cores=2, extra=extra)
    args = parse_args(argv[3:])
    assert args.orient == (30.0, 0.0, -45.0)
    assert args.support and args.support_thickness == 1.2
    assert "-v" not in argv and argv[-1] == "--progress-stream"


def test_the_view_rests_the_part_on_the_bed_and_marks_what_overhangs():
    import numpy as np

    from latticegen2 import orient
    from latticegen2.gui import view3d

    # A unit cube far from its own origin, as CAD parts usually are.
    verts = np.array([[x, y, z] for x in (100, 101) for y in (50, 51)
                      for z in (20, 21)], dtype=float)
    for angles in ((0.0, 0.0, 0.0), (30.0, 20.0, -110.0)):
        placed, _shift = view3d.place_on_bed(verts, orient.rotation_matrix(angles))
        assert placed[:, 2].min() == pytest.approx(0.0, abs=1e-9)
        assert placed[:, :2].min(axis=0) == pytest.approx(-placed[:, :2].max(axis=0))

    normals = np.array([[0, 0, 1], [0, 0, -1], [1, 0, 0]], dtype=float)
    # Upright: the lid is a ceiling of the volume, the floor and a wall are not.
    assert view3d.overhang_mask(normals, (0, 0, 0)).tolist() == [True, False, False]
    # Turned upside down, the old floor is the ceiling.
    assert view3d.overhang_mask(normals, (180, 0, 0)).tolist() == [False, True, False]
    # On its side (+X up), the wall is.
    assert view3d.overhang_mask(normals, (0, -90, 0)).tolist() == [False, False, True]
    # The bar is strict 60 deg: a ceiling tilted 30 deg from horizontal (exactly
    # 60 from vertical) is left bare, one tilted 29 deg is supported.
    tilted = np.array([[0.0, np.sin(np.radians(a)), np.cos(np.radians(a))]
                       for a in (30.0, 29.0)])
    assert view3d.overhang_mask(tilted, (0, 0, 0)).tolist() == [False, True]


def test_dragging_turns_the_part_and_reports_angles_that_reproduce_it():
    import numpy as np

    from latticegen2 import orient
    from latticegen2.gui import view3d

    R = orient.rotation_matrix((10.0, 20.0, 30.0))
    for dx, dy in ((40, 0), (0, -25), (13, 31)):
        R = view3d.drag_rotation(R, dx, dy)
        assert np.allclose(R @ R.T, np.eye(3), atol=1e-12)
        again = orient.rotation_matrix(view3d.snapped_angles(R))
        assert np.allclose(again, R, atol=5e-3)        # tenth-of-a-degree boxes
    # A horizontal drag is a turn about the bed's normal: +Z stays +Z.
    turned = view3d.drag_rotation(np.eye(3), 50, 0)
    assert np.allclose(turned @ [0, 0, 1], [0, 0, 1], atol=1e-12)


def test_the_orientation_panel_is_always_there_and_support_fields_grey_out(tmp_path, tk_root):
    from latticegen2.gui.app import App

    step = tmp_path / "part.step"
    step.write_text("x")
    app = App(tk_root)
    try:
        tk_root.update()
        # Always shown: the orientation boxes and the view are gridded at once.
        assert app.view.winfo_manager() == "grid"
        assert all(s.winfo_manager() == "grid" for s in app.orient_spins)
        assert app.thickness_spin.winfo_manager() == "grid"
        assert not hasattr(app, "overhang_spin")
        # Greyed out, not hidden, while support is off.
        assert app.thickness_spin.instate(["disabled"])

        app.input_var.set(str(step))
        app.cc_var.set("10")
        app.t_var.set("1.5")
        app.rx_var.set("30")
        tk_root.update()
        assert "--support" not in app._pending_argv()
        assert app.derived_var.get().endswith("part-lattice-cc10t1.5-rotx30y0z0.step")

        app.support_var.set(True)
        tk_root.update()
        assert app.thickness_spin.instate(["!disabled"])
        # The thickness starts at the command line's 1 mm default, stays put
        # when t changes, and is always passed with --support.
        assert app.thickness_var.get() == "1"
        assert app._pending_argv()[-3:] == ["--support", "--support-thickness", "1"]
        app.t_var.set("1.2")
        tk_root.update()
        assert app.thickness_var.get() == "1"
        assert app.derived_var.get().endswith("-rotx30y0z0-sup1.step")
        assert app.start.instate(["!disabled"])
        app.thickness_var.set("0.8")
        app.t_var.set("1.5")
        tk_root.update()
        assert app._pending_argv()[-3:] == ["--support", "--support-thickness", "0.8"]

        # Beyond a - t no support grid can hold it: the parser refuses it, its
        # message is shown, and Start! greys out.
        app.cc_var.set("10")
        app.thickness_var.set("6")
        tk_root.update()
        assert app.start.instate(["disabled"])
        assert app.message_var.get() == ("To start, lower Min thickness to at "
                                         "most 5.57 mm, or raise cc to at least 10.61 mm.")
        app.thickness_var.set("1.5")
        tk_root.update()
        assert app.start.instate(["!disabled"])

        # A drag reports angles straight into the three boxes.
        app._dragged((12.5, -3.0, 90.0))
        assert (app.rx_var.get(), app.ry_var.get(), app.rz_var.get()) == ("12.5", "-3", "90")
    finally:
        if app._preview is not None:
            app._preview.abandon()

import numpy as np


# -- the bad-geometry viewer (specification.md §3.1) ---------------------------


def _mark_event(severity, title, message, places):
    from latticegen2.runlog import mark_fields

    return {"ev": progress.MARK, "t": 1.0, **mark_fields(severity, title, message, places)}


def _cube_mesh(lo=(0.0, 0.0, 0.0), size=10.0):
    """A cube as the preview helper writes one: six faces of two triangles,
    each triangle wound outward and tagged with its face, and each face's box."""
    lo = np.asarray(lo, dtype=float)
    corners = np.array([[x, y, z] for z in (0, 1) for y in (0, 1) for x in (0, 1)],
                       dtype=float) * size + lo
    quads = [(0, 2, 3, 1), (4, 5, 7, 6), (0, 1, 5, 4), (2, 6, 7, 3), (0, 4, 6, 2), (1, 3, 7, 5)]
    tris, face = [], []
    for k, (a, b, c, d) in enumerate(quads):
        tris += [(a, b, c), (a, c, d)]
        face += [k, k]
    tris = np.array(tris, dtype=np.int64)
    face = np.array(face, dtype=np.int64)
    face_lo = np.array([corners[list(q)].min(axis=0) for q in quads])
    face_hi = np.array([corners[list(q)].max(axis=0) for q in quads])
    return corners, tris, face, face_lo, face_hi


def test_a_place_prints_exactly_the_line_it_replaces():
    """A warning built from places reads byte for byte as one built from strings,
    so the .log does not change for the viewer's sake."""
    from latticegen2.runlog import Place, box_text, point_text, warning_text

    p, lo, hi = (1, 2, 3), (0, 0, 0), (1, 2, 3)
    as_text = warning_text("t", "w", "r", [point_text(p), box_text(lo, hi),
                                           f"junction (1, 2, 3) at {point_text(p)}",
                                           "curve ends at A, vertex at B"])
    as_places = warning_text("t", "w", "r", [
        Place.point(p), Place.box(lo, hi),
        Place.point(p, f"junction (1, 2, 3) at {point_text(p)}"),
        Place.segment(p, lo, "curve ends at A, vertex at B")])
    assert str(as_places) == str(as_text)
    assert [pl.kind for pl in as_places.places] == ["point", "box", "point", "segment"]
    assert as_text.places == ()       # plain lines point nowhere


def test_a_warning_and_a_failure_keep_their_places_across_a_process_boundary():
    import pickle

    from latticegen2.errors import InputGeometryError, ParamError
    from latticegen2.runlog import Place, warning_text

    w = warning_text("t", "w", "r", [Place.point((1, 2, 3))])
    back = pickle.loads(pickle.dumps(w))
    assert back == w and back.title == "t" and back.places == w.places
    e = pickle.loads(pickle.dumps(InputGeometryError("bad", places=[Place.point((0, 0, 1))])))
    assert str(e) == "bad" and e.places[0].coords == ((0.0, 0.0, 1.0),)
    p = pickle.loads(pickle.dumps(ParamError("no", fix="raise {cc}")))
    assert p.fix == "raise {cc}"


def test_a_mark_is_read_from_its_event_and_hostile_places_are_dropped():
    from latticegen2.gui import marks
    from latticegen2.runlog import Place

    ev = _mark_event("warning", "the input file is imprecise",
                     "WARNING: the input file is imprecise: a vertex is loose. It may "
                     "fail. Inspect, in the part's coordinates:\n  at [1.000, 2.000, 3.000] mm",
                     [Place.point((1, 2, 3), "at [1.000, 2.000, 3.000] mm")])
    ev["places"] += [{"kind": "point", "p": [1, 2]}, {"kind": "blob"}, "junk",
                     {"kind": "box", "lo": [0, 0, "x"], "hi": [1, 1, 1]},
                     {"kind": "point", "p": [float("nan"), 0, 0]}]
    mark = marks.from_event(ev)
    assert mark is not None and len(mark.places) == 1
    assert mark.sentence() == "a vertex is loose. It may fail."
    assert mark.colour == marks.WARNING_COLOUR and not mark.is_error
    assert marks.from_event({"ev": "mark", "severity": "info", "places": ev["places"]}) is None
    assert marks.from_event({**ev, "places": [{"kind": "blob"}]}) is None
    text = marks.detail_text(mark, mark.places[0])
    assert text.startswith("Warning: the input file is imprecise")
    assert "at [1.000, 2.000, 3.000] mm" in text and "a vertex is loose" in text


def test_marks_reach_the_run_state_and_say_there_is_bad_geometry():
    from latticegen2.runlog import Place

    s = model.RunState()
    s.apply(json.loads(event(progress.HELLO, stages=list(STAGES))))
    assert not s.has_bad_geometry
    s.apply(_mark_event("error", "the run failed", "boom", [Place.point((0, 0, 0))]))
    s.apply({"ev": progress.MARK, "t": 2.0, "severity": "warning", "title": "x",
             "message": "y", "places": [{"kind": "nonsense"}]})
    assert s.has_bad_geometry and len(s.marks) == 1 and s.marks[0].is_error


def test_a_failures_face_is_found_among_the_preview_faces_by_box():
    from latticegen2.gui import marks

    _v, _t, _f, face_lo, face_hi = _cube_mesh()
    # The top face, its box tightened a little as OCCT does once it is meshed.
    assert marks.match_face((0.001, 0.0, 10.0), (9.999, 10.0, 10.0), face_lo, face_hi) == 1
    assert marks.match_face((40, 40, 40), (50, 50, 50), face_lo, face_hi) is None
    assert marks.match_face((0, 0, 0), (1, 1, 1), None, None) is None


def test_thinning_a_preview_keeps_every_face():
    from latticegen2.preview import thin_keeping_faces

    owner = np.repeat(np.arange(50), [1] + [99] * 49)
    keep = thin_keeping_faces(owner, 200)
    assert len(keep) <= 200 and set(owner[keep]) == set(range(50))
    assert np.all(np.diff(keep) > 0)
    assert np.array_equal(thin_keeping_faces(owner[:10], 200), np.arange(10))


def test_the_viewer_camera_zooms_about_the_cursor_and_fits_the_part():
    from latticegen2.gui import view3d

    centre, cursor = (400.0, 300.0), (520.0, 260.0)
    scale, pan = 2.0, (15.0, -8.0)
    # The point under the cursor, in the scene's own units.
    xy = ((cursor[0] - centre[0] - pan[0]) / scale, (cursor[1] - centre[1] - pan[1]) / scale)
    new, (px, py) = view3d.zoom_about(scale, pan, 1.6, cursor, centre)
    assert new == pytest.approx(3.2)
    assert centre[0] + px + new * xy[0] == pytest.approx(cursor[0])
    assert centre[1] + py + new * xy[1] == pytest.approx(cursor[1])
    assert view3d.fit_scale(10.0, 800, 600) == pytest.approx(600 * 0.5 / (10.0 * 1.15))
    pts = np.array([[1.0, 2.0, 3.0], [0.0, 0.0, 0.0]])
    assert np.array_equal(view3d.view_points(pts, np.eye(3), np.array([5.0, 5.0, 5.0])), pts)


def test_the_overhang_tint_is_blue_in_both_views():
    from latticegen2.gui import marks, orientview

    r, g, b = orientview.OVERHANG_COLOUR
    assert b > 0.8 and r < 0.3                       # blue
    assert orientview.OVERHANG_COLOUR not in (marks.ERROR_COLOUR, marks.WARNING_COLOUR)


def test_the_viewer_cap_is_the_helpers():
    from latticegen2 import preview
    from latticegen2.gui import app

    assert app.VIEWER_TRIANGLES == preview.VIEWER_TRIANGLES


def _viewer(tk_root, marks_=(), orient=(0, 0, 0)):
    from latticegen2.gui.badgeom import BadGeometryViewer

    v = BadGeometryViewer(tk_root, orient=orient, mesh=_cube_mesh(), marks=marks_)
    v.geometry("800x600+40+40")
    for _ in range(5):
        tk_root.update()
    v.fit()
    tk_root.update()
    return v


def _sample_marks():
    from latticegen2.gui import marks
    from latticegen2.runlog import Place

    warn = marks.from_event(_mark_event(
        "warning", "corner", "WARNING: corner: it moved. It may fail. Inspect, in the "
        "part's coordinates:\n  at [10.000, 10.000, 10.000] mm",
        [Place.point((10, 10, 10), "at [10.000, 10.000, 10.000] mm"),
         Place.box((0, 0, 0), (5, 5, 5))]))
    err = marks.from_event(_mark_event(
        "error", "the run failed", "Surface tessellation for face 1 failed",
        [Place.face((5, 5, 10), (0, 0, 10), (10, 10, 10), "face 1")]))
    return [warn, err]


def test_the_viewer_fills_its_window_and_overlays_its_controls(tk_root):
    v = _viewer(tk_root, _sample_marks())
    try:
        assert v.canvas.winfo_manager() == "pack"
        assert v.overlay.winfo_manager() == "place"
        assert v.overlay.place_info()["anchor"] == "ne"
        # The controls are described in the top-left corner.
        (item,) = v.canvas.find_withtag("controls")
        x0, y0, _x1, _y1 = v.canvas.bbox(item)
        assert x0 < 20 and y0 < 20
        assert "Pan: middle drag" in v.canvas.itemcget(item, "text")
        # No print bed: nothing but the mesh, the marks, the axes and the text.
        tags = {t for i in v.canvas.find_all() for t in v.canvas.gettags(i)}
        assert not any("bed" in t for t in tags)
        # The issue list starts open, folds to its title, and the overhang box
        # above it stays either way.
        assert v.issues_expanded and v.issues_body.winfo_manager() == "pack"
        v.toggle_issues()
        tk_root.update()
        assert not v.issues_expanded and not v.issues_body.winfo_manager()
        assert v.overhangs_box.winfo_manager() == "pack"
        v.toggle_issues()
        assert v.issues_body.winfo_manager() == "pack"
        # Errors are listed first.
        first = v.tree.get_children()[0]
        assert v.tree.item(first, "text").startswith("Error")
        # Resizing re-fits the viewport to the window.
        v.geometry("1000x700")
        for _ in range(5):
            tk_root.update()
        assert abs(v.canvas.winfo_width() - 1000) <= 2
        assert abs(v.canvas.winfo_height() - 700) <= 2
    finally:
        v.close()


def test_the_viewer_marks_in_red_and_yellow_and_opens_details_on_a_click(tk_root):
    from latticegen2.gui import marks
    from latticegen2.gui.orientview import _hex

    v = _viewer(tk_root, _sample_marks())
    try:
        red, yellow = _hex(marks.ERROR_COLOUR), _hex(marks.WARNING_COLOUR)
        fills = {v.canvas.itemcget(i, "fill") for i in v.canvas.find_withtag("place:0:0")}
        assert fills == {yellow}
        # The failure's face is the cube's top face, tinted red, plus its pin.
        assert v._face_match == {(1, 0): 1}
        assert red in {v.canvas.itemcget(i, "fill") for i in v.canvas.find_withtag("place:1:0")}
        # Click the warning's point: its detail box opens, bordered in yellow.
        (sx, sy), = v.screen_of(v.to_bed((10, 10, 10)))[0]
        v.click(sx, sy)
        assert v.detail == (0, 0)
        boxes = [i for i in v.canvas.find_withtag("detail")
                 if v.canvas.type(i) == "rectangle"]
        assert boxes and v.canvas.itemcget(boxes[0], "outline") == yellow
        text = [v.canvas.itemcget(i, "text") for i in v.canvas.find_withtag("detail")
                if v.canvas.type(i) == "text"][0]
        assert "at [10.000, 10.000, 10.000] mm" in text and "it moved" in text
        # A click on nothing closes it.
        v.click(2, 590)
        assert v.detail is None
        # Picking from the list centres that place and opens its box.
        v.tree.selection_set("m1_0")
        tk_root.update()
        assert v.detail == (1, 0)
        (cx, cy), = v.screen_of(v.to_bed((5, 5, 10)))[0]
        assert abs(cx - v.canvas.winfo_width() / 2) < 1.0
        assert abs(cy - v.canvas.winfo_height() / 2) < 1.0
    finally:
        v.close()


def test_the_viewer_turns_pans_zooms_and_fits_without_touching_the_orientation(tk_root):
    from latticegen2.gui.badgeom import CLICK_SLOP

    v = _viewer(tk_root, _sample_marks(), orient=(30, 0, 0))
    try:
        c = v.canvas
        R = v._R.copy()
        # A click that hardly moves is a click, not a rotation.
        c.event_generate("<ButtonPress-1>", x=300, y=300)
        c.event_generate("<B1-Motion>", x=300 + CLICK_SLOP, y=300)
        c.event_generate("<ButtonRelease-1>", x=300 + CLICK_SLOP, y=300)
        assert np.array_equal(v.view_rotation, np.eye(3))
        # A left drag turns the view as a drag turns the part in the small view.
        c.event_generate("<ButtonPress-1>", x=300, y=300)
        c.event_generate("<B1-Motion>", x=340, y=310)
        c.event_generate("<ButtonRelease-1>", x=340, y=310)
        from latticegen2.gui import view3d

        expect = view3d.drag_rotation(np.eye(3), 40, 10)
        assert np.allclose(v.view_rotation, expect)
        assert np.array_equal(v._R, R)               # the run's orientation is untouched
        # A middle drag pans.
        scale, pan = v.camera
        c.event_generate("<ButtonPress-2>", x=100, y=100)
        c.event_generate("<B2-Motion>", x=130, y=80)
        c.event_generate("<ButtonRelease-2>", x=130, y=80)
        assert v.camera[1] == (pan[0] + 30, pan[1] - 20)
        # The wheel zooms, about the cursor.
        c.event_generate("<MouseWheel>", x=200, y=200, delta=120)
        assert v.camera[0] == pytest.approx(scale * 1.15)
        # F fits again: scale back, pan gone, rotation kept.
        c.focus_force()
        tk_root.update()
        c.event_generate("<KeyPress-f>")
        tk_root.update()
        assert v.camera == (pytest.approx(scale), (0.0, 0.0))
        assert np.allclose(v.view_rotation, expect)
    finally:
        v.close()


def test_the_overhang_toggle_starts_off_and_tints_blue(tk_root):
    from latticegen2.gui.orientview import OVERHANG_COLOUR, _hex

    v = _viewer(tk_root, orient=(0, 0, 0))
    try:
        def blues():
            # The cube's top face overhangs (it is a ceiling); shaded blues only.
            return [i for i in v.canvas.find_all()
                    if v.canvas.type(i) == "polygon"
                    and v.canvas.itemcget(i, "fill").startswith("#")
                    and int(v.canvas.itemcget(i, "fill")[5:7], 16)
                    > int(v.canvas.itemcget(i, "fill")[1:3], 16) + 80]
        assert v.overhangs_var.get() is False and not blues()
        v.overhangs_var.set(True)
        v.redraw()
        assert blues()
        assert OVERHANG_COLOUR == (0.18, 0.44, 0.92) and _hex(OVERHANG_COLOUR)
    finally:
        v.close()


def test_show_bad_geometry_appears_with_the_first_mark_and_opens_one_viewer(tk_root):
    from latticegen2.gui.app import App
    from latticegen2.runlog import Place

    class FakeRun:
        def __init__(self, lines):
            import queue

            self.events = queue.Queue()
            for line in lines:
                self.events.put(("event", json.loads(line) if isinstance(line, str) else line))
            self.running = True
            self.stderr_lines = []

        def returncode(self):
            return 0

        def clear_cancel(self):
            pass

    app = App(tk_root)
    try:
        app.run = FakeRun([event(progress.HELLO, stages=list(STAGES), output="o", log="l")])
        app._started = 0.0
        app.panel.grid(row=9, column=0, columnspan=3, sticky="ew")
        app._drain()
        tk_root.update()
        assert not app.bad_button.winfo_manager()
        # Mid-run, the first mark brings the button.
        app.run.events.put(("event", _mark_event("warning", "x", "WARNING: x: y. z.",
                                                 [Place.point((1, 2, 3))])))
        app._drain()
        tk_root.update()
        assert app.bad_button.winfo_manager() == "grid"
        app._show_bad_geometry()
        viewer = app._viewer
        assert viewer is not None and len(viewer._marks) == 1
        app._show_bad_geometry()
        assert app._viewer is viewer                    # raised, not a second one
        # A failure's red mark joins the open viewer live.
        app.run.events.put(("event", _mark_event("error", "the run failed", "boom",
                                                 [Place.point((0, 0, 0))])))
        app._drain()
        assert len(viewer._marks) == 2
        # Closed with the window's x, it is forgotten.
        viewer.protocol("WM_DELETE_WINDOW")
        viewer.close()
        assert app._viewer is None
    finally:
        app.run = None
        app._close_viewer()
        if app._preview is not None:
            app._preview.abandon()


def test_a_parameter_failure_shows_no_bad_geometry_button(tk_root):
    from latticegen2.gui.app import App

    app = App(tk_root)
    try:
        app.state.apply(json.loads(event(progress.RESULT, status="failed", exit=2,
                                         reason="-t too large")))
        app._refresh()
        tk_root.update()
        assert not app.bad_button.winfo_manager()
    finally:
        if app._preview is not None:
            app._preview.abandon()
