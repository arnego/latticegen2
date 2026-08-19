"""Process entry point: failure reporting and exit codes.

specification.md §7 and docs/algorithm.md §10 require *exactly one* human-readable
reason line per non-zero exit. This covers the branch that reports a failure
raised after the log file is already open, which is where the duplication of
issue #6 lived.
"""

import os

from latticegen2.__main__ import main


def _run(tmp_path, capsys, argv):
    code = main(argv)
    captured = capsys.readouterr()
    return code, captured


def test_a_failure_after_the_log_opens_reports_one_reason_line(tmp_path, capsys):
    """One line on stderr, one in the log, and nothing echoed to stdout.

    `RunLog.always` writes the log *and* prints to stdout, so pairing it with a
    `print` to stderr showed the user the same FAILED line twice — the second
    symptom reported in issue #6.
    """
    bad = tmp_path / "not-really.step"
    bad.write_text("this is not a STEP file\n", encoding="utf-8")
    out = tmp_path / "out.step"

    code, captured = _run(
        tmp_path, capsys, ["-i", str(bad), "-o", str(out), "-cc", "10", "-t", "1"]
    )

    assert code == 3  # InputGeometryError
    stderr_reasons = [ln for ln in captured.err.splitlines() if ln.startswith("FAILED:")]
    stdout_reasons = [ln for ln in captured.out.splitlines() if ln.startswith("FAILED:")]
    assert len(stderr_reasons) == 1, captured.err
    assert stdout_reasons == [], captured.out

    log = tmp_path / "out.log"
    assert log.is_file()
    logged = [ln for ln in log.read_text(encoding="utf-8").splitlines() if "FAILED:" in ln]
    assert len(logged) == 1, logged


def test_a_failure_before_the_log_opens_writes_no_files(tmp_path, capsys):
    """A run rejected in preflight must leave no .step and no .log behind."""
    missing = tmp_path / "absent.step"
    out = tmp_path / "out.step"

    code, captured = _run(
        tmp_path, capsys, ["-i", str(missing), "-o", str(out), "-cc", "10", "-t", "1"]
    )

    assert code == 3
    assert len([ln for ln in captured.err.splitlines() if ln.startswith("FAILED:")]) == 1
    assert os.listdir(tmp_path) == []


# --- the graphical front-end's dispatch (specification.md §3.1) --------------


def test_no_arguments_on_a_headless_machine_keeps_the_old_behaviour(monkeypatch, capsys):
    """The gate that stops this being a silent breaking change.

    A bare ``latticegen2`` has always exited 2 with usage. Opening a window
    instead would change that for every script and CI job that relies on it, and
    on a machine with no display it could not open one anyway — so the window is
    offered only where a window can exist.
    """
    monkeypatch.setattr("latticegen2.__main__._display_available", lambda: False)
    code = main([])
    captured = capsys.readouterr()
    assert code == 2
    assert captured.err.count("FAILED:") == 1
    assert "-i, --input" in captured.err


def test_no_arguments_with_a_display_opens_the_window(monkeypatch):
    monkeypatch.setattr("latticegen2.__main__._display_available", lambda: True)
    import latticegen2.gui as gui

    monkeypatch.setattr(gui, "run", lambda: 0)
    assert main([]) == 0


def test_gui_explicitly_tries_even_with_no_display(monkeypatch):
    """``--gui`` is a request, not a hint: it reports why it could not, rather
    than quietly printing usage as if the user had asked for the CLI."""
    monkeypatch.setattr("latticegen2.__main__._display_available", lambda: False)
    import latticegen2.gui as gui

    monkeypatch.setattr(gui, "run", lambda: 0)
    assert main(["--gui"]) == 0


def test_a_window_that_cannot_start_reports_one_line_and_no_traceback(monkeypatch, capsys):
    """Under ``pythonw.exe`` there is no stderr, so this is also the path that
    puts the message in a dialog — the alternative is a double-click that does
    nothing at all."""
    import latticegen2.gui as gui

    def explode():
        raise ImportError("No module named '_tkinter'")

    monkeypatch.setattr(gui, "_report_startup_failure", gui._report_startup_failure)
    monkeypatch.setattr("latticegen2.gui.app.launch", explode, raising=False)
    monkeypatch.setattr("latticegen2.__main__._display_available", lambda: True)
    code = main(["--gui"])
    captured = capsys.readouterr()
    assert code == 1
    assert captured.err.count("FAILED:") == 1
    assert "_tkinter" in captured.err
    assert "Traceback" not in captured.err
