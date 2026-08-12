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
