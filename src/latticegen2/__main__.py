"""Process entry point: argument handling, failure reporting, exit codes.

Every non-zero exit prints exactly one human-readable reason line naming what
failed (specification.md §7). A cancelled run prints ``CANCELLED:`` rather than
``FAILED:`` and no traceback — it did what the user asked, it did not
malfunction.

Under ``--progress-stream`` the same run is additionally reported as
machine-readable events for the graphical front-end (specification.md §3.1,
docs/algorithm.md §10), and this module is also where that run can be cancelled
from — see :func:`_watch_for_cancel`.
"""

from __future__ import annotations

import os
import sys

from . import progress
from .cli import HelpRequested, USAGE, parse_args, preflight_checks
from .errors import CancelledError, LatticeGenError
from .runlog import STAGES, RunLog


#: How often the run looks for its cancel sentinel. Fast enough that Stop feels
#: immediate, and utterly negligible beside a stage measured in minutes.
CANCEL_POLL_S = 0.25


def _watch_for_cancel(log_path: str) -> None:
    """Turn the appearance of the cancel sentinel into ``Ctrl+C``.

    :func:`latticegen2.progress.cancel_sentinel_path` records why the channel is
    a file rather than a signal or a line on stdin; both alternatives were tried
    and one of them was measured breaking the worker pool outright.

    ``_thread.interrupt_main`` raises ``KeyboardInterrupt`` in the main thread,
    so a cancelled run takes the path it has always taken: workers torn down,
    one ``CANCELLED:`` line, exit 130, temp folder kept (specification.md §3).
    Nothing new had to be threaded through the pipeline for it.

    **A sentinel already on disk when the run starts is removed, not obeyed.**
    It is left over from something else — a previous run that was cancelled and
    whose window closed before it could tidy up — and honouring it would kill a
    fresh run instantly for a reason nobody could see.

    One honest limit, which the front-end states rather than hides: the
    interrupt is checked between bytecodes, so inside a long GIL-holding OCCT
    call it fires when that call returns and not before. That is exactly what
    Ctrl+C does on the command line today, so it is a property of the pipeline
    rather than a shortcoming of this channel.
    """
    import _thread
    import threading
    import time

    path = progress.cancel_sentinel_path(log_path)
    try:
        os.remove(path)
    except OSError:
        pass

    def watch() -> None:
        while True:
            if os.path.exists(path):
                try:
                    os.remove(path)
                except OSError:
                    pass
                _thread.interrupt_main()
                return
            time.sleep(CANCEL_POLL_S)

    threading.Thread(target=watch, daemon=True, name="cancel-watch").start()


def _display_available() -> bool:
    """Whether a window could be opened at all."""
    if sys.platform == "win32":
        return True
    return bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)

    # The graphical front-end, dispatched before `parse_args` (which requires
    # `-i`) and before `set_background_priority` below (the window must stay at
    # normal priority; only the pipeline it launches is demoted).
    #
    # **Zero arguments open a window only where one can exist.** Without the
    # display gate this would silently change what a bare `latticegen2` does on
    # a headless machine — today it exits 2 with usage, and a script or a CI job
    # relying on that would instead hang against a window that can never appear.
    # An explicit `--gui` always tries, and says so in one line if it cannot.
    if not argv or argv == ["--gui"]:
        if argv or _display_available():
            from .gui import run  # deferred: `spawn` children never import tkinter

            return run()
        print(f"FAILED: -i/--input is required.\n\n{USAGE}", file=sys.stderr)
        return 2

    try:
        args = parse_args(argv)
    except HelpRequested:
        print(USAGE)
        return 0
    except LatticeGenError as exc:
        print(f"FAILED: {exc}", file=sys.stderr)
        return exc.exit_code

    try:
        preflight_checks(args)
    except LatticeGenError as exc:
        # Deliberately before the log file is opened: a run rejected at this
        # point must leave no .step and no .log behind.
        print(f"FAILED: {exc}", file=sys.stderr)
        return exc.exit_code

    # Unconditional, and here rather than in `parse_args`, so a run rejected for
    # a bad parameter never touches this process's scheduling at all. Each
    # worker sets its own once, from the pool initializer.
    from .parallel import set_background_priority

    set_background_priority()

    from .pipeline import run_pipeline

    rl = RunLog(args.log_path, verbose=args.verbose).open()

    def finish(status: str, code: int, reason: "str | None") -> int:
        if rl.emit is not None:
            rl.emit(
                progress.RESULT, status=status, exit=code, reason=reason,
                output=args.output, log=args.log_path, tmpdir=rl.tmpdir,
            )
        return code

    try:
        if args.progress_stream:
            # Inside the guarded region, not before it. A cancellation can
            # arrive the instant the watcher starts, and raised above this
            # `try` it would escape `main` as an uncaught KeyboardInterrupt —
            # a traceback on a run that did exactly what it was asked, which
            # is the one thing specification.md §3 says a cancelled run must
            # not produce.
            #
            # After the log is open, so the emitter's clock agrees with the
            # log's, and after `preflight_checks`, so a run rejected there
            # still leaves nothing behind and reports itself as it always has.
            rl.emit = progress.NdjsonEmitter(sys.stdout, t0=rl.t0)
            rl.emit(
                progress.HELLO, v=progress.SCHEMA_VERSION, pid=os.getpid(),
                stages=list(STAGES), output=args.output, log=args.log_path,
                workers=args.workers,
            )
            _watch_for_cancel(args.log_path)
        rl.header(args.as_dict())
        stats = run_pipeline(args, rl)
        rl.summary(args.as_dict(), stats)
        return finish(progress.OK, 0, None)
    except KeyboardInterrupt:
        rl.always("CANCELLED: interrupted by user (Ctrl+C). Intermediate files kept.")
        return finish(progress.CANCELLED, CancelledError.exit_code,
                      "interrupted by user. Intermediate files kept.")
    except LatticeGenError as exc:
        # `line(console=False)` rather than `always`: the reason goes to the log
        # here and to stderr below. `always` would echo it to stdout as well and
        # the user would see the same line twice.
        rl.line(f"FAILED: {exc}", console=False)
        print(f"FAILED: {exc}", file=sys.stderr)
        return finish(progress.FAILED, exc.exit_code, str(exc))
    except Exception as exc:  # unexpected: report it, but still as one line first
        rl.line(f"FAILED: unexpected {type(exc).__name__}: {exc}", console=False)
        print(f"FAILED: unexpected {type(exc).__name__}: {exc}", file=sys.stderr)
        finish(progress.FAILED, 1, f"unexpected {type(exc).__name__}: {exc}")
        raise
    finally:
        rl.close()


if __name__ == "__main__":
    sys.exit(main())
