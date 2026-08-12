"""Process entry point: argument handling, failure reporting, exit codes.

Every non-zero exit prints exactly one human-readable reason line naming what
failed (specification.md §7). A cancelled run prints ``CANCELLED:`` rather than
``FAILED:`` and no traceback — it did what the user asked, it did not
malfunction.
"""

from __future__ import annotations

import sys

from .cli import HelpRequested, USAGE, parse_args, preflight_checks
from .errors import CancelledError, LatticeGenError
from .runlog import RunLog


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
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

    if args.background:
        from .boundary import _set_background_priority

        _set_background_priority()

    from .pipeline import run_pipeline

    rl = RunLog(args.log_path, verbose=args.verbose).open()
    try:
        rl.header(args.as_dict())
        stats = run_pipeline(args, rl)
        rl.summary(args.as_dict(), stats)
        return 0
    except KeyboardInterrupt:
        rl.always("CANCELLED: interrupted by user (Ctrl+C). Intermediate files kept.")
        return CancelledError.exit_code
    except LatticeGenError as exc:
        # `line(console=False)` rather than `always`: the reason goes to the log
        # here and to stderr below. `always` would echo it to stdout as well and
        # the user would see the same line twice.
        rl.line(f"FAILED: {exc}", console=False)
        print(f"FAILED: {exc}", file=sys.stderr)
        return exc.exit_code
    except Exception as exc:  # unexpected: report it, but still as one line first
        rl.line(f"FAILED: unexpected {type(exc).__name__}: {exc}", console=False)
        print(f"FAILED: unexpected {type(exc).__name__}: {exc}", file=sys.stderr)
        raise
    finally:
        rl.close()


if __name__ == "__main__":
    sys.exit(main())
