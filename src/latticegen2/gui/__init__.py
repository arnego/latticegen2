"""The graphical front-end (specification.md §3.1).

Started by the launchers with no arguments, or by ``--gui`` explicitly. It adds
no capability the command line lacks: it collects the same parameters, validates
them with the same :mod:`latticegen2.cli` code, and runs the same
``src/main.py`` as a child process. What it adds is the ability to *watch* an
hour-long run — which stage is running, how far into it, and how much memory it
is holding.

Import discipline, which matters more here than the size of the module suggests:
nothing under this package may import :mod:`latticegen2.pipeline`,
:mod:`latticegen2.occ` or ``OCP``. The window has to appear immediately and must
survive a run that crashes, and both follow from the geometry kernel living
entirely in the child. Only :mod:`latticegen2.cli` (and through it ``numpy``) is
loaded here, so validation is shared rather than reimplemented.
"""

from __future__ import annotations


def run() -> int:
    """Open the window and return the process exit code.

    Failures here are reported by the only channel a windowed process has. Under
    ``pythonw.exe`` there is no stderr for a message to land in, so a Tk that
    cannot start would otherwise be a double-click that does nothing at all —
    the least diagnosable failure this tool could have.
    """
    try:
        from .app import launch

        return launch()
    except BaseException as exc:  # noqa: BLE001 - the last line of defence
        _report_startup_failure(exc)
        return 1


def _report_startup_failure(exc: BaseException) -> None:
    """Say what went wrong, on whichever channel this process actually has.

    A dialog box is the *fallback*, not the default, and the condition is what
    makes that safe: under ``pythonw.exe`` — how the launcher starts the window —
    there is no console, and CPython sets ``sys.stderr`` to ``None``. That is
    exactly when a printed message would vanish and a modal dialog is the only
    way to reach anyone.

    When stderr does exist, printing is enough and the dialog must **not**
    appear: it blocks until somebody clicks it, which in a test run or a
    scripted ``--gui`` on a build machine means a process that never exits.
    """
    import sys

    message = (
        f"latticegen2 could not start its window.\n\n"
        f"{type(exc).__name__}: {exc}\n\n"
        f"The command line is unaffected: run latticegen2 with parameters "
        f"(-i, -cc, -t) for the usual behaviour, or with --help for usage."
    )
    stderr = getattr(sys, "stderr", None)
    if stderr is not None:
        print(f"FAILED: {message}", file=stderr)
        return
    if sys.platform == "win32":
        try:
            import ctypes

            ctypes.windll.user32.MessageBoxW(None, message, "latticegen2", 0x10)
        except Exception:
            pass
