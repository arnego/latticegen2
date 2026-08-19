"""Launching, watching and stopping a run. No Tk in this module.

The pipeline runs as a **child process**, never in the window's own process, and
the reason that is not negotiable is the GIL: G7 measured OCP holding it across
``ShapeUpgrade_UnifySameDomain`` (`tools/prototypes/RESULTS.md`, quoted at
:func:`latticegen2.pipeline._unify`), and ``simplify`` plus ``export`` are 47 %
of a rehearsal run. In-process, the Tk event loop could not acquire the GIL for
that long and the window would be "Not Responding" for half the run.

Three more consequences fall out of the same choice, all of them wanted: the
child's ``__main__`` stays ``src/main.py``, so ``spawn``'s re-import contract is
exactly what it has always been; a 19 GB peak and any OCCT access violation stay
on the far side of a process boundary; and ``set_background_priority`` demotes
the pipeline without demoting the user interface.
"""

from __future__ import annotations

import os
import queue
import subprocess
import sys
import threading

from .. import progress

#: How long Stop waits for the child to shut down cleanly before killing it.
#: specification.md §3 asks for exactly this shape — orderly first, forced only
#: if the run does not respond "within a short grace period".
GRACE_SECONDS = 10.0


def child_python() -> str:
    """The interpreter to run the pipeline with.

    Under ``pythonw.exe`` — which is how the launcher starts the window, so no
    console lingers behind it — ``sys.executable`` is the *windowed* interpreter.
    The child needs the console one: it is not drawing anything, and a windowed
    interpreter has no usable standard streams, which is the whole channel this
    module depends on.
    """
    exe = sys.executable
    base = os.path.basename(exe).lower()
    if base.startswith("pythonw"):
        console = os.path.join(os.path.dirname(exe),
                               os.path.basename(exe).replace("pythonw", "python", 1))
        if os.path.isfile(console):
            return console
    return exe


class LaunchError(RuntimeError):
    """The child could not be started. Carries a line fit to show the user."""


def main_script() -> str:
    """``src/main.py``, resolved from this package's own location.

    Deliberately not the installed console script and not ``-m latticegen2``:
    boundary workers re-import ``__main__`` under ``spawn``, and ``src/main.py``
    is what puts ``src/`` on their path (latticegen2.bat, specification.md §2).
    Launching anything else would be a different program to the one every other
    invocation of this tool runs.
    """
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    path = os.path.join(here, "main.py")
    if not os.path.isfile(path):
        # Reachable, and worth naming rather than letting `Popen` raise a bare
        # FileNotFoundError into a Tk callback where nothing would show it: a
        # `pip install .` puts this package under site-packages, where there is
        # no `src/` above it, and `pyproject.toml`'s console script reaches the
        # window with no arguments.
        raise LaunchError(
            f"cannot find {path}. The graphical front-end runs the pipeline from "
            f"a checkout or a release bundle; an installed package has no "
            f"src/main.py beside it. Use the launcher in the bundle, or run the "
            f"command line instead."
        )
    return path


def build_argv(*, input_path: str, output_path: str, cc: str, t: str,
               cores: int) -> list[str]:
    """The exact command line the child is given.

    ``-v`` is never passed. Every log line crosses as an event carrying whether
    a CLI run would have shown it, so verbosity is something the window filters
    on its own side and can change while the run is going.
    """
    return [
        child_python(), "-u", main_script(),
        "-i", input_path,
        "-cc", cc,
        "-t", t,
        "-o", output_path,
        "--cores", str(cores),
        "--progress-stream",
    ]


class Run:
    """One child process, its two reader threads, and the queue they feed.

    Nothing here touches a widget. The window drains :attr:`events` from its own
    event loop, which is what keeps every Tk call on the main thread — the
    classic way a Tk application corrupts itself is a background thread updating
    a widget directly.
    """

    def __init__(self, argv: list[str], log_path: str, cwd: str | None = None):
        self.argv = argv
        self.cancel_path = progress.cancel_sentinel_path(log_path)
        self.events: queue.Queue = queue.Queue()
        self.stderr_lines: list[str] = []
        self._cancelling = False
        creation = 0
        popen_extra: dict = {}
        if os.name == "nt":
            # No console window for the child. The parent is a GUI process and
            # a flashing console on every run would be the most visible thing
            # about it.
            creation = getattr(subprocess, "CREATE_NO_WINDOW", 0)
        else:
            # Its own session, so a force-stop can take the whole worker tree
            # rather than orphaning the pool.
            popen_extra["start_new_session"] = True

        # **stdin must be the null device, not a pipe.** Handing the child an
        # anonymous pipe for file descriptor 0 stalls its worker pool the first
        # time `spawn` runs: measured on this machine, `classify` took 0.2 s
        # with the null device and never finished with a pipe — with or without
        # `CREATE_NO_WINDOW`, and whether or not this end was kept open. That is
        # why Stop is a file rather than a line on stdin
        # (`progress.cancel_sentinel_path`).
        self.proc = subprocess.Popen(
            argv, cwd=cwd,
            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, encoding="utf-8", errors="replace", bufsize=1,
            creationflags=creation, **popen_extra,
        )
        self._reader(self.proc.stdout, self._on_stdout)
        self._reader(self.proc.stderr, self._on_stderr)

    # -- plumbing ----------------------------------------------------------

    def _reader(self, stream, handle) -> None:
        def pump():
            try:
                for line in stream:
                    handle(line)
            except Exception:
                pass
            finally:
                self.events.put(("closed", None))

        threading.Thread(target=pump, daemon=True).start()

    def _on_stdout(self, line: str) -> None:
        event = progress.parse_line(line)
        if event is not None:
            self.events.put(("event", event))
        elif line.strip():
            # Not an event, and not a reason to stop reading: OCCT writes to
            # file descriptor 1 below `sys.stdout`, so kernel chatter can land
            # here mid-stream.
            self.events.put(("raw", line.rstrip("\n")))

    def _on_stderr(self, line: str) -> None:
        if line.strip():
            self.stderr_lines.append(line.rstrip("\n"))
            self.events.put(("stderr", line.rstrip("\n")))

    # -- control -----------------------------------------------------------

    @property
    def running(self) -> bool:
        return self.proc.poll() is None

    @property
    def cancelling(self) -> bool:
        return self._cancelling

    def cancel(self) -> None:
        """Ask the run to stop the way Ctrl+C does.

        The child turns this into ``KeyboardInterrupt``, so it takes the whole
        existing shutdown: workers stopped in order, one ``CANCELLED`` line,
        exit 130, and the temp folder **left in place** for analysis
        (specification.md §3, §4.4).

        It is not instant, and the window says so rather than pretending: the
        interrupt is checked between bytecodes, so inside a long kernel call it
        lands when that call returns. Exactly what Ctrl+C does on the command
        line — a property of the pipeline, not of this channel.
        """
        self._cancelling = True
        try:
            with open(self.cancel_path, "w", encoding="utf-8") as fh:
                fh.write("cancelled from the latticegen2 window\n")
        except OSError:
            # Unwritable output folder: the run could not have started, and the
            # grace period below force-stops it either way.
            pass

    def force_stop(self) -> None:
        """Kill the child *and its workers* after the grace period.

        The tree matters. Under ``spawn`` the pool's workers are children of the
        master, and on Windows killing the master leaves them running — up to
        ``--cores`` orphaned OCCT processes, each potentially holding gigabytes.
        ``taskkill /T`` and a POSIX process-group kill are what make Stop mean
        stopped.
        """
        try:
            if os.name == "nt":
                subprocess.run(
                    ["taskkill", "/F", "/T", "/PID", str(self.proc.pid)],
                    capture_output=True,
                    creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
                )
            else:
                import signal
                os.killpg(os.getpgid(self.proc.pid), signal.SIGKILL)
        except Exception:
            try:
                self.proc.kill()
            except Exception:
                pass

    def returncode(self) -> int | None:
        return self.proc.poll()

    def clear_cancel(self) -> None:
        """Remove the sentinel once the run is over.

        The child removes it as it acts on it, so this only matters when the run
        ended some other way first — a failure, or a force-stop — and it stops a
        stale file being there for the next run to trip over.
        """
        try:
            os.remove(self.cancel_path)
        except OSError:
            pass
