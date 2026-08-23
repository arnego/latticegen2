"""The window itself (specification.md §3.1).

Small and utilitarian on purpose: it is meant to sit in a corner of the screen
for the hour a large part takes, so it is fixed-width, grows only while a run is
in flight, and shows nothing it cannot actually measure.

The split with :mod:`latticegen2.gui.model` is deliberate — that module decides
*what* to show and is fully unit-tested without a display; this one only draws
it and drives the event loop. The only rule that matters here is that **no Tk
call is ever made off the main thread**: the reader threads in
:mod:`latticegen2.gui.runner` push into a queue, and this module drains it from
an ``after`` callback.
"""

from __future__ import annotations

import os
import subprocess
import sys
import time
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from .. import __version__, progress
from ..cli import (
    CC_RANGE,
    CORES_RANGE,
    T_RANGE,
    USAGE,
    parse_args,
    resolve_output_paths,
)
from ..errors import LatticeGenError
from ..sysinfo import logical_core_count
from . import model
from .model import RunState
from .runner import GRACE_SECONDS, LaunchError, Run, build_argv

POLL_MS = 100
BAR_HEIGHT = 16
WIDTH = 460
#: Height of the log pane, in text lines. Deep enough to read a stage's worth
#: of output at a glance and shallow enough that the window still sits in a
#: corner of the screen, which is what §3.1 asks of it.
LOG_ROWS = 8

WHAT_IT_MAKES = """\
latticegen2 fills the solid body of a STEP file with a diamond-strut lattice \
and writes the result as a single AP214 STEP, in the same coordinate system and \
trimmed exactly to the input's surfaces. Struts stand on their tip, three to a \
node, reclined about 55° from vertical.

The output sits beside a .log recording the parameters, the per-stage timings, \
the run's characteristics and its peak memory — written every run, whether or \
not this window is open.
"""


class Bar(tk.Canvas):
    """A progress bar with its own text drawn over it.

    A canvas rather than :class:`ttk.Progressbar` because ttk has no way to put
    a label inside a bar — the alternatives are a ``place``d widget on top or a
    custom style layout, and both fight whatever theme the platform supplies.
    Thirty lines of canvas look the same on Windows and Linux and make the
    indeterminate sweep trivial.
    """

    def __init__(self, master, **kw):
        super().__init__(master, height=BAR_HEIGHT, highlightthickness=0,
                         bd=0, **kw)
        self._fraction = 0.0
        self._text = ""
        self._sweep = None          # None when determinate
        #: Identifies the current sweep, so a stale rescheduled callback from a
        #: previous one retires instead of doubling the animation's speed. A
        #: pending `after` cannot be un-scheduled by setting `_sweep` to None:
        #: it still fires, and if the bar has gone indeterminate again by then
        #: it would see a live sweep and reschedule itself alongside the new
        #: loop.
        self._sweep_id = 0
        self.bind("<Configure>", lambda _e: self._redraw())

    def set(self, fraction: float | None, text: str = "") -> None:
        self._text = text
        if fraction is None:
            if self._sweep is None:
                self._sweep = 0.0
                self._sweep_id += 1
                self._animate(self._sweep_id)
        else:
            self._sweep = None
            self._fraction = min(max(fraction, 0.0), 1.0)
        self._redraw()

    def clear(self) -> None:
        """Empty the bar and stop any sweep.

        What a finished run leaves behind. An indeterminate bar means "something
        is happening and it cannot be counted", so continuing to sweep once the
        run has ended says the opposite of the truth — and it is the one part of
        the window still moving, which reads as a run that has not stopped.
        """
        self._sweep = None
        self._sweep_id += 1
        self._fraction = 0.0
        self._text = ""
        self._redraw()

    def _animate(self, sweep_id: int) -> None:
        if self._sweep is None or sweep_id != self._sweep_id:
            return
        self._sweep = (self._sweep + 0.02) % 1.0
        self._redraw()
        self.after(40, self._animate, sweep_id)

    def _redraw(self) -> None:
        self.delete("all")
        width = max(self.winfo_width(), 1)
        height = BAR_HEIGHT
        self.create_rectangle(0, 0, width, height, fill="#e6e6e6", outline="#bdbdbd")
        if self._sweep is None:
            self.create_rectangle(0, 0, int(width * self._fraction), height,
                                  fill="#4a90d9", outline="")
        else:
            # A block travelling left to right: honest about not knowing the
            # fraction, while still showing the run is alive.
            span = max(int(width * 0.25), 24)
            left = int((width + span) * self._sweep) - span
            self.create_rectangle(max(left, 0), 0, min(left + span, width), height,
                                  fill="#4a90d9", outline="")
        if self._text:
            self.create_text(width // 2, height // 2, text=self._text,
                             font=("TkDefaultFont", 8))


class App:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.run: Run | None = None
        self.state = RunState()
        #: Wall-clock deadline for force-stopping, and the run's own start.
        #: **Wall clock, not the event stream.** Both were derived from
        #: `state.elapsed` at first, which only advances when the child emits
        #: something — so the grace period could not expire during a long silent
        #: kernel call, which is precisely the case it exists for, and the
        #: elapsed display froze for the whole of `simplify`.
        self._force_at: float | None = None
        self._started: float | None = None
        #: How many of the child's two streams have reached end of file.
        self._closed_seen = 0

        root.title(f"latticegen2 {__version__}")
        root.resizable(False, False)

        self.input_var = tk.StringVar()
        self.outdir_var = tk.StringVar()
        self.cc_var = tk.StringVar(value="10")
        self.t_var = tk.StringVar(value="1.5")
        self.cores_var = tk.StringVar(value=str(logical_core_count()))
        self.verbose_var = tk.BooleanVar(value=False)
        self.derived_var = tk.StringVar(value="—")
        self.message_var = tk.StringVar(value="")
        #: What the log pane was last drawn from — the number of lines that have
        #: arrived and the filter in force. Redrawing only when that changes
        #: keeps the 10 Hz poll from rebuilding a 400-line widget for nothing.
        self._log_key: tuple | None = None

        self._build(root)
        for var in (self.input_var, self.cc_var, self.t_var, self.cores_var,
                    self.outdir_var):
            var.trace_add("write", lambda *_a: self._revalidate())
        self.verbose_var.trace_add("write", lambda *_a: self._render_log())
        self._revalidate()
        root.protocol("WM_DELETE_WINDOW", self._on_exit)
        root.after(POLL_MS, self._drain)

    # -- layout ------------------------------------------------------------

    def _build(self, root: tk.Tk) -> None:
        outer = ttk.Frame(root, padding=8)
        outer.grid(sticky="nsew")
        outer.columnconfigure(1, weight=1)

        header = ttk.Frame(outer)
        header.grid(row=0, column=0, columnspan=3, sticky="ew")
        header.columnconfigure(0, weight=1)
        ttk.Label(header, text="Generate a lattice into a STEP volume").grid(
            row=0, column=0, sticky="w")
        ttk.Button(header, text="?", width=3, command=self._show_help).grid(
            row=0, column=1, sticky="e")

        self.fields: list[ttk.Widget] = []

        ttk.Label(outer, text="Input").grid(row=1, column=0, sticky="w", pady=(8, 2))
        entry = ttk.Entry(outer, textvariable=self.input_var, width=44)
        entry.grid(row=1, column=1, sticky="ew", pady=(8, 2))
        browse = ttk.Button(outer, text="Browse…", command=self._pick_input)
        browse.grid(row=1, column=2, sticky="e", padx=(4, 0), pady=(8, 2))
        self.fields += [entry, browse]

        ttk.Label(outer, text="Output").grid(row=2, column=0, sticky="w", pady=2)
        outdir = ttk.Entry(outer, textvariable=self.outdir_var, width=44)
        outdir.grid(row=2, column=1, sticky="ew", pady=2)
        outbrowse = ttk.Button(outer, text="Browse…", command=self._pick_outdir)
        outbrowse.grid(row=2, column=2, sticky="e", padx=(4, 0), pady=2)
        self.fields += [outdir, outbrowse]

        ttk.Label(outer, textvariable=self.derived_var, foreground="#555").grid(
            row=3, column=1, columnspan=2, sticky="w")

        params = ttk.Frame(outer)
        params.grid(row=4, column=0, columnspan=3, sticky="ew", pady=(6, 2))
        self.fields += [
            self._spin(params, "cc", self.cc_var, CC_RANGE, 0.1, 0),
            self._spin(params, "t", self.t_var, T_RANGE, 0.1, 2),
            self._spin(params, "cores", self.cores_var, CORES_RANGE, 1, 4),
        ]
        # **Deliberately not in `self.fields`.** Everything else is a parameter
        # of the run and is frozen once it starts; this is a filter over output
        # that has already arrived. Every line crosses the event stream whichever
        # way the box is set (docs/algorithm.md §10), so it can be changed
        # *during* a run — or after a failure, to read what led to it — and the
        # pane redraws from what was already received. That is the whole reason
        # the child is never given `-v` itself.
        self.verbose_check = ttk.Checkbutton(params, text="verbose",
                                             variable=self.verbose_var)
        self.verbose_check.grid(row=0, column=6, sticky="w", padx=(12, 0))

        self.message = ttk.Label(outer, textvariable=self.message_var,
                                 foreground="#a11", wraplength=WIDTH - 24)
        self.message.grid(row=5, column=0, columnspan=3, sticky="w", pady=(2, 4))

        buttons = ttk.Frame(outer)
        buttons.grid(row=6, column=0, columnspan=3, sticky="ew")
        buttons.columnconfigure(0, weight=1)
        buttons.columnconfigure(1, weight=1)
        self.start = ttk.Button(buttons, text="Start!", command=self._on_start)
        self.start.grid(row=0, column=0, sticky="ew", padx=(0, 4))
        ttk.Button(buttons, text="Exit", command=self._on_exit).grid(
            row=0, column=1, sticky="ew", padx=(4, 0))

        # --- the panel that appears while a run is in flight ---------------
        self.panel = ttk.Frame(outer, padding=(0, 8, 0, 0))
        self.panel.columnconfigure(0, weight=1)

        self.stage_label = ttk.Label(self.panel, text="")
        self.stage_label.grid(row=0, column=0, sticky="w")
        self.stage_bar = Bar(self.panel)
        self.stage_bar.grid(row=1, column=0, sticky="ew", pady=(1, 3))
        self.sub_bar = Bar(self.panel)
        self.sub_bar.grid(row=2, column=0, sticky="ew")
        self.resource_label = ttk.Label(self.panel, text="", foreground="#555")
        self.resource_label.grid(row=3, column=0, sticky="w", pady=(3, 0))

        # Created but not gridded. `_render_log` shows it only when it has
        # something in it, so an ordinary run does not carry an empty box
        # around for the hour it takes — and the box appearing is itself the
        # signal that the child said something unexpected.
        self.log_frame = ttk.Frame(self.panel)
        log = self.log_frame
        log.columnconfigure(0, weight=1)
        # `width=1` rather than a character count, so the pane takes the width
        # the rest of the window already has instead of setting it: the root is
        # not resizable, so a Text asking for 60 columns would widen everything.
        self.log_text = tk.Text(log, height=LOG_ROWS, width=1, wrap="word",
                                font=("TkFixedFont", 8), state="disabled",
                                background="#fbfbfb", relief="solid",
                                borderwidth=1, highlightthickness=0)
        self.log_text.grid(row=0, column=0, sticky="ew")
        scroll = ttk.Scrollbar(log, orient="vertical",
                               command=self.log_text.yview)
        scroll.grid(row=0, column=1, sticky="ns")
        self.log_text.configure(yscrollcommand=scroll.set)

        self.result_label = ttk.Label(self.panel, text="", wraplength=WIDTH - 24)
        self.result_label.grid(row=5, column=0, sticky="w", pady=(6, 0))
        self.open_button = ttk.Button(self.panel, text="Open export folder",
                                      command=self._open_output_folder)

    def _spin(self, parent, label, var, bounds, step, column):
        ttk.Label(parent, text=label).grid(row=0, column=column, sticky="w",
                                           padx=(0 if column == 0 else 10, 4))
        spin = ttk.Spinbox(parent, textvariable=var, width=8, increment=step,
                           from_=bounds[0], to=bounds[1])
        spin.grid(row=0, column=column + 1, sticky="w")
        return spin

    # -- input handling ----------------------------------------------------

    def _pick_input(self) -> None:
        path = filedialog.askopenfilename(
            title="Input STEP file",
            filetypes=[("STEP files", "*.step *.stp *.STEP *.STP"), ("All files", "*.*")],
        )
        if path:
            self.input_var.set(os.path.abspath(path))
            if not self.outdir_var.get().strip():
                self.outdir_var.set(os.path.dirname(os.path.abspath(path)))

    def _pick_outdir(self) -> None:
        current = self.outdir_var.get().strip() or os.path.dirname(
            self.input_var.get().strip())
        path = filedialog.askdirectory(title="Output folder",
                                       initialdir=current or os.getcwd())
        if path:
            self.outdir_var.set(os.path.abspath(path))

    def _pending_argv(self) -> list[str]:
        """The command line the current field values describe."""
        return [
            "-i", self.input_var.get().strip(),
            "-cc", self.cc_var.get().strip(),
            "-t", self.t_var.get().strip(),
            "--cores", self.cores_var.get().strip(),
        ]

    def _revalidate(self) -> None:
        """Re-check the fields with the *parser*, not with a copy of its rules.

        Every bound, the ``t < cc/sqrt(2)`` cross-constraint and the derived
        output name all come from :mod:`latticegen2.cli`. Reimplementing any of
        them here would mean the window eventually accepting something the run
        would refuse, and the user meeting the difference as a failed run rather
        than as a disabled button.
        """
        if self.run is not None:
            return
        self.message_var.set("")
        try:
            args = parse_args(self._pending_argv())
        except LatticeGenError as exc:
            self.derived_var.set("—")
            first = str(exc).split("\n", 1)[0]
            self.message_var.set(first if self.input_var.get().strip() else "")
            self.start.state(["disabled"])
            return

        name = os.path.basename(args.output)
        self.derived_var.set(f"→ {name}")
        self.start.state(["!disabled"])

    def _resolved_output(self) -> str:
        """Where the .step goes: the chosen folder, with the derived filename.

        The name is never editable, and comes from
        :func:`latticegen2.cli.resolve_output_paths` rather than being formatted
        here — ``format_param`` renders ``20.0`` as ``20`` and ``1.5`` as
        ``1.5``, and a second implementation of that would drift.
        """
        args = parse_args(self._pending_argv())
        default_step, _log = resolve_output_paths(
            args.input, None, args.cc, args.t)
        folder = self.outdir_var.get().strip() or os.path.dirname(args.input)
        return os.path.join(os.path.abspath(folder), os.path.basename(default_step))

    # -- running -----------------------------------------------------------

    def _on_start(self) -> None:
        if self.run is not None:
            self._on_stop()
            return
        try:
            args = parse_args(self._pending_argv())
            output = self._resolved_output()
            # The same filesystem checks the CLI runs, run here so an unreadable
            # input or an unwritable folder is a message in the window rather
            # than a child that exits before it has drawn anything.
            from ..cli import preflight_checks

            preflight_checks(parse_args(self._pending_argv() + ["-o", output]))
        except LatticeGenError as exc:
            self.message_var.set(str(exc).split("\n", 1)[0])
            return

        self.state = RunState()
        try:
            self.run = Run(
                build_argv(
                    input_path=args.input, output_path=output,
                    cc=self.cc_var.get().strip(), t=self.t_var.get().strip(),
                    cores=int(self.cores_var.get().strip()),
                ),
                log_path=os.path.splitext(output)[0] + ".log",
            )
        except (LaunchError, OSError) as exc:
            # Anything that stops the child existing at all. Without this the
            # exception escapes into the Tk callback, where a windowed process
            # has nowhere to print it: the button would simply do nothing.
            self.message_var.set(str(exc))
            return
        self._force_at = None
        self._started = time.monotonic()
        self._closed_seen = 0
        self._log_key = None
        self.start.configure(text="Stop!")
        for widget in self.fields:
            widget.state(["disabled"])
        self.result_label.configure(text="")
        self.open_button.grid_forget()
        self.panel.grid(row=7, column=0, columnspan=3, sticky="ew")
        self._refresh()

    def _on_stop(self) -> None:
        if self.run is None or not self.run.running:
            return
        self.run.cancel()
        self._force_at = time.monotonic() + GRACE_SECONDS
        self.start.state(["disabled"])
        self.start.configure(text="Cancelling…")
        # Said plainly rather than hidden: the interrupt is checked between
        # bytecodes, so inside a long kernel call it lands when that call
        # returns. Exactly what Ctrl+C does on the command line.
        self.result_label.configure(
            text="Cancelling — this can take a moment if a kernel call is running.")

    def _on_exit(self) -> None:
        if self.run is not None and self.run.running:
            if not messagebox.askokcancel(
                    "latticegen2", "A run is in progress. Stop it and exit?"):
                return
            self.run.cancel()
            self.run.force_stop()
        self.root.destroy()

    # -- the event loop ----------------------------------------------------

    def _drain(self) -> None:
        """Move everything the reader threads collected into the widgets.

        The only place events reach Tk, and it runs on the main thread by
        construction because ``after`` schedules it there.
        """
        run = self.run
        if run is not None:
            drained = False
            while True:
                try:
                    kind, payload = run.events.get_nowait()
                except Exception:
                    break
                drained = True
                if kind == "event":
                    self.state.apply(payload)
                elif kind in ("raw", "stderr"):
                    self.state.add_line(payload, model.OUT)
                elif kind == "closed":
                    self._closed_seen += 1
            if drained or run.running:
                self._refresh()
            if self._force_at is not None and time.monotonic() > self._force_at \
                    and run.running:
                run.force_stop()
                self._force_at = None
            # **Both streams must have reached EOF, not merely "the queue looks
            # empty".** The child's exit code is available the instant it dies,
            # while the pipes can still hold buffered lines and a reader thread
            # can be between its read and its `put` — so finishing on an empty
            # queue alone can miss the terminal `result` event and report a
            # perfectly good run as a failure.
            if not run.running and self._closed_seen >= 2 and run.events.empty():
                self._finish(run)
        self.root.after(POLL_MS, self._drain)

    def _refresh(self) -> None:
        state = self.state
        if self.run is not None and self._started is not None:
            # Advanced here rather than only from events, so the clock keeps
            # moving through a stage that reports nothing for nineteen minutes.
            # `max` so the summary's own authoritative duration still wins at
            # the end.
            state.elapsed = max(state.elapsed, time.monotonic() - self._started)
        self.stage_label.configure(text=state.stage_text())
        self.stage_bar.set(state.permille / 1000.0, f"{state.permille // 10}%")
        if state.finished:
            # Success, cancellation and failure alike: there is no sub-stage
            # left to report, and `sub_fraction` is `None` at this point for the
            # same reason it is during `export` — which would otherwise leave the
            # bar sweeping for as long as the window stayed open.
            self.sub_bar.clear()
        else:
            self.sub_bar.set(state.sub_fraction, state.sub_text())
        self.resource_label.configure(text=state.resource_text())
        self._render_log()

    def _render_log(self) -> None:
        """Draw the lines the tick box admits, and hide the pane when there
        are none.

        Keyed on ``lines_seen`` rather than ``len(lines)`` because the backlog
        is a bounded deque: once it is full its length stops changing, and a
        length-keyed redraw would freeze the pane at exactly the point in a long
        run where there is most to read.
        """
        verbose = bool(self.verbose_var.get())
        key = (self.state.lines_seen, verbose)
        if key == self._log_key:
            return
        self._log_key = key
        visible = self.state.visible_lines(verbose)
        if visible:
            self.log_frame.grid(row=4, column=0, sticky="ew", pady=(6, 0))
        else:
            self.log_frame.grid_forget()
        self.log_text.configure(state="normal")
        self.log_text.delete("1.0", "end")
        if visible:
            self.log_text.insert("1.0", "\n".join(visible))
        self.log_text.configure(state="disabled")
        self.log_text.see("end")

    def _finish(self, run: Run) -> None:
        code = run.returncode()
        state = self.state
        if state.status is None:
            # The child died without reporting itself — killed, or a crash
            # below Python. Its stderr is the only account of it there is.
            tail = " ".join(run.stderr_lines[-2:])[:300]
            state.status = progress.CANCELLED if code == 130 else progress.FAILED
            state.reason = tail or f"the run ended with exit code {code}"
        run.clear_cancel()
        self.run = None
        self._force_at = None
        self._started = None
        self.start.state(["!disabled"])
        self.start.configure(text="Start!")
        for widget in self.fields:
            widget.state(["!disabled"])
        self._refresh()

        text = state.result_text()
        if state.tmpdir:
            text += f"\nIntermediate files kept in: {state.tmpdir}"
        self.result_label.configure(
            text=("✔ " if state.status == progress.OK else "✘ ") + text)
        self.open_button.grid(row=6, column=0, sticky="w", pady=(4, 0))
        self._revalidate()

    # -- odds and ends -----------------------------------------------------

    def _open_output_folder(self) -> None:
        folder = os.path.dirname(self.state.output or "") or \
            self.outdir_var.get().strip()
        if not folder or not os.path.isdir(folder):
            return
        try:
            if sys.platform == "win32":
                os.startfile(folder)  # noqa: S606 - a directory this app chose
            elif sys.platform == "darwin":
                subprocess.Popen(["open", folder])
            else:
                subprocess.Popen(["xdg-open", folder])
        except Exception as exc:  # noqa: BLE001
            messagebox.showwarning("latticegen2", f"Could not open {folder}\n\n{exc}")

    def _show_help(self) -> None:
        top = tk.Toplevel(self.root)
        top.title("latticegen2 — help")
        top.transient(self.root)
        frame = ttk.Frame(top, padding=10)
        frame.grid(sticky="nsew")
        ttk.Label(frame, text="What this produces", font=("TkDefaultFont", 9, "bold")
                  ).grid(sticky="w")
        ttk.Label(frame, text=WHAT_IT_MAKES, wraplength=520, justify="left").grid(
            sticky="w", pady=(2, 8))
        ttk.Label(frame, text="Parameters", font=("TkDefaultFont", 9, "bold")).grid(
            sticky="w")
        # The command line's own usage text, verbatim. One description of the
        # parameters, so the two can never disagree about a range.
        text = tk.Text(frame, width=76, height=22, wrap="none")
        text.insert("1.0", USAGE)
        text.configure(state="disabled")
        text.grid(sticky="nsew", pady=(2, 0))
        ttk.Button(frame, text="Close", command=top.destroy).grid(
            sticky="e", pady=(8, 0))


def launch() -> int:
    root = tk.Tk()
    App(root)
    root.mainloop()
    return 0
