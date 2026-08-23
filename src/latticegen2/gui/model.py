"""What the window shows, derived from the event stream. No Tk in this module.

The front-end is two halves: this one decides *what* to display, and
:mod:`latticegen2.gui.app` decides how to draw it. They are separated because
only one of them can be tested — a Tk widget's state is awkward to assert and
needs a display, while a reduction from a list of events to a dataclass needs
neither.

Everything here tolerates a stream that is incomplete or out of order. It has
to: the producer is a separate process that can be killed at any point, so a
missing ``stage_end``, a truncated final line, or an interpreter warning printed
before the first event are all ordinary rather than exceptional.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field

from .. import progress
from ..runlog import STAGES, format_bytes, format_duration
from .weights import overall_permille

#: How many lines of the run's own output to keep for the failure panel. The
#: rehearsal writes a few hundred; a bounded deque means an hour-long run cannot
#: grow the window's memory without bound.
LOG_BACKLOG = 400


@dataclass
class RunState:
    """Everything the window draws, and nothing about how it is drawn."""

    stage: str | None = None
    stage_number: int = 0
    #: Sub-stage label, e.g. ``"boundary trim"`` or ``"writing STEP"``.
    label: str = ""
    done: int | None = None
    total: int | None = None
    #: Top bar, 0-1000. Monotone by construction — see :meth:`apply`.
    permille: int = 0
    peak_rss: int = 0
    elapsed: float = 0.0
    #: ``None`` while the run is in flight, then one of `progress.STATUSES`.
    status: str | None = None
    reason: str | None = None
    output: str | None = None
    log: str | None = None
    tmpdir: str | None = None
    summary: dict = field(default_factory=dict)
    lines: deque = field(default_factory=lambda: deque(maxlen=LOG_BACKLOG))
    #: How many lines have *ever* arrived, which ``len(lines)`` cannot say once
    #: the backlog is full and the deque starts dropping from the left. The
    #: window redraws its log pane only when this changes, so without it a long
    #: run would stop updating the pane at exactly the point there is most to
    #: read.
    lines_seen: int = 0

    def add_line(self, text: str, console: bool) -> None:
        """Record one line of the run's own output.

        ``console`` is whether a command-line run would have printed it, which
        is what the window's **verbose** tick box filters on (specification.md
        §3.1). Every line crosses regardless, so the filter can be changed
        while the run is going.
        """
        self.lines.append((text, console))
        self.lines_seen += 1

    # -- derived, for the widgets -----------------------------------------

    @property
    def indeterminate(self) -> bool:
        """True when the sub-bar has no fraction to show and should sweep."""
        return self.total is None

    @property
    def sub_fraction(self) -> float | None:
        if self.total is None or self.total <= 0:
            return None
        return min(max((self.done or 0) / self.total, 0.0), 1.0)

    @property
    def finished(self) -> bool:
        return self.status is not None

    def stage_text(self) -> str:
        if self.stage is None:
            return "starting…" if not self.finished else ""
        return f"{self.stage}  ({self.stage_number} of {len(STAGES)})"

    def sub_text(self) -> str:
        """The line drawn over the sub-bar."""
        if not self.label:
            return ""
        if self.total is None:
            return self.label
        return f"{self.label}: {self.done:,} / {self.total:,}"

    def resource_text(self) -> str:
        """The one instrumented figure this tool actually gathers, plus a clock.

        Peak memory is genuinely all there is: `RunLog` tracks ``max_rss`` and
        folds each worker's own reading into it (specification.md §3), and
        nothing else is measured continuously. Showing the two together is
        honest about that rather than padding the line out.
        """
        rss = format_bytes(self.peak_rss) if self.peak_rss else "—"
        return f"peak memory {rss}   ·   elapsed {format_duration(self.elapsed)}"

    def result_text(self) -> str:
        if self.status == progress.OK:
            solids = self.summary.get("solids_written")
            size = self.summary.get("output_size")
            detail = ", ".join(str(x) for x in (
                f"{solids} solid(s)" if solids else None, size) if x)
            return f"Completed in {format_duration(self.elapsed)}" + (
                f" — {detail}" if detail else "")
        if self.status == progress.CANCELLED:
            return f"Cancelled after {format_duration(self.elapsed)}"
        if self.status == progress.FAILED:
            return f"Failed: {self.reason or 'see the log'}"
        return ""

    # -- the reduction -----------------------------------------------------

    def apply(self, event: dict) -> None:
        """Fold one event in. Unknown events are ignored, never fatal."""
        ev = event.get("ev")
        self.elapsed = max(self.elapsed, float(event.get("t") or 0.0))

        if ev == progress.HELLO:
            self.output = event["output"]
            self.log = event["log"]

        elif ev == progress.STAGE_BEGIN:
            self.stage = event["name"]
            self.stage_number = event["i"] + 1
            self.label = ""
            self.done = self.total = None
            self._advance(overall_permille(self.stage, 0.0))

        elif ev == progress.STAGE_END:
            self.peak_rss = max(self.peak_rss, int(event.get("max_rss") or 0))
            # To the *end* of the finished stage's band, which is the start of
            # the next one's. Done here rather than on the next `stage_begin`
            # so the bar reflects completed work even if the run stops between
            # two stages.
            self._advance(overall_permille(event["name"], 1.0))
            self.label = ""
            self.done = self.total = None

        elif ev == progress.PROGRESS:
            # `stage` from the event rather than from our own state: a tick that
            # arrives after its stage ended still belongs to that stage.
            self.label = event["label"] or ""
            self.done = event["done"]
            self.total = event["total"]
            self._advance(overall_permille(event["stage"] or self.stage,
                                           self.sub_fraction))

        elif ev == progress.LOG:
            self.add_line(event["msg"], bool(event["console"]))

        elif ev == progress.SUMMARY:
            self.summary = dict(event["stats"])
            self.peak_rss = max(self.peak_rss, int(event.get("peak_rss") or 0))
            self.elapsed = max(self.elapsed, float(event.get("duration") or 0.0))

        elif ev == progress.RESULT:
            self.status = event["status"]
            self.reason = event["reason"]
            self.output = event["output"] or self.output
            self.log = event["log"] or self.log
            self.tmpdir = event["tmpdir"]
            self.stage = None
            self.label = ""
            self.done = self.total = None
            if self.status == progress.OK:
                self.permille = 1000

    def _advance(self, value: int) -> None:
        """Move the bar forward only.

        Ordered ``imap`` means a completion count never goes backwards, but the
        *bar* has a second source — ``stage_end`` fills its band completely —
        and a late tick from the stage that just ended would otherwise pull it
        back. A progress bar that retreats reads as something having gone wrong,
        which here would be false.
        """
        self.permille = max(self.permille, min(int(value), 1000))


def reduce_stream(lines) -> RunState:
    """Fold a whole stream of stdout lines into one state. Used by the tests."""
    state = RunState()
    for line in lines:
        event = progress.parse_line(line)
        if event is not None:
            state.apply(event)
        elif line.strip():
            state.add_line(line.rstrip("\n"), True)
    return state
