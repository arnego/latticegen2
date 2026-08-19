"""Machine-readable run events: the NDJSON stream the graphical front-end reads.

A run reports itself two ways. The ``.log`` file and the console text are for a
person and are documented in docs/algorithm.md §10; this module is the second
channel, for :mod:`latticegen2.gui`, and it exists because a bar cannot be driven
by prose. It is enabled only by ``--progress-stream`` and is inert otherwise, so
a plain CLI run pays nothing and produces exactly the bytes it always did.

**Nothing here may change what a run computes or writes.** The emitter is handed
to :class:`latticegen2.runlog.RunLog` and observes it; the ``.log`` file is
written by the same code either way, and ``test/test_runlog_events.py`` pins that
as a regression rather than leaving it as an intention.

The schema lives here rather than in either of the two modules that use it, so
the producer (`runlog`, `pipeline`, `__main__`) and the consumer (`gui.runner`)
cannot drift apart: :data:`FIELDS` is the single declaration of what every event
carries, and :class:`NdjsonEmitter` rejects an event whose fields do not match it
exactly. That is deliberately strict — the alternative is a consumer full of
``.get(...)`` calls, each of which silently tolerates a field that stopped being
emitted.

**Stdlib only, and it must stay that way.** The GUI process imports this module
and must never load OCP; see docs/algorithm.md §13.
"""

from __future__ import annotations

import json
import sys
import time

SCHEMA_VERSION = 1

#: The sequence a run announces once, before any stage.
HELLO = "hello"
#: A stage started. There is no counterpart in the ``.log`` — see `RunLog.stage_begin`.
STAGE_BEGIN = "stage_begin"
#: A stage finished, carrying the same elapsed time the ``.log`` records.
STAGE_END = "stage_end"
#: Sub-stage work. ``total`` is ``None`` where the stage has no countable unit.
PROGRESS = "progress"
#: One line of what a CLI run would have written, with whether it would have shown.
LOG = "log"
#: The end-of-run summary's contents, structured rather than formatted.
SUMMARY = "summary"
#: How the run ended. Always the last event.
RESULT = "result"

#: Every event's fields, excluding ``ev`` and ``t`` which all of them carry.
#:
#: An emitted event must supply exactly these keys — no more, and none missing.
#: ``None`` is how "not known" is expressed, so the consumer can index rather
#: than probe.
FIELDS: dict[str, tuple[str, ...]] = {
    HELLO: ("v", "pid", "stages", "output", "log", "workers"),
    STAGE_BEGIN: ("name", "i", "n"),
    STAGE_END: ("name", "i", "n", "elapsed", "max_rss"),
    PROGRESS: ("stage", "label", "done", "total"),
    LOG: ("msg", "console"),
    SUMMARY: ("duration", "peak_rss", "params", "stats", "stages"),
    RESULT: ("status", "exit", "reason", "output", "log", "tmpdir"),
}

EVENTS: tuple[str, ...] = tuple(FIELDS)

#: Suffix of the file whose appearance means "stop this run".
CANCEL_SUFFIX = ".cancel"


def cancel_sentinel_path(log_path: str) -> str:
    """Where the front-end asks a run to stop, given that run's ``.log`` path.

    **A file, because neither of the two obvious channels works here.**

    Signals do not: the window runs under ``pythonw.exe`` and therefore has no
    console, and ``GenerateConsoleCtrlEvent`` can only signal process groups on
    the *caller's* console — a Windows GUI has no way to send a child Ctrl+Break
    at all.

    A line on the child's stdin does not either, and that one had to be measured
    rather than reasoned about. Giving the child an anonymous pipe for stdin
    stalls the worker pool the first time it spawns: on this machine
    ``classify`` ran for 0.2 s with ``stdin`` on the null device and hung
    indefinitely with a pipe, and it hung whether or not the parent kept its end
    open and whether or not ``CREATE_NO_WINDOW`` was set. So it is the child's
    file descriptor 0 *being a pipe* that ``spawn`` cannot live with, not
    anything about how the parent uses it — which rules the whole channel out
    rather than suggesting a way to hold it correctly.

    A file has neither problem: no inherited handle, no signal, and it works the
    same on both platforms. It is derived from the log path rather than passed as
    a flag, so there is nothing extra on the command line and both sides compute
    it from something they already agree on.
    """
    import os

    return os.path.splitext(log_path)[0] + CANCEL_SUFFIX


#: ``result.status`` values.
OK = "ok"
FAILED = "failed"
CANCELLED = "cancelled"
STATUSES = (OK, FAILED, CANCELLED)


class NdjsonEmitter:
    """Write one JSON object per line, to a stream the GUI is reading.

    ``ensure_ascii`` is on and non-negotiable: the stream crosses a pipe whose
    encoding is the child's stdout encoding, which on Windows is whatever the
    console code page happens to be. An input path containing a non-ASCII
    character would otherwise be mojibake or a ``UnicodeEncodeError`` in the
    middle of a run — escaping it to an ASCII unicode escape keeps the line pure
    ASCII and the decoding lossless.

    Every line is flushed. The consumer is a progress bar, so a buffered event
    is a bar that jumps; and the volume is bounded by design (docs/algorithm.md
    §12 — events are emitted at stage and batch boundaries, never inside an
    inner loop).

    **A dead consumer must not kill the run.** If the GUI vanishes, its end of
    the pipe closes and the next write raises. That is not a reason to abandon
    an hour of geometry: the emitter goes inert and the pipeline finishes and
    writes its STEP exactly as a CLI run would. (The child does learn about it,
    separately and deliberately, through stdin reaching EOF — which *is* treated
    as a cancellation, because it is one. See :mod:`latticegen2.__main__`.)
    """

    def __init__(self, stream=None, t0: float | None = None, clock=time.perf_counter):
        self._stream = sys.stdout if stream is None else stream
        self._clock = clock
        self._t0 = clock() if t0 is None else t0
        self.dead = False

    def __call__(self, ev: str, **fields) -> None:
        expected = FIELDS.get(ev)
        if expected is None:
            raise ValueError(f"unknown progress event {ev!r}")
        if set(fields) != set(expected):
            missing = sorted(set(expected) - set(fields))
            extra = sorted(set(fields) - set(expected))
            raise ValueError(
                f"progress event {ev!r} field mismatch: missing {missing}, extra {extra}"
            )
        if self.dead:
            return
        payload = {"ev": ev, "t": round(self._clock() - self._t0, 3)}
        payload.update(fields)
        try:
            self._stream.write(
                json.dumps(payload, ensure_ascii=True, separators=(",", ":")) + "\n"
            )
            self._stream.flush()
        except Exception:
            # Broken pipe, closed stream, or an encoding the stream cannot take.
            # See the class docstring: the run continues, silently.
            self.dead = True


def parse_line(text: str) -> dict | None:
    """Decode one stdout line, or ``None`` if it is not an event.

    Returning ``None`` rather than raising is the whole contract. Three things
    legitimately arrive on this stream that are not events: OCCT writes to file
    descriptor 1 directly, below ``sys.stdout``, so kernel chatter can appear
    mid-stream; a run killed during a write leaves a truncated final line; and
    an interpreter can emit a warning before anything else. None of those is a
    reason to stop reading, so the consumer shows them as raw text and carries
    on.
    """
    text = text.strip()
    if not text or not text.startswith("{"):
        return None
    try:
        obj = json.loads(text)
    except ValueError:
        return None
    if not isinstance(obj, dict) or obj.get("ev") not in FIELDS:
        return None
    return obj
