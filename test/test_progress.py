"""The NDJSON progress schema (docs/algorithm.md §10).

Two properties matter here and neither is obvious from reading the module.

The schema is a *contract between two processes*, so the producer and the
consumer cannot be checked against each other by any type system — only by
agreeing on :data:`latticegen2.progress.FIELDS`. These tests pin that every
event supplies exactly the keys it declares, which is what lets the front-end
index rather than probe.

And the reader has to survive a stream that is not entirely events. OCCT writes
to file descriptor 1 below ``sys.stdout``, and a killed process leaves a
truncated final line — so ``parse_line`` returning ``None`` on junk is a
requirement, not a convenience.
"""

import io
import json

import pytest

from latticegen2 import progress
from latticegen2.runlog import STAGES


def emitter():
    stream = io.StringIO()
    return stream, progress.NdjsonEmitter(stream, t0=0.0, clock=lambda: 0.5)


def test_stage_names_are_the_twelve_pipeline_stages_in_order():
    # Pinned as a literal rather than derived, because the whole point of
    # runlog.STAGES is to be the one declaration everything else agrees with.
    assert STAGES == (
        "template", "import", "tessellate", "classify", "boundary", "connect",
        "stitch", "instance", "assemble", "simplify", "validate", "export",
    )


def test_every_event_writes_exactly_one_json_line():
    stream, emit = emitter()
    emit(progress.STAGE_BEGIN, name="classify", i=3, n=12)
    emit(progress.PROGRESS, stage="classify", label="slices", done=1, total=24)
    text = stream.getvalue()
    assert text.count("\n") == 2
    assert text.endswith("\n")
    for line in text.splitlines():
        json.loads(line)


def test_each_event_carries_its_full_declared_key_set():
    for ev, fields in progress.FIELDS.items():
        stream, emit = emitter()
        emit(ev, **{k: None for k in fields})
        payload = json.loads(stream.getvalue())
        assert set(payload) == {"ev", "t"} | set(fields), ev
        assert payload["ev"] == ev


@pytest.mark.parametrize("supplied", [{}, {"name": "classify"}, {"name": "x", "i": 1, "n": 12, "extra": 1}])
def test_a_field_mismatch_raises_rather_than_emitting_a_partial_event(supplied):
    stream, emit = emitter()
    with pytest.raises(ValueError):
        emit(progress.STAGE_BEGIN, **supplied)
    assert stream.getvalue() == ""


def test_non_ascii_paths_survive_the_stream():
    # The pipe's encoding is the child's stdout encoding, which on Windows is a
    # console code page. ensure_ascii is what stops a Norwegian path being
    # mojibake or a UnicodeEncodeError in the middle of an hour-long run.
    stream, emit = emitter()
    path = "C:\\Bruker\\Ø-deler\\køle-rør.step"
    emit(progress.RESULT, status=progress.OK, exit=0, reason=None,
         output=path, log=None, tmpdir=None)
    raw = stream.getvalue()
    assert raw.isascii()
    assert progress.parse_line(raw)["output"] == path


def test_timestamps_do_not_go_backwards():
    stream = io.StringIO()
    ticks = iter([0.0, 1.0, 1.0, 7.5])
    emit = progress.NdjsonEmitter(stream, t0=0.0, clock=lambda: next(ticks))
    for _ in range(3):
        emit(progress.LOG, msg="x", console=False)
    stamps = [json.loads(line)["t"] for line in stream.getvalue().splitlines()]
    assert stamps == sorted(stamps)


def test_parse_line_round_trips_and_rejects_everything_else():
    stream, emit = emitter()
    emit(progress.LOG, msg="boundary trim: 1/2", console=True)
    assert progress.parse_line(stream.getvalue()) == {
        "ev": "log", "t": 0.5, "msg": "boundary trim: 1/2", "console": True,
    }
    for junk in [
        "",
        "   ",
        "Intel oneMKL FATAL ERROR",          # a native library writing to fd 1
        '{"ev":"stage_begin","name":"cla',   # a run killed mid-write
        '{"ev":"not_an_event","x":1}',       # a future schema, read by an old GUI
        '["ev"]',                            # valid JSON, wrong shape
        '{"no_ev":1}',
    ]:
        assert progress.parse_line(junk) is None, junk


def test_a_dead_consumer_does_not_kill_the_run():
    # The GUI vanishing must not abandon an hour of geometry mid-pipeline. The
    # child learns about it through stdin reaching EOF instead, which is a
    # cancellation it can act on cleanly.
    class Closed(io.StringIO):
        def write(self, _):
            raise BrokenPipeError

    emit = progress.NdjsonEmitter(Closed(), t0=0.0)
    emit(progress.LOG, msg="x", console=False)
    assert emit.dead
    emit(progress.LOG, msg="y", console=False)  # still silent, still no raise
