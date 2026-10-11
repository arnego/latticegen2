"""How much of the progress bar each pipeline stage is worth.

The stages differ by three orders of magnitude in cost, so a bar that advanced
one step per stage would sit still for the ten minutes ``boundary`` takes and
then jump. The shares here come from measured runs instead.

**There are two layouts, and which one a run gets is decided by the stages it
announces.** A run with print supports goes through three more stages than one
without — ``overhang``, ``plan`` and ``support`` — and they are a quarter of
it. They also change what the *shared* stages cost: a supported run has
several times the interfaces to repair in ``connect`` and more faces to carry
through ``stitch``, ``simplify`` and ``validate``. So the supported layout is
not the plain one with three bands squeezed in; each is laid out from its own
measurement, :data:`SECONDS_PLAIN` and :data:`SECONDS_SUPPORT`, and both
tables are computed once, at import.

**Provenance.** One controlled pair on the 6-core development workstation, on
2026-10-11: `TD_HX_rehearsal_test` at ``cc=7, t=1.4``, ``--orient 180 0 0``,
six cores, without and with ``--support --overhang 60 --support-thickness 1``.
Same part, same parameters, same commit; the only difference is the flag. The
supported run is the ``rehearsal-support`` e2e scenario (docs/specification.md
§6.1), so re-measuring it is one command.

**What these tables are not.** They are one part's shape, on one machine. A
small part is boundary- and export-dominated — on ``80mm-test-ball``
``simplify`` is trivial — so the bar will visibly jump there, and a part with
little overhanging surface will run through the three support stages almost at
once. ``validate`` is dominated by the export-truth gate (docs/algorithm.md
§9), which writes, re-reads and tessellates every output solid; at
``cc=5, t=1`` on this same part that gate alone is half an hour and
``validate`` is the largest stage of the run, where here it is a sixth of it.
That is the accepted trade: a fixed weighting is monotone and never lies about
what has finished, where a time-based estimate would sit at 100 % of a stage
that had not ended. Read the bar as "how far through the work", never as an
ETA.

**And it rots.** If either run is re-measured, these two tables are downstream
consumers of those numbers; docs/profiling-reports.md says so beside the
entry.
"""

from __future__ import annotations

from ..runlog import STAGES, SUPPORT_STAGES, stages_for

#: Wall-clock seconds per stage of the run **without** ``--support``.
SECONDS_PLAIN: dict[str, float] = {
    "template": 0.05,
    "import": 1.0,
    "tessellate": 2.3,
    "classify": 14.7,
    "boundary": 668.2,
    "connect": 11.0,
    "stitch": 258.7,
    "instance": 13.1,
    "assemble": 8.3,
    "simplify": 234.9,
    "validate": 434.4,
    "export": 90.3,
}  # 29 m 00 s end to end

#: Wall-clock seconds per stage of the same run **with** ``--support``.
SECONDS_SUPPORT: dict[str, float] = {
    "template": 0.05,
    "import": 1.3,
    "tessellate": 2.3,
    "classify": 12.9,
    "boundary": 622.3,
    "overhang": 7.0,
    "plan": 52.1,
    "support": 740.7,
    "connect": 208.3,
    "stitch": 395.6,
    "instance": 11.9,
    "assemble": 25.4,
    "simplify": 319.7,
    "validate": 502.0,
    "export": 108.8,
}  # 50 m 16 s end to end


def _permille(seconds: dict[str, float], stages: tuple) -> dict[str, int]:
    """Shares of 1000 for ``stages``, in proportion to ``seconds``.

    Every stage gets at least 1, so the bar acknowledges each of them: a stage
    that measures below one part per thousand still has a band to begin and
    end. The rest is apportioned by largest remainder, which is what makes the
    result sum to **exactly** 1000 without a correction step — nothing has to
    be normalised at runtime, and the last band always ends where the bar does.
    A stage with no measurement is worth the floor and nothing more.
    """
    spare = 1000 - len(stages)
    total = sum(max(float(seconds.get(name, 0.0)), 0.0) for name in stages)
    if total <= 0.0:
        exact = {name: spare / len(stages) for name in stages}
    else:
        exact = {name: spare * max(float(seconds.get(name, 0.0)), 0.0) / total
                 for name in stages}
    share = {name: 1 + int(exact[name]) for name in stages}
    left = 1000 - sum(share.values())
    # Ties broken by pipeline order, so the table is the same on every machine.
    by_remainder = sorted(stages, key=lambda n: (-(exact[n] - int(exact[n])),
                                                  stages.index(n)))
    for name in by_remainder[:left]:
        share[name] += 1
    return share


def _starts(share: dict[str, int], stages: tuple) -> dict[str, int]:
    """Where each stage's band begins, in the same units."""
    start, running = {}, 0
    for name in stages:
        start[name] = running
        running += share[name]
    return start


_SUPPORTED: tuple[str, ...] = stages_for(True)

#: Each stage's share of a run without print supports, in parts per thousand.
STAGE_PERMILLE: dict[str, int] = _permille(SECONDS_PLAIN, STAGES)
#: Where each of those bands begins.
BAND_START: dict[str, int] = _starts(STAGE_PERMILLE, STAGES)

#: The same two tables for a run with print supports.
SUPPORT_PERMILLE: dict[str, int] = _permille(SECONDS_SUPPORT, _SUPPORTED)
SUPPORT_BAND_START: dict[str, int] = _starts(SUPPORT_PERMILLE, _SUPPORTED)

assert set(STAGE_PERMILLE) == set(STAGES), "a stage has no share of the bar"
assert set(SUPPORT_PERMILLE) == set(STAGES) | set(SUPPORT_STAGES)
assert sum(STAGE_PERMILLE.values()) == sum(SUPPORT_PERMILLE.values()) == 1000

_LAYOUTS: dict[tuple, tuple[dict, dict]] = {
    STAGES: (STAGE_PERMILLE, BAND_START),
    _SUPPORTED: (SUPPORT_PERMILLE, SUPPORT_BAND_START),
}


def layout(stages=None) -> tuple[dict[str, int], dict[str, int]]:
    """``(share, start)`` per stage, for a run that announced ``stages``.

    The two layouts a run can actually announce are precomputed. Anything else
    — a list from a newer or older pipeline than this window — is laid out on
    the spot from the supported run's measurements, which cover every stage
    there is; a name they do not know is worth the floor. The bar then still
    tiles from 0 to 1000 over exactly the stages it was told about, which is
    the property the window relies on.
    """
    if stages is None:
        return _LAYOUTS[STAGES]
    stages = tuple(stages)
    known = _LAYOUTS.get(stages)
    if known is None:
        share = _permille(SECONDS_SUPPORT, stages)
        known = _LAYOUTS[stages] = (share, _starts(share, stages))
    return known


def overall_permille(stage: str | None, fraction: float | None,
                     stages=None) -> int:
    """Where the top bar should sit while ``stage`` is ``fraction`` complete.

    ``stages`` is the list the run announced in its ``hello``, and it is what
    selects the layout: the plain pipeline's twelve, or the fifteen of a run
    with print supports. Omitted, it is the plain pipeline.

    ``fraction`` is ``None`` for a stage with no countable work — ``export``'s
    single ``STEPControl_Writer`` call, ``connect`` while it repairs
    interfaces, or ``simplify`` while its one dominant solid is being unified.
    The bar then **holds at the stage's band start** rather than creeping, and
    the sub-bar says what is happening instead. A bar that moved on a timer
    would be reporting a number that measures nothing; holding still while a
    label names the stage is the honest version of not knowing.
    """
    share, start = layout(stages)
    if stage is None or stage not in share:
        return 0
    if fraction is None:
        return start[stage]
    return start[stage] + int(round(share[stage] * min(max(fraction, 0.0), 1.0)))
