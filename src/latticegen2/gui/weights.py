"""How much of the progress bar each pipeline stage is worth.

The twelve stages differ by three orders of magnitude in cost, so a bar that
advanced one twelfth per stage would sit at 42 % for the twenty minutes
``simplify`` takes and then jump. These shares come from a real run instead.

**Provenance.** The 2026-08-18 controlled pair in docs/profiling-reports.md —
`TD_HX_rehearsal_test` at ``cc=5, t=1`` on the 6-core development workstation,
the two halves run on the same machine with every untouched stage agreeing to
within 1.5 %, which is what makes averaging their columns legitimate rather than
optimistic. ``template`` and ``import`` measure 0 s there (the profiler samples
every 2 s) and are given a floor of one part per thousand, which the rounding
absorbs exactly.

**What this table is not.** It is one part's shape, on one machine. A small part
is boundary- and export-dominated — on ``80mm-test-ball`` ``simplify`` is
trivial — so the bar will visibly jump there. That is the accepted trade: a
fixed weighting is monotone and never lies about what has finished, where a
time-based estimate would sit at 100 % of a stage that had not ended. Read the
bar as "how far through the work", never as an ETA.

**And it rots.** If the rehearsal is ever re-profiled, this table is a downstream
consumer of those numbers; docs/profiling-reports.md says so beside the entry.
"""

from __future__ import annotations

from ..runlog import STAGES

#: Each stage's share of a run, in parts per thousand. Sums to exactly 1000, so
#: nothing has to be normalised at runtime.
STAGE_PERMILLE: dict[str, int] = {
    "template": 1,
    "import": 1,
    "tessellate": 1,
    "classify": 13,
    "boundary": 250,
    "connect": 4,
    "stitch": 207,
    "instance": 14,
    "assemble": 7,
    "simplify": 360,
    "validate": 35,
    "export": 107,
}

#: Where each stage's band begins, in the same units.
BAND_START: dict[str, int] = {}
_running = 0
for _name in STAGES:
    BAND_START[_name] = _running
    _running += STAGE_PERMILLE[_name]
del _running, _name

assert set(STAGE_PERMILLE) == set(STAGES), "a stage has no share of the bar"
assert sum(STAGE_PERMILLE.values()) == 1000


def overall_permille(stage: str | None, fraction: float | None) -> int:
    """Where the top bar should sit while ``stage`` is ``fraction`` complete.

    ``fraction`` is ``None`` for a stage with no countable work — ``export``'s
    single ``STEPControl_Writer`` call, or ``simplify`` while its one dominant
    solid is being unified. The bar then **holds at the stage's band start**
    rather than creeping, and the sub-bar says what is happening instead. A bar
    that moved on a timer would be reporting a number that measures nothing;
    holding still while a label names the stage is the honest version of not
    knowing.
    """
    if stage is None or stage not in STAGE_PERMILLE:
        return 0
    start = BAND_START[stage]
    if fraction is None:
        return start
    return start + int(round(STAGE_PERMILLE[stage] * min(max(fraction, 0.0), 1.0)))
