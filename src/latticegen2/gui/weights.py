"""How much of the progress bar each pipeline stage is worth.

The twelve stages differ by three orders of magnitude in cost, so a bar that
advanced one twelfth per stage would sit at 42 % for the twenty minutes
``simplify`` takes and then jump. These shares come from a real run instead.

**Provenance.** The 2026-08-26 entry in docs/profiling-reports.md —
`TD_HX_rehearsal_test` at ``cc=5, t=1`` on the 6-core development workstation,
90.9 minutes end to end. ``template`` measures 0 s (the profiler samples every
2 s) and ``import``/``tessellate`` round below one part per thousand; all three
are given a floor of 1, and the two largest stages absorb the two parts that
costs.

**This one is a single run, not a controlled pair, and that is a weaker basis
than the table it replaces** — the 2026-08-18 pair could be averaged precisely
because its untouched stages agreed to within 1.5 %. It is used anyway because
the alternative is worse: the previous shares gave ``validate`` **35** of 1000
when it is now the largest stage in the run at **385**, so the bar would sit
still for thirty-five minutes and then jump. A weighting measured once is wrong
about the margins; the one it replaces was wrong about the shape.

**What this table is not.** It is one part's shape, on one machine. A small part
is boundary- and export-dominated — on ``80mm-test-ball`` ``simplify`` is
trivial — so the bar will visibly jump there. ``validate``'s 385 is dominated by
the export-truth gate (docs/algorithm.md §9), which writes, re-reads and
tessellates **every** output solid: 2,003 s of that stage's 2,105 s here. On a
part whose bodies are small that gate is seconds, and the bar will run through
``validate`` almost at once. That is the accepted trade: a
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
    "classify": 5,
    "boundary": 211,
    "connect": 6,
    "stitch": 108,
    "instance": 8,
    "assemble": 4,
    "simplify": 204,
    "validate": 385,
    "export": 66,
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

#: Shares of the three stages a run with print supports adds, in the same units
#: as :data:`STAGE_PERMILLE` *before* the whole is renormalised to 1000.
#:
#: **Provenance, and it is thin:** one run of `test-cylinder.STEP` at
#: ``cc=10, t=1.5`` on the development workstation, where ``plan`` took 12 s and
#: ``support`` 21 s beside a 7.4 s ``boundary``. That part's whole run is two
#: minutes, so these are a shape, not a calibration — read them the way the
#: table above says to read it.
SUPPORT_SHARES: dict[str, int] = {"overhang": 10, "plan": 150, "support": 400}


def _bands(stages) -> tuple[dict, dict]:
    """``(share, start)`` per stage for a run that announced ``stages``."""
    raw = {name: STAGE_PERMILLE.get(name, SUPPORT_SHARES.get(name, 1))
           for name in stages}
    total = float(sum(raw.values())) or 1.0
    share, start, running = {}, {}, 0.0
    for name in stages:
        start[name] = running
        share[name] = raw[name] * 1000.0 / total
        running += share[name]
    return share, start


def overall_permille(stage: str | None, fraction: float | None,
                     stages=None) -> int:
    """Where the top bar should sit while ``stage`` is ``fraction`` complete.

    ``fraction`` is ``None`` for a stage with no countable work — ``export``'s
    single ``STEPControl_Writer`` call, or ``simplify`` while its one dominant
    solid is being unified. The bar then **holds at the stage's band start**
    rather than creeping, and the sub-bar says what is happening instead. A bar
    that moved on a timer would be reporting a number that measures nothing;
    holding still while a label names the stage is the honest version of not
    knowing.
    """
    if stages is not None and tuple(stages) != STAGES:
        # A run with print supports goes through more stages than the plain
        # pipeline; it says which in its `hello`, and the bar is laid out over
        # exactly those. The plain pipeline keeps the table above untouched.
        if stage is None or stage not in stages:
            return 0
        share, begin = _bands(tuple(stages))
        if fraction is None:
            return int(round(begin[stage]))
        return int(round(begin[stage] + share[stage] * min(max(fraction, 0.0), 1.0)))
    if stage is None or stage not in STAGE_PERMILLE:
        return 0
    start = BAND_START[stage]
    if fraction is None:
        return start
    return start + int(round(STAGE_PERMILLE[stage] * min(max(fraction, 0.0), 1.0)))
