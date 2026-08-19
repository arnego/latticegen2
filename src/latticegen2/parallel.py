"""Shared spawn-process worker pool infrastructure.

Every stage that distributes work across processes — classification
(:mod:`latticegen2.classify`), boundary trimming
(:mod:`latticegen2.boundary`), the boundary-sew tiling and its round-2 merge
(:mod:`latticegen2.weld`), same-domain unification
(:mod:`latticegen2.pipeline`) — follows the same discipline: a ``spawn``-context
pool, ordered ``imap`` so results land in job order regardless of which worker
finishes first (reproducibility matters more than the marginal speed of
consuming results as they arrive — docs/algorithm.md §7, §8), and only file
paths and small plain data crossing the process boundary, never live OCCT
geometry.

That pattern used to be duplicated once per stage, with a fresh pool built and
torn down each time. :class:`WorkerPool` centralises it: one pool is built in
:func:`latticegen2.pipeline._run` for the whole run and handed to whichever
stage wants it, so ``spawn``'s process-creation cost — non-trivial on a part
with many independent tiled components, as `TD_HX_rehearsal_test`'s 14 is
(specification.md §10) — is paid once per run, not once per stage.

A caller whose own job count is too small to be worth parallelising (``workers
<= 1``, or too few jobs to amortise `spawn`'s cost) should not build a
:class:`WorkerPool` at all — that decision is job-shape-specific (a boundary
batch is not a sew tile is not a solid to unify) and stays with the caller,
exactly as it did before this module existed.
"""

from __future__ import annotations

import multiprocessing as mp
import os
from typing import Callable

from OCP.BRep import BRep_Builder
from OCP.BRepTools import BRepTools
from OCP.OSD import OSD_ThreadPool
from OCP.TopoDS import TopoDS_Iterator, TopoDS_Shape

from .errors import ProcessingError


def set_background_priority() -> None:
    """Drop this process to below-normal scheduling priority.

    Unconditional for every process in a run, master and workers alike. It used
    to be the opt-in ``-bg`` flag (specification.md §3); making it the only
    behaviour removed a choice with no interesting second option — the cost is a
    scheduling hint on a machine with nothing else to run, and the benefit is
    that a long run never makes the desktop unusable.
    """
    try:
        if os.name == "nt":
            import ctypes

            BELOW_NORMAL_PRIORITY_CLASS = 0x00004000
            ctypes.windll.kernel32.SetPriorityClass(
                ctypes.windll.kernel32.GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS
            )
        else:
            os.nice(5)
    except Exception:
        pass  # priority is a courtesy, never a reason to fail a run


def set_thread_budget(cores: int) -> None:
    """Cap OCCT's *own* native thread pool to the ``--cores`` budget.

    Distinct from :class:`WorkerPool`, and the distinction is the whole point.
    That pool is processes this code dispatches to; this is the thread pool
    OCCT launches inside a single call when asked — today only
    :func:`latticegen2.occ.is_valid`, via ``BRepCheck_Analyzer``'s
    ``theIsParallel`` (docs/algorithm.md §9). Left at its default, that pool
    sizes itself to the *machine*, which would quietly make ``--cores 2`` on a
    six-core box use six threads and break the plain reading of
    specification.md §3: "Maximum CPU cores this run may use", honoured exactly.

    ``SetNbDefaultThreadsToLaunch`` caps what one launcher may lock rather than
    resizing the pool, which is what makes it safe to call on an already-created
    default pool — ``Init`` would throw if any job were active.

    This is deliberately *not* called in the worker initializer. Nothing
    dispatched to a worker asks OCCT for threads, so there is nothing there to
    bound; and were that to change, the budget would have to be divided by the
    worker count rather than repeated in each one.
    """
    try:
        OSD_ThreadPool.DefaultPool_s().SetNbDefaultThreadsToLaunch(max(1, int(cores)))
    except Exception:
        pass  # a thread-count hint is never a reason to fail a run


def read_brep(path: str) -> TopoDS_Shape:
    """Read one intermediate ``.brep`` file, raising :class:`ProcessingError`."""
    shape = TopoDS_Shape()
    builder = BRep_Builder()
    if not BRepTools.Read_s(shape, path, builder):
        raise ProcessingError(f"Could not read intermediate BREP file: {path}")
    return shape


def write_brep(shape: TopoDS_Shape, path: str) -> None:
    """Write one intermediate ``.brep`` file."""
    BRepTools.Write_s(shape, path)


def compound_children(shape: TopoDS_Shape) -> list[TopoDS_Shape]:
    """Direct children of a compound, in the order they were added."""
    out = []
    it = TopoDS_Iterator(shape)
    while it.More():
        out.append(it.Value())
        it.Next()
    return out


def _pool_initializer() -> None:
    """Run once per worker process, before it takes any jobs.

    Priority is a process-lifetime setting, not a per-job one, so it belongs in
    one-time setup: paying for :func:`set_background_priority` on every job (as
    each stage used to, when it built its own transient pool) was only ever
    wasted work.
    """
    set_background_priority()


class WorkerPool:
    """One ``spawn``-context process pool, shared across every stage of a run.

    ``workers <= 1`` builds no real pool at all — :attr:`active` is ``False``
    and :meth:`run` refuses to be called — so a caller that already checked
    ``workers <= 1`` and took its own sequential path never has to ask this
    class for a pool it does not need.
    """

    def __init__(self, workers: int):
        self.workers = workers
        self._pool = None

    @property
    def active(self) -> bool:
        return self._pool is not None

    def __enter__(self) -> "WorkerPool":
        if self.workers > 1:
            ctx = mp.get_context("spawn")
            self._pool = ctx.Pool(
                processes=self.workers,
                initializer=_pool_initializer,
            )
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:
        if self._pool is not None:
            # Delegates to Pool's own context-manager protocol (unconditional
            # `terminate()`) rather than reimplementing it, so a Ctrl+C mid-run
            # tears the workers down exactly as every transient pool this
            # replaces always did — caught at the top level as CancelledError
            # (docs/algorithm.md §10).
            self._pool.__exit__(exc_type, exc, tb)
            self._pool = None
        return False

    def run(
        self, worker_fn: Callable, jobs: list, *, rss_index: int = -1,
        sort_by: Callable | None = None,
        on_result: Callable[[int, int], None] | None = None,
    ) -> tuple[list, int]:
        """Run ``worker_fn`` over ``jobs``, folding each result's peak-RSS.

        Ordered ``imap`` (``chunksize=1``): jobs still run in parallel, but
        results come back in job order, not completion order, so a caller
        building output from them gets the same order every run regardless of
        which worker happened to finish first — the same reproducibility
        argument :func:`latticegen2.boundary.trim_boundary` and
        :func:`latticegen2.weld._sew_all_tiles` document at their own call
        sites.

        Every worker result is expected to carry its own
        :func:`latticegen2.runlog.peak_rss_bytes` reading at ``rss_index``
        (the last element, by convention shared with every worker function in
        this codebase); this method only extracts that column for the fold —
        it does not otherwise interpret the tuple, so an existing worker
        result shape needs no change to be used through this class.

        ``sort_by``, if given, is a key function used to submit ``jobs`` to
        the pool largest-first (descending) rather than in their given order,
        for callers whose jobs are very unequal in size — same-domain
        unification's 14 wildly different solids on the
        docs/specification.md §10 rehearsal is the motivating case. `imap`
        assigns queued tasks to workers as they idle, so a large job sitting
        late in submission order can start only once an earlier, unrelated
        job frees a worker, stretching the tail of the run past what six busy
        workers should cost; submitting it first lets it start immediately
        instead. The returned list is always back in ``jobs``' own order
        regardless of ``sort_by`` — only *dispatch* order changes, per this
        method's own ordering guarantee above.

        ``on_result``, if given, is called ``(completed, total)`` once per
        finished job, from inside the loop below. It is what lets a stage report
        progress *while* it runs rather than after it: every parallel stage in
        the pipeline dispatches through this one loop, so one hook here serves
        classification, boundary trimming, boundary-sew round 1 and same-domain
        unification alike (docs/algorithm.md §12).

        **It is an observer and nothing else, and the ordering guarantee above
        is what makes that safe.** It receives two integers, is called after the
        result has been stored, and cannot reach ``results``. In particular
        ``imap`` stays *ordered* with ``chunksize=1``: switching to
        ``imap_unordered`` would make progress arrive sooner and would break the
        byte-reproducibility this method exists to provide, so it must not be
        done for the front-end's benefit. One consequence to keep in mind at the
        call site: the count is a *completion count*, not job identity — a slow
        early job holds back the ticks for every job behind it, and a caller
        using ``sort_by`` cannot infer *which* jobs are done.
        """
        if not self.active:
            raise ProcessingError(
                "WorkerPool.run() called with no pool active (workers <= 1). "
                "A caller in that situation should take its own sequential "
                "path instead of asking this class for a pool it does not have."
            )
        if sort_by is not None:
            order = sorted(range(len(jobs)), key=lambda i: sort_by(jobs[i]), reverse=True)
        else:
            order = list(range(len(jobs)))
        dispatch = [jobs[i] for i in order]

        results: list = [None] * len(jobs)
        max_rss = 0
        completed = 0
        for i, out in zip(order, self._pool.imap(worker_fn, dispatch, chunksize=1)):
            results[i] = out
            max_rss = max(max_rss, out[rss_index])
            completed += 1
            if on_result is not None:
                on_result(completed, len(jobs))
        return results, max_rss
