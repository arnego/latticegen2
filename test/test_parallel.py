"""The shared worker pool (docs/algorithm.md §8, §12).

:class:`latticegen2.parallel.WorkerPool` replaced two near-identical transient
`multiprocessing` pools (one in :mod:`latticegen2.boundary`, one in
:mod:`latticegen2.weld`) with one pool built once per run and handed to every
parallel stage. These tests pin its own contract in isolation, with a trivial
worker function rather than real geometry — kernel-free, so they stay in the
fast subset docs/testing.md names for iterative development — plus one test
that a real ``spawn`` pool actually runs a job in a separate process, mirroring
how `test_weld.py`'s tiling test is the one place that crosses a process
boundary for real.
"""

from __future__ import annotations

import pytest

from latticegen2.errors import ProcessingError
from latticegen2.parallel import WorkerPool

# Worker functions used with a real `spawn` pool must be importable at module
# scope — a closure or lambda cannot be pickled to send to a child process.


def _times_two_with_rss(job):
    """``(value, rss)`` -> ``(value * 2, rss)``, mimicking every real worker's
    "small plain data plus a peak-RSS reading" convention (docs/algorithm.md §7)."""
    value, rss = job
    return value * 2, rss


def _boom(job):
    raise ValueError(f"boom {job}")


def _os_pid(_job):
    import os

    return os.getpid(), 0


# --- construction and the workers<=1 contract --------------------------------


def test_workers_le_1_builds_no_real_pool():
    with WorkerPool(1) as pool:
        assert not pool.active


def test_zero_or_negative_workers_also_builds_no_real_pool():
    # `cli.py` already floors the derived worker count at 1, but this class
    # should not rely on that — a caller passing 0 must degrade the same way,
    # not crash.
    with WorkerPool(0) as pool:
        assert not pool.active


def test_run_refuses_when_no_pool_is_active():
    with WorkerPool(1) as pool:
        with pytest.raises(ProcessingError, match="no pool active"):
            pool.run(_times_two_with_rss, [(1, 0)])


def test_workers_greater_than_one_builds_a_real_pool():
    with WorkerPool(2) as pool:
        assert pool.active


def test_pool_is_inert_again_after_the_with_block_exits():
    pool = WorkerPool(2)
    with pool:
        assert pool.active
    assert not pool.active


# --- run(): ordering, RSS folding, real parallelism ---------------------------


def test_results_come_back_in_input_order_not_completion_order():
    jobs = [(i, i) for i in range(8)]  # rss = i, so later jobs "look" heavier
    with WorkerPool(4) as pool:
        results, max_rss = pool.run(_times_two_with_rss, jobs)
    assert results == [(i * 2, i) for i in range(8)]
    assert max_rss == 7


def test_sort_by_changes_dispatch_but_not_the_returned_order():
    """The one thing `sort_by` may never do is reorder what the caller sees.

    Dispatching largest-first is purely a scheduling optimization
    (:meth:`WorkerPool.run`'s own docstring) — a caller that already built its
    output list assuming ``results[i]`` answers ``jobs[i]`` must get the same
    answer whether or not a sort key was given.
    """
    jobs = [(v, 0) for v in [3, 1, 4, 1, 5, 9, 2, 6]]
    with WorkerPool(4) as pool:
        plain, _ = pool.run(_times_two_with_rss, jobs)
        sorted_dispatch, _ = pool.run(_times_two_with_rss, jobs, sort_by=lambda j: j[0])
    assert plain == sorted_dispatch == [(v * 2, 0) for v in [3, 1, 4, 1, 5, 9, 2, 6]]


def test_max_rss_is_the_maximum_not_the_last():
    jobs = [(0, 5), (0, 90), (0, 12)]
    with WorkerPool(3) as pool:
        _results, max_rss = pool.run(_times_two_with_rss, jobs)
    assert max_rss == 90


def test_a_worker_exception_propagates_to_the_caller():
    with WorkerPool(2) as pool:
        with pytest.raises(ValueError, match="boom"):
            pool.run(_boom, [1, 2, 3])


def test_a_job_actually_runs_in_a_different_process():
    """The one real cross-process check, mirroring test_weld.py's own note on
    why this is exercised for real rather than only through the sequential
    branch (spawning processes per test is slow, so this is deliberately the
    single test that pays for it)."""
    import os

    with WorkerPool(2) as pool:
        results, _ = pool.run(_os_pid, [None, None])
    pids = {pid for pid, _rss in results}
    assert os.getpid() not in pids
