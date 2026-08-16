"""Machine resource detection (specification.md §3, :mod:`latticegen2.sysinfo`).

Deliberately unmocked, unlike ``test_cli.py``'s budget tests: those mock this
module out so their assertions are exact regardless of the machine, which means
nothing there would notice if the detection itself stopped working. This file is
the counterpart — it checks the real calls on the real machine, asserting only
what must be true anywhere rather than anything about this particular box.
"""

from latticegen2 import sysinfo


def test_logical_core_count_is_a_usable_worker_count():
    cores = sysinfo.logical_core_count()
    assert isinstance(cores, int)
    assert cores >= 1


def test_memory_detection_returns_sane_figures():
    total = sysinfo.total_ram_gb()
    free = sysinfo.free_ram_gb()
    # A machine running this test has at least a fraction of a GB, and cannot
    # have more free than it has in total. Both are the properties the CLI's
    # bounds actually rest on; the absolute figures are the machine's business.
    assert total > 0.1
    assert 0.0 <= free <= total
