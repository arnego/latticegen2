"""Machine resource detection (specification.md §3, :mod:`latticegen2.sysinfo`).

Deliberately unmocked, unlike ``test_cli.py``'s ``--cores`` tests: those mock
this module out so their assertions are exact regardless of the machine, which
means nothing there would notice if the detection itself stopped working. This
file is the counterpart — it checks the real call on the real machine,
asserting only what must be true anywhere rather than anything about this
particular box.

There used to be a second pair of functions here, ``total_ram_gb`` and
``free_ram_gb``, behind the CLI's ``--ram`` budget. Both were removed along with
the flag (specification.md §11).
"""

from latticegen2 import sysinfo


def test_logical_core_count_is_a_usable_worker_count():
    cores = sysinfo.logical_core_count()
    assert isinstance(cores, int)
    assert cores >= 1
