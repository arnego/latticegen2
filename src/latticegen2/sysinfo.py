"""Machine resource detection for the CLI's optional budgets (specification.md §3).

``--cores`` and ``--ram`` are *budgets*, not hints: each is optional, and when
omitted it resolves to something the machine actually reports rather than to a
static literal. That resolution lives here, apart from :mod:`latticegen2.cli`'s
parsing so the parser stays free of platform-specific calls, and apart from
:mod:`latticegen2.parallel` so that module stays about the worker pool itself.

Core count comes from the stdlib. Memory does not — there is no stdlib call for
physical RAM — so ``psutil`` is a runtime dependency for exactly this, and only
this (see licenses/LICENSES.md).

It is imported at module scope, exactly like ``numpy`` and ``OCP`` elsewhere in
the package, and that is a decision rather than a default. Deferring it into the
functions that need it would keep the process alive past a missing dependency
and report it as a *parameter* error, which is the wrong name for it: nothing is
wrong with the parameter. It would also make the launchers lie — they probe the
interpreter after a non-zero exit and say "the tool never ran" (docs/algorithm.md
§10), which is true of a dependency that fails at import and false of one that
fails on use. Importing here puts psutil in the same class as every other
dependency: absent means the tool does not start, and the launcher says so and
names the interpreter it chose.
"""

from __future__ import annotations

import os

import psutil


def logical_core_count() -> int:
    """Logical cores on this machine, the default worker count.

    Logical rather than physical (specification.md §3): boundary-junction jobs
    are constant-size, independent and process-parallel, so a hyperthread is
    still a place to run one. Falls back to 2 on the platforms where
    ``os.cpu_count()`` declines to answer — one worker plus one is a safe floor,
    never a zero.
    """
    return os.cpu_count() or 2


def total_ram_gb() -> float:
    """Total physical RAM in GB — the upper bound ``--ram`` is validated against.

    A budget above what the machine physically has is not a budget, so this is a
    hard ceiling rather than a warning.
    """
    return psutil.virtual_memory().total / (1024.0 ** 3)


def free_ram_gb() -> float:
    """Available RAM in GB at this moment — the default when ``--ram`` is omitted.

    ``available`` rather than ``free``: it accounts for the cache and buffers the
    OS would reclaim under pressure, which is what "free" means to a process
    about to ask for memory. ``free`` would under-report it badly on a machine
    that has been up for a while, and a budget is more useful slightly generous
    than badly pessimistic.
    """
    return psutil.virtual_memory().available / (1024.0 ** 3)
