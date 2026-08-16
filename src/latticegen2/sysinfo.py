"""Machine resource detection for the CLI's optional budgets (specification.md §3).

``--cores`` and ``--ram`` are *budgets*, not hints: each is optional, and when
omitted it resolves to something the machine actually reports rather than to a
static literal. That resolution lives here, apart from :mod:`latticegen2.cli`'s
parsing so the parser stays free of platform-specific calls, and apart from
:mod:`latticegen2.parallel` so that module stays about the worker pool itself.

Core count comes from the stdlib. Memory does not — there is no stdlib call for
physical RAM — so ``psutil`` is a runtime dependency for exactly this, and only
this (see licenses/libraries.md). It is deliberately imported lazily inside the
two functions that need it: a machine-detection failure should surface as this
tool's own one-line parameter error, not as an import-time traceback out of a
module the user never asked about.
"""

from __future__ import annotations

import os

from .errors import ParamError


def logical_core_count() -> int:
    """Logical cores on this machine, the default worker count.

    Logical rather than physical (specification.md §3): boundary-junction jobs
    are constant-size, independent and process-parallel, so a hyperthread is
    still a place to run one. Falls back to 2 on the platforms where
    ``os.cpu_count()`` declines to answer — one worker plus one is a safe floor,
    never a zero.
    """
    return os.cpu_count() or 2


def _virtual_memory():
    """``psutil.virtual_memory()``, or a :class:`ParamError` naming the failure.

    Every caller here is resolving a ``--ram`` budget, which specification.md §7
    puts squarely before any computation starts — so a detection failure is a
    parameter failure (exit 2), not a resource one. No new exit code is needed
    for a dependency that is present in every supported install.
    """
    try:
        import psutil

        return psutil.virtual_memory()
    except Exception as exc:
        raise ParamError(
            f"Could not determine this machine's memory ({exc}). "
            f"--ram needs it both as its upper bound and as its default."
        ) from None


def total_ram_gb() -> float:
    """Total physical RAM in GB — the upper bound ``--ram`` is validated against.

    A budget above what the machine physically has is not a budget, so this is a
    hard ceiling rather than a warning.
    """
    return _virtual_memory().total / (1024.0 ** 3)


def free_ram_gb() -> float:
    """Available RAM in GB at this moment — the default when ``--ram`` is omitted.

    ``available`` rather than ``free``: it accounts for the cache and buffers the
    OS would reclaim under pressure, which is what "free" means to a process
    about to ask for memory. ``free`` would under-report it badly on a machine
    that has been up for a while, and a budget is more useful slightly generous
    than badly pessimistic.
    """
    return _virtual_memory().available / (1024.0 ** 3)
