"""Machine resource detection for the CLI's optional budget (specification.md §3).

``--cores`` is a *budget*, not a hint: it is optional, and when omitted it
resolves to something the machine actually reports rather than to a static
literal. That resolution lives here, apart from :mod:`latticegen2.cli`'s
parsing so the parser stays free of platform-specific calls, and apart from
:mod:`latticegen2.parallel` so that module stays about the worker pool itself.

Core count comes from the stdlib, which is why this module has no third-party
dependency of its own (specification.md §11 records the removal of the former
``--ram`` budget and the ``psutil`` dependency that existed only to serve it).
"""

from __future__ import annotations

import os


def logical_core_count() -> int:
    """Logical cores on this machine, the default worker count.

    Logical rather than physical (specification.md §3): boundary-junction jobs
    are constant-size, independent and process-parallel, so a hyperthread is
    still a place to run one. Falls back to 2 on the platforms where
    ``os.cpu_count()`` declines to answer — one worker plus one is a safe floor,
    never a zero.
    """
    return os.cpu_count() or 2
