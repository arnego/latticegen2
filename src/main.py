#!/usr/bin/env python3
"""Runnable entry point: ``python src/main.py -i part.step -cc 10 -t 1.5``.

Adds ``src/`` to ``sys.path`` so the tool runs straight from a checkout with no
installation step, which is what the offline deployment target needs
(specification.md §2). An installed copy can use the ``latticegen2`` console
script or ``python -m latticegen2`` instead.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from latticegen2.__main__ import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
