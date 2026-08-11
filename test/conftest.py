"""Make the package importable from a bare checkout, with no install step."""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

TESTDIR = os.path.dirname(os.path.abspath(__file__))
