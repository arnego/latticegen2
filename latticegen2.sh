#!/bin/sh
# Convenience wrapper: ./latticegen2.sh -i part.step -cc 10 -t 1.5 [options]
#
# Picks an interpreter in this order:
#   1. $LATTICEGEN2_PYTHON      explicit override always wins
#   2. ./runtime/bin/python3    portable release bundle, nothing installed
#   3. ./.venv/bin/python       wheels release bundle, after install.sh
#   4. python3                  a plain checkout on a prepared machine
#
# The tool runs straight from this directory with no install step, so an offline
# workstation needs either a release bundle or Python plus the two dependencies
# listed in README.md.
set -e
DIR=$(cd "$(dirname "$0")" && pwd)

if [ -n "$LATTICEGEN2_PYTHON" ]; then
    PY="$LATTICEGEN2_PYTHON"
elif [ -x "$DIR/runtime/bin/python3" ]; then
    PY="$DIR/runtime/bin/python3"
elif [ -x "$DIR/.venv/bin/python" ]; then
    PY="$DIR/.venv/bin/python"
else
    PY=python3
fi

# Always launch src/main.py rather than the installed console script: boundary
# workers use multiprocessing "spawn", which re-imports this __main__ module in
# each child, and main.py is what puts src/ on sys.path for them.
exec "$PY" "$DIR/src/main.py" "$@"
