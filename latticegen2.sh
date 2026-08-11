#!/bin/sh
# Convenience wrapper: ./latticegen2.sh -i part.step -cc 10 -t 1.5 [options]
#
# Uses $LATTICEGEN2_PYTHON if set, otherwise python3. The tool runs straight
# from this checkout, with no install step, so an offline workstation only needs
# Python plus the two dependencies listed in README.md.
set -e
DIR=$(cd "$(dirname "$0")" && pwd)
"${LATTICEGEN2_PYTHON:-python3}" "$DIR/src/main.py" "$@"
