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
# workstation needs either a release bundle or Python plus the dependencies
# listed in README.md.
set -e
DIR=$(cd "$(dirname "$0")" && pwd)

if [ -n "$LATTICEGEN2_PYTHON" ]; then
    PY="$LATTICEGEN2_PYTHON"
elif [ -d "$DIR/runtime" ]; then
    # A portable bundle. Never fall through to a system interpreter from here:
    # that one has none of the dependencies, so the operator would get a confusing
    # ModuleNotFoundError instead of the actual problem, which is almost always
    # an extraction that dropped the executable bit.
    PY="$DIR/runtime/bin/python3"
    if [ ! -x "$PY" ]; then
        echo "FAILED: $PY is not executable." >&2
        echo "  Extracting this bundle lost file permissions. Restore them with:" >&2
        echo "    chmod +x '$DIR/latticegen2.sh' '$DIR/runtime/bin/'*" >&2
        exit 1
    fi
elif [ -x "$DIR/.venv/bin/python" ]; then
    PY="$DIR/.venv/bin/python"
else
    PY=python3
fi

# Always launch src/main.py rather than the installed console script: boundary
# workers use multiprocessing "spawn", which re-imports this __main__ module in
# each child, and main.py is what puts src/ on sys.path for them.
#
# Not exec'd, because a non-zero exit needs one more look. It is normally the
# tool reporting a real failure, which is passed straight through; but it is also
# what an interpreter that cannot load the dependencies produces *before*
# main.py runs, and that case cannot be reported from inside Python -- a native
# library aborting during import kills the process outright (MKL's "cannot load
# mkl_intel_thread.2.dll" is not catchable as an exception) and exits 2,
# colliding with the parameter-error code. Probing only after something has gone
# wrong keeps the cost off every healthy run.
#
# The probe is scoped to codes that could actually be that. 3, 4, 5, 6 and 130
# come from the pipeline's own error paths (docs/algorithm.md §10), so by then
# the tool has demonstrably run and probing is pure delay -- which after a Ctrl+C
# would mean an unexplained pause following a shutdown the spec asks to be clean.
# Everything else can be an interpreter that never started: 1, 2, and the shell's
# own 126/127 for something it could not execute.
# No arguments at all means the graphical front-end (specification.md §3.1),
# but only where a window can exist: on a headless machine a bare invocation
# still exits 2 with usage, so nothing scripted changes meaning. `--gui`
# explicitly always tries, and reports one line if it cannot.
#
# The dependency probe runs *first* here, unlike the command line's, which
# probes only after a non-zero exit: this exec's, so there is no exit code to
# come back and inspect. tkinter is named because the front-end imports it at
# module scope (docs/algorithm.md §10).
if [ $# -eq 0 ] && { [ -n "$DISPLAY" ] || [ -n "$WAYLAND_DISPLAY" ]; }; then
    if "$PY" -c "import numpy, OCP.TopoDS, tkinter" >/dev/null 2>&1; then
        exec "$PY" "$DIR/src/main.py" --gui
    fi
    {
        echo
        echo "FAILED: the interpreter below cannot load latticegen2's dependencies,"
        echo "        so the graphical front-end cannot start."
        echo
        echo "  interpreter: $PY"
        echo
        "$PY" -c "import numpy, OCP.TopoDS, tkinter"
        echo
        echo "  On Debian and Ubuntu, tkinter is the separate python3-tk package."
        echo "  See README.md \"Installation\". A release bundle carries its own"
        echo "  interpreter and needs none of this."
    } >&2
    exit 1
fi

set +e
"$PY" "$DIR/src/main.py" "$@"
rc=$?
[ "$rc" -eq 0 ] && exit 0
case "$rc" in
    3|4|5|6|130) exit "$rc" ;;
esac
"$PY" -c "import numpy, OCP.TopoDS" >/dev/null 2>&1 && exit "$rc"

{
    echo
    echo "FAILED: the interpreter below cannot load latticegen2's dependencies, so"
    echo "        the tool never ran and the exit code above is not its own."
    echo
    echo "  interpreter: $PY"
    echo
    "$PY" -c "import numpy, OCP.TopoDS"
    echo
    echo "  Point LATTICEGEN2_PYTHON at a Python 3.11+ that has numpy and"
    echo "  cadquery-ocp:"
    echo "    export LATTICEGEN2_PYTHON=/path/to/python3"
    echo
    echo "  See README.md \"Installation\". A release bundle carries its own interpreter"
    echo "  and needs none of this."
} >&2
exit 1
