#!/usr/bin/env bash
# Convenience wrapper for Linux. See specification.md §2 (Packaging form)
# and README.md for usage. Forwards all arguments to src/main.jl using the
# Julia project defined by Project.toml/Manifest.toml in this directory.
#
# `-t auto` enables Julia's own thread pool (distinct from the --cores/
# --workers Distributed *process* pool used for tiling, docs/algorithm.md
# §6.2): the classification stage's per-strut loop is threaded
# (docs/algorithm.md §11.2) and is otherwise silently single-threaded, since
# Julia defaults to 1 thread unless told otherwise.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec julia -t auto --project="$SCRIPT_DIR" "$SCRIPT_DIR/src/main.jl" "$@"
