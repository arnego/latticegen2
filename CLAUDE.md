# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**latticegen2** is a geometry generation tool that generates complex 3D lattice structures based on input parameters and a provided volume. 
Read @docs/specification.md
The file @docs/specification.md can be modified by claude following the rules defined in the top section of the document.
The detailed, normative algorithm implementation (exact math, pipeline diagrams, optimization strategy) is in @docs/algorithm.md — source code must match it exactly.

## Key Considerations

Goals from most important to least important:
1. Precision (that the output is correct)
2. Memory usage (that the script does not consume more memory than available and become unstable, other than that the use of all available resources is OK)
3. Speed (that the script completes in the shortest amount of time)

## Development Workflow

### Initial Setup

Before implementing features, establish:
- **Language/runtime**: **Python 3.11+**.
- **Dependencies** for geometry processing (note down selected geometries here)
  - Geometry kernel: **Open CASCADE (OCCT)**, accessed directly through the
    **OCP** bindings (`cadquery-ocp` wheels, which bundle the OCCT binaries).
    Used for STEP I/O, B-rep primitive construction, the single junction fuse,
    per-junction intersection, sewing, meshing, and `BRepCheck_Analyzer` validity
    checking. **NumPy** (array maths) is the only other dependency; core-count
    detection for the `--cores` budget (`src/latticegen2/sysinfo.py`) uses only
    the Python standard library. License texts live in
    `licenses/`, cross-referenced in `licenses/LICENSES.md`. See
    [docs/algorithm.md](docs/algorithm.md) for how the kernel is used in the
    pipeline, and note that the design deliberately keeps boolean operations off
    the hot path.
- **Testing framework** for geometry validation
  - `pytest` for unit tests (`test/test_*.py`), plus standalone scripts in
    `tools/` (`e2e.py`, `verify_geometry.py`) for end-to-end and geometry-validity
    checks. See [docs/testing.md](docs/testing.md).

### Review code

User will run /code-review manually or ask for it explicitly when needed. 

### Testing Procedures

Read @docs/testing.md 

### Releases

Offline release bundles and the tag-driven GitHub Actions workflow that
publishes them are documented in @docs/release.md. The version is single-sourced
from `src/latticegen2/__init__.py`; `pyproject.toml` derives it. Dependency pins
live in `requirements-bundle.txt` and must stay in sync with
`licenses/LICENSES.md`.

### Documentation

Maintain project documentation in README.md, including:
- Project introduction
- Dependencies list
- Installation instructions (if any)
- Parameter reference with valid ranges and defaults
- File format specifications for inputs
- File format specifications for outputs
- Example usage and typical workflows
- Document memory usage patterns and limits for large lattices
- Algorithm overview for the lattice generation method including implemented optimization strategy

## Tone Preference
Keep responses focused, brief, and concise. Lead with the outcome — the first
sentence should answer "what happened" or "what did I find" — with supporting
detail after it for anyone who wants more.

Before your first tool call, say in one sentence what you're about to do.
While working, give an update only when you find something important or
change direction. Skip narration of routine steps.

Only flag a correction to something you said earlier if the error would
change my code, conclusions, or decisions. State it plainly and move on.
For slips that change nothing for me, just fix it silently.

