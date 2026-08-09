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
- **Dependencies** for geometry processing (note down selected geometries here)
  - Geometry kernel: **Gmsh SDK** (bundles Open CASCADE / OCCT), accessed from Julia via the `Gmsh.jl` package. Used for STEP/BREP I/O, B-rep primitive construction, and boolean operations (fuse/common). It is the only third-party Julia dependency; everything else is stdlib. License texts live in `licenses/` (GPL-2+ for Gmsh, LGPL-2.1 for the bundled OCCT), cross-referenced in `licenses/libraries.md`. See [docs/algorithm.md](docs/algorithm.md) for how it's used in the pipeline.
- **Testing framework** for geometry validation
  - Julia stdlib `Test` for unit tests (`test/runtests.jl`), plus standalone scripts in `tools/` (`e2e.jl`, `verify_geometry.jl`) for end-to-end and geometry-validity checks. See [docs/testing.md](docs/testing.md).

### Review code

User will run /code-review manually or ask for it explicitly when needed. 

### Testing Procedures

Read @docs/testing.md 

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

