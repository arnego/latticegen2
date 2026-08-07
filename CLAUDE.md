# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**latticegen2** is a geometry generation tool that generates complex 3D lattice structures based on input parameters and a provided volume. 
Read @docs/specification.md
The file @docs/specification.md can be modified by claude following the rules defined in the top section of the document.

## Key Considerations

Goals from most important to least important:
1. Precision (that the output is correct)
2. Memory usage (that the script does not consume more memory than available and become unstable, other than that the use of all available resources is OK)
3. Speed (that the script completes in the shortest amount of time)

## Development Workflow

### Initial Setup

Before implementing features, establish:
- **Dependencies** for geometry processing (note down selected geometries here)
- **Testing framework** for geometry validation

### Review code

Run /code-review before pushing and fix any issues that come up. 

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

