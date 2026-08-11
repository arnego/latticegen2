"""latticegen2 — parameterised diamond-strut lattice generation.

Generates a strut lattice filling the solid volume of an input STEP file and
writes it back out as an exact B-rep AP214 STEP body.

The architecture is documented in docs/algorithm.md; the short version is that
the lattice is built by *instancing* one pre-fused junction solid at every
lattice node and stitching the instances along their shared mid-strut faces, so
the only boolean operations in the program are one six-operand fuse (the
template) and one single-operand intersection per boundary junction.
"""

__version__ = "2.0.0"
