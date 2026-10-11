"""Print orientation (specification.md §4.5, docs/algorithm.md §2).

The orientation of the lattice is set by the print bed, not by the coordinate
system of the input geometry. Three angles in degrees rotate the *part* relative
to the bed, about the bed's fixed X, then Y, then Z axis::

    R = Rz(rz) @ Ry(ry) @ Rx(rx)        p_bed = R @ p_part

The bed's normal is ``+Z``, so the build direction seen from the part's own
coordinate system is ``b = R.T @ z``. The input body is never transformed: the
lattice frame is rotated into the part's coordinates instead, which keeps the
output in the input's coordinate system (specification.md §1).

Pure NumPy, and imported by the command line, the lattice maths and the window
alike, so there is one definition of what the three angles mean.

**Zero angles give exactly the identity.** ``cos(0) == 1.0`` and
``sin(0) == 0.0`` exactly, and a product of matrices whose entries are exact
ones and zeros is exact. That is what lets the pipeline apply the rotation on
every run with no identity special case and still reproduce the unrotated
lattice (specification.md §4.5) — ``test/test_orient.py`` pins it.
"""

from __future__ import annotations

import numpy as np

Orientation = tuple[float, float, float]

IDENTITY: Orientation = (0.0, 0.0, 0.0)


def rotation_matrix(orient: Orientation) -> np.ndarray:
    """``R`` with ``p_bed = R @ p_part``, for ``(rx, ry, rz)`` in degrees."""
    rx, ry, rz = (float(np.radians(a)) for a in orient)
    cx, sx = float(np.cos(rx)), float(np.sin(rx))
    cy, sy = float(np.cos(ry)), float(np.sin(ry))
    cz, sz = float(np.cos(rz)), float(np.sin(rz))
    Rx = np.array([[1.0, 0.0, 0.0], [0.0, cx, -sx], [0.0, sx, cx]])
    Ry = np.array([[cy, 0.0, sy], [0.0, 1.0, 0.0], [-sy, 0.0, cy]])
    Rz = np.array([[cz, -sz, 0.0], [sz, cz, 0.0], [0.0, 0.0, 1.0]])
    return Rz @ Ry @ Rx


def build_direction(orient: Orientation) -> np.ndarray:
    """The bed's normal in part coordinates, ``b = R.T @ z``."""
    return rotation_matrix(orient).T @ np.array([0.0, 0.0, 1.0])


def euler_from_matrix(R: np.ndarray) -> Orientation:
    """``(rx, ry, rz)`` in degrees with ``rotation_matrix(...) == R``.

    The inverse the window needs when the part is dragged: the drag composes a
    rotation, and the three boxes must then show angles that reproduce it. At
    the gimbal singularity (``ry = ±90``) ``rx`` and ``rz`` are not separately
    determined, so ``rz`` is taken as zero and ``rx`` carries the whole of it.
    """
    R = np.asarray(R, dtype=float)
    sy = -float(R[2, 0])
    sy = min(1.0, max(-1.0, sy))
    ry = float(np.arcsin(sy))
    if abs(sy) < 1.0 - 1e-12:
        rx = float(np.arctan2(R[2, 1], R[2, 2]))
        rz = float(np.arctan2(R[1, 0], R[0, 0]))
    else:
        rz = 0.0
        rx = float(np.arctan2(-R[1, 2], R[1, 1]))
    return (float(np.degrees(rx)), float(np.degrees(ry)), float(np.degrees(rz)))


def is_rotated(orient: Orientation) -> bool:
    """Whether any angle is non-zero — for **naming only**, never for geometry.

    The file and part name carry the rotation only when there is one
    (specification.md §4.5). The geometry path does not ask this question: the
    rotation is applied whether or not it is the identity.
    """
    return any(float(a) != 0.0 for a in orient)
