"""Lattice mathematics.

Implements docs/algorithm.md §2 and §3.1 exactly — strut directions, the node
basis, the diamond profile frame, and candidate enumeration from a world-space
AABB. Every angle and length comes from an expression (``asin(sqrt(2/3))``),
never a decimal literal, so precision is IEEE-754 double throughout.

This module is pure NumPy: it never touches the geometry kernel, which is what
makes the lattice math unit-testable on its own (``test/test_lattice.py``).

The one concept this module adds beyond docs/algorithm.md is the **half-strut**:
the half of a strut adjacent to one node, obtained by cutting every strut at its
midpoint. Each node owns exactly six of them (``±e_k`` for ``k = 0, 1, 2``), and
their union is the junction solid ``J`` that the whole fuse-free architecture is
built on — see docs/algorithm.md §3.2.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# --- Strut recline angle (docs/algorithm.md §2.1) --------------------------

THETA = float(np.arcsin(np.sqrt(2.0 / 3.0)))
"""Exact strut recline angle from +Z, ~54.7356°, computed natively."""

SIN_THETA = float(np.sqrt(2.0 / 3.0))
COS_THETA = float(1.0 / np.sqrt(3.0))


def strut_direction(k: int) -> np.ndarray:
    """Unit direction of strut family ``k in {0,1,2}`` (docs/algorithm.md §2.1).

    Azimuth ``2*pi*k/3`` around Z, reclined from +Z by :data:`THETA`.
    """
    phi = 2.0 * np.pi * k / 3.0
    return np.array(
        [SIN_THETA * np.cos(phi), SIN_THETA * np.sin(phi), COS_THETA], dtype=float
    )


def cube_edge(cc: float) -> float:
    """Cube edge length ``a = cc / sqrt(2)`` (docs/algorithm.md §2.2)."""
    return float(cc) / np.sqrt(2.0)


# --- Precomputed lattice geometry ------------------------------------------


@dataclass(frozen=True)
class LatticeParams:
    """All lattice geometry derived from the CLI parameters ``cc`` and ``t``.

    Attributes are indexed ``0:3`` for strut families ``e0, e1, e2``, matching
    the ``k`` of docs/algorithm.md.

    ``e``, ``u`` and ``v`` are ``(3, 3)`` arrays whose **rows** are the
    per-direction vectors; ``B`` is the ``(3, 3)`` basis matrix whose
    **columns** are ``a * e_k``, so ``node(i,j,k) = B @ [i,j,k]``.
    """

    cc: float
    t: float
    a: float
    """Cube edge length, ``cc / sqrt(2)``."""
    r: float
    """Strut circumradius, ``t / sqrt(2)`` — half the profile's diagonal."""
    e: np.ndarray
    u: np.ndarray
    """Horizontal profile axis per direction (docs/algorithm.md §3.1)."""
    v: np.ndarray
    """Vertical-plane profile axis per direction (docs/algorithm.md §3.1)."""
    B: np.ndarray


def lattice_params(cc: float, t: float) -> LatticeParams:
    """Build a :class:`LatticeParams` from the two CLI geometry parameters."""
    a = cube_edge(cc)
    e = np.array([strut_direction(k) for k in range(3)], dtype=float)
    zhat = np.array([0.0, 0.0, 1.0])
    u = np.array([_normalize(np.cross(zhat, e[k])) for k in range(3)], dtype=float)
    v = np.array([np.cross(e[k], u[k]) for k in range(3)], dtype=float)
    B = a * np.column_stack([e[0], e[1], e[2]])
    return LatticeParams(cc=float(cc), t=float(t), a=a, r=float(t) / np.sqrt(2.0),
                         e=e, u=u, v=v, B=B)


def _normalize(x: np.ndarray) -> np.ndarray:
    return x / np.linalg.norm(x)


def node(lp: LatticeParams, ijk) -> np.ndarray:
    """World-space position of lattice node ``(i,j,k)`` (docs/algorithm.md §2.2)."""
    return lp.B @ np.asarray(ijk, dtype=float)


def nodes(lp: LatticeParams, ijk: np.ndarray) -> np.ndarray:
    """Vectorized :func:`node`: ``(N, 3)`` integer indices -> ``(N, 3)`` positions."""
    return np.asarray(ijk, dtype=float) @ lp.B.T


def profile_vertices(lp: LatticeParams, center: np.ndarray, k: int) -> np.ndarray:
    """The four corners of the diamond cross-section centred at ``center``.

    Order is ``+u, +v, -u, -v`` (docs/algorithm.md §3.1), which traverses the
    square consistently — consecutive entries are adjacent corners, so the list
    is directly usable as a closed polygon.
    """
    hu = lp.r * lp.u[k]
    hv = lp.r * lp.v[k]
    c = np.asarray(center, dtype=float)
    return np.array([c + hu, c + hv, c - hu, c - hv], dtype=float)


# --- Half-struts ------------------------------------------------------------

HALF_STRUTS: tuple[tuple[int, int], ...] = (
    (0, +1), (1, +1), (2, +1), (0, -1), (1, -1), (2, -1),
)
"""The six ``(k, sign)`` half-struts every node owns, in canonical order.

Index ``h`` in this tuple is the half-strut's id throughout the pipeline;
``0..2`` are the outgoing (``+e_k``) halves and ``3..5`` the incoming ones, so
``h`` and ``(h + 3) % 6`` are always the two halves pointing opposite ways along
the same axis.
"""


def half_strut_offset(lp: LatticeParams, h: int) -> np.ndarray:
    """Vector from a node to the far (cap) end of half-strut ``h``."""
    k, s = HALF_STRUTS[h]
    return s * (lp.a / 2.0) * lp.e[k]


def neighbor_step(h: int) -> np.ndarray:
    """Index-space step to the node sharing half-strut ``h``'s cap.

    Half-strut ``h`` of node ``n`` and half-strut ``(h + 3) % 6`` of node
    ``n + neighbor_step(h)`` are the two halves of one strut, and their cap
    quads are the same square in world space.
    """
    k, s = HALF_STRUTS[h]
    step = np.zeros(3, dtype=np.int64)
    step[k] = s
    return step


OPPOSITE_HALF: tuple[int, ...] = tuple((h + 3) % 6 for h in range(6))
"""``OPPOSITE_HALF[h]`` is the half-strut facing ``h`` across the shared cap."""


# --- Candidate index range (docs/algorithm.md §2.4) -------------------------


def index_range(lp: LatticeParams, lo: np.ndarray, hi: np.ndarray):
    """Inclusive ``(lo_idx, hi_idx)`` node-index bounds covering AABB ``[lo, hi]``.

    Ported verbatim from docs/algorithm.md §2.4: map the box's eight corners
    into index space through ``B``, take the floor/ceil hull, and pad by
    ``ceil(r/a) + 1`` cells so no strut that could geometrically touch the
    (circumradius-expanded) volume is missed.

    The padding is comfortably more than the ``a/2`` reach of a junction, so a
    node outside the returned range cannot contribute material inside the box.
    """
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    corners = np.array(
        [[x, y, z] for x in (lo[0], hi[0]) for y in (lo[1], hi[1]) for z in (lo[2], hi[2])],
        dtype=float,
    )
    idx_corners = np.linalg.solve(lp.B, corners.T).T
    pad = int(np.ceil(lp.r / lp.a)) + 1
    lo_idx = np.floor(idx_corners.min(axis=0)).astype(np.int64) - pad
    hi_idx = np.ceil(idx_corners.max(axis=0)).astype(np.int64) + pad
    return lo_idx, hi_idx


def candidate_nodes(lp: LatticeParams, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """Every node index in :func:`index_range`, as an ``(N, 3)`` int array.

    Rows are ordered so the last index varies fastest, which keeps nodes that
    are spatial neighbours close together in memory for the classification
    sweep.
    """
    lo_idx, hi_idx = index_range(lp, lo, hi)
    ii = np.arange(lo_idx[0], hi_idx[0] + 1, dtype=np.int64)
    jj = np.arange(lo_idx[1], hi_idx[1] + 1, dtype=np.int64)
    kk = np.arange(lo_idx[2], hi_idx[2] + 1, dtype=np.int64)
    grid = np.meshgrid(ii, jj, kk, indexing="ij")
    return np.stack([g.ravel() for g in grid], axis=1)


def format_param(x: float) -> str:
    """Format a float for file names and part names without trailing zeros.

    ``20.0 -> "20"``, ``1.5 -> "1.5"`` (docs/algorithm.md §9's naming rules).
    """
    if x == int(x):
        return str(int(x))
    return repr(float(x)).rstrip("0").rstrip(".")


def part_name(input_path: str, cc: float, t: float) -> str:
    """STEP part name ``<input_stem>+cc<cc>+t<t>`` (specification.md §5)."""
    import os

    stem = os.path.splitext(os.path.basename(input_path))[0]
    return f"{stem}+cc{format_param(cc)}+t{format_param(t)}"
