"""What the support stage receives and what it must hand back (docs/algorithm.md §14.2).

The support generator is handed the **completed lattice** — every solid the
pipeline has assembled and unified, before ``validate`` — together with the
input body, the run's parameters and the inner overhangs already detected. It
returns the **support material only**. The support stage then merges it with
the lattice into the final product (:mod:`.merge`), and everything after it
(``validate``, the export-truth gate, ``export``) gates that product unchanged.

Both types are plain frozen dataclasses so the generator can be exercised in a
unit test without a pipeline around it. docs/algorithm.md §14.4 says how to
write a generator against them.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np
from OCP.TopoDS import TopoDS_Shape

from ..classify import TriMesh
from ..lattice import LatticeParams
from ..parallel import WorkerPool
from .overhang import OverhangSamples


@dataclass(frozen=True)
class SupportInput:
    """Everything the support generator may read. Nothing in it is to be mutated."""

    lattice: tuple[TopoDS_Shape, ...]
    """The lattice solids as the pipeline built them: assembled, closed and
    unified (``simplify``), **not yet validated**. One per connected component
    that survived the floating-body rule (specification.md §5), largest
    component not necessarily first. In the input's own coordinate system."""

    body: TopoDS_Shape
    """The input solid, as read from the STEP file. It still carries the
    classification triangulation (``classify.tessellate_surface``), which
    :func:`.overhang.detect` refines from."""

    body_path: str
    """``body`` staged as a ``.brep`` in ``tmpdir``, for worker processes to
    read instead of having the shape pickled to them (docs/algorithm.md §7)."""

    mesh: TriMesh
    """The classification mesh of ``body``, with its **measured** chordal
    deviation ``mesh.deviation`` (docs/algorithm.md §5.1). What kernel-free
    inside/outside questions should be asked of, e.g. ``classify.OutsideProbe``
    with a margin of ``2 * mesh.deviation`` (§7.2)."""

    mesh_path: str
    """``mesh`` staged as an ``.npz`` in ``tmpdir``, for worker processes."""

    lp: LatticeParams
    """The lattice geometry: ``cc``, ``t``, ``a``, ``r``, the strut frame
    ``e``/``u``/``v``, the basis ``B``, the print orientation ``orient`` and the
    build direction ``b`` in part coordinates (docs/algorithm.md §14.0)."""

    overhang_deg: float
    """``--overhang``: maximum overhang angle from the build direction, degrees."""

    thickness_mm: float
    """``--support-thickness``: minimum thickness of a free-standing member, mm."""

    samples: OverhangSamples
    """The inner overhangs of ``body`` at this orientation and limit
    (docs/algorithm.md §14.1). Never empty: a run with nothing to support does
    not call the generator. ``samples.required`` marks the ones that exceed
    the limit itself, which are what coverage must be proven for."""

    kept_nodes: np.ndarray
    """``(N, 3)`` int64 lattice indices of every node whose junction is part of
    ``lattice``, sorted. ``lattice.nodes(lp, kept_nodes)`` gives their
    positions; the struts between them are known analytically from there,
    without reading them back off the B-rep."""

    pool: WorkerPool
    """The run's shared worker pool (``pool.active`` is False when the run has
    one core). Workers already run at background priority."""

    workers: int
    """``--cores`` as resolved: the most processes this run may occupy."""

    tmpdir: str
    """The run's temporary folder, kept on failure for analysis
    (specification.md §4.4). Intermediate files belong here."""

    report: Optional[Callable[[str, int, Optional[int]], None]] = None
    """Progress callback, ``report(label, done, total)``; ``total=None`` for
    work with no countable fraction. Drives the window's second bar."""


@dataclass(frozen=True)
class SupportResult:
    """What the support generator hands back to the pipeline."""

    supports: list
    """The support solids, and nothing else: no lattice. Each must be a valid
    closed solid inside the input body. They may overlap the lattice and each
    other — the support stage fuses every support with the lattice solids and
    supports it touches (:func:`.merge.merge_with_lattice`), unifies the
    result, and checks that no lattice was lost. A support touching nothing
    stays a body of its own, attached to the wall. Empty means no support was
    generated."""

    stats: dict = field(default_factory=dict)
    """Run characteristics for the end-of-run summary (specification.md §3).
    Keys should start with ``support_``."""

    notes: list = field(default_factory=list)
    """Lines printed to the console and the log regardless of ``-v``: things a
    run should say out loud, such as material left out."""

    lines: list = field(default_factory=list)
    """Lines for the log, and the console only with ``-v``."""
