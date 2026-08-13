"""Gate G5 — what makes ``BRepBuilderAPI_Sewing`` expensive at boundary scale.

The scale rehearsal (``TD_HX_Indre_Volum``, ``cc=5 t=1``) spent **4 h 45 m** of a
5 h 04 m run in the ``sew`` stage, single-threaded, stitching 21,697 shells.
docs/specification.md §10 predicted this and sketched two ways out, but they cost
very different amounts to build and the choice turns on *which* term dominates:

  * the **free edges** sewing has to pair up by geometric search — if that is the
    cost, the fix is to stop searching for a pairing the junction graph already
    knows (index welding), or at least to partition the search;
  * the sheer **face count** it indexes on the way in — if that is the cost,
    partitioning helps and welding the interfaces does not.

This gate separates them, on the real pieces the failed run left behind rather
than on a synthetic stand-in:

  1. a scaling ladder over piece count, to get the exponent;
  2. the same ladder with OCCT's ``Cutting`` (free-edge splitting),
     ``Analysis`` and ``SameParameter`` phases turned off — all three are work
     this architecture does not need, since interfaces match by construction;
  3. the delta from adding a large **closed** instanced interior shell, which
     contributes faces and *no* free edges, isolating the face-count term.

Run against the kept temp folder of a failed or cancelled run::

    python tools/prototypes/g5_stitch_scaling.py --pieces path/to/temp/<stamp>

There is no PASS bar here. This is a measurement that decides a design, and it
prints the numbers and the fitted exponent for docs/prototypes/RESULTS.md.
"""

import argparse
import glob
import os
import re
import sys
import time

import _bootstrap  # noqa: F401
import numpy as np
from OCP.BRep import BRep_Builder
from OCP.BRepBuilderAPI import BRepBuilderAPI_Sewing
from OCP.BRepTools import BRepTools
from OCP.TopAbs import TopAbs_ShapeEnum
from OCP.TopoDS import TopoDS_Iterator, TopoDS_Shape, TopoDS_Shell

from latticegen2 import occ
from latticegen2.connect import lattice_interfaces
from latticegen2.interior import build_interior_shell, extract_template_mesh
from latticegen2.junction import build_template
from latticegen2.lattice import lattice_params

SEW_TOL = 1e-6

DEFAULT_PIECES = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "..", "..", "gmsh-occt-memory-write-69e632", "temp", "20260812-192100",
)
"""Where the 5-hour failure left its 21,955 trimmed boundary pieces.

A specific run's leavings, so it will not be there forever — pass ``--pieces``
to point at any kept ``temp/<stamp>`` folder. Nothing else in the repository
depends on it."""

VARIANTS = {
    # (cutting, analysis, same_parameter)
    "default": (True, True, True),
    "nocut": (False, True, True),
    "nocut+nosameparam": (False, True, False),
    "minimal": (False, False, False),
}
"""``default`` is what ``latticegen2.occ.sew`` uses today. The rest switch off
phases this architecture cannot benefit from: ``Cutting`` splits free edges so
they match, but ours already match exactly; ``Analysis`` looks for degenerate
shapes that instanced prisms and their trims do not produce; ``SameParameter``
re-fits merged edges."""


# --- loading the pieces -----------------------------------------------------


def _children(shape):
    out = []
    it = TopoDS_Iterator(shape)
    while it.More():
        out.append(it.Value())
        it.Next()
    return out


def _as_shell(shape: TopoDS_Shape) -> TopoDS_Shape:
    """Accept either layout a run's ``boundary_*.brep`` can be in.

    Before the symmetric-interface rework a worker wrote one already-assembled
    shell per piece; it now writes a compound of that piece's faces, because
    which of them are dropped is decided later on the master. Either way what
    sewing wants is a shell.
    """
    if shape.ShapeType() == TopAbs_ShapeEnum.TopAbs_SHELL:
        return shape
    shell = TopoDS_Shell()
    builder = BRep_Builder()
    builder.MakeShell(shell)
    for face in occ.faces(shape):
        builder.Add(shell, face)
    return shell


def load_pieces(directory: str, limit: int):
    """Up to ``limit`` trimmed boundary pieces, in the order the run produced them."""
    paths = sorted(
        glob.glob(os.path.join(directory, "boundary_*.brep")),
        key=lambda p: int(re.search(r"boundary_(\d+)", p).group(1)),
    )
    if not paths:
        raise SystemExit(
            f"No boundary_*.brep in {directory}.\n"
            f"Point --pieces at a kept temp/<stamp> folder from a failed or "
            f"cancelled run (the folder is left in place for exactly this)."
        )
    out = []
    for path in paths:
        shape = TopoDS_Shape()
        if not BRepTools.Read_s(shape, path, BRep_Builder()):
            raise SystemExit(f"Could not read {path}")
        for child in _children(shape):
            out.append(_as_shell(child))
            if len(out) >= limit:
                return out
    return out


# --- the measurement --------------------------------------------------------


def sew_once(shapes, cutting: bool, analysis: bool, same_parameter: bool):
    """One sewing run, returning wall time and the diagnostics OCCT offers."""
    sewing = BRepBuilderAPI_Sewing(SEW_TOL, True, analysis, cutting, False)
    sewing.SetSameParameterMode(same_parameter)
    for s in shapes:
        sewing.Add(s)
    t0 = time.perf_counter()
    sewing.Perform()
    elapsed = time.perf_counter() - t0
    result = sewing.SewedShape()
    shells = occ.shells(result)
    return {
        "seconds": elapsed,
        "shells": len(shells),
        "open_shells": sum(1 for s in shells if not s.Closed()),
        "free_edges": sewing.NbFreeEdges(),
        "multiple_edges": sewing.NbMultipleEdges(),
        "contiguous_edges": sewing.NbContigousEdges(),
    }


def fit_exponent(sizes, times):
    """Least-squares slope of log(time) against log(size) — the scaling exponent."""
    usable = [(n, t) for n, t in zip(sizes, times) if t > 0]
    if len(usable) < 2:
        return float("nan")
    x = np.log(np.array([n for n, _ in usable], dtype=float))
    y = np.log(np.array([t for _, t in usable], dtype=float))
    return float(np.polyfit(x, y, 1)[0])


def ladder(pieces, counts, variant, budget):
    """Sew growing prefixes until one exceeds ``budget`` seconds."""
    cutting, analysis, same_parameter = VARIANTS[variant]
    sizes, times = [], []
    for n in counts:
        if n > len(pieces):
            break
        got = sew_once(pieces[:n], cutting, analysis, same_parameter)
        sizes.append(n)
        times.append(got["seconds"])
        print(
            f"  {variant:<18} {n:>6} pieces  {got['seconds']:8.2f} s  "
            f"{got['shells']:>5} shell(s) ({got['open_shells']} open)  "
            f"free {got['free_edges']:>6}  multiple {got['multiple_edges']:>5}  "
            f"contiguous {got['contiguous_edges']:>6}",
            flush=True,
        )
        if got["seconds"] > budget:
            print(f"  {variant:<18} stopping: past the {budget:g} s per-measurement budget")
            break
    if len(sizes) >= 2:
        print(f"  {variant:<18} fitted exponent {fit_exponent(sizes, times):.2f}"
              f"  ({sizes[0]} -> {sizes[-1]} pieces)")
    return sizes, times


def face_count_delta(pieces, n, m, budget):
    """Cost of adding a large closed shell: faces, and no free edges at all.

    A closed instanced grid contributes nothing for sewing to pair up, so any
    time it adds is the price of indexing faces rather than of the geometric
    search. That is exactly the term that decides whether partitioning the sew
    would help.
    """
    lp = lattice_params(10.0, 1.5)
    tpl = build_template(lp)
    tmesh = extract_template_mesh(lp, tpl)
    a = np.arange(m, dtype=np.int64)
    g = np.meshgrid(a, a, a, indexing="ij")
    ns = np.stack([x.ravel() for x in g], axis=1)
    kept = {(int(x), int(y), int(z)) for x, y, z in ns}
    t0 = time.perf_counter()
    shell, stats = build_interior_shell(lp, tpl, tmesh, ns, lattice_interfaces(kept))
    print(f"  interior shell: {stats['interior_faces']} faces, "
          f"{stats['interior_open_edges']} free edges, built in "
          f"{time.perf_counter() - t0:.1f} s")
    if not shell.Closed():
        print("  WARNING: the reference interior shell is not closed; delta is not clean")

    cutting, analysis, same_parameter = VARIANTS["nocut"]
    without = sew_once(pieces[:n], cutting, analysis, same_parameter)
    with_shell = sew_once(list(pieces[:n]) + [shell], cutting, analysis, same_parameter)
    print(f"  {n} pieces alone            {without['seconds']:8.2f} s")
    print(f"  {n} pieces + {stats['interior_faces']} faces "
          f"{with_shell['seconds']:8.2f} s  "
          f"(+{with_shell['seconds'] - without['seconds']:.2f} s)")
    return without, with_shell, stats


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pieces", default=DEFAULT_PIECES,
                    help="kept temp/<stamp> folder holding boundary_*.brep")
    ap.add_argument("--budget", type=float, default=600.0,
                    help="seconds; stop a ladder once one measurement exceeds it")
    ap.add_argument("--max", type=int, default=32000, help="most pieces to load")
    ap.add_argument("--grid", type=int, default=20,
                    help="edge of the reference interior grid (m^3 junctions)")
    ap.add_argument("--variants", default=",".join(VARIANTS),
                    help=f"comma-separated subset of {', '.join(VARIANTS)}")
    args = ap.parse_args()

    chosen = [v.strip() for v in args.variants.split(",") if v.strip()]
    unknown = [v for v in chosen if v not in VARIANTS]
    if unknown:
        raise SystemExit(f"Unknown variant(s) {unknown}; choose from {list(VARIANTS)}")

    directory = os.path.abspath(args.pieces)
    print(f"G5: reading trimmed boundary pieces from {directory}")
    t0 = time.perf_counter()
    pieces = load_pieces(directory, args.max)
    faces = sum(len(occ.faces(p)) for p in pieces)
    print(f"G5: {len(pieces)} pieces, {faces} faces, loaded in "
          f"{time.perf_counter() - t0:.1f} s "
          f"({faces / max(len(pieces), 1):.1f} faces/piece)\n")

    counts = [n for n in (500, 1000, 2000, 4000, 8000, 16000, 32000) if n <= len(pieces)]
    if len(pieces) not in counts:
        counts.append(len(pieces))

    print("G5a — scaling of the current configuration and of cheaper ones")
    results = {}
    for variant in chosen:
        results[variant] = ladder(pieces, counts, variant, args.budget)
        print()

    print("G5b — what the face count alone costs")
    reference = min(4000, len(pieces))
    face_count_delta(pieces, reference, args.grid, args.budget)
    return 0


if __name__ == "__main__":
    sys.exit(main())
