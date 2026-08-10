"""The junction template ``J`` — the one solid this whole architecture rests on.

Cut every strut at its midpoint and assign each half to its nearer node. Every
node then owns exactly six half-struts whose union is a solid that is
**congruent at every node of the lattice**, differing only by a translation
(docs/research/perf-rearchitecture-proposal.md §2, Insight 3). Because the three
strut directions are mutually orthogonal, all strut-strut overlap is contained
inside that solid: two adjacent junctions never overlap volumetrically, they
meet exactly on the shared mid-strut cap quad.

So the fused lattice is translated copies of one solid glued on coincident
quads, and the single 6-operand fuse in this module is the **only** general
boolean the program ever performs. Everything downstream is instancing and
topological stitching.

The cap-integrity check in :func:`build_template` is what makes that safe: if the
node-region overlap ever swallowed a cap quad, the interfaces would not exist and
the instancing argument would collapse.

It never does, and the reason is worth stating because it removes a restriction
the design was expected to carry. A half-strut along ``e_j`` extends toward an
orthogonal ``e_k`` only as far as its profile's *support* in that direction. That
support is ``r/sqrt(2) = t/2`` — the square's **inradius**, not its circumradius —
because each diamond profile presents an edge, never a corner, to every
orthogonal strut direction (verified numerically for all six ``(j, k)`` pairs in
``test/test_junction.py``). Caps therefore stay intact while ``t/2 < a/2``, i.e.
for exactly ``t < a``, which is the CLI's existing cross-constraint. There is no
narrower parameter window for the fast path: a swept ``t/cc`` from 0.45 up to the
``t < a`` limit at ``cc`` of 5, 10 and 20 mm keeps all six caps throughout.

The check is kept regardless, as a cheap (~40 ms) geometric verification at the
run's actual parameters, so a future change to the profile or the cell cannot
silently invalidate the argument.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from OCP.BRep import BRep_Builder
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCP.GeomAbs import GeomAbs_SurfaceType
from OCP.TopoDS import TopoDS_Face, TopoDS_Shape, TopoDS_Shell
from OCP.TopTools import TopTools_ListOfShape

from . import occ
from .errors import ParamError, ProcessingError
from .lattice import HALF_STRUTS, LatticeParams, half_strut_offset, profile_vertices


@dataclass
class JunctionTemplate:
    """One node's fused six half-struts, plus everything instancing needs."""

    solid: TopoDS_Shape
    """The fused junction ``J``, centred on the origin."""
    open_shell: TopoDS_Shell
    """``J``'s boundary **minus all six cap quads** — the interior instance.

    Every cap of an INTERIOR node is always dropped (its neighbour is always a
    kept node — see :mod:`latticegen2.classify`), so this one shell serves every
    interior junction in the run and is instanced with an O(1) ``Moved()``.
    """
    cap_faces: list[TopoDS_Face]
    """The six cap quads, indexed by half-strut id (``lattice.HALF_STRUTS``)."""
    cap_centers: np.ndarray
    """``(6, 3)`` local cap-quad centres, indexed by half-strut id."""
    volume: float
    """Exact volume of ``J`` in mm³ — the basis of the interior volume identity."""
    n_faces: int
    cap_area_error: float
    """Worst relative deviation of a cap's area from the ideal ``t²``."""


def _half_strut_solid(lp: LatticeParams, h: int) -> TopoDS_Shape:
    """One half-strut prism at the origin: profile polygon extruded by ``a/2``."""
    k, _ = HALF_STRUTS[h]
    base = occ.polygon_face(profile_vertices(lp, np.zeros(3), k))
    return occ.prism(base, half_strut_offset(lp, h))


def _plane_of(face: TopoDS_Face) -> tuple[np.ndarray, np.ndarray] | None:
    """``(normal, point)`` of a planar face, or ``None`` if it is not planar."""
    surf = BRepAdaptor_Surface(face)
    if surf.GetType() != GeomAbs_SurfaceType.GeomAbs_Plane:
        return None
    pln = surf.Plane()
    d = pln.Axis().Direction()
    p = pln.Location()
    return np.array([d.X(), d.Y(), d.Z()]), np.array([p.X(), p.Y(), p.Z()])


def _identify_caps(lp: LatticeParams, shape: TopoDS_Shape, tol: float):
    """Find the six cap quads of a junction centred at the origin.

    A cap for half-strut ``h`` is the planar face whose normal is parallel to
    ``e_k``, whose plane sits at signed distance ``s * a/2`` from the node along
    ``e_k``, and whose area is ``t²``. Lateral faces cannot be mistaken for one:
    a lateral plane of any half-strut lies at distance ``r = t/sqrt(2)`` from the
    node, and the ``t < cc/2`` regime this template requires puts ``r`` strictly
    below the caps' ``a/2``.
    """
    found: dict[int, TopoDS_Face] = {}
    worst_area_err = 0.0
    ideal = lp.t * lp.t
    for face in occ.faces(shape):
        plane = _plane_of(face)
        if plane is None:
            continue
        normal, point = plane
        for h, (k, s) in enumerate(HALF_STRUTS):
            if abs(abs(float(np.dot(normal, lp.e[k]))) - 1.0) > 1e-9:
                continue
            if abs(float(np.dot(point, lp.e[k])) - s * lp.a / 2.0) > tol:
                continue
            a_face = occ.area(face)
            if abs(a_face - ideal) > 1e-6 * max(ideal, 1.0):
                # Right plane, wrong size: the cap has been eaten into by the
                # node-region overlap. Report it as missing rather than accept a
                # partial interface.
                continue
            worst_area_err = max(worst_area_err, abs(a_face - ideal) / ideal)
            found[h] = face
            break
    return found, worst_area_err


def build_template(lp: LatticeParams, tol: float = 1e-7) -> JunctionTemplate:
    """Build ``J`` and verify its six cap quads survived the fuse intact.

    Raises :class:`ParamError` if a cap is missing. No accepted parameter
    combination is expected to trigger it (see the module docstring), so this is
    a defensive gate rather than a documented limit — it runs in about 40 ms,
    before any real work, so a violated assumption surfaces immediately instead
    of as broken geometry.
    """
    halves = [_half_strut_solid(lp, h) for h in range(6)]

    fuse = BRepAlgoAPI_Fuse()
    args = TopTools_ListOfShape()
    args.Append(halves[0])
    tools = TopTools_ListOfShape()
    for s in halves[1:]:
        tools.Append(s)
    fuse.SetArguments(args)
    fuse.SetTools(tools)
    fuse.Build()
    if not fuse.IsDone():
        raise ProcessingError(
            "Could not fuse the six half-struts of the junction template. This is "
            "the only general boolean the pipeline performs; its failure means the "
            "requested profile/cell combination produces degenerate geometry."
        )
    fuse.SimplifyResult()
    solid = fuse.Shape()

    caps, area_err = _identify_caps(lp, solid, tol)
    missing = [h for h in range(6) if h not in caps]
    if missing:
        raise ParamError(
            f"Junction cap integrity check failed for cc={lp.cc} t={lp.t}: "
            f"{len(missing)} of 6 mid-strut cap quads did not survive the template "
            f"fuse as intact {lp.t}x{lp.t} faces (half-strut ids {missing}). Those "
            f"faces are the interfaces this generator joins junctions along, so it "
            f"cannot build a watertight lattice at these parameters. This is not "
            f"expected for any -cc/-t inside the documented ranges — please report it "
            f"with these values. A smaller -t relative to -cc will work around it."
        )

    cap_faces = [caps[h] for h in range(6)]
    cap_centers = np.array(
        [half_strut_offset(lp, h) for h in range(6)], dtype=float
    )

    all_faces = occ.faces(solid)
    keep = [f for f in all_faces if not any(f.IsSame(c) for c in cap_faces)]
    if len(keep) != len(all_faces) - 6:
        raise ProcessingError(
            "Junction template face bookkeeping is inconsistent: expected to drop "
            f"exactly 6 cap faces from {len(all_faces)}, dropped {len(all_faces) - len(keep)}."
        )

    shell = TopoDS_Shell()
    builder = BRep_Builder()
    builder.MakeShell(shell)
    for f in keep:
        builder.Add(shell, f)

    return JunctionTemplate(
        solid=solid,
        open_shell=shell,
        cap_faces=cap_faces,
        cap_centers=cap_centers,
        volume=occ.volume(solid),
        n_faces=len(all_faces),
        cap_area_error=area_err,
    )


def is_cap_plane_face(lp: LatticeParams, face: TopoDS_Face, node_pos: np.ndarray,
                      tol: float = 1e-6) -> int | None:
    """Which cap plane (half-strut id) a face lies in, or ``None``.

    Used on **trimmed** boundary junctions, where a cap may have been cut down
    to a fragment or removed entirely, so the area test of :func:`_identify_caps`
    does not apply — only the plane does. The plane alone is unambiguous for the
    same reason: lateral planes sit at ``r`` from the node, caps at ``a/2``.
    """
    plane = _plane_of(face)
    if plane is None:
        return None
    normal, point = plane
    rel = point - np.asarray(node_pos, dtype=float)
    for h, (k, s) in enumerate(HALF_STRUTS):
        if abs(abs(float(np.dot(normal, lp.e[k]))) - 1.0) > 1e-9:
            continue
        if abs(float(np.dot(rel, lp.e[k])) - s * lp.a / 2.0) > tol:
            continue
        return h
    return None
