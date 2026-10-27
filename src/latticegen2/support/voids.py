"""Pockets: free space the supports seal off, filled (docs/algorithm.md §14.3 S5).

Powder must not be trapped. A **pocket** is a region of free space -- inside the
body, outside the lattice and the supports -- that the supports have cut off
from the free space it used to belong to. Every new pocket is bounded partly by
support, and below the support that seals it, so it lies in the **support
envelope**: one box per region, the cells' bounds inflated by ``t_s`` and swept
down along ``-b`` to the body's lowest point. Inflated alone, the box reaches
``t_s`` below the cells, and a cup 5 mm deep under a lid would touch its floor
and never read as sealed.

The free space is found by cutting, not by splitting: ``envelope ∩ body`` minus
the lattice is the free space before the supports, and minus the cells as well
is the free space after them, each a list of solids. The cells are disjoint by
then (S6 runs first), so the second cut must remove exactly their volume; a
split by the cells, tried first, lost volume where neighbouring cells share
faces. Every free piece after the supports is placed by a **witness** point --
asked of tessellations, not of the kernel (§7.2) -- in the piece of free space
it came from, and a piece is a pocket when it was cut off:

* where that free space reaches the envelope's own faces, every piece that no
  longer does is cut off;
* where the envelope holds all of it -- a small body, whose free space the body
  alone already encloses -- the largest piece is that free space and every other
  piece is cut off from it.

The fill is the pocket itself: an exact solid each of whose faces lies on the
body, the lattice or a support, shipped as a support body of its own. Having no
free face, it is self-supporting. With the fills in place no piece is cut off
any more (gate G5).
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.gp import gp_Pnt

from .. import occ
from ..errors import ProcessingError

WITNESS_DEFLECTION = 0.01
"""mm. Tessellation the witness tests read pieces through."""

ON_ENVELOPE_TOL = 1e-5
"""mm. A face whose every mesh point lies this close to one of the envelope's
planes is a face of the envelope."""

CUT_FUZZY_LADDER = (0.0, 1e-5, 1e-4, 1e-3)
"""mm. Fuzzy values a cut is tried with until it removes what it should -- the
same ladder a cell's trim uses (``build.FUZZY_LADDER``). Measured on the test
chamber: nine disjoint cells sharing faces with each other, and the cut removes
178.67 of their 205.75 mm^3 at every value up to 1e-4 and all of it at 1e-3."""

CUT_VOLUME_REL_TOL = 1e-5
"""What a cut removes must match what it was asked to remove, to this fraction
of the envelope's volume."""

MIN_POCKET_VOLUME = 1e-3
"""mm^3. Free pieces smaller than this -- a cube a tenth of a millimetre on a
side -- are boolean specks along faces two bodies share, not pockets: on the
cup fixture one of 2e-5 mm^3 had no interior point its own tessellation could
find."""


@dataclass
class Pockets:
    fills: list = field(default_factory=list)
    stats: dict = field(default_factory=dict)


def _solid_mesh(shape, deflection: float = WITNESS_DEFLECTION):
    """A closed, welded ``TriMesh`` of ``shape``'s surface."""
    from ..classify import TriMesh, _weld
    from .build import _mesh_of

    verts, tris, offset = [], [], 0
    for _f, p, _u, t in _mesh_of(shape, deflection):
        verts.append(p)
        tris.append(t + offset)
        offset += len(p)
    if not verts:
        return None
    v, t = _weld(np.vstack(verts), np.vstack(tris))
    return TriMesh(v, t)


class _Inside:
    """Point-in-solid for a list of solids, by ray parity over their meshes."""

    def __init__(self, shapes: list):
        from ..classify import PointInside

        self.tests = []
        for s in shapes:
            m = _solid_mesh(s)
            if m is not None:
                lo, hi = m.bounds()
                self.tests.append((lo, hi, PointInside(m), m))

    def which(self, p: np.ndarray) -> list:
        """Indices of the solids holding ``p``.

        Failing that, the one whose mesh is nearest, within the mesh's own
        deflection: a sliver's deepest point can lie microns inside a curved
        wall the mesh's chords cut inside of -- on the cup fixture, 4 microns
        from a wall meshed to 0.01 mm.
        """
        from ..classify import _point_triangle_dist

        out = []
        for k, (lo, hi, test, _m) in enumerate(self.tests):
            if np.all(p >= lo) and np.all(p <= hi) and bool(test(p[None])[0]):
                out.append(k)
        if out:
            return out
        best, where = WITNESS_DEFLECTION, None
        for k, (lo, hi, _t, m) in enumerate(self.tests):
            if np.any(p < lo - best) or np.any(p > hi + best):
                continue
            a, b, c = m.triangle_points
            d = float(_point_triangle_dist(p, a, b, c).min())
            if d <= best:
                best, where = d, k
        return [] if where is None else [where]


def witness(piece) -> np.ndarray | None:
    """A point well inside ``piece``, by its own tessellation, or ``None``.

    Candidates are the surface centroid and points a little inside every
    triangle; of those the piece's own mesh puts inside, the one furthest from
    its surface wins, so a test made at it is not decided by a mesh's chords.
    """
    from ..classify import PointInside, _point_triangle_dist

    m = _solid_mesh(piece)
    if m is None:
        return None
    a, b, c = m.triangle_points
    n = np.cross(b - a, c - a)
    nl = np.linalg.norm(n, axis=1)
    ok = nl > 1e-14
    a, b, c, n = a[ok], b[ok], c[ok], n[ok] / nl[ok][:, None]
    cen = (a + b + c) / 3.0
    if len(cen) > 600:
        pick = np.linspace(0, len(cen) - 1, 600).astype(int)
        cen, n = cen[pick], n[pick]
    size = float(np.linalg.norm(m.verts.max(axis=0) - m.verts.min(axis=0)))
    cands = [occ.centroid(piece)[None]]
    for eps in (0.1, 0.02, 0.004):
        cands.append(cen - min(eps, size / 20.0) * n)
    cands = np.vstack(cands)
    cands = cands[PointInside(m)(cands)]
    if not len(cands):
        return None
    A, B, C = m.triangle_points
    depth = np.array([_point_triangle_dist(p, A, B, C).min() for p in cands[:200]])
    return cands[int(np.argmax(depth))]


def _box(lo: np.ndarray, hi: np.ndarray):
    return BRepPrimAPI_MakeBox(gp_Pnt(*map(float, lo)), gp_Pnt(*map(float, hi))).Shape()


def _free(region: list, tools: list, removes: float, what: str) -> list:
    """The solids of ``region`` minus ``tools``, checked to remove ``removes`` mm^3."""
    from .build import _boolean

    total = sum(occ.volume(s) for s in region)
    if not tools:
        return list(region)
    for fuzzy in CUT_FUZZY_LADDER:
        pieces = []
        for s in region:
            got = _boolean(BRepAlgoAPI_Cut, s, tools, fuzzy)
            if got is None:
                break
            pieces += occ.solids(got)
        else:
            left = sum(occ.volume(p) for p in pieces)
            if abs((total - left) - removes) <= CUT_VOLUME_REL_TOL * total:
                return pieces
    got = _free_one_by_one(region, tools, total)
    if got is not None and abs((total - sum(occ.volume(p) for p in got)) - removes) \
            <= CUT_VOLUME_REL_TOL * total:
        return got
    raise ProcessingError(
        f"Cutting {what} out of the support envelope did not remove their volume "
        f"({removes:.4f} mm^3 of {total:.4f}) at any fuzzy value "
        f"(docs/algorithm.md §14.3 S5)."
    )


ONE_BY_ONE_PASSES = 3
"""Orders :func:`_free_one_by_one` tries before giving up."""


def _free_one_by_one(region: list, tools: list, total: float):
    """``region`` minus ``tools``, one tool at a time, or ``None``.

    The kernel can take one of two face-sharing neighbours for the inside of
    the other. Measured on the cup fixture once its wall cells reached the
    crease: one cut by all 13 cells, each of which cuts exactly on its own,
    removed 536.05 of 581.50 mm^3 -- short by exactly one cell's 45.44 --
    at every fuzzy value; fusing that cell with its neighbour returned the
    neighbour alone; and cutting the pair one at a time worked in one order
    and not the other.

    So each step is checked against what that tool alone takes out of the
    region, and a tool whose step falls short is moved to the front of the
    order for the next pass.
    """
    from .build import _boolean

    alone = []
    for t in tools:
        got = [_boolean(BRepAlgoAPI_Common, t, [s]) for s in region]
        if any(g is None for g in got):
            return None
        alone.append(sum(occ.volume(x) for g in got for x in occ.solids(g)))
    order = list(range(len(tools)))
    for _pass in range(ONE_BY_ONE_PASSES):
        pieces = list(region)
        short = []
        for k in order:
            before = sum(occ.volume(p) for p in pieces)
            nxt = []
            for s in pieces:
                got = _boolean(BRepAlgoAPI_Cut, s, [tools[k]])
                if got is None:
                    return None
                nxt += occ.solids(got)
            if abs((before - sum(occ.volume(p) for p in nxt)) - alone[k]) \
                    > CUT_VOLUME_REL_TOL * total:
                short.append(k)
                continue
            pieces = nxt
        if not short:
            return pieces
        order = short + [k for k in order if k not in short]
    return None


def _on_envelope(piece, lo: np.ndarray, hi: np.ndarray) -> bool:
    """Whether any face of ``piece`` lies on one of the envelope box's planes."""
    from .build import _mesh_of

    for _f, pts, _u, _t in _mesh_of(piece, 0.1):
        for k in range(3):
            for plane in (lo[k], hi[k]):
                if np.all(np.abs(pts[:, k] - plane) <= ON_ENVELOPE_TOL):
                    return True
    return False


def envelope(cells: list, b: np.ndarray, body_mesh, thickness: float):
    """``(lo, hi)``: the cells' box inflated by ``thickness``, swept down ``-b``
    to the lowest point of the body."""
    lo = np.min([occ.bounding_box(c)[0] for c in cells], axis=0) - thickness
    hi = np.max([occ.bounding_box(c)[1] for c in cells], axis=0) + thickness
    corners = np.array([[x, y, z] for x in (lo[0], hi[0]) for y in (lo[1], hi[1])
                        for z in (lo[2], hi[2])])
    drop = float((corners @ b).min() - (body_mesh.verts @ b).min()) + thickness
    swept = np.vstack([corners, corners - max(drop, 0.0) * b])
    return swept.min(axis=0), swept.max(axis=0)


def find_pockets(body, body_mesh, cells: list, lattice_near, thickness: float,
                 b: np.ndarray) -> Pockets:
    """S5 for one region's ``cells``, which overlap nothing: the fill of every pocket.

    ``lattice_near(box)`` returns the lattice junctions near a box (an empty list
    with no lattice). ``body_mesh`` is the classification mesh, whose measured
    deviation sizes the witness tolerance.
    """
    from ..classify import OutsideProbe

    if not cells:
        return Pockets(stats={"support_pockets": 0})
    lo, hi = envelope(cells, np.asarray(b, dtype=float), body_mesh, thickness)
    region = occ.solids(_boolean(_box(lo, hi), body))
    if not region:
        return Pockets(stats={"support_pockets": 0})
    lattice = list(lattice_near((lo, hi)))
    probe = OutsideProbe(body_mesh, 2.0 * body_mesh.deviation)
    lattice_in = (sum(occ.volume(s) for s in occ.solids(_boolean(_compound(region), *lattice)))
                  if lattice else 0.0)
    before = _free(region, lattice, lattice_in, "the lattice")
    before = [(q, _on_envelope(q, lo, hi)) for q in before
              if occ.volume(q) > MIN_POCKET_VOLUME]
    homes = _Inside([q for q, _t in before])

    def pockets_after(tools: list, removes: float) -> list:
        """``[(piece, witness)]``: the free pieces cut off by ``tools``."""
        groups = defaultdict(list)
        for p in _free(region, tools, removes, "the supports"):
            if occ.volume(p) <= MIN_POCKET_VOLUME:
                continue
            w = witness(p)
            if w is None:
                raise ProcessingError(
                    f"No point inside a {occ.volume(p):.4g} mm^3 piece of free space "
                    f"near {np.round(occ.centroid(p), 3).tolist()} could be found from "
                    f"its own tessellation (docs/algorithm.md §14.3 S5)."
                )
            if probe.worst_outside(w[None]) > 0.0:
                raise ProcessingError(
                    f"A piece of free space cut from the body near "
                    f"{np.round(w, 3).tolist()} lies outside the body by its "
                    f"tessellation: the cut and the witness disagree "
                    f"(docs/algorithm.md §14.3 S5)."
                )
            home = homes.which(w)
            if not home:
                raise ProcessingError(
                    f"A {occ.volume(p):.4g} mm^3 piece of free space near "
                    f"{np.round(w, 3).tolist()} lies in no piece of the free space "
                    f"before the supports: the cuts and the witness disagree "
                    f"(docs/algorithm.md §14.3 S5)."
                )
            groups[home[0]].append((p, w, _on_envelope(p, lo, hi)))
        out = []
        for k, members in groups.items():
            if before[k][1]:
                out += [(p, w) for p, w, touches in members if not touches]
            elif len(members) > 1:
                members = sorted(members, key=lambda m: occ.volume(m[0]), reverse=True)
                out += [(p, w) for p, w, _t in members[1:]]
        return out

    removes = lattice_in + sum(occ.volume(c) for c in cells)
    fills = []
    for p, w in pockets_after(list(cells) + lattice, removes):
        if not occ.is_valid(p, parallel=False):
            raise ProcessingError(
                f"The fill of a pocket near {np.round(w, 3).tolist()} is not a "
                f"valid solid (docs/algorithm.md §14.3 S5)."
            )
        fills.append(p)
    if fills:
        # G5: with the fills in place nothing is cut off any more.
        again = pockets_after(list(cells) + fills + lattice,
                              removes + sum(occ.volume(f) for f in fills))
        if again:
            where = np.round(again[0][1], 3).tolist()
            raise ProcessingError(
                f"A pocket near {where} is still cut off after the pockets were "
                f"filled (docs/algorithm.md §14.3 S5, gate G5)."
            )
    return Pockets(fills, {"support_pockets": len(fills),
                           "support_pocket_volume_mm3":
                               round(sum(occ.volume(f) for f in fills), 4)})


def _compound(shapes: list):
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Compound

    builder = BRep_Builder()
    comp = TopoDS_Compound()
    builder.MakeCompound(comp)
    for s in shapes:
        builder.Add(comp, s)
    return comp


def _boolean(a, *tools):
    from .build import _boolean as run

    got = run(BRepAlgoAPI_Common, a, list(tools))
    if got is None:
        raise ProcessingError(
            "Intersecting the support envelope failed (docs/algorithm.md §14.3 S5)."
        )
    return got
