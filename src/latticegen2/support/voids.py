"""Pockets: free space the supports seal off, filled (docs/algorithm.md §14.3 S7).

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
body, the lattice or a support, joined to the supports that seal it. Having no
free face, it is self-supporting. With the fills in place no piece is cut off
any more (gate G5).

The envelope is searched in **slabs**, one per worker of the run's pool
(:func:`find_pockets`): on the ball, three cuts by all 777 junctions in the
envelope were 132 of the search's 157 s.
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
of the region's volume -- or :data:`CUT_PER_TOOL_TOL` per tool, if larger."""

CUT_PER_TOOL_TOL = 0.05
"""mm^3 per tool. Slivers where neighbouring tools touch: junctions meet on
their cap faces, and the volume a cut by many of them removes and the volume
their intersection measures differ by a little per junction. Measured on the
ball's middle slab: 2.7 mm^3 over 263 junctions, 0.01 each, where the
fraction alone allowed 0.43 -- against which a dropped tool, the case this
check exists for, is a whole support node or junction (45 mm^3 on the cup
fixture)."""

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
    bar = max(CUT_VOLUME_REL_TOL * total, CUT_PER_TOOL_TOL * len(tools))
    for fuzzy in CUT_FUZZY_LADDER:
        pieces = []
        for s in region:
            got = _boolean(BRepAlgoAPI_Cut, s, tools, fuzzy)
            if got is None:
                break
            pieces += occ.solids(got)
        else:
            left = sum(occ.volume(p) for p in pieces)
            if abs((total - left) - removes) <= bar:
                return pieces
    got = _free_one_by_one(region, tools, total)
    if got is not None and abs((total - sum(occ.volume(p) for p in got)) - removes) <= bar:
        return got
    raise ProcessingError(
        f"Cutting {what} out of the support envelope did not remove their volume "
        f"({removes:.4f} mm^3 of {total:.4f}) at any fuzzy value "
        f"(docs/algorithm.md §14.3 S7)."
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


SLAB_MIN_JUNCTIONS = 60
"""Lattice junctions per slab below which the envelope is not divided further:
a slab's cuts cost a fixed overhead whatever they cut, so a slab is only worth
a worker when it has this much lattice to cut."""

SLAB_NUDGE = 0.0137
"""mm. Interior slab planes are moved this far off the even division, so that
none lies on a face of the input by accident -- a plane coinciding with a
planar face makes the box-and-body intersection degenerate."""

PLANE_TOL = 1e-4
"""mm. A face whose box lies within this of a slab plane is on it."""

CONTACT_AREA = 1e-8
"""mm^2. Two pieces in neighbouring slabs whose faces on the plane between
them overlap by more than this are one piece of free space."""


def _slab_count(workers: int, junctions: int) -> int:
    """Slabs to cut the envelope into: one per worker, while each still has
    :data:`SLAB_MIN_JUNCTIONS` of lattice to cut."""
    return max(1, min(workers, junctions // SLAB_MIN_JUNCTIONS))


@dataclass
class _Slab:
    """The free space in one slab of the envelope, before and after the tools."""

    before: list
    after: list
    home: list
    """Per ``after`` piece, the index of the ``before`` piece holding it."""
    before_env: list
    after_env: list


def _touches(shape, lo, hi) -> bool:
    slo, shi = (np.asarray(x, dtype=float) for x in occ.bounding_box(shape))
    return bool(np.all(slo <= hi + 1e-6) and np.all(lo <= shi + 1e-6))


def _slab_free(body, probe, lo, hi, env_lo, env_hi, lattice: list, tools: list) -> _Slab:
    """The free space of the body within the box ``[lo, hi]``, before and after
    ``tools``: the region cut by the lattice, then those pieces cut by the
    tools, every cut checked to remove what it should. Each piece after is
    placed by a witness point in the piece before it came from."""
    region = occ.solids(_boolean(_box(lo, hi), body))
    if not region:
        return _Slab([], [], [], [], [])
    lattice_in = (sum(occ.volume(s) for s in occ.solids(_boolean(_compound(region), *lattice)))
                  if lattice else 0.0)
    before = [q for q in _free(region, lattice, lattice_in, "the lattice")
              if occ.volume(q) > MIN_POCKET_VOLUME]
    # Tools lie in the body and overlap neither the lattice nor each other, so
    # what the second cut must remove is their volume inside this box -- one
    # tool per call: fills share faces with the supports they seal, and one
    # call over touching operands partitions them (docs/algorithm.md §7).
    box = _box(lo, hi)
    tools_in = sum(occ.volume(s) for tool in tools for s in occ.solids(_boolean(tool, box)))
    after, home = [], []
    homes = _Inside(before)
    for p in _free(before, tools, tools_in, "the supports"):
        if occ.volume(p) <= MIN_POCKET_VOLUME:
            continue
        w = witness(p)
        if w is None:
            raise ProcessingError(
                f"No point inside a {occ.volume(p):.4g} mm^3 piece of free space "
                f"near {np.round(occ.centroid(p), 3).tolist()} could be found from "
                f"its own tessellation (docs/algorithm.md §14.3 S7)."
            )
        if probe.worst_outside(w[None]) > 0.0:
            raise ProcessingError(
                f"A piece of free space cut from the body near "
                f"{np.round(w, 3).tolist()} lies outside the body by its "
                f"tessellation: the cut and the witness disagree "
                f"(docs/algorithm.md §14.3 S7)."
            )
        k = homes.which(w)
        if not k:
            raise ProcessingError(
                f"A {occ.volume(p):.4g} mm^3 piece of free space near "
                f"{np.round(w, 3).tolist()} lies in no piece of the free space "
                f"before the supports: the cuts and the witness disagree "
                f"(docs/algorithm.md §14.3 S7)."
            )
        after.append(p)
        home.append(k[0])
    return _Slab(before, after, home,
                 [_on_envelope(q, env_lo, env_hi) for q in before],
                 [_on_envelope(p, env_lo, env_hi) for p in after])


def _slab_worker(job):
    """One slab in a worker process: everything is read from the paths and
    plain arrays the job carries, and the pieces go back as ``.brep`` files."""
    from ..classify import OutsideProbe, load_mesh
    from ..lattice import lattice_params
    from ..parallel import compound_children, read_brep, write_brep
    from ..runlog import peak_rss_bytes
    from .assemble import _LocalLattice

    (body_path, mesh_path, lp_args, kept, tools_path, lo, hi, env_lo, env_hi, out) = job
    body = read_brep(body_path)
    mesh = load_mesh(mesh_path)
    probe = OutsideProbe(mesh, 2.0 * mesh.deviation)
    local = _LocalLattice(lattice_params(*lp_args), kept)
    lattice = list(local.near((lo, hi))) if len(kept) else []
    tools = [s for s in compound_children(read_brep(tools_path)) if _touches(s, lo, hi)]
    got = _slab_free(body, probe, lo, hi, env_lo, env_hi, lattice, tools)
    write_brep(_compound(got.before), out + "-before.brep")
    write_brep(_compound(got.after), out + "-after.brep")
    return (len(got.before), len(got.after), got.home, got.before_env, got.after_env,
            out, peak_rss_bytes())


def _faces_on(piece, ax: int, c: float) -> list:
    out = []
    for f in occ.faces(piece):
        flo, fhi = occ.bounding_box(f)
        if abs(flo[ax] - c) <= PLANE_TOL and abs(fhi[ax] - c) <= PLANE_TOL:
            out.append(f)
    return out


def _joined_across(left: list, right: list, ax: int, c: float) -> list:
    """``(i, j)`` pairs of pieces either side of the plane ``x[ax] = c`` whose
    faces on it overlap: one piece of free space cut in two by the slabs."""
    from .build import _boolean as run

    lf = [_faces_on(p, ax, c) for p in left]
    rf = [_faces_on(p, ax, c) for p in right]
    pairs = []
    for i, a in enumerate(lf):
        for j, b in enumerate(rf):
            if not a or not b:
                continue
            got = run(BRepAlgoAPI_Common, _compound(a), [_compound(b)])
            if got is not None and occ.area(got) > CONTACT_AREA:
                pairs.append((i, j))
    return pairs


def _components(n: int, pairs) -> list:
    parent = list(range(n))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j in pairs:
        parent[root(j)] = root(i)
    return [root(i) for i in range(n)]


def find_pockets(body, body_mesh, cells: list, lattice_near, thickness: float,
                 b: np.ndarray, *, pool=None, staged=None) -> Pockets:
    """S7 for a group of support ``cells``, which overlap nothing: the fill of
    every pocket.

    ``lattice_near(box)`` returns the lattice junctions near a box (an empty
    list with no lattice). ``body_mesh`` is the classification mesh, whose
    measured deviation sizes the witness tolerance.

    **The envelope is cut into slabs** across its longest extent, one per
    worker of ``pool`` while each has :data:`SLAB_MIN_JUNCTIONS` of lattice to
    cut, and every slab's free space is found on its own -- across the pool
    when there is one, in this process otherwise, by the same code. Pieces of
    neighbouring slabs whose faces overlap on the plane between them are one
    piece of free space, so pockets are judged on the whole envelope exactly
    as before. ``staged`` -- ``{"body_path", "mesh_path", "lp_args",
    "kept_nodes", "tmpdir"}`` -- is what the workers read.
    """
    from ..classify import OutsideProbe
    from ..parallel import compound_children, read_brep, write_brep

    if not cells:
        return Pockets(stats={"support_pockets": 0})
    lo, hi = envelope(cells, np.asarray(b, dtype=float), body_mesh, thickness)
    ax = int(np.argmax(hi - lo))
    parallel = pool is not None and pool.active and staged is not None
    n = _slab_count(pool.workers if parallel else 1, len(lattice_near((lo, hi))))
    edges = lo[ax] + (hi[ax] - lo[ax]) * np.arange(n + 1) / n
    edges[1:-1] += SLAB_NUDGE
    boxes = []
    for s in range(n):
        slo, shi = lo.copy(), hi.copy()
        slo[ax], shi[ax] = edges[s], edges[s + 1]
        boxes.append((slo, shi))
    probe = OutsideProbe(body_mesh, 2.0 * body_mesh.deviation)

    def slabs(tools: list) -> list:
        if n == 1 or not parallel:
            return [_slab_free(body, probe, slo, shi, lo, hi, list(lattice_near((slo, shi))),
                               [t for t in tools if _touches(t, slo, shi)])
                    for slo, shi in boxes]
        import os
        import uuid

        tag = os.path.join(staged["tmpdir"], f"pockets-{uuid.uuid4().hex[:8]}")
        write_brep(_compound(tools), tag + "-tools.brep")
        lp_args = staged["lp_args"]
        jobs = [(staged["body_path"], staged["mesh_path"], lp_args, staged["kept_nodes"],
                 tag + "-tools.brep", slo, shi, lo, hi, f"{tag}-{s}")
                for s, (slo, shi) in enumerate(boxes)]
        results, _rss = pool.run(_slab_worker, jobs)
        out = []
        for nb, na, home, benv, aenv, path, _r in results:
            before = compound_children(read_brep(path + "-before.brep")) if nb else []
            after = compound_children(read_brep(path + "-after.brep")) if na else []
            out.append(_Slab(before, after, list(home), list(benv), list(aenv)))
        return out

    def pockets(tools: list) -> list:
        got = slabs(tools)
        # Global indices of every slab's pieces, then the pieces joined across
        # the planes between slabs.
        b0 = np.cumsum([0] + [len(s.before) for s in got])
        a0 = np.cumsum([0] + [len(s.after) for s in got])
        bpairs, apairs = [], []
        for s in range(n - 1):
            c = edges[s + 1]
            bpairs += [(b0[s] + i, b0[s + 1] + j)
                       for i, j in _joined_across(got[s].before, got[s + 1].before, ax, c)]
            apairs += [(a0[s] + i, a0[s + 1] + j)
                       for i, j in _joined_across(got[s].after, got[s + 1].after, ax, c)]
        bcomp = _components(int(b0[-1]), bpairs)
        acomp = _components(int(a0[-1]), apairs)
        before_env = [e for s in got for e in s.before_env]
        after = [(s, i) for s in range(n) for i in range(len(got[s].after))]
        groups = defaultdict(lambda: defaultdict(list))   # before comp -> after comp -> pieces
        for g, (s, i) in enumerate(after):
            groups[bcomp[b0[s] + got[s].home[i]]][acomp[g]].append(g)
        out = []
        for bc, members in groups.items():
            reaches = any(e for k, e in enumerate(before_env) if bcomp[k] == bc)
            comps = list(members.values())
            if reaches:
                cut_off = [m for m in comps
                           if not any(got[after[g][0]].after_env[after[g][1]] for g in m)]
            elif len(comps) > 1:
                vol = [sum(occ.volume(got[after[g][0]].after[after[g][1]]) for g in m)
                       for m in comps]
                largest = int(np.argmax(vol))
                cut_off = [m for k, m in enumerate(comps) if k != largest]
            else:
                cut_off = []
            for m in cut_off:
                out += [got[after[g][0]].after[after[g][1]] for g in m]
        return out

    fills = []
    for p in pockets(list(cells)):
        if not occ.is_valid(p, parallel=False):
            raise ProcessingError(
                f"The fill of a pocket near {np.round(occ.centroid(p), 3).tolist()} is "
                f"not a valid solid (docs/algorithm.md §14.3 S7)."
            )
        fills.append(p)
    if fills:
        # G5: with the fills in place nothing is cut off any more.
        again = pockets(list(cells) + fills)
        if again:
            where = np.round(occ.centroid(again[0]), 3).tolist()
            raise ProcessingError(
                f"A pocket near {where} is still cut off after the pockets were "
                f"filled (docs/algorithm.md §14.3 S7, gate G5)."
            )
    return Pockets(fills, {"support_pockets": len(fills), "support_pocket_slabs": n,
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
            "Intersecting the support envelope failed (docs/algorithm.md §14.3 S7)."
        )
    return got
