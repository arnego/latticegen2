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
any more (gate G5) -- by construction, since each fill is a whole cut-off piece.

The envelope is searched in **slabs** across the run's pool
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
from ..runlog import box_text


class CutRefused(ProcessingError):
    """The kernel would not cut a slab of the support envelope as checked: a
    slab raising this is split and searched again, and left open at last."""

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
    unsearched: list = field(default_factory=list)
    """One line per region the kernel would not cut, with its coordinates."""
    refusals: list = field(default_factory=list)
    """The kernel's reason for each of them, for the log."""


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


def _free(region: list, tools: list, removes: float, what: str,
          one_by_one: bool = True) -> list:
    """The solids of ``region`` minus ``tools``, checked to remove ``removes`` mm^3.

    ``one_by_one=False`` skips the tool-at-a-time fallback, for a caller that
    has a cheaper second route to try first."""
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
    got = _free_one_by_one(region, tools, total) if one_by_one else None
    if got is not None and abs((total - sum(occ.volume(p) for p in got)) - removes) <= bar:
        return got
    raise CutRefused(
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


SLAB_JUNCTIONS = 150
"""Lattice junctions per slab the envelope is divided for. A cut's cost grows
faster than the number of tools it is given, so the slabs follow the lattice,
not the workers: on the rehearsal part six slabs of ~3,000 junctions each were
still cutting after half an hour, where the ball's slabs of ~150 took seconds."""

SLAB_REACH_MULTIPLE = 4.0
"""A slab is at least this many junction reaches thick: a junction is cut by
every slab it reaches into, so thinner slabs share most of their tools."""

SLAB_NUDGE = 0.0137
"""mm. Interior slab planes are moved this far off the even division, so that
none lies on a face of the input by accident -- a plane coinciding with a
planar face makes the box-and-body intersection degenerate."""

PLANE_TOL = 1e-4
"""mm. A face whose box lies within this of a slab plane is on it."""

SLAB_RETRY_SPLITS = 2
"""Times a slab whose cuts the kernel refuses is split in two and searched
again. Measured on the rehearsal part: one slab of 106 was refused by both
routes, the combined cut and the supports alone, at every fuzzy value."""

UNSEARCHED_WARNING = "enclosed pockets may remain unfilled"
"""Title of the warning naming the regions left unsearched."""

def _slab_count(junctions: int, extent: float, reach: float) -> int:
    """Slabs to cut an envelope ``extent`` long into: one per
    :data:`SLAB_JUNCTIONS` of lattice in it, none thinner than
    :data:`SLAB_REACH_MULTIPLE` junction reaches. The same with a pool and
    without, so the free space found does not depend on how it was run."""
    by_lattice = -(-junctions // SLAB_JUNCTIONS)
    by_size = int(extent // (SLAB_REACH_MULTIPLE * reach)) if reach > 0.0 else by_lattice
    return max(1, min(by_lattice, by_size))


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
    ``tools``: the region cut by the lattice, and by the lattice and the tools
    together, every cut checked to remove what it should. Each piece after is
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
    # The region cut by lattice and supports in one call first. Neither route
    # is always taken by the kernel, measured on the rehearsal part: the free
    # space above cut by the supports alone -- whose faces lie on the rod faces
    # that cut has just made -- was refused on slab after slab, and on one slab
    # the combined cut was refused where that route removed exactly the
    # expected 1,234.62 mm^3 (a support beside junctions it touches). So the
    # combined cut is tried without its slow tool-at-a-time fallback, and the
    # second route, with only the supports' few tools, carries the fallback.
    try:
        freed = _free(region, list(lattice) + list(tools), lattice_in + tools_in,
                      "the lattice and the supports", one_by_one=False)
    except CutRefused:
        freed = _free(before, tools, tools_in, "the supports")
    after, home = [], []
    homes = _Inside(before)
    for p in freed:
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

    (body_path, mesh_path, lp_args, kept, lost_path, tools_path, lo, hi, env_lo, env_hi,
     out) = job
    body = read_brep(body_path)
    mesh = load_mesh(mesh_path)
    probe = OutsideProbe(mesh, 2.0 * mesh.deviation)
    lost = compound_children(read_brep(lost_path)) if lost_path else ()
    local = _LocalLattice(lattice_params(*lp_args), kept, lost)
    lattice = list(local.near((lo, hi))) if len(kept) else []
    tools = [s for s in compound_children(read_brep(tools_path)) if _touches(s, lo, hi)]
    try:
        got = _slab_free(body, probe, lo, hi, env_lo, env_hi, lattice, tools)
    except CutRefused as exc:
        return ("refused", str(exc), peak_rss_bytes())
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


def _inner_points(faces: list) -> list:
    """A point inside each of ``faces``: the centroid of its largest triangle."""
    from .build import _mesh_of

    out = []
    for f in faces:
        for _f, pts, _u, tris in _mesh_of(f, 0.05):
            a, b, c = pts[tris[:, 0]], pts[tris[:, 1]], pts[tris[:, 2]]
            k = int(np.argmax(np.linalg.norm(np.cross(b - a, c - a), axis=1)))
            out.append((a[k] + b[k] + c[k]) / 3.0)
    return out


def _in_any(point, faces: list) -> bool:
    from OCP.BRepClass import BRepClass_FaceClassifier
    from OCP.TopAbs import TopAbs_State

    for f in faces:
        cls = BRepClass_FaceClassifier(f, gp_Pnt(*map(float, point)), PLANE_TOL)
        if cls.State() == TopAbs_State.TopAbs_IN:
            return True
    return False


def _joined_across(left: list, right: list, ax: int, c: float) -> list:
    """``(i, j)`` pairs of pieces either side of the plane ``x[ax] = c`` whose
    faces on it overlap: one piece of free space cut in two by the slabs.

    Both sides' faces on the plane are the same cross-section -- the body less
    the same tools, cut from either side -- so a point inside a face of one
    piece lies inside a face of the other exactly where the two are one piece.
    Asked both ways, by the face classifier: a boolean between two such faces,
    each the whole part's section less hundreds of strut holes, costs seconds,
    and there are two per plane."""
    lf = [_faces_on(p, ax, c) for p in left]
    rf = [_faces_on(p, ax, c) for p in right]
    lp = [_inner_points(a) for a in lf]
    rp = [_inner_points(b) for b in rf]
    pairs = []
    for i, a in enumerate(lf):
        for j, b in enumerate(rf):
            if not a or not b:
                continue
            if any(_in_any(q, b) for q in lp[i]) or any(_in_any(q, a) for q in rp[j]):
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
                 b: np.ndarray, *, pool=None, staged=None, report=None,
                 reach: float = 0.0) -> Pockets:
    """S7 for a group of support ``cells``, which overlap nothing: the fill of
    every pocket.

    ``lattice_near(box)`` returns the lattice junctions near a box (an empty
    list with no lattice). ``body_mesh`` is the classification mesh, whose
    measured deviation sizes the witness tolerance.

    **The envelope is cut into slabs** across its longest extent, one per
    :data:`SLAB_JUNCTIONS` of lattice in it and none thinner than
    :data:`SLAB_REACH_MULTIPLE` junction ``reach``-es, and every slab's free
    space is found on its own -- across the pool when there is one, in this
    process otherwise, by the same code. Pieces of
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
    n = _slab_count(len(lattice_near((lo, hi))), float(hi[ax] - lo[ax]), reach)
    edges = lo[ax] + (hi[ax] - lo[ax]) * np.arange(n + 1) / n
    edges[1:-1] += SLAB_NUDGE
    probe = OutsideProbe(body_mesh, 2.0 * body_mesh.deviation)

    def tick(done: int, total: int) -> None:
        if report is not None:
            report("searching the support envelope for sealed pockets, slab", done, total)

    def box(span):
        slo, shi = lo.copy(), hi.copy()
        slo[ax], shi[ax] = span
        return slo, shi

    def search(spans: list, tools: list) -> list:
        """Per span, its :class:`_Slab`, or the reason the kernel refused it."""
        if len(spans) == 1 or not parallel:
            out = []
            for s, span in enumerate(spans):
                tick(s, len(spans))
                slo, shi = box(span)
                try:
                    out.append(_slab_free(body, probe, slo, shi, lo, hi,
                                          list(lattice_near((slo, shi))),
                                          [t for t in tools if _touches(t, slo, shi)]))
                except CutRefused as exc:
                    out.append(str(exc))
            tick(len(spans), len(spans))
            return out
        import os
        import uuid

        tag = os.path.join(staged["tmpdir"], f"pockets-{uuid.uuid4().hex[:8]}")
        write_brep(_compound(tools), tag + "-tools.brep")
        lp_args = staged["lp_args"]
        jobs = [(staged["body_path"], staged["mesh_path"], lp_args, staged["kept_nodes"],
                 staged.get("lost_path"), tag + "-tools.brep", *box(span), lo, hi, f"{tag}-{s}")
                for s, span in enumerate(spans)]
        tick(0, len(spans))
        results, _rss = pool.run(_slab_worker, jobs, on_result=tick)
        out = []
        for r in results:
            if r[0] == "refused":
                out.append(r[1])
                continue
            nb, na, home, benv, aenv, path, _r = r
            before = compound_children(read_brep(path + "-before.brep")) if nb else []
            after = compound_children(read_brep(path + "-after.brep")) if na else []
            out.append(_Slab(before, after, list(home), list(benv), list(aenv)))
        return out

    def slabs(tools: list):
        """The spans along ``ax`` and their slabs. A refused slab is split in
        two and searched again, :data:`SLAB_RETRY_SPLITS` times; one still
        refused stays a reason string."""
        spans = [(float(edges[s]), float(edges[s + 1])) for s in range(n)]
        got = search(spans, tools)
        for _ in range(SLAB_RETRY_SPLITS):
            if not any(isinstance(g, str) for g in got):
                break
            nspans, keep = [], []
            for span, g in zip(spans, got):
                if isinstance(g, str):
                    mid = 0.5 * (span[0] + span[1]) + SLAB_NUDGE
                    nspans += [(span[0], mid), (mid, span[1])]
                    keep += [None, None]
                else:
                    nspans.append(span)
                    keep.append(g)
            todo = [k for k, g in enumerate(keep) if g is None]
            for k, g in zip(todo, search([nspans[k] for k in todo], tools)):
                keep[k] = g
            spans, got = nspans, keep
        return spans, got

    unsearched, refusals, searched = [], [], {}
    body_lo, body_hi = body_mesh.verts.min(axis=0), body_mesh.verts.max(axis=0)

    def pockets(tools: list) -> list:
        spans, got = slabs(tools)
        # A slab still refused is left open: its free space is not known, so
        # every piece beside it counts as reaching the envelope and nothing
        # it could belong to is judged sealed.
        refused = [k for k, g in enumerate(got) if isinstance(g, str)]
        run = []
        for k in refused + [None]:
            # Neighbours left refused by a split read as one region.
            if run and (k is None or k != run[-1] + 1):
                rlo, rhi = box((spans[run[0]][0], spans[run[-1]][1]))
                rlo, rhi = np.maximum(rlo, body_lo), np.minimum(rhi, body_hi)
                unsearched.append(box_text(rlo, rhi))
                refusals.append(f"{unsearched[-1]}: {got[run[0]]}")
                run = []
            if k is not None:
                run.append(k)
        for k in refused:
            got[k] = _Slab([], [], [], [], [])
        for k in refused:
            for s, c in ((k - 1, spans[k][0]), (k + 1, spans[k][1])):
                if 0 <= s < len(got) and s not in refused:
                    g = got[s]
                    g.before_env = [e or bool(_faces_on(q, ax, c))
                                    for q, e in zip(g.before, g.before_env)]
                    g.after_env = [e or bool(_faces_on(q, ax, c))
                                   for q, e in zip(g.after, g.after_env)]
        n = searched["slabs"] = len(got)
        # Global indices of every slab's pieces, then the pieces joined across
        # the planes between slabs.
        b0 = np.cumsum([0] + [len(s.before) for s in got])
        a0 = np.cumsum([0] + [len(s.after) for s in got])
        bpairs, apairs = [], []
        for s in range(n - 1):
            c = spans[s][1]
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
    # G5 -- with the fills in place nothing is cut off -- holds by construction:
    # every fill is a whole cut-off piece of the free space just found, and what
    # remains is exactly the pieces that were not cut off. Re-finding the free
    # space with the fills as tools cut through nothing but coincident faces,
    # and was refused on a slab of the rehearsal part after an hour of cutting.
    return Pockets(fills, {"support_pockets": len(fills), "support_pocket_slabs": searched.get("slabs", n),
                           "support_pocket_slabs_unsearched": len(unsearched),
                           "support_pocket_volume_mm3":
                               round(sum(occ.volume(f) for f in fills), 4)},
                   unsearched, refusals)


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
        raise CutRefused(
            "Intersecting the support envelope failed (docs/algorithm.md §14.3 S7)."
        )
    return got
