"""Putting the supports beside the lattice (docs/algorithm.md §14.2).

The generator returns support material only, and the stage does **not** merge
it with the lattice: the output carries the lattice solids exactly as
``simplify`` left them and, after them, the support bodies, built to share faces
with the lattice and the input's surfaces so the merge in the downstream CAD
tool succeeds (specification.md §4.6).

What a separate body can be held to without a second opinion on its geometry is
that it **overlaps nothing**: not the lattice, and not another support. Volumes
are then additive by construction, which is what the run summary reports.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.TopAbs import TopAbs_ShapeEnum, TopAbs_State
from OCP.TopTools import TopTools_ListOfShape
from OCP.gp import gp_Pnt

from .. import lattice as _lattice
from .. import occ
from ..errors import ProcessingError

SUPPORT_OVERLAP_TOL = 1e-6
"""mm³. The most volume a support may share with the lattice or with another
support. Faces may coincide; volume may not."""

KEY_OFFSET = 1 << 20
KEY_SPAN = 1 << 21
"""Packing of a node index triple into one integer: each index shifted by
:data:`KEY_OFFSET` into ``[0, KEY_SPAN)``, a million cells either way of the
origin -- 0.4 mm cells reach 400 m."""

CONTACT_GAP = 1e-3
"""mm. Bounding boxes closer than this are treated as touching when choosing
which pairs to intersect. Over-selecting costs a boolean, never a wrong answer."""

OVERLAP_GROUP_JUNCTIONS = 32
"""Junction instances per intersection in the overlap check. The junctions do
not overlap each other (docs/algorithm.md §3.2), so a support's overlap with
all of them is the sum of its overlaps with any partition of them, and each
group is one job for the worker pool. Measured on the ball's and the
cylinder's supports: an intersection costs ~45 ms per junction in it, and the
same in groups of 8 to 160 to within 10 %, so small groups cost nothing and
let one large support spread over every worker."""

JUNCTION_FACES = 30
"""Faces of the junction template (docs/algorithm.md §9), used only to put a
group of junctions and a support on one scale when ordering the jobs."""

OVERLAP_GROUP_MIN = 8
"""The smallest group used to give every worker several jobs on a small part:
the bottom of the range measured flat."""


@dataclass
class AssembleResult:
    solids: list
    """The final product: the lattice solids, then the support bodies."""
    lattice_volume: "float | None"
    """``None`` when there were no supports: the volume is then the output's,
    which ``validate`` measures."""
    support_volume: float
    stats: dict = field(default_factory=dict)


def check_result(result, lattice: list) -> None:
    """Fail (exit 4) unless ``result.supports`` are solids and none is the lattice."""
    check_supports(list(result.supports), lattice)


def check_supports(supports: list, lattice: list) -> None:
    """Fail (exit 4) unless every support is a solid and none is the lattice."""
    for i, s in enumerate(supports):
        if s is None or s.IsNull() or s.ShapeType() != TopAbs_ShapeEnum.TopAbs_SOLID:
            raise ProcessingError(
                f"The support generator returned something other than a solid "
                f"at position {i} of {len(supports)} (docs/algorithm.md §14.2)."
            )
        if any(s.IsSame(l) for l in lattice):
            raise ProcessingError(
                f"The support generator returned a lattice solid as support "
                f"at position {i} of {len(supports)}; it must return support "
                f"material only (docs/algorithm.md §14.2)."
            )


def _common_volume(a, tools: list) -> float:
    return float(sum(occ.volume(s) for s in _common_solids(a, tools)))


def _overlap_at(solids: list) -> list:
    """Where the largest of ``solids`` is, for the message naming an overlap."""
    return np.round(occ.centroid(max(solids, key=occ.volume)), 3).tolist()


def _common_solids(a, tools: list) -> list:
    if not tools:
        return []
    algo = BRepAlgoAPI_Common()
    args = TopTools_ListOfShape()
    args.Append(a)
    tl = TopTools_ListOfShape()
    for s in tools:
        tl.Append(s)
    algo.SetArguments(args)
    algo.SetTools(tl)
    algo.Build()
    if not algo.IsDone():
        raise ProcessingError(
            "Checking a support for overlap did not complete: the intersection "
            "failed (docs/algorithm.md §14.2)."
        )
    return occ.solids(algo.Shape())


def _boxes(solids):
    out = []
    for s in solids:
        lo, hi = occ.bounding_box(s)
        out.append((np.asarray(lo) - CONTACT_GAP, np.asarray(hi) + CONTACT_GAP))
    return out


def _touch(box_a, box_b) -> bool:
    return bool(np.all(box_a[0] <= box_b[1]) and np.all(box_b[0] <= box_a[1]))


def _where(box) -> list:
    return np.round((box[0] + box[1]) / 2.0, 3).tolist()


class _LocalLattice:
    """Junction instances near a box: the lattice's material, by construction.

    Built from the junction template at the kept nodes rather than read off the
    lattice solids, whose dominant body has faces in the hundreds of thousands:
    an intersection against it per support would put that body on the support
    path once per cell. The template's faces are the lattice's own faces
    (docs/algorithm.md §3.2), so the two answer the overlap question alike
    inside the body, which is where supports are.

    ``lost`` is the exception: pieces of a kept node's junction that are **not**
    in the lattice, because the floating-body rule dropped them (specification.md
    §5). A trim can cut a strut's tip off from its node; the node stays kept
    through its other arms, and the tip is gone. Measured on
    `TD_HX_rehearsal_test`: 146 ceiling samples under such tips, read as backed
    by a strut that is not there. Neither :meth:`contains` nor :meth:`near`
    counts that material.
    """

    def __init__(self, lp, kept_nodes: np.ndarray, lost=()):
        self.lp = lp
        self._nodes = np.asarray(kept_nodes, dtype=np.int64).reshape(-1, 3)
        self.positions = (_lattice.nodes(lp, kept_nodes) if len(kept_nodes)
                          else np.empty((0, 3)))
        self._template = None
        self._keys = None
        # A junction reaches a/2 along a strut, plus the profile's circumradius.
        self.reach = lp.a / 2.0 + lp.r
        self.lost = list(lost)
        self.lost_boxes = [tuple(np.asarray(x, dtype=float) for x in occ.bounding_box(s))
                           for s in self.lost]
        self._lost_lo = self._lost_hi = None
        self._lost_cls = {}
        self._instances = {}

    def near(self, box) -> list:
        if not len(self.positions):
            return []
        lo, hi = box
        hit = np.all((self.positions >= lo - self.reach)
                     & (self.positions <= hi + self.reach), axis=1)
        return self.instances(np.nonzero(hit)[0])

    def near_solid(self, solid) -> np.ndarray:
        """Indices of the junctions that can share volume with ``solid``, sorted.

        Tighter than :meth:`near` on its box, and exact: a junction lies within
        ``reach`` of its node along every axis, so one whose box meets the box
        of no face of ``solid`` has none of ``solid``'s boundary in it, and is
        wholly inside or wholly outside; its node says which. Inside it is
        kept, and reports its whole volume as overlap."""
        if not len(self.positions):
            return np.zeros(0, dtype=np.int64)
        lo, hi = _boxes([solid])[0]
        R = self.reach
        idx = np.nonzero(np.all((self.positions >= lo - R)
                                & (self.positions <= hi + R), axis=1))[0]
        if not len(idx):
            return idx
        fb = np.array([np.concatenate(occ.bounding_box(f)) for f in occ.faces(solid)])
        flo, fhi = fb[:, :3] - CONTACT_GAP, fb[:, 3:] + CONTACT_GAP
        P = self.positions[idx]
        meets = np.zeros(len(P), dtype=bool)
        for a in range(0, len(P), 512):
            p = P[a:a + 512, None, :]
            meets[a:a + 512] = np.all((p - R <= fhi[None]) & (p + R >= flo[None]),
                                      axis=2).any(axis=1)
        inside = []
        rest = idx[~meets]
        if len(rest):
            cls = BRepClass3d_SolidClassifier(solid)
            for j in rest:
                cls.Perform(gp_Pnt(*map(float, self.positions[j])), 1e-7)
                if cls.State() != TopAbs_State.TopAbs_OUT:
                    inside.append(j)
        return np.sort(np.concatenate([idx[meets], np.asarray(inside, dtype=np.int64)]))

    def groups(self, idx: np.ndarray, size: int = OVERLAP_GROUP_JUNCTIONS) -> list:
        """``idx`` split into groups of about ``size`` junctions near each
        other: sorted by 4-cell block of the lattice index, then cut in runs."""
        if not len(idx):
            return []
        b = self._nodes[idx] // 4
        order = np.lexsort((b[:, 2], b[:, 1], b[:, 0]))
        k = max(1, int(np.ceil(len(idx) / size)))
        return [g for g in np.array_split(idx[order], k) if len(g)]

    def instances(self, idx) -> list:
        """Junction instances ``idx``, less their lost pieces."""
        if not len(idx):
            return []
        if self._template is None:
            from ..junction import build_template

            self._template = build_template(self.lp).solid
        return [self._instance(int(i)) for i in idx]

    def _instance(self, i: int):
        """Junction ``i``'s template, less any lost piece of it.

        A cut that does not take out what it should leaves the template: the
        overlap check then falls back on the lattice solids, and the pocket
        search on its refusal path, so neither is misled."""
        got = self._instances.get(i)
        if got is not None:
            return got
        p = self.positions[i]
        got = self._template.Moved(occ.translation(p))
        for s, (lo, hi) in zip(self.lost, self.lost_boxes):
            if np.any(lo > p + self.reach) or np.any(hi < p - self.reach):
                continue
            v = occ.volume(got)
            cut = BRepAlgoAPI_Cut(got, s)
            if not cut.IsDone():
                continue
            gone = v - occ.volume(cut.Shape())
            if -1e-9 <= gone <= occ.volume(s) * (1.0 + 1e-3):
                got = cut.Shape()
        self._instances[i] = got
        return got

    def contains(self, points) -> np.ndarray:
        """Which ``points`` lie in the lattice's material: the junction
        instances', less what was lost. :meth:`in_template` is the first half."""
        P = np.asarray(points, dtype=float).reshape(-1, 3)
        out = self.in_template(P)
        if not self.lost or not out.any():
            return out
        for j in self._lost_near(P[out]):
            lo, hi = self.lost_boxes[j]
            at = np.nonzero(out & np.all((P >= lo - 1e-6) & (P <= hi + 1e-6), axis=1))[0]
            for i in at:
                if self._in_lost(j, P[i]):
                    out[i] = False
        return out

    def any_kept(self, points) -> bool:
        """Whether any of ``points``, all in a junction's template, is not in a
        lost piece: :meth:`contains` for a caller that needs only the one answer,
        so it stops at the first such point."""
        P = np.asarray(points, dtype=float).reshape(-1, 3)
        near = self._lost_near(P)
        for p in P:
            for j in near:
                lo, hi = self.lost_boxes[j]
                if np.all(p >= lo - 1e-6) and np.all(p <= hi + 1e-6) and self._in_lost(j, p):
                    break
            else:
                return True
        return False

    def _lost_near(self, P: np.ndarray) -> np.ndarray:
        """The lost pieces whose boxes meet the box of ``P``."""
        if self._lost_lo is None:
            self._lost_lo = np.array([b[0] for b in self.lost_boxes]).reshape(-1, 3)
            self._lost_hi = np.array([b[1] for b in self.lost_boxes]).reshape(-1, 3)
        if not len(P):
            return np.zeros(0, dtype=np.int64)
        return np.nonzero(np.all((self._lost_lo <= P.max(axis=0) + 1e-6)
                                 & (self._lost_hi >= P.min(axis=0) - 1e-6), axis=1))[0]

    def _in_lost(self, j: int, p: np.ndarray) -> bool:
        cls = self._lost_cls.get(j)
        if cls is None:
            cls = self._lost_cls[j] = BRepClass3d_SolidClassifier(self.lost[j])
        cls.Perform(gp_Pnt(*map(float, p)), 1e-7)
        return cls.State() in (TopAbs_State.TopAbs_IN, TopAbs_State.TopAbs_ON)

    def is_lost(self, solid) -> bool:
        """Whether ``solid`` is lost material, for one connected solid inside a
        junction's trim: it lies in one connected piece of that trim, and a piece
        is lost or kept whole, so one point inside it decides. No such point reads
        as kept, which fails towards refusing rather than shipping an overlap."""
        from .voids import witness

        w = witness(solid)
        if w is None:
            return False
        return any(self._in_lost(int(j), w) for j in self._lost_near(w.reshape(1, 3)))

    def in_template(self, points) -> np.ndarray:
        """Which ``points`` lie in the junction instances' material, exactly.

        ``J`` is the union of six half-strut prisms (docs/algorithm.md §3.2), so
        a point is inside it iff, for some strut family ``k``, it lies within
        ``a/2`` of the node along ``e_k`` and inside the diamond profile across
        it: ``|x·u_k| + |x·v_k| <= r``. A point within a junction is less than
        one cell from its node in every index coordinate, so the nodes to ask
        are the corners of the index cell holding it.
        """
        P = np.asarray(points, dtype=float).reshape(-1, 3)
        out = np.zeros(len(P), dtype=bool)
        if not len(self.positions) or not len(P):
            return out
        if self._keys is None:
            self._keys = np.sort(self._key(self._nodes))
        lp = self.lp
        q = np.linalg.solve(lp.B, P.T).T
        base = np.floor(q).astype(np.int64)
        for corner in np.array(np.meshgrid([0, 1], [0, 1], [0, 1])).T.reshape(-1, 3):
            n = base + corner
            key = self._key(n)
            at = np.searchsorted(self._keys, key)
            kept = (at < len(self._keys)) & (self._keys[np.minimum(at, len(self._keys) - 1)] == key)
            if not kept.any():
                continue
            d = P[kept] - (lp.B @ n[kept].T).T
            hit = np.zeros(len(d), dtype=bool)
            for k in range(3):
                hit |= ((np.abs(d @ lp.e[k]) <= lp.a / 2.0)
                        & (np.abs(d @ lp.u[k]) + np.abs(d @ lp.v[k]) <= lp.r))
            out[np.nonzero(kept)[0][hit]] = True
        return out

    @staticmethod
    def _key(n: np.ndarray) -> np.ndarray:
        n = np.asarray(n, dtype=np.int64) + KEY_OFFSET
        return (n[:, 0] * KEY_SPAN + n[:, 1]) * KEY_SPAN + n[:, 2]


class _Overlaps:
    """The overlap check's jobs, run alike in this process or in a worker.

    A job is ``("lattice", i, junction indices)`` -- support ``i`` against a
    group of junction instances -- or ``("pair", i, j)``, two supports. Each
    answers ``(volume, where)``: ``where`` is the largest piece's centroid for a
    pair over the bar, the one place a message names, and ``None`` otherwise --
    a lattice reading is confirmed against the lattice solids, which name their
    own."""

    def __init__(self, supports: list, local: _LocalLattice):
        self.supports = supports
        self.local = local

    def run(self, job) -> tuple:
        kind, i, other = job
        if kind == "lattice":
            # Cutting a lost piece out of an instance can fail, so what is left
            # of one here is not an overlap: on the rehearsal part, two tips of
            # 0.008 and 0.002 mm³, which would send a 9,000-face support body to
            # the dominant lattice solid for nothing.
            common = [c for c in _common_solids(self.supports[i],
                                                self.local.instances(other))
                      if not self.local.is_lost(c)]
        else:
            common = _common_solids(self.supports[i], [self.supports[other]])
        vols = [occ.volume(c) for c in common]
        total = float(sum(vols))
        if kind == "pair" and total > SUPPORT_OVERLAP_TOL:
            return total, _overlap_at(common)
        return total, None


_WORKER_OVERLAPS: dict = {}
"""Per worker process: the staged supports and lattice, read once per run."""


def _overlap_worker(job):
    from ..lattice import lattice_params
    from ..parallel import compound_children, read_brep
    from ..runlog import peak_rss_bytes

    (supports_path, kept_path, lost_path, lp_args), task, _cost = job
    chk = _WORKER_OVERLAPS.get(supports_path)
    if chk is None:
        _WORKER_OVERLAPS.clear()
        lost = compound_children(read_brep(lost_path)) if lost_path else ()
        local = _LocalLattice(lattice_params(*lp_args), np.load(kept_path), lost)
        chk = _WORKER_OVERLAPS[supports_path] = _Overlaps(
            compound_children(read_brep(supports_path)), local)
    return chk.run(task) + (peak_rss_bytes(),)


def assemble(lattice: list, supports: list, lp, kept_nodes: np.ndarray,
             report=None, lost=(), pool=None, tmpdir: "str | None" = None) -> AssembleResult:
    """Check ``supports`` overlap nothing, then return lattice + supports.

    Every support is intersected with the junction instances that can reach it
    (:meth:`_LocalLattice.near_solid`), in groups, and with every support whose
    box touches its own. The groups and pairs are independent jobs, run across
    ``pool`` when it is active and ``tmpdir`` is given to stage the shapes in,
    and in this process otherwise, by the same code."""
    if not supports:
        # Nothing to add and nothing to measure: the lattice passes through as
        # the same objects, and `validate` measures its volume anyway.
        return AssembleResult(list(lattice), None, 0.0, {})
    check_supports(supports, lattice)
    t0 = time.perf_counter()

    kept_nodes = np.asarray(kept_nodes)
    local = _LocalLattice(lp, kept_nodes, lost)
    boxes = _boxes(supports)
    n = len(supports)
    parallel = pool is not None and pool.active and tmpdir is not None
    near = [local.near_solid(s) for s in supports]
    n_tools = int(sum(len(x) for x in near))
    # Four jobs a worker where the part is small, within the range measured
    # flat (OVERLAP_GROUP_JUNCTIONS), so no worker waits on one long job.
    workers = pool.workers if parallel else 1
    size = int(np.clip(np.ceil(n_tools / (4 * workers)),
                       OVERLAP_GROUP_MIN, OVERLAP_GROUP_JUNCTIONS))
    # A job's cost is estimated as the support's faces times the tools' faces,
    # one unit for groups and pairs alike, so dispatch order compares them.
    nf = [len(occ.faces(s)) for s in supports]
    jobs, cost = [], []
    for i in range(n):
        for g in local.groups(near[i], size):
            jobs.append(("lattice", i, g))
            cost.append(nf[i] * len(g) * JUNCTION_FACES)
    pairs = 0
    for i in range(n):
        for j in range(i + 1, n):
            if _touch(boxes[i], boxes[j]):
                jobs.append(("pair", i, j))
                cost.append(nf[i] * nf[j])
                pairs += 1
    select_s = time.perf_counter() - t0

    label = "checking supports for overlap"
    if report is not None:
        report(label, 0, len(jobs))
    if parallel and len(jobs) > 1:
        import os
        import uuid

        from ..parallel import write_brep

        tag = os.path.join(tmpdir, f"overlap-{uuid.uuid4().hex[:8]}")
        write_brep(occ.compound(list(supports)), tag + "-supports.brep")
        np.save(tag + "-kept.npy", kept_nodes)
        lost_path = None
        if lost:
            lost_path = tag + "-lost.brep"
            write_brep(occ.compound(list(lost)), lost_path)
        stage = (tag + "-supports.brep", tag + "-kept.npy", lost_path,
                 (lp.cc, lp.t, tuple(map(float, lp.orient))))
        # Largest first, so the one large support's groups start at once and
        # the small jobs fill in behind them.
        got, _rss = pool.run(
            _overlap_worker, [(stage, job, c) for job, c in zip(jobs, cost)],
            sort_by=lambda w: w[2],
            on_result=(lambda done, total: report(label, done, total)) if report else None)
        results = [r[:2] for r in got]
    else:
        chk = _Overlaps(supports, local)
        results = []
        for k, job in enumerate(jobs):
            results.append(chk.run(job))
            if report is not None:
                report(label, k + 1, len(jobs))

    # Junction overlaps sum per support over its groups: the junctions share no
    # volume, so the groups' answers add up to the whole.
    against_lattice = np.zeros(n)
    for job, (vol, _at) in zip(jobs, results):
        if job[0] == "lattice":
            against_lattice[job[1]] += vol
    lattice_boxes = None
    for i, s in enumerate(supports):
        if against_lattice[i] > SUPPORT_OVERLAP_TOL:
            # A junction instance is the template, untrimmed, so a reading here
            # is confirmed against the lattice solids themselves before it is
            # believed.
            if lattice_boxes is None:
                lattice_boxes = _boxes(lattice)
            common = _common_solids(
                s, [l for l, lb in zip(lattice, lattice_boxes) if _touch(boxes[i], lb)])
            real = float(sum(occ.volume(c) for c in common))
            if real > SUPPORT_OVERLAP_TOL:
                raise ProcessingError(
                    f"A support body overlaps the lattice by {real:.3e} mm^3 near "
                    f"{_overlap_at(common)}; supports must share faces with the "
                    f"lattice, never volume (docs/algorithm.md §14.2)."
                )
    for job, (vol, at) in zip(jobs, results):
        if job[0] == "pair" and vol > SUPPORT_OVERLAP_TOL:
            raise ProcessingError(
                f"Two support bodies overlap by {vol:.3e} mm^3 near "
                f"{at}; neighbouring supports must share faces, "
                f"never volume (docs/algorithm.md §14.2)."
            )

    lattice_volume = float(sum(occ.volume(s) for s in lattice))
    support_volume = float(sum(occ.volume(s) for s in supports))
    return AssembleResult(list(lattice) + list(supports), lattice_volume,
                          support_volume, {
                              "support_bodies": n,
                              "support_pairs_checked": pairs,
                              "support_overlap_jobs": len(jobs),
                              "support_overlap_junctions": n_tools,
                              "support_overlap_select_s": round(select_s, 2),
                              "support_overlap_s": round(time.perf_counter() - t0, 2),
                          })
