"""The support generator (docs/algorithm.md §14.3).

Ceiling points (:mod:`.regions`), frustum supports on the lattice grid
(:mod:`.frustum`), the pockets they seal (:mod:`.voids`) and the joining of
bodies that share a face (:mod:`.join`), in that order.

:func:`generate` receives the completed lattice and the input body
(:class:`.contract.SupportInput`) and returns the support solids alone
(:class:`.contract.SupportResult`).

What happens to them is not the generator's job: the support stage checks
they overlap neither the lattice nor each other and appends them beside the
lattice as separate bodies (:mod:`.assemble`); then ``validate`` ->
STEP round-trip check -> ``export`` run on the final product. What an implementation
owes that chain:

* **Inside the body.** No material outside the input solid (specification.md
  §1). A boolean against the body can be wrong while reporting success
  (docs/algorithm.md §7.2): every trim is held to a sampled estimate and to
  ``classify.OutsideProbe`` over ``inp.mesh``.
* **Coverage.** Every ceiling point inside support or lattice material --
  ``tools/verify_geometry.overhang_coverage`` asks this of the written file in
  e2e.
* **Self-supporting.** No support surface or edge overhangs more than the rod
  angle :data:`.params.SUPPORT_OVERHANG_DEG` from ``inp.lp.b`` -- by
  construction, every face lying on a plane of the lattice frame.
* **Separate, manifold bodies.** Separate from the lattice, sharing faces with
  it and with the input's surfaces and overlapping nothing; supports sharing a
  face or volume with each other are one body (:mod:`.join`). No edge of a
  body is used by more than two of its faces. Support bodies are exempt from
  the floating-body rule (specification.md §5).
* **Valid, writable solids.** Closed and ``BRepCheck_Analyzer``-valid, and
  surviving a STEP round trip, which ``validate`` checks.
* **Supports only.** Never a copy of the lattice.
"""

from __future__ import annotations

import os
import tempfile
import time

import numpy as np

from ..runlog import WARNING_PREFIX, point_text, points_text, warning_text
from .contract import SupportInput, SupportResult

UNBACKED_WARNING = "ceiling left without support"
"""Title of the warning naming ceiling points nothing in the output backs."""


def _unbacked_warning(bare: np.ndarray, total: int) -> str:
    return warning_text(
        UNBACKED_WARNING,
        f"{len(bare)} of {total} ceiling point(s) are backed by nothing in the output",
        "Those overhangs may print unsupported",
        [point_text(p) for p in bare],
    )


def _touching(a, b, gap: float = 1e-3) -> bool:
    return bool(np.all(a[0] <= b[1] + gap) and np.all(b[0] <= a[1] + gap))


def remove_micro_caps(solids: list) -> tuple[list, list, list]:
    """``(solids, removed, refused)``: every body with its sub-micron caps
    removed (:func:`latticegen2.occ.remove_micro_caps`).

    ``removed`` holds one position per cap taken out, ``refused`` one per body
    whose repair failed its check and was kept as built. The repair is exact by
    construction, and checked rather than trusted: the result must be one
    valid solid with every edge used by exactly two faces, and its boundary may
    not move by more than :data:`latticegen2.occ.MICRO_CAP_TOL` -- read, as §9
    reads unification, as ``|dV| / surface area``, since the volume is
    re-integrated over the host face and moves by quadrature noise (5.7e-6 mm^3
    on the 6,442 mm^3 body this was found on, for a cap of ~1e-12 mm^3). A
    refused repair leaves the body to ``validate``, whose STEP round-trip check
    still judges it.
    """
    from .. import occ
    from . import frustum

    out, removed, refused = [], [], []
    for s in solids:
        fixed, where = occ.remove_micro_caps(s)
        if not where:
            out.append(s)
            continue
        got = occ.solids(fixed)
        ok = len(got) == 1 and occ.is_valid(got[0], parallel=False)
        if ok:
            ok = frustum.edges_not_used_twice(got[0]) == 0
        if ok:
            moved = abs(occ.volume(got[0]) - occ.volume(s)) / max(occ.area(s), 1e-30)
            ok = moved <= occ.MICRO_CAP_TOL
        if ok:
            out.append(got[0])
            removed += where
        else:
            out.append(s)
            refused += where
    return out, removed, refused


def absorb_lost(solids: list, lost, bare: np.ndarray, added: list | None = None
                ) -> tuple[list, int]:
    """``(solids, n)``: ``solids`` with the ``n`` lost junction pieces that hold a
    ``bare`` point made part of them.

    A fine support node a wall crosses, with a lost strut tip on one side and
    kept lattice on the other, stays lattice, so the tip's share of it stays
    empty, and the ceiling above it is backed only if the supports around happen
    to seal it as a pocket. The tip itself is the exact solid that is missing: it
    lies in the body, leaves the kept lattice along faces, and its faces lie on
    rod planes and the body's surface. What the supports do not already hold of
    it joins them. A piece that would make a body non-manifold is left out, for
    the count at the end to report. What joins is appended to ``added``.
    """
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut

    from .. import occ
    from . import build, frustum, join

    n = 0
    for piece in lost:
        if not frustum.exact_inside(piece, bare).any():
            continue
        box = occ.bounding_box(piece)
        near = [s for s in solids if _touching(occ.bounding_box(s), box)]
        got = build._boolean(BRepAlgoAPI_Cut, piece, near) if near else piece
        rest = [] if got is None else occ.solids(got)
        rest = [r for r in rest if occ.volume(r) > 1e-3 and occ.is_valid(r, parallel=False)]
        if not rest or sum(occ.volume(r) for r in rest) > occ.volume(piece) * (1.0 + 1e-6):
            continue
        joined, _js = join.join(near + rest, overlaps=False)
        if any(frustum.non_manifold_edges(s) for s in joined):
            continue
        keep = [s for s in solids if not any(s is o for o in near)]
        solids = keep + joined
        if added is not None:
            added += rest
        n += 1
    return solids, n


def supports_for(body, mesh, lattice: list, kept_nodes, lp, thickness: float,
                 report=None, tmpdir: str | None = None, pool=None, staged=None,
                 pockets_begin=None, samples=None, lost=()):
    """``(solids, stats, lines)``: every support body for ``body``.

    The generator proper, shared by :func:`generate` and by
    ``tools/support_dev.py``, which calls it with no lattice: every anchor is
    then on a wall. ``pool`` and ``staged`` -- the body's ``.brep`` and the
    mesh's ``.npz`` -- let the pocket search run across the run's workers;
    without them it runs in this process. ``samples`` is what the detector
    found, detected here if not given; ``lost`` the junction pieces the
    floating-body rule took from kept nodes.
    """
    from .. import occ
    from . import build, frustum, join, overhang, regions, voids
    from .assemble import _LocalLattice

    t0 = time.perf_counter()
    local = _LocalLattice(lp, np.asarray(kept_nodes), lost)
    if samples is None:
        samples = overhang.detect(body, lp.b)
    # Ceiling a strut ends on is backed by the lattice already.
    rmesh = regions.build(body, lp.b, report=report,
                          covered=local.contains if len(kept_nodes) else None)
    gates = build.Gates(mesh, lp.b)
    ceiling_s = time.perf_counter() - t0
    if tmpdir is None:
        tmpdir = tempfile.mkdtemp(prefix="latticegen2-support-")
    # Keels and tips have no face to refine: their samples are ceiling points
    # as they stand, asked what backs them straight below (§14.1).
    other = samples.kind != overhang.FACE
    got = frustum.build(body, mesh, rmesh, local, lp, thickness, gates, tmpdir,
                        report=report, extra=(samples.points[other], -samples.inward[other]))
    stats = dict(got.stats)
    stats["support_ceiling_s"] = round(ceiling_s, 2)
    lines = list(got.lines) + list(got.notes)
    covered = got.covered.copy()
    shallow = got.points - frustum.SHALLOW_DEPTH * got.normals
    bodies = list(got.solids)
    # Inside the body these are the bodies, with no face on its surface: what
    # the pocket search cuts with where the bodies themselves are refused.
    cores = list(got.cores)
    if lost and not covered.all():
        bodies, stats["support_lost_pieces_absorbed"] = absorb_lost(bodies, lost,
                                                                    shallow[~covered], cores)
    if not bodies:
        # Nothing to search pockets in or to ask again, but what the grid left
        # bare -- every piece dropped, say -- is still reported.
        stats["support_bodies"] = 0
        if not covered.all():
            lines.append(_unbacked_warning(got.points[~covered], len(covered)))
        return [], stats, lines

    # S7 per group of bodies whose envelopes touch: free space the supports seal
    # off is filled, and a fill joins the bodies whose faces it lies on.
    if pockets_begin is not None:
        pockets_begin()
    t3 = time.perf_counter()
    fills, unsearched, refusals, widenings = [], [], [], 0
    gs = join.groups(bodies, grow=thickness)
    for k, g in enumerate(gs):
        if report is not None:
            report("looking for sealed pockets", k, len(gs))
        mine = set(g)
        pockets = voids.find_pockets(body, mesh, [bodies[i] for i in g], local.near,
                                     thickness,
                                     others=[s for i, s in enumerate(bodies) if i not in mine],
                                     pool=pool, staged=staged, report=report,
                                     reach=local.reach, cores=cores)
        fills += [f for f in pockets.fills if not voids.found_before(f, fills)]
        unsearched += pockets.unsearched
        refusals += pockets.refusals
        widenings += pockets.stats.get("support_pocket_widenings", 0)
    out = list(bodies)
    if fills:
        every = bodies + fills
        out = []
        for g in join.groups(every):
            members = [every[i] for i in g]
            if any(i >= len(bodies) for i in g):
                members, _js = join.join(members, overlaps=False)
            out += members
    # A support node's corner within a micron of the input's wall leaves a cap
    # no mesher can close, which the STEP round-trip check then refuses (§14.3 S6).
    out, caps, caps_refused = remove_micro_caps(out)
    stats["support_micro_caps_removed"] = len(caps)
    if caps:
        lines.append(f"note: {len(caps)} sub-micron cap(s) removed from support bodies "
                     f"where a support node's corner grazed the input's surface, at "
                     f"{points_text(caps)}")
    if caps_refused:
        lines.append(f"note: sub-micron cap repair refused by its check, bodies kept as built "
                     f"for the STEP round-trip check to judge, "
                     f"at {points_text(caps_refused)}")
    for s in out:
        nm = frustum.non_manifold_edges(s)
        if nm:
            from ..errors import ProcessingError

            where = np.round(np.asarray(occ.centroid(s), dtype=float), 3).tolist()
            raise ProcessingError(
                f"A support body near {where} has {nm} edge(s) used by more than two "
                f"faces: it is not a manifold solid (docs/algorithm.md §14.3)."
            )
    stats["support_pockets_s"] = round(time.perf_counter() - t3, 2)
    stats["support_pockets"] = len(fills)
    stats["support_pocket_volume_mm3"] = round(sum(occ.volume(f) for f in fills), 4)
    stats["support_bodies"] = len(out)
    stats["support_pocket_slabs_unsearched"] = len(unsearched)
    stats["support_pocket_widenings"] = widenings
    if unsearched:
        lines.append(warning_text(
            voids.UNSEARCHED_WARNING,
            f"the geometry kernel could not cut {len(unsearched)} region(s) of the "
            f"support envelope, even split in halves, so they were not searched",
            "Powder may be trapped in enclosed pockets there",
            unsearched,
        ))
        lines += [f"pocket search refused: {r}" for r in refusals]
    # The grid's answer was a first one: a pocket fill or an absorbed piece can
    # back what it could not. The final bodies decide, and are few to ask about.
    for s in out:
        if covered.all():
            break
        todo = np.nonzero(~covered)[0]
        covered[todo] = (frustum.exact_inside(s, shallow[todo])
                         | frustum.exact_inside(s, got.points[todo] - frustum.PROBE_DEPTH
                                                * got.normals[todo]))
    stats["support_ceiling_points_uncovered"] = int((~covered).sum())
    if not covered.all():
        lines.append(_unbacked_warning(got.points[~covered], len(covered)))
    lines.append(
        f"support: {len(out)} bod(y/ies), {sum(occ.volume(s) for s in out):.3f} mm^3, "
        f"{len(fills)} sealed pocket(s) filled; "
        f"{int(covered.sum())} of {len(covered)} ceiling point(s) backed"
    )
    lines.append(
        "support timing: " + ", ".join(
            f"{k[len('support_'):-2]} {stats[k]:.1f} s" for k in
            ("support_ceiling_s", "support_grid_s", "support_staircase_s",
             "support_trim_s", "support_backed_s", "support_pockets_s"))
    )
    return out, stats, lines


def generate(inp: SupportInput) -> SupportResult:
    """Support solids backing the ceiling of ``inp.body``, anchored on ``inp.lattice``."""
    staged = {"body_path": inp.body_path, "mesh_path": inp.mesh_path,
              "lp_args": (inp.lp.cc, inp.lp.t, tuple(inp.lp.orient)),
              "kept_nodes": np.asarray(inp.kept_nodes), "tmpdir": inp.tmpdir}
    if inp.lost_lattice:
        from .. import occ
        from ..parallel import write_brep

        staged["lost_path"] = os.path.join(inp.tmpdir, "lost-lattice.brep")
        write_brep(occ.compound(list(inp.lost_lattice)), staged["lost_path"])
    solids, stats, lines = supports_for(inp.body, inp.mesh, list(inp.lattice),
                                        inp.kept_nodes, inp.lp, inp.thickness_mm,
                                        report=inp.report, tmpdir=inp.tmpdir,
                                        pool=inp.pool, staged=staged,
                                        pockets_begin=inp.pockets_begin,
                                        samples=inp.samples, lost=inp.lost_lattice)
    # What a run did not do goes to the console whatever -v says.
    notes = [line for line in lines if line.startswith(WARNING_PREFIX)]
    lines = [line for line in lines if not line.startswith(WARNING_PREFIX)]
    return SupportResult(supports=solids, stats=stats, lines=lines, notes=notes)
