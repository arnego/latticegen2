"""Pipeline steps that must degrade rather than abort.

Same-domain unification (docs/algorithm.md §9) makes the output smaller, not more
correct, so a kernel that refuses to perform it must not end the run — the
failure mode of issue #6's second half.
"""

import os

import numpy as np
import pytest
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.Standard import Standard_Failure

from latticegen2 import occ
from latticegen2.errors import ProcessingError
from latticegen2.parallel import WorkerPool
from latticegen2.pipeline import (
    UNIFY_MAX_DISPLACEMENT,
    UNIFY_VOLUME_TOL,
    _check_unify_result,
    _unify,
    _unify_one,
    _validate,
)


def two_boxes_sharing_a_face():
    """One solid whose boundary carries a redundant, mergeable partition."""
    a = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
    b = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
    from OCP.gp import gp_Trsf, gp_Vec

    move = gp_Trsf()
    move.SetTranslation(gp_Vec(10.0, 0.0, 0.0))
    b = BRepBuilderAPI_Transform(b, move, True).Shape()
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse

    fused = BRepAlgoAPI_Fuse(a, b).Shape()
    return occ.solids(fused)[0]


def test_unification_merges_the_redundant_partition():
    solid = two_boxes_sharing_a_face()
    before, _ = occ.count_subshapes(solid)
    merged, ran, _degraded = _unify_one(solid)
    after, _ = occ.count_subshapes(merged)
    assert ran
    assert after < before
    assert after == 6  # a fused pair of boxes is one 20x10x10 box


def test_split_unification_matches_a_single_combined_call():
    """Running the two passes separately must describe the solid identically.

    `_unify_one` merges faces first and concatenates edges in a second call,
    rather than asking for both at once. The sub-body tiling that originally
    required the split is disproved (docs/specification.md §11, G15), but the
    split is kept for its degradation behaviour — a throwing edge pass no longer
    discards a completed face merge. Either way the reordering is only
    admissible while it produces the same B-rep, which is what this pins: same
    faces *and* same edges as the combined call.
    """
    solid = two_boxes_sharing_a_face()
    combined = occ.unify_same_domain(solid)
    split, ran, _degraded = _unify_one(solid)

    assert ran
    assert occ.count_subshapes(split) == occ.count_subshapes(combined)
    assert occ.volume(split) == pytest.approx(occ.volume(combined))


def test_a_failing_edge_pass_keeps_the_face_merge(monkeypatch):
    """The edge pass is what throws; losing it must not lose the face merge.

    It now runs last and alone, so a refusal costs only the edge concatenation —
    the face merge that carries the value has already happened and is kept.
    """
    solid = two_boxes_sharing_a_face()
    real = occ.unify_same_domain
    seen = []

    def refuse_edge_merging(shape, unify_edges=True, unify_faces=True):
        seen.append((unify_edges, unify_faces))
        if unify_edges:
            raise Standard_Failure("Courbes non jointives")
        return real(shape, unify_edges=False)

    monkeypatch.setattr(occ, "unify_same_domain", refuse_edge_merging)
    merged, ran, _degraded = _unify_one(solid)
    faces, _ = occ.count_subshapes(merged)
    assert seen == [(False, True), (True, False)]
    assert ran  # the kernel did run, just not the edge pass
    assert faces == 6


def test_an_invalid_merge_keeps_the_un_unified_solid(monkeypatch):
    """A kernel that does not throw can still hand back an invalid solid.

    docs/algorithm.md §9 names the validity gate as the *stronger* guard on this
    step -- the volume bars cannot see a bad merge, and on
    `TD_HX_rehearsal_test` at ``cc=7, t=1.4`` they did not: a sound 279,358-face
    body unified to 193,721 faces with one invalid face of 3.426439 mm², at a
    volume drift of 6.20e-06, comfortably inside the 1e-4 pre-filter. The run
    then *failed* at `validate`, which is the one response §11 forbids for a
    step whose only job is to make the output smaller.

    So the result is checked and the input kept. What is asserted here is the
    behaviour, not the plumbing: the solid that comes back is the one that went
    in, and the run is told so.
    """
    solid = two_boxes_sharing_a_face()
    faces_in, edges_in = occ.count_subshapes(solid)

    # Signature-faithful: `_unify_one` asks for `parallel=False` because it
    # runs inside a worker, and a stub that did not accept it would pass
    # while the real call raised.
    monkeypatch.setattr(occ, "is_valid", lambda shape, parallel=True: False)
    merged, ran, degraded = _unify_one(solid)

    assert ran, "the kernel did run -- it is the result that did not hold up"
    assert degraded, "and the run has to be able to say so"
    assert occ.count_subshapes(merged) == (faces_in, edges_in), (
        "the un-unified solid is what is kept, so the output is larger and not "
        "different"
    )


def test_the_workers_validity_check_does_not_ask_occt_for_threads(monkeypatch):
    """`--cores` is honoured exactly, and this call is inside a worker.

    docs/algorithm.md S9 keeps the *gate* on the master precisely because W
    worker processes each launching W OCCT threads is W^2 threads on W cores,
    and `parallel.set_thread_budget` is deliberately not called in the worker
    initializer -- so OCCT's default pool there sizes itself to the machine.
    A degrade check that asked for threads would reinstate exactly that.
    """
    seen = []

    def record(shape, parallel=True):
        seen.append(parallel)
        return True

    monkeypatch.setattr(occ, "is_valid", record)
    _unify_one(two_boxes_sharing_a_face())

    assert seen == [False], "the worker must ask for a single-threaded check"


def test_a_valid_merge_is_not_degraded():
    """The control beside it: an ordinary merge must still merge.

    A degrade that fired on sound geometry would quietly double every output
    file, which is the failure this project would notice last.
    """
    solid = two_boxes_sharing_a_face()
    faces_in, _ = occ.count_subshapes(solid)
    merged, ran, degraded = _unify_one(solid)
    faces_out, _ = occ.count_subshapes(merged)

    assert ran and not degraded
    assert faces_out < faces_in, "the redundant partition is still merged away"


def test_an_already_minimal_solid_is_not_reported_as_a_refusal():
    """Unifying to itself is success, not a kernel that declined.

    docs/algorithm.md §9 records the junction template doing exactly this — 30
    faces to 30 — so a stat inferred from "the face count did not change" would
    make the run claim a failure that never happened.
    """
    minimal = BRepPrimAPI_MakeBox(10.0, 10.0, 10.0).Shape()
    before, _ = occ.count_subshapes(minimal)

    _, stats = _unify([occ.solids(minimal)[0]])

    assert stats["output_faces"] == before  # nothing to merge on a bare box
    assert stats["unmerged_solids"] == 0


def test_unification_that_cannot_run_at_all_returns_the_solid_unchanged(monkeypatch):
    """No merge is a bigger file, never a failed run or a different body."""
    solid = two_boxes_sharing_a_face()

    def always_refuse(shape, unify_edges=True, unify_faces=True):
        raise Standard_Failure("Courbes non jointives")

    monkeypatch.setattr(occ, "unify_same_domain", always_refuse)
    out, stats = _unify([solid])

    assert out[0] is solid
    assert stats["unmerged_solids"] == 1
    assert stats["output_faces"] == stats["faces_before"]
    assert stats["volume_drift"] == pytest.approx(0.0, abs=UNIFY_VOLUME_TOL)
    assert occ.volume(out[0]) == pytest.approx(2000.0)


# --- the shared guard, decoupled from where unification ran ------------------
#
# `_check_unify_result` is the one place both the serial and the parallel path
# (docs/specification.md §10 paths 1 and 4) decide pass/fail, so pinning it on
# its own — no geometry, no worker process — is what keeps the two paths from
# silently drifting apart in what they consider a pass.


def test_check_unify_result_accepts_an_exact_merge():
    assert _check_unify_result(100.0, 100.0, 1) == pytest.approx(0.0)


def test_check_unify_result_accepts_drift_within_tolerance():
    post = 100.0 * (1 + 0.5 * UNIFY_VOLUME_TOL)
    drift = _check_unify_result(100.0, post, 1)
    assert 0.0 < drift < UNIFY_VOLUME_TOL


def test_check_unify_result_rejects_drift_beyond_tolerance():
    post = 100.0 * (1 + 2 * UNIFY_VOLUME_TOL)
    with pytest.raises(ProcessingError, match="changed a solid's volume"):
        _check_unify_result(100.0, post, 1)


# --- and the second stage of it ---------------------------------------------
#
# Relative drift is only the pre-filter. What decides is the displacement it
# implies over the solid's own surface, because the same boundary movement
# shows up as a far larger *fraction* of a small solid's volume than of a big
# one's -- so the relative figure is biased by size rather than by the thing
# the guard is looking for. The two solids below are the real pair from one
# `SpiralTest.step` run, unified by the same call in the same stage.


def test_check_unify_result_accepts_a_small_solid_whose_boundary_barely_moved():
    """The regression: `SpiralTest`'s 4.17 mm^3 island, at 1.8e-03 relative.

    Refused outright by the relative bar, and sound -- the exact symmetric
    difference against the un-unified solid is 0 mm^3 both ways, and the merged
    edges carry recorded tolerances of 2.1e-02 mm, 85x the movement implied
    here. Its displacement is 2.494e-04 mm.
    """
    pre, post, area = 4.170537957, 4.162877342, 30.7213
    drift = _check_unify_result(pre, post, 1, area_of=lambda: area)
    assert drift > UNIFY_VOLUME_TOL, "the pre-filter must still trip, or nothing is tested"
    assert abs(pre - post) / area < UNIFY_MAX_DISPLACEMENT


def test_check_unify_result_still_rejects_a_boundary_that_really_moved():
    """Same solid, same area, a drift an order past what displacement allows."""
    pre, area = 4.170537957, 30.7213
    post = pre - 10 * UNIFY_MAX_DISPLACEMENT * area
    with pytest.raises(ProcessingError, match="boundary displacement"):
        _check_unify_result(pre, post, 1, area_of=lambda: area)


def test_check_unify_result_does_not_measure_area_unless_the_pre_filter_trips():
    """The whole reason `area_of` is a callable: 7.7 s on a 64k-face solid.

    Measuring it on every solid would put ~70 s on every rehearsal-scale run
    for a guard that almost never fires.
    """
    calls = []

    def area_of():
        calls.append(1)
        return 30.7213

    _check_unify_result(100.0, 100.0 * (1 + 0.5 * UNIFY_VOLUME_TOL), 1, area_of=area_of)
    assert calls == [], "a solid inside the pre-filter must not pay for its area"


def test_check_unify_result_without_an_area_falls_back_to_the_pre_filter():
    """A caller with no shape in hand keeps the old, stricter behaviour."""
    post = 100.0 * (1 + 2 * UNIFY_VOLUME_TOL)
    with pytest.raises(ProcessingError, match="changed a solid's volume"):
        _check_unify_result(100.0, post, 1)


def test_check_unify_result_rejects_a_solid_count_change():
    with pytest.raises(ProcessingError, match="turned one solid into 2"):
        _check_unify_result(100.0, 100.0, 2)


def test_check_unify_result_rejects_zero_solids_produced():
    with pytest.raises(ProcessingError, match="turned one solid into 0"):
        _check_unify_result(100.0, float("nan"), 0)


# --- parallel dispatch: same result as the master, via a real worker ---------
#
# `_unify_one`'s own ladder (the fallback tests above) never changes between
# the two paths — only *where* it runs does — so what these need to prove is
# dispatch, not the ladder a second time. Note this cannot exercise the
# "kernel refuses inside a worker" branch specifically: the `monkeypatch`
# fallback tests above patch `occ.unify_same_domain` in this process, and a
# `spawn` worker re-imports the module fresh, so that patch does not cross the
# process boundary. `_unify_one` itself is identical code on both paths, so
# the ladder's correctness is not in question — only IPC is, and that is what
# these test.


def two_solids_that_can_unify():
    return [two_boxes_sharing_a_face(), two_boxes_sharing_a_face()]


def test_unify_parallel_matches_serial_on_the_same_solids(tmp_path):
    serial_out, serial_stats = _unify(two_solids_that_can_unify())
    with WorkerPool(2) as pool:
        parallel_out, parallel_stats = _unify(
            two_solids_that_can_unify(), pool=pool, tmpdir=str(tmp_path)
        )

    assert parallel_stats["faces_before"] == serial_stats["faces_before"]
    assert parallel_stats["output_faces"] == serial_stats["output_faces"]
    assert parallel_stats["output_edges"] == serial_stats["output_edges"]
    assert parallel_stats["unmerged_solids"] == serial_stats["unmerged_solids"]
    assert [occ.volume(s) for s in parallel_out] == pytest.approx(
        [occ.volume(s) for s in serial_out]
    )


def test_unify_only_dispatches_to_the_pool_when_eligible(tmp_path):
    solids = two_solids_that_can_unify()

    # No pool given at all: the serial path, same as every existing test above.
    _, no_pool_stats = _unify(solids)
    assert no_pool_stats["max_worker_rss"] == 0

    # A pool is given and there is more than one solid: dispatched for real —
    # `max_worker_rss` is only ever nonzero if a worker process actually ran
    # and reported its own peak (the same proof test_weld.py's tiling test
    # uses for the same reason).
    with WorkerPool(2) as pool:
        _, pool_stats = _unify(solids, pool=pool, tmpdir=str(tmp_path))
    assert pool_stats["max_worker_rss"] > 0

    # A single solid: not worth a process round-trip regardless of the pool.
    with WorkerPool(2) as pool:
        _, single_stats = _unify([two_boxes_sharing_a_face()], pool=pool, tmpdir=str(tmp_path))
    assert single_stats["max_worker_rss"] == 0


def test_validate_passes_sound_solids_and_sums_their_volume():
    solids = two_solids_that_can_unify()
    invalid, volume = _validate(solids)
    assert invalid == []
    assert volume == pytest.approx(sum(occ.volume(s) for s in solids))


def _open_shell_solid():
    """A ``TopoDS_Solid`` wrapping an open shell — one face short of a box.

    Deliberately invalid: ``BRepCheck_Analyzer`` flags a non-closed shell,
    which is what this needs to prove the gate fires from a worker and not
    only on the master.
    """
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Shell

    box = BRepPrimAPI_MakeBox(5.0, 5.0, 5.0).Shape()
    faces = occ.faces(box)
    shell = TopoDS_Shell()
    builder = BRep_Builder()
    builder.MakeShell(shell)
    for f in faces[:-1]:
        builder.Add(shell, f)
    return occ.make_solid(shell)


def test_validate_catches_an_invalid_solid():
    """The gate must fire, and must name *which* solid failed."""
    bad_solid = _open_shell_solid()
    assert not occ.is_valid(bad_solid)  # sanity: this really is invalid

    invalid, _volume = _validate([two_boxes_sharing_a_face(), bad_solid])
    assert invalid == [1]


def test_occt_parallel_flag_does_not_change_a_validity_verdict():
    """Gate G18's part A, pinned as a regression rather than left in a report.

    :func:`latticegen2.occ.is_valid` runs ``BRepCheck_Analyzer`` with OCCT's own
    ``theIsParallel`` on (docs/algorithm.md §9). That is a threading change to
    the pipeline's exit-4 correctness gate, so what has to hold is not that it
    is faster but that it answers **identically** — and identically on faults
    above all. A check that only ever agrees on sound geometry proves nothing:
    that is G10's lesson in this codebase, where a blind scanner reporting "no
    faults" cost real time.

    ``test_weld.py`` covers the same property on the four *real* invalid faces
    committed from the `cc=5, t=1` rehearsal, by way of the ``occ.is_valid``
    calls in its own fixtures; this covers a whole invalid **solid**, which is
    the shape ``validate`` actually receives.
    """
    from OCP.BRepCheck import BRepCheck_Analyzer

    for shape in (two_boxes_sharing_a_face(), _open_shell_solid()):
        serial = BRepCheck_Analyzer(shape, True, False).IsValid()
        assert occ.is_valid(shape) == serial

    assert occ.is_valid(two_boxes_sharing_a_face())
    assert not occ.is_valid(_open_shell_solid())


def test_thread_budget_is_bounded_and_never_raises():
    """``--cores`` must keep meaning what specification.md §3 says it means.

    OCCT's default pool sizes itself to the machine, so without this a
    ``--cores 2`` run could launch a thread per physical core inside
    ``BRepCheck_Analyzer``. The call is also on the run's startup path, so it
    must never throw — a thread-count hint is not a reason to fail a run.
    """
    from OCP.OSD import OSD_ThreadPool

    from latticegen2.parallel import set_thread_budget

    try:
        set_thread_budget(2)
        assert OSD_ThreadPool.DefaultPool_s().NbDefaultThreadsToLaunch() == 2
        set_thread_budget(0)  # nonsense in, still bounded at 1, still no raise
        assert OSD_ThreadPool.DefaultPool_s().NbDefaultThreadsToLaunch() == 1
    finally:
        set_thread_budget(os.cpu_count() or 1)
