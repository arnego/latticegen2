"""Pipeline steps that must degrade rather than abort.

Same-domain unification (docs/algorithm.md §9) makes the output smaller, not more
correct, so a kernel that refuses to perform it must not end the run — the
failure mode of issue #6's second half.
"""

import numpy as np
import pytest
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.Standard import Standard_Failure

from latticegen2 import occ
from latticegen2.errors import ProcessingError
from latticegen2.parallel import WorkerPool
from latticegen2.pipeline import (
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
    merged, ran = _unify_one(solid)
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
    split, ran = _unify_one(solid)

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
    merged, ran = _unify_one(solid)
    faces, _ = occ.count_subshapes(merged)
    assert seen == [(False, True), (True, False)]
    assert ran  # the kernel did run, just not the edge pass
    assert faces == 6


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


def test_validate_parallel_matches_serial(tmp_path):
    solids = two_solids_that_can_unify()
    serial_invalid, serial_volume, serial_rss = _validate(solids)
    with WorkerPool(2) as pool:
        parallel_invalid, parallel_volume, parallel_rss = _validate(
            solids, pool=pool, tmpdir=str(tmp_path)
        )

    assert parallel_invalid == serial_invalid == []
    assert parallel_volume == pytest.approx(serial_volume)
    assert serial_rss == 0
    assert parallel_rss > 0


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


def test_validate_catches_an_invalid_solid_through_the_pool(tmp_path):
    """The gate must still fire when it runs in a worker, not just on the master."""
    bad_solid = _open_shell_solid()
    assert not occ.is_valid(bad_solid)  # sanity: this really is invalid

    with WorkerPool(2) as pool:
        invalid, _volume, _rss = _validate(
            [two_boxes_sharing_a_face(), bad_solid], pool=pool, tmpdir=str(tmp_path)
        )
    assert invalid == [1]
