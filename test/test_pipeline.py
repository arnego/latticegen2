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
from latticegen2.pipeline import UNIFY_VOLUME_TOL, _unify, _unify_one


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


def test_unification_falls_back_to_face_only_merging(monkeypatch):
    """The edge pass is what throws; dropping it keeps the valuable half."""
    solid = two_boxes_sharing_a_face()
    real = occ.unify_same_domain
    seen = []

    def refuse_edge_merging(shape, unify_edges=True):
        seen.append(unify_edges)
        if unify_edges:
            raise Standard_Failure("Courbes non jointives")
        return real(shape, unify_edges=False)

    monkeypatch.setattr(occ, "unify_same_domain", refuse_edge_merging)
    merged, ran = _unify_one(solid)
    faces, _ = occ.count_subshapes(merged)
    assert seen == [True, False]
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

    def always_refuse(shape, unify_edges=True):
        raise Standard_Failure("Courbes non jointives")

    monkeypatch.setattr(occ, "unify_same_domain", always_refuse)
    out, stats = _unify([solid])

    assert out[0] is solid
    assert stats["unmerged_solids"] == 1
    assert stats["output_faces"] == stats["faces_before"]
    assert stats["volume_drift"] == pytest.approx(0.0, abs=UNIFY_VOLUME_TOL)
    assert occ.volume(out[0]) == pytest.approx(2000.0)
