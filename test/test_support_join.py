"""Supports sharing a face or volume become one body (docs/algorithm.md §14.3 S6)."""

from __future__ import annotations

import pytest


def _box(origin, size):
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt

    return BRepPrimAPI_MakeBox(gp_Pnt(*origin), *size).Shape()


def _solid(shape):
    from latticegen2 import occ

    return occ.solids(shape)[0]


def test_supports_sharing_a_face_are_one_body_and_its_faces_merge():
    from latticegen2 import occ
    from latticegen2.support import join

    a = _solid(_box((0, 0, 0), (2, 2, 2)))
    b = _solid(_box((2, 0, 0), (2, 2, 2)))
    got, stats = join.join([a, b])
    assert len(got) == 1
    assert occ.volume(got[0]) == pytest.approx(16.0, rel=1e-12)
    assert occ.is_valid(got[0])
    # The shared face is gone and the coplanar ones are merged: one box.
    assert len(occ.faces(got[0])) == 6
    assert stats["support_joined_groups"] == 1
    assert stats["support_overlapping_pairs"] == 0


def test_overlapping_supports_are_fused_to_their_union():
    from latticegen2 import occ
    from latticegen2.support import join

    a = _solid(_box((0, 0, 0), (3, 2, 2)))
    b = _solid(_box((2, 0, 0), (3, 2, 2)))                 # 1 mm into `a`
    got, stats = join.join([a, b])
    assert len(got) == 1
    assert occ.volume(got[0]) == pytest.approx(20.0, rel=1e-12)
    assert stats["support_overlapping_pairs"] == 1


def test_supports_meeting_along_an_edge_stay_separate():
    """Neighbouring cells across a parting meet in a V from one edge: joined,
    they would not be a manifold solid."""
    from latticegen2 import occ
    from latticegen2.support import join

    a = _solid(_box((0, 0, 0), (2, 2, 2)))
    b = _solid(_box((2, 2, 0), (2, 2, 2)))
    got, stats = join.join([a, b])
    assert len(got) == 2
    assert sorted(occ.volume(s) for s in got) == pytest.approx([8.0, 8.0], rel=1e-12)
    assert stats["support_joined_groups"] == 0


def test_supports_apart_are_left_as_they_are():
    from latticegen2.support import join

    a = _solid(_box((0, 0, 0), (2, 2, 2)))
    b = _solid(_box((5, 0, 0), (2, 2, 2)))
    got, _stats = join.join([a, b])
    assert len(got) == 2 and got[0].IsSame(a) and got[1].IsSame(b)


def test_a_fuse_that_loses_volume_is_folded_one_at_a_time(monkeypatch):
    from latticegen2 import occ
    from latticegen2.support import join

    real = join._fuse
    calls = []

    def first_call_returns_an_operand(solids):
        calls.append(len(solids))
        if len(calls) == 1:
            return [solids[0]]
        return real(solids)

    monkeypatch.setattr(join, "_fuse", first_call_returns_an_operand)
    boxes = [_solid(_box((2 * k, 0, 0), (2, 2, 2))) for k in range(3)]
    got, _stats = join.join(boxes)
    assert calls == [3, 2, 2]
    assert len(got) == 1 and occ.volume(got[0]) == pytest.approx(24.0, rel=1e-12)


def test_face_sharing_supports_the_kernel_cannot_fuse_stay_as_they_were(monkeypatch):
    from latticegen2.support import join

    monkeypatch.setattr(join, "_fuse", lambda solids: None)
    a = _solid(_box((0, 0, 0), (2, 2, 2)))
    b = _solid(_box((2, 0, 0), (2, 2, 2)))
    got, stats = join.join([a, b])
    assert len(got) == 2 and got[0].IsSame(a) and got[1].IsSame(b)
    assert stats["support_joins_declined"] == 1


def test_overlapping_supports_the_kernel_cannot_fuse_fail_the_run(monkeypatch):
    from latticegen2.errors import ProcessingError
    from latticegen2.support import join

    monkeypatch.setattr(join, "_fuse", lambda solids: None)
    a = _solid(_box((0, 0, 0), (3, 2, 2)))
    b = _solid(_box((2, 0, 0), (3, 2, 2)))
    with pytest.raises(ProcessingError, match="could not be fused"):
        join.join([a, b])
