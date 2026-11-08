"""Sub-micron caps on lattice bodies (docs/algorithm.md §9).

The cap docs/algorithm.md §14.3 S6 removes from support bodies can form on the
lattice too: a strut's corner grazing the input's wall within a micron leaves
its tip as a few faces of ~1e-8 mm^2 closing one inner wire of the wall face,
and §9's export-truth gate refuses the run. `simplify` removes it, on the
finished solid, checking the repair locally and riding the validity check it
already runs on every unified solid.

The fixture is a cube's corner pressed into a cylinder's crest. Lattice struts
run along three mutually perpendicular directions (docs/algorithm.md §2.3), so
the corner a junction presents to a wall is exactly that trihedral.
"""

from __future__ import annotations

import numpy as np

from latticegen2 import occ
from latticegen2.parallel import WorkerPool
from latticegen2.pipeline import _simplify_one, _unify
from latticegen2.support import frustum

from test_support_micro_caps import _dimpled_cylinder

CAP_AT = (20.3, 0.0, 20.0)


def test_simplify_removes_a_cap_from_a_lattice_solid():
    s = _dimpled_cylinder(1e-4)
    out, ran, degraded, removed, refused = _simplify_one(s)
    assert ran and not degraded
    assert refused == []
    assert len(removed) == 1 and np.allclose(removed[0], CAP_AT, atol=1e-3)
    assert len(occ.faces(out)) == 3
    assert occ.is_valid(out)
    assert frustum.edges_not_used_twice(out) == 0
    assert abs(occ.volume(out) - occ.volume(s)) < 1e-6


def test_a_solid_without_a_cap_is_unchanged_by_the_repair():
    s = _dimpled_cylinder(0.5)
    out, ran, degraded, removed, refused = _simplify_one(s)
    assert removed == [] and refused == []
    assert len(occ.faces(out)) == len(occ.faces(s))


def test_the_local_check_accepts_the_repair_and_names_a_wrong_one():
    s = _dimpled_cylinder(1e-4)
    hosts: list = []
    fixed, where = occ.remove_micro_caps(s, hosts=hosts)
    assert len(where) == 1 and len(hosts) == 1
    assert occ.micro_cap_repair_fault(s, fixed, hosts) is None
    # Claiming one cap face fewer than was removed: the count no longer adds up.
    host, wire, n = hosts[0]
    assert "face count" in occ.micro_cap_repair_fault(s, fixed, [(host, wire, n - 1)])
    # The input itself is not the input less a cap.
    assert occ.micro_cap_repair_fault(s, s, hosts) is not None


def test_a_repair_failing_the_local_check_keeps_the_cap(monkeypatch):
    s = _dimpled_cylinder(1e-4)
    monkeypatch.setattr(occ, "micro_cap_repair_fault", lambda *a, **k: "refused")
    out, ran, degraded, removed, refused = _simplify_one(s)
    assert removed == []
    assert len(refused) == 1 and np.allclose(refused[0], CAP_AT, atol=1e-3)
    assert len(occ.faces(out)) == 6
    assert not degraded


def test_a_repair_the_whole_solid_check_refuses_falls_back_to_the_unrepaired(monkeypatch):
    """A repaired solid `BRepCheck_Analyzer` refuses costs the run nothing: the
    unrepaired one is checked in its place and kept if it holds."""
    s = _dimpled_cylinder(1e-4)
    real = occ.is_valid

    def refuses_the_repaired(shape, parallel=True):
        if len(occ.faces(shape)) == 3:
            return False
        return real(shape, parallel=parallel)

    monkeypatch.setattr(occ, "is_valid", refuses_the_repaired)
    out, ran, degraded, removed, refused = _simplify_one(s)
    assert removed == []
    assert len(refused) == 1
    assert len(occ.faces(out)) == 6
    assert not degraded


def test_simplify_reports_caps_on_both_paths(tmp_path):
    def solids():
        return [_dimpled_cylinder(1e-4), _dimpled_cylinder(0.5)]

    serial_out, serial = _unify(solids())
    with WorkerPool(2) as pool:
        parallel_out, parallel = _unify(solids(), pool=pool, tmpdir=str(tmp_path))
    assert parallel["max_worker_rss"] > 0
    for stats in (serial, parallel):
        assert len(stats["caps_removed"]) == 1
        assert np.allclose(stats["caps_removed"][0], CAP_AT, atol=1e-3)
        assert stats["caps_refused"] == []
    assert [len(occ.faces(s)) for s in parallel_out] == [
        len(occ.faces(s)) for s in serial_out]
