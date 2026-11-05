"""The E2E harness's own containment check (docs/testing.md, specification.md §6.2).

`tools/` is normally tested by running it — `tools/e2e.py` *is* the gate for the
harness. This file exists for the one part of it whose failure mode is a
plausible number rather than an exception, and which therefore cannot be caught
by reading the output of a passing run.

`material_outside` reported **354,733 mm³** of material outside the input body —
the entire lattice — on the `cc=12, t=2.5` output of `TD_HX_rehearsal_test`,
whose every junction was *built* by intersecting with that body. The boolean did
not decline: it reported `IsDone`, `HasModified` and `HasGenerated`, and returned
43,672 faces where 43,530 went in. A lattice trimmed from a body has a large
share of its faces lying exactly *on* that body's surface, which is the classic
ill-conditioned case for a cut.

What these tests pin is the shape of the repair rather than that one part: a
trustworthy answer is *exactly* zero, so any remainder is contradicted against a
boolean-free containment check, and a disagreement is reported as unmeasured
rather than as a pass. The real part is too large to carry in a unit test — it is
covered by `tools/e2e.py` and by the measurements in specification.md §11 — so
these use boxes, where "contained" and "sticking out" are known by construction.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tools"))

import verify_geometry as vg  # noqa: E402

from latticegen2 import occ  # noqa: E402


def box(dx: float, dy: float, dz: float, at=(0.0, 0.0, 0.0)):
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt

    return BRepPrimAPI_MakeBox(gp_Pnt(*at), dx, dy, dz).Solid()


@pytest.fixture(scope="module")
def outer():
    return box(100.0, 100.0, 100.0)


def test_a_contained_solid_is_measured_exactly_and_at_zero(outer, tmp_path):
    """The ordinary case must stay exact — the repair must not have turned every
    measurement into the weaker fallback."""
    inner = box(10.0, 10.0, 10.0, at=(20.0, 20.0, 20.0))
    cand, body = tmp_path / "cand.step", tmp_path / "body.step"
    occ.write_step(inner, str(cand), "cand")
    occ.write_step(outer, str(body), "body")

    result = vg.material_outside(str(cand), str(body))

    assert result["exact"] is True
    assert result["unmeasured"] == []
    assert result["volume_mm3"] == 0.0
    assert all(s["trusted"] for s in result["per_solid"])


def test_a_solid_sticking_out_is_reported_not_swallowed(outer, tmp_path):
    """The contradiction path must never be able to hide real material outside.

    A boolean-free containment check sees a protrusion too, so whichever route
    the result takes, the scenario fails rather than passes.
    """
    poking = box(10.0, 10.0, 10.0, at=(95.0, 20.0, 20.0))   # 5 mm of 10 mm is out
    cand, body = tmp_path / "cand.step", tmp_path / "body.step"
    occ.write_step(poking, str(cand), "cand")
    occ.write_step(outer, str(body), "body")

    result = vg.material_outside(str(cand), str(body))

    if result["exact"]:
        assert result["volume_mm3"] == pytest.approx(500.0, rel=1e-6)
    else:
        # Taken as unmeasured, the containment evidence must still condemn it.
        assert len(result["unmeasured"]) == 1
        assert result["unmeasured"][0]["containment"]["contained"] is False


def test_a_negative_cut_remainder_is_not_believed(outer, tmp_path, monkeypatch):
    """A cut can come back with a small negative volume -- measured on a
    support body of the rehearsal part. That is no measurement of material
    outside, so it goes to the containment check like any other non-zero
    remainder, instead of failing a contained solid by its absolute value."""
    inner = box(10.0, 10.0, 10.0, at=(20.0, 20.0, 20.0))
    cand, body = tmp_path / "cand.step", tmp_path / "body.step"
    occ.write_step(inner, str(cand), "cand")
    occ.write_step(outer, str(body), "body")
    monkeypatch.setattr(vg, "_cut_volume", lambda a, b: -1.5e-3)

    result = vg.material_outside(str(cand), str(body))

    assert result["exact"] is False
    assert len(result["unmeasured"]) == 1
    assert result["unmeasured"][0]["containment"]["contained"] is True


def test_a_bare_sliver_is_told_from_a_bare_ceiling(tmp_path):
    """Coverage of a flat roof: material stopping just short of it is a sliver
    and does not fail, stopping well short does, and a sample that probes out of
    the body is never counted against the output."""
    body = box(20.0, 20.0, 20.0)
    body_path = tmp_path / "body.step"
    occ.write_step(body, str(body_path), "body")

    def coverage(gap):
        cand = tmp_path / f"cand-{gap}.step"
        occ.write_step(box(20.0 - gap, 20.0, 20.0), str(cand), "cand")
        return vg.overhang_coverage(str(cand), str(body_path), 10.0, 1.5, (0.0, 0.0, 0.0), 1.5,
                                    pitch=0.5)

    whole = coverage(0.0)
    assert whole["uncovered"] == 0 and whole["slivers"] == 0

    near = coverage(vg.COVERAGE_SLIVER_TOL / 2.0)
    assert near["uncovered"] == 0 and near["slivers"] > 0

    far = coverage(2.0)
    assert far["uncovered"] > 0


def test_containment_is_boolean_free_and_agrees_with_construction(outer):
    """`surface_points_outside` is the arbiter, so it has to be right on its own."""
    inside = vg.surface_points_outside(box(10.0, 10.0, 10.0, at=(20.0, 20.0, 20.0)), outer)
    assert inside["sampled"] > 0
    assert inside["outside"] == 0
    assert inside["contained"] is True
    assert inside["cleared_at_mm"] == vg.CONTAINMENT_SWEEP[0], (
        "a solid well clear of the boundary must clear at the tightest tolerance"
    )

    out = vg.surface_points_outside(box(10.0, 10.0, 10.0, at=(95.0, 20.0, 20.0)), outer)
    assert out["outside"] > 0
    assert out["contained"] is False
    assert out["cleared_at_mm"] is None


def test_a_face_lying_on_the_boundary_is_a_tie_not_a_protrusion(outer):
    """The case the whole tolerance sweep exists for.

    A trimmed lattice's faces lie exactly *on* the input surface, so its points
    there are ties. They must not read as material outside — otherwise the
    fallback would condemn every part it was built to rescue.
    """
    flush = box(10.0, 10.0, 10.0, at=(90.0, 20.0, 20.0))    # right face coincident
    ev = vg.surface_points_outside(flush, outer)

    assert ev["contained"] is True, ev
    assert ev["cleared_at_mm"] <= vg.CONTAINMENT_TOL


def test_a_real_protrusion_is_not_swallowed_by_the_distance_contradiction(outer):
    """The load-bearing half: the contradiction must not become an off switch.

    `surface_points_outside` now contradicts a point the tolerance sweep cannot
    clear by measuring its exact distance to the body, because the sweep tops
    out at 2e-03 mm and a 1e-09 mm tie reads there exactly like a 1 mm
    protrusion. A contradiction that only ever exonerates would be worse than
    no check at all (docs/testing.md, G10), so this pins the other direction:
    a box genuinely sticking out measures a real distance and is still
    condemned.
    """
    out = vg.surface_points_outside(box(10.0, 10.0, 10.0, at=(95.0, 20.0, 20.0)), outer)
    assert out["contained"] is False
    assert out["worst_distance_mm"] > vg.PROTRUSION_MIN_MM
    assert out["worst_distance_at"] is not None, "and it says where"


def test_the_distance_contradiction_exonerates_only_points_actually_on_the_body(outer):
    """And the direction it was built for, stated as the quantity rather than a count.

    A face flush with the boundary has its points *on* the body, so their exact
    distance to it is zero however the classifier votes about which side they
    are on. That is what makes the bar tight enough to be meaningless as a
    calibration: the ties measured on `SpiralTest.step` sit at 0.000000 mm
    (median 9.95e-10) against protrusions of 0.239 to 0.771 mm.
    """
    flush = box(10.0, 10.0, 10.0, at=(90.0, 20.0, 20.0))
    ev = vg.surface_points_outside(flush, outer)
    assert ev["contained"] is True
    assert ev["worst_distance_mm"] <= vg.PROTRUSION_MIN_MM


# --- the final print's free overhangs (tools/verify_geometry.py) -------------
#
# A support node is a cube standing on its corner: its faces lean 35.26 deg
# from vertical and its edges exactly the rod angle, so the only thing about it
# that can hang is its lowest corner. Each case below is one detector's true
# positive in free space and the same solid rooted on a floor, the lattice or
# another support, which must pass.

import numpy as np  # noqa: E402

UP = np.array([0.0, 0.0, 1.0])


def _moved(shape, rotate=None, to=(0.0, 0.0, 0.0)):
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
    from OCP.gp import gp_Ax1, gp_Dir, gp_Pnt, gp_Trsf, gp_Vec

    if rotate is not None:
        axis, angle = rotate
        r = gp_Trsf()
        r.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(*axis)), angle)
        shape = BRepBuilderAPI_Transform(shape, r, True).Shape()
    m = gp_Trsf()
    m.SetTranslation(gp_Vec(*to))
    return occ.solids(BRepBuilderAPI_Transform(shape, m, True).Shape())[0]


def _node(at):
    """A 2 mm cube with its (1,1,1) diagonal vertical, lowest corner at ``at``."""
    axis = np.cross([1.0, 1.0, 1.0], UP)
    return _moved(box(2.0, 2.0, 2.0), (axis, float(np.arccos(1.0 / np.sqrt(3.0)))), at)


def _keel(at):
    """A 2 mm cube rolled 45 deg about x: a horizontal edge at the bottom, its
    two faces 45 deg from vertical, inside the rod angle."""
    return _moved(box(2.0, 2.0, 2.0), ((1.0, 0.0, 0.0), np.pi / 4.0), at)


def _free(supports, lattice=(), body=None):
    got = vg.free_support_overhangs(list(supports), list(lattice),
                                    body if body is not None else box(100.0, 100.0, 100.0), UP)
    return got["candidates"], got["free"]


def test_a_support_node_hanging_in_free_space_is_a_free_tip():
    cand, free = _free([_node((50.0, 50.0, 50.0))])
    assert cand == {"face": 0, "edge": 0, "vertex": 1}
    assert free == {"face": 0, "edge": 0, "vertex": 1}


def test_a_support_node_standing_on_a_floor_the_lattice_or_a_support_is_rooted():
    on_floor = _node((50.0, 50.0, 0.0))
    assert _free([on_floor])[1] == {"face": 0, "edge": 0, "vertex": 0}
    strut = box(4.0, 4.0, 4.0, at=(48.0, 48.0, 46.0))          # top at z = 50
    assert _free([_node((50.0, 50.0, 50.0))], [strut])[1]["vertex"] == 0
    below = box(4.0, 4.0, 4.0, at=(20.0, 20.0, 26.0))
    cand, free = _free([_node((22.0, 22.0, 30.0)), below])
    assert cand["vertex"] == 1 and free["vertex"] == 0         # the slab's own corners rise flat


def test_a_horizontal_edge_in_free_space_is_a_free_keel_and_rooted_on_a_floor():
    cand, free = _free([_keel((50.0, 50.0, 50.0))])
    assert free["face"] == 0 and free["vertex"] == 0
    assert cand["edge"] >= 5 and free["edge"] == cand["edge"]   # 2 mm at 0.5 mm
    assert _free([_keel((50.0, 50.0, 0.0))])[1] == {"face": 0, "edge": 0, "vertex": 0}


def test_a_flat_underside_hangs_unless_it_rests_on_something():
    cand, free = _free([box(4.0, 4.0, 1.0, at=(50.0, 50.0, 50.0))])
    assert cand["face"] > 0 and free["face"] == cand["face"]
    assert free["edge"] == 0 and free["vertex"] == 0
    cand, free = _free([box(4.0, 4.0, 1.0, at=(50.0, 50.0, 0.0))])
    assert cand["face"] > 0 and free == {"face": 0, "edge": 0, "vertex": 0}
    strut = box(10.0, 10.0, 2.0, at=(47.0, 47.0, 48.0))
    cand, free = _free([box(4.0, 4.0, 1.0, at=(50.0, 50.0, 50.0))], [strut])
    assert free == {"face": 0, "edge": 0, "vertex": 0}


def test_a_support_against_a_tilted_wall_is_rooted_on_it():
    """Where a support's rod-plane face meets the wall, their edge can run
    flatter than the rod angle; it lies on the wall and is rooted there. The
    floor here is tilted 30 deg and cuts the node's lowest corner off, so the
    trimmed node rests on it with a face, keels and corners, none of them
    free."""
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Common

    body = _moved(box(100.0, 100.0, 100.0), ((0.0, 1.0, 0.0), np.radians(30.0)), (0, 0, 0))
    node = _moved(box(20.0, 20.0, 20.0), (np.cross([1.0, 1.0, 1.0], UP),
                                          float(np.arccos(1.0 / np.sqrt(3.0)))), (30.0, 50.0, -20.0))
    trimmed = occ.solids(BRepAlgoAPI_Common(node, body).Shape())[0]
    cand, free = _free([trimmed], body=body)
    assert cand["face"] > 0
    assert free == {"face": 0, "edge": 0, "vertex": 0}
    # The control: the same solid with the floor somewhere else hangs whole.
    assert _free([trimmed])[1] == cand


def test_a_saddle_whose_edges_all_rise_is_not_a_tip():
    """Three unit boxes in the octants (+,-,-), (+,-,+) and (+,+,-) round a
    point, with up along (1,1,1): every edge at that point runs along a +e_k,
    yet the faces there hold material below it. The one tip is the union's
    lowest corner, (0,-1,-1) from it."""
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse
    from OCP.TopTools import TopTools_ListOfShape

    o = np.array([50.0, 50.0, 50.0])
    parts = [box(1.0, 1.0, 1.0, at=tuple(o + d)) for d in
             ((0, -1, -1), (0, -1, 0), (0, 0, -1))]
    fuse = BRepAlgoAPI_Fuse()
    args, tools = TopTools_ListOfShape(), TopTools_ListOfShape()
    args.Append(parts[0])
    for p in parts[1:]:
        tools.Append(p)
    fuse.SetArguments(args)
    fuse.SetTools(tools)
    fuse.Build()
    # Unified, as the generator's staircases are: split coplanar faces would
    # leave edges running down from the saddle, and the edge test alone would
    # pass.
    solid = occ.unify_same_domain(occ.solids(fuse.Shape())[0])
    up = np.ones(3) / np.sqrt(3.0)
    occ.mesh_shape(solid, vg.SUPPORT_MESH_DEFLECTION)
    tips = vg._material_tip_samples(solid, -up)
    assert len(tips) == 1
    assert np.allclose(tips[0], o + (0.0, -1.0, -1.0))
