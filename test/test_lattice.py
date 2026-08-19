"""Lattice mathematics — docs/algorithm.md §2.3's verified identities.

Includes the mutual-orthogonality identity the whole instancing architecture
depends on.
"""

import os

import numpy as np
import pytest

from latticegen2.lattice import (
    HALF_STRUTS,
    OPPOSITE_HALF,
    THETA,
    candidate_nodes,
    cube_edge,
    format_param,
    half_strut_offset,
    index_range,
    lattice_params,
    neighbor_step,
    node,
    part_name,
    profile_vertices,
    strut_direction,
)

ATOL = 1e-9


def test_strut_directions_are_unit_vectors():
    for k in range(3):
        assert np.linalg.norm(strut_direction(k)) == pytest.approx(1.0, abs=ATOL)


def test_recline_angle_matches_theta():
    z = np.array([0.0, 0.0, 1.0])
    for k in range(3):
        assert np.arccos(np.dot(strut_direction(k), z)) == pytest.approx(THETA, abs=ATOL)


def test_theta_is_about_54_7_degrees():
    assert np.degrees(THETA) == pytest.approx(54.7356103172, abs=1e-9)


def test_horizontal_projections_are_120_degrees_apart():
    proj = [strut_direction(k)[:2] for k in range(3)]
    for a, b in ((0, 1), (1, 2), (2, 0)):
        cos = np.dot(proj[a], proj[b]) / (np.linalg.norm(proj[a]) * np.linalg.norm(proj[b]))
        assert np.degrees(np.arccos(cos)) == pytest.approx(120.0, abs=1e-9)


def test_strut_directions_are_mutually_orthogonal():
    """The identity the whole junction-instancing architecture rests on.

    Mutual orthogonality is what confines strut-strut overlap to inside a single
    junction, so two adjacent junctions meet on a shared face instead of
    overlapping in volume.
    """
    for a, b in ((0, 1), (0, 2), (1, 2)):
        assert np.dot(strut_direction(a), strut_direction(b)) == pytest.approx(0.0, abs=ATOL)


def test_in_plane_neighbour_spacing_equals_cc():
    lp = lattice_params(10.0, 1.5)
    p = node(lp, (1, 0, 0))
    q = node(lp, (0, 1, 0))
    assert np.linalg.norm(p - q) == pytest.approx(lp.cc, abs=ATOL)
    assert p[2] == pytest.approx(q[2], abs=ATOL), "cc must be a pure XY-plane distance"


def test_vertical_space_diagonal():
    lp = lattice_params(10.0, 1.5)
    d = node(lp, (1, 1, 1)) - node(lp, (0, 0, 0))
    assert d[0] == pytest.approx(0.0, abs=ATOL)
    assert d[1] == pytest.approx(0.0, abs=ATOL)
    assert d[2] == pytest.approx(lp.a * np.sqrt(3), abs=ATOL)


def test_cube_edge_and_circumradius():
    lp = lattice_params(10.0, 1.5)
    assert lp.a == pytest.approx(cube_edge(10.0))
    assert lp.a == pytest.approx(10.0 / np.sqrt(2))
    assert lp.r == pytest.approx(1.5 / np.sqrt(2))


def test_profile_frame_orientation():
    """u is horizontal; v lies in the vertical plane containing the strut axis."""
    lp = lattice_params(10.0, 1.5)
    for k in range(3):
        assert lp.u[k][2] == pytest.approx(0.0, abs=ATOL)
        assert np.dot(lp.u[k], lp.e[k]) == pytest.approx(0.0, abs=ATOL)
        assert np.dot(lp.v[k], lp.e[k]) == pytest.approx(0.0, abs=ATOL)
        assert np.linalg.norm(np.cross(lp.u[k], lp.v[k]) - lp.e[k]) == pytest.approx(0.0, abs=ATOL)


def test_profile_is_a_square_of_side_t():
    lp = lattice_params(10.0, 1.5)
    verts = profile_vertices(lp, np.zeros(3), 0)
    sides = [np.linalg.norm(verts[i] - verts[(i + 1) % 4]) for i in range(4)]
    for s in sides:
        assert s == pytest.approx(lp.t, abs=ATOL)
    # The diagonals are t*sqrt(2) and one of them is horizontal.
    assert np.linalg.norm(verts[0] - verts[2]) == pytest.approx(lp.t * np.sqrt(2), abs=ATOL)
    assert (verts[0] - verts[2])[2] == pytest.approx(0.0, abs=ATOL)


def test_half_struts_pair_across_shared_caps():
    """Half-strut h of node n and its opposite at n+step(h) share one cap quad."""
    lp = lattice_params(10.0, 1.5)
    n = np.array([2, -1, 3])
    for h in range(6):
        step = neighbor_step(h)
        opp = OPPOSITE_HALF[h]
        cap_here = node(lp, n) + half_strut_offset(lp, h)
        cap_there = node(lp, n + step) + half_strut_offset(lp, opp)
        assert np.linalg.norm(cap_here - cap_there) == pytest.approx(0.0, abs=ATOL)


def test_half_strut_offsets_are_half_a_cell():
    lp = lattice_params(10.0, 1.5)
    for h in range(6):
        assert np.linalg.norm(half_strut_offset(lp, h)) == pytest.approx(lp.a / 2, abs=ATOL)
    assert len(HALF_STRUTS) == 6
    assert len(set(HALF_STRUTS)) == 6


def test_index_range_covers_the_box_with_padding():
    lp = lattice_params(10.0, 1.5)
    lo = np.array([-20.0, -20.0, -20.0])
    hi = np.array([20.0, 20.0, 20.0])
    lo_idx, hi_idx = index_range(lp, lo, hi)
    assert np.all(lo_idx < hi_idx)
    # Every corner of the box maps inside the index range.
    for x in (lo[0], hi[0]):
        for y in (lo[1], hi[1]):
            for z in (lo[2], hi[2]):
                idx = np.linalg.solve(lp.B, np.array([x, y, z]))
                assert np.all(idx >= lo_idx) and np.all(idx <= hi_idx)


def test_candidate_nodes_matches_the_index_range():
    lp = lattice_params(20.0, 4.0)
    lo, hi = np.array([-10.0, -10.0, -10.0]), np.array([10.0, 10.0, 10.0])
    lo_idx, hi_idx = index_range(lp, lo, hi)
    cand = candidate_nodes(lp, lo, hi)
    expected = np.prod(hi_idx - lo_idx + 1)
    assert len(cand) == expected
    assert len(np.unique(cand, axis=0)) == expected


@pytest.mark.parametrize(
    "value,expected", [(20.0, "20"), (1.5, "1.5"), (0.4, "0.4"), (10.0, "10"), (2.25, "2.25")]
)
def test_format_param_drops_trailing_zeros(value, expected):
    assert format_param(value) == expected


def test_part_name_uses_the_plus_convention():
    assert part_name("/x/y/ball.step", 2.5, 0.4) == "ball+lattice+cc2.5+t0.4"


def test_part_name_carries_the_same_components_as_the_default_file_name():
    """The two names differ in punctuation alone (specification.md §5)."""
    from latticegen2.cli import resolve_output_paths

    step_path, _ = resolve_output_paths("/x/y/ball.step", None, 2.5, 0.4)
    stem = os.path.splitext(os.path.basename(step_path))[0]
    assert stem == "ball-lattice-cc2.5t0.4"
    assert part_name("/x/y/ball.step", 2.5, 0.4).split("+") == [
        "ball",
        "lattice",
        "cc2.5",
        "t0.4",
    ]
