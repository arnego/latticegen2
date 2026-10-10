"""Print orientation (specification.md §4.5): the three angles, and the frame.

Kernel-free. The property everything else rests on is at the top: zero angles
give *exactly* the identity, which is what lets the pipeline apply the rotation
on every run with no identity special case.
"""

from __future__ import annotations

import numpy as np
import pytest

from latticegen2 import orient
from latticegen2.lattice import (
    COS_THETA,
    THETA,
    candidate_nodes,
    candidate_nodes_of_points,
    cube_edge,
    lattice_params,
    nodes,
    part_name,
    strut_direction,
)

ATOL = 1e-9
TILTS = [(30.0, 0.0, 0.0), (0.0, 45.0, 0.0), (0.0, 0.0, 70.0),
         (35.0, 20.0, -110.0), (180.0, 0.0, 0.0), (90.0, 90.0, 90.0)]


def test_zero_angles_are_exactly_the_identity():
    R = orient.rotation_matrix((0.0, 0.0, 0.0))
    assert np.array_equal(R, np.eye(3))
    assert np.array_equal(orient.build_direction((0.0, 0.0, 0.0)), [0.0, 0.0, 1.0])


def test_unrotated_frame_equals_the_canonical_expressions_element_for_element():
    """No bypass exists, so this is what guarantees the unrotated lattice.

    The expected values are written out the way ``lattice_params`` computed
    them before the orientation existed, and compared with ``==``.
    """
    cc, t = 10.0, 1.5
    lp = lattice_params(cc, t)
    a = cube_edge(cc)
    e = np.array([strut_direction(k) for k in range(3)], dtype=float)
    zhat = np.array([0.0, 0.0, 1.0])
    u = np.array([np.cross(zhat, e[k]) / np.linalg.norm(np.cross(zhat, e[k]))
                  for k in range(3)])
    v = np.array([np.cross(e[k], u[k]) for k in range(3)])
    B = a * np.column_stack([e[0], e[1], e[2]])
    assert np.array_equal(lp.e, e)
    assert np.array_equal(lp.u, u)
    assert np.array_equal(lp.v, v)
    assert np.array_equal(lp.B, B)
    assert np.array_equal(lattice_params(cc, t, (0, 0, 0)).B, B)


@pytest.mark.parametrize("angles", TILTS)
def test_rotation_is_orthonormal_and_proper(angles):
    R = orient.rotation_matrix(angles)
    assert np.allclose(R @ R.T, np.eye(3), atol=1e-12)
    assert np.linalg.det(R) == pytest.approx(1.0, abs=1e-12)


def test_convention_is_extrinsic_x_then_y_then_z():
    # 90 deg about bed X takes the part's +Y to the bed's +Z ...
    R = orient.rotation_matrix((90.0, 0.0, 0.0))
    assert np.allclose(R @ [0, 1, 0], [0, 0, 1], atol=1e-12)
    # ... so the build direction, seen from the part, is its +Y.
    assert np.allclose(orient.build_direction((90.0, 0.0, 0.0)), [0, 1, 0], atol=1e-12)
    # X first, then Z: the part's +Y goes up, and stays up under the Z turn.
    R = orient.rotation_matrix((90.0, 0.0, 90.0))
    assert np.allclose(R @ [0, 1, 0], [0, 0, 1], atol=1e-12)
    assert np.allclose(R @ [1, 0, 0], [0, 1, 0], atol=1e-12)


@pytest.mark.parametrize("angles", TILTS + [(10.0, 90.0, 25.0), (10.0, -90.0, 25.0)])
def test_euler_round_trip_reproduces_the_rotation(angles):
    R = orient.rotation_matrix(angles)
    again = orient.rotation_matrix(orient.euler_from_matrix(R))
    assert np.allclose(again, R, atol=1e-9)


@pytest.mark.parametrize("angles", TILTS)
def test_lattice_identities_hold_about_the_build_direction(angles):
    """docs/algorithm.md §2.3, with the bed's normal standing where Z stood."""
    lp = lattice_params(10.0, 1.5, angles)
    b = lp.b
    assert np.allclose(np.linalg.norm(lp.e, axis=1), 1.0, atol=ATOL)
    assert np.allclose(lp.e @ b, COS_THETA, atol=ATOL)                 # recline
    assert np.allclose(np.degrees(np.arccos(lp.e @ b)), np.degrees(THETA), atol=1e-7)
    assert np.allclose(lp.e @ lp.e.T, np.eye(3), atol=ATOL)            # orthogonal
    assert np.allclose(lp.u @ b, 0.0, atol=ATOL)                       # u horizontal
    for k in range(3):
        # v lies in the vertical plane containing the strut axis.
        assert abs(float(np.dot(lp.v[k], np.cross(lp.e[k], b)))) < ATOL
    # In-plane neighbours are cc apart and level; the body diagonal is vertical.
    p, q = nodes(lp, np.array([[1, 0, 0], [0, 1, 0]]))
    assert np.linalg.norm(p - q) == pytest.approx(lp.cc, abs=ATOL)
    assert float((p - q) @ b) == pytest.approx(0.0, abs=ATOL)
    diag = nodes(lp, np.array([[1, 1, 1]]))[0]
    assert np.allclose(diag, lp.a * np.sqrt(3.0) * b, atol=ATOL)


@pytest.mark.parametrize("angles", [(0.0, 0.0, 0.0)] + TILTS)
def test_tight_candidates_are_a_subset_that_keeps_every_reachable_node(angles):
    """Bounding the surface points drops only nodes that cannot reach the body.

    A sphere's surface points stand in for the mesh. Every node whose junction
    could touch the sphere (within ``a/2`` of it) must survive, in the same
    relative order the box enumeration lists it in.
    """
    lp = lattice_params(10.0, 1.5, angles)
    rng = np.random.default_rng(7)
    x = rng.normal(size=(4000, 3))
    radius = 23.0
    pts = radius * x / np.linalg.norm(x, axis=1, keepdims=True) + [5.0, -3.0, 11.0]
    tight = candidate_nodes_of_points(lp, pts)
    box = candidate_nodes(lp, pts.min(axis=0), pts.max(axis=0))

    tight_set = {tuple(r) for r in tight.tolist()}
    kept = [tuple(r) for r in box.tolist() if tuple(r) in tight_set]
    reach = np.linalg.norm(nodes(lp, box) - [5.0, -3.0, 11.0], axis=1) <= radius + lp.a / 2
    assert all(tuple(r) in tight_set for r in box[reach].tolist())
    # Same relative order: the survivors of the box list, read in order, are a
    # subsequence of the tight list.
    order = {n: i for i, n in enumerate(map(tuple, tight.tolist()))}
    positions = [order[n] for n in kept]
    assert positions == sorted(positions)


def test_names_carry_rotation_and_support_only_when_they_apply():
    assert part_name("/x/ball.step", 20, 4) == "ball+lattice+cc20+t4"
    assert part_name("/x/ball.step", 20, 4, (0.0, 0.0, 0.0)) == "ball+lattice+cc20+t4"
    assert (part_name("/x/ball.step", 20, 4, (30.0, 0.0, -22.5))
            == "ball+lattice+cc20+t4+rotx30y0z-22.5")
    assert (part_name("/x/ball.step", 20, 4, (30.0, 0.0, 0.0), 60.0)
            == "ball+lattice+cc20+t4+rotx30y0z0+sup60")
    assert part_name("/x/ball.step", 20, 4, overhang=45.0) == "ball+lattice+cc20+t4+sup45"
