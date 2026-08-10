"""The junction graph and the floating-body rule (specification.md §5).

The rule is deliberately asymmetric: a body is only ever dropped once it is
*proven* disconnected. A small fragment still attached to the rest is kept, and
there is no case where connectivity is unknown — that was the whole point of
answering it by graph rather than by attempting a boolean.
"""

import numpy as np
import pytest

from latticegen2.connect import build_components
from latticegen2.lattice import OPPOSITE_HALF


class FakePiece:
    """Stands in for a trimmed boundary junction piece."""

    def __init__(self, node, caps, volume):
        self.node = node
        self.caps = frozenset(caps)
        self.volume = volume


def nodes(*rows):
    return np.array(rows, dtype=np.int64) if rows else np.empty((0, 3), dtype=np.int64)


def test_a_chain_of_interior_nodes_is_one_component():
    comps = build_components(nodes((0, 0, 0), (1, 0, 0), (2, 0, 0)), 10.0, [])
    assert len(comps.volumes) == 1
    assert comps.volumes[0] == pytest.approx(30.0)


def test_two_separated_clusters_stay_separate():
    comps = build_components(nodes((0, 0, 0), (1, 0, 0), (10, 0, 0), (11, 0, 0)), 10.0, [])
    assert len(comps.volumes) == 2
    assert sorted(comps.volumes.values()) == [pytest.approx(20.0), pytest.approx(20.0)]


def test_diagonal_neighbours_are_not_connected():
    """Junctions join only along a shared strut axis, never across a diagonal."""
    comps = build_components(nodes((0, 0, 0), (1, 1, 0)), 10.0, [])
    assert len(comps.volumes) == 2


def test_a_boundary_piece_joins_the_interior_through_its_cap():
    piece = FakePiece((1, 0, 0), caps=[OPPOSITE_HALF[0]], volume=2.0)
    comps = build_components(nodes((0, 0, 0)), 10.0, [piece])
    assert len(comps.volumes) == 1
    assert comps.volumes[0] == pytest.approx(12.0)


def test_a_boundary_piece_without_the_matching_cap_stays_separate():
    """A cap the trim removed cannot carry a connection."""
    piece = FakePiece((1, 0, 0), caps=[], volume=2.0)
    comps = build_components(nodes((0, 0, 0)), 10.0, [piece])
    assert len(comps.volumes) == 2


def test_two_boundary_pieces_join_across_a_shared_cap():
    a = FakePiece((0, 0, 0), caps=[0], volume=3.0)
    b = FakePiece((1, 0, 0), caps=[OPPOSITE_HALF[0]], volume=4.0)
    comps = build_components(nodes(), 0.0, [a, b])
    assert len(comps.volumes) == 1
    assert comps.volumes[0] == pytest.approx(7.0)


def test_one_junction_split_into_pieces_contributes_one_vertex_each():
    a = FakePiece((5, 5, 5), caps=[0], volume=1.0)
    b = FakePiece((5, 5, 5), caps=[1], volume=2.0)
    comps = build_components(nodes(), 0.0, [a, b])
    assert len(comps.vertices) == 2
    assert len(comps.volumes) == 2, "pieces of one junction are not joined by sharing a node"


def test_small_floating_component_is_dropped_and_a_large_one_is_not():
    big = FakePiece((0, 0, 0), caps=[], volume=100.0)
    small = FakePiece((9, 9, 9), caps=[], volume=0.5)
    comps = build_components(nodes(), 0.0, [big, small])
    dropped = comps.dropped(threshold=1.0)
    assert len(dropped) == 1
    assert comps.volumes[next(iter(dropped))] == pytest.approx(0.5)


def test_a_small_but_connected_piece_is_never_dropped():
    """specification.md §5: only *floating* sub-threshold bodies may be removed."""
    tiny = FakePiece((1, 0, 0), caps=[OPPOSITE_HALF[0]], volume=0.01)
    comps = build_components(nodes((0, 0, 0)), 100.0, [tiny])
    assert comps.dropped(threshold=1.0) == set()
    assert len(comps.volumes) == 1


def test_component_membership_covers_every_vertex():
    pieces = [FakePiece((0, 0, 0), caps=[0], volume=1.0),
              FakePiece((1, 0, 0), caps=[OPPOSITE_HALF[0]], volume=1.0)]
    comps = build_components(nodes((5, 5, 5)), 2.0, pieces)
    assert sum(len(m) for m in comps.members.values()) == len(comps.vertices)
    assert set(comps.labels) == set(comps.volumes)
