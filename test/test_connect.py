"""The junction graph and the floating-body rule (specification.md §5).

The rule is deliberately asymmetric: a body is only ever dropped once it is
*proven* disconnected. A small fragment still attached to the rest is kept, and
there is no case where connectivity is unknown — that was the whole point of
answering it by graph rather than by attempting a boolean.
"""

import numpy as np
import pytest

from latticegen2.connect import build_components
from latticegen2.errors import ProcessingError
from latticegen2.lattice import OPPOSITE_HALF, neighbor_step


class FakePiece:
    """Stands in for a trimmed boundary junction piece.

    ``caps`` is given as bare half-strut ids, relative to ``node`` — the
    ordinary case, and all these tests need — but stored the same way
    production ``BoundaryPiece.caps`` now is, as ``(node, h)`` pairs, since a
    piece :func:`latticegen2.boundary.fuse_disagreeing_pairs` has merged can
    hold caps belonging to more than one node (docs/algorithm.md §7.1).
    """

    def __init__(self, node, caps, volume):
        self.node = node
        self.caps = frozenset((node, h) for h in caps)
        self.volume = volume


def nodes(*rows):
    return np.array(rows, dtype=np.int64) if rows else np.empty((0, 3), dtype=np.int64)


def across(node, h):
    step = neighbor_step(h)
    return (node[0] + int(step[0]), node[1] + int(step[1]), node[2] + int(step[2]))


def interfaces_for(interior, pieces):
    """What ``boundary.resolve_interfaces`` would produce for this fixture.

    Interior nodes present all six caps; a piece presents the caps it gives up.
    A cap is an interface exactly when both sides present it — the same
    symmetric rule, expressed without the geometry kernel.
    """
    present = {(tuple(int(x) for x in n), h) for n in interior for h in range(6)}
    present |= {key for p in pieces for key in p.caps}
    out = set()
    for node, h in present:
        partner = (across(node, h), OPPOSITE_HALF[h])
        if partner in present:
            out.add((node, h))
            out.add(partner)
    return out


def components(interior, volume, pieces):
    """``build_components`` with the interface set its inputs imply."""
    return build_components(interior, volume, pieces, interfaces_for(interior, pieces))


def test_a_chain_of_interior_nodes_is_one_component():
    comps = components(nodes((0, 0, 0), (1, 0, 0), (2, 0, 0)), 10.0, [])
    assert len(comps.volumes) == 1
    assert comps.volumes[0] == pytest.approx(30.0)


def test_two_separated_clusters_stay_separate():
    comps = components(nodes((0, 0, 0), (1, 0, 0), (10, 0, 0), (11, 0, 0)), 10.0, [])
    assert len(comps.volumes) == 2
    assert sorted(comps.volumes.values()) == [pytest.approx(20.0), pytest.approx(20.0)]


def test_diagonal_neighbours_are_not_connected():
    """Junctions join only along a shared strut axis, never across a diagonal."""
    comps = components(nodes((0, 0, 0), (1, 1, 0)), 10.0, [])
    assert len(comps.volumes) == 2


def test_a_boundary_piece_joins_the_interior_through_its_cap():
    piece = FakePiece((1, 0, 0), caps=[OPPOSITE_HALF[0]], volume=2.0)
    comps = components(nodes((0, 0, 0)), 10.0, [piece])
    assert len(comps.volumes) == 1
    assert comps.volumes[0] == pytest.approx(12.0)


def test_a_boundary_piece_without_the_matching_cap_stays_separate():
    """A cap the trim removed cannot carry a connection."""
    piece = FakePiece((1, 0, 0), caps=[], volume=2.0)
    comps = components(nodes((0, 0, 0)), 10.0, [piece])
    assert len(comps.volumes) == 2


def test_two_boundary_pieces_join_across_a_shared_cap():
    a = FakePiece((0, 0, 0), caps=[0], volume=3.0)
    b = FakePiece((1, 0, 0), caps=[OPPOSITE_HALF[0]], volume=4.0)
    comps = components(nodes(), 0.0, [a, b])
    assert len(comps.volumes) == 1
    assert comps.volumes[0] == pytest.approx(7.0)


def test_one_junction_split_into_pieces_contributes_one_vertex_each():
    """Two pieces of the same trimmed junction are not joined by sharing a node.

    Each reaches its own neighbour along a different strut axis, so they stay in
    separate components even though both sit at (5, 5, 5).
    """
    a = FakePiece((5, 5, 5), caps=[0], volume=1.0)
    b = FakePiece((5, 5, 5), caps=[1], volume=2.0)
    comps = components(nodes(across((5, 5, 5), 0), across((5, 5, 5), 1)), 4.0, [a, b])
    assert len(comps.vertices) == 4
    assert len(comps.volumes) == 2, "pieces of one junction are not joined by sharing a node"
    assert sorted(comps.volumes.values()) == [pytest.approx(5.0), pytest.approx(6.0)]


def test_a_cap_given_up_with_no_partner_is_a_hard_failure():
    """The regression this whole interface rework exists for.

    A piece that gives up a cap the other side never gave up has punched a hole
    with nothing behind it. Silently skipping it — what the previous
    implementation did — deferred the consequence to the stitcher, which
    discovered it hours later as "N of M stitched shells are not closed" with no
    indication of where. It has to fail here, naming the junction.
    """
    orphan = FakePiece((1, 0, 0), caps=[OPPOSITE_HALF[0]], volume=2.0)
    with pytest.raises(ProcessingError) as excinfo:
        build_components(nodes(), 0.0, [orphan], interfaces={((1, 0, 0), OPPOSITE_HALF[0])})
    assert "(1, 0, 0)" in str(excinfo.value)
    assert "(0, 0, 0)" in str(excinfo.value)


def test_an_interior_node_only_registers_caps_that_are_interfaces():
    """An interior cap facing nothing is exterior surface, not a connection.

    The interior shell keeps that cap face, so the graph must not claim a
    connection across it either — the two have to agree or the assembly does not
    close.
    """
    lone = nodes((0, 0, 0))
    comps = build_components(lone, 10.0, [], interfaces=set())
    assert len(comps.volumes) == 1
    assert comps.volumes[0] == pytest.approx(10.0)


def test_small_floating_component_is_dropped_and_a_large_one_is_not():
    big = FakePiece((0, 0, 0), caps=[], volume=100.0)
    small = FakePiece((9, 9, 9), caps=[], volume=0.5)
    comps = components(nodes(), 0.0, [big, small])
    dropped = comps.dropped(threshold=1.0)
    assert len(dropped) == 1
    assert comps.volumes[next(iter(dropped))] == pytest.approx(0.5)


def test_a_small_but_connected_piece_is_never_dropped():
    """specification.md §5: only *floating* sub-threshold bodies may be removed."""
    tiny = FakePiece((1, 0, 0), caps=[OPPOSITE_HALF[0]], volume=0.01)
    comps = components(nodes((0, 0, 0)), 100.0, [tiny])
    assert comps.dropped(threshold=1.0) == set()
    assert len(comps.volumes) == 1


def test_component_membership_covers_every_vertex():
    pieces = [FakePiece((0, 0, 0), caps=[0], volume=1.0),
              FakePiece((1, 0, 0), caps=[OPPOSITE_HALF[0]], volume=1.0)]
    comps = components(nodes((5, 5, 5)), 2.0, pieces)
    assert sum(len(m) for m in comps.members.values()) == len(comps.vertices)
    assert set(comps.labels) == set(comps.volumes)
