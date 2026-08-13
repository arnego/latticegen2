"""Connectivity and the floating-body rule, answered by graph, not by geometry.

specification.md §5 forbids emitting a floating body smaller than ``t³``, and is
emphatic that "floating" means *provably disconnected* — a small fragment that is
still attached to the rest must never be deleted. Establishing that by experiment,
attempting a boolean and seeing what merges, is both slow and inconclusive when
the boolean itself fails on degenerate input.

Here connectivity is known by construction. Two junctions are joined exactly
when they share a mid-strut cap interface, and which caps are interfaces is
decided once, symmetrically, by :func:`latticegen2.boundary.resolve_interfaces`.
So the question reduces to connected components of a small graph: vertices are
instanced junctions (a trimmed boundary junction contributes one vertex per
connected piece it split into), edges are interfaces. No boolean is involved and
there is no unresolvable case.

This module deliberately stays free of the geometry kernel — it is pure graph
work over the interface set it is handed, which is what keeps
``test/test_connect.py`` runnable without OCCT.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .errors import ProcessingError
from .lattice import OPPOSITE_HALF, neighbor_step

NodeKey = tuple[int, int, int]


def lattice_interfaces(node_set) -> set[tuple[NodeKey, int]]:
    """Interfaces of a set of nodes that all present their six caps whole.

    The synthetic all-interior case: two nodes are joined exactly when they are
    neighbours and both are present. The real pipeline uses
    :func:`latticegen2.boundary.resolve_interfaces` instead, which also accounts
    for caps a boundary trim cut down or removed.
    """
    steps = [tuple(int(x) for x in neighbor_step(h)) for h in range(6)]
    out: set[tuple[NodeKey, int]] = set()
    for node in node_set:
        for h in range(6):
            step = steps[h]
            other = (node[0] + step[0], node[1] + step[1], node[2] + step[2])
            if other in node_set:
                out.add((node, h))
    return out


@dataclass
class Vertex:
    """One instanced junction (or one piece of a trimmed one)."""

    node: NodeKey
    volume: float
    is_interior: bool
    piece_index: int = -1
    """Index into the boundary piece list, or ``-1`` for an interior node."""


class _UnionFind:
    def __init__(self, n: int):
        self.parent = list(range(n))

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


@dataclass
class Components:
    vertices: list[Vertex]
    labels: np.ndarray
    """``labels[i]`` is the component id of ``vertices[i]``."""
    volumes: dict[int, float]
    members: dict[int, list[int]]

    def dropped(self, threshold: float) -> set[int]:
        """Component ids that are floating bodies below ``threshold`` volume."""
        return {cid for cid, vol in self.volumes.items() if vol < threshold}


def build_components(
    interior_nodes: np.ndarray,
    interior_volume: float,
    boundary_pieces,
    interfaces: set[tuple[NodeKey, int]],
) -> Components:
    """Connected components of the junction graph.

    ``interfaces`` holds both sides of every cap that is stitched across, so a
    node registers cap ``h`` exactly when ``(node, h)`` is in it. For a boundary
    piece that is the same thing as its own ``caps`` — the caps it gave up — which
    is how a junction split into several pieces attributes each interface to the
    piece that actually carries it.
    """
    vertices: list[Vertex] = []
    cap_owner: dict[tuple[NodeKey, int], list[int]] = {}

    for row in interior_nodes:
        node = (int(row[0]), int(row[1]), int(row[2]))
        vid = len(vertices)
        vertices.append(Vertex(node=node, volume=interior_volume, is_interior=True))
        for h in range(6):
            if (node, h) in interfaces:
                cap_owner.setdefault((node, h), []).append(vid)

    for pi, piece in enumerate(boundary_pieces):
        vid = len(vertices)
        vertices.append(
            Vertex(node=piece.node, volume=piece.volume, is_interior=False, piece_index=pi)
        )
        for h in piece.caps:
            cap_owner.setdefault((piece.node, int(h)), []).append(vid)

    uf = _UnionFind(len(vertices))
    steps = [neighbor_step(h) for h in range(6)]
    for (node, h), owners in cap_owner.items():
        step = steps[h]
        other = (node[0] + int(step[0]), node[1] + int(step[1]), node[2] + int(step[2]))
        partners = cap_owner.get((other, OPPOSITE_HALF[h]))
        if not partners:
            # Both sides of an interface give up their cap, so both sides are
            # always registered here. A missing partner means this junction has
            # punched a hole with nothing behind it, and the output would not be
            # watertight — the failure `resolve_interfaces` exists to prevent.
            # Reported here, in seconds, rather than discovered hours later as an
            # unclosed shell coming out of the stitcher.
            raise ProcessingError(
                f"Junction {node} gives up cap {h} as an interface, but the "
                f"junction across it at {other} does not give up cap "
                f"{OPPOSITE_HALF[h]}. That is a hole with no matching hole to "
                f"stitch to, so the output could not be watertight."
            )
        if h >= 3:
            continue  # each interface is unioned once, from its outgoing side
        # When a trim splits a cap region across several pieces, joining every
        # pair over-connects rather than under-connects. That is the safe
        # direction: it can only keep a body that might have been droppable,
        # never delete one that was attached (specification.md §5).
        for a in owners:
            for b in partners:
                uf.union(a, b)

    roots = {}
    labels = np.empty(len(vertices), dtype=np.int64)
    for i in range(len(vertices)):
        root = uf.find(i)
        if root not in roots:
            roots[root] = len(roots)
        labels[i] = roots[root]

    volumes: dict[int, float] = {}
    members: dict[int, list[int]] = {}
    for i, v in enumerate(vertices):
        cid = int(labels[i])
        volumes[cid] = volumes.get(cid, 0.0) + v.volume
        members.setdefault(cid, []).append(i)

    return Components(vertices=vertices, labels=labels, volumes=volumes, members=members)
