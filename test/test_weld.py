"""Assembly by adoption (docs/algorithm.md §8).

The boundary layer is sewn to itself, and the interior is then built *onto* the
topology that came out — adopting its vertices and edges rather than making its
own and reconciling afterwards. These tests pin the two properties that took the
longest to get right, both of which fail silently when they are wrong:

* the assembled surface must be closed **and orientable**, because a shell whose
  faces are joined back-to-front still gives every edge exactly two users;
* an adopted edge's own direction is not ours to choose, so the index has to
  observe it rather than assume it.
"""

import numpy as np
import pytest

from OCP.BRep import BRep_Builder
from OCP.TopAbs import TopAbs_ShapeEnum
from OCP.TopoDS import TopoDS, TopoDS_Shell

from latticegen2 import occ, weld
from latticegen2.boundary import trim_junction
from latticegen2.connect import lattice_interfaces
from latticegen2.interior import build_interior_shell, extract_template_mesh
from latticegen2.junction import build_template, is_cap_plane_face
from latticegen2.lattice import OPPOSITE_HALF, lattice_params, neighbor_step, nodes

CC, T = 10.0, 1.5


@pytest.fixture(scope="module")
def template():
    lp = lattice_params(CC, T)
    tpl = build_template(lp)
    return lp, tpl, extract_template_mesh(lp, tpl)


def big_box():
    return occ.prism(
        occ.polygon_face(np.array([[-50.0, -50.0, -50.0], [50.0, -50.0, -50.0],
                                   [50.0, 50.0, -50.0], [-50.0, 50.0, -50.0]])),
        np.array([0.0, 0.0, 100.0]),
    )


# --- ring matching ----------------------------------------------------------


def test_a_trimmed_cap_matches_the_ring_predicted_from_lattice_maths(template):
    """The interior's side of an interface is known before any geometry exists.

    That is what lets an interface that will not join be rejected while both
    sides still have their cap face to fall back on.
    """
    lp, tpl, tmesh = template
    node = (0, 0, 0)
    faces, tags, _ = trim_junction(lp, tpl, np.zeros(3), big_box())[0]
    cap = next(f for f, t in zip(faces, tags) if t == 0)
    predicted = weld.template_cap_ring(lp, tmesh, node, 0)
    assert weld.match_rings(predicted, weld.ring_of_face(cap)) is not None


def test_rings_of_different_size_never_match(template):
    lp, _, tmesh = template
    four = weld.template_cap_ring(lp, tmesh, (0, 0, 0), 0)
    three = weld.Ring(edges=four.edges[:3], verts=four.verts[:3], points=four.points[:3])
    assert weld.match_rings(four, three) is None


def test_matching_is_indifferent_to_which_way_each_edge_runs(template):
    """The two sides are free to have parametrised the same segment either way."""
    lp, _, tmesh = template
    ring = weld.template_cap_ring(lp, tmesh, (0, 0, 0), 0)
    flipped = weld.Ring(
        edges=list(ring.edges),
        verts=list(ring.verts),
        points=[p[::-1] for p in ring.points],
    )
    assert weld.match_rings(ring, flipped) is not None


# --- the watertightness proof ----------------------------------------------


def shell_of(faces):
    shell = TopoDS_Shell()
    builder = BRep_Builder()
    builder.MakeShell(shell)
    for f in faces:
        builder.Add(shell, f)
    return shell


def test_a_closed_shell_has_no_defects(template):
    lp, tpl, _ = template
    faces, _, _ = trim_junction(lp, tpl, np.zeros(3), big_box())[0]
    assert weld.shell_defects(shell_of(faces)) [:2] == (0, 0)


def test_a_face_joined_back_to_front_is_caught(template):
    """Counting uses alone would pass this: every edge still has two faces.

    Measured on the 80 mm ball while this was being built — 0 open edges, 954
    edges traversed the same way twice, and a volume of 29,111 mm³ against a true
    51,393 mm³, with nothing else looking wrong.
    """
    lp, tpl, _ = template
    faces, _, _ = trim_junction(lp, tpl, np.zeros(3), big_box())[0]
    flipped = [faces[0].Reversed()] + list(faces[1:])
    open_edges, misoriented, _ = weld.shell_defects(shell_of(flipped))
    assert open_edges == 0, "the flip does not open the shell — that is the point"
    assert misoriented > 0


def test_a_missing_face_is_caught(template):
    lp, tpl, _ = template
    faces, _, _ = trim_junction(lp, tpl, np.zeros(3), big_box())[0]
    open_edges, _, samples = weld.shell_defects(shell_of(faces[1:]))
    assert open_edges > 0
    assert samples, "an open edge must be located, not just counted"


# --- adoption ---------------------------------------------------------------


def test_the_interior_joins_a_trimmed_junction_by_adopting_its_topology(template):
    """The whole mechanism, at two nodes.

    Node (0,0,0) is instanced by the index; node (1,0,0) comes out of a boolean,
    exactly as a real boundary piece does. The index adopts the boundary's ring,
    so the two share four edges rather than owning two coincident copies — and
    the assembled solid's volume is the exact sum, which is only true if nothing
    was double-counted or lost.
    """
    lp, tpl, tmesh = template
    a, b = (0, 0, 0), (1, 0, 0)
    h = 0
    interfaces = {(a, h), (b, OPPOSITE_HALF[h])}

    pos_b = nodes(lp, np.array([b], dtype=np.int64))[0]
    faces, tags, vol_b = trim_junction(lp, tpl, pos_b, big_box())[0]
    assert vol_b == pytest.approx(tpl.volume, rel=1e-12), "the box contains it whole"
    boundary_faces = [f for f, t in zip(faces, tags) if t != OPPOSITE_HALF[h]]

    rings = weld.interface_rings(lp, tmesh, {0: boundary_faces}, {(a, h): 0})
    assert set(rings) == {(a, h)}

    built = build_interior_shell(
        lp, tpl, tmesh, np.array([a], dtype=np.int64), interfaces,
        groups={a: 0}, adopted=rings,
    )
    shells = weld.assemble(built.shells, {0: boundary_faces})
    solids, _ = weld.close_shells(shells)

    assert len(solids) == 1
    assert occ.is_valid(solids[0])
    assert occ.volume(solids[0]) == pytest.approx(2 * tpl.volume, rel=1e-9)


def test_adoption_shares_the_edges_rather_than_duplicating_them(template):
    """Four edges are shared, so the assembly has four fewer than the two parts."""
    lp, tpl, tmesh = template
    a, b = (0, 0, 0), (1, 0, 0)
    interfaces = {(a, 0), (b, OPPOSITE_HALF[0])}
    pos_b = nodes(lp, np.array([b], dtype=np.int64))[0]
    faces, tags, _ = trim_junction(lp, tpl, pos_b, big_box())[0]
    boundary_faces = [f for f, t in zip(faces, tags) if t != OPPOSITE_HALF[0]]

    rings = weld.interface_rings(lp, tmesh, {0: boundary_faces}, {(a, 0): 0})
    built = build_interior_shell(
        lp, tpl, tmesh, np.array([a], dtype=np.int64), interfaces,
        groups={a: 0}, adopted=rings,
    )
    assembled = weld.assemble(built.shells, {0: boundary_faces})[0]

    separate = len(occ.faces(shell_of(boundary_faces)))
    assert len(occ.faces(assembled)) == separate + built.stats["interior_faces"]
    # Every edge of the assembly is used twice: nothing is coincident-but-distinct.
    assert weld.shell_defects(assembled)[:2] == (0, 0)


def test_a_ring_the_boundary_does_not_present_is_a_named_failure(template):
    """Silence here would become an unclosed shell much later, with no location."""
    lp, tpl, tmesh = template
    faces, tags, _ = trim_junction(lp, tpl, np.zeros(3), big_box())[0]
    with pytest.raises(Exception, match="no unique free edge"):
        weld.interface_rings(lp, tmesh, {0: faces}, {((7, 7, 7), 0): 0})
