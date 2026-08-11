"""The junction template and the interior shell built by instancing it.

These cover the two claims the whole architecture rests on: that the six
mid-strut cap quads survive the template fuse intact, and that instancing the
template produces a watertight solid whose volume is exactly
``n_nodes * volume(J)``.

That volume identity is an exact, independent check rather than an approximate
one. Half-struts partition the struts, so all strut-strut overlap is internal to
a single junction and adjacent junctions meet with zero-volume contact — the
union's volume is therefore the plain sum, with no inclusion-exclusion term.
"""

import numpy as np
import pytest

from latticegen2 import occ
from latticegen2.interior import build_interior_shell, extract_template_mesh
from latticegen2.junction import build_template, is_cap_plane_face
from latticegen2.lattice import HALF_STRUTS, half_strut_offset, lattice_params


@pytest.fixture(scope="module")
def template():
    lp = lattice_params(10.0, 1.5)
    tpl = build_template(lp)
    return lp, tpl, extract_template_mesh(lp, tpl)


def grid(m: int) -> np.ndarray:
    a = np.arange(m, dtype=np.int64)
    g = np.meshgrid(a, a, a, indexing="ij")
    return np.stack([x.ravel() for x in g], axis=1)


def test_template_has_six_intact_caps(template):
    lp, tpl, _ = template
    assert len(tpl.cap_faces) == 6
    for h, face in enumerate(tpl.cap_faces):
        assert occ.area(face) == pytest.approx(lp.t ** 2, rel=1e-9)


def test_template_cap_centres_are_at_half_a_cell(template):
    lp, tpl, _ = template
    for h in range(6):
        assert np.allclose(tpl.cap_centers[h], half_strut_offset(lp, h))
        assert np.linalg.norm(tpl.cap_centers[h]) == pytest.approx(lp.a / 2, abs=1e-12)


def test_template_is_a_valid_solid(template):
    _, tpl, _ = template
    assert occ.is_valid(tpl.solid)
    assert tpl.volume > 0


def test_template_volume_is_below_six_disjoint_half_struts(template):
    """The junction is a union, so overlap near the node must reduce the total."""
    lp, tpl, _ = template
    disjoint = 6 * (lp.t ** 2) * (lp.a / 2)
    assert tpl.volume < disjoint
    assert tpl.volume > 0.5 * disjoint


@pytest.mark.parametrize("ratio", [0.05, 0.25, 0.5, 0.6, 0.7])
def test_caps_stay_intact_across_the_whole_valid_parameter_range(ratio):
    """There is no narrower window than the CLI's own ``t < a`` constraint.

    The design was expected to need roughly ``t < cc/2``; it does not, because a
    diamond profile presents its edge rather than its corner to every orthogonal
    strut direction (see :func:`test_profile_support_is_the_inradius`).
    """
    cc = 10.0
    lp = lattice_params(cc, cc * ratio)
    assert lp.t < lp.a, "the case must be inside the CLI's valid range"
    tpl = build_template(lp)
    assert len(tpl.cap_faces) == 6
    for face in tpl.cap_faces:
        assert occ.area(face) == pytest.approx(lp.t ** 2, rel=1e-9)


def test_profile_support_is_the_inradius(template):
    """Why caps survive: the support of profile j along any e_k is t/2, not r.

    ``t/2 < a/2`` is exactly ``t < a``, so cap integrity adds no restriction
    beyond the cross-constraint the CLI already enforces.
    """
    lp, _, _ = template
    for j in range(3):
        for k in range(3):
            if j == k:
                continue
            support = lp.r * max(abs(np.dot(lp.u[j], lp.e[k])), abs(np.dot(lp.v[j], lp.e[k])))
            assert support == pytest.approx(lp.t / 2, abs=1e-12)
            assert support < lp.a / 2


def test_template_mesh_pairs_every_incoming_cap_vertex(template):
    lp, tpl, tmesh = template
    incoming = [i for i, c in enumerate(tmesh.vertex_cap) if c >= 3]
    assert len(incoming) == 12, "three incoming caps of four corners each"
    for i in incoming:
        j = tmesh.cap_partner[i]
        k, _ = HALF_STRUTS[tmesh.vertex_cap[i]]
        shift = lp.a * lp.e[k]
        assert np.allclose(tmesh.verts[i] + shift, tmesh.verts[j], atol=1e-9)


def test_template_mesh_marks_exactly_six_cap_faces(template):
    _, _, tmesh = template
    assert sorted(c for c in tmesh.face_cap if c >= 0) == list(range(6))


def test_cap_plane_detection_finds_caps_and_rejects_lateral_faces(template):
    lp, tpl, _ = template
    origin = np.zeros(3)
    found = {is_cap_plane_face(lp, f, origin) for f in tpl.cap_faces}
    assert found == set(range(6))
    laterals = [f for f in occ.faces(tpl.solid) if not any(f.IsSame(c) for c in tpl.cap_faces)]
    assert laterals, "the template must have non-cap faces"
    assert all(is_cap_plane_face(lp, f, origin) is None for f in laterals)


@pytest.mark.parametrize("m", [1, 2, 3])
def test_instanced_grid_is_watertight_with_exact_volume(template, m):
    lp, tpl, tmesh = template
    ns = grid(m)
    kept = {(int(a), int(b), int(c)) for a, b, c in ns}
    shell, stats = build_interior_shell(lp, tpl, tmesh, ns, kept)
    assert stats["interior_open_edges"] == 0
    assert shell.Closed()
    solid = occ.make_solid(shell)
    assert occ.is_valid(solid)
    assert occ.volume(solid) == pytest.approx(len(ns) * tpl.volume, rel=1e-9)


def test_unification_is_a_no_op_on_the_template(template):
    """The template is already minimal — the redundancy comes from instancing.

    Pins the root-cause finding behind the export-time unification step: within
    one junction the +k and -k lateral faces are coplanar but *not* adjacent,
    because the other four half-struts cut them apart at the node. If a future
    change to the profile or the fuse made the template itself reducible, the
    reasoning in docs/algorithm.md §9 would need revisiting.
    """
    _, tpl, _ = template
    merged = occ.unify_same_domain(tpl.solid)
    assert len(occ.faces(merged)) == tpl.n_faces
    assert occ.volume(merged) == pytest.approx(tpl.volume, rel=1e-12)


def test_unification_halves_an_instanced_grid_without_moving_it(template):
    """Across every shared cap, two coplanar lateral faces become one."""
    lp, tpl, tmesh = template
    ns = grid(3)
    kept = {(int(a), int(b), int(c)) for a, b, c in ns}
    shell, _ = build_interior_shell(lp, tpl, tmesh, ns, kept)
    solid = occ.make_solid(shell)

    before_faces, before_edges = occ.count_subshapes(solid)
    merged = occ.unify_same_domain(solid)
    after_faces, after_edges = occ.count_subshapes(merged)

    # The exact count is OCCT's result, not ours; assert the effect, not a number.
    assert after_faces < 0.75 * before_faces, (before_faces, after_faces)
    assert after_edges < before_edges
    # ...and that it only re-described the boundary, never moved it.
    assert len(occ.solids(merged)) == 1
    assert occ.volume(merged) == pytest.approx(len(ns) * tpl.volume, rel=1e-9)
    assert occ.is_valid(merged)


def test_adjacent_instances_share_topology(template):
    """Neighbours must reference one shared vertex set, not two coincident ones."""
    lp, tpl, tmesh = template
    single, s1 = build_interior_shell(lp, tpl, tmesh, grid(1), {(0, 0, 0)})
    pair_nodes = np.array([[0, 0, 0], [1, 0, 0]], dtype=np.int64)
    _, s2 = build_interior_shell(lp, tpl, tmesh, pair_nodes, {(0, 0, 0), (1, 0, 0)})
    # Two independent junctions would have 2x the vertices; sharing one cap's
    # four corners must leave strictly fewer.
    assert s2["interior_vertices"] == 2 * s1["interior_vertices"] - 4


def test_unpaired_caps_are_kept_as_exterior_surface(template):
    """A junction with no neighbours must still close, using all six caps."""
    lp, tpl, tmesh = template
    shell, stats = build_interior_shell(lp, tpl, tmesh, grid(1), {(0, 0, 0)})
    assert shell.Closed()
    assert stats["interior_faces"] == tpl.n_faces
    assert occ.volume(occ.make_solid(shell)) == pytest.approx(tpl.volume, rel=1e-9)


def test_interior_shell_opens_exactly_where_a_neighbour_is_absent(template):
    """Dropping a neighbour from the kept set must leave one square hole."""
    lp, tpl, tmesh = template
    ns = np.array([[0, 0, 0]], dtype=np.int64)
    # Claim a neighbour exists along +e0 without instancing it: the cap toward it
    # is dropped, leaving that quad's four edges free.
    _, stats = build_interior_shell(lp, tpl, tmesh, ns, {(0, 0, 0), (1, 0, 0)})
    assert stats["interior_open_edges"] == 4
    assert stats["interior_faces"] == tpl.n_faces - 1
