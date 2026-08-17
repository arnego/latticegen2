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
from latticegen2.connect import lattice_interfaces
from latticegen2.interior import _newell, build_interior_shell, extract_template_mesh

def _one_shell(built):
    """The single shell a synthetic all-interior grid produces."""
    shells = list(built.shells.values())
    assert len(shells) == 1, shells
    shells[0].Closed(built.stats["interior_open_edges"] == 0)
    return shells[0]

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
    built = build_interior_shell(lp, tpl, tmesh, ns, lattice_interfaces(kept))
    shell, stats = _one_shell(built), built.stats
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


def test_the_instanced_grid_arrives_already_merged(template):
    """Across every shared cap, two coplanar lateral faces are built as one.

    This merge used to be left to `simplify` (docs/algorithm.md §9), which had to
    rediscover by search over the whole solid a pairing that is known here by
    construction — one per surviving mid-strut interface. Building it merged
    instead is what docs/specification.md §10 Phase 2 does, so what this pins is
    the *effect*: the shell comes out at well under the un-merged face count, and
    the kernel then finds nothing left to do.
    """
    lp, tpl, tmesh = template
    ns = grid(3)
    kept = {(int(a), int(b), int(c)) for a, b, c in ns}
    ifaces = lattice_interfaces(kept)
    built = build_interior_shell(lp, tpl, tmesh, ns, ifaces)
    shell = _one_shell(built)
    solid = occ.make_solid(shell)

    # What the build would have emitted without merging: every template face at
    # every node, less the one cap each interface drops.
    unmerged = len(ns) * tpl.n_faces - len(ifaces)
    assert built.stats["interior_faces"] < 0.75 * unmerged

    before_faces, _ = occ.count_subshapes(solid)
    merged = occ.unify_same_domain(solid)
    after_faces, _ = occ.count_subshapes(merged)

    assert after_faces == before_faces, "the build left the kernel nothing to merge"
    # ...and the geometry is the same one the merge used to produce.
    assert len(occ.solids(merged)) == 1
    assert occ.volume(merged) == pytest.approx(len(ns) * tpl.volume, rel=1e-9)
    assert occ.is_valid(merged)


def test_adjacent_instances_share_topology(template):
    """Neighbours must reference one shared vertex set, not two coincident ones.

    A shared cap's four corners are not merely shared but *gone*: in the merged
    full-strut lateral face the two edges meeting at such a corner both run along
    the strut axis, so it stops being a corner and is dropped. Two junctions
    therefore lose eight vertices between them — four to sharing, four to the
    merge — which is the edge and vertex reduction Phase 2 delivers by
    construction rather than paying the kernel's edge pass for.
    """
    lp, tpl, tmesh = template
    s1 = build_interior_shell(lp, tpl, tmesh, grid(1), lattice_interfaces({(0, 0, 0)})).stats
    pair_nodes = np.array([[0, 0, 0], [1, 0, 0]], dtype=np.int64)
    s2 = build_interior_shell(lp, tpl, tmesh, pair_nodes, lattice_interfaces({(0, 0, 0), (1, 0, 0)})).stats
    assert s2["interior_vertices"] == 2 * s1["interior_vertices"] - 8


@pytest.mark.parametrize(
    "cc,t",
    [(cc, t) for cc in (5.0, 10.0, 20.0, 50.0) for t in (0.4, 1.0, 1.5, 4.0, 20.0)
     if t < cc / np.sqrt(2.0)],
)
def test_merged_lateral_faces_are_sound_across_the_parameter_range(cc, t):
    """docs/specification.md §10 Phase 2, risk R4, over G1's own sweep.

    The spliced wire has to be planar, correctly wound and exactly area-
    preserving for every ``(cc, t)`` the CLI accepts — a merged loop that is any
    of those things wrongly would produce a face that is not the union of the two
    it replaces. `_splice_lateral` enforces all three and returns ``None`` on
    failure, so the first assertion is that nothing had to fall back; the rest
    re-derive the area identity here rather than trusting the code that built it.
    """
    lp = lattice_params(cc, t)
    tpl = build_template(lp)
    tmesh = extract_template_mesh(lp, tpl)

    assert tmesh.merged_halves == {0, 1, 2}, "a strut family fell back to half-faces"
    assert len(tmesh.merged) == 12  # four lateral faces per outgoing half-strut

    for fi, mf in tmesh.merged.items():
        h = tmesh.face_half[fi]
        k, _ = HALF_STRUTS[h]
        shift = lp.a * lp.e[k]
        pts = np.array([tmesh.verts[v] + (shift if side else 0.0) for side, v in mf.loop])

        # Planar, and on the plane the template says this face lies in.
        assert np.abs((pts - pts[0]) @ mf.normal).max() < 1e-9
        # Wound so the right-hand normal is the template's outward one.
        newell = _newell(pts)
        assert np.dot(newell / np.linalg.norm(newell), mf.normal) == pytest.approx(1.0, abs=1e-9)
        # Exactly the two halves it replaces, no more and no less.
        half = 0.5 * np.linalg.norm(_newell(tmesh.verts[tmesh.loops[fi]]))
        area = 0.5 * np.linalg.norm(newell)
        assert area == pytest.approx(2.0 * half, rel=1e-12)


def test_a_merged_grid_is_watertight_and_exact_across_the_parameter_range():
    """The merged build must still close, and still hit ``N x volume(J)``.

    `build_interior_shell` raises on any edge used more than twice or twice the
    same way, so reaching a closed shell at all is the malformed-edge check R4
    asks for; the volume identity is what proves the merge did not move anything.
    """
    for cc, t in ((5.0, 1.0), (10.0, 1.5), (20.0, 4.0), (50.0, 20.0)):
        lp = lattice_params(cc, t)
        tpl = build_template(lp)
        tmesh = extract_template_mesh(lp, tpl)
        ns = grid(2)
        kept = {(int(a), int(b), int(c)) for a, b, c in ns}
        built = build_interior_shell(lp, tpl, tmesh, ns, lattice_interfaces(kept))
        assert built.stats["interior_open_edges"] == 0, (cc, t)
        solid = occ.make_solid(_one_shell(built))
        assert occ.is_valid(solid), (cc, t)
        assert occ.volume(solid) == pytest.approx(len(ns) * tpl.volume, rel=1e-9), (cc, t)


def test_unpaired_caps_are_kept_as_exterior_surface(template):
    """A junction with no neighbours must still close, using all six caps."""
    lp, tpl, tmesh = template
    built = build_interior_shell(lp, tpl, tmesh, grid(1), lattice_interfaces({(0, 0, 0)}))
    shell, stats = _one_shell(built), built.stats
    assert shell.Closed()
    assert stats["interior_faces"] == tpl.n_faces
    assert occ.volume(occ.make_solid(shell)) == pytest.approx(tpl.volume, rel=1e-9)


def test_interior_shell_opens_exactly_where_a_neighbour_is_absent(template):
    """Dropping a neighbour from the kept set must leave one square hole."""
    lp, tpl, tmesh = template
    ns = np.array([[0, 0, 0]], dtype=np.int64)
    # Claim a neighbour exists along +e0 without instancing it: the cap toward it
    # is dropped, leaving that quad's four edges free.
    stats = build_interior_shell(lp, tpl, tmesh, ns, lattice_interfaces({(0, 0, 0), (1, 0, 0)})).stats
    assert stats["interior_open_edges"] == 4
    assert stats["interior_faces"] == tpl.n_faces - 1
