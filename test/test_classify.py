"""Classification: distance primitives, the mesh gate, and node classes.

The pure algorithm is tested against a synthetic analytic mesh so the maths is
covered without needing a CAD file, alongside the mesh-coverage gate cases from
docs/algorithm.md §5.1.
"""

import os

import numpy as np
import pytest

from latticegen2 import classify as classify_mod
from latticegen2 import occ
from latticegen2.classify import (
    Klass,
    Occupancy,
    PointInside,
    SpatialHash,
    TriMesh,
    _max_distance_to_mesh,
    _point_triangle_dist,
    _segment_segment_dist,
    check_surface_mesh_coverage,
    chordal_target,
    classify_nodes,
    measure_mesh_deviation,
    segment_triangle_dist,
    surface_samples,
    tessellate_surface,
)
from latticegen2.errors import InputGeometryError
from latticegen2.lattice import candidate_nodes, lattice_params
from latticegen2.parallel import WorkerPool

TESTDIR = os.path.dirname(os.path.abspath(__file__))
CYLINDER = os.path.join(TESTDIR, "test-cylinder.STEP")
BALL = os.path.join(TESTDIR, "80mm-test-ball.step")


def icosphere(radius: float, subdivisions: int = 2) -> TriMesh:
    """A closed analytic sphere mesh — no CAD kernel involved."""
    phi = (1 + 5 ** 0.5) / 2
    verts = np.array(
        [[-1, phi, 0], [1, phi, 0], [-1, -phi, 0], [1, -phi, 0],
         [0, -1, phi], [0, 1, phi], [0, -1, -phi], [0, 1, -phi],
         [phi, 0, -1], [phi, 0, 1], [-phi, 0, -1], [-phi, 0, 1]],
        dtype=float,
    )
    tris = np.array(
        [[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
         [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
         [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
         [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]],
        dtype=np.int64,
    )
    for _ in range(subdivisions):
        new_tris = []
        cache: dict[tuple[int, int], int] = {}
        verts = list(verts)

        def midpoint(a, b):
            key = (min(a, b), max(a, b))
            if key not in cache:
                cache[key] = len(verts)
                verts.append((verts[a] + verts[b]) / 2.0)
            return cache[key]

        for a, b, c in tris:
            ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
            new_tris += [[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]
        verts = np.array(verts)
        tris = np.array(new_tris, dtype=np.int64)
    verts = radius * verts / np.linalg.norm(verts, axis=1)[:, None]
    return TriMesh(verts=verts, tris=tris, deviation=0.0)


# --- distance primitives ----------------------------------------------------


def test_point_triangle_distance_interior_and_edge_cases():
    A = np.array([[0.0, 0, 0]])
    B = np.array([[1.0, 0, 0]])
    C = np.array([[0.0, 1, 0]])
    # Above the interior: perpendicular distance.
    assert _point_triangle_dist(np.array([0.25, 0.25, 2.0]), A, B, C)[0] == pytest.approx(2.0)
    # Beyond a vertex: distance to that vertex.
    assert _point_triangle_dist(np.array([-1.0, -1.0, 0.0]), A, B, C)[0] == pytest.approx(np.sqrt(2))
    # In-plane, outside one edge.
    assert _point_triangle_dist(np.array([2.0, 0.0, 0.0]), A, B, C)[0] == pytest.approx(1.0)


def test_segment_segment_distance():
    p, q = np.array([0.0, 0, 0]), np.array([1.0, 0, 0])
    # Parallel, offset in z.
    d = _segment_segment_dist(p, q, np.array([[0.0, 0, 3]]), np.array([[1.0, 0, 3]]))
    assert d[0] == pytest.approx(3.0)
    # Crossing in projection but separated in z.
    d = _segment_segment_dist(p, q, np.array([[0.5, -1, 2]]), np.array([[0.5, 1, 2]]))
    assert d[0] == pytest.approx(2.0)


def test_segment_triangle_distance_is_zero_when_piercing():
    A = np.array([[0.0, 0, 0]])
    B = np.array([[1.0, 0, 0]])
    C = np.array([[0.0, 1, 0]])
    d = segment_triangle_dist(np.array([0.2, 0.2, -1.0]), np.array([0.2, 0.2, 1.0]), A, B, C)
    assert d[0] == pytest.approx(0.0)


def test_segment_triangle_distance_when_clear():
    A = np.array([[0.0, 0, 0]])
    B = np.array([[1.0, 0, 0]])
    C = np.array([[0.0, 1, 0]])
    d = segment_triangle_dist(np.array([0.2, 0.2, 5.0]), np.array([0.2, 0.2, 7.0]), A, B, C)
    assert d[0] == pytest.approx(5.0)


# --- spatial structures -----------------------------------------------------


def test_spatial_hash_finds_the_triangles_near_a_box():
    mesh = icosphere(10.0, subdivisions=1)
    sh = SpatialHash(mesh)
    lo = np.array([-11.0, -11.0, -11.0])
    hi = np.array([11.0, 11.0, 11.0])
    assert len(sh.query(lo, hi)) == len(mesh.tris)
    far = sh.query(np.array([100.0, 100.0, 100.0]), np.array([101.0, 101.0, 101.0]))
    assert len(far) == 0


def test_spatial_hash_honours_a_minimum_cell_size():
    """A fine mesh must not force a fine grid — that is what bounds query cost."""
    mesh = icosphere(10.0, subdivisions=3)
    fine = SpatialHash(mesh)
    coarse = SpatialHash(mesh, min_cell=5.0)
    assert coarse.cell >= 5.0
    assert coarse.cell > fine.cell
    # Coarsening must not lose anything: a query still finds every triangle it
    # would have found, it just returns more candidates.
    lo = np.array([-11.0, -11.0, -11.0])
    hi = np.array([11.0, 11.0, 11.0])
    assert set(fine.query(lo, hi).tolist()) <= set(coarse.query(lo, hi).tolist())


def test_classification_is_unchanged_by_a_much_finer_mesh():
    """Mesh density must not move nodes between classes.

    The grid's cell size is derived from the mesh, so this also guards the
    minimum-cell floor above: if coarsening the grid dropped candidate
    triangles, near-surface nodes would silently stop being BOUNDARY.
    """
    lp = lattice_params(10.0, 1.5)
    cand = candidate_nodes(lp, np.array([-20.0] * 3), np.array([20.0] * 3))
    coarse = icosphere(18.0, subdivisions=2)
    fine = icosphere(18.0, subdivisions=4)
    coarse.deviation = fine.deviation = 1.0
    a = classify_nodes(lp, coarse, cand)
    b = classify_nodes(lp, fine, cand)
    differing = int(np.count_nonzero(a.node_class != b.node_class))
    assert differing <= 2, f"{differing} nodes changed class with mesh density"


def test_occupancy_is_conservative_but_excludes_the_far_field():
    mesh = icosphere(10.0, subdivisions=2)
    occupancy = Occupancy(mesh, 2.0)
    on_surface = np.array([[10.0, 0.0, 0.0]])
    deep_inside = np.array([[0.0, 0.0, 0.0]])
    far_away = np.array([[500.0, 0.0, 0.0]])
    assert occupancy.near(on_surface)[0]
    assert not occupancy.near(deep_inside)[0]
    assert not occupancy.near(far_away)[0]


def test_point_inside_on_an_analytic_sphere():
    mesh = icosphere(10.0, subdivisions=3)
    inside = PointInside(mesh)
    pts = np.array([[0.0, 0, 0], [5.0, 0, 0], [0.0, 0, 9.0], [20.0, 0, 0], [0.0, 15.0, 0]])
    got = inside(pts)
    assert list(got) == [True, True, True, False, False]


# --- node classification ----------------------------------------------------


def test_classification_of_a_sphere_is_layered_correctly():
    """Deep nodes INTERIOR, shell nodes BOUNDARY, far nodes OUTSIDE."""
    lp = lattice_params(10.0, 1.5)
    mesh = icosphere(30.0, subdivisions=3)
    mesh.deviation = 0.2
    cand = candidate_nodes(lp, np.array([-30.0] * 3), np.array([30.0] * 3))
    result = classify_nodes(lp, mesh, cand)
    counts = result.counts()
    assert counts["interior"] > 0 and counts["boundary"] > 0 and counts["outside"] > 0

    from latticegen2.lattice import nodes as node_positions

    pos = node_positions(lp, result.node_index)
    radius = np.linalg.norm(pos, axis=1)
    margin = lp.r + mesh.deviation
    interior = result.node_class == int(Klass.INTERIOR)
    outside = result.node_class == int(Klass.OUTSIDE)
    # Every INTERIOR node must be far enough inside that its junction is clear.
    assert np.all(radius[interior] < 30.0 - margin)
    # Every OUTSIDE node must be beyond the sphere.
    assert np.all(radius[outside] > 30.0 - lp.a / 2 - margin)


def test_worker_reuses_one_classify_index_across_the_slices_it_is_given(tmp_path):
    """Rebuilt *per worker* has to mean per worker, not per slice.

    ``_classify_parallel`` dispatches ``workers * 4`` slices and the pool hands
    them out one at a time, so a worker is called once per slice. With no cache
    it rebuilt the mesh and its three indices on every one of those calls —
    four times over on the default budget — while ``_worker_classify``'s own
    docstring and docs/algorithm.md §5.4 both priced the rebuild as paid once.

    Reuse is safe because :class:`~latticegen2.classify._ClassifyIndex` is a
    function of ``lp`` and the mesh alone and is read-only once built, which is
    the same property that makes the sweep divisible. This asserts object
    identity rather than a timing, since the saving is wall clock only — the
    classification is bit-identical either way, which the sibling test pins.
    """
    lp = lattice_params(10.0, 1.5)
    mesh = icosphere(30.0, subdivisions=3)
    mesh.deviation = 0.2
    cand = candidate_nodes(lp, np.array([-30.0] * 3), np.array([30.0] * 3))

    mesh_path = str(tmp_path / "mesh.npz")
    nodes_path = str(tmp_path / "nodes.npz")
    classify_mod.save_mesh(mesh, mesh_path)
    np.savez(nodes_path, candidates=np.asarray(cand, dtype=np.int64))

    classify_mod._WORKER_INDEX.clear()
    try:
        stride = 4
        parts = [
            classify_mod._worker_classify(
                (mesh_path, nodes_path, lp.cc, lp.t, offset, stride)
            )[0]
            for offset in range(stride)
        ]
        assert len(classify_mod._WORKER_INDEX) == 1, "one index, not one per slice"
        cached = classify_mod._WORKER_INDEX[(mesh_path, lp.cc, lp.t)]

        # And the cache is what the slices actually ran against, not a spare
        # copy built beside them: re-running a slice must hand back that object.
        classify_mod._worker_classify((mesh_path, nodes_path, lp.cc, lp.t, 0, stride))
        assert classify_mod._WORKER_INDEX[(mesh_path, lp.cc, lp.t)] is cached

        # Reuse must not have changed the answer.
        serial = classify_nodes(lp, mesh, cand).node_class
        reassembled = np.empty(len(cand), dtype=np.int64)
        for offset, part in enumerate(parts):
            reassembled[offset::stride] = part
        assert np.array_equal(serial, reassembled)
    finally:
        classify_mod._WORKER_INDEX.clear()


def test_parallel_classification_is_identical_to_the_serial_sweep(tmp_path):
    """The one classification path that crosses a process boundary.

    :func:`latticegen2.classify.classify_slice` divides the sweep on the claim
    that every node is decided independently of every other, so the parallel
    path must return not "equivalent" results but **the same** ones — the
    strided slicing and reassembly leave no room for a near-miss, and a
    tolerance here would hide exactly the off-by-one that stride arithmetic
    invites.

    This also exercises the IPC plumbing for real, the way
    ``test_weld.py``'s worker-path test does: a bad job tuple, a mesh that does
    not survive the ``.npz`` round trip, or a worker exception swallowed rather
    than propagated fails here rather than at rehearsal scale.
    """
    lp = lattice_params(10.0, 1.5)
    mesh = icosphere(30.0, subdivisions=3)
    mesh.deviation = 0.2
    cand = candidate_nodes(lp, np.array([-30.0] * 3), np.array([30.0] * 3))

    serial = classify_nodes(lp, mesh, cand)
    with WorkerPool(4) as pool:
        parallel = classify_nodes(lp, mesh, cand, tmpdir=str(tmp_path), pool=pool)

    counts = serial.counts()
    assert counts["interior"] > 0 and counts["boundary"] > 0 and counts["outside"] > 0, \
        "a sweep with an empty class would not test the reassembly at all"
    assert np.array_equal(serial.node_class, parallel.node_class)
    assert np.array_equal(serial.node_index, parallel.node_index)
    assert parallel.counts() == counts
    assert parallel.max_worker_rss > 0, "a worker actually ran and reported its RSS"


def test_interior_nodes_never_neighbour_an_outside_node():
    """The property :mod:`latticegen2.interior` relies on to drop every cap."""
    lp = lattice_params(10.0, 1.5)
    mesh = icosphere(30.0, subdivisions=3)
    mesh.deviation = 0.2
    cand = candidate_nodes(lp, np.array([-30.0] * 3), np.array([30.0] * 3))
    result = classify_nodes(lp, mesh, cand)
    klass = {tuple(int(x) for x in row): int(c)
             for row, c in zip(result.node_index, result.node_class)}
    from latticegen2.lattice import neighbor_step

    steps = [neighbor_step(h) for h in range(6)]
    for nodekey, c in klass.items():
        if c != int(Klass.INTERIOR):
            continue
        for step in steps:
            nb = (nodekey[0] + int(step[0]), nodekey[1] + int(step[1]), nodekey[2] + int(step[2]))
            assert klass.get(nb, int(Klass.BOUNDARY)) != int(Klass.OUTSIDE)


# --- the mesh gate ----------------------------------------------------------


@pytest.mark.parametrize("part,cc,t", [(CYLINDER, 10.0, 1.5), (BALL, 20.0, 4.0)])
def test_real_cad_tessellates_and_passes_the_coverage_gate(part, cc, t):
    """Real CAD must pass the coverage gate, not merely fail to crash it."""
    lp = lattice_params(cc, t)
    shape = occ.read_step(part)
    mesh = tessellate_surface(shape, lp)
    assert len(mesh.tris) > 0
    assert mesh.deviation >= chordal_target(lp)


def test_every_face_of_the_cylinder_is_meshed_to_its_exact_area():
    """Every face must mesh to its exact trimmed area, the tight quantity the
    coverage gate compares against (docs/algorithm.md §5.1)."""
    lp = lattice_params(10.0, 1.5)
    shape = occ.read_step(CYLINDER)
    occ.mesh_shape(shape, chordal_target(lp))
    for face in occ.faces(shape):
        exact = occ.area(face)
        verts, tris = occ.face_triangulation(face)
        a, b, c = verts[tris[:, 0]], verts[tris[:, 1]], verts[tris[:, 2]]
        meshed = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1).sum()
        assert meshed == pytest.approx(exact, rel=0.02)


def test_coverage_gate_rejects_a_genuinely_undercovered_face(monkeypatch):
    """A face reporting far more exact area than it meshes must be an error."""
    lp = lattice_params(10.0, 1.5)
    shape = occ.read_step(CYLINDER)
    occ.mesh_shape(shape, chordal_target(lp))
    real_area = occ.area

    def inflated(target):
        return real_area(target) * 100.0

    monkeypatch.setattr(occ, "area", inflated)
    with pytest.raises(InputGeometryError, match="covers only"):
        check_surface_mesh_coverage(shape, lp)


def test_coverage_gate_tolerates_a_sliver_shortfall():
    """A shortfall below one mesh element is noise, whatever its ratio.

    A pure relative test trips on the thousands of tiny trimmed faces a lattice
    has; the absolute floor is what keeps this function usable on generated
    output as well as on input CAD (docs/algorithm.md §5.1).
    """
    lp = lattice_params(10.0, 1.5)
    shape = occ.read_step(CYLINDER)
    occ.mesh_shape(shape, chordal_target(lp))
    # min(t, a)^2 = 2.25 mm^2 here, so a shortfall must exceed that to count.
    assert min(lp.t, lp.a) ** 2 == pytest.approx(2.25)
    check_surface_mesh_coverage(shape, lp)


def test_planar_faces_are_not_sampled_for_deviation():
    """A plane's triangulation is exact, so it contributes no samples at all."""
    lp = lattice_params(10.0, 1.5)
    shape = occ.read_step(CYLINDER)
    occ.mesh_shape(shape, chordal_target(lp))
    planar = [f for f in occ.faces(shape) if occ.is_planar(f)]
    assert planar
    assert all(surface_samples(f) is None for f in planar)


def test_sphere_deviation_is_the_real_sagitta_not_a_pole_artifact():
    """Regression for issue #6.

    The ball is two hemispherical faces, so it has a degenerate parametric edge
    at each pole where every cap triangle owns a vertex at ``v = ±π/2``.
    Comparing a sample against the triangle that owns its *parameters* — rather
    than against the nearest one — put the worst sample 90° of longitude away
    from the triangle it was measured against and reported 2.1988 mm on a mesh
    whose true worst sagitta is 0.0786 mm. At ``cc=10`` that cleared ``a/4`` and
    failed valid input (docs/algorithm.md §5.1).
    """
    lp = lattice_params(10.0, 1.0)
    shape = occ.read_step(BALL)
    mesh = tessellate_surface(shape, lp)

    assert mesh.deviation < lp.a / 4.0
    # Requested is min(t, a)/10 = 0.1 mm; the measured part must be of that
    # order, not the 2.2 mm the parametric pairing produced.
    assert mesh.deviation < 3.0 * chordal_target(lp)


def test_measured_deviation_is_never_below_the_true_sagitta():
    """The margin's guarantee needs ``d`` to be an upper bound on mesh error.

    Over-estimating only sends more junctions down the boundary path; under-
    estimating could let a strut that touches the surface pass as clear of it.
    The ball is an exact sphere of radius 40 centred on the origin, so the true
    sagitta at every triangle centroid is known analytically.

    The sagitta is taken from the mesh `tessellate_surface` actually returned,
    not from a second tessellation: meshing again with different parameters
    would compare a deviation measured on one mesh against a reference measured
    on another, and the assertion would only hold by coincidence.
    """
    lp = lattice_params(10.0, 1.0)
    shape = occ.read_step(BALL)
    mesh = tessellate_surface(shape, lp)

    centroids = mesh.verts[mesh.tris].mean(axis=1)
    true_sagitta = float((40.0 - np.linalg.norm(centroids, axis=1)).max())
    assert true_sagitta > 0.0
    assert mesh.deviation >= true_sagitta


def test_deviation_search_matches_an_exhaustive_nearest_triangle_scan():
    """The local search must agree with a brute-force scan over every triangle.

    Restricting the search to a grid neighbourhood is only allowed to
    over-estimate, and this pins that it does not need to.
    """
    lp = lattice_params(10.0, 1.5)
    shape = occ.read_step(CYLINDER)
    mesh = tessellate_surface(shape, lp)
    samples = np.vstack(
        [s for s in (surface_samples(f) for f in occ.faces(shape)) if s is not None]
    )
    subset = samples[::7]
    a, b, c = mesh.triangle_points
    brute = max(
        float(_point_triangle_dist(np.repeat(p[None, :], len(a), axis=0), a, b, c).min())
        for p in subset
    )
    assert _max_distance_to_mesh(subset, mesh) == pytest.approx(brute)
    # And the full measurement is at least as large, since it sees every sample.
    assert measure_mesh_deviation(shape, mesh) >= brute - 1e-12
