# Phase 0 — de-risking measurements

The five load-bearing assumptions behind the fuse-free architecture
([../../docs/research/perf-rearchitecture-proposal.md](../../docs/research/perf-rearchitecture-proposal.md)),
measured before the implementation was built out. Each script in this directory
prints PASS/FAIL against its own bar and can be re-run:

```
python tools/prototypes/g1_cap_integrity.py
python tools/prototypes/g2_instancing_join.py
python tools/prototypes/g3_g4_boundary_export.py
```

Machine: Windows 11, 6-core CPU, 32 GB RAM, Python 3.11.13, OCP 7.9.3.1.1
(OCCT 7.9.3), NumPy 2.4.6.

---

## G1 — junction template cap integrity ✅ PASS

Builds the junction template `J` for every `(cc, t)` across the documented CLI
range and checks that all six mid-strut cap quads survive the fuse as intact
`t × t` faces, cross-checking `volume(J)` against an inclusion-exclusion
prediction computed from independent pairwise/triple intersections.

| Result | |
|---|---|
| Combinations tested | 16 valid `(cc, t)` pairs across `cc ∈ {5,10,20,50}`, `t ∈ {0.4,1,1.5,4,20}` |
| Caps intact | **6 of 6, every combination** |
| Template size | 30 faces, 32 vertices |
| `volume(J)` vs inclusion-exclusion | agreement to **≤ 4.4 × 10⁻¹⁶** relative, every combination |
| Build time | 38–51 ms |

**The expected parameter restriction does not exist.** The proposal anticipated
needing roughly `t < cc/2`. A sweep of `t/cc` from 0.45 up to the `t < a` limit,
at `cc ∈ {5, 10, 20}` mm, keeps all six caps throughout — the caps only stop
existing where the CLI already rejects the parameters. The reason is exact: a
half-strut's reach toward an orthogonal strut direction is the profile's
**inradius** `t/2`, not its circumradius `r = t/√2`, because a diamond profile
presents an edge rather than a corner to every orthogonal direction. So caps are
intact for exactly `t < a`. See docs/algorithm.md §3.3; asserted in
`test/test_junction.py`.

**Consequence:** the "fall back to the legacy pipeline outside the fast path's
parameter window" branch the guide called for was not needed and does not exist.

---

## G2 — instancing and join throughput ✅ PASS (with a redesign)

Three candidate mechanisms for joining instanced junctions, over an `m × m × m`
grid.

| Nodes | Mechanism | Build | Join | Solids | Faces | Volume error | Valid |
|---|---|---|---|---|---|---|---|
| 1,000 | `BOPAlgo_GlueFull` | 0.01 s | 9.54 s | **1,000** | 30,000 | — | — |
| 1,000 | `BRepBuilderAPI_Sewing` | 0.01 s | 14.46 s | 1 | 24,600 | 2.8 × 10⁻¹⁴ | yes |
| 8,000 | `BRepBuilderAPI_Sewing` | — | **> 250 s CPU, killed** | — | — | — | — |
| 64 | indexed shared topology | 0.15 s | — | 1 | 1,632 | 2.8 × 10⁻¹⁵ | yes |
| 512 | indexed shared topology | 1.20 s | — | 1 | 12,672 | 2.4 × 10⁻¹⁴ | yes |
| 1,000 | indexed shared topology | 3.03 s | — | 1 | 24,600 | — | yes |
| 8,000 | indexed shared topology | 23.79 s | — | 1 | 194,400 | — | yes |
| 64,000 | indexed shared topology | **198.52 s** | — | 1 | 1,545,600 | — | yes |

**Glue mode does not merge.** `BOPAlgo_Builder` with
`SetGlue(BOPAlgo_GlueFull)` is OCCT's dedicated mode for operands meeting only on
coincident faces, which is exactly this contact pattern. It returned 1,000 solids
from 1,000 inputs — no merging at all. Rejected.

**Sewing works but does not scale.** It has to *discover* by geometric search a
pairing that is already known exactly, and the cost grows clearly worse than
linearly: 14.9 s at 1,000 junctions, still running after 250 s of CPU at 8,000.

**The indexed build was written in response** and became the interior path
(docs/algorithm.md §6): global vertices keyed by `(owning node, local template
vertex)` through a precomputed cap correspondence, so neighbouring instances
reference the same `TopoDS_Vertex`/`TopoDS_Edge` objects and the shell is
watertight by construction. It is linear — **64,000 junctions / 1.55 M faces in
199 s**, inside the gate's 5-minute bar — and its volume matches
`N × volume(J)` exactly, which is a sound identity here because adjacent
junctions have zero-volume contact.

Two things this cost to get right, both recorded in the code:

* `BRepBuilderAPI_MakeFace(wire, onlyPlane=True)` infers an *arbitrary* plane
  normal, producing a shell that is closed but encloses zero volume. Faces must
  be built on a plane whose normal is stated from the template's outward normal.
* `BRep_Builder` leaves a hand-built shell's `Closed` flag false regardless of
  the geometry, and `BRepBuilderAPI_MakeSolid` reads the flag. Closure is now
  computed from the edge-use tally (every edge exactly twice, opposite
  directions) and the flag set from that.

Sewing is still used, but only where the pairing genuinely is not known — joining
trimmed boundary junctions, which scales with surface area rather than volume.

---

## G3 — per-junction intersection latency ✅ PASS

One single-operand `BRepAlgoAPI_Common` per boundary junction, against
`test/test-cylinder.STEP` at `cc=10, t=1.5`, over 200 of the run's 968 boundary
junctions.

| | Bar | Measured |
|---|---|---|
| median | < 50 ms | **24.0 ms** |
| p95 | < 250 ms | **31.7 ms** |
| max | — | 39.6 ms |

The tight spread is the point: these are constant-size, independent jobs, which
is what makes the boundary stage trivially parallel. Corroborated end-to-end —
the `dense-lattice` run trims all 968 boundary junctions in **6.1 s** across 5
worker processes.

The localisation optimisation the proposal held in reserve (intersect the input
body once per spatial bucket, then trim junctions against the local piece) is
**not needed** at this input complexity and was not built.

---

## G4 — STEP writer throughput and round trip ✅ PASS

`STEPControl_Writer` on a 12×12×12 instanced interior shell, then re-read.

| | Bar | Measured |
|---|---|---|
| write throughput | ≥ 3,000 faces/s | **3,569 faces/s** (42,336 faces, 131.9 MB, 11.9 s) |
| round trip | solid count and volume preserved | 1 solid, volume relative error **5.8 × 10⁻¹³** |

Re-reading took 28.1 s — reading is slower than writing, and at the projected
64× scale the file itself (multi-GB) becomes the dominant cost, exactly as the
proposal's §4 predicted. That remains an open item for the scale rehearsal
(specification.md §10).

---

## G5 — does COMMON preserve untouched faces bit-exactly? — **not needed**

The gate existed to decide whether trimmed boundary junctions could be attached
by exact index pairing, with tolerance-based sewing as the fallback. The
implementation takes the fallback **by design**, for a reason the gate did not
anticipate: a trimmed junction's faces come back from the boolean as new
topology that cannot be re-indexed against the template at all, whatever the
vertices happen to be. Boundary interfaces are therefore always stitched with a
tolerance (1e-6 mm, against a ~1e-14 mm real discrepancy and a ≥ 0.4 mm smallest
feature), and exact index pairing is used only on the interior path, where both
sides are template instances and no boolean has intervened.

This is not a weakening: the interior is where the scaling is, and it is exact.

---

## What Phase 0 changed about the plan

1. **No parameter-window restriction, and no legacy fallback path** (G1).
2. **The interior join is a custom indexed build, not sewing or glue** (G2) — the
   guide listed this as a "v2, only if needed" escalation; it was needed
   immediately.
3. **G5 became moot** once the boundary attachment strategy was settled.
