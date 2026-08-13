# Design measurements

The load-bearing assumptions behind the lattice-generation architecture
([../../docs/algorithm.md](../../docs/algorithm.md)), measured rather than
assumed. Each script in this directory prints PASS/FAIL against its own bar and
can be re-run:

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

**Cap integrity imposes no restriction of its own.** It was expected to need
roughly `t < cc/2`. A sweep of `t/cc` from 0.45 up to the `t < a` limit, at
`cc ∈ {5, 10, 20}` mm, keeps all six caps throughout — the caps only stop
existing where the CLI already rejects the parameters. The reason is exact: a
half-strut's reach toward an orthogonal strut direction is the profile's
**inradius** `t/2`, not its circumradius `r = t/√2`, because a diamond profile
presents an edge rather than a corner to every orthogonal direction. So caps are
intact for exactly `t < a`. See docs/algorithm.md §3.3; asserted in
`test/test_junction.py`.

---

## G2 — how the interior gets joined ✅ PASS

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

**The indexed build is what the interior path uses** (docs/algorithm.md §6):
global vertices keyed by `(owning node, local template vertex)` through a
precomputed cap correspondence, so neighbouring instances reference the same
`TopoDS_Vertex`/`TopoDS_Edge` objects and the shell is watertight by
construction. It is linear — **64,000 junctions / 1.55 M faces in 199 s** — and
its volume matches `N × volume(J)` exactly, which is a sound identity here
because adjacent junctions have zero-volume contact.

Two things this cost to get right, both recorded in the code:

* `BRepBuilderAPI_MakeFace(wire, onlyPlane=True)` infers an *arbitrary* plane
  normal, producing a shell that is closed but encloses zero volume. Faces must
  be built on a plane whose normal is stated from the template's outward normal.
* `BRep_Builder` leaves a hand-built shell's `Closed` flag false regardless of
  the geometry, and `BRepBuilderAPI_MakeSolid` reads the flag. Closure is
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

The localisation optimisation held in reserve (intersect the input body once per
spatial bucket, then trim junctions against the local piece) is **not needed** at
this input complexity and was not built.

---

## G4 — STEP writer throughput and round trip ✅ PASS

`STEPControl_Writer` on a 12×12×12 instanced interior shell, then re-read.

| | Bar | Measured |
|---|---|---|
| write throughput | ≥ 3,000 faces/s | **3,569 faces/s** (42,336 faces, 131.9 MB, 11.9 s) |
| round trip | solid count and volume preserved | 1 solid, volume relative error **5.8 × 10⁻¹³** |

Re-reading took 28.1 s — reading is slower than writing, and at large scale the
file itself becomes the dominant cost. That remains an open item for the scale
rehearsal (specification.md §10).

---

## Same-domain unification

Instancing merges nothing, so across every shared mid-strut interface two
coplanar lateral faces sit unmerged — every strut carrying eight lateral faces
where four suffice. `ShapeUpgrade_UnifySameDomain` before export recovers it, on
`dense-lattice`:

| | Before | After |
|---|---|---|
| Faces | 29,974 | **15,966** |
| Edges | 122,556 | **67,898** |
| Triangles (re-tessellated) | 85,832 | **62,152** |
| File size | 98.9 MB | **52.6 MB** |
| Volume (mm³) | 43574.0966 | 43574.0965 |

It costs ~8 s and pays for itself: export drops 9.3 s → 5.7 s and the round-trip
check 23.5 s → 13.5 s. The template itself is already minimal — unifying it
leaves 30 faces — so the redundancy is created entirely by instancing. See
docs/algorithm.md §9.

---

## G5 — what makes boundary stitching expensive ⚠️ ARCHITECTURAL FINDING

Run on the 21,955 trimmed boundary pieces (302,577 faces) the failed `cc=5, t=1`
scale rehearsal left in its temp folder, so these are the real shapes that took
4 h 45 m to sew, not a synthetic stand-in.

```
python tools/prototypes/g5_stitch_scaling.py --pieces path/to/temp/<stamp>
```

### G5a — scaling, and whether any configuration rescues it

| Pieces | `default` | `minimal` | Free edges | Merged edges |
|---|---|---|---|---|
| 500 | 1.46 s | 1.44 s | 2,000 | 2,575 |
| 1,000 | 5.99 s | 5.47 s | 4,378 | 5,263 |
| 2,000 | 24.34 s | 23.81 s | 8,961 | 10,812 |
| 4,000 | 77.16 s | 76.34 s | 16,235 | 20,587 |
| 8,000 | 195.77 s | 194.29 s | 29,550 | 36,273 |
| fitted exponent | **1.78** | **1.80** | | |

`minimal` turns off every optional phase — `Cutting` (splitting free edges so they
match), `Analysis` (degenerate-shape detection) and `SameParameter` (re-fitting
merged edges). **All three together account for under 2 %.** There is no
configuration of `BRepBuilderAPI_Sewing` that makes this affordable; the cost is
the merging itself.

Free edges grow linearly with piece count while time grows at `n^1.8`, so the
cost is roughly *quadratic in the number of interfaces to pair up* — the
signature of a search, which is exactly what it is. The pairing is already known
exactly from the junction graph (docs/algorithm.md §7.1), and sewing rediscovers
it.

### G5b — the face count alone, with nothing to merge

A **closed** instanced interior shell contributes faces and **zero** free edges,
so anything it costs is the price of indexing faces, not of the search:

| | Seconds |
|---|---|
| 4,000 pieces | **76.51** |
| 4,000 pieces + one closed 194,400-face shell | **716.61** |
| delta, for zero extra free edges | **+640.10** |

**This refutes the assumption the pipeline was built on.** `latticegen2.occ.sew`
documented that "sewing only works on free edges, so a shell whose interior edges
are already shared contributes only its holes to the workload", and that is what
justified handing `BRepBuilderAPI_Sewing` the whole interior shell. It is false:
9.4× the cost for nothing to merge. The interior shell at rehearsal scale is
705,000 faces — 3.6× the shell measured here — which is where the missing hours
of the 5 h 04 m run went.

### What this decides

* **Tuning the sewing call is dead.** Under 2 % is available (G5a).
* **Partitioning the boundary sew alone is not enough.** Even a perfectly tiled
  boundary stitch still has to join the interior shell, and that is the dominant
  term (G5b).
* **The interior shell must never enter sewing.** Which means stitching has to
  join known partners by index rather than by geometric search — the direction
  [specification.md](../../docs/specification.md) §10 sketches, now measured
  rather than assumed. Every interface is known from both sides after
  docs/algorithm.md §7.1, and an interior↔boundary cap is provably the whole
  template quad (§5.3(b)), so both sides' topology is addressable exactly.

Note that the shell counts in G5a are an artefact of sewing a *prefix* of the
boundary layer — an arbitrary subset has many genuinely unmatched free edges.
Only the timings and edge counts are meaningful as scale indicators.

---

## G6 — does tiling the boundary sew actually save wall time?

docs/specification.md §10 "Tile the boundary sew" proposes splitting a
component's pieces into spatial tiles, sewing each on its own, then sewing the
tile results together. G5b is a reason to be suspicious of that plan rather than
to trust it on faith: it already found that sewing pays a face-count cost even
where there is *nothing* to merge, and round 2 of tiling feeds sewing several
tile shells whose *combined* face count equals the untiled input's — only the
tile-boundary edges are genuinely free. If G5b's finding generalises, round 2
could cost close to what one monolithic sew costs, and tiling would only add
round 1's cost on top for nothing.

Measured directly, on the same real trimmed boundary pieces G5 used:

```
python tools/prototypes/g6_tile_stitch.py --pieces path/to/temp/<stamp>
```

Tiling here is by contiguous chunks of the loaded piece order, not by lattice
index — the saved `.brep` files carry geometry only, not the node each piece
belongs to. `boundary._split_batches`'s own docstring records that a worker's
batch already has spatial locality ("a worker's junctions share input-body
regions"), and pieces load in that same batch order, so contiguous chunks stand
in reasonably for real spatial tiles for the purpose of this measurement, which
is about the sewing call's cost model rather than the exact partition.

### 4,000 pieces

| | round 1 | round 2 | total | vs. baseline |
|---|---|---|---|---|
| baseline (1 call) | — | — | **91.3 s** | — |
| 4 tiles (1,000 each) | 21.5 s | 48.0 s | 69.4 s | **1.31×** |
| 8 tiles (500 each) | 11.0 s | 51.8 s | 62.8 s | **1.45×** |
| 16 tiles (250 each) | 6.3 s | 59.2 s | 65.5 s | **1.39×** |

### 8,000 pieces

| | round 1 | round 2 | total | vs. baseline |
|---|---|---|---|---|
| baseline (1 call) | — | — | **228.3 s** | — |
| 8 tiles (1,000 each) | 27.8 s | 131.8 s | 159.6 s | **1.43×** |
| 16 tiles (500 each) | 17.2 s | 148.3 s | 165.5 s | **1.38×** |
| 32 tiles (250 each) | 10.9 s | 159.3 s | 170.3 s | **1.34×** |

Same shape at twice the scale: round 1 keeps shrinking with tile count, round 2
keeps *growing* slightly instead, and the best total is again around a few
hundred to ~1,000 pieces per tile rather than at the finest split tried.

Round 1 shrinks with tile count roughly as G5a's `n^1.8` model predicts. Round 2
does **not** shrink to match — it *grows* slightly as tiles get smaller — which
confirms the G5b-based suspicion above rather than refuting it: sewing several
already-mostly-closed shells still pays close to a monolithic sew's face-count
cost, so round 2 puts a floor under how much tiling alone can save. The result is
a real but bounded win with a shallow optimum around a few hundred to ~1,000
pieces per tile, not a runaway improvement from finer tiling — which is why
`latticegen2.weld.TILE_TARGET_PIECES` is set at 500, inside that plateau at both
scales measured, rather than pushed smaller.

This measurement sums round 1 serially; the real pipeline runs it across worker
processes, so production's wall-clock saving should exceed the 1.3–1.45× measured
here by roughly the parallel speedup on round 1's share of the total —
confirming by how much is what the `TD_HX_Indre_Volum` rehearsal re-run is for
(docs/specification.md §10), since this gate only reached 8,000 pieces per
component and the rehearsal runs at 21,955.

### What this decides

* **Tiling the boundary sew is worth doing, but it is not the whole answer.**
  1.3–1.45× measured at 4,000 and 8,000 pieces, likely more in production once
  round 1's parallelism is credited — real money on a multi-hour stage, but
  nowhere near cancelling the `n^1.8` term outright, because round 2 does not
  scale down with tile count.
* **The tile target is a measured choice, not a guess.** 500 pieces per tile
  matches G5a's cheapest-measured single-tile size and sits inside G6's observed
  optimum plateau at both scales measured; smaller tiles trade a shrinking
  round 1 for a growing round 2 for no net gain.
* **The next lever, if one is needed, is round 2 itself** — recursive tiling, or
  finding a way to make the second sew pay only for the tile-boundary free edges
  rather than for every face it is handed — but that is future work, not
  something this gate's data justifies building yet.
