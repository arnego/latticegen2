# Design measurements

The load-bearing assumptions behind the lattice-generation architecture
([../../docs/algorithm.md](../../docs/algorithm.md)), measured rather than
assumed. Each script in this directory prints PASS/FAIL against its own bar and
can be re-run:

```
python tools/prototypes/g1_cap_integrity.py
python tools/prototypes/g2_instancing_join.py
python tools/prototypes/g3_g4_boundary_export.py
python tools/prototypes/g5_stitch_scaling.py
python tools/prototypes/g6_tile_stitch.py
python tools/prototypes/g7_thread_scaling.py
python tools/prototypes/g8_seam_only_sew.py
python tools/prototypes/g10_pinhole_wires.py
python tools/prototypes/g12_self_intersecting_wire.py
python tools/prototypes/g13_unify_scaling.py
python tools/prototypes/g14_tiled_unify_trimmed.py
python tools/prototypes/g18_validate_decomposition.py
python tools/prototypes/g19_pinhole_volume_bias.py
python tools/prototypes/g21_straddling_seam_split.py
python tools/prototypes/g22_parallel_face_scan.py
```

G9, G11 and G20 have no standalone script: all three were found and fixed
against the real production rehearsal rather than a synthetic prototype (see
their sections).

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
confirming by how much is what the `TD_HX_rehearsal_test` rehearsal re-run is for
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

---

## G7 — threads or processes for `simplify` / `validate`? ❌ THREADS REJECTED

specification.md §10 ranks parallelising `simplify` (same-domain unification)
and `validate` (`BRepCheck_Analyzer`) as optimization paths 1 and 4 — together
16.5 min of the `cc=5, t=1` rehearsal's 73.1. Every other parallel stage in
this codebase pays a `multiprocessing` spawn-context pool's IPC cost (a `.brep`
file round-trip); for these two stages the payload is the run's entire output,
up to 2 GB, so that round-trip adds a serial master-side read-back that did not
exist before. A `ThreadPoolExecutor` would avoid that entirely — worth checking
before building the heavier path, per docs/algorithm.md §11 ("measured, not
assumed").

Method: 8 unequal solids (702–8,526 faces, echoing the rehearsal's 14
very-unequal ones), built from real instanced-lattice grids of increasing size,
`occ.unify_same_domain` and `occ.is_valid` run over them serially versus
through a 6-worker `ThreadPoolExecutor`.

    python tools/prototypes/g7_thread_scaling.py

| | Serial | Threaded (6 workers) | Ratio |
|---|---|---|---|
| `unify_same_domain` | 5.138 s | 4.697 s | 0.91 |
| `is_valid` | 5.510 s | 5.539 s | 1.01 |

**OCP holds the GIL around both calls.** `is_valid` shows no speedup at all
(ratio 1.01); `unify_same_domain`'s 9 % is well short of anything six real
cores would produce and is consistent with ordinary scheduling noise, not
released-GIL parallelism. Neither clears the 0.6 bar this gate set in advance
for "threads win outright".

### What this decides

Both stages are parallelised with the established process-pool-plus-`.brep`
pattern (:class:`latticegen2.parallel.WorkerPool`), the same mechanism boundary
trimming and the boundary sew already use, not `ThreadPoolExecutor`. The
master-side read-back this implies (each unified solid is read back once, to
run `validate` and `export` against a live shape) is new serial cost that did
not exist before; docs/specification.md §10's re-measurement of the rehearsal
reports what it actually costs.

---

## G8 — can boundary-sew round 2 skip faces that are already fully sewn? ✅ PASS

specification.md §10 ranks a cheaper round 2 (optimization path 3) after
rejecting a hierarchical tree reduction on paper (docs/algorithm.md §8: round
2's cost tracks total face count almost flatly, not shape count, so a tree
would pay that flat cost once per level for nothing). The lever that survives
that argument: only faces bearing a free edge after round 1 can possibly be
affected by round 2 — everything else is already fully joined within its own
tile and can be carried into the final shell unchanged, by reference.

    python tools/prototypes/g8_seam_only_sew.py

Method: a chain of trimmed junctions closed at both ends (so a correct result
is fully watertight either way), tiled, round 1 run, then round 2 computed two
ways — a full sew of every tile's result, and sewing only the free-edge-bearing
subset while carrying the rest through unsewn — and compared for exact
identity, not just plausibility.

| Scale | Seam fraction | Full round 2 | Seam-only round 2 | Face count | Volume diff |
|---|---|---|---|---|---|
| 40 pieces / 20 tiles | 13.5 % | 0.065–0.071 s | 0.015–0.017 s | identical | 2.3e-12 |
| 150 pieces / 75 tiles | 14.1 % | 0.361 s | 0.061 s | identical | 2.5e-11 |

**PASS at both scales**: same total face count, both fully closed (0 free
edges), volume identical to machine precision. Adopted into
`latticegen2.weld._split_seam_interior` / `_sew_round_two`.

**A production-scale bug was found and fixed after this gate passed, and is
recorded here because the gate itself did not catch it.** The first
implementation of the split (committed, then caught before merge) computed
free edges as a plain Python list and tested every face's every edge against
it with `.IsSame()` — `O(faces × edges_per_face × free_edges)`, fine at this
gate's scale (hundreds to low thousands of faces) but not at the rehearsal's:
run against the real `cc=5, t=1` part, the `stitch` stage went from **8 m 57 s
to 51 m 07 s** — a 5.7× *regression*, not the intended improvement. The fix
uses `TopTools_IndexedMapOfShape`, OCCT's own shape-identity map, so every
membership test is the map's own near-`O(1)` lookup instead of a Python-level
scan; re-measured on a synthetic 600-piece / 200-tile chain (16,802 faces),
the split itself costs 0.159 s. The lesson this leaves for future gates at
this project's scale: **a correctness gate on a few hundred faces is not a
performance gate**, and a design meant to run at hundred-thousand-face scale
needs at least one measurement taken there before being trusted, exactly as
docs/algorithm.md §11 already says about correctness bars — the same applies
to complexity, not only to tolerances.

---

## G9 — does the seam-only split hold on real, heavily trimmed production geometry? ❌ FAILED, then repaired

G8 passed on every prototype scale tried, all of them lightly trimmed junctions
(a box far larger than the chain). None of that geometry is heavily trimmed the
way a real part's boundary layer is, and that turned out to be load-bearing.

No standalone script — this was found and fixed against the real `cc=5, t=1`
production rehearsal (docs/specification.md §10) rather than a synthetic
prototype, because reproducing it needs the exact pathology heavy, near-tangential
trimming produces, which no synthetic setup so far has recreated at any scale
tried. Reproduction:

    python src/main.py -i test/TD_HX_rehearsal_test.step -cc 5 -t 1 \
      -o <out>.step --cores 6 --ram 20 -v

**FAILED, first measured 2026-08-16 on commit `f0fa0bb`.** A controlled pair,
identical in every respect but `_split_seam_interior` forced on or off:

| `_split_seam_interior` | open edges at `assemble` | `stitch` |
|---|---|---|
| on (chapter 10 default) | 118,760 | 1 m 20 s |
| off | 10 | 10 m 22 s |

Mechanism: after round 1, a tile's own sewn result can still have a "straddling"
edge — shared between a face round 2's seam-only split keeps (because it also
bears a genuine free edge) and a face the split carries through unchanged
(because within the tile it has none) — measured at 144–720 such edges in every
synthetic 3-D block tried once the check for them existed. Sewing the seam-only
subset without its carried neighbour present lets `BRepBuilderAPI_Sewing`
rebuild the straddling edge onto a new `TopoDS_Edge` object while the carried
face keeps the original one — one shared edge becomes two, each used once, and
`shell_defects`'s shape-identity map (correctly) counts them as two opens.

**Repaired, not reverted.** `_split_seam_interior` itself is untouched — the
guarantee needed is downstream, at `_sew_round_two`. Every interior interface a
correctly sewn boundary layer presents is the whole template cap quad (four
edges, docs/algorithm.md §5.3(b)); every other cap a boundary piece carries
stays a closed face (§7.1). So a component's post-round-2 free-edge count must
equal exactly `4 × its interior interfaces` — no more, no less — and
`weld.sew_boundary` now checks that against `want_rings` (the same dict its
caller already computes for `interface_rings`) and redoes any component that
fails the check on the unsplit tile results, reporting the count as
`SewStats.repaired_components`.

**Re-run after the fix, same commit history plus the repair, 2026-08-16:**
`assemble` reported exactly **2** open edges — at `[1874.836, 60.370, 970.121]`
and `[1874.836, 59.912, 969.775]`, the pre-existing, separately tracked
micro-loop defect (docs/specification.md §10's "residual defect", out of scope
here) — not the seam-split failure. The run log confirms the mechanism fired
exactly where needed: `1 tiled component(s)' boundary-sew round 2 ... left a
free-edge count other than its interior interfaces ... redone with a full
unsplit sew`. `stitch` cost 9 m 54 s — near the no-split control's 10 m 22 s,
since the one affected component holds nearly all of this part's 21,955
pieces, so repairing it pays close to what disabling the split everywhere
would have, while the other 13 small, untiled components (which were never
split in the first place, tiling requiring `MIN_PIECES_TO_TILE`) are
unaffected either way. `test_weld.py`'s
`test_round_two_repairs_a_component_the_seam_split_got_wrong` pins the repair
mechanism itself as a permanent regression, since no scale short of a full
rehearsal reproduces the real defect it stands in for.

---

## G10 — the rehearsal's last 2 open edges: pinhole wires, not debris edges ✅ PASS

With G9's seam-split repair in place, `assemble` on `TD_HX_Indre_Volum` at
`cc=5, t=1` was left with exactly 2 open edges, which
docs/specification.md §10 had tracked since 2026-08-14 as "micron-scale debris
edges" and proposed to fix with `ShapeFix` small-edge removal — the last of its
three candidates, and the only one not yet disproved.

    python tools/prototypes/g10_pinhole_wires.py

Run against the real part at the real node (`(591,-46,-70)`), deliberately with
no synthetic stand-in. Both the diagnosis and the proposed fix were wrong.

### What they actually are

| property | recorded in §10 | measured |
|---|---|---|
| length | 3.171690e-06 / 5.808982e-06 mm | same |
| non-degenerate | yes | yes |
| owning faces | (not recorded) | **1** |
| endpoints | (not recorded) | **do not meet** — 3.17e-06 / 5.81e-06 mm apart |
| position in face | "a planar face carrying 8 edges where 7 would do" | a **1-edge INNER wire** |

Each is an inner wire of a planar face (1.19 and 1.25 mm²) consisting of one
edge whose two endpoints do not meet. It bounds **no area**. It is a pinhole,
not a sliver — which is exactly why it is used by one face, and why
`BRepCheck_Analyzer` calls the solid **valid** while `weld.shell_defects`
rejects it. The defect is visible on the trimmed piece alone, before any
sewing: `shell_defects` on the raw boolean output reports `(2, 0)`.

### Why OCCT's own repairs cannot help, measured

| tool | result |
|---|---|
| `ShapeFix_Wireframe` (small **edges**), prec 1e-5 … 1e-2 | `CheckSmallEdges` finds **0** candidates; 134 edges in, 134 out |
| `ShapeFix_Face.FixSmallAreaWire` (small **wires**), prec 1e-4, 1e-2 | **0** wires removed |

Both want a well-formed closed wire — one to merge an edge into its neighbours,
the other to integrate an area. A single non-closing edge is neither.

### The fix, and what it costs

`occ.remove_pinhole_wires` drops a non-outer wire when every edge in it is
below `PINHOLE_WIRE_TOL` (3e-5 mm) **and** every edge is already used exactly
once. Measured on the piece above:

| quantity | result |
|---|---|
| pinhole wires removed | 2 |
| faces | 32 → 32 (none lost) |
| edges | 134 → 132 |
| surface area drift | **0.000e+00** — bit-identical |
| volume drift | 2.684e-15 |
| cap areas (`resolve_interfaces`' quantity) | unchanged, drift **0.000e+00** |
| `BRepCheck_Analyzer` | still valid |
| `shell_defects` | **(2, 0) → (0, 0)** |

The area result is the load-bearing one and is why `PINHOLE_AREA_TOL` is 1e-12
rather than a comfortable margin: a wire bounding no area *cannot* change the
area of the face carrying it, so unlike same-domain unification (§9) there is no
re-integration and no quadrature noise to allow for. Anything but zero means a
wire that bounded something was removed.

### The safety property, and why the length bar is not load-bearing

Run against the **template solid** — a closed solid, every edge paired — at a
bar of 35.36 mm, larger than every edge in it, so length alone would condemn
all of them: **0 wires removed.** The "already used exactly once" condition is
what does the work. A real feature is paired and is therefore out of reach of
this repair at any threshold, which is a stronger guarantee than a tolerance
gap.

### The lesson, which is the same one G8's postscript records

The pre-fix diagnosis was not vague — it was specific, quantitative, and wrong
in the one property that mattered. An earlier attempt built a synthetic
near-tangential graze reproducing the symptom's scale to four significant
figures (a genuine, non-degenerate ~3e-06 mm edge from a real boolean), swept a
threshold, measured cap-area drift across 25 configurations, and passed. It was
repairing ordinary **two-owner small edges**, which `ShapeFix` removes happily
and which were never the problem. Checking `owners` and endpoint coincidence on
the real part — one inspection, minutes — would have ruled that approach out
before any of it was built.

**A gate must reproduce the defect's mechanism, not its symptom at the right
order of magnitude.** Where the failing geometry is committed to the repo, as
it is here, test against it directly; `test/test_boundary.py`'s pinhole tests do.

---

## G11 — vertices recorded off their edge's curve after sewing ✅ 30 of 34 (rest in G12)

With G10's pinhole repair in place the rehearsal assembles 14 watertight solids
and then fails `validate`: 1 of 14 solids carries **34 individually invalid
faces**. This gate establishes what they are and what fixes them.

Method: dump the faces out of the failed run's staged `unify_0.brep` and
inspect them alone — no synthetic stand-in, per G10's lesson.

### It is 17 bad edges, not 34 bad faces

Each invalid face has exactly **one** individually invalid edge, and the faces
come in pairs because the two faces of a pair **share** that edge (verified by
identical endpoints). 17 edges x 2 owning faces = 34.

The edges are ordinary — spans of 0.183 and 0.356 mm, real `Geom_Ellipse` and
`Geom_BSplineCurve` curves, pcurves present, not degenerate — and every
face-level check passes on planar, cylindrical and B-spline faces alike.

**The fault is a vertex sitting off its edge's 3D curve**, by 2.474044e-05 mm
(ellipse) and 3.316370e-04 mm (B-spline), with that vertex's tolerance inflated
to *exactly* that distance. The validity test therefore sits on the knife edge.

**Not created by the trim.** Every trimmed junction around both sample
locations comes out of `trim_junction` with zero invalid faces and zero invalid
edges, so sewing introduces it — which is why the repair runs on the sewn
boundary rather than in the worker beside the pinhole repair.

### Repairs measured, on the four dumped faces

| repair | fixes | area drift |
|---|---|---|
| `BRepLib.UpdateTolerances_s(shape, True)` | 3 of 4 — fails the B-spline face | 0.000e+00 |
| **`ShapeFix_Edge.FixVertexTolerance(edge, face)`** | **4 of 4** | **0.000e+00** |
| widening the offending vertex by 1 % | 4 of 4 | 0.000e+00 |

`FixVertexTolerance` is adopted: it is OCCT's own repair for this fault, fixes
every case, and adjusts a recorded tolerance rather than geometry — the area
comes back bit-identical, which is the bound `occ.fix_vertex_tolerances`
enforces as a hard failure. `UpdateTolerances` is the obvious alternative and
is rejected on the B-spline case; both fixtures are committed to `test/`
precisely because they discriminate between the two.

### At production scale: 34 -> 4, and the remainder is a different fault

The full rehearsal reports **15 faces retoleranced** (a shared edge's fix
validates its partner face too, so 15 repairs cover ~30 of the 34) and **4
faces remaining**. Those 4 have **no** individually invalid edge or vertex and
pass every face-level check, so their fault is *contextual* — an edge invalid
only in the context of its face, i.e. pcurve against 3D curve. Measured
against them: `FixVertexTolerance` and `FixSameParameter` fix none,
`BRepLib.SameParameter` fixes 3 of 4 at <= 1.7e-09 drift, and
`ShapeFix_Shape.Perform` fixes all 4 but moves geometry by up to **6.4e-04**
relative area and rebuilds faces — which on an already-proven-watertight shell
is the very mechanism behind G9's regression. Not adopted.

**PARTIAL when written**: the fault this gate names is fixed and bounded, and
the rehearsal still failed `validate` on the remaining 4.

**Superseded by G12, which closes them** — and which found the "contextual
pcurve against 3D curve" reading above to be **wrong**. The deviation it
measures is real but is not what the analyzer rejects; the actual fault is a
falsely self-intersecting wire. Read G12 before trusting the last paragraph of
this section.

---

## G12 — the residual 4 faces are falsely self-intersecting wires ✅ FIXED

G11 closed 30 of the 34 invalid faces and named the remaining 4 a *contextual*
fault — "an edge invalid only in the context of its face, i.e. pcurve against
3D curve". **That reading was wrong**, and it is the fourth diagnosis in this
family to be wrong in a way that matched the symptom convincingly. The pcurve
deviation is real, but it is not what `BRepCheck_Analyzer` is rejecting.

Script: [`g12_self_intersecting_wire.py`](g12_self_intersecting_wire.py), run
against the four faces lifted out of the failed run's `unify_0.brep`.

### What the symptom looked like, and why it misled

Measuring every edge's pcurve against its 3D curve does find something: on each
of the four faces exactly one edge — always the fat-tolerance `Geom_BSplineCurve`
the boolean fitted to the strut/input-surface intersection — deviates by
**98–100 % of its own recorded tolerance**:

| face | culprit edge tol | pcurve deviation | dev/tol |
|---|---|---|---|
| `residual_0` (B-spline surface) | 1.539546e-03 | 1.509950e-03 | 0.9808 |
| `residual_1` (cylinder) | 8.741113e-04 | 8.639311e-04 | 0.9884 |
| `residual_2` (cylinder) | 8.743809e-04 | 8.744186e-04 | **1.0000** |
| `residual_3` (cylinder) | 8.744153e-04 | 8.744211e-04 | **1.0000** |

That is the same knife-edge shape as G11's fault, on a different quantity, and
it is exactly the sort of quantitative near-miss that reads as a root cause.
Two measurements disprove it:

* **Widening the culprit edge's tolerance fixes nothing** — not at 1.05x, and
  not at 5x, on any of the four. A knife-edge fault would clear at 1.05x.
* `BRepLib.SameParameter` **does** fix three of them, and on `residual_2` it
  genuinely improves the deviation (8.744e-04 -> 8.296e-04) and lowers the
  tolerance to match — and the face stays invalid. So the deviation is not the
  thing being rejected.

### The mechanism, established with a controlled probe

`BRepCheck_Analyzer.IsValid(subshape)` — the subshape overload, which checks
*in the context of* the analyzed shape, unlike a fresh analyzer on a lone
subshape — points at the **wire**, on all four. The probe is controlled: it
fires before the repair on all four and clears afterwards on the three
`SameParameter` fixes, so a "found nothing" on `residual_2` could not have been
the probe being blind (G10's lesson).

`BRepCheck_Wire::InContext` runs three checks. `Closed` and `Orientation`
return `BRepCheck_NoError`; the fault is **`BRepCheck_SelfIntersectingWire`**,
reported for a pair of edges that is **adjacent in the wire** in every case —
a short, tight-tolerance (1e-07 to 5e-06) trim edge against the fat-tolerance
B-spline intersection edge.

**It is not a real self-intersection.** The two pcurves cross at exactly one
point, and that point lies *at the shared vertex, inside its tolerance*:

| face | pair | intersections | dist to shared vertex | shared vertex tol |
|---|---|---|---|---|
| `residual_0` | e7 ↔ e8 | 1 point, 0 segments | 1.229058e-03 | 1.539646e-03 |
| `residual_1` | e3 ↔ e4 | 1 point, 0 segments | 2.409343e-04 | 8.742113e-04 |
| `residual_2` | e0 ↔ e5 | 1 point, 0 segments | 3.521887e-04 | 8.743809e-04 |
| `residual_3` | e0 ↔ e4 | 1 point, 0 segments | 3.151464e-04 | 8.744153e-04 |

The edges meet where they are supposed to. The shared vertex's recorded
tolerance is simply left a little too tight for OCCT's check to swallow the
crossing — the same *class* of fault as G11 (a tolerance recorded wrong, no
geometry wrong), on a different quantity.

### Which tolerance the check keys on

Sweeping each candidate independently. `c` = check clean, `S` = still
self-intersecting; `V` = face valid. No configuration moved any area.

| target | x1.05 | x1.1 | x1.25 | x1.5 | x2 | x5 |
|---|---|---|---|---|---|---|
| `residual_0` shared vertex | cV | cV | cV | cV | cV | cV |
| `residual_1` shared vertex | cV | cV | cV | cV | cV | cV |
| `residual_2` shared vertex | S. | S. | **cV** | cV | cV | cV |
| `residual_3` shared vertex | S. | **cV** | cV | cV | cV | cV |
| *all four*, fat edge instead | S. | S. | S. | S. | S. | S. |

The shared vertex is the knob; the fat edge is not, at any factor. **This also
explains `residual_2`, which is what the whole investigation turned on:**
`SameParameter` fixes a face only incidentally, by *raising* the culprit edge's
tolerance to 1.05x its deviation, which propagates to the shared vertex. On
`residual_1` and `residual_3` that raised the vertex tolerance (to 9.071e-04
and 9.163e-04) and cleared the check as a side effect. On `residual_2` the
re-fit *lowered* the edge tolerance instead, so the vertex was never widened
and the check never cleared.

### Repairs measured, on all four faces

| repair | fixes | area drift | topology replaced |
|---|---|---|---|
| `ShapeFix_Wire.FixSelfIntersection` | 0 of 4 | 0.000e+00 | none — a no-op |
| `BRepLib.SameParameter` | 3 of 4 — not `residual_2` | ≤ 1.7e-09 | none |
| `ShapeFix_Shape.Perform` | 4 of 4 | up to **6.382e-04** | **2–3 edges/face** |
| **widen the shared vertex** | **4 of 4** | **0.000e+00** | **none** |

`ShapeFix_Wire.FixSelfIntersection` is the tool the fault's name points at and
is a complete no-op here, at either of its modes — the fourth time in this
family that OCCT's named repair for the named symptom does not touch the actual
defect. `ShapeFix_Shape` works and is rejected on two independent grounds: it
moves geometry, and it rebuilds faces, keeping only 6 of 9 and 4 of 6 edge
objects. Minting new edges on a shell `assemble` has already proven watertight
is precisely G9's regression mechanism.

Widening the shared vertex is adopted. It is metadata-only, so surface area
comes back **bit-identical** and every `TopoDS_Edge` and `TopoDS_Vertex` object
survives — both asserted as permanent regressions in `test_weld.py`, the second
being the property that makes this safe on a proven shell. It is also
monotonically *permissive*: every check reading a vertex tolerance is a "within
tolerance" test, so a neighbouring face sharing the vertex can only become more
valid. Measured rather than argued — widening *every* vertex of each repaired
face to 1e-01 mm, twenty-five times the repair's own absolute cap, leaves all
four valid with area still bit-identical.

`occ._widen_self_intersection_vertices` searches for the smallest widening that
satisfies **OCCT's own predicate** — it asks `BRepCheck_Wire::SelfIntersect`
rather than re-deriving the rule OCCT applies — bounded at
`SELF_INTERSECT_TOL_GROWTH` (4x what the kernel itself recorded) and
`SELF_INTERSECT_MAX_VERTEX_TOL` (4e-3 mm, a hundredfold below the CLI's
smallest legal strut). Measured need is at most 1.25x / 1.093e-03 mm, so both
bounds clear it comfortably; a face that exhausts them is left in
`still_invalid` for `validate` to report, never widened without limit.

`residual_2` and `residual_0` are committed as `test/self-intersecting-wire-
{cylinder,bspline}.brep`. Both are kept for the same reason G11 kept two:
they discriminate between the candidates — `SameParameter` fixes the B-spline
one and not the cylinder one.

---

## G13 — can same-domain unification be tiled *below* the solid? ⚠️ CONDITIONAL PASS

docs/specification.md §10 records optimization path 1 — parallelising
`simplify` across the shared pool — as implemented, correct, and **not a
wall-clock win** on the `cc=5, t=1` rehearsal (17 m 17 s → 18 m 39 s, still
0.99 cores), because that part's 14 solids are one dominant body plus 13
scraps and "the largest single solid is the floor, not the sum". Dispatching
body-for-body cannot lower that floor; tiling *within* a solid is the only
lever that can. Two questions decide whether it is worth building.

    python tools/prototypes/g13_unify_scaling.py --tiles 2 3 4

> **Re-running this gate now gives different numbers, by design.** Every figure
> below was measured *before* docs/specification.md §10 Phase 2, when the
> interior was instanced with one lateral face per half-strut. Phase 2 builds
> the merged full-strut face directly, so the grids this script constructs now
> arrive at roughly the face count unification used to produce — which is
> exactly the change Phase 2 made, and it means part A's inputs are no longer
> the ones tabulated here. The scaling exponent and the tiling identity results
> still stand as measurements of `ShapeUpgrade_UnifySameDomain` itself.

Method: closed all-interior `m × m × m` grids at `cc=10, t=1.5` (the same
family G7 used), `m` from 4 to 16. Part A times a whole-solid unification at
each scale. Part B partitions the largest two solids' faces into an `n³`
spatial grid **by face centroid**, unifies each tile's face set on its own,
concatenates the results with `BRep_Builder.Add` and counts edges used by
exactly one face. The input is closed, so "did the tile seams survive" is a
count against zero, not a judgement.

### A — cost against face count: mildly superlinear, and worsening

| m | faces in | unify s | ms/face | faces out |
|---|---|---|---|---|
| 4 | 1,632 | 0.209 | 0.128 | 1,056 |
| 8 | 12,672 | 2.197 | 0.173 | 7,296 |
| 12 | 42,336 | 7.996 | 0.189 | 23,328 |
| 14 | 67,032 | 13.831 | 0.206 | 36,456 |
| 16 | 99,840 | 22.968 | 0.230 | 53,760 |

Overall log-log slope **1.135**, just under this gate's 1.15 "superlinear"
bar — but the local slopes climb with scale (1.044 at 12 k→25 k faces, 1.193
at 42 k→67 k, **1.273** at 67 k→100 k) and cost per face rises monotonically
from 0.128 to 0.230 ms. So: **plan on tiling being worth about `W` and no
more.** The serial sums in part B agree — 8 tiles cost 6.494 s against the
whole solid's 7.054 s, an 8 % saving before any parallelism — so essentially
all of the available win is parallelism, not a smaller `n^k` term. The
upward trend does mean the extrapolation to a 1 M-face solid is not linear,
and it is the reason the per-face rates in G7 (0.17 ms/face at ≤8.5 k faces)
and the rehearsal (1.11 ms/face at 1 M) cannot be reconciled by assuming one.

### B — tiles reassemble without sewing, but only with `unify_edges=False`

| Scale | `unify_edges` | tiles | whole | slowest tile | faces | free edges | valid |
|---|---|---|---|---|---|---|---|
| m=16 | True | 8 | 23.007 s | 2.361 s | +4.94 % | **3,632** | — |
| m=16 | True | 27 | 23.007 s | 2.270 s | +2.14 % | **11,232** | — |
| m=16 | False | 8 | 7.054 s | 0.908 s | +4.94 % | **0** | yes |
| m=16 | False | 27 | 7.054 s | 0.866 s | +2.14 % | **0** | yes |
| m=16 | False | 64 | 7.054 s | 0.374 s | +7.17 % | **0** | yes |

**With edge unification on, tiling is impossible.** Thousands of seam edges
come back as two distinct objects and the reassembled shell is full of holes.
That is the edge pass doing exactly what it is for: concatenating the
collinear pairs left inside a merged wire rewrites edges *on* the tile
boundary, so the two sides stop being the same `TShape`.

**With it off, identity is exact.** Zero free edges at every tile count and
both scales, `BRepCheck_Analyzer` valid, volume preserved to ~1e-13 relative.
`BRep_Builder.Add` alone suffices — no sewing ever touches the
volume-scaling face set, which is what docs/algorithm.md §6 and §8 require.

### The finding that was not being looked for

**`unify_edges=False` is 3.1–3.3× faster and merges exactly the same faces.**
At m=16: 23.007 s → 7.054 s, both producing 53,760 faces. At m=14: 14.183 s →
4.512 s, both 36,456 faces. The edge pass is ~70 % of the stage's cost and
contributes nothing to face merging.

It is **not** free to drop, and docs/algorithm.md §9's "worth almost nothing"
(4 edges out of 81,816, measured on the 80 mm ball) does not hold at lattice
scale: at m=16 the edge pass takes 307,200 edges down to 215,040, a **30 %
reduction**.

### Follow-up on a real part: neither the drop nor the split is a speed win

The trade above was measured on `dense-lattice` before being taken, and **both
candidate forms failed** (docs/specification.md §10, Phase 1):

| `simplify` on `dense-lattice` | time | output |
|---|---|---|
| combined, one call (baseline) | 13.21 s, 16.18 s | 67,898 edges, 52.80 MB |
| edge pass dropped | 9.45 s | 94,476 edges, 71.29 MB |
| split: faces, then edges alone | 13.87 s, 13.94 s | 67,898 edges, 52.80 MB |

**Dropping the edge pass fails on the downstream cost it creates.** The 3.8 s
it saves is handed straight back to `validate` (6.24 → 8.21 s) and `export`
(6.25 → 10.50 s), which scale with edge count as well as face count, and the
file grows 35 %. Net run time does not change (57.28 → 57.57 s).

**Splitting the passes is neutral, not cheap.** The hypothesis that the edge
pass would be far cheaper over an already-face-merged solid is **disproved**:
it costs ~4.4 s either way. The 3.1–3.3× above therefore does not transfer from
the synthetic grid to a real trimmed solid, where the edge pass is ~33 % of the
stage rather than ~70 %.

The split was nonetheless adopted, on structural grounds rather than speed: it
is what allows the face merge to run with edge merging **off**, which is a
precondition for the tiling below, and it improves degradation (a throwing edge
pass no longer discards a completed face merge). `test/test_pipeline.py` pins
the split's B-rep as identical to the combined call's.

**The lesson for the rest of this gate.** Part A's exponent and the tiling
ceilings below are measured on synthetic all-planar instanced grids. This
follow-up shows that a ratio measured there can be off by 2× on a real part
carrying trimmed curved faces. Treat every projection derived from part A as an
upper bound until it is confirmed on a real one.

### Merge loss, and why the partition should not be generic

A merge group straddling a tile seam merges partially in each tile, so the
tiled result carries more faces than a whole-solid unification: **+2.1 % to
+11.4 %** here, non-monotonic in tile count because it depends on where the
centroid grid falls relative to the lattice. That is a "the output is slightly
larger" failure mode, which docs/algorithm.md §11 accepts — but it need not be
paid at all. The centroid partition used here is deliberately the *generic*
one, so its loss is an upper bound: the merge pairs are known by construction
(one per surviving mid-strut interface), so bucketing by **strut** rather than
by face centroid gives a partition no merge group can straddle, and a loss of
zero.

### What this decides

* Sub-body tiling is viable, and `unify_edges=False` is a **precondition**,
  not an option.
* Expect ~`W`, not more. The mechanism must therefore not add serial cost
  comparable to what it saves — reassembly and `.brep` IPC included.
* Tile by strut, not by centroid, so the merge loss is zero by construction.
* The edge-merge trade above is a separate decision from tiling, worth more
  than tiling is on its own, and needs one measurement on a real part.

---

## G14 — does tiled unification survive *trimmed, curved* faces? ✅ PASS

G13 proved that unifying a solid's faces in spatial tiles and concatenating the
results with `BRep_Builder.Add` reproduces a closed valid solid with zero free
edges, provided edge merging is off — the property docs/specification.md §10
Phase 3 rests on. But it measured that on **all-planar instanced grids only**,
which Phase 3 records as blocking risk R2: real output carries faces trimmed by
`BRepAlgoAPI_Common` against the input body, where `TShape` identity across a
tile seam was unproven.

Phase 2 sharpened the question rather than retiring it. With the interior now
built pre-merged (docs/algorithm.md §6), the faces still entering `simplify` are
mostly boundary-derived — on `dense-lattice`, 15,718 of 25,234 — so the region
Phase 3 would tile is exactly the region G13 never tested.

    python tools/prototypes/g14_tiled_unify_trimmed.py --tiles 2 3 4 --pieces path/to/temp/<stamp>

Two inputs, both carrying genuinely curved trimmed faces:

**(a) A closed solid: instanced grid ∩ sphere** — 1,804 faces, **76 curved**,
valid, 4590.670414 mm³.

| tiles | free edges | faces | valid | volume drift |
|---|---|---|---|---|
| whole solid | 0 | 1,800 | yes | 6.06e-12 |
| 8 | **0** | 1,804 (+0.22 %) | yes | 1.19e-15 |
| 27 | **0** | 1,800 (+0.00 %) | yes | 6.06e-12 |
| 64 | **0** | 1,804 (+0.22 %) | yes | 1.39e-15 |

**(b) Real trimmed boundary pieces** from a kept `temp/<stamp>` of an 80 mm ball
run, sewn as `weld` sews them — 2,342 faces, **176 curved**. This shell is open
at every interface hole, so the bar is that tiling introduces no *new* free
edge, not that there are none.

| tiles | free edges (baseline 0) | faces | area drift |
|---|---|---|---|
| 8 | **0** | 2,342 | 5.31e-10 |
| 26 | **0** | 2,328 | 3.69e-16 |
| 60 | **0** | 2,342 | 5.31e-10 |

**PASS on both.** Curved trimmed faces do not break tile identity, so Phase 3
need not fall back to tiling only the interior-derived faces.

> **Superseded by G15, and this verdict must not be read as a green light.**
> Every measurement in this gate — like G13's — reassembles tiles *in one
> process*, where identity is pointer identity. Phase 3 dispatches tiles to
> separate worker processes as `.brep` files, and G15 measures the seam not
> surviving that. What G14 establishes is narrower than it reads: curved
> trimmed faces are not what breaks tiling. Something else is.

### Two things this gate does *not* establish

**The amount of merging in it is small.** Whole-solid unification took (a) from
1,804 faces to 1,800 and (b) from 2,342 to 2,328 — because Phase 2 already
merged the interior before this gate ever sees it. Identity is only really
stressed where unification *rebuilds* faces, so the strongest evidence for
"identity survives heavy merging" remains G13's planar case (99,840 → 53,760).
Together they cover the two axes separately rather than jointly; a part that
merges heavily *and* is mostly curved is not represented here.

**The bar had to be corrected mid-gate, and the first version was wrong in a
familiar way.** Volume/area preservation was initially barred at 1e-12, carried
over from G13's planar measurements, and it failed a *correct* tiling at an area
drift of 5.31e-10. Quadrature over curved trimmed faces is not exact —
docs/algorithm.md §9 already records 2.4e-7 relative drift from the same cause
on `dense-lattice` — so the bar now reuses `pipeline.UNIFY_VOLUME_TOL` (1e-5),
the figure production already applies to this exact operation. This is the same
mistake docs/specification.md §10 records against the boundary-unification
attempt, made again; the tolerance-free part of the gate is the free-edge
count, and that is what the verdict should rest on.

A first run was also **vacuous** and said so: the sphere radius was scaled from
the bounding box *diagonal*, and the lattice box is elongated along Z, so the
sphere contained the whole grid and trimmed nothing. The curved-face count is
now checked before the result is believed.

---

## G15 — does tile identity survive the `.brep` round trip? ❌ FAIL

G13 and G14 both concluded that a tiled same-domain unification reassembles with
`BRep_Builder.Add` alone — 0 free edges, valid, volume preserved — because
`ShapeUpgrade_UnifySameDomain` leaves the edges it did not merge as *the same
objects*, so a tile seam stays one `TShape`. That is the property
docs/specification.md §10 Phase 3 is built on.

**Both gates measured it in one process, where identity is pointer identity.**
Phase 3 does not run in one process. Its entire purpose is to spread tiles
across `latticegen2.parallel.WorkerPool`, and every stage that does so moves
geometry as `.brep` files (docs/algorithm.md §7, §8). A `.brep` is a
serialization: sharing is preserved *within* one file and cannot be preserved
*between* two, because each file writes its own copy of every edge it
references. Two tiles sharing a seam edge are written by two different workers
into two different files.

    python tools/prototypes/g15_tiled_unify_ipc.py

Input: G14's instanced grid ∩ sphere — 1,804 faces, **76 curved**, closed, so
the whole-solid unification has 0 free edges by construction and the test is a
count rather than a judgement. Three reassembly routes, the first two controls:

| route | 8 tiles | 27 tiles | |
|---|---|---|---|
| in-process (G14's route) | **0** | **0** | control: reproduces G14 |
| one file, all tiles as one compound | **0** | **0** | control: serialization itself is fine |
| **one file per tile** | **864** | **1,760** | what Phase 3 actually does |

**The middle row is what makes this readable.** Serialization does not destroy
sharing — a compound written to one file and read back keeps every seam. What
destroys it is the *file boundary*, which is precisely the boundary Phase 3
puts between workers. So this is not a `.brep` defect to work around; it is
what one-process-per-tile means.

### Why the obvious repairs do not rescue it

* **Re-identify the duplicates on the master.** The geometry is unchanged, so
  matching the pairs is a cheap exact lookup. *Merging* them is not: it needs
  `BRepTools_ReShape` to replace edges **and** their vertices, and
  docs/algorithm.md §8 already measured that — replacing an edge's vertices
  leaves the neighbouring edges pointing at the old ones, the wire comes apart
  (`BRepCheck_NotConnected`), the solid is invalid and the volume wrong, while
  every edge still has two faces and the shell still "closes".
* **Sew only the seam-bearing faces.** This is G8's split, and its failure mode
  is documented at production scale: on the rehearsal's real trimmed pieces it
  produced 118,760 open edges where 10 were expected, and the checked fallback
  is a full unsplit sew. A full sew over a volume-scaling face set is the one
  thing docs/algorithm.md §6 and §8 exist to prevent (it cost 4 h 45 m of a
  5 h 04 m run). The failure mode here would be "catastrophically more work",
  not §11's acceptable "a little more work".
* **Tile inside a single worker.** Sound — the middle row proves it — but then
  the tiles are serial, and G13 measured serial tiling at an 8 % saving
  (6.494 s against 7.054 s at 8 tiles). G13's own conclusion was "plan on tiling
  being worth about `W` and no more; the win is parallelism". Without the
  parallelism there is close to nothing left.

**The lesson is G8's, in a new place.** A property proven at prototype scale is
proven *under the prototype's conditions*, and "runs in one process" was a
condition neither G13 nor G14 stated because neither had a reason to. Phase 3's
risk list carried R2 (curved faces), R3 (merge loss), R5 (memory/IPC volume) and
R6 (the `IsSame` trap) — five gates' worth of scepticism — and none of them
asked whether the mechanism survives the process boundary the plan's first
sentence puts it across.

---

## G16 — is unification priced by the size of its input? ✅ YES (and that is the surprise)

docs/specification.md §11 records two attempts to make `simplify` cheaper by
giving `ShapeUpgrade_UnifySameDomain` less to look at. Phase 2 cut its input
~31 % and the stage fell 12.4 %. The restricted face merge cut it **46 %** and
the stage did not fall at all, on a byte-identical output.

The natural conclusion — "the call is priced by what it emits, not what it
consumes" — was written into this project's docs twice before being tested.
**It is wrong.**

    python tools/prototypes/g16_unify_elasticity.py

On G14's trimmed test solid (4,350 faces, 270 curved), timing the face merge
over spatially coherent subsets, with `elasticity = (relative time saved) /
(relative input removed)`:

| input | share | time | saved | elasticity |
|---|---|---|---|---|
| 4,350 | 1.00 | 0.156 s | — | — |
| 3,915 | 0.90 | 0.141 s | 9.5 % | 0.95 |
| 3,480 | 0.80 | 0.127 s | 18.9 % | 0.94 |
| 2,610 | 0.60 | 0.094 s | 39.7 % | 0.99 |
| 1,740 | 0.40 | 0.057 s | 63.3 % | 1.06 |

**Mean 0.98 — very nearly linear.** Remove a representative 40 % of the faces
and you save 40 % of the time.

### What that changes

The pipeline's own measurement stands: on a real `dense-lattice` solid,
removing 20 % of the faces saved **6 %** of the time, an elasticity near 0.3.
Against 0.98 for a generic removal of the same size, the gap is not the
kernel's pricing — it is **which** faces a correct restriction removes.

A correct restriction skips exactly the faces unification would have returned
unchanged, because those are the only ones it is *allowed* to skip. Those are
also the cheap ones. It keeps exactly the faces that merge, which are the
expensive ones. **The property that makes the restriction correct is the
property that makes it worthless**, so no implementation improves on 0.3, and
the ~0.045 ms/face of bookkeeping any restriction needs is never recovered.

That is a stronger and more transferable claim than "the kernel is priced by
its output", and it is the one to carry forward: *any* future proposal to feed
this stage less inherits the same 0.3, whatever mechanism it uses to identify
what to skip.

### What this gate does not measure

The **targeted** elasticity itself. Doing that needs geometry that merges
heavily, and the trimmed test solid here merges 4,350 → 4,344 — six faces —
now that Phase 2 builds the interior pre-merged (docs/algorithm.md §6). The
generic figure is the honest half of the comparison, and the targeted one is
the pipeline measurement above, taken on `dense-lattice` where the real merge
is 25,234 → 15,966. Quoting this gate's 0.98 as if it were what a restriction
would achieve would repeat exactly the mistake it exists to correct.

---

## G17 — threads instead of processes for a tiled unification? ❌ NO

G15 killed sub-body tiling on the *transport*: unified tiles reassemble by
shared topology only inside one process, and shipping them to workers as
`.brep` files gives every seam edge two copies. Threads would sidestep that
entirely — one heap, tiles pointing at literally the same `TShape`, no
serialization anywhere — so the honest follow-up is whether the transport was
ever the real obstacle.

    python tools/prototypes/g17_thread_tiled_unify.py --m 12 --threads 6

Input: a closed 12×12×12 instanced grid, 23,328 faces, 0 free edges.

**Probe A — is the GIL released during the call?** A Python counter spins in a
background thread while one `unify_same_domain` runs on the main thread, and its
rate is compared against the same counter with nothing competing:

| | iterations/s |
|---|---|
| counter alone | 8,369,062 |
| counter during a 0.87 s OCCT call | 312,934 |

**3.7 % retained — the GIL is held for essentially the whole call.**

**Probe B — does threading go faster anyway?** 27 tiles, serial against 6
threads:

| | time | faces | free edges after reassembly |
|---|---|---|---|
| serial | 0.935 s | 23,328 | — |
| 6 threads | 0.896 s | 23,328 | **0** |

**1.04×.** No parallelism, exactly as probe A predicts.

### The finding is the pair of results, not either one

**Threads fix the thing processes break, and break the thing processes fix.**
Reassembly after threaded tiling leaves **0 free edges** — identity is perfect,
because nothing was ever serialized. And the tiles run one at a time regardless
of how many threads dispatch them. Processes give real parallelism (`boundary`
reaches 5.2 of 6 cores) and destroy tile identity. There is no third option in
this architecture, so sub-body parallelism in `simplify` is closed, not merely
unbuilt.

This also upgrades G7's finding from symptom to mechanism. G7 inferred a held
GIL from a 0.91–1.01× speedup; probe A observes the interpreter stalling
directly, which is the thing a scaling measurement can only imply.

### The escape hatches, and why none is available

* **No internal parallel mode.** `ShapeUpgrade_UnifySameDomain` exposes no
  `SetRunParallel` or thread-pool hook (checked: its whole surface is
  `AllowInternalEdges`, `Build`, `History`, `Initialize`, `KeepShape(s)`,
  `SetAngularTolerance`, `SetLinearTolerance`, `SetSafeInputMode`, `Shape`), so
  unlike OCCT's booleans and mesher there is nothing to switch on.
* **A C extension releasing the GIL around the call** would work in principle
  and is ruled out by packaging: specification.md §2 requires ordinary wheels,
  no build step and no compiler on the target.
* **Free-threaded CPython (PEP 703)** is the only thing that would change the
  answer. The project pins Python 3.11 with pinned wheels, and no free-threaded
  `cadquery-ocp` build exists.
* **Shared memory instead of files** does not help. It addresses I/O, which was
  never the cost — `simplify`'s round trip measured 940 MB against an 18-minute
  stage. A `TopoDS_Shape` is a graph of pointer-linked C++ objects with
  vtables; another process cannot use it at a different base address, which is
  why `.brep` serialization is the only supported transfer and why the identity
  loss is structural rather than an artefact of choosing files.

### What this gate does *not* establish

**That threading would be safe if it were fast.** Tiles sharing a seam edge
share a `TShape`, so two threads unifying neighbouring tiles touch the same
object and the same reference count concurrently — exactly the case OCCT's
thread-safety notes do not cover. Probe B ran without crashing, and that is not
evidence: a data race that does not fire is indistinguishable from no race.
The question simply never becomes worth asking, because probe A closes it first.

---

## G18 — can `validate` use more than one core, and at what risk?

`tools/prototypes/g18_validate_decomposition.py`. Run on a real trimmed lattice
solid (`test-cylinder` at `cc=10, t=1.5`) plus the four real invalid faces
committed from the `cc=5, t=1` rehearsal.

**Asked because `validate` and `simplify` are not the same problem, though
docs/testing.md long described them together.** Per-*body* dispatch failed for
both for one reason — one dominant solid sets the floor. Going *below* the body
is where they diverge: `simplify` must hand back geometry, so its tiles have to
reassemble by shared topology and G15's file boundary destroys that.
`_worker_validate` handed back `(bool, float)`. Nothing reassembles, so G15 has
nothing to attach to.

### The answer was not the split this gate was written to evaluate

`BRepCheck_Analyzer` has a parallel flag of its own:

    BRepCheck_Analyzer(S, GeomControls=True, theIsParallel=False, theIsExact=False)

G17 recorded that `ShapeUpgrade_UnifySameDomain` "exposes no `SetRunParallel` or
thread-pool hook, so unlike OCCT's booleans and mesher there is nothing to
switch on". **This one has the hook.** The threads are OCCT's own native ones,
so the GIL result does not bind them, and the verdict stays OCCT's rather than a
conjunction assembled here.

| solid | serial | `theIsParallel` | speedup | cores |
|---|---|---|---|---|
| 14,790 faces | 3.791 s | 2.370 s | **1.60×** | 0.96 → **3.43** |
| 1,176 faces | 0.283 s | 0.179 s | **1.58×** | 0.94 → **3.93** |

**A — verdict safety, the part that decides it.** All four committed invalid
faces (`invalid-vertex-tolerance-{ellipse,bspline}.brep`,
`self-intersecting-wire-{cylinder,bspline}.brep` — the real geometry behind
docs/algorithm.md §8's two repair rungs) read `False` on **both** sides, and
both valid solids read `True` on both. This control is the point of the gate,
not the timings: `validate` is the exit-4 gate before a 2 GB file ships, and per
G10 a scanner that only ever agrees on sound geometry proves nothing. Pinned as
a permanent regression in `test_pipeline.py`.

### The manual split would go further, and was not taken

| | share of the serial check |
|---|---|
| per-face checks, summed | **94.4 %** |
| structural only (shell closure + orientation) | 4.8 % |

So ~94 % of the work is local to a face and would divide across `W` workers —
roughly 5× against the flag's 1.6×. **Not built**, by decision, and the reason
is not effort: `BRepCheck_Analyzer(solid)` checks subshapes **in context**, and
a standalone per-face check is a different predicate — that difference is
precisely what G12 exploited to find rung 2's mechanism. Replacing this
project's final correctness gate with a hand-assembled conjunction needs a
control proving it catches an in-context-only fault, and no such control exists
cheaply. Under docs/algorithm.md §11 the failure mode of an optimization must be
"do more work", never "produce a wrong result"; this one's would be the latter.
The 94.4 % figure is recorded so the option stays visible if the stage ever
justifies building that control.

### What shipped

The flag, plus moving `validate` off the worker pool onto the master
(docs/algorithm.md §9). The two are one decision, not two: `--cores` is
"honoured exactly" (specification.md §3), and `W` worker processes each
launching `W` OCCT threads is `W²` threads on `W` cores.
`latticegen2.parallel.set_thread_budget` caps OCCT's own pool to the budget.
Running on the master also deletes a `.brep` round trip that existed only to
reach the workers — 464 MB each way on the rehearsal, to compute two scalars per
solid.

Measured as a controlled pair on `dense-lattice`, back to back:
`validate` **5.49 s → 2.11 s (−61.6 %)**, better than the flag alone because the
round trip goes too. Output byte-identical.

**What this gives up**, stated plainly: on a part whose components *are* evenly
sized, per-solid dispatch across `W` processes would beat 1.6×. That was
docs/specification.md §10 path 4's justification for keeping a change that
measured *slower* (2 m 59 s → 3 m 29.6 s) on the only part it was ever run on.
No evenly-sized part has been measured, so what is traded away here is a
projected benefit against a measured one.
## G19 — the pinhole repair's volume guard refused a valid run ✅ FIXED

`-cc 12 -t 2.5` on `TD_HX_rehearsal_test` failed with

    FAILED: Removing 1 pinhole wire(s) from the junction at (2036.89, 60.0,
    969.998) changed its volume from 77.3988208 to 77.3988207 mm^3
    (1.235e-09 relative, tolerance 1e-09).

Both parameters are inside the CLI's ranges and satisfy `t < a`, so this is
docs/algorithm.md §11's worst failure mode — refusing correct input — and it had
been there since the repair was written.

    python tools/prototypes/g19_pinhole_volume_bias.py

### The three readings that fit the symptom, and why each is wrong

The absolute change is 9.6e-08 mm³ on a 77.4 mm³ piece, which reads as noise of
some kind. It is not.

| reading | test | result |
|---|---|---|
| quadrature noise | adaptive Gauss-Kronrod at eps 1e-3 … **1e-11** | drift **9.555971e-08** at every eps — the two shapes integrate to genuinely different values |
| a coordinate-magnitude artifact (the part sits 2,257 mm from the origin) | properties taken about the origin / the node / 100× the node; piece translated to the origin and 10× further out | **9.55597e-08 throughout** — location-invariant |
| a scaling the bar got wrong for junction size or `t` | compare against the `cc=5, t=1` junction the bar was calibrated on | the cc=12 wire is **shorter** (1.010641e-06 mm against 3.171691e-06 and 5.808982e-06) and drifts **seven orders further** (1.235e-09 against 3.097e-15) |

The third is the one worth dwelling on, because "the bar is mis-scaled for a
larger junction" is the natural diagnosis and it is quantitatively hopeless: the
piece is 9× the volume and the drift is 4×10⁵ times larger, on a shorter wire.

### What it actually is

`BRepGProp::VolumeProperties` states that its shape "must be exempt of any free
boundary … results will be false if [the conditions] are not respected". **A
pinhole wire is a free boundary** — an edge used by exactly one face is the
definition of the defect (G10). So the volume of the *unrepaired* piece is
computed outside OCCT's stated precondition and the repaired one inside it. The
guard was comparing a figure OCCT does not promise against one it does, and
calling the difference geometry loss.

Surface area is untouched by this, at 172.14799496928404 mm² **bit-identical**
before and after — which is exactly what docs/algorithm.md §7 argues must
happen, and it is why the area check was always the one doing the work.

### The control: add a free boundary to a clean face

If a free boundary biases the integral, adding one to a repaired piece must
reproduce the bias. A synthetic open inner wire — one edge, with a pcurve, on a
planar face, nothing else changed — was added to four faces of the repaired
piece at four lengths each:

| face | area mm² | dV at L=1e-6 | at L=1e-5 | at L=1e-4 | at L=1e-3 | dA |
|---|---|---|---|---|---|---|
| 0 | 7.3207 | +9.665055e-08 | +9.665057e-08 | +9.665058e-08 | +9.665089e-08 | **0** |
| 1 | 6.4368 | +2.488447e-07 | +2.488447e-07 | +2.488446e-07 | +2.488443e-07 | **0** |
| 2 | 7.1440 | −2.413901e-08 | −2.413948e-08 | −2.414404e-08 | −2.418949e-08 | **0** |
| 3 | 3.4968 | +2.549071e-08 | +2.549102e-08 | +2.549426e-08 | +2.552665e-08 | **0** |

and, on the face that carried the real pinhole, +9.471e-08 against the real
defect's 9.556e-08 — **the same footprint to within 1 %**.

Three things follow, and together they close the question of what bar to use:

* **Area is exactly blind to it**, at every length and position tried (a further
  15-point sweep of the wire's position over ±0.5 in UV: `dA = 0.0` throughout).
* **The magnitude is not the wire's.** A thousandfold change in length moves it
  by at most 0.21 %, where a length-proportional term would move it by 1000×.
* **The magnitude is the face's**, varying 10.3× across four faces of one piece
  (3.12e-10 to 3.22e-09 relative) and changing sign between them.

So there is nothing the repair controls to express a tolerance in terms of, and
the 2.7e-15 the original bar was calibrated on was not headroom — it was two
wires that happened to sit on faces contributing almost nothing.

### The fix: prove it structurally instead

`occ.only_inner_wires_dropped(before, after, n_removed)` replaces the bar and
returns a reason or `None`. It requires the face count unchanged; every face's
**outer wire to be the same object** as some face's before, one for one; the
paired faces to have the same orientation and bit-identical area; every wire
kept to have been a wire of its counterpart; and the wires lost to number
exactly `n_removed`. That pins the enclosed region object-for-object, which is
a question with an exact answer where volume is one OCCT can only answer with a
bias while the defect is present.

It is cheap because `BRepTools_ReShape` rebuilds only what it touches: **27 of
28** faces come back `IsSame`, so this is a walk over objects plus one area
evaluation per face.

| | before the fix | after |
|---|---|---|
| `cc=5, t=1`, node (591,−46,−70) | passes | passes, `only_inner_wires_dropped` clean |
| `cc=12, t=2.5`, node (262,−27,−37) | **exit 4** | passes, `only_inner_wires_dropped` clean |
| a deliberately wrong `n_removed` | — | rejected |
| the repair run backwards (a wire *gained*) | — | rejected |

A serial census of all **2,907** boundary junctions at `cc=12, t=2.5` finds
exactly **one** pinhole-bearing piece — the one that failed — so this defect is
rare, and a run either clears it or loses everything.

### The lesson, which is the fourth of its kind here

G9–G12 each turned on a repair being aimed at the wrong object. This one turns
on a *check* being aimed at a quantity that is not defined on the shape it is
handed. docs/algorithm.md §5.1 already says a gate is only as trustworthy as
the tightness of the quantity it compares; this adds that it is only as
trustworthy as whether that quantity means anything on the input. **Where a
check can be made structural — the same objects, in the same places — it should
be.**

---

## G20 — the unification volume bar refused a valid run too ✅ FIXED

Found immediately after G19, on the same command, because clearing the pinhole
guard let `-cc 12 -t 2.5` reach `simplify` for the first time:

    FAILED: Same-domain unification changed a solid's volume from 181.005368 to
    181.007867 mm^3 (1.38e-05 relative, tolerance 1e-05).

No standalone script — like G9 and G11, this was measured against the real run's
kept temp directory (`unify_*.brep` and `unify_*_out.brep`, the exact shapes the
workers unified), which is more faithful than any reproduction.

### Nothing moved, and that is measured three ways

| test | result |
|---|---|
| exact symmetric difference, `BRepAlgoAPI_Cut` both ways | **0.000000000e+00 mm³** each way |
| `BRepCheck_Analyzer` | valid before **and** after |
| whole-run output | 9 solids, unchanged everywhere else |

### But it is not integrator noise either, and that is the surprise

The obvious reading — quadrature over larger merged regions — predicts that a
better integrator converges the two figures. It does not:

| | before | after | drift |
|---|---|---|---|
| default integrator | 181.00536767629785 | 181.0078668508196 | 1.3807e-05 |
| Gauss-Kronrod, eps 1e-6 | 181.0053536917804 | 181.0078507197767 | 1.3795e-05 |
| Gauss-Kronrod, eps **1e-11** | 181.00535226444575 | 181.00784910601723 | **1.3794e-05** |

Surface area moves too — 359.32946624725827 → 359.3283313326181, 3.16e-06
relative. So unification really does re-describe the boundary slightly; it does
not merely re-integrate it. (Contrast G19, where area was bit-identical and the
volume difference was an artifact of the *defect*, not of the repair.)

### Sized as a displacement, it is 7 nanometres

`|ΔV| / surface area` is the boundary movement the drift implies, and it is the
quantity that makes the nine solids of this run comparable:

| solid | volume mm³ | area mm² | A/V | rel. volume drift | `|ΔV|/A` mm |
|---|---|---|---|---|---|
| 0 (the body) | 354,735.33 | 535,623.19 | 1.51 | 2.584e-06 | 1.711e-06 |
| **3** | **181.01** | **359.33** | **1.99** | **1.381e-05** | **6.955e-06** |
| 4 (3's mirror twin) | 181.16 | 360.11 | 1.99 | 4.787e-07 | 2.408e-07 |
| 1, 2, 5, 6, 7, 8 | 16.8 – 42.9 | 61.7 – 97.8 | 2.26 – 3.67 | 1.1e-11 – 1.6e-09 | 2.9e-12 – 4.3e-10 |

Against that, the faces being merged carry recorded tolerances of **8.7e-04 to
1.5e-03 mm** (G12) — two orders larger. The two descriptions are the same
surface to well inside the kernel's own idea of one, which is also why the
boolean above finds nothing to report.

### Why no tighter bar is defensible

**Solid 3 and solid 4 are mirror images** — same volume and area to 0.1 %, both
99 → 68 faces, both 436 → 340 edges — and they drift **29x** apart. The
magnitude belongs to which merge the kernel happened to perform, not to the
geometry, so it cannot be predicted from the part, and a bar set just above
today's worst case is a bar waiting to refuse tomorrow's.

`UNIFY_VOLUME_TOL` is therefore **1e-5 → 1e-4**. Across the 1.5–3.7 per mm
surface-to-volume ratios this run spans, that admits at most ~3e-05 to 7e-05 mm
of boundary movement — still over an order of magnitude inside the face
tolerances above, so a merge across genuinely different surfaces cannot hide
under it. The guards that actually catch such a merge are unchanged and are the
ones docs/algorithm.md §9 always named as the stronger pair: the solid-count
check beside it, and `BRepCheck_Analyzer` immediately after.

### Both halves of the session's lesson

G19 and G20 are the same failure at one remove from each other: a scalar proxy
standing in for "did the geometry change", with its bar expressed in units that
do not match what it varies with. G19's answer was to replace the proxy with an
exact structural test, because one existed. G20's is to keep the proxy and size
it by the physical quantity it stands for, because no cheap exact test does — the
symmetric-difference boolean that settles it took seconds on a 99-face island and
does not finish quickly on the 69,305-face body beside it.

---

## G22 — one parallel analyzer instead of 301,505 serial ones ✅ PASS

`occ.fix_vertex_tolerances` (docs/algorithm.md §8) finds its work by asking, per
boundary face, whether `BRepCheck_Analyzer` rejects it. On the `cc=5, t=1`
rehearsal that scan is **44.6 s** over 301,505 faces — 5.8 % of `stitch`, ~1.2 %
of the run — serial on the master, and docs/specification.md §10 carried it as
the smaller of that chapter's two items along with two reasons the obvious
speed-ups do not work:

* dispatching it to workers has no free I/O route (~380 MB of writes for a
  ~37 s saving), and
* one analyzer over the assembled **solid** is a *different predicate* —
  `IsValid(subshape)` there is the in-context overload, and the standalone check
  is what this repair is calibrated against. That difference is exactly what
  G12 used to diagnose repair rung 2.

```
python tools/prototypes/g22_parallel_face_scan.py [extra.step ...]
```

**What this gate tests is neither of those.** A `TopoDS_Compound` of loose
faces has no shell and no solid in it, so a face's parent in the traversal is
the compound itself and there is no context to be checked in — while
`theIsParallel` (G18) still divides the scan across OCCT's own native threads.
If that reading is right, the predicate is unchanged and the scan needs neither
dispatch nor extra I/O.

### A, B — predicate identity, with the real faults in the corpus

Corpus: both committed golden samples (a real lattice output whose faces share
edges exactly as the sewn boundary layer's do) plus the four real invalid faces
committed from the rehearsal — 17,308 faces, of which 4 are genuinely invalid.

| Comparison against the per-face predicate | Disagreements |
|---|---|
| compound scan, `theIsParallel=False` | **0** of 17,308 |
| compound scan, `theIsParallel=True` | **0** of 17,308 |
| chunked at 1,000 (window boundaries between edge-sharing faces) | **0** of 17,308 |
| after `fix_vertex_tolerances` repairs all four | **0** of 4 |

All four committed faults read invalid on both sides, before repair and valid
on both sides after it (`repaired=4, still_invalid=0`). **That is the
load-bearing half**: per G10, a scanner that only ever agrees on sound geometry
proves nothing, and this scan decides which faces the repair ever looks at.

### C — what it buys

17,308 faces, 6-core machine:

| Scan | Wall | CPU | Core-equivalents | Speedup |
|---|---|---|---|---|
| serial, one analyzer per face | 4.067 s | 3.953 s | 0.97 | 1.00 |
| compound, `theIsParallel=False` | 4.238 s | 3.938 s | 0.93 | **0.96** |
| compound, `theIsParallel=True` | 2.454 s | 8.766 s | 3.57 | **1.66** |
| chunked at 1,000 / 5,000 / 20,000 | 2.450 / 2.443 / 2.458 s | | 3.35–3.54 | 1.65–1.66 |

The two compound rows are reported separately on purpose: making **one call
instead of N is worth nothing** (0.96×), so the whole 1.66× is OCCT's threads —
the same figure G18 measured for the whole-solid check (1.60×), which is what
one would expect if it is the same mechanism. Chunking costs nothing measurable
at any window size from 1,000 up.

### D — why it is chunked anyway

The analyzer holds its per-subshape results for as long as it is alive:
**~14 kB per face**, measured (RSS 333 → 589 MB over 17,308 faces). One call
over the rehearsal's 301,505 faces would hold **~4.2 GB** at a point in `stitch`
already carrying ~3 GB. `occ.FACE_SCAN_CHUNK` is therefore 20,000 faces
(~0.3 GB), sized by memory rather than by speed — which costs nothing, since
the speedup is flat from 1,000 upward.

**Adopted** into `occ.invalid_faces`, used by `fix_vertex_tolerances`.
Everything downstream of the scan — the repair itself and every re-check that
decides `repaired` against `still_invalid` — still runs one face at a time,
standalone, unchanged. `test_weld.py::test_the_batch_scan_is_the_same_predicate_as_the_per_face_one`
pins the result, fault fixtures included, as a permanent regression.

### E — what the corpus could not have caught, and what actually broke

**The change was wrong on first landing and the rehearsal said so within the
hour**: 19 faces corrected and **15 remaining invalid**, where the serial scan
had always reported 19 and none. The same run's validity gate then passed all
14 solids, which is what identified the 15 as phantoms.

**The first explanation fitted perfectly and was wrong.** Every fault in this
gate's corpus is a **loose** face — the four committed `.brep` fixtures — while
a sewn layer's faults have neighbours, and a compound analyzer holds one result
per subshape shared between the faces using it, so a face could be rejected in
batch for the fault beside it. A confirmation stage was added on that reading
and **the count stayed at 15**, which settles it: all 34 candidates are invalid
standalone, and the two predicates have still never been observed to differ.

**The real cause is ordering, in the caller.** Both repair rungs widen
tolerances on vertices and edges that neighbouring faces share, and widening is
monotonically permissive — so repairing one face can make the next one valid
before the loop reaches it. The serial loop asked at the moment it arrived at
each face and so never counted them; a scan-then-repair pass asks before any
repair has happened. `fix_vertex_tolerances` now re-checks each candidate as it
reaches it, which restores the original set exactly and cannot fail the other
way, since no repair here invalidates a face.

The confirmation stage stays, relabelled: **insurance on the case this corpus
cannot reach, not a measured necessity.** It costs tens of checks against
hundreds of thousands, and whatever the two stages might disagree about can only
run one way — the batch scan records at least the statuses the standalone one
does, so it can over-report and cannot hide a fault.

Two lessons. **A control whose faults are isolated cannot test a predicate whose
difference is about neighbours** — G10's, in a new place. And an explanation
that fits the symptom is not the mechanism: that is the fifth time in this
file's history (G9, G10, G11, G12), and what settled it here was the cheapest
available test — ship the fix the explanation implies and see whether the number
moves.

---

## G21 — can the seam-only round 2 be made correct on real trimmed parts? ❌ NO

docs/specification.md §10's one remaining live item, and the only lever it
listed as worth more than ~1 %. Boundary-sew round 2 sews only the
free-edge-bearing subset of each tile (`weld._split_seam_interior`), which G9
found wrong on the `cc=5, t=1` rehearsal; the `expected_rings` check catches it
and redoes the component with a full unsplit sew, costing **545 s of that run's
769 s `stitch`**. The proposal was to carry a seam face's *straddling*
neighbours — edges used twice in the tile but once in the subset — into the
sewn subset, so sewing cannot rebuild one while the carried face keeps the
original.

```
python tools/prototypes/g21_straddling_seam_split.py --cores 6
python tools/prototypes/g21_straddling_seam_split.py --tiles <dir>
```

Real geometry throughout: `TD_HX_rehearsal_test` at `cc=5, t=1`, classified and
trimmed for real (19,552 junctions → 21,955 pieces, 754 s), tiled into 33
tiles, round 1 run through worker `.brep` files exactly as production does.
Round-1 results are cached per tile so the analysis can be replayed without
re-trimming.

### The closure is the whole tile — measured, on five tiles

| Tile | Faces | \|S0\| (today) | \|S1\| | \|S2\| | Fixpoint | Hops |
|---|---|---|---|---|---|---|
| 0 | 23,546 | 7,555 | 12,491 | 19,954 | **23,523** | 8 |
| 1 | 22,872 | 6,514 | 10,966 | 17,769 | **22,837** | 8 |
| 2 | 19,201 | 6,275 | 10,303 | 16,339 | **19,189** | 6 |
| 3 | 16,831 | 4,116 | 6,983 | 11,635 | **16,723** | 18 |
| 4 | 16,172 | 4,581 | 7,740 | 12,663 | **16,160** | 7 |

**A subset with no straddling edge is a union of connected components of the
tile's face-adjacency graph**, and a tile's round-1 result is one connected
shell but for a few strays — so the fixpoint is the tile itself, less the 12–108
faces that are genuinely separate. Every hop short of it leaves a fresh
frontier of exactly the edges the closure exists to remove.

The cost table says the same thing in wall clock, on the full 302,576-face set:

| Round 2 variant | Time |
|---|---|
| full, unsplit | 520.89 s |
| seam-only (today) | 9.34 s |
| closure, 1 hop | 43.42 s |
| closure, 2 hops | 605.16 s |
| **closure to fixpoint** | **523.20 s** |

Closing to correctness costs what the fallback already costs, because it *is*
the fallback. **The proposal is disproved: do not build it.**

### What the production check actually reports, and the 10 edges hiding in it

`_sew_round_two` now records `(component, want, got_split, got_unsplit)` for
every repaired component — free, since the unsplit sew has to run before either
number is known. On the rehearsal:

    component 0: expected 73984 free edge(s), seam-only split gave 192692,
                 full unsplit sew gives 73994

Two separate findings in one line.

* **The split really is wrong** at production scale — 192,692 against 73,994 —
  so G9's mechanism is live, and the repair earns its cost on this part.
* **The correct sew misses the expectation too, by exactly 10.** `free_edges`
  counted **degenerate** edges, which are not holes: an edge with no extent is
  a parametric artefact whose one owning face uses it once by construction.
  `shell_defects` has always skipped them and records finding exactly 10 on
  this part (3.0e-9 to 8.3e-8 mm). So the check fired *whatever* the split did,
  and could not have reported a correct split as correct. Fixed by excluding
  degenerate edges in `free_edges`, which is where the two counts should have
  agreed all along; `test_weld.py` pins it with a sphere, whose poles the
  kernel gives degenerate edges of its own.

### The gate's own limitation, which nearly inverted its verdict

Part A of this gate — sewing the same tiles both ways and comparing — reported
the two **identical** (302,576 faces, 74,049 free edges, area to six decimals,
60,228 straddling edges and *none* rebuilt). That is not a contradiction of the
production numbers above; it is the gate reconstructing the pipeline's front
half and stopping one call short. `finalize_pieces` is what drops interface cap
faces and folds the kept ones back in, so without it the sewn face set has
holes where production has closed caps, and the defect does not arise on it.

**A prototype that rebuilds part of the pipeline has to rebuild all of it**, or
it measures a shape the pipeline never sees. The cross-check that caught it was
the production run's own evidence line — which is the general form of the
lesson G8 and G9 already left twice: the gate's scale was not the problem this
time, its *inputs* were.

## G23 — what are the 26 bad edges on the rehearsal's dominant body? ✅ MEASURED

`TD_HX_rehearsal_test.step` at `cc=5, t=1` is refused at `export truth`: solid
0, the 583,894-face dominant body, tessellates into 1,427,670 triangles carrying
26 edges not used by exactly two of them. Two exonerations were already on
record — the STEP write (the same 26 appear in memory) and `simplify` (the
pre-unification solid gives the same 26 in 1,509,507 triangles). So the defect
is already in the shell `weld.assemble` proved watertight.

**That is the tension this gate resolves.** `weld.shell_defects` proves every
*B-rep* edge is used exactly twice, keyed on `TopoDS_Shape::IsSame`.
`occ.exported_mesh_defects` counts *mesh* edges keyed on coordinates rounded to
1e-7 mm. `tools/prototypes/g23_bad_edge_provenance.py` records, per bad mesh
edge, which face and which B-rep edge produced it — via
`BRep_Tool.PolygonOnTriangulation_s`, which maps every boundary segment to its
edge exactly, with no tolerance and no search — and classifies it:

| kind | meaning |
|---|---|
| **T** | two coincident but *distinct* `TopoDS_Edge`s, each used by two faces — duplicate material |
| **M** | one *shared* edge whose two faces discretize it differently — no hole exists in the B-rep |
| **one-face** | a B-rep edge genuinely used by one face — a real hole, and a contradiction of `assemble` |
| **interior** | a segment no B-rep edge claims, i.e. inside one face's own triangulation |

### The result, and it refutes the documented reading

    1x/M 9    1x/interior 4    4x/T 2    4x/interior 10      one-face 0

**There are no holes.** docs/specification.md §10 recorded the 14 used-once
edges as holes, including "a 2.2 mm free edge is a real hole in a `t = 1` mm
lattice rather than a rounding artefact". Every one of those edges is `M`: the
B-rep edge carries **two** faces and only one of them contributed that segment.
The inference was reasonable and it was wrong, which is the sixth time in this
project's history (G9–G12, the pcurve chapter) that the first convincing reading
did not survive being measured.

**Only 2 of the 26 are a generator defect.** Both are `T`, at
[2055.604, −88.617, 968.389] and [2055.885, −88.396, 968.858], and they
implicate faces 121454 and 531940 — **both `Plane`, both area 7.368093e-01 mm²,
identical**. Those are the two coincident cap sheets left by `weld.unweldable`
declining the cap at [2056.81, −87.5, 967.548]: `InterfaceSet.decline` withdraws
*both* keys and `finalize_pieces` then keeps *both* sides' cap faces. The run
assembled **14** solids, not 15, so the two pieces stayed connected through
their other interfaces and the duplicate sheets ended up *inside* one solid —
which is why "declining costs an extra solid, never a hole" did not hold here.

### The controls, because an instrument that only agrees on sound geometry proves nothing (G10)

* **C1 negative** — the 13 sibling solids the same run measured clean come back
  at **0**, with triangle counts matching the production log exactly (348, 84,
  68 …).
* **C2 positive** — removing one face from a sound solid yields exactly its 4
  non-degenerate B-rep edges as 7 used-once mesh segments, **all classified
  `one-face`**. So a genuine hole is known to read as `one-face` — and solid 0
  has none.
* **Polygon coverage** — 2,516,915 edge/face pairs carry a polygon and **0** do
  not, so `interior` means "no B-rep edge claims this", not "the mapping failed".

### The faces are all valid, and the tolerances are tight

`BRepCheck_Analyzer` passes **all 18** implicated faces. Face tolerances are
1e-07 throughout and max edge tolerances run 1e-07 to 8.7e-04. **This is not a
fat-tolerance defect**, which rules out the family §8's two repair rungs address.

### The `interior` folds are a cone-apex degeneracy, and they reproduce alone

Meshed on its own, in a compound with no neighbour, **every implicated
`ConicalSurface` face reproduces its own over-used edge**; the planar and
cylindrical ones do not. All eight folding cones share one property: their `v`
range begins at exactly **−1.7321 = −√3**, the cone's apex, where the whole
u-range collapses to a single 3D point. The same three faces recur at
lattice-periodic positions (504823 / 558350 / 581848, identical to seven
figures), which is why the coordinates repeat at Δy = 45.0 and 95.0.

### Nine of the 25 were the gate's own arithmetic — fixed

The first reading of this was that the apex is where bad edges *sit*, and that
excluding edges incident to a pole node accounts for 9 of the 25. That is the
right count for the wrong reason, and the mechanism is sharper than the
correlation.

**A triangle with two vertices at one point contributes its collapsed key once
*and the real edge twice***, since that edge appears as both `(A,P)` and `(P,A)`
within the same triangle. With two sound neighbours also using it, a perfectly
closed fan reads as four uses. Mesh points are interned by rounded coordinate,
so at a pole — a cone apex, a sphere pole — a whole parametric range collapses
to one id and exactly that happens.

**The smallest statement of the bug is not the rehearsal at all:**

| shape | triangles | degenerate triangles | bad before | `by_use` before | after |
|---|---|---|---|---|---|
| `MakeSphere(10.0)` | 2022 | 2 | **4** | `{4: 2, 1: 2}` | **0** |
| `MakeCone(5, 0, 10)` | 445 | 1 | **2** | `{4: 1, 1: 1}` | **0** |

Two closed, sound solids the export-truth gate would have refused outright, with
one defect per degenerate triangle. On the rehearsal's dominant body the count
falls **25 → 16**, exactly **9** triangles are skipped, every pole-incident
reading disappears (`bad edges touching a pole: 0 of 16`) and **no `M` or `T`
edge does** — the rule does not reach into real defects.

**The blind-spot control is what chooses the rule.** Excluding edges *near* a
pole would hide a genuine hole whose boundary reached one; skipping triangles
that bound nothing cannot. Measured on a cone with its base disk removed — a
real hole on the boundary of the very face carrying the apex:

| | bad | `by_use` |
|---|---|---|
| before | 65 | `{1: 64, 4: 1}` |
| after | **63** | `{1: 63}` |

63 is exactly the base polyline. The fix makes the reading *cleaner*, not
weaker: the old rule piled two apex artefacts on top of a real hole.

The test is **interned-id equality rather than geometric area**, deliberately —
the defect is created by the interning, so it belongs at the interning, and a
thin triangle with three distinct ids still bounds material and is still
counted. It also subsumes the collapsed `lo == hi` key production used to count,
so production and this script now share one rule.

**This is the same counting bug this project has now found three times.**
`shell_defects` has always skipped degenerate edges; `free_edges` counted them
until G21 (10 phantom holes on this very part); `exported_mesh_defects` never
had the test. The controls are clean: the 13 sibling solids read 0 both ways,
`--control-c2` still reports a removed face's edges as `one-face`, and
`spiral-island-unwritable.brep` — the one body known to be genuinely unwritable
— meshes with **0** degenerate triangles, so all 10 of its defects are untouched.

**It does not make the part ship**, and that is the point to carry: 16 readings
remain — 9 `M`, 4 `interior`, 2 `T` and 1 `4x/interior` — on faces
`BRepCheck_Analyzer` calls valid.

### What this does and does not settle

Settled: there are no holes; the "holes" reading in §10 is withdrawn; the
declined cap really does leave duplicate material inside one solid; and the
faces carrying the rest are valid geometry that OCCT's own mesher discretizes
inconsistently.

Not settled: the 14 readings that are neither `T` nor pole-adjacent. They sit on
valid faces with tight tolerances, and `spiral-island-unwritable.brep` — a body
known to be genuinely broken — shows the **same** `M` and `interior` classes
(8 `M`, 2 `interior`, 0 `T`). So neither class can be dismissed as noise on the
strength of the rehearsal alone, and no narrowing of the gate is defensible on
two data points.

**Consequence for planning:** repairing the declined cap removes 2 of 26. The
part does not ship on the strength of that fix, and the remaining readings are
not something `boundary` or `weld` can repair.
