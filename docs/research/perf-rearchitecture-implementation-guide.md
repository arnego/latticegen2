# Implementation Guide: Fuse-Free Lattice Synthesis

Companion to [perf-rearchitecture-proposal.md](perf-rearchitecture-proposal.md).
This document is written so a less capable model (or a human unfamiliar with the
history) can implement the proposal **without re-deriving any analysis**. Follow it
in order. Where it says MUST, deviation is a bug. Where it gives a numeric pass bar,
measure before proceeding.

**Do not begin implementation until the user has approved the proposal** (the
project's specification rules require user sign-off on Claude-proposed features).

---

## 0. Ground rules (read first, apply always)

1. **Priorities are, in order: correctness, memory stability, speed**
   (specification.md "Key Considerations"). Every shortcut below is justified only
   because its failure mode is "do more work," never "wrong output." Preserve that
   property in anything you add.
2. **Do not modify the existing Julia pipeline.** It stays on `main` as the reference
   implementation. The new tool lives in a new top-level directory `pyfast/` on this
   research branch until golden-sample parity is proven.
3. **Never re-derive the lattice math.** Copy the formulas from
   [docs/algorithm.md](../algorithm.md) §2–§3 exactly (they are restated in §2 below
   with the additions this design needs). All angles/lengths from expressions
   (`asin(sqrt(2/3))`), never decimal literals.
4. **Never call a general (non-glued) boolean fuse on more than ~10 solids.** The
   entire point of this architecture is that large fuses are unnecessary. If you find
   yourself needing one, you have made an error — stop and re-read §4–§6 here.
5. **One object operand per COMMON call.** A COMMON with two overlapping object
   operands *partitions* them instead of trimming them (docs/algorithm.md §6.3, the
   operand-disjointness invariant). In this design every COMMON has exactly one
   object (a single already-fused junction instance), so the failure mode is
   unreachable — keep it that way.
6. **CLI, exit codes, log format, STEP metadata are contracts** — implement them
   exactly per specification.md §3/§5/§7 and docs/algorithm.md §8/§9. The e2e harness
   and the user's muscle memory both depend on them.
7. Target platform: Python 3.11, `cadquery-ocp` (OCP) for OCCT, NumPy. Vendor wheels
   for offline install; add license texts to `licenses/` and update
   `licenses/libraries.md` (OCP: Apache-2.0; NumPy: BSD-3; OCCT: LGPL-2.1 already
   present). Zero network at runtime.

---

## 1. Vocabulary

| Term | Meaning |
|---|---|
| `cc`, `t`, `a = cc/√2`, `r = t/√2`, `θ`, `e_k`, `B`, `node(i,j,k)` | Exactly as docs/algorithm.md §1–§2. |
| half-strut `(n, k, s)` | The half of a strut adjacent to node `n`, along direction `s·e_k`, `s ∈ {+1,−1}`, length `a/2`. |
| junction solid `J` | Union of the 6 half-struts owned by one node. Congruent at every node. |
| cap quad | The `t×t` diamond cross-section face at distance `a/2` from the node along `±e_k` — the interface where two adjacent junctions meet exactly. |
| junction graph | Graph with one vertex per instantiated node and one edge per surviving cap-quad interface. Connectivity queries run here, never in the geometry kernel. |

---

## 2. Exact geometry definitions

All from docs/algorithm.md §2–§3; the only new object is the half-strut.

```python
import numpy as np
theta = np.arcsin(np.sqrt(2/3))
def e(k):                       # unit strut direction, k = 0,1,2
    phi = 2*np.pi*k/3
    s, c = np.sqrt(2/3), 1/np.sqrt(3)
    return np.array([s*np.cos(phi), s*np.sin(phi), c])
# FACT (verified): e(0), e(1), e(2) are mutually orthogonal. Unit-test this.
a  = cc / np.sqrt(2)
r  = t / np.sqrt(2)
B  = a * np.column_stack([e(0), e(1), e(2)])
node = lambda i, j, k: B @ np.array([i, j, k])

def frame(k):                   # profile frame, docs/algorithm.md §3.1
    u = np.cross([0,0,1], e(k)); u /= np.linalg.norm(u)
    v = np.cross(e(k), u)
    return u, v

def profile_verts(center, k):   # diamond square, side t
    u, v = frame(k); h = t/np.sqrt(2)
    return [center + h*u, center + h*v, center - h*u, center - h*v]
```

**Half-strut solid** `(origin_node_pos p, k, s)`: the 4-vertex profile polygon at
`p` (same `frame(k)` for both signs — a strut's profile is constant along its axis),
extruded by `s * (a/2) * e(k)`. Build as: closed polygon wire → planar face →
`BRepPrimAPI_MakePrism` with vector `s*(a/2)*e(k)`.

**Junction template `J`:** at the origin (`p = (0,0,0)`), build all 6 half-struts and
fuse them with **one** `BRepAlgoAPI_Fuse`/`BOPAlgo_Builder` call (6 operands — this is
the ONLY general fuse in the entire program, executed exactly once per run,
milliseconds). Then run the **cap-integrity check** (§3, gate G1).

**Candidate index range:** port docs/algorithm.md §2.4 verbatim (`B \ corners`,
`pad = ceil(r/a) + 1`).

---

## 3. Phase 0 — de-risking prototypes

Standalone scripts under `pyfast/prototypes/`, each printing PASS/FAIL against its
bar. **Do not write product code until all five pass or their fallbacks are chosen.**

### G1 — junction template cap integrity
For every `(cc, t)` in the grid {cc ∈ 0.4, 5, 10, 20, 50} × {t ∈ 0.4, 1, 1.5, 4, 20,
and `t = 0.49·cc`, skipping invalid combos where `t ≥ a`}: build `J`, then verify each
of the 6 cap quads survives as a face of the fused solid — find a face of `J` whose 4
vertices match `profile_verts(p + s*(a/2)*e(k), k)` within `1e-9` and whose area is
`t²` within `1e-9` relative. Also verify `volume(J)` equals the inclusion-exclusion
prediction from independent pairwise/triple COMMONs (the §11.1 technique) to `1e-9`
relative.
**Pass bar:** all 6 caps intact for all combos with `t < cc/2`; record precisely
where (if anywhere) in `cc/2 ≤ t < a` caps start failing.
**Fallback if a needed combo fails:** for that parameter range only, the runtime MUST
refuse the fast interior path and fall back to the legacy Julia pipeline (print a
clear message); never silently produce geometry with broken interfaces.

### G2 — instancing + join throughput
Build `J` once; instance it at an `m×m×m` grid of nodes via `TopLoc_Location`
translation (no geometry copy); join with **both** candidate mechanisms and time
them at m³ ∈ {10³, 8·10³, 6.4·10⁴}:
  (a) `BRepBuilderAPI_Sewing` over all faces minus paired caps;
  (b) `BOPAlgo_Builder` with `SetGlue(BOPAlgo_GlueFull)` over the instanced solids.
**Pass bar:** either mechanism joins 6.4·10⁴ junctions (~1.6M faces) in < 5 min and
< 16 GB RSS. Record which is faster; the winner becomes the v1 join. (v2, only if
needed for the 64× target: direct shared-topology construction, §5 step 4.)

### G3 — per-junction COMMON latency
Load `test/test-cylinder.STEP`; take 200 junction instances straddling its surface;
run one single-object COMMON each against the body.
**Pass bar:** median < 50 ms, p95 < 250 ms. If input-body complexity ever breaks
this, implement the localization option (proposal §3) — not before.

### G4 — STEP writer throughput/memory
Write the G2 result via `STEPControl_Writer` (AP214, mm).
**Pass bar:** ≥ 3,000 faces/s and memory < 2× the model's in-core size. Also re-read
the file and confirm solid count and total volume round-trip.

### G5 — COMMON preserves untouched faces bit-exactly
For 50 G3 results whose interior-facing cap was classified strictly inside: check the
cap's 4 vertices are bit-identical (or ≤ 1e-12) to the pre-trim instance's.
**Pass bar:** all 50. **Fallback:** glue boundary junctions with sewing at tolerance
`1e-7` instead of exact index pairing (bounded O(surface) cost — acceptable).

---

## 4. Phase 1 — interior-only vertical slice

Deliverable: `pyfast/latticegen.py` generating a lattice for a **box** input, interior
path only, end-to-end to STEP.

1. **CLI** (`pyfast/cli.py`): port specification.md §3 exactly — flags, ranges,
   mutual-exclusion rules, exit code 2 semantics, derived output/log naming. Port
   `test/test_cli.jl` cases as pytest.
2. **Classification** (`pyfast/classify.py`): port docs/algorithm.md §5 with one
   change — classify **half-struts** (segments `p → p + s*(a/2)*e(k)`), then derive
   node classes:
   - half-strut INTERIOR = segment-to-mesh distance > `r + d` AND midpoint inside
     (3-ray parity, fixed directions — copy the constants from `src/classify.jl`);
   - half-strut OUTSIDE = distance > `r + d` AND midpoint outside;
   - else half-strut BOUNDARY.
   - **node INTERIOR** iff all 6 incident half-struts INTERIOR; **node OUTSIDE** iff
     all 6 OUTSIDE; else **node BOUNDARY**.
   Implement the spatial hash + segment-triangle distance with NumPy in batches
   (vectorize over struts, not triangles); port `check_surface_mesh_coverage`
   unchanged (both bars, §5.1 — it is a hard-won correctness gate, docs/algorithm.md
   §11.3). Unit-test against known sphere cases ported from `test/test_classify.jl`.
   For Phase 1's box input, an analytic inside-test may stub the mesh path, but the
   real mesh path must exist before Phase 2.
3. **Template** (`pyfast/junction.py`): §2 construction + G1 check at startup —
   G1 failure at the run's actual `(cc,t)` → exit 2 with a human-readable message
   naming the legacy fallback.
4. **Interior build** (`pyfast/interior.py`):
   - Extract `J`'s topology ONCE into an indexed structure: vertices (coords), edges
     (vertex-id pairs), faces (wire of edge ids + orientation), the 6 cap face ids,
     and for each cap the **precomputed correspondence map** cap`(+k)` vertex/edge ids
     ↔ cap`(−k)` ids of the neighbor template (match by coordinates once, at template
     build time, `1e-9`; store as index pairs — runtime does integer lookups only,
     never coordinate matching).
   - v1 (default): instance solids by `TopLoc_Location`, join with the G2 winner.
   - v2 (only if G2 shows v1 misses the 64× target): build one shell directly with
     `BRep_Builder`, reusing shared `TopoDS_Vertex`/`TopoDS_Edge` objects across
     neighboring instances via the correspondence map, omitting both members of every
     paired cap. Each mesh edge must be used exactly twice with opposite orientation —
     assert this.
5. **Export** (`pyfast/stepout.py`): `STEPControl_Writer`, AP214, mm; then port
   `src/stepmeta.jl`'s header rewrite verbatim (FILE_NAME first field = part name
   `<stem>+cc<cc>+t<t>`; append params to FILE_DESCRIPTION; fill FILE_SCHEMA **only
   if blank**). Round-trip re-read gate before exit 0.
6. **Phase-1 acceptance:** box `L×W×H`, `cc=10, t=1.5`: `BRepCheck_Analyzer` valid;
   and the exact volume identity holds: because half-struts partition the struts,
   all strut-strut overlap is internal to each `J`, and adjacent junctions meet only
   on shared faces (zero-volume contact), the union's volume is **exactly
   `N_instanced_nodes × volume(J)`** — measure `volume(J)` once with `GProp` and
   assert the assembled solid matches within `1e-6` relative. Also run
   `tools/verify_geometry.jl`'s manifold + self-intersection checks on the output
   (they are stack-independent).

## 5. Phase 2 — boundary + connectivity

1. **Boundary workers** (`pyfast/boundary.py`): for each BOUNDARY node, instance `J`,
   translate, one single-object `BRepAlgoAPI_Common` against the input body.
   Distribute across `multiprocessing` workers (respect `--workers`/`--cores`
   semantics, §3; `-bg` priority per docs/algorithm.md §7.3). Workers return
   serialized shapes or `.brep` paths + stats — mirror the existing small-IPC design
   (§6.2). Auto-tuning simplifies: no tile sizing needed (junction jobs are
   constant-size); keep the RSS watchdog with its bounded pause (§7.2).
2. **Attach:** for each surviving trimmed junction, glue to neighbors at intact cap
   quads (exact pairing if G5 passed; sewing fallback otherwise). A cap that was cut
   by the trim is simply exterior surface now — no action.
3. **Connectivity / floating rule** (`pyfast/connect.py`): build the junction graph —
   vertices = instantiated junctions (interior + surviving boundary pieces); edge iff
   the shared cap interface exists on both sides post-trim. A trimmed junction that
   COMMON split into multiple solids contributes one graph vertex per piece; a piece
   connects through a cap only if that cap face belongs to that piece (face-membership
   check, no booleans). Then apply specification.md §5 exactly: BFS components; a
   component is dropped iff total volume < `t³` AND it shares no interface with the
   rest — **connectivity is now proof by construction, so the exit-4 "unresolvable"
   path of the old `filter_floating!` should be unreachable; keep it as a defensive
   assertion, and keep the single aggregate removal-log line format (§8).**
4. **Logging** (`pyfast/runlog.py`): port docs/algorithm.md §9 — always-on `.log`,
   `-v` console verbosity, per-stage lines, per-junction-batch stats, end-of-run
   summary (all spec §3 fields), exit codes 0/2/3/4/5/6/130, Ctrl+C graceful shutdown
   semantics (§9.1: catchable interrupt, orderly worker stop with 2 s grace,
   `CANCELLED` line, temp kept).

## 6. Phase 3 — parity and verification

1. Port `tools/e2e.jl` scenarios (specification.md §6.1) to drive the new CLI;
   keep using the **Julia** `tools/verify_geometry.jl` + `golden_sample_volume_diff`
   as the independent checker (different stack than the generator — that independence
   is deliberate).
2. **Golden-sample gates:** `smoke-verified` vs
   `test/80mm-test-ball-cc20t4-golden-sample.step`; `dense-lattice` vs
   `test/test-cylinder-cc10t1.5-golden-sample.step`. Volume diff near zero both ways.
   Any mismatch is a stop-the-line correctness bug — never adjust tolerances to pass.
3. **New gate the old pipeline couldn't have:** `BRepCheck_Analyzer` on the final
   shape; also run it on the golden samples and record results (this finally
   resolves docs/algorithm.md §11.1's open self-intersection question with an exact
   check).
4. **Performance acceptance:** `dense-lattice` < 10 min end-to-end on the 6-core /
   32 GB dev machine; log the stage table for comparison with §1 of the proposal.

## 7. Phase 4 — scale rehearsal

Scale a synthetic input (e.g. the cylinder scaled 2× linearly) at `cc=5, t=1`:
measure classify / instancing / boundary / join / write / RSS. Compare against the
proposal §4 projections; if STEP write dominates, report file size + write time +
import feasibility to the user for the downstream-format decision — do not
unilaterally change the output contract.

---

## 8. Guardrails — common failure modes to avoid

- Do NOT fuse tiles, batches, or "just this group of 50" — see rule 4. The only
  legal kernel booleans are: the one 6-operand template fuse, and one single-object
  COMMON per boundary junction.
- Do NOT resolve connectivity questions with trial fuses. The junction graph answers
  them exactly and instantly.
- Do NOT match interface geometry by floating-point search at runtime; use the
  precomputed template index correspondence (or the sewing fallback wholesale).
- Do NOT weaken a verification gate to make a test pass (project history:
  docs/algorithm.md §11.1/§11.3 — gates were fixed, never loosened).
- Do NOT introduce per-item synchronize/flush patterns in loops (the 202× lesson,
  §11.2 finding #2). OCP has no gmsh model layer, but the same class of mistake
  exists anywhere an O(model) refresh hides inside an O(1)-looking call.
- Do NOT delete any solid the graph hasn't proven disconnected (spec §5).
- Any RAM estimate must respect the watchdog pattern: bounded pauses (§7.2), never
  unbounded waits on a monotonic high-water mark.
- Every temp artifact goes in `temp/<date><time>/` next to the output; keep on
  failure, delete on success (spec §4.4).

## 9. Definition of done

All of: Phase 0 bars recorded in `pyfast/prototypes/RESULTS.md`; Phases 1–3 tests
green; both golden-sample diffs ≈ 0; `BRepCheck` clean; `dense-lattice` < 10 min;
logging/exit-code parity spot-checked against a `main`-branch run; licenses folder
updated; README updated (deps, install, usage, memory notes, algorithm overview
pointing at the proposal doc). Then hand back to the user for the adoption decision
(merge strategy, deprecation of the Julia pipeline, spec updates with proper TODO
tags) — those are user decisions, not implementer decisions.
