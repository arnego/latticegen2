# Implementation Guide: Fuse-Free Lattice Synthesis

Companion to [perf-rearchitecture-proposal.md](perf-rearchitecture-proposal.md).
This document is written so a less capable model (or a human unfamiliar with the
history) can implement the proposal **without re-deriving any analysis**. Follow it
in order. Where it says MUST, deviation is a bug. Where it gives a numeric pass bar,
measure before proceeding.

**Status: approved by the user on 2026-08-10**, with two modifications folded into
ground rules 2 and 6 below (the Julia implementation moves aside into `old-julia/` and
the new one takes its place; the CLI/exit-code/log/STEP-metadata surface is a
suggestion to improve on, not a contract to reproduce).

**Implemented on 2026-08-10.** This document is kept as the record of what was
planned; where the plan and the result diverge, the *code* and
[../algorithm.md](../algorithm.md) are authoritative. The three material
divergences, all driven by Phase-0 measurement
([../../tools/prototypes/RESULTS.md](../../tools/prototypes/RESULTS.md)):

1. **No parameter-window restriction and no legacy fallback path.** The `t < cc/2`
   limit §3's G1 was written to map does not exist — caps stay intact for the whole
   of `t < a` (algorithm.md §3.3).
2. **The interior join is a custom indexed build, not sewing.** §4 step 4 listed
   this as "v2, only if needed"; G2 showed sewing does not scale and glued booleans
   do not merge, so v2 was needed immediately and v1 was never shipped.
3. **G5 became moot.** Boundary junctions are attached by sewing at tolerance by
   design, so whether COMMON preserves untouched faces bit-exactly does not arise.

---

## 0. Ground rules (read first, apply always)

1. **Priorities are, in order: correctness, memory stability, speed**
   (specification.md "Key Considerations"). Every shortcut below is justified only
   because its failure mode is "do more work," never "wrong output." Preserve that
   property in anything you add.
2. **The new implementation is the product; the Julia pipeline is moved aside as
   reference.** (User decision, 2026-08-10 — this supersedes the earlier "do not touch
   the Julia pipeline, keep it on `main`" rule.) Concretely:
   - The Julia sources move **unmodified** into a new top-level `old-julia/` directory
     (`src/`, the Julia-only `tools/`, `Project.toml`/`Manifest.toml`, the
     `latticegen2.bat`/`.sh` wrappers, `tools/build/`). They are reference material,
     not a maintained parallel product — do not fix bugs or add features there.
   - The new implementation takes the front-and-center position in the repo: top-level
     `src/` (Python package `latticegen2`), top-level `tests/`, top-level `tools/`,
     and the root `README.md` documents it as *the* tool.
   - **Clean up everything that only existed to serve the old implementation** before
     committing: `Project.toml`/`Manifest.toml` at the repo root, the Julia wrapper
     scripts, `tools/build_app.jl` + `tools/build/` (PackageCompiler), and the
     `licenses/` entries that no longer apply (Julia, Gmsh, the Gmsh.jl binding,
     PackageCompiler). Keep `licenses/occt-LICENSE_LGPL_21.txt` (OCCT is still the
     kernel) and add the new dependencies' licenses. `licenses/libraries.md` must end
     up describing exactly the dependency set the new tool actually has.
   - **Ported tests are kept and used.** Every Julia test that covers behaviour the new
     implementation still has (lattice math, CLI validation, STEP header rewrite,
     classification, cleanup rules) is ported to pytest and must pass. Julia tests that
     only covered removed machinery (tiling, balanced fuse, the gmsh geometry kernel
     wrapper) are dropped with the machinery.
   - The old implementation's end of life is the `v1.0` tag on `main`. If it is ever
     needed again it is branched from there; `old-julia/` is a convenience copy for
     the transition, not the archive of record. Say so in `old-julia/README.md`.
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
6. **CLI, exit codes, log format, and STEP metadata are suggestions, not contracts.**
   (User decision, 2026-08-10 — this supersedes the earlier "these are contracts,
   implement them exactly" rule. There is no muscle memory or downstream automation
   built around them.) Treat specification.md §3/§5/§7 and docs/algorithm.md §8/§9 as
   the well-considered starting point they are, and improve on them where the new
   architecture makes something obsolete, misleading, or needlessly awkward. Two
   obligations come with that freedom:
   - **Surface every change.** Each deviation goes in a single "Changes from the Julia
     implementation" section of `README.md`, saying what changed and why. Do not let a
     change reach the user only as a surprise in a log file.
   - **Do not lose information.** specification.md §3's required end-of-run summary
     content (input parameters, start timestamp, duration, run characteristics, peak
     memory, output path) must still be reported, and every failure must still exit
     nonzero with one human-readable reason line. The *shape* is yours; the
     *information* is not optional.
   The e2e harness is expected to need adapting to whatever surface you settle on —
   adapt it deliberately, and never by weakening a check (rule 8's guardrail still
   applies: a gate may be *re-expressed*, never *loosened*, to make a test pass).
   Parameters that only existed to tune the old tiling/fusion machinery
   (`--tile-cells`, and `--workers`' tile-stage meaning) are prime candidates for
   removal or redefinition, since the stages they controlled no longer exist.
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

Standalone scripts under `tools/prototypes/`, each printing PASS/FAIL against its
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

Deliverable: `src/latticegen2/` generating a lattice for a **box** input, interior
path only, end-to-end to STEP.

1. **CLI** (`src/latticegen2/cli.py`): start from specification.md §3 — flags, ranges,
   the `t < a` cross-constraint, exit code 2 semantics, derived output/log naming —
   and simplify it per ground rule 6 where the old flags described machinery that no
   longer exists. Port `test/test_cli.jl`'s cases as pytest for whatever surface you
   settle on.
2. **Classification** (`src/latticegen2/classify.py`): port docs/algorithm.md §5 with one
   change — classify **half-struts** (segments `p → p + s*(a/2)*e(k)`), then derive
   node classes:
   - half-strut INTERIOR = segment-to-mesh distance > `r + d` AND midpoint inside
     (3-ray parity, fixed directions — copy the constants from `old-julia/src/classify.jl`);
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
3. **Template** (`src/latticegen2/junction.py`): §2 construction + G1 check at startup —
   G1 failure at the run's actual `(cc,t)` → exit 2 with a human-readable message
   naming the legacy fallback in `old-julia/`.
4. **Interior build** (`src/latticegen2/interior.py`):
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
5. **Export** (`src/latticegen2/stepout.py`): `STEPControl_Writer`, AP214, mm; then port
   `old-julia/src/stepmeta.jl`'s header rewrite (FILE_NAME first field = part name
   `<stem>+cc<cc>+t<t>`; append params to FILE_DESCRIPTION; fill FILE_SCHEMA **only
   if blank**). Round-trip re-read gate before exit 0.
6. **Phase-1 acceptance:** box `L×W×H`, `cc=10, t=1.5`: `BRepCheck_Analyzer` valid;
   and the exact volume identity holds: because half-struts partition the struts,
   all strut-strut overlap is internal to each `J`, and adjacent junctions meet only
   on shared faces (zero-volume contact), the union's volume is **exactly
   `N_instanced_nodes × volume(J)`** — measure `volume(J)` once with `GProp` and
   assert the assembled solid matches within `1e-6` relative. Also run
   `old-julia/tools/verify_geometry.jl`'s manifold + self-intersection checks on the output
   (they are stack-independent).

## 5. Phase 2 — boundary + connectivity

1. **Boundary workers** (`src/latticegen2/boundary.py`): for each BOUNDARY node, instance `J`,
   translate, one single-object `BRepAlgoAPI_Common` against the input body.
   Distribute across `multiprocessing` workers (respect `--workers`/`--cores`
   semantics, §3; `-bg` priority per docs/algorithm.md §7.3). Workers return
   serialized shapes or `.brep` paths + stats — mirror the existing small-IPC design
   (§6.2). Auto-tuning simplifies: no tile sizing needed (junction jobs are
   constant-size); keep the RSS watchdog with its bounded pause (§7.2).
2. **Attach:** for each surviving trimmed junction, glue to neighbors at intact cap
   quads (exact pairing if G5 passed; sewing fallback otherwise). A cap that was cut
   by the trim is simply exterior surface now — no action.
3. **Connectivity / floating rule** (`src/latticegen2/connect.py`): build the junction graph —
   vertices = instantiated junctions (interior + surviving boundary pieces); edge iff
   the shared cap interface exists on both sides post-trim. A trimmed junction that
   COMMON split into multiple solids contributes one graph vertex per piece; a piece
   connects through a cap only if that cap face belongs to that piece (face-membership
   check, no booleans). Then apply specification.md §5 exactly: BFS components; a
   component is dropped iff total volume < `t³` AND it shares no interface with the
   rest — **connectivity is now proof by construction, so the exit-4 "unresolvable"
   path of the old `filter_floating!` should be unreachable; keep it as a defensive
   assertion, and keep the single aggregate removal-log line format (§8).**
4. **Logging** (`src/latticegen2/runlog.py`): follow docs/algorithm.md §9's *intent* —
   always-on `.log`, `-v` console verbosity, per-stage lines, per-junction-batch stats,
   end-of-run summary carrying every spec §3 field, distinct nonzero exit codes with
   one human-readable reason line each, Ctrl+C graceful shutdown (§9.1: catchable
   interrupt, orderly worker stop with a short grace period, `CANCELLED` line, temp
   kept). Per ground rule 6 the exact codes and line formats are yours to set; retire
   codes whose failure mode no longer exists rather than reserving them for nothing,
   and record the final table in `README.md`.

## 6. Phase 3 — parity and verification

1. Port `tools/e2e.jl`'s scenarios (specification.md §6.1) to drive the new CLI, and
   port `old-julia/tools/verify_geometry.jl`'s checks (`manifold_check`,
   `triangles_properly_cross`/`self_intersection_check`, `golden_sample_volume_diff`)
   to Python alongside them — ground rule 2 retires the Julia toolchain, so keeping a
   Julia checker would mean keeping Julia, gmsh, and their licenses for one script.
   The checks stay *algorithmically* independent of the generator (they re-tessellate
   the finished STEP and reason only about the resulting triangles), and the
   independence lost by sharing a runtime is more than repaid by item 3, which the old
   stack could not do at all.
2. **Golden-sample gates:** `smoke-verified` vs
   `test/80mm-test-ball-cc20t4-golden-sample.step`; `dense-lattice` vs
   `test/test-cylinder-cc10t1.5-golden-sample.step`. Volume diff near zero both ways.
   Never adjust a tolerance to pass. **Note (user, 2026-08-10): both golden samples
   were produced by the old pipeline and are now considered suspect** — a mismatch is
   not automatically a bug in the new implementation. If a golden gate fails, stop and
   report the numbers to the user for manual verification rather than either
   "fixing" the generator to match or loosening the gate.
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

All of: Phase 0 bars recorded in `tools/prototypes/RESULTS.md`; Phases 1–3 tests
green; `BRepCheck` clean; `dense-lattice` < 10 min; golden-sample diffs either ≈ 0 or
escalated to the user per §6 item 2; repo restructured per ground rule 2 (`old-julia/`
populated, obsolete build/wrapper/license files removed, `licenses/libraries.md`
matching the real dependency set); ported tests green; README updated (deps, install,
usage, memory notes, algorithm overview pointing at the proposal doc, **and the
"Changes from the Julia implementation" section ground rule 6 requires**). Then hand
back to the user for the remaining decisions (merge strategy, the `v1.0` tag that ends
the Julia implementation's life, spec updates with proper TODO tags) — those are user
decisions, not implementer decisions.
