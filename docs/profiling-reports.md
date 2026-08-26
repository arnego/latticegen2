# Profiling reports

Per-stage CPU, memory and I/O for whole runs of the **scale rehearsal**, kept
verbatim as [`tools/profile_report.py`](../tools/profile_report.py) printed
them. Wall-clock-only stage timings live in each run's own `.log`; this file is
for the resource picture that the log cannot give — how many cores a stage
actually used, how much memory it held, and how much it moved to and from disk.

Entries are in chronological order of the code they measure, oldest first, with
a note between consecutive entries saying what changed in between that is
relevant to performance. Add to the end rather than rewriting: the value of this
file is the series, and a figure is only comparable to its neighbours if the
method that produced it is stated.

## How these were produced

```
python tools/profile_run.py --out run.csv -- python src/main.py -i test/TD_HX_rehearsal_test.step -cc 5 -t 1 -o out.step --cores 6 -v
```

```
python tools/profile_report.py run.csv out.log --cores 6
```

`profile_run.py` samples the whole process tree — master plus every worker —
every 2 s; `profile_report.py` joins those samples to the stage boundaries in
the `.log`. **Cores used** is the number to read: 1.00 is one core fully busy,
6.00 is the whole machine. `profile_run.py` needs `psutil` for sampling the
process tree; that is a development-only dependency of `tools/`, installed
separately (`python -m pip install psutil`) — it is not a runtime dependency of
latticegen2 and does not ship in a bundle (docs/testing.md,
docs/specification.md §11).

Entries dated before 2026-08-17 quote a command line carrying `--ram`, which
was a real flag when they were run and has since been removed
(docs/specification.md §11). They are left as they were: they are records of
what was run, not instructions.

Machine for every entry below: Windows 11, 6-core CPU, 32 GB RAM,
Python 3.11.13, OCP 7.9.3.1.1 (OCCT 7.9.3), NumPy 2.4.6.

**Read a pair only when both halves were run on the same machine, close
together.** docs/specification.md §10 records a comparison across two sessions
whose *untouched* stages moved 25–36 % on ordinary machine load, which made its
touched-stage deltas unreadable. Where an entry below is one half of a
controlled pair, it says so.

---

## 2026-08-17 — `82adbb1` (pre-Phase-2 baseline)

Commit `82adbb1`, "Fix the last invalid boundary faces: falsely self-intersecting
wires (4 -> 0) (#18)". The first state of the tool in which this part completes
end to end: passes `validate` and writes its STEP unaided, with no bypass of the
`assemble` watertightness gate.

Run duration from the `.log`: **55 m 17.9 s**, peak memory 18.65 GB as the master
reports it. Output 2.01 GB, 584,028 faces, 2,517,881 edges, 14 solids, lattice
volume 330,354.002 mm³. Interior faces 705,000.

```
run started : 2026-08-17 03:36:40
samples     : 1620  (03:36:39 -> 04:31:57)
stages      : 12
cores       : 6 physical (100% CPU == 1 core)
```

| Stage | Duration | CPU mean | CPU peak | Cores used | RSS mean | RSS peak | Procs | Read | Written |
|---|---|---|---|---|---|---|---|---|---|
| template | 0.0s | 0% | 0% | 0.00 | 0 MB | 0 MB | 0 | 0 MB | 0 MB |
| import | 0.0s | 0% | 0% | 0.00 | 0 MB | 0 MB | 0 | 0 MB | 0 MB |
| tessellate | 3.0s | 103% | 108% | 1.03 | 266 MB | 299 MB | 1 | 0 MB | 0 MB |
| classify | 126.0s | 98% | 101% | 0.98 | 288 MB | 308 MB | 1 | 0 MB | 0 MB |
| boundary | 715.0s | 523% | 561% | 5.23 | 1,705 MB | 2,170 MB | 7 | 184 MB | 156 MB |
| connect | 11.0s | 99% | 100% | 0.99 | 2,235 MB | 2,241 MB | 7 | 0 MB | 0 MB |
| stitch | 629.0s | 109% | 550% | 1.09 | 3,039 MB | 3,139 MB | 7 | 382 MB | 365 MB |
| instance | 75.0s | 95% | 100% | 0.95 | 3,924 MB | 4,947 MB | 7 | 0 MB | 0 MB |
| assemble | 31.0s | 97% | 100% | 0.97 | 4,447 MB | 4,449 MB | 7 | 0 MB | 0 MB |
| simplify | 1,241.0s | 98% | 107% | 0.98 | 8,658 MB | 9,956 MB | 7 | 1,189 MB | 1,188 MB |
| validate | 243.0s | 98% | 102% | 0.98 | 10,089 MB | 15,394 MB | 7 | 469 MB | 469 MB |
| export | 229.0s | 98% | 101% | 0.98 | 11,779 MB | 19,827 MB | 7 | 2,017 MB | 2,017 MB |

```
total to last stage : 3,303.0s (55.0 min)
peak tree RSS       : 19,827 MB
```

---

## What changed between these two entries

Two changes to what the `simplify` stage has to do
(docs/specification.md §10, docs/algorithm.md §6 and §9). **Nothing else in the
pipeline was touched**, which is why the untouched stages below serve as the
control.

**The interior is now built pre-merged, and this is the change that matters
here.** `interior.py` used to emit four lateral faces per half-strut, so every
strut with two interior ends carried eight where four suffice, and `simplify`
spent the stage merging them back — rediscovering by geometric search over the
whole solid a pairing that is known by construction, one per surviving mid-strut
interface. The two half-faces are now spliced at template-build time into one
full-strut face, with the shared cap corners dropped because in the merged face
the edges meeting there both run along the strut axis. Fewer faces and edges are
therefore *created*, rather than created and then merged away. On this part that
takes interior faces from **705,000 to 389,492 (−44.8 %)**.

**Same-domain unification's two passes are now two calls** — faces first with
edge merging off, then edge merging alone over the result. This is measured
**neutral** for performance and is kept for structural reasons (it is what would
let a tiled unification reassemble by shared topology, docs/specification.md §10
Phase 3). Dropping the edge pass outright was tried and rejected: it made
`simplify` faster and handed all of it back to `validate` and `export`, which
scale with edge count too, for a 35 % larger file.

**What to expect in the numbers, therefore:** `instance` and `assemble` fall
roughly with the face count, since they build and collect those faces;
`simplify`, `validate` and `export` fall by less, since they are driven by the
*output*, which is unchanged by construction; `classify`, `boundary`, `connect`
and `stitch` should not move at all, and are the control.

---

## 2026-08-17 — `928cc57` (post-Phase-2)

Commit `928cc57`, "Build the interior pre-merged; split unification's two
passes". **The second half of a controlled pair with the entry above**: same
machine, same input and parameters, the two runs back to back. Note the
timestamps run the other way — this run executed first (02:44–03:36) and the
baseline immediately after it (03:36–04:31); nothing else was running on the
machine during either.

Run duration from the `.log`: **51 m 43.3 s**, peak memory 18.12 GB as the master
reports it. Output 2.00 GB, 584,028 faces, 2,517,881 edges, 14 solids, lattice
volume 330,354.002 mm³ — **identical to the baseline in every one of those
figures**, which is the check that this changed how the result is built and not
what it is. Interior faces 389,492.

```
run started : 2026-08-17 02:44:25
samples     : 1515  (02:44:26 -> 03:36:08)
stages      : 12
cores       : 6 physical (100% CPU == 1 core)
```

| Stage | Duration | CPU mean | CPU peak | Cores used | RSS mean | RSS peak | Procs | Read | Written |
|---|---|---|---|---|---|---|---|---|---|
| template | 0.0s | 0% | 0% | 0.00 | 0 MB | 0 MB | 0 | 0 MB | 0 MB |
| import | 0.0s | 0% | 0% | 0.00 | 0 MB | 0 MB | 0 | 0 MB | 0 MB |
| tessellate | 3.0s | 49% | 98% | 0.49 | 253 MB | 279 MB | 1 | 0 MB | 0 MB |
| classify | 126.0s | 98% | 102% | 0.98 | 288 MB | 306 MB | 1 | 0 MB | 0 MB |
| boundary | 721.0s | 520% | 562% | 5.20 | 1,704 MB | 2,132 MB | 7 | 182 MB | 156 MB |
| connect | 11.0s | 98% | 99% | 0.98 | 2,233 MB | 2,242 MB | 7 | 0 MB | 0 MB |
| stitch | 623.0s | 109% | 542% | 1.09 | 3,094 MB | 3,205 MB | 7 | 382 MB | 371 MB |
| instance | 42.0s | 99% | 101% | 0.99 | 3,556 MB | 4,046 MB | 7 | 0 MB | 0 MB |
| assemble | 21.0s | 98% | 101% | 0.98 | 3,822 MB | 3,822 MB | 7 | 0 MB | 0 MB |
| simplify | 1,088.0s | 98% | 108% | 0.98 | 6,884 MB | 7,339 MB | 7 | 941 MB | 940 MB |
| validate | 225.0s | 98% | 107% | 0.98 | 10,124 MB | 15,742 MB | 7 | 464 MB | 464 MB |
| export | 232.0s | 98% | 101% | 0.98 | 11,436 MB | 19,291 MB | 7 | 2,007 MB | 2,007 MB |

```
total to last stage : 3,092.0s (51.5 min)
peak tree RSS       : 19,291 MB
```

---

## Reading the pair

**The control held.** `classify`, `boundary`, `connect` and `stitch` — untouched
code, same input — agree to within 1 % across the pair (−0.3 %, +0.9 %, −0.6 %,
−1.0 %). That is what makes the rest of the deltas readable, and it is the
property the 2026-08-14/15 comparison lacked.

**Stage deltas**, from the `.log` timings rather than the rounded seconds above:
`instance` 1 m 14.2 s → 41.7 s (−43.8 %), `assemble` 30.9 s → 21.0 s (−31.9 %),
`simplify` 20 m 41.3 s → 18 m 07.9 s (−12.4 %), `validate` 4 m 02.8 s →
3 m 44.9 s (−7.4 %), `export` 3 m 49.3 s → 3 m 52.5 s (+1.4 %); total
55 m 17.9 s → 51 m 43.3 s (**−6.5 %**).

**Two thirds of the machine is idle for most of the run.** `boundary` is the only
stage that uses it (5.2 of 6 cores). Every other stage sits at 0.98 — one core.
Summing cores × duration gives roughly 6,340 CPU-seconds before and 6,140 after
against 3,303 s and 3,092 s of wall clock, i.e. **about 1.9–2.0 cores of 6, ~33 %
utilisation** (derived from the tables, not sampled directly). Phase 2 does not
change that and was not meant to — it removes work rather than spreading it, so
the CPU-seconds saved and the wall-clock saved match, as they must for stages
pinned to one core.

**`stitch`'s 109 % mean against a 542–550 % peak** is docs/specification.md §10's
open item visible in the profile: round 1 tiles across workers and genuinely
reaches ~5.5 cores, then round 2's unsplit redo runs serially and dominates the
mean.

**Memory climbs monotonically and peaks in `export`, not `simplify`:** 288 MB →
2.2 GB (connect) → 3.2 GB (stitch) → 4 GB (instance) → 7.3 GB (simplify) →
15.7 GB (validate) → 19.3 GB (export). The run's peak is the master holding the
whole 2 GB result while serialising it. So the headline "peak RSS −2.7 %"
understates the change: it barely touches the export peak that *sets* that
number, while cutting `simplify`'s own peak 9,956 → 7,339 MB (−26.3 %) and its
mean 8,658 → 6,884 MB (−20.5 %).

**I/O is almost entirely `.brep` worker round-trips plus the STEP write.**
`simplify`'s round-trip traffic falls 1,189 → 941 MB (−20.9 %), tracking the face
count; `export`'s ~2 GB written is the output file and is unchanged.

**One figure not to read anything into:** `tessellate` shows 1.03 cores before
and 0.49 after. It is a 3-second stage sampled every 2 s, so that is one or two
samples — noise, not a change.

**Where the remaining headroom is.** `simplify`, `validate` and `export` together
are 1,545 s, half the run, and all three are pinned at one core. `simplify` and
`validate` already dispatch across the shared pool *per body*;
docs/specification.md §10 records why that bought nothing on this part's very
unequal 14 solids. Sub-body tiling was aimed exactly at that gap and is now
**disproved** (docs/specification.md §11, G15) — as is the input-side
alternative tried in its place, which the entry below measures. This pair's own
finding is what predicted the second one: `simplify` scales with the output it
must produce (unchanged at 584,028 faces) more than with the input it consumes
(down ~31 %).

---

## 2026-08-17 — the restricted face merge (**reverted, kept as the measurement**)

Not a state of the tool. This run carries an uncommitted change that restricted
same-domain unification's *face* merge to the boundary layer plus one hop,
carrying the interior into the result by reference (docs/specification.md §11).
It is recorded because the number is the whole reason the change was reverted,
and a negative result nobody can see is a negative result somebody repeats.

**Read this against the entry above with more caution than that pair deserves.**
It was run the following morning rather than back to back, and the untouched
stages show it: `classify` agrees to **+0.1 %** and `connect` to +0.5 %, but
`boundary` is +4.0 % and `export` **+17.6 %** on identical code. So the
`simplify` figure below supports "no measurable win" and not a precise delta.

Run duration **54 m 10.9 s** against 51 m 43.3 s. Output **identical to both
entries above** in every figure — 584,028 faces, 2,517,881 edges, 14 solids,
330,354.002 mm³, 2.00 GB, volume drift 1.60e-07 — with 375,489 of 690,997 faces
reaching the kernel (−46 %) and 0 reassembly fallbacks. The restriction loses no
merge; it simply does not pay.

```
run started : 2026-08-17 11:41:30
samples     : 1588  (11:41:31 -> 12:35:41)
stages      : 12
cores       : 6 physical (100% CPU == 1 core)
```

| Stage | Duration | CPU mean | CPU peak | Cores used | RSS mean | RSS peak | Procs | Read | Written |
|---|---|---|---|---|---|---|---|---|---|
| template | 0.0s | 0% | 0% | 0.00 | 0 MB | 0 MB | 0 | 0 MB | 0 MB |
| import | 0.0s | 0% | 0% | 0.00 | 0 MB | 0 MB | 0 | 0 MB | 0 MB |
| tessellate | 3.0s | 47% | 95% | 0.47 | 247 MB | 268 MB | 1 | 0 MB | 0 MB |
| classify | 126.0s | 95% | 98% | 0.95 | 287 MB | 308 MB | 1 | 0 MB | 0 MB |
| boundary | 751.0s | 507% | 544% | 5.07 | 1,703 MB | 2,056 MB | 7 | 165 MB | 156 MB |
| connect | 11.0s | 98% | 100% | 0.98 | 2,231 MB | 2,254 MB | 7 | 0 MB | 0 MB |
| stitch | 636.0s | 109% | 526% | 1.09 | 3,071 MB | 3,178 MB | 7 | 382 MB | 360 MB |
| instance | 43.0s | 97% | 100% | 0.97 | 3,540 MB | 4,038 MB | 7 | 0 MB | 0 MB |
| assemble | 21.0s | 97% | 99% | 0.97 | 3,803 MB | 3,803 MB | 7 | 0 MB | 0 MB |
| simplify | 1,146.0s | 97% | 105% | 0.97 | 6,657 MB | 7,604 MB | 7 | 940 MB | 939 MB |
| validate | 228.0s | 97% | 104% | 0.97 | 9,811 MB | 15,128 MB | 7 | 464 MB | 464 MB |
| export | 274.0s | 96% | 100% | 0.96 | 11,675 MB | 19,266 MB | 7 | 2,007 MB | 2,007 MB |

```
total to last stage : 3,239.0s (54.0 min)
peak tree RSS       : 19,266 MB
```

**What this rules out, beyond the change itself.** `simplify`'s I/O is
**940 MB read / 939 MB written against the previous entry's 941 / 940** — so
the extra time is not the IPC of shipping the split, which was the obvious
suspect and is now excluded. Its RSS peak moved 7,339 → 7,604 MB (+3.6 %), the
bookkeeping's own footprint. Cores stayed at 0.97: this lever never touched
parallelism, and the stage is exactly as single-threaded as before.

The mechanism was then measured directly rather than inferred from this delta —
cutting the face merge's input 20 % cuts its time 6 %, an elasticity near 0.3,
where a *generic* subset of the same size gives 0.98 (G16) — which is what makes
the conclusion safe despite the imperfect pair, and what identifies the cause as
the restriction's own selection rather than the kernel. See
docs/specification.md §11.

---

## What changed between the previous entry and the next

Two stages were parallelised and one measurement gate closed a third direction.

* **`classify` dispatches across the shared pool** (docs/algorithm.md §5.4).
  Strided slices, an `.npz` mesh, pure NumPy — the one parallel stage that moves
  no geometry, so neither G7/G17's GIL result nor G15's identity result applies.
  The `WorkerPool` is consequently built *before* `classify` rather than after.
* **`validate` came off the pool and onto the master**, with
  `BRepCheck_Analyzer`'s own `theIsParallel` flag (docs/algorithm.md §9, G18).
  Native OCCT threads, so the GIL result does not bind them; and since the stage
  returns a scalar rather than geometry, G15 has nothing to attach to. This
  replaces specification.md §10's path 4, which cannot coexist with it under
  `--cores` (`W` processes × `W` threads).
* **`stitch` gained per-phase timers** — round 1, split, round 2, repair,
  retolerance, rings — because the stage's 0.96 mean against a 4.38 peak was
  visible but not attributable.

## 2026-08-17 — parallel `classify`, master-side `validate` (**partially contaminated**)

**Read the caveats before the table.** This is *not* a controlled pair, and two
separate things limit it:

1. **Another process started at 21:36:22**, three seconds before `stitch` ended.
   `template` through `stitch` are clean; `instance`, `assemble`, `simplify`,
   `validate` and `export` ran under competition and their wall times are
   inflated by an unknown amount.
2. **`boundary` landed in this project's known slow band.** It measured 889 s at
   **4.28 cores** against the previous entry's 721 s at 5.20 — but the
   2026-08-15 profile recorded 14 m 51 s at **4.22 cores** on untouched code,
   labelled "no — variance" in specification.md §10's own table. This run
   reproduces that point to within 1 %. So the whole-run total below is **not**
   comparable to the 3,092 s of the previous entry, and no conclusion should be
   drawn from it.

| Stage | Duration | CPU mean | CPU peak | Cores used | RSS mean | RSS peak | Procs | Read | Written |
|---|---|---|---|---|---|---|---|---|---|
| template | 0.0s | 0% | 0% | 0.00 | 0 MB | 0 MB | 0 | 0 MB | 0 MB |
| import | 0.0s | 0% | 0% | 0.00 | 0 MB | 0 MB | 0 | 0 MB | 0 MB |
| tessellate | 3.0s | 49% | 98% | 0.49 | 251 MB | 269 MB | 1 | 1 MB | 0 MB |
| **classify** | **47.0s** | **402%** | **492%** | **4.02** | 1,837 MB | 1,867 MB | 7 | 238 MB | 3 MB |
| boundary | 889.0s | 428% | 467% | 4.28 | 1,822 MB | 2,173 MB | 7 | 170 MB | 156 MB |
| connect | 13.0s | 94% | 98% | 0.94 | 2,332 MB | 2,347 MB | 7 | 0 MB | 0 MB |
| stitch | 770.0s | 96% | 438% | 0.96 | 3,203 MB | 3,361 MB | 7 | 382 MB | 353 MB |
| instance | 46.0s | 95% | 99% | 0.95 | 3,687 MB | 4,168 MB | 7 | 0 MB | 0 MB |
| assemble | 22.0s | 97% | 100% | 0.97 | 3,928 MB | 3,928 MB | 7 | 0 MB | 0 MB |
| simplify | 1,155.0s | 96% | 106% | 0.96 | 6,469 MB | 6,974 MB | 7 | 940 MB | 940 MB |
| **validate** | **114.0s** | **159%** | **464%** | **1.59** | 9,845 MB | 14,476 MB | 7 | **0 MB** | **0 MB** |
| export | 347.0s | 96% | 101% | 0.96 | 11,512 MB | 19,387 MB | 7 | 2,007 MB | 2,007 MB |

```
total to last stage : 3,406.0s (56.8 min)   <- not comparable, see caveats
peak tree RSS       : 19,387 MB
```

### The two readable results

**`classify`: 126 s → 47.0 s, 0.98 → 4.02 cores.** Measured in a clean window.
The classification is unchanged — 527,425 candidates → 29,375 interior, 19,552
boundary, exactly the previous entries' figures. The stage does not reach 6
cores because each worker rebuilds the mesh-derived indices (0.37 s measured
directly) and reads the staged `.npz` — visible as this stage's 238 MB read,
which did not exist before.

Note it undershoots the naive projection: the sweep alone measures 122.6 s
serial off the committed mesh, so 6 workers "should" give ~21 s. It gives 47 s.
The difference is staging, dispatch and per-worker index rebuild — the reason to
quote the stage rather than the kernel.

**`validate`: 225 s → 114.0 s, 0.98 → 1.59 cores, and I/O to exactly zero.**
This ran *under competition*, so the true figure is better than shown. The
`0 MB read / 0 MB written` is the most legible part of this entry: the previous
entry moved 464 MB each way to compute two scalars per solid, and running on the
master deleted that round trip outright. The 1.59 mean against a 4.64 peak is
the shape G18 predicts — OCCT's threads fully engaged on the dominant solid and
idle across the 13 scraps.

### What `stitch`'s new phase timers say, and what it costs

    round1 49.1s   split 3.0s   round2 15.1s   repair 651.2s
    retolerance 44.6s   rings 6.7s

**The repair is 85 % of the stage**, and this is the number that closes a
proposal rather than opening one. Speculatively running the unsplit sew
alongside the seam-only one — dispatching both and keeping whichever passes the
free-edge check, to take the discarded attempt off the critical path — can
recover at most the **15.1 s** the attempt costs. That is 0.5 % of the run, not
the minutes it was worth building for. The seam subset is small, so computing it
and throwing it away is nearly free.

What remains expensive is exactly what specification.md §10 already said: the
651 s full unsplit sew, reachable only by making the seam-only split correct in
the presence of straddling edges. That is an algorithmic fix, not a parallelism
lever, and no scheduling change touches it.

`retolerance` at 44.6 s is the other newly-visible figure — 5.8 % of `stitch`,
1.3 % of the run, and the one remaining embarrassingly-parallel item in this
stage.

### Correctness at production scale

Unaffected by contention, and the reason this run was still worth finishing:
**all 14 solids pass `BRepCheck_Analyzer`**, 19 vertex tolerances corrected with
no residual, 122,180 interfaces, 584,028 faces, 330,354.002 mm³, 2.00 GB — every
figure matching the previous entries.

One exception worth recording rather than smoothing over: **edges came back at
2,517,853 against 2,517,881, a difference of 28 (0.001 %)**. Nothing in this
change can plausibly cause it — the classification is identical, `validate` is
read-only, and `set_thread_budget(6)` is a no-op when OCCT's default on this
machine is already 6. It sits inside what docs/algorithm.md §9 describes as
same-domain unification's own representation choice. Flagged so a future run can
confirm it is variance rather than drift.

## What changed between the previous entry and the next

Two sessions' worth. `--ram` was removed (specification.md §11) and three guards
that refused valid input were fixed, neither of which touches per-stage cost.
Then the pair below, which closes specification.md §10:

* `occ.fix_vertex_tolerances` finds the faces it repairs with one parallel
  `BRepCheck_Analyzer` per 20,000-face window instead of one per face
  (docs/algorithm.md §8, G22). This is the only change expected to move a
  number, and it moves one *phase* of `stitch`, not a stage.
* `weld.free_edges` no longer counts degenerate edges, which corrects the
  round-2 check's own arithmetic. It changes what the log *says*, not what the
  run does on this part — the repair fires either way here.
* `_sew_round_two` records `(component, want, got_split, got_unsplit)` when it
  repairs, computed from a sew that had to run regardless.

## 2026-08-18 — `7e82e2a` (before), controlled pair

```
| Stage | Duration | CPU mean | CPU peak | Cores used | RSS mean | RSS peak | Procs | Read | Written |
|---|---|---|---|---|---|---|---|---|---|
| template | 0.0s | 0% | 0% | 0.00 | 0 MB | 0 MB | 0 | 0 MB | 0 MB |
| import | 1.0s | 0% | 0% | 0.00 | 233 MB | 233 MB | 1 | 0 MB | 0 MB |
| tessellate | 3.0s | 91% | 91% | 0.91 | 274 MB | 274 MB | 1 | 0 MB | 0 MB |
| classify | 42.0s | 443% | 533% | 4.43 | 1,790 MB | 1,853 MB | 7 | 333 MB | 4 MB |
| boundary | 786.0s | 486% | 543% | 4.86 | 1,808 MB | 2,255 MB | 7 | 181 MB | 156 MB |
| connect | 12.0s | 97% | 100% | 0.97 | 2,323 MB | 2,330 MB | 7 | 0 MB | 0 MB |
| stitch | 656.0s | 108% | 532% | 1.08 | 3,220 MB | 3,387 MB | 7 | 382 MB | 356 MB |
| instance | 45.0s | 97% | 99% | 0.97 | 3,727 MB | 4,192 MB | 7 | 0 MB | 0 MB |
| assemble | 23.0s | 96% | 100% | 0.96 | 3,961 MB | 3,961 MB | 7 | 0 MB | 0 MB |
| simplify | 1,129.0s | 97% | 103% | 0.97 | 6,611 MB | 7,472 MB | 7 | 940 MB | 940 MB |
| validate | 111.0s | 164% | 540% | 1.64 | 9,700 MB | 14,362 MB | 7 | 0 MB | 0 MB |
| export | 340.0s | 96% | 100% | 0.96 | 11,999 MB | 19,742 MB | 7 | 2,007 MB | 2,007 MB |

total to last stage : 3,148.0s (52.5 min)
peak tree RSS       : 19,742 MB
```

Run log: 52m 43.0s, 14 valid solids, 584,028 faces, 2,517,853 edges, 2.00 GB.
`stitch` phases: `round1 45.3s, split 2.8s, round2 14.7s, repair 542.0s,
retolerance 44.1s, rings 6.7s`.

## 2026-08-18 — the branch (after), same machine, 1 h 47 m later

```
| Stage | Duration | CPU mean | CPU peak | Cores used | RSS mean | RSS peak | Procs | Read | Written |
|---|---|---|---|---|---|---|---|---|---|
| template | 0.0s | 0% | 0% | 0.00 | 0 MB | 0 MB | 0 | 0 MB | 0 MB |
| import | 0.0s | 0% | 0% | 0.00 | 0 MB | 0 MB | 0 | 0 MB | 0 MB |
| tessellate | 3.0s | 48% | 95% | 0.48 | 250 MB | 267 MB | 1 | 1 MB | 0 MB |
| classify | 43.0s | 438% | 526% | 4.38 | 1,806 MB | 1,850 MB | 7 | 282 MB | 4 MB |
| boundary | 798.0s | 480% | 542% | 4.80 | 1,806 MB | 2,122 MB | 7 | 148 MB | 156 MB |
| connect | 13.0s | 96% | 98% | 0.96 | 2,320 MB | 2,336 MB | 7 | 0 MB | 0 MB |
| stitch | 652.0s | 112% | 527% | 1.12 | 3,233 MB | 3,813 MB | 7 | 382 MB | 352 MB |
| instance | 46.0s | 96% | 101% | 0.96 | 3,775 MB | 4,190 MB | 7 | 0 MB | 0 MB |
| assemble | 22.0s | 97% | 99% | 0.97 | 3,951 MB | 3,954 MB | 7 | 0 MB | 0 MB |
| simplify | 1,148.0s | 97% | 107% | 0.97 | 7,022 MB | 7,485 MB | 7 | 940 MB | 940 MB |
| validate | 112.0s | 160% | 524% | 1.60 | 10,358 MB | 15,078 MB | 7 | 0 MB | 0 MB |
| export | 339.0s | 95% | 101% | 0.95 | 11,720 MB | 19,050 MB | 7 | 2,007 MB | 4,056 MB |

total to last stage : 3,176.0s (52.9 min)
peak tree RSS       : 19,050 MB
```

Run log: 53m 11.3s, 14 valid solids, and **the same output byte for byte** —
584,028 faces, 2,517,853 edges, 2,148,818,507 bytes, identical to the before
run's file size and every reported figure. `stitch` phases: `round1 44.3s,
split 2.8s, round2 14.3s, repair 559.1s, retolerance 22.6s, rings 8.6s`.

## Reading the pair

**These two tables have a downstream consumer.** The stage weights the
graphical front-end's progress bar uses
([`src/latticegen2/gui/weights.py`](../src/latticegen2/gui/weights.py)) are the
per-mille shares of the mean of these two columns — this pair specifically,
because its untouched stages agree to within 1.5 %, which is what makes
averaging them legitimate. If this part is ever re-profiled, that table is
what needs updating with it; a test pins that its keys match
`runlog.STAGES` and that it sums to 1000, but nothing can tell it the numbers
have gone stale.

**Read the phase, not the stage.** The only thing this pair changes is one
phase of `stitch`:

| | before | after | |
|---|---|---|---|
| `retolerance` | 44.1 s | **22.6 s** | **−48.8 %** |
| `stitch` | 656 s | 652 s | −0.6 % |
| whole run | 52m 43s | 53m 11s | +0.9 % |

The stage moves by less than its own repair phase varies between runs (542.0 s
against 559.1 s, +3.2 %, on identical code paths and identical input), and the
whole run is inside ordinary variance. **That is the honest result**: 21.5 s
off a 3,175-second run is 0.7 %, real but invisible at stage granularity, and
the only reason it is measurable at all is that `stitch` reports its phases
separately. A pair like this one is worth running anyway, and this one paid for
itself twice over — see below.

The untouched stages agree to within 1.5 % (`boundary` +1.5 %, `simplify`
+1.7 %, `classify` +2.4 %, `validate` +0.9 %, `export` −0.3 %), which is what
makes the phase delta readable. `peak tree RSS` falls 3.5 % and
`min system available` rises, both incidental.

**A third run sits between these two and is not reported as a pair half.** The
branch's first version scanned every face up front and then repaired, and its
run reported 19 faces corrected with **15 "still invalid"** where the serial
scan had always reported none — while its validity gate passed all 14 solids,
which is what identified them as phantoms. The cause was ordering rather than
the predicate (specification.md §11): repairs widen shared tolerances, so a
neighbour can be fixed for free before the loop reaches it. The run above is
after that fix, and reports the original 19 with no residual.

**Two things this pair establishes that the timings do not.** The evidence line
`component 0: expected 73984 free edge(s), seam-only split gave 192682, full
unsplit sew gives 73984` says the seam-only split is genuinely wrong on this
part by a factor of 2.6 — and that the unsplit sew now meets the expectation
exactly, where before the degenerate-edge fix it read 73,994 and could not.

## What changed between the previous entry and the next

Two things, and only the first is a performance change at all.

**`validate` now carries the export-truth gate on every solid, with no size
bound** (docs/algorithm.md §9). The 2026-08-18 pair predates that gate
entirely, so the `validate` row below is not comparable with the one above it,
and neither is the whole-run total.

**The gate gained a second pass** (`occ.refine_until_manifold`): a body that
produces readings at 0.05 mm is re-measured at finer deflections over a
neighbourhood of the faces carrying them, and clears only on an exact zero. It
runs on the failing path alone. Here it runs once, on solid 0, over a 53-face
neighbourhood, and costs nothing measurable against the 2,003 s the gate's first
pass spends writing, re-reading and tessellating a 583,892-face body.

Also in this branch, and visible in the run log rather than the table:
same-domain unification now keeps the un-unified solid when its result is
invalid, and the pinhole repair's per-face area test compares an integral with a
relative bar instead of `!=` (docs/specification.md §11).

## 2026-08-26 — the branch, `TD_HX_rehearsal_test` at `cc=5, t=1`, single run

**This is one run, not a controlled pair.** Every performance claim in this file
that matters is made from a pair run back to back on the same machine; this is
not one, and it is not offered as a performance measurement. What it is for is
pricing the part end to end now that it *completes*, and supplying the shares in
`src/latticegen2/gui/weights.py`.

```
run started : 2026-08-26 09:45:13
samples     : 2683  (09:45:14 -> 11:16:24)
stages      : 12
cores       : 6 physical (100% CPU == 1 core)

| Stage | Duration | CPU mean | CPU peak | Cores used | RSS mean | RSS peak | Procs | Read | Written |
|---|---|---|---|---|---|---|---|---|---|
| template | 0.0s | 0% | 0% | 0.00 | 0 MB | 0 MB | 0 | 0 MB | 0 MB |
| import | 1.0s | 0% | 0% | 0.00 | 233 MB | 233 MB | 1 | 0 MB | 0 MB |
| tessellate | 3.0s | 114% | 114% | 1.14 | 272 MB | 272 MB | 1 | 0 MB | 0 MB |
| classify | 30.0s | 532% | 589% | 5.32 | 1,804 MB | 1,823 MB | 7 | 294 MB | 4 MB |
| boundary | 1,158.0s | 557% | 592% | 5.57 | 2,069 MB | 2,410 MB | 7 | 154 MB | 158 MB |
| connect | 33.0s | 100% | 101% | 1.00 | 2,592 MB | 2,621 MB | 7 | 0 MB | 0 MB |
| stitch | 588.0s | 116% | 584% | 1.16 | 3,727 MB | 4,318 MB | 7 | 383 MB | 358 MB |
| instance | 41.0s | 102% | 151% | 1.02 | 4,246 MB | 4,654 MB | 7 | 0 MB | 0 MB |
| assemble | 20.0s | 99% | 100% | 0.99 | 4,430 MB | 4,430 MB | 7 | 0 MB | 0 MB |
| simplify | 1,115.0s | 106% | 582% | 1.06 | 7,970 MB | 16,006 MB | 7 | 941 MB | 940 MB |
| validate | 2,105.0s | 109% | 586% | 1.09 | 14,910 MB | 21,123 MB | 7 | 2,007 MB | 2,007 MB |
| export | 362.0s | 99% | 101% | 0.99 | 13,260 MB | 21,182 MB | 7 | 2,007 MB | 4,056 MB |

total to last stage : 5,456.0s (90.9 min)
peak tree RSS       : 21,182 MB
min system available: 335 MB
total written       : 7,563 MB
total read          : 5,886 MB
```

**The part writes its STEP for the first time.** Run log: 1 h 31 m 10.5 s,
2,148,943,726 bytes, 14 solids, 584,114 faces, 2,518,001 edges, lattice volume
330,346.9858 mm³, peak 18.65 GB. `export_truth_s: 2003.07` — **95 % of the
`validate` stage**, and 37 % of the whole run.

The line this entry exists for:

```
solid 0: 583892 faces -> 1427664 triangle(s) after a STEP round trip,
         13 non-manifold edge(s) (9 degenerate triangle(s) skipped)
  by use count {1: 13}; positions [...]
  on 11 face(s); re-measured over a 53-face neighbourhood at
  0.01:8 0.002:4 0.0005:4 0.0001:0 -- resolved
```

**Those four numbers are G24's whole-solid sweep, reproduced by the cheap
route.** G24 measured the entire 583,892-face body at the same four deflections
and got 8, 4, 4, 0. The neighbourhood is derived automatically here — 11 core
faces rather than the 3 G24 picked by hand, 53 faces in the extract rather than
68 — and it returns the same counts at every rung. That is the oracle agreeing
on the *counts*, not merely on the verdict, on the one production body this has
ever been asked about.

**`validate` at 1.09 cores is not a parallelism opportunity**, and the
"Parallelization candidates" section at the foot of the report is wrong about it
for a reason worth stating: the stage is one `write -> read -> tessellate` per
solid and the run's 14 solids are one dominant body plus 13 scraps, so the floor
is that body however the work is spread. It is the same shape of finding as
`simplify`'s, which docs/specification.md §11 records at length.

**Reading this against the 2026-08-18 pair is a mistake**, and the numbers are
far enough apart to invite it. That pair totals 3,176 s against this run's
5,456 s, but 2,003 s of the difference is a correctness gate that did not exist
then, and the rest is one run's worth of machine variance against a pair's. The
untouched stages here (`boundary` 1,158 s against 798 s, `stitch` 588 s against
652 s) move in *both* directions, which is the signature of variance rather than
of a change.

**One caveat on `simplify`'s 1,115 s.** This run predates the code-review fix
that stops the worker-side validity check asking OCCT for a thread pool: at the
time, six workers each launched six threads on six cores. The verdict is
unaffected (G18 measured `BRepCheck_Analyzer` returning the same answer either
way), so the geometry and every other figure here stand — but that stage's time
was measured under an over-subscription that no longer happens, and would need
re-measuring before it is quoted as a cost. It is 204 of 1000 in the bar, so
being somewhat wrong about it moves the bar by less than the single-run basis
already does.

**Downstream consumer.** `src/latticegen2/gui/weights.py` takes its per-mille
shares from this table's Duration column. The previous shares gave `validate`
**35** of 1000; it is **385** here, so the bar would otherwise have sat still for
thirty-five minutes and then jumped. If the rehearsal is re-profiled again, that
table has to be recomputed with it.
