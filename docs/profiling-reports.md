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
python tools/profile_run.py --out run.csv -- python src/main.py -i test/TD_HX_rehearsal_test.step -cc 5 -t 1 -o out.step --cores 6 --ram 20 -v
```

```
python tools/profile_report.py run.csv out.log --cores 6
```

`profile_run.py` samples the whole process tree — master plus every worker —
every 2 s; `profile_report.py` joins those samples to the stage boundaries in
the `.log`. **Cores used** is the number to read: 1.00 is one core fully busy,
6.00 is the whole machine. `profile_run.py` needs `psutil` for sampling the
process tree; `psutil` is a runtime dependency of the tool itself, so it is
already present by any of README.md's routes (docs/testing.md).

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
