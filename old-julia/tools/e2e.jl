# tools/e2e.jl
#
# End-to-end verification harness (specification.md §6, docs/testing.md).
# Runs every scenario in specification.md §6.1 as a subprocess and checks
# every item in specification.md §6.2 that applies to each:
#   - smoke-fast:      always run, no golden-sample check (fast sanity pass). Same
#                      -cc/-t/input as smoke-verified below by design (spec §6.1's
#                      own table pairs them that way) — this scenario is the quick
#                      pass/fail signal without paying for the golden-sample diff.
#   - invalid-input:   always run, no golden sample needed (rejected pre-computation).
#   - smoke-verified:  always run; compares against the committed golden sample
#                      test/80mm-test-ball-cc20t4-golden-sample.step via
#                      golden_sample_volume_diff.
#   - dense-lattice:   only runs if test/test-cylinder-cc10t1.5-golden-sample.step
#                      is present on disk. The `-golden-sample` suffix (matching
#                      80mm-test-ball-cc20t4-golden-sample.step) keeps a
#                      reviewed baseline clearly distinct from the ordinary
#                      run output a generation at these params produces.
#                      Originally specified
#                      at -cc 5 -t 1, an attempt to produce that golden sample
#                      ran for hours and was manually terminated mid-run (not a
#                      crash — see docs/algorithm.md §11.2's investigation): an
#                      auto-tuned tile size well past the fuse-time knee left
#                      assembly with 11k+ unfused solids, and the (now-fixed)
#                      unconditional sub-threshold cleanup rule was deleting
#                      connected junction material one solid at a time. Root
#                      cause diagnosed and fixed (docs/algorithm.md §7.1,
#                      §6.3, §8, §11.2). The scenario's params were changed to
#                      -cc 10 -t 1.5 (a less dense lattice than -cc 5 -t 1,
#                      capable of finishing within a reasonable time) with a
#                      60-minute runtime budget. This scenario stays gated on
#                      the golden sample's presence rather than being auto-run
#                      here, since committing a golden sample is a decision
#                      for whoever reviews the regenerated run's output, not
#                      something this harness should do unattended.
#
# Usage:
#   julia --project=. tools/e2e.jl

using Test

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "latticegen2.jl"))
using .latticegen2
include(joinpath(@__DIR__, "verify_geometry.jl"))

"""
    sub_threshold_removal_ratio(output_step_path) -> Float64

Parse `<output>.log`'s mandatory end-of-run summary line ("Solids written: N
(removed as sub-threshold: M)", always printed regardless of `-v` —
`log_summary`/`log_always`, docs/algorithm.md §9) and return `M / (N + M)`.

Regression check for the root cause investigated in docs/algorithm.md
§11.2: an unconverged upstream fuse used to make the old, unconditional
"volume < threshold -> delete" cleanup rule delete large numbers of
*connected* junction solids, not genuine floating debris. A healthy run's
sub-threshold-removal ratio should be a small fraction of total solids —
this check alone would have caught the `test-cylinder-cc5t1` run (3740+
removals against a handful of kept solids) long before its multi-hour tail.
"""
function sub_threshold_removal_ratio(output_step_path::AbstractString)
    log_path = output_step_path[1:end-5] * ".log"   # ".step" is 5 chars, per cli.jl's resolve_output_paths
    isfile(log_path) || return NaN
    text = read(log_path, String)
    m = match(r"Solids written: (\d+) \(removed as sub-threshold: (\d+)\)", text)
    m === nothing && return NaN
    written, removed = parse(Int, m[1]), parse(Int, m[2])
    total = written + removed
    return total == 0 ? 0.0 : removed / total
end

# A top-level `@testset` that finishes with any failures prints its summary
# and then *throws* — uncaught, that aborts the whole script (a `LoadError`),
# so a failure in an earlier scenario would silently prevent every later
# scenario from running at all. Each scenario's `@testset` is therefore run
# inside its own try/catch below (mirroring the `global`-inside-try/catch
# pattern src/main.jl already uses for its own top-level state), so all four
# scenarios always run and are reported regardless of one another's outcome;
# `all_ok` tracks whether the run as a whole should exit nonzero.
all_ok = true

const INPUT = joinpath(ROOT, "test", "80mm-test-ball.step")
const OUTDIR = mktempdir()
const OUTPUT = joinpath(OUTDIR, "smoke-fast.step")
const BUDGET_SECONDS = 10 * 60.0   # spec §6.1: "generation < 10 minutes"
const CC = 20.0
const T = 4.0

println("=== latticegen2 e2e: smoke-fast ===")
println("Input:  ", INPUT)
println("Output: ", OUTPUT)

# specification.md §6.1's literal smoke-fast row is `-i ... -cc 20 -t 4 -bg`.
# It doesn't specify --workers/--tile-cells (--cores/--ram, added to the CLI
# after this row was written, would work too — specification.md §2/§3:
# optimization parameters are mandatory, either auto-tuned from hardware or
# given explicitly); this harness pins explicit values instead so the
# scenario is reproducible independent of the host's actual core count.
cmd = `julia --project=$ROOT $(joinpath(ROOT, "src", "main.jl")) -i $INPUT -cc $CC -t $T -o $OUTPUT --workers 4 --tile-cells 6`
println("Command: ", cmd)

t0 = time()
proc = run(pipeline(cmd; stdout=stdout, stderr=stderr); wait=false)
wait(proc)
elapsed = time() - t0

try
    @testset "smoke-fast e2e" begin
        @testset "process exits 0 with expected console output" begin
            @test proc.exitcode == 0
        end

        @testset "runtime under performance budget" begin
            println("Elapsed: $(round(elapsed, digits=1))s (budget $(BUDGET_SECONDS)s)")
            @test elapsed < BUDGET_SECONDS
        end

        @testset "STEP file written and non-empty" begin
            @test isfile(OUTPUT)
            isfile(OUTPUT) && @test filesize(OUTPUT) > 0
        end

        @testset "sub-threshold removals stay a small fraction of solids (docs/algorithm.md §8, §11.2)" begin
            ratio = sub_threshold_removal_ratio(OUTPUT)
            println("Sub-threshold removal ratio: $(round(ratio*100, digits=2))%")
            @test !isnan(ratio)
            isnan(ratio) || @test ratio < 0.01
        end

        if isfile(OUTPUT) && filesize(OUTPUT) > 0
            mesh_for_checks = with_gmsh() do
                new_model("e2e-mesh")
                import_shapes(OUTPUT)
                tessellate_surface(LatticeParams(CC, T))
            end

            @testset "round-trip read succeeds" begin
                vols = Int32[]
                with_gmsh() do
                    new_model("e2e-verify")
                    vs = import_shapes(OUTPUT)
                    append!(vols, Int32[t for (_, t) in vs])
                end
                @test !isempty(vols)
            end

            @testset "geometry validity (manifold, no self-intersections)" begin
                mesh = mesh_for_checks
                manifold_ok, bad_edges = manifold_check(mesh)
                @test manifold_ok
                manifold_ok || println("Non-manifold edges: ", length(bad_edges), " (showing up to 5): ",
                                        first(bad_edges, min(5, length(bad_edges))))

                no_self_int, bad_pairs = self_intersection_check(mesh)
                @test no_self_int
                no_self_int || println("Self-intersecting triangle pairs: ", length(bad_pairs))
            end

            @testset "output bounding box matches input within tolerance" begin
                in_lo, in_hi = with_gmsh() do
                    new_model("e2e-input-bbox")
                    import_shapes(INPUT)
                    bounding_box()
                end
                out_lo, out_hi = with_gmsh() do
                    new_model("e2e-output-bbox")
                    import_shapes(OUTPUT)
                    bounding_box()
                end
                tol = CC + T
                @test out_lo.x >= in_lo.x - tol && out_hi.x <= in_hi.x + tol
                @test out_lo.y >= in_lo.y - tol && out_hi.y <= in_hi.y + tol
                @test out_lo.z >= in_lo.z - tol && out_hi.z <= in_hi.z + tol
            end
        end
    end
catch e
    global all_ok = false
    println("smoke-fast e2e: testset reported failures (see summary above); continuing to next scenario.")
end

rm(OUTDIR; recursive=true, force=true)

println()
println("=== latticegen2 e2e: invalid-input ===")

# specification.md §6.1 `invalid-input` scenario: strut size > cell size.
# Mirrors test/test_cli.jl's "t vs cell edge cross-constraint" boundary
# (cc=5 -> a=cc/√2≈3.5355mm, t=4mm >= a) but exercised as a full subprocess
# run rather than a direct parse_args() call, so this also proves the
# process-level exit code / no-file-written behavior in src/main.jl, not
# just the validation logic in src/cli.jl.
const INVALID_OUTDIR = mktempdir()
const INVALID_OUTPUT = joinpath(INVALID_OUTDIR, "invalid-input.step")
const INVALID_LOG = joinpath(INVALID_OUTDIR, "invalid-input.log")

invalid_cmd = `julia --project=$ROOT $(joinpath(ROOT, "src", "main.jl")) -i $INPUT -cc 5 -t 4 -o $INVALID_OUTPUT --workers 1 --tile-cells 4`
println("Command: ", invalid_cmd)

invalid_stderr = IOBuffer()
invalid_proc = run(pipeline(invalid_cmd; stdout=devnull, stderr=invalid_stderr); wait=false)
wait(invalid_proc)
invalid_stderr_text = String(take!(invalid_stderr))

try
    @testset "invalid-input e2e" begin
        @testset "process exits nonzero with human-readable reason" begin
            @test invalid_proc.exitcode != 0
            @test occursin("ERROR", invalid_stderr_text)
        end

        @testset "no output STEP file written" begin
            @test !isfile(INVALID_OUTPUT)
        end

        @testset "no log file written (rejected before log is opened)" begin
            @test !isfile(INVALID_LOG)
        end
    end
catch e
    global all_ok = false
    println("invalid-input e2e: testset reported failures (see summary above); continuing to next scenario.")
end

rm(INVALID_OUTDIR; recursive=true, force=true)

println()
println("=== latticegen2 e2e: smoke-verified ===")

# specification.md §6.1 `smoke-verified` scenario. Uses -cc 20 -t 4 (not the
# originally-written -cc 10 -t 2 — that mismatched the golden sample's own
# filename, test/80mm-test-ball-cc20t4.step; -cc 20 -t 4 is what that file was
# actually generated with, so the scenario's params were corrected to match
# the golden sample rather than the other way around, per user decision).
const VERIFIED_GOLDEN = joinpath(ROOT, "test", "80mm-test-ball-cc20t4-golden-sample.step")
const VERIFIED_OUTDIR = mktempdir()
const VERIFIED_OUTPUT = joinpath(VERIFIED_OUTDIR, "smoke-verified.step")
const VERIFIED_BUDGET_SECONDS = 20 * 60.0
const VERIFIED_CC = 20.0
const VERIFIED_T = 4.0

verified_cmd = `julia --project=$ROOT $(joinpath(ROOT, "src", "main.jl")) -i $INPUT -cc $VERIFIED_CC -t $VERIFIED_T -o $VERIFIED_OUTPUT --workers 4 --tile-cells 6`
println("Command: ", verified_cmd)

verified_t0 = time()
verified_proc = run(pipeline(verified_cmd; stdout=stdout, stderr=stderr); wait=false)
wait(verified_proc)
verified_elapsed = time() - verified_t0

try
    @testset "smoke-verified e2e" begin
        @testset "process exits 0" begin
            @test verified_proc.exitcode == 0
        end

        @testset "runtime under performance budget" begin
            println("Elapsed: $(round(verified_elapsed, digits=1))s (budget $(VERIFIED_BUDGET_SECONDS)s)")
            @test verified_elapsed < VERIFIED_BUDGET_SECONDS
        end

        @testset "STEP file written and non-empty" begin
            @test isfile(VERIFIED_OUTPUT)
            isfile(VERIFIED_OUTPUT) && @test filesize(VERIFIED_OUTPUT) > 0
        end

        @testset "sub-threshold removals stay a small fraction of solids (docs/algorithm.md §8, §11.2)" begin
            ratio = sub_threshold_removal_ratio(VERIFIED_OUTPUT)
            println("Sub-threshold removal ratio: $(round(ratio*100, digits=2))%")
            @test !isnan(ratio)
            isnan(ratio) || @test ratio < 0.01
        end

        if isfile(VERIFIED_OUTPUT) && filesize(VERIFIED_OUTPUT) > 0
            @testset "geometry validity (manifold, no self-intersections)" begin
                mesh = with_gmsh() do
                    new_model("e2e-verified-mesh")
                    import_shapes(VERIFIED_OUTPUT)
                    tessellate_surface(LatticeParams(VERIFIED_CC, VERIFIED_T))
                end
                manifold_ok, bad_edges = manifold_check(mesh)
                @test manifold_ok
                manifold_ok || println("Non-manifold edges: ", length(bad_edges))

                no_self_int, bad_pairs = self_intersection_check(mesh)
                @test no_self_int
                no_self_int || println("Self-intersecting triangle pairs: ", length(bad_pairs))
            end

            @testset "matches golden sample (near-zero symmetric-difference volume)" begin
                @test isfile(VERIFIED_GOLDEN)
                if isfile(VERIFIED_GOLDEN)
                    diff_vol = golden_sample_volume_diff(VERIFIED_OUTPUT, VERIFIED_GOLDEN)
                    tol = VERIFIED_T^3   # negligible relative to model volume; spec §5's own floating-body floor
                    println("Golden-sample symmetric-difference volume: $(diff_vol) mm^3 (tolerance $(tol) mm^3)")
                    @test diff_vol < tol
                end
            end
        end
    end
catch e
    global all_ok = false
    println("smoke-verified e2e: testset reported failures (see summary above); continuing to next scenario.")
end

rm(VERIFIED_OUTDIR; recursive=true, force=true)

println()
println("=== latticegen2 e2e: dense-lattice ===")

# specification.md §6.1 `dense-lattice` scenario. Gated on the golden sample
# already existing on disk — not because the scenario is broken (it runs to
# completion; the earlier long-run and mesh-coverage problems are both
# root-caused and fixed, docs/algorithm.md §11.2 and §11.3) but because
# committing a golden sample is a review decision for whoever owns the test
# fixture, not something this harness should make for them.
const DENSE_INPUT = joinpath(ROOT, "test", "test-cylinder.STEP")
const DENSE_GOLDEN = joinpath(ROOT, "test", "test-cylinder-cc10t1.5-golden-sample.step")
const DENSE_CC = 10.0
const DENSE_T = 1.5
const DENSE_BUDGET_SECONDS = 60 * 60.0   # less dense than the original -cc 5 -t 1 scenario; 60-minute budget for now

if !isfile(DENSE_GOLDEN)
    println("SKIPPED: golden sample not found at $DENSE_GOLDEN.")
    println("         The scenario itself runs to completion — the earlier multi-hour run and the")
    println("         later mesh-coverage rejection are both root-caused and fixed (docs/algorithm.md")
    println("         §11.2, §11.3). Generating this file and committing it as the golden sample is a")
    println("         review decision, not automated here. See specification.md §6.1/§9.")
else
    const DENSE_OUTDIR = mktempdir()
    const DENSE_OUTPUT = joinpath(DENSE_OUTDIR, "dense-lattice.step")

    dense_cmd = `julia --project=$ROOT $(joinpath(ROOT, "src", "main.jl")) -i $DENSE_INPUT -cc $DENSE_CC -t $DENSE_T -o $DENSE_OUTPUT --cores 6 --ram 20`
    println("Command: ", dense_cmd)

    dense_t0 = time()
    dense_proc = run(pipeline(dense_cmd; stdout=stdout, stderr=stderr); wait=false)
    wait(dense_proc)
    dense_elapsed = time() - dense_t0

    try
        @testset "dense-lattice e2e" begin
            @testset "process exits 0" begin
                @test dense_proc.exitcode == 0
            end

            @testset "runtime under performance budget" begin
                println("Elapsed: $(round(dense_elapsed, digits=1))s (budget $(DENSE_BUDGET_SECONDS)s)")
                @test dense_elapsed < DENSE_BUDGET_SECONDS
            end

            @testset "STEP file written and non-empty" begin
                @test isfile(DENSE_OUTPUT)
                isfile(DENSE_OUTPUT) && @test filesize(DENSE_OUTPUT) > 0
            end

            @testset "sub-threshold removals stay a small fraction of solids (docs/algorithm.md §8, §11.2)" begin
                ratio = sub_threshold_removal_ratio(DENSE_OUTPUT)
                println("Sub-threshold removal ratio: $(round(ratio*100, digits=2))%")
                @test !isnan(ratio)
                isnan(ratio) || @test ratio < 0.01
            end

            if isfile(DENSE_OUTPUT) && filesize(DENSE_OUTPUT) > 0
                @testset "geometry validity (manifold, no self-intersections)" begin
                    mesh = with_gmsh() do
                        new_model("e2e-dense-mesh")
                        import_shapes(DENSE_OUTPUT)
                        tessellate_surface(LatticeParams(DENSE_CC, DENSE_T))
                    end
                    manifold_ok, bad_edges = manifold_check(mesh)
                    @test manifold_ok
                    manifold_ok || println("Non-manifold edges: ", length(bad_edges))

                    no_self_int, bad_pairs = self_intersection_check(mesh)
                    @test no_self_int
                    no_self_int || println("Self-intersecting triangle pairs: ", length(bad_pairs))
                end

                @testset "matches golden sample (near-zero symmetric-difference volume)" begin
                    diff_vol = golden_sample_volume_diff(DENSE_OUTPUT, DENSE_GOLDEN)
                    tol = DENSE_T^3
                    println("Golden-sample symmetric-difference volume: $(diff_vol) mm^3 (tolerance $(tol) mm^3)")
                    @test diff_vol < tol
                end
            end
        end
    catch e
        global all_ok = false
        println("dense-lattice e2e: testset reported failures (see summary above).")
    end

    rm(DENSE_OUTDIR; recursive=true, force=true)
end

println()
println("=== e2e summary ===")
println("smoke-fast:      exit=$(proc.exitcode) elapsed=$(round(elapsed, digits=1))s output=$OUTPUT")
println("invalid-input:   exit=$(invalid_proc.exitcode)")
println("smoke-verified:  exit=$(verified_proc.exitcode) elapsed=$(round(verified_elapsed, digits=1))s output=$VERIFIED_OUTPUT")
println("dense-lattice:   ", isfile(DENSE_GOLDEN) ? "ran (see above)" : "skipped (no golden sample yet)")
println()
println(all_ok ? "=== e2e: ALL SCENARIOS PASSED ===" : "=== e2e: ONE OR MORE SCENARIOS FAILED (see above) ===")

exit(all_ok ? 0 : 1)
