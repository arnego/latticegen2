# Unit tests for src/cli.jl — pure Julia, no gmsh required.

@testset "cli.jl" begin
    @testset "valid explicit args" begin
        a = parse_args(["-i", "test/80mm-test-ball.step", "-cc", "10", "-t", "2",
                         "--workers", "4", "--tile-cells", "8"])
        @test a.input == "test/80mm-test-ball.step"
        @test a.cc == 10.0
        @test a.t == 2.0
        @test a.workers == 4
        @test a.tile_cells == 8
        @test a.background == false
        @test a.verbose == false
        @test endswith(a.output, "80mm-test-ball-cc10t2.step")
        @test endswith(a.log_path, "80mm-test-ball-cc10t2.log")
    end

    @testset "valid auto-tune args" begin
        a = parse_args(["-i", "x.step", "-cc", "5", "-t", "1", "--cores", "6", "--ram", "24"])
        @test a.cores == 6
        @test a.ram == 24.0
    end

    @testset "custom output path" begin
        # A user-supplied -o path is preserved verbatim (plus the .step/.log
        # extension) rather than being re-normalized through joinpath.
        a = parse_args(["-i", "x.step", "-cc", "5", "-t", "1", "--cores", "6", "--ram", "24",
                         "-o", "out/result"])
        @test a.output == "out/result.step"
        @test a.log_path == "out/result.log"
    end

    @testset "flags" begin
        a = parse_args(["-i", "x.step", "-cc", "5", "-t", "1", "--cores", "6", "--ram", "24",
                         "-bg", "-v"])
        @test a.background
        @test a.verbose
    end

    @testset "missing required" begin
        @test_throws ParamError parse_args(["-cc", "5", "-t", "1", "--cores", "6", "--ram", "24"])
        @test_throws ParamError parse_args(["-i", "x.step", "-t", "1", "--cores", "6", "--ram", "24"])
        @test_throws ParamError parse_args(["-i", "x.step", "-cc", "5", "--cores", "6", "--ram", "24"])
    end

    @testset "out-of-range" begin
        @test_throws ParamError parse_args(["-i", "x.step", "-cc", "0.1", "-t", "1", "--cores", "6", "--ram", "24"])
        @test_throws ParamError parse_args(["-i", "x.step", "-cc", "5", "-t", "50", "--cores", "6", "--ram", "24"])
        @test_throws ParamError parse_args(["-i", "x.step", "-cc", "5", "-t", "1", "--cores", "200", "--ram", "24"])
    end

    @testset "t vs cell edge cross-constraint" begin
        # a = cc/sqrt2; for cc=5, a ≈ 3.5355; t=4 should fail (t >= a)
        @test_throws ParamError parse_args(["-i", "x.step", "-cc", "5", "-t", "4", "--cores", "6", "--ram", "24"])
        a = parse_args(["-i", "x.step", "-cc", "5", "-t", "3", "--cores", "6", "--ram", "24"])
        @test a.t == 3.0
    end

    @testset "mutual exclusivity of optimization params" begin
        @test_throws ParamError parse_args(["-i", "x.step", "-cc", "5", "-t", "1"])
        @test_throws ParamError parse_args(["-i", "x.step", "-cc", "5", "-t", "1",
                                             "--cores", "6", "--ram", "24", "--workers", "2", "--tile-cells", "8"])
        @test_throws ParamError parse_args(["-i", "x.step", "-cc", "5", "-t", "1", "--cores", "6"])
        @test_throws ParamError parse_args(["-i", "x.step", "-cc", "5", "-t", "1", "--workers", "4"])
    end

    @testset "unknown flag" begin
        @test_throws ParamError parse_args(["-i", "x.step", "-cc", "5", "-t", "1",
                                             "--cores", "6", "--ram", "24", "--bogus"])
    end

    @testset "resolve_output_paths formatting" begin
        step, log = resolve_output_paths("dir/name.step", nothing, 2.5, 0.4)
        @test step == joinpath("dir", "name-cc2.5t0.4.step")
        @test log == joinpath("dir", "name-cc2.5t0.4.log")

        step2, _ = resolve_output_paths("name.step", nothing, 10.0, 1.0)
        @test step2 == joinpath(".", "name-cc10t1.step")
    end

    @testset "preflight_checks" begin
        mktempdir() do dir
            infile = joinpath(dir, "in.step")
            write(infile, "ISO-10303-21;\nENDSEC;\n")
            args = parse_args(["-i", infile, "-cc", "5", "-t", "1", "--cores", "2", "--ram", "8",
                                "-o", joinpath(dir, "out")])
            @test preflight_checks(args) === nothing

            missing_args = parse_args(["-i", joinpath(dir, "nope.step"), "-cc", "5", "-t", "1",
                                        "--cores", "2", "--ram", "8"])
            @test_throws InputGeometryError preflight_checks(missing_args)

            baddir_args = parse_args(["-i", infile, "-cc", "5", "-t", "1", "--cores", "2", "--ram", "8",
                                       "-o", joinpath(dir, "no_such_dir", "out")])
            @test_throws OutputError preflight_checks(baddir_args)
        end
    end
end
