# Unit tests for src/runlog.jl — pure Julia, no gmsh required.

@testset "runlog.jl" begin
    @testset "format_param" begin
        @test format_param(10.0) == "10"
        @test format_param(2.5) == "2.5"
        @test format_param(0.4) == "0.4"
        @test format_param(1) == "1"
        @test format_param(0.0) == "0"
    end

    @testset "format_bytes" begin
        @test format_bytes(500 * 1024^2) == "500.0 MB"
        @test occursin("GB", format_bytes(2 * 1024^3))
    end

    @testset "format_duration" begin
        @test occursin("s", format_duration(5.5))
        @test occursin("m", format_duration(90.0))
        @test occursin("h", format_duration(3700.0))
    end

    @testset "exit codes" begin
        @test exit_code(ParamError("x")) == 2
        @test exit_code(InputGeometryError("x")) == 3
        @test exit_code(ProcessingError("x")) == 4
        @test exit_code(ResourceError("x")) == 5
        @test exit_code(OutputError("x")) == 6
        @test exit_code(CancelledError("x")) == 130
        @test exit_code(ErrorException("x")) == 4
    end

    @testset "is_interrupt unwraps the shapes Ctrl+C arrives in" begin
        @test is_interrupt(InterruptException())
        @test is_interrupt(CancelledError("cancelled"))

        # Raised inside a remote call: RemoteException(pid, CapturedException).
        # Constructed structurally here (Distributed is not loaded in the fast
        # test subset) — is_interrupt matches on the `captured` field for the
        # same reason.
        struct_like_remote = (captured = CapturedException(InterruptException(), Any[]),)
        @test is_interrupt(struct_like_remote)
        @test is_interrupt(CapturedException(InterruptException(), Any[]))

        # Raised inside an @async task under @sync (dispatch_tiles' pool):
        # CompositeException wrapping a TaskFailedException.
        composite = try
            @sync @async throw(InterruptException())
        catch e
            e
        end
        @test is_interrupt(composite)

        failed_task = try
            wait(@async throw(InterruptException()))
        catch e
            e
        end
        @test is_interrupt(failed_task)

        # Genuine failures must not be mistaken for a cancellation, including
        # a worker dying for a non-interrupt reason.
        @test !is_interrupt(ErrorException("boom"))
        @test !is_interrupt(ProcessingError("fuse failed"))
        @test !is_interrupt(CapturedException(ErrorException("boom"), Any[]))
        @test !is_interrupt((captured = CapturedException(ErrorException("oom"), Any[]),))
    end

    @testset "RunLog lifecycle" begin
        mktempdir() do dir
            path = joinpath(dir, "test.log")
            rl = open_runlog(path; verbose=false)
            log_line(rl, "hello")
            log_stage(rl, "stage1", 1.23)
            log_run_header(rl, Dict("cc" => 5.0, "t" => 1.0))
            log_summary(rl, Dict("cc" => 5.0), Dict{Symbol,Any}(:output_path => "out.step"))
            close_runlog(rl)
            @test isfile(path)
            content = read(path, String)
            @test occursin("hello", content)
            @test occursin("stage1", content)
            @test occursin("run summary", content)
            @test occursin("out.step", content)
        end
    end
end
