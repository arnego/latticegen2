# Unit tests for src/stepmeta.jl — pure Julia, no gmsh required (operates
# on hand-written STEP-header-shaped text).

const _SAMPLE_HEADER = """
ISO-10303-21;
HEADER;
FILE_NAME('Open CASCADE Shape Model','2026-08-08T10:06:07',('Author'),(
    'Open CASCADE'),'Open CASCADE STEP processor 7.9','Open CASCADE 7.9'
  ,'Unknown');
FILE_DESCRIPTION(('Open CASCADE Model'),'2;1');
FILE_SCHEMA((''));
ENDSEC;
DATA;
#1 = APPLICATION_PROTOCOL_DEFINITION('international standard',
  'automotive_design',2000,#2);
ENDSEC;
END-ISO-10303-21;
"""

@testset "stepmeta.jl" begin
    @testset "_find_entity_span locates entities across newlines" begin
        span = latticegen2._find_entity_span(_SAMPLE_HEADER, "FILE_NAME")
        @test span !== nothing
        text = _SAMPLE_HEADER[span]
        @test startswith(text, "FILE_NAME(")
        @test endswith(text, ";")
        @test occursin("Open CASCADE Shape Model", text)

        span2 = latticegen2._find_entity_span(_SAMPLE_HEADER, "FILE_DESCRIPTION")
        @test span2 !== nothing
        @test occursin("FILE_DESCRIPTION(", _SAMPLE_HEADER[span2])

        @test latticegen2._find_entity_span(_SAMPLE_HEADER, "NO_SUCH_ENTITY") === nothing
    end

    @testset "_replace_first_string_field" begin
        entity = "FILE_NAME('old name','ts',('Author'),('Org'),'pre','orig','auth');"
        replaced = latticegen2._replace_first_string_field(entity, "new+part+name")
        @test occursin("'new+part+name'", replaced)
        @test !occursin("old name", replaced)
        @test occursin("'ts'", replaced)   # other fields untouched

        @testset "escapes embedded quotes" begin
            r2 = latticegen2._replace_first_string_field(entity, "O'Brien")
            @test occursin("'O''Brien'", r2)
        end
    end

    @testset "_append_to_first_list" begin
        entity = "FILE_DESCRIPTION(('Open CASCADE Model'),'2;1');"
        appended = latticegen2._append_to_first_list(entity, "extra text")
        @test occursin("'Open CASCADE Model'", appended)
        @test occursin("'extra text'", appended)
        @test endswith(appended, ";")
    end

    @testset "_is_blank_schema" begin
        @test latticegen2._is_blank_schema("FILE_SCHEMA((''));")
        @test !latticegen2._is_blank_schema("FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));")
    end

    @testset "generation_params_text" begin
        txt = generation_params_text("test/80mm-test-ball.step", 10.0, 2.0)
        @test occursin("cc=10 mm", txt)
        @test occursin("t=2 mm", txt)
        @test occursin("80mm-test-ball.step", txt)
    end

    @testset "rewrite_step_header! end-to-end on a real file" begin
        mktempdir() do dir
            path = joinpath(dir, "sample.step")
            write(path, _SAMPLE_HEADER)
            rewrite_step_header!(path, "ball+cc10+t2", "latticegen2 test params cc=10 t=2")
            text = read(path, String)

            @test occursin("ball+cc10+t2", text)
            @test !occursin("Open CASCADE Shape Model", text)
            @test occursin("latticegen2 test params cc=10 t=2", text)
            @test occursin("Open CASCADE Model", text)     # original description entry preserved
            @test occursin("'AUTOMOTIVE_DESIGN'", text)     # blank schema filled in
            # Untouched sections still present and well-formed:
            @test occursin("APPLICATION_PROTOCOL_DEFINITION", text)
            @test occursin("END-ISO-10303-21;", text)
        end
    end

    @testset "rewrite_step_header! does not overwrite a non-blank schema" begin
        mktempdir() do dir
            path = joinpath(dir, "sample2.step")
            populated = replace(_SAMPLE_HEADER, "FILE_SCHEMA((''));" => "FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));")
            write(path, populated)
            rewrite_step_header!(path, "part", "params")
            text = read(path, String)
            @test occursin("'CONFIG_CONTROL_DESIGN'", text)
            @test !occursin("AUTOMOTIVE_DESIGN", text)
        end
    end
end
