# src/stepmeta.jl
#
# Post-export STEP header text rewrite and round-trip verification.
# Implements docs/algorithm.md §8.
#
# STEP (ISO-10303-21) header entities are plain ASCII text of the form
# `NAME(arg1, arg2, ...);`, where arguments can themselves be nested lists
# and quoted strings use `''` to escape a literal `'`. The helpers below are
# small, quote-aware scanners — deliberately not a general STEP parser, just
# enough to locate and surgically edit the three header entities this
# module touches, without disturbing anything else in the file (in
# particular, FILE_SCHEMA's *value* is only ever filled in when blank, never
# overwritten — see `rewrite_step_header!`).

using Dates

"""
    _find_entity_span(text, name) -> Union{UnitRange{Int},Nothing}

Character range of the `name(...)` entity call in `text`, from `name`'s
first character through the terminating `;`, respecting STEP's quoted
string syntax so parentheses/semicolons inside string literals are not
mistaken for structural delimiters. `nothing` if `name(` doesn't occur.
"""
function _find_entity_span(text::String, name::String)
    startpos = findfirst(name * "(", text)
    startpos === nothing && return nothing
    open_paren = last(startpos)
    i = nextind(text, open_paren)
    depth = 1
    in_string = false
    n = lastindex(text)
    while i <= n
        c = text[i]
        if in_string
            if c == '\''
                if i < n && text[nextind(text, i)] == '\''
                    i = nextind(text, i)   # doubled '' -> literal quote, skip both
                else
                    in_string = false
                end
            end
        else
            if c == '\''
                in_string = true
            elseif c == '('
                depth += 1
            elseif c == ')'
                depth -= 1
                if depth == 0
                    j = findnext(';', text, i)
                    j === nothing && return nothing
                    return first(startpos):j
                end
            end
        end
        i = nextind(text, i)
    end
    return nothing
end

"""Find the closing quote index matching the opening quote at `text[q1]`
(quote-aware: doubled `''` is a literal quote, not a terminator)."""
function _matching_close_quote(text::String, q1::Int)
    i = nextind(text, q1)
    n = lastindex(text)
    while i <= n
        if text[i] == '\''
            if i < n && text[nextind(text, i)] == '\''
                i = nextind(text, nextind(text, i))
                continue
            else
                return i
            end
        end
        i = nextind(text, i)
    end
    return nothing
end

"""Replace the *first* quoted string field in `entity_text`'s argument list
(e.g. FILE_NAME's `name` field) with `newvalue` (escaping embedded `'` as
`''` per STEP convention)."""
function _replace_first_string_field(entity_text::String, newvalue::AbstractString)
    open_paren = findfirst('(', entity_text)
    open_paren === nothing && throw(OutputError("Malformed STEP entity (no '('): $entity_text"))
    q1 = findnext('\'', entity_text, open_paren)
    q1 === nothing && throw(OutputError("Malformed STEP entity (no string field): $entity_text"))
    q2 = _matching_close_quote(entity_text, q1)
    q2 === nothing && throw(OutputError("Malformed STEP entity (unterminated string): $entity_text"))
    escaped = replace(String(newvalue), "'" => "''")
    return entity_text[1:q1] * escaped * entity_text[q2:end]
end

"""Append `extra` as a new quoted string to the first nested list argument
of `entity_text` (e.g. FILE_DESCRIPTION's `(description...)` list)."""
function _append_to_first_list(entity_text::String, extra::AbstractString)
    open1 = findfirst('(', entity_text)
    open1 === nothing && throw(OutputError("Malformed STEP entity (no '('): $entity_text"))
    open2 = findnext('(', entity_text, nextind(entity_text, open1))
    open2 === nothing && throw(OutputError("Malformed STEP entity (no nested list): $entity_text"))
    i = nextind(entity_text, open2)
    depth = 1
    in_string = false
    n = lastindex(entity_text)
    closepos = nothing
    while i <= n
        c = entity_text[i]
        if in_string
            if c == '\''
                if i < n && entity_text[nextind(entity_text, i)] == '\''
                    i = nextind(entity_text, i)
                else
                    in_string = false
                end
            end
        else
            if c == '\''
                in_string = true
            elseif c == '('
                depth += 1
            elseif c == ')'
                depth -= 1
                if depth == 0
                    closepos = i
                    break
                end
            end
        end
        i = nextind(entity_text, i)
    end
    closepos === nothing && throw(OutputError("Malformed STEP entity (unterminated list): $entity_text"))
    escaped = replace(String(extra), "'" => "''")
    insertion = ",'" * escaped * "'"
    return entity_text[1:prevind(entity_text, closepos)] * insertion * entity_text[closepos:end]
end

"""True if a `FILE_SCHEMA` entity's argument list contains only a single
blank string (gmsh/OCCT's default STEP output — see `rewrite_step_header!`)."""
_is_blank_schema(entity_text::String) = occursin(r"\(\s*'\s*'\s*\)", entity_text)

"""Generation-parameters text embedded into FILE_DESCRIPTION
(specification.md §5)."""
function generation_params_text(input::AbstractString, cc::Real, t::Real)
    return "latticegen2 generation parameters: input=$(input), cc=$(format_param(cc)) mm, " *
           "t=$(format_param(t)) mm, generated=$(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))"
end

"""
    rewrite_step_header!(path, part_name, params_text)

Post-export text rewrite of a gmsh/OCCT-written STEP file's HEADER section
(docs/algorithm.md §8): sets `FILE_NAME`'s `name` field to `part_name`,
appends `params_text` as an extra `FILE_DESCRIPTION` entry, and fills in
`FILE_SCHEMA` with `'AUTOMOTIVE_DESIGN'` **only if it is currently blank**
(gmsh's default STEP export leaves it blank even though the file's
`APPLICATION_PROTOCOL_DEFINITION` already correctly identifies AP214 — see
docs/algorithm.md §8). A non-blank `FILE_SCHEMA` is never overwritten.
Rewrites the file in place.
"""
function rewrite_step_header!(path::AbstractString, part_name::AbstractString, params_text::AbstractString)
    text = read(path, String)

    fn_span = _find_entity_span(text, "FILE_NAME")
    fn_span === nothing && throw(OutputError("STEP file has no FILE_NAME entity: $path"))
    new_fn = _replace_first_string_field(text[fn_span], part_name)
    text = text[1:prevind(text, first(fn_span))] * new_fn * text[nextind(text, last(fn_span)):end]

    fd_span = _find_entity_span(text, "FILE_DESCRIPTION")
    fd_span === nothing && throw(OutputError("STEP file has no FILE_DESCRIPTION entity: $path"))
    new_fd = _append_to_first_list(text[fd_span], params_text)
    text = text[1:prevind(text, first(fd_span))] * new_fd * text[nextind(text, last(fd_span)):end]

    fs_span = _find_entity_span(text, "FILE_SCHEMA")
    if fs_span !== nothing
        current = text[fs_span]
        if _is_blank_schema(current)
            new_fs = _replace_first_string_field(current, "AUTOMOTIVE_DESIGN")
            text = text[1:prevind(text, first(fs_span))] * new_fs * text[nextind(text, last(fs_span)):end]
        end
    end

    write(path, text)
    return nothing
end

"""
    round_trip_check(path)

Re-import `path` into a fresh gmsh session and confirm at least one volume
is present, throwing `OutputError` otherwise (docs/algorithm.md §8's
mandatory success gate — a run is only reported as `0`/success after this
passes).
"""
function round_trip_check(path::AbstractString)
    with_gmsh() do
        new_model("verify")
        vols = import_shapes(path)
        isempty(vols) && throw(OutputError("Round-trip re-import of $path found no volumes"))
    end
    return true
end
