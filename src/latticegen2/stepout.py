"""STEP export: header metadata rewrite and the round-trip success gate.

Ported from the Julia implementation's ``src/stepmeta.jl`` (docs/algorithm.md
§8). STEP (ISO-10303-21) headers are plain ASCII entities of the form
``NAME(arg, ...);`` where arguments nest and quoted strings escape a literal
``'`` by doubling it. The helpers below are small quote-aware scanners —
deliberately not a general STEP parser, just enough to edit the three header
entities this module touches without disturbing anything else.

``FILE_SCHEMA``'s value is only ever filled in when blank, never overwritten:
that is what keeps the output a clean standard AP214 document rather than a
hand-patched hybrid.
"""

from __future__ import annotations

import datetime as _dt

from . import occ
from .errors import OutputError
from .lattice import format_param


def _find_entity_span(text: str, name: str) -> tuple[int, int] | None:
    """``(start, end)`` character span of ``name(...);``, or ``None``.

    Respects STEP's quoted-string syntax so parentheses and semicolons inside
    string literals are not mistaken for structural delimiters.
    """
    start = text.find(name + "(")
    if start < 0:
        return None
    i = text.index("(", start) + 1
    depth = 1
    in_string = False
    n = len(text)
    while i < n:
        ch = text[i]
        if in_string:
            if ch == "'":
                if i + 1 < n and text[i + 1] == "'":
                    i += 1
                else:
                    in_string = False
        elif ch == "'":
            in_string = True
        elif ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                semi = text.find(";", i)
                if semi < 0:
                    return None
                return start, semi + 1
        i += 1
    return None


def _matching_close_quote(text: str, q1: int) -> int | None:
    i = q1 + 1
    n = len(text)
    while i < n:
        if text[i] == "'":
            if i + 1 < n and text[i + 1] == "'":
                i += 2
                continue
            return i
        i += 1
    return None


def _replace_first_string_field(entity: str, new_value: str) -> str:
    """Replace the first quoted string in an entity's argument list."""
    open_paren = entity.find("(")
    if open_paren < 0:
        raise OutputError(f"Malformed STEP entity (no '('): {entity}")
    q1 = entity.find("'", open_paren)
    if q1 < 0:
        raise OutputError(f"Malformed STEP entity (no string field): {entity}")
    q2 = _matching_close_quote(entity, q1)
    if q2 is None:
        raise OutputError(f"Malformed STEP entity (unterminated string): {entity}")
    return entity[: q1 + 1] + new_value.replace("'", "''") + entity[q2:]


def _append_to_first_list(entity: str, extra: str) -> str:
    """Append ``extra`` as a new quoted string to the first nested list argument."""
    open1 = entity.find("(")
    if open1 < 0:
        raise OutputError(f"Malformed STEP entity (no '('): {entity}")
    open2 = entity.find("(", open1 + 1)
    if open2 < 0:
        raise OutputError(f"Malformed STEP entity (no nested list): {entity}")
    i = open2 + 1
    depth = 1
    in_string = False
    n = len(entity)
    close = None
    while i < n:
        ch = entity[i]
        if in_string:
            if ch == "'":
                if i + 1 < n and entity[i + 1] == "'":
                    i += 1
                else:
                    in_string = False
        elif ch == "'":
            in_string = True
        elif ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                close = i
                break
        i += 1
    if close is None:
        raise OutputError(f"Malformed STEP entity (unterminated list): {entity}")
    return entity[:close] + ",'" + extra.replace("'", "''") + "'" + entity[close:]


def _is_blank_schema(entity: str) -> bool:
    import re

    return re.search(r"\(\s*'\s*'\s*\)", entity) is not None


def generation_params_text(input_path: str, cc: float, t: float) -> str:
    """The parameter string embedded into FILE_DESCRIPTION (specification.md §5)."""
    stamp = _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return (
        f"latticegen2 generation parameters: input={input_path}, "
        f"cc={format_param(cc)} mm, t={format_param(t)} mm, generated={stamp}"
    )


def rewrite_step_header(path: str, part_name: str, params_text: str) -> None:
    """Set FILE_NAME's name, append the parameters to FILE_DESCRIPTION, fill FILE_SCHEMA.

    ``FILE_SCHEMA`` is filled with ``'AUTOMOTIVE_DESIGN'`` only when it is blank
    — a populated value is never overwritten.
    """
    with open(path, "r", encoding="utf-8", errors="surrogateescape") as fh:
        text = fh.read()

    span = _find_entity_span(text, "FILE_NAME")
    if span is None:
        raise OutputError(f"STEP file has no FILE_NAME entity: {path}")
    start, end = span
    text = text[:start] + _replace_first_string_field(text[start:end], part_name) + text[end:]

    span = _find_entity_span(text, "FILE_DESCRIPTION")
    if span is None:
        raise OutputError(f"STEP file has no FILE_DESCRIPTION entity: {path}")
    start, end = span
    text = text[:start] + _append_to_first_list(text[start:end], params_text) + text[end:]

    span = _find_entity_span(text, "FILE_SCHEMA")
    if span is not None:
        start, end = span
        current = text[start:end]
        if _is_blank_schema(current):
            text = text[:start] + _replace_first_string_field(current, "AUTOMOTIVE_DESIGN") + text[end:]

    with open(path, "w", encoding="utf-8", errors="surrogateescape") as fh:
        fh.write(text)


def round_trip_check(path: str) -> int:
    """Re-read the written file and confirm it still contains solids.

    The mandatory success gate: a run only exits 0 after the file it just wrote
    has been parsed back independently (docs/algorithm.md §8).
    """
    shape = occ.read_step(path)
    found = occ.solids(shape)
    if not found:
        raise OutputError(f"Round-trip re-import of {path} found no solids.")
    return len(found)
