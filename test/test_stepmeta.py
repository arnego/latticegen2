"""STEP header rewriting (docs/algorithm.md §9)."""

import pytest

from latticegen2.errors import OutputError
from latticegen2.stepout import (
    _append_to_first_list,
    _find_entity_span,
    _is_blank_schema,
    _replace_first_string_field,
    generation_params_text,
    rewrite_step_header,
)

HEADER = """ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('Open CASCADE Model'),'2;1');
FILE_NAME('Open CASCADE Shape Model','2026-08-10T12:00:00',('Author'),(
    'Open CASCADE'),'Open CASCADE STEP processor 7.9','','Unknown');
FILE_SCHEMA(('
'));
ENDSEC;
DATA;
#1=SOME_ENTITY('x;y(z)');
ENDSEC;
END-ISO-10303-21;
"""


def test_find_entity_span_locates_a_whole_entity():
    span = _find_entity_span(HEADER, "FILE_DESCRIPTION")
    assert span is not None
    text = HEADER[span[0]:span[1]]
    assert text.startswith("FILE_DESCRIPTION(")
    assert text.endswith(";")


def test_find_entity_span_ignores_delimiters_inside_strings():
    text = "FILE_NAME('a;b(c)','d');NEXT;"
    span = _find_entity_span(text, "FILE_NAME")
    assert text[span[0]:span[1]] == "FILE_NAME('a;b(c)','d');"


def test_find_entity_span_returns_none_when_absent():
    assert _find_entity_span("HEADER;", "FILE_NAME") is None


def test_replace_first_string_field():
    out = _replace_first_string_field("FILE_NAME('old','rest');", "new")
    assert out == "FILE_NAME('new','rest');"


def test_replace_first_string_field_escapes_quotes():
    out = _replace_first_string_field("FILE_NAME('old');", "it's")
    assert out == "FILE_NAME('it''s');"


def test_replace_first_string_field_rejects_malformed_input():
    with pytest.raises(OutputError):
        _replace_first_string_field("FILE_NAME(1,2);", "x")


def test_append_to_first_list():
    out = _append_to_first_list("FILE_DESCRIPTION(('a'),'2;1');", "b")
    assert out == "FILE_DESCRIPTION(('a','b'),'2;1');"


def test_blank_schema_detection():
    assert _is_blank_schema("FILE_SCHEMA(('\n'));")
    assert _is_blank_schema("FILE_SCHEMA((''));")
    assert not _is_blank_schema("FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));")


def test_rewrite_sets_name_appends_params_and_fills_blank_schema(tmp_path):
    path = tmp_path / "out.step"
    path.write_text(HEADER, encoding="utf-8")
    rewrite_step_header(str(path), "ball+cc20+t4", "params here")
    text = path.read_text(encoding="utf-8")
    assert "FILE_NAME('ball+cc20+t4'" in text
    assert "'params here'" in text
    assert "AUTOMOTIVE_DESIGN" in text
    assert "#1=SOME_ENTITY('x;y(z)');" in text, "the DATA section must be untouched"


def test_rewrite_never_overwrites_a_populated_schema(tmp_path):
    path = tmp_path / "out.step"
    path.write_text(HEADER.replace("FILE_SCHEMA(('\n'));", "FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));"),
                    encoding="utf-8")
    rewrite_step_header(str(path), "n", "p")
    assert "CONFIG_CONTROL_DESIGN" in path.read_text(encoding="utf-8")
    assert "AUTOMOTIVE_DESIGN" not in path.read_text(encoding="utf-8")


def test_rewrite_rejects_a_file_without_a_header(tmp_path):
    path = tmp_path / "bad.step"
    path.write_text("ISO-10303-21;\nDATA;\nENDSEC;\n", encoding="utf-8")
    with pytest.raises(OutputError, match="FILE_NAME"):
        rewrite_step_header(str(path), "n", "p")


def test_generation_params_text_carries_the_parameters():
    text = generation_params_text("/in/part.step", 20.0, 4.0)
    assert "cc=20 mm" in text and "t=4 mm" in text and "/in/part.step" in text
