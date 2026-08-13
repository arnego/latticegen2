"""CLI parsing and validation (specification.md §3).

Covers the ranges, the ``t < a`` cross-constraint, path resolution, and the rule
that a rejected run touches no files.
"""

import os

import pytest

from latticegen2.cli import (
    HelpRequested,
    default_workers,
    parse_args,
    preflight_checks,
    resolve_output_paths,
)
from latticegen2.errors import InputGeometryError, OutputError, ParamError

BASE = ["-i", "part.step", "-cc", "10", "-t", "1.5"]


def test_parses_the_required_arguments():
    args = parse_args(BASE)
    assert args.input == "part.step"
    assert args.cc == 10.0
    assert args.t == 1.5
    assert args.workers >= 1


def test_missing_input_is_rejected():
    with pytest.raises(ParamError, match="input"):
        parse_args(["-cc", "10", "-t", "1.5"])


def test_missing_cc_is_rejected():
    with pytest.raises(ParamError, match="-cc"):
        parse_args(["-i", "part.step", "-t", "1.5"])


def test_missing_t_is_rejected():
    with pytest.raises(ParamError, match="-t"):
        parse_args(["-i", "part.step", "-cc", "10"])


def test_unknown_argument_is_rejected():
    with pytest.raises(ParamError, match="Unknown argument"):
        parse_args(BASE + ["--nope"])


def test_missing_value_is_rejected():
    with pytest.raises(ParamError, match="Missing value"):
        parse_args(["-i", "part.step", "-cc", "10", "-t"])


def test_non_numeric_value_is_rejected():
    with pytest.raises(ParamError, match="Invalid numeric"):
        parse_args(["-i", "part.step", "-cc", "ten", "-t", "1.5"])


@pytest.mark.parametrize("cc", [0.39, 50.1])
def test_cc_out_of_range(cc):
    with pytest.raises(ParamError, match="-cc"):
        parse_args(["-i", "p.step", "-cc", str(cc), "-t", "0.4"])


@pytest.mark.parametrize("cc,t", [(50.0, 0.39), (50.0, 20.1)])
def test_t_out_of_range(cc, t):
    with pytest.raises(ParamError, match="-t"):
        parse_args(["-i", "p.step", "-cc", str(cc), "-t", str(t)])


def test_cc_and_t_at_their_documented_bounds_are_accepted():
    args = parse_args(["-i", "p.step", "-cc", "50", "-t", "20"])
    assert (args.cc, args.t) == (50.0, 20.0)
    assert parse_args(["-i", "p.step", "-cc", "0.6", "-t", "0.4"]).t == 0.4


def test_t_at_or_above_the_cell_edge_is_rejected():
    # cc = 5 gives a = 3.5355 mm, so t = 4 cannot fit inside one cell.
    with pytest.raises(ParamError, match="cube edge length"):
        parse_args(["-i", "p.step", "-cc", "5", "-t", "4"])


def test_t_below_the_cell_edge_is_accepted_right_up_to_the_limit():
    """``t < a`` is the only cross-constraint; the fast path adds none of its own.

    An earlier revision also rejected ``t >= cc/2`` on the assumption that the
    junction's mid-strut interfaces would be swallowed by the node overlap. They
    are not (test_junction.py), and rejecting valid input is not an acceptable
    failure mode — see docs/algorithm.md's history of the same mistake.
    """
    a = 10.0 / (2 ** 0.5)
    assert parse_args(["-i", "p.step", "-cc", "10", "-t", "7"]).t == 7.0
    assert 7.0 < a
    with pytest.raises(ParamError, match="cube edge length"):
        parse_args(["-i", "p.step", "-cc", "10", "-t", "7.1"])


def test_help_is_signalled_not_an_error():
    with pytest.raises(HelpRequested):
        parse_args(["-h"])


def test_flags_are_recognised():
    args = parse_args(BASE + ["-bg", "-v"])
    assert args.background and args.verbose


def test_worker_count_can_be_given_explicitly():
    assert parse_args(BASE + ["--workers", "3"]).workers == 3


def test_worker_count_is_derived_from_cores_when_not_given():
    """One worker per core: the master is blocked on ``Pool.imap`` throughout.

    It used to reserve a core for the master, which cost a fifth of the boundary
    stage's throughput on a 6-core machine to protect a process that spends over
    99 % of that stage waiting.
    """
    assert parse_args(BASE + ["--cores", "6"]).workers == 6
    assert parse_args(BASE + ["--cores", "1"]).workers == 1


def test_worker_count_is_capped():
    assert parse_args(BASE + ["--cores", "128"]).workers == 8
    assert default_workers(None) >= 1


def test_ram_hint_is_recorded():
    assert parse_args(BASE + ["--ram", "20"]).ram == 20.0


@pytest.mark.parametrize("flag,value", [("--workers", "0"), ("--cores", "0"), ("--ram", "0.5")])
def test_hint_ranges_are_validated(flag, value):
    with pytest.raises(ParamError, match=flag):
        parse_args(BASE + [flag, value])


def test_default_output_path_is_derived_from_the_input():
    step, log = resolve_output_paths(os.path.join("dir", "ball.step"), None, 20.0, 4.0)
    assert os.path.basename(step) == "ball-cc20t4.step"
    assert os.path.basename(log) == "ball-cc20t4.log"
    assert os.path.dirname(step) == "dir"


def test_explicit_output_gets_a_step_extension_and_a_matching_log():
    step, log = resolve_output_paths("in.step", os.path.join("out", "thing"), 10.0, 1.5)
    assert step.endswith("thing.step")
    assert log.endswith("thing.log")


def test_log_never_becomes_step_dot_log():
    _, log = resolve_output_paths("in.step", "out.step", 10.0, 1.5)
    assert log == "out.log"


# Built from the platform's own separators rather than hard-coded. A backslash
# is a directory separator on Windows but an ordinary filename character on
# POSIX, so `out\` must be rejected on one and accepted on the other — listing
# it unconditionally fails the Linux leg of the release workflow.
DIRECTORY_LIKE = ["." + os.sep, ".", "..", "out" + os.sep, "", "   "]
if os.altsep:
    DIRECTORY_LIKE += ["." + os.altsep, "out" + os.altsep]


@pytest.mark.parametrize("bad", DIRECTORY_LIKE)
def test_output_that_names_a_directory_is_rejected(bad):
    """``-o`` names a file. A directory used to be silently turned into one.

    ``-o .\\`` appended the extension to nothing and wrote ``.\\.step`` beside
    ``.\\.step.log`` — a basename that is all extension has no stem for the log
    to share, so both files landed under names nobody typed.
    """
    with pytest.raises(ParamError, match="must name the output .step file"):
        resolve_output_paths("in.step", bad, 10.0, 1.5)


@pytest.mark.skipif(os.altsep is not None, reason="backslash is a separator here")
def test_a_backslash_is_an_ordinary_filename_character_on_posix():
    """The rejection is about separators, and which characters those are is
    platform-dependent. On POSIX ``out\\`` is a legal filename, not a directory."""
    step, _ = resolve_output_paths("in.step", "out\\", 10.0, 1.5)
    assert step == "out\\.step"


@pytest.mark.parametrize(
    "good,stem",
    [
        ("hx", "hx"),
        ("hx.step", "hx"),
        (".hidden.step", ".hidden"),
        (os.path.join("a.b", "hx"), "hx"),
    ],
)
def test_ordinary_output_names_are_still_accepted(good, stem):
    """The rejection is syntactic, so it must not catch a real filename.

    Every valid basename has at least one non-dot character, including leading-
    dot names and paths whose *directory* carries a dot.
    """
    step, log = resolve_output_paths("in.step", good, 10.0, 1.5)
    assert os.path.basename(step) == stem + ".step"
    assert os.path.basename(log) == stem + ".log"


def test_preflight_rejects_an_output_that_is_an_existing_directory(tmp_path):
    """The syntactic rule cannot see the filesystem; this case needs a look."""
    target = tmp_path / "already.step"
    target.mkdir()
    args = parse_args(["-i", "in.step", "-cc", "10", "-t", "1.5", "-o", str(target)])
    with pytest.raises(ParamError, match="existing directory"):
        preflight_checks(args)


def test_preflight_rejects_a_missing_input(tmp_path):
    args = parse_args(["-i", str(tmp_path / "nope.step"), "-cc", "10", "-t", "1.5",
                       "-o", str(tmp_path / "out.step")])
    with pytest.raises(InputGeometryError, match="not found"):
        preflight_checks(args)


def test_preflight_rejects_a_missing_output_directory(tmp_path):
    src = tmp_path / "in.step"
    src.write_text("ISO-10303-21;")
    args = parse_args(["-i", str(src), "-cc", "10", "-t", "1.5",
                       "-o", str(tmp_path / "nodir" / "out.step")])
    with pytest.raises(OutputError, match="does not exist"):
        preflight_checks(args)


def test_preflight_accepts_a_writable_target(tmp_path):
    src = tmp_path / "in.step"
    src.write_text("ISO-10303-21;")
    args = parse_args(["-i", str(src), "-cc", "10", "-t", "1.5",
                       "-o", str(tmp_path / "out.step")])
    preflight_checks(args)
