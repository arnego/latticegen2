"""CLI parsing and validation (specification.md §3).

Covers the ranges, the ``t < a`` cross-constraint, path resolution, and the rule
that a rejected run touches no files.
"""

import math
import os
import re

import pytest

from latticegen2 import cli, sysinfo
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
    assert parse_args(BASE + ["-v"]).verbose


@pytest.mark.parametrize(
    "argv", [["-bg"], ["--background"], ["--workers", "4"], ["--ram", "20"]]
)
def test_removed_flags_are_rejected(argv):
    """``-bg``, ``--workers`` and ``--ram`` are gone, and silence would be wrong.

    Every run is below-normal priority now, so ``-bg`` has nothing left to
    request; ``--cores`` is the only worker-count control; ``--ram`` was removed
    outright, being validated and recorded but never enforced (specification.md
    §11). A removed flag has to fail loudly rather than be ignored — a
    script still passing ``--ram 20`` would otherwise look like it was being
    honoured.
    """
    with pytest.raises(ParamError, match="Unknown argument"):
        parse_args(BASE + argv)


def test_worker_count_is_derived_from_cores_when_not_given():
    """One worker per core: the master is blocked on ``Pool.imap`` throughout.

    It used to reserve a core for the master, which cost a fifth of the boundary
    stage's throughput on a 6-core machine to protect a process that spends over
    99 % of that stage waiting.
    """
    assert parse_args(BASE + ["--cores", "6"]).workers == 6
    assert parse_args(BASE + ["--cores", "1"]).workers == 1


def test_worker_count_is_not_capped():
    """``--cores`` is a budget the user set, so it is honoured exactly.

    An earlier revision capped the derived worker count at 8, which silently
    contradicted an explicit ``--cores`` above that.
    """
    assert parse_args(BASE + ["--cores", "128"]).workers == 128


def test_worker_count_defaults_to_the_machines_logical_cores(monkeypatch):
    monkeypatch.setattr(sysinfo, "logical_core_count", lambda: 12)
    assert parse_args(BASE).workers == 12
    assert default_workers(None) == 12


def test_budget_ranges_are_validated():
    with pytest.raises(ParamError, match="--cores"):
        parse_args(BASE + ["--cores", "0"])


def test_default_output_path_is_derived_from_the_input():
    step, log = resolve_output_paths(os.path.join("dir", "ball.step"), None, 20.0, 4.0)
    assert os.path.basename(step) == "ball-lattice-cc20t4.step"
    assert os.path.basename(log) == "ball-lattice-cc20t4.log"
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


# --- name and path lengths (specification.md §3, -o) -----------------------

PATH_LIMIT = cli.WINDOWS_PATH_MAX if os.name == "nt" else cli.POSIX_PATH_MAX


def folder_of_length(base, n):
    """An absolute folder path exactly ``n`` characters long, below ``base``.

    Nothing is created: the check is pure syntax, so the folders need not
    exist. Segments stay well inside a component's 255.
    """
    folder = os.path.abspath(str(base))
    assert len(folder) + 2 <= n
    while len(folder) < n:
        room = n - len(folder) - 1                 # one separator per segment
        seg = min(room, 200)
        if room - seg in (1,):                     # never leave a gap of one
            seg -= 1
        folder = os.path.join(folder, "d" * seg)
    assert len(folder) == n
    return folder


def test_a_file_name_one_past_the_limit_is_refused_naming_the_excess():
    # The cancel sentinel `<stem>.cancel` is the longest name a run writes.
    stem = "p" * (cli.NAME_MAX - len(cli.LONGEST_SIBLING_SUFFIX) + 1)
    with pytest.raises(ParamError, match="1 character") as err:
        cli.check_path_lengths(stem + ".step")
    assert err.value.fix_text({"input": "input file"}) == (
        "shorten the input file name by 1 character(s)")


def test_a_file_name_at_the_limit_is_not_refused_for_its_name():
    stem = "p" * (cli.NAME_MAX - len(cli.LONGEST_SIBLING_SUFFIX))
    try:
        cli.check_path_lengths(stem + ".step")
    except ParamError as exc:  # the folder can still make the *path* too long
        assert "file name" not in str(exc)


def test_a_folder_too_deep_for_the_temp_files_is_refused(tmp_path):
    # Short stem: the files staged under temp/<stamp>/ are the longest paths.
    n = PATH_LIMIT - cli.TEMP_PATH_RESERVE - 1 + 5
    step = os.path.join(folder_of_length(tmp_path, n), "p.step")
    with pytest.raises(ParamError, match="5 character") as err:
        cli.check_path_lengths(step)
    assert err.value.fix_text({"input": "input file", "output": "output folder"}) == (
        "choose an output folder with a path 5 character(s) shorter")


def test_a_folder_exactly_deep_enough_is_accepted(tmp_path):
    n = PATH_LIMIT - cli.TEMP_PATH_RESERVE - 1
    cli.check_path_lengths(os.path.join(folder_of_length(tmp_path, n), "p.step"))


def test_a_long_input_name_making_the_path_too_long_names_the_input(tmp_path):
    stem = "p" * 100        # longer than the temp reserve, so it decides
    name_len = len(stem) + len(cli.LONGEST_SIBLING_SUFFIX)
    n = PATH_LIMIT - name_len - 1 + 3
    step = os.path.join(folder_of_length(tmp_path, n), stem + ".step")
    with pytest.raises(ParamError, match="3 character") as err:
        cli.check_path_lengths(step)
    assert err.value.fix_text({"input": "input file", "output": "output folder"}) == (
        "choose an output folder with a path 3 character(s) shorter, "
        "or shorten the input file name")


def test_path_length_is_counted_the_way_the_platform_counts_it():
    # One UTF-16 unit on Windows; two UTF-8 bytes on Linux.
    assert cli._path_units("é") == (1 if os.name == "nt" else 2)


def test_preflight_refuses_a_path_too_long_before_touching_the_disk(tmp_path):
    n = PATH_LIMIT - cli.TEMP_PATH_RESERVE + 10
    out = os.path.join(folder_of_length(tmp_path, n), "p.step")
    args = parse_args(["-i", "in.step", "-cc", "10", "-t", "1.5", "-o", out])
    # The input does not exist either: the length is refused first, as a
    # parameter (exit 2), not as a missing input.
    with pytest.raises(ParamError, match="longer than"):
        preflight_checks(args)


# --- print orientation and support (specification.md §4.5, §4.6) -----------


def _args(*extra):
    return parse_args(["-i", "part.step", "-cc", "10", "-t", "1.5", *extra])


def _args_at(cc, t, *extra):
    return parse_args(["-i", "part.step", "-cc", cc, "-t", t, *extra])


def test_orientation_and_support_default_to_off():
    args = _args()
    assert args.orient == (0.0, 0.0, 0.0)
    assert args.support is False
    assert os.path.basename(args.output) == "part-lattice-cc10t1.5.step"
    assert "orient" not in args.as_dict() and "support" not in args.as_dict()


def test_orientation_is_accepted_without_support_and_named():
    args = _args("--orient", "30", "0", "-22.5")
    assert args.orient == (30.0, 0.0, -22.5)
    assert args.support is False
    assert os.path.basename(args.output) == "part-lattice-cc10t1.5-rotx30y0z-22.5.step"
    assert args.as_dict()["orient"] == "30 0 -22.5 deg"


def test_support_defaults_and_names():
    args = _args("--support")
    # A fixed 1 mm unless set, whatever the strut size.
    assert args.support and args.support_thickness == 1.0
    assert _args_at("20", "4", "--support").support_thickness == 1.0
    assert os.path.basename(args.output) == "part-lattice-cc10t1.5-sup1.step"
    assert args.as_dict()["support_thickness"] == "1 mm"
    both = _args("--orient", "30", "0", "0", "--support", "--support-thickness", "2")
    assert both.support_thickness == 2.0
    assert os.path.basename(both.output) == "part-lattice-cc10t1.5-rotx30y0z0-sup2.step"
    assert os.path.basename(both.log_path) == "part-lattice-cc10t1.5-rotx30y0z0-sup2.log"
    thin = _args("--support", "--support-thickness", "0.8")
    assert os.path.basename(thin.output) == "part-lattice-cc10t1.5-sup0.8.step"


@pytest.mark.parametrize("extra", [
    ("--support-thickness", "2"),
    ("--orient", "10", "20"),
    ("--orient", "400", "0", "0"),
    ("--support", "--support-thickness", "0.1"),
])
def test_bad_orientation_and_support_arguments_are_rejected(extra):
    with pytest.raises(ParamError):
        _args(*extra)


@pytest.mark.parametrize("cc, t, thickness, m, thinnest", [
    (10.0, 1.5, 1.5, 3, 1.607),
    (7.0, 1.4, 1.0, 2, 1.775),
    (20.0, 4.0, 4.0, 2, 5.071),
    (5.0, 1.0, 1.0, 2, 1.268),
    (10.0, 1.5, 0.2, 7, 0.260),
])
def test_the_support_grid_is_the_finest_whose_thinnest_member_holds(cc, t, thickness, m, thinnest):
    from latticegen2.cli import support_grid

    a = cc / math.sqrt(2.0)
    got_m, got = support_grid(a, t, thickness)
    assert (got_m, round(got, 3)) == (m, thinnest)
    assert got >= thickness
    # One division finer would be too thin.
    assert a / (m + 1) - t / 2 < thickness


def test_a_support_thickness_beyond_a_minus_t_is_refused_naming_the_largest():
    from latticegen2.cli import support_grid

    # cc 5, t 1: a - t = 2.5355 mm, delivered at m = 1; the message names it
    # rounded down, so the figure it names is itself accepted.
    a = 5.0 / math.sqrt(2.0)
    assert support_grid(a, 1.0, a - 1.0) == (1, pytest.approx(a - 1.0))
    args = _args_at("5", "1", "--support", "--support-thickness", "2.535")
    assert args.support_thickness == 2.535
    with pytest.raises(ParamError, match=r"at most a - t = 2\.535 mm"):
        _args_at("5", "1", "--support", "--support-thickness", "2.536")
    # Asked only of a supported run: without --support the default reaches
    # nothing, so a legal lattice whose a - t is under 1 mm is not refused.
    assert _args_at("5", "3").support is False
    with pytest.raises(ParamError, match="at most a - t"):
        _args_at("5", "3", "--support")


FLAG_LABELS = {"cc": "-cc", "t": "-t", "cores": "--cores", "thickness": "--support-thickness",
               "rx": "rx", "ry": "ry", "rz": "rz", "orient": "--orient", "input": "input"}


def _fix(*argv) -> str:
    with pytest.raises(ParamError) as info:
        cli.parse_args(["-i", "part.step", *argv])
    return info.value.fix_text(FLAG_LABELS)


@pytest.mark.parametrize("argv, flags", [
    # A strut too thick for its cell: either change is enough on its own.
    (["-cc", "10", "-t", "18"], ("-t", "-cc")),
    (["-cc", "5", "-t", "3.6"], ("-t", "-cc")),
    # A support thickness no grid can hold.
    (["-cc", "10", "-t", "1.5", "--support", "--support-thickness", "6"],
     ("--support-thickness", "-cc")),
    (["-cc", "5", "-t", "3", "--support"], ("--support-thickness", "-cc")),
    # Neither a thinner support nor a larger cell is in range: a thinner strut.
    (["-cc", "28.43", "-t", "20", "--support", "--support-thickness", "19"], ("-t",)),
])
def test_every_change_a_refusal_suggests_is_accepted(argv, flags):
    """The window shows the fix while Start! is greyed out (specification.md
    §3.1), so each change it offers must be enough on its own -- bounds are
    rounded in the direction that makes the named value itself pass."""
    fix = _fix(*argv)
    found = re.findall(r"(-cc|-t|--support-thickness) (?:smaller|larger) than ([\d.]+)"
                       r"|(-cc|-t|--support-thickness) to at (?:most|least) ([\d.]+)", fix)
    changes = [(a or c, b or d) for a, b, c, d in found]
    assert tuple(f for f, _ in changes) == flags, fix
    for flag, value in changes:
        changed = list(argv)
        if flag in changed:
            changed[changed.index(flag) + 1] = value
        else:
            changed += [flag, value]
        cli.parse_args(["-i", "part.step", *changed])


def test_a_refusal_names_the_parameter_and_its_range():
    assert _fix("-cc", "100", "-t", "1") == "set -cc between 0.4 and 50 mm"
    assert _fix("-cc", "abc", "-t", "1") == "enter a number for -cc"
    assert _fix("-cc", "10", "-t", "1", "--cores", "1.5") == "enter a whole number for --cores"
    assert _fix("-cc", "10", "-t", "1", "--orient", "0", "400", "0") == "set ry between -360 and 360°"
    with pytest.raises(ParamError) as info:
        cli.parse_args(["-cc", "10", "-t", "1"])
    assert info.value.fix_text(FLAG_LABELS) == "choose an input"
    # No single change in range is enough: all three are named, none bounded.
    assert (_fix("-cc", "28.43", "-t", "20", "--support", "--support-thickness", "20")
            == "lower --support-thickness or -t, or raise -cc")


def test_a_fix_naming_a_field_without_a_label_keeps_the_field_name():
    assert ParamError("x", fix="raise {cc} or {budget}").fix_text({"cc": "cc"}) == (
        "raise cc or budget")


@pytest.mark.parametrize("extra", [("--overhang", "60"),
                                   ("--support", "--overhang", "45")])
def test_the_removed_overhang_flag_is_refused_by_name(extra):
    """Not "unknown argument": a script written against the skeleton is told why."""
    with pytest.raises(ParamError, match="--overhang has been removed"):
        _args(*extra)


def test_the_support_angles_are_the_rod_angle_and_sixty_degrees():
    """Fixed, not parameters (specification.md §4.6); the rod angle is the
    expression, not a typed decimal."""
    import math

    from latticegen2.cli import SUPPORT_OVERHANG_DEG, SUPPORT_START_DEG

    assert SUPPORT_OVERHANG_DEG == math.degrees(math.asin(math.sqrt(2.0 / 3.0)))
    assert SUPPORT_OVERHANG_DEG == pytest.approx(54.7356, abs=1e-4)
    assert SUPPORT_START_DEG == 60.0
