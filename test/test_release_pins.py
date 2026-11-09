"""The bundle's pins and its licence table agree (licenses/LICENSES.md).

Every wheel a release bundle holds is pinned in requirements-bundle.txt, which
tools/build_release.py enforces at build time. This holds each pin to a row of
the licence table at the same version, with every licence file it links to
present, so a pin cannot be added or moved without the table following.
"""

import os
import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))

import build_release  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
LICENSES = ROOT / "licenses"


def _table_rows() -> list[list[str]]:
    rows = []
    for line in (LICENSES / "LICENSES.md").read_text(encoding="utf-8").splitlines():
        if line.startswith("| **"):
            rows.append([c.strip() for c in line.strip().strip("|").split("|")])
    return rows


def _row_for(name: str, rows: list[list[str]]) -> list[str] | None:
    """The row naming ``name`` as ``**name**`` or as ``(`name`)``, exactly."""
    want = build_release._canonical(name)
    for row in rows:
        names = re.findall(r"\*\*([^*]+)\*\*|`([^`]+)`", row[0])
        if any(build_release._canonical(a or b) == want for a, b in names):
            return row
    return None


PINS = build_release.pinned_requirements()


@pytest.mark.parametrize("name", sorted(PINS))
def test_every_pin_has_a_licence_row_at_its_version(name):
    row = _row_for(name, _table_rows())
    assert row is not None, f"{name} is pinned but has no row in LICENSES.md"
    assert row[1] == PINS[name], (
        f"{name}: pinned {PINS[name]}, LICENSES.md records {row[1]}")
    files = re.findall(r"\]\(([^)]+)\)", row[-1])
    assert files, f"{name}: its row links no licence text"
    for f in files:
        path = LICENSES / f
        assert path.is_file() or (path.is_dir() and any(path.rglob("*"))), \
            f"{name}: {f} is missing"


def _linked_files() -> set[Path]:
    text = (LICENSES / "LICENSES.md").read_text(encoding="utf-8")
    out = set()
    for target in re.findall(r"\]\(([^)#]+)\)", text):
        if "://" in target or target.startswith(".."):
            continue
        path = (LICENSES / target).resolve()
        out |= set(path.rglob("*")) if path.is_dir() else {path}
    return {p for p in out if p.is_file()}


def test_every_licence_link_resolves():
    text = (LICENSES / "LICENSES.md").read_text(encoding="utf-8")
    for target in re.findall(r"\]\(([^)#]+)\)", text):
        if "://" in target:
            continue
        assert (LICENSES / target).exists(), f"LICENSES.md links {target}, which is missing"


def test_every_licence_text_is_linked_from_the_table():
    """A text no row links is one nobody can tell the purpose of."""
    texts = {p.resolve() for p in LICENSES.rglob("*")
             if p.is_file() and p.name != "LICENSES.md" and "runtime" not in p.parts}
    orphans = sorted(str(p.relative_to(LICENSES.resolve())) for p in texts - _linked_files())
    assert not orphans, f"licence texts linked from nowhere: {orphans}"


def test_vtk_and_what_it_pulls_in_are_pinned():
    """The ten packages vtk 9.6.2 brings in were once resolved unpinned."""
    for name in ("vtk", "matplotlib", "contourpy", "cycler", "fonttools",
                 "kiwisolver", "packaging", "pillow", "pyparsing",
                 "python-dateutil", "six"):
        assert name in PINS


def test_an_unpinned_or_mismatched_wheel_is_refused():
    pins = {"numpy": "2.4.6", "python-dateutil": "2.9.0.post0"}
    wheels = [Path("numpy-2.4.6-cp311-cp311-win_amd64.whl"),
              Path("python_dateutil-2.9.0.post0-py2.py3-none-any.whl"),
              Path("numpy-2.4.5-cp311-cp311-win_amd64.whl"),
              Path("six-1.17.0-py2.py3-none-any.whl")]
    assert build_release.unpinned_wheels(wheels, pins) == [
        "numpy-2.4.5-cp311-cp311-win_amd64.whl",
        "six-1.17.0-py2.py3-none-any.whl"]
