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
        assert (LICENSES / f).is_file(), f"{name}: {f} is missing"


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
