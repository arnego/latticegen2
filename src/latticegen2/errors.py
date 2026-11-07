"""Typed failure modes and their process exit codes.

Every non-zero exit prints exactly one human-readable reason line naming what
failed and why (specification.md §7). The full table is in docs/algorithm.md §10.
"""

from __future__ import annotations


class LatticeGenError(Exception):
    """Base class for every failure this tool reports with a specific exit code."""

    exit_code = 1


class ParamError(LatticeGenError):
    """Invalid, missing, out-of-range or conflicting parameter (exit 2).

    Raised before any computation starts (specification.md §7).

    ``fix`` says what to change for the value to be accepted, as a template
    whose ``{cc}``, ``{t}``, ``{cores}``, ``{thickness}``, ``{rx}``, ``{ry}``,
    ``{rz}``, ``{orient}`` and ``{input}`` fields name the parameters; the
    window fills them with its own labels (specification.md §3.1). ``None``
    where no field of the window can cause the refusal.
    """

    exit_code = 2

    def __init__(self, message: str, fix: str | None = None):
        super().__init__(message)
        self.fix = fix

    def fix_text(self, labels: dict[str, str]) -> str | None:
        """``fix`` with each parameter named by ``labels``."""
        return self.fix.format(**labels) if self.fix else None


class InputGeometryError(LatticeGenError):
    """The input STEP could not be read, parsed, or faithfully meshed (exit 3)."""

    exit_code = 3


class ProcessingError(LatticeGenError):
    """Geometry processing failed (exit 4).

    Covers the boundary trim, the shell build, and the connectivity gate's
    refusal to guess about an unresolved fragment (docs/algorithm.md §8).
    """

    exit_code = 4


class ResourceError(LatticeGenError):
    """A resource limit could not be respected — memory, disk (exit 5)."""

    exit_code = 5


class OutputError(LatticeGenError):
    """The output file could not be written (exit 6)."""

    exit_code = 6


class CancelledError(LatticeGenError):
    """The user cancelled the run with Ctrl+C (exit 130 = 128 + SIGINT)."""

    exit_code = 130
