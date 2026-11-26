"""What the bad-geometry viewer marks, and nothing about how (specification.md §3.1).

A run that warns, or fails, at a place in the part sends a ``mark`` event
(:mod:`latticegen2.progress`) carrying that place as numbers. This module turns
such an event into :class:`Mark` objects and answers the questions the viewer
asks of them -- where to centre, which preview face a failure names, what the
detail box says -- without a display, like :mod:`latticegen2.gui.model`.

The event arrives from another process, so it is checked rather than trusted:
a malformed place is skipped, and a mark left with none is no mark at all.
"""

from __future__ import annotations

import textwrap
from dataclasses import dataclass

import numpy as np

from .. import progress
from ..runlog import WARNING_PREFIX

#: Marks of a failure are red, of a warning yellow.
ERROR_COLOUR = (0.86, 0.12, 0.12)
WARNING_COLOUR = (0.98, 0.80, 0.10)

#: How many coordinates each place kind carries, in the order the event names them.
_KEYS = {"point": ("p",), "segment": ("a", "b"), "box": ("lo", "hi"),
         "face": ("centroid", "lo", "hi")}

#: How long a failure's message may run in a detail box before it is cut.
DETAIL_CHARS = 700


@dataclass(frozen=True)
class MarkPlace:
    """One place a message points at: its line of text and its points."""

    kind: str
    text: str
    #: ``(k, 3)``: the point; a segment's two ends; a box's ``lo`` and ``hi``;
    #: a face's centroid, ``lo`` and ``hi``.
    points: np.ndarray

    def anchor(self) -> np.ndarray:
        """The one point the view centres on and the pin is drawn at."""
        if self.kind == "point" or self.kind == "face":
            return self.points[0]
        return self.points.mean(axis=0)

    def box(self) -> "tuple[np.ndarray, np.ndarray] | None":
        if self.kind == "box":
            return self.points[0], self.points[1]
        if self.kind == "face":
            return self.points[1], self.points[2]
        return None


@dataclass(frozen=True)
class Mark:
    severity: str
    title: str
    message: str
    places: tuple

    @property
    def is_error(self) -> bool:
        return self.severity == progress.ERROR_MARK

    @property
    def colour(self) -> tuple:
        return ERROR_COLOUR if self.is_error else WARNING_COLOUR

    def headline(self) -> str:
        """``Warning: <title>`` or ``Error: <title>``, for the issue list."""
        return ("Error: " if self.is_error else "Warning: ") + (self.title or "the run failed")

    def sentence(self) -> str:
        """What went wrong and what it may cause, without the places.

        A warning's first line, less its prefix, title and the closing
        "Inspect, in the part's coordinates:"; a failure's whole reason, cut at
        :data:`DETAIL_CHARS`.
        """
        if self.is_error:
            text = " ".join(self.message.split())
            return text if len(text) <= DETAIL_CHARS else text[:DETAIL_CHARS - 1] + "…"
        first = self.message.split("\n", 1)[0]
        if first.startswith(WARNING_PREFIX):
            first = first[len(WARNING_PREFIX):]
        if self.title and first.startswith(self.title + ": "):
            first = first[len(self.title) + 2:]
        suffix = " Inspect, in the part's coordinates:"
        if first.endswith(suffix):
            first = first[:-len(suffix)]
        return first


def from_event(event: dict) -> "Mark | None":
    """The mark a ``mark`` event describes, or ``None`` if it describes none."""
    severity = event.get("severity")
    if severity not in progress.SEVERITIES:
        return None
    places = []
    for raw in event.get("places") or ():
        place = _place(raw)
        if place is not None:
            places.append(place)
    if not places:
        return None
    return Mark(severity=severity, title=str(event.get("title") or ""),
                message=str(event.get("message") or ""), places=tuple(places))


def _place(raw) -> "MarkPlace | None":
    if not isinstance(raw, dict):
        return None
    kind = raw.get("kind")
    keys = _KEYS.get(kind)
    if keys is None:
        return None
    try:
        points = np.array([[float(x) for x in raw[k]] for k in keys], dtype=float)
    except (KeyError, TypeError, ValueError):
        return None
    if points.shape != (len(keys), 3) or not np.all(np.isfinite(points)):
        return None
    return MarkPlace(kind=kind, text=str(raw.get("text") or ""), points=points)


def detail_text(mark: Mark, place: MarkPlace, width: int = 52) -> str:
    """What the detail box beside a clicked mark says."""
    rows = [mark.headline(), ""]
    rows += textwrap.wrap(place.text, width) or [""]
    rows.append("")
    rows += textwrap.wrap(mark.sentence(), width)
    return "\n".join(rows)


def match_face(lo, hi, face_lo: np.ndarray, face_hi: np.ndarray) -> "int | None":
    """The preview face a failure's face is, by bounding box.

    By box rather than by index: the pipeline repairs the input before it
    names a face (docs/algorithm.md §7.5), and a face's place in the traversal
    order is not something either side promises. Two boxes of one face can
    still differ a little -- OCCT tightens a face's box to its triangulation
    once it has one -- so the nearest box is taken, if it is within a tenth of
    the named box's diagonal (and a tenth of a millimetre).
    """
    if face_lo is None or face_hi is None or len(face_lo) == 0:
        return None
    lo, hi = np.asarray(lo, dtype=float), np.asarray(hi, dtype=float)
    face_lo, face_hi = np.asarray(face_lo, dtype=float), np.asarray(face_hi, dtype=float)
    gap = np.maximum(np.abs(face_lo - lo).max(axis=1), np.abs(face_hi - hi).max(axis=1))
    best = int(np.argmin(gap))
    bar = 0.1 * float(np.linalg.norm(hi - lo)) + 0.1
    return best if gap[best] <= bar else None
