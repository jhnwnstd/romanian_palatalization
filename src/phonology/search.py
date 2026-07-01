"""SEARCH primitive: one mechanism, terminator-breadth as its parameter.

The paper (Dabbous et al. 2024) frames phonological locality as a
search procedure with three parameters: direction, terminator class,
and condition on the terminator. The terminator class is the
discriminating bit:

  - Empty / wildcard terminator (``terminator = {}``) is **broad**:
    the scan halts on the *first* following segment. The condition is
    then applied as a yes/no filter on that segment. Anything that
    isn't a member of the condition class blocks the rule. This is
    dorsal palatalization (latex.tex e:velar-pal): a /T/ between /K/
    and /e/ stops the scan and the dorsal stays unpalatalized.

  - Non-empty terminator is **narrow**: the scan steps past any
    segment that isn't a member of the terminator class. Non-members
    are transparent. This is S-palatalization and assibilation
    (latex.tex e:s-pal-rev, e:assib): /S/ in ``prost`` palatalizes
    across the stop because the stop isn't in the postalveolar place
    class so the scan keeps going.

  - Leftward + broad = strict left-adjacency. This is the bleed rule
    (latex.tex e:bleed-rev): the postalveolar must sit immediately
    before the target. The same Search class expresses it without a
    separate "adjacent" knob.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Final

from .segments import FeatureMap, Word


class Direction(StrEnum):
    LEFT = "L"
    RIGHT = "R"


class SearchOutcome(StrEnum):
    """Why a :meth:`Search.locate` call returned the value it did.

    Used by :mod:`phonology.diagnostics.explain` to narrate rule
    behaviour. Not needed for the fast :meth:`locate` path.
    """

    LICENSED = "LICENSED"  # trigger found and condition holds
    BROAD_TERMINATOR_BLOCK = (
        "BROAD_TERMINATOR_BLOCK"  # first seg failed condition
    )
    NARROW_SCAN_EXHAUSTED = "NARROW_SCAN_EXHAUSTED"  # scanned to word edge
    EMPTY_SCAN = "EMPTY_SCAN"  # target at edge; no segments to look at


@dataclass(frozen=True, slots=True)
class LocateReport:
    """Detailed outcome of a :meth:`Search.locate_verbose` call."""

    trigger_index: int | None
    outcome: SearchOutcome
    inspected_index: int | None = None  # segment the scan halted on


@dataclass(frozen=True, slots=True)
class Search:
    """Locate a trigger segment for a rule with target at ``anchor``.

    ``locate`` returns the index of the matched segment, or ``None`` if
    no match. The target segment itself (``word[anchor]``) is skipped
    — searches always look at *other* segments.
    """

    direction: Direction
    terminator: FeatureMap  # ``{}`` for broad next-segment
    condition: FeatureMap  # what the terminator must satisfy

    def locate(self, word: Word, anchor: int) -> int | None:
        step = 1 if self.direction is Direction.RIGHT else -1
        i = anchor + step
        n = len(word)
        while 0 <= i < n:
            seg = word[i]
            if seg.matches(self.terminator):
                if seg.matches(self.condition):
                    return i
                return None  # broad-terminator block
            i += step
        return None

    def locate_verbose(self, word: Word, anchor: int) -> LocateReport:
        """Same as :meth:`locate` but returns a :class:`LocateReport`
        naming the outcome and (when relevant) the segment that
        halted the scan.

        Explain-mode uses this to say why a rule didn't fire —
        distinguishing "broad terminator halted on non-matching X"
        from "narrow terminator scanned to word edge without finding
        anything".
        """
        step = 1 if self.direction is Direction.RIGHT else -1
        i = anchor + step
        n = len(word)
        if not (0 <= i < n):
            return LocateReport(None, SearchOutcome.EMPTY_SCAN)
        while 0 <= i < n:
            seg = word[i]
            if seg.matches(self.terminator):
                if seg.matches(self.condition):
                    return LocateReport(i, SearchOutcome.LICENSED, i)
                return LocateReport(
                    None,
                    SearchOutcome.BROAD_TERMINATOR_BLOCK,
                    i,
                )
            i += step
        return LocateReport(None, SearchOutcome.NARROW_SCAN_EXHAUSTED)


__all__: Final[tuple[str, ...]] = (
    "Direction",
    "LocateReport",
    "Search",
    "SearchOutcome",
)
