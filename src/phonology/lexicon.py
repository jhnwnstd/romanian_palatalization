"""Lexicon layer: streaming row iterator + neutral row dataclass.

Historically the driver and ``scripts/build_contingency_table.py``
each carried their own CSV filter chain, their own dataclass, and
their own cluster-detection helpers — three sources of truth for the
same schema that drifted (validator's ``TARGET_STEMS`` included /z/,
contingency's did not). This module is the single answer:

  - :class:`LexRow` — one immutable row shape everyone consumes.
  - :class:`ClusterTag` — orthographic cluster classification computed
    ONCE at row construction, replacing the four ``has_*_cluster``
    properties that were duplicated across scripts.
  - :func:`iter_lexicon_rows` — the single streaming iterator. Every
    consumer of the lexicon should call this, not open the CSV
    directly.

Neither this module nor its consumers should hold Romanian-specific
constants (``TARGET_STEMS``, etc.) — those live with the analysis at
:mod:`phonology.analyses.romanian_palatalization` and are passed in.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Final, Iterator


class ClusterTag(StrEnum):
    """Orthographic cluster classification of a lemma.

    Computed once from ``lemma.lower()`` at :meth:`LexRow.from_dict`.
    Values roughly follow the paper's stated cluster categories.
    """

    NONE = "NONE"
    CT = "CT"  # -ct, -ctă  (Latinate cluster, dorsal blocked)
    SC = "SC"  # -sc, -scă  (palatalises via /sk/ → [ʃt])
    SHCA = "SHCA"  # -șcă, -şcă (already-palatalized /ʃ/ in the sg)
    ST = "ST"  # -st, -stă, -ște, -şte (S-pal + bleed)


def _tag_cluster(lemma: str) -> ClusterTag:
    """Classify the lemma's trailing consonant cluster.

    Order matters: SHCA must be checked before SC (since -șcă
    endswith 'că' but is a different phonological beast — its
    singular already carries surface /ʃ/), and ST-family suffixes
    include the palatalized -ște/-şte forms that end in 'e'.
    """
    lem = lemma.lower()
    if lem.endswith(("șcă", "şcă")):
        return ClusterTag.SHCA
    if lem.endswith(("ct", "ctă")):
        return ClusterTag.CT
    if lem.endswith(("sc", "scă")):
        return ClusterTag.SC
    if lem.endswith(("st", "stă", "ște", "şte")):
        return ClusterTag.ST
    return ClusterTag.NONE


@dataclass(frozen=True, slots=True)
class LexRow:
    """One lexicon CSV row, narrowed to the fields all consumers need.

    Immutable. ``cluster`` is derived from ``lemma`` at construction
    time; every reader that used to call a ``has_*_cluster`` property
    now just checks ``row.cluster is ClusterTag.X``.
    """

    lemma: str
    plural: str
    pos: str
    stem_final: str
    opportunity: str
    mutation: bool
    ipa_lemma: str
    ipa_pl: str
    exception_reason: str
    cluster: ClusterTag

    @classmethod
    def from_dict(cls, r: dict[str, str]) -> LexRow:
        """Build a LexRow from a raw csv.DictReader row.

        Missing fields become empty strings; ``mutation`` is parsed
        case-insensitively as True iff the value is 'TRUE'.
        """
        lemma = (r.get("lemma") or "").strip()
        return cls(
            lemma=lemma,
            plural=(r.get("plural") or "").strip(),
            pos=(r.get("pos") or "").strip(),
            stem_final=(r.get("stem_final") or "").strip(),
            opportunity=(r.get("opportunity") or "").strip(),
            mutation=(r.get("mutation") or "").strip().upper() == "TRUE",
            ipa_lemma=(r.get("ipa_normalized_lemma") or "").strip(),
            ipa_pl=(r.get("ipa_normalized_pl") or "").strip(),
            exception_reason=(r.get("exception_reason") or "").strip(),
            cluster=_tag_cluster(lemma),
        )


def iter_lexicon_rows(
    csv_path: Path,
    *,
    target_stems: frozenset[str],
    require_plural: bool = True,
    skip_opportunity_none: bool = True,
    skip_nde_exceptions: bool = True,
    pos_filter: str | None = "N",
    limit: int | None = None,
) -> Iterator[LexRow]:
    """Stream a lexicon CSV as :class:`LexRow` under the given filters.

    Every filter is a keyword arg so callers can subset without
    editing shared code. The default filter chain matches what the
    palatalization validator has always used: nouns only, target
    stem-final consonant present, non-empty plural, opportunity !=
    'none', NDE exceptions excluded. Toggle knobs to relax any of
    these for other analyses.
    """
    yielded = 0
    with csv_path.open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            row = LexRow.from_dict(r)
            if pos_filter is not None and row.pos != pos_filter:
                continue
            if row.stem_final not in target_stems:
                continue
            if require_plural and not row.plural:
                continue
            if skip_opportunity_none and row.opportunity == "none":
                continue
            if skip_nde_exceptions and row.exception_reason.startswith("nde:"):
                continue
            yield row
            yielded += 1
            if limit is not None and yielded >= limit:
                return


__all__: Final[tuple[str, ...]] = (
    "ClusterTag",
    "LexRow",
    "iter_lexicon_rows",
)
