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

:class:`LexRow` is the row contract every analysis shares. The
Romanian lexicon's specifics also live here, as the defaults of two
:class:`~phonology.profile.AnalysisProfile` hooks:
:meth:`LexRow.from_dict` reads the Romanian CSV schema and tags
clusters with the Romanian orthographic tagger (``row_factory``), and
:func:`romanian_cluster_key` buckets rows into the paper's cluster
cells (``cluster_key``). A lexicon with another schema passes its own
``row_factory`` to :func:`iter_lexicon_rows`; rows from such a factory
carry ``ClusterTag.NONE`` unless the analysis has a use for the tag.
Analysis constants (``TARGET_STEMS``, etc.) live with the analysis
at :mod:`phonology.analyses.romanian_palatalization` and are passed
in.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Final, Iterator, Mapping, Protocol


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
    def from_dict(cls, r: Mapping[str, str]) -> LexRow:
        """Build a LexRow from a raw csv.DictReader row.

        Reads the Romanian lexicon's column names; the default
        ``row_factory``. Missing fields become empty strings;
        ``mutation`` is parsed case-insensitively as True iff the
        value is 'TRUE'.
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


class RowFactoryFn(Protocol):
    """Turn one ``csv.DictReader`` row into a :class:`LexRow`.

    One CSV row yields one LexRow, so a lexicon with several cells
    per lexeme is written in long format (one row per cell).
    """

    def __call__(self, r: Mapping[str, str], /) -> LexRow: ...


# Report buckets for the paper's cluster cells: the validator's
# "Cluster-specific verification" table and the mismatch CSV's
# ``cluster`` column. Other (cluster, opportunity) pairs get no bucket.
_ROMANIAN_CLUSTER_KEYS: Final[Mapping[tuple[ClusterTag, str], str]] = {
    (ClusterTag.SC, "e"): "sc(ă)+e",
    (ClusterTag.SC, "i"): "sc(ă)+i",
    (ClusterTag.CT, "e"): "ct(ă)+e",
    (ClusterTag.CT, "i"): "ct(ă)+i",
    (ClusterTag.ST, "i"): "st(ă)+i",
}


def romanian_cluster_key(row: LexRow) -> str | None:
    """The row's cluster-report bucket, or ``None`` if it has none.

    The default :attr:`~phonology.profile.AnalysisProfile.cluster_key`.
    """
    return _ROMANIAN_CLUSTER_KEYS.get((row.cluster, row.opportunity))


def iter_lexicon_rows(
    csv_path: Path,
    *,
    target_stems: frozenset[str],
    require_plural: bool = True,
    skip_opportunity_none: bool = True,
    skip_nde_exceptions: bool = True,
    pos_filter: str | None = "N",
    limit: int | None = None,
    row_factory: RowFactoryFn = LexRow.from_dict,
) -> Iterator[LexRow]:
    """Stream a lexicon CSV as :class:`LexRow` under the given filters.

    Every filter is a keyword arg so callers can subset without
    editing shared code. The default filter chain matches what the
    palatalization validator has always used: nouns only, target
    stem-final consonant present, non-empty plural, opportunity !=
    'none', NDE exceptions excluded. Toggle knobs to relax any of
    these for other analyses.

    ``row_factory`` builds each LexRow (default: the Romanian schema;
    pass ``profile.row_factory``). The filters test LexRow fields, so
    a factory for another schema maps into the same vocabulary (pos
    ``"N"``, opportunity ``"none"``, ``"nde:"`` exception prefix) or
    the caller relaxes the filters.
    """
    yielded = 0
    with csv_path.open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            row = row_factory(r)
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
    "RowFactoryFn",
    "iter_lexicon_rows",
    "romanian_cluster_key",
)
