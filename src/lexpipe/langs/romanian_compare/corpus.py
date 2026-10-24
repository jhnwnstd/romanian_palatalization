"""The Romanian Leipzig word-frequency table as a core attestation table.

``data/leipzig/freq/ron_news_2024_1M_freq.csv`` (``word,freq``,
lowercased by ``scripts/download_leipzig.py``; read-only, R2) is the
corpus of the legacy ``freq_ron_news_2024_1M`` column. Keys are folded
as ``scripts/attach_frequencies.py`` folds them (NFC, cedilla ş/ţ →
comma-below ș/ț, lowercase; counts of keys that fold together are
summed), which reproduces that column for every lexicon row (a test
checks it). Without the cedilla fold 1,617 rows would differ.

The table is lowercased, so Romanian has no case-sensitive counts:
:class:`lexpipe.attest.AttestationTable` gets the same mapping for both
of its tables, and ``total_tokens`` is the sum of the table's counts
(the cleaned word tokens of ron_news_2024_1M, 18,628,329). Thresholds
τ are then decided by the core code (:data:`lexpipe.attest.TAUS`).
"""

from __future__ import annotations

import csv
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

from lexpipe.attest import TAUS, AttestationTable
from lexpipe.langs.romanian_compare.lexicon import fold_key


@dataclass(frozen=True, slots=True)
class RomanianCounts:
    """Folded corpus counts and the corpus size in word tokens."""

    counts: Mapping[str, int]
    total_tokens: int

    def count(self, form: str) -> int:
        """Tokens of *form* (0 if absent or empty)."""
        key = fold_key(form)
        return self.counts.get(key, 0) if key else 0


def load_counts(path: Path) -> RomanianCounts:
    """Read the ``word,freq`` table (one streaming pass)."""
    counts: dict[str, int] = {}
    total = 0
    with path.open("r", encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            key = fold_key(row.get("word") or "")
            try:
                freq = int((row.get("freq") or "").strip())
            except ValueError:
                continue
            if not key or freq <= 0:
                continue
            counts[key] = counts.get(key, 0) + freq
            total += freq
    if total <= 0:
        raise ValueError(f"{path}: no usable rows")
    return RomanianCounts(counts, total)


def attestation_table(
    counts: Mapping[str, int], total_tokens: int
) -> AttestationTable:
    """A core :class:`AttestationTable` over *counts* (folded keys)."""
    return AttestationTable(
        counts_cs=counts, counts_folded=counts, total_tokens=total_tokens
    )


def attested_by_tau(
    table: AttestationTable, form: str, *, in_dictionary: bool = True
) -> dict[str, bool]:
    """``{tau: attested}`` for every τ of :data:`lexpipe.attest.TAUS`.

    A form that is not a dictionary form (a corpus-only candidate) is
    never attested at ``dictionary``.
    """
    result = table.attested_by_tau(fold_key(form), TAUS)
    if not in_dictionary:
        result["dictionary"] = False
    return result


def subset_counts(
    counts: RomanianCounts, forms: Iterable[str]
) -> dict[str, int]:
    """Counts of *forms* only (keys folded, zero counts kept), sorted."""
    keys = sorted({fold_key(form) for form in forms if fold_key(form)})
    return {key: counts.counts.get(key, 0) for key in keys}
