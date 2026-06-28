#!/usr/bin/env python3
"""Build the inflection-dependence contingency table.

Canonical script for the table the abstract reports. Replicates the
Steriade (2008) §10.2.3.5 lexical counts under the appropriate coding
for our data:

ROW
    ``plural alternates`` yes/no, derived from the ``mutation`` field
    of the lexicon CSV (which is the direct singular-vs-plural surface
    comparison: t→ț, d→z, s→ș for coronals; orthographic ci/ce/gi/ge
    for dorsals).

COLUMN
    Verbalizer allomorph (``-i`` vs ``-a``). This is Steriade's affix-
    avoidance tabulation in her (10.20). We do not filter the column by
    whether the verb's conjugated paradigm independently surfaces the
    alternation — that is a separate (manifestation b) question, and
    Steriade is explicit that it does not arise with derived verbs.

EXCLUSIONS
    - Rows whose ``exception_reason`` begins with ``nde:`` (the three
      non-derived-environment-blocking classes Steriade also excludes).
    - Derived verbs that don't pass the DEX etymology audit. See
      :mod:`romanian_palatalization.dex_etymology`.
    - Explicit USER_KEEP allowlist for back-formations the etym audit
      mis-rejects: clint, dinte, vargă, verigă.

The output is two panels — TOP-1000 most-frequent noun lemmas and FULL
DEX (no frequency filter) — each split into dorsal (c, g), coronal
(t, d, s), and pooled subpanels.

Important methodological note
-----------------------------
The row coding here is ``mutation``-based, NOT the historic
``opportunity``-based allomorph proxy ``-i/-e`` vs ``-uri``. For
coronals, the surface allomorph "-e" does not trigger assibilation
(e.g. *casă/case* does not alternate), so coding by allomorph
mis-classifies coronal -e nouns as alternating. The ``mutation`` field
records the actual sg-vs-pl comparison, which is what Steriade's K/TS
vs K/K split is about.
"""

from __future__ import annotations

import csv
import sys
from collections import Counter
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Callable, Final

from scipy.stats import fisher_exact

# Add src to path so we can import the shared etymology + cache helpers.
_PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PROJECT_ROOT / "src"))

from dex_cache_reader import (  # noqa: E402
    DEFAULT_CACHE_PATH,
    get_html,
    load_cache_subset,
)
from dex_etymology import (  # noqa: E402
    USER_KEEP,
    categorize_etymology,
    get_verb_etymology,
    is_accepted,
)

# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

CSV_PATH: Final[Path] = (
    _PROJECT_ROOT / "data" / "romanian_lexicon_with_freq.csv"
)
CACHE_PATH: Final[Path] = DEFAULT_CACHE_PATH

TOP_N_FREQUENCY: Final[int] = 1_000

DORSAL_STEMS: Final[frozenset[str]] = frozenset({"c", "g"})
CORONAL_STEMS: Final[frozenset[str]] = frozenset({"t", "d", "s"})
TARGET_STEMS: Final[frozenset[str]] = DORSAL_STEMS | CORONAL_STEMS
PLURAL_FRONT_OPPS: Final[frozenset[str]] = frozenset({"i", "e"})
ACCEPTED_VERBALIZERS: Final[frozenset[str]] = frozenset({"-i", "-a"})


# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------


class StemClass(StrEnum):
    """Articulatory class of the noun stem-final consonant."""

    DORSAL = "dorsal"
    CORONAL = "coronal"


class PluralAlt(StrEnum):
    """Whether the noun's plural form palatalizes the stem-final."""

    ALTERNATES = "alt"
    NON_ALTERNATING = "nonalt"


class Verbalizer(StrEnum):
    """Verbalizer allomorph the derived verb selects."""

    FRONT_I = "i"
    BACK_A = "a"


@dataclass(frozen=True, slots=True)
class NounRecord:
    """A single noun row that passes the contingency-table row filters.

    Frozen + slots: this object flows between layers and shouldn't
    mutate. ``derived_verbs`` is a tuple so the dataclass is hashable
    if the caller needs it.
    """

    lemma: str
    derived_verbs: tuple[str, ...]
    stem_final: str
    opportunity: str
    deriv_suffix: str  # "-i" or "-a"
    mutation: bool
    frequency: float

    @property
    def stem_class(self) -> StemClass:
        return (
            StemClass.DORSAL
            if self.stem_final in DORSAL_STEMS
            else StemClass.CORONAL
        )

    @property
    def plural_alt(self) -> PluralAlt:
        return (
            PluralAlt.ALTERNATES
            if self.mutation
            else PluralAlt.NON_ALTERNATING
        )

    @property
    def verbalizer(self) -> Verbalizer:
        return (
            Verbalizer.FRONT_I
            if self.deriv_suffix == "-i"
            else Verbalizer.BACK_A
        )


@dataclass(frozen=True, slots=True)
class PanelStats:
    """Fisher 2x2 result for one stem-class panel."""

    alt_i: int
    alt_a: int
    nonalt_i: int
    nonalt_a: int

    @property
    def n(self) -> int:
        return self.alt_i + self.alt_a + self.nonalt_i + self.nonalt_a

    @property
    def is_complete(self) -> bool:
        """All four marginals non-zero (required by Fisher's exact)."""
        return (
            self.alt_i + self.alt_a > 0
            and self.nonalt_i + self.nonalt_a > 0
            and self.alt_i + self.nonalt_i > 0
            and self.alt_a + self.nonalt_a > 0
        )

    def fisher(self) -> tuple[float, float] | None:
        if not self.is_complete:
            return None
        odds, p = fisher_exact(
            [[self.alt_i, self.alt_a], [self.nonalt_i, self.nonalt_a]]
        )
        return float(odds), float(p)

    def format_line(self, name: str) -> str:
        head = (
            f"  {name:8s}:  alt-pl  {self.alt_i:3d}  -i / "
            f"{self.alt_a:3d}  -a    "
            f"nonalt-pl  {self.nonalt_i:3d}  -i / "
            f"{self.nonalt_a:3d}  -a    n = {self.n:4d}"
        )
        fisher = self.fisher()
        if fisher is None:
            return head
        odds, p = fisher
        return f"{head}    OR = {odds:.2f}    p = {p:.3f}"


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def _load_records() -> tuple[list[NounRecord], dict[str, float]]:
    """Read the lexicon CSV and return:
    - the NounRecord list passing the row filters
    - a {lemma → frequency} dict for ALL nouns (needed for top-N ranking
      independently of whether the noun is in the contingency table)

    Filters applied here (the trust boundary):
      pos == "N"
      stem_final in TARGET_STEMS
      mutation field is non-empty
      opportunity in {i, e, uri}
      deriv_suffix in {-i, -a}
      derived_verbs non-empty
      exception_reason not starting with "nde:"
    """
    records: list[NounRecord] = []
    seen: set[str] = set()
    frequency_by_lemma: dict[str, float] = {}

    with CSV_PATH.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["pos"] != "N":
                continue
            lemma = row["lemma"]
            try:
                freq = float(row.get("freq_ron_news_2024_1M") or 0)
            except ValueError:
                freq = 0.0
            frequency_by_lemma[lemma] = freq

            stem_final = row["stem_final"]
            if stem_final not in TARGET_STEMS:
                continue
            if row["mutation"].strip() == "":
                continue
            opp = row["opportunity"].strip()
            if opp not in PLURAL_FRONT_OPPS and opp != "uri":
                continue
            deriv_suffix = row["deriv_suffixes"].strip()
            if deriv_suffix not in ACCEPTED_VERBALIZERS:
                continue
            derived_verbs_field = (row.get("derived_verbs") or "").strip()
            if not derived_verbs_field:
                continue
            if (row.get("exception_reason") or "").startswith("nde:"):
                continue
            if lemma in seen:
                continue
            seen.add(lemma)
            verbs = tuple(
                v.strip() for v in derived_verbs_field.split("|") if v.strip()
            )
            records.append(
                NounRecord(
                    lemma=lemma,
                    derived_verbs=verbs,
                    stem_final=stem_final,
                    opportunity=opp,
                    deriv_suffix=deriv_suffix,
                    mutation=row["mutation"].strip().upper() == "TRUE",
                    frequency=freq,
                )
            )
    return records, frequency_by_lemma


# ---------------------------------------------------------------------------
# Etymology audit
# ---------------------------------------------------------------------------


def _passes_etymology_audit(record: NounRecord, cache: dict[str, str]) -> bool:
    """At least one derived verb must classify as A, B, or D, OR the
    lemma must appear in USER_KEEP."""
    if record.lemma in USER_KEEP:
        return True
    for verb in record.derived_verbs:
        html = get_html(cache, verb)
        if html is None:
            continue
        etymology = get_verb_etymology(html, verb)
        if is_accepted(categorize_etymology(record.lemma, etymology)):
            return True
    return False


# ---------------------------------------------------------------------------
# Table building
# ---------------------------------------------------------------------------


_ScopeFilter = Callable[[NounRecord], bool]


def _build_panels(
    records: list[NounRecord],
    cache: dict[str, str],
    scope: _ScopeFilter,
) -> tuple[dict[StemClass | str, PanelStats], int]:
    """Build dorsal, coronal, and pooled panels for the records that
    pass ``scope`` and the etymology audit.

    Returns the three panels keyed by ``StemClass`` plus a sentinel
    "pooled" string, and the count of kept records.
    """
    cells: dict[StemClass | str, Counter[tuple[PluralAlt, Verbalizer]]] = {
        StemClass.DORSAL: Counter(),
        StemClass.CORONAL: Counter(),
        "pooled": Counter(),
    }
    n_kept = 0
    for record in records:
        if not scope(record):
            continue
        if not _passes_etymology_audit(record, cache):
            continue
        n_kept += 1
        key = (record.plural_alt, record.verbalizer)
        cells[record.stem_class][key] += 1
        cells["pooled"][key] += 1

    def pack(c: Counter[tuple[PluralAlt, Verbalizer]]) -> PanelStats:
        return PanelStats(
            alt_i=c[(PluralAlt.ALTERNATES, Verbalizer.FRONT_I)],
            alt_a=c[(PluralAlt.ALTERNATES, Verbalizer.BACK_A)],
            nonalt_i=c[(PluralAlt.NON_ALTERNATING, Verbalizer.FRONT_I)],
            nonalt_a=c[(PluralAlt.NON_ALTERNATING, Verbalizer.BACK_A)],
        )

    return (
        {k: pack(v) for k, v in cells.items()},
        n_kept,
    )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    records, frequency_by_lemma = _load_records()
    print(
        f"Loaded {len(records)} (lemma, dvs, ...) tuples "
        f"that pass row filters."
    )

    unique_verbs = {v for r in records for v in r.derived_verbs}
    cache = load_cache_subset(unique_verbs, cache_path=CACHE_PATH)
    print(
        f"Cached DEX pages: {len(cache)} "
        f"(for {len(unique_verbs)} unique verbs)"
    )

    ranked = sorted(frequency_by_lemma.items(), key=lambda kv: -kv[1])
    top_lemmas: frozenset[str] = frozenset(
        lemma for lemma, _ in ranked[:TOP_N_FREQUENCY]
    )

    scopes: list[tuple[str, _ScopeFilter]] = [
        (f"TOP-{TOP_N_FREQUENCY}", lambda r: r.lemma in top_lemmas),
        ("FULL DEX", lambda r: True),
    ]

    print("\n" + "=" * 78)
    print(
        "Contingency table: row = plural alternates (mutation),"
        " column = verbalizer allomorph"
    )
    print("=" * 78)
    for label, scope in scopes:
        panels, n_kept = _build_panels(records, cache, scope)
        print(f"\n{label}  (n_kept = {n_kept})")
        for name in (StemClass.DORSAL, StemClass.CORONAL, "pooled"):
            stats = panels[name]
            display_name = name.value if isinstance(name, StemClass) else name
            print(stats.format_line(display_name))


if __name__ == "__main__":
    main()
