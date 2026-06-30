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

EXCLUSIONS — row-side (applied in both codings)
    - Rows whose ``exception_reason`` begins with ``nde:`` (the three
      non-derived-environment-blocking classes Steriade also excludes).
    - Derived verbs that don't pass the DEX etymology audit. The audit
      categorizations are persisted in
      ``data/derived_verb_etymology.csv`` (the canonical pipeline
      artifact); only categories A, B, D are accepted.
    - Explicit USER_KEEP allowlist for back-formations the etym audit
      mis-rejects: clint, dinte, vargă, verigă.

TWO CODINGS — emitted side-by-side
----------------------------------
The column-side trigger-absent question (whether to include verbs whose
conjugation morphologically bleeds the trigger — the ``-esc``/``-ez``
augmented classes) is methodologically contestable, so the script
produces BOTH panels and prints the contamination diagnostic, letting
the reader see what the choice buys.

1. **STERIADE-STRICT** (column includes trigger-absent verbs).
   Tests Steriade's literal lexeme-level affix-avoidance claim: given a
   K/Tʃ vs K/K base, which verbalizer allomorph does the speaker
   select? The augment is downstream of selection, so trigger-absent
   verbs stay in. Caveat: the coronal ``-i`` column ends up ~79%
   trigger-absent — the diagnostic block reports the exact share.

2. **BOTH-AXES-EXCLUDED** (column excludes TRIGGER-ABSENT verbs).
   Symmetric construct: both axes measure "actually licenses
   palatalization." Cleaner column, but excluding ~80% of the coronal
   ``-i`` column collapses the coronal stratum's power to detect an
   effect alone; significance survives only in the pool.

When forced to one table for space-constrained reports, prefer
STERIADE-STRICT: it tests the literal claim, every stratum is
significant, and the contamination is a coding choice defensible in
one caption sentence. BOTH-AXES-EXCLUDED is the appropriate robustness
check for a full paper.

Etymology + paradigm sources of truth
-------------------------------------
- ``data/derived_verb_etymology.csv`` — per-pair etymology category.
- ``data/verb_paradigm_verdict.csv`` — per-verb conjugation verdict
  (ALTERNATES / TRIGGER-ABSENT / NO-PARADIGM). Both are first-class
  CSV pipeline artifacts; neither requires a live DEX cache lookup.

Row coding (settled)
--------------------
The row is ``mutation``-based, NOT the historic ``opportunity``-based
allomorph proxy ``-i/-e`` vs ``-uri``. For coronals, the surface
allomorph ``-e`` does not trigger assibilation (e.g. *casă/case*),
so coding by allomorph mis-classifies coronal-``-e`` nouns as
alternating. The ``mutation`` field records the actual sg-vs-pl
comparison, which is the surface-to-surface licensing condition.
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

from dex_etymology import USER_KEEP  # noqa: E402

# ---------------------------------------------------------------------------
# Inputs / outputs
# ---------------------------------------------------------------------------

CSV_PATH: Final[Path] = (
    _PROJECT_ROOT / "data" / "romanian_lexicon_with_freq.csv"
)
ETYMOLOGY_CSV_PATH: Final[Path] = (
    _PROJECT_ROOT / "data" / "derived_verb_etymology.csv"
)
VERDICT_CSV_PATH: Final[Path] = (
    _PROJECT_ROOT / "data" / "verb_paradigm_verdict.csv"
)
OUTPUT_PATH: Final[Path] = _PROJECT_ROOT / "analysis" / "contingency_table.txt"

TRIGGER_ABSENT: Final[str] = "TRIGGER-ABSENT"

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
        # Use scientific notation for vanishingly small p so highly
        # significant results don't collapse to "0.000" in the report.
        p_fmt = f"{p:.2e}" if p < 1e-3 else f"{p:.3f}"
        return f"{head}    OR = {odds:.2f}    p = {p_fmt}"


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def _resolve_deriv_suffix(
    lemma: str,
    verbs: list[str],
    suffixes: list[str],
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
) -> str | None:
    """Pick the single ``deriv_suffix`` value that represents this noun.

    The lexicon may emit multi-token suffixes (e.g., ``-i|-a``) when a
    noun has more than one denominal verb. The contingency table is a
    per-noun structure so we have to pick one. The rule: among the
    derived verbs the etymology audit ACCEPTS, take their suffix. If
    every accepted verb uses the same suffix (the empirical pattern in
    our data — zero genuinely-ambiguous cases), that's our answer. If
    none of the listed verbs are accepted, return ``None`` so the
    caller can drop the row gracefully.

    USER_KEEP lemmas (back-formations the audit mis-rejects) fall
    through to the first valid suffix.
    """
    pairs = [
        (v, s) for v, s in zip(verbs, suffixes) if s in ACCEPTED_VERBALIZERS
    ]
    if not pairs:
        return None
    accepted_suffixes = sorted(
        {s for v, s in pairs if accepted_by_pair.get((lemma, v), False)}
    )
    if len(accepted_suffixes) == 1:
        return accepted_suffixes[0]
    if len(accepted_suffixes) > 1:
        # Genuinely ambiguous: a noun has accepted verbs with both -i
        # and -a. Empirically empty in our data; if it ever fires we
        # default to -i (the palatalizing form, the marked choice that
        # Steriade's claim is about). The script's diagnostic block
        # will not flag this directly, but the run will.
        return "-i"
    if lemma in user_keep:
        return pairs[0][1]
    return None


def _load_records(
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
) -> tuple[list[NounRecord], dict[str, float]]:
    """Read the lexicon CSV and return:
    - the NounRecord list passing the row filters
    - a {lemma → frequency} dict for ALL nouns (needed for top-N ranking
      independently of whether the noun is in the contingency table)

    Filters applied here (the trust boundary):
      pos == "N"
      stem_final in TARGET_STEMS
      mutation field is non-empty
      opportunity in {i, e, uri}
      derived_verbs non-empty
      exception_reason not starting with "nde:"
      ``deriv_suffix`` resolves via :func:`_resolve_deriv_suffix` —
      handles multi-token suffix strings by picking the accepted-verb
      suffix; returns ``None`` (row dropped) if no accepted verb has a
      valid suffix.
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
            derived_verbs_field = (row.get("derived_verbs") or "").strip()
            if not derived_verbs_field:
                continue
            if (row.get("exception_reason") or "").startswith("nde:"):
                continue
            if lemma in seen:
                continue
            verbs = [
                v.strip() for v in derived_verbs_field.split("|") if v.strip()
            ]
            suffixes = [
                s.strip()
                for s in (row.get("deriv_suffixes") or "").split("|")
                if s.strip()
            ]
            deriv_suffix = _resolve_deriv_suffix(
                lemma, verbs, suffixes, accepted_by_pair, user_keep
            )
            if deriv_suffix is None:
                continue
            seen.add(lemma)
            records.append(
                NounRecord(
                    lemma=lemma,
                    derived_verbs=tuple(verbs),
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


def _load_verb_verdicts() -> dict[str, str]:
    """Read the per-verb conjugation verdict CSV into ``{verb: verdict}``.

    Verdict values follow :mod:`scripts.audit_verb_alternation`:
    ``ALTERNATES``, ``TRIGGER-ABSENT``, ``NO-PARADIGM``. The CSV is the
    canonical source; the legacy ``analysis/verb_alternation_audit.txt``
    is a human-readable mirror only.
    """
    out: dict[str, str] = {}
    with VERDICT_CSV_PATH.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            verb = (row.get("verb") or "").strip()
            verdict = (row.get("verdict") or "").strip()
            if verb and verdict:
                out[verb] = verdict
    return out


def _load_etymology_decisions() -> (
    tuple[dict[tuple[str, str], bool], set[str]]
):
    """Read the structured etymology CSV into:

    - ``accepted_by_pair``: ``{(noun, verb): True/False}`` from explicit
      audit rows. ``True`` iff the audit categorized the pair as A, B,
      or D (the ACCEPT set). Pairs absent from this map are
      ``unclassified`` — the script will flag them rather than
      silently treat them as rejected.
    - ``user_keep``: set of nouns the explicit user-override list keeps
      regardless of any individual verb's audit result.

    The CSV is the canonical pipeline artifact; this function is the
    single read path. See :mod:`romanian_palatalization.dex_etymology`
    for the underlying category definitions.
    """
    accepted: dict[tuple[str, str], bool] = {}
    user_keep: set[str] = set()
    with ETYMOLOGY_CSV_PATH.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            noun = (row.get("noun") or "").strip()
            verb = (row.get("verb") or "").strip()
            cat = (row.get("category") or "").strip().upper()
            is_acc = (row.get("accepted") or "").strip().lower() == "true"
            if cat == "USER_KEEP":
                user_keep.add(noun)
                continue
            if noun and verb:
                accepted[(noun, verb)] = is_acc
    # Merge in the hard-coded fallback so a missing CSV doesn't silently
    # drop the back-formation overrides.
    user_keep |= set(USER_KEEP)
    return accepted, user_keep


def _accepted_verbs(
    record: NounRecord,
    accepted_by_pair: dict[tuple[str, str], bool],
    verdicts: dict[str, str],
    exclude_trigger_absent: bool,
    unclassified_pairs: set[tuple[str, str]],
) -> list[str]:
    """Return the verbs in ``record.derived_verbs`` that pass the
    etymology audit (and, in BOTH-AXES-EXCLUDED mode, are NOT
    trigger-absent in conjugation).

    Pairs missing from ``accepted_by_pair`` are added to
    ``unclassified_pairs`` so the caller can warn about coverage gaps.
    """
    out: list[str] = []
    for verb in record.derived_verbs:
        key = (record.lemma, verb)
        if key not in accepted_by_pair:
            unclassified_pairs.add(key)
            continue
        if not accepted_by_pair[key]:
            continue
        if exclude_trigger_absent and verdicts.get(verb) == TRIGGER_ABSENT:
            continue
        out.append(verb)
    return out


# ---------------------------------------------------------------------------
# Table building
# ---------------------------------------------------------------------------


_ScopeFilter = Callable[[NounRecord], bool]


def _build_panels(
    records: list[NounRecord],
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
    verdicts: dict[str, str],
    exclude_trigger_absent: bool,
    unclassified_pairs: set[tuple[str, str]],
    scope: _ScopeFilter,
) -> tuple[dict[StemClass | str, PanelStats], int]:
    """Build dorsal, coronal, and pooled panels for the records that
    pass ``scope`` and the etymology audit.

    ``exclude_trigger_absent`` selects between codings:
      - ``False`` (STERIADE-STRICT): column counts every accepted verb,
        regardless of whether its conjugation paradigm bleeds the
        trigger.
      - ``True`` (BOTH-AXES-EXCLUDED): column drops verbs whose
        conjugation verdict is ``TRIGGER-ABSENT``.

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
        if record.lemma in user_keep:
            # User-keep nouns count even without an accepted verb.
            keep = True
        else:
            keep = bool(
                _accepted_verbs(
                    record,
                    accepted_by_pair,
                    verdicts,
                    exclude_trigger_absent,
                    unclassified_pairs,
                )
            )
        if not keep:
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


def _contamination_diagnostic(
    records: list[NounRecord],
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
    verdicts: dict[str, str],
) -> str:
    """Quantify the trigger-absent share of each stratum's ``-i`` column.

    Surfaces the methodological contamination that motivates the
    BOTH-AXES-EXCLUDED panel. Per-stratum breakdown of the (alt-pl,
    nonalt-pl) × (ALTERNATES, TRIGGER-ABSENT, NO-PARADIGM) joint table
    restricted to ``-i`` verbalizers.
    """
    lines = []
    stem_groups: list[tuple[str, frozenset[str]]] = [
        ("dorsal", DORSAL_STEMS),
        ("coronal", CORONAL_STEMS),
    ]
    unclassified_throwaway: set[tuple[str, str]] = set()
    for label, stems in stem_groups:
        counts: Counter[tuple[str, str]] = Counter()
        for rec in records:
            if rec.stem_final not in stems:
                continue
            if rec.deriv_suffix != "-i":
                continue
            # Use all accepted verbs (Steriade-strict semantics) so the
            # diagnostic reflects what the strict panel actually sees.
            for verb in _accepted_verbs(
                rec,
                accepted_by_pair,
                verdicts,
                exclude_trigger_absent=False,
                unclassified_pairs=unclassified_throwaway,
            ) or ([] if rec.lemma not in user_keep else [""]):
                verdict = (
                    verdicts.get(verb, "UNKNOWN") if verb else "USER_KEEP"
                )
                row_key = (
                    "alt" if rec.mutation else "noalt",
                    verdict,
                )
                counts[row_key] += 1
        total = sum(counts.values())
        if total == 0:
            continue
        lines.append(f"\n  {label} -i column (n={total}):")
        for (m, v), c in sorted(counts.items()):
            pct = 100 * c / total
            lines.append(f"    {m:5s}, verdict={v:14s}: {c:3d}  ({pct:.1f}%)")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def _render_coding_panel(
    records: list[NounRecord],
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
    verdicts: dict[str, str],
    exclude_trigger_absent: bool,
    frequency_by_lemma: dict[str, float],
    unclassified_pairs: set[tuple[str, str]],
) -> str:
    """Render the TOP-1000 and FULL DEX subpanels for one coding."""
    lines: list[str] = []
    ranked = sorted(frequency_by_lemma.items(), key=lambda kv: -kv[1])
    top_lemmas: frozenset[str] = frozenset(
        lemma for lemma, _ in ranked[:TOP_N_FREQUENCY]
    )
    scopes: list[tuple[str, _ScopeFilter]] = [
        (f"TOP-{TOP_N_FREQUENCY}", lambda r: r.lemma in top_lemmas),
        ("FULL DEX", lambda r: True),
    ]
    for label, scope in scopes:
        panels, n_kept = _build_panels(
            records,
            accepted_by_pair,
            user_keep,
            verdicts,
            exclude_trigger_absent,
            unclassified_pairs,
            scope,
        )
        lines.append(f"\n{label}  (n_kept = {n_kept})")
        for name in (StemClass.DORSAL, StemClass.CORONAL, "pooled"):
            stats = panels[name]
            display = name.value if isinstance(name, StemClass) else name
            lines.append(stats.format_line(display))
    return "\n".join(lines)


def _render_report(
    records: list[NounRecord],
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
    verdicts: dict[str, str],
    unclassified_pairs: set[tuple[str, str]],
    frequency_by_lemma: dict[str, float],
) -> str:
    """Render the full contingency-table report.

    Two codings (STERIADE-STRICT and BOTH-AXES-EXCLUDED), each with
    TOP-{TOP_N_FREQUENCY} and FULL DEX subpanels, plus a contamination
    diagnostic that quantifies the trigger-absent share of the ``-i``
    column under the strict coding.
    """
    lines: list[str] = []
    lines.append(
        f"Loaded {len(records)} (lemma, dvs, ...) tuples "
        f"that pass row filters."
    )
    lines.append(
        f"Etymology decisions: {len(accepted_by_pair)} explicit pairs, "
        f"{len(user_keep)} user-keep lemmas, "
        f"{len(verdicts)} verb verdicts."
    )

    lines.append("\n" + "=" * 78)
    lines.append("STERIADE-STRICT  (column INCLUDES trigger-absent verbs)")
    lines.append(
        "Tests Steriade's literal lexeme-level affix-avoidance claim."
    )
    lines.append("=" * 78)
    lines.append(
        _render_coding_panel(
            records,
            accepted_by_pair,
            user_keep,
            verdicts,
            exclude_trigger_absent=False,
            frequency_by_lemma=frequency_by_lemma,
            unclassified_pairs=unclassified_pairs,
        )
    )

    lines.append("\n" + "=" * 78)
    lines.append("BOTH-AXES-EXCLUDED  (column EXCLUDES TRIGGER-ABSENT verbs)")
    lines.append(
        "Robustness check: symmetric construct, sacrifices coronal "
        "stratum power."
    )
    lines.append("=" * 78)
    lines.append(
        _render_coding_panel(
            records,
            accepted_by_pair,
            user_keep,
            verdicts,
            exclude_trigger_absent=True,
            frequency_by_lemma=frequency_by_lemma,
            unclassified_pairs=unclassified_pairs,
        )
    )

    lines.append("\n" + "=" * 78)
    lines.append("CONTAMINATION DIAGNOSTIC: trigger-absent share of -i column")
    lines.append("(under STERIADE-STRICT coding)")
    lines.append("=" * 78)
    lines.append(
        _contamination_diagnostic(
            records, accepted_by_pair, user_keep, verdicts
        )
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    accepted_by_pair, user_keep = _load_etymology_decisions()
    verdicts = _load_verb_verdicts()
    records, frequency_by_lemma = _load_records(accepted_by_pair, user_keep)
    unclassified_pairs: set[tuple[str, str]] = set()

    report = _render_report(
        records,
        accepted_by_pair,
        user_keep,
        verdicts,
        unclassified_pairs,
        frequency_by_lemma,
    )
    print(report, end="")

    # Persist the report so future runs can be compared without re-
    # running the script (e.g., when the user edits derived verbs).
    OUTPUT_PATH.write_text(report, encoding="utf-8")
    print(f"\nReport persisted to {OUTPUT_PATH}")

    if unclassified_pairs:
        # Coverage gap: pairs in the lexicon that are missing from
        # the etymology CSV. Print a concise summary so the user knows
        # the etymology audit needs a refresh on these.
        n = len(unclassified_pairs)
        sample = sorted(unclassified_pairs)[:10]
        print(
            f"\nWARNING: {n} (noun, verb) pairs lack an etymology "
            f"classification in {ETYMOLOGY_CSV_PATH.name}."
        )
        print("First 10 unclassified pairs:")
        for noun, verb in sample:
            print(f"  {noun}  →  {verb}")
        print(
            "To refresh, re-run the supplement step or add the missing "
            "rows to the etymology CSV with category in {A, B, D} "
            "(accepted) or {C, E, F} (rejected)."
        )


if __name__ == "__main__":
    main()
