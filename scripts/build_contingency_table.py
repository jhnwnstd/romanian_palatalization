#!/usr/bin/env python3
"""Build the inflection-dependence contingency table.

Canonical script for the table Steriade 2008 §10.2.3.5 offers as
support for her paradigm-consultation account. Replicates her lexical
counts on the modern lexicon, and reports the same test under three
additional codings so the reader can see how sensitive the effect is
to methodological choices.

The claim under test
--------------------
For a noun with a target obstruent stem-final, the plural allomorph
and the verbalizer allomorph co-vary: nouns whose plural palatalizes
the stem-final tend to select the palatalizing verbalizer ``-i``;
nouns whose plural doesn't palatalize tend to select ``-a``.

Axes and coding
---------------
Two axes, two coding choices each → four panels.

ROW = "does the plural palatalize?"
  (a) MUTATION-FLAG — the lexicon's ``mutation`` boolean, a direct
      surface-orth comparison of singular vs plural. Cheap and
      Steriade-original.
  (b) IPA-SG-vs-PL — compare the last consonantal segment of the
      singular's normalised IPA to the same position in the plural's,
      after stripping known plural-suffix material. Three empirical
      states:
        BASE→BASE   → K/K   (non-alternating)
        BASE→PAL    → K/T͡ʃ ALTERNATING (mutation)
        PAL→PAL     → K/T͡ʃ LEXICAL (base already palatal — arici,
                       pace, elogiu; behaviourally K/T͡ʃ under both
                       theories, so pooled with ALTERNATING; the split
                       is reported for transparency)
        PAL→BASE    → inconsistent, dropped
      Catches (i) cluster cases where the orthographic comparison
      conflates the sc→ʃt cascade into a single ``True`` flag; (ii)
      lexical-palatal bases that the orthographic comparison
      correctly marks non-mutating but that behave as K/T͡ʃ.

COL = "does the derivation palatalize?"
  (a) SUFFIX-ONLY — ``-i`` vs ``-a``. Steriade's literal claim is
      about affix avoidance, which is what this measures.
  (b) ALTERNATES-ONLY — a ``-i`` verb only counts as PAL evidence
      if its paradigm verdict is ``ALTERNATES`` (verb actually surfaces
      the palatalization in a diagnostic cell). ``TRIGGER-ABSENT``
      verbs (whose ``-esc``/``-ez`` augment bleeds the phonological
      trigger) test morphology, not phonology, and are dropped —
      they can't distinguish "the base can palatalize" from "the base
      can't". ``NO-PARADIGM`` verbs are also dropped (no evidence).

The four panels
---------------
    row = MUTATION,  col = SUFFIX-ONLY       → LITERAL-STERIADE
    row = MUTATION,  col = ALTERNATES-ONLY   → MUTATION-ALT
    row = IPA,       col = SUFFIX-ONLY       → IPA-SUFFIX
    row = IPA,       col = ALTERNATES-ONLY   → PRIMARY  (recommended)

The PRIMARY panel is the honest test: both axes measure surface
palatalization, neither is contaminated by allomorph choice masquerading
as phonology or by augment-bled morphology masquerading as inertness.
LITERAL-STERIADE reproduces the paper's own coding exactly; the other
two are robustness cells of the 2×2 factorial. If all four agree in
direction and significance, the effect is armoured against coding.

Exclusions applied everywhere
-----------------------------
- ``pos == 'N'``; ``stem_final in {c, g, t, d, s, z}``
- ``exception_reason`` doesn't start with ``nde:`` (Steriade's NDEB
  excludes)
- ``derived_verbs`` field non-empty
- Nouns whose etymology-accepted verbs use *both* ``-i`` and ``-a``
  emit two records (one per suffix), each carrying only the verbs of
  that suffix. This keeps the ambiguous noun's evidence honest instead
  of collapsing it to whichever suffix happens to be preferred.
- Column-side: etymology audit accepts the (noun, verb) pair, or the
  noun is in the ``USER_KEEP`` back-formation allowlist
- Column-side: verbalizer restricted to ``{-i, -a}``. ``-ui`` is a
  third suffix Steriade treats separately; reported in a small sidebar,
  not in the main table
- Row-side: rows whose IPA lemma or plural is missing/foreign-only get
  a UNCLEAR classification and are dropped from the IPA-row panels
  (but retained in the MUTATION-row panels, which don't need IPA)

Two frequency scopes
--------------------
- TOP-1000: the top-1000 nouns by frequency, matching the size of
  Steriade's own hand-collected sample.
- FULL DEX: the modern-lexicon replication.

Statistics
----------
Fisher's exact 2×2 per stem-class panel (dorsal / coronal / pooled).
Fisher over chi-square because per-stem-class cells can be small.
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
from scipy.stats.contingency import odds_ratio

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PROJECT_ROOT / "src"))

from dex_etymology import USER_KEEP  # noqa: E402

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
ALTERNATES: Final[str] = "ALTERNATES"

TOP_N_FREQUENCY: Final[int] = 1_000

DORSAL_STEMS: Final[frozenset[str]] = frozenset({"c", "g"})
CORONAL_STEMS: Final[frozenset[str]] = frozenset({"t", "d", "s", "z"})
TARGET_STEMS: Final[frozenset[str]] = DORSAL_STEMS | CORONAL_STEMS
PLURAL_FRONT_OPPS: Final[frozenset[str]] = frozenset({"i", "e"})
ACCEPTED_VERBALIZERS: Final[frozenset[str]] = frozenset({"-i", "-a"})


# ---------------------------------------------------------------------------
# IPA classifier for the row axis
# ---------------------------------------------------------------------------
#
# The row axis under IPA coding is a three-state comparison of the
# last stem-final consonant in the singular's normalised IPA to the
# same position in the plural's. Both sides use the shared tokeniser
# below so tie-barred affricates count as one segment.

BASE_IPA: Final[dict[str, str]] = {
    "c": "k",
    "g": "ɡ",
    "t": "t",
    "d": "d",
    "s": "s",
    "z": "z",
}
PAL_IPA: Final[dict[str, str]] = {
    "c": "t͡ʃ",
    "g": "d͡ʒ",
    "t": "t͡s",
    "d": "z",
    "s": "ʃ",
    "z": "ʒ",
}

_PLURAL_SUFFIXES: Final[tuple[str, ...]] = (
    "uri",
    "urj",
    "urʲ",
    "ele",
    "ile",
    "j",
    "ʲ",
    "i",
    "e",
)
# Singular strippable material — only vowels that participate in
# inflectional class (-ă/-ə for feminine, -e for -e class). Other
# vowels (o, ɨ) at word end are treated as stem-internal.
_SING_SUFFIXES: Final[tuple[str, ...]] = ("ə", "e")

_VOWELS: Final[set[str]] = set("aeiouəɨɛɔ")
_GLIDES: Final[set[str]] = set("jw")
_MARKS: Final[set[str]] = set("ˈˌ.'` ")


def _tokens(ipa: str) -> list[str]:
    """Split IPA into segment atoms; tie-barred affricates (t͡ʃ etc.)
    are one atom of three code points."""
    out: list[str] = []
    i = 0
    while i < len(ipa):
        ch = ipa[i]
        if ch in _MARKS:
            i += 1
            continue
        if i + 2 < len(ipa) and ipa[i + 1] == "͡":
            out.append(ipa[i : i + 3])
            i += 3
        else:
            out.append(ch)
            i += 1
    return out


def _first_variant(field_val: str) -> str:
    """Return the first non-empty, non-foreign-language variant of a
    pipe-separated IPA field. ``fr:absorbant``-style rows are dropped."""
    for v in field_val.split(" | "):
        v = v.strip()
        if not v:
            continue
        if len(v) >= 3 and v[2] == ":" and v[:2].isalpha() and v[:2].islower():
            continue
        return v
    return ""


def _strip_suffix(ipa: str, suffixes: tuple[str, ...]) -> str:
    """Strip the longest matching suffix once."""
    for s in sorted(suffixes, key=len, reverse=True):
        if ipa.endswith(s):
            return ipa[: -len(s)]
    return ipa


def _last_consonantal(ipa: str) -> str | None:
    """Walk reverse-tokens, skipping trailing vowels/glides/ʲ,
    return the first consonantal segment we find."""
    for tok in reversed(_tokens(ipa)):
        head = tok[0]
        if head in _VOWELS or head in _GLIDES or tok == "ʲ":
            continue
        return tok
    return None


def _form_class(ipa: str, sf: str, suffixes: tuple[str, ...]) -> str:
    """Classify one form ('pal' / 'base' / 'other') for stem-final sf.

    Cluster carve-out: for s-stems and c-stems, if /ʃ/ appears within
    the last 4 chars of the stripped stem, the form counts as
    palatalised (the sc→ʃt cascade puts ʃ before a coronal stop, not
    at the very end).
    """
    ipa = _first_variant(ipa)
    if not ipa:
        return "other"
    stem = _strip_suffix(ipa, suffixes)
    if not stem:
        return "other"
    last = _last_consonantal(stem)
    if last is None:
        return "other"
    if last == PAL_IPA[sf]:
        return "pal"
    if last == BASE_IPA[sf]:
        return "base"
    if sf in ("s", "c") and "ʃ" in stem[-4:]:
        return "pal"
    return "other"


class IPARowClass(StrEnum):
    """Three-state outcome of the sg-vs-pl IPA comparison."""

    NONALT = "nonalt"  # sg base, pl base
    ALT = "alt"  # sg base, pl palatal (a real mutation)
    LEX = "lex"  # sg palatal, pl palatal (already palatal)
    UNCLEAR = "unclear"  # missing IPA or PAL→BASE inconsistency


def _classify_ipa_row(sg_ipa: str, pl_ipa: str, sf: str) -> IPARowClass:
    """Three-state classification of the plural mutation via sg-vs-pl."""
    sg = _form_class(sg_ipa, sf, _SING_SUFFIXES)
    pl = _form_class(pl_ipa, sf, _PLURAL_SUFFIXES)
    if sg == "other" or pl == "other":
        return IPARowClass.UNCLEAR
    if sg == "base" and pl == "pal":
        return IPARowClass.ALT
    if sg == "base" and pl == "base":
        return IPARowClass.NONALT
    if sg == "pal" and pl == "pal":
        return IPARowClass.LEX
    return IPARowClass.UNCLEAR  # pal→base is inconsistent


# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------


class StemClass(StrEnum):
    DORSAL = "dorsal"
    CORONAL = "coronal"


class RowCoding(StrEnum):
    """How the row axis (plural class) is coded."""

    MUTATION = "mutation"  # orthographic ``mutation`` field
    IPA = "ipa"  # sg-vs-pl IPA comparison


class ColCoding(StrEnum):
    """How the column axis (verbalizer class) is coded."""

    SUFFIX_ONLY = "suffix"  # -i vs -a, no paradigm filter
    ALTERNATES_ONLY = "alternates"  # -i verbs must be ALTERNATES verdict


@dataclass(frozen=True, slots=True)
class NounRecord:
    """A single noun row that passes the row-side filters."""

    lemma: str
    derived_verbs: tuple[str, ...]
    stem_final: str
    opportunity: str
    deriv_suffix: str
    mutation: bool
    frequency: float
    ipa_row_class: IPARowClass

    @property
    def stem_class(self) -> StemClass:
        return (
            StemClass.DORSAL
            if self.stem_final in DORSAL_STEMS
            else StemClass.CORONAL
        )

    def plural_alt(self, coding: RowCoding) -> bool | None:
        """Return True (alt) / False (nonalt) / None (drop) under the
        given row coding."""
        if coding is RowCoding.MUTATION:
            return self.mutation
        if self.ipa_row_class is IPARowClass.ALT:
            return True
        if self.ipa_row_class is IPARowClass.LEX:
            return True  # pooled with ALT — both are K/T͡ʃ class
        if self.ipa_row_class is IPARowClass.NONALT:
            return False
        return None  # UNCLEAR


@dataclass(frozen=True, slots=True)
class PanelStats:
    """Fisher 2×2 result for one panel."""

    alt_i: int
    alt_a: int
    nonalt_i: int
    nonalt_a: int

    @property
    def n(self) -> int:
        return self.alt_i + self.alt_a + self.nonalt_i + self.nonalt_a

    @property
    def is_complete(self) -> bool:
        return (
            self.alt_i + self.alt_a > 0
            and self.nonalt_i + self.nonalt_a > 0
            and self.alt_i + self.nonalt_i > 0
            and self.alt_a + self.nonalt_a > 0
        )

    def fisher(self) -> tuple[float, float, float, float] | None:
        """Return (odds ratio, 95% CI low, 95% CI high, p) or None if
        any margin is zero.

        OR and CI come from :func:`scipy.stats.contingency.odds_ratio`
        with ``kind='sample'`` (matches Fisher's plain ``a·d/(b·c)``).
        The CI is the exact conditional Fisher CI. p-value is Fisher's
        two-sided exact test.
        """
        if not self.is_complete:
            return None
        table = [
            [self.alt_i, self.alt_a],
            [self.nonalt_i, self.nonalt_a],
        ]
        _, p = fisher_exact(table)
        or_result = odds_ratio(table, kind="sample")
        ci = or_result.confidence_interval(confidence_level=0.95)
        return (
            float(or_result.statistic),
            float(ci.low),
            float(ci.high),
            float(p),
        )

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
        odds, ci_low, ci_high, p = fisher
        p_fmt = f"{p:.2e}" if p < 1e-3 else f"{p:.3f}"

        def _or_fmt(x: float) -> str:
            if x == float("inf"):
                return "  inf"
            if x >= 100:
                return f"{x:5.0f}"
            return f"{x:5.2f}"

        return (
            f"{head}    OR = {_or_fmt(odds)}  "
            f"[{_or_fmt(ci_low)}, {_or_fmt(ci_high)}]    p = {p_fmt}"
        )


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def _resolve_records(
    lemma: str,
    verbs: list[str],
    suffixes: list[str],
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
) -> list[tuple[str, tuple[str, ...]]]:
    """Return one (suffix, verbs-for-that-suffix) pair per emit.

    A noun with etymology-accepted verbs under one suffix returns a
    single-element list. A noun with accepted verbs under *both* ``-i``
    and ``-a`` returns two elements — the analysis then counts the
    noun once in each column, since the language demonstrably accepts
    both derivations from the same base. A USER_KEEP noun with no
    audited pair falls back to the first ACCEPTED_VERBALIZERS suffix
    in the field. Any other configuration returns an empty list (drop).
    """
    pairs = [
        (v, s) for v, s in zip(verbs, suffixes) if s in ACCEPTED_VERBALIZERS
    ]
    if not pairs:
        return []
    by_suffix: dict[str, list[str]] = {}
    for v, s in pairs:
        if accepted_by_pair.get((lemma, v), False):
            by_suffix.setdefault(s, []).append(v)
    if by_suffix:
        return [(s, tuple(vs)) for s, vs in sorted(by_suffix.items())]
    if lemma in user_keep:
        first_v, first_s = pairs[0]
        return [(first_s, (first_v,))]
    return []


def _load_records(
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
) -> tuple[list[NounRecord], dict[str, float]]:
    """Read the lexicon CSV and return the list of NounRecords that
    pass the row-side filters, plus a {lemma → frequency} map."""
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
            resolutions = _resolve_records(
                lemma, verbs, suffixes, accepted_by_pair, user_keep
            )
            if not resolutions:
                continue
            seen.add(lemma)
            ipa_row_class = _classify_ipa_row(
                row.get("ipa_normalized_lemma", ""),
                row.get("ipa_normalized_pl", ""),
                stem_final,
            )
            mutation = row["mutation"].strip().upper() == "TRUE"
            for deriv_suffix, resolved_verbs in resolutions:
                records.append(
                    NounRecord(
                        lemma=lemma,
                        derived_verbs=resolved_verbs,
                        stem_final=stem_final,
                        opportunity=opp,
                        deriv_suffix=deriv_suffix,
                        mutation=mutation,
                        frequency=freq,
                        ipa_row_class=ipa_row_class,
                    )
                )
    return records, frequency_by_lemma


# ---------------------------------------------------------------------------
# Etymology + paradigm audits
# ---------------------------------------------------------------------------


def _load_verb_verdicts() -> dict[str, str]:
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
    user_keep |= set(USER_KEEP)
    return accepted, user_keep


def _accepted_verbs(
    record: NounRecord,
    accepted_by_pair: dict[tuple[str, str], bool],
    verdicts: dict[str, str],
    col_coding: ColCoding,
    unclassified_pairs: set[tuple[str, str]],
) -> list[str]:
    """Return the verbs on this noun that count as evidence under the
    column coding.

    Under SUFFIX-ONLY, all etymology-accepted verbs pass through.
    Under ALTERNATES-ONLY, ``-i`` verbs whose paradigm verdict isn't
    ``ALTERNATES`` are dropped (their morphology bleeds the trigger,
    or there's no paradigm data — either way, no phonological
    evidence). ``-a`` verbs are always kept.
    """
    out: list[str] = []
    for verb in record.derived_verbs:
        key = (record.lemma, verb)
        if key not in accepted_by_pair:
            unclassified_pairs.add(key)
            continue
        if not accepted_by_pair[key]:
            continue
        if col_coding is ColCoding.ALTERNATES_ONLY:
            # -a verbs: automatic. -i verbs: only if verdict ALTERNATES.
            if record.deriv_suffix == "-i":
                if verdicts.get(verb) != ALTERNATES:
                    continue
        out.append(verb)
    return out


# ---------------------------------------------------------------------------
# Panel building
# ---------------------------------------------------------------------------


_ScopeFilter = Callable[[NounRecord], bool]


def _build_panels(
    records: list[NounRecord],
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
    verdicts: dict[str, str],
    row_coding: RowCoding,
    col_coding: ColCoding,
    unclassified_pairs: set[tuple[str, str]],
    scope: _ScopeFilter,
) -> tuple[dict[StemClass | str, PanelStats], int]:
    """Build dorsal / coronal / pooled panels under a specific coding.

    Returns the three panels + the count of kept records.
    """
    cells: dict[StemClass | str, Counter[tuple[bool, str]]] = {
        StemClass.DORSAL: Counter(),
        StemClass.CORONAL: Counter(),
        "pooled": Counter(),
    }
    n_kept = 0
    for record in records:
        if not scope(record):
            continue
        row_alt = record.plural_alt(row_coding)
        if row_alt is None:
            continue  # UNCLEAR IPA
        if record.lemma in user_keep:
            keep = True
        else:
            keep = bool(
                _accepted_verbs(
                    record,
                    accepted_by_pair,
                    verdicts,
                    col_coding,
                    unclassified_pairs,
                )
            )
        if not keep:
            continue
        n_kept += 1
        key = (row_alt, record.deriv_suffix)
        cells[record.stem_class][key] += 1
        cells["pooled"][key] += 1

    def pack(c: Counter[tuple[bool, str]]) -> PanelStats:
        return PanelStats(
            alt_i=c[(True, "-i")],
            alt_a=c[(True, "-a")],
            nonalt_i=c[(False, "-i")],
            nonalt_a=c[(False, "-a")],
        )

    return ({k: pack(v) for k, v in cells.items()}, n_kept)


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------


def _contamination_diagnostic(
    records: list[NounRecord],
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
    verdicts: dict[str, str],
) -> str:
    """Under the SUFFIX-ONLY column coding, what share of the ``-i``
    column has a trigger-absent paradigm verdict?

    High shares mean the SUFFIX-ONLY panel is measurably contaminated
    by morphological-block cases where the -i suffix was selected but
    palatalisation never applied. This is what motivates promoting
    the ALTERNATES-ONLY column to primary.
    """
    lines: list[str] = []
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
            for verb in _accepted_verbs(
                rec,
                accepted_by_pair,
                verdicts,
                ColCoding.SUFFIX_ONLY,
                unclassified_throwaway,
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


def _ipa_row_diagnostic(records: list[NounRecord]) -> str:
    """Show the three-state IPA-row breakdown, especially the LEX
    (already-palatal-in-singular) column that's silently pooled with
    ALT.

    Deduped by lemma so that ambiguous nouns (two records under the
    multi-suffix fix) count once toward the row-state distribution.
    """
    lines: list[str] = []
    lines.append("\n  Distribution of IPA-row states (deduped by lemma):")
    seen: set[str] = set()
    tab: Counter[tuple[str, IPARowClass]] = Counter()
    for rec in records:
        if rec.lemma in seen:
            continue
        seen.add(rec.lemma)
        tab[(rec.stem_final, rec.ipa_row_class)] += 1

    lines.append(
        f"    {'sf':3s}  {'ALT':>5s}  {'LEX':>5s}  {'NONALT':>6s}  "
        f"{'UNCLEAR':>7s}"
    )
    for sf in sorted(TARGET_STEMS):
        alt = tab.get((sf, IPARowClass.ALT), 0)
        lex = tab.get((sf, IPARowClass.LEX), 0)
        non = tab.get((sf, IPARowClass.NONALT), 0)
        unc = tab.get((sf, IPARowClass.UNCLEAR), 0)
        lines.append(f"    {sf:3s}  {alt:5d}  {lex:5d}  {non:6d}  {unc:7d}")

    # Lex breakdown: are the LEX rows selecting -i or -a? Under both
    # theories they should behave like ALT (K/T͡ʃ class). If they don't,
    # something's wrong.
    lex_verbalizer: Counter[str] = Counter()
    for rec in records:
        if rec.ipa_row_class is IPARowClass.LEX:
            lex_verbalizer[rec.deriv_suffix] += 1
    if lex_verbalizer:
        lines.append(
            "\n  LEX rows (base already palatal in singular) — "
            "verbalizer choice:"
        )
        for suf, n in sorted(lex_verbalizer.items()):
            lines.append(f"    {suf:3s}: {n}")
    return "\n".join(lines)


def _ui_sidebar(
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
) -> str:
    """Report ``-ui`` verbs that the main panels drop.

    Steriade treats ``-ui`` separately: like ``-i`` it ends in a front
    vocoid but unlike ``-i`` it does not palatalise the base. If the
    inflection-dependence hypothesis is right, ``-ui`` verbs should
    pattern with ``-a`` (nonalt bases), not with ``-i`` (alt bases).
    """
    lines: list[str] = []
    lines.append(
        "\n  ``-ui`` verbs on target-stem nouns (etymology-audited, "
        "deduped by lemma):"
    )
    lines.append(
        "  Prediction under Steriade: if -ui doesn't palatalise, it "
        "should pattern like -a — i.e. attach to nonalt (K/K) bases."
    )
    counts: Counter[tuple[str, bool]] = Counter()
    total = 0
    kept_lemmas: set[str] = set()
    with CSV_PATH.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["pos"] != "N":
                continue
            sf = row["stem_final"]
            if sf not in TARGET_STEMS:
                continue
            if row["mutation"].strip() == "":
                continue
            opp = row["opportunity"].strip()
            if opp not in PLURAL_FRONT_OPPS and opp != "uri":
                continue
            if (row.get("exception_reason") or "").startswith("nde:"):
                continue
            derived_field = (row.get("derived_verbs") or "").strip()
            if not derived_field:
                continue
            lemma = row["lemma"]
            if lemma in kept_lemmas:
                continue
            verbs = [v.strip() for v in derived_field.split("|") if v.strip()]
            suffixes = [
                s.strip()
                for s in (row.get("deriv_suffixes") or "").split("|")
                if s.strip()
            ]
            ui_verbs = [v for v, s in zip(verbs, suffixes) if s == "-ui"]
            if not ui_verbs:
                continue
            # Etymology gate: at least one accepted -ui pair, or lemma
            # is USER_KEEP.
            accepted = any(
                accepted_by_pair.get((lemma, v), False) for v in ui_verbs
            )
            if not accepted and lemma not in user_keep:
                continue
            kept_lemmas.add(lemma)
            total += 1
            mutation = row["mutation"].strip().upper() == "TRUE"
            counts[(sf, mutation)] += 1

    if total == 0:
        lines.append("    (no audited -ui pairs)")
        return "\n".join(lines)
    lines.append(
        f"\n    {'sf':3s}  {'alt-pl':>7s}  {'nonalt-pl':>10s}  {'n':>4s}"
    )
    for sf in sorted(TARGET_STEMS):
        a = counts.get((sf, True), 0)
        n = counts.get((sf, False), 0)
        if a + n == 0:
            continue
        lines.append(f"    {sf:3s}  {a:7d}  {n:10d}  {a + n:4d}")
    lines.append(f"    {'':3s}  {'':7s}  {'':10s}  {total:4d}  (total)")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


PANEL_MATRIX: Final[tuple[tuple[str, RowCoding, ColCoding, str], ...]] = (
    (
        "PRIMARY",
        RowCoding.IPA,
        ColCoding.ALTERNATES_ONLY,
        "Row: sg-vs-pl IPA (3-state, pooled). "
        "Col: -i verbs must be ALTERNATES; -a always kept.",
    ),
    (
        "MUTATION-ALT",
        RowCoding.MUTATION,
        ColCoding.ALTERNATES_ONLY,
        "Row: mutation flag (orth sg-vs-pl). "
        "Col: -i verbs must be ALTERNATES; -a always kept.",
    ),
    (
        "IPA-SUFFIX",
        RowCoding.IPA,
        ColCoding.SUFFIX_ONLY,
        "Row: sg-vs-pl IPA. Col: -i vs -a suffix, no paradigm filter.",
    ),
    (
        "LITERAL-STERIADE",
        RowCoding.MUTATION,
        ColCoding.SUFFIX_ONLY,
        "Row: mutation flag. Col: -i vs -a suffix. "
        "Literal Steriade coding — no ALTERNATES filter applied.",
    ),
)


def _render_panel(
    name: str,
    row_coding: RowCoding,
    col_coding: ColCoding,
    caption: str,
    records: list[NounRecord],
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
    verdicts: dict[str, str],
    frequency_by_lemma: dict[str, float],
    unclassified_pairs: set[tuple[str, str]],
) -> str:
    lines: list[str] = []
    lines.append("\n" + "=" * 78)
    lines.append(f"{name}")
    lines.append(caption)
    lines.append("=" * 78)

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
            row_coding,
            col_coding,
            unclassified_pairs,
            scope,
        )
        lines.append(f"\n{label}  (n_kept = {n_kept})")
        for stem in (StemClass.DORSAL, StemClass.CORONAL, "pooled"):
            stats = panels[stem]
            disp = stem.value if isinstance(stem, StemClass) else stem
            lines.append(stats.format_line(disp))
    return "\n".join(lines)


def _render_report(
    records: list[NounRecord],
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
    verdicts: dict[str, str],
    unclassified_pairs: set[tuple[str, str]],
    frequency_by_lemma: dict[str, float],
) -> str:
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

    # IPA-row diagnostic — show the LEX split so it can't hide
    lines.append("\n" + "=" * 78)
    lines.append("IPA-ROW DIAGNOSTIC: three-state sg-vs-pl breakdown")
    lines.append("=" * 78)
    lines.append(_ipa_row_diagnostic(records))

    # Four panels, PRIMARY first
    for name, row_coding, col_coding, caption in PANEL_MATRIX:
        lines.append(
            _render_panel(
                name,
                row_coding,
                col_coding,
                caption,
                records,
                accepted_by_pair,
                user_keep,
                verdicts,
                frequency_by_lemma,
                unclassified_pairs,
            )
        )

    # Contamination diagnostic under SUFFIX-ONLY column coding
    lines.append("\n" + "=" * 78)
    lines.append("CONTAMINATION DIAGNOSTIC: trigger-absent share of -i column")
    lines.append(
        "(under SUFFIX-ONLY column coding; motivates ALTERNATES-ONLY)"
    )
    lines.append("=" * 78)
    lines.append(
        _contamination_diagnostic(
            records, accepted_by_pair, user_keep, verdicts
        )
    )

    # -ui sidebar — Steriade treats -ui separately from -i/-a
    lines.append("\n" + "=" * 78)
    lines.append("SIDEBAR: -ui verbs (excluded from the main panels)")
    lines.append("=" * 78)
    lines.append(_ui_sidebar(accepted_by_pair, user_keep))
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
    OUTPUT_PATH.write_text(report, encoding="utf-8")
    print(f"\nReport persisted to {OUTPUT_PATH}")

    if unclassified_pairs:
        n = len(unclassified_pairs)
        sample = sorted(unclassified_pairs)[:10]
        print(
            f"\nWARNING: {n} (noun, verb) pairs lack an etymology "
            f"classification in {ETYMOLOGY_CSV_PATH.name}."
        )
        print("First 10 unclassified pairs:")
        for noun, verb in sample:
            print(f"  {noun}  →  {verb}")


if __name__ == "__main__":
    main()
