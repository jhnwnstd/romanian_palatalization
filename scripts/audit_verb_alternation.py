#!/usr/bin/env python3
"""Diagnostic audit of derived-verb conjugational morphology.

Reads each derived verb's DEX conjugation paradigm from the cache and
classifies it on two axes:

- AUGMENT
    Whether the present-indicative paradigm shows an augment
    (``-esc/-ești/-ește`` for 4th-conjugation, ``-ez/-ezi/-ează`` for
    1st-conjugation). The augment inserts a vowel between the stem and
    any front-vowel inflectional ending, bleeding the phonological
    trigger for assibilation/palatalization at the stem-final position.

- VERDICT
    A three-way classification of whether the stem-final alternation
    surfaces:

    - :attr:`VerbVerdict.ALTERNATES` — the palatalized counterpart of
      the noun's stem-final appears at the stem-final position of some
      diagnostic cell (2sg/3sg present indicative, 3sg present
      subjunctive).
    - :attr:`VerbVerdict.NO_ALTERNATION` — diagnostic cells exist with
      the phonological trigger present, but no cell shows the
      alternation. Empirically near-empty in our dataset.
    - :attr:`VerbVerdict.TRIGGER_ABSENT` — no diagnostic cell offers a
      bare front-vowel trigger at the stem-final position; the augment
      morphologically blocks the alternation regardless of the
      underlying segment.
    - :attr:`VerbVerdict.NO_PARADIGM` — cache hit, but the verb tree on
      the page has no parseable lexeme paradigm table.
    - :attr:`VerbVerdict.NO_PAGE` — the verb was never cached.

Purpose
-------
This audit is **diagnostic, not corrective**. The verdicts here
characterize Romanian morphology of these verbs; they do NOT justify
excluding verbs from the affix-avoidance contingency table.

Steriade's lexical counts in (10.20) tabulate the verbalizer allomorph
distribution across all K/TS and K/K bases without filtering by
whether the verb's conjugated paradigm independently surfaces the
alternation. She is explicit (p. 9, after 10.15) that *"only the former
[affix avoidance] arises with derived verbs"* — the phonotactic-
blocking manifestation applies to single-form suffixes like ``-ist``,
not to verbs.

For the abstract's canonical contingency table see
``scripts/build_contingency_table.py``.

Detector anchoring
------------------
The alternation detector is anchored to the stem-final POSITION (the
end of the stem after stripping a known inflectional ending), not
"is the palatalized counterpart anywhere in the form". This avoids
false positives from suffix-internal ``ș`` in forms like ``rostești``
where the ``ș`` comes from the ``-ești`` augment ending, not from
palatalization of the stem-``s``.
"""

from __future__ import annotations

import csv
import re
import sys
from collections import Counter
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import Final, Iterable

# Add src to path so we can import the shared cache reader.
_PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PROJECT_ROOT / "src"))

from dex_cache_reader import (  # noqa: E402
    DEFAULT_CACHE_PATH,
    get_html,
    load_cache_subset,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

CSV_PATH: Final[Path] = (
    _PROJECT_ROOT / "data" / "romanian_lexicon_with_freq.csv"
)
CACHE_PATH: Final[Path] = DEFAULT_CACHE_PATH
OUT_PATH: Final[Path] = (
    _PROJECT_ROOT / "analysis" / "verb_alternation_audit.txt"
)

DORSAL_STEMS: Final[frozenset[str]] = frozenset({"c", "g"})
CORONAL_STEMS: Final[frozenset[str]] = frozenset({"t", "d", "s"})
TARGET_STEMS: Final[frozenset[str]] = DORSAL_STEMS | CORONAL_STEMS

# Front-vowel triggers in Romanian orthography. Coronal assibilation
# applies before /i/ only; dorsal palatalization applies before /i/ or
# /e/.
FRONT_VOWELS: Final[str] = "ei"

# Orthographic palatalized counterparts of the coronal stem-finals.
PALATAL_COUNTERPART: Final[dict[str, str]] = {
    "t": "ț",
    "d": "z",
    "s": "ș",
}


# ---------------------------------------------------------------------------
# Enums and result types
# ---------------------------------------------------------------------------


class VerbVerdict(StrEnum):
    """Three-way verdict on whether the verb shows the alternation."""

    ALTERNATES = "ALTERNATES"
    NO_ALTERNATION = "NO-ALTERNATION"
    TRIGGER_ABSENT = "TRIGGER-ABSENT"
    NO_PARADIGM = "NO-PARADIGM"
    NO_PAGE = "NO-PAGE"


class PresentColumn(StrEnum):
    """Columns of the present-tense section of a DEX paradigm table."""

    INDICATIVE = "prez_ind"
    SUBJUNCTIVE = "prez_subj"
    IMPERFECT = "imperfect"
    PERF_SIMPLU = "perf_simplu"
    MMCP = "mmcp"


# DEX uses these literal labels for present-tense person rows.
class PersonRow(StrEnum):
    SG_1 = "I (eu)"
    SG_2 = "a II-a (tu)"
    SG_3 = "a III-a (el, ea)"
    PL_1 = "I (noi)"
    PL_2 = "a II-a (voi)"
    PL_3 = "a III-a (ei, ele)"


DIAGNOSTIC_PERSONS: Final[tuple[PersonRow, ...]] = (
    PersonRow.SG_2,
    PersonRow.SG_3,
)
AUGMENT_PROBE_PERSONS: Final[tuple[PersonRow, ...]] = (
    PersonRow.SG_1,
    PersonRow.SG_2,
    PersonRow.SG_3,
    PersonRow.PL_3,
)

# Ordered columns of the present-tense section, as they appear in the
# DEX table. Used to map raw form cells to named columns.
_PRESENT_COLUMN_ORDER: Final[tuple[PresentColumn, ...]] = (
    PresentColumn.INDICATIVE,
    PresentColumn.SUBJUNCTIVE,
    PresentColumn.IMPERFECT,
    PresentColumn.PERF_SIMPLU,
    PresentColumn.MMCP,
)


# A paradigm is {person → {column → [surface forms]}}.
Paradigm = dict[PersonRow, dict[PresentColumn, list[str]]]


@dataclass(frozen=True, slots=True)
class VerbAudit:
    """Per-verb audit result.

    ``verdict`` is the headline classification. The remaining fields
    are kept around so the saved audit file is self-describing.
    """

    verdict: VerbVerdict
    augment: bool | None
    conjugation_class: str  # "IV", "II", etc. — DEX's Roman-numeral tag
    paradigm: Paradigm = field(default_factory=dict)


# ---------------------------------------------------------------------------
# DEX paradigm-table parser
# ---------------------------------------------------------------------------

_LEXEME_TABLE_RE: Final[re.Pattern[str]] = re.compile(
    r'<table[^>]*class="[^"]*lexeme[^"]*"[^>]*>(.*?)</table>',
    re.DOTALL | re.IGNORECASE,
)
_PERSON_ROW_RE: Final[re.Pattern[str]] = re.compile(
    r'<td[^>]*class="[^"]*inflection person[^"]*"[^>]*>' r"([^<]+)</td>",
    re.IGNORECASE,
)
_FORM_CELL_RE: Final[re.Pattern[str]] = re.compile(
    r'<td[^>]*class="[^"]*form[^"]*"[^>]*>(.*?)</td>',
    re.DOTALL | re.IGNORECASE,
)
_TONIC_SPAN_RE: Final[re.Pattern[str]] = re.compile(
    r'<span[^>]*class="[^"]*tonic-accent[^"]*"[^>]*>([^<]+)</span>',
    re.IGNORECASE,
)
_ANY_TAG_RE: Final[re.Pattern[str]] = re.compile(r"<[^>]+>")
_SE_PREFIX_RE: Final[re.Pattern[str]] = re.compile(r"^\s*\(s[ăa]\)\s*")


def _clean_form(td_content: str) -> str:
    """Reduce a ``<td class="form">`` cell to a single surface form."""
    x = _TONIC_SPAN_RE.sub(r"\1", td_content)
    x = _ANY_TAG_RE.sub(" ", x)
    x = x.replace("&#x2011", "").replace("‑", "")
    x = _SE_PREFIX_RE.sub("", x)
    parts = [p.strip() for p in re.split(r"\s{2,}|\n", x) if p.strip()]
    return parts[0] if parts else x.strip()


def _all_forms(td_content: str) -> list[str]:
    """Return all forms in a cell — DEX often lists historical variants."""
    x = _TONIC_SPAN_RE.sub(r"\1", td_content)
    lis = re.findall(r"<li[^>]*>(.*?)</li>", x, re.DOTALL)
    if not lis:
        return [_clean_form(td_content)]
    out: list[str] = []
    for li in lis:
        s = (
            _ANY_TAG_RE.sub("", li)
            .replace("&#x2011", "")
            .replace("‑", "")
            .strip()
        )
        if s:
            out.append(s)
    return out


def _extract_present_paradigm(html: str) -> Paradigm:
    """Parse the DEX lexeme table and return the present-tense block."""
    tables = _LEXEME_TABLE_RE.findall(html)
    if not tables:
        return {}
    table = tables[0]
    persons = list(_PERSON_ROW_RE.finditer(table))
    if not persons:
        return {}
    paradigm: Paradigm = {}
    for i, match in enumerate(persons):
        label = match.group(1).strip()
        try:
            person = PersonRow(label)
        except ValueError:
            continue
        end = persons[i + 1].start() if i + 1 < len(persons) else len(table)
        row_html = table[match.end() : end]
        cells = _FORM_CELL_RE.findall(row_html)
        row: dict[PresentColumn, list[str]] = {}
        for col_i, col_name in enumerate(_PRESENT_COLUMN_ORDER):
            if col_i < len(cells):
                row[col_name] = _all_forms(cells[col_i])
        paradigm[person] = row
    return paradigm


def _conjugation_class(html: str) -> str:
    """DEX prints ``conjugarea a X-a`` near the lexeme header."""
    match = re.search(r"conjugarea a ([IVX]+)-a", html)
    return match.group(1) if match else ""


# ---------------------------------------------------------------------------
# Augment + alternation detection
# ---------------------------------------------------------------------------

_AUGMENT_PATTERNS: Final[tuple[re.Pattern[str], ...]] = (
    re.compile(r"(?:esc|ești|ește)\b", re.IGNORECASE),
    re.compile(r"(?:ez|ezi|ează)\b", re.IGNORECASE),
)


# Inflectional endings we know how to strip. Longest first so a longer
# ending is preferred when there's overlap. This list captures the
# Romanian present-tense + a handful of other common forms; we don't
# need a complete morphological analyzer — we just need to isolate the
# stem-final position for the diagnostic cells.
_INFLECTION_ENDINGS: Final[tuple[str, ...]] = (
    "ească",
    "ești",
    "ește",
    "ează",
    "ezi",
    "eze",
    "esc",
    "ind",
    "ăm",
    "ați",
    "âm",
    "are",
    "ere",
    "ire",
    "ăsc",
    "i",
    "e",
    "a",
    "ă",
    "u",
)


def _strip_ending(form: str) -> str:
    """Strip the longest matching inflectional ending."""
    for ending in _INFLECTION_ENDINGS:
        if form.endswith(ending) and len(form) > len(ending):
            return form[: -len(ending)]
    return form


def _has_augment(paradigm: Paradigm) -> bool:
    """True iff any present-paradigm cell carries an augment ending.

    Paradigm-wide rather than 1sg-only: an irregular or syncretic 1sg
    can mislead a single-cell classifier, but the augment is a property
    of the inflectional class and shows up in multiple cells.
    """
    for person in AUGMENT_PROBE_PERSONS:
        row = paradigm.get(person, {})
        for column in (
            PresentColumn.INDICATIVE,
            PresentColumn.SUBJUNCTIVE,
        ):
            for form in row.get(column, []):
                for pattern in _AUGMENT_PATTERNS:
                    if pattern.search(form):
                        return True
    return False


def _diagnostic_forms(paradigm: Paradigm) -> Iterable[str]:
    """Yield every form in the cells where the alternation is testable."""
    for person in DIAGNOSTIC_PERSONS:
        row = paradigm.get(person, {})
        yield from row.get(PresentColumn.INDICATIVE, [])
        yield from row.get(PresentColumn.SUBJUNCTIVE, [])
    # 1sg present indicative shows the augment but rarely the
    # alternation; we include it so the detector sees augment material.
    yield from paradigm.get(PersonRow.SG_1, {}).get(
        PresentColumn.INDICATIVE, []
    )


def _stem_alternates(stem_final: str, form: str) -> bool:
    """Does ``form`` palatalize the stem-final at the stem position?"""
    if not form:
        return False
    sf = stem_final.lower()
    stem = _strip_ending(form.lower())
    if not stem:
        return False
    if sf in PALATAL_COUNTERPART:
        return stem.endswith(PALATAL_COUNTERPART[sf])
    if sf in DORSAL_STEMS:
        ending = form.lower()[len(stem) :]
        if stem.endswith(sf) and ending and ending[0] in FRONT_VOWELS:
            return True
        for i, ch in enumerate(stem[:-1]):
            if ch == sf and stem[i + 1] in FRONT_VOWELS:
                return True
        return False
    return False


def _trigger_present(stem_final: str, paradigm: Paradigm) -> bool:
    """Does any diagnostic form put the stem-final before a trigger?"""
    sf = stem_final.lower()
    for form in _diagnostic_forms(paradigm):
        f = form.lower()
        stem = _strip_ending(f)
        if not stem:
            continue
        ending = f[len(stem) :]
        if sf in PALATAL_COUNTERPART:
            counterpart = PALATAL_COUNTERPART[sf]
            if (
                (stem.endswith(sf) or stem.endswith(counterpart))
                and ending
                and ending[0] == "i"
            ):
                return True
        elif sf in DORSAL_STEMS:
            if stem.endswith(sf) and ending and ending[0] in FRONT_VOWELS:
                return True
            for i, ch in enumerate(stem[:-1]):
                if ch == sf and stem[i + 1] in FRONT_VOWELS:
                    return True
    return False


def classify_verb(html: str, stem_final: str) -> VerbAudit:
    """Score a derived verb against its conjugation paradigm.

    See module docstring for the verdict semantics. ``html`` is the
    cached DEX page for the verb's infinitive; ``stem_final`` is the
    noun's stem-final consonant (so we know which palatalized
    counterpart to look for).
    """
    paradigm = _extract_present_paradigm(html)
    if not paradigm:
        return VerbAudit(
            verdict=VerbVerdict.NO_PARADIGM,
            augment=None,
            conjugation_class=_conjugation_class(html),
        )
    alternates = any(
        _stem_alternates(stem_final, form)
        for form in _diagnostic_forms(paradigm)
    )
    if alternates:
        verdict = VerbVerdict.ALTERNATES
    elif not _trigger_present(stem_final, paradigm):
        verdict = VerbVerdict.TRIGGER_ABSENT
    elif _has_augment(paradigm):
        # If we got here, _trigger_present saw something — but an
        # augmented verb still effectively bleeds the alternation,
        # so we report TRIGGER_ABSENT to keep the category clean.
        verdict = VerbVerdict.TRIGGER_ABSENT
    else:
        verdict = VerbVerdict.NO_ALTERNATION

    return VerbAudit(
        verdict=verdict,
        augment=_has_augment(paradigm),
        conjugation_class=_conjugation_class(html),
        paradigm=paradigm,
    )


# ---------------------------------------------------------------------------
# Validation harness
# ---------------------------------------------------------------------------

# (verb, stem_final, expected verdict, prose justification)
_VALIDATION_ANCHORS: Final[tuple[tuple[str, str, VerbVerdict, str], ...]] = (
    (
        "plăti",
        "t",
        VerbVerdict.TRIGGER_ABSENT,
        "augmented -esc; t never meets bare /i/",
    ),
    (
        "lupta",
        "t",
        VerbVerdict.ALTERNATES,
        "2sg lupți shows t→ț",
    ),
    (
        "munci",
        "c",
        VerbVerdict.ALTERNATES,
        "augmented but dorsal: c+e is palatalized too",
    ),
    (
        "fabrica",
        "c",
        VerbVerdict.ALTERNATES,
        "2sg fabrici shows c+i palatalization",
    ),
    (
        "ataca",
        "c",
        VerbVerdict.ALTERNATES,
        "2sg ataci shows c+i palatalization",
    ),
    (
        "scuza",
        "z",
        VerbVerdict.TRIGGER_ABSENT,
        "stem already z; no further alternation",
    ),
)


def _run_validation(cache: dict[str, str]) -> None:
    print("=" * 70, file=sys.stderr)
    print("VALIDATION", file=sys.stderr)
    print("=" * 70, file=sys.stderr)
    for verb, sf, expected, note in _VALIDATION_ANCHORS:
        html = get_html(cache, verb)
        if html is None:
            print(f"  {verb:14s}  NOT IN CACHE", file=sys.stderr)
            continue
        result = classify_verb(html, sf)
        status = "OK" if result.verdict == expected else "❌"
        print(
            f"  {status} {verb:14s} stem={sf}  "
            f"verdict={result.verdict:14s}  "
            f"expected={expected}  ({note})",
            file=sys.stderr,
        )


# ---------------------------------------------------------------------------
# Audit row loading
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class AuditRow:
    noun: str
    verb: str
    stem_final: str
    opportunity: str
    deriv_suffix: str


def _load_audit_rows() -> list[AuditRow]:
    """Same row filters as the contingency table; one row per (noun,verb)."""
    out: list[AuditRow] = []
    with CSV_PATH.open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["pos"] != "N":
                continue
            sf = r["stem_final"]
            if sf not in TARGET_STEMS:
                continue
            if r["mutation"].strip() == "":
                continue
            opp = r["opportunity"].strip()
            if opp not in {"i", "e", "uri"}:
                continue
            deriv = r["deriv_suffixes"].strip()
            if deriv not in {"-i", "-a"}:
                continue
            dv = (r.get("derived_verbs") or "").strip()
            if not dv:
                continue
            if (r.get("exception_reason") or "").startswith("nde:"):
                continue
            for verb in (v.strip() for v in dv.split("|") if v.strip()):
                out.append(
                    AuditRow(
                        noun=r["lemma"],
                        verb=verb,
                        stem_final=sf,
                        opportunity=opp,
                        deriv_suffix=deriv,
                    )
                )
    return out


# ---------------------------------------------------------------------------
# Reporting and persistence
# ---------------------------------------------------------------------------


def _stem_class(stem_final: str) -> str:
    return "dorsal" if stem_final in DORSAL_STEMS else "coronal"


def _report_distribution(
    results: list[tuple[AuditRow, VerbAudit]],
) -> None:
    counts: Counter[tuple[str, str, str]] = Counter()
    for row, result in results:
        counts[
            (
                _stem_class(row.stem_final),
                row.deriv_suffix,
                str(result.verdict),
            )
        ] += 1
    print("\nVerdict distribution by stem-class × verbalizer:")
    for seg in ("dorsal", "coronal"):
        print(f"\n  {seg.upper()}:")
        for verbalizer in ("-i", "-a"):
            print(f"    {verbalizer}:")
            for verdict in (
                VerbVerdict.ALTERNATES,
                VerbVerdict.NO_ALTERNATION,
                VerbVerdict.TRIGGER_ABSENT,
                VerbVerdict.NO_PARADIGM,
                VerbVerdict.NO_PAGE,
            ):
                n = counts[(seg, verbalizer, verdict.value)]
                if n:
                    print(f"      {verdict.value:18s}: {n:4d}")


def _write_audit_file(
    results: list[tuple[AuditRow, VerbAudit]],
) -> None:
    OUT_PATH.parent.mkdir(exist_ok=True, parents=True)
    with OUT_PATH.open("w", encoding="utf-8") as f:
        f.write("# Verb-alternation column audit\n")
        f.write(
            "# noun\tderived_verb\tstem_final\topp\tds\tverdict"
            "\taugment\tconjugare\t2sg.pres\t3sg.pres\t3sg.subj\n\n"
        )
        for row, result in sorted(
            results,
            key=lambda pair: (
                pair[0].stem_final,
                pair[0].deriv_suffix,
                pair[0].noun,
            ),
        ):
            paradigm = result.paradigm
            sg2 = "|".join(
                paradigm.get(PersonRow.SG_2, {}).get(
                    PresentColumn.INDICATIVE, []
                )
            )
            sg3 = "|".join(
                paradigm.get(PersonRow.SG_3, {}).get(
                    PresentColumn.INDICATIVE, []
                )
            )
            sj3 = "|".join(
                paradigm.get(PersonRow.SG_3, {}).get(
                    PresentColumn.SUBJUNCTIVE, []
                )
            )
            f.write(
                f"{row.noun}\t{row.verb}\t{row.stem_final}\t"
                f"{row.opportunity}\t{row.deriv_suffix}\t"
                f"{result.verdict}\t{result.augment}\t"
                f"{result.conjugation_class}\t{sg2}\t{sg3}\t{sj3}\n"
            )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    validation_verbs = {anchor[0] for anchor in _VALIDATION_ANCHORS}
    audit_rows = _load_audit_rows()
    print(
        f"\n(noun, verb, sf) tuples: {len(audit_rows)}",
        file=sys.stderr,
    )

    unique_verbs = {row.verb for row in audit_rows} | validation_verbs
    print(
        f"  Unique verbs to look up: {len(unique_verbs)}",
        file=sys.stderr,
    )
    cache = load_cache_subset(unique_verbs, cache_path=CACHE_PATH)
    print(
        f"  Pages matched in cache: {len(cache)}",
        file=sys.stderr,
    )

    _run_validation(cache)

    print("\n" + "=" * 70, file=sys.stderr)
    print("FULL AUDIT", file=sys.stderr)
    print("=" * 70, file=sys.stderr)
    results: list[tuple[AuditRow, VerbAudit]] = []
    for row in audit_rows:
        html = get_html(cache, row.verb)
        if html is None:
            results.append(
                (
                    row,
                    VerbAudit(
                        verdict=VerbVerdict.NO_PAGE,
                        augment=None,
                        conjugation_class="",
                    ),
                )
            )
            continue
        results.append((row, classify_verb(html, row.stem_final)))

    _report_distribution(results)
    _write_audit_file(results)
    print(f"\nFull audit written to {OUT_PATH}", file=sys.stderr)


if __name__ == "__main__":
    main()
