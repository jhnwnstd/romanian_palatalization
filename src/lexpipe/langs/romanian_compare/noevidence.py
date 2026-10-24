"""No-evidence roots and derivative palatalization (R21).

The plural-less prediction: a root the learner never hears alternate
is stored prespecified, so its derivatives should not palatalize.
"No alternation evidence" is decided from the corpus, not from the
dictionary paradigm (R21): the base class of a root at τ comes from
its plural cell (nominal inflection, E0):

    alternating          the plural is attested at τ and shows the
                         palatal reflex;
    non-alternating      the plural is attested at τ before one of the
                         class's triggers (dorsal: -i, -e; coronal, z:
                         -i) and keeps the plain consonant;
    non-alternating -uri the plural is attested at τ with -uri (the
                         spell-out account's sign of prespecification);
    no evidence          anything else: no plural, a plural below τ, a
                         coronal -e plural (not a trigger), a plural the
                         legacy alignment could not parse.

The derivative side is coded at the base × suffix level (R23): every
derived verb or adjective of the lexicon that is attested at τ and
whose suffix starts with a front vowel gives ``pal`` (reflex) or
``nonpal`` (plain); a base with both outcomes for one suffix is
``both``; derivatives with a back-vowel suffix (-a, -ui, -os, ...) are
``avoid``. The 3×2 tables collapse the two non-alternating rows (as
``analysis/learnability_checks/check_claims.py`` does) and are given
per suffix group, per class and per base stratum (R22).

Candidates (R21): every root's -i verb is generated in both spellings
(foc → foci / fochi; vad → vazi / vadi) and looked up among all
derived-verb headwords of the lexicon and in the corpus. Only a
corpus hit that is not a lexicon derivative is new, and such hits are
reported separately because a corpus form can be a homograph (foci
is also a plural of focar/focus).
"""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from lexpipe.langs.romanian_compare import morphology as mo
from lexpipe.langs.romanian_compare.lexicon import STRATA, fold_key
from lexpipe.langs.romanian_compare.records import Derivative, RootRecord

Attest = Callable[[str], Mapping[str, bool]]
BASE_CLASSES: Final[tuple[str, ...]] = (
    "alternating",
    "non-alternating",
    "non-alternating -uri",
    "no evidence",
)
OUTCOMES: Final[tuple[str, ...]] = ("pal", "nonpal", "both", "avoid")
SUFFIX_GROUPS: Final[tuple[str, ...]] = (
    "verb -i",
    "verb -iza/-ifica",
    "adjective, front suffix",
    "all front suffixes",
)
CLASSES: Final[Mapping[str, tuple[str, ...]]] = {
    "dorsal": ("dorsal",),
    "coronal": ("coronal",),
    "z": ("z",),
}
REPORT_TAUS: Final[tuple[str, ...]] = ("dictionary", "ge1", "ge2pm")
_LEXEME: Final[str] = "#lexeme"


def base_class(record: RootRecord, tau: str, attest: Attest) -> str:
    """The root's base class at τ (module docstring)."""
    if not record.plural or not attest(fold_key(record.plural)).get(
        tau, False
    ):
        return "no evidence"
    if record.plural_allomorph == "-uri":
        return "non-alternating -uri"
    vowel = {"-i": "i", "-e": "e"}.get(record.plural_allomorph)
    if vowel is None:
        return "no evidence"
    if record.plural_palatal is True:
        return "alternating"
    if mo.triggers_alternation(record.seg, vowel):
        return "non-alternating"
    return "no evidence"


def suffix_group(derivative: Derivative) -> str | None:
    """Which front-suffix group a derivative belongs to (None: back)."""
    if derivative.vowel not in ("i", "e"):
        return None
    if derivative.kind == "verb":
        if derivative.suffix == "-i":
            return "verb -i"
        if derivative.suffix in ("-iza", "-ifica"):
            return "verb -iza/-ifica"
        return None
    return "adjective, front suffix"


def _attested(derivative: Derivative, tau: str, attest: Attest) -> bool:
    key = fold_key(derivative.form) + (
        _LEXEME if derivative.kind == "verb" else ""
    )
    return attest(key).get(tau, False)


@dataclass(frozen=True, slots=True)
class BaseSuffix:
    """One base × suffix unit."""

    root_id: str
    lemma: str
    target_class: str
    stratum: str
    base: str
    group: str
    suffix: str
    outcome: str


def units(
    records: Sequence[RootRecord], tau: str, attest: Attest
) -> list[BaseSuffix]:
    """Base × suffix units of target roots with an attested derivative."""
    out: list[BaseSuffix] = []
    for record in records:
        if not record.is_target or record.nde:
            continue
        base = base_class(record, tau, attest)
        by_suffix: dict[tuple[str, str], list[bool]] = defaultdict(list)
        avoid: set[str] = set()
        for derivative in record.derivatives:
            if derivative.vowel is None or not _attested(
                derivative, tau, attest
            ):
                continue
            group = suffix_group(derivative)
            label = derivative.suffix or derivative.kind
            if group is None:
                avoid.add(label)
                continue
            by_suffix[(group, label)].append(bool(derivative.palatal))
        for (group, label), outcomes in sorted(by_suffix.items()):
            kinds = set(outcomes)
            outcome = (
                "both"
                if len(kinds) == 2
                else ("pal" if True in kinds else "nonpal")
            )
            out.append(
                BaseSuffix(
                    record.root_id,
                    record.lemma,
                    record.target_class,
                    record.stratum,
                    base,
                    group,
                    label,
                    outcome,
                )
            )
        for label in sorted(avoid):
            out.append(
                BaseSuffix(
                    record.root_id,
                    record.lemma,
                    record.target_class,
                    record.stratum,
                    base,
                    "back suffix",
                    label,
                    "avoid",
                )
            )
    return out


def table(items: Iterable[BaseSuffix]) -> dict[str, dict[str, Counter[str]]]:
    """Counts by base class → suffix group → outcome.

    Front-suffix units count under their group and under
    ``all front suffixes``; back-suffix units under ``avoid``.
    """
    out: dict[str, dict[str, Counter[str]]] = {
        base: {group: Counter() for group in (*SUFFIX_GROUPS, "back suffix")}
        for base in BASE_CLASSES
    }
    for item in items:
        if item.outcome == "avoid":
            out[item.base]["back suffix"]["avoid"] += 1
            continue
        out[item.base][item.group][item.outcome] += 1
        out[item.base]["all front suffixes"][item.outcome] += 1
    return out


# ---------------------------------------------------------------------------
# Candidates in both spellings
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class CandidateHit:
    """A generated -i verb found as a headword or in the corpus."""

    lemma: str
    base: str
    form: str
    palatal: bool
    in_lexicon: bool


def nominal_forms(records: Iterable[RootRecord]) -> frozenset[str]:
    """Every noun lemma and plural of the lexicon (lowercased)."""
    forms: set[str] = set()
    for record in records:
        forms.add(record.lemma.lower())
        if record.plural:
            forms.add(record.plural.lower())
    return frozenset(forms)


def candidate_hits(
    records: Sequence[RootRecord],
    headwords: frozenset[str],
    tau: str,
    attest: Attest,
) -> list[CandidateHit]:
    """Generated candidates found as a derived-verb headword of the
    lexicon or with at least one corpus token.

    Skipped: the root's own derivatives, and any form that is a noun
    lemma or plural of the lexicon (coți is cot's plural, not a verb):
    the corpus count of such a form says nothing about a verb.
    """
    nominal = nominal_forms(records)
    hits: list[CandidateHit] = []
    for record in records:
        if not record.is_target or record.nde:
            continue
        own = {d.form for d in record.derivatives}
        base = base_class(record, tau, attest)
        for candidate in record.candidates:
            form = candidate.form
            if form in own or form in nominal:
                continue
            in_lexicon = form in headwords
            if in_lexicon or candidate.count >= 1:
                hits.append(
                    CandidateHit(
                        record.lemma, base, form, candidate.palatal, in_lexicon
                    )
                )
    return hits


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


_GROUP_LABELS: Final[tuple[tuple[str, str], ...]] = (
    ("verb -i", "verb -i"),
    ("verb -iza/-ifica", "-iza/-ifica"),
    ("adjective, front suffix", "adj. front"),
    ("all front suffixes", "all front"),
)


def _table_lines(tab: Mapping[str, Mapping[str, Counter[str]]]) -> list[str]:
    """One row per base class: pal/nonpal per group, both, avoid."""
    head = "".join(f"{label:>14s}" for _, label in _GROUP_LABELS)
    lines = [f"      {'base':22s}{head} {'both':>5s} {'avoid':>6s}"]
    sub = "".join(f"{'pal/nonpal':>14s}" for _ in _GROUP_LABELS)
    lines.append(f"      {'':22s}{sub}")
    pooled: dict[str, Counter[str]] = {g: Counter() for g, _ in _GROUP_LABELS}
    pooled["back suffix"] = Counter()
    rows: list[tuple[str, Mapping[str, Counter[str]]]] = [
        (base, tab[base]) for base in BASE_CLASSES
    ]
    for base in ("non-alternating", "non-alternating -uri"):
        for group, counter in tab[base].items():
            pooled[group].update(counter)
    rows.append(("(non-alt. pooled)", pooled))
    for name, groups in rows:
        cells = "".join(
            f"{groups[g]['pal']:>7d}/{groups[g]['nonpal']:<6d}"
            for g, _ in _GROUP_LABELS
        )
        both = groups["all front suffixes"]["both"]
        avoid = groups["back suffix"]["avoid"]
        lines.append(f"      {name:22s}{cells} {both:>5d} {avoid:>6d}")
    return lines


def legacy_prediction_lines(path: Path) -> list[str]:
    """The Prediction 1 block of check_claims_output.txt, verbatim."""
    if not path.is_file():
        return [f"  ({path.name} not found)"]
    text = path.read_text(encoding="utf-8").splitlines()
    out: list[str] = []
    keep = False
    for line in text:
        if line.startswith("PREDICTION 1"):
            keep = True
            continue
        if keep and line.startswith("TOP-1000"):
            break
        if keep and not re.fullmatch(r"=+", line.strip()):
            out.append("  | " + line)
    while out and out[-1].strip() in ("|", ""):
        out.pop()
    return out


def render_report(
    records: Sequence[RootRecord],
    attest: Attest,
    check_claims_output: Path,
) -> str:
    """The text of ``analysis/comparison/romanian_noevidence.txt``."""
    lines: list[str] = []
    add = lines.append
    add("Romanian no-evidence roots × derivative palatalization (design R21)")
    add("=" * 72)
    add("")
    add("Base class from the corpus (plural cell, E0) at tau: alternating /")
    add(
        "non-alternating (attested plural before a trigger, plain consonant) /"
    )
    add("non-alternating -uri / no evidence (no plural, plural below tau,")
    add("coronal -e plural, unparsed plural). Derivatives: lexicon verbs and")
    add("adjectives attested at tau (verbs: infinitive + participle tokens),")
    add("coded per base x suffix: pal / nonpal / both (front-vowel suffix),")
    add("avoid (back-vowel suffix). Structural NDE roots are excluded.")
    add("")
    headwords = frozenset(
        d.form for r in records for d in r.derivatives if d.kind == "verb"
    )
    for tau in REPORT_TAUS:
        items = units(records, tau, attest)
        add(f"tau = {tau}")
        add("-" * 72)
        bases = Counter(
            base_class(r, tau, attest)
            for r in records
            if r.is_target and not r.nde
        )
        add(
            "  target roots by base class: "
            + ", ".join(f"{b} {bases[b]}" for b in BASE_CLASSES)
        )
        for name, classes in CLASSES.items():
            add(f"    class {name} (base x suffix units)")
            lines.extend(
                _table_lines(
                    table(i for i in items if i.target_class in classes)
                )
            )
        add("  By base stratum (all classes):")
        for stratum in STRATA:
            add(f"    stratum {stratum}")
            lines.extend(
                _table_lines(table(i for i in items if i.stratum == stratum))
            )
        add("")
    hits = candidate_hits(records, headwords, "dictionary", attest)
    in_lex = [h for h in hits if h.in_lexicon]
    corpus_only = [h for h in hits if not h.in_lexicon]
    add(
        "Generated -i verb candidates (both spellings), not listed for the root"
    )
    add("-" * 72)
    add(
        f"  found among lexicon derived-verb headwords: {len(in_lex)}; "
        f"corpus-only (>= 1 token): {len(corpus_only)}"
    )
    add("  Forms equal to a noun lemma or plural of the lexicon are skipped.")
    add("  Matches are unverified: a headword can be an unrelated verb (minți")
    add("  'to lie' for mintă 'mint') and a corpus form any homograph (a 2sg,")
    add("  an adjective plural). They are listed for checking, not counted in")
    add("  the tables above.")
    by_base = Counter((h.base, h.palatal, h.in_lexicon) for h in hits)
    for base in BASE_CLASSES:
        head_pal = by_base[(base, True, True)]
        head_non = by_base[(base, False, True)]
        corp_pal = by_base[(base, True, False)]
        corp_non = by_base[(base, False, False)]
        add(
            f"    {base:22s} headword pal {head_pal:>4d} nonpal {head_non:>4d}"
            f" | corpus-only pal {corp_pal:>4d} nonpal {corp_non:>4d}"
        )
    no_ev = sorted(
        f"{h.lemma}→{h.form}" for h in in_lex if h.base == "no evidence"
    )
    if no_ev:
        add(
            "  headword candidates of no-evidence roots: "
            + ", ".join(no_ev[:40])
        )
    add("")
    add("Comparison with analysis/learnability_checks/check_claims_output.txt")
    add("-" * 72)
    add("  check_claims' Prediction 1 (quoted below) differs by definition:")
    add(
        "  its 'no plural (confirmed)' roots are lexicon nouns without a plural"
    )
    add(
        "  that DEX or Wiktionary confirms as uncountable (21 of 75 candidates);"
    )
    add(
        "  its unit is the verb (pair), its derivatives are verbs only, and its"
    )
    add(
        "  'non-alternating plural' pools -uri, non-alternating -i/-e and coronal"
    )
    add("  -e plurals. Here 'no evidence' is corpus-based (no attested plural")
    add(
        "  form showing or resisting the trigger), units are base x suffix, and"
    )
    add(
        "  adjectives are included; DEX is not consulted (R1 import isolation)."
    )
    lines.extend(legacy_prediction_lines(check_claims_output))
    return "\n".join(lines) + "\n"
