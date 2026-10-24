"""Spell-out TP on certain-specification roots (R20).

Gorman & Winstead's spell-out rules give the non-palatalizing
verbalizer (-a; -ui shields the root with /u/) only to roots
prespecified as inalterable, and -i elsewhere. A root seen
alternating in its plural is underspecified, so the account predicts
-i for it. The spell-out TP counts:

    roots   alternating-plural roots: the plural (-i, or -e for
            dorsals and sc) is attested at τ and shows the target's
            palatal reflex;
    N       those with an attested derived verb from the competing
            pair -i vs -a/-ui (lexicon ``derived_verbs``;
            infinitive + past participle tokens at τ);
    e       roots that take only -a/-ui;
            roots with both are doublets, coded by corpus share (R19).

It is computed with :func:`lexpipe.tp.two_outcome_rows`: in the cell
``verbalizer`` every derived verb of the pair is a ``P`` observation
whose ``alternates`` means "chose -i" (so ``palatalize`` is the
spell-out default, ``retain`` the -a/-ui rule), stratified by base
stratum (R22) and in the R24 V = 1000 frame.

The critique's figures (60 non-palatalizing pairs vs θ(164) = 32.2 on
the full lexicon; 17 of 42 vs θ(42) = 11.2 in a top-1000) come from
``analysis/learnability_checks/check_claims.py``: they count verbs
(pairs), not roots; their n includes every suffix the lexicon lists
(-ări, -iza, -ăni, -arisi as well as -i/-a/-ui); they cover dorsals
only; and they use no attestation threshold. :func:`critique_counts`
reproduces them from the same columns (its own plural-evidence
classes and its own top-1000 frame) so the report can show both.
"""

from __future__ import annotations

import math
import re
from collections import Counter
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Final

from lexpipe import tp
from lexpipe.attest import TAUS
from lexpipe.langs.romanian_compare import morphology as mo
from lexpipe.langs.romanian_compare.harmonized import (
    LANG,
    finish,
    frame,
)
from lexpipe.langs.romanian_compare.lexicon import STRATA, fold_key
from lexpipe.langs.romanian_compare.records import Derivative, RootRecord

TRIGGER: Final[str] = "verbalizer"
TABLE: Final[str] = "spellout"
PAIR: Final[tuple[str, ...]] = (
    mo.PALATALIZING_VERBALIZERS + mo.NONPALATALIZING_VERBALIZERS
)
CLASSES: Final[Mapping[str, tuple[str, ...]]] = {
    "dorsal": ("dorsal",),
    "coronal": ("coronal",),
    "z": ("z",),
    "all": ("dorsal", "coronal", "z"),
}
Attest = Callable[[str], Mapping[str, bool]]
_LEXEME: Final[str] = "#lexeme"


def pair_verbs(record: RootRecord) -> list[Derivative]:
    """The root's derived verbs from the competing pair -i vs -a/-ui."""
    return [
        d for d in record.derivatives if d.kind == "verb" and d.suffix in PAIR
    ]


def plural_attested(record: RootRecord, tau: str, attest: Attest) -> bool:
    return bool(record.plural) and attest(fold_key(record.plural)).get(
        tau, False
    )


def certain(record: RootRecord, tau: str, attest: Attest) -> bool:
    """Alternating-plural root (underspecified with certainty) at τ."""
    return (
        record.is_target
        and not record.nde
        and record.plural_palatal is True
        and plural_attested(record, tau, attest)
    )


def prespecified(record: RootRecord, tau: str, attest: Attest) -> bool:
    """Contrast set: an attested plural that avoids or resists the
    trigger (-uri, or a non-alternating plural before a trigger)."""
    if not record.is_target or record.nde:
        return False
    if not plural_attested(record, tau, attest):
        return False
    if record.plural_allomorph == "-uri":
        return True
    vowel = {"-i": "i", "-e": "e"}.get(record.plural_allomorph)
    return (
        vowel is not None
        and mo.triggers_alternation(record.seg, vowel)
        and record.plural_palatal is False
    )


def spellout_root(record: RootRecord, attest: Attest) -> tp.RootObservation:
    """The root with one ``verbalizer`` cell (``alternates`` = -i)."""
    observations = tuple(
        tp.Observation(
            alternates=d.suffix in mo.PALATALIZING_VERBALIZERS,
            attested_by_tau=attest(fold_key(d.form) + _LEXEME),
            allomorph="P",
            form=d.form,
            count=d.attest_count,
        )
        for d in pair_verbs(record)
    )
    return tp.RootObservation(
        record.root_id,
        record.target_class,
        float(record.lemma_count),
        {TRIGGER: observations},
    )


def _cell(
    target: str, tau: str, coding: str, V: int | None = None
) -> tp.CellSpec:
    return tp.CellSpec(
        lang=LANG,
        table=TABLE,
        target_class=target,
        trigger=TRIGGER,
        tau=tau,
        evidence_set="E0",
        coding=coding,
        V=V,
    )


def _subsets() -> dict[str, Callable[[RootRecord], bool]]:
    subsets: dict[str, Callable[[RootRecord], bool]] = {
        "all": lambda r: True,
        "core": lambda r: not r.flags,
    }
    for name in STRATA:
        subsets[f"stratum:{name}"] = (lambda s: lambda r: r.stratum == s)(name)
    return subsets


def spellout_rows(
    records: Sequence[RootRecord], attest: Attest
) -> list[dict[str, Any]]:
    """TP long rows of the spell-out TP (both rules, three codings)."""
    rows: list[dict[str, Any]] = []
    top1000 = frame(records, 1000, "lemma")
    for target, classes in CLASSES.items():
        for subset, member in _subsets().items():
            for tau in TAUS:
                chosen = [
                    r
                    for r in records
                    if r.target_class in classes
                    and member(r)
                    and certain(r, tau, attest)
                ]
                roots = [spellout_root(r, attest) for r in chosen]
                memo = "spell-out -i vs -a/-ui"
                for coding in tp.CODINGS:
                    rows += [
                        finish(row)
                        for row in tp.two_outcome_rows(
                            roots,
                            _cell(target, tau, coding),
                            subset=subset,
                            memo=memo,
                        )
                    ]
                if subset == "all":
                    inside = [
                        spellout_root(r, attest)
                        for r in chosen
                        if r.root_id in top1000
                    ]
                    rows += [
                        finish(row)
                        for row in tp.two_outcome_rows(
                            inside,
                            _cell(target, tau, "dominant", 1000),
                            subset=subset,
                            memo=memo,
                        )
                    ]
    return rows


# ---------------------------------------------------------------------------
# Reproduction of the critique (check_claims.py)
# ---------------------------------------------------------------------------


def _post_palatal(stem: str) -> bool:
    return bool(re.search(r"[șj][td]$", stem))


def check_claims_evidence(record: RootRecord) -> str:
    """``plural_evidence`` of check_claims.py from the legacy columns."""
    seg = record.seg
    legacy = record.legacy
    nde = (
        legacy.exception_reason[4:]
        if legacy.exception_reason.startswith("nde:")
        else ""
    )
    mutated = legacy.mutation is True
    if nde == "gimpe":
        return "lexical_palatal"
    if record.lemma.endswith("ică") and record.plural.endswith("ele"):
        return "suppletion"
    if seg in ("t", "d") and _post_palatal(
        re.sub("[ăeiu]$", "", record.lemma)
    ):
        return "post_palatal"
    if (
        seg in ("c", "g")
        and not mutated
        and re.search(r"ș[cg]", record.lemma[-4:])
        and "t" in legacy.orth_change
    ):
        return "alt"
    if nde in ("ochi", "paduchi"):
        return "nonalt_front"
    if legacy.opportunity == "uri":
        return "uri"
    if legacy.opportunity in ("i", "e"):
        triggers = ("i", "e") if seg in ("c", "g") else ("i",)
        if legacy.opportunity not in triggers:
            return "no_trigger"
        return "alt" if mutated else "nonalt_front"
    return "none"


@dataclass(frozen=True, slots=True)
class CritiqueCounts:
    """Verbalizer pairs on check_claims' underspecified dorsals."""

    frame: str
    pairs: int
    by_suffix: tuple[tuple[str, int], ...]
    nonpalatalizing: int
    theta: float
    roots_in_pair: int
    roots_i_only: int
    roots_nonpal_only: int
    roots_both: int


def critique_counts(records: Sequence[RootRecord]) -> list[CritiqueCounts]:
    """The critique's full-lexicon and top-1000 numbers (pairs).

    Frames as in check_claims.py: nouns with a target ``stem_final``,
    one per lemma, by descending frequency (ties by lemma here; the
    script's pandas sort leaves ties unordered); the top-1000 keeps
    nouns with an -i/-e/-uri opportunity and frequency > 0.
    """
    nouns = sorted(
        (r for r in records if r.seg),
        key=lambda r: (-r.lemma_count, r.root_id),
    )
    full = [r for r in nouns if check_claims_evidence(r) != "none"]
    top = [
        r
        for r in nouns
        if r.legacy.opportunity in ("i", "e", "uri") and r.lemma_count > 0
    ][:1000]
    out: list[CritiqueCounts] = []
    for name, sample in (("full lexicon", full), ("top-1000", top)):
        under = [
            r
            for r in sample
            if r.seg in ("c", "g") and check_claims_evidence(r) == "alt"
        ]
        suffixes = Counter(
            d.suffix for r in under for d in r.derivatives if d.kind == "verb"
        )
        pairs = sum(suffixes.values())
        nonpal = suffixes["-a"] + suffixes["-ui"]
        kinds = Counter(_choice(r) for r in under)
        out.append(
            CritiqueCounts(
                frame=name,
                pairs=pairs,
                by_suffix=tuple(
                    sorted(suffixes.items(), key=lambda kv: (-kv[1], kv[0]))
                ),
                nonpalatalizing=nonpal,
                theta=pairs / math.log(pairs) if pairs > 1 else 0.0,
                roots_in_pair=kinds["i_only"]
                + kinds["nonpal_only"]
                + kinds["both"],
                roots_i_only=kinds["i_only"],
                roots_nonpal_only=kinds["nonpal_only"],
                roots_both=kinds["both"],
            )
        )
    return out


def critique_only_roots(
    records: Sequence[RootRecord], attest: Attest
) -> list[str]:
    """Dorsals check_claims calls alternating (with a pair verb) that
    :func:`certain` does not, as ``lemma/plural``."""
    out: list[str] = []
    for record in records:
        if record.seg not in ("c", "g") or not pair_verbs(record):
            continue
        if check_claims_evidence(record) != "alt":
            continue
        if not certain(record, "dictionary", attest):
            out.append(f"{record.lemma}/{record.plural}")
    return sorted(out)


def _choice(record: RootRecord) -> str:
    chosen = {d.suffix for d in pair_verbs(record)}
    if not chosen:
        return "none"
    if chosen <= set(mo.PALATALIZING_VERBALIZERS):
        return "i_only"
    if not chosen & set(mo.PALATALIZING_VERBALIZERS):
        return "nonpal_only"
    return "both"


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


def _row(rows: Iterable[Mapping[str, Any]], **match: Any) -> Mapping[str, Any]:
    found = [r for r in rows if all(r[k] == v for k, v in match.items())]
    if len(found) != 1:
        raise LookupError(f"{len(found)} rows match {match}")
    return found[0]


def _fmt_theta(value: Any) -> str:
    return "-" if value is None else f"{float(value):.1f}"


def _fmt_bool(value: Any) -> str:
    return {True: "yes", False: "no", None: "n/a"}[value]


def render_report(
    records: Sequence[RootRecord],
    rows: Sequence[Mapping[str, Any]],
    attest: Attest,
) -> str:
    """The text of ``analysis/comparison/romanian_spellout.txt``."""
    lines: list[str] = []
    add = lines.append
    add("Romanian spell-out TP (design R20) — harmonized comparison")
    add("=" * 72)
    add("")
    add("Question: do roots whose plural alternates (underspecified with")
    add("certainty) take the palatalizing verbalizer -i rather than -a/-ui,")
    add("as the spell-out account predicts? Unit: root types. N = certain")
    add("roots with an attested derived verb from the pair -i vs -a/-ui;")
    add("e = roots taking only -a/-ui; roots with both are doublets coded by")
    add("corpus share (dominant >= 90% decides; R19). theta = N/ln N.")
    add(
        "Attestation: plural and verb (infinitive + participle tokens) at tau."
    )
    add("")
    add("1. The critique's numbers, reproduced (check_claims.py definitions)")
    add("-" * 72)
    for item in critique_counts(records):
        suffixes = ", ".join(f"{s} {n}" for s, n in item.by_suffix)
        add(
            f"  {item.frame}: {item.pairs} verbalizer pairs on alternating-plural"
        )
        add(f"    dorsals ({suffixes});")
        add(
            f"    non-palatalizing pairs (-a, -ui) = {item.nonpalatalizing} vs "
            f"theta({item.pairs}) = {item.theta:.1f}"
        )
        add(
            f"    as root types (-i/-a/-ui only): {item.roots_in_pair} roots = "
            f"{item.roots_i_only} -i only + {item.roots_nonpal_only} -a/-ui only"
            f" + {item.roots_both} both"
        )
    add("")
    add("  Differences from the harmonized definition below: pairs vs roots;")
    add("  n includes suffixes outside the competing pair (-ări, -iza, -ăni,")
    add(
        "  -arisi); dorsals only; no attestation threshold; the top-1000 frame"
    )
    add("  is check_claims' (nouns with -i/-e/-uri opportunity), not R24's.")
    only = critique_only_roots(records, attest)
    if only:
        add(
            "  Alternating-plural dorsals with a pair verb that check_claims"
            " counts and the"
        )
        add(
            "  harmonized set does not (plural not marked as mutated in the "
            "lexicon): " + ", ".join(only)
        )
    add("")
    add("2. Harmonized spell-out TP (palatalize = -i default; e as defined)")
    add("-" * 72)
    header = (
        f"  {'class':8s} {'subset':24s} {'tau':10s} {'coding':18s} "
        f"{'N_target':>8s} {'N':>5s} {'e_only':>6s} {'e':>5s} {'theta':>6s} "
        f"{'tol?':>4s} {'frag':>5s}"
    )
    add(header)
    for target in CLASSES:
        for subset in ("all", "core"):
            for tau in TAUS:
                for coding in tp.CODINGS:
                    row = _row(
                        rows,
                        target_class=target,
                        subset=subset,
                        tau=tau,
                        coding=coding,
                        rule="palatalize",
                        V=None,
                    )
                    add(
                        f"  {target:8s} {subset:24s} {tau:10s} {coding:18s} "
                        f"{row['N_target']:>8d} {row['N_trigger']:>5d} "
                        f"{row['e_nonalt']:>6d} {row['e']:>5d} "
                        f"{_fmt_theta(row['theta_trigger']):>6s} "
                        f"{_fmt_bool(row['tolerated_trigger']):>4s} "
                        f"{'-' if row['fragility'] is None else row['fragility']:>5}"
                    )
    add("")
    add("  By base stratum (R22; dominant coding):")
    for target in CLASSES:
        for name in STRATA:
            for tau in ("dictionary", "ge1"):
                row = _row(
                    rows,
                    target_class=target,
                    subset=f"stratum:{name}",
                    tau=tau,
                    coding="dominant",
                    rule="palatalize",
                    V=None,
                )
                add(
                    f"  {target:8s} {name:24s} {tau:10s} N={row['N_trigger']:<5d} "
                    f"e={row['e']:<5d} theta={_fmt_theta(row['theta_trigger'])} "
                    f"tolerated={_fmt_bool(row['tolerated_trigger'])}"
                )
    add("")
    add("  R24 frame, V = 1000 most frequent nouns (dominant coding):")
    for target in CLASSES:
        for tau in ("dictionary", "ge1"):
            row = _row(
                rows,
                target_class=target,
                subset="all",
                tau=tau,
                coding="dominant",
                rule="palatalize",
                V=1000,
            )
            add(
                f"  {target:8s} {tau:10s} N={row['N_trigger']:<5d} e={row['e']:<5d} "
                f"theta={_fmt_theta(row['theta_trigger'])} "
                f"tolerated={_fmt_bool(row['tolerated_trigger'])}"
            )
    add("")
    add("3. Contrast and residue (dictionary tau)")
    add("-" * 72)
    for target, classes in CLASSES.items():
        kinds = Counter(
            _choice(r)
            for r in records
            if r.target_class in classes
            and prespecified(r, "dictionary", attest)
        )
        add(
            f"  {target:8s} non-alternating-plural roots (-uri or resisting): "
            f"-i only {kinds['i_only']}, -a/-ui only {kinds['nonpal_only']}, "
            f"both {kinds['both']}"
        )
    add("")
    unpal = sorted(
        {
            f"{r.lemma}→{d.form}"
            for r in records
            if certain(r, "dictionary", attest)
            for d in pair_verbs(r)
            if d.suffix == "-i" and d.palatal is False
        }
    )
    add(
        f"  -i verbs of alternating-plural roots that keep the plain consonant "
        f"({len(unpal)}; rule-following for spell-out, phonological exceptions"
    )
    add("  for the root):")
    for start in range(0, len(unpal), 6):
        add("    " + ", ".join(unpal[start : start + 6]))
    add("")
    add("Unit and data: noun roots of data/romanian_lexicon_with_freq.csv;")
    add("derived verbs from its derived_verbs/deriv_suffixes columns; counts")
    add("from ron_news_2024_1M (aggregates only here, R9).")
    return "\n".join(lines) + "\n"
