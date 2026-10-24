"""The -log / -fag / -urg table (R27).

One suffixoid, two plural behaviours: person nouns take -i and
palatalize (psiholog/psihologi), thing nouns take -uri and do not
(dialog/dialoguri) or take -e and palatalize (catalog/cataloage).
The same table is built for every language with the same person/thing
classifier, cell outcomes, attested frequencies and derivatives.

Learned suffixoid or lookalike
------------------------------
Many Romanian nouns end in -log or -urg without the Greek suffixoid:
folk words (bârlog 'den', milog 'beggar', zdârlog 'thistle', polog
'swath', ...) and simplexes (amurg 'dusk', burg, foburg). Both
classifiers are the shared ones of :mod:`lexpipe.suffixoids` (the same
for Italian and Czech, R27): a noun counts as the learned suffixoid
(:func:`learned_suffixoid`) when its English gloss carries the
suffixoid's reflex (-log: "log" as in psychologist, dialogue,
catalogue; -fag: "phag"; -urg: "urg" as in surgeon, dramaturge,
metallurgist), when an -olog noun is glossed as a learned profession,
or when a curated decision says so (empty or opaque glosses: dialog,
necrolog 'obituary', liturg 'public official'). The Romanian lexicon
has no etymology templates, so no plugin evidence is passed.
Lookalikes stay in the per-noun table, marked ``learned_suffixoid =
False``, and are left out of the aggregates.

Classifier
----------
:func:`semantic_class` is :func:`lexpipe.suffixoids.semantic_class` on
the first English gloss (rules, in order: a curated decision, animate
non-person patterns, thing patterns, person patterns) and records
which rule decided (``class_basis``). Gender is *not* used: Romanian
masculine/neuter gender is close to the plural allomorph itself (masc.
-i, neut. -uri/-e), so it would make the classifier circular for the
allomorph question. ``unknown`` is kept separate.

Curated decisions
-----------------
``curated/romanian_compare/suffixoid_classes.csv``, keyed by lemma
(the Romanian lexicon has no lexeme_ids) with the guard column
``guard_gloss`` (the first gloss the curator saw); read with
:func:`load_decisions`, checked by :func:`check_decisions`.

Outputs
-------
The tracked ``analysis/comparison/romanian_logtable.csv`` has one row
per noun with lexicon facts only (no per-word Leipzig counts, R9);
``data/romanian_compare/romanian_compare_logtable_full.csv`` adds the
corpus counts and τ flags. τ-dependent aggregates are TP long rows
(:func:`logtable_rows`, table ``logtable``).

A lexicon row whose ``source`` ends in ``| dex_qc`` had no plural in
Wiktionary: ``scripts/dex_qc_main_csv.py`` adds that mark only when it
fills (or replaces) the plural from dexonline (DEX), whose terms nobody
has checked (``dex_qc``, class ``unverified``: not redistributable,
R9). Such nouns (``lexicon_source`` ``…+dex_qc``) are left out of the
tracked table, whose rows must rest on Wiktionary alone, and stay in
the untracked full table; the aggregates count them like every other
root of the harvested lexicon.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any, Final

from lexpipe import curated, suffixoids, tp
from lexpipe.attest import TAUS
from lexpipe.langs.romanian_compare import morphology as mo
from lexpipe.langs.romanian_compare.harmonized import (
    LANG,
    LICENCE_CLASS,
    ObservationCache,
    finish,
)
from lexpipe.langs.romanian_compare.lexicon import fold_key
from lexpipe.langs.romanian_compare.records import RootRecord

SUFFIXOIDS: Final[tuple[str, ...]] = ("log", "fag", "urg")
CLASSES: Final[tuple[str, ...]] = suffixoids.CLASSES
KEY_FIELD: Final[str] = "lemma"
Attest = Callable[[str], Mapping[str, bool]]
Decisions = Mapping[str, suffixoids.Decision]
NO_DECISIONS: Final[Decisions] = MappingProxyType({})
# Row sources by the lexicon's edition (lexicon.lexicon_source labels).
# A "+dex_qc" row's plural was filled from DEX, not confirmed by it, so
# such rows stay out of the tracked table (module docstring).
DEX_QC_LABEL: Final[str] = "dex_qc"
EDITION_SOURCES: Final[Mapping[str, str]] = {
    "en.wiktionary": "wiktionary_en",
    "ro.wiktionary": "wiktionary_ro",
}

LOG_FIELDS: Final[tuple[str, ...]] = (
    "lemma",
    "suffixoid",
    "learned_suffixoid",
    "suffixoid_basis",
    "semantic_class",
    "class_basis",
    "gender",
    "gloss",
    "plural",
    "plural_allomorph",
    "plural_palatalized",
    "abstract_noun",
    "derivatives",
    "lexicon_source",
    "source",
    "licence_class",
)
FULL_EXTRA: Final[tuple[str, ...]] = (
    "lemma_count",
    "plural_count",
    *(f"plural_attested_{tau}" for tau in TAUS),
)


def decisions_path(curated_dir: Path) -> Path:
    """The curated decisions file under ``curated/romanian_compare/``."""
    return curated_dir / suffixoids.DECISIONS_FILE


def load_decisions(
    path: Path,
) -> tuple[dict[str, curated.CuratedRow], dict[str, suffixoids.Decision]]:
    """Curated rows and decisions keyed by lemma (none without a file)."""
    return suffixoids.load_decisions(path, KEY_FIELD)


def check_decisions(
    rows: Mapping[str, curated.CuratedRow], records: Sequence[RootRecord]
) -> None:
    """Stop (:class:`lexpipe.curated.CuratedGuardError`) when a decided
    lemma is gone or its first gloss changed."""
    current: dict[str, dict[str, Any]] = {}
    for record in records:
        if record.lemma in rows and record.lemma not in current:
            current[record.lemma] = {"gloss": first_gloss(record.gloss)}
    curated.check_guards(rows, current)


def suffixoid(lemma: str) -> str:
    """``log`` / ``fag`` / ``urg`` if the noun ends in one, else ""."""
    for item in SUFFIXOIDS:
        if lemma.endswith(item):
            return item
    return ""


def learned_suffixoid(
    lemma: str, gloss: str, decisions: Decisions = NO_DECISIONS
) -> tuple[bool, str]:
    """(is the learned suffixoid, basis) for a noun ending in one
    (:func:`lexpipe.suffixoids.learned_suffixoid`).

    -olog nouns glossed as a learned profession without the reflex
    (africanolog 'africanist', politolog 'political scientist',
    merceolog 'an expert in commodities') count as learned too; the
    folk -olog nouns are glossed as kinds of person (ciotolog 'small
    and ugly man', motolog 'stupid person').
    """
    decision = decisions.get(lemma)
    if decision is not None and decision.learned is not None:
        return decision.learned, "manual"
    ending = suffixoid(lemma)
    if not ending:
        return False, "no suffixoid"
    return suffixoids.learned_suffixoid(
        lemma, ending, gloss, decision=decision
    )


def first_gloss(gloss: str) -> str:
    """The first gloss (before ``;``), trimmed."""
    return suffixoids.first_gloss(gloss)


def semantic_class(
    lemma: str, gloss: str, decisions: Decisions = NO_DECISIONS
) -> tuple[str, str]:
    """(class, basis) for a -log/-fag/-urg noun (module docstring)."""
    return suffixoids.semantic_class(gloss, decision=decisions.get(lemma))


def from_dex(record: RootRecord) -> bool:
    """Did DEX QC fill the row's plural (module docstring)?"""
    return DEX_QC_LABEL in record.lexicon_source.split("+")


def row_source(record: RootRecord) -> str:
    """The Wiktionary edition a lexicon row came from (R9)."""
    edition = record.lexicon_source.split("+")[0]
    try:
        return EDITION_SOURCES[edition]
    except KeyError:
        raise ValueError(
            f"{record.lemma}: unknown lexicon source "
            f"{record.lexicon_source!r}"
        ) from None


def log_records(records: Sequence[RootRecord]) -> list[RootRecord]:
    """g-final noun roots ending in -log, -fag or -urg (NDE included:
    -log nouns are never NDE in this lexicon, but the filter is by
    lemma shape only)."""
    return [r for r in records if r.seg == "g" and suffixoid(r.lemma)]


def _derivative_text(record: RootRecord) -> list[str]:
    out: list[str] = []
    for d in record.derivatives:
        if d.vowel is None:
            outcome = "unread"
        elif d.vowel == "back":
            outcome = "back"
        else:
            outcome = "pal" if d.palatal else "nonpal"
        out.append(f"{d.form}:{outcome}")
    return out


def _abstract(record: RootRecord, nouns: frozenset[str]) -> str:
    candidate = record.lemma + "ie"
    return candidate if candidate in nouns else ""


def table_rows(
    records: Sequence[RootRecord],
    *,
    full: bool,
    decisions: Decisions = NO_DECISIONS,
    attest: Attest | None = None,
) -> list[dict[str, Any]]:
    """One row per -log/-fag/-urg noun, sorted by lemma.

    ``full=True`` adds the corpus columns (untracked file only); without
    it, nouns whose plural DEX filled are left out (module docstring).
    """
    nouns = frozenset(r.lemma for r in records)
    rows: list[dict[str, Any]] = []
    for record in sorted(log_records(records), key=lambda r: r.lemma):
        if not full and from_dex(record):
            continue
        klass, basis = semantic_class(record.lemma, record.gloss, decisions)
        learned, learned_basis = learned_suffixoid(
            record.lemma, record.gloss, decisions
        )
        row: dict[str, Any] = {
            "lemma": record.lemma,
            "suffixoid": "-" + suffixoid(record.lemma),
            "learned_suffixoid": learned,
            "suffixoid_basis": learned_basis,
            "semantic_class": klass,
            "class_basis": basis,
            "gender": record.gender,
            "gloss": first_gloss(record.gloss),
            "plural": record.plural,
            "plural_allomorph": record.plural_allomorph,
            "plural_palatalized": record.plural_palatal,
            "abstract_noun": _abstract(record, nouns),
            "derivatives": _derivative_text(record),
            "lexicon_source": record.lexicon_source,
            "source": row_source(record),
            "licence_class": LICENCE_CLASS,
        }
        if full:
            if attest is None:
                raise ValueError("full rows need an attester")
            row["lemma_count"] = record.lemma_count
            row["plural_count"] = record.plural_count
            flags = attest(fold_key(record.plural)) if record.plural else {}
            for tau in TAUS:
                row[f"plural_attested_{tau}"] = bool(flags.get(tau, False))
        rows.append(row)
    return rows


def logtable_rows(
    records: Sequence[RootRecord],
    cache: ObservationCache,
    decisions: Decisions = NO_DECISIONS,
) -> list[dict[str, Any]]:
    """TP long aggregates by semantic class (table ``logtable``, E0),
    learned suffixoid nouns only."""
    observed = cache.get("E0")
    chosen = [
        r
        for r in log_records(records)
        if learned_suffixoid(r.lemma, r.gloss, decisions)[0]
    ]
    classes = {
        r.root_id: semantic_class(r.lemma, r.gloss, decisions)[0]
        for r in chosen
    }
    rows: list[dict[str, Any]] = []
    for subset in ("all", *CLASSES):
        roots = [
            observed[r.root_id]
            for r in chosen
            if subset == "all" or classes[r.root_id] == subset
        ]
        for trigger in mo.TRIGGERS:
            for tau in TAUS:
                cell = tp.CellSpec(
                    lang=LANG,
                    table="logtable",
                    target_class="g:-log/-fag/-urg",
                    trigger=trigger,
                    tau=tau,
                    evidence_set="E0",
                )
                rows += [
                    finish(row)
                    for row in tp.two_outcome_rows(
                        roots,
                        cell,
                        subset=subset,
                        memo="R27 person/thing classifier",
                    )
                ]
    return rows


@dataclass(frozen=True, slots=True)
class CrossTab:
    """semantic class × plural allomorph × palatalized counts."""

    counts: Mapping[tuple[str, str, str], int]


def crosstab(
    records: Sequence[RootRecord], decisions: Decisions = NO_DECISIONS
) -> CrossTab:
    counts: dict[tuple[str, str, str], int] = {}
    for record in log_records(records):
        if not learned_suffixoid(record.lemma, record.gloss, decisions)[0]:
            continue
        klass = semantic_class(record.lemma, record.gloss, decisions)[0]
        palatal = {True: "pal", False: "nonpal", None: "-"}[
            record.plural_palatal
        ]
        key = (klass, record.plural_allomorph, palatal)
        counts[key] = counts.get(key, 0) + 1
    return CrossTab(dict(sorted(counts.items())))
