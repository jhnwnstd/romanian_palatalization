"""Stage ``logtable``: the -logo/-fago/-urgo table of Italian (R27).

The same learned suffixoid splits by meaning in all three languages:
person nouns take -logi/-fagi (teologi, antropofagi), things and
adjectives -loghi/-faghi (cataloghi, analoghi), with doublets
(sarcofaghi ~ sarcofagi, astrologi ~ astrologhi); the Greek -urgo
person nouns palatalize or vary (taumaturgi; chirurghi ~ chirurgi).
One table per language with the same person/thing classifier, cell
outcomes, attested frequencies and derivatives (R27).

Rows: every in-domain lexeme of ``italian_lexicon_complete.csv`` whose
``suffixoid`` is -logo, -fago or -urgo (the ``features`` stage's
endings with evidence: masculine endings, feminines with evidence,
deverbal nouns left out). ``learned_suffixoid`` and ``suffixoid_basis``
come from the shared detector of :mod:`lexpipe.suffixoids` (Italian
evidence first; a curated ``False`` in
``curated/italian/suffixoid_classes.csv`` marks a lookalike such as
spurgo, listed here but left out of the aggregates);
``semantic_class`` and ``class_basis`` from the shared person/thing
classifier on the first gloss (adjectives: ``adjective``; curated
classes win). Romanian and Czech use the same two classifiers (R27).

``analysis/italian/italian_logtable.csv`` (tracked,
:data:`LOGTABLE_FIELDS`): Wiktionary-derived columns only — identity,
suffixoid, semantic class, gloss, stress, stratum, the dictionary
plurals of the -i and -e cells with their outcomes, the abstract noun
(-logia/-fagia/-urgia, when it is a headword) and the derivatives
linked by the dictionary (``named`` links of the ``derivations``
stage). Sorted by suffixoid, lemma, lexeme_id.

``data/italian/italian_logtable_attestation.csv`` (gitignored:
Leipzig counts, R9; :data:`ATTESTATION_FIELDS`): every plural form of
those lexemes with its corpus counts, τ flags and the cell's doublet
class (R19). The aggregate TP rows by suffixoid × class and τ are in
``italian_tp_long.csv`` (table ``logtable``, built by the ``tp``
stage).
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, Final

from lexpipe import attest, csvio
from lexpipe.context import RunContext, run_stage
from lexpipe.langs.italian import derivations, etymology, features, keys
from lexpipe.langs.italian.layout import Layout

STAGE: Final[str] = "logtable"
# 2: shared R27 classifiers (lexpipe.suffixoids); lookalikes listed.
VERSION: Final[int] = 2
SOURCE: Final[str] = "wiktionary_en"
LICENCE_CLASS: Final[str] = "by-sa"
_ABSTRACT: Final[Mapping[str, str]] = {
    "-logo": "logia",
    "-fago": "fagia",
    "-urgo": "urgia",
}

LOGTABLE_FIELDS: Final[tuple[str, ...]] = (
    "lexeme_id",
    "root_id",
    "lemma",
    "pos",
    "gender",
    "suffixoid",
    "learned_suffixoid",
    "suffixoid_basis",
    "semantic_class",
    "class_basis",
    "gloss",
    "stress_class",
    "stratum",
    "plural",
    "plural_outcome",
    "plural_e",
    "plural_e_outcome",
    "plural_status",
    "abstract_noun",
    "derivatives",
    "source",
    "licence_class",
)
ATTESTATION_FIELDS: Final[tuple[str, ...]] = (
    "lexeme_id",
    "lemma",
    "suffixoid",
    "semantic_class",
    "cell",
    "form",
    "opportunity",
    "mutation",
    "count_cs",
    "count_folded",
    *(f"attested_{tau}" for tau in attest.TAUS),
    "doublet_class",
    "source",
)


def suffixoid_rows(path: Path) -> list[dict[str, str]]:
    """In-domain lexemes ending in a -logo/-fago/-urgo suffixoid
    (learned or lookalike)."""
    return [
        row
        for row in csvio.read_rows(path)
        if row["in_domain"] == "True"
        and row["suffixoid"] in etymology.SUFFIXOIDS
    ]


def named_derivatives(path: Path) -> dict[str, list[str]]:
    """Dictionary-linked derivatives (named tier, accepted) by root."""
    out: dict[str, set[str]] = defaultdict(set)
    for row in csvio.read_rows(path):
        if (
            row["link_status"] == derivations.ACCEPTED
            and row["link_tier"] == derivations.TIER_NAMED
        ):
            out[row["base_root_id"]].add(row["derivative"])
    return {root: sorted(words) for root, words in out.items()}


def abstract_noun(
    lemma: str, suffixoid: str, lemma_keys: frozenset[str]
) -> str:
    """The -logia/-fagia/-urgia headword of a -logo/-fago/-urgo lemma
    (teologo → teologia), or ""."""
    key = keys.match_key(lemma)
    if len(key) < 3 or key[-1] not in "oa":
        return ""
    candidate = key[:-1] + "ia"
    if candidate.endswith(_ABSTRACT[suffixoid]) and candidate in lemma_keys:
        return candidate
    return ""


def table_row(
    row: Mapping[str, str],
    derived: Mapping[str, Sequence[str]],
    lemma_keys: frozenset[str],
) -> dict[str, Any]:
    """One :data:`LOGTABLE_FIELDS` row (Wiktionary-derived only)."""
    return {
        "lexeme_id": row["lexeme_id"],
        "root_id": row["root_id"],
        "lemma": row["word"],
        "pos": row["pos"],
        "gender": row["gender"],
        "suffixoid": row["suffixoid"],
        "learned_suffixoid": csvio.parse_bool(row["learned_suffixoid"]),
        "suffixoid_basis": row["suffixoid_basis"],
        "semantic_class": row["semantic_class"],
        "class_basis": row["class_basis"],
        "gloss": row["first_gloss"],
        "stress_class": row["stress_class"],
        "stratum": row["stratum"],
        "plural": csvio.split_list(row["forms_i"]),
        "plural_outcome": row["outcome_i"],
        "plural_e": csvio.split_list(row["forms_e"]),
        "plural_e_outcome": row["outcome_e"],
        "plural_status": row["plural_status"],
        "abstract_noun": abstract_noun(
            row["word"], row["suffixoid"], lemma_keys
        ),
        "derivatives": list(derived.get(row["root_id"], ())),
        "source": SOURCE,
        "licence_class": LICENCE_CLASS,
    }


def attestation_rows(
    lexemes: Sequence[Mapping[str, str]], cells: Path
) -> list[dict[str, Any]]:
    """Plural forms of the suffixoid lexemes with their counts."""
    by_id = {row["lexeme_id"]: row for row in lexemes}
    out: list[dict[str, Any]] = []
    for row in csvio.read_rows(cells):
        lexeme = by_id.get(row["lexeme_id"])
        if lexeme is None:
            continue
        item: dict[str, Any] = {
            "lexeme_id": row["lexeme_id"],
            "lemma": row["word"],
            "suffixoid": lexeme["suffixoid"],
            "semantic_class": lexeme["semantic_class"],
            "cell": row["cell"],
            "form": row["form"],
            "opportunity": row["opportunity"],
            "mutation": row["mutation"],
            "count_cs": row["count_cs"],
            "count_folded": row["count_folded"],
            "doublet_class": row["doublet_class"],
            "source": "wiktionary_en|leipzig",
        }
        for tau in attest.TAUS:
            item[f"attested_{tau}"] = row[f"attested_{tau}"]
        out.append(item)
    return out


def _body(ctx: RunContext, layout: Layout) -> int:
    lexemes = suffixoid_rows(layout.lexicon_complete)
    derived = named_derivatives(layout.derivations)
    lemma_keys = features.load_lemma_keys(layout.form_index)
    rows = [table_row(row, derived, lemma_keys) for row in lexemes]
    rows.sort(key=lambda r: (r["suffixoid"], r["lemma"], r["lexeme_id"]))
    count = csvio.write_tracked(
        layout.logtable, rows, LOGTABLE_FIELDS, ctx=ctx
    )
    detail = attestation_rows(lexemes, layout.cells)
    csvio.write_rows(
        layout.logtable_attestation, detail, ATTESTATION_FIELDS, ctx=ctx
    )
    print(
        f"  {count:,} -logo/-fago/-urgo lexemes, {len(detail):,} plural forms"
    )
    return 0


def outputs(layout: Layout) -> list[Path]:
    return [layout.logtable, layout.logtable_attestation]


def run(ctx: RunContext) -> int | None:
    """Run the ``logtable`` stage (skipped when its manifest matches)."""
    layout = Layout(ctx)
    return run_stage(
        ctx,
        STAGE,
        inputs=[
            layout.lexicon_complete,
            layout.cells,
            layout.derivations,
            layout.form_index,
        ],
        params={"version": VERSION},
        outputs=outputs(layout),
        body=lambda: _body(ctx, layout),
    )
