"""Stage ``logtable``: the -log/-urg/-graf/-fil table (R27).

One learned suffixoid, two behaviours: person nouns (biolog, chirurg)
take NOM.PL -ové and LOC.PL -lozích/-rzích in the dictionary; thing
nouns (katalog, dialog) are inanimate, NOM.PL -y, LOC.SG -u and LOC.PL
-lozích. -graf and -fil (and -fob) contain no target and serve as the
control: their person nouns take -ové too (fotografové, bibliofilové),
so -ové goes with the suffixoid class, not with the alternation.

Classifiers (the same for every language, R27)
----------------------------------------------
Both classifiers are the shared ones of :mod:`lexpipe.suffixoids`
(Italian and Romanian use them too):

- :func:`learned_suffixoid`: the learned suffixoid or a lookalike
  (profil 'profile' ends like -fil). Czech evidence first: the R26
  suffix evidence of the lexicon (an etymology template or a category;
  basis ``etymology``/``category``; a spelling match is no evidence),
  for a cs.wiktionary lexeme (Czech gloss) the Greek combining vowel
  -o- before the suffixoid (:func:`czech_evidence`), then the English
  gloss's reflex (log, urg, graph, phil), then an
  -olog noun glossed as a learned profession (English, or the Czech
  :data:`CZECH_PROFESSION` for cs.wiktionary glosses). Lookalikes are
  listed (``learned_suffixoid = False``) and left out of the
  aggregates.
- :func:`semantic_class`: person/animate/thing/unknown from the first
  gloss. English glosses go to the shared rules with their
  parenthesised explanation cut (:func:`gloss_head`: paragraf
  'section (part of a legal document …)' is no agent noun in
  -ent), then the Czech layer decides an ``unknown``; Czech glosses
  of cs.wiktionary lexemes go to the Czech layer only
  (:data:`CZECH_PERSON`, :data:`CZECH_THING`; the English -or rule
  read soubor in dekalog's gloss as an agent noun). ``class_basis``
  records which rule decided. Animacy is not used: in Czech it decides
  the paradigm itself (NOM.PL -ové vs -y), so it would make the
  classifier circular.

Curated decisions: ``curated/czech/suffixoid_classes.csv`` (keyed by
lexeme_id; guards ``lemma``, ``readings``, ``first_gloss``), e.g. the
inanimate homographs biograf 'cinema' and telegraf, and the things
dekalog and paragraf.

Outputs
-------
``analysis/czech/czech_logtable.csv``  one row per noun with lexicon
    facts only (forms and allomorphs from Wiktionary; no corpus counts,
    R9), written through the licence gate
``data/czech/czech_logtable_full.csv``  the same rows with corpus
    counts and τ flags of the NOM.PL and LOC.PL forms (untracked)
``analysis/czech/czech_logtable_summary.txt``  class × allomorph
    aggregates at every τ
"""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

from lexpipe import csvio, curated, suffixoids
from lexpipe.attest import TAUS
from lexpipe.context import RunContext, run_stage
from lexpipe.langs.czech import jsonio, ortho
from lexpipe.langs.czech.analysis_layout import AnalysisPaths, analysis_paths
from lexpipe.langs.czech.attestation import CORPUS_TAUS
from lexpipe.langs.czech.layout import CzechPaths, paths_for

STAGE: Final[str] = "logtable"
# 3: the shared classifiers of lexpipe.suffixoids (learned suffixoid
# or lookalike, person/thing) and curated decisions in curated/czech/.
# 4: English rules on the gloss head only, Czech glosses by the Czech
# layer only (audit 2026-10-23).
VERSION: Final[str] = "4"
SUFFIXOIDS: Final[tuple[str, ...]] = ("log", "urg", "graf", "fil")
TARGET_SUFFIXOIDS: Final[frozenset[str]] = frozenset({"log", "urg"})
CLASSES: Final[tuple[str, ...]] = suffixoids.CLASSES
# R26 suffix evidence of the lexicon that counts as evidence of the
# learned suffixoid (a spelling match does not).
EVIDENCE: Final[frozenset[str]] = frozenset({"etymology", "category"})
CS_EDITION: Final[str] = "cs"
COMBINING_VOWEL: Final[str] = "o"
CS_COMPOUND: Final[str] = "combining -o- (cs gloss)"
SOURCE: Final[list[str]] = ["wiktionary_en", "wiktionary_cs"]
LICENCE_CLASS: Final[str] = "by-sa"
CELLS: Final[tuple[str, ...]] = ("NOM.PL", "LOC.SG", "LOC.PL")

# Czech glosses of cs.wiktionary lexemes (decide only after the shared
# classifier returned "unknown").
CZECH_PERSON: Final[re.Pattern[str]] = re.compile(
    r"\b(člověk|muž|žena|odborník|lékař|kdo|ten,? kdo|obdivovatel|"
    r"stoupenec|milovník|vědec|specialista)\b"
)
CZECH_THING: Final[re.Pattern[str]] = re.compile(
    r"\b(přístroj|zařízení|soubor|nádor|písemnost|záznam|text|"
    r"projev|seznam|graf)\b"
)
# -olog nouns of cs.wiktionary glossed as a learned profession.
CZECH_PROFESSION: Final[re.Pattern[str]] = re.compile(
    r"\b(lékař|odborník|vědec|specialista|znalec|ten,? kdo)\b"
)
CZECH_GLOSSES: Final[suffixoids.GlossPatterns] = suffixoids.GlossPatterns(
    label="cs",
    person=CZECH_PERSON,
    thing=CZECH_THING,
    profession=CZECH_PROFESSION,
)
Decisions = Mapping[str, suffixoids.Decision]

FIELDS: Final[tuple[str, ...]] = (
    "lexeme_id",
    "lemma",
    "suffixoid",
    "learned_suffixoid",
    "suffixoid_basis",
    "semantic_class",
    "class_basis",
    "readings",
    "gloss",
    "edition",
    "nom_pl",
    "nom_pl_allomorph",
    "loc_sg",
    "loc_sg_allomorph",
    "loc_pl",
    "loc_pl_allomorph",
    "loc_pl_palatalized",
    "derivatives",
    "source",
    "licence_class",
)
FULL_EXTRA: Final[tuple[str, ...]] = (
    "lemma_count",
    *(f"{cell}_count" for cell in ("nom_pl", "loc_sg", "loc_pl")),
    *(
        f"{cell}_attested_{tau}"
        for cell in ("nom_pl", "loc_sg", "loc_pl")
        for tau in CORPUS_TAUS
    ),
)


def suffixoid(lemma: str) -> str:
    """``log`` / ``urg`` / ``graf`` / ``fil`` if the noun ends in one."""
    for item in SUFFIXOIDS:
        if lemma.endswith(item) and len(lemma) > len(item) + 1:
            return item
    return ""


_PARENTHESES: Final[re.Pattern[str]] = re.compile(r"\s*\([^()]*\)")


def gloss_head(gloss: str) -> str:
    """The gloss without its parenthesised parts: the English person
    rules then match the head noun, not a word of its explanation
    (paragraf 'section (part of a legal document …)': document is no
    agent noun)."""
    return _PARENTHESES.sub("", gloss).strip() or gloss


def semantic_class(
    gloss: str,
    decision: suffixoids.Decision | None = None,
    edition: str = "",
) -> tuple[str, str]:
    """(class, basis): a curated decision, then for an English gloss
    the shared English rules on its head (:func:`gloss_head`) and the
    Czech layer (:func:`lexpipe.suffixoids.semantic_class`); a Czech
    gloss of a cs.wiktionary lexeme gets the Czech layer only (the
    English agent-noun rule would read soubor as -or, audit
    2026-10-23)."""
    if edition != CS_EDITION:
        return suffixoids.semantic_class(
            gloss_head(gloss), decision=decision, extra=CZECH_GLOSSES
        )
    if decision is not None and decision.semantic_class:
        return decision.semantic_class, "manual"
    text = suffixoids.first_gloss(gloss).lower()
    if not text:
        return suffixoids.UNKNOWN_CLASS, "no gloss"
    if CZECH_PERSON.search(text):
        return suffixoids.PERSON_CLASS, f"gloss ({CZECH_GLOSSES.label})"
    if CZECH_THING.search(text):
        return suffixoids.THING_CLASS, f"gloss ({CZECH_GLOSSES.label})"
    return suffixoids.UNKNOWN_CLASS, "gloss"


def czech_evidence(lemma: str, suffix_evidence: str, edition: str) -> str:
    """The Czech evidence that *lemma* ends in the learned suffixoid:
    the lexicon's R26 ``suffix_evidence`` when it is an etymology
    template or a category (a spelling match is none), or, for a
    cs.wiktionary lexeme (whose Czech gloss the English reflex rule
    cannot read), the Greek combining vowel -o- before the suffixoid
    (germano-fil, seismo-graf, bio-log); "" otherwise."""
    if suffix_evidence in EVIDENCE:
        return suffix_evidence
    ending = suffixoid(lemma)
    if (
        edition == CS_EDITION
        and ending
        and lemma[: -len(ending)].endswith(COMBINING_VOWEL)
    ):
        return CS_COMPOUND
    return ""


def learned_suffixoid(
    lemma: str,
    gloss: str,
    evidence: str = "",
    decision: suffixoids.Decision | None = None,
) -> tuple[bool, str]:
    """(learned suffixoid?, basis) of a -log/-urg/-graf/-fil noun
    (:func:`lexpipe.suffixoids.learned_suffixoid`); *evidence* is
    :func:`czech_evidence` ("" when there is none)."""
    return suffixoids.learned_suffixoid(
        lemma,
        suffixoid(lemma),
        gloss,
        decision=decision,
        evidence=evidence,
        extra=CZECH_GLOSSES,
    )


def load_decisions(
    path: Path, nouns: Sequence[Mapping[str, str]]
) -> dict[str, suffixoids.Decision]:
    """Curated decisions, guards checked against the suffixoid nouns
    (a decided lexeme that is gone or changed stops the stage)."""
    rows, decisions = suffixoids.load_decisions(path)
    by_id = {noun["lexeme_id"]: noun for noun in nouns}
    curated.check_guards(
        rows,
        {
            lexeme_id: {
                "lemma": by_id[lexeme_id]["lemma"],
                "readings": csvio.split_list(by_id[lexeme_id]["readings"]),
                "first_gloss": by_id[lexeme_id]["first_gloss"],
            }
            for lexeme_id in rows
            if lexeme_id in by_id
        },
    )
    return decisions


@dataclass
class CellForms:
    """Primary forms of one cell with their allomorph classes."""

    forms: list[str]
    allomorphs: list[str]
    palatalized: bool | None
    count: int
    attested: dict[str, bool]


def _allomorph(row: Mapping[str, str]) -> str:
    """Ending label of a row: the ending itself (-ové, -i, -u, -ích …)."""
    ending = row["ending"]
    return f"-{ending}" if ending else "-0"


def cell_forms(rows: Sequence[Mapping[str, str]], cell: str) -> CellForms:
    chosen = [r for r in rows if r["cell"] == cell and r["primary"] == "True"]
    forms: list[str] = []
    allomorphs: list[str] = []
    attested = dict.fromkeys(CORPUS_TAUS, False)
    palatalized: bool | None = None
    count = 0
    for row in chosen:
        if row["form"] in forms:
            continue
        forms.append(row["form"])
        label = _allomorph(row)
        if label not in allomorphs:
            allomorphs.append(label)
        count += int(row["count"] or 0)
        for tau in CORPUS_TAUS:
            attested[tau] = attested[tau] or row[f"attested_{tau}"] == "True"
        if row["allomorph_class"] == ortho.P:
            palatalized = bool(palatalized) or row["mutation"] == "True"
    return CellForms(forms, allomorphs, palatalized, count, attested)


# ---------------------------------------------------------------------------
# Stage
# ---------------------------------------------------------------------------


def run(ctx: RunContext) -> int | None:
    """Stage entry point (see the module docstring)."""
    paths = paths_for(ctx)
    apaths = analysis_paths(ctx)
    inputs = [paths.lexicon_complete, paths.cells, apaths.derivations]
    if apaths.suffixoid_classes.is_file():
        inputs.append(apaths.suffixoid_classes)
    return run_stage(
        ctx,
        STAGE,
        inputs=inputs,
        params={"version": VERSION},
        outputs=[
            apaths.logtable,
            apaths.logtable_full,
            apaths.logtable_summary,
        ],
        body=lambda: _logtable(ctx, paths, apaths),
    )


def load(
    paths: CzechPaths, apaths: AnalysisPaths
) -> tuple[
    list[dict[str, str]], dict[str, list[dict[str, str]]], dict[str, list[str]]
]:
    """Suffixoid nouns, their cell rows and their derivatives."""
    nouns = [
        row
        for row in jsonio.iter_csv(paths.lexicon_complete)
        if row["pos"] == "noun"
        and row["domain"] == "in"
        and suffixoid(row["lemma"])
    ]
    wanted = {row["lexeme_id"] for row in nouns}
    cells: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in jsonio.iter_csv(paths.cells):
        if row["lexeme_id"] in wanted:
            cells[row["lexeme_id"]].append(row)
    derivatives: dict[str, list[str]] = defaultdict(list)
    for row in jsonio.iter_csv(apaths.derivations):
        if row["base_lexeme_id"] in wanted and row["basis"] in (
            "linked",
            "wikt",
        ):
            text = f"{row['derivative']}:{row['suffix']}:{row['outcome']}"
            if text not in derivatives[row["base_lexeme_id"]]:
                derivatives[row["base_lexeme_id"]].append(text)
    return nouns, cells, derivatives


def table_rows(
    nouns: Sequence[Mapping[str, str]],
    cells: Mapping[str, Sequence[Mapping[str, str]]],
    derivatives: Mapping[str, Sequence[str]],
    *,
    full: bool,
    decisions: Decisions | None = None,
) -> list[dict[str, Any]]:
    """One row per noun, sorted by lemma then lexeme_id."""
    decided = decisions or {}
    out: list[dict[str, Any]] = []
    for noun in sorted(nouns, key=lambda r: (r["lemma"], r["lexeme_id"])):
        decision = decided.get(noun["lexeme_id"])
        klass, basis = semantic_class(
            noun["first_gloss"], decision, noun.get("edition", "")
        )
        learned, learned_basis = learned_suffixoid(
            noun["lemma"],
            noun["first_gloss"],
            czech_evidence(
                noun["lemma"],
                noun.get("suffix_evidence", ""),
                noun.get("edition", ""),
            ),
            decision,
        )
        rows = cells.get(noun["lexeme_id"], ())
        found = {cell: cell_forms(rows, cell) for cell in CELLS}
        row: dict[str, Any] = {
            "lexeme_id": noun["lexeme_id"],
            "lemma": noun["lemma"],
            "suffixoid": f"-{suffixoid(noun['lemma'])}",
            "learned_suffixoid": learned,
            "suffixoid_basis": learned_basis,
            "semantic_class": klass,
            "class_basis": basis,
            "readings": csvio.split_list(noun["readings"]),
            "gloss": suffixoids.first_gloss(noun["first_gloss"]),
            "edition": noun["edition"],
            "nom_pl": found["NOM.PL"].forms,
            "nom_pl_allomorph": found["NOM.PL"].allomorphs,
            "loc_sg": found["LOC.SG"].forms,
            "loc_sg_allomorph": found["LOC.SG"].allomorphs,
            "loc_pl": found["LOC.PL"].forms,
            "loc_pl_allomorph": found["LOC.PL"].allomorphs,
            "loc_pl_palatalized": found["LOC.PL"].palatalized,
            "derivatives": list(derivatives.get(noun["lexeme_id"], ())),
            "source": list(SOURCE),
            "licence_class": LICENCE_CLASS,
        }
        if full:
            row["lemma_count"] = noun["lemma_count"]
            for cell, name in zip(
                CELLS, ("nom_pl", "loc_sg", "loc_pl"), strict=True
            ):
                row[f"{name}_count"] = found[cell].count
                for tau in CORPUS_TAUS:
                    row[f"{name}_attested_{tau}"] = found[cell].attested[tau]
        out.append(row)
    return out


def _logtable(
    ctx: RunContext, paths: CzechPaths, apaths: AnalysisPaths
) -> int:
    nouns, cells, derivatives = load(paths, apaths)
    decisions = load_decisions(apaths.suffixoid_classes, nouns)
    csvio.write_tracked(
        apaths.logtable,
        table_rows(nouns, cells, derivatives, full=False, decisions=decisions),
        FIELDS,
        ctx=ctx,
    )
    full = table_rows(
        nouns, cells, derivatives, full=True, decisions=decisions
    )
    csvio.write_rows(apaths.logtable_full, full, FIELDS + FULL_EXTRA, ctx=ctx)
    text = summary(full)
    ctx.write_tracked_text(apaths.logtable_summary, text, sources=SOURCE)
    print(text, flush=True)
    return 0


def _attested_cell(row: Mapping[str, Any], name: str, tau: str) -> bool:
    if tau == "dictionary":
        return bool(row[name])
    return bool(row[f"{name}_attested_{tau}"])


def summary(rows: Sequence[Mapping[str, Any]]) -> str:
    """Aggregates: suffixoid × class × allomorph at every τ, over the
    learned suffixoid nouns (lookalikes counted apart)."""
    learned = [r for r in rows if r["learned_suffixoid"] is True]
    lookalikes = sorted(
        str(r["lemma"]) for r in rows if r["learned_suffixoid"] is not True
    )
    lines = [
        "Czech -log/-urg/-graf/-fil table (design R27): aggregates",
        "=" * 72,
        "",
        "Classifiers: the shared ones of lexpipe.suffixoids (the same for",
        "Italian and Romanian): learned suffixoid or lookalike (R26",
        "etymology/category evidence, the English gloss's reflex, -olog",
        "professions), person/thing on the first gloss (Czech gloss rules",
        "when the English ones say unknown), curated decisions in",
        "curated/czech/suffixoid_classes.csv. A noun counts at tau when the",
        "cell has a form attested at tau (Leipzig ces_news_2024_1M;",
        "dictionary = any form).",
        "",
        f"Lookalikes left out of the aggregates ({len(lookalikes)}): "
        + (", ".join(lookalikes) or "-"),
        "",
        "Learned suffixoid nouns by suffixoid and class:",
    ]
    counts = Counter((r["suffixoid"], r["semantic_class"]) for r in learned)
    for (item, klass), count in sorted(counts.items()):
        lines.append(f"  {item:<6} {klass:<8} {count:>5}")
    for tau in TAUS:
        lines += ["", f"tau = {tau}: allomorphs of the attested cells"]
        for name in ("nom_pl", "loc_sg", "loc_pl"):
            table: Counter[tuple[str, str, str]] = Counter()
            for row in learned:
                if not _attested_cell(row, name, tau):
                    continue
                label = "|".join(row[f"{name}_allomorph"]) or "-"
                table[(row["suffixoid"], row["semantic_class"], label)] += 1
            if not table:
                continue
            lines.append(f"  {name}")
            for (item, klass, label), count in sorted(table.items()):
                lines.append(
                    f"    {item:<6} {klass:<8} {label:<18} {count:>5}"
                )
    lines += [
        "",
        "LOC.PL palatalized (-lozích, -rzích) among -log/-urg nouns:",
    ]
    for tau in TAUS:
        table2: Counter[tuple[str, str]] = Counter()
        for row in learned:
            if row["suffixoid"] not in ("-log", "-urg"):
                continue
            if not _attested_cell(row, "loc_pl", tau):
                continue
            table2[
                (row["semantic_class"], str(row["loc_pl_palatalized"]))
            ] += 1
        parts = ", ".join(
            f"{k}:{v} {n}" for (k, v), n in sorted(table2.items())
        )
        lines.append(f"  {tau:<10} {parts or '-'}")
    return "\n".join(lines) + "\n"
