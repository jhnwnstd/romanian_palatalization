"""Stage ``lexicon``: Czech noun and adjective lexemes and their cells.

Reads the projected editions written by :mod:`.extract` and writes

``czech_lexicon_raw.csv`` (tracked, Wiktionary data only)
    The lean raw lexicon, :data:`RAW_FIELDS`: what was extracted from
    the lexeme's primary entry and nothing derived — ``lexeme_id``,
    ``lemma``, ``pos``, ``gender`` and ``animacy`` (the readings of the
    head code), ``head_args`` (the head templates as written,
    :func:`lexpipe.kaikki.head_args`; the name of the standard
    cs-noun/cs-adj template left out), ``decl_spec`` (the raw cs-ndecl
    spec), ``ipa`` (raw transcriptions; phrase transcriptions with a
    ``|`` prosodic bar included, the bar backslash-escaped by
    :mod:`lexpipe.csvio`), ``etymology_summary`` (the deciding
    etymology template, :func:`lexpipe.etymology.summary`; the
    etymology *number* is part of ``lexeme_id``), ``first_gloss``, the
    forms of the cells the analyses use (``cell_<CELL>``, as extracted)
    and ``source`` (the primary edition). One row per lexeme, sorted by
    ``lexeme_id``, except proper names and multiword expressions (out
    of domain, R25), which only the full table holds.
    Written through :func:`lexpipe.csvio.write_tracked`.
``czech_lexicon_full.csv`` (gitignored)
    Every lexeme with every column of the stage (:data:`LEXICON_FIELDS`):
    identity (R4), the parsed head code (R13), the raw cs-ndecl spec
    (audit only), the domain decision (R13, R25), IPA, head-line
    derivatives and the forms of the cells the analysis uses;
    ``source`` names the editions the row draws on and
    ``licence_class`` is ``by-sa``. Later stages read this table.
``czech_cells_raw.csv`` (gitignored)
    One row per (lexeme, cell, form, source): every cell form of every
    lexeme except proper names, from ``forms`` entries with exactly one
    case and one number tag (:mod:`.cells`), variants unranked,
    syncretic forms recorded once per cell, with ``form_source`` (R11).
``czech_lexicon_stats.tsv``
    Counts behind the stage report.

Population and identity
-----------------------
English Wiktionary lemma entries with ``pos`` noun, adj or name are
lexemes. A Czech Wiktionary noun or adjective lemma entry is matched to
the English lexemes with the same word and part of speech whose gender
and animacy readings are compatible; its cells are recorded as a
**second source** for them (``primary = False``). A Czech entry whose
word has no English lexeme of that part of speech becomes a lexeme of
its own, unless it is a proper name. Hence ``paradigm_source``:

    en    the English entry has a declension table (primary cells)
    cs    the lexeme exists only in Czech Wiktionary (primary cells)
    none  the English entry has no table (plech); second-source cells
          may still exist and are kept, flagged non-primary. R13 allows
          filling such gaps from MorfFlex only in NC files.

``lexeme_id`` comes from :func:`lexpipe.ids.assign_lexeme_ids` over the
whole population. ``root_id`` (:func:`lexpipe.ids.root_id`) keys on
(word, etymology number): entries of one page and etymology share a
root (noha 'leg' and noha 'leg of furniture'; POS conversions), other
etymologies do not. Czech adds one refinement: when one etymology has
both an animate-only and an inanimate-only masculine noun (šach
'shah' vs šach 'check'), animacy joins the key. These are different
lexical items with different paradigms (šaších vs šachách), and
MorfFlex and DeriNet keep them apart too (R14).

R26 suffix evidence
-------------------
``lemma_suffix`` is the R26 suffix the lemma ends in and
``suffix_evidence`` the first evidence that it is that suffix
(:func:`suffix_evidence`; :mod:`.morphology`): an etymology template or
suffix category, a derivational form-of sense naming a lexeme
(``derivational_formof``: klouček 'diminutive of kluk'), another
lexeme's headword line (Čech f=Češka), a -ský/-cký adjective, or a base
lexeme the word decomposes against when the word is not inherited or
borrowed whole (sedlák ← sedlo; not banka, matka); the -log/-urg/-graf
suffixoids keep their spelling evidence.

Domain (R13, R25)
-----------------
``domain`` is ``out`` with an ``ood_class`` for, in this priority:
``name:<surname|given_name|toponym|other_name>`` (:mod:`.names`;
ethnonyms such as Čech stay in), ``multiword``, ``abbreviation``
(all-capital lemmas), ``hyphenated``, ``dual_remnant``
(cs-decl-noun-dual: ruka, oko; cs duál tables), ``irregular_stem``
(wiktextract class ``irreg-stem``: ucho 'ear', člověk), ``tstem`` (the
kuře type, whose oblique stem adds -et-), ``indeclinable``. Everything
else is ``in``; soft and other non-target stems stay in the domain
and are classified by the features stage.
"""

from __future__ import annotations

import re
import unicodedata
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Final

from lexpipe import csvio
from lexpipe import etymology as core_etymology
from lexpipe import ids, kaikki
from lexpipe.context import RunContext, run_stage
from lexpipe.langs.czech import cells as cellmod
from lexpipe.langs.czech import declspec
from lexpipe.langs.czech import etymology as czech_etymology
from lexpipe.langs.czech import extract, heads, jsonio
from lexpipe.langs.czech import morphology as morph
from lexpipe.langs.czech import names, ortho
from lexpipe.langs.czech.layout import CzechPaths, paths_for
from lexpipe.normalize import normalizer_for

STAGE: Final[str] = "lexicon"
# 2: lean tracked raw lexicon + gitignored full table; IPA phrase
# transcriptions and head links with '|' kept (csvio escapes them).
# 3: derivational form-of lexemes (derivational_formof); R26 suffix
# evidence identified morphologically (:mod:`.morphology`); Czech
# place/name glosses and constellation homographs are names (R25).
VERSION: Final[str] = "3"

EDITION_EN: Final[str] = "en"
EDITION_CS: Final[str] = "cs"
SOURCE_EN: Final[str] = "wiktionary_en"
SOURCE_CS: Final[str] = "wiktionary_cs"
SOURCE_OF: Final[Mapping[str, str]] = {
    EDITION_EN: SOURCE_EN,
    EDITION_CS: SOURCE_CS,
}
LICENCE_CLASS: Final[str] = "by-sa"
PARADIGM_NONE: Final[str] = "none"

LEXEME_POS: Final[tuple[str, ...]] = ("noun", "adj", "name")
CS_POS: Final[tuple[str, ...]] = ("noun", "adj")

DOMAIN_IN: Final[str] = "in"
DOMAIN_OUT: Final[str] = "out"
OOD_MULTIWORD: Final[str] = "multiword"
OOD_ABBREVIATION: Final[str] = "abbreviation"
OOD_HYPHENATED: Final[str] = "hyphenated"
OOD_DUAL: Final[str] = "dual_remnant"
OOD_IRREGULAR: Final[str] = "irregular_stem"
OOD_TSTEM: Final[str] = "tstem"
OOD_INDECLINABLE: Final[str] = "indeclinable"

DUAL_TEMPLATE: Final[str] = "cs-decl-noun-dual"
NDECL_TEMPLATE: Final[str] = "cs-ndecl"
MANUAL_TEMPLATES: Final[frozenset[str]] = frozenset(
    {DUAL_TEMPLATE, "cs-adecl-manual"}
)
FORM_EXPLICIT: Final[str] = "explicit"
FORM_DEFAULT: Final[str] = "module_default"

# Cells written as wide columns into the tracked lexicon: the cells the
# P/N classification of design §6 uses, plus GEN.SG (the stem cell).
WIDE_CELLS: Final[tuple[str, ...]] = (
    "GEN.SG",
    "DAT.SG",
    "VOC.SG",
    "LOC.SG",
    "INS.SG",
    "NOM.PL",
    "LOC.PL",
    cellmod.HARD_ADJECTIVE_CELL,
)
HEAD_LINKS: Final[tuple[tuple[str, str], ...]] = (
    ("f", "head_feminine"),
    ("m", "head_masculine"),
    ("dim", "head_diminutive"),
    ("adj", "head_adjective"),
)

LEXICON_FIELDS: Final[tuple[str, ...]] = (
    "lexeme_id",
    "root_id",
    "lemma",
    "pos",
    "etymology",
    "edition",
    "paradigm_source",
    "domain",
    "ood_class",
    "gender",
    "animacy",
    "readings",
    "gender_ambiguous",
    "animacy_ambiguous",
    "pluralia_tantum",
    "by_sense",
    "head_template",
    "head_code",
    "decl_templates",
    "decl_spec",
    "spec_flags",
    "explicit_cells",
    "stem_class",
    "lemma_suffix",
    "suffix_evidence",
    "derivational_formof",
    "table_animacy",
    "capitalized",
    "name_evidence",
    "first_gloss",
    "ipa",
    "ipa_kind",
    *(column for _, column in HEAD_LINKS),
    "n_tables",
    "n_cells",
    "n_cells_second",
    "second_source_entries",
    *(f"cell_{cell}" for cell in WIDE_CELLS),
    "source",
    "licence_class",
)
# The tracked raw lexicon: extracted columns only (module docstring).
RAW_FIELDS: Final[tuple[str, ...]] = (
    "lexeme_id",
    "lemma",
    "pos",
    "gender",
    "animacy",
    "head_args",
    "decl_spec",
    "ipa",
    "etymology_summary",
    "first_gloss",
    *(f"cell_{cell}" for cell in WIDE_CELLS),
    "source",
)
# The standard head template of each part of speech (head_args leaves
# its name out).
DEFAULT_HEADS: Final[Mapping[str, str]] = {"noun": "cs-noun", "adj": "cs-adj"}
CELL_FIELDS: Final[tuple[str, ...]] = (
    "lexeme_id",
    "cell",
    "form",
    "source",
    "primary",
    "table",
    "table_animacy",
    "order",
    "form_source",
    "qualifiers",
    "syncretic_cells",
    "in_other_edition",
    "tag_error",
)

_ABBREVIATION_TAGS: Final[frozenset[str]] = frozenset(
    {"abbreviation", "acronym", "initialism"}
)
_CS_INDECLINABLE: Final[str] = "nesklonné"
_CS_DUAL: Final[re.Pattern[str]] = re.compile(r"^duálov")


# ---------------------------------------------------------------------------
# Lexeme records
# ---------------------------------------------------------------------------


@dataclass
class Lexeme:
    """One lexeme while the stage builds it."""

    record: Mapping[str, Any]
    edition: str
    lemma: str
    pos: str
    etymology: str
    head: heads.HeadInfo
    head_template: str = ""
    head_args: Mapping[str, str] = field(default_factory=dict)
    templates: tuple[str, ...] = ()
    specs: tuple[declspec.DeclSpec | None, ...] = ()
    paradigm: cellmod.Paradigm = field(default_factory=cellmod.Paradigm)
    second: list[cellmod.CellForm] = field(default_factory=list)
    second_entries: int = 0
    name: names.NameVerdict = names.COMMON
    domain: str = DOMAIN_IN
    ood_class: str = ""
    stem_class: str = ""
    cs_indeclinable: bool = False
    cs_dual: bool = False
    lexeme_id: str = ""
    root_id: str = ""

    @property
    def is_name(self) -> bool:
        return self.name.is_name

    @property
    def has_paradigm(self) -> bool:
        return bool(self.paradigm.cells)

    @property
    def paradigm_source(self) -> str:
        return self.edition if self.has_paradigm else PARADIGM_NONE


def _nfc(text: str) -> str:
    return unicodedata.normalize("NFC", text).strip()


def _templates(record: Mapping[str, Any]) -> list[tuple[str, dict[str, str]]]:
    out: list[tuple[str, dict[str, str]]] = []
    for template in record.get("inflection_templates") or ():
        if not isinstance(template, Mapping):
            continue
        args = template.get("args")
        out.append(
            (
                str(template.get("name") or ""),
                (
                    {str(k): str(v) for k, v in args.items()}
                    if isinstance(args, Mapping)
                    else {}
                ),
            )
        )
    return out


def en_lexeme(record: Mapping[str, Any]) -> Lexeme:
    """An English Wiktionary lexeme from its projected record."""
    pos = str(record.get("pos") or "")
    found = heads.head_template(record)
    templates = _templates(record)
    specs = tuple(
        (
            declspec.parse_decl_spec(args.get("1", ""))
            if name == NDECL_TEMPLATE
            else None
        )
        for name, args in templates
    )
    paradigm = cellmod.en_cells(record, adjective=pos == "adj")
    stem_class = next((stem for _, stem in paradigm.tables if stem), "")
    return Lexeme(
        record=record,
        edition=EDITION_EN,
        lemma=_nfc(str(record.get("word") or "")),
        pos=pos,
        etymology=kaikki.etymology_label(record),
        head=heads.head_info(record),
        head_template=found[0] if found else "",
        head_args=found[1] if found else {},
        templates=tuple(name for name, _ in templates),
        specs=specs,
        paradigm=paradigm,
        name=names.classify_en(record),
        stem_class=stem_class,
    )


def cs_lexeme(record: Mapping[str, Any]) -> Lexeme:
    """A Czech Wiktionary lexeme (or match candidate) from its record."""
    pos = str(record.get("pos") or "")
    tags = [t for t in record.get("tags") or () if isinstance(t, str)]
    raw_tags = [t for t in record.get("raw_tags") or () if isinstance(t, str)]
    readings = heads.readings_from_tags(tags)
    animacies = {r.animacy for r in readings if r.gender == "m"}
    animacy = animacies.pop() if len(animacies) == 1 else ""
    return Lexeme(
        record=record,
        edition=EDITION_CS,
        lemma=_nfc(str(record.get("word") or "")),
        pos=pos,
        etymology=kaikki.etymology_label(record),
        head=heads.HeadInfo(code="", readings=readings),
        paradigm=cellmod.cs_cells(
            record, adjective=pos == "adj", animacy=animacy
        ),
        name=names.classify_cs(record),
        cs_indeclinable=_CS_INDECLINABLE in raw_tags,
        cs_dual=any(_CS_DUAL.match(t) for t in raw_tags),
    )


# ---------------------------------------------------------------------------
# Domain
# ---------------------------------------------------------------------------


def _is_abbreviation(lexeme: Lexeme) -> bool:
    """All-capital lemmas (UFO, AIDS, B) or entries whose every sense is
    tagged abbreviation, acronym or initialism."""
    letters = [ch for ch in lexeme.lemma if ch.isalpha()]
    if letters and all(ch.isupper() for ch in letters):
        if len(letters) > 1 or len(lexeme.lemma) == 1:
            return True
    senses = [
        s for s in lexeme.record.get("senses") or () if isinstance(s, Mapping)
    ]
    return bool(senses) and all(
        _ABBREVIATION_TAGS.intersection(s.get("tags") or ()) for s in senses
    )


def _is_indeclinable(lexeme: Lexeme) -> bool:
    if lexeme.cs_indeclinable:
        return True
    if lexeme.head_args.get("indecl", "").strip() not in ("", "0"):
        return True
    specs = [spec for spec in lexeme.specs if spec is not None]
    return bool(specs) and all(
        spec.has_flag("indecl") or spec.has_flag("mostlyindecl")
        for spec in specs
    )


def domain_of(lexeme: Lexeme) -> tuple[str, str]:
    """``(domain, ood_class)`` (see the module docstring)."""
    if lexeme.is_name:
        return DOMAIN_OUT, f"name:{lexeme.name.name_class}"
    if " " in lexeme.lemma or any(
        spec is not None and spec.multiword for spec in lexeme.specs
    ):
        return DOMAIN_OUT, OOD_MULTIWORD
    if _is_abbreviation(lexeme):
        return DOMAIN_OUT, OOD_ABBREVIATION
    if "-" in lexeme.lemma.strip("-"):
        return DOMAIN_OUT, OOD_HYPHENATED
    if (
        DUAL_TEMPLATE in lexeme.templates
        or lexeme.cs_dual
        or lexeme.paradigm.dropped.get(cellmod.DROP_DUAL, 0)
    ):
        return DOMAIN_OUT, OOD_DUAL
    if lexeme.stem_class == "irreg-stem":
        return DOMAIN_OUT, OOD_IRREGULAR
    if lexeme.stem_class == "t-stem" or any(
        spec is not None and spec.has_flag("tstem") for spec in lexeme.specs
    ):
        return DOMAIN_OUT, OOD_TSTEM
    if _is_indeclinable(lexeme):
        return DOMAIN_OUT, OOD_INDECLINABLE
    return DOMAIN_IN, ""


# ---------------------------------------------------------------------------
# Matching the Czech edition
# ---------------------------------------------------------------------------


def readings_compatible(
    left: Sequence[heads.GenderReading], right: Sequence[heads.GenderReading]
) -> bool:
    """Do two reading lists share a gender (and, for masculines, an
    animacy, an unknown animacy matching either)?"""
    if not left or not right:
        return not left and not right
    for a in left:
        for b in right:
            if a.gender != b.gender:
                continue
            if a.gender != "m" or not a.animacy or not b.animacy:
                return True
            if a.animacy == b.animacy:
                return True
    return False


@dataclass(frozen=True, slots=True)
class MatchStats:
    """How the Czech edition's entries were used."""

    matched: int
    cs_only: int
    names: int
    unmatched: int


def match_cs(
    en: Sequence[Lexeme], cs: Sequence[Lexeme]
) -> tuple[list[Lexeme], MatchStats]:
    """Attach Czech cells to English lexemes; return the Czech-only ones.

    Czech entries are visited in record order (the projection is
    sorted), so the order of second-source cells is deterministic.
    """
    by_key: dict[tuple[str, str], list[Lexeme]] = defaultdict(list)
    for lexeme in en:
        if not lexeme.is_name:
            by_key[(lexeme.lemma, lexeme.pos)].append(lexeme)
    cs_only: list[Lexeme] = []
    matched = names_count = unmatched = 0
    for candidate in cs:
        if candidate.is_name:
            names_count += 1
            continue
        targets = by_key.get((candidate.lemma, candidate.pos))
        if not targets:
            cs_only.append(candidate)
            continue
        compatible = [
            t
            for t in targets
            if readings_compatible(t.head.readings, candidate.head.readings)
        ]
        if not compatible:
            unmatched += 1
            continue
        matched += 1
        for target in compatible:
            offset = max((c.order for c in target.second), default=-1) + 1
            target.second.extend(
                cellmod.CellForm(
                    cell=c.cell,
                    form=c.form,
                    table=target.second_entries,
                    table_animacy=c.table_animacy,
                    order=offset + c.order,
                    qualifiers=c.qualifiers,
                )
                for c in candidate.paradigm.cells
            )
            target.second_entries += 1
    return cs_only, MatchStats(matched, len(cs_only), names_count, unmatched)


# ---------------------------------------------------------------------------
# Identifiers
# ---------------------------------------------------------------------------


def _single_animacy(lexeme: Lexeme) -> str:
    """'an' or 'in' if the noun's masculine readings have one animacy."""
    if lexeme.pos != "noun" or lexeme.head.genders != ("m",):
        return ""
    animacies = lexeme.head.animacies
    return animacies[0] if len(animacies) == 1 else ""


def _root_key(parts: Sequence[str]) -> ids.RootKeyFn:
    frozen = tuple(parts)
    return lambda _entry: frozen


def assign_identifiers(lexemes: Sequence[Lexeme]) -> None:
    """Set ``lexeme_id`` and ``root_id`` on every lexeme (R4)."""
    lexeme_ids = ids.assign_lexeme_ids("cs", [lx.record for lx in lexemes])
    animacies: dict[tuple[str, str], set[str]] = defaultdict(set)
    for lexeme in lexemes:
        animacy = _single_animacy(lexeme)
        if animacy:
            animacies[(lexeme.lemma, lexeme.etymology)].add(animacy)
    split = {key for key, found in animacies.items() if len(found) > 1}
    for lexeme, lexeme_id in zip(lexemes, lexeme_ids, strict=True):
        lexeme.lexeme_id = lexeme_id
        key = (lexeme.lemma, lexeme.etymology)
        animacy = _single_animacy(lexeme) if key in split else ""
        lexeme.root_id = ids.root_id(
            "cs",
            lexeme.record,
            _root_key([lexeme.lemma, lexeme.etymology, animacy]),
        )


# ---------------------------------------------------------------------------
# Rows
# ---------------------------------------------------------------------------


def _table_spec(lexeme: Lexeme, table: int) -> declspec.DeclSpec | None:
    """The spec of table *table* when tables and templates align."""
    if len(lexeme.paradigm.tables) == len(lexeme.templates) and table < len(
        lexeme.specs
    ):
        return lexeme.specs[table]
    return None


def _union_spec(lexeme: Lexeme) -> declspec.DeclSpec:
    specs = [spec for spec in lexeme.specs if spec is not None]
    overrides: set[tuple[str, str]] = set()
    flags: set[str] = set()
    for spec in specs:
        overrides |= spec.overrides
        flags |= set(spec.flags)
    return declspec.DeclSpec(
        raw="",
        overrides=frozenset(overrides),
        flags=tuple(sorted(flags)),
        manual=any(spec.manual for spec in specs),
    )


def _table_template(lexeme: Lexeme, table: int) -> str:
    if len(lexeme.paradigm.tables) == len(lexeme.templates):
        return lexeme.templates[table] if table < len(lexeme.templates) else ""
    return lexeme.templates[0] if len(set(lexeme.templates)) == 1 else ""


def form_source(lexeme: Lexeme, cell: cellmod.CellForm) -> str:
    """``explicit`` or ``module_default`` for a primary cell form (R11)."""
    if lexeme.edition == EDITION_CS:
        return FORM_EXPLICIT  # cs.wiktionary tables are literal
    template = _table_template(lexeme, cell.table)
    if template in MANUAL_TEMPLATES:
        return FORM_EXPLICIT
    if lexeme.pos == "adj" and template != NDECL_TEMPLATE:
        return FORM_DEFAULT
    spec = _table_spec(lexeme, cell.table) or _union_spec(lexeme)
    if spec.manual:
        return FORM_EXPLICIT
    animate = cell.table_animacy == "an" or (
        not cell.table_animacy and lexeme.head.labels == ("m.an",)
    )
    case, number = cell.cell.split(".")[:2]
    explicit = declspec.explicit_cells(spec, masculine_animate=animate)
    return FORM_EXPLICIT if (case, number) in explicit else FORM_DEFAULT


def _forms_by_cell(cells: Iterable[cellmod.CellForm]) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for cell in cells:
        forms = out.setdefault(cell.cell, [])
        if cell.form not in forms:
            forms.append(cell.form)
    return out


def cell_rows(lexeme: Lexeme) -> list[dict[str, Any]]:
    """Rows of ``czech_cells_raw.csv`` for one lexeme."""
    primary = list(lexeme.paradigm.cells)
    second = list(lexeme.second)
    primary_forms = _forms_by_cell(primary)
    second_forms = _forms_by_cell(second)
    rows: list[dict[str, Any]] = []
    for group, is_primary, other in (
        (primary, True, second_forms),
        (second, False, primary_forms),
    ):
        source_edition = lexeme.edition if is_primary else EDITION_CS
        sync = cellmod.syncretism(group)
        for cell in group:
            others = other.get(cell.cell)
            rows.append(
                {
                    "lexeme_id": lexeme.lexeme_id,
                    "cell": cell.cell,
                    "form": cell.form,
                    "source": SOURCE_OF[source_edition],
                    "primary": is_primary,
                    "table": cell.table,
                    "table_animacy": cell.table_animacy,
                    "order": cell.order,
                    "form_source": (
                        form_source(lexeme, cell)
                        if is_primary
                        else FORM_EXPLICIT
                    ),
                    "qualifiers": [q for q in cell.qualifiers if "|" not in q],
                    "syncretic_cells": [
                        c
                        for c in sync[(cell.table, cell.form)]
                        if c != cell.cell
                    ],
                    "in_other_edition": (
                        None if others is None else cell.form in others
                    ),
                    "tag_error": cell.tag_error,
                }
            )
    rows.sort(key=_cell_row_key)
    return rows


def _cell_row_key(row: Mapping[str, Any]) -> tuple[Any, ...]:
    return (
        cellmod.cell_sort_key(str(row["cell"])),
        not row["primary"],
        int(row["table"]),
        int(row["order"]),
        str(row["form"]),
    )


def _ipa(record: Mapping[str, Any]) -> tuple[list[str], str]:
    """Raw IPA strings and the kind of the first parsable one.

    Phrase transcriptions holding ``|`` (the prosodic bar of multiword
    terms) are kept: :mod:`lexpipe.csvio` escapes a ``|`` inside a list
    item.
    """
    raw = list(kaikki.ipa_list(record))
    normalizer = normalizer_for("cs")
    kinds = [
        parsed.kind for r in raw if (parsed := normalizer.ipa(r)) is not None
    ]
    return raw, (str(kinds[0]) if kinds else "")


def derivational_formof(record: Mapping[str, Any]) -> list[str]:
    """``relation:target`` of the entry's derivational form-of senses
    (:func:`.extract.derivational_relation`: ``diminutive:kluk``,
    ``female_equivalent:Čech``), in source order."""
    out: list[str] = []
    for sense in record.get("senses") or ():
        if not isinstance(sense, Mapping):
            continue
        relation = extract.derivational_relation(sense)
        if not relation:
            continue
        for target in sense.get("form_of") or sense.get("alt_of") or ():
            word = target.get("word") if isinstance(target, Mapping) else ""
            if isinstance(word, str) and word.strip():
                item = f"{relation}:{_nfc(word)}"
                if item not in out:
                    out.append(item)
    return out


def wiktionary_stratum(record: Mapping[str, Any]) -> str:
    """The entry's R22 Wiktionary stratum (:mod:`.etymology`)."""
    found = czech_etymology.steps(record.get("etymology_templates") or ())
    text = str(record.get("etymology_text") or "")
    word = str(record.get("word") or "")
    return czech_etymology.wiktionary_stratum(found, text, word=word)[0]


def suffix_evidence(
    lexeme: Lexeme,
    suffix: str,
    formof: Sequence[str],
    context: morph.SuffixContext,
) -> str:
    """R26 evidence kind for the lexeme's suffix (:mod:`.morphology`):
    etymology or category first, then the morphological kinds, then
    the orthographic evidence of the -log/-urg/-graf suffixoids."""
    evidence = ortho.suffix_evidence(
        suffix,
        [
            t
            for t in lexeme.record.get("etymology_templates") or ()
            if isinstance(t, Mapping)
        ],
        names.czech_categories(lexeme.record),
    )
    if evidence and evidence != morph.ORTHOGRAPHIC:
        return evidence
    found = morph.morphological_evidence(
        lexeme.lemma,
        lexeme.pos,
        suffix,
        formof_targets=[item.split(":", 1)[1] for item in formof],
        stratum=(
            wiktionary_stratum(lexeme.record)
            if lexeme.edition == EDITION_EN
            else czech_etymology.UNKNOWN
        ),
        context=context,
    )
    return found or evidence


def lexicon_row(
    lexeme: Lexeme, context: morph.SuffixContext = morph.EMPTY_CONTEXT
) -> dict[str, Any]:
    """The full lexicon row of one lexeme (*context*: the lexicon-wide
    facts the R26 suffix evidence needs, :func:`suffix_context`)."""
    head = lexeme.head
    forms = _forms_by_cell(lexeme.paradigm.cells)
    ipa, ipa_kind = _ipa(lexeme.record)
    specs = [spec for spec in lexeme.specs if spec is not None]
    explicit = sorted(
        {
            f"{case}.{number}"
            for spec in specs
            for case, number in spec.overrides
        },
        key=cellmod.cell_sort_key,
    )
    sources = [SOURCE_OF[lexeme.edition]]
    if lexeme.second:
        sources.append(SOURCE_CS)
    suffix = ortho.lemma_suffix(lexeme.lemma)
    formof = derivational_formof(lexeme.record)
    evidence = suffix_evidence(lexeme, suffix, formof, context)
    row: dict[str, Any] = {
        "lexeme_id": lexeme.lexeme_id,
        "root_id": lexeme.root_id,
        "lemma": lexeme.lemma,
        "pos": lexeme.pos,
        "etymology": lexeme.etymology,
        "edition": lexeme.edition,
        "paradigm_source": lexeme.paradigm_source,
        "domain": lexeme.domain,
        "ood_class": lexeme.ood_class,
        "gender": list(head.genders),
        "animacy": list(head.animacies),
        "readings": list(head.labels),
        "gender_ambiguous": head.gender_ambiguous,
        "animacy_ambiguous": head.animacy_ambiguous,
        "pluralia_tantum": head.pluralia_tantum,
        "by_sense": head.by_sense,
        "head_template": lexeme.head_template,
        "head_code": head.code,
        "decl_templates": list(lexeme.templates),
        "decl_spec": [spec.raw for spec in specs],
        "spec_flags": sorted({flag for spec in specs for flag in spec.flags}),
        "explicit_cells": explicit,
        "stem_class": lexeme.stem_class,
        "lemma_suffix": suffix,
        "suffix_evidence": evidence,
        "derivational_formof": formof,
        "table_animacy": [animacy for animacy, _ in lexeme.paradigm.tables],
        "capitalized": lexeme.lemma[:1].isupper(),
        "name_evidence": lexeme.name.evidence,
        "first_gloss": names.first_gloss(lexeme.record),
        "ipa": ipa,
        "ipa_kind": ipa_kind,
        "n_tables": len(lexeme.paradigm.tables),
        "n_cells": len(forms),
        "n_cells_second": len(_forms_by_cell(lexeme.second)),
        "second_source_entries": lexeme.second_entries,
        "source": sources,
        "licence_class": LICENCE_CLASS,
    }
    for arg, column in HEAD_LINKS:
        row[column] = _split_links(lexeme.head_args.get(arg, ""))
    for cell in WIDE_CELLS:
        row[f"cell_{cell}"] = forms.get(cell, [])
    return row


def _split_links(value: str) -> list[str]:
    """Head-template link lists ``žačka,žákyně`` (inline modifiers cut;
    a ``|`` in an item is escaped by :mod:`lexpipe.csvio`)."""
    out: list[str] = []
    for item in value.split(","):
        word = re.sub(r"<[^<>]*>", "", item).strip()
        if word and word not in out:
            out.append(word)
    return out


def tracked(lexeme: Lexeme) -> bool:
    """Does the tracked raw lexicon hold *lexeme*? Not proper names or
    multiword expressions (R25; the full table holds them)."""
    return not lexeme.is_name and lexeme.ood_class != OOD_MULTIWORD


def raw_row(lexeme: Lexeme, full: Mapping[str, Any]) -> dict[str, Any]:
    """The tracked raw lexicon row of one lexeme (:data:`RAW_FIELDS`),
    from its full row *full*."""
    found = core_etymology.steps(
        lexeme.record.get("etymology_templates") or (), "cs"
    )
    row: dict[str, Any] = {
        "lexeme_id": lexeme.lexeme_id,
        "lemma": lexeme.lemma,
        "pos": lexeme.pos,
        "gender": full["gender"],
        "animacy": full["animacy"],
        "head_args": kaikki.head_args(
            lexeme.record, DEFAULT_HEADS.get(lexeme.pos, "")
        ),
        "decl_spec": full["decl_spec"],
        "ipa": full["ipa"],
        "etymology_summary": core_etymology.summary(found, "cs"),
        "first_gloss": full["first_gloss"],
        "source": SOURCE_OF[lexeme.edition],
    }
    for cell in WIDE_CELLS:
        row[f"cell_{cell}"] = full[f"cell_{cell}"]
    return row


# ---------------------------------------------------------------------------
# Stage
# ---------------------------------------------------------------------------


def run(ctx: RunContext) -> int | None:
    """Stage entry point (see the module docstring)."""
    paths = paths_for(ctx)
    return run_stage(
        ctx,
        STAGE,
        inputs=[paths.en_entries, paths.cs_entries],
        params={"version": VERSION},
        outputs=[
            paths.lexicon_full,
            paths.lexicon_raw,
            paths.cells_raw,
            paths.lexicon_stats,
        ],
        body=lambda: _build(ctx, paths),
    )


def load_lexemes(paths: CzechPaths) -> tuple[list[Lexeme], MatchStats]:
    """All lexemes, identified, with domain decisions (no I/O writes)."""
    en = [
        en_lexeme(record)
        for record in jsonio.read_jsonl(paths.en_entries)
        if record.get("pos") in LEXEME_POS
    ]
    cs = [
        cs_lexeme(record)
        for record in jsonio.read_jsonl(paths.cs_entries)
        if record.get("pos") in CS_POS
    ]
    cs_only, stats = match_cs(en, cs)
    lexemes = en + cs_only
    for lexeme in lexemes:
        lexeme.domain, lexeme.ood_class = domain_of(lexeme)
    assign_identifiers(lexemes)
    lexemes.sort(key=lambda lx: lx.lexeme_id)
    return lexemes, stats


# Headword-line links that name a derivative of the entry (R26
# ``headword`` evidence): feminine, diminutive, relational adjective.
# m= names the masculine base of a feminine, not a derivative.
DERIVATIVE_HEAD_ARGS: Final[tuple[str, ...]] = ("f", "dim", "adj")


def verb_lemmas(paths: CzechPaths) -> list[str]:
    """Verb lemmas of both editions (bases of -ák, -ník … nouns)."""
    out: set[str] = set()
    for path in (paths.en_entries, paths.cs_entries):
        for record in jsonio.read_jsonl(path):
            if record.get("pos") == "verb":
                word = _nfc(str(record.get("word") or ""))
                if word and " " not in word:
                    out.add(word)
    return sorted(out)


def suffix_context(
    lexemes: Sequence[Lexeme], verbs: Iterable[str] = ()
) -> morph.SuffixContext:
    """The lexicon-wide facts of the R26 evidence (:mod:`.morphology`):
    every lexeme lemma, the derivatives named on headword lines and
    the base index over nouns, adjectives and *verbs*."""
    index = morph.BaseIndex()
    named: set[str] = set()
    for lexeme in lexemes:
        if " " not in lexeme.lemma:
            index.add(lexeme.lemma, lexeme.pos)
        for arg in DERIVATIVE_HEAD_ARGS:
            named.update(_split_links(lexeme.head_args.get(arg, "")))
    for verb in verbs:
        index.add(verb, "verb")
    return morph.SuffixContext(
        lemmas=frozenset(lx.lemma for lx in lexemes),
        head_named=frozenset(named),
        index=index,
    )


def _build(ctx: RunContext, paths: CzechPaths) -> int:
    lexemes, stats = load_lexemes(paths)
    context = suffix_context(lexemes, verb_lemmas(paths))
    full = [lexicon_row(lexeme, context) for lexeme in lexemes]
    csvio.write_rows(paths.lexicon_full, full, LEXICON_FIELDS, ctx=ctx)
    csvio.write_tracked(
        paths.lexicon_raw,
        (
            raw_row(lexeme, row)
            for lexeme, row in zip(lexemes, full, strict=True)
            if tracked(lexeme)
        ),
        RAW_FIELDS,
        ctx=ctx,
    )
    counts: Counter[tuple[str, str]] = Counter()

    def all_cell_rows() -> Iterable[dict[str, Any]]:
        for lexeme in lexemes:
            if lexeme.is_name:
                continue  # R25: proper names are out; their cells are not kept
            for row in cell_rows(lexeme):
                if row["primary"]:
                    counts[("form_source", str(row["form_source"]))] += 1
                counts[("cell_rows", str(row["source"]))] += 1
                yield row

    csvio.write_rows(paths.cells_raw, all_cell_rows(), CELL_FIELDS, ctx=ctx)
    _write_stats(ctx, paths, lexemes, stats, counts)
    return 0


def _write_stats(
    ctx: RunContext,
    paths: CzechPaths,
    lexemes: Sequence[Lexeme],
    stats: MatchStats,
    counts: Counter[tuple[str, str]],
) -> None:
    table: Counter[tuple[str, str]] = Counter(counts)
    for lexeme in lexemes:
        table[("lexemes", f"{lexeme.edition}:{lexeme.pos}")] += 1
        table[
            ("domain", f"{lexeme.pos}:{lexeme.domain}:{lexeme.ood_class}")
        ] += 1
        table[
            ("paradigm_source", f"{lexeme.pos}:{lexeme.paradigm_source}")
        ] += 1
        if lexeme.pos == "noun" and lexeme.domain == DOMAIN_IN:
            label = "|".join(lexeme.head.labels) or "?"
            table[("noun_in_domain_readings", label)] += 1
            if lexeme.head.animacy_ambiguous:
                table[("noun_in_domain_flags", "animacy_ambiguous")] += 1
            if lexeme.head.pluralia_tantum:
                table[("noun_in_domain_flags", "pluralia_tantum")] += 1
            if lexeme.head.by_sense:
                table[("noun_in_domain_flags", "by_sense")] += 1
            if lexeme.second:
                table[("noun_in_domain_flags", "has_second_source")] += 1
        for reason, count in lexeme.paradigm.dropped.items():
            table[("dropped_forms", reason)] += count
        if lexeme.domain == DOMAIN_IN and lexeme.pos in CS_POS:
            for cell in {c.cell for c in lexeme.paradigm.cells}:
                table[("cell_coverage", f"{lexeme.pos}:{cell}")] += 1
    table[("cs_entries", "matched")] = stats.matched
    table[("cs_entries", "cs_only")] = stats.cs_only
    table[("cs_entries", "names")] = stats.names
    table[("cs_entries", "unmatched")] = stats.unmatched
    lines = ["section\tkey\tvalue"]
    lines.extend(f"{s}\t{k}\t{v}" for (s, k), v in sorted(table.items()))
    ctx.write_text(paths.lexicon_stats, "\n".join(lines) + "\n")
