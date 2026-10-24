"""Paradigm cells of Czech nouns and adjectives (design rule R13).

A cell is a case × number slot, ``NOM.SG`` … ``INS.PL``; adjectives add
the gender (and masculine animacy) their forms are tagged with, so the
hard-adjective cell of §6 is ``NOM.PL.M.AN``. Cells are read only from
wiktextract ``forms`` that carry **exactly one case tag and exactly one
number tag**:

  - a dual form (``ruce`` ["dual", "nominative"]) has a number tag but
    no SG/PL cell, so it is dropped and counted as ``dual``;
  - merged adjective cells that lost their number (``české``
    ["feminine", "nominative"]) are dropped as ``no_number``;
  - forms with two case tags (``alternative, dative, locative``) are
    dropped as ``multi_case``.

Every other tag of a kept form (``rare``, ``archaic``, ``dialectal`` …)
and its raw tags are kept as ``qualifiers``; nothing is ranked
("variants unranked"). A form is recorded once per (table, cell): a
repeat (cs entries that print an animate and an inanimate table
without a separator, alternant declensions) only adds qualifiers.
``order`` is the form's position in the entry, recorded only so that
the source order can be reproduced.

English Wiktionary tables
-------------------------
An en entry may hold several declension tables (agent: an animate and
an inanimate one; týden: a regular and a manual one). wiktextract
starts each table with two pseudo-forms, ``table-tags`` (whose form is
``animate``, ``inanimate`` or ``no-table-tags``) and
``inflection-template``; a ``class`` pseudo-form (``i-stem``,
``t-stem``, ``irreg-stem`` …) may follow. :func:`en_cells` numbers the
tables from 0 in that order, so table *i* belongs to
``inflection_templates[i]`` when the counts agree. Only forms from a
table (``source == "declension"``) are cells: forms on the headword
line (feminine equivalent, diminutive, relational adjective) are other
lexemes.

Czech Wiktionary tables
-----------------------
The cs edition (:func:`cs_cells`) has one literal table per entry and
omits cells identical to the headword (NOM.SG; ACC.SG of inanimates).
Its forms carry no ``source``.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Final

from lexpipe import kaikki
from lexpipe.normalize import normalize_orthography

CASES: Final[tuple[tuple[str, str], ...]] = (
    ("nominative", "NOM"),
    ("genitive", "GEN"),
    ("dative", "DAT"),
    ("accusative", "ACC"),
    ("vocative", "VOC"),
    ("locative", "LOC"),
    ("instrumental", "INS"),
)
CASE_TAGS: Final[Mapping[str, str]] = dict(CASES)
CASE_ORDER: Final[Mapping[str, int]] = {
    code: index for index, (_, code) in enumerate(CASES)
}
NUMBER_TAGS: Final[Mapping[str, str]] = {"singular": "SG", "plural": "PL"}
DUAL_TAG: Final[str] = "dual"
NUMBER_ORDER: Final[Mapping[str, int]] = {"SG": 0, "PL": 1}
GENDER_TAGS: Final[Mapping[str, str]] = {
    "masculine": "M",
    "feminine": "F",
    "neuter": "N",
}
ANIMACY_TAGS: Final[Mapping[str, str]] = {"animate": "AN", "inanimate": "IN"}
_GENDER_SUFFIX_ORDER: Final[tuple[str, ...]] = (
    "",
    ".M.AN",
    ".M.IN",
    ".M",
    ".F",
    ".N",
)

NOUN_CELLS: Final[tuple[str, ...]] = tuple(
    f"{case}.{number}" for number in ("SG", "PL") for _, case in CASES
)
HARD_ADJECTIVE_CELL: Final[str] = "NOM.PL.M.AN"

TABLE_TAGS: Final[str] = "table-tags"
TEMPLATE_TAGS: Final[str] = "inflection-template"
CLASS_TAGS: Final[str] = "class"
_TABLE_ANIMACY: Final[Mapping[str, str]] = {
    "animate": "an",
    "inanimate": "in",
}
DECLENSION_SOURCE: Final[str] = "declension"

DROP_DUAL: Final[str] = "dual"
DROP_NO_CASE: Final[str] = "no_case"
DROP_MULTI_CASE: Final[str] = "multi_case"
DROP_NO_NUMBER: Final[str] = "no_number"
DROP_MULTI_NUMBER: Final[str] = "multi_number"
DROP_NOT_TABLE: Final[str] = "not_table"


@dataclass(frozen=True, slots=True)
class CellForm:
    """One form filling one cell (see the module docstring)."""

    cell: str
    form: str
    table: int = 0
    table_animacy: str = ""
    order: int = 0
    qualifiers: tuple[str, ...] = ()
    tag_error: str = ""


@dataclass(frozen=True, slots=True)
class Paradigm:
    """The cells of one entry plus what was dropped and the tables seen.

    ``tables`` lists ``(animacy, stem_class)`` per table, in order.
    """

    cells: tuple[CellForm, ...] = ()
    dropped: Mapping[str, int] = field(default_factory=dict)
    tables: tuple[tuple[str, str], ...] = ()


def cell_label(
    tags: Iterable[str], *, adjective: bool
) -> tuple[str | None, str]:
    """The cell of a form's tags, or ``(None, reason)``.

    Nouns ignore gender tags; adjectives append ``.M``/``.F``/``.N`` and,
    for masculines, ``.AN``/``.IN`` when tagged.
    """
    tag_set = set(tags)
    cases = [CASE_TAGS[t] for t in tag_set if t in CASE_TAGS]
    numbers = [NUMBER_TAGS[t] for t in tag_set if t in NUMBER_TAGS]
    dual = DUAL_TAG in tag_set
    if not cases:
        return None, DROP_NO_CASE
    if len(cases) > 1:
        return None, DROP_MULTI_CASE
    if dual and not numbers:
        return None, DROP_DUAL
    if not numbers:
        return None, DROP_NO_NUMBER
    if len(numbers) > 1 or dual:
        return None, DROP_MULTI_NUMBER
    cell = f"{cases[0]}.{numbers[0]}"
    if adjective:
        genders = [GENDER_TAGS[t] for t in GENDER_TAGS if t in tag_set]
        if len(genders) == 1:
            cell += f".{genders[0]}"
            animacies = [ANIMACY_TAGS[t] for t in ANIMACY_TAGS if t in tag_set]
            if genders[0] == "M" and len(animacies) == 1:
                cell += f".{animacies[0]}"
        elif genders:
            # A merged cell (two genders share the form) keeps both.
            cell += "." + "+".join(genders)
    return cell, ""


def cell_sort_key(cell: str) -> tuple[int, int, int, str]:
    """Canonical order: case order within number (SG first), then the
    adjective gender suffix."""
    parts = cell.split(".")
    case = CASE_ORDER.get(parts[0], len(CASE_ORDER))
    number = NUMBER_ORDER.get(parts[1], 2) if len(parts) > 1 else 2
    suffix = "".join(f".{p}" for p in parts[2:])
    rank = (
        _GENDER_SUFFIX_ORDER.index(suffix)
        if suffix in _GENDER_SUFFIX_ORDER
        else len(_GENDER_SUFFIX_ORDER)
    )
    return (number, case, rank, cell)


def _strings(value: Any) -> tuple[str, ...]:
    return tuple(item for item in value or () if isinstance(item, str))


def _qualifiers(
    tags: Sequence[str], raw_tags: Sequence[str]
) -> tuple[str, ...]:
    structural = (
        set(CASE_TAGS)
        | set(NUMBER_TAGS)
        | set(GENDER_TAGS)
        | set(ANIMACY_TAGS)
        | {DUAL_TAG}
        | kaikki.TAG_ERROR_TAGS
    )
    found = {t for t in tags if t not in structural}
    found.update(t.strip() for t in raw_tags if t.strip())
    return tuple(sorted(found))


def _form_text(form: Mapping[str, Any]) -> str | None:
    text = form.get("form")
    if not isinstance(text, str) or text.strip() in kaikki.PLACEHOLDER_FORMS:
        return None
    return normalize_orthography(text)


def en_cells(entry: Mapping[str, Any], *, adjective: bool) -> Paradigm:
    """Cells of an English Wiktionary entry (see the module docstring)."""
    cells: list[CellForm] = []
    seen: dict[tuple[int, str, str], int] = {}
    dropped: Counter[str] = Counter()
    tables: list[list[str]] = []
    for order, form in enumerate(entry.get("forms") or ()):
        if not isinstance(form, Mapping):
            continue
        tags = _strings(form.get("tags"))
        text = _form_text(form)
        if TABLE_TAGS in tags:
            tables.append([_TABLE_ANIMACY.get(text or "", ""), ""])
            continue
        if CLASS_TAGS in tags:
            if tables and text:
                tables[-1][1] = text
            continue
        if kaikki.PSEUDO_FORM_TAGS.intersection(tags) or text is None:
            continue
        if form.get("source") != DECLENSION_SOURCE:
            if set(tags) & set(CASE_TAGS):
                dropped[DROP_NOT_TABLE] += 1
            continue
        cell, reason = cell_label(tags, adjective=adjective)
        if cell is None:
            dropped[reason] += 1
            continue
        index = max(len(tables) - 1, 0)
        _add(
            cells,
            seen,
            CellForm(
                cell=cell,
                form=text,
                table=index,
                table_animacy=tables[index][0] if tables else "",
                order=order,
                qualifiers=_qualifiers(tags, _strings(form.get("raw_tags"))),
                tag_error="|".join(
                    sorted(kaikki.TAG_ERROR_TAGS.intersection(tags))
                ),
            ),
        )
    return Paradigm(
        cells=tuple(cells),
        dropped=dict(sorted(dropped.items())),
        tables=tuple((animacy, stem) for animacy, stem in tables),
    )


def cs_cells(
    entry: Mapping[str, Any], *, adjective: bool, animacy: str = ""
) -> Paradigm:
    """Cells of a Czech Wiktionary entry (one literal table)."""
    cells: list[CellForm] = []
    seen: dict[tuple[int, str, str], int] = {}
    dropped: Counter[str] = Counter()
    for order, form in enumerate(entry.get("forms") or ()):
        if not isinstance(form, Mapping):
            continue
        tags = _strings(form.get("tags"))
        text = _form_text(form)
        if text is None or kaikki.PSEUDO_FORM_TAGS.intersection(tags):
            continue
        if not set(tags) & set(CASE_TAGS):
            continue
        cell, reason = cell_label(tags, adjective=adjective)
        if cell is None:
            dropped[reason] += 1
            continue
        _add(
            cells,
            seen,
            CellForm(
                cell=cell,
                form=text,
                table=0,
                table_animacy=animacy,
                order=order,
                qualifiers=_qualifiers(tags, _strings(form.get("raw_tags"))),
            ),
        )
    return Paradigm(
        cells=tuple(cells),
        dropped=dict(sorted(dropped.items())),
        tables=((animacy, ""),),
    )


def _add(
    cells: list[CellForm], seen: dict[tuple[int, str, str], int], new: CellForm
) -> None:
    """Append *new* unless its (table, cell, form) is already there; a
    repeat (cs tables printed twice, alternant declensions) only adds
    its qualifiers."""
    key = (new.table, new.cell, new.form)
    index = seen.get(key)
    if index is None:
        seen[key] = len(cells)
        cells.append(new)
        return
    old = cells[index]
    merged = tuple(sorted(set(old.qualifiers) | set(new.qualifiers)))
    if merged != old.qualifiers:
        cells[index] = CellForm(
            old.cell,
            old.form,
            old.table,
            old.table_animacy,
            old.order,
            merged,
            old.tag_error,
        )


def syncretism(
    cells: Iterable[CellForm],
) -> dict[tuple[int, str], tuple[str, ...]]:
    """``(table, form) -> cells`` holding that form, in canonical order.

    A form is recorded once per cell it fills; this map lets a row name
    the other cells of its table that share the form (kluku: DAT.SG,
    VOC.SG, LOC.SG).
    """
    found: dict[tuple[int, str], set[str]] = {}
    for cell in cells:
        found.setdefault((cell.table, cell.form), set()).add(cell.cell)
    return {
        key: tuple(sorted(value, key=cell_sort_key))
        for key, value in found.items()
    }
