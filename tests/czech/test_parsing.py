"""Unit tests: head codes, cs-ndecl specs, R13 cells and R25 names."""

from __future__ import annotations

import json
from typing import Any

import pytest

from lexpipe.langs.czech import cells, declspec, heads, names
from tests.czech import treebuilder

pytestmark = pytest.mark.cs


def _entries(
    path: Any, word: str, pos: str | None = None
) -> list[dict[str, Any]]:
    out = []
    for line in path.read_text(encoding="utf-8").splitlines():
        entry = json.loads(line)
        if entry["word"] == word and (pos is None or entry["pos"] == pos):
            out.append(entry)
    return out


def en(word: str, pos: str | None = None) -> list[dict[str, Any]]:
    return _entries(treebuilder.EN_FIXTURE, word, pos)


def cs(word: str, pos: str | None = None) -> list[dict[str, Any]]:
    return _entries(treebuilder.CS_FIXTURE, word, pos)


# -- head codes ---------------------------------------------------------------


def test_head_code_single_reading() -> None:
    info = heads.parse_head_code("m-an")
    assert info.labels == ("m.an",)
    assert info.genders == ("m",) and info.animacies == ("an",)
    assert not info.animacy_ambiguous and not info.pluralia_tantum


def test_head_code_comma_list_gives_both_readings() -> None:
    info = heads.parse_head_code("m-an,m-in")
    assert info.labels == ("m.an", "m.in")
    assert info.animacy_ambiguous and not info.gender_ambiguous


def test_head_code_qualifier_commas_do_not_split() -> None:
    info = heads.parse_head_code("m-in,m-an<l:archaic,or,literary>")
    assert info.labels == ("m.in", "m.an")
    assert info.readings[1].qualifier == "<l:archaic,or,literary>"


@pytest.mark.parametrize(
    ("code", "tantum", "gender_ambiguous"),
    [
        ("f-p", True, False),
        ("m-in-p,f-p", True, True),
        ("f,m-in<l:rare>", False, True),
    ],
)
def test_head_code_pluralia_tantum(
    code: str, tantum: bool, gender_ambiguous: bool
) -> None:
    info = heads.parse_head_code(code)
    assert info.pluralia_tantum is tantum
    assert info.gender_ambiguous is gender_ambiguous


def test_head_code_by_sense_and_unknown() -> None:
    info = heads.parse_head_code("mfbysense-an")
    assert info.by_sense and info.genders == ("m", "f")
    assert info.animacies == ("an",)
    unknown = heads.parse_head_code("?-p")
    assert unknown.genders == ("?",) and unknown.pluralia_tantum


def test_readings_from_cs_tags() -> None:
    assert (
        heads.readings_from_tags(["animate", "masculine"])[0].label == "m.an"
    )
    assert heads.readings_from_tags(["masculine"])[0].label == "m"
    assert heads.readings_from_tags(["feminine"])[0].label == "f"


def test_morfflex_gender_letters() -> None:
    animate = heads.morfflex_genders(heads.GenderReading("m", "an"))
    assert animate == {"M"}
    assert heads.morfflex_gender_matches("Y", animate)
    assert not heads.morfflex_gender_matches("I", animate)
    assert not heads.morfflex_gender_matches("H", animate)
    assert heads.morfflex_genders(heads.GenderReading("m")) == {"M", "I"}


def test_head_info_of_fixture_entries() -> None:
    (agent_ambiguous,) = [
        e for e in en("agent") if heads.head_info(e).animacy_ambiguous
    ]
    assert heads.head_info(agent_ambiguous).code == "m-an,m-in"
    (shah,) = [e for e in en("šach") if heads.head_info(e).labels == ("m.an",)]
    assert shah["senses"][0]["glosses"] == ["shah"]


# -- cs-ndecl specs ------------------------------------------------------------


@pytest.mark.parametrize(
    ("raw", "overrides"),
    [
        ("m.an.nompli:ové", {("NOM", "PL")}),
        ("m.gena.locu:e", {("GEN", "SG"), ("LOC", "SG")}),
        (
            "f.genplou:-.locplou:ách.insplama",
            {("GEN", "PL"), ("LOC", "PL"), ("INS", "PL")},
        ),
        ("m.locplách", {("LOC", "PL")}),
        ("m.an.nomplé:i", {("NOM", "PL")}),
        ("m.an", set()),
        ("n.sg.+.insým.[rare]", {("INS", "SG")}),
    ],
)
def test_spec_overrides(raw: str, overrides: set[tuple[str, str]]) -> None:
    assert declspec.parse_decl_spec(raw).overrides == overrides


def test_spec_flags_and_special_forms() -> None:
    surname = declspec.parse_decl_spec("m.an.surname")
    assert surname.has_flag("surname") and surname.overrides == frozenset()
    combined = declspec.parse_decl_spec("Kačer<m.an.dat+locu.nompli>")
    assert combined.overrides == {("DAT", "SG"), ("LOC", "SG"), ("NOM", "PL")}
    ucho = declspec.parse_decl_spec(
        "((n.sg,uši<n.pl.decllemma:uší.inspl:ušima>))"
    )
    assert ucho.alternants == 2
    assert ucho.overrides == {("INS", "PL")}
    assert {"altlemma", "decllemma", "pl", "sg"} <= set(ucho.flags)
    assert declspec.parse_decl_spec("Česká<+> republika<f.sg>").multiword
    manual = declspec.parse_decl_spec(
        "m.!\n .gen:téhodne\n .nompl:tédni:tídnové"
    )
    assert manual.manual and ("GEN", "SG") in manual.overrides
    assert (
        "[rare]" not in declspec.parse_decl_spec("n.sg.+.insým.[rare]").flags
    )
    assert declspec.parse_decl_spec("n.tstem").has_flag("tstem")


def test_explicit_cells_follow_module_copies() -> None:
    nompl = declspec.parse_decl_spec("m.locplách.nomplové")
    inanimate = declspec.explicit_cells(nompl, masculine_animate=False)
    assert {("VOC", "PL"), ("ACC", "PL"), ("LOC", "PL")} <= inanimate
    animate = declspec.explicit_cells(
        declspec.parse_decl_spec("m.an.gena.nomplové"), masculine_animate=True
    )
    assert ("ACC", "SG") in animate and ("ACC", "PL") not in animate


# -- cells -------------------------------------------------------------------


@pytest.mark.parametrize(
    ("tags", "adjective", "expected"),
    [
        (["nominative", "plural"], False, ("NOM.PL", "")),
        (["dual", "nominative"], False, (None, cells.DROP_DUAL)),
        (
            ["alternative", "dative", "locative", "singular"],
            False,
            (None, cells.DROP_MULTI_CASE),
        ),
        (["feminine", "nominative"], True, (None, cells.DROP_NO_NUMBER)),
        (
            ["animate", "masculine", "nominative", "plural"],
            True,
            ("NOM.PL.M.AN", ""),
        ),
        (["genitive", "plural"], True, ("GEN.PL", "")),
        (["genitive", "plural", "uncommon"], False, ("GEN.PL", "")),
        (["feminine"], False, (None, cells.DROP_NO_CASE)),
    ],
)
def test_cell_label(
    tags: list[str], adjective: bool, expected: tuple[str | None, str]
) -> None:
    assert cells.cell_label(tags, adjective=adjective) == expected


def test_cell_order() -> None:
    shuffled = ["LOC.PL", "NOM.SG", "INS.SG", "NOM.PL", "GEN.SG"]
    assert sorted(shuffled, key=cells.cell_sort_key) == [
        "NOM.SG",
        "GEN.SG",
        "INS.SG",
        "NOM.PL",
        "LOC.PL",
    ]
    assert cells.NOUN_CELLS[0] == "NOM.SG" and cells.NOUN_CELLS[-1] == "INS.PL"


def test_en_cells_segment_tables() -> None:
    (two_tables,) = [
        e for e in en("agent") if len(e["inflection_templates"]) == 2
    ]
    paradigm = cells.en_cells(two_tables, adjective=False)
    assert [animacy for animacy, _ in paradigm.tables] == ["an", "in"]
    nominative = [
        (c.table, c.form) for c in paradigm.cells if c.cell == "NOM.PL"
    ]
    assert nominative == [(0, "agenti"), (1, "agenty")]
    assert all(
        c.table_animacy == ("an" if c.table == 0 else "in")
        for c in paradigm.cells
    )


def test_en_cells_drop_dual_and_headline_forms() -> None:
    (ruka,) = en("ruka", "noun")
    paradigm = cells.en_cells(ruka, adjective=False)
    assert paradigm.dropped[cells.DROP_DUAL] > 0
    assert {c.form for c in paradigm.cells if c.cell == "NOM.PL"} == {"ruce"}
    assert "ruční" not in {c.form for c in paradigm.cells}  # head-line form


def test_syncretic_forms_once_per_cell() -> None:
    (kluk,) = en("kluk", "noun")
    paradigm = cells.en_cells(kluk, adjective=False)
    sync = cells.syncretism(paradigm.cells)
    assert sync[(0, "kluku")] == ("DAT.SG", "VOC.SG", "LOC.SG")
    kluku = [c.cell for c in paradigm.cells if c.form == "kluku"]
    assert sorted(kluku) == sorted(set(kluku))


def test_cs_cells_keep_qualifiers() -> None:
    (kluk,) = cs("kluk", "noun")
    paradigm = cells.cs_cells(kluk, adjective=False, animacy="an")
    nominative = {
        c.form: c.qualifiers for c in paradigm.cells if c.cell == "NOM.PL"
    }
    assert nominative == {"kluci": (), "klucí": ("dialectal",)}
    assert "NOM.SG" not in {c.cell for c in paradigm.cells}  # omitted by cs


# -- names (R25) ---------------------------------------------------------------


def test_ethnonym_stays_and_surname_goes() -> None:
    (noun,) = en("Čech", "noun")
    (name,) = en("Čech", "name")
    assert names.classify_en(noun) == names.COMMON
    verdict = names.classify_en(name)
    assert verdict.is_name and verdict.name_class == names.SURNAME
    assert verdict.evidence == "pos"


def test_given_names_under_a_noun_header() -> None:
    verdicts = sorted(
        (names.classify_en(e).is_name, names.classify_en(e).name_class)
        for e in en("Inka")
    )
    assert verdicts == [
        (False, ""),
        (True, "given_name"),
        (True, "given_name"),
    ]


def test_toponym_from_gloss() -> None:
    entry = {
        "word": "Egypt",
        "pos": "name",
        "senses": [
            {"glosses": ["Egypt (a country in North Africa and West Asia)"]}
        ],
    }
    assert names.classify_en(entry).name_class == names.TOPONYM


def test_cs_edition_names() -> None:
    classes = sorted(names.classify_cs(e).name_class for e in cs("Čech"))
    # ethnonym (common), surname, mythical founder, form-of stub
    assert classes.count(names.SURNAME) == 1
    assert "" in classes and names.OTHER_NAME in classes
