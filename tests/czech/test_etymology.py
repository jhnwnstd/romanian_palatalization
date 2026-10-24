"""Etymological strata (R22): Wiktionary templates, czetyl, combination."""

from __future__ import annotations

from pathlib import Path

import pytest

from lexpipe.langs.czech import etymology as et

pytestmark = pytest.mark.cs


def ety(*args: str, **named: str) -> dict[str, object]:
    values = {"1": "cs", **{str(i): a for i, a in enumerate(args, 2)}}
    values.update(named)
    return {"name": "ety", "args": values}


def test_affix_step_with_modifiers() -> None:
    found = et.steps([ety(":af", "kniha<t:book>", "-ní")])
    assert [(s.kind, s.lang, s.terms) for s in found] == [
        ("af", "cs", ("kniha", "-ní"))
    ]
    assert et.affix_parts(found) == [("kniha", "-ní")]
    assert et.wiktionary_stratum(found, "") == (et.INTERNAL, "af")


def test_inheritance_with_a_surface_analysis_after_it() -> None:
    found = et.steps(
        [ety(":inh", "sla-pro:*strašiti", ":af<surf>", "strach", "-it")]
    )
    assert [s.kind for s in found] == ["inh", "surf"]
    assert et.wiktionary_stratum(found, "")[0] == et.INHERITED
    assert et.affix_parts(found) == [("strach", "-it")]


def test_surface_analysis_first_is_borrowed_whole() -> None:
    found = et.steps([ety(":af<surf>", "Egypt", "-an")])
    assert et.wiktionary_stratum(found, "")[0] == et.BORROWED


def test_duplicate_template_artifacts_are_skipped() -> None:
    clean = ety(":inh", "zlw-ocs:kniha")
    artifact = ety(":inh", "['zlw-ocs:', <LINK(['kniha']){} >]")
    found = et.steps([clean, artifact, clean])
    assert [(s.kind, s.terms) for s in found] == [("inh", ("kniha",))]


@pytest.mark.parametrize(
    ("template", "stratum"),
    [
        (
            {"name": "bor", "args": {"1": "cs", "2": "it", "3": "banca"}},
            et.BORROWED,
        ),
        (
            {
                "name": "inh+",
                "args": {"1": "cs", "2": "zlw-ocs", "3": "prach"},
            },
            et.INHERITED,
        ),
        (
            {"name": "der", "args": {"1": "cs", "2": "ar", "3": "-"}},
            et.BORROWED,
        ),
        (
            {"name": "der", "args": {"1": "cs", "2": "sla-pro", "3": "*x"}},
            et.INHERITED,
        ),
        (
            {"name": "der", "args": {"1": "cs", "2": "cs", "3": "x"}},
            et.INTERNAL,
        ),
        (
            {"name": "suffix", "args": {"1": "cs", "2": "Bělov", "3": "ský"}},
            et.INTERNAL,
        ),
    ],
)
def test_classic_templates(template: dict[str, object], stratum: str) -> None:
    assert et.wiktionary_stratum(et.steps([template]), "")[0] == stratum


@pytest.mark.parametrize(
    ("text", "stratum"),
    [
        ("From sníh (“snow”, noun).", et.INTERNAL),
        ("Inherited from Old Czech prach.", et.INHERITED),
        ("Borrowed from German Liga.", et.BORROWED),
        ("From Latin liga.", et.BORROWED),
        ("", et.UNKNOWN),
    ],
)
def test_text_fallback(text: str, stratum: str) -> None:
    assert et.wiktionary_stratum([], text)[0] == stratum


def test_non_mapping_template_items_are_skipped() -> None:
    assert et.steps(["not a template", 3, None]) == []


def test_czetyl(tmp_path: Path) -> None:
    path = tmp_path / "czetyl.tsv"
    path.write_text(
        "word\tetymology\tis_loan\n"
        "a\t\tnative\n"
        "a-2\tell\tloan\n"
        "liga\tdeu\tloan\n"
        "kniha\t\tnative\n",
        encoding="utf-8",
    )
    table = et.load_czetyl(path)
    assert table == {"kniha": et.INHERITED, "liga": et.BORROWED}


def test_czetyl_rejects_another_format(tmp_path: Path) -> None:
    path = tmp_path / "czetyl.tsv"
    path.write_text("lemma\tx\n", encoding="utf-8")
    with pytest.raises(ValueError):
        et.load_czetyl(path)


def test_combine_precedence() -> None:
    assert et.combine((et.INHERITED, "inh"), et.BORROWED, True) == et.Stratum(
        et.INHERITED, et.SOURCE_WIKTIONARY, "inh"
    )
    assert et.combine((et.UNKNOWN, ""), et.BORROWED, None).source == (
        et.SOURCE_CZETYL
    )
    assert et.combine((et.UNKNOWN, ""), None, True) == et.Stratum(
        et.BORROWED, et.SOURCE_DERINET, "Loanword"
    )
    assert et.combine((et.UNKNOWN, ""), None, None).stratum == et.UNKNOWN
