"""Per-form plural provenance (lexpipe.langs.italian.plurals, R12)."""

from __future__ import annotations

from typing import Any

import pytest

from lexpipe.langs.italian import plurals

pytestmark = pytest.mark.it


def _noun(args: dict[str, str]) -> list[dict[str, Any]]:
    return [{"name": "it-noun", "args": args}]


def _adj(args: dict[str, str]) -> list[dict[str, Any]]:
    return [{"name": "it-adj", "args": args}]


def _rows(info: plurals.PluralInfo) -> list[tuple[str, str, str, str]]:
    return [(r.cell, r.form, r.provenance, r.code) for r in info.rows]


def test_pure_default_noun() -> None:
    info = plurals.resolve_plurals("lago", _noun({"1": "m"}))
    assert _rows(info) == [("pl", "laghi", "module_default", "default")]
    assert info.pure_default and info.provenance == "module_default"
    assert info.status == "regular"


def test_mago_default_plus_labelled_variant() -> None:
    info = plurals.resolve_plurals(
        "mago",
        _noun(
            {
                "1": "m",
                "2": "+,magi<q:archaic, literary, dialectal or "
                "Zoroastrianism>",
            }
        ),
    )
    assert _rows(info) == [
        ("pl", "maghi", "module_default", "+"),
        ("pl", "magi", "explicit", ""),
    ]
    assert info.rows[1].labels == ("archaic", "dialectal", "literary")
    assert info.rows[1].raw_arg.startswith("2=magi<q:")
    assert info.provenance == "mixed"


def test_fuoco_old_italian_and_tag_labels() -> None:
    info = plurals.resolve_plurals(
        "fuoco",
        _noun({"1": "m", "2": "fuochi,fòcora<l:Old Italian>,fuòcora"}),
    )
    assert info.forms == ("fuochi", "fòcora", "fuòcora")
    assert info.rows[1].labels == ("old italian",)
    assert {r.provenance for r in info.rows} == {"explicit"}


def test_valigia_explicit_then_default() -> None:
    info = plurals.resolve_plurals(
        "valigia", _noun({"1": "f", "2": "valigie,+"})
    )
    assert _rows(info) == [
        ("pl", "valigie", "explicit", ""),
        ("pl", "valige", "module_default", "+"),
    ]


def test_belga_gendered_cells() -> None:
    info = plurals.resolve_plurals(
        "belga",
        _noun({"1": "mfbysense", "mpl": "belgi", "fpl": "belghe"}),
    )
    assert _rows(info) == [
        ("m.pl", "belgi", "explicit", ""),
        ("f.pl", "belghe", "explicit", ""),
    ]
    collega = plurals.resolve_plurals("collega", _noun({"1": "mfbysense"}))
    assert _rows(collega) == [
        ("m.pl", "colleghi", "module_default", "default"),
        ("f.pl", "colleghe", "module_default", "default"),
    ]


def test_codes() -> None:
    greco = plurals.resolve_plurals("greco", _noun({"1": "m", "2": "-"}))
    assert greco.status == "uncountable" and greco.forms == ()
    specie = plurals.resolve_plurals("specie", _noun({"1": "f", "2": "#"}))
    assert specie.forms == ("specie",) and specie.invariable
    peso = plurals.resolve_plurals("peso", _noun({"1": "m", "2": "#s"}))
    assert peso.forms == ("pesos",)
    latte = plurals.resolve_plurals("latte", _noun({"1": "m", "2": "~"}))
    assert latte.forms == ("latti",)
    assert latte.countability == "countable_and_uncountable"
    uva = plurals.resolve_plurals("uva", _noun({"1": "f", "2": "-,+"}))
    assert uva.forms == ("uve",) and "uncountable" in uva.countability
    letto = plurals.resolve_plurals(
        "divano-letto", _noun({"1": "m", "2": "++"})
    )
    assert letto.forms == ("divani-letti",)
    first = plurals.resolve_plurals(
        "anno luce", _noun({"1": "m", "2": "+first"})
    )
    assert first.rows[0].tag_error == "not_ported"


def test_capo_by_sense_gives_both_cells() -> None:
    info = plurals.resolve_plurals(
        "capocronista", _noun({"1": "mfbysense", "2": "cap*+"})
    )
    assert _rows(info) == [
        ("m.pl", "capicronisti", "module_default", "cap*+"),
        ("f.pl", "capocroniste", "module_default", "cap*+"),
    ]


def test_eco_mixed_gender_items() -> None:
    info = plurals.resolve_plurals(
        "eco",
        _noun(
            {
                "1": "f,m<l:colloquial,or,literary>",
                "2": "+<g:m>,#<g:f><l:<<proscribed>>>",
            }
        ),
    )
    assert [(r.form, r.gender, r.code) for r in info.rows] == [
        ("echi", "m", "+"),
        ("eco", "f", "#"),
    ]
    assert [g.code for g in info.genders] == ["f", "m"]


def test_plurale_tantum_has_no_rows() -> None:
    info = plurals.resolve_plurals("forbici", _noun({"1": "f-p"}))
    assert info.rows == () and info.status == "plurale_tantum"


def test_adjectives() -> None:
    carico = plurals.resolve_plurals("carico", _adj({"mpl": "carichi"}))
    assert _rows(carico) == [
        ("m.pl", "carichi", "explicit", ""),
        ("f.pl", "cariche", "module_default", "default"),
    ]
    antico = plurals.resolve_plurals("antico", _adj({}))
    assert antico.cell_forms("m.pl") == ("antici",)
    inv = plurals.resolve_plurals("blu", _adj({"inv": "1"}))
    assert _rows(inv) == [("pl", "blu", "explicit", "#")]
    assert inv.rows[0].raw_arg == "inv=1"


def test_several_heads_are_merged() -> None:
    info = plurals.resolve_plurals(
        "otorinolaringoiatra",
        _noun({"1": "m"}) + _noun({"1": "f", "2": "otorinolaringoiatre"}),
    )
    assert info.forms == ("otorinolaringoiatri", "otorinolaringoiatre")
    assert info.rows[1].raw_arg == "h2:2=otorinolaringoiatre"
    assert not info.pure_default


def test_no_head_template() -> None:
    info = plurals.resolve_plurals("x", [{"name": "head", "args": {}}])
    assert info.template == "" and info.status == "no_head_template"


def _entry(word: str, forms: list[tuple[str, list[str]]]) -> dict[str, Any]:
    return {
        "word": word,
        "forms": [{"form": f, "tags": tags} for f, tags in forms],
    }


def test_drift_statuses() -> None:
    lago = plurals.resolve_plurals("lago", _noun({"1": "m"}))
    match = plurals.drift(lago, _entry("lago", [("laghi", ["plural"])]))
    assert match is not None
    assert (match.scope, match.status) == ("pure_default_noun", "match")
    wrong = plurals.drift(lago, _entry("lago", [("lagi", ["plural"])]))
    assert wrong is not None and wrong.status == "mismatch"
    assert wrong.resolved == ("pl:laghi",) and wrong.rendered == ("pl:lagi",)
    # kaikki renders an invariable plural as "(invariable)", not a form.
    specie = plurals.resolve_plurals("specie", _noun({"1": "f", "2": "#"}))
    invariable = plurals.drift(specie, _entry("specie", []))
    assert invariable is not None and invariable.status == "match"
    multi = plurals.resolve_plurals("abete rosso", _noun({"1": "m"}))
    not_ported = plurals.drift(
        multi, _entry("abete rosso", [("abeti rossi", ["plural"])])
    )
    assert not_ported is not None and not_ported.status == "not_ported"
    tantum = plurals.resolve_plurals("forbici", _noun({"1": "f-p"}))
    assert plurals.drift(tantum, _entry("forbici", [])) is None
