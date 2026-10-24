"""Head-template argument grammar (lexpipe.langs.italian.heads, R12)."""

from __future__ import annotations

import pytest

from lexpipe.langs.italian import heads

pytestmark = pytest.mark.it


def test_split_items_respects_brackets() -> None:
    value = (
        "braccia<g:f><q:in most meanings 1, 2, 4>,"
        'bracci<g:m><q:in the figurative meanings "branch, wing" 6, 7>'
    )
    assert [item.split("<")[0] for item in heads.split_items(value)] == [
        "braccia",
        "bracci",
    ]
    assert heads.split_items(" a , ,b ") == ["a", "b"]
    assert heads.split_items("") == []


def test_parse_item_modifiers_and_codes() -> None:
    item = heads.parse_item(
        "+,magi<q:archaic, literary, dialectal or Zoroastrianism>".split(
            ",", 1
        )[1]
    )
    assert item.text == "magi"
    assert item.code == ""
    assert item.qualifiers == (
        "q:archaic, literary, dialectal or Zoroastrianism",
    )
    assert item.labels == ("archaic", "dialectal", "literary")
    plus = heads.parse_item("+<g:m>")
    assert plus.code == "+" and plus.gender == "m"
    assert heads.parse_item("#s").code == ""
    assert heads.parse_item("cap*+").code == "cap*+"
    assert heads.parse_item("+first").code == "+first"


def test_nested_reference_modifiers() -> None:
    raw = (
        "#<g:f><l:<<proscribed>><ref name=dop/>><ref:“eco” in "
        "garzantilinguistica.it – Garzanti Linguistica, De Agostini>"
    )
    item = heads.parse_item(raw)
    assert item.error == ""
    assert item.code == "#"
    assert item.gender == "f"
    assert item.labels == ("proscribed",)
    # References never reach the qualifiers.
    assert item.qualifiers == ("l:proscribed",)


def test_unbalanced_item_is_flagged() -> None:
    item = heads.parse_item("+<l:chiefly,archaic,or,idiomatic>>")
    assert item.error == heads.ERROR_UNBALANCED
    assert item.code == "+"
    assert item.labels == ("archaic", "chiefly", "idiomatic")
    assert heads.parse_item("x<nokey>").error == heads.ERROR_BAD_MODIFIER


def test_parse_gender() -> None:
    codes = heads.parse_gender(
        "f,m<l:colloquial,or,literary><ref:<<name:dop>>>"
    )
    assert [c.code for c in codes] == ["f", "m"]
    assert codes[1].labels == ("colloquial", "literary")
    tantum = heads.parse_gender("m-p")
    assert tantum[0].base == "m" and tantum[0].plural_only
    assert heads.parse_gender("mfbysense")[0].base == "mfbysense"
    assert heads.parse_gender("") == ()


def test_plural_slots_order_and_cells() -> None:
    args = {"1": "m", "fpl": "dee", "2": "dei", "pl2": "dii", "mpl": "x"}
    slots = heads.plural_slots(heads.NOUN_TEMPLATE, args)
    assert [(s.key, s.cell) for s in slots] == [
        ("2", "pl"),
        ("pl2", "pl"),
        ("mpl", "m.pl"),
        ("fpl", "f.pl"),
    ]
    adjective = heads.plural_slots(
        heads.ADJ_TEMPLATE, {"mpl": "antichi", "sup": "x", "pl": "y"}
    )
    assert [(s.key, s.cell) for s in adjective] == [
        ("pl", "pl"),
        ("mpl", "m.pl"),
    ]
    # An empty argument is absent, as in the module.
    assert heads.plural_slots(heads.NOUN_TEMPLATE, {"1": "m", "2": ""}) == []


def test_head_template_picks_italian_heads() -> None:
    entry = {
        "head_templates": [
            {"name": "head", "args": {"1": "it"}},
            {"name": "it-noun", "args": {"1": "m", "2": "greci"}},
        ]
    }
    assert heads.head_template(entry) == ("it-noun", {"1": "m", "2": "greci"})
    assert heads.head_template({}) == ("", {})
