"""The port of Module:it-common (revision 79032535) and its defaults."""

from __future__ import annotations

import hashlib

import pytest

from lexpipe.langs.italian import itcommon
from tests.italian import helpers

pytestmark = pytest.mark.it

UPSTREAM = helpers.PROJECT_ROOT / "data" / "italian" / "upstream"


@pytest.mark.parametrize(
    ("term", "gender", "plural"),
    [
        ("bacio", "m", "baci"),  # -io → -i
        ("teologo", "m", "teologi"),  # -ologo → -ologi
        ("amico", "m", "amici"),  # -ico → -ici
        ("maniaco", "m", "maniaci"),  # -aco → -aci
        ("medico", "m", "medici"),
        ("lago", "m", "laghi"),  # other -go → -ghi
        ("parco", "m", "parchi"),
        ("catalogo", "m", "cataloghi"),  # -logo but not -ologo
        ("fuoco", "m", "fuochi"),
        ("gatto", "m", "gatti"),
        ("amica", "f", "amiche"),  # -ca f → -che
        ("patriarca", "m", "patriarchi"),  # -ca m → -chi
        ("collega", "m", "colleghi"),
        ("biologia", "f", "biologie"),  # -logia f
        ("valigia", "f", "valige"),  # -gia f → -ge (the explicit valigie wins)
        ("camicia", "m", "camici"),
        ("casa", "f", "case"),
        ("poeta", "m", "poeti"),
        ("luce", "f", "luci"),
        ("psiche", "f", "psichi"),
        ("città", "f", None),  # final accent: no default
        ("sport", "m", None),
    ],
)
def test_make_plural_rules(term: str, gender: str, plural: str | None) -> None:
    assert itcommon.make_plural(term, gender) == plural


def test_gender_specific_defaults_raise_for_mf() -> None:
    with pytest.raises(itcommon.DefaultPluralError):
        itcommon.make_plural("collega", "mfbysense")
    with pytest.raises(itcommon.DefaultPluralError):
        itcommon.make_plural("taggiasca", "?")
    with pytest.raises(itcommon.DefaultPluralError):
        itcommon.make_plural("biologia", "m")
    # The -e and -o rules are not gender-specific.
    assert itcommon.make_plural("cantante", "mfbysense") == "cantanti"


def test_capo_specials() -> None:
    assert itcommon.make_plural("capocuoco", "m", "cap*+") == "capicuochi"
    assert itcommon.make_plural("capotribù", "m", "cap*") == "capitribù"
    assert itcommon.make_plural("capoarea", "f", "cap*") == "capoarea"
    with pytest.raises(itcommon.DefaultPluralError):
        itcommon.make_plural("bosco", "m", "cap*")


def test_multiword_is_not_ported() -> None:
    with pytest.raises(itcommon.NotPortedError):
        itcommon.make_plural("abito da sposa", "m")
    with pytest.raises(itcommon.NotPortedError):
        itcommon.make_plural("divano-letto", "m")
    with pytest.raises(itcommon.NotPortedError):
        itcommon.make_plural("amico", "m", "first")


def test_feminine_and_masculine() -> None:
    assert itcommon.make_feminine("amico") == "amica"
    assert itcommon.make_feminine("editore") == "editrice"
    assert itcommon.make_feminine("birbone") == "birbona"
    assert itcommon.make_feminine("cantante") == "cantante"
    assert itcommon.make_masculine("amica") == "amico"
    assert itcommon.make_masculine("attrice") == "attore"
    # Like the module, multiword masculines fall through to the
    # single-word rules on the whole term.
    assert itcommon.make_masculine("buona amica") == "buona amico"


def test_noun_default_cells() -> None:
    assert itcommon.noun_default_cells("lago", "m") == [("pl", "laghi")]
    assert itcommon.noun_default_cells("collega", "mfbysense") == [
        ("m.pl", "colleghi"),
        ("f.pl", "colleghe"),
    ]
    assert itcommon.noun_default_cells("cantante", "mfbysense") == [
        ("pl", "cantanti")
    ]
    assert itcommon.noun_default_cells("sport", "m") == []


def test_adjective_default_cells() -> None:
    assert itcommon.adjective_default_cells("medico") == [
        ("m.pl", "medici"),
        ("f.pl", "mediche"),
    ]
    assert itcommon.adjective_default_cells("mobile") == [("pl", "mobili")]
    assert itcommon.adjective_default_cells("comunista") == [
        ("m.pl", "comunisti"),
        ("f.pl", "comuniste"),
    ]
    assert itcommon.adjective_default_cells("editore") == [
        ("m.pl", "editori"),
        ("f.pl", "editrici"),
    ]
    assert itcommon.adjective_default_cells(
        "invasore", ["invaditrice", "invasora"]
    ) == [("m.pl", "invasori"), ("f.pl", "invaditrici"), ("f.pl", "invasore")]
    assert itcommon.adjective_default_cells("aurica", feminine_only=True) == [
        ("f.pl", "auriche")
    ]
    assert itcommon.adjective_default_cells("blu") == []


def test_revision_constants() -> None:
    assert itcommon.MODULE_REVISION == 79032535
    assert str(itcommon.MODULE_REVISION) in itcommon.MODULE_URL
    assert str(itcommon.MODULE_REVISION) in (itcommon.__doc__ or "")
    assert itcommon.MODULE_SHA1 in (itcommon.__doc__ or "")


@pytest.mark.needs_data
def test_ported_revision_matches_the_fetched_source() -> None:
    """The fetched module text is the revision the port cites."""
    source = UPSTREAM / "Module_it-common.lua"
    record = UPSTREAM / "Module_it-common.revision.json"
    if not source.is_file() or not record.is_file():
        pytest.skip("fetched Module:it-common copy not present")
    data = source.read_bytes()
    digest = hashlib.sha1(data, usedforsecurity=False).hexdigest()
    assert digest == itcommon.MODULE_SHA1
    assert len(data) == itcommon.MODULE_SIZE
    assert str(itcommon.MODULE_REVISION) in record.read_text("utf-8")
    assert "ex.make_plural" in data.decode("utf-8")
