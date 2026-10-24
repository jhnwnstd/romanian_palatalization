"""Velar-stem verbs: infinitives, cells, form parsing, paradigm types."""

from __future__ import annotations

import pytest

from lexpipe.langs.italian import verbs

pytestmark = pytest.mark.it


@pytest.mark.parametrize(
    ("key", "expected"),
    [
        ("cercare", ("are", "cerc", False)),
        ("amicarsi", ("are", "amic", True)),
        ("vincere", ("ere", "vinc", False)),
        ("bianchire", ("ire", "bianch", False)),
        ("cerca", None),
    ],
)
def test_infinitive_parts(
    key: str, expected: tuple[str, str, bool] | None
) -> None:
    assert verbs.infinitive_parts(key) == expected


def test_velar_stem_and_target_letters() -> None:
    assert verbs.velar_stem("bianch") == "bianc"
    assert verbs.velar_stem("cerc") == "cerc"
    assert verbs.velar_stem("lanci") == ""  # palatal spelling
    assert verbs.velar_stem("macchi") == ""  # velar + glide
    assert verbs.target_letters("conosc") == "sc"
    assert verbs.target_letters("legg") == "gg"
    assert verbs.target_letters("vinc") == "c"


def test_cell_of_and_strip_clitic() -> None:
    assert (
        verbs.cell_of(["indicative", "present", "second-person", "singular"])
        == "2sg"
    )
    assert verbs.cell_of(["imperative", "second-person", "singular"]) == ""
    assert (
        verbs.cell_of(["present", "subjunctive", "third-person", "singular"])
        == "sbjv3sg"
    )
    assert (
        verbs.cell_of(
            [
                "proscribed",
                "present",
                "subjunctive",
                "third-person",
                "singular",
            ]
        )
        == ""
    )
    assert verbs.strip_clitic("ti amìchi") == "amìchi"
    assert verbs.strip_clitic("vìnci") == "vìnci"


@pytest.mark.parametrize(
    ("cell", "form", "vowel", "edge"),
    [
        ("2sg", "vinci", "i", verbs.PALATAL),
        ("3sg", "vince", "e", verbs.PALATAL),
        ("1pl", "vinciamo", "i", verbs.PALATAL),
        ("1sg", "vinco", "back", verbs.VELAR),
        ("sbjv3sg", "vinca", "back", verbs.VELAR),
        ("2sg", "cerchi", "i", verbs.VELAR),
        ("1pl", "cerchiamo", "i", verbs.VELAR),
        ("3sg", "cerca", "back", verbs.VELAR),
        ("1sg", "taccio", "back", verbs.PALATAL),  # over-application
        ("past1sg", "tacqui", "back", verbs.VELAR),  # velar before /w/
        ("2sg", "bianchisci", "i", verbs.VELAR),
        ("1sg", "agisco", "i", verbs.PALATAL),
        ("1pl", "taciamo", "i", verbs.PALATAL),
    ],
)
def test_parse_form(cell: str, form: str, vowel: str, edge: str) -> None:
    parsed = verbs.parse_form(cell, form)
    assert (parsed.vowel, parsed.edge) == (vowel, edge)


def _cells(*items: tuple[str, str]) -> list[tuple[str, verbs.Parsed]]:
    return [(cell, verbs.parse_form(cell, form)) for cell, form in items]


def test_paradigm_types() -> None:
    vincere = _cells(("1sg", "vinco"), ("2sg", "vinci"), ("3sg", "vince"))
    assert verbs.paradigm_type(vincere) == (verbs.ALTERNATING, False)
    cercare = _cells(("1sg", "cerco"), ("2sg", "cerchi"), ("3sg", "cerca"))
    assert verbs.paradigm_type(cercare) == (verbs.VELAR_THROUGHOUT, False)
    tacere = _cells(("1sg", "taccio"), ("2sg", "taci"), ("past1sg", "tacqui"))
    assert verbs.paradigm_type(tacere) == (verbs.OVERAPPLYING, True)
    cucire = _cells(("1sg", "cucio"), ("2sg", "cuci"), ("3sg", "cuce"))
    assert verbs.paradigm_type(cucire) == (verbs.PALATAL_THROUGHOUT, True)
    assert verbs.paradigm_type([]) == (verbs.NO_CELLS, False)


def _lemma(
    word: str, canonical: str, category: bool = False, ipa: int | None = None
) -> verbs.VerbLemma:
    parts = verbs.infinitive_parts(word)
    assert parts is not None
    ending, written, reflexive = parts
    return verbs.VerbLemma(
        verb_id=f"it:{word}:verb:0:00000000",
        word=word,
        verb_key=written + ending,
        reflexive=reflexive,
        canonical=canonical,
        ending=ending,
        stem_written=written,
        stem_key=verbs.velar_stem(written),
        root_stressed_category=category,
        stratum="inherited",
        ipa_stress=ipa,
    )


def test_conjugation_classes() -> None:
    vincere = _lemma("vincere", "vìncere")
    stress, source = verbs.infinitive_stress(vincere, {})
    assert (stress, source) == ("root", "canonical")
    assert verbs.conj_class(vincere, {}, stress) == verbs.ERE_ROOT
    tacere = _lemma("tacere", "tacére")
    stress, _ = verbs.infinitive_stress(tacere, {})
    assert verbs.conj_class(tacere, {}, stress) == verbs.ERE_STRESSED
    # the IPA decides first (/ˈvin.t͡ʃe.re/: three syllables from the end)
    from_ipa = _lemma("vincere", "vincére", ipa=3)
    assert verbs.infinitive_stress(from_ipa, {}) == ("root", "ipa")
    assert verbs.canonical_stress(from_ipa) == "ending"
    assert verbs.ipa_stress({"sounds": [{"ipa": "/taˈt͡ʃe.re/"}]}) == 2
    noted = _lemma("dispiacere", "dispiacére [auxiliary essere]")
    assert verbs.infinitive_stress(noted, {}) == ("ending", "canonical")
    bianchire = _lemma("bianchire", "bianchìre")
    forms = {"1sg": ["bianchìsco"]}
    assert verbs.conj_class(bianchire, forms, "") == verbs.ISC
    fuggire = _lemma("fuggire", "fuggìre")
    assert verbs.conj_class(fuggire, {"1sg": ["fùggo"]}, "") == verbs.IRE
    cercare = _lemma("cercare", "cercàre")
    assert verbs.conj_class(cercare, {}, "") == verbs.ARE
    assert verbs.dire_stem({"1sg": ["benedìco"]}) == "benedic"
