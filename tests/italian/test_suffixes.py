"""The suffix inventory and the stem-edge reading (design §5, R21)."""

from __future__ import annotations

import pytest

from lexpipe.langs.italian import suffixes as sx

pytestmark = pytest.mark.it

S = sx.BY_NAME


@pytest.mark.parametrize(
    ("stem", "suffix", "key", "edge"),
    [
        ("grec", "-ista", "grecista", sx.Edge.PALATAL),
        ("leg", "-ista", "leghista", sx.Edge.VELAR),
        ("fuoc", "-ista", "fochista", sx.Edge.VELAR),  # mobile diphthong
        ("ciec", "-ità", "cecità", sx.Edge.PALATAL),
        ("antic", "-ità", "antichità", sx.Edge.VELAR),
        ("monac", "-esimo", "monachesimo", sx.Edge.VELAR),
        ("chirurg", "-ia", "chirurgia", sx.Edge.PALATAL),
        ("banc", "-iere", "banchiere", sx.Edge.VELAR),
        ("barc", "-etto", "barchetta", sx.Edge.VELAR),  # gender variant
        ("fuoc", "-erello", "fuocherello", sx.Edge.VELAR),
        ("caric", "-are", "caricare", sx.Edge.BACK),
        ("arcaic", "-ismo", "arcaismo", sx.Edge.TRUNCATED),
        ("ceramic", "-ista", "ceramista", sx.Edge.TRUNCATED),
        ("acromatic", "-izzare", "acromatizzare", sx.Edge.TRUNCATED),
        ("barc", "-are", "barcollare", sx.Edge.MISMATCH),
    ],
)
def test_read_edge(stem: str, suffix: str, key: str, edge: sx.Edge) -> None:
    assert sx.read_edge(stem, S[suffix], key)[0] == edge


def test_prefixes_only_when_allowed() -> None:
    assert sx.read_edge("bosc", S["-ire"], "imboschire")[0] == sx.Edge.MISMATCH
    found = sx.read_edge("bosc", S["-ire"], "imboschire", allow_prefix=True)
    assert found == (sx.Edge.VELAR, "ire")


def test_candidates_cover_both_spellings() -> None:
    forms = [c.form for c in sx.candidates("leg", S["-ista"])]
    assert forms == ["legista", "leghista"]
    spelled = {c.form: c.spelling for c in sx.candidates("fuoc", S["-etto"])}
    assert spelled["fochetto"] == sx.Edge.VELAR
    assert spelled["focetto"] == sx.Edge.PALATAL
    assert [c.form for c in sx.candidates("caric", S["-are"])] == ["caricare"]
    assert sx.candidates("grec", S["-ia"]) == []  # -ia is not generated


def test_other_spelling_and_canonical() -> None:
    assert sx.other_spelling("leghista", sx.Edge.VELAR, "ista") == "legista"
    assert sx.other_spelling("grecista", sx.Edge.PALATAL, "ista") == (
        "grechista"
    )
    assert sx.other_spelling("caricare", sx.Edge.BACK, "are") == ""
    assert sx.canonical("ezza") is S["-ezza"]  # {{suffix}} drops the hyphen
    assert sx.canonical("-etta") is S["-etto"]
    assert sx.canonical("-iero") is S["-iere"]
    assert sx.canonical("-mento") is None


def test_best_reading_prefers_the_longest_body() -> None:
    found = sx.best_reading("fuoc", "fuocherello")
    assert found is not None and found[0] is S["-erello"]
    found = sx.best_reading("turc", "turchità")
    assert found is not None and found[:2] == (S["-ità"], sx.Edge.VELAR)
    assert sx.best_reading("bianc", "biancaneve") is None


def test_corpus_forms() -> None:
    assert sx.corpus_forms("leghista", S["-ista"]) == [
        "leghista",
        "leghisti",
        "leghiste",
    ]
    # verbs: infinitive and past participles only (the 3sg of an -are
    # verb is often the base noun itself: amica)
    assert sx.corpus_forms("tipicizzare", S["-izzare"]) == [
        "tipicizzare",
        "tipicizzato",
        "tipicizzata",
        "tipicizzati",
        "tipicizzate",
    ]
    assert sx.corpus_forms("antichità", S["-ità"]) == ["antichità"]
    assert sx.corpus_forms("amichevole", S["-evole"]) == [
        "amichevole",
        "amichevoli",
    ]


def test_groups_match_the_design() -> None:
    assert sx.LEARNED_SET == (
        "-issimo",
        "-izzare",
        "-ista",
        "-ismo",
        "-esimo",
        "-ità",
        "-izia",
        "-ia",
    )
    assert set(sx.NATIVE_SET) == {
        "-ezza",
        "-etto",
        "-ino",
        "-erello",
        "-ello",
        "-eria",
        "-essa",
        "-evole",
        "-eggiare",
        "-ire",
        "-ese",
        "-iere",
    }
    assert {s.name for s in sx.SUFFIXES if s.spellout} == {
        "-izzare",
        "-eggiare",
        "-ire",
        "-are",
    }
