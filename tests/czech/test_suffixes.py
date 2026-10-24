"""Derivational suffixes: decomposition, candidates, paradigms (R14, R21)."""

from __future__ import annotations

import pytest

from lexpipe.langs.czech import suffixes as sx

pytestmark = pytest.mark.cs


@pytest.mark.parametrize(
    ("derivative", "stem", "lemma", "cluster", "suffix", "outcome", "change"),
    [
        ("knižní", "knih", "kniha", "", "-ní/-ný", "pal", "h→ž"),
        ("knížka", "knih", "kniha", "", "-ka", "pal", "h→ž"),
        ("ruční", "ruk", "ruka", "", "-ní/-ný", "pal", "k→č"),
        ("ručička", "ruk", "ruka", "", "-ka", "pal", "k→č"),
        ("ručit", "ruk", "ruka", "", "-it", "pal", "k→č"),
        ("ručník", "ruk", "ruka", "", "-ník", "pal", "k→č"),
        ("vzdušný", "vzduch", "vzduch", "", "-ní/-ný", "pal", "ch→š"),
        ("vzduchový", "vzduch", "vzduch", "", "-ový/-ovní", "plain", ""),
        ("bankovní", "bank", "banka", "", "-ový/-ovní", "plain", ""),
        ("sněžit", "sněh", "sníh", "", "-it", "pal", "h→ž"),
        ("klouček", "kluk", "kluk", "", "-ek", "pal", "k→č"),
        ("růžek", "roh", "roh", "", "-ek", "pal", "h→ž"),
        ("říční", "řek", "řeka", "", "-ní/-ný", "pal", "k→č"),
        ("český", "Čech", "Čech", "", "-ský/-cký", "pal", "ch→š"),
        ("otrocký", "otrok", "otrok", "", "-ský/-cký", "pal", "k→c"),
        ("pražský", "Prah", "Praha", "", "-ský/-cký", "pal", "h→ž"),
        ("turecký", "Turk", "Turek", "", "-ský/-cký", "pal", "k→c"),
        ("barokní", "barok", "baroko", "", "-ní/-ný", "exception", ""),
        ("vrchní", "vrch", "vrch", "", "-ní/-ný", "exception", ""),
        ("tankista", "tank", "tank", "", "-ista", "plain", ""),
        ("žákyně", "žák", "žák", "", "-yně", "plain", ""),
        ("žačka", "žák", "žák", "", "-ka", "pal", "k→č"),
        ("bioložka", "biolog", "biolog", "", "-ka", "pal", "g→ž"),
        ("bratřit", "bratr", "bratr", "Cr", "-it", "pal", "r→ř"),
        ("bratrský", "bratr", "bratr", "Cr", "-ský/-cký", "nontarget", ""),
    ],
)
def test_decompose(
    derivative: str,
    stem: str,
    lemma: str,
    cluster: str,
    suffix: str,
    outcome: str,
    change: str,
) -> None:
    found = sx.decompose(derivative, stem, lemma, cluster=cluster)
    assert found is not None
    target = "r" if stem.endswith("r") else "velar"
    assert found.suffix.label == suffix
    assert found.outcome(target, cluster) == outcome
    assert found.orth_change == change


def test_prefixes_only_when_allowed() -> None:
    assert sx.decompose("snožit", "noh", "noha") is None
    found = sx.decompose("snožit", "noh", "noha", allow_prefix=True)
    assert found is not None
    assert (found.prefix, found.suffix.label, found.mutated) == (
        "s",
        "-it",
        True,
    )


def test_unrelated_word_does_not_decompose() -> None:
    assert sx.decompose("knihovna", "knih", "kniha") is None
    assert sx.decompose("parkoviště", "park", "park") is None


def test_suffix_targets() -> None:
    ek = sx.BY_LABEL["-ek"]
    assert ek.palatalizes("r", "Cr")  # vnitro → vnitřek
    assert not ek.palatalizes("r", "")  # večer → večírek
    assert not sx.BY_LABEL["-ský/-cký"].palatalizes("r")
    assert not sx.BY_LABEL["-ový/-ovní"].palatalizes("velar")
    assert {s.pair for s in sx.SUFFIXES if s.pair} == set(sx.PAIRS)


def test_candidates_both_spellings() -> None:
    words = {
        c.word: c for c in sx.candidates("vzduch", "", sx.BY_LABEL["-ní/-ný"])
    }
    assert words["vzdušný"].mutated and words["vzdušný"].size == "big"
    assert not words["vzduchný"].mutated
    plain = [
        c.word for c in sx.candidates("vzduch", "", sx.BY_LABEL["-ový/-ovní"])
    ]
    assert "vzduchový" in plain


def test_candidates_spell_t_before_front_vowels() -> None:
    words = [c.word for c in sx.candidates("česk", "sk", sx.BY_LABEL["-ina"])]
    assert "čeština" in words
    assert "češťina" not in words


def test_diminutive_candidates_lengthen_the_root_vowel() -> None:
    words = [c.word for c in sx.candidates("knih", "", sx.BY_LABEL["-ka"])]
    assert "knížka" in words
    words = [c.word for c in sx.candidates("kluk", "", sx.BY_LABEL["-ek"])]
    assert "klouček" in words


def test_paradigm_avoids_forms_of_other_paradigms() -> None:
    forms = sx.paradigm("klukův", sx.BY_LABEL["-ův"])
    assert "klukovi" not in forms  # DAT.SG of kluk
    assert "klukových" not in forms  # a form of an -ový adjective
    assert sx.paradigm("knížka", sx.BY_LABEL["-ka"])[:3] == (
        "knížka",
        "knížky",
        "knížku",
    )
    assert "ručíme" in sx.paradigm("ručit", sx.BY_LABEL["-it"])


def test_fold() -> None:
    assert sx.fold("Klouček") == "kluček"
    assert sx.fold_ei("říční") == "řični"
    assert sx.final_unit("česk", "sk") == "sk"
    assert sx.final_unit("vzduch") == "ch"
