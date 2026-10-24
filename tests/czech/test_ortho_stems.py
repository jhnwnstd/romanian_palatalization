"""Stems of nouns whose root vowel alternates in length (data-layer fix).

Without length folding, the common prefix of dráhy ~ drahám stops at
the vowel (dr-) and the stem final, hence the target class, is wrong
(dráha r instead of velar, skála velar instead of lateral). Folding
must not run into the plural endings either (Frankám ~ Frankami). The
stem is read off GEN.SG (audit 2026-10-23); the folded common prefix
is its consistency check and the stem of nouns without GEN.SG.
"""

from __future__ import annotations

import pytest

from lexpipe.langs.czech import ortho

pytestmark = pytest.mark.cs


def forms(
    gen: str, dat_pl: str, ins_pl: str, ins: str
) -> dict[str, list[str]]:
    return {
        "GEN.SG": [gen],
        "DAT.PL": [dat_pl],
        "INS.PL": [ins_pl],
        "INS.SG": [ins],
    }


@pytest.mark.parametrize(
    ("lemma", "cells", "stem", "target"),
    [
        (
            "dráha",
            forms("dráhy", "drahám", "drahami", "drahou"),
            "dráh",
            "velar",
        ),
        (
            "skála",
            forms("skály", "skalám", "skalami", "skálou"),
            "skál",
            "lateral",
        ),
        (
            "kráva",
            forms("krávy", "kravám", "kravami", "krávou"),
            "kráv",
            "labial",
        ),
        (
            "Turek",
            forms("Turka", "Turkům", "Turky", "Turkem"),
            "Turk",
            "velar",
        ),
        ("dům", forms("domu", "domům", "domy", "domem"), "dom", "labial"),
        ("kost", forms("kosti", "kostem", "kostmi", "kostí"), "kost", "soft"),
    ],
)
def test_noun_stem(
    lemma: str, cells: dict[str, list[str]], stem: str, target: str
) -> None:
    found, method = ortho.noun_stem(cells, lemma)
    assert (found, method) == (stem, "genitive")
    assert ortho.common_stem(cells) == stem  # the consistency check agrees
    decl = ortho.declension(cells, found)
    assert (
        ortho.target_info(found, declension_type=decl).target_class == target
    )


def test_plural_endings_are_not_folded_into_the_stem() -> None:
    cells = {"DAT.PL": ["Frankám"], "INS.PL": ["Frankami"]}
    assert ortho.noun_stem(cells, "Franky") == ("Frank", "lcp")


def test_nadrazi_keeps_its_stem() -> None:
    cells = forms("nádraží", "nádražím", "nádražími", "nádražím")
    assert ortho.noun_stem(cells, "nádraží")[0] == "nádraž"


def test_align_folds_vowel_length() -> None:
    plain = ortho.align("dráh", "drahách")
    assert (plain.aligned, plain.ending, plain.orth_change) == (
        True,
        "ách",
        "",
    )
    mutated = ortho.align("dráh", "draze")
    assert (mutated.aligned, mutated.ending, mutated.orth_change) == (
        True,
        "e",
        "h→z",
    )


@pytest.mark.parametrize(
    ("lemma", "cells", "stem"),
    [
        (
            "chovankyně",
            {
                "DAT.SG": ["chovankyni"],
                "INS.SG": ["chovankyní"],
                "DAT.PL": ["chovankyním"],
                "INS.PL": ["chovankyněmi"],
            },
            "chovankyn",
        ),
        ("lázně", {"DAT.PL": ["lázním"], "INS.PL": ["lázněmi"]}, "lázn"),
        ("Polabí", {"INS.SG": ["Polabím"], "DAT.SG": ["Polabí"]}, "Polab"),
    ],
)
def test_soft_declension_without_genitive(
    lemma: str, cells: dict[str, list[str]], stem: str
) -> None:
    found, _ = ortho.noun_stem(cells, lemma)
    assert found == stem
    analysis = ortho.analyse_stem("noun", lemma, cells)
    assert analysis.declension == ortho.SOFT_DECL
    assert analysis.target.target_class == "soft"


def test_hard_plurale_tantum_stays_unknown() -> None:
    cells = {"DAT.PL": ["botám"], "INS.PL": ["botami"]}
    analysis = ortho.analyse_stem("noun", "boty", cells)
    assert analysis.declension == ortho.UNKNOWN
    assert analysis.target.target_class == "coronal_inferred"
