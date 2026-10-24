"""The shared etymology parser and R22 stratum classifier (lexpipe.etymology).

Italian and Czech read Wiktionary etymology templates with this one
module; the plugins only supply a profile. The cases below cover the
template tables (the union of the two former copies), the parsing
details, every branch of the classifier and the reconciled behaviours
listed in the module docstring.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lexpipe import etymology as et  # noqa: E402

IT = et.ITALIAN
CS = et.CZECH


def t(name: str, *args: str, **named: str) -> dict[str, Any]:
    values = {str(i): a for i, a in enumerate(args, 1)}
    values.update(named)
    return {"name": name, "args": values}


def kinds(templates: list[dict[str, Any]], code: str) -> list[tuple[Any, ...]]:
    return [(s.kind, s.lang, s.terms) for s in et.steps(templates, code)]


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


def test_classic_templates_and_variants() -> None:
    assert kinds([t("inh", "it", "la", "amīcus")], "it") == [
        ("inh", "la", ("amīcus",))
    ]
    assert kinds([t("inh-lite", "it", "la", "x")], "it")[0][0] == "inh"
    assert kinds([t("bor+", "cs", "it", "banca")], "cs")[0][0] == "bor"
    assert kinds([t("cal", "it", "en", "x")], "it")[0][0] == "calque"
    assert kinds([t("cog", "it", "ro", "amic")], "it") == []


def test_affix_parts_keep_their_prefix_and_record_languages() -> None:
    (step,) = et.steps([t("af", "cs", "la:chronicus", "-ický")], "cs")
    assert step.terms == ("la:chronicus", "-ický")  # as written
    assert step.term_langs == ("la", "")
    (step,) = et.steps(
        [t("af", "it", "a-", "χρῶμα", "-ia", lang2="grc")], "it"
    )
    assert step.term_lang(1) == "grc" and step.term_lang(0) == ""
    assert et.steps([t("con", "it", "x-", "-y")], "it")[0].kind == "af"
    assert et.steps([t("surface analysis", "it", "x", "-y")], "it")[
        0
    ].kind == ("surf")


def test_ety_template_steps() -> None:
    found = et.steps(
        [
            t(
                "ety",
                "cs",
                ":inh",
                "sla-pro:*strašiti",
                ":af<surf>",
                "strach",
                "-it",
            )
        ],
        "cs",
    )
    assert [(s.kind, s.lang) for s in found] == [
        ("inh", "sla-pro"),
        ("surf", "cs"),
    ]
    assert found[0].term_langs == ("sla-pro",)
    af = et.steps([t("ety", "it", ":af", "banca<t:bank>", "-ario")], "it")
    assert [(s.kind, s.lang, s.terms) for s in af] == [
        ("af", "it", ("banca", "-ario"))
    ]
    aliases = et.steps([t("ety", "cs", ":afeq", "x", ":ellipsis", "y")], "cs")
    assert [s.kind for s in aliases] == ["af", "clipping"]


def test_artifacts_parentheses_and_duplicates() -> None:
    clean = t("ety", "cs", ":inh", "zlw-ocs:kniha")
    artifact = t("ety", "cs", ":inh", "['zlw-ocs:', <LINK(['kniha']){} >]")
    found = et.steps([clean, artifact, clean], "cs")
    assert [(s.kind, s.terms) for s in found] == [("inh", ("kniha",))]
    (step,) = et.steps(
        [t("univerbation", "it", "(dipingere)", "a", "fresco")], "it"
    )
    assert step.terms == ("a", "fresco")
    assert et.clean_term("mieti(trice)<t:x>") == "mieti"
    assert et.clean_term("(nux) Abellāna") == "Abellāna"
    assert et.clean_term("indicum (folium)") == "indicum"
    # an optional letter inside a word is kept
    assert et.clean_term("*ga(v)ota") == "*ga(v)ota"
    assert et.clean_term("sla-pro:*o(b)krajь") == "sla-pro:*o(b)krajь"
    assert et.clean_term("*h₂(e)rmos") == "*h₂(e)rmos"
    assert et.steps(["not a template", 3, None], "it") == []


def test_deverbal_templates_of_each_language() -> None:
    assert et.steps([t("it-deverbal", "pesare")], "it")[0].terms == ("pesare",)
    assert et.steps([t("it-deverbal fpp", "x")], "it")[0].kind == "deverbal"
    assert et.steps([t("deverbal", "cs", "zrodit")], "cs")[0].terms == (
        "zrodit",
    )
    # Another language's own template is not this language's deverbal.
    assert et.steps([t("cs-deverbal", "x")], "it") == []
    assert et.steps([t("onom", "it")], "it")[0].kind == "onomatopoeic"


# ---------------------------------------------------------------------------
# Stratum
# ---------------------------------------------------------------------------


def stratum(
    templates: list[dict[str, Any]],
    text: str = "",
    profile: et.EtymologyProfile = IT,
    word: str = "",
) -> tuple[str, str]:
    return et.stratum(
        et.steps(templates, profile.code), text, profile, word=word
    )


@pytest.mark.parametrize(
    ("templates", "text", "profile", "expected"),
    [
        ([t("inh", "it", "la", "amīcus")], "", IT, (et.INHERITED, "inh")),
        ([t("bor", "it", "fr", "x")], "", IT, (et.BORROWED, "bor")),
        ([t("der", "it", "la", "x")], "", IT, (et.BORROWED, "der:la")),
        (
            [t("uder", "it", "la-vul", "*x")],
            "",
            IT,
            (et.INHERITED, "uder:la-vul"),
        ),
        (
            [t("der", "cs", "sla-pro", "*x")],
            "",
            CS,
            (et.INHERITED, "der:sla-pro"),
        ),
        ([t("der", "cs", "cs", "x")], "", CS, (et.INTERNAL, "der")),
        # :from is classified by its language (not always internal)
        ([t("ety", "cs", ":from", "la:x")], "", CS, (et.BORROWED, "from:la")),
        ([t("ety", "it", ":from", "caricare")], "", IT, (et.INTERNAL, "from")),
        (
            [t("af", "it", "sarco-", "-fago")],
            "From sarco- + -fago.",
            IT,
            (et.INTERNAL, "af"),
        ),
        ([t("surf", "it", "x", "-ico")], "", IT, (et.BORROWED, "surf")),
        ([], "", IT, (et.UNKNOWN, "")),
    ],
)
def test_template_strata(
    templates: list[dict[str, Any]],
    text: str,
    profile: et.EtymologyProfile,
    expected: tuple[str, str],
) -> None:
    assert stratum(templates, text, profile) == expected


def test_internal_kinds_whose_text_names_a_borrowing() -> None:
    # Italian rule, applied to Czech too: the text opens with a borrowing.
    assert stratum(
        [t("af", "cs", "la:āctīvus", "-ní")],
        "From Latin āctīvus + Czech -ní.",
        CS,
    ) == (et.BORROWED, "internal_analysis")
    assert stratum(
        [t("af", "it", "sporco", "-izia")],
        "From Latin spurcitia. Internally, from sporco + -izia.",
    ) == (et.BORROWED, "internal_analysis")


def test_surface_analysis_keeps_internal_bases() -> None:
    # A leading affixation of native bases stays internal (antichità,
    # altoatesino, elektromobilita); bases of another language do not.
    assert stratum(
        [t("af", "it", "antico", "-ità"), t("surf", "+bor+", "it", "la", "x")],
        "From antico + -ità. By surface analysis, borrowed from Latin x.",
    ) == (et.INTERNAL, "af")
    assert stratum(
        [t("compound", "it", "alto", "atesino")],
        "From alto + atesino. By surface analysis, Alto Adige + -ino.",
    ) == (et.INTERNAL, "af")
    assert stratum(
        [t("ety", "cs", ":af", "elektro-", "mobilita")],
        "From elektro- + mobilita. By surface analysis, elektromobil + -ita.",
        CS,
    ) == (et.INTERNAL, "af")
    # unikat: u- + Proto-Slavic *nikàti, an ancestor's base
    assert stratum(
        [t("ety", "cs", ":af", "u-", "sla-pro:*nikàti")],
        "From u- + Proto-Slavic *nikàti. By surface analysis, uniknout + -at.",
        CS,
    ) == (et.INTERNAL, "af")
    assert stratum(
        [t("af", "it", "ἀγράμματος", "-ismo", lang1="grc")],
        "From Ancient Greek ἀγράμματος + -ismo. By surface analysis, a- + -ismo.",
    ) == (et.BORROWED, "internal_analysis")
    assert stratum(
        [t("af", "it", "a-", "χρῶμα", "-ia", lang2="grc")],
        "From a- + Ancient Greek χρῶμα + -ia. By surface analysis, a- + -cromia.",
    ) == (et.BORROWED, "surface_analysis")
    # A text that opens with the surface analysis is a surface analysis.
    assert stratum(
        [t("ety", "it", ":af", "ampere", "-metrico")],
        "By surface analysis, ampere + -o- + -metrico.",
    ) == (et.BORROWED, "surface_analysis")


@pytest.mark.parametrize(
    ("text", "profile", "expected"),
    [
        ("Inherited from Old Czech prach.", CS, et.INHERITED),
        ("From Proto-Slavic *x.", CS, et.INHERITED),
        ("From Vulgar Latin *x.", IT, et.INHERITED),
        ("From Old Italian x.", IT, et.INHERITED),
        ("Borrowed from German Liga.", CS, et.BORROWED),
        ("Unadapted borrowing from Piedmontese x.", IT, et.BORROWED),
        ("From Latin liga.", CS, et.BORROWED),
        ("From sníh (“snow”, noun).", CS, et.INTERNAL),
        ("From acciacco.", IT, et.INTERNAL),
        ("From alto + -issimo. By surface analysis, x.", IT, et.INTERNAL),
        ("Compare x. By surface analysis, y + -z.", IT, et.BORROWED),
        ("Of uncertain origin.", IT, et.UNKNOWN),
    ],
)
def test_text_fallback(
    text: str, profile: et.EtymologyProfile, expected: str
) -> None:
    assert stratum([], text, profile)[0] == expected


def test_etymology_trees_are_skipped() -> None:
    tree = (
        "Etymology tree\nProto-Slavic *-ina\nOld Czech -ina\nCzech -ina\n"
        "Czech cenina\nFrom cenný + -ina. By surface analysis, cena + ina."
    )
    assert et.prose(tree, "cenina", CS) == (
        "From cenný + -ina. By surface analysis, cena + ina."
    )
    assert et.prose(tree) == (
        "From cenný + -ina. By surface analysis, cena + ina."
    )
    assert et.prose("From x.", "x", CS) == "From x."
    assert stratum([], tree, CS, "cenina")[0] == et.INTERNAL


def test_summary() -> None:
    found = et.steps([t("inh", "it", "la", "amīcus")], "it")
    assert et.summary(found, "it") == "inh:la:amīcus"
    found = et.steps([t("ety", "it", ":af", "banca", "-ario")], "it")
    assert et.summary(found, "it") == "af"
    found = et.steps([t("ety", "cs", ":from", "la:x")], "cs")
    assert et.summary(found, "cs") == "from:la:x"
    assert et.summary([], "it") == ""


def test_etymon_and_families() -> None:
    found = et.steps([t("ety", "it", ":bor", "fro:parc")], "it")
    assert et.etymon(found, "it") == ("fro", "parc")
    assert (
        et.etymon(et.steps([t("ety", "it", ":from", "x")], "it"), "it") is None
    )
    assert et.language_family("la-med") == "la"
    assert et.language_family("grc-koi") == "grc"
    assert et.fold("Àntico") == "antico"
