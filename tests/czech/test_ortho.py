"""Stems, targets, mutations and the P/N table of design §6."""

from __future__ import annotations

import pytest

from lexpipe.langs.czech import ortho

pytestmark = pytest.mark.cs


def forms(**cells: str) -> dict[str, list[str]]:
    return {
        cell.replace("_", "."): value.split("|")
        for cell, value in cells.items()
    }


# -- stems ---------------------------------------------------------------------


def test_fleeting_e_and_quantity() -> None:
    turek = forms(
        GEN_SG="Turka", DAT_PL="Turkům", INS_PL="Turky", INS_SG="Turkem"
    )
    assert ortho.noun_stem(turek, "Turek") == ("Turk", "genitive")
    assert ortho.common_stem(turek) == "Turk"
    dum = forms(GEN_SG="domu", DAT_PL="domům", INS_PL="domy", INS_SG="domem")
    assert ortho.noun_stem(dum, "dům")[0] == "dom"


def test_shared_long_vowel_is_not_stem() -> None:
    nadrazi = forms(
        GEN_SG="nádraží",
        DAT_PL="nádražím",
        INS_PL="nádražími",
        INS_SG="nádražím",
    )
    stem, _ = ortho.noun_stem(nadrazi, "nádraží")
    assert stem == "nádraž"
    assert ortho.declension(nadrazi, stem) == ortho.SOFT_DECL


def test_stem_fallbacks() -> None:
    assert ortho.noun_stem(forms(GEN_SG="plechu"), "plech") == (
        "plech",
        "genitive",
    )
    assert ortho.noun_stem(forms(DAT_PL="plechům"), "plechy") == (
        "plech",
        "strip",
    )
    assert ortho.noun_stem(
        forms(DAT_PL="plechům", INS_PL="plechy"), "plechy"
    ) == ("plech", "lcp")
    assert ortho.noun_stem({}, "plech") == ("plech", "lemma")


@pytest.mark.parametrize(
    ("genitive", "stem", "expected"),
    [
        ("kluka", "kluk", ortho.HARD),
        ("muže", "muž", ortho.SOFT_DECL),
        ("kosti", "kost", ortho.ISTEM),
        ("vrátného", "vrátn", ortho.ADJECTIVAL),
    ],
)
def test_declension(genitive: str, stem: str, expected: str) -> None:
    assert ortho.declension({"GEN.SG": [genitive]}, stem) == expected


# -- targets -------------------------------------------------------------------


@pytest.mark.parametrize(
    ("stem", "final", "cluster", "target"),
    [
        ("kluk", "k", "", ortho.VELAR),
        ("Čech", "ch", "", ortho.VELAR),
        ("česk", "k", "sk", ortho.VELAR),
        ("německ", "k", "ck", ortho.VELAR),
        ("bratr", "r", "Cr", ortho.R_CLASS),
        ("doktor", "r", "", ortho.R_CLASS),
        ("student", "t", "", ortho.CORONAL),
        ("most", "t", "st", ortho.CORONAL),
        ("hvězd", "d", "zd", ortho.CORONAL),
        ("otc", "c", "", ortho.C_CLASS),
        ("les", "s", "", ortho.SIBILANT),
        ("ryb", "b", "", ortho.LABIAL),
        ("stol", "l", "", ortho.LATERAL),
        ("muž", "ž", "", ortho.SOFT),
        ("rádi", "i", "", ortho.VOWEL),
    ],
)
def test_target_classes(
    stem: str, final: str, cluster: str, target: str
) -> None:
    info = ortho.target_info(stem)
    assert (info.stem_final, info.cluster, info.target_class) == (
        final,
        cluster,
        target,
    )


def test_soft_declension_and_hints_override() -> None:
    assert ortho.target_info(
        "kost", declension_type=ortho.ISTEM
    ).target_class == (ortho.SOFT)
    assert (
        ortho.target_info("absinth", flags=["-velar"]).target_class
        == ortho.OTHER
    )
    assert (
        ortho.target_info("kluk").core
        and not ortho.target_info("student").core
    )


# -- alignment -----------------------------------------------------------------


@pytest.mark.parametrize(
    ("stem", "form", "ending", "change"),
    [
        ("kluk", "kluci", "i", "k→c"),
        ("Čech", "Češi", "i", "ch→š"),
        ("česk", "čeští", "í", "sk→šť"),
        ("německ", "němečtí", "í", "ck→čť"),
        ("Etrusk", "Etruscích", "ích", "sk→sc"),
        ("vrah", "vrazi", "i", "h→z"),
        ("biolog", "biolozích", "ích", "g→z"),
        ("doktor", "doktoři", "i", "r→ř"),
        ("boh", "bože", "e", "h→ž"),
        ("otc", "otče", "e", "c→č"),
        ("potok", "potoku", "u", ""),
        ("student", "studenti", "i", ""),
    ],
)
def test_align(stem: str, form: str, ending: str, change: str) -> None:
    alignment = ortho.align(stem, form)
    assert alignment.aligned
    assert (alignment.ending, alignment.orth_change) == (ending, change)


def test_t_spelling_only_before_front_vowels() -> None:
    # "št" stands for šť only before i/í/ě; before ou it is not sk→šť.
    assert ortho.align("pusk", "puští").orth_change == "sk→šť"
    assert not ortho.align("pusk", "puštou").aligned
    assert not ortho.align("kluk", "dům").aligned


# -- P/N table -----------------------------------------------------------------

VELAR = ortho.target_info("kluk")
CR = ortho.target_info("bratr")
R_AFTER_VOWEL = ortho.target_info("doktor")
CORONAL = ortho.target_info("student")
EC = ortho.TargetInfo("c", "", ortho.C_CLASS, ec_type=True)
CE = ortho.target_info("soudc")


@pytest.mark.parametrize(
    ("cell", "gender", "ending", "target", "expected"),
    [
        ("NOM.PL", "m.an", "i", VELAR, (ortho.P, "")),
        ("NOM.PL", "m.an", "ové", VELAR, (ortho.N, ortho.N_ALLOMORPH)),
        ("NOM.PL", "m.an", "é", CORONAL, (ortho.N, ortho.N_ALLOMORPH)),
        ("NOM.PL", "m.in", "y", VELAR, (ortho.NONE, "")),
        ("LOC.SG", "m.in", "u", VELAR, (ortho.N, ortho.N_ALLOMORPH)),
        ("LOC.SG", "m.in", "e", VELAR, (ortho.P, "")),
        ("LOC.SG", "m.an", "u", VELAR, (ortho.NONE, "")),
        ("DAT.SG", "f", "e", VELAR, (ortho.P, "")),
        ("LOC.PL", "n", "ách", VELAR, (ortho.N, ortho.N_ALLOMORPH)),
        ("LOC.PL", "m.an", "ích", VELAR, (ortho.P, "")),
        ("LOC.PL", "m.an", "ech", CORONAL, (ortho.N, ortho.N_NONTRIGGER)),
        ("LOC.PL", "f", "ách", VELAR, (ortho.NONE, "")),
        ("VOC.SG", "m.an", "e", VELAR, (ortho.P, "")),
        ("VOC.SG", "m.an", "u", VELAR, (ortho.NONE, "")),
        ("VOC.SG", "m.an", "e", CR, (ortho.P, "")),
        ("VOC.SG", "m.an", "e", R_AFTER_VOWEL, (ortho.N, ortho.N_NONTRIGGER)),
        ("VOC.SG", "m.an", "e", CORONAL, (ortho.N, ortho.N_NONTRIGGER)),
        ("VOC.SG", "m.an", "e", EC, (ortho.P, "")),
        ("VOC.SG", "m.an", "e", CE, (ortho.N, ortho.N_NONTRIGGER)),
        ("INS.SG", "m.an", "em", VELAR, (ortho.N, ortho.N_NONTRIGGER)),
        ("NOM.PL.M.AN", "adj", "í", VELAR, (ortho.P, "")),
    ],
)
def test_classify(
    cell: str,
    gender: str,
    ending: str,
    target: ortho.TargetInfo,
    expected: tuple[str, str],
) -> None:
    verdict = ortho.classify(cell, gender, ending, target)
    assert (verdict.allomorph_class, verdict.n_kind) == expected


def test_mutation_kinds() -> None:
    studenti = ortho.mutation_of(
        ortho.align("student", "studenti"), CORONAL, ortho.P
    )
    assert (studenti.mutation, studenti.kind, studenti.palatal_output) == (
        True,
        ortho.MUTATION_INFERRED,
        "c",
    )
    rybe = ortho.mutation_of(
        ortho.align("ryb", "rybě"), ortho.target_info("ryb"), ortho.P
    )
    assert (rybe.mutation, rybe.kind, rybe.palatal_output) == (
        False,
        ortho.MUTATION_GLIDE,
        "bj",
    )
    kluci = ortho.mutation_of(ortho.align("kluk", "kluci"), VELAR, ortho.P)
    assert kluci.palatal_output == "t͡s"


def test_exception_labels() -> None:
    p = ortho.Classification(ortho.P)
    n = ortho.Classification(ortho.N, ortho.N_ALLOMORPH)
    nontrigger = ortho.Classification(ortho.N, ortho.N_NONTRIGGER)

    def label(**kw: object) -> str:
        base: dict[str, object] = {
            "in_domain": True,
            "cell": "NOM.PL",
            "classification": p,
            "target": VELAR,
            "mutation": True,
        }
        base.update(kw)
        return ortho.exception_label(**base)  # type: ignore[arg-type]

    assert label() == ortho.UNDERGOER
    assert label(mutation=False) == ortho.UNEXPLAINED
    assert label(classification=n) == ortho.ALLOMORPH_LABEL
    assert label(classification=nontrigger) == ortho.NON_EXCEPTION
    assert label(target=ortho.target_info("muž")) == ortho.NON_EXCEPTION
    assert label(in_domain=False) == ortho.NON_EXCEPTION
    vowel = ortho.target_info("energi")
    assert label(target=vowel, nde_class="tautomorphemic_dorsal") == (
        "nde:tautomorphemic_dorsal"
    )


# -- suffixes and NDE ------------------------------------------------------------


def test_suffix_targets() -> None:
    biolog = ortho.target_info("biolog")
    templates = [
        {"name": "ety", "args": {"2": ":af", "3": "bio-", "4": "-log"}}
    ]
    assert ortho.lemma_suffix("biolog") == "log"
    assert ortho.suffix_evidence("log", templates, []) == "etymology"
    assert ortho.target_is_suffix("log", "etymology", biolog)
    # -ka without morphological evidence is not a suffix (ruka, louka)
    assert ortho.suffix_evidence(ortho.lemma_suffix("ruka"), [], []) == ""
    assert (
        ortho.suffix_evidence("ek", [], ["Czech terms suffixed with -ek"])
        == "category"
    )
    # -graf contains no target
    assert not ortho.target_is_suffix(
        "graf", "orthographic", ortho.target_info("fotograf")
    )


def test_nde_classes() -> None:
    energie = ortho.target_info("energi")
    assert ortho.nde_class("energi", energie, []) == "tautomorphemic_dorsal"
    radio = ortho.target_info("rádi")
    assert (
        ortho.nde_class("rádi", radio, ["[ˈraːdɪjo]"]) == "loan_hard_coronal"
    )
    assert ortho.nde_class("kluk", VELAR, []) == ""
    assert ortho.ndeb_tokens("kino", ["[ˈkɪno]"]) == ("ki",)
    assert ortho.ndeb_tokens("diplomat", ["[ˈdɪplomat]"]) == ("dɪ",)


def test_opportunity_cells() -> None:
    assert ortho.opportunity_cells("noun", ["m.an"]) == {
        "NOM.PL",
        "VOC.SG",
        "LOC.PL",
    }
    assert ortho.opportunity_cells("noun", ["f"]) == {"DAT.SG", "LOC.SG"}
    assert ortho.opportunity_cells("adj", []) == {"NOM.PL.M.AN"}


def test_analyse_stem_marks_ec_nouns() -> None:
    otec = forms(GEN_SG="otce", DAT_PL="otcům", INS_PL="otci", INS_SG="otcem")
    analysis = ortho.analyse_stem("noun", "otec", otec)
    assert analysis.target.ec_type and analysis.stem == "otc"
    soudce = forms(
        GEN_SG="soudce", DAT_PL="soudcům", INS_PL="soudci", INS_SG="soudcem"
    )
    assert not ortho.analyse_stem("noun", "soudce", soudce).target.ec_type
    assert ortho.analyse_stem("adj", "český", {}).target.cluster == "sk"
    assert (
        ortho.analyse_stem("adj", "jarní", {}).target.target_class
        == ortho.SOFT
    )
