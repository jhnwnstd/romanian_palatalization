"""Etymology steps, stratum (R22), deverbal flags, -ico suffix (R26),
-logo/-fago/-urgo person/thing (R27)."""

from __future__ import annotations

from typing import Any

import pytest

from lexpipe import suffixoids
from lexpipe.langs.italian import etymology

pytestmark = pytest.mark.it


def _t(name: str, **args: str) -> dict[str, Any]:
    return {"name": name, "args": args}


def _steps(*templates: dict[str, Any]) -> list[etymology.Step]:
    return etymology.steps(templates)


def test_classic_and_ety_templates() -> None:
    found = _steps(
        _t("ety", **{"1": "it", "2": ":lbor", "3": "la:antīcus<t:anterior>"}),
        _t("lbor", **{"1": "it", "2": "la", "3": "antīcus"}),
        _t("ety", **{"1": "it", "2": ":lbor", "3": "la:antīcus<t:anterior>"}),
    )
    assert [(s.kind, s.lang, s.terms) for s in found] == [
        ("lbor", "la", ("antīcus",))
    ]
    af = _steps(
        _t("ety", **{"1": "it", "2": ":af", "3": "banca", "4": "-ario"})
    )
    assert [(s.kind, s.lang, s.terms) for s in af] == [
        ("af", "it", ("banca", "-ario"))
    ]
    deverbal = _steps(
        _t(
            "ety",
            **{
                "1": "it",
                "2": ":deverbal",
                "3": "concertare",
                "4": ":af",
                "5": "-o<id:deverbal>",
            },
        )
    )
    assert [s.kind for s in deverbal] == ["deverbal", "af"]
    assert _steps(_t("it-deverbal", **{"1": "pesare"}))[0].terms == ("pesare",)
    assert _steps(_t("cog", **{"1": "ro", "2": "amic"})) == []


def test_etymon_key_keeps_source_words_apart() -> None:
    # The key keeps the etymon's ending: Latin arca 'chest' and arcus
    # 'bow' (arca, arco) are different words; gender pairs (amīca,
    # amīcus) share a root through the feminine's stated counterpart.
    arca = _steps(_t("inh", **{"1": "it", "2": "la", "3": "arca"}))
    arcus = _steps(_t("inh", **{"1": "it", "2": "la", "3": "arcus"}))
    assert etymology.etymon_key(arca) == "la:arca"
    assert etymology.etymon_key(arca) != etymology.etymon_key(arcus)
    amico = _steps(_t("inh", **{"1": "it", "2": "la", "3": "amīcus"}))
    assert etymology.etymon_key(amico) == "la:amicus"
    antico = _steps(_t("inh", **{"1": "it", "2": "la", "3": "antīquus"}))
    anterior = _steps(_t("lbor", **{"1": "it", "2": "la", "3": "antīcus"}))
    assert etymology.etymon_key(antico) != etymology.etymon_key(anterior)
    park = _steps(_t("ety", **{"1": "it", "2": ":bor", "3": "fro:parc"}))
    frugal = _steps(_t("lbor", **{"1": "it", "2": "la", "3": "parcus"}))
    assert etymology.etymon_key(park) == "fro:parc"
    assert etymology.etymon_key(frugal) == "la:parcus"


@pytest.mark.parametrize(
    ("templates", "text", "expected"),
    [
        (
            [_t("inh", **{"1": "it", "2": "la", "3": "amīcus"})],
            "",
            "inherited",
        ),
        (
            [_t("bor+", **{"1": "it", "2": "la", "3": "pūblicus"})],
            "",
            "borrowed whole",
        ),
        (
            [_t("der", **{"1": "it", "2": "la", "3": "ēlectricus"})],
            "",
            "borrowed whole",
        ),
        (
            [_t("uder", **{"1": "it", "2": "la-vul", "3": "*trānsfrīcō"})],
            "",
            "inherited",
        ),
        (
            [_t("ety", **{"1": "it", "2": ":from", "3": "caricare"})],
            "",
            "internal formation",
        ),
        (
            [_t("af", **{"1": "it", "2": "sarco-", "3": "-fago"})],
            "From sarco- + -fago.",
            "internal formation",
        ),
        (
            [_t("af", **{"1": "it", "2": "sporco", "3": "-izia"})],
            "From Latin spurcitia. Internally, from sporco + -izia.",
            "borrowed whole",
        ),
        (
            [
                _t("der", **{"1": "it", "2": "la", "3": "theologus"}),
                _t("surf", **{"1": "it", "2": "teo-", "3": "-logo"}),
            ],
            "From Latin theologus. By surface analysis, teo- + -logo.",
            "borrowed whole",
        ),
        (
            # antichità: the affixation stands; the surface phrase names
            # an alternative borrowing, not an analysis of one
            [
                _t("af", **{"1": "it", "2": "antico", "3": "-ità"}),
                _t(
                    "surf",
                    **{"1": "+bor+", "2": "it", "3": "la", "4": "antiquitas"},
                ),
            ],
            "From antico + -ità. By surface analysis, borrowed from Latin "
            "antīquitātem.",
            "internal formation",
        ),
        (
            [
                _t("bor", **{"1": "it", "2": "fr", "3": "calorifère"}),
                _t("surf", **{"1": "it", "2": "calori-", "3": "-fero"}),
            ],
            "Borrowed from French calorifère. By surface analysis, from "
            "calori- + -fero.",
            "borrowed whole",
        ),
        ([], "", "unknown"),
    ],
)
def test_stratum(
    templates: list[dict[str, Any]], text: str, expected: str
) -> None:
    assert etymology.stratum(_steps(*templates), text)[0] == expected


def test_deverbal() -> None:
    carico = _steps(_t("ety", **{"1": "it", "2": ":from", "3": "caricare"}))
    text = "From the short past participle of caricare (“to load”) in Tuscan"
    assert etymology.deverbal_verb(carico, text) == "caricare"
    assert etymology.is_first_conjugation("caricare")
    assert etymology.deverbal_verb([], "Deverbal from pescare + -a.") == (
        "pescare"
    )
    assert not etymology.is_first_conjugation("salire")
    assert etymology.deverbal_verb([], "From Latin medicus.") == ""


def _ico(key: str, *templates: dict[str, Any], **kwargs: Any) -> Any:
    defaults: dict[str, Any] = {
        "ico_stressed": False,
        "adjectival": True,
        "deverbal": False,
        "italian_lemmas": frozenset({"elettro", "classe"}),
    }
    defaults.update(kwargs)
    return etymology.ico_suffix(key, _steps(*templates), **defaults)


def test_ico_suffix_by_morphology() -> None:
    medico = _ico(
        "medico", _t("inh", **{"1": "it", "2": "la", "3": "medicus"})
    )
    assert not medico.target_is_suffix  # simplex inherited Latin noun
    pubblico = _ico(
        "pubblico", _t("bor+", **{"1": "it", "2": "la", "3": "pūblicus"})
    )
    assert pubblico.target_is_suffix and pubblico.evidence == "latin_icus"
    assert pubblico.base == ""
    elettrico = _ico(
        "elettrico", _t("der", **{"1": "it", "2": "la", "3": "ēlectricus"})
    )
    assert elettrico.target_is_suffix and elettrico.base == "elettro"
    template = _ico(
        "bancario", _t("ety", **{"1": "it", "2": ":af", "3": "x", "4": "-ico"})
    )
    assert not template.target_is_suffix  # does not end in -ico
    tipico = _ico(
        "artistico",
        _t("ety", **{"1": "it", "2": ":af", "3": "arte", "4": "-istico"}),
    )
    assert tipico.target_is_suffix and tipico.evidence == "template"
    assert tipico.raw_suffix == "-istico"
    carico = _ico(
        "carico",
        _t("ety", **{"1": "it", "2": ":from", "3": "caricare"}),
        deverbal=True,
    )
    assert not carico.target_is_suffix
    amico = _ico(
        "amico",
        _t("bor", **{"1": "it", "2": "la", "3": "amīcus"}),
        ico_stressed=True,
    )
    assert not amico.target_is_suffix  # stressed -ìco
    # arcaico /arˈkaj.ko/: the penult is stressed on the a, the i of
    # -ico is a glide; the suffix is unstressed (lexicon.ico_stressed)
    arcaico = _ico(
        "arcaico",
        _t("der", **{"1": "it", "2": "grc", "3": "ἀρχαϊκός"}),
        ico_stressed=False,
    )
    assert arcaico.target_is_suffix and arcaico.evidence == "latin_icus"
    unknown = _ico(
        "pubblico",
        _t("bor+", **{"1": "it", "2": "la", "3": "pūblicus"}),
        ico_stressed=None,
    )
    assert unknown.target_is_suffix  # unknown stress is allowed
    noun_only = _ico(
        "pubblico",
        _t("bor+", **{"1": "it", "2": "la", "3": "pūblicus"}),
        adjectival=False,
    )
    assert not noun_only.target_is_suffix


def test_suffixoids_and_their_evidence() -> None:
    assert etymology.suffixoid("teologo") == "-logo"
    assert etymology.suffixoid("chirurgo") == "-urgo"
    assert etymology.suffixoid("sarcofago") == "-fago"
    assert etymology.suffixoid("spurgo", deverbal=True) == ""
    assert etymology.suffixoid("purga") == ""
    assert etymology.suffixoid("riloga") == ""
    assert etymology.suffixoid("psicologa", masculine_forms=["psicologo"]) == (
        "-logo"
    )
    assert etymology.suffixoid_evidence("teologo") == ("-logo", "ending")
    assert etymology.suffixoid_evidence(
        "psicologa", masculine_forms=["psicologo"]
    ) == ("-logo", "masculine counterpart")
    template = _steps(_t("af", **{"1": "it", "2": "psico", "3": "-loga"}))
    assert etymology.suffixoid_evidence("psicologa", template) == (
        "-logo",
        "etymology (template)",
    )
    assert etymology.SUFFIXOID_KEYS["-fago"] == "fag"


def test_italian_classes_come_from_the_shared_classifier() -> None:
    # The person/thing rules are lexpipe.suffixoids' (the same for
    # Romanian and Czech, R27); the Italian evidence decides learnedness.
    assert suffixoids.semantic_class("theologian")[0] == "person"
    assert suffixoids.semantic_class("surgeon")[0] == "person"
    assert suffixoids.semantic_class("chiromancer")[0] == "person"
    assert suffixoids.semantic_class("catalogue/catalog; list")[0] == "thing"
    assert suffixoids.semantic_class("sarcophagus")[0] == "thing"
    assert suffixoids.semantic_class("cannibal")[0] == "person"
    learned, basis = suffixoids.learned_suffixoid(
        "drammaturgo", "urg", "dramatist, playwright", evidence="ending"
    )
    assert (learned, basis) == (True, "ending")
    assert suffixoids.learned_suffixoid(
        "spurgo",
        "urg",
        "clearing",
        decision=suffixoids.Decision(learned=False),
        evidence="ending",
    ) == (False, "manual")


def test_stratum_keeps_internal_bases_under_a_surface_analysis() -> None:
    # The shared classifier (lexpipe.etymology): a leading affixation of
    # Italian bases stays internal whatever the surface phrase says.
    altoatesino = _steps(
        _t("compound", **{"1": "it", "2": "alto", "3": "atesino"}),
        _t("der", **{"1": "it", "2": "la", "3": "athesīnus"}),
        _t("surf", **{"1": "it", "2": "Alto Adige", "3": "-ino"}),
    )
    text = (
        "From alto + atesino, from Latin athesīnus. By surface analysis, "
        "Alto Adige + -ino."
    )
    assert etymology.stratum(altoatesino, text) == (
        "internal formation",
        "af",
    )
    greek = _steps(
        _t(
            "af",
            **{"1": "it", "2": "ἀγράμματος", "3": "-ismo", "lang1": "grc"},
        )
    )
    text = "From Ancient Greek ἀγράμματος + -ismo. By surface analysis, a- + -ismo."
    assert etymology.stratum(greek, text)[0] == "borrowed whole"
    surface_only = _steps(
        _t("ety", **{"1": "it", "2": ":af", "3": "ampere", "4": "-metrico"})
    )
    assert etymology.stratum(
        surface_only, "By surface analysis, ampere + -o- + -metrico."
    ) == ("borrowed whole", "surface_analysis")


def test_suffix_and_confix_templates_name_the_suffix() -> None:
    # {{suf}} and {{confix}} write the suffix without its hyphen
    feroico = _ico(
        "feroico", _t("suf", **{"1": "it", "2": "Fær Øer", "3": "ico"})
    )
    assert feroico.target_is_suffix and feroico.evidence == "template"
    assert feroico.raw_suffix == "-ico"
    decaedrico = _ico(
        "decaedrico", _t("confix", **{"1": "it", "2": "deca", "3": "edrico"})
    )
    assert decaedrico.target_is_suffix and decaedrico.raw_suffix == "-edrico"
    # a prefix template's last part is a base word, not a suffix
    prefixed = _ico(
        "antiamico", _t("prefix", **{"1": "it", "2": "anti", "3": "amico"})
    )
    assert not prefixed.target_is_suffix
    assert etymology.ico_template(
        _steps(_t("af", **{"1": "it", "2": "bosone", "3": "-ico"}))
    )
    assert not etymology.ico_template(
        _steps(_t("af", **{"1": "it", "2": "fora", "3": "sacco"}))
    )


def test_ico_base() -> None:
    lemmas = frozenset({"elettro", "tallio", "fitoterapia", "classe"})
    assert etymology.ico_base("elettrico", lemmas) == "elettro"
    assert etymology.ico_base("tallico", lemmas) == "tallio"
    assert etymology.ico_base("fitoterapico", lemmas) == "fitoterapia"
    assert etymology.ico_base("geomatica", lemmas) == ""
    assert etymology.ico_base("parco", lemmas) == ""


def test_back_formations_from_nouns_are_not_deverbal() -> None:
    noun = _steps(_t("back-form", **{"1": "it", "2": "metallurgia"}))
    assert (
        etymology.deverbal_verb(noun, "Back-formation from metallurgia") == ""
    )
    verb = _steps(_t("back-form", **{"1": "it", "2": "condurre"}))
    assert etymology.deverbal_verb(verb, "") == "condurre"
    # the deverbal templates name a verb by definition
    assert (
        etymology.deverbal_verb(
            _steps(_t("it-deverbal", **{"1": "caricare"})), ""
        )
        == "caricare"
    )
    # metallurgo is then a -urgo suffixoid (R27), not a deverbal lookalike
    assert etymology.suffixoid("metallurgo", noun, deverbal=False) == "-urgo"


def test_italian_compound_templates_are_internal_formations() -> None:
    found = _steps(_t("it-verb-obj", **{"1": "forare", "2": "sacco"}))
    assert [(s.kind, s.lang, s.terms) for s in found] == [
        ("af", "it", ("forare", "sacco"))
    ]
    text = "Verb-object compound, composed of fora + sacco."
    assert etymology.stratum(found, text)[0] == "internal formation"
    pair = _steps(_t("it-verb-verb", **{"1": "saliscendi", "2": "x"}))
    assert pair[0].kind == "af"
