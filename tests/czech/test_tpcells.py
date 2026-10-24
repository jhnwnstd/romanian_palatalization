"""TP cells and observations (design §6, R15–R19) on hand-made roots."""

from __future__ import annotations

import pytest

from lexpipe import tp
from lexpipe.langs.czech import roots
from lexpipe.langs.czech import tpcells as tc

pytestmark = pytest.mark.cs

ALL: tuple[bool, ...] = (True, True, True)
NONE_: tuple[bool, ...] = (False, False, False)


def real(
    cell: str,
    form: str,
    ending: str,
    gender: str,
    allomorph: str,
    *,
    mutation: bool | None = False,
    n_kind: str = "",
    count: int = 10,
    attested: tuple[bool, ...] = ALL,
) -> roots.Realization:
    return roots.Realization(
        cell=cell,
        form=form,
        ending=ending,
        gender=gender,
        allomorph=allomorph,
        n_kind=n_kind,
        mutation=mutation,
        aligned=True,
        count=count,
        attested=attested,
        ambiguous=False,
        confounded=False,
    )


def record(
    lemma: str,
    target: str,
    realizations: list[roots.Realization],
    *,
    suffix_target: bool = False,
) -> roots.RootRecord:
    return roots.RootRecord(
        root_id=f"cs:{lemma}:x",
        pos="noun",
        lemma=lemma,
        lexeme_ids=[f"cs:{lemma}:noun:0:x"],
        target_class=target,
        stem=lemma,
        stem_final=lemma[-1],
        cluster="",
        readings=["m.an"],
        suffix_target=suffix_target,
        lemma_suffix="",
        nde_class="",
        stratum="inherited",
        lemma_freq=100.0,
        paradigm_freq=200.0,
        realizations=realizations,
    )


BIOLOG = record(
    "biolog",
    "velar",
    [
        real("NOM.PL", "biologové", "ové", "m.an", "N", n_kind="allomorph"),
        real("VOC.PL", "biologové", "ové", "m.an", "none"),
        real(
            "LOC.PL",
            "biolozích",
            "ích",
            "m.an",
            "P",
            mutation=True,
            attested=NONE_,
            count=0,
        ),
        real("INS.SG", "biologem", "em", "m.an", "N", n_kind="nontrigger"),
    ],
)
KLUK = record(
    "kluk",
    "velar",
    [
        real("NOM.PL", "kluci", "i", "m.an", "P", mutation=True, count=870),
        real("LOC.PL", "klucích", "ích", "m.an", "P", mutation=True),
        real("INS.SG", "klukem", "em", "m.an", "N", n_kind="nontrigger"),
    ],
)
STUDENT = record(
    "student",
    "coronal_inferred",
    [
        real("NOM.PL", "studenti", "i", "m.an", "P", mutation=True),
        real("VOC.PL", "studenti", "i", "m.an", "none", mutation=False),
    ],
)


def _status(
    rec: roots.RootRecord, key: str, obs: tuple[tp.Observation, ...], tau: str
) -> tp.RootStatus:
    root = tp.RootObservation(rec.key, rec.target_class, 1.0, {key: obs})
    return tp.root_status(root, key, tau, "dominant")


def trigger(name: str) -> tc.RuleTrigger:
    return next(t for t in tc.RULE_TRIGGERS if t.name == name)


def test_rule_trigger_p_and_shielding() -> None:
    nom = trigger("NOM.PL m.an -i")
    kluci = tc.rule_observations(KLUK, nom)
    assert (
        _status(KLUK, nom.key, kluci, "dictionary") == tp.RootStatus.ALTERNATES
    )
    shielded = tc.rule_observations(BIOLOG, nom)
    assert [o.allomorph for o in shielded] == ["N"]
    assert _status(BIOLOG, nom.key, shielded, "ge1") == tp.RootStatus.SHIELDED


def test_nontrigger_rows_are_all_exceptions() -> None:
    ins = trigger("INS.SG m/n -em")
    observed = tc.rule_observations(KLUK, ins)
    assert [(o.alternates, o.allomorph) for o in observed] == [(False, "P")]
    assert _status(KLUK, ins.key, observed, "ge1") == tp.RootStatus.NONALT


def test_c_is_evaluated_only_in_the_vocative() -> None:
    assert tc.rule_applies("c", trigger("VOC.SG m.an -e"))
    assert not tc.rule_applies("c", trigger("NOM.PL m.an -i"))
    assert tc.rule_applies("velar", trigger("NOM.PL m.an -i"))


def test_competition_slots() -> None:
    comp = tc.COMPETITION_BY_NAME["LOC.PL m.an"]
    assert comp.slots("velar") == (
        frozenset({"ích"}),
        frozenset({"ách", "ech"}),
    )
    assert comp.slots("labial") == (frozenset({"ech"}), frozenset({"ách"}))
    nompl = tc.COMPETITION_BY_NAME["NOM.PL m.an"]
    observed = tc.competition_observations(BIOLOG, nompl)
    assert [o.alternates for o in observed] == [False]  # chose -ové


def test_leave_one_cell_out() -> None:
    family = [BIOLOG]
    loc_pl = tc.SPELLOUT_CELLS["LOC.PL m.an"]
    nom_pl = tc.SPELLOUT_CELLS["NOM.PL m.an"]
    # biolozích is evidence for NOM.PL in the dictionary...
    assert tc.alternation_evidence(family, "dictionary", "E0", nom_pl)
    # ...but not when LOC.PL is the tested cell, nor at τ ≥ 1 (0 tokens)
    assert not tc.alternation_evidence(family, "dictionary", "E0", loc_pl)
    assert not tc.alternation_evidence(family, "ge1", "E0", nom_pl)


def test_derivatives_count_from_e1_and_e2() -> None:
    noun = record("strach", "velar", [])
    noun.derivatives = [
        roots.DerivativeEvidence(
            "strašit",
            "-it",
            "P",
            "it_ovat",
            "pal",
            True,
            135,
            ALL,
            "linked",
            roots.FRONT_HIGH,
            "inherited",
        ),
        roots.DerivativeEvidence(
            "strašný",
            "-ní/-ný",
            "P",
            "ni_ovy",
            "pal",
            True,
            330,
            ALL,
            "linked",
            roots.CONSONANT,
            "inherited",
        ),
    ]
    empty: frozenset[str] = frozenset()
    assert not tc.alternation_evidence([noun], "ge1", "E0", empty)
    assert tc.alternation_evidence([noun], "ge1", "E1", empty)
    noun.derivatives = noun.derivatives[1:]
    assert not tc.alternation_evidence([noun], "ge1", "E1", empty)
    assert tc.alternation_evidence([noun], "ge1", "E2", empty)


def test_r17_uses_one_reflex_per_form() -> None:
    high = frozenset({roots.FRONT_HIGH})
    observed = tc.r17_observations(STUDENT, high, "E0")
    assert [(o.form, o.alternates) for o in observed] == [("studenti", True)]
    mid = frozenset({roots.FRONT_MID})
    assert [o.alternates for o in tc.r17_observations(KLUK, mid, "E0")] == [
        False
    ]


def test_onsets_and_syllables() -> None:
    assert roots.onset("ích") == roots.FRONT_HIGH
    assert roots.onset("em") == roots.FRONT_MID
    assert roots.onset("y") == roots.Y_VOWEL
    assert roots.onset("ové") == roots.BACK
    assert roots.onset("ka") == roots.CONSONANT
    assert roots.onset("") == roots.NO_ENDING
    assert roots.monosyllable("vlk") and roots.monosyllable("mág")
    assert roots.monosyllable("Dán") and not roots.monosyllable("biolog")
    assert roots.monosyllable("houf")


def test_root_observation_has_every_trigger() -> None:
    observed = tc.root_observation(KLUK, "E0")
    assert "rule:NOM.PL m.an -i" in observed.triggers
    assert "comp:LOC.SG m.in" in observed.triggers
    assert any(key.startswith("r17:") for key in observed.triggers)
    assert observed.frequency == 100.0
    assert tc.root_observation(KLUK, "E0", "paradigm").frequency == 200.0


def test_nom_pl_e_is_a_nontrigger_row() -> None:
    diplomat = record(
        "diplomat",
        "coronal_inferred",
        [
            real("NOM.PL", "diplomaté", "é", "m.an", "N", n_kind="allomorph"),
            real("NOM.PL", "diplomati", "i", "m.an", "P", mutation=True),
        ],
    )
    rows = tc.rule_observations(diplomat, trigger("NOM.PL m.an -é"))
    assert [(o.form, o.alternates) for o in rows] == [("diplomaté", False)]


def test_r17c_infers_the_coronal_change_in_every_front_form() -> None:
    kmen = record(
        "kmen",
        "coronal_inferred",
        [
            real("LOC.SG", "kmeni", "i", "m.in", "N", n_kind="soft"),
            real("INS.SG", "kmenem", "em", "m.in", "N", n_kind="nontrigger"),
        ],
    )
    high = frozenset({roots.FRONT_HIGH})
    mid = frozenset({roots.FRONT_MID})
    assert [o.alternates for o in tc.r17_observations(kmen, high, "E0")] == [
        True
    ]
    assert [o.alternates for o in tc.r17_observations(kmen, mid, "E0")] == [
        False
    ]


def test_baseline_is_hard_declension_only() -> None:
    hard = record("chlap", "labial", [])
    soft = record("učitel", "lateral", [])
    soft.declension = "soft"
    baseline = tc.COMPETITION_TARGETS["baseline"]
    assert baseline(hard) and not baseline(soft)
    assert hard.in_baseline and not soft.in_baseline
    assert not tc.COMPETITION_TARGETS["lateral"](soft)
