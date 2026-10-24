"""Unit tests of the Italian analyses on synthetic root records.

Covers the observation semantics of :mod:`lexpipe.langs.italian.roots`
(trigger cells, evidence sets, scopes, dedup), the equivalence of the
learner cells' status-based verdicts (:func:`lexpipe.learners.status_verdict`
over :func:`lexpipe.langs.italian.tplearners.statuses_of`) with
:func:`lexpipe.learners.verdict`, the draw-sharing of
:func:`lexpipe.learners.simulate_all` for the Italian cells, and the
units of the spell-out (R20), co-variation (R23) and no-evidence (R21)
analyses.
"""

from __future__ import annotations

import dataclasses
import random
from collections.abc import Sequence
from pathlib import Path

import pytest

from lexpipe import learners, tp
from lexpipe.attest import TAUS
from lexpipe.langs.italian import (
    covariation,
    docstage,
    noevidence,
    roots,
    spellout,
    tplearners,
    tptables,
)

pytestmark = pytest.mark.it

ALL_TAUS = dict.fromkeys(TAUS, True)
DICT_ONLY = {tau: tau == "dictionary" for tau in TAUS}


def _real(
    key: str,
    vowel: str,
    palatal: bool | None,
    *,
    evidence: str = "E0",
    scope: str = "noun",
    count: int = 10,
    attested: dict[str, bool] | None = None,
    ambiguous: bool = False,
) -> roots.Realization:
    return roots.Realization(
        evidence=evidence,
        scope=scope,
        source="pl",
        key=key,
        vowel=vowel,
        palatal=palatal,
        count_cs=count,
        count_folded=count,
        attested=attested or ALL_TAUS,
        ambiguous=ambiguous,
    )


def _record(
    root_id: str,
    items: Sequence[roots.Realization],
    *,
    post_tonic: bool | None = None,
    ico: bool = False,
    stratum: str = "internal formation",
    frequency: int = 10,
    target: str = "dorsal",
) -> roots.RootRecord:
    return roots.RootRecord(
        root_id=root_id,
        lemma=root_id.split(":")[1] + "o",
        stem=root_id.split(":")[1],
        target_class=target,
        seg="c",
        target_letters="c",
        has_noun=True,
        has_adj=False,
        velar_post_tonic=post_tonic,
        stress_class="",
        ico_suffix=ico,
        iaco_suffix=False,
        suffixoid="",
        semantic_class="",
        orth_ico=ico,
        deverbal_from_are=False,
        shares_1st_conj=False,
        stratum=stratum,
        freq_lemma={"all": frequency, "noun": frequency, "adj": 0},
        freq_paradigm=frequency,
        realizations=list(items),
    )


# ---------------------------------------------------------------------------
# Observations
# ---------------------------------------------------------------------------


def test_trigger_cells_and_statuses() -> None:
    amico = _record(
        "it:amic",
        [_real("amici", "i", True), _real("amiche", "e", False)],
    )
    found = roots.observation(amico, "E0")
    status = {
        t: tp.root_status(found, t, "ge1", "dominant") for t in roots.TRIGGERS
    }
    assert status["-i"] == tp.RootStatus.ALTERNATES
    assert status["-e"] == tp.RootStatus.NONALT
    # both cells pooled (R17a): the -i and -e cells disagree, which is not
    # a doublet (R19 codes variants of one cell) but an exception to both
    # rules, whatever the tokens and the coding
    for coding in tp.CODINGS:
        assert tp.root_status(found, "-i/-e", "ge1", coding) == (
            tp.RootStatus.VARIABLE
        )
    assert [obs.cell for obs in found.triggers["-i/-e"]] == [
        "plural:e",
        "plural:i",
    ]
    banca = _record("it:banc", [_real("banche", "e", False)])
    shielded = roots.observation(banca, "E0")
    assert tp.root_status(shielded, "-i", "ge1", "dominant") == (
        tp.RootStatus.SHIELDED
    )


def test_evidence_sets_and_scope() -> None:
    record = _record(
        "it:traffic",
        [
            _real("traffici", "i", True, count=20),
            _real(
                "traffichi", "i", False, evidence="E1", scope="verb", count=0
            ),
            _real("trafficone", "back", None, evidence="E2", scope="deriv"),
            _real("traffici", "i", True, scope="adj"),
        ],
    )
    e0 = roots.realizations_for(record.realizations, "E0")
    assert [r.key for r in e0] == ["traffici"]  # deduplicated
    e1 = roots.realizations_for(record.realizations, "E1")
    assert {r.key for r in e1} == {"traffici", "traffichi"}
    assert len(roots.realizations_for(record.realizations, "E2")) == 3
    adj = roots.realizations_for(record.realizations, "E0", "adj")
    assert [r.scope for r in adj] == ["adj"]
    # E1 adds a velar verb form to a palatal plural: two cells that
    # disagree, so the root is an exception to both rules under every
    # coding (not a doublet that traffici's tokens would decide)
    found = roots.observation(record, "E1")
    for coding in tp.CODINGS:
        assert tp.root_status(found, "-i", "dictionary", coding) == (
            tp.RootStatus.VARIABLE
        )


def test_copies_of_a_form_merge_their_explicit_flag() -> None:
    # nemico: the noun's plural nemici was typed by an editor, the
    # adjective's is a Module:it-common default. The adjective's copy
    # (cell m.pl) sorts before the noun's (cell pl), so keeping the
    # first copy dropped the typed plural from τ = dictionary_explicit.
    typed = {**DICT_ONLY, roots.EXPLICIT_TAU: True}
    default = {**DICT_ONLY, roots.EXPLICIT_TAU: False}
    adjective = dataclasses.replace(
        _real("nemici", "i", True, scope="adj", attested=default),
        source="m.pl",
    )
    noun = _real("nemici", "i", True, scope="noun", attested=typed)
    record = _record("it:nemic", [noun, adjective], post_tonic=True)
    merged = roots.realizations_for(record.realizations, "E0")
    assert [r.key for r in merged] == ["nemici"]  # still counted once
    assert merged[0].attested[roots.EXPLICIT_TAU] is True
    assert merged[0].attested["ge1"] is False
    found = roots.observation(record, "E0")
    assert tp.root_status(found, "-i", roots.EXPLICIT_TAU, "dominant") == (
        tp.RootStatus.ALTERNATES
    )
    # the noun scope sees the noun's copy alone
    nouns = roots.realizations_for(record.realizations, "E0", "noun")
    assert nouns[0].attested[roots.EXPLICIT_TAU] is True
    # an E1 copy of an E0 form is not merged: earlier evidence wins
    verb = _real(
        "nemici",
        "i",
        True,
        evidence="E1",
        scope="verb",
        attested={**ALL_TAUS, roots.EXPLICIT_TAU: False},
    )
    mixed = _record("it:nemic", [adjective, verb])
    e1 = roots.realizations_for(mixed.realizations, "E1")
    assert [r.evidence for r in e1] == ["E0"]
    assert e1[0].attested["ge1"] is False
    assert tp.root_status(
        roots.observation(record, "E0"), "-i", "dictionary", "dominant"
    ) == (tp.RootStatus.ALTERNATES)


def test_unambiguous_attestation() -> None:
    record = _record(
        "it:caric", [_real("carichi", "i", False, ambiguous=True)]
    )
    plain = roots.observation(record, "E0")
    strict = roots.observation(record, "E0", unambiguous=True)
    assert (
        tp.root_status(plain, "-i", "ge1", "dominant") == tp.RootStatus.NONALT
    )
    assert tp.root_status(strict, "-i", "ge1", "dominant") == (
        tp.RootStatus.NO_EVIDENCE
    )
    assert tp.root_status(strict, "-i", "dictionary", "dominant") == (
        tp.RootStatus.NONALT
    )


# ---------------------------------------------------------------------------
# Learners
# ---------------------------------------------------------------------------


def _random_roots(n: int, seed: int) -> list[roots.RootRecord]:
    rng = random.Random(seed)
    out: list[roots.RootRecord] = []
    for index in range(n):
        items: list[roots.Realization] = []
        for _ in range(rng.randint(0, 2)):
            vowel = rng.choice(["i", "e", "back"])
            palatal = None if vowel == "back" else rng.random() < 0.6
            items.append(
                _real(
                    f"f{index}{len(items)}",
                    vowel,
                    palatal,
                    count=rng.randint(0, 30),
                    attested=ALL_TAUS if rng.random() < 0.7 else DICT_ONLY,
                )
            )
        out.append(
            _record(
                f"it:r{index:03d}",
                items,
                post_tonic=rng.choice([True, False, None]),
                frequency=rng.randint(1, 100),
            )
        )
    return out


@pytest.mark.parametrize("trigger", roots.TRIGGERS)
@pytest.mark.parametrize("tau", ["dictionary", "ge1"])
def test_status_verdict_equals_the_core_verdict(
    trigger: str, tau: str
) -> None:
    records = _random_roots(160, seed=11)
    observed = [roots.observation(r, "E0") for r in records]
    spec_cell = tptables.cell("t", "dorsal", trigger, tau=tau)
    post = frozenset(r.root_id for r in records if r.velar_post_tonic)
    specific_core = [
        tp.NestedRule("post", "retain", lambda root: root.root_id in post)
    ]
    core = learners.verdict(spec_cell, specific=specific_core)
    fast = learners.status_verdict(
        tplearners.statuses_of(observed, spec_cell), [("retain", post)]
    )
    plain_core = learners.verdict(spec_cell)
    plain_fast = learners.status_verdict(
        tplearners.statuses_of(observed, spec_cell)
    )
    rng = random.Random(5)
    for _ in range(300):
        vocabulary = rng.sample(observed, rng.randint(0, 60))
        assert fast(vocabulary) == core(vocabulary)
        assert plain_fast(vocabulary) == plain_core(vocabulary)


def test_simulate_all_equals_separate_calls() -> None:
    records = _random_roots(80, seed=3)
    pool = [roots.observation(r, "E0") for r in records]
    cells = [
        tptables.cell("t", "dorsal", "-i", tau="dictionary"),
        tptables.cell("t", "dorsal", "-e", tau="ge1"),
    ]
    evaluators = {
        str(index): learners.status_verdict(tplearners.statuses_of(pool, c))
        for index, c in enumerate(cells)
    }
    for workers in (1, 2):
        shared = learners.simulate_all(
            pool, evaluators, (10, 20), 30, seed=7, workers=workers
        )
        for name, evaluate in evaluators.items():
            alone = learners.simulate_learners(
                pool, evaluate, (10, 20), 30, seed=7
            )
            assert [s.counts for s in shared[name]] == [
                s.counts for s in alone
            ]


# ---------------------------------------------------------------------------
# Spell-out, co-variation, no-evidence units
# ---------------------------------------------------------------------------


def _derivative(
    root_id: str,
    word: str,
    suffix: str,
    palatalized: bool | None,
    *,
    stratum: str = "internal formation",
    tier: str = "named",
    attested: dict[str, bool] | None = None,
    group: str = "learned",
    feminine: bool = False,
) -> spellout.Derivative:
    return spellout.Derivative(
        base_root_id=root_id,
        word=word,
        key=word,
        suffix=suffix,
        group=group,
        palatalized=palatalized,
        edge="palatal" if palatalized else "velar",
        tier=tier,
        deriv_stratum=stratum,
        base_stratum="borrowed whole",
        count_cs=5,
        attested=attested or ALL_TAUS,
        base_lexeme_id=f"{root_id}:lex",
        base_word=word[:4],
        base_pos="noun",
        base_gender="f" if feminine else "m",
        masc_counterpart="",
        ic_suffix=False,
        fact_check="",
    )


def _records(*records: roots.RootRecord) -> roots.Records:
    return roots.Records({r.root_id: r for r in records}, {}, {})


def test_spellout_counts_certain_roots() -> None:
    tipic = _record("it:tipic", [_real("tipici", "i", True)])
    politic = _record("it:politic", [_real("politici", "i", True)])
    antic = _record("it:antic", [_real("antichi", "i", False)])
    items = [
        _derivative("it:tipic", "tipicizzare", "-izzare", True),
        _derivative("it:politic", "politicheggiare", "-eggiare", False),
        _derivative("it:antic", "antichizzare", "-izzare", False),
    ]
    rows = spellout.spellout_rows(_records(tipic, politic, antic), items)
    pal = next(
        r
        for r in rows
        if r["trigger"] == "verbalizer"
        and r["subset"] == "all"
        and r["tau"] == "dictionary"
        and r["coding"] == "dominant"
        and r["rule"] == "palatalize"
    )
    assert (pal["N_target"], pal["N_trigger"], pal["e"]) == (2, 2, 1)
    contrast = next(
        r
        for r in rows
        if r["trigger"] == "verbalizer"
        and r["subset"] == "contrast:nonalt"
        and r["tau"] == "dictionary"
        and r["rule"] == "palatalize"
    )
    assert (contrast["N_trigger"], contrast["e"]) == (1, 1)


def test_covariation_units_and_tables() -> None:
    tipic = _record("it:tipic", [_real("tipici", "i", True)])
    antic = _record("it:antic", [_real("antichi", "i", False)])
    lega = _record("it:leg", [_real("leghe", "e", False)])
    items = [
        _derivative("it:tipic", "tipicità", "-ità", True),
        _derivative("it:antic", "antichità", "-ità", False),
        _derivative("it:antic", "anticità", "-ità", True, tier="candidate"),
        _derivative("it:leg", "leghista", "-ista", False),
    ]
    records = _records(tipic, antic, lega)
    statuses = spellout.root_statuses(records)
    units = covariation.build_units(
        records,
        items,
        statuses,
        "dictionary",
        covariation.TIERS["named"],
        ("-ità", "-ista"),
    )
    assert [(u.lemma, u.base_class, u.outcome) for u in units] == [
        ("antico", "N", "N"),
        ("tipico", "P", "P"),
    ]  # lega has no -i plural: R21, not R23
    assert covariation.table_of(units) == ((1, 0), (0, 1))
    with_candidates = covariation.build_units(
        records,
        items,
        statuses,
        "dictionary",
        covariation.TIERS["named+candidates"],
        ("-ità",),
    )
    assert [u.outcome for u in with_candidates] == ["both", "P"]
    assert covariation.table_of(with_candidates, both_as="P") == (
        (1, 0),
        (1, 0),
    )


def test_noevidence_classes() -> None:
    attested_dict = _record(
        "it:music", [_real("musici", "i", True, attested=DICT_ONLY)]
    )
    lega = _record("it:leg", [_real("leghe", "e", False)])
    amico = _record("it:amic", [_real("amici", "i", True)])
    assert noevidence.base_class(attested_dict, "dictionary") == (
        noevidence.ALTERNATING
    )
    assert (
        noevidence.base_class(attested_dict, "ge1") == noevidence.NO_EVIDENCE
    )
    assert noevidence.base_class(lega, "ge1") == noevidence.NO_EVIDENCE
    assert noevidence.base_class(amico, "ge1") == noevidence.ALTERNATING
    items = [
        _derivative("it:leg", "leghista", "-ista", False, tier="candidate"),
        _derivative("it:music", "musicista", "-ista", True, feminine=True),
    ]
    units = noevidence.build_units(
        _records(attested_dict, lega, amico),
        items,
        "ge1",
        noevidence.TIERS["named+candidates"],
    )
    assert [(u.lemma, u.klass, u.outcome, u.subtype) for u in units] == [
        ("lego", noevidence.NO_EVIDENCE, "N", "feminine-only (-e plurals)"),
        ("musico", noevidence.NO_EVIDENCE, "P", "-i plural below ge1"),
    ]
    assert units[1].feminine_named


def test_readme_block_replacement() -> None:
    block = f"{docstage.BEGIN}\nnew\n{docstage.END}\n"
    text = f"intro\n\n{docstage.BEGIN}\nold\n{docstage.END}\n\noutro\n"
    assert docstage.update_readme(text, block) == (
        f"intro\n\n{docstage.BEGIN}\nnew\n{docstage.END}\n\noutro\n"
    )
    assert docstage.update_readme("intro\n", block) == "intro\n\n" + block


# ---------------------------------------------------------------------------
# R26 counts once; typed plurals; E2 strata; det plurals; R23 output
# ---------------------------------------------------------------------------


def _ico_record(
    root_id: str, palatal: bool, count: int = 10
) -> roots.RootRecord:
    stem = root_id.split(":")[1]
    return _record(
        root_id,
        [
            _real(
                f"{stem}i" if palatal else f"{stem}hi",
                "i",
                palatal,
                count=count,
            )
        ],
        post_tonic=False,
        ico=True,
    )


def test_suffix_units_count_once_in_the_tables() -> None:
    plain = [
        _record("it:fuoc", [_real("fuochi", "i", False)], post_tonic=True),
        _record("it:parc", [_real("parchi", "i", False)], post_tonic=True),
        _record("it:medic", [_real("medici", "i", True)], post_tonic=False),
        _record("it:bianc", [_real("bianchi", "i", False)], post_tonic=False),
    ]
    suffix = [_ico_record(f"it:tipic{n}", True) for n in range(6)]
    chosen = plain + suffix
    observed = [roots.observation(r, "E0") for r in chosen]
    kept, collapsed = tptables.collapse(chosen, observed, tptables.r26_unit)
    assert [r.root_id for r in kept] == [r.root_id for r in plain]
    assert [r.root_id for r in collapsed][-1] == "suffix:-ico"
    unit = collapsed[-1]
    assert len(unit.triggers["-i"]) == 6
    assert tp.root_status(unit, "-i", "ge1", "dominant") == (
        tp.RootStatus.ALTERNATES
    )
    spec_cell = tptables.cell("it_stress", "dorsal", "-i", tau="ge1")
    by_name = {n.name: n for n in tptables.NESTED_STRESS}
    per_root = tptables.nested_set_rows(
        chosen, observed, spec_cell, by_name["giavazzi"], "m"
    )
    once = tptables.nested_set_rows(
        chosen, observed, spec_cell, by_name["suffix_once"], "m"
    )
    # elsewhere: medico, bianco and six -ico roots, or one -ico unit
    assert [(r["N_trigger"], r["e"]) for r in per_root] == [(2, 0), (8, 1)]
    assert [(r["N_trigger"], r["e"]) for r in once] == [(2, 0), (3, 1)]


def test_suffix_once_learners_pool_the_drawn_members() -> None:
    plain = [
        _record(
            f"it:fuoc{n}", [_real(f"fuoc{n}hi", "i", False)], post_tonic=True
        )
        for n in range(5)
    ] + [
        _record("it:medic", [_real("medici", "i", True)], post_tonic=False),
        _record("it:bianc", [_real("bianchi", "i", False)], post_tonic=False),
    ]
    suffix = [_ico_record(f"it:tipic{n}", True) for n in range(5)]
    spec_cell = tptables.cell("t", "dorsal", "-i", tau="ge1")
    observed = {r.root_id: roots.observation(r, "E0") for r in plain + suffix}
    statuses = tplearners.statuses_of(
        [observed[r.root_id] for r in plain], spec_cell
    )
    units = {r.root_id: "suffix:-ico" for r in suffix}
    members = {r.root_id: observed[r.root_id] for r in suffix}
    post = frozenset(f"it:fuoc{n}" for n in range(5))
    evaluate = tplearners.suffix_once_verdict(
        statuses, units, members, spec_cell, [("retain", post)]
    )
    per_root = learners.status_verdict(
        tplearners.statuses_of(list(observed.values()), spec_cell),
        [("retain", post)],
    )
    every = list(observed.values())
    # Per root: no pooled rule (6 velar, 6 palatal of 12); post-tonic
    # retain (5, e 0) and elsewhere palatalize (medico, bianco and five
    # -ico roots: e 1 of 7) hold. Once: 5 velars, medico, bianco and the
    # -ico unit; retain has 2 exceptions of 8 and wins outright.
    assert per_root(every) == "conditioned"
    assert evaluate(every) == "retain"
    # a vocabulary without -ico roots is the same either way
    no_ico = [observed[r.root_id] for r in plain]
    assert evaluate(no_ico) == per_root(no_ico)


def test_typed_plurals_and_e2_strata(tmp_path: Path) -> None:
    from lexpipe import csvio

    cells = tmp_path / "cells.csv"
    header = (
        "root_id",
        "opportunity",
        "in_domain",
        "aligned",
        "filtered",
        "adjudicated_out",
        "target",
        "nde_class",
        "frontstem",
        "mutation",
        "pos",
        "cell",
        "form",
        "count_cs",
        "count_folded",
        "form_ambiguous",
        "provenance",
        *(f"attested_{tau}" for tau in TAUS),
    )
    common = {
        "opportunity": "i",
        "in_domain": True,
        "aligned": True,
        "filtered": False,
        "adjudicated_out": False,
        "nde_class": "",
        "pos": "noun",
        "cell": "pl",
        "count_cs": 3,
        "count_folded": 3,
        "form_ambiguous": False,
        **{f"attested_{tau}": True for tau in TAUS},
    }
    rows = [
        {
            **common,
            "root_id": "it:dolc",
            "target": "dorsal",
            "frontstem": False,
            "mutation": False,
            "form": "dolchi",
            "provenance": "explicit",
        },
        {  # a lexically palatal homograph merged into the dorsal root
            **common,
            "root_id": "it:dolc",
            "target": "",
            "frontstem": True,
            "mutation": False,
            "form": "dolci",
            "provenance": "module_default",
        },
    ]
    csvio.write_rows(cells, rows, header)
    record = _record("it:dolc", [])
    assert roots.add_cells({"it:dolc": record}, cells) == 1
    (only,) = record.realizations
    assert only.key == "dolchi"
    assert only.attested[roots.EXPLICIT_TAU] is True
    derivatives = tmp_path / "derivations.csv"
    deriv_header = (
        "base_root_id",
        "link_status",
        "link_tier",
        "edge",
        "suffix_vowel",
        "suffix",
        "derivative",
        "deriv_stratum",
        "count_cs",
        "count_folded",
        *(f"attested_{tau}" for tau in TAUS),
    )
    base = {
        "base_root_id": "it:dolc",
        "link_status": "accepted",
        "link_tier": "named",
        "suffix_vowel": "i",
        "count_cs": 1,
        "count_folded": 1,
        **{f"attested_{tau}": True for tau in TAUS},
    }
    csvio.write_rows(
        derivatives,
        [
            {
                **base,
                "edge": "palatal",
                "suffix": "-ità",
                "derivative": "dolcità",
                "deriv_stratum": "borrowed whole",
            },
            {
                **base,
                "edge": "velar",
                "suffix": "-ista",
                "derivative": "dolchista",
                "deriv_stratum": "internal formation",
            },
        ],
        deriv_header,
    )
    assert roots.add_derivatives({"it:dolc": record}, derivatives) == 2
    e2 = roots.realizations_for(record.realizations, "E2")
    e2_all = roots.realizations_for(record.realizations, roots.E2_ALL)
    assert sorted(r.key for r in e2) == ["dolchi", "dolchista"]
    assert sorted(r.key for r in e2_all) == ["dolchi", "dolchista", "dolcità"]
    found = roots.observation(record, "E0")
    assert (
        tp.root_status(found, "-i", roots.EXPLICIT_TAU, "dominant")
        == tp.RootStatus.NONALT
    )


def test_det_plural_gives_a_no_evidence_base_a_class() -> None:
    poco = _record("it:poc", [])
    item = _derivative("it:poc", "pochissimo", "-issimo", False)
    with_det = dataclasses.replace(
        item,
        base_det_plural_palatal=False,
        base_det_plural_taus=frozenset({"dictionary", "ge1"}),
    )
    assert noevidence.det_plural_class([item], "ge1") == ""
    assert noevidence.det_plural_class([with_det], "ge1") == (
        noevidence.NON_ALTERNATING
    )
    assert noevidence.det_plural_class([with_det], "ge2pm") == ""
    units = noevidence.build_units(
        _records(poco), [with_det], "ge1", noevidence.TIERS["named"]
    )
    assert [(u.klass, u.outcome) for u in units] == [
        (noevidence.NON_ALTERNATING, "N")
    ]


def test_separated_logistic_fit_prints_no_p_values() -> None:
    pytest.importorskip("statsmodels")
    units = [
        covariation.Unit(
            root_id=f"it:r{index}",
            lemma=f"r{index}",
            suffix="-ità" if index % 2 else "-ista",
            stratum="internal formation",
            base_class="P" if index < 12 else "N",
            outcome="P" if index < 12 else ("P" if index % 5 == 0 else "N"),
            words=(f"w{index}",),
            feminine_base=False,
            log_freq=float(index % 7),
            base_stratum="borrowed whole" if index % 3 else "inherited",
        )
        for index in range(30)
    ]
    line = covariation.logistic_line(units)
    # every P base has a P derivative: base_class separates the outcome
    assert "quasi-complete separation on base_class" in line
    assert "base_class " in line and "(se NA, p NA)" in line
    assert "se 0" not in line and "p 0" not in line


def test_palatal_base_rejections_form_a_sensitivity_tier(
    tmp_path: Path,
) -> None:
    from lexpipe import csvio

    path = tmp_path / "derivations.csv"
    header = [
        "base_root_id",
        "derivative",
        "suffix",
        "suffix_group",
        "palatalized",
        "edge",
        "link_tier",
        "link_status",
        "deriv_stratum",
        "base_stratum",
        "count_cs",
        *(f"attested_{tau}" for tau in TAUS),
        "base_lexeme_id",
        "base_word",
        "base_pos",
        "base_gender",
        "masc_counterpart",
        "ic_suffix",
        "fact_check",
        "base_det_plural_palatal",
        "base_det_plural_taus",
    ]
    common: dict[str, object] = {
        "base_root_id": "it:calc",
        "suffix": "-ista",
        "suffix_group": "learned",
        "link_tier": "candidate",
        "deriv_stratum": "unknown",
        "base_stratum": "unknown",
        "count_cs": 3,
        **{f"attested_{tau}": True for tau in TAUS},
        "base_lexeme_id": "x",
        "base_word": "calco",
        "base_pos": "noun",
        "base_gender": "m",
        "masc_counterpart": "",
        "ic_suffix": False,
        "fact_check": "",
        "base_det_plural_palatal": None,
        "base_det_plural_taus": [],
    }
    csvio.write_rows(
        path,
        [
            {
                **common,
                "derivative": "calcista",
                "palatalized": True,
                "edge": "palatal",
                "link_status": "rejected:palatal_base_exists",
            },
            {
                **common,
                "derivative": "calchista",
                "palatalized": False,
                "edge": "velar",
                "link_status": "accepted",
            },
            {
                **common,
                "derivative": "calcetto",
                "palatalized": True,
                "edge": "palatal",
                "link_status": "rejected:capitalized_only",
            },
        ],
        header,
    )
    items = spellout.load_derivatives(path)
    assert [(i.word, i.tier) for i in items] == [
        ("calcista", spellout.OTHER_BASE_TIER),
        ("calchista", "candidate"),
    ]
    assert (
        spellout.OTHER_BASE_TIER not in covariation.TIERS["named+candidates"]
    )
    assert (
        spellout.OTHER_BASE_TIER
        in noevidence.TIERS["named+candidates+other_base"]
    )
