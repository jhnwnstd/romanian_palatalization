"""Golden expectations of the audit fixes on the real Czech outputs.

Skipped when the pipeline has not been run
(``python scripts/lang_pipeline.py --lang cs``).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from lexpipe import csvio, spec
from lexpipe.context import RunContext
from lexpipe.langs.czech import analysis_layout, jsonio, layout
from tests.czech import treebuilder

pytestmark = [pytest.mark.cs, pytest.mark.needs_data]

ROOT = treebuilder.PROJECT_ROOT
CTX = RunContext(spec.LANGUAGES["cs"], root=ROOT)
PATHS = layout.paths_for(CTX)
APATHS = analysis_layout.analysis_paths(CTX)


def _require(path: Path) -> None:
    if not path.is_file():
        pytest.skip(f"gitignored output absent: {path}")


def _lexicon(*lemmas: str) -> list[dict[str, str]]:
    _require(PATHS.lexicon_complete)
    wanted = set(lemmas)
    return [
        row
        for row in jsonio.iter_csv(PATHS.lexicon_complete)
        if row["lemma"] in wanted
    ]


def _one(rows: list[dict[str, str]], lemma: str, pos: str = "noun") -> Any:
    found = [r for r in rows if r["lemma"] == lemma and r["pos"] == pos]
    assert found, lemma
    return found[0]


def test_derivational_form_of_entries_are_lexemes() -> None:
    rows = _lexicon("klouček", "Češka", "muška")
    klouček = _one(rows, "klouček")
    assert "diminutive:kluk" in csvio.split_list(
        klouček["derivational_formof"]
    )
    assert klouček["target_is_suffix"] == "True"
    assert _one(rows, "Češka")["suffix_evidence"] in ("formof", "etymology")


def test_suffix_targets_are_morphological() -> None:
    rows = _lexicon("typický", "český", "banka", "ruka", "knížka")
    assert _one(rows, "typický", "adj")["target_is_suffix"] == "True"
    assert _one(rows, "český", "adj")["target_is_suffix"] == "True"
    assert _one(rows, "banka")["target_is_suffix"] == "False"
    assert _one(rows, "knížka")["target_is_suffix"] == "True"


def test_stems_of_the_audit() -> None:
    rows = _lexicon("oheň", "míra", "fleky", "stupeň")
    assert _one(rows, "oheň")["target_class"] == "soft"
    assert _one(rows, "stupeň")["target_class"] == "soft"
    assert _one(rows, "míra")["target_class"] == "r"
    assert _one(rows, "fleky")["target_class"] == "velar"


def test_names_leave_the_domain() -> None:
    rows = _lexicon("Rak", "Býk", "Pha", "rak", "býk")
    for lemma in ("Rak", "Býk", "Pha"):
        assert _one(rows, lemma)["domain"] == "out", lemma
    assert _one(rows, "rak")["domain"] == "in"


def test_homograph_lemma_frequency_is_flagged() -> None:
    rows = _lexicon("jak")
    yak = [r for r in rows if r["pos"] == "noun"]
    assert yak and all(r["lemma_ambiguous"] == "True" for r in yak)
    for row in yak:
        assert int(row["paradigm_sum_unambiguous"] or 0) < int(
            row["lemma_count"] or 0
        )


def _pairs() -> set[tuple[str, str]]:
    _require(APATHS.derivations)
    return {
        (row["base_lemma"], row["derivative"])
        for row in jsonio.iter_csv(APATHS.derivations)
    }


def test_reversed_and_false_pairs_are_dropped() -> None:
    pairs = _pairs()
    for pair in (
        ("obloha", "obložit"),
        ("útok", "útočit"),
        ("tucha", "tušit"),
        ("katolík", "katolicismus"),
        ("ekolog", "ekologista"),
        ("sníh", "Sněžka"),
    ):
        assert pair not in pairs, pair
    assert ("poloha", "polohovat") in pairs  # formed from the noun


def test_homograph_roots_keep_one_copy() -> None:
    _require(APATHS.derivations)
    roots: dict[tuple[str, str, str], set[str]] = {}
    for row in jsonio.iter_csv(APATHS.derivations):
        key = (row["base_lemma"], row["suffix"], row["derivative"])
        roots.setdefault(key, set()).add(row["base_root_id"])
    assert all(len(found) == 1 for found in roots.values())
    sokyne = [
        row
        for row in jsonio.iter_csv(APATHS.derivations)
        if row["derivative"] == "sokyně"
    ]
    assert sokyne and all(
        "m.an" in csvio.split_list(r["base_readings"]) for r in sokyne
    )


def test_strata_split_adapted_loans_and_hybrids() -> None:
    _require(APATHS.lexeme_strata)
    found = {
        row["lemma"]: (row["stratum_wikt"], row["stratum_wikt_basis"])
        for row in jsonio.iter_csv(APATHS.lexeme_strata)
        if row["lemma"] in ("tragický", "fanoušek", "Velšan")
    }
    assert found["tragický"] == ("borrowed whole", "adapted_loan")
    assert found["fanoušek"] == ("internal formation", "hybrid")
    assert found["Velšan"] == ("internal formation", "hybrid")


def test_logtable_things() -> None:
    _require(APATHS.logtable)
    classes = {
        row["lemma"]: row["semantic_class"]
        for row in csvio.read_rows(APATHS.logtable)
    }
    assert classes["dekalog"] == "thing"
    assert classes["paragraf"] == "thing"


def _tp(**want: str) -> dict[str, str]:
    _require(APATHS.tp_long)
    for row in jsonio.iter_csv(APATHS.tp_long):
        if all(row[key] == value for key, value in want.items()):
            return row
    raise AssertionError(want)


def test_adverbs_of_plain_velar_adjectives_palatalize() -> None:
    common = {
        "table": "competition_adj",
        "target_class": "velar",
        "tau": "ge1",
        "evidence_set": "E0",
        "coding": "dominant",
    }
    plain = _tp(
        trigger="ADV -y", subset="class:other", rule="palatalize", **common
    )
    assert plain["tolerated_trigger"] == "True"
    sky = _tp(trigger="ADV -y", subset="class:sky", rule="retain", **common)
    assert sky["tolerated_trigger"] == "True"


def test_strict_rows_exist() -> None:
    row = _tp(
        table="competition",
        target_class="velar",
        trigger="LOC.SG m.in",
        subset="all",
        tau="ge1_strict",
        evidence_set="E0",
        coding="dominant",
        rule="retain",
    )
    assert int(row["N_trigger"]) < 50  # LOC.SG -u is GEN.SG/DAT.SG too


def test_learners_on_three_bases() -> None:
    _require(APATHS.learners)
    rows = csvio.read_rows(APATHS.learners)
    assert {r["freq_basis"] for r in rows} == {
        "lemma",
        "paradigm",
        "paradigm_unambiguous",
    }
    assert {r["tau"] for r in rows} == {"dictionary", "ge1", "ge1_strict"}


def test_counterexample_words_are_untracked() -> None:
    _require(APATHS.noevidence)
    tracked = APATHS.noevidence.read_text(encoding="utf-8")
    assert "Counterexamples at tau = ge1" not in tracked
    assert "Counterexamples at tau = dictionary" in tracked
    _require(APATHS.noevidence_words)
