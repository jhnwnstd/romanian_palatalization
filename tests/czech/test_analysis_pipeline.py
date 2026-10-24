"""The Czech analysis stages on the fixture tree (design §6, R9, R14–R27).

The tree holds the CC BY-SA Wiktionary fixtures plus synthetic DeriNet,
czetyl, MorfFlex, Leipzig and CNK tables (:mod:`tests.czech.analysis_tree`).
"""

from __future__ import annotations

import importlib.util
import subprocess
from collections.abc import Mapping
from pathlib import Path

import pytest

from lexpipe import csvio, licensing
from lexpipe.langs import czech
from lexpipe.langs.czech import analysis_layout, schema
from tests.czech import analysis_tree

pytestmark = pytest.mark.cs


def _paths(root: Path) -> analysis_layout.AnalysisPaths:
    return analysis_layout.analysis_paths(analysis_tree.context(root))


def _pairs(root: Path) -> list[dict[str, str]]:
    return csvio.read_rows(_paths(root).derivations)


def _pair(
    rows: list[dict[str, str]], base: str, derivative: str
) -> dict[str, str]:
    found = [
        r
        for r in rows
        if r["base_lemma"] == base and r["derivative"] == derivative
    ]
    assert found, f"no pair {base} → {derivative}"
    return found[0]


def test_wiktionary_pairs(analysis_root: Path) -> None:
    rows = _pairs(analysis_root)
    knizni = _pair(rows, "kniha", "knižní")
    assert (knizni["suffix"], knizni["outcome"], knizni["basis"]) == (
        "-ní/-ný",
        "pal",
        "linked",
    )
    assert "wikt_ety" in knizni["evidence"]
    assert "derinet" in knizni["evidence"]
    bankovni = _pair(rows, "banka", "bankovní")
    assert (bankovni["suffix_class"], bankovni["outcome"]) == ("N", "plain")
    assert bankovni["suffix_pair"] == "ni_ovy"
    tankista = _pair(rows, "tank", "tankista")
    assert tankista["suffix_class"] == "NT"


def test_derinet_only_and_candidate_pairs(analysis_root: Path) -> None:
    rows = _pairs(analysis_root)
    knihovat = _pair(rows, "kniha", "knihovat")
    assert knihovat["basis"] == "all"  # a DeriNet-only pair
    assert knihovat["derinet_verbalizer"] == "-ovat"
    assert knihovat["verbalizer_consistent"] == "True"
    assert knihovat["conj_class_consistent"] == "True"
    ligovy = _pair(rows, "liga", "ligový")
    assert "cand_headword" in ligovy["evidence"]  # a cs.wiktionary headword


def test_dual_remnant_base(analysis_root: Path) -> None:
    rows = _pairs(analysis_root)
    rucni = _pair(rows, "ruka", "ruční")
    assert rucni["base_domain"] == "out:dual_remnant"
    assert rucni["outcome"] == "pal"


def test_strata(analysis_root: Path) -> None:
    rows = _pairs(analysis_root)
    vzdusny = _pair(rows, "vzduch", "vzdušný")
    assert vzdusny["base_stratum_wikt"] == "borrowed whole"
    assert vzdusny["derivative_stratum_wikt"] == "internal formation"
    liga = _pair(rows, "liga", "ligový")
    # Wiktionary has no etymology for liga; czetyl (synthetic) says loan
    assert liga["base_stratum_wikt"] == "unknown"
    assert (liga["base_stratum"], liga["base_stratum_source"]) == (
        "borrowed whole",
        "czetyl",
    )


def test_stem_with_quantity_alternation(analysis_root: Path) -> None:
    complete = {
        r["lemma"]: r
        for r in csvio.read_rows(
            analysis_root / "data" / "czech" / "czech_lexicon_complete.csv"
        )
        if r["pos"] == "noun"
    }
    assert (complete["dráha"]["stem"], complete["dráha"]["target_class"]) == (
        "dráh",
        "velar",
    )
    assert complete["skála"]["target_class"] == "lateral"


def test_rule_tp_has_no_velar_exceptions(analysis_root: Path) -> None:
    rows = csvio.read_rows(_paths(analysis_root).tp_long)
    velar = [
        r
        for r in rows
        if r["table"] == "rule"
        and r["target_class"] == "velar"
        and r["trigger"] == "NOM.PL m.an -i"
        and r["subset"] == "all"
        and r["rule"] == "palatalize"
    ]
    assert velar and all(r["e"] == "0" for r in velar)
    assert any(int(r["N_trigger"]) > 0 for r in velar)
    nontrigger = [
        r
        for r in rows
        if r["table"] == "rule"
        and r["trigger"] == "INS.SG m/n -em"
        and r["rule"] == "palatalize"
        and r["subset"] == "all"
    ]
    assert nontrigger and all(r["e"] == r["N_trigger"] for r in nontrigger)


def test_log_nouns_form_their_own_nested_class(analysis_root: Path) -> None:
    rows = csvio.read_rows(_paths(analysis_root).tp_long)
    suffixoid = [
        r
        for r in rows
        if r["table"] == "competition"
        and r["target_class"] == "velar"
        and r["trigger"] == "NOM.PL m.an"
        and r["subset"] == "class:suffixoid"
        and r["tau"] == "dictionary"
        and r["coding"] == "dominant"
        and r["rule"] == "retain"
    ]
    # biolog: biologové, the -ové of the -log/-urg loans (IJP 226)
    assert suffixoid and suffixoid[0]["e"] == "0"
    assert int(suffixoid[0]["N_trigger"]) >= 1


def test_tracked_outputs_pass_the_licence_gate(analysis_root: Path) -> None:
    paths = _paths(analysis_root)
    for path in (
        paths.derivations_summary,
        paths.tp_long,
        paths.learners,
        paths.logtable,
    ):
        assert licensing.is_tracked_location(path, analysis_root)
        for row in csvio.read_rows(path):
            licensing.check_redistributable(
                licensing.split_sources(row["source"])
            )
            licensing.check_licence_class(row["licence_class"])


def test_nc_outputs_stay_in_nc(analysis_root: Path) -> None:
    paths = _paths(analysis_root)
    for path in (
        paths.derinet_db,
        paths.derivations,
        paths.lexeme_strata,
        paths.spellout_roots,
        paths.noevidence_pairs,
        paths.covariation_units,
    ):
        assert path.parent.name == "nc" and path.is_file()


def test_tracked_text_reports_have_no_nc_rows(analysis_root: Path) -> None:
    """Reports are aggregates: no DeriNet node ids or czetyl rows."""
    paths = _paths(analysis_root)
    for path in (
        paths.derivations_report,
        paths.tp_summary,
        paths.spellout,
        paths.noevidence,
        paths.covariation,
        paths.logtable_summary,
    ):
        text = path.read_text(encoding="utf-8")
        assert "#NNF" not in text and "is_loan\t" not in text


def test_reports_exist_and_mention_their_rules(analysis_root: Path) -> None:
    paths = _paths(analysis_root)
    assert "R20" in paths.spellout.read_text(encoding="utf-8")
    assert "R21" in paths.noevidence.read_text(encoding="utf-8")
    covariation = paths.covariation.read_text(encoding="utf-8")
    assert "R23" in covariation and "logistic model" in covariation
    if importlib.util.find_spec("statsmodels") is None:
        # The skipped model names the missing optional dependency.
        assert "statsmodels" in covariation
    assert "R18" in paths.tp_summary.read_text(encoding="utf-8")


def test_logtable(analysis_root: Path) -> None:
    rows = {
        r["lemma"]: r for r in csvio.read_rows(_paths(analysis_root).logtable)
    }
    biolog = rows["biolog"]
    assert (biolog["suffixoid"], biolog["semantic_class"]) == (
        "-log",
        "person",
    )
    assert biolog["nom_pl_allomorph"] == "-ové"
    assert biolog["loc_pl_palatalized"] == "True"
    assert "count" not in " ".join(rows["biolog"])  # no corpus columns


def test_headers_match_the_schema(analysis_root: Path) -> None:
    for table in schema.tables():
        path = analysis_root / table.path
        if not path.is_file():
            continue
        header = path.read_text(encoding="utf-8").split("\n", 1)[0]
        assert header.split(",") == list(table.fields), table.name


def test_data_dictionary_is_generated(analysis_root: Path) -> None:
    text = (analysis_root / "docs" / "czech" / "data_dict.md").read_text(
        encoding="utf-8"
    )
    assert text == schema.render()


def test_rerun_is_skipped(analysis_root: Path) -> None:
    """Every analysis stage's manifest matches: a second run skips."""
    ctx = analysis_tree.context(analysis_root)
    messages: list[str] = []
    for name in analysis_tree.ANALYSIS_STAGES:
        before = _digests(analysis_root)
        status = czech.STAGES[name](ctx)
        assert status in (0, None)
        assert _digests(analysis_root) == before
        messages.append(name)
    assert messages == list(analysis_tree.ANALYSIS_STAGES)


def _digests(root: Path) -> Mapping[str, int]:
    out: dict[str, int] = {}
    for path in sorted((root / "analysis" / "czech").glob("*")):
        out[path.name] = path.stat().st_mtime_ns
    return out


def test_tracked_outputs_are_not_gitignored() -> None:
    """The tracked analysis files must not be hidden by .gitignore."""
    root = Path(__file__).resolve().parents[2]
    result = subprocess.run(
        [
            "git",
            "-C",
            str(root),
            "check-ignore",
            "analysis/czech/czech_tp_long.csv",
            "docs/czech/data_dict.md",
            "docs/czech/README.md",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.stdout.strip() == ""


def test_spellout_roots(analysis_root: Path) -> None:
    rows = csvio.read_rows(_paths(analysis_root).spellout_roots)
    feminine = {
        (r["base_lemma"], r["pair"]): r
        for r in rows
        if r["root_set"] == "feminine"
    }
    assert feminine[("kniha", "ni_ovy")]["choice_ge1"] == "P"
    assert feminine[("liga", "ni_ovy")]["choice_ge1"] == "N"
    assert feminine[("banka", "ni_ovy")]["choice_ge1"] == "N"
    # noha has two lexemes of one root: one root type (R4)
    noha = [
        r
        for r in rows
        if r["base_lemma"] == "noha"
        and r["pair"] == "ni_ovy"
        and r["root_set"] == "feminine"
    ]
    assert len(noha) == 1


def test_noevidence_base_classes(analysis_root: Path) -> None:
    rows = csvio.read_rows(_paths(analysis_root).noevidence_pairs)
    classes = {r["base_lemma"]: r for r in rows}
    # biologové is attested, biolozích is not (synthetic Leipzig counts)
    assert classes["biolog"]["base_class_dictionary"] == "alternating"
    assert classes["biolog"]["base_class_ge1"] == "non-alternating"
    assert classes["kniha"]["base_class_ge1"] == "alternating"


def test_covariation_units(analysis_root: Path) -> None:
    rows = csvio.read_rows(_paths(analysis_root).covariation_units)
    jazyk = next(
        r for r in rows if r["lemma"] == "jazyk" and r["pair"] == "ni_ovy"
    )
    # v jazyce (P) but jazykový (the non-palatalizing member)
    assert jazyk["LOC.SG_any_ge1"] == "P"
    assert jazyk["choice_dictionary"] == "N"
