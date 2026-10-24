"""Golden expectations on the real Czech analysis outputs (needs_data).

Skipped when the analysis stages have not been run
(``python scripts/lang_pipeline.py --lang cs``).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lexpipe import csvio, licensing, spec
from lexpipe.context import RunContext
from lexpipe.langs.czech import analysis_layout, jsonio
from tests.czech import treebuilder

pytestmark = [pytest.mark.cs, pytest.mark.needs_data]

ROOT = treebuilder.PROJECT_ROOT
PATHS = analysis_layout.analysis_paths(
    RunContext(spec.LANGUAGES["cs"], root=ROOT)
)


def _require(path: Path) -> None:
    if not path.is_file():
        pytest.skip(f"gitignored output absent: {path}")


def _pairs() -> dict[tuple[str, str], dict[str, str]]:
    _require(PATHS.derivations)
    out: dict[tuple[str, str], dict[str, str]] = {}
    for row in jsonio.iter_csv(PATHS.derivations):
        out.setdefault((row["base_lemma"], row["derivative"]), row)
    return out


@pytest.mark.parametrize(
    ("base", "derivative", "suffix", "outcome"),
    [
        ("kniha", "knižní", "-ní/-ný", "pal"),
        ("ruka", "ruční", "-ní/-ný", "pal"),
        ("řeka", "říční", "-ní/-ný", "pal"),
        ("liga", "ligový", "-ový/-ovní", "plain"),
        ("banka", "bankovní", "-ový/-ovní", "plain"),
        ("droga", "drogový", "-ový/-ovní", "plain"),
        ("vzduch", "vzdušný", "-ní/-ný", "pal"),
        ("vzduch", "vzduchový", "-ový/-ovní", "plain"),
        ("strach", "strašný", "-ní/-ný", "pal"),
        ("sníh", "sněžit", "-it", "pal"),
        ("tank", "tankista", "-ista", "plain"),
        ("baroko", "barokní", "-ní/-ný", "exception"),
        ("biolog", "bioložka", "-ka", "pal"),
    ],
)
def test_golden_pairs(
    base: str, derivative: str, suffix: str, outcome: str
) -> None:
    row = _pairs().get((base, derivative))
    assert row is not None, (base, derivative)
    assert (row["suffix"], row["outcome"]) == (suffix, outcome)


def test_false_candidates_are_rejected() -> None:
    pairs = _pairs()
    assert ("droga", "drozní") not in pairs  # drozní is drozdí
    assert ("prak", "pracovat") not in pairs  # pracovat ← práce


def test_rule_tp_invariant() -> None:
    _require(PATHS.tp_long)
    rows = [
        r
        for r in jsonio.iter_csv(PATHS.tp_long)
        if r["table"] == "rule"
        and r["rule"] == "palatalize"
        and r["subset"] == "all"
        and r["tau"] == "ge1"
        and r["target_class"] in ("velar", "r", "coronal_inferred")
        and not r["memo"].endswith("N = e")
    ]
    assert rows
    for row in rows:
        assert int(row["e"]) <= 1, row  # Czech rule-TP e ≈ 0 (§8)


def test_r17_dorsal_front_vowel_rule_fails() -> None:
    _require(PATHS.tp_long)
    row = next(
        r
        for r in jsonio.iter_csv(PATHS.tp_long)
        if r["table"] == "R17a"
        and r["trigger"] == "-e/-i"
        and r["rule"] == "palatalize"
        and r["subset"] == "all"
        and r["tau"] == "ge1"
        and r["evidence_set"] == "E0"
        and r["coding"] == "dominant"
        and r["V"] == ""
    )
    assert row["tolerated_trigger"] == "False"


def test_tracked_files_pass_the_licence_gate() -> None:
    for path in (
        PATHS.derivations_summary,
        PATHS.tp_long,
        PATHS.learners,
        PATHS.logtable,
    ):
        _require(path)
        for row in csvio.read_rows(path):
            licensing.check_redistributable(
                licensing.split_sources(row["source"])
            )


def test_learners_ran_with_the_full_frame() -> None:
    _require(PATHS.learners)
    rows = csvio.read_rows(PATHS.learners)
    assert {r["V"] for r in rows} == {"500", "1000", "2000", "5000"}
    assert {r["n_samples"] for r in rows} == {"1000"}
