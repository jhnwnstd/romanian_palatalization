"""Golden checks on the full Italian analysis outputs (needs_data).

Skipped when the gitignored outputs under ``data/italian/`` are absent.
The expectations come from the design (§5) and the it_phon.json
research and verification:

- verbs: -are never palatalizes; root-stressed -ere alternates; -isc-
  verbs with a velar keep it; tacere over-applies;
- derivatives: the fact-check cases are coded as the verification says
  (grecizzare is not a derivative of greco; musicista and politichese
  have feminine bases; porcino is inherited);
- TP: the feminine -e cell is N = e for palatalize (only the editing
  errors palatalize); coronals never palatalize before -i; the
  stress-nested grammar is tolerated in both branches at τ =
  dictionary;
- the tracked outputs pass the licence gate (R9).
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import pytest

from lexpipe import licensing
from tests.italian import helpers

pytestmark = [pytest.mark.it, pytest.mark.needs_data]

DATA = helpers.PROJECT_ROOT / "data" / "italian"
ANALYSIS = helpers.PROJECT_ROOT / "analysis" / "italian"


def _read(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        pytest.skip(f"gitignored data absent: {path}")
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def test_verb_classes_on_the_full_data() -> None:
    first: dict[str, dict[str, str]] = {}
    for row in _read(DATA / "italian_verbs.csv"):
        first.setdefault(row["verb_key"], row)
    kinds = Counter(
        (row["conj_class"], row["paradigm_type"])
        for row in first.values()
        if row["in_domain"] == "True"
    )
    are = [k for k in kinds if k[0] == "are"]
    assert are == [("are", "velar_throughout")]  # -are: never
    assert kinds[("ere_root", "alternating")] > 200
    assert kinds[("isc", "velar_throughout")] > 20
    assert first["tacere"]["paradigm_type"] == "overapplying"
    assert first["dire"]["conj_class"] == "dire"


def test_fact_check_cases_on_the_full_data() -> None:
    rows = _read(DATA / "italian_derivations.csv")

    def links(word: str) -> list[dict[str, str]]:
        return [r for r in rows if r["derivative"] == word]

    grecizzare = links("grecizzare")
    assert grecizzare and all(
        r["link_status"] == "rejected:named_other_base" for r in grecizzare
    )
    musicista = [
        r for r in links("musicista") if r["link_status"] == "accepted"
    ]
    assert [r["base_word"] for r in musicista] == ["musica"]
    assert musicista[0]["masc_counterpart"] == "other_root"
    politichese = links("politichese")
    assert [r["base_word"] for r in politichese] == ["politica"]
    porcino = [r for r in links("porcino") if r["link_status"] == "accepted"]
    assert porcino and porcino[0]["deriv_stratum"] == "inherited"
    leghista = [r for r in links("leghista") if r["base_word"] == "lega"]
    assert leghista and leghista[0]["edge"] == "velar"


def _tp(**where: str) -> dict[str, str]:
    defaults = {"evidence_set": "E0", "coding": "dominant", "V": ""}
    defaults.update(where)
    rows = [
        r
        for r in _read(ANALYSIS / "italian_tp_long.csv")
        if all(r[k] == v for k, v in defaults.items())
    ]
    assert len(rows) == 1, (where, len(rows))
    return rows[0]


def test_tp_golden_rows() -> None:
    e_cell = _tp(
        table="R17b",
        trigger="-e",
        subset="all",
        tau="dictionary",
        rule="retain",
    )
    assert int(e_cell["e"]) <= 20  # the -ice/-ige editing errors only
    coronal = _tp(
        table="R17c",
        trigger="-i",
        target_class="coronal",
        subset="all",
        tau="ge1",
        rule="palatalize",
    )
    assert coronal["N_trigger"] == coronal["e"]
    for subset, rule in (
        ("nested:post_tonic", "retain"),
        ("nested:elsewhere", "palatalize"),
    ):
        row = [
            r
            for r in _read(ANALYSIS / "italian_tp_long.csv")
            if r["table"] == "it_stress"
            and r["target_class"] == "dorsal"
            and r["tau"] == "dictionary"
            and r["evidence_set"] == "E0"
            and r["subset"] == subset
            and r["rule"] == rule
            and r["memo"].endswith("giavazzi")
        ]
        assert len(row) == 1 and row[0]["tolerated_trigger"] == "True"


def test_tracked_outputs_pass_the_licence_gate() -> None:
    for name in (
        "italian_tp_long.csv",
        "italian_learners.csv",
        "italian_logtable.csv",
    ):
        for row in _read(ANALYSIS / name):
            licensing.check_redistributable(
                licensing.split_sources(row["source"])
            )
            if row.get("licence_class"):
                licensing.check_licence_class(row["licence_class"])
    learners = _read(ANALYSIS / "italian_learners.csv")
    assert {r["n_samples"] for r in learners} == {"1000"}
