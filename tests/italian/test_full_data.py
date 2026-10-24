"""Checks on the full Italian data (gitignored; skipped when absent).

- lexeme_id uniqueness (R4): the ``extract`` projection's identifiers
  are unique and equal to :func:`lexpipe.ids.assign_lexeme_ids` over
  every noun/adjective lemma entry of the reference split;
- the drift test (R12): the port of Module:it-common reproduces
  kaikki's plural for every pure-default noun of the split (no plural
  argument), cell by cell; mismatches are listed;
- the golden expectations hold on the full ``features`` output.
"""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

import pytest

from lexpipe import ids, kaikki
from lexpipe.langs.italian import plurals, store
from tests.italian import helpers

pytestmark = [pytest.mark.it, pytest.mark.needs_data]

SPLIT = helpers.PROJECT_ROOT / "data" / "sources" / "kaikki"
SPLIT_FILE = kaikki.split_path(SPLIT, "it")
DATA = helpers.PROJECT_ROOT / "data" / "italian"
PROJECTION = DATA / "extract" / "italian_en_lexemes.jsonl.gz"
COMPLETE = DATA / "italian_lexicon_complete.csv"


def _need(path: Path) -> None:
    if not path.is_file():
        pytest.skip(f"gitignored data absent: {path}")


def _lemma_entries() -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    with SPLIT_FILE.open("rb") as fh:
        for line in fh:
            if b'"pos": "noun"' not in line and b'"pos": "adj"' not in line:
                continue
            entry = json.loads(line)
            if entry.get("pos") in ("noun", "adj") and kaikki.is_lemma_entry(
                entry
            ):
                out.append(entry)
    return out


@pytest.mark.slow
def test_lexeme_ids_are_unique_and_match_the_core() -> None:
    _need(SPLIT_FILE)
    _need(PROJECTION)
    _, records = store.read_jsonl(PROJECTION)
    found = [str(record["lexeme_id"]) for record in records]
    assert len(found) > 80_000
    assert len(found) == len(set(found))
    expected = ids.assign_lexeme_ids("it", _lemma_entries())
    assert sorted(found) == sorted(expected)


@pytest.mark.slow
def test_drift_every_pure_default_noun() -> None:
    """port(lemma, gender) == kaikki plural for every pure-default noun."""
    _need(SPLIT_FILE)
    statuses: Counter[str] = Counter()
    mismatches: list[str] = []
    with SPLIT_FILE.open("rb") as fh:
        for line in fh:
            if b'"it-noun"' not in line:
                continue
            entry = json.loads(line)
            if entry.get("pos") != "noun" or not kaikki.is_lemma_entry(entry):
                continue
            info = plurals.resolve_plurals(
                entry["word"], entry.get("head_templates") or ()
            )
            result = plurals.drift(info, entry)
            if result is None or result.scope != "pure_default_noun":
                continue
            statuses[result.status] += 1
            if result.status == plurals.DRIFT_MISMATCH:
                mismatches.append(
                    f"{entry['word']}: port {result.resolved} "
                    f"kaikki {result.rendered}"
                )
    print(f"drift test (pure-default nouns): {dict(statuses)}")
    assert statuses[plurals.DRIFT_MATCH] > 45_000
    assert not mismatches, "\n".join(mismatches[:50])
    # Only the gender '?' noun has no default (the module raises).
    assert statuses[plurals.DRIFT_PORT_ERROR] <= 5


def test_golden_rows_on_the_full_output() -> None:
    _need(COMPLETE)
    with COMPLETE.open(encoding="utf-8", newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if r["pos"] == "noun"]
    by_word: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_word.setdefault(row["word"], []).append(row)

    def only(word: str) -> dict[str, str]:
        found = [r for r in by_word[word] if r["plural_status"] == "regular"]
        assert len(found) == 1, word
        return found[0]

    assert only("amico")["outcome_i"] == "palatalized"
    assert only("medico")["outcome_i"] == "palatalized"
    assert only("parco")["outcome_i"] == "velar"
    assert only("lago")["outcome_i"] == "velar"
    assert only("chirurgo")["outcome_i"] == "variable"
    assert only("amica")["outcome_e"] == "velar"
    assert only("luce")["frontstem"] == "True"
    assert only("occhio")["nde_class"] == "occhio"
    assert only("medico")["target_is_suffix"] == "False"
    assert only("pubblico")["target_is_suffix"] == "True"
