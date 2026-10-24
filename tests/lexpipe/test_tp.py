"""Tests for lexpipe.tp against the Romanian R output.

``tp_row`` must reproduce the rows ``tp_table()`` in
``analysis/analyze_romanian_palatalization.R`` wrote to
``analysis/romanian_tp_summary_{full,downsampled}.csv``: same
columns, same rate, and the same TRUE/FALSE/NA in ``majority?`` and
``tolerated?``.
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path
from typing import Any

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lexpipe import tp  # noqa: E402

TP_FULL = PROJECT_ROOT / "analysis" / "romanian_tp_summary_full.csv"
TP_DOWNSAMPLED = (
    PROJECT_ROOT / "analysis" / "romanian_tp_summary_downsampled.csv"
)
R_LOGICAL = {"TRUE": True, "FALSE": False, "NA": None}


def _r_rows(path: Path) -> list[dict[str, str]]:
    """Data rows of an R summary (section headers have N = NA)."""
    if not path.exists():
        pytest.skip(f"R summary not present: {path}")
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        assert tuple(reader.fieldnames or ()) == tp.TP_FIELDS
        return [row for row in reader if row["N"] != "NA"]


def _r_row(path: Path, label: str) -> dict[str, str]:
    matches = [row for row in _r_rows(path) if row["type"] == label]
    assert len(matches) == 1, f"{label!r} not found once in {path}"
    return matches[0]


def _assert_same_row(row: dict[str, Any], r_row: dict[str, str]) -> None:
    assert row["type"] == r_row["type"]
    assert row["N"] == int(r_row["N"])
    assert row["mutated"] == int(r_row["mutated"])
    assert row["non-mutated"] == int(r_row["non-mutated"])
    if r_row["rate"] == "NA":
        assert row["rate"] is None
    else:
        # readr writes the shortest round-trip repr: exact equality.
        assert row["rate"] == float(r_row["rate"])
    assert row["majority?"] is R_LOGICAL[r_row["majority?"]]
    assert row["tolerated?"] is R_LOGICAL[r_row["tolerated?"]]
    assert row["memo"] == r_row["memo"]


def _from_r(r_row: dict[str, str]) -> dict[str, Any]:
    return tp.tp_row(
        r_row["type"], int(r_row["N"]), int(r_row["mutated"]), r_row["memo"]
    )


# ---------------------------------------------------------------------------
# θ_N
# ---------------------------------------------------------------------------


def test_theta() -> None:
    assert tp.theta(0) == 0.0
    assert tp.theta(1) == 0.0
    assert tp.theta(2) == pytest.approx(2.8853900817779268)
    assert tp.theta(1000) == pytest.approx(144.76482730108395)
    assert tp.theta(2642) == 2642 / math.log(2642)


# ---------------------------------------------------------------------------
# Rows of the Romanian full-lexicon summary
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "label",
    [
        # Productive by majority, not tolerated: e = 601 > θ ≈ 335.
        "<c> + <-i, -e> plural, WITH NDEB",
        # Small cell, tolerated: e = 4 ≤ θ_14 ≈ 5.3.
        "<z> + <-i> plural",
    ],
)
def test_tp_row_matches_romanian_full_summary(label: str) -> None:
    r_row = _r_row(TP_FULL, label)
    row = _from_r(r_row)
    _assert_same_row(row, r_row)
    assert list(row) == list(tp.TP_FIELDS)
    exceptions = row["non-mutated"]
    assert row["tolerated?"] is (exceptions <= tp.theta(row["N"]))


@pytest.mark.parametrize("path", [TP_FULL, TP_DOWNSAMPLED], ids=["full", "ds"])
def test_tp_row_reproduces_every_r_row(path: Path) -> None:
    rows = _r_rows(path)
    assert rows
    for r_row in rows:
        _assert_same_row(_from_r(r_row), r_row)


def test_known_values_without_the_csv() -> None:
    row = tp.tp_row("<c> + <-i, -e> plural, WITH NDEB", 2642, 2041)
    assert row["rate"] == 2041 / 2642
    assert row["majority?"] is True
    assert row["tolerated?"] is False
    assert row["memo"] == ""


# ---------------------------------------------------------------------------
# NA semantics
# ---------------------------------------------------------------------------


def test_empty_cell_is_all_missing() -> None:
    row = tp.tp_row("empty", 0, 0)
    assert (row["rate"], row["majority?"], row["tolerated?"]) == (
        None,
        None,
        None,
    )


def test_single_item_cell_has_no_tolerance_verdict() -> None:
    # R: theta is NA for N <= 1 (cf. "<sc> + <-e> plural (downsampled)").
    for mutated in (0, 1):
        row = tp.tp_row("one", 1, mutated)
        assert row["tolerated?"] is None
        assert row["majority?"] is bool(mutated)


def test_tie_has_no_majority() -> None:
    row = tp.tp_row("tie", 2, 1)
    assert row["majority?"] is None
    assert row["tolerated?"] is True  # e = 1 <= 2 / ln 2


@pytest.mark.parametrize(("n", "mutated"), [(-1, 0), (3, 4), (3, -1)])
def test_invalid_counts(n: int, mutated: int) -> None:
    with pytest.raises(ValueError):
        tp.tp_row("bad", n, mutated)


def test_header_and_blank_rows() -> None:
    header = tp.header_row("Dorsals")
    assert list(header) == list(tp.TP_FIELDS)
    assert header["type"] == "Dorsals" and header["memo"] == ""
    assert all(header[k] is None for k in tp.TP_FIELDS[1:-1])
    assert tp.blank_row()["type"] == ""


def test_small_n_flag() -> None:
    assert tp.SMALL_CELL_THRESHOLD == 20
    assert tp.small_n(19) and not tp.small_n(20)


# ---------------------------------------------------------------------------
# Downsampling
# ---------------------------------------------------------------------------

ROWS = [
    {"lemma": "casa", "plural": "case", "freq": "50"},
    {"lemma": "amico", "plural": "amici", "freq": "70"},
    {"lemma": "bosco", "plural": "boschi", "freq": "70"},
    {"lemma": "casa", "plural": "casi", "freq": "90"},  # max wins
    {"lemma": "rarissimo", "plural": "", "freq": "0"},
    {"lemma": "blank", "plural": "", "freq": ""},
    {"lemma": "dado", "plural": "dadi", "freq": "10"},
]


def test_top_n_ranks_lemmas_by_max_frequency() -> None:
    top = tp.top_n(ROWS, 3, "freq")
    assert [(r["lemma"], r["plural"]) for r in top] == [
        ("casa", "case"),  # every row of a lemma, in input order
        ("casa", "casi"),
        ("amico", "amici"),  # tie at 70 broken by lemma
        ("bosco", "boschi"),
    ]


def test_top_n_never_selects_zero_frequency() -> None:
    lemmas = {r["lemma"] for r in tp.top_n(ROWS, 100, "freq")}
    assert lemmas == {"casa", "amico", "bosco", "dado"}


def test_top_n_with_callable_and_custom_key() -> None:
    rows: list[dict[str, Any]] = [
        {"id": "b", "f": 2.0},
        {"id": "a", "f": 2.0},
        {"id": "c", "f": 5},
    ]
    top = tp.top_n(rows, 2, lambda r: r["f"], lemma_key="id")
    assert [r["id"] for r in top] == ["c", "a"]
    assert tp.top_n(rows, 0, "f", lemma_key="id") == []
