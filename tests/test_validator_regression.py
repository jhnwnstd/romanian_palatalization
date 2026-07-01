#!/usr/bin/env python3
"""Regression pins for the palatalization validator's headline numbers.

If someone edits the rules, patches, or underspecs and drops below
the historical bars, these tests fail loudly. If they IMPROVE the
numbers, bump the constants.

Uses the Python API (``phonology.validation``) — no CLI needed.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from phonology.analyses.romanian_palatalization import build_profile  # noqa: E402
from phonology.lexicon import iter_lexicon_rows  # noqa: E402
from phonology.validation import ValidationReport, run_validation  # noqa: E402


CSV_PATH = PROJECT_ROOT / "data" / "romanian_lexicon_with_freq.csv"


@pytest.fixture(scope="module")
def report() -> ValidationReport:
    if not CSV_PATH.exists():
        pytest.skip(f"lexicon CSV missing at {CSV_PATH}")
    profile = build_profile()
    rows = iter_lexicon_rows(CSV_PATH, target_stems=profile.target_stems)
    return run_validation(profile, rows)


def test_row_count(report: ValidationReport) -> None:
    """The lexicon should still yield roughly 11,641 target stems."""
    assert 11_600 <= report.total <= 11_700


def test_boolean_match_rate(report: ValidationReport) -> None:
    """Overall boolean-mutation match must stay ≥97.25%."""
    rate = 100 * report.boolean_matches / report.total
    assert rate >= 97.25, (
        f"Boolean match rate regressed to {rate:.2f}% (was ≥97.25%)"
    )


def test_surface_ipa_match_rate(report: ValidationReport) -> None:
    """Surface-IPA match must stay ≥93% (after normalisation)."""
    rate = 100 * report.ipa_matches / report.total
    assert rate >= 93.0, (
        f"Surface-IPA match rate regressed to {rate:.2f}% (was ≥93%)"
    )


def test_all_cluster_cells_100_percent(report: ValidationReport) -> None:
    """Every cluster cell (ct+e, ct+i, sc+e, sc+i, st+i) must be
    at 100% — the paper's headline claim (latex.tex:611-623).
    """
    expected_cells = {"ct(ă)+e", "ct(ă)+i", "sc(ă)+e", "sc(ă)+i", "st(ă)+i"}
    for cell in expected_cells:
        counts = report.by_cluster.get(cell)
        assert counts is not None, f"cluster cell {cell!r} has zero rows"
        total = counts["match"] + counts["mismatch"]
        assert total > 0, f"cluster cell {cell!r} has zero rows"
        assert counts["mismatch"] == 0, (
            f"cluster cell {cell!r} regressed: "
            f"{counts['match']} match / {counts['mismatch']} mismatch"
        )
