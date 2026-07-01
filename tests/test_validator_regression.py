#!/usr/bin/env python3
"""Regression pin for the palatalization validator's headline numbers.

The workflow-approved plan sets these bars: overall boolean match
≥97.25%, every cluster cell at 100%. If someone edits the rules,
patches, or underspecs and drops below these numbers this test
fails loudly. If they IMPROVE the numbers, bump the constants.

This test is expensive — it runs the validator on all 11,641 target-
stem lexicon rows — so it lives in its own module and is skipped when
the CSV isn't present (e.g., a fresh clone before the pipeline runs).
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from phonology.analyses.romanian_palatalization import (  # noqa: E402
    RULES_FROM_PAPER,
    load_inventory,
)
from phonology.pipeline import RulePipeline  # noqa: E402
import validate_palatalization_rules as vp  # noqa: E402


CSV_PATH: Path = vp.CSV_PATH


@pytest.fixture(scope="module")
def results():
    if not CSV_PATH.exists():
        pytest.skip(f"lexicon CSV missing at {CSV_PATH}")
    inv = load_inventory()
    pipe = RulePipeline(
        rules=RULES_FROM_PAPER,
        resolve=lambda label: inv.segment(label),
    )
    out: list[vp.RowResult] = []
    for row in vp._load_rows():
        out.append(vp._evaluate(row, pipe, inv, strict_ipa=False))
    return out


def test_row_count(results):
    """The lexicon should still yield the paper's 11,641 target stems."""
    assert len(results) >= 11_600
    assert len(results) <= 11_700


def test_boolean_match_rate(results):
    """Overall boolean-mutation match must stay ≥97.25%."""
    total = len(results)
    matched = sum(1 for r in results if r.boolean_match)
    rate = 100 * matched / total
    assert rate >= 97.25, (
        f"Boolean match rate regressed to {rate:.2f}% "
        f"(was ≥97.25%)"
    )


def test_surface_ipa_match_rate(results):
    """Surface-IPA match must stay ≥93% (after normalisation).

    Lower than boolean because feature-level rules can produce the
    right feature bundle in the wrong Unicode form (e.g. tie-bar
    conventions across dialects). The compare stack folds the known
    noise categories; the residual is the actual rule failure set.
    """
    total = len(results)
    matched = sum(1 for r in results if r.ipa_match)
    rate = 100 * matched / total
    assert rate >= 93.0, (
        f"Surface-IPA match rate regressed to {rate:.2f}% "
        f"(was ≥93%)"
    )


def test_all_cluster_cells_100_percent(results):
    """Every cluster cell (ct+e, ct+i, sc+e, sc+i, st+i) must be at
    100% — the paper's headline claim (latex.tex:611-623)."""
    from collections import Counter, defaultdict
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for r in results:
        key = vp._cluster_key(r.row)
        if key is None:
            continue
        counts[key]["match" if r.boolean_match else "mismatch"] += 1
    expected_cells = {"ct(ă)+e", "ct(ă)+i", "sc(ă)+e", "sc(ă)+i", "st(ă)+i"}
    for cell in expected_cells:
        c = counts[cell]
        total = c["match"] + c["mismatch"]
        assert total > 0, f"cluster cell {cell!r} has zero rows"
        assert c["mismatch"] == 0, (
            f"cluster cell {cell!r} regressed: "
            f"{c['match']} match / {c['mismatch']} mismatch"
        )
