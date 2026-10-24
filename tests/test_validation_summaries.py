#!/usr/bin/env python3
"""The tracked validator summaries agree with the per-variant reports.

``scripts/validate_lang_rules.py --lang it|cs [--variant NAME]`` writes
one report per specification policy
(``analysis/<language>/<language>_rule_validation[_<variant>].txt``),
and ``python -m phonology.analyses.<language>_palatalization`` writes
the cross-variant summary
(``analysis/<language>/<language>_rule_validation_summary.txt``) that
``docs/<language>/phonology.md`` quotes. Neither is a pipeline stage
with a manifest, so nothing re-runs the summary when the cell data
change: the summaries once lagged the reports by a data-layer re-run
(Italian 7,339 rows against 7,338; Czech 46,449 against 48,143), and
the phonology notes quoted the stale numbers while
``analysis/comparison/README.md`` quoted the current ones.

These tests compare tracked files only (no data needed): each
variant's row count and boolean and IPA match rates in the summary
must equal its report's header, and each phonology note must quote
the current row count.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent

ITALIAN_VARIANTS = ("lexical", "all_underspecified", "stress")
CZECH_VARIANTS = ("lexical", "all_underspecified", "lexical_loo", "spellout")

_ROWS_RE = re.compile(r"Target-stem rows validated: ([\d,]+)")
_BOOL_RE = re.compile(r"BOOLEAN mutation match: ([\d,]+)/([\d,]+)")
_IPA_RE = re.compile(r"SURFACE-IPA match: ([\d,]+)/([\d,]+)")


def _int(text: str) -> int:
    return int(text.replace(",", ""))


def _report_header(path: Path) -> tuple[int, str, str]:
    """(rows, boolean %, IPA %) of a validator report, rates formatted
    as the summaries print them."""
    head = "\n".join(path.read_text(encoding="utf-8").splitlines()[:8])
    rows = re.search(_ROWS_RE, head)
    boolean = re.search(_BOOL_RE, head)
    ipa = re.search(_IPA_RE, head)
    assert rows and boolean and ipa, f"{path}: no validator header"
    n = _int(rows.group(1))
    assert _int(boolean.group(2)) == n and _int(ipa.group(2)) == n
    return (
        n,
        f"{100 * _int(boolean.group(1)) / n:.2f}%",
        f"{100 * _int(ipa.group(1)) / n:.2f}%",
    )


def _report_path(language: str, variant: str) -> Path:
    suffix = "" if variant == "lexical" else f"_{variant}"
    return (
        PROJECT_ROOT
        / "analysis"
        / language
        / f"{language}_rule_validation{suffix}.txt"
    )


def _summary(language: str) -> str:
    path = (
        PROJECT_ROOT
        / "analysis"
        / language
        / f"{language}_rule_validation_summary.txt"
    )
    return path.read_text(encoding="utf-8")


def _italian_overall(summary: str) -> dict[str, tuple[int, str, str]]:
    """``{variant: (rows, bool %, IPA %)}`` of the summary's Overall
    section."""
    found: dict[str, tuple[int, str, str]] = {}
    for line in summary.splitlines():
        m = re.fullmatch(
            r"  (\w+)\s+(\d+)\s+([\d.]+%)\s+([\d.]+%)\s+UR not built: \d+",
            line,
        )
        if m and m.group(1) in ITALIAN_VARIANTS:
            found[m.group(1)] = (int(m.group(2)), m.group(3), m.group(4))
    return found


def _czech_policies(summary: str) -> dict[str, tuple[int, str, str]]:
    """``{variant: (rows, bool %, IPA %)}`` of the summary's policy
    sections."""
    found: dict[str, tuple[int, str, str]] = {}
    policy = ""
    for line in summary.splitlines():
        if m := re.fullmatch(r"Policy: (\w+)", line):
            policy = m.group(1)
            continue
        m = re.fullmatch(
            r"rows ([\d,]+); boolean match\s+([\d.]+%); "
            r"IPA match\s+([\d.]+%); both\s+[\d.]+%",
            line,
        )
        if m and policy and policy not in found:
            found[policy] = (_int(m.group(1)), m.group(2), m.group(3))
    return found


@pytest.mark.it
def test_italian_summary_matches_the_variant_reports() -> None:
    overall = _italian_overall(_summary("italian"))
    assert set(overall) == set(ITALIAN_VARIANTS)
    for variant in ITALIAN_VARIANTS:
        report = _report_header(_report_path("italian", variant))
        assert overall[variant] == report, (
            f"italian_rule_validation_summary.txt ({variant}) is stale: "
            "re-run PYTHONPATH=src python -m "
            "phonology.analyses.italian_palatalization"
        )


@pytest.mark.cs
def test_czech_summary_matches_the_variant_reports() -> None:
    policies = _czech_policies(_summary("czech"))
    assert set(policies) == set(CZECH_VARIANTS)
    for variant in CZECH_VARIANTS:
        report = _report_header(_report_path("czech", variant))
        assert policies[variant] == report, (
            f"czech_rule_validation_summary.txt ({variant}) is stale: "
            "re-run PYTHONPATH=src python -m "
            "phonology.analyses.czech_palatalization"
        )


@pytest.mark.parametrize(
    ("language", "pattern"),
    [
        pytest.param(
            "italian",
            r"## 9\. Results \(real data, ([\d,]+) rows\)",
            marks=pytest.mark.it,
        ),
        pytest.param(
            "czech",
            r"That gives ([\d,]+) rows",
            marks=pytest.mark.cs,
        ),
    ],
)
def test_phonology_note_quotes_the_current_row_count(
    language: str, pattern: str
) -> None:
    rows, _, _ = _report_header(_report_path(language, "lexical"))
    note = (PROJECT_ROOT / "docs" / language / "phonology.md").read_text(
        encoding="utf-8"
    )
    quoted = re.search(pattern, note)
    assert quoted is not None, f"docs/{language}/phonology.md: no row count"
    assert _int(quoted.group(1)) == rows
