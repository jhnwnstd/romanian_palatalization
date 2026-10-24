"""Helpers: a temporary repository tree holding the fixture inputs.

The tree has the layout the plugin reads (``data/romanian_lexicon_with_freq.csv``,
``data/verb_paradigm_verdict.csv``, ``data/leipzig/freq/ron_news_2024_1M_freq.csv``,
the legacy TP summaries and ``analysis/learnability_checks/check_claims_output.txt``).
The legacy summaries of the fixture are written from the reproduction of the R
definitions (``parity.legacy_tables``) in the R script's CSV format, so the parity
stage has inputs; their agreement with the real R output is tested on the real data
(``needs_data``).
"""

from __future__ import annotations

import csv
import shutil
import sys
from collections.abc import Mapping
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lexpipe import tp  # noqa: E402
from lexpipe.langs.romanian_compare import lexicon as lx  # noqa: E402
from lexpipe.langs.romanian_compare import parity  # noqa: E402

FIXTURES = Path(__file__).resolve().parent / "fixtures"
REAL_LEXICON = PROJECT_ROOT / "data" / "romanian_lexicon_with_freq.csv"
REAL_VERDICTS = PROJECT_ROOT / "data" / "verb_paradigm_verdict.csv"
REAL_FREQ = (
    PROJECT_ROOT / "data" / "leipzig" / "freq" / "ron_news_2024_1M_freq.csv"
)
REAL_SUMMARY_FULL = PROJECT_ROOT / "analysis" / "romanian_tp_summary_full.csv"
REAL_SUMMARY_DS = (
    PROJECT_ROOT / "analysis" / "romanian_tp_summary_downsampled.csv"
)
# The curated -log decisions (tracked; read-only to the pipeline).
DECISIONS = (
    PROJECT_ROOT / "curated" / "romanian_compare" / "suffixoid_classes.csv"
)
CHECK_CLAIMS = """\
========================================================================
PREDICTION 1: plural-less nouns resist palatalization in derivatives
========================================================================
candidates (no plural in data, has derived verb): 1
========================================================================
TOP-1000 LEARNER SAMPLE (same construction as the R script)
"""


def _flag(value: bool | None) -> str:
    return "NA" if value is None else ("TRUE" if value else "FALSE")


def write_legacy_csv(
    path: Path, counts: Mapping[str, tuple[int, int]]
) -> None:
    """Rows in the R script's summary format (labels sorted)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(tp.TP_FIELDS)
        writer.writerow(["Dorsals", "NA", "NA", "NA", "NA", "NA", "NA", ""])
        for label in sorted(counts):
            mutated, non = counts[label]
            row = tp.tp_row(label, mutated + non, mutated)
            rate = row["rate"]
            writer.writerow(
                [
                    label,
                    row["N"],
                    row["mutated"],
                    row["non-mutated"],
                    "NA" if rate is None else repr(rate),
                    _flag(row["majority?"]),
                    _flag(row["tolerated?"]),
                    "",
                ]
            )


def write_decisions(root: Path, lemmas: set[str]) -> Path:
    """The tracked curated -log decisions restricted to *lemmas* (a
    decision about a lemma the tree does not have would stop the run:
    its guard finds no lexeme)."""
    with DECISIONS.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        fields = list(reader.fieldnames or [])
        rows = [row for row in reader if row["lemma"] in lemmas]
    target = root / "curated" / "romanian_compare" / DECISIONS.name
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return target


def build_tree(root: Path, lexicon: Path = FIXTURES / "lexicon.csv") -> Path:
    """Lay the fixture inputs out under *root*; return *root*."""
    data = root / "data"
    (data / "leipzig" / "freq").mkdir(parents=True, exist_ok=True)
    shutil.copyfile(lexicon, data / "romanian_lexicon_with_freq.csv")
    shutil.copyfile(
        FIXTURES / "verb_paradigm_verdict.csv",
        data / "verb_paradigm_verdict.csv",
    )
    shutil.copyfile(
        FIXTURES / "freq.csv",
        data / "leipzig" / "freq" / "ron_news_2024_1M_freq.csv",
    )
    checks = root / "analysis" / "learnability_checks"
    checks.mkdir(parents=True, exist_ok=True)
    (checks / "check_claims_output.txt").write_text(
        CHECK_CLAIMS, encoding="utf-8"
    )
    rows = lx.read_lexicon(data / "romanian_lexicon_with_freq.csv")
    write_decisions(root, {row.lemma for row in rows if row.pos == "N"})
    full, small = parity.legacy_tables(rows)
    write_legacy_csv(root / "analysis" / "romanian_tp_summary_full.csv", full)
    write_legacy_csv(
        root / "analysis" / "romanian_tp_summary_downsampled.csv", small
    )
    return root


def real_data_or_skip() -> None:
    missing = [
        p
        for p in (
            REAL_LEXICON,
            REAL_VERDICTS,
            REAL_FREQ,
            REAL_SUMMARY_FULL,
            REAL_SUMMARY_DS,
        )
        if not p.is_file()
    ]
    if missing:
        pytest.skip(f"gitignored data absent: {missing[0]}")
