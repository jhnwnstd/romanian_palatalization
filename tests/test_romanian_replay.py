#!/usr/bin/env python3
"""Romanian regression guard: replay stages 4-5 and the validator.

The multilingual extension shares code and data directories with the
Romanian pipeline, whose outputs must not move. This test replays,
into a temporary directory and each in a fresh interpreter:

  - stage 4, ``scripts/romanian_processor_main.py``:
    ``process_csv(romanian_lexicon_raw_dex.csv → complete.csv)``;
  - stage 5, ``scripts/attach_frequencies.py``: ``main()`` with its
    lexicon and output paths redirected (``complete → with_freq``);

and compares SHA-256 against the baseline pinned below. A fresh
interpreter matters: ``romanian_processor_lib`` calibrates its
plural-plausibility thresholds in module state, so a second
``process_csv`` in one process gives different output. The validator
API (``phonology.validation``, called as
``scripts/validate_palatalization_rules.py`` does) then runs on the
replayed ``with_freq`` CSV, and its report text and mismatch CSV must
equal the tracked ``analysis/palatalization_rule_validation.txt`` and
``analysis/palatalization_mismatches.csv`` byte for byte.

Stage 5 attaches every ``*_freq.csv`` in ``data/leipzig/freq/``, so a
stray table there (an Italian one, say) changes the replay;
:func:`test_romanian_freq_dir_is_unchanged` names it.

The inputs are gitignored: the tests skip when they are absent. If an
input changes on purpose (a new DEX QC run, a new Leipzig release),
re-run the pipeline, review the outputs and re-pin the constants.

Marked ``replay``; it runs by default (deselect with
``-m "not replay"``).
"""

from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from phonology.analyses.romanian_palatalization import (  # noqa: E402
    build_profile,
)
from phonology.lexicon import iter_lexicon_rows  # noqa: E402
from phonology.validation import (  # noqa: E402
    render_report,
    run_validation,
    write_mismatch_csv,
)

pytestmark = pytest.mark.replay

RAW_DEX = PROJECT_ROOT / "data" / "romanian_lexicon_raw_dex.csv"
FREQ_DIR = PROJECT_ROOT / "data" / "leipzig" / "freq"
FREQ_CSV = FREQ_DIR / "ron_news_2024_1M_freq.csv"
STAGE4_SCRIPT = PROJECT_ROOT / "scripts" / "romanian_processor_main.py"
STAGE5_SCRIPT = PROJECT_ROOT / "scripts" / "attach_frequencies.py"
REPORT_PATH = PROJECT_ROOT / "analysis" / "palatalization_rule_validation.txt"
MISMATCH_CSV = PROJECT_ROOT / "analysis" / "palatalization_mismatches.csv"

# Baseline confirmed by replay on 2026-10-23.
RAW_DEX_SHA256 = (
    "f61f4fd6aa6389e75e17d9327351716c0c92e2b98604e091c0c206e09a6fd8fe"
)
FREQ_CSV_SHA256 = (
    "790b8078a525485c0a631ae8bcc315e4a148f48ea6f0309c9dd2eaeb4a7a98e8"
)
COMPLETE_SHA256 = (
    "b8d57673289620b0214bd507524f46f3707bef3d229c53f5264eae06d35e7690"
)
WITH_FREQ_SHA256 = (
    "fb51334926d12eaffae85dd946efcc0b88149f03109e2175542706f50c3250ad"
)
PINNED_FREQ_FILES = ["ron_news_2024_1M_freq.csv"]

_LOAD_SCRIPT = """
import importlib.util
import sys
from pathlib import Path

spec = importlib.util.spec_from_file_location("replayed", sys.argv[1])
module = importlib.util.module_from_spec(spec)
sys.modules["replayed"] = module
spec.loader.exec_module(module)
"""
_STAGE4 = _LOAD_SCRIPT + "module.process_csv(sys.argv[2], sys.argv[3])\n"
_STAGE5 = _LOAD_SCRIPT + (
    "module.LEXICON_PATH = Path(sys.argv[2])\n"
    "module.OUTPUT_PATH = Path(sys.argv[3])\n"
    "sys.exit(module.main())\n"
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def _run_stage(code: str, script: Path, src: Path, dest: Path) -> None:
    """Run one stage in a fresh interpreter, cwd = the output directory."""
    result = subprocess.run(
        [sys.executable, "-c", code, str(script), str(src), str(dest)],
        cwd=dest.parent,
        capture_output=True,
        text=True,
        timeout=1800,
        check=False,
    )
    assert result.returncode == 0, (
        f"{script.name} failed ({result.returncode}):\n"
        f"{result.stdout[-2000:]}\n{result.stderr[-4000:]}"
    )


def _repin_hint(what: str) -> str:
    return (
        f"{what} changed. If intended, review the new outputs and re-pin "
        "the constants in tests/test_romanian_replay.py."
    )


@pytest.fixture(scope="module")
def inputs() -> None:
    for path in (RAW_DEX, FREQ_CSV):
        if not path.exists():
            pytest.skip(f"gitignored replay input absent: {path}")


@pytest.fixture(scope="module")
def replay(
    inputs: None, tmp_path_factory: pytest.TempPathFactory
) -> tuple[Path, Path]:
    out = tmp_path_factory.mktemp("romanian_replay")
    complete = out / "romanian_lexicon_complete.csv"
    with_freq = out / "romanian_lexicon_with_freq.csv"
    _run_stage(_STAGE4, STAGE4_SCRIPT, RAW_DEX, complete)
    _run_stage(_STAGE5, STAGE5_SCRIPT, complete, with_freq)
    return complete, with_freq


@pytest.fixture(scope="module")
def validation(replay: tuple[Path, Path]) -> tuple[str, Path]:
    _, with_freq = replay
    profile = build_profile()
    rows = iter_lexicon_rows(with_freq, target_stems=profile.target_stems)
    report = run_validation(profile, rows)
    text = render_report(report, include_distance_summary=True)
    mismatch_csv = with_freq.parent / "palatalization_mismatches.csv"
    write_mismatch_csv(report, mismatch_csv)
    return text, mismatch_csv


def test_inputs_match_pinned_baseline(inputs: None) -> None:
    assert _sha256(RAW_DEX) == RAW_DEX_SHA256, _repin_hint(RAW_DEX.name)
    assert _sha256(FREQ_CSV) == FREQ_CSV_SHA256, _repin_hint(FREQ_CSV.name)


def test_romanian_freq_dir_is_unchanged(inputs: None) -> None:
    found = sorted(p.name for p in FREQ_DIR.glob("*_freq.csv"))
    assert found == PINNED_FREQ_FILES, (
        f"data/leipzig/freq/ holds {found}: attach_frequencies.py adds a "
        "column per file there. Non-Romanian tables belong in "
        "data/<language>/leipzig/freq/."
    )


def test_stage4_replay_matches_baseline(replay: tuple[Path, Path]) -> None:
    complete, _ = replay
    assert _sha256(complete) == COMPLETE_SHA256, _repin_hint(
        "Stage 4 output (romanian_lexicon_complete.csv)"
    )


def test_stage5_replay_matches_baseline(replay: tuple[Path, Path]) -> None:
    _, with_freq = replay
    assert _sha256(with_freq) == WITH_FREQ_SHA256, _repin_hint(
        "Stage 5 output (romanian_lexicon_with_freq.csv)"
    )


def test_validator_report_matches_tracked_file(
    validation: tuple[str, Path],
) -> None:
    text, _ = validation
    replayed = hashlib.sha256(text.encode("utf-8")).hexdigest()
    assert replayed == _sha256(REPORT_PATH), (
        "render_report() on the replayed lexicon differs from "
        f"{REPORT_PATH.relative_to(PROJECT_ROOT)}"
    )


def test_mismatch_csv_matches_tracked_file(
    validation: tuple[str, Path],
) -> None:
    _, mismatch_csv = validation
    assert _sha256(mismatch_csv) == _sha256(MISMATCH_CSV), (
        "write_mismatch_csv() on the replayed lexicon differs from "
        f"{MISMATCH_CSV.relative_to(PROJECT_ROOT)}"
    )
