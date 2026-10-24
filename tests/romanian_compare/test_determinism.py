"""Determinism of every stage (R5).

The whole pipeline runs three times on the fixture in fresh
interpreters: with ``PYTHONHASHSEED=0``, with ``PYTHONHASHSEED=1`` and
the learner simulation in one process (``--param workers=1``; the
others use worker processes), and with every input file's data rows
shuffled. Every output must be byte-identical, except two input stamps
that legitimately record the inputs' SHA-256 (the roots JSONL header
line and the parity report's "Inputs" block), which are compared after
removing those stamps.
"""

from __future__ import annotations

import csv
import hashlib
import os
import random
import re
import subprocess
import sys
from pathlib import Path

import pytest

from tests.romanian_compare.helpers import PROJECT_ROOT, build_tree

pytestmark = [pytest.mark.ro, pytest.mark.slow]

STAMPED = {"romanian_compare_roots.jsonl", "romanian_parity.txt"}


def shuffle_rows(path: Path, seed: int) -> None:
    """Shuffle a CSV's data rows (quoted newlines kept in their row)."""
    with path.open("r", encoding="utf-8", newline="") as fh:
        rows = list(csv.reader(fh))
    header, body = rows[:1], rows[1:]
    random.Random(seed).shuffle(body)
    with path.open("w", encoding="utf-8", newline="") as fh:
        csv.writer(fh, lineterminator="\n").writerows(header + body)


def run_pipeline(root: Path, hash_seed: str, *extra: str) -> None:
    env = dict(os.environ)
    env["PYTHONHASHSEED"] = hash_seed
    env["PYTHONPATH"] = str(PROJECT_ROOT / "src")
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "lexpipe.langs.romanian_compare",
            "--root",
            str(root),
            "--learner-samples",
            "40",
            "--learner-V",
            "5,10,20",
            *extra,
        ],
        cwd=root,
        env=env,
        capture_output=True,
        text=True,
        timeout=600,
        check=False,
    )
    assert result.returncode == 0, (
        result.stdout[-3000:] + result.stderr[-3000:]
    )


def unstamped(path: Path) -> bytes:
    text = path.read_text(encoding="utf-8")
    if path.name == "romanian_compare_roots.jsonl":
        return text.split("\n", 1)[1].encode("utf-8")
    if path.name == "romanian_parity.txt":
        return re.sub(r"Inputs \(SHA-256\):\n(  .*\n)+", "", text).encode(
            "utf-8"
        )
    return path.read_bytes()


def outputs(root: Path) -> dict[str, str]:
    found: dict[str, str] = {}
    for base in (
        root / "analysis" / "comparison",
        root / "data" / "romanian_compare",
    ):
        for path in sorted(base.rglob("*")):
            if not path.is_file() or "manifests" in path.parts:
                continue
            data = (
                unstamped(path) if path.name in STAMPED else path.read_bytes()
            )
            found[path.relative_to(root).as_posix()] = hashlib.sha256(
                data
            ).hexdigest()
    return found


def test_outputs_do_not_depend_on_hash_seed_or_input_order(
    tmp_path: Path,
) -> None:
    first = build_tree(tmp_path / "a")
    second = build_tree(tmp_path / "b")
    third = build_tree(tmp_path / "c")
    for index, name in enumerate(
        (
            "data/romanian_lexicon_with_freq.csv",
            "data/verb_paradigm_verdict.csv",
            "data/leipzig/freq/ron_news_2024_1M_freq.csv",
        )
    ):
        shuffle_rows(third / name, seed=index + 1)
    run_pipeline(first, "0", "--param", "workers=3")
    run_pipeline(second, "1", "--param", "workers=1")
    run_pipeline(third, "0")
    reference = outputs(first)
    assert len(reference) == 12
    assert outputs(second) == reference
    assert outputs(third) == reference
