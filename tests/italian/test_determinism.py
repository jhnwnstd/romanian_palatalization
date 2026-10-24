"""Determinism of every Italian stage (design rule R5).

The fixture pipeline runs three times in fresh interpreters:

  A  PYTHONHASHSEED=0, inputs in file order;
  B  PYTHONHASHSEED=1, inputs in file order;
  C  PYTHONHASHSEED=1, the lines of every input shuffled (en split,
     it.wiktionary extract, Morph-it!, the words table).

A and B must give byte-identical outputs. C must give the same outputs
except for the source stamps, which record the SHA-256 of the (now
different) input files: the ``_meta`` line of the projections, the
``#`` line of the form index and the input-hash lines of the QC report
are left out of that comparison.
"""

from __future__ import annotations

import gzip
import hashlib
import os
import re
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path

import pytest

from lexpipe.langs.italian import attestation, extract, features, lexicon, qc
from lexpipe.langs.italian.layout import Layout
from tests.italian import helpers

pytestmark = [pytest.mark.it, pytest.mark.slow]

STAGE_OUTPUTS: dict[str, Callable[[Layout], list[Path]]] = {
    "extract": extract.outputs,
    "lexicon": lexicon.outputs,
    "qc": qc.outputs,
    "attest": attestation.outputs,
    "features": features.outputs,
}
_SHA256_LINE = re.compile(r"[0-9a-f]{64}")


def _run(root: Path, hash_seed: str, shuffle: int | None) -> Layout:
    env = dict(os.environ, PYTHONHASHSEED=hash_seed)
    command = [sys.executable, str(helpers.DRIVER), str(root)]
    if shuffle is not None:
        command.append(str(shuffle))
    result = subprocess.run(
        command,
        cwd=helpers.PROJECT_ROOT,
        env=env,
        capture_output=True,
        text=True,
        timeout=600,
        check=False,
    )
    assert result.returncode == 0, result.stderr[-4000:]
    return Layout(helpers.make_context(root))


@pytest.fixture(scope="module")
def runs(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Layout]:
    base = tmp_path_factory.mktemp("determinism")
    return {
        "A": _run(base / "a", "0", None),
        "B": _run(base / "b", "1", None),
        "C": _run(base / "c", "1", 7),
    }


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _unstamped(path: Path) -> bytes:
    """File content without its source stamps (module docstring)."""
    if path.name.endswith(".gz"):
        lines = gzip.decompress(path.read_bytes()).splitlines(keepends=True)
        return b"".join(lines[1:])
    data = path.read_bytes()
    if path.suffix == ".txt":
        kept = [
            line
            for line in data.splitlines(keepends=True)
            if not _SHA256_LINE.search(line.decode("utf-8"))
        ]
        return b"".join(kept)
    return data


@pytest.mark.parametrize("stage", list(STAGE_OUTPUTS))
def test_stage_outputs_do_not_depend_on_hash_seed_or_input_order(
    stage: str, runs: dict[str, Layout]
) -> None:
    outputs = {
        name: STAGE_OUTPUTS[stage](layout) for name, layout in runs.items()
    }
    for a, b, c in zip(outputs["A"], outputs["B"], outputs["C"], strict=True):
        assert a.name == b.name == c.name
        assert _digest(a.read_bytes()) == _digest(b.read_bytes()), a.name
        assert _digest(_unstamped(a)) == _digest(_unstamped(c)), a.name


def test_shuffled_inputs_really_differ(runs: dict[str, Layout]) -> None:
    """Negative control: the shuffled run did read reordered inputs."""
    a = runs["A"].en_split.read_bytes()
    c = runs["C"].en_split.read_bytes()
    assert a != c and sorted(a.splitlines()) == sorted(c.splitlines())
