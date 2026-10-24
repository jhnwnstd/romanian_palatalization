"""Determinism of the Italian analysis stages (design rule R5).

The analysis fixture pipeline (data layer + derivations … docs) runs
three times in fresh interpreters:

  A  PYTHONHASHSEED=0, inputs in file order;
  B  PYTHONHASHSEED=1, inputs in file order;
  C  PYTHONHASHSEED=1, the lines of every input shuffled (en split,
     it.wiktionary extract, Morph-it!, the words table).

Every output of the analysis stages must be byte-identical in A, B and
C: none of them records input hashes, and the learners draw from the
run context's seed only.
"""

from __future__ import annotations

import hashlib
import os
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path

import pytest

from lexpipe.langs.italian import (
    covariation,
    derivations,
    docstage,
    logtable,
    noevidence,
    spellout,
    tpstage,
    verbs,
)
from lexpipe.langs.italian.layout import Layout
from tests.italian import analysis_helpers, helpers

pytestmark = [pytest.mark.it, pytest.mark.slow]

STAGE_OUTPUTS: dict[str, Callable[[Layout], list[Path]]] = {
    "derivations": derivations.outputs,
    "verbs": verbs.outputs,
    "tp": tpstage.outputs,
    "covariation": covariation.outputs,
    "spellout": spellout.outputs,
    "noevidence": noevidence.outputs,
    "logtable": logtable.outputs,
    "docs": docstage.outputs,
}


def _run(root: Path, hash_seed: str, shuffle: int | None) -> Layout:
    env = dict(os.environ, PYTHONHASHSEED=hash_seed)
    command = [sys.executable, str(analysis_helpers.DRIVER), str(root)]
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
    return Layout(analysis_helpers.make_context(root))


@pytest.fixture(scope="module")
def runs(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Layout]:
    base = tmp_path_factory.mktemp("analysis_determinism")
    return {
        "A": _run(base / "a", "0", None),
        "B": _run(base / "b", "1", None),
        "C": _run(base / "c", "1", 7),
    }


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.mark.parametrize("stage", list(STAGE_OUTPUTS))
def test_outputs_do_not_depend_on_hash_seed_or_input_order(
    stage: str, runs: dict[str, Layout]
) -> None:
    outputs = {
        name: STAGE_OUTPUTS[stage](layout) for name, layout in runs.items()
    }
    for a, b, c in zip(outputs["A"], outputs["B"], outputs["C"], strict=True):
        assert a.name == b.name == c.name
        assert _digest(a) == _digest(b), a.name
        assert _digest(a) == _digest(c), a.name


def test_shuffled_inputs_really_differ(runs: dict[str, Layout]) -> None:
    """Negative control: the shuffled run read reordered inputs."""
    a = runs["A"].en_split.read_bytes()
    c = runs["C"].en_split.read_bytes()
    assert a != c and sorted(a.splitlines()) == sorted(c.splitlines())


def test_learners_ran(runs: dict[str, Layout]) -> None:
    """The learner rows are not empty (so their determinism is tested)."""
    text = runs["A"].learners.read_text(encoding="utf-8")
    assert len(text.splitlines()) > 1


def test_readme_has_one_generated_block(runs: dict[str, Layout]) -> None:
    readme = runs["C"].readme.read_text(encoding="utf-8")
    assert readme.count(docstage.BEGIN) == readme.count(docstage.END) == 1
