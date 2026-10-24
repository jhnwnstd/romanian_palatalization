"""Determinism of the Czech analysis stages (design rule R5).

The data-layer and analysis stages run on the analysis fixture tree in
fresh interpreters with ``PYTHONHASHSEED`` 0 and 1, and once more with
every input shuffled (kaikki lines, the synthetic DeriNet trees, czetyl
rows, Leipzig rows, MorfFlex and CNK rows). Every output of the
analysis stages must have the same SHA-256 across hash seeds; with
shuffled inputs, the outputs stamped with a source SHA-256 (the DeriNet
projection's ``meta`` table, R7) are compared with the stamp left out.
Learners use the seed of the run context, so their shares repeat too.
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import os
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

from lexpipe.langs.czech import analysis_layout
from tests.czech import analysis_tree

pytestmark = [pytest.mark.cs, pytest.mark.slow]

_DRIVER = """
import json, sys
from pathlib import Path
sys.path.insert(0, sys.argv[1])
from tests.czech import analysis_tree
root = analysis_tree.build_tree(
    Path(sys.argv[2]), shuffle_seed=json.loads(sys.argv[3])
)
analysis_tree.run_all(root)
"""


def outputs(root: Path) -> list[Path]:
    """Every file the analysis stages write."""
    paths = analysis_layout.analysis_paths(analysis_tree.context(root))
    return [
        paths.derinet_db,
        paths.derivations,
        paths.lexeme_strata,
        paths.derivations_summary,
        paths.derivations_report,
        paths.tp_long,
        paths.learners,
        paths.tp_summary,
        paths.tp_roots,
        paths.spellout,
        paths.spellout_roots,
        paths.noevidence,
        paths.noevidence_pairs,
        paths.covariation,
        paths.covariation_units,
        paths.logtable,
        paths.logtable_full,
        paths.logtable_summary,
        paths.data_dict,
    ]


def unstamped(path: Path) -> str:
    """SHA-256 of a file without its source-hash stamps."""
    digest = hashlib.sha256()
    if path.suffix == ".sqlite":
        with contextlib.closing(sqlite3.connect(path)) as con:
            tables = [
                name
                for (name,) in con.execute(
                    "SELECT name FROM sqlite_master WHERE type = 'table' "
                    "ORDER BY name"
                )
            ]
            for table in tables:
                for row in con.execute(
                    f"SELECT * FROM {table} ORDER BY rowid"
                ):
                    if table == "meta" and "sha256" in str(row[0]):
                        continue
                    digest.update(repr(row).encode())
        return digest.hexdigest()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def _run(
    root: Path, hash_seed: str, shuffle_seed: int | None
) -> dict[str, str]:
    env = {**os.environ, "PYTHONHASHSEED": hash_seed}
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            _DRIVER,
            str(analysis_tree.PROJECT_ROOT),
            str(root),
            json.dumps(shuffle_seed),
        ],
        cwd=analysis_tree.PROJECT_ROOT,
        env=env,
        capture_output=True,
        text=True,
        timeout=900,
        check=False,
    )
    assert result.returncode == 0, result.stderr[-4000:]
    return {
        path.relative_to(root).as_posix(): (
            unstamped(path)
            if shuffle_seed is not None
            else hashlib.sha256(path.read_bytes()).hexdigest()
        )
        for path in outputs(root)
    }


def test_analysis_outputs_identical_across_seeds_and_input_order(
    tmp_path: Path,
) -> None:
    first = _run(tmp_path / "seed0", "0", None)
    second = _run(tmp_path / "seed1", "1", None)
    assert len(first) == len(outputs(tmp_path / "seed0"))
    assert [k for k in first if first[k] != second[k]] == []
    base = {key: unstamped(tmp_path / "seed0" / key) for key in first}
    shuffled = _run(tmp_path / "shuffled", "1", 20080715)
    assert [k for k in base if base[k] != shuffled[k]] == []
