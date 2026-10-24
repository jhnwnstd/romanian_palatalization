"""Determinism of every Czech stage (design rule R5).

The five stages run on the fixture tree in fresh interpreters with
``PYTHONHASHSEED`` 0 and 1, and once more with every input shuffled
(kaikki lines of both editions, MorfFlex rows, Leipzig rows, CNK rows).
Every output of every stage must have the same SHA-256 in the two
hash-seed runs. With shuffled inputs the outputs that are stamped with
their source's SHA-256 (R7: the extract stats, the SQLite derivatives'
``meta`` tables, the QC report) necessarily differ in that stamp, so
they are compared with the stamp left out. Manifests are not compared
(they record input hashes and the git state).
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

from tests.czech import treebuilder

pytestmark = [pytest.mark.cs, pytest.mark.slow]

_DRIVER = """
import json, sys
from pathlib import Path
sys.path.insert(0, sys.argv[1])
from tests.czech import treebuilder
root = treebuilder.build_tree(Path(sys.argv[2]), shuffle_seed=json.loads(sys.argv[3]))
treebuilder.run_stages(root)
"""


def _unstamped(path: Path) -> str:
    """SHA-256 of a file's content without source-hash stamps."""
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
    for line in path.read_bytes().splitlines():
        if b"sha256" not in line:
            digest.update(line + b"\n")
    return digest.hexdigest()


def _run(
    root: Path,
    hash_seed: str,
    shuffle_seed: int | None,
    *,
    stamped: bool = False,
) -> dict[str, str]:
    env = {**os.environ, "PYTHONHASHSEED": hash_seed}
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            _DRIVER,
            str(treebuilder.PROJECT_ROOT),
            str(root),
            json.dumps(shuffle_seed),
        ],
        cwd=treebuilder.PROJECT_ROOT,
        env=env,
        capture_output=True,
        text=True,
        timeout=600,
        check=False,
    )
    assert result.returncode == 0, result.stderr[-4000:]
    digests: dict[str, str] = {}
    for stage, outputs in treebuilder.stage_outputs(root).items():
        for path in outputs:
            digest = (
                _unstamped(path)
                if stamped
                else hashlib.sha256(path.read_bytes()).hexdigest()
            )
            digests[f"{stage}:{path.relative_to(root).as_posix()}"] = digest
    return digests


def test_outputs_identical_across_hash_seeds_and_input_order(
    tmp_path: Path,
) -> None:
    first = _run(tmp_path / "seed0", "0", None)
    second = _run(tmp_path / "seed1", "1", None)
    assert len(first) == sum(
        len(v) for v in treebuilder.stage_outputs(tmp_path / "seed0").values()
    )
    assert [k for k in first if first[k] != second[k]] == []
    base = {
        key: _unstamped(tmp_path / "seed0" / key.split(":", 1)[1])
        for key in first
    }
    shuffled = _run(tmp_path / "shuffled", "1", 20080715, stamped=True)
    assert [k for k in base if base[k] != shuffled[k]] == []
