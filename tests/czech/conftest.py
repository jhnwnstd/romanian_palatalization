"""Shared fixtures of the Czech tests: one pipeline run on the fixtures."""

from __future__ import annotations

import sys
from collections.abc import Iterator
from pathlib import Path

import pytest

# Test modules import lexpipe directly: put src/ on sys.path first, as
# the other test packages do.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from lexpipe import csvio  # noqa: E402
from lexpipe.langs.czech import layout  # noqa: E402
from tests.czech import analysis_tree, treebuilder  # noqa: E402


@pytest.fixture(scope="session")
def fixture_root(tmp_path_factory: pytest.TempPathFactory) -> Iterator[Path]:
    """A temp tree on which all five stages have run once."""
    root = treebuilder.build_tree(tmp_path_factory.mktemp("czech_tree"))
    treebuilder.run_stages(root)
    yield root


@pytest.fixture(scope="session")
def paths(fixture_root: Path) -> layout.CzechPaths:
    return layout.paths_for(treebuilder.context(fixture_root))


@pytest.fixture(scope="session")
def lexicon_rows(paths: layout.CzechPaths) -> list[dict[str, str]]:
    """Every lexeme with every column (the gitignored full table)."""
    return csvio.read_rows(paths.lexicon_full)


@pytest.fixture(scope="session")
def raw_rows(paths: layout.CzechPaths) -> list[dict[str, str]]:
    """The tracked, lean raw lexicon."""
    return csvio.read_rows(paths.lexicon_raw)


@pytest.fixture(scope="session")
def complete_rows(paths: layout.CzechPaths) -> list[dict[str, str]]:
    return csvio.read_rows(paths.lexicon_complete)


@pytest.fixture(scope="session")
def cell_rows(paths: layout.CzechPaths) -> list[dict[str, str]]:
    return csvio.read_rows(paths.cells)


@pytest.fixture(scope="session")
def analysis_root(tmp_path_factory: pytest.TempPathFactory) -> Iterator[Path]:
    """A temp tree on which the data-layer and analysis stages have run."""
    root = analysis_tree.build_tree(tmp_path_factory.mktemp("czech_analysis"))
    analysis_tree.run_all(root)
    yield root
