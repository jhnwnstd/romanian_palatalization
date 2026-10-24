"""Shared pytest fixtures of the Romanian comparison tests."""

from __future__ import annotations

from pathlib import Path

import pytest

from tests.romanian_compare.helpers import build_tree


@pytest.fixture
def tree(tmp_path: Path) -> Path:
    """A temporary repository tree with the fixture inputs."""
    return build_tree(tmp_path / "repo")
