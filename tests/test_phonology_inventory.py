#!/usr/bin/env python3
"""Tests for inventory files: lookup, loading, format and git tracking.

Pins:

- **Lookup.** :func:`inventory_path` resolves ``<name>_features.json``
  in the package data directory ``phonology/inventories/``, falls back
  to the repository root for ``romanian`` only, and with ``directory=``
  looks in that one directory. :func:`load_inventory_file` loads what
  it resolves, with patches and underspecified segments.
- **Romanian unchanged.** The Romanian analysis still loads the
  repository-root ``romanian_features.json``.
- **Format.** Every inventory follows the format documented in
  ``src/phonology/inventories/README.md``.
- **Tracking.** ``*_features.json`` files in the inventory directory
  are not gitignored (skipped outside a git checkout).

Temporary inventories live in pytest's tmp directories; the repository
is never written to.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
import unicodedata
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import phonology.analyses.romanian_palatalization as romanian  # noqa: E402
import phonology.inventory as inventory_module  # noqa: E402
from phonology.inventory import (  # noqa: E402
    INVENTORY_DIR,
    FeatureInventory,
    FeaturePatch,
    UnderspecifiedSegment,
    inventory_path,
    load_inventory_file,
)

ROMANIAN_JSON = PROJECT_ROOT / "romanian_features.json"
README = INVENTORY_DIR / "README.md"

TOY_FEATURES = ("Syllabic", "Consonantal", "CORONAL", "DORSAL")
TOY_SEGMENTS = {
    "a": dict(zip(TOY_FEATURES, ("+", "-", "-", "+"), strict=True)),
    "k": dict(zip(TOY_FEATURES, ("-", "+", "-", "+"), strict=True)),
    "t": dict(zip(TOY_FEATURES, ("-", "+", "+", "-"), strict=True)),
}


def _write_toy(directory: Path, name: str = "toy") -> Path:
    path = directory / f"{name}_features.json"
    inventory = {
        "metadata": {"name": "Toy inventory", "notes": "For tests."},
        "features": list(TOY_FEATURES),
        "segments": TOY_SEGMENTS,
    }
    path.write_text(
        json.dumps(inventory, ensure_ascii=False), encoding="utf-8"
    )
    return path


@pytest.fixture
def package_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Stand an empty temporary directory in for INVENTORY_DIR."""
    directory = tmp_path / "inventories"
    directory.mkdir()
    monkeypatch.setattr(inventory_module, "INVENTORY_DIR", directory)
    return directory


# ---------------------------------------------------------------------------
# Lookup and loading
# ---------------------------------------------------------------------------


class TestLookup:
    def test_temporary_inventory_through_directory(
        self, tmp_path: Path
    ) -> None:
        path = _write_toy(tmp_path)
        assert inventory_path("toy", directory=tmp_path) == path
        inventory = load_inventory_file("toy", directory=tmp_path)
        assert inventory.features == TOY_FEATURES
        assert dict(inventory.base) == TOY_SEGMENTS

    def test_load_applies_patches_and_underspec(self, tmp_path: Path) -> None:
        path = _write_toy(tmp_path)
        patches = (
            FeaturePatch(segment="k", updates={"CORONAL": "0"}, reason="t"),
        )
        underspec = (
            UnderspecifiedSegment(
                label="K", base="k", clear=frozenset({"DORSAL"}), reason="t"
            ),
        )
        inventory = load_inventory_file(
            "toy", patches, underspec, directory=tmp_path
        )
        assert inventory.base["k"]["CORONAL"] == "0"
        assert inventory.features_of("K")["DORSAL"] == "0"
        assert inventory == FeatureInventory.load(path, patches, underspec)

    def test_directory_is_the_only_place_looked(self, tmp_path: Path) -> None:
        # No fallback to the repository root, even for romanian.
        with pytest.raises(FileNotFoundError) as excinfo:
            inventory_path("romanian", directory=tmp_path)
        assert str(tmp_path / "romanian_features.json") in str(excinfo.value)
        assert str(ROMANIAN_JSON) not in str(excinfo.value)

    def test_default_lookup_reads_the_package_directory(
        self, package_dir: Path
    ) -> None:
        path = _write_toy(package_dir)
        assert inventory_path("toy") == path
        assert load_inventory_file("toy").features == TOY_FEATURES

    def test_romanian_falls_back_to_the_repository_root(
        self, package_dir: Path
    ) -> None:
        assert inventory_path("romanian") == ROMANIAN_JSON

    def test_package_directory_comes_first(self, package_dir: Path) -> None:
        copy = package_dir / "romanian_features.json"
        shutil.copyfile(ROMANIAN_JSON, copy)
        assert inventory_path("romanian") == copy

    def test_only_romanian_has_a_fallback(self, package_dir: Path) -> None:
        with pytest.raises(FileNotFoundError) as excinfo:
            inventory_path("italian")
        message = str(excinfo.value)
        assert str(package_dir / "italian_features.json") in message
        assert str(PROJECT_ROOT / "italian_features.json") not in message

    @pytest.mark.parametrize(
        "name",
        [
            "",
            "Romanian",
            "../romanian",
            "czech/../romanian",
            "romanian_features.json",
            "it-IT",
            "1czech",
        ],
    )
    def test_rejects_names_that_are_not_identifiers(self, name: str) -> None:
        with pytest.raises(ValueError):
            inventory_path(name)


# ---------------------------------------------------------------------------
# The Romanian analysis
# ---------------------------------------------------------------------------


class TestRomanianUnchanged:
    def test_analysis_uses_the_repository_root_file(self) -> None:
        assert inventory_path("romanian") == ROMANIAN_JSON
        assert romanian.build_profile().inventory_json == ROMANIAN_JSON

    def test_loader_builds_the_analysis_inventory(self) -> None:
        loaded = load_inventory_file(
            "romanian", romanian.PATCHES, romanian.UNDERSPEC
        )
        assert loaded == romanian.load_inventory()
        assert loaded == FeatureInventory.load(
            ROMANIAN_JSON, romanian.PATCHES, romanian.UNDERSPEC
        )


# ---------------------------------------------------------------------------
# The package data directory and the file format
# ---------------------------------------------------------------------------


def _inventory_files() -> list[Path]:
    """The reference inventory plus every inventory in INVENTORY_DIR."""
    return [ROMANIAN_JSON, *sorted(INVENTORY_DIR.glob("*_features.json"))]


def _reference_features() -> list[str]:
    raw = json.loads(ROMANIAN_JSON.read_text(encoding="utf-8"))
    features: list[str] = raw["features"]
    return features


class TestPackageDirectory:
    def test_directory_sits_in_the_package_with_its_readme(self) -> None:
        assert (
            INVENTORY_DIR == PROJECT_ROOT / "src" / "phonology" / "inventories"
        )
        assert README.is_file()

    def test_json_files_carry_the_tracked_suffix(self) -> None:
        for path in INVENTORY_DIR.glob("*.json"):
            assert path.name.endswith("_features.json"), path.name

    def test_readme_lists_the_reference_features_in_order(self) -> None:
        cells = re.findall(
            r"\|\s*(\d+)\s*\|\s*`(\w+)`", README.read_text("utf-8")
        )
        listed = [name for _, name in sorted(cells, key=lambda c: int(c[0]))]
        assert listed == _reference_features()
        assert len(listed) == 27 == len(set(listed))


class TestFormat:
    @pytest.mark.parametrize("path", _inventory_files(), ids=lambda p: p.name)
    def test_follows_the_documented_format(self, path: Path) -> None:
        raw = json.loads(path.read_text(encoding="utf-8"))
        features = _reference_features()
        assert list(raw) == ["metadata", "features", "segments"]
        assert isinstance(raw["metadata"]["name"], str)
        assert raw["metadata"]["name"]
        assert isinstance(raw["metadata"]["notes"], str)
        assert raw["features"] == features
        assert raw["segments"]
        for label, values in raw["segments"].items():
            assert label and unicodedata.is_normalized("NFC", label), label
            assert list(values) == features, label
            assert set(values.values()) <= {"+", "-", "0"}, label


# ---------------------------------------------------------------------------
# Git tracking
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def git() -> str:
    exe = shutil.which("git")
    if exe is None or not (PROJECT_ROOT / ".git").exists():
        pytest.skip("not a git checkout")
    return exe


def _ignored(git: str, relpath: str) -> bool:
    """True iff the gitignore rules exclude ``relpath`` (read-only)."""
    result = subprocess.run(
        [git, "check-ignore", "--quiet", "--no-index", relpath],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode in (0, 1), result.stderr
    return result.returncode == 0


class TestGitTracking:
    def test_inventory_files_are_not_ignored(self, git: str) -> None:
        directory = INVENTORY_DIR.relative_to(PROJECT_ROOT).as_posix()
        for name in ("italian", "czech"):
            assert not _ignored(git, f"{directory}/{name}_features.json")
        assert not _ignored(git, f"{directory}/README.md")
        assert not _ignored(git, "romanian_features.json")

    def test_other_json_names_would_be_ignored(self, git: str) -> None:
        directory = INVENTORY_DIR.relative_to(PROJECT_ROOT).as_posix()
        assert _ignored(git, f"{directory}/italian.json")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
