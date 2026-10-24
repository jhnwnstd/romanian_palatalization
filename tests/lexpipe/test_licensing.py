"""Tests for lexpipe.licensing and the csvio.write_tracked gate (R9).

The last test scans every CSV git would commit (tracked, or untracked
and not ignored) and fails if a ``source`` or ``licence_class`` column
holds a value whose data may not be redistributed in the repository.
"""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lexpipe import csvio, licensing, sources, spec  # noqa: E402
from lexpipe.licensing import LicenceClass, LicenceError  # noqa: E402
from lexpipe.spec import LanguageInfo  # noqa: E402

REDISTRIBUTABLE = {
    "wiktionary_en": True,
    "wiktionary_it": True,
    "wiktionary_cs": True,
    "wiktionary_ro": True,
    "morphit": True,
    "cnk": True,
    "morfflex": False,
    "derinet": False,
    "czetyl": False,
    "leipzig": False,
    "dex_qc": False,
    "treccani": True,
    "ijp": True,
}
FIELDS = ["lexeme_id", "plural", "source"]


def test_registry() -> None:
    assert {
        name: entry.redistributable
        for name, entry in licensing.LICENCES.items()
    } == REDISTRIBUTABLE
    assert licensing.licence_of("morfflex").licence == "CC BY-NC-SA 4.0"
    assert licensing.licence_of("morphit").licence == "CC BY-SA 2.0 / LGPL"
    assert licensing.licence_of("cnk").licence_class is LicenceClass.BY
    assert licensing.licence_of("leipzig").licence == "terms unverified"
    for name in ("treccani", "ijp"):
        entry = licensing.licence_of(name)
        assert entry.licence_class is LicenceClass.CITATION
        assert "citation" in entry.licence


def test_every_download_has_a_row_source() -> None:
    assert set(licensing.DOWNLOADS) == set(sources.SOURCES)
    for name, code in licensing.DOWNLOADS.items():
        entry = licensing.licence_of(code)
        download = sources.get_source(name)
        non_commercial = "NC" in download.licence
        unverified = "unverified" in download.licence
        assert entry.redistributable is not (non_commercial or unverified)


def test_unknown_source() -> None:
    with pytest.raises(LicenceError, match="unknown source 'dex'"):
        licensing.licence_of("dex")


def test_split_sources() -> None:
    assert licensing.split_sources("wiktionary_it|morphit") == [
        "wiktionary_it",
        "morphit",
    ]
    assert licensing.split_sources(["cnk", "ijp"]) == ["cnk", "ijp"]
    assert licensing.split_sources("") == []


@pytest.mark.parametrize(
    "sources",
    [["wiktionary_cs", "morfflex"], ["leipzig"], [], [""], ["x"]],
)
def test_check_redistributable_refuses(sources: list[str]) -> None:
    with pytest.raises(LicenceError):
        licensing.check_redistributable(sources)


def test_check_redistributable_accepts() -> None:
    licensing.check_redistributable(["wiktionary_it", "morphit", "treccani"])


def test_licence_classes() -> None:
    licensing.check_licence_class("by-sa|by|citation")
    for bad in ("by-nc-sa", "by-sa|unverified", "public-domain"):
        with pytest.raises(LicenceError):
            licensing.check_licence_class(bad)


@pytest.mark.parametrize(
    ("relative", "expected"),
    [
        ("analysis/italian/tp_long.csv", True),
        ("analysis/czech/sub/dir/x.csv", True),
        ("analysis/comparison/tp_long.csv", True),
        ("docs/czech/data_dict.md", True),
        ("data/italian/italian_lexicon_raw.csv", True),
        ("data/czech/czech_lexicon_raw.csv", True),
        ("data/czech/italian_lexicon_raw.csv", False),
        ("data/italian/italian_plural_forms.csv", False),
        ("data/italian/nc/worklist.csv", False),
        ("analysis/romanian_tp_summary_full.csv", False),  # Romanian, R1
        ("analysis/romanian/x.csv", False),
        ("curated/italian/overrides.csv", False),  # read-only to stages
        ("docs/sources_lock.tsv", False),
        ("x.csv", False),
        # comparison-only (Romanian): analysis/comparison/ only
        ("analysis/romanian_compare/x.csv", False),
        ("docs/romanian_compare/data_dict.md", False),
        ("data/romanian_compare/romanian_compare_lexicon_raw.csv", False),
    ],
)
def test_tracked_locations(
    tmp_path: Path, relative: str, expected: bool
) -> None:
    path = tmp_path / relative
    assert licensing.is_tracked_location(path, tmp_path) is expected
    assert licensing.is_tracked_location(Path("/elsewhere/x.csv")) is False


def test_romanian_sources() -> None:
    ro = licensing.licence_of("wiktionary_ro")
    assert ro.licence == "CC BY-SA 4.0 + GFDL"
    assert ro.licence_class is LicenceClass.BY_SA and ro.redistributable
    dex = licensing.licence_of("dex_qc")
    assert dex.licence_class is LicenceClass.UNVERIFIED
    assert not dex.redistributable and "dexonline" in dex.licence
    with pytest.raises(LicenceError, match="dex_qc"):
        licensing.check_redistributable(["wiktionary_ro", "dex_qc"])


def test_derived_subsets_have_the_dump_licence() -> None:
    for name in ("kaikki_raw_it", "kaikki_raw_cs"):
        assert licensing.DOWNLOADS[name] == licensing.DOWNLOADS["kaikki_raw"]


def test_aggregate_sources() -> None:
    entries = licensing.check_aggregate_sources(
        ["wiktionary_cs", "morfflex", "derinet", "morfflex", "ijp"]
    )
    assert [e.source for e in entries] == [
        "derinet",
        "ijp",
        "morfflex",
        "wiktionary_cs",
    ]
    assert licensing.AGGREGATE_CLASSES == (
        set(LicenceClass) - {LicenceClass.UNVERIFIED}
    )
    for bad in (["leipzig"], ["morfflex", "dex_qc"], ["nope"], [], [""]):
        with pytest.raises(LicenceError):
            licensing.check_aggregate_sources(bad)


def test_aggregate_header() -> None:
    header = licensing.aggregate_header(["morfflex", "wiktionary_cs"])
    lines = header.splitlines()
    assert lines[:2] == [
        licensing.AGGREGATE_HEADER_TITLE,
        "-" * len(licensing.AGGREGATE_HEADER_TITLE),
    ]
    assert header.endswith("\n\n")
    assert "Aggregate counts only" in header
    listed = [line for line in lines if line.startswith("  ")]
    assert [line.split()[0] for line in listed] == [
        "morfflex",
        "wiktionary_cs",
    ]
    assert "CC BY-NC-SA 4.0" in listed[0]
    assert "DATA_LICENSE.md" in header and "commercially" in header
    # Without an NC source there is no NC paragraph; order is canonical.
    plain = licensing.aggregate_header(["wiktionary_it", "cnk"])
    assert "DATA_LICENSE.md" not in plain
    assert plain == licensing.aggregate_header(["cnk", "wiktionary_it"])
    with pytest.raises(LicenceError):
        licensing.aggregate_header(["leipzig"])


def test_with_language() -> None:
    toy = LanguageInfo("zz", "toy", "m")
    merged = licensing.with_language(toy)
    assert merged["zz"] is toy and set(merged) == {"it", "cs", "ro", "zz"}
    # A registered code keeps its registry row.
    unflagged = LanguageInfo("ro", "romanian_compare", "m")
    assert licensing.with_language(unflagged)["ro"].comparison_only


def test_comparison_only_languages_add_no_tracked_location(
    tmp_path: Path,
) -> None:
    registry = {
        "it": spec.LANGUAGES["it"],
        "zz": LanguageInfo(
            "zz", "toy_compare", "m", file_prefix="toy_", comparison_only=True
        ),
    }
    for relative, expected in (
        ("analysis/toy_compare/x.csv", False),
        ("docs/toy_compare/x.md", False),
        ("data/toy_compare/toy_lexicon_raw.csv", False),
        ("analysis/comparison/x.csv", True),
        ("analysis/italian/x.csv", True),
    ):
        path = tmp_path / relative
        assert licensing.is_tracked_location(path, tmp_path, registry) is (
            expected
        )


# ---------------------------------------------------------------------------
# csvio.write_tracked
# ---------------------------------------------------------------------------


def _rows(*sources: Any) -> list[dict[str, Any]]:
    return [
        {"lexeme_id": f"it:x{i}:noun:0:00000000", "plural": "x", "source": s}
        for i, s in enumerate(sources)
    ]


def test_write_tracked_writes_redistributable_rows(tmp_path: Path) -> None:
    path = tmp_path / "analysis" / "italian" / "plurals.csv"
    rows = _rows("wiktionary_en", "wiktionary_en|morphit", ["cnk", "ijp"])
    assert csvio.write_tracked(path, rows, FIELDS, root=tmp_path) == 3
    written = csvio.read_rows(path)
    assert [row["source"] for row in written] == [
        "wiktionary_en",
        "wiktionary_en|morphit",
        "cnk|ijp",
    ]


@pytest.mark.parametrize(
    "bad",
    ["morfflex", "wiktionary_cs|derinet", ["czetyl"], "leipzig", "", None],
)
def test_write_tracked_refuses_rows(tmp_path: Path, bad: Any) -> None:
    path = tmp_path / "analysis" / "czech" / "plurals.csv"
    csvio.write_tracked(path, _rows("wiktionary_cs"), FIELDS, root=tmp_path)
    before = path.read_bytes()
    with pytest.raises(LicenceError, match="row 2"):
        csvio.write_tracked(
            path, _rows("wiktionary_cs", bad), FIELDS, root=tmp_path
        )
    assert path.read_bytes() == before
    assert not list(path.parent.glob("*.tmp"))


def test_write_tracked_refuses_locations(tmp_path: Path) -> None:
    rows = _rows("wiktionary_it")
    for relative in ("data/italian/nc/x.csv", "analysis/x.csv", "x.csv"):
        with pytest.raises(LicenceError, match="not a tracked location"):
            csvio.write_tracked(
                tmp_path / relative, rows, FIELDS, root=tmp_path
            )
        assert not (tmp_path / relative).exists()


def test_write_tracked_needs_a_source_column(tmp_path: Path) -> None:
    path = tmp_path / "analysis" / "italian" / "x.csv"
    with pytest.raises(LicenceError, match="'source' column"):
        csvio.write_tracked(path, [{"a": 1}], ["a"], root=tmp_path)
    rows = [{"a": 1, "origin": "morphit"}]
    assert (
        csvio.write_tracked(
            path, rows, ["a", "origin"], source_field="origin", root=tmp_path
        )
        == 1
    )


def test_write_tracked_checks_licence_class(tmp_path: Path) -> None:
    path = tmp_path / "analysis" / "italian" / "x.csv"
    fields = [*FIELDS, "licence_class"]
    rows = _rows("wiktionary_it")
    rows[0]["licence_class"] = LicenceClass.BY_SA
    assert csvio.write_tracked(path, rows, fields, root=tmp_path) == 1
    rows[0]["licence_class"] = "by-nc-sa"
    with pytest.raises(LicenceError, match="by-nc-sa"):
        csvio.write_tracked(path, rows, fields, root=tmp_path)


# ---------------------------------------------------------------------------
# What git would commit
# ---------------------------------------------------------------------------


def _committable_csvs() -> list[Path]:
    try:
        result = subprocess.run(
            [
                "git",
                "-C",
                str(PROJECT_ROOT),
                "ls-files",
                "-z",
                "--cached",
                "--others",
                "--exclude-standard",
                "--",
                "*.csv",
            ],
            capture_output=True,
            timeout=120,
            check=False,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        pytest.skip(f"git unavailable: {exc}")
    if result.returncode != 0:
        pytest.skip("not a git checkout")
    names = result.stdout.decode("utf-8").split("\0")
    return [PROJECT_ROOT / name for name in names if name]


def _violations(path: Path) -> list[str]:
    found: list[str] = []
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        columns = [
            c
            for c in ("source", "licence_class")
            if c in (reader.fieldnames or ())
        ]
        if not columns:
            return found
        for row in reader:
            for item in licensing.split_sources(row.get("source") or ""):
                entry = licensing.LICENCES.get(item)
                if entry is not None and not entry.redistributable:
                    found.append(f"{path.name}:{reader.line_num}: {item}")
            for item in licensing.split_sources(
                row.get("licence_class") or ""
            ):
                if item in {
                    c.value for c in licensing.NON_REDISTRIBUTABLE_CLASSES
                }:
                    found.append(f"{path.name}:{reader.line_num}: {item}")
    return found


def test_committable_csvs_hold_no_restricted_data() -> None:
    csvs = [path for path in _committable_csvs() if path.is_file()]
    assert csvs, "git lists no CSV at all"
    violations = [v for path in csvs for v in _violations(path)]
    assert not violations, violations[:20]
