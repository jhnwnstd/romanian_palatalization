"""Tests for lexpipe.csvio output conventions."""

from __future__ import annotations

import sys
from enum import StrEnum
from pathlib import Path
from typing import Any

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lexpipe import csvio  # noqa: E402

FIELDS = ["lemma", "mutation", "plurals", "rate", "note"]


class Label(StrEnum):
    UNDERGOER = "undergoer"


def test_file_layout(tmp_path: Path) -> None:
    path = tmp_path / "out" / "lexicon.csv"
    rows: list[dict[str, Any]] = [
        {
            "lemma": "chirurgo",
            "mutation": None,
            "plurals": ["chirurghi", "chirurgi"],
            "rate": 0.5,
            "note": 'variable, "both"',
        },
        {"lemma": "amico", "mutation": True, "plurals": ("amici",)},
        {"lemma": "bosco", "mutation": False, "plurals": [], "rate": 1.0},
        {"lemma": "x", "rate": float("nan"), "note": Label.UNDERGOER},
    ]
    assert csvio.write_rows(path, rows, FIELDS) == 4
    assert path.read_bytes().decode("utf-8") == (
        "lemma,mutation,plurals,rate,note\n"
        'chirurgo,,chirurghi|chirurgi,0.5,"variable, ""both"""\n'
        "amico,True,amici,,\n"
        "bosco,False,,1.0,\n"
        "x,,,,undergoer\n"
    )
    assert not list(path.parent.glob("*.tmp"))


def test_round_trip(tmp_path: Path) -> None:
    path = tmp_path / "rt.csv"
    csvio.write_rows(
        path,
        [{"lemma": "biolog", "mutation": True, "plurals": ["biologové"]}],
        FIELDS,
    )
    (row,) = csvio.read_rows(path)
    assert row == {
        "lemma": "biolog",
        "mutation": "True",
        "plurals": "biologové",
        "rate": "",
        "note": "",
    }
    assert csvio.parse_bool(row["mutation"]) is True
    assert csvio.parse_bool(row["rate"]) is None
    assert csvio.split_list(row["plurals"]) == ["biologové"]
    assert csvio.split_list("") == []
    assert csvio.split_list("a|b") == ["a", "b"]


def test_format_value() -> None:
    assert csvio.format_value(True) == "True"
    assert csvio.format_value(False) == "False"
    assert csvio.format_value(0) == "0"
    assert csvio.format_value(1 / 3) == "0.3333333333333333"
    assert csvio.format_value(["a", None, True]) == "a||True"
    # A "|" or backslash inside an item is escaped, never refused.
    assert csvio.format_value(["a|b"]) == "a\\|b"
    assert csvio.format_value(["a\\b", "c"]) == "a\\\\b|c"
    # Plain strings are never escaped.
    assert csvio.format_value("a|b") == "a|b"
    with pytest.raises(TypeError, match="source order"):
        csvio.format_value({"a", "b"})


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("True", True),  # the Romanian and lexpipe dialect (R10)
        ("true", True),
        ("TRUE", True),  # readr
        (" False ", False),
        ("false", False),
        ("FALSE", False),
        ("", None),
    ],
)
def test_parse_bool_ignores_case(text: str, expected: bool | None) -> None:
    assert csvio.parse_bool(text) is expected


@pytest.mark.parametrize("text", ["yes", "1", "T", "nan"])
def test_parse_bool_rejects_other_words(text: str) -> None:
    with pytest.raises(ValueError, match="not a boolean"):
        csvio.parse_bool(text)


def test_unknown_columns(tmp_path: Path) -> None:
    path = tmp_path / "x.csv"
    rows = [{"lemma": "casa", "stray": 1}]
    with pytest.raises(ValueError, match="stray"):
        csvio.write_rows(path, rows, FIELDS)
    assert not path.exists()
    assert csvio.write_rows(path, rows, FIELDS, extra="ignore") == 1
    assert csvio.read_rows(path)[0]["lemma"] == "casa"
    with pytest.raises(ValueError, match="duplicate"):
        csvio.write_rows(path, rows, ["lemma", "lemma"])


def test_failed_write_keeps_previous_file(tmp_path: Path) -> None:
    path = tmp_path / "keep.csv"
    csvio.write_rows(path, [{"lemma": "casa"}], ["lemma"])
    before = path.read_bytes()
    with pytest.raises(ValueError):
        csvio.write_rows(path, [{"lemma": "x"}, {"stray": 1}], ["lemma"])
    assert path.read_bytes() == before
    with pytest.raises(TypeError):
        csvio.write_rows(path, [{"lemma": {"a", "b"}}], ["lemma"])
    assert path.read_bytes() == before
    assert not list(tmp_path.glob("*.tmp"))


def test_read_rows_rejects_overlong_rows(tmp_path: Path) -> None:
    path = tmp_path / "long.csv"
    path.write_text("a,b\n1,2\n1,2,3\n", encoding="utf-8")
    with pytest.raises(ValueError, match="more cells"):
        csvio.read_rows(path)
    path.write_text("a,b\n1\n", encoding="utf-8")
    assert csvio.read_rows(path) == [{"a": "1", "b": ""}]


# ---------------------------------------------------------------------------
# Escaped list items
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "items",
    [
        ["ˈmɛdiko | ˈmɛditʃi", "x"],  # IPA phrase with a prosodic bar
        ["a\\b", "c\\"],  # backslashes, one at the end of an item
        ["|", "\\|", "\\\\", ""],
        ["", "x"],
        ["only"],
    ],
)
def test_escaped_items_round_trip(tmp_path: Path, items: list[str]) -> None:
    cell = csvio.join_list(items)
    assert csvio.split_list(cell) == items
    path = tmp_path / "rt.csv"
    csvio.write_rows(path, [{"lemma": "x", "plurals": items}], FIELDS)
    (row,) = csvio.read_rows(path)
    assert csvio.split_list(row["plurals"]) == items


def test_plain_cells_read_and_write_as_before() -> None:
    # Items without "|" or a backslash are written exactly as before.
    assert csvio.join_list(["chirurghi", "chirurgi"]) == "chirurghi|chirurgi"
    assert csvio.join_list(("x", 1, 0.5, False)) == "x|1|0.5|False"
    # Cells without a backslash split exactly as before.
    for cell in ("a|b", "a||b", "|", "a"):
        assert csvio.split_list(cell) == cell.split("|")
    # A backslash before anything but "|" or a backslash is literal.
    assert csvio.split_list("C:\\x|y") == ["C:\\x", "y"]
    assert csvio.split_list("dangling\\") == ["dangling\\"]
    assert csvio.escape_list_item("a|b\\") == "a\\|b\\\\"


def test_curated_list_guards_use_the_escaped_form(tmp_path: Path) -> None:
    from lexpipe import curated

    path = tmp_path / "decisions.csv"
    csvio.write_rows(
        path,
        [{"lexeme_id": "x", "guard_ipa": ["a | b", "c"]}],
        ["lexeme_id", "guard_ipa"],
    )
    rows = curated.load_curated(path)
    assert curated.guard_mismatches(rows, {"x": {"ipa": ["a | b", "c"]}}) == []


# ---------------------------------------------------------------------------
# Guarded writes (ctx)
# ---------------------------------------------------------------------------


def _ctx(root: Path, code: str = "it") -> Any:
    from lexpipe import spec
    from lexpipe.context import RunContext

    return RunContext(spec.LANGUAGES[code], root=root)


def test_write_rows_with_a_context(tmp_path: Path) -> None:
    ctx = _ctx(tmp_path)
    target = ctx.data_dir / "nc" / "worklist.csv"
    assert csvio.write_rows(target, [{"lemma": "x"}], ["lemma"], ctx=ctx) == 1
    assert csvio.read_rows(target) == [{"lemma": "x"}]
    for relative in (
        "data/sources/kaikki/x.csv",
        "curated/italian/x.csv",
        "data/leipzig/freq/x_freq.csv",
        "data/romanian_lexicon_with_freq.csv",
        "analysis/romanian_tp_summary_full.csv",
    ):
        with pytest.raises(PermissionError):
            csvio.write_rows(
                tmp_path / relative, [{"lemma": "x"}], ["lemma"], ctx=ctx
            )
        assert not (tmp_path / relative).exists()
    # Without a context nothing is guarded (unchanged behaviour).
    free = tmp_path / "curated" / "italian" / "x.csv"
    assert csvio.write_rows(free, [{"lemma": "x"}], ["lemma"]) == 1


def test_write_rows_follows_the_context_into_symlinks(tmp_path: Path) -> None:
    ctx = _ctx(tmp_path, "ro")
    outside = tmp_path / "analysis" / "romanian_tp_summary_full.csv"
    outside.parent.mkdir(parents=True)
    outside.write_text("legacy\n", encoding="utf-8")
    link = ctx.comparison_dir / "innocent.csv"
    link.parent.mkdir(parents=True)
    link.symlink_to(outside)
    with pytest.raises(PermissionError):
        csvio.write_rows(link, [{"lemma": "x"}], ["lemma"], ctx=ctx)
    assert outside.read_text(encoding="utf-8") == "legacy\n"


def test_write_tracked_with_a_context(tmp_path: Path) -> None:
    ctx = _ctx(tmp_path)
    rows = [{"lemma": "x", "source": "wiktionary_en"}]
    fields = ["lemma", "source"]
    path = ctx.analysis_dir / "x.csv"
    # root comes from the context (the repository is not consulted).
    assert csvio.write_tracked(path, rows, fields, ctx=ctx) == 1
    assert csvio.write_tracked(path, rows, fields, ctx=ctx, root=tmp_path) == 1
    with pytest.raises(ValueError, match="differs"):
        csvio.write_tracked(
            path, rows, fields, ctx=ctx, root=tmp_path / "other"
        )
    with pytest.raises(PermissionError):
        csvio.write_tracked(
            tmp_path / "analysis" / "romanian_x.csv", rows, fields, ctx=ctx
        )
    from lexpipe import licensing

    with pytest.raises(licensing.LicenceError):
        csvio.write_tracked(
            path, [{"lemma": "x", "source": "morfflex"}], fields, ctx=ctx
        )
    with pytest.raises(licensing.LicenceError, match="tracked location"):
        csvio.write_tracked(ctx.nc_dir / "x.csv", rows, fields, ctx=ctx)


def test_write_tracked_for_a_comparison_language(tmp_path: Path) -> None:
    ctx = _ctx(tmp_path, "ro")
    rows = [{"lemma": "x", "source": "wiktionary_ro"}]
    fields = ["lemma", "source"]
    assert csvio.write_tracked(
        ctx.comparison_dir / "x.csv", rows, fields, ctx=ctx
    )
    with pytest.raises(PermissionError):  # confined (context guard)
        csvio.write_tracked(
            tmp_path / "analysis" / "italian" / "x.csv", rows, fields, ctx=ctx
        )
    from lexpipe import licensing

    with pytest.raises(licensing.LicenceError):  # writable, not tracked
        csvio.write_tracked(
            ctx.data_dir / "romanian_compare_lexicon_raw.csv",
            rows,
            fields,
            ctx=ctx,
        )
