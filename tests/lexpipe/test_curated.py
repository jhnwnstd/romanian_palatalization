"""Tests for lexpipe.curated: loading and guard checking."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lexpipe import csvio  # noqa: E402
from lexpipe.curated import (  # noqa: E402
    CuratedGuardError,
    check_guards,
    guard_mismatches,
    load_curated,
)

FIELDS = ("lexeme_id", "plural_adjudicated", "guard_gender", "guard_plurals")


def _write(path: Path, rows: list[dict[str, object]]) -> Path:
    csvio.write_rows(path, rows, FIELDS)
    return path


def test_missing_file_is_empty(tmp_path: Path) -> None:
    assert load_curated(tmp_path / "absent.csv") == {}


def test_load_splits_values_and_guards(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "adjudications.csv",
        [
            {
                "lexeme_id": "it:chirurgo:noun:0:abcd1234",
                "plural_adjudicated": "chirurghi",
                "guard_gender": "m",
                "guard_plurals": ["chirurghi", "chirurgi"],
            }
        ],
    )
    rows = load_curated(path)
    row = rows["it:chirurgo:noun:0:abcd1234"]
    assert dict(row.values) == {"plural_adjudicated": "chirurghi"}
    assert dict(row.guards) == {"gender": "m", "plurals": "chirurghi|chirurgi"}


def test_duplicate_id_rejected(tmp_path: Path) -> None:
    row = {"lexeme_id": "it:x:noun:0:00000000", "guard_gender": "m"}
    path = _write(tmp_path / "dup.csv", [row, row])
    with pytest.raises(ValueError, match="duplicate"):
        load_curated(path)


def test_guards_pass_and_fail(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "adjudications.csv",
        [
            {
                "lexeme_id": "it:mago:noun:0:11111111",
                "plural_adjudicated": "maghi",
                "guard_gender": "m",
                "guard_plurals": ["maghi", "magi"],
            }
        ],
    )
    curated = load_curated(path)
    same = {
        "it:mago:noun:0:11111111": {
            "gender": "m",
            "plurals": ["maghi", "magi"],
        }
    }
    assert guard_mismatches(curated, same) == []
    check_guards(curated, same)

    changed = {
        "it:mago:noun:0:11111111": {"gender": "m", "plurals": ["maghi"]}
    }
    with pytest.raises(CuratedGuardError, match="guard 'plurals'"):
        check_guards(curated, changed)

    with pytest.raises(CuratedGuardError, match="lexeme not found"):
        check_guards(curated, {})


def test_default_rows_carry_their_key(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "a.csv",
        [{"lexeme_id": "it:x:noun:0:00000000", "guard_gender": "m"}],
    )
    (row,) = load_curated(path).values()
    assert row.key == ("it:x:noun:0:00000000",)
    assert row.lexeme_id == "it:x:noun:0:00000000"
    assert load_curated(path, "lexeme_id") == load_curated(path)


# ---------------------------------------------------------------------------
# Composite keys (cell-level decisions)
# ---------------------------------------------------------------------------

CELL_FIELDS = ("lexeme_id", "cell", "confirmed", "guard_forms")
JAZYK = "cs:jazyk:noun:0:12345678"


def _cells(path: Path, rows: list[dict[str, object]]) -> Path:
    csvio.write_rows(path, rows, CELL_FIELDS)
    return path


def test_composite_keys(tmp_path: Path) -> None:
    path = _cells(
        tmp_path / "cell_adjudications.csv",
        [
            {
                "lexeme_id": JAZYK,
                "cell": "LOC.SG",
                "confirmed": ["jazyce", "jazyku"],
                "guard_forms": ["jazyce", "jazyku"],
            },
            {
                "lexeme_id": JAZYK,
                "cell": "NOM.PL",
                "confirmed": "jazyky",
                "guard_forms": "jazyky",
            },
        ],
    )
    rows = load_curated(path, key_fields=("lexeme_id", "cell"))
    assert list(rows) == [(JAZYK, "LOC.SG"), (JAZYK, "NOM.PL")]
    loc = rows[(JAZYK, "LOC.SG")]
    assert loc.key == (JAZYK, "LOC.SG") and loc.lexeme_id == JAZYK
    assert dict(loc.values) == {"confirmed": "jazyce|jazyku"}  # no key
    assert dict(loc.guards) == {"forms": "jazyce|jazyku"}
    # The same file read by lexeme_id alone holds two rows for one id.
    with pytest.raises(ValueError, match="duplicate lexeme_id"):
        load_curated(path)

    current = {
        (JAZYK, "LOC.SG"): {"forms": ["jazyce", "jazyku"]},
        (JAZYK, "NOM.PL"): {"forms": "jazyky"},
    }
    assert guard_mismatches(rows, current) == []
    check_guards(rows, current)
    changed = dict(current)
    changed[(JAZYK, "NOM.PL")] = {"forms": "jazyky|jazykové"}
    with pytest.raises(CuratedGuardError, match="NOM.PL: guard 'forms'"):
        check_guards(rows, changed)
    del changed[(JAZYK, "LOC.SG")]
    problems = guard_mismatches(rows, changed)
    assert problems[0] == f"{JAZYK}/LOC.SG: lexeme not found"


def test_composite_key_errors(tmp_path: Path) -> None:
    row = {"lexeme_id": JAZYK, "cell": "LOC.SG", "guard_forms": "x"}
    path = _cells(tmp_path / "dup.csv", [row, row])
    with pytest.raises(ValueError, match="duplicate key"):
        load_curated(path, key_fields=("lexeme_id", "cell"))
    with pytest.raises(ValueError, match="no 'number' column"):
        load_curated(path, key_fields=("lexeme_id", "number"))
    for bad in ((), ("cell", "cell"), ("lexeme_id", "guard_forms")):
        with pytest.raises(ValueError):
            load_curated(path, key_fields=bad)
    assert load_curated(tmp_path / "absent.csv", key_fields=("a", "b")) == {}
