"""Schema of the Czech CSVs and the generated data dictionary (design §8)."""

from __future__ import annotations

import dataclasses
from pathlib import Path

import pytest

from lexpipe import licensing, spec
from lexpipe.context import RunContext
from lexpipe.langs.czech import analysis_layout, layout, schema
from tests.czech import treebuilder

pytestmark = pytest.mark.cs

ROOT = treebuilder.PROJECT_ROOT


def test_every_column_is_documented() -> None:
    for table in schema.tables():
        for column in table.columns():
            assert column.description, (table.name, column.name)
            assert column.kind in (
                "str",
                "int",
                "float",
                "bool",
                "list",
                "enum",
            )


def test_undocumented_column_raises() -> None:
    with pytest.raises(KeyError):
        schema.describe("no_such_column")


def test_patterned_columns() -> None:
    assert schema.describe("cell_LOC.PL").kind == "list"
    assert "ge2pm" in schema.describe("attested_ge2pm").description
    strict = schema.describe("summary_E0_ge1_strict").description
    assert "strict" in strict


def test_table_paths_match_the_layouts() -> None:
    ctx = RunContext(spec.LANGUAGES["cs"], root=ROOT)
    known = {
        path.relative_to(ROOT).as_posix()
        for paths in (
            layout.paths_for(ctx),
            analysis_layout.analysis_paths(ctx),
        )
        for path in (getattr(paths, f.name) for f in dataclasses.fields(paths))
        if isinstance(path, Path)
    }
    for table in schema.tables():
        assert table.path in known, table.path


def test_tracked_tables_are_tracked_locations() -> None:
    for table in schema.tables():
        assert licensing.is_tracked_location(ROOT / table.path, ROOT) == (
            table.tracked
        ), table.name


def test_committed_data_dictionary_is_up_to_date() -> None:
    path = ROOT / "docs" / "czech" / "data_dict.md"
    assert path.read_text(encoding="utf-8") == schema.render()


@pytest.mark.needs_data
def test_real_headers_match_the_schema() -> None:
    found = 0
    for table in schema.tables():
        path = ROOT / table.path
        if not path.is_file():
            continue
        with path.open(encoding="utf-8") as fh:
            header = fh.readline().rstrip("\n").split(",")
        assert header == list(table.fields), table.name
        found += 1
    if not found:
        pytest.skip("no Czech outputs on disk")


def test_value_validation_flags_values_outside_the_sets() -> None:
    table = next(t for t in schema.tables() if t.name == "tp_long")
    good = {
        "rule": "palatalize",
        "target_class": "velar:k",
        "tau": "ge1_strict",
        "evidence_set": "E0",
        "coding": "dominant",
        "freq_basis": "paradigm_unambiguous",
        "tolerated_trigger": "True",
        "licence_class": "by-sa",
    }
    assert schema.validate_rows(table, [good]) == []
    bad = dict(good, tau="ge3", tolerated_trigger="yes")
    problems = schema.validate_rows(table, [bad])
    assert any("tp_long.tau" in p and "'ge3'" in p for p in problems)
    assert any("tolerated_trigger" in p for p in problems)


def test_list_columns_check_every_item() -> None:
    table = next(t for t in schema.tables() if t.name == "spellout_roots")
    assert schema.validate_rows(table, [{"basis": "linked|wikt"}]) == []
    assert schema.validate_rows(table, [{"basis": "linked|nope"}])


def test_family_enums_have_value_sets() -> None:
    assert "both" in schema.describe("choice_ge1").values
    assert "mixed" in schema.describe("summary_E0_ge1_strict").values


@pytest.mark.needs_data
def test_real_values_are_in_the_documented_sets() -> None:
    found = 0
    for table in schema.tables():
        if not (ROOT / table.path).is_file():
            continue
        assert schema.validate_file(table, ROOT) == [], table.name
        found += 1
    if not found:
        pytest.skip("no Czech outputs on disk")
