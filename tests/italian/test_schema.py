"""The Italian schema (``docs/italian/schema.toml``) against the code.

Every CSV writer's header constant must equal the schema's column list,
in order; the data dictionary renders every file; controlled columns of
the files on disk stay within their allowed values (needs_data).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lexpipe.langs.italian import schema
from tests.italian import helpers

pytestmark = pytest.mark.it

SCHEMA = helpers.SCHEMA


@pytest.fixture(scope="module")
def loaded() -> dict[str, schema.FileSchema]:
    return schema.load_schema(SCHEMA)


def test_every_writer_header_equals_the_schema(
    loaded: dict[str, schema.FileSchema],
) -> None:
    assert schema.header_problems(loaded) == []
    for item in loaded.values():
        assert schema.code_fields(item.fields) == item.names, item.key


def test_every_italian_csv_is_described(
    loaded: dict[str, schema.FileSchema],
) -> None:
    paths = {item.path for item in loaded.values()}
    expected = {
        "data/italian/italian_lexicon_raw.csv",
        "data/italian/italian_lexicon_full.csv",
        "data/italian/italian_plural_forms.csv",
        "data/italian/italian_plural_drift.csv",
        "data/italian/italian_qc.csv",
        "analysis/italian/italian_manual_worklist.csv",
        "data/italian/italian_attestation.csv",
        "data/italian/italian_cells.csv",
        "data/italian/italian_lexicon_complete.csv",
        "data/italian/italian_derivations.csv",
        "data/italian/italian_verbs.csv",
        "analysis/italian/italian_tp_long.csv",
        "analysis/italian/italian_learners.csv",
        "analysis/italian/italian_logtable.csv",
        "data/italian/italian_logtable_attestation.csv",
    }
    assert paths == expected
    for item in loaded.values():
        tracked = item.path.startswith("analysis/") or item.path.endswith(
            "lexicon_raw.csv"
        )
        assert item.tracked == tracked, item.key
        for column in item.columns:
            assert column.doc, (item.key, column.name)
            if item.tracked:  # Leipzig counts never in tracked files (R9)
                assert column.licence != "unverified", (item.key, column.name)


def test_data_dict_renders_every_column(
    loaded: dict[str, schema.FileSchema],
) -> None:
    text = schema.render_data_dict(loaded)
    for item in loaded.values():
        assert f"## `{item.path}`" in text
        for column in item.columns:
            assert f"| `{column.name}` |" in text


def test_malformed_schema_is_refused(tmp_path: Path) -> None:
    bad = tmp_path / "schema.toml"
    bad.write_text(
        '[files.x]\npath = "a.csv"\nstage = "s"\ntracked = false\n'
        'fields = "m:F"\ncolumns = [{ name = "a", type = "text", '
        'provenance = "derived", licence = "by-sa", doc = "d" }]\n',
        encoding="utf-8",
    )
    with pytest.raises(schema.SchemaError):
        schema.load_schema(bad)


@pytest.mark.needs_data
def test_files_on_disk_match_the_schema(
    loaded: dict[str, schema.FileSchema],
) -> None:
    root = helpers.PROJECT_ROOT
    present = [i for i in loaded.values() if (root / i.path).is_file()]
    if not present:
        pytest.skip("no Italian outputs on disk")
    assert schema.header_problems(loaded, root) == []
    for item in present:
        problems = schema.value_problems(item, root / item.path)
        assert problems == [], problems[:5]
