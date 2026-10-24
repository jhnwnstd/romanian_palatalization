"""The schema of every Italian CSV and its data dictionary (design §8).

``docs/italian/schema.toml`` is the source of truth: one table per file
(path, stage, tracked, the code constant holding its header, a
description) with one entry per column (name, type, allowed values,
provenance class, licence class, doc; the file's header comment defines
the vocabularies). This module reads it (:func:`load_schema`), checks it
against the code and the files (:func:`header_problems`,
:func:`value_problems`) and renders ``docs/italian/data_dict.md``
(:func:`render_data_dict`). ``tests/italian/test_schema.py`` asserts that
every writer's header equals the schema; the ``docs`` stage
(:mod:`.docstage`) refuses to write the dictionary while one does not.
"""

from __future__ import annotations

import csv
import importlib
import tomllib
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

TYPES: Final[frozenset[str]] = frozenset(
    {"str", "int", "float", "bool", "list", "enum"}
)
PROVENANCE: Final[frozenset[str]] = frozenset(
    {
        "wiktionary_en",
        "wiktionary_it",
        "morphit",
        "leipzig",
        "module",
        "curated",
        "derived",
    }
)
LICENCES: Final[frozenset[str]] = frozenset(
    {"by-sa", "unverified", "citation"}
)


@dataclass(frozen=True, slots=True)
class Column:
    """One column of a file."""

    name: str
    type: str
    allowed: tuple[str, ...]
    provenance: str
    licence: str
    doc: str


@dataclass(frozen=True, slots=True)
class FileSchema:
    """One CSV of the Italian pipeline."""

    key: str
    path: str
    stage: str
    tracked: bool
    fields: str
    description: str
    columns: tuple[Column, ...]

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(column.name for column in self.columns)


class SchemaError(ValueError):
    """The schema file is malformed or disagrees with the code."""


def _column(raw: Mapping[str, Any], where: str) -> Column:
    name = str(raw.get("name", ""))
    kind = str(raw.get("type", ""))
    provenance = str(raw.get("provenance", ""))
    licence = str(raw.get("licence", ""))
    if not name or kind not in TYPES:
        raise SchemaError(f"{where}: bad column {raw!r}")
    if provenance not in PROVENANCE or licence not in LICENCES:
        raise SchemaError(f"{where}.{name}: bad provenance/licence")
    allowed = raw.get("allowed", [])
    if not isinstance(allowed, list):
        raise SchemaError(f"{where}.{name}: allowed must be a list")
    return Column(
        name=name,
        type=kind,
        allowed=tuple(str(a) for a in allowed),
        provenance=provenance,
        licence=licence,
        doc=str(raw.get("doc", "")),
    )


def load_schema(path: Path) -> dict[str, FileSchema]:
    """Every file of the schema, in file order."""
    with path.open("rb") as fh:
        data = tomllib.load(fh)
    files = data.get("files")
    if not isinstance(files, dict):
        raise SchemaError(f"{path}: no [files] tables")
    out: dict[str, FileSchema] = {}
    for key, raw in files.items():
        columns = tuple(
            _column(item, f"{path}:{key}") for item in raw.get("columns", [])
        )
        names = [c.name for c in columns]
        if len(set(names)) != len(names):
            raise SchemaError(f"{path}:{key}: duplicate column names")
        out[key] = FileSchema(
            key=key,
            path=str(raw["path"]),
            stage=str(raw["stage"]),
            tracked=bool(raw["tracked"]),
            fields=str(raw["fields"]),
            description=str(raw.get("description", "")),
            columns=columns,
        )
    return out


def code_fields(spec: str) -> tuple[str, ...]:
    """The header constant named ``"module:CONSTANT"``."""
    module_name, _, attribute = spec.partition(":")
    value = getattr(importlib.import_module(module_name), attribute)
    return tuple(str(item) for item in value)


def file_header(path: Path) -> tuple[str, ...]:
    """The first line of a CSV, parsed."""
    with path.open("r", encoding="utf-8", newline="") as fh:
        return tuple(next(csv.reader(fh), []))


def header_problems(
    schema: Mapping[str, FileSchema], root: Path | None = None
) -> list[str]:
    """Differences between the schema and the code constants (and the
    files under *root*, when given and present)."""
    problems: list[str] = []
    for item in schema.values():
        expected = item.names
        found = code_fields(item.fields)
        if found != expected:
            problems.append(
                f"{item.key}: code header {item.fields} differs from the "
                f"schema: missing {sorted(set(found) - set(expected))}, "
                f"extra {sorted(set(expected) - set(found))}, or order"
            )
        if root is not None and (root / item.path).is_file():
            on_disk = file_header(root / item.path)
            if on_disk != expected:
                problems.append(f"{item.path}: header differs from the schema")
    return problems


def _allowed_values(column: Column, cell: str) -> Iterable[str]:
    if column.type == "list":
        return cell.split("|") if cell else []
    return [cell]


def value_problems(
    item: FileSchema, path: Path, max_rows: int | None = None
) -> list[str]:
    """Cells of controlled columns whose value is not allowed."""
    constrained = {c.name: c for c in item.columns if c.allowed}
    if not constrained:
        return []
    problems: list[str] = []
    with path.open("r", encoding="utf-8", newline="") as fh:
        for number, row in enumerate(csv.DictReader(fh), 2):
            if max_rows is not None and number > max_rows + 1:
                break
            for name, column in constrained.items():
                for value in _allowed_values(column, row.get(name) or ""):
                    if value not in column.allowed and not (
                        value == "" and column.type == "bool"
                    ):
                        problems.append(
                            f"{path.name}:{number}: {name} = {value!r}"
                        )
    return problems


def _escape(text: str) -> str:
    return text.replace("|", "\\|").replace("\n", " ")


def render_data_dict(schema: Mapping[str, FileSchema]) -> str:
    """``docs/italian/data_dict.md`` (generated; do not edit)."""
    lines = [
        "# Italian data dictionary",
        "",
        "Generated by the `docs` stage from `docs/italian/schema.toml` "
        "(`lexpipe.langs.italian.schema`); edit the schema, not this file.",
        "",
        "Types: `bool` cells are `True`/`False` (R10), `list` cells are "
        "`|`-joined, empty cells are missing values. Provenance: "
        "`wiktionary_en`, `wiktionary_it`, `morphit`, `leipzig` (corpus "
        "counts), `module` (Module:it-common port), `curated`, `derived`. "
        "Licence: `by-sa` (Wiktionary, Morph-it!), `unverified` (Leipzig "
        "counts: never in tracked files, R9), `citation` (hand "
        "adjudication).",
        "",
        "## Files",
        "",
        "| file | stage | tracked | columns |",
        "|---|---|---|---|",
    ]
    for item in schema.values():
        lines.append(
            f"| `{item.path}` | {item.stage} | "
            f"{'yes' if item.tracked else 'no'} | {len(item.columns)} |"
        )
    for item in schema.values():
        lines += [
            "",
            f"## `{item.path}`",
            "",
            _escape(item.description),
            "",
            f"Stage `{item.stage}`; header `{item.fields}`; "
            f"{'tracked' if item.tracked else 'gitignored'}.",
            "",
            "| column | type | allowed | provenance | licence | description |",
            "|---|---|---|---|---|---|",
        ]
        for column in item.columns:
            allowed = ", ".join(f"`{a}`" for a in column.allowed if a)
            lines.append(
                f"| `{column.name}` | {column.type} | {_escape(allowed)} | "
                f"{column.provenance} | {column.licence} | "
                f"{_escape(column.doc)} |"
            )
    return "\n".join(lines) + "\n"


def tracked_files(schema: Mapping[str, FileSchema]) -> Sequence[FileSchema]:
    """Files the schema marks as tracked."""
    return [item for item in schema.values() if item.tracked]
