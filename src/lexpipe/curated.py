"""Curated, hand-maintained inputs keyed by ``lexeme_id`` with guard fields.

Curated files hold decisions a human made about individual lexemes:
plural adjudications looked up by hand in Treccani or IJP, overrides
of a mis-parsed paradigm, person/thing judgements for -logo nouns.
They live under ``curated/<language>/`` (tracked) and are read-only
to the pipeline (design rule R6).

A decision is only valid for the lexeme the curator looked at. Every
row therefore carries ``guard_*`` columns recording what the curator
saw (for example ``guard_gender`` or ``guard_plurals``). When a new
kaikki dump changes that lexeme, the guard no longer matches the
current data and :func:`check_guards` fails loudly, instead of
silently re-attaching the decision to a different homograph (design
rules R4 and R6). This is the failure the Romanian pipeline's
hand-patched ``raw_dex.csv`` rows are exposed to.

File format: a CSV written with :mod:`lexpipe.csvio` conventions and a
``lexeme_id`` column. Columns starting with ``guard_`` are guards; the
remaining columns are the curated values. Guard values are compared as
strings after the same rendering :func:`lexpipe.csvio.format_value`
applies, so list-valued guards are ``|``-joined in source order.

Composite keys
--------------
Some decisions concern one cell of a lexeme (a Czech LOC.SG
adjudication), so a file may hold several rows per lexeme. Pass the
key columns as a tuple, ``load_curated(path, key_fields=("lexeme_id",
"cell"))``: rows are then keyed by the tuple of their key values
(``("cs:jazyk:noun:0:…", "LOC.SG")``), a repeated tuple is refused,
and the key columns are left out of ``values``. The guards of such a
file are checked against current data keyed the same way. The
default, ``key_fields="lexeme_id"``, keeps one row per lexeme keyed
by the plain id.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any, Final, TypeVar, overload

from lexpipe import csvio

ID_FIELD: Final[str] = "lexeme_id"
GUARD_PREFIX: Final[str] = "guard_"

# A curated file's key: the plain lexeme_id, or a tuple of key values.
CuratedKey = TypeVar("CuratedKey", str, tuple[str, ...])


class CuratedGuardError(ValueError):
    """Raised when a curated row no longer matches the lexeme it was made for."""


@dataclass(frozen=True, slots=True)
class CuratedRow:
    """One curated decision: its lexeme, its values and its guards.

    ``key`` holds the values of the file's key columns, in their order
    (``(lexeme_id,)`` for a file keyed by lexeme_id alone).
    """

    lexeme_id: str
    values: Mapping[str, str]
    guards: Mapping[str, str]
    key: tuple[str, ...] = ()


def _key_label(key: str | tuple[str, ...]) -> str:
    return key if isinstance(key, str) else "/".join(key)


@overload
def load_curated(path: Path) -> dict[str, CuratedRow]: ...


@overload
def load_curated(path: Path, key_fields: str) -> dict[str, CuratedRow]: ...


@overload
def load_curated(
    path: Path, key_fields: tuple[str, ...]
) -> dict[tuple[str, ...], CuratedRow]: ...


def load_curated(
    path: Path, key_fields: str | tuple[str, ...] = ID_FIELD
) -> dict[str, CuratedRow] | dict[tuple[str, ...], CuratedRow]:
    """Read a curated CSV into ``{key: CuratedRow}``.

    With a single key column (the default, ``lexeme_id``) the keys are
    its values; with a tuple of key columns they are tuples of values
    (module docstring). A missing file means "no decisions yet" and
    gives an empty dict. Raises :class:`ValueError` when a key column
    is missing or a key occurs twice: a second row for the same key
    would make the decision order-dependent.
    """
    single = isinstance(key_fields, str)
    fields: tuple[str, ...] = (
        (key_fields,) if isinstance(key_fields, str) else tuple(key_fields)
    )
    if not fields or len(set(fields)) != len(fields):
        raise ValueError(f"key fields must be distinct and given: {fields}")
    if any(field.startswith(GUARD_PREFIX) for field in fields):
        raise ValueError(f"a guard column cannot be a key: {fields}")
    if not path.exists():
        return {}
    rows = csvio.read_rows(path)
    by_tuple: dict[tuple[str, ...], CuratedRow] = {}
    for line, row in enumerate(rows, start=2):
        for name in fields:
            if name not in row:
                raise ValueError(f"{path}: no {name!r} column")
        key = tuple(row[name] for name in fields)
        if key in by_tuple:
            what = (
                f"{fields[0]} {key[0]!r}"
                if single
                else f"key {dict(zip(fields, key, strict=True))!r}"
            )
            raise ValueError(f"{path}:{line}: duplicate {what}")
        guards = {
            name[len(GUARD_PREFIX) :]: value
            for name, value in row.items()
            if name.startswith(GUARD_PREFIX)
        }
        values = {
            name: value
            for name, value in row.items()
            if name not in fields
            and name != ID_FIELD
            and not name.startswith(GUARD_PREFIX)
        }
        by_tuple[key] = CuratedRow(
            lexeme_id=row.get(ID_FIELD, ""),
            values=MappingProxyType(values),
            guards=MappingProxyType(guards),
            key=key,
        )
    if single:
        return {key[0]: row for key, row in by_tuple.items()}
    return by_tuple


def guard_mismatches(
    curated: Mapping[CuratedKey, CuratedRow],
    current: Mapping[CuratedKey, Mapping[str, Any]],
) -> list[str]:
    """Describe every curated row whose guards disagree with ``current``.

    ``current`` maps each key (``lexeme_id``, or the tuple of a
    composite key) to the present field values (rendered with
    :func:`lexpipe.csvio.format_value` before comparing). A curated
    key missing from ``current`` is reported too: the lexeme (or its
    cell) disappeared or its id changed.
    """
    problems: list[str] = []
    for key in sorted(curated):
        row = curated[key]
        label = _key_label(key)
        fields = current.get(key)
        if fields is None:
            problems.append(f"{label}: lexeme not found")
            continue
        for name in sorted(row.guards):
            expected = row.guards[name]
            actual = csvio.format_value(fields.get(name))
            if actual != expected:
                problems.append(
                    f"{label}: guard {name!r} expected {expected!r}, "
                    f"found {actual!r}"
                )
    return problems


def check_guards(
    curated: Mapping[CuratedKey, CuratedRow],
    current: Mapping[CuratedKey, Mapping[str, Any]],
) -> None:
    """Raise :class:`CuratedGuardError` listing every guard mismatch."""
    problems = guard_mismatches(curated, current)
    if problems:
        listing = "\n  ".join(problems)
        raise CuratedGuardError(
            f"{len(problems)} curated decision(s) no longer match the data "
            f"(re-check them by hand, then update the guard columns):\n  {listing}"
        )
