"""CSV conventions for the Italian and Czech outputs.

Every CSV the language plugins write goes through :func:`write_rows`
(or :func:`write_tracked` for tracked files), so all of them read the
same way in Python and R:

  - booleans are ``True`` / ``False``, as in the Romanian outputs
    (design rule R10; ``docs/data_dict.md``'s "lowercase" claim is
    wrong), and :func:`parse_bool` reads any capitalisation;
  - None, NaN and the empty string are written as an empty cell;
  - lists and tuples are joined with ``|`` in source order, each item
    rendered by the same rules and then escaped: a backslash in an
    item is written ``\\\\`` and a ``|`` is written ``\\|``
    (:func:`escape_list_item`), so IPA with prosodic bars or glosses
    with a ``|`` survive a round trip;
  - sets are rejected: they have no source order, and the output
    must be byte-identical across runs;
  - floats use Python's shortest round-trip repr; everything else
    is ``str(value)`` (so a StrEnum writes its value; a plain string
    is never escaped);
  - UTF-8, ``\\n`` line endings, minimal quoting, header first.

Files are written to a temporary name and renamed, so a reader never
sees half a table. :func:`read_rows` gives every cell back as a
string; :func:`parse_bool` and :func:`split_list` undo the encoding.
List items without a backslash or ``|`` are written exactly as before
escaping existed, and :func:`split_list` reads a cell without a
backslash exactly as before; a backslash that precedes neither
``\\`` nor ``|`` is kept as it is.

Guarded writes
--------------
:func:`write_rows` and :func:`write_tracked` take an optional
``ctx`` (a :class:`lexpipe.context.RunContext`): the path then goes
through ``ctx.check_writable`` first (the protected locations of R1,
R2 and R6, and the confinement of comparison-only languages), so a
stage need not call it before every write, and the tracked-location
check of :func:`write_tracked` uses ``ctx.root``.

Tracked outputs
---------------
:func:`write_tracked` is the licence gate of R9: it refuses a path
outside the tracked locations, a header without the source column,
and any row whose source cell (``|``-joined allowed) names a source
whose data may not be redistributed (or an unknown one), before the
target file is touched; a ``licence_class`` column is checked the
same way. See :mod:`lexpipe.licensing`. NC-bearing outputs are
written with :func:`write_rows` under ``data/<language>/nc/``.

The Romanian pipeline predates these rules (CRLF line endings, its
own writers) and is deliberately left as it is.
"""

from __future__ import annotations

import csv
import math
import os
from collections.abc import Iterable, Iterator, Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any, Final, Literal

from lexpipe import PROJECT_ROOT, licensing
from lexpipe.spec import LANGUAGES, LanguageInfo

if TYPE_CHECKING:
    from lexpipe.context import RunContext

LIST_SEPARATOR: Final[str] = "|"
LIST_ESCAPE: Final[str] = "\\"
TRUE: Final[str] = "True"
FALSE: Final[str] = "False"
LINE_TERMINATOR: Final[str] = "\n"
LICENCE_CLASS_FIELD: Final[str] = "licence_class"


def escape_list_item(text: str) -> str:
    """Escape one rendered list item: ``\\`` → ``\\\\``, ``|`` → ``\\|``."""
    return text.replace(LIST_ESCAPE, LIST_ESCAPE * 2).replace(
        LIST_SEPARATOR, LIST_ESCAPE + LIST_SEPARATOR
    )


def format_value(value: Any) -> str:
    """Render one cell according to the module conventions."""
    if value is None:
        return ""
    if isinstance(value, bool):
        return TRUE if value else FALSE
    if isinstance(value, float):
        return "" if math.isnan(value) else repr(value)
    if isinstance(value, (list, tuple)):
        return LIST_SEPARATOR.join(
            escape_list_item(format_value(item)) for item in value
        )
    if isinstance(value, (set, frozenset)):
        raise TypeError("sets have no source order; pass a list")
    return str(value)


def join_list(items: Iterable[Any]) -> str:
    """One list cell from *items* (in their order), as :func:`format_value`
    writes a list."""
    return format_value(list(items))


def _guarded(path: Path, ctx: RunContext | None) -> Path:
    return path if ctx is None else ctx.check_writable(path)


def write_rows(
    path: Path,
    rows: Iterable[Mapping[str, Any]],
    fields: Sequence[str],
    *,
    extra: Literal["raise", "ignore"] = "raise",
    ctx: RunContext | None = None,
) -> int:
    """Write *rows* under the header *fields*; return the row count.

    Keys missing from a row become empty cells. Keys not in *fields*
    raise :class:`ValueError` unless ``extra="ignore"``. If iterating
    *rows* raises, the previous file (if any) is left untouched. With
    *ctx*, *path* must pass ``ctx.check_writable`` (else
    :class:`PermissionError`, before anything is written) and the
    resolved path it returns is written.
    """
    path = _guarded(path, ctx)
    columns = list(fields)
    if len(set(columns)) != len(columns):
        raise ValueError(f"duplicate column names: {columns}")
    known = set(columns)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    count = 0
    try:
        with tmp.open("w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh, lineterminator=LINE_TERMINATOR)
            writer.writerow(columns)
            for row in rows:
                unknown = [key for key in row if key not in known]
                if unknown and extra == "raise":
                    raise ValueError(
                        f"row {count + 1} has columns not in fields: "
                        f"{unknown}"
                    )
                writer.writerow(
                    [format_value(row.get(column)) for column in columns]
                )
                count += 1
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    os.replace(tmp, path)
    return count


def _checked(
    rows: Iterable[Mapping[str, Any]],
    source_field: str,
    check_class: bool,
) -> Iterator[Mapping[str, Any]]:
    """Yield *rows*, raising :class:`licensing.LicenceError` on the
    first row that may not be redistributed."""
    for number, row in enumerate(rows, 1):
        value = row.get(source_field)
        if value is None:
            value = ""
        try:
            licensing.check_redistributable(
                licensing.split_sources(
                    value if isinstance(value, (list, tuple)) else str(value)
                )
            )
            if check_class:
                cell = row.get(LICENCE_CLASS_FIELD)
                if cell not in (None, ""):
                    licensing.check_licence_class(format_value(cell))
        except licensing.LicenceError as exc:
            raise licensing.LicenceError(f"row {number}: {exc}") from None
        yield row


def write_tracked(
    path: Path,
    rows: Iterable[Mapping[str, Any]],
    fields: Sequence[str],
    *,
    source_field: str = "source",
    root: Path | None = None,
    extra: Literal["raise", "ignore"] = "raise",
    ctx: RunContext | None = None,
) -> int:
    """:func:`write_rows` for a tracked file, behind the licence gate.

    Raises :class:`lexpipe.licensing.LicenceError` if *path* is not a
    tracked location under *root*, if *fields* lacks *source_field*,
    or if a row's sources (or ``licence_class``, when that column
    exists) are not all redistributable. Nothing is written then: an
    existing file keeps its previous content.

    *root* defaults to ``ctx.root`` with a *ctx* and to the repository
    otherwise; a *root* that differs from ``ctx.root`` raises
    :class:`ValueError`. With *ctx*, *path* must also pass
    ``ctx.check_writable`` (:class:`PermissionError`), and the tracked
    locations are those of the registry plus ``ctx.language`` when it
    is not registered (:func:`lexpipe.licensing.with_language`).
    """
    if ctx is None:
        base = PROJECT_ROOT if root is None else root
        languages: Mapping[str, LanguageInfo] = LANGUAGES
    else:
        if root is not None and Path(root).resolve() != ctx.root:
            raise ValueError(
                f"root {root} differs from the run context's {ctx.root}"
            )
        path = ctx.check_writable(path)
        base = ctx.root
        languages = licensing.with_language(ctx.language)
    licensing.check_tracked_location(path, base, languages)
    if source_field not in fields:
        raise licensing.LicenceError(
            f"tracked tables need a {source_field!r} column (R9)"
        )
    checked = _checked(
        rows, source_field, check_class=LICENCE_CLASS_FIELD in fields
    )
    return write_rows(path, checked, fields, extra=extra)


def read_rows(path: Path) -> list[dict[str, str]]:
    """Read a CSV written by :func:`write_rows` (all cells as strings).

    Short rows are padded with empty cells; a row with more cells than
    the header raises :class:`ValueError`.
    """
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if None in row:
                raise ValueError(
                    f"{path}:{reader.line_num}: more cells than columns"
                )
            rows.append({key: value or "" for key, value in row.items()})
    return rows


def parse_bool(text: str) -> bool | None:
    """``True``/``true``/``TRUE`` → True, likewise False; ``""`` → None.

    Any capitalisation is accepted (R10), so the Romanian ``True`` and
    readr's ``TRUE`` read the same.
    """
    lowered = text.strip().lower()
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    if lowered == "":
        return None
    raise ValueError(f"not a boolean cell: {text!r}")


def split_list(text: str) -> list[str]:
    """Undo the ``|`` join and the item escapes; an empty cell is an
    empty list.

    ``\\|`` is a literal ``|`` and ``\\\\`` a literal backslash; any
    other backslash is kept as it is, so a cell without escapes reads
    exactly as it did before escaping existed.
    """
    if not text:
        return []
    if LIST_ESCAPE not in text:
        return text.split(LIST_SEPARATOR)
    items: list[str] = []
    current: list[str] = []
    index, end = 0, len(text)
    while index < end:
        char = text[index]
        if (
            char == LIST_ESCAPE
            and index + 1 < end
            and text[index + 1] in (LIST_ESCAPE, LIST_SEPARATOR)
        ):
            current.append(text[index + 1])
            index += 2
            continue
        if char == LIST_SEPARATOR:
            items.append("".join(current))
            current = []
        else:
            current.append(char)
        index += 1
    items.append("".join(current))
    return items
