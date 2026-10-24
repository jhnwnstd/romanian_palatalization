"""Deterministic JSONL, gzip and SQLite files for the Czech derivatives.

Every derivative is byte-identical across runs, hash seeds and input
orders (design rule R5), so stage manifests can be compared and the
determinism tests can hash outputs:

  - :func:`canonical_line` serializes one record as canonical JSON
    (sorted keys, no whitespace, non-ASCII kept), so equal records give
    equal lines and the lines themselves are a total sort key;
  - :func:`write_jsonl_gz` writes gzip with a zero mtime and no file
    name in the header, the two fields that otherwise differ per run;
  - :func:`build_sqlite` creates a database from rows that are first
    collected in a staging file and then inserted in a fixed order
    into a fresh file, so the page layout does not depend on the order
    rows arrived in.

All writers go through :meth:`lexpipe.context.RunContext.open_output`
(atomic rename, protected locations refused).
"""

from __future__ import annotations

import contextlib
import csv
import gzip
import json
import os
import sqlite3
from collections.abc import Iterable, Iterator, Sequence
from pathlib import Path
from typing import Any, Final

from lexpipe import kaikki
from lexpipe.context import RunContext

GZIP_LEVEL: Final[int] = 6
_BATCH: Final[int] = 50_000


def canonical_line(record: Any) -> str:
    """Canonical JSON of one record (no trailing newline)."""
    return json.dumps(
        record,
        sort_keys=True,
        ensure_ascii=False,
        separators=(",", ":"),
        allow_nan=False,
    )


def write_jsonl_gz(ctx: RunContext, path: Path, lines: Iterable[str]) -> int:
    """Write *lines* (each without newline) as gzip; return the count."""
    count = 0
    with ctx.open_output(path, "wb") as raw:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            fileobj=raw,
            mtime=0,
            compresslevel=GZIP_LEVEL,
        ) as fh:
            for line in lines:
                fh.write(line.encode("utf-8"))
                fh.write(b"\n")
                count += 1
    return count


def iter_csv(path: Path) -> Iterator[dict[str, str]]:
    """Rows of a CSV written by :mod:`lexpipe.csvio`, streamed.

    Like :func:`lexpipe.csvio.read_rows` (every cell a string, missing
    cells empty) without holding the file in memory.
    """
    csv.field_size_limit(1 << 30)
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if None in row:
                raise ValueError(
                    f"{path}:{reader.line_num}: more cells than columns"
                )
            yield {key: value or "" for key, value in row.items()}


def read_jsonl(path: Path) -> Iterator[dict[str, Any]]:
    """Records of a (gzipped) JSONL file, in file order."""
    return kaikki.iter_jsonl(path)


def build_sqlite(
    ctx: RunContext,
    path: Path,
    *,
    table: str,
    columns: Sequence[tuple[str, str]],
    rows: Iterable[Sequence[Any]],
    order_by: Sequence[str],
    indexes: Sequence[tuple[str, Sequence[str]]] = (),
    meta: Sequence[tuple[str, str]] = (),
    distinct: bool = False,
) -> int:
    """Create ``path`` holding *table* filled from *rows*; return the
    number of rows stored.

    Rows are staged in a temporary database next to the target, then
    copied into a fresh database ``ORDER BY`` *order_by* (with
    ``DISTINCT`` if asked), indexed, and renamed over *path*. A
    ``meta(key, value)`` table records *meta* (source hashes, counts).
    """
    target = ctx.check_writable(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    staging = target.with_name(f".{target.name}.{os.getpid()}.staging")
    final = target.with_name(f".{target.name}.{os.getpid()}.tmp")
    for leftover in (staging, final):
        leftover.unlink(missing_ok=True)
    column_sql = ", ".join(f"{name} {kind}" for name, kind in columns)
    names = ", ".join(name for name, _ in columns)
    marks = ", ".join("?" for _ in columns)
    try:
        with contextlib.closing(sqlite3.connect(staging)) as con:
            con.execute(f"CREATE TABLE {table} ({column_sql})")
            batch: list[Sequence[Any]] = []
            for row in rows:
                batch.append(row)
                if len(batch) >= _BATCH:
                    con.executemany(
                        f"INSERT INTO {table} VALUES ({marks})", batch
                    )
                    batch.clear()
            if batch:
                con.executemany(f"INSERT INTO {table} VALUES ({marks})", batch)
            con.commit()
        with contextlib.closing(sqlite3.connect(final)) as con:
            con.execute("PRAGMA journal_mode = OFF")
            con.execute(f"CREATE TABLE {table} ({column_sql})")
            con.execute("CREATE TABLE meta (key TEXT, value TEXT)")
            con.execute("ATTACH DATABASE ? AS staging", (str(staging),))
            select = "SELECT DISTINCT" if distinct else "SELECT"
            con.execute(
                f"INSERT INTO {table} ({names}) {select} {names} "
                f"FROM staging.{table} ORDER BY {', '.join(order_by)}"
            )
            count = int(
                con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
            )
            con.executemany(
                "INSERT INTO meta VALUES (?, ?)",
                sorted([*meta, ("rows", str(count))]),
            )
            con.commit()
            con.execute("DETACH DATABASE staging")
            for index_name, index_columns in indexes:
                con.execute(
                    f"CREATE INDEX {index_name} ON {table} "
                    f"({', '.join(index_columns)})"
                )
            con.commit()
        os.replace(final, target)
    finally:
        staging.unlink(missing_ok=True)
        final.unlink(missing_ok=True)
    return count
