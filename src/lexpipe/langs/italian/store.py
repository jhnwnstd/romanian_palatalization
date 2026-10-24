"""Compact, source-stamped derivatives of the Italian pipeline (R5, R7).

The big sources are read once by the ``extract`` stage and projected
into small files under ``data/italian/extract/``; later stages read only
those. Two formats, both gzip-compressed:

  - **JSONL** (projections of dictionary entries): the first line is
    ``{"_meta": {...}}``, every other line one record written with
    :func:`dumps` (sorted keys, no spaces, non-ASCII kept);
  - **TSV** (the form index): the first line is ``#`` followed by the
    meta JSON, the second the column names, then one row per line.

The meta object stamps a derivative with the format name and the
SHA-256 of every source file it was built from, so a reader can refuse
a projection made from another dump (:func:`check_stamp`).

Determinism: records are written in the order given (callers sort
them), and the gzip header carries neither a file name nor a time
stamp (``mtime=0``), so identical content gives identical bytes.
Files are written through :meth:`lexpipe.context.RunContext.open_output`
(atomic, guarded).
"""

from __future__ import annotations

import gzip
import io
import json
from collections.abc import Iterable, Iterator, Mapping, Sequence
from pathlib import Path
from typing import IO, Any, Final

from lexpipe.context import RunContext

META_KEY: Final[str] = "_meta"
TSV_META_PREFIX: Final[str] = "#"
COMPRESSLEVEL: Final[int] = 6


class StaleDerivativeError(ValueError):
    """A derivative was built from other source files than expected."""


def dumps(obj: Any) -> str:
    """Canonical one-line JSON (sorted keys, compact, UTF-8 text)."""
    return json.dumps(
        obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")
    )


def _gzip_writer(raw: IO[Any]) -> io.TextIOWrapper:
    gz = gzip.GzipFile(
        filename="",
        mode="wb",
        fileobj=raw,
        mtime=0,
        compresslevel=COMPRESSLEVEL,
    )
    return io.TextIOWrapper(gz, encoding="utf-8", newline="\n")


def write_jsonl(
    ctx: RunContext,
    path: Path,
    meta: Mapping[str, Any],
    lines: Iterable[str],
) -> int:
    """Write a stamped JSONL.gz; *lines* are pre-serialized records.

    Returns the number of records written.
    """
    count = 0
    with ctx.open_output(path, "wb") as raw:
        with _gzip_writer(raw) as fh:
            fh.write(dumps({META_KEY: dict(meta)}) + "\n")
            for line in lines:
                fh.write(line + "\n")
                count += 1
    return count


def read_jsonl(path: Path) -> tuple[dict[str, Any], Iterator[dict[str, Any]]]:
    """``(meta, records)`` of a stamped JSONL.gz (records are lazy)."""
    fh = gzip.open(path, "rt", encoding="utf-8")
    first = fh.readline()
    try:
        head = json.loads(first) if first.strip() else {}
    except json.JSONDecodeError as exc:
        fh.close()
        raise ValueError(f"{path}: first line is not JSON meta") from exc
    if not isinstance(head, dict) or META_KEY not in head:
        fh.close()
        raise ValueError(f"{path}: no {META_KEY!r} line")
    meta = head[META_KEY]

    def records() -> Iterator[dict[str, Any]]:
        with fh:
            for lineno, line in enumerate(fh, 2):
                if not line.strip():
                    continue
                obj = json.loads(line)
                if not isinstance(obj, dict):
                    raise ValueError(f"{path}:{lineno}: not a JSON object")
                yield obj

    return dict(meta), records()


def write_tsv(
    ctx: RunContext,
    path: Path,
    meta: Mapping[str, Any],
    fields: Sequence[str],
    rows: Iterable[str],
) -> int:
    """Write a stamped TSV.gz; *rows* are tab-joined lines."""
    count = 0
    with ctx.open_output(path, "wb") as raw:
        with _gzip_writer(raw) as fh:
            fh.write(TSV_META_PREFIX + dumps(dict(meta)) + "\n")
            fh.write("\t".join(fields) + "\n")
            for row in rows:
                fh.write(row + "\n")
                count += 1
    return count


def read_tsv(
    path: Path,
) -> tuple[dict[str, Any], list[str], Iterator[list[str]]]:
    """``(meta, fields, rows)`` of a stamped TSV.gz (rows are lazy)."""
    fh = gzip.open(path, "rt", encoding="utf-8", newline="\n")
    first = fh.readline()
    if not first.startswith(TSV_META_PREFIX):
        fh.close()
        raise ValueError(f"{path}: no meta line")
    meta = json.loads(first[len(TSV_META_PREFIX) :])
    fields = fh.readline().rstrip("\n").split("\t")

    def rows() -> Iterator[list[str]]:
        with fh:
            for line in fh:
                text = line.rstrip("\n")
                if text:
                    yield text.split("\t")

    return dict(meta), fields, rows()


def check_stamp(
    meta: Mapping[str, Any], fmt: str, sources: Mapping[str, str]
) -> None:
    """Raise :class:`StaleDerivativeError` unless *meta* has format
    *fmt* and exactly the source SHA-256s in *sources*."""
    if meta.get("format") != fmt:
        raise StaleDerivativeError(
            f"expected format {fmt!r}, found {meta.get('format')!r}"
        )
    stamped = meta.get("sources")
    if stamped != dict(sources):
        raise StaleDerivativeError(
            f"derivative built from {stamped!r}, expected {dict(sources)!r}"
        )


def tsv_cell(value: str) -> str:
    """A value safe for a TSV cell (tabs and newlines become spaces)."""
    return value.replace("\t", " ").replace("\n", " ").replace("\r", " ")
