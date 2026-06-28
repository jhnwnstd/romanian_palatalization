#!/usr/bin/env python3
"""One-shot migration: ``dex_cache.json`` (6.9 GB JSON) → ``dex_cache.db``.

The JSON cache is a flat ``{url: html}`` mapping accumulated by the DEX
QC and supplement runs. Loading or streaming it once per analysis is
~30 s of overhead even for a few hundred lookups. Moving to a SQLite
k/v store makes random lookups O(log n) and avoids the multi-GB load
entirely.

Memory and disk safety
----------------------
- The JSON is *streamed* via :mod:`ijson` (``kvitems``). At no point is
  the full file loaded into memory; we hold at most one
  ``BATCH_SIZE``-sized list of ``(url, html)`` pairs.
- The DB uses the default DELETE rollback journal — NOT WAL — so the
  per-transaction journal is small and is unlinked after each commit.
  A WAL approach can grow the ``.wal`` file to the size of all pending
  writes; killing the process leaves the data stranded in the WAL.
- Each batch is its own transaction. If the process is killed
  mid-migration, the already-committed rows persist and the next run
  resumes via ``INSERT OR REPLACE`` (idempotent).

Schema
------
::

    CREATE TABLE cache (
        url   TEXT PRIMARY KEY,
        html  TEXT NOT NULL
    );

Usage
-----
::

    python scripts/migrate_dex_cache_to_sqlite.py [--json PATH] [--db PATH]

After successful migration, the JSON file can be deleted (it is
already gitignored).
"""

from __future__ import annotations

import argparse
import sqlite3
import sys
import time
from pathlib import Path
from typing import Final

import ijson

PROJECT_ROOT: Final[Path] = Path(__file__).resolve().parent.parent
DEFAULT_JSON: Final[Path] = PROJECT_ROOT / "dex_cache.json"
DEFAULT_DB: Final[Path] = PROJECT_ROOT / "dex_cache.db"

# Rows per transaction. Large enough for throughput, small enough that
# the rollback journal stays a few MB instead of multi-GB.
BATCH_SIZE: Final[int] = 1_000

# Progress print cadence (in committed rows).
LOG_EVERY: Final[int] = 10_000


def _open_db(db_path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    # DELETE journal mode (default) — the journal lives only for the
    # duration of a transaction and is unlinked on commit. No risk of
    # a multi-GB orphan .wal file if the process dies.
    conn.execute("PRAGMA journal_mode = DELETE")
    # Default synchronous=FULL is fine here: this is a one-shot ETL,
    # not a hot write path. We'd rather pay a small per-commit fsync
    # cost than risk corruption if the box loses power mid-migration.
    conn.execute(
        "CREATE TABLE IF NOT EXISTS cache ("
        "url TEXT PRIMARY KEY, html TEXT NOT NULL)"
    )
    conn.commit()
    return conn


def migrate(json_path: Path, db_path: Path) -> int:
    """Stream the JSON cache into SQLite. Returns the row count written.

    Returns -1 if the source JSON is missing.
    """
    if not json_path.exists():
        print(f"ERROR: source JSON not found: {json_path}", file=sys.stderr)
        return -1
    print(f"Source: {json_path} ({json_path.stat().st_size / 1e9:.2f} GB)")
    print(f"Target: {db_path}")

    conn = _open_db(db_path)
    try:
        t0 = time.perf_counter()
        n_written = 0
        batch: list[tuple[str, str]] = []
        with json_path.open("rb") as f:
            for key, value in ijson.kvitems(f, ""):
                batch.append((key, value))
                if len(batch) >= BATCH_SIZE:
                    conn.executemany(
                        "INSERT OR REPLACE INTO cache "
                        "(url, html) VALUES (?, ?)",
                        batch,
                    )
                    conn.commit()
                    n_written += len(batch)
                    batch.clear()
                    if n_written % LOG_EVERY == 0:
                        elapsed = time.perf_counter() - t0
                        rate = n_written / elapsed
                        print(
                            f"  wrote {n_written:>8,d} rows  "
                            f"({rate:.0f}/s, {elapsed:.0f}s elapsed)",
                            flush=True,
                        )
        if batch:
            conn.executemany(
                "INSERT OR REPLACE INTO cache (url, html) VALUES (?, ?)",
                batch,
            )
            conn.commit()
            n_written += len(batch)
        elapsed = time.perf_counter() - t0
        print(
            f"\nDone. Wrote {n_written:,} rows in {elapsed:.1f}s "
            f"({n_written / max(elapsed, 0.001):.0f}/s).",
            flush=True,
        )
        print(f"DB size: {db_path.stat().st_size / 1e9:.2f} GB")
        return n_written
    finally:
        conn.close()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--json",
        type=Path,
        default=DEFAULT_JSON,
        help="path to source JSON cache (default: dex_cache.json)",
    )
    parser.add_argument(
        "--db",
        type=Path,
        default=DEFAULT_DB,
        help="path to destination SQLite DB (default: dex_cache.db)",
    )
    args = parser.parse_args()
    return 0 if migrate(args.json, args.db) >= 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
