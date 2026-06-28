"""Random-access reader for the on-disk DEX HTML cache.

The cache lives in a SQLite k/v store (``dex_cache.db``, schema
``cache(url TEXT PRIMARY KEY, html TEXT)``). Lookups are O(log n) per
URL — independent of cache size — so analysis scripts that need a few
hundred pages don't pay a multi-GB scan.

Historical note
---------------
Previously the cache was a single 6.9 GB JSON file. Loading or
streaming it cost ~30 s per analysis run even when we only needed a
few hundred lookups. See ``scripts/migrate_dex_cache_to_sqlite.py`` for
the migration.

API
---
- :func:`url_variants` — every cache-key variant we might have stored
  a DEX page under.
- :func:`load_cache_subset` — return a small in-memory dict for the
  forms an analysis run will touch. Wraps the SQLite query in the
  same shape that the previous streaming-JSON loader returned, so
  callers don't need to change.
- :func:`get_html` — look up a single form against an existing subset.

The functions are read-only. Writes to the cache go through
``polite_get`` in :mod:`dex_utils`.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Iterable
from urllib.parse import quote

DEFAULT_CACHE_PATH: Path = (
    Path(__file__).resolve().parent.parent / "dex_cache.db"
)


def url_variants(form: str) -> tuple[str, ...]:
    """All cache-key variants we might have stored a DEX page under.

    The cache writer historically used both URL-encoded and bare-form
    keys, with and without the ``/definitii`` suffix. Any of the four
    counts as a cache hit for ``form``.
    """
    q = quote(form)
    return (
        f"https://dexonline.ro/definitie/{q}/definitii",
        f"https://dexonline.ro/definitie/{q}",
        f"https://dexonline.ro/definitie/{form}/definitii",
        f"https://dexonline.ro/definitie/{form}",
    )


def load_cache_subset(
    forms: Iterable[str],
    cache_path: Path = DEFAULT_CACHE_PATH,
) -> dict[str, str]:
    """Return the cache pages relevant to ``forms`` as a ``{url: html}``.

    Issues one prepared SELECT per URL variant. Total time is
    proportional to the number of forms, not the cache size.
    """
    needed: list[str] = []
    for form in forms:
        needed.extend(url_variants(form))
    out: dict[str, str] = {}
    with sqlite3.connect(f"file:{cache_path}?mode=ro", uri=True) as conn:
        conn.execute("PRAGMA query_only = ON")
        # Chunk the IN(...) clause to avoid SQLite's parameter limit
        # (default 999 for older builds, 32766 for newer). 500 is a
        # safe round number.
        CHUNK = 500
        cursor = conn.cursor()
        for i in range(0, len(needed), CHUNK):
            chunk = needed[i : i + CHUNK]
            placeholders = ",".join("?" * len(chunk))
            cursor.execute(
                f"SELECT url, html FROM cache WHERE url IN ({placeholders})",
                chunk,
            )
            for url, html in cursor.fetchall():
                # SQLite returns Any-typed rows; narrow at the trust
                # boundary so callers get real str.
                out[str(url)] = str(html)
    return out


def get_html(cache: dict[str, str], form: str) -> str | None:
    """Look up ``form`` in an already-loaded cache subset.

    Tries variants in the order returned by :func:`url_variants`,
    returning the first hit. Returns ``None`` if no variant is present.
    """
    for url in url_variants(form):
        html = cache.get(url)
        if html is not None:
            return html
    return None
