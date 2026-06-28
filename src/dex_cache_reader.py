"""Streaming reader for the on-disk DEX HTML cache.

The cache file is ~7 GB JSON keyed by full DEX URL. Loading it
into memory blows past the available RAM on a typical workstation,
so this module exposes a streaming API that scans the JSON once,
pulls only the keys our analysis asks for, and returns a small
mapping the rest of the code can use synchronously.

The URL schema mirrors what ``polite_get`` in ``dex_utils`` stores,
so cached pages remain shareable across tools.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable
from urllib.parse import quote

import ijson

# Default cache location. Tools that need a different path should pass
# it explicitly; this constant captures the conventional one.
DEFAULT_CACHE_PATH: Path = (
    Path(__file__).resolve().parent.parent / "dex_cache.json"
)


def url_variants(form: str) -> tuple[str, ...]:
    """All cache-key variants we might have stored a DEX page under.

    The cache writer historically used both URL-encoded and bare-form
    keys, and both the ``/definitii`` and bare paths. Any of them
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
    """Stream the cache file and return only the URLs we need.

    ``forms`` is the set of DEX surface lemmas we want pages for.
    Returns a dict keyed by URL (so callers can choose which variant
    matched). The pass is single-shot: the cache is opened in binary
    mode and traversed once via ijson.
    """
    needed: set[str] = set()
    for form in forms:
        needed.update(url_variants(form))
    matched: dict[str, str] = {}
    with cache_path.open("rb") as f:
        for key, value in ijson.kvitems(f, ""):
            if key in needed:
                matched[key] = value
    return matched


def get_html(cache: dict[str, str], form: str) -> str | None:
    """Return the cached HTML for ``form`` if any URL variant matched.

    Tries variants in the same order as ``url_variants``, returning the
    first hit. Returns ``None`` if no variant is in the cache.
    """
    for url in url_variants(form):
        html = cache.get(url)
        if html is not None:
            return html
    return None
