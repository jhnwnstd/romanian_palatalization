#!/usr/bin/env python3
"""Core DEX (Dicționarul Explicativ al Limbii Române) lookup utilities.

Exposes :func:`polite_get` and the persistent ``DISK_CACHE`` k/v store
used by the harvester and the DEX-supplement scripts.

Cache architecture
------------------
- ``HTML_CACHE`` is a small in-process dict (L1). It satisfies repeat
  lookups within a single run without hitting disk.
- ``DISK_CACHE`` is a SQLite-backed k/v store (L2). It persists across
  runs. Every successful HTTP fetch writes through to L2, so worker
  threads see each other's writes immediately and crashes don't lose
  partial work.

The L2 store used to be a single 6.9 GB JSON file loaded into memory at
module-import time. That implementation has been retired; the on-disk
schema is now SQLite (``dex_cache.db``). Random lookups are O(log n)
and there is no module-import I/O cost.

Thread safety
-------------
- ``_CACHE_LOCK`` guards mutations of the in-memory ``HTML_CACHE``.
- The :class:`_SqliteKVStore` carries its own internal lock for its
  SQLite connection (single shared connection used across worker
  threads, per SQLite's threading model with ``check_same_thread=False``).
- The two locks protect disjoint resources, so they can't deadlock.
"""

import random
import re
import sqlite3
import threading
import time
import unicodedata
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from urllib.parse import quote

import requests
from bs4 import BeautifulSoup

# ------------------------------------------------------------
# Config (polite crawling)
# ------------------------------------------------------------
BASE = "https://dexonline.ro"
UA = (
    "RomanianLexiconQC/2.0 "
    "(+academic research on Romanian palatalization; "
    "contact via github.com/anthropics/claude-code)"
)
HEADERS = {"User-Agent": UA}
THROTTLE = (0.05, 0.1)  # jittered sleep per request (min, max)
RETRIES = 4
TIMEOUT = 10

KEEP_POS = {"s.m.", "s.f.", "s.n.", "adj.", "vb."}
DASHES = {"-", "—", "–"}

# Accept simple plural tokens (or dash for none)
PL_OK = re.compile(r"^([-—–]|[A-Za-zăâîșțĂÂÎȘȚ\- ()]+)$")

PAREN_TAIL = re.compile(
    r"\s*\([^)]*\)\s*$"
)  # strip parenthetical notes at end

# DEX spaced-stress pattern: "CONF O RT" -> "CONFORT"
_STRESS_SPACE_RE = re.compile(r"\b([A-ZĂÂÎȘȚ]+) ([A-ZĂÂÎȘȚ]) ([A-ZĂÂÎȘȚ]+)\b")

# Plural extraction from DEX definition text:
#   "Pl. conforturi", "pl. oportunități", etc.
_DEX_PL_RE = re.compile(
    r"\bpl\.\s+([a-zăâîșț][a-zăâîșț\-]+)",
    re.I,
)

# Short POS forms used in older DEX entries: "f.", "m.", "n."
_SHORT_POS_RE = re.compile(
    r"(?<![a-zA-Z])\b([fmn])\.\s",
)
_SHORT_POS_MAP = {"f": "s.f.", "m": "s.m.", "n": "s.n."}


# ------------------------------------------------------------
# Helpers / normalization
# ------------------------------------------------------------
def _nfc(s: str) -> str:
    """Normalize to Unicode NFC."""
    return unicodedata.normalize("NFC", s or "")


def _norm_ws(s: str) -> str:
    """
    Normalize whitespace: replace nbsp/narrow-nbsp with regular space,
    collapse runs.
    """
    s = s.replace("\u00a0", " ").replace("\u202f", " ")
    return re.sub(r"\s+", " ", s)


def norm_lemma(s: str) -> str:
    """Normalize lemma: NFC, strip punctuation/whitespace."""
    s = _nfc(s or "").strip()
    s = s.strip(".,;:!?)(")
    s = _norm_ws(s)
    return s


def _clean_plural_token(pl: str) -> str:
    """Normalize plural token; return '-' if invalid."""
    if not pl or pl in DASHES:
        return "-"
    x = norm_lemma(pl)
    if not PL_OK.match(x):
        return "-"
    return x


def _classify_plural(pl: str) -> str:
    """Classify plural as 'i', 'e', 'uri', 'none', or 'unknown'.

    Returns:
        'i': plural ends in -i or -ii
        'e': plural ends in -e
        'uri': plural ends in -uri
        'none': explicitly no plural (dash markers)
        'unknown': cannot determine or ambiguous
    """
    if not pl:
        return "unknown"
    if pl in DASHES:
        return "none"
    x = PAREN_TAIL.sub("", pl).strip().lower()
    if x.endswith("uri"):
        return "uri"
    if x.endswith(("ii", "i")):
        return "i"
    if x.endswith("e"):
        return "e"
    return "unknown"


def _gender_from_pos(pos_raw: str) -> str:
    """Map s.m./s.f./s.n. to MASC/FEM/NEUT."""
    return {"s.m.": "MASC", "s.f.": "FEM", "s.n.": "NEUT"}.get(pos_raw, "")


# ------------------------------------------------------------
# POS normalization & abbr scanning
# ------------------------------------------------------------
POS_CANON = {
    "sm": "s.m.",
    "s.m": "s.m.",
    "s.m.": "s.m.",
    "sf": "s.f.",
    "s.f": "s.f.",
    "s.f.": "s.f.",
    "sn": "s.n.",
    "s.n": "s.n.",
    "s.n.": "s.n.",
    "adj.": "adj.",
    "vb.": "vb.",
}


def _canonize_pos_token(raw: str) -> str:
    """Canonize POS abbreviation."""
    token = raw or ""
    token = _norm_ws(token)
    token = re.sub(r"\s+", "", token.lower())
    if token in POS_CANON:
        return POS_CANON[token]
    if token in {"sm.", "sf.", "sn."}:
        return POS_CANON.get(token[:-1], "")
    return ""


def _scan_abbr_pos_anywhere(soup: "BeautifulSoup") -> Optional[str]:
    """Scan for <abbr> tags containing noun POS."""
    for ab in soup.find_all("abbr"):
        txt = ab.get_text(" ", strip=True) or ""
        pos = _canonize_pos_token(txt)
        if pos in {"s.m.", "s.f.", "s.n."}:
            return pos
        title = ab.get("title") or ""
        if isinstance(title, str):
            low = title.strip().lower()
            if "substantiv masculin" in low:
                return "s.m."
            if "substantiv feminin" in low:
                return "s.f."
            if "substantiv neutru" in low or "substantiv neuter" in low:
                return "s.n."
    return ""


def _scan_pos_in_text(soup: "BeautifulSoup") -> Optional[str]:
    """Scan for s.m./s.f./s.n. in full text."""
    blob = _norm_ws(soup.get_text(" ", strip=True))
    match = re.search(r"\bs\.\s*[mfn]\.\b", blob, flags=re.I)
    if match:
        return _canonize_pos_token(match.group(0))
    return ""


# ------------------------------------------------------------
# Flexible header matcher (DEX)
# ------------------------------------------------------------
# DEX header: "LEMMA, PLURAL, s.m." — lemma is ALL-CAPS (after
# collapsing stress spaces).  The plural is a single lowercase
# Romanian word, not a sentence fragment.
HEADER_FLEX_RE = re.compile(
    r"""
    ^\s*
    (?P<lemma>[A-ZĂÂÎȘȚ][A-ZĂÂÎȘȚ\- ]{0,30})
    \s*,\s*
    (?:
        -[A-ZĂÂÎȘȚa-zăâîșț]*
        \s*,\s*
    )?
    (?P<plural>[a-zăâîșț][a-zăâîșț\- ]{0,30}|[-—–])
    \s*,\s*
    (?P<posblob>(?:s\.\s*[mfn]\.|adj\.|vb\.)\s*.*)
    \s*$
    """,
    re.X,
)
HEADER_POS_RE = re.compile(
    r"^\s*([A-ZĂÂÎȘȚ][A-ZĂÂÎȘȚ\- ]{0,30})\s*,\s*(adj\.|vb\.)\s*$"
)


def _collapse_stress_spaces(text: str) -> str:
    """Collapse DEX spaced-stress marks: 'CONF O RT' -> 'CONFORT'."""
    return _STRESS_SPACE_RE.sub(
        lambda m: m.group(1) + m.group(2) + m.group(3), text
    )


def _pick_noun_pos_from_blob(posblob: str) -> Optional[str]:
    """Extract noun POS from posblob string."""
    if not posblob:
        return None
    blob = _norm_ws(posblob.lower())
    for tag in ("s. m.", "s.m.", "s. f.", "s.f.", "s. n.", "s.n."):
        if tag in blob:
            return tag.replace(" ", "")
    # Short forms: standalone "f.", "m.", "n." at start of POS blob
    m = _SHORT_POS_RE.search(blob)
    if m:
        return _SHORT_POS_MAP.get(m.group(1))
    return None


# ------------------------------------------------------------
# Data model (DEX)
# ------------------------------------------------------------
@dataclass
class DexEntry:  # pylint: disable=too-many-instance-attributes
    """Structured data extracted from DEX."""

    lemma: str
    pos_raw: str
    plural: str
    gender: str
    url: str
    source: str = "DEX"
    dex_confidence: str = "header"  # 'header' | 'pos_only'
    plural_class: str = ""
    notes: str = ""


# ------------------------------------------------------------
# Session + persistent disk cache
# ------------------------------------------------------------
SESSION = requests.Session()
SESSION.headers.update(HEADERS)
HTML_CACHE: Dict[str, str] = {}  # In-memory L1 cache for this run


class _SqliteKVStore:
    """Dict-like view of a SQLite k/v cache.

    Implements ``__contains__``, ``__getitem__``, ``__setitem__``, and
    ``__len__`` so existing callers that treated ``DISK_CACHE`` as a
    dict keep working unchanged.

    Reads and writes go straight to the on-disk DB — there is no
    in-memory mirror at this layer. Repeated lookups within a run
    should hit the in-process ``HTML_CACHE`` (L1) instead.

    A single shared SQLite connection is reused across worker threads
    (``check_same_thread=False``) and guarded by an internal lock so
    concurrent writes serialize correctly.
    """

    __slots__ = ("_path", "_conn", "_lock")

    def __init__(self, path: Path) -> None:
        self._path = path
        self._lock = threading.Lock()
        self._conn = sqlite3.connect(
            str(path), check_same_thread=False
        )
        # Default DELETE journal mode keeps the per-transaction
        # journal small; commits are durable per-row.
        self._conn.execute(
            "CREATE TABLE IF NOT EXISTS cache "
            "(url TEXT PRIMARY KEY, html TEXT NOT NULL)"
        )
        self._conn.commit()

    def __contains__(self, url: object) -> bool:
        if not isinstance(url, str):
            return False
        with self._lock:
            row = self._conn.execute(
                "SELECT 1 FROM cache WHERE url = ? LIMIT 1", (url,)
            ).fetchone()
        return row is not None

    def __getitem__(self, url: str) -> str:
        with self._lock:
            row = self._conn.execute(
                "SELECT html FROM cache WHERE url = ?", (url,)
            ).fetchone()
        if row is None:
            raise KeyError(url)
        return str(row[0])

    def __setitem__(self, url: str, html: str) -> None:
        with self._lock:
            self._conn.execute(
                "INSERT OR REPLACE INTO cache "
                "(url, html) VALUES (?, ?)",
                (url, html),
            )
            self._conn.commit()

    def __len__(self) -> int:
        with self._lock:
            row = self._conn.execute(
                "SELECT COUNT(*) FROM cache"
            ).fetchone()
        return int(row[0])


# Persistent disk cache (SQLite, in repo root).
DISK_CACHE_PATH = Path(__file__).resolve().parent.parent / "dex_cache.db"
DISK_CACHE: _SqliteKVStore = _SqliteKVStore(DISK_CACHE_PATH)
# Guards mutations of HTML_CACHE only. DISK_CACHE has its own lock.
_CACHE_LOCK = threading.Lock()


def load_disk_cache() -> None:
    """No-op retained for backward compatibility.

    With the SQLite backend, the on-disk cache is queried lazily per
    lookup. Nothing needs to be loaded into memory at module import.
    """


def save_disk_cache() -> None:
    """No-op retained for backward compatibility.

    Every cache write goes through ``polite_get`` → ``DISK_CACHE``
    setitem, which commits to SQLite per row. There is nothing to
    flush.
    """


def _cached_response(url: str, html: str) -> requests.Response:
    """Build a fake Response from cached HTML."""
    resp = requests.Response()
    resp.status_code = 200
    resp._content = html.encode("utf-8")  # noqa: SLF001
    resp.url = url
    resp.encoding = "utf-8"
    return resp


def polite_get(url: str) -> requests.Response:
    """GET with retries, rate limiting, and persistent caching.

    L1 (in-memory ``HTML_CACHE``) is checked first under
    ``_CACHE_LOCK``. L2 (SQLite ``DISK_CACHE``) is checked next; on
    L2 hit the page is promoted to L1 so subsequent requests within
    the same run avoid the SQLite round-trip.

    On HTTP success the page is written to both L1 and L2. The L2
    write is durable (one INSERT-or-REPLACE commit), so a crash
    mid-batch does not lose successful fetches.
    """
    with _CACHE_LOCK:
        if url in HTML_CACHE:
            return _cached_response(url, HTML_CACHE[url])
    if url in DISK_CACHE:
        html = DISK_CACHE[url]
        with _CACHE_LOCK:
            HTML_CACHE[url] = html
        return _cached_response(url, html)
    last_exc: Optional[requests.RequestException] = None
    for attempt in range(1, RETRIES + 1):
        try:
            response = SESSION.get(url, timeout=TIMEOUT)
            if response.status_code in (429, 503):
                time.sleep(min(10, attempt * 1.5))
                continue
            response.raise_for_status()
            text = response.text
            with _CACHE_LOCK:
                HTML_CACHE[url] = text
            DISK_CACHE[url] = text
            time.sleep(random.uniform(*THROTTLE))
            return response
        except requests.RequestException as exc:
            last_exc = exc
            time.sleep(min(10, attempt * 1.2))
    if last_exc is None:
        raise RuntimeError(
            f"Failed to retrieve {url} after {RETRIES} attempts"
        )
    raise last_exc


# ------------------------------------------------------------
# DEX parsing
# ------------------------------------------------------------
def dex_url_for(lemma: str) -> str:
    """Build DEX URL for a given lemma."""
    lemma_norm = norm_lemma(_nfc(lemma))
    return f"{BASE}/definitie/{quote(lemma_norm)}/definitii"


# DEX inline header: "LEMMA s. n. definition..." or "LEMMA f. def..."
_DEX_INLINE_RE = re.compile(
    r"""
    ^\s*
    (?P<lemma>[A-ZĂÂÎȘȚ][A-ZĂÂÎȘȚa-zăâîșț\-]*?)
    \s+
    (?P<posblob>(?:s\.\s*[mfn]\.|[fmn]\.)\s*)
    (?P<rest>.*)
    $
    """,
    re.X,
)


def _try_parse_flex_header(
    line: str, seen_headers: Set[Tuple[str, str, str]]
) -> Optional[DexEntry]:
    """Try to parse a flex header (lemma + plural + POS)."""
    collapsed = _collapse_stress_spaces(line)
    match = HEADER_FLEX_RE.match(collapsed)
    if not match:
        return None
    lemma = norm_lemma(match.group("lemma"))
    plural = _clean_plural_token(match.group("plural"))
    posblob = _nfc(match.group("posblob") or "").strip()
    noun_pos = _pick_noun_pos_from_blob(posblob)
    if not noun_pos or noun_pos not in KEEP_POS:
        return None
    key = (lemma, plural, noun_pos)
    if key in seen_headers:
        return None
    seen_headers.add(key)
    return DexEntry(
        lemma=lemma,
        pos_raw=noun_pos,
        plural=plural,
        gender=_gender_from_pos(noun_pos),
        url="",
        dex_confidence="header",
        plural_class=_classify_plural(plural),
        notes="parsed_header_flex",
    )


def _try_parse_inline_header(
    line: str, seen_headers: Set[Tuple[str, str, str]]
) -> Optional[DexEntry]:
    """Parse DEX inline format: 'CONFORT s. n. Definition...'

    Also extracts plural from 'Pl. ...' in the definition text.
    """
    collapsed = _collapse_stress_spaces(line)
    match = _DEX_INLINE_RE.match(collapsed)
    if not match:
        return None
    lemma = norm_lemma(match.group("lemma"))
    posblob = match.group("posblob").strip()
    rest = match.group("rest") or ""

    noun_pos = _pick_noun_pos_from_blob(posblob)
    if not noun_pos or noun_pos not in KEEP_POS:
        return None

    # Try to extract plural from definition text
    plural = "-"
    pl_match = _DEX_PL_STRICT_RE.search(";" + rest)
    if pl_match:
        plural = _clean_plural_token(pl_match.group(1))
    if plural == "-":
        pl_match = _DEX_PL_SPACED_RE.search(";" + rest)
        if pl_match:
            collapsed = _collapse_spaced_plural(pl_match.group(1))
            plural = _clean_plural_token(collapsed)

    key = (lemma, plural, noun_pos)
    if key in seen_headers:
        return None
    seen_headers.add(key)

    confidence = "header" if plural != "-" else "inline_pos"
    return DexEntry(
        lemma=lemma,
        pos_raw=noun_pos,
        plural=plural,
        gender=_gender_from_pos(noun_pos),
        url="",
        dex_confidence=confidence,
        plural_class=_classify_plural(plural),
        notes="parsed_inline",
    )


def _try_parse_pos_only_header(line: str) -> Optional[DexEntry]:
    """Try to parse a POS-only header (lemma + POS)."""
    collapsed = _collapse_stress_spaces(line)
    match = HEADER_POS_RE.match(collapsed)
    if not match:
        return None
    lemma, pos = match.groups()
    lemma = norm_lemma(lemma)
    pos = _nfc(pos).strip()
    if pos not in KEEP_POS:
        return None
    return DexEntry(
        lemma=lemma,
        pos_raw=pos,
        plural="-",
        gender=_gender_from_pos(pos),
        url="",
        dex_confidence="pos_only",
        plural_class="unknown",
        notes="parsed_pos_only",
    )


# Strict plural extraction: match "pl. WORD" only when WORD
# looks like a Romanian noun (all lowercase Romanian letters,
# at least 3 chars, no spaces after collapsing stress marks).
_DEX_PL_STRICT_RE = re.compile(
    r";\s*pl\.\s+" r"([a-zăâîșț][a-zăâîșț]{2,})" r"(?:\s|;|,|\.|$)",
)

# Spaced-stress variant: "pl. oportunit ă ți"
# Single-char tokens between word parts are stressed vowels.
_DEX_PL_SPACED_RE = re.compile(
    r";\s*pl\.\s+"
    r"([a-zăâîșț]+"
    r"(?:\s[a-zăâîșț]\s[a-zăâîșț]+)*)"
    r"(?:\s|;|,|\.|$)",
)


def _collapse_spaced_plural(raw: str) -> str:
    """Collapse 'oportunit ă ți' -> 'oportunități'."""
    parts = raw.split()
    out: List[str] = []
    for p in parts:
        if len(p) == 1 and p.isalpha():
            out.append(p)
        else:
            out.append(p)
    return "".join(out)


def _extract_plural_from_page(
    soup: "BeautifulSoup",
    lemma: str,
) -> str:
    """Try to extract a plural from the full page text.

    Looks for ``; pl. WORD`` patterns in DEX definition headers,
    handling DEX's spaced-stress format.
    Returns the plural token or ``"-"`` if not found.
    """
    blob = _norm_ws(soup.get_text(" ", strip=True))
    lemma_lower = (lemma or "").lower()

    # Try strict pattern first (no spaced stress)
    for m in _DEX_PL_STRICT_RE.finditer(blob):
        candidate = _clean_plural_token(m.group(1))
        if (
            candidate != "-"
            and candidate.lower() != lemma_lower
            and len(candidate) >= 3
        ):
            return candidate

    # Try spaced-stress pattern
    for m in _DEX_PL_SPACED_RE.finditer(blob):
        collapsed = _collapse_spaced_plural(m.group(1).strip())
        candidate = _clean_plural_token(collapsed)
        if (
            candidate != "-"
            and candidate.lower() != lemma_lower
            and len(candidate) >= 3
        ):
            return candidate

    return "-"


def _parse_header_lines(
    soup: "BeautifulSoup",
) -> List[DexEntry]:
    """
    Parse DEX header lines from BeautifulSoup object into DexEntry
    objects.
    """
    entries: List[DexEntry] = []
    candidates: List[str] = []
    for tag in soup.find_all(["p", "div", "li", "span"]):
        txt = (tag.get_text(" ", strip=True) or "").strip()
        if 2 <= len(txt) <= 220:
            candidates.append(txt)
    seen_headers: Set[Tuple[str, str, str]] = set()
    for line in candidates:
        entry = _try_parse_flex_header(line, seen_headers)
        if entry:
            entries.append(entry)
            continue
        # Try DEX inline format ("CONFORT s. n. definition...")
        entry = _try_parse_inline_header(line, seen_headers)
        if entry:
            entries.append(entry)
            continue
        entry = _try_parse_pos_only_header(line)
        if entry:
            entries.append(entry)

    noun_have = [e for e in entries if e.pos_raw in {"s.m.", "s.f.", "s.n."}]
    if noun_have:
        # If we found noun headers but none have a real plural,
        # try to extract one from the page text.
        has_real_plural = any(
            e.plural and e.plural not in DASHES for e in noun_have
        )
        if not has_real_plural:
            # Use the first noun header's lemma for context
            page_lemma = noun_have[0].lemma or ""
            page_plural = _extract_plural_from_page(soup, page_lemma)
            if page_plural != "-":
                noun_have[0].plural = page_plural
                noun_have[0].plural_class = _classify_plural(page_plural)
                noun_have[0].dex_confidence = "header"
                noun_have[0].notes = "plural_from_page_text"
        return entries

    pos_canon = _scan_abbr_pos_anywhere(soup) or _scan_pos_in_text(soup)
    if pos_canon in {"s.m.", "s.f.", "s.n."}:
        # Try to get plural from page text
        page_plural = _extract_plural_from_page(soup, "")
        confidence = "header" if page_plural != "-" else "pos_only"
        entries.append(
            DexEntry(
                lemma="",
                pos_raw=pos_canon,
                plural=page_plural,
                gender=_gender_from_pos(pos_canon),
                url="",
                dex_confidence=confidence,
                plural_class=_classify_plural(page_plural),
                notes="abbr_fallback",
            )
        )
        return entries
    return []


def fetch_dex_page(lemma: str) -> Tuple[List[DexEntry], str]:
    """
    Fetch DEX page for lemma and parse headers.

    Returns:
        (headers, html): List of DexEntry objects and raw HTML.
    """
    url = dex_url_for(lemma)
    response = polite_get(url)
    soup = BeautifulSoup(response.text, "html.parser")
    headers = _parse_header_lines(soup)
    for header in headers:
        header.url = url
        if not header.lemma:
            header.lemma = norm_lemma(lemma)
    return headers, response.text


def best_dex_noun_header(
    headers: List[DexEntry], lemma: str
) -> Optional[DexEntry]:
    """
    Pick the best noun header for a given lemma.

    Prefers:
    1. Exact case-insensitive match with real plural
       that looks like a plausible inflection of the lemma
    2. Exact case-insensitive match (even if no plural)
    3. Any noun header with real plural
    4. First noun header
    """
    if not headers:
        return None
    lemma_norm_lower = norm_lemma(lemma).lower()
    noun_headers = [
        h for h in headers if h.pos_raw in {"s.m.", "s.f.", "s.n."}
    ]
    if not noun_headers:
        return None

    def _plural_looks_plausible(h: DexEntry) -> bool:
        """Reject plurals that are clearly from a different word.

        E.g., 'nord' getting plural 'nord-americance' from a
        compound entry, or 'engleză' getting '-Ă adj'.
        """
        pl = (h.plural or "").lower().strip()
        if not pl or pl in DASHES:
            return False
        # Reject if plural contains markers of wrong entry
        if " " in pl or "adj" in pl:
            return False
        # Reject if plural doesn't share any prefix with lemma
        # (a plausible inflection should overlap significantly)
        lem = lemma_norm_lower
        shared = 0
        for a, b in zip(lem, pl):
            if a == b:
                shared += 1
            else:
                break
        # Require at least 2 shared chars or 40% of lemma length
        min_shared = max(2, len(lem) * 2 // 5)
        return shared >= min_shared

    # Filter by exact match
    exact = [
        h
        for h in noun_headers
        if norm_lemma(h.lemma).lower() == lemma_norm_lower
    ]
    if exact:
        # Prefer exact matches with plausible plurals
        exact_plausible = [h for h in exact if _plural_looks_plausible(h)]
        if exact_plausible:
            return exact_plausible[0]
        # Exact match with any real plural
        exact_with_plural = [
            h for h in exact if h.plural and h.plural not in DASHES
        ]
        if exact_with_plural:
            return exact_with_plural[0]
        return exact[0]

    # No exact match — prefer plausible plurals only
    plausible = [h for h in noun_headers if _plural_looks_plausible(h)]
    if plausible:
        return plausible[0]

    # No plausible match.  Return None rather than a header from
    # an unrelated word (e.g., ȘUVAR on the page for "munte").
    return None


def dex_has_entry(lemma: str, restrict_to_keep_pos: bool = False) -> bool:
    """
    Lightweight existence check: does DEX have *any* definition for this
    lemma?

    Uses two signals:
      1. _parse_header_lines() producing any entries (works for nouns,
         adjectives via the structured header parsers)
      2. the raw HTML containing definition markers ("sursa:"), which
         is present for verbs and other POS that the header parsers
         do not currently recognize

    If restrict_to_keep_pos is True, requires that at least one parsed
    header has pos_raw in KEEP_POS (s.m., s.f., s.n., adj., vb.). In
    that mode the HTML-level fallback is skipped because we cannot
    distinguish POS without the parse.
    """
    lemma_norm = norm_lemma(lemma)
    if not lemma_norm:
        return False
    try:
        headers, html = fetch_dex_page(lemma_norm)
    except Exception:
        # On network / parsing error, treat as "unknown" – caller can decide
        return False

    if headers:
        if not restrict_to_keep_pos:
            return True
        return any(h.pos_raw in KEEP_POS for h in headers)

    if restrict_to_keep_pos:
        return False

    # No parsed headers: fall back to an HTML-level check. DEX pages
    # for real entries always contain "sursa:" markers (source
    # attributions for definitions). Nonexistent lemma pages have
    # zero such markers. This catches verb entries (e.g., "acidifica",
    # "amurgi") that the header parsers currently skip because the
    # flex/inline parsers only emit noun-POS results.
    return "sursa:" in html
