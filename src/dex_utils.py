#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Core DEX (Dicționarul Explicativ al Limbii Române) lookup utilities.
Extracted from dex_lookup.py to be reusable by QC scripts.
"""

import json
import random
import re
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

TARGET_FINALS = set("cgtdsz")
KEEP_POS = {"s.m.", "s.f.", "s.n.", "adj.", "vb."}
DASHES = {"-", "—", "–"}

# Accept simple plural tokens (or dash for none)
PL_OK = re.compile(r"^([-—–]|[A-Za-zăâîșțĂÂÎȘȚ\- ()]+)$")

PAREN_TAIL = re.compile(
    r"\s*\([^)]*\)\s*$"
)  # strip parenthetical notes at end


# ------------------------------------------------------------
# Helpers / normalization
# ------------------------------------------------------------
def nfc(s: str) -> str:
    """Normalize to Unicode NFC."""
    return unicodedata.normalize("NFC", s or "")


def norm_ws(s: str) -> str:
    """
    Normalize whitespace: replace nbsp/narrow-nbsp with regular space,
    collapse runs.
    """
    s = s.replace("\u00a0", " ").replace("\u202f", " ")
    return re.sub(r"\s+", " ", s)


def norm_lemma(s: str) -> str:
    """Normalize lemma: NFC, strip punctuation/whitespace."""
    s = nfc(s or "").strip()
    s = s.strip(".,;:!?)(")
    s = norm_ws(s)
    return s


def looks_like_target(lemma: str) -> bool:
    """
    True if lemma ends in one of TARGET_FINALS (c/g/t/d/s/z) and is
    single-token.
    """
    lemma_lower = norm_lemma(lemma).lower()
    if not lemma_lower or " " in lemma_lower:
        return False
    return lemma_lower[-1] in TARGET_FINALS


def clean_plural_token(pl: str) -> str:
    """Normalize plural token; return '-' if invalid."""
    if not pl or pl in DASHES:
        return "-"
    x = norm_lemma(pl)
    if not PL_OK.match(x):
        return "-"
    return x


def classify_plural(pl: str) -> str:
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


def gender_from_pos(pos_raw: str) -> str:
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


def canonize_pos_token(raw: str) -> str:
    """Canonize POS abbreviation."""
    token = raw or ""
    token = norm_ws(token)
    token = re.sub(r"\s+", "", token.lower())
    if token in POS_CANON:
        return POS_CANON[token]
    if token in {"sm.", "sf.", "sn."}:
        return POS_CANON.get(token[:-1], "")
    return ""


def scan_abbr_pos_anywhere(soup: "BeautifulSoup") -> Optional[str]:
    """Scan for <abbr> tags containing noun POS."""
    for ab in soup.find_all("abbr"):
        txt = ab.get_text(" ", strip=True) or ""
        pos = canonize_pos_token(txt)
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


def scan_pos_in_text(soup: "BeautifulSoup") -> Optional[str]:
    """Scan for s.m./s.f./s.n. in full text."""
    blob = norm_ws(soup.get_text(" ", strip=True))
    match = re.search(r"\bs\.\s*[mfn]\.\b", blob, flags=re.I)
    if match:
        return canonize_pos_token(match.group(0))
    return ""


# ------------------------------------------------------------
# Flexible header matcher (DEX)
# ------------------------------------------------------------
HEADER_FLEX_RE = re.compile(
    r"""
    ^\s*
    (?P<lemma>[A-ZĂÂÎȘȚ][A-ZĂÂÎȘȚa-zăâîșț\- ]*?)
    \s*,\s*
    (?:
        -[A-ZĂÂÎȘȚa-zăâîșț]*
        \s*,\s*
    )?
    (?P<plural>[^,]*?)
    \s*,\s*
    (?P<posblob>.+?)
    \s*$
    """,
    re.X,
)
HEADER_POS_RE = re.compile(
    r"^\s*([A-ZĂÂÎȘȚ][A-ZĂÂÎȘȚa-zăâîșț\- ]*?)\s*,\s*(adj\.|vb\.)\s*$"
)


def pick_noun_pos_from_blob(posblob: str) -> Optional[str]:
    """Extract noun POS from posblob string."""
    if not posblob:
        return None
    blob = norm_ws(posblob.lower())
    for tag in ("s. m.", "s.m.", "s. f.", "s.f.", "s. n.", "s.n."):
        if tag in blob:
            return tag.replace(" ", "")
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
HTML_CACHE: Dict[str, str] = {}  # In-memory cache for this run

# Persistent disk cache (in repo root)
DISK_CACHE_PATH = Path(__file__).resolve().parent.parent / "dex_cache.json"
DISK_CACHE: Dict[str, str] = {}
cache_dirty = False  # Track if we need to save


def load_disk_cache() -> None:
    """Load persistent cache from disk on startup."""
    global DISK_CACHE  # pylint: disable=global-statement
    if DISK_CACHE_PATH.exists():
        try:
            with open(DISK_CACHE_PATH, "r", encoding="utf-8") as f:
                DISK_CACHE = json.load(f)
            print(
                f"[dex_utils] Loaded {len(DISK_CACHE)} cached DEX "
                "pages from disk"
            )
        except (json.JSONDecodeError, OSError) as e:
            print(f"[dex_utils] Warning: Could not load cache: {e}")
            DISK_CACHE = {}
    else:
        DISK_CACHE = {}


def save_disk_cache() -> None:
    """Save persistent cache to disk."""
    global cache_dirty  # pylint: disable=global-statement
    if not cache_dirty:
        return
    try:
        with open(DISK_CACHE_PATH, "w", encoding="utf-8") as f:
            json.dump(DISK_CACHE, f, ensure_ascii=False, indent=2)
        cache_dirty = False
        print(f"[dex_utils] Saved {len(DISK_CACHE)} DEX pages to disk cache")
    except OSError as exc:
        print(f"[dex_utils] Warning: Could not save cache: {exc}")


# Load cache on module import
load_disk_cache()


def polite_get(url: str) -> requests.Response:
    """GET with retries, rate limiting, and persistent caching."""
    global cache_dirty  # pylint: disable=global-statement
    if url in HTML_CACHE:
        resp = requests.Response()
        resp.status_code = 200
        resp._content = HTML_CACHE[
            url
        ].encode(  # pylint: disable=protected-access
            "utf-8"
        )
        resp.url = url
        resp.encoding = "utf-8"
        return resp
    if url in DISK_CACHE:
        html = DISK_CACHE[url]
        HTML_CACHE[url] = html
        resp = requests.Response()
        resp.status_code = 200
        resp._content = html.encode(
            "utf-8"
        )  # pylint: disable=protected-access
        resp.url = url
        resp.encoding = "utf-8"
        return resp
    last_exc: Optional[requests.RequestException] = None
    for attempt in range(1, RETRIES + 1):
        try:
            response = SESSION.get(url, timeout=TIMEOUT)
            if response.status_code in (429, 503):
                time.sleep(min(10, attempt * 1.5))
                continue
            response.raise_for_status()
            HTML_CACHE[url] = response.text
            DISK_CACHE[url] = response.text
            cache_dirty = True
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
    lemma_norm = norm_lemma(nfc(lemma))
    return f"{BASE}/definitie/{quote(lemma_norm)}/definitii"


def _try_parse_flex_header(
    line: str, seen_headers: Set[Tuple[str, str, str]]
) -> Optional[DexEntry]:
    """Try to parse a flex header (lemma + plural + POS)."""
    match = HEADER_FLEX_RE.match(line)
    if not match:
        return None
    lemma = norm_lemma(match.group("lemma"))
    plural = clean_plural_token(match.group("plural"))
    posblob = nfc(match.group("posblob") or "").strip()
    noun_pos = pick_noun_pos_from_blob(posblob)
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
        gender=gender_from_pos(noun_pos),
        url="",
        dex_confidence="header",
        plural_class=classify_plural(plural),
        notes="parsed_header_flex",
    )


def _try_parse_pos_only_header(line: str) -> Optional[DexEntry]:
    """Try to parse a POS-only header (lemma + POS)."""
    match = HEADER_POS_RE.match(line)
    if not match:
        return None
    lemma, pos = match.groups()
    lemma = norm_lemma(lemma)
    pos = nfc(pos).strip()
    if pos not in KEEP_POS:
        return None
    return DexEntry(
        lemma=lemma,
        pos_raw=pos,
        plural="-",
        gender=gender_from_pos(pos),
        url="",
        dex_confidence="pos_only",
        plural_class="unknown",
        notes="parsed_pos_only",
    )


def parse_header_lines(
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
        entry = _try_parse_pos_only_header(line)
        if entry:
            entries.append(entry)
    noun_have = [e for e in entries if e.pos_raw in {"s.m.", "s.f.", "s.n."}]
    if noun_have:
        return entries
    pos_canon = scan_abbr_pos_anywhere(soup) or scan_pos_in_text(soup)
    if pos_canon in {"s.m.", "s.f.", "s.n."}:
        entries.append(
            DexEntry(
                lemma="",
                pos_raw=pos_canon,
                plural="-",
                gender=gender_from_pos(pos_canon),
                url="",
                dex_confidence="pos_only",
                plural_class="unknown",
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
    headers = parse_header_lines(soup)
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

    # Filter by exact match
    exact = [
        h
        for h in noun_headers
        if norm_lemma(h.lemma).lower() == lemma_norm_lower
    ]
    if exact:
        # Prefer exact matches with real plurals
        exact_with_plural = [
            h for h in exact if h.plural and h.plural not in DASHES
        ]
        if exact_with_plural:
            return exact_with_plural[0]
        return exact[0]

    # No exact match - prefer any header with real plural
    with_plural = [
        h for h in noun_headers if h.plural and h.plural not in DASHES
    ]
    if with_plural:
        return with_plural[0]

    return noun_headers[0]


def dex_has_entry(lemma: str, restrict_to_keep_pos: bool = False) -> bool:
    """
    Lightweight existence check: does DEX have *any* header for this lemma?

    If restrict_to_keep_pos is True, require that at least one header has
    pos_raw in KEEP_POS (s.m., s.f., s.n., adj., vb.).
    """
    lemma_norm = norm_lemma(lemma)
    if not lemma_norm:
        return False
    try:
        headers, _ = fetch_dex_page(lemma_norm)
    except Exception:
        # On network / parsing error, treat as "unknown" – caller can decide
        return False

    if not headers:
        return False
    if not restrict_to_keep_pos:
        return True
    return any(h.pos_raw in KEEP_POS for h in headers)
