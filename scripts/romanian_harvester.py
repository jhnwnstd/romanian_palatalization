#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Romanian Wiktionary Data Harvester - Extraction Only.

Extracts raw data from Romanian Wiktionary entries. Provides clean,
normalized extracted fields for downstream processing by
romanian_processor_main.py.

Output fields:
    - lemma, gloss, pos, gender, plural
    - derived_verbs, derived_adj
    - etym_lang, source, notes
    - ipa_raw_lemma, ipa_raw_pl

Usage:
    python romanian_harvester.py
"""

import json
import os
import random
import re
import sqlite3
import time
import unicodedata
from collections import defaultdict
from collections.abc import Mapping
from functools import lru_cache
from pathlib import Path
from typing import Any, DefaultDict, Optional, Set
from urllib.parse import quote, unquote

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)
from urllib3.util.retry import Retry

# Optional deps (safe fallbacks)
HAS_BS4: bool
try:
    from bs4 import BeautifulSoup

    HAS_BS4 = True
except ImportError:
    BeautifulSoup = None  # type: ignore
    HAS_BS4 = False

# ============================================================================
# CONFIGURATION
# ============================================================================

# API endpoints
WIKI_API_EN = "https://en.wiktionary.org/w/api.php"
WIKI_API_RO = "https://ro.wiktionary.org/w/api.php"
UA = (
    "RomanianLexicon/4.0 "
    "(https://example.edu/your-project-page; your.email@example.edu)"
)

# Discovery limits (EN Wiktionary has ~63k nouns, ~18k adj, ~7k verbs)
VERB_LIMIT = 20000
NOUN_LIMIT = 70000
ADJ_LIMIT = 25000
RO_NOUN_LIMIT = 30000
RO_ADJ_LIMIT = 15000
RO_VERB_LIMIT = 15000

# Rate limiting — Wikimedia allows ~200 req/s for identified bots.
# maxlag=5 handles server-side throttling; this is client-side politeness.
THROTTLE_DELAY = 0.05

# Output CSV
OUTPUT_CSV = Path(__file__).parent.parent / "data" / "romanian_lexicon_raw.csv"

# Disk cache files (in data directory)
_base_dir = Path(__file__).parent.parent / "data"
_base_dir.mkdir(parents=True, exist_ok=True)
CACHE_EN_PATH = str(_base_dir / "wt_cache_en.json")
CACHE_RO_PATH = str(_base_dir / "wt_cache_ro.json")
IPA_CACHE_PATH = str(_base_dir / "ipa_cache.json")
# HTML caches use SQLite for performance (the JSON versions were 2.7GB+
# and caused multi-minute load/save pauses and OOM kills)
HTML_DB_PATH = str(_base_dir / "html_cache.db")

# Enable IPA fetch from HTML
ENABLE_IPA_FETCH = True

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

ROMANIAN_CHARS = set("abcdefghijklmnopqrstuvwxyzăâîșț")

# IPA extraction patterns
IPA_SPAN_SELECTOR = "span.IPA, span.ipa, span.mw-IPA, span.mw-ipa"
IPA_TOK_RE = re.compile(r"[/\[][^\]/\]]{1,64}[/\]]")
IPA_HINT = re.compile(r"(?:ʃ|ʒ|ɨ|ə|ɡ|ɲ|ʎ|ɟ|ç|ʝ|t͡s|d͡z|t͡ʃ|d͡ʒ|ˈ|ˌ)")

# Wikitext section patterns
ROMANIAN_SECTION_RE = re.compile(
    r"^==\s*Romanian\s*==\s*$(.*?)(?=^==\s*[^=]|\Z)",
    re.MULTILINE | re.DOTALL | re.IGNORECASE,
)

# POS headings: support EN + RO labels
POS_HEAD_RE = re.compile(
    r"^===\s*(Noun|Substantiv|Adjective|Adjectiv|Verb)\s*===\s*$"
    r"(.*?)(?=^===|\Z)",
    re.MULTILINE | re.DOTALL | re.IGNORECASE,
)

# Normalize heading labels to our internal POS keys
POS_MAP = {
    # English headings (en.wiktionary)
    "noun": "N",
    "adjective": "ADJ",
    "verb": "V",
    # Romanian headings (ro.wiktionary)
    "substantiv": "N",
    "adjectiv": "ADJ",
}

PROPN_RE = re.compile(r"===\s*Proper noun\s*===", re.I)
GLOSS_LINE_RE = re.compile(r"^#\s*(.+)$", re.MULTILINE)
ETYM_RE = re.compile(
    r"^===\s*Etymology\s*(?:\d+)?\s*===\s*$\n(.*?)(?=^===|\Z)",
    re.MULTILINE | re.DOTALL,
)
IPA_WT_RE = re.compile(
    r"\{\{\s*IPA\s*\|\s*ro\s*\|\s*([^}|]+)",
    re.I,
)

# Gender normalization
GENDER_MAP = {
    "m": "MASC",
    "f": "FEM",
    "n": "NEUT",
    "masc": "MASC",
    "fem": "FEM",
    "neut": "NEUT",
    "masculine": "MASC",
    "feminine": "FEM",
    "neuter": "NEUT",
}

GENDER_TO_DECL = {"MASC": "i", "FEM": "e", "NEUT": "uri"}

# Plural extraction patterns
PLURAL_KEYS = [
    re.compile(r"\|\s*pl\s*=\s*([^\n\|\}]+)", re.I),
    re.compile(r"\|\s*plural\s*=\s*([^\n\|\}]+)", re.I),
]

TABLE_ANY_PL_RE = re.compile(
    r"\{\{\s*ro-noun-table[^}]*\|\s*(?:pl|plural)\s*=\s*([^\n\|}\]]+)",
    re.I,
)

# Etymology language tag extraction
ETYM_LANG_RE = re.compile(
    r"\{\{\s*(?:bor|der|inh|lbor|calque)\s*\|\s*ro\s*\|\s*"
    r"([a-z]{2,3})\s*(?:\||}})",
    re.I,
)

UNCOUNTABLE_RE = re.compile(
    r"\{\{\s*(?:unc|uncountable)\s*(?:\|[^}]*)?\}\}", re.I
)

# Canary examples
CANARY_LEMMAS = [
    "ac",
    "mic",
    "lung",
    "frate",
    "verde",
    "gros",
    "obraz",
    "arhitect",
    "artist",
    "banc",
]

# Steriade (2008) key examples
STERIADE_EXAMPLES = [
    # K-palatalization
    "nuc",
    "sărac",
    "mic",
    "lung",
    "olog",
    "bolnav",
    "alb",
    "mare",
    "foc",
    "loc",
    "tanc",
    "bloc",
    "catalog",
    "pedagog",
    "demolog",
    "teolog",
    "monolog",
    "psiholog",
    "geolog",
    "ideolog",
    "filolog",
    "biolog",
    "dialog",
    "patolog",
    "bogat",
    "curat",
    "artist",
    "romantic",
    "pitic",
    "critic",
    "politic",
    "domestic",
    "etic",
    "cinic",
    "fanatic",
    # Proper names, common nouns
    "Puică",
    "Volga",
    "Olga",
    "Iorga",
    "Coca",
    "Lunca",
    "Mureș",
    "Narcisa",
    "algă",
    "casă",
    "fată",
    "copac",
    "rege",
    "vulpe",
    "om",
    "frate",
    "lemn",
    "lup",
    "trup",
    "apă",
    "aripă",
    "lipsă",
    "nume",
    "lume",
    "vreme",
    "fiică",
    "mag",
    "face",
    "suge",
    "catarg",
    "stângă",
    "covrig",
    "săracă",
    # NDEB
    "ochi",
    "triunghi",
    "kilogram",
    "chestie",
    "ghem",
    "ghinion",
    # Denominal bases
    "pădure",
    "clește",
    "pildă",
    "păianjen",
    "drept",
    "grijă",
    "Franco",
    "Goga",
    "falangă",
    "logică",
    # Assibilation
    "verde",
    "supus",
    "viteaz",
    "gata",
    "lat",
    "aminte",
    "cuminte",
    "popas",
    "pas",
    "șase",
    "sănătos",
    "dovadă",
    "dovedi",
    "suflet",
    "însufleți",
    "perete",
    "pecețe",
    "miere",
    "cărare",
    "lene",
    "pășune",
    "fasole",
    "nimic",
    "vlagă",
    # Chitoran (2002) z~ʒ
    "obraz",
    "arbuz",
    "albigenz",
    "cartaginez",
    "albigenză",
    "aspergiloză",
]

# Enable string-based denominal verb heuristic
ENABLE_DENOMINAL_HEURISTIC = True

# Romanian vowels (for stripping final vowel in noun stems)
RO_VOWELS = set("aăâeiîouAĂÂEIÎOU")

# Highly transparent denominal verb suffixes of interest
DENOMINAL_VERB_SUFFIXES = ("i", "ui", "í", "uí")

# Derivational verbal prefixes that often form denominal verbs
DERIV_PREFIXES: tuple[str, ...] = (
    # Longer ones first for readability; we'll sort programmatically anyway
    "reîn",
    "între",
    "supra",
    "stră",
    "peste",
    "fără",
    "înd",
    "răs",
    "răz",
    "sub",
    "dez",
    "des",
    "pre",
    "ră",
    "în",
    "îm",
    "ne",
    "re",
    "de",
)

# For parsing we want longest-match-first
DERIV_PREFIXES_SORTED: tuple[str, ...] = tuple(
    sorted(DERIV_PREFIXES, key=len, reverse=True)
)


def _sleep_jitter(base_delay: float) -> None:
    """Sleep with ±20% jitter to avoid bursty request patterns."""
    jitter = base_delay * 0.2 * (2 * random.random() - 1)
    time.sleep(base_delay + jitter)


def _wiktionary_page_url(title: str, lang: str = "ro") -> str:
    """Build canonical Wiktionary page URL."""
    base = (
        "https://ro.wiktionary.org/wiki/"
        if lang == "ro"
        else "https://en.wiktionary.org/wiki/"
    )
    norm = unicodedata.normalize("NFC", title).replace(" ", "_")
    return base + quote(norm)


# ============================================================================
# HTTP SESSION & CACHE
# ============================================================================


def _build_session_with_retries() -> requests.Session:
    """Build Session with automatic retry on 5xx, connection errors."""
    session = requests.Session()
    retry_strategy = Retry(
        total=5,
        backoff_factor=1,
        status_forcelist=[500, 502, 503, 504],
        allowed_methods=frozenset({"GET", "POST"}),
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    session.headers.update({"User-Agent": UA})
    return session


_SESSION = _build_session_with_retries()
_WT_CACHE_EN: dict[str, str] = {}
_WT_CACHE_RO: dict[str, str] = {}
_IPA_CACHE: dict[str, list[str]] = {}


class _SqliteCache:
    """SQLite-backed key-value cache for HTML pages.

    Avoids loading multi-GB JSON into memory. Lookups and inserts
    are instant; no bulk load/save needed.
    """

    def __init__(self, db_path: str, table: str) -> None:
        self._table = table
        self._conn = sqlite3.connect(db_path)
        self._conn.execute("PRAGMA journal_mode=WAL")
        self._conn.execute("PRAGMA synchronous=NORMAL")
        self._conn.execute(
            f"CREATE TABLE IF NOT EXISTS {table} "
            f"(key TEXT PRIMARY KEY, value TEXT)"
        )
        self._conn.commit()

    def __contains__(self, key: str) -> bool:
        row = self._conn.execute(
            f"SELECT 1 FROM {self._table} WHERE key=? LIMIT 1",
            (key,),
        ).fetchone()
        return row is not None

    def __getitem__(self, key: str) -> str:
        row = self._conn.execute(
            f"SELECT value FROM {self._table} WHERE key=?",
            (key,),
        ).fetchone()
        if row is None:
            raise KeyError(key)
        # sqlite returns Any-typed rows; the value column is TEXT so
        # narrowing to str is correct here.
        return str(row[0])

    def __setitem__(self, key: str, value: str) -> None:
        self._conn.execute(
            f"INSERT OR REPLACE INTO {self._table} (key, value) "
            f"VALUES (?, ?)",
            (key, value),
        )
        self._conn.commit()

    def get(self, key: str, default: str = "") -> str:
        try:
            return self[key]
        except KeyError:
            return default

    def __len__(self) -> int:
        row = self._conn.execute(
            f"SELECT COUNT(*) FROM {self._table}"
        ).fetchone()
        return row[0] if row else 0

    def close(self) -> None:
        self._conn.close()


_HTML_CACHE_EN: _SqliteCache = None  # type: ignore[assignment]
_HTML_CACHE_RO: _SqliteCache = None  # type: ignore[assignment]

# Avoids repeated BeautifulSoup parsing
_PLURAL_TABLE_CACHE: dict[str, Optional[str]] = {}

DENOMINAL_VERBS: DefaultDict[str, Set[str]] = defaultdict(set)
DEADJECTIVAL_VERBS: DefaultDict[str, Set[str]] = defaultdict(set)
DENOMINAL_ADJS: DefaultDict[str, Set[str]] = defaultdict(
    set
)  # noun -> adjectives
DEADJECTIVAL_ADJS: DefaultDict[str, Set[str]] = defaultdict(
    set
)  # adj  -> adjectives


@retry(
    reraise=True,
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=0.5, min=0.5, max=5),
    retry=retry_if_exception_type(requests.RequestException),
)
def get(url: str, params: Mapping[str, Any]) -> dict[str, Any]:
    """Fetch JSON from Wiktionary API with rate limiting and retry."""
    params = {"maxlag": "5", **params}
    r = _SESSION.get(url, params=params, timeout=10)
    if r.status_code in (429, 503):
        try:
            retry_after = float(r.headers.get("Retry-After", "2"))
        except ValueError:
            retry_after = 2.0
        time.sleep(min(10.0, retry_after + 1.0))
        raise requests.RequestException(f"Rate-limited: {r.status_code}")
    r.raise_for_status()
    # r.json() is typed as Any; narrow at the trust boundary so the
    # caller doesn't have to deal with Any propagation.
    data: dict[str, Any] = r.json()
    if "error" in data and data["error"].get("code") == "maxlag":
        time.sleep(1.0 + random.random())
        raise requests.RequestException("Server under replication lag")
    return data


# ============================================================================
# TEXT NORMALIZATION
# ============================================================================


@lru_cache(maxsize=10000)
def normalize_unicode(s: str) -> str:
    """Normalize Unicode and convert legacy cedilla to comma-below."""
    if not s:
        return ""
    s = unicodedata.normalize("NFC", s)
    s = s.replace("Ş", "Ș").replace("ş", "ș")
    s = s.replace("Ţ", "Ț").replace("ţ", "ț")
    return s


def _normalize_ws(s: str) -> str:
    """Collapse whitespace to single spaces and strip edges."""
    return re.sub(r"\s+", " ", s or "").strip()


def _strip_wiki_markup(s: str) -> str:
    """Remove wiki markup (fast-path: skip regex if no markup chars)."""
    if not s:
        return s
    if not any(c in s for c in "[]{}&"):
        return s.strip()
    s = re.sub(r"\[\[([^\]|]+)(?:\|[^\]]+)?\]\]", r"\1", s)
    s = re.sub(r"\{\{[^\}]+\}\}", "", s)
    s = re.sub(r"&[a-z]+;", "", s)
    s = s.replace("[", "").replace("]", "")
    return s.strip()


def _clean_ipa_raw(ipa_str: str) -> str:
    """Remove IPA delimiters, HTML, percent-encoding; keep stress markers.

    Handles malformed HTML patterns like:
        encoded" title="IPA">IPA<
        c%C3%B2%27t%C9%99%CC%80" title="cò&#39;tə̀">cò'tə̀<
    """
    s = ipa_str.strip()
    if not s or s.lower() == "nan":
        return ""

    # Multiple rounds of URL decoding (some are double-encoded)
    try:
        for _ in range(3):
            decoded = unquote(s)
            if decoded == s:
                break
            s = decoded
    except (ImportError, ValueError, TypeError):
        pass

    # Extract from malformed HTML pattern: encoded" title="IPA">IPA<
    # Try to extract the cleanest version (usually in title or after >)
    title_match = re.search(r'title="([^"]*)"', s)
    content_match = re.search(r">([^<]+)<", s)

    if title_match:
        s = title_match.group(1)
    elif content_match:
        s = content_match.group(1)

    # Decode HTML entities
    s = re.sub(r"&#(\d+);", lambda m: chr(int(m.group(1))), s)
    s = re.sub(r"&#x([0-9a-fA-F]+);", lambda m: chr(int(m.group(1), 16)), s)

    # Remove remaining HTML
    s = re.sub(r"<[^>]+>", "", s)
    s = re.sub(r'\b(?:title|href|class|id|lang|span)="[^"]*"', "", s)
    s = re.sub(r"&[a-zA-Z]+;", " ", s)
    s = s.replace(">", "").replace("<", "").replace('"', "")
    s = s.replace("[", "").replace("]", "").replace("/", "")

    # Normalize whitespace
    s = re.sub(r"\s+", " ", s)
    s = re.sub(r"\s*\|\s*", " | ", s)

    # Deduplicate: sometimes HTML has the same IPA twice (in title and content)
    tokens = s.split()
    seen = set()
    deduped = []
    for tok in tokens:
        if tok not in seen:
            seen.add(tok)
            deduped.append(tok)
    s = " ".join(deduped)

    # Reject strings contaminated with foreign language labels
    _FOREIGN_MARKERS = (
        "français",
        "türkçe",
        "português",
        "español",
        "deutsch",
        "italiano",
        "english",
        "polski",
        "русский",
        "العربية",
        "中文",
        "日本語",
        "한국어",
        "wikipedia",
        "wiktionary",
    )
    s_lower = s.lower()
    if any(marker in s_lower for marker in _FOREIGN_MARKERS):
        return ""

    return s.strip()


def _clean_plural(s: str) -> str:
    """Strip markup, take first variant, reject placeholders."""
    if not s:
        return ""
    s = _strip_wiki_markup(s)
    s = normalize_unicode(s)
    first = re.split(r"\s*[,;/]|(?:\s+or\s+)|(?:\s+aka\s+)", s)[0]
    first = re.sub(r"\s*\([^)]*\)\s*$", "", first).strip()
    cleaned = first.strip(" .;:,")[:100]

    # Reject placeholders
    if cleaned in ("-", "—", "–", "?", "!", "n/a", "N/A", "none", "None"):
        return ""
    # Reject adjectival annotations
    if re.search(r"\badj\b", cleaned, flags=re.IGNORECASE):
        return ""

    return cleaned[:100]


def _normalize_gloss(gloss: str) -> str:
    """Normalize gloss field with enhanced HTML/template cleanup."""
    if not gloss:
        return ""

    # Remove HTML comments
    s = re.sub(r"<!--.*?-->", "", gloss, flags=re.DOTALL)

    # Remove <ref> tags (with or without content)
    s = re.sub(r"<ref[^>]*>.*?</ref>", "", s, flags=re.DOTALL)
    s = re.sub(r"<ref[^>]*/>", "", s)
    s = re.sub(r"<ref></ref>", "", s)

    # Basic cleanup
    s = _strip_wiki_markup(s)
    s = normalize_unicode(s)
    s = _normalize_ws(s)

    # Remove unclosed template markers (artifacts from incomplete parsing)
    s = re.sub(r"\}\}+$", "", s)  # Trailing }}
    s = re.sub(r"^\{\{+", "", s)  # Leading {{
    s = re.sub(r",\s*\}\}+", "", s)  # , followed by }}

    # Clean up leading/trailing punctuation artifacts
    s = re.sub(r"^[,:;\s]+|[,:;\s]+$", "", s)

    # Remove parenthetical notes (but keep if entire gloss is in parens)
    if not re.match(r"^\([^)]+\)$", s):
        s = re.sub(r"\s*\([^)]*\)\s*", " ", s)

    s = _normalize_ws(s)
    return s


def _extract_best_gloss(block: str) -> str:
    """Pick the first non-empty, non-inflectional gloss from a POS block."""
    if not block:
        return ""

    for mg in GLOSS_LINE_RE.finditer(block):
        raw = mg.group(1) or ""
        # Strip simple link brackets early
        raw = _normalize_ws(re.sub(r"\[\[|\]\]", "", raw))
        gloss_norm = _normalize_gloss(raw)
        if not gloss_norm:
            continue

        low = gloss_norm.lower()
        # Skip pure inflectional / meta definitions
        if low.startswith(("wikipedia", "wiktionary")):
            continue
        if low.startswith(
            (
                "plural of",
                "alternative form of",
                "alternative spelling of",
                "inflection of",
            )
        ):
            continue

        return gloss_norm

    return ""


# ============================================================================
# TITLE VALIDATION
# ============================================================================


def is_candidate_title(t: str) -> bool:
    """Reject proper nouns, acronyms, multiword, foreign chars."""
    if not t:
        return False
    if (
        " " in t
        or t.count("-") > 1
        or t.startswith("-")
        or t.endswith("-")
        or "'" in t
    ):
        return False
    if not t[0].islower():
        if not (t[0].isupper() and (len(t) == 1 or t[1:].islower())):
            return False
    low = t.lower()
    return all(ch in ROMANIAN_CHARS for ch in low)


def _has_romanian_section(wt: str) -> bool:
    return bool(ROMANIAN_SECTION_RE.search(wt))


def _is_uncountable(ro_section: str) -> bool:
    """Check if wikitext contains uncountable noun template."""
    return bool(UNCOUNTABLE_RE.search(ro_section))


def has_verb_pos(ro_section: str) -> bool:
    """Return True if the Romanian section has a verb POS heading."""
    for m in POS_HEAD_RE.finditer(ro_section):
        raw_label = m.group(1).lower()
        if POS_MAP.get(raw_label) == "V":
            return True
    return False


# ============================================================================
# TEMPLATE PARSING
# ============================================================================


def _extract_gender_from_template(block: str) -> Optional[str]:
    """Extract gender from {{ro-noun}} template."""
    # Pattern: {{ro-noun|m|... or {{ro-noun|f|... or {{ro-noun|n|...
    m = re.search(r"\{\{ro-noun\s*\|\s*([mfn])\s*(?:\||}})", block, re.I)
    if m:
        g_raw = m.group(1).lower()
        return GENDER_MAP.get(g_raw, "")
    return ""


def _extract_plural_from_template(block: str) -> str:
    """Extract plural from inline template parameters."""
    for pattern in PLURAL_KEYS:
        m = pattern.search(block)
        if m:
            return _normalize_ws(m.group(1))
    return ""


def _extract_from_templates(wikitext: str) -> dict[str, Any]:
    """Extract metadata from ro-noun, ro-adj, ro-verb templates.

    Returns dict with: plural, gender, m_pl, is_propn
    """
    result: dict[str, Any] = {
        "plural": None,
        "gender": None,
        "m_pl": None,
        "is_propn": False,
    }

    # Check for proper noun
    if "{{ro-proper noun" in wikitext.lower():
        result["is_propn"] = True

    # Extract from ro-noun template
    m_noun = re.search(
        r"\{\{ro-noun\s*\|\s*([mfn])\s*\|\s*([^\n\|\}]+)",
        wikitext,
        re.I,
    )
    if m_noun:
        result["gender"] = GENDER_MAP.get(m_noun.group(1).lower(), "")
        result["plural"] = _clean_plural(m_noun.group(2))

    # Extract masculine plural from ro-adj
    m_adj = re.search(
        r"\{\{ro-adj(?:\s*\|\s*([^\n\|\}]+))?(?:\s*\|\s*([^\n\|\}]+))?",
        wikitext,
        re.I,
    )
    if m_adj and m_adj.group(1):
        result["m_pl"] = _clean_plural(m_adj.group(1))

    return result


def _extract_etym_language_tag(wikitext: str) -> Optional[str]:
    """Extract etymology source language tag."""
    m = ETYM_LANG_RE.search(wikitext)
    if m:
        return m.group(1).lower()
    return None


def _strip_final_vowel(lemma: str) -> str:
    """Return lemma minus a final vowel, if present (>= 3 chars)."""
    s = normalize_unicode(lemma)
    if len(s) > 2 and s[-1] in RO_VOWELS:
        return s[:-1]
    return s


def _build_noun_stem_index(
    noun_lemmas: Set[str],
) -> DefaultDict[str, Set[str]]:
    """
    Map noun stems (lowercased, final vowel stripped) -> set of noun lemmas.

    Example:
        pădure -> stem 'pădur'
        drept  -> stem 'drept' (no final vowel to strip)
    """
    index: DefaultDict[str, Set[str]] = defaultdict(set)
    for lemma in noun_lemmas:
        lemma_norm = normalize_unicode(lemma)
        stem = _strip_final_vowel(lemma_norm).lower()
        if stem:
            index[stem].add(lemma_norm)
    return index


def _candidate_base_stems_for_verb(verb_lower: str, suffix: str) -> set[str]:
    """
    Given a lowercased verb and a denominal suffix (e.g. 'i', 'ui'),
    return possible *noun stems* (lowercased) that could underlie it.

    We consider:
      verb_lower = stem + suffix
      verb_lower = PREFIX + stem + suffix, where PREFIX in DERIV_PREFIXES.
    """
    stems: set[str] = set()
    # Remove the verbal suffix to get the "core"
    if not verb_lower.endswith(suffix):
        return stems

    core = verb_lower[: -len(suffix)]
    if len(core) >= 2:
        stems.add(core)  # bare denominal, no prefix

    # Prefixed denominals: PREFIX + noun_stem
    for pref in DERIV_PREFIXES_SORTED:
        if core.startswith(pref):
            base = core[len(pref) :]
            if len(base) >= 2:
                stems.add(base)

    return stems


def _augment_denominal_verbs_with_heuristics(
    noun_lemmas: Set[str],
    verb_lemmas: Set[str],
) -> int:
    """
    Add extra denominal verb links purely from transparent morphology.

    Heuristic:
      - Look for verbs ending in -i / -ui / í / uí (DENOMINAL_VERB_SUFFIXES).
      - Remove the suffix to get a core.
      - Treat core either as:
            core = noun_stem
        or  core = PREFIX + noun_stem, where PREFIX in DERIV_PREFIXES.
      - Match noun_stem against noun stems (noun minus final vowel).
      - Only add (N, V) if the noun lemma actually appears in the
        Romanian section of the verb's Wiktionary entry.

    This fills the global DENOMINAL_VERBS[base_noun] -> {verb, ...}
    and returns the number of new (noun, verb) pairs added.
    """
    noun_index = _build_noun_stem_index(noun_lemmas)
    added_pairs = 0

    for verb in verb_lemmas:
        verb_norm = normalize_unicode(verb)
        verb_lower = verb_norm.lower()

        # Check each allowed denominal suffix
        suffix_matched = None
        for suff in DENOMINAL_VERB_SUFFIXES:
            if verb_lower.endswith(suff):
                suffix_matched = suff
                break
        if suffix_matched is None:
            continue

        # Derive candidate noun stems (with and without prefixes)
        stem_candidates = _candidate_base_stems_for_verb(
            verb_lower, suffix_matched
        )
        if not stem_candidates:
            continue

        # Collect all noun lemmas whose stems match any candidate
        candidate_bases: set[str] = set()
        for stem in stem_candidates:
            bases = noun_index.get(stem)
            if bases:
                candidate_bases.update(bases)

        if not candidate_bases:
            continue

        # Optional precision boost: require that the noun lemma
        # actually occurs somewhere in the Romanian section
        wt = _get_wikitext_cached(verb_norm)
        m = ROMANIAN_SECTION_RE.search(wt)
        ro_section = m.group(1) if m else wt
        ro_lower = normalize_unicode(ro_section).lower()

        for base in candidate_bases:
            base_norm = normalize_unicode(base)
            base_lower = base_norm.lower()

            if base_lower not in ro_lower:
                # Drop this check if you prefer more recall over precision
                continue

            if verb_norm not in DENOMINAL_VERBS[base_norm]:
                DENOMINAL_VERBS[base_norm].add(verb_norm)
                added_pairs += 1

    return added_pairs


def register_verb_derivations(title: str, ro_section: str) -> None:
    """Record denominal/de-adjectival verb links from {{af|ro|...}} templates.

    This fills the global DENOMINAL_VERBS / DEADJECTIVAL_VERBS maps so that
    noun / adjective entries can later expose their derived verbs.
    """
    # Match the whole af-template argument list after `af|ro|`
    for m in re.finditer(
        r"\{\{\s*af\s*\|\s*ro\s*\|([^{}]+)\}\}", ro_section, flags=re.I
    ):
        parts = [p.strip() for p in m.group(1).split("|") if p.strip()]
        if not parts:
            continue

        base = normalize_unicode(parts[0])
        if not base or base == title:
            continue

        # Try to detect base POS (pos1=n/adj). Default to noun if unknown.
        base_pos = ""
        for p in parts[1:]:
            if p.lower().startswith("pos1="):
                base_pos = p.split("=", 1)[1].strip().lower()
                break

        base_pos_norm = base_pos.split(",")[0].strip()  # handle e.g. "adj, tr"
        if base_pos_norm.startswith("adj"):
            DEADJECTIVAL_VERBS[base].add(title)
        else:
            # Treat no pos1, or anything not clearly adjectival, as denominal
            DENOMINAL_VERBS[base].add(title)


def register_adj_derivations(title: str, ro_section: str) -> None:
    """Record denominal and de-adjectival adjective links from {{af|ro|...}}.

    Fills:
      - DENOMINAL_ADJS[base_noun] -> set(derived_adjectives)
      - DEADJECTIVAL_ADJS[base_adj] -> set(derived_adjectives)
    """
    for m in re.finditer(
        r"\{\{\s*af\s*\|\s*ro\s*\|([^{}]+)\}\}", ro_section, flags=re.I
    ):
        parts = [p.strip() for p in m.group(1).split("|") if p.strip()]
        if not parts:
            continue

        base = normalize_unicode(parts[0])
        if not base or base == title:
            continue

        base_pos = ""
        for p in parts[1:]:
            if p.lower().startswith("pos1="):
                base_pos = p.split("=", 1)[1].strip().lower()
                break

        base_pos_norm = base_pos.split(",")[0].strip()
        if base_pos_norm.startswith("adj"):
            DEADJECTIVAL_ADJS[base].add(title)
        else:
            # Assume noun base if not clearly adjectival
            DENOMINAL_ADJS[base].add(title)


# Regex for extracting base nouns from etymology text using {{m|ro|...}}
# Stops at |, }, <, or whitespace to avoid template parameter artifacts
_ETYM_M_RE = re.compile(
    r"\{\{\s*(?:m|mention|l|link)\s*\|\s*ro\s*\|\s*"
    r"([a-zăâîșțA-ZĂÂÎȘȚ][a-zăâîșțA-ZĂÂÎȘȚ\-]*?)(?=[|}<\s])",
    re.I,
)

# Derived terms section pattern
_DERIVED_TERMS_RE = re.compile(
    r"={3,4}\s*(?:Derived terms|Termeni derivați)\s*={3,4}"
    r"\s*(.*?)(?=\n={2,4}\s|\Z)",
    re.DOTALL | re.I,
)

# Suffix categories for derivational mining
DERIV_SUFFIX_CATEGORIES = [
    # Verbalizers (front-vowel → palatalization triggers)
    ("Category:Romanian terms suffixed with -i", "V"),
    ("Category:Romanian terms suffixed with -ui", "V"),
    ("Category:Romanian terms suffixed with -a", "V"),
    # Adjective-forming suffixes
    ("Category:Romanian terms suffixed with -esc", "ADJ"),
    ("Category:Romanian terms suffixed with -ic", "ADJ"),
    ("Category:Romanian terms suffixed with -ică", "ADJ"),
    ("Category:Romanian terms suffixed with -ist", "ADJ"),
]

# Limit for suffix category mining
SUFFIX_CAT_LIMIT = 5000


def _mine_derivation_categories() -> int:
    """Mine Wiktionary suffix categories to discover N→V and N→Adj links.

    For each suffix category, fetches members and then checks their
    etymology to find the base noun/adjective.  Populates the global
    DENOMINAL_VERBS / DENOMINAL_ADJS / DEADJECTIVAL_VERBS maps.

    Returns total number of new pairs added.
    """
    total_added = 0

    for cat_name, derived_pos in DERIV_SUFFIX_CATEGORIES:
        print(f"  Mining {cat_name}...")
        members = _fetch_category_members(
            WIKI_API_EN, cat_name, limit=SUFFIX_CAT_LIMIT
        )
        if not members:
            print("    No members found")
            continue

        # Batch prefetch uncached members for efficiency
        uncached = [
            t
            for t in members
            if t not in _WT_CACHE_EN and t not in _WT_CACHE_RO
        ]
        if uncached:
            print(f"    Batch-prefetching {len(uncached)} uncached...")
            for i in range(0, len(uncached), 800):
                chunk = uncached[i : i + 800]
                _batch_fetch_wikitext(WIKI_API_EN, chunk, _WT_CACHE_EN)

        added = 0
        for title in members:
            title_norm = normalize_unicode(title)
            wt = _get_wikitext_cached(title_norm)
            if not wt:
                continue
            m = ROMANIAN_SECTION_RE.search(wt)
            if not m:
                continue
            ro_section = m.group(1)

            # Extract base words from etymology via multiple patterns
            bases = set()

            # Pattern 1: {{af|ro|base|-suffix}} (already used by
            # register_verb/adj_derivations, but we re-check here
            # for items not in the verb/adj title sets)
            for af_m in re.finditer(
                r"\{\{\s*af\s*\|\s*ro\s*\|([^{}]+)\}\}",
                ro_section,
                flags=re.I,
            ):
                parts = [
                    p.strip() for p in af_m.group(1).split("|") if p.strip()
                ]
                if parts:
                    base = normalize_unicode(parts[0])
                    if base and base != title_norm:
                        bases.add(base)

            # Pattern 2: {{m|ro|base}} in etymology prose
            for m_m in _ETYM_M_RE.finditer(ro_section):
                base = normalize_unicode(m_m.group(1))
                if base and base != title_norm and len(base) >= 2:
                    bases.add(base)

            # Register discovered links
            for base in bases:
                if derived_pos == "V":
                    if title_norm not in DENOMINAL_VERBS.get(base, set()):
                        DENOMINAL_VERBS[base].add(title_norm)
                        added += 1
                elif derived_pos == "ADJ":
                    if title_norm not in DENOMINAL_ADJS.get(base, set()):
                        DENOMINAL_ADJS[base].add(title_norm)
                        added += 1

        print(f"    {len(members)} members, {added} new pairs")
        total_added += added

    return total_added


def extract_derived_terms_from_entry(title: str, ro_section: str) -> None:
    """Extract derivational links from the Derived terms section.

    Looks for ====Derived terms==== sections in the Romanian entry
    and parses linked terms, checking if they are verbs or adjectives.
    Populates DENOMINAL_VERBS and DENOMINAL_ADJS.
    """
    title_norm = normalize_unicode(title)

    m = _DERIVED_TERMS_RE.search(ro_section)
    if not m:
        return

    derived_block = m.group(1)

    # Extract linked terms: [[word]] or {{l|ro|word}}
    terms = set()
    for link_m in re.finditer(r"\[\[([^\]|]+)", derived_block):
        term = normalize_unicode(link_m.group(1).strip())
        if term and term != title_norm:
            terms.add(term)
    for l_m in re.finditer(
        r"\{\{\s*(?:l|m)\s*\|\s*ro\s*\|\s*([^|}]+)",
        derived_block,
        flags=re.I,
    ):
        term = normalize_unicode(l_m.group(1).strip())
        if term and term != title_norm:
            terms.add(term)

    # For each derived term, check its POS from cached wikitext only
    # (avoid triggering API calls for uncached terms)
    for term in terms:
        if term not in _WT_CACHE_EN and term not in _WT_CACHE_RO:
            continue
        wt = _get_wikitext_cached(term)
        if not wt:
            continue
        sec_m = ROMANIAN_SECTION_RE.search(wt)
        if not sec_m:
            continue
        ro = sec_m.group(1)

        is_verb = has_verb_pos(ro)
        is_adj = False
        for pos_m in POS_HEAD_RE.finditer(ro):
            if POS_MAP.get(pos_m.group(1).lower()) == "ADJ":
                is_adj = True
                break

        if is_verb:
            DENOMINAL_VERBS[title_norm].add(term)
        if is_adj:
            DENOMINAL_ADJS[title_norm].add(term)


# ============================================================================
# IPA EXTRACTION
# ============================================================================


def _extract_ipa_from_wikitext(ro_section: str) -> list[str]:
    results = []
    for m in IPA_WT_RE.finditer(ro_section):
        cand = _clean_ipa_raw(m.group(1))
        if cand:
            results.append(cand)
    return list(dict.fromkeys(results))


def _extract_ipa_list_from_html(html: str) -> list[str]:
    """Extract IPA transcriptions from rendered HTML."""
    if not HAS_BS4 or not html or BeautifulSoup is None:
        return []

    soup = BeautifulSoup(html, "html.parser")
    candidates = []

    # Find IPA spans
    for span in soup.select(IPA_SPAN_SELECTOR):
        text = span.get_text()
        if IPA_HINT.search(text):
            candidates.append(text)

    # Extract from brackets/slashes
    if not candidates:
        for m in IPA_TOK_RE.finditer(html):
            tok = m.group(0)
            if IPA_HINT.search(tok):
                candidates.append(tok)

    # Clean and deduplicate
    cleaned = []
    seen = set()
    for c in candidates:
        c_clean = _clean_ipa_raw(c)
        if c_clean and c_clean not in seen:
            cleaned.append(c_clean)
            seen.add(c_clean)

    return cleaned


def _fetch_html_section(api: str, title: str) -> str:
    """Fetch rendered HTML for entire page."""
    # Check HTML cache
    if api == WIKI_API_EN and title in _HTML_CACHE_EN:
        return _HTML_CACHE_EN[title]
    if api == WIKI_API_RO and title in _HTML_CACHE_RO:
        return _HTML_CACHE_RO[title]
    try:
        params = {
            "action": "parse",
            "page": title,
            "prop": "text",
            "format": "json",
        }
        data = get(api, params)
        _sleep_jitter(THROTTLE_DELAY)
        if "parse" in data and "text" in data["parse"]:
            # MediaWiki API responses are Any-typed; narrow to str
            # at the cache-write boundary.
            html = str(data["parse"]["text"].get("*", ""))
            if api == WIKI_API_EN:
                _HTML_CACHE_EN[title] = html
            else:
                _HTML_CACHE_RO[title] = html
            return html
    except (requests.RequestException, KeyError, ValueError):
        pass
    return ""


def _get_ipa_for_form(title: str) -> list[str]:
    """Get IPA for a word form from EN/RO HTML rendering."""
    if not ENABLE_IPA_FETCH:
        return []

    # Check cache
    if title in _IPA_CACHE:
        return _IPA_CACHE[title]

    ipas: list[str] = []

    # EN Wiktionary first
    html_en = _fetch_html_section(WIKI_API_EN, title)
    if html_en:
        ipas = _extract_ipa_list_from_html(html_en)

    # Fallback to RO Wiktionary
    if not ipas:
        html_ro = _fetch_html_section(WIKI_API_RO, title)
        if html_ro:
            ipas = _extract_ipa_list_from_html(html_ro)

    _IPA_CACHE[title] = ipas
    return ipas


# ============================================================================
# PLURAL EXTRACTION FROM HTML TABLES
# ============================================================================


def _sanitize_table_plural(text: str) -> Optional[str]:
    """Sanitize plural extracted from HTML table."""
    if not text:
        return None

    # Remove HTML tags
    text = re.sub(r"<[^>]+>", "", text)
    text = _strip_wiki_markup(text)
    text = _normalize_ws(text)

    # Reject if too short or looks like markup
    if len(text) < 3 or text in ("-", "—", "–", "i", "e", "uri"):
        return None

    return text


def _extract_plural_from_table(html: str) -> Optional[str]:
    """Extract plural from declension table in rendered HTML."""
    if not HAS_BS4 or not html or BeautifulSoup is None:
        return None

    soup = BeautifulSoup(html, "html.parser")

    # Look for table with class "inflection-table"
    for table in soup.find_all("table", class_=re.compile(r"inflection")):
        # Find cells with "plural" or "pl" label
        for cell in table.find_all(["td", "th"]):
            text = cell.get_text().lower()
            if "plural" in text or text.strip() in ("pl", "pl."):
                # Get next cell (might contain plural form)
                next_cell = cell.find_next("td")
                if next_cell:
                    candidate = _sanitize_table_plural(next_cell.get_text())
                    if candidate:
                        return candidate

    return None


def _confirm_plural_via_tables_or_templates(
    title: str,
    block: Optional[str] = None,
    tpl: Optional[dict[str, Any]] = None,
) -> str:
    """
    Return confirmed plural from template params, POS block, or HTML table.
    """

    def _is_valid_plural(s: str) -> bool:
        if not s:
            return False
        s_lower = s.lower()
        if s_lower in ("-", "—", "–", "", "e", "i", "uri"):
            return False
        if len(s) < 3:
            return False
        return True

    # Try template first
    if tpl and tpl.get("plural"):
        cand = _clean_plural(tpl["plural"])
        if _is_valid_plural(cand):
            return cand

    # Try inline template in block
    if block:
        inline_pl = _extract_plural_from_template(block)
        if inline_pl:
            cand = _clean_plural(inline_pl)
            if _is_valid_plural(cand):
                return cand

    # Check cache first to avoid repeated HTML parsing
    if title in _PLURAL_TABLE_CACHE:
        cached = _PLURAL_TABLE_CACHE[title]
        return cached if cached else ""

    # Try HTML tables (EN then RO)
    result = ""
    if HAS_BS4:
        html_en = _fetch_html_section(WIKI_API_EN, title)
        cand_en = _extract_plural_from_table(html_en)
        if cand_en:
            cand = _clean_plural(cand_en)
            if _is_valid_plural(cand):
                result = cand
                _PLURAL_TABLE_CACHE[title] = result
                return result

        html_ro = _fetch_html_section(WIKI_API_RO, title)
        cand_ro = _extract_plural_from_table(html_ro)
        if cand_ro:
            cand = _clean_plural(cand_ro)
            if _is_valid_plural(cand):
                result = cand
                _PLURAL_TABLE_CACHE[title] = result
                return result

    # Cache negative result to avoid repeated lookups
    _PLURAL_TABLE_CACHE[title] = None
    return ""


# ============================================================================
# WIKITEXT CACHE
# ============================================================================


def _migrate_json_to_sqlite(
    json_path: str, db: _SqliteCache, label: str
) -> None:
    """One-time migration: load old JSON cache into SQLite DB."""
    if not os.path.exists(json_path):
        return
    if len(db) > 0:
        # SQLite already has data; skip migration
        return
    print(f"  Migrating {label} JSON → SQLite...")
    try:
        with open(json_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        count = 0
        # Batch insert for speed
        conn = db._conn
        conn.execute("BEGIN")
        for key, value in data.items():
            conn.execute(
                f"INSERT OR IGNORE INTO {db._table} "
                f"(key, value) VALUES (?, ?)",
                (key, value),
            )
            count += 1
        conn.commit()
        print(f"  Migrated {count} {label} entries")
        # Rename old JSON so we don't re-migrate
        os.rename(json_path, json_path + ".migrated")
    except (json.JSONDecodeError, IOError, MemoryError) as e:
        print(f"  WARNING: migration failed for {label}: {e}")


def load_disk_cache() -> None:
    """Load wikitext caches from JSON, HTML caches from SQLite."""
    global _WT_CACHE_EN, _WT_CACHE_RO, _HTML_CACHE_EN, _HTML_CACHE_RO

    # Wikitext EN (JSON — small, ~55MB)
    if os.path.exists(CACHE_EN_PATH):
        try:
            with open(CACHE_EN_PATH, "r", encoding="utf-8") as f:
                _WT_CACHE_EN = json.load(f)
            print(f"Loaded {len(_WT_CACHE_EN)} EN wikitext cache entries")
        except (json.JSONDecodeError, IOError):
            _WT_CACHE_EN = {}

    # Wikitext RO (JSON — small, ~38MB)
    if os.path.exists(CACHE_RO_PATH):
        try:
            with open(CACHE_RO_PATH, "r", encoding="utf-8") as f:
                _WT_CACHE_RO = json.load(f)
            print(f"Loaded {len(_WT_CACHE_RO)} RO wikitext cache entries")
        except (json.JSONDecodeError, IOError):
            _WT_CACHE_RO = {}

    # HTML caches: SQLite (instant open, no bulk load)
    _HTML_CACHE_EN = _SqliteCache(HTML_DB_PATH, "html_en")
    _HTML_CACHE_RO = _SqliteCache(HTML_DB_PATH, "html_ro")

    # Migrate old JSON caches to SQLite if they exist
    old_en = str(_base_dir / "html_cache_en.json")
    old_ro = str(_base_dir / "html_cache_ro.json")
    _migrate_json_to_sqlite(old_en, _HTML_CACHE_EN, "EN HTML")
    _migrate_json_to_sqlite(old_ro, _HTML_CACHE_RO, "RO HTML")

    print(
        f"HTML cache (SQLite): "
        f"{len(_HTML_CACHE_EN)} EN, {len(_HTML_CACHE_RO)} RO"
    )


def save_disk_cache() -> None:
    """Save wikitext caches to JSON. HTML is in SQLite (auto-saved)."""
    # Wikitext EN
    try:
        with open(CACHE_EN_PATH, "w", encoding="utf-8") as f:
            json.dump(_WT_CACHE_EN, f, ensure_ascii=False)
        print(f"Saved {len(_WT_CACHE_EN)} EN wikitext cache entries")
    except IOError as e:
        print(f"Warning: failed to save EN wikitext cache: {e}")

    # Wikitext RO
    try:
        with open(CACHE_RO_PATH, "w", encoding="utf-8") as f:
            json.dump(_WT_CACHE_RO, f, ensure_ascii=False)
        print(f"Saved {len(_WT_CACHE_RO)} RO wikitext cache entries")
    except IOError as e:
        print(f"Warning: failed to save RO wikitext cache: {e}")

    # HTML caches are SQLite — already persisted on each insert
    print(
        f"HTML cache (SQLite): "
        f"{len(_HTML_CACHE_EN)} EN, {len(_HTML_CACHE_RO)} RO"
    )


def load_ipa_cache() -> None:
    """Load IPA cache from disk and clean entries."""
    global _IPA_CACHE

    if os.path.exists(IPA_CACHE_PATH):
        try:
            with open(IPA_CACHE_PATH, "r", encoding="utf-8") as f:
                raw_cache = json.load(f)
        except (json.JSONDecodeError, IOError):
            _IPA_CACHE = {}
            return

        cleaned: dict[str, list[str]] = {}
        for form, ipas in raw_cache.items():
            if not isinstance(ipas, list):
                continue
            new_list: list[str] = []
            for ipa in ipas:
                c = _clean_ipa_raw(ipa)
                if c:
                    new_list.append(c)
            if new_list:
                cleaned[form] = new_list

        _IPA_CACHE = cleaned
        print(f"Loaded {len(_IPA_CACHE)} IPA cache entries (cleaned)")
    else:
        _IPA_CACHE = {}


def _save_ipa_cache() -> None:
    """Save IPA cache to disk."""
    try:
        with open(IPA_CACHE_PATH, "w", encoding="utf-8") as f:
            json.dump(_IPA_CACHE, f, ensure_ascii=False)
        print(f"Saved {len(_IPA_CACHE)} IPA cache entries")
    except IOError as e:
        print(f"Warning: failed to save IPA cache: {e}")


def _batch_fetch_wikitext(
    api: str, titles: list[str], cache: dict[str, str]
) -> None:
    """Batch-fetch wikitext for multiple titles."""
    if not titles:
        return

    # Filter out already cached
    to_fetch = [t for t in titles if t not in cache]
    if not to_fetch:
        return

    # Fetch in batches of 50
    for i in range(0, len(to_fetch), 50):
        batch = to_fetch[i : i + 50]
        try:
            params = {
                "action": "query",
                "titles": "|".join(batch),
                "prop": "revisions",
                "rvprop": "content",
                "rvslots": "main",
                "redirects": "",
                "format": "json",
            }
            data = get(api, params)
            _sleep_jitter(THROTTLE_DELAY)

            if "query" in data and "pages" in data["query"]:
                for page in data["query"]["pages"].values():
                    title = page.get("title", "")
                    if "revisions" in page:
                        content = page["revisions"][0]["slots"]["main"]["*"]
                        cache[title] = content
        except (requests.RequestException, KeyError, ValueError) as exc:
            print(f"  WARNING: batch fetch failed: {exc}")
            continue


def _fetch_wikitext_from_api(api_url: str, title: str) -> Optional[str]:
    """Fetch wikitext from a single Wiktionary API."""
    try:
        params = {
            "action": "query",
            "titles": title,
            "prop": "revisions",
            "rvprop": "content",
            "rvslots": "main",
            "redirects": "",
            "format": "json",
        }
        data = get(api_url, params)
        _sleep_jitter(THROTTLE_DELAY)
        if "query" in data and "pages" in data["query"]:
            for page in data["query"]["pages"].values():
                if "revisions" in page:
                    # MediaWiki revisions slot is Any-typed; narrow to
                    # str at the trust boundary.
                    return str(page["revisions"][0]["slots"]["main"]["*"])
    except (requests.RequestException, KeyError, ValueError):
        pass
    return None


def _get_wikitext_cached(title: str) -> str:
    """Get wikitext for title with caching (EN preferred, RO fallback)."""
    if title in _WT_CACHE_EN:
        return _WT_CACHE_EN[title]
    if title in _WT_CACHE_RO:
        return _WT_CACHE_RO[title]
    content = _fetch_wikitext_from_api(WIKI_API_EN, title)
    if content:
        _WT_CACHE_EN[title] = content
        return content
    content = _fetch_wikitext_from_api(WIKI_API_RO, title)
    if content:
        _WT_CACHE_RO[title] = content
        return content
    return ""


# ============================================================================
# CATEGORY MEMBER FETCHING
# ============================================================================


def _fetch_category_members(
    api: str, category: str, limit: int = 5000
) -> list[str]:
    """Fetch all members of a Wiktionary category."""
    members: list[str] = []
    cmcontinue = None
    while len(members) < limit:
        params: dict[str, Any] = {
            "action": "query",
            "list": "categorymembers",
            "cmtitle": category,
            "cmlimit": min(500, limit - len(members)),
            "format": "json",
        }
        if cmcontinue:
            params["cmcontinue"] = cmcontinue
        try:
            data = get(api, params)
            _sleep_jitter(THROTTLE_DELAY)
            if "query" not in data or "categorymembers" not in data["query"]:
                break
            for item in data["query"]["categorymembers"]:
                members.append(item["title"])
            if "continue" in data and "cmcontinue" in data["continue"]:
                cmcontinue = data["continue"]["cmcontinue"]
            else:
                break
        except (requests.RequestException, KeyError, ValueError):
            break
    return members[:limit]


# ============================================================================
# ENTRY PARSING
# ============================================================================


def parse_romanian_entry(
    title: str, skip_ipa: bool = False
) -> Optional[dict[str, Any]]:
    """Parse a Romanian Wiktionary entry and extract raw fields.

    Returns dict with extracted fields or None if entry invalid.
    """
    if not title:
        return None
    title_normalized = normalize_unicode(title)
    wt = _get_wikitext_cached(title_normalized)
    m = ROMANIAN_SECTION_RE.search(wt)
    if not m:
        return None
    ro = m.group(1)
    # Exclude proper nouns (except Steriade examples)
    if PROPN_RE.search(ro) and title_normalized not in STERIADE_EXAMPLES:
        return None
    found: dict[str, Optional[str]] = {
        "N": None,
        "ADJ": None,
        "V": None,
    }
    for pos_m in POS_HEAD_RE.finditer(ro):
        raw_label = pos_m.group(1).lower()
        kind = POS_MAP.get(raw_label)
        if not kind:
            continue
        block = pos_m.group(2)
        if found.get(kind) is None:
            found[kind] = block
    # Prefer nouns > adjectives > verbs
    pos, block = None, None
    for k in ("N", "ADJ", "V"):
        if found[k]:
            pos, block = k, found[k]
            break
    if pos is None or block is None:
        return None
    unc = _is_uncountable(ro)
    tpl = _extract_from_templates(ro)
    if tpl.get("is_propn"):
        return None
    result = {
        "lemma": title_normalized,
        "pos": pos,
        "gender": "",
        "plural": "",
        "derived_verbs": "",
        "derived_adj": "",
        "gloss": "",
        "etym_lang": "",
        "source": "",
        "notes": "",
        "ipa_raw_lemma": "",
        "ipa_raw_pl": "",
    }
    if pos == "N":
        gender = tpl.get("gender") or _extract_gender_from_template(block) or ""
        result["gender"] = gender
        plural = ""
        if tpl.get("plural"):
            plural = tpl["plural"]
        else:
            for pattern in PLURAL_KEYS:
                m_pl = re.search(pattern, block)
                if m_pl:
                    plural = _normalize_ws(m_pl.group(1))
                    break
            if not plural:
                mtab = TABLE_ANY_PL_RE.search(ro)
                if mtab:
                    plural = _normalize_ws(mtab.group(1))
        if plural:
            plural = _clean_plural(plural)
        # Only fall back to HTML table extraction if wikitext
        # gave us nothing usable.  This avoids expensive HTTP
        # calls for entries where the template already provided
        # a valid plural.
        if not plural or plural in ("-", ""):
            if not unc:
                confirmed = _confirm_plural_via_tables_or_templates(
                    title=title_normalized, block=block, tpl=tpl
                )
                if confirmed:
                    plural = confirmed
                else:
                    plural = ""
            else:
                plural = ""
        # Reject bare suffix forms that aren't real words
        if plural and re.fullmatch(r"(?:i|ii|e|uri)", plural.strip(), re.I):
            plural = ""
        result["plural"] = plural
    gloss = _extract_best_gloss(block)
    if gloss:
        result["gloss"] = gloss
    for em in ETYM_RE.finditer(ro):
        etym_text = em.group(1)
        ety_lang = _extract_etym_language_tag(etym_text)
        if ety_lang:
            result["etym_lang"] = ety_lang
            break
    # Attach denominal / de-adjectival verbs discovered in the first pass.
    if pos == "N":
        # denominal verbs
        verbs = sorted(DENOMINAL_VERBS.get(title_normalized, ()))
        if verbs:
            result["derived_verbs"] = " | ".join(verbs)

        # denominal adjectives
        adjs = sorted(DENOMINAL_ADJS.get(title_normalized, ()))
        if adjs:
            result["derived_adj"] = " | ".join(adjs)

    elif pos == "ADJ":
        # de-adjectival verbs
        verbs = sorted(DEADJECTIVAL_VERBS.get(title_normalized, ()))
        if verbs:
            result["derived_verbs"] = " | ".join(verbs)

        # adjective → adjective derivations
        adjs = sorted(DEADJECTIVAL_ADJS.get(title_normalized, ()))
        if adjs:
            result["derived_adj"] = " | ".join(adjs)
    if title_normalized in _WT_CACHE_RO:
        lang_code = "ro"
    elif title_normalized in _WT_CACHE_EN:
        lang_code = "en"
    else:
        # Detect Romanian by diacritics
        low = title_normalized.lower()
        if any(ch in low for ch in "ăâîșț"):
            lang_code = "ro"
        else:
            lang_code = "en"
    result["source"] = _wiktionary_page_url(title_normalized, lang_code)
    notes = []
    if pos == "N" and unc:
        notes.append("uncountable")
    if pos == "N" and not result["plural"] and not unc:
        notes.append("needs plural confirmation")
    if pos == "N" and not result["gender"]:
        notes.append("needs gender confirmation")
    if notes:
        result["notes"] = " | ".join(notes)

    # IPA: prefer wikitext templates (no HTTP needed).
    # HTML-based IPA fetch only if HTML is already cached (from a
    # plural confirmation fetch).  Uncached HTML IPA fetches are
    # skipped — the downstream pipeline generates IPA via G2P.
    if not skip_ipa and pos in {"N", "ADJ"}:
        try:
            ipas_lemma = _extract_ipa_from_wikitext(ro)
            if not ipas_lemma and title_normalized in _HTML_CACHE_EN:
                ipas_lemma = _extract_ipa_list_from_html(
                    _HTML_CACHE_EN[title_normalized]
                )
            if not ipas_lemma and title_normalized in _HTML_CACHE_RO:
                ipas_lemma = _extract_ipa_list_from_html(
                    _HTML_CACHE_RO[title_normalized]
                )
            if ipas_lemma:
                result["ipa_raw_lemma"] = " | ".join(ipas_lemma)
        except (requests.RequestException, ValueError):
            pass
        if result["plural"]:
            try:
                ipas_pl = _get_ipa_for_form(result["plural"])
                if ipas_pl:
                    result["ipa_raw_pl"] = " | ".join(ipas_pl)
            except (requests.RequestException, ValueError):
                pass
    return result


# ============================================================================
# MAIN HARVEST
# ============================================================================


def harvest_data() -> pd.DataFrame:
    """Harvest raw data from Wiktionary."""
    all_rows = []
    seen = set()
    print("Discovering titles...")
    noun_titles_en = set(
        _fetch_category_members(
            WIKI_API_EN, "Category:Romanian nouns", limit=NOUN_LIMIT
        )
    )
    verb_titles_en = set(
        _fetch_category_members(
            WIKI_API_EN, "Category:Romanian verbs", limit=VERB_LIMIT
        )
    )
    adj_titles_en = set(
        _fetch_category_members(
            WIKI_API_EN, "Category:Romanian adjectives", limit=ADJ_LIMIT
        )
    )
    noun_titles_ro = set(
        _fetch_category_members(
            WIKI_API_RO, "Category:Substantive în română", limit=RO_NOUN_LIMIT
        )
    )
    adj_titles_ro = set(
        _fetch_category_members(
            WIKI_API_RO, "Category:Adjective în română", limit=RO_ADJ_LIMIT
        )
    )
    verb_titles_ro = set(
        _fetch_category_members(
            WIKI_API_RO, "Category:Verbe în română", limit=RO_VERB_LIMIT
        )
    )
    # Preserve order while deduplicating, include test cases
    all_titles = list(
        dict.fromkeys(
            list(noun_titles_en | verb_titles_en | adj_titles_en)
            + list(noun_titles_ro | verb_titles_ro | adj_titles_ro)
            + CANARY_LEMMAS
            + STERIADE_EXAMPLES
        )
    )
    print(f"Discovered {len(all_titles)} candidate titles")
    before = len(all_titles)
    all_titles = [t for t in all_titles if is_candidate_title(t)]
    print(f"Prefiltered: {len(all_titles)}/{before} kept")

    # Split into titles we need for verb-derivation scanning vs
    # titles we actually want rows for (nouns + adjectives).
    title_set = set(all_titles)
    verb_title_set = (verb_titles_en | verb_titles_ro) & title_set
    noun_adj_title_set = (
        (
            (noun_titles_en | noun_titles_ro | adj_titles_en | adj_titles_ro)
            & title_set
        )
        | set(CANARY_LEMMAS)
        | set(STERIADE_EXAMPLES)
    )

    titles_for_verbs = [t for t in all_titles if t in verb_title_set]
    titles_for_rows = [t for t in all_titles if t in noun_adj_title_set]
    adj_title_set = (adj_titles_en | adj_titles_ro) & title_set
    titles_for_adjs = [t for t in all_titles if t in adj_title_set]

    load_disk_cache()
    load_ipa_cache()

    # Batch prefetch to reduce API calls (both verbs and row entries)
    print("Batch-prefetching wikitext...")
    for i in range(0, len(all_titles), 800):
        chunk = all_titles[i : i + 800]
        _batch_fetch_wikitext(WIKI_API_EN, chunk, _WT_CACHE_EN)
        _batch_fetch_wikitext(WIKI_API_RO, chunk, _WT_CACHE_RO)
        if i % 2400 == 0:
            print(f"  Prefetched {i}/{len(all_titles)} titles...")

    # First pass: scan verb entries once to build denominal/de-adjectival maps.
    print("\nScanning verbs for denominal / de-adjectival derivations...")
    for i, title in enumerate(titles_for_verbs, 1):
        if i % 500 == 0:
            print(f"  Scanned {i}/{len(titles_for_verbs)} verb titles...")
        title_norm = normalize_unicode(title)
        wt = _get_wikitext_cached(title_norm)
        if not wt:
            continue
        m = ROMANIAN_SECTION_RE.search(wt)
        if not m:
            continue
        ro_section = m.group(1)
        # Skip if no verb POS heading (handles EN + RO labels via POS_MAP)
        if not has_verb_pos(ro_section):
            continue
        register_verb_derivations(title_norm, ro_section)

    print(
        f"  Found denominal verbs for {len(DENOMINAL_VERBS)} bases "
        f"and de-adjectival verbs for {len(DEADJECTIVAL_VERBS)} bases"
    )
    # Heuristic pass: string-based denominal verb discovery
    if ENABLE_DENOMINAL_HEURISTIC:
        print("\nAdding heuristic denominal verb links (-i / -ui)...")

        # Normalized noun and verb lemma sets, restricted to candidate titles
        noun_lemmas_for_heuristic: Set[str] = {
            normalize_unicode(t)
            for t in (noun_titles_en | noun_titles_ro) & title_set
        }
        verb_lemmas_for_heuristic: Set[str] = {
            normalize_unicode(t) for t in verb_title_set
        }

        heuristic_added = _augment_denominal_verbs_with_heuristics(
            noun_lemmas_for_heuristic,
            verb_lemmas_for_heuristic,
        )
        print(f"  Heuristic denominal verb pairs added: {heuristic_added}")

    # Second pass: scan adjectives for denominal/de-adjectival adj derivations.
    print("\nScanning adjectives for denominal / de-adjectival adjectives...")
    for i, title in enumerate(titles_for_adjs, 1):
        if i % 500 == 0:
            print(f"  Scanned {i}/{len(titles_for_adjs)} adjective titles...")
        title_norm = normalize_unicode(title)
        wt = _get_wikitext_cached(title_norm)
        if not wt:
            continue
        m = ROMANIAN_SECTION_RE.search(wt)
        if not m:
            continue
        ro_section = m.group(1)
        register_adj_derivations(title_norm, ro_section)

    print(
        f"  Found denominal adjectives for {len(DENOMINAL_ADJS)} bases "
        f"and de-adjectival adjectives for {len(DEADJECTIVAL_ADJS)} bases"
    )

    # Derivation category mining: discover N→V and N→Adj from
    # Wiktionary suffix categories (e.g., "terms suffixed with -i")
    print("\nMining Wiktionary suffix categories for derivations...")
    cat_added = _mine_derivation_categories()
    print(f"  Category mining added {cat_added} derivation pairs")
    print(
        f"  Total denominal verbs: "
        f"{sum(len(v) for v in DENOMINAL_VERBS.values())} pairs "
        f"across {len(DENOMINAL_VERBS)} bases"
    )
    print(
        f"  Total denominal adjs: "
        f"{sum(len(v) for v in DENOMINAL_ADJS.values())} pairs "
        f"across {len(DENOMINAL_ADJS)} bases"
    )

    # Derived terms pass: extract derivational links from noun entries
    print("\nExtracting derived terms from noun entries...")
    derived_count = 0
    for i, title in enumerate(titles_for_rows, 1):
        if i % 2000 == 0:
            print(f"  Scanned {i}/{len(titles_for_rows)} entries...")
        title_norm = normalize_unicode(title)
        wt = _get_wikitext_cached(title_norm)
        if not wt:
            continue
        sec_m = ROMANIAN_SECTION_RE.search(wt)
        if not sec_m:
            continue
        before = sum(len(v) for v in DENOMINAL_VERBS.values())
        before += sum(len(v) for v in DENOMINAL_ADJS.values())
        extract_derived_terms_from_entry(title_norm, sec_m.group(1))
        after = sum(len(v) for v in DENOMINAL_VERBS.values())
        after += sum(len(v) for v in DENOMINAL_ADJS.values())
        derived_count += after - before
    print(f"  Derived terms extraction added {derived_count} pairs")

    # Final pass: parse only nouns / adjectives into rows.
    print("\nParsing entries...")
    for i, title in enumerate(titles_for_rows, 1):
        if title in seen:
            continue
        if i % 100 == 0:
            print(f"  Processed {i}/{len(titles_for_rows)} entries...")
        entry = parse_romanian_entry(title)
        # We don't want bare verb lemmas as rows
        if entry and entry["pos"] != "V":
            all_rows.append(entry)
            seen.add(title)
    save_disk_cache()
    _save_ipa_cache()
    df = pd.DataFrame(all_rows)
    return df


def main() -> None:
    """Main entry point."""
    print("=" * 80)
    print("Romanian Wiktionary Harvester")
    print("=" * 80)
    df = harvest_data()
    print(f"\n{len(df)} entries harvested")

    # ===========================================================================
    # De-duplicate by (lemma, pos), preferring better sources
    # ===========================================================================
    print("\nDe-duplicating by (lemma, pos)...")

    # 1. Source priority: prefer Romanian Wiktionary over English
    df["source_rank"] = df["source"].apply(
        lambda s: 0 if "ro.wiktionary" in str(s) else 1
    )

    # 2. Prefer rows that actually have a gloss and a plural
    df["has_gloss"] = df["gloss"].notna() & (df["gloss"].str.strip() != "")
    df["has_plural"] = df["plural"].notna() & (df["plural"].str.strip() != "")

    # 3. Sort so the "best" row per lemma/pos comes first
    df = df.sort_values(
        by=["lemma", "pos", "source_rank", "has_gloss", "has_plural"],
        ascending=[True, True, True, False, False],
    )

    # 4. Drop duplicates, keeping the best for each (lemma, pos)
    duplicates_before = len(df)
    df = df.drop_duplicates(subset=["lemma", "pos"], keep="first")
    duplicates_removed = duplicates_before - len(df)

    # 5. Clean up helper columns
    df = df.drop(columns=["source_rank", "has_gloss", "has_plural"])

    print(f"Removed {duplicates_removed} duplicate (lemma, pos) pairs")
    print(f"{len(df)} unique entries remaining")

    # ===========================================================================

    print(f"\nWriting to {OUTPUT_CSV}...")
    df.to_csv(OUTPUT_CSV, index=False)
    print("\n" + "=" * 80)
    print("STATISTICS")
    print("=" * 80)
    print(f"Total entries: {len(df)}")
    print(f"Nouns: {len(df[df['pos'] == 'N'])}")
    print(f"Adjectives: {len(df[df['pos'] == 'ADJ'])}")
    total_denominal = sum(len(vs) for vs in DENOMINAL_VERBS.values())
    total_deadj = sum(len(vs) for vs in DEADJECTIVAL_VERBS.values())
    print(f"Denominal verbs (distinct pairs): {total_denominal}")
    print(f"De-adjectival verbs (distinct pairs): {total_deadj}")
    print(f"With plural: {len(df[df['plural'] != ''])}")
    print(f"With gender: {len(df[df['gender'] != ''])}")
    print(f"With IPA (lemma): {len(df[df['ipa_raw_lemma'] != ''])}")
    print(f"With IPA (plural): {len(df[df['ipa_raw_pl'] != ''])}")
    print("=" * 80)
    print("Done!")


if __name__ == "__main__":
    main()
