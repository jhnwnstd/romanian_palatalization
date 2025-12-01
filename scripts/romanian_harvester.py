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
import time
import unicodedata
from collections import defaultdict
from collections.abc import Mapping
from functools import lru_cache
from pathlib import Path
from typing import Any, DefaultDict, Optional, Set
from urllib.parse import quote

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

# Discovery limits
VERB_LIMIT = 15000
NOUN_LIMIT = 40000
ADJ_LIMIT = 25000
RO_NOUN_LIMIT = 10000
RO_ADJ_LIMIT = 8000
RO_VERB_LIMIT = 8000

# Rate limiting
THROTTLE_DELAY = 0.20

# Output CSV
OUTPUT_CSV = Path(__file__).parent.parent / "data" / "romanian_lexicon_raw.csv"

# Disk cache files (in data directory)
_base_dir = Path(__file__).parent.parent / "data"
_base_dir.mkdir(parents=True, exist_ok=True)
CACHE_EN_PATH = str(_base_dir / "wt_cache_en.json")
CACHE_RO_PATH = str(_base_dir / "wt_cache_ro.json")
IPA_CACHE_PATH = str(_base_dir / "ipa_cache.json")
HTML_EN_CACHE_PATH = str(_base_dir / "html_cache_en.json")
HTML_RO_CACHE_PATH = str(_base_dir / "html_cache_ro.json")

# Enable IPA fetch from HTML
ENABLE_IPA_FETCH = True

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

ROMANIAN_CHARS = set("abcdefghijklmnopqrstuvwxyzăâîșț")
ASCII_LETTERS = set("abcdefghijklmnopqrstuvwxyz")

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
    r"^===\s*(Noun|Substantiv|Adjective|Adjectiv|Verb)\s*===\s*$" r"(.*?)(?=^===|\Z)",
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

FEM_KEYS = [
    re.compile(r"\|\s*f\s*=\s*([^\n\|\}]+)", re.I),
    re.compile(r"\|\s*fem\s*=\s*([^\n\|\}]+)", re.I),
    re.compile(r"\|\s*feminine\s*=\s*([^\n\|\}]+)", re.I),
]

TABLE_ANY_PL_RE = re.compile(
    r"\{\{\s*ro-noun-table[^}]*\|\s*(?:pl|plural)\s*=\s*([^\n\|}\]]+)",
    re.I,
)

# Etymology language tag extraction
ETYM_LANG_RE = re.compile(
    r"\{\{\s*(?:bor|der|inh|lbor|calque)\s*\|\s*ro\s*\|\s*" r"([a-z]{2,3})\s*(?:\||}})",
    re.I,
)

# Derivational affix templates
AFFIX_GENERIC_RE = re.compile(
    r"\{\{\s*af\s*\|\s*ro\s*\|\s*([a-zăâîșțA-ZĂÂÎȘȚ\-]+)",
    re.I,
)
AFFIX_ADJ_RE = re.compile(
    r"\{\{\s*af\s*\|\s*ro\s*\|\s*([a-zăâîșțA-ZĂÂÎȘȚ\-]+)\s*\|\s*"
    r"\-os\s*\|\s*pos1\s*=\s*n",
    re.I,
)

UNCOUNTABLE_RE = re.compile(r"\{\{\s*(?:unc|uncountable)\s*(?:\|[^}]*)?\}\}", re.I)

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

_HTML_CACHE_EN: dict[str, str] = {}
_HTML_CACHE_RO: dict[str, str] = {}

# Cache for extracted plural forms from HTML tables (avoids repeated BeautifulSoup parsing)
_PLURAL_TABLE_CACHE: dict[str, Optional[str]] = {}

DENOMINAL_VERBS: DefaultDict[str, Set[str]] = defaultdict(set)
DEADJECTIVAL_VERBS: DefaultDict[str, Set[str]] = defaultdict(set)
DENOMINAL_ADJS: DefaultDict[str, Set[str]] = defaultdict(set)  # noun -> adjectives
DEADJECTIVAL_ADJS: DefaultDict[str, Set[str]] = defaultdict(set)  # adj  -> adjectives


@retry(
    reraise=True,
    stop=stop_after_attempt(6),
    wait=wait_exponential(multiplier=0.75, min=0.5, max=12),
    retry=retry_if_exception_type(requests.RequestException),
)
def get(url: str, params: Mapping[str, Any]) -> dict[str, Any]:
    """Fetch JSON from Wiktionary API with rate limiting and retry."""
    params = {"maxlag": "5", **params}
    r = _SESSION.get(url, params=params, timeout=30)
    if r.status_code in (429, 503):
        try:
            retry_after = float(r.headers.get("Retry-After", "2"))
        except ValueError:
            retry_after = 2.0
        time.sleep(min(10.0, retry_after + 1.0))
        raise requests.RequestException(f"Rate-limited: {r.status_code}")
    r.raise_for_status()
    data = r.json()
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


def normalize_ws(s: str) -> str:
    """Collapse whitespace to single spaces and strip edges."""
    return re.sub(r"\s+", " ", s or "").strip()


def strip_wiki_markup(s: str) -> str:
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


def clean_ipa_raw(ipa_str: str) -> str:
    """Remove IPA delimiters, HTML, percent-encoding; keep stress markers."""
    if not isinstance(ipa_str, str):
        return ""
    s = ipa_str.strip()
    if not s or s.lower() == "nan":
        return ""
    # HTML decoding and cleanup
    try:
        from urllib.parse import unquote

        s = unquote(s)
    except (ImportError, ValueError, TypeError):
        pass
    # Remove HTML tags
    s = re.sub(r"<[^>]+>", "", s)
    s = re.sub(r'\b(?:title|href|class|id|lang|span)="[^"]*"', "", s)
    s = re.sub(r"&[a-zA-Z]+;", " ", s)
    s = s.replace(">", "").replace("<", "").replace('"', "")
    s = s.replace("[", "").replace("]", "").replace("/", "")
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

    return s.strip()


def clean_plural(s: str) -> str:
    """Strip markup, take first variant, reject placeholders."""
    if not s:
        return ""
    s = strip_wiki_markup(s)
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


def normalize_gloss(gloss: str) -> tuple[str, list[str], list[str]]:
    """Normalize gloss field with basic cleanup.

    Returns: (normalized_gloss, empty_list, empty_list)
    """
    if not gloss:
        return ("", [], [])

    # Basic cleanup
    s = strip_wiki_markup(gloss)
    s = normalize_unicode(s)
    s = normalize_ws(s)

    # Remove parenthetical notes
    s = re.sub(r"\s*\([^)]*\)\s*", " ", s)
    s = normalize_ws(s)

    return (s, [], [])


def extract_best_gloss(block: str) -> str:
    """Pick the first non-empty, non-inflectional gloss from a POS block."""
    if not block:
        return ""

    for mg in GLOSS_LINE_RE.finditer(block):
        raw = mg.group(1) or ""
        # Strip simple link brackets early
        raw = normalize_ws(re.sub(r"\[\[|\]\]", "", raw))
        gloss_norm, _, _ = normalize_gloss(raw)
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
    if " " in t or t.count("-") > 1 or t.startswith("-") or t.endswith("-") or "'" in t:
        return False
    if not t[0].islower():
        if not (t[0].isupper() and (len(t) == 1 or t[1:].islower())):
            return False
    low = t.lower()
    return all(ch in ROMANIAN_CHARS or ch in ASCII_LETTERS for ch in low)


def has_romanian_section(wt: str) -> bool:
    return bool(ROMANIAN_SECTION_RE.search(wt))


def is_uncountable(ro_section: str) -> bool:
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


def extract_gender_from_template(block: str) -> Optional[str]:
    """Extract gender from {{ro-noun}} template."""
    # Pattern: {{ro-noun|m|... or {{ro-noun|f|... or {{ro-noun|n|...
    m = re.search(r"\{\{ro-noun\s*\|\s*([mfn])\s*(?:\||}})", block, re.I)
    if m:
        g_raw = m.group(1).lower()
        return GENDER_MAP.get(g_raw, "")
    return ""


def extract_plural_from_template(block: str) -> str:
    """Extract plural from inline template parameters."""
    for pattern in PLURAL_KEYS:
        m = pattern.search(block)
        if m:
            return normalize_ws(m.group(1))
    return ""


def extract_from_templates(wikitext: str) -> dict[str, Any]:
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
        result["plural"] = clean_plural(m_noun.group(2))

    # Extract masculine plural from ro-adj
    m_adj = re.search(
        r"\{\{ro-adj(?:\s*\|\s*([^\n\|\}]+))?(?:\s*\|\s*([^\n\|\}]+))?",
        wikitext,
        re.I,
    )
    if m_adj and m_adj.group(1):
        result["m_pl"] = clean_plural(m_adj.group(1))

    return result


def extract_etym_language_tag(wikitext: str) -> Optional[str]:
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


def augment_denominal_verbs_with_heuristics(
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

        # Optional: Restrict heuristic to verbs starting with a
        # derivational prefix. To activate, uncomment the `continue`
        # and delete `pass` below.
        if not any(verb_lower.startswith(p) for p in DERIV_PREFIXES):
            # continue
            pass

        # Check each allowed denominal suffix
        suffix_matched = None
        for suff in DENOMINAL_VERB_SUFFIXES:
            if verb_lower.endswith(suff):
                suffix_matched = suff
                break
        if suffix_matched is None:
            continue

        # Derive candidate noun stems (with and without prefixes)
        stem_candidates = _candidate_base_stems_for_verb(verb_lower, suffix_matched)
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
        wt = get_wikitext_cached(verb_norm)
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


# ============================================================================
# IPA EXTRACTION
# ============================================================================


def extract_ipa_from_wikitext(ro_section: str) -> list[str]:
    results = []
    for m in IPA_WT_RE.finditer(ro_section):
        cand = clean_ipa_raw(m.group(1))
        if cand:
            results.append(cand)
    return list(dict.fromkeys(results))


def extract_ipa_list_from_html(html: str) -> list[str]:
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
        c_clean = clean_ipa_raw(c)
        if c_clean and c_clean not in seen:
            cleaned.append(c_clean)
            seen.add(c_clean)

    return cleaned


def fetch_html_section(api: str, title: str) -> str:
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
            html = data["parse"]["text"].get("*", "")
            if api == WIKI_API_EN:
                _HTML_CACHE_EN[title] = html
            else:
                _HTML_CACHE_RO[title] = html
            return html
    except (requests.RequestException, KeyError, ValueError):
        pass
    return ""


def get_ipa_for_form(title: str) -> list[str]:
    """Get IPA for a word form from EN/RO HTML rendering."""
    if not ENABLE_IPA_FETCH:
        return []

    # Check cache
    if title in _IPA_CACHE:
        return _IPA_CACHE[title]

    ipas: list[str] = []

    # EN Wiktionary first
    html_en = fetch_html_section(WIKI_API_EN, title)
    if html_en:
        ipas = extract_ipa_list_from_html(html_en)

    # Fallback to RO Wiktionary
    if not ipas:
        html_ro = fetch_html_section(WIKI_API_RO, title)
        if html_ro:
            ipas = extract_ipa_list_from_html(html_ro)

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
    text = strip_wiki_markup(text)
    text = normalize_ws(text)

    # Reject if too short or looks like markup
    if len(text) < 3 or text in ("-", "—", "–", "i", "e", "uri"):
        return None

    return text


def extract_plural_from_table(html: str) -> Optional[str]:
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


def confirm_plural_via_tables_or_templates(
    title: str,
    block: Optional[str] = None,
    tpl: Optional[dict] = None,
) -> str:
    """
    Return confirmed plural from template params, POS block, or HTML table.
    """

    def is_valid_plural(s: str) -> bool:
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
        cand = clean_plural(tpl["plural"])
        if is_valid_plural(cand):
            return cand

    # Try inline template in block
    if block:
        inline_pl = extract_plural_from_template(block)
        if inline_pl:
            cand = clean_plural(inline_pl)
            if is_valid_plural(cand):
                return cand

    # Check cache first to avoid repeated HTML parsing
    if title in _PLURAL_TABLE_CACHE:
        cached = _PLURAL_TABLE_CACHE[title]
        return cached if cached else ""

    # Try HTML tables (EN then RO)
    result = ""
    if HAS_BS4:
        html_en = fetch_html_section(WIKI_API_EN, title)
        cand_en = extract_plural_from_table(html_en)
        if cand_en:
            cand = clean_plural(cand_en)
            if is_valid_plural(cand):
                result = cand
                _PLURAL_TABLE_CACHE[title] = result
                return result

        html_ro = fetch_html_section(WIKI_API_RO, title)
        cand_ro = extract_plural_from_table(html_ro)
        if cand_ro:
            cand = clean_plural(cand_ro)
            if is_valid_plural(cand):
                result = cand
                _PLURAL_TABLE_CACHE[title] = result
                return result

    # Cache negative result to avoid repeated lookups
    _PLURAL_TABLE_CACHE[title] = None
    return ""


# ============================================================================
# WIKITEXT CACHE
# ============================================================================


def load_disk_cache():
    """Load wikitext and HTML caches from disk."""
    global _WT_CACHE_EN, _WT_CACHE_RO, _HTML_CACHE_EN, _HTML_CACHE_RO

    # Wikitext EN
    if os.path.exists(CACHE_EN_PATH):
        try:
            with open(CACHE_EN_PATH, "r", encoding="utf-8") as f:
                _WT_CACHE_EN = json.load(f)
            print(f"Loaded {len(_WT_CACHE_EN)} EN wikitext cache entries")
        except (json.JSONDecodeError, IOError):
            _WT_CACHE_EN = {}

    # Wikitext RO
    if os.path.exists(CACHE_RO_PATH):
        try:
            with open(CACHE_RO_PATH, "r", encoding="utf-8") as f:
                _WT_CACHE_RO = json.load(f)
            print(f"Loaded {len(_WT_CACHE_RO)} RO wikitext cache entries")
        except (json.JSONDecodeError, IOError):
            _WT_CACHE_RO = {}

    # HTML EN
    if os.path.exists(HTML_EN_CACHE_PATH):
        try:
            with open(HTML_EN_CACHE_PATH, "r", encoding="utf-8") as f:
                _HTML_CACHE_EN = json.load(f)
            print(f"Loaded {len(_HTML_CACHE_EN)} EN HTML cache entries")
        except (json.JSONDecodeError, IOError):
            _HTML_CACHE_EN = {}

    # HTML RO
    if os.path.exists(HTML_RO_CACHE_PATH):
        try:
            with open(HTML_RO_CACHE_PATH, "r", encoding="utf-8") as f:
                _HTML_CACHE_RO = json.load(f)
            print(f"Loaded {len(_HTML_CACHE_RO)} RO HTML cache entries")
        except (json.JSONDecodeError, IOError):
            _HTML_CACHE_RO = {}


def save_disk_cache():
    """Save wikitext and HTML caches to disk."""
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

    # HTML EN
    try:
        with open(HTML_EN_CACHE_PATH, "w", encoding="utf-8") as f:
            json.dump(_HTML_CACHE_EN, f, ensure_ascii=False)
        print(f"Saved {len(_HTML_CACHE_EN)} EN HTML cache entries")
    except IOError as e:
        print(f"Warning: failed to save EN HTML cache: {e}")

    # HTML RO
    try:
        with open(HTML_RO_CACHE_PATH, "w", encoding="utf-8") as f:
            json.dump(_HTML_CACHE_RO, f, ensure_ascii=False)
        print(f"Saved {len(_HTML_CACHE_RO)} RO HTML cache entries")
    except IOError as e:
        print(f"Warning: failed to save RO HTML cache: {e}")


def load_ipa_cache():
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
                c = clean_ipa_raw(ipa)
                if c:
                    new_list.append(c)
            if new_list:
                cleaned[form] = new_list

        _IPA_CACHE = cleaned
        print(f"Loaded {len(_IPA_CACHE)} IPA cache entries (cleaned)")
    else:
        _IPA_CACHE = {}


def save_ipa_cache():
    """Save IPA cache to disk."""
    try:
        with open(IPA_CACHE_PATH, "w", encoding="utf-8") as f:
            json.dump(_IPA_CACHE, f, ensure_ascii=False)
        print(f"Saved {len(_IPA_CACHE)} IPA cache entries")
    except IOError as e:
        print(f"Warning: failed to save IPA cache: {e}")


def batch_fetch_wikitext(api: str, titles: list[str], cache: dict[str, str]) -> None:
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
        except (requests.RequestException, KeyError, ValueError):
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
            "format": "json",
        }
        data = get(api_url, params)
        _sleep_jitter(THROTTLE_DELAY)
        if "query" in data and "pages" in data["query"]:
            for page in data["query"]["pages"].values():
                if "revisions" in page:
                    return page["revisions"][0]["slots"]["main"]["*"]
    except (requests.RequestException, KeyError, ValueError):
        pass
    return None


def get_wikitext_cached(title: str) -> str:
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


def fetch_category_members(api: str, category: str, limit: int = 5000) -> list[str]:
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


def parse_romanian_entry(title: str, skip_ipa: bool = False) -> Optional[dict]:
    """Parse a Romanian Wiktionary entry and extract raw fields.

    Returns dict with extracted fields or None if entry invalid.
    """
    if not title:
        return None
    title_normalized = normalize_unicode(title)
    wt = get_wikitext_cached(title_normalized)
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
    unc = is_uncountable(ro)
    tpl = extract_from_templates(ro)
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
        gender = tpl.get("gender") or extract_gender_from_template(block) or ""
        result["gender"] = gender
        plural = ""
        if tpl.get("plural"):
            plural = tpl["plural"]
        else:
            for pattern in PLURAL_KEYS:
                m_pl = re.search(pattern, block)
                if m_pl:
                    plural = normalize_ws(m_pl.group(1))
                    break
            if not plural:
                mtab = TABLE_ANY_PL_RE.search(ro)
                if mtab:
                    plural = normalize_ws(mtab.group(1))
        if plural:
            plural = clean_plural(plural)
        # Confirm via tables/templates if plural looks incomplete
        if not plural or plural in ("-", "", "i", "e", "uri") or len(plural) < 3:
            if not unc:
                confirmed = confirm_plural_via_tables_or_templates(
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
    gloss = extract_best_gloss(block)
    if gloss:
        result["gloss"] = gloss
    for em in ETYM_RE.finditer(ro):
        etym_text = em.group(1)
        ety_lang = extract_etym_language_tag(etym_text)
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

    # We only need IPA for entries we actually keep (nouns / adjectives).
    if not skip_ipa and pos in {"N", "ADJ"}:
        try:
            # First try IPA from Romanian wikitext
            ipas_lemma = extract_ipa_from_wikitext(ro)
            # Fallback to HTML-based extraction if none found
            if not ipas_lemma:
                ipas_lemma = get_ipa_for_form(title_normalized)
            if ipas_lemma:
                result["ipa_raw_lemma"] = " | ".join(ipas_lemma)
        except (requests.RequestException, ValueError):
            pass
        if result["plural"]:
            try:
                ipas_pl = get_ipa_for_form(result["plural"])
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
        fetch_category_members(WIKI_API_EN, "Category:Romanian nouns", limit=NOUN_LIMIT)
    )
    verb_titles_en = set(
        fetch_category_members(WIKI_API_EN, "Category:Romanian verbs", limit=VERB_LIMIT)
    )
    adj_titles_en = set(
        fetch_category_members(
            WIKI_API_EN, "Category:Romanian adjectives", limit=ADJ_LIMIT
        )
    )
    noun_titles_ro = set(
        fetch_category_members(
            WIKI_API_RO, "Category:Substantive în română", limit=RO_NOUN_LIMIT
        )
    )
    adj_titles_ro = set(
        fetch_category_members(
            WIKI_API_RO, "Category:Adjective în română", limit=RO_ADJ_LIMIT
        )
    )
    verb_titles_ro = set(
        fetch_category_members(
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
        ((noun_titles_en | noun_titles_ro | adj_titles_en | adj_titles_ro) & title_set)
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
        batch_fetch_wikitext(WIKI_API_EN, chunk, _WT_CACHE_EN)
        batch_fetch_wikitext(WIKI_API_RO, chunk, _WT_CACHE_RO)
        if i % 2400 == 0:
            print(f"  Prefetched {i}/{len(all_titles)} titles...")

    # First pass: scan verb entries once to build denominal/de-adjectival maps.
    print("\nScanning verbs for denominal / de-adjectival derivations...")
    for i, title in enumerate(titles_for_verbs, 1):
        if i % 500 == 0:
            print(f"  Scanned {i}/{len(titles_for_verbs)} verb titles...")
        title_norm = normalize_unicode(title)
        wt = get_wikitext_cached(title_norm)
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
            normalize_unicode(t) for t in (noun_titles_en | noun_titles_ro) & title_set
        }
        verb_lemmas_for_heuristic: Set[str] = {
            normalize_unicode(t) for t in verb_title_set
        }

        heuristic_added = augment_denominal_verbs_with_heuristics(
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
        wt = get_wikitext_cached(title_norm)
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

    # Third pass: parse only nouns / adjectives into rows.
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
    save_ipa_cache()
    df = pd.DataFrame(all_rows)
    return df


def main() -> None:
    """Main entry point."""
    print("=" * 80)
    print("Romanian Wiktionary Harvester")
    print("=" * 80)
    df = harvest_data()
    print(f"\n{len(df)} entries harvested")
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
