"""DEX etymology audit for derived-verb validation.

Trust-boundary module. Given a derived-verb candidate's DEX HTML, decide
whether to accept the noun→verb pair into the inflection-dependence
contingency table.

Why this exists
---------------
Wiktionary's *Derived terms* section and our DEX-supplement detector both
produce candidate (noun, derived_verb) pairs that look right by surface
shape. Many of them are not actually denominal verbs of the noun in
question: DEX may host a homographic verb whose etymology cites a
different word (e.g. *dogi* "to crack" cites *doagă* "stave", not *dog*
"mastiff"), and some candidates are simply pages that happened to URL-
collide with our query string.

The audit classifies each candidate into one of six etymology categories
and accepts only those that name our noun (A, B) or are foreign
borrowings whose Romanian form is the synchronic noun-verb pair (D).
This is invoked once at the row's trust boundary; downstream code can
then assume every (noun, verb) pair is a real denominal pair.

The validation is conservative on purpose. We over-reject ambiguous
necunoscută cases rather than under-reject false positives that would
inflate Steriade-direction effects in the contingency table.
"""

from __future__ import annotations

import re
from enum import StrEnum
from typing import Final

# ---------------------------------------------------------------------------
# Public constants
# ---------------------------------------------------------------------------

# Back-formations whose DEX etymology is "necunoscută" but which behave
# synchronically as denominal pairs (the noun and verb share a root in
# modern Romanian, regardless of which came first historically). Kept by
# explicit allowlist because the etymology field doesn't license them.
USER_KEEP: Final[frozenset[str]] = frozenset(
    {"clint", "dinte", "vargă", "verigă"}
)


class EtymCategory(StrEnum):
    """How a DEX etymology field relates a derived verb to a base noun.

    Wire-stable string values. A, B, and D are accepted into the
    contingency table; C, E, and F are rejected.
    """

    # Strong denominal links — accepted
    CITES_OUR_NOUN = "A"
    CITES_PREFIXED_FORM = "B"
    FOREIGN_SOURCE = "D"

    # Reject categories
    CITES_DIFFERENT_WORD = "C"
    UNKNOWN_ORIGIN = "E"
    NO_ETYMOLOGY = "F"


ACCEPT_CATEGORIES: Final[frozenset[EtymCategory]] = frozenset(
    {
        EtymCategory.CITES_OUR_NOUN,
        EtymCategory.CITES_PREFIXED_FORM,
        EtymCategory.FOREIGN_SOURCE,
    }
)


# ---------------------------------------------------------------------------
# DEX HTML structure regexes
# ---------------------------------------------------------------------------

_VERB_POS_RE: Final[re.Pattern[str]] = re.compile(
    r'<span[^>]*class="[^"]*tree-pos-info[^"]*"[^>]*>verb',
    re.IGNORECASE,
)
_TONIC_RE: Final[re.Pattern[str]] = re.compile(
    r'<span[^>]*class="[^"]*tonic-accent[^"]*"[^>]*>(.*?)</span>',
    re.IGNORECASE | re.DOTALL,
)
_INFLECTED_RE: Final[re.Pattern[str]] = re.compile(
    r'<span[^>]*class="[^"]*tree-inflected-form[^"]*"[^>]*>.*?</span>',
    re.IGNORECASE | re.DOTALL,
)
_POSINFO_RE: Final[re.Pattern[str]] = re.compile(
    r'<span[^>]*class="[^"]*tree-pos-info[^"]*"[^>]*>.*?</span>',
    re.IGNORECASE | re.DOTALL,
)
_ANY_TAG: Final[re.Pattern[str]] = re.compile(r"<[^>]+>")
_TREE_HEADING_RE: Final[re.Pattern[str]] = re.compile(
    r'<h3[^>]*class="[^"]*tree-heading[^"]*"[^>]*>(.*?)</h3>',
    re.DOTALL | re.IGNORECASE,
)
_ETYM_RE: Final[re.Pattern[str]] = re.compile(
    r"etimologie\s*[:\s]\s*(.{1,200}?)"
    r"(?=\s+(?:DEX|MDA|DLRLC|DLRM|DEXI|CADE|NODEX|info|sinonime|"
    r"antonime|format|$))",
    re.IGNORECASE | re.DOTALL,
)
_FOREIGN_RE: Final[re.Pattern[str]] = re.compile(
    r"\b(?:limba\s+)?"
    r"(?:franceză|latină|slav[ăa]|veche|neogreacă|greacă|germană|"
    r"italiană|maghiară|engleză|turcă|spaniolă|rusă|bulgară|sârbă|"
    r"polonă|ucrainean|săsească|albaneză|portugheză|cehă|olandeză|"
    r"catalană|lat\.|fr\.|gr\.|sl\.|magh\.|germ\.|it\.|engl\.|tc\.|"
    r"sb\.|bg\.|rus\.|alb\.)",
    re.IGNORECASE,
)
_UNKNOWN_RE: Final[re.Pattern[str]] = re.compile(
    r"\bnecunoscut|et\.\s*nec\.|etimol\.\s*nec\.",
    re.IGNORECASE,
)
_PARTICLES: Final[frozenset[str]] = frozenset(
    {"cf", "vezi", "din", "de", "la", "se", "și", "sau"}
)


# ---------------------------------------------------------------------------
# Stem-variant generation
# ---------------------------------------------------------------------------


def stem_variants(noun: str) -> frozenset[str]:
    """Return plausible orthographic stem variants of a Romanian noun.

    Generates the bare noun, the gender-stripped stem, and forms
    differing by the standard Romanian vowel alternations (a↔ă, o↔oa,
    e↔ea). These are the variants we accept as "same root" when
    comparing an etymology citation to our base noun.
    """
    base = noun.lower()
    if len(base) > 1 and base[-1] in "ăaeio":
        stems = {base, base[:-1]}
    else:
        stems = {base}
    out: set[str] = set(stems)
    for s in stems:
        if "a" in s:
            out.add(s.replace("a", "ă"))
        if "ă" in s:
            out.add(s.replace("ă", "a"))
        if "o" in s:
            i = s.find("o")
            out.add(s[:i] + "oa" + s[i + 1 :])
        if "oa" in s:
            out.add(s.replace("oa", "o"))
        if "e" in s:
            i = s.find("e")
            out.add(s[:i] + "ea" + s[i + 1 :])
        if "ea" in s:
            out.add(s.replace("ea", "e"))
    return frozenset(out)


# ---------------------------------------------------------------------------
# Etymology extraction
# ---------------------------------------------------------------------------


def _heading_lemma(heading_html: str) -> str | None:
    """Extract the lemma word from a DEX tree-heading block."""
    h = _TONIC_RE.sub(r"\1", heading_html)
    h = _INFLECTED_RE.sub("", h)
    h = _POSINFO_RE.sub("", h)
    h = _ANY_TAG.sub("", h)
    h = re.sub(r"\s+", " ", h).strip()
    m = re.match(r"([^,;\s]+)", h)
    return m.group(1).strip().lower() if m else None


def get_verb_etymology(html: str, candidate: str) -> str | None:
    """Return the etymology field of the DEX verb tree whose lemma
    equals ``candidate``.

    Returns ``None`` if the page has no verb tree with the matching
    canonical lemma (a strong indicator that the candidate is not the
    page's primary lemma — e.g. *fabrici* is actually a noun plural form
    whose canonical entry is *fabrică*). Returns ``""`` if the verb tree
    exists but has no parseable etymology field.
    """
    h3 = list(_TREE_HEADING_RE.finditer(html))
    for i, m in enumerate(h3):
        if not _VERB_POS_RE.search(m.group(1)):
            continue
        if _heading_lemma(m.group(1)) != candidate.lower():
            continue
        body_end = h3[i + 1].start() if i + 1 < len(h3) else len(html)
        body = html[m.end() : body_end]
        text = re.sub(r"\s+", " ", _ANY_TAG.sub(" ", body)).strip()
        ets = _ETYM_RE.findall(text)
        return ets[0].strip() if ets else ""
    return None


# ---------------------------------------------------------------------------
# Categorization
# ---------------------------------------------------------------------------


def _matches_variant(cited: str, variants: frozenset[str]) -> bool:
    """Approximate-equality test between an etymology citation and the
    base noun's stem variants. Tolerates prefix overlap in either
    direction (e.g. *colac* citation matching base *colăci* stem)."""
    c = cited.lower()
    return any(c == v or c.startswith(v) or v.startswith(c) for v in variants)


def categorize_etymology(noun: str, etymology: str | None) -> EtymCategory:
    """Classify a DEX etymology string against the base noun.

    ``etymology`` may be ``None`` (no verb tree at all) or ``""`` (verb
    tree present but no etymology field). Both map to NO_ETYMOLOGY.
    """
    if etymology is None or etymology == "":
        return EtymCategory.NO_ETYMOLOGY

    et = etymology.lower()
    variants = stem_variants(noun)

    if _UNKNOWN_RE.search(et):
        return EtymCategory.UNKNOWN_ORIGIN
    if _FOREIGN_RE.search(et):
        return EtymCategory.FOREIGN_SOURCE

    # "vezi X" = "see X" — cross-reference to another lemma. Accept if
    # the cross-reference matches our noun.
    m = re.search(r"\bvezi\s+([a-zA-ZăâîșțĂÂÎȘȚ]+)", et)
    if m and _matches_variant(m.group(1), variants):
        return EtymCategory.CITES_OUR_NOUN

    # "în/îm/in/im + X" = prefixed denominal form.
    m = re.search(r"\b(?:în|îm|in|im)\s*\+\s*([a-zA-ZăâîșțĂÂÎȘȚ]+)", et)
    if m:
        return (
            EtymCategory.CITES_PREFIXED_FORM
            if _matches_variant(m.group(1), variants)
            else EtymCategory.CITES_DIFFERENT_WORD
        )

    # "A + X" prefix form (rarer, e.g. *avânt* < *A + vânt*).
    m = re.search(r"\bA\s*\+\s*([a-zA-ZăâîșțĂÂÎȘȚ]+)", etymology)
    if m:
        return (
            EtymCategory.CITES_PREFIXED_FORM
            if _matches_variant(m.group(1), variants)
            else EtymCategory.CITES_DIFFERENT_WORD
        )

    # Fall back to the first content word, skipping particles.
    m = re.search(r"\b([a-zA-ZăâîșțĂÂÎȘȚ]{2,})", et)
    if m:
        cited = m.group(1).lower()
        if cited in _PARTICLES:
            m2 = re.search(
                r"\b[a-zăâîșț]{2,}\s+" r"([a-zA-ZăâîșțĂÂÎȘȚ]{2,})",
                et,
            )
            if m2:
                cited = m2.group(1).lower()
        return (
            EtymCategory.CITES_OUR_NOUN
            if _matches_variant(cited, variants)
            else EtymCategory.CITES_DIFFERENT_WORD
        )

    return EtymCategory.NO_ETYMOLOGY


def is_accepted(category: EtymCategory) -> bool:
    """Whether an etymology category passes the audit."""
    return category in ACCEPT_CATEGORIES
