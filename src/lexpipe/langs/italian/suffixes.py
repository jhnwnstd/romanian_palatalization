"""Derivational suffixes of the Italian analyses (design §5; R20–R23).

Everything here is a pure function of lookup keys
(:func:`lexpipe.langs.italian.keys.match_key`: lowercase, respelling
accents stripped). The stages that use it are ``derivations`` (links
between a base and its derivatives), ``covariation`` (R23),
``spellout`` (R20) and ``noevidence`` (R21).

Inventory (:data:`SUFFIXES`)
----------------------------
Each :class:`Suffix` has a canonical name (``-ista``), the vowel at its
left edge (``i``, ``e`` or ``back``), the part of speech of its
derivatives, the surface bodies of its lemmas (gender variants:
``etto``/``etta``) and a group (design §5):

    learned    the Italian-internal learned set whose palatalization is
               said to track the base's masculine plural (the
               co-variation set): -issimo, -izzare, -ista, -ismo,
               -esimo, -ità, -izia, -ia (-ità/-izia/-ia only when the
               derivative is an Italian-internal formation; the stratum
               filter is applied by the analyses, R22)
    native     the never-palatalizing control: -ezza, -etto, -ino,
               -ello, -erello, -eria, -essa, -evole, -eggiare,
               derivational -ire, -ese, -iere
    other      reported, in neither set: -ificare, -ile (porcile)
    verbalizer the back-vowel verbalizer -are (no front vowel after the
               root; the non-palatalizing member of the spell-out pair)

``spellout`` gives the role of a verbalizer in the R20 pair "-izzare vs
-eggiare/-are/-ire": ``P`` for -izzare, ``N`` for the others.

Stem edge (:func:`read_edge`)
-----------------------------
A derivative is read at the junction of the base stem and the suffix
body (it_phon.json, operational definition: "palatalized = the derived
form shows c/g + i/e at the stem edge; retained = ch/gh"):

    palatal    stem + body, body front (grec + ista, bosc + etto)
    velar      stem + h + body (fuoc → foc + h + ista, antic + h + ità)
    back       stem + body for a back-vowel suffix (caric + are)
    truncated  the target is not at the junction: the -ico of the base
               is dropped (arcaico → arca + ismo, ceramica → ceram +
               ista); neither reading applies
    mismatch   the derivative does not end in stem + body

The stem may lose a mobile diphthong (fuoco → fochista, cieco →
cecità: :func:`stem_variants`). With ``allow_prefix`` the derivative
may carry a prefix (in + bosco + -ire → imboschire); candidate
generation never uses it.

Candidates (:func:`candidates`) and corpus forms (:func:`corpus_forms`)
-----------------------------------------------------------------------
R21: derivative candidates are generated in **both** spellings from a
base stem (lega + -ista → legista and leghista) and looked up in the
headwords and the corpus. -ia is not generated (c/g + ia is ambiguous
between /t͡ʃa/ and /t͡ʃi.a/, and Grecia, Lucia are names). A
derivative's corpus attestation sums the tokens of a small, documented
set of surface forms (:func:`corpus_forms`): the lemma and its regular
inflections (nouns and adjectives), or the infinitive and the four past
participles (verbs; present forms are left out because the 3sg of an
-are verb is often the base noun itself: amica).
"""

from __future__ import annotations

import unicodedata
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum
from types import MappingProxyType
from typing import Final

LEARNED: Final[str] = "learned"
NATIVE: Final[str] = "native"
OTHER: Final[str] = "other"
VERBALIZER: Final[str] = "verbalizer"
GROUPS: Final[tuple[str, ...]] = (LEARNED, NATIVE, OTHER, VERBALIZER)

VOWEL_I: Final[str] = "i"
VOWEL_E: Final[str] = "e"
VOWEL_BACK: Final[str] = "back"

# Inflection classes for corpus_forms().
_INFL_O4: Final[str] = "o4"  # -o lemma with -a/-i/-e (adjective-like)
_INFL_O2: Final[str] = "o2"  # -o noun with -i
_INFL_A2: Final[str] = "a2"  # -a feminine noun with -e
_INFL_A3: Final[str] = "a3"  # -ista: -a with -i/-e
_INFL_E2: Final[str] = "e2"  # -e with -i
_INFL_INV: Final[str] = "inv"  # invariable (-ità)
_INFL_VERB_ARE: Final[str] = "verb_are"
_INFL_VERB_ISC: Final[str] = "verb_isc"

# Truncation: at most this many letters of the base stem may be lost
# at the junction (-ic of -ico: arcaico → arcaismo).
_MAX_TRUNCATION: Final[int] = 3
_DIPHTHONGS: Final[tuple[tuple[str, str], ...]] = (("uo", "o"), ("ie", "e"))


class Edge(StrEnum):
    """Reading of a derivative at the stem edge (module docstring)."""

    PALATAL = "palatal"
    VELAR = "velar"
    BACK = "back"
    TRUNCATED = "truncated"
    MISMATCH = "mismatch"


READABLE_EDGES: Final[frozenset[Edge]] = frozenset({Edge.PALATAL, Edge.VELAR})


@dataclass(frozen=True, slots=True)
class Suffix:
    """One derivational suffix (see the module docstring)."""

    name: str
    group: str
    vowel: str
    pos: str
    bodies: tuple[str, ...]
    inflection: str
    generate: bool = True
    spellout: str = ""

    def __post_init__(self) -> None:
        if self.group not in GROUPS:
            raise ValueError(f"{self.name}: unknown group {self.group!r}")
        if self.vowel not in (VOWEL_I, VOWEL_E, VOWEL_BACK):
            raise ValueError(f"{self.name}: unknown vowel {self.vowel!r}")
        if not self.name.startswith("-") or not self.bodies:
            raise ValueError(f"malformed suffix {self.name!r}")

    @property
    def front(self) -> bool:
        return self.vowel != VOWEL_BACK


SUFFIXES: Final[tuple[Suffix, ...]] = (
    # learned /i/-initial set (§5; -esimo is the -ismo variant)
    Suffix("-issimo", LEARNED, VOWEL_I, "adj", ("issimo", "issima"), _INFL_O4),
    Suffix(
        "-izzare",
        LEARNED,
        VOWEL_I,
        "verb",
        ("izzare",),
        _INFL_VERB_ARE,
        spellout="P",
    ),
    Suffix("-ista", LEARNED, VOWEL_I, "noun", ("ista",), _INFL_A3),
    Suffix("-ismo", LEARNED, VOWEL_I, "noun", ("ismo",), _INFL_O2),
    Suffix("-esimo", LEARNED, VOWEL_E, "noun", ("esimo",), _INFL_O2),
    Suffix("-ità", LEARNED, VOWEL_I, "noun", ("ità",), _INFL_INV),
    Suffix("-izia", LEARNED, VOWEL_I, "noun", ("izia",), _INFL_A2),
    Suffix("-ia", LEARNED, VOWEL_I, "noun", ("ia",), _INFL_A2, False),
    # native never-palatalizing control (§5)
    Suffix("-ezza", NATIVE, VOWEL_E, "noun", ("ezza",), _INFL_A2),
    Suffix("-etto", NATIVE, VOWEL_E, "noun", ("etto", "etta"), _INFL_O4),
    Suffix("-ino", NATIVE, VOWEL_I, "noun", ("ino", "ina"), _INFL_O4),
    Suffix("-erello", NATIVE, VOWEL_E, "noun", ("erello", "erella"), _INFL_O4),
    Suffix("-ello", NATIVE, VOWEL_E, "noun", ("ello", "ella"), _INFL_O4),
    Suffix("-eria", NATIVE, VOWEL_E, "noun", ("eria",), _INFL_A2),
    Suffix("-essa", NATIVE, VOWEL_E, "noun", ("essa",), _INFL_A2),
    Suffix("-evole", NATIVE, VOWEL_E, "adj", ("evole",), _INFL_E2),
    Suffix(
        "-eggiare",
        NATIVE,
        VOWEL_E,
        "verb",
        ("eggiare",),
        _INFL_VERB_ARE,
        spellout="N",
    ),
    Suffix(
        "-ire", NATIVE, VOWEL_I, "verb", ("ire",), _INFL_VERB_ISC, spellout="N"
    ),
    Suffix("-ese", NATIVE, VOWEL_E, "noun", ("ese",), _INFL_E2),
    Suffix(
        "-iere", NATIVE, VOWEL_I, "noun", ("iere", "iera", "iero"), _INFL_E2
    ),
    # reported, in neither set
    Suffix("-ificare", OTHER, VOWEL_I, "verb", ("ificare",), _INFL_VERB_ARE),
    Suffix("-ile", OTHER, VOWEL_I, "noun", ("ile",), _INFL_E2),
    # the back-vowel verbalizer of the spell-out pair (R20)
    Suffix(
        "-are",
        VERBALIZER,
        VOWEL_BACK,
        "verb",
        ("are",),
        _INFL_VERB_ARE,
        spellout="N",
    ),
)

BY_NAME: Final[Mapping[str, Suffix]] = MappingProxyType(
    {suffix.name: suffix for suffix in SUFFIXES}
)
LEARNED_SET: Final[tuple[str, ...]] = tuple(
    s.name for s in SUFFIXES if s.group == LEARNED
)
NATIVE_SET: Final[tuple[str, ...]] = tuple(
    s.name for s in SUFFIXES if s.group == NATIVE
)
# Template spellings of gender variants and plurals → canonical name.
_ALIASES: Final[Mapping[str, str]] = MappingProxyType(
    {
        "-issima": "-issimo",
        "-etta": "-etto",
        "-ina": "-ino",
        "-ella": "-ello",
        "-erella": "-erello",
        "-iera": "-iere",
        "-iero": "-iere",
        "-isti": "-ista",
    }
)


def _nfc(text: str) -> str:
    return unicodedata.normalize("NFC", text).strip()


def canonical(raw: str) -> Suffix | None:
    """The suffix a template names (``ezza``, ``-etta``, ``-ità``), or
    None when it is not in the inventory."""
    text = _nfc(raw).lower()
    if not text:
        return None
    name = text if text.startswith("-") else "-" + text
    name = _ALIASES.get(name, name)
    return BY_NAME.get(name)


def stem_variants(stem: str) -> list[str]:
    """*stem* and its form without the last mobile diphthong (fuoc →
    foc, ciec → cec), in that order."""
    out = [stem]
    for long, short in _DIPHTHONGS:
        cut = stem.rfind(long)
        if cut > 0:
            out.append(stem[:cut] + short + stem[cut + len(long) :])
    return list(dict.fromkeys(out))


def _matches(key: str, expected: str, allow_prefix: bool) -> bool:
    if key == expected:
        return True
    return allow_prefix and len(key) > len(expected) and key.endswith(expected)


def _truncated(stem: str, key: str, body: str, allow_prefix: bool) -> bool:
    if not key.endswith(body) or len(key) <= len(body):
        return False
    left = key[: -len(body)]
    for length in range(1, _MAX_TRUNCATION + 1):
        kept = stem[:-length]
        if len(kept) < 2:
            break
        if left == kept or (allow_prefix and left.endswith(kept)):
            return True
        # haplology of an i-initial body after -i- (arcai + ismo)
        if body.startswith("i") and left + "i" == kept:
            return True
    return False


def read_edge(
    stem: str,
    suffix: Suffix,
    key: str,
    *,
    allow_prefix: bool = False,
) -> tuple[Edge, str]:
    """``(edge, body)``: how *key* realizes *stem* + *suffix* (see the
    module docstring); *body* is the matched surface body ("" when the
    edge is truncated or a mismatch)."""
    bodies = sorted(suffix.bodies, key=lambda b: (-len(b), b))
    for body in bodies:
        for variant in stem_variants(stem):
            if not suffix.front:
                if _matches(key, variant + body, allow_prefix):
                    return Edge.BACK, body
                continue
            if _matches(key, variant + "h" + body, allow_prefix):
                return Edge.VELAR, body
            if _matches(key, variant + body, allow_prefix):
                return Edge.PALATAL, body
    for body in bodies:
        for variant in stem_variants(stem):
            if _truncated(variant, key, body, allow_prefix):
                return Edge.TRUNCATED, ""
    return Edge.MISMATCH, ""


def best_reading(
    stem: str,
    key: str,
    suffixes: Iterable[Suffix] = SUFFIXES,
    *,
    allow_prefix: bool = False,
) -> tuple[Suffix, Edge, str] | None:
    """The suffix whose longest body *key* realizes after *stem*
    (palatal, velar or back reading only), or None. Used for derived
    lists, which name a derivative without its suffix."""
    found: list[tuple[int, str, Suffix, Edge, str]] = []
    for suffix in suffixes:
        edge, body = read_edge(stem, suffix, key, allow_prefix=allow_prefix)
        if edge in (Edge.PALATAL, Edge.VELAR, Edge.BACK):
            found.append((-len(body), suffix.name, suffix, edge, body))
    if not found:
        return None
    found.sort(key=lambda item: (item[0], item[1]))
    _, _, suffix, edge, body = found[0]
    return suffix, edge, body


@dataclass(frozen=True, slots=True)
class Candidate:
    """A generated derivative (R21)."""

    form: str
    spelling: Edge  # palatal | velar | back
    suffix: Suffix
    body: str


def candidates(stem: str, suffix: Suffix) -> list[Candidate]:
    """Both spellings of *stem* + *suffix* for every body and stem
    variant (empty for a suffix that is not generated)."""
    if not suffix.generate:
        return []
    out: list[Candidate] = []
    for body in suffix.bodies:
        for variant in stem_variants(stem):
            if not suffix.front:
                out.append(Candidate(variant + body, Edge.BACK, suffix, body))
                continue
            out.append(Candidate(variant + body, Edge.PALATAL, suffix, body))
            out.append(
                Candidate(variant + "h" + body, Edge.VELAR, suffix, body)
            )
    seen: dict[str, Candidate] = {}
    for item in out:
        seen.setdefault(item.form, item)
    return list(seen.values())


def other_spelling(form: str, edge: Edge, body: str) -> str:
    """The form with the other reading at the junction (leghista ↔
    legista); "" when *edge* is not palatal or velar."""
    if edge not in READABLE_EDGES or not body or not form.endswith(body):
        return ""
    left = form[: -len(body)]
    if edge == Edge.VELAR:
        return left[:-1] + body if left.endswith("h") else ""
    return left + "h" + body


def _swap_ending(word: str, old: str, news: Sequence[str]) -> list[str]:
    if not word.endswith(old):
        return [word]
    base = word[: -len(old)]
    return [word] + [base + new for new in news]


def corpus_forms(word: str, suffix: Suffix) -> list[str]:
    """Surface forms whose tokens attest a derivative lemma (module
    docstring), in a fixed order without duplicates."""
    infl = suffix.inflection
    forms: list[str]
    if infl == _INFL_VERB_ARE and word.endswith("are"):
        stem = word[:-3]
        forms = [word] + [stem + end for end in ("ato", "ata", "ati", "ate")]
    elif infl == _INFL_VERB_ISC and word.endswith("ire"):
        stem = word[:-3]
        forms = [word] + [stem + end for end in ("ito", "ita", "iti", "ite")]
    elif infl == _INFL_O4 and word.endswith("o"):
        forms = _swap_ending(word, "o", ("a", "i", "e"))
    elif infl == _INFL_O4 and word.endswith("a"):
        forms = _swap_ending(word, "a", ("e",))
    elif infl == _INFL_O2:
        forms = _swap_ending(word, "o", ("i",))
    elif infl == _INFL_A2:
        forms = _swap_ending(word, "a", ("e",))
    elif infl == _INFL_A3:
        forms = _swap_ending(word, "a", ("i", "e"))
    elif infl == _INFL_E2 and word.endswith("e"):
        forms = _swap_ending(word, "e", ("i",))
    elif infl == _INFL_E2 and word.endswith("a"):
        forms = _swap_ending(word, "a", ("e",))
    elif infl == _INFL_E2 and word.endswith("o"):
        forms = _swap_ending(word, "o", ("i",))
    else:
        forms = [word]
    return list(dict.fromkeys(forms))
