"""Orthography and IPA normalization for the language plugins (R8).

The Romanian helpers (``wiktionary_normalizer``) must not be reused:
they fold accents and vowel length (žák → žak, diplomaté → diplomate,
città → citta), reject every bracketed transcription and strip stress
and syllable dots. This module does the opposite: it keeps every
distinction and only canonicalises spelling variants of the same
thing.

Orthography
-----------
- :func:`normalize_orthography`: NFC, surrounding whitespace removed,
  internal runs of whitespace collapsed to one space, optionally
  lowercased. Nothing is folded: Czech á é í ó ú ů ý stay (diplomaté
  ≠ diplomate), Italian final accents stay (città, perché).
- :func:`strip_respelling_accents` (Italian, **match keys only**):
  English Wiktionary respells Italian forms with a grave or acute
  accent on the stressed vowel (fòcora, àntico, nómora). Such accents
  on a vowel that is not word-final are removed; word-final accents
  are orthographic and kept (città, caffè, perché).
- :meth:`Normalizer.match_key`: the per-language lookup key: the
  lowercased orthography, the typographic apostrophe ’ written ', and
  for Italian the respelling accents stripped.

IPA
---
:func:`parse_ipa` accepts one transcription delimited as ``/…/``
(phonemic) or ``[…]`` (phonetic) and records which in
:class:`IpaKind`; anything else (a bare string, a reconstruction
``*/…/``, two transcriptions in one string) gives None. The body is
kept as written, stress marks ``ˈ ˌ`` and syllable dots included,
except that

  - an ASCII ``:`` becomes the length mark ``ː`` (it.wiktionary
    writes ``/ˈɡrɛ:ko/``);
  - the affricate ligatures ʧ ʤ ʦ ʣ and the undertie t͜ʃ are written
    with the tie bar: t͡ʃ d͡ʒ t͡s d͡z;
  - optionally (:func:`tie_affricates`), directly adjacent t+ʃ, d+ʒ,
    t+s, d+z get the tie bar, as the Romanian normalizer does. Pairs
    separated by anything (a syllable dot: ``/ˈɡrat.sja/``, a stress
    mark, a diacritic) are left alone.

Whether adjacent pairs are tied is a per-language convention
(:data:`NORMALIZERS`). Italian ties them: en.wiktionary's it-pr
syllabifies with dots, so an adjacent pair inside a syllable is an
affricate (``/ski.tso.freˈni.a/``). Czech does not: cs-IPA already
writes c and č as t͡s and t͡ʃ, and its untied pairs are clusters across
a prefix boundary ([ˈpodzɪm] pod+zim, [ˈotstavɛt͡s] od+stavec,
[ˈnatʃɛɲiː] nad+šení), which tying would merge.

Syllables and stress
--------------------
:func:`syllabify` splits a transcription at syllable dots, before
stress marks (which stay on their syllable) and at word boundaries
(spaces, ``‿``). :func:`is_syllabified` checks that every piece has
exactly one nucleus (a run of vowels, or a consonant marked syllabic
with ◌̩/◌̍; a vowel marked non-syllabic with ◌̯ is a glide): Czech
transcriptions carry no dots, so ``[ˈabort]`` is one piece with two
nuclei. :func:`stress_from_end` (1 final, 2 penultimate, 3
antepenultimate, 4 for anything earlier) and :func:`post_tonic` give
None for a transcription that is not syllabified, has no primary
stress or spans several words.
"""

from __future__ import annotations

import re
import unicodedata
from collections.abc import Mapping
from dataclasses import dataclass
from enum import StrEnum
from types import MappingProxyType
from typing import Final

PRIMARY_STRESS: Final[str] = "ˈ"
SECONDARY_STRESS: Final[str] = "ˌ"
SYLLABLE_DOT: Final[str] = "."
LENGTH: Final[str] = "ː"
TIE: Final[str] = "\u0361"
_UNDERTIE: Final[str] = "\u035c"

# Affricates whose two halves get the tie bar.
AFFRICATES: Final[tuple[tuple[str, str], ...]] = (
    ("t", "ʃ"),
    ("d", "ʒ"),
    ("t", "s"),
    ("d", "z"),
)
_LIGATURES: Final[Mapping[str, str]] = MappingProxyType(
    {"ʧ": f"t{TIE}ʃ", "ʤ": f"d{TIE}ʒ", "ʦ": f"t{TIE}s", "ʣ": f"d{TIE}z"}
)
_ADJACENT_RE: Final[re.Pattern[str]] = re.compile(
    "|".join(f"(?<={a}){b}" for a, b in AFFRICATES)
)
_UNDERTIED_RE: Final[re.Pattern[str]] = re.compile(
    "|".join(f"(?<={a}){_UNDERTIE}(?={b})" for a, b in AFFRICATES)
)

# Respelling accents of Italian Wiktionary headwords.
_GRAVE_ACUTE: Final[frozenset[str]] = frozenset({"\u0300", "\u0301"})
_ITALIAN_VOWELS: Final[frozenset[str]] = frozenset("aeiouAEIOU")
_APOSTROPHES: Final[Mapping[int, str]] = MappingProxyType({0x2019: "'"})

# IPA vowel letters (the nuclei of syllable counting).
_VOWELS: Final[frozenset[str]] = frozenset("iyɨʉɯuɪʏʊeøɘɵɤoəɛœɜɞʌɔæɐaɶɑɒɚɝ")
_NON_SYLLABIC: Final[str] = "\u032f"
_SYLLABIC: Final[frozenset[str]] = frozenset({"\u0329", "\u030d"})
_WORD_BREAK_RE: Final[re.Pattern[str]] = re.compile(r"[\s‿]+")
_MAX_STRESS_FROM_END: Final[int] = 4

STRESS_CLASS_NAMES: Final[Mapping[int, str]] = MappingProxyType(
    {
        1: "oxytone",
        2: "paroxytone",
        3: "proparoxytone",
        4: "preproparoxytone",
    }
)


# ---------------------------------------------------------------------------
# Orthography
# ---------------------------------------------------------------------------


def normalize_orthography(text: str, *, lowercase: bool = False) -> str:
    """NFC, whitespace tidied, optionally lowercased; nothing folded."""
    tidy = " ".join(unicodedata.normalize("NFC", text).split())
    return tidy.lower() if lowercase else tidy


def _next_base(chars: list[str], start: int) -> str:
    """The first non-combining character at or after *start* ("" at end)."""
    for ch in chars[start:]:
        if not unicodedata.combining(ch):
            return ch
    return ""


def strip_respelling_accents(text: str) -> str:
    """Remove grave/acute accents from vowels that are not word-final.

    For Italian match keys only (see the module docstring):
    ``fòcora`` → ``focora``, ``àntico`` → ``antico``, while
    ``città``, ``perché`` and ``caffè americano`` are unchanged.
    """
    chars = list(unicodedata.normalize("NFD", text))
    out: list[str] = []
    base = ""
    for i, ch in enumerate(chars):
        if not unicodedata.combining(ch):
            base = ch
        elif (
            ch in _GRAVE_ACUTE
            and base in _ITALIAN_VOWELS
            and _next_base(chars, i + 1).isalpha()
        ):
            continue  # a respelling accent inside the word
        out.append(ch)
    return unicodedata.normalize("NFC", "".join(out))


# ---------------------------------------------------------------------------
# IPA
# ---------------------------------------------------------------------------


class IpaKind(StrEnum):
    """Delimiters of a transcription: ``/…/`` or ``[…]``."""

    PHONEMIC = "phonemic"
    PHONETIC = "phonetic"


_DELIMITERS: Final[Mapping[IpaKind, tuple[str, str]]] = MappingProxyType(
    {IpaKind.PHONEMIC: ("/", "/"), IpaKind.PHONETIC: ("[", "]")}
)


@dataclass(frozen=True, slots=True)
class Ipa:
    """A normalized transcription: its body (no delimiters) and kind."""

    text: str
    kind: IpaKind

    @property
    def delimited(self) -> str:
        opening, closing = _DELIMITERS[self.kind]
        return f"{opening}{self.text}{closing}"


def tie_affricates(ipa: str) -> str:
    """Tie directly adjacent t+ʃ, d+ʒ, t+s, d+z (t͡ʃ d͡ʒ t͡s d͡z).

    Pairs that already carry a tie, or are separated by a syllable
    dot, a stress mark or a diacritic, are left alone.
    """
    return _ADJACENT_RE.sub(lambda m: TIE + m.group(0), ipa)


def _canonical_ties(ipa: str) -> str:
    for ligature, tied in _LIGATURES.items():
        ipa = ipa.replace(ligature, tied)
    return _UNDERTIED_RE.sub(TIE, ipa)


def _kind(text: str) -> IpaKind | None:
    for kind, (opening, closing) in _DELIMITERS.items():
        if len(text) > 2 and text[0] == opening and text[-1] == closing:
            inner = text[1:-1]
            if opening not in inner and closing not in inner:
                return kind
    return None


def parse_ipa(raw: str, *, tie: bool = True) -> Ipa | None:
    """Normalize one delimited transcription (see the module docstring).

    *tie* adds tie bars to adjacent affricate halves
    (:func:`tie_affricates`). Returns None for anything that is not a
    single ``/…/`` or ``[…]`` transcription.
    """
    text = unicodedata.normalize("NFC", raw).strip()
    kind = _kind(text)
    if kind is None:
        return None
    body = _canonical_ties(text[1:-1].strip().replace(":", LENGTH))
    if tie:
        body = tie_affricates(body)
    return Ipa(unicodedata.normalize("NFC", body), kind) if body else None


def _body(ipa: Ipa | str) -> str:
    if isinstance(ipa, Ipa):
        return ipa.text
    text = unicodedata.normalize("NFC", ipa).strip()
    return text[1:-1] if _kind(text) is not None else text


def syllabify(ipa: Ipa | str) -> list[str]:
    """Syllables of a transcription, stress marks kept on their syllable.

    *ipa* is an :class:`Ipa`, a delimited string or a bare body.
    Splits at syllable dots, before ``ˈ``/``ˌ`` and at word
    boundaries; empty pieces are dropped. Check
    :func:`is_syllabified` before trusting the result.
    """
    syllables: list[str] = []
    for word in _WORD_BREAK_RE.split(_body(ipa)):
        current = ""
        for ch in word:
            if ch == SYLLABLE_DOT or ch in (PRIMARY_STRESS, SECONDARY_STRESS):
                if current:
                    syllables.append(current)
                current = "" if ch == SYLLABLE_DOT else ch
            else:
                current += ch
        if current:
            syllables.append(current)
    return syllables


def nuclei(syllable: str) -> int:
    """Number of syllable nuclei in a piece of transcription."""
    count = 0
    in_vowel_run = False
    chars = list(unicodedata.normalize("NFD", syllable))
    for i, ch in enumerate(chars):
        if ch in _SYLLABIC:
            count += 1
            in_vowel_run = False
        elif unicodedata.combining(ch) or ch in (LENGTH, "ˑ"):
            continue
        elif ch in _VOWELS and _next_mark_is_not(chars, i + 1):
            if not in_vowel_run:
                count += 1
            in_vowel_run = True
        else:
            in_vowel_run = False
    return count


def _next_mark_is_not(chars: list[str], start: int) -> bool:
    """False if the combining marks after *start* include ◌̯."""
    for ch in chars[start:]:
        if not unicodedata.combining(ch):
            break
        if ch == _NON_SYLLABIC:
            return False
    return True


def is_syllabified(ipa: Ipa | str) -> bool:
    """True if every piece :func:`syllabify` returns has one nucleus."""
    pieces = syllabify(ipa)
    return bool(pieces) and all(nuclei(piece) == 1 for piece in pieces)


def _stressed_index(ipa: Ipa | str) -> tuple[list[str], int] | None:
    body = _body(ipa)
    if _WORD_BREAK_RE.search(body.strip()) or not is_syllabified(body):
        return None
    syllables = syllabify(body)
    stressed = [
        i for i, s in enumerate(syllables) if s.startswith(PRIMARY_STRESS)
    ]
    if not stressed:
        return None
    return syllables, stressed[-1]


def stress_from_end(ipa: Ipa | str) -> int | None:
    """Position of the primary stress counted from the end (1 = final).

    4 means the fourth syllable from the end or earlier. None when the
    transcription is not syllabified, has no ``ˈ`` or spans several
    words: ``/ˈfwɔ.ko/`` → 2, ``/ˈmɛ.di.ko/`` → 3, ``[ˈkluk]`` → 1.
    """
    found = _stressed_index(ipa)
    if found is None:
        return None
    syllables, index = found
    return min(len(syllables) - index, _MAX_STRESS_FROM_END)


def post_tonic(ipa: Ipa | str, target_syllable_index: int) -> bool | None:
    """Is the target syllable the one right after the stressed one?

    *target_syllable_index* indexes :func:`syllabify`'s list; negative
    values count from the end (-1 is the final syllable, where the
    stem-final velar of -co/-go nouns sits). An index outside the
    word raises :class:`IndexError`; the same None cases as
    :func:`stress_from_end` apply. ``/ˈfwɔ.ko/``, -1 → True;
    ``/ˈmɛ.di.ko/``, -1 → False.
    """
    found = _stressed_index(ipa)
    if found is None:
        return None
    syllables, stressed = found
    count = len(syllables)
    if not -count <= target_syllable_index < count:
        raise IndexError(
            f"syllable {target_syllable_index} of a {count}-syllable word"
        )
    return target_syllable_index % count == stressed + 1


# ---------------------------------------------------------------------------
# Per-language normalizers
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Normalizer:
    """One language's normalization conventions (see the module docstring)."""

    code: str
    strip_respelling: bool = False
    tie: bool = True

    def orthography(self, text: str, *, lowercase: bool = False) -> str:
        return normalize_orthography(text, lowercase=lowercase)

    def match_key(self, text: str) -> str:
        """Lookup key: lowercase, ’ → ', respelling accents stripped
        where the language uses them."""
        key = normalize_orthography(text, lowercase=True)
        key = key.translate(_APOSTROPHES)
        return strip_respelling_accents(key) if self.strip_respelling else key

    def ipa(self, raw: str) -> Ipa | None:
        return parse_ipa(raw, tie=self.tie)


NORMALIZERS: Final[Mapping[str, Normalizer]] = MappingProxyType(
    {
        "it": Normalizer("it", strip_respelling=True, tie=True),
        "cs": Normalizer("cs", strip_respelling=False, tie=False),
    }
)


def normalizer_for(code: str) -> Normalizer:
    """The :class:`Normalizer` of a registered language code."""
    try:
        return NORMALIZERS[code]
    except KeyError:
        known = ", ".join(sorted(NORMALIZERS))
        raise KeyError(f"no normalizer for {code!r}; known: {known}") from None
