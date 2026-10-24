"""Lookup keys for Italian word forms (design rule R8).

Orthography is never folded for storage: forms are kept as written
(NFC), final accents included (città, perché). For *matching* — the
form index, Morph-it!, the it.wiktionary edition and the corpus —
two keys are used, both built with :mod:`lexpipe.normalize`:

  - :func:`match_key`: lowercase, typographic apostrophe ’ written ',
    and English Wiktionary's respelling accents on non-final vowels
    stripped (càrichi → carichi, fòcora → focora);
  - :func:`case_key`: the same with the case kept, for case-sensitive
    lookups (Medici the family vs medici the doctors).
"""

from __future__ import annotations

from typing import Final

from lexpipe.normalize import (
    Normalizer,
    normalize_orthography,
    normalizer_for,
    strip_respelling_accents,
)

LANGUAGE_CODE: Final[str] = "it"
_APOSTROPHES: Final[dict[int, str]] = {0x2019: "'"}


def normalizer() -> Normalizer:
    """The Italian :class:`lexpipe.normalize.Normalizer`."""
    return normalizer_for(LANGUAGE_CODE)


def match_key(text: str) -> str:
    """Lowercase lookup key (see the module docstring)."""
    return normalizer().match_key(text)


def case_key(text: str) -> str:
    """Case-preserving lookup key (see the module docstring)."""
    tidy = normalize_orthography(text).translate(_APOSTROPHES)
    return strip_respelling_accents(tidy)


def form_keys(text: str) -> tuple[str, str]:
    """``(match_key(text), case_key(text))``."""
    return match_key(text), case_key(text)


def strip_final_vowel(key: str) -> str:
    """*key* without one final unaccented vowel (amico → amic), if it
    has more than two letters; used for root keys."""
    if len(key) > 2 and key[-1] in "aeiou":
        return key[:-1]
    return key
