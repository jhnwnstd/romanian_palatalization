"""Stem-final targets and plural alternations of Italian (design §5).

Everything works on lookup keys (:func:`lexpipe.langs.italian.keys.
match_key`: lowercase, respelling accents stripped). Orthography is a
reliable guide for native words: ‹c g› before ‹e i› are /t͡ʃ d͡ʒ/, ‹ch
gh› before ‹e i› are /k ɡ/, ‹sc›+‹e i› is /ʃʃ/. Loans that break this
(chip /t͡ʃip/) are invariable and so never reach the alternation
columns.

Lemma (:func:`analyze_lemma`) → :class:`StemAnalysis`
-----------------------------------------------------
Checked in this order:

  1. -chio/-chia/-ghio/-ghia (occhio, macchia, unghia): velar + root
     glide, ``nde_class = occhio``, cluster ``chi``/``ghi``;
  2. singular -che/-ghe (psiche): ``nde_class = psiche``, cluster
     ``che``/``ghe``;
  3. c/g (also sc, cc, gg) followed by a front vowel letter at the end
     (-ce/-ge/-sce, -cia/-gia/-scia, -cio/-gio/-scio, -cie…: luce,
     valigia, bacio, farmacia, specie): lexically palatal,
     ``frontstem``, not a target;
  4. c/g before the final vowel (amico, parco, belga, bosco, sacco):
     the dorsal target; ``cluster`` is ``sc``, ``cc``, ``gg`` or the
     consonant before the velar plus the velar (``nc``, ``ng``, ``rc``,
     ``rg``, ``lc``, ``lg``…), "" after a vowel;
  5. t/d/s/z before the final vowel (gatto, naso, prezzo): a coronal,
     recorded for the "no coronal rule" rows; ``cluster`` likewise
     (``st``, ``nt``, ``tt``…);
  6. anything else: no target (``stem_final`` "").

``stem_final`` is c/g for 1–4 and t/d/s/z for 5.

Plural form (:func:`analyze_form`) → :class:`FormAnalysis`
---------------------------------------------------------
``opportunity`` from the plural's ending: ``i``, ``e``, ``a`` (back
vowel: uova, Old Italian fòcora) or ``none`` (invariable — the form is
the lemma —, a consonant or accented ending, or a form whose stem edge
does not match the lemma's: suppletive or irregular).

``mutation`` for dorsal targets with opportunity i/e (§5 table):

    palatalized (True): co→ci go→gi ca→ci ga→gi sco→sci cco→cci ggo→ggi
    velar kept (False): co→chi go→ghi ca→che ga→ghe ca→chi ga→ghi
                        sco→schi cco→cchi

with ``palatal_output`` t͡ʃ, d͡ʒ, ʃ (from sc), t͡ʃː (cc) or d͡ʒː (gg).
Coronals before i/e never change (False: to→ti); NDE stems keep the
velar (occhio→occhi, psiche→psichi); frontstem forms keep their palatal
(cia→ce is spelling only): all False. ``orth_change`` records the
lemma and plural endings at the target (``co→chi``).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Final

DORSAL: Final[str] = "dorsal"
CORONAL: Final[str] = "coronal"
NDE_OCCHIO: Final[str] = "occhio"
NDE_PSICHE: Final[str] = "psiche"
NDE_GIMPE: Final[str] = "gimpe_like"

OPPORTUNITY_I: Final[str] = "i"
OPPORTUNITY_E: Final[str] = "e"
OPPORTUNITY_A: Final[str] = "a"
OPPORTUNITY_NONE: Final[str] = "none"

PALATAL_OUTPUT: Final[dict[str, str]] = {
    "c": "t͡ʃ",
    "g": "d͡ʒ",
    "sc": "ʃ",
    "cc": "t͡ʃː",
    "gg": "d͡ʒː",
}

_VOWELS: Final[str] = "aeiou"
_CORONALS: Final[str] = "tdsz"
_OCCHIO_RE: Final[re.Pattern[str]] = re.compile(r"(c|g)hi[ao]\Z")
_PSICHE_RE: Final[re.Pattern[str]] = re.compile(r"(c|g)he\Z")
_FRONTSTEM_RE: Final[re.Pattern[str]] = re.compile(
    r"(sc|cc|gg|c|g)(?:i[aeiou]?|e)\Z"
)
_GIMPE_RE: Final[re.Pattern[str]] = re.compile(r"(?:c|g|ch|gh|sc)(?:i?e|i)\Z")
_ROOT_VELAR_FRONT_RE: Final[re.Pattern[str]] = re.compile(r"[cg]h[ie]")


@dataclass(frozen=True, slots=True)
class StemAnalysis:
    """The stem edge of a lemma (see the module docstring)."""

    stem_final: str
    cluster: str
    target: str
    frontstem: bool
    nde_class: str
    stem: str
    final_vowel: str  # the letters after the stem ("ia" for valigia)

    @property
    def dorsal(self) -> bool:
        return self.target == DORSAL

    @property
    def target_letters(self) -> str:
        """The letters that alternate: sc, cc, gg, or the velar."""
        if self.cluster in ("sc", "cc", "gg"):
            return self.cluster
        return self.stem_final


NO_STEM: Final[StemAnalysis] = StemAnalysis("", "", "", False, "", "", "")


@dataclass(frozen=True, slots=True)
class FormAnalysis:
    """One plural form against its lemma."""

    opportunity: str
    mutation: bool | None
    orth_change: str
    palatal_output: str
    aligned: bool


def _cluster(stem: str, letter: str) -> str:
    """``sc``/``cc``/``gg`` or consonant + target letter before it."""
    if len(stem) < 2:
        return ""
    before = stem[-2]
    if before == "s" and letter == "c":
        return "sc"
    if before == letter:
        return letter * 2
    if before.isalpha() and before not in _VOWELS:
        return before + letter
    return ""


def analyze_lemma(key: str) -> StemAnalysis:
    """Stem analysis of a lemma key (see the module docstring)."""
    if not key or not key.isalpha():
        return NO_STEM
    if (match := _OCCHIO_RE.search(key)) is not None:
        letter = match.group(1)
        return StemAnalysis(
            letter, f"{letter}hi", DORSAL, False, NDE_OCCHIO, key[:-1], key[-1]
        )
    if (match := _PSICHE_RE.search(key)) is not None:
        letter = match.group(1)
        return StemAnalysis(
            letter, f"{letter}he", DORSAL, False, NDE_PSICHE, key[:-1], "e"
        )
    if (match := _FRONTSTEM_RE.search(key)) is not None and len(key) > 2:
        target = match.group(1)
        letter = target[-1]
        stem = key[: match.end(1)]
        cluster = target if len(target) == 2 else _cluster(stem, letter)
        return StemAnalysis(
            letter, cluster, "", True, "", stem, key[match.end(1) :]
        )
    if len(key) < 2 or key[-1] not in _VOWELS:
        return NO_STEM
    stem, final_vowel = key[:-1], key[-1]
    letter = stem[-1]
    if letter in "cg":
        return StemAnalysis(
            letter,
            _cluster(stem, letter),
            DORSAL,
            False,
            "",
            stem,
            final_vowel,
        )
    if letter in _CORONALS:
        return StemAnalysis(
            letter,
            _cluster(stem, letter),
            CORONAL,
            False,
            "",
            stem,
            final_vowel,
        )
    return StemAnalysis("", "", "", False, "", stem, final_vowel)


def opportunity_of(lemma_key: str, form_key: str) -> str:
    """i | e | a | none from a plural's ending (not checking the stem)."""
    if not form_key or form_key == lemma_key:
        return OPPORTUNITY_NONE
    last = form_key[-1]
    if last == "i":
        return OPPORTUNITY_I
    if last == "e":
        return OPPORTUNITY_E
    if last == "a":
        return OPPORTUNITY_A
    return OPPORTUNITY_NONE


def _edge(stem: str, size: int = 3) -> str:
    return stem[-size:] if len(stem) > size else stem


def _aligned(stem: str, body: str) -> bool:
    """Does *body* (the plural minus its ending) end like the lemma's
    stem? Only the last three letters count, so capo- compounds
    (capocuoco → capicuochi) and mobile diphthongs further left align."""
    return bool(stem) and body.endswith(_edge(stem))


def _change(lemma_key: str, form_key: str, site: str) -> str:
    """``lemma ending→form ending`` from the last occurrence of *site*
    (the alternating letters) in each: co→chi, gia→ge, chio→chi."""
    left = lemma_key.rfind(site)
    right = form_key.rfind(site)
    if left < 0 or right < 0:
        return ""
    return f"{lemma_key[left:]}→{form_key[right:]}"


def _dorsal(
    stem: StemAnalysis, lemma_key: str, form_key: str, opportunity: str
) -> FormAnalysis:
    letters = stem.target_letters
    vowel = form_key[-1]
    if form_key.endswith(letters + "h" + vowel) and _aligned(
        stem.stem, form_key[:-2]
    ):
        change = _change(lemma_key, form_key, letters)
        return FormAnalysis(opportunity, False, change, "", True)
    if form_key.endswith(letters + vowel) and _aligned(
        stem.stem, form_key[:-1]
    ):
        output = PALATAL_OUTPUT.get(letters, PALATAL_OUTPUT[stem.stem_final])
        change = _change(lemma_key, form_key, letters)
        return FormAnalysis(opportunity, True, change, output, True)
    return FormAnalysis(OPPORTUNITY_NONE, None, "", "", False)


def _nde(
    stem: StemAnalysis, lemma_key: str, form_key: str, opportunity: str
) -> FormAnalysis:
    """occhio → occhi, macchia → macchie, psiche → psichi: velar kept."""
    site = stem.stem_final + "h"
    edge = _edge(stem.stem)
    aligned = form_key.endswith(edge) or form_key[:-1].endswith(edge)
    if site in form_key[-4:] and aligned:
        change = _change(lemma_key, form_key, site)
        return FormAnalysis(opportunity, False, change, "", True)
    return FormAnalysis(OPPORTUNITY_NONE, None, "", "", False)


def _plain(
    stem: StemAnalysis, lemma_key: str, form_key: str, opportunity: str
) -> FormAnalysis:
    """Coronal and frontstem forms: the consonant does not change
    (gatto → gatti; valigia → valigie/valige, bacio → baci)."""
    letters = stem.target_letters
    right = form_key.rfind(letters)
    if right < 0 or not _aligned(stem.stem, form_key[: right + len(letters)]):
        return FormAnalysis(OPPORTUNITY_NONE, None, "", "", False)
    change = _change(lemma_key, form_key, letters)
    return FormAnalysis(opportunity, False, change, "", True)


def analyze_form(
    lemma_key: str, stem: StemAnalysis, form_key: str
) -> FormAnalysis:
    """A plural form against its lemma (see the module docstring)."""
    opportunity = opportunity_of(lemma_key, form_key)
    if opportunity == OPPORTUNITY_NONE:
        return FormAnalysis(opportunity, None, "", "", form_key == lemma_key)
    if opportunity == OPPORTUNITY_A:
        aligned = _aligned(stem.stem, form_key[:-1]) if stem.stem else False
        return FormAnalysis(opportunity, None, "", "", aligned)
    if stem.nde_class:
        return _nde(stem, lemma_key, form_key, opportunity)
    if stem.dorsal:
        return _dorsal(stem, lemma_key, form_key, opportunity)
    if stem.target == CORONAL or stem.frontstem:
        return _plain(stem, lemma_key, form_key, opportunity)
    aligned = bool(stem.stem) and _aligned(stem.stem, form_key[:-1])
    return FormAnalysis(
        opportunity if aligned else OPPORTUNITY_NONE, None, "", "", aligned
    )


def gimpe_like(lemma_key: str, form_key: str) -> bool:
    """Invariable plural with a velar/palatal + front vowel (specie)."""
    return form_key == lemma_key and bool(_GIMPE_RE.search(lemma_key))


def root_internal_velar_front(stem: StemAnalysis, lemma_key: str) -> int:
    """Count of root-internal ‹chi che ghi ghe› (chilo, chiesa, ghiro):
    the non-derived-environment count. The stem-final NDE sequence of
    occhio/psiche is not counted."""
    body = lemma_key
    if stem.nde_class:
        body = lemma_key[: len(stem.stem) - 2]
    return len(_ROOT_VELAR_FRONT_RE.findall(body))
