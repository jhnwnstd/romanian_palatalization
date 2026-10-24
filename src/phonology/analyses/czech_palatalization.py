"""Czech velar and r palatalization: a provisional Logical Phonology analysis.

This module extends the analysis of ``latex.tex`` (Gorman & Winstead's
Romanian palatalization, :mod:`phonology.analyses.romanian_palatalization`)
to the Czech nominal inflection counted by :mod:`lexpipe.langs.czech`
(kluk/kluci, matka/matce, kniha/knize, Čech/Češi, doktor/doktoři;
biolog/biologové, park/v parku). It is the analysis as data: the
inventory patches, the underspecified targets, the suffix-borne
palatalizer, the ordered rules and the specification policies of four
variants, plus the hooks that feed the generated Czech cell data to the
shared validator (:mod:`phonology.validation`, driven by
``scripts/validate_lang_rules.py --lang cs``). The choices are argued,
and the results discussed, in ``docs/czech/phonology.md``.

The analysis
------------
**Inventory.** ``phonology/inventories/czech_features.json``: the 27
features of the Romanian inventory, the shared segments copied exactly,
and the Czech segments chosen so that the small (second) palatalization
is feature-filling: ť ď ň are t d n with [-ant, +dist] (no dorsal
component), ř is r with that same place plus [-son, -approx] (an
obstruent, JIPA 2012), x is the continuant k and ɦ the voiced x (a
dorsal, so h > z fills place exactly as g > z does).

**The trigger is the suffix, not the vowel** (Cavirani & Vanden
Wyngaerd 2026; Ketner 2005). -y/-ý [ɪ iː], -é and the -e of INS.SG -em,
LOC.PL -ech and the coronal vocative do not palatalize, while the
-e/-ě of LOC/DAT.SG, NOM.PL -i, LOC.PL -ích and the adjectival NOM.PL
-í do. So no rule here mentions a front vowel: the UR builder inserts a
floating palatalizer, a placeholder segment, between the stem and every
P ending (and only there), and the rules look for it:

    PAL1  [-syll, -cons, +cor, -ant, +dist, -strid]   CVW's PAL1 'I':
          case/number endings and the adjectival -í (small pattern)
    PAL3  [-syll, -cons, +cor, -ant, +dist, +strid]   CVW's PAL3 'H.I':
          the velar vocative -e (bože, člověče; big pattern)

No real segment is [-cons, +cor] (j and u̯ are CORONAL -), so the
placeholders are found by natural class and deleted after the
palatalization rules.

**Underspecification** (:data:`UNDERSPEC`). Velars are placeless (CVW
2026: the palatalizer gives a placeless consonant its place): /K G X Ɦ/
are k ɡ x ɦ without CORONAL, DORSAL and Strident (and, for /K G/,
DelRel; /G/ also Continuant, the Romanian /D/ device). /T D N/ are t d
n without Anterior and Distributed, /R/ is r without Anterior,
Distributed, Sonorant and Approximant. /Ḱ/ is the k of the
adjectival suffix -ský/-cký: placeless like /K/ but prespecified
[-strident], so it can take coronal place but not stridency (čeští,
němečtí: šť čť rather than *šc *čc); /S TS/ are the s and c of that
suffix, open for Anterior and Distributed. Prespecified k ɡ x ɦ t d c ɟ
are [-strident] (:data:`PATCHES`), so a stem whose target is
prespecified is inalterable by unification failure, as in Romanian.

**Rules** (:data:`RULES`, in order; every palatalization rule has a
broad, next-segment terminator, so only the consonant immediately
before the palatalizer is affected):

  1. ``big-pal`` [-son, +cons] ⊔ {+cor, -dor, +strid, -ant, +dist} /
     _ PAL3: k > t͡ʃ, g ɦ > ʒ, x > ʃ (bože, hříšníče).
  2. ``x-pal`` [-son, +cons, +cont, -voice] ⊔ {+cor, -dor, +strid,
     -ant, +dist} / _ PAL: ch > š before any palatalizer (Češi, mouše).
  3. ``velar-pal`` [-son, +cons] ⊔ {+cor, -dor, +strid, +ant, -dist}
     / _ PAL1: k > t͡s, g ɦ > z (kluci, lize, knize); it fails on /X/
     (already -ant), on /Ḱ/ (-strid) and on every prespecified segment.
  4. ``coronal-pal`` [+cons] ⊔ {+cor, -dor, -ant, +dist} / _ PAL:
     t d n r > ť ď ň ř and /Ḱ/ > ť (the palatalizer's own place).
  5. ``s-pal`` [+cons, +cor, +strid] ⊔ {-ant, +dist} / _ [+cons, +cor,
     -ant, +dist]: /S TS/ before the derived ť (čeští, němečtí), the
     Romanian S-palatalization fed by the dorsal rule (muscă).
  6. ``pal-deletion`` the placeholder is deleted.
  7. ``dorsal-default`` [-son, +cons] ⊔ {-cor, +dor, -strid};
     ``dorsal-stop-default`` [-son, +dor] ⊔ {-cont, -delrel}:
     unpalatalized /K G X Ɦ Ḱ/ surface as k ɡ x ɦ k.
  8. ``sibilant-continuant`` [+cor, +strid] ⊔ {+cont}: palatalized /G/
     (open for Continuant) is z, not *d͡z; ``affricate-release``
     [+cor, +strid, -cont] ⊔ {+delrel}: palatalized /K/ is t͡s.
  9. ``coronal-place-default`` [+cor] ⊔ {+ant, -dist}: unpalatalized
     /T D N R S TS/.
 10. ``trill-obstruent`` [+trill, -ant] ⊔ {-son, -approx} and
     ``trill-sonorant`` [+trill] ⊔ {+son, +approx}: ř vs r.

Specification policies (the variants)
-------------------------------------
Spell-out is not modelled: the UR is built from the attested form, with
its ending and its own palatalizer, so the palatalization rules are
tested given the morphology. The variants differ only in which
stem-final targets are underspecified:

  - ``lexical`` (:func:`build_profile`; Gorman & Winstead, the Inkelas
    learner): the target is underspecified iff its lexeme shows
    alternation evidence in E0 (nominal inflection) at τ = dictionary,
    the data layer's ``p_cells_E0_dictionary``. This is circular: a
    P form with the change is itself the evidence, so the boolean match
    tests the rule mechanics (trigger, locality, cluster outcomes,
    projection), not a lexical prediction.
  - ``all_underspecified``: every target is underspecified: how far a
    single set of rules gets (the rule TP (a) as a derivation).
  - ``lexical_loo``: the evidence must come from a cell outside the
    tested cell and its copies (NOM.PL = VOC.PL, feminine DAT.SG =
    LOC.SG): the non-circular, leave-one-cell-out version of
    ``lexical`` (the spell-out TP (c) as a derivation).
  - ``spellout``: Gorman & Winstead's spell-out claim as a policy: the
    non-palatalizing allomorphs (-ové/-é, -u, -ách) select prespecified
    roots, so a target is prespecified iff its lexeme takes one of them
    in some E0 cell at τ = dictionary (the data layer's
    ``n_cells_E0_dictionary``), and underspecified otherwise. It reads
    the N choices, never the P outcome, so it is not circular; the
    mixed stems (park: v parku but parcích) are its counterexamples.

The suffix type of an sk/ck cluster (/SḰ/ in -ský/-cký, s + /K/
elsewhere) is morphological and the same in every variant.

Rows and the attested form
--------------------------
:func:`iter_rows` reads ``data/czech/czech_cells.csv`` (one row per
lexeme × cell × form) and keeps the primary, aligned, not rejected rows
of velar, r and coronal_inferred stems with a P or N ending of the hard
paradigm, joined with the lemma IPA and the E0 evidence of
``czech_lexicon_complete.csv``. The data carry no transcription of
inflected forms, so the attested IPA (``ipa_pl``) is
:func:`transcribe` of the spelling: a re-implementation of the cs-IPA
conventions of the lemma transcriptions (final devoicing, voicing
assimilation, ř devoicing, ŋ, syllabic r l, ť ď ň before ě i í, hiatus
j), with the hard di ti ni of a loan taken from its lemma IPA. The UR is
the same transcription of the form with the stem-final change undone
and the palatalization across the stem edge suppressed, so the IPA match
tests the derived stem edge against the spelling; t/d/n rows are
``coronal_inferred`` (their change is invisible in spelling). The
comparison ignores stress, syllable marks, the allophones r̝̊ ŋ ɱ and
the degemination of t + affricate (:data:`CZECH_NORMALISERS`).

Hooks implemented here because the shared ones cannot carry Czech
(reported in docs/czech/phonology.md): :class:`CzechRow` subclasses
:class:`~phonology.lexicon.LexRow`; :func:`categorise` maps onto the
shared :class:`~phonology.diagnostics.categorise.RefinedCategory`,
which has no Czech labels, and :func:`mismatch_class` names the Czech
mismatch types for :func:`render_summary`; the normaliser cascade does
its folding in ``strip``, since :class:`Normalisation` has no Czech
labels either.
"""

from __future__ import annotations

import argparse
import csv
import sys
import unicodedata
from collections import Counter, defaultdict
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum
from functools import lru_cache
from pathlib import Path
from types import MappingProxyType
from typing import Final, cast

from lexpipe import PROJECT_ROOT
from lexpipe.csvio import split_list
from lexpipe.spec import LANGUAGES, LanguageInfo

from ..diagnostics.categorise import Categorisation, RefinedCategory
from ..diagnostics.compare import NormaliserCascade, compare_ipa
from ..inventory import (
    FeatureInventory,
    FeaturePatch,
    UnderspecifiedSegment,
    inventory_path,
    load_inventory_file,
)
from ..lexicon import ClusterTag, LexRow
from ..lp import NaturalClass, natural_class
from ..profile import AnalysisProfile
from ..rules import Rule, SegmentDeletionRule, UnificationRule
from ..search import Direction, Search
from ..segments import (
    NON_SEGMENTAL_MARKS,
    LengthHandling,
    Segment,
    TokenizerConfig,
    Word,
    tokenize_ipa,
)
from ..validation import (
    RowResult,
    TraceResult,
    ValidationReport,
    evaluate_row,
    render_trace,
    run_validation,
)

LANGUAGE: Final[LanguageInfo] = LANGUAGES["cs"]
_INVENTORY_NAME: Final[str] = "czech"

# Target classes of the data layer (lexpipe.langs.czech.ortho).
VELAR: Final[str] = "velar"
R_CLASS: Final[str] = "r"
CORONAL: Final[str] = "coronal_inferred"
TARGET_CLASSES: Final[tuple[str, ...]] = (VELAR, R_CLASS, CORONAL)
P: Final[str] = "P"
N: Final[str] = "N"

# Stem finals the analysis covers (ch is one letter).
DORSAL_STEMS: Final[frozenset[str]] = frozenset({"k", "g", "h", "ch"})
CORONAL_STEMS: Final[frozenset[str]] = frozenset({"t", "d", "n"})
TARGET_STEMS: Final[frozenset[str]] = (
    DORSAL_STEMS | CORONAL_STEMS | frozenset({"r"})
)


def _frozen(d: Mapping[str, str]) -> Mapping[str, str]:
    """A read-only copy of a feature pattern (shared by several rules)."""
    return MappingProxyType(dict(d))


# ---------------------------------------------------------------------------
# Transcription (cs-IPA conventions)
# ---------------------------------------------------------------------------

STRESS: Final[str] = "ˈ"
LENGTH: Final[str] = "ː"

_VOWEL_LETTERS: Final[Mapping[str, str]] = MappingProxyType(
    {
        "a": "a",
        "á": "aː",
        "e": "ɛ",
        "é": "ɛː",
        "ě": "ɛ",
        "i": "ɪ",
        "í": "iː",
        "y": "ɪ",
        "ý": "iː",
        "o": "o",
        "ó": "oː",
        "u": "u",
        "ú": "uː",
        "ů": "uː",
        "ä": "ɛ",
        "ë": "ɛ",
        "ö": "ɛː",
        "ü": "ɪ",
    }
)
_CONSONANT_LETTERS: Final[Mapping[str, str]] = MappingProxyType(
    {
        "b": "b",
        "c": "t͡s",
        "č": "t͡ʃ",
        "d": "d",
        "ď": "ɟ",
        "f": "f",
        "g": "ɡ",
        "h": "ɦ",
        "j": "j",
        "k": "k",
        "l": "l",
        "m": "m",
        "n": "n",
        "ň": "ɲ",
        "p": "p",
        "q": "k",
        "r": "r",
        "ř": "r̝",
        "s": "s",
        "š": "ʃ",
        "t": "t",
        "ť": "c",
        "v": "v",
        "w": "v",
        "z": "z",
        "ž": "ʒ",
        # Spelling units of more than one letter (see _units).
        "ch": "x",
        "dž": "d͡ʒ",
        "tz": "t͡s",
    }
)
_DIGRAPHS: Final[tuple[str, ...]] = ("ch", "dž", "tz")
# t d n read as ť ď ň before ě i í (soft reading).
_SOFTENED: Final[Mapping[str, str]] = MappingProxyType(
    {"t": "c", "d": "ɟ", "n": "ɲ"}
)
_SOFTENING_VOWELS: Final[frozenset[str]] = frozenset({"ě", "i", "í"})
# bě pě vě fě [bjɛ …], mě [mɲɛ].
_GLIDE_BEFORE_E: Final[frozenset[str]] = frozenset({"b", "p", "v", "f", "w"})
_HIATUS_VOWELS: Final[frozenset[str]] = frozenset({"i", "í", "y", "ý"})
_DIPHTHONG_FIRST: Final[frozenset[str]] = frozenset({"o", "a", "e"})

_VOICING_PAIRS: Final[tuple[tuple[str, str], ...]] = (
    ("p", "b"),
    ("t", "d"),
    ("c", "ɟ"),
    ("k", "ɡ"),
    ("f", "v"),
    ("s", "z"),
    ("ʃ", "ʒ"),
    ("x", "ɦ"),
    ("t͡s", "d͡z"),
    ("t͡ʃ", "d͡ʒ"),
    ("r̝̊", "r̝"),
)
_TO_VOICED: Final[Mapping[str, str]] = MappingProxyType(dict(_VOICING_PAIRS))
_TO_VOICELESS: Final[Mapping[str, str]] = MappingProxyType(
    {voiced: voiceless for voiceless, voiced in _VOICING_PAIRS}
)
_OBSTRUENT_SEGMENTS: Final[frozenset[str]] = frozenset(_TO_VOICED) | frozenset(
    _TO_VOICELESS
)
_VOICED_OBSTRUENTS: Final[frozenset[str]] = frozenset(_TO_VOICELESS)
_VOICELESS_OBSTRUENTS: Final[frozenset[str]] = frozenset(_TO_VOICED)
_VOWEL_SEGMENTS: Final[frozenset[str]] = frozenset(_VOWEL_LETTERS.values())
_SONORANT_SEGMENTS: Final[frozenset[str]] = frozenset(
    {"m", "n", "ɲ", "ŋ", "l", "r", "j"}
)
_SYLLABIC: Final[str] = "\N{COMBINING VERTICAL LINE BELOW}"


@dataclass(frozen=True, slots=True)
class Transcription:
    """A spelling's segments with the spelling index each one came from.

    ``origins[i]`` is the start index (in the NFC, lowercased text) of
    the spelling unit segment ``i`` was read from, so a caller that
    knows where the stem ends in the spelling can split the segments
    there (:meth:`split`). Segments are cs-IPA tokens, allophones
    included (r̝̊ ŋ); :func:`phonemic` folds those away.
    """

    segments: tuple[str, ...]
    origins: tuple[int, ...]

    def split(self, boundary: int) -> tuple[tuple[str, ...], tuple[str, ...]]:
        """``(stem segments, ending segments)`` at spelling index
        ``boundary``."""
        stem = tuple(
            s for s, o in zip(self.segments, self.origins) if o < boundary
        )
        rest = tuple(
            s for s, o in zip(self.segments, self.origins) if o >= boundary
        )
        return stem, rest

    def render(self) -> str:
        """The cs-IPA string: initial stress mark plus the segments."""
        return render(self.segments)


def render(segments: Iterable[str]) -> str:
    """cs-IPA string of ``segments``, with the (fixed, initial) stress."""
    return STRESS + "".join(segments)


def _units(text: str) -> list[tuple[str, int]]:
    """Spelling units with their start index (ch, dž, tz are one unit)."""
    out: list[tuple[str, int]] = []
    i = 0
    while i < len(text):
        for digraph in _DIGRAPHS:
            if text.startswith(digraph, i):
                out.append((digraph, i))
                i += len(digraph)
                break
        else:
            out.append((text[i], i))
            i += 1
    return out


@lru_cache(maxsize=1 << 16)
def transcribe(
    text: str,
    *,
    boundary: int | None = None,
    underlying: bool = False,
    hard: frozenset[str] = frozenset(),
) -> Transcription:
    """cs-IPA transcription of a Czech spelling (one word).

    The conventions are those of en.wiktionary's ``{{cs-IPA}}``, which
    transcribes the lemmas of the data: i y [ɪ], í ý [iː], e ě [ɛ],
    ou [ou̯]; ť ď ň before ě i í; bě pě vě fě [Cjɛ], mě [mɲɛ]; a hiatus
    j after i before a vowel (biologie [bɪjoloɡɪjɛ]); then
    :func:`_postlexical` (voicing, degemination, ŋ, syllabic r l, ř
    devoicing).

    ``boundary`` is the spelling index where the ending starts.
    ``hard`` lists the letters among t d n whose di ti ni (dí tí ní)
    inside the stem are hard, as in loans (diplomat [dɪplomat]); the
    stem-final consonant before the ending is not affected.
    ``underlying=True`` suppresses what the ending does to the stem
    final (ť ď ň before ě i í, the labial glide before ě): the
    palatalization at the stem edge is the rules' business, so the UR
    builder transcribes the form this way.
    """
    lowered = unicodedata.normalize("NFC", text).lower()
    units = _units(lowered)
    segs: list[str] = []
    origins: list[int] = []

    def emit(segment: str, origin: int) -> None:
        segs.append(segment)
        origins.append(origin)

    k = 0
    n = len(units)
    while k < n:
        unit, pos = units[k]
        nxt, nxt_pos = units[k + 1] if k + 1 < n else ("", len(lowered))
        at_edge = boundary is not None and pos < boundary <= nxt_pos
        if unit in _VOWEL_LETTERS:
            if unit in _DIPHTHONG_FIRST and nxt == "u" and not at_edge:
                emit(_VOWEL_LETTERS[unit], pos)
                emit("u̯", nxt_pos)
                k += 2
                continue
            emit(_VOWEL_LETTERS[unit], pos)
            if unit in _HIATUS_VOWELS and nxt in _VOWEL_LETTERS:
                if not at_edge:
                    emit("j", pos)
        elif unit == "x":
            # ks; the prefix ex- before a vowel is [ɛɡz] (exil).
            if k == 1 and units[0][0] == "e" and nxt in _VOWEL_LETTERS:
                emit("ɡ", pos)
                emit("z", pos)
            else:
                emit("k", pos)
                emit("s", pos)
        elif unit in _CONSONANT_LETTERS:
            segment = _CONSONANT_LETTERS[unit]
            suppressed = at_edge and underlying
            if nxt == "ě" and not suppressed:
                if unit in _SOFTENED:
                    segment = _SOFTENED[unit]
                elif unit in _GLIDE_BEFORE_E:
                    emit(segment, pos)
                    segment = "j"
                elif unit == "m":
                    emit(segment, pos)
                    segment = "ɲ"
            elif (
                nxt in _SOFTENING_VOWELS
                and unit in _SOFTENED
                and not suppressed
                and not (unit in hard and not at_edge)
            ):
                segment = _SOFTENED[unit]
            emit(segment, pos)
        else:
            emit(unit, pos)  # unknown letter: left for the caller
        k += 1
    protected: frozenset[int] = frozenset()
    if underlying and boundary is not None:
        # The stem final and the consonant before it may change in the
        # derivation (měkk-ý: k + /K/; abchaz-sk-ý: z + /S/), so they
        # are not merged with their neighbours here: degemination is
        # surface phonology, done by the comparison (strip_czech).
        stem_units = [pos for _, pos in units if pos < boundary]
        protected = frozenset(stem_units[-2:])
    final_segs, final_origins = _postlexical(segs, origins, protected)
    return Transcription(tuple(final_segs), tuple(final_origins))


def _postlexical(
    segs: list[str],
    origins: list[int],
    protected: frozenset[int] = frozenset(),
) -> tuple[list[str], list[int]]:
    """Voicing, degemination, ŋ, syllabic r l and ř devoicing.

    A segment whose origin is in ``protected`` is never merged into the
    segment before it (degemination).
    """
    out = list(segs)
    # -ism- (absolutismus [ɪzmus]).
    for i in range(1, len(out) - 1):
        if out[i] == "s" and out[i - 1] == "ɪ" and out[i + 1] == "m":
            out[i] = "z"
    # Final devoicing and regressive voicing assimilation, right to
    # left. v undergoes but never triggers; ř triggers only when it is
    # itself voiceless (debř [dɛpr̝̊]).
    voice_next: bool | None = False
    for i in range(len(out) - 1, -1, -1):
        original = out[i]
        if original not in _OBSTRUENT_SEGMENTS:
            voice_next = None
            continue
        if voice_next is True:
            out[i] = _TO_VOICED.get(original, original)
        elif voice_next is False:
            out[i] = _TO_VOICELESS.get(original, original)
        if original == "v":
            voice_next = None
        elif original == "r̝":
            voice_next = False if out[i] == "r̝̊" else None
        else:
            voice_next = out[i] in _VOICED_OBSTRUENTS
    # Degemination: identical consonants (not affricates: hrdliččí
    # [t͡ʃt͡ʃ]), t + affricate (dcera [t͡sɛra]), n + ň (anglofonní), and an
    # affricate before its own fricative (apačština [t͡ʃc]).
    merged: list[str] = []
    merged_origins: list[int] = []
    for segment, origin in zip(out, origins):
        if merged and origin not in protected:
            last = merged[-1]
            if (
                last == segment
                and segment not in _VOWEL_SEGMENTS
                and "͡" not in segment
            ):
                continue
            if (last == "t" and segment in ("t͡s", "t͡ʃ")) or (
                last == "n" and segment == "ɲ"
            ):
                merged[-1] = segment
                merged_origins[-1] = origin
                continue
            if (last, segment) in (("t͡ʃ", "ʃ"), ("t͡s", "s")):
                continue
        merged.append(segment)
        merged_origins.append(origin)
    out = merged
    # n before a velar is [ŋ].
    for i in range(len(out) - 1):
        if out[i] == "n" and out[i + 1] in ("k", "ɡ"):
            out[i] = "ŋ"
    # Syllabic r l: after a consonant and before a consonant, or
    # word-finally after an obstruent (vlk, bratr; but hamr [ɦamr]).
    nonsyllabic = _VOWEL_SEGMENTS | {"j", "u̯"}
    for i, segment in enumerate(out):
        if segment not in ("r", "l"):
            continue
        left = out[i - 1] if i > 0 else ""
        right = out[i + 1] if i + 1 < len(out) else ""
        if not left or left in nonsyllabic:
            continue
        if (right and right not in nonsyllabic) or (
            not right and left not in _SONORANT_SEGMENTS
        ):
            out[i] = segment + _SYLLABIC
    # ř after a voiceless obstruent is voiceless (tři, sestře).
    for i in range(1, len(out)):
        if out[i] == "r̝" and out[i - 1] in _VOICELESS_OBSTRUENTS:
            out[i] = "r̝̊"
    return out, merged_origins


# ---------------------------------------------------------------------------
# Tokenizer, phonemic folding and the comparison convention
# ---------------------------------------------------------------------------

# Allophones the inventory does not list, folded onto their phonemes:
# voiceless ř (tři), n and m before velars and labiodentals, and the
# voiced x and c of voicing assimilation (cs-IPA writes the voiced ch of
# dechberoucí as ɦ or ɣ; salcburský [sald͡zb…]).
_ALLOPHONES: Final[Mapping[str, str]] = MappingProxyType(
    {"r̝̊": "r̝", "ŋ": "n", "ɱ": "m", "ɣ": "ɦ", "d͡z": "t͡s"}
)
# Other spellings of inventory segments in Wiktionary IPA: the
# cs.wiktionary convention (ɔ ʊ ɔʊ̯ e), loans (g h w), ASCII.
_VARIANT_SYMBOLS: Final[Mapping[str, str]] = MappingProxyType(
    {
        "ɔ": "o",
        "ɔː": "oː",
        "ᴐ": "o",
        "ʊ": "u",
        "ʊː": "uː",
        "ʊ̯": "u̯",
        "e": "ɛ",
        "eː": "ɛː",
        "ε": "ɛ",
        "i": "ɪ",
        "g": "ɡ",
        "h": "ɦ",
        "w": "u̯",
    }
)

CZECH_TOKENIZER: Final[TokenizerConfig] = TokenizerConfig(
    skip=NON_SEGMENTAL_MARKS | frozenset("ʔˀ‿()[]/"),
    fold={**_ALLOPHONES, **_VARIANT_SYMBOLS},
    keep_combining=True,
    length=LengthHandling.ATTACH,
)
"""Czech convention: diacritics stay on their base (r̝ r̩ u̯ are one
token), ː attaches to its segment (length is contrastive), allophones
and variant symbols fold onto inventory labels, the glottal stop is
skipped."""


def phonemic(segment: str) -> str:
    """The inventory label of a transcription segment (allophones and
    variant symbols folded, as :data:`CZECH_TOKENIZER` does)."""
    return CZECH_TOKENIZER.fold.get(segment, segment)


def czech_tokens(ipa: str) -> tuple[str, ...]:
    """Tokens of an IPA string under :data:`CZECH_TOKENIZER`."""
    return tokenize_ipa(ipa, tokenizer=CZECH_TOKENIZER)


_ABSORBED_BEFORE: Final[Mapping[str, frozenset[str]]] = MappingProxyType(
    {
        # t (d) before an affricate, n before ň: the affricate (ň) wins.
        "t͡s": frozenset({"t", "d"}),
        "t͡ʃ": frozenset({"t", "d"}),
        "d͡ʒ": frozenset({"t", "d"}),
        "ɲ": frozenset({"n"}),
    }
)
_ABSORBED_AFTER: Final[Mapping[str, str]] = MappingProxyType(
    # an affricate swallows its own fricative (apačština [t͡ʃc])
    {"ʃ": "t͡ʃ", "s": "t͡s", "ʒ": "d͡ʒ"}
)


def _degeminate(tokens: Sequence[str]) -> list[str]:
    """Merge the geminates the surface transcription merges: identical
    consonants (affricates included), t/d + affricate, n + ň, an
    affricate + its own fricative; repeated until nothing merges."""
    out: list[str] = []
    for token in tokens:
        out.append(token)
        while len(out) > 1:
            last, before = out[-1], out[-2]
            if last == before and last not in _VOWEL_SEGMENTS:
                out.pop()
            elif before in _ABSORBED_BEFORE.get(last, frozenset()):
                del out[-2]
            elif _ABSORBED_AFTER.get(last) == before:
                out.pop()
            else:
                break
    return out


def strip_czech(ipa: str) -> str:
    """The comparison form of a Czech transcription.

    Drops stress, syllable dots, the glottal stop and brackets, folds
    allophones (r̝̊ ŋ ɱ ɣ d͡z) and variant symbols onto inventory labels,
    and degeminates (cs-IPA writes matce as [mat͡sɛ]; the derivation
    gives /matt͡sɛ/). Length is kept: it is contrastive in Czech.
    """
    return "".join(_degeminate(czech_tokens(ipa)))


CZECH_NORMALISERS: Final[NormaliserCascade] = NormaliserCascade(
    strip=strip_czech,
    steps=(),
)
"""Every Czech difference the comparison ignores is in ``strip``; the
shared :class:`~phonology.diagnostics.compare.Normalisation`
vocabulary has no Czech labels, so there are no labelled steps."""


# ---------------------------------------------------------------------------
# Patches and underspecified segments
# ---------------------------------------------------------------------------

_NON_STRIDENT: Final[Mapping[str, str]] = _frozen({"Strident": "-"})

PATCHES: Final[tuple[FeaturePatch, ...]] = tuple(
    FeaturePatch(
        segment=label,
        updates=_NON_STRIDENT,
        reason=reason,
    )
    for label, reason in (
        (
            "k",
            "Velars are [-strident]. Prespecified on k ɡ x ɦ so that a "
            "stem whose velar is not underspecified is inalterable to "
            "velar-pal (assibilation), and so that the -ský/-cký /Ḱ/ "
            "can keep the value while /K/ clears it.",
        ),
        ("ɡ", "See k: velars are prespecified [-strident]."),
        ("x", "See k: velars are prespecified [-strident]."),
        ("ɦ", "See k: velars are prespecified [-strident]."),
        (
            "t",
            "As in the Romanian analysis, prespecified /t d/ are "
            "[-strident], the value assibilation (here velar-pal) would "
            "contradict; /T D/ keep it, so velar-pal never applies to "
            "them and t > ť is the coronal place fill alone.",
        ),
        ("d", "See t."),
        (
            "c",
            "ť ď are non-strident palatal stops: a prespecified ť before "
            "a palatalizer must not turn into t͡ʃ.",
        ),
        ("ɟ", "See c."),
    )
)

_DORSAL_PLACE: Final[frozenset[str]] = frozenset({"CORONAL", "DORSAL"})
_CORONAL_PLACE: Final[frozenset[str]] = frozenset({"Anterior", "Distributed"})

UNDERSPEC: Final[tuple[UnderspecifiedSegment, ...]] = (
    UnderspecifiedSegment(
        label="K",
        base="k",
        clear=_DORSAL_PLACE | {"Strident", "DelRel"},
        reason=(
            "Alternating k: placeless (CVW 2026: velars have no place "
            "element; the palatalizer supplies one) and open for "
            "Strident and DelRel, which velar-pal and affricate-release "
            "fill (t͡s), or the dorsal defaults (k)."
        ),
    ),
    UnderspecifiedSegment(
        label="Ḱ",
        base="k",
        clear=_DORSAL_PLACE,
        reason=(
            "The k of the adjectival suffix -ský/-cký: placeless but "
            "prespecified [-strident] (patch), so velar-pal fails on it "
            "and coronal-pal gives it the palatalizer's own place: ť, "
            "as in čeští, němečtí."
        ),
    ),
    UnderspecifiedSegment(
        label="G",
        base="ɡ",
        clear=_DORSAL_PLACE | {"Strident", "DelRel", "Continuant"},
        reason=(
            "Alternating g: as /K/, and open for Continuant (the "
            "Romanian /D/ device), so the palatalized stop becomes the "
            "fricative z (lize), not *d͡z."
        ),
    ),
    UnderspecifiedSegment(
        label="X",
        base="x",
        clear=_DORSAL_PLACE | {"Strident"},
        reason="Alternating ch: placeless, open for Strident (š).",
    ),
    UnderspecifiedSegment(
        label="Ɦ",
        base="ɦ",
        clear=_DORSAL_PLACE | {"Strident"},
        reason=(
            "Alternating h: placeless, open for Strident; voiced and "
            "continuant, so it becomes z like g (knize, vrazi)."
        ),
    ),
    UnderspecifiedSegment(
        label="R",
        base="r",
        clear=_CORONAL_PLACE | {"Sonorant", "Approximant"},
        reason=(
            "Alternating r: open for coronal place, so the palatalizer's "
            "[-ant, +dist] makes it ř, and for Sonorant and Approximant, "
            "which the trill defaults fill (ř is an obstruent)."
        ),
    ),
    UnderspecifiedSegment(
        label="T",
        base="t",
        clear=_CORONAL_PLACE,
        reason="Alternating t (ť): open for coronal place.",
    ),
    UnderspecifiedSegment(
        label="D",
        base="d",
        clear=_CORONAL_PLACE,
        reason="Alternating d (ď): open for coronal place.",
    ),
    UnderspecifiedSegment(
        label="N",
        base="n",
        clear=_CORONAL_PLACE,
        reason="Alternating n (ň): open for coronal place.",
    ),
    UnderspecifiedSegment(
        label="S",
        base="s",
        clear=_CORONAL_PLACE,
        reason=(
            "The s of the suffix -ský: open for coronal place, so s-pal "
            "makes it š before the derived ť (čeští)."
        ),
    ),
    UnderspecifiedSegment(
        label="TS",
        base="t͡s",
        clear=_CORONAL_PLACE,
        reason=(
            "The c of the suffix -cký: open for coronal place, so s-pal "
            "makes it č before the derived ť (němečtí)."
        ),
    ),
)

# Surface segment → its underspecified twin, for the stem-final target.
UNDERSPEC_FOR: Final[Mapping[str, str]] = MappingProxyType(
    {
        "k": "K",
        "ɡ": "G",
        "x": "X",
        "ɦ": "Ɦ",
        "r": "R",
        "t": "T",
        "d": "D",
        "n": "N",
    }
)
SUFFIX_VELAR: Final[str] = "Ḱ"
SUFFIX_SIBILANT_FOR: Final[Mapping[str, str]] = MappingProxyType(
    {"s": "S", "t͡s": "TS"}
)
# Stem final (spelling) → the segment the transcription gives it.
BASE_SEGMENT: Final[Mapping[str, str]] = MappingProxyType(
    {
        "k": "k",
        "g": "ɡ",
        "h": "ɦ",
        "ch": "x",
        "r": "r",
        "t": "t",
        "d": "d",
        "n": "n",
    }
)


def load_inventory() -> FeatureInventory:
    """The Czech inventory with this analysis's patches applied."""
    return load_inventory_file(
        _INVENTORY_NAME,
        patches=PATCHES,
        underspec=UNDERSPEC,
    )


# ---------------------------------------------------------------------------
# The palatalizers (floating feature bundles)
# ---------------------------------------------------------------------------

PAL1: Final[str] = "PAL1"
PAL3: Final[str] = "PAL3"
_PALATALIZER_PLACE: Final[Mapping[str, str]] = _frozen(
    {
        "Syllabic": "-",
        "Consonantal": "-",
        "CORONAL": "+",
        "Anterior": "-",
        "Distributed": "+",
    }
)
PALATALIZER_FEATURES: Final[Mapping[str, Mapping[str, str]]] = (
    MappingProxyType(
        {
            PAL1: _frozen({**_PALATALIZER_PLACE, "Strident": "-"}),
            PAL3: _frozen({**_PALATALIZER_PLACE, "Strident": "+"}),
        }
    )
)
"""PAL1 is CVW's 'I' (a coronal palatal place); PAL3 'H.I' adds the
stridency that makes the big pattern postalveolar (k > t͡ʃ)."""


def palatalizer(label: str) -> Segment:
    """The placeholder segment of a palatalizer (:data:`PAL1`,
    :data:`PAL3`)."""
    return Segment(label=label, features=PALATALIZER_FEATURES[label])


# ---------------------------------------------------------------------------
# Natural classes (LP [ ])
# ---------------------------------------------------------------------------

CONSONANT: Final[NaturalClass] = natural_class(Consonantal="+")
OBSTRUENT: Final[NaturalClass] = natural_class(Consonantal="+", Sonorant="-")
VOICELESS_FRICATIVE: Final[NaturalClass] = natural_class(
    Consonantal="+", Sonorant="-", Continuant="+", Voice="-"
)
DORSAL_OBSTRUENT: Final[NaturalClass] = natural_class(Sonorant="-", DORSAL="+")
SIBILANT: Final[NaturalClass] = natural_class(
    Consonantal="+", CORONAL="+", Strident="+"
)
STRIDENT_STOP: Final[NaturalClass] = natural_class(
    CORONAL="+", Strident="+", Continuant="-"
)
CORONAL_SEGMENT: Final[NaturalClass] = natural_class(CORONAL="+")
POSTALVEOLAR_CONSONANT: Final[NaturalClass] = natural_class(
    Consonantal="+", CORONAL="+", Anterior="-", Distributed="+"
)
TRILL: Final[NaturalClass] = natural_class(Trill="+")
PALATAL_TRILL: Final[NaturalClass] = natural_class(Trill="+", Anterior="-")
ANY_PALATALIZER: Final[NaturalClass] = natural_class(
    Syllabic="-", Consonantal="-", CORONAL="+"
)
SMALL_PALATALIZER: Final[NaturalClass] = natural_class(
    Consonantal="-", CORONAL="+", Strident="-"
)
BIG_PALATALIZER: Final[NaturalClass] = natural_class(
    Consonantal="-", CORONAL="+", Strident="+"
)

# Feature changes (LP { }).
POSTALVEOLAR_SIBILANT: Final[Mapping[str, str]] = _frozen(
    {
        "CORONAL": "+",
        "DORSAL": "-",
        "Strident": "+",
        "Anterior": "-",
        "Distributed": "+",
    }
)
ALVEOLAR_SIBILANT: Final[Mapping[str, str]] = _frozen(
    {
        "CORONAL": "+",
        "DORSAL": "-",
        "Strident": "+",
        "Anterior": "+",
        "Distributed": "-",
    }
)
PALATAL_PLACE: Final[Mapping[str, str]] = _frozen(
    {"CORONAL": "+", "DORSAL": "-", "Anterior": "-", "Distributed": "+"}
)


def _next_segment(condition: NaturalClass) -> Search:
    """A broad rightward SEARCH: the next segment, which must be in
    ``condition`` (latex.tex: the dorsal rule's terminator)."""
    return Search(
        direction=Direction.RIGHT,
        terminator=NaturalClass.universal().spec,
        condition=condition.spec,
    )


# ---------------------------------------------------------------------------
# Rules
# ---------------------------------------------------------------------------

BIG_PAL: Final[UnificationRule] = UnificationRule(
    name="big-pal",
    target=OBSTRUENT.spec,
    supply=POSTALVEOLAR_SIBILANT,
    search=_next_segment(BIG_PALATALIZER),
)
"""The big pattern before PAL3 (VOC -e of velars): k > t͡ʃ, g ɦ > ʒ, x > ʃ."""

X_PAL: Final[UnificationRule] = UnificationRule(
    name="x-pal",
    target=VOICELESS_FRICATIVE.spec,
    supply=POSTALVEOLAR_SIBILANT,
    search=_next_segment(ANY_PALATALIZER),
)
"""ch > š before any palatalizer. Ordered before velar-pal, which then
fails on the -ant /X/; voiced /Ɦ/ is not in the target and becomes z."""

VELAR_PAL: Final[UnificationRule] = UnificationRule(
    name="velar-pal",
    target=OBSTRUENT.spec,
    supply=ALVEOLAR_SIBILANT,
    search=_next_segment(SMALL_PALATALIZER),
)
"""Small velar palatalization (assibilation): k > t͡s, g ɦ > z.

The supply is coronal stridency with the unmarked alveolar place, CVW's
phonetic interpretation of a placeless stop plus 'I' as [ʦ]. Every
prespecified obstruent resists by unification failure: k ɡ x ɦ (-cor),
t d c ɟ (-strid, patched), the postalveolars and ř (-ant); s z t͡s
already have the bundle. So do /X/ (-ant after x-pal) and /Ḱ/ (-strid).
"""

CORONAL_PAL: Final[UnificationRule] = UnificationRule(
    name="coronal-pal",
    target=CONSONANT.spec,
    supply=PALATAL_PLACE,
    search=_next_segment(ANY_PALATALIZER),
)
"""The palatalizer's own place: t d n r > ť ď ň ř, and /Ḱ/ > ť.

Fails on everything with fixed contrary place: t d n r l s z t͡s
(+ant), labials and prespecified velars (-cor), and on /K G Ɦ/ already
made +ant by velar-pal."""

S_PAL: Final[UnificationRule] = UnificationRule(
    name="s-pal",
    target=SIBILANT.spec,
    supply=_frozen({"Anterior": "-", "Distributed": "+"}),
    search=_next_segment(POSTALVEOLAR_CONSONANT),
)
"""/S TS/ > š č before a postalveolar consonant (the derived ť of /Ḱ/):
-ský > -ští, -cký > -čtí. Latex.tex's S-palatalization fed by the
dorsal rule (muscă), with a strictly adjacent trigger."""

PAL_DELETION: Final[SegmentDeletionRule] = SegmentDeletionRule(
    name="pal-deletion",
    target=ANY_PALATALIZER.spec,
)
"""The floating palatalizer never surfaces on its own here."""

DORSAL_DEFAULT: Final[UnificationRule] = UnificationRule(
    name="dorsal-default",
    target=OBSTRUENT.spec,
    supply=_frozen({"CORONAL": "-", "DORSAL": "+", "Strident": "-"}),
)
"""Unpalatalized /K Ḱ G X Ɦ/ are velars (latex.tex e:default)."""

DORSAL_STOP_DEFAULT: Final[UnificationRule] = UnificationRule(
    name="dorsal-stop-default",
    target=DORSAL_OBSTRUENT.spec,
    supply=_frozen({"Continuant": "-", "DelRel": "-"}),
)
"""Unpalatalized /K G/ are stops (x ɦ, +cont, resist)."""

SIBILANT_CONTINUANT: Final[UnificationRule] = UnificationRule(
    name="sibilant-continuant",
    target=SIBILANT.spec,
    supply=_frozen({"Continuant": "+"}),
)
"""Palatalized /G/ (open for Continuant) is z or ʒ; /K/ (-cont) stays an
affricate. As in Romanian, the language lacks the voiced affricate the
rules would otherwise produce."""

AFFRICATE_RELEASE: Final[UnificationRule] = UnificationRule(
    name="affricate-release",
    target=STRIDENT_STOP.spec,
    supply=_frozen({"DelRel": "+"}),
)
"""A strident stop is an affricate: palatalized /K/ is t͡s or t͡ʃ."""

CORONAL_PLACE_DEFAULT: Final[UnificationRule] = UnificationRule(
    name="coronal-place-default",
    target=CORONAL_SEGMENT.spec,
    supply=_frozen({"Anterior": "+", "Distributed": "-"}),
)
"""Unpalatalized /T D N R S TS/ are alveolar."""

TRILL_OBSTRUENT: Final[UnificationRule] = UnificationRule(
    name="trill-obstruent",
    target=PALATAL_TRILL.spec,
    supply=_frozen({"Sonorant": "-", "Approximant": "-"}),
)
"""The palatalized trill ř is an obstruent."""

TRILL_SONORANT: Final[UnificationRule] = UnificationRule(
    name="trill-sonorant",
    target=TRILL.spec,
    supply=_frozen({"Sonorant": "+", "Approximant": "+"}),
)
"""Unpalatalized /R/ is the sonorant r."""

RULES: Final[tuple[Rule, ...]] = cast(
    "tuple[Rule, ...]",
    (
        BIG_PAL,
        X_PAL,
        VELAR_PAL,
        CORONAL_PAL,
        S_PAL,
        PAL_DELETION,
        DORSAL_DEFAULT,
        DORSAL_STOP_DEFAULT,
        SIBILANT_CONTINUANT,
        AFFRICATE_RELEASE,
        CORONAL_PLACE_DEFAULT,
        TRILL_OBSTRUENT,
        TRILL_SONORANT,
    ),
)
"""The ordered rule set (module docstring)."""

PALATALIZATION_RULES: Final[tuple[Rule, ...]] = cast(
    "tuple[Rule, ...]", (BIG_PAL, X_PAL, VELAR_PAL, CORONAL_PAL, S_PAL)
)
PALATALIZATION_RULE_NAMES: Final[frozenset[str]] = frozenset(
    r.name for r in PALATALIZATION_RULES
)
"""The rules whose firing means the stem alternated (the deletion and
the defaults apply whether or not it did)."""


# ---------------------------------------------------------------------------
# Rows
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class CzechRow(LexRow):
    """One validated cell form: the shared row plus the Czech fields.

    The shared fields are filled as the validator expects: ``plural``
    is the cell form, ``opportunity`` the cell and ending ("NOM.PL -i"),
    ``ipa_lemma`` the lemma's Wiktionary IPA, ``ipa_pl`` the attested
    (transcribed) form, ``exception_reason`` the data layer's
    ``exception_label`` and ``cluster`` the Romanian tag (always NONE;
    the Czech cluster is ``czech_cluster``).
    """

    lexeme_id: str
    cell: str
    ending: str
    allomorph: str
    n_kind: str
    target_class: str
    czech_cluster: str
    orth_change: str
    gender_class: str
    suffix_cluster: bool
    alt_cells: tuple[str, ...]
    resist_cells: tuple[str, ...]
    hard_coronals: frozenset[str]

    @property
    def stem_spelling(self) -> str:
        """The form without its ending (the possibly mutated stem)."""
        if not self.ending:
            return self.plural
        return self.plural[: len(self.plural) - len(self.ending)]


_POS: Final[Mapping[str, str]] = MappingProxyType({"noun": "N", "adj": "ADJ"})
_SUFFIX_ENDINGS: Final[tuple[str, ...]] = ("ský", "cký")


def suffix_cluster(lemma: str, cluster: str) -> bool:
    """Is the sk/ck cluster the adjectival suffix -ský/-cký (/SḰ TSḰ/)?

    Morphology, read off the lemma (adjectives and the adjectival nouns
    hospodský, krámský): the same in every variant.
    """
    return cluster in ("sk", "ck") and lemma.lower().endswith(_SUFFIX_ENDINGS)


def _bool(text: str | None) -> bool:
    return (text or "").strip().lower() == "true"


def _cells(text: str | None) -> tuple[str, ...]:
    return tuple(c for c in split_list((text or "").strip()) if c)


def row_factory(r: Mapping[str, str]) -> LexRow:
    """A :class:`CzechRow` from one ``czech_cells.csv`` row.

    Reads the cells schema (docs/czech/data_dict.md) plus four optional
    lexeme columns that :func:`iter_rows` joins in: ``ipa_lemma`` (the
    lemma's first Wiktionary transcription), ``alt_cells`` and
    ``resist_cells`` (``|``-joined E0 cells with an alternating P form,
    and with a competing N allomorph, at τ = dictionary) and
    ``hard_coronals`` (letters among t d n that are hard in the lemma's
    di ti ni). Without them a row has no evidence (it is prespecified
    under ``lexical``) and soft di ti ni.
    """
    form = (r.get("form") or "").strip()
    ending = (r.get("ending") or "").strip()
    cell = (r.get("cell") or "").strip()
    lemma = (r.get("lemma") or "").strip()
    cluster = (r.get("cluster") or "").strip()
    hard = frozenset((r.get("hard_coronals") or "").strip())
    stem_len = len(form) - len(ending) if ending else len(form)
    attested = (
        transcribe(form, boundary=stem_len, hard=hard).render() if form else ""
    )
    allomorph = (r.get("allomorph_class") or "").strip()
    return CzechRow(
        lemma=lemma,
        plural=form,
        pos=_POS.get((r.get("pos") or "").strip(), ""),
        stem_final=(r.get("stem_final") or "").strip(),
        opportunity=(f"{cell} -{ending}" if allomorph in (P, N) else "none"),
        mutation=_bool(r.get("mutation")),
        ipa_lemma=(r.get("ipa_lemma") or "").strip(),
        ipa_pl=attested,
        exception_reason=(r.get("exception_label") or "").strip(),
        cluster=ClusterTag.NONE,
        lexeme_id=(r.get("lexeme_id") or "").strip(),
        cell=cell,
        ending=ending,
        allomorph=allomorph,
        n_kind=(r.get("n_kind") or "").strip(),
        target_class=(r.get("target_class") or "").strip(),
        czech_cluster=cluster,
        orth_change=(r.get("orth_change") or "").strip(),
        gender_class=(r.get("gender_class") or "").strip(),
        suffix_cluster=suffix_cluster(lemma, cluster),
        alt_cells=_cells(r.get("alt_cells")),
        resist_cells=_cells(r.get("resist_cells")),
        hard_coronals=hard,
    )


def as_czech(row: LexRow) -> CzechRow:
    """``row`` as a :class:`CzechRow` (the hooks need its fields)."""
    if not isinstance(row, CzechRow):
        raise TypeError(
            "the Czech analysis needs CzechRow rows (row_factory / "
            f"iter_rows); got {type(row).__name__}"
        )
    return row


def _in_scope(r: Mapping[str, str]) -> bool:
    """The validated rows: primary, aligned, not rejected, spelled with
    letters only (numeral compounds such as 10letý have no
    transcription), a velar, r or coronal_inferred target, a P ending
    or a hard-paradigm N ending (soft-declension endings, kmeni, are
    not competitors)."""
    return (
        _bool(r.get("primary"))
        and _bool(r.get("aligned"))
        and (r.get("form") or "").isalpha()
        and (r.get("curated_status") or "") != "rejected"
        and r.get("pos") in _POS
        and r.get("target_class") in TARGET_CLASSES
        and r.get("stem_final") in TARGET_STEMS
        and (
            r.get("allomorph_class") == P
            or (r.get("allomorph_class") == N and r.get("n_kind") != "soft")
        )
    )


@dataclass(frozen=True, slots=True)
class LexemeInfo:
    """What a row needs from ``czech_lexicon_complete.csv``."""

    ipa_lemma: str
    alt_cells: tuple[str, ...]
    resist_cells: tuple[str, ...]
    hard_coronals: frozenset[str]


def data_paths(root: Path) -> tuple[Path, Path]:
    """``(czech_cells.csv, czech_lexicon_complete.csv)`` under ``root``."""
    data = LANGUAGE.data_dir(root)
    return (
        data / LANGUAGE.file_name("cells"),
        data / LANGUAGE.file_name("lexicon_complete"),
    )


def first_transcription(ipa_field: str) -> str:
    """The first Wiktionary transcription of a lexicon ``ipa`` cell."""
    for item in split_list(ipa_field):
        for variant in item.split("|"):
            text = variant.strip()
            if text:
                return text
    return ""


def _align(a: Sequence[str], b: Sequence[str]) -> list[tuple[str, str]]:
    """Levenshtein alignment of two token sequences (substitutions and
    matches only; gaps are dropped)."""
    rows, cols = len(a) + 1, len(b) + 1
    cost = [[0] * cols for _ in range(rows)]
    for i in range(rows):
        cost[i][0] = i
    for j in range(cols):
        cost[0][j] = j
    for i in range(1, rows):
        for j in range(1, cols):
            cost[i][j] = min(
                cost[i - 1][j] + 1,
                cost[i][j - 1] + 1,
                cost[i - 1][j - 1] + (a[i - 1] != b[j - 1]),
            )
    pairs: list[tuple[str, str]] = []
    i, j = len(a), len(b)
    while i > 0 and j > 0:
        if cost[i][j] == cost[i - 1][j - 1] + (a[i - 1] != b[j - 1]):
            pairs.append((a[i - 1], b[j - 1]))
            i, j = i - 1, j - 1
        elif cost[i][j] == cost[i - 1][j] + 1:
            i -= 1
        else:
            j -= 1
    pairs.reverse()
    return pairs


# Our soft reading → the hard segment (and letter) Wiktionary may write.
_HARDENED: Final[Mapping[str, str]] = MappingProxyType(
    {"c": "t", "ɟ": "d", "ɲ": "n"}
)


def hard_coronals(lemma: str, ipa: str) -> frozenset[str]:
    """Letters among t d n whose soft reading the lemma IPA contradicts.

    The lemma is transcribed with soft di ti ni and aligned with its
    Wiktionary IPA; where ours has ť ď ň and Wiktionary t d n (a
    ``{{cs-IPA|[dy]}}`` hint: diplomat [dɪplomat], politika), the letter
    is hard in this lexeme's stem.
    """
    if not lemma or not ipa:
        return frozenset()
    ours = [phonemic(s) for s in transcribe(lemma).segments]
    theirs = list(czech_tokens(ipa))
    letters: set[str] = set()
    for mine, wiki in _align(ours, theirs):
        if _HARDENED.get(mine) == wiki:
            letters.add(wiki)
    return frozenset(letters)


def load_lexemes(lexicon_complete: Path) -> dict[str, LexemeInfo]:
    """Lemma IPA, E0 evidence and hard coronals of every in-domain
    noun and adjective with a velar, r or coronal_inferred target."""
    lexemes: dict[str, LexemeInfo] = {}
    with lexicon_complete.open(encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f):
            if r.get("domain") != "in" or r.get("pos") not in _POS:
                continue
            if r.get("target_class") not in TARGET_CLASSES:
                continue
            ipa = first_transcription(r.get("ipa") or "")
            lexemes[r["lexeme_id"]] = LexemeInfo(
                ipa_lemma=ipa,
                alt_cells=_cells(r.get("p_cells_E0_dictionary")),
                resist_cells=_cells(r.get("n_cells_E0_dictionary")),
                hard_coronals=hard_coronals(r.get("lemma") or "", ipa),
            )
    return lexemes


def load_rows(root: Path) -> tuple[CzechRow, ...]:
    """Every validated row of the generated data under ``root``, sorted
    by (lexeme_id, cell, form), one row per (lexeme, cell, form)."""
    cells_path, complete_path = data_paths(root)
    lexemes = load_lexemes(complete_path)
    rows: dict[tuple[str, str, str], CzechRow] = {}
    with cells_path.open(encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f):
            if not _in_scope(r):
                continue
            key = (r["lexeme_id"], r["cell"], r["form"])
            if key in rows:
                continue  # the same form from another table
            info = lexemes.get(r["lexeme_id"])
            enriched = dict(r)
            if info is not None:
                enriched["ipa_lemma"] = info.ipa_lemma
                enriched["alt_cells"] = "|".join(info.alt_cells)
                enriched["resist_cells"] = "|".join(info.resist_cells)
                enriched["hard_coronals"] = "".join(sorted(info.hard_coronals))
            rows[key] = as_czech(row_factory(enriched))
    return tuple(rows[key] for key in sorted(rows))


def iter_rows(root: Path) -> Iterator[LexRow]:
    """The rows to validate (driver contract of
    ``scripts/validate_lang_rules.py``): :func:`load_rows`."""
    yield from load_rows(root)


# ---------------------------------------------------------------------------
# Specification policies
# ---------------------------------------------------------------------------

LEXICAL: Final[str] = "lexical"
ALL_UNDERSPECIFIED: Final[str] = "all_underspecified"
LEXICAL_LOO: Final[str] = "lexical_loo"
SPELLOUT: Final[str] = "spellout"
POLICIES: Final[tuple[str, ...]] = (
    LEXICAL,
    ALL_UNDERSPECIFIED,
    LEXICAL_LOO,
    SPELLOUT,
)

# Cells whose forms copy each other (lexpipe.langs.czech.features).
COPY_GROUPS: Final[tuple[frozenset[str], ...]] = (
    frozenset({"NOM.PL", "VOC.PL"}),
    frozenset({"NOM.PL.M.AN", "VOC.PL.M.AN"}),
)
FEMININE_COPY: Final[frozenset[str]] = frozenset({"DAT.SG", "LOC.SG"})


def copy_group(cell: str, gender_class: str) -> frozenset[str]:
    """The tested cell plus the cells that copy it."""
    if gender_class == "f" and cell in FEMININE_COPY:
        return FEMININE_COPY
    for group in COPY_GROUPS:
        if cell in group:
            return group
    return frozenset({cell})


def is_underspecified(row: LexRow, policy: str) -> bool:
    """Is the row's stem-final target underspecified under ``policy``?"""
    czech = as_czech(row)
    if policy == ALL_UNDERSPECIFIED:
        return True
    if policy == LEXICAL:
        return bool(czech.alt_cells)
    if policy == LEXICAL_LOO:
        excluded = copy_group(czech.cell, czech.gender_class)
        return any(cell not in excluded for cell in czech.alt_cells)
    if policy == SPELLOUT:
        return not czech.resist_cells
    raise ValueError(f"unknown policy {policy!r}; known: {POLICIES}")


# ---------------------------------------------------------------------------
# UR construction
# ---------------------------------------------------------------------------

# Data-layer mutation label → (base spelling, spellings of the change).
_MUTATION_SPELLINGS: Final[Mapping[str, tuple[str, tuple[str, ...]]]] = (
    MappingProxyType(
        {
            "k→c": ("k", ("c",)),
            "g→z": ("g", ("z",)),
            "h→z": ("h", ("z",)),
            "ch→š": ("ch", ("š",)),
            "r→ř": ("r", ("ř",)),
            "sk→sc": ("sk", ("sc",)),
            "sk→šť": ("sk", ("šť", "št")),
            "ck→čť": ("ck", ("čť", "čt")),
            "k→č": ("k", ("č",)),
            "g→ž": ("g", ("ž",)),
            "h→ž": ("h", ("ž",)),
        }
    )
)


def base_stem(row: LexRow) -> str | None:
    """The form's stem with the stem-final change undone (kluc- →
    kluk-, češt- → česk-), or None if the spelling does not show the
    change the data layer recorded."""
    czech = as_czech(row)
    stem = czech.stem_spelling
    if not czech.orth_change:
        return stem
    spelled = _MUTATION_SPELLINGS.get(czech.orth_change)
    if spelled is None:
        return None
    original, spellings = spelled
    lowered = stem.lower()
    for spelling in spellings:
        if lowered.endswith(spelling):
            return stem[: len(stem) - len(spelling)] + original
    return None


def palatalizer_of(row: LexRow) -> str | None:
    """The palatalizer the row's ending carries: PAL3 for the velar
    vocative -e, PAL1 for every other P ending, None for N endings."""
    czech = as_czech(row)
    if czech.allomorph != P:
        return None
    if czech.cell == "VOC.SG" and czech.target_class == VELAR:
        return PAL3
    return PAL1


def ur_labels(row: LexRow, policy: str) -> tuple[str, ...] | None:
    """The UR as inventory/underspecified labels, or None.

    The form with its stem-final change undone is transcribed with the
    stem edge left alone (``underlying=True``); allophones fold onto
    phonemes; the stem-final target (and the s/c of a -ský/-cký
    cluster) is swapped for its underspecified twin if the policy says
    so; the ending's palatalizer goes between stem and ending.
    """
    czech = as_czech(row)
    base = base_stem(czech)
    expected = BASE_SEGMENT.get(czech.stem_final)
    if base is None or expected is None or not base:
        return None
    transcription = transcribe(
        base + czech.ending,
        boundary=len(base),
        underlying=True,
        hard=czech.hard_coronals,
    )
    stem_segs, ending_segs = transcription.split(len(base))
    stem = [phonemic(s) for s in stem_segs]
    ending = [phonemic(s) for s in ending_segs]
    if not stem or stem[-1] != expected:
        return None
    if is_underspecified(czech, policy):
        if czech.suffix_cluster:
            sibilant = (
                SUFFIX_SIBILANT_FOR.get(stem[-2]) if len(stem) > 1 else None
            )
            if expected != "k" or sibilant is None:
                return None
            stem[-2], stem[-1] = sibilant, SUFFIX_VELAR
        else:
            stem[-1] = UNDERSPEC_FOR[expected]
    pal = palatalizer_of(czech)
    return tuple(stem) + ((pal,) if pal is not None else ()) + tuple(ending)


def build_ur(
    row: LexRow, inventory: FeatureInventory, policy: str
) -> Word | None:
    """The UR :data:`Word` of ``row`` under ``policy`` (None when a label
    is not in the inventory)."""
    labels = ur_labels(row, policy)
    if labels is None:
        return None
    segments: list[Segment] = []
    for label in labels:
        if label in PALATALIZER_FEATURES:
            segments.append(palatalizer(label))
        elif label in inventory.base or label in inventory.underspec:
            segments.append(inventory.segment(label))
        else:
            return None
    return tuple(segments)


@dataclass(frozen=True, slots=True)
class UrBuilder:
    """The ``ur_builder`` hook of one policy (pure, picklable)."""

    policy: str

    def __call__(
        self,
        row: LexRow,
        inventory: FeatureInventory,
        /,
    ) -> Word | None:
        return build_ur(row, inventory, self.policy)


def expected_firings(row: LexRow, policy: str) -> frozenset[str]:
    """The palatalization rules that should fire under ``policy``.

    The analysis's claim: an underspecified target before a palatalizer
    alternates by the rule of its class; a prespecified target or an N
    ending fires nothing. Never read off the attested outcome.
    """
    czech = as_czech(row)
    pal = palatalizer_of(czech)
    if pal is None or not is_underspecified(czech, policy):
        return frozenset()
    if pal == PAL3:
        return frozenset({BIG_PAL.name})
    if czech.stem_final == "ch":
        return frozenset({X_PAL.name})
    if czech.stem_final in ("k", "g", "h"):
        if czech.suffix_cluster:
            return frozenset({CORONAL_PAL.name, S_PAL.name})
        return frozenset({VELAR_PAL.name})
    return frozenset({CORONAL_PAL.name})


@dataclass(frozen=True, slots=True)
class ExpectedFirings:
    """The ``expected_firings`` hook of one policy."""

    policy: str

    def __call__(self, row: LexRow) -> frozenset[str]:
        return expected_firings(row, self.policy)


@dataclass(frozen=True, slots=True)
class FallbackPredictor:
    """The ``fallback_predict`` hook: the policy's own prediction for a
    row whose UR cannot be built (never the data)."""

    policy: str

    def __call__(self, row: LexRow) -> bool:
        return bool(expected_firings(row, self.policy))


# ---------------------------------------------------------------------------
# Categoriser, cluster buckets and Czech mismatch classes
# ---------------------------------------------------------------------------


def categorise(
    row: LexRow,
    predicted_sr: str,
    attested_field: str,
    boolean_predicted: bool,
    /,
    *,
    ur_built: bool,
    normalisers: NormaliserCascade,
) -> Categorisation:
    """The ``categoriser`` hook, on the shared vocabulary.

    :class:`RefinedCategory` has no Czech labels: a row is
    UR_BUILD_FAILED or RULE_FAILURE (the cascade has no labelled steps,
    so no NORM_* category arises). :func:`mismatch_class` gives the
    Czech cause.
    """
    if not ur_built:
        return Categorisation(RefinedCategory.UR_BUILD_FAILED, (), 999.0, "")
    cmp = compare_ipa(predicted_sr, attested_field, normalisers=normalisers)
    return Categorisation(
        RefinedCategory.RULE_FAILURE,
        cmp.normalisation_applied,
        cmp.edit_distance,
        cmp.attested_variant,
    )


def row_class(row: LexRow) -> str:
    """``velar``, ``r`` or ``coronal_inferred``, with ``(adj)`` for
    adjectives."""
    czech = as_czech(row)
    suffix = " (adj)" if czech.pos == "ADJ" else ""
    return czech.target_class + suffix


def cluster_key(row: LexRow) -> str | None:
    """The cluster-report bucket: the cluster and ending of a clustered
    stem (sk+í is -ský), the velar vocative, or None."""
    czech = as_czech(row)
    if czech.cell == "VOC.SG" and czech.allomorph == P:
        return f"VOC {czech.target_class}"
    if not czech.czech_cluster:
        return None
    suffix = " (suffix)" if czech.suffix_cluster else ""
    return f"{czech.czech_cluster}+{czech.ending}{suffix}"


class MismatchClass(StrEnum):
    """Why a Czech row failed (the shared categories cannot say)."""

    UR_NOT_BUILT = "ur-not-built"
    UNEXPLAINED_FORM = "unexplained-form"
    OVERAPPLICATION = "rules-change-unchanged-form"
    MIXED_STEM_ALTERNATES = "mixed-stem-alternates"
    NO_EVIDENCE_ALTERNATES = "no-evidence-alternates"
    CLUSTER_OUTCOME = "cluster-outcome"
    CORONAL_INFERRED = "coronal-inferred-ipa"
    OTHER_IPA = "other-ipa"


MISMATCH_GLOSSES: Final[Mapping[MismatchClass, str]] = MappingProxyType(
    {
        MismatchClass.UR_NOT_BUILT: (
            "no UR (the spelling does not show the recorded change, or "
            "a letter is outside the inventory); scored by the fallback"
        ),
        MismatchClass.UNEXPLAINED_FORM: (
            "P ending without the change in the data (the data layer's "
            "'unexplained' forms, mostly cs.wiktionary table errors); "
            "the rules palatalize"
        ),
        MismatchClass.OVERAPPLICATION: (
            "the rules change a form the data shows unchanged (not an "
            "'unexplained' row)"
        ),
        MismatchClass.MIXED_STEM_ALTERNATES: (
            "the policy left the target prespecified, but the form "
            "alternates, and the stem also takes a competing N allomorph "
            "(in another cell, or beside this form: park has v parku "
            "but parcích): choosing an N allomorph does not mark a "
            "non-alternating root"
        ),
        MismatchClass.NO_EVIDENCE_ALTERNATES: (
            "the policy left the target prespecified for want of "
            "evidence (no other cell alternates or competes: the "
            "feminine DAT/LOC.SG, the adjective NOM.PL), but the form "
            "alternates"
        ),
        MismatchClass.CLUSTER_OUTCOME: (
            "boolean right, IPA wrong at an sk/ck/st/zd cluster " "(sc vs šť)"
        ),
        MismatchClass.CORONAL_INFERRED: (
            "boolean right, IPA wrong on a t/d/n row (change inferred "
            "from spelling)"
        ),
        MismatchClass.OTHER_IPA: "boolean right, IPA wrong elsewhere",
    }
)


def mismatch_class(result: RowResult) -> MismatchClass | None:
    """The Czech cause of a failed row (None for a full match)."""
    if result.boolean_match and result.ipa_match:
        return None
    czech = as_czech(result.row)
    if result.ur is None:
        return MismatchClass.UR_NOT_BUILT
    if not result.boolean_match:
        if result.boolean_predicted:
            if czech.exception_reason == "unexplained":
                return MismatchClass.UNEXPLAINED_FORM
            return MismatchClass.OVERAPPLICATION
        if czech.resist_cells:
            return MismatchClass.MIXED_STEM_ALTERNATES
        return MismatchClass.NO_EVIDENCE_ALTERNATES
    if czech.czech_cluster in ("sk", "ck", "st", "zd"):
        return MismatchClass.CLUSTER_OUTCOME
    if czech.target_class == CORONAL:
        return MismatchClass.CORONAL_INFERRED
    return MismatchClass.OTHER_IPA


# ---------------------------------------------------------------------------
# Profiles
# ---------------------------------------------------------------------------


def build_variant(policy: str) -> AnalysisProfile:
    """The :class:`AnalysisProfile` of one specification policy."""
    if policy not in POLICIES:
        raise ValueError(f"unknown policy {policy!r}; known: {POLICIES}")
    return AnalysisProfile(
        name=f"czech-palatalization ({policy})",
        inventory=load_inventory(),
        patches=PATCHES,
        underspec=UNDERSPEC,
        rules=RULES,
        inventory_json=inventory_path(_INVENTORY_NAME),
        target_stems=TARGET_STEMS,
        palatalization_rule_names=PALATALIZATION_RULE_NAMES,
        expected_firings=ExpectedFirings(policy),
        fallback_predict=FallbackPredictor(policy),
        ur_builder=UrBuilder(policy),
        cluster_key=cluster_key,
        normalisers=CZECH_NORMALISERS,
        categoriser=categorise,
        row_factory=row_factory,
        tokenizer=CZECH_TOKENIZER,
    )


def build_profile() -> AnalysisProfile:
    """The default analysis: the ``lexical`` policy (Gorman & Winstead)."""
    return build_variant(LEXICAL)


def _all_underspecified() -> AnalysisProfile:
    return build_variant(ALL_UNDERSPECIFIED)


def _lexical_loo() -> AnalysisProfile:
    return build_variant(LEXICAL_LOO)


def _spellout() -> AnalysisProfile:
    return build_variant(SPELLOUT)


PROFILE_VARIANTS: Final[Mapping[str, Callable[[], AnalysisProfile]]] = (
    MappingProxyType(
        {
            LEXICAL: build_profile,
            ALL_UNDERSPECIFIED: _all_underspecified,
            LEXICAL_LOO: _lexical_loo,
            SPELLOUT: _spellout,
        }
    )
)
"""The variants ``scripts/validate_lang_rules.py --variant`` accepts."""


# ---------------------------------------------------------------------------
# Summary across variants
# ---------------------------------------------------------------------------


def _rate(part: int, whole: int) -> str:
    return f"{100 * part / whole:6.2f}%" if whole else "     -"


@dataclass(frozen=True, slots=True)
class _Tally:
    rows: int
    boolean: int
    ipa: int
    both: int


def _tally(results: Sequence[RowResult]) -> _Tally:
    return _Tally(
        rows=len(results),
        boolean=sum(r.boolean_match for r in results),
        ipa=sum(r.ipa_match for r in results),
        both=sum(r.boolean_match and r.ipa_match for r in results),
    )


def _table(
    title: str,
    groups: Mapping[str, Sequence[RowResult]],
) -> list[str]:
    lines = [title, ""]
    width = max([len(k) for k in groups] + [10])
    lines.append(
        f"  {'':{width}}  {'rows':>6}  {'boolean':>8}  {'IPA':>8}  "
        f"{'both':>8}"
    )
    for key in sorted(groups):
        t = _tally(groups[key])
        lines.append(
            f"  {key:{width}}  {t.rows:6d}  {_rate(t.boolean, t.rows):>8}  "
            f"{_rate(t.ipa, t.rows):>8}  {_rate(t.both, t.rows):>8}"
        )
    lines.append("")
    return lines


def _grouped(
    results: Iterable[RowResult], key: Callable[[RowResult], str]
) -> dict[str, list[RowResult]]:
    groups: dict[str, list[RowResult]] = defaultdict(list)
    for result in results:
        groups[key(result)].append(result)
    return groups


def _cell_key(result: RowResult) -> str:
    czech = as_czech(result.row)
    return f"{row_class(czech):22} {czech.cell:11} {czech.allomorph}"


def _final_key(result: RowResult) -> str:
    czech = as_czech(result.row)
    return f"{czech.stem_final:2} {czech.allomorph}"


def _example(result: RowResult) -> str:
    czech = as_czech(result.row)
    ur = "".join(s.label for s in result.ur) if result.ur else "-"
    return (
        f"{czech.lemma} → {czech.plural} ({czech.cell} -{czech.ending}): "
        f"UR /{ur}/ predicted [{result.predicted_sr}] attested "
        f"[{result.attested_first.lstrip(STRESS)}]"
    )


def render_variant_summary(
    report: ValidationReport, *, examples: int = 6
) -> str:
    """Boolean and IPA match rates by class × cell, by stem final, by
    cluster, and the Czech mismatch classes of one variant."""
    results = report.results
    total = _tally(results)
    lines = [
        f"=== {report.profile_name} ===",
        "",
        f"rows {total.rows:,}; boolean match {_rate(total.boolean, total.rows)}"
        f"; IPA match {_rate(total.ipa, total.rows)}; both "
        f"{_rate(total.both, total.rows)}",
        "",
    ]
    lines += _table(
        "By class × cell × ending class (P palatalizing, N not):",
        _grouped(results, _cell_key),
    )
    lines += _table("By class:", _grouped(results, lambda r: row_class(r.row)))
    lines += _table(
        "By stem final × ending class:", _grouped(results, _final_key)
    )
    clustered = [r for r in results if cluster_key(r.row) is not None]
    lines += _table(
        "Clusters and the vocative (bucket = cluster+ending):",
        _grouped(clustered, lambda r: cluster_key(r.row) or ""),
    )
    classes: dict[MismatchClass, list[RowResult]] = defaultdict(list)
    for result in results:
        cls = mismatch_class(result)
        if cls is not None:
            classes[cls].append(result)
    lines.append("Mismatch classes (rows failing boolean or IPA):")
    lines.append("")
    for cls in sorted(classes, key=lambda c: (-len(classes[c]), c.value)):
        items = classes[cls]
        lines.append(
            f"  [{cls.value}] n = {len(items)}: {MISMATCH_GLOSSES[cls]}"
        )
        by_class = Counter(_cell_key(r) for r in items)
        for key, count in sorted(by_class.items(), key=lambda kv: -kv[1])[:6]:
            lines.append(f"      {count:6d}  {key}")
        for result in items[:examples]:
            lines.append(f"    {_example(result)}")
        lines.append("")
    return "\n".join(lines)


@dataclass(frozen=True, slots=True)
class G2pAgreement:
    """Agreement of :func:`transcribe` with the lemma transcriptions."""

    lemmas: int
    exact: int
    exact_with_hard: int


def g2p_agreement(lexicon_complete: Path) -> G2pAgreement:
    """:func:`transcribe` vs the en.wiktionary lemma IPA of in-domain
    one-word lexemes with a single transcription (stress and syllable
    dots removed; the hard di ti ni of the lemma's own IPA applied in
    the second count)."""
    lemmas = exact = exact_hard = 0
    with lexicon_complete.open(encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f):
            if r.get("edition") != "en" or r.get("domain") != "in":
                continue
            items = split_list(r.get("ipa") or "")
            lemma = r.get("lemma") or ""
            if len(items) != 1 or "|" in items[0] or not lemma.isalpha():
                continue
            wiki = items[0].strip("[]/").replace(STRESS, "").replace(".", "")
            ours = "".join(transcribe(lemma).segments)
            hard = hard_coronals(lemma, items[0])
            ours_hard = "".join(transcribe(lemma, hard=hard).segments)
            lemmas += 1
            exact += ours == wiki
            exact_hard += ours_hard == wiki
    return G2pAgreement(lemmas, exact, exact_hard)


_SUMMARY_HEADER: Final[str] = """\
Czech palatalization: Logical Phonology analysis, all variants
(phonology.analyses.czech_palatalization; docs/czech/phonology.md)

Source: Czech entries of English Wiktionary (and of Czech Wiktionary
for lexemes English Wiktionary lacks) via kaikki.org/wiktextract
(Ylonen 2022), CC BY-SA 4.0. Rows: data/czech/czech_cells.csv, primary
aligned rows of velar, r and coronal_inferred stems with a P ending or
a hard-paradigm N ending (nouns and adjectives), one per lexeme × cell
× form. Evidence for the policies: E0 at τ = dictionary (no corpus
counts are used). Attested IPA: the spelling transcribed with the
cs-IPA conventions; compared modulo stress, syllable marks, the
allophones r̝̊ ŋ ɱ ɣ d͡z and degemination. t/d/n rows
(coronal_inferred) compare with a change inferred from the spelling.
"""


def render_summary(
    reports: Mapping[str, ValidationReport],
    *,
    agreement: G2pAgreement | None = None,
    examples: int = 6,
) -> str:
    """The cross-variant summary text (``render_variant_summary`` of
    every report, after a header and the transcription check)."""
    lines = [_SUMMARY_HEADER]
    if agreement is not None:
        exact = _rate(agreement.exact, agreement.lemmas).strip()
        hard = _rate(agreement.exact_with_hard, agreement.lemmas).strip()
        lines += [
            "Transcription check: transcribe() equals the en.wiktionary",
            f"lemma IPA of {agreement.exact:,} of {agreement.lemmas:,} "
            f"one-word lemmas ({exact});",
            "with the hard di ti ni of the lemma's own IPA, "
            f"{agreement.exact_with_hard:,} ({hard}).",
            "",
        ]
    for policy, report in reports.items():
        lines.append(f"Policy: {policy}")
        lines.append(render_variant_summary(report, examples=examples))
    return "\n".join(lines).rstrip() + "\n"


def trace(lemma: str, rows: Iterable[LexRow], policy: str = LEXICAL) -> str:
    """Derivation tables of every row of ``lemma`` under ``policy``."""
    profile = build_variant(policy)
    out: list[str] = []
    inv_base = profile.inventory.base_segments()
    for row in rows:
        if row.lemma != lemma:
            continue
        result = evaluate_row(row, profile)
        out.append(render_trace(TraceResult(row=row, result=result), inv_base))
    return "\n".join(out) if out else f"{lemma}: no validated rows\n"


def main(argv: Sequence[str] | None = None) -> int:
    """Run every variant and write the cross-variant summary::

        python -m phonology.analyses.czech_palatalization [--root DIR]
            [--out FILE] [--trace LEMMA ...]

    writes ``analysis/czech/czech_rule_validation_summary.txt`` under
    ``--root``; the per-variant reports and mismatch CSVs are
    ``scripts/validate_lang_rules.py --lang cs [--variant NAME]``.
    """
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument("--root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--trace", nargs="*", default=())
    args = parser.parse_args(argv)
    rows = load_rows(args.root)
    if args.trace:
        for lemma in args.trace:
            print(trace(lemma, rows))
        return 0
    reports = {
        policy: run_validation(build_variant(policy), rows)
        for policy in POLICIES
    }
    _, complete = data_paths(args.root)
    text = render_summary(reports, agreement=g2p_agreement(complete))
    out = args.out or (
        LANGUAGE.analysis_dir(args.root)
        / LANGUAGE.file_name("rule_validation_summary", ".txt")
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    print(text, end="")
    print(f"\nSummary written to {out}")
    return 0


__all__: Final[tuple[str, ...]] = (
    "AFFRICATE_RELEASE",
    "ALL_UNDERSPECIFIED",
    "BIG_PAL",
    "CORONAL_PAL",
    "CORONAL_PLACE_DEFAULT",
    "CZECH_NORMALISERS",
    "CZECH_TOKENIZER",
    "CzechRow",
    "DORSAL_DEFAULT",
    "DORSAL_STOP_DEFAULT",
    "LEXICAL",
    "LEXICAL_LOO",
    "MismatchClass",
    "PAL1",
    "PAL3",
    "PAL_DELETION",
    "PALATALIZATION_RULES",
    "PALATALIZATION_RULE_NAMES",
    "PATCHES",
    "POLICIES",
    "PROFILE_VARIANTS",
    "RULES",
    "SIBILANT_CONTINUANT",
    "SPELLOUT",
    "S_PAL",
    "TARGET_STEMS",
    "TRILL_OBSTRUENT",
    "TRILL_SONORANT",
    "Transcription",
    "UNDERSPEC",
    "VELAR_PAL",
    "X_PAL",
    "build_profile",
    "build_ur",
    "build_variant",
    "categorise",
    "cluster_key",
    "expected_firings",
    "hard_coronals",
    "is_underspecified",
    "iter_rows",
    "load_inventory",
    "load_rows",
    "mismatch_class",
    "render_summary",
    "row_factory",
    "strip_czech",
    "transcribe",
    "ur_labels",
)


if __name__ == "__main__":
    sys.exit(main())
