"""Stage ``lexicon``: Italian noun and adjective lexemes (R4, R8, R12, R25).

Inputs: the ``extract`` projections ``italian_en_lexemes.jsonl.gz``
(every noun/adjective lemma entry of the en split, with its
``lexeme_id``) and ``italian_it_entries.jsonl.gz`` (it.wiktionary).

Outputs
-------
``data/italian/italian_lexicon_raw.csv`` (tracked, written through
:func:`lexpipe.csvio.write_tracked`): the lean raw lexicon,
:data:`RAW_FIELDS` — what was extracted from the en.wiktionary entry
and nothing derived:

  lexeme_id    the identifier (R4)
  word, pos    lemma and part of speech
  gender       base of the first gender code (m f mf mfbysense n ?)
  head_args    the head templates as written (:func:`lexpipe.kaikki.head_args`:
               ``1=m;2=amici``; the template name is written only when
               it is not the part of speech's standard it-noun/it-adj)
  ipa_raw      the primary transcription as written
  etymology_summary  the deciding etymology template, summarized
               (:func:`lexpipe.etymology.summary`: ``inh:la:amīcus``, ``af``)
  first_gloss  first gloss of the first sense
  source       ``wiktionary_en``

One row per lexeme, sorted by ``lexeme_id`` (R5), except out-of-domain
proper names (domain ``capitalized``) and multiword expressions, which
only the full table holds (:data:`RAW_EXCLUDED_DOMAINS`).

``data/italian/italian_lexicon_full.csv`` (gitignored): every lexeme
with every column of the stage (:data:`LEXICON_FIELDS`; ``source`` =
``wiktionary_en`` or ``wiktionary_en|wiktionary_it`` when it.wiktionary
columns are filled). Later stages read this table. Column groups:

  identity    lexeme_id, root_id, word, pos, etymology_number,
              head_template, first_gloss
  domain      ``domain`` (in_domain or the first out-of-domain class)
              and ``domain_flags`` (all of them, in this order):
              multiword, hyphenated, affix, capitalized (proper or
              capitalised common nouns), nonalphabetic, abbreviation
              (every lexeme sense tagged abbreviation, acronym or
              initialism), clipping (a ``clipping`` template, an
              etymology starting "Clipping of"/"Short for", or every
              sense tagged clipping)
  gender      gender_raw (argument 1 as written), gender (first code's
              base: m f mf mfbysense n ?), genders (every code, in
              order: f|m for *eco*), gender_labels, by_sense,
              plurale_tantum (a ``-p`` code)
  plurals     plural_status, plural_provenance (explicit |
              module_default | mixed | none), countability, invariable,
              pure_default (no plural argument), plural_drift (the
              drift-test status), pl_forms / mpl_forms / fpl_forms (in
              source order), plural_labels (``form:label/label`` for
              labelled variants), feminine_forms / masculine_forms
              (counterparts stated by the head)
  phonology   ipa_raw (the primary transcription as written: the first
              without a regional/traditional/nonstandard tag), ipa
              (normalised body, stress and syllable dots kept, R8),
              ipa_kind, ipa_count, stress_from_end (1 final … 4 earlier),
              stress_class, stress_source (ipa | hyphenation, the en
              accented hyphenation | it_hyphenation, it.wiktionary's,
              used only when en has neither),
              stress_variants (distinct stresses over all IPAs),
              stress_missing, final_post_tonic (the final syllable — the
              one holding a stem-final -co/-go velar — follows the
              stressed one), hyphenation (en, accented)
  etymology   etymon_lang, etymon, first_transmission (the template kind
              deciding the stratum; :mod:`.etymology`)
  it-edition  it_entries (matching lemma entries: same word and part of
              speech), it_gender, it_forms (``cell:form``), it_invariable,
              it_hyphenation, it_stress_from_end, stress_agrees

``root_id`` (R4) is :func:`lexpipe.ids.root_id` over the key
``[stem, etymon key, stress, (etymology label)]``: the lemma's match key
without its final vowel, :func:`.etymology.etymon_key` (amīcus and
amīca agree, antīquus and antīcus do not), the stress position (antico
/anˈti.ko/ vs àntico /ˈan.ti.ko/), and — only when no etymon is known —
the page's etymology number, so that two etymology sections of one page
without templates stay apart. A feminine noun that names its masculine
counterpart (``m=`` argument, or a "female equivalent of X" sense) takes
the root key of that lexeme, which merges gender pairs (amico/amica).
POS conversions (amico N/ADJ) share the key by construction.

``data/italian/italian_plural_forms.csv`` (R12): one row per plural
form (:mod:`.plurals` documents the columns), sorted by lexeme_id,
cell, item.

``data/italian/italian_plural_drift.csv`` (R12): the drift test. For
every lexeme whose plurals can be compared, the resolved plurals against
the plurals kaikki rendered; ``scope`` pure_default_noun is the test the
design asks for (port(lemma, gender) == kaikki plural for every noun
without a plural argument), pure_default_adj extends it to adjectives
and with_arguments checks the argument parser. ``status`` match |
mismatch | port_error (the module has no default) | not_ported
(multiword).

``data/italian/italian_lexicon_summary.tsv``: counts.
"""

from __future__ import annotations

import dataclasses
import re
import unicodedata
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any, Final

from lexpipe import csvio
from lexpipe import etymology as core_etymology
from lexpipe import ids, kaikki, normalize
from lexpipe.context import RunContext, run_stage
from lexpipe.langs.italian import (
    etymology,
    extract,
    heads,
    itcommon,
    keys,
    plurals,
    store,
)
from lexpipe.langs.italian.layout import Layout

STAGE: Final[str] = "lexicon"
# 2: lean tracked raw lexicon + gitignored full table; shared etymology
# parser (lexpipe.etymology).
# 3: checked hyphenation stress (stress_rejected, ico_stressed), the
#    plurale_tantum and plural_form domains, superlatives out of
#    it_forms, root keys with the final-vowel class, the full etymon and
#    POS-conversion stress; agreement counts in the summary.
VERSION: Final[int] = 3
SOURCE_EN: Final[str] = "wiktionary_en"
SOURCE_IT: Final[str] = "wiktionary_it"
IN_DOMAIN: Final[str] = "in_domain"

LEXICON_FIELDS: Final[tuple[str, ...]] = (
    "lexeme_id",
    "root_id",
    "word",
    "pos",
    "etymology_number",
    "head_template",
    "domain",
    "domain_flags",
    "gender_raw",
    "gender",
    "genders",
    "gender_labels",
    "by_sense",
    "plurale_tantum",
    "plural_status",
    "plural_provenance",
    "countability",
    "invariable",
    "pure_default",
    "plural_drift",
    "pl_forms",
    "mpl_forms",
    "fpl_forms",
    "plural_labels",
    "feminine_forms",
    "masculine_forms",
    "ipa_raw",
    "ipa",
    "ipa_kind",
    "ipa_count",
    "stress_from_end",
    "stress_class",
    "stress_source",
    "stress_variants",
    "stress_missing",
    "final_post_tonic",
    "ico_stressed",
    "hyphenation",
    "etymon_lang",
    "etymon",
    "first_transmission",
    "it_entries",
    "it_gender",
    "it_forms",
    "it_invariable",
    "it_hyphenation",
    "it_stress_from_end",
    "stress_agrees",
    "stress_rejected",
    "first_gloss",
    "source",
)
# The tracked raw lexicon: extracted columns only (module docstring).
RAW_FIELDS: Final[tuple[str, ...]] = (
    "lexeme_id",
    "word",
    "pos",
    "gender",
    "head_args",
    "ipa_raw",
    "etymology_summary",
    "first_gloss",
    "source",
)
# Out-of-domain lexemes the tracked raw lexicon leaves out: proper
# names (capitalized) and multiword expressions.
RAW_EXCLUDED_DOMAINS: Final[frozenset[str]] = frozenset(
    {"capitalized", "multiword"}
)
# The standard head template of each part of speech (head_args leaves
# its name out).
DEFAULT_HEADS: Final[Mapping[str, str]] = {"noun": "it-noun", "adj": "it-adj"}
PLURAL_FIELDS: Final[tuple[str, ...]] = (
    "lexeme_id",
    "cell",
    "form",
    "provenance",
    "code",
    "gender",
    "qualifiers",
    "labels",
    "raw_arg",
    "tag_error",
    "item",
    "source",
)
DRIFT_FIELDS: Final[tuple[str, ...]] = (
    "lexeme_id",
    "word",
    "pos",
    "gender",
    "scope",
    "status",
    "resolved",
    "rendered",
    "module_revision",
)
_CELL_ORDER: Final[Mapping[str, int]] = {
    cell: index for index, cell in enumerate(heads.CELLS)
}
# Sound tags that mark a transcription as not the standard one.
NONSTANDARD_SOUND_TAGS: Final[frozenset[str]] = frozenset(
    {
        "traditional",
        "regional",
        "nonstandard",
        "proscribed",
        "dialectal",
        "colloquial",
        "uncommon",
        "sometimes",
        "Latinate",
        "Northern",
        "Milan",
        "Tuscany",
        "Central-and-Southern-Italy",
        "Switzerland",
        "Italy",
    }
)
_ABBREVIATION_TAGS: Final[frozenset[str]] = frozenset(
    {"abbreviation", "acronym", "initialism"}
)
_CLIPPING_TEXT_RE: Final[re.Pattern[str]] = re.compile(
    r"\A\s*(?:Clipping|Clipped form|Short|Shortening) (?:of|for)\b"
)
_FEMALE_GLOSS_RE: Final[re.Pattern[str]] = re.compile(
    r"\A(?:female|feminine) (?:equivalent )?of\b", re.IGNORECASE
)
# A gloss that makes an entry a plural form (domain plural_form): "plural
# of X", with an optional qualifier before it ("masculine", "Old
# Italian", "The"), or "(plural of X)" / "; plural of X" after a gloss.
_PLURAL_OF_RE: Final[re.Pattern[str]] = re.compile(
    r"\A(?:\([^)]*\)\s*)?(?:(?:the|old|old italian|archaic|obsolete|"
    r"masculine|feminine)\s+)*plural of\b|\(plural of\b|;\s*plural of\b",
    re.IGNORECASE,
)
_HYPHEN_SEPARATORS: Final[str] = "‧"
# Syllable separators found in hyphenations: the hyphenation point, a
# middle dot, hyphens and spaces ("post - so - vi - è - ti - co",
# "ra‧dio‧e‧stè‧si co").
_SYLLABLE_SPLIT_RE: Final[re.Pattern[str]] = re.compile(r"[‧·\-\s]+")
_STRESS_ACCENTS: Final[frozenset[str]] = frozenset({"̀", "́"})
_MAX_STRESS: Final[int] = 4
# Vowels of an Italian transcription (the nucleus of a syllable).
_IPA_VOWELS: Final[frozenset[str]] = frozenset("aeiouɛɔ")
# Lemma endings whose stressed vowel matters for R26 (-ico/-iaco).
_ICO_ENDINGS: Final[tuple[str, ...]] = ("ico", "ica")
# Why a hyphenation's stress was not used (column stress_rejected).
REJECT_OTHER_WORD: Final[str] = "other_word"
REJECT_ICO_ACCENT: Final[str] = "ico_accent_vs_template"


# ---------------------------------------------------------------------------
# Stress
# ---------------------------------------------------------------------------


def hyphenation_syllables(parts: Sequence[str]) -> list[str]:
    """The syllables of a hyphenation given as parts, each part split
    again at the separators editors also use (hyphens, spaces, middle
    dots: :data:`_SYLLABLE_SPLIT_RE`); pieces without a letter (a
    footnote number: ``ri‧le‧gà‧to 1``) are dropped."""
    out: list[str] = []
    for part in parts:
        out.extend(
            piece
            for piece in _SYLLABLE_SPLIT_RE.split(part)
            if any(ch.isalpha() for ch in piece)
        )
    return out


def _letters(text: str) -> str:
    """Lowercase letters without diacritics (apostrophes and separators
    dropped), for comparing a hyphenation with its word."""
    decomposed = unicodedata.normalize("NFD", text.lower())
    return "".join(
        ch
        for ch in decomposed
        if not unicodedata.combining(ch) and ch.isalpha()
    )


def spells(word: str, syllables: Sequence[str]) -> bool:
    """Do the hyphenation's syllables spell *word*? Accents, case and
    apostrophes are ignored, and a different final vowel is allowed
    (the hyphenation of agonistico on agonistica keeps the stress
    position); another word (bi‧o‧ce‧nò‧si on biocenotico) is not."""
    written, target = _letters("".join(syllables)), _letters(word)
    if written == target:
        return True
    return (
        len(written) == len(target) > 1
        and written[:-1] == target[:-1]
        and written[-1] in "aeiou"
        and target[-1] in "aeiou"
    )


def _accented(syllables: Sequence[str]) -> list[int]:
    return [
        index
        for index, syllable in enumerate(syllables)
        if _STRESS_ACCENTS & set(unicodedata.normalize("NFD", syllable))
    ]


def syllable_stress(parts: Sequence[str]) -> int | None:
    """Stress from an accented hyphenation (lù|ce → 2, cit|tà → 1).

    None unless exactly one syllable carries a grave or acute accent
    (a monosyllable needs none). Parts are split again at hyphens and
    spaces (:func:`hyphenation_syllables`).
    """
    syllables = hyphenation_syllables(parts)
    if not syllables:
        return None
    accented = _accented(syllables)
    if len(accented) == 1:
        return min(len(syllables) - accented[0], _MAX_STRESS)
    if not accented and len(syllables) == 1:
        return 1
    return None


def accented_vowel(parts: Sequence[str]) -> str:
    """The base letter that carries the accent of a hyphenation ("" if
    no single syllable is accented): ``i`` for ``bo‧so‧nì‧co``, ``a``
    for ``lài‧co``."""
    syllables = hyphenation_syllables(parts)
    accented = _accented(syllables)
    if len(accented) != 1:
        return ""
    chars = unicodedata.normalize("NFD", syllables[accented[0]])
    for index, ch in enumerate(chars[1:], 1):
        if ch in _STRESS_ACCENTS:
            return chars[index - 1].lower()
    return ""


def _first_alternative(text: str) -> list[str]:
    """The syllables of the first alternative of an en hyphenation."""
    first = re.split(r",\s*", text.strip())[0]
    return hyphenation_syllables(first.split(_HYPHEN_SEPARATORS))


def hyphenation_stress(text: str, word: str = "") -> int | None:
    """:func:`syllable_stress` of an en hyphenation such as ``a‧mì‧co``
    (the first alternative when several are listed); None when *word*
    is given and the hyphenation spells another word (:func:`spells`)."""
    syllables = _first_alternative(text)
    if word and not spells(word, syllables):
        return None
    return syllable_stress(syllables)


def ipa_vowel_of(syllable: str) -> str:
    """The first vowel of an IPA syllable (its nucleus in Italian: a
    following i/u is a glide, written j/w or as a falling diphthong)."""
    for ch in unicodedata.normalize("NFD", syllable):
        if ch in _IPA_VOWELS:
            return ch
    return ""


def is_ico_word(word: str) -> bool:
    """Does the lemma end in -ico/-ica (-iaco/-iaca included)?"""
    return keys.match_key(word).endswith(_ICO_ENDINGS)


def ico_stressed(
    word: str, stress_from_end: int | None, vowel: str
) -> bool | None:
    """Does the i of a final -ico/-ica carry the stress (amìco)?

    *vowel* is the stressed vowel (:func:`ipa_vowel_of` of the stressed
    IPA syllable, or :func:`accented_vowel` of the hyphenation). True
    only for a penultimate stress on that i; a penultimate stress on
    the vowel before it (arcaico /arˈkaj.ko/, laico /ˈlaj.ko/: the i is
    a glide) is False, as is any other position. None when the word
    does not end in -ico/-ica or its stress is unknown.
    """
    if not is_ico_word(word) or stress_from_end is None:
        return None
    if stress_from_end != 2:
        return False
    return vowel == "i" if vowel else None


def ipa_stress(ipa: normalize.Ipa) -> int | None:
    """Stress position of a transcription; 1 for a monosyllable
    written without a stress mark."""
    found = normalize.stress_from_end(ipa)
    if found is not None:
        return found
    pieces = normalize.syllabify(ipa)
    if (
        len(pieces) == 1
        and normalize.is_syllabified(ipa)
        and " " not in ipa.text.strip()
    ):
        return 1
    return None


@dataclass(frozen=True, slots=True)
class Stress:
    """The phonological columns of one lexeme."""

    ipa_raw: str = ""
    ipa: str = ""
    ipa_kind: str = ""
    ipa_count: int = 0
    stress_from_end: int | None = None
    stress_source: str = ""
    stress_variants: tuple[int, ...] = ()
    final_post_tonic: bool | None = None
    hyphenation: str = ""
    # The vowel that carries the stress (ipa_vowel_of / accented_vowel).
    stressed_vowel: str = ""
    # Hyphenation stresses not used, as "source:reason".
    rejected: tuple[str, ...] = ()

    @property
    def stress_class(self) -> str:
        if self.stress_from_end is None:
            return ""
        return normalize.STRESS_CLASS_NAMES[self.stress_from_end]


def _ipa_stressed_vowel(ipa: normalize.Ipa) -> str:
    """The vowel of the stressed syllable of a transcription (the only
    syllable of an unmarked monosyllable)."""
    pieces = normalize.syllabify(ipa)
    marked = [p for p in pieces if p.startswith(normalize.PRIMARY_STRESS)]
    if marked:
        return ipa_vowel_of(marked[-1])
    return ipa_vowel_of(pieces[0]) if len(pieces) == 1 else ""


def _primary_sound(
    sounds: Sequence[Mapping[str, Any]],
) -> tuple[str, normalize.Ipa | None]:
    parsed = [
        (str(s["ipa"]), keys.normalizer().ipa(str(s["ipa"])), s)
        for s in sounds
        if isinstance(s.get("ipa"), str)
    ]
    usable = [(raw, ipa, s) for raw, ipa, s in parsed if ipa is not None]
    for raw, ipa, sound in usable:
        tags = set(sound.get("tags") or ()) | set(sound.get("raw_tags") or ())
        if not tags & NONSTANDARD_SOUND_TAGS:
            return raw, ipa
    if usable:
        return usable[0][0], usable[0][1]
    return (parsed[0][0], None) if parsed else ("", None)


def hyphenation_choice(
    word: str,
    syllables: Sequence[str],
    source: str,
    found: Sequence[etymology.Step] = (),
    lemmas: Mapping[str, bool] = MappingProxyType({}),
) -> tuple[int | None, str, str]:
    """``(stress, stressed vowel, rejection)`` of a hyphenation used as
    the stress source *source* (``hyphenation`` or ``it_hyphenation``).

    The stress is not used (rejection ``<source>:<reason>``) when the
    hyphenation spells another word (:data:`REJECT_OTHER_WORD`:
    bi‧o‧ce‧nò‧si on biocenotico), or when it accents the i of a final
    -ico/-ica although the word is formed with the unstressed -ico
    suffix (:data:`REJECT_ICO_ACCENT`): its etymology says X + -ico
    (iperone + -ico), or the word is an adjective (*lemmas* maps lemma
    keys to whether one of their entries is an adjective) with an
    Italian base among the lemmas (bosone → bo‧so‧nì‧co, tallio →
    tal‧lì‧co, :func:`.etymology.ico_base`). No such formation is
    stressed on the suffix's i (simplex nouns can be: nasìca < Latin
    nāsīca); a hyphenation is the fallback for words without a
    transcription, so a mistyped accent has nothing to override it.
    """
    if not syllables:
        return None, "", ""
    if not spells(word, syllables):
        return None, "", f"{source}:{REJECT_OTHER_WORD}"
    stress = syllable_stress(syllables)
    if stress is None:
        return None, "", ""
    vowel = accented_vowel(syllables)
    key = keys.match_key(word)
    if ico_stressed(word, stress, vowel) is True and (
        etymology.ico_template(found)
        or (lemmas.get(key, False) and bool(etymology.ico_base(key, lemmas)))
    ):
        return None, "", f"{source}:{REJECT_ICO_ACCENT}"
    return stress, vowel, ""


def stress_of(
    record: Mapping[str, Any],
    found: Sequence[etymology.Step] = (),
    lemmas: Mapping[str, bool] = MappingProxyType({}),
) -> Stress:
    """Stress columns from the projected sounds and hyphenation (the
    en hyphenation only when it passes :func:`hyphenation_choice`;
    *found* are the entry's etymology steps, *lemmas* the lemma keys)."""
    sounds = [s for s in record.get("sounds") or () if isinstance(s, Mapping)]
    raw, ipa = _primary_sound(sounds)
    variants = sorted(
        {
            found_stress
            for s in sounds
            if (parsed := keys.normalizer().ipa(str(s.get("ipa", ""))))
            is not None
            and (found_stress := ipa_stress(parsed)) is not None
        }
    )
    word = str(record.get("word") or "")
    hyphenation = next(iter(record.get("hyphenation") or ()), "")
    stress = ipa_stress(ipa) if ipa is not None else None
    source = "ipa" if stress is not None else ""
    vowel = (
        _ipa_stressed_vowel(ipa)
        if stress is not None and ipa is not None
        else ""
    )
    rejected: list[str] = []
    if stress is None and hyphenation:
        stress, vowel, problem = hyphenation_choice(
            word, _first_alternative(hyphenation), "hyphenation", found, lemmas
        )
        source = "hyphenation" if stress is not None else ""
        if problem:
            rejected.append(problem)
    post_tonic: bool | None = None
    if source == "ipa" and ipa is not None:
        post_tonic = normalize.post_tonic(ipa, -1)
    elif stress is not None:
        post_tonic = stress == 2
    return Stress(
        ipa_raw=raw,
        ipa=ipa.text if ipa is not None else "",
        ipa_kind=str(ipa.kind) if ipa is not None else "",
        ipa_count=len(sounds),
        stress_from_end=stress,
        stress_source=source,
        stress_variants=tuple(variants),
        final_post_tonic=post_tonic,
        hyphenation=hyphenation,
        stressed_vowel=vowel,
        rejected=tuple(rejected),
    )


# ---------------------------------------------------------------------------
# Domain, counterparts, it-edition
# ---------------------------------------------------------------------------


def _lexeme_senses(record: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    return [
        s
        for s in record.get("senses") or ()
        if isinstance(s, Mapping)
        and not s.get("form_of")
        and not s.get("alt_of")
        and not {"form-of", "alt-of"} & set(s.get("tags") or ())
    ]


def _first_glosses(senses: Sequence[Mapping[str, Any]]) -> list[str]:
    return [
        str(glosses[0])
        for sense in senses
        if (glosses := sense.get("glosses")) and isinstance(glosses, list)
    ]


def is_plural_form_entry(record: Mapping[str, Any]) -> bool:
    """Is the entry a plural of another word? Its first sense says so
    ("plural of alimento", "(Old Italian) plural of casa"; the senses
    after it are the lexicalized plural's meanings: alimenti 'alimony',
    casalinghi 'housewares'), or every lexeme sense does ("crostini
    (plural of crostino)"). Such entries are plurals, not lemmas
    (domain ``plural_form``): their headword is not a singular."""
    senses = [s for s in record.get("senses") or () if isinstance(s, Mapping)]
    first = _first_glosses(senses[:1])
    if first and _PLURAL_OF_RE.search(first[0]):
        return True
    glosses = _first_glosses(_lexeme_senses(record))
    return bool(glosses) and all(
        _PLURAL_OF_RE.search(gloss) for gloss in glosses
    )


def domain_flags(
    record: Mapping[str, Any],
    found: Sequence[etymology.Step],
    info: plurals.PluralInfo | None = None,
) -> list[str]:
    """Out-of-domain classes of a lexeme (module docstring order).
    ``plurale_tantum`` needs the head's gender codes (*info*)."""
    word = str(record.get("word") or "")
    flags: list[str] = []
    if " " in word:
        flags.append("multiword")
    if "-" in word.strip("-"):
        flags.append("hyphenated")
    if word.startswith("-") or word.endswith("-"):
        flags.append("affix")
    if any(ch.isupper() for ch in word):
        flags.append("capitalized")
    if not all(ch.isalpha() or ch in " -'’" for ch in word):
        flags.append("nonalphabetic")
    senses = _lexeme_senses(record)
    if senses and all(
        _ABBREVIATION_TAGS & set(s.get("tags") or ()) for s in senses
    ):
        flags.append("abbreviation")
    text = str(record.get("etymology_text") or "")
    clipped = (
        any(step.kind == "clipping" for step in found)
        or bool(_CLIPPING_TEXT_RE.search(text))
        or (
            senses and all("clipping" in (s.get("tags") or ()) for s in senses)
        )
    )
    if clipped:
        flags.append("clipping")
    if (
        info is not None
        and info.genders
        and all(code.plural_only for code in info.genders)
    ):
        flags.append("plurale_tantum")
    if is_plural_form_entry(record):
        flags.append("plural_form")
    return flags


def female_counterparts(
    record: Mapping[str, Any], info: plurals.PluralInfo
) -> tuple[str, ...]:
    """Masculine words a feminine noun names as its counterpart."""
    found: dict[str, None] = dict.fromkeys(info.masculines)
    for sense in record.get("senses") or ():
        if not isinstance(sense, Mapping):
            continue
        glosses = [g for g in sense.get("glosses") or () if isinstance(g, str)]
        if glosses and _FEMALE_GLOSS_RE.search(glosses[0]):
            for word in sense.get("form_of") or ():
                found.setdefault(str(word), None)
    return tuple(found)


@dataclass(frozen=True, slots=True)
class ItEdition:
    """What it.wiktionary says about a word and part of speech."""

    entries: int = 0
    gender: str = ""
    forms: tuple[str, ...] = ()
    invariable: bool = False
    hyphenation: str = ""
    stress_from_end: int | None = None
    ipa: str = ""
    # The hyphenation's syllables (split again at hyphens and spaces).
    syllables: tuple[str, ...] = ()


# Forms of an it.wiktionary paradigm that are not plurals of the lemma
# although tagged plural (timidissimi: an absolute superlative).
_NOT_PLURAL_TAGS: Final[frozenset[str]] = frozenset(
    {"superlative", "comparative"}
)
_PL_RAW_TAG_RE: Final[re.Pattern[str]] = re.compile(r"\Apl\.?:\s*(.+)\Z")


def _it_cell(tags: set[str]) -> str:
    if "plural" not in tags or tags & _NOT_PLURAL_TAGS:
        return ""
    if "masculine" in tags and "feminine" not in tags:
        return heads.CELL_MPL
    if "feminine" in tags and "masculine" not in tags:
        return heads.CELL_FPL
    return heads.CELL_PL


def it_edition(entries: Sequence[Mapping[str, Any]]) -> ItEdition:
    """Merge the it.wiktionary lemma entries of one word and POS."""
    if not entries:
        return ItEdition()
    genders: dict[str, None] = {}
    forms: dict[str, None] = {}
    invariable = False
    for entry in entries:
        tags = set(entry.get("tags") or ())
        if "masculine" in tags:
            genders.setdefault("m", None)
        if "feminine" in tags:
            genders.setdefault("f", None)
        invariable = invariable or "invariable" in tags
        for form in entry.get("forms") or ():
            cell = _it_cell(set(form.get("tags") or ()))
            if cell and form.get("form"):
                forms.setdefault(f"{cell}:{form['form']}", None)
        for sense in entry.get("senses") or ():
            for raw in sense.get("raw_tags") or ():
                match = _PL_RAW_TAG_RE.match(str(raw).strip())
                if match:
                    for text in match.group(1).split(","):
                        if text.strip():
                            forms.setdefault(f"pl:{text.strip()}", None)
    first = entries[0]
    parts: list[str] = [
        str(part) for part in next(iter(first.get("hyphenations") or ()), [])
    ]
    syllables = hyphenation_syllables(parts)
    return ItEdition(
        entries=len(entries),
        gender="|".join(genders),
        forms=tuple(forms),
        invariable=invariable,
        hyphenation="‧".join(parts),
        stress_from_end=syllable_stress(syllables) if syllables else None,
        ipa=next(iter(first.get("ipa") or ()), ""),
        syllables=tuple(syllables),
    )


def load_it_edition(path: Path) -> dict[tuple[str, str], list[dict[str, Any]]]:
    """it.wiktionary lemma entries keyed by (word, pos)."""
    meta, records = store.read_jsonl(path)
    if meta.get("format") != extract.IT_ENTRIES_FORMAT:
        raise store.StaleDerivativeError(f"{path}: unexpected format {meta}")
    out: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        if "form-of" in (record.get("tags") or ()):
            continue
        out[(record["word"], record["pos"])].append(record)
    return dict(out)


def with_it_fallback(
    stress: Stress,
    edition: ItEdition,
    word: str = "",
    found: Sequence[etymology.Step] = (),
    lemmas: Mapping[str, bool] = MappingProxyType({}),
) -> Stress:
    """Use it.wiktionary's accented hyphenation when en gives no stress
    (``stress_source = it_hyphenation``), if it passes
    :func:`hyphenation_choice` for *word* (a rejection is recorded in
    ``rejected``)."""
    if stress.stress_from_end is not None or not edition.syllables:
        return stress
    found_stress, vowel, problem = hyphenation_choice(
        word, edition.syllables, "it_hyphenation", found, lemmas
    )
    rejected = stress.rejected + ((problem,) if problem else ())
    if found_stress is None:
        return dataclasses.replace(stress, rejected=rejected)
    return dataclasses.replace(
        stress,
        stress_from_end=found_stress,
        stress_source="it_hyphenation",
        final_post_tonic=found_stress == 2,
        stressed_vowel=vowel,
        rejected=rejected,
    )


# ---------------------------------------------------------------------------
# One lexeme
# ---------------------------------------------------------------------------


@dataclass
class Lexeme:
    """One lexicon row in the making, plus what the root key needs."""

    row: dict[str, Any]
    plural_rows: list[dict[str, Any]]
    drift: plurals.Drift | None
    stem_key: str
    etymon_key: str
    stress: str
    ety_label: str
    counterparts: tuple[str, ...]
    gender: str
    raw: dict[str, Any] = field(default_factory=dict)
    root_parts: list[str] = field(default_factory=list)
    vowel_class: str = ""

    @property
    def tracked(self) -> bool:
        """Does the tracked raw lexicon hold this lexeme?"""
        return self.row["domain"] not in RAW_EXCLUDED_DOMAINS


def _labels(info: plurals.PluralInfo) -> list[str]:
    out: dict[str, None] = {}
    for row in info.rows:
        if row.form and row.labels:
            out.setdefault(
                f"{row.form}:{'/'.join(row.labels)}".replace("|", "/"), None
            )
    return list(out)


def _gender_columns(info: plurals.PluralInfo, raw: str) -> dict[str, Any]:
    genders = info.genders
    labels = [
        f"{code.code}:{'/'.join(code.labels)}"
        for code in genders
        if code.labels
    ]
    return {
        "gender_raw": raw,
        "gender": genders[0].base if genders else "",
        "genders": [code.code for code in genders],
        "gender_labels": labels,
        "by_sense": any(code.base == "mfbysense" for code in genders),
        "plurale_tantum": info.plurale_tantum,
    }


def _plural_columns(
    info: plurals.PluralInfo, drift: plurals.Drift | None
) -> dict[str, Any]:
    return {
        "plural_status": info.status,
        "plural_provenance": info.provenance,
        "countability": info.countability,
        "invariable": info.invariable,
        "pure_default": info.pure_default,
        "plural_drift": drift.status if drift else "",
        "pl_forms": list(info.cell_forms(heads.CELL_PL)),
        "mpl_forms": list(info.cell_forms(heads.CELL_MPL)),
        "fpl_forms": list(info.cell_forms(heads.CELL_FPL)),
        "plural_labels": _labels(info),
        "feminine_forms": list(info.feminines),
        "masculine_forms": list(info.masculines),
    }


def _stress_columns(stress: Stress, word: str) -> dict[str, Any]:
    return {
        "ipa_raw": stress.ipa_raw,
        "ipa": stress.ipa,
        "ipa_kind": stress.ipa_kind,
        "ipa_count": stress.ipa_count,
        "stress_from_end": stress.stress_from_end,
        "stress_class": stress.stress_class,
        "stress_source": stress.stress_source,
        "stress_variants": (
            list(stress.stress_variants)
            if len(stress.stress_variants) > 1
            else []
        ),
        "stress_missing": stress.stress_from_end is None,
        "final_post_tonic": stress.final_post_tonic,
        "ico_stressed": ico_stressed(
            word, stress.stress_from_end, stress.stressed_vowel
        ),
        "hyphenation": stress.hyphenation,
        "stress_rejected": list(stress.rejected),
    }


def _it_columns(edition: ItEdition, stress: Stress) -> dict[str, Any]:
    agrees: bool | None = None
    if (
        edition.stress_from_end is not None
        and stress.stress_from_end is not None
        and stress.stress_source != "it_hyphenation"
    ):
        agrees = edition.stress_from_end == stress.stress_from_end
    return {
        "it_entries": edition.entries,
        "it_gender": edition.gender,
        "it_forms": list(edition.forms),
        "it_invariable": edition.invariable if edition.entries else None,
        "it_hyphenation": edition.hyphenation,
        "it_stress_from_end": edition.stress_from_end,
        "stress_agrees": agrees,
    }


def _plural_row(lexeme_id: str, row: plurals.PluralRow) -> dict[str, Any]:
    return {
        "lexeme_id": lexeme_id,
        "cell": row.cell,
        "form": row.form,
        "provenance": row.provenance,
        "code": row.code,
        "gender": row.gender,
        "qualifiers": list(row.qualifiers),
        "labels": list(row.labels),
        "raw_arg": row.raw_arg,
        "tag_error": row.tag_error,
        "item": row.item,
        "source": SOURCE_EN,
    }


def build_lexeme(
    record: Mapping[str, Any],
    it_index: Mapping[tuple[str, str], Sequence[Mapping[str, Any]]],
    lemmas: Mapping[str, bool] = MappingProxyType({}),
) -> Lexeme:
    """Everything the lexicon knows about one projected entry (*lemmas*:
    every noun/adjective lemma key → has an adjective entry, for the
    stress check of :func:`hyphenation_choice`)."""
    word = str(record["word"])
    pos = str(record["pos"])
    lexeme_id = str(record["lexeme_id"])
    head_name, args = heads.head_template(record)
    info = plurals.resolve_plurals(word, record.get("head_templates") or ())
    drift = plurals.drift(info, record)
    found = etymology.steps(record.get("etymology_templates") or ())
    flags = domain_flags(record, found, info)
    source_etymon = etymology.etymon(found)
    first = etymology.first_transmission(found)
    edition = it_edition(it_index.get((word, pos), ()))
    stress = with_it_fallback(
        stress_of(record, found, lemmas),
        edition,
        word=word,
        found=found,
        lemmas=lemmas,
    )
    row: dict[str, Any] = {
        "lexeme_id": lexeme_id,
        "word": word,
        "pos": pos,
        "etymology_number": record.get("etymology_number", "0"),
        "head_template": head_name,
        "domain": flags[0] if flags else IN_DOMAIN,
        "domain_flags": flags,
        **_gender_columns(info, args.get("1", "") if head_name else ""),
        **_plural_columns(info, drift),
        **_stress_columns(stress, word),
        "etymon_lang": source_etymon[0] if source_etymon else "",
        "etymon": source_etymon[1] if source_etymon else "",
        "first_transmission": first.kind if first else "",
        **_it_columns(edition, stress),
        "first_gloss": ids.first_gloss(record),
        "source": (
            SOURCE_EN if not edition.entries else f"{SOURCE_EN}|{SOURCE_IT}"
        ),
    }
    gender = row["gender"]
    raw = {
        "lexeme_id": lexeme_id,
        "word": word,
        "pos": pos,
        "gender": gender,
        "head_args": kaikki.head_args(record, DEFAULT_HEADS.get(pos, "")),
        "ipa_raw": row["ipa_raw"],
        "etymology_summary": core_etymology.summary(
            found, etymology.PROFILE.code
        ),
        "first_gloss": row["first_gloss"],
        "source": SOURCE_EN,
    }
    return Lexeme(
        row=row,
        raw=raw,
        plural_rows=[_plural_row(lexeme_id, r) for r in info.rows],
        drift=drift,
        stem_key=keys.strip_final_vowel(keys.match_key(word)),
        etymon_key=etymology.etymon_key(found),
        stress=(
            ""
            if stress.stress_from_end is None
            else str(stress.stress_from_end)
        ),
        ety_label=str(record.get("etymology_number", "0")),
        counterparts=(
            female_counterparts(record, info) if gender == "f" else ()
        ),
        gender=gender,
        vowel_class=final_vowel_class(word),
    )


# ---------------------------------------------------------------------------
# Root identifiers
# ---------------------------------------------------------------------------


def final_vowel_class(word: str) -> str:
    """The root-key class of a lemma's final vowel: ``o/a`` for -o and
    -a (the gender pairs and the -ica nouns of -ico adjectives share a
    root), the vowel itself for -e, -i and -u (lattice and lattico,
    dolce and dolco are different words), "" without a final vowel."""
    key = keys.match_key(word)
    if keys.strip_final_vowel(key) == key:
        return ""
    return "o/a" if key[-1] in "oa" else key[-1]


def own_root_parts(lexeme: Lexeme, stress: str | None = None) -> list[str]:
    """``[stem, final-vowel class, etymon key, stress, (ety label)]``
    (module docstring); *stress* replaces the lexeme's own (a POS
    conversion's, :func:`conversion_stresses`)."""
    parts = [
        lexeme.stem_key,
        lexeme.vowel_class,
        lexeme.etymon_key,
        lexeme.stress if stress is None else stress,
    ]
    if not lexeme.etymon_key and lexeme.ety_label not in ("", "0"):
        parts.append(f"ety{lexeme.ety_label}")
    return parts


def conversion_stresses(lexemes: Sequence[Lexeme]) -> dict[int, str]:
    """The stress a lexeme without one takes for its root key: the one
    stress of the other entries of the same word and page etymology
    (a POS conversion: psicotico N without a transcription, psicotico
    ADJ /psiˈkɔ.ti.ko/), keyed by ``id(lexeme)``. Words whose entries
    have several stresses (antico, àntico) fill nothing."""
    known: dict[tuple[str, str], set[str]] = defaultdict(set)
    for lexeme in lexemes:
        if lexeme.stress:
            known[(lexeme.row["word"], lexeme.ety_label)].add(lexeme.stress)
    out: dict[int, str] = {}
    for lexeme in lexemes:
        found = known.get((lexeme.row["word"], lexeme.ety_label), set())
        if not lexeme.stress and len(found) == 1:
            out[id(lexeme)] = next(iter(found))
    return out


def assign_root_ids(lang: str, lexemes: Sequence[Lexeme]) -> None:
    """Fill ``root_id`` of every lexeme (gender pairs merged)."""
    by_word: dict[str, list[Lexeme]] = defaultdict(list)
    for lexeme in lexemes:
        by_word[lexeme.row["word"]].append(lexeme)
    filled = conversion_stresses(lexemes)
    for lexeme in lexemes:
        lexeme.root_parts = own_root_parts(lexeme, filled.get(id(lexeme)))
    for lexeme in lexemes:
        partner = _counterpart(lexeme, by_word)
        if partner is not None:
            lexeme.root_parts = own_root_parts(
                partner, filled.get(id(partner))
            )
    for lexeme in lexemes:
        lexeme.row["root_id"] = ids.root_id(
            lang, {"word": lexeme.row["word"]}, fixed_key(lexeme.root_parts)
        )


def fixed_key(parts: Sequence[str]) -> ids.RootKeyFn:
    """A :data:`lexpipe.ids.RootKeyFn` that returns *parts* whatever the
    entry (the parts are computed across lexemes beforehand)."""
    frozen = list(parts)

    def key(_entry: Mapping[str, Any]) -> Sequence[str]:
        return frozen

    return key


def _counterpart(
    lexeme: Lexeme, by_word: Mapping[str, Sequence[Lexeme]]
) -> Lexeme | None:
    """The masculine lexeme a feminine noun points to, if unique enough:
    candidates with the same stress first, then nouns, then the
    smallest lexeme_id."""
    candidates = [
        other
        for word in lexeme.counterparts
        for other in by_word.get(word, ())
        if other is not lexeme
    ]
    if not candidates:
        return None

    def rank(other: Lexeme) -> tuple[int, int, str]:
        same_stress = 0 if other.stress == lexeme.stress else 1
        noun = 0 if other.row["pos"] == "noun" else 1
        return same_stress, noun, str(other.row["lexeme_id"])

    return min(candidates, key=rank)


# ---------------------------------------------------------------------------
# The stage
# ---------------------------------------------------------------------------


def load_lexemes(
    layout: Layout,
) -> list[Lexeme]:
    """Build every lexeme from the projections (sorted by lexeme_id)."""
    it_index = load_it_edition(layout.it_entries)
    meta, records = store.read_jsonl(layout.en_lexemes)
    if meta.get("format") != extract.EN_LEXEMES_FORMAT:
        raise store.StaleDerivativeError(
            f"{layout.en_lexemes}: unexpected format {meta}"
        )
    _, first_pass = store.read_jsonl(layout.en_lexemes)
    lemmas: dict[str, bool] = {}
    for item in first_pass:
        key = keys.match_key(str(item["word"]))
        lemmas[key] = lemmas.get(key, False) or item.get("pos") == "adj"
    lexemes = [build_lexeme(record, it_index, lemmas) for record in records]
    lexemes.sort(key=lambda lexeme: str(lexeme.row["lexeme_id"]))
    return lexemes


def _drift_rows(lexemes: Iterable[Lexeme]) -> Iterable[dict[str, Any]]:
    for lexeme in lexemes:
        drift = lexeme.drift
        if drift is None:
            continue
        yield {
            "lexeme_id": lexeme.row["lexeme_id"],
            "word": lexeme.row["word"],
            "pos": lexeme.row["pos"],
            "gender": lexeme.row["gender"],
            "scope": drift.scope,
            "status": drift.status,
            "resolved": list(drift.resolved),
            "rendered": list(drift.rendered),
            "module_revision": itcommon.MODULE_REVISION,
        }


def _sorted_plural_rows(lexemes: Iterable[Lexeme]) -> list[dict[str, Any]]:
    rows = [row for lexeme in lexemes for row in lexeme.plural_rows]
    rows.sort(
        key=lambda r: (r["lexeme_id"], _CELL_ORDER[r["cell"]], r["item"])
    )
    return rows


def summary_counts(lexemes: Sequence[Lexeme]) -> Counter[str]:
    """Aggregate counts for ``italian_lexicon_summary.tsv``."""
    counts: Counter[str] = Counter()
    for lexeme in lexemes:
        row = lexeme.row
        counts["lexemes"] += 1
        counts["tracked_raw"] += lexeme.tracked
        counts[f"pos:{row['pos']}"] += 1
        counts[f"domain:{row['domain']}"] += 1
        counts[f"plural_status:{row['plural_status']}"] += 1
        counts[f"plural_provenance:{row['plural_provenance']}"] += 1
        counts[f"stress_source:{row['stress_source'] or 'missing'}"] += 1
        if row["stress_agrees"] is not None:
            # en stress (ipa or hyphenation) vs it.wiktionary's, by source
            agree = "agree" if row["stress_agrees"] else "disagree"
            counts[f"stress_vs_it:{row['stress_source']}:{agree}"] += 1
        for reason in row["stress_rejected"]:
            counts[f"stress_rejected:{reason}"] += 1
        counts[f"it_matched:{bool(row['it_entries'])}"] += 1
        counts[f"gender:{row['gender'] or 'none'}"] += 1
        if lexeme.drift is not None:
            counts[f"drift:{lexeme.drift.scope}:{lexeme.drift.status}"] += 1
        for plural in lexeme.plural_rows:
            counts[f"plural_rows:{plural['provenance']}"] += 1
            if plural["tag_error"]:
                counts[f"plural_tag_error:{plural['tag_error']}"] += 1
    counts["roots"] = len({lexeme.row["root_id"] for lexeme in lexemes})
    return counts


def _write_summary(ctx: RunContext, path: Path, counts: Counter[str]) -> None:
    lines = ["key\tvalue"]
    lines.extend(f"{key}\t{counts[key]}" for key in sorted(counts))
    ctx.write_text(path, "\n".join(lines) + "\n")


def _body(ctx: RunContext, layout: Layout) -> int:
    lexemes = load_lexemes(layout)
    assign_root_ids(ctx.language.code, lexemes)
    csvio.write_rows(
        layout.lexicon_full,
        (lexeme.row for lexeme in lexemes),
        LEXICON_FIELDS,
        ctx=ctx,
    )
    csvio.write_tracked(
        layout.lexicon_raw,
        (lexeme.raw for lexeme in lexemes if lexeme.tracked),
        RAW_FIELDS,
        ctx=ctx,
    )
    csvio.write_rows(
        layout.plural_forms,
        _sorted_plural_rows(lexemes),
        PLURAL_FIELDS,
        ctx=ctx,
    )
    csvio.write_rows(
        layout.plural_drift, _drift_rows(lexemes), DRIFT_FIELDS, ctx=ctx
    )
    counts = summary_counts(lexemes)
    _write_summary(ctx, layout.lexicon_summary, counts)
    drift = {
        key.split(":", 1)[1]: value
        for key, value in counts.items()
        if key.startswith("drift:pure_default_noun:")
    }
    print(
        f"  {counts['lexemes']:,} lexemes ({counts['pos:noun']:,} nouns, "
        f"{counts['pos:adj']:,} adjectives), {counts['roots']:,} roots; "
        f"drift test (pure-default nouns): {drift}"
    )
    return 0


def outputs(layout: Layout) -> list[Path]:
    return [
        layout.lexicon_full,
        layout.lexicon_raw,
        layout.plural_forms,
        layout.plural_drift,
        layout.lexicon_summary,
    ]


def run(ctx: RunContext) -> int | None:
    """Run the ``lexicon`` stage (skipped when its manifest matches)."""
    layout = Layout(ctx)
    params = {
        "version": VERSION,
        "module_it_common_revision": itcommon.MODULE_REVISION,
    }
    return run_stage(
        ctx,
        STAGE,
        inputs=[layout.en_lexemes, layout.it_entries],
        params=params,
        outputs=outputs(layout),
        body=lambda: _body(ctx, layout),
    )
