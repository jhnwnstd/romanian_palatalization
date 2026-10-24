"""Morphological evidence for R26 suffix targets (design rule R26).

R26 asks for suffix-internal targets (-ník, -ík, -ák, -ek, -ka, -ička,
-ský/-cký, -log, -urg, -graf) to be identified **morphologically**,
not orthographically: knížka's k is the diminutive suffix's, ruka's
and banka's k belong to the root. The lexicon stage records the
evidence kind in ``suffix_evidence`` (:data:`EVIDENCE_KINDS`, in
priority order) and the features stage sets ``target_is_suffix`` when
the stem-final target lies inside the suffix
(:func:`lexpipe.langs.czech.ortho.target_is_suffix`):

``etymology``  an English Wiktionary etymology template names the
               suffix (``{{af|cs|kniha|-ka}}``);
``category``   "Czech terms suffixed with -X";
``formof``     a derivational form-of sense names a base that is a
               lexeme (klouček 'diminutive of kluk', Češka 'female
               equivalent of Čech'; :mod:`.extract`);
``headword``   another lexeme's headword line names the word as its
               feminine, masculine or diminutive (Čech f=Češka,
               kniha dim=knížka);
``adjective``  a hard adjective in -ský/-cký (-ický included): Czech
               has no simplex adjective root in -sk-/-ck-, so the k is
               always the adjectival suffix's (český, německý,
               typický);
``base``       the word minus the suffix is the stem of another lexeme
               (noun, adjective or verb), plain or with the first
               palatalization undone and vowel length, ě/e and the
               soft coronals folded (sedlák ← sedlo, ručník ← ruka,
               klouček ← kluk), the base stem has at least
               :data:`MIN_BASE_STEM` letters, and the word's own
               Wiktionary stratum is not ``inherited`` or ``borrowed
               whole`` (banka is borrowed from Italian banca, not
               ban + -ka; matka is inherited whole).

-log, -urg and -graf keep their orthographic evidence (R27: the
suffixoids of whole-word loans). Everything here is a pure function of
spellings and of the lexeme sets passed in.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import Final

ETYMOLOGY: Final[str] = "etymology"
CATEGORY: Final[str] = "category"
FORMOF: Final[str] = "formof"
HEADWORD: Final[str] = "headword"
ADJECTIVE: Final[str] = "adjective"
BASE: Final[str] = "base"
ORTHOGRAPHIC: Final[str] = "orthographic"
EVIDENCE_KINDS: Final[tuple[str, ...]] = (
    ETYMOLOGY,
    CATEGORY,
    FORMOF,
    HEADWORD,
    ADJECTIVE,
    BASE,
    ORTHOGRAPHIC,
)
MIN_BASE_STEM: Final[int] = 3
# Suffixes whose base may be found by decomposition (``base``), with
# the letters the suffix takes off the word. -ička is also tried as -ka
# on a base in -ice (lavička ← lavic-e).
NOMINAL_SUFFIX_CUTS: Final[Mapping[str, tuple[int, ...]]] = {
    "ička": (4, 2),
    "ník": (3,),
    "ík": (2,),
    "ák": (2,),
    "ek": (2,),
    "ka": (2,),
}
ADJECTIVE_SUFFIXES: Final[frozenset[str]] = frozenset({"ský", "cký"})
# Strata that make a word opaque: inherited whole or borrowed whole.
OPAQUE_STRATA: Final[frozenset[str]] = frozenset(
    {"inherited", "borrowed whole"}
)

_FOLD: Final[dict[int, int]] = str.maketrans("áéěíóúůýňťď", "aeeiouoyntd")
# First-palatalization outputs and the finals they come from (the
# derivational alternations of R26 suffixes: kluk → klouček, ruka →
# ručník, Čech → Češka, kniha → knížka, bratr → bratříček).
_UNDO: Final[tuple[tuple[str, tuple[str, ...]], ...]] = (
    ("šť", ("sk", "st")),
    ("čť", ("ck",)),
    ("č", ("k", "c", "t")),
    ("ž", ("h", "g", "z", "d")),
    ("š", ("ch", "s")),
    ("ř", ("r",)),
    ("c", ("k",)),
    ("z", ("h", "g")),
)
_VERB_ENDINGS: Final[tuple[str, ...]] = (
    "ovat",
    "nout",
    "at",
    "át",
    "it",
    "ít",
    "et",
    "ět",
    "out",
)
_NOUN_ENDINGS: Final[tuple[str, ...]] = (
    "ium",
    "um",
    "us",
    "os",
    "on",
    "a",
    "e",
    "ě",
    "o",
    "í",
    "y",
    "i",
    "u",
)


def fold(text: str) -> str:
    """Comparison key: lowercase, vowel length and ě folded, ou → u,
    ň ť ď → n t d (Czech spells them n t d before ě/i)."""
    return text.lower().translate(_FOLD).replace("ou", "u")


def _strip(word: str, endings: Iterable[str]) -> str:
    lowered = word.lower()
    for ending in endings:
        if lowered.endswith(ending) and len(lowered) > len(ending) + 1:
            return lowered[: -len(ending)]
    return lowered


def base_stems(lemma: str, pos: str) -> set[str]:
    """Folded stems a derivative of *lemma* can build on.

    Nouns lose a final vowel (or -um/-us/-on loan endings) and a
    fleeting e (Turek → turk); adjectives -ý/-í; verbs their
    infinitive suffix.
    """
    lowered = lemma.lower()
    out: set[str] = set()
    if pos == "adj":
        out.add(_strip(lowered, ("ý", "í")))
    elif pos == "verb":
        out.add(_strip(lowered, _VERB_ENDINGS))
    else:
        stem = _strip(lowered, _NOUN_ENDINGS)
        out.add(stem)
        if stem.endswith(("ek", "ec", "el", "en", "er")) and len(stem) > 3:
            out.add(stem[:-2] + stem[-1])  # fleeting e: Turek → turk
    return {fold(stem) for stem in out if stem}


def undone(rest: str) -> set[str]:
    """*rest* as written and with a final first-palatalization output
    undone (klouč → klouk; kníž → knih, knig, kniz, knid), folded."""
    lowered = rest.lower()
    out = {fold(lowered)}
    for output, finals in _UNDO:
        if lowered.endswith(output):
            head = lowered[: -len(output)]
            out.update(fold(head + final) for final in finals)
            break
    return out


@dataclass
class BaseIndex:
    """Folded base stems of the lexicon (:func:`base_stems`)."""

    stems: dict[str, set[str]] = field(
        default_factory=lambda: defaultdict(set)
    )

    def add(self, lemma: str, pos: str) -> None:
        for stem in base_stems(lemma, pos):
            if len(stem) >= MIN_BASE_STEM:
                self.stems[stem].add(lemma)

    def bases(self, word: str, suffix: str) -> list[str]:
        """Lemmas whose stem is *word* minus *suffix* (sorted; the word
        itself excluded)."""
        cuts = NOMINAL_SUFFIX_CUTS.get(suffix)
        if not cuts:
            return []
        lowered = word.lower()
        found: set[str] = set()
        for cut in cuts:
            if len(lowered) <= cut:
                continue
            rest = lowered[:-cut]
            for key in undone(rest):
                if len(key) < MIN_BASE_STEM:
                    continue
                found.update(self.stems.get(key, ()))
        found.discard(word)
        return sorted(found)


def adjective_suffix(lemma: str, pos: str, suffix: str) -> bool:
    """A hard -ský/-cký adjective (``adjective`` evidence)."""
    return (
        pos == "adj" and suffix in ADJECTIVE_SUFFIXES and lemma.endswith("ý")
    )


@dataclass(frozen=True, slots=True)
class SuffixContext:
    """What the lexicon knows beyond one entry (module docstring):
    the lexeme lemmas, the words other lexemes' headword lines name
    (f=, m=, dim=), and the base index."""

    lemmas: frozenset[str] = frozenset()
    head_named: frozenset[str] = frozenset()
    index: BaseIndex = field(default_factory=BaseIndex)


EMPTY_CONTEXT: Final[SuffixContext] = SuffixContext()


def morphological_evidence(
    lemma: str,
    pos: str,
    suffix: str,
    *,
    formof_targets: Iterable[str],
    stratum: str,
    context: SuffixContext,
) -> str:
    """The first of ``formof``, ``headword``, ``adjective``, ``base``
    that holds for *lemma* with its R26 *suffix* ("" if none; module
    docstring). Etymology and category evidence come first and are
    read by :func:`lexpipe.langs.czech.ortho.suffix_evidence`."""
    if not suffix:
        return ""
    targets = [t for t in formof_targets if t and t != lemma]
    if any(target in context.lemmas for target in targets):
        return FORMOF
    if lemma in context.head_named:
        return HEADWORD
    if adjective_suffix(lemma, pos, suffix):
        return ADJECTIVE
    if (
        pos == "noun"
        and stratum not in OPAQUE_STRATA
        and context.index.bases(lemma, suffix)
    ):
        return BASE
    return ""
