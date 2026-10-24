"""Orthographic reading of verb cells and derivatives (E1, E2; R15).

A root-final target consonant meets a suffix vowel in a verb cell
(E1) or a derivative (E2). What a learner can read off such a form is
(a) which vowel follows the root, and (b) whether the consonant shows
its palatal reflex. Both are read from the spelling, which is
transparent in Romanian:

    target  palatalized before i/e   plain before i/e   plain before a/o/u
    c       ci, ce                   chi, che           ca, co, cu
    g       gi, ge                   ghi, ghe           ga, go, gu
    t       ți, țe                   ti, te             ta, ...
    d       zi, ze                   di, de             da, ...
    s       și, șe                   si, se             sa, ...
    z       ji, je                   zi, ze             za, ...
    st, sc  ști, ște                 sti, ste / sci     sta, sca, ...
    ct      cți, cțe                 cti, cte           cta, ...

A form is split as ``root part + ending``; :func:`reflex` reads the
end of the root part. The vowel class of the ending is ``i`` or ``e``
(front) or ``back`` (a, ă, â, î, o, u); an ending that starts with a
consonant, or a root part that does not end in the target or its
reflex (picur-a from pic, a stem extension), gives no reading (None).

Verb cells (E1)
---------------
The design's E1 cells are the 2sg present and the 3sg subjunctive.
Endings, longest first: 2sg ``-ezi``/``-ești`` (augment: the root
meets the augment's /e/), ``-i``; 3sg subjunctive ``-eze``/``-ească``
(augment), ``-e``, ``-ă`` (back). The augments are read by the vowel
that follows the root, so a dorsal before -ez- (blochezi: ch + e) is a
non-palatalized dorsal before a front vowel and a dorsal before -esc-
(muncești: c + e) a palatalized one. ``scripts/audit_verb_alternation.py``
(the source of ``data/verb_paradigm_verdict.csv``) instead treats every
augmented cell as TRIGGER-ABSENT; :func:`audit_like_verdict` summarizes
a pair this module's way so the two can be compared.

Derivatives (E2)
----------------
Verbs carry their suffix in ``deriv_suffixes`` (-a, -i, -ui, -iza,
-ifica, -ări, -ăni, -arisi). Adjectives do not, so their suffix is
the longest entry of :data:`ADJECTIVE_SUFFIXES` that leaves a root
part ending in the target or its reflex.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Final

from lexpipe.langs.romanian_compare.lexicon import (
    PALATAL_LETTER,
    VOWELS,
    strip_final_vowel,
)

FRONT: Final[str] = "front"
VOWEL_CLASSES: Final[tuple[str, ...]] = ("i", "e", "back")
TRIGGER_VOWELS: Final[dict[str, tuple[str, ...]]] = {
    "-i": ("i",),
    "-e": ("e",),
    "-i/-e": ("i", "e"),
}
TRIGGERS: Final[tuple[str, ...]] = ("-i", "-e", "-i/-e")

VERB_ENDINGS_2SG: Final[tuple[tuple[str, bool], ...]] = (
    ("ești", True),
    ("ezi", True),
    ("i", False),
)
VERB_ENDINGS_SUBJ3SG: Final[tuple[tuple[str, bool], ...]] = (
    ("ească", True),
    ("eze", True),
    ("e", False),
    ("ă", False),
)
VERB_SUFFIXES: Final[tuple[str, ...]] = (
    "-a",
    "-i",
    "-ui",
    "-iza",
    "-ifica",
    "-ări",
    "-ăni",
    "-arisi",
)
# Competing verbalizers of the spell-out TP (R20): -i vs -a/-ui.
PALATALIZING_VERBALIZERS: Final[tuple[str, ...]] = ("-i",)
NONPALATALIZING_VERBALIZERS: Final[tuple[str, ...]] = ("-a", "-ui")

# Adjectival (and a few nominal) suffixes seen in derived_adj, longest
# match first. Feminine forms are listed because the column holds some.
ADJECTIVE_SUFFIXES: Final[tuple[str, ...]] = tuple(
    sorted(
        {
            "icesc",
            "iceasca",
            "esc",
            "ească",
            "ic",
            "ică",
            "ist",
            "istă",
            "ian",
            "iană",
            "ean",
            "eană",
            "iu",
            "ie",
            "ios",
            "ioasă",
            "os",
            "oasă",
            "ar",
            "ară",
            "an",
            "ană",
            "at",
            "ată",
            "it",
            "ită",
            "ut",
            "ută",
            "uit",
            "el",
            "ea",
            "ei",
            "ior",
            "uț",
            "uleț",
            "eț",
            "areț",
            "aș",
            "al",
            "ală",
            "iv",
            "ivă",
            "ent",
            "ant",
            "ui",
            "uie",
            "nic",
            "nică",
            "tor",
        },
        key=lambda suffix: (-len(suffix), suffix),
    )
)


def vowel_class(letter: str) -> str | None:
    """``i`` / ``e`` / ``back`` for a vowel letter, None otherwise."""
    if letter == "i":
        return "i"
    if letter == "e":
        return "e"
    if letter in VOWELS:
        return "back"
    return None


def reflex(root_part: str, seg: str, cluster: str, vowel: str) -> bool | None:
    """Does *root_part* end in the palatal reflex of the target?

    True for the reflex, False for the plain consonant, None when the
    root part ends in neither. Dorsals spelled c/g are palatal only
    before a front vowel (c + a is plain [k]); ch/gh is always plain.
    """
    text = root_part.lower()
    front = vowel in ("i", "e")
    if seg in ("c", "g"):
        if text.endswith(seg + "h"):
            return False
        if text.endswith(seg):
            return front
        return None
    if seg == "s" and cluster.startswith(("st", "sc")):
        if text.endswith("șt"):
            return True
        if text.endswith(cluster[:2]):
            return False
        return None
    if seg == "t" and cluster.startswith("ct"):
        if text.endswith("cț"):
            return True
        if text.endswith("ct"):
            return False
        return None
    if seg in PALATAL_LETTER:
        if text.endswith(PALATAL_LETTER[seg]):
            return True
        if text.endswith(seg):
            return False
    return None


@dataclass(frozen=True, slots=True)
class Reading:
    """What a form shows about its root-final target."""

    form: str
    vowel: str
    palatal: bool
    augment: bool = False
    suffix: str = ""


def read_verb_cell(
    form: str,
    endings: Sequence[tuple[str, bool]],
    seg: str,
    cluster: str,
) -> Reading | None:
    """Split a verb cell by its ending and read the target (E1)."""
    text = form.strip().lower()
    for ending, augment in endings:
        if not text.endswith(ending) or len(text) <= len(ending):
            continue
        vowel = vowel_class(ending[0])
        if vowel is None:
            continue
        palatal = reflex(text[: -len(ending)], seg, cluster, vowel)
        if palatal is None:
            return None
        return Reading(text, vowel, palatal, augment)
    return None


def first_variant(cell: str, verb: str) -> str:
    """The first ``|``-separated variant that looks like a form of
    *verb* (same first letter); ``&mdash;`` and empty cells give ""."""
    for variant in cell.split("|"):
        text = variant.strip()
        if not text or text.startswith("&"):
            continue
        if verb and text[0].lower() != verb[0].lower():
            continue
        return text
    return ""


def read_derived_verb(
    verb: str, suffix: str, seg: str, cluster: str
) -> Reading | None:
    """Read a derived verb's infinitive with its known suffix (E2)."""
    ending = suffix.lstrip("-")
    text = verb.strip().lower()
    if not ending or not text.endswith(ending) or len(text) <= len(ending):
        return None
    vowel = vowel_class(ending[0])
    if vowel is None:
        return None
    palatal = reflex(text[: -len(ending)], seg, cluster, vowel)
    if palatal is None:
        return None
    return Reading(text, vowel, palatal, suffix=suffix)


def read_derived_adjective(
    adjective: str, seg: str, cluster: str
) -> Reading | None:
    """Find the suffix of a derived adjective and read the target (E2)."""
    text = adjective.strip().lower()
    for suffix in ADJECTIVE_SUFFIXES:
        if not text.endswith(suffix) or len(text) - len(suffix) < 2:
            continue
        vowel = vowel_class(suffix[0])
        if vowel is None:
            continue
        palatal = reflex(text[: -len(suffix)], seg, cluster, vowel)
        if palatal is not None:
            return Reading(text, vowel, palatal, suffix="-" + suffix)
    return None


def participle(verb: str, suffix: str) -> str:
    """Past participle (masc. sg.) of a derived verb, for attestation:
    -a → -at, -ui → -uit, other -i verbs → -it; "" if unknown."""
    text = verb.strip().lower()
    if suffix == "-ui" and text.endswith("ui"):
        return text + "t"
    if text.endswith("a"):
        return text[:-1] + "at"
    if text.endswith("i"):
        return text[:-1] + "it"
    return ""


def triggers_alternation(seg: str, vowel: str) -> bool:
    """Is *vowel* a trigger for the target under the language's own
    rules (dorsals before i/e, coronals and z before i)?"""
    if seg in ("c", "g"):
        return vowel in ("i", "e")
    return vowel == "i"


def audit_like_verdict(seg: str, readings: Sequence[Reading]) -> str:
    """Summarize a verb's cells as the audit's verdict labels.

    ALTERNATES if a trigger cell shows the reflex, NO-ALTERNATION if
    trigger cells exist and none does, TRIGGER-ABSENT if no cell puts
    the target before a trigger, NO-READING when no cell could be read
    (no paradigm, or a stem change: picur-a from pic, strunj-i from
    strung).
    """
    if not readings:
        return "NO-READING"
    triggering = [r for r in readings if triggers_alternation(seg, r.vowel)]
    if any(r.palatal for r in triggering):
        return "ALTERNATES"
    if triggering:
        return "NO-ALTERNATION"
    return "TRIGGER-ABSENT"


def allomorph_for(trigger: str, vowel: str | None) -> str | None:
    """lexpipe.tp allomorph of a realization in a trigger cell.

    ``P`` when the vowel after the root is one of the cell's triggers,
    ``N`` when it is another front vowel (the cell is realized without
    this trigger), ``nonfront`` for a back vowel; None (no
    observation) when the vowel is unknown.
    """
    if vowel is None:
        return None
    if vowel in TRIGGER_VOWELS[trigger]:
        return "P"
    if vowel in ("i", "e"):
        return "N"
    return "nonfront"


def candidate_verbs(
    lemma: str, seg: str, cluster: str
) -> tuple[str, str] | None:
    """(palatalized, plain) -i infinitive candidates of a target root
    (R21: derivative candidates in both spellings): foc → foci / fochi,
    vad → vazi / vadi, prost → proști / prosti. None when the stem does
    not end in the target."""
    stem = strip_final_vowel(lemma.lower())
    if seg in ("c", "g"):
        if not stem.endswith(seg):
            return None
        return stem + "i", stem + "hi"
    if seg == "s" and cluster.startswith(("st", "sc")):
        if not stem.endswith(cluster[:2]):
            return None
        return stem[:-2] + "ști", stem + "i"
    if seg == "t" and cluster.startswith("ct"):
        if not stem.endswith("ct"):
            return None
        return stem[:-1] + "ți", stem + "i"
    if seg in PALATAL_LETTER and stem.endswith(seg):
        return stem[:-1] + PALATAL_LETTER[seg] + "i", stem + "i"
    return None
