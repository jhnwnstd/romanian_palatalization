"""Czech derivational suffixes: inventory, decomposition, candidates.

Pure functions of spellings, shared by the ``derivations`` stage and the
analyses that read its pairs (design §6; rules R14, R20–R23, R26).

Inventory (:data:`SUFFIXES`)
----------------------------
Each :class:`Suffix` has a label, a class and the target classes it
palatalizes (``targets``; Połomská 2025, IJP "Změny hlásek při tvoření
slov", Cavirani & Vanden Wyngaerd 2026)::

    P   -it, -et/-ět        verbs      velar, r, c (big pattern)
        -ní/-ný             adjectives velar, c   (r and t/d/n stay)
        -ský/-cký           adjectives velar      (h→ž, ch+s→s, k+s→c)
        -ina                nouns      velar, c
        -í                  adjectives velar, c   (ptačí, ovčí)
        -ek, -ka, -ík       nouns      velar, c, r after a consonant
        -an                 nouns      velar, r, c
        -ník, -nice         nouns      velar, c   (ručník, knižnice)
    N   -ovat               verbs      never
        -ový/-ovní          adjectives never
        -yně                nouns      never
        -ův                 adjectives never
    NT  -ista, -ismus       nouns      never: front vowel, no trigger
                                       (tankista, šachista; R22)

The pairs of design R20 (spell-out) and R23 (co-variation) are
``ni_ovy`` (-ní/-ný vs -ový/-ovní), ``it_ovat`` (-it vs -ovat) and
``ka_yne`` (-ka vs -yně).

Decomposition (:func:`decompose`)
---------------------------------
A derivative is split as optional prefix + base stem (plain or with a
mutated final) + an ending of the suffix. Stems are compared after
:func:`fold` (lowercase, vowel length removed, ``ou`` → ``u``), which
undoes the quantity alternations of word formation (kluk → klouček,
kniha → knížka, roh → růžek). The base is tried with its oblique stem
(features stage: Turek → Turk-) and with its lemma minus a final vowel
(Němec → němec-ký). -ský merges with a velar: h/g → ž + ský (pražský),
ch + ský → ský (Čech → český), k + ský → cký (otrok → otrocký). The
outcome is ``pal`` (the final is mutated), ``plain`` (it is not), and the
suffix-specific verdict ``exception`` (a P suffix that palatalizes the
target but did not: barokní) or ``nontarget`` (the suffix does not
palatalize that target, or the target is not a palatalization target).
Coronal t/d/n changes before i/ě are invisible in spelling, so coronal
bases get ``nontarget`` unless the spelling changed (R16 coronal rule).
The non-trigger suffixes -ista and -ismus decompose only against the
plain stem: katolicismus, mysticismus, stoicismus are the learned
suffix -icismus (Latin -icismus) of a learned -ik noun, not katolík +
-ismus with k → c.

Candidates (R21)
----------------
:func:`candidates` generates, for a base and a suffix, the derivative in
both spellings, plain and palatalized (vzduchový/vzdušný,
parkový/parčový), with the long-vowel variants of diminutives; the
derivations stage looks them up among Wiktionary headwords and in the
corpus. Prefixes are never generated. :func:`paradigm` gives a
derivative's inflected forms by suffix class, for corpus counts of
candidates without a Wiktionary paradigm.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from typing import Final

P_CLASS: Final[str] = "P"
N_CLASS: Final[str] = "N"
NT_CLASS: Final[str] = "NT"

PAL: Final[str] = "pal"
PLAIN: Final[str] = "plain"
EXCEPTION: Final[str] = "exception"
NONTARGET: Final[str] = "nontarget"
OUTCOMES: Final[tuple[str, ...]] = (PAL, PLAIN, EXCEPTION, NONTARGET)

VELAR: Final[str] = "velar"
R_CLASS: Final[str] = "r"
C_CLASS: Final[str] = "c"
CR_CLUSTER: Final[str] = "Cr"
# Target classes whose palatalization is visible in spelling.
VISIBLE_TARGETS: Final[frozenset[str]] = frozenset({VELAR, R_CLASS, C_CLASS})

PAIR_NI_OVY: Final[str] = "ni_ovy"
PAIR_IT_OVAT: Final[str] = "it_ovat"
PAIR_KA_YNE: Final[str] = "ka_yne"
PAIRS: Final[tuple[str, ...]] = (PAIR_NI_OVY, PAIR_IT_OVAT, PAIR_KA_YNE)
PAIR_LABELS: Final[Mapping[str, str]] = {
    PAIR_NI_OVY: "-ní/-ný vs -ový/-ovní",
    PAIR_IT_OVAT: "-it vs -ovat",
    PAIR_KA_YNE: "-ka vs -yně",
}

_VELAR_C: Final[frozenset[str]] = frozenset({VELAR, C_CLASS})
_VELAR_R_C: Final[frozenset[str]] = frozenset({VELAR, R_CLASS, C_CLASS})


@dataclass(frozen=True, slots=True)
class Suffix:
    """One derivational suffix (see the module docstring)."""

    label: str
    klass: str
    pos: str
    endings: tuple[str, ...]
    targets: frozenset[str] = frozenset()
    pair: str = ""
    r_after_consonant_only: bool = False
    diminutive: bool = False

    @property
    def palatalizing(self) -> bool:
        return self.klass == P_CLASS

    def palatalizes(self, target_class: str, cluster: str = "") -> bool:
        """Does this suffix palatalize a stem of *target_class*?"""
        if not self.palatalizing or target_class not in self.targets:
            return False
        if target_class == R_CLASS and self.r_after_consonant_only:
            return cluster == CR_CLUSTER
        return True


SUFFIXES: Final[tuple[Suffix, ...]] = (
    Suffix("-it", P_CLASS, "verb", ("it",), _VELAR_R_C, PAIR_IT_OVAT),
    Suffix("-et", P_CLASS, "verb", ("et", "ět"), _VELAR_R_C),
    Suffix("-ní/-ný", P_CLASS, "adj", ("ní", "ný"), _VELAR_C, PAIR_NI_OVY),
    Suffix("-ský/-cký", P_CLASS, "adj", ("ský", "cký"), frozenset({VELAR})),
    Suffix("-ina", P_CLASS, "noun", ("ina",), _VELAR_C),
    Suffix("-í", P_CLASS, "adj", ("í",), _VELAR_C),
    Suffix(
        "-ek",
        P_CLASS,
        "noun",
        ("ek",),
        _VELAR_R_C,
        r_after_consonant_only=True,
        diminutive=True,
    ),
    Suffix(
        "-ka",
        P_CLASS,
        "noun",
        ("ička", "ka"),
        _VELAR_R_C,
        PAIR_KA_YNE,
        r_after_consonant_only=True,
        diminutive=True,
    ),
    Suffix(
        "-ík",
        P_CLASS,
        "noun",
        ("íček", "ík"),
        _VELAR_R_C,
        r_after_consonant_only=True,
        diminutive=True,
    ),
    Suffix("-an", P_CLASS, "noun", ("an",), _VELAR_R_C),
    Suffix("-ník", P_CLASS, "noun", ("ník",), _VELAR_C),
    Suffix("-nice", P_CLASS, "noun", ("nice",), _VELAR_C),
    Suffix("-ovat", N_CLASS, "verb", ("ovat",), pair=PAIR_IT_OVAT),
    Suffix("-ový/-ovní", N_CLASS, "adj", ("ovní", "ový"), pair=PAIR_NI_OVY),
    Suffix("-yně", N_CLASS, "noun", ("yně",), pair=PAIR_KA_YNE),
    Suffix("-ův", N_CLASS, "adj", ("ův",)),
    Suffix("-ista", NT_CLASS, "noun", ("ista",)),
    Suffix("-ismus", NT_CLASS, "noun", ("ismus",)),
)
BY_LABEL: Final[Mapping[str, Suffix]] = {s.label: s for s in SUFFIXES}
# Decomposition order: longer endings first, so -ovní wins over -ní,
# -ička over -ka, -ník over -ík and -ista over -a stems.
_ENDING_ORDER: Final[tuple[tuple[str, Suffix], ...]] = tuple(
    sorted(
        ((ending, s) for s in SUFFIXES for ending in s.endings),
        key=lambda item: (-len(item[0]), item[1].label, item[0]),
    )
)

# Mutations of a stem-final unit: first (big) and second (small)
# palatalization outputs, plus the cluster outputs (sk→šť/sc, ck→čť).
# Before i/í/ě the ť of šť/čť is spelled t (čeština, němčina).
MUTATIONS: Final[Mapping[str, tuple[tuple[str, str], ...]]] = {
    "k": (("č", "big"), ("c", "small")),
    "g": (("ž", "big"), ("z", "small")),
    "h": (("ž", "big"), ("z", "small")),
    "ch": (("š", "big"),),
    "r": (("ř", "big"),),
    "c": (("č", "big"),),
    "sk": (("šť", "big"), ("št", "big"), ("sc", "small")),
    "ck": (("čť", "big"), ("čt", "big"), ("č", "big")),
    "t": (("c", "big"), ("ť", "big")),
    "d": (("z", "big"), ("ď", "big")),
    "n": (("ň", "big"),),
    "s": (("š", "big"),),
    "z": (("ž", "big"),),
    "st": (("šť", "big"), ("št", "big")),
    "zd": (("žď", "big"), ("žd", "big")),
}

# Prefixes allowed in front of a linked (never a generated) derivative:
# parasynthetic verbs (s-nož-it) and prefixed adjectives (ná-rož-ní).
PREFIXES: Final[tuple[str, ...]] = (
    "před",
    "pod",
    "nad",
    "roz",
    "pře",
    "při",
    "pro",
    "na",
    "ná",
    "za",
    "zá",
    "do",
    "po",
    "vy",
    "vý",
    "od",
    "ob",
    "vz",
    "se",
    "ze",
    "ve",
    "s",
    "z",
    "u",
    "ú",
    "v",
    "o",
)

_FOLD_TABLE: Final[dict[int, int]] = str.maketrans("áéěíóúůý", "aeeiouoy")
_VOWELS: Final[frozenset[str]] = frozenset("aáeéěiíoóuúůyý")
_LENGTHEN: Final[Mapping[str, tuple[str, ...]]] = {
    "a": ("á",),
    "e": ("é",),
    "i": ("í",),
    "o": ("ů", "ó"),
    "u": ("ou", "ú"),
    "y": ("ý",),
}


def fold(text: str) -> str:
    """Comparison key of a stem: lowercase, no length, ``ou`` → ``u``."""
    return text.lower().translate(_FOLD_TABLE).replace("ou", "u")


def final_unit(stem: str, cluster: str = "") -> str:
    """The unit a mutation replaces: a cluster (sk, ck, st, zd), ``ch``
    or the last letter (lowercase)."""
    lowered = stem.lower()
    if cluster in ("sk", "ck", "st", "zd") and lowered.endswith(cluster):
        return cluster
    if lowered.endswith("ch"):
        return "ch"
    return lowered[-1:]


@dataclass(frozen=True, slots=True)
class Decomposition:
    """A derivative split against its base (see the module docstring)."""

    suffix: Suffix
    ending: str
    prefix: str
    mutated: bool
    orth_change: str
    size: str
    stem_used: str

    def outcome(self, target_class: str, cluster: str = "") -> str:
        """``pal``, ``plain``, ``exception`` or ``nontarget``."""
        if self.mutated:
            return PAL
        if self.suffix.klass != P_CLASS:
            return PLAIN
        if target_class in VISIBLE_TARGETS and self.suffix.palatalizes(
            target_class, cluster
        ):
            return EXCEPTION
        return NONTARGET


def base_stems(stem: str, lemma: str) -> tuple[str, ...]:
    """Stems a derivative may build on: the oblique stem, then the lemma
    minus a final vowel (fleeting e: Němec → němec-ký)."""
    lowered = lemma.lower()
    alternative = lowered
    while alternative and alternative[-1] in _VOWELS:
        alternative = alternative[:-1]
    out: list[str] = []
    for item in (stem.lower(), alternative):
        if item and item not in out:
            out.append(item)
    return tuple(out)


def _sky_match(
    rest: str, ending: str, stem: str, unit: str
) -> tuple[bool, str] | None:
    """-ský/-cký against one stem: ``(mutated, change)`` or None."""
    folded, base = fold(rest), fold(stem)
    head = fold(stem[: -len(unit)]) if unit else base
    if ending == "ský":
        if folded == base:
            return False, ""
        if unit in ("h", "g") and folded == head + "ž":
            return True, f"{unit}→ž"
        if unit == "ch" and folded in (head, head + "š"):
            return True, "ch→š"
        if unit == "sk" and folded == head:
            return False, "sk+s→s"
    if ending == "cký":
        if unit == "k" and folded == head:
            return True, "k→c"
        if unit == "c" and folded == head:
            return False, "c+s→c"
        if unit == "ck" and folded == head:
            return False, "ck+s→c"
    return None


def fold_ei(text: str) -> str:
    """:func:`fold` with e and i merged: the old ie alternation of a
    root vowel (řeka → říční, světlo → svítit)."""
    return fold(text).replace("e", "i")


def _stem_match(
    rest: str, stem: str, unit: str
) -> tuple[bool, str, str] | None:
    """``(mutated, change, size)`` if *rest* is *stem*, plain or with its
    final unit mutated; None otherwise. The comparison is made with
    :func:`fold`, then with :func:`fold_ei`."""
    for key in (fold, fold_ei):
        folded = key(rest)
        if folded == key(stem):
            return False, "", ""
        head = key(stem[: -len(unit)]) if unit else ""
        if not unit or not folded.startswith(head):
            continue
        for output, size in MUTATIONS.get(unit, ()):
            if folded == head + key(output):
                return True, f"{unit}→{output}", size
    return None


def _match_suffix(
    word: str,
    ending: str,
    suffix: Suffix,
    stems: Sequence[str],
    cluster: str,
    allow_prefix: bool,
) -> Decomposition | None:
    rest_all = word[: -len(ending)]
    prefixes = ("",) + (PREFIXES if allow_prefix else ())
    for prefix in prefixes:
        if prefix and not rest_all.startswith(prefix):
            continue
        rest = rest_all[len(prefix) :]
        if len(fold(rest)) < 2:
            continue
        for stem in stems:
            unit = final_unit(stem, cluster)
            if suffix.label == "-ský/-cký":
                sky = _sky_match(rest, ending, stem, unit)
                if sky is not None:
                    mutated, change = sky
                    return Decomposition(
                        suffix, ending, prefix, mutated, change, "", stem
                    )
                continue
            found = _stem_match(rest, stem, unit)
            if found is not None:
                mutated, change, size = found
                if mutated and suffix.klass == NT_CLASS:
                    # -ista/-ismus never palatalize: a mutated stem is
                    # the learned -icismus of an -ik base (katolík →
                    # katolicismus, Latin -icismus), not k + -ismus.
                    continue
                return Decomposition(
                    suffix, ending, prefix, mutated, change, size, stem
                )
    return None


def decompose(
    derivative: str,
    stem: str,
    lemma: str,
    *,
    cluster: str = "",
    allow_prefix: bool = False,
    suffixes: Sequence[Suffix] | None = None,
) -> Decomposition | None:
    """Split *derivative* against a base (see the module docstring).

    *stem* is the base's oblique stem, *lemma* its citation form;
    *suffixes* restricts the inventory (default: all). Returns None
    when no suffix fits.
    """
    word = derivative.lower()
    allowed = None if suffixes is None else {s.label for s in suffixes}
    stems = base_stems(stem, lemma)
    for ending, suffix in _ENDING_ORDER:
        if allowed is not None and suffix.label not in allowed:
            continue
        if not word.endswith(ending) or len(word) <= len(ending) + 1:
            continue
        found = _match_suffix(
            word, ending, suffix, stems, cluster, allow_prefix
        )
        if found is not None:
            return found
    return None


# ---------------------------------------------------------------------------
# Candidates (R21)
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Candidate:
    """A generated derivative: word, suffix, spelling and change
    (``size``: big or small palatalization, "" for -ský merges)."""

    word: str
    suffix: Suffix
    mutated: bool
    orth_change: str
    size: str = ""


def _lengthened(stem: str) -> list[str]:
    """*stem* with its last vowel lengthened (knih → kníh, kluk → klouk)."""
    for index in range(len(stem) - 1, -1, -1):
        vowel = stem[index]
        if vowel in _LENGTHEN:
            return [
                stem[:index] + long + stem[index + 1 :]
                for long in _LENGTHEN[vowel]
            ]
        if vowel in _VOWELS:
            return []
    return []


def _sky_candidates(stem: str, unit: str) -> Iterator[tuple[str, bool, str]]:
    yield stem + "ský", False, ""
    head = stem[: -len(unit)]
    if unit in ("h", "g"):
        yield head + "žský", True, f"{unit}→ž"
    elif unit == "ch":
        yield head + "ský", True, "ch→š"
    elif unit == "k":
        yield head + "cký", True, "k→c"


def candidates(stem: str, cluster: str, suffix: Suffix) -> list[Candidate]:
    """Derivatives of *stem* with *suffix* in both spellings (R21).

    Diminutives (-ek, -ka, -ík) also get the stem with its last vowel
    lengthened (knížka, klouček, růžek). Duplicates are dropped; the
    order is deterministic.
    """
    lowered = stem.lower()
    unit = final_unit(lowered, cluster)
    out: dict[str, Candidate] = {}

    def add(word: str, mutated: bool, change: str, size: str = "") -> None:
        if word not in out:
            out[word] = Candidate(word, suffix, mutated, change, size)

    if suffix.label == "-ský/-cký":
        for word, mutated, change in _sky_candidates(lowered, unit):
            add(word, mutated, change)
        return list(out.values())
    head = lowered[: -len(unit)] if unit else lowered
    heads = [head] + (_lengthened(head) if suffix.diminutive else [])
    for current in heads:
        for ending in suffix.endings:
            add(current + unit + ending, False, "")
            for output, size in MUTATIONS.get(unit, ()):
                if (
                    not _spelling_ok(output, ending)
                    or (
                        unit,
                        output,
                    )
                    in _NOT_GENERATED
                ):
                    continue
                add(current + output + ending, True, f"{unit}→{output}", size)
    return list(out.values())


# Mutations recognized in linked derivatives but never generated: ck→č
# needs a lost vowel (německý → němčina); as a candidate it only finds
# other words (micka → míček, racek → ráček).
_NOT_GENERATED: Final[frozenset[tuple[str, str]]] = frozenset({("ck", "č")})


def _spelling_ok(output: str, ending: str) -> bool:
    """Czech spells ť as t before i/í/ě and nowhere else (čeština,
    but sněžťový is not a spelling)."""
    front = ending[:1] in ("i", "í", "ě")
    if output.endswith("ť"):
        return not front
    if output.endswith("t") and output != "t":
        return front
    return True


# ---------------------------------------------------------------------------
# Generated paradigms (corpus counts of candidates)
# ---------------------------------------------------------------------------

_HARD_ADJ: Final[tuple[str, ...]] = (
    "ý",
    "ého",
    "ému",
    "ém",
    "ým",
    "á",
    "é",
    "ou",
    "ých",
    "ými",
)
_SOFT_ADJ: Final[tuple[str, ...]] = ("í", "ího", "ímu", "ím", "ích", "ími")
# Possessive -ův forms that no other paradigm shares: -ovi is a noun's
# DAT.SG (klukovi), -ově/-ovým/-ových/-ovými are forms of -ový
# adjectives (systémových), so they are not generated.
_POSSESSIVE: Final[tuple[str, ...]] = ("ův", "ova", "ovo", "ovy", "ovu")
_VERB_IT: Final[tuple[str, ...]] = (
    "it",
    "ím",
    "íš",
    "í",
    "íme",
    "íte",
    "il",
    "ila",
    "ilo",
    "ili",
    "ily",
    "en",
    "ena",
    "eno",
    "ení",
)
_VERB_OVAT: Final[tuple[str, ...]] = (
    "ovat",
    "uji",
    "uju",
    "uješ",
    "uje",
    "ujeme",
    "ujete",
    "ují",
    "oval",
    "ovala",
    "ovalo",
    "ovali",
    "ovaly",
    "ován",
    "ována",
    "ování",
)
_FEM_A: Final[tuple[str, ...]] = ("a", "y", "u", "ou", "o", "ám", "ách", "ami")
_MASC_AN: Final[tuple[str, ...]] = (
    "",
    "a",
    "u",
    "em",
    "ovi",
    "ové",
    "ů",
    "ům",
    "y",
    "ech",
)
_SOFT_FEM: Final[tuple[str, ...]] = ("e", "i", "í", "ím", "ích", "emi")


def paradigm(word: str, suffix: Suffix) -> tuple[str, ...]:
    """Inflected forms of a derivative, by its suffix class (lowercase,
    lemma first, duplicates dropped). An approximation for corpus
    counts: forms with a palatalized stem alternant (knížce, kloučci)
    are generated only where the alternation is regular."""
    lowered = word.lower()
    label = suffix.label
    forms: list[str] = [lowered]
    if suffix.pos == "adj":
        if lowered.endswith("ý"):
            forms += [lowered[:-1] + e for e in _HARD_ADJ]
        elif lowered.endswith("ův"):
            forms += [lowered[:-2] + e for e in _POSSESSIVE]
        elif lowered.endswith("í"):
            forms += [lowered[:-1] + e for e in _SOFT_ADJ]
    elif label in ("-it", "-et"):
        vowel = lowered[-2]
        stem = lowered[:-2]
        endings = _VERB_IT if vowel == "i" else _VERB_IT[1:6]
        forms += [stem + e for e in endings]
        if vowel != "i":
            forms += [stem + vowel + e for e in ("l", "la", "lo", "li", "ly")]
    elif label == "-ovat":
        forms += [lowered[:-4] + e for e in _VERB_OVAT]
    elif label in ("-ka", "-ina"):
        stem = lowered[:-1]
        forms += [stem + e for e in _FEM_A]
        if label == "-ka":
            forms += [stem[:-1] + "ce", stem[:-1] + "ek"]
        else:
            forms.append(stem + "ě")
    elif label == "-an":
        forms += [lowered + e for e in _MASC_AN[1:] if e != "ové"]
        forms.append(lowered + "é")
    elif label in ("-ek", "-ík", "-ník"):
        stem = lowered[:-2] + "k" if label == "-ek" else lowered
        forms += [stem + e for e in _MASC_AN[1:]]
        forms += [stem[:-1] + "ci", stem[:-1] + "cích"]
    elif label in ("-nice", "-yně"):
        stem = lowered[:-1]
        forms += [stem + e for e in _SOFT_FEM]
    elif label == "-ista":
        stem = lowered[:-1]
        forms += [stem + e for e in ("y", "ovi", "u", "o", "ou", "é", "i")]
        forms += [stem + e for e in ("ů", "ům", "ech", "ech", "ami")]
    elif label == "-ismus":
        stem = lowered[:-2]
        forms += [stem + e for e in ("u", "em", "y", "ů", "ech")]
    return tuple(dict.fromkeys(forms))
