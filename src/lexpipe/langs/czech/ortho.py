"""Czech stems, palatalization targets, mutations and P/N endings (design §6).

Everything here is a pure function of spellings, so the qc and
features stages and the tests share one definition.

Stems
-----
The stem is the first GEN.SG form minus its ending
(:func:`genitive_stem`): the oblique stem, which undoes the fleeting
e (Turek: Turka → Turk-) and the quantity alternation (dům: domu →
dom-). Czech writes ň ť ď as n t d before ě (ohně = ohň-e), so n/t/d
+ ě counts as ň/ť/ď + e there, and a lemma in ň/ť/ď (loď: lodi) gives
its palatal back to the stem. The longest common prefix of the first
forms of the oblique cells whose endings never palatalize
(:data:`STEM_CELLS`: GEN.SG, DAT.PL, INS.PL, INS.SG; compared with
vowel length folded and the palatal spelling) is a consistency check:
``stem_method`` is ``genitive`` when it agrees and ``genitive_check``
when it does not (míra: měrou/mírou, měrám; a bad cs.wiktionary cell,
kabátníci in INS.PL), which the features report counts for manual
review. Without a GEN.SG form the common prefix decides (``lcp``); one
form alone is stripped of its cell's ending (``strip``: the pluralia
tantum fleky, whose only stem cell is DAT.PL flekům → flek-); with no
oblique form the lemma's final vowel is removed (``lemma``). A long
í/ý that every ending shares (nádraží, nádražím) is not part of the
stem, and a stem that contains the whole -í lemma (Polabí: INS.SG
Polabím alone) is cut back to the lemma minus that vowel.

The declension type is read off the GEN.SG ending (:func:`declension`):
-a/-u/-y hard, -e/-ě/-í soft, -i i-stem, -ého/-ího adjectival, and -é
of an -á noun adjectival too (vysoká, pokojská, dovolená). Without a
GEN.SG form (cs.wiktionary tables that leave it out, pluralia tantum)
an INS.SG in -í, or a lemma that is the stem plus -e/-ě/-í
(chovankyně, lázně, Polabí), marks the soft declension.

Targets (:func:`target_info`)
-----------------------------
The stem-final unit (``ch`` is one letter) gives the target class::

    velar             k g h ch   (clusters sk ck)
    r                 r          (cluster Cr: r after a consonant)
    coronal_inferred  t d n      (clusters st zd); palatalization is
                                 invisible in spelling, so it is
                                 inferred from a following ě/i/í
    c                 c          first palatalization only (VOC otče)
    sibilant          s z        no palatalization in inflection
    labial            p b f v m w  glide before ě (rybě [bjɛ])
    lateral           l          never a target
    soft              ž š č ř j ď ť ň, and any stem of a soft or
                      i-stem declension
    vowel             a vowel-final stem (rádi-o, muze-um)
    other             anything else (x, q; ``-velar``/``c_as_k`` hints)

Core targets are velar and r; coronal_inferred rows are kept apart
(R17/§6: excluded from rule TP and learner evidence).

Mutations (:func:`align`)
-------------------------
A cell form is aligned with the stem as stem + ending, trying the
unchanged stem first and then the orthographic mutations of the stem
final: small (second) palatalization k→c, g/h→z, ch→š, r→ř, sk→sc
(before e, ích) / sk→šť (before i, í), ck→čť; big (first, vocative)
k→č, g/h→ž, c→č. Before i, í and ě the ť of šť/čť is spelled t
(čeští, němečtí), and those spellings are matched there only.

P and N endings (:func:`classify`)
----------------------------------
Per cell and gender class, following the table of design §6:

    NOM.PL m.an       P -i (-í of adjectival nouns)
                      N -ové, -é (allomorph)
    LOC.SG m.in, n    P -e/-ě         N -u (allomorph), -i (soft)
    LOC.SG f, DAT.SG f P -e/-ě        (no N in the hard paradigm)
    LOC.PL m, n       P -ích          N -ách (allomorph), -ech (non-trigger)
    VOC.SG m.an       P -e after a velar (člověče), r after a consonant
                      (bratře) or the c of an -ec noun (otče); N -e
                      elsewhere (soudce, studente,
                      doktore: non-trigger); -u is the regular velar
                      vocative and no opportunity
    INS.SG m, n       N -em (non-trigger front vowel)
    NOM.PL.M.AN adj   P -í

``n_kind`` separates the N allomorphs that compete with a P ending in
that cell (``allomorph``: -ové, -é, -u, -ách) from front-vowel endings
that never trigger (``nontrigger``: -em, -ech, coronal VOC -e; R16's
N = e rows) and soft-paradigm endings (``soft``).
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Final

VOWELS: Final[frozenset[str]] = frozenset("aáäeéěëiíoóöuúůüyý")
FRONT_VOWELS: Final[frozenset[str]] = frozenset("eéěií")
VELARS: Final[frozenset[str]] = frozenset({"k", "g", "h", "ch"})
CORONALS: Final[frozenset[str]] = frozenset({"t", "d", "n"})
LABIALS: Final[frozenset[str]] = frozenset("pbfvmw")
SIBILANTS: Final[frozenset[str]] = frozenset("sz")
SOFT_LETTERS: Final[frozenset[str]] = frozenset("žščřjďťň")

VELAR: Final[str] = "velar"
R_CLASS: Final[str] = "r"
CORONAL: Final[str] = "coronal_inferred"
C_CLASS: Final[str] = "c"
SIBILANT: Final[str] = "sibilant"
LABIAL: Final[str] = "labial"
LATERAL: Final[str] = "lateral"
SOFT: Final[str] = "soft"
VOWEL: Final[str] = "vowel"
OTHER: Final[str] = "other"
NONE: Final[str] = "none"
CORE_TARGETS: Final[frozenset[str]] = frozenset({VELAR, R_CLASS})
MUTATING_TARGETS: Final[frozenset[str]] = frozenset({VELAR, R_CLASS, CORONAL})

STEM_CELLS: Final[tuple[str, ...]] = ("GEN.SG", "DAT.PL", "INS.PL", "INS.SG")
_STRIP_ENDINGS: Final[tuple[str, ...]] = (
    "ého",
    "ího",
    "ete",
    "ěte",
    "ena",
    "ou",
    "a",
    "u",
    "y",
    "e",
    "ě",
    "i",
    "í",
    "o",
)

HARD: Final[str] = "hard"
SOFT_DECL: Final[str] = "soft"
ISTEM: Final[str] = "i-stem"
ADJECTIVAL: Final[str] = "adjectival"
UNKNOWN: Final[str] = "unknown"

P: Final[str] = "P"
N: Final[str] = "N"
N_ALLOMORPH: Final[str] = "allomorph"
N_NONTRIGGER: Final[str] = "nontrigger"
N_SOFT: Final[str] = "soft"

# (original, spelling, label, kind): kind "small" (second
# palatalization, case endings) or "big" (first palatalization:
# vocative, derivation). Before i/í/ě Czech spells ť as t, so the
# clusters šť and čť appear as št and čt there (čeští, němečtí); both
# spellings carry the label sk→šť / ck→čť.
MUTATIONS: Final[tuple[tuple[str, str, str, str], ...]] = (
    ("sk", "sc", "sk→sc", "small"),
    ("sk", "šť", "sk→šť", "small"),
    ("sk", "št", "sk→šť", "small"),
    ("ck", "čť", "ck→čť", "small"),
    ("ck", "čt", "ck→čť", "small"),
    ("ch", "š", "ch→š", "small"),
    ("k", "c", "k→c", "small"),
    ("g", "z", "g→z", "small"),
    ("h", "z", "h→z", "small"),
    ("r", "ř", "r→ř", "small"),
    ("k", "č", "k→č", "big"),
    ("g", "ž", "g→ž", "big"),
    ("h", "ž", "h→ž", "big"),
    ("c", "č", "c→č", "big"),
)
# Spellings that stand for ť only before these letters.
_T_BEFORE: Final[frozenset[str]] = frozenset("iíě")
PALATAL_OUTPUT: Final[Mapping[str, str]] = {
    "k→c": "t͡s",
    "g→z": "z",
    "h→z": "z",
    "ch→š": "ʃ",
    "r→ř": "r̝",
    "sk→sc": "st͡s",
    "sk→šť": "ʃc",
    "ck→čť": "t͡ʃc",
    "k→č": "t͡ʃ",
    "g→ž": "ʒ",
    "h→ž": "ʒ",
    "c→č": "t͡ʃ",
}
INFERRED_OUTPUT: Final[Mapping[str, str]] = {
    "t": "c",
    "d": "ɟ",
    "n": "ɲ",
    "st": "sc",
    "zd": "zɟ",
}
GLIDE_OUTPUT: Final[Mapping[str, str]] = {
    "p": "pj",
    "b": "bj",
    "f": "fj",
    "v": "vj",
    "w": "vj",
    "m": "mɲ",
}
MUTATION_ORTHOGRAPHIC: Final[str] = "orthographic"
MUTATION_INFERRED: Final[str] = "inferred"
MUTATION_GLIDE: Final[str] = "labial_glide"


def final_unit(stem: str) -> str:
    """The stem's last letter, with ``ch`` as one letter (lowercase)."""
    lowered = stem.lower()
    if lowered.endswith("ch"):
        return "ch"
    return lowered[-1:] if lowered else ""


def common_prefix(forms: Sequence[str]) -> str:
    """Longest common prefix of *forms* (character-wise, NFC input)."""
    if not forms:
        return ""
    first = forms[0]
    end = len(first)
    for form in forms[1:]:
        end = min(end, len(form))
        for index in range(end):
            if first[index] != form[index]:
                end = index
                break
    return first[:end]


def strip_ending(form: str) -> str:
    """*form* minus the longest known nominal ending (keeps ≥ 1 letter)."""
    lowered = form.lower()
    for ending in _STRIP_ENDINGS:
        if lowered.endswith(ending) and len(form) > len(ending):
            return form[: -len(ending)]
    return form


# Endings of the stem cells (:data:`STEM_CELLS`), longest first, for a
# stem read off one cell form (the pluralia tantum fleky: DAT.PL
# flekům only). GEN.SG -é is the adjectival feminine (vysoké).
CELL_ENDINGS: Final[Mapping[str, tuple[str, ...]]] = {
    "GEN.SG": (
        "ého",
        "ího",
        "ete",
        "ěte",
        "a",
        "u",
        "y",
        "e",
        "ě",
        "i",
        "í",
        "é",
    ),
    "DAT.PL": ("ům", "ám", "ím", "em", "ým"),
    "INS.PL": ("ami", "emi", "ími", "ými", "mi", "y", "i"),
    "INS.SG": ("em", "ěm", "ím", "ým", "ou", "í"),
}
_PALATAL_OF: Final[Mapping[str, str]] = {"n": "ň", "t": "ť", "d": "ď"}
_PLAIN_OF: Final[Mapping[str, str]] = {"ň": "n", "ť": "t", "ď": "d"}


def strip_cell_ending(cell: str, form: str) -> str:
    """*form* minus the longest ending of *cell* (keeps ≥ 1 letter;
    :data:`CELL_ENDINGS`, else :func:`strip_ending`)."""
    endings = CELL_ENDINGS.get(cell)
    if endings is None:
        return strip_ending(form)
    lowered = form.lower()
    for ending in sorted(endings, key=len, reverse=True):
        if lowered.endswith(ending) and len(form) > len(ending):
            return form[: -len(ending)]
    return form


def _usable(forms: Sequence[str]) -> list[str]:
    """Cell forms as single words: a cs.wiktionary cell that holds
    several variants in one string ("pestřce / pestřece") gives them
    one by one; other forms with a space are left out."""
    out: list[str] = []
    for form in forms:
        for part in form.split("/"):
            word = part.strip()
            if word and " " not in word and word not in out:
                out.append(word)
    return out


def genitive_stem(forms: Mapping[str, Sequence[str]], lemma: str) -> str:
    """The first GEN.SG form minus its ending (module docstring): n/t/d
    before the ending -ě is ň/ť/ď (ohně → ohň-, koně → koň-), and a
    lemma in ň/ť/ď gives its palatal back to the stem (loď: lodi →
    loď-); "" without a usable GEN.SG form or a known ending."""
    genitives = _usable(forms.get("GEN.SG") or ())
    if not genitives:
        return ""
    form = genitives[0]
    stem = strip_cell_ending("GEN.SG", form)
    if stem == form:
        return ""
    ending = form[len(stem) :]
    final = stem[-1:]
    if ending.startswith("ě") and final.lower() in _PALATAL_OF:
        palatal = _PALATAL_OF[final.lower()]
        stem = stem[:-1] + (palatal.upper() if final.isupper() else palatal)
    lemma_final = lemma[-1:].lower()
    if (
        lemma_final in _PLAIN_OF
        and stem[-1:].lower() == _PLAIN_OF[lemma_final]
    ):
        stem = stem[:-1] + lemma_final  # loď: lodi → loď-, zeď: zdi → zď-
    return stem


def _palatal_agrees(common: str, genitive: str) -> bool:
    """Does the common prefix agree with a GEN.SG stem in ň/ť/ď? The
    prefix stops before the palatal (oheň: ohně, ohňům → oh-) or keeps
    its plain spelling (kůň: koně, koním → kon-)."""
    final = genitive[-1:].lower()
    if final not in _PLAIN_OF:
        return False
    head = fold_length(genitive[:-1])
    folded = fold_length(common)
    return folded in (head, head + _PLAIN_OF[final])


_LENGTH_FOLD: Final[dict[int, int]] = str.maketrans("áéíóúůý", "aeiouoy")


def fold_length(text: str) -> str:
    """*text* lowercased with vowel length removed (á→a … ů→o).

    One character in, one out, so prefix lengths carry over to the
    unfolded string. Stems alternate in quantity (dráha: dráhy ~
    drahám, skála ~ skalám); compared unfolded, their common prefix
    would stop at the vowel (dr-, sk-) and misplace the stem final.
    """
    return text.lower().translate(_LENGTH_FOLD)


def stem_length(
    forms: Mapping[str, Sequence[str]], firsts: Sequence[str]
) -> int:
    """Length of the common stem of *firsts* (the first forms of the
    stem cells).

    Compared with vowel length folded, a stem that alternates in
    quantity keeps its final consonant (dráhy ~ drahám → drah); but
    folding must not run into the endings, where -ám and -ami differ
    only in length (Franky: Frankám, Frankami). The DAT.PL/INS.PL
    endings (-ům/-y, -ám/-ami, -ím/-ími …) differ in their first
    letter as written, so the folded prefix is capped by their
    unfolded common prefix; without that pair the unfolded prefix of
    *firsts* is used, as before.
    """
    plural = [
        found[0] for cell in ("DAT.PL", "INS.PL") if (found := forms.get(cell))
    ]
    if len(plural) == 2 and plural[0] != plural[1]:
        folded = len(common_prefix([fold_length(f) for f in firsts]))
        return min(folded, len(common_prefix(plural)))
    return len(common_prefix(firsts))


def common_stem(forms: Mapping[str, Sequence[str]]) -> str:
    """The common prefix of the stem cells' first usable forms (cut
    from the first with the length of :func:`stem_length`, so dráha
    gives dráh-); "" with fewer than two distinct forms."""
    usable: dict[str, list[str]] = {
        cell: found
        for cell in STEM_CELLS
        if (found := _usable(forms.get(cell) or ()))
    }
    firsts = [found[0] for found in usable.values()]
    if len(set(firsts)) < 2:
        return ""
    return firsts[0][: stem_length(usable, firsts)]


def noun_stem(
    forms: Mapping[str, Sequence[str]], lemma: str
) -> tuple[str, str]:
    """``(stem, method)`` of a noun from its cells (module docstring)."""
    genitive = genitive_stem(forms, lemma)
    common = common_stem(forms)
    firsts = [
        (cell, found[0])
        for cell in STEM_CELLS
        if (found := _usable(forms.get(cell) or ()))
    ]
    if len(common) > 1 and common[-1] in "íý":
        common = common[:-1]  # nádraží, nádražím: the shared í
    if genitive:
        agrees = (
            not common
            or fold_length(common) == fold_length(genitive)
            or _palatal_agrees(common, genitive)
        )
        stem, method = genitive, "genitive" if agrees else "genitive_check"
    elif common:
        stem, method = common, "lcp"
    elif firsts:
        cell, form = firsts[0]
        stem, method = strip_cell_ending(cell, form), "strip"
    else:
        stem, method = strip_ending(lemma), "lemma"
    if not stem:
        stem, method = strip_ending(lemma), "lemma"
    if _runs_past_lemma(stem, lemma):
        stem = lemma[:-1]
    if len(stem) > 1 and stem[-1] in "íý":
        stem = stem[:-1]
    return stem, method


def _runs_past_lemma(stem: str, lemma: str) -> bool:
    """A stem that contains the whole -í lemma has run into an ending
    (Polabí: INS.SG Polabím alone, or Krušnohořím = DAT.PL). Other
    vowel-final lemmas can have longer stems (drama: dramat-)."""
    return (
        len(lemma) > 1
        and lemma.endswith("í")
        and stem.lower().startswith(lemma.lower())
    )


SOFT_LEMMA_ENDINGS: Final[frozenset[str]] = frozenset({"e", "ě", "í"})


def declension(
    forms: Mapping[str, Sequence[str]], stem: str, lemma: str = ""
) -> str:
    """Declension type from the GEN.SG ending, else from the INS.SG or
    the lemma (module docstring)."""
    genitives = _usable(forms.get("GEN.SG") or ())
    if not genitives:
        return _declension_without_genitive(forms, stem, lemma)
    genitive = genitives[0]
    plain_stem = stem
    if stem[-1:] in _PLAIN_OF:  # ohň- ~ ohně, loď- ~ lodi
        plain_stem = stem[:-1] + _PLAIN_OF[stem[-1]]
    if not (genitive.startswith(stem) or genitive.startswith(plain_stem)):
        return UNKNOWN
    ending = genitive[len(stem) :]
    if ending in ("ého", "ího"):
        return ADJECTIVAL
    if ending == "é" and lemma.endswith("á"):
        return ADJECTIVAL  # adjectival feminine: vysoká, vysoké
    if ending in ("a", "u", "y"):
        return HARD
    if ending in ("e", "ě", "í", "ete", "ěte"):
        return SOFT_DECL
    if ending == "i":
        return ISTEM
    return UNKNOWN


def _declension_without_genitive(
    forms: Mapping[str, Sequence[str]], stem: str, lemma: str
) -> str:
    """Soft when the INS.SG ends in -í (růží, kostí) or the lemma is the
    stem plus -e/-ě/-í (chovankyně, lázně, Polabí); else unknown."""
    if not stem:
        return UNKNOWN
    instrumentals = forms.get("INS.SG") or ()
    if instrumentals and instrumentals[0] == stem + "í":
        return SOFT_DECL
    if lemma.startswith(stem) and lemma[len(stem) :] in SOFT_LEMMA_ENDINGS:
        return SOFT_DECL
    return UNKNOWN


def adjective_stem(lemma: str) -> tuple[str, str]:
    """``(stem, type)``: hard -ý, soft -í, possessive -ův/-in, other."""
    if lemma.endswith("ý") and len(lemma) > 1:
        return lemma[:-1], HARD
    if lemma.endswith("í") and len(lemma) > 1:
        return lemma[:-1], SOFT_DECL
    if lemma.endswith("ův") and len(lemma) > 2:
        return lemma[:-2] + "ov", "possessive"
    if lemma.endswith("in"):
        return lemma, "possessive"
    return lemma, UNKNOWN


@dataclass(frozen=True, slots=True)
class TargetInfo:
    """The stem-final unit, its cluster and its target class.

    ``ec_type`` marks the -ec nouns (otec: stem otc-, VOC otče), the
    only c-stems whose vocative palatalizes; the soudce type keeps
    VOC = NOM.
    """

    stem_final: str
    cluster: str
    target_class: str
    ec_type: bool = False

    @property
    def core(self) -> bool:
        return self.target_class in CORE_TARGETS

    @property
    def coronal_inferred(self) -> bool:
        return self.target_class == CORONAL


def _cluster(stem: str, unit: str) -> str:
    lowered = stem.lower()
    for cluster in ("sk", "ck", "st", "zd"):
        if lowered.endswith(cluster):
            return cluster
    if unit == "r" and len(lowered) > 1 and lowered[-2] not in VOWELS:
        return "Cr"
    return ""


def target_info(
    stem: str,
    *,
    declension_type: str = HARD,
    flags: Sequence[str] = (),
) -> TargetInfo:
    """Classify the stem final (see the module docstring)."""
    unit = final_unit(stem)
    cluster = _cluster(stem, unit)
    if not unit:
        return TargetInfo("", "", NONE)
    if "-velar" in flags or "c_as_k" in flags:
        return TargetInfo(unit, cluster, OTHER)
    if unit in VOWELS:
        target = VOWEL
    elif unit in VELARS:
        target = VELAR
    elif unit == "r":
        target = R_CLASS
    elif unit in CORONALS:
        target = CORONAL
    elif unit == "c":
        target = C_CLASS
    elif unit in SIBILANTS:
        target = SIBILANT
    elif unit in LABIALS:
        target = LABIAL
    elif unit == "l":
        target = LATERAL
    elif unit in SOFT_LETTERS:
        target = SOFT
    else:
        target = OTHER
    if declension_type in (SOFT_DECL, ISTEM) and target in (
        VELAR,
        R_CLASS,
        CORONAL,
        LABIAL,
    ):
        target = SOFT
    return TargetInfo(unit, cluster, target)


@dataclass(frozen=True, slots=True)
class Alignment:
    """A cell form split into (possibly mutated) stem and ending."""

    aligned: bool
    ending: str = ""
    orth_change: str = ""
    mutation_size: str = ""

    @property
    def mutated(self) -> bool:
        return bool(self.orth_change)


def _palatal_aligned(stem: str, form: str) -> bool:
    """A stem in ň/ť/ď against a form that spells it n/t/d before
    ě/i/í (ohň-: ohně, ohni) or as written (ohňům)."""
    final = stem[-1:].lower()
    if final not in _PLAIN_OF or len(form) <= len(stem):
        return False
    head = stem[:-1]
    if not fold_length(form).startswith(fold_length(head)):
        return False
    letter = form[len(head)].lower()
    after = form[len(stem) : len(stem) + 1].lower()
    return letter == final or (letter == _PLAIN_OF[final] and after in "ěií")


def align(stem: str, form: str) -> Alignment:
    """Split *form* as stem + ending (see the module docstring).

    Stem and form are compared as written, then with vowel length
    folded (:func:`fold_length`: dráh- aligns with drahách).
    """
    if _palatal_aligned(stem, form):
        return Alignment(True, form[len(stem) :])  # ohně, ohni ~ ohň-
    lowered_form = form.lower()
    lowered_stem = stem.lower()
    folded_form = fold_length(form)
    if lowered_stem and (
        lowered_form.startswith(lowered_stem)
        or folded_form.startswith(fold_length(stem))
    ):
        return Alignment(True, form[len(stem) :])
    for original, spelling, label, size in MUTATIONS:
        if not lowered_stem.endswith(original):
            continue
        mutated = lowered_stem[: -len(original)] + spelling
        if not (
            lowered_form.startswith(mutated)
            or folded_form.startswith(fold_length(mutated))
        ):
            continue
        ending = form[len(mutated) :]
        if spelling.endswith("t") and ending[:1].lower() not in _T_BEFORE:
            continue  # št/čt is šť/čť only before i, í, ě
        return Alignment(True, ending, label, size)
    return Alignment(False)


def gender_class(readings: Sequence[str], table_animacy: str) -> str:
    """``m.an`` ``m.in`` ``m`` ``f`` ``n`` or ``mixed`` for a noun row.

    A masculine table's animacy wins; otherwise the lexeme's readings
    decide (one gender needed).
    """
    if table_animacy in ("an", "in"):
        return f"m.{table_animacy}"
    distinct = list(dict.fromkeys(readings))
    if len(distinct) == 1:
        return distinct[0]
    genders = {label.split(".")[0] for label in distinct}
    if len(genders) == 1:
        return "m" if genders == {"m"} else distinct[0]
    return "mixed"


@dataclass(frozen=True, slots=True)
class Classification:
    """P/N verdict of one ending in one cell."""

    allomorph_class: str
    n_kind: str = ""

    @property
    def opportunity(self) -> bool:
        return self.allomorph_class in (P, N)


NO_OPPORTUNITY: Final[Classification] = Classification(NONE)
_MASCULINE: Final[frozenset[str]] = frozenset({"m.an", "m.in", "m"})


def _voc_class(ending: str, target: TargetInfo) -> Classification:
    if ending not in ("e", "ě"):
        return NO_OPPORTUNITY
    if target.target_class == VELAR:
        return Classification(P)
    if target.target_class == C_CLASS and target.ec_type:
        return Classification(P)
    if target.target_class == R_CLASS and target.cluster == "Cr":
        return Classification(P)
    return Classification(N, N_NONTRIGGER)


def classify(
    cell: str, gender: str, ending: str, target: TargetInfo
) -> Classification:
    """P/N class of *ending* in *cell* (table in the module docstring)."""
    if gender == "adj":
        if cell == "NOM.PL.M.AN" and ending == "í":
            return Classification(P)
        return NO_OPPORTUNITY
    if cell == "NOM.PL" and gender in ("m.an", "m"):
        if ending in ("i", "í"):  # -í: adjectival nouns (hajný → hajní)
            return Classification(P)
        if ending in ("ové", "é"):
            return Classification(N, N_ALLOMORPH)
        return NO_OPPORTUNITY
    if cell == "LOC.SG":
        if ending in ("e", "ě"):
            return Classification(P)
        if gender in ("m.in", "n", "m"):
            if ending == "u":
                return Classification(N, N_ALLOMORPH)
            if ending == "i":
                return Classification(N, N_SOFT)
        return NO_OPPORTUNITY
    if cell == "DAT.SG" and gender == "f":
        return Classification(P) if ending in ("e", "ě") else NO_OPPORTUNITY
    if cell == "LOC.PL" and (gender in _MASCULINE or gender == "n"):
        if ending == "ích":
            return Classification(P)
        if ending == "ách":
            return Classification(N, N_ALLOMORPH)
        if ending == "ech":
            return Classification(N, N_NONTRIGGER)
        return NO_OPPORTUNITY
    if cell == "VOC.SG" and gender in ("m.an", "m"):
        return _voc_class(ending, target)
    if cell == "INS.SG" and (gender in _MASCULINE or gender == "n"):
        return (
            Classification(N, N_NONTRIGGER)
            if ending == "em"
            else NO_OPPORTUNITY
        )
    return NO_OPPORTUNITY


@dataclass(frozen=True, slots=True)
class Mutation:
    """What happened to the stem final in one cell form."""

    mutation: bool | None
    kind: str = ""
    palatal_output: str = ""


def mutation_of(
    alignment: Alignment, target: TargetInfo, allomorph_class: str
) -> Mutation:
    """Orthographic, inferred (t/d/n + ě/i/í) or glide change, if any."""
    if not alignment.aligned:
        return Mutation(None)
    if alignment.mutated:
        return Mutation(
            True,
            MUTATION_ORTHOGRAPHIC,
            PALATAL_OUTPUT.get(alignment.orth_change, ""),
        )
    first = alignment.ending[:1]
    if target.target_class == CORONAL and allomorph_class == P:
        if first in ("ě", "i", "í"):
            key = (
                target.cluster
                if target.cluster in ("st", "zd")
                else target.stem_final
            )
            return Mutation(
                True, MUTATION_INFERRED, INFERRED_OUTPUT.get(key, "")
            )
        return Mutation(False)
    if target.target_class == LABIAL and first == "ě":
        return Mutation(
            False, MUTATION_GLIDE, GLIDE_OUTPUT.get(target.stem_final, "")
        )
    return Mutation(False)


UNDERGOER: Final[str] = "undergoer"
UNEXPLAINED: Final[str] = "unexplained"
ALLOMORPH_LABEL: Final[str] = "non exception (allomorph)"
NON_EXCEPTION: Final[str] = "non exception"


def exception_label(
    *,
    in_domain: bool,
    cell: str,
    classification: Classification,
    target: TargetInfo,
    mutation: bool | None,
    nde_class: str = "",
) -> str:
    """Row label: undergoer, unexplained, non exception (allomorph),
    nde:<class> or non exception (design §6 and cs_phon (6))."""
    if not in_domain:
        return NON_EXCEPTION
    mutable = target.target_class in MUTATING_TARGETS or (
        target.target_class == C_CLASS and cell == "VOC.SG"
    )
    if classification.allomorph_class == P and mutable:
        return UNDERGOER if mutation else UNEXPLAINED
    if (
        classification.allomorph_class == N
        and classification.n_kind == N_ALLOMORPH
        and target.target_class in MUTATING_TARGETS
    ):
        return ALLOMORPH_LABEL
    if nde_class:
        return f"nde:{nde_class}"
    return NON_EXCEPTION


# ---------------------------------------------------------------------------
# Suffixes (R26) and non-derived environments
# ---------------------------------------------------------------------------

R26_SUFFIXES: Final[tuple[str, ...]] = (
    "ička",
    "ník",
    "ský",
    "cký",
    "graf",
    "log",
    "urg",
    "ík",
    "ák",
    "ek",
    "ka",
)
# Suffixoids of whole-word loans, accepted on spelling alone (R27).
ORTHOGRAPHIC_SUFFIXES: Final[frozenset[str]] = frozenset(
    {"log", "urg", "graf"}
)
# The target consonant each suffix contains (None: no target).
SUFFIX_TARGET: Final[Mapping[str, str | None]] = {
    "ička": "k",
    "ník": "k",
    "ský": "k",
    "cký": "k",
    "graf": None,
    "log": "g",
    "urg": "g",
    "ík": "k",
    "ák": "k",
    "ek": "k",
    "ka": "k",
}
_SUFFIX_TEMPLATES: Final[frozenset[str]] = frozenset(
    {"af", "affix", "suf", "suffix", "surf", "ety", "compound", "con"}
)


def lemma_suffix(lemma: str) -> str:
    """The longest R26 suffix the lemma ends in ("" if none)."""
    lowered = lemma.lower()
    for suffix in R26_SUFFIXES:
        if lowered.endswith(suffix) and len(lowered) > len(suffix) + 1:
            return suffix
    return ""


def suffix_evidence(
    suffix: str,
    etymology_templates: Sequence[Mapping[str, object]],
    categories: Sequence[str],
) -> str:
    """``etymology``, ``category``, ``orthographic`` or "" (R26)."""
    if not suffix:
        return ""
    wanted = "-" + suffix
    for template in etymology_templates:
        if template.get("name") not in _SUFFIX_TEMPLATES:
            continue
        args = template.get("args")
        if not isinstance(args, Mapping):
            continue
        for value in args.values():
            text = re.sub(r"<[^<>]*>", "", str(value)).strip()
            if text == wanted or text.endswith(":" + wanted):
                return "etymology"
    if f"Czech terms suffixed with {wanted}" in categories:
        return "category"
    if suffix in ORTHOGRAPHIC_SUFFIXES:
        return "orthographic"
    return ""


def target_is_suffix(suffix: str, evidence: str, target: TargetInfo) -> bool:
    """Does the stem-final target lie inside a supported R26 suffix?"""
    consonant = SUFFIX_TARGET.get(suffix)
    return bool(evidence and consonant and target.stem_final == consonant)


_IPA_HARD_CORONAL: Final[re.Pattern[str]] = re.compile(r"[tdn](?=[ɪi])")
# Root-internal NDEB sequences (:func:`ndeb_tokens`): y/ý is [ɪ iː].
_DORSAL_FRONT: Final[re.Pattern[str]] = re.compile(r"(ch|k|g|h)([eéěiíyý])")
_IPA_HARD_CORONAL_TOKEN: Final[re.Pattern[str]] = re.compile(
    r"[tdn](?=([ɪi]))"
)
_IPA_TAIL: Final[str] = "ː.ˈˌ "


def nde_class(stem: str, target: TargetInfo, ipa: Sequence[str]) -> str:
    """NDE class of a vowel-final stem: the velar or coronal before the
    final front vowel belongs to the root (energi-e, rádi-o)."""
    if target.target_class != VOWEL:
        return ""
    lowered = stem.lower()
    core = lowered.rstrip("".join(VOWELS))
    vowel = lowered[len(core) : len(core) + 1]
    if not core or vowel not in FRONT_VOWELS:
        return ""
    unit = final_unit(core)
    if unit in VELARS:
        return "tautomorphemic_dorsal"
    if unit in CORONALS and vowel in ("i", "í"):
        if any(_IPA_HARD_CORONAL.search(text) for text in ipa):
            return "loan_hard_coronal"
    return ""


def ndeb_tokens(lemma: str, ipa: Sequence[str]) -> tuple[str, ...]:
    """Root-internal NDEB sequences: velar + front vowel letters, y/ý
    included (kino: ki; chemie: che; gymnázium: gy) and hard t/d/n +
    [ɪ i] in the IPA (diplomat [ˈdɪplomat]: dɪ). A vowel that ends the
    lemma (or its transcription, length mark aside) is the inflectional
    ending, not part of the root: mladý [ˈmladiː] and vysoký have no
    token."""
    lowered = lemma.lower()
    found = [
        m.group(0)
        for m in _DORSAL_FRONT.finditer(lowered)
        if m.end() < len(lowered)
    ]
    text = (ipa[0] if ipa else "").strip().strip("[]/").rstrip(_IPA_TAIL)
    found.extend(
        m.group(0) + "ɪ"
        for m in _IPA_HARD_CORONAL_TOKEN.finditer(text)
        if m.end(1) < len(text)
    )
    return tuple(found)


# ---------------------------------------------------------------------------
# One call for a lexeme
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class StemAnalysis:
    """Stem, how it was found, declension type and target class."""

    stem: str
    method: str
    declension: str
    target: TargetInfo


def analyse_stem(
    pos: str,
    lemma: str,
    forms: Mapping[str, Sequence[str]],
    flags: Sequence[str] = (),
) -> StemAnalysis:
    """Stem and target of a noun (from its cells) or an adjective (from
    its lemma: hard -ý adjectives are targets; soft -í ones are soft;
    possessives and the rest are ``other``)."""
    if pos == "adj":
        stem, adjective_type = adjective_stem(lemma)
        if adjective_type == HARD:
            target = target_info(stem, flags=flags)
        else:
            target = TargetInfo(
                final_unit(stem),
                "",
                SOFT if adjective_type == SOFT_DECL else OTHER,
            )
        return StemAnalysis(stem, "adjective", adjective_type, target)
    stem, method = noun_stem(forms, lemma)
    decl = declension(forms, stem, lemma)
    target = target_info(stem, declension_type=decl, flags=flags)
    if target.target_class == C_CLASS and is_ec_noun(lemma, stem):
        target = TargetInfo(target.stem_final, target.cluster, C_CLASS, True)
    return StemAnalysis(stem, method, decl, target)


def is_ec_noun(lemma: str, stem: str) -> bool:
    """otec/otc-, chlapec/chlapc-: -ec with a fleeting e."""
    return (
        lemma.lower().endswith("ec")
        and stem.lower().endswith("c")
        and len(stem) == len(lemma) - 1
    )


# Cells with a P/N choice per gender class (the §6 table).
OPPORTUNITY_CELLS: Final[Mapping[str, frozenset[str]]] = {
    "adj": frozenset({"NOM.PL.M.AN"}),
    "m.an": frozenset({"NOM.PL", "VOC.SG", "LOC.PL"}),
    "m.in": frozenset({"LOC.SG", "LOC.PL"}),
    "m": frozenset({"NOM.PL", "VOC.SG", "LOC.SG", "LOC.PL"}),
    "n": frozenset({"LOC.SG", "LOC.PL"}),
    "f": frozenset({"DAT.SG", "LOC.SG"}),
}


def opportunity_cells(pos: str, readings: Sequence[str]) -> frozenset[str]:
    """Cells where the lexeme's readings offer a P/N choice."""
    if pos == "adj":
        return OPPORTUNITY_CELLS["adj"]
    found: set[str] = set()
    for label in readings or ("m",):
        found |= OPPORTUNITY_CELLS.get(label, OPPORTUNITY_CELLS["m"])
    return frozenset(found)
