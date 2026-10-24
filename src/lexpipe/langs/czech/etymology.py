"""Etymological stratum of Czech bases and derivatives (design rule R22).

Every base and derivative gets a stratum ∈ {``inherited``, ``borrowed
whole``, ``internal formation``, ``unknown``} (:data:`STRATA`, the labels
of every language, so cross-language tables join). Three sources, in
precedence order (:func:`combine`):

1. **English Wiktionary** etymology templates (CC BY-SA; the only source
   allowed into tracked tables, R9), read with the shared parser and
   classifier of :mod:`lexpipe.etymology` (the same code as Italian)
   and the Czech profile :data:`lexpipe.etymology.CZECH`: ``der``/
   ``uder`` from Old Czech, Proto-West Slavic, Proto-Slavic or older is
   inheritance; a term without a language prefix is Czech; "Inherited
   from …" / "From <Czech word> …" sentences decide without a template
   (sněžit: "From sníh"). :func:`wiktionary_stratum` is
   :func:`lexpipe.etymology.stratum` with that profile.
2. **czetyl** (Czech Etymological Lexicon 1.0, Rejzek et al. 2025;
   CC BY-NC-SA): ``loan`` → borrowed whole, ``native`` → inherited.
   Homonyms with both verdicts give no verdict.
3. **DeriNet 2.3** (CC BY-NC-SA): ``Loanword=True`` on a node without a
   derivational parent → borrowed whole. DeriNet is otherwise used for
   structure only: a derived node's flag is inherited from its tree's
   root and says nothing about the derivative's own history.

``stratum_wikt`` (source 1 only) feeds tracked tables; ``stratum`` (all
three) feeds the NC files and the NC-labelled report sections.

Adapted loans and hybrids (Czech convention, audit 2026-10-23)
--------------------------------------------------------------
The shared classifier makes an internal kind whose text opens with a
borrowing "borrowed whole" (basis ``internal_analysis``: "From Latin
tragicus + Czech -ický"). That is right for an **adapted loan**, a
foreign word given the Czech suffix that is the equivalent of its own
ending (-icus/-ικός/-ique → -ický; -tiō → -ace/-uce/-ce; a foreign
adjective → -ní/-ný/-ální; -ism → -ismus, -ist → -ista, -in → -ýn;
-ics → -ika), but not for a **hybrid**, a Czech formation on a
borrowed stem or name (fan + -ouš + -ek → fanoušek, Welsh + -an →
Velšan, Hobel + -ovka → hoblovka, mainstream + -ový, Provençal + -ec,
filou + -ta). :func:`refine_stratum` keeps ``borrowed whole`` with the
basis ``adapted_loan`` for the first and makes the second ``internal
formation`` with the basis ``hybrid``; a foreign term that exists as
a Czech headword (mainstream) also makes a hybrid. Steps without a
Czech suffix (compounds: hokejbal, orogeneze) keep the shared verdict.

:func:`affix_parts` reads the (base, suffix) pairs that Czech
affixation steps state, for the derivations stage.
"""

from __future__ import annotations

import re
import unicodedata
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

from lexpipe import etymology as core
from lexpipe.etymology import (
    BORROWED,
    INHERITED,
    INTERNAL,
    STRATA,
    UNKNOWN,
    Step,
    clean_term,
    first_transmission,
)

__all__ = [
    "BORROWED",
    "INHERITED",
    "INTERNAL",
    "STRATA",
    "UNKNOWN",
    "Step",
    "clean_term",
    "first_transmission",
]

PROFILE: Final[core.EtymologyProfile] = core.CZECH
CZECH: Final[str] = PROFILE.code
INHERITED_SOURCES: Final[frozenset[str]] = PROFILE.inherited_sources

SOURCE_WIKTIONARY: Final[str] = "wiktionary"
SOURCE_CZETYL: Final[str] = "czetyl"
SOURCE_DERINET: Final[str] = "derinet"


def template_steps(template: Mapping[str, Any]) -> list[Step]:
    """The steps one etymology template states
    (:func:`lexpipe.etymology.template_steps`)."""
    return core.template_steps(template, CZECH)


def steps(templates: Iterable[Any]) -> list[Step]:
    """The steps of a Czech entry (:func:`lexpipe.etymology.steps`)."""
    return core.steps(templates, CZECH)


def wiktionary_stratum(
    found: Sequence[Step], text: str, *, word: str = ""
) -> tuple[str, str]:
    """``(stratum, basis)`` from Wiktionary
    (:func:`lexpipe.etymology.stratum` with the Czech profile; *word*
    locates the prose of an etymology tree)."""
    return core.stratum(found, text, PROFILE, word=word)


ADAPTED_LOAN: Final[str] = "adapted_loan"
HYBRID: Final[str] = "hybrid"
INTERNAL_ANALYSIS: Final[str] = "internal_analysis"
# Czech suffixes that adapt a foreign word of the same kind whatever
# its ending (learned adjectives and -tiō nouns, -ics nouns).
ADAPTING_SUFFIXES: Final[frozenset[str]] = frozenset(
    {
        "-ický",
        "-ace",
        "-uce",
        "-ce",
        "-ze",
        "-ance",
        "-ence",
        "-ika",
        "-ní",
        "-ný",
        "-ální",
        "-ární",
        "-ivní",
        "-ilní",
        "-ózní",
    }
)
# Czech suffixes that adapt a foreign word only when it already ends in
# the equivalent (virilismus ← virīlis + -ismus is a hybrid).
EQUIVALENT_ENDINGS: Final[Mapping[str, tuple[str, ...]]] = {
    "-ismus": ("ism", "ismus", "isme", "ismo", "ισμός"),
    "-izmus": ("ism", "ismus", "isme", "ismo", "ισμός"),
    "-ista": ("ist", "ista", "iste", "ιστής"),
    "-ita": ("itas", "itās", "ity", "ité", "ität", "ità"),
    "-ýn": ("in", "ine", "ino", "ina"),
    "-ín": ("in", "ine", "ino", "ina"),
}
_LANG_PREFIX: Final[re.Pattern[str]] = re.compile(r"^[a-z]{2,3}(?:-[a-z]+)*:")


def adaptation_key(found: Sequence[Step]) -> str:
    """``suffix|term|…`` of the first transmission step when it is an
    internal kind with a Czech suffix (the last ``-x`` term) and
    non-affix terms; "" otherwise. :func:`refine_stratum` reads it."""
    first = first_transmission(found)
    if first is None or first.kind not in core.INTERNAL_KINDS:
        return ""
    suffixes = [t for t in first.terms if t.startswith("-")]
    terms = [
        _LANG_PREFIX.sub("", t)
        for t in first.terms
        if t and not t.startswith("-") and not t.endswith("-")
    ]
    if not suffixes or not terms:
        return ""
    return "|".join([suffixes[-1], *terms])


def refine_stratum(
    stratum: tuple[str, str], key: str, headwords: frozenset[str]
) -> tuple[str, str]:
    """Adapted loan or hybrid for a ``borrowed whole`` /
    ``internal_analysis`` verdict (module docstring); any other verdict
    is returned unchanged. *key* is :func:`adaptation_key` of the
    entry's steps, *headwords* the words of both editions."""
    if stratum != (BORROWED, INTERNAL_ANALYSIS) or not key:
        return stratum
    suffix, *terms = key.split("|")
    folded = {t.lower() for t in headwords}
    if any(term.lower() in folded for term in terms):
        return INTERNAL, HYBRID
    if suffix in ADAPTING_SUFFIXES:
        return BORROWED, ADAPTED_LOAN
    endings = EQUIVALENT_ENDINGS.get(suffix, ())
    if endings and any(
        core.fold(term).endswith(tuple(core.fold(e) for e in endings))
        for term in terms
    ):
        return BORROWED, ADAPTED_LOAN
    return INTERNAL, HYBRID


def affix_parts(found: Sequence[Step]) -> list[tuple[str, str]]:
    """``(base, suffix)`` pairs stated by Czech affixation steps
    (``:af``/``af``/``suffix``/``surf``): every non-affix term before
    the last ``-suffix`` term. ``{{ety|cs|:af|kniha<t:book>|-ní}}`` →
    ``[("kniha", "-ní")]``."""
    out: list[tuple[str, str]] = []
    for step in found:
        if step.kind not in ("af", "surf") or step.lang not in ("", CZECH):
            continue
        suffixes = [t for t in step.terms if t.startswith("-")]
        if not suffixes:
            continue
        suffix = suffixes[-1]
        for term in step.terms:
            if term.startswith("-") or term.endswith("-") or ":" in term:
                continue
            pair = (term, suffix)
            if pair not in out:
                out.append(pair)
    return out


# ---------------------------------------------------------------------------
# czetyl and the combination
# ---------------------------------------------------------------------------

_HOMONYM_RE: Final[re.Pattern[str]] = re.compile(r"-\d+\Z")


def load_czetyl(path: Path) -> dict[str, str]:
    """``word -> stratum`` from czetyl.tsv (``word``, ``etymology``,
    ``is_loan``). Homonym indices (``a-2``) are dropped; a word whose
    homonyms disagree is left out."""
    verdicts: dict[str, set[str]] = defaultdict(set)
    with path.open("r", encoding="utf-8", newline="") as fh:
        header = fh.readline().rstrip("\r\n").split("\t")
        if header[:3] != ["word", "etymology", "is_loan"]:
            raise ValueError(f"{path}: unexpected header {header!r}")
        for line in fh:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 3 or not parts[0]:
                continue
            word = unicodedata.normalize(
                "NFC", _HOMONYM_RE.sub("", parts[0].strip())
            )
            value = parts[2].strip()
            if value == "loan":
                verdicts[word].add(BORROWED)
            elif value == "native":
                verdicts[word].add(INHERITED)
    return {
        w: next(iter(v)) for w, v in sorted(verdicts.items()) if len(v) == 1
    }


@dataclass(frozen=True, slots=True)
class Stratum:
    """A stratum with the source that decided it."""

    stratum: str
    source: str
    basis: str = ""


def combine(
    wiktionary: tuple[str, str],
    czetyl: str | None,
    derinet_loan_root: bool | None,
) -> Stratum:
    """Precedence Wiktionary > czetyl > DeriNet (module docstring)."""
    stratum, basis = wiktionary
    if stratum != UNKNOWN:
        return Stratum(stratum, SOURCE_WIKTIONARY, basis)
    if czetyl:
        return Stratum(czetyl, SOURCE_CZETYL, "is_loan")
    if derinet_loan_root:
        return Stratum(BORROWED, SOURCE_DERINET, "Loanword")
    return Stratum(UNKNOWN, "", "")
