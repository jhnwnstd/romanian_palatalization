"""Etymology of Italian lexemes from Wiktionary templates (R22, R26, R27).

The template parser and the R22 stratum classifier are the shared ones
of :mod:`lexpipe.etymology` (the same code as Czech), with the Italian
profile :data:`lexpipe.etymology.ITALIAN`: ``der``/``uder`` from Vulgar
Latin or Old Italian is inheritance, terms without a language prefix
are Italian, the ``it-deverbal`` templates are deverbals. This module
adds what is Italian:

- :func:`steps`: the shared parser, after the Italian compound
  templates ``it-verb-obj`` (fora + sacco) and ``it-verb-verb``, which
  name their parts from argument 1, are read as ``compound``
  (:func:`italian_templates`): internal formations (R22);
- :func:`etymon_key` — the root key of an etymon (R4);
- :func:`deverbal_verb`: the verb a noun is formed from, from the
  deverbal templates, a back-formation from a verb, ``{{ety|:from|X}}``
  with an Italian verb, or the text "Deverbal from X" / "(short) past
  participle of X"; ``deverbal_from_are`` is true when that verb is a
  first-conjugation verb (-are, -arsi).
- :func:`ico_suffix` (R26): the -ico/-iaco suffix evidence of a lemma
  (see its docstring).
- :func:`suffixoid_evidence` (R27): -logo/-fago/-urgo lemmas and the
  evidence that they carry the learned suffixoid, for the shared
  learned-suffixoid detector of :mod:`lexpipe.suffixoids` (which also
  holds the person/thing classifier).
"""

from __future__ import annotations

import re
from collections.abc import Container, Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Final

from lexpipe import etymology as core
from lexpipe.etymology import (
    BORROWED,
    INHERITED,
    INTERNAL,
    STRATA,
    UNKNOWN,
    Step,
    first_transmission,
    fold,
    language_family,
)

__all__ = [
    "BORROWED",
    "INHERITED",
    "INTERNAL",
    "STRATA",
    "UNKNOWN",
    "Step",
    "first_transmission",
    "fold",
    "language_family",
]

PROFILE: Final[core.EtymologyProfile] = core.ITALIAN
LATIN_GREEK: Final[tuple[str, ...]] = core.LATIN_GREEK
SOURCE_KINDS: Final[frozenset[str]] = core.SOURCE_KINDS
_VERB_RE: Final[re.Pattern[str]] = re.compile(
    r"(?:are|ere|ire|arsi|ersi|irsi)\Z"
)
# A back-formation's source counts as a verb only with an infinitive
# ending (-rre included: porre); metallurgo is a back-formation from the
# noun metallurgia, not a deverbal noun.
_BF_VERB_RE: Final[re.Pattern[str]] = re.compile(
    r"(?:are|ere|ire|rre|arsi|ersi|irsi|rsi)\Z"
)
_DEVERBAL_TEXT_RE: Final[re.Pattern[str]] = re.compile(
    r"\b(?:[Dd]everbal|[Bb]ack-formation)\s+(?:from|of)\s+([^\s,.;:()“”\"]+)"
)
_PARTICIPLE_TEXT_RE: Final[re.Pattern[str]] = re.compile(
    r"\b(?:short\s+)?past participle of (?:the verb\s+)?" r"([^\s,.;:()“”\"]+)"
)
_FIRST_CONJUGATION_RE: Final[re.Pattern[str]] = re.compile(r"(?:are|arsi)\Z")


# Italian compound templates that name their parts from argument 1
# (no language argument): {{it-verb-obj|forare|sacco}}.
COMPOUND_TEMPLATES: Final[frozenset[str]] = frozenset(
    {"it-verb-obj", "it-verb-verb"}
)


def italian_templates(
    templates: Iterable[Any],
) -> Iterator[Any]:
    """*templates* with each :data:`COMPOUND_TEMPLATES` template written
    as the ``compound`` template the shared parser reads (language code
    first, the parts shifted by one; other arguments dropped). Other
    items are passed through unchanged."""
    for template in templates:
        if (
            isinstance(template, Mapping)
            and template.get("name") in COMPOUND_TEMPLATES
        ):
            raw = template.get("args")
            args = raw if isinstance(raw, Mapping) else {}
            parts = sorted(
                (int(key), str(value))
                for key, value in args.items()
                if str(key).isdigit() and int(key) >= 1
            )
            shifted = {str(index + 1): value for index, value in parts}
            yield {
                "name": "compound",
                "args": {"1": PROFILE.code, **shifted},
            }
        else:
            yield template


def steps(templates: Iterable[Mapping[str, Any]]) -> list[Step]:
    """The etymological steps of an Italian entry
    (:func:`lexpipe.etymology.steps` over :func:`italian_templates`)."""
    return core.steps(italian_templates(templates), PROFILE.code)


def etymon(found: Sequence[Step]) -> tuple[str, str] | None:
    """``(language, term)`` of the first source step
    (:func:`lexpipe.etymology.etymon`)."""
    return core.etymon(found, PROFILE.code)


def etymon_key(found: Sequence[Step]) -> str:
    """``family:term`` of :func:`etymon` for root keys ("" if none).

    The term is folded (case, diacritics, a leading ``*``) but keeps
    its ending: different source words stay apart even when they differ
    only in their inflection class (Latin arca 'chest' and arcus 'bow',
    orca and Orcus, dulcis and dulcō), as do different etymologies
    (antīquus vs antīcus; Latin parcus vs Old French parc). Gender pairs
    whose etymons differ (amīcus, amīca) share a root through the
    feminine's stated counterpart (:mod:`.lexicon`), not through their
    etymons.
    """
    found_etymon = etymon(found)
    if found_etymon is None:
        return ""
    lang, term = found_etymon
    key = fold(term).lstrip("*").strip()
    return f"{language_family(lang)}:{key}"


def stratum(
    found: Sequence[Step], text: str, *, word: str = ""
) -> tuple[str, str]:
    """``(stratum, basis)`` per R22 (:func:`lexpipe.etymology.stratum`
    with the Italian profile; *word* locates the prose of an etymology
    tree)."""
    return core.stratum(found, text, PROFILE, word=word)


def deverbal_verb(found: Sequence[Step], text: str) -> str:
    """The verb a lexeme is formed from ("" if none is stated)."""
    for step in found:
        if step.kind == "deverbal" and step.terms:
            return step.terms[0]
        if (
            step.kind == "bf"
            and step.terms
            and _BF_VERB_RE.search(fold(step.terms[0]))
        ):
            return step.terms[0]
        if (
            step.kind == "from"
            and step.lang in ("", "it")
            and step.terms
            and _VERB_RE.search(fold(step.terms[0]))
        ):
            return step.terms[0]
    for pattern in (_DEVERBAL_TEXT_RE, _PARTICIPLE_TEXT_RE):
        match = pattern.search(text or "")
        if match is not None and _VERB_RE.search(fold(match.group(1))):
            return match.group(1)
    return ""


def is_first_conjugation(verb: str) -> bool:
    """-are or -arsi (accents ignored)."""
    return bool(_FIRST_CONJUGATION_RE.search(fold(verb)))


# ---------------------------------------------------------------------------
# -ico suffix (R26)
# ---------------------------------------------------------------------------

ICO: Final[str] = "-ico"
IACO: Final[str] = "-iaco"
_ICO_ENDINGS: Final[tuple[str, ...]] = ("ico", "ica")
_IACO_ENDINGS: Final[tuple[str, ...]] = ("iaco", "iaca")
_LATIN_ICUS_RE: Final[re.Pattern[str]] = re.compile(
    r"(?:icus|icum|ica|ικος|ικον|ικη|iacus|iacum|ιακος)\Z"
)


@dataclass(frozen=True, slots=True)
class SuffixEvidence:
    """Result of :func:`ico_suffix`."""

    lemma_suffix: str
    target_is_suffix: bool
    evidence: str
    raw_suffix: str
    base: str


NO_SUFFIX: Final[SuffixEvidence] = SuffixEvidence("", False, "", "", "")


def ico_family(key: str) -> str:
    """``-iaco``/``-ico`` if the lemma key ends like one (else "")."""
    if key.endswith(_IACO_ENDINGS):
        return IACO
    if key.endswith(_ICO_ENDINGS):
        return ICO
    return ""


# Templates whose last part is a suffix whether or not it is written with
# a hyphen ({{suf|it|Fær Øer|ico}}, {{confix|it|deca|edrico}}: the
# templates add the hyphen themselves).
SUFFIX_LAST_TEMPLATES: Final[frozenset[str]] = frozenset(
    {"suffix", "suf", "confix", "con"}
)


def _suffix_part(step: Step) -> str:
    """The step's last part as a suffix (with its hyphen), or "" when
    the template does not mark it as one."""
    last = step.terms[-1]
    if last.startswith("-"):
        return last
    if step.template.rstrip("+") in SUFFIX_LAST_TEMPLATES and last:
        return "-" + last
    return ""


def _template_suffix(found: Sequence[Step]) -> tuple[str, str, str]:
    """``(kind, suffix, base)`` of the first af/surf step whose last part
    is an -ico-family suffix (:func:`_suffix_part`)."""
    for step in found:
        if step.kind not in ("af", "surf") or len(step.terms) < 2:
            continue
        suffix = _suffix_part(step)
        if suffix and fold(suffix).endswith(_ICO_ENDINGS + _IACO_ENDINGS):
            return step.kind, suffix, step.terms[0]
    return "", "", ""


def ico_template(found: Sequence[Step]) -> bool:
    """Does an affix or surface template form the word with an -ico-
    family suffix (X + -ico, X + -istico)?"""
    return bool(_template_suffix(found)[0])


def ico_suffix(
    key: str,
    found: Sequence[Step],
    *,
    ico_stressed: bool | None,
    adjectival: bool,
    deverbal: bool,
    italian_lemmas: frozenset[str] = frozenset(),
) -> SuffixEvidence:
    """-ico suffix status of a lemma, by morphology (R26).

    A lemma whose key ends in -ico/-ica (-iaco/-iaca) is a suffix target
    when the suffix is unstressed (*ico_stressed* is not True: the i of
    -ico does not carry the stress, :func:`.lexicon.ico_stressed`;
    unknown stress is allowed, and a penultimate stress on the vowel
    before -ico — arcaico /arˈkaj.ko/ — leaves the suffix unstressed)
    and either

      - ``template``/``surface``: an affix or surface-analysis template
        ends in an -ico-family suffix (X + -ico, X + -istico, …; the
        last part of a ``suffix``/``confix`` template is a suffix
        without its hyphen too: feroico, decaedrico); or
      - ``latin_icus``: its first source is a Latin or Greek -icus/-ικός
        word that was *not* inherited (borrowed or derived), and the
        lexeme is adjectival (an adjective, or a noun whose root has an
        adjective: *adjectival*).

    Never a target: deverbal nouns (carico, from caricare) and inherited
    words (the simplex Latin nouns: medico < medicus, inherited).

    ``base`` is the template's base, or — for ``latin_icus`` — an
    Italian lemma (in *italian_lemmas*) equal to the stem plus o/a/e/ia
    (elettrico → elettro, classico → classe); R26's "attested Italian
    base". It is recorded, not required: Wiktionary rarely states the
    base of a borrowed -icus adjective (pubblico < pūblicus has none),
    so requiring it would drop the bulk of the learned -ico class.
    """
    family = ico_family(key)
    if not family or deverbal or ico_stressed is True:
        return NO_SUFFIX
    first = first_transmission(found)
    if first is not None and first.kind == "inh":
        return NO_SUFFIX
    kind, raw, base = _template_suffix(found)
    if kind:
        evidence = "template" if kind == "af" else "surface"
        return SuffixEvidence(family, True, evidence, raw, base)
    source = etymon(found)
    if source is None or not adjectival:
        return NO_SUFFIX
    lang, term = source
    if language_family(lang) not in LATIN_GREEK:
        return NO_SUFFIX
    if not _LATIN_ICUS_RE.search(fold(term)):
        return NO_SUFFIX
    return SuffixEvidence(
        family, True, "latin_icus", "", ico_base(key, italian_lemmas)
    )


# Endings of an attested Italian base of an -ico adjective (elettro →
# elettrico, classe → classico, fitoterapia → fitoterapico, tallio →
# tallico).
ICO_BASE_ENDINGS: Final[tuple[str, ...]] = ("o", "a", "e", "ia", "io")


def ico_base(key: str, italian_lemmas: Container[str]) -> str:
    """An Italian lemma (in *italian_lemmas*) equal to the stem of an
    -ico/-iaco word plus :data:`ICO_BASE_ENDINGS` ("" if none): the
    attested base R26 asks for."""
    family = ico_family(key)
    if not family:
        return ""
    stem = key[: -len(family) + 1]
    return next(
        (
            stem + ending
            for ending in ICO_BASE_ENDINGS
            if stem + ending in italian_lemmas and stem + ending != key
        ),
        "",
    )


# ---------------------------------------------------------------------------
# -logo / -fago / -urgo (R27)
# ---------------------------------------------------------------------------

SUFFIXOIDS: Final[tuple[str, ...]] = ("-logo", "-fago", "-urgo")
# Suffixoid → the key of lexpipe.suffixoids.GLOSS_REFLEX.
SUFFIXOID_KEYS: Final[Mapping[str, str]] = {
    "-logo": "log",
    "-fago": "fag",
    "-urgo": "urg",
}
_MASCULINE_ENDINGS: Final[Mapping[str, str]] = {
    "logo": "-logo",
    "fago": "-fago",
    "urgo": "-urgo",
}
_FEMININE_ENDINGS: Final[Mapping[str, str]] = {
    "loga": "-logo",
    "faga": "-fago",
    "urga": "-urgo",
}
_SUFFIXOID_PART_RE: Final[re.Pattern[str]] = re.compile(
    r"-(?:log|fag|urg)[oa]\Z"
)
_SUFFIXOID_ETYMON_RE: Final[re.Pattern[str]] = re.compile(
    r"(?:logus|logos|λογος|phagus|phagos|φαγος|urgus|ourgos|ουργος)\Z"
)
# Evidence labels (the basis lexpipe.suffixoids.learned_suffixoid
# reports for an Italian lexeme).
EVIDENCE_ENDING: Final[str] = "ending"
EVIDENCE_TEMPLATE: Final[str] = "etymology (template)"
EVIDENCE_ETYMON: Final[str] = "etymology (etymon)"
EVIDENCE_MASCULINE: Final[str] = "masculine counterpart"


def suffixoid_evidence(
    key: str,
    found: Sequence[Step] = (),
    *,
    deverbal: bool = False,
    masculine_forms: Sequence[str] = (),
) -> tuple[str, str]:
    """``(suffixoid, evidence)`` of a -logo/-fago/-urgo lemma, or
    ``("", "")``.

    A masculine lemma in -logo/-fago/-urgo qualifies by its ending
    (teologo, catalogo, chirurgo; evidence ``ending``: Italian has no
    folk lookalikes ending so, apart from deverbal nouns), unless it is
    deverbal (spurgo < spurgare). A feminine in -loga/-faga/-urga needs
    evidence: an affix template with the suffixoid, a Greek/Latin
    etymon in -logos/-phagos/-ourgos, or a stated masculine
    counterpart that qualifies (psicologa → psicologo); this keeps out
    purga, riloga, egloga.
    """
    if deverbal:
        return "", ""
    for ending, name in _MASCULINE_ENDINGS.items():
        if key.endswith(ending) and len(key) > len(ending) + 1:
            return name, EVIDENCE_ENDING
    for ending, name in _FEMININE_ENDINGS.items():
        if not key.endswith(ending) or len(key) <= len(ending) + 1:
            continue
        parts = (
            fold(term)
            for step in found
            if step.kind in ("af", "surf")
            for term in step.terms
        )
        if any(_SUFFIXOID_PART_RE.search(part) for part in parts):
            return name, EVIDENCE_TEMPLATE
        source = etymon(found)
        if source is not None and _SUFFIXOID_ETYMON_RE.search(fold(source[1])):
            return name, EVIDENCE_ETYMON
        if any(
            fold(form).endswith(ending[:-1] + "o") for form in masculine_forms
        ):
            return name, EVIDENCE_MASCULINE
    return "", ""


def suffixoid(
    key: str,
    found: Sequence[Step] = (),
    *,
    deverbal: bool = False,
    masculine_forms: Sequence[str] = (),
) -> str:
    """``-logo``/``-fago``/``-urgo`` (R27) or "" (the suffixoid of
    :func:`suffixoid_evidence`)."""
    return suffixoid_evidence(
        key, found, deverbal=deverbal, masculine_forms=masculine_forms
    )[0]
