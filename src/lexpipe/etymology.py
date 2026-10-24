"""Wiktionary etymology templates and the R22 stratum, for every language.

English Wiktionary states etymologies with templates that wiktextract
keeps in ``etymology_templates``. One parser (:func:`steps`) turns
them into :class:`Step`\\s and one classifier (:func:`stratum`) gives
the R22 stratum of an entry, so Italian and Czech (and any later
plugin) classify etymologies identically. A plugin supplies only an
:class:`EtymologyProfile`: its language code and name, the ancestors
whose ``der``/``uder`` count as inheritance, and the letters of a
native word in "From <word> …" sentences (:data:`ITALIAN`,
:data:`CZECH`).

Steps
-----
Each step is a transmission kind, a source language and the terms
(:class:`Step`):

  - classic templates ``inh`` (and ``inh-lite``) ``bor`` ``lbor``
    ``slbor`` ``obor`` ``ubor`` ``der`` ``uder`` ``calque`` (``cal``,
    ``clq``) ``psm`` ``semantic loan`` (``+`` variants included):
    language in argument 2, term in argument 3;
  - affix templates ``af``/``affix``/``suffix``/``prefix``/``confix``/
    ``compound`` (and their short names; kind ``af``), ``surf`` and
    ``surface analysis`` ("By surface analysis", kind ``surf``),
    ``blend``, ``clipping``,
    ``univerbation``: the parts in order. A part keeps an inline
    ``lang:`` prefix as written; its language (the prefix or the
    template's ``langN`` argument) is recorded in ``term_langs``;
  - ``deverbal`` and the language's own ``<code>-deverbal`` (``…
    fpp``) templates (kind ``deverbal``; the term is argument 1, or
    argument 2 after the language code), ``back-form``/``bf`` (kind
    ``bf``);
  - ``onomatopoeic``/``onom``, ``named-after``, ``pseudo-loan``;
  - ``{{ety}}``/``{{etymon}}``: ``:kind`` tokens open a step (with the
    aliases of :data:`ETY_KIND_ALIASES`; a ``<surf>`` modifier on the
    kind makes it ``surf``), the tokens after them are its terms,
    written ``lang:term<modifiers>``; a term without a prefix is in
    the entry's language (``{{ety|it|:af|banca|-ario}}``,
    ``{{ety|cs|:inh|sla-pro:*strašiti|:af<surf>|strach|-it}}``). The
    step's language is that of its first prefixed term.

Terms lose their ``<modifiers>`` and parenthesised parts at a word
edge (``mieti(trice)`` → ``mieti``, ``(nux) Abellāna`` →
``Abellāna``); a parenthesis inside a word marks an optional letter
and is kept (``*ga(v)ota``, ``cal(i)dus``, ``*o(b)krajь``). Wikitext
artifacts that wiktextract keeps in duplicated templates
(``['zlw-ocs:', <LINK…>]``) are skipped; identical consecutive steps
(a wrapper and its expansion) are kept once. Cognates, doublets, roots
and mentions are not steps.

Stratum (R22)
-------------
:func:`stratum` returns ``(stratum, basis)`` with stratum ∈
:data:`STRATA`. The first transmission step decides
(:func:`first_transmission`):

  - ``inh`` → inherited;
  - a borrowing (``bor``, ``lbor``, …, ``calque``) or ``surf`` →
    borrowed whole;
  - ``der``/``uder``/``:from`` by language: the language itself →
    internal formation; an ancestor
    (:attr:`EtymologyProfile.inherited_sources`) → inherited, basis
    ``kind:lang``; any other language → borrowed whole, basis
    ``kind:lang``;
  - an internal kind (affixation, compound, deverbal, back-formation,
    clipping, blend, univerbation, onomatopoeia, named-after):
    borrowed whole (basis ``internal_analysis``) when the etymology
    text says "Internally, from …" or *opens* with a borrowing ("From
    Latin X + -ní", "Borrowed from …", "Learned borrowing from …",
    "Calque of …": the word came from elsewhere and is only analysed
    internally); borrowed whole (basis ``surface_analysis``) when the
    text opens with "By surface analysis" (the affixation is only the
    surface analysis), or has it later and the step names no base of
    the language itself (or of an ancestor) — the surface phrase then
    analyses a borrowing; otherwise internal formation. A step that
    names an internal base keeps its affixation whatever a later
    surface phrase says ("From antico + -ità. By surface analysis,
    borrowed from Latin …"; "From elektro- + mobilita. By surface
    analysis, elektromobil + -ita.").

Without a transmission template the text decides, from its first
sentence: "Inherited from …" or "From <ancestor>" → inherited; a
borrowing phrase ("Borrowed from …", "From <foreign language> …",
"By surface analysis …") → borrowed whole; "Internally, from …"
anywhere, or "By surface analysis" after a sentence that does not open
with an affixation of a native word → borrowed whole
(``surface_analysis``); "From <native word> (…", "From <word> +
-suffix", "From <word>." → internal formation; otherwise unknown.
wiktextract prefixes the texts of ``{{etymon}}`` trees with the tree
("Etymology tree / Proto-Slavic *X / … / Czech <word>"); only the prose
after the tree's last node is read (:func:`prose`).

:func:`summary` writes the deciding step compactly for tracked raw
tables (``inh:la:amīcus``, ``af``).

Reconciliation
--------------
The Italian and Czech plugins had separate copies of this code. The
merged behaviour (``docs/multilingual_design.md``, "Corrections since
adoption") is:

  - the union of both template tables (Italian gains ``inh-lite``,
    ``onom``, ``afeq``, ``ellipsis``; Czech gains ``con``), the Czech
    artifact skipping, ``<surf>`` kinds and parenthesis stripping
    (Czech stripped every parenthesis, Italian none; the merged rule
    strips them at word edges and keeps optional letters inside a
    word);
  - ``:from`` steps are classified by language (Italian), not always
    internal (Czech);
  - internal kinds whose text opens with a borrowing or says
    "Internally, from" are borrowed whole (Italian; new for Czech:
    "From Latin āctīvus + Czech -ní");
  - "By surface analysis" makes an internal kind borrowed whole only
    when the step names no native base (generalizing the Italian
    antichità exception, which kept only "…, borrowed from" phrases
    internal: altoatesino, antropometrico are internal now);
  - der/uder bases are written ``kind:lang`` (Czech) for both;
  - the text fallback reads inheritance, borrowing, surface analysis
    and native "From X" sentences (Czech), with language names for
    both, and strips etymology trees.
"""

from __future__ import annotations

import functools
import re
import unicodedata
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Final

INHERITED: Final[str] = "inherited"
BORROWED: Final[str] = "borrowed whole"
INTERNAL: Final[str] = "internal formation"
UNKNOWN: Final[str] = "unknown"
STRATA: Final[tuple[str, ...]] = (INHERITED, BORROWED, INTERNAL, UNKNOWN)

CLASSIC_TEMPLATES: Final[Mapping[str, str]] = {
    "inh": "inh",
    "inh-lite": "inh",
    "bor": "bor",
    "lbor": "lbor",
    "slbor": "slbor",
    "obor": "obor",
    "ubor": "ubor",
    "der": "der",
    "uder": "uder",
    "calque": "calque",
    "cal": "calque",
    "clq": "calque",
    "psm": "psm",
    "semantic loan": "semantic_loan",
    "sl": "semantic_loan",
}
AFFIX_TEMPLATES: Final[Mapping[str, str]] = {
    "af": "af",
    "affix": "af",
    "suffix": "af",
    "suf": "af",
    "prefix": "af",
    "pre": "af",
    "confix": "af",
    "con": "af",
    "compound": "af",
    "com": "af",
    "surf": "surf",
    "surface analysis": "surf",
    "blend": "blend",
    "clipping": "clipping",
    "univerbation": "univerbation",
}
DEVERBAL_TEMPLATES: Final[Mapping[str, str]] = {
    "deverbal": "deverbal",
    "back-form": "bf",
    "back-formation": "bf",
    "bf": "bf",
}
# <code>-deverbal templates (it-deverbal, it-deverbal fpp).
_LANGUAGE_DEVERBAL_SUFFIXES: Final[tuple[str, ...]] = (
    "-deverbal",
    "-deverbal fpp",
)
OTHER_TEMPLATES: Final[Mapping[str, str]] = {
    "onomatopoeic": "onomatopoeic",
    "onom": "onomatopoeic",
    "named-after": "named_after",
    "pseudo-loan": "pseudo_loan",
}
ETY_KIND_ALIASES: Final[Mapping[str, str]] = {
    "affix": "af",
    "afeq": "af",
    "suffix": "af",
    "prefix": "af",
    "compound": "af",
    "back-form": "bf",
    "cal": "calque",
    "clq": "calque",
    "univ": "univerbation",
    "ellipsis": "clipping",
}

BORROWING_KINDS: Final[frozenset[str]] = frozenset(
    {
        "bor",
        "lbor",
        "slbor",
        "obor",
        "ubor",
        "calque",
        "psm",
        "semantic_loan",
        "pseudo_loan",
    }
)
SOURCE_KINDS: Final[frozenset[str]] = BORROWING_KINDS | {"inh", "der", "uder"}
INTERNAL_KINDS: Final[frozenset[str]] = frozenset(
    {
        "af",
        "deverbal",
        "bf",
        "from",
        "blend",
        "clipping",
        "univerbation",
        "onomatopoeic",
        "named_after",
    }
)
TRANSMISSION_KINDS: Final[frozenset[str]] = (
    SOURCE_KINDS | INTERNAL_KINDS | {"surf"}
)
LATIN_GREEK: Final[tuple[str, ...]] = ("la", "grc")

# English names of source languages in etymology sentences ("From Latin
# …", "From Old French …"), and the modifiers that may precede them.
FOREIGN_LANGUAGE_NAMES: Final[tuple[str, ...]] = (
    "Latin",
    "Greek",
    "French",
    "English",
    "Italian",
    "Spanish",
    "Portuguese",
    "German",
    "Arabic",
    "Occitan",
    "Provençal",
    "Catalan",
    "Persian",
    "Turkish",
    "Hebrew",
    "Russian",
    "Polish",
    "Hungarian",
    "Dutch",
    "Church Slavonic",
    "Slovak",
    "Yiddish",
    "Czech",
)
LANGUAGE_MODIFIERS: Final[tuple[str, ...]] = (
    "Late",
    "Medieval",
    "New",
    "Vulgar",
    "Ecclesiastical",
    "Classical",
    "Old",
    "Middle",
    "Ancient",
    "Modern",
    "Byzantine",
)
TREE_HEADER: Final[str] = "Etymology tree"

_LANG_TERM_RE: Final[re.Pattern[str]] = re.compile(
    r"([a-z]{2,3}(?:-[a-z]{2,4})*):(.+)", re.DOTALL
)
# A parenthesised part at a word edge; one inside a word (an optional
# letter: *ga(v)ota) is kept.
_PARENTHESIS_RE: Final[re.Pattern[str]] = re.compile(
    r"(?<!\w)\([^()]*\)|\([^()]*\)(?!\w)"
)
_SURFACE_RE: Final[re.Pattern[str]] = re.compile(r"By surface analysis")
_INTERNALLY_RE: Final[re.Pattern[str]] = re.compile(r"\bInternally,? from\b")
_SURFACE_START_RE: Final[re.Pattern[str]] = re.compile(
    r"\s*By surface analysis"
)


@dataclass(frozen=True, slots=True)
class EtymologyProfile:
    """What the parser and the classifier need to know of a language.

    ``code`` is the Wiktionary code of the language (a term without a
    ``lang:`` prefix is in it; ``der`` from it is an internal
    formation); ``name`` its English name in etymology sentences
    ("Czech -ický", "Etymology tree … Czech kniha");
    ``inherited_sources`` the ancestor codes whose ``der``/``uder``
    are inheritance; ``inherited_names`` the ancestors' English names
    ("From Vulgar Latin …" is inheritance); ``native_letters`` a
    regex character class of the lowercase letters of a native word
    ("From sníh (…)" is an internal formation).
    """

    code: str
    name: str
    inherited_sources: frozenset[str]
    inherited_names: tuple[str, ...]
    native_letters: str

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "inherited_sources", frozenset(self.inherited_sources)
        )
        object.__setattr__(
            self, "inherited_names", tuple(self.inherited_names)
        )


ITALIAN: Final[EtymologyProfile] = EtymologyProfile(
    code="it",
    name="Italian",
    inherited_sources=frozenset({"la-vul", "roa-oit"}),
    inherited_names=("Vulgar Latin", "Old Italian"),
    native_letters="a-zàèéìíòóù",
)
CZECH: Final[EtymologyProfile] = EtymologyProfile(
    code="cs",
    name="Czech",
    inherited_sources=frozenset(
        {"zlw-ocs", "zlw-pro", "sla-pro", "ine-bsl-pro", "ine-pro"}
    ),
    inherited_names=(
        "Old Czech",
        "Proto-Slavic",
        "Proto-West Slavic",
        "Proto-Balto-Slavic",
        "Proto-Indo-European",
    ),
    native_letters="a-záčďéěíňóřšťúůýž",
)


@dataclass(frozen=True, slots=True)
class _Patterns:
    inherited: re.Pattern[str]
    borrowed: re.Pattern[str]
    native_from: re.Pattern[str]
    leading_affixation: re.Pattern[str]


@functools.lru_cache(maxsize=None)
def _patterns(profile: EtymologyProfile) -> _Patterns:
    foreign = sorted(
        (name for name in FOREIGN_LANGUAGE_NAMES if name != profile.name),
        key=lambda name: (-len(name), name),
    )
    languages = "|".join(re.escape(name) for name in foreign)
    modifiers = "|".join(LANGUAGE_MODIFIERS)
    ancestors = "|".join(re.escape(name) for name in profile.inherited_names)
    inherited = r"Inherited from" + (
        rf"|From (?:{ancestors})" if ancestors else ""
    )
    letters = profile.native_letters
    return _Patterns(
        inherited=re.compile(rf"\s*(?:{inherited})\b"),
        borrowed=re.compile(
            r"\s*(?:Borrowed from|Learned borrowing from|Semi-learned "
            r"borrowing from|Unadapted borrowing from|Calque of|By surface "
            rf"analysis|From (?:(?:{modifiers})\s+)*(?:{languages}))\b"
        ),
        native_from=re.compile(
            rf"\s*From [{letters}]+(?:\s+\(|\s*\+\s*-|\.|,)"
        ),
        leading_affixation=re.compile(
            rf"\s*From [{letters}]+-?(?:\s+\([^)]*\))?\s*\+"
        ),
    )


@dataclass(frozen=True, slots=True)
class Step:
    """One etymological step: kind, source language, terms.

    ``term_langs`` gives the language of each term when the template
    says ("" = the step's language): an ``{{ety}}`` prefix
    (``grc:δαίμων``), an inline prefix of an affix part
    (``la:chronicus``) or the template's ``langN`` argument.
    """

    kind: str
    lang: str
    terms: tuple[str, ...]
    template: str
    term_langs: tuple[str, ...] = ()

    def term_lang(self, index: int) -> str:
        """The language of term *index* ("" = the step's language)."""
        if index < len(self.term_langs):
            return self.term_langs[index]
        return ""


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


def clean_term(text: str) -> str:
    """A term without its ``<modifiers>`` and the parenthesised parts at
    a word edge (NFC, stripped; module docstring); a ``lang:`` prefix
    is kept as written."""
    cut = text.find("<")
    term = text[:cut] if cut >= 0 else text
    term = _PARENTHESIS_RE.sub("", term)
    return unicodedata.normalize("NFC", term.strip())


def _numbered(args: Mapping[str, str], start: int) -> list[tuple[int, str]]:
    return sorted(
        (int(key), value)
        for key, value in args.items()
        if key.isdigit() and int(key) >= start
    )


def _args(template: Mapping[str, Any]) -> dict[str, str]:
    raw = template.get("args")
    if not isinstance(raw, Mapping):
        return {}
    return {str(k): str(v) for k, v in raw.items()}


def _kind_token(token: str) -> str:
    raw = token[1:]
    if "<surf>" in raw:
        return "surf"
    name = clean_term(raw).lower()
    return ETY_KIND_ALIASES.get(name, name)


def _is_artifact(token: str) -> bool:
    return token.startswith("[") or "<LINK" in token


def _ety_steps(name: str, args: Mapping[str, str], code: str) -> list[Step]:
    out: list[Step] = []
    kind = ""
    lang = ""
    terms: list[str] = []
    langs: list[str] = []
    artifact = False

    def flush() -> None:
        # A step whose only term was a wikitext artifact says nothing.
        if terms or (kind and not artifact):
            out.append(
                Step(
                    kind or "from",
                    lang or code,
                    tuple(terms),
                    name,
                    tuple(langs),
                )
            )

    for _, token in _numbered(args, 2):
        stripped = token.strip()
        if _is_artifact(stripped):
            artifact = True
            continue
        if stripped.startswith(":"):
            flush()
            kind = _kind_token(stripped)
            lang, terms, langs, artifact = "", [], [], False
            continue
        match = _LANG_TERM_RE.fullmatch(stripped)
        term_lang = ""
        if match is not None and "<" not in match.group(1):
            term_lang = match.group(1)
            if not lang:
                lang = term_lang
            term = clean_term(match.group(2))
        else:
            term = clean_term(stripped)
        if term:
            terms.append(term)
            langs.append(term_lang)
    flush()
    return out


def _parts(args: Mapping[str, str]) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """The parts of an affix template (argument 2 on) and their
    languages (``langN`` of part N, or an inline ``lang:`` prefix)."""
    parts: list[str] = []
    langs: list[str] = []
    for index, raw in _numbered(args, 2):
        part_lang = args.get(f"lang{index - 1}", "")
        match = _LANG_TERM_RE.fullmatch(raw.strip())
        if match is not None and "<" not in match.group(1) and not part_lang:
            part_lang = match.group(1)
        term = clean_term(raw)
        if term:
            parts.append(term)
            langs.append(part_lang)
    return tuple(parts), tuple(langs)


def _is_deverbal(base: str, code: str) -> bool:
    if base in DEVERBAL_TEMPLATES:
        return True
    prefix = code + "-"
    return base.startswith(prefix) and (
        "-" + base[len(prefix) :] in _LANGUAGE_DEVERBAL_SUFFIXES
    )


def template_steps(template: Mapping[str, Any], code: str) -> list[Step]:
    """The steps one etymology template states (module docstring);
    *code* is the entry's language (:attr:`EtymologyProfile.code`)."""
    name = str(template.get("name") or "")
    args = _args(template)
    base = name[:-1] if name.endswith("+") else name
    if base in CLASSIC_TEMPLATES:
        term = clean_term(args.get("3", ""))
        return [
            Step(
                CLASSIC_TEMPLATES[base],
                args.get("2", ""),
                (term,) if term else (),
                name,
            )
        ]
    if base in AFFIX_TEMPLATES:
        parts, langs = _parts(args)
        return [
            Step(AFFIX_TEMPLATES[base], args.get("1", ""), parts, name, langs)
        ]
    if _is_deverbal(base, code):
        first = args.get("1", "")
        term = clean_term(args.get("2", "") if first == code else first)
        kind = DEVERBAL_TEMPLATES.get(base, "deverbal")
        return [Step(kind, code, (term,) if term else (), name)]
    if base in OTHER_TEMPLATES:
        parts, langs = _parts(args)
        return [
            Step(OTHER_TEMPLATES[base], args.get("1", ""), parts, name, langs)
        ]
    if base in ("ety", "etymon"):
        return _ety_steps(base, args, code)
    return []


def steps(templates: Iterable[Any], code: str) -> list[Step]:
    """The steps of an entry, in source order (module docstring).

    Consecutive duplicates (same kind, language and terms) are dropped
    and non-mapping items of the untrusted input are skipped.
    """
    out: list[Step] = []
    for template in templates:
        if not isinstance(template, Mapping):
            continue
        for step in template_steps(template, code):
            if out and (out[-1].kind, out[-1].lang, out[-1].terms) == (
                step.kind,
                step.lang,
                step.terms,
            ):
                continue
            out.append(step)
    return out


# ---------------------------------------------------------------------------
# Derived values
# ---------------------------------------------------------------------------


def first_transmission(found: Sequence[Step]) -> Step | None:
    """The first step that says how the word came about."""
    for step in found:
        if step.kind in TRANSMISSION_KINDS:
            return step
    return None


def language_family(lang: str) -> str:
    """``la`` for every Latin variety (la, la-lat, la-med, la-vul…),
    ``grc`` for Ancient Greek varieties, else the code itself."""
    for family in LATIN_GREEK:
        if lang == family or lang.startswith(family + "-"):
            return family
    return lang


def fold(text: str) -> str:
    """Lowercase, diacritics removed (comparison only)."""
    decomposed = unicodedata.normalize("NFKD", text.lower())
    return "".join(ch for ch in decomposed if not unicodedata.combining(ch))


def etymon(found: Sequence[Step], code: str) -> tuple[str, str] | None:
    """``(language, term)`` of the first source step (inherited,
    borrowed or derived from another language, or ``:from`` a word of
    another language than *code*)."""
    for step in found:
        if step.kind in SOURCE_KINDS and step.terms and step.lang:
            return step.lang, step.terms[0]
        if step.kind == "from" and step.terms and step.lang not in ("", code):
            return step.lang, step.terms[0]
    return None


def summary(found: Sequence[Step], code: str) -> str:
    """A compact summary of the first transmission step, for tracked
    raw tables: ``kind:lang:term`` for a source step (``inh:la:amīcus``,
    ``bor:fr:parc``, ``der:la:ēlectricus``) or a ``:from`` step from
    another language than *code*; the kind alone for an internal one
    (``af``, ``deverbal``, ``surf``); "" without a transmission
    template."""
    first = first_transmission(found)
    if first is None:
        return ""
    if first.kind in SOURCE_KINDS or (
        first.kind == "from" and first.lang not in ("", code)
    ):
        term = first.terms[0] if first.terms else ""
        return f"{first.kind}:{first.lang}:{term}"
    return first.kind


def prose(
    text: str, word: str = "", profile: EtymologyProfile | None = None
) -> str:
    """The etymology sentences of a wiktextract ``etymology_text``.

    A text that starts with :data:`TREE_HEADER` lists the tree of an
    ``{{etymon}}`` first, one node per line, ending with the entry's
    own node ("Czech kniha"); the prose follows it. Without the
    entry's node (no *word* or *profile*), the last line is the prose.
    Other texts are returned unchanged.
    """
    if not text.startswith(TREE_HEADER):
        return text
    lines = text.split("\n")
    if word and profile is not None:
        node = f"{profile.name} {word}"
        found = [index for index, line in enumerate(lines) if line == node]
        if found:
            return "\n".join(lines[found[-1] + 1 :])
    return lines[-1] if len(lines) > 1 else ""


def _native(lang: str, profile: EtymologyProfile) -> bool:
    return lang in ("", profile.code) or lang in profile.inherited_sources


def names_native_base(step: Step, profile: EtymologyProfile) -> bool:
    """Does an internal step name a base of the language itself (or of
    an ancestor)? Affixes (``-ità``, ``a-``) are not bases; a part
    with a foreign language (``grc:χρῶμα``, ``lang2=grc``) is not
    native."""
    for index, term in enumerate(step.terms):
        if not term or term.startswith("-") or term.endswith("-"):
            continue
        if _native(step.term_lang(index), profile):
            return True
    return False


def _no_template_stratum(
    text: str, profile: EtymologyProfile, patterns: _Patterns
) -> tuple[str, str]:
    if patterns.inherited.match(text):
        return INHERITED, "text"
    if patterns.borrowed.match(text):
        return BORROWED, "text"
    if _INTERNALLY_RE.search(text):
        return BORROWED, "surface_analysis"
    if _SURFACE_RE.search(text) and not patterns.leading_affixation.match(
        text
    ):
        return BORROWED, "surface_analysis"
    if patterns.native_from.match(text):
        return INTERNAL, "text"
    return UNKNOWN, ""


def stratum(
    found: Sequence[Step],
    text: str,
    profile: EtymologyProfile,
    *,
    word: str = "",
) -> tuple[str, str]:
    """``(stratum, basis)`` per R22 (module docstring).

    *found* are the entry's :func:`steps`, *text* its
    ``etymology_text`` (an etymology tree is skipped; *word* finds the
    tree's last node, :func:`prose`).
    """
    sentences = prose(text or "", word, profile)
    patterns = _patterns(profile)
    first = first_transmission(found)
    if first is None:
        return _no_template_stratum(sentences, profile, patterns)
    kind = first.kind
    if kind == "inh":
        return INHERITED, kind
    if kind in BORROWING_KINDS or kind == "surf":
        return BORROWED, kind
    if kind in ("der", "uder", "from"):
        if language_family(first.lang) in ("", profile.code):
            return INTERNAL, kind
        if first.lang in profile.inherited_sources:
            return INHERITED, f"{kind}:{first.lang}"
        return BORROWED, f"{kind}:{first.lang}"
    if kind in INTERNAL_KINDS:
        surface_first = bool(_SURFACE_START_RE.match(sentences))
        if _INTERNALLY_RE.search(sentences) or (
            patterns.borrowed.match(sentences) and not surface_first
        ):
            return BORROWED, "internal_analysis"
        if surface_first or (
            _SURFACE_RE.search(sentences)
            and not names_native_base(first, profile)
        ):
            return BORROWED, "surface_analysis"
        return INTERNAL, kind
    return UNKNOWN, kind
