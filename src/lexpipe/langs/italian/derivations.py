"""Stage ``derivations``: dorsal-final bases and their derivatives.

Design §5 ("Derivatives"), R21 (candidates in both spellings), R22
(stratum of base and derivative), R26 (morphological suffix targets).

Inputs: ``italian_lexicon_complete.csv`` (``features``: the bases), the
``extract`` projections ``italian_en_lexemes.jsonl.gz`` and
``italian_en_verbs.jsonl.gz`` (etymologies and derived lists), the form
index (every headword of the en split, any part of speech) and the
Leipzig words table (``attest``).

Bases
-----
A base is an in-domain noun or adjective lexeme whose stem ends in a
dorsal target (``dorsal_target``: -co/-ca/-go/-ga, -sco, -cco …; not
the NDE classes). Bases are grouped by root (R4); the stem is the lemma
key without its final vowel (:func:`.morphology.analyze_lemma`).

Links (``link_source``, best first)
-----------------------------------
template            the derivative's en.wiktionary etymology has an
                    ``af``/``affix``/``suffix``/``surf`` step (or
                    ``{{ety|it|:af|…}}``) "base + -suffix" or "prefix- +
                    base + -suffix" (:mod:`.etymology`)
text                the etymology text says "From base + -suffix"
                    without a template
derived_list        the base's entry lists the derivative under
                    "Derived terms"
head_form           the base adjective's headword line lists the
                    superlative (amico: "superlative amicissimo")
candidate_headword  generated from the base stem in **both** spellings
                    (:func:`.suffixes.candidates`, R21) and found among
                    the headwords of the en split
candidate_corpus    generated and attested (≥ 1 token over its corpus
                    forms) in ita_news_2024_1M but not a headword

``link_tier`` is ``named`` for the first four (the dictionary states the
link) and ``candidate`` for the last two (generated, unverified). Stems
shorter than :data:`MIN_CANDIDATE_STEM` letters generate no candidates
(ago → agire, eco → ecista). ``link_status`` is ``accepted`` or one of:

    rejected:named_other_base     the derivative's etymology names a
                                  base that is not a lexeme of the root
                                  (grecizzare names Grecia, not greco:
                                  the it_phon.json fact-check;
                                  ponticello names ponte)
    rejected:etymology_unrelated  a candidate headword whose etymology
                                  is a transmission that does not
                                  involve the base (pagina, from Latin)
    rejected:form_of_other_lemma  a corpus-only candidate one of whose
                                  corpus forms is a form of another
                                  lemma of the en split (accese, a form
                                  of accendere)
    rejected:palatal_base_exists  a palatal candidate whose stem is also
                                  the stem of a frontstem lemma (calcio
                                  → calcista, calcetto). This rule can
                                  only remove palatal candidates (no
                                  other lemma explains a velar ch/gh
                                  before the suffix: a -chio base gives
                                  -chi- derivatives), so it is not
                                  outcome-blind (R16): the analyses
                                  report the candidate tier with these
                                  links and without them
    rejected:capitalized_only     a corpus-only candidate none of whose
                                  tokens is lowercase: a name (R25;
                                  Cecchetto, Legino)

The rejection rules apply to the candidate tier only, except that a
derived-list derivative whose own etymology names another base for the
same suffix is rejected too. A named base is resolved to the dorsal
roots of its lemma; a candidate to the dorsal roots of its stem. With
several roots (antico 'ancient' and àntico 'anterior'), the most
frequent root is chosen (paradigm sum, ties by root_id) and
``base_ambiguous`` says whether the roots' dictionary -i outcomes
differ.

Readings
--------
``edge`` is the stem-edge reading of :func:`.suffixes.read_edge`
(palatal, velar, back, truncated, mismatch); ``palatalized`` is True
for palatal and False for velar. ``deriv_stratum`` comes from the
derivative's own etymology (:func:`.etymology.stratum`; "unknown" when
it has no en entry), ``base_stratum`` from the base root (the first
known stratum of its lexemes, masculine -o lemmas first).

Base descriptors (outcome-blind, for R21's separate reports):
``base_feminine_only`` (every lexeme of the root is a feminine noun),
``base_has_i_cell`` (some lexeme has an -i plural in the dictionary),
``masc_counterpart`` (``other_root``: another dorsal root on the same
stem has an -i cell, musica/musico; ``same_root``: the named base is a
feminine lexeme whose root gets its -i cell from a masculine lexeme,
politica/politico) and ``ic_suffix`` (the base lemma ends in -ico/-ica
whose i is not stressed — ``ico_stressed`` of :mod:`.lexicon`: not
amìco, but arcaico — or a lexeme of its root is an R26 -ico suffix
target). ``base_det_plural`` is the masculine plural of a
determiner or pronoun entry of a base lemma (poco → pochi: the noun is
invariable, the determiner is not), with its stem-edge reading
(``base_det_plural_palatal``) and the thresholds at which it is
attested (``base_det_plural_taus``, :func:`.leipzig.attested_flags`):
inflectional evidence outside the noun/adjective entries, which R21's
no-evidence classes read.

Counts: a derivative's tokens are summed over
:func:`.suffixes.corpus_forms` (folded and case-sensitive); the
``attested_<τ>`` flags attest those forms as one generated paradigm
(:func:`lexpipe.attest.attest_paradigm`: the τ thresholds applied to
the folded sum), except ``attested_dictionary``, which is
True for derivatives listed in the dictionary (headwords, templates,
derived lists, head forms) and False for corpus-only candidates. ``other_spelling``
is the same derivative with the other reading at the junction (leghista
↔ legista), with its headword status and token count.

``fact_check`` annotates the derivatives the it_phon.json verification
corrected (grecizzare, politichese, musicista; porcino, porcello and
porcile are inherited); the coding itself comes from the data.

Outputs (gitignored: Leipzig counts, R9)
---------------------------------------
``data/italian/italian_derivations.csv`` (:data:`DERIVATION_FIELDS`), one
row per (base root, suffix, derivative), sorted by base_root_id, suffix,
derivative; ``data/italian/italian_derivations_summary.tsv``: counts.
"""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

from lexpipe import attest, csvio
from lexpipe.attest import AttestationTable
from lexpipe.context import RunContext, run_stage
from lexpipe.langs.italian import (
    etymology,
    extract,
    features,
    keys,
    leipzig,
    morphology,
    store,
)
from lexpipe.langs.italian import suffixes as sx
from lexpipe.langs.italian.layout import Layout

STAGE: Final[str] = "derivations"
# 5: shared etymology classifier; paradigm attestation from the core
# (lexpipe.attest.attest_paradigm).
# 6: names filter on the τ flags (leipzig.paradigm_flags); the
#    determiner/pronoun plurals of the bases (R21).
VERSION: Final[int] = 7
SOURCE: Final[str] = "wiktionary_en|leipzig"

TEMPLATE: Final[str] = "template"
TEXT: Final[str] = "text"
DERIVED_LIST: Final[str] = "derived_list"
HEAD_FORM: Final[str] = "head_form"
CANDIDATE_HEADWORD: Final[str] = "candidate_headword"
CANDIDATE_CORPUS: Final[str] = "candidate_corpus"
LINK_SOURCES: Final[tuple[str, ...]] = (
    TEMPLATE,
    TEXT,
    DERIVED_LIST,
    HEAD_FORM,
    CANDIDATE_HEADWORD,
    CANDIDATE_CORPUS,
)
NAMED_SOURCES: Final[frozenset[str]] = frozenset(
    {TEMPLATE, TEXT, DERIVED_LIST, HEAD_FORM}
)
DICTIONARY_SOURCES: Final[frozenset[str]] = NAMED_SOURCES | {
    CANDIDATE_HEADWORD
}
TIER_NAMED: Final[str] = "named"
TIER_CANDIDATE: Final[str] = "candidate"
MIN_CANDIDATE_STEM: Final[int] = 3
ACCEPTED: Final[str] = "accepted"
REJECTED_OTHER_BASE: Final[str] = "rejected:named_other_base"
REJECTED_ETYMOLOGY: Final[str] = "rejected:etymology_unrelated"
REJECTED_OTHER_LEMMA: Final[str] = "rejected:form_of_other_lemma"
REJECTED_PALATAL_BASE: Final[str] = "rejected:palatal_base_exists"
REJECTED_CAPITALIZED: Final[str] = "rejected:capitalized_only"
# Frontstem lemma endings whose stem would explain a palatal candidate.
_FRONTSTEM_ENDINGS: Final[tuple[str, ...]] = ("io", "ia", "e", "i")
SAME_ROOT: Final[str] = "same_root"
OTHER_ROOT: Final[str] = "other_root"

# The it_phon.json verification's corrections (annotation only).
FACT_CHECK: Final[Mapping[str, str]] = {
    "grecizzare": (
        "fact-check: from Grecia + -izzare, not greco; not plural-tracking "
        "evidence"
    ),
    "politichese": (
        "fact-check: from feminine politica + -ese; not evidence about a "
        "-ci base"
    ),
    "musicista": (
        "fact-check: from feminine musica + -ista (musiche); masculine "
        "musico/musici exists"
    ),
    "porcino": "fact-check: inherited (Latin porcinus)",
    "porcello": "fact-check: inherited (Late Latin porcellus)",
    "porcile": "fact-check: inherited (Vulgar Latin *porcile)",
}

DERIVATION_FIELDS: Final[tuple[str, ...]] = (
    "base_root_id",
    "base_word",
    "base_lexeme_id",
    "base_pos",
    "base_gender",
    "base_stem",
    "target_letters",
    "base_roots",
    "base_ambiguous",
    "base_outcome_i",
    "base_stratum",
    "base_feminine_only",
    "base_has_i_cell",
    "masc_counterpart",
    "ic_suffix",
    "base_det_plural",
    "base_det_plural_palatal",
    "base_det_plural_taus",
    "derivative",
    "deriv_id",
    "deriv_pos",
    "deriv_headword",
    "suffix",
    "suffix_group",
    "suffix_vowel",
    "spellout_role",
    "edge",
    "palatalized",
    "link_source",
    "link_sources",
    "link_tier",
    "link_status",
    "named_bases",
    "deriv_stratum",
    "deriv_stratum_basis",
    "corpus_forms",
    "count_folded",
    "count_cs",
    *(f"attested_{tau}" for tau in attest.TAUS),
    "other_spelling",
    "other_headword",
    "other_count_folded",
    "fact_check",
    "source",
)

_TEXT_LINK_RE: Final[re.Pattern[str]] = re.compile(
    r"\bFrom ([^\s(),.;:+]+)(?:\s*\([^)]*\))?\s*\+\s*(-[^\s(),.;:]+)"
)


# ---------------------------------------------------------------------------
# Bases
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class BaseLexeme:
    """One dorsal noun/adjective lexeme (lexicon_complete columns)."""

    lexeme_id: str
    root_id: str
    word: str
    key: str
    pos: str
    gender: str
    stem: str
    target_letters: str
    stratum: str
    has_i_cell: bool
    outcome_i: str
    stress_from_end: str
    suffix_target: bool
    ico_stressed: str = ""


@dataclass
class BaseRoot:
    """The dorsal lexemes of one root."""

    root_id: str
    stem: str
    frequency: int
    lexemes: list[BaseLexeme] = field(default_factory=list)

    @property
    def keys(self) -> frozenset[str]:
        return frozenset(lex.key for lex in self.lexemes)

    @property
    def has_i_cell(self) -> bool:
        return any(lex.has_i_cell for lex in self.lexemes)

    @property
    def feminine_only(self) -> bool:
        return all(
            lex.pos == "noun" and lex.gender == "f" for lex in self.lexemes
        )

    @property
    def outcome_i(self) -> str:
        """Dictionary -i outcome of the root (palatalized, velar,
        variable or none)."""
        found = {lex.outcome_i for lex in self.lexemes if lex.outcome_i}
        found.discard(features.NONE)
        if not found:
            return features.NONE
        if len(found) == 1:
            return next(iter(found))
        return features.VARIABLE

    def primary(self) -> BaseLexeme:
        """Masculine -o lemmas first, then nouns, then by lexeme_id."""
        return min(
            self.lexemes,
            key=lambda lex: (
                not lex.key.endswith("o"),
                lex.pos != "noun",
                lex.lexeme_id,
            ),
        )

    @property
    def stratum(self) -> str:
        ordered = sorted(
            self.lexemes,
            key=lambda lex: (not lex.key.endswith("o"), lex.lexeme_id),
        )
        for lex in ordered:
            if lex.stratum != etymology.UNKNOWN:
                return lex.stratum
        return etymology.UNKNOWN

    @property
    def ic_suffix(self) -> bool:
        """An R26 target among the lexemes, or a lemma in -ico/-ica whose
        i is not stressed (not amìco; arcaico counts)."""
        return any(
            lex.suffix_target
            or (
                lex.key.endswith(("ico", "ica")) and lex.ico_stressed != "True"
            )
            for lex in self.lexemes
        )


@dataclass
class BaseIndex:
    """Dorsal bases by lemma key, root and stem."""

    by_key: dict[str, list[BaseLexeme]]
    roots: dict[str, BaseRoot]
    by_stem: dict[str, list[str]]


def _bool(text: str) -> bool:
    return csvio.parse_bool(text) is True


def load_bases(path: Path) -> BaseIndex:
    """Dorsal in-domain noun/adjective lexemes of lexicon_complete."""
    by_key: dict[str, list[BaseLexeme]] = defaultdict(list)
    roots: dict[str, BaseRoot] = {}
    for row in csvio.read_rows(path):
        if not (_bool(row["in_domain"]) and _bool(row["dorsal_target"])):
            continue
        key = keys.match_key(row["word"])
        stem = morphology.analyze_lemma(key)
        lexeme = BaseLexeme(
            lexeme_id=row["lexeme_id"],
            root_id=row["root_id"],
            word=row["word"],
            key=key,
            pos=row["pos"],
            gender=row["gender"],
            stem=stem.stem,
            target_letters=stem.target_letters,
            stratum=row["stratum"] or etymology.UNKNOWN,
            has_i_cell=bool(row["forms_i"]),
            outcome_i=row["outcome_i"],
            stress_from_end=row["stress_from_end"],
            ico_stressed=row["ico_stressed"],
            suffix_target=_bool(row["target_is_suffix"]),
        )
        by_key[key].append(lexeme)
        root = roots.get(row["root_id"])
        if root is None:
            root = BaseRoot(
                row["root_id"], stem.stem, int(row["freq_root_folded"] or 0)
            )
            roots[row["root_id"]] = root
        root.lexemes.append(lexeme)
    by_stem: dict[str, list[str]] = defaultdict(list)
    for root_id in sorted(roots):
        by_stem[roots[root_id].stem].append(root_id)
    return BaseIndex(dict(by_key), roots, dict(by_stem))


def choose_root(index: BaseIndex, root_ids: Iterable[str]) -> str:
    """The most frequent root (ties by root_id)."""
    ordered = sorted(
        set(root_ids), key=lambda r: (-index.roots[r].frequency, r)
    )
    if not ordered:
        raise ValueError("no root to choose from")
    return ordered[0]


# ---------------------------------------------------------------------------
# Headwords and named links
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Headword:
    """An en noun, adjective or verb entry, reduced to what links need.

    ``named`` holds the (base key, suffix) pairs of inventory suffixes;
    ``named_any`` every base key an affix step names (any affix).
    """

    ident: str
    word: str
    key: str
    pos: str
    stratum: str
    basis: str
    named: tuple[tuple[str, str], ...]
    named_any: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class RawLink:
    """A derivative linked to a base lemma or stem before resolution."""

    derivative: str
    key: str
    suffix: sx.Suffix
    source: str
    base_key: str = ""
    stem: str = ""


def _plain_term(term: str) -> bool:
    return bool(term) and not term.startswith("-") and not term.endswith("-")


def affix_bases(steps: Sequence[etymology.Step]) -> list[str]:
    """Every plain term (no leading or trailing hyphen) of an affix or
    surface-analysis step: the bases an etymology names."""
    return [
        term
        for step in steps
        if step.kind in ("af", "surf")
        for term in step.terms
        if _plain_term(term)
    ]


def named_links(
    steps: Sequence[etymology.Step], text: str
) -> list[tuple[str, sx.Suffix, str]]:
    """``(base word, suffix, source)`` the etymology names (module
    docstring): 2-part and prefix + base + suffix affix steps, else the
    "From X + -suffix" text."""
    out: list[tuple[str, sx.Suffix, str]] = []
    for step in steps:
        if step.kind not in ("af", "surf"):
            continue
        terms = [term for term in step.terms if term]
        base = ""
        if len(terms) == 2 and _plain_term(terms[0]):
            base = terms[0]
        elif (
            len(terms) == 3
            and terms[0].endswith("-")
            and _plain_term(terms[1])
        ):
            base = terms[1]
        suffix = sx.canonical(terms[-1]) if base else None
        if suffix is not None:
            out.append((base, suffix, TEMPLATE))
    if not out:
        for match in _TEXT_LINK_RE.finditer(text or ""):
            suffix = sx.canonical(match.group(2))
            if suffix is not None:
                out.append((match.group(1), suffix, TEXT))
    seen: dict[tuple[str, str], tuple[str, sx.Suffix, str]] = {}
    for base, suffix, source in out:
        seen.setdefault(
            (keys.match_key(base), suffix.name), (base, suffix, source)
        )
    return list(seen.values())


@dataclass
class Scan:
    """Everything one pass over the projections yields."""

    headwords: dict[str, list[Headword]]
    links: list[RawLink]


def _strings(value: Any) -> list[str]:
    if not isinstance(value, list):
        return []
    return [item for item in value if isinstance(item, str)]


def _scan_record(
    record: Mapping[str, Any],
    ident_field: str,
    bases: BaseIndex,
    scan: Scan,
) -> None:
    word = str(record.get("word") or "")
    key = keys.match_key(word)
    if not key:
        return
    steps = etymology.steps(record.get("etymology_templates") or ())
    text = str(record.get("etymology_text") or "")
    stratum, basis = etymology.stratum(steps, text, word=word)
    named = named_links(steps, text)
    scan.headwords.setdefault(key, []).append(
        Headword(
            ident=str(record[ident_field]),
            word=word,
            key=key,
            pos=str(record.get("pos") or ""),
            stratum=stratum,
            basis=basis,
            named=tuple(
                sorted((keys.match_key(b), s.name) for b, s, _ in named)
            ),
            named_any=tuple(
                sorted(
                    {keys.match_key(b) for b in affix_bases(steps)}
                    | {keys.match_key(b) for b, _, _ in named}
                )
            ),
        )
    )
    for base, suffix, source in named:
        base_key = keys.match_key(base)
        if base_key in bases.by_key:
            scan.links.append(
                RawLink(word, key, suffix, source, base_key=base_key)
            )
    if key in bases.by_key:
        for derived in _strings(record.get("derived")):
            _derived_link(derived, key, bases, scan)


def _derived_link(
    derived: str, base_key: str, bases: BaseIndex, scan: Scan
) -> None:
    derived_key = keys.match_key(derived)
    if not derived_key.isalpha():
        return
    stem = bases.by_key[base_key][0].stem
    reading = sx.best_reading(stem, derived_key, allow_prefix=True)
    if reading is not None:
        scan.links.append(
            RawLink(
                derived,
                derived_key,
                reading[0],
                DERIVED_LIST,
                base_key=base_key,
            )
        )


def scan_projections(
    lexemes_path: Path, verbs_path: Path, bases: BaseIndex
) -> Scan:
    """One pass over each projection (noun/adjective, verb)."""
    scan = Scan({}, [])
    for path, fmt, ident in (
        (lexemes_path, extract.EN_LEXEMES_FORMAT, "lexeme_id"),
        (verbs_path, extract.EN_VERBS_FORMAT, "verb_id"),
    ):
        meta, records = store.read_jsonl(path)
        if meta.get("format") != fmt:
            raise store.StaleDerivativeError(f"{path}: unexpected {meta}")
        for record in records:
            _scan_record(record, ident, bases, scan)
    for entries in scan.headwords.values():
        entries.sort(key=lambda h: h.ident)
    return scan


# ---------------------------------------------------------------------------
# The form index
# ---------------------------------------------------------------------------


@dataclass
class FormIndexScan:
    """Headword keys of every part of speech, the superlatives the
    adjective heads of dorsal bases list, and the masculine plurals of
    determiner/pronoun entries of base lemmas (pochi of poco)."""

    lemma_keys: frozenset[str]
    head_links: list[RawLink]
    det_plurals: Mapping[str, str] = field(default_factory=dict)


# Parts of speech whose plural is inflectional evidence for a noun or
# adjective base of the same lemma (R21: poco N is invariable, pochi DET).
DET_POS: Final[frozenset[str]] = frozenset({"det", "pron"})


def _index_rows(path: Path) -> Iterator[dict[str, str]]:
    meta, fields, rows = store.read_tsv(path)
    if meta.get("format") != extract.FORM_INDEX_FORMAT:
        raise store.StaleDerivativeError(f"{path}: unexpected {meta}")
    for row in rows:
        yield dict(zip(fields, row, strict=False))


def scan_form_index(path: Path, bases: BaseIndex) -> FormIndexScan:
    """One pass over the form index (module docstring)."""
    superlative = sx.BY_NAME["-issimo"]
    lemma_keys: set[str] = set()
    links: list[RawLink] = []
    det_plurals: dict[str, set[str]] = defaultdict(set)
    for row in _index_rows(path):
        if row["source"] == "lemma":
            lemma_keys.add(row["key"])
            continue
        lemma_key = keys.match_key(row["lemma"])
        tags = set(row["tags"].split("|"))
        if (
            row["pos"] in DET_POS
            and "plural" in tags
            and "feminine" not in tags
            and lemma_key in bases.by_key
            and row["key"].isalpha()
        ):
            det_plurals[lemma_key].add(row["key"])
        if (
            row["source"] == "head"
            and row["pos"] == "adj"
            and "superlative" in row["tags"].split("|")
            and lemma_key in bases.by_key
            and row["key"].isalpha()
        ):
            links.append(
                RawLink(
                    row["form"],
                    row["key"],
                    superlative,
                    HEAD_FORM,
                    base_key=lemma_key,
                )
            )
    return FormIndexScan(
        frozenset(lemma_keys),
        links,
        {key: min(forms) for key, forms in det_plurals.items()},
    )


def lemmas_of_forms(path: Path, wanted: frozenset[str]) -> dict[str, set[str]]:
    """Lemma keys of every form-index row whose key is in *wanted*."""
    out: dict[str, set[str]] = defaultdict(set)
    for row in _index_rows(path):
        if row["key"] in wanted:
            out[row["key"]].add(keys.match_key(row["lemma"]))
    return dict(out)


# ---------------------------------------------------------------------------
# Candidates (R21)
# ---------------------------------------------------------------------------


def corpus_count(
    table: AttestationTable, forms: Iterable[str], case_sensitive: bool
) -> int:
    """Tokens summed over distinct *forms*."""
    return sum(
        table.count(form, case_sensitive=case_sensitive)
        for form in dict.fromkeys(forms)
    )


def candidate_links(
    bases: BaseIndex,
    lemma_keys: frozenset[str],
    table: AttestationTable,
) -> Iterator[RawLink]:
    """Generated derivatives found among the headwords or in the corpus
    (stems of at least :data:`MIN_CANDIDATE_STEM` letters)."""
    for stem in sorted(bases.by_stem):
        if len(stem) < MIN_CANDIDATE_STEM:
            continue
        for suffix in sx.SUFFIXES:
            for candidate in sx.candidates(stem, suffix):
                if candidate.form in lemma_keys:
                    source = CANDIDATE_HEADWORD
                elif corpus_count(
                    table, sx.corpus_forms(candidate.form, suffix), False
                ):
                    source = CANDIDATE_CORPUS
                else:
                    continue
                yield RawLink(
                    candidate.form, candidate.form, suffix, source, stem=stem
                )


# ---------------------------------------------------------------------------
# Resolution
# ---------------------------------------------------------------------------


@dataclass
class Resolved:
    """A link resolved to a base root (one output row)."""

    link: RawLink
    root_id: str
    roots: tuple[str, ...]
    base_lexeme: BaseLexeme
    sources: set[str] = field(default_factory=set)


def _roots_of(link: RawLink, bases: BaseIndex) -> tuple[str, ...]:
    if link.base_key:
        found = {lex.root_id for lex in bases.by_key.get(link.base_key, ())}
    else:
        found = set(bases.by_stem.get(link.stem, ()))
    return tuple(sorted(found))


def _base_lexeme(link: RawLink, root: BaseRoot) -> BaseLexeme:
    """The named lexeme (one with an -i plural first: greco the person,
    not greco the language), else the root's primary lexeme."""
    if link.base_key:
        named = [lex for lex in root.lexemes if lex.key == link.base_key]
        if named:
            return min(
                named,
                key=lambda lex: (
                    not lex.has_i_cell,
                    lex.pos != "noun",
                    lex.lexeme_id,
                ),
            )
    return root.primary()


def resolve(links: Iterable[RawLink], bases: BaseIndex) -> list[Resolved]:
    """Group raw links by (derivative, suffix, chosen root), keeping the
    best link source (module docstring)."""
    merged: dict[tuple[str, str, str], Resolved] = {}
    rank = {source: i for i, source in enumerate(LINK_SOURCES)}
    for link in links:
        roots = _roots_of(link, bases)
        if not roots:
            continue
        root_id = choose_root(bases, roots)
        group = (link.key, link.suffix.name, root_id)
        current = merged.get(group)
        if current is None:
            current = Resolved(
                link, root_id, roots, _base_lexeme(link, bases.roots[root_id])
            )
            merged[group] = current
        elif rank[link.source] < rank[current.link.source]:
            current.link = link
            current.base_lexeme = _base_lexeme(link, bases.roots[root_id])
        current.sources.add(link.source)
    return [merged[group] for group in sorted(merged)]


# ---------------------------------------------------------------------------
# Rows
# ---------------------------------------------------------------------------


@dataclass
class Context:
    """Shared inputs of the row builder."""

    bases: BaseIndex
    headwords: Mapping[str, Sequence[Headword]]
    lemma_keys: frozenset[str]
    table: AttestationTable
    form_lemmas: Mapping[str, set[str]] = field(default_factory=dict)
    sums: Mapping[str, Mapping[str, bool]] | None = None
    det_plurals: Mapping[str, str] = field(default_factory=dict)


def det_plural_columns(root: BaseRoot, ctx: Context) -> dict[str, Any]:
    """``base_det_plural`` and its reading and thresholds (module
    docstring) for the root's first lexeme (by key) that has one."""
    for lexeme in sorted(root.lexemes, key=lambda lex: lex.key):
        form = ctx.det_plurals.get(lexeme.key)
        if not form:
            continue
        stem = morphology.analyze_lemma(lexeme.key)
        reading = morphology.analyze_form(lexeme.key, stem, form)
        flags = leipzig.attested_flags(ctx.table, form)
        return {
            "base_det_plural": form,
            "base_det_plural_palatal": reading.mutation,
            "base_det_plural_taus": [t for t in attest.TAUS if flags[t]],
        }
    return {
        "base_det_plural": "",
        "base_det_plural_palatal": None,
        "base_det_plural_taus": [],
    }


def _headword_for(
    key: str, suffix: sx.Suffix, base_key: str, ctx: Context
) -> Headword | None:
    entries = list(ctx.headwords.get(key, ()))
    if not entries:
        return None
    return min(
        entries,
        key=lambda h: (
            not any(n == (base_key, suffix.name) for n in h.named),
            h.pos != suffix.pos,
            h.ident,
        ),
    )


def named_bases(key: str, suffix: sx.Suffix, ctx: Context) -> list[str]:
    """Base keys the derivative's entries name for *suffix* (sorted)."""
    return sorted(
        {
            base
            for entry in ctx.headwords.get(key, ())
            for base, name in entry.named
            if name == suffix.name
        }
    )


def _palatal_base_exists(item: Resolved, edge: sx.Edge, ctx: Context) -> bool:
    if edge != sx.Edge.PALATAL:
        return False
    stem = item.base_lexeme.stem
    return any(stem + end in ctx.lemma_keys for end in _FRONTSTEM_ENDINGS)


def _candidate_status(item: Resolved, edge: sx.Edge, ctx: Context) -> str:
    root = ctx.bases.roots[item.root_id]
    link = item.link
    entries = ctx.headwords.get(link.key, ())
    if link.source == CANDIDATE_HEADWORD and entries:
        mentioned = {base for entry in entries for base in entry.named_any}
        if mentioned and root.keys.isdisjoint(mentioned):
            return REJECTED_OTHER_BASE
        if not mentioned and any(
            entry.stratum != etymology.UNKNOWN for entry in entries
        ):
            return REJECTED_ETYMOLOGY
    if link.source == CANDIDATE_CORPUS:
        forms = sx.corpus_forms(link.key, link.suffix)
        allowed = root.keys | {link.key}
        for form in forms:
            if ctx.form_lemmas.get(form, set()) - allowed:
                return REJECTED_OTHER_LEMMA
        if not corpus_count(ctx.table, forms, True):
            return REJECTED_CAPITALIZED
    if _palatal_base_exists(item, edge, ctx):
        return REJECTED_PALATAL_BASE
    return ACCEPTED


def _status(item: Resolved, edge: sx.Edge, ctx: Context) -> str:
    """``link_status`` (module docstring)."""
    source = item.link.source
    if source in (TEMPLATE, TEXT, HEAD_FORM):
        return ACCEPTED
    if source == DERIVED_LIST:
        named = named_bases(item.link.key, item.link.suffix, ctx)
        root = ctx.bases.roots[item.root_id]
        if named and root.keys.isdisjoint(named):
            return REJECTED_OTHER_BASE
        return ACCEPTED
    return _candidate_status(item, edge, ctx)


def _masc_counterpart(item: Resolved, ctx: Context) -> str:
    root = ctx.bases.roots[item.root_id]
    base = item.base_lexeme
    if any(lex.has_i_cell for lex in root.lexemes if lex.key == base.key):
        return ""
    if root.has_i_cell:
        return SAME_ROOT
    others = [
        ctx.bases.roots[r]
        for r in ctx.bases.by_stem.get(root.stem, ())
        if r != item.root_id
    ]
    return OTHER_ROOT if any(o.has_i_cell for o in others) else ""


def _edge(item: Resolved) -> tuple[sx.Edge, str]:
    link = item.link
    stem = item.base_lexeme.stem
    allow_prefix = link.source in NAMED_SOURCES
    return sx.read_edge(stem, link.suffix, link.key, allow_prefix=allow_prefix)


def _palatalized(edge: sx.Edge) -> bool | None:
    if edge == sx.Edge.PALATAL:
        return True
    if edge == sx.Edge.VELAR:
        return False
    return None


def build_row(item: Resolved, ctx: Context) -> dict[str, Any]:
    """One :data:`DERIVATION_FIELDS` row."""
    link = item.link
    suffix = link.suffix
    root = ctx.bases.roots[item.root_id]
    base = item.base_lexeme
    edge, body = _edge(item)
    headword = _headword_for(link.key, suffix, base.key, ctx)
    forms = sx.corpus_forms(link.key, suffix)
    other = sx.other_spelling(link.key, edge, body)
    outcomes = {ctx.bases.roots[r].outcome_i for r in item.roots}
    attested = (
        dict(ctx.sums.get(link.key, dict.fromkeys(attest.TAUS, False)))
        if ctx.sums is not None
        else dict.fromkeys(attest.TAUS, False)
    )
    attested["dictionary"] = link.source in DICTIONARY_SOURCES or bool(
        headword
    )
    row: dict[str, Any] = {
        "base_root_id": item.root_id,
        "base_word": base.word,
        "base_lexeme_id": base.lexeme_id,
        "base_pos": base.pos,
        "base_gender": base.gender,
        "base_stem": base.stem,
        "target_letters": base.target_letters,
        "base_roots": list(item.roots),
        "base_ambiguous": len(outcomes) > 1,
        "base_outcome_i": root.outcome_i,
        "base_stratum": root.stratum,
        "base_feminine_only": root.feminine_only,
        "base_has_i_cell": root.has_i_cell,
        "masc_counterpart": _masc_counterpart(item, ctx),
        "ic_suffix": root.ic_suffix,
        **det_plural_columns(root, ctx),
        "derivative": headword.word if headword else link.derivative,
        "deriv_id": headword.ident if headword else "",
        "deriv_pos": headword.pos if headword else suffix.pos,
        "deriv_headword": link.key in ctx.lemma_keys,
        "suffix": suffix.name,
        "suffix_group": suffix.group,
        "suffix_vowel": suffix.vowel,
        "spellout_role": suffix.spellout,
        "edge": edge.value,
        "palatalized": _palatalized(edge),
        "link_source": link.source,
        "link_sources": [s for s in LINK_SOURCES if s in item.sources],
        "link_tier": (
            TIER_NAMED if link.source in NAMED_SOURCES else TIER_CANDIDATE
        ),
        "link_status": _status(item, edge, ctx),
        "named_bases": named_bases(link.key, suffix, ctx),
        "deriv_stratum": headword.stratum if headword else etymology.UNKNOWN,
        "deriv_stratum_basis": headword.basis if headword else "no_entry",
        "corpus_forms": forms,
        "count_folded": corpus_count(ctx.table, forms, False),
        "count_cs": corpus_count(ctx.table, forms, True),
        "other_spelling": other,
        "other_headword": bool(other) and other in ctx.lemma_keys,
        "other_count_folded": (
            corpus_count(ctx.table, sx.corpus_forms(other, suffix), False)
            if other
            else None
        ),
        "fact_check": FACT_CHECK.get(link.key, ""),
        "source": SOURCE,
    }
    for tau in attest.TAUS:
        row[f"attested_{tau}"] = attested[tau]
    return row


def attested_sums(
    resolved: Sequence[Resolved], table: AttestationTable
) -> dict[str, Mapping[str, bool]]:
    """``{derivative key: {τ: attested}}``: the derivative's corpus forms
    attested as one generated paradigm (:func:`.leipzig.paradigm_flags`:
    folded counts summed over distinct forms, with a token in its own
    case among them, R25)."""
    return {
        item.link.key: leipzig.paradigm_flags(
            table, sx.corpus_forms(item.link.key, item.link.suffix)
        )
        for item in resolved
    }


def summary(rows: Sequence[Mapping[str, Any]], scan: Scan) -> str:
    """Counts by link source, status, suffix and edge."""
    counts: Counter[str] = Counter()
    counts["headword_keys"] = len(scan.headwords)
    counts["raw_links"] = len(scan.links)
    for row in rows:
        counts["rows"] += 1
        counts[f"link_source:{row['link_source']}"] += 1
        counts[f"link_status:{row['link_status']}"] += 1
        counts[f"link_tier:{row['link_tier']}"] += 1
        counts[f"edge:{row['edge']}"] += 1
        if row["link_status"] != ACCEPTED:
            continue
        counts[
            f"accepted:{row['link_tier']}:{row['suffix']}:{row['edge']}"
        ] += 1
        counts[f"accepted_deriv_stratum:{row['deriv_stratum']}"] += 1
    lines = ["key\tvalue"]
    lines.extend(f"{key}\t{counts[key]}" for key in sorted(counts))
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# The stage
# ---------------------------------------------------------------------------


def _candidate_forms(resolved: Sequence[Resolved]) -> frozenset[str]:
    return frozenset(
        form
        for item in resolved
        if item.link.source == CANDIDATE_CORPUS
        for form in sx.corpus_forms(item.link.key, item.link.suffix)
    )


def _body(ctx: RunContext, layout: Layout, words: Path) -> int:
    bases = load_bases(layout.lexicon_complete)
    scan = scan_projections(layout.en_lexemes, layout.en_verbs, bases)
    index = scan_form_index(layout.form_index, bases)
    table = leipzig.attestation_table(ctx, words)
    links = list(scan.links) + index.head_links
    links.extend(candidate_links(bases, index.lemma_keys, table))
    resolved = resolve(links, bases)
    shared = Context(
        bases,
        scan.headwords,
        index.lemma_keys,
        table,
        det_plurals=index.det_plurals,
    )
    shared.form_lemmas = lemmas_of_forms(
        layout.form_index, _candidate_forms(resolved)
    )
    shared.sums = attested_sums(resolved, table)
    rows = [build_row(item, shared) for item in resolved]
    rows.sort(
        key=lambda r: (
            r["base_root_id"],
            r["suffix"],
            keys.match_key(str(r["derivative"])),
            r["derivative"],
        )
    )
    csvio.write_rows(layout.derivations, rows, DERIVATION_FIELDS, ctx=ctx)
    ctx.write_text(layout.derivations_summary, summary(rows, scan))
    accepted = Counter(
        str(r["link_tier"])
        for r in rows
        if r["link_status"] == ACCEPTED and r["edge"] in ("palatal", "velar")
    )
    print(
        f"  {len(bases.roots):,} dorsal roots; {len(rows):,} derivative "
        f"links; accepted with a palatal/velar edge: "
        f"{accepted[TIER_NAMED]:,} named, "
        f"{accepted[TIER_CANDIDATE]:,} candidates"
    )
    return 0


def outputs(layout: Layout) -> list[Path]:
    return [layout.derivations, layout.derivations_summary]


def run(ctx: RunContext) -> int | None:
    """Run the ``derivations`` stage (skipped when its manifest matches)."""
    layout = Layout(ctx)
    words = leipzig.words_file(ctx, layout)
    return run_stage(
        ctx,
        STAGE,
        inputs=[
            layout.lexicon_complete,
            layout.en_lexemes,
            layout.en_verbs,
            layout.form_index,
            words,
        ],
        params={
            "version": VERSION,
            "suffixes": [s.name for s in sx.SUFFIXES],
            "min_candidate_stem": MIN_CANDIDATE_STEM,
            "elision_clitics": sorted(ctx.language.elision_clitics),
        },
        outputs=outputs(layout),
        body=lambda: _body(ctx, layout, words),
    )
