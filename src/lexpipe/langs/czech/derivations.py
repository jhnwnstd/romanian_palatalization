"""Stage ``derivations``: base × derivative pairs (design §6; R14, R20–R22).

Every in-domain noun of the lexicon is a possible **base**, and so are
the dual remnants (ruka, oko), whose singular is regular. A pair
(base, derivative) is found by any of these kinds of evidence,
recorded in ``evidence`` (sorted list):

``wikt_ety``      the derivative's English Wiktionary etymology states
                  base + suffix (``{{ety|cs|:af|kniha|-ní}}``,
                  :func:`lexpipe.langs.czech.etymology.affix_parts`);
``wikt_head``     the base's headword line names it (adj=, dim=, f=, m=);
``wikt_derived``  the base's "Derived terms" list;
``wikt_related``  the base's "Related terms" list (noisy, so only words
                  that decompose against the base are kept);
``derinet``       a DeriNet 2.3 ``Type=Derivation`` edge whose parent is
                  the base (R14: joined on lemma, gender and animacy;
                  Fictitious nodes dropped; projection in
                  :mod:`.derinet`);
``cand_headword`` a generated candidate (R21, :func:`.suffixes.candidates`)
                  that is an English or Czech Wiktionary headword of the
                  suffix's part of speech;
``cand_corpus``   a generated candidate that is no headword at all but
                  has corpus tokens of its own.

Candidates (R21). Both spellings are generated for P and NT suffixes
(vzdušný/vzduchný, tankista/tancista), only the plain one for N
suffixes (a mutated -ový/-ovat/-yně/-ův spelling is never a Czech
derivative: its hits are other words, pracovat ← práce, not prak). A
candidate is kept only if (:func:`candidate_problem`) its spelling is
the expected one (the first palatalization for a P suffix that
palatalizes the base's target: drozní is drozdí, not from droga; plain
otherwise), its stem is not shared by bases with other lemmas (deka /
deko → dečka), its own Wiktionary etymology does not tell another
story (inherited, borrowed, or formed from another base), an -ista or
-ismus candidate is not the derivative of an -ie noun of the same stem
(ekologista ← ekologie, monarchista ← monarchie), and, for a
corpus-only candidate, its stem has at least :data:`MIN_CORPUS_STEM`
letters, its lemma form occurs and its own forms have at least
:data:`MIN_CORPUS_TOKENS` tokens. The spelling check would make a P
suffix's outcome ``pal`` by construction (R16: no outcome-conditioned
acceptance), so a **headword** candidate in the plain spelling of a P
suffix that palatalizes the base (pantofágní) is accepted as
``cand_headword_plain`` evidence (basis ``wikt``) when it passes the
other checks; the spelling check stays for corpus-only candidates and
for mutated spellings of the wrong palatalization (drozní). An
unexpected spelling also corroborates a pair a Wiktionary link already
found (barokní).

Only pairs whose derivative decomposes into (prefix +) base stem
(plain or mutated) + a suffix of :data:`.suffixes.SUFFIXES` are kept;
undecomposable links are counted in the report. ``basis`` is the
narrowest evidence level of a pair: ``linked`` (a Wiktionary link)
⊂ ``wikt`` (+ Wiktionary headword candidates; the main basis of the
analyses) ⊂ ``wikt+corpus`` (+ corpus-only candidates) ⊂ ``all``
(+ DeriNet-only pairs: NC). Tracked tables use the first three.

Which base a link belongs to (audit 2026-10-23):

- **Homographs.** Derived and related terms belong to the entry that
  lists them (keyed by its lexeme id), not to every noun with that
  lemma (drak 'dragon' lists dračí, drak 'glider' does not). A
  derivative found for several homograph roots of one lemma (an
  etymology template or a candidate names the lemma only; lemmas
  compared case-insensitively: Buddha, buddha) is kept on one root
  (:func:`dedupe_homographs`): a feminine agent (-yně, or -ka glossed
  "female …", "woman", or a "female equivalent of" sense) on an
  animate masculine base, then the base whose own entry links it, then
  the one with more evidence kinds, then the first lexeme_id.
- **Related terms** are noisy: a related-terms link is dropped
  (``wikt_related`` evidence removed) when the derivative's own
  etymology names another base (rozkročit is roz- + krok + -it,
  vystrašit vy- + strašit), when the derivative is borrowed whole
  (papoušek ← papouch), or when it is an inherited verb: an inherited
  noun and an inherited verb are related, the direction unknown
  (křik – křičet, lék – léčit, dluh – dlužit). An inherited adjective
  or noun of the base keeps its link (strach – strašný).
- **Direction.** A base whose own etymology derives it from a verb
  (``deverbal``/``back-formation`` templates, "Deverbal from X",
  "Formed from X", "From the verb X") records those verbs in
  ``base_deverbal_from``. A verb derivative of that verb family
  (obloha → obložit, útok → útočit, leh → ležet ← lehat, tucha →
  tušit) is the base's source, not its derivative: the pair is
  dropped (``reversed``).
- **Names.** A derivative that is only a proper name (Sněžka,
  Kameníček) is dropped (R25).

``derivative_feminine_agent`` marks the feminizing derivatives that
the -ka vs -yně pair compares (R20, R23 use it with animate
masculine bases only).

Each pair records the suffix, its class (P / N / NT), the competing
pair it belongs to (R20, R23), the outcome (``pal`` / ``plain`` /
``exception`` / ``nontarget``, :class:`.suffixes.Decomposition`), the
DeriNet verbalizer of verbs and whether ConjugClass agrees (R14), the
strata of base and derivative (R22, :mod:`.etymology`), and the
derivative's corpus attestation: the summed tokens of its forms
(Wiktionary paradigm if it is a headword, else
:func:`.suffixes.paradigm`), each form counted case-sensitively plus
its sentence-initial capitalized spelling, and forms that the form
index assigns to another entry left out (klukovi is kluk's DAT.SG, not
a form of klukův). τ flags use the thresholds of :mod:`lexpipe.attest`.

Outputs
-------
``data/czech/nc/czech_derinet.sqlite``      DeriNet projection (NC)
``data/czech/nc/czech_derivations.csv``     every pair (NC: DeriNet,
                                            czetyl and per-word Leipzig
                                            counts)
``data/czech/nc/czech_lexeme_strata.csv``   R22 strata of every noun and
                                            adjective (Wiktionary-only and
                                            all-source columns)
``analysis/czech/czech_derivations_summary.csv``  tracked aggregates of
                                            the ``linked``, ``wikt`` and
                                            ``wikt+corpus`` bases, through
                                            the licence gate
``analysis/czech/czech_derivations_report.txt``   tracked aggregate
                                            report (DeriNet appears as
                                            counts only)
"""

from __future__ import annotations

import contextlib
import re
import sqlite3
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Final

from lexpipe import attest, csvio
from lexpipe import etymology as core_etymology
from lexpipe import ids, kaikki
from lexpipe.attest import TAUS, AttestationTable
from lexpipe.context import RunContext, run_stage
from lexpipe.langs.czech import derinet, etymology, extract, jsonio, leipzig
from lexpipe.langs.czech import morphology as morph
from lexpipe.langs.czech import suffixes as sx
from lexpipe.langs.czech.analysis_layout import AnalysisPaths, analysis_paths
from lexpipe.langs.czech.attestation import (
    CORPUS_TAUS,
    base_id,
    counts_as_owner,
)
from lexpipe.langs.czech.layout import CzechPaths, paths_for

STAGE: Final[str] = "derivations"
# 5: shared etymology classifier (lexpipe.etymology); core paradigm
# attestation (lexpipe.attest.attest_paradigm); licence-gated report.
# 6: entry-keyed links, homograph dedupe, reversed (deverbal) pairs,
# related-terms pruning, name derivatives, feminine agents, plain
# headword candidates, -icismus and -ie checks; Czech strata refined
# (adapted loan vs hybrid; audit 2026-10-23).
VERSION: Final[str] = "6"
MIN_STEM: Final[int] = 3
MIN_CORPUS_STEM: Final[int] = 4
MIN_CORPUS_TOKENS: Final[int] = 2
# Target classes whose bases get generated candidates: the visible
# targets, the inferred coronals and the non-target baseline (R23).
CANDIDATE_CLASSES: Final[frozenset[str]] = frozenset(
    {
        "velar",
        "r",
        "c",
        "coronal_inferred",
        "labial",
        "lateral",
        "sibilant",
    }
)
LINK_EVIDENCE: Final[frozenset[str]] = frozenset(
    {"wikt_ety", "wikt_head", "wikt_derived", "wikt_related"}
)
# Evidence that comes from the base entry itself (not from a lemma
# match): the homograph dedupe prefers it.
OWN_EVIDENCE: Final[frozenset[str]] = frozenset(
    {"wikt_head", "wikt_derived", "wikt_related", "derinet"}
)
HEADWORD_CANDIDATES: Final[frozenset[str]] = frozenset(
    {"cand_headword", "cand_headword_plain"}
)
VERB_SUFFIXES: Final[frozenset[str]] = frozenset({"-it", "-et", "-ovat"})
FEMININE_AGENT_SUFFIXES: Final[frozenset[str]] = frozenset({"-ka", "-yně"})
_FEMALE_GLOSS: Final[re.Pattern[str]] = re.compile(
    r"\b(?:female|woman|women|girl|feminine)\b", re.IGNORECASE
)
# Etymology sentences that derive a noun from a verb (module docstring).
_DEVERBAL_TEXT: Final[re.Pattern[str]] = re.compile(
    r"(?:Deverbal from|Back-formation from|Backformation from|"
    r"Formed from|From the verb)\s+([^\W\d_]+)(?:[^.]*?\bof\s+"
    r"([^\W\d_]+))?"
)
_DEVERBAL_KINDS: Final[frozenset[str]] = frozenset(
    {"deverbal", "backform", "back-formation", "back-form"}
)
BASES: Final[tuple[str, ...]] = ("linked", "wikt", "wikt+corpus", "all")
SOURCE_EN: Final[str] = "wiktionary_en"
SOURCE_CS: Final[str] = "wiktionary_cs"
TRACKED_SOURCE: Final[list[str]] = [SOURCE_EN, SOURCE_CS]
# The tracked report counts DeriNet pairs and czetyl/DeriNet strata:
# aggregates of CC BY-NC-SA sources (R9; DATA_LICENSE.md draft policy).
REPORT_SOURCES: Final[tuple[str, ...]] = (
    SOURCE_EN,
    SOURCE_CS,
    "derinet",
    "czetyl",
)
# Bases whose aggregates may be tracked: no DeriNet-only pair (R9).
TRACKED_BASES: Final[tuple[str, ...]] = ("linked", "wikt", "wikt+corpus")
LICENCE_BY_SA: Final[str] = "by-sa"

DERIVATION_FIELDS: Final[tuple[str, ...]] = (
    "base_lexeme_id",
    "base_root_id",
    "base_lemma",
    "base_readings",
    "base_stem",
    "base_target_class",
    "base_stem_final",
    "base_cluster",
    "base_suffix_target",
    "base_domain",
    "base_stratum_wikt",
    "base_stratum",
    "base_stratum_source",
    "base_deverbal_from",
    "derivative",
    "derivative_pos",
    "derivative_feminine_agent",
    "suffix",
    "suffix_class",
    "suffix_pair",
    "suffix_palatalizes_base",
    "ending",
    "prefix",
    "outcome",
    "orth_change",
    "mutation_size",
    "evidence",
    "basis",
    "headword_editions",
    "derinet_child_id",
    "derinet_verbalizer",
    "derinet_conj_class",
    "verbalizer_consistent",
    "conj_class_consistent",
    "derivative_stratum_wikt",
    "derivative_stratum",
    "derivative_stratum_source",
    "count",
    "lemma_count",
    "n_forms_counted",
    *(f"attested_{tau}" for tau in CORPUS_TAUS),
    "source",
    "licence_class",
)
STRATA_FIELDS: Final[tuple[str, ...]] = (
    "lexeme_id",
    "root_id",
    "lemma",
    "pos",
    "stratum_wikt",
    "stratum_wikt_basis",
    "stratum",
    "stratum_source",
    "source",
    "licence_class",
)
SUMMARY_FIELDS: Final[tuple[str, ...]] = (
    "basis",
    "tau",
    "base_target_class",
    "suffix",
    "suffix_class",
    "outcome",
    "base_stratum_wikt",
    "n_pairs",
    "n_bases",
    "source",
    "licence_class",
)


# ---------------------------------------------------------------------------
# Bases
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Base:
    """An in-domain noun lexeme that derivatives may build on."""

    lexeme_id: str
    root_id: str
    lemma: str
    readings: tuple[str, ...]
    stem: str
    target_class: str
    stem_final: str
    cluster: str
    suffix_target: bool
    etymology: str
    head_links: tuple[tuple[str, str], ...]
    domain: str = "in"

    @property
    def entry_id(self) -> str:
        """The id of the base's English Wiktionary entry
        (:func:`lexpipe.ids.lexeme_id`, before tie-breaking)."""
        return base_id(self.lexeme_id, self.lemma, "noun", self.etymology)

    @property
    def animate_masculine(self) -> bool:
        return "m.an" in self.readings


_HEAD_COLUMNS: Final[tuple[tuple[str, str], ...]] = (
    ("head_adjective", "adj"),
    ("head_diminutive", "dim"),
    ("head_feminine", "f"),
    ("head_masculine", "m"),
)
# Out-of-domain classes that still serve as derivation bases: a dual
# remnant (ruka, oko) inflects regularly in the singular (ruce), so its
# derivatives (ruční, ručit) and its DAT/LOC.SG count for R20.
BASE_OOD_CLASSES: Final[frozenset[str]] = frozenset({"dual_remnant"})


def load_bases(paths: CzechPaths) -> list[Base]:
    """Nouns of ``czech_lexicon_complete.csv`` with a stem: the
    in-domain ones and the dual remnants (``domain`` records which).

    Affix entries (-slav, -ule) are left out: a lemma that starts or
    ends with a hyphen is not a noun of the language.
    """
    out: list[Base] = []
    for row in jsonio.iter_csv(paths.lexicon_complete):
        lemma = row["lemma"]
        domain = (
            "in"
            if row["domain"] == "in"
            else (
                f"out:{row['ood_class']}"
                if row["ood_class"] in BASE_OOD_CLASSES
                else ""
            )
        )
        if (
            row["pos"] != "noun"
            or not domain
            or not row["stem"]
            or lemma.startswith("-")
            or lemma.endswith("-")
        ):
            continue
        links = tuple(
            (kind, word)
            for column, kind in _HEAD_COLUMNS
            for word in csvio.split_list(row[column])
        )
        out.append(
            Base(
                lexeme_id=row["lexeme_id"],
                root_id=row["root_id"],
                lemma=lemma,
                readings=tuple(csvio.split_list(row["readings"])),
                stem=row["stem"],
                target_class=row["target_class"],
                stem_final=row["stem_final"],
                cluster=row["cluster"],
                suffix_target=csvio.parse_bool(row["target_is_suffix"])
                is True,
                etymology=row["etymology"],
                head_links=links,
                domain=domain,
            )
        )
    return out


# ---------------------------------------------------------------------------
# Wiktionary
# ---------------------------------------------------------------------------


@dataclass
class WiktionaryIndex:
    """What the derivations stage reads from the two editions."""

    headwords: dict[str, set[tuple[str, str]]] = field(
        default_factory=lambda: defaultdict(set)
    )
    entry_strata: dict[tuple[str, str, str], set[tuple[str, str]]] = field(
        default_factory=lambda: defaultdict(set)
    )
    word_strata: dict[tuple[str, str], set[tuple[str, str]]] = field(
        default_factory=lambda: defaultdict(set)
    )
    affix_links: list[tuple[str, str, str]] = field(default_factory=list)
    # Derived/related terms keyed by the listing entry's id
    # (:attr:`Base.entry_id`).
    noun_links: dict[str, list[tuple[str, str]]] = field(
        default_factory=lambda: defaultdict(list)
    )
    derivative_bases: dict[str, set[str]] = field(
        default_factory=lambda: defaultdict(set)
    )
    # Verbs a noun entry is derived from (entry id -> verbs).
    deverbal: dict[str, tuple[str, ...]] = field(default_factory=dict)
    # Nouns glossed as female agents or "female equivalent of".
    feminine_agents: set[str] = field(default_factory=set)
    # Raw (stratum, basis) per entry and word, refined into the Czech
    # strata (adapted loan / hybrid) once every headword is known.
    raw_strata: list[tuple[str, str, str, tuple[str, str], str]] = field(
        default_factory=list
    )

    def contradicts(self, word: str, pos: str, base_lemma: str) -> bool:
        """Does *word*'s own etymology tell another story? It does when
        it derives the word from another base, or states an origin
        (inherited, borrowed) without naming this base (drožka is
        borrowed from Russian, not formed from droga)."""
        stated = {b.lower() for b in self.derivative_bases.get(word, ())}
        if base_lemma.lower() in stated:
            return False
        if stated:
            return True
        return self.stratum(word, pos)[0] != etymology.UNKNOWN

    def related_conflict(
        self, word: str, pos: str, base_lemma: str, verb: bool
    ) -> bool:
        """Does a related-terms link's derivative tell another story
        (module docstring)? Its etymology names another base
        (vystrašit: vy- + strašit), or it is borrowed whole (papoušek),
        or it is an inherited verb: a noun and a verb both inherited
        are related, the direction unknown (křik – křičet, lék –
        léčit). An inherited adjective or noun of the base stays
        (strach – strašný, vlk – vlčí)."""
        stated = {b.lower() for b in self.derivative_bases.get(word, ())}
        if base_lemma.lower() in stated:
            return False
        if stated:
            return True
        found = self.stratum(word, pos)[0]
        return found == etymology.BORROWED or (
            verb and found == etymology.INHERITED
        )

    def has_headword(self, word: str, pos: str) -> bool:
        return any(p == pos for _, p in self.headwords.get(word, ()))

    def editions(self, word: str, pos: str) -> list[str]:
        return sorted({e for e, p in self.headwords.get(word, ()) if p == pos})

    def add_stratum(
        self, word: str, pos: str, ety: str, value: tuple[str, str]
    ) -> None:
        self.entry_strata[(word, pos, ety)].add(value)
        self.word_strata[(word, pos)].add(value)

    def stratum(self, word: str, pos: str, ety: str = "") -> tuple[str, str]:
        """Wiktionary stratum of an entry; homographs that disagree (and
        no exact etymology match) give ``unknown``."""
        exact = self.entry_strata.get((word, pos, ety))
        if exact and len(exact) == 1:
            return next(iter(exact))
        found = self.word_strata.get((word, pos), set())
        strata = {s for s, _ in found}
        if len(strata) == 1:
            return sorted(found)[0]
        if not found:
            return etymology.UNKNOWN, ""
        return etymology.UNKNOWN, "conflict"


_LINK_KINDS: Final[tuple[tuple[str, str], ...]] = (
    ("derived", "wikt_derived"),
    ("related", "wikt_related"),
)


def deverbal_sources(
    found: Sequence[etymology.Step], text: str, word: str = ""
) -> tuple[str, ...]:
    """Verbs a noun's etymology derives it from (module docstring):
    the terms of ``deverbal``/back-formation steps and the words of
    "Deverbal from X", "Formed from X (… of Y)", "From the verb X"."""
    out: list[str] = []
    for step in found:
        if step.kind in _DEVERBAL_KINDS and step.lang in ("", "cs"):
            for term in step.terms:
                if term and not term.startswith("-") and term not in out:
                    out.append(term)
    prose = core_etymology.prose(text, word, etymology.PROFILE)
    for match in _DEVERBAL_TEXT.finditer(prose):
        for verb in match.groups():
            if verb and verb not in out:
                out.append(verb)
    return tuple(out)


def feminine_agent_entry(record: Mapping[str, Any]) -> bool:
    """A noun entry glossed as a female agent ("female rival", "woman
    from …") or with a "female equivalent of" sense."""
    for sense in record.get("senses") or ():
        if not isinstance(sense, Mapping):
            continue
        if extract.derivational_relation(sense) == "female_equivalent":
            return True
        glosses = [g for g in sense.get("glosses") or () if isinstance(g, str)]
        if glosses and _FEMALE_GLOSS.search(glosses[0]):
            return True
    return False


def scan_wiktionary(paths: CzechPaths) -> WiktionaryIndex:
    """One pass over each projected edition (module docstring)."""
    index = WiktionaryIndex()
    for record in jsonio.read_jsonl(paths.en_entries):
        word = str(record.get("word") or "")
        pos = str(record.get("pos") or "")
        if not word or not pos:
            continue
        index.headwords[word].add(("en", pos))
        found = etymology.steps(record.get("etymology_templates") or ())
        text = str(record.get("etymology_text") or "")
        index.raw_strata.append(
            (
                word,
                pos,
                kaikki.etymology_label(record),
                etymology.wiktionary_stratum(found, text, word=word),
                etymology.adaptation_key(found),
            )
        )
        for base, _suffix in etymology.affix_parts(found):
            index.affix_links.append((base, word, pos))
            index.derivative_bases[word].add(base)
        if pos == "noun":
            entry = ids.lexeme_id("cs", record)
            for key, kind in _LINK_KINDS:
                for item in record.get(key) or ():
                    linked = str(item.get("word") or "").strip()
                    if linked and " " not in linked:
                        index.noun_links[entry].append((kind, linked))
            sources = deverbal_sources(found, text, word)
            if sources:
                index.deverbal[entry] = sources
            if feminine_agent_entry(record):
                index.feminine_agents.add(word)
    for record in jsonio.read_jsonl(paths.cs_entries):
        word = str(record.get("word") or "")
        pos = str(record.get("pos") or "")
        if word and pos:
            index.headwords[word].add(("cs", pos))
    lemmas = frozenset(index.headwords)
    for word, pos, label, raw, key in index.raw_strata:
        index.add_stratum(
            word, pos, label, etymology.refine_stratum(raw, key, lemmas)
        )
    index.raw_strata.clear()
    return index


def wiktionary_forms(
    paths: CzechPaths, wanted: set[tuple[str, str]]
) -> dict[tuple[str, str], list[str]]:
    """Paradigm forms (table forms, lemma included) of the wanted
    ``(word, pos)`` entries, English edition first, then Czech."""
    out: dict[tuple[str, str], list[str]] = defaultdict(list)
    for edition_path, needs_source in (
        (paths.en_entries, True),
        (paths.cs_entries, False),
    ):
        for record in jsonio.read_jsonl(edition_path):
            key = (str(record.get("word") or ""), str(record.get("pos") or ""))
            if key not in wanted:
                continue
            forms = out[key]
            if key[0] not in forms:
                forms.append(key[0])
            for form in kaikki.inflected_forms(record):
                if needs_source and not form.source:
                    continue
                if " " in form.form or form.form in forms:
                    continue
                forms.append(form.form)
    return out


# ---------------------------------------------------------------------------
# Pairs
# ---------------------------------------------------------------------------


@dataclass
class Pair:
    """One base × derivative pair while the stage builds it."""

    base: Base
    derivative: str
    decomposition: sx.Decomposition
    evidence: set[str] = field(default_factory=set)
    edge: derinet.DerinetEdge | None = None

    @property
    def suffix(self) -> sx.Suffix:
        return self.decomposition.suffix

    @property
    def key(self) -> tuple[str, str, str]:
        return (self.base.lexeme_id, self.suffix.label, self.derivative)

    def basis(self) -> str:
        """The narrowest basis that contains the pair (module docstring)."""
        if self.evidence & LINK_EVIDENCE:
            return "linked"
        if self.evidence & HEADWORD_CANDIDATES:
            return "wikt"
        if "cand_corpus" in self.evidence:
            return "wikt+corpus"
        return "all"


class PairSet:
    """Pairs keyed by (base lexeme, suffix, derivative)."""

    def __init__(self) -> None:
        self.pairs: dict[tuple[str, str, str], Pair] = {}
        self.unmatched: Counter[str] = Counter()

    def add(
        self,
        base: Base,
        derivative: str,
        evidence: str,
        *,
        allow_prefix: bool,
        suffixes: Sequence[sx.Suffix] | None = None,
        edge: derinet.DerinetEdge | None = None,
    ) -> bool:
        if derivative.lower() == base.lemma.lower():
            return False
        found = sx.decompose(
            derivative,
            base.stem,
            base.lemma,
            cluster=base.cluster,
            allow_prefix=allow_prefix,
            suffixes=suffixes,
        )
        if found is None:
            self.unmatched[evidence] += 1
            return False
        pair = Pair(base, derivative, found)
        current = self.pairs.setdefault(pair.key, pair)
        current.evidence.add(evidence)
        if edge is not None and current.edge is None:
            current.edge = edge
        return True


def _bases_by_lemma(bases: Sequence[Base]) -> dict[str, list[Base]]:
    out: dict[str, list[Base]] = defaultdict(list)
    for base in bases:
        out[base.lemma].append(base)
    return out


def add_wiktionary_pairs(
    pairs: PairSet, bases: Sequence[Base], index: WiktionaryIndex
) -> None:
    """``wikt_ety``, ``wikt_head``, ``wikt_derived``, ``wikt_related``
    (derived and related terms from the base's own entry)."""
    by_lemma = _bases_by_lemma(bases)
    for base_word, derivative, _pos in sorted(index.affix_links):
        for base in by_lemma.get(base_word, ()):
            pairs.add(base, derivative, "wikt_ety", allow_prefix=True)
    for base in bases:
        for _kind, word in base.head_links:
            pairs.add(base, word, "wikt_head", allow_prefix=True)
        for kind, word in index.noun_links.get(base.entry_id, ()):
            pairs.add(base, word, kind, allow_prefix=True)


# ---------------------------------------------------------------------------
# Pruning (module docstring: homographs, related terms, direction, names)
# ---------------------------------------------------------------------------


def prune_related(pairs: PairSet, index: WiktionaryIndex) -> int:
    """Drop ``wikt_related`` evidence whose derivative tells another
    story (:meth:`WiktionaryIndex.related_conflict`); a pair left
    without evidence is dropped. Returns the pairs dropped."""
    dropped = 0
    for key in sorted(pairs.pairs):
        pair = pairs.pairs[key]
        if "wikt_related" not in pair.evidence:
            continue
        if not index.related_conflict(
            pair.derivative,
            pair.suffix.pos,
            pair.base.lemma,
            pair.suffix.label in VERB_SUFFIXES,
        ):
            continue
        pair.evidence.discard("wikt_related")
        if not pair.evidence:
            del pairs.pairs[key]
            dropped += 1
    return dropped


def _verb_core(verb: str) -> set[str]:
    """A verb's root keys for the direction check: infinitive suffix
    and l-participle ending off, folded, plus the stem with a final
    first-palatalization output undone (ležet ~ lehat: leh)."""
    lowered = verb.lower()
    for ending in (
        "ovat",
        "nout",
        "ila",
        "ilo",
        "il",
        "ěl",
        "el",
        "al",
        "at",
        "át",
        "it",
        "ít",
        "et",
        "ět",
        "out",
    ):
        if lowered.endswith(ending) and len(lowered) > len(ending) + 1:
            lowered = lowered[: -len(ending)]
            break
    return morph.undone(lowered)


def _verb_cores(verb: str) -> set[str]:
    """:func:`_verb_core` of the verb and of the verb without a prefix
    (obložit, ložit)."""
    lowered = verb.lower()
    out = set(_verb_core(lowered))
    for prefix in sx.PREFIXES:
        if lowered.startswith(prefix) and len(lowered) > len(prefix) + 3:
            out |= _verb_core(lowered[len(prefix) :])
    return out


def reversed_pair(pair: Pair, sources: Sequence[str]) -> bool:
    """Is a verb derivative the source of a deverbal base (module
    docstring)? The derivative is one of the source verbs, or an -it/-et
    verb (or the -ovat imperfective of one, with the changed stem:
    vyloučit, vylučovat ← výluka) of the same root (prefix and
    infinitive suffix off, palatalization undone). An -ovat verb on the
    plain stem is formed from the noun (poloha → polohovat, tlak →
    tlakovat), and a derivative whose own etymology names the base
    (``wikt_ety``) is never reversed."""
    if pair.suffix.label not in VERB_SUFFIXES or not sources:
        return False
    if "wikt_ety" in pair.evidence:
        return False
    word = pair.derivative.lower()
    if any(verb.lower() == word for verb in sources):
        return True
    if pair.suffix.label == "-ovat" and not pair.decomposition.mutated:
        return False
    derivative = _verb_cores(pair.derivative)
    return any(derivative & _verb_cores(verb) for verb in sources)


def drop_reversed(
    pairs: PairSet, deverbal: Mapping[str, tuple[str, ...]]
) -> int:
    """Drop the reversed pairs (:func:`reversed_pair`); return how
    many."""
    dropped = 0
    for key in sorted(pairs.pairs):
        pair = pairs.pairs[key]
        if reversed_pair(pair, deverbal.get(pair.base.entry_id, ())):
            del pairs.pairs[key]
            dropped += 1
    return dropped


def drop_names(pairs: PairSet, names_only: frozenset[str]) -> int:
    """Drop pairs whose derivative is only a proper name (R25)."""
    dropped = 0
    for key in sorted(pairs.pairs):
        if pairs.pairs[key].derivative in names_only:
            del pairs.pairs[key]
            dropped += 1
    return dropped


def feminine_agent(pair: Pair, index: WiktionaryIndex) -> bool:
    """A feminizing derivative: -yně, or -ka named on the base's
    headword line as its feminine (f=) or glossed as a female agent
    (module docstring)."""
    if pair.suffix.label == "-yně":
        return True
    if pair.suffix.label != "-ka":
        return False
    if ("f", pair.derivative) in pair.base.head_links:
        return True
    return pair.derivative in index.feminine_agents


def _dedupe_rank(pair: Pair, index: WiktionaryIndex) -> tuple[Any, ...]:
    agent = feminine_agent(pair, index)
    return (
        0 if not agent or pair.base.animate_masculine else 1,
        0 if pair.evidence & OWN_EVIDENCE else 1,
        -len(pair.evidence),
        pair.base.lexeme_id,
    )


def dedupe_homographs(pairs: PairSet, index: WiktionaryIndex) -> int:
    """Keep a derivative on one root among homograph roots of one
    lemma (module docstring: sokyně on sok 'rival', not on the
    inanimate sok; dračí on drak 'dragon'); return the pairs
    dropped."""
    groups: dict[tuple[str, str, str], list[tuple[str, str, str]]] = (
        defaultdict(list)
    )
    for key, pair in pairs.pairs.items():
        # case-insensitive: Buddha and buddha both give buddhista
        lemma = pair.base.lemma.lower()
        groups[(lemma, pair.suffix.label, pair.derivative.lower())].append(key)
    dropped = 0
    for group in groups.values():
        roots = {pairs.pairs[key].base.root_id for key in group}
        if len(roots) < 2:
            continue
        ranked = sorted(
            group, key=lambda k: _dedupe_rank(pairs.pairs[k], index)
        )
        keep_root = pairs.pairs[ranked[0]].base.root_id
        for key in group:
            if pairs.pairs[key].base.root_id != keep_root:
                del pairs.pairs[key]
                dropped += 1
    return dropped


def add_derinet_pairs(
    pairs: PairSet,
    bases: Sequence[Base],
    nouns: Sequence[derinet.DerinetNoun],
    edges: Sequence[derinet.DerinetEdge],
) -> dict[str, list[derinet.DerinetNoun]]:
    """``derinet`` evidence; returns the DeriNet nouns joined to each
    base lexeme (R14: lemma + compatible gender/animacy)."""
    by_lemma = _bases_by_lemma(bases)
    node_bases: dict[str, list[Base]] = defaultdict(list)
    joined: dict[str, list[derinet.DerinetNoun]] = defaultdict(list)
    for noun in nouns:
        for base in by_lemma.get(noun.lemma, ()):
            if derinet.compatible(base.readings, noun.gender, noun.animacy):
                node_bases[noun.node_id].append(base)
                joined[base.lexeme_id].append(noun)
    for edge in edges:
        for base in node_bases.get(edge.parent_id, ()):
            pairs.add(
                base,
                edge.child_lemma,
                "derinet",
                allow_prefix=True,
                edge=edge,
            )
    return joined


# ---------------------------------------------------------------------------
# Candidates (R21) and corpus counts
# ---------------------------------------------------------------------------


class Counts:
    """Corpus counts of word forms: case-sensitive tokens plus the
    sentence-initial capitalized spelling of a lowercase form.

    Generated forms are concatenations of NFC strings, so the
    case-sensitive table is read directly (no re-normalization).
    """

    def __init__(self, table: AttestationTable) -> None:
        self.table = table
        self._cs = table.counts_cs

    def raw(self, form: str) -> int:
        count = self._cs.get(form, 0)
        if form[:1].islower():
            count += self._cs.get(form[:1].upper() + form[1:], 0)
        return count

    def any_tokens(self, forms: Iterable[str]) -> bool:
        return any(self.raw(form) for form in forms)


def _spellings(word: str, base: Base) -> list[str]:
    """A candidate as generated (lowercase) and, for a capitalized base
    (Čech → Češka), capitalized."""
    out = [word]
    if base.lemma[:1].isupper():
        out.append(word[:1].upper() + word[1:])
    return out


@dataclass(frozen=True, slots=True)
class Found:
    """A generated candidate that was found (R21)."""

    base: Base
    word: str
    suffix: sx.Suffix
    evidence: str
    expected: bool
    mutated: bool = False

    @property
    def plain_palatalizing(self) -> bool:
        """The plain spelling of a P suffix that palatalizes the base's
        target (pantofágní): a possible exception (module docstring)."""
        return not self.mutated and self.suffix.palatalizes(
            self.base.target_class, self.base.cluster
        )


def expected_spelling(base: Base, cand: sx.Candidate) -> bool:
    """Palatalized with the first (big) palatalization for a suffix that
    palatalizes the base's target (derivation never takes the second:
    drozní is drozdí, not from droga; -ský merges are size ""), plain
    otherwise."""
    if not cand.suffix.palatalizes(base.target_class, base.cluster):
        return not cand.mutated
    return cand.mutated and cand.size in ("big", "")


def _generated(base: Base) -> Iterator[sx.Candidate]:
    """Candidates of one base: both spellings for P and NT suffixes,
    the plain spelling only for N suffixes (a mutated -ový, -ovat,
    -yně or -ův spelling is not a Czech derivative; its hits are other
    words: pracovat ← práce, not prak)."""
    for suffix in sx.SUFFIXES:
        for cand in sx.candidates(base.stem, base.cluster, suffix):
            if cand.mutated and suffix.klass == sx.N_CLASS:
                continue
            if cand.word != base.lemma.lower():
                yield cand


def candidate_words(
    bases: Sequence[Base], index: WiktionaryIndex, counts: Counts
) -> list[Found]:
    """Generated derivatives found as headwords of the suffix's part of
    speech, or in the corpus without being any headword."""
    found: list[Found] = []
    for base in bases:
        if base.target_class not in CANDIDATE_CLASSES:
            continue
        if len(sx.fold(base.stem)) < MIN_STEM:
            continue
        for cand in _generated(base):
            suffix = cand.suffix
            expected = expected_spelling(base, cand)
            spellings = _spellings(cand.word, base)
            headword = next(
                (w for w in spellings if index.has_headword(w, suffix.pos)),
                None,
            )
            if headword is not None:
                found.append(
                    Found(
                        base,
                        headword,
                        suffix,
                        "cand_headword",
                        expected,
                        cand.mutated,
                    )
                )
            elif any(w in index.headwords for w in spellings):
                continue  # a headword of another part of speech
            elif counts.any_tokens(sx.paradigm(cand.word, suffix)):
                found.append(
                    Found(
                        base,
                        cand.word,
                        suffix,
                        "cand_corpus",
                        expected,
                        cand.mutated,
                    )
                )
    return found


def shared_stems(found: Sequence[Found]) -> set[tuple[str, str]]:
    """``(word, suffix)`` candidates generated from bases with different
    lemmas (deka and deko both give dečka): which base the word belongs
    to cannot be told from its shape."""
    lemmas: dict[tuple[str, str], set[str]] = defaultdict(set)
    for item in found:
        lemmas[(item.word.lower(), item.suffix.label)].add(
            item.base.lemma.lower()
        )
    return {key for key, names in lemmas.items() if len(names) > 1}


def candidate_problem(
    item: Found, index: WiktionaryIndex, shared: set[tuple[str, str]]
) -> str:
    """Why a candidate fails the R21 checks ("" if it passes). A
    headword candidate in the plain spelling of a palatalizing suffix
    passes the spelling check (module docstring)."""
    if not item.expected and not (
        item.evidence == "cand_headword" and item.plain_palatalizing
    ):
        return "unexpected spelling"
    if (item.word.lower(), item.suffix.label) in shared:
        return "stem shared by several bases"
    if index.contradicts(item.word, item.suffix.pos, item.base.lemma):
        return "own etymology, other base"
    if item.suffix.klass == sx.NT_CLASS and _ie_noun(item.base, index):
        return "-ista/-ismus of an -ie noun"
    if (
        item.evidence == "cand_corpus"
        and len(sx.fold(item.base.stem)) < MIN_CORPUS_STEM
    ):
        return "short stem"
    return ""


def _ie_noun(base: Base, index: WiktionaryIndex) -> bool:
    """Is there a noun "<base stem>ie" (ekolog: ekologie)? Its -ista,
    -ismus derivatives are the -ie noun's (module docstring)."""
    stem = base.stem
    return any(
        index.has_headword(word, "noun")
        for word in {stem + "ie", stem.lower() + "ie"}
    )


def add_candidates(
    pairs: PairSet,
    found: Sequence[Found],
    index: WiktionaryIndex,
    rejected: Counter[str],
) -> None:
    """Add the candidates that pass :func:`candidate_problem` (to a new
    pair or as evidence of an existing one). A candidate that fails
    only for its spelling still corroborates a pair that a Wiktionary
    link already found (barokní); nothing else is added, so a
    DeriNet-only pair never gains Wiktionary-basis evidence from a
    candidate that fails the checks. Corpus-only candidates are also
    checked for tokens of their own once the form index is read."""
    shared = shared_stems(found)
    for item in found:
        key = (item.base.lexeme_id, item.suffix.label, item.word)
        problem = candidate_problem(item, index, shared)
        existing = pairs.pairs.get(key)
        if problem == "unexpected spelling" and existing is not None:
            if existing.evidence & LINK_EVIDENCE:
                existing.evidence.add(item.evidence)
                continue
        if problem:
            rejected[f"{item.evidence}: {problem}"] += 1
            continue
        evidence = item.evidence
        if not item.expected:
            evidence = "cand_headword_plain"  # accepted, flagged
        if existing is not None:
            existing.evidence.add(evidence)
            continue
        pairs.add(
            item.base,
            item.word,
            evidence,
            allow_prefix=False,
            suffixes=(item.suffix,),
        )


def form_owners(paths: CzechPaths, forms: set[str]) -> dict[str, set[str]]:
    """``form -> owner words`` from the form index (lemma, paradigm and
    form-of rows; :func:`.attestation.counts_as_owner`), for the given
    exact forms."""
    owners: dict[str, set[str]] = defaultdict(set)
    if not forms:
        return owners
    uri = paths.form_index.resolve().as_uri() + "?mode=ro"
    with contextlib.closing(sqlite3.connect(uri, uri=True)) as con:
        con.execute("CREATE TEMP TABLE wanted (form TEXT PRIMARY KEY)")
        con.executemany(
            "INSERT INTO wanted VALUES (?)", ((f,) for f in sorted(forms))
        )
        rows = con.execute(
            "SELECT f.form, f.kind, f.owner_word, f.target FROM forms f "
            "JOIN wanted w ON f.form = w.form"
        )
        for form, kind, owner_word, target in rows:
            if not counts_as_owner(kind, target):
                continue  # spelling variants, forms of phrases
            owners[form].add(
                owner_word if kind in ("lemma", "paradigm") else target
            )
    return owners


@dataclass(frozen=True, slots=True)
class Attestation:
    """Corpus attestation of one derivative."""

    count: int
    lemma_count: int
    n_forms: int
    flags: Mapping[str, bool]


def _capitalized(form: str) -> str:
    return form[:1].upper() + form[1:]


def attest_derivative(
    word: str,
    forms: Sequence[str],
    owners: Mapping[str, set[str]],
    counts: Counts,
) -> Attestation:
    """Summed tokens of the forms that belong to *word* only.

    The paradigm is attested with the core's
    :func:`lexpipe.attest.attest_paradigm` (case-sensitive): every form
    and, for a lowercase form, its capitalized spelling (sentence-initial
    uses), each distinct spelling counted once. A form that is also a
    form of another word is left out with its capitalized spelling
    (unless that spelling is one of the forms itself), and a
    capitalized spelling that belongs to another word (Poláci) is left
    out alone. ``n_forms`` counts the forms (not their capitalized
    spellings) that enter the sum.
    """
    own = {word, word.lower()}
    distinct = list(dict.fromkeys(forms))
    listed = set(distinct)
    excluded: set[str] = set()
    spellings: list[str] = []
    for form in distinct:
        owned = bool(owners.get(form, set()) - own)
        if owned:
            excluded.add(form)
        spellings.append(form)
        if form[:1].islower():
            cap = _capitalized(form)
            spellings.append(cap)
            if cap not in listed and (owned or owners.get(cap, set()) - own):
                excluded.add(cap)
    result = attest.attest_paradigm(
        counts.table,
        spellings,
        taus=CORPUS_TAUS,
        case_sensitive=True,
        exclude=excluded,
    )
    return Attestation(
        count=result.total,
        lemma_count=counts.raw(word),
        n_forms=sum(1 for form in distinct if form not in excluded),
        flags=dict(result.attested_by_tau),
    )


# ---------------------------------------------------------------------------
# Strata
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class StrataSources:
    """Everything the stratum of a base or derivative is read from."""

    index: WiktionaryIndex
    czetyl: Mapping[str, str]
    derinet_nouns: Mapping[str, list[derinet.DerinetNoun]]

    def base(self, base: Base) -> tuple[str, etymology.Stratum]:
        wikt = self.index.stratum(base.lemma, "noun", base.etymology)
        loan_root = any(
            noun.loan_root
            for noun in self.derinet_nouns.get(base.lexeme_id, ())
        )
        combined = etymology.combine(
            wikt, self.czetyl.get(base.lemma), loan_root or None
        )
        return wikt[0], combined

    def derivative(self, word: str, pos: str) -> tuple[str, etymology.Stratum]:
        wikt = self.index.stratum(word, pos)
        combined = etymology.combine(wikt, self.czetyl.get(word), None)
        return wikt[0], combined


def strata_rows(
    paths: CzechPaths, strata: StrataSources, bases: Sequence[Base]
) -> list[dict[str, Any]]:
    """R22 strata of every noun and adjective lexeme of the lexicon
    (``czech_lexeme_strata.csv``, NC: czetyl and DeriNet columns)."""
    by_id = {base.lexeme_id: base for base in bases}
    rows: list[dict[str, Any]] = []
    for row in jsonio.iter_csv(paths.lexicon_complete):
        pos = row["pos"]
        if pos not in ("noun", "adj"):
            continue
        base = by_id.get(row["lexeme_id"])
        if base is not None:
            wikt_label, combined = strata.base(base)
            basis = strata.index.stratum(base.lemma, "noun", base.etymology)[1]
        else:
            wikt = strata.index.stratum(row["lemma"], pos, row["etymology"])
            wikt_label, basis = wikt
            combined = etymology.combine(
                wikt, strata.czetyl.get(row["lemma"]), None
            )
        sources = [SOURCE_EN]
        if combined.source == etymology.SOURCE_CZETYL:
            sources.append("czetyl")
        elif combined.source == etymology.SOURCE_DERINET:
            sources.append("derinet")
        rows.append(
            {
                "lexeme_id": row["lexeme_id"],
                "root_id": row["root_id"],
                "lemma": row["lemma"],
                "pos": pos,
                "stratum_wikt": wikt_label,
                "stratum_wikt_basis": basis,
                "stratum": combined.stratum,
                "stratum_source": combined.source,
                "source": sources,
                "licence_class": (
                    "by-sa" if sources == [SOURCE_EN] else "by-nc-sa"
                ),
            }
        )
    return rows


# ---------------------------------------------------------------------------
# Rows
# ---------------------------------------------------------------------------


def pair_row(
    pair: Pair,
    strata: StrataSources,
    attestation: Attestation,
) -> dict[str, Any]:
    """One row of ``czech_derivations.csv``."""
    base = pair.base
    index = strata.index
    found = pair.decomposition
    suffix = found.suffix
    base_wikt, base_all = strata.base(base)
    deriv_wikt, deriv_all = strata.derivative(pair.derivative, suffix.pos)
    edge = pair.edge
    verbalizer_ok: bool | None = None
    conj_ok: bool | None = None
    if edge is not None and edge.verbalizer:
        verbalizer_ok = edge.verbalizer == suffix.label
        conj_ok = derinet.conj_consistent(edge.verbalizer, edge.conj_class)
    sources = [SOURCE_EN, SOURCE_CS]
    if "derinet" in pair.evidence:
        sources.append("derinet")
    if etymology.SOURCE_CZETYL in (base_all.source, deriv_all.source):
        sources.append("czetyl")
    sources.append("leipzig")
    return {
        "base_lexeme_id": base.lexeme_id,
        "base_root_id": base.root_id,
        "base_lemma": base.lemma,
        "base_readings": list(base.readings),
        "base_stem": base.stem,
        "base_target_class": base.target_class,
        "base_stem_final": base.stem_final,
        "base_cluster": base.cluster,
        "base_suffix_target": base.suffix_target,
        "base_domain": base.domain,
        "base_stratum_wikt": base_wikt,
        "base_stratum": base_all.stratum,
        "base_stratum_source": base_all.source,
        "base_deverbal_from": list(index.deverbal.get(base.entry_id, ())),
        "derivative": pair.derivative,
        "derivative_pos": suffix.pos,
        "derivative_feminine_agent": feminine_agent(pair, index),
        "suffix": suffix.label,
        "suffix_class": suffix.klass,
        "suffix_pair": suffix.pair,
        "suffix_palatalizes_base": suffix.palatalizes(
            base.target_class, base.cluster
        ),
        "ending": found.ending,
        "prefix": found.prefix,
        "outcome": found.outcome(base.target_class, base.cluster),
        "orth_change": found.orth_change,
        "mutation_size": found.size,
        "evidence": sorted(pair.evidence),
        "basis": pair.basis(),
        "headword_editions": strata.index.editions(
            pair.derivative, suffix.pos
        ),
        "derinet_child_id": edge.child_id if edge else "",
        "derinet_verbalizer": edge.verbalizer if edge else "",
        "derinet_conj_class": edge.conj_class if edge else "",
        "verbalizer_consistent": verbalizer_ok,
        "conj_class_consistent": conj_ok,
        "derivative_stratum_wikt": deriv_wikt,
        "derivative_stratum": deriv_all.stratum,
        "derivative_stratum_source": deriv_all.source,
        "count": attestation.count,
        "lemma_count": attestation.lemma_count,
        "n_forms_counted": attestation.n_forms,
        **{f"attested_{tau}": attestation.flags[tau] for tau in CORPUS_TAUS},
        "source": sources,
        "licence_class": "by-nc-sa",
    }


def in_basis(row: Mapping[str, Any], basis: str) -> bool:
    """Is a pair row inside *basis* (``linked`` ⊂ ``wikt`` ⊂
    ``wikt+corpus`` ⊂ ``all``)? Works on written rows and built rows."""
    return BASES.index(str(row["basis"])) <= BASES.index(basis)


def row_attested(row: Mapping[str, Any], tau: str) -> bool:
    """τ flag of a pair row (dictionary: always)."""
    if tau == "dictionary":
        return True
    value = row[f"attested_{tau}"]
    if isinstance(value, bool):
        return value
    return csvio.parse_bool(str(value)) is True


def summary_rows(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Tracked aggregates (:data:`TRACKED_BASES` only)."""
    pairs: Counter[tuple[str, ...]] = Counter()
    bases: dict[tuple[str, ...], set[str]] = defaultdict(set)
    for row in rows:
        for basis in TRACKED_BASES:
            if not in_basis(row, basis):
                continue
            for tau in TAUS:
                if not row_attested(row, tau):
                    continue
                key = (
                    basis,
                    tau,
                    str(row["base_target_class"]),
                    str(row["suffix"]),
                    str(row["suffix_class"]),
                    str(row["outcome"]),
                    str(row["base_stratum_wikt"]),
                )
                pairs[key] += 1
                bases[key].add(str(row["base_root_id"]))
    out: list[dict[str, Any]] = []
    for label in sorted(pairs):
        out.append(
            {
                "basis": label[0],
                "tau": label[1],
                "base_target_class": label[2],
                "suffix": label[3],
                "suffix_class": label[4],
                "outcome": label[5],
                "base_stratum_wikt": label[6],
                "n_pairs": pairs[label],
                "n_bases": len(bases[label]),
                "source": TRACKED_SOURCE,
                "licence_class": LICENCE_BY_SA,
            }
        )
    return out


# ---------------------------------------------------------------------------
# Stage
# ---------------------------------------------------------------------------


def run(ctx: RunContext) -> int | None:
    """Stage entry point (see the module docstring)."""
    paths = paths_for(ctx)
    apaths = analysis_paths(ctx)
    words = leipzig.ensure_words_file(ctx, paths)
    return run_stage(
        ctx,
        STAGE,
        inputs=[
            paths.lexicon_complete,
            paths.lexicon_full,
            paths.en_entries,
            paths.cs_entries,
            paths.form_index,
            words,
            apaths.derinet_source,
            apaths.czetyl_source,
        ],
        params={"version": VERSION, "min_stem": MIN_STEM},
        outputs=[
            apaths.derinet_db,
            apaths.derivations,
            apaths.lexeme_strata,
            apaths.derivations_summary,
            apaths.derivations_report,
        ],
        body=lambda: _derivations(ctx, paths, apaths),
    )


@dataclass
class _Stats:
    bases: int = 0
    projection: derinet.ProjectionCounts = field(
        default_factory=derinet.ProjectionCounts
    )
    unmatched: Counter[str] = field(default_factory=Counter)
    candidates: Counter[str] = field(default_factory=Counter)
    rejected: Counter[str] = field(default_factory=Counter)
    pruned: Counter[str] = field(default_factory=Counter)


def names_only(paths: CzechPaths) -> frozenset[str]:
    """Lemmas whose every lexeme is a proper name (R25)."""
    kinds: dict[str, set[bool]] = defaultdict(set)
    for row in jsonio.iter_csv(paths.lexicon_full):
        kinds[row["lemma"]].add(row["ood_class"].startswith("name:"))
    return frozenset(
        lemma for lemma, found in kinds.items() if found == {True}
    )


def _derivations(
    ctx: RunContext, paths: CzechPaths, apaths: AnalysisPaths
) -> int:
    stats = _Stats()
    bases = load_bases(paths)
    stats.bases = len(bases)
    print(f"  bases: {len(bases):,} nouns", flush=True)
    noun_lemmas = frozenset(base.lemma for base in bases)
    stats.projection = derinet.project(
        ctx, apaths.derinet_source, apaths.derinet_db, noun_lemmas
    )
    nouns, edges = derinet.read_projection(apaths.derinet_db)
    print(
        f"  DeriNet: {stats.projection.nodes:,} nodes, "
        f"{len(nouns):,} lexicon nouns, {len(edges):,} noun-parent edges",
        flush=True,
    )
    index = scan_wiktionary(paths)
    words_file = leipzig.ensure_words_file(ctx, paths)
    counts = Counts(leipzig.attestation_table(words_file))
    pairs = PairSet()
    add_wiktionary_pairs(pairs, bases, index)
    joined = add_derinet_pairs(pairs, bases, nouns, edges)
    found = candidate_words(bases, index, counts)
    stats.candidates = Counter(item.evidence for item in found)
    add_candidates(pairs, found, index, stats.rejected)
    stats.unmatched = pairs.unmatched
    stats.pruned["related terms: own etymology, other base"] = prune_related(
        pairs, index
    )
    stats.pruned["reversed: base derived from the verb"] = drop_reversed(
        pairs, index.deverbal
    )
    stats.pruned["derivative is a proper name (R25)"] = drop_names(
        pairs, names_only(paths)
    )
    stats.pruned["homograph roots: kept on one root"] = dedupe_homographs(
        pairs, index
    )
    strata = StrataSources(
        index, etymology.load_czetyl(apaths.czetyl_source), joined
    )
    rows = _rows(paths, pairs, index, strata, counts, stats)
    csvio.write_rows(apaths.derivations, rows, DERIVATION_FIELDS, ctx=ctx)
    csvio.write_rows(
        apaths.lexeme_strata,
        strata_rows(paths, strata, bases),
        STRATA_FIELDS,
        ctx=ctx,
    )
    csvio.write_tracked(
        apaths.derivations_summary,
        summary_rows(rows),
        SUMMARY_FIELDS,
        ctx=ctx,
    )
    report = render_report(rows, stats)
    ctx.write_tracked_text(
        apaths.derivations_report,
        report,
        sources=REPORT_SOURCES,
        nc_aggregates=True,
    )
    print(report, flush=True)
    return 0


def _rows(
    paths: CzechPaths,
    pairs: PairSet,
    index: WiktionaryIndex,
    strata: StrataSources,
    counts: Counts,
    stats: _Stats,
) -> list[dict[str, Any]]:
    ordered = [pairs.pairs[key] for key in sorted(pairs.pairs)]
    wanted = {
        (pair.derivative, pair.suffix.pos)
        for pair in ordered
        if index.has_headword(pair.derivative, pair.suffix.pos)
    }
    wikt_forms = wiktionary_forms(paths, wanted)
    forms_of: dict[tuple[str, str], list[str]] = {}
    all_forms: set[str] = set()
    for pair in ordered:
        key = (pair.derivative, pair.suffix.pos)
        if key in forms_of:
            continue
        forms = wikt_forms.get(key) or list(
            sx.paradigm(pair.derivative, pair.suffix)
        )
        forms_of[key] = forms
        for form in forms:
            all_forms.add(form)
            all_forms.add(form[:1].upper() + form[1:])
    owners = form_owners(paths, all_forms)
    cache: dict[tuple[str, str], Attestation] = {}
    rows: list[dict[str, Any]] = []
    for pair in ordered:
        key = (pair.derivative, pair.suffix.pos)
        if key not in cache:
            cache[key] = attest_derivative(
                pair.derivative, forms_of[key], owners, counts
            )
        if pair.evidence == {"cand_corpus"} and (
            cache[key].count < MIN_CORPUS_TOKENS or cache[key].lemma_count == 0
        ):
            # Its tokens belong to other entries' forms, or the lemma
            # form itself never occurs (R21).
            stats.rejected["cand_corpus: too few tokens of its own"] += 1
            continue
        rows.append(pair_row(pair, strata, cache[key]))
    return rows


def _pct(part: int, whole: int) -> str:
    return f"{100 * part / whole:5.1f}%" if whole else "   n/a"


def render_report(rows: Sequence[Mapping[str, Any]], stats: _Stats) -> str:
    """Tracked aggregate report of the stage (counts only, R9)."""
    lines = [
        "Czech derivations: base x derivative pairs (design R14, R20-R22)",
        "=" * 72,
        "",
        "Aggregates only (R9). DeriNet 2.3 and czetyl are CC BY-NC-SA: they",
        "appear here as counts; the pairs themselves are in",
        "data/czech/nc/czech_derivations.csv (gitignored).",
        "",
        f"noun bases (in-domain nouns and dual remnants): {stats.bases:,}",
        f"DeriNet: {stats.projection.nodes:,} nodes in "
        f"{stats.projection.trees:,} trees; {stats.projection.fictitious:,} "
        "Fictitious nodes dropped;",
        f"  {stats.projection.nouns_kept:,} nodes of lexicon noun lemmas, "
        f"{stats.projection.edges_kept:,} Type=Derivation edges with such a "
        "noun parent",
        "",
        "Links that do not decompose into base + a listed suffix (dropped):",
    ]
    for kind, count in sorted(stats.unmatched.items()):
        lines.append(f"  {kind:<14} {count:>8,}")
    lines += ["", "Generated candidates found (R21), before the checks:"]
    for kind, count in sorted(stats.candidates.items()):
        lines.append(f"  {kind:<14} {count:>8,}")
    lines += ["", "Candidates rejected (not corroborated by a link):"]
    for reason, count in sorted(stats.rejected.items()):
        lines.append(f"  {reason:<44} {count:>8,}")
    lines += ["", "Pairs dropped after linking (module docstring):"]
    for reason, count in sorted(stats.pruned.items()):
        lines.append(f"  {reason:<44} {count:>8,}")
    evidence: Counter[str] = Counter()
    basis: Counter[str] = Counter()
    for row in rows:
        basis[str(row["basis"])] += 1
        for item in row["evidence"]:
            evidence[str(item)] += 1
    lines += ["", f"pairs: {len(rows):,}; by narrowest basis:"]
    for name in BASES:
        lines.append(f"  {name:<12} {basis[name]:>8,}")
    lines += ["", "pairs by evidence (a pair can have several):"]
    for kind, count in sorted(evidence.items()):
        lines.append(f"  {kind:<14} {count:>8,}")
    lines += _verbalizer_lines(rows)
    lines += _outcome_lines(rows)
    lines += _strata_lines(rows)
    return "\n".join(lines) + "\n"


def _verbalizer_lines(rows: Sequence[Mapping[str, Any]]) -> list[str]:
    agree: Counter[tuple[str, str]] = Counter()
    for row in rows:
        verbalizer = str(row["derinet_verbalizer"])
        if not verbalizer or row["derivative_pos"] != "verb":
            continue
        agree[
            (
                f"{row['suffix']} vs DeriNet {verbalizer}",
                f"ConjugClass {row['conj_class_consistent']}",
            )
        ] += 1
    lines = [
        "",
        "Verbs with a DeriNet edge: suffix read from the lemma vs the",
        "verbalizer read from DeriNet's last D morph (R14), and whether",
        "ConjugClass agrees with that verbalizer:",
    ]
    for (pair, conj), count in sorted(agree.items()):
        lines.append(f"  {pair:<28} {conj:<20} {count:>6,}")
    return lines


def _outcome_lines(rows: Sequence[Mapping[str, Any]]) -> list[str]:
    lines = [
        "",
        "Velar bases, basis wikt, tau = dictionary: pairs by suffix and",
        "outcome (pal / plain / exception / nontarget)",
    ]
    table: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        if row["base_target_class"] != "velar" or not in_basis(row, "wikt"):
            continue
        table[str(row["suffix"])][str(row["outcome"])] += 1
    for suffix in [s.label for s in sx.SUFFIXES]:
        counter = table.get(suffix)
        if not counter:
            continue
        cells = "  ".join(f"{o} {counter[o]:>5,}" for o in sx.OUTCOMES)
        lines.append(f"  {suffix:<12} {cells}")
    return lines


def _strata_lines(rows: Sequence[Mapping[str, Any]]) -> list[str]:
    lines = ["", "Base strata (R22) of the pairs' bases, one count per base:"]
    seen: dict[str, tuple[str, str, str]] = {}
    for row in rows:
        seen.setdefault(
            str(row["base_lexeme_id"]),
            (
                str(row["base_stratum_wikt"]),
                str(row["base_stratum"]),
                str(row["base_stratum_source"]) or "-",
            ),
        )
    wikt = Counter(value[0] for value in seen.values())
    combined = Counter((value[1], value[2]) for value in seen.values())
    for name in etymology.STRATA:
        lines.append(f"  Wiktionary only   {name:<20} {wikt[name]:>7,}")
    for (name, source), count in sorted(combined.items()):
        lines.append(
            f"  all sources (NC)  {name:<20} {source:<10} {count:>7,}"
        )
    return lines


def iter_pair_rows(path: Any) -> Iterator[dict[str, str]]:
    """Stream ``czech_derivations.csv`` (for the analysis stages)."""
    return jsonio.iter_csv(path)
