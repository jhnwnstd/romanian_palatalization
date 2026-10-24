"""Corpus attestation of word forms (design rules R11 and R24).

Dictionary forms are not attestation: Wiktionary and MorfFlex generate
paradigm cells nobody uses (Czech biolozích, snězích: 0 tokens in
ces_news_2024_1M). A form is *attested* by its token count in the
language's Leipzig corpus, read from the corpus's words table with
:func:`lexpipe.frequency.count_words`.

Counts
------
:class:`AttestationTable` keeps both a case-sensitive table (``Medici``
and ``medici`` apart) and a folded one, and the corpus size used for
rates: ``total_tokens`` is the number of **word** tokens (tokens
:func:`lexpipe.frequency.clean_token` keeps: no punctuation, numbers
or URLs), after clitic splitting. A sentence-initial correction of
the case-sensitive counts needs the sentences file and is not made
here.

Italian elided clitics
----------------------
Leipzig tokenizes ``dell'amico`` as one token, which hides 27–73% of
the singular tokens of vowel-initial nouns. For a language whose
registry row lists ``elision_clitics`` (Italian:
:data:`lexpipe.spec.ITALIAN_ELISION_CLITICS`; Czech and Romanian have
none), :func:`split_elided` counts the part after the first
apostrophe (``'`` or ``’``) as a token of its own when the part before
it is a listed clitic (in any case) and the part after it is
alphabetic: ``l'amico`` → ``l`` + ``amico``, ``L’Italia`` → ``L`` +
``Italia``, ``c'è`` → ``c`` + ``è``; ``po'`` and ``dev'essere`` stay
whole. Both pieces count, so a split adds one token to the total.

Thresholds
----------
Every learner-evidence and TP table is reported over the thresholds
:data:`TAUS` (R11): ``dictionary`` (every dictionary form counts),
``ge1`` (at least one token), ``ge0.5pm`` and ``ge2pm`` (at least 0.5
and 2 tokens per million word tokens). Rates are compared exactly
(rational arithmetic), so a form exactly at a threshold is attested.
:func:`attested_count` applies a threshold to any precomputed count
(a sum over forms, a count with a sentence-initial correction): it is
the arithmetic :meth:`AttestationTable.attested` uses, so the two can
never disagree.

Paradigms
---------
Dictionaries generate paradigm cells nobody uses, so a generated
paradigm (a derivative's forms, a lexeme's cells) is attested form by
form. :func:`attest_paradigm` looks up every distinct form (forms
that share a lookup key, such as ``Medici``/``medici`` when counts are
folded, count once), tests each against every τ, and sums the counts
of the forms the caller does not exclude (forms that are also forms
of another lexeme, R11's ``form_ambiguous``): the paradigm is attested
at τ when that sum reaches τ.
"""

from __future__ import annotations

import unicodedata
from collections.abc import Callable, Collection, Iterable, Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from types import MappingProxyType
from typing import Final

from lexpipe import PROJECT_ROOT
from lexpipe.frequency import (
    TOKEN_TRIM,
    count_words,
    find_words_file,
    language_layout,
    normalize_key,
)
from lexpipe.spec import LanguageInfo

TAUS: Final[tuple[str, ...]] = ("dictionary", "ge1", "ge0.5pm", "ge2pm")
APOSTROPHES: Final[str] = "'’"
_MILLION: Final[int] = 1_000_000
_PER_MILLION: Final[Mapping[str, Fraction]] = MappingProxyType(
    {"ge0.5pm": Fraction(1, 2), "ge2pm": Fraction(2)}
)


def split_elided(token: str, clitics: Collection[str]) -> tuple[str, ...]:
    """Split an elided clitic off *token* (see the module docstring).

    Returns ``(clitic, word)`` or ``(token,)`` unchanged. *clitics*
    are lowercase and written without the apostrophe.
    """
    if not clitics:
        return (token,)
    text = unicodedata.normalize("NFC", token).strip().strip(TOKEN_TRIM)
    cut = next((i for i, ch in enumerate(text) if ch in APOSTROPHES), -1)
    if cut <= 0:
        return (token,)
    prefix, rest = text[:cut], text[cut + 1 :]
    if prefix.lower() in clitics and rest.isalpha():
        return (prefix, rest)
    return (token,)


def _check_tau(tau: str) -> str:
    if tau not in TAUS:
        raise ValueError(f"unknown threshold {tau!r}; known: {TAUS}")
    return tau


def attested_count(count: int, tau: str, total_tokens: int) -> bool:
    """Does a token *count* reach threshold *tau* in a corpus of
    *total_tokens* word tokens?

    ``dictionary`` is always true and ``ge1`` needs one token; the
    per-million thresholds compare ``count / total_tokens`` with 0.5
    or 2 per million exactly (``count · 10⁶ ≥ rate · total``), so a
    count exactly at a threshold is attested. Raises
    :class:`ValueError` for an unknown *tau*, a negative *count* or a
    *total_tokens* that is not positive.
    """
    _check_tau(tau)
    if count < 0:
        raise ValueError(f"negative token count: {count}")
    if total_tokens <= 0:
        raise ValueError(f"the corpus size must be positive: {total_tokens}")
    if tau == "dictionary":
        return True
    if tau == "ge1":
        return count >= 1
    return count * _MILLION >= _PER_MILLION[tau] * total_tokens


def attested_count_by_tau(
    count: int, total_tokens: int, taus: Sequence[str] = TAUS
) -> dict[str, bool]:
    """``{tau: attested_count(count, tau, total_tokens)}`` for *taus*."""
    return {tau: attested_count(count, tau, total_tokens) for tau in taus}


@dataclass(frozen=True, slots=True)
class AttestationTable:
    """Token counts of one corpus (see the module docstring).

    ``counts_cs`` is keyed by NFC forms with their case,
    ``counts_folded`` by :func:`lexpipe.frequency.normalize_key`
    forms; ``total_tokens`` is the corpus size in word tokens.
    """

    counts_cs: Mapping[str, int]
    counts_folded: Mapping[str, int]
    total_tokens: int

    def __post_init__(self) -> None:
        if self.total_tokens <= 0:
            raise ValueError("an attestation table needs a non-empty corpus")
        for name in ("counts_cs", "counts_folded"):
            frozen = MappingProxyType(dict(getattr(self, name)))
            object.__setattr__(self, name, frozen)

    @classmethod
    def from_words_file(
        cls, path: Path, *, clitics: Collection[str] = frozenset()
    ) -> AttestationTable:
        """Read a Leipzig words table, splitting *clitics* if given."""
        counts = count_words(
            path, split=lambda token: split_elided(token, clitics)
        )
        return cls(counts.case_sensitive, counts.folded, counts.word_tokens)

    @classmethod
    def for_language(
        cls, info: LanguageInfo, *, root: Path = PROJECT_ROOT
    ) -> AttestationTable:
        """The table of the language's Leipzig corpus, with its token
        rules (the extracted words table under
        ``data/<language>/leipzig/corpora/``)."""
        layout = language_layout(info.name, root)
        words = find_words_file(
            layout.corpus_dir(info.leipzig_id), info.leipzig_id
        )
        return cls.from_words_file(words, clitics=info.elision_clitics)

    def key(self, form: str, case_sensitive: bool = False) -> str:
        """The lookup key of *form*: NFC with its case, or folded."""
        if case_sensitive:
            return unicodedata.normalize("NFC", form).strip()
        return normalize_key(form)

    def count(self, form: str, case_sensitive: bool = False) -> int:
        """Tokens of *form* (0 if unattested)."""
        counts = self.counts_cs if case_sensitive else self.counts_folded
        return counts.get(self.key(form, case_sensitive), 0)

    def per_million(self, form: str, case_sensitive: bool = False) -> float:
        """Tokens of *form* per million word tokens."""
        return self.count(form, case_sensitive) * _MILLION / self.total_tokens

    def attested(
        self, form: str, tau: str, *, case_sensitive: bool = False
    ) -> bool:
        """Does *form* reach threshold *tau* (one of :data:`TAUS`)?"""
        if _check_tau(tau) == "dictionary":
            return True
        return attested_count(
            self.count(form, case_sensitive), tau, self.total_tokens
        )

    def attested_count(self, count: int, tau: str) -> bool:
        """:func:`attested_count` of *count* in this corpus."""
        return attested_count(count, tau, self.total_tokens)

    def attested_by_tau(
        self,
        form: str,
        taus: Sequence[str] = TAUS,
        *,
        case_sensitive: bool = False,
    ) -> dict[str, bool]:
        """``{tau: attested(form, tau)}`` for every threshold in *taus*."""
        return {
            tau: self.attested(form, tau, case_sensitive=case_sensitive)
            for tau in taus
        }


# ---------------------------------------------------------------------------
# Paradigms (see the module docstring)
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class FormAttestation:
    """One distinct form of a paradigm.

    ``key`` is the corpus lookup key (:meth:`AttestationTable.key`);
    ``counted`` is False for a form the caller excluded from the sum
    (its own ``count`` and ``attested_by_tau`` are still reported).
    """

    form: str
    key: str
    count: int
    attested_by_tau: Mapping[str, bool]
    counted: bool = True

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "attested_by_tau",
            MappingProxyType(dict(self.attested_by_tau)),
        )


@dataclass(frozen=True, slots=True)
class ParadigmAttestation:
    """The forms of a generated paradigm, attested one by one.

    ``forms`` holds one entry per distinct lookup key, in input order
    (the first spelling seen). ``total`` sums the counts of the
    counted forms; ``attested_by_tau`` tests that sum against each τ.
    """

    forms: tuple[FormAttestation, ...]
    total: int
    total_tokens: int
    attested_by_tau: Mapping[str, bool]

    def __post_init__(self) -> None:
        object.__setattr__(self, "forms", tuple(self.forms))
        object.__setattr__(
            self,
            "attested_by_tau",
            MappingProxyType(dict(self.attested_by_tau)),
        )

    @property
    def n_counted(self) -> int:
        """Forms that enter the sum."""
        return sum(1 for item in self.forms if item.counted)

    @property
    def excluded(self) -> tuple[str, ...]:
        """Forms left out of the sum, in input order."""
        return tuple(item.form for item in self.forms if not item.counted)

    def attested_forms(self, tau: str) -> tuple[str, ...]:
        """Counted forms that reach *tau* on their own, in input order."""
        _check_tau(tau)
        return tuple(
            item.form
            for item in self.forms
            if item.counted and item.attested_by_tau.get(tau, False)
        )


def _excluder(
    exclude: Callable[[str], bool] | Collection[str],
) -> Callable[[str], bool]:
    if callable(exclude):
        return exclude
    excluded = frozenset(exclude)
    return excluded.__contains__


def attest_paradigm(
    table: AttestationTable,
    forms: Iterable[str],
    *,
    taus: Sequence[str] = TAUS,
    case_sensitive: bool = False,
    exclude: Callable[[str], bool] | Collection[str] = (),
) -> ParadigmAttestation:
    """Attest every form of a generated paradigm (module docstring).

    Forms are looked up in *table* (folded, or with their case when
    *case_sensitive*); empty forms are skipped and forms sharing a
    lookup key count once. *exclude* names forms (a collection of
    spellings, or a predicate on the spelling) whose tokens must not
    be counted for this paradigm, typically forms that also belong to
    another lexeme. Every τ in *taus* is checked against the table's
    corpus size with :func:`attested_count`.
    """
    checked = tuple(_check_tau(tau) for tau in taus)
    skip = _excluder(exclude)
    seen: dict[str, FormAttestation] = {}
    for form in forms:
        if not form:
            continue
        key = table.key(form, case_sensitive)
        if key in seen:
            continue
        count = table.count(form, case_sensitive)
        seen[key] = FormAttestation(
            form=form,
            key=key,
            count=count,
            attested_by_tau=attested_count_by_tau(
                count, table.total_tokens, checked
            ),
            counted=not skip(form),
        )
    items = tuple(seen.values())
    total = sum(item.count for item in items if item.counted)
    return ParadigmAttestation(
        forms=items,
        total=total,
        total_tokens=table.total_tokens,
        attested_by_tau=attested_count_by_tau(
            total, table.total_tokens, checked
        ),
    )
