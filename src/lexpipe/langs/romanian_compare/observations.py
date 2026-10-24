"""Root records → :class:`lexpipe.tp.RootObservation` (R15, R16).

Trigger cells
-------------
A root's realizations are the forms in which its final target meets
a suffix vowel:

    E0   the plural (-i, -e, or -uri, whose /u/ follows the root);
    E1   E0 + the 2sg present and 3sg subjunctive of the noun's verb
         (``data/verb_paradigm_verdict.csv``);
    E2   E1 + the derived verbs and adjectives of the lexicon.

Each harmonized row is about one trigger: ``-i`` (/i/), ``-e`` (/e/)
or ``-i/-e`` (either front vowel). A realization becomes an
:class:`lexpipe.tp.Observation` of that trigger cell
(:func:`lexpipe.langs.romanian_compare.morphology.allomorph_for`):

    ``P``         the vowel after the root is one of the cell's
                  triggers; ``alternates`` = the form shows the
                  target's palatal reflex;
    ``N``         another front vowel (the cell is realized without
                  this trigger: a coronal's -e plural in the -i cell);
    ``nonfront``  a back vowel (-uri, -a, -ui, -os, ...).

So a root is *shielded* in a cell when every attested realization
avoids the cell's trigger, and *no_evidence* when none is attested
at τ (a plural-less noun, an unattested plural, a plural the legacy
alignment could not parse).

Attestation (R11): a form is attested at τ when the core
:class:`lexpipe.attest.AttestationTable` says so for its corpus count.
A derived verb counts its infinitive and past participle together.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Final

from lexpipe.attest import AttestationTable
from lexpipe.langs.romanian_compare import morphology as mo
from lexpipe.langs.romanian_compare.corpus import attestation_table
from lexpipe.langs.romanian_compare.lexicon import fold_key
from lexpipe.langs.romanian_compare.records import RootRecord, RootsFile
from lexpipe.tp import Observation, RootObservation

EVIDENCE_SETS: Final[tuple[str, ...]] = ("E0", "E1", "E2")
_PLURAL_VOWEL: Final[Mapping[str, str]] = {
    "-i": "i",
    "-e": "e",
    "-uri": "back",
}
_LEXEME_SUFFIX: Final[str] = "#lexeme"


@dataclass(frozen=True, slots=True)
class Realization:
    """A form in which the root's target meets a suffix vowel."""

    source: str
    form: str
    vowel: str
    palatal: bool | None
    count: int
    key: str


def realizations(record: RootRecord, evidence: str) -> list[Realization]:
    """The root's realizations under evidence set E0, E1 or E2."""
    if evidence not in EVIDENCE_SETS:
        raise ValueError(f"unknown evidence set {evidence!r}")
    out: list[Realization] = []
    vowel = _PLURAL_VOWEL.get(record.plural_allomorph)
    if vowel is not None:
        out.append(
            Realization(
                "plural",
                record.plural,
                vowel,
                record.plural_palatal,
                record.plural_count,
                fold_key(record.plural),
            )
        )
    if evidence in ("E1", "E2"):
        for cell in record.verb_cells:
            out.append(
                Realization(
                    f"verb:{cell.cell}",
                    cell.form,
                    cell.vowel,
                    cell.palatal,
                    cell.count,
                    fold_key(cell.form),
                )
            )
    if evidence == "E2":
        for derivative in record.derivatives:
            if derivative.vowel is None:
                continue
            key = fold_key(derivative.form)
            if derivative.kind == "verb":
                key += _LEXEME_SUFFIX
            out.append(
                Realization(
                    f"deriv:{derivative.kind}",
                    derivative.form,
                    derivative.vowel,
                    derivative.palatal,
                    derivative.attest_count,
                    key,
                )
            )
    return out


def attestation_counts(records: Iterable[RootRecord]) -> dict[str, int]:
    """Corpus count of every realization key of every evidence set."""
    counts: dict[str, int] = {}
    for record in records:
        for item in realizations(record, "E2"):
            counts[item.key] = item.count
        counts[fold_key(record.lemma)] = record.lemma_count
    return counts


class Attester:
    """``{tau: attested}`` of realizations, by the core attestation table."""

    def __init__(self, roots: RootsFile) -> None:
        self._table: AttestationTable = attestation_table(
            attestation_counts(roots.records), roots.total_tokens
        )
        self._cache: dict[str, Mapping[str, bool]] = {}

    @property
    def table(self) -> AttestationTable:
        return self._table

    def __call__(self, key: str) -> Mapping[str, bool]:
        found = self._cache.get(key)
        if found is None:
            found = self._table.attested_by_tau(key)
            self._cache[key] = found
        return found


def observations(
    items: Sequence[Realization],
    trigger: str,
    attest: Callable[[str], Mapping[str, bool]],
) -> tuple[Observation, ...]:
    """The realizations as observations of one trigger cell."""
    out: list[Observation] = []
    for item in items:
        allomorph = mo.allomorph_for(trigger, item.vowel)
        if allomorph is None:
            continue
        out.append(
            Observation(
                alternates=item.palatal if allomorph == "P" else None,
                attested_by_tau=attest(item.key),
                allomorph=allomorph,
                form=item.form,
                count=item.count,
            )
        )
    return tuple(out)


def root_observation(
    record: RootRecord,
    evidence: str,
    attest: Callable[[str], Mapping[str, bool]],
    *,
    basis: str = "lemma",
    with_cells: bool = True,
) -> RootObservation:
    """A root with its three trigger cells under an evidence set.

    ``with_cells=False`` gives a root with no cells at all (a noun of
    the sampling frame that is not in the counted class).
    """
    items = realizations(record, evidence) if with_cells else []
    return RootObservation(
        root_id=record.root_id,
        target_class=record.target_class or record.nde or "none",
        frequency=float(record.frequency(basis)),
        triggers={
            trigger: observations(items, trigger, attest)
            for trigger in mo.TRIGGERS
        },
    )
