"""Tolerance Principle rows with the semantics of the Romanian R analysis.

Yang's Tolerance Principle: a rule over N items is productive iff its
exceptions e satisfy e ≤ θ_N = N / ln N. The Romanian tables
(``analysis/romanian_tp_summary_*.csv``) come from
``tp_table()`` in ``analysis/analyze_romanian_palatalization.R``;
:func:`tp_row` reproduces one of their rows, column for column, so
Italian and Czech tables are counted the same way:

    type          the cell label
    N             mutated + non-mutated
    mutated       items that undergo the rule
    non-mutated   items that do not; always the exceptions e
    rate          mutated / N; missing when N = 0
    majority?     mutated > non-mutated; missing when N = 0 or on a tie
    tolerated?    e ≤ θ_N; missing when N ≤ 1 (θ_N undefined in R)
    memo          free text

Missing values are None here and empty cells once written with
:mod:`lexpipe.csvio` (the R files say ``NA``).

θ for small N: R's ``tp_threshold`` is NA for N ≤ 1, which is what
makes ``tolerated?`` missing there. :func:`theta` returns 0.0 for
N ≤ 1 (the convention of ``analysis/learnability_checks``) so that it
is always a number; decide tolerance with :func:`tolerated` or
:func:`tp_row`, not by comparing against ``theta(n)`` yourself.

Downsampling (:func:`top_n`) follows ``build_downsampled_lexicon()``
in the R script: a lemma's frequency is the maximum over its rows,
lemmas with zero frequency are dropped, and the n most frequent are
kept with ties broken by lemma (code-point order).

Harmonized long rows (R16–R19)
------------------------------
The cross-language tables count **root types** per target class ×
trigger cell with one definition for all three languages, from
:class:`RootObservation` records: a root, its target class, its token
frequency, and per trigger the :class:`Observation`\\s of its
realizations (does it alternate there; is the form attested at each
threshold τ; which allomorph: ``P`` palatalizing/front, ``N``
non-palatalizing, ``nonfront``). For a cell (:class:`CellSpec`: the
trigger, τ, doublet coding, …) every root gets one :class:`RootStatus`
from its realizations attested at τ:

    no_evidence  none attested in the cell
    shielded     all attested realizations use an N or non-front
                 allomorph (Romanian -uri, Czech -u/-ové/-ách)
    unknown      attested before the trigger, alternation unknown
    alternates   attested before the trigger, alternates
    nonalt       attested before the trigger, never alternates
    variable     a true doublet counted as an exception (below)
    excluded     a true doublet excluded from the domain (below)

and :func:`harmonized_row` reports, in :data:`TP_LONG_FIELDS`:

    N_target    target-final roots (minus excluded doublets)
    N_trigger   alternates + nonalt + variable
    e_nonalt    nonalt
    n_shielded  shielded
    e           exceptions to the row's rule: nonalt + variable for
                ``palatalize``, alternates + variable for ``retain``
    theta_trigger / theta_target   θ of N_trigger / N_target (empty
                when N ≤ 1)
    tolerated_trigger              e ≤ θ(N_trigger)
    tolerated_counterfactual       the check_claims.py counterfactual,
                generalized: if every target root met the trigger and
                behaved as the learner has evidence for, the rule's
                exceptions among N_target would be N_target −
                alternates (palatalize) or alternates + variable
                (retain)
    majority    N_trigger − e > e (missing on a tie or N_trigger = 0)

``e`` is not in the design's first list of columns: with doublets and
two rules, the exceptions of a row cannot be recovered from
``e_nonalt`` alone.

Both outcomes are evaluated (:func:`two_outcome_rows`, R18), and
nested rules in Yang's order (:func:`nested_rows`): each specific
subclass is evaluated with its own rule over the roots it contains
(and earlier subclasses did not claim), then the elsewhere rule over
the remainder. Roots are partitioned by their properties, never by
their outcome.

Doublets (R19). A root whose attested P realizations both alternate
and do not (chirurghi ~ chirurgi) is a doublet. :func:`classify_doublet`
codes it by corpus counts: a variant with ≥ 90% of the tokens is
dominant and decides the root under every coding; otherwise it is a
true doublet, and the cell's coding decides: ``dominant`` takes the
more frequent outcome (a tie counts as an exception, ``variable``),
``doublet_exception`` counts it as an exception to either rule
(``variable``), ``doublet_excluded`` removes it from the domain.

Two optional refinements, off by default (every row is then counted
as described above):

  - :attr:`Observation.cell` names the paradigm cell a realization
    belongs to. A doublet is competing variants of **one** cell, so
    the doublet coding runs within each cell; a root whose cells
    disagree (amici before -i but amiche before -e, pooled in one
    ``-i/-e`` trigger) is an exception to both rules, ``variable``,
    under every coding (a cell that is itself a true doublet keeps its
    coding: ``excluded`` wins, then ``variable``). Observations
    without a cell all count as one cell.
  - :attr:`CellSpec.doublet_min_tokens`: a doublet whose variants have
    fewer tokens in all is never dominant (:func:`classify_doublet`'s
    *min_total*); it is a true doublet, coded by the cell's coding.

Tracked rows
------------
Every language writes the long rows with three extra columns,
:data:`TP_LONG_EXTRA_FIELDS` (the header is
:data:`TP_LONG_TRACKED_FIELDS`): ``fragility`` (:func:`fragility` of
N_trigger and e: the fewest item changes that flip
``tolerated_trigger``), ``source`` and ``licence_class`` (the R9 gate
of :func:`lexpipe.csvio.write_tracked`). :func:`tracked_row` adds
them to a row.

Precomputed statuses
--------------------
A root's status depends only on its observations, the trigger, τ and
the doublet coding, not on the rule, the subset, the nested rule set
or the vocabulary it is counted in. Compute it once and reuse it:

  - :func:`cell_statuses` → :class:`CellStatuses`, every root's status
    in one cell (trigger × τ × coding × evidence set);
    :func:`cell_statuses_by_coding` does all codings in one pass
    (:func:`root_status_by_coding`).
  - :func:`harmonized_row`, :func:`two_outcome_rows`,
    :func:`nested_rows` and :func:`domain_counts` take
    ``statuses=``: the roots' statuses are then looked up, not
    recomputed (a :class:`CellStatuses` for another trigger, τ, coding
    or evidence set is refused).
  - :func:`status_counts`, :func:`status_row` and :func:`status_rows`
    evaluate a plain list of statuses (no root objects at all);
    :func:`counts_row` evaluates :class:`DomainCounts`. They share
    :func:`harmonized_row`'s arithmetic.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum
from types import MappingProxyType
from typing import Any, Final, TypeAlias, TypeVar

# R: SMALL_CELL_THRESHOLD; the TP binary is unreliable below it.
SMALL_CELL_THRESHOLD: Final[int] = 20

TP_FIELDS: Final[tuple[str, ...]] = (
    "type",
    "N",
    "mutated",
    "non-mutated",
    "rate",
    "majority?",
    "tolerated?",
    "memo",
)

TPRow: TypeAlias = dict[str, Any]
Row = TypeVar("Row", bound=Mapping[str, Any])


def theta(n: int) -> float:
    """θ_N = N / ln N, and 0.0 for N ≤ 1 (see the module docstring)."""
    return n / math.log(n) if n > 1 else 0.0


def majority(mutated: int, non_mutated: int) -> bool | None:
    """Does the rule apply to most items? None on a tie or an empty cell."""
    if mutated + non_mutated == 0 or mutated == non_mutated:
        return None
    return mutated > non_mutated


def tolerated(n: int, exceptions: int) -> bool | None:
    """e ≤ θ_N, or None when N ≤ 1 (R leaves θ_N undefined there)."""
    if n <= 1:
        return None
    return exceptions <= theta(n)


def small_n(n: int) -> bool:
    """True for cells below :data:`SMALL_CELL_THRESHOLD` (R's ``small_n``)."""
    return n < SMALL_CELL_THRESHOLD


def tp_row(label: str, n: int, mutated: int, memo: str = "") -> TPRow:
    """One TP summary row (columns :data:`TP_FIELDS`).

    *n* is the cell size and *mutated* the number of items that
    undergo the rule; the rest are the exceptions.
    """
    if n < 0 or not 0 <= mutated <= n:
        raise ValueError(
            f"{label!r}: need 0 <= mutated <= N, got {mutated}/{n}"
        )
    non_mutated = n - mutated
    return {
        "type": label,
        "N": n,
        "mutated": mutated,
        "non-mutated": non_mutated,
        "rate": mutated / n if n > 0 else None,
        "majority?": majority(mutated, non_mutated),
        "tolerated?": tolerated(n, non_mutated),
        "memo": memo,
    }


def header_row(label: str) -> TPRow:
    """Section header row (R's ``header_row``): a label, nothing else."""
    row: TPRow = dict.fromkeys(TP_FIELDS)
    row["type"] = label
    row["memo"] = ""
    return row


def blank_row() -> TPRow:
    """Separator row between sections (R's ``blank_row``)."""
    return header_row("")


def _column(name: str) -> Callable[[Mapping[str, Any]], Any]:
    def get(row: Mapping[str, Any]) -> Any:
        return row.get(name)

    return get


def _as_frequency(value: Any) -> float:
    """Frequency cell → number; blank, None and non-finite become 0."""
    if value is None or (isinstance(value, str) and not value.strip()):
        return 0.0
    number = float(value)
    return number if math.isfinite(number) else 0.0


def top_n(
    rows: Iterable[Row],
    n: int,
    freq_key: str | Callable[[Row], float],
    *,
    lemma_key: str = "lemma",
) -> list[Row]:
    """Rows of the *n* most frequent lemmas, most frequent first.

    *freq_key* is a column name (blank cells count as 0) or a function
    of the row. A lemma's frequency is the maximum over its rows;
    lemmas with frequency 0 are never selected; ties are broken by
    lemma. Every row of a selected lemma is returned, in input order
    within the lemma (the R script then keeps one row per lemma;
    deduplicate afterwards if you need that).
    """
    materialised = list(rows)
    if n <= 0:
        return []
    frequency_of: Callable[[Row], Any] = (
        _column(freq_key) if isinstance(freq_key, str) else freq_key
    )
    best: dict[str, float] = {}
    for row in materialised:
        lemma = str(row.get(lemma_key) or "")
        freq = _as_frequency(frequency_of(row))
        best[lemma] = max(best.get(lemma, 0.0), freq)
    ranked = sorted(
        (lemma for lemma, freq in best.items() if freq > 0),
        key=lambda lemma: (-best[lemma], lemma),
    )
    rank = {lemma: i for i, lemma in enumerate(ranked[:n])}
    chosen = [
        row for row in materialised if str(row.get(lemma_key) or "") in rank
    ]
    chosen.sort(key=lambda row: rank[str(row.get(lemma_key) or "")])
    return chosen


# ---------------------------------------------------------------------------
# Harmonized long rows (R16-R19); see the module docstring
# ---------------------------------------------------------------------------

TP_LONG_FIELDS: Final[tuple[str, ...]] = (
    "lang",
    "table",
    "rule",
    "target_class",
    "trigger",
    "subset",
    "tau",
    "evidence_set",
    "coding",
    "freq_basis",
    "V",
    "N_target",
    "N_trigger",
    "e_nonalt",
    "n_shielded",
    "e",
    "theta_trigger",
    "theta_target",
    "tolerated_trigger",
    "tolerated_counterfactual",
    "majority",
    "memo",
)
# The extra columns every tracked long table carries, in this order
# (the Italian, Czech and Romanian plugins' EXTRA_FIELDS).
TP_LONG_EXTRA_FIELDS: Final[tuple[str, ...]] = (
    "fragility",
    "source",
    "licence_class",
)
TP_LONG_TRACKED_FIELDS: Final[tuple[str, ...]] = (
    TP_LONG_FIELDS + TP_LONG_EXTRA_FIELDS
)

RULES: Final[tuple[str, ...]] = ("palatalize", "retain")
ALLOMORPHS: Final[tuple[str, ...]] = ("P", "N", "nonfront")
CODINGS: Final[tuple[str, ...]] = (
    "dominant",
    "doublet_exception",
    "doublet_excluded",
)
TRUE_DOUBLET: Final[str] = "true_doublet"
DOMINANT_PREFIX: Final[str] = "dominant:"
# A variant with at least this share of a doublet's tokens is dominant.
DOUBLET_DOMINANCE: Final[tuple[int, int]] = (9, 10)

_ALTERNATING: Final[str] = "alternating"
_NONALTERNATING: Final[str] = "nonalternating"


def _frozen_map(mapping: Mapping[str, Any]) -> Mapping[str, Any]:
    return MappingProxyType(dict(mapping))


@dataclass(frozen=True, slots=True)
class Observation:
    """One realization of a root in a trigger cell.

    ``alternates`` is None when it cannot be told (or does not apply:
    N and non-front allomorphs). ``attested_by_tau`` maps a threshold
    (:data:`lexpipe.attest.TAUS`) to whether the form reaches it; a
    missing threshold counts as not attested. ``form`` and ``count``
    (corpus tokens) are only needed to code doublets; ``cell`` (the
    paradigm cell, optional) keeps realizations of different cells
    from being coded as one doublet (module docstring).
    """

    alternates: bool | None
    attested_by_tau: Mapping[str, bool]
    allomorph: str
    form: str = ""
    count: int = 0
    # The paradigm cell of the realization ("" = the trigger's one cell;
    # see the module docstring on doublets).
    cell: str = ""

    def __post_init__(self) -> None:
        if self.allomorph not in ALLOMORPHS:
            raise ValueError(
                f"allomorph must be one of {ALLOMORPHS}: {self.allomorph!r}"
            )
        if self.count < 0:
            raise ValueError(f"negative token count: {self.count}")
        object.__setattr__(
            self, "attested_by_tau", _frozen_map(self.attested_by_tau)
        )


@dataclass(frozen=True, slots=True)
class RootObservation:
    """A root type with its realizations, keyed by trigger cell."""

    root_id: str
    target_class: str
    frequency: float
    triggers: Mapping[str, Sequence[Observation]]

    def __post_init__(self) -> None:
        if not math.isfinite(self.frequency) or self.frequency < 0:
            raise ValueError(f"{self.root_id}: bad frequency {self.frequency}")
        frozen = {name: tuple(obs) for name, obs in self.triggers.items()}
        object.__setattr__(self, "triggers", MappingProxyType(frozen))


@dataclass(frozen=True, slots=True)
class CellSpec:
    """The labels of one harmonized cell; ``tau`` and ``coding`` also
    decide the counting."""

    lang: str
    table: str
    target_class: str
    trigger: str
    tau: str = "dictionary"
    evidence_set: str = "E0"
    coding: str = "dominant"
    freq_basis: str = "lemma"
    V: int | None = None
    # A doublet with fewer tokens than this is never dominant (R19;
    # 0 = no minimum, the default).
    doublet_min_tokens: int = 0

    def __post_init__(self) -> None:
        if self.coding not in CODINGS:
            raise ValueError(
                f"coding must be one of {CODINGS}: {self.coding!r}"
            )
        if self.doublet_min_tokens < 0:
            raise ValueError(
                f"doublet_min_tokens must be >= 0: {self.doublet_min_tokens}"
            )


class RootStatus(StrEnum):
    """What a cell knows about one root (see the module docstring)."""

    NO_EVIDENCE = "no_evidence"
    SHIELDED = "shielded"
    UNKNOWN = "unknown"
    ALTERNATES = "alternates"
    NONALT = "nonalt"
    VARIABLE = "variable"
    EXCLUDED = "excluded"


def classify_doublet(
    counts_by_variant: Mapping[str, int], min_total: int = 0
) -> str:
    """``"dominant:<variant>"`` if one variant has ≥ 90% of the tokens,
    else ``"true_doublet"`` (also when no variant has a token, or when
    the variants have fewer than *min_total* tokens in all)."""
    if not counts_by_variant:
        raise ValueError("a doublet needs at least one variant")
    if any(count < 0 for count in counts_by_variant.values()):
        raise ValueError(f"negative count in {dict(counts_by_variant)}")
    total = sum(counts_by_variant.values())
    variant, best = max(counts_by_variant.items(), key=lambda kv: kv[1])
    share_num, share_den = DOUBLET_DOMINANCE
    if (
        total > 0
        and total >= min_total
        and best * share_den >= share_num * total
    ):
        return DOMINANT_PREFIX + variant
    return TRUE_DOUBLET


def _doublet_status(
    known: Sequence[Observation], coding: str, min_tokens: int = 0
) -> RootStatus:
    counts = {_ALTERNATING: 0, _NONALTERNATING: 0}
    seen: set[tuple[str, bool]] = set()
    for obs in known:
        key = (obs.form, bool(obs.alternates))
        if obs.form and key in seen:
            continue  # the same form listed twice counts once
        seen.add(key)
        counts[
            _ALTERNATING if obs.alternates else _NONALTERNATING
        ] += obs.count
    verdict = classify_doublet(counts, min_tokens)
    if verdict == DOMINANT_PREFIX + _ALTERNATING:
        return RootStatus.ALTERNATES
    if verdict == DOMINANT_PREFIX + _NONALTERNATING:
        return RootStatus.NONALT
    if coding == "doublet_excluded":
        return RootStatus.EXCLUDED
    if (
        coding == "dominant"
        and counts[_ALTERNATING] != counts[_NONALTERNATING]
    ):
        if counts[_ALTERNATING] > counts[_NONALTERNATING]:
            return RootStatus.ALTERNATES
        return RootStatus.NONALT
    return RootStatus.VARIABLE


def _decided_or_doublet(
    root: RootObservation, trigger: str, tau: str
) -> RootStatus | list[list[Observation]]:
    """The status if the coding cannot change it; else the root's known
    front observations grouped by cell (both outcomes: a doublet, or
    cells that disagree)."""
    attested = [
        obs
        for obs in root.triggers.get(trigger, ())
        if obs.attested_by_tau.get(tau, False)
    ]
    if not attested:
        return RootStatus.NO_EVIDENCE
    front = [obs for obs in attested if obs.allomorph == "P"]
    if not front:
        return RootStatus.SHIELDED
    known = [obs for obs in front if obs.alternates is not None]
    if not known:
        return RootStatus.UNKNOWN
    outcomes = {obs.alternates for obs in known}
    if outcomes == {True}:
        return RootStatus.ALTERNATES
    if outcomes == {False}:
        return RootStatus.NONALT
    groups: dict[str, list[Observation]] = {}
    for obs in known:
        groups.setdefault(obs.cell, []).append(obs)
    if len(groups) > 1 and all(
        len({obs.alternates for obs in group}) == 1
        for group in groups.values()
    ):
        return RootStatus.VARIABLE  # cells disagree, none is a doublet
    return list(groups.values())


def _group_status(
    group: Sequence[Observation], coding: str, min_tokens: int
) -> RootStatus:
    outcomes = {obs.alternates for obs in group}
    if outcomes == {True}:
        return RootStatus.ALTERNATES
    if outcomes == {False}:
        return RootStatus.NONALT
    return _doublet_status(group, coding, min_tokens)


def _resolve(
    groups: Sequence[Sequence[Observation]], coding: str, min_tokens: int
) -> RootStatus:
    """The status of a root's cell groups under *coding* (module
    docstring: a doublet is coded within its cell; cells that disagree
    make the root variable)."""
    statuses = {_group_status(group, coding, min_tokens) for group in groups}
    if len(statuses) == 1:
        return next(iter(statuses))
    if RootStatus.EXCLUDED in statuses:
        return RootStatus.EXCLUDED
    return RootStatus.VARIABLE


def root_status(
    root: RootObservation,
    trigger: str,
    tau: str,
    coding: str,
    *,
    min_tokens: int = 0,
) -> RootStatus:
    """The status of *root* in the *trigger* cell at threshold *tau*
    (*min_tokens*: :attr:`CellSpec.doublet_min_tokens`)."""
    decided = _decided_or_doublet(root, trigger, tau)
    if isinstance(decided, RootStatus):
        return decided
    return _resolve(decided, coding, min_tokens)


def _check_coding(coding: str) -> str:
    if coding not in CODINGS:
        raise ValueError(f"coding must be one of {CODINGS}: {coding!r}")
    return coding


def root_status_by_coding(
    root: RootObservation,
    trigger: str,
    tau: str,
    codings: Sequence[str] = CODINGS,
    *,
    min_tokens: int = 0,
) -> dict[str, RootStatus]:
    """``{coding: root_status(root, trigger, tau, coding)}`` in one pass.

    Only a true doublet's status depends on the coding, so the
    attested realizations are filtered once.
    """
    wanted = [_check_coding(coding) for coding in codings]
    decided = _decided_or_doublet(root, trigger, tau)
    if isinstance(decided, RootStatus):
        return dict.fromkeys(wanted, decided)
    return {coding: _resolve(decided, coding, min_tokens) for coding in wanted}


@dataclass(frozen=True, slots=True)
class DomainCounts:
    """Root counts of one cell (see the module docstring)."""

    n_target: int
    n_alternating: int
    n_nonalt: int
    n_variable: int
    n_shielded: int

    @property
    def n_trigger(self) -> int:
        return self.n_alternating + self.n_nonalt + self.n_variable

    def exceptions(self, rule: str) -> int:
        """Exceptions to *rule* among the N_trigger roots."""
        if rule == "palatalize":
            return self.n_nonalt + self.n_variable
        if rule == "retain":
            return self.n_alternating + self.n_variable
        raise ValueError(f"rule must be one of {RULES}: {rule!r}")

    def counterfactual_exceptions(self, rule: str) -> int:
        """Exceptions to *rule* among N_target in the counterfactual."""
        if rule == "palatalize":
            return self.n_target - self.n_alternating
        return self.exceptions(rule)


def _domain(statuses: Mapping[RootStatus, int], n_roots: int) -> DomainCounts:
    """DomainCounts from status tallies over *n_roots* roots."""
    return DomainCounts(
        n_target=n_roots - statuses.get(RootStatus.EXCLUDED, 0),
        n_alternating=statuses.get(RootStatus.ALTERNATES, 0),
        n_nonalt=statuses.get(RootStatus.NONALT, 0),
        n_variable=statuses.get(RootStatus.VARIABLE, 0),
        n_shielded=statuses.get(RootStatus.SHIELDED, 0),
    )


def status_counts(statuses: Iterable[RootStatus | str]) -> DomainCounts:
    """DomainCounts of one cell from one status per root.

    Statuses may be :class:`RootStatus` members or their string values
    (as read back from a CSV); anything else raises
    :class:`ValueError`. Each item is one root: a duplicate cannot be
    detected here, so pass every root once.
    """
    tallies: Counter[RootStatus] = Counter()
    n_roots = 0
    for status in statuses:
        tallies[RootStatus(status)] += 1
        n_roots += 1
    return _domain(tallies, n_roots)


@dataclass(frozen=True, slots=True)
class CellStatuses:
    """Every root's status in one cell, computed once and reused.

    The status of a root depends on its observations (and so on the
    evidence set they were built for), the trigger, τ, the doublet
    coding and the doublet minimum: those are recorded, and
    :meth:`check` refuses a :class:`CellSpec` that differs in any of
    them. ``statuses`` maps root_id to status for every root given to
    :func:`cell_statuses`, ``no_evidence`` ones included.
    """

    trigger: str
    tau: str
    coding: str
    evidence_set: str
    statuses: Mapping[str, RootStatus]
    doublet_min_tokens: int = 0

    def __post_init__(self) -> None:
        _check_coding(self.coding)
        object.__setattr__(
            self,
            "statuses",
            MappingProxyType(
                {
                    key: RootStatus(value)
                    for key, value in self.statuses.items()
                }
            ),
        )

    def matches(self, cell: CellSpec) -> bool:
        """Same trigger, τ, coding, evidence set and doublet minimum as
        *cell*?"""
        return (
            self.trigger,
            self.tau,
            self.coding,
            self.evidence_set,
            self.doublet_min_tokens,
        ) == (
            cell.trigger,
            cell.tau,
            cell.coding,
            cell.evidence_set,
            cell.doublet_min_tokens,
        )

    def check(self, cell: CellSpec) -> None:
        """Raise :class:`ValueError` unless :meth:`matches` *cell*."""
        if not self.matches(cell):
            raise ValueError(
                "statuses computed for trigger/tau/coding/evidence "
                f"{(self.trigger, self.tau, self.coding, self.evidence_set)}"
                " cannot count the cell "
                f"{(cell.trigger, cell.tau, cell.coding, cell.evidence_set)}"
            )

    def status(self, root_id: str) -> RootStatus:
        """The status of one root; :class:`ValueError` if it was not
        among the roots the statuses were computed for."""
        try:
            return self.statuses[root_id]
        except KeyError:
            raise ValueError(
                f"no precomputed status for root {root_id!r}"
            ) from None

    def of(self, root_ids: Iterable[str]) -> list[RootStatus]:
        """The statuses of *root_ids*, in their order."""
        return [self.status(root_id) for root_id in root_ids]

    def counts(self, root_ids: Iterable[str]) -> DomainCounts:
        """DomainCounts over *root_ids* (each counted once; a repeated
        id raises :class:`ValueError`)."""
        ids = list(root_ids)
        if len(set(ids)) != len(ids):
            raise ValueError("a root is listed twice")
        return status_counts(self.of(ids))


def _unique(roots: Iterable[RootObservation]) -> list[RootObservation]:
    listed = list(roots)
    seen: set[str] = set()
    for root in listed:
        if root.root_id in seen:
            raise ValueError(f"root {root.root_id!r} listed twice")
        seen.add(root.root_id)
    return listed


def cell_statuses(
    roots: Iterable[RootObservation], cell: CellSpec
) -> CellStatuses:
    """Every root's status in *cell*. Duplicate root_ids raise."""
    return CellStatuses(
        trigger=cell.trigger,
        tau=cell.tau,
        coding=cell.coding,
        evidence_set=cell.evidence_set,
        statuses={
            root.root_id: root_status(
                root,
                cell.trigger,
                cell.tau,
                cell.coding,
                min_tokens=cell.doublet_min_tokens,
            )
            for root in _unique(roots)
        },
        doublet_min_tokens=cell.doublet_min_tokens,
    )


def cell_statuses_by_coding(
    roots: Iterable[RootObservation],
    cell: CellSpec,
    codings: Sequence[str] = CODINGS,
) -> dict[str, CellStatuses]:
    """:func:`cell_statuses` for each coding (``cell.coding`` is
    ignored), each root's realizations filtered once."""
    wanted = [_check_coding(coding) for coding in dict.fromkeys(codings)]
    per_root = {
        root.root_id: root_status_by_coding(
            root,
            cell.trigger,
            cell.tau,
            wanted,
            min_tokens=cell.doublet_min_tokens,
        )
        for root in _unique(roots)
    }
    return {
        coding: CellStatuses(
            trigger=cell.trigger,
            tau=cell.tau,
            coding=coding,
            evidence_set=cell.evidence_set,
            statuses={key: value[coding] for key, value in per_root.items()},
            doublet_min_tokens=cell.doublet_min_tokens,
        )
        for coding in wanted
    }


def domain_counts(
    roots: Iterable[RootObservation],
    cell: CellSpec,
    *,
    statuses: CellStatuses | None = None,
) -> DomainCounts:
    """Count the roots of one cell. Duplicate root_ids raise.

    With *statuses* (computed for this cell's trigger, τ, coding and
    evidence set) each root's status is looked up, not recomputed.
    """
    if statuses is not None:
        statuses.check(cell)
    tallies: Counter[RootStatus] = Counter()
    seen: set[str] = set()
    for root in roots:
        if root.root_id in seen:
            raise ValueError(f"root {root.root_id!r} listed twice")
        seen.add(root.root_id)
        if statuses is None:
            status = root_status(
                root,
                cell.trigger,
                cell.tau,
                cell.coding,
                min_tokens=cell.doublet_min_tokens,
            )
        else:
            status = statuses.status(root.root_id)
        tallies[status] += 1
    return _domain(tallies, len(seen))


def _theta_or_none(n: int) -> float | None:
    return theta(n) if n > 1 else None


def counts_row(
    counts: DomainCounts,
    cell: CellSpec,
    *,
    rule: str = "palatalize",
    subset: str = "all",
    memo: str = "",
) -> TPRow:
    """One :data:`TP_LONG_FIELDS` row for *rule* from precomputed counts
    (the arithmetic of :func:`harmonized_row`)."""
    e = counts.exceptions(rule)
    n_trigger = counts.n_trigger
    return {
        "lang": cell.lang,
        "table": cell.table,
        "rule": rule,
        "target_class": cell.target_class,
        "trigger": cell.trigger,
        "subset": subset,
        "tau": cell.tau,
        "evidence_set": cell.evidence_set,
        "coding": cell.coding,
        "freq_basis": cell.freq_basis,
        "V": cell.V,
        "N_target": counts.n_target,
        "N_trigger": n_trigger,
        "e_nonalt": counts.n_nonalt,
        "n_shielded": counts.n_shielded,
        "e": e,
        "theta_trigger": _theta_or_none(n_trigger),
        "theta_target": _theta_or_none(counts.n_target),
        "tolerated_trigger": tolerated(n_trigger, e),
        "tolerated_counterfactual": tolerated(
            counts.n_target, counts.counterfactual_exceptions(rule)
        ),
        "majority": majority(n_trigger - e, e),
        "memo": memo,
    }


def status_row(
    statuses: Iterable[RootStatus | str],
    cell: CellSpec,
    *,
    rule: str = "palatalize",
    subset: str = "all",
    memo: str = "",
) -> TPRow:
    """:func:`harmonized_row` from one precomputed status per root
    (:func:`status_counts`); *cell* only labels the row."""
    return counts_row(
        status_counts(statuses), cell, rule=rule, subset=subset, memo=memo
    )


def status_rows(
    statuses: Iterable[RootStatus | str],
    cell: CellSpec,
    *,
    subset: str = "all",
    memo: str = "",
) -> list[TPRow]:
    """:func:`two_outcome_rows` from precomputed statuses."""
    counts = status_counts(statuses)
    return [
        counts_row(counts, cell, rule=rule, subset=subset, memo=memo)
        for rule in RULES
    ]


def harmonized_row(
    roots: Iterable[RootObservation],
    cell: CellSpec,
    *,
    rule: str = "palatalize",
    subset: str = "all",
    memo: str = "",
    statuses: CellStatuses | None = None,
) -> TPRow:
    """One :data:`TP_LONG_FIELDS` row for *rule* over *roots*.

    *roots* are the target-final roots of the cell's target class (or
    of a subset of it); the caller selects them. *statuses* (from
    :func:`cell_statuses`) saves recomputing their statuses.
    """
    counts = domain_counts(roots, cell, statuses=statuses)
    return counts_row(counts, cell, rule=rule, subset=subset, memo=memo)


def two_outcome_rows(
    roots: Iterable[RootObservation],
    cell: CellSpec,
    *,
    subset: str = "all",
    memo: str = "",
    statuses: CellStatuses | None = None,
) -> list[TPRow]:
    """The ``palatalize`` and the ``retain`` row over the same roots."""
    counts = domain_counts(roots, cell, statuses=statuses)
    return [
        counts_row(counts, cell, rule=rule, subset=subset, memo=memo)
        for rule in RULES
    ]


@dataclass(frozen=True, slots=True)
class NestedRule:
    """A rule specific to a subclass, for :func:`nested_rows`."""

    subset: str
    rule: str
    member: Callable[[RootObservation], bool]

    def __post_init__(self) -> None:
        if self.rule not in RULES:
            raise ValueError(f"rule must be one of {RULES}: {self.rule!r}")


def nested_rows(
    roots: Iterable[RootObservation],
    cell: CellSpec,
    specific: Sequence[NestedRule],
    *,
    elsewhere: str = "palatalize",
    elsewhere_subset: str = "elsewhere",
    statuses: CellStatuses | None = None,
) -> list[TPRow]:
    """Rows of a nested rule set in Yang's order.

    Each specific rule is evaluated over the remaining roots it
    contains, in the order given (most specific first); the elsewhere
    rule over the roots no specific rule contains. *statuses* (from
    :func:`cell_statuses`) saves recomputing the roots' statuses.
    """
    remaining = list(roots)
    rows: list[TPRow] = []
    for nested in specific:
        inside = [nested.member(root) for root in remaining]
        members = [r for r, yes in zip(remaining, inside, strict=True) if yes]
        remaining = [
            r for r, yes in zip(remaining, inside, strict=True) if not yes
        ]
        rows.append(
            harmonized_row(
                members,
                cell,
                rule=nested.rule,
                subset=nested.subset,
                statuses=statuses,
            )
        )
    rows.append(
        harmonized_row(
            remaining,
            cell,
            rule=elsewhere,
            subset=elsewhere_subset,
            statuses=statuses,
        )
    )
    return rows


# ---------------------------------------------------------------------------
# Fragility and the tracked columns
# ---------------------------------------------------------------------------


def fragility(n: int, e: int) -> int | None:
    """Fewest item changes that flip the verdict e ≤ θ_N (N fixed).

    A tolerated rule flips when e exceeds ⌊θ_N⌋; an untolerated one
    when e drops to ⌊θ_N⌋. None when there is no verdict (N ≤ 1) or
    no number of changes flips it (θ_N ≥ N, i.e. N = 2).
    """
    if n < 0 or not 0 <= e <= n:
        raise ValueError(f"need 0 <= e <= N, got e={e}, N={n}")
    verdict_now = tolerated(n, e)
    if verdict_now is None:
        return None
    limit = math.floor(theta(n))
    if verdict_now:
        return limit + 1 - e if limit + 1 <= n else None
    return e - limit


def tracked_row(
    row: Mapping[str, Any],
    *,
    source: str | Sequence[str],
    licence_class: str,
) -> TPRow:
    """*row* (a :data:`TP_LONG_FIELDS` row) with the columns of
    :data:`TP_LONG_EXTRA_FIELDS`: ``fragility`` of its N_trigger and e,
    *source* (one code, or a list written ``|``-joined) and
    *licence_class*. The result has every :data:`TP_LONG_TRACKED_FIELDS`
    column."""
    out = dict(row)
    out["fragility"] = fragility(int(row["N_trigger"]), int(row["e"]))
    out["source"] = source if isinstance(source, str) else list(source)
    out["licence_class"] = licence_class
    return out
