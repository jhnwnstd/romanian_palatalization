"""Simulated learners and fragility of TP verdicts (design rule R18).

A single top-V lexicon cannot say whether speakers would split; cells
are small (the Romanian top-1000 c + -i cell has N = 28), so one or
two items can flip a verdict. This module asks two questions of a
set of :class:`lexpipe.tp.RootObservation` roots.

Simulated learners
------------------
:func:`simulate_learners` draws ``n_samples`` vocabularies of V root
types for every V in ``V_values`` and records which rule each
simulated learner would find productive (:data:`LEARNER_OUTCOMES`):
``palatalize``, ``retain``, ``neither`` or ``conditioned`` (no
pooled rule, but a nested rule set where every rule is productive).

Sampling scheme. Roots are drawn by token frequency without
replacement: each learner gives every root with positive frequency f
an independent key E/f with E ~ Exp(1) (``random.Random.expovariate``)
and takes the V roots with the smallest keys. This is the
exponential-clock form of Efraimidis & Spirakis's weighted sampling:
the order of the keys is the order of drawing roots one at a time,
each with probability proportional to its frequency among the roots
not yet drawn. One learner's vocabularies are nested (its V = 500
vocabulary lies inside its V = 1000 one: a lexicon that grows), and
roots of frequency 0 are never drawn. Keys are generated from
``random.Random(seed)`` with the roots in root_id order, so the result
depends only on the seed and the set of roots, not on their order.

The rule evaluation is a function from a vocabulary to an outcome;
:func:`verdict` builds the standard one from a cell
(:class:`lexpipe.tp.CellSpec`) and optional nested rules, and
:func:`status_verdict` the same one from precomputed root statuses
(:class:`lexpipe.tp.CellStatuses` or a ``{root_id: status}``
mapping), which is much faster: a root's status does not depend on
the vocabulary it is drawn into.

Several rules on the same draws
-------------------------------
:func:`simulate_all` evaluates several rule evaluations, by name, on
the **same** vocabularies: every evaluator sees exactly the draws
:func:`simulate_learners` would give it with the same seed, so its
shares equal a separate :func:`simulate_learners` call, while the
vocabularies are drawn once instead of once per evaluator. With
``workers > 1`` the evaluators are split into contiguous groups that
run in forked worker processes; each process redraws the same
vocabularies from the seed, and results are collected in submission
order, so the output is the same for any number of workers (R5).
Without the ``fork`` start method (Windows) it runs in-process.

Fragility
---------
:func:`fragility` is the minimum number of items whose status must
change (a rule-following item becoming an exception or the reverse,
N fixed) to flip the tolerance verdict e ≤ θ_N. It lives in
:mod:`lexpipe.tp` (:func:`lexpipe.tp.fragility`) and is re-exported
here.
"""

from __future__ import annotations

import concurrent.futures
import heapq
import multiprocessing
import random
from collections import Counter
from collections.abc import (
    Callable,
    Collection,
    Iterator,
    Mapping,
    Sequence,
)
from dataclasses import dataclass
from types import MappingProxyType
from typing import Any, Final, TypeAlias

from lexpipe import tp
from lexpipe.tp import CellSpec, NestedRule, RootObservation
from lexpipe.tp import fragility as fragility  # noqa: F401 (re-export)
from lexpipe.tp import tolerated

LEARNER_OUTCOMES: Final[tuple[str, ...]] = (
    "palatalize",
    "retain",
    "neither",
    "conditioned",
)
DEFAULT_V_VALUES: Final[tuple[int, ...]] = (500, 1000, 2000, 5000)
DEFAULT_SAMPLES: Final[int] = 1000

RuleEval: TypeAlias = Callable[[Sequence[RootObservation]], str]
# (rule, member root_ids) of one specific rule of status_verdict.
StatusRule: TypeAlias = tuple[str, Collection[str]]


@dataclass(frozen=True, slots=True)
class LearnerShares:
    """Outcome counts of the simulated learners for one V."""

    V: int
    n_samples: int
    counts: Mapping[str, int]

    @property
    def shares(self) -> dict[str, float]:
        """``counts`` divided by ``n_samples``, in outcome order."""
        return {
            outcome: self.counts[outcome] / self.n_samples
            for outcome in LEARNER_OUTCOMES
        }


# ---------------------------------------------------------------------------
# Drawing vocabularies
# ---------------------------------------------------------------------------


def _vocabularies(
    roots: Sequence[RootObservation],
    v_values: Sequence[int],
    rng: random.Random,
) -> dict[int, list[RootObservation]]:
    """One learner's nested vocabularies (see the module docstring)."""
    keyed = [(rng.expovariate(root.frequency), root) for root in roots]
    drawn = heapq.nsmallest(max(v_values), keyed, key=lambda kv: kv[0])
    order = [root for _, root in drawn]
    return {v: order[:v] for v in v_values}


@dataclass(frozen=True, slots=True)
class _Plan:
    """Checked arguments of a simulation: the pool, V values, samples."""

    pool: tuple[RootObservation, ...]
    v_values: tuple[int, ...]
    n_samples: int
    seed: int


def _plan(
    roots: Sequence[RootObservation],
    V_values: Sequence[int],
    n_samples: int,
    seed: int,
) -> _Plan:
    if n_samples <= 0:
        raise ValueError("n_samples must be positive")
    v_values = tuple(dict.fromkeys(V_values))
    if not v_values or min(v_values) <= 0:
        raise ValueError(f"V values must be positive: {list(V_values)}")
    ids = [root.root_id for root in roots]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate root_ids")
    pool = tuple(
        sorted(
            (root for root in roots if root.frequency > 0),
            key=lambda root: root.root_id,
        )
    )
    if max(v_values) > len(pool):
        raise ValueError(
            f"V = {max(v_values)} exceeds the {len(pool)} roots with "
            "positive frequency"
        )
    return _Plan(pool, v_values, n_samples, seed)


def _draws(plan: _Plan) -> Iterator[dict[int, list[RootObservation]]]:
    """Every learner's vocabularies, in sampling order."""
    rng = random.Random(plan.seed)
    for _ in range(plan.n_samples):
        yield _vocabularies(plan.pool, plan.v_values, rng)


Tallies: TypeAlias = dict[str, dict[int, dict[str, int]]]


def _tally(plan: _Plan, evaluators: Sequence[tuple[str, RuleEval]]) -> Tallies:
    """Outcome counts per evaluator name and V over the plan's draws."""
    tallies: Tallies = {
        name: {v: dict.fromkeys(LEARNER_OUTCOMES, 0) for v in plan.v_values}
        for name, _ in evaluators
    }
    for vocabularies in _draws(plan):
        for v, vocabulary in vocabularies.items():
            for name, evaluate in evaluators:
                outcome = evaluate(vocabulary)
                counts = tallies[name][v]
                if outcome not in counts:
                    raise ValueError(
                        f"evaluator {name!r} returned {outcome!r}, not one "
                        f"of {LEARNER_OUTCOMES}"
                    )
                counts[outcome] += 1
    return tallies


def _shares(
    plan: _Plan, tally: Mapping[int, Mapping[str, int]]
) -> list[LearnerShares]:
    return [
        LearnerShares(v, plan.n_samples, MappingProxyType(dict(tally[v])))
        for v in plan.v_values
    ]


def simulate_learners(
    roots: Sequence[RootObservation],
    rule_eval: RuleEval,
    V_values: Sequence[int] = DEFAULT_V_VALUES,
    n_samples: int = DEFAULT_SAMPLES,
    *,
    seed: int,
) -> list[LearnerShares]:
    """Outcome shares of simulated learners, one entry per V (in order).

    *rule_eval* returns one of :data:`LEARNER_OUTCOMES` for a
    vocabulary. *seed* has no default: pass the run context's
    (``ctx.seed``, R5). A V larger than the number of roots with
    positive frequency raises :class:`ValueError`, as do duplicate
    root_ids.
    """
    plan = _plan(roots, V_values, n_samples, seed)
    tallies = _tally(plan, [("rule_eval", rule_eval)])
    return _shares(plan, tallies["rule_eval"])


# Worker-process state of simulate_all (set by _init_worker in each
# forked child; the parent never sets it).
_WORKER: dict[str, Any] = {}


def _init_worker(
    plan: _Plan, evaluators: Sequence[tuple[str, RuleEval]]
) -> None:
    _WORKER["plan"] = plan
    _WORKER["evaluators"] = dict(evaluators)


def _worker_job(names: tuple[str, ...]) -> Tallies:
    evaluators: dict[str, RuleEval] = _WORKER["evaluators"]
    return _tally(
        _WORKER["plan"], [(name, evaluators[name]) for name in names]
    )


def _groups(names: Sequence[str], n_groups: int) -> list[tuple[str, ...]]:
    """*names* in *n_groups* contiguous groups of near-equal size."""
    size, extra = divmod(len(names), n_groups)
    groups: list[tuple[str, ...]] = []
    start = 0
    for index in range(n_groups):
        end = start + size + (1 if index < extra else 0)
        groups.append(tuple(names[start:end]))
        start = end
    return [group for group in groups if group]


def _fork_context() -> multiprocessing.context.BaseContext | None:
    try:
        return multiprocessing.get_context("fork")
    except ValueError:
        return None


def simulate_all(
    roots: Sequence[RootObservation],
    evaluators: Mapping[str, RuleEval],
    V_values: Sequence[int] = DEFAULT_V_VALUES,
    n_samples: int = DEFAULT_SAMPLES,
    *,
    seed: int,
    workers: int = 1,
) -> dict[str, list[LearnerShares]]:
    """Shares of every evaluator on the **same** simulated vocabularies.

    Returns ``{name: [LearnerShares per V]}`` in the order of
    *evaluators*; each list equals ``simulate_learners(roots,
    evaluators[name], V_values, n_samples, seed=seed)``. With
    *workers* > 1 the evaluators run in up to *workers* forked
    processes (module docstring); the result does not depend on
    *workers*. Arguments are checked as by :func:`simulate_learners`;
    an evaluator returning something other than a
    :data:`LEARNER_OUTCOMES` value raises :class:`ValueError`.
    """
    if workers < 1:
        raise ValueError(f"workers must be at least 1: {workers}")
    plan = _plan(roots, V_values, n_samples, seed)
    named = list(evaluators.items())
    if not named:
        return {}
    groups = _groups([name for name, _ in named], min(workers, len(named)))
    context = _fork_context() if len(groups) > 1 else None
    if context is None:
        tallies = _tally(plan, named)
    else:
        tallies = {}
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=len(groups),
            mp_context=context,
            initializer=_init_worker,
            initargs=(plan, named),
        ) as pool:
            futures = [pool.submit(_worker_job, group) for group in groups]
            for future in futures:
                tallies.update(future.result())
    return {name: _shares(plan, tallies[name]) for name, _ in named}


# ---------------------------------------------------------------------------
# Rule evaluations
# ---------------------------------------------------------------------------


def _pooled(counts: tp.DomainCounts) -> str | None:
    """The pooled verdict, or None when neither pooled rule holds."""
    n = counts.n_trigger
    e_pal = counts.exceptions("palatalize")
    e_ret = counts.exceptions("retain")
    pal = tolerated(n, e_pal) is True
    ret = tolerated(n, e_ret) is True
    if pal and ret:
        if e_pal == e_ret:
            return "neither"
        return "palatalize" if e_pal < e_ret else "retain"
    if pal:
        return "palatalize"
    if ret:
        return "retain"
    return None


def verdict(
    cell: CellSpec,
    *,
    specific: Sequence[NestedRule] = (),
    elsewhere: str = "palatalize",
    statuses: tp.CellStatuses | None = None,
) -> RuleEval:
    """The standard rule evaluation for :func:`simulate_learners`.

    The pooled rules come first: ``palatalize`` or ``retain`` if
    exactly one is tolerated (both can be, for N ≤ 7: then the one
    with fewer exceptions, ``neither`` on a tie). Otherwise, with
    *specific* nested rules, ``conditioned`` if every row of
    :func:`lexpipe.tp.nested_rows` is tolerated (a row with N ≤ 1 is
    not), else ``neither``. *statuses* (:func:`lexpipe.tp.cell_statuses`
    over every root that can be drawn) saves recomputing statuses.
    """
    if statuses is not None:
        statuses.check(cell)

    def evaluate(vocabulary: Sequence[RootObservation]) -> str:
        counts = tp.domain_counts(vocabulary, cell, statuses=statuses)
        pooled = _pooled(counts)
        if pooled is not None:
            return pooled
        if specific:
            rows = tp.nested_rows(
                vocabulary,
                cell,
                specific,
                elsewhere=elsewhere,
                statuses=statuses,
            )
            if all(row["tolerated_trigger"] is True for row in rows):
                return "conditioned"
        return "neither"

    return evaluate


def _check_rule(rule: str) -> str:
    if rule not in tp.RULES:
        raise ValueError(f"rule must be one of {tp.RULES}: {rule!r}")
    return rule


def status_verdict(
    statuses: tp.CellStatuses | Mapping[str, tp.RootStatus],
    specific: Sequence[StatusRule] = (),
    elsewhere: str = "palatalize",
) -> RuleEval:
    """:func:`verdict` from precomputed root statuses.

    *statuses* maps root_id to status (a :class:`lexpipe.tp.CellStatuses`
    or a plain mapping); a root missing from it counts as
    ``no_evidence``, so the mapping may leave those out. *specific*
    lists ``(rule, member root_ids)`` in Yang's order, the analogue of
    :class:`lexpipe.tp.NestedRule`; *elsewhere* is the rule of the
    remaining roots. Each vocabulary gives the same outcome as
    :func:`verdict` over the same roots.
    """
    table = (
        statuses.statuses
        if isinstance(statuses, tp.CellStatuses)
        else statuses
    )
    lookup: dict[str, tp.RootStatus] = {
        root_id: tp.RootStatus(status) for root_id, status in table.items()
    }
    rules = [
        (_check_rule(rule), frozenset(members)) for rule, members in specific
    ]
    _check_rule(elsewhere)

    def counts_of(ids: Sequence[str]) -> tp.DomainCounts:
        tallies = Counter(map(lookup.get, ids))
        return tp.DomainCounts(
            n_target=len(ids) - tallies[tp.RootStatus.EXCLUDED],
            n_alternating=tallies[tp.RootStatus.ALTERNATES],
            n_nonalt=tallies[tp.RootStatus.NONALT],
            n_variable=tallies[tp.RootStatus.VARIABLE],
            n_shielded=tallies[tp.RootStatus.SHIELDED],
        )

    def nested_holds(ids: list[str]) -> bool:
        remaining = ids
        for rule, members in rules:
            inside = [root_id for root_id in remaining if root_id in members]
            remaining = [
                root_id for root_id in remaining if root_id not in members
            ]
            part = counts_of(inside)
            if tolerated(part.n_trigger, part.exceptions(rule)) is not True:
                return False
        rest = counts_of(remaining)
        return tolerated(rest.n_trigger, rest.exceptions(elsewhere)) is True

    def evaluate(vocabulary: Sequence[RootObservation]) -> str:
        ids = [root.root_id for root in vocabulary]
        pooled = _pooled(counts_of(ids))
        if pooled is not None:
            return pooled
        if rules and nested_holds(ids):
            return "conditioned"
        return "neither"

    return evaluate
