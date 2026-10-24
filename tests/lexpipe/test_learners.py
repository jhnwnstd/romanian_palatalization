"""Tests for lexpipe.learners on synthetic roots (R18)."""

from __future__ import annotations

import math
import random
import sys
from collections import Counter
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lexpipe import learners, tp  # noqa: E402
from lexpipe.tp import (  # noqa: E402
    CellSpec,
    NestedRule,
    Observation,
    RootObservation,
)

CORPUS = {"dictionary": True, "ge1": True, "ge0.5pm": True, "ge2pm": True}
CELL = CellSpec(lang="it", table="nouns", target_class="k", trigger="-i")


def root(root_id: str, alternates: bool, frequency: float) -> RootObservation:
    observation = Observation(alternates, CORPUS, "P")
    return RootObservation(root_id, "k", frequency, {"-i": (observation,)})


def pool(alternating: int, non_alternating: int) -> list[RootObservation]:
    roots = [root(f"a{i:03d}", True, 1.0 + i) for i in range(alternating)]
    roots += [
        root(f"n{i:03d}", False, 1.0 + i) for i in range(non_alternating)
    ]
    return roots


def never(_: Sequence[RootObservation]) -> str:
    return "neither"


# ---------------------------------------------------------------------------
# Fragility
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("n", "e", "expected"),
    [
        (100, 21, 1),  # θ = 21.7: one more exception breaks it
        (100, 0, 22),
        (100, 22, 1),  # not tolerated: one fewer exception fixes it
        (100, 25, 4),
        (28, 0, 9),  # the Romanian top-1000 c + -i cell, θ = 8.4
        (2, 0, None),  # θ_2 > 2: nothing flips it
        (1, 0, None),  # no verdict for N <= 1
        (0, 0, None),
    ],
)
def test_fragility(n: int, e: int, expected: int | None) -> None:
    assert learners.fragility(n, e) == expected
    if expected is not None:
        before = tp.tolerated(n, e)
        changed = e + expected if before else e - expected
        assert tp.tolerated(n, changed) is (not before)
        nearly = e + expected - 1 if before else e - expected + 1
        assert tp.tolerated(n, nearly) is before


@pytest.mark.parametrize(("n", "e"), [(-1, 0), (3, 4), (3, -1)])
def test_fragility_rejects_impossible_cells(n: int, e: int) -> None:
    with pytest.raises(ValueError):
        learners.fragility(n, e)


# ---------------------------------------------------------------------------
# Sampling
# ---------------------------------------------------------------------------


def test_sampling_follows_token_frequency() -> None:
    roots = [
        root("heavy", True, 1000.0),
        root("light", True, 1.0),
        root("never", True, 0.0),
    ]
    drawn: Counter[str] = Counter()

    def record(vocabulary: Sequence[RootObservation]) -> str:
        drawn.update(r.root_id for r in vocabulary)
        return "neither"

    learners.simulate_learners(
        roots, record, V_values=(1,), n_samples=400, seed=3
    )
    assert drawn["never"] == 0
    assert drawn["heavy"] + drawn["light"] == 400
    assert drawn["heavy"] > 380  # p = 1000 / 1001


def test_vocabularies_are_nested_and_without_replacement() -> None:
    roots = pool(6, 6)
    seen: list[list[str]] = []

    def record(vocabulary: Sequence[RootObservation]) -> str:
        ids = [r.root_id for r in vocabulary]
        assert len(set(ids)) == len(ids)
        seen.append(ids)
        return "neither"

    learners.simulate_learners(
        roots, record, V_values=(3, 8), n_samples=25, seed=3
    )
    assert len(seen) == 50
    for small, large in zip(seen[::2], seen[1::2], strict=True):
        assert (len(small), len(large)) == (3, 8)
        assert large[:3] == small


def test_simulation_is_deterministic_and_order_free() -> None:
    roots = pool(30, 10)
    evaluate = learners.verdict(CELL)
    args = {"V_values": (10, 20), "n_samples": 200}
    first = learners.simulate_learners(roots, evaluate, seed=1, **args)
    again = learners.simulate_learners(roots[::-1], evaluate, seed=1, **args)
    other = learners.simulate_learners(roots, evaluate, seed=2, **args)
    assert first == again
    assert first != other
    for shares in first:
        assert sum(shares.counts.values()) == shares.n_samples == 200
        assert list(shares.shares) == list(learners.LEARNER_OUTCOMES)
        assert math.isclose(sum(shares.shares.values()), 1.0)


# ---------------------------------------------------------------------------
# Verdicts
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("alternating", "non_alternating", "outcome"),
    [
        (20, 0, "palatalize"),
        (18, 2, "palatalize"),  # e = 2 <= θ_20 = 6.7
        (0, 20, "retain"),
        (10, 10, "neither"),
        (3, 2, "palatalize"),  # both tolerated (N = 5): fewer exceptions
        (2, 2, "neither"),  # both tolerated, tied
    ],
)
def test_verdict(alternating: int, non_alternating: int, outcome: str) -> None:
    roots = pool(alternating, non_alternating)
    assert learners.verdict(CELL)(roots) == outcome
    shares = learners.simulate_learners(
        roots,
        learners.verdict(CELL),
        V_values=(len(roots),),
        n_samples=5,
        seed=0,
    )
    # V = the whole pool: every learner sees every root.
    assert shares[0].counts[outcome] == 5


def test_conditioned_rule() -> None:
    roots = pool(10, 10)
    specific = [
        NestedRule("n-roots", "retain", lambda r: r.root_id.startswith("n"))
    ]
    evaluate = learners.verdict(CELL, specific=specific)
    assert evaluate(roots) == "conditioned"
    # A nested set with an unproductive member is no conditioned rule.
    wrong = [NestedRule("n-roots", "palatalize", specific[0].member)]
    assert learners.verdict(CELL, specific=wrong)(roots) == "neither"


def test_simulation_arguments_are_checked() -> None:
    roots = pool(3, 0)
    with pytest.raises(ValueError, match="exceeds"):
        learners.simulate_learners(roots, never, (4,), 1, seed=0)
    with pytest.raises(ValueError, match="positive"):
        learners.simulate_learners(roots, never, (0,), 1, seed=0)
    with pytest.raises(ValueError, match="n_samples"):
        learners.simulate_learners(roots, never, (1,), 0, seed=0)
    with pytest.raises(ValueError, match="duplicate"):
        learners.simulate_learners([*roots, roots[0]], never, (1,), 1, seed=0)
    with pytest.raises(ValueError, match="not one of"):
        learners.simulate_learners(roots, lambda v: "maybe", (1,), 1, seed=0)
    with pytest.raises(TypeError):
        learners.simulate_learners(roots, never, (1,), 1)  # type: ignore[call-arg]


def test_default_sweep() -> None:
    assert learners.DEFAULT_V_VALUES == (500, 1000, 2000, 5000)
    assert learners.DEFAULT_SAMPLES == 1000


# ---------------------------------------------------------------------------
# Several evaluators on the same draws
# ---------------------------------------------------------------------------


def mixed_pool(n: int = 60, seed: int = 11) -> list[RootObservation]:
    """Roots of every status at once, with skewed frequencies."""
    rng = random.Random(seed)
    roots = []
    for i in range(n):
        kind = rng.choice(["alt", "alt", "non", "doublet", "shield", "none"])
        if kind == "alt":
            observations = (Observation(True, CORPUS, "P"),)
        elif kind == "non":
            observations = (Observation(False, CORPUS, "P"),)
        elif kind == "doublet":
            observations = (
                Observation(True, CORPUS, "P", f"f{i}a", rng.randint(0, 9)),
                Observation(False, CORPUS, "P", f"f{i}b", rng.randint(0, 9)),
            )
        elif kind == "shield":
            observations = (Observation(None, CORPUS, "N"),)
        else:
            observations = ()
        triggers = {"-i": observations} if observations else {}
        roots.append(
            RootObservation(f"r{i:03d}", "k", rng.paretovariate(1.2), triggers)
        )
    return roots


def nested_by_id(prefix_digits: str) -> NestedRule:
    return NestedRule(
        "odd", "retain", lambda r: r.root_id[-1] in prefix_digits
    )


EVALUATORS = {
    "dominant": learners.verdict(CELL),
    "exception": learners.verdict(
        CellSpec("it", "nouns", "k", "-i", coding="doublet_exception")
    ),
    "nested": learners.verdict(CELL, specific=[nested_by_id("13579")]),
    "never": never,
}


def test_simulate_all_equals_one_run_per_evaluator() -> None:
    roots = mixed_pool()
    args: dict[str, Any] = {"V_values": (10, 30, 10), "n_samples": 60}
    together = learners.simulate_all(roots, EVALUATORS, seed=4, **args)
    assert list(together) == list(EVALUATORS)
    for name, evaluate in EVALUATORS.items():
        alone = learners.simulate_learners(roots, evaluate, seed=4, **args)
        assert together[name] == alone
        assert [s.V for s in alone] == [10, 30]


@pytest.mark.parametrize("workers", [2, 3, 8])
def test_simulate_all_does_not_depend_on_workers(workers: int) -> None:
    roots = mixed_pool()
    args: dict[str, Any] = {"V_values": (10, 25), "n_samples": 40, "seed": 9}
    serial = learners.simulate_all(roots, EVALUATORS, **args)
    parallel = learners.simulate_all(
        roots[::-1], EVALUATORS, workers=workers, **args
    )
    assert parallel == serial
    assert list(parallel) == list(EVALUATORS)


def test_simulate_all_arguments() -> None:
    roots = pool(3, 0)
    assert learners.simulate_all(roots, {}, (1,), 1, seed=0) == {}
    with pytest.raises(ValueError, match="workers"):
        learners.simulate_all(roots, {"x": never}, (1,), 1, seed=0, workers=0)
    with pytest.raises(ValueError, match="exceeds"):
        learners.simulate_all(roots, {}, (4,), 1, seed=0)
    with pytest.raises(ValueError, match="'bad' returned 'maybe'"):
        learners.simulate_all(
            roots, {"x": never, "bad": lambda v: "maybe"}, (1,), 2, seed=0
        )
    # The error of a worker process reaches the caller.
    with pytest.raises(ValueError, match="not one of"):
        learners.simulate_all(
            roots,
            {"x": never, "bad": lambda v: "maybe"},
            (1,),
            2,
            seed=0,
            workers=2,
        )


def test_fragility_is_reexported() -> None:
    assert learners.fragility is tp.fragility


# ---------------------------------------------------------------------------
# Verdicts from precomputed statuses
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("coding", tp.CODINGS)
def test_status_verdict_equals_verdict(coding: str) -> None:
    roots = mixed_pool(80, seed=coding.__len__())
    spec_cell = CellSpec("it", "nouns", "k", "-i", coding=coding)
    statuses = tp.cell_statuses(roots, spec_cell)
    odd = nested_by_id("13579")
    low = NestedRule("low", "palatalize", lambda r: r.root_id < "r030")
    members = {
        "odd": frozenset(r.root_id for r in roots if odd.member(r)),
        "low": frozenset(r.root_id for r in roots if low.member(r)),
    }
    cases = [
        ((), [], "palatalize"),
        ((odd,), [("retain", members["odd"])], "palatalize"),
        (
            (odd, low),
            [("retain", members["odd"]), ("palatalize", members["low"])],
            "retain",
        ),
    ]
    rng = random.Random(3)
    outcomes: Counter[str] = Counter()
    for specific, status_specific, elsewhere in cases:
        slow = learners.verdict(
            spec_cell, specific=specific, elsewhere=elsewhere
        )
        cached = learners.verdict(
            spec_cell,
            specific=specific,
            elsewhere=elsewhere,
            statuses=statuses,
        )
        fast = learners.status_verdict(statuses, status_specific, elsewhere)
        # A plain mapping without the no_evidence roots works the same.
        sparse = {
            root_id: status
            for root_id, status in statuses.statuses.items()
            if status is not tp.RootStatus.NO_EVIDENCE
        }
        plain = learners.status_verdict(sparse, status_specific, elsewhere)
        for _ in range(150):
            vocabulary = rng.sample(roots, rng.randint(0, len(roots)))
            expected = slow(vocabulary)
            outcomes[expected] += 1
            assert cached(vocabulary) == expected
            assert fast(vocabulary) == expected
            assert plain(vocabulary) == expected
    # The comparison exercised more than one outcome.
    assert len(outcomes) >= 2


def test_status_verdict_simulates_like_verdict() -> None:
    roots = mixed_pool()
    statuses = tp.cell_statuses(roots, CELL)
    args: dict[str, Any] = {"V_values": (15, 40), "n_samples": 50, "seed": 2}
    assert learners.simulate_learners(
        roots, learners.status_verdict(statuses), **args
    ) == learners.simulate_learners(roots, learners.verdict(CELL), **args)


def test_status_verdict_checks_rules_and_cells() -> None:
    with pytest.raises(ValueError, match="rule"):
        learners.status_verdict({}, [("mutate", frozenset())])
    with pytest.raises(ValueError, match="rule"):
        learners.status_verdict({}, elsewhere="mutate")
    with pytest.raises(ValueError):
        learners.status_verdict({"r1": "sometimes"})  # type: ignore[dict-item]
    statuses = tp.cell_statuses(pool(2, 2), CELL)
    other = CellSpec("it", "nouns", "k", "-i", tau="ge1")
    with pytest.raises(ValueError, match="cannot count"):
        learners.verdict(other, statuses=statuses)
