"""Rule-ordering search.

When a rule chain fails on a lemma, this module asks: is there a
DIFFERENT ordering of the same rules that would succeed? Answering
"yes" points at an ordering bug (the rule set is right but wrong
sequence). Answering "no" points at a rule-content bug (no
permutation of these rules produces the attested SR).

The search is parametrised by a strategy because full permutation
of nine rules is 362,880 — usable but wasteful. Strategies:

  - ``DECLARED`` — try only the declared ordering. Baseline check.
  - ``ADJACENT_SWAP`` — try every single adjacent swap of the
    declared ordering (N-1 permutations). Cheap and often
    diagnostic — if a rule needs to move one position, this finds it.
  - ``PAIRWISE_SWAP`` — try every pairwise swap (N(N-1)/2). Still
    fast at N=9 (36 permutations).
  - ``MOVABLE_PERM`` — permute only the rules whose position isn't
    logically constrained (declared via a ``movable_indices`` tuple
    on the analysis). Good middle-ground: the analyst names the
    subset they don't trust the order of.
  - ``FULL`` — every permutation of the rule tuple. Prohibitive at
    N > 8 unless bounded by ``limit``.

The success predicate is :func:`phonology.diagnostics.compare.compare`,
so an ordering "works" iff the SR it produces matches the attested
IPA under the shared normalisation stack — the same standard the
driver uses for the mismatch report.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from itertools import permutations
from typing import Callable, Final, Sequence

from ..inventory import FeatureInventory
from ..pipeline import RulePipeline
from ..rules import Rule
from ..segments import Segment, Word, segments_to_ipa
from .compare import compare


class OrderingStrategy(StrEnum):
    DECLARED = "DECLARED"
    ADJACENT_SWAP = "ADJACENT_SWAP"
    PAIRWISE_SWAP = "PAIRWISE_SWAP"
    MOVABLE_PERM = "MOVABLE_PERM"
    FULL = "FULL"


@dataclass(frozen=True, slots=True)
class OrderingCandidate:
    """One ordering + the SR it produced."""

    permutation: tuple[int, ...]     # index-into-original-rules order
    rule_names: tuple[str, ...]      # names in the produced order
    sr: str                          # rendered IPA
    matched: bool
    edit_distance: int


@dataclass(frozen=True, slots=True)
class OrderingSearchResult:
    """Result of an ordering search."""

    strategy: OrderingStrategy
    orderings_tried: int
    baseline_matched: bool
    baseline_sr: str
    successful: tuple[OrderingCandidate, ...]
    closest_unsuccessful: OrderingCandidate | None


def _permutations_for(
    strategy: OrderingStrategy,
    n: int,
    movable_indices: Sequence[int] | None,
    limit: int | None,
):
    """Yield permutation-index tuples according to the strategy."""
    identity = tuple(range(n))
    yield identity  # always try baseline first
    if strategy is OrderingStrategy.DECLARED:
        return
    seen = {identity}
    count = 1
    if strategy is OrderingStrategy.ADJACENT_SWAP:
        for i in range(n - 1):
            perm = list(identity)
            perm[i], perm[i + 1] = perm[i + 1], perm[i]
            t = tuple(perm)
            if t in seen:
                continue
            seen.add(t)
            yield t
            count += 1
            if limit is not None and count >= limit:
                return
        return
    if strategy is OrderingStrategy.PAIRWISE_SWAP:
        for i in range(n):
            for j in range(i + 1, n):
                perm = list(identity)
                perm[i], perm[j] = perm[j], perm[i]
                t = tuple(perm)
                if t in seen:
                    continue
                seen.add(t)
                yield t
                count += 1
                if limit is not None and count >= limit:
                    return
        return
    if strategy is OrderingStrategy.MOVABLE_PERM:
        if not movable_indices:
            return
        movable = tuple(movable_indices)
        for perm_of_movable in permutations(movable):
            perm = list(identity)
            for original, target in zip(movable, perm_of_movable):
                perm[original] = target
            t = tuple(perm)
            if t in seen:
                continue
            seen.add(t)
            yield t
            count += 1
            if limit is not None and count >= limit:
                return
        return
    if strategy is OrderingStrategy.FULL:
        for perm in permutations(range(n)):
            if perm in seen:
                continue
            seen.add(perm)
            yield perm
            count += 1
            if limit is not None and count >= limit:
                return
        return


def find_valid_orderings(
    rules: Sequence[Rule],
    ur: Word,
    expected_field: str,
    resolve: Callable[[str], Segment],
    inventory: FeatureInventory,
    strategy: OrderingStrategy = OrderingStrategy.ADJACENT_SWAP,
    movable_indices: Sequence[int] | None = None,
    limit: int | None = None,
) -> OrderingSearchResult:
    """Search rule orderings for ones that produce ``expected_field``.

    ``rules`` is the declared rule tuple; ``ur`` is the underlying
    representation to feed through the pipeline; ``expected_field`` is
    the attested IPA (pipe-separated variants ok). Success is defined
    by :func:`phonology.diagnostics.compare.compare` — the same
    matcher the driver's report uses.
    """
    rules_tuple = tuple(rules)
    n = len(rules_tuple)
    base = inventory.base_segments()
    baseline_sr: str = ""
    baseline_matched = False
    successful: list[OrderingCandidate] = []
    closest: OrderingCandidate | None = None
    closest_dist = 10**9
    tried = 0

    for perm in _permutations_for(strategy, n, movable_indices, limit):
        permuted = tuple(rules_tuple[i] for i in perm)
        pipeline = RulePipeline(rules=permuted, resolve=resolve)
        try:
            deriv = pipeline.derive(ur)
        except Exception:
            continue
        sr = segments_to_ipa(deriv.sr, base)
        cmp = compare(sr, expected_field)
        cand = OrderingCandidate(
            permutation=perm,
            rule_names=tuple(r.name for r in permuted),
            sr=sr,
            matched=cmp.matched,
            edit_distance=cmp.edit_distance,
        )
        if perm == tuple(range(n)):
            baseline_sr = sr
            baseline_matched = cmp.matched
        if cmp.matched:
            successful.append(cand)
        else:
            if cmp.edit_distance < closest_dist:
                closest = cand
                closest_dist = cmp.edit_distance
        tried += 1

    return OrderingSearchResult(
        strategy=strategy,
        orderings_tried=tried,
        baseline_matched=baseline_matched,
        baseline_sr=baseline_sr,
        successful=tuple(successful),
        closest_unsuccessful=None if baseline_matched else closest,
    )


def format_result(result: OrderingSearchResult) -> str:
    """Human-readable rendering of an OrderingSearchResult."""
    lines: list[str] = []
    lines.append(
        f"Strategy: {result.strategy.value}  "
        f"tried: {result.orderings_tried}  "
        f"baseline: /{result.baseline_sr}/ "
        f"{'✓' if result.baseline_matched else '✗'}"
    )
    if result.baseline_matched:
        lines.append("  baseline ordering already succeeds — no re-order needed")
        return "\n".join(lines)
    lines.append(f"  {len(result.successful)} alternative ordering(s) succeeded:")
    for cand in result.successful[:6]:
        lines.append(f"    {' → '.join(cand.rule_names)}")
    if len(result.successful) > 6:
        lines.append(f"    ... and {len(result.successful) - 6} more")
    if not result.successful and result.closest_unsuccessful:
        c = result.closest_unsuccessful
        lines.append(
            f"  no ordering matched; closest was distance "
            f"{c.edit_distance} → /{c.sr}/"
        )
        lines.append(f"    (via {' → '.join(c.rule_names)})")
    return "\n".join(lines)


__all__: Final[tuple[str, ...]] = (
    "OrderingCandidate",
    "OrderingSearchResult",
    "OrderingStrategy",
    "find_valid_orderings",
    "format_result",
)
