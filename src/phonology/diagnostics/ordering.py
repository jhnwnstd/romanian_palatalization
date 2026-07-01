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
from typing import TYPE_CHECKING, Final, Iterator, Sequence

from ..inventory import FeatureInventory
from ..pipeline import RulePipeline
from ..rules import ResolverFn, Rule
from ..segments import Word, segments_to_ipa
from .compare import compare_ipa

if TYPE_CHECKING:
    from ..analyses.romanian_palatalization import AnalysisProfile


class OrderingStrategy(StrEnum):
    DECLARED = "DECLARED"
    ADJACENT_SWAP = "ADJACENT_SWAP"
    PAIRWISE_SWAP = "PAIRWISE_SWAP"
    MOVABLE_PERM = "MOVABLE_PERM"
    FULL = "FULL"


@dataclass(frozen=True, slots=True)
class OrderingCandidate:
    """One ordering + the SR it produced."""

    permutation: tuple[int, ...]  # index-into-original-rules order
    rule_names: tuple[str, ...]  # names in the produced order
    sr: str  # rendered IPA
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
) -> Iterator[tuple[int, ...]]:
    """Yield permutation-index tuples according to the strategy."""
    identity: tuple[int, ...] = tuple(range(n))
    yield identity  # always try baseline first
    if strategy is OrderingStrategy.DECLARED:
        return
    seen: set[tuple[int, ...]] = {identity}
    count = 1
    if strategy is OrderingStrategy.ADJACENT_SWAP:
        for i in range(n - 1):
            swap = list(identity)
            swap[i], swap[i + 1] = swap[i + 1], swap[i]
            t = tuple(swap)
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
                swap = list(identity)
                swap[i], swap[j] = swap[j], swap[i]
                t = tuple(swap)
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
            swap = list(identity)
            for original, target in zip(movable, perm_of_movable):
                swap[original] = target
            t = tuple(swap)
            if t in seen:
                continue
            seen.add(t)
            yield t
            count += 1
            if limit is not None and count >= limit:
                return
        return
    if strategy is OrderingStrategy.FULL:
        for full_perm in permutations(range(n)):
            if full_perm in seen:
                continue
            seen.add(full_perm)
            yield full_perm
            count += 1
            if limit is not None and count >= limit:
                return
        return


def find_valid_orderings(
    rules: Sequence[Rule],
    ur: Word,
    expected_field: str | Sequence[str],
    resolve: ResolverFn,
    inventory: FeatureInventory,
    strategy: OrderingStrategy = OrderingStrategy.ADJACENT_SWAP,
    movable_indices: Sequence[int] | None = None,
    limit: int | None = None,
) -> OrderingSearchResult:
    """Search rule orderings for ones that produce ``expected_field``.

    Low-level API taking all pipeline pieces separately. Most callers
    should prefer :func:`search_orderings`, which takes an
    :class:`AnalysisProfile` and unpacks it for you.

    ``expected_field`` accepts a pipe-separated raw string OR a tuple
    of variants (see :func:`compare_ipa`).
    """
    rules_tuple = tuple(rules)
    n = len(rules_tuple)
    if strategy is OrderingStrategy.MOVABLE_PERM and not movable_indices:
        raise ValueError(
            "OrderingStrategy.MOVABLE_PERM requires a non-empty "
            "movable_indices argument (or use ADJACENT_SWAP / "
            "PAIRWISE_SWAP / FULL instead)"
        )
    base = inventory.base_segments()
    baseline_sr: str = ""
    baseline_matched = False
    successful: list[OrderingCandidate] = []
    closest: OrderingCandidate | None = None
    closest_dist = float("inf")
    tried = 0

    for perm in _permutations_for(strategy, n, movable_indices, limit):
        permuted = tuple(rules_tuple[i] for i in perm)
        pipeline = RulePipeline(rules=permuted, resolve=resolve)
        try:
            deriv = pipeline.derive(ur)
        except Exception:
            continue
        sr = segments_to_ipa(deriv.sr, base)
        cmp = compare_ipa(sr, expected_field)
        cand = OrderingCandidate(
            permutation=perm,
            rule_names=tuple(r.name for r in permuted),
            sr=sr,
            matched=cmp.matched,
            edit_distance=(
                int(cmp.edit_distance)
                if cmp.edit_distance != float("inf")
                else 10**9
            ),
        )
        if perm == tuple(range(n)):
            baseline_sr = sr
            baseline_matched = cmp.matched
        if cmp.matched:
            successful.append(cand)
        else:
            if cand.edit_distance < closest_dist:
                closest = cand
                closest_dist = cand.edit_distance
        tried += 1

    return OrderingSearchResult(
        strategy=strategy,
        orderings_tried=tried,
        baseline_matched=baseline_matched,
        baseline_sr=baseline_sr,
        successful=tuple(successful),
        closest_unsuccessful=None if baseline_matched else closest,
    )


def search_orderings(
    profile: "AnalysisProfile",
    ur: Word,
    attested: str | Sequence[str],
    *,
    strategy: OrderingStrategy = OrderingStrategy.ADJACENT_SWAP,
    movable_indices: Sequence[int] | None = None,
    limit: int | None = None,
) -> OrderingSearchResult:
    """Ergonomic wrapper: takes an :class:`AnalysisProfile` and
    unpacks its rules/inventory/resolver for you.

    Prefer this over :func:`find_valid_orderings` in new code.
    """
    inv = profile.inventory
    return find_valid_orderings(
        rules=profile.rules,
        ur=ur,
        expected_field=attested,
        resolve=lambda label: inv.segment(label),
        inventory=inv,
        strategy=strategy,
        movable_indices=movable_indices,
        limit=limit,
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
        lines.append(
            "  baseline ordering already succeeds — no re-order needed"
        )
        return "\n".join(lines)
    lines.append(
        f"  {len(result.successful)} alternative ordering(s) succeeded:"
    )
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
    "search_orderings",
)
