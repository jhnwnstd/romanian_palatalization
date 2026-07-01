"""Distance-to-working metric for mismatches.

Given a predicted SR and the attested plural IPA, produce a numeric
score that quantifies how close the derivation is to matching. Small
distances mean "one edit away from correct"; large distances mean
"the rules are fundamentally wrong for this stem". Sorting the
mismatch report by distance surfaces the near-misses — the rows most
likely to benefit from a small feature/rule tweak.

The score has three components (see :class:`DistanceWeights`):

  - **ipa_edit** — Levenshtein distance between predicted and the
    closest attested variant. Captures gross string-level distance.
  - **feature_edit** — number of feature values that differ at the
    stem-final position between predicted and attested (approximated
    by counting differing explicit values in the corresponding
    inventory segments). Captures how many featural facts are wrong.
  - **rule_edit** — cardinality of the symmetric difference between
    the set of palatalisation rules that fired vs. the set expected
    to fire. Captures whether the RIGHT rules are firing.

``total`` is a weighted sum. Featural edits are weighted highest
because they're the most actionable diagnostic (a patch or clear-set
edit); rule firings next; IPA edit distance last as a tiebreaker.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final, Iterable

from ..inventory import FeatureInventory
from ..pipeline import Derivation
from ..segments import Segment, Word, tokenize_ipa
from .compare import CompareStatus, compare_ipa, split_variants


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class DistanceWeights:
    """Weights combining the three distance components into ``total``.

    Featural distance weighted highest because it maps directly to
    an actionable fix (edit a patch, add/remove a clear-set entry).
    IPA edit distance is a tiebreaker only — pure transcription
    differences shouldn't dominate the ranking.
    """

    feature: float = 0.5
    rule: float = 0.3
    ipa: float = 0.2
    ipa_clamp: int = 20                    # cap raw IPA distance's contribution


DEFAULT_WEIGHTS: Final[DistanceWeights] = DistanceWeights()

# Sentinel for "attested data missing / unparseable"; consumers
# should treat this as "cannot be ranked", not as "very far".
UNKNOWN_DISTANCE: Final[float] = math.inf


# ---------------------------------------------------------------------------
# Score
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class DistanceScore:
    """Composite distance between a predicted and attested SR."""

    ipa_edit: float
    feature_edit: float
    rule_edit: int
    total: float


def _last_stem_consonant(word: Word) -> Segment | None:
    """Right-to-left search for a +Consonantal segment.

    Suffixes appended by g2p (i, e, u+r+i) are vowels or -Consonantal,
    so this walk lands on the last consonant of the stem — the
    position the palatalization rules were supposed to change.
    """
    for seg in reversed(word):
        if seg.features.get("Consonantal") == "+":
            return seg
    return None


def _attested_last_consonant_features(
    attested_ipa: str,
    inventory: FeatureInventory,
):
    """Return the feature map of the last +Consonantal inventory
    segment in the first attested variant, or None if none is
    findable in the inventory.
    """
    variants = split_variants(attested_ipa)
    if not variants:
        return None
    tokens = tokenize_ipa(variants[0])
    for tok in reversed(tokens):
        feats = inventory.base.get(tok)
        if feats is None:
            continue
        if feats.get("Consonantal") == "+":
            return feats
    return None


def feature_edit_distance(
    predicted_last_cons: Segment,
    attested_ipa: str,
    inventory: FeatureInventory,
) -> float:
    """Count feature-value differences between the predicted last-
    stem consonant and the last consonant of the attested plural.

    Returns ``UNKNOWN_DISTANCE`` if the attested last-consonant can't
    be located in the inventory (unusual characters, foreign IPA
    only, etc.). Callers must check for infinity before treating the
    score as comparable.
    """
    attested_feats = _attested_last_consonant_features(
        attested_ipa, inventory,
    )
    if attested_feats is None:
        return UNKNOWN_DISTANCE
    diff = 0
    for feat, val in attested_feats.items():
        if val == "0":
            continue
        pred_val = predicted_last_cons.features.get(feat, "0")
        if pred_val != "0" and pred_val != val:
            diff += 1
    return float(diff)


def rule_firing_distance(
    derivation: Derivation,
    expected_fired: Iterable[str],
    palatalization_rules: frozenset[str],
) -> int:
    """Symmetric difference between actually-fired palatalization
    rules and the set expected to fire.

    ``palatalization_rules`` is passed in (typically sourced from
    the analysis module's ``PALATALIZATION_RULE_NAMES``) so this
    module holds no hard-coded copy that could drift on a rename.
    """
    fired = {
        s.rule for s in derivation.steps
        if s.fired and s.rule in palatalization_rules
    }
    expected = set(expected_fired) & palatalization_rules
    return len(fired ^ expected)


def score(
    predicted_sr: str,
    attested_field: str,
    derivation: Derivation | None,
    expected_fired: Iterable[str],
    inventory: FeatureInventory,
    palatalization_rules: frozenset[str],
    weights: DistanceWeights = DEFAULT_WEIGHTS,
) -> DistanceScore:
    """Compute the composite distance score for one row.

    ``expected_fired`` is the rule-firing oracle for this row (the
    analysis module's ``expected_firings(row)`` — passed in so this
    module stays analysis-agnostic).
    """
    cmp = compare_ipa(predicted_sr, attested_field)
    if cmp.status is CompareStatus.MATCH:
        ipa: float = 0.0
    elif cmp.status is CompareStatus.MISMATCH:
        ipa = cmp.edit_distance
    else:
        ipa = UNKNOWN_DISTANCE

    feat: float = 0.0
    rule: int = 0
    if derivation is not None:
        last = _last_stem_consonant(derivation.sr)
        if last is not None and attested_field:
            feat = feature_edit_distance(last, attested_field, inventory)
        rule = rule_firing_distance(
            derivation, expected_fired, palatalization_rules,
        )

    if math.isinf(ipa) or math.isinf(feat):
        total = UNKNOWN_DISTANCE
    else:
        total = round(
            weights.feature * feat
            + weights.rule * rule
            + weights.ipa * min(ipa, weights.ipa_clamp),
            3,
        )
    return DistanceScore(
        ipa_edit=ipa, feature_edit=feat, rule_edit=rule, total=total,
    )


__all__: Final[tuple[str, ...]] = (
    "DEFAULT_WEIGHTS",
    "DistanceScore",
    "DistanceWeights",
    "UNKNOWN_DISTANCE",
    "feature_edit_distance",
    "rule_firing_distance",
    "score",
)
