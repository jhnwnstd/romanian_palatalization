"""Distance-to-working metric for mismatches.

Given a predicted SR and the attested plural IPA, produce a numeric
score that quantifies how close the derivation is to matching. Small
distances mean "one edit away from correct"; large distances mean
"the rules are fundamentally wrong for this stem". Sorting the
mismatch report by distance surfaces the near-misses — the rows most
likely to benefit from a small feature/rule tweak.

The score has three components:

  - **ipa_edit** — Levenshtein distance between predicted and the
    closest attested variant. Captures gross string-level distance.
  - **feature_edit** — number of feature values that differ at the
    stem-final position between predicted and attested (approximated
    by counting differing explicit values in the corresponding
    inventory segments). Captures how many featural facts are wrong.
  - **rule_edit** — cardinality of the symmetric difference between
    the set of palatalisation rules that fired vs. the set expected
    to fire (inferred from ``stem_final × opportunity × cluster``).
    Captures whether the RIGHT rules are firing.

``total`` is a weighted sum. The weights are tuned to
prioritise featural corrections (which are actionable via a
patch/clear-set edit) over pure IPA edit distance (which is often
transcription noise).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final, Iterable

from ..inventory import FeatureInventory
from ..pipeline import Derivation
from ..segments import Segment, Word
from .compare import compare, split_variants


@dataclass(frozen=True, slots=True)
class DistanceScore:
    """Composite distance between a predicted and attested SR."""

    ipa_edit: int
    feature_edit: int
    rule_edit: int
    total: float


def ipa_edit_distance(predicted: str, attested_field: str) -> int:
    """Minimum Levenshtein between predicted and any attested variant.

    Returns 0 on match, positive integer otherwise. Uses the same
    normalisation stack as :func:`compare` — trailing glide / stress
    marks / etc. — so pure transcription differences don't inflate
    the score.
    """
    result = compare(predicted, attested_field)
    if result.matched:
        return 0
    return result.edit_distance


def feature_edit_distance(
    predicted_last_cons: Segment,
    attested_ipa: str,
    inventory: FeatureInventory,
) -> int:
    """Count feature-value differences between the predicted last-
    stem consonant and the last consonant of the attested plural.

    Uses inventory features when the attested consonant is a known
    inventory label; falls back to ``inf`` (encoded as a large int)
    for consonants not in the inventory.
    """
    from ..segments import tokenize_ipa
    variants = split_variants(attested_ipa)
    if not variants:
        return 999
    attested = variants[0]
    tokens = tokenize_ipa(attested)
    # Walk from right, find first +Consonantal segment.
    for tok in reversed(tokens):
        if tok not in inventory.base:
            continue
        feats = inventory.base[tok]
        if feats.get("Consonantal") != "+":
            continue
        # Count explicit-value differences.
        diff = 0
        for feat, val in feats.items():
            if val == "0":
                continue
            pred_val = predicted_last_cons.features.get(feat, "0")
            if pred_val != "0" and pred_val != val:
                diff += 1
        return diff
    return 999


def _last_stem_consonant(word: Word) -> Segment | None:
    """Right-to-left search for a +Consonantal segment (skip suffixes
    that are vowels or glides)."""
    for seg in reversed(word):
        if seg.features.get("Consonantal") == "+":
            return seg
    return None


PALATALIZATION_RULE_NAMES: Final[frozenset[str]] = frozenset({
    "dorsal-pal", "s-pal-rev", "assibilation", "bleed-rev",
})


def rule_firing_distance(
    derivation: Derivation,
    expected_fired: Iterable[str],
) -> int:
    """Symmetric difference between actually-fired palatalization
    rules and the set expected to fire.

    ``expected_fired`` is typically derived from ``mutation ×
    opportunity × cluster`` — e.g. for a t-final -i lemma we expect
    ``assibilation``; for a c-final -i we expect ``dorsal-pal``.
    """
    fired = {
        s.rule for s in derivation.steps
        if s.fired and s.rule in PALATALIZATION_RULE_NAMES
    }
    expected = set(expected_fired) & PALATALIZATION_RULE_NAMES
    return len(fired ^ expected)


def score(
    predicted_sr: str,
    attested_field: str,
    derivation: Derivation | None,
    expected_fired: Iterable[str],
    inventory: FeatureInventory,
) -> DistanceScore:
    """Compute the composite distance score for one row."""
    ipa = ipa_edit_distance(predicted_sr, attested_field)
    feat = 0
    rule = 0
    if derivation is not None:
        last = _last_stem_consonant(derivation.sr)
        if last is not None and attested_field:
            feat = feature_edit_distance(last, attested_field, inventory)
        rule = rule_firing_distance(derivation, expected_fired)
    # Weights: featural edits are the most actionable (patch or
    # clear-set), so weight them highest; rule firings next; ipa
    # distance as tiebreaker only.
    total = 0.5 * feat + 0.3 * rule + 0.2 * min(ipa, 20)
    return DistanceScore(
        ipa_edit=ipa,
        feature_edit=feat,
        rule_edit=rule,
        total=round(total, 3),
    )


__all__: Final[tuple[str, ...]] = (
    "DistanceScore",
    "PALATALIZATION_RULE_NAMES",
    "feature_edit_distance",
    "ipa_edit_distance",
    "rule_firing_distance",
    "score",
)
