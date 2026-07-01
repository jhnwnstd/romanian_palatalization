#!/usr/bin/env python3
"""Tests for the Logical Phonology layer.

Exercises the LP primitives independently of the Romanian analysis
so a failure here points at the framework's LP compliance, not at
whether Romanian happens to work.

Covers:
  - Natural-class subsumption (section 3 of the LP notes)
  - The "cannot target absence" invariant (section 2)
  - Unification / subtraction semantics (section 6)
  - Circle notation via NaturalClass.from_segment (section 5)
  - Bracket-type distinction: NaturalClass vs FeatureChange (section 4)
  - Rule application: LP totality (all matches from input state)
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from phonology.lp import (  # noqa: E402
    FeatureChange,
    NaturalClass,
    feature_change,
    is_consistent,
    natural_class,
    subtract,
    unify,
    valued_features,
)
from phonology.rules import (  # noqa: E402
    DeletionRule,
    SegmentDeletionRule,
    UnificationRule,
)
from phonology.search import Direction, Search  # noqa: E402
from phonology.segments import Segment  # noqa: E402


def _seg(label: str, **feats: str) -> Segment:
    return Segment(label=label, features=dict(feats))


# ---------------------------------------------------------------------------
# Section 3 — natural classes and subsumption
# ---------------------------------------------------------------------------


class TestNaturalClass:
    """N(C) = { σ | σ ⊇ C } as sets of valued features."""

    def test_universal_class_contains_everything(self) -> None:
        # N(∅) — the universal class — must contain any segment.
        universe = NaturalClass.universal()
        assert universe.contains(_seg("k", CORONAL="-", DORSAL="+"))
        assert universe.contains(_seg("i", Syllabic="+"))
        assert universe.contains(_seg("empty"))

    def test_more_features_smaller_class(self) -> None:
        # More features in the spec → smaller class. The [+Cor, +Ant]
        # class is contained in [+Cor].
        broad = natural_class(CORONAL="+")
        narrow = natural_class(CORONAL="+", Anterior="+")
        assert broad.subsumes(narrow)      # N(broad) ⊇ N(narrow)
        assert not narrow.subsumes(broad)

    def test_subsumption_via_contains(self) -> None:
        t = _seg("t", CORONAL="+", Anterior="+", Strident="-")
        s = _seg("s", CORONAL="+", Anterior="+", Strident="+")
        cor = natural_class(CORONAL="+")
        strident = natural_class(Strident="+")
        # Both are in [+Cor]; only /s/ is in [+Strident].
        assert cor.contains(t)
        assert cor.contains(s)
        assert not strident.contains(t)
        assert strident.contains(s)


class TestCannotTargetAbsence:
    """LP §2: absence is not a property; you cannot target it.

    An underspecified segment (F = "0") is NOT in the class of
    segments explicitly specified for F. This is the mechanism
    forcing rules to aim broadly and let the operator sort them out.
    """

    def test_underspecified_segment_not_in_explicit_class(self) -> None:
        # A segment open for Anterior is NOT in the [+Anterior] class.
        # This is the "cannot target absence" invariant: even though
        # /T/ might LATER receive +Anterior via unification, the class
        # membership check on the CURRENT feature set says no.
        T = _seg("T", CORONAL="+", Anterior="0", Strident="0")
        plus_ant = natural_class(CORONAL="+", Anterior="+")
        minus_ant = natural_class(CORONAL="+", Anterior="-")
        assert not plus_ant.contains(T)
        assert not minus_ant.contains(T)

    def test_you_cannot_write_a_targets_absence_pattern(self) -> None:
        # There's no NaturalClass constructor that says "match iff F is
        # absent". Trying to reach the underspecified /T/ specifically
        # while excluding /t/ requires a natural-class specification
        # that doesn't exist — this test enforces the LP prediction
        # that you have to aim broadly and let the operation sort.
        T = _seg("T", CORONAL="+", Anterior="0")
        t = _seg("t", CORONAL="+", Anterior="+")
        # Any class containing T also contains t under subsumption:
        for spec in ({"CORONAL": "+"}, {}):
            cls = NaturalClass(spec=spec)
            if cls.contains(T):
                assert cls.contains(t), (
                    "no natural class can carve out /T/ without also "
                    "including /t/"
                )


# ---------------------------------------------------------------------------
# Section 6 — operators from Ω
# ---------------------------------------------------------------------------


class TestUnification:
    """LP §6b: three outcomes — fill, no-op, refuse."""

    def test_fill_open_slot(self) -> None:
        # {} ⊔ {+F} = {+F}
        s = _seg("T", Strident="0")
        result = unify(s, {"Strident": "+"})
        assert result is not None
        assert result.features["Strident"] == "+"

    def test_vacuous_when_already_present(self) -> None:
        # {+F} ⊔ {+F} = {+F}
        s = _seg("s", Strident="+")
        result = unify(s, {"Strident": "+"})
        assert result is not None
        assert result.features["Strident"] == "+"

    def test_refuse_on_contradiction(self) -> None:
        # {-F} ⊔ {+F} = ∅ (unification fails, LP: identity on segment)
        s = _seg("t", Strident="-")
        result = unify(s, {"Strident": "+"})
        assert result is None

    def test_inalterability_is_unification_failure(self) -> None:
        # /t/ prespecified -Strident is inalterable to a rule that
        # supplies +Strident. No annotation on /t/ says "exception"—
        # the operator itself refuses because there's no room.
        t = _seg("t", Strident="-")
        T = _seg("T", Strident="0")
        assibilate = {"Strident": "+"}
        assert unify(t, assibilate) is None       # prespecified: no
        assert unify(T, assibilate) is not None   # underspec: yes


class TestSubtraction:
    """LP §6a: A \\ {+F, -F} strips F at any value."""

    def test_strip_feature_at_any_value(self) -> None:
        # Subtract at "any value it bears" — the bleed's semantics.
        pos = _seg("x", Anterior="+", Distributed="-", Strident="+")
        result = subtract(pos, frozenset({"Anterior"}))
        assert result.features["Anterior"] == "0"
        assert result.features["Distributed"] == "-"

    def test_derived_underspecification(self) -> None:
        # Subtraction creates underspecification — this is what feeds
        # the default fill-in. LP predicts feature-changing rules can't
        # bypass an underspec because they must pass THROUGH one.
        s = _seg("s", Anterior="+", Distributed="-", Strident="+")
        stripped = subtract(
            s, frozenset({"Anterior", "Distributed", "Strident"}),
        )
        assert stripped.features["Anterior"] == "0"
        assert stripped.features["Distributed"] == "0"
        assert stripped.features["Strident"] == "0"

    def test_total_operation_never_fails(self) -> None:
        # Subtracting a feature the segment doesn't have is a no-op.
        s = _seg("bare", A="0")
        result = subtract(s, frozenset({"NonexistentFeature"}))
        assert result.features.get("NonexistentFeature", "0") == "0"


# ---------------------------------------------------------------------------
# Section 5 — circle notation
# ---------------------------------------------------------------------------


class TestCircleNotation:
    """[Ⓢ] = N(σ.features): the natural class defined by σ's bundle."""

    def test_from_fully_specified_segment_gives_singleton_style_class(
        self,
    ) -> None:
        # A fully specified segment's [Ⓢ] contains only that segment
        # (and any that carry MORE features, but not fewer).
        sigma = _seg(
            "t", CORONAL="+", Anterior="+", Distributed="-",
            Strident="-", Continuant="-",
        )
        cs = NaturalClass.from_segment(sigma)
        # sigma itself is in [Ⓢ]:
        assert cs.contains(sigma)
        # a segment lacking one of those features is NOT:
        lax = _seg("half", CORONAL="+", Anterior="+")
        assert not cs.contains(lax)

    def test_from_partially_specified_segment_contains_supersets(self) -> None:
        # Ⓢ = {+Cor, -Ant, +Dist}. The class contains any segment
        # subsuming those three features — including /tʃ/ (with the
        # extra Strident:+) and /ʃ/ (Strident:+, Continuant:+).
        S_derived = _seg(
            "postalv", CORONAL="+", Anterior="-", Distributed="+",
        )
        cs = NaturalClass.from_segment(S_derived)
        tsh = _seg(
            "tʃ", CORONAL="+", Anterior="-", Distributed="+",
            Strident="+", Continuant="-",
        )
        sh = _seg(
            "ʃ", CORONAL="+", Anterior="-", Distributed="+",
            Strident="+", Continuant="+",
        )
        # Both derived postalveolars are in the class defined by Ⓢ.
        assert cs.contains(tsh)
        assert cs.contains(sh)

    def test_absence_dropped_from_defining_bundle(self) -> None:
        # If σ is underspecified for F, [Ⓢ] doesn't constrain F.
        # This matches LP: only *valued* features participate.
        partial = _seg("x", CORONAL="+", Anterior="0")
        cs = NaturalClass.from_segment(partial)
        assert "Anterior" not in cs.spec
        assert cs.spec["CORONAL"] == "+"


# ---------------------------------------------------------------------------
# Section 4 — bracket type distinction
# ---------------------------------------------------------------------------


class TestBracketTypes:
    """[ ] and { } are different types; the API should distinguish."""

    def test_natural_class_and_feature_change_are_distinct(self) -> None:
        nc = natural_class(CORONAL="+", Strident="+")
        fc = feature_change(add={"Anterior": "-", "Distributed": "+"})
        # No accidental interchangeability:
        assert isinstance(nc, NaturalClass)
        assert isinstance(fc, FeatureChange)
        assert not isinstance(nc, FeatureChange)
        assert not isinstance(fc, NaturalClass)


# ---------------------------------------------------------------------------
# Consistency
# ---------------------------------------------------------------------------


class TestConsistency:
    def test_segment_with_only_plus_minus_zero_is_consistent(self) -> None:
        s = _seg("x", A="+", B="-", C="0")
        assert is_consistent(s)

    def test_valued_features_projects_out_absence(self) -> None:
        s = _seg("x", A="+", B="-", C="0")
        fs = valued_features(s)
        assert fs == frozenset({("+", "A"), ("-", "B")})


# ---------------------------------------------------------------------------
# LP totality — rules compute all firings from the input state
# ---------------------------------------------------------------------------


class TestLPTotality:
    """LP §9: within one rule, all licensed positions fire from the
    ORIGINAL input; nothing feeds or bleeds within one rule."""

    def test_unification_fires_on_all_matching_targets(self) -> None:
        # Two targets both matching, both licensed: BOTH fire in one
        # rule application (not just the leftmost).
        rule = UnificationRule(
            name="test-fill",
            target={"A": "0"},           # anything with A=0 will match
            supply={"A": "+"},
            search=None,                 # env-free default-fill
        )
        w = (
            _seg("x", A="0"),
            _seg("y", A="0"),
            _seg("z", A="+"),
        )
        result = rule.apply(w)
        assert result.word[0].features["A"] == "+"
        assert result.word[1].features["A"] == "+"
        assert result.word[2].features["A"] == "+"   # was already +
        # Two edits (z was already + so no edit for it):
        assert len(result.edits) == 2

    def test_search_looks_at_input_state_not_intermediate(self) -> None:
        # A search-driven rule where an earlier firing might have
        # changed the trigger IF the rule looked at intermediate
        # state. Under LP totality, the search must see the input.
        rule = UnificationRule(
            name="test",
            target={"A": "0"},
            supply={"A": "+"},
            search=Search(
                direction=Direction.RIGHT,
                terminator={},
                condition={"A": "0"},    # trigger licensed on OPEN A
            ),
        )
        # Two adjacent open-A segments: seg0's search sees seg1 open,
        # seg1's search would find nothing (end of word). Under LP
        # totality both fire because seg0's mutation doesn't retro-
        # actively close seg1's trigger.
        w = (_seg("a", A="0"), _seg("b", A="0"))
        result = rule.apply(w)
        # seg0 fires (search finds seg1's open A); seg1's search
        # finds nothing (no segment to its right).
        assert result.word[0].features["A"] == "+"
        # seg1 stays 0 because its search failed, not because seg0
        # fired first — a sequential impl might have made seg1's
        # trigger unavailable, but LP-simultaneous doesn't.
        assert result.word[1].features["A"] == "0"


# ---------------------------------------------------------------------------
# Segment deletion (LP § 7)
# ---------------------------------------------------------------------------


class TestSegmentDeletion:
    """LP §7: [ ... ] ↦ ε removes segments from the string."""

    def test_deletes_all_licensed_matches(self) -> None:
        rule = SegmentDeletionRule(
            name="delete-a",
            target={"Tag": "+"},
        )
        w = (
            _seg("a", Tag="+"),
            _seg("b", Tag="-"),
            _seg("c", Tag="+"),
            _seg("d", Tag="-"),
        )
        result = rule.apply(w)
        labels = [s.label for s in result.word]
        assert labels == ["b", "d"]
        assert len(result.edits) == 2

    def test_no_matches_leaves_word_unchanged(self) -> None:
        rule = SegmentDeletionRule(
            name="noop",
            target={"Tag": "+"},
        )
        w = (_seg("a", Tag="-"), _seg("b", Tag="-"))
        result = rule.apply(w)
        assert result.word == w
        assert result.edits == ()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
