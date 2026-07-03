#!/usr/bin/env python3
"""Unit tests for the phonology framework primitives.

These tests exercise the SEARCH primitive, unification / deletion
semantics, glide formation, and the underspecification mechanism —
independent of the Romanian analysis. If a test here fails, the bug
is in the framework, not in the ruleset.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from phonology.inventory import FeatureInventory, UnderspecifiedSegment
from phonology.rules import (
    DeletionRule,
    GlideFormationRule,
    UnificationRule,
)
from phonology.search import Direction, Search
from phonology.segments import Segment, tokenize_ipa
from phonology.analyses.romanian_palatalization import load_inventory


# ---------------------------------------------------------------------------
# Unification semantics
# ---------------------------------------------------------------------------


class TestUnificationSemantics:
    """Segment.unify: 0 accepts anything; explicit accepts same only."""

    def _seg(self, **feats: str) -> Segment:
        return Segment(label="x", features=dict(feats))

    def test_unify_into_zero_writes_value(self) -> None:
        s = self._seg(F="0")
        out = s.unify({"F": "+"})
        assert out is not None
        assert out.features["F"] == "+"

    def test_unify_matching_value_is_noop(self) -> None:
        s = self._seg(F="+")
        out = s.unify({"F": "+"})
        assert out is not None
        assert out.features["F"] == "+"

    def test_unify_conflicting_value_fails(self) -> None:
        s = self._seg(F="+")
        out = s.unify({"F": "-"})
        assert out is None

    def test_unify_supply_zero_is_wildcard(self) -> None:
        s = self._seg(F="+")
        out = s.unify({"F": "0"})
        assert out is not None
        assert out.features["F"] == "+"

    def test_delete_resets_named_features(self) -> None:
        s = self._seg(A="+", B="-", C="+")
        out = s.delete(frozenset({"A", "B"}))
        assert out.features["A"] == "0"
        assert out.features["B"] == "0"
        assert out.features["C"] == "+"


# ---------------------------------------------------------------------------
# Search / terminator breadth
# ---------------------------------------------------------------------------


class TestSearchTerminatorBreadth:
    """Broad terminator halts on next segment; narrow scans past
    non-members."""

    def _word(self, *specs: dict[str, str]) -> tuple[Segment, ...]:
        return tuple(
            Segment(label=str(i), features=spec)
            for i, spec in enumerate(specs)
        )

    def test_broad_terminator_halts_on_next(self) -> None:
        # word: target=0, blocker=1, real-trigger=2
        w = self._word(
            {"T": "+"},
            {"F": "+"},           # not in condition class
            {"C": "+"},           # would satisfy condition, but blocked
        )
        s = Search(
            direction=Direction.RIGHT,
            terminator={},
            condition={"C": "+"},
        )
        assert s.locate(w, 0) is None    # blocker halts scan and fails

    def test_broad_terminator_licenses_next_when_condition_holds(
        self,
    ) -> None:
        w = self._word(
            {"T": "+"},
            {"C": "+"},           # immediately satisfies condition
        )
        s = Search(
            direction=Direction.RIGHT,
            terminator={},
            condition={"C": "+"},
        )
        assert s.locate(w, 0) == 1

    def test_narrow_terminator_scans_past_non_members(self) -> None:
        w = self._word(
            {"T": "+"},
            {"X": "+"},           # not in terminator class — transparent
            {"C": "+"},           # in terminator class
        )
        s = Search(
            direction=Direction.RIGHT,
            terminator={"C": "+"},
            condition={"C": "+"},
        )
        assert s.locate(w, 0) == 2

    def test_leftward_broad_is_strict_left_adjacency(self) -> None:
        w = self._word(
            {"C": "+"},           # would satisfy condition
            {"X": "+"},           # intervener — halts leftward broad scan
            {"T": "+"},
        )
        s = Search(
            direction=Direction.LEFT,
            terminator={},
            condition={"C": "+"},
        )
        assert s.locate(w, 2) is None


# ---------------------------------------------------------------------------
# Rule kinds
# ---------------------------------------------------------------------------


class TestRuleKinds:
    """Each rule kind behaves as advertised."""

    def _seg(self, label: str, **feats: str) -> Segment:
        return Segment(label=label, features=dict(feats))

    def test_unification_rule_fires_when_target_and_trigger_match(
        self,
    ) -> None:
        target_seg = self._seg("t", A="0", B="+")
        trigger_seg = self._seg("g", C="+")
        rule = UnificationRule(
            name="r",
            target={"B": "+"},
            supply={"A": "+"},
            search=Search(
                direction=Direction.RIGHT,
                terminator={},
                condition={"C": "+"},
            ),
        )
        result = rule.apply((target_seg, trigger_seg))
        assert result.edits
        assert result.word[0].features["A"] == "+"

    def test_unification_no_op_on_prespecified_conflict(self) -> None:
        target_seg = self._seg("t", A="-", B="+")   # A fixed to -
        trigger_seg = self._seg("g", C="+")
        rule = UnificationRule(
            name="r",
            target={"B": "+"},
            supply={"A": "+"},        # supply +A conflicts with -A
            search=Search(
                direction=Direction.RIGHT,
                terminator={},
                condition={"C": "+"},
            ),
        )
        result = rule.apply((target_seg, trigger_seg))
        assert not result.edits
        assert result.word == (target_seg, trigger_seg)

    def test_deletion_rule_clears_features(self) -> None:
        target_seg = self._seg("t", A="+", B="+")
        cond_seg = self._seg("s", C="+")
        rule = DeletionRule(
            name="r",
            target={"A": "+"},
            clear=frozenset({"A", "B"}),
            search=Search(
                direction=Direction.LEFT,
                terminator={},
                condition={"C": "+"},
            ),
        )
        result = rule.apply((cond_seg, target_seg))
        assert result.edits
        assert result.word[1].features["A"] == "0"
        assert result.word[1].features["B"] == "0"

    def test_default_fill_rule_writes_into_open_features(self) -> None:
        # search=None means environment-free default fill; iterate over
        # word, unify supply into every matching target.
        rule = UnificationRule(
            name="default",
            target={"C": "+"},
            supply={"A": "-"},
            search=None,
        )
        word = (
            self._seg("x", C="+", A="0"),   # will get A=-
            self._seg("y", C="+", A="+"),   # conflict, unchanged
            self._seg("z", C="-", A="0"),   # target fails, unchanged
        )
        result = rule.apply(word)
        assert result.word[0].features["A"] == "-"
        assert result.word[1].features["A"] == "+"
        assert result.word[2].features["A"] == "0"

    def test_glide_formation_desyllabifies_word_final_i(self) -> None:
        c = self._seg("t", Consonantal="+", Syllabic="-")
        i = self._seg(
            "i",
            Syllabic="+",
            Consonantal="-",
            CORONAL="+",
            Anterior="-",
            Distributed="+",
        )
        rule = GlideFormationRule(
            name="glide",
            requires_preceding={"Consonantal": "+"},
        )
        result = rule.apply((c, i))
        assert result.edits
        assert result.word[1].label == "j"
        assert result.word[1].features["Syllabic"] == "-"

    def test_glide_formation_needs_preceding_consonant(self) -> None:
        v = self._seg("a", Syllabic="+", Consonantal="-")
        i = self._seg("i", Syllabic="+", Consonantal="-")
        rule = GlideFormationRule(
            name="glide",
            requires_preceding={"Consonantal": "+"},
        )
        result = rule.apply((v, i))
        assert not result.edits


# ---------------------------------------------------------------------------
# Underspecification / inventory patches
# ---------------------------------------------------------------------------


class TestInventoryAndUnderspec:
    """Patches and underspec declarations plug into the inventory
    without editing the JSON."""

    def test_underspec_clears_declared_features(self) -> None:
        inv = load_inventory()
        k = inv.segment("k")
        K = inv.segment("K")
        assert k.features["CORONAL"] == "-"
        assert K.features["CORONAL"] == "0"
        # Paper e:dor-inventory: /K/ = /k/ minus only Coronal. Dorsal
        # is *not* cleared — it stays +Dorsal from /k/, and
        # palatalisation leaves Dorsal untouched, so the derived
        # affricate is +Cor +Dor.
        assert K.features["DORSAL"] == "+"
        # Continuant is NOT cleared for /K/ — this is the reason the
        # palatalised result is [t͡ʃ] and not [ʃ].
        assert K.features["Continuant"] == "-"

    def test_patches_apply_to_base_segments(self) -> None:
        inv = load_inventory()
        # /i/ patched to postalveolar coronal
        i = inv.segment("i")
        assert i.features["CORONAL"] == "+"
        assert i.features["Anterior"] == "-"
        assert i.features["Distributed"] == "+"
        # /t/ patched to prespecified -Strident
        t = inv.segment("t")
        assert t.features["Strident"] == "-"

    def test_underspec_has_extra_labels_beyond_base(self) -> None:
        inv = load_inventory()
        # /Z/ intentionally absent: the paper treats /z/ as fully
        # specified and inalterable (latex.tex:170).
        for label in ("K", "G", "S", "T", "D"):
            assert label in inv.underspec
            assert label not in inv.base
        assert "Z" not in inv.underspec


# ---------------------------------------------------------------------------
# tokenize_ipa
# ---------------------------------------------------------------------------


class TestTokenizeIpa:
    """The IPA tokeniser recognises tie bars, diphthong ties, and
    single-codepoint segments."""

    def test_plain_string_splits_by_codepoint(self) -> None:
        assert tokenize_ipa("kat") == ("k", "a", "t")

    def test_tie_bar_binds_affricate(self) -> None:
        # t͡ʃ is three codepoints (t, tie, ʃ) but one segment.
        assert tokenize_ipa("t͡ʃe") == ("t͡ʃ", "e")

    def test_diphthong_tie_binds_to_base(self) -> None:
        # o̯ is two codepoints (o, U+032F).
        assert tokenize_ipa("bro̯aska") == (
            "b", "r", "o̯", "a", "s", "k", "a",
        )

    def test_stress_and_syllable_marks_are_dropped(self) -> None:
        assert tokenize_ipa("ˈka.t") == ("k", "a", "t")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
