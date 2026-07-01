#!/usr/bin/env python3
"""Tests for the diagnostics package.

Exercises compare, distance, explain, ordering, and perturbation
against controlled inputs so failures point at the diagnostic layer
rather than the underlying analysis.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from phonology.analyses.romanian_palatalization import (
    PALATALIZATION_RULE_NAMES,
    PATCHES,
    RULES_FROM_PAPER,
    UNDERSPEC,
    load_inventory,
)
from phonology.diagnostics.categorise import (
    RefinedCategory,
    refine_category,
)
from phonology.diagnostics.compare import Normalisation, compare
from phonology.diagnostics.distance import score as distance_score
from phonology.diagnostics.explain import (
    explain_derivation,
    format_explanation,
)
from phonology.diagnostics.ordering import (
    OrderingStrategy,
    find_valid_orderings,
)
from phonology.diagnostics.perturb import (
    PerturbationKind,
    search_perturbations,
)
from phonology.g2p import build_ur
from phonology.pipeline import RulePipeline
from phonology.search import SearchOutcome
from phonology.segments import segments_to_ipa


INVENTORY_JSON = PROJECT_ROOT / "romanian_features.json"


@pytest.fixture(scope="module")
def inventory():
    return load_inventory()


@pytest.fixture(scope="module")
def pipeline(inventory):
    return RulePipeline(
        rules=RULES_FROM_PAPER,
        resolve=lambda label: inventory.segment(label),
    )


# ---------------------------------------------------------------------------
# compare
# ---------------------------------------------------------------------------


class TestCompare:
    """The canonical predicted-vs-attested match stack."""

    def test_exact_match(self) -> None:
        r = compare("kat", "kat")
        assert r.matched
        assert r.normalisation_applied == ()
        assert r.edit_distance == 0

    def test_trailing_glide_j_i(self) -> None:
        r = compare("ruʃj", "ruʃi")
        assert r.matched
        assert Normalisation.TRAILING_GLIDE in r.normalisation_applied

    def test_unstressed_schwa_bidir(self) -> None:
        assert compare("abilitat͡si", "abilitət͡si").matched
        assert compare("abilitət͡si", "abilitat͡si").matched

    def test_ea_monophthong(self) -> None:
        r = compare("kɨʃlərease", "kɨʃlərese")
        assert r.matched

    def test_pipe_variants(self) -> None:
        r = compare("ruʃi", "fr:russe | ruʃi | ruʃʲ")
        assert r.matched
        assert r.attested_variant == "ruʃi"

    def test_far_string_no_match(self) -> None:
        r = compare("hello", "kat")
        assert not r.matched
        assert r.edit_distance >= 3


# ---------------------------------------------------------------------------
# distance
# ---------------------------------------------------------------------------


class TestDistance:
    def test_zero_on_exact_match(self, inventory, pipeline) -> None:
        ur = build_ur("muskə", "e", inventory, stem_final="s")
        deriv = pipeline.derive(ur)
        pred = segments_to_ipa(deriv.sr, inventory.base_segments())
        s = distance_score(
            predicted_sr=pred, attested_field="muʃte",
            derivation=deriv, expected_fired=frozenset({"dorsal-pal"}),
            inventory=inventory,
            palatalization_rules=PALATALIZATION_RULE_NAMES,
        )
        assert s.ipa_edit == 0
        assert s.total < 1.0

    def test_non_zero_on_mismatch(self, inventory, pipeline) -> None:
        ur = build_ur("muskə", "e", inventory, stem_final="s")
        deriv = pipeline.derive(ur)
        pred = segments_to_ipa(deriv.sr, inventory.base_segments())
        s = distance_score(
            predicted_sr=pred, attested_field="xyz",
            derivation=deriv, expected_fired=frozenset(),
            inventory=inventory,
            palatalization_rules=PALATALIZATION_RULE_NAMES,
        )
        assert s.ipa_edit > 0
        assert s.total > 0


# ---------------------------------------------------------------------------
# explain
# ---------------------------------------------------------------------------


class TestExplain:
    def test_explanation_covers_all_rules(self, inventory, pipeline) -> None:
        ur = build_ur("brad", "i", inventory, stem_final="d")
        exp = explain_derivation(pipeline, ur)
        # One trace per rule in the pipeline.
        assert len(exp.rule_traces) == len(RULES_FROM_PAPER)
        names = [t.rule_name for t in exp.rule_traces]
        assert "assibilation" in names
        assert "cor-default-continuant" in names

    def test_assibilation_fires_on_brad(self, inventory, pipeline) -> None:
        ur = build_ur("brad", "i", inventory, stem_final="d")
        exp = explain_derivation(pipeline, ur)
        assib = next(t for t in exp.rule_traces if t.rule_name == "assibilation")
        assert assib.fired
        assert assib.search_outcome is SearchOutcome.LICENSED

    def test_dorsal_pal_blocks_in_punct(self, inventory, pipeline) -> None:
        ur = build_ur("punkt", "e", inventory, stem_final="t")
        exp = explain_derivation(pipeline, ur)
        dp = next(t for t in exp.rule_traces if t.rule_name == "dorsal-pal")
        # /K/ blocked by broad terminator halting on /T/ (fails +Syll,+Front).
        assert not dp.fired
        # The trace should record either a blocked search or a unify conflict
        # for the /K/ target.
        assert (
            dp.search_outcome is SearchOutcome.BROAD_TERMINATOR_BLOCK
            or dp.unify_conflict is not None
        )

    def test_format_explanation_returns_string(
        self, inventory, pipeline,
    ) -> None:
        ur = build_ur("prost", "i", inventory, stem_final="s")
        text = format_explanation(explain_derivation(pipeline, ur))
        assert "s-pal-rev" in text
        assert "bleed-rev" in text


# ---------------------------------------------------------------------------
# ordering
# ---------------------------------------------------------------------------


class TestOrdering:
    def test_declared_ordering_succeeds_on_working_case(
        self, inventory, pipeline,
    ) -> None:
        ur = build_ur("muskə", "e", inventory, stem_final="s")
        result = find_valid_orderings(
            rules=RULES_FROM_PAPER,
            ur=ur,
            expected_field="muʃte",
            resolve=lambda l: inventory.segment(l),
            inventory=inventory,
            strategy=OrderingStrategy.DECLARED,
        )
        assert result.baseline_matched
        assert len(result.successful) >= 1

    def test_adjacent_swap_stays_tractable(self, inventory) -> None:
        ur = build_ur("rus", "i", inventory, stem_final="s")
        result = find_valid_orderings(
            rules=RULES_FROM_PAPER,
            ur=ur,
            expected_field="ruʃi",
            resolve=lambda l: inventory.segment(l),
            inventory=inventory,
            strategy=OrderingStrategy.ADJACENT_SWAP,
        )
        # N-1 adjacent swaps + baseline = 9 permutations
        assert result.orderings_tried <= len(RULES_FROM_PAPER)


# ---------------------------------------------------------------------------
# perturbation
# ---------------------------------------------------------------------------


class TestPerturbation:
    def test_baseline_working_case_reports_success(
        self, inventory,
    ) -> None:
        # A case that already works should report baseline matched.
        def builder(inv):
            return build_ur("rus", "i", inv, stem_final="s")
        report = search_perturbations(
            ur_builder=builder,
            expected_field="ruʃi",
            inventory_json=INVENTORY_JSON,
            baseline_patches=PATCHES,
            baseline_underspec=UNDERSPEC,
            baseline_rules=RULES_FROM_PAPER,
            kinds=(PerturbationKind.PATCH_VALUE,),
            limit=20,
        )
        assert report.baseline_matched
        assert report.tried > 0


# ---------------------------------------------------------------------------
# categorise
# ---------------------------------------------------------------------------


class _FakeRow:
    def __init__(self, lemma, plural, stem_final, opp, mutation):
        self.lemma = lemma
        self.plural = plural
        self.stem_final = stem_final
        self.opportunity = opp
        self.mutation = mutation


class TestCategorise:
    def test_trailing_glide_is_data_side(self) -> None:
        row = _FakeRow("prost", "proști", "s", "i", True)
        cat = refine_category(row, "proʃtj", "proʃti", True, ur_built=True)
        assert cat.category is RefinedCategory.NORM_TRAILING_GLIDE

    def test_ez_ethnonym_survives(self) -> None:
        row = _FakeRow("albanez", "albanezi", "z", "i", False)
        cat = refine_category(
            row, "albaneʒj", "albanezi", True, ur_built=True,
        )
        # trailing j vs i doesn't match under normalization because
        # the medial ʒ vs z differs — should NOT be NORM_TRAILING_GLIDE.
        # Should land on EZ_ETHNONYM.
        assert cat.category is RefinedCategory.EZ_ETHNONYM

    def test_ur_build_failed(self) -> None:
        row = _FakeRow("foo", "foos", "s", "i", False)
        cat = refine_category(row, "", "foos", False, ur_built=False)
        assert cat.category is RefinedCategory.UR_BUILD_FAILED


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
