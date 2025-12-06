#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test suite for Romanian palatalization pipeline invariants.

These tests validate the conceptual correctness of the data pipeline:
- TP domain membership (what counts as N vs exceptions)
- NDE classification (gimpe/ochi/paduchi)
- Mutation vs opportunity consistency
- Exclusion of NDEB from main TP exception counts

Run with: pytest tests/test_pipeline_invariants.py -v
"""

import re
from pathlib import Path

import pandas as pd
import pytest

# Constants matching the R analysis script
SEGMENTS_OF_INTEREST = {"c", "g", "t", "d", "s", "z"}
PLURAL_OPPORTUNITIES = {"i", "e"}
NDE_CLASSES = {"gimpe", "ochi", "paduchi"}
NDE_OBSERVABLE = {"ochi", "paduchi"}  # Observable as DE in plural domain

# Data paths
PROJECT_ROOT = Path(__file__).parent.parent
LEX_PATH = PROJECT_ROOT / "data" / "romanian_lexicon_with_freq.csv"


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture(scope="module")
def lex():
    """Load the processed lexicon once for all tests."""
    assert LEX_PATH.exists(), f"Lexicon not found at {LEX_PATH}"
    df = pd.read_csv(LEX_PATH)

    # Derive nde_class from exception_reason for tests
    # exception_reason values: "nde:gimpe", "nde:ochi", "nde:paduchi", "undergoer", etc.
    df["nde_class"] = df["exception_reason"].apply(
        lambda x: x.replace("nde:", "") if isinstance(x, str) and x.startswith("nde:") else ""
    )

    # Create opportunity_tp from opportunity for compatibility with tests
    df["opportunity_tp"] = df["opportunity"]

    return df


@pytest.fixture(scope="module")
def nouns(lex):
    """All noun rows."""
    return lex[lex["pos"].str.upper() == "N"].copy()


@pytest.fixture(scope="module")
def nouns_opp(nouns):
    """
    Full i/e domain (including NDEB).
    Matches R script's nouns_opp definition.
    """
    return nouns[
        nouns["opportunity_tp"].isin(PLURAL_OPPORTUNITIES)
        & nouns["stem_final"].isin(SEGMENTS_OF_INTEREST)
    ].copy()


@pytest.fixture(scope="module")
def nouns_tp(nouns):
    """
    PRODUCTIVE TP domain (excluding NDEB and structural non-alternators).
    Matches R script's nouns_tp definition using tp_in_domain flag.
    """
    # tp_in_domain is stored as boolean True/False in CSV
    return nouns[nouns["tp_in_domain"] == True].copy()  # noqa: E712


@pytest.fixture(scope="module")
def nouns_opp_no_ndeb(nouns_opp):
    """
    Legacy: i/e domain excluding NDEB (manual filtering).
    Should match nouns_tp closely.
    """
    return nouns_opp[~nouns_opp["nde_class"].isin(NDE_CLASSES)].copy()


# =============================================================================
# Helper Functions
# =============================================================================


def has_front_env(lemma: str) -> bool:
    """
    Check if lemma contains C+{i,e} tautomorphemic sequence OR palatal graphemes.

    GIMPE = tautomorphemic C+front-vowel, which includes:
    - Dorsals: [cg] before i/e (excluding hC sequences)
    - Coronals: [sztd] before i/e
    - Palatal graphemes already present: [țșjğ]

    Examples:
        aciu → True (ci)
        abagiu → True (gi)
        abazie → True (zie)
        aboliționistă → True (has ș)
        batuc → False (no C+{i,e} or palatals)
    """
    if not isinstance(lemma, str):
        return False

    # Dorsals: non-h [cg] before i/e
    if re.search(r"(?<!h)[cg][ie]", lemma):
        return True

    # Coronals: [sztd] before i/e
    if re.search(r"[sztd][ie]", lemma):
        return True

    # Palatal graphemes (outputs of palatalization)
    if re.search(r"[țșjğ]", lemma):
        return True

    return False


# =============================================================================
# Structural Tests: TP Domain Definition
# =============================================================================


def test_tp_domain_exists(lex):
    """tp_in_domain field exists and has correct type."""
    assert "tp_in_domain" in lex.columns
    assert lex["tp_in_domain"].dtype == bool


def test_tp_domain_nonempty(nouns_tp):
    """Productive TP domain is non-empty."""
    assert len(nouns_tp) > 0, "TP domain is empty!"


def test_tp_domain_matches_definition(nouns_tp):
    """All tp_in_domain=True items meet basic criteria."""
    # All should be nouns
    assert (nouns_tp["pos"].str.upper() == "N").all()

    # All should have i/e opportunity
    assert nouns_tp["opportunity"].isin(PLURAL_OPPORTUNITIES).all()

    # All should have target segments
    assert nouns_tp["stem_final"].isin(SEGMENTS_OF_INTEREST).all()

    # No NDE classes
    assert not nouns_tp["nde_class"].isin(NDE_CLASSES).any()

    # No "non exception" items
    assert not (nouns_tp["exception_reason"] == "non exception").any()


def test_tp_domain_excludes_ndeb(nouns_tp):
    """TP domain should contain NO NDEB items."""
    assert not nouns_tp["nde_class"].isin(NDE_CLASSES).any()


def test_tp_domain_vs_manual_filtering(nouns_tp, nouns_opp_no_ndeb):
    """
    tp_in_domain should closely match manual NDEB filtering.

    Difference should be only structural non-alternators
    (frontstem dorsals, underlying palatals, suffix-internal).
    """
    # nouns_opp_no_ndeb is the "old" way (just exclude NDE classes)
    # nouns_tp is the "new" way (exclude NDE + underlying palatals, etc.)

    diff = len(nouns_opp_no_ndeb) - len(nouns_tp)

    # We expect ~1200-1300 items removed (frontstem dorsals, etc.)
    assert 1000 <= diff <= 1500, (
        f"Unexpected difference between manual and tp_in_domain filtering: "
        f"{diff} items"
    )


def test_no_ndeb_in_tp_exceptions(nouns_tp):
    """
    Main TP exception counts should contain NO NDEB items.

    This is the critical test: NDEB should not inflate exception counts.
    """
    exceptions = nouns_tp[nouns_tp["mutation"] == False]  # noqa: E712

    assert not exceptions["nde_class"].isin(NDE_CLASSES).any(), (
        "Found NDEB items in TP domain exceptions! "
        "These should be excluded by tp_in_domain filter."
    )


# =============================================================================
# NDE Classification Tests
# =============================================================================


def test_all_gimpe_have_front_env(lex):
    """
    Most gimpe-type items should have C+{i,e} or palatals in lemma.

    GIMPE = tautomorphemic C+front-vowel sequences OR vowel-only changes.
    Some GIMPE items are vowel-only (ă→e) without C+front sequences.
    We just check that the majority have the expected pattern.
    """
    gimpe = lex[lex["nde_class"] == "gimpe"]

    if gimpe.empty:
        pytest.skip("No gimpe items in lexicon")

    with_pattern = gimpe[gimpe["lemma"].apply(has_front_env)]

    # At least 85% should have the pattern
    # (remaining are vowel-only changes like acuarelistă→acuareliste)
    proportion = len(with_pattern) / len(gimpe)
    assert proportion >= 0.85, (
        f"Only {proportion:.1%} of GIMPE items have C+{{i,e}} or palatals in lemma "
        f"(expected >= 85%)"
    )


def test_gimpe_items_have_correct_exception_reason(lex):
    """All gimpe items should have exception_reason = 'nde:gimpe'."""
    gimpe = lex[lex["nde_class"] == "gimpe"]

    if gimpe.empty:
        pytest.skip("No gimpe items in lexicon")

    assert (gimpe["exception_reason"] == "nde:gimpe").all()


def test_ochi_items_have_lemma_equals_plural(lex):
    """
    OCHI-type items should have lemma = plural.

    These are ambiguous: is the final -i part of root or plural suffix?
    """
    ochi = lex[lex["nde_class"] == "ochi"]

    if ochi.empty:
        pytest.skip("No ochi items in lexicon")

    # All ochi items should have lemma = plural
    assert (ochi["lemma"] == ochi["plural"]).all(), (
        "Found OCHI items where lemma != plural"
    )


def test_paduchi_dorsal_items_have_h_blocking(lex):
    """
    PADUCHI dorsal items should have che/ghe → chi/ghi pattern.

    This is h-insertion blocking palatalization.
    """
    paduchi_dorsal = lex[
        (lex["nde_class"] == "paduchi") & (lex["stem_final"].isin({"c", "g"}))
    ]

    if paduchi_dorsal.empty:
        pytest.skip("No paduchi dorsal items in lexicon")

    # Check a few have che→chi or ghe→ghi pattern
    # (Not all will, some might be chi→chi OCHI-style)
    che_chi = paduchi_dorsal[
        paduchi_dorsal["lemma"].str.endswith(("che", "ghe"))
        & paduchi_dorsal["plural"].str.endswith(("chi", "ghi"))
    ]

    assert len(che_chi) > 0, (
        "Expected at least some che→chi or ghe→ghi PADUCHI items"
    )


# =============================================================================
# Specific Lemma Tests (Golden Examples)
# =============================================================================


def expect_gimpe(lex, lemma: str):
    """Helper: expect a specific lemma to be classified as gimpe."""
    rows = lex[lex["lemma"] == lemma]
    assert not rows.empty, f"Missing lemma {lemma!r} in lexicon"
    assert (
        rows["nde_class"] == "gimpe"
    ).any(), f"{lemma} is not labeled gimpe"

    # If in i/e domain, mutation should be False
    in_ie = rows["opportunity_tp"].isin(PLURAL_OPPORTUNITIES)
    if in_ie.any():
        assert not rows.loc[in_ie, "mutation"].any(), (
            f"{lemma} is gimpe but mutation=True"
        )


def expect_paduchi(lex, lemma: str):
    """Helper: expect a specific lemma to be classified as paduchi."""
    rows = lex[lex["lemma"] == lemma]
    assert not rows.empty, f"Missing lemma {lemma!r} in lexicon"
    assert (
        rows["nde_class"] == "paduchi"
    ).any(), f"{lemma} is not labeled paduchi"
    assert (
        rows["exception_reason"] == "nde:paduchi"
    ).any()


def expect_undergoes(lex, lemma: str):
    """Helper: expect a specific lemma to palatalize."""
    rows = lex[lex["lemma"] == lemma]
    assert not rows.empty, f"Missing lemma {lemma!r} in lexicon"
    assert rows["mutation"].any(), f"{lemma} does not palatalize (mutation=False)"
    assert (
        rows["exception_reason"] == "undergoer"
    ).any()


@pytest.mark.parametrize(
    "lemma",
    [
        "abagiu",  # gi in lemma → gii in plural
        "aciu",  # ci in lemma → cii in plural
        "acacia",  # ci in lemma → cii in plural
        "abazie",  # zie in lemma → zii in plural (vowel-only)
    ],
)
def test_specific_gimpe_examples(lex, lemma):
    """Specific lemmas that should be gimpe-type."""
    expect_gimpe(lex, lemma)


@pytest.mark.parametrize(
    "lemma",
    [
        "genunche",  # che → chi (h-blocking)
        "genuche",  # che → chi (h-blocking)
    ],
)
def test_specific_paduchi_examples(lex, lemma):
    """Specific lemmas that should be paduchi-type."""
    expect_paduchi(lex, lemma)


@pytest.mark.parametrize(
    "lemma,plural",
    [
        ("ciurechi", "ciurechi"),  # lemma = plural, chi pattern
        ("genunchi", "genunchi"),  # lemma = plural, chi pattern
    ],
)
def test_specific_ochi_examples(lex, lemma, plural):
    """Specific lemmas that should be ochi-type."""
    rows = lex[lex["lemma"] == lemma]
    assert not rows.empty, f"Missing lemma {lemma!r}"
    assert (rows["nde_class"] == "ochi").any()
    assert (rows["lemma"] == rows["plural"]).all()


@pytest.mark.parametrize(
    "lemma",
    [
        "abulic",  # c+i → č+i (regular palatalization: abulic → abulici)
        "acrobatic",  # c+i → č+i (regular palatalization)
    ],
)
def test_specific_undergoes_examples(lex, lemma):
    """Specific lemmas that should palatalize regularly."""
    expect_undergoes(lex, lemma)


# =============================================================================
# Mutation vs Opportunity Consistency
# =============================================================================


def test_mutation_only_in_ie_domain(nouns):
    """
    Items with mutation=True should be in i/e opportunity domain.

    This catches cases where mutation is set incorrectly.
    """
    mutated = nouns[nouns["mutation"] == True]  # noqa: E712

    bad = mutated[~mutated["opportunity"].isin(PLURAL_OPPORTUNITIES)]

    assert bad.empty, (
        f"Found {len(bad)} mutated items outside i/e opportunity:\n"
        f"{bad[['lemma', 'plural', 'opportunity', 'stem_final']].head(10)}"
    )


def test_no_na_stem_final_in_tp_domain(nouns_tp):
    """TP domain should have no NA stem_final values."""
    assert not nouns_tp["stem_final"].isna().any()


# =============================================================================
# Dorsal vs Coronal Behavior
# =============================================================================


def test_dorsals_in_tp_domain_all_palatalize(nouns_tp):
    """
    All dorsals (c/g) in TP domain should palatalize.

    This is a KEY test: after excluding NDEB and frontstem dorsals,
    all remaining dorsals should be 100% undergoers.
    """
    dorsals = nouns_tp[nouns_tp["stem_final"].isin({"c", "g"})]

    if dorsals.empty:
        pytest.skip("No dorsals in TP domain")

    # All should palatalize
    non_palatalizing = dorsals[dorsals["mutation"] == False]  # noqa: E712

    assert non_palatalizing.empty, (
        f"Found {len(non_palatalizing)} non-palatalizing dorsals in TP domain!\n"
        f"These should have been excluded as NDEB or frontstem dorsals:\n"
        f"{non_palatalizing[['lemma', 'plural', 'stem_final', 'nde_class', 'exception_reason']].head(10)}"
    )


def test_coronals_show_variation(nouns_tp):
    """
    Coronals (t/d/s/z) in TP domain should show variation.

    Unlike dorsals, coronals don't all palatalize.
    We just check that we have both palatalizing and non-palatalizing cases.
    """
    coronals = nouns_tp[nouns_tp["stem_final"].isin({"t", "d", "s", "z"})]

    if coronals.empty:
        pytest.skip("No coronals in TP domain")

    # Should have both mutation=True and mutation=False
    palatalizing = coronals[coronals["mutation"] == True]  # noqa: E712
    non_palatalizing = coronals[coronals["mutation"] == False]  # noqa: E712

    assert len(palatalizing) > 0, "No palatalizing coronals found"
    assert len(non_palatalizing) > 0, "No non-palatalizing coronals found"


# =============================================================================
# Count Stability Tests
# =============================================================================


def test_total_ndeb_count(lex):
    """Total NDEB count should be stable across runs."""
    ndeb = lex[lex["nde_class"].isin(NDE_CLASSES)]

    # Expected: ~2200-2300 NDEB items
    assert 2200 <= len(ndeb) <= 2400, (
        f"Unexpected NDEB count: {len(ndeb)} "
        f"(expected ~2200-2300)"
    )


def test_ndeb_breakdown(lex):
    """NDEB breakdown by class should be stable."""
    ndeb = lex[lex["nde_class"].isin(NDE_CLASSES)]

    breakdown = ndeb["nde_class"].value_counts().to_dict()

    # Expected approximate counts (allow some variation)
    assert 2000 <= breakdown.get("gimpe", 0) <= 2400, (
        f"Unexpected gimpe count: {breakdown.get('gimpe', 0)}"
    )
    assert 10 <= breakdown.get("ochi", 0) <= 20, (
        f"Unexpected ochi count: {breakdown.get('ochi', 0)}"
    )
    assert 50 <= breakdown.get("paduchi", 0) <= 80, (
        f"Unexpected paduchi count: {breakdown.get('paduchi', 0)}"
    )


def test_tp_domain_size(nouns_tp):
    """TP domain size should be stable (~4800-5000 items)."""
    assert 4500 <= len(nouns_tp) <= 5200, (
        f"Unexpected TP domain size: {len(nouns_tp)} "
        f"(expected ~4800-5000)"
    )


def test_tp_domain_exception_count(nouns_tp):
    """Exception count in TP domain should be stable (~2600-2800)."""
    exceptions = nouns_tp[nouns_tp["mutation"] == False]  # noqa: E712

    assert 2400 <= len(exceptions) <= 3000, (
        f"Unexpected exception count: {len(exceptions)} "
        f"(expected ~2600-2800)"
    )


# =============================================================================
# Run Tests
# =============================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
