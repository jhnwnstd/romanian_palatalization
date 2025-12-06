#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Snapshot tests for TP summary tables.

These tests validate that key TP numbers remain stable across pipeline changes.
If these fail, it means the pipeline logic has changed in a way that affects
the final TP calculations - which may or may not be intentional.

Run with: pytest tests/test_tp_snapshot.py -v
"""

from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).parent.parent
TP_FULL_PATH = PROJECT_ROOT / "analysis" / "romanian_tp_summary_full.csv"
TP_DS_PATH = (
    PROJECT_ROOT / "analysis" / "romanian_tp_summary_downsampled.csv"
)


@pytest.fixture(scope="module")
def tp_full():
    """Load full-lexicon TP summary."""
    assert TP_FULL_PATH.exists(), f"TP summary not found at {TP_FULL_PATH}"
    return pd.read_csv(TP_FULL_PATH)


@pytest.fixture(scope="module")
def tp_downsampled():
    """Load downsampled TP summary."""
    assert TP_DS_PATH.exists(), f"TP summary not found at {TP_DS_PATH}"
    return pd.read_csv(TP_DS_PATH)


def get_row(df, type_label):
    """Helper to get a single row by type label."""
    rows = df[df["type"] == type_label]
    assert len(rows) == 1, f"Expected 1 row for '{type_label}', got {len(rows)}"
    return rows.iloc[0]


# =============================================================================
# Full Lexicon Snapshot Tests
# =============================================================================


def test_full_gimpe_type(tp_full):
    """GIMPE type counts (full lexicon)."""
    row = get_row(tp_full, "gimpe type")

    # These should be very stable
    assert 2150 <= row["N"] <= 2250, f"Unexpected gimpe N: {row['N']}"
    assert row["mutated"] == 0, "GIMPE should not mutate"
    assert row["non-mutated"] == row["N"]
    assert row["rate"] == 0.0
    assert row["majority?"] is False
    assert row["tolerated?"] is True


def test_full_ochi_type(tp_full):
    """OCHI type counts (full lexicon)."""
    row = get_row(tp_full, "ochi-ochi type")

    # Small stable count
    assert 10 <= row["N"] <= 20, f"Unexpected ochi N: {row['N']}"
    assert row["mutated"] == 0
    assert row["rate"] == 0.0


def test_full_paduchi_type(tp_full):
    """PADUCHI type counts (full lexicon)."""
    row = get_row(tp_full, "paduche-paduchi type")

    # Small stable count
    assert 50 <= row["N"] <= 80, f"Unexpected paduchi N: {row['N']}"
    assert row["mutated"] == 0
    assert row["rate"] == 0.0


def test_full_dorsals_c_i_e(tp_full):
    """Dorsal c + i,e combined (full lexicon)."""
    row = get_row(tp_full, "<c> + <-i, -e> plural")

    # With NDEB included, we expect ~1400-1500
    assert 1300 <= row["N"] <= 1600, f"Unexpected c i+e N: {row['N']}"

    # Should have majority mutating
    assert row["majority?"] is True


def test_full_dorsals_g_i_e(tp_full):
    """Dorsal g + i,e combined (full lexicon)."""
    row = get_row(tp_full, "<g> + <-i, -e> plural")

    # With NDEB, expect ~800-900
    assert 750 <= row["N"] <= 1000, f"Unexpected g i+e N: {row['N']}"


def test_full_coronals_t_i(tp_full):
    """Coronal t + i (full lexicon)."""
    row = get_row(tp_full, "<t> + <-i> plural")

    # Should have good N
    assert 1100 <= row["N"] <= 1400, f"Unexpected t+i N: {row['N']}"

    # Should have high mutation rate
    assert row["rate"] > 0.75, f"Unexpected t+i rate: {row['rate']}"
    assert row["majority?"] is True


def test_full_coronals_z_e(tp_full):
    """Coronal z + e (full lexicon) - almost never palatalizes."""
    row = get_row(tp_full, "<z> + <-e> plural")

    # z+e should have very low mutation rate
    assert row["rate"] < 0.05, f"Unexpected z+e rate: {row['rate']}"
    assert row["majority?"] is False


def test_full_cluster_st_i(tp_full):
    """Cluster st + i (full lexicon) - should palatalize well."""
    row = get_row(tp_full, "<st> + <-i> plural")

    # st + i should have very high mutation rate
    assert row["N"] > 500, f"Unexpected st+i N: {row['N']}"
    assert row["rate"] > 0.95, f"Unexpected st+i rate: {row['rate']}"
    assert row["tolerated?"] is True


# =============================================================================
# Downsampled Lexicon Snapshot Tests
# =============================================================================


def test_downsampled_gimpe_stable(tp_downsampled):
    """GIMPE counts in downsampled lexicon."""
    row = get_row(tp_downsampled, "gimpe type")

    # Downsampled should have ~180-200 gimpe
    assert 150 <= row["N"] <= 220, f"Unexpected downsampled gimpe N: {row['N']}"
    assert row["mutated"] == 0
    assert row["tolerated?"] is True


def test_downsampled_c_i_stable(tp_downsampled):
    """Dorsal c + i in downsampled lexicon."""
    row = get_row(tp_downsampled, "<c> + <-i> plural (downsampled)")

    # Should have ~80-120 items
    assert 70 <= row["N"] <= 130, f"Unexpected downsampled c+i N: {row['N']}"

    # Should have high mutation rate
    assert row["rate"] > 0.70, f"Unexpected downsampled c+i rate: {row['rate']}"
    assert row["majority?"] is True


def test_downsampled_g_i_stable(tp_downsampled):
    """Dorsal g + i in downsampled lexicon."""
    row = get_row(tp_downsampled, "<g> + <-i> plural (downsampled)")

    # Smaller count
    assert 30 <= row["N"] <= 60, f"Unexpected downsampled g+i N: {row['N']}"


def test_downsampled_t_i_stable(tp_downsampled):
    """Coronal t + i in downsampled lexicon."""
    row = get_row(tp_downsampled, "<t> + <-i> plural (downsampled)")

    # Should have reasonable N
    assert 150 <= row["N"] <= 250, f"Unexpected downsampled t+i N: {row['N']}"

    # Should have high mutation rate
    assert row["rate"] > 0.75, f"Unexpected downsampled t+i rate: {row['rate']}"
    assert row["majority?"] is True


def test_downsampled_st_i_stable(tp_downsampled):
    """Cluster st + i in downsampled lexicon."""
    row = get_row(tp_downsampled, "<st> + <-i> plural (downsampled)")

    # Should have some items (downsampling reduces count)
    assert row["N"] > 30, f"Unexpected downsampled st+i N: {row['N']}"

    # Should have very high rate
    assert row["rate"] > 0.90, f"Unexpected downsampled st+i rate: {row['rate']}"


# =============================================================================
# Derivational Snapshot Tests
# =============================================================================


def test_full_nv_derivation_all(tp_full):
    """N→V derivation ALL (full lexicon)."""
    row = get_row(tp_full, "N→V derivation, ALL")

    # Should have some items
    assert row["N"] > 30, f"Unexpected N→V ALL N: {row['N']}"

    # Should show some palatalization
    assert row["mutated"] > 0, "Expected some N→V palatalization"


def test_full_nadj_derivation_all(tp_full):
    """N→Adj derivation ALL (full lexicon)."""
    row = get_row(tp_full, "N→Adj derivation, ALL")

    # Should have some items
    assert row["N"] > 10, f"Unexpected N→Adj ALL N: {row['N']}"


# =============================================================================
# Consistency Tests Across Tables
# =============================================================================


def test_combined_equals_sum_of_parts(tp_full):
    """
    Combined i+e counts should equal sum of i and e counts.

    This validates internal consistency of the TP tables.
    """
    # Test for c
    c_i = get_row(tp_full, "<c> + <-i> plural")
    c_e = get_row(tp_full, "<c> + <-e> plural")
    c_ie = get_row(tp_full, "<c> + <-i, -e> plural")

    assert c_ie["N"] == c_i["N"] + c_e["N"], (
        f"c i+e N ({c_ie['N']}) != c i N ({c_i['N']}) + c e N ({c_e['N']})"
    )
    assert c_ie["mutated"] == c_i["mutated"] + c_e["mutated"]
    assert c_ie["non-mutated"] == c_i["non-mutated"] + c_e["non-mutated"]

    # Test for g
    g_i = get_row(tp_full, "<g> + <-i> plural")
    g_e = get_row(tp_full, "<g> + <-e> plural")
    g_ie = get_row(tp_full, "<g> + <-i, -e> plural")

    assert g_ie["N"] == g_i["N"] + g_e["N"]
    assert g_ie["mutated"] == g_i["mutated"] + g_e["mutated"]


def test_full_vs_downsampled_rates_similar(tp_full, tp_downsampled):
    """
    Mutation rates in downsampled lexicon should be similar to full lexicon.

    This checks that downsampling preserves the grammar structure.
    """
    # Compare c + i rates
    c_i_full = get_row(tp_full, "<c> + <-i> plural")
    c_i_ds = get_row(tp_downsampled, "<c> + <-i> plural (downsampled)")

    # Rates should be within ~10% of each other
    assert abs(c_i_full["rate"] - c_i_ds["rate"]) < 0.15, (
        f"c+i rates differ too much: full={c_i_full['rate']:.3f}, "
        f"ds={c_i_ds['rate']:.3f}"
    )

    # Compare t + i rates
    t_i_full = get_row(tp_full, "<t> + <-i> plural")
    t_i_ds = get_row(tp_downsampled, "<t> + <-i> plural (downsampled)")

    assert abs(t_i_full["rate"] - t_i_ds["rate"]) < 0.15


# =============================================================================
# Run Tests
# =============================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
