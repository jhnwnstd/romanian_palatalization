#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Romanian Palatalization Data Processor - Main Script

Complete pipeline for deriving palatalization-related fields from
Romanian Wiktionary data.

Usage:
    python romanian_processor_main.py

Input:
    romanian_lexicon_raw_dex.csv

Output:
    romanian_lexicon_complete.csv
"""

import csv
import sys
from pathlib import Path
from typing import Dict, Optional

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

# Import all derivation functions from library
# Part 1: Basic derivations
# Part 2: Alignment-based derivations
# Part 3: Suffix and G2P
# Part 4: NDE and exceptions
# pylint: disable=wrong-import-position
from romanian_processor_lib import (  # noqa: E402
    derive_derived_adj_fields,
    derive_derived_verbs_fields,
    derive_exception_reason,
    derive_is_true_exception,
    derive_lemma_suffix,
    derive_mutation_and_orth_change,
    derive_nde_class,
    derive_opportunity,
    derive_palatal_consonant_pl,
    derive_stem_final_and_cluster,
    derive_suffix_triggers_plural_mutation,
    derive_target_is_suffix,
    explode_derived_verbs_row,
    set_ipa_normalizer,
    to_ipa,
    tweak_nominal_ipa,
    validate_plural_quality,
)

# Import normalization functions
from wiktionary_normalizer import (  # noqa: E402
    normalize_ipa,
    normalize_orthography,
)

# pylint: enable=wrong-import-position

# Configure the library to use IPA normalizer
# This keeps dependencies linear: main → lib and main → normalizer (not lib → normalizer)
set_ipa_normalizer(lambda ipa: normalize_ipa(ipa, remove_stress=True))


def process_row(row: Dict[str, str]) -> Optional[Dict[str, str]]:
    """
    Process a single CSV row through the complete pipeline.

    Args:
        row: Dictionary representing a CSV row

    Returns:
        Processed row with all derived fields added, or None if row
        should be skipped
    """
    result = row.copy()
    for field in ("lemma", "plural", "gloss", "etym_lang"):
        value = result.get(field)
        if value:
            result[field] = normalize_orthography(value)

    def _ensure_ipa(
        orth_key: str,
        raw_key: str,
        norm_key: str,
        tweak_fn=None,
    ) -> None:
        raw_val = result.get(raw_key)
        if raw_val:
            result[norm_key] = normalize_ipa(raw_val, remove_stress=True)
            return
        orth_val = result.get(orth_key, "")
        if orth_val:
            ipa = to_ipa(orth_val)
            if tweak_fn is not None:
                ipa = tweak_fn(orth_val, ipa)
            result[norm_key] = normalize_ipa(ipa, remove_stress=True)

    _ensure_ipa(
        orth_key="lemma",
        raw_key="ipa_raw_lemma",
        norm_key="ipa_normalized_lemma",
        tweak_fn=tweak_nominal_ipa,
    )
    _ensure_ipa(
        orth_key="plural",
        raw_key="ipa_raw_pl",
        norm_key="ipa_normalized_pl",
        tweak_fn=None,
    )
    result["pos"] = (result.get("pos") or "").strip().upper()
    if result["pos"] == "N" and not result.get("gender"):
        return None
    derive_stem_final_and_cluster(result)
    validate_plural_quality(result)
    derive_mutation_and_orth_change(result)
    derive_opportunity(result)
    derive_palatal_consonant_pl(result)
    derive_lemma_suffix(result)
    derive_target_is_suffix(result)
    derive_suffix_triggers_plural_mutation(result)
    derive_derived_verbs_fields(result)
    derive_derived_adj_fields(result)
    derive_nde_class(result)
    derive_exception_reason(result)
    derive_is_true_exception(result)
    return result


OUTPUT_FIELDS = [
    "lemma",
    "gloss",
    "pos",
    "gender",
    "stem_final",
    "cluster",
    "plural",
    "mutation",
    "orth_change",
    "opportunity",
    "palatal_consonant_pl",
    "ipa_normalized_lemma",
    "ipa_normalized_pl",
    "derived_verbs",
    "deriv_suffixes",
    "ipa_derived_verbs",
    "derived_adj",
    "ipa_derived_adj",
    "etym_lang",
    "exception_reason",
    "nde_class",
    "is_true_exception",
    "lemma_suffix",
    "target_is_suffix",
    "suffix_triggers_plural_mutation",
    "source",
    "notes",
    "ipa_raw_lemma",
    "ipa_raw_pl",
]


def process_csv(input_path: str, output_path: str) -> None:
    """
    Read input CSV, process all rows, write output CSV with only OUTPUT_FIELDS.
    """
    print(f"Reading {input_path}...")

    with open(input_path, "r", encoding="utf-8") as f_in:
        reader = csv.DictReader(f_in)

        rows = []
        for i, row in enumerate(reader, start=1):
            if i % 1000 == 0:
                print(f"  Processed {i} rows...")
            processed = process_row(row)
            if processed is None:
                continue
            exploded = explode_derived_verbs_row(processed)
            rows.extend(exploded)

    print(f"\nWriting {len(rows)} rows to {output_path}...")

    with open(output_path, "w", encoding="utf-8", newline="") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=OUTPUT_FIELDS)
        writer.writeheader()

        for r in rows:
            filtered_row = {field: r.get(field, "") for field in OUTPUT_FIELDS}
            writer.writerow(filtered_row)

    print("Done!")
    print("\n" + "=" * 80)
    print("BASIC SUMMARY")
    print("=" * 80)

    total = len(rows)
    nouns = sum(1 for r in rows if r.get("pos") == "N")

    print(f"\nTotal entries: {total}")
    print(f"Nouns: {nouns}")

    mutations_true = sum(1 for r in rows if r.get("mutation") == "True")
    mutations_false = sum(1 for r in rows if r.get("mutation") == "False")

    print("\nMutation:")
    print(f"  True: {mutations_true}")
    print(f"  False: {mutations_false}")

    opp_counts: Dict[str, int] = {}
    for r in rows:
        opp = r.get("opportunity", "none")
        opp_counts[opp] = opp_counts.get(opp, 0) + 1

    print("\nOpportunity distribution:")
    for opp in ["i", "e", "uri", "none"]:
        print(f"  {opp}: {opp_counts.get(opp, 0)}")

    nde_gimpe = sum(1 for r in rows if r.get("nde_class") == "gimpe")
    nde_ochi = sum(1 for r in rows if r.get("nde_class") == "ochi")
    nde_paduchi = sum(1 for r in rows if r.get("nde_class") == "paduchi")
    true_exceptions = sum(
        1 for r in rows if r.get("is_true_exception") == "True"
    )

    print("\nException classification:")
    print(f"  NDE gimpe (tautomorphemic): {nde_gimpe}")
    print(f"  NDE ochi (singular=plural): {nde_ochi}")
    print(f"  NDE paduchi (che/ghe→chi/ghi): {nde_paduchi}")
    print(f"  True exceptions (unexplained): {true_exceptions}")

    print("\n" + "=" * 80)


if __name__ == "__main__":
    base_dir = Path(__file__).parent.parent
    INPUT_CSV = base_dir / "data" / "romanian_lexicon_raw_dex.csv"
    OUTPUT_CSV = base_dir / "data" / "romanian_lexicon_complete.csv"
    process_csv(str(INPUT_CSV), str(OUTPUT_CSV))
