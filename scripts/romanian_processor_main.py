#!/usr/bin/env python3
"""Romanian palatalization processor — pipeline stage 4.

Reads ``romanian_lexicon_raw_dex.csv`` (the DEX-QC'd raw lexicon) and
writes ``romanian_lexicon_complete.csv`` with all derived columns:
``stem_final``, ``mutation``, ``opportunity``, ``derived_verbs``,
``deriv_suffixes``, IPA, etc.

The actual derivation logic lives in ``src/romanian_processor_lib.py``;
this script orchestrates the row-level pass, dedup, and writing.

Invariants
----------
- Each ``(lemma, pos, plural)`` triple appears at most once in the
  output. Duplicates from Wiktionary merge by preferring the row with
  more derivation data.
- A row whose ``lemma`` equals the plural of a different row with the
  same pos is dropped (it's a stray plural form mis-listed as a lemma)
  — except for the ochi/genunchi case where lemma == plural for the
  same entry.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

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
    derive_frontstem_dorsal,
    derive_lemma_suffix,
    derive_mutation_and_orth_change,
    derive_nde_class,
    derive_opportunity,
    derive_palatal_consonant_pl,
    derive_stem_final_and_cluster,
    derive_target_is_suffix,
    ensure_ipa_fields,
    explode_derived_verbs_row,
    fix_nde_mutations,
    fix_underlying_palatals,
    mark_tp_domain,
    set_ipa_normalizer,
    should_process_row,
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
# This keeps dependencies linear: main → lib and main → normalizer
# (not lib → normalizer)
set_ipa_normalizer(lambda ipa: normalize_ipa(ipa, remove_stress=True))


def process_row(row: dict[str, str]) -> dict[str, str] | None:
    """
    Process a single CSV row through the complete pipeline.

    Args:
        row: Dictionary representing a CSV row

    Returns:
        Processed row with all derived fields added, or None if row
        should be skipped
    """
    # Early filter: skip rows without target consonants
    if not should_process_row(row):
        return None

    result = row.copy()

    # Normalize orthographic fields
    for field in ("lemma", "plural", "gloss", "etym_lang"):
        value = result.get(field)
        if value:
            result[field] = normalize_orthography(value)

    # Ensure IPA fields are populated
    ensure_ipa_fields(
        result,
        orth_key="lemma",
        raw_key="ipa_raw_lemma",
        norm_key="ipa_normalized_lemma",
        tweak_fn=tweak_nominal_ipa,
    )
    ensure_ipa_fields(
        result,
        orth_key="plural",
        raw_key="ipa_raw_pl",
        norm_key="ipa_normalized_pl",
        tweak_fn=None,
    )

    # G2P fallback for plural IPA if still missing
    if not result.get("ipa_normalized_pl"):
        plural = result.get("plural", "")
        if plural:
            result["ipa_normalized_pl"] = normalize_ipa(
                to_ipa(plural), remove_stress=True
            )
            notes = result.get("notes", "")
            result["notes"] = (notes + " [IPA from G2P]").strip()

    # Validate POS and gender
    result["pos"] = (result.get("pos") or "").strip().upper()
    if result["pos"] == "N" and not result.get("gender"):
        return None

    # Run derivation pipeline
    derive_stem_final_and_cluster(result)
    derive_frontstem_dorsal(result)  # Mark gi/ge/ci/ce stems
    validate_plural_quality(result)
    derive_mutation_and_orth_change(result)
    derive_opportunity(result)
    derive_palatal_consonant_pl(result)
    derive_lemma_suffix(result)
    derive_target_is_suffix(result)
    derive_derived_verbs_fields(result)
    derive_derived_adj_fields(result)
    derive_nde_class(result)
    fix_nde_mutations(result)  # Fix mutation status for NDE items
    fix_underlying_palatals(result)  # Remove frontstem dorsals unless NDE
    derive_exception_reason(result)  # Must be after fix_underlying_palatals
    mark_tp_domain(
        result
    )  # Mark TP domain membership after all filters applied

    return result


OUTPUT_FIELDS = [
    "lemma",
    "gloss",
    "pos",
    "gender",
    "stem_final",
    "cluster",
    "frontstem_dorsal",
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
    "tp_in_domain",
    "target_is_suffix",
    "lemma_suffix",
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

        rows: list[dict[str, str]] = []
        # Track (lemma, pos, plural) -> index in rows list
        # so we can replace an entry if a later duplicate has
        # richer derivation data.
        seen_entries: dict[tuple[str, str, str], int] = {}
        duplicates_skipped = 0

        for i, row in enumerate(reader, start=1):
            if i % 1000 == 0:
                print(f"  Processed {i} rows...")
            processed = process_row(row)
            if processed is None:
                continue
            exploded = explode_derived_verbs_row(processed)

            for exp_row in exploded:
                key = (
                    exp_row.get("lemma", ""),
                    exp_row.get("pos", ""),
                    exp_row.get("plural", ""),
                )
                if key in seen_entries:
                    # Merge derivation data into the existing entry.
                    # Two sources of duplicate keys:
                    #   1. explode_derived_verbs_row split a multi-verb
                    #      input row into one row per verb (the common
                    #      case for nouns with several denominal verbs).
                    #   2. The raw CSV genuinely had two rows for the
                    #      same (lemma, pos, plural) — rare but possible.
                    # In both cases we want the UNION of the verbs, not
                    # the first one to arrive. Dropping the second was a
                    # latent bug that systematically lost ~80% of
                    # multi-verb data (e.g., fabrică keeping the plural
                    # noun "fabrici" while dropping the real infinitive
                    # "fabrica").
                    existing_idx = seen_entries[key]
                    existing = rows[existing_idx]
                    old_dv = (existing.get("derived_verbs") or "").strip()
                    new_dv = (exp_row.get("derived_verbs") or "").strip()
                    if new_dv:
                        old_verbs = [
                            v.strip() for v in old_dv.split("|") if v.strip()
                        ]
                        new_verbs = [
                            v.strip() for v in new_dv.split("|") if v.strip()
                        ]
                        merged = old_verbs + [
                            v for v in new_verbs if v not in old_verbs
                        ]
                        if merged != old_verbs:
                            existing["derived_verbs"] = "|".join(merged)
                            # Re-derive companion fields from the union.
                            derive_derived_verbs_fields(existing)
                    new_has_da = bool(exp_row.get("derived_adj"))
                    old_has_da = bool(existing.get("derived_adj"))
                    if new_has_da and not old_has_da:
                        existing["derived_adj"] = exp_row["derived_adj"]
                        existing["ipa_derived_adj"] = exp_row.get(
                            "ipa_derived_adj", ""
                        )
                    duplicates_skipped += 1
                    continue
                seen_entries[key] = len(rows)
                rows.append(exp_row)

        if duplicates_skipped > 0:
            print(f"  Skipped {duplicates_skipped} duplicate entries")

    # ── Post-processing: remove plural-as-lemma duplicates ──
    # Wiktionary sometimes lists plural forms as standalone lemmas
    # (e.g., "componente" listed separately from "componentă").
    # These get assigned their own plurals, creating phantom
    # mutations. Remove any lemma that equals the plural of a
    # DIFFERENT entry with the same pos.  Ochi-type items where
    # lemma == plural (e.g., "genunchi") must be kept.
    plural_to_lemma: dict[tuple[str, str], set[str]] = {}
    for r in rows:
        pl = (r.get("plural", "") or "").strip().lower()
        lem = (r.get("lemma", "") or "").strip().lower()
        pos = (r.get("pos", "") or "").strip()
        if pl:
            plural_to_lemma.setdefault((pl, pos), set()).add(lem)

    # Build set of lemmas that have their own entries
    own_entries: set[tuple[str, str]] = set()
    for r in rows:
        lem = (r.get("lemma", "") or "").strip().lower()
        pl = (r.get("plural", "") or "").strip().lower()
        pos = (r.get("pos", "") or "").strip()
        own_entries.add((lem, pos))

    def _is_plural_of_other(r: dict[str, str]) -> bool:
        lem = (r.get("lemma", "") or "").strip().lower()
        pos = (r.get("pos", "") or "").strip()
        key = (lem, pos)
        if key not in plural_to_lemma:
            return False
        # The lemma appears as a plural of some other entry.
        sources = plural_to_lemma[key]
        if sources == {lem}:
            # Only plural of itself (ochi-type) — keep.
            return False
        # It IS the plural of a different lemma.  But does
        # that different lemma ALSO have its own entry?  If
        # the base exists independently, this entry is
        # redundant.  Only remove if the base form is in the
        # dataset AND this lemma's own plural is empty or
        # looks like a double-plural (different from itself).
        own_pl = (r.get("plural", "") or "").strip().lower()
        # Keep if this lemma has a normal plural that differs
        # from both itself and the base's plural — it may be
        # a legitimate homophone.
        if own_pl and own_pl != lem:
            # Its plural is something else; check if ANY of
            # the base lemmas are truly different words.
            for src in sources:
                if src != lem and (src, pos) in own_entries:
                    return True
        elif not own_pl:
            # No plural at all — likely a Wiktionary plural
            # form page with no content.
            for src in sources:
                if src != lem:
                    return True
        return False

    pre_filter = len(rows)
    rows = [r for r in rows if not _is_plural_of_other(r)]
    plural_as_lemma = pre_filter - len(rows)
    if plural_as_lemma > 0:
        print(f"  Removed {plural_as_lemma} plural-as-lemma" " duplicates")

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

    opp_counts: dict[str, int] = {}
    for r in rows:
        opp = r.get("opportunity", "none")
        opp_counts[opp] = opp_counts.get(opp, 0) + 1

    print("\nOpportunity distribution:")
    for opp in ["i", "e", "uri", "none"]:
        print(f"  {opp}: {opp_counts.get(opp, 0)}")

    nde_gimpe = sum(
        1 for r in rows if r.get("exception_reason") == "nde:gimpe"
    )
    nde_ochi = sum(1 for r in rows if r.get("exception_reason") == "nde:ochi")
    nde_paduchi = sum(
        1 for r in rows if r.get("exception_reason") == "nde:paduchi"
    )
    true_exceptions = sum(
        1 for r in rows if r.get("exception_reason") == "unexplained"
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
