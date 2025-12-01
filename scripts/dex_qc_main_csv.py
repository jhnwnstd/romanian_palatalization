#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Surgical QC assistant for romanian_lexicon_raw.csv

Uses DEX as an external oracle to fill missing data and resolve uncertainties.
Only queries DEX for rows that need help; doesn't overwrite good data.

Runs BEFORE romanian_processor_main.py, so does not derive
palatalization fields. The processor will compute mutation, orth_change,
opportunity, etc. from scratch.

Usage:
    python dex_qc_main_csv.py
"""

import csv
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from dex_utils import (  # noqa: E402  # pylint: disable=wrong-import-position
    DexEntry,
    best_dex_noun_header,
    dex_has_entry,
    fetch_dex_page,
    norm_lemma,
    save_disk_cache,
)

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

# Base project dir = repo root (parent of this file's directory)
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"

# Default locations for pipeline
INPUT_CSV = str(DATA_DIR / "romanian_lexicon_raw.csv")
OUTPUT_CSV = str(DATA_DIR / "romanian_lexicon_raw_dex.csv")
AUDIT_CSV = str(DATA_DIR / "dex_qc_audit.csv")
DISAGREEMENTS_CSV = str(DATA_DIR / "dex_qc_disagreements.csv")


# ------------------------------------------------------------
# Row classification
# ------------------------------------------------------------
@dataclass
class QCTarget:
    """Identifies what needs fixing in a row."""

    row_idx: int
    lemma: str
    needs_gender: bool = False
    needs_plural: bool = False
    has_uncertainty_note: bool = False
    reason: str = ""


def row_needs_qc(row: Dict[str, str], idx: int) -> Optional[QCTarget]:
    """
    Determine if a row needs QC assistance from DEX.

    Criteria (tailored to our dataset):
    - Missing gender or plural
    - Harvester flags in notes:
        * "needs plural confirmation"
        * "needs gender confirmation"
    - Plural is malformed or obviously harvested meta-text
      (placeholders, Wiktionary gloss fragments, truncated -ț, etc.)
    """
    lemma = norm_lemma(row.get("lemma") or "")
    if not lemma:
        return None
    pos = (row.get("pos") or "").strip().upper()
    if not pos:
        return None
    notes_raw = row.get("notes") or ""
    source_raw = row.get("source") or ""
    notes = notes_raw.lower()
    source = source_raw.lower()
    gender = (row.get("gender") or "").strip()
    plural = (row.get("plural") or "").strip()
    plural_lower = plural.lower()
    if "dex-confirmed" in notes or "dex_qc" in source:
        if gender and plural:
            if (
                plural not in ("-", "—", "–", "?", "n/a")
                and "adj" not in plural_lower
            ):
                return None
    has_plural_note = "needs plural confirmation" in notes
    has_gender_note = "needs gender confirmation" in notes
    has_uncertainty_note = has_plural_note or has_gender_note
    needs_gender = not gender or has_gender_note
    needs_plural = not plural or has_plural_note
    # Allow uncountable nouns to be flagged but not treated as errors
    invariant_hint = any(
        hint in notes for hint in ("invariant", "indeclinable", "uncountable")
    )
    is_invariant_candidate = (not plural) and invariant_hint
    if plural:
        if plural in ("-", "—", "–", "?", "n/a"):
            needs_plural = True
        if plural_lower.endswith(("ț", "ţ")) and not plural_lower.endswith(
            ("ți", "ţi")
        ):
            needs_plural = True
        if any(ch.isspace() for ch in plural) or any(
            marker in plural_lower
            for marker in (
                " adj",
                "verb",
                "obsolete",
                "archaic",
                "[",
                "]",
                "(",
                ")",
            )
        ):
            needs_plural = True
        if len(plural) <= 2 and plural_lower not in ("e", "i", "ii"):
            needs_plural = True
    if (
        needs_gender
        or needs_plural
        or has_uncertainty_note
        or is_invariant_candidate
    ):
        reasons: List[str] = []
        if needs_gender:
            reasons.append("missing_gender")
        if needs_plural:
            reasons.append("missing_plural")
        if has_uncertainty_note:
            reasons.append("uncertainty_note")
        if is_invariant_candidate:
            reasons.append("invariant_candidate")
        return QCTarget(
            row_idx=idx,
            lemma=lemma,
            needs_gender=needs_gender,
            needs_plural=needs_plural,
            has_uncertainty_note=has_uncertainty_note,
            reason=" | ".join(reasons),
        )
    return None


# ------------------------------------------------------------
# DEX repair logic
# ------------------------------------------------------------
def is_safe_orthographic_repair(csv_plural: str, dex_plural: str) -> bool:
    """
    Detect if DEX plural is just a safe orthographic fix
    (e.g., truncated -ț → -ți).

    Safe repairs are ones where the phonological content is clearly the same.

    Examples:
    - actanț → actanți (missing final i)
    - israeliț → israeliți (missing final i)
    """
    if not csv_plural or not dex_plural:
        return False
    csv_lower = csv_plural.lower()
    dex_lower = dex_plural.lower()
    # Truncated -ț → -ți pattern
    if csv_lower.endswith(("ț", "ţ")) and not csv_lower.endswith(("ți", "ţi")):
        if dex_lower.endswith(("ți", "ţi")):
            dex_candidates = (
                csv_lower + "i",
                csv_lower.replace("ţ", "ț") + "i",
            )
            return dex_lower in dex_candidates
    return False


def apply_dex_repair(
    row: Dict[str, str], target: QCTarget, dex_entry: DexEntry
) -> Dict[str, Any]:
    """
    Apply DEX data to fill missing or flagged-as-uncertain fields.

    Conservative approach:
    - Fill if field is empty/whitespace OR flagged in notes
    - Log all changes for audit

    This version does NOT derive palatalization fields - the processor
    will handle that.
    """
    repaired: Dict[str, Any] = row.copy()
    changes: List[str] = []
    filled_plural = False
    is_safe_repair = False
    if target.needs_gender and dex_entry.gender:
        old_gender = (row.get("gender") or "").strip()
        if dex_entry.gender != old_gender:
            repaired["gender"] = dex_entry.gender
            changes.append(f"gender:{dex_entry.gender}")
    if target.needs_plural and dex_entry.plural and dex_entry.plural != "-":
        old_plural = (row.get("plural") or "").strip()
        is_safe_repair = is_safe_orthographic_repair(
            old_plural, dex_entry.plural
        )
        if dex_entry.plural != old_plural:
            repaired["plural"] = dex_entry.plural
            changes.append(f"plural:{dex_entry.plural}")
            filled_plural = True
            if is_safe_repair:
                changes.append("safe-orthographic-repair")
    if changes:
        existing_notes = (row.get("notes") or "").strip()
        for change in changes:
            if change.startswith("gender:"):
                existing_notes = existing_notes.replace(
                    "needs gender confirmation", ""
                )
            if change.startswith("plural:"):
                existing_notes = existing_notes.replace(
                    "needs plural confirmation", ""
                )
        parts = [p.strip() for p in existing_notes.split("|") if p.strip()]
        cleaned_notes = " | ".join(parts) if parts else ""
        dex_note = f"DEX-confirmed: {', '.join(changes)}"
        if cleaned_notes:
            repaired["notes"] = f"{cleaned_notes} | {dex_note}"
        else:
            repaired["notes"] = dex_note
        if filled_plural:
            src = (row.get("source") or "").strip()
            if src:
                if "dex_qc" not in src.lower():
                    repaired["source"] = f"{src} | dex_qc"
            else:
                repaired["source"] = "dex_qc"
    return repaired


def detect_disagreement(
    row: Dict[str, str], dex_entry: DexEntry
) -> Optional[Dict[str, Any]]:
    """
    Detect disagreements between CSV and DEX for non-empty fields.

    Returns dict with disagreement details if found, else None.
    """
    disagreements = []
    csv_gender = (row.get("gender") or "").strip()
    if csv_gender and dex_entry.gender and csv_gender != dex_entry.gender:
        disagreements.append(
            {
                "field": "gender",
                "csv_value": csv_gender,
                "dex_value": dex_entry.gender,
            }
        )
    csv_plural = norm_lemma(row.get("plural") or "")
    dex_plural = norm_lemma(dex_entry.plural)
    if csv_plural and dex_plural != "-" and csv_plural != dex_plural:
        disagreements.append(
            {
                "field": "plural",
                "csv_value": csv_plural,
                "dex_value": dex_plural,
            }
        )
    if disagreements:
        return {
            "lemma": row.get("lemma"),
            "disagreements": disagreements,
            "dex_url": dex_entry.url,
        }
    return None


# ------------------------------------------------------------
# Main QC pipeline helpers
# ------------------------------------------------------------
def _create_audit_record(
    target: QCTarget, status: str, changes: str = "", dex_url: str = ""
) -> Dict[str, Any]:
    """Create a single audit record."""
    return {
        "lemma": target.lemma,
        "status": status,
        "reason": target.reason,
        "changes": changes,
        "dex_url": dex_url,
    }


def _process_no_dex_entry(
    targets: List[QCTarget],
) -> tuple[List[Dict[str, Any]], int]:
    """Handle case where DEX has no entry for lemma."""
    audit_records = []
    for target in targets:
        audit_records.append(_create_audit_record(target, "no_dex_entry"))
    return audit_records, len(targets)


def _process_targets_for_lemma(
    targets: List[QCTarget],
    rows: List[Dict[str, str]],
    dex_entry: DexEntry,
) -> tuple[List[Dict[str, Any]], List[Dict[str, Any]], int]:
    """Process all targets for a given lemma."""
    audit_records = []
    disagreement_records = []
    repaired_count = 0
    for target in targets:
        disagreement = detect_disagreement(rows[target.row_idx], dex_entry)
        if disagreement:
            disagreement_records.append(disagreement)
        original_row = rows[target.row_idx].copy()
        repaired_row = apply_dex_repair(
            rows[target.row_idx], target, dex_entry
        )
        changes = []
        for field in ["gender", "plural"]:
            old_val = original_row.get(field, "").strip()
            new_val = repaired_row.get(field, "").strip()
            if old_val != new_val:
                changes.append(f"{field}: '{old_val}' → '{new_val}'")
        if changes:
            rows[target.row_idx] = repaired_row
            repaired_count += 1
            audit_records.append(
                _create_audit_record(
                    target, "repaired", " | ".join(changes), dex_entry.url
                )
            )
        else:
            audit_records.append(
                _create_audit_record(
                    target, "no_changes_needed", "", dex_entry.url
                )
            )
    return audit_records, disagreement_records, repaired_count


def _process_error(
    targets: List[QCTarget], exc: Exception
) -> tuple[List[Dict[str, Any]], int]:
    """Handle errors during DEX query."""
    audit_records = []
    for target in targets:
        audit_records.append(
            _create_audit_record(target, f"error:{type(exc).__name__}")
        )
    return audit_records, len(targets)


def _qc_derivational_fields(
    rows: List[Dict[str, str]],
    audit_records: List[Dict[str, Any]],
    max_extra_queries: Optional[int] = None,
) -> None:
    """
    Light-weight QC on derived_verbs / derived_adj using DEX.

    Strategy:
      - For each lemma with non-empty derived_verbs / derived_adj:
          * Split the pipe-separated list.
          * For each derived lemma, check whether DEX has any entry.
          * If DEX clearly has no entry for that lemma, drop it.
      - Log any pruning into the existing audit log.

    We never drop data just because of network errors or exhausted budget.
    """

    # Cache: derived lemma -> bool (True = DEX has entry, False = clearly no)
    dex_presence_cache: Dict[str, bool] = {}
    queries_used = 0

    def _has_dex(lemma: str) -> bool:
        nonlocal queries_used
        lemma_norm = norm_lemma(lemma)
        if not lemma_norm:
            return False

        if lemma_norm in dex_presence_cache:
            return dex_presence_cache[lemma_norm]

        if max_extra_queries is not None and queries_used >= max_extra_queries:
            # Budget exhausted: be conservative and KEEP it
            dex_presence_cache[lemma_norm] = True
            return True

        try:
            present = dex_has_entry(lemma_norm, restrict_to_keep_pos=False)
        except Exception:
            # On any error, don't nuke the data
            present = True

        dex_presence_cache[lemma_norm] = present
        queries_used += 1
        return present

    for row in rows:
        base_lemma = row.get("lemma") or ""
        for field in ("derived_verbs", "derived_adj"):
            raw = (row.get(field) or "").strip()
            if not raw:
                continue

            parts = [p.strip() for p in raw.split("|")]
            lemmas = [norm_lemma(p) for p in parts if norm_lemma(p)]
            if not lemmas:
                continue

            kept: List[str] = []
            dropped: List[str] = []

            for dlemma in lemmas:
                if _has_dex(dlemma):
                    kept.append(dlemma)
                else:
                    dropped.append(dlemma)

            # Nothing pruned → nothing to log
            if not dropped:
                continue

            row[field] = " | ".join(kept) if kept else ""

            audit_records.append(
                {
                    "lemma": base_lemma,
                    "status": "deriv_qc",
                    "reason": f"{field}_dex_filter",
                    "changes": (
                        f"{field}: kept [{', '.join(kept)}], "
                        f"dropped [{', '.join(dropped)}]"
                    ),
                    "dex_url": "",
                }
            )


# ------------------------------------------------------------
# Main QC pipeline
# ------------------------------------------------------------
def run_qc(  # pylint: disable=too-many-locals
    input_csv: str = INPUT_CSV,
    output_csv: str = OUTPUT_CSV,
    audit_csv: str = AUDIT_CSV,
    disagreements_csv: str = DISAGREEMENTS_CSV,
    max_queries: Optional[int] = None,
) -> None:
    """
    Run surgical QC on romanian_lexicon_raw.csv.

    Steps:
    1. Load CSV and normalize POS field
    2. Identify rows needing help
    3. Query DEX only for those rows
    4. Fill missing data conservatively
    5. Log all changes and disagreements
    6. Save cleaned CSV
    """
    with open(input_csv, "r", encoding="utf-8") as f_in:
        reader = csv.DictReader(f_in)
        fieldnames = reader.fieldnames or []
        rows = list(reader)
    print(f"Loaded {len(rows)} rows from {input_csv}")
    for row in rows:
        row["pos"] = (row.get("pos") or "").strip().upper()
    targets: List[QCTarget] = [
        target
        for idx, row in enumerate(rows)
        if (target := row_needs_qc(row, idx))
    ]
    print(f"Found {len(targets)} rows needing QC")
    if max_queries:
        targets = targets[:max_queries]
        print(f"Limited to {len(targets)} queries")
    # Deduplicate by lemma to minimize expensive DEX API calls
    targets_by_lemma: Dict[str, List[QCTarget]] = defaultdict(list)
    for target in targets:
        targets_by_lemma[target.lemma].append(target)
    unique_lemmas = list(targets_by_lemma.keys())
    print(
        f"Unique lemmas to query: {len(unique_lemmas)} "
        f"(down from {len(targets)} rows)"
    )
    audit_records: List[Dict[str, Any]] = []
    disagreement_records: List[Dict[str, Any]] = []
    repaired_count = 0
    failed_count = 0
    for i, lemma in enumerate(unique_lemmas, 1):
        if i % 25 == 0:
            print(f"Processing {i}/{len(unique_lemmas)}...")
        try:
            headers, _ = fetch_dex_page(lemma)
            dex_entry = best_dex_noun_header(headers, lemma)
            if not dex_entry:
                new_audit, fail_count = _process_no_dex_entry(
                    targets_by_lemma[lemma]
                )
                audit_records.extend(new_audit)
                failed_count += fail_count
                continue
            new_audit, new_disagreements, repairs = _process_targets_for_lemma(
                targets_by_lemma[lemma], rows, dex_entry
            )
            audit_records.extend(new_audit)
            disagreement_records.extend(new_disagreements)
            repaired_count += repairs
        except Exception as exc:  # pylint: disable=broad-exception-caught
            new_audit, fail_count = _process_error(
                targets_by_lemma[lemma], exc
            )
            audit_records.extend(new_audit)
            failed_count += fail_count
    _qc_derivational_fields(
        rows,
        audit_records,
        max_extra_queries=None,  # or an int if you want a separate budget
    )
    with open(output_csv, "w", encoding="utf-8", newline="") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    with open(audit_csv, "w", encoding="utf-8", newline="") as f_audit:
        audit_fields = ["lemma", "status", "reason", "changes", "dex_url"]
        writer = csv.DictWriter(f_audit, fieldnames=audit_fields)
        writer.writeheader()
        for rec in audit_records:
            writer.writerow({k: rec.get(k, "") for k in audit_fields})
    if disagreement_records:
        with open(
            disagreements_csv, "w", encoding="utf-8", newline=""
        ) as f_dis:
            csv_writer = csv.writer(f_dis)
            csv_writer.writerow(
                ["lemma", "field", "csv_value", "dex_value", "dex_url"]
            )
            for rec in disagreement_records:
                for dis in rec["disagreements"]:
                    csv_writer.writerow(
                        [
                            rec["lemma"],
                            dis["field"],  # type: ignore[index]
                            dis["csv_value"],  # type: ignore[index]
                            dis["dex_value"],  # type: ignore[index]
                            rec.get("dex_url", ""),
                        ]
                    )
    save_disk_cache()
    print("\n=== DEX QC Summary ===")
    print(f"Total rows:            {len(rows)}")
    print(f"Rows needing QC:       {len(targets)}")
    print(f"Successfully repaired: {repaired_count}")
    print(f"Failed/no data:        {failed_count}")
    print(f"Disagreements found:   {len(disagreement_records)}")
    print("\nOutput files:")
    print(f"  Cleaned CSV:         {output_csv}")
    print(f"  Audit log:           {audit_csv}")
    if disagreement_records:
        print(f"  Disagreements:       {disagreements_csv}")


# ------------------------------------------------------------
# CLI entry point
# ------------------------------------------------------------
if __name__ == "__main__":
    max_queries_arg: Optional[int] = None
    if len(sys.argv) > 1:
        try:
            max_queries_arg = int(sys.argv[1])
            print(f"Limiting to {max_queries_arg} DEX queries")
        except ValueError:
            print(f"Usage: {sys.argv[0]} [max_queries]")
            sys.exit(1)
    run_qc(
        input_csv=INPUT_CSV,
        output_csv=OUTPUT_CSV,
        audit_csv=AUDIT_CSV,
        disagreements_csv=DISAGREEMENTS_CSV,
        max_queries=max_queries_arg,
    )
