#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Attach Leipzig frequency counts to the Romanian lexicon CSV.

Usage:

  1. Set LEXICON_PATH, OUTPUT_PATH, FREQ_DIR as needed.
  2. Optionally set CORPUS_IDS to a specific list of corpus IDs,
     or leave it as None to use all *_freq.csv files in FREQ_DIR.
  3. Run:

       python scripts/attach_frequencies.py

Assumptions:
- Leipzig freq CSVs have columns: word, freq (or legacy: orth, freq).
- Lemmas in the lexicon are already project-normalized apart from NFC.
"""

import csv
import sys
import unicodedata
from pathlib import Path
from typing import Dict, List, Optional

HERE = Path(__file__).resolve()
PROJECT_ROOT = HERE.parents[1]

LEXICON_PATH: Path = PROJECT_ROOT / "data" / "romanian_lexicon_complete.csv"
OUTPUT_PATH: Path = PROJECT_ROOT / "data" / "romanian_lexicon_with_freq.csv"
FREQ_DIR: Path = PROJECT_ROOT / "data" / "leipzig" / "freq"

# If None auto-discover; otherwise set explicitly, e.g.:
# CORPUS_IDS = ["ron_wikipedia_2021_1M"]
CORPUS_IDS: Optional[List[str]] = None

# Keep imports aligned with rest of project layout.
sys.path.insert(0, str(PROJECT_ROOT / "src"))
try:
    from wiktionary_normalizer import normalize_orthography
except ImportError as exc:  # pragma: no cover
    print(
        "ERROR: Could not import normalize_orthography "
        "from wiktionary_normalizer.py"
    )
    print("PROJECT_ROOT:", PROJECT_ROOT)
    raise exc


def discover_corpora(freq_dir: Path) -> List[str]:
    """
    Discover available corpus IDs from freq_dir by listing *_freq.csv files.
    """
    return sorted(
        p.stem.removesuffix("_freq") for p in freq_dir.glob("*_freq.csv")
    )


def normalize_freq_key(token: str) -> str:
    """
    Normalize a token from the Leipzig frequency table for lookup.
    """
    t = unicodedata.normalize("NFC", token)
    return t.strip().lower()


def load_freq_table(freq_path: Path) -> Dict[str, int]:
    """
    Load a Leipzig frequency table from a CSV file.
    Returns a mapping from normalized word to frequency.
    """
    freq_map: Dict[str, int] = {}
    with freq_path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            raw = row.get("word") or ""
            word = normalize_freq_key(raw)
            if not word:
                continue
            try:
                freq = int((row.get("freq") or "").strip())
            except ValueError:
                continue
            if freq > 0:
                freq_map[word] = freq_map.get(word, 0) + freq
    return freq_map


def build_freq_maps(
    freq_dir: Path, corpus_ids: List[str]
) -> Dict[str, Dict[str, int]]:
    """
    Load frequency tables for the specified corpus IDs.
    Returns a mapping from corpus ID to its frequency table (word -> freq).
    """
    freq_maps: Dict[str, Dict[str, int]] = {}
    for corpus_id in corpus_ids:
        freq_path = freq_dir / f"{corpus_id}_freq.csv"
        if not freq_path.exists():
            print(
                f"WARNING: freq file not found for {corpus_id}: "
                f"{freq_path}, skipping."
            )
            continue
        print(f"Loading frequency table: {freq_path}")
        table = load_freq_table(freq_path)
        if not table:
            print(
                f"  WARNING: no entries parsed from {freq_path}, "
                "skipping corpus."
            )
            continue
        freq_maps[corpus_id] = table
    return freq_maps


def normalize_lemma(lemma: str) -> str:
    """Normalize lemma for frequency lookup using project normalization."""
    pre = unicodedata.normalize("NFC", lemma)
    norm = normalize_orthography(pre) or ""
    norm = unicodedata.normalize("NFC", norm)
    return norm.strip().lower()


def main():
    """Main entry point for attaching frequency data to the lexicon."""
    for label, path in (("lexicon", LEXICON_PATH), ("freq-dir", FREQ_DIR)):
        if not path.exists():
            print(f"ERROR: {label} not found: {path}")
            return

    requested_ids = (
        CORPUS_IDS if CORPUS_IDS is not None else discover_corpora(FREQ_DIR)
    )
    if not requested_ids:
        print(
            "No corpus frequency files found in freq-dir. Nothing to attach."
        )
        return

    freq_maps = build_freq_maps(FREQ_DIR, requested_ids)
    if not freq_maps:
        print("No usable frequency maps loaded. Exiting.")
        return

    corpus_ids = sorted(freq_maps)
    print("Using corpora:")
    print(*(f"  - {cid}" for cid in corpus_ids), sep="\n")
    print(f"\nReading lexicon: {LEXICON_PATH}")
    n_rows = 0
    with LEXICON_PATH.open("r", encoding="utf-8") as f_in, OUTPUT_PATH.open(
        "w", encoding="utf-8", newline=""
    ) as f_out:
        reader = csv.DictReader(f_in)
        base_fieldnames: List[str] = list(reader.fieldnames or [])
        extra_cols = [
            col
            for cid in corpus_ids
            for col in [f"freq_{cid}"]
            if col not in base_fieldnames
        ]
        out_fieldnames = base_fieldnames + extra_cols
        print(
            f"Writing output with {len(out_fieldnames)} columns "
            f"to: {OUTPUT_PATH}"
        )
        writer = csv.DictWriter(f_out, fieldnames=out_fieldnames)
        writer.writeheader()
        for row in reader:
            n_rows += 1
            key = normalize_lemma(row.get("lemma") or "")
            for cid in corpus_ids:
                # Frequencies stored as strings: "0" for missing,
                # else integer as string
                row[f"freq_{cid}"] = str(freq_maps[cid].get(key, 0))
            writer.writerow(row)
    print(f"Done. Wrote {n_rows} rows to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
