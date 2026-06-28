#!/usr/bin/env python3
"""Attach Leipzig frequency counts to the Romanian lexicon CSV.

Reads each Leipzig corpus `*_freq.csv` in :data:`FREQ_DIR`, normalizes
the word keys to match our lexicon convention (NFC + cedilla → comma-
below), and writes a copy of the lexicon CSV with one extra
``freq_{corpus_id}`` column per corpus.

Contracts
---------
- Leipzig freq CSVs have columns ``word`` and ``freq`` (integer).
- Lemmas in the lexicon are already project-normalized apart from NFC;
  we run :func:`normalize_orthography` on the lookup key.
- Missing entries are written as ``"0"`` (string), not blank, so the
  downstream R/Python loaders can parse the column as integer
  uniformly.
"""

from __future__ import annotations

import csv
import sys
import unicodedata
from pathlib import Path
from typing import Final

HERE: Final[Path] = Path(__file__).resolve()
PROJECT_ROOT: Final[Path] = HERE.parents[1]

LEXICON_PATH: Final[Path] = (
    PROJECT_ROOT / "data" / "romanian_lexicon_complete.csv"
)
OUTPUT_PATH: Final[Path] = (
    PROJECT_ROOT / "data" / "romanian_lexicon_with_freq.csv"
)
FREQ_DIR: Final[Path] = PROJECT_ROOT / "data" / "leipzig" / "freq"

# If None: auto-discover from FREQ_DIR. Set explicitly to restrict to a
# specific subset, e.g. CORPUS_IDS = ["ron_wikipedia_2021_1M"].
CORPUS_IDS: Final[list[str] | None] = None

# Keep imports aligned with rest of project layout.
sys.path.insert(0, str(PROJECT_ROOT / "src"))
try:
    from wiktionary_normalizer import normalize_cedilla, normalize_orthography
except ImportError as exc:  # pragma: no cover
    print(
        "ERROR: Could not import normalize_orthography "
        "from wiktionary_normalizer.py"
    )
    print("PROJECT_ROOT:", PROJECT_ROOT)
    raise exc


def discover_corpora(freq_dir: Path) -> list[str]:
    """Return sorted corpus IDs by listing ``*_freq.csv`` in *freq_dir*."""
    return sorted(
        p.stem.removesuffix("_freq") for p in freq_dir.glob("*_freq.csv")
    )


def normalize_freq_key(token: str) -> str:
    """
    Normalize a token from the Leipzig frequency table for lookup.

    Applies cedilla-to-comma normalization to match the lexicon side,
    which uses normalize_orthography (ş→ș, ţ→ț).
    """
    t = unicodedata.normalize("NFC", token)
    t = normalize_cedilla(t)
    return t.strip().lower()


def load_freq_table(freq_path: Path) -> dict[str, int]:
    """Load a Leipzig frequency table, summing any duplicate keys.

    Returns a mapping from normalized word → frequency. Frequencies of
    0 are dropped (keeps the map small). Rows whose ``freq`` field
    doesn't parse as int are skipped silently.
    """
    freq_map: dict[str, int] = {}
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
                # Sum (rather than overwrite) so duplicate keys in a
                # source corpus don't silently lose counts.
                freq_map[word] = freq_map.get(word, 0) + freq
    return freq_map


def build_freq_maps(
    freq_dir: Path, corpus_ids: list[str]
) -> dict[str, dict[str, int]]:
    """Load frequency tables for the specified corpus IDs."""
    freq_maps: dict[str, dict[str, int]] = {}
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
    norm = normalize_orthography(lemma) or ""
    return norm.strip()


def main() -> int:
    """Run the attach step. Returns 0 on success, 1 on input error."""
    for label, path in (("lexicon", LEXICON_PATH), ("freq-dir", FREQ_DIR)):
        if not path.exists():
            print(f"ERROR: {label} not found: {path}")
            return 1

    requested_ids = (
        CORPUS_IDS if CORPUS_IDS is not None else discover_corpora(FREQ_DIR)
    )
    if not requested_ids:
        print(
            "No corpus frequency files found in freq-dir. "
            "Nothing to attach."
        )
        return 1

    freq_maps = build_freq_maps(FREQ_DIR, requested_ids)
    if not freq_maps:
        print("No usable frequency maps loaded. Exiting.")
        return 1

    corpus_ids = sorted(freq_maps)
    print("Using corpora:")
    print(*(f"  - {cid}" for cid in corpus_ids), sep="\n")
    print(f"\nReading lexicon: {LEXICON_PATH}")
    n_rows = 0
    with (
        LEXICON_PATH.open("r", encoding="utf-8") as f_in,
        OUTPUT_PATH.open("w", encoding="utf-8", newline="") as f_out,
    ):
        reader = csv.DictReader(f_in)
        base_fieldnames: list[str] = list(reader.fieldnames or [])
        extra_cols = [
            f"freq_{cid}"
            for cid in corpus_ids
            if f"freq_{cid}" not in base_fieldnames
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
                # Missing → "0" so downstream loaders can parse the
                # column uniformly as integer.
                row[f"freq_{cid}"] = str(freq_maps[cid].get(key, 0))
            writer.writerow(row)
    print(f"Done. Wrote {n_rows} rows to {OUTPUT_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
