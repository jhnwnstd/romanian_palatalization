#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Download and preprocess Leipzig Corpora frequency data for Romanian.

- Downloads one or more fixed Romanian corpora (news/newscrawl/web/Wikipedia)
- Extracts archives under: data/leipzig/corpora/<corpus_id>/
- Finds the *-words*.txt / *-words_pos_base*.txt file
- Builds a cleaned frequency table (word -> frequency) as CSV

Usage:
    python scripts/download_leipzig.py
"""

import csv
import subprocess
import tarfile
import unicodedata
from pathlib import Path
from typing import Dict

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data" / "leipzig"

LEIPZIG_BASE_URL = "https://downloads.wortschatz-leipzig.de/corpora"

RAW_DIR = DATA_DIR / "raw"
CORPORA_DIR = DATA_DIR / "corpora"
FREQ_DIR = DATA_DIR / "freq"

DEFAULT_CORPORA = [
    # News
    # "ron_news_2015_1M",
    # "ron_news_2019_1M",
    # "ron_news_2020_1M",
    # "ron_news_2022_1M",
    # "ron_news_2024_1M",
    # Newscrawl
    # "ron_newscrawl_2011_1M",
    # "ron_newscrawl_2015_1M",
    # Web (2015 Romania + older web)
    # "ron_web_2011_1M",
    # "ron-md_web_2013_1M",
    # "ron-md_web_2014_1M",
    # "ron-ro_web_2015_1M",
    # Wikipedia
    # "ron_wikipedia_2010_1M",
    # "ron_wikipedia_2011_1M",
    # "ron_wikipedia_2014_1M",
    # "ron_wikipedia_2018_1M",
    "ron_wikipedia_2021_1M",
]

# Original download page (for reference):
# https://wortschatz-leipzig.de/de/download/ron


def ensure_dir(path: Path) -> None:
    """Ensure a directory exists, creating it if necessary."""
    path.mkdir(parents=True, exist_ok=True)


def run_curl(url: str, dest: Path) -> None:
    """Download a file using curl."""
    print(f"  Downloading {url} -> {dest}")
    subprocess.run(["curl", "-L", "-o", str(dest), url], check=True)


def ensure_downloaded(corpus_id: str) -> Path:
    """Download a corpus archive if not already cached."""
    ensure_dir(RAW_DIR)
    archive = RAW_DIR / f"{corpus_id}.tar.gz"
    if archive.exists():
        print(f"  Archive already present: {archive}")
        return archive
    url = f"{LEIPZIG_BASE_URL}/{corpus_id}.tar.gz"
    run_curl(url, archive)
    return archive


def ensure_extracted(corpus_id: str, archive: Path) -> Path:
    """Extract a corpus archive if not already extracted."""
    dest_dir = CORPORA_DIR / corpus_id
    ensure_dir(CORPORA_DIR)
    if dest_dir.exists() and any(dest_dir.iterdir()):
        print(f"  Using existing extracted directory: {dest_dir}")
        return dest_dir
    print(f"  Extracting {archive} -> {dest_dir}")
    ensure_dir(dest_dir)
    with tarfile.open(archive, "r:gz") as tf:
        tf.extractall(dest_dir)
    return dest_dir


def find_words_file(corpus_dir: Path, corpus_id: str) -> Path:
    """Find the words frequency file in an extracted corpus."""
    # Leipzig archives currently unpack as:
    # <corpus_dir>/<corpus_id>/<corpus_id>-words*.txt
    inner = corpus_dir / corpus_id
    base = inner if inner.is_dir() else corpus_dir
    for suffix in ("-words_pos_base.txt", "-words.txt"):
        candidate = base / f"{corpus_id}{suffix}"
        if candidate.exists():
            print(
                f"  Found words file: " f"{candidate.relative_to(corpus_dir)}"
            )
            return candidate
    # Fallback in case the layout shifts slightly in future
    # releases.
    patterns = (
        f"{corpus_id}-words_pos_base.txt",
        f"{corpus_id}-words.txt",
        "*-words_pos_base.txt",
        "*-words.txt",
    )
    for pattern in patterns:
        matches = list(base.rglob(pattern))
        if matches:
            chosen = min(
                matches, key=lambda p: len(p.relative_to(corpus_dir).parts)
            )
            print(
                f"  Found words file (fallback): "
                f"{chosen.relative_to(corpus_dir)}"
            )
            return chosen
    raise FileNotFoundError(
        f"No words or words_pos_base file found for "
        f"{corpus_id} in {corpus_dir}"
    )


def clean_token(token: str) -> str:
    """Clean and normalize a token, filtering out junk and noise."""
    # Filters out punctuation-anchored junk, multi-word strings, and
    # obvious noise like URLs.
    t = unicodedata.normalize("NFC", token).strip()
    t = t.strip(".,;:!?\"“”„'«»()[]{}<>")
    if not t:
        return ""
    if any(ch.isspace() for ch in t):
        return ""
    lower = t.lower()
    if any(x in lower for x in ("http", "www.", "://", "@")):
        return ""
    if not any(ch.isalpha() for ch in t):
        return ""
    return lower


def parse_words_file(path: Path) -> Dict[str, int]:
    """Parse a Leipzig words file and build a frequency table."""
    freq_map: Dict[str, int] = {}
    print(f"  Parsing words file: {path}")
    with path.open("r", encoding="utf-8") as f:
        for i, line in enumerate(f, 1):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            if i == 1 and not parts[0].isdigit():
                continue
            freq_str = parts[2]
            if not freq_str.isdigit():
                continue
            token = clean_token(parts[1])
            if not token:
                continue
            freq = int(freq_str)
            if freq <= 0:
                continue
            freq_map[token] = freq_map.get(token, 0) + freq
    print(f"  Parsed {len(freq_map):,} word types from {path.name}")
    return freq_map


def write_freq_csv(corpus_id: str, freq_map: Dict[str, int]) -> Path:
    """Write a frequency table to CSV format."""
    ensure_dir(FREQ_DIR)
    out_path = FREQ_DIR / f"{corpus_id}_freq.csv"
    print(f"  Writing frequency CSV: {out_path}")
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["word", "freq"])
        for word, freq in sorted(
            freq_map.items(), key=lambda x: (-x[1], x[0])
        ):
            # freq is int here, csv.writer converts to string automatically
            writer.writerow([word, freq])
    return out_path


def main() -> None:
    """Main entry point for downloading and processing Leipzig corpora."""
    print("Selected corpora:")
    print(*(f"  - {cid}" for cid in DEFAULT_CORPORA), sep="\n")
    print()
    for corpus_id in DEFAULT_CORPORA:
        print(f"Processing corpus: {corpus_id}")
        try:
            archive = ensure_downloaded(corpus_id)
            corpus_dir = ensure_extracted(corpus_id, archive)
            words_path = find_words_file(corpus_dir, corpus_id)
            freq_map = parse_words_file(words_path)
            write_freq_csv(corpus_id, freq_map)
        except (
            Exception
        ) as exc:  # keep broad; we want to continue over partial failures
            print(f"ERROR while processing {corpus_id}: {exc}")
        print()
    print("Done.")


if __name__ == "__main__":
    main()
