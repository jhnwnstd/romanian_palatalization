#!/usr/bin/env python3
"""Download and preprocess Leipzig Corpora frequency data for Romanian.

Pipeline per corpus:

1. Download ``{corpus_id}.tar.gz`` from the Leipzig server if not cached.
2. Extract into ``data/leipzig/corpora/{corpus_id}/`` using the safe
   ``filter='data'`` extraction policy.
3. Locate the ``*-words*.txt`` / ``*-words_pos_base*.txt`` table.
4. Build a cleaned ``word → freq`` map (filtering URLs, punctuation,
   multi-word strings) and write it as ``{corpus_id}_freq.csv``.

Contracts
---------
- The Leipzig archive layout puts the words file under
  ``<archive>/<corpus_id>/<corpus_id>-words*.txt``. We accept the
  alternative layout where the inner directory is missing.
- Tarballs are extracted with ``filter='data'`` (PEP 706) to refuse
  any member that would write outside the destination — required since
  Python 3.12 and a good idea before that.
- The output frequency CSV always has ``word,freq`` columns even if
  empty, so downstream loaders can rely on the header.
"""

from __future__ import annotations

import csv
import subprocess
import tarfile
import unicodedata
from pathlib import Path
from typing import Final

PROJECT_ROOT: Final[Path] = Path(__file__).resolve().parents[1]
DATA_DIR: Final[Path] = PROJECT_ROOT / "data" / "leipzig"

LEIPZIG_BASE_URL: Final[str] = (
    "https://downloads.wortschatz-leipzig.de/corpora"
)

RAW_DIR: Final[Path] = DATA_DIR / "raw"
CORPORA_DIR: Final[Path] = DATA_DIR / "corpora"
FREQ_DIR: Final[Path] = DATA_DIR / "freq"

DEFAULT_CORPORA: Final[tuple[str, ...]] = (
    # News
    # "ron_news_2015_1M",
    # "ron_news_2019_1M",
    # "ron_news_2020_1M",
    # "ron_news_2022_1M",
    "ron_news_2024_1M",
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
    # "ron_wikipedia_2021_1M",
)

# Original download page (for reference):
# https://wortschatz-leipzig.de/de/download/ron

# Punctuation we strip from the ends of tokens before normalising. Kept
# as a module constant rather than rebuilt per token.
_TOKEN_TRIM: Final[str] = ".,;:!?\"“”„'«»()[]{}<>"


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
    """Extract a corpus archive if not already extracted.

    Uses the ``filter='data'`` extraction policy (PEP 706), which
    refuses members containing absolute paths, ``..`` traversal, or
    device/named-pipe nodes. The behavior is mandatory under Python
    3.12+ and a sensible default before that for untrusted archives.
    """
    dest_dir = CORPORA_DIR / corpus_id
    ensure_dir(CORPORA_DIR)
    if dest_dir.exists() and any(dest_dir.iterdir()):
        print(f"  Using existing extracted directory: {dest_dir}")
        return dest_dir
    print(f"  Extracting {archive} -> {dest_dir}")
    ensure_dir(dest_dir)
    with tarfile.open(archive, "r:gz") as tf:
        tf.extractall(dest_dir, filter="data")
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
    """Normalize and filter a Leipzig token.

    Returns the empty string for tokens we drop:

    - empty after trimming punctuation
    - containing whitespace (multi-word strings)
    - URLs / email-like (``http``, ``www.``, ``://``, ``@``)
    - no alphabetic character

    Otherwise returns the NFC-normalized, lowercased form.
    """
    t = unicodedata.normalize("NFC", token).strip()
    t = t.strip(_TOKEN_TRIM)
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


def parse_words_file(path: Path) -> dict[str, int]:
    """Parse a Leipzig words file and build a frequency table."""
    freq_map: dict[str, int] = {}
    print(f"  Parsing words file: {path}")
    with path.open("r", encoding="utf-8") as f:
        for i, line in enumerate(f, 1):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            # Skip a header row (Leipzig files sometimes ship one).
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


def write_freq_csv(corpus_id: str, freq_map: dict[str, int]) -> Path:
    """Write a frequency table to CSV (``word,freq``)."""
    ensure_dir(FREQ_DIR)
    out_path = FREQ_DIR / f"{corpus_id}_freq.csv"
    print(f"  Writing frequency CSV: {out_path}")
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["word", "freq"])
        # Sort by descending frequency, then alphabetically for ties,
        # so the output file is deterministic across runs.
        for word, freq in sorted(
            freq_map.items(), key=lambda kv: (-kv[1], kv[0])
        ):
            writer.writerow([word, freq])
    return out_path


def main() -> int:
    """Download, extract, and build frequency tables for each corpus.

    Returns 0 if every selected corpus completed; 1 if any failed.
    Continues across failures so a single broken corpus doesn't abort
    the others.
    """
    print("Selected corpora:")
    print(*(f"  - {cid}" for cid in DEFAULT_CORPORA), sep="\n")
    print()
    failures = 0
    for corpus_id in DEFAULT_CORPORA:
        print(f"Processing corpus: {corpus_id}")
        try:
            archive = ensure_downloaded(corpus_id)
            corpus_dir = ensure_extracted(corpus_id, archive)
            words_path = find_words_file(corpus_dir, corpus_id)
            freq_map = parse_words_file(words_path)
            write_freq_csv(corpus_id, freq_map)
        except Exception as exc:  # noqa: B902
            # Broad on purpose: a single corpus failure shouldn't
            # abort processing of the others.
            print(f"ERROR while processing {corpus_id}: {exc}")
            failures += 1
        print()
    print(f"Done. ({failures} failure(s))")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
