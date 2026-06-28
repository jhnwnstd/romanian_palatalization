#!/usr/bin/env python3
"""Supplement ``derived_verbs`` in the DEX-QC'd lexicon.

For each row whose ``derived_verbs`` is empty AND ``stem_final`` is one
of the target consonants, generate plausible Romanian denominal-verb
infinitives by morphological rule and ask DEX whether each candidate
actually exists. Confirmed verbs get appended to the row's
``derived_verbs`` column for the downstream processor to consume.

Inputs / outputs
----------------
- read  ``data/romanian_lexicon_with_freq.csv`` (frequency + stem info)
- read  ``data/romanian_lexicon_raw_dex.csv`` (the row we'll mutate)
- write ``data/romanian_lexicon_raw_dex_supplemented.csv``

After this script completes, the pipeline runs stages 4–6 on the
supplemented CSV.

Concurrency
-----------
Workers are bounded by ``N_WORKERS`` and share the DEX disk cache.
Cache mutation is guarded by ``_CACHE_LOCK`` in ``dex_utils``;
``save_disk_cache`` writes atomically via tmp-file + rename. The
periodic ``save_every`` checkpoint keeps the cache crash-safe across
long runs.
"""

from __future__ import annotations

import csv
import re
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Final

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from dex_utils import (  # noqa: E402
    DISK_CACHE,
    HTML_CACHE,
    norm_lemma,
    polite_get,
    save_disk_cache,
)

REPO: Final[Path] = Path(__file__).parent.parent
# Input: with_freq has stem_final / opportunity computed; we use it to
# *identify* target lemmas. The actual write-back is to raw_dex
# (upstream of stage 4) so the pipeline can recompute everything from
# the new derived_verbs.
WITH_FREQ_CSV: Final[Path] = REPO / "data" / "romanian_lexicon_with_freq.csv"
RAW_DEX_CSV: Final[Path] = REPO / "data" / "romanian_lexicon_raw_dex.csv"
OUTPUT_CSV: Final[Path] = (
    REPO / "data" / "romanian_lexicon_raw_dex_supplemented.csv"
)

TARGET_STEMS: Final[frozenset[str]] = frozenset({"c", "g", "t", "d", "s", "z"})

# Limit to the N highest-frequency eligible lemmas. None = process all.
TOP_N_BY_FREQUENCY: Final[int | None] = None

# Concurrent workers for DEX queries.
N_WORKERS: int = 8

# ---------------------------------------------------------------------------
# Candidate generation
# ---------------------------------------------------------------------------

VOWELS = set("aăâeiîou")


def strip_gender_suffix(lemma: str) -> str:
    """Strip a single final vowel (gender/declension marker)."""
    if len(lemma) > 1 and lemma[-1] in VOWELS:
        return lemma[:-1]
    return lemma


def mutate_a_to_aa(stem: str) -> str | None:
    """Generate the 'a -> ă everywhere in stem' variant. Returns None if no
    'a' is present. Romanian unstressed-syllable mutation is lexical, so we
    just produce the variant and let DEX confirm or reject."""
    if "a" not in stem:
        return None
    return stem.replace("a", "ă")


def generate_candidates(lemma: str, stem_final: str) -> list[tuple[str, str]]:
    """Return a list of (candidate_verb, suffix_label) tuples to test against DEX.

    Conservative set: 2 stems x 3 suffixes x 2 prefixes = up to 12 candidates.
    """
    stem = strip_gender_suffix(lemma)
    stems = {stem}
    mut = mutate_a_to_aa(stem)
    if mut:
        stems.add(mut)

    # The three productive Romanian verbalizers
    SUFFIXES: list[tuple[str, str]] = [
        ("i", "-i"),
        ("a", "-a"),
        ("ui", "-ui"),
    ]

    # Either no prefix or in- (with morphophonemic îm-/în- alternation)
    PREFIXES: list[str] = ["", "în"]

    cands: set[tuple[str, str]] = set()
    for s in stems:
        if not s:
            continue
        for sfx, label in SUFFIXES:
            for pre in PREFIXES:
                cand = pre + s + sfx
                # in/im alternation: îm before labials, în elsewhere
                if (
                    cand.startswith("în")
                    and len(cand) > 2
                    and cand[2] in "bpm"
                ):
                    cand = "îm" + cand[2:]
                cands.add((cand, label))

    return [(c, l) for c, l in cands if len(c) >= 3]


_TITLE_RE = re.compile(
    r"<title>\s*([^|<-]+?)\s*-\s*defini[țt]ie", re.IGNORECASE
)
_VB_LEXEM_RE = re.compile(r"<span[^>]*>\s*(?:vb\.|verb)", re.IGNORECASE)
_NORM_TR = str.maketrans(
    {"â": "î"}
)  # DEX uses both â and î; normalize for compare


def _title_lemma(html: str) -> str | None:
    m = _TITLE_RE.search(html or "")
    if not m:
        return None
    return m.group(1).strip().lower()


def _fast_get_html(candidate: str) -> str | None:
    """Return cached HTML for `candidate` without BeautifulSoup parsing.
    Falls back to polite_get only if not in cache. Uses URL-encoded keys
    matching what polite_get would store (urllib.parse.quote)."""
    from urllib.parse import quote

    cand = norm_lemma(candidate)
    if not cand:
        return None
    cand_quoted = quote(cand)
    urls = [
        f"https://dexonline.ro/definitie/{cand_quoted}/definitii",
        f"https://dexonline.ro/definitie/{cand_quoted}",
        # Also try the unquoted form in case any cached entries used it
        f"https://dexonline.ro/definitie/{cand}/definitii",
        f"https://dexonline.ro/definitie/{cand}",
    ]
    for url in urls:
        if url in HTML_CACHE:
            return HTML_CACHE[url]
        if url in DISK_CACHE:
            HTML_CACHE[url] = DISK_CACHE[url]
            return DISK_CACHE[url]
    # Cache miss — fetch via polite_get (handles throttle, cache write, retries).
    try:
        resp = polite_get(urls[0])
        return resp.text
    except Exception:
        return None


def is_attested_verb(candidate: str) -> bool:
    """Return True iff DEX attests this exact form as a verb infinitive.
    Three checks (all must pass):
      (a) page exists and its canonical title lemma equals `candidate`
          (rules out conjugated-form redirects like ataci -> ataca title)
      (b) page contains a vb./verb POS marker
      (c) page contains the literal "a {candidate}" infinitive phrase, which
          distinguishes a real verb lemma from a page that merely shares its
          URL with verb content (e.g. câștigi as 2sg of câștiga: page has
          vb. content for câștiga but doesn't say "a câștigi").
    """
    html = _fast_get_html(candidate)
    if not html:
        return False
    if not _VB_LEXEM_RE.search(html):
        return False
    title_lemma = _title_lemma(html)
    if not title_lemma:
        return False
    cand_norm = norm_lemma(candidate)
    if title_lemma.translate(_NORM_TR) != cand_norm.translate(_NORM_TR):
        return False
    # Final: require the literal "a {cand}" infinitive phrase
    return f"a {cand_norm}" in html.lower()


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------


def main() -> None:
    if not WITH_FREQ_CSV.exists():
        print(f"ERROR: not found: {WITH_FREQ_CSV}", file=sys.stderr)
        sys.exit(1)
    if not RAW_DEX_CSV.exists():
        print(f"ERROR: not found: {RAW_DEX_CSV}", file=sys.stderr)
        sys.exit(1)

    # Read with_freq for stem_final and freq lookup, but use raw_dex to determine
    # eligibility (empty derived_verbs) since raw_dex is the upstream source of
    # truth and avoids stale state from previous supplement runs.
    print(f"Reading {WITH_FREQ_CSV} for stem_final/freq lookup ...")
    stem_by_lemma = {}
    freq_by_lemma = {}
    with WITH_FREQ_CSV.open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            lemma = (r.get("lemma") or "").strip().lower()
            pos = (r.get("pos") or "").strip().upper()
            stem_final = (r.get("stem_final") or "").strip()
            if pos != "N" or not lemma:
                continue
            stem_by_lemma[lemma] = stem_final
            try:
                freq_by_lemma[lemma] = float(
                    r.get("freq_ron_news_2024_1M") or 0
                )
            except ValueError:
                freq_by_lemma[lemma] = 0.0

    print(
        f"Reading {RAW_DEX_CSV} to find lemmas with currently-empty derived_verbs ..."
    )
    target_records = []
    seen = set()
    with RAW_DEX_CSV.open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            lemma = (r.get("lemma") or "").strip().lower()
            pos = (r.get("pos") or "").strip().upper()
            existing = (r.get("derived_verbs") or "").strip()
            if pos != "N" or not lemma:
                continue
            if lemma in seen:
                continue
            if existing:
                continue
            sf = stem_by_lemma.get(lemma, "")
            if sf not in TARGET_STEMS:
                continue
            seen.add(lemma)
            target_records.append((lemma, sf, freq_by_lemma.get(lemma, 0.0)))

    # Sort by frequency (high first) and apply top-N filter
    target_records.sort(key=lambda x: -x[2])
    if TOP_N_BY_FREQUENCY is not None:
        target_records = target_records[:TOP_N_BY_FREQUENCY]
        if target_records:
            print(
                f"  Restricted to top-{TOP_N_BY_FREQUENCY} by frequency "
                f"(min freq in set: {target_records[-1][2]:.0f})"
            )
    targets_by_lemma = {lemma: sf for lemma, sf, _ in target_records}
    print(f"  Lemmas eligible for supplementation: {len(targets_by_lemma):,}")

    # Read raw_dex (upstream of stage 4) — we'll write supplemented derived_verbs
    # here so the processor will recompute everything cleanly.
    print(f"\nReading {RAW_DEX_CSV} ...")
    with RAW_DEX_CSV.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    print(f"  Loaded {len(rows):,} raw_dex rows")

    # Build the row list to supplement: only those whose lemma appears in targets_by_lemma.
    targets = []
    for r in rows:
        lemma = (r.get("lemma") or "").strip().lower()
        if lemma in targets_by_lemma:
            r["__stem_final"] = targets_by_lemma[lemma]
            targets.append(r)
    print(f"  Matched in raw_dex: {len(targets):,}")
    print(f"\nNouns eligible for supplementation: {len(targets):,}")
    print(
        f"  (target stems = {sorted(TARGET_STEMS)}, currently empty derived_verbs)"
    )

    queries = 0
    confirmed = 0
    confirmed_lemmas = []
    progress_lock = threading.Lock()
    state = {"idx": 0}

    t0 = time.time()
    save_every = 500

    def process_row(row: dict[str, str]) -> None:
        nonlocal queries, confirmed
        lemma = row["lemma"].strip().lower()
        stem_final = row["__stem_final"]
        cands = generate_candidates(lemma, stem_final)
        kept_verbs: list[str] = []
        kept_suffixes: list[str] = []
        local_q = 0
        for cand, label in cands:
            local_q += 1
            if is_attested_verb(cand):
                kept_verbs.append(cand)
                kept_suffixes.append(label)
        with progress_lock:
            queries += local_q
            state["idx"] += 1
            idx = state["idx"]
            if kept_verbs:
                row["derived_verbs"] = " | ".join(kept_verbs)
                # Don't write deriv_suffixes here — that column is not in
                # raw_dex.csv (it's computed by stage 4 / romanian_processor_main.py
                # from derived_verbs, using its own canonical labelling).
                confirmed += 1
                confirmed_lemmas.append((lemma, row["derived_verbs"]))
            if idx % 25 == 0:
                elapsed = time.time() - t0
                rate = idx / elapsed if elapsed > 0 else 0
                est_total = (len(targets) / max(rate, 1e-6)) / 60
                print(
                    f"  [{idx:5}/{len(targets)}]  confirmed={confirmed:4d}  "
                    f"queries={queries:7d}  rate={rate:.2f} lemmas/s  "
                    f"elapsed={elapsed/60:.1f}min  eta_total={est_total:.0f}min",
                    flush=True,
                )
            if idx % save_every == 0:
                save_disk_cache()

    print(f"\nLaunching {N_WORKERS} worker threads ...")
    with ThreadPoolExecutor(max_workers=N_WORKERS) as pool:
        futures = [pool.submit(process_row, r) for r in targets]
        for fut in futures:
            try:
                fut.result()
            except Exception as e:
                print(f"  worker error: {e!r}", flush=True)

    save_disk_cache()
    elapsed = time.time() - t0

    print("\n=== Done ===")
    print(f"  Eligible nouns processed: {len(targets):,}")
    print(f"  Nouns supplemented with at least one verb: {confirmed:,}")
    print(f"  Total candidate queries: {queries:,}")
    print(f"  Elapsed: {elapsed/60:.1f} min", flush=True)

    # Show a sample of supplemented entries
    print("\n=== Sample supplemented (first 30) ===")
    for lemma, verbs in confirmed_lemmas[:30]:
        print(f"  {lemma:18}  derived_verbs = {verbs}")

    # Write output (raw_dex format — drop the temp __stem_final column we added)
    print(f"\nWriting {OUTPUT_CSV} ...")
    for r in rows:
        r.pop("__stem_final", None)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Wrote {len(rows):,} rows")
    print("\nNext steps to integrate into the pipeline:")
    print(
        f"  cp {OUTPUT_CSV} {RAW_DEX_CSV}   # backup the original first if you like"
    )
    print(
        "  ./run_pipeline.sh                # re-runs stages 4-6 with the new derived_verbs"
    )


if __name__ == "__main__":
    main()
