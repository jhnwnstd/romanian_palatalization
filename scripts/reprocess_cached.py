#!/usr/bin/env python3
"""Re-parse all cached wikitext using the improved harvester logic.

Runs entirely from cache — zero HTTP calls. Uses the existing
wikitext and HTML caches to re-extract plurals, IPA, derivations,
glosses, and gender with all the latest parsing improvements.

Usage:
    python scripts/reprocess_cached.py
"""

import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

import romanian_harvester as rh  # noqa: E402
from romanian_harvester import (  # noqa: E402
    CANARY_LEMMAS,
    DENOMINAL_ADJS,
    DENOMINAL_VERBS,
    OUTPUT_CSV,
    POS_HEAD_RE,
    POS_MAP,
    ROMANIAN_SECTION_RE,
    STERIADE_EXAMPLES,
    extract_derived_terms_from_entry,
    has_verb_pos,
    is_candidate_title,
    load_disk_cache,
    load_ipa_cache,
    normalize_unicode,
    parse_romanian_entry,
    register_adj_derivations,
    register_verb_derivations,
)

# Disable network calls — cache only
rh.ENABLE_IPA_FETCH = False


def main():
    print("=" * 70)
    print("REPROCESS FROM CACHE (no HTTP calls)")
    print("=" * 70)

    load_disk_cache()
    load_ipa_cache()

    # Collect all titles from both caches (use module-level refs)
    all_titles_set = set(rh._WT_CACHE_EN.keys()) | set(rh._WT_CACHE_RO.keys())
    all_titles = [t for t in all_titles_set if is_candidate_title(t)]
    # Add canary/Steriade examples that might be in cache
    for t in CANARY_LEMMAS + STERIADE_EXAMPLES:
        if t not in all_titles_set and (
            t in rh._WT_CACHE_EN or t in rh._WT_CACHE_RO
        ):
            all_titles.append(t)

    print(f"Cached titles: {len(all_titles)}")

    # Pass 1: Scan verbs for derivations
    print("\nScanning verbs for derivations...")
    verb_count = 0
    for title in all_titles:
        title_norm = normalize_unicode(title)
        wt = rh._WT_CACHE_EN.get(title_norm) or rh._WT_CACHE_RO.get(
            title_norm, ""
        )
        if not wt:
            continue
        m = ROMANIAN_SECTION_RE.search(wt)
        if not m:
            continue
        ro = m.group(1)
        if not has_verb_pos(ro):
            continue
        verb_count += 1
        register_verb_derivations(title_norm, ro)

    print(
        f"  {verb_count} verbs scanned, "
        f"{len(DENOMINAL_VERBS)} denominal bases"
    )

    # Pass 2: Scan adjectives for derivations
    print("Scanning adjectives for derivations...")
    adj_count = 0
    for title in all_titles:
        title_norm = normalize_unicode(title)
        wt = rh._WT_CACHE_EN.get(title_norm) or rh._WT_CACHE_RO.get(
            title_norm, ""
        )
        if not wt:
            continue
        m = ROMANIAN_SECTION_RE.search(wt)
        if not m:
            continue
        ro = m.group(1)
        has_adj = False
        for pos_m in POS_HEAD_RE.finditer(ro):
            if POS_MAP.get(pos_m.group(1).lower()) == "ADJ":
                has_adj = True
                break
        if not has_adj:
            continue
        adj_count += 1
        register_adj_derivations(title_norm, ro)

    print(
        f"  {adj_count} adjectives scanned, "
        f"{len(DENOMINAL_ADJS)} denominal adj bases"
    )

    # Pass 3: Extract derived terms from noun entries
    print("Extracting derived terms from entries...")
    for title in all_titles:
        title_norm = normalize_unicode(title)
        wt = rh._WT_CACHE_EN.get(title_norm) or rh._WT_CACHE_RO.get(
            title_norm, ""
        )
        if not wt:
            continue
        m = ROMANIAN_SECTION_RE.search(wt)
        if not m:
            continue
        extract_derived_terms_from_entry(title_norm, m.group(1))

    total_v = sum(len(v) for v in DENOMINAL_VERBS.values())
    total_a = sum(len(v) for v in DENOMINAL_ADJS.values())
    print(f"  Total: {total_v} N→V pairs, {total_a} N→Adj pairs")

    # Pass 4: Parse entries (cache only, skip_ipa for uncached)
    print(f"\nParsing {len(all_titles)} entries (cache only)...")
    all_rows = []
    seen = set()
    for i, title in enumerate(all_titles, 1):
        if title in seen:
            continue
        if i % 1000 == 0:
            print(f"  {i}/{len(all_titles)}...")
        entry = parse_romanian_entry(title, skip_ipa=False)
        if entry and entry["pos"] != "V":
            all_rows.append(entry)
            seen.add(title)

    print(f"\n{len(all_rows)} entries parsed")

    df = pd.DataFrame(all_rows)

    # Dedup
    df["source_rank"] = df["source"].apply(
        lambda s: 0 if "ro.wiktionary" in str(s) else 1
    )
    df["has_gloss"] = df["gloss"].notna() & (df["gloss"].str.strip() != "")
    df["has_plural"] = df["plural"].notna() & (df["plural"].str.strip() != "")
    df = df.sort_values(
        by=[
            "lemma",
            "pos",
            "source_rank",
            "has_gloss",
            "has_plural",
        ],
        ascending=[True, True, True, False, False],
    )
    before = len(df)
    df = df.drop_duplicates(subset=["lemma", "pos"], keep="first")
    df = df.drop(columns=["source_rank", "has_gloss", "has_plural"])
    print(f"Deduped: {before} → {len(df)}")

    # Write to the standard output path
    # (safe: the background harvester writes to a temp file first)
    print(f"Writing to {OUTPUT_CSV}...")
    df.to_csv(OUTPUT_CSV, index=False)

    print("\n" + "=" * 70)
    print("STATISTICS")
    print("=" * 70)
    print(f"Total: {len(df)}")
    print(f"Nouns: {len(df[df['pos'] == 'N'])}")
    print(f"Adjectives: {len(df[df['pos'] == 'ADJ'])}")
    print(f"With plural: {(df['plural'] != '').sum()}")
    print(f"With IPA: {(df['ipa_raw_lemma'] != '').sum()}")
    dv = df[df["derived_verbs"] != ""]
    da = df[df["derived_adj"] != ""]
    print(f"With derived verbs: {len(dv)}")
    print(f"With derived adj: {len(da)}")
    print("=" * 70)

    print("Done!")


if __name__ == "__main__":
    main()
