#!/usr/bin/env python3
"""Build the inflection-dependence contingency table the right way.

This is the canonical script for the table that goes into the abstract.
Replicates Steriade (2008) §10.2.3.5 lexical counts under the appropriate
coding for our data:

  ROW    — plural alternates yes/no  (using the `mutation` field, which is
           the direct sg-vs-pl surface comparison: t→ț, d→z, s→ș for
           coronals; c/g → orthographic ci/ge for dorsals).
  COLUMN — verbalizer allomorph  (-i / -a)  — Steriade's affix-avoidance
           tabulation in her (10.20).
  EXCLUSIONS — derived verbs whose DEX page does not pass etymology
           validation (verb tree's lemma matches the candidate; etymology
           cites our noun, an "În + our noun" form, or a foreign source).

  USER-OVERRIDE KEEPS — back-formations the etym detector mis-rejects:
           clint, dinte, vargă, verigă.

Key methodological points:

  1. The row is coded by ACTUAL alternation (`mutation`), not by the
     surface allomorph (-i/-e vs -uri). Coronal -e plurals don't
     palatalize, so they belong in the non-alternating row even though
     their suffix is phonologically "front". This is the row recode
     that fixed the original measurement error.

  2. The column is NOT filtered by whether the conjugated verb
     independently surfaces the alternation. Augmented verbs (plăti,
     bloca, etc.) are KEPT even though their paradigms morphologically
     bleed the trigger. This matches Steriade's own tabulation in
     (10.20), which is an affix-avoidance test — what allomorph did the
     speaker select?  Not — does the verb's paradigm show palatalization?

     The orthogonal question of augment-driven trigger absence is
     reported separately by scripts/audit_verb_alternation.py.

Outputs two tables:
  - TOP-1000 most-frequent noun lemmas
  - FULL DEX (no frequency filter)

Each split into three panels: dorsal (c, g), coronal (t, d, s), pooled.
"""
from __future__ import annotations

import csv
import ijson
import re
from collections import Counter
from pathlib import Path
from urllib.parse import quote

from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parent.parent
CSV = ROOT / "data" / "romanian_lexicon_with_freq.csv"
CACHE = ROOT / "dex_cache.json"

PLURAL_FRONT = {"i", "e"}
SEGS_DORSAL = {"c", "g"}
SEGS_CORONAL = {"t", "d", "s"}

# Back-formations whose DEX etymology field is "necunoscută" but which
# we keep as valid noun-verb pairs because synchronically they look like
# typical denominal pairs (the noun and verb share the same root, even
# if the historical direction is verb→noun).
USER_KEEP = {"clint", "dinte", "vargă", "verigă"}


# ===========================================================================
# DEX etymology validation
# ===========================================================================

VERB_POS_RE = re.compile(
    r'<span[^>]*class="[^"]*tree-pos-info[^"]*"[^>]*>verb', re.IGNORECASE)
TONIC_RE = re.compile(
    r'<span[^>]*class="[^"]*tonic-accent[^"]*"[^>]*>(.*?)</span>',
    re.IGNORECASE | re.DOTALL)
INFLECTED_RE = re.compile(
    r'<span[^>]*class="[^"]*tree-inflected-form[^"]*"[^>]*>.*?</span>',
    re.IGNORECASE | re.DOTALL)
POSINFO_RE = re.compile(
    r'<span[^>]*class="[^"]*tree-pos-info[^"]*"[^>]*>.*?</span>',
    re.IGNORECASE | re.DOTALL)
ANY_TAG = re.compile(r'<[^>]+>')
ETYM_RE = re.compile(
    r'etimologie\s*[:\s]\s*(.{1,200}?)'
    r'(?=\s+(?:DEX|MDA|DLRLC|DLRM|DEXI|CADE|NODEX|info|sinonime|antonime|format|$))',
    re.IGNORECASE | re.DOTALL)
FOREIGN_RE = re.compile(
    r'\b(?:limba\s+)?(?:franceză|latină|slav[ăa]|veche|neogreacă|greacă|'
    r'germană|italiană|maghiară|engleză|turcă|spaniolă|rusă|bulgară|sârbă|'
    r'polonă|ucrainean|săsească|albaneză|portugheză|cehă|olandeză|catalană|'
    r'lat\.|fr\.|gr\.|sl\.|magh\.|germ\.|it\.|engl\.|tc\.|sb\.|bg\.|rus\.|alb\.)',
    re.IGNORECASE)
UNKNOWN_RE = re.compile(
    r'\bnecunoscut|et\.\s*nec\.|etimol\.\s*nec\.', re.IGNORECASE)


def heading_lemma(h: str) -> str | None:
    h = TONIC_RE.sub(r"\1", h)
    h = INFLECTED_RE.sub("", h)
    h = POSINFO_RE.sub("", h)
    h = ANY_TAG.sub("", h)
    m = re.match(r"([^,;\s]+)", re.sub(r"\s+", " ", h).strip())
    return m.group(1).strip().lower() if m else None


def get_etym_for_verb(html: str, candidate: str) -> str | None:
    """Return the etymology field of the DEX verb tree whose lemma matches
    `candidate`, or None if no such tree exists."""
    h3 = list(re.finditer(
        r'<h3[^>]*class="[^"]*tree-heading[^"]*"[^>]*>(.*?)</h3>',
        html, re.DOTALL | re.IGNORECASE))
    for i, m in enumerate(h3):
        if not VERB_POS_RE.search(m.group(1)):
            continue
        if heading_lemma(m.group(1)) != candidate.lower():
            continue
        body = html[m.end(): h3[i + 1].start() if i + 1 < len(h3) else len(html)]
        text = re.sub(r"\s+", " ", ANY_TAG.sub(" ", body)).strip()
        ets = ETYM_RE.findall(text)
        return ets[0].strip() if ets else ""
    return None


def stem_variants(noun: str) -> set[str]:
    b = noun.lower()
    stems = {b, b[:-1]} if (len(b) > 1 and b[-1] in "ăaeio") else {b}
    out = set(stems)
    for s in stems:
        if "a" in s: out.add(s.replace("a", "ă"))
        if "ă" in s: out.add(s.replace("ă", "a"))
        if "o" in s: out.add(s[:s.find("o")] + "oa" + s[s.find("o") + 1:])
        if "oa" in s: out.add(s.replace("oa", "o"))
        if "e" in s: out.add(s[:s.find("e")] + "ea" + s[s.find("e") + 1:])
        if "ea" in s: out.add(s.replace("ea", "e"))
    return out


def categorize(noun: str, etym: str | None) -> str:
    """Return A/B/C/D/E/F per the etymology category schema:
      A — cites our noun
      B — cites 'În + our noun' (prefixed denominal)
      C — cites a different Romanian word
      D — cites a foreign source (synchronic pair)
      E — necunoscută (DEX disclaims any source)
      F — no etymology field on page or no verb tree
    """
    if etym is None or etym == "":
        return "F"
    et = etym.lower()
    nv = stem_variants(noun)
    if UNKNOWN_RE.search(et):
        return "E"
    if FOREIGN_RE.search(et):
        return "D"
    m = re.search(r'\bvezi\s+([a-zA-ZăâîșțĂÂÎȘȚ]+)', et)
    if m and any(
        m.group(1).lower() == v or m.group(1).lower().startswith(v)
        or v.startswith(m.group(1).lower()) for v in nv
    ):
        return "A"
    m = re.search(r'\b(?:în|îm|in|im)\s*\+\s*([a-zA-ZăâîșțĂÂÎȘȚ]+)', et)
    if m:
        c = m.group(1).lower()
        return "B" if any(c == v or c.startswith(v) or v.startswith(c)
                          for v in nv) else "C"
    m = re.search(r'\bA\s*\+\s*([a-zA-ZăâîșțĂÂÎȘȚ]+)', etym)
    if m:
        c = m.group(1).lower()
        return "B" if any(c == v or c.startswith(v) or v.startswith(c)
                          for v in nv) else "C"
    m = re.search(r'\b([a-zA-ZăâîșțĂÂÎȘȚ]{2,})', et)
    if m:
        c = m.group(1).lower()
        if c in {"cf", "vezi", "din", "de", "la", "se", "și", "sau"}:
            m2 = re.search(
                r'\b[a-zăâîșț]{2,}\s+([a-zA-ZăâîșțĂÂÎȘȚ]{2,})', et)
            if m2:
                c = m2.group(1).lower()
        return "A" if any(c == v or c.startswith(v) or v.startswith(c)
                          for v in nv) else "C"
    return "F"


def url_variants(form: str) -> list[str]:
    q = quote(form)
    return [
        f"https://dexonline.ro/definitie/{q}/definitii",
        f"https://dexonline.ro/definitie/{q}",
        f"https://dexonline.ro/definitie/{form}/definitii",
        f"https://dexonline.ro/definitie/{form}",
    ]


def load_cache_subset(verbs: set[str]) -> dict[str, str]:
    needed = set()
    for v in verbs:
        needed.update(url_variants(v))
    cache: dict[str, str] = {}
    with CACHE.open("rb") as f:
        for k, v in ijson.kvitems(f, ""):
            if k in needed:
                cache[k] = v
    return cache


def get_html(cache: dict[str, str], form: str) -> str | None:
    for u in url_variants(form):
        if u in cache:
            return cache[u]
    return None


def keep_pair(noun: str, dvs: list[str], cache: dict[str, str]) -> bool:
    """Etymology validation: at least one of the derived verbs must
    categorize as A, B, or D."""
    if noun in USER_KEEP:
        return True
    cats = set()
    for dv in dvs:
        html = get_html(cache, dv)
        if html is None:
            cats.add("?")
            continue
        cats.add(categorize(noun, get_etym_for_verb(html, dv)))
    return bool(cats & {"A", "B", "D"})


# ===========================================================================
# Build the table
# ===========================================================================

def load_data() -> tuple[list[tuple], dict[str, float]]:
    """Load all (lemma, dvs, stem_final, opportunity, deriv_suffix, mutation,
    freq) tuples that survive the basic row filters."""
    work = []
    seen = set()
    freq = {}
    with CSV.open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["pos"] != "N":
                continue
            try:
                freq[r["lemma"]] = float(r.get("freq_ron_news_2024_1M") or 0)
            except ValueError:
                freq[r["lemma"]] = 0.0
            sf = r["stem_final"]
            if sf not in SEGS_DORSAL | SEGS_CORONAL:
                continue
            mut_str = r["mutation"].strip()
            if mut_str == "":
                continue
            opp = r["opportunity"].strip()
            if opp not in PLURAL_FRONT and opp != "uri":
                continue
            ds = r["deriv_suffixes"].strip()
            if ds not in ("-i", "-a"):
                continue
            dv_str = (r.get("derived_verbs") or "").strip()
            if not dv_str:
                continue
            if (r.get("exception_reason") or "").startswith("nde:"):
                continue
            if r["lemma"] in seen:
                continue
            seen.add(r["lemma"])
            mutation = mut_str.upper() == "TRUE"
            dvs = [d.strip() for d in dv_str.split("|") if d.strip()]
            work.append((r["lemma"], dvs, sf, opp, ds, mutation,
                         freq[r["lemma"]]))
    return work, freq


def fisher_panel(name: str, cells: Counter) -> str:
    fi = cells[("alt", "i")]
    fa = cells[("alt", "a")]
    bi = cells[("nonalt", "i")]
    ba = cells[("nonalt", "a")]
    n = fi + fa + bi + ba
    out = (f"  {name:8s}:  alt-pl  {fi:3d}  -i / {fa:3d}  -a    "
           f"nonalt-pl  {bi:3d}  -i / {ba:3d}  -a    n = {n:4d}")
    if (fi + fa) > 0 and (bi + ba) > 0 and (fi + bi) > 0 and (fa + ba) > 0:
        odds, p = fisher_exact([[fi, fa], [bi, ba]])
        out += f"    OR = {odds:.2f}    p = {p:.3f}"
    return out


def build(work: list[tuple], cache: dict[str, str],
          freq_filter) -> tuple[dict[str, Counter], int]:
    """Build the three panels (dorsal / coronal / pooled) over the rows
    that pass etymology validation and the frequency filter."""
    cells = {"dorsal": Counter(), "coronal": Counter(),
             "pooled": Counter()}
    n_kept = 0
    for noun, dvs, sf, opp, ds, mutation, fr in work:
        if not freq_filter(noun, fr):
            continue
        if not keep_pair(noun, dvs, cache):
            continue
        n_kept += 1
        seg = "dorsal" if sf in SEGS_DORSAL else "coronal"
        row = "alt" if mutation else "nonalt"
        col = "i" if ds == "-i" else "a"
        cells[seg][(row, col)] += 1
        cells["pooled"][(row, col)] += 1
    return cells, n_kept


def main():
    work, freq = load_data()
    print(f"Loaded {len(work)} (lemma, dvs, ...) tuples that pass row filters.")
    all_dvs = {dv for _, dvs, *_ in work for dv in dvs}
    cache = load_cache_subset(all_dvs)
    print(f"Cached DEX pages: {len(cache)} (for {len(all_dvs)} unique verbs)")

    ranked = sorted(freq.items(), key=lambda x: -x[1])
    top1000 = {lem for lem, _ in ranked[:1000]}

    scopes = [
        ("TOP-1000", lambda n, fr: n in top1000),
        ("FULL DEX", lambda n, fr: True),
    ]

    print("\n" + "=" * 78)
    print("Contingency table: row = plural alternates (mutation),"
          " column = verbalizer allomorph")
    print("=" * 78)
    for label, fn in scopes:
        cells, n = build(work, cache, fn)
        print(f"\n{label}  (n_kept = {n})")
        for s in ("dorsal", "coronal", "pooled"):
            print(fisher_panel(s, cells[s]))


if __name__ == "__main__":
    main()
