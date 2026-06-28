#!/usr/bin/env python3
"""Diagnostic audit of derived-verb conjugational morphology.

For each derived verb in the inflection-dependence contingency table:

  1. Pull its DEX conjugation table from the cache (no new fetches).
  2. Classify its present-paradigm augment status (paradigm-wide).
  3. Check whether the stem-final consonant surfaces in its palatalized
     form in any diagnostic cell (2sg/3sg pres indicative, 3sg pres
     subjunctive).
  4. Emit a three-way verdict per verb:
        ALTERNATES     — surface form shows the alternation
        NO-ALTERNATION — diagnostic cells exist but no alternation surfaces
        TRIGGER-ABSENT — augment/conjugation bleeds the trigger entirely

**Purpose**: characterize Romanian morphology of these verbs. The output
is DIAGNOSTIC, not CORRECTIVE. The verdicts here do NOT justify excluding
verbs from the affix-avoidance contingency table.

**Why**: Steriade's lexical counts in (10.20) tabulate verbalizer
allomorph distribution across all K/TS and K/K bases without filtering
by whether the verb's conjugated paradigm independently shows the
alternation. She is explicit (p. 9, after 10.15) that *"only the former
[affix avoidance] arises with derived verbs"* — the phonotactic-blocking
manifestation applies to single-form suffixes like -ist, not to verbs.

So a verb like plăti (augmented, TRIGGER-ABSENT) is still valid data
for the affix-avoidance test: plată chose -i as its verbalizer, and
that choice is what Steriade tabulates.

For the actual contingency table, see scripts/build_contingency_table.py.
"""
from __future__ import annotations

import csv
import ijson
import re
import sys
from collections import Counter
from pathlib import Path
from urllib.parse import quote

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

ROOT = Path(__file__).resolve().parent.parent
CSV_IN = ROOT / "data" / "romanian_lexicon_with_freq.csv"
CACHE = ROOT / "dex_cache.json"
OUT_TXT = ROOT / "analysis" / "verb_alternation_audit.txt"

# ---------------------------------------------------------------------------
# Conjugation-table parser
# ---------------------------------------------------------------------------

LEXEME_TABLE_RE = re.compile(
    r'<table[^>]*class="[^"]*lexeme[^"]*"[^>]*>(.*?)</table>',
    re.DOTALL | re.IGNORECASE,
)
PERSON_ROW_RE = re.compile(
    r'<td[^>]*class="[^"]*inflection person[^"]*"[^>]*>([^<]+)</td>',
    re.IGNORECASE,
)
FORM_CELL_RE = re.compile(
    r'<td[^>]*class="[^"]*form[^"]*"[^>]*>(.*?)</td>',
    re.DOTALL | re.IGNORECASE,
)
TONIC_RE = re.compile(
    r'<span[^>]*class="[^"]*tonic-accent[^"]*"[^>]*>([^<]+)</span>',
    re.IGNORECASE,
)
ANY_TAG = re.compile(r'<[^>]+>')


def clean_form(td_content: str) -> str:
    """Strip HTML cruft from a <td class="form"> cell and return its text."""
    x = TONIC_RE.sub(r'\1', td_content)
    x = ANY_TAG.sub(' ', x)
    x = x.replace('&#x2011', '').replace('‑', '')  # non-breaking hyphen
    # The subjunctive cells are prefixed with "(să) ". Strip it; the form
    # itself is what we want.
    x = re.sub(r'^\s*\(s[ăa]\)\s*', '', x)
    # Some cells have multiple forms separated by line breaks; take the first.
    parts = [p.strip() for p in re.split(r'\s{2,}|\n', x) if p.strip()]
    return parts[0] if parts else x.strip()


def all_forms(td_content: str) -> list[str]:
    """Return all surface forms in a form cell (handles commaList variants)."""
    x = TONIC_RE.sub(r'\1', td_content)
    # Each form sits inside a <li>...</li>
    lis = re.findall(r'<li[^>]*>(.*?)</li>', x, re.DOTALL)
    if not lis:
        return [clean_form(td_content)]
    out = []
    for li in lis:
        s = ANY_TAG.sub('', li).replace('&#x2011', '').replace('‑', '').strip()
        if s:
            out.append(s)
    return out


def extract_present_paradigm(html: str) -> dict[str, dict[str, list[str]]]:
    """Return {person: {column: [forms]}} for the first lexeme table.

    Columns in the present-tense section, in order:
        prez_ind, prez_subj, imperfect, perfect_simplu, mmcp.
    """
    tables = LEXEME_TABLE_RE.findall(html)
    if not tables:
        return {}
    table = tables[0]
    # Split on person markers to get one chunk per person row.
    # We find every person <td> and slice up to the next one.
    persons = list(PERSON_ROW_RE.finditer(table))
    if not persons:
        return {}
    out: dict[str, dict[str, list[str]]] = {}
    columns = ["prez_ind", "prez_subj", "imperfect", "perf_simplu", "mmcp"]
    for i, m in enumerate(persons):
        person = m.group(1).strip()
        end = persons[i + 1].start() if i + 1 < len(persons) else len(table)
        row_html = table[m.end():end]
        # Get the form cells in order.
        cells = FORM_CELL_RE.findall(row_html)
        row: dict[str, list[str]] = {}
        for col_i, col_name in enumerate(columns):
            if col_i < len(cells):
                row[col_name] = all_forms(cells[col_i])
        out[person] = row
    return out


def get_conjugare_class(html: str) -> str:
    """Return the conjugarea group tag (e.g. 'I', 'IV') if present."""
    # DEX tags verbs with "<a ...>conjugarea a IV-a</a>" or similar.
    m = re.search(r'conjugarea a ([IVX]+)-a', html)
    if m:
        return m.group(1)
    return ""


# ---------------------------------------------------------------------------
# Augment classifier (paradigm-wide)
# ---------------------------------------------------------------------------

# Augment patterns. The augment sits between the stem and the inflectional
# ending, inserting a vowel that bleeds the front-vowel trigger.
#
# 4th-conj augment: -esc, -ești, -ește. Appears in 1sg/2sg/3sg pres ind
# (and 3pl pres ind, which is identical to 1sg).
#
# 1st-conj augment: -ez, -ezi, -ează. Same positions.
AUGMENT_RES = [
    re.compile(r"(?:esc|ești|ește)\b", re.IGNORECASE),  # -esc class
    re.compile(r"(?:ez|ezi|ează)\b", re.IGNORECASE),    # -ez class
]


def has_augment(paradigm: dict) -> bool:
    """Return True if any singular-present cell shows an augment ending.

    Looks at the WHOLE present indicative (singular and plural) plus the
    present subjunctive — not just 1sg — to avoid misclassifying verbs
    whose 1sg is syncretic or irregular.
    """
    diagnostic_persons = ["I (eu)", "a II-a (tu)", "a III-a (el, ea)",
                          "a III-a (ei, ele)"]
    for person, row in paradigm.items():
        if person not in diagnostic_persons:
            continue
        for col in ("prez_ind", "prez_subj"):
            for form in row.get(col, []):
                for rgx in AUGMENT_RES:
                    if rgx.search(form):
                        return True
    return False


# ---------------------------------------------------------------------------
# Alternation detector
# ---------------------------------------------------------------------------

# Palatalized counterparts of stem-final consonants in Romanian orthography.
# For dorsals, the palatalized form is signalled by the FOLLOWING vowel:
# c+i, c+e, g+i, g+e are pronounced [tʃi tʃe dʒi dʒe]. So for dorsals the
# diagnostic is "is the stem-final c/g immediately followed by e or i in the
# verb stem?" rather than a consonant substitution.
#
# For coronals, the alternation is segmental: t→ț, d→z, s→ș. We look for
# the palatalized counterpart appearing in the position the stem-final
# occupies.

PAL_CORONAL = {"t": "ț", "d": "z", "s": "ș"}

FRONT_VOWELS = "ei"

# Known inflectional endings to strip when locating the stem-final position
# in a conjugated form. Listed longest-first so longer ones strip preferentially.
INFLECTION_ENDINGS = [
    "ească", "ești", "ește", "ească",
    "ează", "ezi", "eze", "esc",
    "ești", "ește", "ești",
    "ăm", "ați", "âm", "ind",
    "are", "ere", "ire",
    "ăsc",
    "i", "e", "a", "ă", "u",
]


def strip_ending(form: str) -> str:
    """Strip the longest known inflectional ending and return the stem prefix.

    This isolates the stem so we can check whether the stem-final consonant
    is the palatalized or unpalatalized counterpart at the *correct position*
    — not anywhere in the form, which would over-detect palatalization from
    suffix material like -ești (whose ș is the augment's own).
    """
    for end in INFLECTION_ENDINGS:
        if form.endswith(end) and len(form) > len(end):
            return form[: -len(end)]
    return form


def stem_alternates_in_form(noun_stem: str, verb_stem_final: str,
                             form: str) -> bool:
    """Does this surface form show palatalization of the noun's stem-final
    consonant *at the stem-final position*?

    We isolate the stem by stripping a known inflectional ending, then check
    whether the last consonant of that stem is the palatalized counterpart
    (for coronals: t→ț, d→z, s→ș) or sits in a palatalizing position (for
    dorsals: c/g + front vowel inside the stem itself).
    """
    if not form:
        return False
    form_l = form.lower()
    sf = verb_stem_final.lower()
    stem = strip_ending(form_l)
    if not stem:
        return False
    if sf in PAL_CORONAL:
        return stem.endswith(PAL_CORONAL[sf])
    if sf in ("c", "g"):
        # Dorsal: stem-final orthographically c/g is palatalized iff followed
        # by a front vowel orthographically *inside the stem material* — but
        # in many dorsal forms the stem-final ends the stem and the front-
        # vowel ending follows (e.g. fabrici = fabric+i, where c+i is the
        # palatalized pronunciation [tʃi]). So we accept either:
        #   (a) the stem itself ends in c/g and the ORIGINAL FULL FORM has
        #       that c/g immediately followed by a front vowel in the
        #       stripped position, OR
        #   (b) somewhere inside the stem, c/g is followed by a front vowel.
        # Implementation: check if c/g is immediately followed by e or i
        # at the boundary between stem and ending OR inside the stem.
        # Combine by checking the original form for c/g + front-V where the
        # position aligns with the end of the stem.
        if stem.endswith(sf):
            # stem ends in c/g; check the original form: does the c/g sit
            # right before a front vowel?
            pos = len(stem) - 1
            if pos + 1 < len(form_l) and form_l[pos + 1] in FRONT_VOWELS:
                return True
        # Also accept palatalization triggered inside the stem itself
        for i, ch in enumerate(stem[:-1]):
            if ch == sf and stem[i + 1] in FRONT_VOWELS:
                return True
        return False
    return False


def diagnostic_cells(paradigm: dict) -> list[str]:
    """The cells where assibilation/palatalization is most diagnostic:
    2sg/3sg present indicative and 3sg present subjunctive (Steriade's
    "să sece" form).
    """
    out = []
    for person in ("a II-a (tu)", "a III-a (el, ea)"):
        row = paradigm.get(person, {})
        for form in row.get("prez_ind", []):
            out.append(form)
        for form in row.get("prez_subj", []):
            out.append(form)
    # Also include 1sg pres ind for completeness (it shows the augment).
    row = paradigm.get("I (eu)", {})
    for form in row.get("prez_ind", []):
        out.append(form)
    return out


def trigger_present(stem_final: str, paradigm: dict) -> bool:
    """Does any diagnostic cell put the stem-final consonant in a position
    where the alternation COULD surface?

    The diagnostic position is the stem-final position in the conjugated
    form, isolated by stripping the inflectional ending. The trigger is
    present if, *at that stem-final position*, the consonant is followed
    by a front vowel (the form's stem-final immediately precedes an
    inflectional vowel-initial ending).
    """
    sf = stem_final.lower()
    for form in diagnostic_cells(paradigm):
        f = form.lower()
        stem = strip_ending(f)
        if not stem:
            continue
        end = f[len(stem):]
        # The trigger is present iff the stem ends in our target consonant
        # (or its palatalized counterpart, meaning it already alternated)
        # AND the ending starts with a front vowel.
        if sf in PAL_CORONAL:
            pal = PAL_CORONAL[sf]
            if stem.endswith(sf) or stem.endswith(pal):
                if end and end[0] == "i":
                    return True
        elif sf in ("c", "g"):
            if stem.endswith(sf):
                if end and end[0] in FRONT_VOWELS:
                    return True
            # Or palatalization already realized inside the stem
            for i, ch in enumerate(stem[:-1]):
                if ch == sf and stem[i + 1] in FRONT_VOWELS:
                    return True
    return False


# ---------------------------------------------------------------------------
# Main audit
# ---------------------------------------------------------------------------

def url_variants(form: str) -> list[str]:
    q = quote(form)
    return [
        f"https://dexonline.ro/definitie/{q}/definitii",
        f"https://dexonline.ro/definitie/{q}",
        f"https://dexonline.ro/definitie/{form}/definitii",
        f"https://dexonline.ro/definitie/{form}",
    ]


def classify_verb(html: str, stem_final: str) -> dict:
    """Return a dict with augment status, alternation status, and verdict."""
    paradigm = extract_present_paradigm(html)
    if not paradigm:
        return {"verdict": "NO-PARADIGM", "augment": None, "alt": None,
                "paradigm": None}
    aug = has_augment(paradigm)
    # Did any cell actually show alternation?
    alt = any(stem_alternates_in_form("", stem_final, f)
              for f in diagnostic_cells(paradigm))
    trig = trigger_present(stem_final, paradigm)
    if alt:
        verdict = "ALTERNATES"
    elif not trig:
        verdict = "TRIGGER-ABSENT"
    elif aug:
        verdict = "TRIGGER-ABSENT"  # augment blocks even if some cell looks trigger-like
    else:
        verdict = "NO-ALTERNATION"
    return {
        "verdict": verdict,
        "augment": aug,
        "alt": alt,
        "trigger": trig,
        "conjugare": get_conjugare_class(html),
        "paradigm": paradigm,
    }


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


def main():
    # ----- validation pass first -----
    print("=" * 70, file=sys.stderr)
    print("VALIDATION", file=sys.stderr)
    print("=" * 70, file=sys.stderr)
    print("Loading validation anchors from cache...", file=sys.stderr)
    val_cache = load_cache_subset({"plăti", "lupta", "munci", "fabrica",
                                    "ataca", "scuza"})
    for verb, stem_final, expected in [
        ("plăti",   "t", "TRIGGER-ABSENT"),    # augmented -esc, no alternation
        ("lupta",   "t", "ALTERNATES"),         # unaugmented -a, 2sg lupți
        ("munci",   "c", "ALTERNATES"),         # dorsal -i
        ("fabrica", "c", "ALTERNATES"),         # dorsal -a (1sg fabric, 2sg fabrici)
        ("ataca",   "c", "ALTERNATES"),         # dorsal -a
        ("scuza",   "z", "TRIGGER-ABSENT"),     # has -ez augment; stem already z
    ]:
        html = get_html(val_cache, verb)
        if html is None:
            print(f"  {verb:14s}: NOT IN CACHE", file=sys.stderr)
            continue
        r = classify_verb(html, stem_final)
        status = "OK" if r["verdict"] == expected else "❌"
        print(f"  {status} {verb:14s} stem={stem_final}  verdict={r['verdict']:14s}  "
              f"expected={expected}  conjug={r['conjugare']}  augment={r['augment']}",
              file=sys.stderr)
        if r.get("paradigm"):
            sg_pres = r["paradigm"].get("a II-a (tu)", {}).get("prez_ind", [])
            sg3_pres = r["paradigm"].get("a III-a (el, ea)", {}).get("prez_ind", [])
            print(f"      2sg.pres={sg_pres}  3sg.pres={sg3_pres}", file=sys.stderr)

    # ----- full audit pass -----
    print("\n" + "=" * 70, file=sys.stderr)
    print("FULL AUDIT", file=sys.stderr)
    print("=" * 70, file=sys.stderr)

    # Load every (noun, derived_verb, stem_final) from the with_freq.csv that
    # would land in the contingency table (the etym-validated 388-row pool).
    # For this audit we don't re-run etym validation; we just read all rows
    # that match the basic filters.
    rows = []
    with CSV_IN.open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["pos"] != "N": continue
            sf = r["stem_final"]
            if sf not in {"c","g","t","d","s"}: continue
            if r["mutation"].strip() == "": continue
            opp = r["opportunity"].strip()
            if opp not in {"i","e","uri"}: continue
            ds = r["deriv_suffixes"].strip()
            if ds not in ("-i","-a"): continue
            dv_str = (r.get("derived_verbs") or "").strip()
            if not dv_str: continue
            if (r.get("exception_reason") or "").startswith("nde:"): continue
            verbs = [d.strip() for d in dv_str.split("|") if d.strip()]
            for v in verbs:
                rows.append((r["lemma"], v, sf, opp, ds))

    print(f"  (noun, verb, sf) tuples: {len(rows)}", file=sys.stderr)
    all_verbs = {v for _, v, *_ in rows}
    print(f"  Unique verbs to look up: {len(all_verbs)}", file=sys.stderr)
    cache = load_cache_subset(all_verbs)
    print(f"  Pages matched in cache: {len(cache)}", file=sys.stderr)

    results = []
    for noun, verb, sf, opp, ds in rows:
        html = get_html(cache, verb)
        if html is None:
            results.append((noun, verb, sf, opp, ds, {"verdict": "NO-PAGE"}))
            continue
        r = classify_verb(html, sf)
        results.append((noun, verb, sf, opp, ds, r))

    # Tally
    by_seg_verb_verdict = Counter()
    for noun, verb, sf, opp, ds, r in results:
        seg = "dorsal" if sf in {"c","g"} else "coronal"
        by_seg_verb_verdict[(seg, ds, r["verdict"])] += 1
    print("\nVerdict distribution by stem-class × verbalizer:", file=sys.stderr)
    for seg in ("dorsal","coronal"):
        print(f"\n  {seg.upper()}:", file=sys.stderr)
        for verbalizer in ("-i","-a"):
            print(f"    {verbalizer}:", file=sys.stderr)
            for verdict in ("ALTERNATES","NO-ALTERNATION","TRIGGER-ABSENT",
                            "NO-PARADIGM","NO-PAGE"):
                n = by_seg_verb_verdict[(seg, verbalizer, verdict)]
                if n:
                    print(f"      {verdict:18s}: {n:4d}", file=sys.stderr)

    # Write full audit
    OUT_TXT.parent.mkdir(exist_ok=True, parents=True)
    with OUT_TXT.open("w", encoding="utf-8") as f:
        f.write("# Verb-alternation column audit\n")
        f.write("# noun TAB derived_verb TAB stem_final TAB opp TAB ds TAB verdict TAB augment TAB conjugare TAB 2sg.pres TAB 3sg.pres TAB 3sg.subj\n\n")
        for noun, verb, sf, opp, ds, r in sorted(results,
                                                  key=lambda x: (x[2], x[4], x[0])):
            p = r.get("paradigm") or {}
            sg2 = "|".join(p.get("a II-a (tu)", {}).get("prez_ind", []) or [])
            sg3 = "|".join(p.get("a III-a (el, ea)", {}).get("prez_ind", []) or [])
            sj3 = "|".join(p.get("a III-a (el, ea)", {}).get("prez_subj", []) or [])
            aug = r.get("augment")
            conj = r.get("conjugare", "")
            f.write(f"{noun}\t{verb}\t{sf}\t{opp}\t{ds}\t{r['verdict']}\t{aug}\t{conj}\t{sg2}\t{sg3}\t{sj3}\n")
    print(f"\nFull audit written to {OUT_TXT}", file=sys.stderr)


if __name__ == "__main__":
    main()
