"""Check the data-checkable claims of the learnability paragraph.

Definitions (following the paragraph):
  - A stem is "underspecified" iff the learner sees it alternate; otherwise
    it is stored prespecified.
  - TP exceptions = prespecified stems that occur before a palatalizing
    suffix without alternating.
  - Dorsal rule: c, g before -i/-e. Coronal rule: t, d, s before -i.
Lexical palatals (repo's nde:gimpe class: papugiu, beci, potecea) are
underlying /tʃ dʒ/, not /k g/ or /K G/, so they are out of the domain.

Run from the repository root:
    python3 analysis/learnability_checks/check_claims.py
"""

import math
import re
import sys
from collections import Counter

import pandas as pd

sys.path.insert(0, "src")
from bs4 import BeautifulSoup  # noqa: E402

import dex_utils as du  # noqa: E402
from dex_cache_reader import get_html, load_cache_subset  # noqa: E402

DORS, COR = set("cg"), set("tds")
PAL_ORTH = {"c": ("c", "ch"), "g": ("g", "gh"), "t": ("ț", "t"), "d": ("z", "d"), "s": ("ș", "s")}

d = pd.read_csv("data/romanian_lexicon_with_freq.csv", dtype=str, keep_default_na=False)
d["freq"] = pd.to_numeric(d.freq_ron_news_2024_1M, errors="coerce").fillna(0)
n = d[(d.pos == "N") & d.stem_final.isin(list("cgtdsz"))].copy()
n = n.sort_values("freq", ascending=False).drop_duplicates("lemma")
n["nde"] = n.exception_reason.str.extract(r"nde:(\w+)")[0].fillna("")
n["mut"] = n.mutation.str.lower() == "true"
n["suffix_target"] = n.target_is_suffix.str.lower() == "true"


def verbs(row):
    vs = [v for v in row.derived_verbs.split("|") if v]
    ss = row.deriv_suffixes.split("|")
    return list(zip(vs, ss + [""] * (len(vs) - len(ss))))


def verb_outcome(verb, suffix, sf):
    """'pal' / 'nonpal' for -i verbs, 'avoid' for -a/-ui/other."""
    if suffix != "-i":
        return "avoid"
    if not verb.endswith("i"):
        return None
    stem = verb[:-1]
    if sf in "td" and post_palatal(stem):
        return None
    pal, plain = PAL_ORTH[sf]
    if sf in DORS:
        if stem.endswith(plain):  # chi / ghi
            return "nonpal"
        return "pal" if stem.endswith(pal) else None
    if stem.endswith(pal) or stem.endswith(pal + "t"):  # ț, z, ș, șt
        return "pal"
    if stem.endswith(plain) or stem.endswith(plain + "t"):
        return "nonpal"
    return None


def post_palatal(stem):
    # stop after ș/j: protected by the bleed rule (proʃt), rule-conforming
    return bool(re.search(r"[șj][td]$", stem))


def plural_evidence(r):
    sf = r.stem_final
    if r.nde == "gimpe":
        return "lexical_palatal"
    if r.lemma.endswith("ică") and r.plural.endswith("ele"):
        return "suppletion"  # -ică → -ele replaces the suffix
    if sf in "td" and post_palatal(re.sub("[ăeiu]$", "", r.lemma)):
        return "post_palatal"
    if sf in DORS and not r.mut and re.search(r"ș[cg]", r.lemma[-4:]) and "t" in r.orth_change:
        return "alt"  # șc → șt is the sc alternation
    if r.nde in ("ochi", "paduchi"):
        return "nonalt_front"  # sg/pl -i/-e without alternation
    if r.opportunity == "uri":
        return "uri"
    if r.opportunity in ("i", "e"):
        trig = ("i", "e") if sf in DORS else ("i",)
        if r.opportunity not in trig:
            return "no_trigger"  # coronal -e: not a palatalizing env
        return "alt" if r.mut else "nonalt_front"
    return "none"


n["pl_ev"] = n.apply(plural_evidence, axis=1)
out = []
P = lambda *a: out.append(" ".join(str(x) for x in a))  # noqa: E731


def theta(N):
    return N / math.log(N) if N > 1 else 0


# ---------------------------------------------------------------- Claim 1
P("=" * 72)
P("PREDICTION 1: plural-less nouns resist palatalization in derivatives")
P("=" * 72)
cand = n[(n.pl_ev == "none") & (n.derived_verbs != "") & (n.nde == "")]
cache = load_cache_subset(cand.lemma.tolist())
conf = {}
for lem in cand.lemma:
    h = get_html(cache, lem)
    if h is None:
        continue
    b = du.best_dex_noun_header(du._parse_header_lines(BeautifulSoup(h, "html.parser")), lem)
    if b is not None:
        conf[lem] = b.plural_class
wikt_unc = set(n[n.notes.str.contains("uncountable")].lemma)
st = cand[cand.lemma.map(lambda l: conf.get(l) == "none" or l in wikt_unc)]
P(f"candidates (no plural in data, has derived verb): {len(cand)}; "
  f"DEX/Wiktionary-confirmed no plural: {len(st)}; "
  f"DEX actually lists a plural: {sum(v not in ('none',) for v in conf.values())}")
P(f"proper names (capitalized) among candidates: {sum(cand.lemma.str[0].str.isupper())}")


def tally(df):
    c = Counter()
    ex = {"pal": [], "nonpal": [], "avoid": []}
    for r in df.itertuples():
        if r.stem_final not in PAL_ORTH:
            continue
        for v, s in verbs(r):
            o = verb_outcome(v, s, r.stem_final)
            if o:
                c[o] += 1
                ex[o].append(f"{r.lemma}→{v}")
    return c, ex


groups = {
    "no plural (confirmed)": st,
    "alternating plural": n[n.pl_ev == "alt"],
    "non-alternating plural": n[n.pl_ev.isin(["uri", "nonalt_front", "no_trigger"])],
}
P(f"\n  {'base class':26s} {'cls':8s} {'pal -i':>7s} {'nonpal -i':>9s} {'-a/-ui/etc':>10s} {'%pal':>6s}")
for name, df in groups.items():
    for cls, segs in (("dorsal", DORS), ("coronal", COR)):
        c, ex = tally(df[df.stem_final.isin(segs)])
        tot = sum(c.values())
        P(f"  {name:26s} {cls:8s} {c['pal']:7d} {c['nonpal']:9d} {c['avoid']:10d} "
          f"{(100 * c['pal'] / tot if tot else float('nan')):5.1f}%")
        if name.startswith("no plural"):
            for k in ("pal", "nonpal", "avoid"):
                if ex[k]:
                    P(f"      {k:7s}: {', '.join(ex[k])}")

# ------------------------------------------------------- Claims 2-4 (TP)
top = n[n.opportunity.isin(["i", "e", "uri"]) & (n.freq > 0)].head(1000)
P("\n" + "=" * 72)
P("TOP-1000 LEARNER SAMPLE (same construction as the R script)")
P("=" * 72)
P(f"  {len(top)} nouns; dorsal plural evidence: "
  + str(dict(Counter(top[top.stem_final.isin(DORS)].pl_ev))))


def tp_cell(df, segs, use_verbs, label, no_suffix_targets=False):
    """Stem-level TP: in domain if seen before a palatalizing suffix;
    exception if seen there but never alternates (= prespecified)."""
    N = e = irregular = 0
    exs = []
    sub = df[df.stem_final.isin(segs)]
    if no_suffix_targets:
        sub = sub[~sub.suffix_target]
    for r in sub.itertuples():
        envs = []
        if r.pl_ev in ("alt", "nonalt_front"):
            envs.append(r.pl_ev == "alt")
        if use_verbs:
            for v, s in verbs(r):
                o = verb_outcome(v, s, r.stem_final)
                if o in ("pal", "nonpal"):
                    envs.append(o == "pal")
        if not envs:
            continue
        N += 1
        if not any(envs):
            e += 1
            exs.append(r.lemma)
        elif not all(envs):
            irregular += 1
    th = theta(N)
    ok = "PRODUCTIVE" if N and e <= th else "not productive"
    P(f"  {label:44s} N={N:5d} e={e:4d} θ={th:6.1f} → {ok}"
      + (f"   (+{irregular} alternating stems with a non-alternating env)" if irregular else ""))
    if exs and e <= 45:
        P(f"      exceptions: {', '.join(exs)}")
    return N, e


for sample_name, df in (("top-1000", top), ("full lexicon", n)):
    P(f"\n--- {sample_name} ---")
    for use_verbs in (False, True):
        tag = "plurals + -i verbs" if use_verbs else "plurals only"
        tp_cell(df, DORS, use_verbs, f"dorsal c,g before -i/-e  [{tag}]")
        for sf in "tds":
            tp_cell(df, {sf}, use_verbs, f"coronal {sf} before -i     [{tag}]")
        tp_cell(df, COR, use_verbs, f"coronal t,d,s pooled      [{tag}]")
    tp_cell(df, DORS, True, "dorsal, excl. -ic/-ică suffix targets [pl+verbs]", True)

# Alternative rule carvings (plurals only, top-1000)
P("\n--- alternative rule carvings, top-1000, plurals only ---")
t = top.copy()


def carve(segs, opps, label):
    sub = t[t.stem_final.isin(segs) & t.opportunity.isin(opps) & (t.nde == "")]
    N, e = len(sub), int((~sub.mut).sum())
    P(f"  {label:44s} N={N:5d} e={e:4d} θ={theta(N):6.1f} → "
      + ("PRODUCTIVE" if e <= theta(N) else "not productive"))


carve(DORS | COR, ["i", "e"], "one rule: c,g,t,d,s before -i/-e")
carve(DORS | COR, ["i"], "one rule: c,g,t,d,s before -i")
carve(DORS, ["e"], "  + dorsal-only rule before -e")
carve(COR, ["e"], "coronal t,d,s before -e")

# --------------------------------------------- Claim 3: allomorph selection
P("\n" + "=" * 72)
P("CLAIM: most prespecified dorsals select -uri and -a")
P("=" * 72)
for sample_name, df in (("top-1000", top), ("full lexicon", n[n.pl_ev != "none"])):
    dd = df[df.stem_final.isin(DORS) & (df.pl_ev != "lexical_palatal")]
    pre = dd[dd.pl_ev.isin(["uri", "nonalt_front"])]
    und = dd[dd.pl_ev == "alt"]
    w = lambda x: x.freq.sum()  # noqa: E731
    P(f"\n  [{sample_name}] prespecified dorsals (no plural alternation): {len(pre)}")
    P(f"    plural -uri: {sum(pre.pl_ev == 'uri')}/{len(pre)} types "
      f"({100 * sum(pre.pl_ev == 'uri') / len(pre):.1f}%), "
      f"{100 * w(pre[pre.pl_ev == 'uri']) / w(pre):.1f}% of tokens")
    nf = pre[pre.pl_ev == "nonalt_front"]
    if len(nf):
        P(f"    non-alternating -i/-e: {', '.join(nf.lemma.head(25))}")
    for nm, g in (("prespecified", pre), ("underspecified", und)):
        c = Counter(s for r in g.itertuples() for _, s in verbs(r))
        tot = sum(c.values())
        P(f"    verbalizers on {nm:14s}: n={tot:3d}  "
          + "  ".join(f"{k}:{v}" for k, v in c.most_common(5))
          + (f"   (-a {100 * c['-a'] / tot:.0f}%, -a∪-ui {100 * (c['-a'] + c['-ui']) / tot:.0f}%)" if tot else ""))

# Counterfactual: what if prespecified -uri dorsals had taken -i?
dd = top[top.stem_final.isin(DORS) & (top.pl_ev != "lexical_palatal")]
N_dom = int(dd.pl_ev.isin(["alt", "nonalt_front"]).sum())
e_dom = int((dd.pl_ev == "nonalt_front").sum())
n_uri = int((dd.pl_ev == "uri").sum())
P(f"\n  Counterfactual (top-1000): actual dorsal domain N={N_dom}, e={e_dom}, θ={theta(N_dom):.1f}")
N2, e2 = N_dom + n_uri, e_dom + n_uri
P(f"  if the {n_uri} -uri dorsals had taken -i/-e without alternating: "
  f"N={N2}, e={e2}, θ={theta(N2):.1f} → " + ("PRODUCTIVE" if e2 <= theta(N2) else "NOT productive"))
need = math.floor(theta(N2))
P(f"  break-even: palatalization survives up to {math.floor(theta(N_dom + 1))}–{need} such exceptions")

print("\n".join(out))
