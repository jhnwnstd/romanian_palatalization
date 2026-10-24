"""Follow-up checks on check_claims.py results (verb lists, Fisher tests).

Run from the repository root:
    python3 analysis/learnability_checks/followup.py analysis/learnability_checks/check_claims.py
"""
import runpy, io, contextlib, sys
from scipy.stats import fisher_exact
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    g = runpy.run_path(sys.argv[1])
n, top, verbs, verb_outcome = g["n"], g["top"], g["verbs"], g["verb_outcome"]
for name, df in (("top-1000", top), ("full", n)):
    uri = df[df.stem_final.isin(list("cg")) & (df.pl_ev == "uri")]
    rows = [(r.lemma, r.plural, v, verb_outcome(v, s, r.stem_final)) for r in uri.itertuples() for v, s in verbs(r) if s == "-i"]
    print(name, "-uri dorsals with -i verbs:", rows)
print()
for sf in "tds":
    alt = top[(top.stem_final == sf) & (top.pl_ev == "alt")]
    irr = [(r.lemma, r.plural, v) for r in alt.itertuples() for v, s in verbs(r) if verb_outcome(v, s, sf) == "nonpal"]
    print(f"top-1000 {sf}: alternating-plural stems with unpalatalized -i verb:", irr)
print()
print("Fisher, confirmed no-plural vs alternating-plural, dorsal pal-vs-not:", fisher_exact([[1, 6], [90, 71]]))
print("Fisher, confirmed no-plural vs non-alternating, dorsal:", fisher_exact([[1, 6], [8, 41]]))
print("Fisher, coronal no-plural vs alternating:", fisher_exact([[0, 15], [7, 58]]))
