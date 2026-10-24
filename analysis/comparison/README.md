# Cross-language comparison: Romanian, Italian, Czech

This directory holds the cross-language results of the multilingual
extension (design: [`docs/multilingual_design.md`](../../docs/multilingual_design.md),
rules R1–R27). Romanian (recounted read-only from the legacy pipeline's
outputs), Italian and Czech are counted by one shared core
(`src/lexpipe/`), and [`crosslang_summary.csv`](crosslang_summary.csv)
joins their tracked outputs at matched settings. This README is the
report: what was measured and how, the comparison tables, and a check of
the learnability paragraph of the Gorman & Winstead draft (`latex.tex`)
against all three languages.

Romanian numbers are **harmonized** (counted by
`lexpipe.langs.romanian_compare` with the shared definitions) unless
they are labelled **legacy** (the R analysis of the Romanian pipeline,
`analysis/romanian_tp_summary_*.csv`, whose top-1000 frame is the
draft's productivity table).

## Summary

- **The rule family separates the languages.** Romanian dorsals
  palatalize before -i and -e, and t, d, s before -i, each rule well
  inside the Tolerance threshold (dorsals before -i/-e: 2,874 noun roots,
  180 exceptions, θ = 360.9). Czech dorsals palatalize before -i (3,531
  roots, 3 exceptions), as do r and t, d, n, but no rule stated on the
  vowel covers -e, because Czech -e triggers in some endings and not in
  others. Italian
  dorsal nouns have no productive rule at all before -i (≥ 1 token: 520
  roots, 264 keep the velar and 256 palatalize, θ = 83.1), and Italian
  coronals never alternate (Table 1).
- **Each language's best fit is a grammar of its own** (Table 3):
  Romanian, a dorsal rule over front vowels plus a coronal rule over -i;
  Italian, -i nested by stress (post-tonic → retain, elsewhere →
  palatalize); Czech, a trigger indexed by suffix, with no exception
  once a palatalizing allomorph is chosen, so that all the variation is in
  which allomorph a stem takes.
- **The learnability paragraph** (claims (a)–(g) below): one claim is
  supported outright (f, the Romanian counts) and one as a way of
  counting (d), which makes Romanian and Czech pass and Italian fail.
  Three are mixed (c, e, g) and two are not supported (a, as a default
  for stems not yet seen alternating; b, plural-less nouns resisting
  palatalization). The decisive evidence is that bases without
  alternation evidence palatalize in their derivatives as often as
  alternating ones: in Czech, 155 of 156 units at ≥ 1 token against 480
  of 483; among Romanian dorsals, 17 of 17 against 59 of 59. A Czech
  learner who underspecifies a stem only after seeing it alternate would
  mispredict 48% of the 48,143 inflected forms of target stems (Table 8).
- **Five results contradict the draft directly:**
  - the Romanian spell-out TP fails: 51 of 128 dorsal roots whose plural
    alternates take only -a/-ui (θ = 26.4);
  - Romanian -uri dorsals palatalize whenever they take a front-vowel
    derivational suffix (10 of 10 units);
  - Czech mixed stems: -log persons take NOM.PL -ové but palatalize in
    LOC.PL -ích (107 of 109);
  - Czech LOC.SG -u is the default of every hard inanimate masculine,
    velar or not;
  - Italian first-conjugation verbs never palatalize (940 of 940), so a
    root that alternates in its plural (pratici) keeps the velar in its
    verb (pratichi).

| claim of the paragraph | Romanian | Italian | Czech | overall |
|---|---|---|---|---|
| (a) underspecify only after seeing an alternation | not supported | mixed | not supported | **not supported** as a default |
| (b) plural-less nouns resist palatalization in derivatives | not supported (dorsals) | mixed, small N | not supported | **not supported** |
| (c) the inflection-dependence pattern follows from learner evidence | mixed | mixed (-ic- driven) | not supported | **mixed** |
| (d) TP exceptions = attested non-alternating stems | rule passes | no single rule | rule passes | **supported** as a way of counting |
| (e) allomorph selection structures the evidence | mixed | consistent (one plural exponent) | mixed | **mixed** |
| (f) Romanian: dorsals before -i/-e, t/d/s before -i | supported | – | – | **supported** |
| (g) Italian fails → no single rule, as in Krämer's speakers | – | mixed | – | **mixed** |

Contents: [Methods in brief](#methods-in-brief) ·
[Tables](#tables) ·
[The learnability paragraph, checked against three languages](#the-learnability-paragraph-checked-against-three-languages) ·
[Where the data contradict the draft](#where-the-data-contradict-the-draft) ·
[Open points](#open-points) ·
[Romanian (harmonized, read-only)](#romanian-harmonized-read-only)

## Methods in brief

### Running and files

```bash
python scripts/lang_pipeline.py --lang it        # Italian, every stage
python scripts/lang_pipeline.py --lang cs        # Czech
python scripts/lang_pipeline.py --lang ro        # the Romanian recount (read-only)
PYTHONPATH=src python -m lexpipe.compare         # crosslang_summary.csv + the tables below
PYTHONPATH=src python -m lexpipe.compare --check # exit status 1 if they are stale
```

`lexpipe.compare` reads only tracked files (under `analysis/`), never
the gitignored `data/`, so a checkout can rebuild the summary. It runs in
about a second, writes only under `analysis/comparison/` and is skipped
while its manifest (`data/romanian_compare/manifests/crosslang.manifest.tsv`:
inputs, version and the module's own SHA-256) matches. Every count is
copied from a language's own table or report. Two things are computed
here:

- the odds ratio and Fisher p of a 2×2 summed over pairs or suffixes,
  with `lexpipe.stats`; for Italian they equal the Italian README's;
- θ and the verdicts where a source gives only counts: the pooled
  Italian -logo/-fago/-urgo rows, the Czech -log rows, and the opposite
  rule of the Czech spell-out and legacy Romanian rows.

Files:

- [`crosslang_summary.csv`](crosslang_summary.csv): one row per
  comparison cell, with the columns defined in the module docstring of
  `src/lexpipe/compare.py`. `section` groups the rows (r17, learners,
  model_fit, spellout, noevidence, covariation, logtable, phonology,
  legacy), `matched` marks the shared settings, `input` names the tracked
  file each row comes from, and `source`/`licence_class` pass the R9 gate.
- `romanian_*`: the harmonized Romanian recount and its parity with the
  legacy R tables (section at the end, generated by
  `lexpipe.langs.romanian_compare`).

### Romanian

The lexicon is the legacy pipeline's (`run_pipeline.sh`: Wiktionary
noun pages, about 70% from ro.wiktionary, checked against DEX, which also
filled the plurals Wiktionary lacked for 837 nouns; one plural
per lemma; Leipzig `ron_news_2024_1M`, 18.6 million tokens, lowercased
counts). The legacy R analysis counts lemma rows per stem-final ×
plural cell. Outcome-conditioned exclusions remove the z→zi "foreign"
plurals and items already palatal. Its top-1000 frame (the most frequent
nouns with a target final and an -i or -e plural) is the draft's table.
The harmonized recount
(`lexpipe.langs.romanian_compare`, design §7) counts root types with the
shared R16–R19 definitions. It replaces the exclusions with outcome-blind
flags (subset `core` drops them), uses the R24 frame of the most frequent
nouns of any final, and adds evidence sets E1 (verb cells) and E2
(derivatives) and the R20, R21, R23 and R27 analyses. It reproduces every
legacy row exactly (parity: 37/37 numeric rows of each legacy table).
Where the definitions differ, the counts differ: legacy c + -i is 473
roots with 0 exceptions, harmonized dorsal before -i is 1,616 roots with
49 exceptions in `all` and 852 with 0 in `core`. The legacy fields are
defined in [`docs/data_dict.md`](../../docs/data_dict.md), and the
reconciliation, root by root, is in
[`romanian_parity.txt`](romanian_parity.txt). The Logical Phonology
analysis is `src/phonology/analyses/romanian_palatalization.py`.

### Italian

[`docs/italian/README.md`](../../docs/italian/README.md) and
[`docs/italian/phonology.md`](../../docs/italian/phonology.md). The
lexicon comes from kaikki.org's English and Italian Wiktionary extracts,
with Morph-it! as QC and Leipzig `ita_news_2024_1M` (19.8 million
tokens). A root groups POS conversions and gender pairs (amico N/ADJ,
amica). The cells are -i (masculine plural: the trigger), -e (feminine
plural: a non-trigger cell) and -a. Stress comes from the IPA (missing
for 37% of lexemes), -ico/-iaco suffix targets are identified
morphologically (R26), velar-stem verbs enter as E1 and derivatives as
E2. Read the ≥ 1 token rows first. 96.5% of the dictionary's dorsal -i
plurals (6,890 of 7,137) are defaults generated by `Module:it-common`
(-ico/-aco/-ologo palatal, every other -co/-go velar), so rows at τ =
dictionary partly restate a spelling rule. Of the plurals attested at
≥ 1 token, 95.5% are defaults too, but they are attested.

### Czech

[`docs/czech/README.md`](../../docs/czech/README.md) and
[`docs/czech/phonology.md`](../../docs/czech/phonology.md). The lexicon
comes from kaikki.org's English and Czech Wiktionary extracts, with
MorfFlex as QC, DeriNet for derivations and czetyl for etymology (all
three CC BY-NC-SA: never in a tracked table; the Czech text reports show
their aggregates in sections marked NC, and this report uses none of
them). The corpus is Leipzig `ces_news_2024_1M`: 15.0 million tokens,
case-sensitive counts with a sentence-initial correction, plus a
**strict** variant that does not count forms of other lexemes or forms
shared with another cell.

Cells are case × number. Each palatalizing (P) allomorph has its
competitor (N): NOM.PL -i vs -ové/-é, LOC.SG -e vs -u, LOC.PL -ích vs
-ách/-ech, VOC -e vs -u. There are three TP analyses: the rule TP (does
the stem change before a P allomorph?), the competition TP (which
allomorph does a stem take?) and spell-out leave-one-cell-out.
Palatalization of t/d/n is invisible in spelling and inferred (e = 0 by
construction). R26 suffix targets (-ský/-cký, -ka, -ek, …) are
identified morphologically.

### Shared definitions and matched settings

Unit: root types (R4). N = roots attested before the trigger; e =
exceptions; θ = N/ln N; a rule is productive when e ≤ θ. Both rules are
evaluated (R18): *palatalize* (e = roots that keep the plain consonant)
and *retain* (e = roots that change). Roots whose attested forms all
avoid the trigger (Romanian -uri, Czech -ové/-u/-ách, the Italian
feminine -e plural in the -i cell) are **shielded**: they are neither N
nor exceptions (R16). Doublets are coded by their dominant variant (R19).

Simulated learners draw 1,000 vocabularies by token frequency (lemma
basis, seed 20080715) and report the share of learners for whom each
rule is productive. Nested rules are evaluated in Yang's order.
`matched` rows share E0 (nominal inflection), τ ∈ {dictionary, ≥ 1
token}, dominant coding, the lemma frequency basis, nouns (the Italian
`nouns` subset, the Czech noun rows, the Romanian noun lexicon) and the
full lexicon (learners: V ∈ {1000, 5000}).

### What the matched settings do not match

| issue | Romanian | Italian | Czech |
|---|---|---|---|
| dictionary forms | one listed plural per lemma (no doublets at E0) | 96.5% of dorsal -i plurals are module defaults | paradigm tables are template-generated (biolozích has 0 tokens) |
| corpus counts | lowercased only | folded, plus ≥ 1 lowercase token (names out) | case-sensitive, sentence-initial correction, strict variant |
| syncretism and homographs | not assessed | -are verb forms = noun plurals (càrichi) | LOC.SG -u = GEN/DAT; VOC -u = DAT/LOC; homograph lemmas (jak 'yak') |
| corpus tagging | none | none | none |
| stress | not used | 37% of lexemes without stress; vowel + -ico syllabification undecided | fixed initial stress, not used |
| suffix targets (R26) | counted as roots | counted as roots (`suffix_once`: once per suffix) | counted as roots |
| doublets across cells (R17a) | n/a (one plural) | a root palatal before -i, velar before -e is an exception to both rules | coded by token majority across cells |
| strata (R22) | etymon language only (no internal formations) | Wiktionary templates | Wiktionary (czetyl NC, not used here) |
| coronals | t, d, s; z apart | coded unchanged by construction | t/d/n inferred from spelling; r, s/z, l apart |
| verbs (E1) | 2sg and 3sg subjunctive | 2sg, 3sg, 1pl, 3sg subjunctive | -it/-et verbs |

## Tables

The tables below are generated from `crosslang_summary.csv` by
`python -m lexpipe.compare`; edit the code, not the block.

<!-- BEGIN crosslang: generated by lexpipe.compare from crosslang_summary.csv; edit the code, not this block -->
Unit: root types. θ = N/ln N; e (pal) = exceptions to *palatalize* (roots attested before the trigger that keep the plain consonant); e (ret) = exceptions to *retain* (roots that change); productive = the rule(s) with e ≤ θ (**both** only for very small N). τ = dictionary counts every listed form; τ ≥ 1 token counts forms with at least one corpus token. Shielded = roots attested in the cell only with a non-front or non-palatalizing allomorph. Matched settings: E0, dominant doublet coding, lemma frequency basis, nouns; every column of every row is in the CSV.

### Table 1. The R17 rule family at matched settings

| rule-family row | language | target | τ = dictionary: N: e (pal) / e (ret), θ → productive | τ ≥ 1 token | shielded (dict. / ≥ 1) |
|---|---|---|---|---|---|
| R17a dorsal/_[+front V] | Romanian | dorsal | 2,874: 180 / 2,694, θ 360.9 → **palatalize** | 601: 39 / 562, θ 93.9 → **palatalize** | 501 / 93 |
| R17a dorsal/_[+front V] | Italian | dorsal | 2,550: 1,829 / 776, θ 325.1 → neither | 825: 592 / 256, θ 122.9 → neither | 8 / 0 |
| R17a dorsal/_[+front V] | Czech | dorsal | 7,760: 3,367 / 6,492, θ 866.4 → neither | 2,666: 1,006 / 1,719, θ 338.0 → neither | 0 / 0 |
| R17b dorsal/_[+high,+front] | Romanian | dorsal | 1,616: 49 / 1,567, θ 218.7 → **palatalize** | 397: 9 / 388, θ 66.3 → **palatalize** | 1,759 / 297 |
| R17b dorsal/_[+high,+front] | Italian | dorsal | 1,481: 739 / 758, θ 202.9 → neither | 520: 264 / 256, θ 83.1 → neither | 1,077 / 305 |
| R17b dorsal/_[+high,+front] | Czech | dorsal | 3,531: 3 / 3,528, θ 432.2 → **palatalize** | 837: 0 / 837, θ 124.4 → **palatalize** | 0 / 0 |
| R17b dorsal/_[+high,+front], non-trigger -e | Romanian | dorsal | 1,258: 131 / 1,127, θ 176.3 → **palatalize** | 204: 30 / 174, θ 38.4 → **palatalize** | 2,117 / 490 |
| R17b dorsal/_[+high,+front], non-trigger -e | Italian | dorsal | 1,132: 1,114 / 18, θ 161.0 → **retain** | 338: 338 / 0, θ 58.0 → **retain** | 1,426 / 487 |
| R17b dorsal/_[+high,+front], non-trigger -e | Czech | dorsal | 7,696: 3,884 / 3,827, θ 860.0 → neither | 2,379: 1,245 / 1,134, θ 306.0 → neither | 0 / 0 |
| R17c coronal/_[+high,+front] | Romanian | t, d, s | 2,898: 124 / 2,774, θ 363.5 → **palatalize** | 871: 10 / 861, θ 128.7 → **palatalize** | 4,384 / 1,321 |
| R17c coronal/_[+high,+front] | Romanian | z | 126: 116 / 10, θ 26.1 → **retain** | 48: 46 / 2, θ 12.4 → **retain** | 941 / 171 |
| R17c coronal/_[+high,+front] | Italian | t, d, s, z | 10,277: 10,277 / 0, θ 1112.5 → **retain** | 3,277: 3,277 / 0, θ 404.8 → **retain** | 3,413 / 1,258 |
| R17c coronal/_[+high,+front] | Czech | t, d, n (inferred) | 1,156: 0 / 1,156, θ 163.9 → **palatalize** | 322: 0 / 322, θ 55.8 → **palatalize** | 0 / 0 |
| R17c coronal/_[+high,+front] | Czech | r | 376: 5 / 371, θ 63.4 → **palatalize** | 186: 1 / 185, θ 35.6 → **palatalize** | 0 / 0 |
| R17c coronal/_[+high,+front] | Czech | s, z | 239: 239 / 0, θ 43.6 → **retain** | 114: 114 / 0, θ 24.1 → **retain** | 0 / 0 |
| R17c coronal/_[+high,+front] | Czech | l | 449: 449 / 0, θ 73.5 → **retain** | 220: 220 / 0, θ 40.8 → **retain** | 0 / 0 |
| R17c coronal/_[+high,+front], non-trigger -e | Romanian | t, d, s | 3,246: 3,232 / 14, θ 401.5 → **retain** | 1,087: 1,083 / 4, θ 155.5 → **retain** | 4,036 / 1,105 |
| R17c coronal/_[+high,+front], non-trigger -e | Romanian | z | 864: 864 / 0, θ 127.8 → **retain** | 155: 155 / 0, θ 30.7 → **retain** | 203 / 64 |
| R17c coronal/_[+high,+front], non-trigger -e | Italian | t, d, s, z | 5,039: 5,039 / 0, θ 591.1 → **retain** | 1,452: 1,452 / 0, θ 199.4 → **retain** | 8,651 / 3,083 |
| R17c coronal/_[+high,+front], non-trigger -e | Czech | t, d, n (inferred) | 5,758: 3,599 / 2,354, θ 665.0 → neither | 2,227: 1,384 / 846, θ 288.9 → neither | 0 / 0 |
| R17c coronal/_[+high,+front], non-trigger -e | Czech | r | 1,497: 1,160 / 467, θ 204.8 → neither | 645: 504 / 141, θ 99.7 → neither | 0 / 0 |
| R17c coronal/_[+high,+front], non-trigger -e | Czech | s, z | 1,040: 1,040 / 0, θ 149.7 → **retain** | 451: 451 / 0, θ 73.8 → **retain** | 0 / 0 |
| R17c coronal/_[+high,+front], non-trigger -e | Czech | l | 1,442: 1,442 / 0, θ 198.2 → **retain** | 659: 659 / 0, θ 101.5 → **retain** | 0 / 0 |

### Table 1b. The same in the R24 frame of the 1000 most frequent noun roots

| rule-family row | language | target | τ = dictionary: N: e (pal) / e (ret), θ → productive | τ ≥ 1 token | shielded (dict. / ≥ 1) |
|---|---|---|---|---|---|
| R17a dorsal/_[+front V] | Romanian | dorsal | 25: 0 / 25, θ 7.8 → **palatalize** | 22: 0 / 22, θ 7.1 → **palatalize** | 14 / 10 |
| R17a dorsal/_[+front V] | Italian | dorsal | 42: 35 / 12, θ 11.2 → neither | 40: 33 / 12, θ 10.8 → neither | 0 / 0 |
| R17a dorsal/_[+front V] | Czech | dorsal | 113: 61 / 52, θ 23.9 → neither | 111: 59 / 52, θ 23.6 → neither | 0 / 0 |
| R17b dorsal/_[+high,+front] | Romanian | dorsal | 19: 0 / 19, θ 6.5 → **palatalize** | 17: 0 / 17, θ 6.0 → **palatalize** | 20 / 15 |
| R17b dorsal/_[+high,+front] | Italian | dorsal | 28: 16 / 12, θ 8.4 → neither | 27: 15 / 12, θ 8.2 → neither | 14 / 13 |
| R17b dorsal/_[+high,+front] | Czech | dorsal | 77: 0 / 77, θ 17.7 → **palatalize** | 59: 0 / 59, θ 14.5 → **palatalize** | 0 / 0 |
| R17b dorsal/_[+high,+front], non-trigger -e | Romanian | dorsal | 6: 0 / 6, θ 3.3 → **palatalize** | 5: 0 / 5, θ 3.1 → **palatalize** | 33 / 27 |
| R17b dorsal/_[+high,+front], non-trigger -e | Italian | dorsal | 22: 22 / 0, θ 7.1 → **retain** | 21: 21 / 0, θ 6.9 → **retain** | 20 / 19 |
| R17b dorsal/_[+high,+front], non-trigger -e | Czech | dorsal | 111: 77 / 34, θ 23.6 → neither | 109: 75 / 34, θ 23.2 → neither | 0 / 0 |
| R17c coronal/_[+high,+front] | Romanian | t, d, s | 93: 1 / 92, θ 20.5 → **palatalize** | 89: 1 / 88, θ 19.8 → **palatalize** | 202 / 173 |
| R17c coronal/_[+high,+front] | Romanian | z | 1: 1 / 0, θ – → – | 1: 1 / 0, θ – → – | 17 / 16 |
| R17c coronal/_[+high,+front] | Italian | t, d, s, z | 218: 218 / 0, θ 40.5 → **retain** | 217: 217 / 0, θ 40.3 → **retain** | 109 / 106 |
| R17c coronal/_[+high,+front] | Czech | t, d, n (inferred) | 17: 0 / 17, θ 6.0 → **palatalize** | 16: 0 / 16, θ 5.8 → **palatalize** | 0 / 0 |
| R17c coronal/_[+high,+front] | Czech | r | 13: 0 / 13, θ 5.1 → **palatalize** | 13: 0 / 13, θ 5.1 → **palatalize** | 0 / 0 |
| R17c coronal/_[+high,+front] | Czech | s, z | 15: 15 / 0, θ 5.5 → **retain** | 14: 14 / 0, θ 5.3 → **retain** | 0 / 0 |
| R17c coronal/_[+high,+front] | Czech | l | 27: 27 / 0, θ 8.2 → **retain** | 27: 27 / 0, θ 8.2 → **retain** | 0 / 0 |
| R17c coronal/_[+high,+front], non-trigger -e | Romanian | t, d, s | 115: 115 / 0, θ 24.2 → **retain** | 113: 113 / 0, θ 23.9 → **retain** | 180 / 149 |
| R17c coronal/_[+high,+front], non-trigger -e | Romanian | z | 12: 12 / 0, θ 4.8 → **retain** | 12: 12 / 0, θ 4.8 → **retain** | 6 / 5 |
| R17c coronal/_[+high,+front], non-trigger -e | Italian | t, d, s, z | 119: 119 / 0, θ 24.9 → **retain** | 116: 116 / 0, θ 24.4 → **retain** | 208 / 207 |
| R17c coronal/_[+high,+front], non-trigger -e | Czech | t, d, n (inferred) | 177: 123 / 54, θ 34.2 → neither | 172: 118 / 54, θ 33.4 → neither | 0 / 0 |
| R17c coronal/_[+high,+front], non-trigger -e | Czech | r | 46: 38 / 8, θ 12.0 → **retain** | 45: 37 / 8, θ 11.8 → **retain** | 0 / 0 |
| R17c coronal/_[+high,+front], non-trigger -e | Czech | s, z | 44: 44 / 0, θ 11.6 → **retain** | 44: 44 / 0, θ 11.6 → **retain** | 0 / 0 |
| R17c coronal/_[+high,+front], non-trigger -e | Czech | l | 63: 63 / 0, θ 15.2 → **retain** | 59: 59 / 0, θ 14.5 → **retain** | 0 / 0 |

### Table 1c. Italian nouns and adjectives together (the Italian headline subset)

| rule-family row | language | target | τ = dictionary: N: e (pal) / e (ret), θ → productive | τ ≥ 1 token | shielded (dict. / ≥ 1) |
|---|---|---|---|---|---|
| R17a dorsal/_[+front V] | Italian | dorsal | 7,537: 7,169 / 5,584, θ 844.2 → neither | 1,842: 1,494 / 980, θ 245.0 → neither | 8 / 0 |
| R17b dorsal/_[+high,+front] | Italian | dorsal | 6,612: 1,061 / 5,571, θ 751.7 → neither | 1,302: 322 / 980, θ 181.5 → neither | 933 / 540 |
| R17b dorsal/_[+high,+front], non-trigger -e | Italian | dorsal | 6,566: 6,553 / 17, θ 747.0 → **retain** | 1,248: 1,248 / 0, θ 175.1 → **retain** | 979 / 594 |
| R17c coronal/_[+high,+front] | Italian | t, d, s, z | 15,935: 15,935 / 0, θ 1646.8 → **retain** | 5,601: 5,601 / 0, θ 649.0 → **retain** | 3,380 / 1,497 |
| R17c coronal/_[+high,+front], non-trigger -e | Italian | t, d, s, z | 10,063: 10,063 / 0, θ 1091.8 → **retain** | 3,424: 3,424 / 0, θ 420.7 → **retain** | 9,252 / 3,674 |

### Table 1d. Without suffix targets (subset core; Italian nouns and adjectives)

| rule-family row | language | target | τ = dictionary: N: e (pal) / e (ret), θ → productive | τ ≥ 1 token | shielded (dict. / ≥ 1) |
|---|---|---|---|---|---|
| R17a dorsal/_[+front V] | Romanian | dorsal | 1,556: 0 / 1,556, θ 211.7 → **palatalize** | 324: 0 / 324, θ 56.0 → **palatalize** | 442 / 83 |
| R17a dorsal/_[+front V] | Italian | dorsal | 5,298: 5,178 / 3,446, θ 617.8 → neither | 1,167: 992 / 439, θ 165.2 → neither | 8 / 0 |
| R17a dorsal/_[+front V] | Czech | dorsal | 3,293: 1,699 / 2,485, θ 406.6 → neither | 1,125: 522 / 627, θ 160.1 → neither | 0 / 0 |
| R17b dorsal/_[+high,+front] | Romanian | dorsal | 852: 0 / 852, θ 126.3 → **palatalize** | 239: 0 / 239, θ 43.6 → **palatalize** | 1,146 / 168 |
| R17b dorsal/_[+high,+front] | Italian | dorsal | 4,447: 1,025 / 3,436, θ 529.4 → neither | 747: 308 / 439, θ 112.9 → neither | 859 / 420 |
| R17b dorsal/_[+high,+front] | Czech | dorsal | 1,559: 1 / 1,558, θ 212.1 → **palatalize** | 369: 0 / 369, θ 62.4 → **palatalize** | 0 / 0 |
| R17b dorsal/_[+high,+front], non-trigger -e | Romanian | dorsal | 704: 0 / 704, θ 107.4 → **palatalize** | 85: 0 / 85, θ 19.1 → **palatalize** | 1,294 / 322 |
| R17b dorsal/_[+high,+front], non-trigger -e | Italian | dorsal | 4,599: 4,589 / 13, θ 545.3 → **retain** | 758: 758 / 0, θ 114.3 → **retain** | 707 / 409 |
| R17b dorsal/_[+high,+front], non-trigger -e | Czech | dorsal | 3,237: 1,913 / 1,336, θ 400.5 → neither | 1,001: 621 / 380, θ 144.9 → neither | 0 / 0 |
| R17c coronal/_[+high,+front] | Romanian | t, d, s | 1,975: 6 / 1,969, θ 260.3 → **palatalize** | 663: 0 / 663, θ 102.1 → **palatalize** | 4,355 / 1,309 |
| R17c coronal/_[+high,+front] | Romanian | z | 126: 116 / 10, θ 26.1 → **retain** | 48: 46 / 2, θ 12.4 → **retain** | 941 / 171 |
| R17c coronal/_[+high,+front] | Czech | t, d, n (inferred) | 1,156: 0 / 1,156, θ 163.9 → **palatalize** | 322: 0 / 322, θ 55.8 → **palatalize** | 0 / 0 |
| R17c coronal/_[+high,+front] | Czech | r | 376: 5 / 371, θ 63.4 → **palatalize** | 186: 1 / 185, θ 35.6 → **palatalize** | 0 / 0 |
| R17c coronal/_[+high,+front] | Czech | s, z | 239: 239 / 0, θ 43.6 → **retain** | 114: 114 / 0, θ 24.1 → **retain** | 0 / 0 |
| R17c coronal/_[+high,+front] | Czech | l | 449: 449 / 0, θ 73.5 → **retain** | 220: 220 / 0, θ 40.8 → **retain** | 0 / 0 |
| R17c coronal/_[+high,+front], non-trigger -e | Romanian | t, d, s | 3,231: 3,217 / 14, θ 399.8 → **retain** | 1,080: 1,076 / 4, θ 154.6 → **retain** | 3,099 / 892 |
| R17c coronal/_[+high,+front], non-trigger -e | Romanian | z | 864: 864 / 0, θ 127.8 → **retain** | 155: 155 / 0, θ 30.7 → **retain** | 203 / 64 |
| R17c coronal/_[+high,+front], non-trigger -e | Czech | t, d, n (inferred) | 5,758: 3,599 / 2,354, θ 665.0 → neither | 2,227: 1,384 / 846, θ 288.9 → neither | 0 / 0 |
| R17c coronal/_[+high,+front], non-trigger -e | Czech | r | 1,497: 1,160 / 467, θ 204.8 → neither | 645: 504 / 141, θ 99.7 → neither | 0 / 0 |
| R17c coronal/_[+high,+front], non-trigger -e | Czech | s, z | 1,040: 1,040 / 0, θ 149.7 → **retain** | 451: 451 / 0, θ 73.8 → **retain** | 0 / 0 |
| R17c coronal/_[+high,+front], non-trigger -e | Czech | l | 1,442: 1,442 / 0, θ 198.2 → **retain** | 659: 659 / 0, θ 101.5 → **retain** | 0 / 0 |

### Table 2. Simulated learners (1000 vocabularies by token frequency, E0, lemma basis): shares of learners by productive rule

| language | cell | τ | V = 1000 | V = 5000 |
|---|---|---|---|---|
| Romanian | R17a dorsal -i/-e [with_nde] (nested: velar_front_lemma) | dictionary | palatalize 1.000 | palatalize 1.000 |
| Romanian | R17a dorsal -i/-e [with_nde] (nested: velar_front_lemma) | ge1 | palatalize 1.000 | palatalize 1.000 |
| Romanian | R17b dorsal -i | dictionary | palatalize 1.000 | palatalize 1.000 |
| Romanian | R17b dorsal -i | ge1 | palatalize 1.000 | palatalize 1.000 |
| Romanian | R17b dorsal -e | dictionary | palatalize 1.000 | palatalize 1.000 |
| Romanian | R17b dorsal -e | ge1 | palatalize 1.000 | palatalize 1.000 |
| Romanian | R17c coronal -i (nested: palatal_context) | dictionary | palatalize 1.000 | palatalize 1.000 |
| Romanian | R17c coronal -i (nested: palatal_context) | ge1 | palatalize 1.000 | palatalize 1.000 |
| Romanian | R17c coronal -e | dictionary | retain 1.000 | retain 1.000 |
| Romanian | R17c coronal -e | ge1 | retain 1.000 | retain 1.000 |
| Romanian | R17c coronal -i/-e (nested: e_plural) | dictionary | conditioned 1.000 | conditioned 1.000 |
| Romanian | R17c coronal -i/-e (nested: e_plural) | ge1 | conditioned 1.000 | conditioned 1.000 |
| Romanian | R17c z -i (nested: borrowed) | dictionary | retain 0.688, neither 0.308 | retain 1.000 |
| Romanian | R17c z -i (nested: borrowed) | ge1 | retain 0.632, neither 0.368 | retain 1.000 |
| Romanian | R17c z -e | dictionary | retain 1.000 | retain 1.000 |
| Romanian | R17c z -e | ge1 | retain 1.000 | retain 1.000 |
| Italian | R17a dorsal -i/-e [with_nde] (nested: velar_front_lemma) | dictionary | retain 0.513, neither 0.487 | neither 1.000 |
| Italian | R17a dorsal -i/-e [with_nde] (nested: velar_front_lemma) | ge1 | neither 0.528, retain 0.472 | neither 1.000 |
| Italian | R17b dorsal -i | dictionary | neither 0.966, retain 0.034 | neither 1.000 |
| Italian | R17b dorsal -i | ge1 | neither 0.977, retain 0.023 | neither 1.000 |
| Italian | R17b dorsal -e | dictionary | retain 1.000 | retain 1.000 |
| Italian | R17b dorsal -e | ge1 | retain 1.000 | retain 1.000 |
| Italian | R17c coronal -i | dictionary | retain 1.000 | retain 1.000 |
| Italian | R17c coronal -i | ge1 | retain 1.000 | retain 1.000 |
| Italian | it_stress dorsal -i (nested: post_tonic) | dictionary | conditioned 0.963, retain 0.034 | conditioned 0.993, neither 0.007 |
| Italian | it_stress dorsal -i (nested: post_tonic) | ge1 | conditioned 0.975, retain 0.023 | conditioned 1.000 |
| Italian | R17b dorsal -i [suffix_once] | dictionary | retain 0.843, neither 0.157 | neither 1.000 |
| Italian | R17b dorsal -i [suffix_once] | ge1 | retain 0.811, neither 0.189 | neither 1.000 |
| Italian | it_stress dorsal -i [suffix_once] (nested: post_tonic) | dictionary | retain 0.843, conditioned 0.157 | neither 0.653, conditioned 0.347 |
| Italian | it_stress dorsal -i [suffix_once] (nested: post_tonic) | ge1 | retain 0.811, conditioned 0.189 | conditioned 0.669, neither 0.331 |
| Italian | it_stress dorsal -i [orth_ico_once] (nested: post_tonic) | dictionary | retain 0.990, conditioned 0.010 | retain 0.999 |
| Italian | it_stress dorsal -i [orth_ico_once] (nested: post_tonic) | ge1 | retain 0.988, conditioned 0.012 | retain 0.998 |
| Czech | R17a dorsal -e/-i | dictionary | neither 1.000 | neither 1.000 |
| Czech | R17a dorsal -e/-i | ge1 | neither 1.000 | neither 1.000 |
| Czech | R17b dorsal -i | dictionary | palatalize 1.000 | palatalize 1.000 |
| Czech | R17b dorsal -i | ge1 | palatalize 1.000 | palatalize 1.000 |
| Czech | R17b dorsal -e | dictionary | neither 1.000 | neither 1.000 |
| Czech | R17b dorsal -e | ge1 | neither 1.000 | neither 1.000 |
| Czech | R17c coronal_inferred -i | dictionary | palatalize 1.000 | palatalize 1.000 |
| Czech | R17c coronal_inferred -i | ge1 | palatalize 1.000 | palatalize 1.000 |
| Czech | R17c r -i | dictionary | palatalize 1.000 | palatalize 1.000 |
| Czech | R17c r -i | ge1 | palatalize 1.000 | palatalize 1.000 |
| Czech | R17c r -e | dictionary | retain 0.955, neither 0.045 | neither 0.779, retain 0.221 |
| Czech | R17c r -e | ge1 | retain 0.960, neither 0.040 | retain 0.540, neither 0.460 |
| Czech | R17c sibilant -i | dictionary | retain 1.000 | retain 1.000 |
| Czech | R17c sibilant -i | ge1 | retain 1.000 | retain 1.000 |
| Czech | R17c lateral -i | dictionary | retain 1.000 | retain 1.000 |
| Czech | R17c lateral -i | ge1 | retain 1.000 | retain 1.000 |
| Czech | competition velar LOC.SG m.in | dictionary | retain 1.000 | retain 1.000 |
| Czech | competition velar LOC.SG m.in | ge1 | retain 1.000 | retain 1.000 |
| Czech | competition baseline LOC.SG m.in | dictionary | retain 1.000 | retain 1.000 |
| Czech | competition baseline LOC.SG m.in | ge1 | retain 1.000 | retain 1.000 |
| Czech | competition velar NOM.PL m.an | dictionary | palatalize 0.950, neither 0.050 | neither 0.958, palatalize 0.042 |
| Czech | competition velar NOM.PL m.an | ge1 | palatalize 0.999 | palatalize 0.915, neither 0.085 |
| Czech | competition baseline NOM.PL m.an | dictionary | retain 0.994, neither 0.005 | retain 0.765, neither 0.235 |
| Czech | competition baseline NOM.PL m.an | ge1 | retain 0.988, neither 0.008 | retain 0.944, neither 0.056 |
| Czech | competition velar LOC.PL m.in | dictionary | palatalize 1.000 | palatalize 1.000 |
| Czech | competition velar LOC.PL m.in | ge1 | palatalize 1.000 | palatalize 1.000 |
| Czech | competition baseline LOC.PL m.in | dictionary | palatalize 1.000 | palatalize 1.000 |
| Czech | competition baseline LOC.PL m.in | ge1 | palatalize 1.000 | palatalize 1.000 |
| Czech | spellout_loo velar LOC.SG m.in [underspecified] | dictionary | retain 1.000 | retain 1.000 |
| Czech | spellout_loo velar LOC.SG m.in [underspecified] | ge1 | retain 1.000 | retain 1.000 |
| Czech | spellout_loo velar NOM.PL m.an [underspecified] | dictionary | palatalize 0.951, neither 0.049 | neither 0.939, palatalize 0.061 |
| Czech | spellout_loo velar NOM.PL m.an [underspecified] | ge1 | palatalize 0.999 | palatalize 1.000 |

### Table 2b. Simulated learners with verb cells added (E1)

| language | cell | τ | V = 1000 | V = 5000 |
|---|---|---|---|---|
| Romanian | R17a dorsal -i/-e [with_nde] (nested: velar_front_lemma) | dictionary | palatalize 1.000 | palatalize 1.000 |
| Romanian | R17a dorsal -i/-e [with_nde] (nested: velar_front_lemma) | ge1 | palatalize 1.000 | palatalize 1.000 |
| Romanian | R17b dorsal -i | dictionary | palatalize 1.000 | palatalize 1.000 |
| Romanian | R17b dorsal -i | ge1 | palatalize 1.000 | palatalize 1.000 |
| Romanian | R17b dorsal -e | dictionary | palatalize 1.000 | palatalize 1.000 |
| Romanian | R17b dorsal -e | ge1 | palatalize 1.000 | palatalize 1.000 |
| Romanian | R17c coronal -i (nested: palatal_context) | dictionary | palatalize 1.000 | palatalize 1.000 |
| Romanian | R17c coronal -i (nested: palatal_context) | ge1 | palatalize 1.000 | palatalize 1.000 |
| Romanian | R17c coronal -e | dictionary | retain 1.000 | retain 1.000 |
| Romanian | R17c coronal -e | ge1 | retain 1.000 | retain 1.000 |
| Romanian | R17c coronal -i/-e (nested: e_plural) | dictionary | conditioned 0.980, neither 0.020 | conditioned 1.000 |
| Romanian | R17c coronal -i/-e (nested: e_plural) | ge1 | conditioned 0.985, neither 0.015 | conditioned 1.000 |
| Romanian | R17c z -i (nested: borrowed) | dictionary | retain 0.688, neither 0.308 | retain 1.000 |
| Romanian | R17c z -i (nested: borrowed) | ge1 | retain 0.632, neither 0.368 | retain 1.000 |
| Romanian | R17c z -e | dictionary | retain 1.000 | retain 1.000 |
| Romanian | R17c z -e | ge1 | retain 1.000 | retain 1.000 |
| Italian | R17a dorsal -i/-e [with_nde] (nested: velar_front_lemma) | dictionary | neither 0.506, retain 0.494 | neither 1.000 |
| Italian | R17a dorsal -i/-e [with_nde] (nested: velar_front_lemma) | ge1 | neither 0.548, retain 0.452 | neither 1.000 |
| Italian | R17b dorsal -i | dictionary | neither 0.864, retain 0.136 | neither 1.000 |
| Italian | R17b dorsal -i | ge1 | neither 0.886, retain 0.114 | neither 1.000 |
| Italian | R17b dorsal -e | dictionary | retain 1.000 | retain 1.000 |
| Italian | R17b dorsal -e | ge1 | retain 1.000 | retain 1.000 |
| Italian | it_stress dorsal -i (nested: post_tonic) | dictionary | neither 0.848, retain 0.136, conditioned 0.016 | neither 1.000 |
| Italian | it_stress dorsal -i (nested: post_tonic) | ge1 | neither 0.546, conditioned 0.340, retain 0.114 | neither 1.000 |
| Italian | R17b dorsal -i [suffix_once] | dictionary | retain 0.955, neither 0.045 | neither 0.958, retain 0.042 |
| Italian | R17b dorsal -i [suffix_once] | ge1 | retain 0.949, neither 0.051 | neither 0.968, retain 0.032 |
| Italian | it_stress dorsal -i [suffix_once] (nested: post_tonic) | dictionary | retain 0.955, neither 0.044 | neither 0.958, retain 0.042 |
| Italian | it_stress dorsal -i [suffix_once] (nested: post_tonic) | ge1 | retain 0.949, neither 0.034, conditioned 0.017 | neither 0.968, retain 0.032 |
| Italian | it_stress dorsal -i [orth_ico_once] (nested: post_tonic) | dictionary | retain 0.997 | retain 1.000 |
| Italian | it_stress dorsal -i [orth_ico_once] (nested: post_tonic) | ge1 | retain 0.998 | retain 1.000 |

### Table 3. Language-specific model fits

Nested rule sets, Romanian (Yang's order; each branch: N, e (θ), ✓ = tolerated):

| rule set | τ | evidence | branches | grammar |
|---|---|---|---|---|
| coronal rule (t, d, s) | dictionary | E0 | -e plural → retain: 3,246, e 14 (θ 401.5) ✓; elsewhere (-i) → palatalize: 2,898, e 124 (θ 363.5) ✓ | **conditioned** |
| coronal rule (t, d, s) | ge1 | E0 | -e plural → retain: 1,087, e 4 (θ 155.5) ✓; elsewhere (-i) → palatalize: 871, e 10 (θ 128.7) ✓ | **conditioned** |

Nested rule sets, Italian (Yang's order; each branch: N, e (θ), ✓ = tolerated):

| rule set | τ | evidence | branches | grammar |
|---|---|---|---|---|
| stress-nested -i (giavazzi) | dictionary | E0 | velar post-tonic → retain: 845, e 81 (θ 125.4) ✓; elsewhere → palatalize: 3,405, e 91 (θ 418.7) ✓ | **conditioned** |
| stress-nested -i (suffix_once) | dictionary | E0 | velar post-tonic → retain: 818, e 54 (θ 122.0) ✓; elsewhere → palatalize: 1,851, e 91 (θ 246.0) ✓ | **conditioned** |
| stress-nested -i (orth_ico_once) | dictionary | E0 | velar post-tonic → retain: 799, e 35 (θ 119.6) ✓; elsewhere → palatalize: 275, e 64 (θ 49.0) ✗ | not conditioned |
| stress-nested -i (giavazzi) | ge1 | E0 | velar post-tonic → retain: 299, e 22 (θ 52.5) ✓; elsewhere → palatalize: 901, e 32 (θ 132.4) ✓ | **conditioned** |
| stress-nested -i (suffix_once) | ge1 | E0 | velar post-tonic → retain: 292, e 15 (θ 51.4) ✓; elsewhere → palatalize: 449, e 32 (θ 73.5) ✓ | **conditioned** |
| stress-nested -i (orth_ico_once) | ge1 | E0 | velar post-tonic → retain: 286, e 9 (θ 50.6) ✓; elsewhere → palatalize: 104, e 21 (θ 22.4) ✓ | **conditioned** |
| stress-nested -i (giavazzi) | dictionary_explicit | E0 | velar post-tonic → retain: 70, e 25 (θ 16.5) ✗; elsewhere → palatalize: 109, e 57 (θ 23.2) ✗ | not conditioned |
| stress-nested -i (suffix_once) | dictionary_explicit | E0 | velar post-tonic → retain: 70, e 25 (θ 16.5) ✗; elsewhere → palatalize: 109, e 57 (θ 23.2) ✗ | not conditioned |
| stress-nested -i (orth_ico_once) | dictionary_explicit | E0 | velar post-tonic → retain: 70, e 25 (θ 16.5) ✗; elsewhere → palatalize: 82, e 30 (θ 18.6) ✗ | not conditioned |
| category_stress (noun + verb) | dictionary | E1 | -are verbs → retain: 940, e 0 (θ 137.3) ✓; -isc- verbs → retain: 36, e 0 (θ 10.0) ✓; post-tonic nouns → retain: 845, e 81 (θ 125.4) ✓; elsewhere → palatalize: 6,048, e 291 (θ 694.6) ✓ | **conditioned** |
| category_stress (noun + verb) | ge1 | E1 | -are verbs → retain: 227, e 0 (θ 41.8) ✓; -isc- verbs → retain: 3, e 0 (θ 2.7) ✓; post-tonic nouns → retain: 299, e 22 (θ 52.5) ✓; elsewhere → palatalize: 1,068, e 45 (θ 153.2) ✓ | **conditioned** |

Nested rule sets, Czech (Yang's order; each branch: N, e (θ), ✓ = tolerated):

| rule set | τ | evidence | branches | grammar |
|---|---|---|---|---|
| NOM.PL m.an nested | dictionary | E0 | -log/-urg/-graf → retain: 109, e 3 (θ 23.2) ✓; -ský/-cký → palatalize: 4, e 0 (θ 2.9) ✓; other R26 suffixes (-ák, -ík, -ník, -ek) → palatalize: 779, e 214 (θ 117.0) ✗; monosyllables → palatalize: 62, e 30 (θ 15.0) ✗; loans → palatalize: 37, e 12 (θ 10.2) ✗; elsewhere → palatalize: 458, e 139 (θ 74.8) ✗ | not conditioned |
| NOM.PL m.an nested | ge1 | E0 | -log/-urg/-graf → retain: 44, e 0 (θ 11.6) ✓; -ský/-cký → palatalize: 3, e 0 (θ 2.7) ✓; other R26 suffixes (-ák, -ík, -ník, -ek) → palatalize: 283, e 20 (θ 50.1) ✓; monosyllables → palatalize: 37, e 12 (θ 10.2) ✗; loans → palatalize: 16, e 4 (θ 5.8) ✓; elsewhere → palatalize: 150, e 20 (θ 29.9) ✓ | not conditioned |

Single rules and pooled rows (a two-rule row reads N: e (pal) / e (ret), θ → productive). In a Czech competition row *palatalize* is the rule 'take the palatalizing allomorph' (e = roots taking the competitor) and *retain* the rule 'take the competitor'; the baseline is the hard labial, l, s and z stems, where no alternation is involved:

| language | rule | target | τ | evidence | N, e (θ) | shielded | verdict |
|---|---|---|---|---|---|---|---|
| Romanian | dorsal rule: dorsal/_[+front V] → palatalize | dorsal | dictionary | E0 | 2,874, e 180 (θ 360.9) ✓ | 501 | **palatalize** |
| Romanian | single coronal rule over -i and -e (t, d, s) | t, d, s | dictionary | E0 | 6,144: 3,356 / 2,788, θ 704.3 → neither | 1,138 | neither |
| Romanian | z before -i → retain (prespecified /z/) | z | dictionary | E0 | 126, e 10 (θ 26.1) ✓ | 941 | **retain** |
| Romanian | dorsal rule: dorsal/_[+front V] → palatalize | dorsal | ge1 | E0 | 601, e 39 (θ 93.9) ✓ | 93 | **palatalize** |
| Romanian | single coronal rule over -i and -e (t, d, s) | t, d, s | ge1 | E0 | 1,958: 1,093 / 865, θ 258.3 → neither | 234 | neither |
| Romanian | z before -i → retain (prespecified /z/) | z | ge1 | E0 | 48, e 2 (θ 12.4) ✓ | 171 | **retain** |
| Italian | pooled noun + verb units, -i | dorsal | dictionary | E1 | 7,869: 2,037 / 5,852, θ 877.2 → neither | 936 | neither |
| Italian | verbs -are (2sg, 3sg sbjv, 1pl) | dorsal | dictionary | E1 | 940: 940 / 0, θ 137.3 → **retain** | 0 | **retain** |
| Italian | verbs root-stressed -ere | dorsal | dictionary | E1 | 249: 0 / 249, θ 45.1 → **palatalize** | 3 | **palatalize** |
| Italian | verbs -isc- | dorsal | dictionary | E1 | 36: 36 / 0, θ 10.0 → **retain** | 0 | **retain** |
| Italian | all velar-stem verbs | dorsal | dictionary | E1 | 1,257: 976 / 281, θ 176.1 → neither | 3 | neither |
| Italian | pooled noun + verb units, -i | dorsal | ge1 | E1 | 1,597: 552 / 1,045, θ 216.5 → neither | 810 | neither |
| Italian | verbs -are (2sg, 3sg sbjv, 1pl) | dorsal | ge1 | E1 | 227: 227 / 0, θ 41.8 → **retain** | 196 | **retain** |
| Italian | verbs root-stressed -ere | dorsal | ge1 | E1 | 57: 0 / 57, θ 14.1 → **palatalize** | 61 | **palatalize** |
| Italian | verbs -isc- | dorsal | ge1 | E1 | 3: 3 / 0, θ 2.7 → **retain** | 0 | **retain** |
| Italian | all velar-stem verbs | dorsal | ge1 | E1 | 295: 230 / 65, θ 51.9 → neither | 270 | neither |
| Czech | rule TP (a): velar before NOM.PL m.an -i | velar | dictionary | E0 | 1,272, e 2 (θ 177.9) ✓ | 206 | **palatalize** |
| Czech | rule TP (a): velar before DAT/LOC.SG f -e | velar | dictionary | E0 | 3,786, e 1 (θ 459.5) ✓ | 0 | **palatalize** |
| Czech | rule TP (a): velar before LOC.PL m -ích | velar | dictionary | E0 | 3,480, e 0 (θ 426.7) ✓ | 52 | **palatalize** |
| Czech | rule TP (a): velar before NOM.PL.M.AN -í | velar | dictionary | E0 | 2,813, e 12 (θ 354.2) ✓ | 0 | **palatalize** |
| Czech | rule TP (a): r before NOM.PL m.an -i | r | dictionary | E0 | 367, e 1 (θ 62.1) ✓ | 37 | **palatalize** |
| Czech | rule TP (a): r before DAT/LOC.SG f -e | r | dictionary | E0 | 328, e 0 (θ 56.6) ✓ | 0 | **palatalize** |
| Czech | competition NOM.PL m.an: -i vs -ové/-é | velar | dictionary | E0 | 1,449: 504 / 1,235, θ 199.1 → neither | 0 | neither |
| Czech | competition NOM.PL m.an: -i vs -ové/-é | baseline | dictionary | E0 | 413: 363 / 149, θ 68.6 → neither | 0 | neither |
| Czech | competition LOC.SG m.in: -e vs -u | velar | dictionary | E0 | 2,045: 2,039 / 16, θ 268.3 → **retain** | 0 | **retain** |
| Czech | competition LOC.SG m.in: -e vs -u | baseline | dictionary | E0 | 2,237: 2,210 / 146, θ 290.0 → **retain** | 0 | **retain** |
| Czech | competition LOC.SG n: -e vs -u | velar | dictionary | E0 | 402: 399 / 6, θ 67.0 → **retain** | 0 | **retain** |
| Czech | competition LOC.SG n: -e vs -u | baseline | dictionary | E0 | 450: 378 / 209, θ 73.7 → neither | 0 | neither |
| Czech | competition LOC.PL m.in: -ích vs -ách/-ech (baseline: -ech vs -ách) | velar | dictionary | E0 | 2,070: 351 / 2,011, θ 271.1 → neither | 0 | neither |
| Czech | competition LOC.PL m.in: -ích vs -ách/-ech (baseline: -ech vs -ách) | baseline | dictionary | E0 | 2,220: 0 / 2,220, θ 288.1 → **palatalize** | 0 | **palatalize** |
| Czech | competition VOC.SG m.an: -e vs -u | velar | dictionary | E0 | 1,425: 1,415 / 13, θ 196.2 → **retain** | 0 | **retain** |
| Czech | competition VOC.SG m.an: -e vs -u | baseline | dictionary | E0 | 391: 5 / 386, θ 65.5 → **palatalize** | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -it | velar | dictionary | E2 | 52, e 0 (θ 13.2) ✓ | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -et | velar | dictionary | E2 | 7, e 1 (θ 3.6) ✓ | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -ní/-ný | velar | dictionary | E2 | 131, e 4 (θ 26.9) ✓ | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -ský/-cký | velar | dictionary | E2 | 303, e 3 (θ 53.0) ✓ | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -ina | velar | dictionary | E2 | 26, e 0 (θ 8.0) ✓ | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -í | velar | dictionary | E2 | 100, e 0 (θ 21.7) ✓ | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -ek | velar | dictionary | E2 | 482, e 0 (θ 78.0) ✓ | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -ka | velar | dictionary | E2 | 370, e 0 (θ 62.6) ✓ | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -ík | velar | dictionary | E2 | 17, e 0 (θ 6.0) ✓ | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -an | velar | dictionary | E2 | 9, e 3 (θ 4.1) ✓ | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -ník | velar | dictionary | E2 | 34, e 0 (θ 9.6) ✓ | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -nice | velar | dictionary | E2 | 11, e 0 (θ 4.6) ✓ | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -ista | velar | dictionary | E2 | 4, e 4 (θ 2.9) ✗ | 0 | non-trigger suffix: N = e |
| Czech | R22 derivation: velar base + -ismus | velar | dictionary | E2 | 3, e 3 (θ 2.7) ✗ | 0 | non-trigger suffix: N = e |
| Czech | rule TP (a): velar before NOM.PL m.an -i | velar | ge1 | E0 | 442, e 0 (θ 72.6) ✓ | 99 | **palatalize** |
| Czech | rule TP (a): velar before DAT/LOC.SG f -e | velar | ge1 | E0 | 1,124, e 0 (θ 160.0) ✓ | 0 | **palatalize** |
| Czech | rule TP (a): velar before LOC.PL m -ích | velar | ge1 | E0 | 484, e 0 (θ 78.3) ✓ | 13 | **palatalize** |
| Czech | rule TP (a): velar before NOM.PL.M.AN -í | velar | ge1 | E0 | 456, e 0 (θ 74.5) ✓ | 0 | **palatalize** |
| Czech | rule TP (a): r before NOM.PL m.an -i | r | ge1 | E0 | 184, e 0 (θ 35.3) ✓ | 11 | **palatalize** |
| Czech | rule TP (a): r before DAT/LOC.SG f -e | r | ge1 | E0 | 131, e 0 (θ 26.9) ✓ | 0 | **palatalize** |
| Czech | competition NOM.PL m.an: -i vs -ové/-é | velar | ge1 | E0 | 533: 100 / 436, θ 84.9 → neither | 0 | neither |
| Czech | competition NOM.PL m.an: -i vs -ové/-é | baseline | ge1 | E0 | 124: 92 / 34, θ 25.7 → neither | 0 | neither |
| Czech | competition LOC.SG m.in: -e vs -u | velar | ge1 | E0 | 1,007: 1,002 / 5, θ 145.6 → **retain** | 0 | **retain** |
| Czech | competition LOC.SG m.in: -e vs -u | baseline | ge1 | E0 | 972: 946 / 29, θ 141.3 → **retain** | 0 | **retain** |
| Czech | competition LOC.SG n: -e vs -u | velar | ge1 | E0 | 126: 123 / 3, θ 26.1 → **retain** | 0 | **retain** |
| Czech | competition LOC.SG n: -e vs -u | baseline | ge1 | E0 | 134: 82 / 54, θ 27.4 → neither | 0 | neither |
| Czech | competition LOC.PL m.in: -ích vs -ách/-ech (baseline: -ech vs -ách) | velar | ge1 | E0 | 372: 14 / 358, θ 62.8 → **palatalize** | 0 | **palatalize** |
| Czech | competition LOC.PL m.in: -ích vs -ách/-ech (baseline: -ech vs -ách) | baseline | ge1 | E0 | 381: 0 / 381, θ 64.1 → **palatalize** | 0 | **palatalize** |
| Czech | competition VOC.SG m.an: -e vs -u | velar | ge1 | E0 | 166: 165 / 1, θ 32.5 → **retain** | 0 | **retain** |
| Czech | competition VOC.SG m.an: -e vs -u | baseline | ge1 | E0 | 41: 0 / 41, θ 11.0 → **palatalize** | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -it | velar | ge1 | E2 | 40, e 0 (θ 10.8) ✓ | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -et | velar | ge1 | E2 | 5, e 0 (θ 3.1) ✓ | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -ní/-ný | velar | ge1 | E2 | 116, e 3 (θ 24.4) ✓ | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -ský/-cký | velar | ge1 | E2 | 188, e 2 (θ 35.9) ✓ | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -ina | velar | ge1 | E2 | 17, e 0 (θ 6.0) ✓ | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -í | velar | ge1 | E2 | 53, e 0 (θ 13.3) ✓ | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -ek | velar | ge1 | E2 | 165, e 0 (θ 32.3) ✓ | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -ka | velar | ge1 | E2 | 169, e 0 (θ 32.9) ✓ | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -ík | velar | ge1 | E2 | 6, e 0 (θ 3.3) ✓ | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -an | velar | ge1 | E2 | 5, e 1 (θ 3.1) ✓ | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -ník | velar | ge1 | E2 | 21, e 0 (θ 6.9) ✓ | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -nice | velar | ge1 | E2 | 7, e 0 (θ 3.6) ✓ | 0 | **palatalize** |
| Czech | R22 derivation: velar base + -ista | velar | ge1 | E2 | 4, e 4 (θ 2.9) ✗ | 0 | non-trigger suffix: N = e |
| Czech | R22 derivation: velar base + -ismus | velar | ge1 | E2 | 2, e 2 (θ 2.9) ✓ | 0 | non-trigger suffix: N = e |

### Table 4. Spell-out on certain-specification roots (R20)

e (pal) = roots taking only the non-palatalizing member (doublets by the dominant variant); e (ret) = roots taking the palatalizing member. Contrast rows: roots without alternation (prespecified for the draft).

| language | pair | target | roots | τ = dictionary: N: e (pal) / e (ret), θ → productive | τ ≥ 1 token |
|---|---|---|---|---|---|
| Romanian | -i vs -a/-ui verbalizer, alternating-plural roots | dorsal | all | 128: 51 / 78, θ 26.4 → neither | 77: 33 / 44, θ 17.7 → neither |
| Romanian | -i vs -a/-ui verbalizer, alternating-plural roots | t, d, s | all | 57: 34 / 23, θ 14.1 → neither | 36: 21 / 15, θ 10.0 → neither |
| Romanian | contrast: non-alternating-plural roots | dorsal | contrast:nonalt | 39 roots: non-palatalizing member only 28, -i only 9, both 2 | – |
| Romanian | contrast: non-alternating-plural roots | t, d, s | contrast:nonalt | 110 roots: non-palatalizing member only 80, -i only 27, both 3 | – |
| Italian | -izzare vs -eggiare/-ire/-are | dorsal | all | 44: 25 / 20, θ 11.6 → neither | 17: 11 / 6, θ 6.0 → **retain** |
| Italian | contrast: -izzare vs -eggiare/-ire/-are | dorsal | contrast:nonalt | 88: 87 / 1, θ 19.7 → **retain** | 39: 38 / 1, θ 10.6 → **retain** |
| Italian | -izzare vs -eggiare/-ire | dorsal | all | 26: 7 / 20, θ 8.0 → **palatalize** | 8: 2 / 6, θ 3.8 → **palatalize** |
| Italian | contrast: -izzare vs -eggiare/-ire | dorsal | contrast:nonalt | 27: 25 / 2, θ 8.2 → **retain** | 8: 7 / 1, θ 3.8 → **retain** |
| Italian | palatal vs velar -ista/-ismo/-ità | dorsal | all | 231: 0 / 231, θ 42.4 → **palatalize** | 73: 0 / 73, θ 17.0 → **palatalize** |
| Italian | contrast: palatal vs velar -ista/-ismo/-ità | dorsal | contrast:nonalt | 23: 18 / 5, θ 7.3 → **retain** | 3: 2 / 1, θ 2.7 → **both** |
| Czech | -ní/-ný vs -ový/-ovní | velar | feminine roots; all | 301: 273 / 30, θ 52.7 → **retain** | 109: 85 / 24, θ 23.2 → neither |
| Czech | -ní/-ný vs -ový/-ovní | velar | feminine roots; base:root-final | 153: 125 / 30, θ 30.4 → **retain** | 59: 35 / 24, θ 14.5 → neither |
| Czech | -ní/-ný vs -ový/-ovní | velar | feminine roots; base:suffix | 148: 148 / 0, θ 29.6 → **retain** | 50: 50 / 0, θ 12.8 → **retain** |
| Czech | -it vs -ovat | velar | feminine roots; all | 34: 22 / 12, θ 9.6 → neither | 20: 12 / 8, θ 6.7 → neither |
| Czech | -it vs -ovat | velar | feminine roots; base:root-final | 22: 11 / 11, θ 7.1 → neither | 12: 5 / 7, θ 4.8 → neither |
| Czech | -it vs -ovat | velar | feminine roots; base:suffix | 12: 11 / 1, θ 4.8 → **retain** | 8: 7 / 1, θ 3.8 → **retain** |
| Czech | -ka vs -yně | velar | feminine roots; all | 0: 0 / 0, θ – → – | 0: 0 / 0, θ – → – |
| Czech | -ka vs -yně | velar | E0-certain roots; all | 250: 19 / 233, θ 45.3 → **palatalize** | 77: 10 / 67, θ 17.7 → **palatalize** |
| Czech | any pair | velar | feminine roots; all | 321: 287 / 36, θ 55.6 → **retain** | 123: 94 / 29, θ 25.6 → neither |

### Table 5. No alternation evidence (R21)

Roots by base class (alternation evidence in inflection, E0):

| language | target | τ | alternating | non-alternating | no evidence |
|---|---|---|---|---|---|
| Romanian | all | dictionary | 5,492 | 2,136 | 8,106 (52%) |
| Romanian | all | ge1 | 1,429 | 438 | 13,867 (88%) |
| Italian | dorsal | dictionary | 5,575 | 1,037 | 1,120 (14%) |
| Italian | dorsal | ge1 | 980 | 322 | 6,430 (83%) |
| Czech | velar | dictionary | 7,330 | 447 | 518 (6%) |
| Czech | r | dictionary | 730 | 747 | 135 (8%) |
| Czech | velar | ge1 | 1,973 | 912 | 5,410 (65%) |
| Czech | r | ge1 | 329 | 360 | 923 (57%) |
| Czech | velar | ge1_strict | 1,758 | 130 | 6,407 (77%) |
| Czech | r | ge1_strict | 302 | 11 | 1,299 (81%) |
| Czech | velar | ge2pm | 280 | 250 | 7,765 (94%) |
| Czech | r | ge2pm | 79 | 114 | 1,419 (88%) |

Derivative units by base class (palatalized / not, % palatalized):

| language | units | target | τ | alternating | non-alternating | no evidence | both | back-vowel suffix instead (alt. / non-alt. / no ev.) |
|---|---|---|---|---|---|---|---|---|
| Romanian | front-vowel derivatives (verbs -i, -iza/-ifica, adjectives) | dorsal | dictionary | 129 / 0 (100%) | 10 / 0 (100%) | 11 / 0 (100%) | 0 | 92 / 46 / 19 |
| Romanian | front-vowel derivatives (verbs -i, -iza/-ifica, adjectives) | t, d, s | dictionary | 20 / 64 (24%) | 3 / 36 (8%) | 9 / 92 (9%) | 0 | 56 / 96 / 251 |
| Romanian | front-vowel derivatives (verbs -i, -iza/-ifica, adjectives) | dorsal | ge1 | 59 / 0 (100%) | 2 / 0 (100%) | 17 / 0 (100%) | 0 | 46 / 19 / 23 |
| Romanian | front-vowel derivatives (verbs -i, -iza/-ifica, adjectives) | t, d, s | ge1 | 7 / 22 (24%) | 0 / 16 (0%) | 5 / 73 (6%) | 0 | 28 / 42 / 222 |
| Italian | learned set | dorsal | dictionary | 280 / 0 (100%) | 9 / 40 (18%) | 10 / 4 (71%) | 0 | – |
| Italian | learned set, internal formations | dorsal | dictionary | 258 / 0 (100%) | 5 / 23 (18%) | 9 / 3 (75%) | 0 | – |
| Italian | learned set, with generated candidates | dorsal | dictionary | 322 / 1 (100%) | 19 / 47 (29%) | 12 / 4 (75%) | 2 | – |
| Italian | native control suffixes | dorsal | dictionary | 1 / 31 (3%) | 4 / 120 (3%) | 6 / 77 (7%) | 0 | – |
| Italian | learned set | dorsal | ge1 | 89 / 0 (100%) | 3 / 14 (18%) | 8 / 1 (89%) | 0 | – |
| Italian | learned set, internal formations | dorsal | ge1 | 78 / 0 (100%) | 1 / 6 (14%) | 6 / 1 (86%) | 0 | – |
| Italian | learned set, with generated candidates | dorsal | ge1 | 120 / 1 (99%) | 4 / 25 (14%) | 17 / 3 (85%) | 0 | – |
| Italian | native control suffixes | dorsal | ge1 | 1 / 8 (11%) | 3 / 58 (5%) | 4 / 41 (9%) | 0 | – |
| Czech | P suffixes (velar bases) | velar | dictionary | 1,472 / 8 (99%) | 32 / 2 (94%) | 27 / 0 (100%) | 1 | – |
| Czech | P suffixes, without learned -ik(a) → -ický | velar | dictionary | 1,340 / 8 (99%) | 16 / 2 (89%) | 27 / 0 (100%) | 1 | – |
| Czech | P suffixes, derivatives internal formations only | velar | dictionary | 357 / 4 (99%) | 3 / 1 (75%) | 5 / 0 (100%) | 0 | – |
| Czech | P suffixes (r bases) | r | dictionary | 75 / 3 (96%) | 28 / 1 (97%) | 2 / 0 (100%) | 0 | – |
| Czech | P suffixes (velar bases) | velar | ge1 | 480 / 3 (99%) | 151 / 2 (99%) | 155 / 1 (99%) | 0 | – |
| Czech | P suffixes, without learned -ik(a) → -ický | velar | ge1 | 432 / 3 (99%) | 148 / 2 (99%) | 113 / 1 (99%) | 0 | – |
| Czech | P suffixes, derivatives internal formations only | velar | ge1 | 172 / 3 (98%) | 69 / 1 (99%) | 33 / 1 (97%) | 0 | – |
| Czech | P suffixes (r bases) | r | ge1 | 31 / 0 (100%) | 13 / 1 (93%) | 12 / 2 (86%) | 0 | – |
| Czech | P suffixes (velar bases) | velar | ge1_strict | 429 / 2 (100%) | 41 / 0 (100%) | 316 / 4 (99%) | 0 | – |
| Czech | P suffixes, without learned -ik(a) → -ický | velar | ge1_strict | 382 / 2 (99%) | 39 / 0 (100%) | 272 / 4 (99%) | 0 | – |
| Czech | P suffixes, derivatives internal formations only | velar | ge1_strict | 158 / 2 (99%) | 22 / 0 (100%) | 94 / 3 (97%) | 0 | – |
| Czech | P suffixes (r bases) | r | ge1_strict | 23 / 0 (100%) | 0 / 1 (0%) | 33 / 2 (94%) | 0 | – |
| Czech | P suffixes (velar bases) | velar | ge2pm | 75 / 0 (100%) | 62 / 4 (94%) | 113 / 1 (99%) | 0 | – |
| Czech | P suffixes, without learned -ik(a) → -ický | velar | ge2pm | 64 / 0 (100%) | 62 / 4 (94%) | 74 / 1 (99%) | 0 | – |
| Czech | P suffixes, derivatives internal formations only | velar | ge2pm | 38 / 0 (100%) | 35 / 4 (90%) | 36 / 0 (100%) | 0 | – |
| Czech | P suffixes (r bases) | r | ge2pm | 4 / 0 (100%) | 3 / 0 (100%) | 6 / 0 (100%) | 0 | – |

Italian no-evidence units (learned set, named links) by what the base is:

| base | τ | palatalized / velar |
|---|---|---|
| -ic- suffix base | dictionary | 7 / 0 (100%) |
| masculine counterpart (other root) | dictionary | 5 / 0 (100%) |
| masculine counterpart (same root) | dictionary | 0 / 0 |
| neither | dictionary | 3 / 4 (43%) |
| feminine-only (-e plurals) | dictionary | 7 / 4 (64%) |
| no plural | dictionary | 3 / 0 (100%) |
| -ic- suffix base | ge1 | 7 / 0 (100%) |
| masculine counterpart (other root) | ge1 | 2 / 0 (100%) |
| masculine counterpart (same root) | ge1 | 0 / 0 |
| neither | ge1 | 1 / 1 (50%) |
| -i plural below ge1 | ge1 | 5 / 0 (100%) |
| feminine-only (-e plurals) | ge1 | 3 / 1 (75%) |

### Table 6. Co-variation of inflectional class and derivative choice (R23)

P base = the base alternates (takes the palatalizing allomorph); P/N = the derivative palatalizes (takes the palatalizing member) or not; units with both outcomes left out. OR and Fisher p of the 2×2 (summed over pairs or suffixes); MH = Mantel–Haenszel over each language's strata; perm. p = stratified permutation (10,000 shuffles); logit = base_class in deriv_pal ~ base_class + log_lemma_freq + base_stratum + deriv_stratum, base-clustered errors.

| language | test | base | τ | P base: P / N | N base: P / N | OR [95% CI] | Fisher p | MH OR [95% CI] (strata) | perm. p | logit base_class |
|---|---|---|---|---|---|---|---|---|---|---|
| Romanian | plural class × verbalizer (-i vs -a/-ui) | dorsal | dictionary | 72 / 42 | 9 / 28 | 5.33 [2.3, 12.4] | 5.01e-05 | 5.34 [2.29, 12.5] (3) | < 0.0001 | 1.57 (se 0.459, p 6.41e-04) |
| Romanian | plural class × verbalizer (-i vs -a/-ui) | dorsal | ge1 | 41 / 25 | 2 / 13 | 10.7 [2.22, 51.2] | 0.00104 | 13.7 [2.51, 74.5] (3) | 0.0003 | 2.25 (se 0.907, p 0.0131) |
| Romanian | plural class × verbalizer (-i vs -a/-ui) | t, d, s | dictionary | 22 / 33 | 27 / 80 | 1.98 [0.987, 3.95] | 0.0705 | 2.13 [1.03, 4.41] (3) | 0.0188 | 0.745 (se 0.374, p 0.0465) |
| Romanian | plural class × verbalizer (-i vs -a/-ui) | t, d, s | ge1 | 14 / 20 | 12 / 42 | 2.45 [0.96, 6.25] | 0.0919 | 2.31 [0.874, 6.11] (3) | 0.048 | 0.824 (se 0.516, p 0.11) |
| Romanian | plural class × verbalizer (-i vs -a/-ui) | all | dictionary | 95 / 75 | 45 / 119 | 3.35 [2.12, 5.29] | 1.57e-07 | 3.3 [1.93, 5.64] (9) | < 0.0001 | 1.16 (se 0.244, p 2.07e-06) |
| Romanian | plural class × verbalizer (-i vs -a/-ui) | all | ge1 | 55 / 45 | 18 / 59 | 4.01 [2.07, 7.74] | 2.78e-05 | 4.09 [1.85, 9.06] (8) | < 0.0001 | 1.19 (se 0.369, p 0.00122) |
| Italian | learned suffixes, internal formations | dorsal | dictionary | 258 / 0 | 5 / 22 | ∞ [–, ∞] | 2.05e-28 | ∞ (6) | < 0.0001 | 45.7 (separated) |
| Italian | learned suffixes, internal formations | dorsal | ge1 | 78 / 0 | 1 / 5 | ∞ [–, ∞] | 1.94e-07 | ∞ (6) | 0.0009 | 122 (separated) |
| Italian | learned suffixes, internal formations, no -ic- bases | dorsal | dictionary | 7 / 0 | 5 / 22 | ∞ [–, ∞] | 1.47e-04 | ∞ (6) | – | – |
| Italian | learned suffixes, internal formations, no -ic- bases | dorsal | ge1 | 1 / 0 | 1 / 5 | ∞ [–, ∞] | 0.286 | ∞ (6) | – | – |
| Italian | native control suffixes, internal formations | dorsal | dictionary | 0 / 27 | 2 / 102 | 0 [0, –] | 1 | 0 (11) | 1 | -20 (separated) |
| Italian | native control suffixes, internal formations | dorsal | ge1 | 0 / 6 | 2 / 52 | 0 [0, –] | 1 | – | – | -23 (separated) |
| Czech | LOC.SG class (any) × suffix choice (all pairs) | velar | dictionary | 4 / 6 | 68 / 384 | 3.76 [1.04, 13.7] | 0.0544 | 1.96 [0.5, 7.74] (8) | 0.22 | -0.282 (se 1.05, p 0.789) |
| Czech | LOC.SG class (any) × suffix choice (all pairs) | velar | ge1 | 4 / 2 | 60 / 164 | 5.47 [0.976, 30.6] | 0.0522 | 2.25 [0.44, 11.4] (8) | 0.25 | 0.274 (se 1.01, p 0.785) |
| Czech | LOC.SG class (any) × suffix choice (all pairs) | baseline | dictionary | 96 / 105 | 125 / 258 | 1.89 [1.33, 2.68] | 4.48e-04 | 1.63 [1.12, 2.36] (8) | 0.01 | – |
| Czech | LOC.SG class (any) × suffix choice (all pairs) | baseline | ge1 | 40 / 41 | 112 / 160 | 1.39 [0.847, 2.29] | 0.203 | 1.21 [0.71, 2.06] (8) | 0.26 | – |
| Czech | NOM.PL class (dominant) × suffix choice (all pairs) | velar | dictionary | 175 / 21 | 62 / 8 | 1.08 [0.453, 2.55] | 0.827 | 1.15 [0.36, 3.6] (12) | 0.51 | – |
| Czech | NOM.PL class (dominant) × suffix choice (all pairs) | velar | ge1 | 69 / 9 | 31 / 6 | 1.48 [0.486, 4.53] | 0.557 | 1.48 [0.35, 6.24] (11) | 0.4 | – |
| Czech | NOM.PL class (any) × suffix choice (all pairs) | velar | dictionary | 192 / 29 | 57 / 4 | 0.465 [0.157, 1.38] | 0.183 | 0.78 [0.19, 3.15] (12) | 0.55 | -0.0495 (se 0.604, p 0.935) |
| Czech | NOM.PL class (any) × suffix choice (all pairs) | velar | ge1 | 74 / 11 | 30 / 5 | 1.12 [0.359, 3.5] | 1 | 2.11 [0.46, 9.61] (11) | 0.24 | 0.332 (se 0.633, p 0.6) |
| Czech | NOM.PL class (dominant) × suffix choice (all pairs) | baseline | dictionary | 3 / 9 | 85 / 4 | 0.0157 [0.00302, 0.0814] | 6.67e-08 | 0.03 [0, 0.46] (11) | 1 | – |
| Czech | NOM.PL class (dominant) × suffix choice (all pairs) | baseline | ge1 | 0 / 4 | 37 / 1 | 0 [0, –] | 4.47e-05 | – | – | – |

### Table 7. -log/-fag/-urg nouns (R27)

| language | nouns | cell | N | τ = dictionary: front (e pal / e ret), shielded | τ ≥ 1 token |
|---|---|---|---|---|---|
| Romanian | -log/-fag/-urg person | -i/-e | 194 | 194 before a front vowel (0 / 194), 0 shielded → **palatalize** | 45 before a front vowel (0 / 45), 0 shielded → **palatalize** |
| Romanian | -log/-fag/-urg thing | -i/-e | 14 | 8 before a front vowel (0 / 8), 6 shielded → **palatalize** | 5 before a front vowel (0 / 5), 2 shielded → **palatalize** |
| Italian | -logo/-fago/-urgo person | -i | 251 | 251 before a front vowel (14 / 242), 0 shielded → **palatalize** | 64 before a front vowel (5 / 59), 0 shielded → **palatalize** |
| Italian | -logo/-fago/-urgo thing | -i | 21 | 21 before a front vowel (10 / 11), 0 shielded → neither | 8 before a front vowel (5 / 3), 0 shielded → **retain** |
| Czech | -log/-urg person | NOM.PL m.an | 111 | 3 before a front vowel (1 / 2), 106 shielded → **both** | 0 before a front vowel (0 / 0), 43 shielded → – |
| Czech | -log/-urg person | LOC.PL | 111 | 109 before a front vowel (2 / 107), 0 shielded → **palatalize** | 2 before a front vowel (0 / 2), 0 shielded → **both** |
| Czech | -log/-urg thing | NOM.PL m.an | 9 | 0 before a front vowel (0 / 0), 9 shielded → – | 0 before a front vowel (0 / 0), 5 shielded → – |
| Czech | -log/-urg thing | LOC.PL | 9 | 9 before a front vowel (0 / 9), 0 shielded → **palatalize** | 4 before a front vowel (0 / 4), 0 shielded → **palatalize** |

### Table 8. Logical Phonology validators (share of target-stem rows whose derived surface form matches)

| language | policy | rows | boolean match | IPA match | policy |
|---|---|---|---|---|---|
| Romanian | class | 11,641 | 97.65% | 93.99% | every dorsal, and every coronal before -i, underspecified given the attested plural suffix (the draft's rules; legacy validator) |
| Italian | lexical | 7,338 | 99.52% | 99.52% | underspecified iff the root alternates (reads the -i outcome: circular for the -i cell) |
| Italian | stress | 7,338 | 97.63% | 97.63% | underspecified iff the velar is not post-tonic (no lexical listing) |
| Italian | all_underspecified | 7,338 | 85.28% | 85.28% | every dorsal underspecified |
| Czech | lexical | 48,143 | 99.98% | 99.95% | underspecified iff the lexeme alternates in E0 (the tested cell is its own evidence: circular) |
| Czech | all_underspecified | 48,143 | 99.92% | 99.89% | every target underspecified |
| Czech | lexical_loo | 48,143 | 52.39% | 52.39% | Inkelas learner: underspecified only with alternation evidence outside the tested cell |
| Czech | spellout | 48,143 | 90.63% | 90.60% | prespecified iff the lexeme takes a non-palatalizing ending somewhere (G&W's spell-out claim as a policy) |

### Table 9. Romanian legacy tables (the draft's productivity table and the full lexicon): N: e (pal) / e (ret) → productive

| stem final | -i, top-1000 (draft) | -i, full lexicon | -e, top-1000 (draft) | -e, full lexicon |
|---|---|---|---|---|
| <c> | 28: 0 / 28 → **palatalize** | 473: 0 / 473 → **palatalize** | 14: 0 / 14 → **palatalize** | 434: 0 / 434 → **palatalize** |
| <g> | 17: 0 / 17 → **palatalize** | 380: 0 / 380 → **palatalize** | 3: 0 / 3 → **palatalize** | 273: 0 / 273 → **palatalize** |
| <t> | 207: 0 / 207 → **palatalize** | 1,522: 0 / 1,522 → **palatalize** | 229: 229 / 0 → **retain** | 2,049: 2,048 / 1 → **retain** |
| <d> | 23: 0 / 23 → **palatalize** | 261: 0 / 261 → **palatalize** | 36: 36 / 0 → **retain** | 608: 608 / 0 → **retain** |
| <s> | 11: 0 / 11 → **palatalize** | 138: 0 / 138 → **palatalize** | 44: 44 / 0 → **retain** | 493: 493 / 0 → **retain** |
| <z> | 2: 1 / 1 → **both** | 14: 4 / 10 → **palatalize** | 48: 48 / 0 → **retain** | 862: 862 / 0 → **retain** |
| <st> | 4: 0 / 4 → **palatalize** | 26: 0 / 26 → **palatalize** | – | – |
| <ct> | 3: 0 / 3 → **palatalize** | 9: 0 / 9 → **palatalize** | 22: 22 / 0 → **retain** | 61: 61 / 0 → **retain** |
<!-- END crosslang -->

## The learnability paragraph, checked against three languages

> Beyond Romanian, I will extend the analysis to Italian and ask whether
> the underspecified class is learnable from the data available to
> learners in each language. Italian masculine plurals show the same
> lexical split for dorsals ([amiko, amitʃi] 'friend(s)' but [fwɔko,
> fwɔki] 'fire(s)') but no corresponding coronal alternation, and Kramer
> (2009) finds that speakers differ in how they extend this split to
> novel words. Since learners never observe underlying representations
> directly, each generation must reconstruct the underspecified class
> from surface alternations (cf. Hale & Reiss 2008). Following the
> underspecification criterion of Inkelas (1995), a learner
> underspecifies a segment only when it alternates and otherwise stores
> it fully specified. Thus, a learner can classify a stem only after
> observing it in an alternating environment. This account makes two
> predictions. First, nouns lacking the relevant plural, e.g., singularia
> tantum mass nouns and proper names, should be stored with prespecified
> obstruents and resist palatalization in derivatives. On this account,
> the pattern Steriade (2008) attributes to inflection dependence follows
> instead from the evidence available to the learner. Second, the
> Tolerance Principle (Yang 2005, 2016) supplies a psychologically
> plausible threshold for when a learner generalizes palatalization
> productively. Under the proposed underspecification analysis, the
> relevant exceptions are prespecified stems that occur before a
> palatalizing suffix without alternating. Allomorph selection therefore
> structures the learner's evidence. Because most prespecified dorsals
> select -uri and -a, they rarely enter the environments in which
> palatalization is evaluated, so palatalization faces few exceptions.
> Romanian counts show productive palatalization of dorsals before -i and
> -e and of t, d, and s before -i, motivating separate dorsal and coronal
> rules. If the Italian dorsal pattern fails the threshold, the Tolerance
> Principle predicts no single productive palatalization rule, consistent
> with the lack of a uniform generalization in the speaker variation
> reported by Kramer.

Every number below is a row of `crosslang_summary.csv` (the table
numbers point to the tables above). "Units" are base × suffix pairs;
"≥ 1" is τ ≥ 1 token. Four caveats apply throughout:

- the evidence is a news corpus of 15–20 million tokens per language,
  not a child's input;
- dictionary paradigms contain generated forms nobody uses (Italian
  module defaults, Czech template tables);
- the corpora are untagged, so syncretic forms (Czech LOC.SG -u)
  over-attest a cell;
- no number here rests on a CC BY-NC-SA source (MorfFlex, DeriNet,
  czetyl). The NC sections of the Czech reports (DeriNet pairs, czetyl
  strata) point the same way for R20 and R21 but are not used.

**The premise.** Italian masculine plurals do show the dorsal split with
no coronal counterpart. Masculine dorsal nouns divide almost evenly: 758
palatalize and 739 keep the velar at τ = dictionary, 256 and 264 at ≥ 1
(Table 1). No coronal plural changes (10,277 noun roots before -i). That
part is by construction, since the data layer codes every coronal plural
as unchanged, but it is the descriptive fact (gatto, gatti). On Krämer
see (g).

### (a) Learners reconstruct the underspecified class from alternations

*Claim:* a segment is underspecified only once it has been seen to
alternate; otherwise it is stored fully specified.

**Verdict: not supported** as a default for stems not yet seen
alternating (Romanian, Czech); **mixed** for Italian. The procedure
itself cannot be observed in corpus counts; its consequences can.

- **Few stems are ever seen alternating.** At ≥ 1 token (Table 5):
  - Romanian: 1,429 of 15,734 target roots alternate (9%) and 13,867
    (88%) show no evidence;
  - Italian: 980 of 7,732 dorsal roots alternate (13%) and 6,430 (83%)
    have no attested -i form;
  - Czech: 1,973 of 8,295 velar roots alternate (24%) and 5,410 (65%)
    show nothing.

  A learner with this much text classifies at most a quarter of the
  target stems as underspecified and stores the rest fully specified.
- **Stems stored that way do not behave as prespecified when they are
  used.** Czech gives the direct test (Table 8). Underspecifying only
  lexemes that alternate outside the tested cell (`lexical_loo`, the
  Inkelas learner) matches 52.39% of 48,143 target-stem forms, while
  underspecifying every target matches 99.92%. For most lexemes the
  palatalizing cell (feminine DAT/LOC.SG -e, adjective NOM.PL -í) is
  their only alternating environment, and it never fails: velar
  DAT/LOC.SG -e has 1,124 roots and 0 exceptions at ≥ 1 (Table 3).
  Romanian's validator, which underspecifies every dorsal and every
  coronal before -i, matches 97.65% of 11,641 rows. Bases without
  evidence also palatalize in derivatives (claim (b)).
- **Counting those stems as exceptions sinks every rule.** If every
  target-final root not seen alternating counted against the rule
  (`counterfactual_palatalize`, θ of N_target), no palatalization rule
  in any language would be productive. For Romanian dorsals before
  -i/-e, for example, that is 1,369 exceptions among 4,063 roots against
  θ = 488.9.
- **Italian needs lexical information, but not necessarily
  alternation.** Underspecifying every dorsal matches only 85.28% of
  7,338 rows: the fuoco type keeps its velar. Underspecifying by stress
  alone, with no lexical listing, matches 97.63%, and by attested
  alternation 99.52%; the last policy reads the outcome and is circular
  for the -i cell.

*Caveats:* the Czech test runs on dictionary forms; the lexical policies
are circular by construction; whether Italian learners use stress or
stored alternation is a question for nonce studies.

*Reading:* the data favour a learner who, once the rule is productive,
treats a stem not yet classified as an undergoer and stores only stems
seen not to alternate (in Italian, a rule conditioned on stress). That
is the reverse of the default the paragraph adopts.

### (b) Nouns without the relevant plural resist palatalization in derivatives

**Verdict: not supported** (Czech, Romanian dorsals); **mixed, small N**
(Italian); proper names cannot be tested, since names are out of the
domain in all three languages (R25). Counts are units palatalized of all
units (Table 5).

- **Czech**, velar bases, palatalizing suffixes, Wiktionary basis. Bases
  seen only with competing allomorphs and bases with no evidence
  palatalize as alternating bases do:

  | base class | τ = dictionary | ≥ 1 token | ≥ 1, strict |
  |---|---|---|---|
  | alternating | 1,472 of 1,480 | 480 of 483 | 429 of 431 |
  | competing allomorphs only | 32 of 34 | 151 of 153 | 41 of 41 |
  | no evidence | 27 of 28 | 155 of 156 | 316 of 320 |

  Derivatives that are Czech internal formations, coined in Czech and
  not stored from Proto-Slavic (the decisive test), give 33 of 34 for
  no-evidence bases at ≥ 1.
- **Romanian dorsals.** Front-vowel derivatives (verbs in -i,
  -iza/-ifica, front-vowel adjective suffixes) of no-evidence bases are
  palatalized 11 of 11 at τ = dictionary and 17 of 17 at ≥ 1, like
  alternating bases (129 of 129; 59 of 59) and -uri bases (10 of 10).
  What differs is the choice of suffix. At τ = dictionary, no-evidence
  bases take a back-vowel suffix (-a/-ui, back adjective suffixes) in 19
  of 30 units, alternating bases in 92 of 221 and non-alternating bases
  (almost all -uri) in 46 of 56.
- **Romanian coronals** do not separate the accounts. The derivatives
  rarely assibilate whatever the base: alternating 20 of 84, non-
  alternating 3 of 39, no evidence 9 of 101. Derived verbs are read from
  the infinitive and participle (`romanian_noevidence.txt`), and the
  infinitival -i does not assibilate, so these rows say little about the
  root's specification.
- **Italian**, learned suffixes, named links. No-evidence bases
  palatalize in 10 of 14 units at τ = dictionary and 8 of 9 at ≥ 1.
  That is between alternating bases (280 of 280) and bases with an
  attested velar plural (9 of 49).
  - Seven of the 10 palatalized units have an -ic- suffix base, five of
    them also a masculine counterpart (musica → musicista, logica →
    logicista).
  - The 7 units with neither divide 3 palatalized (elicità, imagismo,
    imagista) and 4 velar (cronachista, droghista, bisarchista,
    ribechista).
  - By why the evidence is missing: feminine-only -ca/-ga nouns 7 and 4;
    no plural at all 3 and 0.

*Caveats:*
- "No evidence" means no attested form in a news corpus (at τ =
  dictionary, no listed plural), not that a learner never heard one.
- Derivatives may be stored or borrowed whole; the Czech
  internal-formation rows address this, the Romanian strata cannot.
- The Italian units are few and mostly suffix-driven.
- Singularia tantum are a subset of the no-evidence classes (Italian
  "no plural": 3 of 3 palatalized).

### (c) The inflection-dependence pattern follows from the learner's evidence

**Verdict: mixed.** The co-variation Steriade describes is robust in
Romanian, suffix-driven in Italian and absent in Czech beyond what
non-alternating stems show too. Where it exists, it does not follow the
learner's alternation evidence (Table 6).

- **Romanian**, plural class × verbalizer -i vs -a/-ui:

  | | dorsal, dictionary | dorsal, ≥ 1 | coronal (t, d, s), dictionary |
  |---|---|---|---|
  | units, base P (P/N) and base N (P/N) | 72/42, 9/28 | 41/25, 2/13 | 22/33, 27/80 |
  | odds ratio [95% CI] | 5.33 [2.30, 12.4] | 10.7 [2.22, 51.2] | 1.98 [0.99, 3.95] |
  | Fisher p | 5.0e-05 | 0.001 | 0.07 |
  | Mantel–Haenszel OR | 5.34 [2.29, 12.5] | 13.7 | 2.13 [1.03, 4.41] |
  | permutation p | < 0.0001 | 0.0003 | 0.019 |
  | logistic base_class | 1.57 (se 0.459, p 6.4e-04) | 2.25 (p 0.013) | 0.745 (p 0.047) |

  The logistic model controls lemma frequency and stratum.
- **Italian**, learned suffixes, Italian-internal formations: 258/0
  against 5/22 (Fisher p 2.05e-28, permutation p < 0.0001). But 251 of
  the 258 palatal units have an -ic- base: plural and derivative both
  contain the suffix (abulico → abulicità).
  - Without -ic- bases: 7/0 against 5/22 (p 1.5e-04) at τ = dictionary;
    1/0 against 1/5 (p 0.29) at ≥ 1.
  - The logistic model is quasi-completely separated (no estimate).
  - Native suffixes never follow the base: 0/27 against 2/102.
- **Czech**, LOC.SG class × derivative suffix choice, velar.
  - The -ní/-ný vs -ový/-ovní pair alone gives OR 5.86 [1.52, 22.6]
    (Fisher p 0.018) from 4 + 5 units with a P base.
  - Stratified it does not hold: Mantel–Haenszel 1.96 [0.50, 7.74],
    permutation p 0.22, logistic base_class −0.28 (p 0.79).
  - The same association appears in hard stems with no alternation at
    all: Mantel–Haenszel 1.63 [1.12, 2.36], permutation p 0.01.
  - The NOM.PL class does not predict the choice: Mantel–Haenszel 1.15
    [0.36, 3.60], p 0.51 at τ = dictionary; 1.48, p 0.40 at ≥ 1.
- **Not the learner's evidence.** If the pattern came from the evidence,
  bases without evidence would pattern with non-alternating ones. They do
  not palatalize less (claim (b)). In Romanian they take back-vowel
  suffixes more often than alternating bases (63% against 42% of units),
  the one place where they lean toward the non-alternating class. Yet
  they palatalize every time they take a front suffix (11 of 11), which
  a prespecified root could not do.

*Caveats:*
- Few R23 units in Czech (4 + 5 P-base units for the pair) and Italian
  (27 units with a velar plural).
- Romanian strata are coarse (borrowed / latin_unresolved / unknown).
- The Romanian E2 reading is orthographic.

### (d) The threshold counts prespecified stems that meet a palatalizing suffix without alternating

**Verdict: supported as a way of counting.** It is the harmonized
definition (R16): N counts the roots attested before the trigger and e
those that do not alternate; roots whose attested forms all avoid the
trigger are shielded, not exceptions. Counted this way (Table 1):

| rule | language | τ = dictionary: N, e, θ | ≥ 1 token: N, e, θ | verdict |
|---|---|---|---|---|
| dorsals before -i/-e | Romanian | 2,874, 180, 360.9 | 601, 39, 93.9 | productive |
| dorsals before -i | Czech | 3,531, 3, 432.2 | 837, 0, 124.4 | productive |
| dorsals before -i | Italian nouns | 1,481, 739, 202.9 | 520, 264, 83.1 | not productive |

*Caveats:*
- **The definition does the work.** Counting every target-final root not
  seen alternating as an exception, no rule passes in any language
  (claim (a)). The TP verdict therefore assumes that unseen stems do not
  count against the rule. That holds only if they undergo it when they
  turn up before a trigger, which they do (claim (b)); but then they were
  not stored as prespecified.
- **Fragility** (fewest changes that flip the verdict, τ = dictionary /
  ≥ 1): Romanian 181 / 55, Czech 430 / 125; the Italian failure is 537 /
  181 exceptions beyond θ.
- **Small frames.** In the frame of the 1,000 most frequent nouns (Table
  1b), Romanian dorsals before -i have 19 roots and 0 exceptions (θ =
  6.5); Italian 28 roots, 16 / 12 (θ = 8.4).

### (e) Allomorph selection structures the learner's evidence

**Verdict: mixed.** The structural point holds: where a language lets
stems choose a non-palatalizing allomorph, the stems that do not
alternate are shielded and the rule faces few exceptions (Romanian,
Czech). Italian masculines have one plural exponent, and their rule
fails. But the shielding allomorphs are not selected by root
specification (Tables 3 and 4).

- **Romanian.** Of the dorsal roots attested in the plural without
  alternating, 501 take -uri (shielded) and 180 appear unchanged before
  -i/-e (the exceptions): 74% shielded at τ = dictionary, 70% (93 and 39)
  at ≥ 1. Among dorsals without an alternating plural that have a verb
  from the pair, 28 take only -a/-ui, 9 only -i, 2 both. So most
  non-alternating dorsals do select -uri and -a.
- **Czech.** Wherever a palatalizing allomorph is attested, the stem
  changes. At ≥ 1, e = 0 for velar NOM.PL -i (442 roots), DAT/LOC.SG f
  -e (1,124), LOC.PL -ích (484) and adjective NOM.PL -í (456). Not
  alternating shows only as the competing allomorph: 99 velar animates
  attested only with -ové/-é, 1,001 inanimates only with LOC.SG -u, 13
  only with LOC.PL -ách/-ech.
- **Italian.** Masculine dorsals that do not alternate cannot avoid the
  -i plural: 739 of the 1,481 nouns attested before -i keep the velar
  (264 of 520 at ≥ 1). The only shielding is gender: 1,077 feminine
  -ca/-ga nouns whose -e plural is no trigger. This is the mechanism of
  (e) seen from the other side.
- **Selection does not follow specification:**
  - roots whose plural certainly alternates take only the
    non-palatalizing verbalizer in 51 of 128 Romanian dorsals (θ = 26.4)
    and 34 of 57 coronals (θ = 14.1). In Italian, with -are in the set,
    25 of 44 certain roots take only -eggiare/-ire/-are (θ = 11.6); the
    contrast set, roots with -chi/-ghi plurals, takes them 87 of 88
    times;
  - Czech LOC.SG -u is the default of every hard inanimate masculine:
    1,002 of 1,007 velars and 946 of 972 baseline (labial, l, s, z)
    stems take it at ≥ 1. In NOM.PL four fifths of velar animates take
    -i (436 of 533), while three quarters of the baseline take -ové/-é
    (92 of 124);
  - Czech mixed stems: of 111 -log/-urg person nouns, 106 of the 109
    with a NOM.PL take only -ové (or -y), yet 107 of the 109 with a
    LOC.PL take -ích with the change (biologové, biolozích). The same
    stems are shielded in one cell and undergo in another;
  - Romanian -uri dorsals palatalize whenever they take a front-vowel
    derivational suffix (10 of 10 units at τ = dictionary; the draft's
    buluc → buluci is one of the -uri bases choosing -i in
    `romanian_covariation.txt`), although -uri is the allomorph the
    draft makes prespecified roots select.

*Caveat:* "prespecified" is operationalized as "not alternating", the
paragraph's own learner criterion. Czech LOC.SG and VOC counts at ≥ 1
are upper bounds: strictly counted, velar LOC.SG m.in has 4 roots.

### (f) Romanian: dorsals before -i and -e, t, d, s before -i; separate dorsal and coronal rules

**Verdict: supported**, in the harmonized recount and in the legacy
tables.

**Harmonized** (nouns, E0; Tables 1, 3 and 2):

| rule | τ = dictionary: N, e, θ | ≥ 1 token: N, e, θ | |
|---|---|---|---|
| dorsals before -i and -e | 2,874, 180, 360.9 | 601, 39, 93.9 | productive |
| dorsals before -e alone (retain) | 1,258, 1,127, 176.3 | | fails: dorsals change before -e too |
| t, d, s before -i | 2,898, 124, 363.5 | 871, 10, 128.7 | productive |
| t, d, s before -e (retain) | 3,246, 14, 401.5 | | productive; 13 of the 14 are sc → ște, which the account derives through the dorsal rule |
| one coronal rule over -i and -e | 6,144, 3,356, 704.3 | | fails |

- The coronal rule nested as "-e plural → retain, elsewhere →
  palatalize" passes, and 100% of simulated learners reach it.
- Without the outcome-blind flags (`core`): dorsals 1,556 roots, 0
  exceptions; before -i 852, 0; t, d, s before -i 1,975, 6.
- Learners: palatalize productive for 100% of learners in every dorsal
  cell and for t, d, s before -i, at V = 1000 and 5000.
- z is no part of the coronal rule: 116 of 126 z-final roots keep z
  before -i; retain is tolerated and productive for 69% of learners at V
  = 1000.

**Legacy** (Table 9). The draft's top-1000 frame: c 28/0 before -i and
14/0 before -e; g 17/0 and 3/0; t 207/0, d 23/0 and s 11/0 before -i;
t, d, s never change before -e (229, 36, 44 roots). Full lexicon: c
473/0 and 434/0, t 1,522/0.

*Caveats:*
- The legacy and harmonized cells are not the same counts:
  outcome-conditioned exclusions against outcome-blind flags, lemma rows
  against root types, the legacy frame against the R24 frame. Legacy c +
  -i in the top 1,000 has 28 roots; harmonized dorsals before -i among
  the 1,000 most frequent nouns of any final have 19.
- One plural per lemma means no doublets at E0.

*Cross-language:*
- Italian has the dorsal rule only before -i, with no coronal rule: 1,114
  of 1,132 dorsal nouns keep the velar before -e, and none of the 18
  exceptions is attested (338 of 338 keep it at ≥ 1).
- Czech triggers are suffixes, not vowels: a dorsal rule stated on -e
  fails (3,884 and 3,827 exceptions to the two rules among 7,696 roots),
  while r (376 roots, 5 exceptions) and t/d/n (inferred) palatalize
  before -i and s, z, l never do.

So each language needs its own trigger: front vowels for Romanian
dorsals, high front vowels for Romanian coronals and Italian dorsals,
palatalizing allomorphs for Czech.

### (g) If the Italian dorsal pattern fails the threshold, there is no single productive rule, consistent with Krämer

**Verdict: mixed.** The failure is robust; whether the TP then predicts
variation among speakers depends on how the -ico suffix is counted.

- **No single rule:**
  - dorsal nouns before -i: 1,481 roots, 739 and 758 exceptions to the
    two rules, θ = 202.9; at ≥ 1, 520 roots, 264 / 256, θ = 83.1;
  - nouns and adjectives together: 6,612 roots, 1,061 / 5,571; at ≥ 1,
    1,302 roots, 322 / 980;
  - among the 1,000 most frequent nouns: 28 roots, 16 / 12;
  - 96.6% of simulated learners (V = 1000) and 100% (V = 5000) find
    neither rule productive.
- **A conditioned grammar passes.** Post-tonic velar → retain (845
  roots, 81 exceptions, θ = 125.4; at ≥ 1, 299, 22, 52.5) and elsewhere →
  palatalize (3,405, 91, 418.7; 901, 32, 132.4). Counting every -ico
  adjective as a root of its own, 96.3% of learners (V = 1000) and 99.3%
  (V = 5000) end with it. That is a uniform generalization, which would
  not explain variation among speakers.
- **Counting each R26 suffix once**, as R26 says a suffix counts in
  learner evidence, the grammar still passes over root types (818 / 54
  and 1,851 / 91), but the learners divide:

  | | V = 1000 | V = 5000 |
  |---|---|---|
  | τ = dictionary | retain 84.3%, conditioned 15.7% | neither 65.3%, conditioned 34.7% |
  | ≥ 1 token | retain 81.1%, conditioned 18.9% | conditioned 66.9%, neither 33.1% |

  Counting every unstressed -ico/-ica root as one unit (the upper bound),
  99% find "retain" productive. A population split between "retain" and
  "conditioned" resembles Krämer (2009a, CILT 304): one group of speakers
  palatalized every nonce plural, the other none. The 4/1/5 split
  attributed to Krämer (2009b) is unverified.
- **Typed plurals only.** Over plurals an editor typed, not module
  defaults, both stress rules fail: post-tonic → retain 70 roots, 25
  exceptions, θ = 16.5; elsewhere → palatalize 109, 57, θ = 23.2.
- **Verbs.** Adding the first-conjugation verb cells (E1) breaks the
  stress grammar for learners: conditioned 1.6% at V = 1000, 0% at V =
  5000 (Table 2b).

*Caveats:*
- The TP here models type-level evidence, not generalization to nonce
  words.
- The nonce studies are known here only from verified abstracts and
  search excerpts (`docs/italian/README.md`, research notes). In
  Giavazzi's production study the speakers chose the stress of the nonce
  items themselves, so stress and palatalization were not manipulated
  independently.
- Stress is missing for 37% of Italian lexemes. The stress-nested rows
  count only the roots whose stress is known; the learners put the
  others in the elsewhere branch. Whether vowel + -ico adjectives
  (arcaico, laico) are paroxytones moves dozens of roots between the
  post-tonic exceptions and the suffix row.
- The R26 counting choice is a decision for the authors.

## Where the data contradict the draft

1. **The Romanian spell-out TP fails (R20).**
   - The draft makes -i the elsewhere verbalizer of roots that alternate
     (underspecified) and -a/-ui the choice of prespecified roots.
   - Among roots whose plural certainly alternates, 51 of 128 dorsals
     take only -a/-ui (θ = 26.4; 33 of 77 at ≥ 1, θ = 17.7) and 34 of 57
     t, d, s roots (θ = 14.1).
   - It fails under every doublet coding and in the frame of the 1,000
     most frequent nouns (15 roots, 10 exceptions; `romanian_spellout.txt`).
   - The co-variation itself is real (claim (c)). Allomorph choice is a
     strong tendency of root class, not a rule conditioned by the root's
     specification.
2. **Romanian -uri dorsals palatalize before front-vowel suffixes.**
   Prespecified roots select -uri and are inalterable in the draft. Yet
   all 10 front-vowel derivative units of -uri dorsal bases palatalize
   at τ = dictionary: the draft's exception class (c), buluc → buluci, is
   the regular outcome when such a root takes -i at all. It does so in
   10 of 54 units.
3. **Czech mixed stems.** The -log persons choose -ové in NOM.PL but
   palatalize in LOC.PL -ích and in their -ka feminines (bioložka).
   - The Czech `spellout` policy (prespecified iff the lexeme takes a
     non-palatalizing ending anywhere) matches 90.63% of forms, against
     99.92% with every target underspecified.
   - Spell-out leave-one-cell-out holds in one direction per cell, never
     in both (`docs/czech/README.md`).
   - Allomorph choice here is a property of the cell and the lexical
     class, not of a fixed root specification.
4. **Czech LOC.SG -u is not chosen by velars as such.** It is the
   default of every hard inanimate masculine (velar 99.5%, baseline 97%
   of roots at ≥ 1). The palatalizing allomorph is therefore not
   everywhere the elsewhere case, as the draft has -i and -e for
   Romanian.
5. **Italian first-conjugation verbs never palatalize.**
   - -are verbs keep the velar in all 940 velar-stem units (227 at ≥ 1),
     root-stressed -ere verbs always palatalize (249 of 249).
   - An underspecified root (pratici, pubblici, amici) keeps the velar
     in its -are verb (pratichi, pubblichi, amichi): 18 such roots in the
     5,000-noun frame.
   - Adding verb cells drops the stress grammar from 96.3% of learners to
     1.6%.
   - A single root specification needs a category or theme-vowel
     condition (the `category_stress` set passes: Table 3).
6. **Italian -e never triggers.**
   - The draft's dorsal rule fires on any [+syllabic, +front] vowel, so
     it would palatalize amiche and pubbliche.
   - 1,114 of 1,132 Italian dorsal nouns keep the velar before -e; none
     of the 18 exceptions is attested (0 of 338 at ≥ 1).
   - The Italian analysis restricts the trigger to the coronal vocoids /i
     j/ (`docs/italian/phonology.md`): the Romanian trigger does not carry
     over.
7. **Italian doublets** (chirurghi ~ chirurgi, manici ~ manichi) are
   the only failures of the lexical Italian derivation (32 of 35
   mismatched rows). The draft needs two listed underlying forms or a
   variable choice for them.

## Open points

- **Which Krämer 2009?** The CILT chapter (two groups of speakers) or
  the book (4/1/5, unverified); design §10.
- **R26.** Should the -ico suffix count once in learner evidence? The
  Italian learner result turns on it (claim (g)), and Romanian and Czech
  count suffix targets as roots too.
- **A tagged corpus** (UDPipe or MorphoDiTa over the Leipzig sentences)
  would give each syncretic form its cell (Czech LOC.SG/VOC -u, Italian
  -are forms) and each homograph its tokens.
- **Nonce tests** that manipulate stress and suffixhood independently
  would tell an Italian stress default from stored alternation (claims
  (a), (g)).
- **No-evidence derivatives.** Hand-check the R21 counterexamples by
  stratum: the Czech internal-formation rows are the cleanest test, and
  the Romanian strata cannot separate inherited from learned Latin.
- **Coronal derivatives in Romanian.** The infinitive's -i does not
  assibilate. A reading of the finite cells of the derived verbs would
  test the draft's coronal spell-out (klopoT-v → klopotsi).

<!-- BEGIN romanian_compare: generated by lexpipe.langs.romanian_compare; edit the code, not this section -->
## Romanian (harmonized, read-only)

Computed by `lexpipe.langs.romanian_compare` from `data/romanian_lexicon_with_freq.csv`, `data/verb_paradigm_verdict.csv` and the Leipzig table `data/leipzig/freq/ron_news_2024_1M_freq.csv` (all read-only) with the same core code as Italian and Czech (`lexpipe.tp`, `lexpipe.attest`, `lexpipe.learners`). The legacy Romanian outputs are untouched.

Files:

- `romanian_harmonized_tp_long.csv`: `lexpipe.tp.TP_LONG_FIELDS` + `fragility`, `source`, `licence_class`; tables `cells`, `R17a`/`R17b`/`R17c` (the R17 rule family, with τ, E0–E2, three doublet codings, V sweep on lemma and paradigm frequency, strata, nested rules), `nde`, `logtable`, `spellout`.
- `romanian_learners.csv`: simulated learners (R18), 1000 vocabularies per V, seed from the run context.
- `romanian_spellout.txt` (R20), `romanian_noevidence.txt` (R21), `romanian_covariation.txt` (R23: plural class x verbalizer choice), `romanian_logtable.csv` (R27), `romanian_parity.txt`.

Definitions: a root is a noun lemma; N_target = target-final roots of the class (structural NDE classes are separate rows); N_trigger = roots with a realization before the trigger attested at τ; shielded = every attested realization avoids the trigger (-uri, or an -e plural in an -i cell). Subset `all` keeps every outcome-blind flag (suffix targets, palatal context, suppletive stems, -ică suppletion) and keeps the legacy outcome-conditioned exclusions (z→zi 'foreign' plurals, fată/pată) as exceptions; `core` drops the flags. E1 adds the 2sg and 3sg-subjunctive cells of `data/verb_paradigm_verdict.csv`, E2 the lexicon's derived verbs and adjectives, read orthographically (the vowel after the root decides the cell).

Coverage: 51,731 noun roots, 15,734 target-final outside the 2,867 structural NDE roots; E1: 422 nouns with a verdict row, 370 with a readable cell; E2: 1,088 of 2,317 listed derivatives read. One plural per lemma (4 lemmas have two rows, one without a plural), so the three doublet codings differ only at E1/E2.

R17 rows (E0, full lexicon, dominant coding; a non-trigger -e cell shows the `retain` row, whose exceptions are the alternations):

| table | target | trigger | rule | subset | τ | N_target | N_trigger | e | shielded | θ | tolerated | fragility |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| R17a | dorsal | -i/-e | palatalize | all | dictionary | 4063 | 2874 | 180 | 501 | 360.9 | yes | 181 |
| R17a | dorsal | -i/-e | palatalize | all | ge1 | 4063 | 601 | 39 | 93 | 93.9 | yes | 55 |
| R17a | dorsal | -i/-e | palatalize | core | dictionary | 2360 | 1556 | 0 | 442 | 211.7 | yes | 212 |
| R17b | dorsal | -i | palatalize | all | dictionary | 4063 | 1616 | 49 | 1759 | 218.7 | yes | 170 |
| R17b | dorsal | -i | palatalize | all | ge1 | 4063 | 397 | 9 | 297 | 66.3 | yes | 58 |
| R17b | dorsal | -i | palatalize | core | dictionary | 2360 | 852 | 0 | 1146 | 126.3 | yes | 127 |
| R17b | dorsal | -e | retain | all | dictionary | 4063 | 1258 | 1127 | 2117 | 176.3 | no | 951 |
| R17b | dorsal | -e | retain | all | ge1 | 4063 | 204 | 174 | 490 | 38.4 | no | 136 |
| R17b | dorsal | -e | retain | core | dictionary | 2360 | 704 | 704 | 1294 | 107.4 | no | 597 |
| R17c | coronal | -i | palatalize | all | dictionary | 10377 | 2898 | 124 | 4384 | 363.5 | yes | 240 |
| R17c | coronal | -i | palatalize | all | ge1 | 10377 | 871 | 10 | 1321 | 128.7 | yes | 119 |
| R17c | coronal | -i | palatalize | core | dictionary | 9383 | 1975 | 6 | 4355 | 260.3 | yes | 255 |
| R17c | coronal | -e | retain | all | dictionary | 10377 | 3246 | 14 | 4036 | 401.5 | yes | 388 |
| R17c | coronal | -e | retain | all | ge1 | 10377 | 1087 | 4 | 1105 | 155.5 | yes | 152 |
| R17c | coronal | -e | retain | core | dictionary | 9383 | 3231 | 14 | 3099 | 399.8 | yes | 386 |
| R17c | coronal | -i/-e | palatalize | all | dictionary | 10377 | 6144 | 3356 | 1138 | 704.3 | no | 2652 |
| R17c | coronal | -i/-e | palatalize | all | ge1 | 10377 | 1958 | 1093 | 234 | 258.3 | no | 835 |
| R17c | coronal | -i/-e | palatalize | core | dictionary | 9383 | 5206 | 3223 | 1124 | 608.4 | no | 2615 |
| R17c | z | -i | palatalize | all | dictionary | 1294 | 126 | 116 | 941 | 26.1 | no | 90 |
| R17c | z | -i | palatalize | all | ge1 | 1294 | 48 | 46 | 171 | 12.4 | no | 34 |
| R17c | z | -i | palatalize | core | dictionary | 1294 | 126 | 116 | 941 | 26.1 | no | 90 |
| R17c | z | -e | retain | all | dictionary | 1294 | 864 | 0 | 203 | 127.8 | yes | 128 |
| R17c | z | -e | retain | all | ge1 | 1294 | 155 | 0 | 64 | 30.7 | yes | 31 |
| R17c | z | -e | retain | core | dictionary | 1294 | 864 | 0 | 203 | 127.8 | yes | 128 |

Simulated learners at V = 1000 (E0, τ = dictionary):

| table | target | trigger | palatalize | retain | neither | conditioned |
|---|---|---|---|---|---|---|
| R17a | dorsal | -i/-e | 1.0 | 0.0 | 0.0 | 0.0 |
| R17b | dorsal | -i | 1.0 | 0.0 | 0.0 | 0.0 |
| R17b | dorsal | -e | 1.0 | 0.0 | 0.0 | 0.0 |
| R17c | coronal | -i | 1.0 | 0.0 | 0.0 | 0.0 |
| R17c | coronal | -e | 0.0 | 1.0 | 0.0 | 0.0 |
| R17c | coronal | -i/-e | 0.0 | 0.0 | 0.0 | 1.0 |
| R17c | z | -i | 0.0 | 0.7 | 0.3 | 0.0 |
| R17c | z | -e | 0.0 | 1.0 | 0.0 | 0.0 |

Spell-out TP (R20; -i vs -a/-ui on alternating-plural roots, e = roots taking only -a/-ui): dorsal (dictionary, doublets as exceptions: N=128, e=54, θ=26.4; dominant: N=128, e=51, θ=26.4; τ ≥ 1: N=77, e=33, θ=17.7); coronal (dictionary, doublets as exceptions: N=57, e=35, θ=14.1; dominant: N=57, e=34, θ=14.1; τ ≥ 1: N=36, e=21, θ=10.0). The critique's 60 vs θ(164) = 32.2 counts verb pairs (and -ări, -iza, -ăni, -arisi), not roots; see `romanian_spellout.txt`.

-log/-fag/-urg nouns with the learned suffixoid, by class and plural: animate -i 5, person -e 1, person -i 191, thing -e 7, thing -i 1, thing -uri 5; 21 lookalikes (amurg, milog, ...) listed apart.

Parity with the legacy R tables: analysis/romanian_tp_summary_full.csv: 37/37 numeric rows identical; analysis/romanian_tp_summary_downsampled.csv: 37/37 numeric rows identical.

Sources and caveats: the harvested lexicon's noun rows come from ro.wiktionary (about 70%) and en.wiktionary (CC BY-SA 4.0 + GFDL), so aggregate rows name `wiktionary_ro|wiktionary_en` and the -log table names each row's edition. DEX QC (dexonline, terms unverified) filled the plurals Wiktionary lacked for some nouns (lexicon source `| dex_qc`): the aggregates count them, as the legacy tables do, while the -log table leaves those nouns out (they are in the untracked full table). The Leipzig table is lowercased (no case-sensitive counts); stratum `latin_unresolved` mixes inherited and learned Latin. `ro` is a comparison-only row of `lexpipe.spec.LANGUAGES`; run the stages with `python scripts/lang_pipeline.py --lang ro`.
<!-- END romanian_compare -->
