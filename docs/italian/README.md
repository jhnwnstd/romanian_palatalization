# Italian palatalization data

Italian dorsal palatalization (amico/amici, medico/medici vs fuoco/fuochi,
lago/laghi; amica/amiche) counted with the shared multilingual core
(`src/lexpipe/`), so the Italian, Czech and harmonized Romanian tables come
from the same code (design: `docs/multilingual_design.md`, rules R1–R27).
The plugin is `src/lexpipe/langs/italian/`; the data dictionary of every CSV
is `docs/italian/data_dict.md` (generated from `docs/italian/schema.toml`).
The data-layer notes are in `docs/italian/data_layer_notes.md`.

```
python scripts/lang_pipeline.py --lang it            # every stage
python scripts/lang_pipeline.py --lang it tp docs    # some stages
python scripts/lang_pipeline.py --lang it --force tp tp
```

Stages run in order and are skipped while their manifest
(`data/italian/manifests/`) matches (R6):

| stage | module | writes |
|---|---|---|
| extract | `extract` | `data/italian/extract/*` (projections of the kaikki sources, form index) |
| lexicon | `lexicon` | `data/italian/italian_lexicon_raw.csv` (tracked, lean), `italian_lexicon_full.csv` (every column), plural forms, drift test |
| qc | `qc` | `italian_qc.csv`, `analysis/italian/italian_qc_report.txt`, the manual worklist |
| attest | `attestation` | Leipzig counts of every form (`italian_attestation.csv`) |
| features | `features` | `italian_cells.csv`, `italian_lexicon_complete.csv` |
| derivations | `derivations` | `italian_derivations.csv` (derivatives of dorsal bases) |
| verbs | `verbs` | `italian_verbs.csv` (velar-stem verbs and their cells) |
| tp | `tpstage` | `analysis/italian/italian_tp_long.csv`, `italian_learners.csv`, `italian_tp_summary.txt` |
| covariation | `covariation` | `analysis/italian/italian_covariation.txt` |
| spellout | `spellout` | `analysis/italian/italian_spellout.txt` |
| noevidence | `noevidence` | `analysis/italian/italian_noevidence.txt` |
| logtable | `logtable` | `analysis/italian/italian_logtable.csv` (+ gitignored attestation) |
| docs | `docstage` | `docs/italian/data_dict.md`, the generated block below |

On the full data the analysis stages take about 2½ minutes (the
simulated learners, every cell on the same draws in 8 worker processes,
under a minute of it); the data layer about 2½ minutes.

## Sources, provenance and licences

| role | source | licence (R9) |
|---|---|---|
| primary lexicon | kaikki.org raw wiktextract dump of English Wiktionary, top-level `lang_code == "it"` (624,122 entries; 2026-10-01 enwiktionary dump, extraction 2026-10-21) | CC BY-SA 4.0 + GFDL |
| second edition | kaikki.org it.wiktionary extract (Wikizionario) | CC BY-SA 4.0 + GFDL |
| automated QC | Morph-it! 0.4.8 (Zanchetta & Baroni 2005) | CC BY-SA 2.0 / LGPL |
| plural defaults | port of `Module:it-common` (revision 79032535), checked by a drift test | CC BY-SA (module code) |
| frequency, attestation | Leipzig Corpora Collection `ita_news_2024_1M` (Goldhahn et al. 2012), elided clitics split | terms unverified: per-word counts never tracked |
| manual adjudication | Treccani, De Mauro: by hand only (Treccani opts out of text and data mining) | coded value + citation |

Every row records its `source`; tracked CSVs pass the licence gate of
`lexpipe.csvio.write_tracked`. Tracked analysis files carry aggregates and
Wiktionary-derived columns only; per-word Leipzig counts stay in
`data/italian/` (gitignored). Cite Wiktionary's contributors and
wiktextract: Tatu Ylonen, *Wiktextract: Wiktionary as Machine-Readable
Structured Data*, LREC 2022, 1317–1325.

**Provenance classes** (column `provenance` of the schema): `wiktionary_en`,
`wiktionary_it`, `morphit`, `leipzig`, `module`, `curated`, `derived`.
Dictionary-generated forms are not attestation (R11): plurals record
whether an editor typed them (`explicit`) or `Module:it-common` generated
them (`module_default`), and every form carries corpus counts and the four
thresholds τ ∈ {dictionary, ≥ 1 token, ≥ 0.5 per million, ≥ 2 per million}.
The corpus thresholds are tested on folded counts (comparable with the
Romanian tables, design §11), and a form needs at least one token in its
own (lowercase) spelling: a form seen only capitalized is a surname or a
firm name (Bocchi, Monte dei Paschi), which R25 keeps out
(`lexpipe.langs.italian.leipzig.attested_flags`; without this filter 51
velar and 18 palatalizing dorsal roots would count as attested at ≥ 1
token on capitalized tokens alone).

**Dictionary-level rows restate Module:it-common.** 6,890 of the 7,137
dorsal -i plurals listed in the dictionary (96.5%) are Module:it-common
defaults, and every default follows the module's spelling rule: -ico, -aco
and -ologo palatal, every other -co/-go velar. At τ = dictionary the stress,
-ico and -logo rows therefore mostly restate that rule
(Module:it-common-predicted); the rows at τ ≥ 1 rest on attestation (still
95.5% defaults, but used) and are listed first in the tables below. The
`it_stress` rows are also given at τ = `dictionary_explicit` (only plurals
an editor typed: a sample biased toward the module's exceptions). For
each default the attest stage looks up the competing spelling
(`competitor`: amichi for amici); only 10 defaults have one with a
lowercase token that is not a form of another word (sarcofagi,
bolscevichi, loci), 4 of them more often than the default: the corpus
rarely contradicts a default, but it attests only 1,517 of the 7,137
plurals.

**QC tiers.** Tier 1 (stage `qc`) compares every Wiktionary plural with
Morph-it! noun lines and it.wiktionary (`plural_qc`: confirmed, variant,
module_default_only, conflict). Tier 2 is a worklist with Treccani and De
Mauro URLs for checking by hand; decisions go to
`curated/italian/plural_adjudications.csv` (keyed by lexeme_id with guard
columns; a stale guard stops the `features` stage).

## Definitions

- **Lexeme and root (R4).** `lexeme_id` identifies a Wiktionary entry;
  `root_id` groups POS conversions and gender pairs of one etymology
  (amico N/ADJ, amica) and keeps etymologies apart (antico 'ancient' vs
  àntico 'anterior'). TP counts root types; a root's evidence pools all its
  forms. The root key is the stem, the class of the final vowel (-o and -a
  together; -e, -i apart: lattice ≠ lattico, dolce ≠ dolco), the etymon
  with its ending (Latin arca ≠ arcus, orca ≠ Orcus) and the stress; a
  feminine noun that names its masculine counterpart joins its root
  (amica → amico), and an entry without a stress of its own takes the one
  stress of the other entries of its word and page etymology (psicotico N
  joins psicotico ADJ). A dorsal root takes only the plurals of its
  dorsal-target lexemes. Lemma entries that are plurals (pluralia tantum:
  brache, mezzemaniche; "plural of X" entries: alimenti, eliche) are out
  of the domain.
- **Targets (design §5).** Dorsal: stems in -co/-ca/-go/-ga (also -sco,
  -cco, -ggo). Coronal: t/d/s/z (no coronal rule in Italian). Lexically
  palatal (`frontstem`: luce, valigia) and the NDE classes (`occhio`:
  -chio/-chia; `psiche`: singular -che/-ghe) are separate rows.
- **Cells.** `-i` (masculine plurals; masculine -a nouns), `-e`
  (feminine plurals: the non-trigger cell), `-a` (irregular plurals:
  back vowel). For verbs: 2sg, 3sg, 1pl present indicative and the 3sg
  present subjunctive.
- **Stress.** `velar_post_tonic`: the stem-final dorsal sits in the
  syllable right after the stressed one (from the IPA syllable dots and
  ˈ; Giavazzi; Shao et al.). Without a transcription the accented
  hyphenation (en, then it.wiktionary) is used, split at hyphens and
  spaces too, and refused when it spells another word (bi‧o‧ce‧nò‧si on
  biocenotico: 27 lexemes) or accents the i of an -ico formed with the
  suffix (bo‧so‧nì‧co: 8 lexemes); `stress_rejected` records why. Stress
  is missing for 37% of the in-domain lexemes (31,107 of 83,812); those
  roots form the `stress:unknown` rows.
- **Suffix targets (R26).** -ico/-iaco by morphology (an "X + -ico"
  template — `suffix` and `confix` templates write the suffix without its
  hyphen —, or a non-inherited Latin/Greek -icus adjective), never the
  deverbal carico type (a back-formation counts as deverbal only from a
  verb: metallurgo, from metallurgia, is a -urgo suffixoid) and never a
  stressed -ìco (`ico_stressed`: the i of the ending carries the stress;
  arcaico /arˈkaj.ko/, stressed on its a, is a target); the orthographic
  lookalikes are a separate row (`orth:-ico`). R26's "each suffix counts
  once as a morpheme in learner evidence" is the `suffix_once` subset
  (below).
- **Strata (R22).** inherited, borrowed whole, internal formation,
  unknown, from the Wiktionary etymology templates with the classifier
  every language shares (`lexpipe.etymology`): "By surface analysis"
  marks a borrowing analysed as X + -suffix, unless a leading
  affixation names an Italian base (antichità, altoatesino: internal);
  a text that opens with a borrowing ("From Latin …") or says
  "Internally, from" makes an affixation a borrowing.
- **Evidence sets (R15).** E0 nominal inflection; E1 + the cells of
  velar-stem verbs on the same stem (amico ↔ amicare, traffico ↔
  trafficare, bianco ↔ bianchire; of several homograph roots, the ones
  with etymological evidence: pescare goes with pesca 'fishing', not pesca
  'peach'); E2 + derivatives (named links) that are internal formations or
  of unknown stratum (`E2_all` adds the borrowed-whole and inherited ones:
  cecità, sporcizia, porcino).
- **TP rows (R16, R18).** For each target class × trigger cell: N_target
  (target-final roots), N_trigger (roots attested before the trigger at τ),
  shielded roots, e for both rules (palatalize: roots that keep the velar;
  retain: roots that palatalize), θ = N/ln N, the counterfactual over
  N_target, fragility (fewest changes that flip the verdict). Doublets
  (chirurghi ~ chirurgi) are coded three ways (R19): the dominant variant
  (≥ 90% of at least 10 case-sensitive tokens, else the majority), as an
  exception, or excluded; sarcofagi 2 vs sarcofaghi 0 is a true doublet,
  not a dominant variant. A doublet is competing variants of one cell:
  a root whose cells disagree — amici before -i but amiche before -e in
  R17a, the plural medici but the verb form medichi in E1 — is an
  exception to both rules under every coding (`lexpipe.tp`,
  `Observation.cell`), not a doublet decided by tokens.
- **Rule family (R17)** (the cross-language comparison): R17a dorsal /
  _[+front V] over -i and -e; R17b dorsal / _[+high,+front] over -i, with
  -e as its non-trigger cell; R17c coronal / _[+high,+front]. The Italian
  analyses (stress-nested, suffix rows, verbs by conjugation, the pooled
  noun + verb and category-indexed rows) are model fits in separate tables.
- **R26 counts once.** The stress sets and learners are given three ways:
  every -ico adjective a root of its own (`giavazzi`, learner subset
  `all`, as in Romanian and Czech); each R26 suffix counted once
  (`suffix_once`: the -ico targets of a vocabulary are one unit whose
  status pools their plurals, likewise -iaco; a unit is never post-tonic);
  and every root in unstressed -ico/-ica as one unit (`orth_ico_once`, the
  upper bound, since R26 detection is a lower bound).
- **Sampling frame (R24).** The V most frequent noun roots of any final
  segment, V ∈ {500, 1000, 2000, 5000}, by lemma frequency and by paradigm
  sum. Simulated learners draw 1000 vocabularies per V by token frequency
  (`lexpipe.learners`; seed from the run context).
- **Derivatives (stage `derivations`).** Links from etymology templates
  ("X + -suffix"), derived-term lists and the superlatives adjective heads
  list (the `named` tier), plus candidates generated from each dorsal stem
  in both spellings (leghista/legista) and found among the headwords or in
  the corpus (the `candidate` tier, R21). Candidates are rejected when the
  headword's etymology names another base (grecizzare names Grecia), is a
  transmission unrelated to the base (pagina), when a corpus form belongs to
  another lemma (accese), when a frontstem lemma explains the palatal
  spelling (calcio → calcista) or when no token is lowercase (names).
  The frontstem rule can only remove palatal candidates (nothing but the
  base explains a velar ch/gh before the suffix), so the candidate
  tables are also given with those candidates (`…+other_base`). The
  stem edge is read as palatal, velar, back, truncated (-ico dropped:
  arcaismo) or mismatch. Learned set (the co-variation set): -issimo,
  -izzare, -ista, -ismo, -esimo, -ità, -izia, -ia; native control: -ezza,
  -etto, -ino, -ello, -erello, -eria, -essa, -evole, -eggiare, -ire, -ese,
  -iere.
- **Verbs (stage `verbs`).** Velar-stem infinitives (stem-final c/g,
  ch/gh) with their conjugation class: -are; root-stressed -ere (stress
  from the infinitive's IPA, /ˈvin.t͡ʃe.re/ vs /taˈt͡ʃe.re/; the category
  "verbs with root-stressed infinitive" and the accented canonical
  infinitive as fallbacks, with an agreement count in
  `italian_verbs_summary.tsv`); -ére; -ire; -isc-; the dire family. Paradigm types: alternating (vinco/vinci), velar throughout
  (cerchi, bianchisci), over-applying (taccio but tacqui), palatal
  throughout (cucio, cuocio: no velar evidence, out of the TP domain;
  cuocere's past is còssi, so it is palatal throughout, while nuocere,
  with nòcqui, over-applies).
- **Spell-out (R20).** Certain roots (attested -ci/-gi) × an attested
  derivative of a competing pair; e = roots choosing only the
  non-palatalizing member. A derivative whose stem edge is truncated
  (tipizzare from tipico: the -ic- dropped) shows neither member and is
  not counted (28 -izzare derivatives of certain roots at τ =
  dictionary; `italian_spellout.txt` lists them).
- **Co-variation (R23).** Unit: base root × suffix × derivative stratum,
  coded P / N / both; base class from the attested -i plural; per-suffix
  OR and Fisher p, Mantel–Haenszel OR over suffix × stratum, stratified
  permutation; the logistic model every language fits
  (`lexpipe.stats.r23_logistic`: deriv_pal ~ base_class + log lemma
  frequency + base stratum + derivative stratum, base-clustered errors;
  needs statsmodels). Italian bases with -ci/-gi never have a velar
  learned derivative, so the base_class effect is quasi-completely
  separated: the maximum-likelihood estimate does not exist, and the
  report prints where the fit stopped with NA for every standard error
  and p-value; the Fisher, Mantel–Haenszel and permutation results are
  the ones to read. Each table is also given without the -ic- bases
  (abulico → abulicità: plural and derivative are palatal because both
  carry the suffix) and with -ic- as a Mantel–Haenszel stratum.
- **No evidence (R21).** Roots with no -i form attested at τ (feminine-only
  -ca/-ga, no plural, plural below τ) × derivative palatalized; the -ic-
  suffix bases and the masculine-counterpart cases (musica/musico,
  politica/politico) are reported apart. A base whose determiner or
  pronoun entry has an attested masculine plural takes its class from it
  (poco N is invariable, but pochi DET is attested: pochissimo has a
  non-alternating base).

## Results

The block below is regenerated by the `docs` stage from the analysis
outputs. In words:

- **Read τ ≥ 1 first.** 96.5% of the dorsal -i plurals at τ = dictionary
  are Module:it-common defaults (all following its -ico/-aco/-ologo rule),
  so the dictionary rows mostly restate the module's spelling rule; the
  τ ≥ 1 rows count plurals the corpus attests. Over typed plurals only
  (τ = `dictionary_explicit`, biased toward the module's exceptions) both
  stress rules fail: post-tonic → retain N 70, e 25 > θ 16.5; elsewhere
  → palatalize N 109, e 57 > θ 23.2.
- **No single rule over the -i cell is productive.** Pooled over every
  dorsal root (R17b), both palatalize and retain fail at every threshold
  (τ ≥ 1: N 1,302, e 322 and 980, θ 181.5), in the nouns-only table too;
  the R17a rule over -i and -e fails worse: the feminine -e cell never
  palatalizes (only the editing errors flagged by QC do), and every root
  with -ci/-gi before -i but -che/-ghe before -e — every -ico adjective:
  pubblici, pubbliche — is an exception to both rules. Coronals never
  palatalize (N = e by construction).
- **The stress grammar holds over root types; learners reach it only if
  every -ico adjective counts as a root.** Nested in Yang's order, post-
  tonic → retain and elsewhere → palatalize are both tolerated (τ ≥ 1:
  299 roots, e 22; 901, e 32), and still with each R26 suffix counted
  once (292, e 15; 449, e 32). With every root in unstressed -ico/-ica as
  one unit (the upper bound), elsewhere → palatalize fails at τ =
  dictionary (275, e 64 > θ 49.0) and passes at τ ≥ 1 with fragility 2
  (104, e 21, θ 22.4). The simulated learners depend on the same
  choice: counting every -ico root, 96% (V = 1000) and 99% (V = 5000) end
  up with the conditioned grammar (E0, τ = dictionary; 98% and 100% at
  τ ≥ 1); counting each R26 suffix once, as R26 asks, only 16% and 35% do
  (19% and 67% at τ ≥ 1) — at V = 1000, 84% find "retain before -i"
  productive outright (81% at τ ≥ 1), and so do 84% in the plain R17b -i
  cell with the suffix counted once; with every -ico/-ica root as one
  unit, about 99% find retain productive. The -ico suffix targets
  themselves are categorical (τ ≥ 1: 473 roots, no velar). g-stems fail
  in the elsewhere branch (Flor's voicing effect: profughi, naufraghi),
  k-stems pass.
- **Verbs split by conjugation.** -are verbs never palatalize, root-stressed
  -ere verbs always do, -isc- verbs with a velar keep it; the pooled noun +
  verb rule fails and the category-indexed rules pass only when the noun
  stress condition is added. Adding the verb cells to the nouns' evidence
  (E1) breaks the stress grammar for the simulated learners (conditioned:
  2% at V = 1000, none at V = 5000, τ = dictionary): 55 dorsal roots of the
  5000-noun frame change status, 35 feminine -ca/-ga nouns shielded in E0
  that meet -i in their -are verb (classifichi, verifichi), which keeps
  the velar (18 not post-tonic: exceptions to "elsewhere palatalize"), and
  18 roots whose palatal plural meets a velar verb form (pratici but
  pratichi, pubblici but pubblichi, amici but amichi), exceptions to both
  rules.
- **Derivation tracks the plural for the learned suffixes, but mostly
  through the -ic- suffix.** Among Italian-internal learned derivatives
  every base with a dictionary-listed -ci/-gi plural has a palatalized
  derivative (258 units at τ = dictionary, 78 at τ ≥ 1); bases with
  -chi/-ghi mostly keep the velar (exceptions: dialogizzare, magismo,
  opacizzare, panicità, valgismo). But 251 of the 258 palatal-base units
  have an -ic- base (abulico → abulicità: both carry the suffix);
  without them the table is 7/0 vs 5/22 (Fisher p 1.5e-4) at τ =
  dictionary and 1/0 vs 1/5 (p 0.29) at τ ≥ 1, where greco → grecista is
  the only palatal base left. The native suffixes keep the velar but for
  two -ino derivatives of velar bases (fungino, mancino). The spell-out
  TP passes for the learned suffixes and for -izzare vs -eggiare/-ire —
  with fragility 1 at τ = dictionary (26 roots, 7 exceptions, θ 8.0): one
  more exception flips it — but fails once -are joins the competing set
  (-are is Italian's default verbalizer, and most such verbs on these
  stems are Latin borrowings).
- **No-evidence roots do not resist.** The few learned derivatives of roots
  without an attested alternating plural mostly palatalize, and nearly all
  of those have an -ic- suffix base (armonicista, rubricista, logicista,
  and at τ ≥ 1 the -ità nouns of adjectives whose plural is unattested)
  or a masculine counterpart (musica/musico → musicista), as the critique
  anticipated; the velar cases are -ca/-ga nouns without the suffix
  (cronachista, droghista, bisarchista, ribechista; leghista among the
  generated candidates). pochissimo is not one: poco's plural pochi is a
  determiner/pronoun form, attested and velar, so its base is
  non-alternating.

<!-- BEGIN generated results: written by the docs stage (lexpipe.langs.italian.docstage); edit the code, not this block -->

### Headline numbers (generated)

Unit: root types (R4); E0 unless stated; doublets coded by the dominant variant (R19); θ = N/ln N. Full rows: `analysis/italian/italian_tp_long.csv`; digest: `analysis/italian/italian_tp_summary.txt`.

**Cross-language rule family (R17).**

| row | τ | N_target | N | e (palatalize) | e (retain) | θ | productive |
|---|---|---|---|---|---|---|---|
| R17a dorsal/_[+front V], -i/-e | ge1 | 7732 | 1842 | 1494 | 980 | 245.0 | neither |
| R17a dorsal/_[+front V], -i/-e | dictionary | 7732 | 7537 | 7169 | 5584 | 844.2 | neither |
| R17b dorsal/_[+high,+front], -i | ge1 | 7732 | 1302 | 322 | 980 | 181.5 | neither |
| R17b dorsal/_[+high,+front], -i | dictionary | 7732 | 6612 | 1061 | 5571 | 751.7 | neither |
| R17b -i, nouns only | ge1 | 2775 | 520 | 264 | 256 | 83.1 | neither |
| R17b -i, nouns only | dictionary | 2775 | 1481 | 739 | 758 | 202.9 | neither |
| R17b -i, core (no -ico, no -logo/-fago/-urgo) | ge1 | 5489 | 747 | 308 | 439 | 112.9 | neither |
| R17b -i, core (no -ico, no -logo/-fago/-urgo) | dictionary | 5489 | 4447 | 1025 | 3436 | 529.4 | neither |
| R17b non-trigger -e (feminine plurals) | ge1 | 7732 | 1248 | 1248 | 0 | 175.1 | retain |
| R17b non-trigger -e (feminine plurals) | dictionary | 7732 | 6566 | 6553 | 17 | 747.0 | retain |
| R17c coronal/_[+high,+front], -i | ge1 | 21161 | 5601 | 5601 | 0 | 649.0 | retain |
| R17c coronal/_[+high,+front], -i | dictionary | 21161 | 15935 | 15935 | 0 | 1646.8 | retain |

**Italian model fits, suffix, verb and spell-out rows.**

| row | τ | N_target | N | e (palatalize) | e (retain) | θ | productive |
|---|---|---|---|---|---|---|---|
| velar post-tonic (paroxytone type) | ge1 | 1269 | 299 | 277 | 22 | 52.5 | retain |
| velar post-tonic (paroxytone type) | dictionary | 1269 | 845 | 770 | 81 | 125.4 | retain |
| velar not post-tonic (proparoxytone type) | ge1 | 3762 | 901 | 32 | 869 | 132.4 | palatalize |
| velar not post-tonic (proparoxytone type) | dictionary | 3762 | 3405 | 91 | 3325 | 418.7 | palatalize |
| velar post-tonic, typed plurals only | dictionary_explicit | 1269 | 70 | 46 | 25 | 16.5 | neither |
| velar not post-tonic, typed plurals only | dictionary_explicit | 3762 | 109 | 57 | 58 | 23.2 | neither |
| -ico suffix targets (R26) | ge1 | 1926 | 473 | 0 | 473 | 76.8 | palatalize |
| -ico suffix targets (R26) | dictionary | 1926 | 1848 | 1 | 1847 | 245.7 | palatalize |
| -logo person nouns | ge1 | 232 | 58 | 0 | 58 | 14.3 | palatalize |
| -logo person nouns | dictionary | 232 | 232 | 2 | 232 | 42.6 | palatalize |
| -logo thing nouns | ge1 | 10 | 5 | 5 | 0 | 3.1 | retain |
| -logo thing nouns | dictionary | 10 | 10 | 9 | 1 | 4.3 | retain |
| deverbal from -are (carico type) | ge1 | 76 | 37 | 34 | 3 | 10.2 | retain |
| deverbal from -are (carico type) | dictionary | 76 | 61 | 57 | 4 | 14.8 | retain |
| verbs -are, 2sg/1pl/sbjv | ge1 | 940 | 227 | 227 | 0 | 41.8 | retain |
| verbs -are, 2sg/1pl/sbjv | dictionary | 940 | 940 | 940 | 0 | 137.3 | retain |
| verbs root-stressed -ere | ge1 | 252 | 57 | 0 | 57 | 14.1 | palatalize |
| verbs root-stressed -ere | dictionary | 252 | 249 | 0 | 249 | 45.1 | palatalize |
| verbs -isc- | ge1 | 36 | 3 | 3 | 0 | 2.7 | retain |
| verbs -isc- | dictionary | 36 | 36 | 36 | 0 | 10.0 | retain |
| pooled noun + verb | ge1 | 8992 | 1597 | 552 | 1045 | 216.5 | neither |
| pooled noun + verb | dictionary | 8992 | 7869 | 2037 | 5852 | 877.2 | neither |
| spell-out: -izzare vs -eggiare/-ire/-are | ge1 | 980 | 17 | 11 | 6 | 6.0 | retain |
| spell-out: -izzare vs -eggiare/-ire/-are | dictionary | 5551 | 44 | 25 | 20 | 11.6 | neither |
| spell-out: same without -are | ge1 | 980 | 8 | 2 | 6 | 3.8 | palatalize |
| spell-out: same without -are | dictionary | 5551 | 26 | 7 | 20 | 8.0 | palatalize |
| spell-out: palatal vs velar -ista/-ismo/-ità | ge1 | 980 | 73 | 0 | 73 | 17.0 | palatalize |
| spell-out: palatal vs velar -ista/-ismo/-ità | dictionary | 5551 | 231 | 0 | 231 | 42.4 | palatalize |

**Nested rule sets (Yang's order).**

| nested set | τ | rule (subset) | N | e | θ | tolerated |
|---|---|---|---|---|---|---|
| it_stress dorsal giavazzi | ge1 | retain (nested:post_tonic) | 299 | 22 | 52.5 | True |
| it_stress dorsal giavazzi | ge1 | palatalize (nested:elsewhere) | 901 | 32 | 132.4 | True |
| it_stress dorsal giavazzi | dictionary | retain (nested:post_tonic) | 845 | 81 | 125.4 | True |
| it_stress dorsal giavazzi | dictionary | palatalize (nested:elsewhere) | 3405 | 91 | 418.7 | True |
| it_stress dorsal giavazzi | dictionary_explicit | retain (nested:post_tonic) | 70 | 25 | 16.5 | False |
| it_stress dorsal giavazzi | dictionary_explicit | palatalize (nested:elsewhere) | 109 | 57 | 23.2 | False |
| it_stress dorsal suffix_once | ge1 | retain (nested:post_tonic) | 292 | 15 | 51.4 | True |
| it_stress dorsal suffix_once | ge1 | palatalize (nested:elsewhere) | 449 | 32 | 73.5 | True |
| it_stress dorsal suffix_once | dictionary | retain (nested:post_tonic) | 818 | 54 | 122.0 | True |
| it_stress dorsal suffix_once | dictionary | palatalize (nested:elsewhere) | 1851 | 91 | 246.0 | True |
| it_stress dorsal suffix_once | dictionary_explicit | retain (nested:post_tonic) | 70 | 25 | 16.5 | False |
| it_stress dorsal suffix_once | dictionary_explicit | palatalize (nested:elsewhere) | 109 | 57 | 23.2 | False |
| it_stress dorsal orth_ico_once | ge1 | retain (nested:post_tonic) | 286 | 9 | 50.6 | True |
| it_stress dorsal orth_ico_once | ge1 | palatalize (nested:elsewhere) | 104 | 21 | 22.4 | True |
| it_stress dorsal orth_ico_once | dictionary | retain (nested:post_tonic) | 799 | 35 | 119.6 | True |
| it_stress dorsal orth_ico_once | dictionary | palatalize (nested:elsewhere) | 275 | 64 | 49.0 | False |
| it_stress dorsal orth_ico_once | dictionary_explicit | retain (nested:post_tonic) | 70 | 25 | 16.5 | False |
| it_stress dorsal orth_ico_once | dictionary_explicit | palatalize (nested:elsewhere) | 82 | 30 | 18.6 | False |
| it_stress g giavazzi | ge1 | retain (nested:post_tonic) | 60 | 1 | 14.7 | True |
| it_stress g giavazzi | ge1 | palatalize (nested:elsewhere) | 76 | 21 | 17.5 | False |
| it_stress g giavazzi | dictionary | retain (nested:post_tonic) | 168 | 5 | 32.8 | True |
| it_stress g giavazzi | dictionary | palatalize (nested:elsewhere) | 206 | 48 | 38.7 | False |
| it_nounverb dorsal category_stress | ge1 | retain (nested:verb_are) | 227 | 0 | 41.8 | True |
| it_nounverb dorsal category_stress | ge1 | retain (nested:verb_isc) | 3 | 0 | 2.7 | True |
| it_nounverb dorsal category_stress | ge1 | retain (nested:noun_post_tonic) | 299 | 22 | 52.5 | True |
| it_nounverb dorsal category_stress | ge1 | palatalize (nested:elsewhere) | 1068 | 45 | 153.2 | True |
| it_nounverb dorsal category_stress | dictionary | retain (nested:verb_are) | 940 | 0 | 137.3 | True |
| it_nounverb dorsal category_stress | dictionary | retain (nested:verb_isc) | 36 | 0 | 10.0 | True |
| it_nounverb dorsal category_stress | dictionary | retain (nested:noun_post_tonic) | 845 | 81 | 125.4 | True |
| it_nounverb dorsal category_stress | dictionary | palatalize (nested:elsewhere) | 6048 | 291 | 694.6 | True |

**Simulated learners (1000 vocabularies per V, R18).**

| cell | subset | E | τ | V | palatalize | retain | neither | conditioned |
|---|---|---|---|---|---|---|---|---|
| R17a -i/-e | with_nde | E0 | dictionary | 1000 | 0.000 | 0.513 | 0.487 | 0.000 |
| R17a -i/-e | with_nde | E0 | dictionary | 5000 | 0.000 | 0.000 | 1.000 | 0.000 |
| R17a -i/-e | with_nde | E1 | dictionary | 1000 | 0.000 | 0.494 | 0.506 | 0.000 |
| R17a -i/-e | with_nde | E1 | dictionary | 5000 | 0.000 | 0.000 | 1.000 | 0.000 |
| R17b -i | all | E0 | dictionary | 1000 | 0.000 | 0.034 | 0.966 | 0.000 |
| R17b -i | all | E0 | dictionary | 5000 | 0.000 | 0.000 | 1.000 | 0.000 |
| R17b -i | all | E0 | ge1 | 1000 | 0.000 | 0.023 | 0.977 | 0.000 |
| R17b -i | all | E0 | ge1 | 5000 | 0.000 | 0.000 | 1.000 | 0.000 |
| R17b -i | all | E1 | dictionary | 1000 | 0.000 | 0.136 | 0.864 | 0.000 |
| R17b -i | all | E1 | dictionary | 5000 | 0.000 | 0.000 | 1.000 | 0.000 |
| R17b -i | all | E1 | ge1 | 1000 | 0.000 | 0.114 | 0.886 | 0.000 |
| R17b -i | all | E1 | ge1 | 5000 | 0.000 | 0.000 | 1.000 | 0.000 |
| R17b -e | all | E0 | dictionary | 1000 | 0.000 | 1.000 | 0.000 | 0.000 |
| R17b -e | all | E0 | dictionary | 5000 | 0.000 | 1.000 | 0.000 | 0.000 |
| R17b -e | all | E1 | dictionary | 1000 | 0.000 | 1.000 | 0.000 | 0.000 |
| R17b -e | all | E1 | dictionary | 5000 | 0.000 | 1.000 | 0.000 | 0.000 |
| R17c -i | all | E0 | dictionary | 1000 | 0.000 | 1.000 | 0.000 | 0.000 |
| R17c -i | all | E0 | dictionary | 5000 | 0.000 | 1.000 | 0.000 | 0.000 |
| it_stress -i | all | E0 | dictionary | 1000 | 0.000 | 0.034 | 0.003 | 0.963 |
| it_stress -i | all | E0 | dictionary | 5000 | 0.000 | 0.000 | 0.007 | 0.993 |
| it_stress -i | all | E0 | ge1 | 1000 | 0.000 | 0.023 | 0.002 | 0.975 |
| it_stress -i | all | E0 | ge1 | 5000 | 0.000 | 0.000 | 0.000 | 1.000 |
| it_stress -i | all | E1 | dictionary | 1000 | 0.000 | 0.136 | 0.848 | 0.016 |
| it_stress -i | all | E1 | dictionary | 5000 | 0.000 | 0.000 | 1.000 | 0.000 |
| it_stress -i | all | E1 | ge1 | 1000 | 0.000 | 0.114 | 0.546 | 0.340 |
| it_stress -i | all | E1 | ge1 | 5000 | 0.000 | 0.000 | 1.000 | 0.000 |
| R17b -i | suffix_once | E0 | dictionary | 1000 | 0.000 | 0.843 | 0.157 | 0.000 |
| R17b -i | suffix_once | E0 | dictionary | 5000 | 0.000 | 0.000 | 1.000 | 0.000 |
| R17b -i | suffix_once | E0 | ge1 | 1000 | 0.000 | 0.811 | 0.189 | 0.000 |
| R17b -i | suffix_once | E0 | ge1 | 5000 | 0.000 | 0.000 | 1.000 | 0.000 |
| R17b -i | suffix_once | E1 | dictionary | 1000 | 0.000 | 0.955 | 0.045 | 0.000 |
| R17b -i | suffix_once | E1 | dictionary | 5000 | 0.000 | 0.042 | 0.958 | 0.000 |
| R17b -i | suffix_once | E1 | ge1 | 1000 | 0.000 | 0.949 | 0.051 | 0.000 |
| R17b -i | suffix_once | E1 | ge1 | 5000 | 0.000 | 0.032 | 0.968 | 0.000 |
| it_stress -i | suffix_once | E0 | dictionary | 1000 | 0.000 | 0.843 | 0.000 | 0.157 |
| it_stress -i | suffix_once | E0 | dictionary | 5000 | 0.000 | 0.000 | 0.653 | 0.347 |
| it_stress -i | suffix_once | E0 | ge1 | 1000 | 0.000 | 0.811 | 0.000 | 0.189 |
| it_stress -i | suffix_once | E0 | ge1 | 5000 | 0.000 | 0.000 | 0.331 | 0.669 |
| it_stress -i | suffix_once | E1 | dictionary | 1000 | 0.000 | 0.955 | 0.044 | 0.001 |
| it_stress -i | suffix_once | E1 | dictionary | 5000 | 0.000 | 0.042 | 0.958 | 0.000 |
| it_stress -i | suffix_once | E1 | ge1 | 1000 | 0.000 | 0.949 | 0.034 | 0.017 |
| it_stress -i | suffix_once | E1 | ge1 | 5000 | 0.000 | 0.032 | 0.968 | 0.000 |
| it_stress -i | orth_ico_once | E0 | dictionary | 1000 | 0.000 | 0.990 | 0.000 | 0.010 |
| it_stress -i | orth_ico_once | E0 | dictionary | 5000 | 0.000 | 0.999 | 0.000 | 0.001 |
| it_stress -i | orth_ico_once | E0 | ge1 | 1000 | 0.000 | 0.988 | 0.000 | 0.012 |
| it_stress -i | orth_ico_once | E0 | ge1 | 5000 | 0.000 | 0.998 | 0.000 | 0.002 |
| it_stress -i | orth_ico_once | E1 | dictionary | 1000 | 0.000 | 0.997 | 0.002 | 0.001 |
| it_stress -i | orth_ico_once | E1 | dictionary | 5000 | 0.000 | 1.000 | 0.000 | 0.000 |
| it_stress -i | orth_ico_once | E1 | ge1 | 1000 | 0.000 | 0.998 | 0.000 | 0.002 |
| it_stress -i | orth_ico_once | E1 | ge1 | 5000 | 0.000 | 1.000 | 0.000 | 0.000 |

**Co-variation (R23): base -i plural × learned derivative, Italian-internal formations, named links.**

| τ | suffix set | P base: P / N | N base: P / N | both | Fisher p | MH OR (strata) |
|---|---|---|---|---|---|---|
| dictionary | learned | 258 / 0 | 5 / 22 | 0 | 2.05e-28 | inf (6) |
| dictionary | learned, no -ic- bases | 7 / 0 | 5 / 22 | 0 | 0.000147 | inf (6) |
| dictionary | native | 0 / 27 | 2 / 102 | 0 | 1 | 0 (11) |
| ge1 | learned | 78 / 0 | 1 / 5 | 0 | 1.94e-07 | inf (6) |
| ge1 | learned, no -ic- bases | 1 / 0 | 1 / 5 | 0 | 0.286 | inf (6) |
| ge1 | native | 0 / 6 | 2 / 52 | 0 | 1 | undefined |

**No-evidence roots (R21): learned-set units by base class.**

| τ | links | alternating P / N | non-alternating P / N | no evidence P / N |
|---|---|---|---|---|
| dictionary | named | 280 / 0 | 9 / 40 | 10 / 4 |
| dictionary | named+candidates | 322 / 1 | 19 / 47 | 12 / 4 |
| dictionary | named+candidates+other_base | 331 / 1 | 23 / 47 | 12 / 4 |
| ge1 | named | 89 / 0 | 3 / 14 | 8 / 1 |
| ge1 | named+candidates | 120 / 1 | 4 / 25 | 17 / 3 |
| ge1 | named+candidates+other_base | 123 / 1 | 4 / 25 | 19 / 3 |

<!-- END generated results -->

## Research notes

Sources verified in the research phase (`data/research_v2/it_phon.json`,
whose `verification` section overrides its `research` section):

- **Krämer (2009a)**, "Crypto-Variation in Italian Velar Palatalisation",
  in Masullo, O'Rourke & Huang (eds.), *Romance Linguistics 2007* (CILT
  304), 193–208, doi:10.1075/cilt.304.13kra. Abstract (verified via
  Crossref): speakers given nonce singulars split into a group that
  palatalized every new plural and a group that palatalized none;
  support for Pater's (2006) lexically indexed constraints. Steddy cites
  this chapter as "Krämer 2007".
- **Krämer (2009b)**, *The Phonology of Italian*, OUP, ch. 4 "Segmental
  Phonology", pp. 44–126, doi:10.1093/oso/9780199290796.001.0001. Kyle
  Gorman's *Wellformedness* post "Italian palatalization" summarizes §4.2.1
  as 10 speakers (4 never palatalized, 1 always, 5 mixed), much higher type
  frequency for non-palatalization, and a "personal decision" quote
  attributed to Krämer (2012:125). These figures were seen only in search
  excerpts and are **unverified**.
- **Krämer (2016)**, "Variation and change in Italian phonology",
  Linguistik Aktuell 234, 205–236, doi:10.1075/la.234.08kra.
- **Giavazzi.** The MIT dissertation (2010, hdl 1721.1/62408) and a
  separate manuscript "Stress-conditioned palatalization in Italian" (Flor
  dates it 2010, Steddy 2012), the source of the nonce results Flor cites;
  her review of Krämer's book, *Phonology* 27(2), 332–341,
  doi:10.1017/S0952675710000151. The WCCFL 27 poster's production study (7
  speakers palatalizing only antepenultimate-stressed forms) is known from
  search excerpts only; its speakers chose the stress themselves, so stress
  and palatalization were not manipulated independently.
- **Flor**, "A unified account of Velar Palatalization in Italian", NELS 50:
  the stress, height (-e never triggers) and voicing effects; NDEB pairs;
  the superlative "follows the base" — stated over base stress, which
  amicissimo contradicts; "no /k/ optionality", which the /k/ doublets
  (manici ~ manichi, stomaci ~ stomachi) contradict.
- **Steddy**, "Palatalisation across the Italian lexicon", ROA 1229:
  -izzare tracks the masculine plural, derivational -ire never palatalizes,
  fuggire is the one exception to his verb account; his 45 and 63 counts
  come from corpus searches over three resources, not specifically
  Morph-it!.
- **Shao, Buech, Hermes & Giavazzi (2023)**, Interspeech 1833–1837,
  doi:10.21437/Interspeech.2023-1456 (comico/comici vs caduco/caduchi);
  their ICPhS 2023 paper (cardiaci vs ubriachi) is at
  https://www.internationalphoneticassociation.org/icphs-proceedings/ICPhS2023/full_papers/814.pdf;
  the 2025 follow-up is *J. Phonetics* 108, 101371,
  doi:10.1016/j.wocn.2024.101371.
- **Celata & Bertinetto (2005)**, "Lexical access in Italian: words with
  and without palatalization", *Lingue e linguaggio* 4(2), 293–318:
  palatalization largely unpredictable, whole-word access.
- **Treccani** (*La grammatica italiana*, 2012): "-co, -go" (proparoxytones
  -ci/-gi, paroxytones -chi/-ghi, "sia pure non in tutti i casi"); "-fago,
  -logo" (persons -logi/-fagi, things and adjectives -loghi/-faghi;
  -logh- person plurals "da evitarsi"); "Sarcofaghi o sarcofagi?" (both
  correct, sarcofaghi far more common). **Accademia della Crusca**:
  consulenza 39 (2002; chirurghi ~ chirurgi, farmaci ~ farmachi, manici ~
  manichi, asparagi; esofagi, sarcofagi), consulenza 272 (Setti 2009/2010,
  quoting Serianni's *Italiano*, p. 100: only abachi, amici and greci are
  stress-rule exceptions among her examples), consulenza 202 (Migliorini's
  -cia/-gia spelling rule).
- **Other**: Steriade (2008), "A pseudo-cyclic effect in Romanian
  morphophonology", in Bachrach & Nevins (eds.), *Inflectional Identity*,
  313–360; Yang (2005), *LIVY* 5, 265–302, and Yang (2016), *The Price of
  Linguistic Productivity*; Inkelas & Cho (1993), *Language* 69(3),
  529–574; Wilson (2006), *Cognitive Science* 30(5), 945–982 (the i vs e
  trigger asymmetry). A 2025 Charles University bachelor's thesis (48
  speakers; stress the strongest predictor) is known from a search snippet
  only.

**Fact-check corrections applied** (it_phon.json verification):
grecizzare is formed from Grecia, politichese from politica and musicista
from musica — none is evidence that a derivative follows a -ci base
(`italian_derivations.csv`: grecizzare is rejected as a derivative of
greco; politichese and musicista keep their feminine bases, with the
masculine counterpart recorded); porcino, porcello and porcile are
inherited (stratum filter); the heavy-penult generalization is refuted
(the -urgo persons and belga palatalize); nuocere over-applies (nòccio
but nòcqui) and dire is a suppletive -ire verb (the `verbs` classes),
while cuocere, which the design lists with nuocere, shows no velar
anywhere in its paradigm (cuòcio, còssi) and is palatal throughout, out
of the TP domain; carico's etymology reads "short past participle of
caricare" (the deverbal rule matches it).

## Data findings for the authors

- The data layer's stratum rule read "From antico + -ità. By surface
  analysis, borrowed from Latin antīquitātem." as a surface-analysed
  borrowing; antichità became an internal formation. The shared
  classifier (2026-10-23) generalizes that fix: 61 in-domain
  affixations whose leading "From X + -suffix" names Italian bases
  (altoatesino, antropometrico) are internal formations now, and 856
  in-domain lexemes without templates whose text reads "From <Italian
  word>." are internal formations instead of unknown
  (`data/research_v2/migration_diff.md`). The Italian compound templates
  `it-verb-obj` and `it-verb-verb` (forasacco: fora + sacco) are read as
  compounds: about 1,150 in-domain lexemes are internal formations
  instead of unknown.
- Most post-tonic palatalizers are the amico type (amico, nemico, greco,
  porco, belga) plus vowel + -ico adjectives (arcaico, ebraico, eroico,
  stoico, laico) whose IPA writes the vowel sequence as a diphthong, which
  makes them paroxytones: 47 of the 81 post-tonic palatalizers at the
  dictionary threshold (row `orth:-Vico post-tonic` in `it_suffix`), 26
  of them R26 suffix targets now that the R26 check reads the stressed
  vowel (arcaico /arˈkaj.ko/ is stressed on its a: the -ico is
  unstressed). The laico syllabification question decides more than half
  of the post-tonic exceptions; with each R26 suffix counted once the
  post-tonic branch keeps 54 (τ = dictionary).
- Hyphenation stress can be wrong in ways the checks do not catch: an
  accent on the i of an -ico adjective without a template or an Italian
  base (geotermoelettrico 'trì', geomatica 'tì') still makes it a
  post-tonic palatalizer.
- Candidate derivatives are noisy even after the filters (legista, from
  Latin legista, is accepted as a palatal derivative of lega because its
  entry has no etymology); the main tables use the `named` tier, and the
  candidate tier is a sensitivity analysis.
- -are verbs on these stems collide with nouns in the corpus (the 3sg
  càrica, mèdica; the 2sg càrichi = the plural carichi): verb attestation is
  inflated, which is why the verb tables also report
  `all_unambiguous`.

## Open questions

- Which "Krämer 2009" is meant (the CILT chapter's two speaker groups or
  the book's 4/1/5 split, the latter unverified)?
- Should the first conjugation count against a single dorsal rule (the
  pooled noun + verb row fails) or is the rule category-indexed (it passes
  with the stress condition)? Both are reported; E1 evidence shows the
  conflict for a learner who links nouns and -are verbs.
- Are vowel + -ico adjectives paroxytones (as Wiktionary's IPA has them) or
  proparoxytones with an -ico suffix? The answer moves dozens of roots
  between the post-tonic exception list and the suffix row.
- Is -are part of the spell-out competition? The TP passes without it and
  fails with it.
- Stress is missing for 37% of lexemes; the `stress:unknown` rows are
  mostly -ico adjectives and palatalize. The stress-nested TP rows are
  computed over the stress-known roots; the simulated learners and the
  `it_nounverb` category_stress set put stress-unknown roots in the
  elsewhere branch (leaving them out moves the learners' conditioned
  share from 0.963 to 0.963 at V = 1000 and from 0.993 to 0.999 at
  V = 5000; E0, τ = dictionary, every root counted).
- R26 asks each suffix to count once in learner evidence. The Italian
  learners report both ways (`all`, `suffix_once`); the Romanian and
  Czech learner cells count every suffix target as its own root, so the
  cross-language R17 learner rows compare like with like, but none of
  them applies R26's rule.
- The R23 logistic model is fitted (statsmodels, optional extra) but is
  quasi-completely separated in every Italian section (the report prints
  NA for its standard errors and p-values); a penalized (Firth) or
  Bayesian fit, or the mixed model with base and suffix intercepts (R
  lme4), was not attempted.
- Leipzig's terms have not been checked by hand: per-word counts stay out
  of tracked files; the tracked tables carry aggregates over τ only.
