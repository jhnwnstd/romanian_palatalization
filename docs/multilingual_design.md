# Multilingual extension: Italian and Czech

Status: design v2, adopted 2026-10-23 after a scientific and an engineering critique. Implemented for Italian, Czech and the read-only Romanian comparison, audited and fixed (2026-10-23). The cross-language synthesis is `analysis/comparison/README.md`, with its tables in `analysis/comparison/crosslang_summary.csv`, built by `lexpipe.compare`. Departures from the rules below are listed in §11.

Research basis (local, gitignored): `data/research_v2/*.json`. It contains the code map, the Italian and Czech phonology and resource reports, their adversarial fact-checks, the stage-1 build reports and both critiques.

Decision rules are numbered (R1, R2, …) so code reviews can cite them.

## 1. Decision: one repository, language plugins, Romanian untouched

We do **not** fork. Italian and Czech are language plugins on a shared core. The legacy Romanian pipeline keeps producing byte-identical outputs.

Why:

- **The engine is already shared.** The Logical Phonology engine (`src/phonology/`) is language-neutral apart from hooks that now exist (`phonology.profile.AnalysisProfile`: `ur_builder`, `cluster_key`, `normalisers`, `categoriser`, `row_factory`).
- **Comparison needs one codebase.** The payoff is a cross-linguistic comparison (Tolerance Principle thresholds, allomorph selection, inflection dependence, learnability), and that requires one set of counting code. A harmonized Romanian table is therefore computed by the same code (§7). The legacy R and Python Romanian outputs stay untouched.
- **Romanian has a pinned baseline.** It replays byte-identically, and `tests/test_romanian_replay.py` pins it: stage 4, stage 5 and the validator report and mismatch CSV.
- **The new languages don't need the Romanian harvester.** They ingest kaikki.org's structured wiktextract extraction of Wiktionary instead of scraping HTML.

**R1. Protected Romanian files.** Romanian files are never modified:

- `scripts/romanian_harvester.py`, `dex_qc_main_csv.py`, `dex_supplement_derived_verbs.py`, `romanian_processor_main.py`, `attach_frequencies.py`, `download_leipzig.py`, `build_contingency_table.py`, `audit_verb_alternation.py`;
- `src/romanian_processor_lib.py`, `src/wiktionary_normalizer.py`, `src/dex_*.py`;
- `analysis/analyze_romanian_palatalization.R`, `run_pipeline.sh`;
- tracked Romanian analysis outputs.

Nothing under `lexpipe/` or the language plugins may import `wiktionary_normalizer`, `romanian_processor_lib`, `dex_utils` or `romanian_harvester`; an import-isolation test enforces this.

**R2. Leipzig files.** No file may be added to `data/leipzig/freq/`, because `attach_frequencies.py` attaches every `*_freq.csv` there. Non-Romanian corpora live under `data/<language>/leipzig/`. A guard test asserts that `data/leipzig/freq/` holds exactly `ron_news_2024_1M_freq.csv`.

## 2. Layout and identity

```
src/lexpipe/                    shared core (language-agnostic modules)
  kaikki.py frequency.py tp.py stats.py csvio.py sources.py spec.py      (stage 1, built)
  normalize.py licensing.py ids.py context.py attest.py learners.py    (stage 2)
  langs/italian/  langs/czech/  plugins (spec, lexicon, qc, features, derivations, analyses)
  langs/romanian_compare.py     read-only adapter that builds the harmonized Romanian table
src/phonology/inventories/{italian,czech}_features.json     feature inventories (tracked)
src/phonology/analyses/{italian,czech}_palatalization.py    provisional LP analyses
scripts/lang_pipeline.py        python scripts/lang_pipeline.py --lang it|cs [stage ...]
curated/<language>/*.csv        tracked hand decisions (overrides, adjudications), keyed by lexeme_id with guard fields
data/<language>/                generated (gitignored; NC-licensed derivatives under data/<language>/nc/)
analysis/<language>/            tracked summaries (only redistributable content, R9)
analysis/comparison/            cross-language tables (Romanian harmonized, Italian, Czech)
docs/<language>/                data dictionaries generated from schema; research notes
```

**R3. Plugin location.** Plugins live under `lexpipe.langs`, not a top-level `languages` package, to avoid namespace collisions. Language identity has one registry entry per ISO 639-1 code:

| Field | Italian | Czech |
|---|---|---|
| name | italian | czech |
| kaikki lang_code | it | cs |
| kaikki language name | Italian | Czech |
| Leipzig ID | ita_news_2024_1M | ces_news_2024_1M |
| file prefix | italian_ | czech_ |

All paths derive from this entry.

**R4. Identifiers.**

- `lexeme_id = {lang}:{word}:{pos}:{ety or 0}:{fp}`, where `fp` is the first 8 hex digits of a SHA-1 over the canonical JSON of `head_templates[*].args` plus `inflection_templates[*].args`. If two entries still tie, add the first gloss, then an ordinal. A uniqueness test runs over the full subsets.
- `(word, etymology_number, pos)` is **not** unique (490 Italian and 245 Czech collisions), and `etymology_number` is mostly null and page-local.
- `root_id` merges POS conversions and gender pairs of one etymology (amico N/ADJ, amico/amica), but keeps different etymologies apart (antico 'ancient' vs àntico 'anterior').
- TP counts **root types per cell**. Learner evidence pools all forms of a `root_id`.
- Every row carries identifying fields (gender/animacy, first gloss, raw IPA).

**R5. Determinism.**

- No first-wins dedupe; merge by an explicit priority.
- No calibration from the first N rows.
- Every pipe-joined list is `sorted()` or kept in source order (documented per field).
- Output rows are sorted by `lexeme_id` (then cell).
- Seeds come from the run context.
- Each stage has a determinism test: two `PYTHONHASHSEED` values plus shuffled input give identical output hashes.

**R6. Stage contracts.**

- Each stage is a function of a `RunContext` (language spec, paths, root, seed, source manifest). No module globals.
- Each stage writes outputs atomically, with a sidecar manifest recording input SHA-256s, parameters, a code fingerprint (SHA-256 over the shared core and the language's plugin source), git SHA plus dirty flag, and output SHA-256s. A stage is skipped when its manifest matches, so editing pipeline code re-runs the affected language's stages.
- `--force STAGE` re-runs derived stages only. It never touches `data/sources/` or `curated/`.
- `--root DIR` lets tests run in a temp tree.

**R7. Big sources.** Each big source gets exactly one streaming pass that writes a compact, indexed derivative stamped with the source SHA-256:

- kaikki: per-language projected JSONL or SQLite holding only the used fields;
- MorfFlex: a SQLite table restricted to NN rows plus the AAMP1 cell, Wiktionary lemmas, and forms in the Leipzig vocabulary;
- DeriNet: noun-parented derivation edges.

Later stages read only these derivatives. Measured costs: a full MorfFlex pass is about 164 s (single xz block); holding the full Italian subset in memory costs about 4 GB.

**R8. Normalization** (`lexpipe/normalize.py`, per language).

- Orthography is NFC with **no accent or length folding**: Czech á/é/í/ó/ú/ů/ý are kept, so diplomaté ≠ diplomate, and Italian final accents (città) are kept. Wiktionary respelling accents on non-final vowels (fòcora) are stripped for matching keys only.
- IPA accepts both `/…/` and `[…]` and records `ipa_kind`. ASCII `:` maps to `ː`. The raw IPA is kept with stress and syllable dots, and stress is computed from it.
- The phonology tokenizer and compare-strip sets come from the profile. Romanian defaults stay unchanged.

**R9. Licensing.**

- Every row/cell records `source` and `licence_class`. `lexpipe/licensing.py` maps each source to its licence and to whether it may be redistributed in the repo.
- Tracked outputs are written only through a gate that rejects NC sources: MorfFlex, DeriNet and czetyl are CC BY-NC-SA.
- NC-bearing outputs, including worklists that copy MorfFlex forms, go to `data/<language>/nc/` (gitignored).
- Leipzig-derived counts are not tracked until Leipzig's terms are checked by hand. Tracked summaries carry aggregates only.
- Treccani and IJP adjudications are stored as coded values plus a citation, never as copied entry text.
- `DATA_LICENSE.md` is a draft for the authors to review: code MIT; tracked data CC BY-SA 4.0 from Wiktionary via wiktextract (Ylonen 2022); other sources as used.

**R10. Booleans.** CSV booleans are written `True`/`False`, matching the actual Romanian output; `docs/data_dict.md`'s 'lowercase' claim is wrong. Readers are case-insensitive.

## 3. Sources and the "DEX" substitutes

| Role | Italian | Czech |
|---|---|---|
| Primary lexicon | kaikki raw en.wiktionary dump, top-level `lang_code == "it"` | same, `lang_code == "cs"` |
| Second edition | kaikki it.wiktionary extract (Wikizionario: hand-entered plurals, gender) | kaikki cs.wiktionary extract (more nouns; labelled variants) |
| Automated QC | Morph-it! 0.4.8 (CC BY-SA 2.0 / LGPL; from the GitHub mirror, because the official download is behind a human-verification page) | MorfFlex CZ 2.1 (CC BY-NC-SA 4.0) |
| Manual adjudication | Treccani, De Mauro: by hand only (Treccani opts out of TDM and plants trap URLs) | IJP: by hand only (robots.txt disallows `/?slovo=`; copying needs ÚJČ consent) |
| Derivations | Wiktionary etymology templates (`af`/`suffix`/`surf`) reverse-indexed + derived lists | DeriNet 2.3 (structure; verbalizer from SEGMENTATION) + Wiktionary templates |
| Etymology / stratum | Wiktionary `inh`/`bor`/`lbor`/`der`/`af`; "By surface analysis" = borrowed whole | Wiktionary + czetyl + Czech Wiktionary inheritance notes |
| Frequency | Leipzig `ita_news_2024_1M` (elided clitics split) | Leipzig `ces_news_2024_1M`; CNK SYN2015 lemma list as a cross-check |
| IPA / stress | lemma IPA (`it-pr`) from kaikki; plural IPA derived from the lemma; it.wiktionary accented hyphenation as a stress cross-check | `cs-IPA` (template-generated: not independent evidence) + rule fallback |

**R11. Form provenance.** Dictionary-generated forms are not attestation.

- Every cell form records `form_source ∈ {explicit, module_default, cs_ndecl_default, …}`.
- Attestation is a token count in the language's Leipzig corpus, both case-sensitive and folded. Italian elided articles and prepositions are split first (l'amico, dell'amico; both ' and ’).
- Forms that are also forms of another lexeme (carichi, also a verb) are flagged `form_ambiguous`.
- Every learner-evidence and TP table is reported over **τ ∈ {dictionary, ≥1 token, ≥0.5 per million, ≥2 per million}**.
- Morph-it! ADJ plurals and MorfFlex cells are generated, so agreement with them is QC, never confirmation of use.
- A POS-tagged count (UDPipe/Stanza over the Leipzig sentences) is a documented future improvement.

**R12. Italian plural provenance.**

- Build `italian_plural_forms.csv` from `head_templates` arguments, not from `forms`. Its columns: `lexeme_id`, `form`, `cell` (pl, m.pl, f.pl), `provenance` (explicit | module_default), `qualifiers`, `raw_arg`, `tag_error`.
- Parse the arguments 2/pl/pl2…/mpl/fpl, the codes `+ # ? ! - ~`, and the inline modifiers `<q:> <qq:> <l:> <ll:> <ref:>`.
- Port `Module:it-common` `make_plural` and run a drift test over every pure-default noun.
- `plural_qc ∈ {confirmed, variant, module_default_only, conflict}` comes from the per-form table plus Morph-it! and the it-edition.

**R13. Czech cells.**

- Cells come only from `forms[]` entries with exactly one case tag and one number tag, mapped to NOM.SG … INS.PL. Variants stay unranked.
- The cs-ndecl spec is kept for audit only, because it omits module defaults.
- The head code is parsed into gender, animacy, pluralia tantum and by-sense; comma lists give both readings plus an `animacy_ambiguous` flag.
- Missing paradigms (plech) get `paradigm_source = none`. They may be filled from MorfFlex only in NC files.
- Dual remnants (ruka) and the kuře type are out of domain. Syncretic forms are recorded once per cell.

**R14. MorfFlex and DeriNet.**

- `parse_pdt_lemma(raw)` follows TR-64 (base, homonym number, name labels, style labels `_,x`, comment). The raw lemma is the join key.
- Keep `NN` rows with tag[10] = 'A'. Var '-' and 1–4 are standard; 5–9 are non-standard; a/b/c are abbreviations and dropped. Style-labelled and proper-name lemmas are flagged.
- DeriNet: the verbalizer comes from the last derivational morph in SEGMENTATION (knih|ova|t → -ovat), with ConjugClass only as a consistency check. Keep `Type=Derivation` edges with a noun parent. Join on (lemma, gender, animacy), because šach and jazyk are each two nodes.

## 4. Shared analytic definitions (all three languages)

**R15. Evidence sets.**

- **E0:** nominal inflection.
- **E1:** E0 plus verb cells with a front-vowel ending:
  - Romanian: 2sg and 3sg subjunctive of -a and -i verbs;
  - Italian: 2sg, 3sg and 1pl of every conjugation;
  - Czech: -it/-et verbs and participles.
- **E2:** E1 plus derivatives.

A root counts as "seen alternating" relative to a stated evidence set and τ.

**R16. Harmonized TP domain** (`lexpipe/tp.py`). For each target class × trigger cell, report:

- `N_target`: all target-final roots;
- `N_trigger`: roots attested before the trigger;
- `e_nonalt`: N_trigger roots that never alternate;
- `n_shielded`: target-final roots whose every attested trigger-cell realization uses a non-front or N allomorph;
- θ(N_trigger), plus θ(N_target) as the counterfactual.

Further rules:

- Outcome-conditioned exclusions are replaced by outcome-blind flags (stratum, suppletion type).
- NDE classes are separate rows.
- Non-trigger front-vowel cells are emitted as N = e rows in every language (Romanian t+-e; Czech -é, -em, -ech, coronal VOC -e; Italian feminine -e).

**R17. Cross-language rule family.** One fixed rule family, applied at a morpheme boundary, is the comparison set:

- dorsal / _[+front V];
- dorsal / _[+high,+front];
- coronal / _[+high,+front].

Language-specific analyses (Italian height-indexed, Czech suffix-indexed) are reported separately as model fits, never as the comparison rows.

**R18. Two-outcome evaluation.** When both outcomes occur in a cell, report e and θ for **both** rules (palatalize and retain), pooled and in nested subclasses evaluated in Yang's order (specific rule first, then the elsewhere rule). Also report:

- **simulated learners:** about 1,000 vocabularies drawn by token frequency for V ∈ {500, 1000, 2000, 5000}, reporting the share of simulated learners with palatalize productive, retain productive, neither, or the conditioned rule;
- **fragility:** the minimum number of item changes that flips each verdict.

**R19. Doublets.** Classify each doublet by corpus ratio. At ≥ 90% for one variant, code the dominant variant; otherwise it is a true doublet. Report the TP under three codings: dominant variant, doublet counted as an exception, doublet excluded.

**R20. Spell-out TP on certain-specification roots.**

- Roots: those with unambiguous alternation evidence (Czech hard feminine velars; Romanian alternating-plural roots; Italian roots with attested -ci/-gi).
- N: roots with an attested derivative from a competing pair:
  - Czech: -ní/-ný vs -ový/-ovní; -it vs -ovat; -ka vs -yně;
  - Romanian: -i vs -a/-ui;
  - Italian: -izzare vs -eggiare/-are/-ire; palatalized vs velar -ista/-ismo/-ità.
- e: roots choosing the non-palatalizing member.
- Stratify by base stratum.

**R21. No-evidence roots** (the plural-less prediction).

- "No alternation evidence" means no attested (≥ τ) alternating form in any inflectional cell. It is computed from the corpus, not from dictionary paradigms, because dictionaries generate unused forms (vzduších, snězích, biolozích: 0 tokens).
- Derivative candidates are generated in both spellings (leghista/legista; vzdušný/vzduchový) and looked up in headwords and the corpus.
- Report 3×2 tables per suffix and stratum: base alternating / non-alternating / no evidence × derivative palatalized or not.
- Italian: feminine-only -ca/-ga stems and uncountable -co/-go stems, with the -ic- suffix and masculine counterparts reported separately (musica / musico).

**R22. Stratum.** Both base and derivative get a stratum ∈ {inherited, borrowed whole, internal formation, unknown}. Every co-variation, spell-out and no-evidence table is reported by stratum, with unknown kept separate. Derivational suffixes get their own trigger rows; Czech -ista/-ismus are non-trigger rows, parallel to Romanian -ist.

**R23. Co-variation statistics.**

- The outcome is coded at the base × suffix level: only P / only N / both.
- Report per-suffix odds ratios, a Mantel–Haenszel pooled OR, and a permutation of the base class within suffix × stratum strata (seeded).
- Report a logistic model `deriv_pal ~ base_class + log_lemma_freq + base_stratum + deriv_stratum (+ stress)` with base-clustered errors. A mixed model with base and suffix intercepts is reported if available (R lme4).

**R24. Sampling frame.**

- One frame for all languages: all nouns of any final segment and plural type, ranked by lemma frequency.
- Report both frequency bases: lemma form, and paradigm sum over attested forms.
- Counts are case-sensitive, with a sentence-initial correction where possible.
- Sweep V ∈ {500, 1000, 2000, 5000}.

**R25. Homonyms and names.**

- Ethnonyms and demonyms are common nouns and stay in (Češi/Čechové carry much of the -i/-ové variation).
- Given names, surnames and toponyms are out.
- Different etymologies are separate roots.
- Sense-conditioned allomorphy within one etymology (jazyk 'language' v jazyce vs 'tongue' na jazyku) is variation within one root, with a sense-split sensitivity check through MorfFlex lemma splits.

**R26. Suffix-internal targets** are identified morphologically, not orthographically:

- Italian -ico: an "X + -ico" template, or a Latin -icus adjective with an attested Italian base and adjectival category. Never deverbal -ico nouns (carico) or simplex Latin nouns (medico).
- Czech: -ník, -ík, -ák, -ek, -ka, -ička, -ský/-cký, -log, -urg, -graf.

Each suffix is its own row and counts once as a morpheme in learner evidence. Carico-type items are listed in both root and suffix tables.

**R27. -log/-fag/-urg table.** One -log/-fag/-urg table per language, with the same person/thing classifier, cell outcomes, attested frequencies and derivatives.

## 5. Italian specifics

**Opportunity.**

- `i`: masculine plural, feminine -o → -i, masculine -a → -i.
- `e`: feminine plural, the non-trigger cell; categorically non-palatalizing, so never pooled with `i`.
- `a`: irregular -a plural (back vowel).
- `none`: invariable, missing or suppletive.

**Mutation.**

- Palatalizing (undergoer): `co→ci go→gi ca→ci ga→gi sco→sci cco→cci`.
- Velar kept: `co→chi go→ghi ca→che ga→ghe sco→schi cco→cchi`.
- Variable when attested plurals compete (R19).

**Classes.**

- `frontstem` (lexically palatal): -ce/-ge/-cia/-gia/-cio/-gio/-scia/-scio.
- `nde:occhio`: -chio/-chia/-ghio/-ghia.
- `nde:psiche`: singular -che/-ghe.
- `nde:gimpe_like`: singular = plural with C + front vowel.
- Root-internal chi/che/ghi/ghe: an NDEB count.
- Invariables and unadapted loans.
- Archaic, Old Italian, literary and proscribed variants are filtered by tag or qualifier and reported.

**Stress.**

- `velar_post_tonic` (the variable that matters: Giavazzi; Shao et al.) and `stress_class` are computed from the IPA syllable dots and ˈ.
- Heavy penults are a tendency, not a rule: the -urgo persons and belga palatalize, and polizza and mandorla are antepenultimate.
- Report missing-stress counts per cell. Keep one row per stress only when the plural differs (antico/àntico).

**Verbs** (TP rows by conjugation; E1 evidence).

- -are: never palatalizes.
- Root-stressed -ere: alternates.
- -ére plus cuocere and nuocere: over-apply.
- -ire: mixed.
- -isc-: never.
- dire: suppletive.
- Report a pooled noun + verb row and a category-indexed row.
- `deverbal_from_are` (the etymology says "Deverbal from" or "short past participle of X-are") is kept separate from `shares_root_with_1st_conj_verb` (synchronic).

**Derivatives.**

- Italian-internal learned formations (-issimo, -izzare, -ista, -ismo/-esimo, Italian-internal -ità) are the co-variation set.
- Native suffixes (-ezza, -etto, -ino, -eria, -essa, -evole, -eggiare, derivational -ire) are the never-palatalizing control.
- Learned borrowings are excluded by stratum.
- Fact-check corrections: grecizzare (from Grecia), politichese (from politica) and musicista (from musica) are not plural-tracking evidence; porcino, porcello and porcile are inherited.

## 6. Czech specifics

**Targets.**

- Core targets: k, h, g, ch, r (`ch` is one letter), with clusters sk, ck, st, zd, Cr.
- t/d/n: their palatalization is invisible in spelling, so they are flagged `coronal_inferred`. They are excluded from rule TP (a) and from learner evidence (c), and analyzed in (b) as allomorph competition (-i vs -é; -ě vs -u; -ech).
- Recorded separately: c (first palatalization only), s/z (big pattern only), labials (glide insertion).
- Never a target: l.

**Palatalizing (P) vs non-palatalizing (N) endings, by cell:**

| Cell | P | N |
|---|---|---|
| NOM.PL masc. animate | -i | -ové, -é |
| LOC.SG masc. inanimate / neuter | -e/-ě | -u, -i (soft) |
| DAT/LOC.SG feminine | -e/-ě | (none in hard paradigm) |
| LOC.PL | -ích | -ách (the default for velar neuters), -ech |
| VOC.SG masc. | -e for velars (člověče, bože), r after a consonant (bratře), -ec → -če | -u, coronal -e, r after a vowel |
| INS.SG | (none) | -em (non-trigger front vowel) |
| Hard adjective NOM.PL.M.AN | -í | (separate domain) |

**The three TP analyses.**

- **(a) Rule TP.** N = roots with ≥ 1 attested P-ending form in the cell; e = those without mutation. Velars, r and clusters only. Non-trigger rows (-é, -em, -ech, coronal VOC -e) are N = e.
- **(b) Competition TP**, per cell: LOC.SG -e/-u; NOM.PL -i/-ové/-é; LOC.PL -ích/-ách, split by gender; VOC; adverb -e/-y.
  - Evaluate both the P-default and the N-default rule within nested classes: gender × animacy × final class × suffix/monosyllable/loan.
  - Use the dominant variant by token frequency.
  - Include a baseline of labial-, l-, s- and z-final hard stems.
- **(c) Spell-out TP, leave-one-cell-out.** Specification is inferred from attested forms in a pre-registered evidence set that excludes the tested cell. N = roots inferred underspecified that are attested in the tested cell; e = those taking the N allomorph.

**Lexeme summary.** `undergoer`, `resister` or `mixed`, relative to the evidence set and τ. Mixed stems (biologové / biolozích / bioložka) are the decisive cases for spell-out by root specification vs lexical class.

## 7. Romanian harmonized comparison (read-only)

`lexpipe.langs.romanian_compare` reads `data/romanian_lexicon_with_freq.csv` (and `data/verb_paradigm_verdict.csv` for E1) and computes, with the same code as Italian and Czech: the harmonized TP rows (R16–R19), the spell-out TP (R20), the no-evidence tables (R21) and the -log table (R27). It writes only to `analysis/comparison/`.

## 8. Testing and documentation

- **Fixtures:** small CC BY-SA kaikki fixtures (with an attribution README) under `tests/fixtures/`. Italian: amico, antico ×4, fuoco, mago, chirurgo, greco ×3, belga, valigia, luce, medico, parco, lago, teologo. Czech: kluk, žák, vrah, Čech (noun and name), biolog, jazyk, park, šach ×2, plech, ruka, diplomat, student, kniha, noha.
- **Golden expectations:** amici / medici / parchi / laghi; chirurghi ~ chirurgi is variable; amiche is the non-trigger -e cell; luce is frontstem; kluci is an undergoer; biologové + biolozích is mixed (dictionary) and resister (τ ≥ 1); plech has no paradigm.
- **Invariants:** lexeme_id uniqueness; Czech rule-TP e ≈ 0; coronal_inferred rows are 100%.
- **Other tests:** parity with the Romanian R and contingency numbers, determinism, licence gate, import isolation, tracked-path guard (`git check-ignore`), the Leipzig freq-directory guard, and schema headers.
- **Markers:** `replay`, `needs_data`, `slow`, `ro`, `it`, `cs`.
- **Schema:** a per-language `schema.toml` (column, type, allowed values, provenance class, licence class). It validates writers and generates `docs/<language>/data_dict.md`.
- **Runbooks** in `docs/runbooks.md`: re-baselining Romanian; moving to a new kaikki dump.
- **README and CITATION.cff** gain data attributions: Ylonen 2022; Zanchetta & Baroni 2005; MorfFlex (Hlaváčová et al.); DeriNet (Vidra et al.); Křen (SYN2015); Goldhahn et al. 2012.

## 9. Actions for the authors (cannot be done by the pipeline)

- **Commit the re-baseline.** The working-tree analysis outputs (`analysis/palatalization_*`, `analysis/romanian_tp_summary_*`) are the reproducible baseline the replay test pins. HEAD's copies are stale.
- **Commit `romanian_features.json`.** It is needed by the tests and was hidden by `*.json`.
- **Review `DATA_LICENSE.md`** and the Leipzig terms before tracking any Leipzig-derived numbers.
- **Morph-it! official tarball (optional):** download in a browser to `data/sources/morphit/morph-it.tgz` to replace the GitHub mirror.

## 10. Open questions for the authors

- **Kramer 2009:** the CILT 304 chapter (two groups of speakers) or *The Phonology of Italian* §4.2.1 (10 speakers: 4/1/5, per the summary on Kyle Gorman's *Wellformedness* blog)? The full texts were unreachable in this session.
- **Italian 1st conjugation:** do these verbs count against a single dorsal rule (pooled row), or does a category/theme-vowel condition apply (category-indexed row)? Both are reported.
- **Mixed Czech stems:** if common, should -ové/-u spell-out refer to lexical class rather than root specification?
- **Czech loan-boundary coronals** (studenti [ɲcɪ]?) lack orthoepic confirmation (Hála 1967; Hůrková 1995).
- **Paper vs code:** /Z/ and /T/'s Continuant value differ between `latex.tex` and `analyses/romanian_palatalization.py`. Pin these before cloning the analysis.

## 11. Corrections since adoption

The implementation (stage-1 core, the Italian, Czech and Romanian plugins, the core follow-up and the migration onto the shared core, 2026-10-23) deviates from the rules above in the places listed here. Each entry names the rule, what the code does and whether it stands (a documented convention) or was corrected by the migration. Code reviews should cite this list together with the rules. The migration's effect on every tracked output is in `data/research_v2/migration_diff.md` (local, gitignored).

**Conventions that stand (documented in the code):**

- **R16, extra `e` column.** The long TP rows carry `e`, the exceptions to the row's own rule, besides `e_nonalt`: with doublets and two rules (R18) the exceptions cannot be recovered from `e_nonalt` alone (`lexpipe.tp`).
- **R16/R18, tracked TP columns.** Every tracked long table (Italian, Czech, Romanian) adds `fragility`, `source` and `licence_class` to the R16 columns: `lexpipe.tp.TP_LONG_TRACKED_FIELDS`, written by `lexpipe.tp.tracked_row`.
- **R11, the `quest` clitic.** `ITALIAN_ELISION_CLITICS` adds *quest* to the elided clitics the design lists (quest'anno: 7,023 tokens in ita_news_2024_1M).
- **R8, Czech IPA stays untied.** The Czech normalizer does not tie adjacent t+ʃ, d+ʒ, t+s, d+z (Italian and Romanian do): cs-IPA writes the affricates c and č tied already, so its untied pairs are clusters across a prefix boundary (pod+zim, od+stavec) that tying would merge (`lexpipe.normalize`).
- **R26, Italian -ico strictness.** An -ico suffix target needs an affix or surface template, or a borrowed Latin/Greek -icus etymon on an adjective. R26 also asks for an attested Italian base; the base is recorded (`suffix_base`, `target_is_suffix_strict`) but not required, because Wiktionary rarely states it (requiring it would drop pubblico, a golden target). The suffix class is a lower bound.
- **R24, Romanian frequencies are lowercased.** The Romanian Leipzig table (R2, protected) has folded counts only: no case-sensitive counts and no sentence-initial correction; its paradigm basis is lemma + indefinite plural. Italian τ flags also use folded counts, for comparability with Romanian (case-sensitive counts only code doublets); Czech uses case-sensitive counts plus the sentence-initial capital correction. Cross-language τ rows therefore differ in this respect.
- **R11/R18, learner thresholds.** Simulated learners run at τ ∈ {dictionary, ≥1 token} only (Italian and Romanian E0/E1, Czech E0), not at all four thresholds; the TP rows themselves cover every τ.
- **R22, Romanian strata.** The Romanian lexicon records only an etymon language: its strata are `borrowed`, `latin_unresolved` (inherited or learned Latin) and `unknown`; internal formations cannot be detected.
- **R15, Romanian E1.** E1 reads the verb cells orthographically (the vowel after the root) rather than the audit's verdict column; they disagree on 63 of 422 pairs (`romanian_parity.txt` §6).
- **R3, Romanian registry row.** `ro` is a registry row with `comparison_only=True` (the R3 table lists Italian and Czech only); comparison-only rows have no normalizer, docs, curated or raw-lexicon paths and their contexts may write only under `analysis/comparison/` and `data/romanian_compare/`. Its curated -log decisions live in `curated/romanian_compare/`.
- **R9, aggregates of NC sources in tracked text.** Czech text reports that count DeriNet pairs, czetyl/DeriNet strata or MorfFlex agreement (`czech_qc_report.txt`, `czech_derivations_report.txt`, `czech_spellout.txt`, `czech_noevidence.txt`) are tracked as aggregates under the draft policy of `DATA_LICENSE.md`: written with `RunContext.write_tracked_text(..., nc_aggregates=True)`, which refuses unverified sources and prepends a licence and attribution header. Tables never carry NC rows. Pending the authors' decision.
- **R9, Leipzig as a threshold.** Tracked aggregates use Leipzig counts as attestation thresholds without naming `leipzig` as a source (the gate refuses unverified terms): the reports say so in their text. The Italian worklist's `freq_rank` is an ordinal derived from Leipzig counts.
- **R9, the Czech worklist.** Which lexeme/cell pairs the tracked Czech worklist lists is chosen by comparison with MorfFlex (CC BY-NC-SA); its rows hold Wiktionary forms and IJP links only.
- **R12, plural table columns.** `italian_plural_forms.csv` adds `code`, `gender`, `labels`, `item` and `source` to the R12 columns; `plural_qc` keeps its four values, `qc_basis` tells an independent confirmation from an explicit-only form.
- **§2, tracked raw lexicons are lean.** `data/<language>/<prefix>lexicon_raw.csv` holds extracted columns only (identifier, lemma, part of speech, gender and animacy, first gloss, raw IPA, the head-template arguments as written, the cells as extracted for Czech, a summary of the deciding etymology template, source) and leaves out out-of-domain proper names and multiword expressions (under 10 MB each); the full stage table with every derived column is `data/<language>/<prefix>lexicon_full.csv` (gitignored), which later stages read.
- **R5, list cells.** A `|` or backslash inside a list item is backslash-escaped (`lexpipe.csvio`); cells without one read as before.

**Corrected by the migration onto the shared core:**

- **R22, one stratum classifier.** Italian and Czech each had a copy of the template parser and classifier. `lexpipe.etymology` now holds one (union of both template tables; Czech artifact skipping and `<surf>` kinds; `:from` steps classified by language; internal kinds whose text opens with a borrowing are borrowed whole, now for Czech too; "By surface analysis" makes an internal kind borrowed whole only when its step names no native base, generalizing the Italian antichità fix; etymology trees stripped from the text; parenthesised parts of a term stripped at a word edge, `(nux) Abellāna` → `Abellāna`, and kept inside a word, `*ga(v)ota`, where Czech stripped every parenthesis and Italian none).
- **R27, one classifier pair.** The person/thing classifier and the learned-suffixoid detector are `lexpipe.suffixoids` for all three languages, with the manual lists moved to `curated/<language>/suffixoid_classes.csv`. Italian had its own gloss rule (classes person, thing, mixed, adjective); it now uses the shared one (adjectives keep the Italian class `adjective`). Czech now flags lookalikes (profil). The shared thing rule matches *cave* as a whole word only (as a substring it made *scavenger* a thing), and *scavenger* is animate.
- **R23, one model.** The logistic model was `base_class + log frequency + deriv_internal + base_borrowed` (Italian, never estimable: a constant predictor) and `base_class` alone (Czech). All three languages now fit `lexpipe.stats.r23_logistic` (`deriv_pal ~ base_class + log_lemma_freq + base_stratum + deriv_stratum`, strata treatment-coded, aliased predictors dropped as R does, a stratum factor with one level in the units reported as having no term, quasi-complete separation reported). Romanian gained the R23 analysis (plural class × verbalizer choice).
- **R23, permutation ties.** The stratified permutation p-value depended on whether statsmodels computed the observed MH OR (floating-point ties); null values within a relative 1e-9 of the observed one now count as ties.
- **R9, gated text.** Every tracked text report goes through `RunContext.write_tracked_text`; the Czech reports with NC numbers were written ungated.
- **R9, Romanian sources.** Romanian rows were filed as `wiktionary_en`; aggregates now name `wiktionary_ro|wiktionary_en` and the -log table each row's edition. The -log table had also filed as Wiktionary the nouns whose plural DEX QC filled (lexicon `source` ending in `| dex_qc`: dexonline, terms unverified); it leaves them out now (final review, 2026-10-23). The aggregates still count them, as the legacy outputs do (`DATA_LICENSE.md`).
- **R5, '|' workarounds removed.** Czech IPA phrase transcriptions (dropped) and glosses or head links containing `|` (replaced or dropped) are kept, escaped.
- **R6, curated cell decisions.** The Czech cell adjudications were a wide template (one row per lexeme); they are keyed by (lexeme_id, cell) now.
- **R18, one draw per learner set.** Each language evaluated its learner cells on separately redrawn vocabularies (Italian and Czech had local forks of the core); every cell is now evaluated on the same draws by `lexpipe.learners.simulate_all` (identical shares, Romanian 15 min → 1 min).

**Added by the audits of the Italian and Czech datasets (2026-10-23).** These are conventions that stand; each language's README and data-layer notes give the details.

- **R19, doublets within one cell (core option).** `lexpipe.tp` has two opt-in refinements:
  - `Observation.cell`: a doublet is competing variants of one paradigm cell, so a root whose cells disagree is an exception to both rules under every coding. An example is amici before -i but amiche before -e in R17a, or the plural medici but the verb form medichi in E1.
  - `CellSpec.doublet_min_tokens`: a doublet below the minimum is never "dominant".

  Italian uses both (cell keys, 10 tokens). Czech and Romanian use neither, so their R17a rows code a root whose cells disagree by token majority. Romanian has one plural per lemma, so its E0 rows are unaffected.
- **R11/R25, Italian names filter.** τ thresholds stay on folded counts (comparable with Romanian), but a form counts as attested only with at least one lowercase token. Surnames and firm names (Bocchi, Monte dei Paschi) no longer attest velar plurals.
- **R11, Italian `dictionary_explicit`.** 96.5% of the dictionary's dorsal -i plurals are `Module:it-common` defaults. The stress rows are also given over the plurals an editor typed. A form listed for two lexemes of a root (nemici, of nemico N and ADJ) counts once with its τ flags merged, so a plural typed for one lexeme counts even when the other's copy is a module default (final review, 2026-10-23: keeping the first copy had dropped nemico from the typed post-tonic row, N 69 → 70).
- **R15/R22, Italian E2.** E2 holds derivatives that are internal formations or of unknown stratum (§5: learned borrowings are excluded). `E2_all` adds the borrowed-whole and inherited ones.
- **R26, "counts once".** The Italian learners and nested stress sets are also given with each R26 suffix counted once (`suffix_once`) and with every unstressed -ico/-ica root as one unit (`orth_ico_once`, an upper bound). The cross-language learner cells of all three languages count every suffix target as a root of its own, so R26's rule is applied only in these Italian sensitivity cells. The Italian learner result depends on it (`analysis/comparison/README.md`, claim (g)).
- **R11, Czech strict thresholds.** `ge1_strict`, `ge0.5pm_strict` and `ge2pm_strict` count a form's tokens only when the form is neither another lexeme's (`form_ambiguous`) nor shared with a cell outside its copy group (LOC.SG parku = GEN.SG = DAT.SG). The rule, competition, spell-out, learner and R21 rows have them. Without a tagged corpus, the ≥ 1-token LOC.SG and VOC -u rows are upper bounds.
- **R24, Czech frequency bases.** There are three bases: lemma, paradigm sum, and paradigm sum without forms of other lexemes. Homograph lemmas (jak 'yak') carry another word's tokens. The Czech learners and V sweep run on all three; Italian and Romanian learners use the lemma basis.
- **R26, Czech suffix targets.** They are identified morphologically: every -ský/-cký adjective, derivational form-of senses, headword lines, and decomposition against an existing base. A -ský/-cký nested class comes before the etymological ones.
- **R22, Czech adapted loans and hybrids.** "From Latin tragicus + Czech -ický" stays borrowed whole (`adapted_loan`). A Czech formation on a borrowed stem (fanoušek, Velšan) is an internal formation (`hybrid`).
- **R4, Czech derivational form-of senses.** Diminutives, female equivalents and verbal nouns are lexemes (`derivational_formof`), since they have paradigms of their own.
- **R9, word lists.** Counterexample lists in tracked Czech reports are printed at τ = dictionary. Lists filtered by Leipzig attestation go to the untracked `data/czech/`. The examples of `italian_qc_report.txt` and the Italian worklist are ordered by Leipzig frequency rank, an ordinal only.

**Cross-language synthesis (2026-10-23).**

- **Where it runs.** `lexpipe.compare` (`PYTHONPATH=src python -m lexpipe.compare`) builds `analysis/comparison/crosslang_summary.csv` and the generated block of `analysis/comparison/README.md`. That block is delimited by markers like the Romanian section's; the prose around both is kept byte for byte. It is a separate entry point, not a stage of the `ro` plugin: a stage there would need Italian and Czech outputs that the Romanian fixture tree lacks.
  - It runs in the `ro` context, so it may write only under `analysis/comparison/` and `data/romanian_compare/`.
  - Its manifest is `crosslang.manifest.tsv`; its parameters carry the module's SHA-256, so a code change re-runs it.
  - `--check` exits 1 when the tracked outputs are stale.
- **Inputs.** It reads tracked outputs only: the three TP long tables, the learner tables, the Czech -log table and the text reports. Report parsing is strict (`ReportFormatError`), and the NC sections of Czech reports are never read.
- **What it computes.** Every count is copied. The step computes only:
  - the sample odds ratio and Fisher p of a 2×2 summed over pairs or suffixes (`lexpipe.stats`; for Italian the same numbers as the Italian README);
  - θ and verdicts where a source gives counts only.
- **Matched settings** (`matched`): E0; τ ∈ {dictionary, ≥ 1 token}; dominant coding; lemma basis; nouns; the full lexicon; learners at V ∈ {1000, 5000}. "Nouns" is the Italian `nouns` subset; the Italian headline subset `all` adds the -ico adjectives and is reported as a sensitivity.
- **Comparability notes.**
  - The Czech R17 rows count every front-vowel realization whatever its allomorph, so they carry no shielded roots; Czech shielding shows in the rule TP (a) rows.
  - The Romanian `all` subset keeps the outcome-blind flags, which `core` drops.

**Final hardening (2026-10-23).**

- **R6, code fingerprint in every manifest.** Stage manifests used to compare inputs and parameters only, so a code change without a hand-bumped `VERSION` parameter left stale outputs in place (the final review found the Czech QC worklist stale this way). `lexpipe.context.code_fingerprint` now hashes `src/lexpipe/*.py` plus the language's plugin package into a `code` manifest row, and a manifest without that row counts as stale. After the change every stage of it, cs and ro was re-run from scratch: all 57 tracked outputs came out byte-identical, and a second plain run skips every stage.
- **R9, the validation driver.** `scripts/validate_lang_rules.py` checks that its sources (English Wiktionary and the language's own edition) are redistributable and writes through the language's `RunContext`, so its write guard applies.
