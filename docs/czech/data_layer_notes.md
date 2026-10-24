# Czech data layer: notes

Notes on `src/lexpipe/langs/czech/` (design v2, `docs/multilingual_design.md`
§2, §3, §6). The module docstrings are the reference; this file records
the decisions taken while building it and the numbers of the first full
run (2026-10-23).

Run: `python scripts/lang_pipeline.py --lang cs extract lexicon qc attest features`.

## Stages and files

| Stage | Writes | Tracked? |
|---|---|---|
| extract | `data/czech/czech_{en,cs}_entries.jsonl.gz` (projected lemma entries), `czech_form_index.sqlite`, `czech_extract_stats.tsv` | no |
| lexicon | `data/czech/czech_lexicon_raw.csv` (lean: extracted columns only, no proper names or multiword expressions; `write_tracked`), `czech_lexicon_full.csv` (every lexeme and column; read by later stages), `czech_cells_raw.csv`, `czech_lexicon_stats.tsv` | lean raw lexicon only |
| qc | `data/czech/nc/czech_morfflex.sqlite`, `nc/czech_morfflex_qc.csv` (NC); `analysis/czech/czech_qc_report.txt`, `analysis/czech/czech_manual_worklist.csv` | report and worklist |
| attest | `data/czech/leipzig/…`, `czech_attestation.csv`, `czech_lexeme_frequency.csv`, `czech_attestation_report.txt` | no (R9: Leipzig) |
| features | `data/czech/czech_cells.csv`, `czech_lexicon_complete.csv`, `czech_features_report.txt` | no |

Hand decisions: `curated/czech/cell_adjudications.csv` (header-only
template, one row per (lexeme_id, cell); contract in `adjudication.py`,
procedure in `curated/czech/README.md`) and
`curated/czech/suffixoid_classes.csv` (R27).

## Decisions

- **Population.** English Wiktionary lemma entries with pos noun, adj,
  name are lexemes, and so are entries whose sense is a
  **derivational** form-of (wiktextract tags "diminutive of kluk",
  "female equivalent of Čech", "verbal noun of kvasit",
  "augmentative of …" as form-of; `extract.derivational_relation`):
  klouček, Češka, kvašení have paradigms of their own, and many R26
  suffix targets (-ek, -ka, -ička, -ík) are among them. The relation is
  kept in `derivational_formof` (`diminutive:kluk`). Inflectional
  form-of senses, spelling variants, possessives of given names and the
  adjective degrees (comparative, superlative: forms of the adjective's
  own paradigm) stay non-lexemes. Czech Wiktionary noun/adjective
  entries are matched to English lexemes by word, part of speech and
  compatible gender and animacy; their cells become second-source rows
  (`primary = False`). Czech entries whose word has no English lexeme of
  that part of speech become lexemes of their own (cs-only), unless they
  are names.
- **`paradigm_source`.** `en` (English table), `cs` (cs-only lexeme),
  `none` (English entry without a table). plech is `none` although
  cs.wiktionary has its table: that table is kept as second-source rows,
  so a lexeme's evidence never mixes editions.
- **`root_id`.** Key (word, etymology number), plus the animacy when one
  etymology has both an animate-only and an inanimate-only masculine noun
  (šach 'shah' / šach 'check'; MorfFlex and DeriNet split them too).
- **Names (R25).** Proper names are pos `name`, a `cs-proper noun`
  head template (Inka 'a female given name' sits under a Noun header) or
  name-only glosses. Categories are page-level (the noun Čech carries
  "Czech proper nouns") and only classify a known name. cs.wiktionary
  tags ethnonyms `vlastní jméno` too; a demonym gloss (obyvatel…,
  …pocházející z…, …národnosti) keeps them in. A **capitalized** entry
  is a name too when its first gloss names a place, a sky object or a
  name type and is no demonym gloss (`names.capitalized_place_en`,
  `names.capitalized_name_cs`; evidence `capitalized_gloss`): the
  constellations Rak, Býk, Lev (sense category cs:Constellations in the
  zodiac), Franky 'Franconia', Pha 'Prague', Jesenicko; cs-only
  Studénka 'město v okrese …', Anuška '… ženské křestní jméno',
  Jindřišská 'ulice …', Azor 'jméno pro psa', Jupiter 'pátá planeta …'.
- **Cells (R13).** Exactly one case and one number tag; forms on the
  headword line are other lexemes; tables are segmented by wiktextract's
  `table-tags` pseudo-forms (agent: an + in), and a repeated
  (table, cell, form) is recorded once.
- **`form_source` (R11).** `explicit` when the cell's slot is overridden
  in the cs-ndecl spec (or copied from one: VOC.PL from NOM.PL), the
  table is manual (`!`, cs-decl-noun-dual) or the edition is cs
  (literal tables); else `module_default`.
- **Out of domain.** names, multiword, abbreviation, hyphenated,
  dual_remnant (ruka, oko, sto), irregular_stem (wiktextract class
  `irreg-stem`: ucho 'ear', člověk), tstem (kuře type), indeclinable.
- **MorfFlex projection (R7, R14).** Rows kept: affirmative NN rows and
  `AAMP1----1A----`, whose base is an in-domain Wiktionary lemma **or**
  whose form is in the Leipzig vocabulary. QC statuses: `present`,
  `styled_lemma` (standard form of a style-labelled lemma: `_,h`
  colloquial, `_,s` spelling variant), `nonstandard_only` (Var 5–9),
  `absent`, `cell_missing`, `no_lemma`. The worklist lists only `absent`
  and `nonstandard_only` forms in cells with a P/N choice, for lexemes
  whose stem ends in a velar, r or t/d/n. It holds Wiktionary forms and
  IJP links only; the qc status itself stays in the NC file.
- **Attestation.** `count` is the case-sensitive count plus the
  capitalized form's count for lowercase forms (sentence-initial
  tokens), unless that capitalized string belongs to another lexeme
  (Poláci). `form_ambiguous` uses the form index: owners other than the
  lexeme or in-domain lexemes of its root; spelling variants (`alt_of`)
  and cs.wiktionary's form-of rows of multiword phrases (kapr under
  "vánoční kapr") are no owners.
- **Frequency bases (R24).** `lemma_count` (lemma form) and
  `paradigm_sum` (attested forms) as R24 asks, plus `lemma_ambiguous`
  (the lemma string is another lexeme's form: jak 'yak' carries the
  conjunction's 33,277 tokens) and `paradigm_sum_unambiguous` (forms
  that are not `form_ambiguous`); the TP stage runs the learners and
  the V sweep on all three.
- **Summaries.** E0 only (nominal inflection); `*_strict` variants drop
  ambiguous forms and forms that also fill a cell outside their copy
  group (LOC.SG parku = GEN.SG = DAT.SG), whose tokens do not attest the
  cell's allomorph choice.
- **Orthography.** Before i/í/ě, Czech spells ť as t, so sk→šť and ck→čť
  show up as `št`/`čt` (čeští, němečtí); both spellings are matched.

## Full run after the audit fixes (2026-10-23)

Lexicon: 57,233 lexemes (English: 23,992 nouns, 7,931 adjectives, 9,740
names; cs-only: 9,452 nouns, 6,118 adjectives). In the domain: 30,430
nouns, 13,852 adjectives. Out: names 9,759 (surname 7,145, given name
944, toponym 1,305, other 365), multiword 2,650, indeclinable 373,
tstem 75, irregular_stem 28, abbreviation 17, hyphenated 46, dual 3.
1,162 English noun entries and 74 adjective entries with a derivational
form-of sense became lexemes (2,107 lexemes carry
`derivational_formof`: 1,208 diminutive, 718 verbal noun, 169 female
equivalent, 12 augmentative; some also have a lemma sense).

The first full run (before the audit fixes) had 56,442 lexemes (the
earlier text said 56,430 in error): English 22,830 nouns, 7,857
adjectives, 9,740 names; cs-only 9,870 nouns, 6,145 adjectives; in the
domain 29,696 nouns and 13,805 adjectives.

Run times: extract 36 s, lexicon 56 s, qc 16 s (plus 110 s when the
MorfFlex projection is rebuilt), attest 49 s, features 59 s.

The stage reports (`czech_lexicon_stats.tsv`, `czech_qc_report.txt`,
`czech_attestation_report.txt`, `czech_features_report.txt`) hold the
full tables.

## Open points

- E1 (verbs) and E2 (derivatives) evidence is left to the derivation
  stages; head-line derivatives (feminine, diminutive, relational
  adjective) are exported as `head_*` columns for them.
- LOC.SG -u counts are confounded with GEN.SG/DAT.SG -u (see `*_strict`).
  Attributing tokens to cells needs a tagged corpus (R11 future work).
- Coronal mutation is inferred (100% undergoer by construction).
- 39 primary velar/r/coronal rows (64 rows in all) are `unexplained`,
  nearly all errors in cs.wiktionary tables (metodickí, nepsychologi):
  QC flags, not exceptions.
- 39 in-domain nouns have `stem_method = genitive_check` (the common
  prefix of their stem cells disagrees with GEN.SG minus its ending):
  for manual review.
- `form_ambiguous` still counts spelling variants listed in another
  entry's table (rasizmus lists rasismus) as other owners.
- cs-only name detection relies on Czech glosses; 360 cs entries share a
  word with an English lexeme but not its gender, and are left out.

## Data-layer fixes made during the analyses (2026-10-23)

Both are in `ortho.py` (features `VERSION` 5), with tests in
`tests/czech/test_ortho_stems.py`; the analyses were rerun on the
corrected cells.

- **Stems with a quantity alternation.** The stem was the unfolded
  common prefix of the oblique forms, so dráhy ~ drahám gave dr- and
  the stem final was wrong (dráha r instead of velar, skála velar
  instead of lateral, kráva r instead of labial). The prefix is now
  taken with vowel length folded, capped by the unfolded common prefix
  of DAT.PL and INS.PL so that folding never runs into the endings
  (Franky: Frankám, Frankami). 23 lexemes changed class.
- **Soft declension without a GEN.SG form.** The declension was read
  only off the GEN.SG ending; cs.wiktionary tables without that cell
  and pluralia tantum came out `unknown`, and their stem final put
  soft nouns into the target classes (chovankyně, the -yně feminines,
  lázně, kleště: coronal_inferred). Without a GEN.SG form an INS.SG in
  -í or a lemma equal to the stem plus -e/-ě/-í now marks the soft
  declension. A stem that contains the whole -í lemma (Polabí, whose
  only stem cell is INS.SG Polabím) is cut back to the lemma minus -í;
  it had made 481 -í neuters "labial" (final m). 735 lexemes changed
  class, 650 of them in-domain nouns (481 labial, 144 coronal, 7 velar
  and 3 r → soft; 15 to lateral, c or sibilant). The velar and r
  classes lost only the loans cache, image, cringe, bridge, deepfake,
  mecheche, Čerokí and kure, sárí, tvore.

Still open: the 64 `unexplained` rows (P forms without the change; 39
primary after the audit fixes, 41 before) are all cs.wiktionary table
forms with no corpus tokens (epistemologe, metodickí, psychiatre); the
TP stage counts them as dictionary-level exceptions and the README says
so.

## Fixes after the adversarial audit (2026-10-23)

Stage versions: extract 2, lexicon 3, attest 4, features 7 (the
analysis stages are listed in `docs/czech/README.md`). Tests:
`tests/czech/test_audit_fixes.py` (units) and
`tests/czech/test_audit_full_data.py` (real data).

- **Stems** (`ortho.noun_stem`). The stem is the first GEN.SG form
  minus its ending; n/t/d before the ending -ě is ň/ť/ď (ohně → ohň-,
  koně → koň-), and a lemma in ň/ť/ď gives its palatal back (loď:
  lodi → loď-). The folded common prefix of the stem cells is a
  consistency check: `stem_method = genitive_check` when it disagrees
  (39 in-domain nouns: míra měrou/mírou, kabátník with a bad INS.PL
  kabátníci, nepřítel; for manual review). Without GEN.SG the common
  prefix decides (`lcp`, 1,648), one cell form is stripped of its
  cell's ending (`strip`, 51: the pluralia tantum fleky, okuláry,
  počty, whose only stem cell is DAT.PL -ům), else the lemma (1,549).
  cs.wiktionary cells that hold several variants in one string
  ("pestřce / pestřece") are split. GEN.SG -é of an -á noun is the
  adjectival declension (vysoká, pokojská, dovolená). 506 stems of
  existing in-domain nouns changed and 130 changed target class: 58
  vowel → soft (paměť, oběť, odpověď), 12 labial → velar and 10
  labial → r (fleky, šusťáky, okuláry, manévry), 10 labial →
  coronal (počty), 7 coronal → velar (kabátník, noviník), velar →
  soft (oheň, kachně), labial → soft (peň, stupeň), labial → r (míra),
  vowel → velar (máslovník).
- **R26 suffix targets** (`lexpipe.langs.czech.morphology`).
  `suffix_evidence` now also takes a derivational form-of sense naming
  a lexeme (`formof`), another lexeme's headword line f=/dim=/adj=
  (`headword`), the -ský/-cký adjectives (`adjective`) and a base
  lexeme the word decomposes against when the word is neither
  inherited nor borrowed whole (`base`: sedlák ← sedlo, ručník ←
  ruka; not banka, matka). In-domain velar targets inside a suffix:
  nouns 1,807 → 4,702 (of 8,342), adjectives 515 → 2,698 (of 2,892;
  every -ský/-cký adjective). A seeded sample of 60 `base` nouns read
  by hand: the k is a suffix's in nearly all (jebačka, fenečka,
  hrošík, smeťák, kapelník); doubtful ones are soudka, ficka, pepka,
  and diminutives of given names that are no lexical nouns (Lenička).
- **Names** (R25): see Decisions. 13 English nouns left the domain
  (the nine zodiac constellations, Franky, Jesenicko, Pha, Prigl) and
  78 cs-only entries (73 of them in the domain before) are now names,
  which are never cs-only lexemes.
- **Ambiguity owners** (`attestation.counts_as_owner`): spelling
  variants (`alt_of`) and cs.wiktionary's form-of rows of multiword
  phrases no longer make a form ambiguous: `form_ambiguous` primary
  rows 50,182 → 44,525; 4,794 in-domain noun lemmas are
  `lemma_ambiguous`.
- **Derivational form-of entries** (see Decisions): 1,159 English nouns
  and 74 adjectives joined the domain; 367 cs-only lexemes (Pražanka,
  Angolanka) became second-source rows of the new English lexemes.
- **NDEB tokens** count root-internal sequences only: the lemma's
  final vowel and the transcription's final vowel are endings (mladý,
  zlatý, vysoký have none), and y/ý count after a velar (gymnázium:
  gy).
