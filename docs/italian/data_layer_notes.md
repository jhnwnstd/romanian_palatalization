# Italian data layer: notes

Notes on `lexpipe.langs.italian` (stages `extract`, `lexicon`, `qc`,
`attest`, `features`) for the people and agents who build on it. Each
module's docstring is the authoritative description of its inputs,
outputs, columns and rules; this file records the decisions that span
modules, the data findings behind them and what is left open. Rule
numbers refer to [`../multilingual_design.md`](../multilingual_design.md).
The data dictionary of the outputs is written separately.

```bash
python scripts/lang_pipeline.py --lang it extract lexicon qc attest features
python -m pytest -q tests/italian            # needs_data tests skip without data/
```

A full run takes about 2¼ minutes (extract ≈ 63 s, features ≈ 25 s) and
at most ≈ 1.5 GB of memory (features). A stage is skipped while its
manifest matches; after a code change bump the module's `VERSION` or run
with `--force <stage>`.

## Outputs at a glance

| stage | output | tracked |
|---|---|---|
| extract | `data/italian/extract/italian_{en_lexemes,en_verbs,it_entries}.jsonl.gz`, `italian_form_index.tsv.gz` (source-stamped, `store.py`) | no |
| lexicon | `data/italian/italian_lexicon_raw.csv` (lean: extracted columns only, no proper names or multiword expressions; via `csvio.write_tracked`) | yes |
| lexicon | `data/italian/italian_lexicon_full.csv` (every lexeme, every column; what later stages read) | no |
| lexicon | `italian_plural_forms.csv` (R12), `italian_plural_drift.csv`, `italian_lexicon_summary.tsv` | no |
| qc | `data/italian/italian_qc.csv` | no |
| qc | `analysis/italian/italian_qc_report.txt`, `analysis/italian/italian_manual_worklist.csv` | yes |
| attest | `data/italian/leipzig/freq/ita_news_2024_1M_freq.csv`, `italian_attestation.csv`, `italian_attestation_summary.tsv` | no |
| features | `data/italian/italian_cells.csv`, `italian_lexicon_complete.csv` | no |

New stages are added in `src/lexpipe/langs/italian/stages.py` (one line
each; see its docstring) with their paths in `layout.py`.

## Identity (R4)

- **Lexemes** are the en entries with `pos` noun/adjective that
  `kaikki.is_lemma_entry` accepts (90,404 in the 2026-10-01 dump): the
  set `tests/lexpipe/test_ids.py` checks. `lexeme_id`s are exactly
  `ids.assign_lexeme_ids` over that set; `extract` calls it only on the
  colliding base identifiers (re-read by file offset), which gives the
  same result without holding every entry in memory
  (`tests/italian/test_full_data.py` compares against the full call).
- **Roots**: key = stem (match key minus final vowel) + the class of
  the final vowel (-o and -a together, -e and -i apart) + etymon key
  (the first source term folded, its ending kept: Latin arca and arcus
  differ) + stress position + (page etymology number when no etymon is
  known). A feminine noun that names its masculine counterpart (`m=`
  argument or a "female equivalent of X" sense) takes that lexeme's
  key, so amico N, amico ADJ and amica share one root; antico
  'ancient', antico 'anterior' (/anˈti.ko/, lbor) and àntico
  (/ˈan.ti.ko/) are three. An entry without a stress takes, for its key,
  the one stress of the other entries of its word and page etymology
  (psicotico N joins psicotico ADJ; 282 words, 27 of them dorsal, were
  split over two roots before). Before 2026-10-23 the etymon's Latin
  ending was folded (amīcus = amīca), which merged 127 groups of words
  with different etymons and no stated counterpart, most of them
  unrelated (arca/arco, casa/caso, porta/porto, blocco 'block'/blocco
  'blockade'), and the final vowel was not part of the key (dolce/dolco,
  lattice/lattico).
- **Domain**: besides multiword, hyphenated, affix, capitalized,
  nonalphabetic, abbreviation and clipping entries, the pluralia tantum
  (every gender code plural-only: brache, mezzemaniche; 191 entries) and
  the plural entries (first sense, or every sense, "plural of X":
  alimenti, eliche, crostini; 129) are out of the domain: their headword
  is not a singular.

## Plurals (R12)

- Plurals come from the head arguments (`heads.py`), never from kaikki's
  `forms`, which lose some (braccio has no plural form at all) and
  invent others (the "form" `invariable` of capotribù; `is feminine`).
- The argument grammar is larger than the design lists: besides
  `+ # - ~ ? !` there are `++` (pluralise every part), `cap*`/`cap*+`
  (capo- compounds), `#s` (lemma + s), the multiword specials
  `+first`, `+last`, `+each`, `<g:>` genders on items, and references
  nested inside labels (`<l:<<proscribed>><ref name=dop/>>`). One head
  (the noun *mangiare*) has unbalanced brackets and is flagged.
- `itcommon.py` ports `make_plural`, `make_feminine` and
  `make_masculine` from Module:it-common revision **79032535**
  (2024-04-24; content SHA-1 64837177…, fetched once on 2026-10-23 into
  `data/italian/upstream/`). Module:it-headword was not fetched; how it
  calls the port is reconstructed and validated by the drift test:
  mf/mfbysense nouns get m.pl and f.pl when the two defaults differ
  (collega → colleghi/colleghe), `+`/`cap*` on such a noun's argument 2
  likewise (capocronista → capicronisti/capocroniste), adjectives take
  f.pl from their feminine singular (editore → editrici), `fonly`
  adjectives only f.pl, `~` alone also gets the default plural.
- **Drift test** (full dump): 49,384 pure-default nouns match kaikki
  cell for cell, 0 mismatch, 1 port error (*taggiasca*, gender `?`,
  for which the module has no default), 3,960 multiword/hyphenated nouns
  not ported (Module:romance utilities). Adjectives: 24,037 match,
  1 mismatch (*grande*: kaikki renders no plural), 90 not ported. Heads
  with arguments (parser check): 11,447 match, 27 mismatch (kaikki
  losses such as braccio, and editing errors such as `2=inv`), 320 not
  ported.
- Provenance per lexeme: 81.0% module_default, 11.1% explicit, 0.3%
  mixed, 7.6% none (uncountable, invariable without rows, plurale
  tantum, no head template).

## QC (tier 1) and the worklist (tier 2)

- `plural_qc` is computed per lexeme × cell against Morph-it! noun lines
  and it.wiktionary; Morph-it! adjective lines are generated and are
  only reported (R11). `confirmed` covers two bases (`qc_basis`): an
  independent attestation, or `explicit_only` (typed by an editor,
  nothing independent either way); a module default without
  independent evidence is `module_default_only`.
- Cells: 47.8% confirmed, 50.9% module_default_only, 0.7% conflict,
  0.5% variant. Real conflicts surface (àntico's *antici* vs Morph-it!
  *antichi*, which conflates the homographs); so do it.wiktionary typos
  (*mosaicge*, *onmcologiche*). it.wiktionary's absolute superlatives
  (*timidissimi*) are tagged plural but are not plurals of the lemma and
  are left out of `it_forms` (they made 77 cells `variant`).
- `suspect = e_palatal`: 18 en.wiktionary plurals palatalize before -e
  (*fluidodinamice*, *atomistice*, …), which Italian never does; they
  are on the worklist as probable editing errors.
- The worklist lists conflicts, velar/palatal variants, suspects, and
  module defaults on dorsal stems for lemmas attested in the corpus,
  ordered by the lemma's frequency rank (rank only). Decisions go to
  `curated/italian/plural_adjudications.csv` (template; see
  `curated/italian/README.md`). Never fetch Treccani or De Mauro
  automatically.

## Features (design §5)

- `target`: dorsal (stem-final c/g, incl. sc, cc, gg and C+velar
  clusters), coronal (t/d/s/z: the "no coronal rule" rows), or none;
  `frontstem` for lexically palatal stems (-ce, -cia, -cio, -scia…);
  `nde_class` occhio (-chio/-chia/-ghio/-ghia), psiche (singular
  -che/-ghe pronounced with a velar: the IPA, or without one a velar
  plural, says so — brioche, fiche, ceviche, cloche are not; 14 lexemes,
  57 before), gimpe_like (invariable C + front vowel: specie).
- Labelled plural variants (archaic, obsolete, Old Italian, literary,
  poetic, dialectal, nonstandard, proscribed) are kept as rows with
  `filtered = True` and left out of `cell_outcome` (magi, fòcora).
- `exception_reason`: amiche-type e cells are `unexplained` (the
  non-trigger N = e rows of R16), coronals before -i too (no coronal
  rule); coronals before -e are `non exception`.
- Doublets (R19): `doublet_class` from the forms' case-sensitive corpus
  counts (`tp.classify_doublet`); dominance needs at least
  `DOUBLET_MIN_TOKENS` (10) tokens in all: sarcofagi 2 vs sarcofaghi 0,
  astrologi 3 vs astrologhi 0 are true doublets (11 of the 13 dominant
  cells were decided on 1–4 tokens).
- Attestation (`attest`): the τ flags test the folded count and need a
  token in the form's own (lowercase) spelling: forms seen only
  capitalized are names (Bocchi, Paschi, Toschi; R25). The competing
  spelling of every dorsal -i plural (`competitor`: amichi for amici,
  parci for parchi) is counted, and flagged when it is a form of another
  word (`competitor_other_lemma`: dolci, luci, merci).
- **Suffix targets (R26)**: `target_is_suffix` requires morphology: an
  af/surf template ending in an -ico-family suffix (-ico, -istico,
  -metrico, …; the last part of a `suffix`/`suf`/`confix` template is a
  suffix even without its hyphen: feroico, decaedrico; 1,669 + 57
  in-domain lexemes), or a non-inherited Latin/Greek -icus etymon on an
  adjectival lexeme (388). Inherited etymons (medico < medicus) and
  deverbal nouns (carico) are never targets; a stressed -ìco never is
  — `ico_stressed` (lexicon stage) asks whether the i of the ending
  carries the stress, so the vowel + -ico adjectives stressed on the
  vowel before it (arcaico /arˈkaj.ko/, laico, ebraico) are targets. A
  back-formation is deverbal only from a verb (metallurgo, from the noun
  metallurgia, is not). Most -ico adjectives have neither kind of
  evidence in Wiktionary and stay root targets: the suffix class is a
  lower bound.
  `suffix_base` (and `target_is_suffix_strict`) records an attested
  Italian base found by spelling (stem + o/a/e/ia); it is a noisy
  heuristic (àntico → *anta*, politico → *polito*) and is not required,
  since requiring it would drop pubblico (golden: a target).
- **Stratum (R22)**: the shared classifier of `lexpipe.etymology` (the
  same code as Czech; this plugin passes the Italian profile). The
  first transmission template decides: `inh` → inherited; `der`/`uder`
  (and `:from`) from Vulgar Latin or Old Italian → inherited, from
  another language → borrowed whole, from Italian → internal
  formation; `surf` → borrowed whole; affixation, deverbal → internal
  formation, unless the text says "Internally, from …" or opens with a
  borrowing ("From Latin …": sporcizia) → borrowed whole, or has "By
  surface analysis" while the affixation names no Italian base (a
  leading affixation of Italian bases stays internal: antichità,
  altoatesino). Without a template the text decides ("From acciacco."
  → internal formation). Etymology trees are skipped. Before the
  migration (2026-10-23) every af-first entry with "By surface
  analysis" was borrowed whole except antichità; see
  `data/research_v2/migration_diff.md` for the shifts.
- **-logo/-fago/-urgo (R27)**: masculine endings qualify by spelling
  (unless deverbal); feminine ones need a template, a Greek/Latin
  etymon in -logos/-phagos/-ourgos or a qualifying masculine
  counterpart (keeps out purga, riloga, egloga): that is the Italian
  evidence (`suffixoid`, `suffixoid_basis`) the shared detector of
  `lexpipe.suffixoids` takes, so `learned_suffixoid` is False only by a
  curated decision (spurgo). `semantic_class` is the shared
  person/thing classifier on the first gloss (the same for Czech and
  Romanian; adjectives are `adjective`); curated decisions are in
  `curated/italian/suffixoid_classes.csv` (formerly `person_thing.csv`).
- `deverbal_from_are` (etymology: it-deverbal, `:deverbal`, `:from` an
  -are verb, "Deverbal from", "(short) past participle of") is kept
  apart from `shares_root_with_1st_conj_verb` (an -are/-arsi verb
  lemma on the same stem exists: medico/medicare).
- Stress: from the primary IPA (first transcription without a regional,
  traditional or nonstandard tag), else en's accented hyphenation,
  else it.wiktionary's (`stress_source`); `velar_post_tonic` =
  the final syllable (holding the stem-final velar) follows the
  stressed one. A hyphenation is split at hyphens and spaces as well as
  hyphenation points (post - so - vi - è - ti - co) and is not used when
  it spells another word (bi‧o‧ce‧nò‧si on biocenotico; another final
  vowel is allowed) or accents the i of an -ico formed with the suffix
  (a template, or an adjective with an Italian base: bo‧so‧nì‧co);
  `stress_rejected` records the reason (27 + 8 lexemes). Where en's
  IPA and it.wiktionary's hyphenation both give a stress they agree for
  22,950 of 23,262 lexemes (312 disagree; `stress_vs_it:<source>` and
  `stress_rejected:<reason>` in `italian_lexicon_summary.tsv`).

## List order in pipe-joined columns (R5)

Source order: plural forms (`pl_forms`, `mpl_forms`, `fpl_forms`,
`forms_i`, `forms_e`, `wiktionary_forms`), `genders`, `qualifiers`,
`feminine_forms`, `masculine_forms`, `it_forms` (entries in projection
order, then forms in entry order), `provenance` (first occurrence).
Sorted: `labels`, `stress_variants`, `morphit_forms`,
`morphit_adj_forms`, `ambiguous_with`, `opportunities`. Fixed order:
`domain_flags` (multiword, hyphenated, affix, capitalized,
nonalphabetic, abbreviation, clipping), `source`.

## Open points

- No sentence-initial correction of case-sensitive counts (needs the
  Leipzig sentences file); `form_ambiguous` is lexical (form index),
  not a POS-tagged corpus count (R11's future improvement).
- Multiword and hyphenated plurals are not ported (out of domain).
- Verb paradigms are available only as form-index rows (source
  `table`, with tags) and the `italian_en_verbs` projection; the verb
  stage builds its own cells from them (`verbs.py`, `italian_verbs.csv`).
  Since extract version 2 the verb projection also carries the
  infinitive's IPA (`sounds`), which decides the root-stressed -ere
  class.
- Stress is missing for 37% of in-domain lexemes (no IPA or accented
  hyphenation in either edition); DiPI/De Mauro by hand for top-N items
  is the documented remedy.
- The tracked raw lexicon was large (≈ 26 MB, 90,404 rows × 49 columns);
  since 2026-10-23 it is lean (`lexicon.RAW_FIELDS`: lexeme_id, word,
  pos, gender, head_args as written, ipa_raw, a summary of the deciding
  etymology template, first_gloss, source; no proper names or
  multiword expressions) and the full table is the gitignored
  `italian_lexicon_full.csv`.
