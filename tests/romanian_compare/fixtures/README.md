# Fixtures of the harmonized Romanian comparison

Small inputs for `tests/romanian_compare/`, laid out as the files the
plugin reads (`lexpipe.langs.romanian_compare`):

| file | stands for | content |
|---|---|---|
| `lexicon.csv` | `data/romanian_lexicon_with_freq.csv` | 55 rows (54 nouns, one adjective) chosen to cover every class the comparison distinguishes: c/g/t/d/s/z targets with -i, -e and -uri plurals, the st/sc/ct clusters, the structural NDE classes (genunchi, păduche, ureche, cruce, anestezie), suffix targets (critic, turist), a palatal context (aluniște), a suppletive stem (oaste), -ică suppletion (alunică), the legacy "foreign plural" (albanez) and "suppletion" (fată, pată) exclusions, -log nouns (psiholog, dialog, catalog, monolog) and a lemma listed twice (esperanto). |
| `verb_paradigm_verdict.csv` | `data/verb_paradigm_verdict.csv` | 15 noun–verb pairs with 2sg, 3sg and 3sg-subjunctive forms. |
| `freq.csv` | `data/leipzig/freq/ron_news_2024_1M_freq.csv` | `word,freq` counts. |

## Provenance and licence

- `lexicon.csv` copies rows of the Romanian lexicon, which
  `scripts/romanian_harvester.py` harvested from English and Romanian
  Wiktionary pages (the `source` column links each page). Text and
  data from Wiktionary are © Wiktionary contributors and available
  under CC BY-SA 4.0 and the GFDL
  (<https://en.wiktionary.org/>, <https://ro.wiktionary.org/>). The
  derived columns (`stem_final`, `mutation`, `opportunity`, ...) are
  the Romanian pipeline's output for those rows. These rows were *not*
  extracted with wiktextract (Ylonen 2022), unlike the Italian and
  Czech fixtures; DEX notes and `dex_qc` source tags were removed.
- **All frequencies are invented.** The `freq_ron_news_2024_1M`
  column and `freq.csv` hold synthetic counts (a hash of the form,
  0–60, about 30% zeros) plus a filler token `de` with 1,000,000
  tokens that sets the corpus size; `urși` / `urşi` test the cedilla
  fold. No count comes from the Leipzig Corpora Collection, whose terms
  are unverified (R9).
- The verb forms in `verb_paradigm_verdict.csv` were written by hand
  from the regular Romanian conjugation patterns; the verdict column
  mimics `scripts/audit_verb_alternation.py` (including its reading of
  `legendezi` as ALTERNATES). Nothing is copied from DEX.

Ylonen, Tatu. 2022. Wiktextract: Wiktionary as Machine-Readable
Structured Data. In *Proceedings of LREC 2022*, 1317–1325.
