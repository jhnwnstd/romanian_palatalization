# Italian test fixtures

Small excerpts used by `tests/italian/` to run the Italian data layer
(`lexpipe.langs.italian`: extract, lexicon, qc, attest, features) on a
fixture-sized tree. Lines are copied byte for byte from the pinned
sources, so the pipeline reads exactly what it reads from the full
files.

## Files and licences

| file | content | source and licence |
|---|---|---|
| `en_wiktionary_it_sample.jsonl` | 62 raw wiktextract records of the English-Wiktionary Italian split (amico, amica, antico ×4, fuoco, mago, chirurgo, greco ×3, belga, valigia, luce, medico, parco, lago, teologo, carico, occhio, pubblico, elettrico, psiche ×2, specie, catalogo, stomaco, form-of pages carichi/medici/stomaci/stomachi, the verbs caricare, medicare, parcare, amicare, amicarsi, pubblicare, the name Medici) | English Wiktionary contributors, CC BY-SA 4.0 and GFDL; extracted with wiktextract and distributed by kaikki.org (raw dump of 2026-10-21, enwiktionary dump of 2026-10-01) |
| `it_extract_sample.jsonl` | 37 records of the Italian-Wiktionary (Wikizionario) extract for the same words | Wikizionario contributors, CC BY-SA 4.0 and GFDL; extracted with wiktextract and distributed by kaikki.org |
| `morphit_sample.txt` | 86 noun and positive-adjective lines of Morph-it! 0.4.8 for the same lemmas | Eros Zanchetta and Marco Baroni, Morph-it!, dual-licensed CC BY-SA 2.0 / LGPL |
| `en_wiktionary_it_analysis_sample.jsonl` | 58 further raw wiktextract records for the analysis stages (derivations … docs): bases unico, tipico, lega ×2, musica, musico ×2, politica, politico ×2, porco ×2, bianco ×3, cieco ×2, opaco, traffico; derivatives grecista, grecismo, grecizzare, grecheggiare, fochista, antichità, antichizzare, unicità, tipicità, tipicizzare (alt-of), leghista, legista, musicista, politichese, amichevole, amichetto, medichessa, porcino ×2, porcello, porcile, bianchezza, bianchire, biancheggiare, cecità, opacità, teologia, astrologo, dialogo; verbs vincere, cercare, tacere, cuocere, cucire, fuggire, dire, trafficare, legare ×2 | English Wiktionary contributors, CC BY-SA 4.0 and GFDL; extracted with wiktextract and distributed by kaikki.org (same dump as above) |

The Leipzig corpus is **not** excerpted: its terms have not been checked
by hand (R9). The tests build a tiny words table with made-up counts in
the Leipzig file format instead (`tests/italian/helpers.py`; the
analysis tests add the made-up counts of `tests/italian/analysis_helpers.py`).

## Citations

- Tatu Ylonen. 2022. Wiktextract: Wiktionary as Machine-Readable
  Structured Data. In *Proceedings of the Thirteenth Language Resources
  and Evaluation Conference (LREC 2022)*, pages 1317–1325, Marseille.
  ELRA.
- Wiktionary, the free dictionary: https://en.wiktionary.org/ and
  https://it.wiktionary.org/ (page histories list the contributors of
  each entry).
- Eros Zanchetta and Marco Baroni. 2005. Morph-it! A free corpus-based
  morphological resource for the Italian language. In *Proceedings of
  Corpus Linguistics 2005*, University of Birmingham.

## Regenerating

The excerpts were cut from `data/sources/kaikki/en_wiktionary_it.jsonl`,
`data/sources/kaikki/it-extract.jsonl.gz` and
`data/sources/morphit/morph-it_048_utf8.txt` by selecting whole lines on
the `word` field (en, it) or the lemma column (Morph-it!); the analysis
excerpt selects whole lines on (`word`, `pos`) and leaves out the lines
already in `en_wiktionary_it_sample.jsonl`. Re-cut them only together
with the golden expectations in `tests/italian/`.
