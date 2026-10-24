# Czech test fixtures

Small excerpts of Wiktionary data, used by `tests/czech/`.

| File | Content | Source |
|---|---|---|
| `en_wiktionary_cs_sample.jsonl` | 33 English Wiktionary entries of the Czech section (`lang_code == "cs"`) | kaikki.org raw wiktextract dump `raw-wiktextract-data.jsonl.gz`, extracted 2026-10-21 from the enwiktionary dump of 2026-10-01 |
| `cs_extract_sample.jsonl` | 32 Czech Wiktionary (Wikislovník) entries with `lang_code == "cs"` | kaikki.org `cs-extract.jsonl.gz`, extracted 2026-10-21 from the cswiktionary dump of 2026-10-01 |
| `en_wiktionary_cs_derivations.jsonl` | 23 English Wiktionary entries for the analysis stages (derivatives and their bases) | as `en_wiktionary_cs_sample.jsonl` |
| `cs_extract_derivations.jsonl` | 1 Czech Wiktionary entry (ligový) | as `cs_extract_sample.jsonl` |

Lines are copied **verbatim** (byte for byte) from the pinned files
`data/sources/kaikki/en_wiktionary_cs.jsonl` and
`data/sources/kaikki/cs-extract.jsonl.gz`; nothing was edited.

Entries: kluk, žák, vrah, Čech (common noun and surname), biolog, jazyk,
park, šach (inanimate 'check' and animate 'shah'), plech, ruka,
diplomat, student, kniha, noha (×2), liga, potok, rok (the words the
design's golden expectations name, `docs/multilingual_design.md` §8),
plus entries exercising specific code paths: český and šachový
(adjectives; šachové is also a form of šach 'shah'), otec (VOC otče),
bratr (VOC bratře), agent (animate and inanimate tables), Inka (an
ethnonym and two given names filed under a Noun header), kuře (t-stem,
out of domain) and the form-of entries kluci, Češi, šachy, roce.

The two `*_derivations.jsonl` files serve the analysis stages
(`tests/czech/analysis_tree.py`): bases banka, droga, dráha, skála
(stems with a long/short root vowel), strach, tank, vzduch and
derivatives bankovní, drogový, klouček, knižní, knížka, parkovat,
ručit, ruční, strašný, tankista, vzduchový, vzdušný, šachista, žačka,
žákyně, bioložka (en) and ligový (cs). They are appended to the sample
files when the analysis tree is built.

## Licence and attribution

The text of these entries comes from Wiktionary and is available under
the Creative Commons Attribution-ShareAlike 4.0 licence (CC BY-SA 4.0)
and the GNU Free Documentation License (GFDL). Authors: the
contributors to English Wiktionary (<https://en.wiktionary.org/>) and
Czech Wiktionary (<https://cs.wiktionary.org/>); the edit history of
each page lists them.

The entries were extracted with wiktextract and distributed by
kaikki.org:

> Tatu Ylonen. 2022. Wiktextract: Wiktionary as Machine-Readable
> Structured Data. In Proceedings of the Thirteenth Language Resources
> and Evaluation Conference (LREC 2022), pages 1317–1325, Marseille.
> ELRA.

The other inputs the Czech tests need (MorfFlex rows, a Leipzig words
table, a CNK lemma list, DeriNet trees and czetyl rows) are
**synthetic**: small made-up tables written by the tests themselves.
No MorfFlex, DeriNet or czetyl (CC BY-NC-SA), Leipzig or CNK data is
stored here.
