# Data licence: DRAFT for the authors

> **Draft, not yet adopted.** This file was written while building the
> multilingual pipeline (design rule R9 in
> [`docs/multilingual_design.md`](docs/multilingual_design.md)). The
> authors must review it, and in particular check the Leipzig terms,
> before relying on it or publishing data. It is not legal advice.

## Code

All code in this repository is under the MIT License ([`LICENSE`](LICENSE)).
The MIT License covers the code only, not the data below.

## Data tracked in this repository

| Path | Derived from | Licence of the data |
|---|---|---|
| `data/romanian_lexicon_raw.csv` | English and Romanian Wiktionary pages, harvested by `scripts/romanian_harvester.py` (the `source` column holds each page's URL) | CC BY-SA 4.0 + GFDL (Wiktionary's dual licence) |
| `data/italian/italian_lexicon_raw.csv`, `data/czech/czech_lexicon_raw.csv` | English Wiktionary (and, for Czech-only lexemes, Czech Wiktionary), via the wiktextract extraction distributed by kaikki.org. Lean raw tables: extracted columns only (identifier, lemma, part of speech, gender/animacy, first gloss, raw IPA, head-template arguments, Czech cell forms, an etymology-template summary, source); out-of-domain proper names and multiword expressions are left out. The full stage tables (`*_lexicon_full.csv`) stay gitignored. | CC BY-SA 4.0 + GFDL |
| `analysis/italian/`, `analysis/czech/`, `analysis/comparison/` | Only sources that may be redistributed (table below). Every CSV row records `source` and may record `licence_class`; `lexpipe.csvio.write_tracked` refuses anything else. Text reports go through `RunContext.write_tracked_text`, which names their sources (see the draft policy below for the four Czech reports with NC aggregates). Leipzig counts appear only as attestation thresholds behind aggregates, never as counts. | CC BY-SA 4.0 when Wiktionary or Morph-it! data is included |
| `curated/<language>/` | Hand decisions (plural and cell adjudications, -log/-fag/-urg classes). Treccani, De Mauro and IJP are consulted by hand and recorded as coded values with a citation, never as copied entry text. | CC BY-SA 4.0 (proposed) |
| `analysis/*.csv`, `analysis/*.txt` (Romanian) | Aggregates and rows computed from the Romanian lexicon | CC BY-SA 4.0 (Wiktionary-derived) |
| `docs/sources_lock.tsv` | Sizes, SHA-256 hashes and HTTP headers of the downloaded files | No third-party content |

Attribution for Wiktionary-derived data: "Wiktionary contributors", with
a link to the page or to <https://www.wiktionary.org/>. Data obtained
through wiktextract should also cite Ylonen (2022), below.

Share-alike: a file that mixes CC BY-SA data with other data must be
released under CC BY-SA 4.0. CC BY-NC-SA data cannot be mixed into such
a file, because the two share-alike licences are incompatible.

## External sources

| Source | Licence | Derived data in this repository? | Notes |
|---|---|---|---|
| Wiktionary via kaikki.org/wiktextract (English, Italian and Czech editions) | CC BY-SA 4.0 + GFDL | Yes | kaikki.org regenerates its files in place, so `docs/sources_lock.tsv` records each file's SHA-256, ETag and Last-Modified. The per-language subsets split from the English dump (`kaikki_raw_it`, `kaikki_raw_cs`) are pinned the same way. |
| Romanian Wiktionary (harvested by `scripts/romanian_harvester.py`) | CC BY-SA 4.0 + GFDL | Yes | Source code `wiktionary_ro`: most noun rows of the harvested Romanian lexicon come from it. The harmonized Romanian tables in `analysis/comparison/` name `wiktionary_ro\|wiktionary_en` (aggregates) or each row's edition (the -log table). |
| Morph-it! 0.4.8 | CC BY-SA 2.0 / LGPL (dual) | Yes | Fetched from a GitHub mirror and pinned by SHA-256. The official tarball is behind a human-verification page. |
| ČNK SYN2015 comparative frequency lists | CC BY 4.0 | Yes | Cross-check for the Czech frequencies. |
| MorfFlex CZ 2.1 | CC BY-NC-SA 4.0 | **No** | Derived files go to `data/czech/nc/` (gitignored). Commercial licences are available from ÚFAL on request. |
| DeriNet 2.3 | CC BY-NC-SA 4.0 | **No** | As for MorfFlex. |
| Czech Etymological Lexicon 1.0 | CC BY-NC-SA 4.0 | **No** | As for MorfFlex. |
| Leipzig Corpora Collection (`ita_news_2024_1M`, `ces_news_2024_1M`, `ron_news_2024_1M`) | **Terms unverified** | **No**, until a person has read the download terms | Frequency tables stay in gitignored `data/` directories. The tracked Romanian outputs hold no Leipzig counts, only results computed with them (e.g. the top-1000 TP summary). |
| Treccani, De Mauro (Italian); IJP (Czech) | Consulted by hand | Coded values and citations only | Treccani opts out of text and data mining; IJP's robots.txt disallows dictionary queries, and copying its content needs ÚJČ's consent. |
| DEX (Romanian QC, dexonline.ro) | **Terms unverified** | **No**, until the authors have checked dexonline's terms | Source code `dex_qc` (class `unverified`). DEX-derived files (`data/romanian_lexicon_raw_dex.csv`, the DEX cache) are gitignored. DEX QC also filled the plural of 837 nouns that Wiktionary lacked (lexicon `source` ending in `\| dex_qc`). The harmonized per-noun -log table leaves those nouns out; the aggregates in `analysis/comparison/` count them, as the legacy Romanian outputs (R1) do, and the legacy per-row files (`analysis/palatalization_mismatches.csv`) list some of their plurals: the authors should decide this together with dexonline's terms. |

`src/lexpipe/licensing.py` holds the machine-readable form of this table.
Change it and this file together.

## Aggregate counts from non-commercial sources (DRAFT, for the authors' decision)

> **Not adopted.** The pipeline implements this policy behind an explicit
> switch, so that the Czech reports can be generated. The authors must
> accept, change or reject it before any such report is published. It is
> not legal advice: whether aggregate counts fall under the NC terms at
> all is a question this draft does not settle.

Some tracked text reports summarise analyses that rest partly on
CC BY-NC-SA sources (MorfFlex CZ 2.1, DeriNet 2.3, the Czech
Etymological Lexicon). An example is the number of base–derivative
pairs DeriNet links, or an odds ratio computed over them. CSV tables
never carry such data: `lexpipe.csvio.write_tracked` refuses NC rows
outright, and this policy does not change that. The reports written
under it are `analysis/czech/czech_qc_report.txt` (MorfFlex agreement
counts), `czech_derivations_report.txt`, `czech_spellout.txt` and
`czech_noevidence.txt` (DeriNet pairs, czetyl and DeriNet strata);
each starts with the standard header.

Proposed policy for tracked **text reports**:

1. **Aggregates only.** A report may contain counts, proportions and
   test statistics computed from NC sources. It may not contain rows,
   word forms, lemma lists, examples or any other text copied from them.
   Each plugin is responsible for this content rule.
2. **Known terms only.** Every source behind the report must have a
   known licence. Sources whose terms nobody has checked yet (Leipzig,
   DEX QC: class `unverified`) are refused even as aggregates.
3. **Standard header.** The report starts with the block produced by
   `lexpipe.licensing.aggregate_header`. The block names every source
   with its licence and states that the report holds aggregate counts
   only. When an NC source is among them, it adds that the counts are
   tracked under this draft policy and may not be used commercially.
4. **One gate.** Reports are written with
   `RunContext.write_tracked_text(path, text, sources=[...], nc_aggregates=True)`,
   which enforces rules 2 and 3 and accepts tracked locations only.
   Without `nc_aggregates`, every source must be redistributable, as
   for CSV tables.

Options for the authors:

- **Adopt** the policy as written.
- **Reject** it: remove `LicenceClass.BY_NC_SA` from
  `AGGREGATE_CLASSES` in `src/lexpipe/licensing.py`. Every report that
  names an NC source then fails to write until its NC sections move to
  `data/<language>/nc/` (gitignored).
- **Ask the licensors**: ÚFAL (Charles University) for MorfFlex and
  DeriNet, which offers commercial licences on request, and the rights
  holders of the Czech Etymological Lexicon, which is extracted from
  Rejzek's *Český etymologický slovník*. Written permission would
  settle the question for that source.

## References

- Ylonen, Tatu. 2022. Wiktextract: Wiktionary as Machine-Readable
  Structured Data. In *Proceedings of the Thirteenth Language Resources
  and Evaluation Conference (LREC 2022)*, pages 1317–1325, Marseille.
  ELRA. <https://aclanthology.org/2022.lrec-1.140/>
- Zanchetta, Eros and Marco Baroni. 2005. Morph-it! A free corpus-based
  morphological resource for the Italian language. In *Proceedings of
  Corpus Linguistics 2005*. University of Birmingham.
- Hajič, Jan, Jaroslava Hlaváčová, Marie Mikulová, Milan Straka and
  Barbora Štěpánková. 2024. *MorfFlex CZ 2.1*. LINDAT/CLARIAH-CZ
  digital library, Charles University.
  <http://hdl.handle.net/11234/1-5833>. Reference paper: Hlaváčová,
  Jaroslava, Marie Mikulová, Barbora Štěpánková, Milan Straka and Jan
  Hajič. 2026. MorfFlex: Handling Rich Morphology. In *Proceedings of
  the Fifteenth Language Resources and Evaluation Conference
  (LREC 2026)*, pages 11495–11505.
- Olbrich, Michal, Viktória Brezinová, Šárka Dohnalová, Vojtěch John,
  Lukáš Kyjánek, Aleš Papáček, Emil Svoboda, Magda Ševčíková, Jonáš
  Vidra and Zdeněk Žabokrtský. 2025. *DeriNet 2.3*. LINDAT/CLARIAH-CZ
  digital library, Charles University.
  <http://hdl.handle.net/11234/1-5846>. Format: Vidra, Jonáš, Zdeněk
  Žabokrtský, Magda Ševčíková and Lukáš Kyjánek. 2019. DeriNet 2.0:
  Towards an All-in-One Word-Formation Resource. In *Proceedings of the
  Second International Workshop on Resources and Tools for Derivational
  Morphology (DeriMo 2019)*, pages 81–89, Prague.
- Rejzek, Jiří, Aleš Papáček, Viktória Brezinová and Zdeněk Žabokrtský.
  2025. *Czech Etymological Lexicon 1.0*. LINDAT/CLARIAH-CZ digital
  library, Charles University. <http://hdl.handle.net/11234/1-5845>.
  It was extracted from Rejzek, Jiří. 2015. *Český etymologický
  slovník*. LEDA.
- Křen, M., V. Cvrček, T. Čapka, A. Čermáková, M. Hnátková, L.
  Chlumská, T. Jelínek, D. Kováříková, V. Petkevič, P. Procházka, H.
  Skoumalová, M. Škrabal, P. Truneček, P. Vondřička and A. Zasina.
  2015. *SYN2015: reprezentativní korpus psané češtiny*. Ústav Českého
  národního korpusu FF UK, Praha. The lists themselves: *Český národní
  korpus: Srovnávací frekvenční seznamy*. Ústav Českého národního
  korpusu FF UK, Praha 2016.
- Goldhahn, Dirk, Thomas Eckart and Uwe Quasthoff. 2012. Building Large
  Monolingual Dictionaries at the Leipzig Corpora Collection: From 100
  to 200 Languages. In *Proceedings of the Eighth International
  Conference on Language Resources and Evaluation (LREC'12)*, pages
  759–765, Istanbul. ELRA. <https://aclanthology.org/L12-1154/>

Author lists and years were checked against the ACL Anthology, LINDAT
and korpus.cz records on 2026-10-23. The Morph-it! entry was checked
against the University of Bologna's record.

## Open points for the authors

1. Read the Leipzig Corpora Collection's download terms. If they allow
   it, set `leipzig` to a redistributable class in
   `src/lexpipe/licensing.py` and update the table above. Do the same
   for dexonline's terms and the `dex_qc` source code.
2. Decide the draft policy on aggregate counts from non-commercial
   sources (section above).
3. Decide the licence statement for `curated/` and the tracked analysis
   outputs (CC BY-SA 4.0 is proposed above).
4. Add a copy of the CC BY-SA 4.0 legal code (for example
   `LICENSES/CC-BY-SA-4.0.txt`), and add the data attributions to
   `CITATION.cff`.
