# Romanian Palatalization 

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## Overview

This repository builds a Romanian lexical dataset with IPA transcriptions and derives features for analyzing stem-final palatalization. It tests whether stem-final consonants (`c, g, t, d, s, z`) systematically avoid mutation before front vowels in plural inflection (`-i/-e`) and denominal derivation, and compares dorsal vs. coronal mutation rates (with suffix class and corpus frequency) in the spirit of Steriade (2008).

The pipeline:

1. Harvest lexical data from Wiktionary → `romanian_lexicon_raw.csv`
2. Attach Leipzig corpus frequency data → `data/leipzig/freq/*.csv`
3. Quality control with DEX → `romanian_lexicon_raw_dex.csv`
4. Derive palatalization features → `romanian_lexicon_complete.csv`
5. Attach corpus frequencies → `romanian_lexicon_with_freq.csv`
6. Run statistical analysis (R)

See `docs/data_dict.md` for full field definitions.

## Quick Start

### 1. Install Python dependencies

```bash
pip install .
```

For development (includes linting, type checking, and testing tools):

```bash
pip install -e .[dev]
```

For the auxiliary analysis and migration scripts as well:

```bash
pip install -e .[dev,tools]
```

Optional sanity check:

```bash
python3 verify_dependencies.py
```

### 2. Run the pipeline

Recommended:

```bash
./run_pipeline.sh          # skips data collection stages if data exists
./run_pipeline.sh --force  # re-run all stages 
```

The pipeline resumes from existing data:
- If raw harvest data exists → skips harvest (saves hours/days)
- If DEX QC data exists → skips DEX QC
- Always runs processing stages (stages 4-6)

**Safety note**: `--force` re-runs stages by overwriting outputs, it does not delete existing data first. 

Manual step-by-step (equivalent to the shell script):

```bash
# 1. Harvest from Wiktionary (slow)
python scripts/romanian_harvester.py

# 2. Download corpus frequency data
python scripts/download_leipzig.py

# 3. QC with DEX
python scripts/dex_qc_main_csv.py

# 4. Derive palatalization features
python scripts/romanian_processor_main.py

# 5. Attach corpus frequencies
python scripts/attach_frequencies.py

# 6. Run analysis (R)
Rscript analysis/analyze_romanian_palatalization.R
```

## Multilingual pipelines

A second pipeline (`src/lexpipe/`, design in
[`docs/multilingual_design.md`](docs/multilingual_design.md)) builds
Italian and Czech lexicons from kaikki.org's wiktextract dumps and
recounts Romanian (read-only) with the same code. The three languages'
Tolerance Principle tables, learner simulations, co-variation
statistics and -log tables therefore come from one implementation, and
a final step joins them into one cross-language report.

**Start with the results:**
[`analysis/comparison/README.md`](analysis/comparison/README.md) holds
the cross-language report. It has the comparison tables and a check of
the draft's learnability paragraph (Inkelas criterion, plural-less
nouns, the Tolerance Principle, allomorph selection, Krämer's
speakers) against all three languages.

### Running

```bash
python scripts/lang_pipeline.py --lang it            # Italian: all stages
python scripts/lang_pipeline.py --lang cs            # Czech: all stages
python scripts/lang_pipeline.py --lang ro            # Romanian comparison only
PYTHONPATH=src python -m lexpipe.compare             # cross-language summary
PYTHONPATH=src python -m lexpipe.compare --check     # exit 1 if it is stale
python scripts/lang_pipeline.py --lang cs --list     # list a language's stages
python scripts/lang_pipeline.py --lang cs --force tp --param workers=4 tp
```

The Romanian pipeline above is untouched: the comparison only reads its
outputs and writes `analysis/comparison/` and `data/romanian_compare/`.

- **Stages.** Each stage is skipped while its manifest (inputs,
  parameters, outputs) matches; `--force STAGE` re-runs one. Run the
  three languages first and `lexpipe.compare` last: it reads only their
  tracked outputs and takes about a second.
- **Sources.** The Italian and Czech data layers need the pinned source
  files under `data/sources/` (`fetch`, `verify-sources`; see
  [`docs/runbooks.md`](docs/runbooks.md)). The analysis stages take a
  few minutes per language.
- **Statistics.** The R23 logistic models need the optional statistics
  extra (`pip install -e '.[stats]'`, statsmodels); without it the
  reports say so.

### Outputs and where the results are

| what | where |
|---|---|
| cross-language report and tables | `analysis/comparison/README.md`, `analysis/comparison/crosslang_summary.csv` |
| Romanian harmonized recount | `analysis/comparison/romanian_*` (section at the end of the comparison README) |
| Italian results | `analysis/italian/` (TP long table, learners, summaries, reports); method and findings in `docs/italian/README.md` |
| Czech results | `analysis/czech/` (same layout); `docs/czech/README.md` |
| Logical Phonology analyses | `src/phonology/analyses/{italian,czech}_palatalization.py`; `docs/<language>/phonology.md`; validator reports `analysis/<language>/*_rule_validation*.txt` |
| column definitions | `docs/italian/data_dict.md`, `docs/czech/data_dict.md`, the docstring of `src/lexpipe/compare.py` |
| tracked raw lexicons (lean) | `data/italian/italian_lexicon_raw.csv`, `data/czech/czech_lexicon_raw.csv` |
| hand decisions | `curated/<language>/` |

Tracked files hold aggregates and Wiktionary-derived columns only. The
full stage tables, per-word corpus counts and everything derived from
CC BY-NC-SA sources stay in the gitignored `data/<language>/`
(`data/<language>/nc/` for the NC ones).

Shared modules every language uses:

- `lexpipe.tp`: TP rows;
- `lexpipe.learners`: simulated learners, every cell on the same draws;
- `lexpipe.attest`: corpus thresholds τ;
- `lexpipe.stats`: odds ratios, Mantel–Haenszel, permutations, the R23
  model;
- `lexpipe.etymology`: Wiktionary etymology templates and the R22
  stratum;
- `lexpipe.suffixoids`: the R27 -log/-fag/-urg classifiers;
- `lexpipe.compare`: the cross-language join.

## Repository Structure

```text
romanian_palatalization/
├── run_pipeline.sh              # pipeline orchestrator
├── pyproject.toml               # package metadata and dependencies
├── scripts/                     # individual pipeline stages
│   ├── romanian_harvester.py
│   ├── download_leipzig.py
│   ├── dex_qc_main_csv.py
│   ├── romanian_processor_main.py
│   ├── attach_frequencies.py
│   └── lang_pipeline.py         # Italian / Czech / Romanian comparison
├── src/                         # library modules
│   ├── romanian_processor_lib.py
│   ├── wiktionary_normalizer.py
│   ├── dex_utils.py
│   └── lexpipe/                 # shared multilingual core + plugins
│       ├── compare.py           # cross-language summary (python -m lexpipe.compare)
│       └── langs/{italian,czech,romanian_compare}/
├── analysis/                    # R scripts for statistical analysis
│   ├── analyze_romanian_palatalization.R
│   └── {italian,czech,comparison}/   # tracked multilingual aggregates
├── curated/                     # hand decisions (per language)
├── docs/
│   ├── data_dict.md             # field definitions & derivation logic
│   ├── multilingual_design.md   # design of the multilingual extension
│   └── {italian,czech}/         # data dictionaries and notes
└── data/                        # CSVs and caches (gitignored, except
                                 #   the raw lexicons)
```

## Requirements

* Python 3.11 or higher
* R 4.x (for statistical analysis)
* R packages used in `analysis/analyze_romanian_palatalization.R`
  (including `cmdstanr` and other Bayesian modeling dependencies)

Python dependencies are defined in [`pyproject.toml`](pyproject.toml).

## License

This project is licensed under the MIT License. See [`LICENSE`](LICENSE) for details.

## Data sources and licences

The MIT License covers the code only. The data come from the sources
below. [`DATA_LICENSE.md`](DATA_LICENSE.md) (a draft for the authors)
gives the details, and `src/lexpipe/licensing.py` keeps data that may
not be redistributed out of tracked files.

| Source | Used for | Licence | In tracked files |
|---|---|---|---|
| Wiktionary (en, ro, it, cs) | lexicons, IPA, etymology | CC BY-SA 4.0 + GFDL | yes |
| Morph-it! 0.4.8 | Italian QC | CC BY-SA 2.0 / LGPL | yes |
| ČNK SYN2015 frequency lists | Czech frequency cross-check | CC BY 4.0 | yes |
| MorfFlex CZ 2.1, DeriNet 2.3, Czech Etymological Lexicon 1.0 | Czech QC, derivations, etymology | CC BY-NC-SA 4.0 | no (`data/czech/nc/`) |
| Leipzig Corpora Collection | corpus frequencies | terms not yet verified | no |

If you use the data, please also cite:

- Ylonen (2022), *Wiktextract: Wiktionary as Machine-Readable Structured
  Data*, LREC 2022: the Italian and Czech lexicons (via kaikki.org).
- Zanchetta & Baroni (2005), *Morph-it! A free corpus-based
  morphological resource for the Italian language*, Corpus Linguistics
  2005.
- Hajič, Hlaváčová, Mikulová, Straka & Štěpánková (2024), *MorfFlex CZ
  2.1*, LINDAT/CLARIAH-CZ (reference paper: Hlaváčová et al. 2026,
  LREC 2026).
- Olbrich et al. (2025), *DeriNet 2.3*, LINDAT/CLARIAH-CZ (format paper:
  Vidra et al. 2019, DeriMo 2019).
- Rejzek, Papáček, Brezinová & Žabokrtský (2025), *Czech Etymological
  Lexicon 1.0*, LINDAT/CLARIAH-CZ.
- Křen et al. (2015), *SYN2015: reprezentativní korpus psané češtiny*,
  and *Srovnávací frekvenční seznamy* (ÚČNK, 2016).
- Goldhahn, Eckart & Quasthoff (2012), *Building Large Monolingual
  Dictionaries at the Leipzig Corpora Collection: From 100 to 200
  Languages*, LREC 2012.

## Citation

If you use this code or dataset in your research, please cite:

```bibtex
@software{romanian_palatalization,
  title  = {Romanian Palatalization Analysis},
  author = {John Winstead},
  year   = {2025},
  url    = {https://github.com/jhnwnstd/romanian_palatalization}
}
```

Primary theoretical reference:

> Steriade, D. (2008). A pseudo-cyclic effect in Romanian morphophonology.
> In A. Bachrach & A. Nevins (Eds.), *Inflectional Identity* (pp. 313–360). Oxford University Press.

## Acknowledgments

* Kyle Gorman for guidance on the analysis and data design
* Wiktionary contributors for lexical data
* Leipzig Corpora Collection for frequency data
* DEX (Dicționarul explicativ al limbii române) for quality control
