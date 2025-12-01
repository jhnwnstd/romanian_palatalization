# Romanian Palatalization 

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
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

See `docs/data_dictionary.md` for full field definitions.

## Quick Start

### 1. Install dependencies

```bash
pip install -r requirements.txt
```

Optional sanity check:

```bash
python3 verify_dependencies.py
```

### 2. Run the pipeline

Recommended:

```bash
./run_pipeline.sh          # runs all stages, skips if harvest or DEX outputs exist
./run_pipeline.sh --force  # re-run everything from scratch (slow)
```

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

## Repository Structure

```text
romanian_palatalization/
├── run_pipeline.sh              # pipeline orchestrator
├── requirements.txt             # Python dependencies
├── scripts/                     # individual pipeline stages
│   ├── romanian_harvester.py
│   ├── download_leipzig.py
│   ├── dex_qc_main_csv.py
│   ├── romanian_processor_main.py
│   └── attach_frequencies.py
├── src/                         # library modules
│   ├── romanian_processor_lib.py
│   ├── wiktionary_normalizer.py
│   └── dex_utils.py
├── analysis/                    # R scripts for statistical analysis
│   └── analyze_romanian_palatalization.R
├── docs/
│   └── data_dict.md             # field definitions & derivation logic
└── data/                        # CSVs and caches (gitignored)
```

## Requirements

* Python 3.9 or higher
* R 4.x (for statistical analysis)
* R packages used in `analysis/analyze_romanian_palatalization.R`
  (including `cmdstanr` and other Bayesian modeling dependencies)

Python dependencies are listed in [`requirements.txt`](requirements.txt).

## License

This project is licensed under the MIT License. See [`LICENSE`](LICENSE) for details.

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