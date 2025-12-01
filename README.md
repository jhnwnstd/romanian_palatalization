# Romanian Palatalization 
## Overview
Builds a Romanian lexical dataset with IPA transcriptions to analyze palatalization patterns. Tests whether stem-final consonants (c, g, t, d, s, z) systematically avoid mutation before front vowels in plural inflection (-i/-e) and denominal derivation. Compares dorsal vs. coronal mutation rates while controlling for suffix class and corpus frequency to evaluate Steriade (2008) avoidance hypothesis.

The pipeline:
1. Harvests lexical data from Wiktionary → `romanian_lexicon_raw.csv`
2. Downloads Leipzig corpus data → `data/leipzig/freq/*.csv`
3. Quality control with DEX (Romanian dictionary) → `romanian_lexicon_raw_dex.csv`
4. Derives palatalization features → `romanian_lexicon_complete.csv`
5. Attaches corpus frequency data → `romanian_lexicon_with_freq.csv`
6. Statistical analysis and modeling

## Quick Start

### Automated Pipeline (Recommended)
```bash
# Install dependencies
pip install -r requirements.txt

# Verify installation (optional but recommended)
python3 verify_dependencies.py

# Run pipeline
./run_pipeline.sh              # Run entire pipeline (skips stages if outputs exist)
./run_pipeline.sh --force      # Re-harvest from scratch (slow)
```

### Manual Pipeline (Step by Step)
```bash
# Install Python dependencies
pip install -r requirements.txt

# 1. Extract from Wiktionary (takes several hours)
python scripts/romanian_harvester.py

# 2. Download corpus frequency data
python scripts/download_leipzig.py

# 3. Quality control with DEX
python scripts/dex_qc_main_csv.py

# 4. Process and derive palatalization features
python scripts/romanian_processor_main.py

# 5. Attach corpus frequencies
python scripts/attach_frequencies.py

# 6. Analyze (R)
Rscript analysis/analyze_romanian_palatalization.R
```

## Repository Structure
```
romanian_palatalization/
├── run_pipeline.sh   Automated pipeline orchestration
├── scripts/          Individual pipeline scripts
├── src/              Library modules and utilities
├── analysis/         Statistical analysis (R)
└── data/             CSV files and caches (gitignored)
```

**Pipeline:**
- [run_pipeline.sh](run_pipeline.sh) - Complete automated pipeline

**Scripts:**
- [scripts/romanian_harvester.py](scripts/romanian_harvester.py) - Wiktionary extraction
- [scripts/download_leipzig.py](scripts/download_leipzig.py) - Download Leipzig corpora
- [scripts/dex_qc_main_csv.py](scripts/dex_qc_main_csv.py) - QC with DEX oracle
- [scripts/romanian_processor_main.py](scripts/romanian_processor_main.py) - Feature derivation
- [scripts/attach_frequencies.py](scripts/attach_frequencies.py) - Attach frequency data

**Modules:**
- [src/romanian_processor_lib.py](src/romanian_processor_lib.py) - Core derivation functions
- [src/wiktionary_normalizer.py](src/wiktionary_normalizer.py) - Text normalization
- [src/dex_utils.py](src/dex_utils.py) - DEX utilities

**Analysis:**
- [analysis/analyze_romanian_palatalization.R](analysis/analyze_romanian_palatalization.R) - Statistical analysis

**Documentation:**
- [docs/data_dictionary.md](docs/data_dictionary.md) - Complete field definitions and derivation logic

## Citation
Steriade, D. (2008). A pseudo-cyclic effect in Romanian morphophonology. In A. Bachrach & A. Nevins (Eds.), Inflectional identity (pp. 313–358). Oxford University Press.
