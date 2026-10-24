"""Language-agnostic pipeline layer shared by the Italian and Czech plugins.

The multilingual extension (``docs/multilingual_design.md``) adds
language plugins under :mod:`lexpipe.langs` that all count, store and
verify things the same way. That shared machinery lives here, so a
cross-linguistic comparison runs literally the same code for every
language:

  - :mod:`lexpipe.spec` — the language registry (one row per ISO
    639-1 code; every path derives from it; Romanian is a
    comparison-only row), the :class:`~lexpipe.spec.LanguagePlugin`
    contract (with optional sub-stages) and the lazy plugin loader
    used by ``scripts/lang_pipeline.py``.
  - :mod:`lexpipe.context` — the :class:`~lexpipe.context.RunContext`
    every stage receives (paths, seed, ``--force``, ``--param``),
    guarded atomic writers, the licence-gated tracked text writer and
    the stage manifests that let a stage be skipped.
  - :mod:`lexpipe.sources` — pinned external resources: URL, size,
    licence, sha256 lock file (with the ETag of rolling URLs),
    resumable download.
  - :mod:`lexpipe.licensing` — which sources may be redistributed in
    tracked files.
  - :mod:`lexpipe.kaikki` — stream and filter wiktextract (kaikki.org)
    JSONL; typed entry accessors (forms, plurals, categories, IPA,
    etymology templates).
  - :mod:`lexpipe.ids` — fingerprinted lexeme_ids and root_ids.
  - :mod:`lexpipe.normalize` — orthography and IPA normalization that
    keeps accents, length, stress and syllable dots.
  - :mod:`lexpipe.frequency` — Leipzig word-frequency tables stored
    under ``data/<language>/leipzig/``; lemma and paradigm-sum lookups,
    case-sensitive counts and token totals.
  - :mod:`lexpipe.attest` — corpus attestation of forms, counts and
    generated paradigms over the thresholds τ (Italian elided clitics
    split).
  - :mod:`lexpipe.tp` — Tolerance Principle rows (θ = N/ln N) with the
    semantics of ``analysis/analyze_romanian_palatalization.R``,
    top-N downsampling, the harmonized long rows of R16–R19 (also from
    precomputed root statuses) and their tracked columns.
  - :mod:`lexpipe.learners` — simulated learners (several rules on the
    same draws, optionally in worker processes) and fragility of TP
    verdicts.
  - :mod:`lexpipe.stats` — odds ratio, Fisher exact and the seeded
    permutation null of ``scripts/build_contingency_table.py``;
    Mantel–Haenszel, stratified permutation, clustered logistic model.
  - :mod:`lexpipe.csvio` — output conventions (``True``/``False``,
    pipe-joined lists with escaped items, empty for missing), writers
    guarded by the run context and the licence gate for tracked files.
  - :mod:`lexpipe.curated` — hand decisions keyed by lexeme_id (or a
    composite key such as lexeme_id × cell) with guard fields.

Romanian keeps its own pipeline (``scripts/romanian_*.py``,
``src/romanian_processor_lib.py``, ``scripts/attach_frequencies.py``).
Nothing in that pipeline imports this package, nothing here imports
it (``tests/test_repo_guards.py``), and ``tests/test_romanian_replay.py``
pins its outputs byte for byte.

Submodules are not imported here: :mod:`lexpipe.stats` pulls in scipy,
which the CLI's source-management stages do not need.
"""

from __future__ import annotations

from pathlib import Path
from typing import Final

# src/lexpipe/__init__.py → repository root. The pipeline runs from a
# checkout (scripts put src/ on sys.path), never from an installed
# wheel, so the data directories are resolved relative to it.
PROJECT_ROOT: Final[Path] = Path(__file__).resolve().parents[2]
DATA_DIR: Final[Path] = PROJECT_ROOT / "data"
SOURCES_DIR: Final[Path] = DATA_DIR / "sources"
ANALYSIS_DIR: Final[Path] = PROJECT_ROOT / "analysis"
DOCS_DIR: Final[Path] = PROJECT_ROOT / "docs"
