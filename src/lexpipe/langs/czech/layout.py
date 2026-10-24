"""Every file the Czech plugin reads or writes, derived from a RunContext.

Paths follow design §2 and R3: inputs come from the pinned sources
under ``data/sources/`` (:mod:`lexpipe.sources`), generated files go to
``data/czech/`` (NC-licensed ones to ``data/czech/nc/``, R9), tracked
summaries to ``analysis/czech/`` and hand decisions are read from
``curated/czech/``. File names carry the registry's ``czech_`` prefix
(:meth:`lexpipe.spec.LanguageInfo.file_name`), and nothing here is
spelled for a specific root, so ``--root DIR`` relocates everything.

Stage outputs, by stage (the plugin's ``__init__`` lists the order)::

    extract     czech_en_entries.jsonl.gz   projected en.wiktionary lemmas
                czech_cs_entries.jsonl.gz   projected cs.wiktionary lemmas
                czech_form_index.sqlite     form -> owning lexemes
                czech_extract_stats.tsv
    lexicon     czech_lexicon_raw.csv       tracked, lean (extracted
                                            columns only, no names)
                czech_lexicon_full.csv      every lexeme and column
                                            (read by later stages)
                czech_cells_raw.csv         one row per cell form
                czech_lexicon_stats.tsv
    qc          nc/czech_morfflex.sqlite    MorfFlex projection (NC)
                nc/czech_morfflex_qc.csv    per-form agreement (NC)
                analysis/czech/czech_qc_report.txt        aggregates
                analysis/czech/czech_manual_worklist.csv  Wiktionary only
    attest      leipzig/freq/ces_news_2024_1M_freq.csv
                czech_attestation.csv, czech_lexeme_frequency.csv,
                czech_attestation_report.txt
    features    czech_cells.csv, czech_lexicon_complete.csv,
                czech_features_report.txt

The English Wiktionary subset ``data/sources/kaikki/en_wiktionary_cs.jsonl``
was split from the raw dump before the pipeline existed and is not in
the source registry. When it is absent, the extract stage splits the
registered raw dump itself into ``data/czech/kaikki/``
(:attr:`CzechPaths.split_subset`), since ``data/sources/`` is read-only
to derived stages.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Final

from lexpipe import kaikki, sources
from lexpipe.context import RunContext
from lexpipe.frequency import LeipzigLayout, language_layout

# Source-registry names of the inputs (lexpipe.sources.SOURCES).
RAW_DUMP_SOURCE: Final[str] = "kaikki_raw"
CS_EXTRACT_SOURCE: Final[str] = "kaikki_cs_extract"
MORFFLEX_SOURCE: Final[str] = "morfflex"
CNK_LEMMA_SOURCE: Final[str] = "cnk_syn2015_lemma"
LEIPZIG_SOURCE: Final[str] = "leipzig_ces_news_2024_1M"

# The template with the column layout of curated cell decisions.
ADJUDICATIONS_FILE: Final[str] = "cell_adjudications.csv"


def _source_path(ctx: RunContext, name: str) -> Path:
    return sources.local_path(sources.get_source(name), ctx.root)


@dataclass(frozen=True, slots=True)
class CzechPaths:
    """Resolved paths of one run (see the module docstring)."""

    # inputs
    en_subset: Path
    split_subset: Path
    raw_dump: Path
    cs_extract: Path
    morfflex: Path
    cnk_lemma: Path
    leipzig_archive: Path
    leipzig: LeipzigLayout
    adjudications: Path
    # extract
    en_entries: Path
    cs_entries: Path
    form_index: Path
    extract_stats: Path
    # lexicon
    lexicon_raw: Path
    lexicon_full: Path
    cells_raw: Path
    lexicon_stats: Path
    # qc
    morfflex_db: Path
    morfflex_lemmas: Path
    qc_csv: Path
    qc_report: Path
    worklist: Path
    # attest
    attestation: Path
    lexeme_frequency: Path
    attestation_report: Path
    # features
    cells: Path
    lexicon_complete: Path
    features_report: Path

    def en_input(self) -> Path:
        """The English Wiktionary subset to read (see module docstring)."""
        return (
            self.en_subset if self.en_subset.is_file() else self.split_subset
        )


def paths_for(ctx: RunContext) -> CzechPaths:
    """All paths of the Czech plugin under ``ctx.root``."""
    info = ctx.language
    data, nc = ctx.data_dir, ctx.nc_dir
    analysis = ctx.analysis_dir

    def name(stem: str, suffix: str = ".csv") -> str:
        return info.file_name(stem, suffix)

    subset_name = kaikki.SPLIT_FILE_TEMPLATE.format(code=info.kaikki_code)
    return CzechPaths(
        en_subset=ctx.sources_dir / "kaikki" / subset_name,
        split_subset=data / "kaikki" / subset_name,
        raw_dump=_source_path(ctx, RAW_DUMP_SOURCE),
        cs_extract=_source_path(ctx, CS_EXTRACT_SOURCE),
        morfflex=_source_path(ctx, MORFFLEX_SOURCE),
        cnk_lemma=_source_path(ctx, CNK_LEMMA_SOURCE),
        leipzig_archive=_source_path(ctx, LEIPZIG_SOURCE),
        leipzig=language_layout(info.name, ctx.root),
        adjudications=ctx.curated_dir / ADJUDICATIONS_FILE,
        en_entries=data / name("en_entries", ".jsonl.gz"),
        cs_entries=data / name("cs_entries", ".jsonl.gz"),
        form_index=data / name("form_index", ".sqlite"),
        extract_stats=data / name("extract_stats", ".tsv"),
        lexicon_raw=data / name("lexicon_raw"),
        lexicon_full=data / name("lexicon_full"),
        cells_raw=data / name("cells_raw"),
        lexicon_stats=data / name("lexicon_stats", ".tsv"),
        morfflex_db=nc / name("morfflex", ".sqlite"),
        morfflex_lemmas=nc / name("morfflex_lemmas", ".txt"),
        qc_csv=nc / name("morfflex_qc"),
        qc_report=analysis / name("qc_report", ".txt"),
        worklist=analysis / name("manual_worklist"),
        attestation=data / name("attestation"),
        lexeme_frequency=data / name("lexeme_frequency"),
        attestation_report=data / name("attestation_report", ".txt"),
        cells=data / name("cells"),
        lexicon_complete=data / name("lexicon_complete"),
        features_report=data / name("features_report", ".txt"),
    )


def check_source_sizes(ctx: RunContext) -> bool:
    """Whether registered byte sizes are enforced (default True).

    Tests run the stages on fixture-sized stand-ins of the sources and
    pass ``params={"check_source_sizes": False}``.
    """
    return bool(ctx.params.get("check_source_sizes", True))
