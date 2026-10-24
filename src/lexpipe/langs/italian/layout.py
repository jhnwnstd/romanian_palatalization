"""Every file the Italian stages read or write, derived from the context.

Paths come from the :class:`lexpipe.context.RunContext` (root, registry
row) and the source registry (:mod:`lexpipe.sources`), never from
literals elsewhere (R3), so ``--root DIR`` moves the whole pipeline::

    sources (read-only)
      data/sources/kaikki/en_wiktionary_it.jsonl   reference split (R7)
      data/sources/kaikki/raw-wiktextract-data.jsonl.gz   if no split
      data/sources/kaikki/it-extract.jsonl.gz      Wikizionario
      data/sources/morphit/morph-it_048_utf8.txt   Morph-it! 0.4.8
      data/sources/leipzig/ita_news_2024_1M.tar.gz Leipzig corpus
    extract     data/italian/extract/italian_*.{jsonl,tsv}.gz
    lexicon     data/italian/italian_lexicon_raw.csv   (tracked, lean)
                data/italian/italian_lexicon_full.csv  (every column;
                                                        read by later stages)
                data/italian/italian_plural_forms.csv
                data/italian/italian_plural_drift.csv
    qc          data/italian/italian_qc.csv
                analysis/italian/italian_qc_report.txt (tracked)
                analysis/italian/italian_manual_worklist.csv (tracked)
    attest      data/italian/leipzig/…  (lexpipe.frequency layout)
                data/italian/italian_attestation.csv
                data/italian/italian_attestation_summary.tsv
    features    data/italian/italian_cells.csv
                data/italian/italian_lexicon_complete.csv
    derivations data/italian/italian_derivations.csv
                data/italian/italian_derivations_summary.tsv
    verbs       data/italian/italian_verbs.csv
                data/italian/italian_verbs_summary.tsv
    tp          analysis/italian/italian_tp_long.csv   (tracked)
                analysis/italian/italian_learners.csv  (tracked)
                analysis/italian/italian_tp_summary.txt (tracked)
    covariation analysis/italian/italian_covariation.txt (tracked)
    spellout    analysis/italian/italian_spellout.txt  (tracked)
    noevidence  analysis/italian/italian_noevidence.txt (tracked)
    logtable    analysis/italian/italian_logtable.csv  (tracked)
                data/italian/italian_logtable_attestation.csv
    docs        docs/italian/data_dict.md              (tracked)
                docs/italian/README.md (generated block; tracked)
    schema      docs/italian/schema.toml               (tracked, read)
    curated     curated/italian/plural_adjudications.csv
                curated/italian/suffixoid_classes.csv  (read-only)

If the reference split is absent, the ``extract`` stage splits the raw
dump into ``data/italian/kaikki/`` (never into ``data/sources/``, R1).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Final

from lexpipe import kaikki, sources, suffixoids
from lexpipe.context import RunContext
from lexpipe.frequency import LeipzigLayout, language_layout

RAW_DUMP_SOURCE: Final[str] = "kaikki_raw"
IT_EXTRACT_SOURCE: Final[str] = "kaikki_it_extract"
MORPHIT_SOURCE: Final[str] = "morphit"
LEIPZIG_SOURCE_PREFIX: Final[str] = "leipzig_"


@dataclass(frozen=True, slots=True)
class Layout:
    """Paths of the Italian pipeline for one run context."""

    ctx: RunContext

    # -- helpers -----------------------------------------------------------

    def _data(self, stem: str, suffix: str = ".csv") -> Path:
        return self.ctx.data_dir / self.ctx.language.file_name(stem, suffix)

    def _extract(self, stem: str, suffix: str) -> Path:
        name = self.ctx.language.file_name(stem, suffix)
        return self.extract_dir / name

    def _source(self, name: str) -> Path:
        return self.ctx.root / sources.SOURCES[name].path

    # -- sources -----------------------------------------------------------

    @property
    def raw_dump(self) -> Path:
        return self._source(RAW_DUMP_SOURCE)

    @property
    def reference_split(self) -> Path:
        """The per-language split under ``data/sources/kaikki/``."""
        directory = self.ctx.sources_dir / "kaikki"
        return kaikki.split_path(directory, self.ctx.language.kaikki_code)

    @property
    def local_split_dir(self) -> Path:
        return self.ctx.data_dir / "kaikki"

    @property
    def local_split(self) -> Path:
        """Where ``extract`` splits the raw dump if no reference split."""
        return kaikki.split_path(
            self.local_split_dir, self.ctx.language.kaikki_code
        )

    @property
    def en_split(self) -> Path:
        """The reference split if present, else the local one."""
        if self.reference_split.is_file():
            return self.reference_split
        return self.local_split

    @property
    def it_extract(self) -> Path:
        return self._source(IT_EXTRACT_SOURCE)

    @property
    def morphit(self) -> Path:
        return self._source(MORPHIT_SOURCE)

    @property
    def leipzig_source(self) -> sources.Source:
        return sources.SOURCES[
            LEIPZIG_SOURCE_PREFIX + self.ctx.language.leipzig_id
        ]

    @property
    def leipzig_archive(self) -> Path:
        return self.ctx.root / self.leipzig_source.path

    @property
    def leipzig(self) -> LeipzigLayout:
        return language_layout(self.ctx.language.name, self.ctx.root)

    @property
    def leipzig_freq(self) -> Path:
        return self.leipzig.freq_csv(self.ctx.language.leipzig_id)

    # -- extract -----------------------------------------------------------

    @property
    def extract_dir(self) -> Path:
        return self.ctx.data_dir / "extract"

    @property
    def en_lexemes(self) -> Path:
        return self._extract("en_lexemes", ".jsonl.gz")

    @property
    def en_verbs(self) -> Path:
        return self._extract("en_verbs", ".jsonl.gz")

    @property
    def it_entries(self) -> Path:
        return self._extract("it_entries", ".jsonl.gz")

    @property
    def form_index(self) -> Path:
        return self._extract("form_index", ".tsv.gz")

    @property
    def extract_summary(self) -> Path:
        return self._extract("extract_summary", ".tsv")

    # -- lexicon -----------------------------------------------------------

    @property
    def lexicon_raw(self) -> Path:
        """The tracked, lean raw lexicon (``lexicon.RAW_FIELDS``)."""
        return self._data("lexicon_raw")

    @property
    def lexicon_full(self) -> Path:
        """Every lexeme with every column of the stage
        (``lexicon.LEXICON_FIELDS``; gitignored), read by later stages."""
        return self._data("lexicon_full")

    @property
    def plural_forms(self) -> Path:
        return self._data("plural_forms")

    @property
    def plural_drift(self) -> Path:
        return self._data("plural_drift")

    @property
    def lexicon_summary(self) -> Path:
        return self._data("lexicon_summary", ".tsv")

    # -- qc ----------------------------------------------------------------

    @property
    def qc(self) -> Path:
        return self._data("qc")

    @property
    def qc_report(self) -> Path:
        return self.ctx.analysis_dir / self.ctx.language.file_name(
            "qc_report", ".txt"
        )

    @property
    def worklist(self) -> Path:
        return self.ctx.analysis_dir / self.ctx.language.file_name(
            "manual_worklist"
        )

    # -- attest ------------------------------------------------------------

    @property
    def attestation(self) -> Path:
        return self._data("attestation")

    @property
    def attestation_summary(self) -> Path:
        return self._data("attestation_summary", ".tsv")

    # -- features ----------------------------------------------------------

    @property
    def cells(self) -> Path:
        return self._data("cells")

    @property
    def lexicon_complete(self) -> Path:
        return self._data("lexicon_complete")

    # -- analyses (stages after features) ----------------------------------

    def _analysis(self, stem: str, suffix: str = ".csv") -> Path:
        return self.ctx.analysis_dir / self.ctx.language.file_name(
            stem, suffix
        )

    @property
    def derivations(self) -> Path:
        return self._data("derivations")

    @property
    def derivations_summary(self) -> Path:
        return self._data("derivations_summary", ".tsv")

    @property
    def verbs(self) -> Path:
        return self._data("verbs")

    @property
    def verbs_summary(self) -> Path:
        return self._data("verbs_summary", ".tsv")

    @property
    def tp_long(self) -> Path:
        return self._analysis("tp_long")

    @property
    def learners(self) -> Path:
        return self._analysis("learners")

    @property
    def tp_summary(self) -> Path:
        return self._analysis("tp_summary", ".txt")

    @property
    def covariation(self) -> Path:
        return self._analysis("covariation", ".txt")

    @property
    def spellout(self) -> Path:
        return self._analysis("spellout", ".txt")

    @property
    def noevidence(self) -> Path:
        return self._analysis("noevidence", ".txt")

    @property
    def logtable(self) -> Path:
        return self._analysis("logtable")

    @property
    def logtable_attestation(self) -> Path:
        return self._data("logtable_attestation")

    @property
    def schema(self) -> Path:
        return self.ctx.docs_dir / "schema.toml"

    @property
    def data_dict(self) -> Path:
        return self.ctx.docs_dir / "data_dict.md"

    @property
    def readme(self) -> Path:
        return self.ctx.docs_dir / "README.md"

    # -- curated (read-only) -----------------------------------------------

    @property
    def plural_adjudications(self) -> Path:
        return self.ctx.curated_dir / "plural_adjudications.csv"

    @property
    def suffixoid_classes(self) -> Path:
        """Curated -logo/-fago/-urgo decisions (R27; :mod:`lexpipe.suffixoids`)."""
        return self.ctx.curated_dir / suffixoids.DECISIONS_FILE
