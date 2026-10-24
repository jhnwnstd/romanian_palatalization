"""Paths of the Czech analysis stages, derived from a RunContext.

The data layer's paths live in :mod:`.layout`; the analysis stages
(design §4, §6; R14–R27) add these. Licensing decides the directory
(R9):

``data/czech/nc/`` (gitignored, CC BY-NC-SA or mixed with it)
    czech_derinet.sqlite          DeriNet projection (R7, R14)
    czech_derivations.csv         one row per base × derivative pair
    czech_lexeme_strata.csv       R22 strata of every noun/adjective
    czech_spellout_roots.csv      spell-out roots and their choices
    czech_noevidence_pairs.csv    no-evidence bases × derivatives
    czech_covariation_units.csv   base × suffix-pair units (R23)
``data/czech/`` (gitignored: per-word Leipzig counts, R9)
    czech_tp_roots.csv            per-root statuses behind the TP rows
    czech_logtable_full.csv       the -log table with corpus counts
    czech_noevidence_counterexamples.txt
                                  R21 counterexample words filtered by
                                  corpus attestation (≥ 1 token)
``analysis/czech/`` (tracked; aggregates only, written through the
licence gate where they are CSV)
    czech_derivations_summary.csv, czech_derivations_report.txt,
    czech_tp_long.csv, czech_learners.csv, czech_tp_summary.txt,
    czech_spellout.txt, czech_noevidence.txt, czech_covariation.txt,
    czech_logtable.csv, czech_logtable_summary.txt
``docs/czech/`` (tracked)
    data_dict.md                  generated from :mod:`.schema`
``curated/czech/`` (tracked, read-only)
    suffixoid_classes.csv         -log/-urg/-graf/-fil decisions (R27)

Source inputs come from the registry (:mod:`lexpipe.sources`): DeriNet
2.3 and the Czech Etymological Lexicon (czetyl).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Final

from lexpipe import sources, suffixoids
from lexpipe.context import RunContext

DERINET_SOURCE: Final[str] = "derinet"
CZETYL_SOURCE: Final[str] = "czetyl"


@dataclass(frozen=True, slots=True)
class AnalysisPaths:
    """Resolved paths of the analysis stages (module docstring)."""

    derinet_source: Path
    czetyl_source: Path
    derinet_db: Path
    derivations: Path
    lexeme_strata: Path
    derivations_summary: Path
    derivations_report: Path
    tp_long: Path
    learners: Path
    tp_summary: Path
    tp_roots: Path
    spellout: Path
    spellout_roots: Path
    noevidence: Path
    noevidence_pairs: Path
    noevidence_words: Path
    covariation: Path
    covariation_units: Path
    logtable: Path
    logtable_full: Path
    logtable_summary: Path
    suffixoid_classes: Path
    data_dict: Path


def analysis_paths(ctx: RunContext) -> AnalysisPaths:
    """All analysis paths under ``ctx.root``."""
    info = ctx.language
    data, nc, analysis = ctx.data_dir, ctx.nc_dir, ctx.analysis_dir

    def name(stem: str, suffix: str = ".csv") -> str:
        return info.file_name(stem, suffix)

    def source(registry_name: str) -> Path:
        return sources.local_path(sources.get_source(registry_name), ctx.root)

    return AnalysisPaths(
        derinet_source=source(DERINET_SOURCE),
        czetyl_source=source(CZETYL_SOURCE),
        derinet_db=nc / name("derinet", ".sqlite"),
        derivations=nc / name("derivations"),
        lexeme_strata=nc / name("lexeme_strata"),
        derivations_summary=analysis / name("derivations_summary"),
        derivations_report=analysis / name("derivations_report", ".txt"),
        tp_long=analysis / name("tp_long"),
        learners=analysis / name("learners"),
        tp_summary=analysis / name("tp_summary", ".txt"),
        tp_roots=data / name("tp_roots"),
        spellout=analysis / name("spellout", ".txt"),
        spellout_roots=nc / name("spellout_roots"),
        noevidence=analysis / name("noevidence", ".txt"),
        noevidence_pairs=nc / name("noevidence_pairs"),
        noevidence_words=data / name("noevidence_counterexamples", ".txt"),
        covariation=analysis / name("covariation", ".txt"),
        covariation_units=nc / name("covariation_units"),
        logtable=analysis / name("logtable"),
        logtable_full=data / name("logtable_full"),
        logtable_summary=analysis / name("logtable_summary", ".txt"),
        suffixoid_classes=ctx.curated_dir / suffixoids.DECISIONS_FILE,
        data_dict=ctx.docs_dir / "data_dict.md",
    )
