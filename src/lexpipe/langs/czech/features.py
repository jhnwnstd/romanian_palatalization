"""Stage ``features``: P/N endings, mutations and lexeme summaries (design §6).

Joins the cells (:mod:`.lexicon`), their attestation
(:mod:`.attestation`) and the curated cell decisions
(:mod:`.adjudication`, guards checked) and classifies every cell form of
every in-domain lexeme with the pure functions of :mod:`.ortho`:

    stem, stem_final, cluster, target_class   the oblique stem and its
                        final unit (Turek → Turk-, dům → dom-); velars,
                        r and the clusters sk ck Cr are core targets;
                        t/d/n are ``coronal_inferred``; c, s/z, labials
                        are recorded; l never is a target
    ending, orth_change, mutation, mutation_kind, palatal_output
                        the form split as (mutated) stem + ending
                        (k→c, h/g→z, ch→š, r→ř, sk→sc/šť, ck→čť; VOC
                        k→č, h/g→ž, c→č); t/d/n + ě/i/í is ``inferred``
    allomorph_class, n_kind, front_vowel_ending, opportunity
                        the §6 P/N table (:func:`lexpipe.langs.czech.ortho.classify`)
    exception_label     undergoer / unexplained / non exception
                        (allomorph) / nde:<class> / non exception
    target_is_suffix    R26: the target lies in -ník -ík -ák -ek -ka
                        -ička -ský -cký -log -urg (-graf has none),
                        identified morphologically (the lexicon's
                        ``suffix_evidence``, :mod:`.morphology`:
                        etymology template, suffix category,
                        derivational form-of sense, another lexeme's
                        headword line, the -ský/-cký adjectives, a base
                        lexeme; -log/-urg/-graf: spelling)
    count, attested_*, form_ambiguous   from the attest stage
    curated_status      confirmed / rejected by a hand decision

Lexeme summary
--------------
Per evidence set and τ, a lexeme whose stem ends in a target is an
``undergoer`` (it alternates in some attested P cell: a P ending with
the mutation), a ``resister`` (only N allomorphs that compete with a P
ending in their cell: -ové/-é, -u, -ách), ``mixed`` (both:
biologové + biolozích) or ``no_evidence``. Non-trigger endings (-em,
-ech, coronal VOC -e) are not evidence either way, nor is the
regular velar vocative -u. For c-final stems only the vocative counts
(otče). Only primary rows count, and not curated-rejected ones.
``summary_E0_<τ>_strict`` (corpus thresholds) also ignores rows whose
form is ambiguous with another lexeme (šachové: also the adjective
šachový) or fills a cell outside its copy group (``cell_confounded``:
LOC.SG parku is also GEN.SG and DAT.SG, so its tokens do not attest
the LOC.SG choice); NOM.PL/VOC.PL and feminine DAT.SG/LOC.SG are copy
groups. The
data layer computes evidence set E0 (nominal inflection); verbs (E1)
and derivatives (E2) belong to the derivation stages. Lexemes with
other finals get ``n/a``; out-of-domain ones ``out``.

Outputs
-------
``czech_cells.csv``           one row per in-domain cell form
``czech_lexicon_complete.csv`` every lexeme: raw columns + features,
                              frequencies and summaries
``czech_features_report.txt`` counts behind the summaries (untracked:
                              it depends on Leipzig counts)
"""

from __future__ import annotations

import itertools
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Final

from lexpipe import csvio, curated
from lexpipe.attest import TAUS
from lexpipe.context import RunContext, run_stage
from lexpipe.langs.czech import adjudication, jsonio, lexicon, ortho
from lexpipe.langs.czech.attestation import CORPUS_TAUS, FREQUENCY_FIELDS
from lexpipe.langs.czech.layout import CzechPaths, paths_for

STAGE: Final[str] = "features"
# 6: reads czech_lexicon_full.csv; cell adjudications keyed by
# (lexeme_id, cell).
# 7: GEN.SG-based stems with the palatal spelling (ohně → ohň-),
# plural-only strip endings, adjectival -é feminines; root-internal
# NDEB tokens; R26 evidence identified morphologically (audit
# 2026-10-23).
VERSION: Final[str] = "7"
EVIDENCE_SET: Final[str] = "E0"
SUMMARY_TARGETS: Final[frozenset[str]] = frozenset(
    {ortho.VELAR, ortho.R_CLASS, ortho.CORONAL, ortho.C_CLASS}
)
UNDERGOER: Final[str] = "undergoer"
RESISTER: Final[str] = "resister"
MIXED: Final[str] = "mixed"
NO_EVIDENCE: Final[str] = "no_evidence"
NOT_APPLICABLE: Final[str] = "n/a"
# Target classes whose per-cell P/N shares the report lists.
SHARE_TARGETS: Final[frozenset[str]] = frozenset({ortho.VELAR, ortho.R_CLASS})
OUT: Final[str] = "out"

CELL_FIELDS: Final[tuple[str, ...]] = (
    "lexeme_id",
    "root_id",
    "lemma",
    "pos",
    "cell",
    "form",
    "source",
    "primary",
    "table",
    "table_animacy",
    "gender_class",
    "form_source",
    "qualifiers",
    "syncretic_cells",
    "cell_confounded",
    "in_other_edition",
    "stem",
    "stem_final",
    "cluster",
    "target_class",
    "coronal_inferred",
    "target_is_suffix",
    "aligned",
    "ending",
    "orth_change",
    "mutation",
    "mutation_kind",
    "palatal_output",
    "allomorph_class",
    "n_kind",
    "front_vowel_ending",
    "opportunity",
    "exception_label",
    "count",
    "count_cs",
    "count_folded",
    *(f"attested_{tau}" for tau in CORPUS_TAUS),
    "form_ambiguous",
    "curated_status",
)
SUMMARY_FIELDS: Final[tuple[str, ...]] = (
    *itertools.chain.from_iterable(
        (
            f"summary_{EVIDENCE_SET}_{tau}",
            f"p_cells_{EVIDENCE_SET}_{tau}",
            f"n_cells_{EVIDENCE_SET}_{tau}",
        )
        for tau in TAUS
    ),
    *(f"summary_{EVIDENCE_SET}_{tau}_strict" for tau in CORPUS_TAUS),
)
# Cells whose forms are systematic copies of each other: a token of a
# form shared only within its group still attests the cell's choice.
COPY_GROUPS: Final[tuple[frozenset[str], ...]] = (
    frozenset({"NOM.PL", "VOC.PL"}),
    frozenset({"NOM.PL.M.AN", "VOC.PL.M.AN"}),
)
FEMININE_COPY: Final[frozenset[str]] = frozenset({"DAT.SG", "LOC.SG"})
FEATURE_FIELDS: Final[tuple[str, ...]] = (
    "stem",
    "stem_method",
    "declension",
    "stem_final",
    "cluster",
    "target_class",
    "coronal_inferred",
    "core_target",
    "target_is_suffix",
    "nde_class",
    "ndeb_tokens",
    "n_opportunity_rows",
    "n_unexplained",
    "n_unaligned",
    *(f for f in FREQUENCY_FIELDS if f not in ("lexeme_id", "lemma", "pos")),
    *SUMMARY_FIELDS,
)
COMPLETE_FIELDS: Final[tuple[str, ...]] = (
    *(
        f
        for f in lexicon.LEXICON_FIELDS
        if f not in ("source", "licence_class")
    ),
    *FEATURE_FIELDS,
    "source",
    "licence_class",
)


# ---------------------------------------------------------------------------
# Input streams (all sorted by lexeme_id)
# ---------------------------------------------------------------------------


def _grouped(
    rows: Iterable[dict[str, str]],
) -> Iterator[tuple[str, list[dict[str, str]]]]:
    for key, group in itertools.groupby(rows, key=lambda r: r["lexeme_id"]):
        yield key, list(group)


class _Aligned:
    """Look up the rows of a lexeme in a stream sorted by lexeme_id."""

    def __init__(self, rows: Iterable[dict[str, str]]) -> None:
        self._groups = _grouped(rows)
        self._current: tuple[str, list[dict[str, str]]] | None = next(
            self._groups, None
        )

    def take(self, lexeme_id: str) -> list[dict[str, str]]:
        while self._current is not None and self._current[0] < lexeme_id:
            self._current = next(self._groups, None)
        if self._current is not None and self._current[0] == lexeme_id:
            rows = self._current[1]
            self._current = next(self._groups, None)
            return rows
        return []


def _attestation_key(row: Mapping[str, str]) -> tuple[str, str, str]:
    return (row["cell"], row["form"], row["source"])


# ---------------------------------------------------------------------------
# One lexeme
# ---------------------------------------------------------------------------


@dataclass
class LexemeFeatures:
    """Feature rows and summary columns of one lexeme."""

    rows: list[dict[str, Any]] = field(default_factory=list)
    columns: dict[str, Any] = field(default_factory=dict)


def _primary_forms(
    cell_rows: Sequence[Mapping[str, str]],
) -> dict[str, list[str]]:
    forms: dict[str, list[str]] = {}
    for row in cell_rows:
        if csvio.parse_bool(row["primary"]) is True:
            found = forms.setdefault(row["cell"], [])
            if row["form"] not in found:
                found.append(row["form"])
    return forms


def _evidence(row: Mapping[str, Any], target: ortho.TargetInfo) -> str:
    """``P`` (alternating P), ``N`` (competing N allomorph) or ""."""
    if target.target_class == ortho.C_CLASS and row["cell"] != "VOC.SG":
        return ""
    if row["allomorph_class"] == ortho.P and row["mutation"] is True:
        return "P"
    if (
        row["allomorph_class"] == ortho.N
        and row["n_kind"] == ortho.N_ALLOMORPH
    ):
        return "N"
    return ""


def cell_confounded(cell: str, gender: str, syncretic: Sequence[str]) -> bool:
    """Does the form also fill a cell outside its copy group?

    Then its tokens do not attest this cell's choice: LOC.SG parku is
    also GEN.SG and DAT.SG. NOM.PL/VOC.PL copies and the feminine
    DAT.SG = LOC.SG pair do not count.
    """
    group: frozenset[str] = frozenset({cell})
    for copies in COPY_GROUPS:
        if cell in copies:
            group = copies
    if gender == "f" and cell in FEMININE_COPY:
        group = FEMININE_COPY
    return any(other not in group for other in syncretic)


def summarize(
    rows: Sequence[Mapping[str, Any]],
    target: ortho.TargetInfo,
    tau: str,
    *,
    strict: bool = False,
) -> tuple[str, list[str], list[str]]:
    """``(summary, P cells, N cells)`` of one lexeme at threshold *tau*.

    *strict* (corpus thresholds only) ignores rows whose form is
    ambiguous with another lexeme or confounded with another cell.
    """
    p_cells: list[str] = []
    n_cells: list[str] = []
    for row in rows:
        if (
            not row["primary"]
            or row["curated_status"] == adjudication.REJECTED
        ):
            continue
        if tau != "dictionary" and not row.get(f"attested_{tau}"):
            continue
        if strict and (row.get("form_ambiguous") or row["cell_confounded"]):
            continue
        kind = _evidence(row, target)
        cells = p_cells if kind == "P" else n_cells if kind == "N" else None
        if cells is not None and row["cell"] not in cells:
            cells.append(row["cell"])
    if p_cells and n_cells:
        summary = MIXED
    elif p_cells:
        summary = UNDERGOER
    elif n_cells:
        summary = RESISTER
    else:
        summary = NO_EVIDENCE
    return summary, p_cells, n_cells


def feature_row(
    lexeme: Mapping[str, str],
    cell_row: Mapping[str, str],
    analysis: ortho.StemAnalysis,
    suffix_inside: bool,
    nde: str,
    attestation: Mapping[str, str] | None,
    decided: Mapping[str, list[str] | None],
) -> dict[str, Any]:
    """One row of czech_cells.csv."""
    target = analysis.target
    pos = lexeme["pos"]
    gender = (
        "adj"
        if pos == "adj"
        else ortho.gender_class(
            csvio.split_list(lexeme["readings"]), cell_row["table_animacy"]
        )
    )
    alignment = ortho.align(analysis.stem, cell_row["form"])
    verdict = (
        ortho.classify(cell_row["cell"], gender, alignment.ending, target)
        if alignment.aligned
        else ortho.NO_OPPORTUNITY
    )
    change = ortho.mutation_of(alignment, target, verdict.allomorph_class)
    label = ortho.exception_label(
        in_domain=True,
        cell=cell_row["cell"],
        classification=verdict,
        target=target,
        mutation=change.mutation,
        nde_class=nde,
    )
    primary = csvio.parse_bool(cell_row["primary"]) is True
    row: dict[str, Any] = {
        "lexeme_id": lexeme["lexeme_id"],
        "root_id": lexeme["root_id"],
        "lemma": lexeme["lemma"],
        "pos": pos,
        "cell": cell_row["cell"],
        "form": cell_row["form"],
        "source": cell_row["source"],
        "primary": primary,
        "table": cell_row["table"],
        "table_animacy": cell_row["table_animacy"],
        "gender_class": gender,
        "form_source": cell_row["form_source"],
        "qualifiers": csvio.split_list(cell_row["qualifiers"]),
        "syncretic_cells": csvio.split_list(cell_row["syncretic_cells"]),
        "cell_confounded": cell_confounded(
            cell_row["cell"],
            gender,
            csvio.split_list(cell_row["syncretic_cells"]),
        ),
        "in_other_edition": csvio.parse_bool(cell_row["in_other_edition"]),
        "stem": analysis.stem,
        "stem_final": target.stem_final,
        "cluster": target.cluster,
        "target_class": target.target_class,
        "coronal_inferred": target.coronal_inferred,
        "target_is_suffix": suffix_inside,
        "aligned": alignment.aligned,
        "ending": alignment.ending,
        "orth_change": alignment.orth_change,
        "mutation": change.mutation,
        "mutation_kind": change.kind,
        "palatal_output": change.palatal_output,
        "allomorph_class": verdict.allomorph_class,
        "n_kind": verdict.n_kind,
        "front_vowel_ending": alignment.ending[:1] in ortho.FRONT_VOWELS,
        "opportunity": verdict.opportunity,
        "exception_label": label,
        "curated_status": (
            adjudication.curated_status(
                decided, cell_row["cell"], cell_row["form"]
            )
            if primary
            else ""
        ),
    }
    if attestation is not None:
        row["count"] = int(attestation["count"])
        row["count_cs"] = int(attestation["count_cs"])
        row["count_folded"] = int(attestation["count_folded"])
        for tau in CORPUS_TAUS:
            row[f"attested_{tau}"] = csvio.parse_bool(
                attestation[f"attested_{tau}"]
            )
        row["form_ambiguous"] = csvio.parse_bool(attestation["form_ambiguous"])
    return row


def lexeme_features(
    lexeme: Mapping[str, str],
    cell_rows: Sequence[Mapping[str, str]],
    attestations: Mapping[tuple[str, str, str], Mapping[str, str]],
    decided: Mapping[str, list[str] | None],
) -> LexemeFeatures:
    """Rows and summary columns of one lexeme."""
    forms = _primary_forms(cell_rows)
    flags = csvio.split_list(lexeme["spec_flags"])
    analysis = ortho.analyse_stem(lexeme["pos"], lexeme["lemma"], forms, flags)
    target = analysis.target
    ipa = csvio.split_list(lexeme["ipa"])
    suffix_inside = ortho.target_is_suffix(
        lexeme["lemma_suffix"], lexeme["suffix_evidence"], target
    )
    nde = ortho.nde_class(analysis.stem, target, ipa)
    out = LexemeFeatures()
    in_domain = lexeme["domain"] == lexicon.DOMAIN_IN
    if in_domain:
        out.rows = [
            feature_row(
                lexeme,
                row,
                analysis,
                suffix_inside,
                nde,
                attestations.get(_attestation_key(row)),
                decided,
            )
            for row in cell_rows
        ]
    columns: dict[str, Any] = {
        "stem": analysis.stem,
        "stem_method": analysis.method,
        "declension": analysis.declension,
        "stem_final": target.stem_final,
        "cluster": target.cluster,
        "target_class": target.target_class,
        "coronal_inferred": target.coronal_inferred,
        "core_target": target.core,
        "target_is_suffix": suffix_inside,
        "nde_class": nde,
        "ndeb_tokens": list(ortho.ndeb_tokens(lexeme["lemma"], ipa)),
        "n_opportunity_rows": sum(
            1 for r in out.rows if r["primary"] and r["opportunity"]
        ),
        "n_unexplained": sum(
            1 for r in out.rows if r["exception_label"] == ortho.UNEXPLAINED
        ),
        "n_unaligned": sum(1 for r in out.rows if not r["aligned"]),
    }
    for tau in TAUS:
        key = f"{EVIDENCE_SET}_{tau}"
        fallback = (
            OUT
            if not in_domain
            else (
                NOT_APPLICABLE
                if target.target_class not in SUMMARY_TARGETS
                else ""
            )
        )
        if fallback:
            columns[f"summary_{key}"] = fallback
            if tau != "dictionary":
                columns[f"summary_{key}_strict"] = fallback
        else:
            summary, p_cells, n_cells = summarize(out.rows, target, tau)
            columns[f"summary_{key}"] = summary
            columns[f"p_cells_{key}"] = p_cells
            columns[f"n_cells_{key}"] = n_cells
            if tau != "dictionary":
                columns[f"summary_{key}_strict"] = summarize(
                    out.rows, target, tau, strict=True
                )[0]
    out.columns = columns
    return out


# ---------------------------------------------------------------------------
# Curated decisions
# ---------------------------------------------------------------------------


def _guard_values(
    lexicon_rows: Mapping[str, Mapping[str, str]],
    cells_by_lexeme: Mapping[str, Mapping[str, list[str]]],
    wanted: Iterable[tuple[str, ...]],
) -> dict[tuple[str, ...], dict[str, Any]]:
    current: dict[tuple[str, ...], dict[str, Any]] = {}
    for key in wanted:
        lexeme_id, cell = key
        row = lexicon_rows.get(lexeme_id)
        if row is None:
            continue
        current[key] = adjudication.guard_fields(
            row["lemma"],
            csvio.split_list(row["gender"]),
            csvio.split_list(row["animacy"]),
            cells_by_lexeme.get(lexeme_id, {}).get(cell, []),
        )
    return current


def load_decisions(
    paths: CzechPaths, lexicon_rows: Mapping[str, Mapping[str, str]]
) -> dict[str, dict[str, list[str] | None]]:
    """Validated curated decisions; raises on stale guards."""
    rows = adjudication.load(paths.adjudications)
    if not rows:
        return {}
    for row in rows.values():
        adjudication.validate(row)
    wanted = {lexeme_id for lexeme_id, _ in rows}
    forms: dict[str, dict[str, list[str]]] = defaultdict(dict)
    for cell_row in jsonio.iter_csv(paths.cells_raw):
        if cell_row["lexeme_id"] in wanted:
            if csvio.parse_bool(cell_row["primary"]) is True:
                found = forms[cell_row["lexeme_id"]].setdefault(
                    cell_row["cell"], []
                )
                if cell_row["form"] not in found:
                    found.append(cell_row["form"])
    curated.check_guards(rows, _guard_values(lexicon_rows, forms, rows))
    return adjudication.decisions(rows)


# ---------------------------------------------------------------------------
# Stage
# ---------------------------------------------------------------------------


def run(ctx: RunContext) -> int | None:
    """Stage entry point (see the module docstring)."""
    paths = paths_for(ctx)
    inputs = [
        paths.lexicon_full,
        paths.cells_raw,
        paths.attestation,
        paths.lexeme_frequency,
    ]
    if paths.adjudications.is_file():
        inputs.append(paths.adjudications)
    return run_stage(
        ctx,
        STAGE,
        inputs=inputs,
        params={"version": VERSION, "evidence_set": EVIDENCE_SET},
        outputs=[paths.cells, paths.lexicon_complete, paths.features_report],
        body=lambda: _features(ctx, paths),
    )


@dataclass
class _Report:
    """Counts for the features report."""

    labels: Counter[tuple[str, str]] = field(default_factory=Counter)
    form_sources: Counter[tuple[str, str]] = field(default_factory=Counter)
    summaries: Counter[tuple[str, str, str]] = field(default_factory=Counter)
    targets: Counter[tuple[str, str]] = field(default_factory=Counter)
    readings: Counter[tuple[str, str]] = field(default_factory=Counter)
    share_cells: Counter[tuple[str, str, str, str, str]] = field(
        default_factory=Counter
    )
    coronal_p: Counter[str] = field(default_factory=Counter)
    unaligned: int = 0

    def add(self, lexeme: Mapping[str, str], result: LexemeFeatures) -> None:
        if lexeme["domain"] != lexicon.DOMAIN_IN or lexeme["pos"] not in (
            "noun",
            "adj",
        ):
            return
        pos = lexeme["pos"]
        target = str(result.columns["target_class"])
        self.targets[(pos, target)] += 1
        self.readings[(pos, lexeme["readings"] or "?")] += 1
        for tau in TAUS:
            summary = str(result.columns[f"summary_{EVIDENCE_SET}_{tau}"])
            self.summaries[(target, tau, summary)] += 1
            if tau != "dictionary":
                strict = str(
                    result.columns[f"summary_{EVIDENCE_SET}_{tau}_strict"]
                )
                self.summaries[(target, f"{tau}_strict", strict)] += 1
        per_cell: dict[tuple[str, str], dict[str, set[str]]] = defaultdict(
            lambda: defaultdict(set)
        )
        for row in result.rows:
            if not row["primary"]:
                continue
            self.labels[(target, str(row["exception_label"]))] += 1
            self.form_sources[(row["source"], str(row["form_source"]))] += 1
            self.unaligned += int(not row["aligned"])
            if target == ortho.CORONAL and row["allomorph_class"] == ortho.P:
                self.coronal_p[str(row["exception_label"])] += 1
            competing = row["allomorph_class"] == ortho.P or (
                row["n_kind"] == ortho.N_ALLOMORPH
            )
            if target not in SHARE_TARGETS or not competing:
                continue
            key = (row["cell"], row["gender_class"])
            for tau in ("dictionary", "ge1"):
                if tau == "ge1" and not row.get("attested_ge1"):
                    continue
                per_cell[key][tau].add(str(row["allomorph_class"]))
        for (cell, gender), by_tau in per_cell.items():
            for tau, classes in by_tau.items():
                kind = (
                    "both"
                    if classes == {ortho.P, ortho.N}
                    else f"{''.join(sorted(classes))} only"
                )
                self.share_cells[(target, cell, gender, tau, kind)] += 1

    def render(self) -> str:
        lines = [
            "Czech features (evidence set E0); untracked: depends on Leipzig (R9)",
            "",
            "In-domain lexemes by part of speech and reading",
        ]
        for (pos, reading), count in sorted(self.readings.items()):
            lines.append(f"  {pos:<5} {reading:<14} {count:>7,}")
        lines += ["", "Target class of the stem final"]
        for (pos, target), count in sorted(self.targets.items()):
            lines.append(f"  {pos:<5} {target:<18} {count:>7,}")
        lines += ["", "Primary rows by source and form_source (R11)"]
        for (source, origin), count in sorted(self.form_sources.items()):
            lines.append(f"  {source:<15} {origin:<15} {count:>9,}")
        lines += ["", "Exception labels of primary rows by target class"]
        for (target, label), count in sorted(self.labels.items()):
            lines.append(f"  {target:<18} {label:<28} {count:>9,}")
        lines += [
            "",
            "Coronal P rows (inferred; 100% undergoer by construction)",
            *(
                f"  {label:<28} {count:>9,}"
                for label, count in sorted(self.coronal_p.items())
            ),
            "",
            f"Unaligned primary rows (form not stem + ending): {self.unaligned:,}",
            "",
            "Lexemes per cell by the competing allomorphs present (P, and N",
            "allomorphs -ové/-é/-u/-ách; non-trigger endings left out)",
        ]
        keys = sorted({(t, c, g) for t, c, g, _, _ in self.share_cells})
        for target, cell, gender in keys:
            for tau in ("dictionary", "ge1"):
                parts = {
                    kind: self.share_cells[(target, cell, gender, tau, kind)]
                    for kind in ("P only", "N only", "both")
                }
                total = sum(parts.values())
                if not total:
                    continue
                shares = "  ".join(
                    f"{kind} {count:>6,} ({100 * count / total:5.1f}%)"
                    for kind, count in parts.items()
                )
                lines.append(
                    f"  {target:<6} {cell:<12} {gender:<5} {tau:<10} "
                    f"n={total:>6,}  {shares}"
                )
        lines += ["", "Lexeme summaries (E0) by target class and threshold"]
        taus = [*TAUS, *(f"{tau}_strict" for tau in CORPUS_TAUS)]
        for target in sorted({t for t, _, _ in self.summaries}):
            for tau in taus:
                counts = {
                    s: self.summaries[(target, tau, s)]
                    for s in (
                        UNDERGOER,
                        RESISTER,
                        MIXED,
                        NO_EVIDENCE,
                        NOT_APPLICABLE,
                    )
                    if self.summaries[(target, tau, s)]
                }
                text = "  ".join(f"{s} {c:,}" for s, c in counts.items())
                lines.append(f"  {target:<18} {tau:<10} {text}")
        return "\n".join(lines) + "\n"


def _features(ctx: RunContext, paths: CzechPaths) -> int:
    lexicon_rows = {
        row["lexeme_id"]: row for row in jsonio.iter_csv(paths.lexicon_full)
    }
    decisions = load_decisions(paths, lexicon_rows)
    frequencies = {
        row["lexeme_id"]: row
        for row in jsonio.iter_csv(paths.lexeme_frequency)
    }
    cells_stream = _Aligned(jsonio.iter_csv(paths.cells_raw))
    attest_stream = _Aligned(jsonio.iter_csv(paths.attestation))
    report = _Report()
    complete: list[dict[str, Any]] = []

    def rows() -> Iterator[dict[str, Any]]:
        for lexeme_id in sorted(lexicon_rows):
            lexeme = lexicon_rows[lexeme_id]
            attestations = {
                _attestation_key(r): r for r in attest_stream.take(lexeme_id)
            }
            result = lexeme_features(
                lexeme,
                cells_stream.take(lexeme_id),
                attestations,
                decisions.get(lexeme_id, {}),
            )
            report.add(lexeme, result)
            complete.append(
                _complete_row(lexeme, result, frequencies.get(lexeme_id))
            )
            yield from result.rows

    csvio.write_rows(paths.cells, rows(), CELL_FIELDS, ctx=ctx)
    csvio.write_rows(
        paths.lexicon_complete, complete, COMPLETE_FIELDS, ctx=ctx
    )
    text = report.render()
    ctx.write_text(paths.features_report, text)
    print(text, flush=True)
    return 0


def _complete_row(
    lexeme: Mapping[str, str],
    result: LexemeFeatures,
    frequency: Mapping[str, str] | None,
) -> dict[str, Any]:
    row: dict[str, Any] = {key: lexeme[key] for key in lexicon.LEXICON_FIELDS}
    row.update(result.columns)
    if frequency is not None:
        for key in FREQUENCY_FIELDS:
            if key not in ("lexeme_id", "lemma", "pos"):
                row[key] = frequency[key]
    return row
