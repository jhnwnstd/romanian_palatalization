"""The validator, as a Python API.

Formerly the driver held all the business logic behind an argparse
CLI. This module extracts it into pure functions that return
immutable dataclasses — the caller controls loading (profile,
lexicon rows) and output (rendering, writing).

Public surface:

  - :func:`evaluate_row` — run the pipeline on one row.
  - :func:`run_validation` — evaluate an iterable of rows and
    aggregate into a :class:`ValidationReport`.
  - :func:`trace_lemmas` — find specific lemmas in the lexicon and
    return :class:`TraceResult` records with opt-in explanation,
    ordering search, and perturbation search.
  - :func:`render_report`, :func:`render_trace` — pure ``str``-returning
    renderers.
  - :func:`write_mismatch_csv`, :func:`render_distance_summary` —
    post-processing helpers.

Nothing here calls ``print``; nothing here reads ``sys.argv``. A
notebook or another script gets exactly the same data as the demo
driver.
"""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Final, Iterable, Iterator, Mapping, Sequence

from .analyses.romanian_palatalization import AnalysisProfile
from .diagnostics.categorise import (
    Categorisation,
    RefinedCategory,
    refine_category,
)
from .diagnostics.compare import (
    CompareStatus,
    compare_ipa,
    split_variants,
)
from .diagnostics.distance import DistanceScore, score as distance_score
from .diagnostics.explain import Explanation, explain_derivation
from .diagnostics.ordering import (
    OrderingSearchResult,
    OrderingStrategy,
    search_orderings,
)
from .diagnostics.perturb import (
    PerturbationKind,
    PerturbationReport,
    search_perturbations_for,
)
from .g2p import build_ur
from .lexicon import ClusterTag, LexRow, iter_lexicon_rows
from .pipeline import Derivation, RulePipeline
from .segments import Word, segments_to_ipa


# ---------------------------------------------------------------------------
# Result types
# ---------------------------------------------------------------------------


class MatchMode(StrEnum):
    """How strict predicted-vs-attested comparison should be.

    ``NORMALISED`` — apply the compare.py cascade (trailing glide,
    schwa reduction, monophthongisation…). Default: mismatches are
    real phonological differences, not transcription noise.

    ``STRICT`` — byte-for-byte match against pipe-separated variants.
    Useful for A/B comparison of two analyses where you want to see
    every character-level difference.
    """

    NORMALISED = "NORMALISED"
    STRICT = "STRICT"


@dataclass(frozen=True, slots=True)
class RowResult:
    """Everything the validator learned about one row."""

    row: LexRow
    ur: Word | None                        # None when build_ur failed
    derivation: Derivation | None
    predicted_sr: str                      # rendered IPA; "" if UR failed
    attested_first: str                    # first non-foreign variant
    boolean_predicted: bool                # rules predict alternation
    boolean_match: bool                    # prediction == row.mutation
    ipa_match: bool                        # SR matches any attested variant
    categorisation: Categorisation
    distance: DistanceScore | None         # None when derivation is None


@dataclass(frozen=True, slots=True)
class TraceResult:
    """A row's evaluation + optional diagnostic add-ons."""

    row: LexRow
    result: RowResult
    explanation: Explanation | None = None
    ordering: OrderingSearchResult | None = None
    perturbation: PerturbationReport | None = None
    found: bool = True                     # False if the lemma isn't in the lexicon


@dataclass(frozen=True, slots=True)
class ValidationReport:
    """Aggregate outcome of :func:`run_validation`.

    The full ``results`` list is retained so callers can slice, plot,
    or filter. Derived tables are precomputed so common queries
    (per-cell match rate, cluster breakdown) don't require re-scanning.
    """

    profile_name: str
    results: tuple[RowResult, ...]
    total: int
    boolean_matches: int
    ipa_matches: int
    by_cell: Mapping[tuple[str, str], Counter[str]]
    by_cluster: Mapping[str, Counter[str]]
    by_category: Mapping[RefinedCategory, tuple[RowResult, ...]]


# ---------------------------------------------------------------------------
# Cluster → report-bucket key
# ---------------------------------------------------------------------------


_CLUSTER_KEYS: Final[
    Mapping[tuple[ClusterTag, str], str]
] = {
    (ClusterTag.SC, "e"): "sc(ă)+e",
    (ClusterTag.SC, "i"): "sc(ă)+i",
    (ClusterTag.CT, "e"): "ct(ă)+e",
    (ClusterTag.CT, "i"): "ct(ă)+i",
    (ClusterTag.ST, "i"): "st(ă)+i",
}


def _cluster_key(row: LexRow) -> str | None:
    return _CLUSTER_KEYS.get((row.cluster, row.opportunity))


# ---------------------------------------------------------------------------
# Boolean-mutation flag
# ---------------------------------------------------------------------------


def _boolean_mutation(
    deriv: Derivation, palatalization_names: frozenset[str],
) -> bool:
    """True iff any palatalization rule fired in the derivation."""
    return any(
        step.fired and step.rule in palatalization_names
        for step in deriv.steps
    )


# ---------------------------------------------------------------------------
# Evaluate one row
# ---------------------------------------------------------------------------


def evaluate_row(
    row: LexRow,
    profile: AnalysisProfile,
    *,
    match_mode: MatchMode = MatchMode.NORMALISED,
) -> RowResult:
    """Run the pipeline on one row and score the outcome.

    Pure: no printing, no file I/O. Rebuilds the resolver from
    ``profile.inventory`` per call — cheap because the inventory is
    already loaded on the profile.
    """
    inv = profile.inventory
    ur = build_ur(row.ipa_lemma, row.opportunity, inv,
                  stem_final=row.stem_final)

    variants = split_variants(row.ipa_pl)
    attested_first = variants[0] if variants else ""

    if ur is None:
        # Can't derive an SR — use the analysis-declared oracle.
        predicted = profile.fallback_predict(row)
        boolean_match = predicted == row.mutation
        cat = refine_category(
            row, "", row.ipa_pl, predicted, ur_built=False,
        )
        return RowResult(
            row=row, ur=None, derivation=None,
            predicted_sr="", attested_first=attested_first,
            boolean_predicted=predicted,
            boolean_match=boolean_match,
            ipa_match=boolean_match,       # can't compare IPA — accept boolean
            categorisation=cat,
            distance=None,
        )

    pipeline = RulePipeline(
        rules=profile.rules,
        resolve=lambda label: inv.segment(label),
    )
    deriv = pipeline.derive(ur)
    inv_base = inv.base_segments()
    predicted_sr = segments_to_ipa(deriv.sr, inv_base)
    boolean_predicted = _boolean_mutation(
        deriv, profile.palatalization_rule_names,
    )
    boolean_match = boolean_predicted == row.mutation

    if match_mode is MatchMode.STRICT:
        cmp = compare_ipa(predicted_sr, row.ipa_pl, strict=True)
    else:
        cmp = compare_ipa(predicted_sr, row.ipa_pl)
    if cmp.status is CompareStatus.EMPTY_ATTESTED:
        ipa_match = boolean_match       # no data to compare against
    else:
        ipa_match = cmp.matched

    cat = refine_category(
        row, predicted_sr, row.ipa_pl,
        boolean_predicted, ur_built=True,
    )
    dist = distance_score(
        predicted_sr=predicted_sr,
        attested_field=row.ipa_pl,
        derivation=deriv,
        expected_fired=profile.expected_firings(row),
        inventory=inv,
        palatalization_rules=profile.palatalization_rule_names,
    )
    return RowResult(
        row=row, ur=ur, derivation=deriv,
        predicted_sr=predicted_sr, attested_first=attested_first,
        boolean_predicted=boolean_predicted,
        boolean_match=boolean_match,
        ipa_match=ipa_match,
        categorisation=cat,
        distance=dist,
    )


# ---------------------------------------------------------------------------
# Run over a batch
# ---------------------------------------------------------------------------


def run_validation(
    profile: AnalysisProfile,
    rows: Iterable[LexRow],
    *,
    match_mode: MatchMode = MatchMode.NORMALISED,
) -> ValidationReport:
    """Evaluate every row and return an aggregate report.

    ``rows`` is any iterable — typically ``iter_lexicon_rows(...)``.
    Streaming: the row iterator can be lazy, but the returned report
    materialises all ``RowResult``s so callers can slice.
    """
    results: list[RowResult] = []
    per_cell: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    per_cluster: dict[str, Counter[str]] = defaultdict(Counter)
    by_cat: dict[RefinedCategory, list[RowResult]] = defaultdict(list)
    total = 0
    bool_m = 0
    ipa_m = 0

    for row in rows:
        res = evaluate_row(row, profile, match_mode=match_mode)
        results.append(res)
        total += 1
        cell = (row.stem_final, row.opportunity)
        per_cell[cell]["match" if res.boolean_match else "mismatch"] += 1
        if res.boolean_match:
            bool_m += 1
        if res.ipa_match:
            ipa_m += 1
        if not (res.boolean_match and res.ipa_match):
            by_cat[res.categorisation.category].append(res)
        ck = _cluster_key(row)
        if ck is not None:
            per_cluster[ck]["match" if res.boolean_match else "mismatch"] += 1

    return ValidationReport(
        profile_name=profile.name,
        results=tuple(results),
        total=total,
        boolean_matches=bool_m,
        ipa_matches=ipa_m,
        by_cell={k: v for k, v in per_cell.items()},
        by_cluster={k: v for k, v in per_cluster.items()},
        by_category={k: tuple(v) for k, v in by_cat.items()},
    )


# ---------------------------------------------------------------------------
# Trace one or more lemmas
# ---------------------------------------------------------------------------


def trace_lemmas(
    lemmas: Iterable[str],
    profile: AnalysisProfile,
    csv_path: Path,
    *,
    include_explanation: bool = False,
    include_ordering: OrderingStrategy | None = None,
    include_perturbation: Sequence[PerturbationKind] | None = None,
    perturbation_limit: int | None = None,
    match_mode: MatchMode = MatchMode.NORMALISED,
) -> tuple[TraceResult, ...]:
    """Find each named lemma in the lexicon and return TraceResults.

    Optional diagnostic add-ons: pass ``include_explanation=True`` to
    attach a per-rule :class:`Explanation`; pass a strategy to
    ``include_ordering`` to attach an :class:`OrderingSearchResult`;
    pass a tuple of :class:`PerturbationKind` to
    ``include_perturbation`` to attach a :class:`PerturbationReport`.

    Lemmas not found in the target-stem-filtered lexicon get a
    :class:`TraceResult` with ``found=False`` and default fields.
    """
    wanted = list(lemmas)
    wanted_set = set(wanted)
    found: dict[str, TraceResult] = {}
    inv = profile.inventory
    pipeline = RulePipeline(
        rules=profile.rules,
        resolve=lambda label: inv.segment(label),
    )

    for row in iter_lexicon_rows(
        csv_path, target_stems=profile.target_stems,
    ):
        if row.lemma not in wanted_set:
            continue
        if row.lemma in found:
            continue      # first hit wins
        res = evaluate_row(row, profile, match_mode=match_mode)
        exp = None
        if include_explanation and res.ur is not None:
            exp = explain_derivation(pipeline, res.ur)
        ord_res = None
        if include_ordering is not None and res.ur is not None:
            ord_res = search_orderings(
                profile=profile,
                ur=res.ur,
                attested=row.ipa_pl,
                strategy=include_ordering,
            )
        pert = None
        if include_perturbation is not None and res.ur is not None:
            pert = search_perturbations_for(
                profile=profile, row=row,
                kinds=tuple(include_perturbation),
                limit=perturbation_limit,
            )
        found[row.lemma] = TraceResult(
            row=row, result=res,
            explanation=exp, ordering=ord_res, perturbation=pert,
            found=True,
        )
        if len(found) == len(wanted_set):
            break

    return tuple(
        found.get(lem, _missing_trace(lem)) for lem in wanted
    )


def _missing_trace(lemma: str) -> TraceResult:
    """Placeholder for a lemma the lexicon iterator didn't yield."""
    empty_row = LexRow(
        lemma=lemma, plural="", pos="", stem_final="",
        opportunity="", mutation=False,
        ipa_lemma="", ipa_pl="", exception_reason="",
        cluster=ClusterTag.NONE,
    )
    empty_res = RowResult(
        row=empty_row, ur=None, derivation=None,
        predicted_sr="", attested_first="",
        boolean_predicted=False, boolean_match=False,
        ipa_match=False,
        categorisation=Categorisation(
            RefinedCategory.UR_BUILD_FAILED, (), 999, "",
        ),
        distance=None,
    )
    return TraceResult(
        row=empty_row, result=empty_res, found=False,
    )


# ---------------------------------------------------------------------------
# Report rendering (pure)
# ---------------------------------------------------------------------------


def render_report(
    report: ValidationReport,
    *,
    include_distance_summary: bool = False,
    max_examples_per_category: int = 6,
) -> str:
    """Render a :class:`ValidationReport` as a multi-line string.

    Pure: no printing, no file I/O. The caller decides whether to
    ``print()`` or ``Path.write_text()`` the result.
    """
    lines: list[str] = []
    total = report.total
    lines.append(f"Analysis: {report.profile_name}")
    lines.append(f"Target-stem rows validated: {total:,}")
    if total:
        lines.append(
            f"BOOLEAN mutation match: {report.boolean_matches:,}/{total:,}  "
            f"({100 * report.boolean_matches / total:.2f}%)"
        )
        lines.append(
            f"SURFACE-IPA match: {report.ipa_matches:,}/{total:,}  "
            f"({100 * report.ipa_matches / total:.2f}%)"
        )
    lines.append("")

    lines.append("=" * 72)
    lines.append("Per (stem_final × opportunity)  bool match / mismatch / rate")
    lines.append("=" * 72)
    for key in sorted(report.by_cell.keys()):
        sf, opp = key
        m = report.by_cell[key]["match"]
        n = report.by_cell[key]["mismatch"]
        rate = 100 * m / (m + n) if (m + n) else 0.0
        flag = "  " if rate >= 95 else (" ?" if rate >= 80 else " !")
        lines.append(
            f"  ({sf!r}, {opp!r:6}){flag}  {m:6} / {n:6}   {rate:6.2f}%"
        )
    lines.append("")

    lines.append("=" * 72)
    lines.append("Cluster-specific verification")
    lines.append("=" * 72)
    for c in sorted(report.by_cluster.keys()):
        m = report.by_cluster[c]["match"]
        n = report.by_cluster[c]["mismatch"]
        rate = 100 * m / (m + n) if (m + n) else 0.0
        lines.append(f"  {c:8} {m:5} / {n:5}   {rate:6.2f}%")
    lines.append("")

    lines.append("=" * 72)
    lines.append(
        "Refined mismatch breakdown"
        " (NORM_* = data-side transcription noise; RULE_FAILURE = residual)"
    )
    lines.append("=" * 72)
    for cat in sorted(
        report.by_category.keys(),
        key=lambda c: -len(report.by_category[c]),
    ):
        items = report.by_category[cat]
        lines.append(f"\n  [{cat.value}]  n = {len(items)}")
        for res in items[:max_examples_per_category]:
            r = res.row
            verdict = (
                "rules→T, data→F" if res.boolean_predicted
                else "rules→F, data→T"
            )
            lines.append(
                f"    {r.lemma:18s} → {r.plural:18s}  "
                f"({r.stem_final}, opp={r.opportunity})  {verdict}"
            )
            if res.predicted_sr and res.attested_first:
                d = res.distance.total if res.distance else "-"
                lines.append(
                    f"      predicted={res.predicted_sr!r:22s} "
                    f"attested={res.attested_first!r:22s} d={d}"
                )
        if len(items) > max_examples_per_category:
            lines.append(
                f"    ... and {len(items) - max_examples_per_category} more"
            )

    if include_distance_summary:
        lines.append("")
        lines.append(render_distance_summary(report))

    return "\n".join(lines) + "\n"


def render_distance_summary(report: ValidationReport) -> str:
    """Histogram of total distance across mismatches."""
    mismatches = [
        r for r in report.results
        if not (r.boolean_match and r.ipa_match) and r.distance is not None
    ]
    if not mismatches:
        return "No mismatches to score."
    buckets: dict[str, int] = defaultdict(int)
    examples: dict[str, str] = {}
    for res in sorted(mismatches, key=lambda r: r.distance.total):
        d = res.distance.total
        if d == 0:
            key = "0.0 (data-side noise)"
        elif d <= 0.5:
            key = "0.0–0.5 (very close)"
        elif d <= 1.0:
            key = "0.5–1.0 (one feature off)"
        elif d <= 2.0:
            key = "1.0–2.0 (a few features off)"
        elif d <= 4.0:
            key = "2.0–4.0 (structural gap)"
        else:
            key = ">4.0 (far)"
        buckets[key] += 1
        examples.setdefault(key, res.row.lemma)
    lines = [
        "=" * 72,
        "Distance-to-working histogram",
        "=" * 72,
    ]
    order = [
        "0.0 (data-side noise)",
        "0.0–0.5 (very close)",
        "0.5–1.0 (one feature off)",
        "1.0–2.0 (a few features off)",
        "2.0–4.0 (structural gap)",
        ">4.0 (far)",
    ]
    total = sum(buckets.values())
    for key in order:
        n = buckets.get(key, 0)
        if not n:
            continue
        share = 100 * n / total
        ex = examples.get(key, "-")
        lines.append(
            f"  {key:34s}  n={n:5d}  ({share:5.1f}%)  ex: {ex}"
        )
    return "\n".join(lines)


def render_trace(
    trace: TraceResult,
    inv_base: Mapping,
) -> str:
    """Render a :class:`TraceResult` as a multi-line string.

    Includes the paper-style column derivation always; the
    explanation, ordering, and perturbation blocks only when the
    trace carries those add-ons.
    """
    if not trace.found:
        return f"{trace.row.lemma}: not found in filtered lexicon rows\n"
    result = trace.result
    if result.ur is None or result.derivation is None:
        return (
            f"{result.row.lemma}: could not build UR from IPA "
            f"{result.row.row.ipa_lemma!r}\n"
        )

    lines: list[str] = []
    r = result.row
    lines.append(
        f"{r.lemma} → {r.plural}   "
        f"(stem_final={r.stem_final}, opp={r.opportunity})"
    )
    lines.append(f"attested plural IPA: {result.attested_first!r}")
    if result.distance is not None:
        lines.append(
            f"distance: total={result.distance.total} "
            f"ipa_edit={result.distance.ipa_edit} "
            f"feature_edit={result.distance.feature_edit} "
            f"rule_edit={result.distance.rule_edit}"
        )
    lines.append("")

    def _row(label: str, word: Word) -> str:
        segs = " ".join(s.label for s in word)
        surface = segments_to_ipa(word, inv_base)
        return f"  {label:26s}  {segs:30s}  = /{surface}/"

    lines.append(_row("UR", result.derivation.ur))
    for step in result.derivation.steps:
        marker = "" if step.fired else "  (no-op)"
        lines.append(_row(step.rule + marker, step.word))
    lines.append(_row("SR", result.derivation.sr))
    lines.append("")
    matched = result.boolean_match and result.ipa_match
    verdict = "✓" if matched else "✗"
    tag = (
        "" if matched
        else f"   category={result.categorisation.category.value}"
    )
    lines.append(
        f"  match: boolean={result.boolean_match}  "
        f"ipa={result.ipa_match}  {verdict}{tag}"
    )

    if trace.explanation is not None:
        from .diagnostics.explain import format_explanation
        lines.append("")
        lines.append("--- per-rule explanation ---")
        lines.append(format_explanation(trace.explanation))

    if trace.ordering is not None:
        from .diagnostics.ordering import format_result as _fmt_ord
        lines.append("")
        lines.append("--- ordering search ---")
        lines.append(_fmt_ord(trace.ordering))

    if trace.perturbation is not None:
        from .diagnostics.perturb import format_report as _fmt_pert
        lines.append("")
        lines.append("--- perturbation search ---")
        lines.append(_fmt_pert(trace.perturbation))

    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Post-processing: mismatch CSV
# ---------------------------------------------------------------------------


def write_mismatch_csv(
    report: ValidationReport,
    path: Path,
    *,
    include_matches: bool = False,
) -> int:
    """Write mismatches (or every row) to CSV. Returns rows written.

    The only side-effect helper in this module. Every other function
    is pure. Callers who want the mismatch data in memory should
    just walk ``report.results``.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "lemma", "plural", "stem_final", "opportunity",
            "ur", "predicted_sr", "attested_sr",
            "category", "normalisations",
            "ipa_edit", "feature_edit", "rule_edit", "total_distance",
            "cluster",
        ])
        for res in report.results:
            if not include_matches and res.boolean_match and res.ipa_match:
                continue
            ur_repr = (
                "".join(s.label for s in res.ur)
                if res.ur is not None else ""
            )
            d = res.distance
            w.writerow([
                res.row.lemma,
                res.row.plural,
                res.row.stem_final,
                res.row.opportunity,
                ur_repr,
                res.predicted_sr,
                res.attested_first,
                res.categorisation.category.value,
                ",".join(
                    n.value for n in res.categorisation.normalisations
                ),
                d.ipa_edit if d else "",
                d.feature_edit if d else "",
                d.rule_edit if d else "",
                d.total if d else "",
                _cluster_key(res.row) or "",
            ])
            n += 1
    return n


__all__: Final[tuple[str, ...]] = (
    "MatchMode",
    "RowResult",
    "TraceResult",
    "ValidationReport",
    "evaluate_row",
    "render_distance_summary",
    "render_report",
    "render_trace",
    "run_validation",
    "trace_lemmas",
    "write_mismatch_csv",
)
