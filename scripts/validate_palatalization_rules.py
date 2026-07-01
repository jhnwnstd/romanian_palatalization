#!/usr/bin/env python3
"""Validate the paper's palatalization rules against the full lexicon.

Runs the analysis in ``src/phonology/analyses/romanian_palatalization``
— the paper's revised rule set (place-class S-palatalization, feature-
deleting bleed, broad-terminator dorsal palatalization, multiply-
articulated /j/) — over every noun in
``data/romanian_lexicon_with_freq.csv`` with a target stem-final
consonant. For each row the driver builds an underlying representation
from the lemma's IPA, feeds it through the rule pipeline, and compares
the resulting surface form to the attested plural IPA using the shared
:mod:`phonology.diagnostics.compare` matcher.

Output
------
Under ``analysis/``:

  - ``palatalization_rule_validation.txt`` — per-cell / cluster /
    refined-mismatch report.
  - ``palatalization_mismatches.csv`` — every mismatch with predicted
    vs attested IPA, refined category, distance score, and applied
    normalisations for downstream inspection.

Trace-mode diagnostics (see ``--trace``, ``--trace-mode``):

  - ``explain`` — per-rule firing/blocking narrative (why did rule X
    fire or not fire?).
  - ``order`` — search rule orderings for one that produces the
    attested SR; distinguishes ordering bugs from content bugs.
  - ``perturb`` — search single-feature perturbations (patches,
    clear sets, rule supplies) that would flip failure to success.
  - ``all`` — all three, plus the paper-style column table.

Usage
-----
  python scripts/validate_palatalization_rules.py
  python scripts/validate_palatalization_rules.py --limit 500
  python scripts/validate_palatalization_rules.py --trace muscă,prost
  python scripts/validate_palatalization_rules.py --trace amant --trace-mode explain
  python scripts/validate_palatalization_rules.py --trace amant --trace-mode order
  python scripts/validate_palatalization_rules.py --trace amant --trace-mode perturb
  python scripts/validate_palatalization_rules.py --distance-summary
  python scripts/validate_palatalization_rules.py --strict-ipa    # disable IPA normalization
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Final, Iterator

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PROJECT_ROOT / "src"))

from phonology.analyses.romanian_palatalization import (  # noqa: E402
    PATCHES,
    RULES_FROM_PAPER,
    UNDERSPEC,
    load_inventory,
)
from phonology.diagnostics.categorise import (  # noqa: E402
    Categorisation,
    RefinedCategory,
    refine_category,
)
from phonology.diagnostics.compare import compare, split_variants  # noqa: E402
from phonology.diagnostics.distance import (  # noqa: E402
    PALATALIZATION_RULE_NAMES,
    DistanceScore,
    score as distance_score,
)
from phonology.diagnostics.explain import (  # noqa: E402
    explain_derivation,
    format_explanation,
)
from phonology.diagnostics.ordering import (  # noqa: E402
    OrderingStrategy,
    find_valid_orderings,
    format_result as format_ordering_result,
)
from phonology.diagnostics.perturb import (  # noqa: E402
    PerturbationKind,
    format_report as format_perturb_report,
    search_perturbations,
)
from phonology.g2p import build_ur, first_variant  # noqa: E402
from phonology.pipeline import Derivation, RulePipeline  # noqa: E402
from phonology.segments import Word, segments_to_ipa  # noqa: E402


CSV_PATH: Final[Path] = (
    _PROJECT_ROOT / "data" / "romanian_lexicon_with_freq.csv"
)
REPORT_PATH: Final[Path] = (
    _PROJECT_ROOT / "analysis" / "palatalization_rule_validation.txt"
)
DEFAULT_MISMATCH_CSV: Final[Path] = (
    _PROJECT_ROOT / "analysis" / "palatalization_mismatches.csv"
)
INVENTORY_JSON: Final[Path] = (
    _PROJECT_ROOT / "romanian_features.json"
)

TARGET_STEMS: Final[frozenset[str]] = frozenset(
    {"c", "g", "t", "d", "s", "z"}
)
DORSAL_STEMS: Final[frozenset[str]] = frozenset({"c", "g"})
# Coronals excluding /z/: the paper treats /z/ as fully specified and
# inalterable (latex.tex:170), so it never triggers the rules.
CORONAL_STEMS: Final[frozenset[str]] = frozenset({"t", "d", "s"})


# ---------------------------------------------------------------------------
# Row and result types
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class NounRow:
    """One lexicon row narrowed to what validation needs."""

    lemma: str
    plural: str
    stem_final: str
    opportunity: str
    mutation: bool
    ipa_lemma: str
    ipa_pl: str

    @property
    def has_ct_cluster(self) -> bool:
        lem = self.lemma.lower()
        return lem.endswith(("ct", "ctă"))

    @property
    def has_sc_cluster(self) -> bool:
        lem = self.lemma.lower()
        return lem.endswith(("sc", "scă")) and not lem.endswith(
            ("șcă", "şcă")
        )

    @property
    def has_shca_cluster(self) -> bool:
        lem = self.lemma.lower()
        return lem.endswith(("șcă", "şcă"))

    @property
    def has_st_cluster(self) -> bool:
        lem = self.lemma.lower()
        return lem.endswith(("st", "stă", "ște", "şte"))


@dataclass(frozen=True, slots=True)
class RowResult:
    """Outcome of running the pipeline on one row."""

    row: NounRow
    ur: Word | None
    derivation: Derivation | None
    predicted_sr: str
    attested_first: str
    boolean_predicted: bool
    boolean_match: bool
    ipa_match: bool
    categorisation: Categorisation
    distance: DistanceScore | None


# ---------------------------------------------------------------------------
# Expected-firing table
# ---------------------------------------------------------------------------


def _expected_firings(row: NounRow) -> frozenset[str]:
    """What palatalisation rules SHOULD fire for this row, given the
    paper's productivity claims. Used by the distance metric to score
    "rule firing" divergence.
    """
    if row.opportunity in {"uri", "none"}:
        return frozenset()
    if row.has_ct_cluster:
        return frozenset({"assibilation"}) if row.opportunity == "i" else frozenset()
    if row.has_sc_cluster:
        return frozenset({"dorsal-pal", "s-pal-rev", "bleed-rev"})
    if row.has_st_cluster and row.opportunity == "i":
        return frozenset({"s-pal-rev", "bleed-rev"})
    if row.stem_final in DORSAL_STEMS:
        if row.opportunity in {"i", "e"}:
            return frozenset({"dorsal-pal"})
    if row.stem_final in CORONAL_STEMS and row.opportunity == "i":
        if row.stem_final in {"t", "d"}:
            return frozenset({"assibilation"})
        return frozenset({"s-pal-rev"})
    return frozenset()


def _fallback_predict(row: NounRow) -> bool:
    """What the rules WOULD predict for boolean mutation when we can't
    build a real UR (missing / unparseable IPA).
    """
    if row.opportunity in {"uri", "none"}:
        return False
    if row.has_ct_cluster:
        return row.opportunity == "i"
    if row.has_sc_cluster or row.has_shca_cluster:
        return True
    if row.has_st_cluster and row.opportunity == "i":
        return True
    if row.stem_final in DORSAL_STEMS:
        return row.opportunity in {"i", "e"}
    if row.stem_final in CORONAL_STEMS:
        return row.opportunity == "i"
    return False


def _boolean_mutation(deriv: Derivation) -> bool:
    """Did any palatalisation rule fire during the derivation?"""
    for step in deriv.steps:
        if step.fired and step.rule in PALATALIZATION_RULE_NAMES:
            return True
    return False


# ---------------------------------------------------------------------------
# CSV streaming
# ---------------------------------------------------------------------------


def _load_rows(limit: int | None = None) -> Iterator[NounRow]:
    """Stream the lexicon and yield ``NounRow`` for validation."""
    yielded = 0
    with CSV_PATH.open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r.get("pos", "").strip() != "N":
                continue
            sf = (r.get("stem_final") or "").strip()
            if sf not in TARGET_STEMS:
                continue
            plural = (r.get("plural") or "").strip()
            if not plural:
                continue
            opp = (r.get("opportunity") or "").strip()
            if opp == "none":
                continue
            if (r.get("exception_reason") or "").startswith("nde:"):
                continue
            yield NounRow(
                lemma=r["lemma"].strip(),
                plural=plural,
                stem_final=sf,
                opportunity=opp,
                mutation=(r.get("mutation") or "").strip().upper() == "TRUE",
                ipa_lemma=(r.get("ipa_normalized_lemma") or "").strip(),
                ipa_pl=(r.get("ipa_normalized_pl") or "").strip(),
            )
            yielded += 1
            if limit is not None and yielded >= limit:
                return


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------


def _ipa_match_strict(predicted: str, attested_field: str) -> bool:
    """Byte-for-byte match against any pipe-separated variant.

    Used when ``--strict-ipa`` is set: no normalisation, no glide
    folding, no vowel-reduction pass. Any transcription difference
    counts as a mismatch.
    """
    if not predicted or not attested_field:
        return False
    return predicted in split_variants(attested_field)


def _evaluate(
    row: NounRow,
    pipeline: RulePipeline,
    inventory,
    strict_ipa: bool = False,
) -> RowResult:
    ur = build_ur(
        row.ipa_lemma, row.opportunity, inventory,
        stem_final=row.stem_final,
    )
    attested_first = first_variant(row.ipa_pl) if row.ipa_pl else ""

    if ur is None:
        fallback_predicted = _fallback_predict(row)
        boolean_match = fallback_predicted == row.mutation
        cat = refine_category(
            row, "", row.ipa_pl, fallback_predicted, ur_built=False,
        )
        return RowResult(
            row=row, ur=None, derivation=None,
            predicted_sr="", attested_first=attested_first,
            boolean_predicted=fallback_predicted,
            boolean_match=boolean_match,
            ipa_match=boolean_match,
            categorisation=cat,
            distance=None,
        )

    deriv = pipeline.derive(ur)
    inv_base = inventory.base_segments()
    predicted_sr = segments_to_ipa(deriv.sr, inv_base)
    boolean_predicted = _boolean_mutation(deriv)
    boolean_match = boolean_predicted == row.mutation

    if strict_ipa:
        ipa_match = _ipa_match_strict(predicted_sr, row.ipa_pl)
    else:
        ipa_match = (
            compare(predicted_sr, row.ipa_pl).matched
            if row.ipa_pl else boolean_match
        )

    cat = refine_category(
        row, predicted_sr, row.ipa_pl,
        boolean_predicted, ur_built=True,
    )
    dist = distance_score(
        predicted_sr=predicted_sr,
        attested_field=row.ipa_pl,
        derivation=deriv,
        expected_fired=_expected_firings(row),
        inventory=inventory,
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
# Cluster keys
# ---------------------------------------------------------------------------


def _cluster_key(row: NounRow) -> str | None:
    if row.has_sc_cluster and row.opportunity == "e":
        return "sc(ă)+e"
    if row.has_sc_cluster and row.opportunity == "i":
        return "sc(ă)+i"
    if row.has_ct_cluster and row.opportunity == "e":
        return "ct(ă)+e"
    if row.has_ct_cluster and row.opportunity == "i":
        return "ct(ă)+i"
    if row.has_st_cluster and row.opportunity == "i":
        return "st(ă)+i"
    return None


# ---------------------------------------------------------------------------
# Report rendering
# ---------------------------------------------------------------------------


def _render_report(results: list[RowResult]) -> str:
    lines: list[str] = []
    per_cell: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    cluster_counts: dict[str, Counter[str]] = defaultdict(Counter)
    by_cat: dict[RefinedCategory, list[RowResult]] = defaultdict(list)
    bool_match = 0
    ipa_match = 0
    total = 0

    for res in results:
        total += 1
        cell = (res.row.stem_final, res.row.opportunity)
        per_cell[cell]["match" if res.boolean_match else "mismatch"] += 1
        if res.boolean_match:
            bool_match += 1
        if res.ipa_match:
            ipa_match += 1
        if not (res.boolean_match and res.ipa_match):
            by_cat[res.categorisation.category].append(res)
        c_key = _cluster_key(res.row)
        if c_key is not None:
            cluster_counts[c_key][
                "match" if res.boolean_match else "mismatch"
            ] += 1

    lines.append(f"Target-stem rows validated: {total:,}")
    if total:
        lines.append(
            f"BOOLEAN mutation match: {bool_match:,}/{total:,}  "
            f"({100 * bool_match / total:.2f}%)"
        )
        lines.append(
            f"SURFACE-IPA match: {ipa_match:,}/{total:,}  "
            f"({100 * ipa_match / total:.2f}%)"
        )
    lines.append("")

    lines.append("=" * 72)
    lines.append("Per (stem_final × opportunity)  bool match / mismatch / rate")
    lines.append("=" * 72)
    for key in sorted(per_cell.keys()):
        sf, opp = key
        m = per_cell[key]["match"]
        n = per_cell[key]["mismatch"]
        rate = 100 * m / (m + n) if (m + n) else 0.0
        flag = "  " if rate >= 95 else (" ?" if rate >= 80 else " !")
        lines.append(
            f"  ({sf!r}, {opp!r:6}){flag}  {m:6} / {n:6}   {rate:6.2f}%"
        )
    lines.append("")

    lines.append("=" * 72)
    lines.append("Cluster-specific verification")
    lines.append("=" * 72)
    for c in sorted(cluster_counts.keys()):
        m = cluster_counts[c]["match"]
        n = cluster_counts[c]["mismatch"]
        rate = 100 * m / (m + n) if (m + n) else 0.0
        lines.append(f"  {c:8} {m:5} / {n:5}   {rate:6.2f}%")
    lines.append("")

    lines.append("=" * 72)
    lines.append(
        "Refined mismatch breakdown"
        " (NORM_* = data-side transcription noise; RULE_FAILURE = residual)"
    )
    lines.append("=" * 72)
    for cat in sorted(by_cat.keys(), key=lambda c: -len(by_cat[c])):
        items = by_cat[cat]
        lines.append(f"\n  [{cat.value}]  n = {len(items)}")
        for res in items[:6]:
            r = res.row
            verdict = (
                "rules→T, data→F"
                if res.boolean_predicted
                else "rules→F, data→T"
            )
            lines.append(
                f"    {r.lemma:18s} → {r.plural:18s}  "
                f"({r.stem_final}, opp={r.opportunity})  {verdict}"
            )
            if res.predicted_sr and res.attested_first:
                dist = res.distance.total if res.distance else "-"
                lines.append(
                    f"      predicted={res.predicted_sr!r:22s} "
                    f"attested={res.attested_first!r:22s} d={dist}"
                )
        if len(items) > 6:
            lines.append(f"    ... and {len(items) - 6} more")

    return "\n".join(lines) + "\n"


def _render_distance_summary(results: list[RowResult]) -> str:
    """Histogram of total distance across mismatches, plus mean stats."""
    mismatches = [
        r for r in results
        if not (r.boolean_match and r.ipa_match) and r.distance is not None
    ]
    if not mismatches:
        return "\nNo mismatches to score.\n"
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
    lines = ["", "=" * 72, "Distance-to-working histogram", "=" * 72]
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
        lines.append(f"  {key:34s}  n={n:5d}  ({share:5.1f}%)  ex: {ex}")
    return "\n".join(lines) + "\n"


def _write_mismatch_csv(
    results: list[RowResult],
    path: Path,
) -> int:
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
        for res in results:
            if res.boolean_match and res.ipa_match:
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


# ---------------------------------------------------------------------------
# Trace mode + diagnostics
# ---------------------------------------------------------------------------


def _render_trace(result: RowResult, inv_base) -> str:
    if result.ur is None or result.derivation is None:
        return (
            f"{result.row.lemma}: could not build UR from IPA "
            f"{result.row.ipa_lemma!r}\n"
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
    return "\n".join(lines) + "\n"


def _explain_trace(result: RowResult, pipeline: RulePipeline) -> str:
    if result.ur is None:
        return "  (no UR to explain)\n"
    exp = explain_derivation(pipeline, result.ur)
    return (
        "\n--- per-rule explanation ---\n" + format_explanation(exp) + "\n"
    )


def _ordering_trace(
    result: RowResult, inventory, strategy: OrderingStrategy,
) -> str:
    if result.ur is None:
        return "  (no UR — cannot search orderings)\n"
    res = find_valid_orderings(
        rules=RULES_FROM_PAPER,
        ur=result.ur,
        expected_field=result.row.ipa_pl,
        resolve=lambda l: inventory.segment(l),
        inventory=inventory,
        strategy=strategy,
    )
    return "\n--- ordering search ---\n" + format_ordering_result(res) + "\n"


def _perturb_trace(result: RowResult, opportunity: str, stem_final: str) -> str:
    if result.ur is None:
        return "  (no UR — cannot perturb)\n"
    ipa_lemma = result.row.ipa_lemma

    def builder(inv):
        return build_ur(ipa_lemma, opportunity, inv, stem_final=stem_final)

    report = search_perturbations(
        ur_builder=builder,
        expected_field=result.row.ipa_pl,
        inventory_json=INVENTORY_JSON,
        baseline_patches=PATCHES,
        baseline_underspec=UNDERSPEC,
        baseline_rules=RULES_FROM_PAPER,
        kinds=tuple(PerturbationKind),
        limit=None,
    )
    return "\n--- perturbation search ---\n" + format_perturb_report(report) + "\n"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--limit", type=int, default=None,
                   help="Stop after N rows.")
    p.add_argument("--trace", type=str, default="",
                   help="Comma-separated lemmas to trace.")
    p.add_argument(
        "--trace-mode", type=str, default="table",
        choices=("table", "explain", "order", "perturb", "all"),
        help=(
            "What to include in --trace output. "
            "table = paper-style column derivation (default); "
            "explain = per-rule firing/blocking diagnostic; "
            "order = search rule orderings; "
            "perturb = search single-feature perturbations; "
            "all = everything."
        ),
    )
    p.add_argument(
        "--order-strategy", type=str, default="pairwise",
        choices=("declared", "adjacent", "pairwise", "full"),
        help="Ordering search strategy (used only in --trace-mode order/all).",
    )
    p.add_argument(
        "--strict-ipa", action="store_true",
        help=(
            "Disable IPA normalisation in match-checking (no trailing "
            "glide fold, no schwa/monophthong tolerance)."
        ),
    )
    p.add_argument(
        "--distance-summary", action="store_true",
        help="Append a distance-to-working histogram to the report.",
    )
    p.add_argument(
        "--mismatch-csv", type=Path, default=DEFAULT_MISMATCH_CSV,
        help="Where to write the mismatch CSV.",
    )
    p.add_argument(
        "--no-mismatch-csv", action="store_true",
        help="Skip writing the mismatch CSV.",
    )
    return p.parse_args(argv)


_STRATEGY_MAP: Final[dict[str, OrderingStrategy]] = {
    "declared": OrderingStrategy.DECLARED,
    "adjacent": OrderingStrategy.ADJACENT_SWAP,
    "pairwise": OrderingStrategy.PAIRWISE_SWAP,
    "full": OrderingStrategy.FULL,
}


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)
    inventory = load_inventory()
    pipeline = RulePipeline(
        rules=RULES_FROM_PAPER,
        resolve=lambda label: inventory.segment(label),
    )

    if args.trace:
        wanted = {t.strip() for t in args.trace.split(",") if t.strip()}
        remaining = set(wanted)
        strategy = _STRATEGY_MAP[args.order_strategy]
        for row in _load_rows():
            if row.lemma not in remaining:
                continue
            remaining.discard(row.lemma)
            res = _evaluate(row, pipeline, inventory, args.strict_ipa)
            if args.trace_mode in ("table", "all"):
                print(_render_trace(res, inventory.base_segments()))
            if args.trace_mode in ("explain", "all"):
                print(_explain_trace(res, pipeline))
            if args.trace_mode in ("order", "all"):
                print(_ordering_trace(res, inventory, strategy))
            if args.trace_mode in ("perturb", "all"):
                print(_perturb_trace(res, row.opportunity, row.stem_final))
            if not remaining:
                break
        for lemma in remaining:
            print(f"{lemma}: not found in filtered lexicon rows")
        return

    results: list[RowResult] = []
    for row in _load_rows(limit=args.limit):
        results.append(_evaluate(row, pipeline, inventory, args.strict_ipa))

    report = _render_report(results)
    if args.distance_summary:
        report += _render_distance_summary(results)
    print(report, end="")
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(report, encoding="utf-8")
    print(f"\nReport persisted to {REPORT_PATH}")

    if not args.no_mismatch_csv:
        n = _write_mismatch_csv(results, args.mismatch_csv)
        print(f"Mismatch CSV: {n:,} rows → {args.mismatch_csv}")


if __name__ == "__main__":
    main()
