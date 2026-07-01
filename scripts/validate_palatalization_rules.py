#!/usr/bin/env python3
"""Validate the paper's palatalization rules against the full lexicon.

Runs the analysis in ``src/phonology/analyses/romanian_palatalization``
— the paper's revised rule set (place-class S-palatalization, feature-
deleting bleed, broad-terminator dorsal palatalization, multiply-
articulated /j/) — over every noun in
``data/romanian_lexicon_with_freq.csv`` with a target stem-final
consonant. For each row the driver builds an underlying representation
from the lemma's IPA, feeds it through the rule pipeline, and compares
the resulting surface form to the attested plural IPA.

Output
------
Two files under ``analysis/``:

  - ``palatalization_rule_validation.txt`` — the per-cell / cluster /
    mismatch-category report (same shape as before, now with both a
    boolean and a surface-IPA match column).
  - ``palatalization_mismatches.csv`` — every mismatch with predicted
    vs attested IPA and its residual-category tag, for downstream
    inspection.

The framework is pluggable — swap ``RULES_FROM_PAPER`` for any other
tuple of rules in the analysis module and re-run to see how the
lexicon-wide numbers move. Underspecification and inventory patches
are equally pluggable via the same module.

Usage
-----
  python scripts/validate_palatalization_rules.py            # full run
  python scripts/validate_palatalization_rules.py --limit 500  # smoke
  python scripts/validate_palatalization_rules.py --trace muscă,prost
  python scripts/validate_palatalization_rules.py --mismatch-csv out.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Final, Iterator

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PROJECT_ROOT / "src"))

from phonology.analyses.romanian_palatalization import (  # noqa: E402
    RULES_FROM_PAPER,
    load_inventory,
)
from phonology.g2p import build_ur, first_variant  # noqa: E402
from phonology.pipeline import Derivation, RulePipeline  # noqa: E402
from phonology.segments import (  # noqa: E402
    Word,
    segments_to_ipa,
)


CSV_PATH: Final[Path] = (
    _PROJECT_ROOT / "data" / "romanian_lexicon_with_freq.csv"
)
REPORT_PATH: Final[Path] = (
    _PROJECT_ROOT / "analysis" / "palatalization_rule_validation.txt"
)
DEFAULT_MISMATCH_CSV: Final[Path] = (
    _PROJECT_ROOT / "analysis" / "palatalization_mismatches.csv"
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


class MismatchCategory(StrEnum):
    """Residual category a mismatch belongs to.

    Categories separate genuine rule failures from known lexical /
    data-side exceptions. See docstring on
    :func:`categorize_mismatch` for what each one means.
    """

    EZ_ETHNONYM = "EZ_ETHNONYM"
    ICA_SUPPLETION = "ICA_SUPPLETION"
    SHCA_CLUSTER_DATA = "SHCA_CLUSTER_DATA"
    SHT_LOANWORD = "SHT_LOANWORD"
    IPA_SEGMENT_MISMATCH = "IPA_SEGMENT_MISMATCH"
    UNKNOWN = "UNKNOWN"


@dataclass(frozen=True, slots=True)
class NounRow:
    """One lexicon row narrowed to the fields validation needs."""

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
        """True iff the lemma ends in a phonologically ``/sk/`` cluster.

        Excludes ``-șcă`` (with s-cedilla): its singular already has
        surface ``/ʃ/``, so the plural formation is a different beast
        that :meth:`has_shca_cluster` catches instead.
        """
        lem = self.lemma.lower()
        return lem.endswith(("sc", "scă")) and not lem.endswith(
            ("șcă", "şcă")
        )

    @property
    def has_shca_cluster(self) -> bool:
        """True iff the lemma ends in ``-șcă`` (already-palatalized)."""
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
    ur: Word | None                         # None if build_ur failed
    derivation: Derivation | None
    predicted_sr: str                       # IPA, "" if UR failed
    attested_sr: str                        # first variant, cleaned
    boolean_predicted: bool                 # rules → mutation Y/N
    boolean_match: bool                     # prediction vs row.mutation
    ipa_match: bool                         # predicted vs any attested variant
    category: MismatchCategory | None       # only set on mismatch


# ---------------------------------------------------------------------------
# IPA comparison
# ---------------------------------------------------------------------------


def _normalize_glide(ipa: str) -> str:
    """Treat trailing /j/ or /ʲ/ as equivalent to /i/.

    Our pipeline outputs /j/ after glide formation; the lexicon writes
    the same fact as either /ʲ/ (palatalisation on the preceding
    consonant) or the vowel /i/. Normalising to /i/ on both sides
    makes the comparison byte-equal.
    """
    if ipa.endswith(("j", "ʲ")):
        return ipa[:-1] + "i"
    return ipa


def _all_variants(ipa_field: str) -> list[str]:
    """Split a pipe-separated IPA field into individual variants."""
    if not ipa_field:
        return []
    return [v.strip() for v in ipa_field.split(" | ") if v.strip()]


def _ipa_matches(predicted: str, attested_field: str) -> bool:
    """True iff ``predicted`` matches any variant in ``attested_field``.

    Both sides normalised for the trailing-glide convention.
    """
    predicted_norm = _normalize_glide(predicted)
    for variant in _all_variants(attested_field):
        if _normalize_glide(variant) == predicted_norm:
            return True
    return False


# ---------------------------------------------------------------------------
# Derivation → boolean mutation
# ---------------------------------------------------------------------------


_PALATALIZATION_RULES: Final[frozenset[str]] = frozenset({
    "dorsal-pal",
    "s-pal-rev",
    "assibilation",
    "bleed-rev",
})


def _fallback_predict(row: NounRow) -> bool:
    """What the rules WOULD predict for boolean mutation on this row.

    Used when we can't build an actual UR (missing / unparseable IPA).
    Mirrors the paper's rule chain at the level of "does the stem
    alternate?": dorsals palatalize before -i and -e; coronals only
    before -i; clusters (sc/st) palatalize whatever the following
    vowel, except ct blocks before both.
    """
    if row.opportunity in {"uri", "none"}:
        return False
    if row.has_ct_cluster:
        # ct+e: broad terminator halts on /T/, /K/ blocked (paper).
        # ct+i: assibilation on /T/ still fires → mutation.
        return row.opportunity == "i"
    if row.has_sc_cluster:
        return True
    if row.has_st_cluster and row.opportunity == "i":
        return True
    if row.stem_final in DORSAL_STEMS:
        return row.opportunity in {"i", "e"}
    if row.stem_final in CORONAL_STEMS:
        return row.opportunity == "i"
    return False


def _boolean_mutation(deriv: Derivation) -> bool:
    """Did any palatalisation rule fire during the derivation?

    ``mutation`` in the lexicon flags whether the stem-final
    orthographic character alternated between singular and plural.
    For cluster stems the "stem-final" character can be the *first*
    of the trailing cluster (e.g. -st- with stem_final='s'), so
    checking a single UR/SR position misses cluster alternations.
    Any firing of a palatalisation-producing rule signals that
    something in the stem alternated, which is what the boolean
    check is really asking about.
    """
    for step in deriv.steps:
        if step.fired and step.rule in _PALATALIZATION_RULES:
            return True
    return False


# ---------------------------------------------------------------------------
# Mismatch categorisation
# ---------------------------------------------------------------------------


def categorize_mismatch(
    row: NounRow,
    boolean_predicted: bool,
    ipa_match: bool,
) -> MismatchCategory:
    """Assign a residual-category tag to one mismatch row.

    - EZ_ETHNONYM: ``-ez`` agentive/ethnonym class; prespecified /z/
      that lexically selects /-i/ against the spell-out rule.
    - ICA_SUPPLETION: ``-ică`` diminutive with suppletive ``-ele``
      plural (alunică → alunele). Not a phonological alternation.
    - SHCA_CLUSTER_DATA: ``-șcă`` cluster the data underclassifies:
      the rules predict alternation (correctly) but the lexicon's
      ``mutation`` field is False.
    - SHT_LOANWORD: ``-șt`` loanwords like *crenvirșt* where /ʃt/ is
      prespecified in the lemma.
    - IPA_SEGMENT_MISMATCH: boolean prediction agrees with the data,
      but the specific surface IPA does not (a segment-level rule
      failure worth inspecting).
    - UNKNOWN: everything else — a genuine isolated exception.
    """
    lem = row.lemma.lower()
    sf = row.stem_final

    if (
        boolean_predicted == row.mutation
        and not ipa_match
    ):
        return MismatchCategory.IPA_SEGMENT_MISMATCH

    if (
        boolean_predicted
        and not row.mutation
        and sf == "z"
        and row.opportunity == "i"
        and (lem.endswith("ez") or lem.endswith("eză"))
    ):
        return MismatchCategory.EZ_ETHNONYM

    if (
        boolean_predicted
        and not row.mutation
        and sf == "c"
        and row.opportunity == "e"
        and lem.endswith("ică")
        and row.plural.lower().endswith("ele")
    ):
        return MismatchCategory.ICA_SUPPLETION

    if (
        boolean_predicted
        and not row.mutation
        and sf == "c"
        and row.opportunity == "i"
        and lem.endswith(("șcă", "şcă"))
    ):
        return MismatchCategory.SHCA_CLUSTER_DATA

    if (
        boolean_predicted
        and not row.mutation
        and sf == "t"
        and row.opportunity == "i"
        and lem.endswith(("șt", "şt"))
    ):
        return MismatchCategory.SHT_LOANWORD

    return MismatchCategory.UNKNOWN


# ---------------------------------------------------------------------------
# CSV streaming
# ---------------------------------------------------------------------------


def _load_rows(limit: int | None = None) -> Iterator[NounRow]:
    """Stream the lexicon and yield ``NounRow`` for validation.

    Filters (identical to the previous validator so the report is
    comparable): POS=N, stem_final in target set, plural non-empty,
    opportunity != "none", not NDE-blocked.
    """
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


def _evaluate(
    row: NounRow,
    pipeline: RulePipeline,
    inventory,
) -> RowResult:
    """Run the pipeline on one row and score the outcome."""
    ur = build_ur(
        row.ipa_lemma,
        row.opportunity,
        inventory,
        stem_final=row.stem_final,
    )
    attested_first = first_variant(row.ipa_pl) if row.ipa_pl else ""
    if ur is None:
        # IPA missing or unparsable (e.g., foreign-language variants,
        # unknown segments). Fall back to what the rules *would* have
        # predicted based on stem_final × opportunity × cluster —
        # mirroring the paper's asymmetries so the fallback tracks the
        # rules rather than a naive "TARGET_STEMS palatalize always"
        # heuristic.
        fallback_predicted = _fallback_predict(row)
        boolean_match = fallback_predicted == row.mutation
        return RowResult(
            row=row,
            ur=None,
            derivation=None,
            predicted_sr="",
            attested_sr=attested_first,
            boolean_predicted=fallback_predicted,
            boolean_match=boolean_match,
            ipa_match=boolean_match,   # can't check IPA — accept the boolean
            category=(
                None
                if boolean_match
                else categorize_mismatch(row, fallback_predicted, True)
            ),
        )

    deriv = pipeline.derive(ur)
    inv_base = inventory.base_segments()
    predicted_sr = segments_to_ipa(deriv.sr, inv_base)
    boolean_predicted = _boolean_mutation(deriv)
    boolean_match = boolean_predicted == row.mutation
    ipa_match = (
        _ipa_matches(predicted_sr, row.ipa_pl) if row.ipa_pl else boolean_match
    )
    category: MismatchCategory | None = None
    if not (boolean_match and ipa_match):
        category = categorize_mismatch(row, boolean_predicted, ipa_match)
    return RowResult(
        row=row,
        ur=ur,
        derivation=deriv,
        predicted_sr=predicted_sr,
        attested_sr=attested_first,
        boolean_predicted=boolean_predicted,
        boolean_match=boolean_match,
        ipa_match=ipa_match,
        category=category,
    )


# ---------------------------------------------------------------------------
# Report rendering
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


def _render_report(results: list[RowResult]) -> str:
    lines: list[str] = []
    per_cell: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    cluster_counts: dict[str, Counter[str]] = defaultdict(Counter)
    mismatches: dict[MismatchCategory, list[RowResult]] = defaultdict(list)
    bool_match = 0
    ipa_match = 0
    total = 0

    for res in results:
        total += 1
        cell = (res.row.stem_final, res.row.opportunity)
        cell_key = "match" if res.boolean_match else "mismatch"
        per_cell[cell][cell_key] += 1
        per_cell[cell]["ipa_match" if res.ipa_match else "ipa_mismatch"] += 1
        if res.boolean_match:
            bool_match += 1
        if res.ipa_match:
            ipa_match += 1
        if res.category is not None:
            mismatches[res.category].append(res)
        c_key = _cluster_key(res.row)
        if c_key is not None:
            cluster_counts[c_key][
                "match" if res.boolean_match else "mismatch"
            ] += 1

    lines.append(f"Target-stem rows validated: {total:,}")
    if total:
        lines.append(
            f"BOOLEAN mutation  match: {bool_match:,}/{total:,}  "
            f"({100 * bool_match / total:.2f}%)"
        )
        lines.append(
            f"SURFACE-IPA match: {ipa_match:,}/{total:,}  "
            f"({100 * ipa_match / total:.2f}%)"
        )
    lines.append("")

    lines.append("=" * 72)
    lines.append(
        "Per (stem_final × opportunity)  bool_match / mismatch / rate"
    )
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
    lines.append("Mismatch breakdown by residual category")
    lines.append("=" * 72)
    for cat in sorted(
        mismatches.keys(), key=lambda c: -len(mismatches[c])
    ):
        items = mismatches[cat]
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
            if res.predicted_sr and res.attested_sr:
                lines.append(
                    f"      predicted={res.predicted_sr!r:20s} "
                    f"attested={res.attested_sr!r}"
                )
        if len(items) > 6:
            lines.append(f"    ... and {len(items) - 6} more")

    return "\n".join(lines) + "\n"


def _write_mismatch_csv(
    results: list[RowResult],
    path: Path,
) -> int:
    """Write every mismatch row to CSV. Returns count written."""
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "lemma", "plural", "stem_final", "opportunity",
            "ur", "predicted_sr", "attested_sr",
            "category", "cluster",
        ])
        for res in results:
            if res.category is None:
                continue
            if res.boolean_match and res.ipa_match:
                continue
            ur_repr = (
                "".join(s.label for s in res.ur)
                if res.ur is not None
                else ""
            )
            w.writerow([
                res.row.lemma,
                res.row.plural,
                res.row.stem_final,
                res.row.opportunity,
                ur_repr,
                res.predicted_sr,
                res.attested_sr,
                res.category.value,
                _cluster_key(res.row) or "",
            ])
            n += 1
    return n


# ---------------------------------------------------------------------------
# Trace mode
# ---------------------------------------------------------------------------


def _render_trace(
    result: RowResult,
    inv_base,
) -> str:
    """Column-by-column derivation table for one lemma."""
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
    lines.append(f"attested plural IPA: {result.attested_sr!r}")
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
    verdict = "✓" if (result.boolean_match and result.ipa_match) else "✗"
    lines.append(
        f"  match: boolean={result.boolean_match}  "
        f"ipa={result.ipa_match}  {verdict}"
    )
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--limit", type=int, default=None,
        help="Stop after N rows (for quick smoke runs).",
    )
    p.add_argument(
        "--trace", type=str, default="",
        help=(
            "Comma-separated lemmas: print the derivation table for "
            "each and exit."
        ),
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
        for row in _load_rows():
            if row.lemma not in remaining:
                continue
            remaining.discard(row.lemma)
            res = _evaluate(row, pipeline, inventory)
            print(_render_trace(res, inventory.base_segments()))
            if not remaining:
                break
        for lemma in remaining:
            print(f"{lemma}: not found in filtered lexicon rows")
        return

    results: list[RowResult] = []
    for row in _load_rows(limit=args.limit):
        results.append(_evaluate(row, pipeline, inventory))

    report = _render_report(results)
    print(report, end="")
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(report, encoding="utf-8")
    print(f"\nReport persisted to {REPORT_PATH}")

    if not args.no_mismatch_csv:
        n = _write_mismatch_csv(results, args.mismatch_csv)
        print(f"Mismatch CSV: {n:,} rows → {args.mismatch_csv}")


if __name__ == "__main__":
    main()
