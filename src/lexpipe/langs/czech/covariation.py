"""Stage ``covariation``: inflectional class × derivative suffix choice (R23).

Steriade-style co-variation asks whether a base's choice of the
non-palatalizing inflectional allomorph goes with the non-palatalizing
member of a derivational pair (trh: na trhu, trhový but tržní, tržit).

Units (R23: coded at the base × suffix level)
    one per noun root × competing pair (``-ní/-ný`` vs ``-ový/-ovní``,
    ``-it`` vs ``-ovat``, ``-ka`` vs ``-yně``), with the derivative
    choice ``P`` (only the palatalizing member attested at τ), ``N``
    (only the other) or ``both`` (:func:`.spellout.member_side`: the
    -ka/-yně pair counts feminine agents of animate masculine bases
    only, the -it/-ovat pair leaves out bases derived from a verb).
Base class (row variable; attested at τ)
    ``LOC.SG``: m.in and n roots that take -e/-ě (``P``) or only -u
    (``N``); ``NOM.PL``: m.an roots that take -i (``P``) or only
    -ové/-é (``N``). Coding ``any``: a root with both is ``P``; coding
    ``dominant``: the variant with ≥ 90% of the tokens decides
    (:func:`lexpipe.tp.classify_doublet`), true doublets are dropped.
2×2 tables (:mod:`lexpipe.stats` orientation)
    rows: base class P / N; columns: derivative P member / N member.
    ``both`` units are left out (main), or counted as P or as N
    (sensitivity).

Statistics, per base cell × target group × τ: per-pair sample odds
ratios with Woolf intervals and Fisher's exact p; the Mantel–Haenszel
pooled OR over strata = pair × base stratum (R22, Wiktionary); a
seeded permutation of the base class within those strata (seed from
the run context); and, for the velar group, the R23 logistic model of
every language (:func:`lexpipe.stats.r23_line`: deriv_pal ~ base_class
+ log lemma frequency + base stratum + derivative stratum, the strata
treatment-coded, aliased predictors dropped; base-clustered errors),
which needs statsmodels (the report says when it is missing). A unit's
derivative stratum is the Wiktionary stratum its counted derivatives
share, ``mixed`` when they differ.

The non-target baseline (hard labial, l, s, z stems) is a control: the
-ní/-ový choice involves no palatalization there, so an association of
the same size would point to lexical age or learnedness rather than to
the alternation (critique, R23).

Outputs: ``analysis/czech/czech_covariation.txt`` (aggregates) and
``data/czech/nc/czech_covariation_units.csv`` (units).
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from collections.abc import Callable, Hashable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Final

from lexpipe import csvio, stats, tp
from lexpipe.context import RunContext, run_stage
from lexpipe.langs.czech import jsonio, roots
from lexpipe.langs.czech import suffixes as sx
from lexpipe.langs.czech import tpcells as tc
from lexpipe.langs.czech.analysis_layout import AnalysisPaths, analysis_paths
from lexpipe.langs.czech.etymology import STRATA, UNKNOWN
from lexpipe.langs.czech.layout import CzechPaths, paths_for
from lexpipe.langs.czech.spellout import member_side

STAGE: Final[str] = "covariation"
# 3: the shared R23 model (lexpipe.stats.r23_line); gated report writer.
# 4: pair members as in R20 (.spellout.member_side: ka_yne feminine
# agents of m.an bases, it_ovat without deverbal bases); NOM.PL class
# without -a masculines; known-animacy cells (audit 2026-10-23).
VERSION: Final[str] = "4"
SOURCE: Final[tuple[str, ...]] = ("wiktionary_en", "wiktionary_cs")
MIXED: Final[str] = "mixed"
BASIS: Final[str] = "wikt"
ORDER: Final[tuple[str, ...]] = ("linked", "wikt", "wikt+corpus", "all")
REPORT_TAUS: Final[tuple[str, ...]] = ("dictionary", "ge1")
GROUPS: Final[Mapping[str, frozenset[str]]] = {
    "velar": frozenset({"velar"}),
    "r": frozenset({"r"}),
    "coronal_inferred": frozenset({"coronal_inferred"}),
    "baseline": frozenset({"labial", "lateral", "sibilant"}),
}
BASE_CELLS: Final[Mapping[str, tuple[str, ...]]] = {
    "LOC.SG": ("LOC.SG m.in", "LOC.SG n"),
    "NOM.PL": ("NOM.PL m.an",),
}
CODINGS: Final[tuple[str, ...]] = ("any", "dominant")
BOTH_CODINGS: Final[tuple[str, ...]] = ("excluded", "as P", "as N")
UNIT_FIELDS: Final[tuple[str, ...]] = (
    "root",
    "lemma",
    "target_class",
    "base_stratum_wikt",
    "pair",
    "basis",
    *(
        f"{cell}_{coding}_{tau}"
        for cell in BASE_CELLS
        for coding in CODINGS
        for tau in REPORT_TAUS
    ),
    *(f"choice_{tau}" for tau in REPORT_TAUS),
    "p_derivatives",
    "n_derivatives",
    "source",
    "licence_class",
)


# ---------------------------------------------------------------------------
# Classes
# ---------------------------------------------------------------------------


def inflection_class(
    record: roots.RootRecord, cell: str, tau: str, coding: str
) -> str:
    """``P``, ``N`` or "" (no attested slot) for a base cell."""
    counts = {True: 0, False: 0}
    seen: set[bool] = set()
    for name in BASE_CELLS[cell]:
        comp = tc.COMPETITION_BY_NAME[name]
        for obs in tc.competition_observations(record, comp):
            if obs.attested_by_tau.get(tau, False):
                side = bool(obs.alternates)
                seen.add(side)
                counts[side] += obs.count
    if not seen:
        return ""
    if seen == {True}:
        return "P"
    if seen == {False}:
        return "N"
    if coding == "any":
        return "P"
    verdict = tp.classify_doublet({"P": counts[True], "N": counts[False]})
    if verdict.startswith(tp.DOMINANT_PREFIX):
        return verdict[len(tp.DOMINANT_PREFIX) :]
    return ""


@dataclass(frozen=True, slots=True)
class Unit:
    """A base × pair unit at one τ."""

    root: str
    target_class: str
    stratum: str
    pair: str
    base: str
    choice: str
    log_freq: float = 0.0
    deriv_stratum: str = UNKNOWN


def _counted(
    rows: Sequence[Mapping[str, str]], pair: str, tau: str
) -> list[tuple[bool, Mapping[str, str]]]:
    """``(member side, row)`` of the pair's derivatives counted at τ."""
    out: list[tuple[bool, Mapping[str, str]]] = []
    for row in rows:
        if row["suffix_pair"] != pair:
            continue
        side = member_side(row)
        if side is None:
            continue
        if (
            tau != "dictionary"
            and csvio.parse_bool(row[f"attested_{tau}"]) is not True
        ):
            continue
        out.append((side, row))
    return out


def derivative_stratum(
    rows: Sequence[Mapping[str, str]], pair: str, tau: str
) -> str:
    """The Wiktionary stratum the counted derivatives share, ``mixed``
    when they differ, ``unknown`` without one."""
    strata = {
        row["derivative_stratum_wikt"] or UNKNOWN
        for _, row in _counted(rows, pair, tau)
    }
    if len(strata) == 1:
        return next(iter(strata))
    return MIXED if strata else UNKNOWN


def derivative_choice(
    rows: Sequence[Mapping[str, str]], pair: str, tau: str
) -> str:
    """``P``, ``N``, ``both`` or "" for one root and pair at τ."""
    sides = {side for side, _ in _counted(rows, pair, tau)}
    if sides == {True, False}:
        return "both"
    if sides == {True}:
        return "P"
    if sides == {False}:
        return "N"
    return ""


def build_units(
    records: Mapping[str, roots.RootRecord],
    derivatives: Mapping[str, Sequence[Mapping[str, str]]],
    cell: str,
    tau: str,
    coding: str,
) -> list[Unit]:
    """Units with a base class and a derivative choice (sorted)."""
    out: list[Unit] = []
    for key in sorted(derivatives):
        record = records[key]
        if (
            record.target_class in roots.BASELINE_CLASSES
            and not record.in_baseline
        ):
            continue  # soft or adjectival paradigm: not the baseline
        base = inflection_class(record, cell, tau, coding)
        if not base:
            continue
        for pair in sx.PAIRS:
            choice = derivative_choice(derivatives[key], pair, tau)
            if choice:
                out.append(
                    Unit(
                        key,
                        record.target_class,
                        record.stratum or UNKNOWN,
                        pair,
                        base,
                        choice,
                        log_freq=math.log1p(record.lemma_freq),
                        deriv_stratum=derivative_stratum(
                            derivatives[key], pair, tau
                        ),
                    )
                )
    return out


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------


def observation(unit: Unit, both: str) -> tuple[bool, bool] | None:
    """``(base P, derivative P)`` or None (``both`` left out)."""
    if unit.choice == "both":
        if both == "excluded":
            return None
        return unit.base == "P", both == "as P"
    return unit.base == "P", unit.choice == "P"


def table(items: Sequence[tuple[bool, bool]]) -> tuple[tuple[int, int], ...]:
    return stats.table_from_pairs(list(items))


def _fmt(value: float) -> str:
    if math.isnan(value):
        return "nan"
    if math.isinf(value):
        return "inf"
    return f"{value:.2f}"


def pair_lines(units: Sequence[Unit], both: str) -> list[str]:
    lines: list[str] = []
    for pair in sx.PAIRS:
        items = [
            obs
            for u in units
            if u.pair == pair and (obs := observation(u, both)) is not None
        ]
        if not items:
            continue
        tab = table(items)
        (a, b), (c, d) = tab
        ratio = stats.odds_ratio(tab)
        text = (
            "OR undefined (empty margin)"
            if ratio is None
            else f"OR {_fmt(ratio.estimate)} [{_fmt(ratio.ci_low)}, "
            f"{_fmt(ratio.ci_high)}] Fisher p={stats.fisher_exact(tab):.3g}"
        )
        lines.append(
            f"      {sx.PAIR_LABELS[pair]:<24} P:P {a:>4} P:N {b:>4} "
            f"N:P {c:>4} N:N {d:>4}  {text}"
        )
    return lines


def stratified(
    units: Sequence[Unit], both: str, *, seed: int, n_iter: int
) -> list[str]:
    rows: list[tuple[Hashable, bool, bool]] = []
    tables: dict[tuple[str, str], list[tuple[bool, bool]]] = defaultdict(list)
    for unit in units:
        obs = observation(unit, both)
        if obs is None:
            continue
        stratum = (unit.pair, unit.stratum)
        rows.append((stratum, obs[0], obs[1]))
        tables[stratum].append(obs)
    pooled = stats.mantel_haenszel(
        [table(tables[key]) for key in sorted(tables)]
    )
    if pooled is None:
        return ["      Mantel-Haenszel: undefined (no informative stratum)"]
    lines = [
        f"      Mantel-Haenszel OR {_fmt(pooled.estimate)} "
        f"[{_fmt(pooled.ci_low)}, {_fmt(pooled.ci_high)}] over "
        f"{pooled.n_strata} strata (pair x base stratum; {pooled.backend})"
    ]
    if math.isfinite(pooled.estimate) and pooled.estimate > 0:
        null = stats.stratified_permutation(
            rows, n_iter=n_iter, seed=seed, observed=pooled.estimate
        )
        lines.append(
            f"      stratified permutation (seed {null.seed}, {null.n_valid:,}"
            f"/{null.n_iter:,} valid): p(OR_MH >= observed) = "
            f"{_fmt(null.p_value)}; null median {_fmt(null.median_null)}, "
            f"97.5% {_fmt(null.upper_975)}"
        )
    else:
        lines.append("      stratified permutation: skipped (OR_MH 0 or inf)")
    return lines


def r23_units(units: Sequence[Unit]) -> list[stats.R23Unit]:
    """The units of the R23 model (``both`` units left out)."""
    return [
        stats.R23Unit(
            cluster=u.root,
            deriv_pal=u.choice == "P",
            base_class=u.base == "P",
            log_lemma_freq=u.log_freq,
            base_stratum=u.stratum,
            deriv_stratum=u.deriv_stratum,
        )
        for u in units
        if u.choice != "both"
    ]


def logistic_note(units: Sequence[Unit]) -> str:
    """The R23 model (:func:`lexpipe.stats.r23_line`), or why it was
    not fitted."""
    return stats.r23_line(r23_units(units), indent="  ")


# ---------------------------------------------------------------------------
# Stage
# ---------------------------------------------------------------------------


def run(ctx: RunContext) -> int | None:
    """Stage entry point (see the module docstring)."""
    paths = paths_for(ctx)
    apaths = analysis_paths(ctx)
    n_iter = int(ctx.params.get("n_permutations", stats.DEFAULT_PERMUTATIONS))
    return run_stage(
        ctx,
        STAGE,
        inputs=[
            paths.lexicon_complete,
            paths.cells,
            apaths.derivations,
            apaths.lexeme_strata,
        ],
        params={
            "version": VERSION,
            "seed": ctx.seed,
            "n_permutations": n_iter,
            "basis": BASIS,
        },
        outputs=[apaths.covariation, apaths.covariation_units],
        body=lambda: _covariation(ctx, paths, apaths, n_iter),
    )


def load(
    paths: CzechPaths, apaths: AnalysisPaths
) -> tuple[dict[str, roots.RootRecord], dict[str, list[dict[str, str]]]]:
    """Noun records with cells, and pair rows on the basis per root."""
    strata = roots.load_strata(apaths.lexeme_strata)
    records, lexeme_key = roots.load_records(paths.lexicon_complete, strata)
    roots.add_cells(records, lexeme_key, paths.cells)
    limit = ORDER.index(BASIS)
    derivatives: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in jsonio.iter_csv(apaths.derivations):
        key = lexeme_key.get(row["base_lexeme_id"])
        if key is None or not row["suffix_pair"]:
            continue
        if (
            ORDER.index(row["basis"]) <= limit
            and records[key].pos == roots.NOUN
        ):
            derivatives[key].append(row)
    return records, derivatives


def _covariation(
    ctx: RunContext, paths: CzechPaths, apaths: AnalysisPaths, n_iter: int
) -> int:
    records, derivatives = load(paths, apaths)
    text = render(records, derivatives, seed=ctx.seed, n_iter=n_iter)
    ctx.write_tracked_text(apaths.covariation, text, sources=SOURCE)
    csvio.write_rows(
        apaths.covariation_units,
        unit_rows(records, derivatives),
        UNIT_FIELDS,
        ctx=ctx,
    )
    print(text, flush=True)
    return 0


def unit_rows(
    records: Mapping[str, roots.RootRecord],
    derivatives: Mapping[str, Sequence[Mapping[str, str]]],
) -> list[dict[str, Any]]:
    """NC rows: every root × pair with its classes and choices."""
    out: list[dict[str, Any]] = []
    for key in sorted(derivatives):
        record = records[key]
        for pair in sx.PAIRS:
            rows = [r for r in derivatives[key] if r["suffix_pair"] == pair]
            if not rows:
                continue
            item: dict[str, Any] = {
                "root": key,
                "lemma": record.lemma,
                "target_class": record.target_class,
                "base_stratum_wikt": record.stratum,
                "pair": pair,
                "basis": sorted({r["basis"] for r in rows}),
                "p_derivatives": sorted(
                    r["derivative"] for r in rows if member_side(r)
                ),
                "n_derivatives": sorted(
                    r["derivative"] for r in rows if member_side(r) is False
                ),
                "source": ["wiktionary_en", "wiktionary_cs", "leipzig"],
                "licence_class": "by-nc-sa",
            }
            for cell in BASE_CELLS:
                for coding in CODINGS:
                    for tau in REPORT_TAUS:
                        item[f"{cell}_{coding}_{tau}"] = inflection_class(
                            record, cell, tau, coding
                        )
            for tau in REPORT_TAUS:
                item[f"choice_{tau}"] = derivative_choice(rows, pair, tau)
            out.append(item)
    return out


def render(
    records: Mapping[str, roots.RootRecord],
    derivatives: Mapping[str, Sequence[Mapping[str, str]]],
    *,
    seed: int,
    n_iter: int,
) -> str:
    """The tracked report (module docstring)."""
    lines = [
        "Czech co-variation: inflectional class x derivative suffix choice (R23)",
        "=" * 72,
        "",
        "Units: noun root x competing pair; base class P (takes the",
        "palatalizing inflectional allomorph) / N (only the competitor);",
        "derivative choice P / N / both. 2x2: rows base P / N, columns",
        "derivative P member / N member (P:P, P:N, N:P, N:N). Derivative",
        f"basis {BASIS} (Wiktionary links + headword candidates); strata",
        "Wiktionary (R22). OR = sample odds ratio, Woolf 95% interval.",
        f"Permutations: {n_iter:,}, seed {seed}.",
    ]
    for cell in BASE_CELLS:
        for tau in REPORT_TAUS:
            for coding in CODINGS:
                units = build_units(records, derivatives, cell, tau, coding)
                lines += [
                    "",
                    f"{cell} class, tau = {tau}, base coding {coding}",
                    "-" * 72,
                ]
                for group, classes in GROUPS.items():
                    chosen = [u for u in units if u.target_class in classes]
                    if not chosen:
                        continue
                    counts = Counter(u.choice for u in chosen)
                    shares = ", ".join(
                        f"{k} {v}" for k, v in sorted(counts.items())
                    )
                    lines.append(
                        f"  {group}: {len(chosen):,} units ({shares})"
                    )
                    for both in BOTH_CODINGS:
                        if coding == "dominant" and both != "excluded":
                            continue
                        lines.append(f"    both {both}:")
                        lines += pair_lines(chosen, both)
                        lines += stratified(
                            chosen, both, seed=seed, n_iter=n_iter
                        )
                    if group == "velar" and coding == "any":
                        lines.append(logistic_note(chosen))
    lines += _strata_counts(records, derivatives)
    return "\n".join(lines) + "\n"


def _strata_counts(
    records: Mapping[str, roots.RootRecord],
    derivatives: Mapping[str, Sequence[Mapping[str, str]]],
) -> list[str]:
    counts = Counter(
        (records[k].target_class, records[k].stratum) for k in derivatives
    )
    lines = ["", "Roots with a pair derivative, by target class and stratum:"]
    for group, classes in GROUPS.items():
        parts = ", ".join(
            f"{s} {sum(counts[(c, s)] for c in classes)}" for s in STRATA
        )
        lines.append(f"  {group}: {parts}")
    return lines


Check = Callable[[Unit], bool]
