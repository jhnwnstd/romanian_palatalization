"""Stage ``covariation``: base plural × derivative palatalization (R23).

Steriade's inflection dependence, translated to Italian (it_phon.json,
theory mapping): a derivative with an /i/-initial learned suffix
palatalizes iff its base's masculine plural does (amici/amicissimo vs
pochi/pochissimo; tipici/tipicizzare vs antichi/antichizzare).

Units (R23: "the outcome is coded at the base × suffix level")
---------------------------------------------------------------
A unit is (base root, suffix, derivative stratum): the base's
derivatives with one suffix and one stratum (R22) that are attested at
τ (``italian_derivations.csv``: accepted links with a palatal or velar
edge; dictionary-listed at τ = dictionary). Its outcome is ``P`` (only
palatalized derivatives), ``N`` (only velar ones) or ``both``.

The base class is the root's -i status at τ (E0, dominant coding;
:func:`.spellout.root_statuses`): ``P`` = alternates (attested
-ci/-gi), ``N`` = nonalt (attested -chi/-ghi). Roots without an
attested -i plural (the no-evidence roots, R21) are left out here and
analysed by the ``noevidence`` stage.

Statistics (:mod:`lexpipe.stats`)
---------------------------------
- per suffix: the 2×2 table (rows base P / N; columns derivative only P
  / only N; ``both`` units reported, not tabulated), the sample odds
  ratio with Woolf's interval and Fisher's exact p;
- the Mantel–Haenszel pooled OR over suffix × derivative-stratum
  strata, for the internal formations (the design's co-variation set:
  learned borrowings are excluded by stratum, §5) and for all strata;
- the stratified permutation of the base class within suffix ×
  stratum (seeded from the run context, :data:`PERMUTATIONS`
  shuffles);
- the R23 logistic model of every language
  (:func:`lexpipe.stats.r23_logistic`: deriv_pal ~ base_class + log
  lemma frequency + base stratum + derivative stratum, strata
  treatment-coded, aliased predictors dropped; base-clustered errors)
  when statsmodels is importable; otherwise the report says it was
  skipped. Under quasi-complete separation the maximum-likelihood
  estimate does not exist: the report then prints where the fit
  stopped, with NA for every standard error and p-value
  (:func:`logistic_line`).

The -ic- suffix (``ic_suffix``: the base ends in an R26 or unstressed
-ico/-ica, abulico → abulicità) makes plural and derivative palatal
together, since both carry the suffix; every learned-set section is
therefore also given with the -ic- bases left out, and with ``-ic-``
as an extra Mantel–Haenszel stratum.

Every table is reported for the learned set and for the native control
set, at τ = dictionary and ≥ 1 token, for the ``named`` links and with
the generated ``candidate`` links added; sensitivity rows code
``both`` as P or as N and drop derivatives whose named base is a
feminine lexeme (politichese, musicista: the it_phon.json fact-check).

Output: ``analysis/italian/italian_covariation.txt`` (tracked;
aggregates and example words, no per-word counts, R9).
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from lexpipe import stats, tp
from lexpipe.context import RunContext, run_stage
from lexpipe.langs.italian import derivations, etymology, roots, spellout
from lexpipe.langs.italian import suffixes as sx
from lexpipe.langs.italian.layout import Layout

STAGE: Final[str] = "covariation"
# 2: the shared R23 model (lexpipe.stats.r23_line); gated report writer.
# 3: -ic- bases apart; NA standard errors under separation.
# 4: the candidate tier with the palatal-base rejections.
VERSION: Final[int] = 4
PERMUTATIONS: Final[int] = stats.DEFAULT_PERMUTATIONS
REPORT_TAUS: Final[tuple[str, ...]] = ("dictionary", "ge1")
TIERS: Final[Mapping[str, frozenset[str]]] = {
    "named": frozenset({derivations.TIER_NAMED}),
    "named+candidates": frozenset(
        {derivations.TIER_NAMED, derivations.TIER_CANDIDATE}
    ),
    # + the candidates only the palatal-base rule rejected (R16: that
    # rule removes palatal candidates only; spellout.OTHER_BASE_TIER)
    "named+candidates+other_base": frozenset(
        {
            derivations.TIER_NAMED,
            derivations.TIER_CANDIDATE,
            spellout.OTHER_BASE_TIER,
        }
    ),
}
INTERNAL: Final[str] = etymology.INTERNAL
SOURCE: Final[str] = "wiktionary_en"
P: Final[str] = "P"
N: Final[str] = "N"
BOTH: Final[str] = "both"
EXAMPLES: Final[int] = 8


@dataclass(frozen=True, slots=True)
class Unit:
    """One base × suffix × derivative-stratum unit."""

    root_id: str
    lemma: str
    suffix: str
    stratum: str
    base_class: str  # P | N
    outcome: str  # P | N | both
    words: tuple[str, ...]
    feminine_base: bool
    log_freq: float
    base_stratum: str
    ic_suffix: bool = False


def _feminine(item: spellout.Derivative) -> bool:
    return item.base_pos == "noun" and item.base_gender == "f"


def build_units(
    records: roots.Records,
    items: Sequence[spellout.Derivative],
    statuses: Mapping[tuple[str, str], tp.RootStatus],
    tau: str,
    tiers: frozenset[str],
    suffixes: Iterable[str],
) -> list[Unit]:
    """The units at τ (module docstring), sorted."""
    wanted = frozenset(suffixes)
    grouped: dict[tuple[str, str, str], list[spellout.Derivative]] = (
        defaultdict(list)
    )
    for item in items:
        if (
            item.suffix in wanted
            and item.tier in tiers
            and item.palatalized is not None
            and item.attested[tau]
        ):
            grouped[
                (item.base_root_id, item.suffix, item.deriv_stratum)
            ].append(item)
    units: list[Unit] = []
    for (root_id, suffix, stratum), found in sorted(grouped.items()):
        status = statuses.get((root_id, tau))
        if status == tp.RootStatus.ALTERNATES:
            base_class = P
        elif status == tp.RootStatus.NONALT:
            base_class = N
        else:
            continue
        kinds = {item.palatalized for item in found}
        outcome = BOTH if len(kinds) == 2 else (P if True in kinds else N)
        record = records.roots[root_id]
        units.append(
            Unit(
                root_id=root_id,
                lemma=record.lemma,
                suffix=suffix,
                stratum=stratum,
                base_class=base_class,
                outcome=outcome,
                words=tuple(sorted({item.word for item in found})),
                feminine_base=any(_feminine(item) for item in found),
                log_freq=math.log1p(record.frequency("lemma")),
                base_stratum=record.stratum,
                ic_suffix=any(item.ic_suffix for item in found),
            )
        )
    return units


def table_of(
    units: Iterable[Unit], both_as: str = ""
) -> tuple[tuple[int, int], tuple[int, int]]:
    """The 2×2 (base P/N × derivative P/N); ``both`` units are left out
    unless *both_as* codes them."""
    counts: Counter[tuple[bool, bool]] = Counter()
    for unit in units:
        outcome = both_as if unit.outcome == BOTH else unit.outcome
        if outcome not in (P, N):
            continue
        counts[(unit.base_class == P, outcome == P)] += 1
    return (
        (counts[(True, True)], counts[(True, False)]),
        (counts[(False, True)], counts[(False, False)]),
    )


def _fmt_or(result: stats.OddsRatio | None) -> str:
    if result is None:
        return "OR undefined (empty row or column)"
    return (
        f"OR {_num(result.estimate)} "
        f"[{_num(result.ci_low)}, {_num(result.ci_high)}]"
    )


def _num(value: float) -> str:
    if math.isnan(value):
        return "nan"
    if math.isinf(value):
        return "inf"
    return f"{value:.3g}"


def suffix_lines(units: Sequence[Unit], suffixes: Sequence[str]) -> list[str]:
    """Per-suffix 2×2 tables with OR and Fisher p."""
    lines = [
        f"    {'suffix':<10} {'P-base:P/N':>11} {'N-base:P/N':>11} "
        f"{'both':>5}  statistics"
    ]
    for name in suffixes:
        chosen = [u for u in units if u.suffix == name]
        if not chosen:
            continue
        (a, b), (c, d) = table_of(chosen)
        both = sum(1 for u in chosen if u.outcome == BOTH)
        table = ((a, b), (c, d))
        detail = "-"
        if a + b + c + d:
            detail = (
                f"{_fmt_or(stats.odds_ratio(table))}; Fisher p = "
                f"{_num(stats.fisher_exact(table))}"
            )
        lines.append(
            f"    {name:<10} {f'{a}/{b}':>11} {f'{c}/{d}':>11} {both:>5}  "
            f"{detail}"
        )
    return lines


def strata_tables(
    units: Sequence[Unit], both_as: str = "", by_ic: bool = False
) -> list[tuple[tuple[int, int], tuple[int, int]]]:
    """2×2 tables per suffix × stratum (× -ic- base, *by_ic*)."""

    def key(unit: Unit) -> tuple[str, str, bool]:
        return unit.suffix, unit.stratum, by_ic and unit.ic_suffix

    keys = sorted({key(u) for u in units})
    return [
        table_of([u for u in units if key(u) == k], both_as=both_as)
        for k in keys
    ]


def _mh_line(
    label: str, units: Sequence[Unit], both_as: str = "", by_ic: bool = False
) -> str:
    pooled = stats.mantel_haenszel(strata_tables(units, both_as, by_ic))
    if pooled is None:
        return f"    {label}: MH OR undefined"
    return (
        f"    {label}: MH OR {_num(pooled.estimate)} "
        f"[{_num(pooled.ci_low)}, {_num(pooled.ci_high)}] over "
        f"{pooled.n_strata} strata ({pooled.backend})"
    )


def _permutation_line(units: Sequence[Unit], seed: int, n_iter: int) -> str:
    rows = [
        ((u.suffix, u.stratum), u.base_class == P, u.outcome == P)
        for u in units
        if u.outcome != BOTH
    ]
    try:
        result = stats.stratified_permutation(rows, n_iter=n_iter, seed=seed)
    except ValueError as exc:
        return f"    stratified permutation: not computed ({exc})"
    p_value = _num(result.p_value)
    if result.n_valid and result.p_value == 0:
        p_value = f"0 (< {1 / result.n_valid:.2g}: no shuffle reached it)"
    return (
        f"    stratified permutation (base class shuffled within suffix x "
        f"stratum; {result.n_valid}/{result.n_iter} valid shuffles, seed "
        f"{result.seed}): observed MH OR {_num(result.observed)}, "
        f"p = {p_value}, null median "
        f"{_num(result.median_null)}, 97.5% {_num(result.upper_975)}"
    )


def r23_units(units: Sequence[Unit]) -> list[stats.R23Unit]:
    """The units of the R23 model (``both`` units left out)."""
    return [
        stats.R23Unit(
            cluster=unit.root_id,
            deriv_pal=unit.outcome == P,
            base_class=unit.base_class == P,
            log_lemma_freq=unit.log_freq,
            base_stratum=unit.base_stratum,
            deriv_stratum=unit.stratum,
        )
        for unit in units
        if unit.outcome != BOTH
    ]


def logistic_line(units: Sequence[Unit], indent: str = "    ") -> str:
    """The R23 model as report text (the layout of
    :func:`lexpipe.stats.r23_line`), or why it was not fitted. Under
    quasi-complete separation every standard error and p-value is NA:
    the estimates are where the fit stopped, not estimates."""
    found = r23_units(units)
    head = f"{indent}logistic model (R23: {stats.R23_FORMULA}"
    if not found:
        return f"{head}): no units"
    try:
        fit = stats.r23_logistic(found)
    except ImportError as exc:
        return f"{head}): skipped: {exc}"
    except ValueError as exc:
        return f"{head}): not estimable: {exc}"
    _, _, references = stats.r23_design(found)
    separated = bool(fit.separation)

    def term(name: str) -> str:
        if separated:
            return f"{name} {_num(fit.params[name])} (se NA, p NA)"
        return (
            f"{name} {_num(fit.params[name])} (se {_num(fit.bse[name])}, "
            f"p {_num(fit.pvalues[name])})"
        )

    lines = [
        f"{head}; {fit.n_obs} units, {fit.n_clusters} bases, reference "
        f"strata base={references['base_stratum']}, "
        f"deriv={references['deriv_stratum']}):",
        f"{indent}  " + "; ".join(term(name) for name in fit.terms),
    ]
    single = [
        f"{name} (all {references[name]})"
        for name, values in (
            ("base_stratum", {u.base_stratum for u in found}),
            ("deriv_stratum", {u.deriv_stratum for u in found}),
        )
        if len(values) == 1
    ]
    if single:
        lines.append(f"{indent}  one level, no term: " + ", ".join(single))
    if fit.aliased:
        lines.append(
            f"{indent}  aliased (constant or collinear, NA): "
            + ", ".join(fit.aliased)
        )
    if separated:
        lines.append(
            f"{indent}  quasi-complete separation on "
            + ", ".join(fit.separation)
            + ": the maximum-likelihood estimate does not exist; the "
            "values are where the fit stopped, and no standard error or "
            "p-value is reported (NA)"
        )
    if fit.warnings:
        lines.append(f"{indent}  fit warnings: " + "; ".join(fit.warnings))
    return "\n".join(lines)


def _examples(
    units: Sequence[Unit], label: str, test: Callable[[Unit], bool]
) -> str:
    chosen = sorted(
        f"{u.lemma} → {'/'.join(u.words)}" for u in units if test(u)
    )
    shown = "; ".join(chosen[:EXAMPLES])
    more = f"; … (+{len(chosen) - EXAMPLES})" if len(chosen) > EXAMPLES else ""
    return f"    {label} [{len(chosen)}]: {shown}{more}"


def section(
    title: str,
    units: Sequence[Unit],
    suffixes: Sequence[str],
    seed: int,
    n_iter: int,
) -> list[str]:
    """Tables and statistics of one suffix set at one τ and tier."""
    internal = [u for u in units if u.stratum == INTERNAL]
    plain = [u for u in internal if not u.ic_suffix]
    lines = [title, "  " + "-" * 70]
    by_stratum = Counter(u.stratum for u in units)
    lines.append(
        "  units by derivative stratum: "
        + ", ".join(f"{k} {by_stratum[k]}" for k in etymology.STRATA)
    )
    lines.append("  internal formations (the co-variation set, §5):")
    lines += suffix_lines(internal, suffixes)
    ic_units = sum(1 for u in internal if u.ic_suffix)
    lines.append(
        f"  internal formations without -ic- bases ({ic_units} of "
        f"{len(internal)} units have an -ic- base, abulico → abulicità):"
    )
    lines += suffix_lines(plain, suffixes)
    lines.append("  all strata pooled (sensitivity):")
    lines += suffix_lines(units, suffixes)
    for stratum in etymology.STRATA:
        if stratum == INTERNAL:
            continue
        chosen = [u for u in units if u.stratum == stratum]
        if chosen:
            lines.append(f"  stratum {stratum}:")
            lines += suffix_lines(chosen, suffixes)
    lines.append("  pooled:")
    lines.append(_mh_line("internal formations", internal))
    lines.append(_mh_line("internal, -ic- bases dropped", plain))
    lines.append(
        _mh_line("internal, -ic- base as a stratum", internal, by_ic=True)
    )
    lines.append(_mh_line("all strata", units))
    lines.append(_mh_line("internal, both = P", internal, both_as=P))
    lines.append(_mh_line("internal, both = N", internal, both_as=N))
    lines.append(
        _mh_line(
            "internal, feminine named bases dropped",
            [u for u in internal if not u.feminine_base],
        )
    )
    lines.append(_permutation_line(internal, seed, n_iter))
    lines.append(logistic_line(internal))
    lines.append("  examples (internal formations; words only):")
    lines.append(
        _examples(
            internal,
            "P base, velar derivative",
            lambda u: u.base_class == P and u.outcome == N,
        )
    )
    lines.append(
        _examples(
            internal,
            "N base, palatal derivative",
            lambda u: u.base_class == N and u.outcome == P,
        )
    )
    lines.append(
        _examples(
            internal,
            "P base, palatal derivative",
            lambda u: u.base_class == P and u.outcome == P,
        )
    )
    lines.append(
        _examples(
            internal,
            "N base, velar derivative",
            lambda u: u.base_class == N and u.outcome == N,
        )
    )
    lines.append(_examples(internal, "both", lambda u: u.outcome == BOTH))
    lines.append("")
    return lines


FACT_CHECK_WORDS: Final[tuple[str, ...]] = tuple(derivations.FACT_CHECK)


def fact_check_lines(path: Path) -> list[str]:
    """How the it_phon.json fact-check cases are coded in the data."""
    from lexpipe import csvio

    found: dict[str, list[str]] = defaultdict(list)
    for row in csvio.read_rows(path):
        word = row["derivative"].lower()
        if word in FACT_CHECK_WORDS:
            found[word].append(
                f"base {row['base_word']} ({row['base_pos']} "
                f"{row['base_gender'] or '-'}), {row['suffix']}, "
                f"{row['edge']}, {row['link_source']}, {row['link_status']}, "
                f"derivative stratum {row['deriv_stratum']}"
            )
    lines = ["Fact-check cases (it_phon.json verification) as coded:"]
    for word in FACT_CHECK_WORDS:
        rows = found.get(word, [])
        detail = (
            "; ".join(rows)
            if rows
            else "no link (not a derivative of a dorsal base)"
        )
        lines.append(f"  {word}: {detail}")
    return lines


def report(
    records: roots.Records,
    items: Sequence[spellout.Derivative],
    derivations_csv: Path,
    seed: int,
    n_iter: int = PERMUTATIONS,
) -> str:
    """The text report (module docstring)."""
    statuses = spellout.root_statuses(records)
    lines = [
        "Italian co-variation: base masculine plural x derivative "
        "palatalization (design R23)",
        "=" * 72,
        "",
        "Unit: base root x suffix x derivative stratum (outcome P = only",
        "palatalized derivatives, N = only velar, both). Base class: the",
        "root's -i plural at tau (P = attested -ci/-gi, N = attested",
        "-chi/-ghi; E0, dominant coding). Table rows: base P, base N;",
        "columns: derivative P/N. OR = sample odds ratio (Woolf interval),",
        "Fisher two-sided p, Mantel-Haenszel OR over suffix x stratum.",
        "",
    ]
    for tau in REPORT_TAUS:
        for tier, tiers in TIERS.items():
            for group, names in (
                ("learned set", sx.LEARNED_SET),
                ("native control", sx.NATIVE_SET),
            ):
                units = build_units(
                    records, items, statuses, tau, tiers, names
                )
                lines += section(
                    f"tau = {tau}; links: {tier}; {group}",
                    units,
                    names,
                    seed,
                    n_iter,
                )
    lines += fact_check_lines(derivations_csv)
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# The stage
# ---------------------------------------------------------------------------


def _body(ctx: RunContext, layout: Layout, n_iter: int) -> int:
    records = spellout.load_records(layout)
    items = spellout.load_derivatives(layout.derivations)
    text = report(records, items, layout.derivations, ctx.seed, n_iter)
    ctx.write_tracked_text(layout.covariation, text, sources=SOURCE)
    print(f"  co-variation report: {len(text.splitlines()):,} lines")
    return 0


def outputs(layout: Layout) -> list[Path]:
    return [layout.covariation]


def run(ctx: RunContext) -> int | None:
    """Run the ``covariation`` stage (skipped when its manifest matches)."""
    layout = Layout(ctx)
    n_iter = int(ctx.params.get("permutations", PERMUTATIONS))
    return run_stage(
        ctx,
        STAGE,
        inputs=spellout.inputs(layout),
        params={"version": VERSION, "seed": ctx.seed, "permutations": n_iter},
        outputs=outputs(layout),
        body=lambda: _body(ctx, layout, n_iter),
    )
