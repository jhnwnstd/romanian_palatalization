"""Co-variation of the plural class and the verbalizer choice (R23).

Steriade's inflection dependence for Romanian verbs: does a noun whose
plural alternates (colac/colaci) take the palatalizing verbalizer -i
more often than one whose plural avoids or resists the trigger (-uri,
or a non-alternating plural)? The spell-out
TP (:mod:`.spellout`, R20) asks the same question as a tolerance count
over the alternating roots; this stage reports the R23 statistics the
Italian and Czech co-variation stages report, with the same code
(:mod:`lexpipe.stats`).

Units (R23: coded at the base × suffix level)
---------------------------------------------
One per noun root with a derived verb from the competing pair -i vs
-a/-ui (:func:`.spellout.pair_verbs`) whose plural class is known at τ:

    base class P   the plural alternates and is attested at τ
                   (:func:`.spellout.certain`);
    base class N   the attested plural avoids the trigger (-uri) or
                   resists it (:func:`.spellout.prespecified`);
    choice         P (only -i verbs attested at τ), N (only -a/-ui) or
                   both.

Statistics, per target class (dorsal, coronal, z, all) × τ
(dictionary, ≥ 1 token):

  - the 2×2 table (rows base P / N; columns derivative P / N; ``both``
    left out, or counted as P or as N), sample odds ratio with Woolf's
    interval and Fisher's exact p;
  - the Mantel–Haenszel pooled OR over strata = target class × base
    stratum (R22; Romanian strata come from the etymon language only),
    and the seeded stratified permutation of the base class within
    those strata (:func:`lexpipe.stats.stratified_permutation`, seed
    from the run context);
  - the R23 logistic model (:func:`lexpipe.stats.r23_line`) with
    base-clustered errors. The derived verbs carry no etymology, so
    ``deriv_stratum`` is ``unknown`` throughout: the model has no
    ``deriv_stratum`` term, and the report says so ("one level, no
    term").

Outputs: ``analysis/comparison/romanian_covariation.txt`` (tracked;
aggregates only) and ``data/romanian_compare/romanian_compare_covariation_units.csv``
(untracked: the units with their words).
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from collections.abc import Callable, Hashable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Final

from lexpipe import stats
from lexpipe.langs.romanian_compare import morphology as mo
from lexpipe.langs.romanian_compare.lexicon import fold_key
from lexpipe.langs.romanian_compare.records import RootRecord
from lexpipe.langs.romanian_compare.spellout import (
    CLASSES,
    certain,
    pair_verbs,
    prespecified,
)

REPORT_TAUS: Final[tuple[str, ...]] = ("dictionary", "ge1")
BOTH_CODINGS: Final[tuple[str, ...]] = ("excluded", "as P", "as N")
P: Final[str] = "P"
N: Final[str] = "N"
BOTH: Final[str] = "both"
DERIV_STRATUM: Final[str] = "unknown"
_LEXEME: Final[str] = "#lexeme"
UNIT_FIELDS: Final[tuple[str, ...]] = (
    "tau",
    "root_id",
    "lemma",
    "target_class",
    "stratum",
    "base_class",
    "choice",
    "verbs_i",
    "verbs_nonpal",
    "lemma_count",
)
Attest = Callable[[str], Mapping[str, bool]]


@dataclass(frozen=True, slots=True)
class Unit:
    """One root × verbalizer pair at one τ."""

    root_id: str
    lemma: str
    target_class: str
    stratum: str
    base_class: str  # P | N
    choice: str  # P | N | both
    verbs_i: tuple[str, ...]
    verbs_nonpal: tuple[str, ...]
    lemma_count: int


def base_class(record: RootRecord, tau: str, attest: Attest) -> str:
    """``P`` (alternating plural), ``N`` (avoids or resists) or ""."""
    if certain(record, tau, attest):
        return P
    if prespecified(record, tau, attest):
        return N
    return ""


def verb_choice(
    record: RootRecord, tau: str, attest: Attest
) -> tuple[str, tuple[str, ...], tuple[str, ...]]:
    """``(choice, -i verbs, -a/-ui verbs)`` among the pair's verbs
    attested at τ (infinitive + participle tokens); choice "" when
    none is."""
    palatalizing: list[str] = []
    other: list[str] = []
    for verb in pair_verbs(record):
        if not attest(fold_key(verb.form) + _LEXEME).get(tau, False):
            continue
        if verb.suffix in mo.PALATALIZING_VERBALIZERS:
            palatalizing.append(verb.form)
        else:
            other.append(verb.form)
    if palatalizing and other:
        choice = BOTH
    elif palatalizing:
        choice = P
    elif other:
        choice = N
    else:
        choice = ""
    return choice, tuple(sorted(palatalizing)), tuple(sorted(other))


def build_units(
    records: Sequence[RootRecord], tau: str, attest: Attest
) -> list[Unit]:
    """Every root with a known base class and an attested pair verb at
    τ, sorted by root_id."""
    out: list[Unit] = []
    for record in sorted(records, key=lambda r: r.root_id):
        if not record.is_target or record.nde or not pair_verbs(record):
            continue
        klass = base_class(record, tau, attest)
        if not klass:
            continue
        choice, verbs_i, verbs_nonpal = verb_choice(record, tau, attest)
        if not choice:
            continue
        out.append(
            Unit(
                root_id=record.root_id,
                lemma=record.lemma,
                target_class=record.target_class,
                stratum=record.stratum,
                base_class=klass,
                choice=choice,
                verbs_i=verbs_i,
                verbs_nonpal=verbs_nonpal,
                lemma_count=record.lemma_count,
            )
        )
    return out


def observation(unit: Unit, both: str) -> tuple[bool, bool] | None:
    """``(base P, derivative P)``, or None for a ``both`` unit left out."""
    if unit.choice == BOTH:
        if both == "excluded":
            return None
        return unit.base_class == P, both == "as P"
    return unit.base_class == P, unit.choice == P


def _fmt(value: float) -> str:
    if math.isnan(value):
        return "nan"
    if math.isinf(value):
        return "inf"
    return f"{value:.3g}"


def table_line(units: Sequence[Unit], both: str) -> str:
    """The 2×2 table with its odds ratio and Fisher p."""
    items = [obs for u in units if (obs := observation(u, both)) is not None]
    if not items:
        return "no units"
    table = stats.table_from_pairs(items)
    (a, b), (c, d) = table
    ratio = stats.odds_ratio(table)
    text = (
        "OR undefined (empty margin)"
        if ratio is None
        else f"OR {_fmt(ratio.estimate)} [{_fmt(ratio.ci_low)}, "
        f"{_fmt(ratio.ci_high)}], Fisher p = {stats.fisher_exact(table):.3g}"
    )
    return f"P:P {a} P:N {b} N:P {c} N:N {d}; {text}"


def stratified_lines(
    units: Sequence[Unit], both: str, *, seed: int, n_iter: int
) -> list[str]:
    """Mantel–Haenszel over target class × base stratum, and the
    stratified permutation of the base class."""
    rows: list[tuple[Hashable, bool, bool]] = []
    tables: dict[tuple[str, str], list[tuple[bool, bool]]] = defaultdict(list)
    for unit in units:
        obs = observation(unit, both)
        if obs is None:
            continue
        stratum = (unit.target_class, unit.stratum)
        rows.append((stratum, obs[0], obs[1]))
        tables[stratum].append(obs)
    pooled = stats.mantel_haenszel(
        [stats.table_from_pairs(tables[key]) for key in sorted(tables)]
    )
    if pooled is None:
        return ["      Mantel-Haenszel: undefined (no informative stratum)"]
    lines = [
        f"      Mantel-Haenszel OR {_fmt(pooled.estimate)} "
        f"[{_fmt(pooled.ci_low)}, {_fmt(pooled.ci_high)}] over "
        f"{pooled.n_strata} strata (class x base stratum; {pooled.backend})"
    ]
    if math.isfinite(pooled.estimate) and pooled.estimate > 0:
        null = stats.stratified_permutation(
            rows, n_iter=n_iter, seed=seed, observed=pooled.estimate
        )
        lines.append(
            f"      stratified permutation (seed {null.seed}, "
            f"{null.n_valid:,}/{null.n_iter:,} valid): p(OR_MH >= "
            f"observed) = {_fmt(null.p_value)}; null median "
            f"{_fmt(null.median_null)}, 97.5% {_fmt(null.upper_975)}"
        )
    else:
        lines.append("      stratified permutation: skipped (OR_MH 0 or inf)")
    return lines


def r23_units(units: Sequence[Unit]) -> list[stats.R23Unit]:
    """The units of the R23 model (``both`` units left out)."""
    return [
        stats.R23Unit(
            cluster=unit.root_id,
            deriv_pal=unit.choice == P,
            base_class=unit.base_class == P,
            log_lemma_freq=math.log1p(unit.lemma_count),
            base_stratum=unit.stratum,
            deriv_stratum=DERIV_STRATUM,
        )
        for unit in units
        if unit.choice != BOTH
    ]


def _examples(units: Sequence[Unit], test: Callable[[Unit], bool]) -> str:
    chosen = sorted(
        f"{u.lemma} → {'/'.join(u.verbs_i + u.verbs_nonpal)}"
        for u in units
        if test(u)
    )
    shown = "; ".join(chosen[:8])
    more = f"; … (+{len(chosen) - 8})" if len(chosen) > 8 else ""
    return f"[{len(chosen)}] {shown}{more}"


def render_report(
    records: Sequence[RootRecord],
    attest: Attest,
    *,
    seed: int,
    n_iter: int,
) -> str:
    """The text of ``analysis/comparison/romanian_covariation.txt``."""
    lines = [
        "Romanian co-variation: plural class x verbalizer choice (design R23)",
        "=" * 72,
        "",
        "Unit: noun root x the verbalizer pair -i vs -a/-ui. Base class P:",
        "the plural alternates (attested at tau); N: the attested plural",
        "avoids (-uri) or resists the trigger. Choice: P (only -i verbs",
        "attested at tau), N (only -a/-ui), both. 2x2: rows base P / N,",
        "columns choice P / N. OR = sample odds ratio (Woolf interval).",
        f"Permutations: {n_iter:,}, seed {seed}. Statistics: lexpipe.stats,",
        "the code of the Italian and Czech co-variation reports.",
    ]
    for tau in REPORT_TAUS:
        units_all = build_units(records, tau, attest)
        for target, classes in CLASSES.items():
            units = [u for u in units_all if u.target_class in classes]
            lines += ["", f"{target}, tau = {tau}", "-" * 72]
            if not units:
                lines.append("  no units")
                continue
            counts = Counter((u.base_class, u.choice) for u in units)
            lines.append(
                f"  {len(units):,} units: "
                + ", ".join(
                    f"base {b} choice {c} {n}"
                    for (b, c), n in sorted(counts.items())
                )
            )
            for both in BOTH_CODINGS:
                lines.append(f"    both {both}: {table_line(units, both)}")
                lines += stratified_lines(
                    units, both, seed=seed, n_iter=n_iter
                )
            lines.append(stats.r23_line(r23_units(units)))
            lines.append(
                "    base P choosing -a/-ui: "
                + _examples(
                    units, lambda u: u.base_class == P and u.choice == N
                )
            )
            lines.append(
                "    base N choosing -i: "
                + _examples(
                    units, lambda u: u.base_class == N and u.choice == P
                )
            )
    lines += [
        "",
        "Data: noun roots of data/romanian_lexicon_with_freq.csv; derived",
        "verbs from its derived_verbs/deriv_suffixes columns; attestation in",
        "ron_news_2024_1M (lowercased table; aggregates only here, R9).",
    ]
    return "\n".join(lines) + "\n"


def unit_rows(
    records: Sequence[RootRecord], attest: Attest
) -> list[dict[str, Any]]:
    """Every unit at every reported τ (the untracked units file)."""
    rows: list[dict[str, Any]] = []
    for tau in REPORT_TAUS:
        for unit in build_units(records, tau, attest):
            rows.append(
                {
                    "tau": tau,
                    "root_id": unit.root_id,
                    "lemma": unit.lemma,
                    "target_class": unit.target_class,
                    "stratum": unit.stratum,
                    "base_class": unit.base_class,
                    "choice": unit.choice,
                    "verbs_i": list(unit.verbs_i),
                    "verbs_nonpal": list(unit.verbs_nonpal),
                    "lemma_count": unit.lemma_count,
                }
            )
    return rows
