"""2×2 contingency statistics, as computed for the Romanian tables.

The inflection-dependence tables (``analysis/contingency_table.txt``)
come from ``scripts/build_contingency_table.py``. This module exposes
the same computations so the Italian and Czech co-variation tables
are tested identically:

  - :func:`odds_ratio` — ``scipy.stats.contingency.odds_ratio`` with
    ``kind="sample"`` (the plain a·d / b·c) and its 95% interval.
    For the sample kind scipy computes Woolf's normal-approximation
    interval on log OR, not the exact conditional one (the Romanian
    script's docstring says "exact conditional"; the numbers are
    Woolf's). A zero cell gives an estimate of 0 or inf and a NaN
    bound, exactly as in the script.
  - :func:`fisher_exact` — two-sided p of ``scipy.stats.fisher_exact``.
  - :func:`permutation_null` — the script's ``_permutation_test``:
    shuffle the row labels with ``random.Random(seed)`` (in place,
    cumulatively), recompute the sample OR, skip shuffles where it is
    undefined (b or c = 0), and report the share of shuffles with
    OR ≥ observed, the median and the 97.5th-percentile null OR.

Table orientation, as in the script::

                    col positive   col negative
    row positive         a              b
    row negative         c              d

e.g. rows = base plural alternates (K/T͡ʃ) or not (K/K), columns =
derivative takes a palatalizing suffix or not. :func:`odds_ratio`
returns None when a whole row or column is empty (the script's
``is_complete``).

Stratified statistics (R23)
---------------------------
Co-variation pairs are not independent and pooling suffixes invites
Simpson's paradox, so the cross-language tables add:

  - :func:`mantel_haenszel` — the Mantel–Haenszel pooled OR over 2×2
    strata (per suffix × stratum) with the Robins–Breslow–Greenland
    95% interval, the estimator of R's ``mantelhaen.test``. It uses
    statsmodels' ``StratifiedTable`` when statsmodels is installed
    and the same formulas implemented here otherwise; strata with no
    observation are dropped either way.
  - :func:`stratified_permutation` — the seeded permutation null of
    :func:`permutation_null` with the row labels shuffled **within**
    strata and the MH OR as the statistic.
  - :func:`logistic_clustered` — a binomial GLM with cluster-robust
    (base-clustered) standard errors. It needs statsmodels, an
    optional dependency (``pip install -e '.[stats]'``), and says so
    when it is missing (:class:`ImportError`). A model that cannot be
    estimated (too few observations or clusters, a rank-deficient
    design, a numerically failed fit) raises :class:`ValueError`;
    with ``drop_aliased`` constant and collinear predictors are left
    out instead (R's ``NA`` coefficients), and quasi-complete
    separation is reported, not hidden.
  - :func:`r23_logistic` / :func:`r23_line` — the R23 model
    ``deriv_pal ~ base_class + log_lemma_freq + base_stratum +
    deriv_stratum`` (strata treatment-coded against their most
    frequent level; base-clustered errors) from :class:`R23Unit`
    records: the one model every language's co-variation report fits.

The permutation p-value of :func:`stratified_permutation` counts null
statistics within a relative :data:`TIE_TOLERANCE` of the observed one
as ties, so it does not depend on which backend computed the observed
MH OR.
"""

from __future__ import annotations

import importlib
import math
import random
import warnings
from collections import Counter, defaultdict
from collections.abc import Hashable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from statistics import NormalDist
from typing import Any, Final, TypeAlias

from scipy.stats import fisher_exact as _scipy_fisher_exact
from scipy.stats.contingency import odds_ratio as _scipy_odds_ratio

# Seed of the Romanian permutation test (build_contingency_table.py).
DEFAULT_PERMUTATION_SEED: Final[int] = 20080715
DEFAULT_PERMUTATIONS: Final[int] = 10_000

Table: TypeAlias = Sequence[Sequence[int]]


@dataclass(frozen=True, slots=True)
class OddsRatio:
    """Sample odds ratio with its confidence interval."""

    estimate: float
    ci_low: float
    ci_high: float
    confidence_level: float = 0.95


@dataclass(frozen=True, slots=True)
class PermutationNull:
    """Result of :func:`permutation_null`.

    ``p_value`` is the share of valid shuffles with a null OR at least
    the observed one; ``n_valid`` counts the shuffles whose OR was
    defined. All three statistics are NaN if none was.
    """

    observed: float
    p_value: float
    median_null: float
    upper_975: float
    n_valid: int
    n_iter: int
    seed: int


def cells(table: Table) -> tuple[int, int, int, int]:
    """``(a, b, c, d)`` of a 2×2 table of non-negative integer counts."""
    if len(table) != 2 or any(len(row) != 2 for row in table):
        raise ValueError(f"expected a 2x2 table, got {table!r}")
    flat = (table[0][0], table[0][1], table[1][0], table[1][1])
    for value in flat:
        if isinstance(value, bool) or int(value) != value or value < 0:
            raise ValueError(f"counts must be non-negative ints: {table!r}")
    a, b, c, d = (int(v) for v in flat)
    return a, b, c, d


def is_complete(table: Table) -> bool:
    """True when no row or column of the table is empty."""
    a, b, c, d = cells(table)
    return a + b > 0 and c + d > 0 and a + c > 0 and b + d > 0


def odds_ratio(
    table: Table, confidence_level: float = 0.95
) -> OddsRatio | None:
    """Sample OR and interval, or None if a row or column is empty."""
    if not is_complete(table):
        return None
    a, b, c, d = cells(table)
    # A zero cell makes scipy divide by zero inside log(); the values
    # (0, inf, nan) are the intended ones, so silence numpy's warning.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        result = _scipy_odds_ratio([[a, b], [c, d]], kind="sample")
        ci = result.confidence_interval(confidence_level=confidence_level)
    return OddsRatio(
        estimate=float(result.statistic),
        ci_low=float(ci.low),
        ci_high=float(ci.high),
        confidence_level=confidence_level,
    )


def fisher_exact(table: Table, alternative: str = "two-sided") -> float:
    """p-value of Fisher's exact test (two-sided by default)."""
    a, b, c, d = cells(table)
    _, p = _scipy_fisher_exact([[a, b], [c, d]], alternative=alternative)
    return float(p)


def sample_odds_ratio(a: int, b: int, c: int, d: int) -> float | None:
    """a·d / b·c, or None when b or c is 0 (the script's ``_sample_or``)."""
    if b == 0 or c == 0:
        return None
    return (a * d) / (b * c)


def table_from_pairs(
    pairs: Sequence[tuple[bool, bool]],
) -> tuple[tuple[int, int], tuple[int, int]]:
    """Count ``(row_positive, col_positive)`` observations into a 2×2."""
    counts = Counter(pairs)
    return (
        (counts[(True, True)], counts[(True, False)]),
        (counts[(False, True)], counts[(False, False)]),
    )


def permutation_null(
    rows: Sequence[tuple[bool, bool]],
    n_iter: int = DEFAULT_PERMUTATIONS,
    seed: int = DEFAULT_PERMUTATION_SEED,
    *,
    observed: float | None = None,
) -> PermutationNull:
    """Seeded permutation null for the OR of ``(row, col)`` observations.

    *rows* holds one ``(row_positive, col_positive)`` pair per kept
    item. The row labels are shuffled against fixed column labels, so
    both margins are preserved. *observed* defaults to the sample OR
    of :func:`odds_ratio`; the table must then have no empty row or
    column.
    """
    if observed is None:
        result = odds_ratio(table_from_pairs(rows))
        if result is None:
            raise ValueError("observed odds ratio undefined: empty margin")
        observed = result.estimate
    rng = random.Random(seed)
    shuffled = [row for row, _ in rows]
    col_labels = [col for _, col in rows]
    null: list[float] = []
    at_least = 0
    for _ in range(n_iter):
        rng.shuffle(shuffled)
        counts = Counter(zip(shuffled, col_labels, strict=True))
        value = sample_odds_ratio(
            counts[(True, True)],
            counts[(True, False)],
            counts[(False, True)],
            counts[(False, False)],
        )
        if value is None:
            continue
        null.append(value)
        if value >= observed:
            at_least += 1
    if not null:
        nan = float("nan")
        return PermutationNull(observed, nan, nan, nan, 0, n_iter, seed)
    null.sort()
    return PermutationNull(
        observed=observed,
        p_value=at_least / len(null),
        median_null=null[len(null) // 2],
        upper_975=null[int(0.975 * len(null))],
        n_valid=len(null),
        n_iter=n_iter,
        seed=seed,
    )


# ---------------------------------------------------------------------------
# Stratified statistics (R23)
# ---------------------------------------------------------------------------

STATS_EXTRA_HINT: Final[str] = (
    "install the optional statistics dependencies: pip install -e '.[stats]'"
)


@dataclass(frozen=True, slots=True)
class PooledOddsRatio:
    """Mantel–Haenszel pooled OR with its RBG confidence interval.

    ``estimate`` is 0 or inf when every stratum lacks discordant
    products on one side; the bounds are NaN then. ``n_strata``
    counts the strata with at least one observation.
    """

    estimate: float
    ci_low: float
    ci_high: float
    n_strata: int
    confidence_level: float = 0.95
    backend: str = "lexpipe"


def _strata(tables: Sequence[Table]) -> list[tuple[int, int, int, int]]:
    out = [cells(table) for table in tables]
    return [abcd for abcd in out if sum(abcd) > 0]


def _mh_sums(
    strata: Sequence[tuple[int, int, int, int]],
) -> tuple[float, float]:
    """(R, S): sums of a·d/n and b·c/n over the strata."""
    r = sum(a * d / (a + b + c + d) for a, b, c, d in strata)
    s = sum(b * c / (a + b + c + d) for a, b, c, d in strata)
    return r, s


def _rbg_variance(strata: Sequence[tuple[int, int, int, int]]) -> float:
    """Robins–Breslow–Greenland variance of log OR_MH (R and S > 0)."""
    r_total, s_total = _mh_sums(strata)
    pr = ps_qr = qs = 0.0
    for a, b, c, d in strata:
        n = a + b + c + d
        p, q = (a + d) / n, (b + c) / n
        r, s = a * d / n, b * c / n
        pr += p * r
        ps_qr += p * s + q * r
        qs += q * s
    return (
        pr / (2 * r_total**2)
        + ps_qr / (2 * r_total * s_total)
        + qs / (2 * s_total**2)
    )


def _mh_builtin(
    strata: Sequence[tuple[int, int, int, int]], confidence_level: float
) -> PooledOddsRatio:
    r_total, s_total = _mh_sums(strata)
    nan = float("nan")
    if s_total == 0:
        return PooledOddsRatio(
            math.inf, nan, nan, len(strata), confidence_level
        )
    estimate = r_total / s_total
    if r_total == 0:
        return PooledOddsRatio(0.0, nan, nan, len(strata), confidence_level)
    z = NormalDist().inv_cdf(0.5 + confidence_level / 2)
    half = z * math.sqrt(_rbg_variance(strata))
    return PooledOddsRatio(
        estimate,
        estimate * math.exp(-half),
        estimate * math.exp(half),
        len(strata),
        confidence_level,
    )


def _mh_statsmodels(
    strata: Sequence[tuple[int, int, int, int]], confidence_level: float
) -> PooledOddsRatio | None:
    """The same estimate from statsmodels, or None if it is missing."""
    try:
        from statsmodels.stats.contingency_tables import StratifiedTable
    except ImportError:
        return None
    table = StratifiedTable([[[a, b], [c, d]] for a, b, c, d in strata])
    low, high = table.oddsratio_pooled_confint(alpha=1 - confidence_level)
    return PooledOddsRatio(
        float(table.oddsratio_pooled),
        float(low),
        float(high),
        len(strata),
        confidence_level,
        backend="statsmodels",
    )


def mantel_haenszel(
    tables: Sequence[Table], confidence_level: float = 0.95
) -> PooledOddsRatio | None:
    """Pooled OR over 2×2 *tables* (same orientation as :func:`cells`).

    None when no table has an observation, or no stratum has a
    defined OR on either side (R = S = 0).
    """
    strata = _strata(tables)
    if not strata:
        return None
    r_total, s_total = _mh_sums(strata)
    if r_total == 0 and s_total == 0:
        return None
    if r_total > 0 and s_total > 0:
        external = _mh_statsmodels(strata, confidence_level)
        if external is not None:
            return external
    return _mh_builtin(strata, confidence_level)


def _table_cells(
    counts: Iterable[tuple[tuple[bool, bool], int]],
) -> tuple[int, int, int, int]:
    """(a, b, c, d) from ``((row, col), count)`` items."""
    table = dict(counts)
    return (
        table.get((True, True), 0),
        table.get((True, False), 0),
        table.get((False, True), 0),
        table.get((False, False), 0),
    )


def _mh_estimate(strata: Sequence[tuple[int, int, int, int]]) -> float | None:
    r_total, s_total = _mh_sums(strata)
    if s_total == 0:
        return None if r_total == 0 else math.inf
    return r_total / s_total


# Relative tolerance of "null OR >= observed OR" in
# stratified_permutation: MH estimates of the same tables computed with
# different floating-point operation orders differ in the last bits.
TIE_TOLERANCE: Final[float] = 1e-9


def tie_threshold(observed: float) -> float:
    """The smallest null statistic counted as "at least *observed*":
    *observed* minus :data:`TIE_TOLERANCE` of its size (inf and 0 are
    returned unchanged)."""
    if math.isinf(observed) or observed == 0:
        return observed
    return observed - TIE_TOLERANCE * abs(observed)


def stratified_permutation(
    rows: Sequence[tuple[Hashable, bool, bool]],
    n_iter: int = DEFAULT_PERMUTATIONS,
    seed: int = DEFAULT_PERMUTATION_SEED,
    *,
    observed: float | None = None,
) -> PermutationNull:
    """Seeded permutation null of the MH OR, shuffling within strata.

    *rows* holds one ``(stratum, row_positive, col_positive)`` triple
    per item (e.g. stratum = (suffix, base stratum)). In every
    iteration the row labels are shuffled within each stratum (in
    place and cumulatively, as in :func:`permutation_null`), keeping
    every stratum's margins; shuffles with an undefined MH OR are
    skipped. Strata are visited in ``repr`` order and their items in a
    canonical order, so the result depends on the seed and the
    multiset of rows only. *observed* defaults to the MH OR of the
    data, which must then be defined.

    A null OR counts as "at least the observed" when it reaches
    :func:`tie_threshold` of *observed*: a shuffle that reproduces the
    observed tables gives the same OR up to floating-point rounding,
    whichever backend computed *observed* (statsmodels' or this
    module's :func:`mantel_haenszel`), so the p-value does not depend
    on whether statsmodels is installed.
    """
    groups: dict[Hashable, list[tuple[bool, bool]]] = defaultdict(list)
    for stratum, row, col in rows:
        groups[stratum].append((bool(row), bool(col)))
    ordered = [sorted(groups[key]) for key in sorted(groups, key=repr)]
    if observed is None:
        counted = [_table_cells(Counter(pairs).items()) for pairs in ordered]
        value = _mh_estimate([abcd for abcd in counted if sum(abcd) > 0])
        if value is None:
            raise ValueError("observed MH odds ratio undefined")
        observed = value
    row_labels = [[row for row, _ in pairs] for pairs in ordered]
    col_labels = [[col for _, col in pairs] for pairs in ordered]
    threshold = tie_threshold(observed)
    rng = random.Random(seed)
    null: list[float] = []
    at_least = 0
    for _ in range(n_iter):
        strata: list[tuple[int, int, int, int]] = []
        for labels, cols in zip(row_labels, col_labels, strict=True):
            rng.shuffle(labels)
            strata.append(
                _table_cells(Counter(zip(labels, cols, strict=True)).items())
            )
        value = _mh_estimate(strata)
        if value is None:
            continue
        null.append(value)
        if value >= threshold:
            at_least += 1
    if not null:
        nan = float("nan")
        return PermutationNull(observed, nan, nan, nan, 0, n_iter, seed)
    null.sort()
    return PermutationNull(
        observed=observed,
        p_value=at_least / len(null),
        median_null=null[len(null) // 2],
        upper_975=null[int(0.975 * len(null))],
        n_valid=len(null),
        n_iter=n_iter,
        seed=seed,
    )


@dataclass(frozen=True, slots=True)
class ClusteredLogit:
    """Coefficients of :func:`logistic_clustered`, keyed by term.

    ``aliased`` lists the predictors left out because they are constant
    or a linear combination of earlier columns (R's ``NA``
    coefficients; only with ``drop_aliased``); ``separation`` the 0/1
    predictors on one of whose levels the outcome never varies
    (quasi-complete separation: the maximum-likelihood estimate of that
    coefficient does not exist, and the reported value is where the
    fit stopped); ``warnings`` the fit's warnings, in order.
    """

    terms: tuple[str, ...]
    params: Mapping[str, float]
    bse: Mapping[str, float]
    pvalues: Mapping[str, float]
    n_obs: int
    n_clusters: int
    aliased: tuple[str, ...] = ()
    separation: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()


def _not_estimable() -> tuple[type[BaseException], ...]:
    """Exceptions of a statsmodels fit that mean "not estimable"
    (numpy's LinAlgError is a ValueError)."""
    errors: tuple[type[BaseException], ...] = (ArithmeticError, ValueError)
    try:
        module = importlib.import_module("statsmodels.tools.sm_exceptions")
    except ImportError:
        return errors
    separation = getattr(module, "PerfectSeparationError", None)
    if isinstance(separation, type) and issubclass(separation, BaseException):
        errors += (separation,)
    return errors


def _separated(
    rows: Sequence[Mapping[str, Any]], outcome: str, predictor: str
) -> bool:
    """Is *predictor* a 0/1 column with a level on which *outcome* is
    constant (quasi-complete separation)?"""
    levels: dict[float, set[float]] = defaultdict(set)
    for row in rows:
        levels[float(row[predictor])].add(float(row[outcome]))
    if set(levels) != {0.0, 1.0}:
        return False
    return any(len(values) == 1 for values in levels.values())


def _independent_columns(
    np: Any, matrix: Any, names: Sequence[str], fixed: int
) -> list[int]:
    """Indices of the columns kept when each column, in order, is kept
    only if it raises the rank (the first *fixed* are always kept)."""
    kept: list[int] = list(range(fixed))
    rank = int(np.linalg.matrix_rank(matrix[:, kept])) if kept else 0
    for index in range(fixed, len(names)):
        trial = kept + [index]
        trial_rank = int(np.linalg.matrix_rank(matrix[:, trial]))
        if trial_rank > rank:
            kept, rank = trial, trial_rank
    return kept


def logistic_clustered(
    rows: Sequence[Mapping[str, Any]],
    outcome: str,
    predictors: Sequence[str],
    cluster: str,
    *,
    intercept: bool = True,
    drop_aliased: bool = False,
) -> ClusteredLogit:
    """Binomial GLM of *outcome* on *predictors*, clustered by *cluster*.

    Columns must be numeric (0/1 outcome; dummy-code categorical
    predictors such as base_class or stratum yourself). Standard
    errors are statsmodels' ``cov_type="cluster"`` with the cluster
    labels as groups (e.g. the base, R23). Raises :class:`ImportError`
    naming the optional extra when statsmodels is not installed.

    With *drop_aliased*, a predictor that is constant or a linear
    combination of the intercept and earlier predictors is left out
    and listed in ``aliased``, as R's ``glm`` reports an ``NA``
    coefficient for it; without it such a design raises.

    Raises :class:`ValueError` when the model cannot be estimated, so
    callers handle every such case alike: no observations; no more
    observations than parameters (the small-sample factor of the
    cluster-robust covariance divides by n − k); fewer than two
    clusters; a rank-deficient design (a constant or collinear
    predictor, unless *drop_aliased*); or a fit that fails numerically
    (statsmodels' singular matrices, divisions by zero, perfect
    separation errors). Quasi-complete separation by a 0/1 predictor
    does not raise: statsmodels warns and returns where the fit
    stopped, and the predictor is listed in ``separation``.
    """
    try:
        # Imported by name: numpy's own stubs need a newer Python than
        # the 3.11 the type checker targets, and statsmodels is optional.
        np = importlib.import_module("numpy")
        sm = importlib.import_module("statsmodels.api")
    except ImportError as exc:
        raise ImportError(
            f"logistic_clustered needs statsmodels; {STATS_EXTRA_HINT}"
        ) from exc
    if not rows:
        raise ValueError("no observations")
    names = (("Intercept",) if intercept else ()) + tuple(predictors)
    full = np.asarray(
        [
            ([1.0] if intercept else []) + [float(row[p]) for p in predictors]
            for row in rows
        ],
        dtype=float,
    )
    aliased: tuple[str, ...] = ()
    if drop_aliased:
        kept = _independent_columns(np, full, names, 1 if intercept else 0)
        aliased = tuple(n for i, n in enumerate(names) if i not in kept)
        terms = tuple(names[i] for i in kept)
        matrix = full[:, kept]
    else:
        terms = names
        matrix = full
    labels = sorted({row[cluster] for row in rows}, key=repr)
    if len(rows) <= len(terms):
        raise ValueError(
            f"{len(rows)} observation(s) for {len(terms)} parameter(s): "
            "the cluster-robust covariance needs more observations than "
            "parameters"
        )
    if len(labels) < 2:
        raise ValueError(
            "cluster-robust standard errors need at least two clusters, "
            f"got {len(labels)}"
        )
    rank = int(np.linalg.matrix_rank(matrix))
    if rank < len(terms):
        raise ValueError(
            f"the design matrix has rank {rank} for {len(terms)} terms "
            f"{list(terms)}: a predictor is constant or collinear"
        )
    index = {label: i for i, label in enumerate(labels)}
    groups = np.asarray([index[row[cluster]] for row in rows])
    y = np.asarray([float(row[outcome]) for row in rows])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            fit = sm.GLM(y, matrix, family=sm.families.Binomial()).fit(
                cov_type="cluster", cov_kwds={"groups": groups}
            )
        except _not_estimable() as exc:
            raise ValueError(
                f"the fit failed ({type(exc).__name__}: {exc})"
            ) from exc
    messages = tuple(
        dict.fromkeys(
            f"{w.category.__name__}: {str(w.message).splitlines()[0]}"
            for w in caught
        )
    )
    separation = tuple(
        name
        for name in terms
        if name != "Intercept" and _separated(rows, outcome, name)
    )
    return ClusteredLogit(
        terms=terms,
        params=dict(zip(terms, map(float, fit.params), strict=True)),
        bse=dict(zip(terms, map(float, fit.bse), strict=True)),
        pvalues=dict(zip(terms, map(float, fit.pvalues), strict=True)),
        n_obs=len(rows),
        n_clusters=len(labels),
        aliased=aliased,
        separation=separation,
        warnings=messages,
    )


# ---------------------------------------------------------------------------
# The R23 model, the same for every language
# ---------------------------------------------------------------------------

R23_FORMULA: Final[str] = (
    "deriv_pal ~ base_class + log_lemma_freq + base_stratum + deriv_stratum"
)


@dataclass(frozen=True, slots=True)
class R23Unit:
    """One base × suffix unit of the R23 model.

    ``cluster`` is the base (clustered errors); ``deriv_pal`` whether
    the unit's derivatives palatalize (units with both outcomes are
    left out by the caller); ``base_class`` whether the base takes the
    palatalizing inflectional allomorph; ``log_lemma_freq`` the log of
    the base's lemma frequency (log(1 + count)); the strata are R22
    labels; ``stress`` an optional 0/1 stress predictor.
    """

    cluster: str
    deriv_pal: bool
    base_class: bool
    log_lemma_freq: float
    base_stratum: str
    deriv_stratum: str
    stress: float | None = None


def _reference(levels: Sequence[str]) -> str:
    """The reference level: the most frequent, ties by name."""
    counts = Counter(levels)
    return min(counts, key=lambda level: (-counts[level], level))


def r23_design(
    units: Sequence[R23Unit],
) -> tuple[list[dict[str, Any]], list[str], dict[str, str]]:
    """Rows, predictor names and reference levels of the R23 model.

    Strata are treatment-coded against their most frequent level
    (``base_stratum[borrowed whole]`` = 1 for a borrowed base); a
    ``stress`` predictor is added when every unit has one.
    """
    references = {
        "base_stratum": _reference([u.base_stratum for u in units]),
        "deriv_stratum": _reference([u.deriv_stratum for u in units]),
    }
    levels = {
        "base_stratum": sorted({u.base_stratum for u in units}),
        "deriv_stratum": sorted({u.deriv_stratum for u in units}),
    }
    stress = bool(units) and all(u.stress is not None for u in units)
    predictors = ["base_class", "log_lemma_freq"]
    for name in ("base_stratum", "deriv_stratum"):
        predictors += [
            f"{name}[{level}]"
            for level in levels[name]
            if level != references[name]
        ]
    if stress:
        predictors.append("stress")
    rows: list[dict[str, Any]] = []
    for unit in units:
        row: dict[str, Any] = {
            "deriv_pal": 1.0 if unit.deriv_pal else 0.0,
            "base_class": 1.0 if unit.base_class else 0.0,
            "log_lemma_freq": unit.log_lemma_freq,
            "cluster": unit.cluster,
        }
        for name, value in (
            ("base_stratum", unit.base_stratum),
            ("deriv_stratum", unit.deriv_stratum),
        ):
            for level in levels[name]:
                if level != references[name]:
                    row[f"{name}[{level}]"] = 1.0 if value == level else 0.0
        if stress:
            row["stress"] = float(unit.stress or 0.0)
        rows.append(row)
    return rows, predictors, references


def r23_logistic(units: Sequence[R23Unit]) -> ClusteredLogit:
    """The R23 model (:data:`R23_FORMULA`) with base-clustered errors,
    aliased predictors dropped (:func:`logistic_clustered`)."""
    rows, predictors, _ = r23_design(units)
    return logistic_clustered(
        rows, "deriv_pal", predictors, "cluster", drop_aliased=True
    )


def _g(value: float) -> str:
    return f"{value:.3g}"


def r23_line(units: Sequence[R23Unit], indent: str = "    ") -> str:
    """The R23 model as report text: every coefficient with its
    clustered SE and p, the stratum factors with a single level (they
    have no term), the aliased predictors and any separation; or why
    the model was not fitted."""
    head = f"{indent}logistic model (R23: {R23_FORMULA}"
    if not units:
        return f"{head}): no units"
    try:
        fit = r23_logistic(units)
    except ImportError as exc:
        return f"{head}): skipped: {exc}"
    except ValueError as exc:
        return f"{head}): not estimable: {exc}"
    _, _, references = r23_design(units)
    terms = "; ".join(
        f"{term} {_g(fit.params[term])} (se {_g(fit.bse[term])}, "
        f"p {_g(fit.pvalues[term])})"
        for term in fit.terms
    )
    lines = [
        f"{head}; {fit.n_obs} units, {fit.n_clusters} bases, reference "
        f"strata base={references['base_stratum']}, "
        f"deriv={references['deriv_stratum']}):",
        f"{indent}  {terms}",
    ]
    single = [
        f"{name} (all {references[name]})"
        for name, values in (
            ("base_stratum", {u.base_stratum for u in units}),
            ("deriv_stratum", {u.deriv_stratum for u in units}),
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
    if fit.separation:
        lines.append(
            f"{indent}  quasi-complete separation on "
            + ", ".join(fit.separation)
            + ": the maximum-likelihood estimate does not exist; the "
            "value is where the fit stopped"
        )
    if fit.warnings:
        lines.append(f"{indent}  fit warnings: " + "; ".join(fit.warnings))
    return "\n".join(lines)
