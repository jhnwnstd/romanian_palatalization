"""Tests for lexpipe.stats, cross-checked with build_contingency_table.py.

Known values use Fisher's tea-tasting table [[3, 1], [1, 3]]: sample
OR 9, two-sided p = 34/70, Woolf 95% interval exp(ln 9 ± z·√(8/3)).
"""

from __future__ import annotations

import importlib.util
import math
import random
import sys
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

pytest.importorskip("scipy")

from lexpipe import stats  # noqa: E402

TEA = [[3, 1], [1, 3]]
Z_975 = 1.959963984540054


@pytest.fixture(scope="module")
def script() -> Iterator[ModuleType]:
    """scripts/build_contingency_table.py, imported by path.

    It is registered in sys.modules while loading (its dataclasses
    resolve string annotations through it) and removed afterwards.
    """
    name = "_lexpipe_test_contingency_script"
    path = PROJECT_ROOT / "scripts" / "build_contingency_table.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
        yield module
    finally:
        sys.modules.pop(name, None)


def _pairs(seed: int, n: int = 120) -> list[tuple[bool, bool]]:
    """Synthetic (row, col) observations with a positive association."""
    rng = random.Random(seed)
    pairs = []
    for _ in range(n):
        row = rng.random() < 0.5
        col = rng.random() < (0.7 if row else 0.3)
        pairs.append((row, col))
    return pairs


def test_odds_ratio_tea_tasting() -> None:
    result = stats.odds_ratio(TEA)
    assert result is not None
    assert result.estimate == 9.0
    half_width = Z_975 * math.sqrt(1 / 3 + 1 + 1 + 1 / 3)
    assert result.ci_low == pytest.approx(math.exp(math.log(9) - half_width))
    assert result.ci_high == pytest.approx(math.exp(math.log(9) + half_width))
    assert result.ci_low == pytest.approx(0.36663693192554553)
    assert result.ci_high == pytest.approx(220.92700692915736)
    assert result.confidence_level == 0.95


def test_fisher_exact_tea_tasting() -> None:
    assert stats.fisher_exact(TEA) == pytest.approx(34 / 70)
    assert stats.fisher_exact(TEA, alternative="greater") == pytest.approx(
        17 / 70
    )


def test_incomplete_table_has_no_odds_ratio() -> None:
    assert stats.odds_ratio([[0, 0], [3, 4]]) is None
    assert stats.odds_ratio([[3, 0], [4, 0]]) is None
    assert not stats.is_complete([[0, 0], [3, 4]])
    assert stats.is_complete([[0, 15], [7, 58]])


def test_zero_cell_keeps_scipy_values() -> None:
    result = stats.odds_ratio([[0, 15], [7, 58]])
    assert result is not None
    assert result.estimate == 0.0
    assert result.ci_low == 0.0
    assert math.isnan(result.ci_high)


@pytest.mark.parametrize(
    "table",
    [[[1, 2, 3], [4, 5, 6]], [[1, 2]], [[1, -2], [3, 4]], [[1.5, 2], [3, 4]]],
)
def test_tables_are_validated(table: list[list[float]]) -> None:
    with pytest.raises(ValueError):
        stats.cells(table)  # type: ignore[arg-type]


def test_matches_contingency_script_panel_stats(script: ModuleType) -> None:
    for a, b, c, d in [(3, 1, 1, 3), (90, 71, 1, 6), (12, 30, 4, 41)]:
        panel = script.PanelStats(
            alt_pal=a, alt_nonpal=b, nonalt_pal=c, nonalt_nonpal=d
        )
        expected = panel.fisher()
        result = stats.odds_ratio([[a, b], [c, d]])
        assert result is not None
        assert (result.estimate, result.ci_low, result.ci_high) == (
            expected[0],
            expected[1],
            expected[2],
        )
        assert stats.fisher_exact([[a, b], [c, d]]) == expected[3]


def test_sample_odds_ratio() -> None:
    assert stats.sample_odds_ratio(3, 1, 1, 3) == 9.0
    assert stats.sample_odds_ratio(3, 0, 1, 3) is None
    assert stats.sample_odds_ratio(0, 1, 1, 3) == 0.0


def test_table_from_pairs() -> None:
    pairs = [(True, True)] * 3 + [(True, False), (False, True)]
    assert stats.table_from_pairs(pairs) == ((3, 1), (1, 0))


def test_permutation_null_is_seeded() -> None:
    pairs = _pairs(1)
    first = stats.permutation_null(pairs, n_iter=500, seed=7)
    again = stats.permutation_null(pairs, n_iter=500, seed=7)
    other = stats.permutation_null(pairs, n_iter=500, seed=8)
    assert first == again
    assert first != other
    observed = stats.odds_ratio(stats.table_from_pairs(pairs))
    assert observed is not None
    assert first.observed == observed.estimate
    assert 0.0 <= first.p_value <= 1.0
    assert first.n_valid == 500 and first.seed == 7
    assert first.median_null <= first.upper_975


def test_permutation_null_matches_contingency_script(
    script: ModuleType,
) -> None:
    pairs = _pairs(2)
    script_pairs = [(row, "pal" if col else "nonpal") for row, col in pairs]
    observed = stats.odds_ratio(stats.table_from_pairs(pairs))
    assert observed is not None
    for seed in (20080715, 3):
        expected = script._permutation_test(
            script_pairs, observed.estimate, n_iter=1000, seed=seed
        )
        result = stats.permutation_null(pairs, n_iter=1000, seed=seed)
        assert (result.p_value, result.median_null, result.upper_975) == (
            expected
        )


def test_permutation_defaults_match_the_script() -> None:
    assert stats.DEFAULT_PERMUTATION_SEED == 20080715
    assert stats.DEFAULT_PERMUTATIONS == 10_000


def test_permutation_null_needs_a_defined_observed_or() -> None:
    one_column = [(True, True), (False, True), (True, True)]
    with pytest.raises(ValueError, match="empty margin"):
        stats.permutation_null(one_column, n_iter=10)
    result = stats.permutation_null(one_column, n_iter=10, observed=1.0)
    assert result.n_valid == 0 and math.isnan(result.p_value)


# ---------------------------------------------------------------------------
# Stratified statistics (R23)
# ---------------------------------------------------------------------------

# R's datasets::UCBAdmissions as Admit x Gender tables, one per
# department. mantelhaen.test(UCBAdmissions) reports a common odds
# ratio of 0.9046968 with 95% interval 0.7719074 - 1.0603298 (the
# Robins-Breslow-Greenland variance).
UCB = [
    [[512, 89], [313, 19]],
    [[353, 17], [207, 8]],
    [[120, 202], [205, 391]],
    [[138, 131], [279, 244]],
    [[53, 94], [138, 299]],
    [[22, 24], [351, 317]],
]


def test_mantel_haenszel_matches_r() -> None:
    result = stats.mantel_haenszel(UCB)
    assert result is not None
    assert result.estimate == pytest.approx(0.9046968, abs=5e-8)
    assert result.ci_low == pytest.approx(0.7719074, abs=5e-8)
    assert result.ci_high == pytest.approx(1.0603298, abs=5e-8)
    assert result.n_strata == 6 and result.confidence_level == 0.95


def _hide_statsmodels(monkeypatch: pytest.MonkeyPatch) -> None:
    """Make every statsmodels import fail, even of submodules another
    test has already imported (they stay in sys.modules otherwise)."""
    names = {"statsmodels", "statsmodels.api"}
    names |= {name for name in sys.modules if name.startswith("statsmodels.")}
    for name in names:
        monkeypatch.setitem(sys.modules, name, None)


def test_mantel_haenszel_builtin_backend(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _hide_statsmodels(monkeypatch)
    result = stats.mantel_haenszel(UCB, confidence_level=0.9)
    assert result is not None and result.backend == "lexpipe"
    assert result.ci_low > 0.7719074 and result.ci_high < 1.0603298


def test_mantel_haenszel_backends_agree() -> None:
    pytest.importorskip("statsmodels")
    external = stats.mantel_haenszel(UCB)
    builtin = stats._mh_builtin(stats._strata(UCB), 0.95)
    assert external is not None and external.backend == "statsmodels"
    assert external.estimate == pytest.approx(builtin.estimate)
    assert external.ci_low == pytest.approx(builtin.ci_low)
    assert external.ci_high == pytest.approx(builtin.ci_high)


def test_one_stratum_is_the_plain_odds_ratio() -> None:
    pooled = stats.mantel_haenszel([TEA])
    single = stats.odds_ratio(TEA)
    assert pooled is not None and single is not None
    assert pooled.estimate == pytest.approx(single.estimate)


def test_mantel_haenszel_edge_cases() -> None:
    assert stats.mantel_haenszel([]) is None
    assert stats.mantel_haenszel([[[0, 0], [0, 0]]]) is None
    assert stats.mantel_haenszel([[[3, 0], [4, 0]]]) is None  # R = S = 0
    empty_dropped = stats.mantel_haenszel([TEA, [[0, 0], [0, 0]]])
    assert empty_dropped is not None and empty_dropped.n_strata == 1
    infinite = stats.mantel_haenszel([[[3, 0], [0, 2]]])
    assert infinite is not None and math.isinf(infinite.estimate)
    assert math.isnan(infinite.ci_low) and math.isnan(infinite.ci_high)
    zero = stats.mantel_haenszel([[[0, 3], [2, 0]]])
    assert zero is not None and zero.estimate == 0.0


def _stratified(seed: int) -> list[tuple[str, bool, bool]]:
    rng = random.Random(seed)
    rows = []
    for stratum, base_rate in (("-ista", 0.2), ("-ità", 0.6), ("-ezza", 0.4)):
        for _ in range(60):
            row = rng.random() < base_rate
            col = rng.random() < (0.75 if row else 0.35)
            rows.append((stratum, row, col))
    return rows


def test_stratified_permutation_is_seeded_and_order_free() -> None:
    rows = _stratified(4)
    first = stats.stratified_permutation(rows, n_iter=300, seed=11)
    again = stats.stratified_permutation(rows[::-1], n_iter=300, seed=11)
    other = stats.stratified_permutation(rows, n_iter=300, seed=12)
    assert first == again
    assert first != other
    pooled = stats.mantel_haenszel(
        [
            stats.table_from_pairs([(r, c) for s, r, c in rows if s == name])
            for name in ("-ista", "-ità", "-ezza")
        ]
    )
    assert pooled is not None
    assert first.observed == pytest.approx(pooled.estimate)
    assert first.observed > 2 and first.p_value < 0.05
    assert first.n_valid == 300 and first.seed == 11


def test_stratified_permutation_keeps_strata_margins() -> None:
    # Within each stratum every row is positive, so no shuffle can
    # change the table: the null is the observed value every time.
    rows = [("a", True, True)] * 3 + [("a", True, False)] * 2
    rows += [("b", False, True)] * 2 + [("b", False, False)] * 4
    with pytest.raises(ValueError, match="undefined"):
        stats.stratified_permutation(rows, n_iter=10)
    result = stats.stratified_permutation(rows, n_iter=10, observed=1.0)
    assert result.n_valid == 0 and math.isnan(result.p_value)


def test_logistic_clustered_names_the_missing_extra(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _hide_statsmodels(monkeypatch)
    rows = [{"y": 1, "x": 0.5, "base": "a"}]
    with pytest.raises(ImportError, match=r"\[stats\]"):
        stats.logistic_clustered(rows, "y", ["x"], "base")


def test_logistic_clustered_fits() -> None:
    pytest.importorskip("statsmodels")
    rng = random.Random(5)
    rows = []
    for base in range(40):
        x = rng.random()
        for _ in range(5):
            p = 1 / (1 + math.exp(-(-1.0 + 2.0 * x)))
            rows.append({"y": int(rng.random() < p), "x": x, "base": base})
    result = stats.logistic_clustered(rows, "y", ["x"], "base")
    assert result.terms == ("Intercept", "x")
    assert result.n_obs == 200 and result.n_clusters == 40
    assert result.params["x"] > 0


def _clustered_rows(n_bases: int = 40) -> list[dict[str, float | int]]:
    rng = random.Random(5)
    rows: list[dict[str, float | int]] = []
    for base in range(n_bases):
        x = rng.random()
        for _ in range(5):
            p = 1 / (1 + math.exp(-(-1.0 + 2.0 * x)))
            rows.append({"y": int(rng.random() < p), "x": x, "base": base})
    return rows


@pytest.mark.parametrize(
    ("rows", "message"),
    [
        ([], "no observations"),
        (  # n = k: statsmodels' small-sample factor divides by n - k
            [{"y": 1, "x": 0.0, "base": "a"}, {"y": 0, "x": 1.0, "base": "b"}],
            "more observations than parameters",
        ),
        (  # one cluster: the cluster factor divides by G - 1
            [{"y": i % 2, "x": i / 4, "base": "a"} for i in range(5)],
            "at least two clusters",
        ),
        (  # a constant predictor is collinear with the intercept
            [{"y": i % 2, "x": 1.0, "base": i} for i in range(6)],
            "rank 1 for 2 terms",
        ),
    ],
)
def test_logistic_clustered_refuses_unestimable_models(
    rows: list[dict[str, float | int | str]], message: str
) -> None:
    pytest.importorskip("statsmodels")
    with pytest.raises(ValueError, match=message):
        stats.logistic_clustered(rows, "y", ["x"], "base")


def test_logistic_clustered_wraps_numerical_failures(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    sm = pytest.importorskip("statsmodels.api")

    class Broken:
        def __init__(self, *args: object, **kwargs: object) -> None:
            pass

        def fit(self, **kwargs: object) -> object:
            raise ZeroDivisionError("float division by zero")

    monkeypatch.setattr(sm, "GLM", Broken)
    with pytest.raises(ValueError, match="ZeroDivisionError"):
        stats.logistic_clustered(_clustered_rows(), "y", ["x"], "base")


def test_logistic_clustered_errors_are_value_errors_or_import_errors(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Without statsmodels every call raises ImportError (the extra);
    the not-estimable checks never mask it."""
    _hide_statsmodels(monkeypatch)
    with pytest.raises(ImportError, match=r"\[stats\]"):
        stats.logistic_clustered([], "y", ["x"], "base")


# ---------------------------------------------------------------------------
# Ties in the stratified permutation; the R23 model
# ---------------------------------------------------------------------------


def test_stratified_permutation_does_not_depend_on_the_observed_backend() -> (
    None
):
    # The observed MH OR from another backend differs from this
    # module's estimate in the last bits; shuffles that reproduce the
    # observed tables must still count as ties (p independent of
    # whether statsmodels is installed).
    rows = _stratified(9)
    exact = stats.stratified_permutation(rows, n_iter=400, seed=3)
    noisy = [
        stats.stratified_permutation(
            rows, n_iter=400, seed=3, observed=exact.observed * factor
        )
        for factor in (1 + 4e-16, 1 - 4e-16, 1 + 1e-12)
    ]
    assert {result.p_value for result in noisy} == {exact.p_value}
    assert stats.tie_threshold(math.inf) == math.inf
    assert stats.tie_threshold(0.0) == 0.0
    assert stats.tie_threshold(2.0) < 2.0


def _r23_units(seed: int, n: int = 160) -> list[stats.R23Unit]:
    rng = random.Random(seed)
    units = []
    for index in range(n):
        base = rng.random() < 0.5
        units.append(
            stats.R23Unit(
                cluster=f"b{index // 2}",
                deriv_pal=rng.random() < (0.8 if base else 0.3),
                base_class=base,
                log_lemma_freq=rng.random() * 6,
                base_stratum=rng.choice(["inherited", "borrowed whole"]),
                deriv_stratum="internal formation",
            )
        )
    return units


def test_r23_design_codes_strata_against_the_most_frequent_level() -> None:
    units = _r23_units(1)
    rows, predictors, references = stats.r23_design(units)
    counts = {
        level: sum(u.base_stratum == level for u in units)
        for level in ("inherited", "borrowed whole")
    }
    reference = max(sorted(counts), key=lambda level: counts[level])
    assert references["base_stratum"] == reference
    other = "inherited" if reference == "borrowed whole" else "borrowed whole"
    assert predictors == [
        "base_class",
        "log_lemma_freq",
        f"base_stratum[{other}]",
    ]
    assert len(rows) == len(units)
    assert {row["cluster"] for row in rows} == {u.cluster for u in units}


def test_logistic_drops_aliased_predictors() -> None:
    pytest.importorskip("statsmodels")
    rng = random.Random(4)
    rows = []
    for index in range(80):
        x = rng.random()
        rows.append(
            {
                "y": 1.0 if rng.random() < 0.3 + 0.4 * x else 0.0,
                "x": x,
                "constant": 1.0,
                "double": 2 * x,
                "base": f"b{index // 2}",
            }
        )
    with pytest.raises(ValueError, match="rank"):
        stats.logistic_clustered(rows, "y", ["x", "constant"], "base")
    fit = stats.logistic_clustered(
        rows, "y", ["x", "constant", "double"], "base", drop_aliased=True
    )
    assert fit.terms == ("Intercept", "x")
    assert fit.aliased == ("constant", "double")


def test_r23_line_reports_coefficients_and_separation() -> None:
    pytest.importorskip("statsmodels")
    text = stats.r23_line(_r23_units(2))
    assert text.startswith("    logistic model (R23: " + stats.R23_FORMULA)
    assert "base_class" in text and "(se " in text
    assert "aliased" not in text  # one deriv stratum: no dummy at all
    assert "one level, no term: deriv_stratum (all internal formation)" in (
        text
    )
    assert "base_stratum (all" not in text
    separated = [
        stats.R23Unit(
            f"b{i}",
            True if i % 2 == 0 else i % 3 == 0,
            i % 2 == 0,
            float(i % 7),
            "inherited",
            "internal formation",
        )
        for i in range(60)
    ]
    assert "quasi-complete separation on base_class" in stats.r23_line(
        separated
    )
    assert "not estimable" in stats.r23_line(separated[:2])
    assert stats.r23_line([]).endswith("no units")


def test_r23_line_without_statsmodels(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _hide_statsmodels(monkeypatch)
    assert "skipped" in stats.r23_line(_r23_units(3))
    assert "statsmodels" in stats.r23_line(_r23_units(3))
