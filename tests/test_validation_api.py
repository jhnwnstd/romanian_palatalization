#!/usr/bin/env python3
"""End-to-end tests of the CLI-free Python API.

Every function callable by an end user is exercised here — the tests
double as usage documentation. If a caller opens a REPL, they'd write
code like these tests.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from phonology.analyses.romanian_palatalization import (  # noqa: E402
    AnalysisProfile,
    build_profile,
    expected_firings,
    fallback_predict,
)
from phonology.diagnostics.ordering import (  # noqa: E402
    OrderingStrategy,
    search_orderings,
)
from phonology.diagnostics.perturb import (  # noqa: E402
    PerturbationKind,
    search_perturbations_for,
)
from phonology.g2p import build_ur  # noqa: E402
from phonology.lexicon import (  # noqa: E402
    ClusterTag,
    LexRow,
    iter_lexicon_rows,
)
from phonology.validation import (  # noqa: E402
    MatchMode,
    RowResult,
    TraceResult,
    ValidationReport,
    evaluate_row,
    render_report,
    render_trace,
    run_validation,
    trace_lemmas,
    write_mismatch_csv,
)


CSV_PATH = PROJECT_ROOT / "data" / "romanian_lexicon_with_freq.csv"


@pytest.fixture(scope="module")
def profile() -> AnalysisProfile:
    return build_profile()


@pytest.fixture(scope="module")
def sample_rows(profile) -> tuple[LexRow, ...]:
    if not CSV_PATH.exists():
        pytest.skip("lexicon CSV missing")
    return tuple(iter_lexicon_rows(
        CSV_PATH, target_stems=profile.target_stems, limit=200,
    ))


# ---------------------------------------------------------------------------
# AnalysisProfile: everything a caller needs, bundled
# ---------------------------------------------------------------------------


class TestProfile:
    def test_profile_carries_all_needed_pieces(self, profile) -> None:
        assert profile.name == "romanian-palatalization"
        assert profile.inventory is not None
        assert len(profile.rules) == 9
        # /z/ included for lexicon filtering even though it can't
        # trigger palatalization — z-final rows should predict False
        # and match the data's False.
        assert profile.target_stems == frozenset("cgtdsz")
        assert "dorsal-pal" in profile.palatalization_rule_names
        assert callable(profile.expected_firings)
        assert callable(profile.fallback_predict)

    def test_immutability(self, profile) -> None:
        with pytest.raises(Exception):
            profile.name = "hack"


# ---------------------------------------------------------------------------
# Oracle: expected_firings ↔ fallback_predict
# ---------------------------------------------------------------------------


class TestOracle:
    """The oracle is one function; fallback_predict is a view of it."""

    def _row(
        self, lemma: str, stem_final: str, opp: str, mutation: bool,
    ) -> LexRow:
        from phonology.lexicon import _tag_cluster
        return LexRow(
            lemma=lemma, plural=lemma + "i", pos="N",
            stem_final=stem_final, opportunity=opp,
            mutation=mutation, ipa_lemma="", ipa_pl="",
            exception_reason="", cluster=_tag_cluster(lemma),
        )

    def test_dorsal_i_expects_dorsal_pal(self) -> None:
        row = self._row("kolak", "c", "i", True)
        firings = expected_firings(row)
        assert "dorsal-pal" in firings

    def test_uri_expects_nothing(self) -> None:
        row = self._row("fok", "c", "uri", False)
        assert expected_firings(row) == frozenset()

    def test_ct_e_blocked(self) -> None:
        row = self._row("punkt", "t", "e", False)
        # /K/ blocked by broad terminator; /T/ can't assibilate
        # before /e/ (not postalveolar). Nothing should fire.
        assert expected_firings(row) == frozenset()

    def test_shca_is_reachable(self) -> None:
        row = self._row("bumașcă", "c", "i", False)
        assert row.cluster is ClusterTag.SHCA
        assert expected_firings(row)   # non-empty

    def test_fallback_matches_expected_boolean(self) -> None:
        row = self._row("brad", "d", "i", True)
        assert fallback_predict(row) is bool(expected_firings(row))


# ---------------------------------------------------------------------------
# LexRow + iter_lexicon_rows
# ---------------------------------------------------------------------------


class TestLexicon:
    def test_iter_yields_lexrow(self, profile) -> None:
        if not CSV_PATH.exists():
            pytest.skip("lexicon CSV missing")
        one = next(iter(iter_lexicon_rows(
            CSV_PATH, target_stems=profile.target_stems, limit=1,
        )))
        assert isinstance(one, LexRow)
        assert one.stem_final in profile.target_stems

    def test_limit_honored(self, sample_rows) -> None:
        assert len(sample_rows) == 200

    def test_cluster_detected_once_at_construction(self) -> None:
        from phonology.lexicon import _tag_cluster
        assert _tag_cluster("muscă") is ClusterTag.SC
        assert _tag_cluster("bumașcă") is ClusterTag.SHCA
        assert _tag_cluster("punct") is ClusterTag.CT
        assert _tag_cluster("prost") is ClusterTag.ST
        assert _tag_cluster("rus") is ClusterTag.NONE


# ---------------------------------------------------------------------------
# evaluate_row + run_validation
# ---------------------------------------------------------------------------


class TestEvaluateAndRun:
    def test_evaluate_returns_rowresult(self, profile, sample_rows) -> None:
        res = evaluate_row(sample_rows[0], profile)
        assert isinstance(res, RowResult)
        assert res.row is sample_rows[0]

    def test_strict_mode_narrower(self, profile, sample_rows) -> None:
        # Strict should never match MORE than normalised on the same row.
        norm = sum(
            1 for r in sample_rows
            if evaluate_row(r, profile, match_mode=MatchMode.NORMALISED).ipa_match
        )
        strict = sum(
            1 for r in sample_rows
            if evaluate_row(r, profile, match_mode=MatchMode.STRICT).ipa_match
        )
        assert strict <= norm

    def test_run_validation_returns_immutable_report(
        self, profile, sample_rows,
    ) -> None:
        report = run_validation(profile, sample_rows)
        assert isinstance(report, ValidationReport)
        assert report.total == len(sample_rows)
        assert isinstance(report.results, tuple)
        with pytest.raises(Exception):
            report.total = 0


# ---------------------------------------------------------------------------
# trace_lemmas
# ---------------------------------------------------------------------------


class TestTraceLemmas:
    def test_returns_one_per_input_lemma(self, profile) -> None:
        if not CSV_PATH.exists():
            pytest.skip("lexicon CSV missing")
        traces = trace_lemmas(
            ["muscă", "prost", "not-a-real-lemma"],
            profile, CSV_PATH,
        )
        assert len(traces) == 3
        assert isinstance(traces[0], TraceResult)
        assert traces[0].found is True
        assert traces[1].found is True
        assert traces[2].found is False

    def test_include_explanation(self, profile) -> None:
        if not CSV_PATH.exists():
            pytest.skip("lexicon CSV missing")
        traces = trace_lemmas(
            ["muscă"], profile, CSV_PATH, include_explanation=True,
        )
        assert traces[0].explanation is not None

    def test_include_ordering(self, profile) -> None:
        if not CSV_PATH.exists():
            pytest.skip("lexicon CSV missing")
        traces = trace_lemmas(
            ["muscă"], profile, CSV_PATH,
            include_ordering=OrderingStrategy.ADJACENT_SWAP,
        )
        assert traces[0].ordering is not None
        assert traces[0].ordering.baseline_matched


# ---------------------------------------------------------------------------
# search_orderings + search_perturbations_for
# ---------------------------------------------------------------------------


class TestSearchWrappers:
    def test_search_orderings_from_profile(self, profile) -> None:
        ur = build_ur("rus", "i", profile.inventory, stem_final="s")
        res = search_orderings(profile, ur, "ruʃi")
        assert res.baseline_matched

    def test_search_orderings_movable_perm_requires_indices(
        self, profile,
    ) -> None:
        ur = build_ur("rus", "i", profile.inventory, stem_final="s")
        with pytest.raises(ValueError):
            search_orderings(
                profile, ur, "ruʃi",
                strategy=OrderingStrategy.MOVABLE_PERM,
            )

    def test_search_perturbations_for(self, profile) -> None:
        if not CSV_PATH.exists():
            pytest.skip("lexicon CSV missing")
        rows = list(iter_lexicon_rows(
            CSV_PATH, target_stems=profile.target_stems, limit=1,
        ))
        report = search_perturbations_for(
            profile, rows[0],
            kinds=(PerturbationKind.PATCH_VALUE,),
            limit=10,
        )
        assert report.tried > 0


# ---------------------------------------------------------------------------
# Renderers (pure)
# ---------------------------------------------------------------------------


class TestRenderers:
    def test_render_report_returns_string(
        self, profile, sample_rows,
    ) -> None:
        report = run_validation(profile, sample_rows)
        text = render_report(report)
        assert "BOOLEAN mutation match" in text
        assert "Cluster-specific verification" in text

    def test_render_report_with_distance_summary(
        self, profile, sample_rows,
    ) -> None:
        report = run_validation(profile, sample_rows)
        text = render_report(report, include_distance_summary=True)
        assert "Distance-to-working histogram" in text

    def test_render_trace_returns_string(self, profile) -> None:
        if not CSV_PATH.exists():
            pytest.skip("lexicon CSV missing")
        traces = trace_lemmas(["muscă"], profile, CSV_PATH)
        text = render_trace(traces[0], profile.inventory.base_segments())
        assert "muscă" in text
        assert "muʃte" in text


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------


class TestMismatchCSV:
    def test_write_mismatch_csv(
        self, profile, sample_rows, tmp_path,
    ) -> None:
        report = run_validation(profile, sample_rows)
        out = tmp_path / "mm.csv"
        n = write_mismatch_csv(report, out)
        assert n >= 0
        assert out.exists()
        header = out.read_text().splitlines()[0]
        assert "predicted_sr" in header
        assert "total_distance" in header
