"""Parity of the shared core with the legacy Romanian R tables (§7, §8).

Fixture tests check the reproduction and the root-by-root
reconciliation mechanics; ``needs_data`` tests run them on the real
lexicon against the tracked ``analysis/romanian_tp_summary_*.csv`` and
pin the headline numbers of ``analysis/comparison/romanian_parity.txt``.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lexpipe import tp
from lexpipe.langs.romanian_compare import harmonized as hz
from lexpipe.langs.romanian_compare import lexicon as lx
from lexpipe.langs.romanian_compare import parity
from lexpipe.langs.romanian_compare.corpus import load_counts
from lexpipe.langs.romanian_compare.observations import Attester
from lexpipe.langs.romanian_compare.records import (
    RootsFile,
    build_records,
    read_verdicts,
)
from tests.romanian_compare.helpers import (
    FIXTURES,
    REAL_FREQ,
    REAL_LEXICON,
    REAL_SUMMARY_DS,
    REAL_SUMMARY_FULL,
    REAL_VERDICTS,
    real_data_or_skip,
    write_legacy_csv,
)

pytestmark = pytest.mark.ro


def load(
    lexicon: Path, verdicts: Path, freq: Path
) -> tuple[list[lx.LexiconRow], RootsFile]:
    rows = lx.read_lexicon(lexicon)
    counts = load_counts(freq)
    records = build_records(rows, read_verdicts(verdicts), counts)
    return rows, RootsFile(counts.total_tokens, {}, tuple(records))


@pytest.fixture(scope="module")
def fixture_data() -> tuple[list[lx.LexiconRow], RootsFile]:
    return load(
        FIXTURES / "lexicon.csv",
        FIXTURES / "verb_paradigm_verdict.csv",
        FIXTURES / "freq.csv",
    )


@pytest.fixture(scope="module")
def real_data() -> tuple[list[lx.LexiconRow], RootsFile]:
    real_data_or_skip()
    return load(REAL_LEXICON, REAL_VERDICTS, REAL_FREQ)


def test_detector_matches_the_r_function() -> None:
    detect = parity.detect_palatal_from_orth
    assert detect("colac", "încolăci|colăci|colaci", "c") is True
    assert detect("bloc", "bloca|debloca", "c") is False
    assert detect("verde", "înverzi", "d") is True
    assert detect("obraz", "obrazi", "z") is False
    assert detect("x", "", "c") is None
    assert detect("x", "y", "r") is None


def test_opportunity_tp_patch() -> None:
    base = lx.parse_row(
        dict.fromkeys(lx.REQUIRED_COLUMNS, "")
        | {
            "lemma": "ureche",
            "pos": "N",
            "plural": "urechi",
            "opportunity": "none",
            "exception_reason": "nde:paduchi",
            "mutation": "False",
        }
    )
    assert parity.opportunity_tp(base) == "i"


def test_fixture_reproduction_round_trip(
    fixture_data: tuple[list[lx.LexiconRow], RootsFile], tmp_path: Path
) -> None:
    rows, _ = fixture_data
    full, small = parity.legacy_tables(rows)
    path = tmp_path / "full.csv"
    write_legacy_csv(path, full)
    lines = parity.compare_csv(path, full)
    assert lines and all(line.ok for line in lines)
    assert full["<z> + <-i> plural"] == (
        1,
        0,
    )  # albanez is excluded by the R rules
    assert full["<t> + <-e> plural"][0] == 0
    assert all(
        label.endswith("(downsampled)")
        or "WITH NDEB" in label
        or label.endswith("type")
        for label in small
    )


def test_fixture_reconciliation_balances(
    fixture_data: tuple[list[lx.LexiconRow], RootsFile],
) -> None:
    _, roots = fixture_data
    cache = hz.ObservationCache(roots.records, Attester(roots))
    for seg in parity.SEGMENTS:
        for opp in ("i", "e"):
            assert parity.reconcile(roots.records, cache, seg, opp).balanced, (
                seg,
                opp,
            )
    z = parity.reconcile(roots.records, cache, "z", "i")
    assert z.added == {parity.ADD_REASONS[7]: (1, 1)}  # albanez


@pytest.mark.needs_data
def test_legacy_rows_are_reproduced_exactly(
    real_data: tuple[list[lx.LexiconRow], RootsFile],
) -> None:
    rows, _ = real_data
    full, small = parity.legacy_tables(rows)
    for path, computed in (
        (REAL_SUMMARY_FULL, full),
        (REAL_SUMMARY_DS, small),
    ):
        lines = parity.compare_csv(path, computed)
        assert len(lines) == 37
        assert [line.label for line in lines if not line.ok] == []
    assert full["<c> + <-i> plural"] == (473, 0)
    assert full["<t> + <-e> plural"] == (1, 2048)
    assert small["<c> + <-i> plural (downsampled)"] == (28, 0)


@pytest.mark.needs_data
def test_lemma_frequencies_match_the_lexicon_column(
    real_data: tuple[list[lx.LexiconRow], RootsFile],
) -> None:
    rows, _ = real_data
    counts = load_counts(REAL_FREQ)
    assert counts.total_tokens == 18_628_329
    mismatches = [r.lemma for r in rows if counts.count(r.lemma) != r.freq]
    assert mismatches == []


@pytest.mark.needs_data
def test_real_reconciliation(
    real_data: tuple[list[lx.LexiconRow], RootsFile],
) -> None:
    _, roots = real_data
    cache = hz.ObservationCache(roots.records, Attester(roots))
    results = {
        (seg, opp): parity.reconcile(roots.records, cache, seg, opp)
        for seg in parity.SEGMENTS
        for opp in ("i", "e")
    }
    assert all(r.balanced for r in results.values())
    c_i = results[("c", "i")]
    assert (c_i.legacy_n, c_i.legacy_e, c_i.harm_n, c_i.harm_e) == (
        473,
        0,
        1236,
        49,
    )
    z_i = results[("z", "i")]
    assert (z_i.legacy_n, z_i.legacy_e, z_i.harm_n, z_i.harm_e) == (
        14,
        4,
        126,
        116,
    )
    assert z_i.added[parity.ADD_REASONS[7]] == (
        112,
        112,
    )  # z→zi "foreign" plurals
    t_e = results[("t", "e")]
    assert t_e.added[parity.ADD_REASONS[6]] == (
        6,
        6,
    )  # fată, pată, ...: suppletion rule
    # The core subset, which drops outcome-blind flags, is close to the
    # legacy domain: compare <c> + <-i> and <g> + <-i, -e>.
    observed = cache.get("E0")
    core_c = [
        observed[r.root_id]
        for r in hz.select(roots.records, hz.TARGETS["c"], hz.SUBSETS["core"])
    ]
    row = tp.harmonized_row(core_c, tp.CellSpec("ro", "cells", "c", "-i"))
    assert (row["N_trigger"], row["e"]) == (472, 0)
