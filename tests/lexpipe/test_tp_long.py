"""Tests for the harmonized long TP rows of lexpipe.tp (R16-R19).

One synthetic k + -i cell exercises every root status: alternating,
non-alternating, dominant and true doublets, a shielded root (only an
N allomorph attested), a root without evidence, one whose alternation
is unknown, and one attested only in the dictionary.
"""

from __future__ import annotations

import importlib
import math
import sys
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lexpipe import csvio, tp  # noqa: E402
from lexpipe.tp import (  # noqa: E402
    CellSpec,
    NestedRule,
    Observation,
    RootObservation,
    RootStatus,
)

CORPUS = {"dictionary": True, "ge1": True, "ge0.5pm": True, "ge2pm": True}
DICTIONARY = {"dictionary": True, "ge1": False, "ge0.5pm": False}


def obs(
    alternates: bool | None,
    allomorph: str = "P",
    attested: Mapping[str, bool] = CORPUS,
    form: str = "",
    count: int = 0,
) -> Observation:
    return Observation(alternates, attested, allomorph, form, count)


def root(
    root_id: str, *observations: Observation, trigger: str = "-i"
) -> RootObservation:
    return RootObservation(root_id, "k", 1.0, {trigger: observations})


ROOTS = [
    root("amico", obs(True)),
    root("medico", obs(True)),
    root("parco", obs(False)),
    root("fuoco", obs(False)),
    root(  # 74 : 0 tokens -> dominant non-alternating variant
        "chirurgo",
        obs(False, form="chirurghi", count=74),
        obs(True, form="chirurgi", count=0),
    ),
    root(  # 15 : 3 = 83% -> a true doublet
        "manico",
        obs(True, form="manici", count=15),
        obs(False, form="manichi", count=3),
    ),
    root(  # no tokens at all -> a true doublet, tied
        "pari",
        obs(True, attested=DICTIONARY, form="a"),
        obs(False, attested=DICTIONARY, form="b"),
    ),
    root("bosco", obs(None, allomorph="N")),
    root("lago", obs(True), trigger="-e"),
    root("ignoto", obs(None)),
    root("raro", obs(True, attested=DICTIONARY)),
]


def cell(tau: str = "dictionary", coding: str = "dominant") -> CellSpec:
    return CellSpec(
        lang="it",
        table="nouns",
        target_class="k",
        trigger="-i",
        tau=tau,
        coding=coding,
        V=None,
    )


def statuses(tau: str, coding: str) -> dict[str, RootStatus]:
    return {r.root_id: tp.root_status(r, "-i", tau, coding) for r in ROOTS}


def test_root_status_dictionary_dominant() -> None:
    S = RootStatus
    assert statuses("dictionary", "dominant") == {
        "amico": S.ALTERNATES,
        "medico": S.ALTERNATES,
        "parco": S.NONALT,
        "fuoco": S.NONALT,
        "chirurgo": S.NONALT,
        "manico": S.ALTERNATES,  # the more frequent outcome
        "pari": S.VARIABLE,  # tie: counted as an exception
        "bosco": S.SHIELDED,
        "lago": S.NO_EVIDENCE,
        "ignoto": S.UNKNOWN,
        "raro": S.ALTERNATES,
    }


def test_codings_differ_only_on_true_doublets() -> None:
    S = RootStatus
    exception = statuses("dictionary", "doublet_exception")
    excluded = statuses("dictionary", "doublet_excluded")
    assert exception["manico"] is S.VARIABLE
    assert exception["pari"] is S.VARIABLE
    assert excluded["manico"] is excluded["pari"] is S.EXCLUDED
    for name in ("chirurgo", "amico", "parco", "bosco"):
        dominant = statuses("dictionary", "dominant")[name]
        assert exception[name] is excluded[name] is dominant


def test_attestation_threshold() -> None:
    S = RootStatus
    at_one_token = statuses("ge1", "dominant")
    assert at_one_token["raro"] is S.NO_EVIDENCE
    assert at_one_token["pari"] is S.NO_EVIDENCE
    assert at_one_token["amico"] is S.ALTERNATES
    # A threshold missing from attested_by_tau is not attested.
    assert statuses("ge2pm", "dominant")["raro"] is S.NO_EVIDENCE


def test_two_outcome_rows() -> None:
    pal, ret = tp.two_outcome_rows(ROOTS, cell())
    assert list(pal) == list(tp.TP_LONG_FIELDS)
    assert (pal["rule"], ret["rule"]) == ("palatalize", "retain")
    for row in (pal, ret):
        assert row["N_target"] == 11
        assert row["N_trigger"] == 8  # 4 alternate, 3 do not, 1 variable
        assert row["e_nonalt"] == 3
        assert row["n_shielded"] == 1
        assert row["theta_trigger"] == pytest.approx(8 / math.log(8))
        assert row["theta_target"] == pytest.approx(11 / math.log(11))
        assert row["subset"] == "all" and row["V"] is None
    assert pal["e"] == 4  # nonalt + variable
    assert pal["tolerated_trigger"] is False  # 4 > 3.85
    assert pal["tolerated_counterfactual"] is False  # 11 - 4 = 7 > 4.59
    assert pal["majority"] is None  # 4 vs 4
    assert ret["e"] == 5  # alternates + variable
    assert ret["tolerated_trigger"] is False
    assert ret["tolerated_counterfactual"] is False  # 5 > 4.59
    assert ret["majority"] is False


def test_excluded_doublets_leave_the_domain() -> None:
    pal, ret = tp.two_outcome_rows(ROOTS, cell(coding="doublet_excluded"))
    assert (pal["N_target"], pal["N_trigger"], pal["e"]) == (9, 6, 3)
    assert pal["tolerated_trigger"] is True  # 3 <= 6 / ln 6
    assert pal["tolerated_counterfactual"] is False  # 9 - 3 = 6 > 4.10
    assert ret["e"] == 3 and ret["tolerated_trigger"] is True


def test_doublet_exception_counts_against_both_rules() -> None:
    pal, ret = tp.two_outcome_rows(ROOTS, cell(coding="doublet_exception"))
    assert pal["N_trigger"] == ret["N_trigger"] == 8
    assert pal["e_nonalt"] == 3
    assert pal["e"] == 3 + 2 and ret["e"] == 3 + 2


def test_corpus_threshold_row() -> None:
    row = tp.harmonized_row(ROOTS, cell(tau="ge1"))
    assert (row["N_target"], row["N_trigger"], row["e"]) == (11, 6, 3)
    assert row["tolerated_trigger"] is True
    assert row["tau"] == "ge1" and row["coding"] == "dominant"


def test_small_cells_have_no_theta() -> None:
    row = tp.harmonized_row(ROOTS[:1], cell())
    assert row["N_trigger"] == 1
    assert row["theta_trigger"] is None and row["theta_target"] is None
    assert row["tolerated_trigger"] is None
    empty = tp.harmonized_row([], cell())
    assert (empty["N_target"], empty["e"], empty["majority"]) == (0, 0, None)


def test_nested_rows_in_yangs_order() -> None:
    proparoxytone = {"medico", "manico", "parco"}
    specific = [
        NestedRule(
            "proparoxytone", "palatalize", lambda r: r.root_id in proparoxytone
        ),
        # Claims what the first rule left: medico is not counted twice.
        NestedRule(
            "velar-final",
            "retain",
            lambda r: r.root_id in {"medico", "fuoco", "chirurgo"},
        ),
    ]
    rows = tp.nested_rows(ROOTS, cell(), specific)
    assert [(r["subset"], r["rule"]) for r in rows] == [
        ("proparoxytone", "palatalize"),
        ("velar-final", "retain"),
        ("elsewhere", "palatalize"),
    ]
    assert [r["N_target"] for r in rows] == [3, 2, 6]
    assert rows[0]["e"] == 1  # parco
    assert rows[1]["e"] == 0  # fuoco and chirurgo both retain
    assert rows[1]["tolerated_trigger"] is True


def test_rows_are_written_with_csvio(tmp_path: Path) -> None:
    rows = tp.two_outcome_rows(ROOTS, cell(), memo="synthetic")
    path = tmp_path / "tp_long.csv"
    assert csvio.write_rows(path, rows, tp.TP_LONG_FIELDS) == 2
    written = csvio.read_rows(path)
    assert written[0]["V"] == "" and written[0]["majority"] == ""
    assert written[1]["tolerated_trigger"] == "False"
    assert written[0]["memo"] == "synthetic"


def test_counts_do_not_depend_on_order() -> None:
    forward = tp.two_outcome_rows(ROOTS, cell())
    assert tp.two_outcome_rows(ROOTS[::-1], cell()) == forward


def test_duplicate_roots_are_refused() -> None:
    with pytest.raises(ValueError, match="listed twice"):
        tp.harmonized_row([ROOTS[0], ROOTS[0]], cell())


@pytest.mark.parametrize(
    ("counts", "expected"),
    [
        ({"chirurghi": 74, "chirurgi": 0}, "dominant:chirurghi"),
        ({"a": 9, "b": 1}, "dominant:a"),  # exactly 90%
        ({"a": 8, "b": 1}, "true_doublet"),  # 88.9%
        ({"a": 15, "b": 3}, "true_doublet"),
        ({"a": 0, "b": 0}, "true_doublet"),
        ({"only": 3}, "dominant:only"),
    ],
)
def test_classify_doublet(counts: dict[str, int], expected: str) -> None:
    assert tp.classify_doublet(counts) == expected


@pytest.mark.parametrize("counts", [{}, {"a": -1, "b": 2}])
def test_classify_doublet_rejects(counts: dict[str, int]) -> None:
    with pytest.raises(ValueError):
        tp.classify_doublet(counts)


@pytest.mark.parametrize(
    "build",
    [
        lambda: Observation(True, CORPUS, "X"),
        lambda: Observation(True, CORPUS, "P", count=-1),
        lambda: RootObservation("x", "k", -1.0, {}),
        lambda: RootObservation("x", "k", math.nan, {}),
        lambda: CellSpec("it", "t", "k", "-i", coding="excluded"),
        lambda: NestedRule("s", "mutate", lambda r: True),
    ],
)
def test_validation(build: Any) -> None:
    with pytest.raises(ValueError):
        build()


def test_records_are_frozen() -> None:
    attested = dict(CORPUS)
    observation = Observation(True, attested, "P")
    attested["ge1"] = False
    assert observation.attested_by_tau["ge1"] is True
    record = root("x", observation)
    with pytest.raises(TypeError):
        record.triggers["-e"] = ()  # type: ignore[index]


# ---------------------------------------------------------------------------
# Tracked columns
# ---------------------------------------------------------------------------


def test_tracked_fields() -> None:
    assert tp.TP_LONG_EXTRA_FIELDS == ("fragility", "source", "licence_class")
    assert tp.TP_LONG_TRACKED_FIELDS == (
        tp.TP_LONG_FIELDS + tp.TP_LONG_EXTRA_FIELDS
    )
    assert len(set(tp.TP_LONG_TRACKED_FIELDS)) == len(
        tp.TP_LONG_TRACKED_FIELDS
    )


@pytest.mark.parametrize(
    "module",
    [
        "lexpipe.langs.italian.tptables",
        "lexpipe.langs.czech.tpstage",
        "lexpipe.langs.romanian_compare.harmonized",
    ],
)
def test_plugin_tp_long_headers_are_the_tracked_fields(module: str) -> None:
    """The three plugins' long tables already use this header."""
    fields = importlib.import_module(module).FIELDS
    assert tuple(fields) == tp.TP_LONG_TRACKED_FIELDS


def test_tracked_row(tmp_path: Path) -> None:
    pal, _ = tp.two_outcome_rows(ROOTS, cell())
    row = tp.tracked_row(pal, source="wiktionary_en", licence_class="by-sa")
    assert list(row) == list(tp.TP_LONG_TRACKED_FIELDS)
    assert row["fragility"] == tp.fragility(pal["N_trigger"], pal["e"]) == 1
    assert (row["source"], row["licence_class"]) == ("wiktionary_en", "by-sa")
    assert "fragility" not in pal  # the input row is not modified
    both = tp.tracked_row(
        pal,
        source=("wiktionary_en", "wiktionary_cs"),
        licence_class="by-sa",
    )
    assert both["source"] == ["wiktionary_en", "wiktionary_cs"]
    path = tmp_path / "analysis" / "italian" / "tp_long.csv"
    written = csvio.write_tracked(
        path, [row, both], tp.TP_LONG_TRACKED_FIELDS, root=tmp_path
    )
    assert written == 2
    assert csvio.read_rows(path)[1]["source"] == "wiktionary_en|wiktionary_cs"


# ---------------------------------------------------------------------------
# Precomputed statuses
# ---------------------------------------------------------------------------

TAUS_AND_CODINGS = [
    (tau, coding) for tau in ("dictionary", "ge1") for coding in tp.CODINGS
]


@pytest.mark.parametrize(("tau", "coding"), TAUS_AND_CODINGS)
def test_status_rows_equal_root_rows(tau: str, coding: str) -> None:
    spec_cell = cell(tau=tau, coding=coding)
    statuses = tp.cell_statuses(ROOTS, spec_cell)
    expected = tp.two_outcome_rows(ROOTS, spec_cell, memo="m")
    assert (
        tp.status_rows(
            statuses.of(r.root_id for r in ROOTS), spec_cell, memo="m"
        )
        == expected
    )
    assert (
        tp.two_outcome_rows(ROOTS, spec_cell, memo="m", statuses=statuses)
        == expected
    )
    for rule in tp.RULES:
        single = tp.harmonized_row(ROOTS, spec_cell, rule=rule, subset="s")
        assert (
            tp.status_row(
                list(statuses.statuses.values()),
                spec_cell,
                rule=rule,
                subset="s",
            )
            == single
        )
        assert (
            tp.harmonized_row(
                ROOTS, spec_cell, rule=rule, subset="s", statuses=statuses
            )
            == single
        )
        assert (
            tp.counts_row(
                tp.domain_counts(ROOTS, spec_cell),
                spec_cell,
                rule=rule,
                subset="s",
            )
            == single
        )


def test_statuses_cover_subsets_and_nested_rules() -> None:
    spec_cell = cell()
    statuses = tp.cell_statuses(ROOTS, spec_cell)
    subset = ROOTS[2:7]
    assert tp.domain_counts(subset, spec_cell, statuses=statuses) == (
        tp.domain_counts(subset, spec_cell)
    )
    assert statuses.counts(r.root_id for r in subset) == tp.domain_counts(
        subset, spec_cell
    )
    specific = [
        NestedRule("pa", "palatalize", lambda r: r.root_id.startswith("p")),
        NestedRule("m", "retain", lambda r: r.root_id.startswith("m")),
    ]
    assert tp.nested_rows(ROOTS, spec_cell, specific, statuses=statuses) == (
        tp.nested_rows(ROOTS, spec_cell, specific)
    )


def test_status_counts_accept_strings_and_check_them() -> None:
    S = RootStatus
    counts = tp.status_counts(
        ["alternates", S.NONALT, "variable", "excluded", "shielded", "unknown"]
    )
    assert counts == tp.DomainCounts(
        n_target=5, n_alternating=1, n_nonalt=1, n_variable=1, n_shielded=1
    )
    assert tp.status_counts([]) == tp.DomainCounts(0, 0, 0, 0, 0)
    with pytest.raises(ValueError):
        tp.status_counts(["mutates"])


def test_cell_statuses_refuse_other_cells_and_roots() -> None:
    spec_cell = cell()
    statuses = tp.cell_statuses(ROOTS, spec_cell)
    assert statuses.matches(spec_cell)
    # Labels that do not change a status may differ ...
    relabelled = CellSpec("it", "other_table", "dorsal", "-i", V=500)
    assert statuses.matches(relabelled)
    tp.harmonized_row(ROOTS, relabelled, statuses=statuses)
    # ... trigger, tau, coding and evidence set may not.
    for other in (
        cell(tau="ge1"),
        cell(coding="doublet_excluded"),
        CellSpec("it", "nouns", "k", "-e"),
        CellSpec("it", "nouns", "k", "-i", evidence_set="E1"),
    ):
        with pytest.raises(ValueError, match="cannot count"):
            tp.harmonized_row(ROOTS, other, statuses=statuses)
    stranger = root("straniero", obs(True))
    with pytest.raises(ValueError, match="no precomputed status"):
        tp.harmonized_row([stranger], spec_cell, statuses=statuses)
    with pytest.raises(ValueError, match="listed twice"):
        tp.cell_statuses([ROOTS[0], ROOTS[0]], spec_cell)
    with pytest.raises(ValueError, match="listed twice"):
        statuses.counts(["amico", "amico"])
    with pytest.raises(TypeError):
        statuses.statuses["amico"] = RootStatus.NONALT  # type: ignore[index]


def test_statuses_by_coding_in_one_pass() -> None:
    by_coding = tp.cell_statuses_by_coding(ROOTS, cell(tau="dictionary"))
    assert list(by_coding) == list(tp.CODINGS)
    for coding, statuses in by_coding.items():
        assert statuses.coding == coding
        assert dict(statuses.statuses) == statuses_of(coding)
    for item in ROOTS:
        assert tp.root_status_by_coding(item, "-i", "ge1") == {
            coding: tp.root_status(item, "-i", "ge1", coding)
            for coding in tp.CODINGS
        }
    with pytest.raises(ValueError, match="coding"):
        tp.root_status_by_coding(ROOTS[0], "-i", "ge1", ("best",))
    with pytest.raises(ValueError, match="coding"):
        tp.CellStatuses("-i", "ge1", "best", "E0", {})


def statuses_of(coding: str) -> dict[str, RootStatus]:
    return statuses("dictionary", coding)


# ---------------------------------------------------------------------------
# Optional refinements: cells (Observation.cell) and a doublet minimum
# (CellSpec.doublet_min_tokens); both off by default
# ---------------------------------------------------------------------------


def cell_obs(
    alternates: bool, cell_name: str, form: str, count: int
) -> Observation:
    return Observation(alternates, CORPUS, "P", form, count, cell_name)


def test_cells_that_disagree_are_variable_under_every_coding() -> None:
    # amici (-i, 3012 tokens) and amiche (-e, 400): one root, two cells.
    amico = root(
        "amico",
        cell_obs(True, "i", "amici", 3012),
        cell_obs(False, "e", "amiche", 400),
        trigger="-i/-e",
    )
    for coding in tp.CODINGS:
        assert tp.root_status(amico, "-i/-e", "ge1", coding) is (
            RootStatus.VARIABLE
        )
    # Without cells the two forms are one doublet, coded by tokens
    # (88% is not dominant; the dominant coding takes the majority).
    pooled = root(
        "amico",
        obs(True, form="amici", count=3012),
        obs(False, form="amiche", count=400),
        trigger="-i/-e",
    )
    assert tp.root_status(pooled, "-i/-e", "ge1", "dominant") is (
        RootStatus.ALTERNATES
    )
    assert tp.root_status(pooled, "-i/-e", "ge1", "doublet_excluded") is (
        RootStatus.EXCLUDED
    )


def test_a_doublet_is_coded_within_its_cell() -> None:
    S = RootStatus
    # -i cell: dominant palatal doublet; -e cell: palatal too.
    agree = root(
        "x",
        cell_obs(True, "i", "a", 95),
        cell_obs(False, "i", "b", 5),
        cell_obs(True, "e", "c", 10),
        trigger="-i/-e",
    )
    for coding in tp.CODINGS:
        assert tp.root_status(agree, "-i/-e", "ge1", coding) is S.ALTERNATES
    # A true doublet in the -i cell next to a velar -e cell.
    mixed = root(
        "y",
        cell_obs(True, "i", "a", 6),
        cell_obs(False, "i", "b", 4),
        cell_obs(False, "e", "c", 10),
        trigger="-i/-e",
    )
    assert tp.root_status_by_coding(mixed, "-i/-e", "ge1") == {
        "dominant": S.VARIABLE,  # -i resolves to alternates, -e does not
        "doublet_exception": S.VARIABLE,
        "doublet_excluded": S.EXCLUDED,
    }
    same = root(
        "z",
        cell_obs(True, "i", "a", 6),
        cell_obs(False, "i", "b", 4),
        cell_obs(True, "e", "c", 10),
        trigger="-i/-e",
    )
    assert tp.root_status(same, "-i/-e", "ge1", "dominant") is S.ALTERNATES


def test_default_rows_ignore_the_refinements() -> None:
    # The fixture roots have no cells: rows equal the unrefined ones.
    spec_cell = cell()
    assert spec_cell.doublet_min_tokens == 0
    statuses = tp.cell_statuses(ROOTS, spec_cell)
    assert statuses.doublet_min_tokens == 0
    assert dict(statuses.statuses) == statuses_of("dominant")


def test_doublet_minimum_tokens() -> None:
    S = RootStatus
    assert tp.classify_doublet({"a": 3, "b": 0}) == "dominant:a"
    assert tp.classify_doublet({"a": 3, "b": 0}, min_total=10) == (
        "true_doublet"
    )
    assert tp.classify_doublet({"a": 74, "b": 0}, min_total=10) == (
        "dominant:a"
    )
    astrologo = root(
        "astrologo",
        obs(True, form="astrologi", count=3),
        obs(False, form="astrologhi", count=0),
    )
    for coding in tp.CODINGS:
        assert tp.root_status(astrologo, "-i", "ge1", coding) is (S.ALTERNATES)
    assert tp.root_status_by_coding(astrologo, "-i", "ge1", min_tokens=10) == {
        "dominant": S.ALTERNATES,  # a true doublet: the majority
        "doublet_exception": S.VARIABLE,
        "doublet_excluded": S.EXCLUDED,
    }
    strict = CellSpec(
        "it", "nouns", "k", "-i", tau="ge1", doublet_min_tokens=10
    )
    by_coding = tp.cell_statuses_by_coding([astrologo], strict)
    assert by_coding["doublet_excluded"].status("astrologo") is S.EXCLUDED
    assert by_coding["doublet_excluded"].doublet_min_tokens == 10
    row = tp.harmonized_row(
        [astrologo],
        CellSpec(
            "it",
            "nouns",
            "k",
            "-i",
            tau="ge1",
            coding="doublet_exception",
            doublet_min_tokens=10,
        ),
    )
    assert (row["N_trigger"], row["e"]) == (1, 1)
    # Statuses computed with another minimum are refused.
    plain = tp.cell_statuses([astrologo], CellSpec("it", "n", "k", "-i"))
    with pytest.raises(ValueError, match="cannot count"):
        tp.harmonized_row([astrologo], strict, statuses=plain)
    with pytest.raises(ValueError, match="doublet_min_tokens"):
        CellSpec("it", "n", "k", "-i", doublet_min_tokens=-1)
