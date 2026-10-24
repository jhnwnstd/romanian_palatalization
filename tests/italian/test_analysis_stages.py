"""The Italian analysis stages end to end on the analysis fixture tree.

Runs the whole pipeline once on ``tests/fixtures/italian`` plus the
analysis excerpt (:mod:`tests.italian.analysis_helpers`) and checks:
the derivative links and the fact-check cases (grecizzare, politichese,
musicista, porcino), the verb classes, the TP long table (R17 rows,
both outcomes, every τ, the non-trigger and coronal N = e rows), the
learners, the reports, the -logo table, the docs, the licence gate of
the tracked outputs (R9) and the stage manifests.
"""

from __future__ import annotations

import csv
from collections.abc import Iterator
from pathlib import Path

import pytest

from lexpipe import csvio, licensing, tp
from lexpipe.langs.italian import (
    docstage,
    logtable,
    schema,
    tplearners,
    tptables,
)
from lexpipe.langs.italian.layout import Layout
from tests.italian import analysis_helpers, helpers

pytestmark = pytest.mark.it


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


@pytest.fixture(scope="module")
def tree(tmp_path_factory: pytest.TempPathFactory) -> Iterator[Path]:
    root = analysis_helpers.build_analysis_tree(
        tmp_path_factory.mktemp("italian_analysis")
    )
    helpers.run_pipeline(analysis_helpers.make_context(root))
    yield root


@pytest.fixture(scope="module")
def layout(tree: Path) -> Layout:
    return Layout(analysis_helpers.make_context(tree))


@pytest.fixture(scope="module")
def derivations(layout: Layout) -> list[dict[str, str]]:
    return _read(layout.derivations)


def _link(
    rows: list[dict[str, str]], derivative: str, base: str
) -> dict[str, str]:
    found = [
        r
        for r in rows
        if r["derivative"] == derivative and r["base_word"] == base
    ]
    assert len(found) == 1, (derivative, base, found)
    return found[0]


# ---------------------------------------------------------------------------
# derivations
# ---------------------------------------------------------------------------


def test_template_links_and_readings(
    derivations: list[dict[str, str]],
) -> None:
    grecista = _link(derivations, "grecista", "greco")
    assert (grecista["edge"], grecista["link_source"]) == (
        "palatal",
        "template",
    )
    assert grecista["deriv_stratum"] == "internal formation"
    fochista = _link(derivations, "fochista", "fuoco")
    assert fochista["edge"] == "velar"  # fuoc → foc (mobile diphthong)
    antichita = _link(derivations, "antichità", "antico")
    assert antichita["edge"] == "velar"
    # the data-layer fix: "From antico + -ità. By surface analysis,
    # borrowed from Latin …" stays an internal formation
    assert antichita["deriv_stratum"] == "internal formation"
    amichevole = _link(derivations, "amichevole", "amico")
    assert (amichevole["edge"], amichevole["deriv_stratum"]) == (
        "velar",
        "inherited",
    )
    amicissimo = _link(derivations, "amicissimo", "amico")
    assert (amicissimo["link_source"], amicissimo["edge"]) == (
        "head_form",
        "palatal",
    )


def test_fact_check_cases(derivations: list[dict[str, str]]) -> None:
    grecizzare = _link(derivations, "grecizzare", "greco")
    assert grecizzare["link_status"] == "rejected:named_other_base"
    assert grecizzare["named_bases"] == "grecia"
    politichese = _link(derivations, "politichese", "politica")
    assert politichese["masc_counterpart"] == "same_root"
    musicista = _link(derivations, "musicista", "musica")
    assert musicista["link_status"] == "accepted"
    assert musicista["masc_counterpart"] == "other_root"
    assert musicista["base_feminine_only"] == "True"
    assert "fact-check" in musicista["fact_check"]
    porcino = _link(derivations, "porcino", "porco")
    assert (porcino["edge"], porcino["deriv_stratum"]) == (
        "palatal",
        "inherited",
    )


def test_candidates_in_both_spellings(
    derivations: list[dict[str, str]],
) -> None:
    leghista = _link(derivations, "leghista", "lega")
    assert leghista["link_tier"] == "candidate"
    assert (leghista["edge"], leghista["link_status"]) == ("velar", "accepted")
    assert leghista["other_spelling"] == "legista"
    assert leghista["attested_ge1"] == "True"
    legista = _link(derivations, "legista", "lega")
    assert (legista["edge"], legista["attested_ge1"]) == ("palatal", "False")
    caricare = _link(derivations, "caricare", "carico")
    assert caricare["link_status"] == "rejected:etymology_unrelated"


# ---------------------------------------------------------------------------
# verbs
# ---------------------------------------------------------------------------


def test_verb_classes(layout: Layout) -> None:
    first: dict[str, dict[str, str]] = {}
    for row in _read(layout.verbs):
        first.setdefault(row["verb"], row)
    expected = {
        "vincere": ("ere_root", "alternating", "True"),
        "cercare": ("are", "velar_throughout", "True"),
        "tacere": ("ere_stressed", "overapplying", "True"),
        "cuocere": ("ere_root", "palatal_throughout", "False"),
        "cucire": ("ire", "palatal_throughout", "False"),
        "fuggire": ("ire", "alternating", "True"),
        "bianchire": ("isc", "velar_throughout", "True"),
        "dire": ("dire", "alternating", "True"),
    }
    for verb, values in expected.items():
        row = first[verb]
        assert (
            row["conj_class"],
            row["paradigm_type"],
            row["in_domain"],
        ) == values, verb
    assert first["trafficare"]["noun_roots"].startswith("it:traffic:")
    assert first["amicarsi"]["verb_key"] == "amicare"


# ---------------------------------------------------------------------------
# tp
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def tp_rows(layout: Layout) -> list[dict[str, str]]:
    return _read(layout.tp_long)


def _row(rows: list[dict[str, str]], **where: str) -> dict[str, str]:
    found = [r for r in rows if all(r[k] == v for k, v in where.items())]
    assert len(found) == 1, (where, len(found))
    return found[0]


def test_tp_long_shape(layout: Layout, tp_rows: list[dict[str, str]]) -> None:
    assert tuple(tp_rows[0]) == tptables.FIELDS
    tables = {r["table"] for r in tp_rows}
    assert tables == {
        "R17a",
        "R17b",
        "R17c",
        "cells",
        "nde",
        "it_stress",
        "it_suffix",
        "it_verbs",
        "it_nounverb",
        "logtable",
        "spellout",
    }
    assert {r["rule"] for r in tp_rows} == {"palatalize", "retain"}
    assert {r["tau"] for r in tp_rows} == {
        "dictionary",
        "ge1",
        "ge0.5pm",
        "ge2pm",
        "dictionary_explicit",
    }
    explicit = {
        r["table"] for r in tp_rows if r["tau"] == "dictionary_explicit"
    }
    assert explicit == {"it_stress"}  # typed plurals only: a sensitivity
    assert {r["evidence_set"] for r in tp_rows} == {"E0", "E1", "E2", "E2_all"}
    for row in tp_rows:  # licence gate of R9: aggregates, by-sa source
        licensing.check_redistributable(licensing.split_sources(row["source"]))


def test_non_trigger_and_coronal_rows(tp_rows: list[dict[str, str]]) -> None:
    common = {
        "evidence_set": "E0",
        "coding": "dominant",
        "V": "",
        "subset": "all",
        "tau": "dictionary",
    }
    e_cell = _row(
        tp_rows, table="R17b", trigger="-e", rule="palatalize", **common
    )
    assert e_cell["N_trigger"] == e_cell["e"]  # amiche, pubbliche, …
    coronal = _row(
        tp_rows,
        table="R17c",
        trigger="-i",
        target_class="coronal",
        rule="palatalize",
        **common,
    )
    assert int(coronal["N_trigger"]) > 0
    assert coronal["N_trigger"] == coronal["e"]  # no coronal rule


def test_dorsal_i_cell_counts(tp_rows: list[dict[str, str]]) -> None:
    row = _row(
        tp_rows,
        table="R17b",
        trigger="-i",
        rule="palatalize",
        subset="all",
        tau="dictionary",
        evidence_set="E0",
        coding="dominant",
        V="",
    )
    retain = _row(
        tp_rows,
        table="R17b",
        trigger="-i",
        rule="retain",
        subset="all",
        tau="dictionary",
        evidence_set="E0",
        coding="dominant",
        V="",
    )
    n = int(row["N_trigger"])
    # palatalize and retain split the same roots; variable roots count
    # against both
    assert int(row["e"]) + int(retain["e"]) >= n
    assert row["fragility"] != ""
    expected = tp.tolerated(n, int(row["e"]))
    assert row["tolerated_trigger"] == str(expected)


def test_learners(layout: Layout) -> None:
    rows = _read(layout.learners)
    assert tuple(rows[0]) == tplearners.LEARNER_FIELDS
    assert {int(r["V"]) for r in rows} == set(analysis_helpers.LEARNER_V)
    for row in rows:
        shares = sum(
            float(row[f"share_{o}"])
            for o in ("palatalize", "retain", "neither", "conditioned")
        )
        assert shares == pytest.approx(1.0)
        assert int(row["n_samples"]) == analysis_helpers.LEARNER_SAMPLES
        licensing.check_redistributable(licensing.split_sources(row["source"]))


# ---------------------------------------------------------------------------
# reports, logtable, docs
# ---------------------------------------------------------------------------


def test_reports(layout: Layout) -> None:
    summary = layout.tp_summary.read_text(encoding="utf-8")
    assert "Cross-language rule family (R17)" in summary
    assert "nested set giavazzi" in summary
    covariation = layout.covariation.read_text(encoding="utf-8")
    assert "Mantel-Haenszel" in covariation
    # Skipped without statsmodels; with it, fitted or (on a fixture this
    # small) not estimable, which lexpipe.stats reports as ValueError.
    assert any(
        line in covariation
        for line in (
            "logistic model skipped",
            "logistic model (",
            "logistic model failed",
        )
    )
    assert "grecizzare" in covariation  # fact-check section
    assert "Pair set verbalizer" in layout.spellout.read_text(encoding="utf-8")
    noev = layout.noevidence.read_text(encoding="utf-8")
    assert "no-evidence units by base descriptor" in noev
    assert "musica → musicista" in noev


def test_logtable(layout: Layout) -> None:
    rows = {r["lemma"]: r for r in _read(layout.logtable)}
    assert tuple(next(iter(rows.values()))) == logtable.LOGTABLE_FIELDS
    assert (rows["teologo"]["semantic_class"], rows["teologo"]["plural"]) == (
        "person",
        "teologi",
    )
    assert rows["teologo"]["abstract_noun"] == "teologia"
    assert rows["catalogo"]["plural_outcome"] == "velar"
    assert rows["chirurgo"]["plural_outcome"] == "variable"
    for row in rows.values():  # no Leipzig column in the tracked table
        licensing.check_redistributable(licensing.split_sources(row["source"]))
    detail = _read(layout.logtable_attestation)
    assert {"count_cs", "count_folded"} <= set(detail[0])


def test_docs(layout: Layout) -> None:
    loaded = schema.load_schema(layout.schema)
    assert schema.header_problems(loaded, layout.ctx.root) == []
    data_dict = layout.data_dict.read_text(encoding="utf-8")
    for item in loaded.values():
        assert f"`{item.path}`" in data_dict
    readme = layout.readme.read_text(encoding="utf-8")
    assert docstage.BEGIN in readme and docstage.END in readme
    assert "Headline numbers" in readme


def test_tracked_csvs_pass_the_licence_gate(layout: Layout) -> None:
    """The checks of :func:`lexpipe.csvio.write_tracked`, re-applied to
    every tracked CSV the schema lists (R9)."""
    loaded = schema.load_schema(layout.schema)
    for item in schema.tracked_files(loaded):
        path = layout.ctx.root / item.path
        assert licensing.is_tracked_location(path, layout.ctx.root), path
        for row in csvio.read_rows(path):
            sources = licensing.split_sources(row["source"])
            licensing.check_redistributable(sources)
            if row.get("licence_class"):
                licensing.check_licence_class(row["licence_class"])


def test_stages_are_skipped_when_up_to_date(
    tree: Path, layout: Layout, capsys: pytest.CaptureFixture[str]
) -> None:
    before = {p: p.read_bytes() for p in layout.ctx.manifest_dir.iterdir()}
    assert {p.name for p in before} >= {
        f"{stage}.manifest.tsv" for stage in analysis_helpers.ANALYSIS_STAGES
    }
    helpers.run_pipeline(analysis_helpers.make_context(tree))
    out = capsys.readouterr().out
    assert out.count("up to date") == len(before)
    after = {p: p.read_bytes() for p in layout.ctx.manifest_dir.iterdir()}
    assert after == before
