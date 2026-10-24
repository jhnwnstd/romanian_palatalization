"""End-to-end run of every stage on the fixture tree (``--root``)."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from lexpipe import csvio, tp
from lexpipe.context import RunContext, read_manifest
from lexpipe.langs.romanian_compare import PLUGIN, RO_INFO
from lexpipe.langs.romanian_compare import __main__ as cli
from lexpipe.langs.romanian_compare import harmonized as hz
from lexpipe.langs.romanian_compare import readme
from lexpipe.langs.romanian_compare import stages as st
from lexpipe.langs.romanian_compare.records import read_roots

pytestmark = pytest.mark.ro

PARAMS = {"learner_samples": 25, "learner_V": [5, 10]}


def run_all(root: Path, **params: object) -> RunContext:
    ctx = RunContext(RO_INFO, root=root, params={**PARAMS, **params})
    for name, stage in PLUGIN.stages.items():
        assert stage(ctx) in (0, None), name
    return ctx


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_all_stages_write_only_their_outputs(tree: Path) -> None:
    before = {p: digest(p) for p in tree.rglob("*") if p.is_file()}
    ctx = run_all(tree)
    p = st.paths(ctx)
    for path in (
        p.roots,
        p.tp_long,
        p.learners,
        p.spellout,
        p.noevidence,
        p.noevidence_units,
        p.covariation,
        p.covariation_units,
        p.logtable,
        p.logtable_full,
        p.parity,
        p.readme,
    ):
        assert path.is_file(), path
    # Inputs untouched; nothing written outside the two directories.
    for path, sha in before.items():
        assert digest(path) == sha, path
    allowed = (
        tree / "analysis" / "comparison",
        tree / "data" / "romanian_compare",
    )
    for path in tree.rglob("*"):
        if path.is_file() and path not in before:
            assert any(base in path.parents for base in allowed), path
    for name in PLUGIN.stages:
        manifest = read_manifest(ctx.manifest_path(name))
        assert manifest is not None and manifest.stage == name


def test_tracked_outputs_pass_the_gate_and_hold_no_counts(tree: Path) -> None:
    ctx = run_all(tree)
    p = st.paths(ctx)
    rows = csvio.read_rows(p.tp_long)
    assert list(rows[0]) == list(hz.FIELDS) == list(tp.TP_LONG_TRACKED_FIELDS)
    # The harvested lexicon draws on both Wiktionary editions (R9).
    assert {r["source"] for r in rows} == {"wiktionary_ro|wiktionary_en"}
    assert {r["licence_class"] for r in rows} == {"by-sa"}
    tables = {r["table"] for r in rows}
    assert tables == {
        "cells",
        "R17a",
        "R17b",
        "R17c",
        "nde",
        "logtable",
        "spellout",
    }
    assert {r["rule"] for r in rows} == set(tp.RULES)
    # Per-word Leipzig counts only in the untracked full -log table;
    # each row names its own edition.
    log = csvio.read_rows(p.logtable)
    assert log and not any("count" in column for column in log[0])
    assert {r["source"] for r in log} <= {"wiktionary_ro", "wiktionary_en"}
    assert all(
        r["source"] == "wiktionary_" + r["lexicon_source"][:2] for r in log
    )
    # Text reports carry no licence header (redistributable sources only).
    for path in (p.spellout, p.noevidence, p.covariation, p.parity):
        assert "Licence and attribution" not in path.read_text("utf-8")
    full = csvio.read_rows(p.logtable_full)
    assert {"lemma_count", "plural_count"} <= set(full[0])
    learners_rows = csvio.read_rows(p.learners)
    assert {r["V"] for r in learners_rows} == {"5", "10"}
    for r in learners_rows:
        total = sum(
            float(r[f"share_{k}"])
            for k in ("palatalize", "retain", "neither", "conditioned")
        )
        assert total == pytest.approx(1.0)


def test_key_fixture_cells(tree: Path) -> None:
    ctx = run_all(tree)
    rows = csvio.read_rows(st.paths(ctx).tp_long)

    def row(**match: str) -> dict[str, str]:
        found = [r for r in rows if all(r[k] == v for k, v in match.items())]
        assert len(found) == 1, match
        return found[0]

    base = {
        "tau": "dictionary",
        "evidence_set": "E0",
        "coding": "dominant",
        "V": "",
        "rule": "palatalize",
    }
    z = row(table="R17c", target_class="z", trigger="-i", subset="all", **base)
    # obraz/obraji alternates; albanez/albanezi (the legacy "foreign
    # plural" exclusion) is kept, as an exception.
    assert (z["N_trigger"], z["e"]) == ("2", "1")
    dorsal = row(
        table="R17a",
        target_class="dorsal",
        trigger="-i/-e",
        subset="all",
        **base,
    )
    assert (
        dorsal["e"] == "1"
    )  # alunică/alunele (-ică suppletion), not in 'core'
    core = row(
        table="R17a",
        target_class="dorsal",
        trigger="-i/-e",
        subset="core",
        **base,
    )
    assert core["e"] == "0"
    coronal_e = row(
        table="R17c",
        target_class="coronal",
        trigger="-e",
        subset="all",
        **base,
    )
    # The non-trigger cell: N = e, except muscă/muște (sc → ște).
    assert int(coronal_e["N_trigger"]) - int(coronal_e["e"]) == 1


def test_rerun_is_skipped_and_force_reruns(
    tree: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    ctx = run_all(tree)
    capsys.readouterr()
    assert st.stage_tp(ctx) == 0
    assert "up to date" in capsys.readouterr().out
    forced = RunContext(
        RO_INFO, root=tree, params=PARAMS, force=frozenset({"tp"})
    )
    sha = digest(st.paths(ctx).tp_long)
    assert st.stage_tp(forced) == 0
    assert "up to date" not in capsys.readouterr().out
    assert digest(st.paths(ctx).tp_long) == sha


def test_readme_section_is_merged(tree: Path) -> None:
    out = tree / "analysis" / "comparison"
    out.mkdir(parents=True)
    other = (
        "# Comparison\n\n<!-- BEGIN italian -->\nItalian text\n"
        "<!-- END italian -->\n"
    )
    (out / "README.md").write_text(other, encoding="utf-8")
    ctx = run_all(tree)
    text = (out / "README.md").read_text(encoding="utf-8")
    assert text.startswith(other)
    assert text.count(readme.BEGIN) == 1 and text.count(readme.END) == 1
    assert "## Romanian (harmonized, read-only)" in text
    # A second run replaces the section in place.
    forced = RunContext(
        RO_INFO, root=tree, params=PARAMS, force=frozenset({"report"})
    )
    assert st.stage_report(forced) == 0
    assert (out / "README.md").read_text(encoding="utf-8") == text
    assert ctx.root == tree.resolve()


def test_merge_section_cases() -> None:
    created = readme.merge_section(None, "body")
    assert created.startswith(readme.TITLE) and created.endswith(
        readme.END + "\n"
    )
    replaced = readme.merge_section(created.replace("body", "old"), "new")
    assert "new" in replaced and "old" not in replaced
    appended = readme.merge_section("# Other\n\ntext\n", "body")
    assert appended.startswith("# Other\n\ntext\n\n" + readme.BEGIN)
    middle = "a\n" + readme.BEGIN + "\nold\n" + readme.END + "\nz\n"
    assert (
        readme.merge_section(middle, "new")
        == "a\n" + readme.BEGIN + "\nnew\n" + readme.END + "\nz\n"
    )


def test_roots_file_is_stamped_and_sorted(tree: Path) -> None:
    ctx = run_all(tree)
    roots = read_roots(st.paths(ctx).roots)
    ids = [r.root_id for r in roots.records]
    assert ids == sorted(ids) and len(ids) == len(set(ids))
    assert set(roots.inputs) == {
        "data/romanian_lexicon_with_freq.csv",
        "data/verb_paradigm_verdict.csv",
        "data/leipzig/freq/ron_news_2024_1M_freq.csv",
    }
    by_lemma = {r.lemma: r for r in roots.records}
    assert "fizic" not in by_lemma  # adjectives are not roots
    assert by_lemma["esperanto"].duplicates == 2
    assert by_lemma["esperanto"].plural == "esperanto"
    assert by_lemma["urs"].plural_count == 7  # urși + cedilla urşi
    assert (
        by_lemma["cruce"].nde == "frontstem"
        and by_lemma["cruce"].target_class == ""
    )


def test_cli_lists_and_runs(
    tree: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    assert cli.main(["--list"]) == 0
    assert "parity" in capsys.readouterr().out
    status = cli.main(
        [
            "--root",
            str(tree),
            "--learner-samples",
            "10",
            "--learner-V",
            "5",
            "roots",
            "tp",
        ]
    )
    assert status == 0
    assert (
        tree / "analysis" / "comparison" / "romanian_harmonized_tp_long.csv"
    ).is_file()
    with pytest.raises(SystemExit):
        cli.main(["--root", str(tree), "nosuchstage"])
