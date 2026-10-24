"""The 'ro' run context may write only under analysis/comparison/ and
data/romanian_compare/ (design §7; lexpipe.context.WRITE_ONLY_UNDER).

Every other guard of RunContext.check_writable still applies, and the
guard does not change what other languages may write.
"""

from __future__ import annotations

import dataclasses
import os
from pathlib import Path

import pytest

from lexpipe import context, csvio, licensing, spec
from lexpipe.context import RunContext
from lexpipe.langs.romanian_compare import PLUGIN, RO_INFO
from lexpipe.langs.romanian_compare import stages as st

pytestmark = pytest.mark.ro

ALLOWED = [
    "analysis/comparison/romanian_harmonized_tp_long.csv",
    "analysis/comparison/README.md",
    "analysis/comparison/sub/dir/file.txt",
    "data/romanian_compare/romanian_compare_roots.jsonl",
    "data/romanian_compare/manifests/tp.manifest.tsv",
]
REFUSED = [
    # Romanian tracked outputs and inputs (R1)
    "analysis/romanian_tp_summary_full.csv",
    "analysis/romanian_tp_summary_downsampled.csv",
    "analysis/romanian_new_table.csv",
    "analysis/palatalization_mismatches.csv",
    "analysis/learnability_checks/check_claims_output.txt",
    "data/romanian_lexicon_with_freq.csv",
    "data/romanian_lexicon_raw.csv",
    "data/romanian_anything.csv",
    "data/verb_paradigm_verdict.csv",
    # Leipzig (R2) and sources
    "data/leipzig/freq/ron_news_2024_1M_freq.csv",
    "data/leipzig/freq/new_freq.csv",
    "data/leipzig/corpora/x.txt",
    "data/sources/kaikki/x.jsonl",
    # other languages and the plugin's own analysis_dir
    "analysis/italian/tp_long.csv",
    "analysis/czech/tp_long.csv",
    "analysis/romanian_compare/tp_long.csv",
    "data/italian/italian_lexicon_raw.csv",
    "data/czech/nc/worklist.csv",
    # everything else in the repository
    "curated/romanian_compare/overrides.csv",
    "docs/romanian_compare/data_dict.md",
    "docs/runbooks.md",
    "README.md",
    "src/lexpipe/tp.py",
    "tests/x.py",
    "analysis/comparison/../romanian_tp_summary_full.csv",
    "data/romanian_compare/../romanian_lexicon_with_freq.csv",
]


@pytest.fixture
def ctx(tmp_path: Path) -> RunContext:
    return RunContext(RO_INFO, root=tmp_path)


def test_registry_row_and_plugin() -> None:
    assert RO_INFO.code == "ro" and RO_INFO.name == "romanian_compare"
    assert RO_INFO.module == "lexpipe.langs.romanian_compare"
    assert (PLUGIN.code, PLUGIN.name) == (RO_INFO.code, RO_INFO.name)
    assert list(PLUGIN.stages) == list(st.STAGES)
    assert context.WRITE_ONLY_UNDER == {
        "ro": ("analysis/comparison", "data/romanian_compare")
    }
    # RO_INFO is the registry's comparison-only row, so registry-wide
    # tests exempt it (no normalizer, docs, curated files or raw
    # lexicon) and the context confines it; a hand-made copy without the
    # flag is still confined by its code (WRITE_ONLY_UNDER).
    assert RO_INFO is spec.LANGUAGES["ro"] and RO_INFO.comparison_only
    unflagged = dataclasses.replace(RO_INFO, comparison_only=False)
    assert context.write_confinements(unflagged) == (
        ("analysis/comparison", "data/romanian_compare"),
    )
    assert spec.load_plugin("ro") is PLUGIN


def test_paths_point_into_the_allowed_directories(ctx: RunContext) -> None:
    assert ctx.data_dir == ctx.root / "data" / "romanian_compare"
    assert (
        ctx.manifest_dir
        == ctx.root / "data" / "romanian_compare" / "manifests"
    )
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
        assert ctx.check_writable(path) == path.resolve()
    for path in (
        p.decisions,
        p.lexicon,
        p.verdicts,
        p.leipzig,
        p.legacy_full,
        p.legacy_downsampled,
        p.check_claims,
    ):
        with pytest.raises(PermissionError):
            ctx.check_writable(path)


@pytest.mark.parametrize("relative", ALLOWED)
def test_allowed(ctx: RunContext, relative: str) -> None:
    path = ctx.write_text(ctx.root / relative, "x\n")
    assert path.read_text(encoding="utf-8") == "x\n"


@pytest.mark.parametrize("relative", REFUSED)
def test_refused(ctx: RunContext, relative: str) -> None:
    target = ctx.root / relative
    with pytest.raises(PermissionError):
        ctx.check_writable(target)
    with pytest.raises(PermissionError):
        ctx.write_text(target, "x")
    with pytest.raises(PermissionError):
        ctx.write_bytes(target, b"x")
    with pytest.raises(PermissionError):
        with ctx.open_output(target) as fh:
            fh.write("x")
    assert not target.resolve().exists()


def test_symlink_out_of_the_allowed_tree_is_refused(ctx: RunContext) -> None:
    outside = ctx.root / "analysis" / "romanian_tp_summary_full.csv"
    outside.parent.mkdir(parents=True)
    outside.write_text("legacy\n", encoding="utf-8")
    link_dir = ctx.root / "analysis" / "comparison"
    link_dir.mkdir(parents=True)
    link = link_dir / "innocent.csv"
    os.symlink(outside, link)
    with pytest.raises(PermissionError):
        ctx.write_text(link, "overwritten\n")
    assert outside.read_text(encoding="utf-8") == "legacy\n"


def test_stage_writers_are_guarded(ctx: RunContext) -> None:
    rows = [{"a": 1, "source": "wiktionary_en"}]
    with pytest.raises(PermissionError):
        st.write_tracked(
            ctx,
            ctx.root / "analysis" / "romanian_x.csv",
            rows,
            ["a", "source"],
        )
    with pytest.raises(PermissionError):
        st.write_untracked(
            ctx, ctx.root / "data" / "romanian_x.csv", rows, ["a", "source"]
        )
    # The licence gate still applies inside the allowed tree (R9).
    with pytest.raises(licensing.LicenceError):
        st.write_tracked(
            ctx,
            ctx.comparison_dir / "x.csv",
            [{"a": 1, "source": "leipzig"}],
            ["a", "source"],
        )
    st.write_tracked(ctx, ctx.comparison_dir / "x.csv", rows, ["a", "source"])
    assert csvio.read_rows(ctx.comparison_dir / "x.csv") == [
        {"a": "1", "source": "wiktionary_en"}
    ]


def test_other_languages_are_not_confined(tmp_path: Path) -> None:
    italian = RunContext(spec.LANGUAGES["it"], root=tmp_path)
    for relative in (
        "analysis/italian/x.csv",
        "data/italian/x.csv",
        "docs/italian/x.md",
    ):
        italian.write_text(tmp_path / relative, "x")
    for relative in ("analysis/romanian_x.csv", "data/leipzig/x.csv"):
        with pytest.raises(PermissionError):
            italian.check_writable(tmp_path / relative)
