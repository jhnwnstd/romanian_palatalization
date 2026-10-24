"""Tests for lexpipe.context: paths, guarded atomic writes, manifests.

Everything runs in a temporary root, as ``--root DIR`` would.
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lexpipe import context, licensing, spec  # noqa: E402
from lexpipe.context import RunContext  # noqa: E402

ITALIAN = spec.LANGUAGES["it"]


@pytest.fixture
def ctx(tmp_path: Path) -> RunContext:
    return RunContext(ITALIAN, root=tmp_path, params={"window": 3})


def test_paths(ctx: RunContext, tmp_path: Path) -> None:
    root = tmp_path.resolve()
    assert ctx.root == root
    assert ctx.data_dir == root / "data" / "italian"
    assert ctx.nc_dir == root / "data" / "italian" / "nc"
    assert ctx.manifest_dir == root / "data" / "italian" / "manifests"
    assert ctx.analysis_dir == root / "analysis" / "italian"
    assert ctx.comparison_dir == root / "analysis" / "comparison"
    assert ctx.curated_dir == root / "curated" / "italian"
    assert ctx.sources_dir == root / "data" / "sources"
    assert ctx.docs_dir == root / "docs" / "italian"
    assert ctx.lock_path == root / "docs" / "sources_lock.tsv"
    assert ctx.manifest_path("tp") == ctx.manifest_dir / "tp.manifest.tsv"
    assert ctx.relative(ctx.data_dir / "x.csv") == "data/italian/x.csv"
    assert ctx.relative(Path("/elsewhere/x")) == "/elsewhere/x"


def test_context_is_frozen(tmp_path: Path) -> None:
    params = {"window": 3}
    ctx = RunContext(ITALIAN, root=tmp_path, force={"tp"}, params=params)
    params["window"] = 4
    assert ctx.params["window"] == 3
    with pytest.raises(TypeError):
        ctx.params["window"] = 5  # type: ignore[index]
    assert ctx.force == frozenset({"tp"}) and ctx.forced("tp")
    assert ctx.seed == context.DEFAULT_SEED == 20080715
    with pytest.raises(TypeError, match="seed"):
        RunContext(ITALIAN, root=tmp_path, seed="7")  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "relative",
    [
        "data/sources/kaikki/x.jsonl",
        "data/sources",
        "curated/italian/overrides.csv",
        "data/leipzig/freq/ita_news_2024_1M_freq.csv",
        "data/romanian_lexicon_with_freq.csv",
        "analysis/romanian_tp_summary_full.csv",
    ],
)
def test_protected_locations(ctx: RunContext, relative: str) -> None:
    with pytest.raises(PermissionError):
        ctx.check_writable(ctx.root / relative)
    with pytest.raises(PermissionError):
        ctx.write_text(ctx.root / relative, "x")
    assert not (ctx.root / relative).is_file()


@pytest.mark.parametrize(
    "relative",
    [
        "data/italian/italian_lexicon_raw.csv",
        "data/italian/nc/worklist.csv",
        "analysis/italian/tp_long.csv",
        "analysis/comparison/tp_long.csv",
        "docs/italian/data_dict.md",
    ],
)
def test_writable_locations(ctx: RunContext, relative: str) -> None:
    path = ctx.write_text(ctx.root / relative, "a,b\n")
    assert path.read_text(encoding="utf-8") == "a,b\n"


def test_atomic_writes(tmp_path: Path) -> None:
    target = tmp_path / "out" / "table.csv"
    context.atomic_write_text(target, "first\r\n")
    assert target.read_bytes() == b"first\r\n"  # line endings untouched
    with pytest.raises(RuntimeError):
        with context.atomic_open(target) as fh:
            fh.write("half")
            raise RuntimeError("stage crashed")
    assert target.read_bytes() == b"first\r\n"
    assert not [p for p in target.parent.iterdir() if p.suffix == ".tmp"]
    context.atomic_write_bytes(target, b"\x00\x01")
    assert target.read_bytes() == b"\x00\x01"
    with pytest.raises(ValueError, match="writes only"):
        with context.atomic_open(target, "a"):
            pass


def test_canonical_params() -> None:
    assert context.canonical_param({"b": [1, 2], "a": "è"}) == (
        '{"a":"è","b":[1,2]}'
    )
    with pytest.raises(TypeError):
        context.canonical_param({1, 2})
    with pytest.raises(ValueError):
        context.canonical_param(float("nan"))


def test_git_state(tmp_path: Path) -> None:
    assert context.git_state(tmp_path) == context.GitState(None, None)
    if shutil.which("git") is None:
        pytest.skip("git is not installed")
    state = context.git_state(PROJECT_ROOT)
    if state.sha is None:
        pytest.skip("not a git checkout")
    assert re.fullmatch(r"[0-9a-f]{40}", state.sha)
    assert state.dirty in (True, False)


# ---------------------------------------------------------------------------
# Manifests and skipping
# ---------------------------------------------------------------------------


class Stage:
    """A stage that copies its input upper-cased; counts its runs."""

    def __init__(self, ctx: RunContext, source: Path, target: Path) -> None:
        self.ctx, self.source, self.target = ctx, source, target
        self.runs = 0

    def body(self) -> int:
        self.runs += 1
        text = self.source.read_text(encoding="utf-8").upper()
        self.ctx.write_text(self.target, text)
        return 0

    def run(self, ctx: RunContext | None = None, **params: object) -> object:
        return context.run_stage(
            ctx or self.ctx,
            "upper",
            inputs=[self.source],
            params=params or {"version": 1},
            outputs=[self.target],
            body=self.body,
            log=lambda _: None,
        )


@pytest.fixture
def stage(ctx: RunContext) -> Stage:
    source = ctx.sources_dir / "toy" / "in.txt"
    source.parent.mkdir(parents=True)
    source.write_text("amico\n", encoding="utf-8")
    return Stage(ctx, source, ctx.data_dir / "upper.txt")


def test_manifest_round_trip(stage: Stage, ctx: RunContext) -> None:
    assert stage.run() == 0 and stage.runs == 1
    manifest = context.read_manifest(ctx.manifest_path("upper"))
    assert manifest is not None
    assert manifest.stage == "upper"
    assert dict(manifest.params) == {"version": "1"}
    assert list(manifest.inputs) == ["data/sources/toy/in.txt"]
    assert list(manifest.outputs) == ["data/italian/upper.txt"]
    assert manifest.outputs["data/italian/upper.txt"] == (
        context.sha256_file(stage.target)
    )
    assert manifest.git == context.GitState(None, None)  # not a repo
    lines = ctx.manifest_path("upper").read_text(encoding="utf-8")
    assert lines.splitlines()[:4] == [
        "kind\tkey\tvalue",
        "stage\tname\tupper",
        "git\tdirty\t",
        "git\tsha\t",
    ]


def test_skip_when_up_to_date(stage: Stage, ctx: RunContext) -> None:
    stage.run()
    assert context.should_skip(ctx, "upper", [stage.source], {"version": 1})
    assert stage.run() == 0 and stage.runs == 1


@pytest.mark.parametrize(
    "change",
    ["params", "input", "output", "missing-output", "force", "manifest"],
)
def test_rerun_when_anything_changed(
    stage: Stage, ctx: RunContext, change: str
) -> None:
    stage.run()
    run_ctx = ctx
    params: dict[str, object] = {}
    if change == "params":
        params = {"version": 2}
    elif change == "input":
        stage.source.write_text("nemico\n", encoding="utf-8")
    elif change == "output":
        stage.target.write_text("edited by hand\n", encoding="utf-8")
    elif change == "missing-output":
        stage.target.unlink()
    elif change == "force":
        run_ctx = RunContext(ITALIAN, root=ctx.root, force={"upper"})
    else:
        ctx.manifest_path("upper").write_text("garbage\n", encoding="utf-8")
    stage.run(run_ctx, **params)
    assert stage.runs == 2


def test_failed_stage_writes_no_manifest(
    ctx: RunContext, stage: Stage
) -> None:
    status = context.run_stage(
        ctx,
        "upper",
        inputs=[stage.source],
        params={},
        outputs=[],
        body=lambda: 2,
        log=lambda _: None,
    )
    assert status == 2
    assert not ctx.manifest_path("upper").exists()


def test_declared_outputs_must_exist(ctx: RunContext, stage: Stage) -> None:
    with pytest.raises(FileNotFoundError, match="did not write"):
        context.run_stage(
            ctx,
            "upper",
            inputs=[stage.source],
            params={},
            outputs=[ctx.data_dir / "never.csv"],
            body=lambda: None,
            log=lambda _: None,
        )


def test_missing_input(ctx: RunContext) -> None:
    with pytest.raises(FileNotFoundError):
        context.should_skip(ctx, "upper", [ctx.root / "nope"], {})


def test_manifest_values_are_written_verbatim(
    stage: Stage, ctx: RunContext
) -> None:
    stage.run(ctx, version=1, label="dorsal\tvelar", cells=["NOM.PL"])
    text = ctx.manifest_path("upper").read_text(encoding="utf-8")
    assert 'param\tcells\t["NOM.PL"]\n' in text
    assert 'param\tlabel\t"dorsal\\tvelar"\n' in text  # JSON escapes tabs
    manifest = context.read_manifest(ctx.manifest_path("upper"))
    assert manifest is not None
    assert manifest.params["label"] == '"dorsal\\tvelar"'
    assert context.should_skip(
        ctx,
        "upper",
        [stage.source],
        {"version": 1, "label": "dorsal\tvelar", "cells": ["NOM.PL"]},
    )


# ---------------------------------------------------------------------------
# Comparison-only languages
# ---------------------------------------------------------------------------

ROMANIAN = spec.LANGUAGES["ro"]
# A comparison-only language the registry does not know.
TOY_COMPARE = spec.LanguageInfo(
    "zz", "toy_compare", "lexpipe_toy", comparison_only=True
)
# A hand-made "ro" row without the flag (like the plugin's RO_INFO).
RO_UNFLAGGED = spec.LanguageInfo("ro", "romanian_compare", "m")


def test_write_only_under_is_derived_from_the_registry() -> None:
    assert dict(context.WRITE_ONLY_UNDER) == {
        "ro": ("analysis/comparison", "data/romanian_compare")
    }
    assert context.comparison_write_dirs(TOY_COMPARE) == (
        "analysis/comparison",
        "data/toy_compare",
    )
    assert context.write_confinements(ITALIAN) == ()
    assert context.write_confinements(ROMANIAN) == (
        ("analysis/comparison", "data/romanian_compare"),
    )
    assert context.write_confinements(RO_UNFLAGGED) == (
        ("analysis/comparison", "data/romanian_compare"),
    )
    assert context.write_confinements(TOY_COMPARE) == (
        ("analysis/comparison", "data/toy_compare"),
    )


@pytest.mark.parametrize("info", [ROMANIAN, RO_UNFLAGGED, TOY_COMPARE])
def test_comparison_only_contexts_are_confined(
    tmp_path: Path, info: spec.LanguageInfo
) -> None:
    ctx = RunContext(info, root=tmp_path)
    for relative in (
        "analysis/comparison/x.csv",
        f"data/{info.name}/x.jsonl",
        f"data/{info.name}/manifests/tp.manifest.tsv",
    ):
        ctx.write_text(tmp_path / relative, "x")
    for relative in (
        f"analysis/{info.name}/x.csv",
        f"docs/{info.name}/x.md",
        "analysis/italian/x.csv",
        "data/italian/x.csv",
        "README.md",
        "analysis/comparison",  # the directory itself
        "analysis/romanian_tp_summary_full.csv",  # R1
        "data/leipzig/freq/x.csv",  # R2
    ):
        with pytest.raises(PermissionError):
            ctx.write_text(tmp_path / relative, "x")
        assert not (tmp_path / relative).is_file()


def test_full_languages_are_not_confined(tmp_path: Path) -> None:
    ctx = RunContext(ITALIAN, root=tmp_path)
    for relative in ("docs/italian/x.md", "analysis/italian/x.csv"):
        ctx.write_text(tmp_path / relative, "x")


# ---------------------------------------------------------------------------
# write_tracked_text (R9)
# ---------------------------------------------------------------------------


def test_tracked_text_from_redistributable_sources(ctx: RunContext) -> None:
    path = ctx.analysis_dir / "summary.txt"
    written = ctx.write_tracked_text(
        path, "N = 12\n", sources=["wiktionary_en", "morphit"]
    )
    assert written == path.resolve()
    assert path.read_text(encoding="utf-8") == "N = 12\n"  # no header
    ctx.write_tracked_text(path, "N = 13\n", sources="wiktionary_it|cnk")
    assert path.read_text(encoding="utf-8") == "N = 13\n"


@pytest.mark.parametrize(
    "sources",
    [["morfflex"], ["wiktionary_cs", "derinet"], ["leipzig"], [], [""], ["x"]],
)
def test_tracked_text_refuses_without_the_nc_flag(
    ctx: RunContext, sources: list[str]
) -> None:
    path = ctx.analysis_dir / "summary.txt"
    with pytest.raises(licensing.LicenceError):
        ctx.write_tracked_text(path, "x\n", sources=sources)
    assert not path.exists()


def test_tracked_text_with_nc_aggregates(tmp_path: Path) -> None:
    ctx = RunContext(spec.LANGUAGES["cs"], root=tmp_path)
    path = ctx.analysis_dir / "czech_covariation.txt"
    ctx.write_tracked_text(
        path,
        "pairs: 120\n",
        sources=["wiktionary_cs", "derinet", "czetyl"],
        nc_aggregates=True,
    )
    text = path.read_text(encoding="utf-8")
    header = licensing.aggregate_header(["czetyl", "derinet", "wiktionary_cs"])
    assert text == header + "pairs: 120\n"
    assert text.startswith(licensing.AGGREGATE_HEADER_TITLE + "\n")
    assert "Aggregate counts only" in text
    assert "CC BY-NC-SA 4.0" in text and "DATA_LICENSE.md" in text


@pytest.mark.parametrize(
    "sources",
    [["leipzig"], ["derinet", "dex_qc"], ["unknown_source"], [], [""]],
)
def test_nc_aggregates_still_refuse_unverified_and_unknown(
    ctx: RunContext, sources: list[str]
) -> None:
    path = ctx.analysis_dir / "summary.txt"
    with pytest.raises(licensing.LicenceError):
        ctx.write_tracked_text(
            path, "x\n", sources=sources, nc_aggregates=True
        )
    assert not path.exists()


@pytest.mark.parametrize(
    "relative",
    [
        "data/italian/nc/report.txt",  # NC derivatives are untracked
        "data/italian/report.txt",
        "docs/report.md",
        "analysis/romanian/report.txt",  # no such language
        "analysis/romanian_compare/report.txt",  # comparison-only: none
    ],
)
def test_tracked_text_needs_a_tracked_location(
    ctx: RunContext, relative: str
) -> None:
    with pytest.raises(licensing.LicenceError, match="tracked location"):
        ctx.write_tracked_text(
            ctx.root / relative, "x\n", sources=["wiktionary_en"]
        )
    assert not (ctx.root / relative).exists()


def test_tracked_text_keeps_the_write_guards(ctx: RunContext) -> None:
    with pytest.raises(PermissionError):
        ctx.write_tracked_text(
            ctx.root / "analysis" / "romanian_new.txt",
            "x\n",
            sources=["wiktionary_en"],
        )


@pytest.mark.parametrize("info", [ROMANIAN, RO_UNFLAGGED, TOY_COMPARE])
def test_tracked_text_of_a_comparison_language(
    tmp_path: Path, info: spec.LanguageInfo
) -> None:
    ctx = RunContext(info, root=tmp_path)
    ctx.write_tracked_text(
        ctx.comparison_dir / "x.txt", "x\n", sources=["wiktionary_ro"]
    )
    # Its data directory is writable but never a tracked location.
    target = ctx.data_dir / info.file_name("lexicon_raw")
    with pytest.raises(licensing.LicenceError):
        ctx.write_tracked_text(target, "x\n", sources=["wiktionary_ro"])
    assert not target.exists()


def test_tracked_text_of_an_unregistered_full_language(
    tmp_path: Path,
) -> None:
    toy = spec.LanguageInfo("zz", "toy", "m", file_prefix="toy_")
    ctx = RunContext(toy, root=tmp_path)
    ctx.write_tracked_text(
        ctx.analysis_dir / "x.txt", "x\n", sources=["wiktionary_en"]
    )
    assert (tmp_path / "analysis" / "toy" / "x.txt").is_file()


def test_manifest_records_the_code_fingerprint(
    stage: Stage, ctx: RunContext
) -> None:
    stage.run()
    manifest = context.read_manifest(ctx.manifest_path("upper"))
    assert manifest is not None
    assert manifest.code == context.code_fingerprint(ctx.language)
    assert manifest.code is not None and len(manifest.code) == 64


def test_rerun_when_code_changed(
    stage: Stage, ctx: RunContext, monkeypatch: pytest.MonkeyPatch
) -> None:
    stage.run()
    assert stage.runs == 1
    monkeypatch.setattr(context, "code_fingerprint", lambda _: "0" * 64)
    assert not context.should_skip(
        ctx, "upper", [stage.source], {"version": 1}
    )
    stage.run()
    assert stage.runs == 2
    stage.run()  # the new fingerprint is recorded now
    assert stage.runs == 2


def test_manifest_without_code_row_is_stale(
    stage: Stage, ctx: RunContext
) -> None:
    stage.run()
    path = ctx.manifest_path("upper")
    kept = [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if not line.startswith("code\t")
    ]
    path.write_text("\n".join(kept) + "\n", encoding="utf-8")
    manifest = context.read_manifest(path)
    assert manifest is not None and manifest.code is None
    assert not context.should_skip(
        ctx, "upper", [stage.source], {"version": 1}
    )


def test_code_fingerprint_follows_the_plugin_source(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    package = tmp_path / "toy_plugin_pkg"
    package.mkdir()
    (package / "__init__.py").write_text("", encoding="utf-8")
    module = package / "stages.py"
    module.write_text("VERSION = 1\n", encoding="utf-8")
    monkeypatch.syspath_prepend(str(tmp_path))
    toy = spec.LanguageInfo("zz", "toy", "toy_plugin_pkg")
    first = context.code_fingerprint(toy)
    assert context.code_fingerprint(toy) == first
    module.write_text("VERSION = 2\n", encoding="utf-8")
    assert context.code_fingerprint(toy) != first
    missing = spec.LanguageInfo("zz", "toy", "no_such_plugin_pkg")
    assert context.code_fingerprint(missing) == context.code_fingerprint(
        missing
    )
