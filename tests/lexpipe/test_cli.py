"""Tests for scripts/lang_pipeline.py with a fake plugin.

``spec.load_plugin`` and ``sources.sources_for`` are monkeypatched,
so nothing here depends on the real plugins or touches the network.
"""

from __future__ import annotations

import importlib.util
import sys
from collections.abc import Callable, Iterator, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from types import ModuleType

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lexpipe import context, sources, spec  # noqa: E402
from lexpipe.context import RunContext  # noqa: E402


@pytest.fixture(scope="module")
def cli() -> Iterator[ModuleType]:
    name = "_lexpipe_test_lang_pipeline"
    path = PROJECT_ROOT / "scripts" / "lang_pipeline.py"
    module_spec = importlib.util.spec_from_file_location(name, path)
    assert module_spec is not None and module_spec.loader is not None
    module = importlib.util.module_from_spec(module_spec)
    sys.modules[name] = module
    try:
        module_spec.loader.exec_module(module)
        yield module
    finally:
        sys.modules.pop(name, None)


@dataclass
class Recorder:
    """Stage calls, with the context each stage received."""

    calls: list[str] = field(default_factory=list)
    contexts: list[RunContext] = field(default_factory=list)

    def stage(
        self, name: str, status: int | None = None
    ) -> Callable[[RunContext], int | None]:
        def run(ctx: RunContext) -> int | None:
            self.calls.append(name)
            self.contexts.append(ctx)
            return status

        return run


@pytest.fixture
def calls(monkeypatch: pytest.MonkeyPatch) -> Recorder:
    """Install a fake 'it' plugin whose stages record their calls."""
    record = Recorder()
    plugin = spec.PluginSpec(
        code="it",
        name="italian",
        stages={
            "build": record.stage("build"),
            "check": record.stage("check", 0),
            "analyse": record.stage("analyse"),
        },
    )
    monkeypatch.setattr(spec, "load_plugin", lambda code: plugin)
    monkeypatch.setattr(sources, "sources_for", lambda code: ())
    return record


def test_list(
    cli: ModuleType, calls: Recorder, capsys: pytest.CaptureFixture[str]
) -> None:
    assert cli.main(["--lang", "it", "--list"]) == 0
    lines = capsys.readouterr().out.splitlines()
    assert lines == [
        "Stages for italian (it), in pipeline order:",
        "  fetch  (common)",
        "  verify-sources  (common)",
        "  build",
        "  check",
        "  analyse",
    ]
    assert calls.calls == []


def test_named_stages_run_in_pipeline_order(
    cli: ModuleType, calls: Recorder, capsys: pytest.CaptureFixture[str]
) -> None:
    assert cli.main(["--lang", "it", "analyse", "build"]) == 0
    assert calls.calls == ["build", "analyse"]
    out = capsys.readouterr().out
    assert "[it] build: done in" in out
    assert "[it] stage timings:" in out and "total" in out


def test_default_runs_everything(cli: ModuleType, calls: Recorder) -> None:
    assert cli.main(["--lang", "it"]) == 0
    assert calls.calls == ["build", "check", "analyse"]


def test_stages_receive_one_run_context(
    cli: ModuleType, calls: Recorder, tmp_path: Path
) -> None:
    argv = ["--lang", "it", "--root", str(tmp_path), "--seed", "7"]
    assert cli.main([*argv, "--force", "check", "build", "check"]) == 0
    assert calls.calls == ["build", "check"]
    first, second = calls.contexts
    assert first is second
    assert first.language == spec.LANGUAGES["it"]
    assert first.root == tmp_path.resolve()
    assert first.seed == 7
    assert first.force == frozenset({"check"})
    assert first.forced("check") and not first.forced("build")
    assert first.data_dir == tmp_path.resolve() / "data" / "italian"


def test_defaults_of_the_context(cli: ModuleType, calls: Recorder) -> None:
    assert cli.main(["--lang", "it", "build"]) == 0
    (ctx,) = calls.contexts
    assert ctx.root == PROJECT_ROOT
    assert ctx.seed == context.DEFAULT_SEED
    assert ctx.force == frozenset()


@pytest.mark.parametrize(
    "argv",
    [
        ["--force", "fetch"],  # common stages manage data/sources/
        ["--force", "verify-sources", "verify-sources"],
        ["--force", "nonesuch"],
        ["--force", "check", "build"],  # forced stage is not run
    ],
)
def test_bad_force(cli: ModuleType, calls: Recorder, argv: list[str]) -> None:
    with pytest.raises(SystemExit) as excinfo:
        cli.main(["--lang", "it", *argv])
    assert excinfo.value.code == 2
    assert calls.calls == []


def test_forcing_a_common_stage_is_refused_before_the_plugin_loads(
    cli: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    def forbidden(code: str) -> spec.LanguagePlugin:
        raise AssertionError("the plugin must not be imported")

    monkeypatch.setattr(spec, "load_plugin", forbidden)
    with pytest.raises(SystemExit) as excinfo:
        cli.main(["--lang", "cs", "--force", "fetch"])
    assert excinfo.value.code == 2


def test_root_must_exist(
    cli: ModuleType, calls: Recorder, tmp_path: Path
) -> None:
    with pytest.raises(SystemExit) as excinfo:
        cli.main(["--lang", "it", "--root", str(tmp_path / "nope")])
    assert excinfo.value.code == 2


def test_failed_stage_stops_the_run(
    cli: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    record: list[str] = []

    def failing(ctx: RunContext) -> int:
        record.append("build")
        return 3

    def later(ctx: RunContext) -> None:
        record.append("analyse")

    plugin = spec.PluginSpec(
        code="it", name="italian", stages={"build": failing, "analyse": later}
    )
    monkeypatch.setattr(spec, "load_plugin", lambda code: plugin)
    assert cli.main(["--lang", "it", "build", "analyse"]) == 1
    assert record == ["build"]


def test_unknown_stage(cli: ModuleType, calls: Recorder) -> None:
    with pytest.raises(SystemExit) as excinfo:
        cli.main(["--lang", "it", "nonesuch"])
    assert excinfo.value.code == 2


def test_missing_plugin(
    cli: ModuleType,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    def missing(code: str) -> spec.LanguagePlugin:
        raise spec.PluginNotFoundError("the czech plugin is not installed")

    monkeypatch.setattr(spec, "load_plugin", missing)
    assert cli.main(["--lang", "cs"]) == 1
    assert "czech plugin is not installed" in capsys.readouterr().err
    assert cli.main(["--lang", "cs", "--list"]) == 1
    out = capsys.readouterr().out
    assert "fetch  (common)" in out and "verify-sources" in out


def test_common_stages_need_no_plugin(
    cli: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    def forbidden(code: str) -> spec.LanguagePlugin:
        raise AssertionError("the plugin must not be imported")

    monkeypatch.setattr(spec, "load_plugin", forbidden)
    monkeypatch.setattr(sources, "sources_for", lambda code: ())
    assert cli.main(["--lang", "cs", "verify-sources", "fetch"]) == 0


# ---------------------------------------------------------------------------
# verify-sources with a rolling source, network stubbed
# ---------------------------------------------------------------------------

PAYLOAD = b"kaikki" * 100
ROLLING = sources.Source(
    name="toy_rolling",
    url="https://example.org/raw.jsonl.gz",
    path="data/sources/toy/raw.jsonl.gz",
    size=len(PAYLOAD),
    licence="CC0",
    citation="none",
    languages=("it",),
    rolling=True,
)


@dataclass
class FakeHeadResponse:
    status_code: int
    headers: Mapping[str, str]

    def close(self) -> None:
        pass


@dataclass
class FakeHeadSession:
    etag: str
    size: int = len(PAYLOAD)
    urls: list[str] = field(default_factory=list)

    def head(
        self,
        url: str,
        *,
        headers: Mapping[str, str],
        timeout: tuple[float, float],
        allow_redirects: bool,
    ) -> FakeHeadResponse:
        self.urls.append(url)
        return FakeHeadResponse(
            200,
            {
                "etag": self.etag,  # lowercase on purpose
                "Last-Modified": "Wed, 21 Oct 2026 08:00:00 GMT",
                "Content-Length": str(self.size),
            },
        )


def test_verify_sources_records_upstream_metadata(
    cli: ModuleType,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.setattr(sources, "sources_for", lambda code: (ROLLING,))
    target = sources.local_path(ROLLING, tmp_path)
    target.parent.mkdir(parents=True)
    target.write_bytes(PAYLOAD)
    ctx = RunContext(spec.LANGUAGES["it"], root=tmp_path)

    session = FakeHeadSession('"abc-1"')
    assert cli.verify_sources(ctx, session=session) == 0
    assert session.urls == [ROLLING.url]
    entry = sources.read_lock(ctx.lock_path)["toy_rolling"]
    assert entry.etag == '"abc-1"'
    assert entry.last_modified == "Wed, 21 Oct 2026 08:00:00 GMT"
    assert "recorded" in capsys.readouterr().out

    assert cli.verify_sources(ctx, session=FakeHeadSession('"abc-1"')) == 0
    assert "upstream unchanged" in capsys.readouterr().out

    # kaikki regenerated the file: still ok locally, but reported.
    assert cli.verify_sources(ctx, session=FakeHeadSession('"abc-2"')) == 0
    assert "upstream changed" in capsys.readouterr().out
    assert sources.read_lock(ctx.lock_path)["toy_rolling"].etag == '"abc-1"'


# ---------------------------------------------------------------------------
# --param
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("200", 200),
        ("-3", -3),
        ("0.5", 0.5),
        ("1e-3", 0.001),
        (".5", 0.5),
        ("True", True),
        ("false", False),
        ("[500, 1000]", [500, 1000]),
        ('["a", 1.5, true]', ["a", 1.5, True]),
        ("velar", "velar"),
        ("nan", "nan"),  # non-finite numbers stay strings
        ("inf", "inf"),
        ("1_000", "1_000"),
        ("", ""),
        ("a=b", "a=b"),
    ],
)
def test_param_values_are_typed(
    cli: ModuleType, text: str, expected: object
) -> None:
    value = cli.parse_param_value(text)
    assert value == expected and type(value) is type(expected)


@pytest.mark.parametrize("text", ["[500,", "[[1]]", '[{"a": 1}]', "[NaN]"])
def test_bad_param_lists(cli: ModuleType, text: str) -> None:
    with pytest.raises(ValueError):
        cli.parse_param_value(text)


def test_params_reach_the_run_context(
    cli: ModuleType, calls: Recorder
) -> None:
    argv = [
        "--lang",
        "it",
        "--param",
        "learner_samples=40",
        "--param",
        "learner_V=[5,10]",
        "--param",
        "verify_source_size=false",
        "--param",
        "label=a b",
        "build",
    ]
    assert cli.main(argv) == 0
    (ctx,) = calls.contexts
    assert dict(ctx.params) == {
        "learner_samples": 40,
        "learner_V": [5, 10],
        "verify_source_size": False,
        "label": "a b",
    }
    # Every value is a valid manifest parameter.
    for value in ctx.params.values():
        context.canonical_param(value)


def test_no_params_by_default(cli: ModuleType, calls: Recorder) -> None:
    assert cli.main(["--lang", "it", "build"]) == 0
    assert dict(calls.contexts[0].params) == {}


@pytest.mark.parametrize(
    "param",
    ["noequals", "=5", "bad key=1", "a\tb=1", "[x]=1"],
)
def test_bad_params(cli: ModuleType, calls: Recorder, param: str) -> None:
    with pytest.raises(SystemExit) as excinfo:
        cli.main(["--lang", "it", "--param", param, "build"])
    assert excinfo.value.code == 2
    assert calls.calls == []


def test_repeated_param(cli: ModuleType, calls: Recorder) -> None:
    with pytest.raises(SystemExit) as excinfo:
        cli.main(["--lang", "it", "--param", "a=1", "--param", "a=2"])
    assert excinfo.value.code == 2
    assert calls.calls == []


# ---------------------------------------------------------------------------
# Sub-stages
# ---------------------------------------------------------------------------


@pytest.fixture
def nested(monkeypatch: pytest.MonkeyPatch) -> Recorder:
    """A fake 'cs' plugin whose qc stage runs a qc-morfflex sub-stage."""
    record = Recorder()
    plugin = spec.PluginSpec(
        code="cs",
        name="czech",
        stages={
            "lexicon": record.stage("lexicon"),
            "qc": record.stage("qc"),
            "tp": record.stage("tp"),
        },
        substages={"qc-morfflex": "qc"},
    )
    monkeypatch.setattr(spec, "load_plugin", lambda code: plugin)
    monkeypatch.setattr(sources, "sources_for", lambda code: ())
    return record


def test_force_a_substage(cli: ModuleType, nested: Recorder) -> None:
    assert cli.main(["--lang", "cs", "--force", "qc-morfflex", "qc"]) == 0
    assert nested.calls == ["qc"]
    (ctx,) = nested.contexts
    assert ctx.forced("qc-morfflex") and not ctx.forced("qc")
    assert cli.main(["--lang", "cs", "--force", "qc-morfflex"]) == 0
    assert nested.calls == ["qc", "lexicon", "qc", "tp"]


def test_substage_of_a_stage_that_is_not_run(
    cli: ModuleType, nested: Recorder, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(SystemExit) as excinfo:
        cli.main(["--lang", "cs", "--force", "qc-morfflex", "tp"])
    assert excinfo.value.code == 2
    assert "qc-morfflex (sub-stage of qc)" in capsys.readouterr().err
    assert nested.calls == []


def test_substages_are_not_run_on_their_own(
    cli: ModuleType, nested: Recorder, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(SystemExit) as excinfo:
        cli.main(["--lang", "cs", "qc-morfflex"])
    assert excinfo.value.code == 2
    assert "sub-stage of qc" in capsys.readouterr().err
    assert nested.calls == []


def test_list_shows_substages(
    cli: ModuleType, nested: Recorder, capsys: pytest.CaptureFixture[str]
) -> None:
    assert cli.main(["--lang", "cs", "--list"]) == 0
    assert capsys.readouterr().out.splitlines() == [
        "Stages for czech (cs), in pipeline order:",
        "  fetch  (common)",
        "  verify-sources  (common)",
        "  lexicon",
        "  qc",
        "    qc-morfflex  (sub-stage of qc; --force only)",
        "  tp",
    ]


# ---------------------------------------------------------------------------
# Derived sources and the Romanian comparison
# ---------------------------------------------------------------------------

SUBSET = sources.Source(
    name="toy_subset",
    url="",
    path="data/sources/toy/subset.jsonl",
    size=10,
    licence="CC0",
    citation="none",
    languages=("it",),
    sha256="0" * 64,
    kind=sources.SourceKind.DERIVED,
    derived_from="toy_rolling",
    derived_from_sha256="1" * 64,
    recipe="python make_subset.py it",
)


def test_fetch_reports_how_to_make_a_derived_source(
    cli: ModuleType,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.setattr(sources, "sources_for", lambda code: (SUBSET,))
    ctx = RunContext(spec.LANGUAGES["it"], root=tmp_path)
    assert cli.fetch_sources(ctx) == 1
    out = capsys.readouterr().out
    assert "ERROR toy_subset" in out and "python make_subset.py it" in out
    assert not sources.local_path(SUBSET, tmp_path).exists()


def test_romanian_comparison_is_a_cli_language(
    cli: ModuleType, capsys: pytest.CaptureFixture[str]
) -> None:
    assert "ro" in cli.build_parser().parse_args(["--lang", "ro"]).lang
    assert cli.main(["--lang", "ro", "--list"]) == 0
    out = capsys.readouterr().out.splitlines()
    assert out[0] == "Stages for romanian_compare (ro), in pipeline order:"
    assert "  tp" in out and "  learners" in out
