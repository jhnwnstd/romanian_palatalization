"""Tests for lexpipe.spec: the plugin contract and lazy registry.

Plugins are faked with throwaway packages written to tmp_path and
registered under test-only codes, so these tests do not depend on
whether lexpipe.langs.italian / lexpipe.langs.czech exist yet.
"""

from __future__ import annotations

import sys
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lexpipe import spec  # noqa: E402
from lexpipe.context import RunContext  # noqa: E402
from lexpipe.spec import (  # noqa: E402
    LanguageInfo,
    PluginNotFoundError,
    PluginSpec,
)

FAKE_PACKAGES = ("lexpipe_fake_ok", "lexpipe_fake_broken", "lexpipe_fake_bare")


def _noop(ctx: RunContext) -> None:
    return None


@pytest.fixture
def plugin_dir(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> Iterator[Path]:
    """A sys.path entry for fake plugin packages, cleaned up after."""
    monkeypatch.syspath_prepend(str(tmp_path))
    yield tmp_path
    for name in FAKE_PACKAGES:
        sys.modules.pop(name, None)


def _package(root: Path, name: str, source: str) -> None:
    (root / name).mkdir()
    (root / name / "__init__.py").write_text(source, encoding="utf-8")


def test_registry_rows() -> None:
    assert spec.LANGUAGES["it"] == LanguageInfo(
        code="it",
        name="italian",
        module="lexpipe.langs.italian",
        kaikki_code="it",
        kaikki_name="Italian",
        leipzig_id="ita_news_2024_1M",
        file_prefix="italian_",
        analysis_module="phonology.analyses.italian_palatalization",
        elision_clitics=spec.ITALIAN_ELISION_CLITICS,
    )
    assert spec.LANGUAGES["cs"] == LanguageInfo(
        code="cs",
        name="czech",
        module="lexpipe.langs.czech",
        kaikki_code="cs",
        kaikki_name="Czech",
        leipzig_id="ces_news_2024_1M",
        file_prefix="czech_",
        analysis_module="phonology.analyses.czech_palatalization",
    )
    info = spec.language_info("cs")
    assert info.data_dir() == PROJECT_ROOT / "data" / "czech"
    assert info.analysis_dir() == PROJECT_ROOT / "analysis" / "czech"
    assert info.file_name("lexicon_raw") == "czech_lexicon_raw.csv"
    assert info.file_name("plural_forms", ".tsv") == "czech_plural_forms.tsv"


def test_registry_is_consistent() -> None:
    for code, info in spec.LANGUAGES.items():
        assert info.code == info.kaikki_code == code
        assert info.module == f"lexpipe.langs.{info.name}"
        assert info.file_prefix == f"{info.name}_"
        assert info.leipzig_id.endswith("_news_2024_1M")
        if info.comparison_only:
            # Named after what it does (romanian_compare); kaikki_name
            # is still the language's.
            assert info.name.startswith(info.kaikki_name.lower() + "_")
            assert not info.analysis_module
        else:
            assert info.kaikki_name == info.name.capitalize()
    # Only Italian splits elided clitics.
    assert spec.LANGUAGES["it"].elision_clitics >= {"l", "dell", "un", "c"}
    assert not spec.LANGUAGES["cs"].elision_clitics
    assert not spec.LANGUAGES["ro"].elision_clitics
    assert all(
        clitic == clitic.lower() and "'" not in clitic
        for clitic in spec.ITALIAN_ELISION_CLITICS
    )


@pytest.mark.parametrize(
    "kwargs",
    [{"name": "Italian"}, {"name": "it/x"}, {"file_prefix": "italian"}],
)
def test_language_info_validates(kwargs: dict[str, str]) -> None:
    fields = {"code": "zz", "name": "zz", "module": "m", **kwargs}
    with pytest.raises(ValueError):
        LanguageInfo(**fields)


def test_unknown_code() -> None:
    with pytest.raises(PluginNotFoundError, match=r"'xx'.*known: .*cs.*it"):
        spec.load_plugin("xx")


def test_missing_plugin_package_has_a_clear_error() -> None:
    registry = {
        "kl": LanguageInfo("kl", "klingon", "lexpipe_missing_lang.klingon")
    }
    with pytest.raises(PluginNotFoundError) as excinfo:
        spec.load_plugin("kl", registry)
    message = str(excinfo.value)
    assert "klingon plugin" in message and "not installed" in message
    assert "src/lexpipe_missing_lang/klingon/__init__.py" in message


def test_plugins_live_under_lexpipe_langs() -> None:
    registry = {"zz": LanguageInfo("zz", "zzlang", "lexpipe.langs.zzlang")}
    with pytest.raises(PluginNotFoundError) as excinfo:
        spec.load_plugin("zz", registry)
    assert "src/lexpipe/langs/zzlang/__init__.py" in str(excinfo.value)


def test_broken_plugin_is_not_reported_as_missing(plugin_dir: Path) -> None:
    _package(
        plugin_dir,
        "lexpipe_fake_broken",
        "import lexpipe_dependency_that_is_not_installed\n",
    )
    registry = {"zz": LanguageInfo("zz", "broken", "lexpipe_fake_broken")}
    with pytest.raises(ModuleNotFoundError) as excinfo:
        spec.load_plugin("zz", registry)
    assert excinfo.type is ModuleNotFoundError
    assert excinfo.value.name == "lexpipe_dependency_that_is_not_installed"


def test_plugin_module_without_plugin(plugin_dir: Path) -> None:
    _package(plugin_dir, "lexpipe_fake_bare", "VALUE = 1\n")
    registry = {"zz": LanguageInfo("zz", "bare", "lexpipe_fake_bare")}
    with pytest.raises(PluginNotFoundError, match="does not define PLUGIN"):
        spec.load_plugin("zz", registry)


def test_working_plugin(plugin_dir: Path) -> None:
    _package(
        plugin_dir,
        "lexpipe_fake_ok",
        "from lexpipe.spec import PluginSpec\n"
        "def build(ctx):\n"
        "    return 0\n"
        "def analyse(ctx):\n"
        "    return None\n"
        "PLUGIN = PluginSpec(code='zz', name='zzlang',\n"
        "                    stages={'build': build, 'analyse': analyse})\n",
    )
    registry = {"zz": LanguageInfo("zz", "zzlang", "lexpipe_fake_ok")}
    plugin = spec.load_plugin("zz", registry)
    assert isinstance(plugin, spec.LanguagePlugin)
    assert list(plugin.stages) == ["build", "analyse"]
    assert plugin.data_dir == PROJECT_ROOT / "data" / "zzlang"
    assert plugin.analysis_dir == PROJECT_ROOT / "analysis" / "zzlang"

    wrong = {"zz": LanguageInfo("zz", "otherlang", "lexpipe_fake_ok")}
    with pytest.raises(ValueError, match="registry expects"):
        spec.load_plugin("zz", wrong)


def test_plugin_spec_validates_stages() -> None:
    with pytest.raises(ValueError, match="no stages"):
        PluginSpec(code="zz", name="zz", stages={})
    with pytest.raises(ValueError, match="every language"):
        PluginSpec(code="zz", name="zz", stages={"fetch": _noop})
    with pytest.raises(ValueError, match="invalid stage name"):
        PluginSpec(code="zz", name="zz", stages={"Build Lexicon": _noop})
    with pytest.raises(TypeError, match="not callable"):
        PluginSpec(
            code="zz", name="zz", stages={"build": "x"}  # type: ignore[dict-item]
        )


def test_plugin_spec_freezes_stage_order(tmp_path: Path) -> None:
    stages = {"b-stage": _noop, "a-stage": _noop}
    plugin = PluginSpec(code="zz", name="zz", stages=stages, root=tmp_path)
    stages["c-stage"] = _noop
    assert list(plugin.stages) == ["b-stage", "a-stage"]
    with pytest.raises(TypeError):
        plugin.stages["d-stage"] = _noop  # type: ignore[index]
    assert plugin.data_dir == tmp_path / "data" / "zz"


# ---------------------------------------------------------------------------
# Comparison-only rows
# ---------------------------------------------------------------------------


def test_romanian_comparison_row() -> None:
    assert spec.LANGUAGES["ro"] == LanguageInfo(
        code="ro",
        name="romanian_compare",
        module="lexpipe.langs.romanian_compare",
        kaikki_code="ro",
        kaikki_name="Romanian",
        leipzig_id="ron_news_2024_1M",
        file_prefix="romanian_compare_",
        comparison_only=True,
    )
    info = spec.language_info("ro")
    assert info.data_dir() == PROJECT_ROOT / "data" / "romanian_compare"
    assert info.file_name("roots", ".jsonl") == "romanian_compare_roots.jsonl"
    # The flag defaults to False, so one-line test rows stay full ones.
    assert LanguageInfo("zz", "zz", "m").comparison_only is False


def test_full_languages_leave_out_comparison_rows() -> None:
    assert list(spec.full_languages()) == ["it", "cs"]
    registry = {
        "aa": LanguageInfo("aa", "alpha", "m"),
        "bb": LanguageInfo("bb", "beta_compare", "m", comparison_only=True),
    }
    assert list(spec.full_languages(registry)) == ["aa"]


def test_the_romanian_plugin_loads_through_the_registry() -> None:
    plugin = spec.load_plugin("ro")
    assert (plugin.code, plugin.name) == ("ro", "romanian_compare")
    assert "tp" in plugin.stages
    assert spec.plugin_substages(plugin) == {}


# ---------------------------------------------------------------------------
# Sub-stages
# ---------------------------------------------------------------------------


def test_plugin_spec_substages() -> None:
    stages = {"lexicon": _noop, "qc": _noop}
    plugin = PluginSpec(
        code="zz",
        name="zz",
        stages=stages,
        substages={"qc-morfflex": "qc"},
    )
    assert dict(plugin.substages) == {"qc-morfflex": "qc"}
    assert spec.plugin_substages(plugin) == {"qc-morfflex": "qc"}
    with pytest.raises(TypeError):
        plugin.substages["x"] = "qc"  # type: ignore[index]
    # Without substages a plugin has none.
    assert PluginSpec(code="zz", name="zz", stages=stages).substages == {}


@pytest.mark.parametrize(
    ("substages", "message"),
    [
        ({"qc-morfflex": "nonesuch"}, "not a stage"),
        ({"qc": "lexicon"}, "name of a stage"),
        ({"fetch": "qc"}, "every language"),
        ({"Bad Name": "qc"}, "invalid stage name"),
    ],
)
def test_plugin_spec_validates_substages(
    substages: dict[str, str], message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        PluginSpec(
            code="zz",
            name="zz",
            stages={"lexicon": _noop, "qc": _noop},
            substages=substages,
        )


def test_plugin_substages_of_duck_typed_plugins() -> None:
    @dataclass
    class Plain:
        code: str = "zz"
        name: str = "zz"
        data_dir: Path = PROJECT_ROOT
        analysis_dir: Path = PROJECT_ROOT
        stages: dict[str, spec.Stage] = field(
            default_factory=lambda: {"qc": _noop}
        )

    plain = Plain()
    assert isinstance(plain, spec.LanguagePlugin)
    assert spec.plugin_substages(plain) == {}

    @dataclass
    class WithSubstages(Plain):
        substages: object = None

    good = WithSubstages(substages={"qc-db": "qc"})
    assert spec.plugin_substages(good) == {"qc-db": "qc"}
    with pytest.raises(ValueError, match="not a stage"):
        spec.plugin_substages(WithSubstages(substages={"x": "nope"}))
    with pytest.raises(TypeError, match="mapping"):
        spec.plugin_substages(WithSubstages(substages=["qc-db"]))
