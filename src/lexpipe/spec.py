"""Language registry, the language-plugin contract and the plugin loader.

Language identity (design rule R3)
----------------------------------
Each language has exactly one registry row, a :class:`LanguageInfo` in
:data:`LANGUAGES`, keyed by its ISO 639-1 code. Paths, file names and
source identifiers are derived from that row, never spelled out
elsewhere:

    field            Italian                  Czech                 Romanian
    code             it                       cs                    ro
    name             italian                  czech                 romanian_compare
    module           lexpipe.langs.italian    lexpipe.langs.czech   lexpipe.langs.
                                                                    romanian_compare
    kaikki_code      it                       cs                    ro
    kaikki_name      Italian                  Czech                 Romanian
    leipzig_id       ita_news_2024_1M         ces_news_2024_1M      ron_news_2024_1M
    file_prefix      italian_                 czech_                romanian_compare_
    analysis_module  phonology.analyses.      phonology.analyses.   (none)
                     italian_palatalization   czech_palatalization
    elision_clitics  ITALIAN_ELISION_CLITICS  (none)                (none)
    comparison_only  False                    False                 True

``kaikki_code`` is the top-level ``lang_code`` of wiktextract entries
and ``kaikki_name`` their ``lang`` (the prefix of their category
names); ``elision_clitics`` drives the token rules of
:mod:`lexpipe.attest`. Directories are resolved against a root by
:class:`lexpipe.context.RunContext`.

Comparison-only rows
--------------------
A row with ``comparison_only=True`` names a plugin that builds no
lexicon of its own: it reads another pipeline's outputs and writes
cross-language tables. Romanian (``ro``, design §7) is one: its
pipeline is the legacy one (``run_pipeline.sh``, R1), and
:mod:`lexpipe.langs.romanian_compare` only recounts its outputs with
the shared core. Such a language has no normalizer, feature
inventory, ``docs/<name>/``, ``curated/<name>/``,
``analysis/<name>/`` or raw lexicon, and its run context may write
only under ``analysis/comparison/`` and ``data/<name>/``
(:meth:`lexpipe.context.RunContext.check_writable`).
:func:`full_languages` lists the other rows.

Plugin contract
---------------
A language plugin is the package named by ``module`` (under
:mod:`lexpipe.langs`); it exposes ``PLUGIN``, an object satisfying
:class:`LanguagePlugin`:

  - ``code`` / ``name`` — equal to the registry row's;
  - ``data_dir`` / ``analysis_dir`` — ``data/<name>/`` and
    ``analysis/<name>/`` under the plugin's root;
  - ``stages`` — an ordered mapping from stage name to a
    :data:`Stage`: a function of a :class:`lexpipe.context.RunContext`
    returning None or 0 on success and a non-zero int to stop the
    pipeline, mirroring the ``main() -> int`` convention of the
    scripts.

:class:`PluginSpec` is a ready-made implementation; a plugin's
``__init__.py`` can simply say::

    PLUGIN = PluginSpec(code="it", name="italian", stages={...})

``scripts/lang_pipeline.py`` adds the stages in :data:`COMMON_STAGES`
(``fetch``, ``verify-sources``) in front of every plugin's own, so a
plugin must not reuse those names.

Sub-stages. A stage may run nested steps that keep manifests of their
own (Czech ``qc`` runs ``qc-morfflex`` through
:func:`lexpipe.context.run_stage`). A plugin may declare them,
optionally, as ``substages``: a mapping from sub-stage name to the
stage that runs it (``PluginSpec(..., substages={"qc-morfflex":
"qc"})``). ``--force`` then accepts the sub-stage's name: it is put
into ``RunContext.force``, so the nested ``run_stage`` call re-runs.
Sub-stages are never run on their own. :func:`plugin_substages`
reads the mapping of any plugin (empty when it declares none), so
the attribute is not part of :class:`LanguagePlugin`.

Plugins are imported only when asked for (:func:`load_plugin`), so a
missing or broken plugin never affects another language, and a
missing one fails with :class:`PluginNotFoundError` naming the package
that should exist. A plugin that exists but fails to import (say, a
missing third-party dependency) raises its own error unchanged.
"""

from __future__ import annotations

import importlib
import re
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING, Final, Protocol, TypeAlias, runtime_checkable

from lexpipe import PROJECT_ROOT

if TYPE_CHECKING:
    from lexpipe.context import RunContext

Stage: TypeAlias = Callable[["RunContext"], int | None]

# Stages scripts/lang_pipeline.py provides for every language.
COMMON_STAGES: Final[tuple[str, ...]] = ("fetch", "verify-sources")

# Stage names appear on the command line: lowercase, digits, hyphens.
_STAGE_NAME_RE: Final[re.Pattern[str]] = re.compile(r"[a-z][a-z0-9-]*")
# Codes and directory names end up in paths and module names.
_CODE_RE: Final[re.Pattern[str]] = re.compile(r"[a-z]{2,3}")
_NAME_RE: Final[re.Pattern[str]] = re.compile(r"[a-z][a-z_]*")

# Italian words elided before a vowel with an apostrophe (l'amico,
# dell'anno, un'altra, quest'anno, c'è, d'Italia; ' or ’). The word
# after the apostrophe is counted as a token of its own (R11).
#   articles and articulated prepositions: l, un, dell, dall, nell,
#     sull, all, coll, pell;
#   determiners and adjectives: quell, quest, bell, sant, tutt, buon;
#   clitic pronouns and particles, and elided di: c, d, m, n, s, t, v,
#     gl.
# "quest" is not in the design's first list; it is added because
# ita_news_2024_1M has 7,023 tokens starting with quest' or quest’
# (quest’anno alone: 2,506).
ITALIAN_ELISION_CLITICS: Final[frozenset[str]] = frozenset(
    {
        "l",
        "un",
        "dell",
        "dall",
        "nell",
        "sull",
        "all",
        "coll",
        "pell",
        "quell",
        "quest",
        "bell",
        "sant",
        "tutt",
        "buon",
        "c",
        "d",
        "m",
        "n",
        "s",
        "t",
        "v",
        "gl",
    }
)


@runtime_checkable
class LanguagePlugin(Protocol):
    """What ``scripts/lang_pipeline.py`` needs from a language plugin."""

    @property
    def code(self) -> str: ...

    @property
    def name(self) -> str: ...

    @property
    def data_dir(self) -> Path: ...

    @property
    def analysis_dir(self) -> Path: ...

    @property
    def stages(self) -> Mapping[str, Stage]: ...


@dataclass(frozen=True, slots=True)
class LanguageInfo:
    """One registry row (see the module docstring for the fields).

    Only ``code``, ``name`` and ``module`` are required, so test
    registries can describe a throwaway plugin in one line.
    ``comparison_only`` marks a plugin that only builds cross-language
    tables from another pipeline's outputs (module docstring).
    """

    code: str
    name: str
    module: str
    kaikki_code: str = ""
    kaikki_name: str = ""
    leipzig_id: str = ""
    file_prefix: str = ""
    analysis_module: str = ""
    elision_clitics: frozenset[str] = frozenset()
    comparison_only: bool = False

    def __post_init__(self) -> None:
        if not _NAME_RE.fullmatch(self.name):
            raise ValueError(f"invalid language name: {self.name!r}")
        if self.file_prefix and not self.file_prefix.endswith("_"):
            raise ValueError(
                f"{self.code}: file_prefix must end with '_': "
                f"{self.file_prefix!r}"
            )

    def data_dir(self, root: Path = PROJECT_ROOT) -> Path:
        return root / "data" / self.name

    def analysis_dir(self, root: Path = PROJECT_ROOT) -> Path:
        return root / "analysis" / self.name

    def file_name(self, stem: str, suffix: str = ".csv") -> str:
        """Prefixed file name: ``file_name("lexicon_raw")`` →
        ``italian_lexicon_raw.csv``."""
        return f"{self.file_prefix}{stem}{suffix}"


def _row(
    code: str,
    name: str,
    kaikki_name: str,
    leipzig_id: str,
    elision_clitics: frozenset[str] = frozenset(),
) -> LanguageInfo:
    if not _CODE_RE.fullmatch(code):
        raise ValueError(f"invalid language code: {code!r}")
    return LanguageInfo(
        code=code,
        name=name,
        module=f"lexpipe.langs.{name}",
        kaikki_code=code,
        kaikki_name=kaikki_name,
        leipzig_id=leipzig_id,
        file_prefix=f"{name}_",
        analysis_module=f"phonology.analyses.{name}_palatalization",
        elision_clitics=elision_clitics,
    )


LANGUAGES: Final[Mapping[str, LanguageInfo]] = MappingProxyType(
    {
        "it": _row(
            "it",
            "italian",
            "Italian",
            "ita_news_2024_1M",
            ITALIAN_ELISION_CLITICS,
        ),
        "cs": _row("cs", "czech", "Czech", "ces_news_2024_1M"),
        # Comparison only (design §7): the harmonized recount of the
        # legacy Romanian pipeline's outputs. Name and prefix keep its
        # files apart from the Romanian pipeline's own (R1).
        "ro": LanguageInfo(
            code="ro",
            name="romanian_compare",
            module="lexpipe.langs.romanian_compare",
            kaikki_code="ro",
            kaikki_name="Romanian",
            leipzig_id="ron_news_2024_1M",
            file_prefix="romanian_compare_",
            comparison_only=True,
        ),
    }
)


def full_languages(
    registry: Mapping[str, LanguageInfo] = LANGUAGES,
) -> dict[str, LanguageInfo]:
    """The rows that are not comparison-only, in registry order.

    These languages have the whole per-language layout: a normalizer,
    a feature inventory, ``docs/<name>/``, ``curated/<name>/``,
    ``analysis/<name>/`` and a tracked raw lexicon.
    """
    return {
        code: info
        for code, info in registry.items()
        if not info.comparison_only
    }


class PluginNotFoundError(LookupError):
    """No plugin is registered or importable for a language code."""


def _check_stage_name(stage_name: str) -> None:
    if not _STAGE_NAME_RE.fullmatch(stage_name):
        raise ValueError(f"invalid stage name: {stage_name!r}")
    if stage_name in COMMON_STAGES:
        raise ValueError(
            f"stage {stage_name!r} is provided for every language "
            "by scripts/lang_pipeline.py; choose another name"
        )


def _check_substages(
    substages: Mapping[str, str], stages: Mapping[str, Stage]
) -> None:
    """Sub-stage names are stage-like, new, and attached to a stage."""
    for sub, parent in substages.items():
        _check_stage_name(sub)
        if sub in stages:
            raise ValueError(
                f"sub-stage {sub!r} has the name of a stage; choose another"
            )
        if parent not in stages:
            raise ValueError(
                f"sub-stage {sub!r} belongs to {parent!r}, which is not a "
                "stage of the plugin"
            )


@dataclass(frozen=True)
class PluginSpec:
    """A plain :class:`LanguagePlugin` built from code, name and stages.

    Stage names are validated (command-line safe, not one of
    :data:`COMMON_STAGES`) and the mapping is frozen in its order.
    ``substages`` maps the name of a nested step with its own manifest
    to the stage that runs it (module docstring); it is frozen too.
    """

    code: str
    name: str
    stages: Mapping[str, Stage]
    root: Path = field(default=PROJECT_ROOT, compare=False)
    substages: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.stages:
            raise ValueError(f"plugin {self.code!r} defines no stages")
        for stage_name, stage in self.stages.items():
            _check_stage_name(stage_name)
            if not callable(stage):
                raise TypeError(f"stage {stage_name!r} is not callable")
        _check_substages(self.substages, self.stages)
        object.__setattr__(self, "stages", MappingProxyType(dict(self.stages)))
        object.__setattr__(
            self, "substages", MappingProxyType(dict(self.substages))
        )

    @property
    def data_dir(self) -> Path:
        return self.root / "data" / self.name

    @property
    def analysis_dir(self) -> Path:
        return self.root / "analysis" / self.name


def plugin_substages(plugin: LanguagePlugin) -> Mapping[str, str]:
    """``{sub-stage: parent stage}`` of *plugin*; empty if it has none.

    Plugins that are not a :class:`PluginSpec` may define a
    ``substages`` mapping too; it is validated like PluginSpec's.
    """
    substages = getattr(plugin, "substages", None)
    if substages is None:
        return MappingProxyType({})
    if not isinstance(substages, Mapping):
        raise TypeError(
            f"plugin {plugin.code!r}: substages must be a mapping, "
            f"got {type(substages).__name__}"
        )
    checked = {str(sub): str(parent) for sub, parent in substages.items()}
    _check_substages(checked, plugin.stages)
    return MappingProxyType(checked)


def language_info(
    code: str, registry: Mapping[str, LanguageInfo] = LANGUAGES
) -> LanguageInfo:
    """Registry row for *code*; :class:`PluginNotFoundError` if unknown."""
    try:
        return registry[code]
    except KeyError:
        known = ", ".join(sorted(registry)) or "none"
        raise PluginNotFoundError(
            f"no language plugin registered for {code!r} (known: {known})"
        ) from None


def _module_is_missing(exc: ModuleNotFoundError, module: str) -> bool:
    """True if *exc* is about *module* or one of its parent packages."""
    missing = exc.name or ""
    return bool(missing) and (
        module == missing or module.startswith(missing + ".")
    )


def load_plugin(
    code: str, registry: Mapping[str, LanguageInfo] = LANGUAGES
) -> LanguagePlugin:
    """Import the plugin for *code* and return its checked ``PLUGIN``."""
    info = language_info(code, registry)
    try:
        module = importlib.import_module(info.module)
    except ModuleNotFoundError as exc:
        if not _module_is_missing(exc, info.module):
            raise
        path = "src/" + info.module.replace(".", "/") + "/__init__.py"
        raise PluginNotFoundError(
            f"the {info.name} plugin ({info.module}) is not installed: "
            f"expected {path} defining PLUGIN"
        ) from exc
    plugin = getattr(module, "PLUGIN", None)
    if plugin is None:
        raise PluginNotFoundError(f"{info.module} does not define PLUGIN")
    if not isinstance(plugin, LanguagePlugin):
        raise TypeError(
            f"{info.module}.PLUGIN does not satisfy LanguagePlugin "
            "(code, name, data_dir, analysis_dir, stages)"
        )
    if (plugin.code, plugin.name) != (info.code, info.name):
        raise ValueError(
            f"{info.module}.PLUGIN is ({plugin.code!r}, {plugin.name!r}); "
            f"the registry expects ({info.code!r}, {info.name!r})"
        )
    clash = [name for name in plugin.stages if name in COMMON_STAGES]
    if clash:
        raise ValueError(
            f"{info.module}.PLUGIN redefines common stages: {clash}"
        )
    return plugin
