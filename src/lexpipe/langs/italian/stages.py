"""The Italian stages in pipeline order, and how to add one.

:data:`STAGE_TABLE` lists ``(stage name, "module:function")`` pairs in
the order ``scripts/lang_pipeline.py`` runs them (after the common
``fetch`` and ``verify-sources``). Each target is a function of a
:class:`lexpipe.context.RunContext` returning None or 0 on success and
a non-zero int to stop the pipeline (:data:`lexpipe.spec.Stage`).

Targets are imported only when their stage runs (:func:`stage_table`),
so importing the plugin stays cheap and a stage module can import
another stage's module without cycles. ``tests/italian`` resolves every
target, so a typo fails a test rather than a run.

Adding a stage (the analysis stages — derivations, verbs, tp,
covariation, spellout, noevidence, logtable, docs — are added this way):

  1. write ``src/lexpipe/langs/italian/<module>.py`` with
     ``run(ctx: RunContext) -> int | None`` (use
     :func:`lexpipe.context.run_stage` so the stage is skipped when its
     manifest matches, and :class:`lexpipe.langs.italian.layout.Layout`
     for paths, adding the new files there);
  2. append ``("<name>", "lexpipe.langs.italian.<module>:run")`` to
     :data:`STAGE_TABLE` at its place in the pipeline order.

Stages built by the data layer: ``extract``, ``lexicon``, ``qc``,
``attest``, ``features``. Analysis stages: ``derivations``, ``verbs``,
``tp``, ``covariation``, ``spellout``, ``noevidence``, ``logtable``,
``docs`` (see the package docstring).
"""

from __future__ import annotations

import importlib
from collections.abc import Callable
from typing import TYPE_CHECKING, Final

from lexpipe.spec import Stage

if TYPE_CHECKING:
    from lexpipe.context import RunContext

STAGE_TABLE: Final[tuple[tuple[str, str], ...]] = (
    ("extract", "lexpipe.langs.italian.extract:run"),
    ("lexicon", "lexpipe.langs.italian.lexicon:run"),
    ("qc", "lexpipe.langs.italian.qc:run"),
    ("attest", "lexpipe.langs.italian.attestation:run"),
    ("features", "lexpipe.langs.italian.features:run"),
    ("derivations", "lexpipe.langs.italian.derivations:run"),
    ("verbs", "lexpipe.langs.italian.verbs:run"),
    ("tp", "lexpipe.langs.italian.tpstage:run"),
    ("covariation", "lexpipe.langs.italian.covariation:run"),
    ("spellout", "lexpipe.langs.italian.spellout:run"),
    ("noevidence", "lexpipe.langs.italian.noevidence:run"),
    ("logtable", "lexpipe.langs.italian.logtable:run"),
    ("docs", "lexpipe.langs.italian.docstage:run"),
)


def resolve(target: str) -> Stage:
    """Import ``"package.module:function"`` and return the function."""
    module_name, sep, attribute = target.partition(":")
    if not sep or not module_name or not attribute:
        raise ValueError(f"stage target must be 'module:function': {target}")
    function: object = getattr(importlib.import_module(module_name), attribute)
    if not callable(function):
        raise TypeError(f"stage target {target} is not callable")
    return _as_stage(function)


def _as_stage(function: Callable[..., object]) -> Stage:
    def stage(ctx: RunContext) -> int | None:
        result = function(ctx)
        if result is None or isinstance(result, int):
            return result
        raise TypeError(
            f"stage returned {type(result).__name__}, expected int or None"
        )

    return stage


def _lazy(target: str) -> Stage:
    def stage(ctx: RunContext) -> int | None:
        return resolve(target)(ctx)

    stage.__name__ = target.rpartition(".")[2].replace(":", "_")
    stage.__qualname__ = stage.__name__
    stage.__doc__ = f"Run {target} (imported when the stage runs)."
    return stage


def stage_table() -> dict[str, Stage]:
    """``{stage name: stage}`` in pipeline order."""
    return {name: _lazy(target) for name, target in STAGE_TABLE}
