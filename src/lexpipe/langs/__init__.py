"""Language plugins of the multilingual pipeline (design rule R3).

One package per language, named after the registry row's ``name``:
``lexpipe.langs.italian`` and ``lexpipe.langs.czech``. Each package's
``__init__.py`` defines ``PLUGIN``, an object satisfying
:class:`lexpipe.spec.LanguagePlugin`; ``scripts/lang_pipeline.py``
imports it only when one of its stages is needed
(:func:`lexpipe.spec.load_plugin`).

Plugins live here rather than in a top-level ``languages`` package so
that they share the project's namespace and cannot shadow, or be
shadowed by, an unrelated installed module.

Contracts for plugin code:

  - language identity, paths and source identifiers come from the
    plugin's :class:`lexpipe.spec.LanguageInfo` row, never from
    literals;
  - stages are functions of a :class:`lexpipe.context.RunContext`
    and keep no module-level state;
  - nothing here imports the Romanian pipeline
    (``wiktionary_normalizer``, ``romanian_processor_lib``,
    ``dex_utils``, ``romanian_harvester``); ``tests/test_repo_guards.py``
    enforces it.
"""
