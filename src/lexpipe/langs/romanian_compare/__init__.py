"""Read-only harmonized Romanian comparison (design §7; R16–R23, R27).

The Romanian pipeline keeps producing its own outputs byte for byte.
This plugin only *reads* them (``data/romanian_lexicon_with_freq.csv``,
``data/verb_paradigm_verdict.csv``, ``data/leipzig/freq/`` and the
legacy TP summaries) and recounts Romanian with the shared core
(:mod:`lexpipe.tp`, :mod:`lexpipe.attest`, :mod:`lexpipe.learners`,
:mod:`lexpipe.stats`, :mod:`lexpipe.suffixoids`), so Romanian, Italian
and Czech come from one counting code. It writes only under
``analysis/comparison/`` (tracked aggregates) and
``data/romanian_compare/`` (untracked derivatives); the run context
enforces it (:func:`lexpipe.context.write_confinements`).

Modules:

  - :mod:`.lexicon` — lexicon rows, outcome-blind root attributes;
  - :mod:`.morphology` — orthographic reading of verb cells and
    derivatives (E1, E2);
  - :mod:`.corpus` — the Romanian Leipzig table as a core
    :class:`lexpipe.attest.AttestationTable`;
  - :mod:`.records` — the per-root derivative (roots JSONL);
  - :mod:`.observations` — records → :class:`lexpipe.tp.RootObservation`;
  - :mod:`.harmonized` — TP long rows, V sweep, nesting, learners;
  - :mod:`.spellout`, :mod:`.noevidence`, :mod:`.covariation`,
    :mod:`.logtable` — R20, R21, R23, R27;
  - :mod:`.parity` — the legacy tables reproduced and reconciled;
  - :mod:`.readme`, :mod:`.stages` — outputs and stage functions.

Registry. :data:`RO_INFO` is the comparison-only registry row
``lexpipe.spec.LANGUAGES["ro"]`` (code ``ro``, name
``romanian_compare``, so ``data_dir`` is ``data/romanian_compare/``).
Run the stages with ``python scripts/lang_pipeline.py --lang ro
[--root DIR] [--seed N] [--force STAGE] [--param KEY=VALUE] [stage
...]``; ``python -m lexpipe.langs.romanian_compare`` is kept as a
shortcut to the same command.
"""

from __future__ import annotations

from typing import Final

from lexpipe.langs.romanian_compare.stages import STAGES
from lexpipe.spec import LANGUAGES, LanguageInfo, PluginSpec

RO_INFO: Final[LanguageInfo] = LANGUAGES["ro"]

PLUGIN: Final[PluginSpec] = PluginSpec(
    code=RO_INFO.code, name=RO_INFO.name, stages=STAGES
)

__all__ = ["PLUGIN", "RO_INFO"]
