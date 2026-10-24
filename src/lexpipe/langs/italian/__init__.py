"""The Italian language plugin (``lexpipe.langs.italian``; design §5).

``PLUGIN`` registers the Italian stages with
``scripts/lang_pipeline.py`` (contract: :mod:`lexpipe.spec`)::

    python scripts/lang_pipeline.py --lang it            # every stage
    python scripts/lang_pipeline.py --lang it lexicon qc

Stages of the data layer, in order (each module's docstring states its
inputs, outputs and rules):

  extract   :mod:`.extract`      one pass over each kaikki source:
                                 projections and the form index (R7)
  lexicon   :mod:`.lexicon`      noun/adjective lexemes, lexeme_id and
                                 root_id (R4), plural provenance (R12),
                                 the Module:it-common drift test
  qc        :mod:`.qc`           plural_qc against Morph-it! and the
                                 it.wiktionary edition; QC report and
                                 manual worklist (tier 2)
  attest    :mod:`.attestation`  Leipzig attestation of every form over
                                 the thresholds τ, form_ambiguous (R11)
  features  :mod:`.features`     cells and lexeme summary (design §5)

Analysis stages (design §4–§5, R15–R27), in order:

  derivations :mod:`.derivations`  derivatives of dorsal bases: template,
                                   text, derived-list and head-form links
                                   plus candidates in both spellings (R21)
  verbs       :mod:`.verbs`        velar-stem verbs, conjugation class,
                                   front-vowel cells (E1, R15)
  tp          :mod:`.tpstage`      harmonized TP long table, model fits,
                                   simulated learners, digest (R16–R19,
                                   R24; rows: :mod:`.tptables`, learners:
                                   :mod:`.tplearners`, digest:
                                   :mod:`.tpsummary`)
  covariation :mod:`.covariation`  base plural × derivative (R23)
  spellout    :mod:`.spellout`     spell-out TP on certain roots (R20)
  noevidence  :mod:`.noevidence`   roots without alternation evidence (R21)
  logtable    :mod:`.logtable`     -logo/-fago/-urgo table (R27)
  docs        :mod:`.docstage`     schema check, data dictionary, README
                                   results (:mod:`.schema`)

Supporting modules: :mod:`.layout` (paths), :mod:`.store` (stamped
derivatives), :mod:`.keys` (match keys), :mod:`.heads` (head-template
grammar), :mod:`.itcommon` (port of Module:it-common), :mod:`.plurals`
(per-form plural table), :mod:`.etymology`, :mod:`.morphology`,
:mod:`.morphit`, :mod:`.leipzig`, :mod:`.suffixes` (derivational
suffixes and stem-edge readings), :mod:`.roots` (root records and their
observations for the core TP code).
"""

from __future__ import annotations

from typing import Final

from lexpipe.langs.italian.stages import stage_table
from lexpipe.spec import LANGUAGES, PluginSpec

_LANGUAGE_CODE: Final[str] = "it"

PLUGIN: Final[PluginSpec] = PluginSpec(
    code=LANGUAGES[_LANGUAGE_CODE].code,
    name=LANGUAGES[_LANGUAGE_CODE].name,
    stages=stage_table(),
)
