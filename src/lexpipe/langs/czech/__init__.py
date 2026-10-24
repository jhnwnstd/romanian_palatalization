"""Czech language plugin (design §2, §3, §6).

The data layer of the Czech pipeline, run with
``python scripts/lang_pipeline.py --lang cs [stage ...]``:

    extract    compact projections of the English and Czech Wiktionary
               editions and a form index (:mod:`.extract`)
    lexicon    noun and adjective lexemes, R13 cells, domain (R25)
               (:mod:`.lexicon`)
    qc         MorfFlex agreement, report and manual worklist (:mod:`.qc`);
               its sub-stage ``qc-morfflex`` builds the MorfFlex projection
               (forced with ``--force qc-morfflex``)
    attest     Leipzig attestation and the CNK cross-check
               (:mod:`.attestation`)
    features   stems, targets, P/N endings, mutations, exception labels
               and lexeme summaries (:mod:`.features`)

and the analyses (design §4, §6; R14–R27):

    derivations  base × derivative pairs from Wiktionary, DeriNet 2.3
                 and generated candidates; R22 strata (:mod:`.derivations`)
    tp           (a) rule TP, (b) competition TP, (c) leave-one-cell-out
                 spell-out TP, the R17 rule family, simulated learners
                 (:mod:`.tpstage`)
    spellout     spell-out TP on certain-specification roots (R20,
                 :mod:`.spellout`)
    noevidence   roots without corpus alternation evidence × their
                 derivatives (R21, :mod:`.noevidence`)
    covariation  inflectional class × derivative suffix choice (R23,
                 :mod:`.covariation`)
    logtable     the -log/-urg/-graf/-fil table (R27, :mod:`.logtable`)
    docs         data dictionary generated from :mod:`.schema`

Paths are listed in :mod:`.layout` and :mod:`.analysis_layout`.
"""

from __future__ import annotations

from typing import Final

from lexpipe.langs.czech import (
    attestation,
    covariation,
    derivations,
    extract,
    features,
    lexicon,
    logtable,
    noevidence,
    qc,
    schema,
    spellout,
    tpstage,
)
from lexpipe.spec import PluginSpec, Stage

STAGES: Final[dict[str, Stage]] = {
    "extract": extract.run,
    "lexicon": lexicon.run,
    "qc": qc.run,
    "attest": attestation.run,
    "features": features.run,
    "derivations": derivations.run,
    "tp": tpstage.run,
    "spellout": spellout.run,
    "noevidence": noevidence.run,
    "covariation": covariation.run,
    "logtable": logtable.run,
    "docs": schema.run,
}

PLUGIN: Final[PluginSpec] = PluginSpec(
    code="cs",
    name="czech",
    stages=STAGES,
    # qc runs the MorfFlex projection with a manifest of its own
    # (python scripts/lang_pipeline.py --lang cs --force qc-morfflex qc).
    substages={qc.DB_STAGE: "qc"},
)
