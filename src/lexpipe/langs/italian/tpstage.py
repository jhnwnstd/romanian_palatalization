"""Stage ``tp``: the Italian TP long table, learners and digest (R16–R19).

Inputs: ``italian_lexicon_complete.csv`` and ``italian_cells.csv``
(``features``), ``italian_verbs.csv`` (``verbs``, E1) and
``italian_derivations.csv`` (``derivations``, E2), loaded into root
records by :mod:`.roots`.

Outputs (tracked, aggregates only, R9):

  ``analysis/italian/italian_tp_long.csv``
      :data:`lexpipe.tp.TP_LONG_TRACKED_FIELDS` (the columns of every
      language's long table, e.g. the Romanian harmonized table
      ``analysis/comparison/romanian_harmonized_tp_long.csv``); tables
      R17a/R17b/R17c, cells, nde, it_stress, it_suffix, it_verbs,
      it_nounverb, logtable (:mod:`.tptables`) and spellout
      (:mod:`.spellout`), in that order.
  ``analysis/italian/italian_learners.csv``
      simulated learner shares (:mod:`.tplearners`,
      :data:`.tplearners.LEARNER_FIELDS`).
  ``analysis/italian/italian_tp_summary.txt``
      the readable digest (:mod:`.tpsummary`).

Both CSVs go through :func:`lexpipe.csvio.write_tracked` and the
digest through :meth:`lexpipe.context.RunContext.write_tracked_text`
(the licence gate); the τ rows use Leipzig counts only as attestation
thresholds and carry aggregates, never per-word counts. The seed comes
from the run context, the number of simulated vocabularies from its
``learner_samples`` parameter (default 1000), the vocabulary sizes from
``learner_V`` (default {500, 1000, 2000, 5000}; tests pass fewer and
smaller) and the worker processes of the simulation from ``workers``
(default min(8, CPUs); the result does not depend on it).
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from pathlib import Path
from typing import Final

from lexpipe import csvio, learners
from lexpipe.context import RunContext, run_stage
from lexpipe.langs.italian import (
    roots,
    spellout,
    tplearners,
    tpsummary,
    tptables,
)
from lexpipe.langs.italian.layout import Layout

STAGE: Final[str] = "tp"
# 2: core simulate_all/status_verdict, tracked_row, shared R27 classes.
# 3: cells keep doublets apart (amici/amiche: variable), doublet minimum,
#    R26 suffix units counted once, dictionary_explicit and E2_all rows,
#    the provenance section of the digest.
# 4: the spell-out subset all+candidates+other_base.
# 5: copies of a plural form (noun and adjective) merge their τ flags, so
#    a typed plural counts at dictionary_explicit (roots.realizations_for).
VERSION: Final[int] = 5


def _body(
    ctx: RunContext,
    layout: Layout,
    n_samples: int,
    v_values: Sequence[int],
    workers: int,
) -> int:
    records = spellout.load_records(layout)
    items = spellout.load_derivatives(layout.derivations)
    rows = tptables.all_rows(records, spellout.spellout_rows(records, items))
    learner_rows = tplearners.learner_rows(
        records,
        seed=ctx.seed,
        n_samples=n_samples,
        V_values=v_values,
        workers=workers,
    )
    csvio.write_tracked(layout.tp_long, rows, tptables.FIELDS, ctx=ctx)
    csvio.write_tracked(
        layout.learners, learner_rows, tplearners.LEARNER_FIELDS, ctx=ctx
    )
    ctx.write_tracked_text(
        layout.tp_summary,
        tpsummary.summary(
            rows,
            learner_rows,
            records.counts,
            n_samples,
            ctx.seed,
            provenance=roots.dictionary_provenance(layout.cells),
        ),
        sources=tptables.SOURCE,
    )
    print(
        f"  {len(rows):,} TP rows, {len(learner_rows):,} learner rows "
        f"({n_samples} vocabularies per V)"
    )
    return 0


def outputs(layout: Layout) -> list[Path]:
    return [layout.tp_long, layout.learners, layout.tp_summary]


def run(ctx: RunContext) -> int | None:
    """Run the ``tp`` stage (skipped when its manifest matches)."""
    layout = Layout(ctx)
    for path in outputs(layout):
        ctx.check_writable(path)
    workers = int(
        ctx.params.get(tplearners.WORKERS_PARAM, min(8, os.cpu_count() or 1))
    )
    n_samples = int(
        ctx.params.get(tplearners.SAMPLES_PARAM, learners.DEFAULT_SAMPLES)
    )
    v_values = [
        int(v)
        for v in ctx.params.get(tplearners.V_PARAM, learners.DEFAULT_V_VALUES)
    ]
    return run_stage(
        ctx,
        STAGE,
        inputs=spellout.inputs(layout),
        params={
            "version": VERSION,
            "seed": ctx.seed,
            "learner_samples": n_samples,
            "V_values": v_values,
        },
        outputs=outputs(layout),
        body=lambda: _body(ctx, layout, n_samples, v_values, workers),
    )
