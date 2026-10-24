"""Stages of the read-only harmonized Romanian comparison (design §7).

Every stage is a function of a :class:`lexpipe.context.RunContext`
whose language is the comparison-only registry row ``ro``
(:data:`lexpipe.langs.romanian_compare.RO_INFO`, which is
``lexpipe.spec.LANGUAGES["ro"]``). The context confines it to
``analysis/comparison/`` and ``data/romanian_compare/``
(:func:`lexpipe.context.write_confinements`); the Romanian inputs are
only read (R1, R2). Each stage runs through
:func:`lexpipe.context.run_stage`, so it is skipped when its manifest
(``data/romanian_compare/manifests/<stage>.manifest.tsv``) still
matches, and its parameters carry :data:`VERSION`.

    roots        lexicon + verb verdicts + Leipzig → roots JSONL (R7)
    tp           harmonized TP long rows (cells, R17, NDE, -log, spell-out)
    learners     simulated learners for the R17 cells (R18; seeded)
    spellout     spell-out TP report (R20)
    noevidence   no-evidence report (R21)
    covariation  plural class × verbalizer choice (R23)
    logtable     -log/-fag/-urg table (R27)
    parity       legacy reproduction and side-by-side report
    report       the Romanian section of analysis/comparison/README.md

Run them with ``python scripts/lang_pipeline.py --lang ro [stage ...]``.

Tracked outputs (``analysis/comparison/``) carry aggregates or lexicon
facts only: CSVs go through :func:`lexpipe.csvio.write_tracked` and
text reports through :meth:`lexpipe.context.RunContext.write_tracked_text`
(R9), both naming the lexicon's editions (:data:`.harmonized.SOURCE`);
everything with per-word Leipzig counts stays in the gitignored
``data/romanian_compare/``. The -log table reads the curated decisions
``curated/romanian_compare/suffixoid_classes.csv`` (guards checked).
Run parameters in ``ctx.params``: ``learner_samples`` (default 1000),
``learner_V`` (default 500, 1000, 2000, 5000; a fixture needs smaller
values), ``workers`` (processes of the learner simulation, default
min(8, CPUs); the result does not depend on it) and
``n_permutations`` (default 10,000); the seed is ``ctx.seed``.
"""

from __future__ import annotations

import os
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

from lexpipe import csvio, learners, stats, suffixoids
from lexpipe.attest import TAUS
from lexpipe.context import RunContext, run_stage, sha256_file
from lexpipe.langs.romanian_compare import (
    covariation,
)
from lexpipe.langs.romanian_compare import harmonized as hz
from lexpipe.langs.romanian_compare import lexicon as lx
from lexpipe.langs.romanian_compare import (
    logtable,
    noevidence,
    parity,
    readme,
    spellout,
)
from lexpipe.langs.romanian_compare.corpus import load_counts
from lexpipe.langs.romanian_compare.observations import Attester
from lexpipe.langs.romanian_compare.records import (
    RootsFile,
    build_records,
    read_roots,
    read_verdicts,
    roots_lines,
)
from lexpipe.spec import Stage

VERSION: Final[str] = "romanian_compare/1"
# Stages whose code changed with the shared-core migration (shared
# suffixoid classifiers and curated decisions, simulate_all, source
# columns, the licence-gated text writer) carry this version instead.
VERSION_2: Final[str] = "romanian_compare/2"
# The -log table without the nouns whose plural DEX filled (and the
# README section that reads it) carry this version.
VERSION_3: Final[str] = "romanian_compare/3"
LOG: Final[Callable[[str], None]] = print


@dataclass(frozen=True, slots=True)
class Paths:
    """Inputs (read-only) and outputs of the stages."""

    lexicon: Path
    verdicts: Path
    leipzig: Path
    legacy_full: Path
    legacy_downsampled: Path
    check_claims: Path
    roots: Path
    tp_long: Path
    learners: Path
    spellout: Path
    noevidence: Path
    noevidence_units: Path
    covariation: Path
    covariation_units: Path
    logtable: Path
    logtable_full: Path
    decisions: Path
    parity: Path
    readme: Path


def paths(ctx: RunContext) -> Paths:
    """Every path, derived from the root and the registry row (R3)."""
    data = ctx.root / "data"
    analysis = ctx.root / "analysis"
    info = ctx.language
    out = ctx.comparison_dir
    return Paths(
        lexicon=data / "romanian_lexicon_with_freq.csv",
        verdicts=data / "verb_paradigm_verdict.csv",
        leipzig=data / "leipzig" / "freq" / f"{info.leipzig_id}_freq.csv",
        legacy_full=analysis / "romanian_tp_summary_full.csv",
        legacy_downsampled=analysis / "romanian_tp_summary_downsampled.csv",
        check_claims=analysis
        / "learnability_checks"
        / "check_claims_output.txt",
        roots=ctx.data_dir / info.file_name("roots", ".jsonl"),
        tp_long=out / "romanian_harmonized_tp_long.csv",
        learners=out / "romanian_learners.csv",
        spellout=out / "romanian_spellout.txt",
        noevidence=out / "romanian_noevidence.txt",
        noevidence_units=ctx.data_dir / info.file_name("noevidence_units"),
        covariation=out / "romanian_covariation.txt",
        covariation_units=ctx.data_dir / info.file_name("covariation_units"),
        logtable=out / "romanian_logtable.csv",
        logtable_full=ctx.data_dir / info.file_name("logtable_full"),
        decisions=logtable.decisions_path(ctx.curated_dir),
        parity=out / "romanian_parity.txt",
        readme=out / "README.md",
    )


def write_tracked(
    ctx: RunContext,
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    fields: Sequence[str],
) -> int:
    """Guarded (context) and licence-gated (R9) tracked CSV writer."""
    return csvio.write_tracked(path, rows, fields, ctx=ctx)


def write_untracked(
    ctx: RunContext,
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    fields: Sequence[str],
) -> int:
    """Guarded CSV writer for ``data/romanian_compare/``."""
    return csvio.write_rows(path, rows, fields, ctx=ctx)


def write_report(ctx: RunContext, path: Path, text: str) -> Path:
    """A tracked text report behind the licence gate (R9): aggregates
    of the lexicon's two Wiktionary editions."""
    return ctx.write_tracked_text(path, text, sources=hz.SOURCE)


def _load(ctx: RunContext) -> tuple[RootsFile, Attester, hz.ObservationCache]:
    roots = read_roots(paths(ctx).roots)
    attest = Attester(roots)
    return roots, attest, hz.ObservationCache(roots.records, attest)


def load_decisions(
    ctx: RunContext, roots: RootsFile
) -> dict[str, suffixoids.Decision]:
    """The curated -log decisions, guards checked against the roots."""
    rows, decisions = logtable.load_decisions(paths(ctx).decisions)
    logtable.check_decisions(rows, roots.records)
    return decisions


def _with_decisions(ctx: RunContext, inputs: list[Path]) -> list[Path]:
    """*inputs* plus the curated decisions file when it exists."""
    path = paths(ctx).decisions
    return inputs + [path] if path.is_file() else inputs


# ---------------------------------------------------------------------------
# Stages
# ---------------------------------------------------------------------------


def stage_roots(ctx: RunContext) -> int | None:
    """Project the Romanian inputs into the roots JSONL (one pass each)."""
    p = paths(ctx)
    inputs = [p.lexicon, p.verdicts, p.leipzig]

    def body() -> int:
        counts = load_counts(p.leipzig)
        records = build_records(
            lx.read_lexicon(p.lexicon), read_verdicts(p.verdicts), counts
        )
        digests = {ctx.relative(path): sha256_file(path) for path in inputs}
        with ctx.open_output(p.roots, "w", newline="") as fh:
            fh.writelines(roots_lines(records, counts.total_tokens, digests))
        LOG(f"  roots: {len(records):,} noun roots → {ctx.relative(p.roots)}")
        return 0

    return run_stage(
        ctx,
        "roots",
        inputs=inputs,
        params={"version": VERSION},
        outputs=[p.roots],
        body=body,
        log=LOG,
    )


def tp_long_rows(
    roots: RootsFile,
    attest: Attester,
    cache: hz.ObservationCache,
    decisions: logtable.Decisions = logtable.NO_DECISIONS,
) -> list[dict[str, Any]]:
    """Every TP long row of the comparison, in table order."""
    records = roots.records
    rows = hz.cells_rows(records, cache)
    rows += hz.r17_rows(records, cache)
    rows += hz.nde_rows(records, cache)
    rows += logtable.logtable_rows(records, cache, decisions)
    rows += spellout.spellout_rows(records, attest)
    return rows


def stage_tp(ctx: RunContext) -> int | None:
    p = paths(ctx)

    def body() -> int:
        roots, attest, cache = _load(ctx)
        rows = tp_long_rows(roots, attest, cache, load_decisions(ctx, roots))
        written = write_tracked(ctx, p.tp_long, rows, hz.FIELDS)
        LOG(f"  tp: {written:,} rows → {ctx.relative(p.tp_long)}")
        return 0

    return run_stage(
        ctx,
        "tp",
        inputs=_with_decisions(ctx, [p.roots]),
        params={"version": VERSION_2},
        outputs=[p.tp_long],
        body=body,
        log=LOG,
    )


def _workers(ctx: RunContext) -> int:
    """Worker processes of the learner simulation (``workers``)."""
    return int(ctx.params.get("workers", min(8, os.cpu_count() or 1)))


def _learner_params(ctx: RunContext) -> dict[str, Any]:
    samples = ctx.params.get("learner_samples", learners.DEFAULT_SAMPLES)
    v_values = ctx.params.get("learner_V", hz.V_VALUES)
    return {
        "version": VERSION_2,
        "seed": ctx.seed,
        "n_samples": int(samples),
        "V": [int(v) for v in v_values],
        "taus": list(hz.LEARNER_TAUS),
        "evidence": list(hz.LEARNER_EVIDENCE),
    }


def stage_learners(ctx: RunContext) -> int | None:
    p = paths(ctx)
    params = _learner_params(ctx)

    def body() -> int:
        roots, _, cache = _load(ctx)
        rows = hz.learner_rows(
            roots.records,
            cache,
            seed=ctx.seed,
            n_samples=params["n_samples"],
            V_values=params["V"],
            workers=_workers(ctx),
            log=LOG,
        )
        written = write_tracked(ctx, p.learners, rows, hz.LEARNER_FIELDS)
        LOG(f"  learners: {written:,} rows → {ctx.relative(p.learners)}")
        return 0

    return run_stage(
        ctx,
        "learners",
        inputs=[p.roots],
        params=params,
        outputs=[p.learners],
        body=body,
        log=LOG,
    )


def stage_spellout(ctx: RunContext) -> int | None:
    p = paths(ctx)

    def body() -> int:
        roots, attest, _ = _load(ctx)
        rows = spellout.spellout_rows(roots.records, attest)
        write_report(
            ctx,
            p.spellout,
            spellout.render_report(roots.records, rows, attest),
        )
        return 0

    return run_stage(
        ctx,
        "spellout",
        inputs=[p.roots],
        params={"version": VERSION_2},
        outputs=[p.spellout],
        body=body,
        log=LOG,
    )


NOEVIDENCE_FIELDS: Final[tuple[str, ...]] = (
    "tau",
    "root_id",
    "lemma",
    "target_class",
    "stratum",
    "base",
    "group",
    "suffix",
    "outcome",
)


def stage_noevidence(ctx: RunContext) -> int | None:
    p = paths(ctx)

    def body() -> int:
        roots, attest, _ = _load(ctx)
        text = noevidence.render_report(roots.records, attest, p.check_claims)
        write_report(ctx, p.noevidence, text)
        units = [
            {"tau": tau, **_unit_dict(unit)}
            for tau in TAUS
            for unit in noevidence.units(roots.records, tau, attest)
        ]
        write_untracked(ctx, p.noevidence_units, units, NOEVIDENCE_FIELDS)
        return 0

    return run_stage(
        ctx,
        "noevidence",
        inputs=[p.roots, p.check_claims],
        params={"version": VERSION_2},
        outputs=[p.noevidence, p.noevidence_units],
        body=body,
        log=LOG,
    )


def _unit_dict(unit: noevidence.BaseSuffix) -> dict[str, str]:
    return {
        "root_id": unit.root_id,
        "lemma": unit.lemma,
        "target_class": unit.target_class,
        "stratum": unit.stratum,
        "base": unit.base,
        "group": unit.group,
        "suffix": unit.suffix,
        "outcome": unit.outcome,
    }


def stage_covariation(ctx: RunContext) -> int | None:
    """Plural class × verbalizer choice (R23; :mod:`.covariation`)."""
    p = paths(ctx)
    n_iter = int(ctx.params.get("n_permutations", stats.DEFAULT_PERMUTATIONS))

    def body() -> int:
        roots, attest, _ = _load(ctx)
        text = covariation.render_report(
            roots.records, attest, seed=ctx.seed, n_iter=n_iter
        )
        write_report(ctx, p.covariation, text)
        write_untracked(
            ctx,
            p.covariation_units,
            covariation.unit_rows(roots.records, attest),
            covariation.UNIT_FIELDS,
        )
        return 0

    return run_stage(
        ctx,
        "covariation",
        inputs=[p.roots],
        params={
            "version": VERSION_2,
            "seed": ctx.seed,
            "n_permutations": n_iter,
        },
        outputs=[p.covariation, p.covariation_units],
        body=body,
        log=LOG,
    )


def stage_logtable(ctx: RunContext) -> int | None:
    p = paths(ctx)

    def body() -> int:
        roots, attest, _ = _load(ctx)
        decisions = load_decisions(ctx, roots)
        tracked = logtable.table_rows(
            roots.records, full=False, decisions=decisions
        )
        written = write_tracked(ctx, p.logtable, tracked, logtable.LOG_FIELDS)
        full = logtable.table_rows(
            roots.records, full=True, decisions=decisions, attest=attest
        )
        write_untracked(
            ctx,
            p.logtable_full,
            full,
            logtable.LOG_FIELDS + logtable.FULL_EXTRA,
        )
        LOG(f"  logtable: {written} nouns → {ctx.relative(p.logtable)}")
        return 0

    return run_stage(
        ctx,
        "logtable",
        inputs=_with_decisions(ctx, [p.roots]),
        params={"version": VERSION_3},
        outputs=[p.logtable, p.logtable_full],
        body=body,
        log=LOG,
    )


def stage_parity(ctx: RunContext) -> int | None:
    p = paths(ctx)
    inputs = [
        p.lexicon,
        p.roots,
        p.legacy_full,
        p.legacy_downsampled,
        p.tp_long,
    ]

    def body() -> int:
        roots, _, cache = _load(ctx)
        text, ok = parity.render_report(
            lexicon_rows=lx.read_lexicon(p.lexicon),
            records=roots.records,
            cache=cache,
            tp_long=_typed(csvio.read_rows(p.tp_long)),
            full_csv=p.legacy_full,
            downsampled_csv=p.legacy_downsampled,
            inputs={
                ctx.relative(path): sha256_file(path)
                for path in (p.lexicon, p.legacy_full, p.legacy_downsampled)
            },
        )
        write_report(ctx, p.parity, text)
        if not ok:
            LOG(
                "  parity: legacy rows NOT reproduced or reconciliation "
                "unbalanced; see the report"
            )
            return 1
        return 0

    return run_stage(
        ctx,
        "parity",
        inputs=inputs,
        params={"version": VERSION_2},
        outputs=[p.parity],
        body=body,
        log=LOG,
    )


def _typed(rows: Sequence[Mapping[str, str]]) -> list[dict[str, Any]]:
    """TP long CSV rows with the counting columns as ints."""
    out: list[dict[str, Any]] = []
    for row in rows:
        typed: dict[str, Any] = dict(row)
        for name in ("N_target", "N_trigger", "e_nonalt", "n_shielded", "e"):
            typed[name] = int(row[name])
        typed["V"] = int(row["V"]) if row["V"] else None
        out.append(typed)
    return out


def stage_report(ctx: RunContext) -> int | None:
    p = paths(ctx)
    inputs = [
        p.roots,
        p.tp_long,
        p.learners,
        p.parity,
        p.spellout,
        p.noevidence,
        p.covariation,
        p.logtable,
    ]

    def body() -> int:
        body_text = render_readme(
            csvio.read_rows(p.tp_long),
            csvio.read_rows(p.learners),
            p.parity.read_text(encoding="utf-8"),
            csvio.read_rows(p.logtable),
            coverage(read_roots(p.roots)),
        )
        existing = (
            p.readme.read_text(encoding="utf-8")
            if p.readme.is_file()
            else None
        )
        write_report(ctx, p.readme, readme.merge_section(existing, body_text))
        return 0

    return run_stage(
        ctx,
        "report",
        inputs=inputs,
        params={"version": VERSION_3},
        outputs=[p.readme],
        body=body,
        log=LOG,
    )


STAGES: Final[Mapping[str, Stage]] = {
    "roots": stage_roots,
    "tp": stage_tp,
    "learners": stage_learners,
    "spellout": stage_spellout,
    "noevidence": stage_noevidence,
    "covariation": stage_covariation,
    "logtable": stage_logtable,
    "parity": stage_parity,
    "report": stage_report,
}


# ---------------------------------------------------------------------------
# README section
# ---------------------------------------------------------------------------


def _r17_table(rows: Sequence[Mapping[str, str]]) -> list[str]:
    """R17 rows at E0: the rule's own row (``retain`` for a non-trigger
    cell, where the rule predicts no change; ``palatalize`` otherwise)."""
    body: list[list[str]] = []
    for spec in hz.R17_CELLS:
        rule = "retain" if "non-trigger" in spec.memo else "palatalize"
        for subset, tau in (
            ("all", "dictionary"),
            ("all", "ge1"),
            ("core", "dictionary"),
        ):
            row = readme.pick(
                rows,
                table=spec.table,
                target_class=spec.target,
                trigger=spec.trigger,
                subset=subset,
                tau=tau,
                evidence_set="E0",
                coding="dominant",
                V="",
                rule=rule,
            )
            if row is None:
                continue
            body.append(
                [
                    spec.table,
                    spec.target,
                    spec.trigger,
                    rule,
                    subset,
                    tau,
                    readme.fmt(row["N_target"]),
                    readme.fmt(row["N_trigger"]),
                    readme.fmt(row["e"]),
                    readme.fmt(row["n_shielded"]),
                    readme.fmt(row["theta_trigger"]),
                    readme.fmt(row["tolerated_trigger"]),
                    readme.fmt(row["fragility"]),
                ]
            )
    header = [
        "table",
        "target",
        "trigger",
        "rule",
        "subset",
        "τ",
        "N_target",
        "N_trigger",
        "e",
        "shielded",
        "θ",
        "tolerated",
        "fragility",
    ]
    return readme.table(header, body)


def _learner_table(rows: Sequence[Mapping[str, str]]) -> list[str]:
    body: list[list[str]] = []
    for row in rows:
        if (
            row["evidence_set"] != "E0"
            or row["tau"] != "dictionary"
            or row["V"] != "1000"
        ):
            continue
        body.append(
            [
                row["table"],
                row["target_class"],
                row["trigger"],
                readme.fmt(row["share_palatalize"]),
                readme.fmt(row["share_retain"]),
                readme.fmt(row["share_neither"]),
                readme.fmt(row["share_conditioned"]),
            ]
        )
    return readme.table(
        [
            "table",
            "target",
            "trigger",
            "palatalize",
            "retain",
            "neither",
            "conditioned",
        ],
        body,
    )


@dataclass(frozen=True, slots=True)
class Coverage:
    """How much of each evidence set the roots file covers."""

    roots: int
    target_roots: int
    nde_roots: int
    verdict_nouns: int
    readable_verbs: int
    derivatives: int
    derivatives_read: int
    duplicate_lemmas: int


def coverage(roots: RootsFile) -> Coverage:
    records = roots.records
    derivatives = [d for r in records for d in r.derivatives]
    return Coverage(
        roots=len(records),
        target_roots=sum(1 for r in records if r.is_target and not r.nde),
        nde_roots=sum(1 for r in records if r.nde),
        verdict_nouns=sum(1 for r in records if r.verb),
        readable_verbs=sum(1 for r in records if r.verb_cells),
        derivatives=len(derivatives),
        derivatives_read=sum(1 for d in derivatives if d.vowel is not None),
        duplicate_lemmas=sum(1 for r in records if r.duplicates > 1),
    )


def _spellout_line(rows: Sequence[Mapping[str, str]]) -> str:
    parts: list[str] = []
    for target in ("dorsal", "coronal"):
        found = [
            readme.pick(
                rows,
                table="spellout",
                target_class=target,
                subset="all",
                tau=tau,
                coding=coding,
                V="",
                rule="palatalize",
            )
            for tau, coding in (
                ("dictionary", "doublet_exception"),
                ("dictionary", "dominant"),
                ("ge1", "dominant"),
            )
        ]
        if any(row is None for row in found):
            continue
        texts = [
            f"{label}: N={row['N_trigger']}, e={row['e']}, "
            f"θ={readme.fmt(row['theta_trigger'])}"
            for label, row in zip(
                ("dictionary, doublets as exceptions", "dominant", "τ ≥ 1"),
                found,
                strict=True,
            )
            if row is not None
        ]
        parts.append(f"{target} ({'; '.join(texts)})")
    return "; ".join(parts)


def render_readme(
    tp_rows: Sequence[Mapping[str, str]],
    learner_rows: Sequence[Mapping[str, str]],
    parity_text: str,
    log_rows: Sequence[Mapping[str, str]],
    cover: Coverage,
) -> str:
    """The Romanian README section (Markdown)."""
    parity_lines = [
        line.strip()
        for line in parity_text.splitlines()
        if "numeric rows identical" in line
    ]
    classes: dict[str, int] = {}
    lookalikes = 0
    for row in log_rows:
        if row["learned_suffixoid"] != "True":
            lookalikes += 1
            continue
        key = f"{row['semantic_class']} {row['plural_allomorph']}"
        classes[key] = classes.get(key, 0) + 1
    samples = learner_rows[0]["n_samples"] if learner_rows else "0"
    lines = [
        "## Romanian (harmonized, read-only)",
        "",
        "Computed by `lexpipe.langs.romanian_compare` from "
        "`data/romanian_lexicon_with_freq.csv`, "
        "`data/verb_paradigm_verdict.csv` and the Leipzig table "
        "`data/leipzig/freq/ron_news_2024_1M_freq.csv` (all read-only) with "
        "the same core code as Italian and Czech (`lexpipe.tp`, "
        "`lexpipe.attest`, `lexpipe.learners`). The legacy Romanian outputs "
        "are untouched.",
        "",
        "Files:",
        "",
        "- `romanian_harmonized_tp_long.csv`: `lexpipe.tp.TP_LONG_FIELDS` + "
        "`fragility`, `source`, `licence_class`; tables `cells`, `R17a`/"
        "`R17b`/`R17c` (the R17 rule family, with τ, E0–E2, three doublet "
        "codings, V sweep on lemma and paradigm frequency, strata, nested "
        "rules), `nde`, `logtable`, `spellout`.",
        "- `romanian_learners.csv`: simulated learners (R18), "
        f"{samples} vocabularies per V, seed from the run context.",
        "- `romanian_spellout.txt` (R20), `romanian_noevidence.txt` (R21), "
        "`romanian_covariation.txt` (R23: plural class x verbalizer choice), "
        "`romanian_logtable.csv` (R27), `romanian_parity.txt`.",
        "",
        "Definitions: a root is a noun lemma; N_target = target-final roots "
        "of the class (structural NDE classes are separate rows); N_trigger "
        "= roots with a realization before the trigger attested at τ; "
        "shielded = every attested realization avoids the trigger (-uri, or "
        "an -e plural in an -i cell). Subset `all` keeps every outcome-blind "
        "flag (suffix targets, palatal context, suppletive stems, -ică "
        "suppletion) and keeps the legacy outcome-conditioned exclusions "
        "(z→zi 'foreign' plurals, fată/pată) as exceptions; `core` drops the "
        "flags. E1 adds the 2sg and 3sg-subjunctive cells of "
        "`data/verb_paradigm_verdict.csv`, E2 the lexicon's derived verbs "
        "and adjectives, read orthographically (the vowel after the root "
        "decides the cell).",
        "",
        f"Coverage: {cover.roots:,} noun roots, {cover.target_roots:,} "
        f"target-final outside the {cover.nde_roots:,} structural NDE roots; "
        f"E1: {cover.verdict_nouns} nouns with a verdict row, "
        f"{cover.readable_verbs} with a readable cell; E2: "
        f"{cover.derivatives_read:,} of {cover.derivatives:,} listed "
        "derivatives read. One plural per lemma "
        f"({cover.duplicate_lemmas} lemmas have two rows, one without a "
        "plural), so the three doublet codings differ only at E1/E2.",
        "",
        "R17 rows (E0, full lexicon, dominant coding; a non-trigger -e cell "
        "shows the `retain` row, whose exceptions are the alternations):",
        "",
        *_r17_table(tp_rows),
        "",
        "Simulated learners at V = 1000 (E0, τ = dictionary):",
        "",
        *_learner_table(learner_rows),
        "",
        "Spell-out TP (R20; -i vs -a/-ui on alternating-plural roots, "
        "e = roots taking only -a/-ui): " + _spellout_line(tp_rows) + ". "
        "The critique's 60 vs θ(164) = 32.2 counts verb pairs (and -ări, "
        "-iza, -ăni, -arisi), not roots; see `romanian_spellout.txt`.",
        "",
        "-log/-fag/-urg nouns with the learned suffixoid, by class and "
        "plural: "
        + ", ".join(f"{k} {v}" for k, v in sorted(classes.items()))
        + f"; {lookalikes} lookalikes (amurg, milog, ...) listed apart.",
        "",
        "Parity with the legacy R tables: " + "; ".join(parity_lines) + ".",
        "",
        "Sources and caveats: the harvested lexicon's noun rows come from "
        "ro.wiktionary (about 70%) and en.wiktionary (CC BY-SA 4.0 + GFDL), "
        "so aggregate rows name `wiktionary_ro|wiktionary_en` and the -log "
        "table names each row's edition. DEX QC (dexonline, terms "
        "unverified) filled the plurals Wiktionary lacked for some nouns "
        "(lexicon source `| dex_qc`): the aggregates count them, as the "
        "legacy tables do, while the -log table leaves those nouns out "
        "(they are in the untracked full table). The "
        "Leipzig table is lowercased (no case-sensitive counts); stratum "
        "`latin_unresolved` mixes inherited and learned Latin. `ro` is a "
        "comparison-only row of `lexpipe.spec.LANGUAGES`; run the stages "
        "with `python scripts/lang_pipeline.py --lang ro`.",
    ]
    return "\n".join(lines) + "\n"
