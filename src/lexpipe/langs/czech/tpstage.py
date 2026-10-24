"""Stage ``tp``: the three Czech TP analyses and the R17 rows (design §6).

Every number is :mod:`lexpipe.tp`'s (:func:`~lexpipe.tp.two_outcome_rows`,
:func:`~lexpipe.tp.nested_rows`, :func:`~lexpipe.tp.top_n`) or
:mod:`lexpipe.learners`' (simulated learners, fragility): this module
only chooses which roots go into which cell (:mod:`.tpcells`). Roots
are :class:`.roots.RootRecord` types (R4); derivative evidence (E1, E2)
is read from the derivations stage on the ``wikt`` basis (Wiktionary
links and headword candidates: no DeriNet-only pair enters a tracked
table, R9).

Tables (the ``table`` column of ``czech_tp_long.csv``)
-------------------------------------------------------
``rule`` / ``rule_adj``
    (a) rule TP per §6 trigger cell for velars, the velar segments and
    clusters, r (and r after a consonant / a vowel), c, and t/d/n
    (``coronal_inferred``: change inferred from spelling, reported
    apart); the R16 non-trigger rows -é, -em, -ech and the
    non-trigger VOC -e (N = e). Subsets ``all``, ``core`` (R26 suffix
    targets out), ``suffix`` and ``stratum:<s>`` (R22, Wiktionary);
    every τ (R11); both rules (R18).
``competition`` / ``competition_adj``
    (b) per cell and gender (LOC.SG m.in/n, NOM.PL m.an, LOC.PL
    m.an/m.in/n, VOC.SG m.an, adverbs: pooled ``ADV`` and its -y and -o
    competitors apart, ``ADV -y`` / ``ADV -o``): ``palatalize`` =
    P-default, ``retain`` = N-default, for velars, r, t/d/n and the
    baseline (labial, l, s, z; and each apart); every τ × the three
    doublet codings (R19); nested classes ``class:suffixoid`` (the
    whole-word loans in -log/-urg/-graf of R27, which IJP 226 puts with
    -ové), ``class:sky`` (the -ský/-cký adjectives, -ický included, whose
    k is the suffix's: R26, identified morphologically, before the
    etymological classes, since their adverb -y and NOM.PL -čtí/-ští are
    the suffix's behaviour), ``class:suffix`` (the other R26 suffix
    targets: -ák, -ík, -ník, -ek …), ``class:monosyllable``,
    ``class:loan``, ``class:other`` and the Yang-order nested set
    ``nested:*`` whose rule per class is the one with fewer exceptions
    (dominant coding).

    Thresholds (rule, competition and spell-out rows): the R11 τ and
    their strict variants ``ge1_strict``, ``ge0.5pm_strict``,
    ``ge2pm_strict``, which count a cell form's tokens only when the
    form is neither ambiguous with another lexeme nor shared with a
    cell outside its copy group (:data:`.roots.STRICT_TAUS`).
``spellout_loo``
    (c) leave-one-cell-out spell-out TP: subset ``underspecified``
    (alternation evidence outside the tested cell) and
    ``prespecified`` (none), per tested cell × E0/E1/E2 × τ × coding,
    for velars, r and both.
``R17a`` / ``R17b`` / ``R17c`` (and ``_adj``)
    the fixed cross-language rule family (R17) over every realization
    whose ending starts with a vowel of the class (see
    :mod:`.tpcells`); E0/E1/E2 × τ × codings × ``all``/``core``; for
    nouns also the V sweep (R24: V ∈ {500, 1000, 2000, 5000} on the
    lemma and paradigm frequency bases, dominant, τ dictionary/≥1)
    and strata (E0).
``nde``
    NDE classes as separate rows (R16): tautomorphemic dorsals
    (energie, chemie) and hard loan coronals; N = e.
``derivation``
    one trigger row per derivational suffix (R22; evidence set E2):
    P suffixes that palatalize the target, the non-trigger -ista and
    -ismus (N = e) and the shielding N suffixes; velar, r and c bases.

Columns: :data:`lexpipe.tp.TP_LONG_TRACKED_FIELDS` (``fragility``,
``source`` and ``licence_class`` added by :func:`lexpipe.tp.tracked_row`;
the R9 gate).

Other outputs
-------------
``analysis/czech/czech_learners.csv``  simulated learners (R18): 1000
    vocabularies per V by token frequency, seed from the run context;
    every cell on the same draws (:func:`lexpipe.learners.simulate_all`,
    ``workers`` processes, the output independent of their number),
    at τ dictionary, ≥ 1 token and ≥ 1 token strict, on each R24
    frequency basis (``freq_basis``: lemma, paradigm sum,
    paradigm sum over unambiguous forms; :data:`.roots.FREQ_BASES`).
``analysis/czech/czech_tp_summary.txt``  headline tables.
``data/czech/czech_tp_roots.csv``  per-root slots and evidence behind
    the rows (untracked: it rests on per-word Leipzig counts).
"""

from __future__ import annotations

import os
import time
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Final

from lexpipe import csvio, learners, tp
from lexpipe.attest import TAUS
from lexpipe.context import RunContext, run_stage
from lexpipe.langs.czech import leipzig, ortho, roots
from lexpipe.langs.czech import suffixes as sx
from lexpipe.langs.czech import tpcells as tc
from lexpipe.langs.czech import tpsummary
from lexpipe.langs.czech.analysis_layout import AnalysisPaths, analysis_paths
from lexpipe.langs.czech.layout import CzechPaths, paths_for
from lexpipe.langs.czech.roots import ADJ, NOUN, RootRecord

STAGE: Final[str] = "tp"
# 6: core simulate_all/status_verdict and tracked_row; gated writers.
# 7: strict thresholds; ADV -y/-o; class:sky; animacy-specific cells
# without 'm'; -a masculines out of NOM.PL; learners on three bases
# (audit 2026-10-23).
VERSION: Final[str] = "7"
LANG: Final[str] = "cs"
SOURCE: Final[list[str]] = ["wiktionary_en", "wiktionary_cs"]
LICENCE_CLASS: Final[str] = "by-sa"
FIELDS: Final[tuple[str, ...]] = tp.TP_LONG_TRACKED_FIELDS
DERIVATIVE_BASES: Final[tuple[str, ...]] = ("linked", "wikt")
V_VALUES: Final[tuple[int, ...]] = learners.DEFAULT_V_VALUES
FREQ_BASES: Final[tuple[str, ...]] = roots.FREQ_BASES
V_TAUS: Final[tuple[str, ...]] = ("dictionary", "ge1")
LEARNER_TAUS: Final[tuple[str, ...]] = ("dictionary", "ge1", "ge1_strict")
# Thresholds of the rule, competition and spell-out rows.
ANALYSIS_TAUS: Final[tuple[str, ...]] = roots.ALL_TAUS

LEARNER_FIELDS: Final[tuple[str, ...]] = (
    "lang",
    "table",
    "target_class",
    "trigger",
    "subset",
    "tau",
    "evidence_set",
    "coding",
    "freq_basis",
    "V",
    "n_samples",
    "seed",
    "nested",
    "share_palatalize",
    "share_retain",
    "share_neither",
    "share_conditioned",
    "source",
    "licence_class",
)
ROOT_FIELDS: Final[tuple[str, ...]] = (
    "root",
    "lemma",
    "pos",
    "target_class",
    "stem_final",
    "cluster",
    "genders",
    "stratum_wikt",
    "suffix_target",
    "monosyllable",
    "lemma_freq",
    "paradigm_freq",
    *(
        f"slots_{name.replace(' ', '_')}_{tau}"
        for name in tc.SPELLOUT_CELLS
        for tau in V_TAUS
    ),
    *(f"evidence_{e}_{tau}" for e in tc.EVIDENCE_SETS for tau in V_TAUS),
)

Rows = list[dict[str, Any]]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def finish(row: Mapping[str, Any]) -> dict[str, Any]:
    """A core row with fragility, source and licence class
    (:func:`lexpipe.tp.tracked_row`)."""
    return tp.tracked_row(row, source=SOURCE, licence_class=LICENCE_CLASS)


def cell(
    table: str,
    target: str,
    trigger: str,
    *,
    tau: str,
    evidence: str = "E0",
    coding: str = "dominant",
    basis: str = "lemma",
    V: int | None = None,
) -> tp.CellSpec:
    return tp.CellSpec(
        lang=LANG,
        table=table,
        target_class=target,
        trigger=trigger,
        tau=tau,
        evidence_set=evidence,
        coding=coding,
        freq_basis=basis,
        V=V,
    )


def two(
    roots_: Sequence[tp.RootObservation],
    spec: tp.CellSpec,
    subset: str,
    memo: str,
) -> Rows:
    """Both rules over the same roots; empty when no root is in it."""
    if not roots_:
        return []
    return [
        finish(row)
        for row in tp.two_outcome_rows(roots_, spec, subset=subset, memo=memo)
    ]


def observe(
    records: Iterable[RootRecord],
    triggers: Callable[[RootRecord], Mapping[str, Sequence[tp.Observation]]],
    basis: str = "lemma",
) -> list[tp.RootObservation]:
    """RootObservations with only the given triggers."""
    return [
        tp.RootObservation(
            root_id=r.key,
            target_class=r.target_class,
            frequency=float(r.frequency(basis)),
            triggers=triggers(r),
        )
        for r in records
    ]


TriggerFn = Callable[[RootRecord], Mapping[str, Sequence[tp.Observation]]]


def rule_triggers(trigger: tc.RuleTrigger) -> TriggerFn:
    """The (a) observations of one trigger, keyed for the core."""

    def build(record: RootRecord) -> Mapping[str, Sequence[tp.Observation]]:
        return {trigger.key: tc.rule_observations(record, trigger)}

    return build


def competition_triggers(comp: tc.Competition) -> TriggerFn:
    """The (b) observations of one competition, keyed for the core."""

    def build(record: RootRecord) -> Mapping[str, Sequence[tp.Observation]]:
        return {comp.key: tc.competition_observations(record, comp)}

    return build


def table_name(base: str, pos: str) -> str:
    return base if pos == NOUN else f"{base}_adj"


# ---------------------------------------------------------------------------
# (a) Rule TP
# ---------------------------------------------------------------------------


def rule_rows(records: Sequence[RootRecord]) -> Rows:
    """(a) and the R16 non-trigger rows (module docstring)."""
    rows: Rows = []
    subsets = tc.subsets()
    for trigger in tc.RULE_TRIGGERS:
        pool = [r for r in records if r.pos == trigger.pos]
        memo = "(a) R16 non-trigger: N = e" if trigger.nontrigger else "(a)"
        for target, member in tc.RULE_TARGETS.items():
            if not tc.rule_applies(target, trigger):
                continue
            chosen = [r for r in pool if member(r)]
            if not chosen:
                continue
            names = (
                subsets
                if target in tc.MAJOR_TARGETS
                else {"all": subsets["all"]}
            )
            for subset, test in names.items():
                inside = [r for r in chosen if test(r)]
                observed = observe(inside, rule_triggers(trigger))
                for tau in ANALYSIS_TAUS:
                    spec = cell(
                        table_name("rule", trigger.pos),
                        target,
                        trigger.key,
                        tau=tau,
                    )
                    for row in two(observed, spec, subset, memo):
                        row["trigger"] = trigger.name
                        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# (b) Competition TP
# ---------------------------------------------------------------------------

CLASS_ORDER: Final[tuple[tuple[str, Callable[[RootRecord], bool]], ...]] = (
    ("suffixoid", lambda r: r.lemma_suffix in ortho.ORTHOGRAPHIC_SUFFIXES),
    (
        "sky",
        lambda r: r.suffix_target and r.lemma_suffix in ("ský", "cký"),
    ),
    ("suffix", lambda r: r.suffix_target),
    ("monosyllable", lambda r: r.monosyllable),
    ("loan", lambda r: r.loan),
)


def _best_rule(roots_: Sequence[tp.RootObservation], spec: tp.CellSpec) -> str:
    """The rule with fewer exceptions (palatalize on a tie)."""
    counts = tp.domain_counts(roots_, spec)
    pal = counts.exceptions("palatalize")
    ret = counts.exceptions("retain")
    return "retain" if ret < pal else "palatalize"


def nested_competition(
    chosen: Sequence[RootRecord],
    observed: Mapping[str, tp.RootObservation],
    spec: tp.CellSpec,
    memo: str,
) -> Rows:
    """Class rows and the Yang-order nested set (module docstring)."""
    rows: Rows = []
    remaining = list(chosen)
    specific: list[tp.NestedRule] = []
    for name, test in CLASS_ORDER:
        inside = [r for r in remaining if test(r)]
        roots_ = [observed[r.key] for r in inside]
        rows += two(roots_, spec, f"class:{name}", memo)
        if roots_:
            keys = frozenset(r.key for r in inside)
            specific.append(
                tp.NestedRule(
                    f"nested:{name}",
                    _best_rule(roots_, spec),
                    (lambda k: lambda root: root.root_id in k)(keys),
                )
            )
        remaining = [r for r in remaining if not test(r)]
    rest = [observed[r.key] for r in remaining]
    rows += two(rest, spec, "class:other", memo)
    if not chosen:
        return rows
    elsewhere = _best_rule(rest, spec) if rest else "palatalize"
    nested = tp.nested_rows(
        [observed[r.key] for r in chosen],
        spec,
        specific,
        elsewhere=elsewhere,
        elsewhere_subset="nested:elsewhere",
    )
    rows += [finish({**row, "memo": f"{memo}; nested"}) for row in nested]
    return rows


def competition_rows(records: Sequence[RootRecord]) -> Rows:
    """(b) rows (module docstring)."""
    rows: Rows = []
    for comp in tc.COMPETITIONS:
        pool = [r for r in records if r.pos == comp.pos]
        for target, member in tc.COMPETITION_TARGETS.items():
            chosen = [r for r in pool if member(r)]
            if not chosen:
                continue
            a, b = comp.slots(chosen[0].target_class)
            memo = f"(b) A={'/'.join(sorted(a))} B={'/'.join(sorted(b))}"
            if target in ("baseline", "labial", "lateral", "sibilant"):
                memo += " baseline"
            observed_list = observe(chosen, competition_triggers(comp))
            observed = {o.root_id: o for o in observed_list}
            table = table_name("competition", comp.pos)
            for tau in ANALYSIS_TAUS:
                for coding in tp.CODINGS:
                    # The core reads observations by ``trigger``: the
                    # observation key goes in, the cell name comes out.
                    spec = cell(
                        table, target, comp.key, tau=tau, coding=coding
                    )
                    rows += _relabel(
                        two(observed_list, spec, "all", memo), comp
                    )
                spec = cell(table, target, comp.key, tau=tau)
                rows += _relabel(
                    nested_competition(chosen, observed, spec, memo), comp
                )
    return rows


def _relabel(rows: Rows, comp: tc.Competition) -> Rows:
    """Show the cell name, not the observation key, in ``trigger``."""
    for row in rows:
        row["trigger"] = comp.name
    return rows


# ---------------------------------------------------------------------------
# (c) Spell-out, leave one cell out
# ---------------------------------------------------------------------------

SPELLOUT_TARGETS: Final[Mapping[str, frozenset[str]]] = {
    "velar": frozenset({"velar"}),
    "r": frozenset({"r"}),
    "velar+r": frozenset({"velar", "r"}),
}


def root_family(
    records: Mapping[str, RootRecord],
) -> dict[str, list[RootRecord]]:
    """``record key -> [the root's noun record, its adjective record]``
    (learner evidence pools all forms of a root_id, R4)."""
    out: dict[str, list[RootRecord]] = {}
    for key, record in records.items():
        family = [record]
        other = (
            records.get(f"{record.root_id}#adj")
            if record.pos == NOUN
            else records.get(record.root_id)
        )
        if other is not None:
            family.append(other)
        out[key] = family
    return out


def spellout_partition(
    chosen: Sequence[RootRecord],
    family: Mapping[str, list[RootRecord]],
    tau: str,
    evidence: str,
    left_out: frozenset[str],
) -> tuple[list[RootRecord], list[RootRecord]]:
    """(underspecified, prespecified) roots of a tested cell."""
    under: list[RootRecord] = []
    pre: list[RootRecord] = []
    for record in chosen:
        if tc.alternation_evidence(
            family[record.key], tau, evidence, left_out
        ):
            under.append(record)
        else:
            pre.append(record)
    return under, pre


def spellout_rows(
    records: Sequence[RootRecord],
    family: Mapping[str, list[RootRecord]],
) -> Rows:
    """(c) rows (module docstring)."""
    rows: Rows = []
    memo = "(c) LOO"
    for name, left_out in tc.SPELLOUT_CELLS.items():
        comp = tc.COMPETITION_BY_NAME[name]
        for target, classes in SPELLOUT_TARGETS.items():
            chosen = [
                r
                for r in records
                if r.pos == NOUN and r.target_class in classes
            ]
            observed = {
                o.root_id: o
                for o in observe(chosen, competition_triggers(comp))
            }
            for evidence in tc.EVIDENCE_SETS:
                for tau in ANALYSIS_TAUS:
                    under, pre = spellout_partition(
                        chosen, family, tau, evidence, left_out
                    )
                    for coding in tp.CODINGS:
                        spec = tp.CellSpec(
                            lang=LANG,
                            table="spellout_loo",
                            target_class=target,
                            trigger=comp.key,
                            tau=tau,
                            evidence_set=evidence,
                            coding=coding,
                        )
                        for subset, part in (
                            ("underspecified", under),
                            ("prespecified", pre),
                        ):
                            rows += _relabel(
                                two(
                                    [observed[r.key] for r in part],
                                    spec,
                                    subset,
                                    memo,
                                ),
                                comp,
                            )
    return rows


# ---------------------------------------------------------------------------
# R17
# ---------------------------------------------------------------------------


def _r17_observed(
    chosen: Sequence[RootRecord],
    definition: tc.R17Def,
    evidence: str,
    basis: str = "lemma",
) -> list[tp.RootObservation]:
    def triggers(record: RootRecord) -> Mapping[str, Sequence[tp.Observation]]:
        if definition.derivatives_only:
            return {
                definition.key: tc.r17_derivative_observations(
                    record, definition.onsets
                )
            }
        return {
            definition.key: tc.r17_observations(
                record, definition.onsets, evidence
            )
        }

    return observe(chosen, triggers, basis)


def _r17_spec(
    definition: tc.R17Def,
    pos: str,
    *,
    tau: str,
    evidence: str,
    coding: str = "dominant",
    basis: str = "lemma",
    V: int | None = None,
) -> tp.CellSpec:
    return tp.CellSpec(
        lang=LANG,
        table=table_name(definition.table, pos),
        target_class=definition.target,
        trigger=definition.key,
        tau=tau,
        evidence_set=evidence,
        coding=coding,
        freq_basis=basis,
        V=V,
    )


def _label(rows: Rows, definition: tc.R17Def) -> Rows:
    for row in rows:
        row["trigger"] = definition.trigger
    return rows


def frame(records: Sequence[RootRecord], V: int, basis: str) -> frozenset[str]:
    """Keys of the V most frequent noun roots of any final (R24)."""
    table = [
        {"lemma": r.key, "freq": r.frequency(basis)}
        for r in records
        if r.pos == NOUN
    ]
    return frozenset(str(row["lemma"]) for row in tp.top_n(table, V, "freq"))


def r17_rows(
    records: Sequence[RootRecord], v_values: Sequence[int] = V_VALUES
) -> Rows:
    """R17 rows (module docstring)."""
    rows: Rows = []
    subsets = tc.subsets()
    frames = {
        (basis, V): frame(records, V, basis)
        for basis in FREQ_BASES
        for V in v_values
    }
    for definition in tc.R17_DEFS:
        member = tc.R17_TARGETS[definition.target]
        for pos in (NOUN, ADJ):
            chosen = [r for r in records if r.pos == pos and member(r)]
            if not chosen:
                continue
            for evidence in tc.EVIDENCE_SETS:
                if definition.derivatives_only and evidence != "E2":
                    continue
                for subset in ("all", "core"):
                    inside = [r for r in chosen if subsets[subset](r)]
                    observed = _r17_observed(inside, definition, evidence)
                    for tau in TAUS:
                        for coding in tp.CODINGS:
                            spec = _r17_spec(
                                definition,
                                pos,
                                tau=tau,
                                evidence=evidence,
                                coding=coding,
                            )
                            rows += _label(
                                two(observed, spec, subset, definition.memo),
                                definition,
                            )
            if pos == NOUN and not definition.derivatives_only:
                rows += _r17_sweeps(
                    chosen, definition, frames, subsets, v_values
                )
    return rows


def _r17_sweeps(
    chosen: Sequence[RootRecord],
    definition: tc.R17Def,
    frames: Mapping[tuple[str, int], frozenset[str]],
    subsets: Mapping[str, Callable[[RootRecord], bool]],
    v_values: Sequence[int],
) -> Rows:
    rows: Rows = []
    for basis in FREQ_BASES:
        observed = _r17_observed(chosen, definition, "E0", basis)
        for V in v_values:
            frame_keys = frames[(basis, V)]
            framed = [o for o in observed if o.root_id in frame_keys]
            for tau in V_TAUS:
                spec = _r17_spec(
                    definition, NOUN, tau=tau, evidence="E0", basis=basis, V=V
                )
                rows += _label(
                    two(framed, spec, "all", definition.memo), definition
                )
    for name in tc.STRATA:
        subset = f"stratum:{name}"
        inside = [r for r in chosen if subsets[subset](r)]
        observed = _r17_observed(inside, definition, "E0")
        for tau in TAUS:
            spec = _r17_spec(definition, NOUN, tau=tau, evidence="E0")
            rows += _label(
                two(observed, spec, subset, definition.memo), definition
            )
    return rows


# ---------------------------------------------------------------------------
# NDE
# ---------------------------------------------------------------------------


DERIVATION_TARGETS: Final[tuple[str, ...]] = ("velar", "r", "c")


def derivation_observations(
    record: RootRecord, suffix: sx.Suffix
) -> tuple[tp.Observation, ...]:
    """One root's derivatives with *suffix* (R22 trigger rows): a P
    suffix that palatalizes the target is a trigger (``alternates`` =
    pal); an NT suffix (-ista, -ismus) is a non-trigger front vowel
    (N = e); an N suffix (back vowel) shields the root."""
    out: list[tp.Observation] = []
    for deriv in record.derivatives:
        if deriv.suffix != suffix.label:
            continue
        if suffix.klass == sx.P_CLASS:
            if not deriv.palatalizes_base:
                continue
            alternates: bool | None = deriv.outcome == sx.PAL
            allomorph = "P"
        elif suffix.klass == sx.NT_CLASS:
            alternates, allomorph = deriv.outcome == sx.PAL, "P"
        else:
            alternates, allomorph = None, "N"
        out.append(
            tp.Observation(
                alternates=alternates,
                attested_by_tau=deriv.attested_by_tau(),
                allomorph=allomorph,
                form=deriv.derivative,
                count=deriv.count,
            )
        )
    return tuple(out)


def derivation_rows(records: Sequence[RootRecord]) -> Rows:
    """One trigger row per derivational suffix (R22): N = bases with
    an attested derivative, e = derivatives without the change; -ista
    and -ismus are non-trigger rows (tankista: N = e)."""
    rows: Rows = []
    nouns = [r for r in records if r.pos == NOUN]
    for suffix in sx.SUFFIXES:
        key = f"deriv:{suffix.label}"
        memo = {
            sx.P_CLASS: "R22 derivational trigger",
            sx.NT_CLASS: "R22 non-trigger suffix: N = e",
            sx.N_CLASS: "R22 non-palatalizing suffix: shields",
        }[suffix.klass]

        def triggers(
            record: RootRecord, s: sx.Suffix = suffix, k: str = key
        ) -> Mapping[str, Sequence[tp.Observation]]:
            return {k: derivation_observations(record, s)}

        for target in DERIVATION_TARGETS:
            chosen = [r for r in nouns if r.target_class == target]
            observed = observe(chosen, triggers)
            for tau in TAUS:
                spec = cell("derivation", target, key, tau=tau, evidence="E2")
                for row in two(observed, spec, "all", memo):
                    row["trigger"] = suffix.label
                    rows.append(row)
    return rows


def nde_rows(records: Sequence[RootRecord]) -> Rows:
    """Structural NDE classes, N = e (R16)."""
    rows: Rows = []
    for pos in (NOUN, ADJ):
        for name in tc.NDE_CLASSES:
            chosen = [
                r for r in records if r.pos == pos and r.nde_class == name
            ]
            observed = observe(
                chosen, lambda r: {"nde": tc.nde_observations(r)}
            )
            for tau in TAUS:
                spec = tp.CellSpec(
                    lang=LANG,
                    table=table_name("nde", pos),
                    target_class=f"nde:{name}",
                    trigger="nde",
                    tau=tau,
                )
                rows += [
                    {**row, "trigger": "root-internal front vowel"}
                    for row in two(observed, spec, "all", "structural NDE")
                ]
    return rows


# ---------------------------------------------------------------------------
# Simulated learners (R18)
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class LearnerCell:
    """One learner configuration: which roots carry which trigger."""

    table: str
    target: str
    trigger: str
    subset: str
    key: str
    member: Callable[[RootRecord], bool]
    observations: Callable[[RootRecord], Sequence[tp.Observation]]


def learner_cells(
    family: Mapping[str, list[RootRecord]], tau: str
) -> list[LearnerCell]:
    """R17 cells, three competitions (velars and baseline) and two
    spell-out cells (underspecified velars at τ, E0)."""
    cells: list[LearnerCell] = []
    for definition in tc.R17_DEFS:
        if definition.derivatives_only or definition.trigger == "-e/-i/-y":
            continue
        if definition.table == "R17a" and definition.trigger != "-e/-i":
            continue
        if (
            definition.table == "R17c"
            and definition.trigger == "-e"
            and definition.target != "r"
        ):
            continue  # s, z, l, t/d/n never change before -e: retain
        cells.append(
            LearnerCell(
                definition.table,
                definition.target,
                definition.trigger,
                "all",
                definition.key,
                tc.R17_TARGETS[definition.target],
                (lambda d: lambda r: tc.r17_observations(r, d.onsets, "E0"))(
                    definition
                ),
            )
        )
    for name in ("LOC.SG m.in", "NOM.PL m.an", "LOC.PL m.in"):
        comp = tc.COMPETITION_BY_NAME[name]
        for target in ("velar", "baseline"):
            cells.append(
                LearnerCell(
                    "competition",
                    target,
                    name,
                    "all",
                    comp.key,
                    tc.COMPETITION_TARGETS[target],
                    (lambda c: lambda r: tc.competition_observations(r, c))(
                        comp
                    ),
                )
            )
    for name in ("LOC.SG m.in", "NOM.PL m.an"):
        comp = tc.COMPETITION_BY_NAME[name]
        left_out = tc.SPELLOUT_CELLS[name]

        def under(record: RootRecord, lo: frozenset[str] = left_out) -> bool:
            return record.target_class == "velar" and tc.alternation_evidence(
                family[record.key], tau, "E0", lo
            )

        cells.append(
            LearnerCell(
                "spellout_loo",
                "velar",
                name,
                "underspecified",
                comp.key,
                under,
                (lambda c: lambda r: tc.competition_observations(r, c))(comp),
            )
        )
    return cells


def _learner_row(
    item: LearnerCell,
    tau: str,
    share: learners.LearnerShares,
    seed: int,
    basis: str = "lemma",
) -> dict[str, Any]:
    values = share.shares
    return {
        "lang": LANG,
        "table": item.table,
        "target_class": item.target,
        "trigger": item.trigger,
        "subset": item.subset,
        "tau": tau,
        "evidence_set": "E0",
        "coding": "dominant",
        "freq_basis": basis,
        "V": share.V,
        "n_samples": share.n_samples,
        "seed": seed,
        "nested": "",
        "share_palatalize": values["palatalize"],
        "share_retain": values["retain"],
        "share_neither": values["neither"],
        "share_conditioned": values["conditioned"],
        "source": list(SOURCE),
        "licence_class": LICENCE_CLASS,
    }


def learner_pool(
    records: Sequence[RootRecord], basis: str = "lemma"
) -> list[tp.RootObservation]:
    """The frame every learner draws from: the noun roots with a
    positive frequency on *basis* (without cells: the cells' statuses
    are computed once, :func:`learner_configs`)."""
    return [
        tp.RootObservation(
            r.key, r.target_class, float(r.frequency(basis)), {}
        )
        for r in records
        if r.pos == NOUN and r.frequency(basis) > 0
    ]


@dataclass(frozen=True, slots=True)
class LearnerConfig:
    """One learner cell at one τ with its rule evaluation."""

    name: str
    item: LearnerCell
    tau: str
    evaluate: learners.RuleEval


def learner_configs(
    records: Sequence[RootRecord],
    family: Mapping[str, list[RootRecord]],
) -> list[LearnerConfig]:
    """Every learner cell × τ, evaluated from precomputed root statuses
    (:func:`lexpipe.learners.status_verdict`; the members' statuses in
    the cell, every other frame root has no evidence there). Statuses
    do not depend on the frequency basis, so every noun root that is a
    member gets one."""
    out: list[LearnerConfig] = []
    for tau in LEARNER_TAUS:
        for index, item in enumerate(learner_cells(family, tau)):
            members = [
                tp.RootObservation(
                    r.key,
                    r.target_class,
                    float(r.lemma_freq),
                    {item.key: item.observations(r)},
                )
                for r in records
                if r.pos == NOUN and item.member(r)
            ]
            spec = tp.CellSpec(
                lang=LANG,
                table=item.table,
                target_class=item.target,
                trigger=item.key,
                tau=tau,
            )
            out.append(
                LearnerConfig(
                    f"{tau} {index}",
                    item,
                    tau,
                    learners.status_verdict(tp.cell_statuses(members, spec)),
                )
            )
    return out


def learner_rows(
    records: Sequence[RootRecord],
    family: Mapping[str, list[RootRecord]],
    *,
    seed: int,
    n_samples: int,
    workers: int = 1,
    v_values: Sequence[int] = V_VALUES,
    log: Callable[[str], None] = print,
) -> Rows:
    """Learner shares per cell, τ and V (lemma frequency, E0).

    Every cell is evaluated on the same simulated vocabularies, drawn
    once from the frame (:func:`lexpipe.learners.simulate_all`; with
    *workers* > 1 the evaluations run in forked processes). Each
    cell's shares equal a separate
    :func:`lexpipe.learners.simulate_learners` call with the seed, and
    the output does not depend on *workers* (R5).
    """
    configs = learner_configs(records, family)
    rows: Rows = []
    for basis in FREQ_BASES:
        pool = learner_pool(records, basis)
        shares = learners.simulate_all(
            pool,
            {config.name: config.evaluate for config in configs},
            v_values,
            n_samples,
            seed=seed,
            workers=workers,
        )
        log(
            f"    learners ({basis}): {len(configs)} configurations on "
            f"{len(pool):,} frame roots, {workers} worker(s)"
        )
        rows += [
            _learner_row(config.item, config.tau, share, seed, basis)
            for config in configs
            for share in shares[config.name]
        ]
    return rows


# ---------------------------------------------------------------------------
# Per-root file
# ---------------------------------------------------------------------------


def _slots(record: RootRecord, comp: tc.Competition, tau: str) -> str:
    chosen = set()
    for obs in tc.competition_observations(record, comp):
        if obs.attested_by_tau.get(tau, False):
            chosen.add("A" if obs.alternates else "B")
    return "".join(sorted(chosen))


def root_rows(
    records: Sequence[RootRecord], family: Mapping[str, list[RootRecord]]
) -> Rows:
    """Per-root slots and evidence of the velar and r noun roots."""
    out: Rows = []
    for record in records:
        if record.pos != NOUN or record.target_class not in ("velar", "r"):
            continue
        row: dict[str, Any] = {
            "root": record.key,
            "lemma": record.lemma,
            "pos": record.pos,
            "target_class": record.target_class,
            "stem_final": record.stem_final,
            "cluster": record.cluster,
            "genders": record.genders(),
            "stratum_wikt": record.stratum,
            "suffix_target": record.suffix_target,
            "monosyllable": record.monosyllable,
            "lemma_freq": record.lemma_freq,
            "paradigm_freq": record.paradigm_freq,
        }
        for name in tc.SPELLOUT_CELLS:
            comp = tc.COMPETITION_BY_NAME[name]
            for tau in V_TAUS:
                row[f"slots_{name.replace(' ', '_')}_{tau}"] = _slots(
                    record, comp, tau
                )
        for evidence in tc.EVIDENCE_SETS:
            for tau in V_TAUS:
                row[f"evidence_{evidence}_{tau}"] = tc.alternation_evidence(
                    family[record.key], tau, evidence, frozenset()
                )
        out.append(row)
    return out


# ---------------------------------------------------------------------------
# Stage
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class TpSettings:
    """Run parameters of the stage (``ctx.params``): learner samples,
    worker processes, and the V values of the sweep and the learners
    (fixture trees have far fewer than 5,000 roots)."""

    n_samples: int
    workers: int
    v_values: tuple[int, ...]

    @classmethod
    def from_context(cls, ctx: RunContext) -> TpSettings:
        params = ctx.params
        return cls(
            n_samples=int(params.get("n_samples", learners.DEFAULT_SAMPLES)),
            workers=int(params.get("workers", min(8, os.cpu_count() or 1))),
            v_values=tuple(int(v) for v in params.get("v_values", V_VALUES)),
        )


def run(ctx: RunContext) -> int | None:
    """Stage entry point (see the module docstring)."""
    paths = paths_for(ctx)
    apaths = analysis_paths(ctx)
    words = leipzig.ensure_words_file(ctx, paths)
    settings = TpSettings.from_context(ctx)
    return run_stage(
        ctx,
        STAGE,
        inputs=[
            paths.lexicon_complete,
            paths.cells,
            paths.en_entries,
            words,
            apaths.lexeme_strata,
            apaths.derivations,
        ],
        params={
            "version": VERSION,
            "seed": ctx.seed,
            "n_samples": settings.n_samples,
            "v_values": list(settings.v_values),
            "derivative_bases": list(DERIVATIVE_BASES),
        },
        outputs=[
            apaths.tp_long,
            apaths.learners,
            apaths.tp_summary,
            apaths.tp_roots,
        ],
        body=lambda: _tp(ctx, paths, apaths, settings),
    )


def load(
    paths: CzechPaths, apaths: AnalysisPaths, ctx: RunContext
) -> dict[str, RootRecord]:
    """Root records with cells, derivatives (wikt basis) and adverbs."""
    strata = roots.load_strata(apaths.lexeme_strata)
    records, lexeme_key = roots.load_records(paths.lexicon_complete, strata)
    roots.add_cells(records, lexeme_key, paths.cells)
    roots.add_derivatives(
        records, lexeme_key, apaths.derivations, DERIVATIVE_BASES
    )
    table = leipzig.attestation_table(leipzig.ensure_words_file(ctx, paths))
    roots.add_adverbs(records, roots.adverb_links(paths.en_entries), table)
    return records


def sort_key(row: Mapping[str, Any]) -> tuple[Any, ...]:
    """Canonical order of the long rows."""
    v = row.get("V")
    return (
        str(row["table"]),
        str(row["target_class"]),
        str(row["trigger"]),
        str(row["subset"]),
        ANALYSIS_TAUS.index(str(row["tau"])),
        str(row["evidence_set"]),
        tp.CODINGS.index(str(row["coding"])),
        str(row["freq_basis"]),
        -1 if v is None else int(v),
        str(row["rule"]),
        str(row["memo"]),
    )


def _tp(
    ctx: RunContext,
    paths: CzechPaths,
    apaths: AnalysisPaths,
    settings: TpSettings,
) -> int:
    records_map = load(paths, apaths, ctx)
    records = roots.ordered(records_map)
    family = root_family(records_map)
    print(f"  root records: {len(records):,}", flush=True)
    rows: Rows = []
    parts: tuple[tuple[str, Callable[[], Rows]], ...] = (
        ("rule", lambda: rule_rows(records)),
        ("competition", lambda: competition_rows(records)),
        ("spellout_loo", lambda: spellout_rows(records, family)),
        ("R17", lambda: r17_rows(records, settings.v_values)),
        ("nde", lambda: nde_rows(records)),
        ("derivation", lambda: derivation_rows(records)),
    )
    for name, build in parts:
        start = time.perf_counter()
        part = build()
        elapsed = time.perf_counter() - start
        print(f"  {name}: {len(part):,} rows ({elapsed:.1f}s)", flush=True)
        rows += part
    rows.sort(key=sort_key)
    csvio.write_tracked(apaths.tp_long, rows, FIELDS, ctx=ctx)
    learner = learner_rows(
        records,
        family,
        seed=ctx.seed,
        n_samples=settings.n_samples,
        workers=settings.workers,
        v_values=settings.v_values,
    )
    csvio.write_tracked(apaths.learners, learner, LEARNER_FIELDS, ctx=ctx)
    csvio.write_rows(
        apaths.tp_roots, root_rows(records, family), ROOT_FIELDS, ctx=ctx
    )
    text = tpsummary.render(rows, learner, records)
    ctx.write_tracked_text(apaths.tp_summary, text, sources=SOURCE)
    print(text, flush=True)
    return 0
