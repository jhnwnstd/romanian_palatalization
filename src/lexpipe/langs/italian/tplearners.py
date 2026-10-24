"""Simulated learners of the Italian TP cells (R18, R24).

:func:`lexpipe.learners.simulate_learners` draws vocabularies of V root
types by token frequency (V ∈ {500, 1000, 2000, 5000}) from the R24
frame: every noun root of any final segment and plural type with a
positive frequency. Roots of the cell's class carry their trigger
cells (E0 or E1, scope ``all``: a noun root's adjective forms are
evidence too, R4); the rest of the frame carries none. For each
vocabulary :func:`lexpipe.learners.verdict` says whether
``palatalize``, ``retain``, ``neither`` or (with a nested rule set) the
``conditioned`` grammar is productive.

Learner cells (:data:`LEARNER_SPECS`): R17a (dorsal, -i/-e; nested
"NDE stems retain" over dorsal + NDE roots), R17b (dorsal -i and the
non-trigger -e cell), R17c (coronal -i, E0 only) and the Italian model
fit ``it_stress`` (dorsal -i nested by stress: post-tonic → retain,
elsewhere palatalize). Each at τ ∈ {dictionary, ge1} on the lemma
frequency basis, E0 and E1. The seed comes from the run context (R5);
``n_samples`` from its ``learner_samples`` parameter (default 1000).
A V larger than the frame is skipped (fixture trees).

R26 ("each suffix … counts once as a morpheme in learner evidence"):
the subset ``all`` counts every root, -ico adjectives included, as
the cross-language R17 cells of Romanian and Czech do; subset
``suffix_once`` (R17b -i and it_stress) counts each R26 suffix of the
vocabulary once — the vocabulary's -ico roots become one unit, its
-iaco roots another (:func:`.tptables.r26_unit`), whose status pools
the forms of the members drawn (:func:`.tptables.pooled_unit`) — and
``orth_ico_once`` (it_stress) does the same with every root in
unstressed -ico/-ica (the upper bound: R26 detection is a lower
bound). A unit is never post-tonic: it falls to the elsewhere rule.

Speed: a root's status in a cell (:func:`lexpipe.tp.root_status`) does
not depend on the vocabulary, so every member's status is computed once
per cell and :func:`lexpipe.learners.status_verdict` evaluates a
vocabulary from the status counts (it reproduces
:func:`lexpipe.learners.verdict` exactly). Every cell is evaluated on
the same vocabularies, drawn once (:func:`lexpipe.learners.simulate_all`,
optionally in ``workers`` forked processes; the result does not depend
on their number).
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any, Final

from lexpipe import learners, tp
from lexpipe.langs.italian.roots import Records, RootRecord, observation
from lexpipe.langs.italian.tptables import (
    LANG,
    LICENCE_CLASS,
    SOURCE,
    TARGETS,
    UnitOf,
    cell,
    orth_ico_unit,
    pooled_unit,
    r26_unit,
)

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
LEARNER_TAUS: Final[tuple[str, ...]] = ("dictionary", "ge1")
LEARNER_EVIDENCE: Final[tuple[str, ...]] = ("E0", "E1")
SAMPLES_PARAM: Final[str] = "learner_samples"
V_PARAM: Final[str] = "learner_V"
WORKERS_PARAM: Final[str] = "workers"
BASIS: Final[str] = "lemma"

Member = Callable[[RootRecord], bool]


def _nde_or_dorsal(record: RootRecord) -> bool:
    return record.dorsal or record.target_class.startswith("nde:")


def _nde(record: RootRecord) -> bool:
    return record.target_class.startswith("nde:")


def _post_tonic(record: RootRecord) -> bool:
    return record.velar_post_tonic is True


@dataclass(frozen=True, slots=True)
class LearnerSpec:
    """One learner cell: members, trigger, an optional nested set and
    optional suffix units (R26, module docstring)."""

    table: str
    target: str
    trigger: str
    subset: str
    members: Member
    nested: tuple[tuple[str, str, Member], ...] = ()
    evidence_sets: tuple[str, ...] = LEARNER_EVIDENCE
    units: UnitOf | None = None


LEARNER_SPECS: Final[tuple[LearnerSpec, ...]] = (
    LearnerSpec(
        "R17a",
        "dorsal",
        "-i/-e",
        "with_nde",
        _nde_or_dorsal,
        (("velar_front_lemma", "retain", _nde),),
    ),
    LearnerSpec("R17b", "dorsal", "-i", "all", TARGETS["dorsal"]),
    LearnerSpec("R17b", "dorsal", "-e", "all", TARGETS["dorsal"]),
    LearnerSpec(
        "R17c",
        "coronal",
        "-i",
        "all",
        TARGETS["coronal"],
        evidence_sets=("E0",),
    ),
    LearnerSpec(
        "it_stress",
        "dorsal",
        "-i",
        "all",
        TARGETS["dorsal"],
        (("post_tonic", "retain", _post_tonic),),
    ),
    LearnerSpec(
        "R17b",
        "dorsal",
        "-i",
        "suffix_once",
        TARGETS["dorsal"],
        units=r26_unit,
    ),
    LearnerSpec(
        "it_stress",
        "dorsal",
        "-i",
        "suffix_once",
        TARGETS["dorsal"],
        (("post_tonic", "retain", _post_tonic),),
        units=r26_unit,
    ),
    LearnerSpec(
        "it_stress",
        "dorsal",
        "-i",
        "orth_ico_once",
        TARGETS["dorsal"],
        (("post_tonic", "retain", _post_tonic),),
        units=orth_ico_unit,
    ),
)


def frame_pool(records: Records) -> list[tp.RootObservation]:
    """The R24 frame: noun roots with a positive lemma frequency (no
    cells; statuses are computed separately, :func:`statuses_for`)."""
    return [
        tp.RootObservation(
            record.root_id,
            record.target_class,
            float(record.frequency(BASIS, "noun")),
            {},
        )
        for record in records.roots.values()
        if record.has_noun and record.frequency(BASIS, "noun") > 0
    ]


def statuses_for(
    records: Records, spec: LearnerSpec, evidence: str, spec_cell: tp.CellSpec
) -> dict[str, tp.RootStatus]:
    """Statuses of the spec's member roots in the frame (scope ``all``)."""
    members = [
        observation(record, evidence, scope="all", basis=BASIS)
        for record in records.roots.values()
        if record.has_noun
        and record.frequency(BASIS, "noun") > 0
        and spec.members(record)
    ]
    return statuses_of(members, spec_cell)


def statuses_of(
    pool: Sequence[tp.RootObservation], spec_cell: tp.CellSpec
) -> dict[str, tp.RootStatus]:
    """Each pool root's status in the cell (no_evidence left out)."""
    out: dict[str, tp.RootStatus] = {}
    for root in pool:
        status = tp.root_status(
            root,
            spec_cell.trigger,
            spec_cell.tau,
            spec_cell.coding,
            min_tokens=spec_cell.doublet_min_tokens,
        )
        if status != tp.RootStatus.NO_EVIDENCE:
            out[root.root_id] = status
    return out


def suffix_once_verdict(
    statuses: dict[str, tp.RootStatus],
    units: dict[str, str],
    unit_members: dict[str, tp.RootObservation],
    spec_cell: tp.CellSpec,
    specific: Sequence[learners.StatusRule],
    elsewhere: str = "palatalize",
) -> learners.RuleEval:
    """:func:`lexpipe.learners.status_verdict` with each unit counted
    once (module docstring): the members of a unit drawn into a
    vocabulary are replaced by the unit, whose status is that of
    :func:`.tptables.pooled_unit` over them. *units* maps a member's
    root_id to its unit, *unit_members* to its observation; *statuses*
    holds the other roots'. Evaluations are cached per combination of
    unit statuses (there are few)."""
    unit_roots = {
        unit: tp.RootObservation(unit, "dorsal", 1.0, {})
        for unit in sorted(set(units.values()))
    }
    cache: dict[tuple[tuple[str, tp.RootStatus], ...], learners.RuleEval] = {}

    def evaluate(vocabulary: Sequence[tp.RootObservation]) -> str:
        kept: list[tp.RootObservation] = []
        drawn: dict[str, list[tp.RootObservation]] = {}
        for root in vocabulary:
            unit = units.get(root.root_id)
            if unit is None:
                kept.append(root)
            else:
                drawn.setdefault(unit, []).append(unit_members[root.root_id])
        key = tuple(
            (
                unit,
                tp.root_status(
                    pooled_unit(unit, drawn[unit]),
                    spec_cell.trigger,
                    spec_cell.tau,
                    spec_cell.coding,
                    min_tokens=spec_cell.doublet_min_tokens,
                ),
            )
            for unit in sorted(drawn)
        )
        found = cache.get(key)
        if found is None:
            found = learners.status_verdict(
                {**statuses, **dict(key)}, specific, elsewhere
            )
            cache[key] = found
        return found(kept + [unit_roots[unit] for unit in sorted(drawn)])

    return evaluate


def _rows_for(
    spec: LearnerSpec,
    evidence: str,
    tau: str,
    shares: Sequence[learners.LearnerShares],
    seed: int,
) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for item in shares:
        share = item.shares
        out.append(
            {
                "lang": LANG,
                "table": spec.table,
                "target_class": spec.target,
                "trigger": spec.trigger,
                "subset": spec.subset,
                "tau": tau,
                "evidence_set": evidence,
                "coding": "dominant",
                "freq_basis": BASIS,
                "V": item.V,
                "n_samples": item.n_samples,
                "seed": seed,
                "nested": [name for name, _, _ in spec.nested],
                "share_palatalize": share["palatalize"],
                "share_retain": share["retain"],
                "share_neither": share["neither"],
                "share_conditioned": share["conditioned"],
                "source": SOURCE,
                "licence_class": LICENCE_CLASS,
            }
        )
    return out


@dataclass(frozen=True, slots=True)
class _Config:
    spec: LearnerSpec
    evidence: str
    tau: str
    evaluate: learners.RuleEval

    @property
    def name(self) -> str:
        return (
            f"{self.spec.table} {self.spec.target} {self.spec.trigger} "
            f"{self.spec.subset} {self.evidence} {self.tau}"
        )


def _evaluation(
    records: Records,
    spec: LearnerSpec,
    evidence: str,
    spec_cell: tp.CellSpec,
    specific: Sequence[learners.StatusRule],
) -> learners.RuleEval:
    """The spec's rule evaluation (with its units, if any)."""
    statuses = statuses_for(records, spec, evidence, spec_cell)
    if spec.units is None:
        return learners.status_verdict(statuses, specific, "palatalize")
    units: dict[str, str] = {}
    members: dict[str, tp.RootObservation] = {}
    for record in records.roots.values():
        unit = spec.units(record)
        if (
            unit
            and record.has_noun
            and record.frequency(BASIS, "noun") > 0
            and spec.members(record)
        ):
            units[record.root_id] = unit
            members[record.root_id] = observation(
                record, evidence, scope="all", basis=BASIS
            )
    plain = {k: v for k, v in statuses.items() if k not in units}
    return suffix_once_verdict(plain, units, members, spec_cell, specific)


def _configs(records: Records) -> list[_Config]:
    out: list[_Config] = []
    for spec in LEARNER_SPECS:
        specific = [
            (
                rule,
                frozenset(
                    r.root_id
                    for r in records.roots.values()
                    if spec.members(r)
                    and member(r)
                    and not (spec.units is not None and spec.units(r))
                ),
            )
            for _, rule, member in spec.nested
        ]
        for evidence in spec.evidence_sets:
            for tau in LEARNER_TAUS:
                spec_cell = cell(
                    spec.table,
                    spec.target,
                    spec.trigger,
                    tau=tau,
                    evidence=evidence,
                    basis=BASIS,
                )
                out.append(
                    _Config(
                        spec,
                        evidence,
                        tau,
                        _evaluation(
                            records, spec, evidence, spec_cell, specific
                        ),
                    )
                )
    return out


def learner_rows(
    records: Records,
    *,
    seed: int,
    n_samples: int = learners.DEFAULT_SAMPLES,
    V_values: Sequence[int] = learners.DEFAULT_V_VALUES,
    workers: int = 1,
    log: Callable[[str], None] = print,
) -> list[dict[str, Any]]:
    """Learner shares for every spec × evidence set × τ (module
    docstring)."""
    pool = frame_pool(records)
    usable = [v for v in V_values if v <= len(pool)]
    if not usable:
        log(f"    learners: frame of {len(pool)} roots is too small")
        return []
    configs = _configs(records)
    log(
        f"    learners: {len(configs)} cells x {len(usable)} V x "
        f"{n_samples} samples over {len(pool):,} frame roots, "
        f"{workers} worker(s)"
    )
    shares = learners.simulate_all(
        pool,
        {config.name: config.evaluate for config in configs},
        usable,
        n_samples,
        seed=seed,
        workers=workers,
    )
    rows: list[dict[str, Any]] = []
    for config in configs:
        rows += _rows_for(
            config.spec, config.evidence, config.tau, shares[config.name], seed
        )
    return rows
