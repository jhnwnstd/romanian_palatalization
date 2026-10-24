"""Harmonized Romanian TP rows, V sweep, nesting, learners (R16–R19, R24).

Every number comes from the core: :func:`lexpipe.tp.two_outcome_rows`
and :func:`lexpipe.tp.nested_rows` (R16, R18), :func:`lexpipe.tp.top_n`
(the V frames), :func:`lexpipe.learners.simulate_all` /
:func:`lexpipe.learners.status_verdict` (simulated learners: every cell
on the same draws, from precomputed root statuses) and
:func:`lexpipe.tp.tracked_row` (fragility and the tracked columns).
This module only chooses which roots go into which cell.

Tables (the ``table`` column)
-----------------------------
``cells``
    the descriptive grid: target c, g, dorsal (c+g), t, d, s,
    coronal (t+d+s) and z × trigger -i, -e, -i/-e × subset
    ``all`` / ``core`` / ``with_nde``, every τ, E0, dominant coding.
``R17a`` / ``R17b`` / ``R17c``
    the fixed cross-language rule family (R17):
    ``R17a`` dorsal / _[+front V] (trigger -i/-e);
    ``R17b`` dorsal / _[+high,+front] (trigger -i), with its
    non-trigger cell -e (Romanian dorsals palatalize there anyway);
    ``R17c`` coronal / _[+high,+front] for t+d+s and for z
    (separately), each with its non-trigger cell -e (N = e expected),
    and the rejected pooled coronal rule (-i/-e). Subsets ``all`` and
    ``core`` × every τ × E0, E1, E2 × the three doublet codings (R19);
    V ∈ {500, 1000, 2000, 5000} on the lemma and paradigm frequency
    bases (R24; subset ``all``, dominant coding); strata (E0); and
    nested rule sets in Yang's order (subset ``nested:<name>``,
    elsewhere ``nested:elsewhere``).
``nde``
    the structural NDE classes as separate rows (R16): ochi, paduchi,
    frontstem (lexical palatals) and coronal_frontstem; E0.

Subsets
-------
``all``       target-final roots of the class, structural NDE classes
              excluded (they have their own rows);
``core``      ``all`` minus every outcome-blind flag (suffix_target,
              palatal_context, suppletive_stem, ica_suppletion): the
              outcome-blind analogue of the legacy TP domain;
``with_nde``  ``all`` plus the class's NDE roots (legacy "WITH NDEB");
``stratum:<s>``  ``all`` restricted to one R22 stratum.

The V frame (R24) is the V most frequent noun roots of any final
segment and plural type (:func:`lexpipe.tp.top_n`, ties by root_id,
zero-frequency roots never chosen); a cell then counts its class's
roots inside the frame. The legacy top-1000 frame (nouns with an
-i/-e plural and a target final only) is reproduced in the parity
report, not here.

Memos
-----
``memo`` names the rule a row evaluates: ``dorsal/_[+front V]`` (R17a),
``dorsal/_[+high,+front]`` (R17b), ``coronal/_[+high,+front]`` (R17c;
``, z`` for z). ``: non-trigger cell`` marks an -e cell under a
[+high,+front] rule: there the rule predicts no change, so its
exceptions are the ``retain`` row's (N = e expected for coronals,
except sc → ște, which the account derives through the dorsal rule;
Romanian dorsals do palatalize before -e). ``coronal/_[+front V]
(pooled; not in R17)`` is the rejected single coronal rule, nested as
"retain before the -e plural, palatalize elsewhere". Nested rows add
``; nested in <subset>``.

Extra columns
-------------
The CSV has :data:`lexpipe.tp.TP_LONG_TRACKED_FIELDS`:
:data:`lexpipe.tp.TP_LONG_FIELDS` plus ``fragility`` (the fewest item
changes that flip ``tolerated_trigger``), ``source`` and
``licence_class`` (the R9 gate of :func:`lexpipe.csvio.write_tracked`).
``source`` is ``wiktionary_ro|wiktionary_en``: the harvested lexicon's
rows come from both editions.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Final

from lexpipe import learners, tp
from lexpipe.attest import TAUS
from lexpipe.langs.romanian_compare import morphology as mo
from lexpipe.langs.romanian_compare.lexicon import STRATA
from lexpipe.langs.romanian_compare.observations import (
    EVIDENCE_SETS,
    root_observation,
)
from lexpipe.langs.romanian_compare.records import RootRecord

LANG: Final[str] = "ro"
# The harvested lexicon's rows come from ro.wiktionary (about 70% of the
# nouns) and en.wiktionary (lexicon.lexicon_source); aggregates rest on
# both (R9). DEX QC filled the plurals Wiktionary lacked for some nouns
# (source "… | dex_qc"; dexonline's terms are unverified): the
# aggregates count those roots like the legacy Romanian tables do, but
# per-noun tables leave them out (logtable).
SOURCE: Final[tuple[str, ...]] = ("wiktionary_ro", "wiktionary_en")
LICENCE_CLASS: Final[str] = "by-sa"
FIELDS: Final[tuple[str, ...]] = tp.TP_LONG_TRACKED_FIELDS
V_VALUES: Final[tuple[int, ...]] = learners.DEFAULT_V_VALUES
FREQ_BASES: Final[tuple[str, ...]] = ("lemma", "paradigm")

Member = Callable[[RootRecord], bool]
Attest = Callable[[str], Mapping[str, bool]]


def _seg(*segments: str) -> Member:
    def member(record: RootRecord) -> bool:
        return record.is_target and record.seg in segments

    return member


def _klass(name: str) -> Member:
    def member(record: RootRecord) -> bool:
        return record.target_class == name

    return member


TARGETS: Final[Mapping[str, Member]] = {
    "c": _seg("c"),
    "g": _seg("g"),
    "dorsal": _klass("dorsal"),
    "t": _seg("t"),
    "d": _seg("d"),
    "s": _seg("s"),
    "coronal": _klass("coronal"),
    "z": _klass("z"),
}


def _no_nde(record: RootRecord) -> bool:
    return not record.nde


def _core(record: RootRecord) -> bool:
    return not record.nde and not record.flags


def _any(record: RootRecord) -> bool:
    return True


def _stratum(name: str) -> Member:
    def member(record: RootRecord) -> bool:
        return not record.nde and record.stratum == name

    return member


SUBSETS: Final[Mapping[str, Member]] = {
    "all": _no_nde,
    "core": _core,
    "with_nde": _any,
    **{f"stratum:{name}": _stratum(name) for name in STRATA},
}


@dataclass(frozen=True, slots=True)
class NestedSpec:
    """A specific rule of a nested set (members by root properties)."""

    name: str
    rule: str
    member: Member


@dataclass(frozen=True, slots=True)
class CellDef:
    """One harmonized cell: target class × trigger, with its memo."""

    table: str
    target: str
    trigger: str
    memo: str
    nested: tuple[NestedSpec, ...] = ()
    nested_subset: str = "all"


def _velar_front(record: RootRecord) -> bool:
    return record.nde in ("ochi", "paduchi")


def _has_flag(flag: str) -> Member:
    def member(record: RootRecord) -> bool:
        return flag in record.flags

    return member


def _e_plural(record: RootRecord) -> bool:
    return record.plural_allomorph == "-e"


def _borrowed(record: RootRecord) -> bool:
    return record.stratum == "borrowed"


R17_CELLS: Final[tuple[CellDef, ...]] = (
    CellDef(
        "R17a",
        "dorsal",
        "-i/-e",
        "dorsal/_[+front V]",
        (NestedSpec("velar_front_lemma", "retain", _velar_front),),
        "with_nde",
    ),
    CellDef(
        "R17b",
        "dorsal",
        "-i",
        "dorsal/_[+high,+front]",
    ),
    CellDef(
        "R17b",
        "dorsal",
        "-e",
        "dorsal/_[+high,+front]: non-trigger cell",
    ),
    CellDef(
        "R17c",
        "coronal",
        "-i",
        "coronal/_[+high,+front]",
        (
            NestedSpec(
                "palatal_context", "retain", _has_flag("palatal_context")
            ),
        ),
    ),
    CellDef(
        "R17c",
        "coronal",
        "-e",
        "coronal/_[+high,+front]: non-trigger cell",
    ),
    CellDef(
        "R17c",
        "coronal",
        "-i/-e",
        "coronal/_[+front V] (pooled; not in R17)",
        (NestedSpec("e_plural", "retain", _e_plural),),
    ),
    CellDef(
        "R17c",
        "z",
        "-i",
        "coronal/_[+high,+front], z",
        (NestedSpec("borrowed", "retain", _borrowed),),
    ),
    CellDef(
        "R17c",
        "z",
        "-e",
        "coronal/_[+high,+front], z: non-trigger cell",
    ),
)


def _nde(name: str) -> Member:
    def member(record: RootRecord) -> bool:
        return record.nde == name

    return member


NDE_TARGETS: Final[Mapping[str, Member]] = {
    f"nde:{name}": _nde(name)
    for name in ("ochi", "paduchi", "frontstem", "coronal_frontstem")
}


# ---------------------------------------------------------------------------
# Row helpers
# ---------------------------------------------------------------------------


def finish(row: Mapping[str, Any]) -> dict[str, Any]:
    """A core TP long row with the tracked columns
    (:func:`lexpipe.tp.tracked_row`)."""
    return tp.tracked_row(row, source=SOURCE, licence_class=LICENCE_CLASS)


class ObservationCache:
    """RootObservations per (evidence set, frequency basis), built once."""

    def __init__(self, records: Sequence[RootRecord], attest: Attest) -> None:
        self._records = records
        self._attest = attest
        self._cache: dict[tuple[str, str], dict[str, tp.RootObservation]] = {}

    def get(
        self, evidence: str, basis: str = "lemma"
    ) -> Mapping[str, tp.RootObservation]:
        key = (evidence, basis)
        if key not in self._cache:
            self._cache[key] = {
                r.root_id: root_observation(
                    r, evidence, self._attest, basis=basis
                )
                for r in self._records
            }
        return self._cache[key]


def select(
    records: Iterable[RootRecord], target: Member, subset: Member
) -> list[RootRecord]:
    """Records of a target class within a subset (input order)."""
    return [r for r in records if target(r) and subset(r)]


def frame(records: Sequence[RootRecord], V: int, basis: str) -> frozenset[str]:
    """root_ids of the V most frequent noun roots (R24; core top_n)."""
    rows = [{"lemma": r.root_id, "freq": r.frequency(basis)} for r in records]
    return frozenset(str(row["lemma"]) for row in tp.top_n(rows, V, "freq"))


def _cell(
    table: str,
    target: str,
    trigger: str,
    *,
    tau: str,
    evidence: str,
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


def _two(
    roots: Sequence[tp.RootObservation],
    cell: tp.CellSpec,
    subset: str,
    memo: str,
) -> list[dict[str, Any]]:
    return [
        finish(row)
        for row in tp.two_outcome_rows(roots, cell, subset=subset, memo=memo)
    ]


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------


def cells_rows(
    records: Sequence[RootRecord], cache: ObservationCache
) -> list[dict[str, Any]]:
    """The ``cells`` grid (E0, dominant, full lexicon)."""
    observed = cache.get("E0")
    rows: list[dict[str, Any]] = []
    for target, member in TARGETS.items():
        for trigger in mo.TRIGGERS:
            for subset in ("all", "core", "with_nde"):
                roots = [
                    observed[r.root_id]
                    for r in select(records, member, SUBSETS[subset])
                ]
                for tau in TAUS:
                    cell = _cell(
                        "cells", target, trigger, tau=tau, evidence="E0"
                    )
                    rows += _two(roots, cell, subset, "")
    return rows


def _nested(
    roots_by_id: Mapping[str, tp.RootObservation],
    records: Sequence[RootRecord],
    spec: CellDef,
    cell: tp.CellSpec,
) -> list[dict[str, Any]]:
    chosen = select(records, TARGETS[spec.target], SUBSETS[spec.nested_subset])
    by_id = {r.root_id: r for r in chosen}
    specific = [
        tp.NestedRule(f"nested:{n.name}", n.rule, _on_roots(n.member, by_id))
        for n in spec.nested
    ]
    rows = tp.nested_rows(
        [roots_by_id[r.root_id] for r in chosen],
        cell,
        specific,
        elsewhere="palatalize",
        elsewhere_subset="nested:elsewhere",
    )
    memo = f"{spec.memo}; nested in {spec.nested_subset}"
    return [finish({**row, "memo": memo}) for row in rows]


def r17_rows(
    records: Sequence[RootRecord], cache: ObservationCache
) -> list[dict[str, Any]]:
    """The R17 comparison rows (see the module docstring)."""
    rows: list[dict[str, Any]] = []
    frames = {
        (basis, V): frame(records, V, basis)
        for basis in FREQ_BASES
        for V in V_VALUES
    }
    for spec in R17_CELLS:
        target = TARGETS[spec.target]
        for evidence in EVIDENCE_SETS:
            observed = cache.get(evidence)
            for subset in ("all", "core"):
                roots = [
                    observed[r.root_id]
                    for r in select(records, target, SUBSETS[subset])
                ]
                for tau in TAUS:
                    for coding in tp.CODINGS:
                        cell = _cell(
                            spec.table,
                            spec.target,
                            spec.trigger,
                            tau=tau,
                            evidence=evidence,
                            coding=coding,
                        )
                        rows += _two(roots, cell, subset, spec.memo)
            chosen = select(records, target, SUBSETS["all"])
            for basis in FREQ_BASES:
                observed_b = cache.get(evidence, basis)
                for V in V_VALUES:
                    inside = frames[(basis, V)]
                    roots = [
                        observed_b[r.root_id]
                        for r in chosen
                        if r.root_id in inside
                    ]
                    for tau in TAUS:
                        cell = _cell(
                            spec.table,
                            spec.target,
                            spec.trigger,
                            tau=tau,
                            evidence=evidence,
                            basis=basis,
                            V=V,
                        )
                        rows += _two(roots, cell, "all", spec.memo)
            if spec.nested:
                for tau in TAUS:
                    cell = _cell(
                        spec.table,
                        spec.target,
                        spec.trigger,
                        tau=tau,
                        evidence=evidence,
                    )
                    rows += _nested(observed, records, spec, cell)
        observed0 = cache.get("E0")
        for name in STRATA:
            subset = f"stratum:{name}"
            roots = [
                observed0[r.root_id]
                for r in select(records, target, SUBSETS[subset])
            ]
            for tau in TAUS:
                cell = _cell(
                    spec.table,
                    spec.target,
                    spec.trigger,
                    tau=tau,
                    evidence="E0",
                )
                rows += _two(roots, cell, subset, spec.memo)
    return rows


def nde_rows(
    records: Sequence[RootRecord], cache: ObservationCache
) -> list[dict[str, Any]]:
    """Structural NDE classes as separate rows (E0, dominant)."""
    observed = cache.get("E0")
    rows: list[dict[str, Any]] = []
    for target, member in NDE_TARGETS.items():
        roots = [observed[r.root_id] for r in records if member(r)]
        for trigger in mo.TRIGGERS:
            for tau in TAUS:
                cell = _cell("nde", target, trigger, tau=tau, evidence="E0")
                rows += _two(roots, cell, "all", "structural NDE class")
    return rows


# ---------------------------------------------------------------------------
# Simulated learners (R18)
# ---------------------------------------------------------------------------

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


def _on_roots(
    member: Member, by_id: Mapping[str, RootRecord]
) -> Callable[[tp.RootObservation], bool]:
    """A record predicate as a predicate on RootObservations; roots
    without a record (the rest of the frame) are never members."""

    def test(root: tp.RootObservation) -> bool:
        record = by_id.get(root.root_id)
        return record is not None and member(record)

    return test


def learner_pool(
    records: Sequence[RootRecord], basis: str
) -> list[tp.RootObservation]:
    """The R24 frame: every noun root with a positive frequency, without
    cells (the cells' statuses are computed once per cell,
    :func:`cell_statuses`)."""
    return [
        tp.RootObservation(
            record.root_id, "none", float(record.frequency(basis)), {}
        )
        for record in records
        if record.frequency(basis) > 0
    ]


def cell_statuses(
    observed: Mapping[str, tp.RootObservation],
    members: Sequence[RootRecord],
    cell: tp.CellSpec,
) -> tp.CellStatuses:
    """The members' statuses in *cell*; every other root of the frame
    has no evidence there."""
    return tp.cell_statuses([observed[r.root_id] for r in members], cell)


@dataclass(frozen=True, slots=True)
class LearnerConfig:
    """One learner cell: an R17 cell × evidence set × τ."""

    spec: CellDef
    subset: str
    evidence: str
    tau: str
    cell: tp.CellSpec
    evaluate: learners.RuleEval

    @property
    def name(self) -> str:
        return (
            f"{self.spec.table} {self.spec.target} {self.spec.trigger} "
            f"{self.evidence} {self.tau}"
        )


def learner_configs(
    records: Sequence[RootRecord],
    cache: ObservationCache,
    *,
    taus: Sequence[str] = LEARNER_TAUS,
    evidence_sets: Sequence[str] = LEARNER_EVIDENCE,
    basis: str = "lemma",
) -> list[LearnerConfig]:
    """Every R17 cell's rule evaluation (subset ``all``; the R17a nested
    set uses ``with_nde``), from precomputed statuses
    (:func:`lexpipe.learners.status_verdict`)."""
    out: list[LearnerConfig] = []
    for spec in R17_CELLS:
        subset = spec.nested_subset if spec.nested else "all"
        chosen = select(records, TARGETS[spec.target], SUBSETS[subset])
        specific = [
            (n.rule, frozenset(r.root_id for r in chosen if n.member(r)))
            for n in spec.nested
        ]
        for evidence in evidence_sets:
            observed = cache.get(evidence, basis)
            for tau in taus:
                cell = _cell(
                    spec.table,
                    spec.target,
                    spec.trigger,
                    tau=tau,
                    evidence=evidence,
                    basis=basis,
                )
                statuses = cell_statuses(observed, chosen, cell)
                out.append(
                    LearnerConfig(
                        spec,
                        subset,
                        evidence,
                        tau,
                        cell,
                        learners.status_verdict(
                            statuses, specific, elsewhere="palatalize"
                        ),
                    )
                )
    return out


def learner_rows(
    records: Sequence[RootRecord],
    cache: ObservationCache,
    *,
    seed: int,
    n_samples: int = learners.DEFAULT_SAMPLES,
    V_values: Sequence[int] = V_VALUES,
    taus: Sequence[str] = LEARNER_TAUS,
    evidence_sets: Sequence[str] = LEARNER_EVIDENCE,
    basis: str = "lemma",
    workers: int = 1,
    log: Callable[[str], None] = print,
) -> list[dict[str, Any]]:
    """Learner shares for every R17 cell (subset ``all``; the R17a
    nested set uses ``with_nde``).

    All cells are evaluated on the same simulated vocabularies, drawn
    once from the R24 frame (:func:`lexpipe.learners.simulate_all`);
    each cell's shares equal a separate
    :func:`lexpipe.learners.simulate_learners` call with the same seed.
    """
    configs = learner_configs(
        records, cache, taus=taus, evidence_sets=evidence_sets, basis=basis
    )
    pool = learner_pool(records, basis)
    log(
        f"    learners: {len(configs)} cells x {len(V_values)} V x "
        f"{n_samples} samples over {len(pool):,} frame roots, "
        f"{workers} worker(s)"
    )
    shares = learners.simulate_all(
        pool,
        {config.name: config.evaluate for config in configs},
        V_values,
        n_samples,
        seed=seed,
        workers=workers,
    )
    rows: list[dict[str, Any]] = []
    for config in configs:
        for item in shares[config.name]:
            share = item.shares
            rows.append(
                {
                    "lang": LANG,
                    "table": config.spec.table,
                    "target_class": config.spec.target,
                    "trigger": config.spec.trigger,
                    "subset": config.subset,
                    "tau": config.tau,
                    "evidence_set": config.evidence,
                    "coding": config.cell.coding,
                    "freq_basis": basis,
                    "V": item.V,
                    "n_samples": item.n_samples,
                    "seed": seed,
                    "nested": "|".join(n.name for n in config.spec.nested),
                    "share_palatalize": share["palatalize"],
                    "share_retain": share["retain"],
                    "share_neither": share["neither"],
                    "share_conditioned": share["conditioned"],
                    "source": list(SOURCE),
                    "licence_class": LICENCE_CLASS,
                }
            )
    return rows
