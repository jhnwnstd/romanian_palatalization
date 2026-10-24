"""Harmonized and model-fit TP rows of Italian (R16–R19, R22, R24, R26).

Every number comes from the core: :func:`lexpipe.tp.two_outcome_rows`
and :func:`lexpipe.tp.nested_rows` (both outcomes, nested rules in
Yang's order), :func:`lexpipe.tp.top_n` (the V frames) and
:func:`lexpipe.tp.tracked_row` (fragility). This module only chooses which
roots go into which cell (:mod:`.roots` builds the observations).

Tables (``table`` column)
-------------------------
``R17a`` / ``R17b`` / ``R17c``
    the fixed cross-language rule family (R17), the comparison rows
    (dorsal rows also at evidence set ``E2_all``: E2 with the
    borrowed-whole and inherited derivatives):
    ``R17a`` dorsal / _[+front V] over the -i and -e cells (trigger
    ``-i/-e``); ``R17b`` dorsal / _[+high,+front] (trigger ``-i``) with
    its non-trigger cell ``-e`` (feminine plurals: N = e expected for
    ``palatalize``); ``R17c`` coronal / _[+high,+front] over coronal -i
    plurals (Italian has no coronal rule: N = e by construction), its
    non-trigger -e cell and the pooled coronal row (not in R17).
    Dorsal rows: subsets ``all`` / ``core`` / ``nouns`` /
    ``adjectives`` × τ × E0/E1/E2 × the three doublet codings (R19);
    V ∈ {500, 1000, 2000, 5000} on the lemma and paradigm frequency
    bases (R24; subsets ``all`` and ``nouns``, E0/E1, dominant
    coding); strata (``stratum:<s>``, E0, R22); R17a also nested
    (``with_nde``: the NDE stems occhio/psiche retain, elsewhere
    palatalize). Coronal rows: E0 only (no coronal verb or derivative
    evidence is collected; E1 = E2 = E0 there).
``cells``
    the descriptive grid: c, g, dorsal, t, d, s, z, coronal × -i, -e,
    -i/-e × ``all`` / ``nouns`` / ``adjectives`` × τ (E0, dominant).
``nde``
    the structural classes as separate rows (R16): nde:occhio,
    nde:psiche and frontstem (lexical palatals), E0.

Italian model fits (never comparison rows, R17)
-----------------------------------------------
``it_stress``
    the height-indexed /i j/ trigger (``-i``) nested by stress
    (``velar_post_tonic``: Giavazzi; Shao et al.): per stress subset
    (``stress:post_tonic`` / ``stress:not_post_tonic`` /
    ``stress:unknown``) and nested sets over the stress-known roots:
    ``giavazzi`` (post-tonic → retain; elsewhere palatalize: every
    -ico root counted as a root of its own), ``suffix_once`` (the same
    rules with each R26 suffix — -ico, -iaco — counted once, as one
    unit whose status pools its roots' forms: R26 "counts once as a
    morpheme in learner evidence"), ``orth_ico_once`` (every root in
    unstressed -ico/-ica, the orthographic class, as one unit: the
    upper bound, since R26 detection is a lower bound),
    ``suffix_then_stress`` (-ico suffix → palatalize; post-tonic →
    retain; elsewhere palatalize) and ``reverse`` (not post-tonic →
    palatalize; elsewhere retain); targets dorsal, c and g; E0/E1.
    The stress subsets and the giavazzi, suffix_once and orth_ico_once
    sets are also counted at τ = ``dictionary_explicit`` (E0, dorsal):
    only plurals an editor typed, not Module:it-common defaults (R11;
    a sample biased toward the plurals the module gets wrong).
``it_suffix``
    suffix rows (R26): each suffix its own row (``suffix:-ico``,
    ``suffix:-iaco``, ``suffixoid:-logo/-fago/-urgo``), the plain
    roots (``root``), the orthographic -ico lookalikes including the
    carico type (``orth:-ico``), vowel + -ico (``orth:-Vico``, all and
    post-tonic: arcaico, laico, whose diphthongal IPA makes them
    paroxytone), ``deverbal_from_are`` and ``shares_1st_conj_verb``;
    E0.
``it_verbs``
    verb units (velar-stem verbs, :class:`.roots.VerbUnit`) by
    conjugation class (``conj:<class>``), the over-appliers
    (``overapplying``), all verbs, and all verbs with ambiguous forms
    unattested above the dictionary threshold (``all_unambiguous``).
``it_nounverb``
    the pooled noun + verb row (units: dorsal noun/adjective roots with
    their E0 cells and verb units with their own cells) and the
    category-indexed rule sets: ``category`` (-are verbs → retain; -isc
    verbs → retain; elsewhere palatalize) and ``category_stress`` (+
    post-tonic nouns → retain).
``logtable``
    R27 rows: -logo/-fago/-urgo × the shared person/thing classes
    (:data:`lexpipe.suffixoids.CLASSES`) and adjective (E0); learned
    suffixoids only.
``spellout``
    R20 rows, built by :mod:`.spellout`.

Memos name the rule a row evaluates (``dorsal/_[+front V]`` …;
``: non-trigger cell``; ``; nested in <subset>``). The columns are
:data:`lexpipe.tp.TP_LONG_TRACKED_FIELDS` (:func:`lexpipe.tp.tracked_row`):
``fragility`` (fewest item changes that flip ``tolerated_trigger``),
``source`` and ``licence_class`` (aggregate rows of Wiktionary-derived
data; τ rows use Leipzig counts only as thresholds, R9).
"""

from __future__ import annotations

import re
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Final

from lexpipe import learners, suffixoids, tp
from lexpipe.attest import TAUS
from lexpipe.langs.italian import etymology, features, keys, verbs
from lexpipe.langs.italian.roots import (
    CORONAL,
    E2_ALL,
    EVIDENCE_SETS,
    EXPLICIT_TAU,
    FREQ_BASES,
    FRONTSTEM,
    TRIGGERS,
    Records,
    RootRecord,
    VerbUnit,
    observation,
    unit_observation,
)

LANG: Final[str] = "it"
SOURCE: Final[str] = "wiktionary_en"
LICENCE_CLASS: Final[str] = "by-sa"
FIELDS: Final[tuple[str, ...]] = tp.TP_LONG_TRACKED_FIELDS
V_VALUES: Final[tuple[int, ...]] = learners.DEFAULT_V_VALUES

Member = Callable[[RootRecord], bool]


# ---------------------------------------------------------------------------
# Selections
# ---------------------------------------------------------------------------


def _dorsal_seg(*segs: str) -> Member:
    def member(record: RootRecord) -> bool:
        return record.dorsal and (not segs or record.seg in segs)

    return member


def _coronal_seg(*segs: str) -> Member:
    def member(record: RootRecord) -> bool:
        return record.target_class == CORONAL and (
            not segs or record.seg in segs
        )

    return member


def _class(name: str) -> Member:
    def member(record: RootRecord) -> bool:
        return record.target_class == name

    return member


TARGETS: Final[Mapping[str, Member]] = {
    "c": _dorsal_seg("c"),
    "g": _dorsal_seg("g"),
    "dorsal": _dorsal_seg(),
    "t": _coronal_seg("t"),
    "d": _coronal_seg("d"),
    "s": _coronal_seg("s"),
    "z": _coronal_seg("z"),
    "coronal": _coronal_seg(),
}
NDE_TARGETS: Final[Mapping[str, Member]] = {
    "nde:occhio": _class("nde:occhio"),
    "nde:psiche": _class("nde:psiche"),
    "frontstem": _class(FRONTSTEM),
}


@dataclass(frozen=True, slots=True)
class Subset:
    """Roots of a subset and the scope of their E0 cells."""

    member: Member
    scope: str


def _any(record: RootRecord) -> bool:
    return True


def _core(record: RootRecord) -> bool:
    return not (record.ico_suffix or record.iaco_suffix or record.suffixoid)


def _stratum(name: str) -> Member:
    def member(record: RootRecord) -> bool:
        return record.stratum == name

    return member


SUBSETS: Final[Mapping[str, Subset]] = {
    "all": Subset(_any, "all"),
    "core": Subset(_core, "all"),
    "nouns": Subset(lambda r: r.has_noun, "noun"),
    "adjectives": Subset(lambda r: r.has_adj, "adj"),
    **{
        f"stratum:{name}": Subset(_stratum(name), "all")
        for name in etymology.STRATA
    },
}


def finish(row: Mapping[str, Any]) -> dict[str, Any]:
    """A core TP long row with the tracked columns
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
    """A :class:`lexpipe.tp.CellSpec` of Italian (doublets need
    :data:`.features.DOUBLET_MIN_TOKENS` tokens to be dominant)."""
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
        doublet_min_tokens=features.DOUBLET_MIN_TOKENS,
    )


def two(
    roots: Sequence[tp.RootObservation],
    spec: tp.CellSpec,
    subset: str,
    memo: str,
) -> list[dict[str, Any]]:
    """Both outcome rows, finished."""
    return [
        finish(row)
        for row in tp.two_outcome_rows(roots, spec, subset=subset, memo=memo)
    ]


class Observations:
    """RootObservations per (root, evidence, scope, basis), built once."""

    def __init__(self, records: Records) -> None:
        self._records = records
        self._cache: dict[tuple[str, str, str, str], tp.RootObservation] = {}

    def get(
        self,
        record: RootRecord,
        evidence: str,
        scope: str = "all",
        basis: str = "lemma",
    ) -> tp.RootObservation:
        key = (record.root_id, evidence, scope, basis)
        found = self._cache.get(key)
        if found is None:
            found = observation(record, evidence, scope=scope, basis=basis)
            self._cache[key] = found
        return found

    def select(
        self,
        target: Member,
        subset: str,
        evidence: str,
        *,
        basis: str = "lemma",
        inside: frozenset[str] | None = None,
    ) -> list[tp.RootObservation]:
        """Observations of a target class within a subset (and a frame)."""
        chosen = SUBSETS[subset]
        return [
            self.get(record, evidence, chosen.scope, basis)
            for record in self._records.roots.values()
            if target(record)
            and chosen.member(record)
            and (inside is None or record.root_id in inside)
        ]


def frame(records: Records, V: int, basis: str) -> frozenset[str]:
    """root_ids of the V most frequent noun roots of any final segment
    (R24; :func:`lexpipe.tp.top_n`, ties by root_id)."""
    rows = [
        {"lemma": r.root_id, "freq": r.frequency(basis, "noun")}
        for r in records.roots.values()
        if r.has_noun
    ]
    return frozenset(str(row["lemma"]) for row in tp.top_n(rows, V, "freq"))


# ---------------------------------------------------------------------------
# R17
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class R17Spec:
    """One R17 cell: target × trigger with its memo."""

    table: str
    target: str
    trigger: str
    memo: str
    evidence_sets: tuple[str, ...]
    subsets: tuple[str, ...]


DORSAL_SUBSETS: Final[tuple[str, ...]] = ("all", "core", "nouns", "adjectives")
CORONAL_SUBSETS: Final[tuple[str, ...]] = ("all", "nouns", "adjectives")
R17_SPECS: Final[tuple[R17Spec, ...]] = (
    R17Spec(
        "R17a",
        "dorsal",
        "-i/-e",
        "dorsal/_[+front V]",
        EVIDENCE_SETS,
        DORSAL_SUBSETS,
    ),
    R17Spec(
        "R17b",
        "dorsal",
        "-i",
        "dorsal/_[+high,+front]",
        EVIDENCE_SETS,
        DORSAL_SUBSETS,
    ),
    R17Spec(
        "R17b",
        "dorsal",
        "-e",
        "dorsal/_[+high,+front]: non-trigger cell",
        EVIDENCE_SETS,
        DORSAL_SUBSETS,
    ),
    R17Spec(
        "R17c",
        "coronal",
        "-i",
        "coronal/_[+high,+front]",
        ("E0",),
        CORONAL_SUBSETS,
    ),
    R17Spec(
        "R17c",
        "coronal",
        "-e",
        "coronal/_[+high,+front]: non-trigger cell",
        ("E0",),
        CORONAL_SUBSETS,
    ),
    R17Spec(
        "R17c",
        "coronal",
        "-i/-e",
        "coronal/_[+front V] (pooled; not in R17)",
        ("E0",),
        CORONAL_SUBSETS,
    ),
)
V_SUBSETS: Final[tuple[str, ...]] = ("all", "nouns")


def _nde_or_dorsal(record: RootRecord) -> bool:
    return record.dorsal or record.target_class.startswith("nde:")


def _on_records(
    member: Member, by_id: Mapping[str, RootRecord]
) -> Callable[[tp.RootObservation], bool]:
    """A record predicate as a predicate on RootObservations."""

    def test(root: tp.RootObservation) -> bool:
        record = by_id.get(root.root_id)
        return record is not None and member(record)

    return test


def _r17_basic(spec: R17Spec, obs: Observations) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    target = TARGETS[spec.target]
    for evidence in spec.evidence_sets:
        for subset in spec.subsets:
            roots = obs.select(target, subset, evidence)
            for tau in TAUS:
                for coding in tp.CODINGS:
                    spec_cell = cell(
                        spec.table,
                        spec.target,
                        spec.trigger,
                        tau=tau,
                        evidence=evidence,
                        coding=coding,
                    )
                    rows += two(roots, spec_cell, subset, spec.memo)
    return rows


def _r17_sweep(
    spec: R17Spec,
    obs: Observations,
    frames: Mapping[tuple[str, int], frozenset[str]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    target = TARGETS[spec.target]
    for evidence in spec.evidence_sets[:2]:
        for basis in FREQ_BASES:
            for V in V_VALUES:
                for subset in V_SUBSETS:
                    roots = obs.select(
                        target,
                        subset,
                        evidence,
                        basis=basis,
                        inside=frames[(basis, V)],
                    )
                    for tau in TAUS:
                        spec_cell = cell(
                            spec.table,
                            spec.target,
                            spec.trigger,
                            tau=tau,
                            evidence=evidence,
                            basis=basis,
                            V=V,
                        )
                        rows += two(roots, spec_cell, subset, spec.memo)
    return rows


def _r17_strata(spec: R17Spec, obs: Observations) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    target = TARGETS[spec.target]
    for name in etymology.STRATA:
        subset = f"stratum:{name}"
        roots = obs.select(target, subset, "E0")
        for tau in TAUS:
            spec_cell = cell(spec.table, spec.target, spec.trigger, tau=tau)
            rows += two(roots, spec_cell, subset, spec.memo)
    return rows


def _r17a_nested(records: Records, obs: Observations) -> list[dict[str, Any]]:
    chosen = [r for r in records.roots.values() if _nde_or_dorsal(r)]
    by_id = {r.root_id: r for r in chosen}
    specific = [
        tp.NestedRule(
            "nested:velar_front_lemma",
            "retain",
            _on_records(lambda r: r.target_class.startswith("nde:"), by_id),
        )
    ]
    rows: list[dict[str, Any]] = []
    for evidence in EVIDENCE_SETS:
        roots = [obs.get(r, evidence) for r in chosen]
        for tau in TAUS:
            spec_cell = cell(
                "R17a", "dorsal", "-i/-e", tau=tau, evidence=evidence
            )
            nested = tp.nested_rows(
                roots,
                spec_cell,
                specific,
                elsewhere="palatalize",
                elsewhere_subset="nested:elsewhere",
            )
            memo = "dorsal/_[+front V]; nested in with_nde"
            rows += [finish({**row, "memo": memo}) for row in nested]
    return rows


def _r17_e2_all(spec: R17Spec, obs: Observations) -> list[dict[str, Any]]:
    """E2 with the borrowed-whole and inherited derivatives (evidence set
    ``E2_all``, subset all, dominant coding): the sensitivity of the E2
    rows to the stratum filter (design §5)."""
    if "E2" not in spec.evidence_sets:
        return []
    rows: list[dict[str, Any]] = []
    roots = obs.select(TARGETS[spec.target], "all", E2_ALL)
    for tau in TAUS:
        spec_cell = cell(
            spec.table, spec.target, spec.trigger, tau=tau, evidence=E2_ALL
        )
        rows += two(roots, spec_cell, "all", spec.memo)
    return rows


def r17_rows(records: Records, obs: Observations) -> list[dict[str, Any]]:
    """The R17 comparison rows (module docstring)."""
    frames = {
        (basis, V): frame(records, V, basis)
        for basis in FREQ_BASES
        for V in V_VALUES
    }
    rows: list[dict[str, Any]] = []
    for spec in R17_SPECS:
        rows += _r17_basic(spec, obs)
        rows += _r17_sweep(spec, obs, frames)
        rows += _r17_strata(spec, obs)
        rows += _r17_e2_all(spec, obs)
    rows += _r17a_nested(records, obs)
    return rows


# ---------------------------------------------------------------------------
# Descriptive grid and NDE rows
# ---------------------------------------------------------------------------


def cells_rows(obs: Observations) -> list[dict[str, Any]]:
    """The ``cells`` grid (E0, dominant)."""
    rows: list[dict[str, Any]] = []
    for target, member in TARGETS.items():
        for trigger in TRIGGERS:
            for subset in ("all", "nouns", "adjectives"):
                roots = obs.select(member, subset, "E0")
                for tau in TAUS:
                    spec_cell = cell("cells", target, trigger, tau=tau)
                    rows += two(roots, spec_cell, subset, "")
    return rows


def nde_rows(obs: Observations) -> list[dict[str, Any]]:
    """Structural classes as separate rows (E0, dominant)."""
    rows: list[dict[str, Any]] = []
    for target, member in NDE_TARGETS.items():
        roots = obs.select(member, "all", "E0")
        for trigger in TRIGGERS:
            for tau in TAUS:
                spec_cell = cell("nde", target, trigger, tau=tau)
                rows += two(roots, spec_cell, "all", "structural class")
    return rows


# ---------------------------------------------------------------------------
# Model fits: stress, suffixes, verbs, noun + verb
# ---------------------------------------------------------------------------


def _post(value: bool | None) -> Member:
    def member(record: RootRecord) -> bool:
        return record.velar_post_tonic is value

    return member


STRESS_SUBSETS: Final[Mapping[str, Member]] = {
    "stress:post_tonic": _post(True),
    "stress:not_post_tonic": _post(False),
    "stress:unknown": _post(None),
}


def r26_unit(record: RootRecord) -> str:
    """The R26 suffix unit of a root ("" for a plain root): each suffix
    counts once (``suffix:-ico``, ``suffix:-iaco``)."""
    if record.ico_suffix:
        return "suffix:-ico"
    if record.iaco_suffix:
        return "suffix:-iaco"
    return ""


def orth_ico_unit(record: RootRecord) -> str:
    """Every root in unstressed -ico/-ica as one unit (``orth:-ico``):
    the upper bound of R26's counts-once rule."""
    return "orth:-ico" if record.orth_ico else ""


UnitOf = Callable[[RootRecord], str]


@dataclass(frozen=True, slots=True)
class NestedSet:
    """A nested rule set over records (most specific rule first);
    *units* collapses suffix roots into one unit each (R26)."""

    name: str
    rules: tuple[tuple[str, str, Member], ...]  # (subset, rule, member)
    elsewhere: str
    units: UnitOf | None = None


def pooled_unit(
    unit_id: str, members: Sequence[tp.RootObservation]
) -> tp.RootObservation:
    """One unit standing for its member roots (R26 "counts once"): its
    realizations in every trigger are the members' (their paradigm
    cells kept), so its status pools their forms — a doublet coded by
    the tokens of all of them."""
    triggers: dict[str, list[tp.Observation]] = {}
    for member in members:
        for trigger, found in member.triggers.items():
            triggers.setdefault(trigger, []).extend(found)
    return tp.RootObservation(
        unit_id,
        DORSAL_UNIT,
        float(sum(member.frequency for member in members)),
        triggers,
    )


DORSAL_UNIT: Final[str] = "dorsal"


def collapse(
    chosen: Sequence[RootRecord],
    roots: Sequence[tp.RootObservation],
    units: UnitOf,
) -> tuple[list[RootRecord], list[tp.RootObservation]]:
    """*chosen* (aligned with *roots*) with every unit's members replaced
    by one pooled unit (:func:`pooled_unit`), units after the plain
    roots in unit order."""
    kept: list[RootRecord] = []
    plain: list[tp.RootObservation] = []
    grouped: dict[str, list[tp.RootObservation]] = {}
    for record, root in zip(chosen, roots, strict=True):
        unit = units(record)
        if unit:
            grouped.setdefault(unit, []).append(root)
        else:
            kept.append(record)
            plain.append(root)
    pooled = [pooled_unit(unit, grouped[unit]) for unit in sorted(grouped)]
    return kept, plain + pooled


NESTED_STRESS: Final[tuple[NestedSet, ...]] = (
    NestedSet(
        "giavazzi",
        (("nested:post_tonic", "retain", _post(True)),),
        "palatalize",
    ),
    NestedSet(
        "suffix_once",
        (("nested:post_tonic", "retain", _post(True)),),
        "palatalize",
        r26_unit,
    ),
    NestedSet(
        "orth_ico_once",
        (("nested:post_tonic", "retain", _post(True)),),
        "palatalize",
        orth_ico_unit,
    ),
    NestedSet(
        "suffix_then_stress",
        (
            ("nested:ico_suffix", "palatalize", lambda r: r.ico_suffix),
            ("nested:post_tonic", "retain", _post(True)),
        ),
        "palatalize",
    ),
    NestedSet(
        "reverse",
        (("nested:not_post_tonic", "palatalize", _post(False)),),
        "retain",
    ),
)


def nested_set_rows(
    chosen: Sequence[RootRecord],
    roots: Sequence[tp.RootObservation],
    spec_cell: tp.CellSpec,
    nested: NestedSet,
    memo: str,
) -> list[dict[str, Any]]:
    """Rows of one nested set over *roots* (aligned with *chosen*); a
    set with units counts each unit once (:func:`collapse`: a unit is
    in no specific rule's subclass, so it falls to the elsewhere
    rule)."""
    if nested.units is not None:
        chosen, roots = collapse(chosen, roots, nested.units)
    by_id = {r.root_id: r for r in chosen}
    specific = [
        tp.NestedRule(subset, rule, _on_records(member, by_id))
        for subset, rule, member in nested.rules
    ]
    rows = tp.nested_rows(
        roots,
        spec_cell,
        specific,
        elsewhere=nested.elsewhere,
        elsewhere_subset="nested:elsewhere",
    )
    note = f"{memo}; nested set {nested.name}"
    return [finish({**row, "memo": note}) for row in rows]


EXPLICIT_SETS: Final[frozenset[str]] = frozenset(
    {"giavazzi", "suffix_once", "orth_ico_once"}
)


def _explicit_rows(
    records: Records, obs: Observations
) -> list[dict[str, Any]]:
    """``it_stress`` dorsal rows at τ = :data:`.roots.EXPLICIT_TAU`
    (E0): the stress subsets and the giavazzi, suffix_once and
    orth_ico_once sets over editor-typed plurals only."""
    rows: list[dict[str, Any]] = []
    memo = "height-indexed /i j/ trigger; explicit plurals only"
    member = TARGETS["dorsal"]
    spec_cell = cell("it_stress", "dorsal", "-i", tau=EXPLICIT_TAU)
    for name, stress in STRESS_SUBSETS.items():
        chosen = [r for r in records.roots.values() if member(r) and stress(r)]
        rows += two([obs.get(r, "E0") for r in chosen], spec_cell, name, memo)
    known = [
        r
        for r in records.roots.values()
        if member(r) and r.velar_post_tonic is not None
    ]
    roots = [obs.get(r, "E0") for r in known]
    for nested in NESTED_STRESS:
        if nested.name in EXPLICIT_SETS:
            rows += nested_set_rows(
                known, roots, spec_cell, nested, f"{memo}; stress-known roots"
            )
    return rows


def stress_rows(records: Records, obs: Observations) -> list[dict[str, Any]]:
    """``it_stress`` (module docstring)."""
    rows: list[dict[str, Any]] = []
    memo = "height-indexed /i j/ trigger"
    for target in ("dorsal", "c", "g"):
        member = TARGETS[target]
        for evidence in ("E0", "E1"):
            for name, stress in STRESS_SUBSETS.items():
                chosen = [
                    r
                    for r in records.roots.values()
                    if member(r) and stress(r)
                ]
                roots = [obs.get(r, evidence) for r in chosen]
                for tau in TAUS:
                    spec_cell = cell(
                        "it_stress", target, "-i", tau=tau, evidence=evidence
                    )
                    rows += two(roots, spec_cell, name, memo)
            known = [
                r
                for r in records.roots.values()
                if member(r) and r.velar_post_tonic is not None
            ]
            roots = [obs.get(r, evidence) for r in known]
            for nested in NESTED_STRESS:
                for tau in TAUS:
                    spec_cell = cell(
                        "it_stress", target, "-i", tau=tau, evidence=evidence
                    )
                    rows += nested_set_rows(
                        known,
                        roots,
                        spec_cell,
                        nested,
                        f"{memo}; stress-known roots",
                    )
    rows += _explicit_rows(records, obs)
    return rows


def _suffixoid(name: str) -> Member:
    def member(record: RootRecord) -> bool:
        return record.suffixoid == name

    return member


_VICO_RE: Final[re.Pattern[str]] = re.compile(r"[aeiou]ic[oa]\Z")


def _vico(record: RootRecord) -> bool:
    """Vowel + -ico/-ica (arcaico, eroico, laico): the IPA often writes
    the vowel sequence as a diphthong, so the stress counts as
    penultimate (it_phon.json open question on laico)."""
    return bool(_VICO_RE.search(keys.match_key(record.lemma)))


SUFFIX_SUBSETS: Final[Mapping[str, Member]] = {
    "suffix:-ico": lambda r: r.ico_suffix,
    "suffix:-iaco": lambda r: r.iaco_suffix,
    "suffixoid:-logo": _suffixoid("-logo"),
    "suffixoid:-fago": _suffixoid("-fago"),
    "suffixoid:-urgo": _suffixoid("-urgo"),
    "root": _core,
    "orth:-ico": lambda r: r.orth_ico,
    "orth:-Vico": _vico,
    "orth:-Vico post-tonic": lambda r: _vico(r) and r.velar_post_tonic is True,
    "deverbal_from_are": lambda r: r.deverbal_from_are,
    "shares_1st_conj_verb": lambda r: r.shares_1st_conj,
}


def suffix_rows(records: Records, obs: Observations) -> list[dict[str, Any]]:
    """``it_suffix`` (module docstring)."""
    rows: list[dict[str, Any]] = []
    for name, member in SUFFIX_SUBSETS.items():
        chosen = [r for r in records.roots.values() if r.dorsal and member(r)]
        roots = [obs.get(r, "E0") for r in chosen]
        for tau in TAUS:
            spec_cell = cell("it_suffix", "dorsal", "-i", tau=tau)
            rows += two(roots, spec_cell, name, "suffix rows (R26)")
    return rows


VERB_SUBSETS: Final[tuple[str, ...]] = (
    *(f"conj:{name}" for name in verbs.CONJ_CLASSES),
    "overapplying",
    "all",
    "all_unambiguous",
)


def _verb_member(subset: str) -> Callable[[VerbUnit], bool]:
    if subset.startswith("conj:"):
        name = subset.split(":", 1)[1]
        return lambda unit: unit.conj_class == name
    if subset == "overapplying":
        return lambda unit: unit.paradigm_type == verbs.OVERAPPLYING
    return lambda unit: True


def verb_rows(records: Records) -> list[dict[str, Any]]:
    """``it_verbs`` (module docstring)."""
    rows: list[dict[str, Any]] = []
    for subset in VERB_SUBSETS:
        member = _verb_member(subset)
        unambiguous = subset == "all_unambiguous"
        roots = [
            unit_observation(unit, unambiguous=unambiguous)
            for unit in records.units.values()
            if member(unit)
        ]
        for trigger in TRIGGERS:
            for tau in TAUS:
                spec_cell = cell(
                    "it_verbs", "dorsal", trigger, tau=tau, evidence="E1"
                )
                rows += two(roots, spec_cell, subset, "verb units")
    return rows


@dataclass(frozen=True, slots=True)
class Unit:
    """A unit of the pooled noun + verb table."""

    root_id: str
    category: str  # noun | verb:<conj class>
    post_tonic: bool | None


def _category(prefix: str) -> Callable[[Unit], bool]:
    def member(unit: Unit) -> bool:
        return unit.category.startswith(prefix)

    return member


def _noun_post_tonic(unit: Unit) -> bool:
    return unit.category == "noun" and unit.post_tonic is True


NOUN_VERB_SETS: Final[
    tuple[tuple[str, tuple[tuple[str, str, Callable[[Unit], bool]], ...]], ...]
] = (
    (
        "category",
        (
            ("nested:verb_are", "retain", _category("verb:are")),
            ("nested:verb_isc", "retain", _category("verb:isc")),
        ),
    ),
    (
        "category_stress",
        (
            ("nested:verb_are", "retain", _category("verb:are")),
            ("nested:verb_isc", "retain", _category("verb:isc")),
            ("nested:noun_post_tonic", "retain", _noun_post_tonic),
        ),
    ),
)


def nounverb_units(
    records: Records, obs: Observations
) -> tuple[list[tp.RootObservation], dict[str, Unit]]:
    """Dorsal noun/adjective roots (E0) and verb units (own cells)."""
    roots: list[tp.RootObservation] = []
    units: dict[str, Unit] = {}
    for record in records.roots.values():
        if record.dorsal:
            roots.append(obs.get(record, "E0"))
            units[record.root_id] = Unit(
                record.root_id, "noun", record.velar_post_tonic
            )
    for verb_unit in records.units.values():
        found = unit_observation(verb_unit)
        roots.append(found)
        units[found.root_id] = Unit(
            found.root_id, f"verb:{verb_unit.conj_class}", None
        )
    return roots, units


def nounverb_rows(records: Records, obs: Observations) -> list[dict[str, Any]]:
    """``it_nounverb`` (module docstring)."""
    roots, units = nounverb_units(records, obs)
    rows: list[dict[str, Any]] = []
    memo = "pooled noun + verb units"
    for trigger in ("-i", "-i/-e"):
        for tau in TAUS:
            spec_cell = cell(
                "it_nounverb", "dorsal", trigger, tau=tau, evidence="E1"
            )
            rows += two(roots, spec_cell, "pooled", memo)
            for name, rules in NOUN_VERB_SETS:
                specific = [
                    tp.NestedRule(
                        subset,
                        rule,
                        _on_units(member, units),
                    )
                    for subset, rule, member in rules
                ]
                nested = tp.nested_rows(
                    roots,
                    spec_cell,
                    specific,
                    elsewhere="palatalize",
                    elsewhere_subset="nested:elsewhere",
                )
                note = f"{memo}; nested set {name}"
                rows += [finish({**row, "memo": note}) for row in nested]
    return rows


def _on_units(
    member: Callable[[Unit], bool], units: Mapping[str, Unit]
) -> Callable[[tp.RootObservation], bool]:
    def test(root: tp.RootObservation) -> bool:
        unit = units.get(root.root_id)
        return unit is not None and member(unit)

    return test


# ---------------------------------------------------------------------------
# -logo/-fago/-urgo (R27)
# ---------------------------------------------------------------------------

# The shared classes (lexpipe.suffixoids) plus Italian adjectives.
SEMANTIC_CLASSES: Final[tuple[str, ...]] = (
    *suffixoids.CLASSES,
    "adjective",
    "any",
)


def logtable_rows(records: Records, obs: Observations) -> list[dict[str, Any]]:
    """``logtable`` rows: suffixoid × semantic class (E0); roots whose
    suffixoid is learned (lookalikes have none, :mod:`.features`)."""
    rows: list[dict[str, Any]] = []
    for suffixoid in etymology.SUFFIXOIDS:
        for klass in SEMANTIC_CLASSES:
            chosen = [
                r
                for r in records.roots.values()
                if r.dorsal
                and r.suffixoid == suffixoid
                and (klass == "any" or r.semantic_class == klass)
            ]
            roots = [obs.get(r, "E0") for r in chosen]
            for tau in TAUS:
                spec_cell = cell("logtable", "dorsal", "-i", tau=tau)
                rows += two(
                    roots, spec_cell, f"{suffixoid}:{klass}", "R27 suffixoid"
                )
    return rows


def all_rows(
    records: Records,
    extra: Iterable[Mapping[str, Any]] = (),
) -> list[dict[str, Any]]:
    """Every table of the module docstring, in table order."""
    obs = Observations(records)
    rows = r17_rows(records, obs)
    rows += cells_rows(obs)
    rows += nde_rows(obs)
    rows += stress_rows(records, obs)
    rows += suffix_rows(records, obs)
    rows += verb_rows(records)
    rows += nounverb_rows(records, obs)
    rows += logtable_rows(records, obs)
    rows += [dict(row) for row in extra]
    return rows
