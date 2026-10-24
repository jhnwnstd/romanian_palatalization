"""Cells of the three Czech TP analyses and their observations (design §6).

Everything here turns a :class:`.roots.RootRecord` into the
:class:`lexpipe.tp.Observation` records the core counts (R16–R19), so
the counting itself is :mod:`lexpipe.tp`'s, identical across languages.
An observation's ``allomorph`` is ``P`` when the realization stands
before the trigger (``alternates`` = it shows the reflex) and ``N``
when it uses a non-palatalizing allomorph instead (it then *shields*
the root, R16).

(a) Rule TP (:data:`RULE_TRIGGERS`)
-----------------------------------
One trigger per §6 P cell: NOM.PL m.an -i, LOC.SG m.in / n -e,
DAT/LOC.SG f -e, LOC.PL m / n -ích, VOC.SG m.an -e, and the hard
adjective NOM.PL.M.AN -í. A P realization is an observation with
``alternates = mutation``; a competing N allomorph (-ové/-é, -u, -ách)
is an ``N`` observation. So N_trigger = roots with an attested P form,
e = those without the change (§6 (a)). The non-trigger rows of R16
(NOM.PL -é, INS.SG -em, LOC.PL -ech, VOC.SG -e of coronals and of r
after a vowel) put every such realization in as ``P`` with
``alternates = False``: N = e by construction.

(b) Competition TP (:data:`COMPETITIONS`)
-----------------------------------------
Per cell, slot A is the palatalizing (or, for non-target stems, the
front default) allomorph and slot B its competitor::

    LOC.SG m.in, n     A -e/-ě    B -u
    NOM.PL m.an        A -i/-í    B -ové, -é   (and A vs -ové alone)
    LOC.PL m.an/m.in/n A -ích     B -ách, -ech   velars
                       A -ech     B -ách         other hard stems
    VOC.SG m.an        A -e       B -u
    ADV (adjectives)   A -e/-ě    B -y, -o
    ADV -y             A -e/-ě    B -y   (the -o adverbs left out)
    ADV -o             A -e/-ě    B -o   (blízko, daleko, těžko: a
                                          formation of its own)

Animacy-specific cells (NOM.PL, VOC.SG m.an; LOC.SG, LOC.PL m.in) take
only rows whose animacy is known (gender class ``m.an``/``m.in``): a
cs-only masculine without animacy (``m``) enters the pooled masculine
cells (LOC.PL m, INS.SG) only. The NOM.PL competitions leave out the
animate masculines in -a (předseda, kolega, monarcha:
:attr:`RootRecord.a_masculine`), whose -ové/-é is their paradigm's.

Every attested realization in A or B is a ``P`` observation whose
``alternates`` means "chose A" (the encoding the Romanian spell-out
rows use), so the core's ``palatalize`` row is the P-default rule
(e = B choosers) and ``retain`` the N-default rule (e = A choosers);
roots with both are doublets coded by token share (R19). The labial,
l, s and z hard stems are the baseline (§6): their A slot involves no
palatalization. Only hard-declension stems count (soft učitel takes
-é and soft endings: :attr:`RootRecord.in_baseline`).

(c) Spell-out leave-one-cell-out
--------------------------------
:func:`alternation_evidence` says whether a root shows the reflex at
τ in its evidence set without the tested cell (and its copies:
NOM.PL with VOC.PL). E0 is nominal inflection (all cells of the root's
noun and adjective records), E1 adds -it/-et verbs (and their
participles, which are forms of the verb), E2 every derivative with a
palatalized stem. Only P realizations with the change count (the R16
non-trigger rows never do). Coronal t/d/n are excluded from (a)'s
learner evidence and (c) because their change is inferred (§6).

R17 rule family (:data:`R17_DEFS`)
----------------------------------
dorsal/_[+front V], dorsal/_[+high,+front] and coronal/_[+high,+front]
over **every** realization whose ending starts with a vowel of the
class, whatever its P/N status: front_mid e/ě/é, front_high i/í, and
y/ý ([ɪ iː]) as a separate sensitivity. E1/E2 add derivatives whose
suffix starts with such a vowel (straš-it, kluč-ina, tank-ista);
E2 also gives the non-trigger rows for consonant-initial
(kniž-ní, kníž-ka) and back-vowel (Praž-an) suffixes, which
palatalize although the rule family predicts no change. Nouns of the
adjectival declension (hajný, vysoká: GEN.SG -ého, -é) are left out of
the noun rows: their -é, -ém, -ého are the adjective's endings, counted
in the adjective rows.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Final

from lexpipe import tp
from lexpipe.langs.czech import ortho
from lexpipe.langs.czech.roots import (
    ADJ,
    BACK,
    CONSONANT,
    FRONT_HIGH,
    FRONT_MID,
    NOUN,
    Y_VOWEL,
    DerivativeEvidence,
    Realization,
    RootRecord,
)

EVIDENCE_SETS: Final[tuple[str, ...]] = ("E0", "E1", "E2")
INFERRED_FRONT: Final[frozenset[str]] = frozenset({"ě", "i", "í"})
VERB_SUFFIXES: Final[frozenset[str]] = frozenset({"-it", "-et"})
MASC: Final[frozenset[str]] = frozenset({"m.an", "m.in", "m"})
# Animacy-specific cells take rows of known animacy only (module
# docstring): ``m`` (a cs-only masculine without animacy) is in MASC.
MASC_AN: Final[frozenset[str]] = frozenset({"m.an"})
MASC_IN: Final[frozenset[str]] = frozenset({"m.in"})
NEUT: Final[frozenset[str]] = frozenset({"n"})
FEM: Final[frozenset[str]] = frozenset({"f"})
ADJ_GENDER: Final[frozenset[str]] = frozenset({ADJ})


# ---------------------------------------------------------------------------
# Target classes and subsets
# ---------------------------------------------------------------------------

Member = Callable[[RootRecord], bool]


def _klass(*names: str) -> Member:
    wanted = frozenset(names)
    return lambda record: record.target_class in wanted


def _segment(name: str) -> Member:
    def member(record: RootRecord) -> bool:
        if record.target_class != ortho.VELAR:
            return False
        if name in ("sk", "ck"):
            return record.cluster == name
        return record.stem_final == name and record.cluster not in (
            "sk",
            "ck",
        )

    return member


def _hard(*classes: str) -> Member:
    """Hard-declension stems of *classes* (:attr:`RootRecord.in_baseline`)."""

    def member(record: RootRecord) -> bool:
        return record.in_baseline and record.target_class in classes

    return member


def _r_cluster(after_consonant: bool) -> Member:
    def member(record: RootRecord) -> bool:
        return record.target_class == ortho.R_CLASS and (
            (record.cluster == "Cr") == after_consonant
        )

    return member


BASELINE: Final[tuple[str, ...]] = (
    ortho.LABIAL,
    ortho.LATERAL,
    ortho.SIBILANT,
)
RULE_TARGETS: Final[Mapping[str, Member]] = {
    "velar": _klass(ortho.VELAR),
    "velar:k": _segment("k"),
    "velar:g": _segment("g"),
    "velar:h": _segment("h"),
    "velar:ch": _segment("ch"),
    "velar:sk": _segment("sk"),
    "velar:ck": _segment("ck"),
    "r": _klass(ortho.R_CLASS),
    "r:Cr": _r_cluster(True),
    "r:Vr": _r_cluster(False),
    "c": _klass(ortho.C_CLASS),
    "coronal_inferred": _klass(ortho.CORONAL),
}
# Classes that get the subset and stratum rows (segments get "all").
MAJOR_TARGETS: Final[tuple[str, ...]] = (
    "velar",
    "r",
    "c",
    "coronal_inferred",
)
COMPETITION_TARGETS: Final[Mapping[str, Member]] = {
    "velar": _klass(ortho.VELAR),
    "r": _klass(ortho.R_CLASS),
    "coronal_inferred": _klass(ortho.CORONAL),
    "baseline": _hard(*BASELINE),
    "labial": _hard(ortho.LABIAL),
    "lateral": _hard(ortho.LATERAL),
    "sibilant": _hard(ortho.SIBILANT),
}
STRATA: Final[tuple[str, ...]] = (
    "inherited",
    "borrowed whole",
    "internal formation",
    "unknown",
)


# c changes only under the first palatalization (VOC otec → otče);
# the second-palatalization cells leave it as it is (otci).
TARGET_TRIGGERS: Final[Mapping[str, frozenset[str]]] = {
    "c": frozenset({"VOC.SG m.an -e"}),
}


def rule_applies(target: str, trigger: "RuleTrigger") -> bool:
    """Is *target* evaluated in the (a) row of *trigger*?"""
    allowed = TARGET_TRIGGERS.get(target)
    return allowed is None or trigger.name in allowed


def subsets() -> dict[str, Member]:
    """``all``, ``core`` (R26 suffix targets out), ``suffix`` and one
    ``stratum:<s>`` subset per R22 stratum (Wiktionary)."""
    out: dict[str, Member] = {
        "all": lambda record: True,
        "core": lambda record: not record.suffix_target,
        "suffix": lambda record: record.suffix_target,
    }
    for name in STRATA:
        out[f"stratum:{name}"] = (lambda s: lambda r: r.stratum == s)(name)
    return out


# ---------------------------------------------------------------------------
# (a) Rule triggers
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class RuleTrigger:
    """A §6 trigger cell (``nontrigger``: an R16 N = e row)."""

    name: str
    cells: frozenset[str]
    genders: frozenset[str]
    pos: str
    nontrigger: bool = False
    endings: frozenset[str] = frozenset()

    @property
    def key(self) -> str:
        return f"rule:{self.name}"


RULE_TRIGGERS: Final[tuple[RuleTrigger, ...]] = (
    RuleTrigger("NOM.PL m.an -i", frozenset({"NOM.PL"}), MASC_AN, NOUN),
    RuleTrigger("LOC.SG m.in -e", frozenset({"LOC.SG"}), MASC_IN, NOUN),
    RuleTrigger("LOC.SG n -e", frozenset({"LOC.SG"}), NEUT, NOUN),
    RuleTrigger("DAT/LOC.SG f -e", frozenset({"DAT.SG"}), FEM, NOUN),
    RuleTrigger("LOC.PL m -ích", frozenset({"LOC.PL"}), MASC, NOUN),
    RuleTrigger("LOC.PL n -ích", frozenset({"LOC.PL"}), NEUT, NOUN),
    RuleTrigger("VOC.SG m.an -e", frozenset({"VOC.SG"}), MASC_AN, NOUN),
    RuleTrigger("NOM.PL.M.AN -í", frozenset({"NOM.PL.M.AN"}), ADJ_GENDER, ADJ),
    RuleTrigger(
        "NOM.PL m.an -é",
        frozenset({"NOM.PL"}),
        MASC_AN,
        NOUN,
        True,
        frozenset({"é"}),
    ),
    RuleTrigger(
        "INS.SG m/n -em",
        frozenset({"INS.SG"}),
        MASC | NEUT,
        NOUN,
        True,
        frozenset({"em"}),
    ),
    RuleTrigger(
        "LOC.PL m/n -ech",
        frozenset({"LOC.PL"}),
        MASC | NEUT,
        NOUN,
        True,
        frozenset({"ech"}),
    ),
    RuleTrigger(
        "VOC.SG m.an -e (non-trigger)",
        frozenset({"VOC.SG"}),
        MASC_AN,
        NOUN,
        True,
        frozenset({"e", "ě"}),
    ),
)


def _observation(
    item: Realization, alternates: bool | None, allomorph: str
) -> tp.Observation:
    return tp.Observation(
        alternates=alternates,
        attested_by_tau=item.attested_by_tau(),
        allomorph=allomorph,
        form=item.form,
        count=item.count,
    )


def rule_observations(
    record: RootRecord, trigger: RuleTrigger
) -> tuple[tp.Observation, ...]:
    """(a) observations of one trigger (module docstring)."""
    out: list[tp.Observation] = []
    for item in record.realizations:
        if (
            item.cell not in trigger.cells
            or item.gender not in trigger.genders
        ):
            continue
        if trigger.nontrigger:
            # -é is a competing allomorph (n_kind allomorph) and a
            # front vowel that does not trigger: both readings hold.
            if (
                item.allomorph == ortho.N
                and item.n_kind in (ortho.N_NONTRIGGER, ortho.N_ALLOMORPH)
                and item.ending in trigger.endings
            ):
                out.append(_observation(item, False, "P"))
            continue
        if item.allomorph == ortho.P:
            out.append(_observation(item, item.mutation, "P"))
        elif item.allomorph == ortho.N and item.n_kind == ortho.N_ALLOMORPH:
            out.append(_observation(item, None, "N"))
    return tuple(out)


# ---------------------------------------------------------------------------
# (b) Competitions
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Competition:
    """A cell with a palatalizing slot A and its competitor B."""

    name: str
    cells: frozenset[str]
    genders: frozenset[str]
    pos: str
    a: frozenset[str]
    b: frozenset[str]
    velar_a: frozenset[str] = frozenset()
    velar_b: frozenset[str] = frozenset()
    label: str = ""

    @property
    def key(self) -> str:
        return f"comp:{self.name}"

    def slots(
        self, target_class: str
    ) -> tuple[frozenset[str], frozenset[str]]:
        if target_class == ortho.VELAR and self.velar_a:
            return self.velar_a, self.velar_b
        return self.a, self.b


_E: Final[frozenset[str]] = frozenset({"e", "ě"})
_U: Final[frozenset[str]] = frozenset({"u"})
_I: Final[frozenset[str]] = frozenset({"i", "í"})
_ICH: Final[frozenset[str]] = frozenset({"ích"})
_ACH: Final[frozenset[str]] = frozenset({"ách"})
_ECH: Final[frozenset[str]] = frozenset({"ech"})

COMPETITIONS: Final[tuple[Competition, ...]] = (
    Competition(
        "LOC.SG m.in",
        frozenset({"LOC.SG"}),
        MASC_IN,
        NOUN,
        _E,
        _U,
        label="-e/-ě vs -u",
    ),
    Competition(
        "LOC.SG n",
        frozenset({"LOC.SG"}),
        NEUT,
        NOUN,
        _E,
        _U,
        label="-e/-ě vs -u",
    ),
    Competition(
        "NOM.PL m.an",
        frozenset({"NOM.PL"}),
        MASC_AN,
        NOUN,
        _I,
        frozenset({"ové", "é"}),
        label="-i vs -ové/-é",
    ),
    Competition(
        "NOM.PL m.an -ové",
        frozenset({"NOM.PL"}),
        MASC_AN,
        NOUN,
        _I,
        frozenset({"ové"}),
        label="-i vs -ové",
    ),
    Competition(
        "LOC.PL m.an",
        frozenset({"LOC.PL"}),
        frozenset({"m.an"}),
        NOUN,
        _ECH,
        _ACH,
        _ICH,
        _ACH | _ECH,
        "-ích vs -ách/-ech (velars); -ech vs -ách",
    ),
    Competition(
        "LOC.PL m.in",
        frozenset({"LOC.PL"}),
        MASC_IN,
        NOUN,
        _ECH,
        _ACH,
        _ICH,
        _ACH | _ECH,
        "-ích vs -ách/-ech (velars); -ech vs -ách",
    ),
    Competition(
        "LOC.PL n",
        frozenset({"LOC.PL"}),
        NEUT,
        NOUN,
        _ECH,
        _ACH,
        _ICH,
        _ACH | _ECH,
        "-ích vs -ách/-ech (velars); -ech vs -ách",
    ),
    Competition(
        "VOC.SG m.an",
        frozenset({"VOC.SG"}),
        MASC_AN,
        NOUN,
        _E,
        _U,
        label="-e vs -u",
    ),
    Competition(
        "ADV",
        frozenset({"ADV"}),
        ADJ_GENDER,
        ADJ,
        _E,
        frozenset({"y", "o"}),
        label="adverb -e/-ě vs -y/-o",
    ),
    Competition(
        "ADV -y",
        frozenset({"ADV"}),
        ADJ_GENDER,
        ADJ,
        _E,
        frozenset({"y"}),
        label="adverb -e/-ě vs -y",
    ),
    Competition(
        "ADV -o",
        frozenset({"ADV"}),
        ADJ_GENDER,
        ADJ,
        _E,
        frozenset({"o"}),
        label="adverb -e/-ě vs -o",
    ),
)
# Competitions whose -ové/-é is the paradigm's for -a masculines.
NOM_PL_COMPETITIONS: Final[frozenset[str]] = frozenset(
    {"NOM.PL m.an", "NOM.PL m.an -ové"}
)
COMPETITION_BY_NAME: Final[Mapping[str, Competition]] = {
    c.name: c for c in COMPETITIONS
}
# The tested cells of (c), with the cells left out with them.
SPELLOUT_CELLS: Final[Mapping[str, frozenset[str]]] = {
    "LOC.SG m.in": frozenset({"LOC.SG"}),
    "LOC.SG n": frozenset({"LOC.SG"}),
    "NOM.PL m.an": frozenset({"NOM.PL", "VOC.PL"}),
    "LOC.PL m.an": frozenset({"LOC.PL"}),
    "LOC.PL m.in": frozenset({"LOC.PL"}),
    "LOC.PL n": frozenset({"LOC.PL"}),
    "VOC.SG m.an": frozenset({"VOC.SG"}),
}


def competition_observations(
    record: RootRecord, competition: Competition
) -> tuple[tp.Observation, ...]:
    """(b) observations: ``alternates`` = the realization is in slot A.
    An -a masculine has none in the NOM.PL competitions (module
    docstring)."""
    if competition.name in NOM_PL_COMPETITIONS and record.a_masculine:
        return ()
    a, b = competition.slots(record.target_class)
    out: list[tp.Observation] = []
    for item in record.realizations:
        if (
            item.cell not in competition.cells
            or item.gender not in competition.genders
        ):
            continue
        if item.ending in a:
            out.append(_observation(item, True, "P"))
        elif item.ending in b:
            out.append(_observation(item, False, "P"))
    return tuple(out)


# ---------------------------------------------------------------------------
# (c) Alternation evidence
# ---------------------------------------------------------------------------


def _derivative_counts(item: DerivativeEvidence, evidence_set: str) -> bool:
    if evidence_set == "E0" or item.outcome != "pal":
        return False
    if evidence_set == "E1":
        return item.suffix in VERB_SUFFIXES
    return True


def alternation_evidence(
    records: Sequence[RootRecord],
    tau: str,
    evidence_set: str,
    left_out: frozenset[str],
) -> bool:
    """Does any of *records* (a root's noun and adjective records) show
    the reflex at τ, outside the cells *left_out* (module docstring)?"""
    for record in records:
        for item in record.realizations:
            if (
                item.cell not in left_out
                and item.allomorph == ortho.P
                and item.mutation is True
                and item.attested_at(tau)
            ):
                return True
        for deriv in record.derivatives:
            if _derivative_counts(deriv, evidence_set) and deriv.attested_at(
                tau
            ):
                return True
    return False


# ---------------------------------------------------------------------------
# R17
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class R17Def:
    """One row definition of the fixed rule family (R17)."""

    table: str
    target: str
    trigger: str
    onsets: frozenset[str]
    memo: str
    derivatives_only: bool = False

    @property
    def key(self) -> str:
        return f"r17:{'|'.join(sorted(self.onsets))}:{self.derivatives_only}"


_FRONT: Final[frozenset[str]] = frozenset({FRONT_MID, FRONT_HIGH})
_HIGH: Final[frozenset[str]] = frozenset({FRONT_HIGH})
_MID: Final[frozenset[str]] = frozenset({FRONT_MID})
DORSAL: Final[str] = "dorsal"
R17_TARGETS: Final[Mapping[str, Member]] = {
    DORSAL: _klass(ortho.VELAR),
    "coronal_inferred": _klass(ortho.CORONAL),
    "r": _klass(ortho.R_CLASS),
    "sibilant": _klass(ortho.SIBILANT),
    "lateral": _klass(ortho.LATERAL),
}
_A: Final[str] = "dorsal/_[+front V]"
_B: Final[str] = "dorsal/_[+high,+front]"
_C: Final[str] = "coronal/_[+high,+front]"
_NT: Final[str] = ": non-trigger cell"

R17_DEFS: Final[tuple[R17Def, ...]] = (
    R17Def("R17a", DORSAL, "-e/-i", _FRONT, _A),
    R17Def("R17a", DORSAL, "-e", _MID, _A),
    R17Def("R17a", DORSAL, "-i", _HIGH, _A),
    R17Def(
        "R17a",
        DORSAL,
        "-e/-i/-y",
        _FRONT | {Y_VOWEL},
        _A + " (y = [ɪ] counted as front)",
    ),
    R17Def("R17b", DORSAL, "-i", _HIGH, _B),
    R17Def("R17b", DORSAL, "-e", _MID, _B + _NT),
    R17Def(
        "R17b",
        DORSAL,
        "C-suffix",
        frozenset({CONSONANT}),
        _B + _NT + " (consonant-initial suffixes: -ní, -ka, -ský)",
        True,
    ),
    R17Def(
        "R17b",
        DORSAL,
        "back-V suffix",
        frozenset({BACK}),
        _B + _NT + " (back-vowel suffixes: -an, -ovat, -ový)",
        True,
    ),
    *(
        item
        for target in ("coronal_inferred", "r", "sibilant", "lateral")
        for item in (
            R17Def("R17c", target, "-i", _HIGH, f"{_C}, {target}"),
            R17Def("R17c", target, "-e", _MID, f"{_C}, {target}{_NT}"),
        )
    ),
)


def _form_reflex(record: RootRecord) -> dict[str, bool | None]:
    """Reflex per form: True when any cell of the form shows it.

    A form that fills several cells is one string with one
    pronunciation; for inferred coronals the data layer marks the
    change only in the P cell (studenti NOM.PL True, the identical
    VOC.PL False), so the form-level value is used. For the same
    reason the change of t/d/n is inferred before ě/i/í in every
    aligned form, not only in P cells (kmen: kmeni [ɲɪ]): R17c's
    coronal rows, like the rule TP's, have e = 0 by construction.
    """
    inferred = record.target_class == ortho.CORONAL
    out: dict[str, bool | None] = {}
    for item in record.realizations:
        value = item.mutation if item.aligned else None
        if inferred and item.aligned and item.ending[:1] in INFERRED_FRONT:
            value = True
        if item.form not in out or value is True:
            out[item.form] = value
        elif out[item.form] is None and value is False:
            out[item.form] = False
    return out


def r17_observations(
    record: RootRecord, onsets: frozenset[str], evidence_set: str
) -> tuple[tp.Observation, ...]:
    """Realizations (and, from E1, derivatives) before a vowel of
    *onsets*, one per form; ``alternates`` = the stem shows the reflex
    (:func:`_form_reflex`). Derivative outcomes ``nontarget``
    (invisible coronal changes) are unknown."""
    if record.pos == NOUN and record.declension == ortho.ADJECTIVAL:
        return ()  # adjectival declension: the adjective rows (docstring)
    reflex = _form_reflex(record)
    out: list[tp.Observation] = []
    seen: set[str] = set()
    for item in record.realizations:
        if item.onset not in onsets or item.cell == "ADV" or item.form in seen:
            continue
        seen.add(item.form)
        out.append(_observation(item, reflex[item.form], "P"))
    out += _r17_derivatives(record, onsets, evidence_set)
    return tuple(out)


def _r17_derivatives(
    record: RootRecord, onsets: frozenset[str], evidence_set: str
) -> list[tp.Observation]:
    if evidence_set == "E0":
        return []
    out: list[tp.Observation] = []
    for deriv in record.derivatives:
        if deriv.onset not in onsets:
            continue
        if evidence_set == "E1" and deriv.suffix not in VERB_SUFFIXES:
            continue
        alternates = (
            None if deriv.outcome == "nontarget" else deriv.outcome == "pal"
        )
        out.append(
            tp.Observation(
                alternates=alternates,
                attested_by_tau=deriv.attested_by_tau(),
                allomorph="P",
                form=deriv.derivative,
                count=deriv.count,
            )
        )
    return out


def r17_derivative_observations(
    record: RootRecord, onsets: frozenset[str]
) -> tuple[tp.Observation, ...]:
    """E2 derivative observations only (the consonant- and
    back-vowel-initial suffix rows)."""
    return tuple(_r17_derivatives(record, onsets, "E2"))


# ---------------------------------------------------------------------------
# NDE
# ---------------------------------------------------------------------------

NDE_CLASSES: Final[tuple[str, ...]] = (
    "tautomorphemic_dorsal",
    "loan_hard_coronal",
)


def nde_observations(record: RootRecord) -> tuple[tp.Observation, ...]:
    """A root-internal velar or hard coronal before a front vowel never
    changes: every attested form is a non-alternation (R16)."""
    return tuple(
        _observation(item, False, "P")
        for item in record.realizations
        if item.cell != "ADV"
    )


# ---------------------------------------------------------------------------
# Root observations
# ---------------------------------------------------------------------------


def root_observation(
    record: RootRecord,
    evidence_set: str,
    basis: str = "lemma",
) -> tp.RootObservation:
    """Every trigger of one root for one evidence set (rule triggers
    and competitions do not depend on it; R17 triggers do)."""
    triggers: dict[str, Sequence[tp.Observation]] = {}
    for trigger in RULE_TRIGGERS:
        if trigger.pos == record.pos:
            triggers[trigger.key] = rule_observations(record, trigger)
    for competition in COMPETITIONS:
        if competition.pos == record.pos:
            triggers[competition.key] = competition_observations(
                record, competition
            )
    for definition in R17_DEFS:
        if definition.key in triggers:
            continue
        if definition.derivatives_only:
            triggers[definition.key] = r17_derivative_observations(
                record, definition.onsets
            )
        else:
            triggers[definition.key] = r17_observations(
                record, definition.onsets, evidence_set
            )
    if record.nde_class:
        triggers["nde"] = nde_observations(record)
    return tp.RootObservation(
        root_id=record.key,
        target_class=record.target_class,
        frequency=float(record.frequency(basis)),
        triggers=triggers,
    )


def members(records: Iterable[RootRecord], *tests: Member) -> list[RootRecord]:
    """Records passing every test, in input order."""
    return [r for r in records if all(test(r) for test in tests)]
