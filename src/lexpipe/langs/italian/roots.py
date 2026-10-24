"""Root records of the Italian analyses (R4, R15, R16, R19, R24).

One :class:`RootRecord` per root_id of the in-domain noun and adjective
lexemes (``italian_lexicon_complete.csv``), with every realization of
its stem-final consonant before a suffix vowel, tagged by evidence set
(R15):

    E0  the plural cells of the root's lexemes (``italian_cells.csv``):
        opportunity i → vowel i, e → e, a → back; filtered variants,
        curator-rejected and unaligned forms are left out, and so are
        the forms of a lexeme that is not a dorsal target in a dorsal
        root (a homograph merged into it)
    E1  E0 + the 2sg, 3sg, 1pl and 3sg-subjunctive cells of the
        in-domain velar-stem verbs linked to the root
        (``italian_verbs.csv``: amico ↔ amicare, traffico ↔
        trafficare, bianco ↔ bianchire)
    E2  E1 + the derivatives linked to the root (``italian_derivations.
        csv``, accepted links with a palatal, velar or back edge; the
        ``named`` tier by default, :data:`DEFAULT_DERIV_TIERS`) that
        are internal formations or of unknown stratum: learned
        borrowings and inherited derivatives (cecità, sporcizia,
        porcino) are excluded by stratum (design §5) and kept only in
        the sensitivity set ``E2_all`` (every stratum)

A realization's ``palatal`` is the alternation the learner sees: the
plural's mutation (dorsal targets), the verb cell's or derivative's
stem-edge reading; always False for coronals, NDE classes and lexically
palatal stems (no change), None for a back vowel.

:func:`observation` turns a record into a
:class:`lexpipe.tp.RootObservation` for an evidence set and a part of
speech scope (``all``; ``noun`` / ``adj``: E0 cells of that part of
speech only, E1/E2 unchanged), with the trigger cells of
:data:`TRIGGERS`:

    -i     P = vowel i; N = vowel e (the root meets the other front
           vowel instead); nonfront = back vowel
    -e     P = vowel e; N = vowel i; nonfront = back vowel
    -i/-e  P = vowel i or e; nonfront = back vowel

so the core's statuses (R16) read: no_evidence (nothing attested in the
cell), shielded (only N or back realizations attested), alternates,
nonalt, or a doublet coded by the case-sensitive token counts of its
forms (R19, as ``features.doublet_class``; dominance needs
:data:`.features.DOUBLET_MIN_TOKENS` tokens). Each observation carries
its paradigm cell (:func:`cell_key`: the -i or -e plural, a verb cell,
a derivative's suffix), so a doublet is coded within one cell and a
root whose cells disagree — amici but amiche under ``-i/-e``, the
plural medici but the verb form medichi in E1 — is an exception to
both rules (``variable``; :mod:`lexpipe.tp`). Equal (form, vowel)
realizations count once (amici of amico N and amico ADJ).

Frequencies (R24): ``lemma`` is the largest folded lemma-form count of
the scope's lexemes, ``paradigm`` the root's folded paradigm sum
(``freq_root_folded``). Lemma counts of homographs are shared (antico
'ancient' and àntico 'anterior' both get the tokens of "antico").

Verb units (:class:`VerbUnit`) are the velar-stem verbs grouped by
infinitive (amicare and amicarsi are one unit), for the verb TP rows;
their observations use the same trigger cells.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from types import MappingProxyType
from typing import Final

from lexpipe import attest, csvio, tp
from lexpipe.langs.italian import derivations, keys, morphology

EVIDENCE_SETS: Final[tuple[str, ...]] = ("E0", "E1", "E2")
# E2 with the borrowed-whole and inherited derivatives (sensitivity).
E2_ALL: Final[str] = "E2_all"
# Derivative strata E2 leaves out (design §5: learned borrowings are
# excluded by stratum; inherited derivatives are Latin formations).
E2_EXCLUDED_STRATA: Final[frozenset[str]] = frozenset(
    {"borrowed whole", "inherited"}
)
TRIGGERS: Final[tuple[str, ...]] = ("-i", "-e", "-i/-e")
SCOPES: Final[tuple[str, ...]] = ("all", "noun", "adj")
FREQ_BASES: Final[tuple[str, ...]] = ("lemma", "paradigm")
DEFAULT_DERIV_TIERS: Final[tuple[str, ...]] = (derivations.TIER_NAMED,)
_TRIGGER_VOWELS: Final[Mapping[str, frozenset[str]]] = MappingProxyType(
    {
        "-i": frozenset({"i"}),
        "-e": frozenset({"e"}),
        "-i/-e": frozenset({"i", "e"}),
    }
)
_FRONT: Final[frozenset[str]] = frozenset({"i", "e"})
_RANK: Final[Mapping[str, int]] = MappingProxyType(
    {name: i for i, name in enumerate((*EVIDENCE_SETS, E2_ALL))}
)
# The dictionary threshold restricted to plurals an editor typed
# (provenance explicit): the Module:it-common defaults (every -ico/-aco/
# -ologo plural palatal, every other -co/-go plural velar) are not
# evidence of their own (R11). E0 only; a sample biased toward the
# plurals the module gets wrong.
EXPLICIT_TAU: Final[str] = "dictionary_explicit"
DORSAL: Final[str] = "dorsal"
CORONAL: Final[str] = "coronal"
FRONTSTEM: Final[str] = "frontstem"
NONE: Final[str] = "none"


@dataclass(frozen=True, slots=True)
class Realization:
    """The root's target before a suffix vowel in one form."""

    evidence: str
    scope: str
    source: str
    key: str
    vowel: str
    palatal: bool | None
    count_cs: int
    count_folded: int
    attested: Mapping[str, bool]
    ambiguous: bool


@dataclass
class RootRecord:
    """One root with its properties and realizations (module docstring)."""

    root_id: str
    lemma: str
    stem: str
    target_class: str
    seg: str
    target_letters: str
    has_noun: bool
    has_adj: bool
    velar_post_tonic: bool | None
    stress_class: str
    ico_suffix: bool
    iaco_suffix: bool
    suffixoid: str
    semantic_class: str
    orth_ico: bool
    deverbal_from_are: bool
    shares_1st_conj: bool
    stratum: str
    freq_lemma: dict[str, int]
    freq_paradigm: int
    realizations: list[Realization] = field(default_factory=list)

    @property
    def dorsal(self) -> bool:
        return self.target_class == DORSAL

    def frequency(self, basis: str, scope: str = "all") -> int:
        """Lemma or paradigm frequency (module docstring)."""
        if basis == "paradigm":
            return self.freq_paradigm
        return self.freq_lemma.get(scope, 0)


def _bool(text: str) -> bool:
    return csvio.parse_bool(text) is True


def _attested(
    row: Mapping[str, str], explicit: bool = False
) -> dict[str, bool]:
    """The row's τ flags, plus :data:`EXPLICIT_TAU` (*explicit*: an
    editor typed the form)."""
    flags = {tau: _bool(row.get(f"attested_{tau}", "")) for tau in attest.TAUS}
    flags[EXPLICIT_TAU] = explicit
    return flags


def _int(text: str) -> int:
    return int(text) if text else 0


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def _target_class(row: Mapping[str, str]) -> str:
    if _bool(row["dorsal_target"]):
        return DORSAL
    if row["nde_class"]:
        return f"nde:{row['nde_class']}"
    if row["target"] == CORONAL:
        return CORONAL
    if _bool(row["frontstem"]):
        return FRONTSTEM
    return NONE


def _primary_order(row: Mapping[str, str]) -> tuple[bool, bool, str]:
    return (
        not row["word"].endswith("o"),
        row["pos"] != "noun",
        row["lexeme_id"],
    )


def _record(root_id: str, rows: Sequence[Mapping[str, str]]) -> RootRecord:
    ordered = sorted(rows, key=_primary_order)
    head = ordered[0]
    key = keys.match_key(head["word"])
    stem = morphology.analyze_lemma(key)
    post = next(
        (r["velar_post_tonic"] for r in ordered if r["velar_post_tonic"]), ""
    )
    suffixoids = [
        r["lemma_suffix"] for r in ordered if r["suffix_kind"] == "suffixoid"
    ]
    strata = [
        r["stratum"]
        for r in ordered
        if r["stratum"] and r["stratum"] != "unknown"
    ]
    freq_lemma = {
        scope: max(
            (
                _int(r["freq_lemma_folded"])
                for r in rows
                if scope == "all" or r["pos"] == scope
            ),
            default=0,
        )
        for scope in SCOPES
    }
    return RootRecord(
        root_id=root_id,
        lemma=head["word"],
        stem=stem.stem,
        target_class=_target_class(head),
        seg=head["stem_final"],
        target_letters=(
            stem.target_letters if stem.target_letters else head["stem_final"]
        ),
        has_noun=any(r["pos"] == "noun" for r in rows),
        has_adj=any(r["pos"] == "adj" for r in rows),
        velar_post_tonic=csvio.parse_bool(post) if post else None,
        stress_class=next(
            (r["stress_class"] for r in ordered if r["stress_class"]), ""
        ),
        ico_suffix=any(
            _bool(r["target_is_suffix"]) and r["lemma_suffix"] == "-ico"
            for r in rows
        ),
        iaco_suffix=any(
            _bool(r["target_is_suffix"]) and r["lemma_suffix"] == "-iaco"
            for r in rows
        ),
        suffixoid=suffixoids[0] if suffixoids else "",
        semantic_class=next(
            (
                r["semantic_class"]
                for r in ordered
                if r["suffix_kind"] == "suffixoid" and r["semantic_class"]
            ),
            "",
        ),
        orth_ico=any(
            keys.match_key(r["word"]).endswith(("ico", "ica"))
            and r["ico_stressed"] != "True"
            for r in rows
        ),
        deverbal_from_are=any(_bool(r["deverbal_from_are"]) for r in rows),
        shares_1st_conj=any(
            _bool(r["shares_root_with_1st_conj_verb"]) for r in rows
        ),
        stratum=strata[0] if strata else "unknown",
        freq_lemma=freq_lemma,
        freq_paradigm=max(_int(r["freq_root_folded"]) for r in rows),
    )


def load_records(lexicon_complete: Path) -> dict[str, RootRecord]:
    """Records (without realizations) of every in-domain noun/adjective
    root, keyed by root_id."""
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in csvio.read_rows(lexicon_complete):
        if _bool(row["in_domain"]):
            grouped[row["root_id"]].append(row)
    return {
        root_id: _record(root_id, grouped[root_id])
        for root_id in sorted(grouped)
    }


def _cell_vowel(opportunity: str) -> str:
    return {"i": "i", "e": "e", "a": "back"}.get(opportunity, "")


def _dorsal_lexeme(row: Mapping[str, str]) -> bool:
    """Is the cell row's own lexeme a dorsal target (not an NDE stem, not
    lexically palatal)?"""
    return (
        row["target"] == DORSAL
        and not row["nde_class"]
        and not _bool(row["frontstem"])
    )


def add_cells(records: Mapping[str, RootRecord], cells: Path) -> int:
    """E0 realizations from the cell table; returns how many. A dorsal
    root takes only the forms of its dorsal-target lexemes: a lexically
    palatal homograph merged into it (dolce in dolco) would otherwise
    read as a velar kept."""
    added = 0
    for row in csvio.read_rows(cells):
        record = records.get(row["root_id"])
        vowel = _cell_vowel(row["opportunity"])
        if (
            record is None
            or not vowel
            or not _bool(row["in_domain"])
            or not _bool(row["aligned"])
            or _bool(row["filtered"])
            or _bool(row["adjudicated_out"])
            or (record.dorsal and not _dorsal_lexeme(row))
        ):
            continue
        palatal: bool | None = None
        if vowel in _FRONT:
            palatal = _bool(row["mutation"]) if record.dorsal else False
        record.realizations.append(
            Realization(
                evidence="E0",
                scope=row["pos"],
                source=row["cell"],
                key=keys.match_key(row["form"]),
                vowel=vowel,
                palatal=palatal,
                count_cs=_int(row["count_cs"]),
                count_folded=_int(row["count_folded"]),
                attested=_attested(row, row["provenance"] == "explicit"),
                ambiguous=_bool(row["form_ambiguous"]),
            )
        )
        added += 1
    return added


# Module:it-common's plural rule for -co/-go (make_plural): -ico/-aco
# (and -ologo) palatalize, every other -co/-go keeps the velar.
_MODULE_PALATAL_ENDINGS: Final[tuple[str, ...]] = ("ico", "aco", "ologo")


def dictionary_provenance(cells: Path) -> Counter[str]:
    """Where the dorsal -i plurals come from (R11; for the digest):
    over the in-domain dorsal-target lexemes' unfiltered, aligned -i
    forms, how many an editor typed (``explicit``) and how many are
    Module:it-common defaults, whether a default follows the module's
    -ico/-aco/-ologo rule, and how often the competing spelling of a
    default (amichi for amici: ``competitor``) has a lowercase token,
    or more tokens than the listed form — counting only competitors
    that are not forms of another word (dolci of dolce for dolco);
    the same over the forms attested at ≥ 1 token (keys ``ge1:…``)."""
    counts: Counter[str] = Counter()
    for row in csvio.read_rows(cells):
        if (
            not _bool(row["in_domain"])
            or not _dorsal_lexeme(row)
            or row["opportunity"] != "i"
            or not _bool(row["aligned"])
            or _bool(row["filtered"])
            or _bool(row["adjudicated_out"])
        ):
            continue
        scopes = [""]
        if _bool(row["attested_ge1"]):
            scopes.append("ge1:")
        provenance = row["provenance"]
        palatal = _bool(row["mutation"])
        module_rule = keys.match_key(row["word"]).endswith(
            _MODULE_PALATAL_ENDINGS
        )
        rival = (
            0
            if _bool(row["competitor_other_lemma"])
            else _int(row["competitor_count_cs"])
        )
        for scope in scopes:
            counts[f"{scope}forms"] += 1
            counts[f"{scope}{provenance}"] += 1
            if provenance != "module_default":
                continue
            counts[f"{scope}module_default:palatal:{palatal}"] += 1
            counts[f"{scope}module_default:follows_rule"] += (
                palatal == module_rule
            )
            if rival > 0:
                counts[f"{scope}module_default:competitor_attested"] += 1
            if rival > _int(row["count_cs"]):
                counts[f"{scope}module_default:competitor_more_frequent"] += 1
    return counts


@dataclass
class VerbUnit:
    """A velar-stem verb (reflexive variants merged) with its cells."""

    verb_key: str
    verb_ids: list[str]
    conj_class: str
    paradigm_type: str
    overapplies: bool
    stem_key: str
    target_letters: str
    noun_roots: list[str]
    stratum: str
    realizations: list[Realization] = field(default_factory=list)

    @property
    def frequency(self) -> int:
        seen: dict[str, int] = {}
        for item in self.realizations:
            seen.setdefault(item.key, item.count_folded)
        return sum(seen.values())


def load_verb_units(verbs: Path) -> dict[str, VerbUnit]:
    """In-domain velar-stem verbs by infinitive key, with their E1 cells."""
    units: dict[str, VerbUnit] = {}
    for row in csvio.read_rows(verbs):
        if not _bool(row["in_domain"]):
            continue
        unit = units.get(row["verb_key"])
        if unit is None:
            unit = VerbUnit(
                verb_key=row["verb_key"],
                verb_ids=[],
                conj_class=row["conj_class"],
                paradigm_type=row["paradigm_type"],
                overapplies=_bool(row["overapplies"]),
                stem_key=row["stem_key"],
                target_letters=row["target_letters"],
                noun_roots=csvio.split_list(row["noun_roots"]),
                stratum=row["stratum"],
            )
            units[row["verb_key"]] = unit
        if row["verb_id"] not in unit.verb_ids:
            unit.verb_ids.append(row["verb_id"])
        if not _bool(row["e1_cell"]) or row["edge"] == "unparsed":
            continue
        vowel = row["vowel"]
        unit.realizations.append(
            Realization(
                evidence="E1",
                scope="verb",
                source=f"verb:{row['cell']}",
                key=row["form_key"],
                vowel=vowel,
                palatal=(
                    (row["edge"] == "palatal") if vowel in _FRONT else None
                ),
                count_cs=_int(row["count_cs"]),
                count_folded=_int(row["count_folded"]),
                attested=_attested(row),
                ambiguous=_bool(row["form_ambiguous"]),
            )
        )
    return {key: units[key] for key in sorted(units)}


def add_verbs(
    records: Mapping[str, RootRecord], units: Mapping[str, VerbUnit]
) -> int:
    """E1 realizations: each unit's cells on its linked noun roots."""
    added = 0
    for unit in units.values():
        for root_id in unit.noun_roots:
            record = records.get(root_id)
            if record is None:
                continue
            record.realizations.extend(unit.realizations)
            added += len(unit.realizations)
    return added


def add_derivatives(
    records: Mapping[str, RootRecord],
    path: Path,
    tiers: Sequence[str] = DEFAULT_DERIV_TIERS,
) -> int:
    """E2 realizations from the accepted derivative links of *tiers*
    (derivatives of :data:`E2_EXCLUDED_STRATA` are tagged ``E2_all``:
    only that sensitivity set reads them)."""
    added = 0
    for row in csvio.read_rows(path):
        record = records.get(row["base_root_id"])
        if (
            record is None
            or row["link_status"] != derivations.ACCEPTED
            or row["link_tier"] not in tiers
            or row["edge"] not in ("palatal", "velar", "back")
        ):
            continue
        vowel = row["suffix_vowel"]
        record.realizations.append(
            Realization(
                evidence=(
                    E2_ALL
                    if row["deriv_stratum"] in E2_EXCLUDED_STRATA
                    else "E2"
                ),
                scope="deriv",
                source=f"deriv:{row['suffix']}",
                key=keys.match_key(row["derivative"]),
                vowel=vowel,
                palatal=(
                    (row["edge"] == "palatal") if vowel in _FRONT else None
                ),
                count_cs=_int(row["count_cs"]),
                count_folded=_int(row["count_folded"]),
                attested=_attested(row),
                ambiguous=False,
            )
        )
        added += 1
    return added


# ---------------------------------------------------------------------------
# Observations
# ---------------------------------------------------------------------------


def realizations_for(
    items: Iterable[Realization], evidence: str, scope: str = "all"
) -> list[Realization]:
    """The realizations of an evidence set and scope, each (form,
    vowel) once, earlier evidence first.

    Copies of a form within one evidence set (amici of amico N and
    amico ADJ) merge their τ flags: the form is attested at τ when any
    copy is, so a plural an editor typed for the noun keeps
    :data:`EXPLICIT_TAU` even when the adjective's copy is a module
    default (the other flags are the same form's and agree anyway).
    """
    limit = _RANK[evidence]
    chosen = [
        r
        for r in items
        if _RANK[r.evidence] <= limit
        and (r.evidence != "E0" or scope == "all" or r.scope == scope)
    ]
    chosen.sort(key=lambda r: (_RANK[r.evidence], r.key, r.vowel, r.source))
    seen: dict[tuple[str, str], Realization] = {}
    for item in chosen:
        key = (item.key, item.vowel)
        kept = seen.get(key)
        if kept is None:
            seen[key] = item
        elif kept.evidence == item.evidence:
            flags = {
                tau: kept.attested.get(tau, False)
                or item.attested.get(tau, False)
                for tau in (*kept.attested, *item.attested)
            }
            if flags != dict(kept.attested):
                seen[key] = replace(kept, attested=MappingProxyType(flags))
    return list(seen.values())


def _allomorph(trigger: str, vowel: str) -> str:
    if vowel not in _FRONT:
        return "nonfront"
    return "P" if vowel in _TRIGGER_VOWELS[trigger] else "N"


def cell_key(item: Realization) -> str:
    """The paradigm cell of a realization (module docstring): the -i or
    -e plural (``plural:i``: pl and m.pl alike), the verb cell
    (``verb:2sg``) or the derivative's suffix (``deriv:-ista``)."""
    if item.evidence == "E0":
        return f"plural:{item.vowel}"
    return item.source


def _attestation(item: Realization, unambiguous: bool) -> Mapping[str, bool]:
    if not unambiguous or not item.ambiguous:
        return item.attested
    return {tau: tau == "dictionary" for tau in attest.TAUS}


def observations(
    items: Sequence[Realization], trigger: str, *, unambiguous: bool = False
) -> tuple[tp.Observation, ...]:
    """The realizations as observations of one trigger cell."""
    out: list[tp.Observation] = []
    for item in items:
        allomorph = _allomorph(trigger, item.vowel)
        out.append(
            tp.Observation(
                alternates=item.palatal if allomorph == "P" else None,
                attested_by_tau=_attestation(item, unambiguous),
                allomorph=allomorph,
                form=item.key,
                count=item.count_cs,
                cell=cell_key(item),
            )
        )
    return tuple(out)


def observation(
    record: RootRecord,
    evidence: str,
    *,
    scope: str = "all",
    basis: str = "lemma",
    unambiguous: bool = False,
    with_cells: bool = True,
) -> tp.RootObservation:
    """A record as a :class:`lexpipe.tp.RootObservation` (module
    docstring); ``with_cells=False`` gives a frame root without cells."""
    items = (
        realizations_for(record.realizations, evidence, scope)
        if with_cells
        else []
    )
    return tp.RootObservation(
        root_id=record.root_id,
        target_class=record.target_class,
        frequency=float(record.frequency(basis, scope)),
        triggers={
            trigger: observations(items, trigger, unambiguous=unambiguous)
            for trigger in TRIGGERS
        },
    )


def unit_observation(
    unit: VerbUnit, *, unambiguous: bool = False
) -> tp.RootObservation:
    """A verb unit as a RootObservation (root_id ``verb:<infinitive>``)."""
    items = realizations_for(unit.realizations, "E1")
    return tp.RootObservation(
        root_id=f"verb:{unit.verb_key}",
        target_class=DORSAL,
        frequency=float(unit.frequency),
        triggers={
            trigger: observations(items, trigger, unambiguous=unambiguous)
            for trigger in TRIGGERS
        },
    )


@dataclass
class Records:
    """Every input of the analyses, loaded once."""

    roots: dict[str, RootRecord]
    units: dict[str, VerbUnit]
    counts: dict[str, int]


def load_all(
    lexicon_complete: Path,
    cells: Path,
    verbs: Path,
    derivations_csv: Path,
    *,
    deriv_tiers: Sequence[str] = DEFAULT_DERIV_TIERS,
) -> Records:
    """Records with E0, E1 and E2 realizations, plus the verb units."""
    roots = load_records(lexicon_complete)
    units = load_verb_units(verbs)
    counts = {
        "roots": len(roots),
        "E0": add_cells(roots, cells),
        "E1": add_verbs(roots, units),
        "E2": add_derivatives(roots, derivations_csv, deriv_tiers),
        "verb_units": len(units),
    }
    return Records(roots, units, counts)
