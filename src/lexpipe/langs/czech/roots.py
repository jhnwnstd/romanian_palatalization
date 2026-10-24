"""Root records for the Czech analyses (design R4, R15, R16, R22, R24).

The TP and the derivational analyses count **root types** (R4: TP
counts root types per cell; learner evidence pools all forms of a
``root_id``). :func:`load_roots` builds one :class:`RootRecord` per
(``root_id``, part of speech) from the data layer and the derivations
stage:

``czech_lexicon_complete.csv``  lexeme features (target class, stem,
                                R26 suffix target, NDE class, both R24
                                frequency bases); in-domain nouns and
                                adjectives only, affix entries left out
``czech_cells.csv``             primary cell forms with their P/N class,
                                mutation and attestation (R11); curated
                                rejections dropped
``nc/czech_lexeme_strata.csv``  the R22 stratum; only ``stratum_wikt``
                                (Wiktionary) is read, so tracked tables
                                stay free of NC data (R9)
``nc/czech_derivations.csv``    derivatives of the root's lexemes (E1,
                                E2), with their basis; the caller picks
                                the basis

Adjective adverbs (hezký → hezky ~ hezce) are read from the English
Wiktionary headword line (``cs-adj|adv=``) and counted in the corpus
(:func:`adverb_realizations`); they are cell ``ADV``.

Realizations keep what the data layer decided: ``allomorph`` is the
§6 P/N class (P palatalizing, N non-palatalizing, ``none`` no
opportunity), ``n_kind`` separates N allomorphs that compete with a P
ending (``allomorph``: -ové, -é, -u, -ách) from non-trigger front
vowels (``nontrigger``: -em, -ech, coronal VOC -e), and ``mutation`` is
True when the form shows the stem final's palatal reflex (inferred for
t/d/n before ě/i/í).

Strict thresholds (:data:`STRICT_TAUS`): ``ge1_strict``,
``ge0.5pm_strict`` and ``ge2pm_strict`` attest a cell form only when
it reaches the corpus threshold **and** is neither ``form_ambiguous``
(válku is válka's, not válek's) nor ``cell_confounded`` (LOC.SG parku
is also GEN.SG and DAT.SG; VOC kluku also DAT/LOC.SG): its tokens then
attest the cell's allomorph choice (the data layer's ``*_strict``
summaries). Adverbs and derivatives have no such confound; their
strict flags equal the plain ones.

Frequency bases (R24, :meth:`RootRecord.frequency`): ``lemma`` (the
lemma form's tokens), ``paradigm`` (the summed tokens of the attested
forms) and ``paradigm_unambiguous`` (the same without forms that
belong to another lexeme: jak 'yak' has the conjunction's lemma
tokens).
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Final

from lexpipe import csvio
from lexpipe.attest import TAUS, AttestationTable, attested_count
from lexpipe.langs.czech import jsonio, ortho
from lexpipe.langs.czech.adjudication import REJECTED
from lexpipe.langs.czech.attestation import CORPUS_TAUS
from lexpipe.langs.czech.etymology import UNKNOWN

NOUN: Final[str] = "noun"
ADJ: Final[str] = "adj"
ADV_CELL: Final[str] = "ADV"
VOWELS: Final[frozenset[str]] = frozenset("aáeéěiíoóuúůyý")

# Onset of an ending or suffix (R17): the class of its first letter.
FRONT_MID: Final[str] = "front_mid"
FRONT_HIGH: Final[str] = "front_high"
Y_VOWEL: Final[str] = "y"
BACK: Final[str] = "back"
CONSONANT: Final[str] = "C"
NO_ENDING: Final[str] = "none"
BASELINE_CLASSES: Final[frozenset[str]] = frozenset(
    {ortho.LABIAL, ortho.LATERAL, ortho.SIBILANT}
)
STRICT_SUFFIX: Final[str] = "_strict"
STRICT_TAUS: Final[tuple[str, ...]] = tuple(
    f"{tau}{STRICT_SUFFIX}" for tau in CORPUS_TAUS
)
# Every threshold of the TP tables: dictionary, the corpus thresholds
# and their strict variants (module docstring).
ALL_TAUS: Final[tuple[str, ...]] = (*TAUS, *STRICT_TAUS)
FREQ_BASES: Final[tuple[str, ...]] = (
    "lemma",
    "paradigm",
    "paradigm_unambiguous",
)


def plain_tau(tau: str) -> str:
    """The corpus threshold of a strict τ (``ge1_strict`` → ``ge1``)."""
    return tau[: -len(STRICT_SUFFIX)] if tau.endswith(STRICT_SUFFIX) else tau


def onset(ending: str) -> str:
    """Onset class of an ending: ``front_mid`` (e ě é), ``front_high``
    (i í), ``y`` (y ý, phonetically [ɪ iː]), ``back``, ``C`` or
    ``none``."""
    first = ending[:1].lower()
    if not first:
        return NO_ENDING
    if first in "eěé":
        return FRONT_MID
    if first in "ií":
        return FRONT_HIGH
    if first in "yý":
        return Y_VOWEL
    if first in VOWELS:
        return BACK
    return CONSONANT


def monosyllable(lemma: str) -> bool:
    """One vowel nucleus (ou/au/eu count once); a lemma without a vowel
    has a syllabic liquid (vlk, krk) and counts as one syllable."""
    lowered = lemma.lower().replace("ou", "o").replace("au", "a")
    lowered = lowered.replace("eu", "e")
    nuclei = 0
    previous = False
    for char in lowered:
        is_vowel = char in VOWELS
        if is_vowel and not previous:
            nuclei += 1
        previous = is_vowel
    return nuclei <= 1


@dataclass(frozen=True, slots=True)
class Realization:
    """One primary cell form of a root (see the module docstring)."""

    cell: str
    form: str
    ending: str
    gender: str
    allomorph: str
    n_kind: str
    mutation: bool | None
    aligned: bool
    count: int
    attested: tuple[bool, ...]
    ambiguous: bool
    confounded: bool

    @property
    def strict_ok(self) -> bool:
        """Its tokens attest this cell (module docstring)."""
        return not (self.ambiguous or self.confounded)

    def attested_by_tau(self) -> dict[str, bool]:
        """``{tau: bool}`` over :data:`ALL_TAUS`."""
        plain = dict(zip(CORPUS_TAUS, self.attested))
        return {
            "dictionary": True,
            **plain,
            **{
                f"{tau}{STRICT_SUFFIX}": value and self.strict_ok
                for tau, value in plain.items()
            },
        }

    def attested_at(self, tau: str) -> bool:
        if tau == "dictionary":
            return True
        base = plain_tau(tau)
        value = self.attested[CORPUS_TAUS.index(base)]
        if base != tau:
            return value and self.strict_ok
        return value

    @property
    def onset(self) -> str:
        return onset(self.ending)


@dataclass(frozen=True, slots=True)
class DerivativeEvidence:
    """A derivative of one of the root's lexemes (derivations stage)."""

    derivative: str
    suffix: str
    suffix_class: str
    pair: str
    outcome: str
    palatalizes_base: bool
    count: int
    attested: tuple[bool, ...]
    basis: str
    onset: str
    stratum_wikt: str

    def attested_by_tau(self) -> dict[str, bool]:
        plain = dict(zip(CORPUS_TAUS, self.attested))
        return {
            "dictionary": True,
            **plain,
            **{f"{tau}{STRICT_SUFFIX}": v for tau, v in plain.items()},
        }

    def attested_at(self, tau: str) -> bool:
        if tau == "dictionary":
            return True
        return self.attested[CORPUS_TAUS.index(plain_tau(tau))]


@dataclass
class RootRecord:
    """A root type of one part of speech (see the module docstring)."""

    root_id: str
    pos: str
    lemma: str
    lexeme_ids: list[str]
    target_class: str
    stem: str
    stem_final: str
    cluster: str
    readings: list[str]
    suffix_target: bool
    lemma_suffix: str
    nde_class: str
    stratum: str
    lemma_freq: float
    paradigm_freq: float
    realizations: list[Realization] = field(default_factory=list)
    derivatives: list[DerivativeEvidence] = field(default_factory=list)
    declension: str = ortho.HARD
    unambiguous_freq: float = 0.0

    @property
    def key(self) -> str:
        """Unique id of the record: the root_id, plus ``#adj`` for the
        adjective record of a root that also has a noun."""
        return self.root_id if self.pos == NOUN else f"{self.root_id}#adj"

    @property
    def monosyllable(self) -> bool:
        return monosyllable(self.lemma)

    @property
    def loan(self) -> bool:
        return self.stratum == "borrowed whole"

    @property
    def in_baseline(self) -> bool:
        """A hard labial, l, s or z stem: the §6 baseline. Soft paradigms
        (učitel: učitelé, učiteli) have endings of their own and stay
        out, as do adjectival nouns."""
        return (
            self.target_class in BASELINE_CLASSES
            and self.declension == ortho.HARD
        )

    @property
    def a_masculine(self) -> bool:
        """An animate masculine noun in -a (předseda, kolega, sluha,
        monarcha): its NOM.PL -ové/-é is its paradigm's ending, not a
        choice between -i and -ové (IJP 226)."""
        return (
            self.pos == NOUN
            and self.lemma.endswith("a")
            and "m.an" in self.readings
            and "f" not in self.readings
        )

    def frequency(self, basis: str) -> float:
        """Token frequency on an R24 basis (:data:`FREQ_BASES`)."""
        if basis == "lemma":
            return self.lemma_freq
        if basis == "paradigm_unambiguous":
            return self.unambiguous_freq
        return self.paradigm_freq

    def genders(self) -> list[str]:
        out: list[str] = []
        for item in self.realizations:
            if item.gender not in out:
                out.append(item.gender)
        return out


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def _float(text: str) -> float:
    return float(text) if text.strip() else 0.0


def _bool(text: str) -> bool:
    return csvio.parse_bool(text) is True


def _flags(row: Mapping[str, str]) -> tuple[bool, ...]:
    return tuple(_bool(row[f"attested_{tau}"]) for tau in CORPUS_TAUS)


def realization_from_row(row: Mapping[str, str]) -> Realization:
    """A :class:`Realization` from a ``czech_cells.csv`` row."""
    mutation = csvio.parse_bool(row["mutation"])
    return Realization(
        cell=row["cell"],
        form=row["form"],
        ending=row["ending"],
        gender=row["gender_class"],
        allomorph=row["allomorph_class"],
        n_kind=row["n_kind"],
        mutation=mutation,
        aligned=_bool(row["aligned"]),
        count=int(row["count"] or 0),
        attested=_flags(row),
        ambiguous=_bool(row["form_ambiguous"]),
        confounded=_bool(row["cell_confounded"]),
    )


def _record(row: Mapping[str, str], stratum: str) -> RootRecord:
    return RootRecord(
        root_id=row["root_id"],
        pos=row["pos"],
        lemma=row["lemma"],
        lexeme_ids=[row["lexeme_id"]],
        target_class=row["target_class"],
        stem=row["stem"],
        stem_final=row["stem_final"],
        cluster=row["cluster"],
        readings=csvio.split_list(row["readings"]),
        suffix_target=_bool(row["target_is_suffix"]),
        lemma_suffix=row["lemma_suffix"],
        nde_class=row["nde_class"],
        stratum=stratum,
        lemma_freq=_float(row["lemma_count"]),
        paradigm_freq=_float(row["paradigm_sum"]),
        declension=row["declension"],
        unambiguous_freq=_float(row.get("paradigm_sum_unambiguous", "")),
    )


def _merge(record: RootRecord, row: Mapping[str, str]) -> None:
    """Add a further lexeme of the same root and part of speech."""
    record.lexeme_ids.append(row["lexeme_id"])
    record.lemma_freq = max(record.lemma_freq, _float(row["lemma_count"]))
    record.paradigm_freq = max(
        record.paradigm_freq, _float(row["paradigm_sum"])
    )
    record.unambiguous_freq = max(
        record.unambiguous_freq,
        _float(row.get("paradigm_sum_unambiguous", "")),
    )
    for reading in csvio.split_list(row["readings"]):
        if reading not in record.readings:
            record.readings.append(reading)


def load_strata(path: Path) -> dict[str, str]:
    """``lexeme_id -> stratum_wikt`` (R22, Wiktionary only)."""
    if not path.is_file():
        return {}
    return {
        row["lexeme_id"]: row["stratum_wikt"] or UNKNOWN
        for row in jsonio.iter_csv(path)
    }


def load_records(
    lexicon_complete: Path, strata: Mapping[str, str]
) -> tuple[dict[str, RootRecord], dict[str, str]]:
    """Records keyed by :attr:`RootRecord.key`, and lexeme → key.

    Lexemes are visited in lexeme_id order, so a root's representative
    lemma and target class are those of its first lexeme.
    """
    records: dict[str, RootRecord] = {}
    lexeme_key: dict[str, str] = {}
    for row in jsonio.iter_csv(lexicon_complete):
        if row["domain"] != "in" or row["pos"] not in (NOUN, ADJ):
            continue
        if row["lemma"].startswith("-") or row["lemma"].endswith("-"):
            continue
        key = row["root_id"] if row["pos"] == NOUN else f"{row['root_id']}#adj"
        stratum = strata.get(row["lexeme_id"], UNKNOWN)
        if key in records:
            _merge(records[key], row)
        else:
            records[key] = _record(row, stratum)
        lexeme_key[row["lexeme_id"]] = key
    return records, lexeme_key


def add_cells(
    records: Mapping[str, RootRecord],
    lexeme_key: Mapping[str, str],
    cells: Path,
) -> int:
    """Attach primary, not curated-rejected cell forms (one per root,
    cell, form and gender); return the number attached."""
    seen: set[tuple[str, str, str, str]] = set()
    attached = 0
    for row in jsonio.iter_csv(cells):
        if row["primary"] != "True" or row["curated_status"] == REJECTED:
            continue
        key = lexeme_key.get(row["lexeme_id"])
        if key is None:
            continue
        ident = (key, row["cell"], row["form"], row["gender_class"])
        if ident in seen:
            continue
        seen.add(ident)
        records[key].realizations.append(realization_from_row(row))
        attached += 1
    return attached


def derivative_from_row(row: Mapping[str, str]) -> DerivativeEvidence:
    """A :class:`DerivativeEvidence` from a ``czech_derivations.csv`` row
    (its onset is that of the ending it was decomposed with)."""
    return DerivativeEvidence(
        derivative=row["derivative"],
        suffix=row["suffix"],
        suffix_class=row["suffix_class"],
        pair=row["suffix_pair"],
        outcome=row["outcome"],
        palatalizes_base=_bool(row["suffix_palatalizes_base"]),
        count=int(row["count"] or 0),
        attested=_flags(row),
        basis=row["basis"],
        onset=onset(row["ending"]),
        stratum_wikt=row["derivative_stratum_wikt"] or UNKNOWN,
    )


def add_derivatives(
    records: Mapping[str, RootRecord],
    lexeme_key: Mapping[str, str],
    derivations: Path,
    bases: Sequence[str],
) -> int:
    """Attach derivatives whose ``basis`` is one of *bases* (dropping a
    derivative repeated through two lexemes of one root)."""
    allowed = frozenset(bases)
    seen: set[tuple[str, str, str]] = set()
    attached = 0
    if not derivations.is_file():
        return 0
    for row in jsonio.iter_csv(derivations):
        if row["basis"] not in allowed:
            continue
        key = lexeme_key.get(row["base_lexeme_id"])
        if key is None:
            continue
        ident = (key, row["suffix"], row["derivative"])
        if ident in seen:
            continue
        seen.add(ident)
        records[key].derivatives.append(derivative_from_row(row))
        attached += 1
    return attached


# ---------------------------------------------------------------------------
# Adverbs (cell ADV)
# ---------------------------------------------------------------------------


def adverb_links(en_entries: Path) -> dict[str, list[str]]:
    """``adjective lemma -> adverbs`` from ``cs-adj|adv=`` (source
    order; ``-`` and multiword values dropped)."""
    out: dict[str, list[str]] = defaultdict(list)
    for record in jsonio.read_jsonl(en_entries):
        if record.get("pos") != ADJ:
            continue
        word = str(record.get("word") or "")
        for template in record.get("head_templates") or ():
            if template.get("name") != "cs-adj":
                continue
            args = template.get("args") or {}
            for item in str(args.get("adv", "")).split(","):
                adverb = item.split("<")[0].strip()
                if adverb and adverb != "-" and " " not in adverb:
                    if adverb not in out[word]:
                        out[word].append(adverb)
    return out


def adverb_realization(
    stem: str,
    target_class: str,
    adverb: str,
    table: AttestationTable,
) -> Realization:
    """The ADV realization of one adverb: -e/-ě is the P slot (with the
    stem final's reflex: hezce, hluboce, staře), -y/-o the N slot
    (hezky, česky, lehko)."""
    alignment = ortho.align(stem, adverb)
    ending = alignment.ending if alignment.aligned else ""
    allomorph = ortho.NONE
    if ending in ("e", "ě"):
        allomorph = ortho.P
    elif ending in ("y", "o"):
        allomorph = ortho.N
    mutation: bool | None = None
    if alignment.aligned:
        mutation = alignment.mutated
        if target_class == ortho.CORONAL and allomorph == ortho.P:
            mutation = ending == "ě"
    count = table.count(adverb, case_sensitive=True)
    count += table.count(adverb[:1].upper() + adverb[1:], case_sensitive=True)
    size = table.total_tokens
    return Realization(
        cell=ADV_CELL,
        form=adverb,
        ending=ending,
        gender=ADJ,
        allomorph=allomorph,
        n_kind=ortho.N_ALLOMORPH if allomorph == ortho.N else "",
        mutation=mutation,
        aligned=alignment.aligned,
        count=count,
        attested=tuple(
            attested_count(count, tau, size) for tau in CORPUS_TAUS
        ),
        ambiguous=False,
        confounded=False,
    )


def add_adverbs(
    records: Mapping[str, RootRecord],
    links: Mapping[str, Sequence[str]],
    table: AttestationTable,
) -> int:
    """Attach ADV realizations to adjective records."""
    attached = 0
    for record in records.values():
        if record.pos != ADJ:
            continue
        for adverb in links.get(record.lemma, ()):
            record.realizations.append(
                adverb_realization(
                    record.stem, record.target_class, adverb, table
                )
            )
            attached += 1
    return attached


def ordered(records: Mapping[str, RootRecord]) -> list[RootRecord]:
    """Records sorted by key (the order every table iterates in)."""
    return [records[key] for key in sorted(records)]


def select(
    records: Iterable[RootRecord], pos: str, classes: Iterable[str]
) -> list[RootRecord]:
    """Records of one part of speech whose target class is listed."""
    wanted = frozenset(classes)
    return [r for r in records if r.pos == pos and r.target_class in wanted]


TAU_LIST: Final[tuple[str, ...]] = TAUS
