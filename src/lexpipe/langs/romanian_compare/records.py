"""Per-root records: the compact derivative of the Romanian inputs (R7).

The ``roots`` stage reads the lexicon, the verb-paradigm verdicts and
the Leipzig table once and writes one JSON object per noun root to
``data/romanian_compare/romanian_compare_roots.jsonl`` (gitignored:
it holds per-word Leipzig counts, R9), sorted by ``root_id``, with
sorted keys (R5). The first line is a header object carrying the
corpus size and the SHA-256 of every input, so the file is stamped
with its sources. Every later stage reads only this file (plus the
legacy outputs it compares with).

A root's scalar fields (gloss, plural, legacy columns) come from its
priority row (:func:`lexpipe.langs.romanian_compare.lexicon.group_nouns`);
derived verbs and adjectives are merged from all its rows. The parity
report counts lemmas whose rows disagree on the plural (none in the
current lexicon), the only case where a second plural would be lost.

Fields of a root (see :class:`RootRecord`): identity and lexical
attributes, the outcome-blind attributes of
:mod:`lexpipe.langs.romanian_compare.lexicon`, the plural cell (E0),
the 2sg / 3sg-subjunctive cells of the verb listed for the noun in
``data/verb_paradigm_verdict.csv`` (E1; 422 nouns have a row there),
the derived verbs and adjectives of the lexicon (E2), corpus counts
of every form, and the legacy columns the parity report needs.
"""

from __future__ import annotations

import csv
import json
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Final

from lexpipe.frequency import paradigm_frequency
from lexpipe.langs.romanian_compare import lexicon as lx
from lexpipe.langs.romanian_compare import morphology as mo
from lexpipe.langs.romanian_compare.corpus import RomanianCounts

SCHEMA: Final[str] = "romanian_compare.roots/1"
VERDICT_COLUMNS: Final[tuple[str, ...]] = (
    "noun",
    "verb",
    "stem_final",
    "opp",
    "ds",
    "verdict",
    "pres_2sg",
    "subj_3sg",
)


@dataclass(frozen=True, slots=True)
class VerbCell:
    """One E1 cell of the noun's verb (2sg present, 3sg subjunctive)."""

    verb: str
    cell: str
    form: str
    vowel: str
    palatal: bool
    augment: bool
    count: int


@dataclass(frozen=True, slots=True)
class Derivative:
    """A derived verb or adjective listed in the lexicon (E2).

    ``vowel`` and ``palatal`` are None when the form could not be read
    (:mod:`lexpipe.langs.romanian_compare.morphology`). For verbs,
    ``attest_count`` adds the past participle's tokens to the
    infinitive's (both forms of one lexeme; infinitives are rare in
    news text).
    """

    kind: str
    form: str
    suffix: str
    vowel: str | None
    palatal: bool | None
    count: int
    attest_count: int


@dataclass(frozen=True, slots=True)
class Candidate:
    """A generated -i verb candidate (R21) and its corpus tokens."""

    form: str
    palatal: bool
    count: int


@dataclass(frozen=True, slots=True)
class Legacy:
    """Legacy columns of the root's primary row (parity report)."""

    tp_in_domain: bool
    exception_reason: str
    opportunity: str
    mutation: bool | None
    orth_change: str
    target_is_suffix: bool
    lemma_suffix: str


@dataclass(frozen=True, slots=True)
class RootRecord:
    """Everything the analyses know about one noun root."""

    root_id: str
    lemma: str
    gloss: str
    gender: str
    seg: str
    cluster: str
    target_class: str
    nde: str
    flags: tuple[str, ...]
    stratum: str
    lexicon_source: str
    plural: str
    plural_allomorph: str
    plural_palatal: bool | None
    lemma_count: int
    plural_count: int
    paradigm_count: int
    verb: str
    audit_verdict: str
    verb_cells: tuple[VerbCell, ...]
    derivatives: tuple[Derivative, ...]
    candidates: tuple[Candidate, ...]
    legacy: Legacy
    duplicates: int = 1

    @property
    def is_target(self) -> bool:
        return self.target_class in ("dorsal", "coronal", "z")

    def frequency(self, basis: str) -> int:
        """Token frequency on a basis: ``lemma`` or ``paradigm``."""
        if basis == "lemma":
            return self.lemma_count
        if basis == "paradigm":
            return self.paradigm_count
        raise ValueError(f"unknown frequency basis {basis!r}")


# ---------------------------------------------------------------------------
# Building
# ---------------------------------------------------------------------------


def read_verdicts(path: Path) -> dict[str, dict[str, str]]:
    """``data/verb_paradigm_verdict.csv`` rows keyed by noun lemma.

    The file has one row per noun; a repeated noun keeps its first
    row in sorted (verb, verdict) order, so the choice is not
    input-order dependent.
    """
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        missing = [
            c for c in VERDICT_COLUMNS if c not in (reader.fieldnames or ())
        ]
        if missing:
            raise ValueError(f"{path}: missing columns {missing}")
        rows = sorted(
            reader, key=lambda r: (r["noun"], r["verb"], r["verdict"])
        )
    out: dict[str, dict[str, str]] = {}
    for row in rows:
        out.setdefault(row["noun"].strip(), dict(row))
    return out


def _verb_cells(
    verdict: Mapping[str, str] | None,
    seg: str,
    cluster: str,
    counts: RomanianCounts,
) -> tuple[str, str, tuple[VerbCell, ...]]:
    if verdict is None:
        return "", "", ()
    verb = verdict["verb"].strip()
    cells: list[VerbCell] = []
    for cell, column, endings in (
        ("2sg", "pres_2sg", mo.VERB_ENDINGS_2SG),
        ("subj3sg", "subj_3sg", mo.VERB_ENDINGS_SUBJ3SG),
    ):
        form = mo.first_variant(verdict.get(column, ""), verb)
        if not form:
            continue
        reading = mo.read_verb_cell(form, endings, seg, cluster)
        if reading is None:
            continue
        cells.append(
            VerbCell(
                verb=verb,
                cell=cell,
                form=reading.form,
                vowel=reading.vowel,
                palatal=reading.palatal,
                augment=reading.augment,
                count=counts.count(reading.form),
            )
        )
    return verb, verdict["verdict"].strip(), tuple(cells)


def _derivatives(
    verbs: Sequence[tuple[str, str]],
    adjectives: Sequence[str],
    seg: str,
    cluster: str,
    counts: RomanianCounts,
) -> tuple[Derivative, ...]:
    out: list[Derivative] = []
    for verb, suffix in verbs:
        reading = mo.read_derived_verb(verb, suffix, seg, cluster)
        count = counts.count(verb)
        part = mo.participle(verb, suffix)
        out.append(
            Derivative(
                kind="verb",
                form=verb.lower(),
                suffix=suffix,
                vowel=reading.vowel if reading else None,
                palatal=reading.palatal if reading else None,
                count=count,
                attest_count=count + (counts.count(part) if part else 0),
            )
        )
    for adjective in adjectives:
        reading = mo.read_derived_adjective(adjective, seg, cluster)
        count = counts.count(adjective)
        out.append(
            Derivative(
                kind="adj",
                form=adjective.lower(),
                suffix=reading.suffix if reading else "",
                vowel=reading.vowel if reading else None,
                palatal=reading.palatal if reading else None,
                count=count,
                attest_count=count,
            )
        )
    return tuple(out)


def _candidates(
    row: lx.LexiconRow, seg: str, counts: RomanianCounts
) -> tuple[Candidate, ...]:
    pair = mo.candidate_verbs(row.lemma, seg, row.cluster) if seg else None
    if pair is None:
        return ()
    return (
        Candidate(pair[0], True, counts.count(pair[0])),
        Candidate(pair[1], False, counts.count(pair[1])),
    )


def _merged_pairs(group: Sequence[lx.LexiconRow]) -> list[tuple[str, str]]:
    """(verb, suffix) pairs of all rows, priority order, no repeats.

    ``derived_verbs`` and ``deriv_suffixes`` are parallel lists; a verb
    without a suffix entry gets "".
    """
    seen: dict[tuple[str, str], None] = {}
    for row in group:
        suffixes = list(row.deriv_suffixes)
        suffixes += [""] * (len(row.derived_verbs) - len(suffixes))
        for verb, suffix in zip(row.derived_verbs, suffixes, strict=False):
            seen.setdefault((verb, suffix), None)
    return list(seen)


def build_record(
    group: Sequence[lx.LexiconRow],
    verdicts: Mapping[str, Mapping[str, str]],
    counts: RomanianCounts,
) -> RootRecord:
    """Merge one lemma's noun rows (priority order) into a record."""
    primary = group[0]
    seg = primary.stem_final if primary.stem_final in lx.TARGETS else ""
    nde = lx.nde_class(primary)
    target_class = lx.TARGET_CLASS.get(seg, "")
    if nde == "frontstem":
        target_class = ""  # lexical palatal: not dorsal-final
    verb_pairs = _merged_pairs(group)
    adjectives = lx.merged_unique(row.derived_adj for row in group)
    verb, audit, cells = _verb_cells(
        verdicts.get(primary.lemma), seg, primary.cluster, counts
    )
    allomorph = lx.plural_allomorph(primary, nde)
    palatal = primary.mutation if allomorph in ("-i", "-e") else None
    return RootRecord(
        root_id=lx.ROOT_PREFIX + primary.lemma,
        lemma=primary.lemma,
        gloss=primary.gloss,
        gender=primary.gender,
        seg=seg,
        cluster=primary.cluster,
        target_class=target_class,
        nde=nde,
        flags=lx.flags(primary) if seg else (),
        stratum=lx.stratum(primary.etym_lang),
        lexicon_source=lx.lexicon_source(primary.source),
        plural=primary.plural,
        plural_allomorph=allomorph,
        plural_palatal=palatal,
        lemma_count=counts.count(primary.lemma),
        plural_count=counts.count(primary.plural) if primary.plural else 0,
        paradigm_count=paradigm_frequency(
            [lx.fold_key(primary.lemma), lx.fold_key(primary.plural)],
            counts.counts,
        ),
        verb=verb,
        audit_verdict=audit,
        verb_cells=cells,
        derivatives=_derivatives(
            verb_pairs, adjectives, seg, primary.cluster, counts
        ),
        candidates=_candidates(primary, seg, counts),
        legacy=Legacy(
            tp_in_domain=primary.tp_in_domain,
            exception_reason=primary.exception_reason,
            opportunity=primary.opportunity,
            mutation=primary.mutation,
            orth_change=primary.orth_change,
            target_is_suffix=primary.target_is_suffix,
            lemma_suffix=primary.lemma_suffix,
        ),
        duplicates=len(group),
    )


def build_records(
    rows: Iterable[lx.LexiconRow],
    verdicts: Mapping[str, Mapping[str, str]],
    counts: RomanianCounts,
) -> list[RootRecord]:
    """One record per noun lemma, sorted by root_id."""
    groups = lx.group_nouns(rows)
    records = [
        build_record(group, verdicts, counts) for group in groups.values()
    ]
    return sorted(records, key=lambda record: record.root_id)


# ---------------------------------------------------------------------------
# JSONL round trip
# ---------------------------------------------------------------------------


def _dump(obj: Any) -> str:
    return json.dumps(
        obj, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    )


def record_to_dict(record: RootRecord) -> dict[str, Any]:
    data = asdict(record)
    data["flags"] = list(record.flags)
    return data


def record_from_dict(data: Mapping[str, Any]) -> RootRecord:
    return RootRecord(
        root_id=data["root_id"],
        lemma=data["lemma"],
        gloss=data["gloss"],
        gender=data["gender"],
        seg=data["seg"],
        cluster=data["cluster"],
        target_class=data["target_class"],
        nde=data["nde"],
        flags=tuple(data["flags"]),
        stratum=data["stratum"],
        lexicon_source=data["lexicon_source"],
        plural=data["plural"],
        plural_allomorph=data["plural_allomorph"],
        plural_palatal=data["plural_palatal"],
        lemma_count=data["lemma_count"],
        plural_count=data["plural_count"],
        paradigm_count=data["paradigm_count"],
        verb=data["verb"],
        audit_verdict=data["audit_verdict"],
        verb_cells=tuple(VerbCell(**cell) for cell in data["verb_cells"]),
        derivatives=tuple(Derivative(**d) for d in data["derivatives"]),
        candidates=tuple(Candidate(**c) for c in data["candidates"]),
        legacy=Legacy(**data["legacy"]),
        duplicates=data["duplicates"],
    )


@dataclass(frozen=True, slots=True)
class RootsFile:
    """A parsed roots JSONL: header fields and the records."""

    total_tokens: int
    inputs: Mapping[str, str]
    records: tuple[RootRecord, ...]


def roots_lines(
    records: Sequence[RootRecord], total_tokens: int, inputs: Mapping[str, str]
) -> Iterator[str]:
    """The JSONL lines (header first), each ending in a newline."""
    header = {
        "schema": SCHEMA,
        "total_tokens": total_tokens,
        "inputs": dict(sorted(inputs.items())),
        "roots": len(records),
    }
    yield _dump(header) + "\n"
    for record in sorted(records, key=lambda r: r.root_id):
        yield _dump(record_to_dict(record)) + "\n"


def read_roots(path: Path) -> RootsFile:
    """Read a file written from :func:`roots_lines`."""
    with path.open("r", encoding="utf-8") as fh:
        header = json.loads(fh.readline())
        if header.get("schema") != SCHEMA:
            raise ValueError(f"{path}: not a {SCHEMA} file")
        records = tuple(
            record_from_dict(json.loads(line)) for line in fh if line.strip()
        )
    if len(records) != header["roots"]:
        raise ValueError(
            f"{path}: {len(records)} roots, header says {header['roots']}"
        )
    return RootsFile(
        total_tokens=int(header["total_tokens"]),
        inputs=dict(header["inputs"]),
        records=records,
    )
