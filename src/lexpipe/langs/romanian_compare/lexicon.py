"""Romanian lexicon rows and outcome-blind root attributes (design §7).

The harmonized Romanian comparison reads the Romanian pipeline's
output ``data/romanian_lexicon_with_freq.csv`` (read-only, R1) and
``data/verb_paradigm_verdict.csv`` (E1, R15). It never imports the
Romanian pipeline (``romanian_processor_lib``, ``wiktionary_normalizer``,
``dex_utils``, ``romanian_harvester``; R1): the few orthographic facts
it needs are restated here and documented where they mirror the
legacy code.

Roots
-----
One root per noun lemma (``pos == "N"``), ``root_id = "ro:<lemma>"``
(R4). The lexicon has one plural per lemma; four lemmas occur in two
rows (an en.wiktionary row without a plural and a ro.wiktionary row
with one). Rows of one lemma are merged by an explicit priority (a
row with a plural first, then by ``source``, then by the whole row),
never first-wins (R5); list fields are merged in that order without
repeats. Adjective rows are not pooled with nouns: the legacy analysis
and this comparison count noun roots only.

Outcome-blind attributes (R16)
------------------------------
The legacy TP domain (``mark_tp_domain``) drops some items only when
they did *not* palatalize: the "foreign plural" rule (111 z→zi nouns
such as albanez/albanezi), the shared-prefix "suppletion" rule
(fată/fete, pată/pete) and the already-palatal check of
``derive_exception_reason``. Here every attribute is computed from
the lemma (or the plural's shape) for every root, whatever its
outcome:

``nde``
    structural non-derived-environment class, from the lemma only:
    ``ochi`` (dorsal + chi/ghi), ``paduchi`` (dorsal + che/ghe),
    ``frontstem`` (lexical palatal: the legacy ``frontstem_dorsal``
    flag, ci/ce/gi/ge at the stem edge; not dorsal-final, so never in
    a dorsal class), ``coronal_frontstem`` (t/d/s/z followed in the
    lemma by a root-internal front vowel: -tie, -zie, -diu, -zi, -ea).
``flags``
    ``suffix_target`` (legacy ``target_is_suffix``: the target sits in
    -ic/-ică/-ist/...), ``palatal_context`` (the legacy already-palatal
    test applied to every root: the stem tail already holds the
    target's palatal letter, or ș/ț/j precedes it), ``suppletive_stem``
    (lemma and plural share too short a prefix once vowels are treated
    as equal and the target may differ only by its palatal reflex), and
    ``ica_suppletion`` (-ică whose plural replaces the suffix, -ele).
``stratum``
    R22 from ``etym_lang``: ``borrowed`` (any etymon language other
    than Latin), ``latin_unresolved`` (``la``: inherited or a learned
    borrowing; the lexicon keeps no template that separates the two)
    and ``unknown`` (empty). "internal formation" cannot be told
    from this column.

The plural cell's allomorph comes from the legacy ``opportunity``
column (``i`` → -i, ``e`` → -e, ``uri`` → -uri; a plural the legacy
alignment could not parse is ``unparsed``). For the structural NDE
classes the R script's ``opportunity_tp`` patch is applied (NDE with
``opportunity == "none"``: the plural's ending, -uri counting as -i).
"""

from __future__ import annotations

import csv
import re
import unicodedata
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from lexpipe import csvio

ROOT_PREFIX: Final[str] = "ro:"
VOWELS: Final[frozenset[str]] = frozenset("aăâeiîou")
FRONT_VOWELS: Final[frozenset[str]] = frozenset("ie")
TARGETS: Final[tuple[str, ...]] = ("c", "g", "t", "d", "s", "z")
TARGET_CLASS: Final[Mapping[str, str]] = {
    "c": "dorsal",
    "g": "dorsal",
    "t": "coronal",
    "d": "coronal",
    "s": "coronal",
    "z": "z",
}
# Letter of each target's palatal reflex (legacy PALATAL_SURFACE plus
# the dorsals, whose reflex is spelled with the same letter before a
# front vowel).
PALATAL_LETTER: Final[Mapping[str, str]] = {
    "c": "c",
    "g": "g",
    "t": "ț",
    "d": "z",
    "s": "ș",
    "z": "j",
}
# Romanian corpora mix cedilla and comma-below letters; the lexicon
# uses comma-below (scripts/attach_frequencies.py folds the same way).
CEDILLA_TO_COMMA: Final[Mapping[int, str]] = str.maketrans(
    {"ş": "ș", "Ş": "Ș", "ţ": "ț", "Ţ": "Ț"}
)
NDE_CLASSES: Final[tuple[str, ...]] = (
    "ochi",
    "paduchi",
    "frontstem",
    "coronal_frontstem",
)
FLAG_NAMES: Final[tuple[str, ...]] = (
    "ica_suppletion",
    "palatal_context",
    "suffix_target",
    "suppletive_stem",
)
STRATA: Final[tuple[str, ...]] = ("borrowed", "latin_unresolved", "unknown")

REQUIRED_COLUMNS: Final[tuple[str, ...]] = (
    "lemma",
    "gloss",
    "pos",
    "gender",
    "stem_final",
    "cluster",
    "frontstem_dorsal",
    "plural",
    "mutation",
    "orth_change",
    "opportunity",
    "derived_verbs",
    "deriv_suffixes",
    "derived_adj",
    "ipa_derived_adj",
    "etym_lang",
    "exception_reason",
    "tp_in_domain",
    "target_is_suffix",
    "lemma_suffix",
    "source",
    "notes",
    "freq_ron_news_2024_1M",
)
FREQ_COLUMN: Final[str] = "freq_ron_news_2024_1M"


def fold_key(word: str) -> str:
    """Corpus lookup key: NFC, cedilla → comma-below, stripped, lower."""
    text = unicodedata.normalize("NFC", word).translate(CEDILLA_TO_COMMA)
    return text.strip().lower()


def strip_final_vowel(word: str) -> str:
    """Drop one final vowel (legacy ``_strip_final_vowel``)."""
    if len(word) > 1 and word[-1] in VOWELS:
        return word[:-1]
    return word


def split_pipe(cell: str) -> tuple[str, ...]:
    """``a|b`` → ``("a", "b")``, items stripped, empty items dropped."""
    return tuple(item.strip() for item in cell.split("|") if item.strip())


@dataclass(frozen=True, slots=True)
class LexiconRow:
    """The columns of one lexicon row that the comparison uses."""

    lemma: str
    gloss: str
    pos: str
    gender: str
    stem_final: str
    cluster: str
    frontstem_dorsal: bool
    plural: str
    mutation: bool | None
    orth_change: str
    opportunity: str
    derived_verbs: tuple[str, ...]
    deriv_suffixes: tuple[str, ...]
    derived_adj: tuple[str, ...]
    ipa_derived_adj: str
    etym_lang: str
    exception_reason: str
    tp_in_domain: bool
    target_is_suffix: bool
    lemma_suffix: str
    source: str
    notes: str
    freq: int


def _bool(cell: str) -> bool | None:
    return csvio.parse_bool(cell)


def _int(cell: str) -> int:
    text = cell.strip()
    return int(text) if text else 0


def parse_row(row: Mapping[str, str]) -> LexiconRow:
    """One CSV row (all cells strings) → :class:`LexiconRow`."""
    return LexiconRow(
        lemma=row["lemma"].strip(),
        gloss=row["gloss"].strip(),
        pos=row["pos"].strip().upper(),
        gender=row["gender"].strip(),
        stem_final=row["stem_final"].strip().lower(),
        cluster=row["cluster"].strip().lower(),
        frontstem_dorsal=_bool(row["frontstem_dorsal"]) is True,
        plural=row["plural"].strip(),
        mutation=_bool(row["mutation"]),
        orth_change=row["orth_change"].strip(),
        opportunity=row["opportunity"].strip(),
        derived_verbs=split_pipe(row["derived_verbs"]),
        deriv_suffixes=(
            tuple(s.strip() for s in row["deriv_suffixes"].split("|"))
            if row["deriv_suffixes"].strip()
            else ()
        ),
        derived_adj=split_pipe(row["derived_adj"]),
        ipa_derived_adj=row["ipa_derived_adj"].strip(),
        etym_lang=row["etym_lang"].strip(),
        exception_reason=row["exception_reason"].strip(),
        tp_in_domain=_bool(row["tp_in_domain"]) is True,
        target_is_suffix=_bool(row["target_is_suffix"]) is True,
        lemma_suffix=row["lemma_suffix"].strip(),
        source=row["source"].strip(),
        notes=row["notes"].strip(),
        freq=_int(row[FREQ_COLUMN]),
    )


def read_lexicon(path: Path) -> list[LexiconRow]:
    """All rows of the lexicon CSV, in file order."""
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        missing = [
            c for c in REQUIRED_COLUMNS if c not in (reader.fieldnames or ())
        ]
        if missing:
            raise ValueError(f"{path}: missing columns {missing}")
        return [parse_row(row) for row in reader]


# ---------------------------------------------------------------------------
# Outcome-blind attributes
# ---------------------------------------------------------------------------


def target_position(lemma: str, seg: str) -> int:
    """Index in *lemma* of its stem-final target *seg*, or -1.

    Mirrors the legacy ``derive_stem_final_and_cluster``: the rightmost
    occurrence of *seg* followed only by vowels (one final vowel
    stripped first, so the vowels after it may include the stem's).
    """
    if not seg:
        return -1
    for index in range(len(lemma) - 1, -1, -1):
        if lemma[index] == seg and all(
            ch in VOWELS for ch in lemma[index + 1 :]
        ):
            return index
    return -1


def nde_class(row: LexiconRow) -> str:
    """Structural NDE class of a noun (see the module docstring)."""
    seg, cluster = row.stem_final, row.cluster
    if seg in ("c", "g") and cluster in ("chi", "ghi"):
        return "ochi"
    if seg in ("c", "g") and cluster in ("che", "ghe"):
        return "paduchi"
    if row.frontstem_dorsal:
        return "frontstem"
    if seg in ("t", "d", "s", "z") and not cluster:
        lemma = row.lemma.lower()
        index = target_position(lemma, seg)
        tail = lemma[index + 1 :] if index >= 0 else ""
        if tail.startswith("i") or (tail.startswith("e") and len(tail) >= 2):
            return "coronal_frontstem"
    return ""


def palatal_context(row: LexiconRow) -> bool:
    """The legacy already-palatal test (``derive_exception_reason``),
    applied to every root, mutated or not."""
    seg = row.stem_final
    if not seg:
        return False
    stem = strip_final_vowel(row.lemma.lower())
    tail = stem[-4:]
    if seg == "d":
        if re.search(r"z[ie]", tail):
            return True
    elif seg in ("t", "s", "z") and PALATAL_LETTER[seg] in tail:
        return True
    if len(stem) >= 2:
        position = stem.rfind(seg)
        if position > 0 and stem[position - 1] in ("ș", "ț", "j"):
            return True
    return False


_REFLEX_EQUIVALENT: Final[frozenset[tuple[str, str]]] = frozenset(
    {("t", "ț"), ("s", "ș"), ("z", "j"), ("d", "z")}
)


def _letters_match(a: str, b: str) -> bool:
    if a == b or (a in VOWELS and b in VOWELS):
        return True
    return (a, b) in _REFLEX_EQUIVALENT or (b, a) in _REFLEX_EQUIVALENT


def suppletive_stem(lemma: str, plural: str) -> bool:
    """Too short a shared prefix for an inflected form of the lemma.

    The legacy threshold (``max(2, 2·len/5)`` letters) with every vowel
    matching every vowel and a target matching its palatal reflex, so a
    vowel alternation (fată/fete, carte/cărți) or the palatalization
    itself never counts as suppletion: outcome-blind.
    """
    if not plural:
        return False
    a, b = lemma.lower(), plural.lower()
    shared = 0
    for x, y in zip(a, b, strict=False):
        if not _letters_match(x, y):
            break
        shared += 1
    return shared < max(2, len(a) * 2 // 5)


def ica_suppletion(lemma: str, plural: str) -> bool:
    """-ică whose plural replaces the suffix (rândunică/rândunele)."""
    lowered, plural_l = lemma.lower(), plural.lower()
    return (
        bool(plural_l)
        and lowered.endswith("ică")
        and not plural_l.endswith(("ici", "ice", "ică"))
    )


def stratum(etym_lang: str) -> str:
    """R22 stratum from the etymon language code (module docstring)."""
    code = etym_lang.strip().lower()
    if not code:
        return "unknown"
    if code == "la":
        return "latin_unresolved"
    return "borrowed"


def flags(row: LexiconRow) -> tuple[str, ...]:
    """Sorted outcome-blind flags of a noun row."""
    found: list[str] = []
    if row.target_is_suffix:
        found.append("suffix_target")
    if palatal_context(row):
        found.append("palatal_context")
    if suppletive_stem(row.lemma, row.plural):
        found.append("suppletive_stem")
    if ica_suppletion(row.lemma, row.plural):
        found.append("ica_suppletion")
    return tuple(sorted(found))


def plural_ending(plural: str) -> str:
    """The R script's ``plural_ending``: uri / i / e / other / none."""
    text = plural.strip()
    if not text:
        return "none"
    if text.endswith("uri"):
        return "uri"
    if text.endswith("i"):
        return "i"
    if text.endswith("e"):
        return "e"
    return "other"


def plural_allomorph(row: LexiconRow, nde: str) -> str:
    """-i / -e / -uri / unparsed / none (see the module docstring)."""
    if not row.plural:
        return "none"
    opportunity = row.opportunity
    if nde in ("ochi", "paduchi", "frontstem", "coronal_frontstem"):
        if opportunity == "none":
            ending = plural_ending(row.plural)
            opportunity = {"i": "i", "e": "e", "uri": "i"}.get(ending, "none")
    return {"i": "-i", "e": "-e", "uri": "-uri"}.get(opportunity, "unparsed")


def lexicon_source(source: str) -> str:
    """Which Wiktionary edition (and DEX QC) a row came from."""
    parts = [part.strip() for part in source.split("|") if part.strip()]
    labels: list[str] = []
    for part in parts:
        if "en.wiktionary" in part:
            labels.append("en.wiktionary")
        elif "ro.wiktionary" in part:
            labels.append("ro.wiktionary")
        elif part == "dex_qc":
            labels.append("dex_qc")
        else:
            labels.append("other")
    return "+".join(labels)


def _priority(row: LexiconRow) -> tuple[int, str, str]:
    return (0 if row.plural else 1, row.source, repr(row))


def group_nouns(rows: Iterable[LexiconRow]) -> dict[str, list[LexiconRow]]:
    """Noun rows by lemma, each group in merge-priority order."""
    groups: dict[str, list[LexiconRow]] = {}
    for row in rows:
        if row.pos == "N" and row.lemma:
            groups.setdefault(row.lemma, []).append(row)
    return {
        lemma: sorted(group, key=_priority)
        for lemma, group in sorted(groups.items())
    }


def merged_unique(values: Iterable[Sequence[str]]) -> tuple[str, ...]:
    """Concatenate sequences, keeping the first of each repeated item."""
    seen: dict[str, None] = {}
    for sequence in values:
        for item in sequence:
            seen.setdefault(item, None)
    return tuple(seen)
