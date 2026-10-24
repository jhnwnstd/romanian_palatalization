"""MorfFlex CZ 2.1: lemma grammar, tags and the compact projection (R7, R14).

MorfFlex is a flat ``lemma<TAB>tag<TAB>form`` list of 126.9 million
rows in one xz block (CC BY-NC-SA 4.0). Everything derived from it is
NC-licensed and stays under ``data/czech/nc/`` (R9).

Lemmas (TR-64, Mikulová et al. 2020)
------------------------------------
A raw lemma is ``<base>[-<homonym>]`` followed by optional segments,
each introduced by ``_``::

    _:T / _:W      aspect (verbs)
    _;X            label: name labels Y (given name), E (ethnic group),
                   G (geographic), U, m (other proper name), o; other
                   letters are semantic labels
    _,X            style label (s spelling variant, h colloquial,
                   a archaic, l, i, e, n, v, t …)
    _^(…)          explanatory comment; ``_^(*3ka)`` is derivation info
    _(…)           derivation info (``_(*2t)``: cut 2 letters, add t)

:func:`parse_pdt_lemma` splits these. Bases may contain digits (Aa1)
and hyphens (Coca-Cola): only ``-<digits>`` directly before the first
segment is a homonym number. The raw lemma stays the join key
(Čech_;E_;Y merges the ethnonym and the surname, which Wiktionary keeps
apart); the base is what is compared with a Wiktionary lemma.

Tags
----
Positions (0-based): 0 POS, 1 SubPOS, 2 gender (M animate, I
inanimate, F, N, or a cover letter), 3 number, 4 case (1–7), 9 degree,
10 negation, 14 variant (Var). Var ``-`` and ``1``–``4`` are standard,
``5``–``9`` non-standard, ``a``/``b``/``c`` abbreviations (dropped).
Variant digits number variants, they do not rank them (a ``1`` can be
archaic), so no order is read into them.

Projection
----------
:func:`iter_rows` streams the xz file in large chunks and lets one byte
regular expression pick the rows R7 keeps: affirmative common-noun rows
(tag ``NN????????A????``) and the hard-adjective cell
``AAMP1----1A----`` (NOM.PL masculine animate, positive, affirmative).
:func:`build_db` stores those whose base is a Wiktionary lemma or whose
form (case-folded) is in the Leipzig vocabulary, in the SQLite table
``morfflex`` (columns :data:`DB_COLUMNS`), stamped in ``meta`` with the
source SHA-256.
"""

from __future__ import annotations

import lzma
import re
from collections.abc import Collection, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Final

from lexpipe.context import RunContext
from lexpipe.langs.czech import jsonio

XZ_MAGIC: Final[bytes] = b"\xfd7zXZ\x00"
CHUNK: Final[int] = 1 << 24
ADJECTIVE_TAG: Final[str] = "AAMP1----1A----"
ROW_RE: Final[re.Pattern[bytes]] = re.compile(
    rb"^([^\t\n]+)\t(NN[^\t\n]{8}A[^\t\n]{4}|AAMP1----1A----)\t([^\t\n]+)$",
    re.MULTILINE,
)
NAME_LABELS: Final[frozenset[str]] = frozenset("YEGUmo")
VAR_STANDARD: Final[frozenset[str]] = frozenset("-1234")
VAR_NONSTANDARD: Final[frozenset[str]] = frozenset("56789")
VAR_ABBREVIATION: Final[frozenset[str]] = frozenset("abc")
STANDARD: Final[str] = "standard"
NONSTANDARD: Final[str] = "nonstandard"
ABBREVIATION: Final[str] = "abbreviation"
CASE_CODES: Final[dict[str, str]] = {
    "1": "NOM",
    "2": "GEN",
    "3": "DAT",
    "4": "ACC",
    "5": "VOC",
    "6": "LOC",
    "7": "INS",
}
NUMBER_CODES: Final[dict[str, str]] = {"S": "SG", "P": "PL"}
TAG_LENGTH: Final[int] = 15
DB_COLUMNS: Final[tuple[tuple[str, str], ...]] = (
    ("raw_lemma", "TEXT"),
    ("base", "TEXT"),
    ("homonym", "INTEGER"),
    ("name_labels", "TEXT"),
    ("style_labels", "TEXT"),
    ("other_labels", "TEXT"),
    ("comment", "TEXT"),
    ("tag", "TEXT"),
    ("cell", "TEXT"),
    ("gender", "TEXT"),
    ("var", "TEXT"),
    ("var_class", "TEXT"),
    ("form", "TEXT"),
    ("form_folded", "TEXT"),
    ("in_lexicon", "INTEGER"),
    ("in_vocab", "INTEGER"),
)
_LEMMA_RE: Final[re.Pattern[str]] = re.compile(
    r"^(?P<base>.+?)(?:-(?P<hom>\d+))?(?P<rest>(?:_[:;,^(].*)?)$", re.DOTALL
)


@dataclass(frozen=True, slots=True)
class PdtLemma:
    """A parsed MorfFlex/PDT lemma (see the module docstring)."""

    raw: str
    base: str
    homonym: int | None = None
    name_labels: tuple[str, ...] = ()
    style_labels: tuple[str, ...] = ()
    other_labels: tuple[str, ...] = ()
    aspect: str = ""
    comment: str = ""
    derivation: str = ""
    unparsed: str = ""

    @property
    def is_proper(self) -> bool:
        return bool(self.name_labels)

    @property
    def is_styled(self) -> bool:
        return bool(self.style_labels)


def _balanced(text: str, start: int) -> tuple[str, int]:
    """Contents of the parenthesized group opening at *start*, and the
    index after its closing parenthesis (end of text if unbalanced)."""
    depth = 0
    for index in range(start, len(text)):
        if text[index] == "(":
            depth += 1
        elif text[index] == ")":
            depth -= 1
            if depth == 0:
                return text[start + 1 : index], index + 1
    return text[start + 1 :], len(text)


def _label(text: str, start: int) -> tuple[str, int]:
    """A label value: everything up to the next ``_`` segment start."""
    end = start
    while end < len(text) and not (
        text[end] == "_" and end + 1 < len(text) and text[end + 1] in ":;,^("
    ):
        end += 1
    return text[start:end], end


def parse_pdt_lemma(raw: str) -> PdtLemma:
    """Split a raw MorfFlex lemma (see the module docstring)."""
    match = _LEMMA_RE.match(raw)
    if match is None:  # empty string
        return PdtLemma(raw=raw, base=raw)
    rest = match.group("rest")
    names: list[str] = []
    styles: list[str] = []
    others: list[str] = []
    aspect = comment = derivation = ""
    index = 0
    while index < len(rest):
        if rest[index] != "_" or index + 1 >= len(rest):
            break
        kind = rest[index + 1]
        if kind in ":;,":
            value, index = _label(rest, index + 2)
            if kind == ":":
                aspect = value
            elif kind == ",":
                styles.append(value)
            elif value in NAME_LABELS:
                names.append(value)
            else:
                others.append(value)
        elif kind == "^" and index + 2 < len(rest) and rest[index + 2] == "(":
            comment, index = _balanced(rest, index + 2)
            if comment.startswith("*"):
                derivation = comment
        elif kind == "(":
            derivation, index = _balanced(rest, index + 1)
        else:
            break
    hom = match.group("hom")
    return PdtLemma(
        raw=raw,
        base=match.group("base"),
        homonym=int(hom) if hom is not None else None,
        name_labels=tuple(names),
        style_labels=tuple(styles),
        other_labels=tuple(others),
        aspect=aspect,
        comment=comment,
        derivation=derivation,
        unparsed=rest[index:],
    )


@dataclass(frozen=True, slots=True)
class PdtTag:
    """The tag positions the plugin uses."""

    raw: str
    pos: str
    subpos: str
    gender: str
    number: str
    case: str
    degree: str
    negation: str
    var: str

    @property
    def cell(self) -> str | None:
        """``NOM.SG`` … for noun tags, ``NOM.PL.M.AN`` for the adjective
        cell; None for other numbers (X, D) or cases."""
        if self.raw == ADJECTIVE_TAG:
            return "NOM.PL.M.AN"
        case = CASE_CODES.get(self.case)
        number = NUMBER_CODES.get(self.number)
        if case is None or number is None:
            return None
        return f"{case}.{number}"


def parse_tag(tag: str) -> PdtTag:
    """Read a 15-position MorfFlex tag; ValueError for another length."""
    if len(tag) != TAG_LENGTH:
        raise ValueError(f"MorfFlex tag of length {len(tag)}: {tag!r}")
    return PdtTag(
        raw=tag,
        pos=tag[0],
        subpos=tag[1],
        gender=tag[2],
        number=tag[3],
        case=tag[4],
        degree=tag[9],
        negation=tag[10],
        var=tag[14],
    )


def var_class(var: str) -> str:
    """standard / nonstandard / abbreviation (TR-64 position 15)."""
    if var in VAR_STANDARD:
        return STANDARD
    if var in VAR_NONSTANDARD:
        return NONSTANDARD
    if var in VAR_ABBREVIATION:
        return ABBREVIATION
    return NONSTANDARD


def _open(path: Path) -> IO[bytes]:
    with path.open("rb") as probe:
        magic = probe.read(len(XZ_MAGIC))
    if magic == XZ_MAGIC:
        return lzma.open(path, "rb")
    return path.open("rb")


def iter_rows(
    path: Path, chunk_size: int = CHUNK
) -> Iterator[tuple[str, str, str]]:
    """``(lemma, tag, form)`` of the kept rows, in file order."""
    tail = b""
    with _open(path) as fh:
        while True:
            chunk = fh.read(chunk_size)
            data = tail + chunk
            if chunk:
                cut = data.rfind(b"\n")
                if cut < 0:
                    tail = data
                    continue
                tail, data = data[cut + 1 :], data[: cut + 1]
            else:
                tail = b""
            for match in ROW_RE.finditer(data):
                yield (
                    match.group(1).decode("utf-8"),
                    match.group(2).decode("ascii"),
                    match.group(3).decode("utf-8"),
                )
            if not chunk:
                return


@dataclass(frozen=True, slots=True)
class BuildStats:
    """Rows seen by the filter and rows stored."""

    matched: int
    stored: int


def _db_rows(
    path: Path,
    lemmas: Collection[str],
    vocabulary: Collection[str],
    counter: list[int],
) -> Iterator[tuple[object, ...]]:
    cache: dict[str, PdtLemma] = {}
    for raw, tag_text, form in iter_rows(path):
        counter[0] += 1
        lemma = cache.get(raw)
        if lemma is None:
            if len(cache) > 100_000:
                cache.clear()
            lemma = cache[raw] = parse_pdt_lemma(raw)
        folded = form.lower()
        in_lexicon = lemma.base in lemmas
        in_vocab = folded in vocabulary
        if not (in_lexicon or in_vocab):
            continue
        tag = parse_tag(tag_text)
        yield (
            raw,
            lemma.base,
            lemma.homonym,
            "|".join(lemma.name_labels),
            "|".join(lemma.style_labels),
            "|".join(lemma.other_labels),
            lemma.comment,
            tag_text,
            tag.cell,
            tag.gender,
            tag.var,
            var_class(tag.var),
            form,
            folded,
            int(in_lexicon),
            int(in_vocab),
        )


def build_db(
    ctx: RunContext,
    source: Path,
    target: Path,
    *,
    lemmas: Collection[str],
    vocabulary: Collection[str],
    meta: tuple[tuple[str, str], ...] = (),
) -> BuildStats:
    """One streaming pass over MorfFlex into *target* (module docstring)."""
    counter = [0]
    stored = jsonio.build_sqlite(
        ctx,
        target,
        table="morfflex",
        columns=DB_COLUMNS,
        rows=_db_rows(source, lemmas, vocabulary, counter),
        order_by=("base", "raw_lemma", "tag", "form"),
        indexes=(
            ("morfflex_base", ("base",)),
            ("morfflex_form", ("form_folded",)),
        ),
        meta=meta,
    )
    return BuildStats(matched=counter[0], stored=stored)
