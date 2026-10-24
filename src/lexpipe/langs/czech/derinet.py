"""DeriNet 2.3: one streaming pass into a compact NC derivative (R7, R14).

DeriNet (Vidra et al. 2019 format; DeriNet 2.3 by Olbrich et al. 2025;
CC BY-NC-SA 4.0) is a 439 MB TSV of derivation trees, one node per
line and a blank line between trees, ten columns::

    ID  LEMID  LEMMA  POS  FEATS  SEGMENTATION  PARENT_ID  RELATION
    OTHER_RELATIONS  OTHER

FEATS is ``&``-separated ``key=value`` (Gender, Animacy, Loanword,
ConjugClass, Aspect, Fictitious …); SEGMENTATION and OTHER are JSON
(morphs typed R/D/I; corpus_stats, techlemma, Etymology, unmotivated).
Lines are parsed, never evaluated; a malformed JSON cell raises
:class:`ValueError` naming the line.

:func:`project` reads the file once, tree by tree, and keeps

``nouns``   every NOUN node whose lemma is a noun lemma of the Czech
            lexicon: its gender and animacy (the R14 join key: šach
            and jazyk are each two nodes), Loanword flag, whether it
            has a derivational parent, DeriNet's Etymology codes and
            corpus count;
``edges``   every ``Type=Derivation`` edge whose parent is such a noun
            (R14: noun-parented derivation edges only), with the
            child's lemma, part of speech, features, segmentation and
            the **verbalizer** read from the last derivational (D)
            morph of the segmentation (knih|ova|t → -ovat;
            s|nož|i|t → -it). ConjugClass, the traditional present-tense
            class (4 holds -it and -et/-ět verbs, 3 mostly -ovat), is
            kept only to check the verbalizer (:func:`conj_consistent`).

Fictitious nodes (affixoid roots such as ``-acid-``) are dropped with
their edges. The derivative is a SQLite file under ``data/czech/nc/``
written by :func:`lexpipe.langs.czech.jsonio.build_sqlite` (fixed row
order, so it is byte-deterministic) with DeriNet's SHA-256 in its
``meta`` table.
"""

from __future__ import annotations

import contextlib
import json
import sqlite3
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

from lexpipe.context import RunContext, sha256_file
from lexpipe.langs.czech import jsonio

UPOS: Final[Mapping[str, str]] = {
    "NOUN": "noun",
    "ADJ": "adj",
    "VERB": "verb",
    "ADV": "adv",
}
GENDERS: Final[Mapping[str, str]] = {"Masc": "m", "Fem": "f", "Neut": "n"}
ANIMACIES: Final[Mapping[str, str]] = {"Anim": "an", "Inan": "in"}
DERIVATION: Final[str] = "Derivation"
N_COLUMNS: Final[int] = 10

# Last D morph of a verb → verbalizer label (R14).
VERBALIZERS: Final[Mapping[str, str]] = {
    "ova": "-ovat",
    "i": "-it",
    "e": "-et",
    "ě": "-et",
    "nou": "-nout",
    "a": "-at",
    "á": "-at",
}
# ConjugClass values a verbalizer is expected to have.
EXPECTED_CLASS: Final[Mapping[str, str]] = {
    "-it": "4",
    "-et": "4",
    "-ovat": "3",
}

NOUN_COLUMNS: Final[tuple[tuple[str, str], ...]] = (
    ("node_id", "TEXT"),
    ("lemma", "TEXT"),
    ("techlemma", "TEXT"),
    ("gender", "TEXT"),
    ("animacy", "TEXT"),
    ("loanword", "TEXT"),
    ("has_parent", "INTEGER"),
    ("etymology", "TEXT"),
    ("corpus_count", "INTEGER"),
)
EDGE_COLUMNS: Final[tuple[tuple[str, str], ...]] = (
    ("parent_id", "TEXT"),
    ("parent_lemma", "TEXT"),
    ("child_id", "TEXT"),
    ("child_lemma", "TEXT"),
    ("child_techlemma", "TEXT"),
    ("child_pos", "TEXT"),
    ("child_gender", "TEXT"),
    ("child_animacy", "TEXT"),
    ("child_loanword", "TEXT"),
    ("conj_class", "TEXT"),
    ("aspect", "TEXT"),
    ("segmentation", "TEXT"),
    ("last_d_morph", "TEXT"),
    ("verbalizer", "TEXT"),
    ("semantic_label", "TEXT"),
    ("corpus_count", "INTEGER"),
)


@dataclass(frozen=True, slots=True)
class Node:
    """One parsed DeriNet line."""

    node_id: str
    lemma: str
    upos: str
    feats: Mapping[str, str]
    morphs: tuple[tuple[str, str], ...]
    parent_id: str
    relation: Mapping[str, str]
    other: Mapping[str, Any]

    @property
    def fictitious(self) -> bool:
        return self.feats.get("Fictitious") == "Yes"

    @property
    def relation_type(self) -> str:
        return self.relation.get("Type", "")

    @property
    def techlemma(self) -> str:
        value = self.other.get("techlemma")
        return value if isinstance(value, str) else ""

    @property
    def corpus_count(self) -> int:
        stats = self.other.get("corpus_stats")
        if isinstance(stats, Mapping):
            count = stats.get("absolute_count")
            if isinstance(count, int) and not isinstance(count, bool):
                return count
        return 0

    @property
    def etymology(self) -> list[str]:
        value = self.other.get("Etymology")
        if isinstance(value, list):
            return [str(code) for code in value]
        return []


def parse_feats(text: str) -> dict[str, str]:
    """``Gender=Fem&Loanword=True`` → ``{"Gender": "Fem", ...}``."""
    out: dict[str, str] = {}
    for item in text.split("&"):
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
    return out


def parse_segmentation(text: str) -> tuple[tuple[str, str], ...]:
    """``[{"morph": "knih", "type": "R"}, …]`` → ``(("knih", "R"), …)``."""
    if not text.strip():
        return ()
    data = json.loads(text)
    if not isinstance(data, list):
        raise ValueError("SEGMENTATION is not a JSON list")
    out: list[tuple[str, str]] = []
    for item in data:
        if isinstance(item, Mapping):
            out.append((str(item.get("morph", "")), str(item.get("type", ""))))
    return tuple(out)


def parse_line(line: str, lineno: int = 0) -> Node:
    """One non-empty DeriNet line (see the module docstring)."""
    parts = line.rstrip("\n").split("\t")
    if len(parts) < N_COLUMNS:
        parts += [""] * (N_COLUMNS - len(parts))
    try:
        morphs = parse_segmentation(parts[5])
        other_raw = json.loads(parts[9]) if parts[9].strip() else {}
    except (ValueError, json.JSONDecodeError) as exc:
        raise ValueError(f"DeriNet line {lineno}: bad JSON ({exc})") from exc
    other = other_raw if isinstance(other_raw, Mapping) else {}
    return Node(
        node_id=parts[0],
        lemma=parts[2],
        upos=parts[3],
        feats=parse_feats(parts[4]),
        morphs=morphs,
        parent_id=parts[6],
        relation=parse_feats(parts[7]),
        other=other,
    )


def iter_trees(path: Path) -> Iterator[list[Node]]:
    """The trees of a DeriNet TSV, in file order."""
    tree: list[Node] = []
    with path.open("r", encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.strip():
                if tree:
                    yield tree
                    tree = []
                continue
            tree.append(parse_line(line, lineno))
    if tree:
        yield tree


def last_d_morph(morphs: Sequence[tuple[str, str]]) -> str:
    """The last morph typed D ("" if none)."""
    for morph, kind in reversed(morphs):
        if kind == "D":
            return morph
    return ""


def verbalizer(node: Node) -> str:
    """R14 verbalizer of a verb node from its segmentation ("" if the
    verb is unsegmented; ``D:<morph>`` for an unlisted morph)."""
    if node.upos != "VERB":
        return ""
    morph = last_d_morph(node.morphs)
    if not morph:
        return ""
    return VERBALIZERS.get(morph, f"D:{morph}")


def conj_consistent(verbalizer_label: str, conj_class: str) -> bool | None:
    """Does ConjugClass agree with the verbalizer? None if either is
    missing or the verbalizer has no expected class."""
    expected = EXPECTED_CLASS.get(verbalizer_label)
    if expected is None or not conj_class:
        return None
    return expected in conj_class.split("#")


def _noun_row(node: Node) -> tuple[Any, ...]:
    return (
        node.node_id,
        node.lemma,
        node.techlemma,
        GENDERS.get(node.feats.get("Gender", ""), ""),
        ANIMACIES.get(node.feats.get("Animacy", ""), ""),
        node.feats.get("Loanword", ""),
        int(bool(node.parent_id) and node.relation_type == DERIVATION),
        "|".join(node.etymology),
        node.corpus_count,
    )


def _edge_row(parent: Node, child: Node) -> tuple[Any, ...]:
    label = verbalizer(child)
    return (
        parent.node_id,
        parent.lemma,
        child.node_id,
        child.lemma,
        child.techlemma,
        UPOS.get(child.upos, child.upos.lower()),
        GENDERS.get(child.feats.get("Gender", ""), ""),
        ANIMACIES.get(child.feats.get("Animacy", ""), ""),
        child.feats.get("Loanword", ""),
        child.feats.get("ConjugClass", ""),
        child.feats.get("Aspect", ""),
        "|".join(f"{m}:{k}" for m, k in child.morphs),
        last_d_morph(child.morphs),
        label,
        child.relation.get("SemanticLabel", ""),
        child.corpus_count,
    )


@dataclass
class ProjectionCounts:
    """What one pass saw and kept."""

    nodes: int = 0
    trees: int = 0
    fictitious: int = 0
    nouns_kept: int = 0
    edges_kept: int = 0


def _tree_rows(
    tree: Sequence[Node], lemmas: frozenset[str], counts: ProjectionCounts
) -> tuple[list[tuple[Any, ...]], list[tuple[Any, ...]]]:
    by_id = {node.node_id: node for node in tree}
    nouns: list[tuple[Any, ...]] = []
    edges: list[tuple[Any, ...]] = []
    for node in tree:
        counts.nodes += 1
        if node.fictitious:
            counts.fictitious += 1
            continue
        if node.upos == "NOUN" and node.lemma in lemmas:
            nouns.append(_noun_row(node))
        if node.relation_type != DERIVATION or not node.parent_id:
            continue
        parent = by_id.get(node.parent_id)
        if (
            parent is None
            or parent.fictitious
            or parent.upos != "NOUN"
            or parent.lemma not in lemmas
        ):
            continue
        edges.append(_edge_row(parent, node))
    return nouns, edges


def project(
    ctx: RunContext,
    source: Path,
    target: Path,
    noun_lemmas: frozenset[str],
) -> ProjectionCounts:
    """Write the compact derivative *target* (module docstring)."""
    counts = ProjectionCounts()
    all_nouns: list[tuple[Any, ...]] = []
    all_edges: list[tuple[Any, ...]] = []
    for tree in iter_trees(source):
        counts.trees += 1
        nouns, edges = _tree_rows(tree, noun_lemmas, counts)
        all_nouns.extend(nouns)
        all_edges.extend(edges)
    counts.nouns_kept = len(all_nouns)
    counts.edges_kept = len(all_edges)
    meta = (("source_sha256:" + source.name, sha256_file(source)),)
    jsonio.build_sqlite(
        ctx,
        target,
        table="nouns",
        columns=NOUN_COLUMNS,
        rows=all_nouns,
        order_by=[name for name, _ in NOUN_COLUMNS],
        indexes=(("nouns_lemma", ("lemma",)),),
        meta=meta,
    )
    _add_table(target, "edges", EDGE_COLUMNS, all_edges)
    return counts


def _add_table(
    path: Path,
    table: str,
    columns: Sequence[tuple[str, str]],
    rows: Sequence[tuple[Any, ...]],
) -> None:
    """Add a second table to a database built by ``build_sqlite``.

    Rows are sorted in Python first, so the file stays deterministic.
    """
    column_sql = ", ".join(f"{name} {kind}" for name, kind in columns)
    marks = ", ".join("?" for _ in columns)
    ordered = sorted(rows, key=lambda row: tuple(str(v) for v in row))
    with contextlib.closing(sqlite3.connect(path)) as con:
        con.execute("PRAGMA journal_mode = OFF")
        con.execute(f"CREATE TABLE {table} ({column_sql})")
        con.executemany(f"INSERT INTO {table} VALUES ({marks})", ordered)
        con.execute(
            f"CREATE INDEX {table}_parent ON {table} ({columns[0][0]})"
        )
        con.execute(
            "INSERT INTO meta VALUES (?, ?)", (f"{table}_rows", str(len(rows)))
        )
        con.commit()


@dataclass(frozen=True, slots=True)
class DerinetNoun:
    """A projected DeriNet noun node."""

    node_id: str
    lemma: str
    gender: str
    animacy: str
    loanword: str
    has_parent: bool
    etymology: str
    corpus_count: int

    @property
    def loan_root(self) -> bool | None:
        """Loanword=True on a node without a derivational parent."""
        if self.has_parent:
            return None
        if self.loanword == "True":
            return True
        return None


@dataclass(frozen=True, slots=True)
class DerinetEdge:
    """A projected noun-parented derivation edge."""

    parent_id: str
    parent_lemma: str
    child_id: str
    child_lemma: str
    child_pos: str
    child_loanword: str
    conj_class: str
    aspect: str
    segmentation: str
    last_d_morph: str
    verbalizer: str
    semantic_label: str
    corpus_count: int


def read_projection(
    path: Path,
) -> tuple[list[DerinetNoun], list[DerinetEdge]]:
    """The two tables of a projection, in stored order."""
    uri = path.resolve().as_uri() + "?mode=ro"
    with contextlib.closing(sqlite3.connect(uri, uri=True)) as con:
        nouns = [
            DerinetNoun(
                node_id=row[0],
                lemma=row[1],
                gender=row[3],
                animacy=row[4],
                loanword=row[5],
                has_parent=bool(row[6]),
                etymology=row[7],
                corpus_count=int(row[8]),
            )
            for row in con.execute("SELECT * FROM nouns ORDER BY rowid")
        ]
        edges = [
            DerinetEdge(
                parent_id=row[0],
                parent_lemma=row[1],
                child_id=row[2],
                child_lemma=row[3],
                child_pos=row[5],
                child_loanword=row[8],
                conj_class=row[9],
                aspect=row[10],
                segmentation=row[11],
                last_d_morph=row[12],
                verbalizer=row[13],
                semantic_label=row[14],
                corpus_count=int(row[15]),
            )
            for row in con.execute("SELECT * FROM edges ORDER BY rowid")
        ]
    return nouns, edges


def compatible(readings: Sequence[str], gender: str, animacy: str) -> bool:
    """Does a DeriNet noun (gender, animacy) fit a lexeme's readings
    (``m.an``, ``f`` …)? Missing values on either side match (R14)."""
    if not gender or not readings:
        return True
    for reading in readings:
        parts = reading.split(".")
        if parts[0] != gender:
            continue
        if gender != "m" or len(parts) < 2 or not animacy:
            return True
        if parts[1] == animacy:
            return True
    return False
