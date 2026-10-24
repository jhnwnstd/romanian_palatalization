"""Fixture trees for the Italian analysis stages (derivations … docs).

:func:`build_analysis_tree` lays out the data-layer fixture tree
(:func:`tests.italian.helpers.build_tree`) and adds

  - the analysis excerpt ``tests/fixtures/italian/
    en_wiktionary_it_analysis_sample.jsonl`` to the en split (derivatives
    such as grecista, fochista, antichità, leghista/legista, musicista,
    politichese, grecizzare, porcino; the verbs vincere, cercare, tacere,
    cuocere, cucire, fuggire, dire, bianchire, trafficare, legare; the
    bases unico, tipico, lega, musica, politica, porco, bianco, cieco);
  - a words table with **made-up** counts for those forms
    (:data:`ANALYSIS_WORDS`, appended to the data layer's
    :data:`tests.italian.helpers.SYNTHETIC_WORDS`; no Leipzig data is
    redistributed, R9).

With a shuffle seed every input's lines (and the words table's rows) are
shuffled, for the determinism tests (R5). :func:`make_context` passes
small learner settings (:data:`LEARNER_V`, :data:`LEARNER_SAMPLES`) and
few permutations, so the stages run in seconds.
"""

from __future__ import annotations

import io
import random
import tarfile
from collections.abc import Sequence
from pathlib import Path
from typing import Final

from lexpipe.context import RunContext
from lexpipe.spec import LANGUAGES
from tests.italian import helpers

ANALYSIS_SAMPLE: Final[Path] = (
    helpers.FIXTURES / "en_wiktionary_it_analysis_sample.jsonl"
)
DRIVER: Final[Path] = Path(__file__).resolve().parent / "analysis_driver.py"
LEARNER_V: Final[tuple[int, ...]] = (5, 10)
LEARNER_SAMPLES: Final[int] = 40
PERMUTATIONS: Final[int] = 200
ANALYSIS_STAGES: Final[tuple[str, ...]] = (
    "derivations",
    "verbs",
    "tp",
    "covariation",
    "spellout",
    "noevidence",
    "logtable",
    "docs",
)

# Invented token counts (Leipzig words-table layout).
ANALYSIS_WORDS: Final[tuple[tuple[str, int], ...]] = (
    ("unico", 300),
    ("unici", 40),
    ("unica", 100),
    ("uniche", 30),
    ("unicità", 30),
    ("tipico", 30),
    ("tipici", 10),
    ("tipica", 20),
    ("tipiche", 8),
    ("tipicità", 5),
    ("tipicizzato", 2),
    ("lega", 200),
    ("leghe", 5),
    ("leghista", 50),
    ("leghisti", 10),
    ("musica", 400),
    ("musiche", 20),
    ("musicista", 40),
    ("musicisti", 12),
    ("musico", 3),
    ("musici", 2),
    ("politica", 600),
    ("politiche", 100),
    ("politico", 150),
    ("politici", 150),
    ("politichese", 2),
    ("grecista", 2),
    ("grecismo", 3),
    ("antichità", 20),
    ("antichizzare", 1),
    ("fochista", 1),
    ("amichevole", 60),
    ("amichevoli", 10),
    ("amichetto", 5),
    ("medichessa", 1),
    ("porco", 3),
    ("porci", 5),
    ("porcino", 8),
    ("porcini", 6),
    ("porcello", 3),
    ("porcile", 4),
    ("bianco", 300),
    ("bianchi", 50),
    ("bianca", 150),
    ("bianche", 40),
    ("bianchezza", 2),
    ("cieco", 10),
    ("ciechi", 5),
    ("cecità", 6),
    ("opaco", 5),
    ("opachi", 2),
    ("opacità", 4),
    ("traffico", 80),
    ("traffici", 15),
    ("vinco", 8),
    ("vinci", 30),
    ("vince", 100),
    ("vinciamo", 5),
    ("vinca", 10),
    ("cerco", 20),
    ("cerchi", 40),
    ("cerca", 200),
    ("cerchiamo", 20),
    ("taci", 5),
    ("tace", 10),
    ("taccio", 2),
    ("cuoce", 3),
    ("cuci", 1),
    ("fuggi", 4),
    ("fugge", 6),
    ("dici", 20),
    ("dice", 400),
    ("diciamo", 50),
    ("dica", 30),
    ("leghi", 3),
    ("astrologo", 3),
    ("astrologi", 2),
    ("astrologhi", 1),
    ("dialogo", 50),
    ("dialoghi", 10),
    ("teologia", 5),
)


def _words_table(seed: int | None) -> bytes:
    merged: dict[str, int] = {}
    for word, count in (*helpers.SYNTHETIC_WORDS, *ANALYSIS_WORDS):
        merged[word] = merged.get(word, 0) + count
    rows = [
        f"{index}\t{word}\t{count}\n".encode()
        for index, (word, count) in enumerate(sorted(merged.items()), 1)
    ]
    if seed is not None:
        random.Random(seed).shuffle(rows)
    return b"".join(rows)


def _archive(seed: int | None) -> bytes:
    data = _words_table(seed)
    buffer = io.BytesIO()
    corpus = helpers.CORPUS_ID
    with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
        info = tarfile.TarInfo(f"{corpus}/{corpus}-words.txt")
        info.size = len(data)
        info.mtime = 0
        archive.addfile(info, io.BytesIO(data))
    return buffer.getvalue()


def _lines(path: Path) -> list[bytes]:
    return [
        line if line.endswith(b"\n") else line + b"\n"
        for line in path.read_bytes().splitlines(keepends=True)
        if line.strip()
    ]


def build_analysis_tree(root: Path, shuffle_seed: int | None = None) -> Path:
    """The data-layer tree plus the analysis excerpt and words."""
    helpers.build_tree(root, shuffle_seed)
    lines = _lines(helpers.EN_SAMPLE) + _lines(ANALYSIS_SAMPLE)
    if shuffle_seed is not None:
        random.Random(shuffle_seed).shuffle(lines)
    split = root / "data" / "sources" / "kaikki" / "en_wiktionary_it.jsonl"
    split.write_bytes(b"".join(lines))
    archive = (
        root / "data" / "sources" / "leipzig" / f"{helpers.CORPUS_ID}.tar.gz"
    )
    archive.write_bytes(_archive(shuffle_seed))
    return root


def make_context(root: Path, force: Sequence[str] = ()) -> RunContext:
    """A run context with small learner and permutation settings."""
    return RunContext(
        language=LANGUAGES["it"],
        root=root,
        force=frozenset(force),
        params={
            "verify_source_size": False,
            "learner_samples": LEARNER_SAMPLES,
            "learner_V": list(LEARNER_V),
            "permutations": PERMUTATIONS,
        },
    )
