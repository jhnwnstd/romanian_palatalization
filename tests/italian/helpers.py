"""Fixture trees for the Italian data-layer tests.

:func:`build_tree` lays out a throwaway repository root holding the
fixture excerpts where the pipeline expects its pinned sources
(``data/sources/…``, see ``tests/fixtures/italian/README.md``), plus a
Leipzig-format corpus archive whose words table has **made-up** counts
(:data:`SYNTHETIC_WORDS`; no Leipzig data is redistributed, R9) and a
copy of ``docs/italian/schema.toml`` (read by the ``docs`` stage). With
a shuffle seed, the lines of every input are shuffled, for the
determinism tests (R5).
"""

from __future__ import annotations

import gzip
import io
import random
import sys
import tarfile
from collections.abc import Sequence
from pathlib import Path
from typing import Final

PROJECT_ROOT: Final[Path] = Path(__file__).resolve().parents[2]
SRC: Final[Path] = PROJECT_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from lexpipe.context import RunContext  # noqa: E402
from lexpipe.spec import LANGUAGES  # noqa: E402

FIXTURES: Final[Path] = PROJECT_ROOT / "tests" / "fixtures" / "italian"
EN_SAMPLE: Final[Path] = FIXTURES / "en_wiktionary_it_sample.jsonl"
IT_SAMPLE: Final[Path] = FIXTURES / "it_extract_sample.jsonl"
MORPHIT_SAMPLE: Final[Path] = FIXTURES / "morphit_sample.txt"
CORPUS_ID: Final[str] = LANGUAGES["it"].leipzig_id
DRIVER: Final[Path] = Path(__file__).resolve().parent / "pipeline_driver.py"
# The schema the ``docs`` stage checks every header against (tracked).
SCHEMA: Final[Path] = PROJECT_ROOT / "docs" / "italian" / "schema.toml"

# Invented token counts in Leipzig's words-table layout (id, word,
# count). Two fillers make the corpus about 1.5 million word tokens, so
# that 0.5 per million is 0.75 tokens and 2 per million 3 tokens.
SYNTHETIC_WORDS: Final[tuple[tuple[str, int], ...]] = (
    ("il", 1_000_000),
    ("di", 500_000),
    (",", 50_000),
    ("amico", 120),
    ("Amico", 3),
    ("l'amico", 25),
    ("amici", 300),
    ("amica", 60),
    ("amiche", 45),
    ("medico", 200),
    ("medici", 150),
    ("Medici", 20),
    ("mediche", 10),
    ("parco", 90),
    ("parchi", 30),
    ("lago", 70),
    ("laghi", 25),
    ("chirurgo", 30),
    ("chirurghi", 40),
    ("carico", 80),
    ("carichi", 15),
    ("stomaco", 25),
    ("stomaci", 3),
    ("stomachi", 2),
    ("mago", 10),
    ("maghi", 8),
    ("magi", 4),
    ("Magi", 6),
    ("pubblico", 500),
    ("pubblici", 200),
    ("pubbliche", 100),
    ("elettrico", 50),
    ("elettrici", 30),
    ("elettriche", 25),
    ("greco", 20),
    ("greci", 30),
    ("greche", 8),
    ("belga", 15),
    ("belgi", 20),
    ("belghe", 2),
    ("luce", 300),
    ("luci", 100),
    ("occhio", 100),
    ("occhi", 400),
    ("teologo", 5),
    ("teologi", 6),
    ("antico", 60),
    ("antichi", 30),
    ("antiche", 30),
    ("valigia", 20),
    ("valigie", 10),
    ("fuoco", 150),
    ("fuochi", 20),
    ("catalogo", 30),
    ("cataloghi", 10),
    ("psiche", 5),
    ("specie", 80),
)


def _lines(path: Path, seed: int | None) -> list[bytes]:
    lines = [
        line if line.endswith(b"\n") else line + b"\n"
        for line in path.read_bytes().splitlines(keepends=True)
        if line.strip()
    ]
    if seed is not None:
        random.Random(seed).shuffle(lines)
    return lines


def words_table(seed: int | None = None) -> bytes:
    """The synthetic words table (optionally in shuffled row order)."""
    rows = [
        f"{index}\t{word}\t{count}\n".encode()
        for index, (word, count) in enumerate(SYNTHETIC_WORDS, 1)
    ]
    if seed is not None:
        random.Random(seed).shuffle(rows)
    return b"".join(rows)


def _corpus_archive(seed: int | None) -> bytes:
    data = words_table(seed)
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
        info = tarfile.TarInfo(f"{CORPUS_ID}/{CORPUS_ID}-words.txt")
        info.size = len(data)
        info.mtime = 0
        archive.addfile(info, io.BytesIO(data))
    return buffer.getvalue()


def build_tree(
    root: Path,
    shuffle_seed: int | None = None,
    *,
    curated: dict[str, str] | None = None,
) -> Path:
    """Lay out the fixture sources under *root*; return *root*.

    *curated* maps file names under ``curated/italian/`` to their text.
    """
    kaikki = root / "data" / "sources" / "kaikki"
    kaikki.mkdir(parents=True)
    split = kaikki / "en_wiktionary_it.jsonl"
    split.write_bytes(b"".join(_lines(EN_SAMPLE, shuffle_seed)))
    extract = kaikki / "it-extract.jsonl.gz"
    with gzip.GzipFile(extract, "wb", mtime=0) as fh:
        fh.write(b"".join(_lines(IT_SAMPLE, shuffle_seed)))
    morphit = root / "data" / "sources" / "morphit"
    morphit.mkdir(parents=True)
    (morphit / "morph-it_048_utf8.txt").write_bytes(
        b"".join(_lines(MORPHIT_SAMPLE, shuffle_seed))
    )
    leipzig = root / "data" / "sources" / "leipzig"
    leipzig.mkdir(parents=True)
    (leipzig / f"{CORPUS_ID}.tar.gz").write_bytes(
        _corpus_archive(shuffle_seed)
    )
    for name, text in (curated or {}).items():
        target = root / "curated" / "italian" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text, encoding="utf-8")
    docs = root / "docs" / "italian"
    docs.mkdir(parents=True, exist_ok=True)
    (docs / "schema.toml").write_bytes(SCHEMA.read_bytes())
    return root


def make_context(root: Path, force: Sequence[str] = ()) -> RunContext:
    """A run context for a fixture tree (no source-size checks)."""
    return RunContext(
        language=LANGUAGES["it"],
        root=root,
        force=frozenset(force),
        params={"verify_source_size": False},
    )


def run_pipeline(ctx: RunContext, stages: Sequence[str] = ()) -> None:
    """Run the plugin's stages (all, or the named ones) in order."""
    from lexpipe.langs.italian import PLUGIN

    for name, stage in PLUGIN.stages.items():
        if stages and name not in stages:
            continue
        status = stage(ctx)
        if status not in (None, 0):
            raise AssertionError(f"stage {name} returned {status}")
