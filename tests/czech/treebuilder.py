"""Build a fixture-sized Czech source tree and run the plugin's stages.

Used by the Czech tests and, in a subprocess, by the determinism test.
The Wiktionary inputs are the tracked CC BY-SA fixtures under
``tests/fixtures/czech/``; MorfFlex, Leipzig and CNK inputs are
**synthetic** tables written here (no NC or unverified data is
tracked). Layout mirrors the real ``data/sources/`` so the stages run
unchanged with ``params={"check_source_sizes": False}``.
"""

from __future__ import annotations

import gzip
import io
import lzma
import random
import shutil
import sys
import tarfile
import zipfile
from collections.abc import Sequence
from pathlib import Path
from typing import Final

PROJECT_ROOT: Final[Path] = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lexpipe import spec  # noqa: E402
from lexpipe.context import RunContext  # noqa: E402
from lexpipe.langs import czech  # noqa: E402
from lexpipe.langs.czech import adjudication, layout  # noqa: E402

FIXTURES: Final[Path] = PROJECT_ROOT / "tests" / "fixtures" / "czech"
EN_FIXTURE: Final[Path] = FIXTURES / "en_wiktionary_cs_sample.jsonl"
CS_FIXTURE: Final[Path] = FIXTURES / "cs_extract_sample.jsonl"
STAGES: Final[tuple[str, ...]] = (
    "extract",
    "lexicon",
    "qc",
    "attest",
    "features",
)
PARAMS: Final[dict[str, object]] = {"check_source_sizes": False}

# Synthetic MorfFlex rows (lemma, tag, form): made up for the tests,
# shaped like MorfFlex CZ 2.1 rows (TR-64 tags).
MORFFLEX_ROWS: Final[tuple[tuple[str, str, str], ...]] = (
    ("kluk", "NNMS1-----A----", "kluk"),
    ("kluk", "NNMS2-----A----", "kluka"),
    ("kluk", "NNMP1-----A----", "kluci"),
    ("kluk", "NNMP1-----A---6", "klucí"),
    ("kluk", "NNMP6-----A----", "klucích"),
    ("kluk", "NNMP1-----N----", "nekluci"),
    ("biolog", "NNMP1-----A----", "biologové"),
    ("biolog", "NNMS2-----A----", "biologa"),
    ("biolog", "NNMP6-----A----", "biologách"),
    ("Čech_;E_;Y", "NNMP1-----A----", "Češi"),
    ("Čech_;E_;Y", "NNMP1-----A---1", "Čechové"),
    ("Čech_;E_;Y", "NNMP6-----A----", "Češích"),
    ("šach", "NNIS6-----A----", "šachu"),
    ("šach", "NNIP6-----A----", "šachách"),
    ("šach-2", "NNMP1-----A----", "šachové"),
    ("šach-2", "NNMP6-----A----", "šaších"),
    ("jazyk", "NNIS6-----A----", "jazyku"),
    ("jazyk", "NNIS6-----A---1", "jazyce"),
    ("jazyk-2", "NNMP1-----A----", "jazykové"),
    ("potok", "NNIS6-----A----", "potoku"),
    ("potok", "NNIS6-----A---1", "potoce"),
    ("park", "NNIS6-----A----", "parku"),
    ("park", "NNIP6-----A----", "parcích"),
    ("rok", "NNIS6-----A----", "roku"),
    ("rok", "NNIS6-----A---6", "roce"),
    ("student", "NNMP1-----A----", "studenti"),
    ("liga_,h", "NNFS3-----A----", "lize"),
    ("český", "AAMP1----1A----", "čeští"),
    ("český", "AAMS1----1A----", "český"),
    ("kupovat_:T", "VB-S---3P-AA---", "kupuje"),
    ("Abacus_;m_,s_^(^DD**Abakus)", "NNMS1-----A----", "Abacus"),
)

# Synthetic Leipzig counts (word, frequency); the filler words make the
# corpus large enough for the per-million thresholds to differ.
LEIPZIG_COUNTS: Final[tuple[tuple[str, int], ...]] = (
    ("a", 1_000_000),
    ("je", 500_000),
    (".", 90_000),
    ("2024", 400),
    ("kluk", 346),
    ("kluci", 870),
    ("Kluci", 30),
    ("klucích", 7),
    ("biolog", 32),
    ("biologové", 7),
    ("Čech", 1300),
    ("Češi", 2000),
    ("Čechové", 9),
    ("Češích", 5),
    ("šachu", 16),
    ("šachové", 7),
    ("Šachové", 3),
    ("šachách", 2),
    ("potok", 53),
    ("potoce", 25),
    ("potoku", 11),
    ("potocích", 3),
    ("student", 298),
    ("studenti", 543),
    ("diplomaté", 200),
    ("diplomati", 20),
    ("rok", 8000),
    ("roce", 18000),
    ("roku", 17000),
    ("rocích", 50),
    ("jazyce", 1600),
    ("jazyku", 1200),
    ("parku", 1400),
    ("parcích", 400),
    ("knize", 300),
    ("noze", 200),
    ("lize", 1984),
    ("vrazi", 259),
    ("vrahové", 60),
    ("žáci", 500),
    ("čeští", 3000),
    ("plechu", 10),
    ("otče", 4),
    ("bratře", 3),
)

CNK_ROWS: Final[tuple[tuple[str, int], ...]] = (
    ("rok", 358211),
    ("kluk", 16297),
    ("student", 15649),
    ("park", 11521),
    ("Čech", 7660),
    ("potok", 3385),
    ("plech", 2344),
    ("šach", 1046),
    ("biolog", 949),
    ("liga", 11783),
)


def _shuffled(lines: list[bytes], seed: int | None) -> list[bytes]:
    if seed is not None:
        random.Random(seed).shuffle(lines)
    return lines


def _write_bytes(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)


def build_tree(root: Path, *, shuffle_seed: int | None = None) -> Path:
    """Create the source tree under *root*; shuffle inputs if asked."""
    sources = root / "data" / "sources"
    en_lines = _shuffled(
        EN_FIXTURE.read_bytes().splitlines(keepends=True), shuffle_seed
    )
    _write_bytes(
        sources / "kaikki" / "en_wiktionary_cs.jsonl", b"".join(en_lines)
    )
    cs_lines = _shuffled(
        CS_FIXTURE.read_bytes().splitlines(keepends=True), shuffle_seed
    )
    _write_bytes(
        sources / "kaikki" / "cs-extract.jsonl.gz",
        gzip.compress(b"".join(cs_lines), mtime=0),
    )
    morf = _shuffled(
        [f"{a}\t{b}\t{c}\n".encode() for a, b, c in MORFFLEX_ROWS],
        shuffle_seed,
    )
    _write_bytes(
        sources / "morfflex" / "czech-morfflex-2.1.tsv.xz",
        lzma.compress(b"".join(morf)),
    )
    _write_leipzig(
        sources / "leipzig" / "ces_news_2024_1M.tar.gz", shuffle_seed
    )
    _write_cnk(sources / "cnk" / "syn2015_lemma_utf8.zip", shuffle_seed)
    curated = root / "curated" / "czech" / layout.ADJUDICATIONS_FILE
    curated.parent.mkdir(parents=True, exist_ok=True)
    curated.write_text(",".join(adjudication.FIELDS) + "\n", encoding="utf-8")
    return root


def _write_leipzig(path: Path, seed: int | None) -> None:
    corpus = "ces_news_2024_1M"
    rows = _shuffled(
        [
            f"{i}\t{w}\t{f}\n".encode()
            for i, (w, f) in enumerate(LEIPZIG_COUNTS, 1)
        ],
        seed,
    )
    members = {
        f"{corpus}/{corpus}-words.txt": b"".join(rows),
        f"{corpus}/{corpus}-meta.txt": b"1\tsynthetic test corpus\n",
    }
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as tar:
        for name, data in sorted(members.items()):
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    _write_bytes(path, buffer.getvalue())


def _write_cnk(path: Path, seed: int | None) -> None:
    rows = _shuffled(
        [
            f"{i}\t{lemma}\t{freq}\t{freq}\t1\t1\t1\t\n".encode()
            for i, (lemma, freq) in enumerate(CNK_ROWS, 1)
        ],
        seed,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("syn2015_lemma_utf8.tsv", b"".join(rows))


def context(root: Path, force: Sequence[str] = ()) -> RunContext:
    """The RunContext of a fixture run."""
    return RunContext(
        spec.LANGUAGES["cs"], root=root, params=PARAMS, force=frozenset(force)
    )


def run_stages(root: Path, stages: Sequence[str] = STAGES) -> None:
    """Run *stages* in order on the tree; raise on a failed stage."""
    ctx = context(root)
    for name in stages:
        status = czech.STAGES[name](ctx)
        if isinstance(status, int) and status != 0:
            raise RuntimeError(f"stage {name} failed with {status}")


def stage_outputs(root: Path) -> dict[str, list[Path]]:
    """Outputs of every stage, as declared by their modules."""
    paths = layout.paths_for(context(root))
    freq = paths.leipzig.freq_csv("ces_news_2024_1M")
    return {
        "extract": [
            paths.en_entries,
            paths.cs_entries,
            paths.form_index,
            paths.extract_stats,
        ],
        "lexicon": [
            paths.lexicon_full,
            paths.lexicon_raw,
            paths.cells_raw,
            paths.lexicon_stats,
        ],
        "qc": [
            paths.morfflex_db,
            paths.morfflex_lemmas,
            paths.qc_csv,
            paths.qc_report,
            paths.worklist,
        ],
        "attest": [
            freq,
            paths.attestation,
            paths.lexeme_frequency,
            paths.attestation_report,
        ],
        "features": [
            paths.cells,
            paths.lexicon_complete,
            paths.features_report,
        ],
    }


def copy_fixture_tree(source: Path, target: Path) -> Path:
    """Copy a built tree (tests that modify inputs work on a copy)."""
    shutil.copytree(source, target)
    return target
