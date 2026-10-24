"""Tests for lexpipe.frequency on a synthetic Leipzig archive.

The fixture archive has the real layout
(``<id>/<id>-words.txt``, ``-sentences.txt``, ``-meta.txt``) and a
words table exercising every cleaning rule. The cleaning is checked
against ``scripts/download_leipzig.py`` itself.
"""

from __future__ import annotations

import importlib.util
import io
import sys
import tarfile
import unicodedata
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lexpipe import frequency  # noqa: E402
from lexpipe.frequency import LeipzigLayout  # noqa: E402

CORPUS = "xyz_news_2024_10K"
DECOMPOSED_CITTA = unicodedata.normalize("NFD", "città")

WORDS_TABLE = "\n".join(
    [
        "id\tword\tfreq",  # header row
        "1\t!\t100",  # no letter
        "2\tcasa\t50",
        "3\tCasa\t7",  # lowercases onto casa
        "4\tcittà\t9",
        "5\tcitta\t4",  # accents are kept: a different key
        f"6\t{DECOMPOSED_CITTA}\t1",  # NFC onto città
        '7\t"perché"\t3',  # quotes trimmed
        "8\thttp://x.it\t5",  # URL
        "9\ta b\t2",  # whitespace inside
        "10\tdell'\t6",  # apostrophe trimmed
        "11\t123\t8",  # digits only
        "12\tzero\t0",  # non-positive count
        "13\tbad\tx",  # non-numeric count
        "14\tshort",  # too few columns
        "15\tŽLUŤOUČKÝ\t2",  # Czech diacritics kept
        "16\tfoo@bar\t3",  # e-mail-like
        "17\tsquare\t²",  # a digit int() cannot parse
    ]
)
EXPECTED_TABLE = {
    "casa": 57,
    "città": 10,
    "dell": 6,
    "citta": 4,
    "perché": 3,
    "žluťoučký": 2,
}


def _add(tf: tarfile.TarFile, name: str, data: bytes) -> None:
    info = tarfile.TarInfo(name)
    info.size = len(data)
    tf.addfile(info, io.BytesIO(data))


def make_archive(path: Path, members: dict[str, str]) -> Path:
    with tarfile.open(path, "w:gz") as tf:
        directory = tarfile.TarInfo(CORPUS)
        directory.type = tarfile.DIRTYPE
        tf.addfile(directory)
        for name, text in members.items():
            _add(tf, name, text.encode("utf-8"))
    return path


@pytest.fixture
def archive(tmp_path: Path) -> Path:
    return make_archive(
        tmp_path / f"{CORPUS}.tar.gz",
        {
            f"{CORPUS}/{CORPUS}-meta.txt": "1\tbuild date\t2025-02-12\n",
            f"{CORPUS}/{CORPUS}-sentences.txt": "1\tUna casa.\n",
            f"{CORPUS}/{CORPUS}-words.txt": WORDS_TABLE + "\n",
        },
    )


@pytest.fixture
def layout(tmp_path: Path) -> LeipzigLayout:
    return frequency.language_layout("xyzlang", tmp_path)


def _quiet(_: str) -> None:
    pass


@pytest.fixture(scope="module")
def romanian() -> Iterator[ModuleType]:
    """scripts/download_leipzig.py, imported by path (nothing runs)."""
    name = "_lexpipe_test_download_leipzig"
    path = PROJECT_ROOT / "scripts" / "download_leipzig.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
        yield module
    finally:
        sys.modules.pop(name, None)


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------


def test_language_layout_paths(tmp_path: Path) -> None:
    layout = frequency.language_layout("italian", tmp_path)
    root = tmp_path / "data" / "italian" / "leipzig"
    assert layout.root == root
    assert layout.freq_csv("ita_news_2024_1M") == (
        root / "freq" / "ita_news_2024_1M_freq.csv"
    )
    assert layout.corpus_dir("ita_news_2024_1M") == (
        root / "corpora" / "ita_news_2024_1M"
    )
    assert layout.archive_path("ita_news_2024_1M") == (
        root / "raw" / "ita_news_2024_1M.tar.gz"
    )


def test_layout_refuses_the_romanian_directory(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="attach_frequencies"):
        LeipzigLayout(tmp_path / "data" / "leipzig")
    with pytest.raises(ValueError, match="attach_frequencies"):
        LeipzigLayout(PROJECT_ROOT / "data" / "leipzig")
    with pytest.raises(ValueError, match="attach_frequencies"):
        LeipzigLayout(PROJECT_ROOT / "data" / "italian" / ".." / "leipzig")


@pytest.mark.parametrize("bad", ["", "../x", "Italian", "it/x"])
def test_language_layout_rejects_odd_names(bad: str) -> None:
    with pytest.raises(ValueError):
        frequency.language_layout(bad)


@pytest.mark.parametrize("bad", ["", "../ita", "ita news", "/abs"])
def test_corpus_ids_are_validated(bad: str) -> None:
    with pytest.raises(ValueError, match="corpus id"):
        frequency.corpus_url(bad)


def test_corpus_url() -> None:
    assert frequency.corpus_url("ces_news_2024_1M") == (
        "https://downloads.wortschatz-leipzig.de/corpora/"
        "ces_news_2024_1M.tar.gz"
    )


# ---------------------------------------------------------------------------
# Building a table
# ---------------------------------------------------------------------------


def test_build_freq_csv(archive: Path, layout: LeipzigLayout) -> None:
    out = frequency.build_freq_csv(
        CORPUS,
        layout,
        archive=archive,
        expected_size=archive.stat().st_size,
        log=_quiet,
    )
    assert out == layout.freq_csv(CORPUS)
    assert out.read_text(encoding="utf-8").splitlines() == [
        "word,freq",
        "casa,57",
        "città,10",
        "dell,6",
        "citta,4",
        "perché,3",
        "žluťoučký,2",
    ]
    corpus_dir = layout.corpus_dir(CORPUS)
    extracted = sorted(
        p.relative_to(corpus_dir).as_posix()
        for p in corpus_dir.rglob("*")
        if p.is_file()
    )
    # Only the words table and the meta file are extracted.
    assert extracted == [
        f"{CORPUS}/{CORPUS}-meta.txt",
        f"{CORPUS}/{CORPUS}-words.txt",
    ]
    assert not list(layout.corpora_dir.glob("*.extracting"))


def test_parse_words_file(archive: Path, layout: LeipzigLayout) -> None:
    words = frequency.extract_words_file(
        archive, layout.corpus_dir(CORPUS), CORPUS
    )
    assert frequency.parse_words_file(words) == EXPECTED_TABLE


def test_cleaning_matches_romanian_script(
    archive: Path,
    layout: LeipzigLayout,
    tmp_path: Path,
    romanian: ModuleType,
) -> None:
    words = frequency.extract_words_file(
        archive, layout.corpus_dir(CORPUS), CORPUS
    )
    # The Romanian parser cannot read the "²" count (int() raises), so
    # compare on the table without that row.
    portable = tmp_path / "portable-words.txt"
    portable.write_text(
        "".join(
            line
            for line in words.read_text(encoding="utf-8").splitlines(True)
            if "²" not in line
        ),
        encoding="utf-8",
    )
    assert frequency.parse_words_file(portable) == (
        romanian.parse_words_file(portable)
    )
    for token in [
        "casa",
        "Casa",
        " «Città» ",
        DECOMPOSED_CITTA,
        "dell'",
        "www.example.org",
        "a@b",
        "1984",
        "e-mail",
        "...",
        "Škoda!",
    ]:
        assert frequency.clean_token(token) == romanian.clean_token(token)


def test_extraction_is_reused(archive: Path, layout: LeipzigLayout) -> None:
    first = frequency.extract_words_file(
        archive, layout.corpus_dir(CORPUS), CORPUS
    )
    archive.unlink()  # a second call must not need the archive
    second = frequency.extract_words_file(
        archive, layout.corpus_dir(CORPUS), CORPUS
    )
    assert first == second


def test_archive_without_inner_directory(
    tmp_path: Path, layout: LeipzigLayout
) -> None:
    flat = make_archive(
        tmp_path / "flat.tar.gz", {f"{CORPUS}-words.txt": "1\tcasa\t3\n"}
    )
    words = frequency.extract_words_file(
        flat, layout.corpus_dir(CORPUS), CORPUS
    )
    assert words == layout.corpus_dir(CORPUS) / f"{CORPUS}-words.txt"


def test_extraction_refuses_path_traversal(
    tmp_path: Path, layout: LeipzigLayout
) -> None:
    evil = make_archive(
        tmp_path / "evil.tar.gz", {f"../{CORPUS}-words.txt": "1\tcasa\t3\n"}
    )
    corpus_dir = layout.corpus_dir(CORPUS)
    with pytest.raises(tarfile.FilterError):
        frequency.extract_words_file(evil, corpus_dir, CORPUS)
    assert not (layout.corpora_dir / f"{CORPUS}-words.txt").exists()
    assert not corpus_dir.exists()
    assert not list(layout.corpora_dir.glob("*.extracting"))


def test_archive_without_words_table(
    tmp_path: Path, layout: LeipzigLayout
) -> None:
    bare = make_archive(
        tmp_path / "bare.tar.gz", {f"{CORPUS}/{CORPUS}-sentences.txt": "x"}
    )
    with pytest.raises(FileNotFoundError, match="words"):
        frequency.extract_words_file(bare, layout.corpus_dir(CORPUS), CORPUS)
    assert not layout.corpus_dir(CORPUS).exists()


def test_empty_table_is_an_error(
    tmp_path: Path, layout: LeipzigLayout
) -> None:
    junk = make_archive(
        tmp_path / "junk.tar.gz",
        {f"{CORPUS}/{CORPUS}-words.txt": "1\t!!!\t5\n2\t42\t3\n"},
    )
    with pytest.raises(ValueError, match="no usable rows"):
        frequency.build_freq_csv(CORPUS, layout, archive=junk, log=_quiet)
    assert not layout.freq_csv(CORPUS).exists()


def test_supplied_archive_size_is_checked(
    archive: Path, layout: LeipzigLayout
) -> None:
    size = archive.stat().st_size
    with pytest.raises(ValueError, match="incomplete download"):
        frequency.ensure_archive(
            CORPUS, layout, archive=archive, expected_size=size + 1
        )
    with pytest.raises(FileNotFoundError):
        frequency.ensure_archive(
            CORPUS, layout, archive=archive.with_name("missing.tar.gz")
        )
    assert (
        frequency.ensure_archive(
            CORPUS, layout, archive=archive, expected_size=size
        )
        == archive
    )


# ---------------------------------------------------------------------------
# Lookups
# ---------------------------------------------------------------------------


def test_load_freq_table_and_lookups(tmp_path: Path) -> None:
    table_path = tmp_path / "x_freq.csv"
    table_path.write_text(
        "word,freq\n"
        "biologa,5\n"
        "Biologa,1\n"
        "biolog,10\n"
        "biolozích,2\n"
        f"{DECOMPOSED_CITTA},4\n"
        "città,1\n"
        ",9\n"
        "zero,0\n"
        "broken,x\n",
        encoding="utf-8",
    )
    table = frequency.load_freq_table(table_path)
    assert table == {"biologa": 6, "biolog": 10, "biolozích": 2, "città": 5}
    assert frequency.lemma_frequency(" Biolog ", table) == 10
    assert frequency.lemma_frequency("unattested", table) == 0
    assert frequency.lemma_frequency("CITTÀ", table) == 5
    assert frequency.lemma_frequency("citta", table) == 0
    # biologa fills GEN.SG and ACC.SG but is counted once.
    forms = ["biolog", "biologa", "biologa", "BIOLOGA", "biolozích", ""]
    assert frequency.paradigm_frequency(forms, table) == 18


def test_freq_csv_round_trips_through_load(
    archive: Path, layout: LeipzigLayout
) -> None:
    out = frequency.build_freq_csv(CORPUS, layout, archive=archive, log=_quiet)
    assert frequency.load_freq_table(out) == EXPECTED_TABLE


# ---------------------------------------------------------------------------
# Case-sensitive counts and totals
# ---------------------------------------------------------------------------


def test_count_words_keeps_case_and_totals(
    archive: Path, layout: LeipzigLayout
) -> None:
    words = frequency.extract_words_file(
        archive, layout.corpus_dir(CORPUS), CORPUS
    )
    counts = frequency.count_words(words)
    assert counts.folded == EXPECTED_TABLE
    assert counts.case_sensitive == {
        "casa": 50,
        "Casa": 7,
        "città": 10,  # NFC and NFD spellings merge
        "citta": 4,
        "dell": 6,
        "perché": 3,
        "ŽLUŤOUČKÝ": 2,
    }
    # Every valid count (punctuation, URLs and numbers included) ...
    assert counts.total_tokens == 200
    # ... and only the kept tokens.
    assert counts.word_tokens == sum(EXPECTED_TABLE.values()) == 82
    assert frequency.parse_words_file(words) == counts.folded


def test_count_words_with_a_split(tmp_path: Path) -> None:
    table = tmp_path / "w.txt"
    table.write_text("1\ta+B\t3\n2\tb\t2\n3\t+\t1\n", encoding="utf-8")

    def split(token: str) -> list[str]:
        return token.split("+")

    counts = frequency.count_words(table, split=split)
    assert counts.case_sensitive == {"a": 3, "B": 3, "b": 2}
    assert counts.folded == {"a": 3, "b": 5}
    assert counts.total_tokens == 3 * 2 + 2 + 1 * 2
    assert counts.word_tokens == 8


def test_clean_token_can_keep_case() -> None:
    assert frequency.clean_token(" «Città» ", fold_case=False) == "Città"
    assert frequency.clean_token(" «Città» ") == "città"
    assert frequency.clean_token("http://X.it", fold_case=False) == ""
