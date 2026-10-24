"""Leipzig Corpora Collection word frequencies for the language plugins.

A language-neutral version of ``scripts/download_leipzig.py`` (build a
``word,freq`` table from a corpus archive) plus the lookup half of
``scripts/attach_frequencies.py`` (lemma and paradigm-summed counts).

Layout
------
Each language keeps its corpora under ``data/<language>/leipzig/``
(:func:`language_layout`)::

    raw/<corpus_id>.tar.gz          archive, unless one is supplied
    corpora/<corpus_id>/...         the words table and meta file only
    freq/<corpus_id>_freq.csv       ``word,freq``

and **never** in ``data/leipzig/freq/``: ``attach_frequencies.py``
attaches every ``*_freq.csv`` it finds there to the Romanian lexicon,
so a stray Italian table would add a column to the Romanian output.
:class:`LeipzigLayout` refuses a ``data/leipzig`` root.

Building a table
----------------
1. Use the supplied archive (for example the copy registered in
   :mod:`lexpipe.sources`) after checking its size, or download
   ``{corpus_id}.tar.gz`` from :data:`LEIPZIG_BASE_URL`.
2. Extract only ``*-words*.txt`` and ``*-meta.txt`` with tarfile's
   ``filter='data'`` policy (PEP 706: no absolute paths, no ``..``,
   no links out of the destination). The sentence and co-occurrence
   files (about 600 MB per corpus) are never needed. Extraction goes
   to a staging directory that is renamed into place only when it
   holds a words table, so an interrupted run is not reused.
3. Clean tokens exactly as ``download_leipzig.py`` does and sum the
   counts of tokens that clean to the same key.
4. Write ``<corpus_id>_freq.csv`` as that script does: header
   ``word,freq``, rows by descending frequency then word, csv-module
   default dialect.

Cleaning and lookup keys
------------------------
A token is NFC-normalised and stripped of surrounding whitespace and
of the punctuation in :data:`TOKEN_TRIM`; it is dropped if it is then
empty, contains whitespace, looks like a URL or e-mail address, or
has no letter; otherwise its lowercase form is the key. Apostrophes
are trimmed like the Romanian corpus's, so Italian elisions such as
``dell'`` count as ``dell``; this does not affect noun and adjective
lemmas. Lookup keys (:func:`normalize_key`) are NFC + strip +
lowercase. Accents and diacritics are kept (``città`` ≠ ``citta``,
``být`` ≠ ``byt``), and there is no Romanian cedilla folding.

Case-sensitive counts and totals
--------------------------------
Leipzig words tables are case-sensitive (``Medici`` and ``medici``
are separate rows). :func:`count_words` reads one table in a single
pass into :class:`WordCounts`: the cleaned tokens with their case
kept (``clean_token(..., fold_case=False)``), the folded table that
:func:`parse_words_file` returns, and two totals: ``total_tokens``,
every count in the table (punctuation and numbers included), and
``word_tokens``, the counts of kept tokens only. An optional
*split* function first breaks a token into several (Italian elided
clitics, :mod:`lexpipe.attest`); every piece counts as a token.
"""

from __future__ import annotations

import csv
import os
import re
import shutil
import tarfile
import unicodedata
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from fnmatch import fnmatch
from pathlib import Path, PurePosixPath
from typing import Final

from lexpipe import PROJECT_ROOT
from lexpipe.sources import HttpSession, download

LEIPZIG_BASE_URL: Final[str] = (
    "https://downloads.wortschatz-leipzig.de/corpora"
)
FREQ_FIELDS: Final[tuple[str, str]] = ("word", "freq")

# Same characters as scripts/download_leipzig.py's _TOKEN_TRIM.
TOKEN_TRIM: Final[str] = ".,;:!?\"“”„'«»()[]{}<>"
_URL_MARKERS: Final[tuple[str, ...]] = ("http", "www.", "://", "@")

# Preference order for the words table inside a corpus directory.
_WORDS_SUFFIXES: Final[tuple[str, ...]] = ("-words_pos_base.txt", "-words.txt")
# Archive members worth extracting (matched against the base name).
_EXTRACT_PATTERNS: Final[tuple[str, ...]] = ("*-words*.txt", "*-meta.txt")

_CORPUS_ID_RE: Final[re.Pattern[str]] = re.compile(r"[A-Za-z0-9][\w-]*")
_LANGUAGE_RE: Final[re.Pattern[str]] = re.compile(r"[a-z][a-z_]*")


def _check_corpus_id(corpus_id: str) -> str:
    if not _CORPUS_ID_RE.fullmatch(corpus_id):
        raise ValueError(f"invalid Leipzig corpus id: {corpus_id!r}")
    return corpus_id


def _is_romanian_location(root: Path) -> bool:
    resolved = root.resolve()
    return resolved.name == "leipzig" and resolved.parent.name == "data"


@dataclass(frozen=True, slots=True)
class LeipzigLayout:
    """Directories holding one language's Leipzig corpora.

    Refuses ``data/leipzig`` itself, the Romanian location (see the
    module docstring).
    """

    root: Path

    def __post_init__(self) -> None:
        if _is_romanian_location(self.root):
            raise ValueError(
                f"{self.root} is the Romanian Leipzig directory: "
                "scripts/attach_frequencies.py attaches every *_freq.csv "
                "in data/leipzig/freq/ to the Romanian lexicon. Use "
                "data/<language>/leipzig/ (language_layout())."
            )

    @property
    def raw_dir(self) -> Path:
        return self.root / "raw"

    @property
    def corpora_dir(self) -> Path:
        return self.root / "corpora"

    @property
    def freq_dir(self) -> Path:
        return self.root / "freq"

    def archive_path(self, corpus_id: str) -> Path:
        return self.raw_dir / f"{_check_corpus_id(corpus_id)}.tar.gz"

    def corpus_dir(self, corpus_id: str) -> Path:
        return self.corpora_dir / _check_corpus_id(corpus_id)

    def freq_csv(self, corpus_id: str) -> Path:
        return self.freq_dir / f"{_check_corpus_id(corpus_id)}_freq.csv"


def language_layout(
    language: str, project_root: Path = PROJECT_ROOT
) -> LeipzigLayout:
    """Layout under ``data/<language>/leipzig/`` (e.g. ``italian``)."""
    if not _LANGUAGE_RE.fullmatch(language):
        raise ValueError(f"invalid language directory name: {language!r}")
    return LeipzigLayout(project_root / "data" / language / "leipzig")


def corpus_url(corpus_id: str) -> str:
    """Download URL of a Leipzig corpus archive."""
    return f"{LEIPZIG_BASE_URL}/{_check_corpus_id(corpus_id)}.tar.gz"


def normalize_key(word: str) -> str:
    """Lookup key: NFC, surrounding whitespace removed, lowercase."""
    return unicodedata.normalize("NFC", word).strip().lower()


def clean_token(token: str, *, fold_case: bool = True) -> str:
    """Key for a Leipzig token, or "" if the token is dropped.

    Same rules as ``scripts/download_leipzig.py``'s ``clean_token``
    (see the module docstring); with ``fold_case=False`` the cleaned
    token keeps its case.
    """
    t = unicodedata.normalize("NFC", token).strip()
    t = t.strip(TOKEN_TRIM)
    if not t:
        return ""
    if any(ch.isspace() for ch in t):
        return ""
    lower = t.lower()
    if any(marker in lower for marker in _URL_MARKERS):
        return ""
    if not any(ch.isalpha() for ch in t):
        return ""
    return lower if fold_case else t


def _is_count(text: str) -> bool:
    # str.isdigit() alone accepts digits like "²" that int() rejects.
    return text.isascii() and text.isdigit()


@dataclass(frozen=True, slots=True)
class WordCounts:
    """One words table read by :func:`count_words`.

    ``case_sensitive`` is keyed by the cleaned token with its case,
    ``folded`` by its lowercase form (as :func:`parse_words_file`).
    ``total_tokens`` sums every valid count in the table, punctuation
    and numbers included; ``word_tokens`` only those of kept tokens.
    Both include the extra tokens a *split* function creates.
    """

    case_sensitive: dict[str, int]
    folded: dict[str, int]
    total_tokens: int
    word_tokens: int


def _no_split(token: str) -> Sequence[str]:
    return (token,)


def _rows(path: Path) -> Iterable[tuple[str, int]]:
    """``(token, count)`` of the valid rows of a Leipzig words table."""
    with path.open("r", encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            if lineno == 1 and not _is_count(parts[0]):
                continue
            if not _is_count(parts[2]):
                continue
            freq = int(parts[2])
            if freq > 0:
                yield parts[1], freq


def count_words(
    path: Path, *, split: Callable[[str], Sequence[str]] | None = None
) -> WordCounts:
    """Case-sensitive and folded counts of a Leipzig words table.

    Rows are filtered as in :func:`parse_words_file`. *split* maps a
    raw token to the tokens it stands for (default: itself); each
    piece is then cleaned and counted on its own.
    """
    splitter = split or _no_split
    cased: dict[str, int] = {}
    folded: dict[str, int] = {}
    total = words = 0
    for raw, freq in _rows(path):
        pieces = splitter(raw)
        total += freq * len(pieces)
        for piece in pieces:
            token = clean_token(piece, fold_case=False)
            if not token:
                continue
            words += freq
            cased[token] = cased.get(token, 0) + freq
            key = token.lower()
            folded[key] = folded.get(key, 0) + freq
    return WordCounts(cased, folded, total, words)


def parse_words_file(path: Path) -> dict[str, int]:
    """Read a Leipzig ``id<TAB>word<TAB>freq`` table into ``key → freq``.

    Skips short rows, a non-numeric header on line 1, non-numeric or
    non-positive counts and tokens :func:`clean_token` drops; counts
    of tokens with the same key are summed.
    """
    return count_words(path).folded


def write_freq_csv(freq_map: Mapping[str, int], out_path: Path) -> Path:
    """Write ``word,freq`` rows by descending frequency, then word.

    Written to a temporary file and renamed, so readers never see a
    partial table.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_name(out_path.name + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(FREQ_FIELDS)
        for word, freq in sorted(
            freq_map.items(), key=lambda kv: (-kv[1], kv[0])
        ):
            writer.writerow([word, freq])
    os.replace(tmp, out_path)
    return out_path


def ensure_archive(
    corpus_id: str,
    layout: LeipzigLayout,
    *,
    archive: Path | None = None,
    expected_size: int | None = None,
    session: HttpSession | None = None,
    log: Callable[[str], None] = print,
) -> Path:
    """Return a complete corpus archive, downloading it if needed.

    With *archive*, that file is used and never downloaded over; with
    *expected_size*, a file of any other size (for example one still
    being downloaded) raises :class:`ValueError`. Without *archive*,
    ``layout.raw_dir/<corpus_id>.tar.gz`` is used or downloaded.
    """
    if archive is not None:
        if not archive.is_file():
            raise FileNotFoundError(f"corpus archive not found: {archive}")
        _check_size(archive, expected_size)
        return archive
    target = layout.archive_path(corpus_id)
    if target.exists():
        _check_size(target, expected_size)
        return target
    log(f"  downloading {corpus_url(corpus_id)}")
    return download(
        corpus_url(corpus_id),
        target,
        expected_size=expected_size,
        session=session,
        log=log,
    )


def _check_size(path: Path, expected_size: int | None) -> None:
    if expected_size is None:
        return
    size = path.stat().st_size
    if size != expected_size:
        raise ValueError(
            f"{path}: {size:,} bytes, expected {expected_size:,} "
            "(incomplete download, or a different release)"
        )


def find_words_file(corpus_dir: Path, corpus_id: str) -> Path:
    """Locate the words table in an extracted corpus directory.

    Archives unpack as ``<corpus_id>/<corpus_id>-words.txt``; a
    missing inner directory is accepted. A ``-words_pos_base.txt``
    table wins over ``-words.txt``; otherwise any ``*-words*.txt``,
    shallowest first.
    """
    inner = corpus_dir / corpus_id
    base = inner if inner.is_dir() else corpus_dir
    for suffix in _WORDS_SUFFIXES:
        candidate = base / f"{corpus_id}{suffix}"
        if candidate.is_file():
            return candidate
    patterns = [f"{corpus_id}{suffix}" for suffix in _WORDS_SUFFIXES]
    patterns.extend(f"*{suffix}" for suffix in _WORDS_SUFFIXES)
    for pattern in patterns:
        matches = [p for p in base.rglob(pattern) if p.is_file()]
        if matches:
            return min(
                matches,
                key=lambda p: (len(p.relative_to(base).parts), str(p)),
            )
    raise FileNotFoundError(
        f"no words or words_pos_base table for {corpus_id} in {corpus_dir}"
    )


def extract_words_file(
    archive: Path, corpus_dir: Path, corpus_id: str
) -> Path:
    """Extract the words table (and meta file) of *archive*; return it.

    An existing words table under *corpus_dir* is reused. Members are
    extracted with ``filter='data'``; an archive member that would
    escape the destination raises :class:`tarfile.FilterError`.
    """
    if corpus_dir.is_dir():
        try:
            return find_words_file(corpus_dir, corpus_id)
        except FileNotFoundError:
            pass
    staging = corpus_dir.with_name(corpus_dir.name + ".extracting")
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True)
    try:
        with tarfile.open(archive, "r:*") as tf:
            for member in tf:
                name = PurePosixPath(member.name).name
                if member.isfile() and any(
                    fnmatch(name, pattern) for pattern in _EXTRACT_PATTERNS
                ):
                    tf.extract(member, staging, filter="data")
        words = find_words_file(staging, corpus_id)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    try:
        os.replace(staging, corpus_dir)
    except OSError as exc:
        shutil.rmtree(staging, ignore_errors=True)
        raise FileExistsError(
            f"{corpus_dir} exists without a words table; remove it and "
            "re-run"
        ) from exc
    return corpus_dir / words.relative_to(staging)


def build_freq_csv(
    corpus_id: str,
    layout: LeipzigLayout,
    *,
    archive: Path | None = None,
    expected_size: int | None = None,
    session: HttpSession | None = None,
    log: Callable[[str], None] = print,
) -> Path:
    """Archive → ``layout.freq_dir/<corpus_id>_freq.csv``; return its path.

    The table is rebuilt on every call (as in the Romanian script);
    downloading and extraction are skipped when already done. An
    archive whose words table yields no usable row raises
    :class:`ValueError` rather than writing an empty table.
    """
    source = ensure_archive(
        corpus_id,
        layout,
        archive=archive,
        expected_size=expected_size,
        session=session,
        log=log,
    )
    words = extract_words_file(source, layout.corpus_dir(corpus_id), corpus_id)
    log(f"  parsing {words}")
    freq_map = parse_words_file(words)
    if not freq_map:
        raise ValueError(f"{words}: no usable rows")
    out = write_freq_csv(freq_map, layout.freq_csv(corpus_id))
    log(f"  wrote {len(freq_map):,} word types to {out}")
    return out


def load_freq_table(path: Path) -> dict[str, int]:
    """Load a ``word,freq`` CSV into ``normalize_key(word) → freq``.

    Keys that normalise alike are summed; rows with a blank word, a
    non-integer or a non-positive count are skipped (as in
    ``scripts/attach_frequencies.py``).
    """
    table: dict[str, int] = {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            key = normalize_key(row.get("word") or "")
            if not key:
                continue
            try:
                freq = int((row.get("freq") or "").strip())
            except ValueError:
                continue
            if freq > 0:
                table[key] = table.get(key, 0) + freq
    return table


def lemma_frequency(lemma: str, table: Mapping[str, int]) -> int:
    """Corpus count of the lemma's own form (0 if unattested)."""
    return table.get(normalize_key(lemma), 0)


def paradigm_frequency(forms: Iterable[str], table: Mapping[str, int]) -> int:
    """Summed count of the *distinct* forms of a paradigm.

    Forms are compared after :func:`normalize_key`, so a form filling
    several cells (Czech GEN.SG = ACC.SG *biologa*) counts once.
    Include the lemma when it is itself a paradigm cell, as it is for
    nouns (NOM.SG) and adjectives (MASC.SG).
    """
    keys = {normalize_key(form) for form in forms}
    keys.discard("")
    return sum(table.get(key, 0) for key in keys)
