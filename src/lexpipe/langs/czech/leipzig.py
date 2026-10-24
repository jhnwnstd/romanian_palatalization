"""Leipzig ``ces_news_2024_1M`` and the CNK SYN2015 lemma list.

Leipzig (Goldhahn et al. 2012) is the Czech attestation corpus (R11).
Its archive is the pinned source ``leipzig_ces_news_2024_1M``; only the
words table and meta file are extracted, into
``data/czech/leipzig/corpora/`` (:mod:`lexpipe.frequency`; never into
``data/leipzig/``, R2). Czech needs no clitic splitting. The qc stage
reads the vocabulary (to restrict the MorfFlex projection), the attest
stage the counts; both call :func:`ensure_words_file`, which extracts
the table once and reuses it afterwards.

Leipzig-derived numbers are not tracked until a person has checked the
Leipzig terms (R9): they stay in ``data/czech/``.

The CNK SYN2015 lemma list (Křen et al. 2015; CC BY 4.0) is a cross-check
of lemma frequencies: a zip holding one TSV with eight columns (rank,
lemma, frequency, normalized frequency, three genre frequencies,
characteristic) for lemmas with frequency ≥ 10. Lemmas keep their case
(Praha, Čech) and homographs of different parts of speech are merged.
"""

from __future__ import annotations

import zipfile
from pathlib import Path
from typing import Final

from lexpipe import sources
from lexpipe.attest import AttestationTable
from lexpipe.context import RunContext
from lexpipe.frequency import (
    build_freq_csv,
    count_words,
    extract_words_file,
    find_words_file,
)
from lexpipe.langs.czech.layout import (
    LEIPZIG_SOURCE,
    CzechPaths,
    check_source_sizes,
)

CNK_LEMMA_COLUMN: Final[int] = 1
CNK_FREQ_COLUMN: Final[int] = 2


def _expected_size(ctx: RunContext) -> int | None:
    if not check_source_sizes(ctx):
        return None
    return sources.get_source(LEIPZIG_SOURCE).size


def ensure_words_file(ctx: RunContext, paths: CzechPaths) -> Path:
    """The extracted words table, extracting it from the archive once."""
    corpus_id = ctx.language.leipzig_id
    corpus_dir = paths.leipzig.corpus_dir(corpus_id)
    try:
        return find_words_file(corpus_dir, corpus_id)
    except FileNotFoundError:
        pass
    archive = paths.leipzig_archive
    if not archive.is_file():
        raise FileNotFoundError(
            f"Leipzig archive missing: {archive} (run the fetch stage)"
        )
    expected = _expected_size(ctx)
    if expected is not None and archive.stat().st_size != expected:
        raise ValueError(
            f"{archive}: {archive.stat().st_size:,} bytes, expected "
            f"{expected:,} (incomplete download or another release)"
        )
    ctx.check_writable(corpus_dir)
    return extract_words_file(archive, corpus_dir, corpus_id)


def ensure_freq_csv(ctx: RunContext, paths: CzechPaths) -> Path:
    """``data/czech/leipzig/freq/<corpus>_freq.csv`` (rebuilt each call)."""
    corpus_id = ctx.language.leipzig_id
    target = paths.leipzig.freq_csv(corpus_id)
    ctx.check_writable(target)
    return build_freq_csv(
        corpus_id,
        paths.leipzig,
        archive=paths.leipzig_archive,
        expected_size=_expected_size(ctx),
        log=lambda message: print(message, flush=True),
    )


def vocabulary(words_file: Path) -> frozenset[str]:
    """Case-folded word types of a Leipzig words table."""
    return frozenset(count_words(words_file).folded)


def attestation_table(words_file: Path) -> AttestationTable:
    """Counts of the Czech corpus; Czech has no elided clitics."""
    return AttestationTable.from_words_file(words_file)


def load_cnk_lemmas(path: Path) -> dict[str, int]:
    """``lemma -> frequency`` from the SYN2015 lemma-list zip.

    Rows whose frequency is not an integer are skipped; a lemma listed
    twice (it should not be) is summed.
    """
    table: dict[str, int] = {}
    with zipfile.ZipFile(path) as archive:
        members = sorted(
            name for name in archive.namelist() if not name.endswith("/")
        )
        if not members:
            raise ValueError(f"{path}: empty zip")
        with archive.open(members[0]) as fh:
            for raw in fh:
                parts = raw.decode("utf-8").rstrip("\r\n").split("\t")
                if len(parts) <= CNK_FREQ_COLUMN:
                    continue
                lemma = parts[CNK_LEMMA_COLUMN].strip()
                count = parts[CNK_FREQ_COLUMN].strip()
                if lemma and count.isascii() and count.isdigit():
                    table[lemma] = table.get(lemma, 0) + int(count)
    return table
