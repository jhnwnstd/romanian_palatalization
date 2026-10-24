"""The Leipzig corpus of the Italian pipeline (R2, R11, R24).

``ita_news_2024_1M`` (the registry row's ``leipzig_id``) is read from
the pinned archive ``data/sources/leipzig/ita_news_2024_1M.tar.gz``;
only its words table is extracted, into
``data/italian/leipzig/corpora/`` (:mod:`lexpipe.frequency`), never into
``data/leipzig/`` (R2). The archive's size is checked against the
source registry unless the run context's parameters set
``verify_source_size`` to False (test trees use a small synthetic
archive).

Counts come from :class:`lexpipe.attest.AttestationTable` with the
Italian elision clitics split off (l'amico → l' + amico, R11), both
case-sensitive and folded.

Attestation (:func:`attested_flags`): the τ thresholds are tested on
the folded count (design §11, for comparability with Romanian), and a
form needs at least one token in its own case — its lowercase spelling
for a common noun — to count as attested at any corpus threshold.
Tokens of a form seen only capitalized are surnames, firm and place
names (Bocchi, Monte dei Paschi), sentence-initial words at most; R25
keeps names out. A generated paradigm (a derivative's forms) needs
one such token among its forms (:func:`paradigm_flags`).
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Final

from lexpipe import attest, frequency
from lexpipe.attest import AttestationTable
from lexpipe.context import RunContext
from lexpipe.langs.italian.layout import Layout

VERIFY_SIZE_PARAM: Final[str] = "verify_source_size"


def expected_size(ctx: RunContext, layout: Layout) -> int | None:
    """The registered archive size, or None when checks are off."""
    if not ctx.params.get(VERIFY_SIZE_PARAM, True):
        return None
    return layout.leipzig_source.size


def words_file(ctx: RunContext, layout: Layout) -> Path:
    """The extracted words table, extracting it once from the archive."""
    corpus_id = ctx.language.leipzig_id
    corpus_dir = layout.leipzig.corpus_dir(corpus_id)
    ctx.check_writable(corpus_dir)
    archive = frequency.ensure_archive(
        corpus_id,
        layout.leipzig,
        archive=layout.leipzig_archive,
        expected_size=expected_size(ctx, layout),
        log=lambda message: print(message),
    )
    return frequency.extract_words_file(archive, corpus_dir, corpus_id)


def attestation_table(ctx: RunContext, words: Path) -> AttestationTable:
    """Counts of *words* with the language's elision clitics split."""
    return AttestationTable.from_words_file(
        words, clitics=ctx.language.elision_clitics
    )


def _names_filtered(flags: dict[str, bool], own_case: int) -> dict[str, bool]:
    if own_case > 0:
        return flags
    return {tau: tau == "dictionary" for tau in flags}


def attested_flags(table: AttestationTable, text: str) -> dict[str, bool]:
    """``{τ: attested}`` of one form (module docstring): the folded
    count reaches τ and the form has a token in its own case."""
    return _names_filtered(
        table.attested_by_tau(text), table.count(text, case_sensitive=True)
    )


def paradigm_flags(
    table: AttestationTable, forms: Iterable[str]
) -> dict[str, bool]:
    """``{τ: attested}`` of a generated paradigm: the folded counts of
    its distinct forms summed (:func:`lexpipe.attest.attest_paradigm`),
    with a token in its own case among them (module docstring)."""
    listed = list(dict.fromkeys(forms))
    own_case = attest.attest_paradigm(table, listed, case_sensitive=True)
    return _names_filtered(
        dict(attest.attest_paradigm(table, listed).attested_by_tau),
        own_case.total,
    )
