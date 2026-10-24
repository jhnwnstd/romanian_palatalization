"""Pinned external resources: where they come from and how to check them.

Every file the Italian and Czech pipelines read from outside the
repository is registered here with its URL, its local path under
``data/sources/``, its exact size, licence and citation. All URLs and
sizes were verified on 2026-10-23.

Two kinds of entry (:class:`SourceKind`):

  - ``download``: a file fetched from its URL;
  - ``derived``: a file made locally from another entry
    (``derived_from``) by a fixed procedure (``recipe``), registered
    so that its provenance is pinned like a download's. The English
    Wiktionary subsets ``data/sources/kaikki/en_wiktionary_{it,cs}.jsonl``
    (``kaikki_raw_it``, ``kaikki_raw_cs``) are the entries of
    ``kaikki_raw`` whose top-level ``lang_code`` is the language's,
    copied byte for byte by :func:`lexpipe.kaikki.split_raw_dump`.
    A derived entry pins its size and SHA-256 in the registry, and the
    SHA-256 of the parent it was made from (``derived_from_sha256``).

Three operations:

  - :func:`fetch` downloads a missing file (resumably, through a
    ``.part`` file) and leaves a file that already has the expected
    size alone. A file with the *wrong* size is never touched: it may
    be another process's download in progress. A derived file is
    never downloaded: when it is missing, :func:`fetch` raises
    :class:`DerivedSourceMissing`, whose message gives the command
    that makes it.
  - :func:`verify` checks the size and SHA-256 of a local file. A
    source whose release hash is known is pinned in the registry
    (``Source.sha256``). Otherwise the first successful verification
    records the hash in ``docs/sources_lock.tsv`` (tracked), and
    every later one must match it, so a silently replaced upstream
    file is caught. A derived file whose parent is now locked at
    another SHA-256 than the one it was made from is ``stale``: the
    parent moved on and the file must be made again.
  - :func:`head` asks the server (one HEAD request) for the ETag,
    Last-Modified and size of a **rolling** source: kaikki.org
    regenerates its files in place at unversioned URLs. The
    ``verify-sources`` stage passes the answer to :func:`verify`,
    which records it in the lock next to the hash, and later reports
    when the server has moved on to another extraction. Only the
    stage does network I/O; :func:`verify` itself never does.

Lock file
---------
``docs/sources_lock.tsv`` has the columns :data:`LOCK_FIELDS`:
``name path bytes sha256 recorded etag last_modified``. ``recorded``
is the date of the first verification; ``etag``/``last_modified``
are what the server said when the local file was first verified
against it (empty for non-rolling sources). A lock written before
the last two columns existed is still read.

Contracts
---------
- Registry paths are relative POSIX paths under ``data/sources/``,
  resolved against a root (the repository by default; tests pass a
  temporary directory).
- Downloads ask for ``Accept-Encoding: identity`` so byte ranges and
  sizes refer to the file itself, never to a compressed transfer.
- The User-Agent names the project and its public repository as the
  contact. Set ``LEXPIPE_CONTACT`` to append another contact (an
  e-mail address, say) when a server's policy asks for one.
- Downloaded files are untrusted data. Nothing here opens them
  beyond hashing.

Licensing: see :mod:`lexpipe.licensing` for what may be redistributed
in the repository. Wiktionary-derived data is CC BY-SA; MorfFlex,
DeriNet and czetyl are CC BY-NC-SA, so columns derived from them stay
in gitignored files (``docs/multilingual_design.md`` R9).
"""

from __future__ import annotations

import contextlib
import csv
import dataclasses
import hashlib
import os
import re
import time
from collections.abc import Callable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from datetime import date
from enum import StrEnum
from pathlib import Path, PurePosixPath
from types import MappingProxyType
from typing import Any, Final, Protocol, TypeAlias

import requests

from lexpipe import DOCS_DIR, PROJECT_ROOT

LOCK_PATH: Final[Path] = DOCS_DIR / "sources_lock.tsv"
LOCK_FIELDS: Final[tuple[str, ...]] = (
    "name",
    "path",
    "bytes",
    "sha256",
    "recorded",
    "etag",
    "last_modified",
)
# Header of locks written before the upstream columns existed.
_LEGACY_LOCK_FIELDS: Final[tuple[str, ...]] = LOCK_FIELDS[:5]
# Plain TSV: ETags are quoted strings ("6ad8a097-b7b68513") and are
# kept verbatim; a field holding a tab or newline raises csv.Error.
_LOCK_DIALECT: Final[Mapping[str, Any]] = MappingProxyType(
    {"delimiter": "\t", "quoting": csv.QUOTE_NONE, "quotechar": None}
)
_SHA256_RE: Final[re.Pattern[str]] = re.compile(r"[0-9a-f]{64}")

CONTACT_URL: Final[str] = "https://github.com/jhnwnstd/romanian_palatalization"
CONTACT_ENV: Final[str] = "LEXPIPE_CONTACT"

# (connect, read) seconds. The read timeout bounds a stalled transfer;
# a re-run resumes from the .part file.
HTTP_TIMEOUT: Final[tuple[float, float]] = (30.0, 120.0)
PROGRESS_INTERVAL_S: Final[float] = 15.0
_CHUNK: Final[int] = 1 << 20

# Citations checked against the ACL Anthology, LINDAT and korpus.cz
# records on 2026-10-23.
_WIKTEXTRACT_CITATION: Final[str] = (
    "Tatu Ylonen. 2022. Wiktextract: Wiktionary as Machine-Readable "
    "Structured Data. In Proceedings of the Thirteenth Language "
    "Resources and Evaluation Conference (LREC 2022), pages "
    "1317-1325, Marseille. ELRA."
)
_LEIPZIG_CITATION: Final[str] = (
    "Dirk Goldhahn, Thomas Eckart and Uwe Quasthoff. 2012. Building "
    "Large Monolingual Dictionaries at the Leipzig Corpora Collection: "
    "From 100 to 200 Languages. In Proceedings of the Eighth "
    "International Conference on Language Resources and Evaluation "
    "(LREC'12), pages 759-765, Istanbul. ELRA."
)
# Nobody has yet read Leipzig's download terms first-hand (the site
# sits behind a bot check); lexpipe.licensing keeps Leipzig-derived
# numbers out of tracked files until someone does.
_LEIPZIG_LICENCE: Final[str] = (
    "terms unverified: Leipzig Corpora Collection download terms, "
    "to be checked by hand"
)
_CNK_CITATION: Final[str] = (
    "Český národní korpus: Srovnávací frekvenční seznamy. Ústav "
    "Českého národního korpusu FF UK, Praha 2016. Lists computed from "
    "SYN2015: Křen, M., Cvrček, V., Čapka, T., Čermáková, A., "
    "Hnátková, M., Chlumská, L., Jelínek, T., Kováříková, D., "
    "Petkevič, V., Procházka, P., Skoumalová, H., Škrabal, M., "
    "Truneček, P., Vondřička, P. and Zasina, A. SYN2015: "
    "reprezentativní korpus psané češtiny. Ústav Českého národního "
    "korpusu FF UK, Praha 2015."
)
_LINDAT: Final[str] = (
    "https://lindat.mff.cuni.cz/repository/server/api/core/bitstreams"
)
_KAIKKI_RAW_PATH: Final[str] = (
    "data/sources/kaikki/raw-wiktextract-data.jsonl.gz"
)
# The kaikki_raw release (docs/sources_lock.tsv, 2026-10-23) the
# English Wiktionary subsets were split from.
_KAIKKI_RAW_SHA256: Final[str] = (
    "244dcdc56868d81dc60e3820657780c06fac2729a90e2c5ffbfa13e6f9ecad1d"
)


def _split_recipe(code: str) -> str:
    """The command (run from the repository root) that splits the
    *code* subset out of kaikki_raw."""
    return (
        "python -c \"import sys; sys.path.insert(0, 'src'); "
        "from lexpipe import kaikki; "
        f"print(kaikki.split_raw_dump('{_KAIKKI_RAW_PATH}', "
        f"'data/sources/kaikki', ['{code}']))\""
    )


class SourceKind(StrEnum):
    """How a registered file comes to exist (module docstring)."""

    DOWNLOAD = "download"
    DERIVED = "derived"


@dataclass(frozen=True, slots=True)
class Source:
    """One pinned external file.

    ``path`` is relative to the repository root and always lies under
    ``data/sources/``; ``size`` is the exact byte count of the release
    the pipeline was built against. ``languages`` lists the plugin
    codes whose ``fetch``/``verify-sources`` stages include it.
    ``sha256`` pins the release in the registry itself when its hash
    is known in advance (empty: the lock file pins it). ``rolling``
    marks a URL that is regenerated in place (kaikki.org), whose
    ETag and Last-Modified ``verify-sources`` records.

    A ``derived`` entry (``kind``) has no URL; it names the entry it is
    made from (``derived_from``), the SHA-256 of the parent file it
    was made from (``derived_from_sha256``), its own ``sha256``, and
    the shell command, run from the repository root, that makes it
    (``recipe``).
    """

    name: str
    url: str
    path: str
    size: int
    licence: str
    citation: str
    notes: str = ""
    languages: tuple[str, ...] = ()
    sha256: str = ""
    rolling: bool = False
    kind: SourceKind = SourceKind.DOWNLOAD
    derived_from: str = ""
    derived_from_sha256: str = ""
    recipe: str = ""

    def __post_init__(self) -> None:
        rel = PurePosixPath(self.path)
        if (
            rel.is_absolute()
            or ".." in rel.parts
            or rel.parts[:2] != ("data", "sources")
            or len(rel.parts) < 3
        ):
            raise ValueError(
                f"{self.name}: path must lie under data/sources/: "
                f"{self.path!r}"
            )
        if self.size <= 0:
            raise ValueError(f"{self.name}: size must be positive")
        object.__setattr__(self, "kind", SourceKind(self.kind))
        for name in ("sha256", "derived_from_sha256"):
            value = getattr(self, name)
            if value and not _SHA256_RE.fullmatch(value):
                raise ValueError(f"{self.name}: {name} must be 64 hex digits")
        if self.kind is SourceKind.DOWNLOAD:
            if not self.url.startswith("https://"):
                raise ValueError(f"{self.name}: url must be https")
            if self.derived_from or self.derived_from_sha256 or self.recipe:
                raise ValueError(
                    f"{self.name}: only derived sources have derived_from, "
                    "derived_from_sha256 and recipe"
                )
            return
        if self.url or self.rolling:
            raise ValueError(
                f"{self.name}: a derived source has no URL and is not rolling"
            )
        missing = [
            name
            for name in ("derived_from", "derived_from_sha256", "sha256")
            if not getattr(self, name)
        ]
        if missing or not self.recipe.strip():
            raise ValueError(
                f"{self.name}: a derived source needs derived_from, "
                "derived_from_sha256, sha256 and recipe"
            )
        if self.derived_from == self.name:
            raise ValueError(
                f"{self.name}: a source cannot derive from itself"
            )

    @property
    def derived(self) -> bool:
        return self.kind is SourceKind.DERIVED


_REGISTRY: Final[tuple[Source, ...]] = (
    Source(
        name="kaikki_raw",
        url="https://kaikki.org/dictionary/raw-wiktextract-data.jsonl.gz",
        path=_KAIKKI_RAW_PATH,
        size=3_082_192_147,
        licence="CC BY-SA 4.0 + GFDL (Wiktionary content)",
        citation=_WIKTEXTRACT_CITATION,
        notes=(
            "All-language wiktextract dump of English Wiktionary, built "
            "from the 2026-10-01 enwiktionary dump. kaikki.org "
            "regenerates it in place, so the lock also records the "
            "server's ETag and Last-Modified. Split it per language "
            "with lexpipe.kaikki.split_raw_dump (top-level lang_code): "
            "kaikki_raw_it and kaikki_raw_cs."
        ),
        languages=("it", "cs"),
        rolling=True,
    ),
    Source(
        name="kaikki_raw_it",
        url="",
        path="data/sources/kaikki/en_wiktionary_it.jsonl",
        size=665_207_473,
        licence="CC BY-SA 4.0 + GFDL (Wiktionary content)",
        citation=_WIKTEXTRACT_CITATION,
        notes=(
            "The 624,122 entries of kaikki_raw whose top-level lang_code "
            "is 'it', copied byte for byte by "
            "lexpipe.kaikki.split_raw_dump: the Italian primary lexicon. "
            "A fresh split of the pinned kaikki_raw reproduced this file "
            "on 2026-10-23. Make it again, and re-pin it here, whenever "
            "kaikki_raw moves (docs/runbooks.md, section 2)."
        ),
        languages=("it",),
        sha256=(
            "ae7de484ba53e3a6a472b9672a69d7e7fdc96cf55b68c5399a5b338adb387f0b"
        ),
        kind=SourceKind.DERIVED,
        derived_from="kaikki_raw",
        derived_from_sha256=_KAIKKI_RAW_SHA256,
        recipe=_split_recipe("it"),
    ),
    Source(
        name="kaikki_raw_cs",
        url="",
        path="data/sources/kaikki/en_wiktionary_cs.jsonl",
        size=223_701_337,
        licence="CC BY-SA 4.0 + GFDL (Wiktionary content)",
        citation=_WIKTEXTRACT_CITATION,
        notes=(
            "The 72,441 entries of kaikki_raw whose top-level lang_code "
            "is 'cs', copied byte for byte by "
            "lexpipe.kaikki.split_raw_dump: the Czech primary lexicon. "
            "A fresh split of the pinned kaikki_raw reproduced this file "
            "on 2026-10-23. Make it again, and re-pin it here, whenever "
            "kaikki_raw moves (docs/runbooks.md, section 2)."
        ),
        languages=("cs",),
        sha256=(
            "674142063ba8dc24be953d56522b08ae5bffd86cbe7f15e7219c44c0fb66c225"
        ),
        kind=SourceKind.DERIVED,
        derived_from="kaikki_raw",
        derived_from_sha256=_KAIKKI_RAW_SHA256,
        recipe=_split_recipe("cs"),
    ),
    Source(
        name="kaikki_it_extract",
        url="https://kaikki.org/dictionary/downloads/it/it-extract.jsonl.gz",
        path="data/sources/kaikki/it-extract.jsonl.gz",
        size=43_849_715,
        licence="CC BY-SA 4.0 + GFDL (Wikizionario content)",
        citation=_WIKTEXTRACT_CITATION,
        notes=(
            "Italian Wiktionary (Wikizionario) extraction: hand-entered "
            "plurals and gender, the second edition used for QC."
        ),
        languages=("it",),
        rolling=True,
    ),
    Source(
        name="kaikki_cs_extract",
        url="https://kaikki.org/dictionary/downloads/cs/cs-extract.jsonl.gz",
        path="data/sources/kaikki/cs-extract.jsonl.gz",
        size=39_633_204,
        licence="CC BY-SA 4.0 + GFDL (Wikislovník content)",
        citation=_WIKTEXTRACT_CITATION,
        notes=(
            "Czech Wiktionary (Wikislovník) extraction: more nouns and "
            "labelled variants, the second edition used for QC."
        ),
        languages=("cs",),
        rolling=True,
    ),
    Source(
        name="morphit",
        url=(
            "https://raw.githubusercontent.com/giodegas/"
            "morphit-lemmatizer/master/master/morph-it_048_utf8.txt"
        ),
        path="data/sources/morphit/morph-it_048_utf8.txt",
        size=19_068_565,
        licence="CC BY-SA 2.0 / LGPL (dual)",
        citation=(
            "Eros Zanchetta and Marco Baroni. 2005. Morph-it! A free "
            "corpus-based morphological resource for the Italian "
            "language. In Proceedings of Corpus Linguistics 2005, "
            "University of Birmingham."
        ),
        notes=(
            "Morph-it! 0.4.8 full-form lexicon (UTF-8 text; gender, "
            "attested plurals): automated Italian QC. Fetched from a "
            "GitHub mirror and pinned by sha256, because the official "
            "download (https://docs.sslmit.unibo.it/lib/exe/fetch.php"
            "?media=resources:morph-it.tgz) sits behind a "
            "human-verification page that the pipeline does not "
            "bypass. The official tarball, downloaded in a browser to "
            "data/sources/morphit/morph-it.tgz (2,617,548 bytes), may "
            "replace the mirror copy."
        ),
        languages=("it",),
        sha256=(
            "7375216351ee59cee4326b455503f3a6c4e9542cd63ecb86ff5cd58ca5c8e800"
        ),
    ),
    Source(
        name="morfflex",
        url=f"{_LINDAT}/cc940076-7bf0-449b-81a4-5e10e023084e/content",
        path="data/sources/morfflex/czech-morfflex-2.1.tsv.xz",
        size=250_486_476,
        licence="CC BY-NC-SA 4.0",
        citation=(
            "Jan Hajič, Jaroslava Hlaváčová, Marie Mikulová, Milan "
            "Straka and Barbora Štěpánková. 2024. MorfFlex CZ 2.1 "
            "(2024-12-23). LINDAT/CLARIAH-CZ digital library, Charles "
            "University, http://hdl.handle.net/11234/1-5833. Reference "
            "paper: Jaroslava Hlaváčová, Marie Mikulová, Barbora "
            "Štěpánková, Milan Straka and Jan Hajič. 2026. MorfFlex: "
            "Handling Rich Morphology. In Proceedings of the Fifteenth "
            "Language Resources and Evaluation Conference (LREC 2026), "
            "pages 11495-11505."
        ),
        notes=(
            "Lemma-tag-form triples for every case x number cell, with "
            "animacy and variant marks; automated Czech QC. Variant "
            "digits are labels, not frequency ranks. NC-SA: derived "
            "columns stay in gitignored files."
        ),
        languages=("cs",),
    ),
    Source(
        name="derinet",
        url=f"{_LINDAT}/6fa2bce5-51cd-42a9-ba70-cb640542cd13/content",
        path="data/sources/derinet/derinet-2-3.tsv",
        size=439_063_096,
        licence="CC BY-NC-SA 4.0",
        citation=(
            "Michal Olbrich, Viktória Brezinová, Šárka Dohnalová, "
            "Vojtěch John, Lukáš Kyjánek, Aleš Papáček, Emil Svoboda, "
            "Magda Ševčíková, Jonáš Vidra and Zdeněk Žabokrtský. 2025. "
            "DeriNet 2.3. LINDAT/CLARIAH-CZ digital library, Charles "
            "University, http://hdl.handle.net/11234/1-5846. Format: "
            "Jonáš Vidra, Zdeněk Žabokrtský, Magda Ševčíková and Lukáš "
            "Kyjánek. 2019. DeriNet 2.0: Towards an All-in-One "
            "Word-Formation Resource. In Proceedings of the Second "
            "International Workshop on Resources and Tools for "
            "Derivational Morphology (DeriMo 2019), pages 81-89."
        ),
        notes=(
            "Czech derivational network: noun->verb/adjective edges "
            "with ConjugClass, loanword flags. NC-SA."
        ),
        languages=("cs",),
    ),
    Source(
        name="czetyl",
        url=f"{_LINDAT}/3cc28c39-b1fa-4bd6-a0e2-4fc489b3d275/content",
        path="data/sources/czetyl/czetyl.tsv",
        size=191_934,
        licence="CC BY-NC-SA 4.0",
        citation=(
            "Jiří Rejzek, Aleš Papáček, Viktória Brezinová and Zdeněk "
            "Žabokrtský. 2025. Czech Etymological Lexicon 1.0. "
            "LINDAT/CLARIAH-CZ digital library, Charles University, "
            "http://hdl.handle.net/11234/1-5845. Extracted from Jiří "
            "Rejzek, Český etymologický slovník, LEDA, 2015."
        ),
        notes=(
            "word, donor-language chain (nearest donor first), "
            "loan|native. NC-SA."
        ),
        languages=("cs",),
    ),
    Source(
        name="cnk_syn2015_lemma",
        url=(
            "https://wiki.korpus.cz/lib/exe/fetch.php/"
            "seznamy:syn2015_lemma_utf8.zip"
        ),
        path="data/sources/cnk/syn2015_lemma_utf8.zip",
        size=1_273_063,
        licence="CC BY 4.0",
        citation=_CNK_CITATION,
        notes=(
            "SYN2015 lemma frequency list (freq >= 10): cross-check "
            "for the Leipzig-derived Czech counts."
        ),
        languages=("cs",),
    ),
    Source(
        name="cnk_syn2015_word",
        url=(
            "https://wiki.korpus.cz/lib/exe/fetch.php/"
            "seznamy:syn2015_word_utf8.zip"
        ),
        path="data/sources/cnk/syn2015_word_utf8.zip",
        size=3_089_732,
        licence="CC BY 4.0",
        citation=_CNK_CITATION,
        notes="SYN2015 word-form frequency list (freq >= 10).",
        languages=("cs",),
    ),
    Source(
        name="leipzig_ita_news_2024_1M",
        url=(
            "https://downloads.wortschatz-leipzig.de/corpora/"
            "ita_news_2024_1M.tar.gz"
        ),
        path="data/sources/leipzig/ita_news_2024_1M.tar.gz",
        size=261_325_052,
        licence=_LEIPZIG_LICENCE,
        citation=_LEIPZIG_CITATION,
        notes=(
            "Italian counterpart of ron_news_2024_1M. Its frequency "
            "table goes to data/italian/leipzig/, never to "
            "data/leipzig/freq/ (see lexpipe.frequency)."
        ),
        languages=("it",),
    ),
    Source(
        name="leipzig_ces_news_2024_1M",
        url=(
            "https://downloads.wortschatz-leipzig.de/corpora/"
            "ces_news_2024_1M.tar.gz"
        ),
        path="data/sources/leipzig/ces_news_2024_1M.tar.gz",
        size=214_278_386,
        licence=_LEIPZIG_LICENCE,
        citation=_LEIPZIG_CITATION,
        notes=(
            "Czech counterpart of ron_news_2024_1M. Its frequency "
            "table goes to data/czech/leipzig/, never to "
            "data/leipzig/freq/ (see lexpipe.frequency)."
        ),
        languages=("cs",),
    ),
)


def check_registry(entries: Sequence[Source]) -> None:
    """Names are unique and a derived entry follows its parent, so
    ``fetch`` and ``verify-sources`` handle the parent first."""
    seen: set[str] = set()
    for source in entries:
        if source.name in seen:
            raise ValueError(f"source {source.name!r} registered twice")
        if source.derived and source.derived_from not in seen:
            raise ValueError(
                f"{source.name}: derived_from {source.derived_from!r} must "
                "be registered before it"
            )
        seen.add(source.name)


check_registry(_REGISTRY)
SOURCES: Final[Mapping[str, Source]] = MappingProxyType(
    {source.name: source for source in _REGISTRY}
)


class DownloadError(RuntimeError):
    """A download failed or produced a file of the wrong size."""


class DerivedSourceMissing(DownloadError):
    """A derived source is absent: :func:`fetch` never downloads it; the
    message gives the command that makes it."""


class VerifyStatus(StrEnum):
    """Outcome of :func:`verify`. Only PINNED and RECORDED are ok."""

    PINNED = "pinned"  # sha256 matches the lock file
    RECORDED = "recorded"  # first verification; sha256 now locked
    MISSING = "missing"
    SIZE_MISMATCH = "size-mismatch"
    SHA256_MISMATCH = "sha256-mismatch"  # vs the lock or the registry
    LOCK_MISMATCH = "lock-mismatch"  # lock row is for another path/size
    STALE = "stale"  # derived from another release of its parent


@dataclass(frozen=True, slots=True)
class VerifyResult:
    """What :func:`verify` found for one source.

    ``upstream_changed`` is True when a rolling source's server now
    reports another ETag, Last-Modified or size than the lock
    recorded: the local file is still the pinned one, but a fresh
    download would not be.
    """

    name: str
    path: Path
    status: VerifyStatus
    size: int | None = None
    sha256: str | None = None
    detail: str = ""
    upstream_changed: bool = False

    @property
    def ok(self) -> bool:
        return self.status in (VerifyStatus.PINNED, VerifyStatus.RECORDED)


@dataclass(frozen=True, slots=True)
class LockEntry:
    """One row of ``docs/sources_lock.tsv``."""

    name: str
    path: str
    size: int
    sha256: str
    recorded: str  # ISO date of the first verification
    etag: str = ""  # server ETag when first verified (rolling sources)
    last_modified: str = ""  # server Last-Modified, likewise


@dataclass(frozen=True, slots=True)
class RemoteMeta:
    """What a HEAD request said about a file (empty when not sent)."""

    etag: str = ""
    last_modified: str = ""
    size: int | None = None


class HttpResponse(Protocol):
    """The part of :class:`requests.Response` that :func:`download` uses."""

    @property
    def status_code(self) -> int: ...

    @property
    def headers(self) -> Mapping[str, str]: ...

    def iter_content(self, chunk_size: int) -> Iterator[bytes]: ...

    def close(self) -> None: ...


class HttpSession(Protocol):
    """The part of :class:`requests.Session` that :func:`download` uses."""

    def get(
        self,
        url: str,
        *,
        headers: Mapping[str, str],
        stream: bool,
        timeout: tuple[float, float],
    ) -> HttpResponse: ...


class HeadResponse(Protocol):
    """The part of :class:`requests.Response` that :func:`head` uses."""

    @property
    def status_code(self) -> int: ...

    @property
    def headers(self) -> Mapping[str, str]: ...

    def close(self) -> None: ...


class HeadSession(Protocol):
    """The part of :class:`requests.Session` that :func:`head` uses."""

    def head(
        self,
        url: str,
        *,
        headers: Mapping[str, str],
        timeout: tuple[float, float],
        allow_redirects: bool,
    ) -> HeadResponse: ...


Log: TypeAlias = Callable[[str], None]


def get_source(name: str, registry: Mapping[str, Source] = SOURCES) -> Source:
    """Look up a registered source, with the known names on failure."""
    try:
        return registry[name]
    except KeyError:
        known = ", ".join(sorted(registry))
        raise KeyError(f"unknown source {name!r}; known: {known}") from None


def sources_for(
    code: str, registry: Mapping[str, Source] = SOURCES
) -> tuple[Source, ...]:
    """Sources used by the plugin with language code *code*, in order."""
    return tuple(s for s in registry.values() if code in s.languages)


def local_path(source: Source, root: Path = PROJECT_ROOT) -> Path:
    """Where *source* lives on disk below *root*."""
    return root.joinpath(*PurePosixPath(source.path).parts)


def user_agent() -> str:
    """Descriptive User-Agent naming the project and a contact."""
    contact = os.environ.get(CONTACT_ENV, "").strip()
    extra = f"; {contact}" if contact else ""
    return (
        "lexpipe/1.0 (romanian-palatalization multilingual research "
        f"pipeline; +{CONTACT_URL}{extra})"
    )


def _progress(written: int, expected: int | None) -> str:
    if not expected:
        return f"  {written:,} bytes"
    share = written / expected
    return f"  {written:,} / {expected:,} bytes ({share:.0%})"


def download(
    url: str,
    dest: Path,
    *,
    expected_size: int | None = None,
    session: HttpSession | None = None,
    log: Log = print,
    chunk_size: int = _CHUNK,
) -> Path:
    """Download *url* to *dest* through ``dest.part``, resuming if possible.

    An existing ``.part`` file is continued with an HTTP Range
    request; a server that ignores the range (200) restarts it. The
    file is renamed to *dest* only when complete, and, with
    *expected_size*, only when the byte count matches exactly.
    Raises :class:`DownloadError` otherwise, keeping the ``.part``
    file so the next call can resume.
    """
    part = dest.with_name(dest.name + ".part")
    dest.parent.mkdir(parents=True, exist_ok=True)
    offset = part.stat().st_size if part.exists() else 0
    if expected_size is not None and offset > expected_size:
        part.unlink()
        offset = 0
    if offset and offset == expected_size:
        # A previous run got every byte but stopped before the rename.
        os.replace(part, dest)
        return dest
    headers = {"User-Agent": user_agent(), "Accept-Encoding": "identity"}
    if offset:
        headers["Range"] = f"bytes={offset}-"
    with contextlib.ExitStack() as stack:
        http = (
            session
            if session is not None
            else stack.enter_context(requests.Session())
        )
        _transfer(
            http,
            url,
            part,
            headers=headers,
            offset=offset,
            expected_size=expected_size,
            log=log,
            chunk_size=chunk_size,
        )
    size = part.stat().st_size
    if expected_size is not None and size != expected_size:
        raise DownloadError(
            f"{url}: got {size:,} bytes, expected {expected_size:,}; "
            f"the partial file {part} is kept for resuming"
        )
    os.replace(part, dest)
    return dest


def _transfer(
    http: HttpSession,
    url: str,
    part: Path,
    *,
    headers: Mapping[str, str],
    offset: int,
    expected_size: int | None,
    log: Log,
    chunk_size: int,
) -> None:
    """One GET: append to the .part file (206) or rewrite it (200)."""
    response = http.get(
        url, headers=headers, stream=True, timeout=HTTP_TIMEOUT
    )
    try:
        status = response.status_code
        if status == 206 and offset:
            content_range = response.headers.get("Content-Range", "")
            if not content_range.startswith(f"bytes {offset}-"):
                raise DownloadError(
                    f"{url}: asked for bytes {offset}-, got Content-Range "
                    f"{content_range!r}"
                )
            mode, written = "ab", offset
            log(f"  resuming at byte {offset:,}")
        elif status == 200:
            mode, written = "wb", 0
        elif status == 416:
            raise DownloadError(
                f"GET {url}: HTTP 416, the server has no byte {offset:,}; "
                f"the upstream file changed, so delete {part} and retry"
            )
        else:
            raise DownloadError(f"GET {url}: HTTP {status}")
        last_report = time.monotonic()
        with part.open(mode) as fh:
            for chunk in response.iter_content(chunk_size):
                if not chunk:
                    continue
                fh.write(chunk)
                written += len(chunk)
                if expected_size is not None and written > expected_size:
                    raise DownloadError(
                        f"{url}: more than the expected {expected_size:,} "
                        "bytes; the upstream file has changed"
                    )
                now = time.monotonic()
                if now - last_report >= PROGRESS_INTERVAL_S:
                    log(_progress(written, expected_size))
                    last_report = now
    finally:
        response.close()


def fetch(
    name: str,
    *,
    root: Path = PROJECT_ROOT,
    registry: Mapping[str, Source] = SOURCES,
    session: HttpSession | None = None,
    log: Log = print,
) -> Path:
    """Make sure source *name* is on disk; download it if absent.

    A file that already has the registered size is left alone. A file
    of any other size raises :class:`DownloadError` without being
    touched: another download may still be writing it. A derived
    source is never downloaded: if it is absent,
    :class:`DerivedSourceMissing` says how to make it
    (:func:`derive_instructions`).
    """
    source = get_source(name, registry)
    dest = local_path(source, root)
    if dest.exists():
        size = dest.stat().st_size
        if size == source.size:
            log(f"{name}: present ({size:,} bytes)")
            return dest
        if source.derived:
            raise DownloadError(
                f"{name}: {dest} has {size:,} bytes, expected "
                f"{source.size:,}. It is derived from "
                f"{source.derived_from}, not downloaded; it was left "
                f"untouched. Delete it and make it again: {source.recipe}"
            )
        raise DownloadError(
            f"{name}: {dest} has {size:,} bytes, expected "
            f"{source.size:,}. Another download may still be writing "
            "it, or it is truncated or a different release. It was "
            "left untouched; delete it to fetch again."
        )
    if source.derived:
        raise DerivedSourceMissing(derive_instructions(source))
    log(f"{name}: downloading {source.url}")
    return download(
        source.url,
        dest,
        expected_size=source.size,
        session=session,
        log=log,
    )


def derive_instructions(source: Source) -> str:
    """What to run to make a derived *source* (fetch never downloads it)."""
    if not source.derived:
        raise ValueError(f"{source.name} is downloaded, not derived")
    return (
        f"{source.path} is derived from {source.derived_from} and never "
        "downloaded. Make it from the repository root with\n"
        f"    {source.recipe}\n"
        "then run verify-sources (docs/runbooks.md, section 2)."
    )


def head(
    url: str,
    *,
    session: HeadSession | None = None,
) -> RemoteMeta:
    """One HEAD request: the server's ETag, Last-Modified and size.

    Redirects are followed. A status other than 200 raises
    :class:`DownloadError`; network errors propagate (they are
    :class:`OSError` subclasses in requests).
    """
    headers = {"User-Agent": user_agent(), "Accept-Encoding": "identity"}
    with contextlib.ExitStack() as stack:
        http: HeadSession = (
            session
            if session is not None
            else stack.enter_context(requests.Session())
        )
        response = http.head(
            url, headers=headers, timeout=HTTP_TIMEOUT, allow_redirects=True
        )
        try:
            if response.status_code != 200:
                raise DownloadError(f"HEAD {url}: HTTP {response.status_code}")
            fields = {k.lower(): v for k, v in response.headers.items()}
        finally:
            response.close()
    length = fields.get("content-length", "").strip()
    return RemoteMeta(
        etag=fields.get("etag", "").strip(),
        last_modified=fields.get("last-modified", "").strip(),
        size=int(length) if length.isascii() and length.isdigit() else None,
    )


def sha256_file(path: Path, chunk_size: int = _CHUNK) -> str:
    """Hex SHA-256 of a file, read in chunks."""
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def read_lock(lock_path: Path = LOCK_PATH) -> dict[str, LockEntry]:
    """Read the sha256 lock file; an absent file is an empty lock.

    Locks with the five columns written before the upstream columns
    existed are accepted (their ETag and Last-Modified are empty).
    """
    if not lock_path.exists():
        return {}
    with lock_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, **_LOCK_DIALECT)
        header = tuple(reader.fieldnames or ())
        if header not in (LOCK_FIELDS, _LEGACY_LOCK_FIELDS):
            raise ValueError(
                f"{lock_path}: expected columns {LOCK_FIELDS}, "
                f"got {reader.fieldnames}"
            )
        entries: dict[str, LockEntry] = {}
        for row in reader:
            entries[row["name"]] = LockEntry(
                name=row["name"],
                path=row["path"],
                size=int(row["bytes"]),
                sha256=row["sha256"],
                recorded=row["recorded"],
                etag=row.get("etag") or "",
                last_modified=row.get("last_modified") or "",
            )
    return entries


def write_lock(
    entries: Mapping[str, LockEntry], lock_path: Path = LOCK_PATH
) -> None:
    """Write the lock file sorted by name (stable diffs), atomically."""
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = lock_path.with_name(lock_path.name + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n", **_LOCK_DIALECT)
        writer.writerow(LOCK_FIELDS)
        for name in sorted(entries):
            e = entries[name]
            writer.writerow(
                [
                    e.name,
                    e.path,
                    e.size,
                    e.sha256,
                    e.recorded,
                    e.etag,
                    e.last_modified,
                ]
            )
    os.replace(tmp, lock_path)


def _upstream_fields(
    remote: RemoteMeta | None, size: int
) -> tuple[str, str, str]:
    """(etag, last_modified, note) worth recording for a local file.

    Nothing is recorded when the server reports another size than the
    local file: it already serves another release.
    """
    if remote is None:
        return "", "", ""
    if remote.size is not None and remote.size != size:
        return (
            "",
            "",
            f"; upstream now serves {remote.size:,} bytes (another "
            "release), so its ETag was not recorded",
        )
    return remote.etag, remote.last_modified, ""


def _upstream_check(
    pinned: LockEntry, remote: RemoteMeta, size: int
) -> tuple[LockEntry, str, bool]:
    """Compare a HEAD answer with the lock row of a pinned file.

    Returns the (possibly completed) lock row, a note for the result
    detail and whether the upstream file changed.
    """
    if not (pinned.etag or pinned.last_modified):
        etag, last_modified, note = _upstream_fields(remote, size)
        if not (etag or last_modified):
            return pinned, note, bool(note)
        completed = dataclasses.replace(
            pinned, etag=etag, last_modified=last_modified
        )
        return completed, "; upstream ETag/Last-Modified recorded", False
    # Only the headers the server sent are compared.
    if not (remote.etag or remote.last_modified or remote.size is not None):
        return pinned, "; upstream sent no ETag, Last-Modified or size", False
    changed = (
        bool(remote.etag and remote.etag != pinned.etag)
        or bool(
            remote.last_modified
            and remote.last_modified != pinned.last_modified
        )
        or (remote.size is not None and remote.size != size)
    )
    if not changed:
        return pinned, "; upstream unchanged", False
    note = (
        f"; upstream changed since {pinned.recorded}: ETag "
        f"{pinned.etag or '?'} -> {remote.etag or '?'}, Last-Modified "
        f"{pinned.last_modified or '?'} -> {remote.last_modified or '?'}"
        " (the local file still matches the lock; a new download "
        "would not)"
    )
    return pinned, note, True


def verify(
    name: str,
    *,
    root: Path = PROJECT_ROOT,
    registry: Mapping[str, Source] = SOURCES,
    lock_path: Path = LOCK_PATH,
    today: date | None = None,
    remote: RemoteMeta | None = None,
) -> VerifyResult:
    """Check size and SHA-256 of source *name* against the lock file.

    The hash of a file with the registered size is recorded in the
    lock on first verification (status RECORDED) and must match on
    every later one (PINNED, else SHA256_MISMATCH). A source pinned
    in the registry (``Source.sha256``) must also match that hash,
    even on first verification. A file with the wrong size is never
    hashed or recorded; a lock row for another path or size is
    LOCK_MISMATCH (the registry entry changed: review, then delete
    the row to re-pin).

    *remote* is the answer to a HEAD request for a rolling source
    (:func:`head`); see the module docstring. Without it the upstream
    columns are left as they are.

    A derived source is checked like a download, and also against its
    parent's lock row: if the parent is locked at another SHA-256 than
    ``derived_from_sha256``, the file was made from another release
    and is STALE (nothing is recorded). A parent without a lock row
    yet is noted, not refused.
    """
    source = get_source(name, registry)
    path = local_path(source, root)
    if not path.is_file():
        detail = "absent"
        if source.derived:
            detail += "; " + derive_instructions(source).replace("\n", " ")
        return VerifyResult(name, path, VerifyStatus.MISSING, detail=detail)
    size = path.stat().st_size
    if size != source.size:
        return VerifyResult(
            name,
            path,
            VerifyStatus.SIZE_MISMATCH,
            size=size,
            detail=f"{size:,} bytes, expected {source.size:,}",
        )
    digest = sha256_file(path)
    if source.sha256 and digest != source.sha256:
        return VerifyResult(
            name,
            path,
            VerifyStatus.SHA256_MISMATCH,
            size=size,
            sha256=digest,
            detail=f"registry pins {source.sha256}",
        )
    lock = read_lock(lock_path)
    stale, parent_note = _parent_check(source, lock)
    if stale:
        return VerifyResult(
            name,
            path,
            VerifyStatus.STALE,
            size=size,
            sha256=digest,
            detail=parent_note,
        )
    pinned = lock.get(name)
    if pinned is None:
        etag, last_modified, note = _upstream_fields(remote, size)
        lock[name] = LockEntry(
            name=name,
            path=source.path,
            size=size,
            sha256=digest,
            recorded=(today or date.today()).isoformat(),
            etag=etag,
            last_modified=last_modified,
        )
        write_lock(lock, lock_path)
        return VerifyResult(
            name,
            path,
            VerifyStatus.RECORDED,
            size=size,
            sha256=digest,
            detail=f"sha256 recorded in {lock_path.name}{note}{parent_note}",
            upstream_changed=bool(note),
        )
    if (pinned.path, pinned.size) != (source.path, source.size):
        return VerifyResult(
            name,
            path,
            VerifyStatus.LOCK_MISMATCH,
            size=size,
            sha256=digest,
            detail=(
                f"lock row is for {pinned.path} ({pinned.size:,} bytes); "
                "the registry entry changed"
            ),
        )
    if pinned.sha256 != digest:
        return VerifyResult(
            name,
            path,
            VerifyStatus.SHA256_MISMATCH,
            size=size,
            sha256=digest,
            detail=f"lock file has {pinned.sha256}",
        )
    note, changed = "", False
    if remote is not None:
        row, note, changed = _upstream_check(pinned, remote, size)
        if row != pinned:
            lock[name] = row
            write_lock(lock, lock_path)
    return VerifyResult(
        name,
        path,
        VerifyStatus.PINNED,
        size=size,
        sha256=digest,
        detail=f"matches lock file{note}{parent_note}",
        upstream_changed=changed,
    )


def _parent_check(
    source: Source, lock: Mapping[str, LockEntry]
) -> tuple[bool, str]:
    """(stale, note) of a derived source against its parent's lock row;
    ``(False, "")`` for a download."""
    if not source.derived:
        return False, ""
    parent = lock.get(source.derived_from)
    if parent is None:
        return False, (
            f"; parent {source.derived_from} has no lock row yet, so the "
            "release it was made from is unchecked"
        )
    if parent.sha256 != source.derived_from_sha256:
        return True, (
            f"made from {source.derived_from} "
            f"{source.derived_from_sha256[:12]}…, but the lock pins "
            f"{source.derived_from} at {parent.sha256[:12]}…: make it "
            "again from the locked release and re-pin it "
            "(docs/runbooks.md, section 2)"
        )
    return False, (
        f"; made from the locked {source.derived_from} "
        f"({parent.sha256[:12]}…)"
    )
