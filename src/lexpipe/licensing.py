"""Which sources' data may be redistributed in this repository (R9).

Every row (or cell) the Italian and Czech pipelines produce records
the ``source`` it comes from, and may record a ``licence_class``. This
module maps each source code to its licence and decides whether data
derived from it may go into a **tracked** file:

    source          licence                       class       tracked
    wiktionary_en   CC BY-SA 4.0 + GFDL           by-sa       yes
    wiktionary_it   CC BY-SA 4.0 + GFDL           by-sa       yes
    wiktionary_cs   CC BY-SA 4.0 + GFDL           by-sa       yes
    wiktionary_ro   CC BY-SA 4.0 + GFDL           by-sa       yes
    morphit         CC BY-SA 2.0 / LGPL           by-sa       yes
    cnk             CC BY 4.0                     by          yes
    morfflex        CC BY-NC-SA 4.0               by-nc-sa    no
    derinet         CC BY-NC-SA 4.0               by-nc-sa    no
    czetyl          CC BY-NC-SA 4.0               by-nc-sa    no
    leipzig         terms unverified              unverified  no
    dex_qc          terms unverified (dexonline)  unverified  no
    treccani        coded value + citation        citation    yes
    ijp             coded value + citation        citation    yes

Notes:

  - CC BY-NC-SA data (MorfFlex, DeriNet, czetyl) and anything copied
    from it, QC worklists included, go to ``data/<language>/nc/``,
    which is gitignored. BY-SA and BY-NC-SA share-alike terms are
    incompatible, so a merged file is not redistributable either.
  - Leipzig-derived numbers stay untracked until a person has read
    the Leipzig download terms; then flip ``leipzig`` here.
  - ``wiktionary_ro``: Romanian Wiktionary pages, which supply most
    noun rows of the harvested Romanian lexicon
    (``scripts/romanian_harvester.py``); same licence as the other
    editions.
  - ``dex_qc``: the DEX quality check of the Romanian lexicon
    (dexonline.ro). Not redistributable until the authors have read
    dexonline's terms; then flip it here and in ``DATA_LICENSE.md``.
  - Treccani and IJP are consulted by hand. What is stored is a coded
    judgement and a citation, never copied entry text, which is why
    those rows may be tracked.

Tracked locations
-----------------
:func:`is_tracked_location` accepts only the places the design tracks
for the language plugins: ``analysis/<language>/``,
``analysis/comparison/``, ``docs/<language>/`` and the raw lexicon
``data/<language>/<prefix>lexicon_raw.csv``. A comparison-only
language (:attr:`lexpipe.spec.LanguageInfo.comparison_only`) has none
of its own: it tracks only ``analysis/comparison/``. ``curated/`` is
tracked too but is read-only to the pipeline, and files directly in
``analysis/`` or ``data/`` belong to the Romanian pipeline (R1), so
neither is accepted. :func:`lexpipe.csvio.write_tracked` enforces both
checks before writing a byte.

Aggregate counts from non-commercial sources (DRAFT policy)
-----------------------------------------------------------
Tables (CSV) never carry NC rows. Text reports may carry **aggregate
counts** (numbers of roots or pairs, odds ratios, test statistics)
computed in part from CC BY-NC-SA sources, under the draft policy in
``DATA_LICENSE.md`` that the authors still have to decide on.
:meth:`lexpipe.context.RunContext.write_tracked_text` is the gate:
with ``nc_aggregates=True`` it accepts every source whose licence is
known (:data:`AGGREGATE_CLASSES`: not ``unverified``, not unknown;
:func:`check_aggregate_sources`) and prepends :func:`aggregate_header`,
which names the sources and their licences and states that the report
holds aggregate counts only. Without the flag it accepts
redistributable sources only (:func:`check_redistributable`).

Source cells may hold several sources joined with ``|`` (a row merged
from two sources); every one of them must be redistributable. Unknown
source codes and empty cells are refused: a gate that cannot tell
must say no. :data:`DOWNLOADS` gives the source code of rows derived
from each file pinned in :mod:`lexpipe.sources`.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from types import MappingProxyType
from typing import Final

from lexpipe import PROJECT_ROOT
from lexpipe.spec import LANGUAGES, LanguageInfo

SOURCE_SEPARATOR: Final[str] = "|"


class LicenceClass(StrEnum):
    """Coarse licence classes, the values of a ``licence_class`` column."""

    BY_SA = "by-sa"  # CC BY-SA (alone or dual-licensed with GFDL/LGPL)
    BY = "by"  # CC BY
    BY_NC_SA = "by-nc-sa"  # non-commercial share-alike: never tracked
    UNVERIFIED = "unverified"  # terms not yet read by a person
    CITATION = "citation"  # coded judgement plus citation, no copied text


# Classes whose data may not appear in tracked files.
NON_REDISTRIBUTABLE_CLASSES: Final[frozenset[LicenceClass]] = frozenset(
    {LicenceClass.BY_NC_SA, LicenceClass.UNVERIFIED}
)
# Classes whose aggregate counts tracked text reports may carry under
# the draft NC policy (module docstring): every known licence, NC
# included; never terms nobody has read.
AGGREGATE_CLASSES: Final[frozenset[LicenceClass]] = frozenset(
    {
        LicenceClass.BY_SA,
        LicenceClass.BY,
        LicenceClass.CITATION,
        LicenceClass.BY_NC_SA,
    }
)
# First line of the header aggregate_header() writes.
AGGREGATE_HEADER_TITLE: Final[str] = "Licence and attribution"


class LicenceError(ValueError):
    """Data that may not be redistributed was headed for a tracked file."""


@dataclass(frozen=True, slots=True)
class SourceLicence:
    """Licence of one source and whether derived data may be tracked."""

    source: str
    licence: str
    licence_class: LicenceClass
    note: str = ""

    @property
    def redistributable(self) -> bool:
        return self.licence_class not in NON_REDISTRIBUTABLE_CLASSES


_WIKTIONARY: Final[str] = "CC BY-SA 4.0 + GFDL"
_NC: Final[str] = "CC BY-NC-SA 4.0"
_CODED: Final[str] = "coded value + citation (no entry text copied)"

_REGISTRY: Final[tuple[SourceLicence, ...]] = (
    SourceLicence(
        "wiktionary_en",
        _WIKTIONARY,
        LicenceClass.BY_SA,
        "English Wiktionary via wiktextract (kaikki.org)",
    ),
    SourceLicence(
        "wiktionary_it",
        _WIKTIONARY,
        LicenceClass.BY_SA,
        "Italian Wiktionary (Wikizionario) via wiktextract",
    ),
    SourceLicence(
        "wiktionary_cs",
        _WIKTIONARY,
        LicenceClass.BY_SA,
        "Czech Wiktionary (Wikislovník) via wiktextract",
    ),
    SourceLicence(
        "wiktionary_ro",
        _WIKTIONARY,
        LicenceClass.BY_SA,
        "Romanian Wiktionary pages harvested by the Romanian pipeline "
        "(scripts/romanian_harvester.py)",
    ),
    SourceLicence(
        "morphit",
        "CC BY-SA 2.0 / LGPL",
        LicenceClass.BY_SA,
        "Morph-it! 0.4.8 (Zanchetta & Baroni 2005), dual-licensed",
    ),
    SourceLicence(
        "cnk",
        "CC BY 4.0",
        LicenceClass.BY,
        "Czech National Corpus SYN2015 frequency lists",
    ),
    SourceLicence("morfflex", _NC, LicenceClass.BY_NC_SA, "MorfFlex CZ 2.1"),
    SourceLicence("derinet", _NC, LicenceClass.BY_NC_SA, "DeriNet 2.3"),
    SourceLicence(
        "czetyl",
        _NC,
        LicenceClass.BY_NC_SA,
        "Czech Etymological Lexicon 1.0",
    ),
    SourceLicence(
        "leipzig",
        "terms unverified",
        LicenceClass.UNVERIFIED,
        "Leipzig Corpora Collection; untracked until a person checks "
        "the download terms",
    ),
    SourceLicence(
        "dex_qc",
        "terms unverified (dexonline)",
        LicenceClass.UNVERIFIED,
        "DEX quality check of the Romanian lexicon (dexonline.ro); "
        "untracked until the authors check dexonline's terms",
    ),
    SourceLicence(
        "treccani",
        _CODED,
        LicenceClass.CITATION,
        "Treccani / De Mauro, consulted by hand (Treccani opts out of "
        "text and data mining)",
    ),
    SourceLicence(
        "ijp",
        _CODED,
        LicenceClass.CITATION,
        "Internetová jazyková příručka, consulted by hand",
    ),
)

LICENCES: Final[Mapping[str, SourceLicence]] = MappingProxyType(
    {entry.source: entry for entry in _REGISTRY}
)

# The source code of rows derived from each pinned download
# (lexpipe.sources registry names). tests/lexpipe/test_licensing.py
# checks that every registered download is listed.
DOWNLOADS: Final[Mapping[str, str]] = MappingProxyType(
    {
        "kaikki_raw": "wiktionary_en",
        "kaikki_raw_it": "wiktionary_en",  # derived: split of kaikki_raw
        "kaikki_raw_cs": "wiktionary_en",  # derived: split of kaikki_raw
        "kaikki_it_extract": "wiktionary_it",
        "kaikki_cs_extract": "wiktionary_cs",
        "morphit": "morphit",
        "morfflex": "morfflex",
        "derinet": "derinet",
        "czetyl": "czetyl",
        "cnk_syn2015_lemma": "cnk",
        "cnk_syn2015_word": "cnk",
        "leipzig_ita_news_2024_1M": "leipzig",
        "leipzig_ces_news_2024_1M": "leipzig",
    }
)


def licence_of(
    source: str, registry: Mapping[str, SourceLicence] = LICENCES
) -> SourceLicence:
    """Registry entry of *source*; :class:`LicenceError` if unknown."""
    try:
        return registry[source]
    except KeyError:
        known = ", ".join(sorted(registry))
        raise LicenceError(
            f"unknown source {source!r}; known: {known}"
        ) from None


def split_sources(value: str | Sequence[str]) -> list[str]:
    """Source codes of one cell: a ``|``-joined string or a sequence."""
    if isinstance(value, str):
        items = value.split(SOURCE_SEPARATOR) if value else []
    else:
        items = [str(item) for item in value]
    return [item.strip() for item in items]


def check_redistributable(
    sources: Iterable[str], registry: Mapping[str, SourceLicence] = LICENCES
) -> None:
    """Raise :class:`LicenceError` unless every source may be tracked.

    Empty and unknown source codes are refused too.
    """
    names = list(sources)
    if not names or any(not name for name in names):
        raise LicenceError(f"empty source in {names!r}")
    refused = [
        f"{name} ({entry.licence})"
        for name in names
        if not (entry := licence_of(name, registry)).redistributable
    ]
    if refused:
        raise LicenceError(
            "not redistributable in the repository: "
            + ", ".join(refused)
            + "; write NC-bearing outputs under data/<language>/nc/"
        )


def check_licence_class(value: str) -> None:
    """Raise :class:`LicenceError` unless every ``|``-joined class in a
    ``licence_class`` cell is known and redistributable."""
    for item in split_sources(value):
        try:
            licence_class = LicenceClass(item)
        except ValueError:
            raise LicenceError(f"unknown licence class {item!r}") from None
        if licence_class in NON_REDISTRIBUTABLE_CLASSES:
            raise LicenceError(
                f"licence class {item!r} may not appear in tracked files"
            )


def check_aggregate_sources(
    sources: Iterable[str], registry: Mapping[str, SourceLicence] = LICENCES
) -> list[SourceLicence]:
    """Entries of *sources* if a tracked report may carry their aggregate
    counts under the draft NC policy (module docstring).

    Every source must be known and of a class in
    :data:`AGGREGATE_CLASSES`: CC BY-NC-SA passes, ``unverified``
    terms (Leipzig, DEX QC) and unknown or empty codes do not
    (:class:`LicenceError`). Returns the entries sorted by source code,
    each once.
    """
    names = list(sources)
    if not names or any(not name for name in names):
        raise LicenceError(f"empty source in {names!r}")
    entries = {name: licence_of(name, registry) for name in names}
    refused = [
        f"{name} ({entry.licence})"
        for name, entry in sorted(entries.items())
        if entry.licence_class not in AGGREGATE_CLASSES
    ]
    if refused:
        raise LicenceError(
            "not even aggregate counts may be tracked from sources whose "
            "terms nobody has checked: " + ", ".join(refused)
        )
    return [entries[name] for name in sorted(entries)]


def aggregate_header(
    sources: Iterable[str], registry: Mapping[str, SourceLicence] = LICENCES
) -> str:
    """The licence and attribution block of a tracked aggregate report.

    Names every source (sorted, once) with its licence, states that
    the report holds aggregate counts only, and, when a CC BY-NC-SA
    source is among them, that those counts are tracked under the
    draft policy of ``DATA_LICENSE.md``. Ends with a blank line, so
    the report text can follow directly. Raises like
    :func:`check_aggregate_sources`.
    """
    entries = check_aggregate_sources(sources, registry)
    width = max(len(entry.source) for entry in entries)
    lines = [
        AGGREGATE_HEADER_TITLE,
        "-" * len(AGGREGATE_HEADER_TITLE),
        "Aggregate counts only: this report holds counts and statistics",
        "computed from the sources below, never their rows, word forms or",
        "text. Sources and licences:",
    ]
    lines += [
        f"  {entry.source:<{width}}  {entry.licence}"
        + (f"; {entry.note}" if entry.note else "")
        for entry in entries
    ]
    if any(e.licence_class is LicenceClass.BY_NC_SA for e in entries):
        lines += [
            "Counts that rest on CC BY-NC-SA sources are tracked under the",
            "draft policy of DATA_LICENSE.md ('Aggregate counts from",
            "non-commercial sources'), pending the authors' decision; they",
            "may not be used commercially.",
        ]
    return "\n".join(lines) + "\n\n"


def with_language(
    info: LanguageInfo, languages: Mapping[str, LanguageInfo] = LANGUAGES
) -> dict[str, LanguageInfo]:
    """*languages* plus *info* when its code is not registered.

    The registry stays authoritative for registered codes (a hand-made
    ``ro`` row without ``comparison_only`` does not make
    ``analysis/romanian_compare/`` tracked); a test language brings its
    own row.
    """
    registry = dict(languages)
    registry.setdefault(info.code, info)
    return registry


def _relative_parts(path: Path, root: Path) -> tuple[str, ...] | None:
    try:
        return Path(path).resolve().relative_to(root.resolve()).parts
    except ValueError:
        return None


def is_tracked_location(
    path: Path,
    root: Path = PROJECT_ROOT,
    languages: Mapping[str, LanguageInfo] = LANGUAGES,
) -> bool:
    """True if *path* is a place where the pipeline writes tracked files.

    See the module docstring; the language names and file prefixes
    come from the registry. Comparison-only languages add no location
    of their own.
    """
    parts = _relative_parts(path, root)
    if parts is None or len(parts) < 3:
        return False
    full = [info for info in languages.values() if not info.comparison_only]
    names = {info.name for info in full}
    top, directory = parts[0], parts[1]
    if top == "analysis":
        return directory in names or directory == "comparison"
    if top == "docs":
        return directory in names
    if top == "data" and len(parts) == 3:
        return any(
            info.name == directory
            and parts[2] == info.file_name("lexicon_raw")
            for info in full
        )
    return False


def check_tracked_location(
    path: Path,
    root: Path = PROJECT_ROOT,
    languages: Mapping[str, LanguageInfo] = LANGUAGES,
) -> None:
    """Raise :class:`LicenceError` unless :func:`is_tracked_location`."""
    if not is_tracked_location(path, root, languages):
        raise LicenceError(
            f"{path} is not a tracked location (analysis/<language>/, "
            "analysis/comparison/, docs/<language>/ or "
            "data/<language>/<prefix>lexicon_raw.csv); NC-bearing "
            "outputs go under data/<language>/nc/"
        )
