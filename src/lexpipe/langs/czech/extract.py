"""Stage ``extract``: compact derivatives of the two Wiktionary editions (R7).

One streaming pass over each kaikki file writes what later stages read
instead of the 224 MB English subset and the 40 MB Czech extract:

``czech_en_entries.jsonl.gz``
    Every English Wiktionary **lexeme** entry of the Czech subset (any
    part of speech; :func:`is_lexeme_entry`) projected to the fields
    the plugin uses (:func:`project_en`). Head and inflection template
    arguments are kept verbatim, so :func:`lexpipe.ids.lexeme_id` gives
    the same fingerprint on the projection as on the raw record. Only
    "Czech …" and topical "cs:…" categories are kept.

    A lexeme entry is a lemma entry (:func:`lexpipe.kaikki.is_lemma_entry`)
    or an entry with a **derivational** form-of sense
    (:func:`derivational_relation`): wiktextract tags "diminutive of
    kluk" (klouček), "female equivalent of Čech" (Češka), "verbal noun
    of kvasit" (kvašení), "augmentative of …" as form-of, but these
    words are lexemes with their own paradigms, and the R26 suffix
    targets (-ek, -ka, -ička, -ík) are mostly among them. Inflectional
    form-of senses (cases, numbers), spelling variants (alternative,
    obsolete, misspelling, abbreviation) and the adjective degrees
    (comparative, superlative: forms of the adjective's own paradigm)
    stay non-lexemes; so do possessive adjectives of given names.
``czech_cs_entries.jsonl.gz``
    The Czech Wiktionary lemma entries with ``lang_code == "cs"``
    (:func:`project_cs`).
``czech_form_index.sqlite``
    Table ``forms``: which entries own a word form, for the R11
    ``form_ambiguous`` flag. One row per (form, owner): the lemma
    itself (kind ``lemma``), every form of its inflection tables
    (``paradigm``; head-line forms such as feminine equivalents are
    other lexemes and are left out), and the targets of form-of and
    alt-of senses (``form_of``/``alt_of``, owner = the target lemma).
    A derivational form-of sense (klouček 'diminutive of kluk') names
    another lexeme, not a form of it, and gives no ``form_of`` row.
    English owners carry their base lexeme_id. Indexed by form and by
    folded form.
``czech_extract_stats.tsv``
    Counts and the SHA-256 of both inputs.

Inputs are untrusted: lines are parsed with :mod:`json` and every
record must be in the raw wiktextract format
(:func:`lexpipe.kaikki.assert_raw_format`); an English record whose
top-level ``lang_code`` is not ``cs`` is an error (the subset was split
wrongly). Projected lines are canonical JSON, sorted by (word, pos,
etymology, line), so the derivatives do not depend on input order
(R5).
"""

from __future__ import annotations

import re
import unicodedata
from collections import Counter
from collections.abc import Callable, Iterator, Mapping
from pathlib import Path
from typing import Any, Final

from lexpipe import ids, kaikki
from lexpipe.context import RunContext, run_stage, sha256_file
from lexpipe.langs.czech import jsonio
from lexpipe.langs.czech.layout import CzechPaths, paths_for

STAGE: Final[str] = "extract"
# 2: derivational form-of entries (diminutives, female equivalents,
# verbal nouns, augmentatives) are lexemes (audit 2026-10-23).
VERSION: Final[str] = "2"

# Derivational relations a form-of gloss can name (module docstring),
# with the relation label recorded by the lexicon stage. The first
# gloss of the sense must start with one of these heads.
DERIVATIONAL_RELATIONS: Final[tuple[tuple[str, str], ...]] = (
    ("diminutive", "diminutive"),
    ("augmentative", "augmentative"),
    ("female equivalent", "female_equivalent"),
    ("feminine equivalent", "female_equivalent"),
    ("feminine", "female_equivalent"),
    ("female", "female_equivalent"),
    ("male equivalent", "male_equivalent"),
    ("masculine equivalent", "male_equivalent"),
    ("verbal noun", "verbal_noun"),
)
_DERIVATIONAL_GLOSS: Final[re.Pattern[str]] = re.compile(
    r"^(?:an?\s+)?("
    + "|".join(re.escape(head) for head, _ in DERIVATIONAL_RELATIONS)
    + r")\s+(?:form\s+)?of\b",
    re.IGNORECASE,
)
_RELATION_OF: Final[Mapping[str, str]] = dict(DERIVATIONAL_RELATIONS)

_SENSE_KEYS: Final[tuple[str, ...]] = (
    "glosses",
    "tags",
    "raw_tags",
    "form_of",
    "alt_of",
)
_FORM_KEYS: Final[tuple[str, ...]] = ("form", "tags", "raw_tags", "source")
_LINK_KEYS: Final[tuple[str, ...]] = ("word", "tags", "raw_tags")
INDEX_COLUMNS: Final[tuple[tuple[str, str], ...]] = (
    ("form", "TEXT"),
    ("form_folded", "TEXT"),
    ("kind", "TEXT"),
    ("edition", "TEXT"),
    ("pos", "TEXT"),
    ("owner_word", "TEXT"),
    ("owner_id", "TEXT"),
    ("target", "TEXT"),
)
EDITION_EN: Final[str] = "en"
EDITION_CS: Final[str] = "cs"


# ---------------------------------------------------------------------------
# Projection
# ---------------------------------------------------------------------------


def derivational_relation(sense: Mapping[str, Any]) -> str:
    """The derivational relation a form-of *sense* names (module
    docstring): ``diminutive``, ``augmentative``, ``female_equivalent``,
    ``male_equivalent`` or ``verbal_noun``; "" for any other sense.

    Only form-of senses (``form_of``/``alt_of`` or the form-of tag)
    qualify, and their first gloss must start with the relation
    ("diminutive of kluk", "female equivalent of Čech").
    """
    tags = sense.get("tags") or ()
    if not (
        "form_of" in sense
        or "alt_of" in sense
        or kaikki.FORM_OF_TAGS.intersection(tags)
    ):
        return ""
    glosses = [g for g in sense.get("glosses") or () if isinstance(g, str)]
    if not glosses:
        return ""
    found = _DERIVATIONAL_GLOSS.match(glosses[0].strip())
    if found is None:
        return ""
    return _RELATION_OF[found.group(1).lower()]


def is_lexeme_entry(entry: Mapping[str, Any]) -> bool:
    """A lemma entry (:func:`lexpipe.kaikki.is_lemma_entry`) or an entry
    with a derivational form-of sense (module docstring)."""
    if kaikki.is_lemma_entry(entry):
        return True
    return any(
        isinstance(sense, Mapping) and derivational_relation(sense)
        for sense in entry.get("senses") or ()
    )


def _subset(obj: Mapping[str, Any], keys: tuple[str, ...]) -> dict[str, Any]:
    return {key: obj[key] for key in keys if key in obj}


def _dicts(value: Any) -> list[Mapping[str, Any]]:
    return [item for item in value or () if isinstance(item, Mapping)]


def _keep_category(name: Any) -> bool:
    return isinstance(name, str) and (
        name.startswith("Czech ") or name.startswith("cs:")
    )


def _project_sense(sense: Mapping[str, Any]) -> dict[str, Any]:
    out = _subset(sense, _SENSE_KEYS)
    for key in ("form_of", "alt_of"):
        if key in out:
            out[key] = [_subset(t, ("word",)) for t in _dicts(out[key])]
    categories = [
        c for c in sense.get("categories") or () if _keep_category(c)
    ]
    if categories:
        out["categories"] = categories
    return out


def _templates(value: Any, keys: tuple[str, ...]) -> list[dict[str, Any]]:
    return [_subset(t, keys) for t in _dicts(value)]


def project_en(entry: Mapping[str, Any]) -> dict[str, Any]:
    """The fields of an English Wiktionary entry the plugin uses."""
    out: dict[str, Any] = _subset(
        entry, ("word", "pos", "lang_code", "etymology_number")
    )
    out["head_templates"] = _templates(
        entry.get("head_templates"), ("name", "args", "expansion")
    )
    out["inflection_templates"] = _templates(
        entry.get("inflection_templates"), ("name", "args")
    )
    out["forms"] = [_subset(f, _FORM_KEYS) for f in _dicts(entry.get("forms"))]
    out["senses"] = [_project_sense(s) for s in _dicts(entry.get("senses"))]
    out["categories"] = [
        c for c in entry.get("categories") or () if _keep_category(c)
    ]
    out["sounds"] = [
        _subset(s, ("ipa", "tags"))
        for s in _dicts(entry.get("sounds"))
        if isinstance(s.get("ipa"), str)
    ]
    out["etymology_templates"] = _templates(
        entry.get("etymology_templates"), ("name", "args")
    )
    if isinstance(entry.get("etymology_text"), str):
        out["etymology_text"] = entry["etymology_text"]
    for key in ("derived", "related"):
        links = [_subset(d, _LINK_KEYS) for d in _dicts(entry.get(key))]
        if links:
            out[key] = links
    return out


def project_cs(entry: Mapping[str, Any]) -> dict[str, Any]:
    """The fields of a Czech Wiktionary entry the plugin uses."""
    out: dict[str, Any] = _subset(
        entry, ("word", "pos", "lang_code", "tags", "raw_tags")
    )
    out["forms"] = [_subset(f, _FORM_KEYS) for f in _dicts(entry.get("forms"))]
    out["senses"] = [_project_sense(s) for s in _dicts(entry.get("senses"))]
    out["sounds"] = [
        _subset(s, ("ipa", "tags"))
        for s in _dicts(entry.get("sounds"))
        if isinstance(s.get("ipa"), str)
    ]
    return out


def sort_key(
    record: Mapping[str, Any], line: str
) -> tuple[str, str, str, str]:
    """Total order of projected records: word, pos, etymology, line."""
    return (
        str(record.get("word", "")),
        str(record.get("pos", "")),
        kaikki.etymology_label(record),
        line,
    )


# ---------------------------------------------------------------------------
# Form index
# ---------------------------------------------------------------------------


def _nfc(text: str) -> str:
    return unicodedata.normalize("NFC", text).strip()


def index_rows(
    entry: Mapping[str, Any], edition: str
) -> Iterator[tuple[str, ...]]:
    """Form-index rows of one entry (see the module docstring)."""
    word = _nfc(str(entry.get("word") or ""))
    pos = str(entry.get("pos") or "")
    if not word or not pos:
        return
    if is_lexeme_entry(entry):
        owner = ids.lexeme_id("cs", entry) if edition == EDITION_EN else ""
        yield _row(word, "lemma", edition, pos, word, owner, word)
        for form in kaikki.inflected_forms(entry):
            if edition == EDITION_EN and not form.source:
                continue  # head-line forms are other lexemes
            yield _row(
                _nfc(form.form), "paradigm", edition, pos, word, owner, word
            )
    for sense in _dicts(entry.get("senses")):
        if derivational_relation(sense):
            continue  # names another lexeme, not a form of it
        for kind in ("form_of", "alt_of"):
            for target in _dicts(sense.get(kind)):
                lemma = target.get("word")
                if isinstance(lemma, str) and lemma.strip():
                    yield _row(word, kind, edition, pos, word, "", _nfc(lemma))


def _row(
    form: str,
    kind: str,
    edition: str,
    pos: str,
    owner_word: str,
    owner_id: str,
    target: str,
) -> tuple[str, ...]:
    return (
        form,
        form.lower(),
        kind,
        edition,
        pos,
        owner_word,
        owner_id,
        target,
    )


# ---------------------------------------------------------------------------
# Stage
# ---------------------------------------------------------------------------


def ensure_en_subset(ctx: RunContext, paths: CzechPaths) -> Path:
    """The English subset to read, splitting the raw dump if needed."""
    found = paths.en_input()
    if found.is_file():
        return found
    if not paths.raw_dump.is_file():
        raise FileNotFoundError(
            f"neither {paths.en_subset} nor the raw dump {paths.raw_dump} "
            "exists; run the fetch stage"
        )
    ctx.check_writable(paths.split_subset)
    print(f"  splitting {paths.raw_dump.name} (lang_code cs)", flush=True)
    kaikki.split_raw_dump(
        paths.raw_dump, paths.split_subset.parent, [ctx.language.kaikki_code]
    )
    return paths.split_subset


def run(ctx: RunContext) -> int | None:
    """Stage entry point (see the module docstring)."""
    paths = paths_for(ctx)
    en_input = ensure_en_subset(ctx, paths)
    return run_stage(
        ctx,
        STAGE,
        inputs=[en_input, paths.cs_extract],
        params={"version": VERSION, "kaikki_code": ctx.language.kaikki_code},
        outputs=[
            paths.en_entries,
            paths.cs_entries,
            paths.form_index,
            paths.extract_stats,
        ],
        body=lambda: _extract(ctx, paths, en_input),
    )


class _Collector:
    """Streams both editions once: yields form-index rows and keeps the
    projected lemma lines and the counts on the side."""

    def __init__(self, code: str) -> None:
        self.code = code
        self.stats: Counter[tuple[str, str]] = Counter()
        self.en_lines: list[tuple[tuple[str, str, str, str], str]] = []
        self.cs_lines: list[tuple[tuple[str, str, str, str], str]] = []

    def rows(
        self, en_input: Path, cs_input: Path
    ) -> Iterator[tuple[str, ...]]:
        for entry in kaikki.iter_jsonl(en_input):
            if entry.get("lang_code") != self.code:
                raise ValueError(
                    f"{en_input}: entry {entry.get('word')!r} has lang_code "
                    f"{entry.get('lang_code')!r}, expected {self.code!r}"
                )
            kaikki.assert_raw_format(entry)
            yield from index_rows(entry, EDITION_EN)
            self._keep(entry, EDITION_EN, project_en, self.en_lines)
        for entry in kaikki.iter_jsonl(cs_input):
            if entry.get("lang_code") != self.code:
                continue
            kaikki.assert_raw_format(entry)
            yield from index_rows(entry, EDITION_CS)
            self._keep(entry, EDITION_CS, project_cs, self.cs_lines)

    def _keep(
        self,
        entry: Mapping[str, Any],
        edition: str,
        project: Callable[[Mapping[str, Any]], dict[str, Any]],
        lines: list[tuple[tuple[str, str, str, str], str]],
    ) -> None:
        pos = str(entry.get("pos") or "")
        if not is_lexeme_entry(entry):
            self.stats[(f"{edition}_formof", pos)] += 1
            return
        record = project(entry)
        line = jsonio.canonical_line(record)
        lines.append((sort_key(record, line), line))
        self.stats[(f"{edition}_lemma", pos)] += 1
        if not kaikki.is_lemma_entry(entry):
            self.stats[(f"{edition}_derivational_formof", pos)] += 1


def _extract(ctx: RunContext, paths: CzechPaths, en_input: Path) -> int:
    digests = (
        ("source_sha256:" + en_input.name, sha256_file(en_input)),
        (
            "source_sha256:" + paths.cs_extract.name,
            sha256_file(paths.cs_extract),
        ),
    )
    collector = _Collector(ctx.language.kaikki_code)
    stored = jsonio.build_sqlite(
        ctx,
        paths.form_index,
        table="forms",
        columns=INDEX_COLUMNS,
        rows=collector.rows(en_input, paths.cs_extract),
        order_by=[name for name, _ in INDEX_COLUMNS[1:]] + ["form"],
        indexes=(
            ("forms_folded", ("form_folded",)),
            ("forms_form", ("form",)),
        ),
        meta=digests,
        distinct=True,
    )
    collector.en_lines.sort()
    collector.cs_lines.sort()
    jsonio.write_jsonl_gz(
        ctx, paths.en_entries, (line for _, line in collector.en_lines)
    )
    jsonio.write_jsonl_gz(
        ctx, paths.cs_entries, (line for _, line in collector.cs_lines)
    )
    _write_stats(ctx, paths, collector.stats, digests, stored)
    print(
        f"  en lemma entries {len(collector.en_lines):,}; cs lemma entries "
        f"{len(collector.cs_lines):,}; form index rows {stored:,}",
        flush=True,
    )
    return 0


def _write_stats(
    ctx: RunContext,
    paths: CzechPaths,
    stats: Counter[tuple[str, str]],
    digests: tuple[tuple[str, str], ...],
    index_rows_stored: int,
) -> None:
    lines = ["section\tkey\tvalue"]
    lines.extend(f"source\t{key}\t{value}" for key, value in digests)
    lines.extend(
        f"{section}\t{pos or '-'}\t{count}"
        for (section, pos), count in sorted(stats.items())
    )
    lines.append(f"form_index\trows\t{index_rows_stored}")
    ctx.write_text(paths.extract_stats, "\n".join(lines) + "\n")
