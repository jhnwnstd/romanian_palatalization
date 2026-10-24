"""Stream and query wiktextract (kaikki.org) JSONL dumps.

The Italian and Czech lexicons are built from the kaikki.org
extraction of English Wiktionary rather than by scraping HTML: one
JSON object per line, one line per (page, language, etymology, part
of speech) entry. This module streams those files and answers the
questions every plugin asks of an entry.

Splitting the raw dump
----------------------
The raw dump (``raw-wiktextract-data.jsonl.gz``, ~3 GB compressed)
holds every language. :func:`split_raw_dump` keeps the entries whose
**top-level** ``lang_code`` is wanted. A substring test alone is
wrong: translation tables and descendant lists nest objects that
carry their own ``lang_code``, so an English entry with an Italian
translation contains ``"lang_code": "it"`` too. A regular expression
(``"lang_code"\\s*:\\s*"<code>"``, so any JSON spacing matches) is
only a fast prefilter; the decision is made on the parsed object.
Kept lines are written byte for byte, so the split files are exact
subsets of the dump.

Entry helpers
-------------
- :func:`entry_key` — ``(word, etymology_number, pos)``. Not unique:
  build identifiers with :mod:`lexpipe.ids`. :func:`etymology_label`
  — the etymology number as a string, sub-etymologies ("1.2")
  included.
- :func:`is_form_of_only` — the entry only points at other entries
  (inflected forms, alternative spellings); :func:`is_lemma_entry`
  — at least one sense is a lexeme sense (not form-of/alt-of, not a
  ``no-gloss`` stub).
- :func:`forms_with_tags` — inflected forms whose tags contain a
  required set and avoid an excluded set; table metadata and dash
  placeholders are never returned.
- :func:`inflected_forms` — every word form as a typed :class:`Form`,
  without table pseudo-forms (:data:`PSEUDO_FORM_TAGS`); wiktextract's
  tag errors (:data:`TAG_ERROR_TAGS`, e.g. fòcora) are kept and
  flagged in ``Form.tag_error`` rather than silently dropped.
- :func:`plural_forms` — forms tagged ``plural`` that are not
  evaluative (pantaloncìni is diminutive *and* plural).
- :func:`head_args` / :func:`args_text` — the head templates with their
  arguments as written, compactly (the tracked raw lexicons'
  ``head_args``); :func:`parse_args_text` reads them back.
- :func:`ipa_list`, :func:`etymology_templates`,
  :func:`head_template_args`, :func:`categories` (entry-level and
  sense-level; two thirds of the noun entries have their categories
  only on the senses; ``lang_name="Italian"`` keeps only that
  language's categories).
- :func:`assert_raw_format` — reject records from kaikki's
  post-processed per-language files (sense ``id``s, dict-valued
  categories): the pipeline reads the raw dump only.

Contracts
---------
- Dump files are untrusted data: they are parsed with :mod:`json`
  and never evaluated. A malformed line raises :class:`ValueError`
  naming the file and line number.
- Helpers never mutate the entry; returned dicts are shallow copies.
- Entries are read as wiktextract writes them: ``etymology_number``
  is a string ("1") in current dumps and an int in older ones;
  categories are strings, or dicts with a ``name`` in some
  post-processed extractions. Both shapes are accepted by the
  helpers; :func:`assert_raw_format` is the place to insist.
"""

from __future__ import annotations

import contextlib
import gzip
import io
import json
import os
import re
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, BinaryIO, Final

GZIP_MAGIC: Final[bytes] = b"\x1f\x8b"

# File name of one language's slice of the raw dump, inside the
# ``out_dir`` given to :func:`split_raw_dump`.
SPLIT_FILE_TEMPLATE: Final[str] = "en_wiktionary_{code}.jsonl"

# Sense tags marking an entry that only points at another entry.
FORM_OF_TAGS: Final[frozenset[str]] = frozenset({"form-of", "alt-of"})

# Pseudo-forms wiktextract emits for inflection-table bookkeeping
# ("animate" as a table header, the template name, the declension
# class). They are not word forms.
META_FORM_TAGS: Final[frozenset[str]] = frozenset(
    {"table-tags", "inflection-template", "class"}
)

# Every pseudo-form :func:`inflected_forms` drops: the table metadata
# above, transliterations, and the headword's canonical spelling
# (diminutìve, "il poi"), which is not an inflected form.
PSEUDO_FORM_TAGS: Final[frozenset[str]] = META_FORM_TAGS | frozenset(
    {"romanization", "canonical"}
)

# Tags wiktextract adds when it could not parse a form's labels.
TAG_ERROR_TAGS: Final[frozenset[str]] = frozenset(
    {"error-unknown-tag", "error-unrecognized-form"}
)

# Evaluative derivatives that wiktextract also tags "plural".
EVALUATIVE_TAGS: Final[frozenset[str]] = frozenset(
    {"diminutive", "augmentative", "pejorative", "endearing", "superlative"}
)

# Table cells that mean "no such form".
PLACEHOLDER_FORMS: Final[frozenset[str]] = frozenset({"", "-", "–", "—"})

# Language codes end up in file names and in a byte needle, so keep
# them to the Wiktionary code alphabet.
_LANG_CODE_RE: Final[re.Pattern[str]] = re.compile(r"[a-z][a-z0-9-]*")


def _open_binary(path: Path) -> io.BufferedIOBase:
    """Open *path* for binary reading, gunzipping if it is gzip data.

    Detection uses the gzip magic number, not the file name, so a
    renamed file still reads correctly.
    """
    with path.open("rb") as probe:
        magic = probe.read(2)
    if magic == GZIP_MAGIC:
        return gzip.open(path, "rb")
    return path.open("rb")


def _parse_line(line: bytes, source: Path, lineno: int) -> Any:
    try:
        return json.loads(line)
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"{source}:{lineno}: invalid JSON ({exc.msg})"
        ) from exc


def iter_jsonl(path: Path | str) -> Iterator[dict[str, Any]]:
    """Yield one dict per non-blank line of a JSONL file (plain or gzip).

    Raises :class:`ValueError` with ``path:line`` for a line that is
    not valid JSON or does not hold a JSON object.
    """
    source = Path(path)
    with _open_binary(source) as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.strip():
                continue
            obj = _parse_line(line, source, lineno)
            if not isinstance(obj, dict):
                raise ValueError(
                    f"{source}:{lineno}: expected a JSON object, "
                    f"got {type(obj).__name__}"
                )
            yield obj


def split_path(out_dir: Path | str, code: str) -> Path:
    """Path of the per-language file :func:`split_raw_dump` writes."""
    return Path(out_dir) / SPLIT_FILE_TEMPLATE.format(code=code)


def lang_code_prefilter(codes: Iterable[str]) -> re.Pattern[bytes]:
    """Byte pattern matching ``"lang_code": "<code>"`` for any of *codes*.

    Whitespace around the colon is optional, so dumps written with
    compact separators match too. A match is necessary but not
    sufficient (nested translations carry ``lang_code`` as well).
    """
    alternatives = b"|".join(re.escape(code.encode("ascii")) for code in codes)
    return re.compile(rb'"lang_code"\s*:\s*"(?:' + alternatives + rb')"')


def split_raw_dump(
    raw_gz: Path | str, out_dir: Path | str, codes: Iterable[str]
) -> dict[str, int]:
    """Copy the entries of the given languages out of the raw dump.

    Writes ``out_dir/en_wiktionary_<code>.jsonl`` for every code (an
    empty file when no entry matches) and returns the number of
    entries written per code. Lines are copied unchanged, in dump
    order.

    Each output is written to a ``.part`` file and renamed only after
    the whole dump has been read, so an interrupted run never leaves
    a truncated file under the final name.
    """
    wanted = tuple(dict.fromkeys(codes))
    if not wanted:
        raise ValueError("split_raw_dump needs at least one language code")
    for code in wanted:
        if not _LANG_CODE_RE.fullmatch(code):
            raise ValueError(f"invalid language code: {code!r}")
    prefilter = lang_code_prefilter(wanted)
    source = Path(raw_gz)
    target_dir = Path(out_dir)
    target_dir.mkdir(parents=True, exist_ok=True)
    final = {code: split_path(target_dir, code) for code in wanted}
    partial = {
        code: path.with_name(path.name + ".part")
        for code, path in final.items()
    }
    counts = dict.fromkeys(wanted, 0)
    try:
        with contextlib.ExitStack() as stack:
            handles: dict[str, BinaryIO] = {
                code: stack.enter_context(partial[code].open("wb"))
                for code in wanted
            }
            dump = stack.enter_context(_open_binary(source))
            for lineno, line in enumerate(dump, 1):
                if not prefilter.search(line):
                    continue
                entry = _parse_line(line, source, lineno)
                lang = (
                    entry.get("lang_code") if isinstance(entry, dict) else None
                )
                if not isinstance(lang, str) or lang not in handles:
                    continue
                if not line.endswith(b"\n"):
                    line += b"\n"
                handles[lang].write(line)
                counts[lang] += 1
    except BaseException:
        for path in partial.values():
            path.unlink(missing_ok=True)
        raise
    for code in wanted:
        os.replace(partial[code], final[code])
    return counts


_ETYMOLOGY_LABEL_RE: Final[re.Pattern[str]] = re.compile(r"\d+(?:\.\d+)*")


def etymology_label(entry: Mapping[str, Any]) -> str:
    """The etymology number as written: ``"0"`` when absent, ``"2"``,
    or a dotted sub-etymology such as ``"1.2"`` (tagliata).

    Raises :class:`ValueError` for any other shape.
    """
    raw = entry.get("etymology_number")
    if raw is None:
        return "0"
    if isinstance(raw, int) and not isinstance(raw, bool):
        return str(raw)
    text = raw.strip() if isinstance(raw, str) else None
    if text == "":
        return "0"
    if text is not None and text.isascii():
        if _ETYMOLOGY_LABEL_RE.fullmatch(text):
            return text
    raise ValueError(f"unexpected etymology_number: {raw!r}")


def entry_key(entry: Mapping[str, Any]) -> tuple[str, int, str]:
    """Return the lexeme key ``(word, etymology_number, pos)``.

    ``etymology_number`` is 0 for pages with a single, unnumbered
    etymology; numbered etymologies start at 1. A dotted
    sub-etymology number ("1.2") raises :class:`ValueError`; use
    :func:`etymology_label` where those can occur.
    """
    raw = entry.get("etymology_number")
    text = raw.strip() if isinstance(raw, str) else None
    if raw is None or text == "":
        number = 0
    elif isinstance(raw, int) and not isinstance(raw, bool):
        number = raw
    elif text is not None and text.isascii() and text.isdigit():
        number = int(text)
    else:
        raise ValueError(f"unexpected etymology_number: {raw!r}")
    return (str(entry.get("word", "")), number, str(entry.get("pos", "")))


def _is_form_of_sense(sense: Mapping[str, Any]) -> bool:
    if FORM_OF_TAGS.intersection(sense.get("tags") or ()):
        return True
    return "form_of" in sense or "alt_of" in sense


def is_form_of_only(entry: Mapping[str, Any]) -> bool:
    """True if every sense is tagged ``form-of`` or ``alt-of``.

    Such entries (``pie`` 'feminine plural of pio') describe forms of
    another lexeme and are not lexemes themselves. An entry without
    senses is not form-of-only: nothing marks it as one.
    """
    senses = [s for s in entry.get("senses") or () if isinstance(s, Mapping)]
    return bool(senses) and all(_is_form_of_sense(s) for s in senses)


def is_lemma_entry(entry: Mapping[str, Any]) -> bool:
    """True if at least one sense is a lexeme sense.

    A lexeme sense is neither form-of/alt-of (see
    :func:`is_form_of_only`) nor a ``no-gloss`` stub. Category names
    are not consulted: "<Language> lemmas" is a page category that
    leaks onto the form-of entries of the same page.
    """
    return any(
        isinstance(sense, Mapping)
        and not _is_form_of_sense(sense)
        and "no-gloss" not in (sense.get("tags") or ())
        for sense in entry.get("senses") or ()
    )


@dataclass(frozen=True, slots=True)
class Form:
    """One inflected form of an entry (:func:`inflected_forms`).

    ``tags`` are wiktextract's tags without the error tags, which are
    reported in ``tag_error`` (``|``-joined, sorted; "" if none).
    ``source`` is where wiktextract found the form ("declension",
    "inflection", or "" for the headword line).
    """

    form: str
    tags: tuple[str, ...]
    raw_tags: tuple[str, ...] = ()
    source: str = ""
    tag_error: str = ""

    def labels(self) -> frozenset[str]:
        """Tags and raw tags together (qualifiers spelled like tags)."""
        return frozenset(self.tags) | frozenset(self.raw_tags)


def _strings(value: Any) -> tuple[str, ...]:
    return tuple(item for item in value or () if isinstance(item, str))


def inflected_forms(entry: Mapping[str, Any]) -> list[Form]:
    """Every word form of the entry, in source order.

    Pseudo-forms (:data:`PSEUDO_FORM_TAGS`) and placeholder cells are
    dropped; forms carrying :data:`TAG_ERROR_TAGS` are kept, flagged.
    """
    out: list[Form] = []
    for form in entry.get("forms") or ():
        if not isinstance(form, Mapping):
            continue
        text = form.get("form")
        if not isinstance(text, str) or text.strip() in PLACEHOLDER_FORMS:
            continue
        tags = _strings(form.get("tags"))
        if PSEUDO_FORM_TAGS.intersection(tags):
            continue
        errors = sorted(TAG_ERROR_TAGS.intersection(tags))
        source = form.get("source")
        out.append(
            Form(
                form=text.strip(),
                tags=tuple(t for t in tags if t not in TAG_ERROR_TAGS),
                raw_tags=_strings(form.get("raw_tags")),
                source=source if isinstance(source, str) else "",
                tag_error="|".join(errors),
            )
        )
    return out


def plural_forms(entry: Mapping[str, Any]) -> list[Form]:
    """Forms tagged ``plural`` without an evaluative label.

    Evaluative labels (:data:`EVALUATIVE_TAGS`) are looked up in tags
    and raw tags alike. Flagged forms (``tag_error``) are included.
    """
    return [
        form
        for form in inflected_forms(entry)
        if "plural" in form.tags and not form.labels() & EVALUATIVE_TAGS
    ]


def forms_with_tags(
    entry: Mapping[str, Any],
    required: Iterable[str],
    excluded: Iterable[str] = (),
) -> list[dict[str, Any]]:
    """Forms whose tags include all of *required* and none of *excluded*.

    ``forms_with_tags(e, {"plural"}, {"archaic"})`` gives the plurals
    of an Italian noun minus archaic variants. Exclusion also matches
    a form's ``raw_tags`` (qualifiers wiktextract could not map to a
    tag, such as "with syntactic gemination after the verb"), so a
    qualifier spelled like a tag is honoured; requirements match
    ``tags`` only.

    Table metadata (:data:`META_FORM_TAGS`) and placeholder cells
    (:data:`PLACEHOLDER_FORMS`) are skipped. The returned dicts are
    copies that keep every key (``tags``, ``raw_tags``, ``source``…),
    in source order.
    """
    need = frozenset(required)
    avoid = frozenset(excluded)
    out: list[dict[str, Any]] = []
    for form in entry.get("forms") or ():
        if not isinstance(form, Mapping):
            continue
        text = form.get("form")
        if not isinstance(text, str) or text.strip() in PLACEHOLDER_FORMS:
            continue
        tags = frozenset(form.get("tags") or ())
        if tags & META_FORM_TAGS:
            continue
        labels = tags | frozenset(form.get("raw_tags") or ())
        if need <= tags and not labels & avoid:
            out.append(dict(form))
    return out


def ipa_list(entry: Mapping[str, Any]) -> list[str]:
    """IPA transcriptions of the entry, in source order, deduplicated.

    Strings are returned as written, delimiters included: ``/…/`` is
    phonemic and ``[…]`` phonetic, which callers may need to know.
    """
    seen: dict[str, None] = {}
    for sound in entry.get("sounds") or ():
        if not isinstance(sound, Mapping):
            continue
        ipa = sound.get("ipa")
        if isinstance(ipa, str) and ipa.strip():
            seen.setdefault(ipa.strip(), None)
    return list(seen)


def etymology_templates(
    entry: Mapping[str, Any], names: Iterable[str] | None = None
) -> list[dict[str, Any]]:
    """Etymology templates (``name``/``args``/``expansion`` dicts).

    With *names*, only templates with one of those names, such as
    ``{"inh", "bor", "lbor", "der"}`` or ``{"af", "suffix"}``. Order
    is source order; wiktextract sometimes repeats a template (an
    ``ety`` wrapper and its expansion), and the repeats are kept.
    """
    wanted = None if names is None else frozenset(names)
    out: list[dict[str, Any]] = []
    for template in entry.get("etymology_templates") or ():
        if not isinstance(template, Mapping):
            continue
        if wanted is None or template.get("name") in wanted:
            out.append(dict(template))
    return out


def head_template_args(
    entry: Mapping[str, Any], name: str
) -> dict[str, str] | None:
    """Arguments of the first head template called *name*, else None.

    ``head_template_args(e, "it-noun")`` returns ``{"1": "m", "2":
    "antichi"}`` for *antico* 'ancient one'; an explicit plural
    argument is what separates hand-entered plurals from the ones
    ``Module:it-common`` generates. ``{}`` means the template was used
    without arguments; None means it was not used.
    """
    for template in entry.get("head_templates") or ():
        if not isinstance(template, Mapping) or template.get("name") != name:
            continue
        args = template.get("args")
        if not isinstance(args, Mapping):
            return {}
        return {str(k): str(v) for k, v in args.items()}
    return None


ARGS_SEPARATOR: Final[str] = ";"
ARGS_ESCAPE: Final[str] = "\\"
HEAD_NAME_SEPARATOR: Final[str] = ":"


def _args_order(key: str) -> tuple[int, int, str]:
    return (0, int(key), "") if key.isdigit() else (1, 0, key)


def args_text(args: Mapping[str, Any]) -> str:
    """Template arguments as one compact string, as written.

    ``key=value`` pairs joined by ``;``: numbered arguments first in
    numeric order, then named ones in code-point order. A backslash
    and a ``;`` inside a value are escaped with a backslash, so
    :func:`parse_args_text` gives the arguments back. ``{"1": "m",
    "2": "amici"}`` → ``1=m;2=amici``; no arguments → "".
    """
    parts: list[str] = []
    for key in sorted((str(k) for k in args), key=_args_order):
        value = str(args[key]).replace(ARGS_ESCAPE, ARGS_ESCAPE * 2)
        value = value.replace(ARGS_SEPARATOR, ARGS_ESCAPE + ARGS_SEPARATOR)
        parts.append(f"{key}={value}")
    return ARGS_SEPARATOR.join(parts)


def parse_args_text(text: str) -> dict[str, str]:
    """The arguments of :func:`args_text` (a value keeps any ``=``)."""
    args: dict[str, str] = {}
    current: list[str] = []
    items: list[str] = []
    index = 0
    while index < len(text):
        char = text[index]
        if char == ARGS_ESCAPE and index + 1 < len(text):
            current.append(text[index : index + 2])
            index += 2
            continue
        if char == ARGS_SEPARATOR:
            items.append("".join(current))
            current = []
        else:
            current.append(char)
        index += 1
    if text:
        items.append("".join(current))
    for item in items:
        key, _, value = item.partition("=")
        unescaped: list[str] = []
        position = 0
        while position < len(value):
            if value[position] == ARGS_ESCAPE and position + 1 < len(value):
                unescaped.append(value[position + 1])
                position += 2
                continue
            unescaped.append(value[position])
            position += 1
        args[key] = "".join(unescaped)
    return args


def head_args(entry: Mapping[str, Any], default: str = "") -> list[str]:
    """Every head template of *entry* with its arguments, as written.

    One item per template, in order: ``name:`` plus :func:`args_text`
    of its arguments, the name left out for the part of speech's
    standard template *default* (``it-noun|1=m`` would repeat what the
    ``pos`` column says): ``["1=m;2=amici"]`` for *amico*,
    ``["head:1=it;2=noun;g=f"]`` for a generic head. An item names its
    template exactly when the text before its first ``:`` holds no
    ``=`` (template names have none; argument keys come first).
    """
    out: list[str] = []
    for template in entry.get("head_templates") or ():
        if not isinstance(template, Mapping):
            continue
        name = str(template.get("name") or "")
        raw = template.get("args")
        text = args_text(raw if isinstance(raw, Mapping) else {})
        out.append(
            text if name == default else f"{name}{HEAD_NAME_SEPARATOR}{text}"
        )
    return out


def _category_groups(entry: Mapping[str, Any]) -> list[Iterable[Any]]:
    groups: list[Iterable[Any]] = [entry.get("categories") or ()]
    groups.extend(
        sense.get("categories") or ()
        for sense in entry.get("senses") or ()
        if isinstance(sense, Mapping)
    )
    return groups


def categories(
    entry: Mapping[str, Any],
    prefix: str | None = None,
    *,
    lang_name: str | None = None,
) -> list[str]:
    """Category names of the entry and its senses, deduplicated.

    Entry-level categories come first, then sense-level ones, each in
    source order. With *prefix*, only names starting with it (for
    example ``"Italian learned borrowings"``). With *lang_name*, only
    names starting with ``lang_name + " "``: that language's own
    categories, without topical ``it:…`` ones or those of another
    language that leak onto the page ("Upper Sorbian …").
    """
    language_prefix = None if lang_name is None else lang_name + " "
    seen: dict[str, None] = {}
    for group in _category_groups(entry):
        for category in group:
            name = (
                category.get("name")
                if isinstance(category, Mapping)
                else category
            )
            if not isinstance(name, str) or not name:
                continue
            if prefix is not None and not name.startswith(prefix):
                continue
            if language_prefix is not None and not name.startswith(
                language_prefix
            ):
                continue
            seen.setdefault(name, None)
    return list(seen)


def assert_raw_format(entry: Mapping[str, Any]) -> None:
    """Raise :class:`ValueError` unless *entry* is a raw wiktextract record.

    kaikki's post-processed per-language files differ from the raw
    dump: their senses carry an ``id`` and their categories are dicts.
    The pipeline is written against the raw format, so such records
    are refused rather than half-understood.
    """
    word = entry.get("word")
    for sense in entry.get("senses") or ():
        if isinstance(sense, Mapping) and "id" in sense:
            raise ValueError(
                f"{word!r}: senses carry 'id' (a post-processed kaikki "
                "file); read the raw dump"
            )
    for group in _category_groups(entry):
        for category in group:
            if not isinstance(category, str):
                raise ValueError(
                    f"{word!r}: non-string category {category!r} (a "
                    "post-processed kaikki file); read the raw dump"
                )
