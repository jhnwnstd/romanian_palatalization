"""Stable identifiers for lexemes and roots (design rule R4).

``(word, etymology_number, pos)`` does not identify a wiktextract
entry: the pinned subsets have 490 Italian and 245 Czech colliding
keys (greco the person vs the language; šach 'check' vs 'shah'), and
``etymology_number`` is mostly null and page-local. Identifiers are
therefore built as::

    lexeme_id = "{lang}:{word}:{pos}:{ety}:{fp}"

  - ``ety`` is the etymology number as written, 0 when absent
    (``1.2`` for a sub-etymology, :func:`lexpipe.kaikki.etymology_label`);
  - ``fp`` is :func:`fingerprint` of the list of every head
    template's ``args`` followed by every inflection template's
    ``args`` (an entry's headword and declension, which differ
    between homographs: ``greco`` ``{2: greci}`` vs ``{2: -}``).

On the pinned subsets this leaves 13 Italian and 1 Czech tie among
noun/adjective lemma entries. :func:`assign_lexeme_ids` resolves them
by appending the slug of the first gloss (``…:a70ddc82:ufo``) and,
for entries still tied, an ordinal (``…:1``, ``…:2``) in the order of
the entries' canonical JSON, so the result does not depend on the
order of the input. ``word`` may contain ``:``; the other fields
never do.

A lexeme_id changes when its entry's templates change or when a new
homograph appears in a later dump; curated files therefore store
guard fields next to it (``docs/runbooks.md``).

``root_id`` groups the lexemes a learner hears as one root: POS
conversions and gender pairs of one etymology (amico N/ADJ,
amico/amica) share it, different etymologies (antico 'ancient' vs
àntico 'anterior') do not. What counts as the same root is language
specific, so :func:`root_id` takes a key function from the plugin.

Canonical JSON (:func:`canonical_json`) is ``json.dumps`` with sorted
keys, no whitespace and non-ASCII characters kept, encoded as UTF-8.
"""

from __future__ import annotations

import hashlib
import json
import re
import unicodedata
from collections import defaultdict
from collections.abc import Callable, Mapping, Sequence
from typing import Any, Final, TypeAlias

from lexpipe.kaikki import etymology_label

FINGERPRINT_LENGTH: Final[int] = 8
GLOSS_SLUG_MAX: Final[int] = 40
NO_GLOSS: Final[str] = "nogloss"
_SLUG_SEPARATORS_RE: Final[re.Pattern[str]] = re.compile(r"[^a-z0-9]+")

RootKeyFn: TypeAlias = Callable[[Mapping[str, Any]], Sequence[str]]


def canonical_json(obj: Any) -> str:
    """JSON with sorted keys and no whitespace (see the module docstring)."""
    return json.dumps(
        obj,
        sort_keys=True,
        ensure_ascii=False,
        separators=(",", ":"),
        allow_nan=False,
    )


def fingerprint(obj: Any) -> str:
    """First 8 hex digits of the SHA-1 of :func:`canonical_json`."""
    digest = hashlib.sha1(
        canonical_json(obj).encode("utf-8"), usedforsecurity=False
    )
    return digest.hexdigest()[:FINGERPRINT_LENGTH]


def _args_of(templates: Any) -> list[Any]:
    out: list[Any] = []
    for template in templates or ():
        if not isinstance(template, Mapping):
            continue
        args = template.get("args")
        out.append(args if isinstance(args, Mapping) else {})
    return out


def template_args(entry: Mapping[str, Any]) -> list[Any]:
    """``head_templates[*].args + inflection_templates[*].args``."""
    return _args_of(entry.get("head_templates")) + _args_of(
        entry.get("inflection_templates")
    )


def lexeme_id(lang: str, entry: Mapping[str, Any]) -> str:
    """``{lang}:{word}:{pos}:{ety or 0}:{fp}`` of a wiktextract entry.

    Ties between entries are not resolved here; use
    :func:`assign_lexeme_ids` for a whole lexicon.
    """
    word = unicodedata.normalize("NFC", str(entry.get("word") or ""))
    pos = str(entry.get("pos") or "")
    if not word or not pos:
        raise ValueError(f"entry without word or pos: {word!r}/{pos!r}")
    ety = etymology_label(entry)
    return f"{lang}:{word}:{pos}:{ety}:{fingerprint(template_args(entry))}"


def first_gloss(entry: Mapping[str, Any]) -> str:
    """The first gloss of the first sense that has one ("" if none)."""
    for sense in entry.get("senses") or ():
        if not isinstance(sense, Mapping):
            continue
        for gloss in sense.get("glosses") or ():
            if isinstance(gloss, str) and gloss.strip():
                return gloss.strip()
    return ""


def gloss_slug(entry: Mapping[str, Any]) -> str:
    """ASCII slug of :func:`first_gloss`, at most 40 characters."""
    folded = unicodedata.normalize("NFKD", first_gloss(entry))
    ascii_text = "".join(ch for ch in folded if not unicodedata.combining(ch))
    slug = _SLUG_SEPARATORS_RE.sub("-", ascii_text.lower()).strip("-")
    return slug[:GLOSS_SLUG_MAX].rstrip("-") or NO_GLOSS


def assign_lexeme_ids(
    lang: str, entries: Sequence[Mapping[str, Any]]
) -> list[str]:
    """Unique lexeme_ids for *entries*, aligned with the input.

    Entries whose :func:`lexeme_id` ties get ``:<gloss slug>``; those
    still tied get ``:<n>``, numbered in the order of their canonical
    JSON (deterministic, independent of input order).
    """
    ids = [lexeme_id(lang, entry) for entry in entries]
    groups: dict[str, list[int]] = defaultdict(list)
    for index, base in enumerate(ids):
        groups[base].append(index)
    for base, members in groups.items():
        if len(members) < 2:
            continue
        by_slug: dict[str, list[int]] = defaultdict(list)
        for index in members:
            by_slug[f"{base}:{gloss_slug(entries[index])}"].append(index)
        for slugged, tied in by_slug.items():
            if len(tied) == 1:
                ids[tied[0]] = slugged
                continue
            ordered = sorted(tied, key=lambda i: canonical_json(entries[i]))
            for ordinal, index in enumerate(ordered, 1):
                ids[index] = f"{slugged}:{ordinal}"
    if len(set(ids)) != len(ids):
        raise ValueError("lexeme_id disambiguation left duplicates")
    return ids


def root_id(lang: str, entry: Mapping[str, Any], key: RootKeyFn) -> str:
    """``{lang}:{stem}:{fp}`` from a language-supplied root key.

    *key* maps an entry to a non-empty sequence of strings: a readable
    stem first, then whatever separates roots with that stem (an
    etymon, say). Entries with equal keys share the root_id; ``fp``
    is :func:`fingerprint` of the whole key.
    """
    parts = list(key(entry))
    if not parts or not all(isinstance(part, str) for part in parts):
        raise ValueError(f"root key must be non-empty strings: {parts!r}")
    stem = unicodedata.normalize("NFC", parts[0])
    if not stem:
        raise ValueError("root key starts with an empty stem")
    return f"{lang}:{stem}:{fingerprint(parts)}"
