"""Stage ``extract``: one streaming pass per kaikki source (R7).

Inputs
  - the English-Wiktionary Italian split, ``en_wiktionary_it.jsonl``
    (top-level ``lang_code == "it"`` lines of the raw dump). The
    reference split under ``data/sources/kaikki/`` is used when present;
    otherwise the raw dump is split once into ``data/italian/kaikki/``
    (:func:`lexpipe.kaikki.split_raw_dump`);
  - the Italian-edition extract ``it-extract.jsonl.gz`` (all languages;
    only ``lang_code == "it"`` is kept).

Every en record must be in raw wiktextract format
(:func:`lexpipe.kaikki.assert_raw_format`); a post-processed record
stops the stage.

Outputs (``data/italian/extract/``, gitignored, source-stamped, see
:mod:`lexpipe.langs.italian.store`)
  - ``italian_en_lexemes.jsonl.gz`` — the noun and adjective lemma
    entries (``pos`` noun/adj and :func:`lexpipe.kaikki.is_lemma_entry`),
    projected to the fields later stages use, each with its
    ``lexeme_id``. The set is the one ``tests/lexpipe/test_ids.py``
    checks, and identifiers are exactly what
    :func:`lexpipe.ids.assign_lexeme_ids` gives for the whole set
    (R4): the function only looks inside groups of equal base
    identifiers, so it is called on the colliding groups, whose raw
    lines are re-read by file offset instead of holding every entry in
    memory. Sorted by ``lexeme_id``.
  - ``italian_en_verbs.jsonl.gz`` — verb lemma entries (head
    arguments, canonical accented infinitive, IPA of the infinitive,
    etymology, first gloss) with a ``verb_id`` built the same way;
    their conjugated forms are in the form index. Sorted by
    ``verb_id``.
  - ``italian_it_entries.jsonl.gz`` — Italian-edition noun and
    adjective entries (lemma and form-of): entry tags (gender, number,
    invariable), forms, IPA, accented hyphenation, sense tags.
    Sorted by (word, pos, record JSON).
  - ``italian_form_index.tsv.gz`` — form → lemma: one row per
    (form, lemma, pos, etymology, tags, source) of the en split:

        key       lowercase match key (respelling accents stripped)
        key_cs    the same with case kept
        form      the form as written (NFC)
        lemma     the lemma it belongs to
        pos       the part of speech of the entry
        ety       the etymology label of the lemma's entry ("" when the
                  row comes from a form-of entry, whose number is
                  page-local and unrelated)
        tags      wiktextract tags, sorted, "|"-joined
        source    lemma (the headword itself) | canonical (accented
                  headword) | head (headword-line form) | table
                  (inflection or conjugation table) | form_of | alt_of

    Used for ``form_ambiguous`` (carichi is also a form of caricare)
    and by later stages for verb paradigms. Sorted, deduplicated.
  - ``italian_extract_summary.tsv`` — entry counts.
"""

from __future__ import annotations

import json
import unicodedata
from collections import Counter, defaultdict
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

from lexpipe import ids, kaikki
from lexpipe.context import RunContext, run_stage, sha256_file
from lexpipe.langs.italian import keys, store
from lexpipe.langs.italian.layout import Layout

STAGE: Final[str] = "extract"
# 2: verb projections carry the infinitive's IPA (sounds), for the
# root-stressed -ere class of the ``verbs`` stage.
VERSION: Final[int] = 2
LEXEME_POS: Final[tuple[str, ...]] = ("noun", "adj")
VERB_POS: Final[str] = "verb"
IT_EDITION_POS: Final[tuple[str, ...]] = ("noun", "adj")

EN_LEXEMES_FORMAT: Final[str] = "italian_en_lexemes/1"
EN_VERBS_FORMAT: Final[str] = "italian_en_verbs/1"
IT_ENTRIES_FORMAT: Final[str] = "italian_it_entries/1"
FORM_INDEX_FORMAT: Final[str] = "italian_form_index/1"
FORM_INDEX_FIELDS: Final[tuple[str, ...]] = (
    "key",
    "key_cs",
    "form",
    "lemma",
    "pos",
    "ety",
    "tags",
    "source",
)
_TABLE_SOURCES: Final[frozenset[str]] = frozenset(
    {"conjugation", "declension", "inflection"}
)
_SKIPPED_FORM_TAGS: Final[frozenset[str]] = kaikki.PSEUDO_FORM_TAGS - {
    "canonical"
}


@dataclass(frozen=True, slots=True)
class _Pending:
    """A projected record waiting for its identifier."""

    base_id: str
    offset: int
    record: str


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


def _nfc(text: str) -> str:
    return unicodedata.normalize("NFC", text)


def _strings(value: Any) -> list[str]:
    if not isinstance(value, Sequence) or isinstance(value, str):
        return []
    return [item for item in value if isinstance(item, str)]


def _mappings(value: Any) -> list[Mapping[str, Any]]:
    if not isinstance(value, Sequence) or isinstance(value, str):
        return []
    return [item for item in value if isinstance(item, Mapping)]


def _args(template: Mapping[str, Any]) -> dict[str, str]:
    args = template.get("args")
    if not isinstance(args, Mapping):
        return {}
    return {str(k): str(v) for k, v in args.items()}


def _templates(entry: Mapping[str, Any], key: str) -> list[dict[str, Any]]:
    return [
        {"name": str(t.get("name") or ""), "args": _args(t)}
        for t in _mappings(entry.get(key))
    ]


def _words(value: Any) -> list[str]:
    return [
        _nfc(str(item["word"]))
        for item in _mappings(value)
        if isinstance(item.get("word"), str) and item["word"].strip()
    ]


# ---------------------------------------------------------------------------
# Projections
# ---------------------------------------------------------------------------


def _project_forms(entry: Mapping[str, Any]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for form in _mappings(entry.get("forms")):
        text = form.get("form")
        if not isinstance(text, str):
            continue
        source = form.get("source")
        out.append(
            {
                "form": _nfc(text.strip()),
                "tags": _strings(form.get("tags")),
                "raw_tags": _strings(form.get("raw_tags")),
                "source": source if isinstance(source, str) else "",
            }
        )
    return out


def _project_sounds(entry: Mapping[str, Any]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for sound in _mappings(entry.get("sounds")):
        ipa = sound.get("ipa")
        if not isinstance(ipa, str) or not ipa.strip():
            continue
        note = sound.get("note")
        out.append(
            {
                "ipa": _nfc(ipa.strip()),
                "tags": _strings(sound.get("tags")),
                "raw_tags": _strings(sound.get("raw_tags")),
                "note": note if isinstance(note, str) else "",
            }
        )
    return out


def _project_senses(entry: Mapping[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "glosses": _strings(sense.get("glosses")),
            "tags": _strings(sense.get("tags")),
            "raw_tags": _strings(sense.get("raw_tags")),
            "topics": _strings(sense.get("topics")),
            "form_of": _words(sense.get("form_of")),
            "alt_of": _words(sense.get("alt_of")),
        }
        for sense in _mappings(entry.get("senses"))
    ]


def project_lexeme(entry: Mapping[str, Any]) -> dict[str, Any]:
    """The fields of an en noun/adjective entry later stages use."""
    text = entry.get("etymology_text")
    return {
        "word": _nfc(str(entry.get("word") or "")),
        "pos": str(entry.get("pos") or ""),
        "etymology_number": kaikki.etymology_label(entry),
        "head_templates": _templates(entry, "head_templates"),
        "inflection_templates": _templates(entry, "inflection_templates"),
        "forms": _project_forms(entry),
        "sounds": _project_sounds(entry),
        "hyphenation": [_nfc(h) for h in _strings(entry.get("hyphenation"))],
        "etymology_text": text if isinstance(text, str) else "",
        "etymology_templates": _templates(entry, "etymology_templates"),
        "senses": _project_senses(entry),
        "categories": kaikki.categories(entry, lang_name="Italian"),
        "derived": _words(entry.get("derived")),
        "related": _words(entry.get("related")),
    }


def project_verb(entry: Mapping[str, Any]) -> dict[str, Any]:
    """The fields of an en verb lemma entry later stages use."""
    text = entry.get("etymology_text")
    canonical = [
        _nfc(str(form["form"]).strip())
        for form in _mappings(entry.get("forms"))
        if isinstance(form.get("form"), str)
        and "canonical" in _strings(form.get("tags"))
    ]
    return {
        "word": _nfc(str(entry.get("word") or "")),
        "pos": str(entry.get("pos") or ""),
        "etymology_number": kaikki.etymology_label(entry),
        "head_templates": _templates(entry, "head_templates"),
        "canonical": canonical,
        "sounds": _project_sounds(entry),
        "first_gloss": ids.first_gloss(entry),
        "etymology_text": text if isinstance(text, str) else "",
        "etymology_templates": _templates(entry, "etymology_templates"),
        "categories": kaikki.categories(entry, lang_name="Italian"),
    }


def project_it_entry(entry: Mapping[str, Any]) -> dict[str, Any]:
    """The fields of an Italian-edition noun/adjective entry."""
    hyphenations = [
        _strings(item.get("parts"))
        for item in _mappings(entry.get("hyphenations"))
    ]
    senses = [
        {
            "tags": _strings(sense.get("tags")),
            "raw_tags": _strings(sense.get("raw_tags")),
            "gloss": next(iter(_strings(sense.get("glosses"))), ""),
        }
        for sense in _mappings(entry.get("senses"))
    ]
    title = entry.get("pos_title")
    return {
        "word": _nfc(str(entry.get("word") or "")),
        "pos": str(entry.get("pos") or ""),
        "pos_title": title if isinstance(title, str) else "",
        "tags": _strings(entry.get("tags")),
        "raw_tags": _strings(entry.get("raw_tags")),
        "forms": [
            {
                "form": _nfc(str(f["form"]).strip()),
                "tags": _strings(f.get("tags")),
                "raw_tags": _strings(f.get("raw_tags")),
            }
            for f in _mappings(entry.get("forms"))
            if isinstance(f.get("form"), str)
        ],
        "ipa": [s["ipa"] for s in _project_sounds(entry)],
        "hyphenations": hyphenations,
        "senses": senses,
        "etymology_texts": _strings(entry.get("etymology_texts")),
    }


# ---------------------------------------------------------------------------
# Form index rows
# ---------------------------------------------------------------------------


def _index_row(
    form: str, lemma: str, pos: str, ety: str, tags: Sequence[str], source: str
) -> str | None:
    text = _nfc(form.strip())
    if text in kaikki.PLACEHOLDER_FORMS:
        return None
    key, key_cs = keys.form_keys(text)
    cells = (key, key_cs, text, lemma, pos, ety, "|".join(sorted(set(tags))))
    return "\t".join(store.tsv_cell(cell) for cell in (*cells, source))


def _form_source(form: Mapping[str, Any], tags: Sequence[str]) -> str:
    if "canonical" in tags:
        return "canonical"
    source = form.get("source")
    return "table" if source in _TABLE_SOURCES else "head"


def index_rows(entry: Mapping[str, Any], lemma_entry: bool) -> Iterator[str]:
    """Form-index rows of one en entry (see the module docstring)."""
    word = _nfc(str(entry.get("word") or ""))
    pos = str(entry.get("pos") or "")
    if not word or not pos:
        return
    ety = kaikki.etymology_label(entry) if lemma_entry else ""
    rows: list[str | None] = []
    if lemma_entry:
        rows.append(_index_row(word, word, pos, ety, (), "lemma"))
        for form in _mappings(entry.get("forms")):
            text = form.get("form")
            tags = _strings(form.get("tags"))
            if not isinstance(text, str) or _SKIPPED_FORM_TAGS & set(tags):
                continue
            source = _form_source(form, tags)
            clean = [t for t in tags if t not in kaikki.TAG_ERROR_TAGS]
            rows.append(_index_row(text, word, pos, ety, clean, source))
    for sense in _mappings(entry.get("senses")):
        tags = [
            t for t in _strings(sense.get("tags")) if t not in ("form-of",)
        ]
        for kind, field in (("form_of", "form_of"), ("alt_of", "alt_of")):
            for target in _words(sense.get(field)):
                rows.append(_index_row(word, target, pos, "", tags, kind))
    yield from (row for row in rows if row is not None)


# ---------------------------------------------------------------------------
# Identifiers
# ---------------------------------------------------------------------------


def _read_entries_at(path: Path, offsets: Sequence[int]) -> list[Any]:
    out: list[Any] = []
    with path.open("rb") as fh:
        for offset in offsets:
            fh.seek(offset)
            out.append(json.loads(fh.readline()))
    return out


def assign_ids(
    lang: str, split: Path, pending: Sequence[_Pending]
) -> list[str]:
    """Identifiers of *pending*, equal to ``ids.assign_lexeme_ids`` over
    the raw entries (see the module docstring)."""
    groups: dict[str, list[int]] = defaultdict(list)
    for index, item in enumerate(pending):
        groups[item.base_id].append(index)
    assigned = [item.base_id for item in pending]
    collided = sorted(
        index
        for members in groups.values()
        if len(members) > 1
        for index in members
    )
    if collided:
        offsets = [pending[index].offset for index in collided]
        entries = _read_entries_at(split, offsets)
        for index, ident in zip(
            collided, ids.assign_lexeme_ids(lang, entries), strict=True
        ):
            assigned[index] = ident
    if len(set(assigned)) != len(assigned):
        raise ValueError("lexeme_id assignment left duplicates")
    return assigned


def _identified(
    lang: str, split: Path, pending: Sequence[_Pending], id_field: str
) -> list[str]:
    """Serialized records with their identifiers, sorted by them."""
    identifiers = assign_ids(lang, split, pending)
    out: list[tuple[str, str]] = []
    for ident, item in zip(identifiers, pending, strict=True):
        record = json.loads(item.record)
        record[id_field] = ident
        out.append((ident, store.dumps(record)))
    out.sort()
    return [line for _, line in out]


# ---------------------------------------------------------------------------
# The passes
# ---------------------------------------------------------------------------


@dataclass
class _EnPass:
    lexemes: list[_Pending]
    verbs: list[_Pending]
    index: list[str]
    counts: Counter[str]


def scan_en_split(lang: str, split: Path) -> _EnPass:
    """Stream the en split once (see the module docstring)."""
    result = _EnPass([], [], [], Counter())
    offset = 0
    with split.open("rb") as fh:
        for lineno, line in enumerate(fh, 1):
            start = offset
            offset += len(line)
            if not line.strip():
                continue
            try:
                entry = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"{split}:{lineno}: invalid JSON") from exc
            if not isinstance(entry, dict):
                raise ValueError(f"{split}:{lineno}: not a JSON object")
            kaikki.assert_raw_format(entry)
            result.counts["entries"] += 1
            if entry.get("lang_code") != lang:
                result.counts["other_lang_code"] += 1
                continue
            _scan_entry(lang, entry, start, result)
    return result


def _scan_entry(
    lang: str, entry: dict[str, Any], offset: int, result: _EnPass
) -> None:
    pos = str(entry.get("pos") or "")
    lemma = kaikki.is_lemma_entry(entry)
    group = pos if pos in (*LEXEME_POS, VERB_POS, "name") else "other"
    result.counts[f"{'lemma' if lemma else 'nonlemma'}:{group}"] += 1
    result.index.extend(index_rows(entry, lemma))
    if not lemma or not entry.get("word"):
        return
    if pos in LEXEME_POS:
        record = store.dumps(project_lexeme(entry))
        base = ids.lexeme_id(lang, entry)
        result.lexemes.append(_Pending(base, offset, record))
    elif pos == VERB_POS:
        record = store.dumps(project_verb(entry))
        base = ids.lexeme_id(lang, entry)
        result.verbs.append(_Pending(base, offset, record))


def scan_it_extract(lang: str, path: Path) -> tuple[list[str], Counter[str]]:
    """Italian-edition noun/adjective projections, sorted."""
    counts: Counter[str] = Counter()
    keyed: list[tuple[str, str, str]] = []
    for entry in kaikki.iter_jsonl(path):
        counts["entries"] += 1
        if entry.get("lang_code") != lang:
            continue
        counts["lang_entries"] += 1
        pos = str(entry.get("pos") or "")
        if pos not in IT_EDITION_POS:
            continue
        projected = project_it_entry(entry)
        line = store.dumps(projected)
        counts[f"kept:{pos}"] += 1
        keyed.append((projected["word"], pos, line))
    keyed.sort()
    return [line for _, _, line in keyed], counts


def _sorted_unique(rows: list[str]) -> Iterator[str]:
    rows.sort()
    previous: str | None = None
    for row in rows:
        if row != previous:
            yield row
        previous = row


def _summary(counts: Mapping[str, int]) -> str:
    lines = ["key\tvalue"]
    lines.extend(f"{key}\t{counts[key]}" for key in sorted(counts))
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# The stage
# ---------------------------------------------------------------------------


def ensure_split(ctx: RunContext, layout: Layout) -> Path:
    """The en split to read, splitting the raw dump if there is none."""
    split = layout.en_split
    if split.is_file():
        return split
    raw = layout.raw_dump
    if not raw.is_file():
        raise FileNotFoundError(
            f"neither {layout.reference_split} nor {raw} exists; run the "
            "fetch stage first"
        )
    ctx.check_writable(layout.local_split)
    print(f"  splitting {raw} (lang_code {ctx.language.kaikki_code!r})")
    kaikki.split_raw_dump(
        raw, layout.local_split_dir, [ctx.language.kaikki_code]
    )
    return layout.local_split


def _body(ctx: RunContext, layout: Layout, split: Path) -> int:
    lang = ctx.language.kaikki_code
    stamps = {
        "en_split": sha256_file(split),
        "it_extract": sha256_file(layout.it_extract),
    }
    en = scan_en_split(lang, split)
    lexemes = _identified(lang, split, en.lexemes, "lexeme_id")
    verbs = _identified(lang, split, en.verbs, "verb_id")
    it_lines, it_counts = scan_it_extract(lang, layout.it_extract)
    counts: Counter[str] = Counter(en.counts)
    counts.update({f"it_extract:{k}": v for k, v in it_counts.items()})
    counts["lexemes"] = len(lexemes)
    counts["verbs"] = len(verbs)

    def meta(fmt: str, *names: str) -> dict[str, Any]:
        return {"format": fmt, "sources": {n: stamps[n] for n in names}}

    store.write_jsonl(
        ctx, layout.en_lexemes, meta(EN_LEXEMES_FORMAT, "en_split"), lexemes
    )
    store.write_jsonl(
        ctx, layout.en_verbs, meta(EN_VERBS_FORMAT, "en_split"), verbs
    )
    store.write_jsonl(
        ctx, layout.it_entries, meta(IT_ENTRIES_FORMAT, "it_extract"), it_lines
    )
    counts["form_index_rows"] = store.write_tsv(
        ctx,
        layout.form_index,
        meta(FORM_INDEX_FORMAT, "en_split"),
        FORM_INDEX_FIELDS,
        _sorted_unique(en.index),
    )
    ctx.write_text(layout.extract_summary, _summary(counts))
    print(
        f"  {counts['lexemes']:,} noun/adjective lexemes, "
        f"{counts['verbs']:,} verbs, "
        f"{counts['form_index_rows']:,} form-index rows, "
        f"{it_counts['lang_entries']:,} it-edition Italian entries"
    )
    return 0


def outputs(layout: Layout) -> list[Path]:
    return [
        layout.en_lexemes,
        layout.en_verbs,
        layout.it_entries,
        layout.form_index,
        layout.extract_summary,
    ]


def run(ctx: RunContext) -> int | None:
    """Run the ``extract`` stage (skipped when its manifest matches)."""
    layout = Layout(ctx)
    split = ensure_split(ctx, layout)
    params = {
        "version": VERSION,
        "kaikki_code": ctx.language.kaikki_code,
        "lexeme_pos": list(LEXEME_POS),
        "it_edition_pos": list(IT_EDITION_POS),
    }
    return run_stage(
        ctx,
        STAGE,
        inputs=[split, layout.it_extract],
        params=params,
        outputs=outputs(layout),
        body=lambda: _body(ctx, layout, split),
    )
