"""Tests for lexpipe.ids (R4).

The last test builds lexeme_ids for every noun and adjective lemma
entry of the reference kaikki splits (gitignored; skipped when
absent) and checks they are unique and independent of input order.
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lexpipe import ids, kaikki  # noqa: E402

SPLIT_DIR = PROJECT_ROOT / "data" / "sources" / "kaikki"

GRECO_PERSON: dict[str, Any] = {
    "word": "greco",
    "pos": "noun",
    "head_templates": [{"name": "it-noun", "args": {"1": "m", "2": "greci"}}],
    "senses": [{"glosses": ["Greek (person)"]}],
}
GRECO_LANGUAGE: dict[str, Any] = {
    "word": "greco",
    "pos": "noun",
    "head_templates": [{"name": "it-noun", "args": {"1": "m", "2": "-"}}],
    "senses": [{"glosses": ["Greek (language)"]}],
}


def test_fingerprint_is_canonical() -> None:
    a = {"b": [1, "è"], "a": {"y": 2, "x": 1}}
    b = {"a": {"x": 1, "y": 2}, "b": [1, "è"]}
    assert ids.canonical_json(a) == '{"a":{"x":1,"y":2},"b":[1,"è"]}'
    assert ids.fingerprint(a) == ids.fingerprint(b)
    expected = hashlib.sha1(ids.canonical_json(a).encode("utf-8"))
    assert ids.fingerprint(a) == expected.hexdigest()[:8]
    assert re.fullmatch(r"[0-9a-f]{8}", ids.fingerprint([]))
    assert ids.fingerprint([{"1": "m"}]) != ids.fingerprint([{"1": "f"}])


def test_lexeme_id_format() -> None:
    entry = {
        **GRECO_PERSON,
        "etymology_number": "2",
        "inflection_templates": [{"name": "x", "args": {"1": "a"}}],
    }
    fp = ids.fingerprint([{"1": "m", "2": "greci"}, {"1": "a"}])
    assert ids.lexeme_id("it", entry) == f"it:greco:noun:2:{fp}"
    assert ids.template_args(entry) == [{"1": "m", "2": "greci"}, {"1": "a"}]


def test_lexeme_id_separates_homographs_by_templates() -> None:
    person = ids.lexeme_id("it", GRECO_PERSON)
    language = ids.lexeme_id("it", GRECO_LANGUAGE)
    assert person != language
    assert person.startswith("it:greco:noun:0:")


def test_template_args_tolerate_odd_shapes() -> None:
    entry = {"head_templates": [{"name": "it-adv"}, "junk", {"args": [1]}]}
    assert ids.template_args(entry) == [{}, {}]
    assert ids.template_args({}) == []


def test_lexeme_id_needs_word_and_pos() -> None:
    with pytest.raises(ValueError, match="word or pos"):
        ids.lexeme_id("it", {"word": "x"})


def test_gloss_slug() -> None:
    assert ids.gloss_slug(GRECO_PERSON) == "greek-person"
    entry = {"senses": [{"tags": ["no-gloss"]}, {"glosses": ["Šach, mat!"]}]}
    assert ids.gloss_slug(entry) == "sach-mat"
    assert ids.gloss_slug({"senses": []}) == ids.NO_GLOSS
    long = {"senses": [{"glosses": ["word " * 30]}]}
    slug = ids.gloss_slug(long)
    assert len(slug) <= ids.GLOSS_SLUG_MAX and not slug.endswith("-")


def _twin(entry: Mapping[str, Any], **changes: Any) -> dict[str, Any]:
    return {**entry, **changes}


def test_assign_lexeme_ids_disambiguates_deterministically() -> None:
    base = {"word": "re", "pos": "noun", "head_templates": []}
    a = _twin(base, senses=[{"glosses": ["king"]}])
    b = _twin(base, senses=[{"glosses": ["re (musical note)"]}])
    c1 = _twin(base, senses=[{"glosses": ["king"]}], forms=[{"form": "re"}])
    unique = _twin(base, word="ra")
    entries = [a, b, c1, unique]
    assigned = ids.assign_lexeme_ids("it", entries)
    stem = ids.lexeme_id("it", a)
    assert assigned[1] == f"{stem}:re-musical-note"
    assert {assigned[0], assigned[2]} == {f"{stem}:king:1", f"{stem}:king:2"}
    assert assigned[3] == ids.lexeme_id("it", unique)
    # The same entries in another order get the same identifiers.
    shuffled = [unique, c1, b, a]
    again = ids.assign_lexeme_ids("it", shuffled)
    assert dict(zip(map(ids.canonical_json, shuffled), again)) == dict(
        zip(map(ids.canonical_json, entries), assigned)
    )


def test_root_id() -> None:
    def italian_key(entry: Mapping[str, Any]) -> list[str]:
        word = entry["word"]
        stem = word[:-1] if word[-1] in "oaie" else word
        return [stem, str(entry.get("etymology_number") or 0)]

    amico = {"word": "amico", "pos": "noun"}
    amica = {"word": "amica", "pos": "noun"}
    amico_adj = {"word": "amico", "pos": "adj"}
    antico_1 = {"word": "antico", "etymology_number": "1"}
    antico_2 = {"word": "antico", "etymology_number": "2"}
    rid = ids.root_id("it", amico, italian_key)
    assert rid == ids.root_id("it", amica, italian_key)
    assert rid == ids.root_id("it", amico_adj, italian_key)
    assert rid.startswith("it:amic:")
    assert ids.root_id("it", antico_1, italian_key) != ids.root_id(
        "it", antico_2, italian_key
    )
    with pytest.raises(ValueError):
        ids.root_id("it", amico, lambda entry: [])
    with pytest.raises(ValueError):
        ids.root_id("it", amico, lambda entry: ["", "x"])


# ---------------------------------------------------------------------------
# Uniqueness over the reference splits (gitignored data)
# ---------------------------------------------------------------------------

_NOUN_OR_ADJ = re.compile(rb'"pos"\s*:\s*"(?:noun|adj)"')


def _lemma_entries(path: Path) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    with path.open("rb") as fh:
        for line in fh:
            if not _NOUN_OR_ADJ.search(line):
                continue
            entry = json.loads(line)
            if entry.get("pos") in ("noun", "adj") and kaikki.is_lemma_entry(
                entry
            ):
                out.append(entry)
    return out


@pytest.mark.needs_data
@pytest.mark.slow
@pytest.mark.parametrize("code", ["it", "cs"])
def test_lexeme_ids_are_unique_on_the_reference_split(code: str) -> None:
    path = kaikki.split_path(SPLIT_DIR, code)
    if not path.is_file():
        pytest.skip(f"reference split not present: {path}")
    entries = _lemma_entries(path)
    assert len(entries) > 10_000
    assigned = ids.assign_lexeme_ids(code, entries)
    assert len(set(assigned)) == len(entries)
    # The fingerprint alone already separates almost every homograph:
    # the critique's 490 Italian / 245 Czech colliding keys.
    bases = [ids.lexeme_id(code, entry) for entry in entries]
    assert len(entries) - len(set(bases)) < 50
    # Independent of input order.
    backwards = ids.assign_lexeme_ids(code, entries[::-1])
    assert backwards[::-1] == assigned


def test_sub_etymology_numbers() -> None:
    entry = {"word": "tagliata", "pos": "noun", "etymology_number": "1.2"}
    assert ids.lexeme_id("it", entry).startswith("it:tagliata:noun:1.2:")
