"""Tests for lexpipe.kaikki on synthetic wiktextract lines.

The raw-dump fixture mimics the real serialisation (``json.dumps``
with ``", "``/``": "`` separators) and contains decoys: entries of
other languages whose nested translations or descendants carry the
wanted ``lang_code``. A split by substring alone would keep them.
"""

from __future__ import annotations

import gzip
import itertools
import json
import sys
from pathlib import Path
from typing import Any

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lexpipe import kaikki  # noqa: E402

SPLIT_DIR = PROJECT_ROOT / "data" / "sources" / "kaikki"

ITALIAN: dict[str, Any] = {
    "word": "antico",
    "pos": "adj",
    "etymology_number": "2",
    "lang": "Italian",
    "lang_code": "it",
    "senses": [{"glosses": ["anterior"]}],
}
ITALIAN_PLURAL_FORM: dict[str, Any] = {
    "word": "pie",
    "pos": "adj",
    "lang": "Italian",
    "lang_code": "it",
    "senses": [{"glosses": ["f pl of pio"], "tags": ["form-of", "plural"]}],
}
CZECH: dict[str, Any] = {
    "word": "biolog",
    "pos": "noun",
    "lang": "Czech",
    "lang_code": "cs",
    "senses": [{"glosses": ["biologist"]}],
}
ENGLISH_DECOY: dict[str, Any] = {
    "word": "friend",
    "pos": "noun",
    "lang": "English",
    "lang_code": "en",
    "translations": [
        {"lang": "Italian", "lang_code": "it", "word": "amico"},
        {"lang": "Czech", "lang_code": "cs", "word": "přítel"},
    ],
}
LATIN_DECOY: dict[str, Any] = {
    "word": "amicus",
    "pos": "noun",
    "lang": "Latin",
    "lang_code": "la",
    "descendants": [{"lang_code": "it", "text": "Italian: amico"}],
}


def _line(entry: dict[str, Any], compact: bool = False) -> bytes:
    separators = (",", ":") if compact else (", ", ": ")
    text = json.dumps(entry, ensure_ascii=False, separators=separators)
    return text.encode("utf-8") + b"\n"


@pytest.fixture
def raw_dump(tmp_path: Path) -> Path:
    """A gzipped all-language dump; the last line has no newline."""
    lines = [
        _line(ENGLISH_DECOY),
        _line(ITALIAN),
        _line(LATIN_DECOY),
        _line(CZECH),
        b"\n",
        _line(ITALIAN_PLURAL_FORM, compact=True).rstrip(b"\n"),
    ]
    path = tmp_path / "raw-wiktextract-data.jsonl.gz"
    with gzip.open(path, "wb") as fh:
        fh.write(b"".join(lines))
    return path


# ---------------------------------------------------------------------------
# Streaming
# ---------------------------------------------------------------------------


def test_iter_jsonl_reads_plain_and_gzip_alike(tmp_path: Path) -> None:
    payload = _line(ITALIAN) + b"\n" + _line(CZECH)
    plain = tmp_path / "entries.jsonl"
    plain.write_bytes(payload)
    packed = tmp_path / "entries.jsonl.gz"
    packed.write_bytes(gzip.compress(payload))
    expected = [ITALIAN, CZECH]
    assert list(kaikki.iter_jsonl(plain)) == expected
    assert list(kaikki.iter_jsonl(packed)) == expected


def test_iter_jsonl_detects_gzip_by_content(tmp_path: Path) -> None:
    misnamed = tmp_path / "entries.jsonl"
    misnamed.write_bytes(gzip.compress(_line(CZECH)))
    assert list(kaikki.iter_jsonl(misnamed)) == [CZECH]


def test_iter_jsonl_names_the_bad_line(tmp_path: Path) -> None:
    path = tmp_path / "broken.jsonl"
    path.write_bytes(_line(ITALIAN) + b'{"word": \n')
    with pytest.raises(ValueError, match=r"broken\.jsonl:2: invalid JSON"):
        list(kaikki.iter_jsonl(path))


def test_iter_jsonl_rejects_non_objects(tmp_path: Path) -> None:
    path = tmp_path / "list.jsonl"
    path.write_bytes(b"[1, 2]\n")
    with pytest.raises(ValueError, match="expected a JSON object"):
        list(kaikki.iter_jsonl(path))


# ---------------------------------------------------------------------------
# Splitting the raw dump
# ---------------------------------------------------------------------------


def test_split_keeps_only_top_level_language(
    raw_dump: Path, tmp_path: Path
) -> None:
    out_dir = tmp_path / "split"
    counts = kaikki.split_raw_dump(raw_dump, out_dir, ["it", "cs"])
    assert counts == {"it": 2, "cs": 1}

    italian = kaikki.split_path(out_dir, "it").read_bytes()
    # Byte-for-byte copies of the dump lines; the unterminated last
    # line gains its newline.
    assert italian == _line(ITALIAN) + _line(ITALIAN_PLURAL_FORM, compact=True)
    czech = list(kaikki.iter_jsonl(kaikki.split_path(out_dir, "cs")))
    assert czech == [CZECH]
    assert not list(out_dir.glob("*.part"))


def test_split_accepts_any_json_spacing(tmp_path: Path) -> None:
    spaced = json.dumps(ITALIAN, ensure_ascii=False).replace(
        '"lang_code": "it"', '"lang_code" :\t "it"'
    )
    assert '"lang_code" :' in spaced
    dump = tmp_path / "raw.jsonl"
    dump.write_bytes(spaced.encode("utf-8") + b"\n" + _line(ENGLISH_DECOY))
    counts = kaikki.split_raw_dump(dump, tmp_path / "out", ["it"])
    assert counts == {"it": 1}


def test_lang_code_prefilter() -> None:
    pattern = kaikki.lang_code_prefilter(["it", "cs"])
    assert pattern.search(b'{"lang_code": "it"}')
    assert pattern.search(b'{"lang_code":"cs"}')
    assert pattern.search(b'{"lang_code"  :  "it", "x": 1}')
    assert not pattern.search(b'{"lang_code": "ita"}')
    assert not pattern.search(b'{"lang_code": "en", "x": "it"}')


def test_split_writes_empty_file_for_absent_language(
    raw_dump: Path, tmp_path: Path
) -> None:
    counts = kaikki.split_raw_dump(raw_dump, tmp_path, ["fr", "it", "fr"])
    assert counts == {"fr": 0, "it": 2}
    assert kaikki.split_path(tmp_path, "fr").read_bytes() == b""


def test_split_rejects_unsafe_codes(raw_dump: Path, tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="invalid language code"):
        kaikki.split_raw_dump(raw_dump, tmp_path, ["../it"])
    with pytest.raises(ValueError, match="at least one"):
        kaikki.split_raw_dump(raw_dump, tmp_path, [])


def test_split_leaves_nothing_behind_on_error(tmp_path: Path) -> None:
    dump = tmp_path / "raw.jsonl.gz"
    dump.write_bytes(
        gzip.compress(_line(ITALIAN) + b'{"lang_code": "it", "word": \n')
    )
    out_dir = tmp_path / "split"
    with pytest.raises(ValueError, match=r"raw\.jsonl\.gz:2"):
        kaikki.split_raw_dump(dump, out_dir, ["it"])
    assert list(out_dir.iterdir()) == []


# ---------------------------------------------------------------------------
# Entry helpers
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("number", "expected"),
    [("2", 2), (" 3 ", 3), (4, 4), (None, 0), ("", 0)],
)
def test_entry_key(number: Any, expected: int) -> None:
    entry = {"word": "antico", "pos": "adj", "etymology_number": number}
    assert kaikki.entry_key(entry) == ("antico", expected, "adj")


def test_entry_key_without_etymology_number() -> None:
    assert kaikki.entry_key(CZECH) == ("biolog", 0, "noun")


@pytest.mark.parametrize("number", ["II", "²", True, 1.5])
def test_entry_key_rejects_odd_numbers(number: Any) -> None:
    with pytest.raises(ValueError, match="etymology_number"):
        kaikki.entry_key(
            {"word": "x", "pos": "noun", "etymology_number": number}
        )


def test_is_form_of_only() -> None:
    assert kaikki.is_form_of_only(ITALIAN_PLURAL_FORM)
    assert not kaikki.is_form_of_only(ITALIAN)
    mixed = {
        "senses": [
            {"tags": ["form-of"], "form_of": [{"word": "pio"}]},
            {"glosses": ["pious woman"]},
        ]
    }
    assert not kaikki.is_form_of_only(mixed)
    alternative = {
        "senses": [{"tags": ["alt-of"]}, {"alt_of": [{"word": "x"}]}]
    }
    assert kaikki.is_form_of_only(alternative)
    assert not kaikki.is_form_of_only({"senses": []})
    assert not kaikki.is_form_of_only({})


CZECH_NOUN: dict[str, Any] = {
    "word": "biolog",
    "forms": [
        {"form": "bioložka", "tags": ["feminine"]},
        {"form": "animate", "source": "declension", "tags": ["table-tags"]},
        {"form": "cs-ndecl", "tags": ["inflection-template"]},
        {"form": "hard", "tags": ["class"]},
        {"form": "biologové", "tags": ["nominative", "plural"]},
        {"form": "biologa", "tags": ["genitive", "singular"]},
        {"form": "-", "tags": ["nominative", "plural"]},
        {
            "form": "biologi",
            "tags": ["nominative", "plural"],
            "raw_tags": ["proscribed"],
        },
        {"form": "biolozích", "tags": ["locative", "plural"]},
        {"form": "biologu", "tags": ["dated", "locative", "singular"]},
    ],
}


def test_forms_with_tags_selects_by_tag_set() -> None:
    nom_pl = kaikki.forms_with_tags(CZECH_NOUN, {"nominative", "plural"})
    assert [f["form"] for f in nom_pl] == ["biologové", "biologi"]
    # Qualifiers survive in the returned dicts.
    assert nom_pl[1]["raw_tags"] == ["proscribed"]

    plurals = kaikki.forms_with_tags(CZECH_NOUN, ["plural"], ["proscribed"])
    assert [f["form"] for f in plurals] == ["biologové", "biolozích"]

    loc_sg = kaikki.forms_with_tags(CZECH_NOUN, ["locative"], ["dated"])
    assert [f["form"] for f in loc_sg] == ["biolozích"]


def test_forms_with_tags_skips_table_metadata() -> None:
    every_form = kaikki.forms_with_tags(CZECH_NOUN, ())
    texts = [f["form"] for f in every_form]
    assert "animate" not in texts
    assert "cs-ndecl" not in texts
    assert "hard" not in texts
    assert "-" not in texts
    assert texts[0] == "bioložka"


def test_forms_with_tags_returns_copies() -> None:
    form = kaikki.forms_with_tags(CZECH_NOUN, {"feminine"})[0]
    form["form"] = "changed"
    assert CZECH_NOUN["forms"][0]["form"] == "bioložka"


def test_ipa_list() -> None:
    entry = {
        "sounds": [
            {"ipa": "/anˈti.ko/"},
            {"rhymes": "-iko"},
            {"ipa": "[anˈtiːko]", "tags": ["phonetic"]},
            {"ipa": "/anˈti.ko/"},
            {"ipa": "  "},
            "junk",
        ]
    }
    assert kaikki.ipa_list(entry) == ["/anˈti.ko/", "[anˈtiːko]"]
    assert kaikki.ipa_list({}) == []


ETYMOLOGY: dict[str, Any] = {
    "etymology_templates": [
        {"name": "ety", "args": {"1": "it", "2": ":lbor"}, "expansion": ""},
        {
            "name": "lbor",
            "args": {"1": "it", "2": "la", "3": "antīcus"},
            "expansion": "Learned borrowing from Latin antīcus",
        },
        {"name": "af", "args": {"1": "it", "2": "antico", "3": "-ità"}},
    ]
}


def test_etymology_templates() -> None:
    names = [t["name"] for t in kaikki.etymology_templates(ETYMOLOGY)]
    assert names == ["ety", "lbor", "af"]
    learned = kaikki.etymology_templates(ETYMOLOGY, {"lbor", "bor"})
    assert len(learned) == 1
    assert learned[0]["args"]["3"] == "antīcus"
    assert kaikki.etymology_templates({}, ["af"]) == []


def test_head_template_args() -> None:
    entry = {
        "head_templates": [
            {"name": "head", "args": {"1": "it", "2": "noun form"}},
            {"name": "it-noun", "args": {"1": "m", "2": "antichi"}},
            {"name": "it-noun", "args": {"1": "m", "2": "antici"}},
            {"name": "it-adv"},
        ]
    }
    assert kaikki.head_template_args(entry, "it-noun") == {
        "1": "m",
        "2": "antichi",
    }
    assert kaikki.head_template_args(entry, "it-adv") == {}
    assert kaikki.head_template_args(entry, "it-adj") is None


def test_categories_merges_entry_and_sense_levels() -> None:
    entry = {
        "categories": [
            "Italian lemmas",
            {"name": "Italian learned borrowings from Latin", "kind": "lang"},
        ],
        "senses": [
            {"categories": ["Italian non-lemma forms", "Italian lemmas"]},
            {"glosses": ["no categories here"]},
        ],
    }
    assert kaikki.categories(entry) == [
        "Italian lemmas",
        "Italian learned borrowings from Latin",
        "Italian non-lemma forms",
    ]
    assert kaikki.categories(entry, "Italian learned") == [
        "Italian learned borrowings from Latin"
    ]
    assert kaikki.categories({}) == []


# ---------------------------------------------------------------------------
# Typed accessors
# ---------------------------------------------------------------------------


def test_is_lemma_entry() -> None:
    assert kaikki.is_lemma_entry(ITALIAN)
    assert not kaikki.is_lemma_entry(ITALIAN_PLURAL_FORM)
    stub = {"senses": [{"tags": ["no-gloss"]}]}
    assert not kaikki.is_lemma_entry(stub)
    mixed = {"senses": [{"tags": ["form-of"]}, {"glosses": ["pious"]}]}
    assert kaikki.is_lemma_entry(mixed)
    assert not kaikki.is_lemma_entry({})
    assert not kaikki.is_lemma_entry({"senses": ["junk"]})


ITALIAN_NOUN: dict[str, Any] = {
    "word": "pantalone",
    "forms": [
        {"form": "pantaloni", "tags": ["plural"]},
        {"form": "pantaloncìni", "tags": ["diminutive", "plural"]},
        {
            "form": "pantalonàcci",
            "tags": ["plural"],
            "raw_tags": ["pejorative"],
        },
        {"form": "pantalóne", "tags": ["canonical"]},
        {"form": "pantalunə", "tags": ["romanization"]},
        {"form": "pantalona", "tags": ["error-unknown-tag", "plural"]},
        {"form": "–", "tags": ["plural"]},
        {"form": " pantalonissimi ", "tags": ["superlative", "plural"]},
        {"form": "pantalonera", "source": "inflection", "tags": ["feminine"]},
    ],
}


def test_inflected_forms_drop_pseudo_forms_and_flag_errors() -> None:
    forms = kaikki.inflected_forms(ITALIAN_NOUN)
    assert [f.form for f in forms] == [
        "pantaloni",
        "pantaloncìni",
        "pantalonàcci",
        "pantalona",
        "pantalonissimi",
        "pantalonera",
    ]
    flagged = forms[3]
    assert flagged.tag_error == "error-unknown-tag"
    assert flagged.tags == ("plural",)
    assert forms[2].raw_tags == ("pejorative",)
    assert forms[5].source == "inflection"
    assert all(f.tag_error == "" for i, f in enumerate(forms) if i != 3)
    czech = [f.form for f in kaikki.inflected_forms(CZECH_NOUN)]
    assert "animate" not in czech and "cs-ndecl" not in czech
    assert "hard" not in czech and "-" not in czech


def test_plural_forms_exclude_evaluatives() -> None:
    plurals = kaikki.plural_forms(ITALIAN_NOUN)
    assert [f.form for f in plurals] == ["pantaloni", "pantalona"]
    assert plurals[1].tag_error == "error-unknown-tag"
    czech = [f.form for f in kaikki.plural_forms(CZECH_NOUN)]
    assert czech == ["biologové", "biologi", "biolozích"]


def test_categories_by_language_name() -> None:
    entry = {
        "categories": ["Czech lemmas", "cs:Chess", "Czechoslovak things"],
        "senses": [
            {"categories": ["Czech masculine nouns", "Upper Sorbian nouns"]}
        ],
    }
    assert kaikki.categories(entry, lang_name="Czech") == [
        "Czech lemmas",
        "Czech masculine nouns",
    ]
    assert kaikki.categories(entry, "Czech m", lang_name="Czech") == [
        "Czech masculine nouns"
    ]
    # The positional prefix keeps its old meaning.
    assert kaikki.categories(entry, "Czech") == [
        "Czech lemmas",
        "Czechoslovak things",
        "Czech masculine nouns",
    ]


def test_assert_raw_format() -> None:
    kaikki.assert_raw_format(ITALIAN)
    kaikki.assert_raw_format({"categories": ["Italian lemmas"]})
    with pytest.raises(ValueError, match="'id'"):
        kaikki.assert_raw_format({"word": "x", "senses": [{"id": "x-1"}]})
    with pytest.raises(ValueError, match="non-string category"):
        kaikki.assert_raw_format({"categories": [{"name": "Italian lemmas"}]})
    with pytest.raises(ValueError, match="non-string category"):
        kaikki.assert_raw_format({"senses": [{"categories": [{"name": "x"}]}]})


# ---------------------------------------------------------------------------
# Smoke test on the reference split (skipped when it is absent)
# ---------------------------------------------------------------------------


@pytest.mark.needs_data
@pytest.mark.parametrize("code", ["it", "cs"])
def test_reference_split_head_is_well_formed(code: str) -> None:
    path = kaikki.split_path(SPLIT_DIR, code)
    if not path.exists():
        pytest.skip(f"reference split not present: {path}")
    head = list(itertools.islice(kaikki.iter_jsonl(path), 2000))
    assert head, f"{path} is empty"
    for entry in head:
        assert entry["lang_code"] == code
        word, number, pos = kaikki.entry_key(entry)
        assert word and pos and number >= 0
        kaikki.is_form_of_only(entry)
        kaikki.forms_with_tags(entry, {"plural"})
        kaikki.ipa_list(entry)
        kaikki.categories(entry)
        kaikki.assert_raw_format(entry)
        kaikki.is_lemma_entry(entry)
        kaikki.plural_forms(entry)
        name = "Italian" if code == "it" else "Czech"
        assert all(
            c.startswith(name + " ")
            for c in kaikki.categories(entry, lang_name=name)
        )


@pytest.mark.parametrize(
    ("number", "expected"),
    [(None, "0"), ("", "0"), (" 2 ", "2"), (3, "3"), ("1.2", "1.2")],
)
def test_etymology_label(number: Any, expected: str) -> None:
    assert kaikki.etymology_label({"etymology_number": number}) == expected


@pytest.mark.parametrize("number", ["II", "1.", ".2", "²", True, 1.5])
def test_etymology_label_rejects_odd_numbers(number: Any) -> None:
    with pytest.raises(ValueError, match="etymology_number"):
        kaikki.etymology_label({"etymology_number": number})


# ---------------------------------------------------------------------------
# Compact head arguments (the tracked raw lexicons' head_args)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "args",
    [
        {},
        {"1": "m", "2": "amici"},
        {"10": "x", "2": "y=z", "1": "f", "dim": "a;b\\c"},
        {"head": "a ; b", "f": "x\\;y"},
    ],
)
def test_args_text_round_trip(args: dict[str, str]) -> None:
    text = kaikki.args_text(args)
    assert kaikki.parse_args_text(text) == args


def test_args_text_order_and_head_args() -> None:
    assert kaikki.args_text({"g": "f", "2": "b", "1": "a", "10": "c"}) == (
        "1=a;2=b;10=c;g=f"
    )
    entry = {
        "head_templates": [
            {"name": "it-noun", "args": {"1": "m", "2": "amici"}},
            {"name": "head", "args": {"1": "it", "2": "noun"}},
            "not a template",
        ]
    }
    assert kaikki.head_args(entry, "it-noun") == [
        "1=m;2=amici",
        "head:1=it;2=noun",
    ]
    assert kaikki.head_args(entry) == [
        "it-noun:1=m;2=amici",
        "head:1=it;2=noun",
    ]
    assert kaikki.head_args({}) == []
