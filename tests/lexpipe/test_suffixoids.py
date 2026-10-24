"""The shared -log/-fag/-urg classifiers (lexpipe.suffixoids, R27)."""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lexpipe import csvio, curated  # noqa: E402
from lexpipe import suffixoids as sx  # noqa: E402


@pytest.mark.parametrize(
    ("gloss", "expected"),
    [
        ("psychologist", ("person", "gloss")),
        ("dramaturge, playwright, dramatist", ("person", "gloss")),
        ("catalogue; list", ("thing", "gloss")),
        ("dialog", ("thing", "gloss")),
        ("ideologue", ("person", "gloss")),
        ("phytophage, herbivore", ("animate", "gloss")),
        ("bacteriophage", ("thing", "gloss")),
        # "cave" is a whole word: scavenger is an animal, not a cave
        ("cave", ("thing", "gloss")),
        ("scavenger", ("animate", "gloss")),
        ("necrophage, scavenger", ("animate", "gloss")),
        ("cannibal", ("person", "gloss")),
        ("", ("unknown", "no gloss")),
        ("zzz", ("unknown", "gloss")),
    ],
)
def test_semantic_class(gloss: str, expected: tuple[str, str]) -> None:
    assert sx.semantic_class(gloss) == expected


def test_only_the_first_gloss_is_read() -> None:
    assert sx.semantic_class("astrologer; synonym of astronomo")[0] == "person"
    assert sx.first_gloss(" a ; b ") == "a"


def test_decisions_and_extra_patterns() -> None:
    decision = sx.Decision(semantic_class="thing")
    assert sx.semantic_class("surgeon", decision=decision) == (
        "thing",
        "manual",
    )
    extra = sx.GlossPatterns(
        label="cs",
        person=re.compile(r"\blékař\b"),
        thing=re.compile(r"\bpřístroj\b"),
    )
    assert sx.semantic_class("lékař", extra=extra) == ("person", "gloss (cs)")
    assert sx.semantic_class("přístroj", extra=extra) == (
        "thing",
        "gloss (cs)",
    )
    # English rules come first.
    assert sx.semantic_class("surgeon", extra=extra) == ("person", "gloss")


@pytest.mark.parametrize(
    ("lemma", "suffixoid", "gloss", "evidence", "expected"),
    [
        ("psiholog", "log", "psychologist", "", (True, "gloss")),
        ("fitofag", "fag", "phytophage", "", (True, "gloss")),
        ("chirurg", "urg", "surgeon", "", (True, "gloss")),
        ("telegraf", "graf", "telegraph", "", (True, "gloss")),
        ("bibliofil", "fil", "bibliophile", "", (True, "gloss")),
        ("bârlog", "log", "cave", "", (False, "gloss")),
        (
            "politolog",
            "log",
            "political scientist",
            "",
            (True, "gloss (profession)"),
        ),
        ("ciotolog", "log", "small and ugly man", "", (False, "gloss")),
        ("drammaturgo", "urg", "playwright", "ending", (True, "ending")),
    ],
)
def test_learned_suffixoid(
    lemma: str,
    suffixoid: str,
    gloss: str,
    evidence: str,
    expected: tuple[bool, str],
) -> None:
    assert sx.learned_suffixoid(
        lemma, suffixoid, gloss, evidence=evidence
    ) == (expected)


def test_learned_suffixoid_order_and_errors() -> None:
    no = sx.Decision(learned=False)
    assert sx.learned_suffixoid(
        "teologo", "log", "theologian", decision=no, evidence="ending"
    ) == (False, "manual")
    only_class = sx.Decision(semantic_class="thing")
    assert sx.learned_suffixoid("dialog", "log", "", decision=only_class) == (
        False,
        "gloss",
    )
    extra = sx.GlossPatterns(label="cs", profession=re.compile(r"\blékař\b"))
    assert sx.learned_suffixoid("rinolog", "log", "lékař", extra=extra) == (
        True,
        "gloss (profession, cs)",
    )
    with pytest.raises(ValueError, match="unknown suffixoid"):
        sx.learned_suffixoid("x", "nom", "y")


def _write(path: Path, rows: list[dict[str, object]]) -> Path:
    fields = (
        "lexeme_id",
        "learned_suffixoid",
        "semantic_class",
        "decided_by",
        "note",
        "guard_word",
    )
    csvio.write_rows(path, rows, fields)
    return path


def test_load_decisions(tmp_path: Path) -> None:
    path = _write(
        tmp_path / sx.DECISIONS_FILE,
        [
            {
                "lexeme_id": "a",
                "learned_suffixoid": False,
                "semantic_class": "",
                "decided_by": "t",
                "note": "lookalike",
                "guard_word": "a",
            },
            {
                "lexeme_id": "b",
                "learned_suffixoid": None,
                "semantic_class": "person",
                "decided_by": "t",
                "note": "",
                "guard_word": "b",
            },
        ],
    )
    rows, decisions = sx.load_decisions(path)
    assert decisions["a"] == sx.Decision(False, "", "lookalike")
    assert decisions["b"] == sx.Decision(None, "person", "")
    curated.check_guards(rows, {"a": {"word": "a"}, "b": {"word": "b"}})
    assert sx.load_decisions(tmp_path / "absent.csv") == ({}, {})


def test_load_decisions_refuses_bad_values(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "bad.csv",
        [
            {
                "lexeme_id": "a",
                "learned_suffixoid": None,
                "semantic_class": "mixed",
                "decided_by": "t",
                "note": "",
                "guard_word": "a",
            }
        ],
    )
    with pytest.raises(ValueError, match="semantic_class"):
        sx.load_decisions(path)
    missing = tmp_path / "missing.csv"
    csvio.write_rows(missing, [{"lexeme_id": "a"}], ("lexeme_id",))
    with pytest.raises(ValueError, match="missing columns"):
        sx.load_decisions(missing)


@pytest.mark.parametrize(
    ("language", "key"),
    [
        ("italian", "lexeme_id"),
        ("czech", "lexeme_id"),
        ("romanian_compare", "lemma"),
    ],
)
def test_every_language_has_a_well_formed_decision_list(
    language: str, key: str
) -> None:
    path = PROJECT_ROOT / "curated" / language / sx.DECISIONS_FILE
    rows, decisions = sx.load_decisions(path, key)
    assert rows and set(rows) == set(decisions)
    for row in rows.values():
        assert row.guards, row.key
        assert row.values["decided_by"], row.key
