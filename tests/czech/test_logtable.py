"""-log/-urg/-graf/-fil classes (R27): the shared classifiers of
lexpipe.suffixoids, the Czech evidence and gloss layer, and the curated
decisions of curated/czech/suffixoid_classes.csv (keyed by lexeme_id,
so a decision follows the reading: biograf 'cinema' is m.in)."""

from __future__ import annotations

from pathlib import Path

import pytest

from lexpipe import curated, suffixoids
from lexpipe.langs.czech import logtable as lt

pytestmark = pytest.mark.cs

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DECISIONS = PROJECT_ROOT / "curated" / "czech" / suffixoids.DECISIONS_FILE
CINEMA = suffixoids.Decision(learned=True, semantic_class="thing")


@pytest.mark.parametrize(
    ("gloss", "decision", "expected"),
    [
        ("biographer", None, ("person", "gloss")),
        ("cinema", CINEMA, ("thing", "manual")),
        ("cinema", None, ("unknown", "gloss")),
        ("biologist", None, ("person", "gloss")),
        ("dialog", None, ("thing", "gloss")),
        ("ideologue", None, ("person", "gloss")),
        (
            "lékař zabývající se léčbou cukrovky",
            None,
            ("person", "gloss (cs)"),
        ),
        (
            "přístroj na zaznamenání zemětřesných vln",
            None,
            ("thing", "gloss (cs)"),
        ),
    ],
)
def test_semantic_class(
    gloss: str,
    decision: suffixoids.Decision | None,
    expected: tuple[str, str],
) -> None:
    assert lt.semantic_class(gloss, decision) == expected


@pytest.mark.parametrize(
    ("lemma", "gloss", "evidence", "edition", "expected"),
    [
        ("biolog", "biologist", "orthographic", "en", (True, "gloss")),
        ("telegraf", "telegraph", "etymology", "en", (True, "etymology")),
        # cs.wiktionary glosses are Czech: the combining -o- decides
        (
            "germanofil",
            "milovník Němců",
            "",
            "cs",
            (True, lt.CS_COMPOUND),
        ),
        (
            "diabetolog",
            "lékař zabývající se léčbou cukrovky",
            "orthographic",
            "en",
            (True, "gloss (profession, cs)"),
        ),
        # profil 'profile' only looks like -fil
        (
            "profil",
            "profile (side view of a head)",
            "",
            "en",
            (False, "gloss"),
        ),
    ],
)
def test_learned_suffixoid(
    lemma: str,
    gloss: str,
    evidence: str,
    edition: str,
    expected: tuple[bool, str],
) -> None:
    found = lt.czech_evidence(lemma, evidence, edition)
    assert lt.learned_suffixoid(lemma, gloss, found) == expected


def test_curated_decisions_are_well_formed_and_guarded() -> None:
    rows, decisions = suffixoids.load_decisions(DECISIONS)
    assert decisions["cs:biograf:noun:0:45ac8625"] == suffixoids.Decision(
        True, "thing", decisions["cs:biograf:noun:0:45ac8625"].note
    )
    biograf = rows["cs:biograf:noun:0:45ac8625"]
    nouns = [
        {
            "lexeme_id": key,
            "lemma": row.guards["lemma"],
            "readings": row.guards["readings"],
            "first_gloss": row.guards["first_gloss"],
        }
        for key, row in rows.items()
    ]
    assert lt.load_decisions(DECISIONS, nouns) == decisions
    changed = [dict(n) for n in nouns]
    for noun in changed:
        if noun["lexeme_id"] == biograf.lexeme_id:
            noun["first_gloss"] = "a movie theater"
    with pytest.raises(curated.CuratedGuardError, match="first_gloss"):
        lt.load_decisions(DECISIONS, changed)
