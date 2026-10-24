"""Tests for lexpipe.attest on synthetic Leipzig words tables."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lexpipe import attest, spec  # noqa: E402
from lexpipe.attest import AttestationTable  # noqa: E402

WORDS = "\n".join(
    [
        "1\tamico\t10",
        "2\tl'amico\t4",
        "3\tL’amico\t2",  # typographic apostrophe, sentence-initial
        "4\tdell'amico\t3",
        "5\tAmico\t1",
        "6\tpo'\t5",  # po' is not a clitic: one token
        "7\tdev'essere\t2",  # nor is dev'
        "8\tc'è\t6",
        "9\tquest'anno\t3",
        "10\tl'ex-ministro\t1",  # the rest is not alphabetic
        "11\t,\t50",
        "12\tmedici\t8",
        "13\tMedici\t2",
        "14\tun'amica\t1",
    ]
)
ITALIAN = spec.LANGUAGES["it"]
CZECH = spec.LANGUAGES["cs"]


@pytest.fixture
def words(tmp_path: Path) -> Path:
    path = tmp_path / "x-words.txt"
    path.write_text(WORDS + "\n", encoding="utf-8")
    return path


@pytest.mark.parametrize(
    ("token", "expected"),
    [
        ("l'amico", ("l", "amico")),
        ("L’Italia", ("L", "Italia")),
        ("DELL'ANNO", ("DELL", "ANNO")),
        ("c'è", ("c", "è")),
        ("«dell'amico»", ("dell", "amico")),
        ("quest’anno", ("quest", "anno")),
        ("po'", ("po'",)),
        ("dev'essere", ("dev'essere",)),
        ("l'ex-ministro", ("l'ex-ministro",)),
        ("l'", ("l'",)),
        ("'amico", ("'amico",)),
        ("amico", ("amico",)),
    ],
)
def test_split_elided(token: str, expected: tuple[str, ...]) -> None:
    assert attest.split_elided(token, ITALIAN.elision_clitics) == expected


def test_no_clitics_no_split() -> None:
    assert attest.split_elided("l'amico", frozenset()) == ("l'amico",)
    assert attest.split_elided("l'amico", CZECH.elision_clitics) == (
        "l'amico",
    )


def test_counts_without_splitting(words: Path) -> None:
    table = AttestationTable.from_words_file(words)
    assert table.count("amico") == 11
    assert table.count("amico", case_sensitive=True) == 10
    assert table.count("Amico", case_sensitive=True) == 1
    assert table.count("l'amico") == 4  # the ASCII and ’ spellings differ
    assert table.total_tokens == 98 - 50


def test_italian_clitics_are_split(words: Path) -> None:
    table = AttestationTable.from_words_file(
        words, clitics=ITALIAN.elision_clitics
    )
    assert table.count("amico") == 10 + 4 + 2 + 3 + 1
    assert table.count("amico", case_sensitive=True) == 19
    assert table.count("AMICO") == 20
    assert table.count("l") == 6 and table.count("L", case_sensitive=True) == 2
    assert table.count("anno") == 3 and table.count("è") == 6
    assert table.count("amica") == 1
    assert table.count("po") == 5  # trimmed, not split
    assert table.count("dev'essere") == 2
    assert table.count("l'ex-ministro") == 1
    assert table.count("medici") == 10
    assert table.count("Medici", case_sensitive=True) == 2
    # Six tokens were split in two: 4 + 2 + 3 + 6 + 3 + 1 more tokens.
    assert table.total_tokens == 48 + 19
    assert table.per_million("amico") == pytest.approx(20e6 / 67)


def test_rates_and_thresholds() -> None:
    table = AttestationTable(
        {"A": 1}, {"a": 1, "b": 5, "c": 20}, total_tokens=10_000_000
    )
    assert attest.TAUS == ("dictionary", "ge1", "ge0.5pm", "ge2pm")
    assert table.per_million("b") == 0.5
    expected = {
        "a": [True, True, False, False],
        "b": [True, True, True, False],  # exactly 0.5 per million
        "c": [True, True, True, True],  # exactly 2 per million
        "zzz": [True, False, False, False],
    }
    for form, flags in expected.items():
        assert list(table.attested_by_tau(form).values()) == flags
    assert table.attested("A", "ge1", case_sensitive=True)
    assert not table.attested("a", "ge1", case_sensitive=True)
    with pytest.raises(ValueError, match="unknown threshold"):
        table.attested("a", "ge5pm")


def test_tables_are_frozen_copies() -> None:
    counts = {"a": 1}
    table = AttestationTable(counts, counts, total_tokens=1)
    counts["a"] = 99
    assert table.count("a") == 1
    with pytest.raises(TypeError):
        table.counts_folded["a"] = 2  # type: ignore[index]
    with pytest.raises(ValueError, match="non-empty"):
        AttestationTable({}, {}, total_tokens=0)


def _install_words(root: Path, info: spec.LanguageInfo) -> None:
    corpus = (
        root / "data" / info.name / "leipzig" / "corpora" / info.leipzig_id
    )
    inner = corpus / info.leipzig_id
    inner.mkdir(parents=True)
    (inner / f"{info.leipzig_id}-words.txt").write_text(
        WORDS + "\n", encoding="utf-8"
    )


def test_for_language_applies_the_token_rules(tmp_path: Path) -> None:
    for info in (ITALIAN, CZECH):
        _install_words(tmp_path, info)
    italian = AttestationTable.for_language(ITALIAN, root=tmp_path)
    czech = AttestationTable.for_language(CZECH, root=tmp_path)
    assert italian.count("amico") == 20
    assert czech.count("amico") == 11
    with pytest.raises(FileNotFoundError):
        AttestationTable.for_language(ITALIAN, root=tmp_path / "empty")


# ---------------------------------------------------------------------------
# attested_count: the threshold arithmetic for precomputed counts
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("count", "expected"),
    [
        (0, [True, False, False, False]),
        (1, [True, True, False, False]),
        (4, [True, True, False, False]),
        (5, [True, True, True, False]),  # exactly 0.5 per million
        (19, [True, True, True, False]),
        (20, [True, True, True, True]),  # exactly 2 per million
    ],
)
def test_attested_count(count: int, expected: list[bool]) -> None:
    total = 10_000_000
    assert [
        attest.attested_count(count, tau, total) for tau in attest.TAUS
    ] == expected
    assert list(attest.attested_count_by_tau(count, total).values()) == (
        expected
    )
    # The table applies the very same arithmetic.
    table = AttestationTable({}, {"form": count}, total_tokens=total)
    assert list(table.attested_by_tau("form").values()) == expected
    assert [table.attested_count(count, tau) for tau in attest.TAUS] == (
        expected
    )


def test_attested_count_is_exact_with_awkward_totals() -> None:
    # 0.5 per million of 3,000,001 tokens is 1.5000005 tokens.
    assert not attest.attested_count(1, "ge0.5pm", 3_000_001)
    assert attest.attested_count(2, "ge0.5pm", 3_000_001)
    # 2 per million of 1,500,000 tokens is exactly 3 tokens.
    assert attest.attested_count(3, "ge2pm", 1_500_000)
    assert not attest.attested_count(2, "ge2pm", 1_500_000)


@pytest.mark.parametrize(
    ("count", "tau", "total"),
    [(1, "ge5pm", 10), (-1, "ge1", 10), (1, "ge1", 0), (1, "dictionary", -5)],
)
def test_attested_count_validates(count: int, tau: str, total: int) -> None:
    with pytest.raises(ValueError):
        attest.attested_count(count, tau, total)


# ---------------------------------------------------------------------------
# attest_paradigm
# ---------------------------------------------------------------------------


@pytest.fixture
def paradigm_table() -> AttestationTable:
    cased = {"biolog": 30, "biologové": 4, "Biologové": 1, "biologa": 2}
    folded = {"biolog": 30, "biologové": 5, "biologa": 2, "biologů": 1}
    return AttestationTable(cased, folded, total_tokens=2_000_000)


def test_attest_paradigm_folded(paradigm_table: AttestationTable) -> None:
    result = attest.attest_paradigm(
        paradigm_table,
        ["biolog", "biologové", "Biologové", "biolozích", "", "biologa"],
    )
    # Biologové folds onto biologové: one form, counted once.
    assert [item.form for item in result.forms] == [
        "biolog",
        "biologové",
        "biolozích",
        "biologa",
    ]
    assert [item.count for item in result.forms] == [30, 5, 0, 2]
    assert result.total == 37 and result.total_tokens == 2_000_000
    biolozich = result.forms[2]
    assert biolozich.attested_by_tau == {
        "dictionary": True,
        "ge1": False,
        "ge0.5pm": False,
        "ge2pm": False,
    }
    # 2 tokens of 2M = 1 per million: ge0.5pm, not ge2pm.
    assert result.forms[3].attested_by_tau["ge0.5pm"] is True
    assert result.forms[3].attested_by_tau["ge2pm"] is False
    assert result.attested_by_tau == {
        "dictionary": True,
        "ge1": True,
        "ge0.5pm": True,
        "ge2pm": True,
    }
    assert result.attested_forms("ge2pm") == ("biolog", "biologové")
    assert result.n_counted == 4 and result.excluded == ()


def test_attest_paradigm_case_sensitive_and_excluded(
    paradigm_table: AttestationTable,
) -> None:
    result = attest.attest_paradigm(
        paradigm_table,
        ["biologové", "Biologové", "biologa"],
        case_sensitive=True,
        exclude={"biologa"},  # also a form of another lexeme
        taus=("ge1", "ge2pm"),
    )
    assert [(i.form, i.count, i.counted) for i in result.forms] == [
        ("biologové", 4, True),
        ("Biologové", 1, True),
        ("biologa", 2, False),
    ]
    assert result.total == 5 and result.excluded == ("biologa",)
    assert result.attested_by_tau == {"ge1": True, "ge2pm": True}
    assert set(result.forms[0].attested_by_tau) == {"ge1", "ge2pm"}
    # The excluded form's own attestation is still reported.
    assert result.forms[2].attested_by_tau["ge1"] is True
    assert result.attested_forms("ge1") == ("biologové", "Biologové")
    # A predicate works as well as a collection.
    same = attest.attest_paradigm(
        paradigm_table,
        ["biologové", "Biologové", "biologa"],
        case_sensitive=True,
        exclude=lambda form: form.endswith("a"),
        taus=("ge1", "ge2pm"),
    )
    assert same == result


def test_attest_paradigm_edge_cases(paradigm_table: AttestationTable) -> None:
    empty = attest.attest_paradigm(paradigm_table, [])
    assert empty.forms == () and empty.total == 0
    # Like any count, a dictionary paradigm is attested at "dictionary".
    assert empty.attested_by_tau["dictionary"] is True
    assert empty.attested_by_tau["ge1"] is False
    unseen = attest.attest_paradigm(paradigm_table, ["snězích"])
    assert unseen.total == 0 and not unseen.attested_by_tau["ge1"]
    with pytest.raises(ValueError, match="unknown threshold"):
        attest.attest_paradigm(paradigm_table, ["biolog"], taus=("ge3pm",))
    with pytest.raises(ValueError, match="unknown threshold"):
        unseen.attested_forms("ge3pm")
