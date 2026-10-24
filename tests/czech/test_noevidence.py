"""No-evidence tables (R21): base classes and the learned -ický split."""

from __future__ import annotations

import pytest

from lexpipe.langs.czech import noevidence as ne

pytestmark = pytest.mark.cs


def unit(
    lemma: str, suffix: str, base_class: str, *words: tuple[str, str]
) -> ne.Unit:
    item = ne.Unit(lemma, lemma, suffix, base_class, "unknown")
    for word, outcome in words:
        item.outcomes.add(outcome)
        item.derivatives.append((word, outcome, "unknown"))
    return item


@pytest.mark.parametrize(
    ("lemma", "suffix", "derivative", "learned"),
    [
        ("kybernetika", "-ský/-cký", "kybernetický", True),
        ("cynik", "-ský/-cký", "cynický", True),
        ("antibiotikum", "-ský/-cký", "antibiotický", True),
        ("lajdák", "-ský/-cký", "lajdácký", False),
        ("Uzbek", "-ský/-cký", "uzbecký", False),
        ("kritika", "-ka", "kritička", False),
    ],
)
def test_learned_ic(
    lemma: str, suffix: str, derivative: str, learned: bool
) -> None:
    assert ne.learned_ic(lemma, suffix, derivative) is learned


def test_tables_without_learned_ic() -> None:
    items = [
        unit("gotika", "-ský/-cký", ne.NO_EVIDENCE, ("gotický", "pal")),
        unit("vzduch", "-ní/-ný", ne.NON_ALTERNATING, ("vzdušný", "pal")),
        unit("kluk", "-ský/-cký", ne.ALTERNATING, ("klukovský", "plain")),
    ]
    pooled = ne.crosstab(items)
    assert pooled[ne.POOLED]["pal"] == 2
    without = ne.crosstab(items, ne._group(ne.WITHOUT_IC))
    assert without[ne.POOLED]["pal"] == 1
    assert without[ne.ALTERNATING]["plain"] == 1
    text = "\n".join(ne.suffix_tables(items))
    assert ne.WITHOUT_IC in text


def test_counterexamples_star_learned_ic() -> None:
    items = [
        unit("gotika", "-ský/-cký", ne.NO_EVIDENCE, ("gotický", "pal")),
        unit("vzduch", "-ní/-ný", ne.NON_ALTERNATING, ("vzdušný", "pal")),
        unit("kluk", "-ek", ne.ALTERNATING, ("klouček", "pal")),
    ]
    found = ne.counterexamples(items)
    assert found == {
        ("unknown", "unknown"): ["gotika → gotický*", "vzduch → vzdušný"]
    }
