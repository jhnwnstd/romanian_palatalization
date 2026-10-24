"""MorfFlex lemma grammar (R14), tags and the streaming filter.

The lemma strings are written inline (TR-64 notation); no MorfFlex
data (CC BY-NC-SA) is stored in the repository.
"""

from __future__ import annotations

import lzma
from pathlib import Path

import pytest

from lexpipe.langs.czech import morfflex

pytestmark = pytest.mark.cs


def test_ethnonym_and_surname_labels() -> None:
    lemma = morfflex.parse_pdt_lemma("Čech_;E_;Y")
    assert lemma.base == "Čech" and lemma.homonym is None
    assert lemma.name_labels == ("E", "Y") and lemma.is_proper
    assert not lemma.is_styled and lemma.unparsed == ""


def test_homonym_number() -> None:
    lemma = morfflex.parse_pdt_lemma("jazyk-2")
    assert (lemma.base, lemma.homonym) == ("jazyk", 2)
    assert morfflex.parse_pdt_lemma("a-33").homonym == 33


def test_style_label_and_comment() -> None:
    lemma = morfflex.parse_pdt_lemma("Abacus_;m_,s_^(^DD**Abakus)")
    assert lemma.base == "Abacus"
    assert lemma.name_labels == ("m",)
    assert lemma.style_labels == ("s",) and lemma.is_styled
    assert lemma.comment == "^DD**Abakus"


def test_digits_belong_to_the_base() -> None:
    lemma = morfflex.parse_pdt_lemma("Aa1_;m_^(stupeň_národního_ratingu)")
    assert lemma.base == "Aa1" and lemma.homonym is None
    assert lemma.comment == "stupeň_národního_ratingu"
    assert lemma.name_labels == ("m",)


def test_derivation_info() -> None:
    lemma = morfflex.parse_pdt_lemma("matčin_^(*3ka)")
    assert lemma.base == "matčin"
    assert lemma.comment == "*3ka" and lemma.derivation == "*3ka"
    other = morfflex.parse_pdt_lemma(
        "automatisovávaný_^(^DD**automatizovávaný)_(*2t)"
    )
    assert other.comment == "^DD**automatizovávaný"
    assert other.derivation == "*2t"


@pytest.mark.parametrize(
    ("raw", "base"),
    [
        ("Coca-Cola", "Coca-Cola"),
        ("kupovat_:T", "kupovat"),
        ("New_York_;G", "New_York"),
    ],
)
def test_hyphens_and_underscores_in_bases(raw: str, base: str) -> None:
    assert morfflex.parse_pdt_lemma(raw).base == base


def test_tags_and_variants() -> None:
    tag = morfflex.parse_tag("NNMP1-----A---1")
    assert (tag.gender, tag.number, tag.case, tag.negation, tag.var) == (
        "M",
        "P",
        "1",
        "A",
        "1",
    )
    assert tag.cell == "NOM.PL"
    assert morfflex.parse_tag(morfflex.ADJECTIVE_TAG).cell == "NOM.PL.M.AN"
    assert morfflex.parse_tag("NNFXX-----A----").cell is None
    with pytest.raises(ValueError, match="length"):
        morfflex.parse_tag("NNMP1")
    assert [morfflex.var_class(v) for v in "-1469a"] == [
        "standard",
        "standard",
        "standard",
        "nonstandard",
        "nonstandard",
        "abbreviation",
    ]


def test_filter_keeps_affirmative_nouns_and_the_adjective_cell(
    tmp_path: Path,
) -> None:
    rows = [
        "kluk\tNNMP1-----A----\tkluci",
        "kluk\tNNMP1-----N----\tnekluci",
        "český\tAAMP1----1A----\tčeští",
        "český\tAAMS1----1A----\tčeský",
        "kupovat_:T\tVB-S---3P-AA---\tkupuje",
        "Čech_;E_;Y\tNNMP1-----A---1\tČechové",
    ]
    path = tmp_path / "m.tsv.xz"
    path.write_bytes(lzma.compress(("\n".join(rows) + "\n").encode()))
    kept = list(morfflex.iter_rows(path, chunk_size=7))
    assert kept == [
        ("kluk", "NNMP1-----A----", "kluci"),
        ("český", "AAMP1----1A----", "čeští"),
        ("Čech_;E_;Y", "NNMP1-----A---1", "Čechové"),
    ]
    plain = tmp_path / "m.tsv"
    plain.write_text("\n".join(rows) + "\n", encoding="utf-8")
    assert list(morfflex.iter_rows(plain)) == kept
