"""DeriNet 2.3 parsing and the compact projection (R7, R14).

Rows are synthetic (DeriNet is CC BY-NC-SA and is not tracked); they
follow the real ten-column format.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lexpipe import spec
from lexpipe.context import RunContext
from lexpipe.langs.czech import derinet
from tests.czech import analysis_tree

pytestmark = pytest.mark.cs


def test_parse_line() -> None:
    line = analysis_tree.DERINET_TREES[1][2]  # knihovat
    node = derinet.parse_line(line)
    assert (node.lemma, node.upos, node.parent_id) == (
        "knihovat",
        "VERB",
        "2.0",
    )
    assert node.feats["ConjugClass"] == "3"
    assert node.morphs == (("knih", "R"), ("ova", "D"), ("t", "I"))
    assert node.relation_type == "Derivation"
    assert derinet.verbalizer(node) == "-ovat"


def test_bad_json_names_the_line() -> None:
    with pytest.raises(ValueError, match="line 7"):
        derinet.parse_line("1.0\tx#N\tx\tNOUN\t\t[{bad\t\t\t\t{}", 7)


@pytest.mark.parametrize(
    ("morphs", "label"),
    [
        ((("knih", "R"), ("ova", "D"), ("t", "I")), "-ovat"),
        ((("s", "D"), ("nož", "R"), ("i", "D"), ("t", "I")), "-it"),
        ((("křič", "R"), ("e", "D"), ("t", "I")), "-et"),
        ((("táh", "R"), ("nou", "D"), ("t", "I")), "-nout"),
        ((("dělat", "R"),), ""),
    ],
)
def test_verbalizer_from_last_d_morph(
    morphs: tuple[tuple[str, str], ...], label: str
) -> None:
    node = derinet.Node("1.1", "x", "VERB", {}, morphs, "1.0", {}, {})
    assert derinet.verbalizer(node) == label


def test_conj_class_is_only_a_check() -> None:
    assert derinet.conj_consistent("-it", "4") is True
    assert derinet.conj_consistent("-ovat", "3#2") is True
    assert derinet.conj_consistent("-ovat", "4") is False
    assert derinet.conj_consistent("-nout", "2") is None
    assert derinet.conj_consistent("-it", "") is None


def test_compatible_join() -> None:
    assert derinet.compatible(["m.in"], "m", "in")
    assert not derinet.compatible(["m.in"], "m", "an")  # šach × 2
    assert not derinet.compatible(["f"], "m", "in")
    assert derinet.compatible(["m"], "m", "an")  # animacy unknown
    assert derinet.compatible(["f"], "", "")  # no gender in DeriNet


def test_projection(tmp_path: Path) -> None:
    source = tmp_path / "derinet.tsv"
    source.write_bytes(analysis_tree._derinet_bytes(None))
    ctx = RunContext(spec.LANGUAGES["cs"], root=tmp_path)
    target = tmp_path / "data" / "czech" / "nc" / "derinet.sqlite"
    lemmas = frozenset({"kniha", "liga", "šach", "strach", "ruka"})
    counts = derinet.project(ctx, source, target, lemmas)
    assert counts.fictitious == 1
    nouns, edges = derinet.read_projection(target)
    pairs = {(e.parent_lemma, e.child_lemma, e.verbalizer) for e in edges}
    assert ("kniha", "knižní", "") in pairs
    assert ("kniha", "knihovat", "-ovat") in pairs
    assert ("ruka", "ručit", "-it") in pairs
    assert ("strach", "strašit", "-it") in pairs
    # strašný hangs under the verb strašit: not a noun-parented edge
    assert not any(e.child_lemma == "strašný" for e in edges)
    # the Fictitious affixoid's child edge is dropped
    assert not any(e.parent_lemma == "-log-" for e in edges)
    sach = sorted((n.gender, n.animacy) for n in nouns if n.lemma == "šach")
    assert sach == [("m", "an"), ("m", "in")]
    liga = next(n for n in nouns if n.lemma == "liga")
    assert liga.loan_root is True and liga.etymology == "deu"


def test_projection_is_deterministic(tmp_path: Path) -> None:
    ctx = RunContext(spec.LANGUAGES["cs"], root=tmp_path)
    digests = []
    for seed in (None, 3):
        source = tmp_path / f"derinet{seed}.tsv"
        source.write_bytes(analysis_tree._derinet_bytes(seed))
        target = tmp_path / "data" / "czech" / "nc" / f"p{seed}.sqlite"
        derinet.project(ctx, source, target, frozenset({"kniha", "liga"}))
        nouns, edges = derinet.read_projection(target)
        digests.append((nouns, edges))
    assert digests[0] == digests[1]
