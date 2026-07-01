#!/usr/bin/env python3
"""Golden derivations from latex.tex.

Every case here corresponds to an entry in the paper's derivation
tables (e:coronal-derivation, latex.tex:472-492, and
e:cluster-derivation, latex.tex:584-596). The predicted surface form
after running the paper's revised rule set must reproduce the paper's
stated SR column-for-column.

If a test here fails, the analysis in
``phonology.analyses.romanian_palatalization`` disagrees with the
paper — either the analysis needs to change, or the paper's stated
SR is wrong, or both.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from phonology.analyses.romanian_palatalization import (
    RULES_FROM_PAPER,
    load_inventory,
)
from phonology.g2p import build_ur
from phonology.pipeline import RulePipeline
from phonology.segments import segments_to_ipa


@pytest.fixture(scope="module")
def pipeline() -> RulePipeline:
    inv = load_inventory()
    return RulePipeline(
        rules=RULES_FROM_PAPER,
        resolve=lambda label: inv.segment(label),
    )


@pytest.fixture(scope="module")
def inventory():
    return load_inventory()


def _normalize_glide(ipa: str) -> str:
    """Treat trailing /j/ as equivalent to /i/ for comparison.

    The paper writes the desyllabified glide as [ʲ] (a palatalisation
    mark on the preceding consonant). The lexicon's normalised IPA
    field writes it as the vowel /i/ instead — the two conventions
    describe the same surface fact. Our pipeline outputs /j/; this
    helper normalises both sides to trailing /i/ for a byte-equal
    comparison.
    """
    if ipa.endswith(("j", "ʲ")):
        return ipa[:-1] + "i"
    return ipa


def _derive(pipeline: RulePipeline, inv, lemma_ipa: str, opp: str) -> str:
    ur = build_ur(lemma_ipa, opp, inv)
    assert ur is not None, f"build_ur failed for {lemma_ipa!r}+{opp!r}"
    deriv = pipeline.derive(ur)
    return segments_to_ipa(deriv.sr, inv.base_segments())


# ---------------------------------------------------------------------------
# e:coronal-derivation (latex.tex:472-492)
# ---------------------------------------------------------------------------


class TestCoronalDerivation:
    """The paper's five worked examples for coronal palatalisation."""

    def test_rus_plus_i_becomes_rushi(self, pipeline, inventory) -> None:
        sr = _derive(pipeline, inventory, "rus", "i")
        # Paper: ruʃ + palatalisation-mark; data: ruʃi
        assert _normalize_glide(sr) == "ruʃi"

    def test_prost_plus_i_becomes_proshti(self, pipeline, inventory) -> None:
        sr = _derive(pipeline, inventory, "prost", "i")
        # Bleed protects /T/, [t] surfaces after /ʃ/.
        assert _normalize_glide(sr) == "proʃti"

    def test_brad_plus_i_becomes_brazi(self, pipeline, inventory) -> None:
        sr = _derive(pipeline, inventory, "brad", "i")
        # Assibilation supplies +Strident to /D/, coronal-default
        # supplies +Continuant → spirantises to [z].
        assert _normalize_glide(sr) == "brazi"

    def test_skut_plus_uri_is_inalterable(self, pipeline, inventory) -> None:
        # Prespecified /t/ + non-palatalising -uri: nothing mutates.
        sr = _derive(pipeline, inventory, "skut", "uri")
        assert _normalize_glide(sr) == "skuturi"

    def test_klopot_plus_i_verb_becomes_klopotsi(
        self, pipeline, inventory,
    ) -> None:
        # Assibilated /T/ stays -Continuant (prespecified) →
        # affricate [t͡s] surfaces, not [s].
        sr = _derive(pipeline, inventory, "klopot", "i")
        assert _normalize_glide(sr) == "klopot͡si"


# ---------------------------------------------------------------------------
# e:cluster-derivation (latex.tex:584-596)
# ---------------------------------------------------------------------------


class TestClusterDerivation:
    """Cluster cases confirm the terminator-breadth split and the
    place-class trigger."""

    def test_muskə_plus_e_becomes_mushte(self, pipeline, inventory) -> None:
        # /K/ palatalises before /e/ (broad terminator finds /e/,
        # condition holds). Derived [tʃ] triggers S-pal on /S/.
        # Bleed clears the derived affricate's place features;
        # coronal-default refills with plain-stop values → [ʃt].
        sr = _derive(pipeline, inventory, "muskə", "e")
        assert sr == "muʃte"

    def test_broaskə_plus_e_becomes_broashte(
        self, pipeline, inventory,
    ) -> None:
        sr = _derive(pipeline, inventory, "broaskə", "e")
        assert sr == "broaʃte"

    def test_punkt_plus_e_stays_punkte(self, pipeline, inventory) -> None:
        # Broad terminator halts on /T/, condition [+Syll,+Front]
        # fails → /K/ blocked; default fill supplies /k/.
        sr = _derive(pipeline, inventory, "punkt", "e")
        assert sr == "punkte"


# ---------------------------------------------------------------------------
# e:noun-plurals paper exemplars (latex.tex:68-76)
# ---------------------------------------------------------------------------


class TestNounPluralExemplars:
    """The paper's headline six palatalisation exemplars."""

    def test_kolak_plus_i_plural_becomes_kolatshi(
        self, pipeline, inventory,
    ) -> None:
        # UR koloK-i; dorsal-pal turns /K/ into [tʃ].
        sr = _derive(pipeline, inventory, "kolak", "i")
        assert _normalize_glide(sr) == "kolat͡ʃi"

    def test_pisikə_plus_i_becomes_pisitshi(
        self, pipeline, inventory,
    ) -> None:
        sr = _derive(pipeline, inventory, "pisikə", "i")
        assert _normalize_glide(sr) == "pisit͡ʃi"

    def test_milog_plus_i_becomes_milodzhi(
        self, pipeline, inventory,
    ) -> None:
        sr = _derive(pipeline, inventory, "miloɡ", "i")
        assert _normalize_glide(sr) == "milod͡ʒi"

    def test_munte_plus_i_becomes_muntsi(self, pipeline, inventory) -> None:
        # -e lemma with plural -i: strip thematic e, assibilate /T/.
        sr = _derive(pipeline, inventory, "munte", "i")
        assert _normalize_glide(sr) == "munt͡si"


# ---------------------------------------------------------------------------
# Inalterability under prespecification
# ---------------------------------------------------------------------------


class TestInalterability:
    """Prespecified obstruents survive palatalising environments
    unchanged (the paper's central claim)."""

    def test_uri_plural_does_not_palatalize_velar(
        self, pipeline, inventory,
    ) -> None:
        # fok-uri: /k/ prespecified -Coronal, cannot palatalise even
        # though the sequence [k]+[u] is present.
        sr = _derive(pipeline, inventory, "fok", "uri")
        assert _normalize_glide(sr) == "fokuri"

    def test_uri_plural_does_not_assibilate_stop(
        self, pipeline, inventory,
    ) -> None:
        sr = _derive(pipeline, inventory, "skut", "uri")
        assert _normalize_glide(sr) == "skuturi"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
