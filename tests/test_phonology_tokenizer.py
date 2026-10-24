#!/usr/bin/env python3
"""Tests for the profile-level IPA tokenizer (design rule R8).

Pins three things:

- **Romanian default.** :data:`ROMANIAN_TOKENIZER` is the default of
  :func:`tokenize_ipa` and of ``AnalysisProfile.tokenizer``, and it
  tokenizes Romanian IPA exactly as ``tokenize_ipa`` did before it took
  a config. The expected tokens below were produced by that version:
  Romanian lexicon strings, plus edge cases whose output is a quirk
  worth keeping.
- **Other conventions.** A non-Romanian :class:`TokenizerConfig` keeps
  ː, ɛ, ɔ and multi-codepoint segments such as Czech r̝ as distinct
  tokens.
- **Threading.** The profile's tokenizer is the one the validator's
  distance metric splits attested IPA with.

Nothing here reads ``data/``.
"""

from __future__ import annotations

import dataclasses
import json
import sys
from pathlib import Path
from typing import Any

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import phonology.analyses.romanian_palatalization as romanian  # noqa: E402
from phonology.diagnostics.distance import feature_edit_distance  # noqa: E402
from phonology.inventory import FeatureInventory  # noqa: E402
from phonology.lexicon import ClusterTag, LexRow  # noqa: E402
from phonology.profile import AnalysisProfile  # noqa: E402
from phonology.segments import (  # noqa: E402
    LENGTH_MARK,
    NON_SEGMENTAL_MARKS,
    ROMANIAN_TOKENIZER,
    LengthHandling,
    TokenizerConfig,
    Word,
    tokenize_ipa,
)
from phonology.validation import evaluate_row  # noqa: E402

TIE = "\N{COMBINING DOUBLE INVERTED BREVE}"  # t͡ʃ
UNDERTIE = "\N{COMBINING DOUBLE BREVE BELOW}"  # t͜ʃ
NON_SYLLABIC = "\N{COMBINING INVERTED BREVE BELOW}"  # o̯, ɪ̯
RAISED = "\N{COMBINING UP TACK BELOW}"  # Czech r̝
SYLLABIC = "\N{COMBINING VERTICAL LINE BELOW}"  # Czech l̩
VOICELESS = "\N{COMBINING RING ABOVE}"  # r̝̊
DENTAL = "\N{COMBINING BRIDGE BELOW}"  # t̪
TILDE = "\N{COMBINING TILDE}"  # ɑ̃
IPA_G = "\N{LATIN SMALL LETTER SCRIPT G}"  # ɡ, not ASCII g
R_RAISED = "r" + RAISED  # Czech ř

# Czech-like convention: diacritics stay on their base, length on its
# segment, nothing folded.
CZECH_LIKE = TokenizerConfig(
    skip=NON_SEGMENTAL_MARKS,
    fold={},
    keep_combining=True,
    length=LengthHandling.ATTACH,
)

# Romanian lexicon IPA (lemma and plural columns), with the tokens the
# pre-config tokenize_ipa gave.
ROMANIAN_STRINGS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("ˈmus.kə", ("m", "u", "s", "k", "ə")),
    ("ˈmuʃ.te", ("m", "u", "ʃ", "t", "e")),
    ("prost", ("p", "r", "o", "s", "t")),
    ("ˈproʃ.tʲ", ("p", "r", "o", "ʃ", "t", "ʲ")),
    ("ˈpunkt", ("p", "u", "n", "k", "t")),
    ("ˈbro̯as.kə", ("b", "r", "o̯", "a", "s", "k", "ə")),
    ("de̯al", ("d", "e̯", "a", "l")),
    ("klɔɨ̯", ("k", "l", "o", "ɨ̯")),
    ("ˈt͡ʃe.ru", ("t͡ʃ", "e", "r", "u")),
    ("ˈa.d͡ʒi", ("a", "d͡ʒ", "i")),
    ("ˈfet͡sɛ", ("f", "e", "t͡s", "e")),
    ("bɔk", ("b", "o", "k")),
    ("greʃ", (IPA_G, "r", "e", "ʃ")),
    ("ˈkaːs", ("k", "a", "s")),
    ("ˌaɡa", ("a", IPA_G, "a")),
    ("vekʲ", ("v", "e", "k", "ʲ")),
    ("ˈruʃʲ", ("r", "u", "ʃ", "ʲ")),
    # Variants of a lexicon field (an undertie; a French nasal vowel).
    ("je" + UNDERTIE + "ofon", ("j", "e", UNDERTIE, "o", "f", "o", "n")),
    ("mə ʁɑ" + TILDE + "ʒ", ("m", "ə", "ʁ", "ɑ", TILDE, "ʒ")),
)

# Edge cases of the Romanian grouping, pinned as they always were.
ROMANIAN_EDGE_CASES: tuple[tuple[str, tuple[str, ...]], ...] = (
    # ɛ folds to e on its own, but not inside a ɛ̯ token.
    ("ɛ" + NON_SYLLABIC + "a", ("ɛ" + NON_SYLLABIC, "a")),
    # A tie bar binds exactly one following code point...
    ("t" + TIE, ("t", TIE)),
    # ...and a mark after the tied pair is a token of its own.
    ("t" + TIE + "ʃ" + NON_SYLLABIC, ("t" + TIE + "ʃ", NON_SYLLABIC)),
    ("t" + DENTAL + TIE + "s" + DENTAL, ("t", DENTAL + TIE + "s", DENTAL)),
    (R_RAISED, ("r", RAISED)),
    ("aːː", ("a",)),
    (LENGTH_MARK, ()),
    ("", ()),
    ("a ˈb\tc`d'e", ("a", "b", "c", "d", "e")),
)


# ---------------------------------------------------------------------------
# The Romanian default
# ---------------------------------------------------------------------------


class TestRomanianDefault:
    def test_romanian_tokenizer_settings(self) -> None:
        assert ROMANIAN_TOKENIZER.skip == NON_SEGMENTAL_MARKS
        assert NON_SEGMENTAL_MARKS == {" ", "\t", ".", "ˈ", "ˌ", "'", "`"}
        assert dict(ROMANIAN_TOKENIZER.fold) == {
            "ɛ": "e",
            "ɔ": "o",
            "g": IPA_G,
        }
        assert ROMANIAN_TOKENIZER.keep_combining is False
        assert ROMANIAN_TOKENIZER.length is LengthHandling.DROP

    @pytest.mark.parametrize(
        ("ipa", "expected"), ROMANIAN_STRINGS + ROMANIAN_EDGE_CASES
    )
    def test_tokens_are_unchanged(
        self, ipa: str, expected: tuple[str, ...]
    ) -> None:
        assert tokenize_ipa(ipa) == expected
        assert tokenize_ipa(ipa, tokenizer=ROMANIAN_TOKENIZER) == expected

    def test_profile_default_is_romanian(self) -> None:
        default = next(
            f.default
            for f in dataclasses.fields(AnalysisProfile)
            if f.name == "tokenizer"
        )
        assert default is ROMANIAN_TOKENIZER
        assert romanian.build_profile().tokenizer is ROMANIAN_TOKENIZER


# ---------------------------------------------------------------------------
# Other conventions
# ---------------------------------------------------------------------------


def _czech(ipa: str) -> tuple[str, ...]:
    return tokenize_ipa(ipa, tokenizer=CZECH_LIKE)


class TestOtherConventions:
    def test_keeps_length_lax_vowels_and_raised_r_distinct(self) -> None:
        ipa = "ˈ" + R_RAISED + "ɛːkɔ"
        assert _czech(ipa) == (R_RAISED, "ɛː", "k", "ɔ")
        # The Romanian default splits r̝, folds ɛ and ɔ, and drops ː.
        assert tokenize_ipa(ipa) == ("r", RAISED, "e", "k", "o")

    @pytest.mark.parametrize(
        ("tense", "lax"),
        [
            ("ˈpe.ska", "ˈpɛ.ska"),  # Italian pésca 'fishing', pèsca 'peach'
            ("ˈbot.te", "ˈbɔt.te"),  # Italian botte 'barrel', bòtte 'blows'
        ],
    )
    def test_lax_vowels_contrast_with_tense_ones(
        self, tense: str, lax: str
    ) -> None:
        kept = _czech(lax)
        assert "ɛ" in kept or "ɔ" in kept
        assert kept != _czech(tense)
        # The Romanian default folds each pair together.
        assert tokenize_ipa(lax) == tokenize_ipa(tense)

    @pytest.mark.parametrize(
        ("length", "expected"),
        [
            (LengthHandling.ATTACH, ("m", "aː", "m", "a")),
            (LengthHandling.SEPARATE, ("m", "a", LENGTH_MARK, "m", "a")),
            (LengthHandling.DROP, ("m", "a", "m", "a")),
        ],
    )
    def test_length_handling(
        self, length: LengthHandling, expected: tuple[str, ...]
    ) -> None:
        config = dataclasses.replace(CZECH_LIKE, length=length)
        assert tokenize_ipa("ˈmaː.ma", tokenizer=config) == expected

    def test_attached_length_closes_the_whole_segment(self) -> None:
        ts = "t" + TIE + "s"
        assert _czech(ts + "ːa") == (ts + "ː", "a")
        assert _czech(R_RAISED + "ː") == (R_RAISED + "ː",)
        assert _czech("aːː") == ("aːː",)
        # A mark that follows no segment is a token of its own.
        assert _czech("ːa") == ("ː", "a")
        assert _czech("a.ː") == ("a", "ː")

    def test_combining_marks_stay_on_their_base(self) -> None:
        # Czech vlk 'wolf' (syllabic l̩) and pout 'pilgrimage' (o + u̯).
        assert _czech("ˈvl" + SYLLABIC + "k") == ("v", "l" + SYLLABIC, "k")
        pout = ("p", "o", "u" + NON_SYLLABIC, "t")
        assert _czech("ˈpou" + NON_SYLLABIC + "t") == pout
        # Several marks on one base; a tie bar between marked bases.
        assert _czech(R_RAISED + VOICELESS) == (R_RAISED + VOICELESS,)
        ts = "t" + DENTAL + TIE + "s" + DENTAL
        assert _czech(ts + "a") == (ts, "a")

    def test_undertie_binds_like_the_tie_bar(self) -> None:
        tsh = "t" + UNDERTIE + "ʃ"
        assert _czech(tsh + "a") == (tsh, "a")
        # A tie bar with nothing after it stays on its base.
        assert _czech("t" + TIE) == ("t" + TIE,)

    def test_fold_maps_whole_tokens(self) -> None:
        tsh = "t" + TIE + "ʃ"
        config = dataclasses.replace(
            CZECH_LIKE,
            fold={
                "u" + NON_SYLLABIC: "w",
                R_RAISED + VOICELESS: R_RAISED,
                "t" + UNDERTIE + "ʃ": tsh,
            },
        )

        def tokens(ipa: str) -> tuple[str, ...]:
            return tokenize_ipa(ipa, tokenizer=config)

        assert tokens("pou" + NON_SYLLABIC + "t") == ("p", "o", "w", "t")
        assert tokens("tr" + RAISED + VOICELESS + "i") == ("t", R_RAISED, "i")
        assert tokens("t" + UNDERTIE + "ʃa") == (tsh, "a")

    def test_modifier_letters_are_tokens_of_their_own(self) -> None:
        for config in (ROMANIAN_TOKENIZER, CZECH_LIKE):
            assert tokenize_ipa("tʲa", tokenizer=config) == ("t", "ʲ", "a")
            assert tokenize_ipa("pʰa", tokenizer=config) == ("p", "ʰ", "a")

    def test_a_stray_mark_is_a_token_of_its_own(self) -> None:
        assert _czech(SYLLABIC + "a") == (SYLLABIC, "a")


# ---------------------------------------------------------------------------
# TokenizerConfig itself
# ---------------------------------------------------------------------------


class TestTokenizerConfig:
    def test_is_immutable(self) -> None:
        with pytest.raises(dataclasses.FrozenInstanceError):
            CZECH_LIKE.keep_combining = False  # type: ignore[misc]
        with pytest.raises(TypeError):
            ROMANIAN_TOKENIZER.fold["ɛ"] = "ɛ"  # type: ignore[index]

    def test_copies_and_coerces_its_fields(self) -> None:
        source = {"g": IPA_G}
        config = TokenizerConfig(
            skip={"."},  # type: ignore[arg-type]
            fold=source,
            keep_combining=False,
            length="ATTACH",  # type: ignore[arg-type]
        )
        source["x"] = "y"
        assert isinstance(config.skip, frozenset)
        assert dict(config.fold) == {"g": IPA_G}
        assert config.length is LengthHandling.ATTACH

    def test_equality_and_hash_are_by_value(self) -> None:
        rebuilt = TokenizerConfig(
            skip=frozenset(NON_SEGMENTAL_MARKS),
            fold={"ɛ": "e", "ɔ": "o", "g": IPA_G},
            keep_combining=False,
            length=LengthHandling.DROP,
        )
        assert rebuilt == ROMANIAN_TOKENIZER
        assert hash(rebuilt) == hash(ROMANIAN_TOKENIZER)
        assert len({rebuilt, ROMANIAN_TOKENIZER, CZECH_LIKE}) == 2
        assert dataclasses.replace(ROMANIAN_TOKENIZER) == ROMANIAN_TOKENIZER
        assert CZECH_LIKE != ROMANIAN_TOKENIZER

    @pytest.mark.parametrize(
        "override",
        [
            {"skip": NON_SEGMENTAL_MARKS | {LENGTH_MARK}},
            {"skip": frozenset({"ˈˌ"})},
            {"fold": {"": "e"}},
            {"fold": {"ɛ": ""}},
        ],
        ids=["length-in-skip", "multichar-skip", "empty-key", "empty-label"],
    )
    def test_rejects_ambiguous_settings(
        self, override: dict[str, Any]
    ) -> None:
        settings: dict[str, Any] = {
            "skip": NON_SEGMENTAL_MARKS,
            "fold": {},
            "keep_combining": True,
            "length": LengthHandling.ATTACH,
        }
        with pytest.raises(ValueError):
            TokenizerConfig(**{**settings, **override})


# ---------------------------------------------------------------------------
# The profile's tokenizer reaches the distance metric
# ---------------------------------------------------------------------------

MINI_FEATURES = ("Syllabic", "Consonantal", "Strident")
MINI_SEGMENTS = {
    "a": dict(zip(MINI_FEATURES, ("+", "-", "0"), strict=True)),
    "k": dict(zip(MINI_FEATURES, ("-", "+", "-"), strict=True)),
    "r": dict(zip(MINI_FEATURES, ("-", "+", "-"), strict=True)),
    R_RAISED: dict(zip(MINI_FEATURES, ("-", "+", "+"), strict=True)),
}

# Czech kař, written with ř = r̝: its last consonant is r̝, not r.
KAR = "ka" + R_RAISED


@pytest.fixture(scope="module")
def mini_json(tmp_path_factory: pytest.TempPathFactory) -> Path:
    path = tmp_path_factory.mktemp("mini") / "mini_features.json"
    inventory = {"features": list(MINI_FEATURES), "segments": MINI_SEGMENTS}
    path.write_text(
        json.dumps(inventory, ensure_ascii=False), encoding="utf-8"
    )
    return path


@pytest.fixture(scope="module")
def mini_inventory(mini_json: Path) -> FeatureInventory:
    return FeatureInventory.load(mini_json)


def _czech_ur(row: LexRow, inventory: FeatureInventory) -> Word | None:
    """Tokenize the lemma the Czech-like way; no suffix, no rules."""
    labels = tokenize_ipa(row.ipa_lemma, tokenizer=CZECH_LIKE)
    if not set(labels) <= set(inventory.base):
        return None
    return tuple(inventory.segment(label) for label in labels)


def _no_firings(row: LexRow) -> frozenset[str]:
    return frozenset()


def _never(row: LexRow) -> bool:
    return False


@pytest.fixture(scope="module")
def mini_profile(
    mini_json: Path, mini_inventory: FeatureInventory
) -> AnalysisProfile:
    """A profile left on the default (Romanian) tokenizer."""
    return AnalysisProfile(
        name="mini",
        inventory=mini_inventory,
        patches=(),
        underspec=(),
        rules=(),
        inventory_json=mini_json,
        target_stems=frozenset({"ř"}),
        palatalization_rule_names=frozenset(),
        expected_firings=_no_firings,
        fallback_predict=_never,
        ur_builder=_czech_ur,
    )


class TestProfileTokenizer:
    def test_feature_edit_distance_splits_with_the_given_tokenizer(
        self, mini_inventory: FeatureInventory
    ) -> None:
        predicted = mini_inventory.segment(R_RAISED)
        # Split as r + U+031D, the attested last consonant is plain /r/,
        # which differs from the predicted /r̝/ in Strident.
        assert feature_edit_distance(predicted, KAR, mini_inventory) == 1.0
        czech = feature_edit_distance(
            predicted, KAR, mini_inventory, tokenizer=CZECH_LIKE
        )
        assert czech == 0.0

    def test_validator_scores_with_the_profile_tokenizer(
        self, mini_profile: AnalysisProfile
    ) -> None:
        row = LexRow(
            lemma="kař",
            plural="kaře",
            pos="N",
            stem_final="ř",
            opportunity="e",
            mutation=False,
            ipa_lemma=KAR,
            ipa_pl=KAR,
            exception_reason="",
            cluster=ClusterTag.NONE,
        )
        default = evaluate_row(row, mini_profile)
        czech = evaluate_row(
            row, dataclasses.replace(mini_profile, tokenizer=CZECH_LIKE)
        )
        for res in (default, czech):
            assert res.predicted_sr == KAR
            assert res.ipa_match and res.boolean_match
        assert default.distance is not None
        assert default.distance.feature_edit == 1.0
        assert czech.distance is not None
        assert czech.distance.feature_edit == 0.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
