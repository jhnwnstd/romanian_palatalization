"""Golden tests for lexpipe.normalize (R8).

The strings are the ones the Romanian normalizers damage (see the
engineering critique in docs/multilingual_design.md's research basis):
Czech length and Italian final accents must survive, bracketed IPA
must parse, stress and syllable dots must be kept.
"""

from __future__ import annotations

import sys
import unicodedata
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lexpipe import normalize, spec  # noqa: E402
from lexpipe.normalize import Ipa, IpaKind  # noqa: E402

ITALIAN = normalize.normalizer_for("it")
CZECH = normalize.normalizer_for("cs")


# ---------------------------------------------------------------------------
# Orthography
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("word", ["žák", "státy", "čeští", "diplomaté"])
def test_czech_length_is_kept(word: str) -> None:
    assert CZECH.orthography(word) == word
    assert CZECH.match_key(word) == word
    assert CZECH.match_key(word.upper()) == word
    decomposed = unicodedata.normalize("NFD", word)
    assert CZECH.orthography(decomposed) == word
    assert CZECH.match_key("diplomaté") != CZECH.match_key("diplomate")


@pytest.mark.parametrize("word", ["città", "perché", "caffè americano"])
def test_italian_final_accents_are_kept(word: str) -> None:
    assert ITALIAN.orthography(word) == word
    assert ITALIAN.match_key(word) == word
    assert normalize.strip_respelling_accents(word) == word


def test_respelling_accents_are_stripped_for_match_keys_only() -> None:
    assert ITALIAN.match_key("fòcora") == "focora"
    assert ITALIAN.match_key("Àntico") == "antico"
    assert ITALIAN.match_key("nómora") == "nomora"
    assert ITALIAN.orthography("fòcora") == "fòcora"
    assert normalize.strip_respelling_accents("pòpò") == "popò"
    # Czech keys never strip anything.
    assert CZECH.match_key("fòcora") == "fòcora"


def test_orthography_tidies_whitespace_and_case() -> None:
    assert normalize.normalize_orthography("  caffè\t americano ") == (
        "caffè americano"
    )
    assert normalize.normalize_orthography("Città", lowercase=True) == "città"
    assert normalize.normalize_orthography("Città") == "Città"
    assert ITALIAN.match_key("po’") == "po'"


# ---------------------------------------------------------------------------
# IPA
# ---------------------------------------------------------------------------


def test_bracketed_ipa_parses_as_phonetic() -> None:
    ipa = CZECH.ipa("[ˈkluk]")
    assert ipa == Ipa("ˈkluk", IpaKind.PHONETIC)
    assert ipa.delimited == "[ˈkluk]"
    assert ITALIAN.ipa("/ˈfwɔ.ko/") == Ipa("ˈfwɔ.ko", IpaKind.PHONEMIC)


@pytest.mark.parametrize(
    "raw", ["", "ˈkluk", "*/tsio/", "/a/, /b/", "/ˈfwɔ.ko", "[]", "//"]
)
def test_anything_else_is_not_a_transcription(raw: str) -> None:
    assert normalize.parse_ipa(raw) is None


def test_ascii_colon_becomes_length() -> None:
    assert ITALIAN.ipa("/ˈɡrɛ:ko/") == Ipa("ˈɡrɛːko", IpaKind.PHONEMIC)


@pytest.mark.parametrize(
    ("raw", "italian", "czech"),
    [
        # adjacent halves: tied in Italian (it-pr syllabifies with dots)
        ("/ski.tso.freˈni.a/", "ski.t͡so.freˈni.a", "ski.tso.freˈni.a"),
        ("/ˈkatʃ.tʃa/", "ˈkat͡ʃ.t͡ʃa", "ˈkatʃ.tʃa"),
        # a syllable dot between the halves: never tied
        ("/ˈɡrat.sja/", "ˈɡrat.sja", "ˈɡrat.sja"),
        # Czech prefix-boundary clusters stay clusters
        ("[ˈpodzɪm]", "ˈpod͡zɪm", "ˈpodzɪm"),
        ("[ˈnatʃɛɲiː]", "ˈnat͡ʃɛɲiː", "ˈnatʃɛɲiː"),
        # already tied: unchanged, never doubled
        ("[ˈt͡sar]", "ˈt͡sar", "ˈt͡sar"),
        ("/ˈd͡zɛ.ro/", "ˈd͡zɛ.ro", "ˈd͡zɛ.ro"),
        # ligatures and the undertie are always written with the tie bar
        ("/ʧao/", "t͡ʃao", "t͡ʃao"),
        ("/ʤʦʣ/", "d͡ʒt͡sd͡z", "d͡ʒt͡sd͡z"),
        ("/at͜sa/", "at͡sa", "at͡sa"),
        # stress, secondary stress, dots, optional segments are kept
        ("/ˌpal.laˈvo.lo/", "ˌpal.laˈvo.lo", "ˌpal.laˈvo.lo"),
        ("/bri.koˈla(d)ʒ/", "bri.koˈla(d)ʒ", "bri.koˈla(d)ʒ"),
    ],
)
def test_affricate_tie_bars(raw: str, italian: str, czech: str) -> None:
    it_ipa, cs_ipa = ITALIAN.ipa(raw), CZECH.ipa(raw)
    assert it_ipa is not None and cs_ipa is not None
    assert it_ipa.text == italian
    assert cs_ipa.text == czech


def test_tie_affricates_is_idempotent() -> None:
    once = normalize.tie_affricates("tʃatsadʒadz")
    assert once == "t͡ʃat͡sad͡ʒad͡z"
    assert normalize.tie_affricates(once) == once


# ---------------------------------------------------------------------------
# Syllables and stress
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("raw", "syllables", "from_end", "post_tonic_final"),
    [
        ("/ˈfwɔ.ko/", ["ˈfwɔ", "ko"], 2, True),  # fuoco
        ("/ˈmɛ.di.ko/", ["ˈmɛ", "di", "ko"], 3, False),  # medico
        ("/aˈmi.ko/", ["a", "ˈmi", "ko"], 2, True),  # amico
        ("/ˈan.ti.ko/", ["ˈan", "ti", "ko"], 3, False),  # àntico
        ("/t͡ʃitˈta/", ["t͡ʃit", "ˈta"], 1, False),  # città
        ("/ˈe.vi.ta.no/", ["ˈe", "vi", "ta", "no"], 4, False),
        ("/ˌpal.laˈvo.lo/", ["ˌpal", "la", "ˈvo", "lo"], 2, True),
        ("[ˈkluk]", ["ˈkluk"], 1, False),
    ],
)
def test_stress(
    raw: str,
    syllables: list[str],
    from_end: int,
    post_tonic_final: bool,
) -> None:
    assert normalize.syllabify(raw) == syllables
    assert normalize.is_syllabified(raw)
    assert normalize.stress_from_end(raw) == from_end
    assert normalize.post_tonic(raw, -1) is post_tonic_final
    ipa = ITALIAN.ipa(raw)
    assert ipa is not None and normalize.stress_from_end(ipa) == from_end


@pytest.mark.parametrize(
    "raw",
    [
        "[ˈabort]",  # Czech: no syllable dots, two nuclei
        "/feudaˈleːzimo/",  # Italian without dots
        "[ˈdiːzl̩]",  # syllabic l is a second nucleus
        "/a/",  # no stress mark
        "/ˈdi ˈvu ˈdi/",  # several words
    ],
)
def test_no_stress_without_syllabified_single_words(raw: str) -> None:
    assert normalize.stress_from_end(raw) is None
    assert normalize.post_tonic(raw, -1) is None


def test_nuclei() -> None:
    assert normalize.nuclei("ˈfwɔ") == 1  # glide + vowel
    assert normalize.nuclei("sou̯") == 1  # non-syllabic u
    assert normalize.nuclei("vl̩k") == 1  # syllabic consonant
    assert normalize.nuclei("ˈabort") == 2
    assert normalize.nuclei("ˈɡrɛːko") == 2
    assert normalize.is_syllabified("[ˈsou̯.sɛt]")


def test_post_tonic_indexes() -> None:
    assert normalize.post_tonic("/ˈmɛ.di.ko/", 1) is True
    assert normalize.post_tonic("/ˈmɛ.di.ko/", -2) is True
    assert normalize.post_tonic("/ˈmɛ.di.ko/", 0) is False
    with pytest.raises(IndexError):
        normalize.post_tonic("/ˈmɛ.di.ko/", 3)
    with pytest.raises(IndexError):
        normalize.post_tonic("/ˈmɛ.di.ko/", -4)


def test_stress_class_names() -> None:
    assert normalize.STRESS_CLASS_NAMES[2] == "paroxytone"
    assert normalize.STRESS_CLASS_NAMES[3] == "proparoxytone"


def test_one_normalizer_per_registered_language() -> None:
    # Comparison-only rows (Romanian) build no lexicon: no normalizer.
    assert set(normalize.NORMALIZERS) == set(spec.full_languages())
    assert spec.LANGUAGES["ro"].comparison_only
    with pytest.raises(KeyError, match="known"):
        normalize.normalizer_for("ro")
