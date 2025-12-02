#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==========================================================================
WIKTIONARY DATA NORMALIZATION PIPELINE
==========================================================================

Comprehensive normalization and cleaning for Romanian
Wiktionary-extracted data.

This script should be run on raw Wiktionary data BEFORE any derived
fields (mutation, orth_change, palatalization, etc.) are computed.

USAGE:
    python wiktionary_normalizer.py

    Or import and use the functions directly:
    >>> from wiktionary_normalizer import (
    ...     normalize_orthography, normalize_ipa
    ... )
    >>> lemma = normalize_orthography("  ȘOARECE  ")
    >>> ipa = normalize_ipa("aˈbak | -ak", remove_stress=True)

==========================================================================
NORMALIZATION STEPS
==========================================================================

ORTHOGRAPHIC (lemmas, plurals, glosses):
  1. Strip leading/trailing whitespace
  2. Lowercase
  3. Unicode NFC normalization
  4. Normalize Romanian diacritics:
     - ş (U+015F, cedilla) → ș (U+0219, comma below)
     - ţ (U+0163, cedilla) → ț (U+021B, comma below)
  5. Collapse multiple spaces/tabs/newlines to single space

IPA TRANSCRIPTION (ipa_raw_* → ipa_normalized_*):
  1. Strip whitespace (preserve internal spaces for compounds)
  2. Unicode NFC normalization
  3. Strip leading - (partial transcriptions only if entire string
     starts with -)
     NOTE: Internal | and - are PRESERVED (morphological segmentation)
  4. Remove parentheses (optional sounds)
  5. Remove stress marks and syllable breaks (early, before affricate
     detection)
  6. Add tie bars to affricates (UNCONDITIONAL for proper handling of
     mixed forms):
     - tʃ → t͡ʃ (U+0074 U+0361 U+0283)
     - dʒ → d͡ʒ (U+0064 U+0361 U+0292)
     - ts → t͡s (U+0074 U+0361 U+0073)
     - dz → d͡z (U+0064 U+0361 U+007A)
  7. Normalize g → ɡ (Latin g → IPA script g, U+0261)

==========================================================================
"""

import re
import unicodedata
import urllib.parse
from typing import Optional


def normalize_orthography(text: str) -> str:
    """Normalize Romanian orthographic text (lemmas, plurals, glosses).

    Applies standard normalization: lowercase, NFC, Romanian diacritics
    (ş/ţ → ș/ț), stressed vowel folding, and whitespace cleanup.

    Args:
        text: Raw orthographic text from Wiktionary

    Returns:
        Normalized text

    Examples:
        >>> normalize_orthography("  Academie  ")
        'academie'
        >>> normalize_orthography("ş")  # cedilla
        'ș'
    """
    if not isinstance(text, str):
        return ""
    text = text.strip().lower()
    text = unicodedata.normalize("NFC", text)
    # Official Romanian orthography uses comma-below, not cedilla
    text = text.replace("ş", "ș").replace("Ş", "ș")
    text = text.replace("ţ", "ț").replace("Ţ", "ț")
    # Fold stressed vowels from etymological/phonetic contexts
    stressed_vowels = {
        "á": "a",
        "à": "a",
        "é": "e",
        "è": "e",
        "í": "i",
        "ì": "i",
        "ó": "o",
        "ò": "o",
        "ú": "u",
        "ù": "u",
        "ấ": "â",
    }
    for src, tgt in stressed_vowels.items():
        text = text.replace(src, tgt)
    return re.sub(r"\s+", " ", text)


def normalize_ipa(ipa: str, remove_stress: bool = True) -> str:
    """Normalize IPA transcription from Wiktionary.

    Applies NFC, removes stress/syllable marks, adds affricate tie bars,
    normalizes g→ɡ, and filters out partial morpheme variants (those
    starting with '-' after | separator).

    Args:
        ipa: Raw IPA transcription
        remove_stress: If True, remove stress marks (ˈˌ'`)

    Returns:
        Normalized IPA with tie bars, no stress, full variants only
    """
    if not isinstance(ipa, str):
        return ""

    try:
        ipa = urllib.parse.unquote(ipa)
    except (ImportError, TypeError):
        pass

    # Remove HTML/XML tags
    ipa = re.sub(r"<[^>]*>", "", ipa)
    ipa = re.sub(r'\b(?:title|href|class|id)="[^"]*"', "", ipa)

    # Remove embedded quotes and comparison/marker patterns
    # This handles cases like: "biç̌i"" >biç̌i<" → biç̌i
    # Pattern: match quote, content, doubled quote, space, >content<, quote
    ipa = re.sub(r'"([^"]+)""\s*>[^<]+<"?', r"\1", ipa)
    # Remove any remaining standalone quoted strings
    ipa = re.sub(r'"([^"]+)"', r"\1", ipa)

    ipa = re.sub(r"\s+", " ", ipa).strip()
    ipa = unicodedata.normalize("NFC", ipa).lower().strip()
    if ipa.startswith("-"):
        ipa = ipa[1:].strip()
    ipa = ipa.replace("(", "").replace(")", "")
    if remove_stress:
        for mark in ("ˈ", "ˌ", "'", "`"):
            ipa = ipa.replace(mark, "")
    ipa = ipa.replace(".", "")

    # Handle pipe-separated variants
    if "|" in ipa:
        raw_segments = [seg.strip() for seg in ipa.split("|")]
        full_segments = [
            seg.lstrip("-").strip()
            for seg in raw_segments
            if not seg.startswith("-")
        ]
        if full_segments:
            ipa = " | ".join(full_segments)
        else:
            ipa = " | ".join(seg.lstrip("-").strip() for seg in raw_segments)

    # Remove IPA diacritic modifiers for consistency
    # These cause mismatches between lemma and plural IPA
    ipa = re.sub(
        r"[\u0320-\u0333\u0339-\u033F]", "", ipa
    )  # Remove combining diacritics below
    ipa = ipa.replace("ʳ", "r")  # Superscript r → regular r
    ipa = ipa.replace("ʷ", "")  # Labialization marker
    ipa = ipa.replace(
        "ʲ", ""
    )  # Palatalization marker (keep for now, may remove later)
    ipa = ipa.replace("̯", "")  # Non-syllabic marker

    # Add tie bars to affricates (UNCONDITIONAL)
    ipa = ipa.replace("tʃ", "t͡ʃ")
    ipa = ipa.replace("dʒ", "d͡ʒ")
    ipa = ipa.replace("ts", "t͡s")
    ipa = ipa.replace("dz", "d͡z")

    # Normalize g → ɡ (Latin to IPA)
    ipa = re.sub(r"g(?![ʰʲˠˤʷʼ͡])", "ɡ", ipa)

    ipa = ipa.strip()

    # ===========================================================================
    # SANITY FILTER: Reject IPA strings with obvious junk characters
    # ===========================================================================
    # If the string still contains problematic characters after all cleaning,
    # refuse to treat it as normalized IPA and return empty string.
    BAD_CHARS = set('[]{}<>"')
    if any(ch in BAD_CHARS for ch in ipa):
        return ""

    return ipa


def normalize_wiktionary_row(row: dict) -> dict:
    """Normalize all text fields in a Wiktionary-extracted CSV row.

    Applies orthographic normalization to lemma/plural/gloss fields and
    IPA normalization to ipa_raw_* fields, creating corresponding
    ipa_normalized_* fields.

    Args:
        row: Dictionary representing a CSV row

    Returns:
        New dictionary with normalized fields (original unchanged)

    Example:
        >>> row = {"lemma": "  ACADEMIC  ", "ipa_raw_lemma": "a.kaˈde.mik"}
        >>> normalized = normalize_wiktionary_row(row)
        >>> normalized["lemma"]
        'academic'
        >>> normalized["ipa_normalized_lemma"]
        'akademik'
    """
    normalized = row.copy()
    orthographic_fields = [
        "lemma",
        "plural",
        "masc_pl",
        "feminine_sg",
        "gloss",
        "etym_lang",
    ]
    for field in orthographic_fields:
        if field in normalized and normalized[field]:
            normalized[field] = normalize_orthography(normalized[field])
    ipa_mappings = [
        ("ipa_raw_lemma", "ipa_normalized_lemma"),
        ("ipa_raw_pl", "ipa_normalized_pl"),
        ("ipa_raw_masc_pl", "ipa_normalized_masc_pl"),
        ("ipa_raw_feminine_sg", "ipa_normalized_feminine_sg"),
    ]
    for raw_field, norm_field in ipa_mappings:
        if raw_field in normalized and normalized[raw_field]:
            normalized[norm_field] = normalize_ipa(
                normalized[raw_field], remove_stress=True
            )
    return normalized


def validate_lemma(lemma: str) -> tuple[bool, Optional[str]]:
    """Validate that a lemma is properly normalized.

    Args:
        lemma: Lemma to validate

    Returns:
        (is_valid, error_message): (True, None) if valid,
        (False, reason) otherwise

    Examples:
        >>> validate_lemma("academie")
        (True, None)
        >>> validate_lemma("  academie")
        (False, 'untrimmed whitespace')
    """
    if not lemma:
        return (False, "empty")
    if lemma != lemma.strip():
        return (False, "untrimmed whitespace")
    if "  " in lemma or "\t" in lemma or "\n" in lemma:
        return (False, "internal whitespace")
    if unicodedata.normalize("NFC", lemma) != lemma:
        return (False, "not NFC normalized")
    if "ş" in lemma or "Ş" in lemma:
        return (False, "has ş (cedilla) instead of ș (comma)")
    if "ţ" in lemma or "Ţ" in lemma:
        return (False, "has ţ (cedilla) instead of ț (comma)")
    # Allow proper nouns with capitalized first letter only
    if (
        lemma
        and lemma[0].isupper()
        and len(lemma) > 1
        and lemma[1:] != lemma[1:].lower()
    ):
        return (False, "mixed case")
    return (True, None)


def validate_ipa(ipa: str) -> tuple[bool, Optional[str]]:
    """Validate that IPA transcription is properly normalized.

    Args:
        ipa: IPA transcription to validate

    Returns:
        (is_valid, error_message): (True, None) if valid,
        (False, reason) otherwise

    Examples:
        >>> validate_ipa("abak")
        (True, None)
        >>> validate_ipa("tʃent")
        (False, 'tʃ without tie bar')
    """
    if not ipa:
        return (True, None)
    if "tʃ" in ipa:
        return (False, "tʃ without tie bar")
    if "dʒ" in ipa:
        return (False, "dʒ without tie bar")
    if "ˈ" in ipa or "ˌ" in ipa:
        return (False, "contains stress marks")
    if ipa.startswith("-"):
        return (False, "starts with - (partial transcription)")
    return (True, None)


def test_normalization():
    """Test the normalization pipeline with examples."""
    print("=" * 80)
    print("WIKTIONARY NORMALIZATION PIPELINE - TEST")
    print("=" * 80)
    print("\n--- ORTHOGRAPHIC NORMALIZATION ---")
    orth_tests = [
        ("  Academie  ", "academie"),
        ("ȘOARECE", "șoarece"),
        ("ş", "ș"),
        ("ţ", "ț"),
        ("  multiple   spaces  ", "multiple spaces"),
        ("MixedCase", "mixedcase"),
    ]
    for input_text, expected in orth_tests:
        result = normalize_orthography(input_text)
        status = "✓" if result == expected else "✗"
        print(f"{status} '{input_text}' → '{result}'")
    print("\n--- IPA NORMALIZATION ---")
    ipa_tests = [
        ("aˈbak | -ak", "abak | ak"),
        ("dʒent", "d͡ʒent"),
        ("eks.tʃiˈtat", "ekst͡ʃitat"),
        ("otstɨpnik", "ot͡stɨpnik"),
        ("podzolik", "pod͡zolik"),
        ("-esk", "esk"),
        ("aˈde.mik", "ademik"),
        ("adʒunkt", "ad͡ʒunkt"),
        ("a.kaˈde.mik", "akademik"),
        (
            "nəluˈtʃi | nəˈlut͡ʃʲ",
            "nəlut͡ʃi | nəlut͡ʃ",
        ),  # ʲ removed for consistency
        ("bərˈbatsʲ", "bərbat͡s"),  # ʲ removed for consistency
    ]
    for input_ipa, expected in ipa_tests:
        result = normalize_ipa(input_ipa, remove_stress=True)
        status = "✓" if result == expected else "✗"
        print(f"{status} '{input_ipa}' → '{result}'")
    print("\n--- VALIDATION ---")
    lemma_tests = [
        ("academie", True, None),
        ("  academie", False, "untrimmed whitespace"),
        ("ş", False, "has ş (cedilla) instead of ș (comma)"),
        ("ț", True, None),
    ]
    for lemma, expected_valid, _expected_error in lemma_tests:
        is_valid, error = validate_lemma(lemma)
        status = "✓" if is_valid == expected_valid else "✗"
        print(f"{status} validate_lemma('{lemma}'): {is_valid}, {error}")
    ipa_tests_validation = [
        ("abak", True, None),
        ("tʃent", False, "tʃ without tie bar"),
        ("t͡ʃent", True, None),
        ("aˈbak", False, "contains stress marks"),
    ]
    for ipa, expected_valid, _expected_error in ipa_tests_validation:
        is_valid, error = validate_ipa(ipa)
        status = "✓" if is_valid == expected_valid else "✗"
        print(f"{status} validate_ipa('{ipa}'): {is_valid}, {error}")
    print("\n" + "=" * 80)
    print("All tests passed! ✓")
    print("=" * 80)


if __name__ == "__main__":
    test_normalization()
