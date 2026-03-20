#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Romanian Palatalization Data Processor - Complete Function Library

Consolidated module containing all derivation functions for Romanian
palatalization analysis.

Sections:
1. Part 1: Basic derivations (stem_final, cluster)
2. Part 2: Alignment and mutation
3. Part 3: Suffix fields and G2P
4. Part 4: NDE and exception handling
"""

import re
import unicodedata
from difflib import SequenceMatcher
from typing import Callable, Dict, List, Set, Tuple

_ipa_normalizer: Callable[[str], str] | None = None


def set_ipa_normalizer(normalizer: Callable[[str], str]) -> None:
    """Set the IPA normalization function to use for G2P outputs.

    The normalizer should be a function that takes an IPA string and returns
    a normalized IPA string (with stress removed, tie bars added, etc.).

    Args:
        normalizer: Function that takes (ipa: str) -> str
    """
    global _ipa_normalizer
    _ipa_normalizer = normalizer


# ============================================================================
# CONSTANTS
# ============================================================================

TARGET_CONSONANTS = {"c", "g", "t", "d", "s", "z"}

NDE_CLASSES = {"gimpe", "ochi", "paduchi"}

PLURAL_OPPORTUNITIES = {"i", "e"}

# Palatal surface forms for coronals (used in NDE alternation detection)
PALATAL_SURFACE = {
    "t": "ț",  # t → ț before i/e
    "d": "z",  # d → z before i/e (NOT when z is underlyingly in lemma)
    "s": "ș",  # s → ș before i/e
    "z": "j",  # z → j before i/e
}

VOWEL_SET = {
    "a",
    "ă",
    "â",
    "e",
    "i",
    "î",
    "o",
    "u",
}

FINAL_CLUSTERS = {
    "st": "s",
    "sc": "s",
    "ct": "t",
}

VELAR_FRONT_SEQUENCES = {
    "chi": "c",
    "che": "c",
    "ghi": "g",
    "ghe": "g",
}

ROMANIAN_CONSONANTS: Set[str] = set("bcdfghjklmnpqrstvwxyzșşţțțţ")

ORTH_TO_PALATAL_IPA = {
    # Core patterns for classification via suffix matching
    # Example: "ate→ăți" recognized via "te→ți" suffix match
    "c→ce": "t͡ʃ",
    "c→ci": "t͡ʃ",
    "ct→cți": "t͡s",
    "ct→cțe": "t͡s",
    "că→ce": "t͡ʃ",
    "că→ci": "t͡ʃ",
    "ca→ce": "t͡ʃ",
    "d→ze": "z",
    "d→zi": "z",
    "d→z": "z",  # bare consonant change
    "de→zi": "z",
    "dă→zi": "z",
    "g→ge": "d͡ʒ",
    "g→gi": "d͡ʒ",
    "gă→ge": "d͡ʒ",
    "gă→gi": "d͡ʒ",
    "ga→ge": "d͡ʒ",
    "go→gi": "d͡ʒ",
    "s→șe": "ʃ",
    "s→și": "ʃ",
    "s→ș": "ʃ",  # bare consonant change
    "sc→ște": "ʃt",
    "sc→ști": "ʃt",
    "scă→ște": "ʃt",
    "scă→ști": "ʃt",
    "st→ște": "ʃt",
    "st→ști": "ʃt",
    "t→țe": "t͡s",
    "t→ți": "t͡s",
    "t→ț": "t͡s",  # bare consonant change (e.g., componente → componențe)
    "te→țe": "t͡s",
    "te→ți": "t͡s",
    "tă→țe": "t͡s",
    "tă→ți": "t͡s",
    "z→je": "ʒ",
    "z→ji": "ʒ",
    "z→j": "ʒ",  # bare consonant change
}

# ============================================================================
# ROW FILTERING
# ============================================================================


def should_process_row(row: Dict[str, str]) -> bool:
    """
    Determine if a row should be processed based on whether it contains
    at least one target consonant in the lemma.

    This filters out rows that cannot possibly undergo palatalization,
    improving processing efficiency.

    Args:
        row: Dictionary representing a CSV row

    Returns:
        True if the row should be processed, False otherwise
    """
    lemma = row.get("lemma", "") or ""
    if not lemma:
        return False
    lemma_lower = lemma.lower()
    for consonant in TARGET_CONSONANTS:
        if consonant in lemma_lower:
            return True
    return False


def orth_change_is_vowels_or_le(orth_change: str) -> bool:
    """Return True if orth_change is vowel-only or includes 'le'.

    Examples: 'e→i', 'a→le', 'ă→ăe'

    This detects cases where the orthographic change is purely vocalic
    and should NOT be treated as introducing a new C+i/e opportunity
    for coronal stems.

    Used to filter out false opportunities like:
    - abazie → abazii (e→i): already has z+i in lemma
    - acadea → acadele (a→le): just vowel/suffix change
    - antiartă → antiartăe (ă→ăe): purely vocalic
    """
    if not orth_change or "→" not in orth_change:
        return False

    for ch in orth_change:
        if ch in {" ", "\t"}:
            continue
        if ch in {"→", "∅"}:
            continue
        if ch in VOWEL_SET:
            continue
        if ch.lower() == "l":
            continue
        return False

    return True


def ensure_ipa_fields(
    row: Dict[str, str],
    orth_key: str,
    raw_key: str,
    norm_key: str,
    tweak_fn: Callable[[str, str], str] | None = None,
) -> None:
    """
    Ensure IPA field is populated in the row, either from raw IPA or via G2P.

    This function modifies the row in place, adding the normalized IPA field.
    It first checks if a raw IPA annotation exists; if so, it normalizes it.
    Otherwise, it generates IPA via grapheme-to-phoneme conversion and
    optionally applies a tweak function.

    Args:
        row: Dictionary representing a CSV row (modified in place)
        orth_key: Key for the orthographic form (e.g., "lemma")
        raw_key: Key for raw IPA annotation (e.g., "ipa_raw_lemma")
        norm_key: Key to store normalized IPA (e.g., "ipa_normalized_lemma")
        tweak_fn: Optional function to adjust G2P output, (orth, ipa) -> ipa

    Requires:
        _ipa_normalizer must be set via set_ipa_normalizer() before calling
    """
    if _ipa_normalizer is None:
        raise RuntimeError(
            "IPA normalizer not set. Call set_ipa_normalizer() first."
        )

    raw_val = row.get(raw_key)
    if raw_val:
        row[norm_key] = _ipa_normalizer(raw_val)
        return

    orth_val = row.get(orth_key, "")
    if orth_val:
        ipa = to_ipa(orth_val)
        if tweak_fn is not None:
            ipa = tweak_fn(orth_val, ipa)
        row[norm_key] = _ipa_normalizer(ipa)


def strip_final_vowel(lemma: str) -> str:
    """Strip a single final vowel from lemma if present."""
    if len(lemma) <= 1:
        return lemma
    if lemma[-1] in VOWEL_SET:
        return lemma[:-1]
    return lemma


def longest_common_substring(a: str, b: str) -> int:
    """Return the length of the longest contiguous matching substring."""
    m = [[0] * (1 + len(b)) for _ in range(1 + len(a))]
    longest = 0
    for i, ca in enumerate(a, 1):
        for j, cb in enumerate(b, 1):
            if ca == cb:
                m[i][j] = m[i - 1][j - 1] + 1
                longest = max(longest, m[i][j])
    return longest


def jaccard_similarity(a: str, b: str) -> float:
    """Compute Jaccard similarity of character bigrams."""
    bigrams_a = {a[i : i + 2] for i in range(len(a) - 1)}
    bigrams_b = {b[i : i + 2] for i in range(len(b) - 1)}
    return len(bigrams_a & bigrams_b) / (len(bigrams_a | bigrams_b) or 1)


def common_prefix_length(a: str, b: str) -> int:
    """Count matching characters at start of strings."""
    for i, (ca, cb) in enumerate(zip(a, b)):
        if ca != cb:
            return i
    return min(len(a), len(b))


def common_suffix_length(a: str, b: str) -> int:
    """Count matching characters at end of strings."""
    for i, (ca, cb) in enumerate(zip(reversed(a), reversed(b))):
        if ca != cb:
            return i
    return min(len(a), len(b))


_PL_PLAUS_SCORES: List[float] = []
_PL_CALIBRATED = False
_PL_THRESHOLDS = {"reject": 0.35, "border": 0.45}
_PL_MIN_CALIBRATION = 500
# Cap calibrated reject threshold.  Romanian has legitimate plural
# patterns that score low on string similarity (short stem + -uri,
# stem-changing plurals like frate→frați, noapte→nopți).  Allowing
# calibration to push the threshold above this cap would discard them.
_PL_MAX_REJECT = 0.40


def _calibrate_plural_thresholds() -> None:
    """Update global plausibility thresholds from observed scores."""
    global _PL_CALIBRATED, _PL_THRESHOLDS  # noqa: F824
    n = len(_PL_PLAUS_SCORES)
    if n < _PL_MIN_CALIBRATION:
        return
    scores = sorted(_PL_PLAUS_SCORES)
    idx_reject = max(0, int(0.05 * (n - 1)))
    idx_border = max(0, int(0.20 * (n - 1)))
    reject_thr = min(scores[idx_reject], _PL_MAX_REJECT)
    border_thr = scores[idx_border]
    if border_thr <= reject_thr:
        border_thr = min(1.0, reject_thr + 0.05)
    _PL_THRESHOLDS["reject"] = reject_thr
    _PL_THRESHOLDS["border"] = border_thr
    _PL_CALIBRATED = True


def compute_plural_plausibility(
    lemma: str, plural: str
) -> Tuple[float, float, float]:
    """
    Compute multi-feature plausibility score for lemma-plural pair.
    Returns (plausibility, seq_sim, lcs_ratio).
    """
    if not lemma or not plural:
        return 0.0, 0.0, 0.0
    max_len = max(len(lemma), len(plural))
    seq_sim = SequenceMatcher(None, lemma, plural).ratio()
    lcs_len = longest_common_substring(lemma, plural)
    lcs_ratio = lcs_len / max_len if max_len > 0 else 0.0
    ngram_sim = jaccard_similarity(lemma, plural)
    edge_len = max(
        common_prefix_length(lemma, plural),
        common_suffix_length(lemma, plural),
    )
    edge_focus = edge_len / max_len if max_len > 0 else 0.0
    plausibility = (seq_sim + lcs_ratio + ngram_sim + edge_focus) / 4
    return plausibility, seq_sim, lcs_ratio


def validate_plural_quality(row: Dict[str, str]) -> None:
    """
    Unsupervised plausibility filter for lemma-plural pairs using
    edit distance, LCS, n-gram similarity, and edge matching.
    """
    pos = (row.get("pos", "") or "").upper()
    lemma = (row.get("lemma", "") or "").strip().lower()
    plural = (row.get("plural", "") or "").strip().lower()
    if pos not in {"N", "ADJ"}:
        return
    if not lemma or not plural:
        return
    if any(ch.isspace() for ch in plural) or (
        "-" in plural and not plural.endswith("uri")
    ):
        row["plural"] = ""
        row["ipa_normalized_pl"] = ""
        return
    len_ratio = len(plural) / len(lemma) if len(lemma) else 0.0
    len_ok = 0.5 <= len_ratio <= 2.0
    plausibility, seq_sim, lcs_ratio = compute_plural_plausibility(
        lemma, plural
    )

    if len_ok:
        _PL_PLAUS_SCORES.append(plausibility)
        if not _PL_CALIBRATED and len(_PL_PLAUS_SCORES) >= _PL_MIN_CALIBRATION:
            _calibrate_plural_thresholds()

    reject_thr = _PL_THRESHOLDS["reject"]
    border_thr = _PL_THRESHOLDS["border"]

    # Stem-change safeguard: Romanian has plurals with major vowel
    # alternation (noapte→nopți, cale→căi, fată→fete, tânăr→tineri)
    # that score very low on string similarity but are valid.
    # If the plural shares a consonant onset with the lemma and
    # has a reasonable length ratio, accept it at a lower threshold.
    stem_change_rescue = False
    if len_ok and plausibility >= 0.15 and len(lemma) >= 3:
        pfx = common_prefix_length(lemma, plural)
        if pfx >= 1 and not lemma[:pfx].replace("-", "").isspace():
            # Shares at least 1 leading character — likely a stem change
            stem_change_rescue = True

    if (
        plausibility < reject_thr
        or not len_ok
        or (seq_sim < 0.1 and lcs_ratio < 0.1)
    ):
        if not stem_change_rescue:
            row["plural"] = ""
            row["ipa_normalized_pl"] = ""
            return

    if _PL_CALIBRATED:
        row["plural_validity"] = (
            "borderline" if plausibility < border_thr else "ok"
        )


def derive_stem_final_and_cluster(row: Dict[str, str]) -> None:
    """
    Identify the target consonant for palatalization and any
    associated cluster (chi/che/ghi/ghe/st/sc/ct).
    """
    lemma = row.get("lemma", "") or ""
    if not lemma:
        row["stem_final"] = ""
        row["cluster"] = ""
        return

    lemma_l = lemma.lower()

    # Velar front clusters at word edge
    for seq, consonant in VELAR_FRONT_SEQUENCES.items():
        if lemma_l.endswith(seq):
            row["stem_final"] = consonant
            row["cluster"] = lemma[-len(seq) :]
            return

    # Final clusters after stripping final vowel
    stem = strip_final_vowel(lemma_l)
    for cluster, consonant in FINAL_CLUSTERS.items():
        if stem.endswith(cluster):
            row["stem_final"] = consonant
            row["cluster"] = lemma[
                -len(cluster) - (len(lemma_l) - len(stem)) : len(lemma)
            ]
            return

    # Bare single consonant at right edge (no consonant to its right)
    for i in range(len(stem) - 1, -1, -1):
        ch = stem[i]
        if ch in TARGET_CONSONANTS:
            trailing = stem[i + 1 :]
            if any(c in ROMANIAN_CONSONANTS for c in trailing):
                continue
            row["stem_final"] = ch
            row["cluster"] = ""
            return

    row["stem_final"] = ""
    row["cluster"] = ""


def derive_frontstem_dorsal(row: Dict[str, str]) -> None:
    """Mark lemmas with frontstem dorsals (ci/ce/gi/ge without h).

    These represent already-palatalized stems where the "target consonant"
    is actually the palatal output, not a plain dorsal.  Only matches
    ci/ce/gi/ge near the stem-final position (last 4 chars of the stem
    after stripping the inflectional vowel), not anywhere in the word.

    Examples: alice, indice, abazie (with 'zie' ending)

    This flag is used post-NDE to remove these from the alternation domain
    unless they're explicitly used as NDEB items.
    """
    lemma = (row.get("lemma", "") or "").strip().lower()
    stem_final = (row.get("stem_final", "") or "").strip().lower()

    row["frontstem_dorsal"] = "False"

    if not lemma or stem_final not in {"c", "g"}:
        return

    # Only check for ci/ce/gi/ge near the stem-final position.
    # Check both the lemma tail (to catch ci/ce at word boundary, e.g.
    # bici, ghiveci) and the stem tail after stripping the inflectional
    # vowel (to catch medial sequences like -gic-, -cie-).
    stem = strip_final_vowel(lemma)
    lemma_tail = lemma[-4:] if len(lemma) >= 4 else lemma
    stem_tail = stem[-4:] if len(stem) >= 4 else stem
    if re.search(r"(?<!h)[cg][ie]", lemma_tail) or re.search(
        r"(?<!h)[cg][ie]", stem_tail
    ):
        row["frontstem_dorsal"] = "True"


def needleman_wunsch(s1: str, s2: str) -> Tuple[str, str]:
    """
    Global alignment using Needleman-Wunsch algorithm.
    Returns (aligned_s1, aligned_s2) with '-' for gaps.
    """
    m, n = len(s1), len(s2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        dp[i][0] = -i
    for j in range(n + 1):
        dp[0][j] = -j
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = dp[i - 1][j - 1] + (1 if s1[i - 1] == s2[j - 1] else -1)
            delete = dp[i - 1][j] - 1
            insert = dp[i][j - 1] - 1
            dp[i][j] = max(match, delete, insert)
    aligned_s1, aligned_s2 = [], []
    i, j = m, n
    while i > 0 or j > 0:
        if (
            i > 0
            and j > 0
            and dp[i][j]
            == dp[i - 1][j - 1] + (1 if s1[i - 1] == s2[j - 1] else -1)
        ):
            aligned_s1.append(s1[i - 1])
            aligned_s2.append(s2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i - 1][j] - 1:
            aligned_s1.append(s1[i - 1])
            aligned_s2.append("-")
            i -= 1
        else:
            aligned_s1.append("-")
            aligned_s2.append(s2[j - 1])
            j -= 1
    return "".join(reversed(aligned_s1)), "".join(reversed(aligned_s2))


def has_tautomorphemic_front_vowel(
    lemma: str, plural: str, stem_final: str
) -> bool:
    """
    Return True iff, around the stem_final consonant, the lemma and plural
    both have the SAME front vowel (i/e) immediately following it.

    Intuition:
      - Align lemma/plural.
      - Locate the aligned position of the *target* consonant in lemma.
      - Look for the next real segment to the right in both lemma and plural.
      - If both share 'i' or 'e' there, that C+front-vowel sequence is
        tautomorphemic (not introduced by a suffix) ⇒ GIMPE-like.

    Used for *all* TARGET_CONSONANTS (c, g, t, d, s, z).
    """
    if not lemma or not plural or not stem_final:
        return False

    lemma = lemma.lower()
    plural = plural.lower()
    stem_final = stem_final.lower()

    if stem_final not in TARGET_CONSONANTS:
        return False

    # Global alignment
    aligned_lemma, aligned_plural = needleman_wunsch(lemma, plural)

    # Map lemma indices → aligned positions
    lemma_to_aligned: Dict[int, int] = {}
    li = 0
    for ai, ch in enumerate(aligned_lemma):
        if ch != "-":
            lemma_to_aligned[li] = ai
            li += 1

    # Take the rightmost occurrence of stem_final in the lemma
    try:
        target_idx = max(i for i, ch in enumerate(lemma) if ch == stem_final)
    except ValueError:
        return False

    if target_idx not in lemma_to_aligned:
        return False

    aidx = lemma_to_aligned[target_idx]

    # Walk right in the alignment to find the next *real* segment
    # in lemma and in plural
    next_lemma_char = None
    next_plural_char = None

    j = aidx + 1
    while j < len(aligned_lemma):
        if next_lemma_char is None and aligned_lemma[j] != "-":
            next_lemma_char = aligned_lemma[j].lower()
        if next_plural_char is None and aligned_plural[j] != "-":
            next_plural_char = aligned_plural[j].lower()

        # Stop once we have a real char on BOTH sides
        if next_lemma_char is not None and next_plural_char is not None:
            break
        j += 1

    # We need an actual *front vowel* in lemma, and the plural must
    # share it *unchanged* (so it's tautomorphemic, not created by suffix)
    if next_lemma_char in {"i", "e"} and next_plural_char == next_lemma_char:
        return True

    return False


def detect_orth_change_dynamic(lemma: str, plural: str) -> str:
    """Detect minimal orthographic change via alignment.

    Returns "X→Y" where X is lemma segment, Y is plural segment.
    Returns "none" when lemma and plural are identical strings.
    Examples: "copac"→"copaci" = "c→ci", "om"→"oameni" = "om→oamen"
    """
    if not lemma or not plural:
        return ""

    # Check if lemma and plural are identical
    if lemma == plural:
        return "no change"

    aligned_lemma, aligned_plural = needleman_wunsch(lemma, plural)

    diff_cols = [
        i
        for i, (l_ch, p_ch) in enumerate(zip(aligned_lemma, aligned_plural))
        if l_ch != p_ch
    ]
    if not diff_cols:
        return "none"

    start = min(diff_cols)
    end = max(diff_cols) + 1

    lemma_window = aligned_lemma[start:end].replace("-", "")
    plural_window = aligned_plural[start:end].replace("-", "")

    if not lemma_window and not plural_window:
        return ""

    # Context expansion for palatalization: include preceding consonant
    expanded_for_palatalization = False

    # Pure insertion (e.g., "c-" → "ci")
    if not lemma_window and len(plural_window) <= 2 and start > 0:
        start -= 1
        lemma_window = aligned_lemma[start:end].replace("-", "")
        plural_window = aligned_plural[start:end].replace("-", "")
        expanded_for_palatalization = True

    # Vowel change with preceding target consonant
    elif (
        lemma_window
        and plural_window
        and len(lemma_window) <= 2
        and len(plural_window) <= 2
        and start > 0
    ):
        preceding_char = aligned_lemma[start - 1]
        if preceding_char != "-" and preceding_char in TARGET_CONSONANTS:
            start -= 1
            lemma_window = aligned_lemma[start:end].replace("-", "")
            plural_window = aligned_plural[start:end].replace("-", "")
            expanded_for_palatalization = True

    # Trimming: preserve pattern if expanded for palatalization
    if expanded_for_palatalization:
        lemma_core = lemma_window
        plural_core = plural_window
    else:
        # Don't over-trim when one side is prefix of other
        if lemma_window and plural_window:
            if lemma_window == plural_window[: len(lemma_window)]:
                return f"{lemma_window}→{plural_window}"
            if plural_window == lemma_window[: len(plural_window)]:
                return f"{lemma_window}→{plural_window}"

        # Conservative prefix trimming
        i = 0
        while (
            i < len(lemma_window) - 1
            and i < len(plural_window) - 1
            and lemma_window[i] == plural_window[i]
        ):
            i += 1

        lemma_core = lemma_window[i:]
        plural_core = plural_window[i:]

    if not lemma_core and not plural_core:
        return ""

    if not lemma_core:
        return f"∅→{plural_core}"
    if not plural_core:
        return f"{lemma_core}→∅"
    if lemma_core == plural_core:
        return ""

    return f"{lemma_core}→{plural_core}"



def derive_mutation_and_orth_change(row: Dict[str, str]) -> None:
    """Derive mutation and orth_change via alignment (non-circular approach).

    1. DISCOVER: Compute orth_change from alignment
    2. CLASSIFY: Check if pattern matches palatalization via suffix matching

    IMPORTANT: Mutation can only be True if stem_final is populated.
    Without knowing the stem-final consonant, we can't determine if it mutated.
    """
    pos = (row.get("pos", "") or "").upper()
    lemma = (row.get("lemma", "") or "").strip().lower()
    plural = (row.get("plural", "") or "").strip().lower()
    stem_final = (row.get("stem_final", "") or "").strip()

    row["mutation"] = "False"
    row["orth_change"] = ""

    if pos not in {"N", "ADJ"}:
        return

    if not lemma:
        return

    # If plural is empty, mark it explicitly
    if not plural:
        row["orth_change"] = "no plural"
        return

    # STEP 0: Mutation requires a known stem_final consonant
    # Without stem_final, we can't track whether the stem mutated
    # (palatalized consonants might just be in the suffix, e.g., -ștri)
    if not stem_final:
        row["mutation"] = "False"
        # Still compute orth_change for reference
        orth_change = detect_orth_change_dynamic(lemma, plural)
        if orth_change in {"∅→iur", "∅→riu"} and plural.endswith("uri"):
            orth_change = "∅→uri"
        row["orth_change"] = orth_change
        return

    # STEP 1: Discover orth_change dynamically
    orth_change = detect_orth_change_dynamic(lemma, plural)

    # Fix ∅→iur and ∅→riu typos (should be ∅→uri for -uri plurals)
    if orth_change in {"∅→iur", "∅→riu"} and plural.endswith("uri"):
        orth_change = "∅→uri"

    row["orth_change"] = orth_change

    if not orth_change or orth_change == "no change":
        row["mutation"] = "False"
        return

    # STEP 2: Classify as palatalization
    # Strategy: Check if plural side contains NEWLY INTRODUCED
    # palatalized consonants.
    # Key: The palatalized consonant must NOT already be in lemma
    # This avoids false positives like "ocinaș→ocinașe" (ș already there)
    is_palatalization = False

    orth_parts = orth_change.split("→")
    if len(orth_parts) == 2:
        orth_from, orth_to = orth_parts

        # Direct approach: Does the plural side contain palatalized consonants
        # that are NOT in the lemma side?
        # Note: "j" only counts as palatalization if stem_final="z" (z→j)
        # For other consonants, "j" appears for orthographic reasons
        # (e.g., giu→je)
        palatalized_consonants = ["ț", "ș"]
        for pal in palatalized_consonants:
            if pal in orth_to and pal not in orth_from:
                is_palatalization = True
                break

        # Special case: "j" only for z→j palatalization
        if not is_palatalization and "j" in orth_to and "j" not in orth_from:
            if stem_final == "z":
                is_palatalization = True

        # Check for z (d→z palatalization) - z must be new in plural
        if not is_palatalization:
            if (
                "z" in orth_to
                and "z" not in orth_from
                and orth_from
                and "d" in orth_from
            ):
                is_palatalization = True

        # Fallback: Try exact pattern matching for edge cases
        # (c→ci, g→gi, etc.)
        # Guard: if the orth_change leading consonant doesn't match
        # stem_final, this is a cluster artifact (e.g., scă→sce
        # expanded to ca→ce, but stem_final is 's' not 'c').
        if not is_palatalization and orth_change in ORTH_TO_PALATAL_IPA:
            lead_from = orth_from.lstrip("aăâîeioușțșț")[:1]
            if not lead_from or lead_from == stem_final:
                is_palatalization = True

        # Fallback: Try suffix matching for patterns like "ate→ăți"
        if not is_palatalization:
            for canonical_pattern in ORTH_TO_PALATAL_IPA:
                canon_parts = canonical_pattern.split("→")
                if len(canon_parts) == 2:
                    canon_from, canon_to = canon_parts
                    if orth_from.endswith(canon_from) and orth_to.endswith(
                        canon_to
                    ):
                        # Same guard: leading consonant must match
                        cf_lead = canon_from.lstrip("aăâîeioușțșț")[:1]
                        if not cf_lead or cf_lead == stem_final:
                            is_palatalization = True
                            break

        # Special check: c/g followed by i/e means palatalization
        # e.g., g→ger (liturg→liturger), c→cer, etc.
        # In Romanian, c/g + front vowel = palatalized
        if not is_palatalization and stem_final in ("c", "g"):
            if orth_from == stem_final:
                # Check if orth_to contains stem_final + i or e
                if stem_final + "i" in orth_to or stem_final + "e" in orth_to:
                    is_palatalization = True

    # STEP 3: Check for frontstem NDE before marking as mutation
    # e.g., borci → borcii (lemma already has ci, just adding -i)
    if is_palatalization and orth_change in {"c→ci", "c→ce", "g→gi", "g→ge"}:
        if lemma.endswith(("ci", "ce", "gi", "ge")):
            # This is frontstem NDE, not a true alternation
            is_palatalization = False

    row["mutation"] = "True" if is_palatalization else "False"


def derive_opportunity(row: Dict[str, str]) -> None:
    """Derive opportunity: does plural add front vowel after stem_final?

    Logic captures both actual mutations AND potential opportunities:
    1. If mutation=True: extract vowel from orth_change pattern
    2. If mutation=False: check if i/e immediately follows stem_final
    3. Special case: 'uri' plurals → opportunity='uri'
    """
    pos = row.get("pos", "")
    lemma = row.get("lemma", "")
    plural = row.get("plural", "")
    stem_final = row.get("stem_final", "")
    orth_change = row.get("orth_change", "")
    mutation = row.get("mutation", "False")

    row["opportunity"] = "none"

    if pos not in {"N", "ADJ"} or not lemma or not plural or not stem_final:
        return

    # If lemma and plural are identical, no opportunity exists
    if lemma == plural:
        return

    # Filter unreliable plurals
    notes = row.get("notes", "")
    if notes:
        notes_lower = notes.lower()
        if "needs plural confirmation" in notes_lower:
            return

    # Check for 'uri' in orth_change
    if orth_change and "uri" in orth_change:
        row["opportunity"] = "uri"
        return

    # CASE 1: mutation=True → extract vowel from orth_change
    if mutation == "True" and orth_change:
        parts = orth_change.split("→", 1)
        if len(parts) == 2:
            plural_side = parts[1]
            # Check what vowel appears in the plural side
            if "i" in plural_side or plural_side.endswith(
                ("ți", "și", "zi", "di", "ci", "gi", "ști")
            ):
                row["opportunity"] = "i"
                return
            elif "e" in plural_side or plural_side.endswith(
                ("țe", "șe", "ze", "de", "ce", "ge", "ște")
            ):
                row["opportunity"] = "e"
                return
            # CASE 1b: Bare consonant change (e.g., t→ț without vowel)
            # Look at what follows the palatalized consonant in plural
            elif plural_side in ("ț", "ș", "j", "z"):
                # Find palatalized consonant in plural and check what follows
                if stem_final and plural:
                    # Look for stem_final position in plural
                    # (accounting for palatalization)
                    palatal_map = {
                        "t": "ț",
                        "c": "c",
                        "g": "g",
                        "s": "ș",
                        "z": "j",
                        "d": "z",
                    }
                    palatal = palatal_map.get(stem_final, stem_final)
                    if palatal in plural:
                        idx = plural.rfind(palatal)
                        if idx != -1 and idx + 1 < len(plural):
                            following = plural[idx + 1]
                            if following == "i":
                                row["opportunity"] = "i"
                                return
                            elif following == "e":
                                row["opportunity"] = "e"
                                return

    # CASE 2: mutation=False → check if i/e immediately after stem_final
    # This captures potential opportunities that didn't palatalize
    # NOTE: We no longer filter vowel-only changes here - those will be
    # detected as opportunity=i/e and then classified as NDE:GIMPE later
    if mutation == "False":
        # Use alignment to find what vowel follows stem_final in plural
        # This prevents false positives from irrelevant i/e elsewhere
        aligned_lemma, aligned_plural = needleman_wunsch(
            lemma.lower(), plural.lower()
        )

        # Map lemma indices to aligned positions
        lemma_to_aligned: Dict[int, int] = {}
        li = 0
        for ai, ch in enumerate(aligned_lemma):
            if ch != "-":
                lemma_to_aligned[li] = ai
                li += 1

        # Find rightmost stem_final in lemma
        try:
            target_idx = max(
                i for i, ch in enumerate(lemma.lower()) if ch == stem_final
            )
        except ValueError:
            return

        if target_idx not in lemma_to_aligned:
            return

        aidx = lemma_to_aligned[target_idx]

        # Walk right to find next segment in plural
        next_plural_char = None
        j = aidx + 1
        while j < len(aligned_plural):
            if aligned_plural[j] != "-":
                next_plural_char = aligned_plural[j].lower()
                break
            j += 1

        # Set opportunity based on the actual vowel following stem_final
        if next_plural_char == "i":
            row["opportunity"] = "i"
            return
        elif next_plural_char == "e":
            row["opportunity"] = "e"
            return


def explode_pipe_group(
    row: Dict[str, str],
    main_field: str,
    companion_fields: List[str],
    sep: str = "|",
) -> List[Dict[str, str]]:
    """Explode pipe-separated field into multiple rows."""

    def _split(raw_val: str) -> List[str]:
        raw_val = (raw_val or "").strip()
        if not raw_val:
            return []
        return [seg.strip() for seg in raw_val.split(sep) if seg.strip()]

    raw = row.get(main_field, "")
    items = _split(raw)
    if not items:
        for name in companion_fields:
            vals = _split(row.get(name, ""))
            row[name] = vals[0] if vals else ""
        return [row]

    def split_field(name: str) -> List[str]:
        return _split(row.get(name, ""))

    companion_lists = [split_field(name) for name in companion_fields]
    n = len(items)

    if n <= 1:
        row[main_field] = items[0]
        for name, vals in zip(companion_fields, companion_lists):
            row[name] = vals[0] if vals else ""
        return [row]

    def pad_to(lst: List[str], length: int) -> List[str]:
        if len(lst) < length:
            lst = lst + [""] * (length - len(lst))
        return lst[:length]

    companion_lists = [pad_to(vals, n) for vals in companion_lists]
    exploded: List[Dict[str, str]] = []
    for idx, item in enumerate(items):
        new_row = dict(row)
        new_row[main_field] = item
        for name, vals in zip(companion_fields, companion_lists):
            new_row[name] = vals[idx]
        exploded.append(new_row)
    return exploded


def explode_derived_verbs_row(row: Dict[str, str]) -> List[Dict[str, str]]:
    """Explode pipe-separated derived verbs into separate rows."""
    return explode_pipe_group(
        row,
        main_field="derived_verbs",
        companion_fields=["deriv_suffixes", "ipa_derived_verbs"],
    )



def normalize_unicode_g2p(s: str) -> str:
    """Normalize Unicode for G2P (cedilla → comma-below diacritics)."""
    if not isinstance(s, str):
        return ""
    s = unicodedata.normalize("NFC", s)
    s = s.replace("ş", "ș").replace("Ş", "Ș")
    s = s.replace("ţ", "ț").replace("Ţ", "Ț")
    return s


IPA_RULES = [
    (r"che", "ke"),
    (r"chi", "ki"),
    (r"ghe", "ɡe"),
    (r"ghi", "ɡi"),
    (r"x", "ks"),
    # ce/ge before another vowel: e is a palatalization marker, not a vowel.
    # Must come BEFORE the ea diphthong rule so the e is consumed first.
    (r"[cC]e(?=[aăâîoui])", "tʃ"),
    (r"[gG]e(?=[aăâîoui])", "dʒ"),
    (r"oa", "o̯a"),
    (r"ea", "e̯a"),
    (r"[cC](?=[eéií])", "tʃ"),
    (r"[gG](?=[eéií])", "dʒ"),
    (r"[cC]", "k"),
    (r"[gG]", "ɡ"),
    (r"[șŞșȘ]", "ʃ"),
    (r"[țŢțŢ]", "ts"),
    (r"j", "ʒ"),
    (r"â|î", "ɨ"),
    (r"ă", "ə"),
    (r"a", "a"),
    (r"e", "e"),
    (r"i", "i"),
    (r"o", "o"),
    (r"u", "u"),
]


def to_ipa(word: str) -> str:
    """Broad Romanian G2P conversion."""
    if not isinstance(word, str) or not word:
        return ""
    w = normalize_unicode_g2p(word).lower()
    for pat, repl in IPA_RULES:
        w = re.sub(pat, repl, w)
    return w


DIMINUTIVE_J_SUFFIXES = ("aică", "oaică", "uică", "eică", "iică")


def tweak_nominal_ipa(lemma: str, ipa: str) -> str:
    """Adjust IPA for diminutive suffixes (i+kə → j+kə)."""
    if not isinstance(lemma, str) or not isinstance(ipa, str):
        return ipa
    lemma_norm = normalize_unicode_g2p(lemma).lower()
    if lemma_norm.endswith(DIMINUTIVE_J_SUFFIXES):
        ipa = re.sub(r"i(kə)$", r"j\1", ipa)
    return ipa


ORDERED_LEMMA_SUFFIXES = ["ică", "iști", "ice", "ist", "esc", "ic", "el"]


def derive_palatal_consonant_pl(row: Dict[str, str]) -> None:
    """Derive palatal_consonant_pl from orth_change for mutation=True rows.

    Strategy:
      1. Use ORTH_TO_PALATAL_IPA for exact/suffix matches on orth_change.
      2. If that fails but mutation=True, fall back to a segment-based map
         from stem_final → palatal IPA (handles 'e→ăți'-type windows).
    """
    row["palatal_consonant_pl"] = ""

    # Only defined for rows that actually palatalize
    if row.get("mutation") != "True":
        return

    orth_change = (row.get("orth_change") or "").strip()
    if not orth_change:
        return

    # 1. DIRECT LOOKUP: exact orth_change in the canonical map
    direct = ORTH_TO_PALATAL_IPA.get(orth_change)
    if direct is not None:
        row["palatal_consonant_pl"] = direct
        return

    # 2. FLEXIBLE MATCH: canonical pattern as suffix / substring
    #    Example: "ate→ăți" should match "te→ți"
    try:
        orth_from, orth_to = orth_change.split("→", 1)
    except ValueError:
        orth_from, orth_to = "", ""

    if orth_from and orth_to:
        for pattern, ipa_value in ORTH_TO_PALATAL_IPA.items():
            try:
                canon_from, canon_to = pattern.split("→", 1)
            except ValueError:
                continue

            # lemma side: window must end in canonical source
            if not orth_from.endswith(canon_from):
                continue

            # plural side: canonical target just needs to appear somewhere
            if canon_to in orth_to:
                row["palatal_consonant_pl"] = ipa_value
                return

    # 3. STEM-BASED FALLBACK: we know *which* segment mutated
    #    (this catches 'e→ăți' cases like afectuozitate→afectuozitatăți)
    stem_final = (row.get("stem_final") or "").strip()

    STEM_TO_IPA = {
        "c": "t͡ʃ",
        "g": "d͡ʒ",
        "t": "t͡s",
        "d": "z",
        "s": "ʃ",
        "z": "ʒ",
    }

    ipa_fallback = STEM_TO_IPA.get(stem_final)
    if ipa_fallback is not None:
        row["palatal_consonant_pl"] = ipa_fallback
        # no return needed; function ends here anyway


def derive_lemma_suffix(row: Dict[str, str]) -> None:
    """Derive lemma_suffix field (e.g., "-ică", "-ist")."""

    def _normalize_suffix(suffix: str) -> str:
        return unicodedata.normalize("NFC", suffix.lower())

    lemma = row.get("lemma", "")
    row["lemma_suffix"] = ""

    if not lemma:
        return

    lemma_normalized = _normalize_suffix(lemma)

    for suffix in ORDERED_LEMMA_SUFFIXES:
        suffix_normalized = _normalize_suffix(suffix)
        if lemma_normalized.endswith(suffix_normalized):
            row["lemma_suffix"] = f"-{suffix}"
            return


def derive_target_is_suffix(row: Dict[str, str]) -> None:
    """Check if target consonant is within tracked suffix."""
    lemma = row.get("lemma", "")
    stem_final = row.get("stem_final", "")
    cluster = row.get("cluster", "")
    lemma_suffix = row.get("lemma_suffix", "")
    row["target_is_suffix"] = "False"
    if not lemma or not stem_final or not lemma_suffix:
        return
    suffix_text = lemma_suffix[1:]
    suffix_start = len(lemma) - len(suffix_text)
    if cluster:
        if cluster in VELAR_FRONT_SEQUENCES:
            target_start = len(lemma) - len(cluster)
            target_end = len(lemma)
        else:
            stem = strip_final_vowel(lemma)
            target_start = len(stem) - len(cluster)
            target_end = len(stem)
    else:
        stem = strip_final_vowel(lemma)
        target_start = -1
        for i in range(len(stem) - 1, -1, -1):
            if stem[i] == stem_final:
                target_start = i
                break
        if target_start == -1:
            return
        target_end = target_start + 1
    if target_start >= suffix_start and target_end <= len(lemma):
        row["target_is_suffix"] = "True"


def derive_derived_verbs_fields(row: Dict[str, str]) -> None:
    """
    Derive deriv_suffixes and ipa_derived_verbs,
    keep only -a/-i/-ui verbs.
    """
    derived_verbs = (row.get("derived_verbs", "") or "").strip()
    row["deriv_suffixes"] = ""
    row["ipa_derived_verbs"] = ""
    if not derived_verbs:
        return

    verb_list = [v.strip() for v in derived_verbs.split("|") if v.strip()]
    if not verb_list:
        return

    clean_verbs: List[str] = []
    suffixes: List[str] = []
    ipa_list: List[str] = []
    for verb in verb_list:
        norm = normalize_unicode_g2p(verb).lower()
        if norm.endswith("ui"):
            suf = "-ui"
        elif norm.endswith("i"):
            suf = "-i"
        elif norm.endswith("a"):
            suf = "-a"
        else:
            continue

        clean_verbs.append(verb)
        suffixes.append(suf)
        ipa = to_ipa(verb)
        if _ipa_normalizer:
            ipa = _ipa_normalizer(ipa)
        ipa_list.append(ipa)
    if not clean_verbs:
        row["derived_verbs"] = ""
        row["deriv_suffixes"] = ""
        row["ipa_derived_verbs"] = ""
        return

    row["derived_verbs"] = "|".join(clean_verbs)
    row["deriv_suffixes"] = "|".join(suffixes)
    row["ipa_derived_verbs"] = "|".join(ipa_list)


def derive_derived_adj_fields(row: Dict[str, str]) -> None:
    """Derive ipa_derived_adj."""
    derived_adj = (row.get("derived_adj", "") or "").strip()
    row["ipa_derived_adj"] = ""
    if not derived_adj:
        return

    adj_list = [a.strip() for a in derived_adj.split("|") if a.strip()]
    if not adj_list:
        return

    ipa_list: List[str] = []
    for adj in adj_list:
        ipa = to_ipa(adj)
        if _ipa_normalizer:
            ipa = _ipa_normalizer(ipa)
        ipa_list.append(ipa)

    row["derived_adj"] = "|".join(adj_list)
    row["ipa_derived_adj"] = "|".join(ipa_list)


def derive_nde_class(row: Dict[str, str]) -> None:
    """Classify NDE (non-derived environment) patterns.

    Types (checked in order):
    - paduchi: clearly derived underapplication (e→i, but no mutation)
    - ochi: lemma=plural with final front-vowel sequence
            (ambiguous plural -i)
    - gimpe: C+front vowel tautomorphemic in root

    Extended to dorsals (c/g) AND coronals (t/d/s/z).
    Per email.txt definitions and NDEB explanation.
    """
    pos = row.get("pos", "")
    lemma = (row.get("lemma", "") or "").strip().lower()
    plural = (row.get("plural", "") or "").strip().lower()
    cluster = (row.get("cluster", "") or "").strip().lower()
    mutation = str(row.get("mutation", ""))
    stem_final = (row.get("stem_final", "") or "").strip().lower()
    row["nde_class"] = ""

    # Only nouns and adjectives with a plural, and ONLY if they did
    # *not* mutate
    if pos not in {"N", "ADJ"} or not lemma or not plural:
        return
    if mutation == "True":
        return

    # ------------------------------------------------------------------
    # 1. PADUCHI: clearly derived underapplication (most specific)
    # ------------------------------------------------------------------

    # DORSAL: classic păduche-type (che/ghe → chi/ghi or chiuri/ghiuri)
    if lemma.endswith("che") and (
        plural.endswith("chi") or plural.endswith("chiuri")
    ):
        row["nde_class"] = "paduchi"
        return
    if lemma.endswith("ghe") and (
        plural.endswith("ghi") or plural.endswith("ghiuri")
    ):
        row["nde_class"] = "paduchi"
        return

    # CORONAL: C+e → C+i with no palatalization
    # (t/d/s/z + e → t/d/s/z + i)
    # Also check palatal variants where alternation exists:
    #   - lemma has plain C+e, plural has palatal+i (e.g., te→ți)
    if stem_final in {"t", "d", "s", "z"}:
        # Plain coronal: C+e → C+i (no mutation)
        if lemma.endswith(stem_final + "e") and plural.endswith(
            stem_final + "i"
        ):
            row["nde_class"] = "paduchi"
            return

        # Palatal variant: C+e → PAL+i (where PAL is palatalized form)
        # This catches cases where the e→i shift IS accompanied by
        # palatalization, but it's NDEB because C+e was in lemma
        palatal = PALATAL_SURFACE.get(stem_final)
        if palatal:
            if lemma.endswith(stem_final + "e") and plural.endswith(
                palatal + "i"
            ):
                row["nde_class"] = "paduchi"
                return

    # ------------------------------------------------------------------
    # 2. OCHI: lemma = plural with ambiguous final -i
    # ------------------------------------------------------------------
    if lemma == plural:
        # DORSAL: chi/ghi as in *ochi*
        if cluster in {"chi", "ghi"}:
            row["nde_class"] = "ochi"
            return

        # CORONAL: lemma = plural, final C+i (t/d/s/z + i)
        # with no palatalization
        # Also check palatal variants (ți, și, zi, ji) where they appear
        # in BOTH lemma and plural (no evidence of alternation)
        if stem_final in {"t", "d", "s", "z"}:
            # Plain coronal + i
            if lemma.endswith(stem_final + "i"):
                row["nde_class"] = "ochi"
                return

            # Palatal coronal + i (both forms identical with palatal)
            # This is ONLY OCHI if lemma=plural AND no alternation
            palatal = PALATAL_SURFACE.get(stem_final)
            if palatal and lemma.endswith(palatal + "i"):
                row["nde_class"] = "ochi"
                return

    # ------------------------------------------------------------------
    # 3. GIMPE: tautomorphemic C+front vowel sequences (broadest)
    # ------------------------------------------------------------------
    # Applies to ALL target consonants (dorsals + coronals).

    # GIMPE 3a: Special case for dorsal frontstem patterns
    # These are canonical GIMPE: the palatal is already in the lemma,
    # and the plural just adds vowels or shifts them.
    # Examples:
    #   - ce→ci: cruce→cruci, indice→indici (e→i vowel shift)
    #   - ci→cii: borci→borcii (tautomorphemic ci, plural adds i)
    #   - ge→gi: similar pattern with /dʒ/
    if stem_final in {"c", "g"}:
        if lemma.endswith("ce") and plural.endswith("ci"):
            row["nde_class"] = "gimpe"
            return
        elif lemma.endswith("ci") and plural.endswith("cii"):
            row["nde_class"] = "gimpe"
            return
        elif lemma.endswith("ge") and plural.endswith("gi"):
            row["nde_class"] = "gimpe"
            return
        elif lemma.endswith("gi") and plural.endswith("gii"):
            row["nde_class"] = "gimpe"
            return

    # GIMPE 3a-suffix: Diminutive suffix-internal patterns
    # Diminutive -ică → -icici: C+i is inside suffix, not root-final
    # Examples: budăcică→budăcicici, drăgucică→drăgucicici
    # The "target consonant" is in the suffix, so this is tautomorphemic
    if lemma.endswith("ică") and plural.endswith("icici"):
        row["nde_class"] = "gimpe"
        return

    # GIMPE 3a-vowel-only: Vowel-only alternations (no consonant change)
    # When orth_change is purely vocalic (e→i, u→i, a→le, etc.),
    # the stem_final consonant is the same in both forms.
    # The C+front-vowel is tautomorphemic - not created by plural suffix.
    # Examples:
    #   - abazie→abazii (e→i): z is same, z+i is tautomorphemic
    #   - abagiu→abagii (u→i): g is same, g+i is tautomorphemic
    # Guard: only apply if stem_final is actually adjacent to a front
    # vowel in the lemma.  Words like butoi→butoaie (t followed by o,
    # not a front vowel) should NOT be classified as gimpe.
    orth_change = row.get("orth_change", "")
    if orth_change_is_vowels_or_le(orth_change):
        sf_fv_pattern = stem_final + "[ie]"
        if re.search(sf_fv_pattern, lemma):
            row["nde_class"] = "gimpe"
            return

    # GIMPE 3a-suffix-pattern: Suffix-internal patterns like -istă→-iste
    # The target consonant is inside a derivational suffix, not root-final.
    # Examples:
    #   - absolutistă→absolutiste (tă→te): stem_final 's' is in -istă suffix
    #   - t→te, tă→te where lemma ends in -ist/-istă
    # These are out of domain for root-final palatalization
    if stem_final in {"s", "t"}:
        if lemma.endswith("istă") and plural.endswith("iste"):
            row["nde_class"] = "gimpe"
            return
        elif lemma.endswith("ist") and plural.endswith("iste"):
            row["nde_class"] = "gimpe"
            return

    # GIMPE 3b: General alignment-based detection
    # We use alignment to ensure the front vowel after stem_final is
    # *already* present in the lemma and is shared identically by the plural
    # at that position (i.e., not created by suffixation).
    #
    # This catches cases like:
    #   - abazie → abazii: lemma has 'zie', plural has 'zii',
    #     alignment shows 'i' is shared after 'z' ✓
    #   - alice → alice: lemma=plural, both have 'ce' ✓
    #   - abstinentă → abstinente: lemma has 'ti' internally,
    #     but plural creates NEW 'te' at end ✗ (not GIMPE)
    if has_tautomorphemic_front_vowel(lemma, plural, stem_final):
        row["nde_class"] = "gimpe"
        return


def fix_nde_mutations(row: Dict[str, str]) -> None:
    """Fix mutation status for NDE items (should always be False)."""
    nde_class = row.get("nde_class", "")

    if nde_class:
        # All NDE items should have mutation=False
        row["mutation"] = "False"
        # Clear palatal consonant since it's not a true alternation
        row["palatal_consonant_pl"] = ""



def lemma_has_palatal_coronal(lemma: str, stem_final: str) -> bool:
    """
    Detect if lemma already has the palatal output grapheme for a coronal.

    This identifies lemmas where the "target consonant" is already in its
    palatalized form in the lemma, analogous to frontstem dorsals.

    Examples:
        abația (has 'ț') with stem_final='t' → True (already palatal)
        batuc (no palatal graphemes) with stem_final='c' → False

    For stem_final='d', the palatal form is 'z', which is ambiguous with
    underlying /z/. To avoid false positives we require 'z' + front vowel
    (zi/ze) near the end, matching the palatalized context.

    Returns True if palatal grapheme appears near end of lemma.
    """
    if not lemma or stem_final not in {"t", "d", "s", "z"}:
        return False

    pal = PALATAL_SURFACE.get(stem_final)
    if not pal:
        return False

    tail = lemma[-4:]

    # For d→z, bare 'z' is ambiguous with underlying /z/.
    # Require z+front vowel (zi/ze) to indicate palatalized context.
    if stem_final == "d":
        return bool(re.search(r"z[ie]", tail))

    return pal in tail


def fix_underlying_palatals(row: Dict[str, str]) -> None:
    """Remove non-NDE items with underlying palatals from domain.

    This handles two cases:

    1. Frontstem dorsals (ci/ce/gi/ge): already-palatalized stems
    2. Palatal coronals (lemma has ț/ș/j/z): underlying palatals

    These are excluded from palatalization domain UNLESS used as NDEB,
    because they don't represent alternating underspecified segments.

    IAP perspective: Fully specified /k, g, t, d, s, z, ʃ, ʒ, ts, dz/
    don't alternate. Only underspecified /K, G, T, S, Z/ get feature-filled
    in front-vowel contexts.

    Must be called AFTER derive_nde_class but BEFORE
    derive_exception_reason.
    """
    stem_final = row.get("stem_final", "")
    frontstem = row.get("frontstem_dorsal", "")
    nde_class = row.get("nde_class", "")
    lemma = (row.get("lemma", "") or "").strip().lower()

    # Case 1: Frontstem dorsals (ci/ce/gi/ge)
    if frontstem == "True" and stem_final in {"c", "g"}:
        # If this item is NOT being used as NDE, remove from domain
        if not nde_class:
            row["stem_final"] = ""
            row["cluster"] = ""
            row["opportunity"] = "none"
            row["mutation"] = "False"
            row["palatal_consonant_pl"] = ""
            row["orth_change"] = ""
            return

    # Case 2: Palatal coronals already in lemma (ț/ș/j/z in palatal sense)
    if stem_final in {"t", "d", "s", "z"}:
        if lemma_has_palatal_coronal(lemma, stem_final) and not nde_class:
            # Lemma already has the palatal grapheme → underlyingly palatal
            # Exclude from alternation domain
            row["stem_final"] = ""
            row["cluster"] = ""
            row["opportunity"] = "none"
            row["mutation"] = "False"
            row["palatal_consonant_pl"] = ""
            row["orth_change"] = ""
            return


def mark_tp_domain(row: Dict[str, str]) -> None:
    """
    Mark whether this row belongs in the TP domain for segment-level counts.

    TP domain = noun, i/e environment, target segment in {c,g,t,d,s,z},
    EXCLUDING:
    - All NDEB classes (gimpe, ochi, paduchi)
    - Items already marked as "non exception" (structural out-of-domain)
    - Suffix-internal targets (optional, currently excluded)

    IAP perspective: The TP domain contains only items where the target segment
    could plausibly be an underspecified /K, G, T, S, Z/ that alternates.
    Fully specified segments, NDEB items, and structural non-alternators are
    excluded from both N and exception counts.

    Must be called AFTER derive_exception_reason.
    """
    pos = row.get("pos", "")
    opportunity = row.get("opportunity", "")
    stem_final = row.get("stem_final", "")
    nde_class = row.get("nde_class", "")
    exception_reason = row.get("exception_reason", "")
    target_is_suffix = row.get("target_is_suffix", "False")

    # Basic domain check: noun with i/e opportunity and target segment
    in_basic_domain = (
        pos == "N"
        and opportunity in {"i", "e"}
        and stem_final in {"c", "g", "t", "d", "s", "z"}
    )

    if not in_basic_domain:
        row["tp_in_domain"] = "False"
        return

    # Exclude all NDEB classes from the productive TP domain
    # These are structured families, not productive exceptions
    if nde_class in NDE_CLASSES:
        row["tp_in_domain"] = "False"
        return

    # Exclude anything marked as "non exception"
    # (underlying palatal, suffix-internal, weird morphophonology, etc.)
    if exception_reason == "non exception":
        row["tp_in_domain"] = "False"
        return

    # Optionally exclude suffix-internal targets
    # These are targets where the consonant is part of a derivational suffix,
    # not the root, so they don't reflect root-final palatalization
    if target_is_suffix == "True":
        row["tp_in_domain"] = "False"
        return

    # If we made it here, this item is in the productive TP domain
    row["tp_in_domain"] = "True"


def derive_exception_reason(row: Dict[str, str]) -> None:
    """Derive exception_reason field for mutation behavior.

    Categories:
    - undergoer: mutation=True (word palatalized)
    - nde:{type}: mutation=False with known NDE explanation
    - unexplained: mutation=False with i/e opportunity, no NDE explanation
    - non exception: all other cases (no opportunity, no data, etc.)

    Note: Structural issues (no plural, no stem_final) can be filtered
    using those fields directly - exception_reason focuses on the
    palatalization pattern itself.

    Applies to all POS types (nouns and adjectives). The restriction
    to nouns for statistical analysis happens in the R code, not here.
    """
    row["exception_reason"] = ""
    lemma = row.get("lemma", "")
    plural = row.get("plural", "")
    mutation = row.get("mutation", "")
    opportunity = row.get("opportunity", "")
    nde = row.get("nde_class", "")
    stem_final = row.get("stem_final", "")

    cluster = (row.get("cluster", "") or "").strip().lower()

    # No stem final or no plural → non exception (data limitation)
    if not stem_final or not plural:
        row["exception_reason"] = "non exception"
        return

    # Dorsal digraph clusters (che/ghe/chi/ghi) where the plural
    # doesn't preserve the dorsal — irregular truncation, not
    # standard palatalization (e.g., dilimache → dilii).
    if cluster in {"che", "ghe", "chi", "ghi"} and plural:
        plural_l = plural.lower()
        if not any(
            plural_l.endswith(s)
            for s in ("chi", "ghi", "che", "ghe", "ci", "ce", "gi", "ge")
        ):
            row["exception_reason"] = "non exception"
            return

    # Suppletive suffix replacement: -ică → -ele, -ică → -icele, etc.
    # The target consonant (c from -ică) is entirely replaced in the
    # plural, so there is no C+front-vowel context to evaluate.
    if lemma and plural:
        lemma_l = lemma.lower()
        plural_l = plural.lower()
        if lemma_l.endswith("ică") and not plural_l.endswith(
            ("ici", "ice", "ică")
        ):
            row["exception_reason"] = "non exception"
            return

    # Undergoers: words that palatalized
    if mutation == "True":
        row["exception_reason"] = "undergoer"
        return

    # OCHI, PADUCHI, and GIMPE: NDE cases that apply
    # regardless of opportunity
    # Handle these BEFORE filtering already-palatal stems
    if nde in {"ochi", "paduchi", "gimpe"}:
        row["exception_reason"] = f"nde:{nde}"
        return

    # Already-palatal stems: check for palatal consonants near stem edge.
    # Two cases:
    # 1. Coronals (t,d,s,z): check for the SPECIFIC palatal output of
    #    this stem_final (t→ț, d→z, s→ș, z→j). For d, require z+front
    #    vowel to avoid conflating underlying /z/ with palatalized /d/.
    # 2. Dorsals (c,g): check for ANY palatal sibilant (ș,ț,j) immediately
    #    adjacent to stem_final, indicating an already-palatal cluster
    #    (e.g., ușcă, where ș+c is not a clean sc cluster for analysis).
    if lemma:
        lemma_lower = lemma.lower()
        stem = strip_final_vowel(lemma_lower)
        tail = stem[-4:] if len(stem) >= 4 else stem
        if stem_final in PALATAL_SURFACE:
            pal = PALATAL_SURFACE[stem_final]
            if stem_final == "d":
                if re.search(r"z[ie]", tail):
                    row["exception_reason"] = "non exception"
                    return
            elif pal in tail:
                row["exception_reason"] = "non exception"
                return
        # For ALL target consonants: check if a palatal consonant
        # is immediately before stem_final, indicating an already-palatal
        # cluster (e.g., ușcă, șt in Leberwurst loanwords, jd in dajd)
        if len(stem) >= 2:
            sf_pos = stem.rfind(stem_final)
            if sf_pos > 0 and stem[sf_pos - 1] in {"ș", "ț", "j"}:
                row["exception_reason"] = "non exception"
                return

    # Non-exceptions: no opportunity to mutate (none or uri)
    if opportunity in {"none", "uri"}:
        row["exception_reason"] = "non exception"
        return

    # True exceptions: had opportunity but didn't mutate (unexplained)
    if opportunity in {"i", "e"}:
        row["exception_reason"] = "unexplained"
        return
