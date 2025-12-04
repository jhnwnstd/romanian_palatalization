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
VOWELS = "aăâeiîou"

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
    if lemma[-1] in VOWELS:
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


def _calibrate_plural_thresholds() -> None:
    """Update global plausibility thresholds from observed scores."""
    global _PL_CALIBRATED, _PL_THRESHOLDS  # noqa: F824
    n = len(_PL_PLAUS_SCORES)
    if n < _PL_MIN_CALIBRATION:
        return
    scores = sorted(_PL_PLAUS_SCORES)
    idx_reject = max(0, int(0.05 * (n - 1)))
    idx_border = max(0, int(0.20 * (n - 1)))
    reject_thr = scores[idx_reject]
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

    if (
        plausibility < reject_thr
        or not len_ok
        or (seq_sim < 0.1 and lcs_ratio < 0.1)
    ):
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


# Mutation patterns: stem_final -> list of (lemma_ending, plural_ending)
MUTATION_PATTERNS = {
    "c": [
        ("c", "ci"),
        ("c", "ce"),
        ("că", "ci"),
        ("că", "ce"),
        ("ca", "ce"),
    ],
    "g": [
        ("g", "gi"),
        ("g", "ge"),
        ("gă", "gi"),
        ("gă", "ge"),
        ("ga", "ge"),
        ("go", "gi"),
    ],
    "t": [
        ("t", "ți"),
        ("t", "țe"),
        ("ct", "cți"),
        ("ct", "cțe"),
        ("te", "ți"),
        ("te", "țe"),
        ("tă", "ți"),
        ("tă", "țe"),
    ],
    "d": [
        ("d", "zi"),
        ("de", "zi"),
        ("d", "ze"),
        ("de", "ze"),
        ("de", "di"),
        ("dă", "zi"),
        ("dă", "ze"),
    ],
    "s": [
        ("s", "și"),
        ("s", "șe"),
        ("st", "ști"),
        ("st", "ște"),
        ("sc", "ști"),
        ("sc", "ște"),
        ("scă", "ști"),
        ("scă", "ște"),
    ],
    "z": [
        ("z", "ji"),
        ("z", "je"),
    ],
}


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


def detect_orth_change_dynamic(lemma: str, plural: str) -> str:
    """Detect minimal orthographic change via alignment (non-circular).

    Returns "X→Y" where X is lemma segment, Y is plural segment.
    Examples: "copac"→"copaci" = "c→ci", "om"→"oameni" = "om→oamen"
    """
    if not lemma or not plural:
        return ""

    aligned_lemma, aligned_plural = needleman_wunsch(lemma, plural)

    diff_cols = [
        i
        for i, (l_ch, p_ch) in enumerate(zip(aligned_lemma, aligned_plural))
        if l_ch != p_ch
    ]
    if not diff_cols:
        return ""

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


def get_change_window(
    lemma: str, plural: str, stem_final: str, cluster: str
) -> Tuple[str, str]:
    """Compute minimal change window around stem_final using alignment."""
    if not lemma or not plural or not stem_final:
        return "", ""
    aligned_lemma, aligned_plural = needleman_wunsch(lemma, plural)

    lemma_to_aligned: Dict[int, int] = {}
    lemma_idx = 0
    for aligned_idx, ch in enumerate(aligned_lemma):
        if ch != "-":
            lemma_to_aligned[lemma_idx] = aligned_idx
            lemma_idx += 1

    if cluster:
        if cluster in VELAR_FRONT_SEQUENCES:
            cluster_start = len(lemma) - len(cluster)
            cluster_end = len(lemma)
        else:
            stem = strip_final_vowel(lemma)
            cluster_start = len(stem) - len(cluster)
            cluster_end = len(stem)
        target_indices = list(range(cluster_start, cluster_end))
    else:
        stem = strip_final_vowel(lemma)
        target_idx = -1
        for i in range(len(stem) - 1, -1, -1):
            if stem[i] == stem_final:
                target_idx = i
                break
        if target_idx == -1:
            return "", ""
        target_indices = [target_idx]

    aligned_target_positions: Set[int] = set()
    for idx in target_indices:
        if idx in lemma_to_aligned:
            aligned_target_positions.add(lemma_to_aligned[idx])
    if not aligned_target_positions:
        return "", ""

    diff_positions: Set[int] = set()
    for i, (char_l, char_p) in enumerate(zip(aligned_lemma, aligned_plural)):
        if char_l != char_p:
            diff_positions.add(i)

    all_positions = aligned_target_positions | diff_positions
    if not all_positions:
        return "", ""
    window_start = min(all_positions)
    last_affected = max(all_positions)
    window_end = last_affected + 1
    if window_end < len(aligned_lemma):
        window_end += 1
    lemma_sub = aligned_lemma[window_start:window_end].replace("-", "")
    plural_sub = aligned_plural[window_start:window_end].replace("-", "")
    return lemma_sub, plural_sub


def derive_mutation_and_orth_change(row: Dict[str, str]) -> None:
    """Derive mutation and orth_change via alignment (non-circular approach).

    1. DISCOVER: Compute orth_change from alignment
    2. CLASSIFY: Check if pattern matches palatalization via suffix matching
    """
    pos = (row.get("pos", "") or "").upper()
    lemma = (row.get("lemma", "") or "").strip().lower()
    plural = (row.get("plural", "") or "").strip().lower()

    row["mutation"] = "False"
    row["orth_change"] = ""

    if pos not in {"N", "ADJ"}:
        return

    if not lemma or not plural:
        return

    # STEP 1: Discover orth_change dynamically
    orth_change = detect_orth_change_dynamic(lemma, plural)

    # Fix ∅→iur and ∅→riu typos (should be ∅→uri for -uri plurals)
    if orth_change in {"∅→iur", "∅→riu"} and plural.endswith("uri"):
        orth_change = "∅→uri"

    row["orth_change"] = orth_change

    if not orth_change:
        row["mutation"] = "False"
        return

    # STEP 2: Classify as palatalization via exact or suffix matching
    # Suffix matching: "ate→ăți" matches "te→ți" if "ate" ends
    # with "te" AND "ăți" ends with "ți"
    is_palatalization = False

    if orth_change in ORTH_TO_PALATAL_IPA:
        is_palatalization = True
    else:
        orth_parts = orth_change.split("→")
        if len(orth_parts) == 2:
            orth_from, orth_to = orth_parts
            for canonical_pattern in ORTH_TO_PALATAL_IPA:
                canon_parts = canonical_pattern.split("→")
                if len(canon_parts) == 2:
                    canon_from, canon_to = canon_parts
                    if orth_from.endswith(canon_from) and orth_to.endswith(
                        canon_to
                    ):
                        is_palatalization = True
                        break

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

    # CASE 2: mutation=False → check if i/e immediately after stem_final
    # This captures potential opportunities that didn't palatalize
    if mutation == "False":
        # Palatalization map for checking palatalized forms too
        palatal_map = {
            "t": "ț",
            "d": "d",
            "s": "ș",
            "z": "z",
            "c": "c",
            "g": "g",
        }

        # Check for stem_final + i (plain or palatalized)
        if stem_final + "i" in plural:
            row["opportunity"] = "i"
            return
        elif stem_final in palatal_map:
            palatalized = palatal_map[stem_final]
            if palatalized + "i" in plural:
                row["opportunity"] = "i"
                return

        # Check for stem_final + e (plain or palatalized)
        if stem_final + "e" in plural:
            row["opportunity"] = "e"
            return
        elif stem_final in palatal_map:
            palatalized = palatal_map[stem_final]
            if palatalized + "e" in plural:
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


def explode_derived_adj_row(row: Dict[str, str]) -> List[Dict[str, str]]:
    """Explode pipe-separated derived adjectives into separate rows."""
    return explode_pipe_group(
        row,
        main_field="derived_adj",
        companion_fields=["ipa_derived_adj"],
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
    "de→di": "dʲ",
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
    "sc→ște": "ʃt",
    "sc→ști": "ʃt",
    "scă→ște": "ʃt",
    "scă→ști": "ʃt",
    "st→ște": "ʃt",
    "st→ști": "ʃt",
    "t→țe": "t͡s",
    "t→ți": "t͡s",
    "te→țe": "t͡s",
    "te→ți": "t͡s",
    "tă→țe": "t͡s",
    "tă→ți": "t͡s",
    "z→je": "ʒ",
    "z→ji": "ʒ",
}

ORDERED_LEMMA_SUFFIXES = ["ică", "iști", "ice", "ist", "esc", "ic", "el"]

LEMMA_SUFFIXES = [
    unicodedata.normalize("NFC", suffix) for suffix in ORDERED_LEMMA_SUFFIXES
]


def derive_palatal_consonant_pl(row: Dict[str, str]) -> None:
    """Derive palatal_consonant_pl from orth_change only for mutation=True."""
    row["palatal_consonant_pl"] = ""

    mutation = row.get("mutation", "")
    if mutation != "True":
        return

    orth_change = row.get("orth_change", "")
    if not orth_change:
        return

    # Exact match
    if orth_change in ORTH_TO_PALATAL_IPA:
        row["palatal_consonant_pl"] = ORTH_TO_PALATAL_IPA[orth_change]
        return

    # Suffix matching (e.g., "ate→ăți" matches "te→ți")
    orth_parts = orth_change.split("→")
    if len(orth_parts) == 2:
        orth_from, orth_to = orth_parts
        for canonical_pattern, ipa_value in ORTH_TO_PALATAL_IPA.items():
            canon_parts = canonical_pattern.split("→")
            if len(canon_parts) == 2:
                canon_from, canon_to = canon_parts
                if orth_from.endswith(canon_from) and orth_to.endswith(
                    canon_to
                ):
                    row["palatal_consonant_pl"] = ipa_value
                    return


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
    - gimpe: C+front vowel tautomorphemic in root (ci/ce/gi/ge in lemma)
    - paduchi: Lemma ends in che/ghe, vowel e→i (derived underapplication)
    - ochi: Lemma=plural with chi/ghi clusters (ambiguous morphology)

    Per email.txt definitions and NDEB explanation.
    """
    pos = row.get("pos", "")
    lemma = row.get("lemma", "") or ""
    plural = row.get("plural", "") or ""
    cluster = row.get("cluster", "") or ""
    mutation = str(row.get("mutation", ""))
    stem_final = row.get("stem_final", "")
    row["nde_class"] = ""

    if pos != "N" or not lemma or not plural:
        return

    if mutation == "True":
        return

    # 1. PADUCHI: lemma ends in che/ghe AND vowel changes to i
    # Checked first because it's most specific
    # "Clearly derived underapplication" - SG has che/ghe, PL has chi/ghi
    if lemma.endswith(("che", "ghe")):
        # Check if plural has chi/ghi (vowel e→i)
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

    # 2. OCHI: lemma=plural with chi/ghi clusters
    # "Ambiguous morphology" - could be /oki/ or /ok-i/
    if lemma == plural and cluster in ("chi", "ghi"):
        row["nde_class"] = "ochi"
        return

    # 3. GIMPE: C+front vowel tautomorphemic in root
    # "Canonical NDEB" - ci/ce/gi/ge already in root, not created by suffix
    # Includes both invariant (alice→alice) and variable (abagiu→abagii)
    if stem_final in ("c", "g"):
        if stem_final + "i" in lemma or stem_final + "e" in lemma:
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


def derive_exception_reason(row: Dict[str, str]) -> None:
    """Derive exception_reason field for mutation behavior.

    Categories:
    - undergoer: mutation=True (word palatalized)
    - nde:{type}: mutation=False with known NDE explanation
    - unexplained: mutation=False with i/e opportunity, no NDE explanation
    - non exception: opportunity=none/uri (no chance to mutate) OR adjectives
    """
    row["exception_reason"] = ""
    pos = row.get("pos", "")
    plural = row.get("plural", "")
    mutation = row.get("mutation", "")
    opportunity = row.get("opportunity", "")
    nde = row.get("nde_class", "")

    # ADJ - mark all adjectives as non exception
    if pos == "ADJ":
        row["exception_reason"] = "non exception"
        return

    # Only process nouns beyond this point
    if pos != "N":
        return

    # Skip if no plural form for remaining analyses
    if not plural:
        # Still mark as non exception if opportunity=none/uri
        if opportunity in {"none", "uri"}:
            row["exception_reason"] = "non exception"
        return

    # Undergoers: words that palatalized
    if mutation == "True":
        row["exception_reason"] = "undergoer"
        return

    # NDE classes explain non-mutation (check before "non exception")
    if nde in {"gimpe", "ochi", "paduchi"}:
        row["exception_reason"] = f"nde:{nde}"
        return

    # True exceptions: had opportunity but didn't mutate (unexplained)
    if opportunity in {"i", "e"}:
        row["exception_reason"] = "unexplained"
        return

    # Non-exceptions: no opportunity to mutate (none or uri)
    if opportunity in {"none", "uri"}:
        row["exception_reason"] = "non exception"
        return
