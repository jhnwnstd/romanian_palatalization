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

# IPA normalizer function to be injected from caller
# Set this before calling derive_derived_verbs_fields or
# derive_derived_adj_fields
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

# Consonants used to sanity-check whether a TARGET_CONSONANT is truly
# stem-final (no other consonant to its right in the stem).
ROMANIAN_CONSONANTS: Set[str] = set("bcdfghjklmnpqrstvwxyzșşţțțţ")


# ============================================================================
# STRING UTILITIES
# ============================================================================


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


# ============================================================================
# PLURAL VALIDATION
# ============================================================================

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


# ============================================================================
# PART 1: DERIVATION FUNCTIONS
# ============================================================================
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
    lemma = row.get("lemma", "")
    if not lemma:
        row["stem_final"] = ""
        row["cluster"] = ""
        return

    # 1) Velar front clusters like "chi/che/ghi/ghe" at word edge
    for seq, consonant in VELAR_FRONT_SEQUENCES.items():
        if lemma.endswith(seq):
            row["stem_final"] = consonant
            row["cluster"] = seq
            return

    # 2) Final clusters like -st, -sc, -ct (after stripping final vowel)
    stem = strip_final_vowel(lemma)
    for cluster, consonant in FINAL_CLUSTERS.items():
        if stem.endswith(cluster):
            row["stem_final"] = consonant
            row["cluster"] = cluster
            return

    # 3) Bare single consonant target at the very right edge of the stem
    for i in range(len(stem) - 1, -1, -1):
        ch = stem[i]
        if ch in TARGET_CONSONANTS:
            trailing = stem[i + 1 :]
            # Only accept this as stem_final if there is no other
            # consonant to its right in the stem.
            if any(c in ROMANIAN_CONSONANTS for c in trailing):
                continue
            row["stem_final"] = ch
            row["cluster"] = ""
            return

    row["stem_final"] = ""
    row["cluster"] = ""


# ============================================================================
# PART 2: CONSTANTS
# ============================================================================

# Mutation patterns: stem_final -> list of (lemma_ending, plural_ending)
MUTATION_PATTERNS = {
    "c": [
        ("c", "ci"),  # expected: core c → ci
        ("c", "ce"),  # expected: core c → ce
        ("că", "ci"),  # expected: -că / -ică → -ci
        ("că", "ce"),  # expected: -că / -ică → -ce
        ("ca", "ce"),  # NEW: abaca → abace
    ],
    "g": [
        ("g", "gi"),  # expected: core g → gi
        ("g", "ge"),  # expected: core g → ge
        ("gă", "gi"),  # expected: -gă → -gi
        ("gă", "ge"),  # expected: -gă → -ge
        ("ga", "ge"),  # NEW: ga → ge
        ("go", "gi"),  # NEW: flamingo → flamingi
    ],
    "t": [
        ("t", "ți"),  # expected: core t → ț(i)
        ("t", "țe"),  # unexpected: theoretical t → ț(e)
        ("ct", "cți"),  # expected: ct → cț(i)
        ("ct", "cțe"),  # unexpected: ct → cț(e)
        ("te", "ți"),  # expected: -tate / -tite → -tăți / -ți
        ("te", "țe"),  # unexpected: theoretical te → ț(e)
        ("tă", "ți"),  # expected (rare): tă → ț(i)
        ("tă", "țe"),  # expected (rare): tă → ț(e)
    ],
    "d": [
        ("d", "zi"),  # expected: cald → calzi, crud → cruzi, surd → surzi
        ("de", "zi"),  # expected: verde → verzi, lespede → lespezi, etc.
        ("d", "ze"),  # unexpected: theoretical
        ("de", "ze"),  # unexpected: theoretical
        ("de", "di"),  # expected (rare): nădejde → nădejdi
        ("dă", "zi"),  # expected: amendă → amenzi, izbândă → izbânzi, etc.
        ("dă", "ze"),  # unexpected: theoretical dă → z(e)
    ],
    "s": [
        ("s", "și"),  # expected: s → ș(i)
        ("s", "șe"),  # unexpected: theoretical s → ș(e)
        ("st", "ști"),  # expected: prost → proști, -ist → -iști
        ("st", "ște"),  # unexpected: theoretical st → șt(e)
        ("sc", "ști"),  # expected: -esc → -ești
        ("sc", "ște"),  # unexpected: theoretical -esc → -ește
        ("scă", "ști"),  # unexpected: theoretical scă → ști
        ("scă", "ște"),  # expected (rare): franciscă → franciște
    ],
    "z": [
        ("z", "ji"),  # expected (rare): obraz → obraji, harbuz → harbuji, etc.
        ("z", "je"),  # unexpected: theoretical z → j(e)
    ],
}


# ============================================================================
# PART 2: NEEDLEMAN-WUNSCH ALIGNMENT
# ============================================================================
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


def get_change_window(
    lemma: str, plural: str, stem_final: str, cluster: str
) -> Tuple[str, str]:
    """
    Compute the minimal change window around stem_final.

    Uses Needleman-Wunsch alignment to identify the region where
    lemma and plural differ, focused on the palatalization site.

    Args:
        lemma: Lemma form
        plural: Plural form
        stem_final: Target consonant
        cluster: Associated cluster (if any)

    Returns:
        Tuple of (lemma_sub, plural_sub): Substrings showing the change
    """
    if not lemma or not plural or not stem_final:
        return "", ""
    aligned_lemma, aligned_plural = needleman_wunsch(lemma, plural)
    # Build index mapping lemma -> aligned
    lemma_to_aligned: Dict[int, int] = {}
    lemma_idx = 0
    for aligned_idx, ch in enumerate(aligned_lemma):
        if ch != "-":
            lemma_to_aligned[lemma_idx] = aligned_idx
            lemma_idx += 1
    # Identify target span in lemma
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
    # Map target indices to aligned positions
    aligned_target_positions: Set[int] = set()
    for idx in target_indices:
        if idx in lemma_to_aligned:
            aligned_target_positions.add(lemma_to_aligned[idx])
    if not aligned_target_positions:
        return "", ""
    # Find all mismatch positions in the alignment
    diff_positions: Set[int] = set()
    for i, (char_l, char_p) in enumerate(zip(aligned_lemma, aligned_plural)):
        if char_l != char_p:
            diff_positions.add(i)
    # Change window: all affected positions (target + diffs)
    all_positions = aligned_target_positions | diff_positions
    if not all_positions:
        return "", ""
    window_start = min(all_positions)
    # Add one extra alignment column of right-hand context if available
    last_affected = max(all_positions)
    window_end = last_affected + 1
    if window_end < len(aligned_lemma):
        window_end += 1
    lemma_sub = aligned_lemma[window_start:window_end].replace("-", "")
    plural_sub = aligned_plural[window_start:window_end].replace("-", "")
    return lemma_sub, plural_sub


# ============================================================================
# PART 2: DERIVATION FUNCTIONS
# ============================================================================
def _canonical_orth_change(
    lemma_sub: str, plural_sub: str, stem_final: str
) -> str:
    """
    Try to collapse noisy alignment windows (e.g. 'ico→ici' vs 'co→ci')
    into a stable orthographic pattern anchored on stem_final.

    Returns a string like 'co→ci', 'che→chi', or '' if no meaningful change.
    """
    if not lemma_sub or not plural_sub or lemma_sub == plural_sub:
        return ""
    lemma_sub = lemma_sub.lower()
    plural_sub = plural_sub.lower()
    stem_final = (stem_final or "").lower()
    idx = lemma_sub.rfind(stem_final) if stem_final else -1
    if idx != -1:
        lemma_tail = lemma_sub[idx:]
        plural_tail = (
            plural_sub[idx:]
            if idx < len(plural_sub)
            else plural_sub[-len(lemma_tail) :]
        )
    else:
        base_len = max(len(stem_final), 1)
        tail_len = min(max(base_len + 1, 2), len(lemma_sub), len(plural_sub))
        lemma_tail = lemma_sub[-tail_len:]
        plural_tail = plural_sub[-tail_len:]
    if lemma_tail == plural_tail:
        return ""
    return f"{lemma_tail}→{plural_tail}"


def derive_mutation_and_orth_change(row: Dict[str, str]) -> None:
    """
    Derive mutation and orth_change fields.

    Modifies row in place, adding:
    - mutation: "True"/"False"/""
    - orth_change:
        * for true palatalization: a clean abstract pattern from
          MUTATION_PATTERNS (e.g. "c→ci", "c→ce", "st→ști", "z→ji")
        * otherwise: a canonicalized local change window from
          _canonical_orth_change(), or "" when no change.

    This keeps orth_change stable and small for real mutations (so it
    lines up with ORTH_TO_PALATAL_IPA), while still using the alignment
    window to describe other orthographic changes.
    """
    pos = row.get("pos", "")
    lemma = row.get("lemma", "")
    plural = row.get("plural", "")
    stem_final = row.get("stem_final", "")
    cluster = row.get("cluster", "")
    row["mutation"] = "False"
    row["orth_change"] = ""
    if pos not in {"N", "ADJ"} or not lemma or not plural or not stem_final:
        return
    # Use alignment to find the local change window
    lemma_sub, plural_sub = get_change_window(
        lemma, plural, stem_final, cluster
    )
    # No detectable change → explicitly no mutation
    if not lemma_sub or not plural_sub or lemma_sub == plural_sub:
        row["mutation"] = "False"
        return
    # Default: canonicalized local window (for non-palatal changes)
    row["orth_change"] = _canonical_orth_change(
        lemma_sub, plural_sub, stem_final
    )
    # If this consonant is not in our mutation inventory, we're done
    if stem_final not in MUTATION_PATTERNS:
        row["mutation"] = "False"
        return
    # Try to recognize one of the abstract palatalization patterns
    for lemma_pattern, plural_pattern in MUTATION_PATTERNS[stem_final]:
        if lemma_sub.endswith(lemma_pattern) and plural_sub.endswith(
            plural_pattern
        ):
            if plural_pattern and plural_pattern[-1] in "ie":
                row["mutation"] = "True"
                row["orth_change"] = f"{lemma_pattern}→{plural_pattern}"
                return
    # No palatalization pattern matched
    row["mutation"] = "False"


def derive_opportunity(row: Dict[str, str]) -> None:
    """
    Derive opportunity field.

    Determines if plural adds a front vowel after stem_final.

    Modifies row in place, adding:
    - opportunity: "i"/"e"/"uri"/"none"

    Args:
        row: Dictionary representing a CSV row
    """
    pos = row.get("pos", "")
    lemma = row.get("lemma", "")
    plural = row.get("plural", "")
    stem_final = row.get("stem_final", "")
    cluster = row.get("cluster", "")
    row["opportunity"] = "none"
    if pos not in {"N", "ADJ"} or not lemma or not plural or not stem_final:
        return

    lemma_l = lemma.lower()
    plural_l = plural.lower()

    # 0) Morphological neuter u/uri class: treat bona fide -uri plurals
    # as "uri" opportunity and do not mix them with i/e.
    if plural_l.endswith("uri"):
        row["opportunity"] = "uri"
        return

    # 1) PREFERRED: local plural change window (for i/e only)
    plural_sub = (row.get("plural_change_window") or "").strip()
    if plural_sub:
        for ch in plural_sub:
            if ch in ("i", "e"):
                row["opportunity"] = ch
                return
        return

    # 2) FALLBACK: canonical palatalization pattern
    orth_change = row.get("orth_change", "")
    if orth_change and orth_change in ORTH_TO_PALATAL_IPA:
        parts = orth_change.split("→", 1)
        if len(parts) == 2:
            plural_side = parts[1]
            if plural_side.endswith("i"):
                row["opportunity"] = "i"
                return
            if plural_side.endswith("e"):
                row["opportunity"] = "e"
                return

    # 3) FALLBACK: alignment-based heuristics
    aligned_lemma, aligned_plural = needleman_wunsch(lemma_l, plural_l)
    # Build index mapping from lemma index → aligned index
    lemma_to_aligned: Dict[int, int] = {}
    lemma_idx = 0
    for aligned_idx, char in enumerate(aligned_lemma):
        if char != "-":
            lemma_to_aligned[lemma_idx] = aligned_idx
            lemma_idx += 1
    # Find target end index in the lemma (orthographic stem-final position)
    if cluster:
        if cluster in VELAR_FRONT_SEQUENCES:
            target_end_idx = len(lemma_l) - 1
        else:
            stem = strip_final_vowel(lemma_l)
            target_end_idx = len(stem) - 1
    else:
        stem = strip_final_vowel(lemma_l)
        target_end_idx = -1
        for i in range(len(stem) - 1, -1, -1):
            if stem[i] == stem_final:
                target_end_idx = i
                break
    if target_end_idx == -1 or target_end_idx not in lemma_to_aligned:
        return
    aligned_target_end = lemma_to_aligned[target_end_idx]
    # Look at following characters in the aligned plural
    following_chars: List[str] = []
    for i in range(aligned_target_end + 1, len(aligned_plural)):
        if aligned_plural[i] != "-":
            following_chars.append(aligned_plural[i])
            if len(following_chars) >= 3:
                break
    if not following_chars:
        return
    first_vowel = following_chars[0]
    if first_vowel == "i":
        if aligned_target_end + 1 < len(aligned_lemma):
            if aligned_lemma[aligned_target_end + 1] != "i":
                row["opportunity"] = "i"
        else:
            row["opportunity"] = "i"
    elif first_vowel == "e":
        if aligned_target_end + 1 < len(aligned_lemma):
            if aligned_lemma[aligned_target_end + 1] != "e":
                row["opportunity"] = "e"
        else:
            row["opportunity"] = "e"


def explode_pipe_group(
    row: Dict[str, str],
    main_field: str,
    companion_fields: List[str],
    sep: str = "|",
) -> List[Dict[str, str]]:
    """
    Explode a row on `main_field` if it's pipe-separated, keeping
    `companion_fields` aligned. All other fields are duplicated.
    """

    def _split(raw_val: str) -> List[str]:
        raw_val = (raw_val or "").strip()
        if not raw_val:
            return []
        # Split on "|" and drop empty segments, normalize whitespace
        return [seg.strip() for seg in raw_val.split(sep) if seg.strip()]

    raw = row.get(main_field, "")
    items = _split(raw)
    if not items:
        # Nothing to explode; just clean spacing on companions
        for name in companion_fields:
            vals = _split(row.get(name, ""))
            row[name] = vals[0] if vals else ""
        return [row]

    def split_field(name: str) -> List[str]:
        return _split(row.get(name, ""))

    companion_lists = [split_field(name) for name in companion_fields]
    n = len(items)

    # Normalize singleton case
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
        new_row = dict(row)  # shallow copy
        new_row[main_field] = item
        for name, vals in zip(companion_fields, companion_lists):
            new_row[name] = vals[idx]
        exploded.append(new_row)
    return exploded


def explode_derived_verbs_row(row: Dict[str, str]) -> List[Dict[str, str]]:
    """
    Make the row tidyverse-friendly with respect to derived verbs.

    If a lemma has multiple derived verbs encoded as a single
    pipe-separated string, e.g.:

        derived_verbs     = "face|înfrunta"
        deriv_suffixes    = "-e|-a"
        ipa_derived_verbs = "fat͡ʃe | ɨnfrunta"

    this function returns one row per derived verb, with all other
    fields duplicated:

        [
          { ..., derived_verbs="face",    deriv_suffixes="-e",
                ipa_derived_verbs="fat͡ʃe" },
          { ..., derived_verbs="înfrunta", deriv_suffixes="-a",
                ipa_derived_verbs="ɨnfrunta" },
        ]

    If there are zero or one derived verbs, it returns a single-element
    list containing a (possibly slightly cleaned) version of the input row.
    """
    return explode_pipe_group(
        row,
        main_field="derived_verbs",
        companion_fields=["deriv_suffixes", "ipa_derived_verbs"],
    )


def explode_derived_adj_row(row: Dict[str, str]) -> List[Dict[str, str]]:
    """
    Make the row tidyverse-friendly with respect to derived adjectives.
    If a lemma has multiple derived adjectives encoded as a single
    pipe-separated string, e.g.:

        derived_adj     = "frumos|frumoasă"
        ipa_derived_adj = "frumos | frumoasă"
    """
    return explode_pipe_group(
        row,
        main_field="derived_adj",
        companion_fields=["ipa_derived_adj"],
    )


# ============================================================================
# PART 3: G2P SETUP
# ============================================================================
def normalize_unicode_g2p(s: str) -> str:
    """
    Normalize Unicode for G2P conversion.

    Args:
        s: Input string

    Returns:
        Normalized string
    """
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
    """
    Broad Romanian grapheme-to-phoneme conversion.

    Args:
        word: Romanian word in standard orthography

    Returns:
        IPA transcription (broad)
    """
    if not isinstance(word, str) or not word:
        return ""
    w = normalize_unicode_g2p(word).lower()
    for pat, repl in IPA_RULES:
        w = re.sub(pat, repl, w)
    return w


DIMINUTIVE_J_SUFFIXES = ("aică", "oaică", "uică", "eică", "iică")


def tweak_nominal_ipa(lemma: str, ipa: str) -> str:
    if not isinstance(lemma, str) or not isinstance(ipa, str):
        return ipa
    lemma_norm = normalize_unicode_g2p(lemma).lower()
    if lemma_norm.endswith(DIMINUTIVE_J_SUFFIXES):
        # Replace final i + kə with j + kə
        ipa = re.sub(r"i(kə)$", r"j\1", ipa)
    return ipa


# ============================================================================
# PART 3: CONSTANTS
# ============================================================================

ORTH_TO_PALATAL_IPA = {
    "c→ce": "t͡ʃ",
    "c→ci": "t͡ʃ",
    "ct→cți": "t͡s",
    "ct→cțe": "t͡s",
    "că→ce": "t͡ʃ",
    "că→ci": "t͡ʃ",
    "ca→ce": "t͡ʃ",  # NEW: abaca → abace
    "d→ze": "z",
    "d→zi": "z",
    "de→di": "dʲ",  # rare allophonic palatalization
    "de→zi": "z",
    "dă→zi": "z",
    "g→ge": "d͡ʒ",
    "g→gi": "d͡ʒ",
    "gă→ge": "d͡ʒ",
    "gă→gi": "d͡ʒ",
    "ga→ge": "d͡ʒ",  # NEW: ga → ge
    "go→gi": "d͡ʒ",  # NEW: flamingo → flamingi
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


# ============================================================================
# PART 3: DERIVATION FUNCTIONS
# ============================================================================
def derive_palatal_consonant_pl(row: Dict[str, str]) -> None:
    """
    Derive palatal_consonant_pl from orth_change.

    Modifies row in place, adding:
    - palatal_consonant_pl: IPA symbol or ""

    Args:
        row: Dictionary representing a CSV row
    """
    orth_change = row.get("orth_change", "")
    row["palatal_consonant_pl"] = ""
    if not orth_change:
        return
    if orth_change in ORTH_TO_PALATAL_IPA:
        row["palatal_consonant_pl"] = ORTH_TO_PALATAL_IPA[orth_change]


def derive_lemma_suffix(row: Dict[str, str]) -> None:
    """
    Derive lemma_suffix field.

    Modifies row in place, adding:
    - lemma_suffix: Suffix like "-ică" or ""

    Args:
        row: Dictionary representing a CSV row
    """

    def _normalize_suffix(suffix: str) -> str:
        """Lowercase + NFC normalization for suffix comparison."""
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
    """
    Derive target_is_suffix field.

    Modifies row in place, adding:
    - target_is_suffix: "True"/"False"

    Args:
        row: Dictionary representing a CSV row
    """
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


def derive_suffix_triggers_plural_mutation(row: Dict[str, str]) -> None:
    """
    Derive suffix_triggers_plural_mutation field.

    Modifies row in place, adding:
    - suffix_triggers_plural_mutation: "True"/"False"

    Args:
        row: Dictionary representing a CSV row
    """
    pos = row.get("pos", "")
    lemma_suffix = row.get("lemma_suffix", "")
    opportunity = row.get("opportunity", "")
    target_is_suffix = row.get("target_is_suffix", "")
    mutation = row.get("mutation", "")
    # Default: suffix is *not* triggering mutation for this lemma
    row["suffix_triggers_plural_mutation"] = "False"
    if (
        pos == "N"
        and opportunity in {"i", "e"}
        and mutation == "True"
        and target_is_suffix == "True"
        and lemma_suffix in {"-ic", "-ist", "-esc", "-ică"}
    ):
        row["suffix_triggers_plural_mutation"] = "True"


# ============================================================================
# PART 4: DERIVATION FUNCTIONS
# ============================================================================
def derive_derived_verbs_fields(row: Dict[str, str]) -> None:
    """
    Derive deriv_suffixes and ipa_derived_verbs.

    Keeps only denominal verbs ending in -a, -i, or -ui.
    Drops other verbs (out of scope for the project).

    Modifies row in place, adding:
    - derived_verbs: cleaned pipe-separated verbs
    - deriv_suffixes: pipe-separated suffixes (-a / -i / -ui)
    - ipa_derived_verbs: pipe-separated normalized IPA
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
        # Normalize orthography for suffix detection
        norm = normalize_unicode_g2p(verb).lower()
        if norm.endswith("ui"):
            suf = "-ui"
        elif norm.endswith("i"):
            suf = "-i"
        elif norm.endswith("a"):
            suf = "-a"
        else:
            # Out-of-domain verb (e.g. ends in -e): drop it
            continue

        clean_verbs.append(verb)
        suffixes.append(suf)
        # Generate IPA, then normalize so affricates get tie bars, etc.
        ipa = to_ipa(verb)
        if _ipa_normalizer:
            ipa = _ipa_normalizer(ipa)
        ipa_list.append(ipa)
    if not clean_verbs:
        # Treat as "no derived verbs" for this project
        row["derived_verbs"] = ""
        row["deriv_suffixes"] = ""
        row["ipa_derived_verbs"] = ""
        return

    row["derived_verbs"] = "|".join(clean_verbs)
    row["deriv_suffixes"] = "|".join(suffixes)
    row["ipa_derived_verbs"] = "|".join(ipa_list)


def derive_derived_adj_fields(row: Dict[str, str]) -> None:
    """
    Derive ipa_derived_adj.

    Modifies row in place, adding:
    - ipa_derived_adj: pipe-separated normalized IPA for each derived adjective
    """
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

    # Normalize spacing in derived_adj itself
    row["derived_adj"] = "|".join(adj_list)
    row["ipa_derived_adj"] = "|".join(ipa_list)


def derive_nde_class(row: Dict[str, str]) -> None:
    """
    Classify NDE patterns.

    Types:
    1. gimpe: Tautomorphemic C+front-vowel
    2. ochi: Singular=plural with chi/ghi
    3. paduchi: che/ghe→chi/ghi under-application
    4. frontstem: ce/ci / ge/gi where plural only toggles e↔i

    Modifies row in place, adding:
    - nde_class: "gimpe"/"ochi"/"paduchi"/"frontstem"/""

    Args:
        row: Dictionary representing a CSV row
    """
    pos = row.get("pos", "")
    lemma = row.get("lemma", "") or ""
    plural = row.get("plural", "") or ""
    stem_final = row.get("stem_final", "") or ""
    cluster = row.get("cluster", "") or ""
    mutation = str(row.get("mutation", ""))
    orth_change = row.get("orth_change", "") or ""
    row["nde_class"] = ""
    # NDEs are non-mutating dorsal nouns with both forms present
    if pos != "N" or not lemma or not plural:
        return

    if mutation == "True":
        return

    if stem_final not in ("c", "g"):
        return

    # (1) ochi-type: sg = pl, chi/ghi
    if lemma == plural and cluster in ("chi", "ghi"):
        row["nde_class"] = "ochi"
        return

    # (2) paduchi-type: che/ghe → chi/ghi (optionally allow chiuri/ghiuri)
    if cluster in ("che", "ghe"):
        expected_cluster_pl = cluster[:-1] + "i"
        if plural.endswith(expected_cluster_pl) or plural.endswith(
            expected_cluster_pl + "uri"
        ):
            row["nde_class"] = "paduchi"
            return

    # (3) gimpe-type: sg = pl, final ci/ce/gi/ge, no velar-front cluster
    if not cluster:
        if lemma == plural and lemma.endswith(("ci", "ce", "gi", "ge")):
            row["nde_class"] = "gimpe"
            return

        # (4) front-stem ce/ci/ge/gi where plural only toggles e↔i
        if orth_change in {"ce→ci", "ci→ce", "ge→gi", "gi→ge"}:
            row["nde_class"] = "frontstem"
            return


def derive_exception_reason(row: Dict[str, str]) -> None:
    """
    Derive exception_reason field.

    Modifies row in place, adding:
    - exception_reason: "nde:gimpe"/"nde:ochi"/"nde:paduchi"/""

    Args:
        row: Dictionary representing a CSV row
    """
    row["exception_reason"] = ""
    pos = row.get("pos", "")
    plural = row.get("plural", "")
    mutation = row.get("mutation", "")
    nde = row.get("nde_class", "")
    # Only nouns with an attested plural
    if pos != "N" or not plural or mutation == "True":
        return
    if nde in {"gimpe", "ochi", "paduchi"}:
        row["exception_reason"] = f"nde:{nde}"


def derive_is_true_exception(row: Dict[str, str]) -> None:
    """
    Derive is_true_exception field.

    A noun is a "true exception" if it has an i/e opportunity,
    does not palatalize, and is not an NDE pattern.

    Modifies row in place, adding:
    - is_true_exception: "True"/"False"

    Args:
        row: Dictionary representing a CSV row
    """
    pos = row.get("pos", "")
    mutation = row.get("mutation", "")
    opportunity = row.get("opportunity", "")
    nde_class = row.get("nde_class", "")
    has_ie_opportunity = opportunity in {"i", "e"}
    # Treat any classified NDE pattern as NDE, including "frontstem"
    is_nde = bool(nde_class)
    if (
        pos == "N"
        and has_ie_opportunity
        and mutation == "False"
        and not is_nde
    ):
        row["is_true_exception"] = "True"
    else:
        row["is_true_exception"] = "False"
