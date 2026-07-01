"""Canonical predicted-vs-attested IPA comparison stack.

The diagnostic modules and the driver's mismatch categoriser all
need to answer "does this predicted SR match the attested plural
IPA?" — but "match" has more than one legitimate meaning in the
lexicon:

  - the same string byte-for-byte,
  - the same string modulo a trailing /j/↔/i/ desyllabification
    convention (the paper writes [ʲ], Wiktionary writes /i/),
  - the same string after normalising unstressed a↔ə, medial j↔i,
    ea→e before palatals, etc.

Routing every check through one function fixes two problems: a
single place to add new normalisations, and consistent "success"
semantics across ordering search, perturbation search, and the
mismatch report. If ordering search calls the SR "matching" but
the driver's CSV calls it "IPA_SEGMENT_MISMATCH", the two
disagree — bugs like that are exactly what this module prevents.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Final


class Normalisation(StrEnum):
    """A specific normalisation applied to make two IPAs equal.

    Reported so the caller can tell WHICH transcription convention
    was smoothed over — a raw-vs-normalised diff is diagnostic.
    """

    NONE = "NONE"
    TRAILING_GLIDE = "TRAILING_GLIDE"     # j/ʲ ↔ i at right edge
    UNSTRESSED_SCHWA = "UNSTRESSED_SCHWA"  # a ↔ ə internal
    MEDIAL_GLIDE = "MEDIAL_GLIDE"          # j ↔ i between vowels
    EA_MONOPHTHONG = "EA_MONOPHTHONG"      # ea ↔ e before palatal
    OA_MONOPHTHONG = "OA_MONOPHTHONG"      # oa ↔ o before palatal
    STRIP_STRESS = "STRIP_STRESS"          # ˈ ˌ ' ` removed
    STRIP_LENGTH = "STRIP_LENGTH"          # ː removed


@dataclass(frozen=True, slots=True)
class CompareResult:
    """Outcome of a predicted-vs-attested comparison."""

    matched: bool
    attested_variant: str                    # which pipe-split variant matched (or "")
    normalisation_applied: tuple[Normalisation, ...] = ()
    edit_distance: int = 0                   # Levenshtein on the closest attested variant


# ---------------------------------------------------------------------------
# Normalisers — each returns a possibly-shorter/simpler string.
# ---------------------------------------------------------------------------


def _strip_marks(s: str) -> str:
    """Remove stress and length marks that shouldn't affect matching."""
    for ch in ("ˈ", "ˌ", "'", "`", "ː"):
        s = s.replace(ch, "")
    return s


def _normalise_trailing_glide(s: str) -> str:
    """Trailing /j/ or /ʲ/ → /i/.

    The paper writes desyllabified plural /i/ as [ʲ] (palatal mark
    on the preceding consonant); Wiktionary writes it as the vowel
    /i/; our pipeline outputs /j/. All three describe the same
    surface fact — normalise to /i/ on both sides.
    """
    while s and s[-1] in ("j", "ʲ"):
        s = s[:-1] + "i"
        break   # only touch the final segment
    return s


def _normalise_unstressed_a(pred: str, attested: str) -> str:
    """Fold predicted ``a``↔``ə`` to match attested at same positions.

    Romanian unstressed /a/ often surfaces as [ə]. The lexicon
    sometimes records the reduction in the plural but not the lemma,
    or vice versa (variant picker might have chosen a schwa variant
    even when the plural has [a]). Fold in either direction so both
    conventions unify.
    """
    if len(pred) != len(attested):
        return pred
    out: list[str] = []
    for p, a in zip(pred, attested):
        if p == "a" and a == "ə":
            out.append("ə")
        elif p == "ə" and a == "a":
            out.append("a")
        else:
            out.append(p)
    return "".join(out)


def _normalise_medial_glide(pred: str, attested: str) -> str:
    """Fold /j/↔/i/ in ``pred`` to match attested at same positions."""
    if len(pred) != len(attested):
        return pred
    out: list[str] = []
    for p, a in zip(pred, attested):
        if p == "j" and a == "i":
            out.append("i")
        elif p == "i" and a == "j":
            out.append("j")
        else:
            out.append(p)
    return "".join(out)


def _normalise_ea_monophthong(pred: str, attested: str) -> str:
    """Fold /ea/ in pred to /e/ where attested has /e/.

    Romanian monophthongises the ``ea`` diphthong to /e/ before
    palatal consonants (``-ancă → -ence``): pred keeps the diphthong
    from the lemma but attested has the monophthongised plural form.
    """
    if "ea" not in pred:
        return pred
    # Alignment-free: try replacing each 'ea' occurrence in pred
    # with 'e' one at a time and return whichever produces attested.
    for i in range(len(pred) - 1):
        if pred[i:i+2] == "ea":
            candidate = pred[:i] + "e" + pred[i+2:]
            if candidate == attested:
                return candidate
    return pred


def _normalise_oa_monophthong(pred: str, attested: str) -> str:
    """Fold /oa/ in pred to /o/ where attested has /o/."""
    if "oa" not in pred:
        return pred
    for i in range(len(pred) - 1):
        if pred[i:i+2] == "oa":
            candidate = pred[:i] + "o" + pred[i+2:]
            if candidate == attested:
                return candidate
    return pred


# ---------------------------------------------------------------------------
# Levenshtein — pure Python fallback so we don't add a hard dep.
# ---------------------------------------------------------------------------


def _levenshtein(a: str, b: str) -> int:
    """Standard iterative edit distance. O(len(a) * len(b))."""
    if a == b:
        return 0
    if not a:
        return len(b)
    if not b:
        return len(a)
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        curr = [i] + [0] * len(b)
        for j, cb in enumerate(b, 1):
            cost = 0 if ca == cb else 1
            curr[j] = min(
                prev[j] + 1,
                curr[j - 1] + 1,
                prev[j - 1] + cost,
            )
        prev = curr
    return prev[-1]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def split_variants(attested_field: str) -> list[str]:
    """Split a pipe-separated IPA field, filter out foreign-language
    annotations, and return the cleaned variants."""
    if not attested_field:
        return []
    parts = [p.strip() for p in attested_field.split(" | ") if p.strip()]
    return [
        p for p in parts
        if not (len(p) >= 3 and p[2] == ":" and p[:2].isalpha()
                and p[:2].islower())
    ]


def compare(predicted: str, attested_field: str) -> CompareResult:
    """Compare ``predicted`` (single IPA string) to a pipe-separated
    ``attested_field`` under successive normalisations.

    Returns the first successful match, the normalisations applied to
    get there, and the minimum edit distance to any attested variant
    (which lets the caller sort mismatches by closeness even when no
    normalisation succeeds).
    """
    if not predicted or not attested_field:
        return CompareResult(matched=False, attested_variant="")

    variants = split_variants(attested_field)
    if not variants:
        return CompareResult(matched=False, attested_variant="")

    # Fast path: exact byte-for-byte match against any variant.
    for v in variants:
        if predicted == v:
            return CompareResult(
                matched=True, attested_variant=v,
                normalisation_applied=(),
                edit_distance=0,
            )

    # Sequential normalisers. We try each layer against every variant
    # and return the first hit; if nothing hits, we track the best
    # (lowest-distance) approximation for the caller.
    pred_stripped = _strip_marks(predicted)
    best_distance = min(_levenshtein(predicted, v) for v in variants)
    best_variant = min(variants, key=lambda v: _levenshtein(predicted, v))

    normalisers: list[tuple[Normalisation, callable]] = [
        (Normalisation.STRIP_STRESS,
            lambda p, a: _strip_marks(p)),
        (Normalisation.TRAILING_GLIDE,
            lambda p, a: _normalise_trailing_glide(p)),
        (Normalisation.UNSTRESSED_SCHWA,
            _normalise_unstressed_a),
        (Normalisation.MEDIAL_GLIDE,
            _normalise_medial_glide),
        (Normalisation.EA_MONOPHTHONG,
            _normalise_ea_monophthong),
        (Normalisation.OA_MONOPHTHONG,
            _normalise_oa_monophthong),
    ]

    for v in variants:
        applied: list[Normalisation] = []
        cur = pred_stripped
        v_stripped = _strip_marks(v)
        # Iteratively apply normalisers; stop early if we hit.
        for label, fn in normalisers:
            if cur == v_stripped:
                return CompareResult(
                    matched=True, attested_variant=v,
                    normalisation_applied=tuple(applied),
                    edit_distance=0,
                )
            new = fn(cur, v_stripped)
            if new != cur:
                applied.append(label)
                cur = new
        if cur == v_stripped:
            return CompareResult(
                matched=True, attested_variant=v,
                normalisation_applied=tuple(applied),
                edit_distance=0,
            )

    return CompareResult(
        matched=False,
        attested_variant=best_variant,
        normalisation_applied=(),
        edit_distance=best_distance,
    )


__all__: Final[tuple[str, ...]] = (
    "CompareResult",
    "Normalisation",
    "compare",
    "split_variants",
)
