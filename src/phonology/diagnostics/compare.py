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

import math
from dataclasses import dataclass
from enum import StrEnum
from typing import Callable, Final, Sequence


class Normalisation(StrEnum):
    """A specific normalisation applied to make two IPAs equal.

    Reported so the caller can tell WHICH transcription convention
    was smoothed over — a raw-vs-normalised diff is diagnostic.
    """

    NONE = "NONE"
    TRAILING_GLIDE = "TRAILING_GLIDE"      # j/ʲ ↔ i at right edge
    UNSTRESSED_SCHWA = "UNSTRESSED_SCHWA"  # a ↔ ə internal
    MEDIAL_GLIDE = "MEDIAL_GLIDE"          # j ↔ i between vowels
    EA_MONOPHTHONG = "EA_MONOPHTHONG"      # ea ↔ e before palatal
    OA_MONOPHTHONG = "OA_MONOPHTHONG"      # oa ↔ o before palatal
    STRIP_STRESS = "STRIP_STRESS"          # ˈ ˌ ' ` and ː removed


class CompareStatus(StrEnum):
    """Terminal state of a compare_ipa() call."""

    MATCH = "MATCH"
    MISMATCH = "MISMATCH"
    EMPTY_ATTESTED = "EMPTY_ATTESTED"      # no variants supplied
    EMPTY_PREDICTED = "EMPTY_PREDICTED"    # nothing to score


@dataclass(frozen=True, slots=True)
class CompareResult:
    """Outcome of a predicted-vs-attested comparison.

    ``edit_distance`` is ``math.inf`` when comparison couldn't happen
    (empty predicted or empty attested set); consumers ranking by
    distance MUST check ``status`` first rather than treating a raw
    integer 0 as "close" — this was a silent-zero bug in an earlier
    version of the function.
    """

    matched: bool
    status: CompareStatus
    attested_variant: str
    normalisation_applied: tuple[Normalisation, ...] = ()
    edit_distance: float = 0.0


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


def _fold_diphthong_bulk(pred: str, attested: str, diphthong: str, mono: str) -> str:
    """Collapse every ``diphthong`` in ``pred`` to ``mono`` if that
    equals ``attested``. Handles the single-occurrence and multi-
    occurrence cases uniformly (the old per-position loop only tried
    one substitution at a time and failed on strings with two
    ``ea``s).
    """
    if diphthong not in pred:
        return pred
    fully_collapsed = pred.replace(diphthong, mono)
    if fully_collapsed == attested:
        return fully_collapsed
    # Also try single-position substitutions in case the caller
    # actually wants to collapse only one of several occurrences.
    for i in range(len(pred) - len(diphthong) + 1):
        if pred[i:i + len(diphthong)] == diphthong:
            candidate = pred[:i] + mono + pred[i + len(diphthong):]
            if candidate == attested:
                return candidate
    return pred


def _normalise_ea_monophthong(pred: str, attested: str) -> str:
    """Fold /ea/ in pred to /e/ where attested has /e/.

    Romanian monophthongises the ``ea`` diphthong to /e/ before
    palatal consonants (``-ancă → -ence``): pred keeps the diphthong
    from the lemma but attested has the monophthongised plural form.
    Handles multiple ``ea`` occurrences via a bulk substitution
    followed by per-position fallbacks.
    """
    return _fold_diphthong_bulk(pred, attested, "ea", "e")


def _normalise_oa_monophthong(pred: str, attested: str) -> str:
    """Fold /oa/ in pred to /o/ where attested has /o/."""
    return _fold_diphthong_bulk(pred, attested, "oa", "o")


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


def is_foreign_annotation(text: str) -> bool:
    """True iff ``text`` starts with a 2-letter lowercase language
    code and a colon (e.g. ``fr:absorbant``, ``tr:emme``).

    Extracted so the driver, compare, and g2p share one predicate
    instead of three lookalike ``p[2] == ':'`` checks that drifted.
    """
    return (
        len(text) >= 3
        and text[2] == ":"
        and text[:2].isalpha()
        and text[:2].islower()
    )


def split_variants(attested_field: str) -> tuple[str, ...]:
    """Split a pipe-separated IPA field, drop foreign-language
    annotations, and return the cleaned variants.

    Returns a tuple (not a list) so callers can rely on immutability
    when passing the result around. Empty input or all-foreign input
    yields an empty tuple — callers should check length rather than
    relying on ``[0]``.
    """
    if not attested_field:
        return ()
    parts = (p.strip() for p in attested_field.split(" | "))
    return tuple(p for p in parts if p and not is_foreign_annotation(p))


# The normaliser cascade — a module-level constant so we don't rebuild
# the list per compare() call. Each entry is (label, fn) where fn takes
# (pred, attested) and returns a possibly-transformed pred. Entries are
# applied in order; the first arrangement that makes pred == attested
# wins.
_NORMALISERS: Final[tuple[
    tuple[Normalisation, Callable[[str, str], str]], ...
]] = (
    (Normalisation.STRIP_STRESS, lambda p, a: _strip_marks(p)),
    (Normalisation.TRAILING_GLIDE,
        lambda p, a: _normalise_trailing_glide(p)),
    (Normalisation.UNSTRESSED_SCHWA, _normalise_unstressed_a),
    (Normalisation.MEDIAL_GLIDE, _normalise_medial_glide),
    (Normalisation.EA_MONOPHTHONG, _normalise_ea_monophthong),
    (Normalisation.OA_MONOPHTHONG, _normalise_oa_monophthong),
)


def compare_ipa(
    predicted: str,
    attested: str | Sequence[str],
    *,
    strict: bool = False,
) -> CompareResult:
    """Compare ``predicted`` to ``attested`` under successive
    normalisations.

    ``attested`` accepts either a raw pipe-separated IPA field (as
    stored in the lexicon CSV) OR an already-split tuple/list of
    individual variants. The former is convenient at the CSV boundary;
    the latter is more ergonomic when the caller already has cleaned
    variants.

    Returns the first successful match, the normalisations applied to
    get there, and (on a genuine mismatch) the minimum edit distance
    to any attested variant. When ``attested`` is empty or all-foreign
    the result's ``status`` is ``EMPTY_ATTESTED`` and ``edit_distance``
    is ``inf`` so distance-ranking callers cannot silently treat "no
    data" as "matches at distance 0".

    Set ``strict=True`` to skip the normalisation cascade — only
    byte-for-byte matches count.
    """
    if not predicted:
        return CompareResult(
            matched=False, status=CompareStatus.EMPTY_PREDICTED,
            attested_variant="", edit_distance=math.inf,
        )

    if isinstance(attested, str):
        variants = split_variants(attested)
    else:
        variants = tuple(v for v in attested if v)
    if not variants:
        return CompareResult(
            matched=False, status=CompareStatus.EMPTY_ATTESTED,
            attested_variant="", edit_distance=math.inf,
        )

    # Fast path: exact byte-for-byte match against any variant.
    for v in variants:
        if predicted == v:
            return CompareResult(
                matched=True, status=CompareStatus.MATCH,
                attested_variant=v, edit_distance=0.0,
            )

    if strict:
        best_variant = min(
            variants, key=lambda v: _levenshtein(predicted, v),
        )
        return CompareResult(
            matched=False, status=CompareStatus.MISMATCH,
            attested_variant=best_variant,
            edit_distance=float(_levenshtein(predicted, best_variant)),
        )

    pred_stripped = _strip_marks(predicted)

    for v in variants:
        v_stripped = _strip_marks(v)
        applied: list[Normalisation] = []
        cur = pred_stripped
        for label, fn in _NORMALISERS:
            if cur == v_stripped:
                return CompareResult(
                    matched=True, status=CompareStatus.MATCH,
                    attested_variant=v,
                    normalisation_applied=tuple(applied),
                    edit_distance=0.0,
                )
            new = fn(cur, v_stripped)
            if new != cur:
                applied.append(label)
                cur = new
        if cur == v_stripped:
            return CompareResult(
                matched=True, status=CompareStatus.MATCH,
                attested_variant=v,
                normalisation_applied=tuple(applied),
                edit_distance=0.0,
            )

    best_variant = min(variants, key=lambda v: _levenshtein(predicted, v))
    return CompareResult(
        matched=False, status=CompareStatus.MISMATCH,
        attested_variant=best_variant,
        edit_distance=float(_levenshtein(predicted, best_variant)),
    )


def compare(predicted: str, attested_field: str) -> CompareResult:
    """Back-compat alias for the pipe-separated-string overload of
    :func:`compare_ipa`. New callers should use ``compare_ipa`` and
    pass either the raw field or a variant tuple as appropriate.
    """
    return compare_ipa(predicted, attested_field)


__all__: Final[tuple[str, ...]] = (
    "CompareResult",
    "CompareStatus",
    "Normalisation",
    "compare",
    "compare_ipa",
    "is_foreign_annotation",
    "split_variants",
)
