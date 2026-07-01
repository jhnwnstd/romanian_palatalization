"""Improved mismatch categorisation.

The driver's original categoriser tagged every mismatch that didn't
fit the four known lexical exception classes as UNKNOWN or
IPA_SEGMENT_MISMATCH — burying real rule failures under a mountain of
transcription-noise cases. This module peels off the noise first, so
the residual RULE_FAILURE bucket is small enough to hand-audit.

Two-pass structure:

1. **Normalisation pass** — try successive normalisers (trailing
   glide, unstressed a↔ə, medial j↔i, ea/oa monophthongisation).
   The first normalisation that makes predicted equal any attested
   variant is recorded as the mismatch's category. Rules were right;
   the mismatch is transcription noise.

2. **Residual pass** — for whatever survives normalisation, apply
   the original lexical-exception classifiers (EZ_ETHNONYM,
   ICA_SUPPLETION, SHCA_CLUSTER_DATA, SHT_LOANWORD) then fall through
   to RULE_FAILURE.

The output CSV gains a ``normalisations`` column so a downstream
reader can filter out transcription-noise rows and focus on the
residual.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Final, Protocol

from .compare import CompareResult, Normalisation, compare


class RefinedCategory(StrEnum):
    """Extended mismatch category vocabulary.

    Overlaps with the driver's :class:`MismatchCategory` but with
    finer distinctions on the transcription-noise side. New entries
    are prefixed with the normalisation that would fix them.
    """

    # Data-side transcription noise (rules RIGHT):
    NORM_TRAILING_GLIDE = "NORM_TRAILING_GLIDE"      # j/ʲ vs i
    NORM_UNSTRESSED_SCHWA = "NORM_UNSTRESSED_SCHWA"  # a vs ə
    NORM_MEDIAL_GLIDE = "NORM_MEDIAL_GLIDE"          # medial j vs i
    NORM_EA_MONOPHTHONG = "NORM_EA_MONOPHTHONG"      # ea vs e
    NORM_OA_MONOPHTHONG = "NORM_OA_MONOPHTHONG"      # oa vs o
    NORM_STRIP_MARKS = "NORM_STRIP_MARKS"            # stress/length marks

    # Known lexical exception classes (rules can't handle):
    EZ_ETHNONYM = "EZ_ETHNONYM"
    ICA_SUPPLETION = "ICA_SUPPLETION"
    SHCA_CLUSTER_DATA = "SHCA_CLUSTER_DATA"
    SHT_LOANWORD = "SHT_LOANWORD"

    # Residual — likely a genuine rule failure worth inspection:
    RULE_FAILURE = "RULE_FAILURE"

    # UR construction failed (bad IPA in the lexicon):
    UR_BUILD_FAILED = "UR_BUILD_FAILED"


_NORM_TO_CATEGORY: Final[dict[Normalisation, RefinedCategory]] = {
    Normalisation.TRAILING_GLIDE: RefinedCategory.NORM_TRAILING_GLIDE,
    Normalisation.UNSTRESSED_SCHWA: RefinedCategory.NORM_UNSTRESSED_SCHWA,
    Normalisation.MEDIAL_GLIDE: RefinedCategory.NORM_MEDIAL_GLIDE,
    Normalisation.EA_MONOPHTHONG: RefinedCategory.NORM_EA_MONOPHTHONG,
    Normalisation.OA_MONOPHTHONG: RefinedCategory.NORM_OA_MONOPHTHONG,
    Normalisation.STRIP_STRESS: RefinedCategory.NORM_STRIP_MARKS,
}


class _RowLike(Protocol):
    """Minimal shape the categoriser needs from a lexicon row."""
    lemma: str
    plural: str
    stem_final: str
    opportunity: str
    mutation: bool


@dataclass(frozen=True, slots=True)
class Categorisation:
    """A categorised mismatch."""

    category: RefinedCategory
    normalisations: tuple[Normalisation, ...]
    edit_distance: int
    attested_variant: str


def refine_category(
    row: _RowLike,
    predicted_sr: str,
    attested_field: str,
    boolean_predicted: bool,
    ur_built: bool,
) -> Categorisation:
    """Assign a :class:`RefinedCategory` to one mismatch row.

    The caller pre-computes whether the boolean prediction is right —
    we don't re-derive here; we just categorise.
    """
    if not ur_built:
        return Categorisation(
            RefinedCategory.UR_BUILD_FAILED, (), 999, ""
        )

    cmp: CompareResult = compare(predicted_sr, attested_field)
    if cmp.matched and cmp.normalisation_applied:
        # Data-side normalisation succeeded — rules were right,
        # transcription differs.
        primary_norm = cmp.normalisation_applied[0]
        cat = _NORM_TO_CATEGORY.get(
            primary_norm, RefinedCategory.RULE_FAILURE,
        )
        return Categorisation(
            cat, cmp.normalisation_applied, 0, cmp.attested_variant,
        )
    if cmp.matched:
        # No normalisation needed — the caller shouldn't be calling
        # us on this row, but categorise defensively.
        return Categorisation(
            RefinedCategory.RULE_FAILURE, (),
            cmp.edit_distance, cmp.attested_variant,
        )

    # Not matched under any normalisation — apply lexical-exception
    # classifiers.
    lem = row.lemma.lower()
    sf = row.stem_final
    if (
        boolean_predicted and not row.mutation
        and sf == "z" and row.opportunity == "i"
        and (lem.endswith("ez") or lem.endswith("eză"))
    ):
        return Categorisation(
            RefinedCategory.EZ_ETHNONYM, (),
            cmp.edit_distance, cmp.attested_variant,
        )
    if (
        boolean_predicted and not row.mutation
        and sf == "c" and row.opportunity == "e"
        and lem.endswith("ică")
        and row.plural.lower().endswith("ele")
    ):
        return Categorisation(
            RefinedCategory.ICA_SUPPLETION, (),
            cmp.edit_distance, cmp.attested_variant,
        )
    if (
        boolean_predicted and not row.mutation
        and sf == "c" and row.opportunity == "i"
        and lem.endswith(("șcă", "şcă"))
    ):
        return Categorisation(
            RefinedCategory.SHCA_CLUSTER_DATA, (),
            cmp.edit_distance, cmp.attested_variant,
        )
    if (
        boolean_predicted and not row.mutation
        and sf == "t" and row.opportunity == "i"
        and lem.endswith(("șt", "şt"))
    ):
        return Categorisation(
            RefinedCategory.SHT_LOANWORD, (),
            cmp.edit_distance, cmp.attested_variant,
        )

    return Categorisation(
        RefinedCategory.RULE_FAILURE, (),
        cmp.edit_distance, cmp.attested_variant,
    )


__all__: Final[tuple[str, ...]] = (
    "Categorisation",
    "RefinedCategory",
    "refine_category",
)
