"""Romanian palatalization analysis — the paper as data.

This module instantiates the framework for the analysis in
``latex.tex``. Editing the analysis is editing this one file; the rest
of ``src/phonology/`` knows nothing about Romanian.

Three things live here:

1. **Patches** to ``romanian_features.json`` that the analysis
   requires. The paper makes claims about /i/ and /j/ being
   postalveolar-coronal, /t/ and /d/ being prespecified -Strident,
   and the affricates' DelRel being underspecified ("assumed for
   free"). Rather than edit the JSON, we declare these patches here
   so the inventory stays auditable.

2. **Underspecified segments** /K, G, S, Z, T, D/ — the inalterable-
   by-prespecification mechanism. Each variant is a clone of a
   surface segment with selected features re-opened to ``"0"``.

3. **The ordered rule tuple** ``RULES_FROM_PAPER`` — the
   palatalization pipeline written rule-by-rule, with cross-
   references to the equation labels in ``latex.tex``.

   Spell-out is not in the pipeline. The validator builds the UR
   with the attested plural's suffix pre-attached (see
   ``phonology.g2p``), which avoids circularity: testing the
   palatalization rules given the morphology, separately from
   testing the morphology itself.

The ``/K/, /G/`` underspecifications clear ``CORONAL`` and ``DORSAL``
but *not* ``Continuant``: dorsal palatalization unifies the
postalveolar place bundle into the segment, but no -Continuant is
supplied, so the segment keeps the ``-Continuant`` it inherited from
``/k/`` or ``/ɡ/`` and projects to the affricate ``/t͡ʃ/`` (resp.
``/d͡ʒ/``), not to the fricative ``/ʃ/``. Clearing ``Continuant``
would let the coronal default ``+Continuant`` fill-in fire and
spirantize the result — exactly the wrong outcome.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Callable, Final, Mapping, cast

from ..inventory import FeatureInventory, FeaturePatch, UnderspecifiedSegment
from ..lexicon import ClusterTag, LexRow
from ..lp import NaturalClass, natural_class
from ..rules import (
    DeletionRule,
    GlideFormationRule,
    Rule,
    UnificationRule,
)
from ..search import Direction, Search
from ..segments import Segment


def _frozen(d: dict[str, str]) -> Mapping[str, str]:
    """Wrap a feature-pattern dict as a read-only mapping so multiple
    rules sharing the same shorthand can't be corrupted by any code
    path that mutates the dict in place."""
    return MappingProxyType(d)


# ---------------------------------------------------------------------------
# Patches to the base inventory
# ---------------------------------------------------------------------------

PATCHES: Final[tuple[FeaturePatch, ...]] = (
    FeaturePatch(
        segment="i",
        updates={"CORONAL": "+", "Anterior": "-", "Distributed": "+"},
        reason=(
            "Paper (latex.tex:529-535): a high front vocoid is "
            "coronal, so /i/ is the postalveolar-place terminator for "
            "S-palatalization and assibilation in the plain -i forms."
        ),
    ),
    FeaturePatch(
        segment="j",
        updates={"CORONAL": "+", "Anterior": "-", "Distributed": "+"},
        reason=(
            "Multiply-articulated /j/: the desyllabified vocoid keeps "
            "the place class of /i/ so it continues to act as a "
            "trigger after glide formation."
        ),
    ),
    FeaturePatch(
        segment="t",
        updates={"Strident": "-"},
        reason=(
            "Prespecified /t/: -Strident is the inalterability "
            "condition for assibilation (latex.tex:339)."
        ),
    ),
    FeaturePatch(
        segment="d",
        updates={"Strident": "-"},
        reason=(
            "Prespecified /d/: -Strident matches /t/'s prespecification "
            "and blocks assibilation on lexical /d/."
        ),
    ),
    FeaturePatch(
        segment="t͡ʃ",
        updates={"DelRel": "0"},
        reason=(
            "Paper: 'we assume DelRel for free' — the affricate "
            "marker is left underspecified so dorsal palatalization "
            "needn't supply it."
        ),
    ),
    FeaturePatch(
        segment="d͡ʒ",
        updates={"DelRel": "0"},
        reason="See /t͡ʃ/ patch: DelRel is assumed free.",
    ),
    FeaturePatch(
        segment="t͡s",
        updates={"DelRel": "0"},
        reason="See /t͡ʃ/ patch: DelRel is assumed free.",
    ),
    # Note: no /d͡z/ patch — Romanian has no [d͡z] in the inventory.
    # This gap is exactly what the analysis predicts (latex.tex:367).
)


# ---------------------------------------------------------------------------
# Underspecified segments
# ---------------------------------------------------------------------------

UNDERSPEC: Final[tuple[UnderspecifiedSegment, ...]] = (
    UnderspecifiedSegment(
        label="K",
        base="k",
        # Closure of every feature dorsal-palatalization writes: if the
        # base /k/ ever gains a concrete Anterior/Distributed/Strident
        # value in the JSON (currently all 0), the analysis would
        # silently break because unification would hit a contrary
        # explicit value. Clearing the full set makes /K/ auditably
        # inalterable-free w.r.t. every dorsal-pal supply feature.
        # Continuant is *not* cleared: it stays -Continuant (from /k/)
        # so the palatalised result is the affricate /t͡ʃ/, not /ʃ/.
        clear=frozenset(
            {"CORONAL", "DORSAL", "Anterior", "Distributed", "Strident"}
        ),
        reason=(
            "/K/ is the palatalisable velar. Clear set = closure of "
            "dorsal-pal's supply features (minus Continuant, which "
            "must stay -Cont to produce the affricate, not the "
            "fricative)."
        ),
    ),
    UnderspecifiedSegment(
        label="G",
        base="ɡ",
        clear=frozenset(
            {"CORONAL", "DORSAL", "Anterior", "Distributed", "Strident"}
        ),
        reason=(
            "Voiced counterpart of /K/. Same closure clear set; "
            "voice is inherited from /ɡ/."
        ),
    ),
    UnderspecifiedSegment(
        label="S",
        base="s",
        clear=frozenset({"Anterior", "Distributed"}),
        reason=(
            "/S/ is the palatalisable strident coronal fricative. "
            "Place features cleared so S-pal can write the postalv "
            "bundle; Strident and Continuant kept (inherited from /s/)."
        ),
    ),
    # Note: no /Z/. The paper explicitly says "fully specified /z/ is
    # inalterable in all contexts" (latex.tex:170), so we do NOT
    # declare an underspecified twin. The rare z-final stems that
    # palatalize (e.g. some hypocoristic verbs) surface as
    # UNKNOWN exceptions and are not covered by the analysis.
    UnderspecifiedSegment(
        label="T",
        base="t",
        clear=frozenset({"Strident"}),
        reason=(
            "/T/ is the assibilable coronal stop. Strident is cleared "
            "so assibilation can write +Strident; Continuant is *not* "
            "cleared. The kept -Continuant is what makes assibilated "
            "/T/ surface as the affricate [t͡s] rather than the "
            "fricative [s] — coronal-default's +Continuant supply "
            "conflicts with /T/'s explicit -Continuant and fails."
        ),
    ),
    UnderspecifiedSegment(
        label="D",
        base="d",
        clear=frozenset({"Strident", "Continuant"}),
        reason=(
            "/D/ asymmetrically clears Continuant in addition to "
            "Strident. After assibilation supplies +Strident, the "
            "coronal-default's +Continuant supply finds /D/'s "
            "Continuant open and fills it — so /D/ spirantizes to "
            "[z]. This is why Romanian has no [d͡z]: the only "
            "segment that could surface as [d͡z] is exactly the one "
            "the default spirantizes (latex.tex:366-367)."
        ),
    ),
)


# ---------------------------------------------------------------------------
# Inventory factory
# ---------------------------------------------------------------------------

_INVENTORY_JSON: Final[Path] = (
    Path(__file__).resolve().parents[3] / "romanian_features.json"
)


def load_inventory() -> FeatureInventory:
    """Load the base JSON inventory with this analysis's patches applied."""
    return FeatureInventory.load(
        json_path=_INVENTORY_JSON,
        patches=PATCHES,
        underspec=UNDERSPEC,
    )


# ---------------------------------------------------------------------------
# Natural classes (LP [ ]) used by multiple rules
# ---------------------------------------------------------------------------
#
# Each of these is a natural class: a specification C, and its
# extension N(C) is every segment whose valued-feature set is a
# superset of C. Rules pattern-match on these; the ``.spec`` dict is
# passed to the rule constructors as the target / terminator /
# condition, all of which are natural classes in the paper's ``[ ]``
# notation.

# [-Syllabic, +Consonantal, -Sonorant, -Approximant, -Nasal] — the
# obstruent class targeted by dorsal palatalization and its default.
OBSTRUENT_CLASS: Final[NaturalClass] = natural_class(
    Syllabic="-",
    Consonantal="+",
    Sonorant="-",
    Approximant="-",
    Nasal="-",
)

# [+Coronal, -Anterior, +Distributed] — the postalveolar place class.
# /i/ and /j/ carry this place after patch; so does the derived [tʃ].
# S-palatalization's terminator (narrow scan) and assibilation's
# terminator both live here.
POSTALVEOLAR_PLACE: Final[NaturalClass] = natural_class(
    CORONAL="+",
    Anterior="-",
    Distributed="+",
)

# [+Syllabic, +Front] — front-vowel condition for dorsal palatalization's
# broad-terminator SEARCH.
FRONT_VOWEL: Final[NaturalClass] = natural_class(
    Syllabic="+",
    Front="+",
)

# The derived postalveolar consonant — the referent of the paper's
# circle-notation ``[Ⓢ]`` in the bleed rule. Its bundle collects
# everything /ʃ/-derived and /tʃ/-derived have in common, plus
# +Consonantal to exclude the vowel /i/ and glide /j/ (which are
# also postalveolar under our patch). NaturalClass.from_segment
# gives us [Ⓢ] directly from this bundle, matching the paper's
# ``[Ⓢ]`` notation.
_DERIVED_POSTALVEOLAR: Final[Segment] = Segment(
    label="Ⓢ",
    features={
        "CORONAL": "+",
        "Anterior": "-",
        "Distributed": "+",
        "Consonantal": "+",
    },
)
CIRCLED_S: Final[NaturalClass] = NaturalClass.from_segment(
    _DERIVED_POSTALVEOLAR,
)


# Backward-compat aliases (the rule constructors below use these as
# their .spec dicts; kept for internal use only).
_OBSTRUENT: Final[Mapping[str, str]] = OBSTRUENT_CLASS.spec
_POSTALVEOLAR_PLACE: Final[Mapping[str, str]] = POSTALVEOLAR_PLACE.spec
_FRONT_VOWEL: Final[Mapping[str, str]] = FRONT_VOWEL.spec
_POSTALVEOLAR_CONSONANT: Final[Mapping[str, str]] = CIRCLED_S.spec


# ---------------------------------------------------------------------------
# Rules (paper-revised set)
# ---------------------------------------------------------------------------

DORSAL_PAL: Final[UnificationRule] = UnificationRule(
    name="dorsal-pal",
    target=_OBSTRUENT,
    supply={
        "CORONAL": "+",
        "DORSAL": "-",
        "Anterior": "-",
        "Distributed": "+",
        "Strident": "+",
    },
    search=Search(
        direction=Direction.RIGHT,
        terminator={},  # broad next-segment terminator
        condition=_FRONT_VOWEL,
    ),
)
"""Dorsal palatalization (latex.tex:188-207 = e:velar-pal).

Target: any obstruent.
Supply: the full postalveolar bundle.
Trigger: broad next-segment terminator + [+Syll, +Front] condition.
The broad terminator is what blocks the rule in -ct- and -kt-: the
next segment after /K/ is the stop /T/, which fails the condition.
"""


DORSAL_DEFAULT: Final[UnificationRule] = UnificationRule(
    name="dorsal-default",
    target=_OBSTRUENT,
    supply={
        "CORONAL": "-",
        "DORSAL": "+",
        "Continuant": "-",
    },
    search=None,  # environment-free default fill
)
"""Dorsal default fill-in (latex.tex:211-225 = e:default).

Supplies the velar feature bundle to any obstruent still open for
those features — namely /K/, /G/ that escaped palatalization. /S/,
/T/, etc. are blocked because they have +Coronal which contradicts
-Coronal.
"""


S_PAL: Final[UnificationRule] = UnificationRule(
    name="s-pal-rev",
    target={
        "CORONAL": "+",
        "Strident": "+",
        "Continuant": "+",
    },
    supply={
        "Anterior": "-",
        "Distributed": "+",
    },
    search=Search(
        direction=Direction.RIGHT,
        terminator=_POSTALVEOLAR_PLACE,
        condition=_POSTALVEOLAR_PLACE,
    ),
)
"""S-palatalization, revised (latex.tex:508-527 = e:s-pal-rev).

Target: strident-coronal fricative class. Supply: postalveolar place
features. Trigger: narrow rightward search for a segment in the
postalveolar place class — /i/, /j/ (patched), or a derived [tʃ].
Non-members are transparent, so /S/ in 'prost' palatalizes across /T/.
"""


BLEED: Final[DeletionRule] = DeletionRule(
    name="bleed-rev",
    target=natural_class(CORONAL="+", Consonantal="+").spec,
    clear=frozenset({"Anterior", "Distributed", "Strident"}),
    # Paper (latex.tex:568-586): environment ``[Ⓢ] __`` — the circle
    # denotes the natural class of any segment subsuming the derived
    # postalveolar's bundle. We build it here via
    # ``NaturalClass.from_segment`` so the code reads the same as
    # the paper's ``[Ⓢ]``. Broad terminator = strict left-adjacency
    # (N(∅) as the terminator = adjacency; see LP notes §9).
    search=Search(
        direction=Direction.LEFT,
        terminator=NaturalClass.universal().spec,
        condition=CIRCLED_S.spec,
    ),
)
"""Bleed, revised (latex.tex:568-586 = e:bleed-rev).

The one feature-changing rule in the coronal set (LP ``\\`` operator).
Deletes place + Strident from a consonantal coronal immediately after
the derived postalveolar. Coronal default refills the unmarked plain-
stop values, giving /ʃt/ for both /st-i/ and /sk-e/.
"""


ASSIBILATION: Final[UnificationRule] = UnificationRule(
    name="assibilation",
    target={
        "CORONAL": "+",
        "Anterior": "+",
        "Sonorant": "-",
    },
    supply={"Strident": "+"},
    search=Search(
        direction=Direction.RIGHT,
        terminator=_POSTALVEOLAR_PLACE,
        condition=_POSTALVEOLAR_PLACE,
    ),
)
"""Assibilation (latex.tex:341-357 = e:assib).

Supplies +Strident to anterior coronal stops before a postalveolar
trigger. Aligned with S-pal on the trigger class for consistency
with the revised analysis (one place class covers both rules per
latex.tex:535).
"""


CORONAL_DEFAULT_STRIDENT: Final[UnificationRule] = UnificationRule(
    name="cor-default-strident",
    target={"CORONAL": "+"},
    supply={"Strident": "-"},
    search=None,
)
"""Coronal default: any coronal still open for Strident → -Strident.

First of three coronal-default sub-rules (latex.tex:370-396 =
e:cor-default). Resolves unmutated /T/ to plain [t], unmutated /D/
to plain [d].
"""


CORONAL_DEFAULT_PLACE: Final[UnificationRule] = UnificationRule(
    name="cor-default-place",
    target={"CORONAL": "+"},
    supply={"Anterior": "+", "Distributed": "-"},
    search=None,
)
"""Coronal default: any coronal still open for place → [+Ant, -Dist].

Second of three sub-rules. Resolves unmutated /S/ to plain [s],
bleached affricate to plain [t]/[d].
"""


CORONAL_DEFAULT_CONT: Final[UnificationRule] = UnificationRule(
    name="cor-default-continuant",
    target={"CORONAL": "+", "Strident": "+"},
    supply={"Continuant": "+"},
    search=None,
)
"""Coronal default: any +Strident coronal still open for Continuant.

Third of three sub-rules. Fires on assibilated /D/ (Continuant open),
giving [z]. Fails on assibilated /T/ (Continuant prespecified -), so
/T/ stays as the affricate [t͡s]. This is the source of the
asymmetric Romanian gap [t͡s] vs *[d͡z].
"""


GLIDE_FORMATION: Final[GlideFormationRule] = GlideFormationRule(
    name="glide-formation",
    requires_preceding={"Consonantal": "+"},
    glide_label="j",
    glide_supply={
        "Syllabic": "-",
        "CORONAL": "+",
        "Anterior": "-",
        "Distributed": "+",
    },
    clear=frozenset({"Syllabic"}),
    target_label="i",
)
"""Glide formation (latex.tex:404-421 = e:glide).

Word-final /i/ after a consonant desyllabifies to multiply-
articulated /j/. The supplied place features unify (no-op) with /i/'s
already-patched postalveolar place, and -Syllabic replaces +Syllabic.
"""


# ---------------------------------------------------------------------------
# Pipelines
# ---------------------------------------------------------------------------

RULES_FROM_PAPER: Final[tuple[Rule, ...]] = cast(
    "tuple[Rule, ...]",
    (
        DORSAL_PAL,
        DORSAL_DEFAULT,
        S_PAL,
        BLEED,
        ASSIBILATION,
        CORONAL_DEFAULT_STRIDENT,
        CORONAL_DEFAULT_PLACE,
        CORONAL_DEFAULT_CONT,
        GLIDE_FORMATION,
    ),
)
"""Paper's revised analysis: place-class S-pal and feature-deleting
bleed. Mirrors the derivation tables at latex.tex:472-492 and
latex.tex:584-596 column-for-column."""


# The four rules whose firing indicates the STEM ALTERNATED (as
# opposed to default fill-in, dorsal-default fill, or glide formation
# — all of which fire even in inalterable cases). Derived from the
# rule tuple itself so a rule rename can't leave stale copies in
# diagnostics/distance.py or the driver.
PALATALIZATION_RULES: Final[tuple[Rule, ...]] = cast(
    "tuple[Rule, ...]",
    (DORSAL_PAL, S_PAL, ASSIBILATION, BLEED),
)
PALATALIZATION_RULE_NAMES: Final[frozenset[str]] = frozenset(
    r.name for r in PALATALIZATION_RULES
)


# ---------------------------------------------------------------------------
# Stem-final classification (shared by driver + contingency table)
# ---------------------------------------------------------------------------

DORSAL_STEMS: Final[frozenset[str]] = frozenset({"c", "g"})
# /z/ excluded from CORONAL_STEMS: paper treats fully specified /z/
# as inalterable (latex.tex:170), so it never triggers palatalization.
# But /z/-final rows are still IN SCOPE for validation (they should
# predict False and match the data's False), so TARGET_STEMS keeps
# /z/ for the lexicon filter.
CORONAL_STEMS: Final[frozenset[str]] = frozenset({"t", "d", "s"})
TARGET_STEMS: Final[frozenset[str]] = (
    DORSAL_STEMS | CORONAL_STEMS | frozenset({"z"})
)


# ---------------------------------------------------------------------------
# Rule-firing oracle — one function, two views
# ---------------------------------------------------------------------------


def expected_firings(row: LexRow) -> frozenset[str]:
    """Which palatalization rules SHOULD fire for this row under the
    paper's productivity claims.

    Single source of truth for the rule-firing oracle. The boolean
    prediction is just ``bool(expected_firings(row))``, so we don't
    keep two shadow tables of the same stem × opportunity × cluster
    truth. Used by :func:`fallback_predict` and by the distance
    metric's rule-firing component.
    """
    if row.opportunity in {"uri", "none"}:
        return frozenset()
    if row.cluster is ClusterTag.CT:
        # ct+e: broad terminator halts on /T/, /K/ blocked.
        # ct+i: /T/ still assibilates.
        if row.opportunity == "i":
            return frozenset({ASSIBILATION.name})
        return frozenset()
    if row.cluster is ClusterTag.SC:
        # -sc(ă): /K/ palatalises, /S/ follows via derived postalv,
        # bleed clears the affricate to plain /t/.
        return frozenset({DORSAL_PAL.name, S_PAL.name, BLEED.name})
    if row.cluster is ClusterTag.SHCA:
        # -șcă already has surface /ʃ/; /K/ palatalises + bleed clears.
        return frozenset({DORSAL_PAL.name, BLEED.name})
    if row.cluster is ClusterTag.ST and row.opportunity == "i":
        return frozenset({S_PAL.name, BLEED.name})
    if row.stem_final in DORSAL_STEMS:
        if row.opportunity in {"i", "e"}:
            return frozenset({DORSAL_PAL.name})
    if row.stem_final in CORONAL_STEMS and row.opportunity == "i":
        if row.stem_final in {"t", "d"}:
            return frozenset({ASSIBILATION.name})
        return frozenset({S_PAL.name})
    return frozenset()


def fallback_predict(row: LexRow) -> bool:
    """Boolean prediction when we can't run the pipeline (missing IPA).

    Derived from :func:`expected_firings` so there's exactly one place
    to update if the rule set changes.
    """
    return bool(expected_firings(row))


# ---------------------------------------------------------------------------
# AnalysisProfile — the "one object to pass around" bundle
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class AnalysisProfile:
    """Everything an analysis run needs, bundled into one immutable
    object so downstream API functions take a single argument.

    Fields:

      - ``name`` — human-readable analysis name.
      - ``inventory`` — the constructed :class:`FeatureInventory`.
      - ``patches`` / ``underspec`` — the source-of-truth tuples used
        to build ``inventory`` (needed by the perturbation searcher
        to reload the inventory with tweaks).
      - ``rules`` — the ordered rule tuple.
      - ``inventory_json`` — path to the JSON that ``patches`` were
        applied on top of (used by perturb to reload).
      - ``target_stems`` — orthographic stem-finals this analysis
        wants to see; passed to :func:`iter_lexicon_rows`.
      - ``palatalization_rule_names`` — which rule names count as
        "mutation-producing"; passed to distance.score.
      - ``expected_firings`` — the rule-firing oracle (per-row).
      - ``fallback_predict`` — the boolean predictor for rows the
        pipeline can't process.
    """

    name: str
    inventory: FeatureInventory
    patches: tuple[FeaturePatch, ...]
    underspec: tuple[UnderspecifiedSegment, ...]
    rules: tuple[Rule, ...]
    inventory_json: Path
    target_stems: frozenset[str]
    palatalization_rule_names: frozenset[str]
    expected_firings: Callable[[LexRow], frozenset[str]]
    fallback_predict: Callable[[LexRow], bool]


def build_profile() -> AnalysisProfile:
    """Construct the default :class:`AnalysisProfile` for Romanian.

    Preferred entry point over :func:`load_inventory` alone — callers
    of the API pass the returned profile to
    :func:`phonology.validation.run_validation` etc.
    """
    inv = load_inventory()
    return AnalysisProfile(
        name="romanian-palatalization",
        inventory=inv,
        patches=PATCHES,
        underspec=UNDERSPEC,
        rules=RULES_FROM_PAPER,
        inventory_json=_INVENTORY_JSON,
        target_stems=TARGET_STEMS,
        palatalization_rule_names=PALATALIZATION_RULE_NAMES,
        expected_firings=expected_firings,
        fallback_predict=fallback_predict,
    )


__all__: Final[tuple[str, ...]] = (
    "ASSIBILATION",
    "AnalysisProfile",
    "BLEED",
    "CORONAL_DEFAULT_CONT",
    "CORONAL_DEFAULT_PLACE",
    "CORONAL_DEFAULT_STRIDENT",
    "CORONAL_STEMS",
    "DORSAL_DEFAULT",
    "DORSAL_PAL",
    "DORSAL_STEMS",
    "GLIDE_FORMATION",
    "PALATALIZATION_RULES",
    "PALATALIZATION_RULE_NAMES",
    "PATCHES",
    "RULES_FROM_PAPER",
    "S_PAL",
    "TARGET_STEMS",
    "UNDERSPEC",
    "build_profile",
    "expected_firings",
    "fallback_predict",
    "load_inventory",
)
