"""Italian velar palatalization: a provisional Logical Phonology analysis.

This module extends the analysis of ``latex.tex`` (Gorman & Winstead's
Romanian palatalization, :mod:`phonology.analyses.romanian_palatalization`)
to the Italian noun and adjective plurals counted by
:mod:`lexpipe.langs.italian` (amico/amici, medico/medici vs fuoco/fuochi,
lago/laghi; amica/amiche). It is the analysis as data: the inventory
patches, the underspecified dorsals, the ordered rules and the
specification policies of three variants, plus the hooks that feed the
generated Italian cell data to the shared validator
(:mod:`phonology.validation`, driven by ``scripts/validate_lang_rules.py
--lang it``). The justifications, the derivations and the results are in
``docs/italian/phonology.md``.

The analysis
------------
**Inventory.** ``phonology/inventories/italian_features.json``: the 27
features of the Romanian inventory, the segments Italian shares with
Romanian copied exactly, close-mid e o (Tense +) beside open-mid ɛ ɔ
(Tense -), d͡z, ɲ, ʎ, and geminates as single segments with Long +
(kː, t͡ʃː, ʃː, …). Wiktionary writes a geminate across the syllable dot
(ˈsak.ko, ˈfat.t͡ʃa); :func:`merge_geminates` rewrites it as the long
segment before :data:`ITALIAN_TOKENIZER` (which attaches ː) splits it.

**Patches** (:data:`PATCHES`). As in Romanian, /i/ and /j/ are
postalveolar-coronal vocoids, and the palatal affricates are dorso-coronal
with DelRel free, so the affricate derived from a velar projects onto
t͡ʃ d͡ʒ (and their geminates). ʃː is dorso-coronal too: it is what the
velar of an sc cluster becomes.

**Underspecification** (:data:`UNDERSPEC`). Only dorsals: /K G/ are k ɡ
minus Coronal (as in Romanian), /Kː Gː/ their geminates, and /Kˢ/ the
velar of an alternating sc unit (osco/osci [ʃː]), open for Coronal,
Continuant and Long. Italian has **no underspecified coronal**: no rule
writes coronal features into t d s z, so the absence of a coronal rule
(gatto/gatti, naso/nasi) is an inventory fact, and the coronal
disjuncts of G&W's spell-out condition P are vacuous.

**Rules** (:data:`RULES`, in order):

  1. ``dorsal-pal`` [obstruent] ⊔ {+Cor, -Ant, +Dist, +Str} / _ [-Cons,
     +Cor, -Ant, +Dist], with a broad (next-segment) terminator. The
     trigger is the coronal vocoid class, i.e. /i j/ under the patches:
     the Height effect (masculine -i triggers, feminine -e never does)
     is the Romanian CORONAL+ patch on /i j/, not a stipulation about
     cells. amica + e keeps its velar although its root is /amiK/.
  2. ``sc-coalescence`` [+Cor, -Ant, +Dist, +Str] ⊔ {+Cont, +Long} /
     [s] _: the palatalized /Kˢ/ of an sc unit becomes [ʃː].
  3. ``s-deletion`` [s] → ∅ / _ [ʃː]: the cluster's s gives way (no
     [sʃ] occurs in Italian), so /ɔsKˢ+i/ → [ɔʃːi], not Romanian [ʃt].
  4. ``dorsal-default`` [obstruent] ⊔ {-Cor, +Dor, -Cont} (latex.tex
     e:default): unpalatalized /K G Kː Gː Kˢ/ surface as velars.
  5. ``length-default`` [+Cons] ⊔ {-Long}: the unpalatalized /Kˢ/.
  6. ``glide-absorption`` /j/ → ∅ / _ /i/: occhio /ɔkːj+i/ → [ɔkːi].

Prespecified segments are inalterable by unification failure, as in
Romanian, so nothing here touches root-internal velars (chilo, ghiro:
NDEB) or the lexical palatals of luce, valigia.

Specification policies (the variants)
-------------------------------------
Spell-out is not modelled: the UR is the lemma's transcription minus
its final vowel plus the cell's ending (-i for the i cell, -e for the e
cell), as in the Romanian validator. What differs between the variants
is only which stem-final dorsals are underspecified:

  - ``lexical`` (:func:`build_profile`, Gorman & Winstead): the target
    is underspecified iff its **root** shows alternation evidence in E0
    (nominal plural cells) at τ = dictionary, computed with the TP
    stage's own definition (:mod:`lexpipe.langs.italian.roots`; R4:
    learner evidence pools all forms of a root, so amica shares amico's
    /K/). This is circular for the i cell: a root is underspecified
    because its -i plural palatalizes, so the boolean match there tests
    the rule mechanics (trigger, locality, projection, the cluster and
    geminate outputs), not a lexical prediction. What it can fail on is
    the e cell (an underspecified root before -e), roots with both
    outcomes (doublets: one UR, one output) and the surface IPA.
  - ``all_underspecified``: every dorsal target, NDE stems included, is
    underspecified: how far a single rule /K/ → [t͡ʃ] / _ i gets (the
    R17b palatalize row as a derivation).
  - ``stress``: underspecified iff the stem-final velar is **not**
    immediately post-tonic (velar_post_tonic of the cell data; unknown
    stress counts as not post-tonic): a learner-style grammar, the
    stress-nested TP rows (retain post-tonic, palatalize elsewhere) as a
    derivation.

The rule-firing oracle (:func:`expected_firings`) is the attested
outcome: Italian palatalization is lexical, so no productivity claim
can stand in for it; the distance metric's rule component then measures
disagreement with the data. ``fallback_predict`` (rows whose UR cannot
be built) applies the variant's own policy, never the data.

Rows and the attested plural
----------------------------
:func:`iter_rows` reads ``data/italian/italian_cells.csv`` (one row per
lexeme × plural cell × form), keeps the in-domain, aligned, unfiltered
dorsal-target rows before -i or -e whose lemma has a transcription, and
joins the lemma IPA (``italian_lexicon_complete.csv``) and the root
evidence. The data carry no plural transcriptions (design §3: plural IPA
is derived from the lemma), so the attested plural (``ipa_pl``) is
reconstructed by :func:`plural_ipa`: the lemma IPA minus its final vowel,
the stem-final velar replaced by the palatal the spelling shows (ci gi
→ t͡ʃ d͡ʒ, cci → t͡ʃː, sci → ʃː) or kept (chi ghi), a root glide absorbed
before -i, plus the ending. The IPA match therefore tests the derived
segments at the stem edge against the orthographic outcome; the rest of
the word is the lemma's own transcription.

Hooks implemented here because the shared ones cannot carry Italian
(reported in docs/italian/phonology.md): :class:`ItalianRow` subclasses
:class:`~phonology.lexicon.LexRow` to carry root evidence, stress and
class; the categoriser maps onto the shared
:class:`~phonology.diagnostics.categorise.RefinedCategory` (which has no
Italian labels) and :func:`mismatch_class` names the Italian mismatch
types for :func:`render_summary`.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from collections.abc import Callable, Collection, Iterator, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Final, cast

from lexpipe import PROJECT_ROOT
from lexpipe.spec import LANGUAGES, LanguageInfo

from ..diagnostics.categorise import Categorisation, RefinedCategory
from ..diagnostics.compare import NormaliserCascade, compare_ipa
from ..inventory import (
    FeatureInventory,
    FeaturePatch,
    UnderspecifiedSegment,
    inventory_path,
    load_inventory_file,
)
from ..lexicon import ClusterTag, LexRow
from ..lp import NaturalClass, natural_class
from ..profile import AnalysisProfile
from ..rules import Rule, SegmentDeletionRule, UnificationRule
from ..search import Direction, Search
from ..segments import (
    NON_SEGMENTAL_MARKS,
    LengthHandling,
    Segment,
    TokenizerConfig,
    Word,
    tokenize_ipa,
)
from ..validation import (
    RowResult,
    TraceResult,
    ValidationReport,
    evaluate_row,
    render_trace,
    run_validation,
)

# ---------------------------------------------------------------------------
# Language, target classes, specification policies
# ---------------------------------------------------------------------------

_INVENTORY_NAME: Final[str] = "italian"
LANGUAGE: Final[LanguageInfo] = LANGUAGES["it"]

DORSAL: Final[str] = "dorsal"
NDE_OCCHIO: Final[str] = "nde:occhio"
NDE_PSICHE: Final[str] = "nde:psiche"
CORONAL: Final[str] = "coronal"
FRONTSTEM: Final[str] = "frontstem"

DORSAL_CLASSES: Final[frozenset[str]] = frozenset(
    {DORSAL, NDE_OCCHIO, NDE_PSICHE}
)
"""Target classes whose stem-final dorsal a policy may underspecify."""

DEFAULT_CLASSES: Final[tuple[str, ...]] = (DORSAL, NDE_OCCHIO, NDE_PSICHE)
"""The rows the validator runs on: dorsal targets and the NDE classes."""

ALL_CLASSES: Final[tuple[str, ...]] = DEFAULT_CLASSES + (CORONAL, FRONTSTEM)
"""Adds the sanity rows: coronals (no coronal rule) and lexical palatals."""

LEXICAL: Final[str] = "lexical"
ALL_UNDERSPECIFIED: Final[str] = "all_underspecified"
STRESS: Final[str] = "stress"
POLICIES: Final[tuple[str, ...]] = (LEXICAL, ALL_UNDERSPECIFIED, STRESS)

DORSAL_STEMS: Final[frozenset[str]] = frozenset({"c", "g"})
CORONAL_STEMS: Final[frozenset[str]] = frozenset({"t", "d", "s", "z"})
TARGET_STEMS: Final[frozenset[str]] = DORSAL_STEMS | CORONAL_STEMS

# ---------------------------------------------------------------------------
# Geminates, tokenizer, comparison convention
# ---------------------------------------------------------------------------

TIE: Final[str] = "\N{COMBINING DOUBLE INVERTED BREVE}"
LENGTH: Final[str] = "ː"
_STRESS_MARKS: Final[frozenset[str]] = frozenset({"ˈ", "ˌ"})
_GEMINATE_SEPARATORS: Final[frozenset[str]] = _STRESS_MARKS | {"."}
_GEMINABLE: Final[frozenset[str]] = frozenset(
    {
        "p",
        "b",
        "t",
        "d",
        "k",
        "ɡ",
        "g",
        "f",
        "v",
        "s",
        "ʃ",
        "m",
        "n",
        "ɲ",
        "l",
        "ʎ",
        "r",
        f"t{TIE}s",
        f"d{TIE}z",
        f"t{TIE}ʃ",
        f"d{TIE}ʒ",
    }
)


def _units(ipa: str) -> list[str]:
    """Split IPA into code points, keeping a tied affricate together."""
    units: list[str] = []
    i = 0
    while i < len(ipa):
        if i + 2 < len(ipa) and ipa[i + 1] == TIE:
            units.append(ipa[i : i + 3])
            i += 3
        else:
            units.append(ipa[i])
            i += 1
    return units


def _second_half(first: str, second: str) -> bool:
    """Is *second* the second half of a geminate starting with *first*?

    The same segment (k.k, ʎ.ʎ), or an affricate whose stop is *first*
    (t.t͡ʃ, d.d͡z: Wiktionary writes the first half as a plain stop).
    """
    if second == first:
        return True
    return len(second) == 3 and second[1] == TIE and second[0] == first


def merge_geminates(ipa: str) -> str:
    """Rewrite each geminate as one long segment: ˈsak.ko → ˈsakːo.

    Two identical consonants (or a stop and the affricate it begins),
    adjacent or separated by syllable dots and stress marks, become the
    second one plus ː: ˈfat.t͡ʃa → ˈfat͡ʃːa, ˈpit.t͡sa → ˈpit͡sːa,
    viʎˈʎak.ko → viˈʎːakːo. A stress mark between the halves is kept
    in front of the long segment; dots between them are dropped.
    Vowels, glides and segments already marked long are left alone.
    Idempotent.
    """
    units = _units(ipa)
    out: list[str] = []
    i = 0
    while i < len(units):
        unit = units[i]
        if unit in _GEMINABLE:
            j = i + 1
            while j < len(units) and units[j] in _GEMINATE_SEPARATORS:
                j += 1
            already_long = j + 1 < len(units) and units[j + 1] == LENGTH
            if (
                j < len(units)
                and _second_half(unit, units[j])
                and not already_long
            ):
                out.extend(m for m in units[i + 1 : j] if m in _STRESS_MARKS)
                out.append(units[j] + LENGTH)
                i = j + 1
                continue
        out.append(unit)
        i += 1
    return "".join(out)


ITALIAN_TOKENIZER: Final[TokenizerConfig] = TokenizerConfig(
    skip=NON_SEGMENTAL_MARKS,
    # Wiktionary IPA uses the IPA letter ɡ (U+0261); an ASCII g in a
    # transcription means the same segment.
    fold={"g": "ɡ"},
    keep_combining=True,
    # Length is contrastive (geminates): ː ends the preceding token.
    length=LengthHandling.ATTACH,
)
"""Italian tokenizer convention (R8): stress marks and syllable dots are
skipped, ː attaches to its segment (kː, t͡ʃː), ɛ and ɔ stay distinct
from e and o (Romanian folds them), and a combining mark stays on its
base. Geminates written across a syllable dot are merged first
(:func:`merge_geminates`); :func:`tokens` does both."""


def tokens(ipa: str) -> tuple[str, ...]:
    """The inventory labels of an Italian transcription."""
    return tokenize_ipa(merge_geminates(ipa), tokenizer=ITALIAN_TOKENIZER)


_COMPARE_STRIP: Final[frozenset[str]] = frozenset(
    {"ˈ", "ˌ", ".", "'", "`", " "}
)


def _strip_italian(ipa: str) -> str:
    """Merge geminates, then drop stress marks, syllable dots and spaces.

    The length mark ː is kept: Italian length is contrastive.
    """
    return "".join(
        ch for ch in merge_geminates(ipa) if ch not in _COMPARE_STRIP
    )


ITALIAN_NORMALISERS: Final[NormaliserCascade] = NormaliserCascade(
    strip=_strip_italian,
    steps=(),
)
"""Italian comparison convention: strip-only. The attested plural is
reconstructed from the same lemma transcription as the UR (module
docstring), so no transcription convention needs smoothing; geminates are
merged on both sides and ː is kept."""


# ---------------------------------------------------------------------------
# Patches to the base inventory
# ---------------------------------------------------------------------------


def _dorso_coronal() -> dict[str, str]:
    """The affricate patch: DelRel free, +Dorsal (a fresh dict each)."""
    return {"DelRel": "0", "DORSAL": "+"}


PATCHES: Final[tuple[FeaturePatch, ...]] = (
    FeaturePatch(
        segment="i",
        updates={"CORONAL": "+", "Anterior": "-", "Distributed": "+"},
        reason=(
            "As in Romanian (latex.tex: a high front vocoid is coronal): "
            "/i/ is a postalveolar-place vocoid. This is the Italian "
            "Height effect: dorsal-pal's trigger is the coronal vocoid "
            "class, which holds /i j/ but not /e/ (CORONAL -), so the "
            "masculine -i triggers and the feminine -e never does."
        ),
    ),
    FeaturePatch(
        segment="j",
        updates={"CORONAL": "+", "Anterior": "-", "Distributed": "+"},
        reason=(
            "/j/ shares the place of /i/ (as in Romanian); a /j/-initial "
            "ending is a trigger too (2nd-conjugation -iamo)."
        ),
    ),
    FeaturePatch(
        segment="t͡ʃ",
        updates=_dorso_coronal(),
        reason=(
            "Romanian patch: palatalization leaves Dorsal untouched, so "
            "the affricate derived from /K/ is +Coronal and +Dorsal; "
            "DelRel is left free so the derived segment (DelRel - from "
            "/k/) projects onto t͡ʃ."
        ),
    ),
    FeaturePatch(
        segment="d͡ʒ",
        updates=_dorso_coronal(),
        reason="See t͡ʃ: the affricate derived from /G/.",
    ),
    FeaturePatch(
        segment="t͡ʃː",
        updates=_dorso_coronal(),
        reason="See t͡ʃ: the geminate derived from /Kː/ (tombacco/tombacci).",
    ),
    FeaturePatch(
        segment="d͡ʒː",
        updates=_dorso_coronal(),
        reason="See t͡ʃ: the geminate derived from /Gː/.",
    ),
    FeaturePatch(
        segment="ʃː",
        updates={"DORSAL": "+"},
        reason=(
            "The velar of an sc unit coalesces into [ʃː] (osco/osci); it "
            "keeps its +Dorsal like the derived affricate, so ʃː is "
            "dorso-coronal. Fricatives have DelRel 0 already. Singleton "
            "ʃ is never derived and keeps its base values."
        ),
    ),
    # No t͡s/d͡z patch: Romanian frees t͡s's DelRel because assibilation
    # derives it; Italian has no assibilation, so no coronal affricate
    # is ever derived. No t/d Strident patch either (Romanian's protects
    # them from assibilation).
)


# ---------------------------------------------------------------------------
# Underspecified segments (dorsals only)
# ---------------------------------------------------------------------------

UNDERSPEC: Final[tuple[UnderspecifiedSegment, ...]] = (
    UnderspecifiedSegment(
        label="K",
        base="k",
        clear=frozenset({"CORONAL"}),
        reason=(
            "As in Romanian (latex.tex): /K/ = /k/ minus Coronal; +Dorsal "
            "and -Continuant survive, so palatalization gives the "
            "affricate t͡ʃ."
        ),
    ),
    UnderspecifiedSegment(
        label="G",
        base="ɡ",
        clear=frozenset({"CORONAL"}),
        reason="Voiced counterpart of /K/ (asparago/asparagi).",
    ),
    UnderspecifiedSegment(
        label="Kː",
        base="kː",
        clear=frozenset({"CORONAL"}),
        reason="Geminate /K/ (tombacco/tombacci [t͡ʃː]); Long stays +.",
    ),
    UnderspecifiedSegment(
        label="Gː",
        base="ɡː",
        clear=frozenset({"CORONAL"}),
        reason="Geminate /G/, the gg counterpart of /Kː/.",
    ),
    UnderspecifiedSegment(
        label="Kˢ",
        base="k",
        clear=frozenset({"CORONAL", "Continuant", "Long"}),
        reason=(
            "The velar of an alternating sc unit (osco/osci [ʃː], "
            "falisco/falisci): open for Continuant and Long as well, "
            "which sc-coalescence fills (+, +) after s and the defaults "
            "fill (-, -) otherwise. The s stays prespecified: Italian has "
            "no underspecified coronal."
        ),
    ),
)

UNDERSPEC_FOR: Final[Mapping[str, str]] = {
    "k": "K",
    "ɡ": "G",
    "kː": "Kː",
    "ɡː": "Gː",
}
"""Velar token → its underspecified twin (sc units take /Kˢ/)."""

SC_VELAR: Final[str] = "Kˢ"


def load_inventory() -> FeatureInventory:
    """The Italian inventory with this analysis's patches applied."""
    return load_inventory_file(
        _INVENTORY_NAME,
        patches=PATCHES,
        underspec=UNDERSPEC,
    )


# ---------------------------------------------------------------------------
# Natural classes (LP [ ])
# ---------------------------------------------------------------------------

OBSTRUENT_CLASS: Final[NaturalClass] = natural_class(
    Syllabic="-",
    Consonantal="+",
    Sonorant="-",
    Approximant="-",
    Nasal="-",
)
"""The target class of dorsal-pal and of the dorsal default."""

POSTALVEOLAR_PLACE: Final[NaturalClass] = natural_class(
    CORONAL="+",
    Anterior="-",
    Distributed="+",
)
"""[+Cor, -Ant, +Dist]: /i j/ under the patches, the palatals."""

TRIGGER_CLASS: Final[NaturalClass] = natural_class(
    **POSTALVEOLAR_PLACE.spec,
    Consonantal="-",
)
"""The coronal vocoids: /i j/ under the patches (never /e/). Over this
inventory [+High, +Front] picks the same two segments; the place
statement is G&W's revised trigger ('/i/ was arguably a place trigger
from the start')."""

POSTALVEOLAR_STRIDENT: Final[NaturalClass] = natural_class(
    **POSTALVEOLAR_PLACE.spec,
    Strident="+",
)
"""Target of sc-coalescence: only a palatalized /Kˢ/ is open for the
Continuant and Long it supplies; lexical t͡ʃ ʃ resist (prespecified)."""

ANTERIOR_SIBILANT: Final[NaturalClass] = natural_class(
    CORONAL="+",
    Anterior="+",
    Strident="+",
    Continuant="+",
    Voice="-",
)
"""[s]: the left context of sc-coalescence and the target of s-deletion."""

LONG_POSTALVEOLAR_FRICATIVE: Final[NaturalClass] = natural_class(
    CORONAL="+",
    Anterior="-",
    Strident="+",
    Continuant="+",
    Long="+",
)
"""[ʃː]: the right context of s-deletion (no [sʃ] occurs in Italian)."""

GLIDE_J: Final[NaturalClass] = natural_class(
    Syllabic="-",
    Consonantal="-",
    High="+",
    Front="+",
)

HIGH_FRONT_VOWEL: Final[NaturalClass] = natural_class(
    Syllabic="+",
    High="+",
    Front="+",
)

_NEXT_SEGMENT: Final[Mapping[str, str]] = NaturalClass.universal().spec


# ---------------------------------------------------------------------------
# Rules
# ---------------------------------------------------------------------------

DORSAL_PAL: Final[UnificationRule] = UnificationRule(
    name="dorsal-pal",
    target=OBSTRUENT_CLASS.spec,
    supply={
        "CORONAL": "+",
        "Anterior": "-",
        "Distributed": "+",
        "Strident": "+",
    },
    search=Search(
        direction=Direction.RIGHT,
        terminator=_NEXT_SEGMENT,
        condition=TRIGGER_CLASS.spec,
    ),
)
"""Dorsal palatalization (latex.tex e:velar-pal) with the Italian trigger.

Supply as in Romanian: +Coronal, -Anterior, +Distributed, +Strident;
Dorsal is untouched, so the output is dorso-coronal. Broad terminator:
only the next segment can trigger, and it must be a coronal vocoid
(/i j/). Prespecified obstruents resist by unification failure (k ɡ are
-Coronal, t d s z are +Anterior)."""

SC_COALESCENCE: Final[UnificationRule] = UnificationRule(
    name="sc-coalescence",
    target=POSTALVEOLAR_STRIDENT.spec,
    supply={"Continuant": "+", "Long": "+"},
    search=Search(
        direction=Direction.LEFT,
        terminator=_NEXT_SEGMENT,
        condition=ANTERIOR_SIBILANT.spec,
    ),
)
"""sc + i → [ʃː], first step: the palatalized /Kˢ/ right after [s] takes
+Continuant and +Long (the two consonants of the unit fuse into one long
fricative). Vacuous or blocked on every prespecified segment."""

S_DELETION: Final[SegmentDeletionRule] = SegmentDeletionRule(
    name="s-deletion",
    target=ANTERIOR_SIBILANT.spec,
    search=Search(
        direction=Direction.RIGHT,
        terminator=_NEXT_SEGMENT,
        condition=LONG_POSTALVEOLAR_FRICATIVE.spec,
    ),
)
"""sc + i → [ʃː], second step: [s] is deleted before [ʃː]. Italian has no
[sʃ] sequence, so the rule never reaches a lexical one; Romanian instead
keeps both ([muʃte]) and bleeds the affricate (latex.tex e:bleed-rev)."""

DORSAL_DEFAULT: Final[UnificationRule] = UnificationRule(
    name="dorsal-default",
    target=OBSTRUENT_CLASS.spec,
    supply={"CORONAL": "-", "DORSAL": "+", "Continuant": "-"},
    search=None,
)
"""Dorsal default fill-in, latex.tex e:default literally ({-Coronal,
+Dorsal, -Continuant}): unpalatalized /K G Kː Gː Kˢ/ become velar stops.
Vacuous on k ɡ and blocked on everything else."""

LENGTH_DEFAULT: Final[UnificationRule] = UnificationRule(
    name="length-default",
    target={"Consonantal": "+"},
    supply={"Long": "-"},
    search=None,
)
"""-Long into a consonant still open for length: only an unpalatalized
/Kˢ/ (bosco-type stems in the variants that underspecify them)."""

GLIDE_ABSORPTION: Final[SegmentDeletionRule] = SegmentDeletionRule(
    name="glide-absorption",
    target=GLIDE_J.spec,
    search=Search(
        direction=Direction.RIGHT,
        terminator=_NEXT_SEGMENT,
        condition=HIGH_FRONT_VOWEL.spec,
    ),
)
"""/j/ → ∅ / _ /i/: the -i ending absorbs a root glide (occhio /ɔkːj+i/
→ [ɔkːi], studio → studi). Italian has no [ji] sequence. Not a
palatalization rule: the occhio class keeps its velar."""

RULES: Final[tuple[Rule, ...]] = cast(
    "tuple[Rule, ...]",
    (
        DORSAL_PAL,
        SC_COALESCENCE,
        S_DELETION,
        DORSAL_DEFAULT,
        LENGTH_DEFAULT,
        GLIDE_ABSORPTION,
    ),
)
"""The ordered rule set; the same for every variant."""

PALATALIZATION_RULES: Final[tuple[Rule, ...]] = cast(
    "tuple[Rule, ...]",
    (DORSAL_PAL, SC_COALESCENCE, S_DELETION),
)
"""The rules whose firing means the stem alternated (the defaults and
glide absorption fire without any alternation)."""

PALATALIZATION_RULE_NAMES: Final[frozenset[str]] = frozenset(
    r.name for r in PALATALIZATION_RULES
)


# ---------------------------------------------------------------------------
# Rows
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class ItalianRow(LexRow):
    """A :class:`LexRow` with the Italian evidence the hooks need.

    The shared row has no place for a root, its alternation evidence or
    the stress of its velar; the hooks of this analysis read these
    fields (and treat a plain :class:`LexRow` as having none).

      - ``target_class`` — dorsal, nde:occhio, nde:psiche, coronal or
        frontstem;
      - ``target_letters`` — the alternating letters: c g sc cc gg for
        dorsals (the data layer's reading), else the stem-final letter;
      - ``velar_post_tonic`` — the stem-final velar is in the syllable
        right after the stress (None when unknown or not a dorsal);
      - ``root_alternates`` — the root has a palatalized E0 form at
        τ = dictionary (the ``lexical`` policy);
      - ``root_variable`` — the root's -i forms include both a
        palatalized and a velar one (a doublet);
      - ``suffix_class`` — -ico or -iaco when the lemma is an R26 suffix
        target; ``suffixoid``/``semantic_class`` — the root's R27
        -logo/-fago/-urgo class; ``deverbal_from_are`` — the carico
        type.
    """

    lexeme_id: str = ""
    root_id: str = ""
    cell: str = ""
    target_class: str = ""
    target_letters: str = ""
    cluster_letters: str = ""
    velar_post_tonic: bool | None = None
    stress_class: str = ""
    root_alternates: bool = False
    root_variable: bool = False
    suffix_class: str = ""
    suffixoid: str = ""
    semantic_class: str = ""
    deverbal_from_are: bool = False


_FINAL_VOWEL_CHARS: Final[frozenset[str]] = frozenset({"a", "e", "o"})
_SUFFIXES: Final[frozenset[str]] = frozenset({"i", "e"})
_VELAR_LETTER: Final[Mapping[str, str]] = {
    "k": "c",
    "kː": "c",
    "ɡ": "g",
    "ɡː": "g",
}
_PALATAL_FOR: Final[Mapping[str, str]] = {
    "k": f"t{TIE}ʃ",
    "kː": f"t{TIE}ʃ{LENGTH}",
    "ɡ": f"d{TIE}ʒ",
    "ɡː": f"d{TIE}ʒ{LENGTH}",
}
_STEM_EDGE_RE: Final[re.Pattern[str]] = re.compile(
    r"(?P<velar>kː|ɡː|k|ɡ)(?P<root>[ji]?)(?P<marks>[.ˈˌ]*)\Z"
)
_S_BEFORE_RE: Final[re.Pattern[str]] = re.compile(r"s[.ˈˌ]*\Z")


def plural_ipa(
    lemma_ipa: str,
    opportunity: str,
    *,
    stem_final: str,
    target_class: str,
    target_letters: str,
    mutation: bool,
) -> str:
    """Reconstruct the attested plural's IPA from the lemma's.

    The data have no plural transcriptions, so the plural is the lemma
    transcription with its final vowel (a e o) replaced by the ending,
    the stem-final velar of a dorsal target replaced by the palatal the
    spelling shows when *mutation* is True (k ɡ kː ɡː → t͡ʃ d͡ʒ t͡ʃː
    d͡ʒː; s + k of an sc unit → ʃː) and kept otherwise, and the root
    glide of an occhio stem absorbed before -i. Stress marks and
    syllable dots are kept; geminates are merged (:func:`merge_geminates`).

    Returns "" when the transcription does not fit the spelling (no
    final vowel, no velar where the stem has c/g, no s before an sc
    velar, a root vocoid outside the occhio class or missing in it).
    """
    body = merge_geminates(lemma_ipa.strip())
    if (
        opportunity not in _SUFFIXES
        or len(body) < 2
        or body[-1] not in _FINAL_VOWEL_CHARS
    ):
        return ""
    stem = body[:-1]
    if target_class not in DORSAL_CLASSES:
        return stem + opportunity
    edge = _STEM_EDGE_RE.search(stem)
    if edge is None or _VELAR_LETTER[edge["velar"]] != stem_final:
        return ""
    velar, root, marks = edge["velar"], edge["root"], edge["marks"]
    if bool(root) != (target_class == NDE_OCCHIO):
        return ""
    head = stem[: edge.start()]
    s_before = _S_BEFORE_RE.search(head)
    if target_letters == "sc" and s_before is None:
        return ""
    if not mutation:
        out = velar
    elif s_before is not None and target_letters == "sc":
        head = head[: s_before.start()]  # s + k fuse into one long [ʃː]
        out = "ʃ" + LENGTH
    else:
        out = _PALATAL_FOR[velar]
    if root == "j" and opportunity == "i":
        root = ""  # occhio → occhi: the ending absorbs the root glide
    return head + out + root + marks + opportunity


def _bool(text: str | None) -> bool:
    return (text or "").strip().lower() == "true"


def _optional_bool(text: str | None) -> bool | None:
    value = (text or "").strip().lower()
    if value == "true":
        return True
    if value == "false":
        return False
    return None


def target_class_of(r: Mapping[str, str]) -> str:
    """The row's target class from the cell table's columns, or ""."""
    nde = (r.get("nde_class") or "").strip()
    if nde in ("occhio", "psiche"):
        return f"nde:{nde}"
    target = (r.get("target") or "").strip()
    if target in (DORSAL, CORONAL):
        return target
    if _bool(r.get("frontstem")):
        return FRONTSTEM
    return ""


def _target_letters(r: Mapping[str, str], target_class: str) -> str:
    stem_final = (r.get("stem_final") or "").strip()
    cluster = (r.get("cluster") or "").strip()
    if target_class == DORSAL and cluster in ("sc", "cc", "gg"):
        return cluster
    return stem_final


_POS: Final[Mapping[str, str]] = {"noun": "N", "adj": "ADJ"}


def row_factory(r: Mapping[str, str]) -> ItalianRow:
    """One cell-table row (plus joined columns) → :class:`ItalianRow`.

    *r* is a row of ``italian_cells.csv`` to which :func:`iter_rows`
    adds ``lemma_ipa`` (the lemma's normalized transcription),
    ``root_alternates``, ``root_variable``, ``suffixoid``,
    ``semantic_class`` and ``deverbal_from_are``; missing columns read
    as empty. pos becomes N/ADJ (the vocabulary of
    :func:`phonology.lexicon.iter_lexicon_rows`), ``ipa_pl`` is
    :func:`plural_ipa` of the lemma transcription.
    """
    target_class = target_class_of(r)
    stem_final = (r.get("stem_final") or "").strip()
    letters = _target_letters(r, target_class)
    opportunity = (r.get("opportunity") or "").strip()
    mutation = _bool(r.get("mutation"))
    lemma_ipa = (r.get("lemma_ipa") or "").strip()
    lemma_suffix = (r.get("lemma_suffix") or "").strip()
    suffix_class = (
        lemma_suffix
        if _bool(r.get("target_is_suffix"))
        and lemma_suffix in ("-ico", "-iaco")
        else ""
    )
    post_tonic = (
        _optional_bool(r.get("velar_post_tonic"))
        if target_class in DORSAL_CLASSES
        else None
    )
    pos = (r.get("pos") or "").strip()
    return ItalianRow(
        lemma=(r.get("word") or "").strip(),
        plural=(r.get("form") or "").strip(),
        pos=_POS.get(pos, pos),
        stem_final=stem_final,
        opportunity=opportunity,
        mutation=mutation,
        ipa_lemma=lemma_ipa,
        ipa_pl=(
            plural_ipa(
                lemma_ipa,
                opportunity,
                stem_final=stem_final,
                target_class=target_class,
                target_letters=letters,
                mutation=mutation,
            )
            if lemma_ipa
            else ""
        ),
        exception_reason=(r.get("exception_reason") or "").strip(),
        cluster=ClusterTag.NONE,
        lexeme_id=(r.get("lexeme_id") or "").strip(),
        root_id=(r.get("root_id") or "").strip(),
        cell=(r.get("cell") or "").strip(),
        target_class=target_class,
        target_letters=letters,
        cluster_letters=(r.get("cluster") or "").strip(),
        velar_post_tonic=post_tonic,
        stress_class=(r.get("stress_class") or "").strip(),
        root_alternates=_bool(r.get("root_alternates")),
        root_variable=_bool(r.get("root_variable")),
        suffix_class=suffix_class,
        suffixoid=(r.get("suffixoid") or "").strip(),
        semantic_class=(r.get("semantic_class") or "").strip(),
        deverbal_from_are=_bool(r.get("deverbal_from_are")),
    )


@dataclass(frozen=True, slots=True)
class RootEvidence:
    """What a root's E0 forms (τ = dictionary) tell a learner."""

    alternates: bool
    variable: bool
    suffixoid: str
    semantic_class: str
    deverbal_from_are: bool


@dataclass(frozen=True, slots=True)
class Coverage:
    """How many in-scope cell rows reach the validator.

    ``in_scope`` counts the in-domain, aligned, unfiltered rows of the
    requested classes before -i or -e; ``no_ipa`` those whose lemma has
    no transcription (not validated); ``by_class`` maps a class to
    (in scope, with a transcription).
    """

    in_scope: int
    no_ipa: int
    by_class: Mapping[str, tuple[int, int]]

    @property
    def validated(self) -> int:
        return self.in_scope - self.no_ipa


def data_paths(root: Path) -> tuple[Path, Path]:
    """``(italian_cells.csv, italian_lexicon_complete.csv)`` under *root*
    (R3: paths derive from the registry row)."""
    data = LANGUAGE.data_dir(root)
    return (
        data / LANGUAGE.file_name("cells"),
        data / LANGUAGE.file_name("lexicon_complete"),
    )


def root_evidence(
    lexicon_complete: Path, cells: Path
) -> dict[str, RootEvidence]:
    """E0 evidence of every in-domain root (``lexical`` policy, R4/R15).

    Computed with the TP stage's own definitions
    (:func:`lexpipe.langs.italian.roots.load_records`, ``add_cells`` and
    ``realizations_for``): a root alternates iff one of its E0
    realizations before a front vowel is palatal and attested at
    τ = dictionary (filtered variants, curator-rejected and unaligned
    forms are not realizations). It is variable iff its -i realizations
    include both a palatal and a velar one.
    """
    from lexpipe.langs.italian import roots

    records = roots.load_records(lexicon_complete)
    roots.add_cells(records, cells)
    evidence: dict[str, RootEvidence] = {}
    for root_id, record in records.items():
        items = [
            item
            for item in roots.realizations_for(record.realizations, "E0")
            if item.attested.get("dictionary", False)
        ]
        before_i = {item.palatal for item in items if item.vowel == "i"}
        evidence[root_id] = RootEvidence(
            alternates=any(item.palatal is True for item in items),
            variable=True in before_i and False in before_i,
            suffixoid=record.suffixoid,
            semantic_class=record.semantic_class,
            deverbal_from_are=record.deverbal_from_are,
        )
    return evidence


def lemma_transcriptions(lexicon_complete: Path) -> dict[str, str]:
    """lexeme_id → the normalized primary transcription of the lemma.

    ``ipa_raw`` normalized by the Italian
    :class:`lexpipe.normalize.Normalizer` (affricates tied, ':' → ː),
    exactly as the lexicon stage computes its ``ipa`` column.
    """
    from lexpipe.langs.italian import keys

    normalizer = keys.normalizer()
    out: dict[str, str] = {}
    with lexicon_complete.open(encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh):
            raw = (r.get("ipa_raw") or "").strip()
            if not raw or not _bool(r.get("in_domain")):
                continue
            ipa = normalizer.ipa(raw)
            if ipa is not None:
                # A few transcriptions hold variants inside one pair of
                # slashes (/podˈd͡ʒin.ko | poˈʒin.ko/): keep the first.
                out[r["lexeme_id"]] = ipa.text.split("|")[0].strip()
    return out


def _in_scope(r: Mapping[str, str], classes: Collection[str]) -> str:
    """The row's class if it is validated at all, else ""."""
    if (
        not _bool(r.get("in_domain"))
        or not _bool(r.get("aligned"))
        or _bool(r.get("filtered"))
        or _bool(r.get("adjudicated_out"))
        or (r.get("opportunity") or "") not in _SUFFIXES
    ):
        return ""
    target_class = target_class_of(r)
    return target_class if target_class in classes else ""


def load_rows(
    root: Path,
    *,
    classes: Collection[str] = DEFAULT_CLASSES,
) -> tuple[tuple[ItalianRow, ...], Coverage]:
    """The validated rows and their coverage (see :func:`iter_rows`)."""
    cells, complete = data_paths(root)
    for path in (cells, complete):
        if not path.is_file():
            raise FileNotFoundError(
                f"{path} is missing: generate it with "
                "`python scripts/lang_pipeline.py --lang it features`"
            )
    evidence = root_evidence(complete, cells)
    ipa = lemma_transcriptions(complete)
    rows: list[ItalianRow] = []
    in_scope: Counter[str] = Counter()
    with_ipa: Counter[str] = Counter()
    with cells.open(encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh):
            target_class = _in_scope(r, classes)
            if not target_class:
                continue
            in_scope[target_class] += 1
            lemma_ipa = ipa.get(r["lexeme_id"], "")
            if not lemma_ipa:
                continue
            with_ipa[target_class] += 1
            found = evidence.get(r["root_id"])
            joined = dict(r)
            joined["lemma_ipa"] = lemma_ipa
            if found is not None:
                joined["root_alternates"] = str(found.alternates)
                joined["root_variable"] = str(found.variable)
                joined["suffixoid"] = found.suffixoid
                joined["semantic_class"] = found.semantic_class
                joined["deverbal_from_are"] = str(found.deverbal_from_are)
            rows.append(row_factory(joined))
    coverage = Coverage(
        in_scope=sum(in_scope.values()),
        no_ipa=sum(in_scope.values()) - sum(with_ipa.values()),
        by_class={
            name: (in_scope[name], with_ipa[name]) for name in sorted(in_scope)
        },
    )
    return tuple(rows), coverage


def iter_rows(
    root: Path,
    *,
    classes: Collection[str] = DEFAULT_CLASSES,
) -> Iterator[LexRow]:
    """The rows to validate (driver contract of
    ``scripts/validate_lang_rules.py``).

    Reads ``root/data/italian/italian_cells.csv`` and keeps the
    in-domain, aligned, unfiltered, not curator-rejected rows before -i
    or -e whose target class is in *classes* (default: dorsal targets
    and the NDE classes; :data:`ALL_CLASSES` adds coronals and lexical
    palatals) and whose lemma has a transcription: a row without one
    cannot be derived, and scoring it with ``fallback_predict`` would
    only restate the policy. Rows come in the cell table's order
    (lexeme_id, cell, item); every in-scope realization is a row, so a
    doublet contributes one row per form and amico N/ADJ one each.
    """
    rows, _ = load_rows(root, classes=classes)
    yield from rows


# ---------------------------------------------------------------------------
# Specification policies and the UR builder
# ---------------------------------------------------------------------------


def is_underspecified(row: LexRow, policy: str) -> bool:
    """Does *policy* underspecify the row's stem-final dorsal?

    Only dorsal targets and the NDE classes have one; coronals and
    lexical palatals are prespecified under every policy.
    """
    if not isinstance(row, ItalianRow):
        return False
    if row.target_class not in DORSAL_CLASSES:
        return False
    if policy == LEXICAL:
        return row.root_alternates
    if policy == ALL_UNDERSPECIFIED:
        return True
    if policy == STRESS:
        return row.velar_post_tonic is not True
    raise ValueError(f"unknown policy {policy!r} (known: {POLICIES})")


def ur_labels(row: LexRow, policy: str) -> tuple[str, ...] | None:
    """The UR as inventory labels: lemma tokens minus the final vowel,
    the stem-final dorsal swapped for its underspecified twin when
    *policy* says so, plus the cell's ending.

    None when the row cannot be built: no transcription, no final
    a/e/o, an ending other than -i/-e, or a dorsal target whose
    transcription has no velar where the spelling has c/g (for sc, no
    s before it; for occhio, no root vocoid after it).
    """
    if not isinstance(row, ItalianRow) or row.opportunity not in _SUFFIXES:
        return None
    toks = list(tokens(row.ipa_lemma))
    if len(toks) < 2 or toks[-1] not in _FINAL_VOWEL_CHARS:
        return None
    stem = toks[:-1]
    if row.target_class in DORSAL_CLASSES:
        index = len(stem) - 1
        if row.target_class == NDE_OCCHIO:
            if stem[index] not in ("j", "i"):
                return None
            index -= 1
        if index < 0 or _VELAR_LETTER.get(stem[index]) != row.stem_final:
            return None
        sc = row.target_letters == "sc"
        if sc and (index == 0 or stem[index - 1] != "s"):
            return None
        if is_underspecified(row, policy):
            stem[index] = SC_VELAR if sc else UNDERSPEC_FOR[stem[index]]
    return tuple(stem) + (row.opportunity,)


def build_ur(
    row: LexRow, inventory: FeatureInventory, policy: str
) -> Word | None:
    """The UR of *row* under *policy* as segments (None: see
    :func:`ur_labels`, or a label the inventory lacks)."""
    labels = ur_labels(row, policy)
    if labels is None:
        return None
    segments: list[Segment] = []
    for label in labels:
        if label not in inventory.base and label not in inventory.underspec:
            return None
        segments.append(inventory.segment(label))
    return tuple(segments)


@dataclass(frozen=True, slots=True)
class UrBuilder:
    """The ``ur_builder`` hook of one policy (equal by policy name)."""

    policy: str

    def __call__(
        self,
        row: LexRow,
        inventory: FeatureInventory,
        /,
    ) -> Word | None:
        return build_ur(row, inventory, self.policy)


@dataclass(frozen=True, slots=True)
class FallbackPredictor:
    """``fallback_predict`` of one policy: the policy's own prediction
    (an underspecified dorsal before -i), never the attested outcome."""

    policy: str

    def __call__(self, row: LexRow) -> bool:
        return row.opportunity == "i" and is_underspecified(row, self.policy)


def expected_firings(row: LexRow) -> frozenset[str]:
    """The palatalization rules that must fire to derive the attested
    plural: dorsal-pal for a palatalized form, plus sc-coalescence and
    s-deletion for an sc unit. Italian palatalization is lexical, so
    the oracle is the data, not a productivity claim."""
    if not row.mutation:
        return frozenset()
    names = {DORSAL_PAL.name}
    if isinstance(row, ItalianRow) and row.target_letters == "sc":
        names |= {SC_COALESCENCE.name, S_DELETION.name}
    return frozenset(names)


# ---------------------------------------------------------------------------
# Categoriser, classes, mismatch types
# ---------------------------------------------------------------------------


def categorise(
    row: LexRow,
    predicted_sr: str,
    attested_field: str,
    boolean_predicted: bool,
    /,
    *,
    ur_built: bool,
    normalisers: NormaliserCascade,
) -> Categorisation:
    """The categoriser hook, on the shared vocabulary only.

    :class:`RefinedCategory` has no Italian labels (adding some would
    mean editing the shared enum), and the Italian cascade has no
    normalisation steps, so a mismatch is UR_BUILD_FAILED or
    RULE_FAILURE, with the edit distance to the attested plural.
    :func:`mismatch_class` names the Italian cause and
    :func:`render_summary` reports it.
    """
    if not ur_built:
        return Categorisation(RefinedCategory.UR_BUILD_FAILED, (), 999.0, "")
    cmp = compare_ipa(predicted_sr, attested_field, normalisers=normalisers)
    return Categorisation(
        RefinedCategory.RULE_FAILURE,
        cmp.normalisation_applied,
        cmp.edit_distance,
        cmp.attested_variant,
    )


def row_class(row: LexRow) -> str:
    """The class a row is reported under (root-level classes, as the TP
    tables use): the NDE classes, coronal, frontstem, the R26 -ico/-iaco
    suffix targets, the R27 -logo/-fago/-urgo roots by person/thing, and
    the other dorsal roots by the stress of the velar."""
    if not isinstance(row, ItalianRow):
        return "other"
    if row.target_class != DORSAL:
        return row.target_class or "other"
    if row.suffix_class:
        return f"suffix {row.suffix_class}"
    if row.suffixoid:
        return f"{row.suffixoid} {row.semantic_class or 'unknown'}"
    if row.velar_post_tonic is True:
        return "root, velar post-tonic"
    if row.velar_post_tonic is False:
        return "root, velar not post-tonic"
    return "root, stress unknown"


def cluster_key(row: LexRow) -> str | None:
    """The cluster-report bucket: ``+<ending> <class>``."""
    if not isinstance(row, ItalianRow):
        return None
    return f"+{row.opportunity} {row_class(row)}"


def mismatch_class(result: RowResult) -> str | None:
    """The Italian name of a row's mismatch (None for a full match).

    Named from the data side (what the row is), so the same row gets
    the same name under every policy; which policies produce it is what
    :func:`render_summary` compares.
    """
    if result.boolean_match and result.ipa_match:
        return None
    row = result.row
    if result.ur is None:
        return "no UR (unknown segment, or IPA does not fit the spelling)"
    if not isinstance(row, ItalianRow):
        return "other"
    if result.boolean_match:
        return "surface IPA only"
    if row.target_class != DORSAL:
        return f"{row.target_class} stem palatalized"
    if row.opportunity == "e":
        if row.mutation:
            return "-e palatalized in the data (QC suspect e_palatal)"
        return "-e palatalized by the rules"
    if row.root_variable:
        return "doublet root (palatal and velar -i plurals)"
    if result.boolean_predicted:
        if row.suffixoid and row.semantic_class != "person":
            return f"{row.suffixoid} {row.semantic_class or 'unknown'} kept"
        if row.velar_post_tonic is True:
            return "post-tonic velar kept (fuoco type)"
        if row.deverbal_from_are:
            return "deverbal from -are kept (carico type)"
        if row.velar_post_tonic is None:
            return "velar kept, stress unknown"
        if row.stem_final == "g":
            return "g-stem not post-tonic kept (profugo type)"
        return "velar not post-tonic kept (other)"
    if row.velar_post_tonic is True:
        return "post-tonic palatalizer (amico type)"
    return "palatalizer, velar not post-tonic"


# ---------------------------------------------------------------------------
# Profiles
# ---------------------------------------------------------------------------


def build_variant(policy: str) -> AnalysisProfile:
    """The analysis under one specification policy (:data:`POLICIES`)."""
    if policy not in POLICIES:
        raise ValueError(f"unknown policy {policy!r} (known: {POLICIES})")
    return AnalysisProfile(
        name=f"italian-palatalization ({policy})",
        inventory=load_inventory(),
        patches=PATCHES,
        underspec=UNDERSPEC,
        rules=RULES,
        inventory_json=inventory_path(_INVENTORY_NAME),
        target_stems=TARGET_STEMS,
        palatalization_rule_names=PALATALIZATION_RULE_NAMES,
        expected_firings=expected_firings,
        fallback_predict=FallbackPredictor(policy),
        ur_builder=UrBuilder(policy),
        cluster_key=cluster_key,
        normalisers=ITALIAN_NORMALISERS,
        categoriser=categorise,
        row_factory=row_factory,
        tokenizer=ITALIAN_TOKENIZER,
    )


def build_profile() -> AnalysisProfile:
    """The default analysis: Gorman & Winstead's lexical specification."""
    return build_variant(LEXICAL)


def _all_underspecified() -> AnalysisProfile:
    return build_variant(ALL_UNDERSPECIFIED)


def _stress() -> AnalysisProfile:
    return build_variant(STRESS)


PROFILE_VARIANTS: Final[Mapping[str, Callable[[], AnalysisProfile]]] = {
    LEXICAL: build_profile,
    ALL_UNDERSPECIFIED: _all_underspecified,
    STRESS: _stress,
}
"""The variants ``scripts/validate_lang_rules.py --variant`` accepts."""


# ---------------------------------------------------------------------------
# Summary across variants
# ---------------------------------------------------------------------------


def _rate(part: int, whole: int) -> str:
    return f"{100 * part / whole:6.2f}%" if whole else "     - "


@dataclass(frozen=True, slots=True)
class _Tally:
    n: int
    boolean: int
    ipa: int

    def cells(self) -> str:
        return (
            f"{self.n:6d} {_rate(self.boolean, self.n)} "
            f"{_rate(self.ipa, self.n)}"
        )


def _tally(results: Sequence[RowResult]) -> _Tally:
    return _Tally(
        n=len(results),
        boolean=sum(r.boolean_match for r in results),
        ipa=sum(r.ipa_match for r in results),
    )


def _grouped(
    report: ValidationReport, key: Callable[[RowResult], str]
) -> dict[str, list[RowResult]]:
    groups: dict[str, list[RowResult]] = defaultdict(list)
    for result in report.results:
        groups[key(result)].append(result)
    return groups


def _cell_key(result: RowResult) -> str:
    return f"{result.row.stem_final}+{result.row.opportunity}"


def _class_key(result: RowResult) -> str:
    return f"+{result.row.opportunity} {row_class(result.row)}"


def _table(
    title: str,
    reports: Mapping[str, ValidationReport],
    key: Callable[[RowResult], str],
) -> list[str]:
    """One line per key: N, then boolean and IPA rates per variant (the
    variants run on the same rows, so N is shared)."""
    names = list(reports)
    grouped = {name: _grouped(reports[name], key) for name in names}
    keys = sorted({k for groups in grouped.values() for k in groups})
    width = max([len(title), *(len(k) for k in keys)]) + 2
    # Each variant's column is wide enough for its name and two rates.
    columns = {name: max(15, len(name)) for name in names}
    head = "".join(f" | {name:^{columns[name]}}" for name in names)
    sub = "".join(
        f" | {'bool':>7} {'IPA':>7}".ljust(columns[name] + 3) for name in names
    )
    lines = [
        f"{title:<{width}}       N{head}".rstrip(),
        f"{'':<{width}}        {sub}".rstrip(),
    ]
    for k in keys:
        first = grouped[names[0]].get(k, [])
        cells = ""
        for name in names:
            tally = _tally(grouped[name].get(k, []))
            cell = (
                f" | {_rate(tally.boolean, tally.n)} "
                f"{_rate(tally.ipa, tally.n)}"
            )
            cells += cell.ljust(columns[name] + 3)
        lines.append(f"{k:<{width}}{len(first):8,}{cells}".rstrip())
    return lines


def _examples(results: Sequence[RowResult], limit: int) -> list[str]:
    """Up to *limit* distinct forms (amico N and ADJ give one line)."""
    lines: list[str] = []
    seen: set[tuple[str, str, str]] = set()
    for result in results:
        row = result.row
        key = (row.lemma, row.plural, row.opportunity)
        if key in seen:
            continue
        seen.add(key)
        lines.append(
            f"      {row.lemma} → {row.plural} (+{row.opportunity}): "
            f"predicted {result.predicted_sr or '-'}, "
            f"attested {result.attested_first or '-'}"
        )
        if len(lines) == limit:
            break
    return lines


def _root_count(results: Sequence[RowResult]) -> int:
    return len(
        {
            row.root_id or row.lemma
            for row in (r.row for r in results)
            if isinstance(row, ItalianRow)
        }
    )


def _height_effect(name: str, report: ValidationReport) -> str:
    """-e rows whose dorsal the variant underspecifies, and how many of
    them the rules palatalize (the Height effect, tested without
    circularity: the policy never looks at the -e cell)."""
    policy = name if name in POLICIES else LEXICAL
    open_e = [
        r
        for r in report.results
        if r.row.opportunity == "e" and is_underspecified(r.row, policy)
    ]
    fired = [r for r in open_e if r.boolean_predicted]
    dorsal = sum(cast(ItalianRow, r.row).target_class == DORSAL for r in fired)
    return (
        f"  {name:20} -e rows with an underspecified dorsal: "
        f"{len(open_e):5,} ({_root_count(open_e):,} roots); palatalized "
        f"by the rules: {len(fired)} (dorsal targets: {dorsal}; NDE "
        f"stems, by their root vocoid: {len(fired) - dorsal})"
    )


def render_summary(
    reports: Mapping[str, ValidationReport],
    coverage: Coverage,
    *,
    sanity: ValidationReport | None = None,
    sanity_coverage: Coverage | None = None,
    examples: int = 4,
) -> str:
    """Boolean and IPA match rates per cell and per class for every
    variant side by side, the Italian mismatch types per variant with
    examples, and the coverage of the cell data."""
    lines = [
        "Italian velar palatalization: provisional LP analysis",
        "(phonology.analyses.italian_palatalization; "
        "docs/italian/phonology.md)",
        "",
        "Rows: in-domain, aligned, unfiltered cell rows of dorsal targets "
        "and NDE classes before -i/-e (τ = dictionary: every Wiktionary "
        "form), lemma transcription required.",
        f"In scope: {coverage.in_scope:,}; without a lemma transcription "
        f"(not validated): {coverage.no_ipa:,}; validated: "
        f"{coverage.validated:,}.",
    ]
    for name, (total, with_ipa) in coverage.by_class.items():
        lines.append(f"  {name:12} {with_ipa:6,} of {total:6,} validated")
    lines.append(
        "Attested plural IPA: the lemma transcription with the stem edge "
        "the spelling shows (module docstring), so the IPA rate tests the "
        "derived segments; a row whose UR cannot be built is scored by the "
        "policy (fallback) and counts as an IPA match iff its boolean "
        "matches (validator convention)."
    )
    lines.append("")
    lines.append("=" * 78)
    lines.append(
        "Overall (bool = boolean mutation match, IPA = surface match)"
    )
    lines.append("=" * 78)
    for name, report in reports.items():
        failed = sum(r.ur is None for r in report.results)
        lines.append(
            f"  {name:20} {_tally(report.results).cells()}   "
            f"UR not built: {failed}"
        )
    lines.append("")
    lines.append("=" * 78)
    lines.append(
        "Height effect: the trigger is /i j/, whatever the specification"
    )
    lines.append("=" * 78)
    lines.extend(
        _height_effect(name, report) for name, report in reports.items()
    )
    lines.append("")
    lines.append("=" * 78)
    lines.append("Per cell (stem_final + ending)")
    lines.append("=" * 78)
    lines.extend(_table("cell", reports, _cell_key))
    lines.append("")
    lines.append("=" * 78)
    lines.append("Per class (+ending, class of the root)")
    lines.append("=" * 78)
    lines.extend(_table("class", reports, _class_key))
    for name, report in reports.items():
        lines.append("")
        lines.append("=" * 78)
        lines.append(f"Mismatch types: {name}")
        lines.append("=" * 78)
        typed: dict[str, list[RowResult]] = defaultdict(list)
        for result in report.results:
            kind = mismatch_class(result)
            if kind is not None:
                typed[kind].append(result)
        if not typed:
            lines.append("  none")
        for kind in sorted(typed, key=lambda k: (-len(typed[k]), k)):
            items = typed[kind]
            lines.append(
                f"  [{kind}]  n = {len(items)} rows, "
                f"{_root_count(items)} roots"
            )
            lines.extend(_examples(items, examples))
    if sanity is not None:
        lines.append("")
        lines.append("=" * 78)
        lines.append(
            "Sanity rows: coronal targets and lexical palatals (no rule "
            "may touch them; the data code coronal plurals as unchanged "
            "by construction)"
        )
        lines.append("=" * 78)
        if sanity_coverage is not None:
            for name, (total, with_ipa) in sanity_coverage.by_class.items():
                if name in (CORONAL, FRONTSTEM):
                    lines.append(
                        f"  {name:12} {with_ipa:6,} of {total:6,} validated"
                    )
        by_class = _grouped(sanity, _class_key)
        for k in sorted(by_class):
            lines.append(f"  {k:30} {_tally(by_class[k]).cells()}")
        misses = [r for r in sanity.results if mismatch_class(r) is not None]
        lines.append(f"  mismatches: {len(misses)}")
        lines.extend(_examples(misses, examples))
    return "\n".join(lines) + "\n"


def trace(
    rows: Sequence[LexRow], lemmas: Sequence[str], profile: AnalysisProfile
) -> str:
    """Paper-style derivation tables of the rows of *lemmas*."""
    out: list[str] = []
    inv_base = profile.inventory.base_segments()
    for lemma in lemmas:
        found = [row for row in rows if row.lemma == lemma]
        if not found:
            out.append(f"{lemma}: no validated row\n")
        for row in found:
            result = evaluate_row(row, profile)
            out.append(render_trace(TraceResult(row, result), inv_base))
    return "\n".join(out)


def main(argv: Sequence[str] | None = None) -> int:
    """Write the cross-variant summary (and optionally print traces)::

        PYTHONPATH=src python -m phonology.analyses.italian_palatalization

    writes ``analysis/italian/italian_rule_validation_summary.txt``. The
    per-variant reports and mismatch CSVs come from
    ``scripts/validate_lang_rules.py --lang it [--variant NAME]``.
    """
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument("--root", type=Path, default=PROJECT_ROOT)
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="summary path (default: analysis/italian/"
        "italian_rule_validation_summary.txt under --root)",
    )
    parser.add_argument(
        "--trace",
        nargs="*",
        default=None,
        metavar="LEMMA",
        help="print the derivations of these lemmas instead",
    )
    parser.add_argument("--variant", default=LEXICAL, choices=POLICIES)
    args = parser.parse_args(argv)
    root: Path = args.root
    rows, coverage = load_rows(root)
    if args.trace is not None:
        print(trace(rows, args.trace, build_variant(args.variant)))
        return 0
    reports = {
        name: run_validation(factory(), rows)
        for name, factory in PROFILE_VARIANTS.items()
    }
    sanity_rows, sanity_coverage = load_rows(
        root, classes=(CORONAL, FRONTSTEM)
    )
    sanity = run_validation(build_profile(), sanity_rows)
    text = render_summary(
        reports,
        coverage,
        sanity=sanity,
        sanity_coverage=sanity_coverage,
    )
    out = args.out or (
        LANGUAGE.analysis_dir(root)
        / LANGUAGE.file_name("rule_validation_summary", ".txt")
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    print(text, end="")
    print(f"\nSummary written to {out}", file=sys.stderr)
    return 0


__all__: Final[tuple[str, ...]] = (
    "ALL_CLASSES",
    "ALL_UNDERSPECIFIED",
    "ANTERIOR_SIBILANT",
    "CORONAL",
    "DEFAULT_CLASSES",
    "DORSAL",
    "DORSAL_CLASSES",
    "DORSAL_DEFAULT",
    "DORSAL_PAL",
    "FRONTSTEM",
    "GLIDE_ABSORPTION",
    "ITALIAN_NORMALISERS",
    "ITALIAN_TOKENIZER",
    "LANGUAGE",
    "LENGTH_DEFAULT",
    "LEXICAL",
    "LONG_POSTALVEOLAR_FRICATIVE",
    "NDE_OCCHIO",
    "NDE_PSICHE",
    "OBSTRUENT_CLASS",
    "PALATALIZATION_RULES",
    "PALATALIZATION_RULE_NAMES",
    "PATCHES",
    "POLICIES",
    "POSTALVEOLAR_PLACE",
    "POSTALVEOLAR_STRIDENT",
    "PROFILE_VARIANTS",
    "RULES",
    "SC_COALESCENCE",
    "SC_VELAR",
    "STRESS",
    "S_DELETION",
    "TARGET_STEMS",
    "TRIGGER_CLASS",
    "UNDERSPEC",
    "UNDERSPEC_FOR",
    "Coverage",
    "FallbackPredictor",
    "ItalianRow",
    "RootEvidence",
    "UrBuilder",
    "build_profile",
    "build_ur",
    "build_variant",
    "categorise",
    "cluster_key",
    "expected_firings",
    "is_underspecified",
    "iter_rows",
    "load_inventory",
    "load_rows",
    "main",
    "merge_geminates",
    "mismatch_class",
    "plural_ipa",
    "render_summary",
    "row_class",
    "row_factory",
    "tokens",
    "trace",
    "ur_labels",
)


if __name__ == "__main__":
    raise SystemExit(main())
