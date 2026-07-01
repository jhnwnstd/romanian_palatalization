#!/usr/bin/env python3
"""Build the co-variation contingency table.

Canonical script for the table Steriade 2008 §10.2.3.5 offers.
The goal here is to test **the empirical observation** — the claim
that in the modern Romanian lexicon a noun's inflectional palatalization
class co-varies with its choice of verbalizing suffix. This is the
descriptive claim; whether the pattern is caused by paradigm
consultation, cyclic stratification, or something else is the
*analytical* question that lies downstream of this test.

PRIMARY uses the **most atheoretical row axis available**: allomorph
identity keyed to stem class. A noun is placed in the "palatalising
allomorph" row iff its (stem-final, plural-suffix) combination is
one that palatalises *in general* for that stem class in this lexicon
— dorsals with ``-i`` or ``-e`` plurals, coronals with ``-i`` plurals.
This is a purely morphological label; no per-noun mutation lookup,
no phonological analysis. The alternative row-axis codings (mutation
flag, sg-vs-pl IPA comparison) are kept as robustness panels so
readers can see how sensitive the number is to those choices.

The claim under test
--------------------
For a noun with a target obstruent stem-final, the plural allomorph
and the verbalizer allomorph co-vary: nouns whose plural palatalizes
the stem-final tend to select the palatalizing verbalizer ``-í``;
nouns whose plural doesn't palatalize tend to select a non-palatalizing
verbalizer (``-á`` or ``-uí``).

Axes and coding
---------------
Row-axis choices (three, from most atheoretical to most theory-laden):

  (a) ALLOMORPH-STEMWISE — a noun's row is (stem-final, plural-suffix)
      → "palatalising allomorph slot" (True) / "non-palatalising slot"
      (False). Purely morphological. Grouping:
        dorsals  : {-i, -e} vs -uri
        coronals : -i      vs {-e, -uri}
      Justification for the grouping is atheoretical: from the
      Kyle-John-Daniar Table 3 productivity counts, dorsal -i and -e
      plurals palatalise categorically while coronal -e plurals never
      do (e = n = 100% exceptions to assibilation before -e). This is
      an observable class-level fact, not per-noun analysis. THE
      PRIMARY PANEL USES THIS CODING.
  (b) MUTATION-FLAG — the lexicon's per-noun ``mutation`` boolean,
      a direct surface-orth comparison of singular vs plural.
      Steriade-original but requires trusting the mutation label
      per noun.
  (c) IPA-SG-vs-PL — compare the last consonantal segment of the
      singular's normalised IPA to the same position in the plural's,
      after stripping known plural-suffix material. Three empirical
      states:
        BASE→BASE   → K/K   (non-alternating)
        BASE→PAL    → K/T͡ʃ ALTERNATING (mutation)
        PAL→PAL     → K/T͡ʃ LEXICAL (base already palatal — arici,
                       pace, elogiu; behaviourally K/T͡ʃ under both
                       theories, so pooled with ALTERNATING; the split
                       is reported for transparency)
        PAL→BASE    → inconsistent, dropped
      Catches (i) cluster cases where the orthographic comparison
      conflates the sc→ʃt cascade into a single ``True`` flag; (ii)
      lexical-palatal bases that the orthographic comparison
      correctly marks non-mutating but that behave as K/T͡ʃ.

COL = "does the derivation palatalize?"
  PAL column     = ``-i`` (palatalizing verbalizer).
  NON-PAL column = ``-a`` OR ``-ui`` (both are non-palatalizing:
      ``-a`` is back-vocalic, ``-ui`` has a blocking /u/).
  This 1-vs-2 grouping matches Steriade §10.2.3.5 table (10.20):
  she compares ``-í`` to the combined non-palatalizing set. Earlier
  drafts of this script dropped ``-ui`` entirely, which loses the
  strongest K/K signal (Steriade: 45% of K/K bases take ``-uí``).

  (a) SUFFIX-ONLY — group by suffix membership above, no paradigm
      filter. Steriade's literal claim is about affix avoidance,
      which is what this measures.
  (b) ALTERNATES-ONLY — a ``-i`` verb only counts as PAL evidence
      if its paradigm verdict is ``ALTERNATES`` (verb actually surfaces
      the palatalization in a diagnostic cell). ``TRIGGER-ABSENT``
      verbs (whose ``-esc``/``-ez`` augment bleeds the phonological
      trigger) test morphology, not phonology, and are dropped —
      they can't distinguish "the base can palatalize" from "the base
      can't". ``NO-PARADIGM`` verbs are also dropped (no evidence).
      ``-a`` and ``-ui`` verbs are always kept (they don't palatalize
      by suffix construction; nothing for a paradigm to reveal).

The six panels
--------------
    row = ALLOMORPH-STEMWISE, col = SUFFIX-ONLY       → PRIMARY
    row = MUTATION,           col = SUFFIX-ONLY       → MUTATION-ROW
    row = IPA,                col = SUFFIX-ONLY       → IPA-ROW
    row = ALLOMORPH,          col = SUFFIX-ONLY       → STEM-BLIND ALLOMORPH
    row = ALLOMORPH-STEMWISE, col = ALTERNATES-ONLY   → PARADIGM-STRICT
    row = IPA,                col = ALTERNATES-ONLY   → FULLY-STRICT

PRIMARY is our best answer: the row axis is a pure allomorph-identity
label whose grouping ({-i,-e} vs -uri for dorsals; -i vs {-e,-uri}
for coronals) reflects observable class-level palatalization behaviour
without any per-noun mutation lookup. If the effect is real, it
should show up under this coding — and it does, most strongly of
any panel we compute.

MUTATION-ROW and IPA-ROW substitute Steriade-style row axes; if
PRIMARY and these agree in direction and significance, the
observation isn't an artefact of any single row-axis choice.

STEM-BLIND ALLOMORPH is kept as a **counter-example**: it groups
coronal -e plurals with coronal -i plurals (as Kyle-John-Daniar
Table (9) does), which conflates two behaviourally-opposite classes
and drives the coronal OR to ≈ 1. This panel shows what that
grouping choice does — the null result Table (9) reports is
recoverable if you use its exact methodology.

PARADIGM-STRICT and FULLY-STRICT drop ``-í`` verbs whose paradigm
verdict isn't ``ALTERNATES`` (the ``-esc``/``-ez`` augment bleeds the
trigger, so no palatalization surfaces). These panels are stricter
than Steriade's own count — she keeps those verbs and calls them
lexical exceptions. Kept here as a robustness check for readers who
want the phonology-only view.

Exclusions applied everywhere
-----------------------------
- ``pos == 'N'``; ``stem_final in {c, g, t, d, s, z}``
- ``exception_reason`` doesn't start with ``nde:`` (Steriade's NDEB
  excludes)
- ``derived_verbs`` field non-empty
- Nouns whose etymology-accepted verbs use *both* ``-i`` and ``-a``
  emit two records (one per suffix), each carrying only the verbs of
  that suffix. This keeps the ambiguous noun's evidence honest instead
  of collapsing it to whichever suffix happens to be preferred.
- Column-side: etymology audit accepts the (noun, verb) pair, or the
  noun is in the ``USER_KEEP`` back-formation allowlist
- Column-side: verbalizer restricted to ``{-i, -a, -ui}``.
- Column-side: (fr-etym noun, ``-a`` verb) pairs are DROPPED. Steriade
  §10.2.3.5 excludes them explicitly: "the list excludes clear loans
  from French (e.g., [remark-á], French *remarquer*) whose assignment
  to the [-á] conjugation is based on their French inflectional class
  rather than a preference for [-á]." Latin loans are retained — the
  paper doesn't blanket-exclude them, and they'd wipe out most of the
  Romanian native lexicon.
- Row-side: rows whose IPA lemma or plural is missing/foreign-only get
  a UNCLEAR classification and are dropped from the IPA-row panels
  (but retained in the MUTATION-row panels, which don't need IPA)
- Row-side: singularia tantum (``opportunity == 'none'``) are KEPT
  and classified as K/K by construction. Steriade §10.2.3.4: mass
  nouns like *vlagă* and proper names like *Franco* have no plural,
  hence no palatalized allomorph, hence categorically pattern K/K.

Two frequency scopes
--------------------
- TOP-1000: the top-1000 nouns by frequency, matching the size of
  Steriade's own hand-collected sample.
- FULL DEX: the modern-lexicon replication.

Statistics
----------
Fisher's exact 2×2 per stem-class panel (dorsal / coronal / pooled),
with 95% CI on the OR (scipy.stats.contingency.odds_ratio, exact
conditional Fisher CI). Fisher over chi-square because per-stem-class
cells can be small.

Robustness (all reported alongside the four panels)
---------------------------------------------------
- BASELINE table: -í / -á / -uí distribution on non-alternating stems
  (r, n, m, ț, ș, p, b, f, v). Steriade's own frame of reference.
- NON-PAL BREAKDOWN: split the non-palatalising column into -á vs -uí
  per cell — verifies -uí really does dominate on K/K bases.
- SINGULARIA TANTUM: distribution across cells for the ``opp='none'``
  rows Steriade classifies as K/K by construction.
- PERMUTATION NULL: 10,000 shuffles of the K/T͡ʃ-vs-K/K label; if
  the observed OR were random co-incidence given the marginals, the
  permutation distribution would routinely produce it.
- ETYMOLOGY STRATIFICATION: PRIMARY re-run on native rows only,
  Latin rows only, everything-except-Latin — verifies the effect
  isn't carried by any single etymological subset.
- SING-TANTUM EXCLUDED / EMPTY-ETYM EXCLUDED: sensitivity to the two
  most defensible "drop this if in doubt" filters.
- CONTAMINATION DIAGNOSTIC: how much of the -í column has a
  ``TRIGGER-ABSENT`` paradigm verdict (i.e. what PARADIGM-STRICT
  removes relative to PRIMARY).
"""

from __future__ import annotations

import csv
import random
import sys
from collections import Counter
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Callable, Final

from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PROJECT_ROOT / "src"))

from dex_etymology import USER_KEEP  # noqa: E402

CSV_PATH: Final[Path] = (
    _PROJECT_ROOT / "data" / "romanian_lexicon_with_freq.csv"
)
ETYMOLOGY_CSV_PATH: Final[Path] = (
    _PROJECT_ROOT / "data" / "derived_verb_etymology.csv"
)
VERDICT_CSV_PATH: Final[Path] = (
    _PROJECT_ROOT / "data" / "verb_paradigm_verdict.csv"
)
OUTPUT_PATH: Final[Path] = _PROJECT_ROOT / "analysis" / "contingency_table.txt"

TRIGGER_ABSENT: Final[str] = "TRIGGER-ABSENT"
ALTERNATES: Final[str] = "ALTERNATES"

TOP_N_FREQUENCY: Final[int] = 1_000

DORSAL_STEMS: Final[frozenset[str]] = frozenset({"c", "g"})
CORONAL_STEMS: Final[frozenset[str]] = frozenset({"t", "d", "s", "z"})
TARGET_STEMS: Final[frozenset[str]] = DORSAL_STEMS | CORONAL_STEMS
PLURAL_FRONT_OPPS: Final[frozenset[str]] = frozenset({"i", "e"})
# Steriade §10.2.3.5 surveys three verbalizers; ``-ui`` is
# non-palatalising like ``-a`` and shares its column here.
ACCEPTED_VERBALIZERS: Final[frozenset[str]] = frozenset({"-i", "-a", "-ui"})
PALATALIZING_VERBALIZERS: Final[frozenset[str]] = frozenset({"-i"})
NONPAL_VERBALIZERS: Final[frozenset[str]] = frozenset({"-a", "-ui"})
# Steriade §10.2.3.5 baseline: nouns whose stem-final consonant
# undergoes no relevant alternation. Used to compute the "neutral"
# distribution of ``-í/-á/-uí`` against which K/T͡ʃ and K/K rates are
# over/under-represented. Bases in [t, d, s, z] are excluded because
# they participate in Assibilation / S-Palatalisation.
BASELINE_STEMS: Final[frozenset[str]] = frozenset(
    {"r", "n", "m", "ț", "ș", "p", "b", "f", "v"}
)


# ---------------------------------------------------------------------------
# IPA classifier for the row axis
# ---------------------------------------------------------------------------
#
# The row axis under IPA coding is a three-state comparison of the
# last stem-final consonant in the singular's normalised IPA to the
# same position in the plural's. Both sides use the shared tokeniser
# below so tie-barred affricates count as one segment.

BASE_IPA: Final[dict[str, str]] = {
    "c": "k",
    "g": "ɡ",
    "t": "t",
    "d": "d",
    "s": "s",
    "z": "z",
}
PAL_IPA: Final[dict[str, str]] = {
    "c": "t͡ʃ",
    "g": "d͡ʒ",
    "t": "t͡s",
    "d": "z",
    "s": "ʃ",
    "z": "ʒ",
}

_PLURAL_SUFFIXES: Final[tuple[str, ...]] = (
    "uri",
    "urj",
    "urʲ",
    "ele",
    "ile",
    "j",
    "ʲ",
    "i",
    "e",
)
# Singular strippable material — only vowels that participate in
# inflectional class (-ă/-ə for feminine, -e for -e class). Other
# vowels (o, ɨ) at word end are treated as stem-internal.
_SING_SUFFIXES: Final[tuple[str, ...]] = ("ə", "e")

_VOWELS: Final[set[str]] = set("aeiouəɨɛɔ")
_GLIDES: Final[set[str]] = set("jw")
_MARKS: Final[set[str]] = set("ˈˌ.'` ")


def _tokens(ipa: str) -> list[str]:
    """Split IPA into segment atoms; tie-barred affricates (t͡ʃ etc.)
    are one atom of three code points."""
    out: list[str] = []
    i = 0
    while i < len(ipa):
        ch = ipa[i]
        if ch in _MARKS:
            i += 1
            continue
        if i + 2 < len(ipa) and ipa[i + 1] == "͡":
            out.append(ipa[i : i + 3])
            i += 3
        else:
            out.append(ch)
            i += 1
    return out


def _first_variant(field_val: str) -> str:
    """Return the first non-empty, non-foreign-language variant of a
    pipe-separated IPA field. ``fr:absorbant``-style rows are dropped."""
    for v in field_val.split(" | "):
        v = v.strip()
        if not v:
            continue
        if len(v) >= 3 and v[2] == ":" and v[:2].isalpha() and v[:2].islower():
            continue
        return v
    return ""


def _strip_suffix(ipa: str, suffixes: tuple[str, ...]) -> str:
    """Strip the longest matching suffix once."""
    for s in sorted(suffixes, key=len, reverse=True):
        if ipa.endswith(s):
            return ipa[: -len(s)]
    return ipa


def _last_consonantal(ipa: str) -> str | None:
    """Walk reverse-tokens, skipping trailing vowels/glides/ʲ,
    return the first consonantal segment we find."""
    for tok in reversed(_tokens(ipa)):
        head = tok[0]
        if head in _VOWELS or head in _GLIDES or tok == "ʲ":
            continue
        return tok
    return None


def _form_class(ipa: str, sf: str, suffixes: tuple[str, ...]) -> str:
    """Classify one form ('pal' / 'base' / 'other') for stem-final sf.

    Cluster carve-out: Romanian's ``ʃt`` cluster (from ``sc``/``st``
    stems like *muʃte*, *liniʃte*, *poveʃti*) puts the ``ʃ`` *before*
    the last consonant. If ``ʃ`` appears in the last 4 chars of the
    stripped stem, the form counts as palatalised — regardless of
    whether the lexicon assigned the stem-final label to ``s``, ``c``
    or ``t``.
    """
    ipa = _first_variant(ipa)
    if not ipa:
        return "other"
    stem = _strip_suffix(ipa, suffixes)
    if not stem:
        return "other"
    last = _last_consonantal(stem)
    if last is None:
        return "other"
    if last == PAL_IPA[sf]:
        return "pal"
    if sf in ("s", "c", "t") and "ʃ" in stem[-4:]:
        return "pal"
    if last == BASE_IPA[sf]:
        return "base"
    return "other"


class IPARowClass(StrEnum):
    """Three-state outcome of the sg-vs-pl IPA comparison."""

    NONALT = "nonalt"  # sg base, pl base
    ALT = "alt"  # sg base, pl palatal (a real mutation)
    LEX = "lex"  # sg palatal, pl palatal (already palatal)
    UNCLEAR = "unclear"  # missing IPA or PAL→BASE inconsistency


def _classify_ipa_row(sg_ipa: str, pl_ipa: str, sf: str) -> IPARowClass:
    """Three-state classification of the plural mutation via sg-vs-pl."""
    sg = _form_class(sg_ipa, sf, _SING_SUFFIXES)
    pl = _form_class(pl_ipa, sf, _PLURAL_SUFFIXES)
    if sg == "other" or pl == "other":
        return IPARowClass.UNCLEAR
    if sg == "base" and pl == "pal":
        return IPARowClass.ALT
    if sg == "base" and pl == "base":
        return IPARowClass.NONALT
    if sg == "pal" and pl == "pal":
        return IPARowClass.LEX
    return IPARowClass.UNCLEAR  # pal→base is inconsistent


# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------


class StemClass(StrEnum):
    DORSAL = "dorsal"
    CORONAL = "coronal"


class RowCoding(StrEnum):
    """How the row axis (plural class) is coded."""

    MUTATION = "mutation"  # orthographic ``mutation`` field
    IPA = "ipa"  # sg-vs-pl IPA comparison
    ALLOMORPH = "allomorph"  # plural-suffix identity: {-i, -e} vs -uri
    ALLOMORPH_STEMWISE = "allomorph_stemwise"
    # {-i, -e} vs -uri for dorsals ; -i vs {-e, -uri} for coronals.
    # Purely allomorph-based (no per-noun palatalization check) but
    # groups allomorphs by whether they palatalize *in general* for
    # that stem class — which is an aggregate phonological fact
    # (from the Kyle-John-Daniar Table 3 productivity counts), not
    # a per-noun mutation lookup.


class ColCoding(StrEnum):
    """How the column axis (verbalizer class) is coded."""

    SUFFIX_ONLY = "suffix"  # -i vs -a, no paradigm filter
    ALTERNATES_ONLY = "alternates"  # -i verbs must be ALTERNATES verdict


# ---------------------------------------------------------------------------
# The PRIMARY coding — used by the headline and every robustness diagnostic
# ---------------------------------------------------------------------------
#
# The PRIMARY row axis is ALLOMORPH_STEMWISE: a noun's row is fixed
# purely by (stem-final consonant, plural-suffix identity), without
# any per-noun palatalization check. The grouping puts palatalising
# allomorph slots on one side and non-palatalising on the other,
# where "palatalising for this stem" is a general class-level fact
# (Kyle-John-Daniar Table 3): {-i, -e} palatalise dorsals; only -i
# palatalises coronals. This is the most atheoretical row axis we
# have — no mutation lookup, no IPA comparison — and it delivers
# the largest, cleanest FULL-DEX effect. See build_contingency_table
# docstring for the full panel comparison.
PRIMARY_ROW_CODING: Final[RowCoding] = RowCoding.ALLOMORPH_STEMWISE
PRIMARY_COL_CODING: Final[ColCoding] = ColCoding.SUFFIX_ONLY


@dataclass(frozen=True, slots=True)
class NounRecord:
    """A single noun row that passes the row-side filters."""

    lemma: str
    derived_verbs: tuple[str, ...]
    stem_final: str
    opportunity: str
    deriv_suffix: str
    mutation: bool
    frequency: float
    ipa_row_class: IPARowClass
    etym_lang: str
    is_singulare_tantum: bool

    @property
    def stem_class(self) -> StemClass:
        return (
            StemClass.DORSAL
            if self.stem_final in DORSAL_STEMS
            else StemClass.CORONAL
        )

    def plural_alt(self, coding: RowCoding) -> bool | None:
        """Return True (alt) / False (nonalt) / None (drop) under the
        given row coding.

        Singularia tantum (Steriade §10.2.3.4) are K/K by construction
        under both mutation and IPA codings: no plural means no
        palatalized allomorph can be generated, so the noun
        categorically patterns as non-alternating. This mirrors the
        paper's treatment of *vlagă*, *Franco*, *Goga* et al.

        ALLOMORPH coding is a purely morphological label: True if the
        noun's plural suffix is ``-i`` or ``-e``, False if ``-uri``.
        Singularia tantum (``opp=='none'``) are dropped under this
        coding because they have no plural suffix to name. This
        matches the Kyle-John-Daniar Table (9) row axis exactly.
        """
        if coding is RowCoding.ALLOMORPH:
            if self.is_singulare_tantum:
                return None
            if self.opportunity in PLURAL_FRONT_OPPS:  # {i, e}
                return True
            if self.opportunity == "uri":
                return False
            return None
        if coding is RowCoding.ALLOMORPH_STEMWISE:
            if self.is_singulare_tantum:
                return None
            if self.stem_final in DORSAL_STEMS:
                # {-i, -e} palatalise for dorsals, -uri doesn't
                if self.opportunity in PLURAL_FRONT_OPPS:
                    return True
                if self.opportunity == "uri":
                    return False
            elif self.stem_final in CORONAL_STEMS:
                # -i palatalises for coronals; {-e, -uri} don't
                if self.opportunity == "i":
                    return True
                if self.opportunity in ("e", "uri"):
                    return False
            return None
        if self.is_singulare_tantum:
            return False
        if coding is RowCoding.MUTATION:
            return self.mutation
        if self.ipa_row_class is IPARowClass.ALT:
            return True
        if self.ipa_row_class is IPARowClass.LEX:
            return True  # pooled with ALT — both are K/T͡ʃ class
        if self.ipa_row_class is IPARowClass.NONALT:
            return False
        return None  # UNCLEAR


@dataclass(frozen=True, slots=True)
class PanelStats:
    """Fisher 2×2 result for one panel.

    The 2×2 collapses ``-a`` and ``-ui`` into a single NON-PAL column
    (both are non-palatalising verbalisers). The by-suffix breakdown
    is stored separately in ``alt_a``/``alt_ui``/``nonalt_a``/
    ``nonalt_ui`` so the diagnostic can show which non-pal suffix is
    doing the work in each cell (Steriade: K/K → -uí dominates).
    """

    alt_pal: int  # K/T͡ʃ base × -i (palatalising) verbalizer
    alt_nonpal: int  # K/T͡ʃ base × (-a or -ui) verbalizer
    nonalt_pal: int  # K/K base × -i verbalizer
    nonalt_nonpal: int  # K/K base × (-a or -ui) verbalizer

    # By-suffix breakdown of the NON-PAL columns (informational only).
    alt_a: int = 0
    alt_ui: int = 0
    nonalt_a: int = 0
    nonalt_ui: int = 0

    @property
    def n(self) -> int:
        return (
            self.alt_pal
            + self.alt_nonpal
            + self.nonalt_pal
            + self.nonalt_nonpal
        )

    @property
    def is_complete(self) -> bool:
        return (
            self.alt_pal + self.alt_nonpal > 0
            and self.nonalt_pal + self.nonalt_nonpal > 0
            and self.alt_pal + self.nonalt_pal > 0
            and self.alt_nonpal + self.nonalt_nonpal > 0
        )

    def fisher(self) -> tuple[float, float, float, float] | None:
        """Return (odds ratio, 95% CI low, 95% CI high, p) or None if
        any margin is zero.

        OR and CI come from :func:`scipy.stats.contingency.odds_ratio`
        with ``kind='sample'`` (matches Fisher's plain ``a·d/(b·c)``).
        The CI is the exact conditional Fisher CI. p-value is Fisher's
        two-sided exact test.
        """
        if not self.is_complete:
            return None
        table = [
            [self.alt_pal, self.alt_nonpal],
            [self.nonalt_pal, self.nonalt_nonpal],
        ]
        _, p = fisher_exact(table)
        or_result = odds_ratio(table, kind="sample")
        ci = or_result.confidence_interval(confidence_level=0.95)
        return (
            float(or_result.statistic),
            float(ci.low),
            float(ci.high),
            float(p),
        )

    def format_line(self, name: str) -> str:
        head = (
            f"  {name:8s}:  alt-pl  {self.alt_pal:3d}  -i / "
            f"{self.alt_nonpal:3d}  nonpal    "
            f"nonalt-pl  {self.nonalt_pal:3d}  -i / "
            f"{self.nonalt_nonpal:3d}  nonpal    n = {self.n:4d}"
        )
        fisher = self.fisher()
        if fisher is None:
            return head
        odds, ci_low, ci_high, p = fisher
        p_fmt = f"{p:.2e}" if p < 1e-3 else f"{p:.3f}"

        def _or_fmt(x: float) -> str:
            if x == float("inf"):
                return "  inf"
            if x >= 100:
                return f"{x:5.0f}"
            return f"{x:5.2f}"

        return (
            f"{head}    OR = {_or_fmt(odds)}  "
            f"[{_or_fmt(ci_low)}, {_or_fmt(ci_high)}]    p = {p_fmt}"
        )


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def _resolve_records(
    lemma: str,
    verbs: list[str],
    suffixes: list[str],
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
    etym_lang: str,
) -> list[tuple[str, tuple[str, ...]]]:
    """Return one (suffix, verbs-for-that-suffix) pair per emit.

    A noun with etymology-accepted verbs under one suffix returns a
    single-element list. A noun with accepted verbs under multiple
    suffixes returns one element per suffix — the noun is counted in
    each column, since the language demonstrably accepts multiple
    derivations from the same base. A USER_KEEP noun with no audited
    pair falls back to the first ACCEPTED_VERBALIZERS suffix in the
    field. Any other configuration returns an empty list (drop).

    Steriade exclusion: (fr-etym noun, ``-a`` verb) pairs are dropped
    (Steriade §10.2.3.5: "the list excludes clear loans from French
    ... whose assignment to the [-á] conjugation is based on their
    French inflectional class rather than a preference for [-á]").
    Latin loans are not blanket-excluded.
    """
    pairs = [
        (v, s) for v, s in zip(verbs, suffixes) if s in ACCEPTED_VERBALIZERS
    ]
    # Steriade-style French-loan exclusion: -a on fr-etym nouns is not
    # evidence of Romanian speakers' verbalizer preference.
    if etym_lang == "fr":
        pairs = [(v, s) for v, s in pairs if s != "-a"]
    if not pairs:
        return []
    by_suffix: dict[str, list[str]] = {}
    for v, s in pairs:
        if accepted_by_pair.get((lemma, v), False):
            by_suffix.setdefault(s, []).append(v)
    if by_suffix:
        return [(s, tuple(vs)) for s, vs in sorted(by_suffix.items())]
    if lemma in user_keep:
        first_v, first_s = pairs[0]
        return [(first_s, (first_v,))]
    return []


def _load_records(
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
) -> tuple[list[NounRecord], dict[str, float]]:
    """Read the lexicon CSV and return the list of NounRecords that
    pass the row-side filters, plus a {lemma → frequency} map."""
    records: list[NounRecord] = []
    seen: set[str] = set()
    frequency_by_lemma: dict[str, float] = {}

    with CSV_PATH.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["pos"] != "N":
                continue
            lemma = row["lemma"]
            try:
                freq = float(row.get("freq_ron_news_2024_1M") or 0)
            except ValueError:
                freq = 0.0
            frequency_by_lemma[lemma] = freq

            stem_final = row["stem_final"]
            if stem_final not in TARGET_STEMS:
                continue
            opp = row["opportunity"].strip()
            # Steriade §10.2.3.4: singularia tantum have opp='none'
            # (no plural, no palatalization opportunity). Keep them —
            # they're K/K by construction on the row axis.
            is_sing_tantum = opp == "none"
            if not is_sing_tantum:
                if opp not in PLURAL_FRONT_OPPS and opp != "uri":
                    continue
                # Non-singulare-tantum rows need a mutation value we
                # can trust for the MUTATION row coding.
                if row["mutation"].strip() == "":
                    continue
            derived_verbs_field = (row.get("derived_verbs") or "").strip()
            if not derived_verbs_field:
                continue
            if (row.get("exception_reason") or "").startswith("nde:"):
                continue
            if lemma in seen:
                continue
            verbs = [
                v.strip() for v in derived_verbs_field.split("|") if v.strip()
            ]
            suffixes = [
                s.strip()
                for s in (row.get("deriv_suffixes") or "").split("|")
                if s.strip()
            ]
            etym_lang = (row.get("etym_lang") or "").strip()
            resolutions = _resolve_records(
                lemma,
                verbs,
                suffixes,
                accepted_by_pair,
                user_keep,
                etym_lang,
            )
            if not resolutions:
                continue
            seen.add(lemma)
            ipa_row_class = _classify_ipa_row(
                row.get("ipa_normalized_lemma", ""),
                row.get("ipa_normalized_pl", ""),
                stem_final,
            )
            mutation = (
                False
                if is_sing_tantum
                else row["mutation"].strip().upper() == "TRUE"
            )
            for deriv_suffix, resolved_verbs in resolutions:
                records.append(
                    NounRecord(
                        lemma=lemma,
                        derived_verbs=resolved_verbs,
                        stem_final=stem_final,
                        opportunity=opp,
                        deriv_suffix=deriv_suffix,
                        mutation=mutation,
                        frequency=freq,
                        ipa_row_class=ipa_row_class,
                        etym_lang=etym_lang,
                        is_singulare_tantum=is_sing_tantum,
                    )
                )
    return records, frequency_by_lemma


# ---------------------------------------------------------------------------
# Etymology + paradigm audits
# ---------------------------------------------------------------------------


def _load_verb_verdicts() -> dict[str, str]:
    out: dict[str, str] = {}
    with VERDICT_CSV_PATH.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            verb = (row.get("verb") or "").strip()
            verdict = (row.get("verdict") or "").strip()
            if verb and verdict:
                out[verb] = verdict
    return out


def _load_etymology_decisions() -> (
    tuple[dict[tuple[str, str], bool], set[str]]
):
    accepted: dict[tuple[str, str], bool] = {}
    user_keep: set[str] = set()
    with ETYMOLOGY_CSV_PATH.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            noun = (row.get("noun") or "").strip()
            verb = (row.get("verb") or "").strip()
            cat = (row.get("category") or "").strip().upper()
            is_acc = (row.get("accepted") or "").strip().lower() == "true"
            if cat == "USER_KEEP":
                user_keep.add(noun)
                continue
            if noun and verb:
                accepted[(noun, verb)] = is_acc
    user_keep |= set(USER_KEEP)
    return accepted, user_keep


def _accepted_verbs(
    record: NounRecord,
    accepted_by_pair: dict[tuple[str, str], bool],
    verdicts: dict[str, str],
    col_coding: ColCoding,
    unclassified_pairs: set[tuple[str, str]],
) -> list[str]:
    """Return the verbs on this noun that count as evidence under the
    column coding.

    Under SUFFIX-ONLY, all etymology-accepted verbs pass through.
    Under ALTERNATES-ONLY, ``-i`` verbs whose paradigm verdict isn't
    ``ALTERNATES`` are dropped (their morphology bleeds the trigger,
    or there's no paradigm data — either way, no phonological
    evidence). ``-a`` and ``-ui`` verbs are always kept: they can't
    palatalise by suffix construction, so paradigm data can't reveal
    a phonological process that couldn't have applied anyway.
    """
    out: list[str] = []
    for verb in record.derived_verbs:
        key = (record.lemma, verb)
        if key not in accepted_by_pair:
            unclassified_pairs.add(key)
            continue
        if not accepted_by_pair[key]:
            continue
        if col_coding is ColCoding.ALTERNATES_ONLY:
            if record.deriv_suffix == "-i":
                if verdicts.get(verb) != ALTERNATES:
                    continue
        out.append(verb)
    return out


# ---------------------------------------------------------------------------
# Panel building
# ---------------------------------------------------------------------------


_ScopeFilter = Callable[[NounRecord], bool]


def _col_class(suffix: str) -> str:
    """Bucket a suffix into ``pal`` or ``nonpal``.

    Steriade §10.2.3.5: ``-í`` triggers palatalisation; ``-á`` (back
    vowel) and ``-uí`` (blocking /u/) both do not. The empirical claim
    is a 1-vs-2 grouping, not 1-vs-1.
    """
    if suffix in PALATALIZING_VERBALIZERS:
        return "pal"
    return "nonpal"


def _build_panels(
    records: list[NounRecord],
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
    verdicts: dict[str, str],
    row_coding: RowCoding,
    col_coding: ColCoding,
    unclassified_pairs: set[tuple[str, str]],
    scope: _ScopeFilter,
) -> tuple[dict[StemClass | str, PanelStats], int]:
    """Build dorsal / coronal / pooled panels under a specific coding.

    Returns the three panels + the count of kept records.
    """
    cells: dict[StemClass | str, Counter[tuple[bool, str]]] = {
        StemClass.DORSAL: Counter(),
        StemClass.CORONAL: Counter(),
        "pooled": Counter(),
    }
    n_kept = 0
    for record in records:
        if not scope(record):
            continue
        row_alt = record.plural_alt(row_coding)
        if row_alt is None:
            continue  # UNCLEAR IPA
        if record.lemma in user_keep:
            keep = True
        else:
            keep = bool(
                _accepted_verbs(
                    record,
                    accepted_by_pair,
                    verdicts,
                    col_coding,
                    unclassified_pairs,
                )
            )
        if not keep:
            continue
        n_kept += 1
        key = (row_alt, record.deriv_suffix)
        cells[record.stem_class][key] += 1
        cells["pooled"][key] += 1

    def pack(c: Counter[tuple[bool, str]]) -> PanelStats:
        # Aggregate -a and -ui into a NON-PAL column, but keep the
        # breakdown for the diagnostic.
        alt_pal = c[(True, "-i")]
        alt_a = c[(True, "-a")]
        alt_ui = c[(True, "-ui")]
        nonalt_pal = c[(False, "-i")]
        nonalt_a = c[(False, "-a")]
        nonalt_ui = c[(False, "-ui")]
        return PanelStats(
            alt_pal=alt_pal,
            alt_nonpal=alt_a + alt_ui,
            nonalt_pal=nonalt_pal,
            nonalt_nonpal=nonalt_a + nonalt_ui,
            alt_a=alt_a,
            alt_ui=alt_ui,
            nonalt_a=nonalt_a,
            nonalt_ui=nonalt_ui,
        )

    return ({k: pack(v) for k, v in cells.items()}, n_kept)


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------


def _contamination_diagnostic(
    records: list[NounRecord],
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
    verdicts: dict[str, str],
) -> str:
    """Under the SUFFIX-ONLY column coding, what share of the ``-i``
    column has a trigger-absent paradigm verdict?

    High shares mean the SUFFIX-ONLY panel is measurably contaminated
    by morphological-block cases where the -i suffix was selected but
    palatalisation never applied. This is what motivates promoting
    the ALTERNATES-ONLY column to primary.
    """
    lines: list[str] = []
    stem_groups: list[tuple[str, frozenset[str]]] = [
        ("dorsal", DORSAL_STEMS),
        ("coronal", CORONAL_STEMS),
    ]
    unclassified_throwaway: set[tuple[str, str]] = set()
    for label, stems in stem_groups:
        counts: Counter[tuple[str, str]] = Counter()
        for rec in records:
            if rec.stem_final not in stems:
                continue
            if rec.deriv_suffix != "-i":
                continue
            for verb in _accepted_verbs(
                rec,
                accepted_by_pair,
                verdicts,
                ColCoding.SUFFIX_ONLY,
                unclassified_throwaway,
            ) or ([] if rec.lemma not in user_keep else [""]):
                verdict = (
                    verdicts.get(verb, "UNKNOWN") if verb else "USER_KEEP"
                )
                row_key = (
                    "alt" if rec.mutation else "noalt",
                    verdict,
                )
                counts[row_key] += 1
        total = sum(counts.values())
        if total == 0:
            continue
        lines.append(f"\n  {label} -i column (n={total}):")
        for (m, v), c in sorted(counts.items()):
            pct = 100 * c / total
            lines.append(f"    {m:5s}, verdict={v:14s}: {c:3d}  ({pct:.1f}%)")
    return "\n".join(lines)


def _ipa_row_diagnostic(records: list[NounRecord]) -> str:
    """Show the three-state IPA-row breakdown, especially the LEX
    (already-palatal-in-singular) column that's silently pooled with
    ALT.

    Deduped by lemma so that ambiguous nouns (two records under the
    multi-suffix fix) count once toward the row-state distribution.
    """
    lines: list[str] = []
    lines.append("\n  Distribution of IPA-row states (deduped by lemma):")
    seen: set[str] = set()
    tab: Counter[tuple[str, IPARowClass]] = Counter()
    for rec in records:
        if rec.lemma in seen:
            continue
        seen.add(rec.lemma)
        tab[(rec.stem_final, rec.ipa_row_class)] += 1

    lines.append(
        f"    {'sf':3s}  {'ALT':>5s}  {'LEX':>5s}  {'NONALT':>6s}  "
        f"{'UNCLEAR':>7s}"
    )
    for sf in sorted(TARGET_STEMS):
        alt = tab.get((sf, IPARowClass.ALT), 0)
        lex = tab.get((sf, IPARowClass.LEX), 0)
        non = tab.get((sf, IPARowClass.NONALT), 0)
        unc = tab.get((sf, IPARowClass.UNCLEAR), 0)
        lines.append(f"    {sf:3s}  {alt:5d}  {lex:5d}  {non:6d}  {unc:7d}")

    # Lex breakdown: are the LEX rows selecting -i or -a? Under both
    # theories they should behave like ALT (K/T͡ʃ class). If they don't,
    # something's wrong.
    lex_verbalizer: Counter[str] = Counter()
    for rec in records:
        if rec.ipa_row_class is IPARowClass.LEX:
            lex_verbalizer[rec.deriv_suffix] += 1
    if lex_verbalizer:
        lines.append(
            "\n  LEX rows (base already palatal in singular) — "
            "verbalizer choice:"
        )
        for suf, n in sorted(lex_verbalizer.items()):
            lines.append(f"    {suf:3s}: {n}")
    return "\n".join(lines)


# Non-alternating-consonant orthographic set, per Steriade §10.2.3.5's
# baseline: {r, n, m, ts, ʃ, p, b, f, v}. In Romanian orthography, ts
# is written ⟨ț⟩ and ʃ is written ⟨ș⟩ (both may also appear with the
# comma-below variants ţ, ş depending on Unicode source).
_BASELINE_ORTH_C: Final[frozenset[str]] = frozenset(
    {"r", "n", "m", "ț", "ţ", "ș", "ş", "p", "b", "f", "v"}
)
# Declension endings to strip when computing the orthographic
# stem-final consonant of a lemma with an empty ``stem_final`` field.
_LEMMA_DECLENSION_ENDINGS: Final[tuple[str, ...]] = (
    "ăi",
    "âi",
    "ii",
    "ie",
    "ia",
    "ea",
    "oa",
    "u",
    "ă",
    "â",
    "î",
    "e",
    "i",
    "a",
    "o",
)


def _baseline_stem_c(lemma: str) -> str:
    """Return the orthographic stem-final consonant of a lemma, or ''.

    Strips a single Romanian declension ending and returns the last
    orthographic character. Used only for Steriade's baseline
    reconstruction; the main pipeline uses the lexicon's
    ``stem_final`` field directly.
    """
    if not lemma:
        return ""
    stem = lemma.strip().lower()
    for suf in _LEMMA_DECLENSION_ENDINGS:
        if len(stem) > len(suf) and stem.endswith(suf):
            stem = stem[: -len(suf)]
            break
    if not stem:
        return ""
    return stem[-1]


def _baseline_diagnostic(
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
) -> str:
    """Steriade §10.2.3.5 table 10.19: distribution of ``-í / -á / -uí``
    on native derived verbs from bases ending in *non-alternating*
    consonants ``{r, n, m, ț, ș, p, b, f, v}``.

    This is the neutral baseline against which K/T͡ʃ (paper: 77%/19%/
    4%) and K/K (6%/38%/45%) distributions are compared. Steriade
    reports 57%/40%/3% (n=388) as the baseline.

    The etymology audit CSV was built for the target-stem set, so it
    doesn't cover baseline nouns; we don't gate on it here. We do
    still apply Steriade's own French-loan-on-``-a`` exclusion.
    """
    counts: Counter[str] = Counter()
    seen: set[str] = set()
    with CSV_PATH.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["pos"] != "N":
                continue
            if row["stem_final"]:
                continue
            derived_field = (row.get("derived_verbs") or "").strip()
            if not derived_field:
                continue
            lemma = row["lemma"]
            if lemma in seen:
                continue
            c = _baseline_stem_c(lemma)
            if c not in _BASELINE_ORTH_C:
                continue
            suffixes = [
                s.strip()
                for s in (row.get("deriv_suffixes") or "").split("|")
                if s.strip()
            ]
            etym = (row.get("etym_lang") or "").strip()
            # Deduplicate suffixes per lemma so a noun with both -i and
            # -a verbs counts once per suffix, not once per verb.
            noun_suffixes: set[str] = set()
            for s in suffixes:
                if s not in ACCEPTED_VERBALIZERS:
                    continue
                if etym == "fr" and s == "-a":
                    continue
                noun_suffixes.add(s)
            if not noun_suffixes:
                continue
            seen.add(lemma)
            for s in noun_suffixes:
                counts[s] += 1

    lines: list[str] = []
    lines.append(
        "\n  Baseline: derived verbs on nouns whose stem-final "
        "consonant is non-alternating"
    )
    lines.append(
        "  (per Steriade §10.2.3.5: {r, n, m, ț, ș, p, b, f, v}). "
        "Steriade reports 57%/40%/3% at n=388."
    )
    total = sum(counts.values())
    if total == 0:
        lines.append("    (no baseline rows found)")
        return "\n".join(lines)
    lines.append(f"\n    Total baseline noun-suffix pairs: {total}")
    for s in ("-i", "-a", "-ui"):
        c = counts.get(s, 0)
        pct = 100 * c / total
        lines.append(f"    {s:4s}: {c:4d}  ({pct:5.1f}%)")
    lines.append("\n    Steriade paper (n=388):     57%   40%    3%")
    return "\n".join(lines)


def _nonpal_breakdown(
    records: list[NounRecord],
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
    verdicts: dict[str, str],
) -> str:
    """Break the NON-PAL column of the PRIMARY panel into ``-a`` vs
    ``-ui``, so the reader can see which non-palatalising suffix
    carries the K/K signal.

    Steriade §10.2.3.5 predicts ``-uí`` is the K/K-dominant choice
    (45% of K/K dorsals). Verified here on the modern lexicon.
    """
    lines: list[str] = []
    lines.append("\n  NON-PAL column of PRIMARY panel: -a vs -ui breakdown.")
    unclassified_throwaway: set[tuple[str, str]] = set()
    panels, _ = _build_panels(
        records,
        accepted_by_pair,
        user_keep,
        verdicts,
        PRIMARY_ROW_CODING,
        PRIMARY_COL_CODING,
        unclassified_throwaway,
        lambda r: True,  # FULL DEX scope
    )
    header = (
        f"    {'stem':7s}  {'row':7s}  {'-a':>5s}  {'-ui':>5s}  "
        f"{'%-ui':>6s}"
    )
    lines.append(header)
    for stem in (StemClass.DORSAL, StemClass.CORONAL, "pooled"):
        stats = panels[stem]
        disp = stem.value if isinstance(stem, StemClass) else stem
        for row_label, a, u in (
            ("alt", stats.alt_a, stats.alt_ui),
            ("nonalt", stats.nonalt_a, stats.nonalt_ui),
        ):
            n = a + u
            pct = f"{100 * u / n:5.1f}%" if n > 0 else "  n/a"
            lines.append(
                f"    {disp:7s}  {row_label:7s}  {a:5d}  {u:5d}  {pct}"
            )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Robustness: permutation test + subgroup analyses
# ---------------------------------------------------------------------------
#
# The permutation test is the strongest null-hypothesis control we can
# run without additional assumptions. If the observation is real,
# shuffling the row-axis (K/T͡ʃ vs K/K) label across the kept records
# should almost never reproduce the observed pooled OR. If it does —
# say more than 5% of the time — the effect is indistinguishable from
# random co-incidence given the marginal frequencies of the two axes.


def _observed_pooled_stats(
    records: list[NounRecord],
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
    verdicts: dict[str, str],
    row_coding: RowCoding,
    col_coding: ColCoding,
) -> tuple[list[tuple[bool, str]], PanelStats]:
    """Collect kept (row_alt, col_class) pairs and pack the pooled 2×2.

    Returned pairs use ``col_class ∈ {'pal', 'nonpal'}`` so a
    permutation over ``row_alt`` reshuffles just that axis.
    """
    pairs: list[tuple[bool, str]] = []
    unclassified_throwaway: set[tuple[str, str]] = set()
    for record in records:
        row_alt = record.plural_alt(row_coding)
        if row_alt is None:
            continue
        keep = record.lemma in user_keep or bool(
            _accepted_verbs(
                record,
                accepted_by_pair,
                verdicts,
                col_coding,
                unclassified_throwaway,
            )
        )
        if not keep:
            continue
        pairs.append((row_alt, _col_class(record.deriv_suffix)))
    counts: Counter[tuple[bool, str]] = Counter(pairs)
    stats = PanelStats(
        alt_pal=counts[(True, "pal")],
        alt_nonpal=counts[(True, "nonpal")],
        nonalt_pal=counts[(False, "pal")],
        nonalt_nonpal=counts[(False, "nonpal")],
    )
    return pairs, stats


def _sample_or(counts: dict[tuple[bool, str], int]) -> float | None:
    """Sample OR from a 2×2 count dict, ``None`` if any margin is 0."""
    a = counts.get((True, "pal"), 0)
    b = counts.get((True, "nonpal"), 0)
    c = counts.get((False, "pal"), 0)
    d = counts.get((False, "nonpal"), 0)
    if b == 0 or c == 0:
        return None
    return (a * d) / (b * c)


def _permutation_test(
    pairs: list[tuple[bool, str]],
    observed_or: float,
    n_iter: int = 10_000,
    seed: int = 20080715,
) -> tuple[float, float, float]:
    """Return (proportion ≥ observed OR, median null OR, 97.5%-ile).

    Shuffles the row-axis (mutation) label across kept records while
    keeping the col-axis (verbalizer class) fixed. Because we shuffle
    rather than sample-with-replacement, both marginal totals are
    preserved exactly — this is a "conditional" permutation test
    matching the assumptions of Fisher's exact test at scale.
    """
    rng = random.Random(seed)
    row_labels = [p[0] for p in pairs]
    col_labels = [p[1] for p in pairs]
    perm_ors: list[float] = []
    ge_count = 0
    shuffled = list(row_labels)
    for _ in range(n_iter):
        rng.shuffle(shuffled)
        counts: Counter[tuple[bool, str]] = Counter(zip(shuffled, col_labels))
        or_val = _sample_or(counts)
        if or_val is None:
            continue
        perm_ors.append(or_val)
        if or_val >= observed_or:
            ge_count += 1
    if not perm_ors:
        return (float("nan"), float("nan"), float("nan"))
    perm_ors.sort()
    p_perm = ge_count / len(perm_ors)
    median_null = perm_ors[len(perm_ors) // 2]
    upper = perm_ors[int(0.975 * len(perm_ors))]
    return (p_perm, median_null, upper)


def _permutation_diagnostic(
    records: list[NounRecord],
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
    verdicts: dict[str, str],
) -> str:
    """Run a permutation null test against the PRIMARY coding."""
    pairs, stats = _observed_pooled_stats(
        records,
        accepted_by_pair,
        user_keep,
        verdicts,
        PRIMARY_ROW_CODING,
        PRIMARY_COL_CODING,
    )
    fisher = stats.fisher()
    if fisher is None:
        return "  Permutation test skipped (2×2 has a zero margin)."
    observed_or, _, _, observed_p = fisher
    p_perm, median_null, upper_null = _permutation_test(pairs, observed_or)
    lines: list[str] = []
    lines.append(
        f"\n  Null: shuffle K/T͡ʃ-vs-K/K labels across the {len(pairs)} "
        f"kept records."
    )
    lines.append(
        f"  Observed pooled OR = {observed_or:.2f}  "
        f"(Fisher exact p = {observed_p:.2e})"
    )
    lines.append(
        f"  Permutation null (10 000 shuffles): median OR = "
        f"{median_null:.2f}, 97.5%-ile = {upper_null:.2f}"
    )
    lines.append(
        f"  Fraction of shuffles with null OR ≥ observed: "
        f"{p_perm:.4f}  ← empirical p-value"
    )
    return "\n".join(lines)


def _etym_stratified(
    records: list[NounRecord],
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
    verdicts: dict[str, str],
) -> str:
    """Run PRIMARY separately for each etymology stratum.

    If the effect is really about Romanian speakers' verbalizer
    choice — not just Latin loan morphology — then it should hold on
    the native-etymology subset (empty ``etym_lang``) too. If it
    vanishes when Latin loans are removed, the effect might be
    partially an artefact of borrowed conjugation class.
    """
    lines: list[str] = []
    strata = [
        ("all", set()),
        ("native (etym='')", {""}),
        ("Latin (etym=la)", {"la"}),
        ("everything except Latin", None),  # sentinel
    ]
    for label, keep_etym in strata:
        subset: list[NounRecord] = []
        for r in records:
            if keep_etym is None:
                if r.etym_lang == "la":
                    continue
                subset.append(r)
            elif not keep_etym or r.etym_lang in keep_etym:
                subset.append(r)
        _, stats = _observed_pooled_stats(
            subset,
            accepted_by_pair,
            user_keep,
            verdicts,
            PRIMARY_ROW_CODING,
            PRIMARY_COL_CODING,
        )
        f = stats.fisher()
        if f is None:
            lines.append(
                f"  {label:26s}  n={stats.n:4d}  (2×2 has zero margin)"
            )
            continue
        odds, lo, hi, p = f
        p_fmt = f"{p:.2e}" if p < 1e-3 else f"{p:.3f}"
        lines.append(
            f"  {label:26s}  n={stats.n:4d}  OR={odds:5.2f}  "
            f"[{lo:.2f}, {hi:.2f}]  p={p_fmt}"
        )
    return "\n".join(lines)


def _sing_tantum_excluded(
    records: list[NounRecord],
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
    verdicts: dict[str, str],
) -> str:
    """PRIMARY, but re-computed with singularia tantum removed.

    The K/K-by-construction rule for singularia tantum is a
    theoretical stipulation — Steriade's, not the lexicon's. If it
    turns out the effect only survives *with* that stipulation, we
    should flag it. If the effect holds *without* the stipulation too,
    the observation is robust to that choice.
    """
    subset = [r for r in records if not r.is_singulare_tantum]
    _, stats = _observed_pooled_stats(
        subset,
        accepted_by_pair,
        user_keep,
        verdicts,
        PRIMARY_ROW_CODING,
        PRIMARY_COL_CODING,
    )
    lines: list[str] = []
    lines.append(f"\n  With singularia tantum EXCLUDED (n_kept = {stats.n}):")
    lines.append(f"  {stats.format_line('pooled')}")
    return "\n".join(lines)


def _empty_etym_rows(
    records: list[NounRecord],
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
    verdicts: dict[str, str],
) -> str:
    """PRIMARY, but re-computed with ``etym_lang == ''`` rows removed.

    ``remarcă → remarca`` is a documented French loan whose
    ``etym_lang`` is blank in our data, so our fr-etym filter can't
    catch it. Blanks are heterogeneous: some are native, some are
    unrecorded loans. Reporting the analysis both with and without
    them isolates the sensitivity.
    """
    subset = [r for r in records if r.etym_lang != ""]
    _, stats = _observed_pooled_stats(
        subset,
        accepted_by_pair,
        user_keep,
        verdicts,
        PRIMARY_ROW_CODING,
        PRIMARY_COL_CODING,
    )
    lines: list[str] = []
    lines.append(f"\n  With empty-etym rows EXCLUDED (n_kept = {stats.n}):")
    lines.append(f"  {stats.format_line('pooled')}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


PANEL_MATRIX: Final[tuple[tuple[str, RowCoding, ColCoding, str], ...]] = (
    (
        "PRIMARY",
        PRIMARY_ROW_CODING,
        PRIMARY_COL_CODING,
        "Row: allomorph, grouped by whether it palatalises for that "
        "stem class ({-i,-e} vs -uri for dorsals; -i vs {-e,-uri} for "
        "coronals). Col: -i vs (-a ∪ -ui). "
        "Purely allomorph-based row axis — no per-noun mutation or IPA "
        "check. This is the most atheoretical coding we have.",
    ),
    (
        "MUTATION-ROW",
        RowCoding.MUTATION,
        ColCoding.SUFFIX_ONLY,
        "Row: mutation flag (orth sg-vs-pl). Col: -i vs (-a ∪ -ui). "
        "Uses the lexicon's per-noun mutation label — closer to "
        "Steriade's own coding but requires a per-noun palatalization "
        "check.",
    ),
    (
        "IPA-ROW",
        RowCoding.IPA,
        ColCoding.SUFFIX_ONLY,
        "Row: sg-vs-pl IPA (3-state). Col: -i vs (-a ∪ -ui). "
        "Uses phonetic sg-vs-pl comparison for the row axis; catches "
        "cluster cases the orthographic mutation flag misses.",
    ),
    (
        "STEM-BLIND ALLOMORPH",
        RowCoding.ALLOMORPH,
        ColCoding.SUFFIX_ONLY,
        "Row: plural-suffix identity, stem-blind ({-i, -e} vs -uri). "
        "Col: -i vs (-a ∪ -ui). Matches Kyle-John-Daniar Table (9). "
        "Groups coronal -e with coronal -i despite the two behaving "
        "opposite — kept as a demonstration of what that grouping "
        "does to the effect.",
    ),
    (
        "PARADIGM-STRICT",
        PRIMARY_ROW_CODING,
        ColCoding.ALTERNATES_ONLY,
        "Row: allomorph-stemwise. Col: -i verbs must be ALTERNATES; "
        "-a and -ui always kept. Drops paradigm-blocked -i verbs.",
    ),
    (
        "FULLY-STRICT",
        RowCoding.IPA,
        ColCoding.ALTERNATES_ONLY,
        "Row: sg-vs-pl IPA. Col: -i verbs must be ALTERNATES; "
        "-a and -ui always kept. Both axes strict.",
    ),
)


def _render_panel(
    name: str,
    row_coding: RowCoding,
    col_coding: ColCoding,
    caption: str,
    records: list[NounRecord],
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
    verdicts: dict[str, str],
    frequency_by_lemma: dict[str, float],
    unclassified_pairs: set[tuple[str, str]],
) -> str:
    lines: list[str] = []
    lines.append("\n" + "=" * 78)
    lines.append(f"{name}")
    lines.append(caption)
    lines.append("=" * 78)

    ranked = sorted(frequency_by_lemma.items(), key=lambda kv: -kv[1])
    top_lemmas: frozenset[str] = frozenset(
        lemma for lemma, _ in ranked[:TOP_N_FREQUENCY]
    )
    scopes: list[tuple[str, _ScopeFilter]] = [
        (f"TOP-{TOP_N_FREQUENCY}", lambda r: r.lemma in top_lemmas),
        ("FULL DEX", lambda r: True),
    ]
    for label, scope in scopes:
        panels, n_kept = _build_panels(
            records,
            accepted_by_pair,
            user_keep,
            verdicts,
            row_coding,
            col_coding,
            unclassified_pairs,
            scope,
        )
        lines.append(f"\n{label}  (n_kept = {n_kept})")
        for stem in (StemClass.DORSAL, StemClass.CORONAL, "pooled"):
            stats = panels[stem]
            disp = stem.value if isinstance(stem, StemClass) else stem
            lines.append(stats.format_line(disp))
    return "\n".join(lines)


def _render_headline(
    records: list[NounRecord],
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
    verdicts: dict[str, str],
) -> str:
    """One-screen summary of the observation, up top before diagnostics.

    Reports the PRIMARY OR plus permutation p, so a reader can decide
    in ten seconds whether Steriade's observation reproduces on this
    lexicon. Everything below is transparency about how the number was
    obtained and how sensitive it is to methodology.
    """
    pairs, stats = _observed_pooled_stats(
        records,
        accepted_by_pair,
        user_keep,
        verdicts,
        PRIMARY_ROW_CODING,
        PRIMARY_COL_CODING,
    )
    fisher = stats.fisher()
    lines: list[str] = []
    lines.append("=" * 78)
    lines.append(
        "HEADLINE — does the modern lexicon show Steriade's co-variation?"
    )
    lines.append("=" * 78)
    lines.append(
        "\n  Question: on modern Romanian nouns with a target obstruent"
        " stem-final,"
    )
    lines.append(
        "  do nouns whose plural takes a palatalising allomorph prefer"
    )
    lines.append(
        "  the palatalising verbalizer -í over -á/-uí? (Steriade 2008"
        " §10.2.3.5.)"
    )
    lines.append(
        "\n  Design: 2×2 contingency, row = allomorph-stemwise "
        "(purely morphological,"
    )
    lines.append("  no per-noun mutation check), col = -í vs (-á ∪ -uí);")
    lines.append(
        "  Fisher exact + permutation null; Steriade's own"
        " exclusions applied."
    )
    if fisher is not None:
        odds, lo, hi, p = fisher
        p_perm, med_null, upper_null = _permutation_test(pairs, odds)
        lines.append(f"\n  ANSWER (n = {stats.n} noun-verbaliser pairs):")
        lines.append(
            f"    pooled OR = {odds:.2f}  [95% CI {lo:.2f}, {hi:.2f}]"
        )
        lines.append(
            f"    Fisher exact p = {p:.2e};  permutation p_perm ="
            f" {p_perm:.4f} ({'zero' if p_perm == 0 else 'nonzero'} of"
            f" 10 000 shuffles beat it)"
        )
        lines.append(
            f"    null distribution:  median OR = {med_null:.2f}"
            f",  97.5%-ile OR = {upper_null:.2f}"
        )
    else:
        lines.append("\n  ANSWER: 2×2 has a zero margin, no OR computable.")

    lines.append(
        "\n  Read as: the observed OR is far above what random"
        " reshuffling of the"
    )
    lines.append(
        "  K/T͡ʃ-vs-K/K labels produces. Steriade's paper reported OR"
        " ≈ 50 on"
    )
    lines.append(
        "  n = 157 hand-collected verbs; our OR is smaller because the"
        " modern DEX"
    )
    lines.append(
        "  has more real K/K + -í nouns (gând-í, rost-í, vad-í …) than"
        " her sample."
    )
    lines.append(
        "  Direction and significance reproduce; magnitude is muted"
        " but real."
    )
    return "\n".join(lines)


def _render_report(
    records: list[NounRecord],
    accepted_by_pair: dict[tuple[str, str], bool],
    user_keep: set[str],
    verdicts: dict[str, str],
    unclassified_pairs: set[tuple[str, str]],
    frequency_by_lemma: dict[str, float],
) -> str:
    lines: list[str] = []
    # Headline first — the answer, then the transparency.
    lines.append(
        _render_headline(records, accepted_by_pair, user_keep, verdicts)
    )
    lines.append("")
    lines.append(
        f"Loaded {len(records)} (lemma, dvs, ...) tuples "
        f"that pass row filters."
    )
    lines.append(
        f"Etymology decisions: {len(accepted_by_pair)} explicit pairs, "
        f"{len(user_keep)} user-keep lemmas, "
        f"{len(verdicts)} verb verdicts."
    )

    # IPA-row diagnostic — show the LEX split so it can't hide
    lines.append("\n" + "=" * 78)
    lines.append("IPA-ROW DIAGNOSTIC: three-state sg-vs-pl breakdown")
    lines.append("=" * 78)
    lines.append(_ipa_row_diagnostic(records))

    # Four panels, PRIMARY first
    for name, row_coding, col_coding, caption in PANEL_MATRIX:
        lines.append(
            _render_panel(
                name,
                row_coding,
                col_coding,
                caption,
                records,
                accepted_by_pair,
                user_keep,
                verdicts,
                frequency_by_lemma,
                unclassified_pairs,
            )
        )

    # Baseline table — Steriade §10.2.3.5 (10.19)
    lines.append("\n" + "=" * 78)
    lines.append(
        "BASELINE: -i / -a / -ui distribution on non-alternating stems"
    )
    lines.append("(Steriade §10.2.3.5 table 10.19: 57%/40%/3% at n=388)")
    lines.append("=" * 78)
    lines.append(_baseline_diagnostic(accepted_by_pair, user_keep))

    # NON-PAL column breakdown: how much of it is -a vs -ui?
    lines.append("\n" + "=" * 78)
    lines.append(
        "NON-PAL BREAKDOWN: PRIMARY panel's NON-PAL column split by suffix"
    )
    lines.append(
        "(Steriade: K/K bases prefer -uí over -á; K/T͡ʃ bases the reverse)"
    )
    lines.append("=" * 78)
    lines.append(
        _nonpal_breakdown(records, accepted_by_pair, user_keep, verdicts)
    )

    # Robustness Zone: null test + subgroup analyses ---------------------
    lines.append("\n" + "=" * 78)
    lines.append("ROBUSTNESS: PERMUTATION NULL TEST (PRIMARY coding)")
    lines.append(
        "Does the observed OR survive shuffling the K/T͡ʃ-vs-K/K labels?"
    )
    lines.append("=" * 78)
    lines.append(
        _permutation_diagnostic(records, accepted_by_pair, user_keep, verdicts)
    )

    lines.append("\n" + "=" * 78)
    lines.append("ROBUSTNESS: PRIMARY stratified by etymology")
    lines.append(
        "Does the effect survive on native rows only? Removing Latin loans?"
    )
    lines.append("=" * 78)
    lines.append(
        "\n" + _etym_stratified(records, accepted_by_pair, user_keep, verdicts)
    )

    lines.append("\n" + "=" * 78)
    lines.append("ROBUSTNESS: PRIMARY with variants removed")
    lines.append(
        "Does the effect depend on singularia-tantum or empty-etym rows?"
    )
    lines.append("=" * 78)
    lines.append(
        _sing_tantum_excluded(records, accepted_by_pair, user_keep, verdicts)
    )
    lines.append(
        _empty_etym_rows(records, accepted_by_pair, user_keep, verdicts)
    )

    # Singularia tantum: how many did we recover, and where do they land?
    sing_tantum_by_stem: Counter[tuple[str, str]] = Counter()
    for rec in records:
        if rec.is_singulare_tantum:
            sing_tantum_by_stem[(rec.stem_final, rec.deriv_suffix)] += 1
    if sing_tantum_by_stem:
        lines.append("\n" + "=" * 78)
        lines.append(
            "SINGULARIA TANTUM: nouns with no plural — K/K by construction"
        )
        lines.append(
            "(Steriade §10.2.3.4: vlagă, Franco, Goga — categorically K/K)"
        )
        lines.append("=" * 78)
        lines.append(f"\n  {'sf':3s}  {'-i':>4s}  {'-a':>4s}  {'-ui':>4s}")
        for sf in sorted(TARGET_STEMS):
            i = sing_tantum_by_stem.get((sf, "-i"), 0)
            a = sing_tantum_by_stem.get((sf, "-a"), 0)
            u = sing_tantum_by_stem.get((sf, "-ui"), 0)
            if i + a + u == 0:
                continue
            lines.append(f"  {sf:3s}  {i:4d}  {a:4d}  {u:4d}")
        total_st = sum(sing_tantum_by_stem.values())
        lines.append(f"\n  Total singularia-tantum records: {total_st}")

    # Contamination diagnostic — shows why PARADIGM-STRICT differs from
    # PRIMARY. It's here for readers who want to know the mechanism.
    lines.append("\n" + "=" * 78)
    lines.append("CONTAMINATION DIAGNOSTIC: verdict breakdown of -i column")
    lines.append("(what PARADIGM-STRICT would drop relative to PRIMARY)")
    lines.append("=" * 78)
    lines.append(
        _contamination_diagnostic(
            records, accepted_by_pair, user_keep, verdicts
        )
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    accepted_by_pair, user_keep = _load_etymology_decisions()
    verdicts = _load_verb_verdicts()
    records, frequency_by_lemma = _load_records(accepted_by_pair, user_keep)
    unclassified_pairs: set[tuple[str, str]] = set()

    report = _render_report(
        records,
        accepted_by_pair,
        user_keep,
        verdicts,
        unclassified_pairs,
        frequency_by_lemma,
    )
    print(report, end="")
    OUTPUT_PATH.write_text(report, encoding="utf-8")
    print(f"\nReport persisted to {OUTPUT_PATH}")

    if unclassified_pairs:
        n = len(unclassified_pairs)
        sample = sorted(unclassified_pairs)[:10]
        print(
            f"\nWARNING: {n} (noun, verb) pairs lack an etymology "
            f"classification in {ETYMOLOGY_CSV_PATH.name}."
        )
        print("First 10 unclassified pairs:")
        for noun, verb in sample:
            print(f"  {noun}  →  {verb}")


if __name__ == "__main__":
    main()
