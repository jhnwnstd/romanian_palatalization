"""Immutable segments and the unification semantics they live in.

A ``Segment`` is the framework's atomic phonological object: a feature
dictionary plus a display ``label`` (for tracing) plus a ``provenance``
tuple recording which rules have touched it. The label is *not* the
phonology — two segments with the same features are phonologically
identical regardless of label. Labels exist so trace output can keep
the underlying ``/K/`` distinct from the surface ``[k]`` even after
default fill-in has resolved them to the same feature bundle.

Unification semantics (Bale, Papillon & Reiss 2014) live here:

  - A feature with value ``"0"`` is *underspecified*: unification freely
    sets it to whatever the rule supplies.
  - A feature with an explicit value (``"+"`` or ``"-"``) is *fixed*:
    unification succeeds only if the rule supplies the same value, and
    fails (returns ``None``) if it supplies the opposite. This is the
    paper's inalterability mechanism — prespecified segments are immune
    to rules that would unify a contradictory value into them.

The diphthong-tie character ``o̯`` is tokenised as a single segment
``o̯`` (matching the inventory JSON key); the affricate tie bar
``͡`` ties the affricate pair into a single segment label like
``t͡ʃ``. The palatalization marker ``ʲ`` is a separate segment of its
own, attached after glide formation.

How :func:`tokenize_ipa` splits IPA into those labels is a per-analysis
convention, a :class:`TokenizerConfig`: which marks to skip, which
tokens to fold onto the inventory's symbols, whether combining
diacritics stay on their base, and what to do with the length mark
``ː``. The default, :data:`ROMANIAN_TOKENIZER`, is the Romanian
convention (drop ``ː``; fold ``ɛ``/``ɔ``/``g`` to ``e``/``o``/``ɡ``),
so every call without a config behaves as it always has. Another
analysis supplies its config as
:attr:`phonology.profile.AnalysisProfile.tokenizer` (design rule R8).
"""

from __future__ import annotations

import unicodedata
from dataclasses import dataclass, field, replace
from enum import StrEnum
from types import MappingProxyType
from typing import Final, Mapping

FeatureValue = str  # one of "+", "-", "0"
FeatureMap = Mapping[str, FeatureValue]

UNSPEC: Final[FeatureValue] = "0"

# IPA punctuation we tokenize as part of a segment, not as boundaries.
_TIE_BAR: Final[str] = "͡"  # combining double inverted breve
_INVERTED_BREVE_BELOW: Final[str] = "̯"  # diphthong tie under o, e, etc.
_PALATALIZATION: Final[str] = "ʲ"
# Tie bars a keep_combining tokenizer honours: the usual one above and
# the undertie U+035C (combining double breve below).
_TIE_BARS: Final[frozenset[str]] = frozenset({_TIE_BAR, "\u035c"})

# The length mark (U+02D0). lexpipe maps ASCII ":" onto it (R8), so it
# is the only length mark a tokenizer config has to handle.
LENGTH_MARK: Final[str] = "ː"

# Marks that never belong to a segment: whitespace, syllable dots,
# and stress marks (IPA primary/secondary and their ASCII stand-ins).
NON_SEGMENTAL_MARKS: Final[frozenset[str]] = frozenset(
    {" ", "\t", ".", "ˈ", "ˌ", "'", "`"}
)

# Unicode general categories of combining marks: non-spacing, spacing
# and enclosing. Every IPA diacritic (U+0300-U+036F) is in "Mn".
_COMBINING_CATEGORIES: Final[frozenset[str]] = frozenset({"Mn", "Mc", "Me"})


@dataclass(frozen=True, slots=True)
class Segment:
    """A feature bundle plus a display label and provenance trail.

    Two segments compare equal iff their labels, features and
    provenance all match — so updating any of those returns a new
    object. Use the helper methods (``unify``, ``delete``,
    ``with_provenance``) rather than constructing copies by hand;
    they preserve the invariant that ``features`` is a fresh mapping
    on every mutation.
    """

    label: str
    features: FeatureMap
    provenance: tuple[str, ...] = field(default_factory=tuple)

    def matches(self, pattern: FeatureMap) -> bool:
        """True iff every fixed value in ``pattern`` agrees with self.

        Pattern entries set to ``"0"`` are wildcards (no constraint).
        Pattern entries set to ``"+"`` or ``"-"`` must equal self's
        value for the same feature. A segment whose own value is
        ``"0"`` does *not* match an explicit pattern value: matching
        is asymmetric, used for rule-target preconditions, not for
        unification.
        """
        for feat, want in pattern.items():
            if want == UNSPEC:
                continue
            if self.features.get(feat, UNSPEC) != want:
                return False
        return True

    def unify(self, supply: FeatureMap) -> Segment | None:
        """Return self with ``supply`` unified in, or None on conflict.

        For each (feat, value) in supply:
          - If self's value is ``"0"`` we adopt ``value``.
          - If self's value equals ``value`` we no-op.
          - Otherwise unification fails: the rule cannot apply.

        ``"0"`` in supply means "supply nothing for this feature"; it
        is treated as a wildcard, not as a fixed value to write.
        """
        new_feats = dict(self.features)
        for feat, value in supply.items():
            if value == UNSPEC:
                continue
            current = new_feats.get(feat, UNSPEC)
            if current == UNSPEC:
                new_feats[feat] = value
            elif current != value:
                return None
        return replace(self, features=new_feats)

    def unify_verbose(
        self,
        supply: FeatureMap,
    ) -> tuple["Segment | None", tuple[str, str, str] | None]:
        """Like :meth:`unify` but also returns diagnostic info.

        On conflict returns ``(None, (feat, wanted, got))`` naming the
        first supply value that contradicted an explicit segment value.
        On success returns ``(new_segment, None)``. Used by
        :mod:`phonology.diagnostics.explain` to explain rule failures.
        """
        new_feats = dict(self.features)
        for feat, value in supply.items():
            if value == UNSPEC:
                continue
            current = new_feats.get(feat, UNSPEC)
            if current == UNSPEC:
                new_feats[feat] = value
            elif current != value:
                return None, (feat, value, current)
        return replace(self, features=new_feats), None

    def delete(self, names: frozenset[str]) -> Segment:
        """Return a copy with every named feature reset to ``"0"``.

        Used by the bleed rule, the only feature-changing rule in the
        Romanian coronal set (see latex.tex e:bleed-rev). All other
        rules in the system are feature-filling and use ``unify``.
        """
        if not names:
            return self
        new_feats = {
            k: (UNSPEC if k in names else v) for k, v in self.features.items()
        }
        for n in names:
            new_feats.setdefault(n, UNSPEC)
        return replace(self, features=new_feats)

    def with_label(self, label: str) -> Segment:
        return replace(self, label=label)

    def with_provenance(self, rule_name: str) -> Segment:
        return replace(self, provenance=self.provenance + (rule_name,))


Word = tuple[Segment, ...]


# ---------------------------------------------------------------------------
# IPA tokenizer
# ---------------------------------------------------------------------------


class LengthHandling(StrEnum):
    """What :func:`tokenize_ipa` does with the length mark ``ː``."""

    DROP = "DROP"  # removed: the inventory has no length contrast
    ATTACH = "ATTACH"  # ends the preceding segment's token: aː ≠ a
    SEPARATE = "SEPARATE"  # a token of its own


@dataclass(frozen=True, slots=True)
class TokenizerConfig:
    """One analysis's convention for splitting IPA into segment labels.

    Fields:

      - ``skip`` — single characters dropped wherever a segment could
        start: whitespace, syllable dots, stress marks
        (:data:`NON_SEGMENTAL_MARKS`). The length mark is not listed
        here; ``length`` decides its fate.
      - ``fold`` — finished token → the inventory's label for it. It
        is looked up on whole tokens, so a key may span several code
        points (``"t͜ʃ": "t͡ʃ"`` under ``keep_combining``); a token
        without an entry is kept as it is.
      - ``keep_combining`` — True: a token is a base character plus
        every combining mark after it, and a tie bar (U+0361, or the
        undertie U+035C) binds the next base, with its marks, into
        the same token, so ``r̝``, ``ɪ̯``, ``r̩`` and ``t̪͡s̪`` are one
        token each. False, the Romanian convention: a U+0361 tie bar
        binds exactly the next code point, a single U+032F glues to
        its base, and any other combining mark is a token of its own.
      - ``length`` — what happens to :data:`LENGTH_MARK`; see
        :class:`LengthHandling`. Under ``ATTACH`` a mark that follows
        no segment (string-initial, or after a skipped character) is
        a token of its own.

    Modifier letters other than the length mark, such as ``ʲ`` and
    ``ʰ``, are not combining marks: under every config they are tokens
    of their own. Input is expected in NFC (lexpipe normalises IPA,
    R8), the form inventory labels are written in.

    Immutable and hashable by value: ``skip`` is stored as a frozenset
    and ``fold`` as a read-only mapping. Construction rejects
    multi-character ``skip`` entries, a ``skip`` that lists the length
    mark, and empty ``fold`` keys or values.
    """

    skip: frozenset[str]
    fold: Mapping[str, str]
    keep_combining: bool
    length: LengthHandling

    def __hash__(self) -> int:
        # The generated hash would fail on the read-only ``fold`` view.
        return hash(
            (
                self.skip,
                frozenset(self.fold.items()),
                self.keep_combining,
                self.length,
            )
        )

    def __post_init__(self) -> None:
        skip = frozenset(self.skip)
        fold = MappingProxyType(dict(self.fold))
        if any(len(ch) != 1 for ch in skip):
            raise ValueError(
                f"skip entries must be single characters: {sorted(skip)!r}"
            )
        if LENGTH_MARK in skip:
            raise ValueError(
                f"{LENGTH_MARK!r} is governed by `length`; drop it with "
                "LengthHandling.DROP rather than listing it in `skip`"
            )
        if any(not k or not v for k, v in fold.items()):
            raise ValueError("fold maps non-empty tokens to non-empty labels")
        object.__setattr__(self, "skip", skip)
        object.__setattr__(self, "fold", fold)
        object.__setattr__(self, "length", LengthHandling(self.length))


ROMANIAN_TOKENIZER: Final[TokenizerConfig] = TokenizerConfig(
    skip=NON_SEGMENTAL_MARKS,
    # Common IPA characters that lax transcriptions use for phonemes
    # our inventory writes with a different symbol. Mapping to the
    # inventory's canonical form lets those transcriptions still
    # tokenise instead of falling through to "unknown segment".
    fold={"ɛ": "e", "ɔ": "o", "g": "ɡ"},
    keep_combining=False,
    # Some Wiktionary transcriptions write ː, but the Romanian
    # inventory has no length contrast.
    length=LengthHandling.DROP,
)
"""The Romanian tokenizer convention: the default of :func:`tokenize_ipa`
and of :attr:`phonology.profile.AnalysisProfile.tokenizer`."""


def _segment_end(ipa: str, start: int, *, keep_combining: bool) -> int:
    """Index just past the segment token that starts at ``ipa[start]``.

    Length marks are not included; :func:`tokenize_ipa` handles them.
    """
    n = len(ipa)
    if not keep_combining:
        # A tie bar joins this base with the next code point.
        if start + 2 < n and ipa[start + 1] == _TIE_BAR:
            return start + 3
        # The diphthong-tie U+032F glues to its base.
        if start + 1 < n and ipa[start + 1] == _INVERTED_BREVE_BELOW:
            return start + 2
        return start + 1
    end = start + 1
    while end < n:
        if ipa[end] in _TIE_BARS and end + 1 < n:
            end += 2  # the tie bar and the base it binds in
        elif unicodedata.category(ipa[end]) in _COMBINING_CATEGORIES:
            end += 1
        else:
            break
    return end


def tokenize_ipa(
    ipa: str,
    *,
    tokenizer: TokenizerConfig = ROMANIAN_TOKENIZER,
) -> tuple[str, ...]:
    """Split an IPA string into segment-label tokens.

    ``tokenizer`` is the analysis's convention (see
    :class:`TokenizerConfig`); the default :data:`ROMANIAN_TOKENIZER`
    gives exactly the tokens this function gave before it took one.

    Tie bars (U+0361) attach to the *preceding* code point and bind it
    with the following one, so ``t͡ʃ`` is a three-codepoint, one-segment
    sequence. The diphthong tie U+032F attaches to its base vowel and
    stays with it (``o̯``). The palatalization marker ``ʲ`` and the
    desyllabified glide carry no special handling — they are separate
    one-codepoint segments. Whitespace, periods, and stress marks
    (which ``normalize_ipa`` should already have stripped) are dropped
    defensively, as is the length mark under ``LengthHandling.DROP``.
    """
    out: list[str] = []
    i = 0
    n = len(ipa)
    while i < n:
        ch = ipa[i]
        if ch in tokenizer.skip:
            i += 1
            continue
        if ch == LENGTH_MARK:
            # Under ATTACH, a mark right after a segment was read with
            # it; one that reaches here follows no segment.
            if tokenizer.length is not LengthHandling.DROP:
                out.append(ch)
            i += 1
            continue
        end = _segment_end(ipa, i, keep_combining=tokenizer.keep_combining)
        if tokenizer.length is LengthHandling.ATTACH:
            while end < n and ipa[end] == LENGTH_MARK:
                end += 1
        token = ipa[i:end]
        out.append(tokenizer.fold.get(token, token))
        i = end
    return tuple(out)


def project(segment: Segment, inventory: Mapping[str, FeatureMap]) -> str:
    """Pick the inventory segment label that best fits ``segment``.

    Used to render an abstract feature bundle back to IPA after rules
    have run. The pick is: the inventory entry whose features are
    consistent with the segment (no contradicting explicit values) and
    that agrees on the largest number of explicit features. Ties are
    broken by inventory-key sort order, making the result
    deterministic.

    If the segment's own label already names an inventory entry whose
    features are consistent, we keep that label — this preserves
    underspecified placeholders like ``/K/`` if no rule has touched
    them yet, instead of silently snapping them to ``/k/``.
    """
    seg_feats = segment.features

    def _consistent(inv_feats: FeatureMap) -> bool:
        for feat, inv_val in inv_feats.items():
            if inv_val == UNSPEC:
                continue
            seg_val = seg_feats.get(feat, UNSPEC)
            if seg_val != UNSPEC and seg_val != inv_val:
                return False
        return True

    def _agreement(inv_feats: FeatureMap) -> int:
        score = 0
        for feat, inv_val in inv_feats.items():
            if inv_val == UNSPEC:
                continue
            if seg_feats.get(feat, UNSPEC) == inv_val:
                score += 1
        return score

    own = inventory.get(segment.label)
    if own is not None and _consistent(own):
        return segment.label

    best: tuple[int, str] | None = None
    for label in sorted(inventory.keys()):
        inv_feats = inventory[label]
        if not _consistent(inv_feats):
            continue
        score = _agreement(inv_feats)
        if best is None or score > best[0]:
            best = (score, label)
    return best[1] if best is not None else segment.label


def segments_to_ipa(word: Word, inventory: Mapping[str, FeatureMap]) -> str:
    """Render a Word back to an IPA string by projecting each segment."""
    return "".join(project(s, inventory) for s in word)
