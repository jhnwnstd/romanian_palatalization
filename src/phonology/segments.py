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
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Final, Mapping

FeatureValue = str  # one of "+", "-", "0"
FeatureMap = Mapping[str, FeatureValue]

UNSPEC: Final[FeatureValue] = "0"

# IPA punctuation we tokenize as part of a segment, not as boundaries.
_TIE_BAR: Final[str] = "͡"  # combining double inverted breve
_INVERTED_BREVE_BELOW: Final[str] = "̯"  # diphthong tie under o, e, etc.
_PALATALIZATION: Final[str] = "ʲ"


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


def tokenize_ipa(ipa: str) -> tuple[str, ...]:
    """Split an IPA string into segment-label tokens.

    Tie bars (U+0361) attach to the *preceding* code point and bind it
    with the following one, so ``t͡ʃ`` is a three-codepoint, one-segment
    sequence. The diphthong tie U+032F attaches to its base vowel and
    stays with it (``o̯``). The palatalization marker ``ʲ`` and the
    desyllabified glide carry no special handling — they are separate
    one-codepoint segments. Whitespace, periods, and stress marks
    (which ``normalize_ipa`` should already have stripped) are dropped
    defensively.
    """
    # Marks we treat as segment-neutral: stress and syllable dots, and
    # the length mark ː that some Wiktionary transcriptions include but
    # our inventory does not distinguish.
    skip = {" ", "\t", ".", "ˈ", "ˌ", "'", "`", "ː"}
    # Common IPA characters that lax transcriptions use for phonemes
    # our inventory writes with a different symbol. Mapping to the
    # inventory's canonical form lets those transcriptions still
    # tokenise instead of falling through to "unknown segment".
    fold = {"ɛ": "e", "ɔ": "o", "g": "ɡ"}
    out: list[str] = []
    i = 0
    n = len(ipa)
    while i < n:
        ch = ipa[i]
        if ch in skip:
            i += 1
            continue
        if ch in fold:
            ch = fold[ch]
        # A tie bar joins this base with the next base into one token.
        if i + 1 < n and ipa[i + 1] == _TIE_BAR and i + 2 < n:
            out.append(ipa[i] + _TIE_BAR + ipa[i + 2])
            i += 3
            continue
        # The diphthong-tie U+032F glues to its base.
        if i + 1 < n and ipa[i + 1] == _INVERTED_BREVE_BELOW:
            out.append(ipa[i] + _INVERTED_BREVE_BELOW)
            i += 2
            continue
        out.append(ch)
        i += 1
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
