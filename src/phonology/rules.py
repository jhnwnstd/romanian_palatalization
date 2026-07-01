"""Rule kinds: unification, deletion, glide formation, allomorph spell-out.

Every rule has a ``name`` and an ``apply(word) -> ApplyResult`` method.
``ApplyResult`` carries the new word plus an edit log so traces can
record which segments changed at which index. Rules that don't fire
return the input unchanged with an empty edit log.

The split into four classes mirrors the paper:

  - :class:`UnificationRule` is the workhorse. Most rules in the paper
    are feature-filling — they unify a feature bundle into a target
    segment when the search finds a licensing trigger.
  - :class:`DeletionRule` is the one feature-changing rule (the bleed).
    Modelling it as its own class makes the "this derivation used the
    only feature-changing rule" question one ``isinstance`` away.
  - :class:`GlideFormationRule` is positional rather than feature-
    based: it fires on a literal /i/ at the right edge, after a
    consonant. Encoding it as its own class avoids contorting the
    Unification machinery with word-position knobs that no other rule
    needs.
  - :class:`AllomorphSelection` is the morphological spell-out. Given a
    morpheme placeholder at the right edge, it appends the chosen
    allomorph based on a Pāṇinian list of (condition, segments) pairs.

Search is always parameterised separately and passed in. Rules don't
know how to find their triggers; they only know what to do once a
trigger is found.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Final, Protocol

from .search import Search
from .segments import FeatureMap, Segment, UNSPEC, Word


@dataclass(frozen=True, slots=True)
class Edit:
    """One segment changed at a known index by a named rule."""

    rule: str
    index: int
    before: Segment
    after: Segment


@dataclass(frozen=True, slots=True)
class ApplyResult:
    """Output of ``Rule.apply``: the new word and what changed."""

    word: Word
    edits: tuple[Edit, ...] = ()


class Rule(Protocol):
    """A phonological rule.

    Implementations must be deterministic and side-effect-free. They
    take a ``Word`` and return an ``ApplyResult`` with the rewritten
    word plus an edit log. If the rule does not fire, return the
    input unchanged with no edits.
    """

    name: str

    def apply(self, word: Word) -> ApplyResult: ...


@dataclass(frozen=True, slots=True)
class UnificationRule:
    """Unify ``supply`` into the first target whose search succeeds.

    Iteration order is left-to-right over the word. For each segment
    matching ``target``, run the search; if the search returns an
    index (trigger licensed), unify ``supply`` into the segment. The
    rule applies at most once per word — multi-target firing is not
    a paper-side commitment we want to bake in. If your rule should
    cascade, run it twice in the pipeline.

    ``search=None`` is an environment-free default-fill rule: unify
    ``supply`` into *every* segment that matches ``target`` (the
    coronal default fill-in is three such rules).
    """

    name: str
    target: FeatureMap
    supply: FeatureMap
    search: Search | None = None

    def apply(self, word: Word) -> ApplyResult:
        if self.search is None:
            return self._apply_default(word)
        return self._apply_search(word)

    def _apply_search(self, word: Word) -> ApplyResult:
        for i, seg in enumerate(word):
            if not seg.matches(self.target):
                continue
            assert self.search is not None
            trigger = self.search.locate(word, i)
            if trigger is None:
                continue
            unified = seg.unify(self.supply)
            if unified is None:
                continue
            if unified.features == seg.features:
                # All supply features were already present; unification
                # is a no-op. Don't record a spurious edit — the rule
                # didn't change anything, so it didn't "fire".
                continue
            unified = unified.with_provenance(self.name)
            new = list(word)
            new[i] = unified
            return ApplyResult(
                word=tuple(new),
                edits=(Edit(self.name, i, seg, unified),),
            )
        return ApplyResult(word=word)

    def _apply_default(self, word: Word) -> ApplyResult:
        new = list(word)
        edits: list[Edit] = []
        for i, seg in enumerate(word):
            if not seg.matches(self.target):
                continue
            unified = seg.unify(self.supply)
            if unified is None or unified is seg:
                continue
            if unified.features == seg.features:
                continue
            unified = unified.with_provenance(self.name)
            new[i] = unified
            edits.append(Edit(self.name, i, seg, unified))
        return ApplyResult(word=tuple(new), edits=tuple(edits))


@dataclass(frozen=True, slots=True)
class DeletionRule:
    """Clear named features on every licensed target.

    Like :class:`UnificationRule`, iteration is left-to-right and the
    rule fires at most once per word. The clear set is applied
    wholesale: every named feature is reset to ``"0"`` whether it was
    explicit or already unspecified. Default-fill rules later in the
    pipeline supply the unmarked values.
    """

    name: str
    target: FeatureMap
    clear: frozenset[str]
    search: Search

    def apply(self, word: Word) -> ApplyResult:
        for i, seg in enumerate(word):
            if not seg.matches(self.target):
                continue
            if self.search.locate(word, i) is None:
                continue
            cleared = seg.delete(self.clear)
            if cleared.features == seg.features:
                continue
            cleared = cleared.with_provenance(self.name)
            new = list(word)
            new[i] = cleared
            return ApplyResult(
                word=tuple(new),
                edits=(Edit(self.name, i, seg, cleared),),
            )
        return ApplyResult(word=word)


@dataclass(frozen=True, slots=True)
class GlideFormationRule:
    """Desyllabify word-final /i/ after a consonant.

    Strips ``+Syllabic`` and supplies the multiply-articulated /j/'s
    feature bundle: ``-Syllabic``, postalveolar place
    (``+Coronal, -Anterior, +Distributed``), retaining ``+High,
    +Front, +Dorsal``. The label is rewritten to ``"j"`` so trace
    output and projection both name the desyllabified vocoid.

    Encoded as its own class because the rule's licensing condition is
    word-position (final), not a feature on a trigger segment, and
    forcing it through ``UnificationRule`` would require a positional
    Search variant no other rule needs.
    """

    name: str
    requires_preceding: FeatureMap
    glide_label: str = "j"
    glide_supply: FeatureMap = field(
        default_factory=lambda: {
            "Syllabic": "-",
            "CORONAL": "+",
            "Anterior": "-",
            "Distributed": "+",
        }
    )
    clear: frozenset[str] = frozenset({"Syllabic"})
    target_label: str = "i"

    def apply(self, word: Word) -> ApplyResult:
        if not word:
            return ApplyResult(word=word)
        last_idx = len(word) - 1
        last = word[last_idx]
        if last.label != self.target_label:
            return ApplyResult(word=word)
        if last_idx == 0:
            return ApplyResult(word=word)
        prev = word[last_idx - 1]
        if not prev.matches(self.requires_preceding):
            return ApplyResult(word=word)
        new_last = last.delete(self.clear)
        unified = new_last.unify(self.glide_supply)
        if unified is None:
            return ApplyResult(word=word)
        unified = unified.with_label(self.glide_label).with_provenance(self.name)
        new = list(word)
        new[last_idx] = unified
        return ApplyResult(
            word=tuple(new),
            edits=(Edit(self.name, last_idx, last, unified),),
        )


@dataclass(frozen=True, slots=True)
class AllomorphChoice:
    """One (condition, segments) clause of a spell-out rule.

    ``condition`` is a feature pattern that must match the last
    *non-morpheme* segment of the stem. ``segments`` is the IPA-side
    label sequence of the allomorph (e.g. ``("u", "r", "i")`` for
    ``-uri``). ``elsewhere`` clauses use ``condition={}`` (always
    matches) and must appear last in the option list.
    """

    label: str
    condition: FeatureMap
    segments: tuple[str, ...]


_MORPHEME_PREFIX: Final[str] = "@"


@dataclass(frozen=True, slots=True)
class AllomorphSelection:
    """Expand a morpheme placeholder into its phonological exponent.

    The pipeline inserts a placeholder segment with ``label =
    "@Num[Pl]"`` (or similar) at the right edge. This rule walks the
    placeholder's options in Pāṇinian order: the first clause whose
    condition matches the preceding stem segment wins. The placeholder
    is replaced by the chosen allomorph's segments.

    Allomorph segments are looked up in ``inventory_resolver`` — a
    callable supplied by the pipeline that maps a label like ``"u"``
    to a :class:`Segment` with the inventory's features.
    """

    name: str
    morpheme: str
    options: tuple[AllomorphChoice, ...]

    def apply_to(
        self,
        word: Word,
        resolve: "ResolverFn",
    ) -> ApplyResult:
        idx = self._placeholder_index(word)
        if idx is None:
            return ApplyResult(word=word)
        if idx == 0:
            return ApplyResult(word=word)
        stem_final = word[idx - 1]
        for choice in self.options:
            if stem_final.matches(choice.condition):
                exponent = tuple(
                    resolve(label).with_provenance(self.name)
                    for label in choice.segments
                )
                new = tuple(word[:idx]) + exponent
                edit = Edit(
                    self.name,
                    idx,
                    word[idx],
                    exponent[0] if exponent else word[idx],
                )
                return ApplyResult(word=new, edits=(edit,))
        return ApplyResult(word=word)

    def _placeholder_index(self, word: Word) -> int | None:
        target = _MORPHEME_PREFIX + self.morpheme
        for i, seg in enumerate(word):
            if seg.label == target:
                return i
        return None


# Type alias for the resolver callable. Keeping it here lets rules.py
# stay self-contained without importing FeatureInventory directly.
class ResolverFn(Protocol):
    def __call__(self, label: str) -> Segment: ...


def placeholder(morpheme: str) -> Segment:
    """Build a morpheme placeholder segment for AllomorphSelection.

    The placeholder carries no phonological features (everything is
    ``"0"``) so no rule can target it; only :class:`AllomorphSelection`
    finds it by label.
    """
    return Segment(label=_MORPHEME_PREFIX + morpheme, features={})


def is_placeholder(seg: Segment) -> bool:
    return seg.label.startswith(_MORPHEME_PREFIX)


__all__: Final[tuple[str, ...]] = (
    "AllomorphChoice",
    "AllomorphSelection",
    "ApplyResult",
    "DeletionRule",
    "Edit",
    "GlideFormationRule",
    "ResolverFn",
    "Rule",
    "UnificationRule",
    "is_placeholder",
    "placeholder",
)
