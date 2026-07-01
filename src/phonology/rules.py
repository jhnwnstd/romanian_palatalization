"""Rule kinds — the Ω of Logical Phonology, plus positional & morph.

Every rule is a *total function* on words: it returns an
:class:`ApplyResult` even when nothing matches. Composition of rules
in an ordered pipeline yields the phonology of a language.

The rule kinds correspond to LP's operator inventory (see
:mod:`phonology.lp` for the primitives themselves):

  - :class:`UnificationRule` — the ``⊔`` operator applied at every
    licensed target. Feature-filling: writes only where consistent,
    silently identity on conflicts. This is the LP inalterability
    mechanism: prespecified segments resist by failed unification.
  - :class:`DeletionRule` — the ``\\`` operator (feature subtraction)
    applied at every licensed target. Strips named features at
    whatever value they bear (LP: ``A \\ {+F, -F}``). Feature-
    changing: creates derived underspecification for the default
    fill-in to resolve.
  - :class:`SegmentDeletionRule` — the ``↦ ε`` operator on whole
    segments. Removes matched segments from the word. Included for
    LP completeness; the Romanian analysis doesn't use it.
  - :class:`GlideFormationRule` — a composed ``\\`` + ``⊔`` at a
    positional (word-final) target. Modelled as its own class
    because its licensing is positional, not feature-based.
  - :class:`AllomorphSelection` — morphological spell-out. Fills a
    morpheme placeholder with the phonological exponent chosen by
    a Pāṇinian list of (condition, segments) clauses.

**LP totality**. Within a single rule, all licensed positions are
computed simultaneously from the input word, so nothing feeds or
bleeds within one rule application (see LP notes §9). Rule
ordering between distinct rules IS extrinsic and does real work —
that's what the :class:`~phonology.pipeline.RulePipeline` handles.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Final, Protocol

from .search import Search
from .segments import FeatureMap, Segment, Word


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
    """LP ``⊔`` applied at every licensed target position.

    Every segment matching ``target`` runs its own SEARCH from the
    ORIGINAL input word; if the search returns a licensed trigger,
    ``supply`` unifies into that segment. All firings are computed
    from the input state and applied in one pass — this is LP
    totality: within one rule, nothing feeds or bleeds. Feeding /
    bleeding between distinct rules is what the pipeline's ordered
    composition handles.

    ``search=None`` is an environment-free default-fill rule: unify
    ``supply`` into every segment that matches ``target`` (the
    coronal default fill-in is three such rules).

    Prespecified segments are inalterable by unification failure:
    if a supply value contradicts an existing explicit value, unify
    returns None and the rule is identity on that position — LP's
    central "the operation, not the target, does the sorting".
    """

    name: str
    target: FeatureMap
    supply: FeatureMap
    search: Search | None = None

    def apply(self, word: Word) -> ApplyResult:
        return self._apply_simultaneous(word)

    def _apply_simultaneous(self, word: Word) -> ApplyResult:
        """All licensed targets fire in parallel from the input state.

        The "computed from the input" invariant matters when a rule
        has multiple candidates in a word: each target's search sees
        the ORIGINAL word, not a partially-updated one, so the rule
        is a well-defined total function on words rather than an
        order-dependent left-fold.
        """
        new = list(word)
        edits: list[Edit] = []
        for i, seg in enumerate(word):
            if not seg.matches(self.target):
                continue
            if self.search is not None:
                # SEARCH against the ORIGINAL word so early firings
                # don't affect later ones' triggers.
                if self.search.locate(word, i) is None:
                    continue
            unified = seg.unify(self.supply)
            if unified is None:
                continue
            if unified.features == seg.features:
                # Vacuous: every supply feature was already present.
                # Not a genuine firing — no edit, no provenance stamp.
                continue
            unified = unified.with_provenance(self.name)
            new[i] = unified
            edits.append(Edit(self.name, i, seg, unified))
        if not edits:
            return ApplyResult(word=word)
        return ApplyResult(word=tuple(new), edits=tuple(edits))


@dataclass(frozen=True, slots=True)
class DeletionRule:
    """LP ``\\`` applied at every licensed target position.

    Every segment matching ``target`` runs its own SEARCH from the
    ORIGINAL input word; each licensed target has the named features
    subtracted (stripped at whatever value they bear). All firings
    computed simultaneously from the input — LP totality.

    Feature-changing rules are ``\\`` then ``⊔`` at the pipeline
    level: this rule creates derived underspecification, a later
    default-fill rule supplies the unmarked values. LP predicts
    "no bypassing" and "no polarity rules" from this decomposition.
    """

    name: str
    target: FeatureMap
    clear: frozenset[str]
    search: Search

    def apply(self, word: Word) -> ApplyResult:
        new = list(word)
        edits: list[Edit] = []
        for i, seg in enumerate(word):
            if not seg.matches(self.target):
                continue
            if self.search.locate(word, i) is None:
                continue
            cleared = seg.delete(self.clear)
            if cleared.features == seg.features:
                continue
            cleared = cleared.with_provenance(self.name)
            new[i] = cleared
            edits.append(Edit(self.name, i, seg, cleared))
        if not edits:
            return ApplyResult(word=word)
        return ApplyResult(word=tuple(new), edits=tuple(edits))


@dataclass(frozen=True, slots=True)
class SegmentDeletionRule:
    """LP ``↦ ε`` on whole segments — removes matched segments from
    the word.

    LP predicts "Delete the Rich": if some but not all X delete in
    a context, the ones that delete are MORE specified than those
    that stay (because only a more-specified class can be carved
    out — see LP notes §7).

    Not used by the Romanian analysis, but included so the framework
    is LP-complete: a caller can implement any analysis translated
    into LP primitives without needing to extend the framework.
    """

    name: str
    target: FeatureMap
    search: Search | None = None

    def apply(self, word: Word) -> ApplyResult:
        # Mark each licensed target; then filter the word in one pass
        # so all deletions are simultaneous.
        to_delete: set[int] = set()
        edits: list[Edit] = []
        for i, seg in enumerate(word):
            if not seg.matches(self.target):
                continue
            if self.search is not None:
                if self.search.locate(word, i) is None:
                    continue
            to_delete.add(i)
        if not to_delete:
            return ApplyResult(word=word)
        new_word = tuple(
            seg for i, seg in enumerate(word) if i not in to_delete
        )
        for i in sorted(to_delete):
            # ``after`` is not meaningful for a deletion; reuse
            # ``before`` for the edit log so traces remain typed.
            edits.append(Edit(self.name, i, word[i], word[i]))
        return ApplyResult(word=new_word, edits=tuple(edits))


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
        unified = unified.with_label(self.glide_label).with_provenance(
            self.name
        )
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
    "SegmentDeletionRule",
    "UnificationRule",
    "is_placeholder",
    "placeholder",
)
