"""RulePipeline: an ordered list of rules + the derivation it produces.

The pipeline walks each rule once in declared order. After each rule
it appends a :class:`Step` to the trace — even no-ops — so the printed
derivation table mirrors the paper's e:coronal-derivation layout
(latex.tex:472-492), one column per rule, even when the rule doesn't
fire on this stem.

The pipeline is generic in the inventory: it takes a resolver callable
that maps a segment label to a :class:`Segment` with features. This
keeps the pipeline free of any Romanian-specific knowledge and lets a
future analysis swap in a different inventory without changing the
core.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

from .rules import AllomorphSelection, ApplyResult, ResolverFn, Rule
from .segments import Word


@dataclass(frozen=True, slots=True)
class Step:
    """One row of the derivation table."""

    rule: str
    word: Word
    fired: bool


@dataclass(frozen=True, slots=True)
class Derivation:
    """Full derivation: UR, ordered steps, SR.

    ``sr`` is the post-last-step word — equivalent to
    ``steps[-1].word`` when at least one step exists, repeated here
    for convenient access.
    """

    ur: Word
    steps: tuple[Step, ...]
    sr: Word


@dataclass(frozen=True, slots=True)
class RulePipeline:
    """Apply an ordered list of rules and record the trace.

    ``rules`` is the declared order from the analysis module. The
    resolver maps an allomorph label like ``"u"`` to a fully-featured
    :class:`Segment`; it is consulted by :class:`AllomorphSelection`
    when expanding morpheme placeholders.
    """

    rules: tuple[Rule, ...]
    resolve: ResolverFn

    def derive(self, ur: Word) -> Derivation:
        steps: list[Step] = []
        current: Word = ur
        for rule in self.rules:
            result = self._apply_one(rule, current)
            steps.append(
                Step(
                    rule=rule.name,
                    word=result.word,
                    fired=bool(result.edits),
                )
            )
            current = result.word
        return Derivation(ur=ur, steps=tuple(steps), sr=current)

    def _apply_one(self, rule: Rule, word: Word) -> ApplyResult:
        if isinstance(rule, AllomorphSelection):
            return rule.apply_to(word, self.resolve)
        # UnificationRule, DeletionRule, GlideFormationRule,
        # SegmentDeletionRule, plus any future Rule-Protocol
        # implementers all take a bare `word` — dispatch through
        # duck-typing keeps the pipeline agnostic to the rule kind.
        return rule.apply(word)


__all__: Final[tuple[str, ...]] = ("Derivation", "RulePipeline", "Step")
