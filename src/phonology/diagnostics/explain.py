"""Per-rule firing/blocking explanations.

For each rule in a derivation, produces a structured record naming
which segments matched the target, whether the search found a
licensing trigger, and — when unification failed — which supply
feature contradicted which segment feature. The intent is to turn
a silent "rule no-op" line into an audit trail an analyst can read
without opening a Python debugger.

The framework's fast path (:meth:`RulePipeline.derive`) stays
unchanged; :func:`explain_derivation` walks the rule tuple again
with :meth:`Search.locate_verbose` and :meth:`Segment.unify_verbose`
to collect the reasons.

Rule-kind coverage:

- :class:`~phonology.rules.UnificationRule` — full diagnosis
  (target match, search outcome, unification conflicts).
- :class:`~phonology.rules.DeletionRule` — target match + search
  outcome (deletion itself cannot fail).
- :class:`~phonology.rules.GlideFormationRule` — checks the
  positional guards separately (word-final? preceding
  +Consonantal?).
- :class:`~phonology.rules.AllomorphSelection` — reports which
  allomorph clause matched (or the elsewhere fallback).
- Other rule kinds — a coarse "fired / didn't fire" summary from
  the before/after word diff.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Final

from ..pipeline import RulePipeline
from ..rules import (
    AllomorphSelection,
    DeletionRule,
    GlideFormationRule,
    Rule,
    UnificationRule,
    is_placeholder,
)
from ..search import SearchOutcome
from ..segments import Word


@dataclass(frozen=True, slots=True)
class TargetMatch:
    """A single (segment index, match/miss) datum for one rule."""

    index: int
    label: str
    matched_target: bool


@dataclass(frozen=True, slots=True)
class RuleTrace:
    """Structured explanation of one rule's step on a word."""

    rule_name: str
    fired: bool
    target_matches: tuple[TargetMatch, ...] = ()
    # For UnificationRule/DeletionRule that ran a search:
    search_outcome: SearchOutcome | None = None
    search_target_index: int | None = None
    search_trigger_index: int | None = None
    search_halted_at: int | None = None
    # For unification conflicts:
    unify_conflict: tuple[str, str, str] | None = None   # (feat, want, got)
    # For AllomorphSelection:
    allomorph_chosen: str | None = None
    allomorph_elsewhere: bool = False
    # Free-text summary for rule kinds that don't fit the fields above.
    notes: tuple[str, ...] = field(default_factory=tuple)


@dataclass(frozen=True, slots=True)
class Explanation:
    """A traced pipeline run: UR + per-rule RuleTrace + SR."""

    ur: Word
    rule_traces: tuple[RuleTrace, ...]
    sr: Word


def explain_derivation(
    pipeline: RulePipeline,
    ur: Word,
) -> Explanation:
    """Run the pipeline in narrating mode and return per-rule traces.

    The pipeline's fast path is unaffected; this function re-walks
    the same rules with verbose helpers to collect the diagnostics.
    """
    traces: list[RuleTrace] = []
    current: Word = ur
    for rule in pipeline.rules:
        trace, current = _explain_one(rule, current, pipeline)
        traces.append(trace)
    return Explanation(ur=ur, rule_traces=tuple(traces), sr=current)


def _explain_one(
    rule: Rule,
    word: Word,
    pipeline: RulePipeline,
) -> tuple[RuleTrace, Word]:
    if isinstance(rule, UnificationRule):
        return _explain_unification(rule, word)
    if isinstance(rule, DeletionRule):
        return _explain_deletion(rule, word)
    if isinstance(rule, GlideFormationRule):
        return _explain_glide(rule, word)
    if isinstance(rule, AllomorphSelection):
        return _explain_allomorph(rule, word, pipeline)
    # Fallback for unrecognised rule kinds — coarse diff.
    result = rule.apply(word)
    trace = RuleTrace(
        rule_name=getattr(rule, "name", type(rule).__name__),
        fired=bool(result.edits),
        notes=("no dedicated explainer for this rule kind",),
    )
    return trace, result.word


def _explain_unification(
    rule: UnificationRule,
    word: Word,
) -> tuple[RuleTrace, Word]:
    matches: list[TargetMatch] = []
    fired = False
    search_outcome: SearchOutcome | None = None
    search_target: int | None = None
    search_trigger: int | None = None
    search_halted: int | None = None
    conflict: tuple[str, str, str] | None = None

    new_word = list(word)

    if rule.search is None:
        # Default-fill: iterate whole word, unify into each match.
        for i, seg in enumerate(word):
            m = seg.matches(rule.target)
            matches.append(TargetMatch(i, seg.label, m))
            if not m:
                continue
            unified, conf = seg.unify_verbose(rule.supply)
            if unified is None:
                if conflict is None:
                    conflict = conf
                continue
            if unified.features == seg.features:
                continue
            new_word[i] = unified.with_provenance(rule.name)
            fired = True
        trace = RuleTrace(
            rule_name=rule.name,
            fired=fired,
            target_matches=tuple(matches),
            unify_conflict=conflict,
        )
        return trace, tuple(new_word)

    # Environment-driven rule: iterate left-to-right, first successful
    # (target + search + unify) fires.
    for i, seg in enumerate(word):
        m = seg.matches(rule.target)
        matches.append(TargetMatch(i, seg.label, m))
        if not m:
            continue
        report = rule.search.locate_verbose(word, i)
        if report.trigger_index is None:
            if search_outcome is None:
                search_outcome = report.outcome
                search_target = i
                search_halted = report.inspected_index
            continue
        unified, conf = seg.unify_verbose(rule.supply)
        if unified is None:
            if conflict is None:
                conflict = conf
                search_outcome = SearchOutcome.LICENSED
                search_target = i
                search_trigger = report.trigger_index
            continue
        if unified.features == seg.features:
            if search_outcome is None:
                search_outcome = SearchOutcome.LICENSED
                search_target = i
                search_trigger = report.trigger_index
            continue
        new_word[i] = unified.with_provenance(rule.name)
        fired = True
        search_outcome = SearchOutcome.LICENSED
        search_target = i
        search_trigger = report.trigger_index
        break   # first-match-wins per UnificationRule._apply_search

    trace = RuleTrace(
        rule_name=rule.name,
        fired=fired,
        target_matches=tuple(matches),
        search_outcome=search_outcome,
        search_target_index=search_target,
        search_trigger_index=search_trigger,
        search_halted_at=search_halted,
        unify_conflict=conflict if not fired else None,
    )
    return trace, tuple(new_word)


def _explain_deletion(
    rule: DeletionRule,
    word: Word,
) -> tuple[RuleTrace, Word]:
    matches: list[TargetMatch] = []
    fired = False
    search_outcome: SearchOutcome | None = None
    search_target: int | None = None
    search_trigger: int | None = None
    search_halted: int | None = None
    new_word = list(word)
    for i, seg in enumerate(word):
        m = seg.matches(rule.target)
        matches.append(TargetMatch(i, seg.label, m))
        if not m:
            continue
        report = rule.search.locate_verbose(word, i)
        if report.trigger_index is None:
            if search_outcome is None:
                search_outcome = report.outcome
                search_target = i
                search_halted = report.inspected_index
            continue
        cleared = seg.delete(rule.clear)
        if cleared.features == seg.features:
            continue
        new_word[i] = cleared.with_provenance(rule.name)
        fired = True
        search_outcome = SearchOutcome.LICENSED
        search_target = i
        search_trigger = report.trigger_index
        break
    trace = RuleTrace(
        rule_name=rule.name,
        fired=fired,
        target_matches=tuple(matches),
        search_outcome=search_outcome,
        search_target_index=search_target,
        search_trigger_index=search_trigger,
        search_halted_at=search_halted,
    )
    return trace, tuple(new_word)


def _explain_glide(
    rule: GlideFormationRule,
    word: Word,
) -> tuple[RuleTrace, Word]:
    notes: list[str] = []
    fired = False
    if not word:
        notes.append("empty word — nothing to desyllabify")
    else:
        last = word[-1]
        if last.label != rule.target_label:
            notes.append(
                f"last segment is /{last.label}/, not /{rule.target_label}/"
            )
        elif len(word) < 2:
            notes.append(
                f"word has only /{last.label}/; no preceding segment"
            )
        else:
            prev = word[-2]
            if not prev.matches(rule.requires_preceding):
                notes.append(
                    f"preceding /{prev.label}/ fails requires_preceding "
                    f"{dict(rule.requires_preceding)}"
                )
            else:
                notes.append(
                    f"preceding /{prev.label}/ matches; desyllabifying /i/"
                )
                fired = True
    result = rule.apply(word)
    return RuleTrace(
        rule_name=rule.name,
        fired=fired,
        notes=tuple(notes),
    ), result.word


def _explain_allomorph(
    rule: AllomorphSelection,
    word: Word,
    pipeline: RulePipeline,
) -> tuple[RuleTrace, Word]:
    notes: list[str] = []
    idx = None
    for i, seg in enumerate(word):
        if is_placeholder(seg) and seg.label == "@" + rule.morpheme:
            idx = i
            break
    chosen: str | None = None
    elsewhere = False
    if idx is None:
        notes.append(f"no placeholder for morpheme {rule.morpheme!r}")
    else:
        stem_final = word[idx - 1] if idx > 0 else None
        for choice in rule.options:
            if stem_final is not None and stem_final.matches(choice.condition):
                chosen = choice.label
                if not choice.condition:
                    elsewhere = True
                notes.append(
                    f"chose {choice.label} on stem-final "
                    f"/{stem_final.label}/"
                )
                break
        else:
            notes.append("no clause matched — allomorph selection failed")
    result = rule.apply_to(word, pipeline.resolve)
    return RuleTrace(
        rule_name=rule.name,
        fired=bool(result.edits),
        allomorph_chosen=chosen,
        allomorph_elsewhere=elsewhere,
        notes=tuple(notes),
    ), result.word


# ---------------------------------------------------------------------------
# Human-readable rendering
# ---------------------------------------------------------------------------


def format_trace(trace: RuleTrace) -> str:
    """One-block plain-text rendering of a RuleTrace."""
    header = (
        f"[{trace.rule_name}]  {'fired' if trace.fired else 'no-op'}"
    )
    body: list[str] = [header]
    matched = [m for m in trace.target_matches if m.matched_target]
    if trace.target_matches:
        summary = (
            f"  target matched {len(matched)} of "
            f"{len(trace.target_matches)} segments"
        )
        if matched:
            summary += (
                ": " + ", ".join(
                    f"idx{m.index}(/{m.label}/)" for m in matched[:6]
                )
            )
            if len(matched) > 6:
                summary += f", +{len(matched)-6} more"
        body.append(summary)
    if trace.search_outcome is not None:
        so = trace.search_outcome
        if so is SearchOutcome.LICENSED:
            body.append(
                f"  search: LICENSED — trigger at idx "
                f"{trace.search_trigger_index}"
            )
        elif so is SearchOutcome.BROAD_TERMINATOR_BLOCK:
            body.append(
                f"  search: BLOCKED — broad terminator halted on "
                f"idx {trace.search_halted_at} (fails condition)"
            )
        elif so is SearchOutcome.NARROW_SCAN_EXHAUSTED:
            body.append(
                f"  search: EXHAUSTED — narrow scan reached word edge "
                f"without finding a match"
            )
        elif so is SearchOutcome.EMPTY_SCAN:
            body.append("  search: EMPTY — target at word edge")
    if trace.unify_conflict is not None:
        feat, want, got = trace.unify_conflict
        body.append(
            f"  unify: CONFLICT on {feat}: supply wants {want!r}, "
            f"segment has {got!r}"
        )
    if trace.allomorph_chosen is not None:
        marker = " (elsewhere)" if trace.allomorph_elsewhere else ""
        body.append(f"  chose allomorph: {trace.allomorph_chosen}{marker}")
    for note in trace.notes:
        body.append(f"  {note}")
    return "\n".join(body)


def format_explanation(explanation: Explanation) -> str:
    return "\n\n".join(
        format_trace(t) for t in explanation.rule_traces
    )


__all__: Final[tuple[str, ...]] = (
    "Explanation",
    "RuleTrace",
    "TargetMatch",
    "explain_derivation",
    "format_explanation",
    "format_trace",
)
