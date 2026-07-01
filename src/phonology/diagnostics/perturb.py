"""Single-feature perturbation search.

If a derivation fails, is there a SMALL change — one feature value
somewhere — that would make it succeed? This module answers the
question by enumerating single-feature edits over four axes and
running the whole pipeline for each perturbation:

  - ``PATCH_VALUE`` — flip the value of one field in one
    :class:`FeaturePatch` (e.g., ``/i/ CORONAL "+" → "-"``).
  - ``CLEAR_ADD`` — add one feature to a
    :class:`UnderspecifiedSegment`'s clear set.
  - ``CLEAR_REMOVE`` — remove one feature from a clear set.
  - ``RULE_SUPPLY`` — flip the value of one supply feature in a
    :class:`UnificationRule`.

A successful perturbation is a smoking-gun diagnostic: "the analysis
would produce the right SR for this lemma if only /X/ had feature F
set to V". Perturbations that succeed on ONE lemma but break many
others are still worth naming — the caller may want to test the
perturbation lexicon-wide.

Only :class:`UnificationRule` supply perturbations are enumerated —
DeletionRule and GlideFormationRule have less structured feature
maps and would need a separate perturbation kind. This is a
deliberate v1 scope choice; extensions can add more kinds later.
"""

from __future__ import annotations

from dataclasses import dataclass
from dataclasses import replace as dc_replace
from enum import StrEnum
from pathlib import Path
from typing import TYPE_CHECKING, Callable, Final, Iterator, Sequence

from ..inventory import FeatureInventory, FeaturePatch, UnderspecifiedSegment
from ..pipeline import RulePipeline
from ..rules import Rule, UnificationRule
from ..segments import Word, segments_to_ipa
from .compare import compare_ipa

if TYPE_CHECKING:
    from ..analyses.romanian_palatalization import AnalysisProfile
    from ..lexicon import LexRow


class PerturbationKind(StrEnum):
    PATCH_VALUE = "PATCH_VALUE"
    CLEAR_ADD = "CLEAR_ADD"
    CLEAR_REMOVE = "CLEAR_REMOVE"
    RULE_SUPPLY = "RULE_SUPPLY"


@dataclass(frozen=True, slots=True)
class Perturbation:
    """One single-feature edit."""

    kind: PerturbationKind
    target: str  # patch segment / underspec label / rule name
    feature: str  # which feature is being edited
    old_value: str  # or "" for CLEAR_ADD
    new_value: str  # or "" for CLEAR_REMOVE

    def describe(self) -> str:
        if self.kind is PerturbationKind.PATCH_VALUE:
            return (
                f"patch {self.target!r} {self.feature}: "
                f"{self.old_value!r} → {self.new_value!r}"
            )
        if self.kind is PerturbationKind.CLEAR_ADD:
            return f"underspec {self.target!r} clear += {{{self.feature}}}"
        if self.kind is PerturbationKind.CLEAR_REMOVE:
            return f"underspec {self.target!r} clear -= {{{self.feature}}}"
        return (
            f"rule {self.target!r} supply {self.feature}: "
            f"{self.old_value!r} → {self.new_value!r}"
        )


@dataclass(frozen=True, slots=True)
class PerturbationResult:
    """Outcome of one perturbation on one lemma."""

    perturbation: Perturbation
    sr: str
    matched: bool
    edit_distance: int


@dataclass(frozen=True, slots=True)
class PerturbationReport:
    """Result of a perturbation search on one lemma."""

    baseline_sr: str
    baseline_matched: bool
    baseline_distance: int
    tried: int
    successful: tuple[PerturbationResult, ...]
    near_misses: tuple[PerturbationResult, ...]  # closer than baseline


_FLIP_VALUE = {"+": "-", "-": "+", "0": "+"}


def _enumerate_perturbations(
    patches: Sequence[FeaturePatch],
    underspec: Sequence[UnderspecifiedSegment],
    rules: Sequence[Rule],
    inventory: FeatureInventory,
    kinds: Sequence[PerturbationKind],
) -> "Iterator[Perturbation]":
    """Yield candidate :class:`Perturbation` objects."""
    all_features = tuple(inventory.features)
    if PerturbationKind.PATCH_VALUE in kinds:
        for patch in patches:
            for feat, val in patch.updates.items():
                new = _FLIP_VALUE.get(val, "+")
                if new != val:
                    yield Perturbation(
                        PerturbationKind.PATCH_VALUE,
                        patch.segment,
                        feat,
                        val,
                        new,
                    )
    if PerturbationKind.CLEAR_ADD in kinds:
        for u in underspec:
            for feat in all_features:
                if feat in u.clear:
                    continue
                yield Perturbation(
                    PerturbationKind.CLEAR_ADD,
                    u.label,
                    feat,
                    "",
                    "",
                )
    if PerturbationKind.CLEAR_REMOVE in kinds:
        for u in underspec:
            for feat in sorted(u.clear):
                yield Perturbation(
                    PerturbationKind.CLEAR_REMOVE,
                    u.label,
                    feat,
                    "",
                    "",
                )
    if PerturbationKind.RULE_SUPPLY in kinds:
        for r in rules:
            if not isinstance(r, UnificationRule):
                continue
            for feat, val in r.supply.items():
                new = _FLIP_VALUE.get(val, "+")
                if new != val:
                    yield Perturbation(
                        PerturbationKind.RULE_SUPPLY,
                        r.name,
                        feat,
                        val,
                        new,
                    )


def _apply_perturbation(
    p: Perturbation,
    patches: tuple[FeaturePatch, ...],
    underspec: tuple[UnderspecifiedSegment, ...],
    rules: tuple[Rule, ...],
) -> tuple[
    tuple[FeaturePatch, ...],
    tuple[UnderspecifiedSegment, ...],
    tuple[Rule, ...],
]:
    """Return new (patches, underspec, rules) tuples with the
    perturbation applied. Originals are not mutated."""
    if p.kind is PerturbationKind.PATCH_VALUE:
        new_patches = list(patches)
        for i, patch in enumerate(new_patches):
            if patch.segment == p.target and p.feature in patch.updates:
                new_updates = dict(patch.updates)
                new_updates[p.feature] = p.new_value
                new_patches[i] = dc_replace(patch, updates=new_updates)
                break
        return tuple(new_patches), underspec, rules
    if p.kind in (PerturbationKind.CLEAR_ADD, PerturbationKind.CLEAR_REMOVE):
        new_underspec = list(underspec)
        for i, u in enumerate(new_underspec):
            if u.label == p.target:
                clear = set(u.clear)
                if p.kind is PerturbationKind.CLEAR_ADD:
                    clear.add(p.feature)
                else:
                    clear.discard(p.feature)
                new_underspec[i] = dc_replace(u, clear=frozenset(clear))
                break
        return patches, tuple(new_underspec), rules
    if p.kind is PerturbationKind.RULE_SUPPLY:
        new_rules = list(rules)
        for i, r in enumerate(new_rules):
            if r.name == p.target and isinstance(r, UnificationRule):
                new_supply = dict(r.supply)
                new_supply[p.feature] = p.new_value
                new_rules[i] = dc_replace(r, supply=new_supply)
                break
        return patches, underspec, tuple(new_rules)
    return patches, underspec, rules


def search_perturbations(
    ur_builder: Callable[[FeatureInventory], Word | None],
    expected_field: str,
    inventory_json: Path,
    baseline_patches: Sequence[FeaturePatch],
    baseline_underspec: Sequence[UnderspecifiedSegment],
    baseline_rules: Sequence[Rule],
    kinds: Sequence[PerturbationKind] = tuple(PerturbationKind),
    limit: int | None = None,
) -> PerturbationReport:
    """Enumerate single-feature perturbations and report which flip the
    derivation from failing to succeeding.

    ``ur_builder`` is a function ``inventory -> Word`` — different
    perturbations rebuild the inventory (patches/underspec changes),
    so the UR is rebuilt each time by re-invoking this callable. This
    lets the caller keep control over how the UR is derived (which
    may itself depend on the inventory features).
    """
    patches = tuple(baseline_patches)
    underspec = tuple(baseline_underspec)
    rules = tuple(baseline_rules)

    def _run(
        p: tuple[FeaturePatch, ...],
        u: tuple[UnderspecifiedSegment, ...],
        r: tuple[Rule, ...],
    ) -> tuple[str, bool, int]:
        try:
            inv = FeatureInventory.load(inventory_json, p, u)
        except Exception:
            return "", False, 999
        ur = ur_builder(inv)
        if ur is None:
            return "", False, 999
        pipe = RulePipeline(
            rules=r,
            resolve=lambda label: inv.segment(label),
        )
        try:
            deriv = pipe.derive(ur)
        except Exception:
            return "", False, 999
        sr = segments_to_ipa(deriv.sr, inv.base_segments())
        cmp = compare_ipa(sr, expected_field)
        dist = (
            int(cmp.edit_distance)
            if cmp.edit_distance != float("inf")
            else 999
        )
        return sr, cmp.matched, dist

    baseline_sr, baseline_matched, baseline_dist = _run(
        patches, underspec, rules
    )
    successful: list[PerturbationResult] = []
    near_misses: list[PerturbationResult] = []
    tried = 0

    for perturbation in _enumerate_perturbations(
        patches,
        underspec,
        rules,
        FeatureInventory.load(
            inventory_json,
            patches,
            underspec,
        ),
        kinds,
    ):
        new_p, new_u, new_r = _apply_perturbation(
            perturbation,
            patches,
            underspec,
            rules,
        )
        sr, matched, dist = _run(new_p, new_u, new_r)
        result = PerturbationResult(
            perturbation=perturbation,
            sr=sr,
            matched=matched,
            edit_distance=dist,
        )
        if matched:
            successful.append(result)
        elif dist < baseline_dist:
            near_misses.append(result)
        tried += 1
        if limit is not None and tried >= limit:
            break

    successful.sort(key=lambda r: r.edit_distance)
    near_misses.sort(key=lambda r: r.edit_distance)

    return PerturbationReport(
        baseline_sr=baseline_sr,
        baseline_matched=baseline_matched,
        baseline_distance=baseline_dist,
        tried=tried,
        successful=tuple(successful),
        near_misses=tuple(near_misses[:20]),
    )


def format_report(report: PerturbationReport) -> str:
    """Human-readable rendering."""
    lines: list[str] = []
    lines.append(
        f"Baseline: /{report.baseline_sr}/ "
        f"{'✓' if report.baseline_matched else '✗'}  "
        f"distance={report.baseline_distance}"
    )
    lines.append(f"Perturbations tried: {report.tried}")
    if report.successful:
        lines.append(f"  ✓ {len(report.successful)} perturbation(s) succeed:")
        for r in report.successful[:10]:
            lines.append(f"     ✓ {r.perturbation.describe()}")
        if len(report.successful) > 10:
            lines.append(f"     ... and {len(report.successful) - 10} more")
    else:
        lines.append("  no single-feature perturbation succeeds")
        if report.near_misses:
            lines.append("  closest near-misses (reduced distance):")
            for r in report.near_misses[:5]:
                lines.append(
                    f"     ▲ {r.perturbation.describe()}  → /{r.sr}/ "
                    f"distance={r.edit_distance}"
                )
    return "\n".join(lines)


def search_perturbations_for(
    profile: "AnalysisProfile",
    row: "LexRow",
    *,
    kinds: Sequence[PerturbationKind] = tuple(PerturbationKind),
    limit: int | None = None,
) -> PerturbationReport:
    """Ergonomic wrapper: takes an :class:`AnalysisProfile` and a
    :class:`LexRow`, unpacks the ur-builder closure and the baseline
    (patches, underspec, rules) tuples for you.

    Prefer this over :func:`search_perturbations` in new code —
    the 7-kwarg surface of the low-level function is unfriendly to
    interactive use.
    """
    from ..g2p import build_ur

    def builder(inv: FeatureInventory) -> Word | None:
        return build_ur(
            row.ipa_lemma,
            row.opportunity,
            inv,
            stem_final=row.stem_final,
        )

    return search_perturbations(
        ur_builder=builder,
        expected_field=row.ipa_pl,
        inventory_json=profile.inventory_json,
        baseline_patches=profile.patches,
        baseline_underspec=profile.underspec,
        baseline_rules=profile.rules,
        kinds=kinds,
        limit=limit,
    )


__all__: Final[tuple[str, ...]] = (
    "Perturbation",
    "PerturbationKind",
    "PerturbationReport",
    "PerturbationResult",
    "format_report",
    "search_perturbations",
    "search_perturbations_for",
)
