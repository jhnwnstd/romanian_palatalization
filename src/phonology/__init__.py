"""Pluggable phonological-rule framework.

A small unification-based phonology engine: feature-bundled
:class:`Segment` s composed into a :class:`Word`, rules that target
segments by feature pattern, and a SEARCH primitive parametrised by
direction, terminator-breadth, and condition. The framework knows
nothing about Romanian; the analysis lives in
``phonology.analyses.romanian_palatalization``.
"""

from __future__ import annotations

from .inventory import FeatureInventory, FeaturePatch, UnderspecifiedSegment
from .pipeline import Derivation, RulePipeline, Step
from .rules import (
    AllomorphChoice,
    AllomorphSelection,
    ApplyResult,
    DeletionRule,
    Edit,
    GlideFormationRule,
    Rule,
    UnificationRule,
    placeholder,
)
from .search import Direction, Search
from .segments import Segment, Word, project, segments_to_ipa, tokenize_ipa

__all__ = (
    "AllomorphChoice",
    "AllomorphSelection",
    "ApplyResult",
    "DeletionRule",
    "Derivation",
    "Direction",
    "Edit",
    "FeatureInventory",
    "FeaturePatch",
    "GlideFormationRule",
    "Rule",
    "RulePipeline",
    "Search",
    "Segment",
    "Step",
    "UnderspecifiedSegment",
    "UnificationRule",
    "Word",
    "placeholder",
    "project",
    "segments_to_ipa",
    "tokenize_ipa",
)
