"""Pluggable phonological-rule framework.

A small unification-based phonology engine: feature-bundled
:class:`Segment` s composed into a :class:`Word`, rules that target
segments by feature pattern, and a SEARCH primitive parametrised by
direction, terminator-breadth, and condition. The engine exported
here is language-neutral: :func:`tokenize_ipa` takes a
:class:`TokenizerConfig`, and only its default,
:data:`ROMANIAN_TOKENIZER`, drops the length mark ː and folds ɛ/ɔ to
e/o for the Romanian inventory. The language-specific steps of a
validation run (UR construction, lexicon row loading, cluster-report
buckets, the predicted-vs-attested comparison convention, mismatch
categorisation, IPA tokenization) are hooks on
``phonology.profile.AnalysisProfile``, whose defaults are the
Romanian implementations in ``g2p``, ``lexicon``, ``diagnostics``
and ``segments``. Feature inventories are found by
:func:`~phonology.inventory.inventory_path` (package data in
``phonology/inventories/``). The analysis itself lives in
``phonology.analyses.romanian_palatalization``.
"""

from __future__ import annotations

from .inventory import (
    FeatureInventory,
    FeaturePatch,
    UnderspecifiedSegment,
    inventory_path,
    load_inventory_file,
)
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
from .segments import (
    ROMANIAN_TOKENIZER,
    LengthHandling,
    Segment,
    TokenizerConfig,
    Word,
    project,
    segments_to_ipa,
    tokenize_ipa,
)

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
    "LengthHandling",
    "ROMANIAN_TOKENIZER",
    "Rule",
    "RulePipeline",
    "Search",
    "Segment",
    "Step",
    "TokenizerConfig",
    "UnderspecifiedSegment",
    "UnificationRule",
    "Word",
    "inventory_path",
    "load_inventory_file",
    "placeholder",
    "project",
    "segments_to_ipa",
    "tokenize_ipa",
)
