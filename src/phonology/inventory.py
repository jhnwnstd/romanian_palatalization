"""The feature inventory: base segments from JSON plus declared patches.

The JSON at ``romanian_features.json`` is the *base*: it contains every
surface segment of Romanian with its feature specifications. The
analysis module declares two layers on top:

  - ``FeaturePatch``: a per-segment feature override. The paper claims
    ``/i/`` is "coronal" because a high front vocoid is coronal — but
    the base JSON has ``/i/`` as ``CORONAL: "-"`` to keep ``DORSAL`` on
    consonantal segments. The patch is how that disagreement is
    encoded *in the analysis* rather than by editing the inventory.
  - ``UnderspecifiedSegment``: an underspecified variant like ``/K/`` —
    a copy of an existing surface segment with selected features
    cleared to ``"0"``. This is the inalterability-as-prespecification
    mechanism: ``/K/`` is ``/k/`` minus the place features the rules
    write, so the rules can write into it; ``/k/`` keeps those
    features fixed and is therefore inalterable.

Patches and underspecifications are pure data — they sit in the
analysis module and the inventory loader applies them without
interpreting them. Swapping one out for another is a one-line edit in
the analysis module, never an edit to this file or to the JSON.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Final, Mapping, Sequence

from .segments import UNSPEC, FeatureMap, FeatureValue, Segment


@dataclass(frozen=True, slots=True)
class FeaturePatch:
    """An override applied to one base inventory segment.

    ``updates`` is unified into the base segment's feature dict by
    overwriting: an explicit ``"+"`` or ``"-"`` in ``updates`` replaces
    whatever was there. Use sparingly — the patch list is the
    auditable record of how the analysis diverges from the inventory.
    """

    segment: str
    updates: Mapping[str, FeatureValue]
    reason: str


@dataclass(frozen=True, slots=True)
class UnderspecifiedSegment:
    """An underspecified variant of an existing inventory segment.

    ``clear`` lists every feature the analysis treats as open on this
    variant. The list must include *every* feature any rule may unify
    into the segment — otherwise the rule will hit a fixed contrary
    value and fail. For ``/K/`` this means clearing CORONAL, DORSAL,
    Anterior, Distributed, Strident, Continuant, and DelRel — the full
    bundle dorsal palatalization writes.
    """

    label: str
    base: str
    clear: frozenset[str]
    reason: str


@dataclass(frozen=True, slots=True)
class FeatureInventory:
    """Loaded inventory plus declared patches and underspecs.

    Build once with :meth:`load`, query with :meth:`segment` /
    :meth:`features_of`. The base mapping is preserved so renderers
    can project feature bundles back to inventory labels.
    """

    features: tuple[str, ...]
    base: Mapping[str, FeatureMap]
    underspec: Mapping[str, UnderspecifiedSegment]

    @classmethod
    def load(
        cls,
        json_path: Path,
        patches: Sequence[FeaturePatch] = (),
        underspec: Sequence[UnderspecifiedSegment] = (),
    ) -> FeatureInventory:
        """Load the JSON inventory and apply the declared patches.

        Patches are applied in order; later patches on the same segment
        overwrite earlier ones (you should never declare overlapping
        patches in the first place, but the rule is documented).
        """
        with json_path.open(encoding="utf-8") as f:
            raw = json.load(f)
        features = tuple(raw["features"])
        base: dict[str, dict[str, FeatureValue]] = {
            label: dict(feats) for label, feats in raw["segments"].items()
        }
        for patch in patches:
            if patch.segment not in base:
                raise KeyError(
                    f"FeaturePatch refers to unknown segment "
                    f"{patch.segment!r}: {patch.reason}"
                )
            base[patch.segment].update(patch.updates)
        # Validate underspec base references; we don't bake the variant
        # into ``base`` because UR construction (g2p) decides when to
        # swap a surface segment for its underspecified twin.
        for u in underspec:
            if u.base not in base:
                raise KeyError(
                    f"UnderspecifiedSegment {u.label!r} references "
                    f"unknown base {u.base!r}"
                )
        underspec_map = {u.label: u for u in underspec}
        return cls(
            features=features,
            base={k: dict(v) for k, v in base.items()},
            underspec=underspec_map,
        )

    def features_of(self, label: str) -> FeatureMap:
        """Return the feature map for a segment label.

        If the label is an underspecified variant (declared via
        :class:`UnderspecifiedSegment`), the base segment's features
        are returned with the variant's ``clear`` set re-opened to
        ``"0"``. Otherwise the base mapping is returned as-is.
        """
        if label in self.underspec:
            u = self.underspec[label]
            base = self.base[u.base]
            return {
                feat: (UNSPEC if feat in u.clear else val)
                for feat, val in base.items()
            }
        return dict(self.base[label])

    def segment(self, label: str) -> Segment:
        """Construct a ``Segment`` for ``label`` with its full features."""
        return Segment(label=label, features=self.features_of(label))

    def base_segments(self) -> Mapping[str, FeatureMap]:
        """Read-only view of the (patched) base inventory.

        Used by :func:`phonology.segments.project` to render abstract
        feature bundles back to IPA labels.
        """
        return self.base


__all__: Final[tuple[str, ...]] = (
    "FeatureInventory",
    "FeaturePatch",
    "UnderspecifiedSegment",
)
