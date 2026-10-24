"""The feature inventory: base segments from JSON plus declared patches.

The JSON (``romanian_features.json`` for Romanian) is the *base*: it
contains every surface segment of the language with its feature
specifications. The analysis module declares two layers on top:

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

Inventory files are found by language name: :func:`inventory_path`
resolves ``<name>_features.json`` in the package data directory
:data:`INVENTORY_DIR` (``phonology/inventories/``, whose README
documents the file format), and :func:`load_inventory_file` loads it
with an analysis's patches. ``romanian_features.json`` predates that
directory and stays at the repository root, where the lookup for
``"romanian"`` falls back to it.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Final, Mapping, Sequence

from .segments import UNSPEC, FeatureMap, FeatureValue, Segment

# ---------------------------------------------------------------------------
# Inventory files
# ---------------------------------------------------------------------------

INVENTORY_DIR: Final[Path] = Path(__file__).resolve().parent / "inventories"
"""Package data directory holding the ``<name>_features.json`` files."""

# ``*.json`` is gitignored but ``!*_features.json`` is re-included, so
# every inventory file name ends in this suffix.
_INVENTORY_SUFFIX: Final[str] = "_features.json"

# A name is a lowercase identifier ("italian", "czech"), so the file it
# names can't lie outside the directory it is looked up in.
_INVENTORY_NAME: Final[re.Pattern[str]] = re.compile(r"[a-z][a-z0-9_]*")

# Inventories that predate INVENTORY_DIR, at their original location.
# The default lookup falls back to these after INVENTORY_DIR.
_LEGACY_INVENTORY_FILES: Final[Mapping[str, Path]] = {
    "romanian": Path(__file__).resolve().parents[2] / "romanian_features.json",
}


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


def inventory_path(name: str, *, directory: Path | None = None) -> Path:
    """Resolve the inventory file ``<name>_features.json``.

    The default lookup tries :data:`INVENTORY_DIR`, then, for an
    inventory that predates it (``"romanian"``), its original location
    at the repository root. ``directory`` replaces that lookup with
    the one directory given (tests pass a temporary one).

    Raises :class:`ValueError` if ``name`` isn't a lowercase
    identifier, and :class:`FileNotFoundError` naming every path tried
    if none of them is a file.
    """
    if not _INVENTORY_NAME.fullmatch(name):
        raise ValueError(
            "inventory names are lowercase identifiers such as "
            f"'italian'; got {name!r}"
        )
    filename = name + _INVENTORY_SUFFIX
    candidates: list[Path] = []
    if directory is not None:
        candidates.append(Path(directory) / filename)
    else:
        candidates.append(INVENTORY_DIR / filename)
        legacy = _LEGACY_INVENTORY_FILES.get(name)
        if legacy is not None:
            candidates.append(legacy)
    for path in candidates:
        if path.is_file():
            return path
    raise FileNotFoundError(
        f"no inventory {filename!r}; looked for "
        + ", ".join(str(p) for p in candidates)
    )


def load_inventory_file(
    name: str,
    patches: Sequence[FeaturePatch] = (),
    underspec: Sequence[UnderspecifiedSegment] = (),
    *,
    directory: Path | None = None,
) -> FeatureInventory:
    """Load the inventory called ``name`` with the given patches.

    Resolves the file with :func:`inventory_path` (same ``directory``
    semantics) and loads it with :meth:`FeatureInventory.load`. An
    analysis that also needs the path, for
    :attr:`~phonology.profile.AnalysisProfile.inventory_json`, calls
    :func:`inventory_path` itself.
    """
    return FeatureInventory.load(
        inventory_path(name, directory=directory),
        patches=patches,
        underspec=underspec,
    )


__all__: Final[tuple[str, ...]] = (
    "INVENTORY_DIR",
    "FeatureInventory",
    "FeaturePatch",
    "UnderspecifiedSegment",
    "inventory_path",
    "load_inventory_file",
)
