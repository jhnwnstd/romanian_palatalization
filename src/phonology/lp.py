"""Logical Phonology primitives.

The Romanian palatalization analysis in ``latex.tex`` is stated in
Logical Phonology (LP) notation (Bale & Reiss and successors). This
module makes the LP formalism explicit as first-class Python types so
the framework is auditably LP-compliant AND general enough to
implement any analysis translated into LP primitives.

The LP ontology, mapped to this module:

  - **Feature**: a name from a universal, finite, innate set. Modelled
    as ``str``.
  - **Value**: ``Value = Literal["+", "-"]`` — exactly two. There is
    no third value: **absence is the lack of a value, not a third
    value**.
  - **Valued feature**: ``(value, feature)`` pair. Written ``+F`` /
    ``-F`` in the paper; here ``("+", "F")``.
  - **Segment** (set formulation): a consistent set of valued
    features. Consistency = no feature appears with both values.
    Represented as :class:`~phonology.segments.Segment` whose
    ``features`` map encodes the function formulation ``σ: F → {∅,
    +, -}`` with ``"0"`` for the absence value ``∅``.
  - **Natural class**: :class:`NaturalClass`. Defined intensionally
    by a feature spec ``C`` and contains every segment whose
    valued-feature set is a **superset** of ``C``. Note the direction:
    the class is defined by the *smaller* spec and contains the
    *larger* segments (subsumption).
  - **Feature change**: :class:`FeatureChange`. The ``{ }`` side of a
    rule — a bare set of valued features to be added or removed by an
    operator. Distinct from a natural class (``[ ]``) at the type
    level so a rule box can't confuse the two.
  - **Operators from Ω**: :func:`unify` (``⊔``), :func:`subtract`
    (``\\``). Segment insertion / deletion is handled by
    :class:`SegmentDeletionRule` in :mod:`phonology.rules`.

Bracket conventions the LP paper insists on and this module preserves:

  - ``[ ]`` in the paper → :class:`NaturalClass` in Python.
  - ``{ }`` in the paper → :class:`FeatureChange` in Python.

Circle notation ``[Ⓢ]`` — the natural class defined by segment σ's
whole feature bundle — is :meth:`NaturalClass.from_segment`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final, Literal, Mapping, cast

from .segments import UNSPEC, FeatureMap, Segment

Value = Literal["+", "-"]
ValuedFeature = tuple[Value, str]


# ---------------------------------------------------------------------------
# Set formulation ↔ function formulation
# ---------------------------------------------------------------------------


def valued_features(segment: Segment) -> frozenset[ValuedFeature]:
    """Return σ as an LP set of valued features.

    The framework internally uses the function formulation (a
    dict mapping feature names to ``"+"``, ``"-"``, or ``"0"``); this
    projects out the absence-marker and returns the equivalent set
    formulation of section 1 of the LP notes.
    """
    return frozenset(
        (cast(Value, v), f)
        for f, v in segment.features.items()
        if v in ("+", "-")
    )


def is_consistent(segment: Segment) -> bool:
    """True iff σ carries no feature with both values.

    LP well-formedness: a segment is a *consistent* set of valued
    features. Internally we encode absence as ``"0"``, so if all
    values in ``features.values()`` are in ``{+, -, 0}`` this holds
    trivially. This helper exists so a caller constructing a Segment
    from raw data can assert LP well-formedness explicitly.
    """
    for v in segment.features.values():
        if v not in ("+", "-", "0"):
            return False
    return True


# ---------------------------------------------------------------------------
# Natural class — the [ ] type
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class NaturalClass:
    """A natural class ``N(C)`` — the set of all segments whose
    valued-feature set is a superset of ``C``.

    Class membership is decided by :meth:`contains`; class-to-class
    subsumption by :meth:`subsumes`. The empty specification gives
    the universal class ``N(∅)`` (see :meth:`universal`) — which is
    also the "adjacency terminator" in SEARCH (any segment halts).

    Circle notation ``[Ⓢ]`` — the class defined by segment σ's whole
    bundle — is :meth:`from_segment`.
    """

    spec: FeatureMap

    def contains(self, segment: Segment) -> bool:
        """True iff σ ∈ N(self.spec).

        Subsumption is asymmetric: for each explicit valued feature
        (F, v) in ``spec``, σ must carry (F, v). Features unspecified
        in ``spec`` (value ``"0"``) impose no constraint. A segment
        that is UNSPECIFIED for a feature the class specifies is NOT
        in the class — LP's central "cannot target absence" rule.
        """
        for feat, want in self.spec.items():
            if want == UNSPEC:
                continue
            if segment.features.get(feat, UNSPEC) != want:
                return False
        return True

    def subsumes(self, other: "NaturalClass") -> bool:
        """True iff N(self) ⊇ N(other) — i.e., every segment in
        ``other``'s class is also in ``self``'s class.

        This holds iff ``self.spec ⊆ other.spec`` (as sets of valued
        features): more features means a smaller class. This is the
        "Phonological Sawzall" corollary — you cannot carve out a
        less-specified segment from under a more-specified one.
        """
        for feat, want in self.spec.items():
            if want == UNSPEC:
                continue
            if other.spec.get(feat, UNSPEC) != want:
                return False
        return True

    def is_consistent(self) -> bool:
        """True iff the spec doesn't require F to be both + and −.

        Since we encode a spec as a mapping from feature name to a
        single value, inconsistency at the spec level is unrepresentable
        by construction — but a caller composing specs with dict
        merge may still produce nonsense values; this method surfaces
        that.
        """
        for v in self.spec.values():
            if v not in ("+", "-", "0"):
                return False
        return True

    @classmethod
    def universal(cls) -> "NaturalClass":
        """The universal class ``N(∅)`` — contains every segment.

        In SEARCH: ``Trm = universal()`` gives adjacency (every
        segment halts the scan). "Adjacency is opaqueness."
        """
        return cls(spec={})

    @classmethod
    def from_segment(cls, segment: Segment) -> "NaturalClass":
        """Circle notation ``[Ⓢ]``: the natural class defined by
        σ's entire feature bundle.

        Every underspecified feature is dropped (absence isn't part of
        the defining spec). For a fully specified segment the class
        is a singleton; for a partially specified one it contains
        every more-specified segment.
        """
        spec = {f: v for f, v in segment.features.items() if v in ("+", "-")}
        return cls(spec=spec)

    def to_dict(self) -> dict[str, str]:
        """Return the spec as a plain dict — useful for legacy rule
        constructors that still take ``Mapping[str, str]``."""
        return dict(self.spec)


# ---------------------------------------------------------------------------
# Feature change — the { } type
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class FeatureChange:
    """A set of valued features to unify into a segment (``⊔``) or a
    set of feature *names* to strip from it (``\\``).

    LP distinguishes the ``{ }`` side of a rule from the ``[ ]`` side
    (natural class). This dataclass makes the distinction explicit at
    the type level — a rule constructor that mixes them up is a
    type error to the reader, matching the paper's non-negotiable
    bracket convention.

    For unification changes, ``add`` carries the valued features.
    For subtraction changes, ``remove`` carries the feature *names*
    (LP's ``A \\ {+F, -F}`` = strip F at whatever value it bears).
    """

    add: Mapping[str, str] | None = None
    remove: frozenset[str] = frozenset()

    def __post_init__(self) -> None:
        # Frozen dataclass — bypass immutability to normalise a
        # None default to an empty mapping so downstream callers
        # can dispatch on ``add`` without a None-check.
        if self.add is None:
            object.__setattr__(self, "add", {})


# ---------------------------------------------------------------------------
# Operators from Ω
# ---------------------------------------------------------------------------


def unify(segment: Segment, change: Mapping[str, str]) -> Segment | None:
    """LP unification ``σ ⊔ change``. Partial operation.

    - If σ is open for F, adopt the change's value.
    - If σ has the same value, no-op.
    - If σ has the opposite value, unification fails: return None.

    A caller building a total rule wraps a None back to identity —
    the rule is a total function even though the operator is not.
    Currently equivalent to :meth:`Segment.unify`; the top-level
    function form matches LP's ``⊔`` operator notation for callers
    reasoning at the ontology level rather than the class level.
    """
    return segment.unify(change)


def subtract(segment: Segment, feature_names: frozenset[str]) -> Segment:
    """LP subtraction ``σ \\ {+F, -F, ...}`` — strip each named
    feature at whatever value it bears.

    Equivalent to projecting σ onto the complement of the named
    feature set. Total: never fails, since removing a feature σ
    doesn't have is a no-op. This is the operator the bleed rule
    uses to delink postalveolar place from a coronal.
    """
    return segment.delete(feature_names)


# ---------------------------------------------------------------------------
# Convenience — the "bracket" constructors, for readability
# ---------------------------------------------------------------------------


def natural_class(**features: str) -> NaturalClass:
    """Ergonomic constructor: ``natural_class(CORONAL='+', Anterior='+')``.

    Reads like the paper's ``[+Coronal, +Anterior]`` notation while
    staying in plain Python — no DSL required.
    """
    return NaturalClass(spec=dict(features))


def feature_change(
    add: Mapping[str, str] | None = None,
    remove: frozenset[str] | None = None,
) -> FeatureChange:
    """Ergonomic constructor for a ``{ }`` change."""
    return FeatureChange(
        add=dict(add) if add else {},
        remove=remove if remove is not None else frozenset(),
    )


__all__: Final[tuple[str, ...]] = (
    "FeatureChange",
    "NaturalClass",
    "Value",
    "ValuedFeature",
    "feature_change",
    "is_consistent",
    "natural_class",
    "subtract",
    "unify",
    "valued_features",
)
