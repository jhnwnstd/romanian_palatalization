"""AnalysisProfile: the one object an analysis run passes around.

An analysis module (e.g.
:mod:`phonology.analyses.romanian_palatalization`) builds one
:class:`AnalysisProfile`; the validator (:mod:`phonology.validation`)
and the profile-aware diagnostics
(:func:`~phonology.diagnostics.ordering.search_orderings`,
:func:`~phonology.diagnostics.perturb.search_perturbations_for`) take
it as their only analysis argument. The class used to live in the
Romanian analysis module, which still re-exports it.

Besides the analysis proper (inventory, rules, oracle), a profile
carries six **hooks**: the language-specific steps of a validation
run.

  - ``ur_builder(row, inventory) -> Word | None`` — build the UR for
    one row. Default :func:`phonology.g2p.build_ur_from_row`.
  - ``cluster_key(row) -> str | None`` — the row's bucket in the
    cluster report and the mismatch CSV's ``cluster`` column.
    Default :func:`phonology.lexicon.romanian_cluster_key`.
  - ``normalisers`` — the :class:`NormaliserCascade` deciding when a
    predicted SR matches an attested one. Every comparison made for
    the profile uses it: the IPA match flag, the categoriser, the
    distance metric, and the ordering and perturbation searches.
    Default :data:`phonology.diagnostics.compare.ROMANIAN_NORMALISERS`.
  - ``categoriser(row, predicted_sr, attested_field,
    boolean_predicted, *, ur_built, normalisers)`` — name the cause of
    a mismatch. Default
    :func:`phonology.diagnostics.categorise.refine_category`.
  - ``row_factory(csv_row) -> LexRow`` — read one row of a lexicon
    CSV, so a lexicon with another column schema can feed the
    validator. Default :meth:`phonology.lexicon.LexRow.from_dict`.
  - ``tokenizer`` — the :class:`~phonology.segments.TokenizerConfig`
    that splits IPA into inventory labels wherever the validator or
    a diagnostic tokenizes for the profile: the attested side of the
    distance metric's feature component. Default
    :data:`phonology.segments.ROMANIAN_TOKENIZER`. The ``ur_builder``
    hook is not handed the profile, so a builder tokenizes the lemma
    with a config of its own choosing; it should be this one (the
    Romanian builder uses ``ROMANIAN_TOKENIZER``).

The defaults are the Romanian implementations, so a profile built
without hooks behaves exactly as profiles did before the hooks
existed. Another analysis passes its own at construction, or derives
a variant with ``dataclasses.replace(profile, ur_builder=...)``. The
low-level functions behind the hooks keep the same Romanian defaults
and take the hook values as keyword arguments (``compare_ipa(...,
normalisers=...)``, ``iter_lexicon_rows(..., row_factory=...)``,
``tokenize_ipa(..., tokenizer=...)``). Category and normalisation
labels are shared vocabularies: a hook that needs a new label adds a
member to :class:`~phonology.diagnostics.categorise.RefinedCategory`
or :class:`~phonology.diagnostics.compare.Normalisation`.

Design rule R8 (``docs/multilingual_design.md``) asks for the
tokenizer and compare-strip sets to come from the profile. The
tokenizer's skip set is ``tokenizer.skip``, with ``tokenizer.length``
deciding the length mark ː. The compare-strip set needs no field of
its own: it is
:attr:`~phonology.diagnostics.compare.NormaliserCascade.strip`, part
of ``normalisers``, which every comparison made for the profile
already uses. The Romanian strip removes the stress marks and ː; a
cascade for a language with contrastive length keeps ː.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Final, Protocol

from .diagnostics.categorise import Categorisation, refine_category
from .diagnostics.compare import ROMANIAN_NORMALISERS, NormaliserCascade
from .g2p import build_ur_from_row
from .inventory import FeatureInventory, FeaturePatch, UnderspecifiedSegment
from .lexicon import LexRow, RowFactoryFn, romanian_cluster_key
from .rules import Rule
from .segments import ROMANIAN_TOKENIZER, TokenizerConfig, Word

# ---------------------------------------------------------------------------
# Hook signatures
# ---------------------------------------------------------------------------


class UrBuilderFn(Protocol):
    """Build the underlying representation for one row.

    Returns ``None`` when the row can't be derived (unusable IPA, an
    unknown segment); the validator then scores the row with
    ``profile.fallback_predict``. Must be pure: the perturbation
    search calls it once per perturbed inventory.
    """

    def __call__(
        self,
        row: LexRow,
        inventory: FeatureInventory,
        /,
    ) -> Word | None: ...


class ClusterKeyFn(Protocol):
    """Map a row to its cluster-report bucket, or ``None``."""

    def __call__(self, row: LexRow, /) -> str | None: ...


class CategoriserFn(Protocol):
    """Categorise one row's mismatch.

    Called for every evaluated row, matched or not. ``ur_built`` is
    False when ``ur_builder`` returned ``None``; ``normalisers`` is
    the profile's cascade, to compare under so NORM_* categories
    agree with the validator's IPA match flag.
    """

    def __call__(
        self,
        row: LexRow,
        predicted_sr: str,
        attested_field: str,
        boolean_predicted: bool,
        /,
        *,
        ur_built: bool,
        normalisers: NormaliserCascade,
    ) -> Categorisation: ...


# ---------------------------------------------------------------------------
# AnalysisProfile — the "one object to pass around" bundle
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class AnalysisProfile:
    """Everything an analysis run needs, bundled into one immutable
    object so downstream API functions take a single argument.

    Fields:

      - ``name`` — human-readable analysis name.
      - ``inventory`` — the constructed :class:`FeatureInventory`.
      - ``patches`` / ``underspec`` — the source-of-truth tuples used
        to build ``inventory`` (needed by the perturbation searcher
        to reload the inventory with tweaks).
      - ``rules`` — the ordered rule tuple.
      - ``inventory_json`` — path to the JSON that ``patches`` were
        applied on top of (used by perturb to reload).
      - ``target_stems`` — orthographic stem-finals this analysis
        wants to see; passed to :func:`iter_lexicon_rows`.
      - ``palatalization_rule_names`` — which rule names count as
        "mutation-producing"; passed to distance.score.
      - ``expected_firings`` — the rule-firing oracle (per-row).
      - ``fallback_predict`` — the boolean predictor for rows the
        pipeline can't process.

    Hooks (optional; see the module docstring): ``ur_builder``,
    ``cluster_key``, ``normalisers``, ``categoriser``,
    ``row_factory``, ``tokenizer``.
    :func:`~phonology.validation.trace_lemmas` reads the CSV with
    ``row_factory``; other callers pass
    ``row_factory=profile.row_factory`` to :func:`iter_lexicon_rows`.
    """

    name: str
    inventory: FeatureInventory
    patches: tuple[FeaturePatch, ...]
    underspec: tuple[UnderspecifiedSegment, ...]
    rules: tuple[Rule, ...]
    inventory_json: Path
    target_stems: frozenset[str]
    palatalization_rule_names: frozenset[str]
    expected_firings: Callable[[LexRow], frozenset[str]]
    fallback_predict: Callable[[LexRow], bool]
    ur_builder: UrBuilderFn = build_ur_from_row
    cluster_key: ClusterKeyFn = romanian_cluster_key
    normalisers: NormaliserCascade = ROMANIAN_NORMALISERS
    categoriser: CategoriserFn = refine_category
    row_factory: RowFactoryFn = LexRow.from_dict
    tokenizer: TokenizerConfig = ROMANIAN_TOKENIZER


__all__: Final[tuple[str, ...]] = (
    "AnalysisProfile",
    "CategoriserFn",
    "ClusterKeyFn",
    "NormaliserCascade",
    "RowFactoryFn",
    "TokenizerConfig",
    "UrBuilderFn",
)
