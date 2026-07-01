"""Build an underlying representation from a lexicon row.

The lemma's IPA (the ``ipa_normalized_lemma`` column) is the input —
already NFC-normalised with tie bars inserted. The output is a
``Word``: a tuple of fully-featured :class:`Segment`s ready to feed
into the rule pipeline.

Three transformations turn the lemma into a UR:

1. **Strip the thematic vowel.** Feminine lemmas in ``-ă`` (IPA
   ``ə``) and consonant-classed lemmas in ``-e`` carry a thematic
   vowel that the plural drops, replacing it with the plural
   allomorph. Everything else is kept verbatim.

2. **Substitute the underspecified twin.** When the lemma's
   ``opportunity`` column says the plural takes ``-i`` or ``-e``,
   the paper's analysis claims the stem-final obstruent is
   underspecified. So we walk from the right edge, swapping each
   palatalisation candidate (``k, ɡ, t, d, s, z``) for its
   underspecified label (``K, G, T, D, S, Z``) until we hit a vowel.
   This handles cluster stems naturally: in ``musk-`` both /s/ and
   /k/ get substituted, in ``punkt-`` only the obstruents do (the
   ``/n/`` has no underspec twin so it is skipped, but the walk
   continues past it).

   When ``opportunity == "uri"`` we leave the obstruent prespecified
   — that is the spell-out rule's claim about which roots take which
   allomorphs, and the analysis predicts inalterability for these.

3. **Append the suffix.** The suffix segments are taken directly
   from the ``opportunity`` column: ``"i" → (i,)``, ``"e" → (e,)``,
   ``"uri" → (u, r, i)``. This sidesteps spell-out in the pipeline;
   we test palatalisation given the morphologically attested suffix.
   Spell-out's choice can be verified separately by comparing the
   stem-final's features to the P condition (see driver).
"""

from __future__ import annotations

from typing import Final, Mapping

from .diagnostics.compare import split_variants
from .inventory import FeatureInventory
from .segments import Segment, Word, tokenize_ipa

# Vowels and vocoids that break a consonantal cluster. /j/, /w/ and
# their diphthong-tie variants count as vocoids here so a preceding
# diphthong (e.g., -aj- before /k/) does not glue itself onto the
# trailing consonant cluster.
_NON_CONSONANTAL: Final[frozenset[str]] = frozenset(
    {
        "a",
        "e",
        "i",
        "o",
        "u",
        "ə",
        "ɨ",
        "j",
        "w",
        "o̯",
        "e̯",
        "u̯",
        "i̯",
        "a̯",
    }
)

# Surface-segment → underspec-label substitutions. Looked up by IPA
# token. Other consonants (n, m, r, l, ...) are left as-is — they
# have no underspec twin in the analysis.
_UNDERSPEC_FOR_SURFACE: Final[Mapping[str, str]] = {
    "k": "K",
    "ɡ": "G",
    "t": "T",
    "d": "D",
    "s": "S",
    # No /z/ → /Z/: the paper treats /z/ as fully specified and
    # inalterable (latex.tex:170).
}

# Dorsal vs coronal split: dorsals underspec even before -e (the -e
# plurals of dorsal stems palatalize categorically per the paper's
# productivity table). Coronals only underspec before -i in the
# single-consonant case; cluster cases underspec both members
# regardless of the suffix vowel.
_DORSAL_SURFACE: Final[frozenset[str]] = frozenset({"k", "ɡ"})
_CORONAL_SURFACE: Final[frozenset[str]] = frozenset({"t", "d", "s"})

# How each opportunity value spells out as a suffix segment list.
SUFFIX_TOKENS: Final[Mapping[str, tuple[str, ...]]] = {
    "i": ("i",),
    "e": ("e",),
    "uri": ("u", "r", "i"),
}

# Orthographic stem_final → canonical IPA segment. Used to correct
# the lemma IPA's final consonant when it disagrees with the
# orthographic classification (e.g. final-devoiced ``brat`` for
# lemma "brad"; ``ð`` in one variant of the same field). Without
# this correction we underspec the wrong segment and derive a wrong
# assibilated / palatalized outcome.
_ORTH_TO_IPA: Final[Mapping[str, str]] = {
    "c": "k",
    "g": "ɡ",
    "t": "t",
    "d": "d",
    "s": "s",
    "z": "z",
}


def _levenshtein(a: str, b: str) -> int:
    """Iterative Levenshtein distance."""
    if a == b:
        return 0
    if not a:
        return len(b)
    if not b:
        return len(a)
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        curr = [i] + [0] * len(b)
        for j, cb in enumerate(b, 1):
            cost = 0 if ca == cb else 1
            curr[j] = min(prev[j] + 1, curr[j - 1] + 1, prev[j - 1] + cost)
        prev = curr
    return prev[-1]


def pick_variant(
    ipa: str,
    stem_final: str | None = None,
    hint: str | None = None,
) -> str:
    """Pick the best variant from a pipe-separated IPA field.

    Wiktionary IPA fields often list several variants that differ on
    both real dialectal facts (schwa-reduction vs. citation-form
    vowels) and on transcription convention (whether internal
    orthographic ``c`` before front vowels is written as ``t͡ʃ`` or
    left as ``s``/``t͡s`` per a non-Romanian transcriber). The first
    variant is not always the most useful for our validator, which
    compares against G2P-derived plural IPA.

    Selection strategy (in order):

    1. Drop foreign-language annotations via
       :func:`~phonology.diagnostics.compare.split_variants`.
    2. If ``stem_final`` is given, restrict candidates to variants
       containing the orthographic stem-final's canonical IPA symbol.
       This eliminates truncated variants like ``əman`` (missing the
       final /t/ of "amant").
    3. If ``hint`` is given (typically ``to_ipa(lemma)``), pick the
       candidate with the smallest edit distance to the hint. This
       prefers variants that agree with our G2P convention on
       internal palatalization (``t͡ʃ`` for ``c`` before front vowels)
       and on segment inventory, which matches what our plural-IPA
       comparison expects.
    4. If no hint is given, pick the longest candidate (a proxy for
       "least abbreviated").

    Empty / all-foreign input returns the empty string.
    """
    parts = split_variants(ipa)
    if not parts:
        return ""
    candidates: tuple[str, ...] = parts
    if stem_final:
        target = _ORTH_TO_IPA.get(stem_final)
        if target is not None:
            filtered = tuple(
                p for p in parts if target in p or (target == "ɡ" and "g" in p)
            )
            if filtered:
                candidates = filtered
    if hint:
        return min(candidates, key=lambda p: _levenshtein(p, hint))
    return max(candidates, key=len)


def build_ur(
    lemma_ipa: str,
    opportunity: str,
    inventory: FeatureInventory,
    stem_final: str | None = None,
    lemma: str | None = None,
) -> Word | None:
    """Build the underlying-representation Word for one lexicon row.

    Returns ``None`` if the IPA is empty, contains a token absent from
    both the inventory and the underspec table, or the opportunity is
    not one of ``i``/``e``/``uri``. The caller should filter such rows
    out of the validation set rather than treat ``None`` as a failure.

    When ``lemma`` (the orthographic form) is supplied, ``pick_variant``
    receives ``to_ipa(lemma)`` as a similarity hint so it can prefer
    Wiktionary variants that agree with our G2P convention on
    internal-palatalization spellings — matching how the plural IPA
    (typically G2P-derived) will look for comparison.
    """
    if not lemma_ipa or opportunity not in SUFFIX_TOKENS:
        return None
    hint: str | None = None
    if lemma:
        # Late import to avoid a heavy load path at module import
        # time; ``romanian_processor_lib`` pulls in the full processor
        # stack.
        from romanian_processor_lib import to_ipa

        hint = to_ipa(lemma)
    tokens = list(tokenize_ipa(pick_variant(lemma_ipa, stem_final, hint=hint)))
    if not tokens:
        return None

    # 1. Strip thematic final vowel. ``ə`` is the feminine -ă class;
    #    ``e`` is the -e class. Other final vowels (``o``, ``ɨ``,
    #    word-final ``a`` in indefinite forms) we leave: they aren't
    #    the inflectional thematic vowel under the paper's analysis.
    if tokens and tokens[-1] in {"ə", "e"} and len(tokens) >= 2:
        tokens.pop()

    # 1b. Reconcile the stem-final consonant with the orthographic
    #     stem_final. Wiktionary occasionally records a final-devoiced
    #     IPA (e.g. ``brat`` for lemma "brad"), which would then get
    #     the wrong underspec twin. Only correct when there is a
    #     single trailing candidate: in cluster stems the last IPA
    #     consonant is not necessarily the one ``stem_final`` names
    #     (e.g. "muscă"'s stem_final is 's' but IPA-last is 'k').
    if stem_final and tokens:
        canonical = _ORTH_TO_IPA.get(stem_final)
        if canonical is not None:
            trailing: list[int] = []
            for i in range(len(tokens) - 1, -1, -1):
                if tokens[i] in _NON_CONSONANTAL:
                    break
                trailing.append(i)
            trailing_candidates = [
                i for i in trailing if tokens[i] in _UNDERSPEC_FOR_SURFACE
            ]
            if len(trailing_candidates) == 1:
                tokens[trailing_candidates[0]] = canonical

    # 2. Underspec substitution. The rule differs by segment class and
    #    by whether the trailing consonant run is a cluster.
    #
    #    * Cluster (two or more adjacent palatalisation candidates at
    #      the right edge — e.g. -sk-, -st-, -kt-): underspec BOTH
    #      candidates. The paper needs this for muskə → muʃte, prost
    #      → proʃti, punkt → punkte etc.
    #    * Non-cluster dorsal (single trailing /k/ or /ɡ/): underspec
    #      regardless of opp — dorsals palatalize categorically
    #      before both -i and -e (paper productivity table
    #      latex.tex:130-133).
    #    * Non-cluster coronal (single trailing /t, d, s, z/):
    #      underspec ONLY when opp == "i". Coronal palatalization
    #      before -i is productive; before -e it is categorically
    #      absent (latex.tex:135-141). Leaving coronals prespecified
    #      before -e makes them inalterable to dorsal-pal (whose
    #      broad terminator would otherwise let /e/ license it) and
    #      to assibilation (which no longer sees a postalveolar
    #      trigger).
    if opportunity in {"i", "e"}:
        cluster_span: list[int] = []
        for i in range(len(tokens) - 1, -1, -1):
            if tokens[i] in _NON_CONSONANTAL:
                break
            cluster_span.append(i)
        cluster_span.reverse()
        candidates = [
            (i, _UNDERSPEC_FOR_SURFACE[tokens[i]])
            for i in cluster_span
            if tokens[i] in _UNDERSPEC_FOR_SURFACE
        ]
        if len(candidates) >= 2:
            for idx, repl in candidates:
                tokens[idx] = repl
        elif len(candidates) == 1:
            idx, repl = candidates[0]
            surface = tokens[idx]
            if surface in _DORSAL_SURFACE:
                tokens[idx] = repl
            elif surface in _CORONAL_SURFACE and opportunity == "i":
                tokens[idx] = repl

    # 3. Materialise segments from inventory (base or underspec).
    segs: list[Segment] = []
    for tok in tokens:
        if tok in inventory.base or tok in inventory.underspec:
            segs.append(inventory.segment(tok))
        else:
            return None  # unknown token: caller treats as skip

    # 4. Append the plural suffix segments.
    for tok in SUFFIX_TOKENS[opportunity]:
        if tok not in inventory.base:
            return None
        segs.append(inventory.segment(tok))

    return tuple(segs)


__all__: Final[tuple[str, ...]] = (
    "SUFFIX_TOKENS",
    "build_ur",
    "pick_variant",
)
