#!/usr/bin/env python3
"""Tests for the Italian palatalization analysis.

Pins (``phonology.analyses.italian_palatalization`` and
``src/phonology/inventories/italian_features.json``):

- **Inventory.** The 27 Romanian features in order; every segment shared
  with Romanian copies its values except the declared deviations (e o
  Tense +); the new segments are the declared ones; geminates differ
  from their singleton in Long only; bundles are unique before and
  after the analysis patches; every label is a token of the Italian
  tokenizer.
- **Tokenizer.** Geminates written across the syllable dot become one
  long segment; ː attaches; ɛ ɔ stay distinct; stress and dots skip.
- **Golden derivations** (lexical policy): amico+i → amit͡ʃi, amica+e →
  amike (Height effect with an underspecified root), fuoco+i → fwɔki
  (prespecified), medico+i → mɛdit͡ʃi, lago+i → laɡi, asparago+i →
  asparad͡ʒi, belga+i → bɛld͡ʒi, occhio+i → no palatalization (NDE),
  plus the cluster and geminate outputs (osco+i → ɔʃːi, tombacco+i →
  tombat͡ʃːi) and the contrasts the other policies make.
- **Rows and the driver contract.** The attested-plural reconstruction,
  the row factory, ``iter_rows`` over a synthetic data tree written in
  the cell and lexeme schemas of ``lexpipe.langs.italian.features``
  (root evidence pooled over amico/amica, rows without IPA counted and
  left out), and ``scripts/validate_lang_rules.py --lang it`` on that
  tree.
- **Profiles.** build_profile() is the lexical variant; the variants
  differ only in their specification policy; hooks are wired.
- **Romanian untouched.** Importing the Italian analysis changes no
  Romanian default (the byte-level replay is
  ``tests/test_romanian_replay.py``).
- **Real data** (``needs_data``): the golden rows on the generated cell
  data, the lexical variant's only mismatches (doublets and the QC's
  e_palatal suspects), the Height effect under every policy, and the
  coronal and lexically palatal sanity rows.

The transcriptions below are English Wiktionary's (CC BY-SA 4.0, via
wiktextract; Ylonen 2022), as they appear in the generated data.
"""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
from collections import Counter
from collections.abc import Mapping, Sequence
from pathlib import Path
from types import ModuleType

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import phonology.analyses.italian_palatalization as it  # noqa: E402
import phonology.analyses.romanian_palatalization as romanian  # noqa: E402
from lexpipe.langs.italian import features  # noqa: E402
from phonology.diagnostics.categorise import RefinedCategory  # noqa: E402
from phonology.diagnostics.compare import ROMANIAN_NORMALISERS  # noqa: E402
from phonology.g2p import build_ur_from_row  # noqa: E402
from phonology.inventory import (  # noqa: E402
    INVENTORY_DIR,
    FeatureInventory,
)
from phonology.lexicon import LexRow  # noqa: E402
from phonology.pipeline import RulePipeline  # noqa: E402
from phonology.profile import AnalysisProfile  # noqa: E402
from phonology.segments import (  # noqa: E402
    ROMANIAN_TOKENIZER,
    LengthHandling,
    segments_to_ipa,
    tokenize_ipa,
)
from phonology.validation import (  # noqa: E402
    RowResult,
    evaluate_row,
    run_validation,
)

pytestmark = pytest.mark.it

ITALIAN_JSON = INVENTORY_DIR / "italian_features.json"
ROMANIAN_JSON = PROJECT_ROOT / "romanian_features.json"
CELLS_CSV = PROJECT_ROOT / "data" / "italian" / "italian_cells.csv"
DRIVER = PROJECT_ROOT / "scripts" / "validate_lang_rules.py"

SHARED = "i u a p b t d k ɡ t͡s t͡ʃ d͡ʒ f v s z ʃ ʒ m n l r j w".split()
DEVIATIONS = {"e": {"Tense": "+"}, "o": {"Tense": "+"}}
NEW = {
    "ɛ": ("e", {"Tense": "-"}),
    "ɔ": ("o", {"Tense": "-"}),
    "d͡z": ("t͡s", {"Voice": "+"}),
    "ɲ": ("n", {"Anterior": "-", "Distributed": "+"}),
    "ʎ": ("l", {"Anterior": "-", "Distributed": "+"}),
}
GEMINABLE = "p b t d k ɡ t͡s d͡z t͡ʃ d͡ʒ f v s ʃ m n ɲ l ʎ r".split()


def _raw(path: Path) -> dict[str, object]:
    data: dict[str, object] = json.loads(path.read_text(encoding="utf-8"))
    return data


def _segments(path: Path) -> dict[str, dict[str, str]]:
    segments = _raw(path)["segments"]
    assert isinstance(segments, dict)
    return segments


def _row(
    word: str,
    form: str,
    ipa: str,
    opportunity: str,
    *,
    stem_final: str,
    mutation: bool,
    target_class: str = it.DORSAL,
    cluster: str = "",
    alternates: bool = False,
    post_tonic: bool | None = None,
    variable: bool = False,
) -> it.ItalianRow:
    """An ItalianRow through the row factory, as iter_rows builds it."""
    nde = target_class[4:] if target_class.startswith("nde:") else ""
    target = it.CORONAL if target_class == it.CORONAL else it.DORSAL
    if target_class == it.FRONTSTEM:
        target = ""
    return it.row_factory(
        {
            "word": word,
            "form": form,
            "pos": "noun",
            "lemma_ipa": ipa,
            "opportunity": opportunity,
            "stem_final": stem_final,
            "mutation": str(mutation),
            "cluster": cluster,
            "target": target,
            "nde_class": nde,
            "frontstem": str(target_class == it.FRONTSTEM),
            "root_alternates": str(alternates),
            "root_variable": str(variable),
            "velar_post_tonic": "" if post_tonic is None else str(post_tonic),
        }
    )


AMICO = _row(
    "amico", "amici", "aˈmi.ko", "i", stem_final="c", mutation=True,
    alternates=True, post_tonic=True,
)  # fmt: skip
AMICA = _row(
    "amica", "amiche", "aˈmi.ka", "e", stem_final="c", mutation=False,
    alternates=True, post_tonic=True,
)  # fmt: skip
FUOCO = _row(
    "fuoco", "fuochi", "ˈfwɔ.ko", "i", stem_final="c", mutation=False,
    post_tonic=True,
)  # fmt: skip
MEDICO = _row(
    "medico", "medici", "ˈmɛ.di.ko", "i", stem_final="c", mutation=True,
    alternates=True, post_tonic=False,
)  # fmt: skip
LAGO = _row(
    "lago", "laghi", "ˈla.ɡo", "i", stem_final="g", mutation=False,
    post_tonic=True,
)  # fmt: skip
ASPARAGO = _row(
    "asparago", "asparagi", "asˈpa.ra.ɡo", "i", stem_final="g",
    mutation=True, alternates=True, post_tonic=False,
)  # fmt: skip
BELGA = _row(
    "belga", "belgi", "ˈbɛl.ɡa", "i", stem_final="g", mutation=True,
    cluster="lg", alternates=True, post_tonic=True,
)  # fmt: skip
OCCHIO = _row(
    "occhio", "occhi", "ˈɔk.kjo", "i", stem_final="c", mutation=False,
    target_class=it.NDE_OCCHIO, cluster="chi", post_tonic=True,
)  # fmt: skip
OSCO = _row(
    "osco", "osci", "ˈɔs.ko", "i", stem_final="c", mutation=True,
    cluster="sc", alternates=True, post_tonic=True,
)  # fmt: skip
TOMBACCO = _row(
    "tombacco", "tombacci", "tomˈbak.ko", "i", stem_final="c",
    mutation=True, cluster="cc", alternates=True, post_tonic=True,
)  # fmt: skip
BOSCO = _row(
    "bosco", "boschi", "ˈbɔs.ko", "i", stem_final="c", mutation=False,
    cluster="sc", post_tonic=True,
)  # fmt: skip
PSICHE = _row(
    "psiche", "psichi", "ˈpsi.ke", "i", stem_final="c", mutation=False,
    target_class=it.NDE_PSICHE, cluster="che", post_tonic=True,
)  # fmt: skip
GATTO = _row(
    "gatto", "gatti", "ˈɡat.to", "i", stem_final="t", mutation=False,
    target_class=it.CORONAL, cluster="tt",
)  # fmt: skip
LUCE = _row(
    "luce", "luci", "ˈlu.t͡ʃe", "i", stem_final="c", mutation=False,
    target_class=it.FRONTSTEM,
)  # fmt: skip


@pytest.fixture(scope="module")
def lexical() -> AnalysisProfile:
    return it.build_profile()


@pytest.fixture(scope="module")
def variants() -> dict[str, AnalysisProfile]:
    return {name: factory() for name, factory in it.PROFILE_VARIANTS.items()}


def _fired(result: RowResult) -> list[str]:
    assert result.derivation is not None
    return [step.rule for step in result.derivation.steps if step.fired]


# ---------------------------------------------------------------------------
# Inventory
# ---------------------------------------------------------------------------


class TestInventory:
    def test_features_are_the_romanian_list(self) -> None:
        assert (
            _raw(ITALIAN_JSON)["features"] == _raw(ROMANIAN_JSON)["features"]
        )
        assert list(_raw(ITALIAN_JSON)) == ["metadata", "features", "segments"]

    def test_shared_segments_copy_romanian_exactly(self) -> None:
        italian, ro = _segments(ITALIAN_JSON), _segments(ROMANIAN_JSON)
        for label in SHARED:
            assert italian[label] == ro[label], label
        for label, updates in DEVIATIONS.items():
            changed = {
                f for f, v in italian[label].items() if ro[label][f] != v
            }
            assert changed == set(updates), label
            for feature, value in updates.items():
                assert italian[label][feature] == value

    def test_romanian_only_segments_are_absent(self) -> None:
        italian = _segments(ITALIAN_JSON)
        assert not {"ɨ", "ə", "h"} & set(italian)

    def test_new_segments_derive_from_their_bases(self) -> None:
        italian = _segments(ITALIAN_JSON)
        for label, (base, updates) in NEW.items():
            expected = dict(italian[base])
            expected.update(updates)
            assert italian[label] == expected, label

    def test_geminates_differ_in_length_only(self) -> None:
        italian = _segments(ITALIAN_JSON)
        for label in GEMINABLE:
            long = italian[label + "ː"]
            assert long["Long"] == "+" and italian[label]["Long"] == "-"
            assert {f: v for f, v in long.items() if f != "Long"} == {
                f: v for f, v in italian[label].items() if f != "Long"
            }
        for label in ("z", "j", "w", "ʒ"):
            assert label + "ː" not in italian
        declared = set(SHARED) | set(DEVIATIONS) | set(NEW)
        declared |= {label + "ː" for label in GEMINABLE}
        assert set(italian) == declared

    def test_mid_vowels_carry_the_only_tense_values(self) -> None:
        italian = _segments(ITALIAN_JSON)
        tense = {
            s: v["Tense"] for s, v in italian.items() if v["Tense"] != "0"
        }
        assert tense == {"e": "+", "o": "+", "ɛ": "-", "ɔ": "-"}
        assert {v["Stress"] for v in italian.values()} == {"-"}

    @pytest.mark.parametrize("patched", [False, True])
    def test_bundles_are_unique(self, patched: bool) -> None:
        inventory = (
            it.load_inventory()
            if patched
            else FeatureInventory.load(ITALIAN_JSON)
        )
        bundles = Counter(
            tuple(sorted(feats.items())) for feats in inventory.base.values()
        )
        assert [b for b, n in bundles.items() if n > 1] == []
        assert len(inventory.base) == 51

    def test_every_label_is_a_token(self) -> None:
        for label in _segments(ITALIAN_JSON):
            assert it.tokens(label) == (label,), label

    def test_underspecified_segments_are_dorsals_only(self) -> None:
        inventory = it.load_inventory()
        for u in it.UNDERSPEC:
            assert inventory.base[u.base]["DORSAL"] == "+", u.label
            assert inventory.base[u.base]["CORONAL"] == "-", u.label
            assert "CORONAL" in u.clear
            assert u.label not in inventory.base
        assert {u.label for u in it.UNDERSPEC} == {"K", "G", "Kː", "Gː", "Kˢ"}

    def test_patches_name_existing_segments(self) -> None:
        italian = _segments(ITALIAN_JSON)
        for patch in it.PATCHES:
            assert patch.segment in italian and patch.reason
        inventory = it.load_inventory()
        for label in ("i", "j"):
            assert inventory.base[label]["CORONAL"] == "+"
        assert inventory.base["e"]["CORONAL"] == "-"  # never a trigger
        for label in ("t͡ʃ", "d͡ʒ", "t͡ʃː", "d͡ʒː", "ʃː"):
            assert inventory.base[label]["DORSAL"] == "+"


# ---------------------------------------------------------------------------
# Tokenizer and geminates
# ---------------------------------------------------------------------------


class TestTokenizer:
    @pytest.mark.parametrize(
        ("ipa", "merged"),
        [
            ("ˈsak.ko", "ˈsakːo"),
            ("ˈfat.t͡ʃa", "ˈfat͡ʃːa"),
            ("ˈpit.t͡sa", "ˈpit͡sːa"),
            ("ˈmɛd.d͡zo", "ˈmɛd͡zːo"),
            ("ˈɔd.d͡ʒi", "ˈɔd͡ʒːi"),
            ("viʎˈʎak.ko", "viˈʎːakːo"),
            ("ˈpeʃ.ʃe", "ˈpeʃːe"),
            ("kavaˈtappi", "kavaˈtapːi"),
            ("ˈt͡si.i", "ˈt͡si.i"),  # vowels never merge
            ("aˈmi.ko", "aˈmi.ko"),
            ("ˈsakːo", "ˈsakːo"),  # idempotent
        ],
    )
    def test_merge_geminates(self, ipa: str, merged: str) -> None:
        assert it.merge_geminates(ipa) == merged
        assert it.merge_geminates(merged) == merged

    @pytest.mark.parametrize(
        ("ipa", "expected"),
        [
            ("ˈsak.ko", ("s", "a", "kː", "o")),
            ("ˈfat.t͡ʃa", ("f", "a", "t͡ʃː", "a")),
            ("ˈfwɔ.ko", ("f", "w", "ɔ", "k", "o")),
            ("ˈbɛl.ɡa", ("b", "ɛ", "l", "ɡ", "a")),
            ("ˈɔk.kjo", ("ɔ", "kː", "j", "o")),
            ("aɲˈɲɛl.lo", ("a", "ɲː", "ɛ", "lː", "o")),
            ("ˈd͡zɛ.ro", ("d͡z", "ɛ", "r", "o")),
        ],
    )
    def test_tokens(self, ipa: str, expected: tuple[str, ...]) -> None:
        assert it.tokens(ipa) == expected

    def test_config(self) -> None:
        config = it.ITALIAN_TOKENIZER
        assert config.length is LengthHandling.ATTACH
        assert "ɛ" not in config.fold and "ɔ" not in config.fold
        # ɛ/ɔ fold to e/o only under the Romanian convention.
        assert tokenize_ipa("ˈpɛ.ʃːe", tokenizer=config) == (
            "p", "ɛ", "ʃː", "e",
        )  # fmt: skip
        assert tokenize_ipa("ˈpɛ.ʃːe") == ("p", "e", "ʃ", "e")

    def test_compare_strip_keeps_length(self) -> None:
        strip = it.ITALIAN_NORMALISERS.strip
        assert strip("ˈsak.ki") == "sakːi" == strip("sakːi")
        assert strip("ˈsa.ki") == "saki" != strip("ˈsak.ki")


# ---------------------------------------------------------------------------
# Golden derivations
# ---------------------------------------------------------------------------

GOLDEN = [
    (AMICO, "amiKi", "amit͡ʃi", ["dorsal-pal"]),
    (AMICA, "amiKe", "amike", ["dorsal-default"]),
    (FUOCO, "fwɔki", "fwɔki", []),
    (MEDICO, "mɛdiKi", "mɛdit͡ʃi", ["dorsal-pal"]),
    (LAGO, "laɡi", "laɡi", []),
    (ASPARAGO, "asparaGi", "asparad͡ʒi", ["dorsal-pal"]),
    (BELGA, "bɛlGi", "bɛld͡ʒi", ["dorsal-pal"]),
    (OCCHIO, "ɔkːji", "ɔkːi", ["glide-absorption"]),
    (OSCO, "ɔsKˢi", "ɔʃːi", ["dorsal-pal", "sc-coalescence", "s-deletion"]),
    (TOMBACCO, "tombaKːi", "tombat͡ʃːi", ["dorsal-pal"]),
    (BOSCO, "bɔski", "bɔski", []),
    (PSICHE, "psiki", "psiki", []),
    (GATTO, "ɡatːi", "ɡatːi", []),
    (LUCE, "lut͡ʃi", "lut͡ʃi", []),
]


class TestGoldenDerivations:
    @pytest.mark.parametrize(
        ("row", "ur", "sr", "fired"),
        GOLDEN,
        ids=[row.lemma for row, *_ in GOLDEN],
    )
    def test_lexical(
        self,
        lexical: AnalysisProfile,
        row: it.ItalianRow,
        ur: str,
        sr: str,
        fired: list[str],
    ) -> None:
        result = evaluate_row(row, lexical)
        assert result.ur is not None
        assert "".join(s.label for s in result.ur) == ur
        assert result.predicted_sr == sr
        assert _fired(result) == fired
        assert result.boolean_predicted is row.mutation
        assert result.boolean_match and result.ipa_match

    def test_occhio_has_no_palatalization(
        self, lexical: AnalysisProfile
    ) -> None:
        result = evaluate_row(OCCHIO, lexical)
        assert not set(_fired(result)) & it.PALATALIZATION_RULE_NAMES
        assert result.boolean_predicted is False

    def test_amica_is_underspecified_but_not_palatalized(
        self, lexical: AnalysisProfile
    ) -> None:
        # The Height effect: /amiK/ before -e reaches the default.
        assert it.is_underspecified(AMICA, it.LEXICAL)
        result = evaluate_row(AMICA, lexical)
        assert result.derivation is not None
        step = result.derivation.steps[0]
        assert step.rule == "dorsal-pal" and not step.fired

    def test_prespecified_k_resists_before_i(self) -> None:
        inventory = it.load_inventory()
        pipeline = RulePipeline(it.RULES, resolve=inventory.segment)
        word = tuple(inventory.segment(x) for x in ("f", "w", "ɔ", "k", "i"))
        derivation = pipeline.derive(word)
        assert segments_to_ipa(derivation.sr, inventory.base) == "fwɔki"
        assert not any(step.fired for step in derivation.steps)

    def test_e_never_triggers(self) -> None:
        inventory = it.load_inventory()
        pipeline = RulePipeline(it.RULES, resolve=inventory.segment)
        for vowel, out in (("i", "t͡ʃ"), ("j", "t͡ʃ"), ("e", "k"), ("a", "k")):
            word = (inventory.segment("K"), inventory.segment(vowel))
            sr = segments_to_ipa(pipeline.derive(word).sr, inventory.base)
            assert sr == out + vowel, vowel

    @pytest.mark.parametrize(
        ("policy", "row", "sr"),
        [
            (it.ALL_UNDERSPECIFIED, FUOCO, "fwɔt͡ʃi"),
            (it.ALL_UNDERSPECIFIED, LAGO, "lad͡ʒi"),
            (it.ALL_UNDERSPECIFIED, OCCHIO, "ɔt͡ʃːi"),
            (it.ALL_UNDERSPECIFIED, PSICHE, "psit͡ʃi"),
            (it.ALL_UNDERSPECIFIED, BOSCO, "bɔʃːi"),
            (it.ALL_UNDERSPECIFIED, AMICA, "amike"),
            (it.STRESS, AMICO, "amiki"),
            (it.STRESS, BELGA, "bɛlɡi"),
            (it.STRESS, MEDICO, "mɛdit͡ʃi"),
            (it.STRESS, FUOCO, "fwɔki"),
            (it.STRESS, OCCHIO, "ɔkːi"),
            (it.STRESS, GATTO, "ɡatːi"),
        ],
    )
    def test_other_policies(
        self,
        variants: dict[str, AnalysisProfile],
        policy: str,
        row: it.ItalianRow,
        sr: str,
    ) -> None:
        assert evaluate_row(row, variants[policy]).predicted_sr == sr

    def test_no_policy_underspecifies_coronals_or_palatals(self) -> None:
        for policy in it.POLICIES:
            assert not it.is_underspecified(GATTO, policy)
            assert not it.is_underspecified(LUCE, policy)
            assert it.is_underspecified(FUOCO, policy) is (
                policy == it.ALL_UNDERSPECIFIED
            )
        with pytest.raises(ValueError):
            it.is_underspecified(FUOCO, "no_such_policy")

    def test_fallback_applies_the_policy_not_the_data(self) -> None:
        stress = it.FallbackPredictor(it.STRESS)
        assert stress(AMICO) is False  # post-tonic, though amici palatalizes
        assert stress(MEDICO) is True
        assert stress(AMICA) is False  # -e never triggers
        assert it.FallbackPredictor(it.LEXICAL)(AMICO) is True

    def test_expected_firings_is_the_attested_outcome(self) -> None:
        assert it.expected_firings(AMICO) == {"dorsal-pal"}
        assert it.expected_firings(FUOCO) == frozenset()
        assert it.expected_firings(OSCO) == {
            "dorsal-pal", "sc-coalescence", "s-deletion",
        }  # fmt: skip


# ---------------------------------------------------------------------------
# Attested plural and rows
# ---------------------------------------------------------------------------


class TestRows:
    @pytest.mark.parametrize(
        ("row", "plural"),
        [
            (AMICO, "aˈmi.t͡ʃi"),
            (AMICA, "aˈmi.ke"),
            (FUOCO, "ˈfwɔ.ki"),
            (BELGA, "ˈbɛl.d͡ʒi"),
            (OCCHIO, "ˈɔkːi"),
            (OSCO, "ˈɔʃːi"),
            (TOMBACCO, "tomˈbat͡ʃːi"),
            (BOSCO, "ˈbɔs.ki"),
            (PSICHE, "ˈpsi.ki"),
            (GATTO, "ˈɡatːi"),
            (LUCE, "ˈlu.t͡ʃi"),
        ],
        ids=lambda x: x.lemma if isinstance(x, LexRow) else x,
    )
    def test_plural_ipa(self, row: it.ItalianRow, plural: str) -> None:
        assert row.ipa_pl == plural

    def test_plural_ipa_refuses_a_transcription_that_misfits(self) -> None:
        def plural(
            ipa: str,
            *,
            stem_final: str = "c",
            target_class: str = it.DORSAL,
            target_letters: str = "c",
        ) -> str:
            return it.plural_ipa(
                ipa,
                "i",
                stem_final=stem_final,
                target_class=target_class,
                target_letters=target_letters,
                mutation=True,
            )

        assert plural("aˈmi.ko") == "aˈmi.t͡ʃi"
        assert plural("seˈvi.t͡ʃe", target_class=it.NDE_PSICHE) == ""
        assert plural("aˈmi.ko", stem_final="g") == ""
        assert plural("aˈmi.ko", target_letters="sc") == ""
        assert plural("ˈɔk.kjo") == ""  # a root glide outside occhio
        assert plural("kafˈfɛ") == ""  # no final a/e/o
        occhio_e = it.plural_ipa(
            "ˈmak.kja",
            "e",
            stem_final="c",
            target_class=it.NDE_OCCHIO,
            target_letters="c",
            mutation=False,
        )
        assert occhio_e == "ˈmakːje"  # -e keeps the root glide

    def test_row_factory_fields(self) -> None:
        assert isinstance(AMICO, LexRow)
        assert AMICO.pos == "N" and AMICO.target_letters == "c"
        assert OSCO.target_letters == "sc" and TOMBACCO.target_letters == "cc"
        assert BELGA.target_letters == "g" and BELGA.cluster_letters == "lg"
        assert OCCHIO.target_class == it.NDE_OCCHIO
        assert GATTO.velar_post_tonic is None  # not a dorsal
        plain = it.row_factory({"word": "x"})
        assert plain.ipa_pl == "" and plain.target_class == ""
        assert it.ur_labels(plain, it.LEXICAL) is None

    def test_plain_lexrow_gets_no_ur(self, lexical: AnalysisProfile) -> None:
        row = LexRow(
            lemma="amico",
            plural="amici",
            pos="N",
            stem_final="c",
            opportunity="i",
            mutation=True,
            ipa_lemma="aˈmi.ko",
            ipa_pl="aˈmi.t͡ʃi",
            exception_reason="",
            cluster=AMICO.cluster,
        )
        assert lexical.ur_builder(row, lexical.inventory) is None
        assert lexical.fallback_predict(row) is False
        assert lexical.cluster_key(row) is None

    def test_classes_and_keys(self) -> None:
        assert it.row_class(FUOCO) == "root, velar post-tonic"
        assert it.row_class(MEDICO) == "root, velar not post-tonic"
        assert it.row_class(OCCHIO) == it.NDE_OCCHIO
        assert it.cluster_key(AMICA) == "+e root, velar post-tonic"


# ---------------------------------------------------------------------------
# iter_rows and the driver on a synthetic data tree
# ---------------------------------------------------------------------------

# (lexeme_id, root_id, word, pos, ipa_raw, stem_final, cluster, target,
#  nde_class, frontstem, velar_post_tonic)
LEXEMES = [
    ("it:amico:noun", "it:amic", "amico", "noun", "/aˈmi.ko/", "c", "",
     "dorsal", "", "False", "True"),
    ("it:amica:noun", "it:amic", "amica", "noun", "/aˈmi.ka/", "c", "",
     "dorsal", "", "False", "True"),
    ("it:fuoco:noun", "it:fuoc", "fuoco", "noun", "/ˈfwɔ.ko/", "c", "",
     "dorsal", "", "False", "True"),
    ("it:medico:noun", "it:medic", "medico", "noun", "/ˈmɛ.di.ko/", "c",
     "", "dorsal", "", "False", "False"),
    ("it:occhio:noun", "it:occhi", "occhio", "noun", "/ˈɔk.kjo/", "c",
     "chi", "dorsal", "occhio", "False", "True"),
    ("it:chirurgo:noun", "it:chirurg", "chirurgo", "noun", "/kiˈrur.ɡo/",
     "g", "rg", "dorsal", "", "False", "True"),
    ("it:lago:noun", "it:lag", "lago", "noun", "", "g", "", "dorsal", "",
     "False", ""),
    ("it:gatto:noun", "it:gatt", "gatto", "noun", "/ˈɡat.to/", "t", "tt",
     "coronal", "", "False", ""),
    ("it:luce:noun", "it:luc", "luce", "noun", "/ˈlu.t͡ʃe/", "c", "", "",
     "", "True", ""),
]  # fmt: skip

# (lexeme_id, cell, form, item, opportunity, mutation, filtered)
CELLS = [
    ("it:amico:noun", "pl", "amici", "0", "i", "True", "False"),
    ("it:amica:noun", "pl", "amiche", "0", "e", "False", "False"),
    ("it:fuoco:noun", "pl", "fuochi", "0", "i", "False", "False"),
    ("it:fuoco:noun", "pl", "fòcora", "1", "a", "", "True"),
    ("it:medico:noun", "pl", "medici", "0", "i", "True", "False"),
    ("it:occhio:noun", "pl", "occhi", "0", "i", "False", "False"),
    ("it:chirurgo:noun", "pl", "chirurghi", "0", "i", "False", "False"),
    ("it:chirurgo:noun", "pl", "chirurgi", "1", "i", "True", "False"),
    ("it:lago:noun", "pl", "laghi", "0", "i", "False", "False"),
    ("it:gatto:noun", "pl", "gatti", "0", "i", "False", "False"),
    ("it:luce:noun", "pl", "luci", "0", "i", "False", "False"),
]


def _write_csv(
    path: Path, fields: tuple[str, ...], rows: Sequence[Mapping[str, str]]
) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(fields), restval="")
        writer.writeheader()
        writer.writerows(rows)


@pytest.fixture(scope="module")
def tree(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """A repository root with synthetic Italian cell and lexeme tables in
    the data layer's schemas (only the columns the analysis reads)."""
    root = tmp_path_factory.mktemp("italian_tree")
    data = root / "data" / "italian"
    data.mkdir(parents=True)
    by_id = {}
    lexeme_rows = []
    for (lid, rid, word, pos, ipa, sf, cluster, target, nde, front,
         post) in LEXEMES:  # fmt: skip
        common = {
            "lexeme_id": lid,
            "root_id": rid,
            "word": word,
            "pos": pos,
            "stem_final": sf,
            "cluster": cluster,
            "target": target,
            "nde_class": nde,
            "frontstem": front,
            "velar_post_tonic": post,
            "in_domain": "True",
        }
        by_id[lid] = common
        lexeme_rows.append(
            {
                **common,
                "ipa_raw": ipa,
                "dorsal_target": str(target == "dorsal" and not nde),
            }
        )
    cell_rows = []
    for lid, cell, form, item, opp, mutation, filtered in CELLS:
        common = by_id[lid]
        cell_rows.append(
            {
                **{k: common[k] for k in common if k != "in_domain"},
                "cell": cell,
                "form": form,
                "item": item,
                "opportunity": opp,
                "mutation": mutation,
                "filtered": filtered,
                "aligned": "True",
                "adjudicated_out": "False",
                "attested_dictionary": "True",
                "in_domain": "True",
            }
        )
    _write_csv(
        data / "italian_lexicon_complete.csv",
        features.LEXEME_FIELDS,
        lexeme_rows,
    )
    _write_csv(data / "italian_cells.csv", features.CELL_FIELDS, cell_rows)
    return root


class TestIterRows:
    def test_rows_and_coverage(self, tree: Path) -> None:
        rows, coverage = it.load_rows(tree)
        words = [(r.lemma, r.plural) for r in rows]
        assert words == [  # the cell table's order
            ("amico", "amici"),
            ("amica", "amiche"),
            ("fuoco", "fuochi"),
            ("medico", "medici"),
            ("occhio", "occhi"),
            ("chirurgo", "chirurghi"),
            ("chirurgo", "chirurgi"),
        ]
        # lago has no transcription; fòcora is a filtered -a plural.
        assert coverage.in_scope == 8 and coverage.no_ipa == 1
        assert coverage.by_class == {
            "dorsal": (7, 6),
            "nde:occhio": (1, 1),
        }
        assert list(it.iter_rows(tree)) == list(rows)

    def test_root_evidence_is_pooled(self, tree: Path) -> None:
        loaded, _ = it.load_rows(tree)
        rows = {(r.lemma, r.plural): r for r in loaded}
        assert rows[("amica", "amiche")].root_alternates  # amici (R4)
        assert not rows[("fuoco", "fuochi")].root_alternates
        chirurgo = rows[("chirurgo", "chirurghi")]
        assert chirurgo.root_alternates and chirurgo.root_variable

    def test_validation_on_the_tree(self, tree: Path) -> None:
        rows, coverage = it.load_rows(tree)
        reports = {
            name: run_validation(factory(), rows)
            for name, factory in it.PROFILE_VARIANTS.items()
        }
        lexical = reports[it.LEXICAL]
        missed = [r.row.plural for r in lexical.results if not r.boolean_match]
        assert missed == ["chirurghi"]  # the doublet's velar variant
        assert all(r.ipa_match == r.boolean_match for r in lexical.results)
        stress = reports[it.STRESS]
        missed = sorted(
            r.row.plural for r in stress.results if not r.boolean_match
        )
        assert missed == ["amici", "chirurgi"]
        text = it.render_summary(reports, coverage)
        assert "doublet root" in text and "amico type" in text
        assert "Height effect" in text

    def test_sanity_classes(self, tree: Path) -> None:
        rows, _ = it.load_rows(tree, classes=(it.CORONAL, it.FRONTSTEM))
        assert [r.lemma for r in rows] == ["gatto", "luce"]
        report = run_validation(it.build_profile(), rows)
        assert report.boolean_matches == report.ipa_matches == 2

    def test_main_writes_the_summary(
        self, tree: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        out = tmp_path / "summary.txt"
        assert it.main(["--root", str(tree), "--out", str(out)]) == 0
        text = out.read_text(encoding="utf-8")
        assert text.startswith("Italian velar palatalization")
        assert "Sanity rows" in text
        assert it.main(["--root", str(tree), "--trace", "amico"]) == 0
        assert "dorsal-pal" in capsys.readouterr().out

    def test_driver_contract(self, tree: Path) -> None:
        spec = importlib.util.spec_from_file_location("driver", DRIVER)
        assert spec is not None and spec.loader is not None
        driver: ModuleType = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(driver)
        assert driver.main(["--lang", "it", "--root", str(tree)]) == 0
        assert (
            driver.main(
                ["--lang", "it", "--root", str(tree), "--variant", "stress"]
            )
            == 0
        )
        out = tree / "analysis" / "italian"
        report = (out / "italian_rule_validation.txt").read_text("utf-8")
        assert report.startswith("Analysis: italian-palatalization (lexical)")
        assert (out / "italian_rule_mismatches_stress.csv").is_file()
        with (out / "italian_rule_mismatches.csv").open(
            encoding="utf-8"
        ) as fh:
            mismatches = list(csv.DictReader(fh))
        assert [m["plural"] for m in mismatches] == ["chirurghi"]
        assert mismatches[0]["category"] == RefinedCategory.RULE_FAILURE


# ---------------------------------------------------------------------------
# Profiles
# ---------------------------------------------------------------------------


class TestProfiles:
    def test_build_profile_is_the_lexical_variant(
        self, lexical: AnalysisProfile
    ) -> None:
        assert set(it.PROFILE_VARIANTS) == {
            "lexical", "all_underspecified", "stress",
        }  # fmt: skip
        assert it.PROFILE_VARIANTS["lexical"]() == lexical
        assert it.build_profile() == lexical
        assert lexical.name == "italian-palatalization (lexical)"

    def test_variants_differ_in_policy_only(
        self, variants: dict[str, AnalysisProfile]
    ) -> None:
        base = variants[it.LEXICAL]
        for name, profile in variants.items():
            assert profile.rules == base.rules == it.RULES
            assert profile.inventory == base.inventory
            assert profile.ur_builder == it.UrBuilder(name)
            assert profile.fallback_predict == it.FallbackPredictor(name)
        with pytest.raises(ValueError):
            it.build_variant("no_such_policy")

    def test_hooks(self, lexical: AnalysisProfile) -> None:
        assert lexical.tokenizer is it.ITALIAN_TOKENIZER
        assert lexical.normalisers is it.ITALIAN_NORMALISERS
        assert lexical.row_factory is it.row_factory
        assert lexical.categoriser is it.categorise
        assert lexical.cluster_key is it.cluster_key
        assert lexical.inventory_json == ITALIAN_JSON
        assert lexical.palatalization_rule_names == {
            "dorsal-pal", "sc-coalescence", "s-deletion",
        }  # fmt: skip
        assert lexical.target_stems == {"c", "g", "t", "d", "s", "z"}

    def test_categoriser_uses_the_shared_vocabulary(
        self, lexical: AnalysisProfile
    ) -> None:
        cat = it.categorise(
            FUOCO,
            "fwɔt͡ʃi",
            FUOCO.ipa_pl,
            True,
            ur_built=True,
            normalisers=it.ITALIAN_NORMALISERS,
        )
        assert cat.category is RefinedCategory.RULE_FAILURE
        assert cat.edit_distance > 0
        failed = it.categorise(
            FUOCO,
            "",
            "",
            False,
            ur_built=False,
            normalisers=it.ITALIAN_NORMALISERS,
        )
        assert failed.category is RefinedCategory.UR_BUILD_FAILED


class TestRomanianUntouched:
    def test_romanian_defaults_survive_the_italian_import(self) -> None:
        profile = romanian.build_profile()
        assert profile.tokenizer is ROMANIAN_TOKENIZER
        assert profile.normalisers is ROMANIAN_NORMALISERS
        assert profile.ur_builder is build_ur_from_row
        assert profile.inventory_json == ROMANIAN_JSON
        assert tokenize_ipa("ˈkaːs") == ("k", "a", "s")
        assert tokenize_ipa("bɔk") == ("b", "o", "k")


# ---------------------------------------------------------------------------
# Real data
# ---------------------------------------------------------------------------

needs_data = pytest.mark.skipif(
    not CELLS_CSV.is_file(), reason=f"gitignored data absent: {CELLS_CSV}"
)


@pytest.fixture(scope="module")
def real() -> tuple[tuple[it.ItalianRow, ...], dict[str, list[RowResult]]]:
    if not CELLS_CSV.is_file():
        pytest.skip(f"gitignored data absent: {CELLS_CSV}")
    rows, _ = it.load_rows(PROJECT_ROOT)
    results = {
        name: list(run_validation(factory(), rows).results)
        for name, factory in it.PROFILE_VARIANTS.items()
    }
    return rows, results


def _find(
    results: list[RowResult], lemma: str, plural: str
) -> list[RowResult]:
    return [
        r for r in results if r.row.lemma == lemma and r.row.plural == plural
    ]


@pytest.mark.needs_data
@needs_data
class TestRealData:
    @pytest.mark.parametrize(
        ("lemma", "plural", "sr", "palatalized"),
        [
            ("amico", "amici", "amit͡ʃi", True),
            ("amica", "amiche", "amike", False),
            ("fuoco", "fuochi", "fwɔki", False),
            ("medico", "medici", "mɛdit͡ʃi", True),
            ("lago", "laghi", "laɡi", False),
            ("asparago", "asparagi", "asparad͡ʒi", True),
            ("belga", "belgi", "bɛld͡ʒi", True),
            ("occhio", "occhi", "ɔkːi", False),
        ],
    )
    def test_golden_rows(
        self,
        real: tuple[tuple[it.ItalianRow, ...], dict[str, list[RowResult]]],
        lemma: str,
        plural: str,
        sr: str,
        palatalized: bool,
    ) -> None:
        found = _find(real[1][it.LEXICAL], lemma, plural)
        assert found, (lemma, plural)
        for result in found:
            assert result.predicted_sr == sr
            assert result.boolean_predicted is palatalized
            assert result.boolean_match and result.ipa_match

    def test_lexical_mismatches_are_doublets_and_qc_suspects(
        self,
        real: tuple[tuple[it.ItalianRow, ...], dict[str, list[RowResult]]],
    ) -> None:
        kinds = Counter(it.mismatch_class(r) for r in real[1][it.LEXICAL])
        kinds.pop(None)
        assert set(kinds) == {
            "doublet root (palatal and velar -i plurals)",
            "-e palatalized in the data (QC suspect e_palatal)",
        }
        total = len(real[1][it.LEXICAL])
        assert sum(kinds.values()) < 0.01 * total

    def test_height_effect_under_every_policy(
        self,
        real: tuple[tuple[it.ItalianRow, ...], dict[str, list[RowResult]]],
    ) -> None:
        for name, results in real[1].items():
            e_dorsal = [
                r
                for r in results
                if r.row.opportunity == "e"
                and isinstance(r.row, it.ItalianRow)
                and r.row.target_class == it.DORSAL
            ]
            assert e_dorsal, name
            assert not any(r.boolean_predicted for r in e_dorsal), name

    def test_policies_order(
        self,
        real: tuple[tuple[it.ItalianRow, ...], dict[str, list[RowResult]]],
    ) -> None:
        rate = {
            name: sum(r.boolean_match for r in results) / len(results)
            for name, results in real[1].items()
        }
        assert rate[it.LEXICAL] > rate[it.STRESS] > rate[it.ALL_UNDERSPECIFIED]

    def test_sanity_rows_are_untouched(self) -> None:
        rows, coverage = it.load_rows(
            PROJECT_ROOT, classes=(it.CORONAL, it.FRONTSTEM)
        )
        assert coverage.validated > 1000
        for factory in it.PROFILE_VARIANTS.values():
            report = run_validation(factory(), rows)
            assert report.boolean_matches == report.total
            assert report.ipa_matches == report.total


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
