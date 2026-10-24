#!/usr/bin/env python3
"""Tests for the Czech palatalization analysis.

Pins (``phonology.analyses.czech_palatalization`` and
``src/phonology/inventories/czech_features.json``):

- **Inventory.** The 27 Romanian features in order; segments shared
  with Romanian copy its values exactly; the new segments derive from
  their bases as documented; bundles are unique before and after the
  patches; every new segment is justified in the notes and in
  ``docs/czech/phonology.md``; the underspecified segments differ from
  their bases in their cleared features only.
- **Tokenizer and transcription.** Diacritics stay on their base, ː
  attaches, allophones (r̝̊ ŋ ɱ ɣ d͡z) and cs.wiktionary symbols fold onto
  inventory labels; ``transcribe`` reproduces the cs-IPA transcription
  of lemmas as English Wiktionary gives them; the hard di ti ni of a
  lemma; the UR mode leaves the stem edge alone.
- **Golden derivations** (lexical policy): kluk+i, matka+e, kniha+e,
  vrah+i, Čech+i, doktor+i, biolog+ové, park+u, potok+e, moucha+e,
  sestra+e, plus the clusters (deska, český, německý, most, hvězda),
  the vocative (bůh, bratr) and the non-trigger front vowels; every
  derived segment is exactly an inventory bundle.
- **Policies and profiles.** lexical, all_underspecified, lexical_loo
  and spellout; build_profile() is lexical; the hooks are wired.
- **Rows.** ``iter_rows`` on a synthetic tree (filters, de-duplication,
  order, the lexeme join) and on the Czech fixture tree built by the
  data-layer tests (``tests/czech/treebuilder``); the driver and the
  summary on that tree.
- **Romanian untouched.** Importing the Czech analysis changes no
  Romanian default (the byte-level replay is
  ``tests/test_romanian_replay.py``).
- **Real data** (``needs_data``, ``slow``): the golden rows in the
  generated cell data.

Transcriptions quoted from English Wiktionary (CC BY-SA 4.0, via
wiktextract; Ylonen 2022) are as they appear in the generated data.
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path
from types import ModuleType

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import phonology.analyses.czech_palatalization as cz  # noqa: E402
import phonology.analyses.romanian_palatalization as romanian  # noqa: E402
from phonology.diagnostics.categorise import RefinedCategory  # noqa: E402
from phonology.diagnostics.compare import (  # noqa: E402
    ROMANIAN_NORMALISERS,
    compare_ipa,
)
from phonology.g2p import build_ur_from_row  # noqa: E402
from phonology.inventory import inventory_path  # noqa: E402
from phonology.lexicon import ClusterTag, LexRow  # noqa: E402
from phonology.profile import AnalysisProfile  # noqa: E402
from phonology.segments import (  # noqa: E402
    ROMANIAN_TOKENIZER,
    LengthHandling,
    project,
    tokenize_ipa,
)
from phonology.validation import (  # noqa: E402
    RowResult,
    evaluate_row,
    run_validation,
)

pytestmark = pytest.mark.cs

CZECH_JSON = PROJECT_ROOT / "src" / "phonology" / "inventories"
CZECH_JSON = CZECH_JSON / "czech_features.json"
ROMANIAN_JSON = PROJECT_ROOT / "romanian_features.json"
PHONOLOGY_DOC = PROJECT_ROOT / "docs" / "czech" / "phonology.md"
DRIVER = PROJECT_ROOT / "scripts" / "validate_lang_rules.py"
CELLS_CSV = PROJECT_ROOT / "data" / "czech" / "czech_cells.csv"

NEW_SEGMENTS = (
    "ɪ",
    "iː",
    "ɛ",
    "ɛː",
    "aː",
    "oː",
    "uː",
    "r̩",
    "l̩",
    "c",
    "ɟ",
    "ɲ",
    "x",
    "ɦ",
    "r̝",
    "u̯",
)


def _load(path: Path) -> dict[str, object]:
    raw: dict[str, object] = json.loads(path.read_text(encoding="utf-8"))
    return raw


def _segments(path: Path) -> dict[str, dict[str, str]]:
    segments = _load(path)["segments"]
    assert isinstance(segments, dict)
    return segments


def row(
    lemma: str,
    form: str,
    cell: str,
    ending: str,
    stem_final: str,
    *,
    orth: str = "",
    cluster: str = "",
    target: str = cz.VELAR,
    allomorph: str = cz.P,
    pos: str = "noun",
    gender: str = "m.an",
    mutation: bool | None = None,
    alt: Sequence[str] = ("NOM.PL", "LOC.PL"),
    resist: Sequence[str] = (),
    hard: str = "",
    exception: str | None = None,
) -> cz.CzechRow:
    """A row as ``iter_rows`` builds it (cells schema + lexeme columns)."""
    changed = bool(orth) or (target == cz.CORONAL and allomorph == cz.P)
    if mutation is None:
        mutation = changed
    if exception is None:
        exception = (
            ("undergoer" if mutation else "unexplained")
            if allomorph == cz.P
            else "non exception"
        )
    return cz.as_czech(
        cz.row_factory(
            {
                "lexeme_id": f"cs:{lemma}:{pos}:0:00000000",
                "lemma": lemma,
                "pos": pos,
                "cell": cell,
                "form": form,
                "ending": ending,
                "stem_final": stem_final,
                "cluster": cluster,
                "target_class": target,
                "allomorph_class": allomorph,
                "n_kind": "" if allomorph == cz.P else "allomorph",
                "orth_change": orth,
                "mutation": str(mutation),
                "gender_class": gender,
                "exception_label": exception,
                "alt_cells": "|".join(alt),
                "resist_cells": "|".join(resist),
                "hard_coronals": hard,
            }
        )
    )


# The golden rows of the task (forms as in the data) and more: the
# lexical UR, the derived SR, the attested (transcribed) form and the
# palatalization rules that fire.
KLUK = row("kluk", "kluci", "NOM.PL", "i", "k", orth="k→c")
MATKA = row("matka", "matce", "DAT.SG", "e", "k", orth="k→c", gender="f")
KNIHA = row("kniha", "knize", "DAT.SG", "e", "h", orth="h→z", gender="f")
VRAH = row("vrah", "vrazi", "NOM.PL", "i", "h", orth="h→z")
CECH = row("Čech", "Češi", "NOM.PL", "i", "ch", orth="ch→š")
DOKTOR = row(
    "doktor", "doktoři", "NOM.PL", "i", "r", orth="r→ř", target=cz.R_CLASS
)
BIOLOG = row(
    "biolog",
    "biologové",
    "NOM.PL",
    "ové",
    "g",
    allomorph=cz.N,
    alt=("LOC.PL",),
    resist=("NOM.PL",),
)
PARK = row(
    "park",
    "parku",
    "LOC.SG",
    "u",
    "k",
    allomorph=cz.N,
    gender="m.in",
    alt=("LOC.PL",),
    resist=("LOC.SG",),
)
POTOK = row("potok", "potoce", "LOC.SG", "e", "k", orth="k→c", gender="m.in")
MOUCHA = row("moucha", "mouše", "DAT.SG", "e", "ch", orth="ch→š", gender="f")
SESTRA = row(
    "sestra",
    "sestře",
    "DAT.SG",
    "e",
    "r",
    orth="r→ř",
    cluster="Cr",
    target=cz.R_CLASS,
    gender="f",
)
DESKA = row("deska", "desce", "DAT.SG", "e", "k", orth="sk→sc", cluster="sk")
CESKY = row(
    "český",
    "čeští",
    "NOM.PL.M.AN",
    "í",
    "k",
    orth="sk→šť",
    cluster="sk",
    pos="adj",
    gender="adj",
)
NEMECKY = row(
    "německý",
    "němečtí",
    "NOM.PL.M.AN",
    "í",
    "k",
    orth="ck→čť",
    cluster="ck",
    pos="adj",
    gender="adj",
)
MOST = row(
    "most",
    "mostě",
    "LOC.SG",
    "ě",
    "t",
    cluster="st",
    target=cz.CORONAL,
    gender="m.in",
)
HVEZDA = row(
    "hvězda",
    "hvězdě",
    "DAT.SG",
    "ě",
    "d",
    cluster="zd",
    target=cz.CORONAL,
    gender="f",
)
STUDENT = row("student", "studenti", "NOM.PL", "i", "t", target=cz.CORONAL)
BUH = row("bůh", "bože", "VOC.SG", "e", "h", orth="h→ž")
BRATR = row(
    "bratr",
    "bratře",
    "VOC.SG",
    "e",
    "r",
    orth="r→ř",
    cluster="Cr",
    target=cz.R_CLASS,
)
LIGA = row("liga", "lize", "DAT.SG", "e", "g", orth="g→z", gender="f")
POLITIKA = row(
    "politika", "politice", "DAT.SG", "e", "k", orth="k→c", hard="t"
)
BLIZOUCKY = row(
    "blizoučký",
    "blizoučcí",
    "NOM.PL.M.AN",
    "í",
    "k",
    orth="k→c",
    pos="adj",
    gender="adj",
)
MEKKY = row(
    "měkký",
    "měkcí",
    "NOM.PL.M.AN",
    "í",
    "k",
    orth="k→c",
    pos="adj",
    gender="adj",
)
KLUKEM = row("kluk", "klukem", "INS.SG", "em", "k", allomorph=cz.N)
STUDENTECH = row(
    "student",
    "studentech",
    "LOC.PL",
    "ech",
    "t",
    target=cz.CORONAL,
    allomorph=cz.N,
)
DOKTORE = row(
    "doktor",
    "doktore",
    "VOC.SG",
    "e",
    "r",
    target=cz.R_CLASS,
    allomorph=cz.N,
)

GOLDEN: tuple[tuple[cz.CzechRow, str, str, str, frozenset[str]], ...] = (
    (KLUK, "kluK PAL1 ɪ", "klut͡sɪ", "ˈklut͡sɪ", frozenset({"velar-pal"})),
    (MATKA, "matK PAL1 ɛ", "matt͡sɛ", "ˈmat͡sɛ", frozenset({"velar-pal"})),
    (KNIHA, "kɲɪꞪ PAL1 ɛ", "kɲɪzɛ", "ˈkɲɪzɛ", frozenset({"velar-pal"})),
    (VRAH, "vraꞪ PAL1 ɪ", "vrazɪ", "ˈvrazɪ", frozenset({"velar-pal"})),
    (CECH, "t͡ʃɛX PAL1 ɪ", "t͡ʃɛʃɪ", "ˈt͡ʃɛʃɪ", frozenset({"x-pal"})),
    (
        DOKTOR,
        "doktoR PAL1 ɪ",
        "doktor̝ɪ",
        "ˈdoktor̝ɪ",
        frozenset({"coronal-pal"}),
    ),
    (BIOLOG, "bɪjoloGovɛː", "bɪjoloɡovɛː", "ˈbɪjoloɡovɛː", frozenset()),
    (PARK, "parKu", "parku", "ˈparku", frozenset()),
    (POTOK, "potoK PAL1 ɛ", "potot͡sɛ", "ˈpotot͡sɛ", frozenset({"velar-pal"})),
    (MOUCHA, "mou̯X PAL1 ɛ", "mou̯ʃɛ", "ˈmou̯ʃɛ", frozenset({"x-pal"})),
    (
        SESTRA,
        "sɛstR PAL1 ɛ",
        "sɛstr̝ɛ",
        "ˈsɛstr̝̊ɛ",
        frozenset({"coronal-pal"}),
    ),
    (DESKA, "dɛsK PAL1 ɛ", "dɛst͡sɛ", "ˈdɛst͡sɛ", frozenset({"velar-pal"})),
    (
        CESKY,
        "t͡ʃɛSḰ PAL1 iː",
        "t͡ʃɛʃciː",
        "ˈt͡ʃɛʃciː",
        frozenset({"coronal-pal", "s-pal"}),
    ),
    (
        NEMECKY,
        "ɲɛmɛTSḰ PAL1 iː",
        "ɲɛmɛt͡ʃciː",
        "ˈɲɛmɛt͡ʃciː",
        frozenset({"coronal-pal", "s-pal"}),
    ),
    (MOST, "mosT PAL1 ɛ", "moscɛ", "ˈmoscɛ", frozenset({"coronal-pal"})),
    (
        HVEZDA,
        "ɦvjɛzD PAL1 ɛ",
        "ɦvjɛzɟɛ",
        "ˈɦvjɛzɟɛ",
        frozenset({"coronal-pal"}),
    ),
    (
        STUDENT,
        "studɛnT PAL1 ɪ",
        "studɛncɪ",
        "ˈstudɛncɪ",
        frozenset({"coronal-pal"}),
    ),
    (BUH, "boꞪ PAL3 ɛ", "boʒɛ", "ˈboʒɛ", frozenset({"big-pal"})),
    (
        BRATR,
        "bratR PAL1 ɛ",
        "bratr̝ɛ",
        "ˈbratr̝̊ɛ",
        frozenset({"coronal-pal"}),
    ),
    (LIGA, "lɪG PAL1 ɛ", "lɪzɛ", "ˈlɪzɛ", frozenset({"velar-pal"})),
    (
        POLITIKA,
        "polɪtɪK PAL1 ɛ",
        "polɪtɪt͡sɛ",
        "ˈpolɪtɪt͡sɛ",
        frozenset({"velar-pal"}),
    ),
    (
        BLIZOUCKY,
        "blɪzou̯t͡ʃK PAL1 iː",
        "blɪzou̯t͡ʃt͡siː",
        "ˈblɪzou̯t͡ʃt͡siː",
        frozenset({"velar-pal"}),
    ),
    (
        MEKKY,
        "mɲɛkK PAL1 iː",
        "mɲɛkt͡siː",
        "ˈmɲɛkt͡siː",
        frozenset({"velar-pal"}),
    ),
    (KLUKEM, "kluKɛm", "klukɛm", "ˈklukɛm", frozenset()),
    (STUDENTECH, "studɛnTɛx", "studɛntɛx", "ˈstudɛntɛx", frozenset()),
    (DOKTORE, "doktoRɛ", "doktorɛ", "ˈdoktorɛ", frozenset()),
)


def _ur_text(labels: Sequence[str]) -> str:
    """Labels joined, the palatalizer set off by spaces (as in GOLDEN)."""
    return "".join(
        f" {x} " if x in cz.PALATALIZER_FEATURES else x for x in labels
    ).strip()


def _fired(result: RowResult) -> frozenset[str]:
    assert result.derivation is not None
    return frozenset(
        step.rule
        for step in result.derivation.steps
        if step.fired and step.rule in cz.PALATALIZATION_RULE_NAMES
    )


@pytest.fixture(scope="module")
def lexical() -> AnalysisProfile:
    return cz.build_profile()


@pytest.fixture(scope="module")
def profiles() -> dict[str, AnalysisProfile]:
    return {name: factory() for name, factory in cz.PROFILE_VARIANTS.items()}


# ---------------------------------------------------------------------------
# Inventory
# ---------------------------------------------------------------------------


class TestInventory:
    def test_file_is_the_registered_inventory(self) -> None:
        assert inventory_path("czech") == CZECH_JSON
        assert cz.build_profile().inventory_json == CZECH_JSON

    def test_features_are_the_romanian_list(self) -> None:
        assert (
            _load(CZECH_JSON)["features"] == _load(ROMANIAN_JSON)["features"]
        )
        assert list(_load(CZECH_JSON)) == ["metadata", "features", "segments"]

    def test_shared_segments_copy_romanian_exactly(self) -> None:
        czech, rom = _segments(CZECH_JSON), _segments(ROMANIAN_JSON)
        shared = set(czech) & set(rom)
        assert shared == set(
            "a o u p b t d k ɡ t͡s t͡ʃ d͡ʒ f v s z ʃ ʒ m n l r j".split()
        )
        for label in shared:
            assert czech[label] == rom[label], label

    def test_new_segments_are_the_declared_ones(self) -> None:
        czech, rom = _segments(CZECH_JSON), _segments(ROMANIAN_JSON)
        assert set(czech) - set(rom) == set(NEW_SEGMENTS)
        assert len(czech) == 39

    @pytest.mark.parametrize(
        ("label", "base", "changes"),
        [
            ("ɪ", "i", {}),
            ("iː", "i", {"Long": "+"}),
            ("ɛ", "e", {}),
            ("ɛː", "e", {"Long": "+"}),
            ("aː", "a", {"Long": "+"}),
            ("oː", "o", {"Long": "+"}),
            ("uː", "u", {"Long": "+"}),
            ("r̩", "r", {"Syllabic": "+"}),
            ("l̩", "l", {"Syllabic": "+"}),
            ("c", "t", {"Anterior": "-", "Distributed": "+"}),
            ("ɟ", "d", {"Anterior": "-", "Distributed": "+"}),
            ("ɲ", "n", {"Anterior": "-", "Distributed": "+"}),
            ("x", "k", {"Continuant": "+", "DelRel": "0"}),
            ("ɦ", "k", {"Continuant": "+", "DelRel": "0", "Voice": "+"}),
            (
                "r̝",
                "r",
                {
                    "Sonorant": "-",
                    "Approximant": "-",
                    "Anterior": "-",
                    "Distributed": "+",
                },
            ),
            ("u̯", "w", {}),
        ],
    )
    def test_new_segments_derive_from_their_bases(
        self, label: str, base: str, changes: Mapping[str, str]
    ) -> None:
        expected = dict(_segments(ROMANIAN_JSON)[base])
        expected.update(changes)
        assert _segments(CZECH_JSON)[label] == expected

    @pytest.mark.parametrize("patched", [False, True])
    def test_bundles_are_unique(self, patched: bool) -> None:
        segments: Mapping[str, Mapping[str, str]] = (
            cz.load_inventory().base if patched else _segments(CZECH_JSON)
        )
        bundles = [tuple(sorted(v.items())) for v in segments.values()]
        assert len(bundles) == len(set(bundles))

    def test_new_segments_are_justified(self) -> None:
        notes = _load(CZECH_JSON)["metadata"]
        assert isinstance(notes, dict)
        doc = PHONOLOGY_DOC.read_text(encoding="utf-8")
        for label in NEW_SEGMENTS:
            assert label in notes["notes"], label
            assert label in doc, label

    def test_every_label_is_a_token(self) -> None:
        for label in _segments(CZECH_JSON):
            assert cz.czech_tokens(label) == (label,), label

    def test_patches(self) -> None:
        inv = cz.load_inventory()
        for label in ("k", "ɡ", "x", "ɦ", "t", "d", "c", "ɟ"):
            assert inv.base[label]["Strident"] == "-", label
        assert {p.segment for p in cz.PATCHES} == {
            "k",
            "ɡ",
            "x",
            "ɦ",
            "t",
            "d",
            "c",
            "ɟ",
        }

    def test_underspecified_segments_clear_only_their_features(self) -> None:
        inv = cz.load_inventory()
        assert not set(inv.underspec) & set(inv.base)
        for u in cz.UNDERSPEC:
            twin, base = inv.features_of(u.label), inv.base[u.base]
            for feature, value in base.items():
                expected = "0" if feature in u.clear else value
                assert twin[feature] == expected, (u.label, feature)
        assert inv.features_of("Ḱ")["Strident"] == "-"
        assert inv.features_of("K")["Strident"] == "0"

    def test_palatalizers_are_no_segment(self) -> None:
        inv = cz.load_inventory()
        for label in cz.PALATALIZER_FEATURES:
            pal = cz.palatalizer(label)
            assert cz.ANY_PALATALIZER.contains(pal)
            for name in list(inv.base) + list(inv.underspec):
                assert not cz.ANY_PALATALIZER.contains(inv.segment(name))
        assert cz.SMALL_PALATALIZER.contains(cz.palatalizer(cz.PAL1))
        assert not cz.SMALL_PALATALIZER.contains(cz.palatalizer(cz.PAL3))
        assert cz.BIG_PALATALIZER.contains(cz.palatalizer(cz.PAL3))


# ---------------------------------------------------------------------------
# Tokenizer, transcription and comparison
# ---------------------------------------------------------------------------


class TestTokenizer:
    def test_config(self) -> None:
        tok = cz.CZECH_TOKENIZER
        assert tok.keep_combining is True
        assert tok.length is LengthHandling.ATTACH
        assert {"ʔ", "ˈ", "."} <= tok.skip

    @pytest.mark.parametrize(
        ("ipa", "expected"),
        [
            ("ˈt͡ʃɛʃciː", ("t͡ʃ", "ɛ", "ʃ", "c", "iː")),
            ("ˈmou̯ʃɛ", ("m", "o", "u̯", "ʃ", "ɛ")),
            ("ˈsɛstr̝̊ɛ", ("s", "ɛ", "s", "t", "r̝", "ɛ")),
            ("[ˈbaŋka]", ("b", "a", "n", "k", "a")),
            ("[ˈvl̩k]", ("v", "l̩", "k")),
            (
                "[ˈʔa.bɪ.gaj.lɪn]",
                ("a", "b", "ɪ", "ɡ", "a", "j", "l", "ɪ", "n"),
            ),
            ("[ˈstɔdvɔʊ̯]", ("s", "t", "o", "d", "v", "o", "u̯")),
            ("[ˈsald͡zbur]", ("s", "a", "l", "t͡s", "b", "u", "r")),
            ("[ˈdɛɣbɛ]", ("d", "ɛ", "ɦ", "b", "ɛ")),
        ],
    )
    def test_tokens(self, ipa: str, expected: tuple[str, ...]) -> None:
        assert cz.czech_tokens(ipa) == expected

    def test_romanian_default_is_unchanged(self) -> None:
        assert tokenize_ipa("ˈmus.kə") == ("m", "u", "s", "k", "ə")
        assert tokenize_ipa("r̝") == ("r", "̝")
        assert dict(ROMANIAN_TOKENIZER.fold) == {"ɛ": "e", "ɔ": "o", "g": "ɡ"}


class TestTranscription:
    @pytest.mark.parametrize(
        ("spelling", "ipa"),
        [
            # English Wiktionary's transcriptions of these lemmas.
            ("kluk", "kluk"),
            ("kniha", "kɲɪɦa"),
            ("moucha", "mou̯xa"),
            ("vrah", "vrax"),
            ("biolog", "bɪjolok"),
            ("sestra", "sɛstra"),
            ("hvězda", "ɦvjɛzda"),
            ("německý", "ɲɛmɛt͡skiː"),
            ("vlk", "vl̩k"),
            ("bratr", "bratr̩"),
            ("hamr", "ɦamr"),
            ("keř", "kɛr̝̊"),
            ("lehký", "lɛxkiː"),
            ("banka", "baŋka"),
            ("chřipka", "xr̝̊ɪpka"),
            ("dcera", "t͡sɛra"),
            ("svatba", "svadba"),
            ("odstavec", "otstavɛt͡s"),
            ("biologie", "bɪjoloɡɪjɛ"),
            ("ticho", "cɪxo"),
            ("anglofonní", "aŋɡlofoɲiː"),
            ("apačština", "apat͡ʃcɪna"),
            ("hrdliččí", "ɦr̩dlɪt͡ʃt͡ʃiː"),
            ("debř", "dɛpr̝̊"),
            ("abchazský", "apxaskiː"),
        ],
    )
    def test_cs_ipa_conventions(self, spelling: str, ipa: str) -> None:
        transcription = cz.transcribe(spelling)
        assert "".join(transcription.segments) == ipa
        assert transcription.render() == "ˈ" + ipa

    @pytest.mark.parametrize(
        ("lemma", "ipa", "hard"),
        [
            ("diplomat", "[ˈdɪplomat]", {"d"}),
            ("politika", "[ˈpolɪtɪka]", {"t"}),
            ("dialektický", "[ˈdɪjalɛktɪt͡skiː]", {"d", "t"}),
            ("kino", "[ˈkɪno]", set()),
            ("ticho", "[ˈcɪxo]", set()),
            ("student", "", set()),
        ],
    )
    def test_hard_coronals(self, lemma: str, ipa: str, hard: set[str]) -> None:
        assert cz.hard_coronals(lemma, ipa) == frozenset(hard)

    def test_hard_letters_spare_the_stem_edge(self) -> None:
        hard = frozenset({"t"})
        assert cz.transcribe("politika", hard=hard).render() == "ˈpolɪtɪka"
        # The stem-final t before the ending -i keeps the soft reading.
        edge = cz.transcribe("diplomati", boundary=8, hard=frozenset("dt"))
        assert "".join(edge.segments) == "dɪplomacɪ"

    def test_underlying_mode_leaves_the_stem_edge(self) -> None:
        assert "".join(cz.transcribe("mostě").segments) == "moscɛ"
        ur = cz.transcribe("mostě", boundary=4, underlying=True)
        assert "".join(ur.segments) == "mostɛ"
        assert ur.split(4) == (("m", "o", "s", "t"), ("ɛ",))
        # ...and does not merge the stem final with its neighbour.
        assert "".join(cz.transcribe("měkkí").segments) == "mɲɛkiː"
        kept = cz.transcribe("měkkí", boundary=4, underlying=True)
        assert kept.split(4)[0] == ("m", "ɲ", "ɛ", "k", "k")

    def test_every_transcribed_token_is_a_label(self) -> None:
        inv = cz.load_inventory()
        words = [r.plural for r, *_ in GOLDEN] + [
            "chřipka",
            "dechberoucí",
            "salcburský",
            "hrdliččí",
            "anglofonní",
            "biologie",
        ]
        for word in words:
            for segment in cz.transcribe(word).segments:
                assert cz.phonemic(segment) in inv.base, (word, segment)

    def test_strip_and_cascade(self) -> None:
        assert cz.strip_czech("matt͡sɛ") == cz.strip_czech("ˈmat͡sɛ")
        assert cz.strip_czech("sɛstr̝ɛ") == cz.strip_czech("[ˈsɛstr̝̊ɛ]")
        assert cz.strip_czech("ˈbaŋka") == "banka"
        assert cz.strip_czech("ɲɛmɛt͡ʃʃciː") == "ɲɛmɛt͡ʃciː"
        # Length is contrastive: never stripped.
        assert cz.strip_czech("dráze") != cz.strip_czech("draze")
        assert cz.strip_czech("ˈdraːzɛ") != cz.strip_czech("ˈdrazɛ")
        assert compare_ipa(
            "matt͡sɛ", "ˈmat͡sɛ", normalisers=cz.CZECH_NORMALISERS
        ).matched
        assert not compare_ipa(
            "kluki", "ˈklut͡sɪ", normalisers=cz.CZECH_NORMALISERS
        ).matched
        assert cz.CZECH_NORMALISERS.steps == ()


# ---------------------------------------------------------------------------
# Golden derivations
# ---------------------------------------------------------------------------


class TestGoldenDerivations:
    @pytest.mark.parametrize(
        ("czech_row", "ur", "sr", "attested", "fired"),
        GOLDEN,
        ids=[f"{g[0].lemma}+{g[0].ending}" for g in GOLDEN],
    )
    def test_lexical(
        self,
        lexical: AnalysisProfile,
        czech_row: cz.CzechRow,
        ur: str,
        sr: str,
        attested: str,
        fired: frozenset[str],
    ) -> None:
        labels = cz.ur_labels(czech_row, cz.LEXICAL)
        assert labels is not None
        assert _ur_text(labels) == ur
        assert czech_row.ipa_pl == attested
        result = evaluate_row(czech_row, lexical)
        assert result.predicted_sr == sr
        assert _fired(result) == fired
        assert result.boolean_predicted is bool(fired)
        assert result.boolean_match and result.ipa_match
        assert cz.expected_firings(czech_row, cz.LEXICAL) == fired

    @pytest.mark.parametrize(
        "czech_row", [g[0] for g in GOLDEN], ids=lambda r: r.plural
    )
    def test_derived_segments_are_inventory_bundles(
        self, lexical: AnalysisProfile, czech_row: cz.CzechRow
    ) -> None:
        inv = lexical.inventory
        result = evaluate_row(czech_row, lexical)
        assert result.derivation is not None
        for segment in result.derivation.sr:
            label = project(segment, inv.base)
            assert dict(segment.features) == dict(inv.base[label]), label

    def test_matka_is_degeminated_by_the_comparison_only(
        self, lexical: AnalysisProfile
    ) -> None:
        result = evaluate_row(MATKA, lexical)
        assert result.predicted_sr == "matt͡sɛ"
        assert result.attested_first == "ˈmat͡sɛ"
        assert result.ipa_match

    def test_non_trigger_front_vowels(self, lexical: AnalysisProfile) -> None:
        # -em, -ech and the coronal vocative -e carry no palatalizer.
        for czech_row in (KLUKEM, STUDENTECH, DOKTORE):
            assert cz.palatalizer_of(czech_row) is None
            result = evaluate_row(czech_row, lexical)
            assert not result.boolean_predicted
            assert result.boolean_match and result.ipa_match

    def test_prespecified_targets_are_inalterable(
        self, lexical: AnalysisProfile
    ) -> None:
        # No alternation evidence: under the lexical policy the target is
        # prespecified, and the palatalizer fails to unify into it.
        resister = row(
            "hater",
            "hateri",
            "NOM.PL",
            "i",
            "r",
            target=cz.R_CLASS,
            alt=(),
            mutation=False,
        )
        labels = cz.ur_labels(resister, cz.LEXICAL)
        assert labels is not None and "R" not in labels and "r" in labels
        result = evaluate_row(resister, lexical)
        assert result.predicted_sr == "ɦatɛrɪ"
        assert not result.boolean_predicted and result.boolean_match
        everything = cz.build_variant(cz.ALL_UNDERSPECIFIED)
        assert evaluate_row(resister, everything).predicted_sr == "ɦatɛr̝ɪ"

    def test_x_and_h_split_on_voicing(self, lexical: AnalysisProfile) -> None:
        assert evaluate_row(CECH, lexical).predicted_sr == "t͡ʃɛʃɪ"
        assert evaluate_row(VRAH, lexical).predicted_sr == "vrazɪ"

    def test_big_palatalization_only_with_pal3(
        self, lexical: AnalysisProfile
    ) -> None:
        assert cz.palatalizer_of(BUH) == cz.PAL3
        assert cz.palatalizer_of(BRATR) == cz.PAL1
        assert cz.palatalizer_of(KLUK) == cz.PAL1
        hrisnik = row("hříšník", "hříšníče", "VOC.SG", "e", "k", orth="k→č")
        result = evaluate_row(hrisnik, lexical)
        assert result.predicted_sr == "ɦr̝iːʃɲiːt͡ʃɛ"
        assert _fired(result) == {"big-pal"}
        assert result.ipa_match


# ---------------------------------------------------------------------------
# Policies
# ---------------------------------------------------------------------------


class TestPolicies:
    def test_names(self) -> None:
        assert cz.POLICIES == (
            "lexical",
            "all_underspecified",
            "lexical_loo",
            "spellout",
        )

    def test_lexical(self) -> None:
        assert cz.is_underspecified(KLUK, cz.LEXICAL)
        assert not cz.is_underspecified(
            row("x", "xi", "NOM.PL", "i", "k", alt=()), cz.LEXICAL
        )

    def test_all_underspecified(self) -> None:
        assert cz.is_underspecified(
            row("x", "xi", "NOM.PL", "i", "k", alt=()), cz.ALL_UNDERSPECIFIED
        )

    def test_leave_one_cell_out(self) -> None:
        # The tested cell and its copies are not evidence.
        feminine = row(
            "matka",
            "matce",
            "DAT.SG",
            "e",
            "k",
            orth="k→c",
            gender="f",
            alt=("DAT.SG", "LOC.SG"),
        )
        assert not cz.is_underspecified(feminine, cz.LEXICAL_LOO)
        nom = row(
            "kluk", "kluci", "NOM.PL", "i", "k", alt=("NOM.PL", "VOC.PL")
        )
        assert not cz.is_underspecified(nom, cz.LEXICAL_LOO)
        assert cz.is_underspecified(KLUK, cz.LEXICAL_LOO)
        assert cz.copy_group("LOC.SG", "m.in") == {"LOC.SG"}
        assert cz.copy_group("LOC.SG", "f") == {"DAT.SG", "LOC.SG"}

    def test_spellout(self) -> None:
        # An N allomorph anywhere marks the root prespecified.
        assert not cz.is_underspecified(BIOLOG, cz.SPELLOUT)
        assert not cz.is_underspecified(PARK, cz.SPELLOUT)
        assert cz.is_underspecified(KLUK, cz.SPELLOUT)

    def test_unknown_policy(self) -> None:
        with pytest.raises(ValueError):
            cz.is_underspecified(KLUK, "stress")
        with pytest.raises(ValueError):
            cz.build_variant("stress")

    def test_mixed_stem_under_spellout(
        self, profiles: dict[str, AnalysisProfile]
    ) -> None:
        parcich = row(
            "park",
            "parcích",
            "LOC.PL",
            "ích",
            "k",
            orth="k→c",
            gender="m.in",
            alt=("LOC.PL",),
            resist=("LOC.SG",),
        )
        assert evaluate_row(parcich, profiles[cz.LEXICAL]).boolean_match
        spellout = evaluate_row(parcich, profiles[cz.SPELLOUT])
        assert spellout.predicted_sr == "parkiːx"
        assert not spellout.boolean_match
        assert (
            cz.mismatch_class(spellout)
            is cz.MismatchClass.MIXED_STEM_ALTERNATES
        )
        loo = evaluate_row(parcich, profiles[cz.LEXICAL_LOO])
        assert cz.mismatch_class(loo) is cz.MismatchClass.MIXED_STEM_ALTERNATES
        feminine = evaluate_row(
            row(
                "kniha",
                "knize",
                "DAT.SG",
                "e",
                "h",
                orth="h→z",
                gender="f",
                alt=("DAT.SG", "LOC.SG"),
            ),
            profiles[cz.LEXICAL_LOO],
        )
        assert (
            cz.mismatch_class(feminine)
            is cz.MismatchClass.NO_EVIDENCE_ALTERNATES
        )

    def test_fallback_applies_the_policy_not_the_data(self) -> None:
        unexplained = row(
            "nepsycholog",
            "nepsychologi",
            "NOM.PL",
            "i",
            "g",
            mutation=False,
        )
        assert cz.FallbackPredictor(cz.LEXICAL)(unexplained)
        assert not cz.FallbackPredictor(cz.LEXICAL)(BIOLOG)
        assert cz.ExpectedFirings(cz.LEXICAL)(DESKA) == {"velar-pal"}

    def test_unexplained_form_class(self, lexical: AnalysisProfile) -> None:
        unexplained = row(
            "nepsycholog",
            "nepsychologi",
            "NOM.PL",
            "i",
            "g",
            mutation=False,
        )
        result = evaluate_row(unexplained, lexical)
        assert result.boolean_predicted and not result.boolean_match
        assert cz.mismatch_class(result) is cz.MismatchClass.UNEXPLAINED_FORM


# ---------------------------------------------------------------------------
# Profiles and hooks
# ---------------------------------------------------------------------------


class TestProfiles:
    def test_build_profile_is_the_lexical_variant(
        self, lexical: AnalysisProfile
    ) -> None:
        assert lexical.name == "czech-palatalization (lexical)"
        assert cz.PROFILE_VARIANTS[cz.LEXICAL] is cz.build_profile
        assert set(cz.PROFILE_VARIANTS) == set(cz.POLICIES)

    def test_variants_differ_in_policy_only(
        self, profiles: dict[str, AnalysisProfile]
    ) -> None:
        for name, profile in profiles.items():
            assert profile.name == f"czech-palatalization ({name})"
            assert profile.rules == cz.RULES
            assert profile.ur_builder == cz.UrBuilder(name)
            assert profile.fallback_predict == cz.FallbackPredictor(name)

    def test_hooks(self, lexical: AnalysisProfile) -> None:
        assert lexical.tokenizer is cz.CZECH_TOKENIZER
        assert lexical.normalisers is cz.CZECH_NORMALISERS
        assert lexical.row_factory is cz.row_factory
        assert lexical.categoriser is cz.categorise
        assert lexical.cluster_key is cz.cluster_key
        assert lexical.target_stems == {
            "k",
            "g",
            "h",
            "ch",
            "r",
            "t",
            "d",
            "n",
        }
        assert lexical.palatalization_rule_names == {
            "big-pal",
            "x-pal",
            "velar-pal",
            "coronal-pal",
            "s-pal",
        }

    def test_rule_order(self) -> None:
        names = [r.name for r in cz.RULES]
        assert names.index("x-pal") < names.index("velar-pal")
        assert names.index("velar-pal") < names.index("coronal-pal")
        assert names.index("coronal-pal") < names.index("s-pal")
        assert names.index("s-pal") < names.index("pal-deletion")
        assert names.index("pal-deletion") < names.index("dorsal-default")

    def test_cluster_key(self) -> None:
        assert cz.cluster_key(CESKY) == "sk+í (suffix)"
        assert cz.cluster_key(DESKA) == "sk+e"
        assert cz.cluster_key(BUH) == "VOC velar"
        assert cz.cluster_key(KLUK) is None

    def test_categoriser_uses_the_shared_vocabulary(
        self, lexical: AnalysisProfile
    ) -> None:
        cat = cz.categorise(
            KLUK,
            "kluki",
            KLUK.ipa_pl,
            False,
            ur_built=True,
            normalisers=cz.CZECH_NORMALISERS,
        )
        assert cat.category is RefinedCategory.RULE_FAILURE
        failed = cz.categorise(
            KLUK,
            "",
            "",
            True,
            ur_built=False,
            normalisers=cz.CZECH_NORMALISERS,
        )
        assert failed.category is RefinedCategory.UR_BUILD_FAILED

    def test_plain_lexrow_is_refused(self) -> None:
        plain = LexRow(
            lemma="kluk",
            plural="kluci",
            pos="N",
            stem_final="k",
            opportunity="i",
            mutation=True,
            ipa_lemma="",
            ipa_pl="",
            exception_reason="",
            cluster=ClusterTag.NONE,
        )
        with pytest.raises(TypeError):
            cz.ur_labels(plain, cz.LEXICAL)

    def test_row_factory_fields(self) -> None:
        assert KLUK.opportunity == "NOM.PL -i"
        assert KLUK.pos == "N" and CESKY.pos == "ADJ"
        assert KLUK.stem_spelling == "kluc"
        assert cz.base_stem(CESKY) == "česk"
        assert cz.base_stem(NEMECKY) == "německ"
        assert CESKY.suffix_cluster and not DESKA.suffix_cluster
        assert KLUK.cluster is ClusterTag.NONE
        broken = row("kluk", "kluki", "NOM.PL", "i", "k", orth="k→c")
        assert cz.base_stem(broken) is None
        assert cz.ur_labels(broken, cz.LEXICAL) is None

    def test_romanian_defaults_survive_the_czech_import(self) -> None:
        profile = romanian.build_profile()
        assert profile.tokenizer is ROMANIAN_TOKENIZER
        assert profile.normalisers is ROMANIAN_NORMALISERS
        assert profile.ur_builder is build_ur_from_row
        assert profile.inventory_json == ROMANIAN_JSON
        assert tokenize_ipa("ˈmuʃ.te") == ("m", "u", "ʃ", "t", "e")


# ---------------------------------------------------------------------------
# Rows: a synthetic tree and the Czech fixture tree
# ---------------------------------------------------------------------------

CELL_FIELDS = (
    "lexeme_id",
    "lemma",
    "pos",
    "cell",
    "form",
    "primary",
    "aligned",
    "curated_status",
    "stem_final",
    "cluster",
    "target_class",
    "gender_class",
    "ending",
    "orth_change",
    "mutation",
    "allomorph_class",
    "n_kind",
    "exception_label",
)
LEXEME_FIELDS = (
    "lexeme_id",
    "lemma",
    "pos",
    "domain",
    "target_class",
    "ipa",
    "p_cells_E0_dictionary",
    "n_cells_E0_dictionary",
)


def _cell(
    lexeme: str, lemma: str, cell: str, form: str, **kw: str
) -> dict[str, str]:
    base = {
        "lexeme_id": lexeme,
        "lemma": lemma,
        "pos": "noun",
        "cell": cell,
        "form": form,
        "primary": "True",
        "aligned": "True",
        "curated_status": "",
        "stem_final": "k",
        "cluster": "",
        "target_class": "velar",
        "gender_class": "m.an",
        "ending": "i",
        "orth_change": "k→c",
        "mutation": "True",
        "allomorph_class": "P",
        "n_kind": "",
        "exception_label": "undergoer",
    }
    base.update(kw)
    return base


def _write(
    path: Path, fields: Sequence[str], rows: Sequence[dict[str, str]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


@pytest.fixture
def synthetic(tmp_path: Path) -> Path:
    cells, complete = cz.data_paths(tmp_path)
    kluk = "cs:kluk:noun:0:aaaaaaaa"
    dipl = "cs:diplomat:noun:0:bbbbbbbb"
    _write(
        cells,
        CELL_FIELDS,
        [
            _cell(kluk, "kluk", "NOM.PL", "kluci"),
            _cell(kluk, "kluk", "NOM.PL", "kluci"),  # another table
            _cell(kluk, "kluk", "LOC.PL", "klucích", ending="ích"),
            _cell(
                kluk,
                "kluk",
                "INS.SG",
                "klukem",
                ending="em",
                orth_change="",
                mutation="False",
                allomorph_class="N",
                n_kind="nontrigger",
                exception_label="non exception",
            ),
            _cell(kluk, "kluk", "NOM.PL", "klucí", primary="False"),
            _cell(kluk, "kluk", "GEN.SG", "kluka", allomorph_class="none"),
            _cell(kluk, "kluk", "DAT.SG", "kluku", aligned="False"),
            _cell(kluk, "kluk", "LOC.SG", "kluci", curated_status="rejected"),
            _cell(
                "cs:kmen:noun:0:cccccccc",
                "kmen",
                "LOC.SG",
                "kmeni",
                stem_final="n",
                target_class="coronal_inferred",
                allomorph_class="N",
                n_kind="soft",
            ),
            _cell(
                "cs:10letý:adj:0:dddddddd",
                "10letý",
                "NOM.PL.M.AN",
                "10letí",
                pos="adj",
                stem_final="t",
                target_class="coronal_inferred",
            ),
            _cell(
                "cs:ryba:noun:0:eeeeeeee",
                "ryba",
                "DAT.SG",
                "rybě",
                stem_final="b",
                target_class="labial",
            ),
            _cell(
                dipl,
                "diplomat",
                "NOM.PL",
                "diplomati",
                stem_final="t",
                target_class="coronal_inferred",
                orth_change="",
            ),
        ],
    )
    _write(
        complete,
        LEXEME_FIELDS,
        [
            {
                "lexeme_id": kluk,
                "lemma": "kluk",
                "pos": "noun",
                "domain": "in",
                "target_class": "velar",
                "ipa": "[ˈkluk]",
                "p_cells_E0_dictionary": "NOM.PL|LOC.PL",
                "n_cells_E0_dictionary": "",
            },
            {
                "lexeme_id": dipl,
                "lemma": "diplomat",
                "pos": "noun",
                "domain": "in",
                "target_class": "coronal_inferred",
                "ipa": "[ˈdɪplomat]",
                "p_cells_E0_dictionary": "NOM.PL",
                "n_cells_E0_dictionary": "NOM.PL",
            },
        ],
    )
    return tmp_path


class TestRows:
    def test_filters_dedupe_order_and_join(self, synthetic: Path) -> None:
        rows = cz.load_rows(synthetic)
        assert [(r.lemma, r.cell, r.plural) for r in rows] == [
            ("diplomat", "NOM.PL", "diplomati"),
            ("kluk", "INS.SG", "klukem"),
            ("kluk", "LOC.PL", "klucích"),
            ("kluk", "NOM.PL", "kluci"),
        ]
        assert list(cz.iter_rows(synthetic)) == list(rows)
        kluci = rows[-1]
        assert kluci.alt_cells == ("NOM.PL", "LOC.PL")
        assert kluci.ipa_lemma == "[ˈkluk]"
        assert kluci.ipa_pl == "ˈklut͡sɪ"
        diplomati = rows[0]
        assert diplomati.hard_coronals == {"d"}
        assert diplomati.resist_cells == ("NOM.PL",)
        assert diplomati.ipa_pl == "ˈdɪplomacɪ"

    def test_validation_on_the_synthetic_tree(self, synthetic: Path) -> None:
        report = run_validation(cz.build_profile(), cz.iter_rows(synthetic))
        assert report.total == 4
        assert report.boolean_matches == report.ipa_matches == 4


@pytest.fixture(scope="module")
def fixture_tree(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The Czech data layer run on the tracked Wiktionary fixtures."""
    from tests.czech import treebuilder

    root = treebuilder.build_tree(tmp_path_factory.mktemp("czech_phonology"))
    treebuilder.run_stages(root)
    return root


def _driver() -> ModuleType:
    spec = importlib.util.spec_from_file_location("driver", DRIVER)
    assert spec is not None and spec.loader is not None
    driver = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(driver)
    return driver


class TestFixtureTree:
    def test_rows(self, fixture_tree: Path) -> None:
        rows = cz.load_rows(fixture_tree)
        keys = [(r.lexeme_id, r.cell, r.plural) for r in rows]
        assert keys == sorted(set(keys))
        lemmas = {r.lemma for r in rows}
        assert {"kluk", "kniha", "biolog", "park", "Čech", "student"} <= lemmas
        assert all(r.target_class in cz.TARGET_CLASSES for r in rows)
        assert all(r.allomorph in (cz.P, cz.N) for r in rows)

    def test_lexical_and_all_underspecified_match(
        self, fixture_tree: Path
    ) -> None:
        rows = cz.load_rows(fixture_tree)
        for policy in (cz.LEXICAL, cz.ALL_UNDERSPECIFIED):
            report = run_validation(cz.build_variant(policy), rows)
            assert report.boolean_matches == report.total, policy
            assert report.ipa_matches == report.total, policy

    def test_leave_one_cell_out_misses_single_cell_lexemes(
        self, fixture_tree: Path
    ) -> None:
        rows = cz.load_rows(fixture_tree)
        report = run_validation(cz.build_variant(cz.LEXICAL_LOO), rows)
        failed = {
            (r.row.lemma, r.row.plural): cz.mismatch_class(r)
            for r in report.results
            if not r.boolean_match
        }
        assert (
            failed[("kniha", "knize")]
            is cz.MismatchClass.NO_EVIDENCE_ALTERNATES
        )
        assert (
            failed[("park", "parcích")]
            is cz.MismatchClass.MIXED_STEM_ALTERNATES
        )
        assert ("kluk", "kluci") not in failed

    def test_driver_contract(self, fixture_tree: Path) -> None:
        driver = _driver()
        assert driver.main(["--lang", "cs", "--root", str(fixture_tree)]) == 0
        argv = ["--lang", "cs", "--root", str(fixture_tree)]
        assert driver.main([*argv, "--variant", "spellout"]) == 0
        out = fixture_tree / "analysis" / "czech"
        report = (out / "czech_rule_validation.txt").read_text("utf-8")
        assert report.startswith("Analysis: czech-palatalization (lexical)")
        assert (out / "czech_rule_mismatches_spellout.csv").is_file()
        with (out / "czech_rule_mismatches.csv").open(encoding="utf-8") as fh:
            assert list(csv.DictReader(fh)) == []

    def test_main_writes_the_summary(
        self,
        fixture_tree: Path,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        out = tmp_path / "summary.txt"
        assert cz.main(["--root", str(fixture_tree), "--out", str(out)]) == 0
        text = out.read_text(encoding="utf-8")
        for policy in cz.POLICIES:
            assert f"Policy: {policy}" in text
        assert "Transcription check" in text
        capsys.readouterr()
        assert cz.main(["--root", str(fixture_tree), "--trace", "kluk"]) == 0
        assert "velar-pal" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# Real data
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def real() -> dict[str, list[RowResult]]:
    if not CELLS_CSV.is_file():
        pytest.skip(f"gitignored data absent: {CELLS_CSV}")
    rows = cz.load_rows(PROJECT_ROOT)
    return {
        policy: list(run_validation(cz.build_variant(policy), rows).results)
        for policy in (cz.LEXICAL, cz.ALL_UNDERSPECIFIED)
    }


@pytest.mark.needs_data
@pytest.mark.slow
class TestRealData:
    @pytest.mark.parametrize(
        ("lemma", "form", "sr"),
        [
            ("kluk", "kluci", "klut͡sɪ"),
            ("matka", "matce", "matt͡sɛ"),
            ("kniha", "knize", "kɲɪzɛ"),
            ("vrah", "vrazi", "vrazɪ"),
            ("Čech", "Češi", "t͡ʃɛʃɪ"),
            ("doktor", "doktoři", "doktor̝ɪ"),
            ("biolog", "biologové", "bɪjoloɡovɛː"),
            ("park", "parku", "parku"),
            ("potok", "potoce", "potot͡sɛ"),
            ("moucha", "mouše", "mou̯ʃɛ"),
            ("sestra", "sestře", "sɛstr̝ɛ"),
            ("český", "čeští", "t͡ʃɛʃciː"),
            ("německý", "němečtí", "ɲɛmɛt͡ʃciː"),
            ("bůh", "bože", "boʒɛ"),
        ],
    )
    def test_golden_rows(
        self, real: dict[str, list[RowResult]], lemma: str, form: str, sr: str
    ) -> None:
        found = [
            r
            for r in real[cz.LEXICAL]
            if r.row.lemma == lemma and r.row.plural == form
        ]
        assert found, (lemma, form)
        for result in found:
            assert result.predicted_sr == sr
            assert result.boolean_match and result.ipa_match

    def test_rates(self, real: dict[str, list[RowResult]]) -> None:
        for policy, results in real.items():
            n = len(results)
            assert n > 40_000
            failed = [r for r in results if not r.boolean_match]
            # The single rule's exceptions are the data layer's
            # 'unexplained' forms: far under the tolerance threshold.
            assert len(failed) < n / math.log(n) / 10, policy
            assert all(
                cz.as_czech(r.row).exception_reason == "unexplained"
                for r in failed
            ), policy
