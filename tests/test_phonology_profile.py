#!/usr/bin/env python3
"""Tests for the AnalysisProfile hooks (:mod:`phonology.profile`).

Pins three things:

- **Import paths.** ``AnalysisProfile`` lives in ``phonology.profile``
  and is still importable from the Romanian analysis module.
- **Romanian defaults.** Every hook defaults to the Romanian
  implementation, so a profile built without hooks behaves as before.
- **Seams.** A toy non-Romanian analysis, with its own inventory, CSV
  schema, UR builder, cluster buckets, normaliser cascade and
  categoriser, runs the validator and the profile-aware diagnostics
  end to end on synthetic rows. Nothing here reads ``data/``.
"""

from __future__ import annotations

import csv
import dataclasses
import json
import sys
from pathlib import Path
from typing import Mapping

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import phonology.analyses.romanian_palatalization as romanian  # noqa: E402
from phonology.diagnostics.categorise import (  # noqa: E402
    Categorisation,
    RefinedCategory,
    refine_category,
)
from phonology.diagnostics.compare import (  # noqa: E402
    ROMANIAN_NORMALISERS,
    Normalisation,
    NormaliserCascade,
    compare_ipa,
)
from phonology.diagnostics.ordering import (  # noqa: E402
    OrderingStrategy,
    search_orderings,
)
from phonology.diagnostics.perturb import (  # noqa: E402
    PerturbationKind,
    search_perturbations_for,
)
from phonology.g2p import build_ur, build_ur_from_row  # noqa: E402
from phonology.inventory import (  # noqa: E402
    FeatureInventory,
    UnderspecifiedSegment,
)
from phonology.lexicon import (  # noqa: E402
    ClusterTag,
    LexRow,
    _tag_cluster,
    iter_lexicon_rows,
    romanian_cluster_key,
)
from phonology.profile import AnalysisProfile  # noqa: E402
from phonology.rules import UnificationRule  # noqa: E402
from phonology.search import Direction, Search  # noqa: E402
from phonology.segments import ROMANIAN_TOKENIZER, Word  # noqa: E402
from phonology.validation import (  # noqa: E402
    RowResult,
    ValidationReport,
    evaluate_row,
    render_trace,
    run_validation,
    trace_lemmas,
    write_mismatch_csv,
)

# ---------------------------------------------------------------------------
# A toy non-Romanian analysis: /K/ palatalizes before a front vowel
# ---------------------------------------------------------------------------

TOY_FEATURES = (
    "Syllabic",
    "Consonantal",
    "High",
    "Front",
    "CORONAL",
    "DORSAL",
    "Continuant",
)


def _bundle(*values: str) -> dict[str, str]:
    return dict(zip(TOY_FEATURES, values, strict=True))


TOY_SEGMENTS = {
    "p": _bundle("-", "+", "0", "0", "-", "-", "-"),
    "k": _bundle("-", "+", "0", "0", "-", "+", "-"),
    "t͡ʃ": _bundle("-", "+", "0", "0", "+", "-", "-"),
    "l": _bundle("-", "+", "0", "0", "+", "-", "+"),
    "a": _bundle("+", "-", "-", "-", "0", "0", "+"),
    "i": _bundle("+", "-", "+", "+", "0", "0", "+"),
    "u": _bundle("+", "-", "+", "-", "0", "0", "+"),
}

TOY_K = UnderspecifiedSegment(
    label="K",
    base="k",
    clear=frozenset({"CORONAL", "DORSAL"}),
    reason="Toy alternating velar: place left open.",
)
_STOP = {"Consonantal": "+", "Continuant": "-"}
TOY_PAL = UnificationRule(
    name="toy-pal",
    target=_STOP,
    supply={"CORONAL": "+", "DORSAL": "-"},
    search=Search(
        direction=Direction.RIGHT,
        terminator={},
        condition={"Syllabic": "+", "Front": "+"},
    ),
)
TOY_DEFAULT = UnificationRule(
    name="toy-default",
    target=_STOP,
    supply={"CORONAL": "-", "DORSAL": "+"},
)

# The lexicon CSV shares no column names with the Romanian schema.
TOY_CSV = "toy_lexicon.csv"
TOY_COLUMNS = (
    "headword",
    "form",
    "ending",
    "alternates",
    "ipa_headword",
    "ipa_form",
)
TOY_ROWS = (
    ("kak", "kachi", "i", "yes", "kak", "kat͡ʃi"),
    ("pak", "paku", "u", "no", "pak", "paku"),  # back ending: no trigger
    ("lak", "laki", "i", "no", "lak", "laki"),  # lexical exception
    ("kuk", "kuchi", "i", "yes", "kuk", "kut͡ʃiː"),  # contrastive length
    ("zak", "zachi", "i", "yes", "zak", "zat͡ʃi"),  # /z/ not in inventory
)

# Length is contrastive in the toy language: strip stress only.
TOY_NORMALISERS = NormaliserCascade(
    strip=lambda s: s.replace("ˈ", ""),
    steps=(),
)


def toy_row(r: Mapping[str, str]) -> LexRow:
    """Map the toy schema onto the LexRow contract."""
    headword = r["headword"]
    return LexRow(
        lemma=headword,
        plural=r["form"],
        pos="N",
        stem_final=headword[-1],
        opportunity=r["ending"],
        mutation=r["alternates"] == "yes",
        ipa_lemma=r["ipa_headword"],
        ipa_pl=r["ipa_form"],
        exception_reason="",
        cluster=ClusterTag.NONE,
    )


def toy_ur(row: LexRow, inventory: FeatureInventory) -> Word | None:
    """Final /k/ is the alternating /K/; the ending is appended as is."""
    labels = list(row.ipa_lemma) + list(row.opportunity)
    if row.ipa_lemma.endswith("k"):
        labels[len(row.ipa_lemma) - 1] = "K"
    known = set(inventory.base) | set(inventory.underspec)
    if not set(labels) <= known:
        return None
    return tuple(inventory.segment(label) for label in labels)


def toy_cluster_key(row: LexRow) -> str | None:
    return f"{row.stem_final}+{row.opportunity}"


def toy_expected_firings(row: LexRow) -> frozenset[str]:
    return frozenset({TOY_PAL.name}) if row.opportunity == "i" else frozenset()


def toy_fallback_predict(row: LexRow) -> bool:
    return bool(toy_expected_firings(row))


class SpyCategoriser:
    """Record each categoriser-hook call, then defer to the Romanian
    categoriser under the cascade the validator handed over."""

    def __init__(self) -> None:
        self.calls: list[tuple[str, bool, NormaliserCascade]] = []

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
    ) -> Categorisation:
        self.calls.append((row.lemma, ur_built, normalisers))
        return refine_category(
            row,
            predicted_sr,
            attested_field,
            boolean_predicted,
            ur_built,
            normalisers=normalisers,
        )


# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def toy_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    d = tmp_path_factory.mktemp("toy")
    inventory = {"features": list(TOY_FEATURES), "segments": TOY_SEGMENTS}
    (d / "toy_features.json").write_text(
        json.dumps(inventory, ensure_ascii=False), encoding="utf-8"
    )
    with (d / TOY_CSV).open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(TOY_COLUMNS)
        writer.writerows(TOY_ROWS)
    return d


@pytest.fixture
def spy() -> SpyCategoriser:
    return SpyCategoriser()


@pytest.fixture
def toy_profile(toy_dir: Path, spy: SpyCategoriser) -> AnalysisProfile:
    json_path = toy_dir / "toy_features.json"
    return AnalysisProfile(
        name="toy-velar-palatalization",
        inventory=FeatureInventory.load(json_path, (), (TOY_K,)),
        patches=(),
        underspec=(TOY_K,),
        rules=(TOY_PAL, TOY_DEFAULT),
        inventory_json=json_path,
        target_stems=frozenset({"k"}),
        palatalization_rule_names=frozenset({TOY_PAL.name}),
        expected_firings=toy_expected_firings,
        fallback_predict=toy_fallback_predict,
        ur_builder=toy_ur,
        cluster_key=toy_cluster_key,
        normalisers=TOY_NORMALISERS,
        categoriser=spy,
        row_factory=toy_row,
    )


@pytest.fixture
def toy_rows(
    toy_dir: Path, toy_profile: AnalysisProfile
) -> tuple[LexRow, ...]:
    return tuple(
        iter_lexicon_rows(
            toy_dir / TOY_CSV,
            target_stems=toy_profile.target_stems,
            row_factory=toy_profile.row_factory,
        )
    )


@pytest.fixture(scope="module")
def legacy_profile() -> AnalysisProfile:
    """The Romanian profile built the pre-hook way: no hook fields."""
    p = romanian.build_profile()
    return AnalysisProfile(
        name=p.name,
        inventory=p.inventory,
        patches=p.patches,
        underspec=p.underspec,
        rules=p.rules,
        inventory_json=p.inventory_json,
        target_stems=p.target_stems,
        palatalization_rule_names=p.palatalization_rule_names,
        expected_firings=p.expected_firings,
        fallback_predict=p.fallback_predict,
    )


def _romanian_row(lemma: str, ipa: str, stem_final: str, opp: str) -> LexRow:
    return LexRow(
        lemma=lemma,
        plural="",
        pos="N",
        stem_final=stem_final,
        opportunity=opp,
        mutation=False,
        ipa_lemma=ipa,
        ipa_pl="",
        exception_reason="",
        cluster=_tag_cluster(lemma),
    )


def _row(rows: tuple[LexRow, ...], lemma: str) -> LexRow:
    return next(r for r in rows if r.lemma == lemma)


def _result(report: ValidationReport, lemma: str) -> RowResult:
    return next(r for r in report.results if r.row.lemma == lemma)


# ---------------------------------------------------------------------------
# Import paths
# ---------------------------------------------------------------------------


class TestImportPaths:
    def test_romanian_module_reexports_the_profile_class(self) -> None:
        assert romanian.AnalysisProfile is AnalysisProfile
        assert AnalysisProfile.__module__ == "phonology.profile"
        assert "AnalysisProfile" in romanian.__all__

    def test_build_profile_returns_a_profile(self) -> None:
        assert isinstance(romanian.build_profile(), AnalysisProfile)


# ---------------------------------------------------------------------------
# Romanian defaults
# ---------------------------------------------------------------------------


class TestRomanianDefaults:
    def test_hook_defaults_are_the_romanian_functions(
        self, legacy_profile: AnalysisProfile
    ) -> None:
        assert legacy_profile.ur_builder is build_ur_from_row
        assert legacy_profile.cluster_key is romanian_cluster_key
        assert legacy_profile.normalisers is ROMANIAN_NORMALISERS
        assert legacy_profile.categoriser is refine_category
        assert legacy_profile.row_factory == LexRow.from_dict
        assert legacy_profile.tokenizer is ROMANIAN_TOKENIZER

    def test_build_profile_wires_the_same_hooks(
        self, legacy_profile: AnalysisProfile
    ) -> None:
        assert romanian.build_profile() == legacy_profile

    def test_build_ur_from_row_is_build_ur_on_the_row_fields(
        self, legacy_profile: AnalysisProfile
    ) -> None:
        inv = legacy_profile.inventory
        # Final-devoiced IPA: stem_final restores /G/ rather than /K/.
        aerolog = _romanian_row("aerolog", "aɛrɔlɔk", "g", "i")
        # Two variants: the lemma's G2P hint picks [aklamat].
        aclamat = _romanian_row("aclamat", "əkləmat | aklamat", "t", "uri")
        musca = _romanian_row("muscă", "muskə", "s", "e")
        for row in (aerolog, aclamat, musca):
            expected = build_ur(
                row.ipa_lemma,
                row.opportunity,
                inv,
                stem_final=row.stem_final,
                lemma=row.lemma,
            )
            assert expected is not None
            assert build_ur_from_row(row, inv) == expected
        # Both optional arguments are load-bearing for those rows.
        no_stem_final = build_ur(aerolog.ipa_lemma, "i", inv, lemma="aerolog")
        assert build_ur_from_row(aerolog, inv) != no_stem_final
        no_lemma = build_ur(aclamat.ipa_lemma, "uri", inv, stem_final="t")
        assert build_ur_from_row(aclamat, inv) != no_lemma

    def test_romanian_cluster_key_buckets(self) -> None:
        expected = {
            ("muscă", "e"): "sc(ă)+e",
            ("muscă", "i"): "sc(ă)+i",
            ("punct", "e"): "ct(ă)+e",
            ("punct", "i"): "ct(ă)+i",
            ("prost", "i"): "st(ă)+i",
            ("prost", "e"): None,
            ("bumașcă", "i"): None,
            ("rus", "i"): None,
        }
        for (lemma, opp), key in expected.items():
            row = _romanian_row(lemma, "", "", opp)
            assert romanian_cluster_key(row) == key

    def test_default_cascade_is_romanian(self) -> None:
        for pred, att in (("ruʃj", "ruʃi"), ("kat", "kaːt"), ("ab", "xy")):
            assert compare_ipa(pred, att) == compare_ipa(
                pred, att, normalisers=ROMANIAN_NORMALISERS
            )
        # Romanian strips stress and the length mark before comparing.
        assert compare_ipa("kat", "ˈkaːt").matched


# ---------------------------------------------------------------------------
# NormaliserCascade
# ---------------------------------------------------------------------------


class TestNormaliserCascade:
    def test_custom_strip_keeps_length_contrast(self) -> None:
        result = compare_ipa("kat", "kaːt", normalisers=TOY_NORMALISERS)
        assert not result.matched

    def test_custom_step_is_applied_and_labelled(self) -> None:
        final_lax_i = NormaliserCascade(
            strip=lambda s: s,
            steps=(
                (
                    Normalisation.TRAILING_GLIDE,
                    lambda p, a: p[:-1] + "i" if p.endswith("ɪ") else p,
                ),
            ),
        )
        result = compare_ipa("kaɪ", "kai", normalisers=final_lax_i)
        assert result.matched
        assert result.normalisation_applied == (Normalisation.TRAILING_GLIDE,)
        assert not compare_ipa("kaɪ", "kai").matched

    def test_categoriser_compares_under_the_given_cascade(self) -> None:
        row = _romanian_row("kat", "kat", "t", "i")
        default = refine_category(row, "kat", "kaːt", False, ur_built=True)
        toy = refine_category(
            row,
            "kat",
            "kaːt",
            False,
            ur_built=True,
            normalisers=TOY_NORMALISERS,
        )
        assert default.edit_distance == 0.0  # ː stripped: a match
        assert toy.edit_distance == 1.0  # ː kept: one edit away


# ---------------------------------------------------------------------------
# A non-Romanian profile, end to end
# ---------------------------------------------------------------------------


class TestToyProfile:
    def test_row_factory_reads_a_foreign_schema(
        self, toy_dir: Path, toy_rows: tuple[LexRow, ...]
    ) -> None:
        assert [r.lemma for r in toy_rows] == [r[0] for r in TOY_ROWS]
        assert _row(toy_rows, "kuk").ipa_pl == "kut͡ʃiː"
        # The Romanian factory finds none of its columns in this CSV.
        default = iter_lexicon_rows(
            toy_dir / TOY_CSV, target_stems=frozenset({"k"})
        )
        assert list(default) == []

    def test_evaluate_row_uses_the_profile_ur_builder(
        self, toy_profile: AnalysisProfile, toy_rows: tuple[LexRow, ...]
    ) -> None:
        pak = _row(toy_rows, "pak")
        # The Romanian builder has no "-u" plural; the toy one does.
        assert build_ur_from_row(pak, toy_profile.inventory) is None
        res = evaluate_row(pak, toy_profile)
        assert res.ur is not None
        assert res.predicted_sr == "paku"
        assert res.boolean_match and res.ipa_match
        kak = evaluate_row(_row(toy_rows, "kak"), toy_profile)
        assert kak.predicted_sr == "kat͡ʃi"
        assert kak.boolean_predicted and kak.boolean_match and kak.ipa_match

    def test_run_validation_end_to_end(
        self,
        toy_profile: AnalysisProfile,
        toy_rows: tuple[LexRow, ...],
        spy: SpyCategoriser,
    ) -> None:
        report = run_validation(toy_profile, toy_rows)
        assert report.profile_name == "toy-velar-palatalization"
        assert report.total == 5
        assert report.boolean_matches == 4  # lak resists
        assert report.ipa_matches == 3  # lak, and kuk's long vowel
        assert dict(report.by_cluster["k+i"]) == {"match": 3, "mismatch": 1}
        assert dict(report.by_cluster["k+u"]) == {"match": 1}
        assert list(report.by_category) == [RefinedCategory.RULE_FAILURE]
        failures = report.by_category[RefinedCategory.RULE_FAILURE]
        assert [r.row.lemma for r in failures] == ["lak", "kuk"]
        zak = _result(report, "zak")
        assert zak.ur is None and zak.boolean_predicted  # fallback_predict
        # The categoriser hook saw every row, under the profile cascade.
        assert spy.calls == [
            ("kak", True, TOY_NORMALISERS),
            ("pak", True, TOY_NORMALISERS),
            ("lak", True, TOY_NORMALISERS),
            ("kuk", True, TOY_NORMALISERS),
            ("zak", False, TOY_NORMALISERS),
        ]

    def test_normalisers_reach_the_match_flag_and_distance(
        self, toy_profile: AnalysisProfile, toy_rows: tuple[LexRow, ...]
    ) -> None:
        kuk = _row(toy_rows, "kuk")
        toy = evaluate_row(kuk, toy_profile)
        assert not toy.ipa_match
        assert toy.distance is not None and toy.distance.ipa_edit == 1.0
        as_romanian = dataclasses.replace(
            toy_profile, normalisers=ROMANIAN_NORMALISERS
        )
        res = evaluate_row(kuk, as_romanian)
        assert res.ipa_match
        assert res.distance is not None and res.distance.ipa_edit == 0.0

    def test_mismatch_csv_uses_the_cluster_key(
        self,
        toy_profile: AnalysisProfile,
        toy_rows: tuple[LexRow, ...],
        tmp_path: Path,
    ) -> None:
        report = run_validation(toy_profile, toy_rows)
        out = tmp_path / "mismatches.csv"
        assert write_mismatch_csv(report, out) == 2
        with out.open(encoding="utf-8", newline="") as f:
            written = {r["lemma"]: r for r in csv.DictReader(f)}
        assert set(written) == {"lak", "kuk"}
        assert {r["cluster"] for r in written.values()} == {"k+i"}
        assert written["kuk"]["ur"] == "kuKi"

    def test_trace_lemmas_reads_with_the_row_factory(
        self, toy_dir: Path, toy_profile: AnalysisProfile
    ) -> None:
        traces = trace_lemmas(
            ["pak", "kak", "nope"],
            toy_profile,
            toy_dir / TOY_CSV,
            include_explanation=True,
            include_ordering=OrderingStrategy.ADJACENT_SWAP,
            include_perturbation=(PerturbationKind.CLEAR_REMOVE,),
        )
        assert [t.found for t in traces] == [True, True, False]
        pak = traces[0]
        assert pak.explanation is not None
        assert pak.ordering is not None and pak.ordering.baseline_matched
        # Rebuilt per perturbation with the toy ur_builder.
        assert pak.perturbation is not None
        assert pak.perturbation.baseline_matched
        text = render_trace(pak, toy_profile.inventory.base_segments())
        assert "= /paku/" in text

    def test_diagnostics_judge_success_under_the_profile_cascade(
        self, toy_profile: AnalysisProfile, toy_rows: tuple[LexRow, ...]
    ) -> None:
        kuk = _row(toy_rows, "kuk")
        ur = toy_profile.ur_builder(kuk, toy_profile.inventory)
        assert ur is not None
        as_romanian = dataclasses.replace(
            toy_profile, normalisers=ROMANIAN_NORMALISERS
        )
        for profile, matched in ((toy_profile, False), (as_romanian, True)):
            ordering = search_orderings(profile, ur, kuk.ipa_pl)
            assert ordering.baseline_matched is matched
            perturbation = search_perturbations_for(
                profile, kuk, kinds=(PerturbationKind.CLEAR_REMOVE,)
            )
            assert perturbation.baseline_matched is matched


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
