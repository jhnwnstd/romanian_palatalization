"""Unit tests: outcome-blind attributes, orthographic readings, cells."""

from __future__ import annotations

from dataclasses import replace

import pytest

from lexpipe import tp
from lexpipe.langs.romanian_compare import lexicon as lx
from lexpipe.langs.romanian_compare import morphology as mo
from lexpipe.langs.romanian_compare import observations as ob
from lexpipe.langs.romanian_compare.records import (
    Candidate,
    Derivative,
    Legacy,
    RootRecord,
    VerbCell,
)

pytestmark = pytest.mark.ro

BASE_ROW = lx.LexiconRow(
    lemma="x",
    gloss="",
    pos="N",
    gender="NEUT",
    stem_final="",
    cluster="",
    frontstem_dorsal=False,
    plural="",
    mutation=False,
    orth_change="",
    opportunity="none",
    derived_verbs=(),
    deriv_suffixes=(),
    derived_adj=(),
    ipa_derived_adj="",
    etym_lang="",
    exception_reason="",
    tp_in_domain=False,
    target_is_suffix=False,
    lemma_suffix="",
    source="https://ro.wiktionary.org/wiki/x",
    notes="",
    freq=0,
)


def row(**fields: object) -> lx.LexiconRow:
    return replace(BASE_ROW, **fields)  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# lexicon: outcome-blind attributes
# ---------------------------------------------------------------------------


def test_fold_key_folds_cedillas_only() -> None:
    assert lx.fold_key(" Urşi ") == "urși"
    assert lx.fold_key("ţară") == "țară"
    assert lx.fold_key("Ăla") == "ăla"  # no accent or breve folding


@pytest.mark.parametrize(
    ("fields", "expected"),
    [
        ({"lemma": "genunchi", "stem_final": "c", "cluster": "chi"}, "ochi"),
        ({"lemma": "ureche", "stem_final": "c", "cluster": "che"}, "paduchi"),
        ({"lemma": "cruce", "frontstem_dorsal": True}, "frontstem"),
        ({"lemma": "anestezie", "stem_final": "z"}, "coronal_frontstem"),
        ({"lemma": "acadea", "stem_final": "d"}, "coronal_frontstem"),
        ({"lemma": "dinte", "stem_final": "t"}, ""),  # inflectional -e
        ({"lemma": "colac", "stem_final": "c"}, ""),
    ],
)
def test_nde_class(fields: dict[str, object], expected: str) -> None:
    assert lx.nde_class(row(**fields)) == expected


def test_nde_class_ignores_the_outcome() -> None:
    plain = row(lemma="ureche", stem_final="c", cluster="che", plural="urechi")
    assert lx.nde_class(plain) == lx.nde_class(replace(plain, mutation=True))


@pytest.mark.parametrize(
    ("lemma", "seg", "expected"),
    [
        ("aluniște", "t", True),  # t after ș: protected
        ("țintă", "t", True),  # ț already in the stem tail
        ("pânză", "d", False),
        ("tată", "t", False),
        ("urs", "s", False),
    ],
)
def test_palatal_context(lemma: str, seg: str, expected: bool) -> None:
    assert lx.palatal_context(row(lemma=lemma, stem_final=seg)) is expected


@pytest.mark.parametrize(
    ("lemma", "plural", "expected"),
    [
        ("fată", "fete", False),  # vowel alternation, legacy "suppletion"
        ("carte", "cărți", False),  # palatal reflex counts as the same letter
        ("tată", "tați", False),
        ("oaste", "oști", True),
        ("oseminte", "oase", True),
        ("x", "", False),
    ],
)
def test_suppletive_stem(lemma: str, plural: str, expected: bool) -> None:
    assert lx.suppletive_stem(lemma, plural) is expected


def test_ica_suppletion_and_flags() -> None:
    assert lx.ica_suppletion("alunică", "alunele")
    assert not lx.ica_suppletion("rândunică", "rândunici")
    flagged = row(
        lemma="alunică",
        stem_final="c",
        plural="alunele",
        target_is_suffix=True,
    )
    assert lx.flags(flagged) == ("ica_suppletion", "suffix_target")


@pytest.mark.parametrize(
    ("code", "expected"),
    [
        ("", "unknown"),
        ("la", "latin_unresolved"),
        ("it", "borrowed"),
        ("FR", "borrowed"),
    ],
)
def test_stratum(code: str, expected: str) -> None:
    assert lx.stratum(code) == expected


@pytest.mark.parametrize(
    ("fields", "nde", "expected"),
    [
        ({"plural": "colaci", "opportunity": "i"}, "", "-i"),
        ({"plural": "abace", "opportunity": "e"}, "", "-e"),
        ({"plural": "focuri", "opportunity": "uri"}, "", "-uri"),
        ({"plural": "esperanto", "opportunity": "none"}, "", "unparsed"),
        ({"plural": "", "opportunity": "none"}, "", "none"),
        ({"plural": "păduchi", "opportunity": "none"}, "paduchi", "-i"),
        ({"plural": "mănunchiuri", "opportunity": "none"}, "ochi", "-i"),
    ],
)
def test_plural_allomorph(
    fields: dict[str, object], nde: str, expected: str
) -> None:
    assert lx.plural_allomorph(row(**fields), nde) == expected


def test_group_nouns_merges_by_priority() -> None:
    without = row(
        lemma="esperanto", source="https://en.wiktionary.org/wiki/Esperanto"
    )
    with_plural = row(
        lemma="esperanto",
        plural="esperanto",
        source="https://ro.wiktionary.org/wiki/esperanto",
    )
    adjective = row(lemma="esperanto", pos="ADJ", plural="zzz")
    for order in (
        [without, with_plural, adjective],
        [adjective, with_plural, without],
    ):
        groups = lx.group_nouns(order)
        assert list(groups) == ["esperanto"]
        assert groups["esperanto"] == [with_plural, without]


def test_lexicon_source_labels() -> None:
    assert (
        lx.lexicon_source("https://ro.wiktionary.org/wiki/a | dex_qc")
        == "ro.wiktionary+dex_qc"
    )
    assert (
        lx.lexicon_source("https://en.wiktionary.org/wiki/A")
        == "en.wiktionary"
    )


# ---------------------------------------------------------------------------
# morphology: orthographic readings
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("root_part", "seg", "cluster", "vowel", "expected"),
    [
        ("atac", "c", "", "i", True),
        ("atac", "c", "", "back", False),
        ("bloch", "c", "", "e", False),
        ("covrig", "g", "", "i", True),
        ("plăț", "t", "", "i", True),
        ("plăt", "t", "", "i", False),
        ("braz", "d", "", "i", True),
        ("brad", "d", "", "i", False),
        ("urș", "s", "", "i", True),
        ("proșt", "s", "st", "i", True),
        ("prost", "s", "st", "i", False),
        ("mușt", "s", "scă", "e", True),
        ("musc", "s", "scă", "e", False),
        ("obraj", "z", "", "i", True),
        ("albanez", "z", "", "i", False),
        ("puncț", "t", "ct", "e", True),
        ("picur", "c", "", "a", None),
    ],
)
def test_reflex(
    root_part: str, seg: str, cluster: str, vowel: str, expected: bool | None
) -> None:
    assert mo.reflex(root_part, seg, cluster, vowel) is expected


@pytest.mark.parametrize(
    ("form", "endings", "seg", "expected"),
    [
        ("ataci", mo.VERB_ENDINGS_2SG, "c", ("i", True, False)),
        ("atace", mo.VERB_ENDINGS_SUBJ3SG, "c", ("e", True, False)),
        ("blochezi", mo.VERB_ENDINGS_2SG, "c", ("e", False, True)),
        ("mijlocești", mo.VERB_ENDINGS_2SG, "c", ("e", True, True)),
        ("legendezi", mo.VERB_ENDINGS_2SG, "d", ("e", False, True)),
        ("îmbogățești", mo.VERB_ENDINGS_2SG, "t", ("e", True, True)),
        ("riști", mo.VERB_ENDINGS_2SG, "s", ("i", True, False)),
    ],
)
def test_read_verb_cell(
    form: str,
    endings: tuple[tuple[str, bool], ...],
    seg: str,
    expected: tuple[str, bool, bool],
) -> None:
    cluster = "sc" if form == "riști" else ""
    reading = mo.read_verb_cell(form, endings, seg, cluster)
    assert reading is not None
    assert (reading.vowel, reading.palatal, reading.augment) == expected


def test_unreadable_cells() -> None:
    assert mo.read_verb_cell("picuri", mo.VERB_ENDINGS_2SG, "c", "") is None
    assert mo.first_variant("&mdash;", "x") == ""
    assert (
        mo.first_variant("îmbogățești|mbogățești", "îmbogăți") == "îmbogățești"
    )
    assert mo.first_variant("nfeudezi|înfeudezi", "înfeuda") == "înfeudezi"


@pytest.mark.parametrize(
    ("verb", "suffix", "seg", "cluster", "expected"),
    [
        ("încolăci", "-i", "c", "", ("i", True)),
        ("ataca", "-a", "c", "", ("back", False)),
        ("târgui", "-ui", "g", "", ("back", False)),
        ("plăti", "-i", "t", "", ("i", False)),
        ("înverzi", "-i", "d", "", ("i", True)),
        ("catalogiza", "-iza", "g", "", ("i", True)),
        ("prosti", "-i", "s", "st", ("i", False)),
    ],
)
def test_read_derived_verb(
    verb: str, suffix: str, seg: str, cluster: str, expected: tuple[str, bool]
) -> None:
    reading = mo.read_derived_verb(verb, suffix, seg, cluster)
    assert reading is not None and (reading.vowel, reading.palatal) == expected


@pytest.mark.parametrize(
    ("adjective", "seg", "expected"),
    [
        ("catolicesc", "c", ("e", True, "-esc")),
        ("focos", "c", ("back", False, "-os")),
        ("turistic", "s", ("i", False, "-ic")),
        ("verzui", "d", ("back", True, "-ui")),
    ],
)
def test_read_derived_adjective(
    adjective: str, seg: str, expected: tuple[str, bool, str]
) -> None:
    cluster = "st" if adjective == "turistic" else ""
    reading = mo.read_derived_adjective(adjective, seg, cluster)
    assert reading is not None
    assert (reading.vowel, reading.palatal, reading.suffix) == expected


def test_participles_and_candidates() -> None:
    assert mo.participle("ataca", "-a") == "atacat"
    assert mo.participle("târgui", "-ui") == "târguit"
    assert mo.participle("buluci", "-i") == "bulucit"
    assert mo.candidate_verbs("foc", "c", "") == ("foci", "fochi")
    assert mo.candidate_verbs("vad", "d", "") == ("vazi", "vadi")
    assert mo.candidate_verbs("prost", "s", "st") == ("proști", "prosti")
    assert mo.candidate_verbs("punct", "t", "ct") == ("puncți", "puncti")
    assert mo.candidate_verbs("tată", "t", "") == ("tați", "tati")


@pytest.mark.parametrize(
    ("trigger", "vowel", "expected"),
    [
        ("-i", "i", "P"),
        ("-i", "e", "N"),
        ("-i", "back", "nonfront"),
        ("-e", "i", "N"),
        ("-i/-e", "e", "P"),
        ("-i/-e", None, None),
    ],
)
def test_allomorph_for(
    trigger: str, vowel: str | None, expected: str | None
) -> None:
    assert mo.allomorph_for(trigger, vowel) == expected


def test_audit_like_verdict() -> None:
    reading = mo.Reading
    assert (
        mo.audit_like_verdict("c", [reading("ataci", "i", True)])
        == "ALTERNATES"
    )
    assert (
        mo.audit_like_verdict("c", [reading("blochezi", "e", False, True)])
        == "NO-ALTERNATION"
    )
    assert (
        mo.audit_like_verdict("t", [reading("plătești", "e", False, True)])
        == "TRIGGER-ABSENT"
    )
    assert mo.audit_like_verdict("c", []) == "NO-READING"


# ---------------------------------------------------------------------------
# observations: realizations become cells
# ---------------------------------------------------------------------------


def record(**fields: object) -> RootRecord:
    base = RootRecord(
        root_id="ro:x",
        lemma="x",
        gloss="",
        gender="NEUT",
        seg="c",
        cluster="",
        target_class="dorsal",
        nde="",
        flags=(),
        stratum="unknown",
        lexicon_source="ro.wiktionary",
        plural="",
        plural_allomorph="none",
        plural_palatal=None,
        lemma_count=5,
        plural_count=0,
        paradigm_count=5,
        verb="",
        audit_verdict="",
        verb_cells=(),
        derivatives=(),
        candidates=(),
        legacy=Legacy(False, "", "none", None, "", False, ""),
    )
    return replace(base, **fields)  # type: ignore[arg-type]


ALL_TAUS = {"dictionary": True, "ge1": True, "ge0.5pm": True, "ge2pm": True}


def attest(key: str) -> dict[str, bool]:
    return dict(ALL_TAUS)


def status(rec: RootRecord, evidence: str, trigger: str) -> tp.RootStatus:
    root = ob.root_observation(rec, evidence, attest)
    return tp.root_status(root, trigger, "dictionary", "dominant")


def test_plural_cells() -> None:
    e_plural = record(
        plural="abace", plural_allomorph="-e", plural_palatal=True
    )
    assert status(e_plural, "E0", "-e") is tp.RootStatus.ALTERNATES
    assert status(e_plural, "E0", "-i/-e") is tp.RootStatus.ALTERNATES
    assert (
        status(e_plural, "E0", "-i") is tp.RootStatus.SHIELDED
    )  # N allomorph
    uri = record(plural="focuri", plural_allomorph="-uri")
    assert status(uri, "E0", "-i/-e") is tp.RootStatus.SHIELDED
    none = record(plural="esperanto", plural_allomorph="unparsed")
    assert status(none, "E0", "-i") is tp.RootStatus.NO_EVIDENCE


def test_evidence_sets_add_cells() -> None:
    rec = record(
        plural="atacuri",
        plural_allomorph="-uri",
        verb="ataca",
        verb_cells=(VerbCell("ataca", "2sg", "ataci", "i", True, False, 3),),
        derivatives=(
            Derivative("verb", "ataca", "-a", "back", False, 9, 12),
            Derivative("adj", "atacic", "-ic", "i", False, 1, 1),
        ),
        candidates=(Candidate("ataci", True, 3),),
    )
    assert status(rec, "E0", "-i") is tp.RootStatus.SHIELDED
    assert status(rec, "E1", "-i") is tp.RootStatus.ALTERNATES
    # E2 adds a non-palatalizing -ic derivative: a doublet, coded by
    # corpus share (3 vs 1 tokens: not dominant, the more frequent wins).
    assert status(rec, "E2", "-i") is tp.RootStatus.ALTERNATES
    root = ob.root_observation(rec, "E2", attest)
    assert (
        tp.root_status(root, "-i", "dictionary", "doublet_exception")
        is tp.RootStatus.VARIABLE
    )
    with pytest.raises(ValueError):
        ob.realizations(rec, "E3")
