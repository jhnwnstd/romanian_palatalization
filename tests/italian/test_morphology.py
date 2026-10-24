"""Stem-final targets and plural alternations (design §5)."""

from __future__ import annotations

import pytest

from lexpipe.langs.italian import keys, morphology

pytestmark = pytest.mark.it


def _analyze(
    lemma: str, form: str
) -> tuple[morphology.StemAnalysis, morphology.FormAnalysis]:
    lemma_key, form_key = keys.match_key(lemma), keys.match_key(form)
    stem = morphology.analyze_lemma(lemma_key)
    return stem, morphology.analyze_form(lemma_key, stem, form_key)


@pytest.mark.parametrize(
    ("lemma", "form", "opportunity", "mutation", "change", "output"),
    [
        ("amico", "amici", "i", True, "co→ci", "t͡ʃ"),
        ("medico", "medici", "i", True, "co→ci", "t͡ʃ"),
        ("belga", "belgi", "i", True, "ga→gi", "d͡ʒ"),
        ("osco", "osci", "i", True, "sco→sci", "ʃ"),
        ("parco", "parchi", "i", False, "co→chi", ""),
        ("lago", "laghi", "i", False, "go→ghi", ""),
        ("bosco", "boschi", "i", False, "sco→schi", ""),
        ("sacco", "sacchi", "i", False, "cco→cchi", ""),
        ("patriarca", "patriarchi", "i", False, "ca→chi", ""),
        ("amica", "amiche", "e", False, "ca→che", ""),
        ("capocuoco", "capicuochi", "i", False, "co→chi", ""),
        ("gatto", "gatti", "i", False, "to→ti", ""),
    ],
)
def test_alternations(
    lemma: str,
    form: str,
    opportunity: str,
    mutation: bool,
    change: str,
    output: str,
) -> None:
    _, analysis = _analyze(lemma, form)
    assert analysis.opportunity == opportunity
    assert analysis.mutation is mutation
    assert analysis.orth_change == change
    assert analysis.palatal_output == output
    assert analysis.aligned


def test_stem_classes() -> None:
    occhio, occhi = _analyze("occhio", "occhi")
    assert (occhio.nde_class, occhio.cluster) == ("occhio", "chi")
    assert occhi.mutation is False and occhi.orth_change == "chio→chi"
    psiche, psichi = _analyze("psiche", "psichi")
    assert psiche.nde_class == "psiche" and psichi.mutation is False
    luce, luci = _analyze("luce", "luci")
    assert luce.frontstem and luce.target == "" and luce.stem_final == "c"
    assert luci.mutation is False and luci.orth_change == "ce→ci"
    valigia, valige = _analyze("valigia", "valige")
    assert valigia.frontstem and valige.orth_change == "gia→ge"
    parco, _ = _analyze("parco", "parchi")
    assert (parco.stem_final, parco.cluster, parco.target) == (
        "c",
        "rc",
        "dorsal",
    )
    gatto, _ = _analyze("gatto", "gatti")
    assert (gatto.stem_final, gatto.target) == ("t", "coronal")
    assert morphology.analyze_lemma("città") == morphology.NO_STEM
    assert morphology.analyze_lemma("abete rosso") == morphology.NO_STEM


def test_no_front_opportunity() -> None:
    _, uova = _analyze("uovo", "uova")
    assert uova.opportunity == "a" and uova.mutation is None
    _, focora = _analyze("fuoco", "fòcora")
    assert focora.opportunity == "a" and focora.mutation is None
    _, uomini = _analyze("uomo", "uomini")
    assert uomini.opportunity == "none" and not uomini.aligned
    _, specie = _analyze("specie", "specie")
    assert specie.opportunity == "none"
    assert morphology.gimpe_like("specie", "specie")
    assert not morphology.gimpe_like("luce", "luci")


def test_root_internal_count() -> None:
    chiesa = morphology.analyze_lemma("chiesa")
    assert morphology.root_internal_velar_front(chiesa, "chiesa") == 1
    occhio = morphology.analyze_lemma("occhio")
    assert morphology.root_internal_velar_front(occhio, "occhio") == 0
