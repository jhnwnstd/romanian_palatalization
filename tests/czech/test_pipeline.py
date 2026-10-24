"""All five stages on the fixture tree: golden expectations (design §8).

Runs extract, lexicon, qc, attest and features once on the CC BY-SA
fixtures plus synthetic MorfFlex, Leipzig and CNK tables
(:mod:`tests.czech.treebuilder`).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lexpipe import csvio, curated, kaikki, licensing
from lexpipe.langs.czech import adjudication, features, layout, lexicon, qc
from tests.czech import treebuilder

pytestmark = pytest.mark.cs


def lexeme(
    rows: list[dict[str, str]], lemma: str, pos: str = "noun", **match: str
) -> dict[str, str]:
    found = [
        r
        for r in rows
        if r["lemma"] == lemma
        and r["pos"] == pos
        and all(r[k] == v for k, v in match.items())
    ]
    assert len(found) == 1, (lemma, pos, match, len(found))
    return found[0]


def cells_of(
    rows: list[dict[str, str]], lemma: str, cell: str, **match: str
) -> dict[str, dict[str, str]]:
    return {
        r["form"]: r
        for r in rows
        if r["lemma"] == lemma
        and r["cell"] == cell
        and r["primary"] == "True"
        and all(r[k] == v for k, v in match.items())
    }


# -- lexicon -----------------------------------------------------------------


def test_lexeme_ids_unique_and_sorted(
    lexicon_rows: list[dict[str, str]],
) -> None:
    ids = [r["lexeme_id"] for r in lexicon_rows]
    assert len(ids) == len(set(ids))
    assert ids == sorted(ids)


def test_headers(paths: layout.CzechPaths) -> None:
    def header(path: Path) -> tuple[str, ...]:
        return tuple(
            path.read_text(encoding="utf-8").split("\n", 1)[0].split(",")
        )

    assert header(paths.lexicon_raw) == lexicon.RAW_FIELDS
    assert header(paths.lexicon_full) == lexicon.LEXICON_FIELDS
    assert header(paths.cells_raw) == lexicon.CELL_FIELDS
    assert header(paths.cells) == features.CELL_FIELDS
    assert header(paths.lexicon_complete) == features.COMPLETE_FIELDS
    assert header(paths.worklist) == qc.WORKLIST_FIELDS
    assert header(paths.adjudications) == adjudication.FIELDS


def test_names_policy(lexicon_rows: list[dict[str, str]]) -> None:
    noun = lexeme(lexicon_rows, "Čech")
    name = lexeme(lexicon_rows, "Čech", "name")
    assert (noun["domain"], noun["capitalized"]) == ("in", "True")
    assert (name["domain"], name["ood_class"]) == ("out", "name:surname")
    assert noun["root_id"] == name["root_id"]
    inka = sorted(r["ood_class"] for r in lexicon_rows if r["lemma"] == "Inka")
    assert inka == ["", "name:given_name", "name:given_name"]


def test_special_types_out_of_domain(
    lexicon_rows: list[dict[str, str]],
) -> None:
    ruka = lexeme(lexicon_rows, "ruka")
    assert (ruka["domain"], ruka["ood_class"]) == ("out", "dual_remnant")
    kure = lexeme(lexicon_rows, "kuře")
    assert (kure["domain"], kure["ood_class"]) == ("out", "tstem")


def test_plech_has_no_paradigm(
    lexicon_rows: list[dict[str, str]], paths: layout.CzechPaths
) -> None:
    plech = lexeme(lexicon_rows, "plech")
    assert plech["paradigm_source"] == "none"
    assert plech["n_cells"] == "0" and int(plech["n_cells_second"]) > 0
    raw = [
        r
        for r in csvio.read_rows(paths.cells_raw)
        if r["lexeme_id"] == plech["lexeme_id"]
    ]
    assert raw and all(
        r["primary"] == "False" and r["source"] == "wiktionary_cs" for r in raw
    )


def test_sach_is_two_lexemes(lexicon_rows: list[dict[str, str]]) -> None:
    check = lexeme(lexicon_rows, "šach", readings="m.in")
    shah = lexeme(lexicon_rows, "šach", readings="m.an")
    assert check["lexeme_id"] != shah["lexeme_id"]
    assert check["root_id"] != shah["root_id"]
    assert (check["cell_LOC.PL"], shah["cell_LOC.PL"]) == ("šachách", "šaších")


def test_head_codes_and_specs(lexicon_rows: list[dict[str, str]]) -> None:
    agent = lexeme(lexicon_rows, "agent", head_code="m-an,m-in")
    assert agent["animacy_ambiguous"] == "True"
    assert agent["table_animacy"] == "an|in"
    potok = lexeme(lexicon_rows, "potok")
    assert potok["decl_spec"] == "m.gena.locu:e"
    assert potok["explicit_cells"] == "GEN.SG|LOC.SG"
    assert potok["cell_LOC.SG"] == "potoku|potoce"  # unranked, source order
    biolog = lexeme(lexicon_rows, "biolog")
    assert biolog["head_feminine"] == "bioložka"
    assert (biolog["lemma_suffix"], biolog["suffix_evidence"]) == (
        "log",
        "etymology",
    )


def test_form_source(paths: layout.CzechPaths) -> None:
    rows = csvio.read_rows(paths.cells_raw)
    by_key = {
        (r["lexeme_id"].split(":")[1], r["cell"], r["form"], r["source"]): r
        for r in rows
    }
    assert (
        by_key[("potok", "LOC.SG", "potoce", "wiktionary_en")]["form_source"]
        == "explicit"
    )
    assert (
        by_key[("kluk", "NOM.PL", "kluci", "wiktionary_en")]["form_source"]
        == "module_default"
    )
    assert (
        by_key[("kluk", "NOM.PL", "klucí", "wiktionary_cs")]["form_source"]
        == "explicit"
    )
    kluku = by_key[("kluk", "LOC.SG", "kluku", "wiktionary_en")]
    assert kluku["syncretic_cells"] == "DAT.SG|VOC.SG"


def test_tracked_lexicon_is_wiktionary_only(
    lexicon_rows: list[dict[str, str]],
    raw_rows: list[dict[str, str]],
    paths: layout.CzechPaths,
    fixture_root: Path,
) -> None:
    for rows in (lexicon_rows, raw_rows):
        sources = {s for r in rows for s in csvio.split_list(r["source"])}
        assert sources <= {"wiktionary_en", "wiktionary_cs"}
        licensing.check_redistributable(sorted(sources))
    for path in (paths.lexicon_raw, paths.worklist, paths.qc_report):
        assert licensing.is_tracked_location(path, fixture_root)
    assert not licensing.is_tracked_location(paths.lexicon_full, fixture_root)


def test_tracked_raw_lexicon_is_lean(
    lexicon_rows: list[dict[str, str]], raw_rows: list[dict[str, str]]
) -> None:
    full = {r["lexeme_id"]: r for r in lexicon_rows}
    kept = [
        r["lexeme_id"]
        for r in lexicon_rows
        if not r["ood_class"].startswith("name")
        and r["ood_class"] != "multiword"
    ]
    assert [r["lexeme_id"] for r in raw_rows] == kept
    kluk = lexeme(raw_rows, "kluk")
    assert kluk["cell_NOM.PL"] == full[kluk["lexeme_id"]]["cell_NOM.PL"]
    assert (kluk["gender"], kluk["animacy"]) == ("m", "an")
    # The head code is the first head argument (cs-noun left unnamed).
    head = csvio.split_list(kluk["head_args"])
    assert (
        kaikki.parse_args_text(head[0])["1"]
        == full[kluk["lexeme_id"]]["head_code"]
    )
    assert lexicon.tracked.__doc__  # names and multiword stay untracked
    assert not any(
        r["lemma"] == "Čech" and r["pos"] == "name" for r in raw_rows
    )


# -- qc -------------------------------------------------------------------


def test_qc_join_and_statuses(paths: layout.CzechPaths) -> None:
    rows = csvio.read_rows(paths.qc_csv)
    status = {(r["lemma"], r["cell"], r["form"]): r["status"] for r in rows}
    assert status[("kluk", "NOM.PL", "kluci")] == qc.PRESENT
    assert status[("rok", "LOC.SG", "roce")] == qc.NONSTANDARD_ONLY
    assert status[("liga", "DAT.SG", "lize")] == qc.STYLED_LEMMA
    assert status[("biolog", "LOC.PL", "biolozích")] == qc.ABSENT
    # plech has no primary paradigm, so it is not compared
    assert not any(lemma == "plech" for lemma, _, _ in status)
    # jazyk m-in joins MorfFlex "jazyk" (I), never "jazyk-2" (M)
    jazyk = {r["morfflex_lemmas"] for r in rows if r["lemma"] == "jazyk"}
    assert jazyk == {"jazyk"}
    sach = {r["morfflex_lemmas"] for r in rows if r["lemma"] == "šach"}
    assert sach == {"šach", "šach-2"}
    cech = {r["name_labels"] for r in rows if r["lemma"] == "Čech"}
    assert cech == {"E|Y"}
    # the negated MorfFlex row never appears
    assert "nekluci" not in {r["form"] for r in rows}


def test_worklist_and_report_hold_no_morfflex_content(
    paths: layout.CzechPaths,
) -> None:
    worklist = csvio.read_rows(paths.worklist)
    keys = {(r["lemma"], r["cell"]) for r in worklist}
    assert ("biolog", "LOC.PL") in keys
    assert all(r["source"].startswith("wiktionary_") for r in worklist)
    assert all(r["ijp_url"].startswith(qc.IJP_URL) for r in worklist)
    morfflex_only = {
        r["form"]
        for r in csvio.read_rows(paths.qc_csv)
        if r["origin"] == "morfflex"
    }
    listed = {
        f for r in worklist for f in csvio.split_list(r["wiktionary_forms"])
    }
    assert not morfflex_only & listed
    report = paths.qc_report.read_text(encoding="utf-8")
    assert "klucí" not in report and "kluci" not in report


# -- attest --------------------------------------------------------------------


def test_attestation_counts(paths: layout.CzechPaths) -> None:
    rows = csvio.read_rows(paths.attestation)
    row = {
        (r["lexeme_id"].split(":")[1], r["cell"], r["form"], r["source"]): r
        for r in rows
    }
    kluci = row[("kluk", "NOM.PL", "kluci", "wiktionary_en")]
    assert (kluci["count_cs"], kluci["count_cap"], kluci["count"]) == (
        "870",
        "30",
        "900",
    )
    cesi = row[("Čech", "NOM.PL", "Češi", "wiktionary_en")]
    assert (cesi["count_cs"], cesi["count"]) == ("2000", "2000")
    sachach = row[("šach", "LOC.PL", "šachách", "wiktionary_en")]
    assert (sachach["attested_ge1"], sachach["attested_ge2pm"]) == (
        "True",
        "False",
    )
    biolozich = row[("biolog", "LOC.PL", "biolozích", "wiktionary_en")]
    assert biolozich["count"] == "0" and biolozich["attested_ge1"] == "False"


def test_form_ambiguity(paths: layout.CzechPaths) -> None:
    rows = csvio.read_rows(paths.attestation)
    flags = {
        (r["form"], r["cell"]): (r["form_ambiguous"], r["ambiguous_with"])
        for r in rows
        if r["source"] == "wiktionary_en"
    }
    assert flags[("šachové", "NOM.PL")] == ("True", "w:adj:šachový")
    assert flags[("šachem", "INS.SG")][0] == "True"  # both šach lexemes
    assert flags[("kluci", "NOM.PL")] == ("False", "")
    # the surname shares the root but is a name: ambiguous
    assert flags[("Čechové", "NOM.PL")][0] == "True"


def test_cnk_cross_check(paths: layout.CzechPaths) -> None:
    rows = {r["lemma"]: r for r in csvio.read_rows(paths.lexeme_frequency)}
    assert rows["rok"]["cnk_lemma_freq"] == "358211"
    assert rows["kluk"]["lemma_count"] == "346"
    report = paths.attestation_report.read_text(encoding="utf-8")
    assert "Spearman" in report


# -- features (golden expectations) ----------------------------------------------


def test_kluci_undergoer(cell_rows: list[dict[str, str]]) -> None:
    kluci = cells_of(cell_rows, "kluk", "NOM.PL")["kluci"]
    assert kluci["exception_label"] == "undergoer"
    assert (kluci["orth_change"], kluci["palatal_output"]) == ("k→c", "t͡s")


def test_biolog_mixed_then_resister(
    complete_rows: list[dict[str, str]], cell_rows: list[dict[str, str]]
) -> None:
    biolog = lexeme(complete_rows, "biolog")
    assert biolog["summary_E0_dictionary"] == "mixed"
    assert biolog["summary_E0_ge1"] == "resister"
    assert (
        biolog["p_cells_E0_dictionary"],
        biolog["n_cells_E0_dictionary"],
    ) == (
        "LOC.PL",
        "NOM.PL",
    )
    assert biolog["target_is_suffix"] == "True"
    nompl = cells_of(cell_rows, "biolog", "NOM.PL")["biologové"]
    assert nompl["exception_label"] == "non exception (allomorph)"


def test_plech_and_ruka(complete_rows: list[dict[str, str]]) -> None:
    plech = lexeme(complete_rows, "plech")
    assert (plech["paradigm_source"], plech["summary_E0_dictionary"]) == (
        "none",
        "no_evidence",
    )
    assert lexeme(complete_rows, "ruka")["summary_E0_dictionary"] == "out"


def test_studenti_coronal_inferred(cell_rows: list[dict[str, str]]) -> None:
    studenti = cells_of(cell_rows, "student", "NOM.PL")["studenti"]
    assert studenti["coronal_inferred"] == "True"
    assert (studenti["mutation_kind"], studenti["exception_label"]) == (
        "inferred",
        "undergoer",
    )
    studente = cells_of(cell_rows, "student", "VOC.SG")["studente"]
    assert (studente["allomorph_class"], studente["n_kind"]) == (
        "N",
        "nontrigger",
    )


def test_potok_locative_doublet(
    cell_rows: list[dict[str, str]], complete_rows: list[dict[str, str]]
) -> None:
    loc = cells_of(cell_rows, "potok", "LOC.SG")
    assert set(loc) == {"potoce", "potoku"}
    assert loc["potoce"]["exception_label"] == "undergoer"
    assert loc["potoku"]["exception_label"] == "non exception (allomorph)"
    assert loc["potoku"]["cell_confounded"] == "True"  # also DAT/VOC.SG
    assert loc["potoce"]["cell_confounded"] == "False"
    assert lexeme(complete_rows, "potok")["summary_E0_ge1"] == "mixed"


def test_sach_lexemes_differ_in_locative_plural(
    cell_rows: list[dict[str, str]], complete_rows: list[dict[str, str]]
) -> None:
    labels = {
        r["form"]: r["exception_label"]
        for r in cell_rows
        if r["lemma"] == "šach"
        and r["cell"] == "LOC.PL"
        and r["primary"] == "True"
    }
    assert labels == {
        "šaších": "undergoer",
        "šachách": "non exception (allomorph)",
    }
    summaries = sorted(
        r["summary_E0_dictionary"]
        for r in complete_rows
        if r["lemma"] == "šach"
    )
    assert summaries == ["mixed", "resister"]


def test_adjective_and_vocatives(cell_rows: list[dict[str, str]]) -> None:
    cesti = cells_of(cell_rows, "český", "NOM.PL.M.AN")["čeští"]
    assert (cesti["orth_change"], cesti["exception_label"]) == (
        "sk→šť",
        "undergoer",
    )
    otce = cells_of(cell_rows, "otec", "VOC.SG")["otče"]
    assert (otce["target_class"], otce["exception_label"]) == (
        "c",
        "undergoer",
    )
    bratre = cells_of(cell_rows, "bratr", "VOC.SG")["bratře"]
    assert (bratre["cluster"], bratre["exception_label"]) == (
        "Cr",
        "undergoer",
    )


def test_invariants(cell_rows: list[dict[str, str]]) -> None:
    coronal_p = [
        r
        for r in cell_rows
        if r["target_class"] == "coronal_inferred"
        and r["allomorph_class"] == "P"
        and r["primary"] == "True"
    ]
    assert coronal_p and all(
        r["exception_label"] == "undergoer" for r in coronal_p
    )
    velar_p = [
        r
        for r in cell_rows
        if r["target_class"] == "velar"
        and r["allomorph_class"] == "P"
        and r["primary"] == "True"
    ]
    assert velar_p and all(r["mutation"] == "True" for r in velar_p)


# -- manifests and curated decisions ---------------------------------------------


def test_rerun_is_skipped(
    fixture_root: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    before = {
        stage: [p.read_bytes() for p in outs]
        for stage, outs in treebuilder.stage_outputs(fixture_root).items()
    }
    treebuilder.run_stages(fixture_root)
    assert capsys.readouterr().out.count("up to date") >= len(
        treebuilder.STAGES
    )
    after = {
        stage: [p.read_bytes() for p in outs]
        for stage, outs in treebuilder.stage_outputs(fixture_root).items()
    }
    assert before == after


def _adjudicate(
    root: Path, guard_forms: str, confirmed: str = "potoce"
) -> None:
    paths = layout.paths_for(treebuilder.context(root))
    potok = lexeme(csvio.read_rows(paths.lexicon_full), "potok")
    row = {
        "lexeme_id": potok["lexeme_id"],
        "cell": "LOC.SG",
        "confirmed": confirmed,
        "citation": "test",
        "guard_lemma": "potok",
        "guard_gender": "m",
        "guard_animacy": "in",
        "guard_forms": guard_forms,
    }
    # One row per decided cell (composite key lexeme_id x cell).
    csvio.write_rows(paths.adjudications, [row], adjudication.FIELDS)


def test_curated_decision_applies(fixture_root: Path, tmp_path: Path) -> None:
    root = treebuilder.copy_fixture_tree(fixture_root, tmp_path / "tree")
    _adjudicate(root, "potoku|potoce")
    treebuilder.run_stages(root, ["features"])
    paths = layout.paths_for(treebuilder.context(root))
    loc = cells_of(csvio.read_rows(paths.cells), "potok", "LOC.SG")
    assert loc["potoce"]["curated_status"] == "confirmed"
    assert loc["potoku"]["curated_status"] == "rejected"
    potok = lexeme(csvio.read_rows(paths.lexicon_complete), "potok")
    assert potok["summary_E0_dictionary"] == "undergoer"


def test_stale_guard_stops_the_run(fixture_root: Path, tmp_path: Path) -> None:
    root = treebuilder.copy_fixture_tree(fixture_root, tmp_path / "tree")
    _adjudicate(root, "potoku|potoce", confirmed="potocem")
    with pytest.raises(ValueError, match="not among the guarded forms"):
        treebuilder.run_stages(root, ["features"])
    # the curator saw other forms than the data now holds
    _adjudicate(root, "potoce|potoku")
    with pytest.raises(curated.CuratedGuardError):
        treebuilder.run_stages(root, ["features"])
