"""The Italian data layer end to end on the fixture tree.

Runs extract, lexicon, qc, attest and features once on
``tests/fixtures/italian`` (plus a synthetic words table) and checks
the golden expectations of the design (§8): amici and medici
palatalize, parchi and laghi keep the velar, chirurghi ~ chirurgi is
variable, amiche is the non-trigger e cell, luce is a frontstem, occhio
is nde:occhio, medico is not a suffix target while pubblico and
elettrico are.
"""

from __future__ import annotations

import csv
import io
import json
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest

from lexpipe import csvio, curated, ids, kaikki, licensing, spec
from lexpipe.langs.italian import PLUGIN, extract, lexicon, stages, store
from lexpipe.langs.italian.layout import Layout
from tests.italian import helpers

pytestmark = pytest.mark.it


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


@pytest.fixture(scope="module")
def tree(tmp_path_factory: pytest.TempPathFactory) -> Iterator[Path]:
    root = helpers.build_tree(tmp_path_factory.mktemp("italian"))
    helpers.run_pipeline(helpers.make_context(root))
    yield root


@pytest.fixture(scope="module")
def layout(tree: Path) -> Layout:
    return Layout(helpers.make_context(tree))


@pytest.fixture(scope="module")
def cells(layout: Layout) -> list[dict[str, str]]:
    return _read(layout.cells)


@pytest.fixture(scope="module")
def lexemes(layout: Layout) -> dict[str, dict[str, str]]:
    return {row["lexeme_id"]: row for row in _read(layout.lexicon_complete)}


def _lexeme(
    lexemes: dict[str, dict[str, str]], word: str, pos: str, ety: str = ""
) -> dict[str, str]:
    found = [
        row
        for row in lexemes.values()
        if row["word"] == word
        and row["pos"] == pos
        and (not ety or row["etymology_number"] == ety)
    ]
    assert len(found) == 1, (word, pos, ety, [r["lexeme_id"] for r in found])
    return found[0]


def _cell(
    cells: list[dict[str, str]], form: str, word: str, pos: str = "noun"
) -> dict[str, str]:
    found = [
        row
        for row in cells
        if row["form"] == form and row["word"] == word and row["pos"] == pos
    ]
    assert found, (form, word, pos)
    return found[0]


# ---------------------------------------------------------------------------
# Plugin
# ---------------------------------------------------------------------------


def test_plugin_registers_the_stages_in_order() -> None:
    plugin = spec.load_plugin("it")
    assert plugin is PLUGIN
    assert list(plugin.stages) == [
        "extract",
        "lexicon",
        "qc",
        "attest",
        "features",
        "derivations",
        "verbs",
        "tp",
        "covariation",
        "spellout",
        "noevidence",
        "logtable",
        "docs",
    ]
    for _, target in stages.STAGE_TABLE:
        assert callable(stages.resolve(target))
    with pytest.raises(ValueError, match="module:function"):
        stages.resolve("lexpipe.langs.italian.extract")


# ---------------------------------------------------------------------------
# Golden expectations (design §8)
# ---------------------------------------------------------------------------


def test_palatalizing_plurals(cells: list[dict[str, str]]) -> None:
    for form, word in (("amici", "amico"), ("medici", "medico")):
        row = _cell(cells, form, word)
        assert row["opportunity"] == "i"
        assert row["mutation"] == "True"
        assert row["exception_reason"] == "undergoer"
        assert row["cell_outcome"] == "palatalized"
        assert row["palatal_output"] == "t͡ʃ"


def test_velar_plurals(cells: list[dict[str, str]]) -> None:
    for form, word, change in (
        ("parchi", "parco", "co→chi"),
        ("laghi", "lago", "go→ghi"),
    ):
        row = _cell(cells, form, word)
        assert row["mutation"] == "False"
        assert row["orth_change"] == change
        assert row["exception_reason"] == "unexplained"
        assert row["cell_outcome"] == "velar"


def test_chirurgo_is_variable(
    cells: list[dict[str, str]], lexemes: dict[str, dict[str, str]]
) -> None:
    velar = _cell(cells, "chirurghi", "chirurgo")
    palatal = _cell(cells, "chirurgi", "chirurgo")
    assert (velar["mutation"], palatal["mutation"]) == ("False", "True")
    assert velar["cell_outcome"] == palatal["cell_outcome"] == "variable"
    row = _lexeme(lexemes, "chirurgo", "noun")
    assert row["outcome_i"] == "variable"
    assert row["mutation_i"] == "variable"
    # Synthetic counts 40 vs 0: one variant has >= 90% (R19).
    assert row["doublet_class"] == "dominant:chirurghi"
    assert row["lemma_suffix"] == "-urgo"
    assert (row["suffixoid"], row["learned_suffixoid"]) == ("-urgo", "True")
    assert row["suffixoid_basis"] == "ending"
    assert (row["semantic_class"], row["class_basis"]) == ("person", "gloss")
    stomaco = _lexeme(lexemes, "stomaco", "noun")
    assert stomaco["doublet_class"] == "true_doublet"  # 3 vs 2 tokens


def test_amiche_is_the_non_trigger_e_cell(
    cells: list[dict[str, str]], lexemes: dict[str, dict[str, str]]
) -> None:
    row = _cell(cells, "amiche", "amica")
    assert row["opportunity"] == "e"
    assert row["mutation"] == "False"
    assert row["orth_change"] == "ca→che"
    assert row["exception_reason"] == "unexplained"
    amica = _lexeme(lexemes, "amica", "noun")
    assert amica["outcome_e"] == "velar" and amica["forms_i"] == ""
    adjective = _cell(cells, "amiche", "amico", pos="adj")
    assert adjective["cell"] == "f.pl" and adjective["opportunity"] == "e"
    # amico and amica are one root (gender pair), as are noun and adj.
    roots = {
        _lexeme(lexemes, "amico", "noun")["root_id"],
        _lexeme(lexemes, "amico", "adj")["root_id"],
        amica["root_id"],
    }
    assert len(roots) == 1


def test_frontstem_and_nde_classes(
    cells: list[dict[str, str]], lexemes: dict[str, dict[str, str]]
) -> None:
    luce = _lexeme(lexemes, "luce", "noun")
    assert luce["frontstem"] == "True" and luce["target"] == ""
    assert luce["exception_reason"] == "non exception"
    assert _cell(cells, "luci", "luce")["cell_outcome"] == "unchanged"
    occhio = _lexeme(lexemes, "occhio", "noun")
    assert occhio["nde_class"] == "occhio"
    assert occhio["exception_reason"] == "nde:occhio"
    assert _cell(cells, "occhi", "occhio")["orth_change"] == "chio→chi"
    psiche = [
        row
        for row in cells
        if row["word"] == "psiche" and row["form"] == "psichi"
    ]
    assert {row["exception_reason"] for row in psiche} == {"nde:psiche"}
    specie = _cell(cells, "specie", "specie")
    assert specie["exception_reason"] == "nde:gimpe_like"
    valigia = _lexeme(lexemes, "valigia", "noun")
    assert valigia["frontstem"] == "True"
    assert valigia["forms_e"] == "valigie|valige"


def test_suffix_targets_by_morphology(
    lexemes: dict[str, dict[str, str]],
) -> None:
    for pos in ("noun", "adj"):
        medico = _lexeme(lexemes, "medico", pos)
        assert medico["target_is_suffix"] == "False"
        assert medico["stratum"] == "inherited"
        pubblico = _lexeme(lexemes, "pubblico", pos)
        assert pubblico["target_is_suffix"] == "True"
        assert pubblico["lemma_suffix"] == "-ico"
    elettrico = _lexeme(lexemes, "elettrico", "adj")
    assert elettrico["target_is_suffix"] == "True"
    assert elettrico["suffix_evidence"] == "latin_icus"
    carico = _lexeme(lexemes, "carico", "noun")
    assert carico["target_is_suffix"] == "False"
    assert carico["deverbal_from_are"] == "True"
    assert carico["deverbal_verb"] == "caricare"
    assert carico["shares_root_with_1st_conj_verb"] == "True"
    medico = _lexeme(lexemes, "medico", "noun")
    # Synchronic only: medicare exists, but medico is not deverbal.
    assert medico["shares_root_with_1st_conj_verb"] == "True"
    assert medico["deverbal_from_are"] == "False"


def test_antico_homographs(lexemes: dict[str, dict[str, str]]) -> None:
    ancient = _lexeme(lexemes, "antico", "adj", "1")
    noun = _lexeme(lexemes, "antico", "noun", "1")
    anterior = _lexeme(lexemes, "antico", "adj", "2")
    stress_shift = _lexeme(lexemes, "antico", "adj", "3")
    assert ancient["root_id"] == noun["root_id"]
    assert len({ancient["root_id"], anterior["root_id"]}) == 2
    assert stress_shift["root_id"] not in (
        ancient["root_id"],
        anterior["root_id"],
    )
    assert ancient["outcome_i"] == "velar"
    assert ancient["velar_post_tonic"] == "True"
    assert stress_shift["outcome_i"] == "palatalized"  # àntico → antici
    assert stress_shift["stress_from_end"] == "3"
    assert stress_shift["velar_post_tonic"] == "False"


def test_labelled_variants_are_filtered(
    cells: list[dict[str, str]], lexemes: dict[str, dict[str, str]]
) -> None:
    magi = _cell(cells, "magi", "mago")
    assert magi["filtered"] == "True"
    assert _lexeme(lexemes, "mago", "noun")["outcome_i"] == "velar"
    focora = _cell(cells, "fòcora", "fuoco")
    assert focora["filtered"] == "True" and focora["opportunity"] == "a"


def test_greco_homographs(
    lexemes: dict[str, dict[str, str]], layout: Layout
) -> None:
    rows = [r for r in lexemes.values() if r["word"] == "greco"]
    assert len(rows) == 3
    assert len({r["root_id"] for r in rows}) == 1
    full = {r["lexeme_id"]: r for r in _read(layout.lexicon_full)}
    statuses = sorted(full[r["lexeme_id"]]["plural_status"] for r in rows)
    assert statuses == ["regular", "regular", "uncountable"]


# ---------------------------------------------------------------------------
# Stage outputs
# ---------------------------------------------------------------------------


def test_lexeme_ids_match_the_core_assignment(layout: Layout) -> None:
    raw = [
        json.loads(line)
        for line in helpers.EN_SAMPLE.read_text("utf-8").splitlines()
    ]
    candidates = [
        e
        for e in raw
        if e.get("pos") in ("noun", "adj") and kaikki.is_lemma_entry(e)
    ]
    expected = set(ids.assign_lexeme_ids("it", candidates))
    _, records = store.read_jsonl(layout.en_lexemes)
    found = [record["lexeme_id"] for record in records]
    assert len(found) == len(set(found))
    assert set(found) == expected


def test_assign_ids_resolves_collisions(tmp_path: Path) -> None:
    def entry(gloss: str) -> dict[str, Any]:
        return {
            "word": "mamma",
            "pos": "noun",
            "lang_code": "it",
            "head_templates": [{"name": "it-noun", "args": {"1": "f"}}],
            "senses": [{"glosses": [gloss]}],
        }

    entries = [entry("mom"), entry("mom"), entry("artichoke"), entry("mom")]
    entries[3]["forms"] = [{"form": "mamme", "tags": ["plural"]}]
    split = tmp_path / "split.jsonl"
    lines = [json.dumps(e, ensure_ascii=False) + "\n" for e in entries]
    split.write_text("".join(lines), encoding="utf-8")
    pending = []
    offset = 0
    for line, item in zip(lines, entries, strict=True):
        base = ids.lexeme_id("it", item)
        pending.append(extract._Pending(base, offset, "{}"))
        offset += len(line.encode("utf-8"))
    assert extract.assign_ids("it", split, pending) == ids.assign_lexeme_ids(
        "it", entries
    )


def test_plural_forms_and_drift(layout: Layout) -> None:
    rows = _read(layout.plural_forms)
    assert set(rows[0]) >= {
        "lexeme_id",
        "form",
        "cell",
        "provenance",
        "qualifiers",
        "raw_arg",
        "tag_error",
    }
    assert {row["provenance"] for row in rows} <= {
        "explicit",
        "module_default",
    }
    magi = [row for row in rows if row["form"] == "magi"]
    assert magi and magi[0]["labels"] == "archaic|dialectal|literary"
    drift = _read(layout.plural_drift)
    nouns = [row for row in drift if row["scope"] == "pure_default_noun"]
    assert nouns and {row["status"] for row in nouns} == {"match"}
    assert {row["module_revision"] for row in drift} == {"79032535"}


def test_lexicon_raw_is_tracked_and_wiktionary_only(layout: Layout) -> None:
    rows = _read(layout.lexicon_raw)
    assert licensing.is_tracked_location(layout.lexicon_raw, layout.ctx.root)
    assert not licensing.is_tracked_location(
        layout.lexicon_full, layout.ctx.root
    )
    assert {row["source"] for row in rows} == {"wiktionary_en"}
    for row in rows:
        licensing.check_redistributable(licensing.split_sources(row["source"]))
    # The full table keeps every column (it.wiktionary ones included).
    full = _read(layout.lexicon_full)
    sources = {row["source"] for row in full}
    assert sources <= {"wiktionary_en", "wiktionary_en|wiktionary_it"}
    assert "wiktionary_en|wiktionary_it" in sources


def test_lexicon_raw_is_lean(layout: Layout) -> None:
    raw = _read(layout.lexicon_raw)
    full = {row["lexeme_id"]: row for row in _read(layout.lexicon_full)}
    assert list(raw[0]) == list(lexicon.RAW_FIELDS)
    assert len(raw) == sum(
        1
        for row in full.values()
        if row["domain"] not in lexicon.RAW_EXCLUDED_DOMAINS
    )
    by_word = {(row["word"], row["pos"]): row for row in raw}
    amico = by_word[("amico", "noun")]
    assert amico["gender"] == full[amico["lexeme_id"]]["gender"]
    assert amico["etymology_summary"] == "inh:la:amīcus"
    assert amico["ipa_raw"] == full[amico["lexeme_id"]]["ipa_raw"]
    # The standard head template keeps only its arguments.
    assert not amico["head_args"].startswith("it-noun:")
    assert kaikki.parse_args_text(csvio.split_list(amico["head_args"])[0])[
        "1"
    ].startswith("m")


def test_attestation_counts(layout: Layout) -> None:
    rows = {
        (row["word"], row["pos"], row["cell"], row["form"]): row
        for row in _read(layout.attestation)
    }
    amico = rows[("amico", "noun", "lemma", "amico")]
    # amico 120 + l'amico 25 (clitic split) + Amico 3 (folded only).
    assert (amico["count_cs"], amico["count_folded"]) == ("145", "148")
    medici = rows[("medico", "noun", "pl", "medici")]
    assert (medici["count_cs"], medici["count_folded"]) == ("150", "170")
    assert medici["form_ambiguous"] == "True"  # Medici, the name
    assert medici["form_ambiguous_cs"] == "False"
    carichi = rows[("carico", "noun", "pl", "carichi")]
    assert carichi["form_ambiguous_cs"] == "True"
    assert "caricare:verb" in carichi["ambiguous_with"]
    stomachi = rows[("stomaco", "noun", "pl", "stomachi")]
    assert stomachi["attested_ge1"] == "True"
    assert stomachi["attested_ge0.5pm"] == "True"
    assert stomachi["attested_ge2pm"] == "False"
    chirurgi = rows[("chirurgo", "noun", "pl", "chirurgi")]
    assert chirurgi["attested_dictionary"] == "True"
    assert chirurgi["attested_ge1"] == "False"
    assert layout.leipzig_freq.is_file()


def test_qc_outputs(layout: Layout) -> None:
    qc_rows = _read(layout.qc)
    assert {row["plural_qc"] for row in qc_rows} <= {
        "confirmed",
        "variant",
        "module_default_only",
        "conflict",
    }
    chirurgo = [row for row in qc_rows if row["word"] == "chirurgo"]
    assert chirurgo[0]["plural_qc"] == "variant"
    assert chirurgo[0]["variant_kind"] == "velar_palatal"
    report = layout.qc_report.read_text("utf-8")
    assert "plural_qc, all cells" in report
    assert "Drift test" in report
    worklist = _read(layout.worklist)
    assert worklist
    for row in worklist:
        assert row["treccani_url"].startswith(
            "https://www.treccani.it/vocabolario/"
        )
        assert row["demauro_url"].startswith(
            "https://dizionario.internazionale.it/parola/"
        )
        licensing.check_redistributable(licensing.split_sources(row["source"]))
    assert "count" not in " ".join(worklist[0])  # ranks only, no counts (R9)


def test_stages_are_skipped_when_up_to_date(
    tree: Path, layout: Layout, capsys: pytest.CaptureFixture[str]
) -> None:
    before = {p: p.read_bytes() for p in layout.ctx.manifest_dir.iterdir()}
    assert len(before) == len(PLUGIN.stages)
    helpers.run_pipeline(helpers.make_context(tree))
    assert capsys.readouterr().out.count("up to date") == len(PLUGIN.stages)
    after = {p: p.read_bytes() for p in layout.ctx.manifest_dir.iterdir()}
    assert after == before


# ---------------------------------------------------------------------------
# Curated decisions (guards)
# ---------------------------------------------------------------------------


def _csv(rows: list[list[str]]) -> str:
    buffer = io.StringIO()
    csv.writer(buffer, lineterminator="\n").writerows(rows)
    return buffer.getvalue()


def _teologo_id(layout: Layout) -> tuple[str, str]:
    for row in _read(layout.lexicon_raw):
        if row["word"] == "teologo":
            return row["lexeme_id"], row["first_gloss"]
    raise AssertionError("teologo missing")


DECISION_HEADER: list[str] = [
    "lexeme_id",
    "learned_suffixoid",
    "semantic_class",
    "decided_by",
    "note",
    "guard_word",
    "guard_pos",
    "guard_first_gloss",
]


def test_curated_decisions_apply(tmp_path: Path, layout: Layout) -> None:
    lexeme_id, gloss = _teologo_id(layout)
    decisions = _csv(
        [
            DECISION_HEADER,
            [lexeme_id, "", "thing", "test", "", "teologo", "noun", gloss],
        ]
    )
    root = helpers.build_tree(
        tmp_path / "t", curated={"suffixoid_classes.csv": decisions}
    )
    helpers.run_pipeline(helpers.make_context(root))
    rows = {
        row["lexeme_id"]: row
        for row in _read(Layout(helpers.make_context(root)).lexicon_complete)
    }
    assert rows[lexeme_id]["semantic_class"] == "thing"
    assert rows[lexeme_id]["class_basis"] == "manual"
    assert rows[lexeme_id]["learned_suffixoid"] == "True"


def test_curated_lookalike_drops_the_suffixoid(
    tmp_path: Path, layout: Layout
) -> None:
    lexeme_id, gloss = _teologo_id(layout)
    decisions = _csv(
        [
            DECISION_HEADER,
            [lexeme_id, "False", "", "test", "", "teologo", "noun", gloss],
        ]
    )
    root = helpers.build_tree(
        tmp_path / "t", curated={"suffixoid_classes.csv": decisions}
    )
    helpers.run_pipeline(helpers.make_context(root))
    row = {
        row["lexeme_id"]: row
        for row in _read(Layout(helpers.make_context(root)).lexicon_complete)
    }[lexeme_id]
    assert (row["suffixoid"], row["learned_suffixoid"]) == ("-logo", "False")
    assert row["suffixoid_basis"] == "manual"
    assert row["lemma_suffix"] == "" and row["suffix_kind"] == ""


def test_curated_guard_mismatch_stops_features(
    tmp_path: Path, layout: Layout
) -> None:
    lexeme_id, _ = _teologo_id(layout)
    decisions = _csv(
        [
            DECISION_HEADER,
            [
                lexeme_id,
                "",
                "thing",
                "test",
                "",
                "teologo",
                "noun",
                "a gloss that changed",
            ],
        ]
    )
    root = helpers.build_tree(
        tmp_path / "t", curated={"suffixoid_classes.csv": decisions}
    )
    with pytest.raises(curated.CuratedGuardError, match="first_gloss"):
        helpers.run_pipeline(helpers.make_context(root))
