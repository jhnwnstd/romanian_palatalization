"""Spell-out TP (R20), no-evidence tables (R21) and the -log table (R27)."""

from __future__ import annotations

import dataclasses
from collections import Counter
from pathlib import Path

import pytest

from lexpipe import tp
from lexpipe.langs.romanian_compare import logtable, noevidence, spellout
from lexpipe.langs.romanian_compare.corpus import load_counts
from lexpipe.langs.romanian_compare.lexicon import read_lexicon
from lexpipe.langs.romanian_compare.observations import Attester
from lexpipe.langs.romanian_compare.records import (
    RootsFile,
    build_records,
    read_verdicts,
)
from tests.romanian_compare.helpers import (
    DECISIONS,
    FIXTURES,
    REAL_FREQ,
    REAL_LEXICON,
    REAL_VERDICTS,
    real_data_or_skip,
)

pytestmark = pytest.mark.ro

# The curated -log decisions of the repository (lemma → decision).
CURATED = logtable.load_decisions(DECISIONS)[1]


def load(lexicon: Path, verdicts: Path, freq: Path) -> RootsFile:
    counts = load_counts(freq)
    records = build_records(
        read_lexicon(lexicon), read_verdicts(verdicts), counts
    )
    return RootsFile(counts.total_tokens, {}, tuple(records))


@pytest.fixture(scope="module")
def fixture_roots() -> RootsFile:
    return load(
        FIXTURES / "lexicon.csv",
        FIXTURES / "verb_paradigm_verdict.csv",
        FIXTURES / "freq.csv",
    )


@pytest.fixture(scope="module")
def real_roots() -> RootsFile:
    real_data_or_skip()
    return load(REAL_LEXICON, REAL_VERDICTS, REAL_FREQ)


# ---------------------------------------------------------------------------
# Spell-out
# ---------------------------------------------------------------------------


def test_spellout_cell_encodes_the_verbalizer_choice(
    fixture_roots: RootsFile,
) -> None:
    attest = Attester(fixture_roots)
    by_lemma = {r.lemma: r for r in fixture_roots.records}
    # urs/urși alternates and lists ursa (-a) and ursi (-i): a doublet.
    urs = spellout.spellout_root(by_lemma["urs"], attest)
    choices = sorted(
        (o.form, o.alternates) for o in urs.triggers[spellout.TRIGGER]
    )
    assert choices == [("ursa", False), ("ursi", True)]
    counts = {o.form: o.count for o in urs.triggers[spellout.TRIGGER]}
    verdict = tp.classify_doublet(
        {"alternating": counts["ursi"], "nonalternating": counts["ursa"]}
    )
    expected = {
        "dominant:alternating": tp.RootStatus.ALTERNATES,
        "dominant:nonalternating": tp.RootStatus.NONALT,
    }.get(verdict, tp.RootStatus.EXCLUDED)
    assert (
        tp.root_status(urs, spellout.TRIGGER, "dictionary", "doublet_excluded")
        is expected
    )
    colac = spellout.spellout_root(by_lemma["colac"], attest)
    assert (
        tp.root_status(colac, spellout.TRIGGER, "dictionary", "dominant")
        is tp.RootStatus.ALTERNATES
    )
    dinte = spellout.spellout_root(by_lemma["dinte"], attest)  # dința (-a)
    assert (
        tp.root_status(dinte, spellout.TRIGGER, "dictionary", "dominant")
        is tp.RootStatus.NONALT
    )
    assert spellout.certain(by_lemma["colac"], "dictionary", attest)
    assert not spellout.certain(by_lemma["atac"], "dictionary", attest)  # -uri
    assert spellout.prespecified(by_lemma["atac"], "dictionary", attest)


def test_spellout_rows_cover_classes_codings_and_frames(
    fixture_roots: RootsFile,
) -> None:
    rows = spellout.spellout_rows(
        fixture_roots.records, Attester(fixture_roots)
    )
    assert {r["table"] for r in rows} == {"spellout"}
    assert {r["target_class"] for r in rows} == set(spellout.CLASSES)
    assert {r["coding"] for r in rows} == set(tp.CODINGS)
    assert {r["V"] for r in rows} == {None, 1000}
    text = spellout.render_report(
        fixture_roots.records, rows, Attester(fixture_roots)
    )
    assert "2. Harmonized spell-out TP" in text


@pytest.mark.needs_data
def test_critique_numbers_are_reproduced(real_roots: RootsFile) -> None:
    full, top = spellout.critique_counts(real_roots.records)
    assert (full.pairs, full.nonpalatalizing, round(full.theta, 1)) == (
        164,
        60,
        32.2,
    )
    assert (
        dict(full.by_suffix)["-i"] == 91 and dict(full.by_suffix)["-a"] == 53
    )
    assert (top.pairs, top.nonpalatalizing, round(top.theta, 1)) == (
        42,
        17,
        11.2,
    )
    assert (full.roots_i_only, full.roots_nonpal_only, full.roots_both) == (
        72,
        43,
        14,
    )


@pytest.mark.needs_data
def test_harmonized_spellout_headline(real_roots: RootsFile) -> None:
    rows = spellout.spellout_rows(real_roots.records, Attester(real_roots))

    def row(**match: object) -> dict[str, object]:
        found = [r for r in rows if all(r[k] == v for k, v in match.items())]
        assert len(found) == 1
        return found[0]

    dorsal = row(
        target_class="dorsal",
        subset="all",
        tau="dictionary",
        coding="doublet_exception",
        rule="palatalize",
        V=None,
    )
    assert dorsal["N_trigger"] == 128 and dorsal["e_nonalt"] == 43
    assert dorsal["tolerated_trigger"] is False


# ---------------------------------------------------------------------------
# No evidence
# ---------------------------------------------------------------------------


def test_base_classes_and_units(fixture_roots: RootsFile) -> None:
    attest = Attester(fixture_roots)
    by_lemma = {r.lemma: r for r in fixture_roots.records}
    assert (
        noevidence.base_class(by_lemma["colac"], "dictionary", attest)
        == "alternating"
    )
    assert (
        noevidence.base_class(by_lemma["atac"], "dictionary", attest)
        == "non-alternating -uri"
    )
    assert (
        noevidence.base_class(by_lemma["albanez"], "dictionary", attest)
        == "non-alternating"
    )
    assert (
        noevidence.base_class(by_lemma["sat"], "dictionary", attest)
        == "no evidence"
    )  # coronal -e
    units = noevidence.units(fixture_roots.records, "dictionary", attest)
    kinds = Counter((u.lemma, u.group, u.outcome) for u in units)
    assert kinds[("colac", "verb -i", "pal")] == 1
    assert kinds[("atac", "back suffix", "avoid")] == 1
    table = noevidence.table(units)
    assert set(table) == set(noevidence.BASE_CLASSES)


def test_noevidence_report(fixture_roots: RootsFile, tmp_path: Path) -> None:
    claims = tmp_path / "check_claims_output.txt"
    claims.write_text(
        "PREDICTION 1: x\nline a\nTOP-1000 y\n", encoding="utf-8"
    )
    text = noevidence.render_report(
        fixture_roots.records, Attester(fixture_roots), claims
    )
    assert (
        "  | line a" in text
        and "TOP-1000" not in text.split("Comparison with")[1]
    )
    missing = noevidence.render_report(
        fixture_roots.records, Attester(fixture_roots), tmp_path / "none.txt"
    )
    assert "not found" in missing


# ---------------------------------------------------------------------------
# -log table
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("lemma", "gloss", "expected"),
    [
        ("psiholog", "psychologist", True),
        ("dialog", "", True),
        ("necrolog", "obituary", True),
        ("chirurg", "surgeon", True),
        ("fitofag", "phytophage, herbivore", True),
        ("bârlog", "cave", False),
        ("amurg", "dusk", False),
        ("milog", "beggar, cadger", False),
        ("x", "anything", False),
    ],
)
def test_learned_suffixoid(lemma: str, gloss: str, expected: bool) -> None:
    assert logtable.learned_suffixoid(lemma, gloss, CURATED)[0] is expected


def test_decisions_come_from_the_curated_file() -> None:
    assert logtable.learned_suffixoid("dialog", "", CURATED) == (
        True,
        "manual",
    )
    assert logtable.learned_suffixoid("dialog", "")[0] is False  # no file
    assert logtable.semantic_class("dialog", "", CURATED) == (
        "thing",
        "manual",
    )
    assert logtable.semantic_class("dialog", "") == ("unknown", "no gloss")
    # The shared rules read 'ideologue' as a person (not as -logue).
    assert logtable.semantic_class("ideolog", "ideologue") == (
        "person",
        "gloss",
    )


@pytest.mark.parametrize(
    ("lemma", "gloss", "expected"),
    [
        ("psiholog", "psychologist", "person"),
        ("dramaturg", "dramaturge, playwright, dramatist", "person"),
        ("dialog", "", "thing"),
        ("catalog", "catalogue", "thing"),
        ("necrolog", "obituary", "thing"),
        ("fitofag", "phytophage, herbivore", "animate"),
        ("ideolog", "ideologue", "person"),
        ("bacteriofag", "bacteriophage", "thing"),
        ("x-log", "", "unknown"),
    ],
)
def test_semantic_class(lemma: str, gloss: str, expected: str) -> None:
    assert logtable.semantic_class(lemma, gloss, CURATED)[0] == expected


def test_log_rows(fixture_roots: RootsFile) -> None:
    rows = logtable.table_rows(
        fixture_roots.records, full=False, decisions=CURATED
    )
    by_lemma = {r["lemma"]: r for r in rows}
    assert sorted(by_lemma) == [
        "amurg",
        "catalog",
        "dialog",
        "monolog",
        "psiholog",
    ]
    assert (
        by_lemma["amurg"]["learned_suffixoid"] is False
    )  # 'dusk': a lookalike
    assert by_lemma["catalog"]["learned_suffixoid"] is True
    assert by_lemma["psiholog"]["semantic_class"] == "person"
    assert by_lemma["psiholog"]["plural_allomorph"] == "-i"
    assert by_lemma["dialog"]["plural_allomorph"] == "-uri"
    assert set(by_lemma["psiholog"]) == set(logtable.LOG_FIELDS)
    assert by_lemma["dialog"]["learned_suffixoid"] is True  # curated
    assert by_lemma["dialog"]["source"] == "wiktionary_ro"
    with pytest.raises(ValueError):
        logtable.table_rows(fixture_roots.records, full=True)
    full = logtable.table_rows(
        fixture_roots.records,
        full=True,
        decisions=CURATED,
        attest=Attester(fixture_roots),
    )
    assert "plural_attested_ge1" in full[0]


def test_log_rows_leave_out_plurals_dex_filled(
    fixture_roots: RootsFile,
) -> None:
    # A "| dex_qc" source means dex_qc_main_csv.py filled the plural from
    # dexonline (terms unverified, not redistributable: R9), so the
    # tracked per-noun table must not carry the row under a Wiktionary
    # source; the untracked full table keeps it.
    records = [
        (
            dataclasses.replace(r, lexicon_source="ro.wiktionary+dex_qc")
            if r.lemma == "psiholog"
            else r
        )
        for r in fixture_roots.records
    ]
    tracked = logtable.table_rows(records, full=False, decisions=CURATED)
    assert "psiholog" not in {r["lemma"] for r in tracked}
    assert all("dex_qc" not in r["lexicon_source"] for r in tracked)
    full = logtable.table_rows(
        records,
        full=True,
        decisions=CURATED,
        attest=Attester(fixture_roots),
    )
    assert "psiholog" in {r["lemma"] for r in full}


@pytest.mark.needs_data
def test_real_log_table(real_roots: RootsFile) -> None:
    full = [
        r
        for r in logtable.table_rows(
            real_roots.records,
            full=True,
            decisions=CURATED,
            attest=Attester(real_roots),
        )
        if r["learned_suffixoid"]
    ]
    kinds = Counter((r["semantic_class"], r["plural_allomorph"]) for r in full)
    assert kinds[("person", "-i")] == 193
    assert (kinds[("thing", "-uri")], kinds[("thing", "-e")]) == (6, 7)
    assert all(
        r["plural_palatalized"] is True
        for r in full
        if r["plural_allomorph"] == "-i"
    )
    # The tracked table leaves out the nouns whose plural DEX filled:
    # dramaturg, subchirurg (person -i) and decalog (thing -uri).
    tracked = [
        r
        for r in logtable.table_rows(
            real_roots.records, full=False, decisions=CURATED
        )
        if r["learned_suffixoid"]
    ]
    kinds = Counter(
        (r["semantic_class"], r["plural_allomorph"]) for r in tracked
    )
    assert kinds[("person", "-i")] == 191
    assert (kinds[("thing", "-uri")], kinds[("thing", "-e")]) == (5, 7)
    assert {r["lemma"] for r in full} - {r["lemma"] for r in tracked} == {
        "dramaturg",
        "subchirurg",
        "decalog",
    }
