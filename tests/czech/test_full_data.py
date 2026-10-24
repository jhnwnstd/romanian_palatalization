"""Checks on the real Czech outputs under data/czech/ (gitignored).

Skipped when the pipeline has not been run
(``python scripts/lang_pipeline.py --lang cs``).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lexpipe import csvio, spec
from lexpipe.context import RunContext
from lexpipe.langs.czech import jsonio, layout
from tests.czech import treebuilder

pytestmark = [pytest.mark.cs, pytest.mark.needs_data]

PATHS = layout.paths_for(
    RunContext(spec.LANGUAGES["cs"], root=treebuilder.PROJECT_ROOT)
)


def _require(path: Path) -> None:
    if not path.is_file():
        pytest.skip(f"gitignored output absent: {path}")


def test_lexeme_ids_unique_over_the_full_lexicon() -> None:
    _require(PATHS.lexicon_full)
    ids = [row["lexeme_id"] for row in jsonio.iter_csv(PATHS.lexicon_full)]
    assert len(ids) > 40_000
    assert len(ids) == len(set(ids))
    assert ids == sorted(ids)


def test_tracked_raw_lexicon_is_lean() -> None:
    _require(PATHS.lexicon_raw)
    assert PATHS.lexicon_raw.stat().st_size < 10_000_000
    ids = [row["lexeme_id"] for row in jsonio.iter_csv(PATHS.lexicon_raw)]
    assert len(ids) > 40_000 and ids == sorted(ids)


def test_cell_rows_unique() -> None:
    _require(PATHS.cells_raw)
    keys = [
        (r["lexeme_id"], r["cell"], r["form"], r["source"], r["table"])
        for r in jsonio.iter_csv(PATHS.cells_raw)
    ]
    assert len(keys) == len(set(keys))


def test_golden_expectations_on_real_data() -> None:
    _require(PATHS.lexicon_complete)
    rows = {
        (r["lemma"], r["pos"], r["readings"]): r
        for r in jsonio.iter_csv(PATHS.lexicon_complete)
        if r["lemma"] in {"kluk", "biolog", "plech", "ruka", "šach", "student"}
    }
    assert (
        rows[("kluk", "noun", "m.an")]["summary_E0_dictionary"] == "undergoer"
    )
    biolog = rows[("biolog", "noun", "m.an")]
    assert biolog["summary_E0_dictionary"] == "mixed"
    assert biolog["summary_E0_ge1"] == "resister"
    assert rows[("plech", "noun", "m.in")]["paradigm_source"] == "none"
    assert rows[("ruka", "noun", "f")]["ood_class"] == "dual_remnant"
    assert (
        rows[("šach", "noun", "m.in")]["lexeme_id"]
        != rows[("šach", "noun", "m.an")]["lexeme_id"]
    )
    assert rows[("student", "noun", "m.an")]["coronal_inferred"] == "True"


def test_tracked_outputs_pass_the_licence_gate() -> None:
    _require(PATHS.lexicon_raw)
    for path in (PATHS.lexicon_raw, PATHS.worklist):
        sources = {
            s
            for r in csvio.read_rows(path)
            for s in csvio.split_list(r["source"])
        }
        assert sources <= {"wiktionary_en", "wiktionary_cs"}
