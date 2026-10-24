"""Stage ``qc``: automated QC of the Wiktionary cells against MorfFlex.

MorfFlex CZ 2.1 plays the role the DEX plays for Romanian (design §3),
with one caveat (R11): it is a generated full-form dictionary, so
agreement with it is quality control, never evidence that a form is
used.

Steps
-----
1. **Projection** (sub-stage ``qc-morfflex``, its own manifest): one
   streaming pass over the MorfFlex xz file
   (:func:`lexpipe.langs.czech.morfflex.build_db`) keeps affirmative
   NN rows and the hard-adjective cell AAMP1----1A---- whose lemma base
   is an in-domain Wiktionary noun or adjective, or whose form is in
   the Leipzig vocabulary (R7). It is rebuilt only when MorfFlex, the
   lemma list or the vocabulary changes (``--force qc`` re-runs the
   comparison, not the pass; delete ``nc/czech_morfflex.sqlite`` to
   force it).
2. **Join** (R14) on (lemma, gender, animacy): a lexeme's MorfFlex
   lemmas are the raw lemmas whose base equals the Wiktionary lemma
   and whose gender letters overlap the lexeme's readings
   (:func:`lexpipe.langs.czech.heads.morfflex_gender_matches`): jazyk
   m-in meets ``jazyk`` (I), not ``jazyk-2`` (M). Adjectives match on
   the base. Proper-name labels are recorded (Čech_;E_;Y).
3. **Agreement** per Wiktionary cell form: ``present`` (a standard
   MorfFlex form of the cell, Var ``-``/1–4, of an unmarked lemma),
   ``styled_lemma`` (a standard form, but of a style-labelled lemma:
   ``_,h`` colloquial Brňák, ``_,s`` spelling variant Afgánec),
   ``nonstandard_only`` (only Var 5–9), ``absent``, ``cell_missing``
   (MorfFlex has the lemma but nothing in that cell) or ``no_lemma``.
   The cell's reference set is its standard forms of unmarked lemmas,
   or of styled lemmas when there are none; reference forms missing
   from Wiktionary are listed as ``morfflex_only``.

Outputs
-------
``data/czech/nc/czech_morfflex_qc.csv`` (NC: copies MorfFlex forms)
    One row per compared form: lexemes by id, cells in canonical
    order, Wiktionary forms in source order, then MorfFlex-only forms.
``analysis/czech/czech_qc_report.txt`` (tracked; aggregates only)
``analysis/czech/czech_manual_worklist.csv`` (tracked)
    In-domain lexemes whose stem ends in a target (velar, r, inferred
    coronal) and whose Wiktionary forms in a cell with a P/N choice for
    their gender (:func:`lexpipe.langs.czech.ortho.opportunity_cells`)
    are absent or non-standard in MorfFlex. Wiktionary forms only, no MorfFlex
    content, with the IJP URL for a **manual** lookup; decisions go to
    ``curated/czech/cell_adjudications.csv``.
"""

from __future__ import annotations

import contextlib
import sqlite3
import urllib.parse
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

from lexpipe import csvio
from lexpipe.context import RunContext, run_stage, sha256_file
from lexpipe.langs.czech import (
    adjudication,
    cells,
    heads,
    jsonio,
    leipzig,
    morfflex,
    ortho,
)
from lexpipe.langs.czech.layout import CzechPaths, paths_for

STAGE: Final[str] = "qc"
DB_STAGE: Final[str] = "qc-morfflex"
# 2: reads czech_lexicon_full.csv; licence-gated report (NC aggregates).
VERSION: Final[str] = "2"
DB_VERSION: Final[str] = "1"
IJP_URL: Final[str] = "https://prirucka.ujc.cas.cz/?slovo="
NC_LICENCE: Final[str] = "by-nc-sa"
WIKT_LICENCE: Final[str] = "by-sa"
SOURCE_MORFFLEX: Final[str] = "morfflex"
# The tracked QC report holds MorfFlex agreement counts: aggregates of
# a CC BY-NC-SA source (R9; DATA_LICENSE.md draft policy).
REPORT_SOURCES: Final[tuple[str, ...]] = (
    "wiktionary_en",
    "wiktionary_cs",
    SOURCE_MORFFLEX,
)

PRESENT: Final[str] = "present"
NONSTANDARD_ONLY: Final[str] = "nonstandard_only"
STYLED_LEMMA: Final[str] = "styled_lemma"
ABSENT: Final[str] = "absent"
CELL_MISSING: Final[str] = "cell_missing"
NO_LEMMA: Final[str] = "no_lemma"
MORFFLEX_ONLY: Final[str] = "morfflex_only"
DISAGREEING: Final[frozenset[str]] = frozenset({ABSENT, NONSTANDARD_ONLY})
WORKLIST_TARGETS: Final[frozenset[str]] = frozenset(
    {ortho.VELAR, ortho.R_CLASS, ortho.CORONAL}
)

QC_FIELDS: Final[tuple[str, ...]] = (
    "lexeme_id",
    "lemma",
    "pos",
    "readings",
    "cell",
    "form",
    "origin",
    "form_source",
    "status",
    "cell_agreement",
    "morfflex_var",
    "morfflex_lemmas",
    "name_labels",
    "source",
    "licence_class",
)
WORKLIST_FIELDS: Final[tuple[str, ...]] = (
    "lexeme_id",
    "lemma",
    "pos",
    "gender",
    "animacy",
    "target_class",
    "cell",
    "wiktionary_forms",
    "form_source",
    "ijp_url",
    "source",
    "licence_class",
)


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------


@dataclass
class QcLexeme:
    """An in-domain lexeme with its primary cells, in Wiktionary order."""

    lexeme_id: str
    lemma: str
    pos: str
    readings: tuple[str, ...]
    gender: tuple[str, ...]
    animacy: tuple[str, ...]
    flags: tuple[str, ...]
    source: str
    forms: dict[str, list[str]] = field(default_factory=dict)
    form_sources: dict[tuple[str, str], str] = field(default_factory=dict)


def load_lexemes(paths: CzechPaths) -> dict[str, QcLexeme]:
    """In-domain nouns and adjectives with a primary paradigm."""
    lexemes: dict[str, QcLexeme] = {}
    for row in jsonio.iter_csv(paths.lexicon_full):
        if row["domain"] != "in" or row["pos"] not in ("noun", "adj"):
            continue
        if row["paradigm_source"] not in ("en", "cs"):
            continue
        lexemes[row["lexeme_id"]] = QcLexeme(
            lexeme_id=row["lexeme_id"],
            lemma=row["lemma"],
            pos=row["pos"],
            readings=tuple(csvio.split_list(row["readings"])),
            gender=tuple(csvio.split_list(row["gender"])),
            animacy=tuple(csvio.split_list(row["animacy"])),
            flags=tuple(csvio.split_list(row["spec_flags"])),
            source="wiktionary_" + row["paradigm_source"],
        )
    for row in jsonio.iter_csv(paths.cells_raw):
        lexeme = lexemes.get(row["lexeme_id"])
        if lexeme is None or csvio.parse_bool(row["primary"]) is not True:
            continue
        forms = lexeme.forms.setdefault(row["cell"], [])
        if row["form"] not in forms:
            forms.append(row["form"])
            lexeme.form_sources[(row["cell"], row["form"])] = row[
                "form_source"
            ]
    return lexemes


def wanted_genders(readings: Sequence[str]) -> frozenset[str]:
    """MorfFlex gender letters of a lexeme's readings (any if unknown)."""
    letters: set[str] = set()
    for label in readings:
        gender, _, animacy = label.partition(".")
        letters |= heads.morfflex_genders(heads.GenderReading(gender, animacy))
    return frozenset(letters) or frozenset("MIFN")


# ---------------------------------------------------------------------------
# MorfFlex side
# ---------------------------------------------------------------------------


@dataclass
class MorfflexCells:
    """The MorfFlex forms of one lexeme, by cell."""

    lemmas: list[str] = field(default_factory=list)
    name_labels: set[str] = field(default_factory=set)
    standard: dict[str, list[str]] = field(default_factory=dict)
    styled: dict[str, list[str]] = field(default_factory=dict)
    nonstandard: dict[str, list[str]] = field(default_factory=dict)
    variants: dict[tuple[str, str], set[str]] = field(default_factory=dict)

    @property
    def found(self) -> bool:
        return bool(self.lemmas)

    def reference(self, cell: str) -> list[str]:
        """Standard forms of unmarked lemmas, else of styled ones."""
        return self.standard.get(cell) or self.styled.get(cell, [])

    def empty(self, cell: str) -> bool:
        return not (
            self.standard.get(cell)
            or self.styled.get(cell)
            or self.nonstandard.get(cell)
        )


def morfflex_cells(con: sqlite3.Connection, lexeme: QcLexeme) -> MorfflexCells:
    """MorfFlex forms of the lemmas joined to *lexeme* (R14)."""
    rows = con.execute(
        "SELECT raw_lemma, name_labels, style_labels, tag, cell, gender, "
        "var, var_class, form FROM morfflex WHERE base = ? AND in_lexicon = 1 "
        "ORDER BY raw_lemma, tag, form",
        (lexeme.lemma,),
    ).fetchall()
    adjective = lexeme.pos == "adj"
    wanted = wanted_genders(lexeme.readings)
    by_lemma: dict[str, list[tuple[Any, ...]]] = defaultdict(list)
    for row in rows:
        is_adjective_row = row[3] == morfflex.ADJECTIVE_TAG
        if adjective != is_adjective_row:
            continue
        by_lemma[row[0]].append(row)
    found = MorfflexCells()
    for raw_lemma in sorted(by_lemma):
        lemma_rows = by_lemma[raw_lemma]
        if not adjective and not any(
            heads.morfflex_gender_matches(r[5], wanted) for r in lemma_rows
        ):
            continue
        found.lemmas.append(raw_lemma)
        for (
            _raw,
            names,
            styles,
            _tag,
            cell,
            gender,
            var,
            klass,
            form,
        ) in lemma_rows:
            if cell is None or klass == morfflex.ABBREVIATION:
                continue
            if not adjective and not heads.morfflex_gender_matches(
                gender, wanted
            ):
                continue
            found.name_labels.update(csvio.split_list(names))
            if klass != morfflex.STANDARD:
                bucket = found.nonstandard
            else:
                bucket = found.styled if styles else found.standard
            forms = bucket.setdefault(cell, [])
            if form not in forms:
                forms.append(form)
            found.variants.setdefault((cell, form), set()).add(var)
    return found


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------


def form_status(form: str, cell: str, found: MorfflexCells) -> str:
    """Agreement of one Wiktionary form (see the module docstring)."""
    if not found.found:
        return NO_LEMMA
    if found.empty(cell):
        return CELL_MISSING
    if form in found.standard.get(cell, []):
        return PRESENT
    if form in found.styled.get(cell, []):
        return STYLED_LEMMA
    if form in found.nonstandard.get(cell, []):
        return NONSTANDARD_ONLY
    return ABSENT


def cell_agreement(
    wiktionary: Sequence[str], cell: str, found: MorfflexCells
) -> str:
    """Set relation of a cell's Wiktionary and standard MorfFlex forms."""
    if not found.found:
        return NO_LEMMA
    if found.empty(cell):
        return CELL_MISSING
    standard = set(found.reference(cell))
    ours = set(wiktionary)
    if ours == standard:
        return "agree"
    if ours < standard:
        return "wikt_subset"
    if ours > standard:
        return "wikt_superset"
    return "overlap" if ours & standard else "disjoint"


def compared_cells(lexeme: QcLexeme) -> list[str]:
    if lexeme.pos == "adj":
        return [c for c in lexeme.forms if c == cells.HARD_ADJECTIVE_CELL]
    return sorted(lexeme.forms, key=cells.cell_sort_key)


def qc_rows(
    lexeme: QcLexeme, found: MorfflexCells
) -> Iterator[dict[str, Any]]:
    """Rows of the NC QC table for one lexeme."""
    common = {
        "lexeme_id": lexeme.lexeme_id,
        "lemma": lexeme.lemma,
        "pos": lexeme.pos,
        "readings": list(lexeme.readings),
        "morfflex_lemmas": found.lemmas,
        "name_labels": sorted(found.name_labels),
    }
    for cell in compared_cells(lexeme):
        wiktionary = lexeme.forms[cell]
        agreement = cell_agreement(wiktionary, cell, found)
        for form in wiktionary:
            yield {
                **common,
                "cell": cell,
                "form": form,
                "origin": "wiktionary",
                "form_source": lexeme.form_sources.get((cell, form), ""),
                "status": form_status(form, cell, found),
                "cell_agreement": agreement,
                "morfflex_var": sorted(found.variants.get((cell, form), ())),
                "source": lexeme.source,
                "licence_class": WIKT_LICENCE,
            }
        for form in found.reference(cell):
            if form in wiktionary:
                continue
            yield {
                **common,
                "cell": cell,
                "form": form,
                "origin": "morfflex",
                "status": MORFFLEX_ONLY,
                "cell_agreement": agreement,
                "morfflex_var": sorted(found.variants.get((cell, form), ())),
                "source": SOURCE_MORFFLEX,
                "licence_class": NC_LICENCE,
            }


def worklist_rows(
    lexeme: QcLexeme, statuses: Mapping[tuple[str, str], str]
) -> Iterator[dict[str, Any]]:
    """Worklist rows of one lexeme: Wiktionary forms only."""
    analysis = ortho.analyse_stem(
        lexeme.pos, lexeme.lemma, lexeme.forms, lexeme.flags
    )
    target = analysis.target.target_class
    if target not in WORKLIST_TARGETS:
        return
    choices = ortho.opportunity_cells(lexeme.pos, lexeme.readings)
    for cell in adjudication.CELLS:
        forms = lexeme.forms.get(cell)
        if not forms or cell not in choices:
            continue
        if not any(
            statuses.get((cell, form)) in DISAGREEING for form in forms
        ):
            continue
        yield {
            "lexeme_id": lexeme.lexeme_id,
            "lemma": lexeme.lemma,
            "pos": lexeme.pos,
            "gender": list(lexeme.gender),
            "animacy": list(lexeme.animacy),
            "target_class": target,
            "cell": cell,
            "wiktionary_forms": forms,
            "form_source": [
                lexeme.form_sources.get((cell, f), "") for f in forms
            ],
            "ijp_url": IJP_URL + urllib.parse.quote(lexeme.lemma, safe=""),
            "source": lexeme.source,
            "licence_class": WIKT_LICENCE,
        }


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


@dataclass
class Tally:
    """Aggregates for the tracked report (counts only)."""

    lexemes: Counter[str] = field(default_factory=Counter)
    joined: Counter[str] = field(default_factory=Counter)
    name_labelled: Counter[str] = field(default_factory=Counter)
    statuses: Counter[tuple[str, str]] = field(default_factory=Counter)
    agreements: Counter[tuple[str, str]] = field(default_factory=Counter)
    by_source: Counter[tuple[str, str]] = field(default_factory=Counter)
    morfflex_only: Counter[str] = field(default_factory=Counter)
    worklist: Counter[str] = field(default_factory=Counter)

    def add(
        self,
        lexeme: QcLexeme,
        found: MorfflexCells,
        rows: Iterable[Mapping[str, Any]],
    ) -> None:
        self.lexemes[lexeme.pos] += 1
        if found.found:
            self.joined[lexeme.pos] += 1
        if found.name_labels:
            self.name_labelled[lexeme.pos] += 1
        seen_cells: set[str] = set()
        for row in rows:
            cell = str(row["cell"])
            if row["origin"] == "morfflex":
                self.morfflex_only[cell] += 1
                continue
            self.statuses[(cell, str(row["status"]))] += 1
            self.by_source[(str(row["form_source"]), str(row["status"]))] += 1
            if cell not in seen_cells:
                seen_cells.add(cell)
                self.agreements[(cell, str(row["cell_agreement"]))] += 1


def _share(part: int, whole: int) -> str:
    return (
        f"{part:>8,} ({100 * part / whole:5.1f}%)" if whole else f"{part:>8,}"
    )


def render_report(tally: Tally, db_meta: Mapping[str, str]) -> str:
    """Plain-text report of aggregates (no forms, no corpus counts)."""
    lines = [
        "Czech QC: Wiktionary cells against MorfFlex CZ 2.1",
        "=" * 52,
        "",
        "MorfFlex is a generated dictionary: agreement is QC, not",
        "attestation (design R11). Aggregates only (R9).",
        "",
        "MorfFlex projection (data/czech/nc/czech_morfflex.sqlite)",
        "  rows matched by the NN/AAMP1 filter: "
        f"{db_meta.get('matched', '?')}",
        "  rows stored (lexicon lemmas or Leipzig forms): "
        f"{db_meta.get('rows', '?')}",
        f"  MorfFlex sha256: {db_meta.get('source_sha256', '?')}",
        "",
        "Lexemes compared (in-domain, with a primary paradigm)",
    ]
    for pos in sorted(tally.lexemes):
        total = tally.lexemes[pos]
        lines.append(
            f"  {pos:<5} {total:>7,}   joined to a MorfFlex lemma: "
            f"{_share(tally.joined[pos], total)}   with name labels: "
            f"{tally.name_labelled[pos]:,}"
        )
    lines += ["", "Wiktionary forms by cell and status"]
    statuses = (
        PRESENT,
        STYLED_LEMMA,
        NONSTANDARD_ONLY,
        ABSENT,
        CELL_MISSING,
        NO_LEMMA,
    )
    lines.append(
        f"  {'cell':<12}"
        + "".join(f"{s:>18}" for s in statuses)
        + f"{'morfflex_only':>15}"
    )
    for cell in sorted(
        {c for c, _ in tally.statuses}, key=cells.cell_sort_key
    ):
        total = sum(tally.statuses[(cell, s)] for s in statuses)
        lines.append(
            f"  {cell:<12}"
            + "".join(
                f"{_share(tally.statuses[(cell, s)], total):>18}"
                for s in statuses
            )
            + f"{tally.morfflex_only[cell]:>15,}"
        )
    lines += [
        "",
        "Cell agreement (set of Wiktionary vs standard MorfFlex forms)",
    ]
    kinds = (
        "agree",
        "wikt_subset",
        "wikt_superset",
        "overlap",
        "disjoint",
        CELL_MISSING,
        NO_LEMMA,
    )
    lines.append(f"  {'cell':<12}" + "".join(f"{k:>15}" for k in kinds))
    for cell in sorted(
        {c for c, _ in tally.agreements}, key=cells.cell_sort_key
    ):
        lines.append(
            f"  {cell:<12}"
            + "".join(f"{tally.agreements[(cell, k)]:>15,}" for k in kinds)
        )
    lines += ["", "Form provenance (R11) by status"]
    for source in sorted({s for s, _ in tally.by_source}):
        total = sum(tally.by_source[(source, s)] for s in statuses)
        lines.append(
            f"  {source or '-':<15}"
            + "".join(
                f"{_share(tally.by_source[(source, s)], total):>18}"
                for s in statuses
            )
        )
    lines += ["", "Manual worklist (analysis/czech/czech_manual_worklist.csv)"]
    lines.append(f"  rows: {sum(tally.worklist.values()):,}")
    for cell in sorted(tally.worklist, key=cells.cell_sort_key):
        lines.append(f"  {cell:<12} {tally.worklist[cell]:>7,}")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Stage
# ---------------------------------------------------------------------------


def in_domain_lemmas(paths: CzechPaths) -> list[str]:
    """Sorted lemmas of the lexemes :func:`load_lexemes` keeps."""
    found: set[str] = set()
    for row in jsonio.iter_csv(paths.lexicon_full):
        if (
            row["domain"] == "in"
            and row["pos"] in ("noun", "adj")
            and row["paradigm_source"] in ("en", "cs")
        ):
            found.add(row["lemma"])
    return sorted(found)


def _write_lemma_list(ctx: RunContext, paths: CzechPaths) -> None:
    text = "\n".join(in_domain_lemmas(paths)) + "\n"
    target = paths.morfflex_lemmas
    if target.is_file() and target.read_text(encoding="utf-8") == text:
        return
    ctx.write_text(target, text)


def ensure_morfflex_db(
    ctx: RunContext, paths: CzechPaths, words: Path
) -> int | None:
    """Build the MorfFlex projection unless its manifest is current."""
    _write_lemma_list(ctx, paths)

    def body() -> int:
        lemmas = frozenset(
            paths.morfflex_lemmas.read_text(encoding="utf-8").split("\n")
        ) - {""}
        print(
            f"  streaming {paths.morfflex.name} ({len(lemmas):,} lemmas)",
            flush=True,
        )
        stats = morfflex.build_db(
            ctx,
            paths.morfflex,
            paths.morfflex_db,
            lemmas=lemmas,
            vocabulary=leipzig.vocabulary(words),
            meta=(("source_sha256", sha256_file(paths.morfflex)),),
        )
        with contextlib.closing(sqlite3.connect(paths.morfflex_db)) as con:
            con.execute(
                "INSERT INTO meta VALUES ('matched', ?)", (str(stats.matched),)
            )
            con.commit()
        print(
            f"  MorfFlex rows matched {stats.matched:,}, stored {stats.stored:,}",
            flush=True,
        )
        return 0

    return run_stage(
        ctx,
        DB_STAGE,
        inputs=[paths.morfflex, paths.morfflex_lemmas, words],
        params={"version": DB_VERSION},
        outputs=[paths.morfflex_db],
        body=body,
    )


def run(ctx: RunContext) -> int | None:
    """Stage entry point (see the module docstring)."""
    paths = paths_for(ctx)
    words = leipzig.ensure_words_file(ctx, paths)
    status = ensure_morfflex_db(ctx, paths, words)
    if isinstance(status, int) and status != 0:
        return status
    return run_stage(
        ctx,
        STAGE,
        inputs=[paths.lexicon_full, paths.cells_raw, paths.morfflex_db],
        params={"version": VERSION},
        outputs=[paths.qc_csv, paths.qc_report, paths.worklist],
        body=lambda: _compare(ctx, paths),
    )


def _db_meta(con: sqlite3.Connection) -> dict[str, str]:
    return {
        str(k): str(v) for k, v in con.execute("SELECT key, value FROM meta")
    }


def _compare(ctx: RunContext, paths: CzechPaths) -> int:
    lexemes = load_lexemes(paths)
    tally = Tally()
    worklist: list[dict[str, Any]] = []
    uri = paths.morfflex_db.resolve().as_uri() + "?mode=ro"
    with contextlib.closing(sqlite3.connect(uri, uri=True)) as con:
        meta = _db_meta(con)

        def all_rows() -> Iterator[dict[str, Any]]:
            for lexeme_id in sorted(lexemes):
                lexeme = lexemes[lexeme_id]
                found = morfflex_cells(con, lexeme)
                rows = list(qc_rows(lexeme, found))
                tally.add(lexeme, found, rows)
                statuses = {
                    (str(r["cell"]), str(r["form"])): str(r["status"])
                    for r in rows
                    if r["origin"] == "wiktionary"
                }
                for item in worklist_rows(lexeme, statuses):
                    tally.worklist[str(item["cell"])] += 1
                    worklist.append(item)
                yield from rows

        csvio.write_rows(paths.qc_csv, all_rows(), QC_FIELDS, ctx=ctx)
    csvio.write_tracked(paths.worklist, worklist, WORKLIST_FIELDS, ctx=ctx)
    ctx.write_tracked_text(
        paths.qc_report,
        render_report(tally, meta),
        sources=REPORT_SOURCES,
        nc_aggregates=True,
    )
    print(
        f"  compared {sum(tally.lexemes.values()):,} lexemes; worklist "
        f"{len(worklist):,} rows",
        flush=True,
    )
    return 0
