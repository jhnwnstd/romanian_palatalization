"""Stage ``attest``: corpus attestation of every Italian form (R11, R24).

Dictionary forms are not attestation (R11): Module:it-common generates
plurals nobody uses. A form is attested by its token count in
``ita_news_2024_1M`` (Leipzig), read with
:class:`lexpipe.attest.AttestationTable` with the Italian elided
clitics split off first (l'amico, dell'amico: both ' and ’).

Inputs: ``italian_lexicon_full.csv``, ``italian_plural_forms.csv``, the
form index of the ``extract`` stage and the corpus archive.

Outputs
-------
``data/italian/leipzig/freq/ita_news_2024_1M_freq.csv``: the folded
``word,freq`` table (:func:`lexpipe.frequency.build_freq_csv`), the
Italian counterpart of the Romanian table, under ``data/italian/``
(R2).

``data/italian/italian_attestation.csv`` (gitignored; Leipzig counts are
not redistributed, R9): one row per (lexeme, cell, form), sorted by
lexeme_id, cell, form:

  cell            lemma (the citation form), f.sg (an adjective's
                  feminine singular), pl | m.pl | f.pl (plurals)
  provenance      lemma | explicit | module_default
  count_cs        tokens of the form with its case (Medici ≠ medici);
                  no sentence-initial correction (it needs the
                  sentences file)
  count_folded    tokens of the lowercased form
  per_million     count_folded per million word tokens
  attested_<τ>    folded count reaches τ ∈ :data:`lexpipe.attest.TAUS`
                  (dictionary, ge1, ge0.5pm, ge2pm) and count_cs ≥ 1:
                  a form seen only capitalized is a name (Bocchi,
                  Paschi; R25), not attested (:func:`.leipzig.attested_flags`)
  attested_cs_<τ> count_cs reaches τ (τ ≠ dictionary)
  competitor      for an -i plural of a dorsal-target lemma (amico,
                  parco; not the NDE and frontstem classes): the
                  competing spelling, velar ↔ palatal (amichi for
                  amici, parci for parchi); "" otherwise
  competitor_count_cs / competitor_count_folded   its tokens (empty
                  without a competitor): a dictionary plural, often a
                  Module:it-common default, whose competitor is used
  competitor_other_lemma  the competitor is also a form or lemma of
                  another word (dolci of dolce for dolco's dolchi, luci
                  of luce): its tokens are not evidence for this plural
  form_ambiguous  the form is also a form or lemma of another root:
                  a verb (carichi of caricare), another part of speech,
                  a proper name, another etymology (antichi of antico
                  'ancient' and of antico 'anterior') — from the form
                  index, case-folded
  form_ambiguous_cs   the same, comparing with the case kept (medici is
                  not ambiguous with the family name Medici)
  ambiguous_with  ``lemma:pos`` of the other lexemes, sorted (at most
                  :data:`MAX_AMBIGUOUS`)

Forms are looked up by their case key (respelling accents stripped:
fòcora → focora), the shape they take in running text.

``data/italian/italian_attestation_summary.tsv``: corpus size and
attestation counts per cell kind and τ.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

from lexpipe import attest, csvio, frequency
from lexpipe.attest import AttestationTable
from lexpipe.context import RunContext, run_stage
from lexpipe.langs.italian import (
    extract,
    itcommon,
    keys,
    leipzig,
    morphology,
    store,
)
from lexpipe.langs.italian.layout import Layout

STAGE: Final[str] = "attest"
# 2: reads italian_lexicon_full.csv (the tracked raw lexicon is lean).
# 3: names filter on the τ flags (R25); competitor spellings (R11).
VERSION: Final[int] = 3
CELL_LEMMA: Final[str] = "lemma"
CELL_FSG: Final[str] = "f.sg"
CELL_ORDER: Final[Mapping[str, int]] = {
    CELL_LEMMA: 0,
    CELL_FSG: 1,
    "pl": 2,
    "m.pl": 3,
    "f.pl": 4,
}
LEXEME_POS: Final[frozenset[str]] = frozenset({"noun", "adj"})
MAX_AMBIGUOUS: Final[int] = 8
CS_TAUS: Final[tuple[str, ...]] = tuple(
    t for t in attest.TAUS if t != "dictionary"
)
ATTESTATION_FIELDS: Final[tuple[str, ...]] = (
    "lexeme_id",
    "root_id",
    "word",
    "pos",
    "cell",
    "form",
    "provenance",
    "count_cs",
    "count_folded",
    "per_million",
    *(f"attested_{tau}" for tau in attest.TAUS),
    *(f"attested_cs_{tau}" for tau in CS_TAUS),
    "form_ambiguous",
    "form_ambiguous_cs",
    "ambiguous_with",
    "competitor",
    "competitor_count_cs",
    "competitor_count_folded",
    "competitor_other_lemma",
    "source",
)
SOURCES: Final[tuple[str, str]] = ("wiktionary_en", "leipzig")


@dataclass(frozen=True, slots=True)
class FormRef:
    """One form of one lexeme to attest."""

    lexeme_id: str
    root_id: str
    word: str
    pos: str
    cell: str
    form: str
    provenance: str


PLURAL_CELLS: Final[frozenset[str]] = frozenset({"pl", "m.pl", "f.pl"})


def competitor(ref: FormRef) -> str:
    """The competing spelling of an -i plural of a dorsal-target lemma
    (amici → amichi, parchi → parci, laghi → lagi, boschi → bosci), or
    "" (other cells and lemmas, forms that do not align)."""
    if ref.cell not in PLURAL_CELLS:
        return ""
    lemma = keys.match_key(ref.word)
    stem = morphology.analyze_lemma(lemma)
    if not stem.dorsal or stem.nde_class:
        return ""
    form = ref.form
    analysis = morphology.analyze_form(lemma, stem, keys.match_key(form))
    if analysis.opportunity != morphology.OPPORTUNITY_I:
        return ""
    if analysis.mutation is False:
        return form[:-2] + form[-1]  # chi → ci
    if analysis.mutation is True:
        return form[:-1] + "h" + form[-1]  # ci → chi
    return ""


@dataclass(frozen=True, slots=True)
class IndexRow:
    """A form-index row reduced to what ambiguity needs."""

    key_cs: str
    lemma: str
    pos: str
    ety: str


def feminine_singulars(row: Mapping[str, str]) -> list[tuple[str, str]]:
    """``(form, provenance)`` of an adjective's feminine singular."""
    if row["pos"] != "adj":
        return []
    stated = csvio.split_list(row["feminine_forms"])
    if stated:
        return [(form, "explicit") for form in stated]
    word = row["word"]
    if " " in word or "-" in word:
        return []
    feminine = itcommon.make_feminine(word)
    return [(feminine, "module_default")] if feminine != word else []


def form_refs(
    lexicon_rows: Sequence[Mapping[str, str]],
    plural_rows: Iterable[Mapping[str, str]],
) -> list[FormRef]:
    """Every (lexeme, cell, form) to attest, sorted."""
    by_id = {row["lexeme_id"]: row for row in lexicon_rows}
    refs: dict[tuple[str, str, str], FormRef] = {}

    def add(row: Mapping[str, str], cell: str, form: str, prov: str) -> None:
        key = (row["lexeme_id"], cell, form)
        if form and key not in refs:
            refs[key] = FormRef(
                row["lexeme_id"],
                row["root_id"],
                row["word"],
                row["pos"],
                cell,
                form,
                prov,
            )

    for row in lexicon_rows:
        add(row, CELL_LEMMA, row["word"], CELL_LEMMA)
        for form, provenance in feminine_singulars(row):
            add(row, CELL_FSG, form, provenance)
    for plural in plural_rows:
        lexeme = by_id.get(plural["lexeme_id"])
        if lexeme is not None and plural["form"]:
            add(lexeme, plural["cell"], plural["form"], plural["provenance"])
    return sorted(
        refs.values(),
        key=lambda r: (r.lexeme_id, CELL_ORDER[r.cell], r.form),
    )


def load_index(path: Path, wanted: set[str]) -> dict[str, list[IndexRow]]:
    """Form-index rows whose match key is in *wanted*."""
    meta, fields, rows = store.read_tsv(path)
    if meta.get("format") != extract.FORM_INDEX_FORMAT:
        raise store.StaleDerivativeError(f"{path}: unexpected format {meta}")
    column = {name: index for index, name in enumerate(fields)}
    out: dict[str, list[IndexRow]] = defaultdict(list)
    for row in rows:
        key = row[column["key"]]
        if key in wanted:
            out[key].append(
                IndexRow(
                    row[column["key_cs"]],
                    row[column["lemma"]],
                    row[column["pos"]],
                    row[column["ety"]],
                )
            )
    return dict(out)


@dataclass(frozen=True, slots=True)
class RootIndex:
    """root_ids of the lexicon by (word, pos) and (word, pos, ety)."""

    by_word: Mapping[tuple[str, str], frozenset[str]]
    by_entry: Mapping[tuple[str, str, str], frozenset[str]]

    @classmethod
    def build(cls, lexicon_rows: Iterable[Mapping[str, str]]) -> RootIndex:
        by_word: dict[tuple[str, str], set[str]] = defaultdict(set)
        by_entry: dict[tuple[str, str, str], set[str]] = defaultdict(set)
        for row in lexicon_rows:
            by_word[(row["word"], row["pos"])].add(row["root_id"])
            key = (row["word"], row["pos"], row["etymology_number"])
            by_entry[key].add(row["root_id"])
        return cls(
            {k: frozenset(v) for k, v in by_word.items()},
            {k: frozenset(v) for k, v in by_entry.items()},
        )

    def same_root(self, row: IndexRow, root_id: str) -> bool:
        if row.pos not in LEXEME_POS:
            return False
        if row.ety:
            roots = self.by_entry.get(
                (row.lemma, row.pos, row.ety), frozenset()
            )
        else:
            roots = self.by_word.get((row.lemma, row.pos), frozenset())
        return root_id in roots


def ambiguity(
    ref: FormRef,
    index: Mapping[str, Sequence[IndexRow]],
    roots: RootIndex,
) -> tuple[bool, bool, list[str]]:
    """``(form_ambiguous, form_ambiguous_cs, ambiguous_with)``."""
    key, key_cs = keys.form_keys(ref.form)
    others = [
        row
        for row in index.get(key, ())
        if not roots.same_root(row, ref.root_id)
    ]
    with_case = [row for row in others if row.key_cs == key_cs]
    names = sorted({f"{row.lemma}:{row.pos}" for row in others})
    return bool(others), bool(with_case), names[:MAX_AMBIGUOUS]


def attestation_row(
    ref: FormRef,
    table: AttestationTable,
    index: Mapping[str, Sequence[IndexRow]],
    roots: RootIndex,
) -> dict[str, Any]:
    """One output row (see the module docstring)."""
    text = keys.case_key(ref.form)
    count_cs = table.count(text, case_sensitive=True)
    count_folded = table.count(text)
    ambiguous, ambiguous_cs, names = ambiguity(ref, index, roots)
    other = competitor(ref)
    other_text = keys.case_key(other) if other else ""
    other_lemma = (
        any(
            not roots.same_root(row, ref.root_id)
            for row in index.get(keys.match_key(other), ())
        )
        if other
        else None
    )
    row: dict[str, Any] = {
        "lexeme_id": ref.lexeme_id,
        "root_id": ref.root_id,
        "word": ref.word,
        "pos": ref.pos,
        "cell": ref.cell,
        "form": ref.form,
        "provenance": ref.provenance,
        "count_cs": count_cs,
        "count_folded": count_folded,
        "per_million": round(table.per_million(text), 6),
        "form_ambiguous": ambiguous,
        "form_ambiguous_cs": ambiguous_cs,
        "ambiguous_with": names,
        "competitor": other,
        "competitor_count_cs": (
            table.count(other_text, case_sensitive=True) if other else None
        ),
        "competitor_count_folded": (
            table.count(other_text) if other else None
        ),
        "competitor_other_lemma": other_lemma,
        "source": list(SOURCES),
    }
    for tau, value in leipzig.attested_flags(table, text).items():
        row[f"attested_{tau}"] = value
    for tau, value in table.attested_by_tau(
        text, CS_TAUS, case_sensitive=True
    ).items():
        row[f"attested_cs_{tau}"] = value
    return row


def summary(rows: Sequence[Mapping[str, Any]], table: AttestationTable) -> str:
    """``italian_attestation_summary.tsv`` (aggregate counts only)."""
    counts: Counter[str] = Counter()
    counts["corpus:word_tokens"] = table.total_tokens
    counts["corpus:types_folded"] = len(table.counts_folded)
    for row in rows:
        kind = (
            row["cell"] if row["cell"] in (CELL_LEMMA, CELL_FSG) else "plural"
        )
        counts[f"forms:{kind}"] += 1
        for tau in attest.TAUS:
            if row[f"attested_{tau}"]:
                counts[f"attested:{kind}:{tau}"] += 1
        if row["form_ambiguous"]:
            counts[f"ambiguous:{kind}"] += 1
        if row["form_ambiguous_cs"]:
            counts[f"ambiguous_cs:{kind}"] += 1
        if row["competitor"]:
            counts[f"competitor:{row['provenance']}"] += 1
            if (
                row["competitor_count_cs"]
                and not row["competitor_other_lemma"]
            ):
                counts[f"competitor_attested:{row['provenance']}"] += 1
        if row["attested_cs_ge1"] is False and row["count_folded"]:
            counts[f"capitalized_only:{kind}"] += 1
    lines = ["key\tvalue"]
    lines.extend(f"{key}\t{counts[key]}" for key in sorted(counts))
    return "\n".join(lines) + "\n"


def _body(ctx: RunContext, layout: Layout, words: Path) -> int:
    corpus_id = ctx.language.leipzig_id
    ctx.check_writable(layout.leipzig_freq)
    frequency.build_freq_csv(
        corpus_id,
        layout.leipzig,
        archive=layout.leipzig_archive,
        expected_size=leipzig.expected_size(ctx, layout),
        log=lambda message: print(message),
    )
    table = leipzig.attestation_table(ctx, words)
    lexicon_rows = csvio.read_rows(layout.lexicon_full)
    refs = form_refs(lexicon_rows, csvio.read_rows(layout.plural_forms))
    wanted = {keys.match_key(ref.form) for ref in refs}
    wanted |= {
        keys.match_key(other) for ref in refs if (other := competitor(ref))
    }
    index = load_index(layout.form_index, wanted)
    roots = RootIndex.build(lexicon_rows)
    rows = [attestation_row(ref, table, index, roots) for ref in refs]
    csvio.write_rows(layout.attestation, rows, ATTESTATION_FIELDS, ctx=ctx)
    ctx.write_text(layout.attestation_summary, summary(rows, table))
    attested = sum(1 for row in rows if row["attested_ge1"])
    print(
        f"  {len(rows):,} forms, {attested:,} attested (≥1 token) in "
        f"{table.total_tokens:,} word tokens"
    )
    return 0


def outputs(layout: Layout) -> list[Path]:
    return [
        layout.leipzig_freq,
        layout.attestation,
        layout.attestation_summary,
    ]


def run(ctx: RunContext) -> int | None:
    """Run the ``attest`` stage (skipped when its manifest matches)."""
    layout = Layout(ctx)
    words = leipzig.words_file(ctx, layout)
    return run_stage(
        ctx,
        STAGE,
        inputs=[
            layout.lexicon_full,
            layout.plural_forms,
            layout.form_index,
            words,
        ],
        params={
            "version": VERSION,
            "elision_clitics": sorted(ctx.language.elision_clitics),
            "taus": list(attest.TAUS),
        },
        outputs=outputs(layout),
        body=lambda: _body(ctx, layout, words),
    )
