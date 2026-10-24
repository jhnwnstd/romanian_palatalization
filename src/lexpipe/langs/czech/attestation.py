"""Stage ``attest``: corpus attestation of every cell form (R11, R24).

Dictionary forms are not attestation: cs-ndecl generates cells nobody
uses (biolozích, snězích: 0 tokens). A form is attested by its token
count in Leipzig ``ces_news_2024_1M`` (:mod:`.leipzig`; Czech has no
elided clitics to split).

Counts
------
Per cell form of an in-domain lexeme (both editions' rows):

    count_cs        tokens of the form exactly as written (Češi, kluci)
    count_folded    case-folded tokens (Kluci + kluci + KLUCI)
    count_cap       tokens of the form with a capital first letter, for
                    a lowercase form: the sentence-initial uses
    cap_claimed     that capitalized string is a form of another
                    lexeme (Poláci 'Poles' for polák 'pochard'), so its
                    tokens cannot be credited
    count           count_cs + count_cap unless cap_claimed: the
                    case-sensitive count with the sentence-initial
                    correction R24 asks for, used for every threshold

Ethnonyms are capitalized (Češi): their count is ``count_cs``. The
thresholds τ (:data:`lexpipe.attest.TAUS`) are reported as
``attested_ge1``, ``attested_ge0.5pm`` and ``attested_ge2pm``: the
core's :func:`lexpipe.attest.attested_count` applied to ``count`` (rates
per million word tokens, exact rational comparison, the arithmetic of
:meth:`lexpipe.attest.AttestationTable.attested`); at τ = dictionary
every row counts.

Ambiguity
---------
``form_ambiguous`` (exact spelling) and ``form_ambiguous_folded`` are
True when the form index (:mod:`.extract`) lists an owner other than
the lexeme itself or an in-domain lexeme of the same root
(``root_id``): a verb form, a proper name, a homograph with another
paradigm (šachu: šach 'check' and šach 'shah'). Owners in English
Wiktionary are compared by base lexeme_id, Czech Wiktionary owners and
form-of targets by (word, part of speech). Two kinds of index rows
are not owners (:func:`counts_as_owner`): an ``alt_of`` row (a spelling
variant of the same word: čůrák, čurák) and a ``form_of`` row whose
target is a multiword phrase (cs.wiktionary files kapr under "vánoční
kapr"): a single-word token is not a form of a phrase.

Lexeme frequencies and the CNK cross-check
------------------------------------------
``czech_lexeme_frequency.csv`` gives per in-domain lexeme the lemma
count (same correction) and the paradigm sum over distinct attested
primary forms (R24: both bases), plus the lemma frequency in the CNK
SYN2015 list (CC BY 4.0; ≥ 10 by construction). The report gives
Spearman and log-Pearson correlations between Leipzig and CNK over the
in-domain nouns found in both.

Homographs inflate both bases: jak 'yak' gets the 33,277 tokens of the
conjunction jak, válek (a diminutive of vál) the tokens of válku and
války, forms of válka 'war'. ``lemma_ambiguous`` flags a lemma string
that the form index assigns to another lexeme (the ``form_ambiguous``
test applied to the lemma), and ``paradigm_sum_unambiguous`` sums the
distinct primary forms that are not ``form_ambiguous`` (the lemma
included when it is not ambiguous): the frequency basis the learners
and the V sweep report as a sensitivity (R24).

Leipzig-derived numbers are untracked until its terms are checked
(R9): every output stays in ``data/czech/``.
"""

from __future__ import annotations

import contextlib
import math
import sqlite3
from collections import defaultdict
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Final

from lexpipe import csvio
from lexpipe.attest import TAUS, AttestationTable, attested_count
from lexpipe.context import RunContext, run_stage
from lexpipe.frequency import paradigm_frequency
from lexpipe.langs.czech import jsonio, leipzig
from lexpipe.langs.czech.layout import CzechPaths, paths_for

STAGE: Final[str] = "attest"
# 2: reads czech_lexicon_full.csv; core attested_count.
# 3: lemma_ambiguous and paradigm_sum_unambiguous (homograph lemmas).
# 4: alt_of rows and form_of rows of multiword phrases are no owners.
VERSION: Final[str] = "4"
MILLION: Final[int] = 1_000_000
CORPUS_TAUS: Final[tuple[str, ...]] = tuple(
    t for t in TAUS if t != "dictionary"
)
AMBIGUITY_SHOWN: Final[int] = 5

ATTESTATION_FIELDS: Final[tuple[str, ...]] = (
    "lexeme_id",
    "cell",
    "form",
    "source",
    "primary",
    "count_cs",
    "count_cap",
    "cap_claimed",
    "count_folded",
    "count",
    "per_million",
    *(f"attested_{tau}" for tau in CORPUS_TAUS),
    "form_ambiguous",
    "form_ambiguous_folded",
    "ambiguous_with",
)
FREQUENCY_FIELDS: Final[tuple[str, ...]] = (
    "lexeme_id",
    "lemma",
    "pos",
    "lemma_count_cs",
    "lemma_count",
    "lemma_count_folded",
    "lemma_ambiguous",
    "paradigm_sum",
    "paradigm_sum_folded",
    "paradigm_sum_unambiguous",
    "n_forms",
    "n_forms_attested",
    "cnk_lemma_freq",
)


# ---------------------------------------------------------------------------
# Lexemes and owners
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class LexemeInfo:
    """What attestation needs to know about one lexeme."""

    lexeme_id: str
    lemma: str
    pos: str
    root_id: str
    in_domain: bool
    base_id: str


def base_id(lexeme_id: str, lemma: str, pos: str, etymology: str) -> str:
    """The lexeme_id before tie-breaking suffixes (R4)."""
    prefix = f"cs:{lemma}:{pos}:{etymology}:"
    if lexeme_id.startswith(prefix):
        return lexeme_id[: len(prefix) + 8]
    return lexeme_id


def load_lexicon(paths: CzechPaths) -> dict[str, LexemeInfo]:
    """Every lexeme of the raw lexicon (names included)."""
    out: dict[str, LexemeInfo] = {}
    for row in jsonio.iter_csv(paths.lexicon_full):
        out[row["lexeme_id"]] = LexemeInfo(
            lexeme_id=row["lexeme_id"],
            lemma=row["lemma"],
            pos=row["pos"],
            root_id=row["root_id"],
            in_domain=row["domain"] == "in",
            base_id=base_id(
                row["lexeme_id"], row["lemma"], row["pos"], row["etymology"]
            ),
        )
    return out


def compatible_owners(
    lexicons: Mapping[str, LexemeInfo],
) -> dict[str, frozenset[str]]:
    """``lexeme_id -> owner keys that do not make a form ambiguous``:
    the lexeme itself and the in-domain lexemes of its root."""
    by_root: dict[str, list[LexemeInfo]] = defaultdict(list)
    for info in lexicons.values():
        if info.in_domain:
            by_root[info.root_id].append(info)
    out: dict[str, frozenset[str]] = {}
    for info in lexicons.values():
        keys = {f"id:{info.base_id}", f"w:{info.pos}:{info.lemma}"}
        if info.in_domain:
            for other in by_root[info.root_id]:
                keys.add(f"id:{other.base_id}")
                keys.add(f"w:{other.pos}:{other.lemma}")
        out[info.lexeme_id] = frozenset(keys)
    return out


def counts_as_owner(kind: str, target: str) -> bool:
    """Does a form-index row make a form ambiguous (module docstring)?
    Not a spelling variant (``alt_of``) nor a form-of row of a
    multiword phrase."""
    if kind == "alt_of":
        return False
    return not (kind == "form_of" and " " in target.strip())


def owner_key(
    kind: str,
    edition: str,
    pos: str,
    owner_word: str,
    owner_id: str,
    target: str,
) -> str:
    """Identity of one form-index owner (see the module docstring)."""
    if kind in ("lemma", "paradigm"):
        if edition == "en" and owner_id:
            return f"id:{owner_id}"
        return f"w:{pos}:{owner_word}"
    return f"w:{pos}:{target}"


@dataclass
class Owners:
    """Owners of the forms of interest, exact and case-folded."""

    exact: dict[str, set[str]] = field(
        default_factory=lambda: defaultdict(set)
    )
    folded: dict[str, set[str]] = field(
        default_factory=lambda: defaultdict(set)
    )


def load_owners(paths: CzechPaths, folded_forms: set[str]) -> Owners:
    """Index rows of the given case-folded forms."""
    owners = Owners()
    uri = paths.form_index.resolve().as_uri() + "?mode=ro"
    with contextlib.closing(sqlite3.connect(uri, uri=True)) as con:
        con.execute("CREATE TEMP TABLE wanted (form_folded TEXT PRIMARY KEY)")
        con.executemany(
            "INSERT INTO wanted VALUES (?)",
            ((f,) for f in sorted(folded_forms)),
        )
        rows = con.execute(
            "SELECT f.form, f.form_folded, f.kind, f.edition, f.pos, "
            "f.owner_word, f.owner_id, f.target FROM forms f "
            "JOIN wanted w ON f.form_folded = w.form_folded"
        )
        for form, folded, kind, edition, pos, word, owner_id, target in rows:
            if not counts_as_owner(kind, target):
                continue
            key = owner_key(kind, edition, pos, word, owner_id, target)
            owners.exact[form].add(key)
            owners.folded[folded].add(key)
    return owners


# ---------------------------------------------------------------------------
# Counting
# ---------------------------------------------------------------------------


def capitalized(form: str) -> str:
    return form[:1].upper() + form[1:]


def attestation_row(
    row: Mapping[str, str],
    table: AttestationTable,
    owners: Owners,
    compatible: frozenset[str],
) -> dict[str, Any]:
    """Counts, thresholds and ambiguity of one cell form."""
    form = row["form"]
    count_cs = table.count(form, case_sensitive=True)
    count_cap = 0
    cap_claimed = False
    if form[:1].islower():
        cap = capitalized(form)
        count_cap = table.count(cap, case_sensitive=True)
        cap_claimed = bool(owners.exact.get(cap, set()) - compatible)
    count = count_cs + (0 if cap_claimed else count_cap)
    others_exact = owners.exact.get(form, set()) - compatible
    others_folded = owners.folded.get(form.lower(), set()) - compatible
    return {
        "lexeme_id": row["lexeme_id"],
        "cell": row["cell"],
        "form": form,
        "source": row["source"],
        "primary": csvio.parse_bool(row["primary"]),
        "count_cs": count_cs,
        "count_cap": count_cap,
        "cap_claimed": cap_claimed,
        "count_folded": table.count(form),
        "count": count,
        "per_million": round(count * MILLION / table.total_tokens, 4),
        **{
            f"attested_{tau}": attested_count(count, tau, table.total_tokens)
            for tau in CORPUS_TAUS
        },
        "form_ambiguous": bool(others_exact),
        "form_ambiguous_folded": bool(others_folded),
        "ambiguous_with": sorted(others_folded)[:AMBIGUITY_SHOWN],
    }


# ---------------------------------------------------------------------------
# Correlation
# ---------------------------------------------------------------------------


def ranks(values: Sequence[float]) -> list[float]:
    """Average ranks (1-based), ties sharing their mean rank."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    out = [0.0] * len(values)
    index = 0
    while index < len(order):
        end = index
        while (
            end + 1 < len(order)
            and values[order[end + 1]] == values[order[index]]
        ):
            end += 1
        mean = (index + end) / 2 + 1
        for k in range(index, end + 1):
            out[order[k]] = mean
        index = end + 1
    return out


def pearson(xs: Sequence[float], ys: Sequence[float]) -> float | None:
    """Pearson's r; None for fewer than 3 points or a constant series."""
    n = len(xs)
    if n < 3:
        return None
    mean_x, mean_y = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mean_x) ** 2 for x in xs)
    syy = sum((y - mean_y) ** 2 for y in ys)
    if sxx == 0 or syy == 0:
        return None
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys, strict=True))
    return sxy / math.sqrt(sxx * syy)


def spearman(xs: Sequence[float], ys: Sequence[float]) -> float | None:
    """Spearman's rho (Pearson on average ranks)."""
    return pearson(ranks(xs), ranks(ys))


def _fmt(value: float | None) -> str:
    return "n/a" if value is None else f"{value:.3f}"


# ---------------------------------------------------------------------------
# Stage
# ---------------------------------------------------------------------------


def run(ctx: RunContext) -> int | None:
    """Stage entry point (see the module docstring)."""
    paths = paths_for(ctx)
    words = leipzig.ensure_words_file(ctx, paths)
    freq_csv = paths.leipzig.freq_csv(ctx.language.leipzig_id)
    return run_stage(
        ctx,
        STAGE,
        inputs=[
            paths.lexicon_full,
            paths.cells_raw,
            paths.form_index,
            words,
            paths.cnk_lemma,
        ],
        params={"version": VERSION, "taus": list(TAUS)},
        outputs=[
            freq_csv,
            paths.attestation,
            paths.lexeme_frequency,
            paths.attestation_report,
        ],
        body=lambda: _attest(ctx, paths),
    )


def _cell_rows(
    paths: CzechPaths, lexicon: Mapping[str, LexemeInfo]
) -> Iterator[dict[str, str]]:
    for row in jsonio.iter_csv(paths.cells_raw):
        info = lexicon.get(row["lexeme_id"])
        if info is not None and info.in_domain:
            yield row


def _attest(ctx: RunContext, paths: CzechPaths) -> int:
    leipzig.ensure_freq_csv(ctx, paths)
    words = leipzig.ensure_words_file(ctx, paths)
    table = leipzig.attestation_table(words)
    lexicon = load_lexicon(paths)
    compatible = compatible_owners(lexicon)
    folded: set[str] = set()
    for row in _cell_rows(paths, lexicon):
        folded.add(row["form"].lower())
        folded.add(capitalized(row["form"]).lower())
    for info in lexicon.values():
        if info.in_domain:
            folded.add(info.lemma.lower())
    owners = load_owners(paths, folded)
    primary_counts: dict[str, dict[str, int]] = defaultdict(dict)
    ambiguous: dict[str, set[str]] = defaultdict(set)
    stats: dict[str, int] = defaultdict(int)

    def rows() -> Iterator[dict[str, Any]]:
        for row in _cell_rows(paths, lexicon):
            out = attestation_row(
                row, table, owners, compatible[row["lexeme_id"]]
            )
            if out["primary"]:
                primary_counts[row["lexeme_id"]][row["form"]] = int(
                    out["count"]
                )
                if out["form_ambiguous"]:
                    ambiguous[row["lexeme_id"]].add(row["form"])
                stats["primary_rows"] += 1
                stats["primary_ambiguous"] += int(bool(out["form_ambiguous"]))
                stats["primary_cap_added"] += int(
                    bool(out["count_cap"]) and not out["cap_claimed"]
                )
                for tau in CORPUS_TAUS:
                    stats[f"primary_{tau}"] += int(
                        bool(out[f"attested_{tau}"])
                    )
            yield out

    csvio.write_rows(paths.attestation, rows(), ATTESTATION_FIELDS, ctx=ctx)
    cnk = leipzig.load_cnk_lemmas(paths.cnk_lemma)
    frequency = [
        _frequency_row(
            info,
            table,
            owners,
            compatible,
            primary_counts,
            cnk,
            ambiguous.get(info.lexeme_id, set()),
        )
        for info in sorted(lexicon.values(), key=lambda i: i.lexeme_id)
        if info.in_domain and info.pos in ("noun", "adj")
    ]
    csvio.write_rows(
        paths.lexeme_frequency, frequency, FREQUENCY_FIELDS, ctx=ctx
    )
    report = _report(table, stats, frequency)
    ctx.write_text(paths.attestation_report, report)
    print(report, flush=True)
    return 0


def _frequency_row(
    info: LexemeInfo,
    table: AttestationTable,
    owners: Owners,
    compatible: Mapping[str, frozenset[str]],
    primary_counts: Mapping[str, Mapping[str, int]],
    cnk: Mapping[str, int],
    ambiguous: frozenset[str] | set[str] = frozenset(),
) -> dict[str, Any]:
    """One row of ``czech_lexeme_frequency.csv`` (module docstring);
    *ambiguous* are the lexeme's primary forms flagged
    ``form_ambiguous``."""
    lemma_row = attestation_row(
        {
            "lexeme_id": info.lexeme_id,
            "cell": "",
            "form": info.lemma,
            "source": "",
            "primary": "True",
        },
        table,
        owners,
        compatible[info.lexeme_id],
    )
    counts = dict(primary_counts.get(info.lexeme_id, {}))
    counts.setdefault(info.lemma, int(lemma_row["count"]))
    lemma_ambiguous = bool(lemma_row["form_ambiguous"])
    unambiguous = {
        form: count
        for form, count in counts.items()
        if form not in ambiguous
        and not (form == info.lemma and lemma_ambiguous)
    }
    return {
        "lexeme_id": info.lexeme_id,
        "lemma": info.lemma,
        "pos": info.pos,
        "lemma_count_cs": lemma_row["count_cs"],
        "lemma_count": lemma_row["count"],
        "lemma_count_folded": lemma_row["count_folded"],
        "lemma_ambiguous": lemma_ambiguous,
        "paradigm_sum": sum(counts.values()),
        "paradigm_sum_folded": paradigm_frequency(counts, table.counts_folded),
        "paradigm_sum_unambiguous": sum(unambiguous.values()),
        "n_forms": len(counts),
        "n_forms_attested": sum(1 for c in counts.values() if c > 0),
        "cnk_lemma_freq": cnk.get(info.lemma),
    }


def _report(
    table: AttestationTable,
    stats: Mapping[str, int],
    frequency: Sequence[Mapping[str, Any]],
) -> str:
    total_rows = stats.get("primary_rows", 0)
    lines = [
        "Czech attestation (Leipzig ces_news_2024_1M); untracked (R9)",
        "",
        f"corpus word tokens: {table.total_tokens:,}",
        f"primary cell-form rows: {total_rows:,}",
    ]
    for tau in CORPUS_TAUS:
        count = stats.get(f"primary_{tau}", 0)
        share = 100 * count / total_rows if total_rows else 0.0
        lines.append(f"  attested {tau:<8} {count:>9,} ({share:5.1f}%)")
    lines.append(f"  form_ambiguous: {stats.get('primary_ambiguous', 0):,}")
    lines.append(
        "  rows credited with sentence-initial capitals: "
        f"{stats.get('primary_cap_added', 0):,}"
    )
    nouns = [r for r in frequency if r["pos"] == "noun"]
    both = [r for r in nouns if r["cnk_lemma_freq"] is not None]
    lines += [
        "",
        "CNK SYN2015 lemma-list cross-check (in-domain nouns)",
        f"  nouns: {len(nouns):,}; in the CNK list (freq >= 10): {len(both):,}",
    ]
    for basis in ("lemma_count", "paradigm_sum"):
        xs = [float(r[basis]) for r in both]
        ys = [float(r["cnk_lemma_freq"]) for r in both]
        log_r = pearson(
            [math.log10(x + 1) for x in xs], [math.log10(y + 1) for y in ys]
        )
        lines.append(
            f"  {basis:<13} Spearman rho {_fmt(spearman(xs, ys))}; "
            f"Pearson r of log10(1+f) {_fmt(log_r)}"
        )
    return "\n".join(lines) + "\n"
