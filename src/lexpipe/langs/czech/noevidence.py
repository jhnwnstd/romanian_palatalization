"""Stage ``noevidence``: roots without alternation evidence × derivatives (R21).

The plural-less prediction (inflection dependence): a root the learner
never hears alternate is stored prespecified, so its derivatives should
not palatalize. "No alternation evidence" is decided **from the corpus**,
not from the dictionary paradigm, because dictionaries generate unused
forms (vzduších, snězích, biolozích: 0 tokens). At τ the base class of
a noun root (its E0 inflection) is

    alternating      an attested P form shows the change (kluci, lize)
    non-alternating  no such form, but an attested competing N
                     allomorph (-u, -ové/-é, -ách) or a P ending without
                     the change: the root is seen avoiding the trigger
    no evidence      neither (only GEN/INS/... forms attested, or none):
                     velar -ko neuters with -u/-ách only fall here or in
                     non-alternating
    (no alternation evidence = non-alternating + no evidence)

Derivatives come from the derivations stage (:mod:`.derivations`;
candidates in both spellings, R21): every derivative whose suffix
palatalizes the base's target (P suffixes: -ní/-ný, -it, -ek, -ka …)
attested at τ, coded at the base × suffix level as ``pal`` (only
palatalized), ``plain`` (only unchanged: barokní) or ``both``.

Reports: 3×2 tables per suffix and pooled, by base stratum and by
derivative stratum (R22, Wiktionary in the tracked report); the pooled
table again without the learned -ik(a)/-ikum → -ický units (kybernetika
→ kybernetický: the adjective may be borrowed whole, Latin -icus,
German -isch, rather than formed as k + -ský); the pooled table over
derivatives that are **internal formations** only (the decisive test:
an inherited derivative was formed in Proto-Slavic and is stored, a
borrowed one was formed elsewhere; a Czech coinage such as bioložka
shows what the grammar does now); the velar neuters on -ko/-go/-ho/-cho
apart; and the counterexamples (no-alternation-evidence bases with a
palatalized derivative) listed by stratum at τ = dictionary, as words
without counts. Thresholds include the strict ``ge1_strict`` and
``ge2pm_strict`` (:data:`.roots.STRICT_TAUS`: a base form's tokens
count only when the form is neither ambiguous with another lexeme nor
shared with a cell outside its copy group). Basis ``wikt`` (Wiktionary
links and headword candidates) is the main one; ``linked`` and
``wikt+corpus`` are sensitivity rows; the NC section adds DeriNet-only
pairs and the czetyl/DeriNet strata (aggregates only).

Outputs: ``analysis/czech/czech_noevidence.txt``,
``data/czech/nc/czech_noevidence_pairs.csv`` and, gitignored because
the lists are filtered by per-word corpus attestation (R9),
``data/czech/czech_noevidence_counterexamples.txt`` (the counterexample
words at τ ≥ 1 token).
"""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Final

from lexpipe import csvio
from lexpipe.attest import TAUS
from lexpipe.context import RunContext, run_stage
from lexpipe.langs.czech import jsonio, ortho, roots
from lexpipe.langs.czech import suffixes as sx
from lexpipe.langs.czech.analysis_layout import AnalysisPaths, analysis_paths
from lexpipe.langs.czech.attestation import CORPUS_TAUS
from lexpipe.langs.czech.etymology import INTERNAL, STRATA, UNKNOWN
from lexpipe.langs.czech.layout import CzechPaths, paths_for

STAGE: Final[str] = "noevidence"
# 3: licence-gated report (its NC section counts DeriNet pairs and
# czetyl/DeriNet strata).
# 4: strict thresholds; internal-formation table; the corpus-filtered
# counterexample words moved to data/czech/ (R9; audit 2026-10-23).
VERSION: Final[str] = "4"
REPORT_SOURCES: Final[tuple[str, ...]] = (
    "wiktionary_en",
    "wiktionary_cs",
    "derinet",
    "czetyl",
)
ALTERNATING: Final[str] = "alternating"
NON_ALTERNATING: Final[str] = "non-alternating"
NO_EVIDENCE: Final[str] = "no evidence"
BASE_CLASSES: Final[tuple[str, ...]] = (
    ALTERNATING,
    NON_ALTERNATING,
    NO_EVIDENCE,
)
POOLED: Final[str] = "(no alternation evidence)"
OUTCOMES: Final[tuple[str, ...]] = ("pal", "plain", "both")
TARGETS: Final[tuple[str, ...]] = (ortho.VELAR, ortho.R_CLASS)
REPORT_TAUS: Final[tuple[str, ...]] = (
    "dictionary",
    "ge1",
    "ge1_strict",
    "ge2pm",
    "ge2pm_strict",
)
INTERNAL_DERIVATIVES: Final[str] = (
    "all P suffixes, derivatives internal formations only"
)
ORDER: Final[tuple[str, ...]] = ("linked", "wikt", "wikt+corpus", "all")
LIST_LIMIT: Final[int] = 40
NEUTER_ENDINGS: Final[tuple[str, ...]] = ("ko", "go", "ho", "cho")
LEARNED_IC_BASES: Final[tuple[str, ...]] = ("ik", "ika", "ikum", "ikus")
ALL_P: Final[str] = "all P suffixes"
WITHOUT_IC: Final[str] = "all P suffixes, without -ik(a) → -ický"

PAIR_FIELDS: Final[tuple[str, ...]] = (
    "root",
    "base_lemma",
    "target_class",
    "genders",
    "base_stratum_wikt",
    "base_stratum",
    *(f"base_class_{tau}" for tau in TAUS),
    "derivative",
    "suffix",
    "outcome",
    "basis",
    "derivative_stratum_wikt",
    "derivative_stratum",
    "count",
    *(f"attested_{tau}" for tau in CORPUS_TAUS),
    "source",
    "licence_class",
)


def base_class(record: roots.RootRecord, tau: str) -> str:
    """Base class of a noun root at τ (module docstring)."""
    seen_n = False
    for item in record.realizations:
        if not item.attested_at(tau):
            continue
        if item.allomorph == ortho.P:
            if item.mutation is True:
                return ALTERNATING
            if item.mutation is False:
                seen_n = True
        elif item.allomorph == ortho.N and item.n_kind == ortho.N_ALLOMORPH:
            seen_n = True
    return NON_ALTERNATING if seen_n else NO_EVIDENCE


def learned_ic(lemma: str, suffix: str, derivative: str) -> bool:
    """-ik/-ika/-ikum + -ický (kritika → kritický, cynik → cynický): the
    learned -ický adjective, which may be borrowed whole or carry the
    suffix -ický rather than show k + -ský → -cký."""
    return (
        suffix == "-ský/-cký"
        and derivative.lower().endswith("ický")
        and lemma.lower().endswith(LEARNED_IC_BASES)
    )


@dataclass
class Unit:
    """One base × suffix unit at one τ and basis."""

    root: str
    lemma: str
    suffix: str
    base_class: str
    base_stratum: str
    outcomes: set[str] = field(default_factory=set)
    derivatives: list[tuple[str, str, str]] = field(default_factory=list)

    @property
    def outcome(self) -> str:
        if self.outcomes == {"pal", "plain"}:
            return "both"
        return next(iter(self.outcomes))

    @property
    def learned(self) -> bool:
        """Every derivative of the unit is a learned -ický adjective."""
        return all(
            learned_ic(self.lemma, self.suffix, word)
            for word, _, _ in self.derivatives
        )


def _attested(row: Mapping[str, str], tau: str) -> bool:
    """A derivative's τ flag (strict τ: the plain one; a derivative's
    forms are already its own, :mod:`.derivations`)."""
    if tau == "dictionary":
        return True
    return csvio.parse_bool(row[f"attested_{roots.plain_tau(tau)}"]) is True


def units(
    records: Mapping[str, roots.RootRecord],
    lexeme_key: Mapping[str, str],
    pairs: Sequence[Mapping[str, str]],
    tau: str,
    basis: str,
    *,
    stratum_field: str = "base_stratum_wikt",
    derivative_stratum_field: str = "derivative_stratum_wikt",
    derivative_strata: frozenset[str] | None = None,
) -> list[Unit]:
    """Base × suffix units at τ on *basis* (P suffixes that palatalize
    the base's target; outcomes pal / plain), from derivatives of
    *derivative_strata* only when given."""
    limit = ORDER.index(basis)
    found: dict[tuple[str, str], Unit] = {}
    for row in pairs:
        key = lexeme_key.get(row["base_lexeme_id"])
        if key is None or ORDER.index(row["basis"]) > limit:
            continue
        if (
            derivative_strata is not None
            and (row[derivative_stratum_field] or UNKNOWN)
            not in derivative_strata
        ):
            continue
        record = records[key]
        if record.target_class not in TARGETS:
            continue
        if row["suffix_class"] != sx.P_CLASS or row["outcome"] not in (
            sx.PAL,
            sx.EXCEPTION,
        ):
            continue
        if not _attested(row, tau):
            continue
        unit = found.setdefault(
            (key, row["suffix"]),
            Unit(
                root=key,
                lemma=record.lemma,
                suffix=row["suffix"],
                base_class=base_class(record, tau),
                base_stratum=row[stratum_field] or UNKNOWN,
            ),
        )
        outcome = "pal" if row["outcome"] == sx.PAL else "plain"
        unit.outcomes.add(outcome)
        unit.derivatives.append(
            (
                row["derivative"],
                outcome,
                row[derivative_stratum_field] or UNKNOWN,
            )
        )
    return [found[k] for k in sorted(found)]


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------


def crosstab(
    items: Sequence[Unit], group: Callable[[Unit], bool] | None = None
) -> dict[str, Counter[str]]:
    """``base class -> Counter(outcome)`` (+ the pooled row)."""
    table: dict[str, Counter[str]] = {
        name: Counter() for name in (*BASE_CLASSES, POOLED)
    }
    for unit in items:
        if group is not None and not group(unit):
            continue
        table[unit.base_class][unit.outcome] += 1
        if unit.base_class != ALTERNATING:
            table[POOLED][unit.outcome] += 1
    return table


def _render_table(table: Mapping[str, Counter[str]], indent: str) -> list[str]:
    lines = [
        f"{indent}{'base class':<28} {'pal':>6} {'plain':>6} {'both':>6} "
        f"{'% pal':>7}"
    ]
    for name in (*BASE_CLASSES, POOLED):
        counts = table[name]
        total = sum(counts.values())
        share = f"{100 * counts['pal'] / total:6.1f}%" if total else "     -"
        lines.append(
            f"{indent}{name:<28} {counts['pal']:>6} {counts['plain']:>6} "
            f"{counts['both']:>6} {share:>7}"
        )
    return lines


def _group(label: str) -> Callable[[Unit], bool] | None:
    if label == ALL_P:
        return None
    if label == WITHOUT_IC:
        return lambda unit: not unit.learned
    return lambda unit: unit.suffix == label


def suffix_tables(items: Sequence[Unit], indent: str = "    ") -> list[str]:
    lines: list[str] = []
    labels = sorted({u.suffix for u in items}, key=_suffix_order)
    for label in [ALL_P, WITHOUT_IC, *labels]:
        table = crosstab(items, _group(label))
        if not any(sum(c.values()) for c in table.values()):
            continue
        lines.append(f"{indent}suffix {label}")
        lines += _render_table(table, indent + "  ")
    return lines


def _suffix_order(label: str) -> int:
    names = [s.label for s in sx.SUFFIXES]
    return names.index(label) if label in names else len(names)


def _in_stratum(stratum: str) -> Callable[[Unit], bool]:
    def test(unit: Unit) -> bool:
        return unit.base_stratum == stratum

    return test


def stratum_tables(items: Sequence[Unit], indent: str = "    ") -> list[str]:
    lines: list[str] = []
    for stratum in STRATA:
        table = crosstab(items, _in_stratum(stratum))
        if not any(sum(c.values()) for c in table.values()):
            continue
        lines.append(f"{indent}base stratum {stratum}")
        lines += _render_table(table, indent + "  ")
    return lines


def counterexamples(items: Sequence[Unit]) -> dict[tuple[str, str], list[str]]:
    """``(base stratum, derivative stratum) -> ["base → derivative"]`` for
    no-alternation-evidence bases with a palatalized derivative; learned
    -ický adjectives are starred."""
    out: dict[tuple[str, str], list[str]] = defaultdict(list)
    for unit in items:
        if unit.base_class == ALTERNATING:
            continue
        for derivative, outcome, stratum in unit.derivatives:
            if outcome == "pal":
                star = (
                    "*"
                    if learned_ic(unit.lemma, unit.suffix, derivative)
                    else ""
                )
                text = f"{unit.lemma} → {derivative}{star}"
                bucket = out[(unit.base_stratum, stratum)]
                if text not in bucket:
                    bucket.append(text)
    return {k: sorted(v) for k, v in sorted(out.items())}


def neuter_velars(
    records: Mapping[str, roots.RootRecord], tau: str
) -> Counter[str]:
    """Base classes of the velar neuters on -ko/-go/-ho/-cho."""
    counts: Counter[str] = Counter()
    for record in records.values():
        if (
            record.pos == roots.NOUN
            and record.target_class == ortho.VELAR
            and record.readings == ["n"]
            and record.lemma.endswith(NEUTER_ENDINGS)
        ):
            counts[base_class(record, tau)] += 1
    return counts


# ---------------------------------------------------------------------------
# Stage
# ---------------------------------------------------------------------------


def run(ctx: RunContext) -> int | None:
    """Stage entry point (see the module docstring)."""
    paths = paths_for(ctx)
    apaths = analysis_paths(ctx)
    return run_stage(
        ctx,
        STAGE,
        inputs=[
            paths.lexicon_complete,
            paths.cells,
            apaths.derivations,
            apaths.lexeme_strata,
        ],
        params={"version": VERSION},
        outputs=[
            apaths.noevidence,
            apaths.noevidence_pairs,
            apaths.noevidence_words,
        ],
        body=lambda: _noevidence(ctx, paths, apaths),
    )


def load(
    paths: CzechPaths, apaths: AnalysisPaths
) -> tuple[dict[str, roots.RootRecord], dict[str, str], list[dict[str, str]]]:
    """Noun root records with cells, lexeme keys and the pair rows."""
    strata = roots.load_strata(apaths.lexeme_strata)
    records, lexeme_key = roots.load_records(paths.lexicon_complete, strata)
    roots.add_cells(records, lexeme_key, paths.cells)
    pairs = list(jsonio.iter_csv(apaths.derivations))
    return records, lexeme_key, pairs


def _noevidence(
    ctx: RunContext, paths: CzechPaths, apaths: AnalysisPaths
) -> int:
    records, lexeme_key, pairs = load(paths, apaths)
    text = render(records, lexeme_key, pairs)
    ctx.write_tracked_text(
        apaths.noevidence, text, sources=REPORT_SOURCES, nc_aggregates=True
    )
    ctx.write_text(
        apaths.noevidence_words,
        counterexample_words(records, lexeme_key, pairs),
    )
    csvio.write_rows(
        apaths.noevidence_pairs,
        pair_rows(records, lexeme_key, pairs),
        PAIR_FIELDS,
        ctx=ctx,
    )
    print(text, flush=True)
    return 0


def pair_rows(
    records: Mapping[str, roots.RootRecord],
    lexeme_key: Mapping[str, str],
    pairs: Sequence[Mapping[str, str]],
) -> list[dict[str, Any]]:
    """NC rows: every P-suffix pair of a velar or r root with its base
    classes at every τ."""
    out: list[dict[str, Any]] = []
    for row in pairs:
        key = lexeme_key.get(row["base_lexeme_id"])
        if key is None:
            continue
        record = records[key]
        if (
            record.target_class not in TARGETS
            or row["suffix_class"] != sx.P_CLASS
            or row["outcome"] not in (sx.PAL, sx.EXCEPTION)
        ):
            continue
        out.append(
            {
                "root": key,
                "base_lemma": record.lemma,
                "target_class": record.target_class,
                "genders": record.genders(),
                "base_stratum_wikt": row["base_stratum_wikt"],
                "base_stratum": row["base_stratum"],
                **{
                    f"base_class_{tau}": base_class(record, tau)
                    for tau in TAUS
                },
                "derivative": row["derivative"],
                "suffix": row["suffix"],
                "outcome": row["outcome"],
                "basis": row["basis"],
                "derivative_stratum_wikt": row["derivative_stratum_wikt"],
                "derivative_stratum": row["derivative_stratum"],
                "count": row["count"],
                **{
                    f"attested_{tau}": row[f"attested_{tau}"]
                    for tau in CORPUS_TAUS
                },
                "source": csvio.split_list(row["source"]),
                "licence_class": "by-nc-sa",
            }
        )
    out.sort(key=lambda r: (r["root"], r["suffix"], r["derivative"]))
    return out


def render(
    records: Mapping[str, roots.RootRecord],
    lexeme_key: Mapping[str, str],
    pairs: Sequence[Mapping[str, str]],
) -> str:
    """The tracked report (module docstring)."""
    lines = [
        "Czech roots without alternation evidence x derivative palatalization (R21)",
        "=" * 72,
        "",
        "Base class from the corpus at tau (E0 inflection): alternating (an",
        "attested P form with the change) / non-alternating (only attested",
        "competing N allomorphs -u, -ove/-e, -ach, or a P ending without the",
        "change) / no evidence. Units: base x P suffix that palatalizes the",
        "base's target; outcome pal / plain / both. Derivatives attested at",
        "tau (summed tokens of their forms). Basis wikt = Wiktionary links",
        "and headword candidates. Strata: Wiktionary (R22). Words without",
        "counts (R9).",
    ]
    for tau in REPORT_TAUS:
        lines += ["", f"tau = {tau}", "-" * 72]
        for target in TARGETS:
            classes = Counter(
                base_class(r, tau)
                for r in records.values()
                if r.pos == roots.NOUN and r.target_class == target
            )
            lines.append(
                f"  {target} noun roots by base class: "
                + ", ".join(f"{c} {classes[c]:,}" for c in BASE_CLASSES)
            )
        main = units(records, lexeme_key, pairs, tau, "wikt")
        for target in TARGETS:
            chosen = [
                u for u in main if records[u.root].target_class == target
            ]
            lines.append(
                f"  target {target} (basis wikt), base x suffix units"
            )
            lines += suffix_tables(chosen)
            lines.append(
                f"  target {target}: by base stratum (all P suffixes)"
            )
            lines += stratum_tables(chosen)
            internal = [
                u
                for u in units(
                    records,
                    lexeme_key,
                    pairs,
                    tau,
                    "wikt",
                    derivative_strata=frozenset({INTERNAL}),
                )
                if records[u.root].target_class == target
            ]
            lines.append(f"  target {target}: {INTERNAL_DERIVATIVES}")
            lines += _render_table(crosstab(internal), "    ")
        lines += _sensitivity(records, lexeme_key, pairs, tau)
        neuters = neuter_velars(records, tau)
        lines.append(
            "  velar neuters on -ko/-go/-ho/-cho by base class: "
            + ", ".join(f"{c} {neuters[c]:,}" for c in BASE_CLASSES)
        )
        neuter_units = [
            u
            for u in main
            if records[u.root].readings == ["n"]
            and records[u.root].lemma.endswith(NEUTER_ENDINGS)
        ]
        lines.append("  velar neuters, all P suffixes:")
        lines += _render_table(crosstab(neuter_units), "    ")
        if tau == "dictionary":
            lines += _counterexample_lines(main, tau)
    lines += _nc_section(records, lexeme_key, pairs)
    return "\n".join(lines) + "\n"


def counterexample_words(
    records: Mapping[str, roots.RootRecord],
    lexeme_key: Mapping[str, str],
    pairs: Sequence[Mapping[str, str]],
) -> str:
    """The counterexample words at τ ≥ 1 token (gitignored: filtered by
    per-word corpus attestation, R9)."""
    lines = [
        "Czech R21 counterexamples at tau >= 1 token (untracked: the lists",
        "are filtered by per-word Leipzig attestation, R9)",
    ]
    for tau in ("ge1", "ge1_strict"):
        main = units(records, lexeme_key, pairs, tau, "wikt")
        lines += _counterexample_lines(main, tau)
    return "\n".join(lines) + "\n"


def _sensitivity(
    records: Mapping[str, roots.RootRecord],
    lexeme_key: Mapping[str, str],
    pairs: Sequence[Mapping[str, str]],
    tau: str,
) -> list[str]:
    lines = ["  sensitivity (velar, all P suffixes):"]
    for basis in ("linked", "wikt+corpus"):
        items = [
            u
            for u in units(records, lexeme_key, pairs, tau, basis)
            if records[u.root].target_class == ortho.VELAR
        ]
        lines.append(f"    basis {basis}")
        lines += _render_table(crosstab(items), "      ")
    return lines


def _counterexample_lines(items: Sequence[Unit], tau: str) -> list[str]:
    lines = [
        "",
        f"  Counterexamples at tau = {tau} (basis wikt): bases without",
        "  alternation evidence that have a palatalized derivative, by",
        "  (base stratum, derivative stratum);",
        "  * = learned -ik(a) → -ický adjective:",
    ]
    for (base, deriv), words in counterexamples(items).items():
        shown = words[:LIST_LIMIT]
        more = (
            f" … (+{len(words) - LIST_LIMIT})"
            if len(words) > LIST_LIMIT
            else ""
        )
        lines.append(f"    base {base} / derivative {deriv} ({len(words)}):")
        lines.append("      " + ", ".join(shown) + more)
    return lines


def _nc_section(
    records: Mapping[str, roots.RootRecord],
    lexeme_key: Mapping[str, str],
    pairs: Sequence[Mapping[str, str]],
) -> list[str]:
    lines = [
        "",
        "NC sources (aggregates only): basis all (DeriNet-only pairs added),",
        "base strata Wiktionary > czetyl > DeriNet (R22)",
        "-" * 72,
    ]
    for tau in ("dictionary", "ge1"):
        items = [
            u
            for u in units(
                records,
                lexeme_key,
                pairs,
                tau,
                "all",
                stratum_field="base_stratum",
                derivative_stratum_field="derivative_stratum",
            )
            if records[u.root].target_class == ortho.VELAR
        ]
        lines.append(f"  tau = {tau}, velar, all P suffixes")
        lines += _render_table(crosstab(items), "    ")
        lines += stratum_tables(items)
    return lines
