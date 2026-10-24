"""Stage ``docs``: data dictionary and the README's generated results.

1. Checks the schema (``docs/italian/schema.toml``) against every
   writer's header constant and the files on disk
   (:func:`.schema.header_problems`); a difference stops the stage.
2. Writes ``docs/italian/data_dict.md`` from the schema.
3. Replaces the block between :data:`BEGIN` and :data:`END` in
   ``docs/italian/README.md`` with the headline numbers of the analyses
   (from ``italian_tp_long.csv``, ``italian_learners.csv`` and the
   co-variation and no-evidence units recomputed from the data), so the
   README never drifts from the outputs. The hand-written text around
   the block is kept; a README without the markers gets the block
   appended (a fixture tree gets a stub).

Aggregates only (R9): both files go through the licence gate of
:meth:`lexpipe.context.RunContext.write_tracked_text`.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, Final

from lexpipe import csvio, stats
from lexpipe.context import RunContext, run_stage
from lexpipe.langs.italian import covariation, noevidence, schema, spellout
from lexpipe.langs.italian import suffixes as sx
from lexpipe.langs.italian.layout import Layout

STAGE: Final[str] = "docs"
# 2: gated writers; the schema of the lean raw lexicon and the shared
# R27 columns.
# 3: R26 suffix-once rows, typed-plural rows, -ic- bases apart.
# 4: the no-evidence table's other_base tier.
VERSION: Final[int] = 4
# What the data dictionary describes and what the results rest on.
DICT_SOURCES: Final[tuple[str, ...]] = (
    "wiktionary_en",
    "wiktionary_it",
    "morphit",
)
RESULTS_SOURCES: Final[tuple[str, ...]] = ("wiktionary_en",)
BEGIN: Final[str] = (
    "<!-- BEGIN generated results: written by the docs stage "
    "(lexpipe.langs.italian.docstage); edit the code, not this block -->"
)
END: Final[str] = "<!-- END generated results -->"
STUB: Final[str] = "# Italian palatalization data\n\n"

Row = Mapping[str, str]


def _find(rows: Sequence[Row], **where: str) -> Row | None:
    defaults = {"evidence_set": "E0", "coding": "dominant", "V": ""}
    defaults.update(where)
    for row in rows:
        if all(row.get(k, "") == v for k, v in defaults.items()):
            return row
    return None


def _theta(row: Row) -> str:
    return (
        f"{float(row['theta_trigger']):.1f}" if row["theta_trigger"] else "-"
    )


def _verdict(pal: Row, ret: Row) -> str:
    good = (
        pal["tolerated_trigger"] == "True",
        ret["tolerated_trigger"] == "True",
    )
    if pal["tolerated_trigger"] == "":
        return "empty"
    return {
        (True, True): "both",
        (True, False): "palatalize",
        (False, True): "retain",
        (False, False): "neither",
    }[good]


def tp_table(
    rows: Sequence[Row],
    items: Sequence[tuple[str, dict[str, str]]],
    taus: Sequence[str] = ("ge1", "dictionary"),
) -> list[str]:
    """Markdown rows: label, τ, N_target, N, e_pal, e_ret, θ, verdict
    (τ ≥ 1 first: the dictionary rows mostly restate Module:it-common's
    defaults)."""
    out = [
        "| row | τ | N_target | N | e (palatalize) | e (retain) | θ | productive |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for label, where in items:
        for tau in (where["tau"],) if "tau" in where else taus:
            found = {k: v for k, v in where.items() if k != "tau"}
            pal = _find(rows, rule="palatalize", tau=tau, **found)
            ret = _find(rows, rule="retain", tau=tau, **found)
            if pal is None or ret is None:
                continue
            out.append(
                f"| {label} | {tau} | {pal['N_target']} | {pal['N_trigger']} "
                f"| {pal['e']} | {ret['e']} | {_theta(pal)} | "
                f"{_verdict(pal, ret)} |"
            )
    return out


R17_ITEMS: Final[tuple[tuple[str, dict[str, str]], ...]] = (
    (
        "R17a dorsal/_[+front V], -i/-e",
        {
            "table": "R17a",
            "target_class": "dorsal",
            "trigger": "-i/-e",
            "subset": "all",
        },
    ),
    (
        "R17b dorsal/_[+high,+front], -i",
        {
            "table": "R17b",
            "target_class": "dorsal",
            "trigger": "-i",
            "subset": "all",
        },
    ),
    (
        "R17b -i, nouns only",
        {
            "table": "R17b",
            "target_class": "dorsal",
            "trigger": "-i",
            "subset": "nouns",
        },
    ),
    (
        "R17b -i, core (no -ico, no -logo/-fago/-urgo)",
        {
            "table": "R17b",
            "target_class": "dorsal",
            "trigger": "-i",
            "subset": "core",
        },
    ),
    (
        "R17b non-trigger -e (feminine plurals)",
        {
            "table": "R17b",
            "target_class": "dorsal",
            "trigger": "-e",
            "subset": "all",
        },
    ),
    (
        "R17c coronal/_[+high,+front], -i",
        {
            "table": "R17c",
            "target_class": "coronal",
            "trigger": "-i",
            "subset": "all",
        },
    ),
)
FIT_ITEMS: Final[tuple[tuple[str, dict[str, str]], ...]] = (
    (
        "velar post-tonic (paroxytone type)",
        {
            "table": "it_stress",
            "target_class": "dorsal",
            "trigger": "-i",
            "subset": "stress:post_tonic",
        },
    ),
    (
        "velar not post-tonic (proparoxytone type)",
        {
            "table": "it_stress",
            "target_class": "dorsal",
            "trigger": "-i",
            "subset": "stress:not_post_tonic",
        },
    ),
    (
        "velar post-tonic, typed plurals only",
        {
            "table": "it_stress",
            "target_class": "dorsal",
            "trigger": "-i",
            "subset": "stress:post_tonic",
            "tau": "dictionary_explicit",
        },
    ),
    (
        "velar not post-tonic, typed plurals only",
        {
            "table": "it_stress",
            "target_class": "dorsal",
            "trigger": "-i",
            "subset": "stress:not_post_tonic",
            "tau": "dictionary_explicit",
        },
    ),
    (
        "-ico suffix targets (R26)",
        {
            "table": "it_suffix",
            "target_class": "dorsal",
            "trigger": "-i",
            "subset": "suffix:-ico",
        },
    ),
    (
        "-logo person nouns",
        {
            "table": "logtable",
            "target_class": "dorsal",
            "trigger": "-i",
            "subset": "-logo:person",
        },
    ),
    (
        "-logo thing nouns",
        {
            "table": "logtable",
            "target_class": "dorsal",
            "trigger": "-i",
            "subset": "-logo:thing",
        },
    ),
    (
        "deverbal from -are (carico type)",
        {
            "table": "it_suffix",
            "target_class": "dorsal",
            "trigger": "-i",
            "subset": "deverbal_from_are",
        },
    ),
    (
        "verbs -are, 2sg/1pl/sbjv",
        {
            "table": "it_verbs",
            "target_class": "dorsal",
            "trigger": "-i",
            "subset": "conj:are",
            "evidence_set": "E1",
        },
    ),
    (
        "verbs root-stressed -ere",
        {
            "table": "it_verbs",
            "target_class": "dorsal",
            "trigger": "-i",
            "subset": "conj:ere_root",
            "evidence_set": "E1",
        },
    ),
    (
        "verbs -isc-",
        {
            "table": "it_verbs",
            "target_class": "dorsal",
            "trigger": "-i",
            "subset": "conj:isc",
            "evidence_set": "E1",
        },
    ),
    (
        "pooled noun + verb",
        {
            "table": "it_nounverb",
            "target_class": "dorsal",
            "trigger": "-i",
            "subset": "pooled",
            "evidence_set": "E1",
        },
    ),
    (
        "spell-out: -izzare vs -eggiare/-ire/-are",
        {
            "table": "spellout",
            "target_class": "dorsal",
            "trigger": "verbalizer",
            "subset": "all",
        },
    ),
    (
        "spell-out: same without -are",
        {
            "table": "spellout",
            "target_class": "dorsal",
            "trigger": "verbalizer_no_are",
            "subset": "all",
        },
    ),
    (
        "spell-out: palatal vs velar -ista/-ismo/-ità",
        {
            "table": "spellout",
            "target_class": "dorsal",
            "trigger": "ista_ismo_ita",
            "subset": "all",
        },
    ),
)


def nested_table(rows: Sequence[Row]) -> list[str]:
    """The nested sets that carry the model-fit verdicts."""
    out = [
        "| nested set | τ | rule (subset) | N | e | θ | tolerated |",
        "|---|---|---|---|---|---|---|",
    ]
    wanted = (
        ("it_stress", "dorsal", "nested set giavazzi", "E0"),
        ("it_stress", "dorsal", "nested set suffix_once", "E0"),
        ("it_stress", "dorsal", "nested set orth_ico_once", "E0"),
        ("it_stress", "g", "nested set giavazzi", "E0"),
        ("it_nounverb", "dorsal", "nested set category_stress", "E1"),
    )
    for table, target, memo, evidence in wanted:
        for tau in ("ge1", "dictionary", "dictionary_explicit"):
            for row in rows:
                if (
                    row["table"] == table
                    and row["target_class"] == target
                    and row["tau"] == tau
                    and row["evidence_set"] == evidence
                    and row["trigger"] == "-i"
                    and row["memo"].endswith(memo)
                ):
                    out.append(
                        f"| {table} {target} {memo.split()[-1]} | {tau} | "
                        f"{row['rule']} ({row['subset']}) | {row['N_trigger']} "
                        f"| {row['e']} | {_theta(row)} | "
                        f"{row['tolerated_trigger']} |"
                    )
    return out


def learner_table(rows: Sequence[Row]) -> list[str]:
    """Learner shares at V = 1000 and 5000 (R17 cells at τ =
    dictionary; the R17b -i and it_stress cells at both τ)."""
    out = [
        "| cell | subset | E | τ | V | palatalize | retain | neither "
        "| conditioned |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for row in rows:
        both_taus = row["trigger"] == "-i" and row["table"] in (
            "R17b",
            "it_stress",
        )
        if row["V"] not in ("1000", "5000") or (
            row["tau"] != "dictionary" and not both_taus
        ):
            continue
        out.append(
            f"| {row['table']} {row['trigger']} | {row['subset']} | "
            f"{row['evidence_set']} | "
            f"{row['tau']} | {row['V']} | {float(row['share_palatalize']):.3f} "
            f"| {float(row['share_retain']):.3f} | "
            f"{float(row['share_neither']):.3f} | "
            f"{float(row['share_conditioned']):.3f} |"
        )
    return out


def covariation_lines(layout: Layout, records: Any, items: Any) -> list[str]:
    """Pooled 2×2 of the learned set (internal formations, named links)."""
    statuses = spellout.root_statuses(records)
    out = [
        "| τ | suffix set | P base: P / N | N base: P / N | both "
        "| Fisher p | MH OR (strata) |",
        "|---|---|---|---|---|---|---|",
    ]
    for tau in ("dictionary", "ge1"):
        for label, names, keep_ic in (
            ("learned", sx.LEARNED_SET, True),
            ("learned, no -ic- bases", sx.LEARNED_SET, False),
            ("native", sx.NATIVE_SET, True),
        ):
            units = [
                u
                for u in covariation.build_units(
                    records,
                    items,
                    statuses,
                    tau,
                    covariation.TIERS["named"],
                    names,
                )
                if u.stratum == covariation.INTERNAL
                and (keep_ic or not u.ic_suffix)
            ]
            (a, b), (c, d) = covariation.table_of(units)
            both = sum(1 for u in units if u.outcome == covariation.BOTH)
            fisher = (
                stats.fisher_exact(((a, b), (c, d)))
                if a + b + c + d
                else float("nan")
            )
            pooled = stats.mantel_haenszel(covariation.strata_tables(units))
            mh = (
                f"{pooled.estimate:.3g} ({pooled.n_strata})"
                if pooled is not None
                else "undefined"
            )
            out.append(
                f"| {tau} | {label} | {a} / {b} | {c} / {d} | {both} | "
                f"{fisher:.3g} | {mh} |"
            )
    return out


def noevidence_lines(records: Any, items: Any) -> list[str]:
    """Learned-set units by R21 base class."""
    out = [
        "| τ | links | alternating P / N | non-alternating P / N "
        "| no evidence P / N |",
        "|---|---|---|---|---|",
    ]
    for tau in ("dictionary", "ge1"):
        for tier, tiers in noevidence.TIERS.items():
            units = [
                u
                for u in noevidence.build_units(records, items, tau, tiers)
                if u.suffix in sx.LEARNED_SET
            ]
            counts = Counter((u.klass, u.outcome) for u in units)
            cells = " | ".join(
                f"{counts[(k, 'P')]} / {counts[(k, 'N')]}"
                for k in noevidence.CLASSES
            )
            out.append(f"| {tau} | {tier} | {cells} |")
    return out


def results_block(layout: Layout) -> str:
    """The generated README block (module docstring)."""
    rows = csvio.read_rows(layout.tp_long)
    learner_rows = csvio.read_rows(layout.learners)
    records = spellout.load_records(layout)
    items = spellout.load_derivatives(layout.derivations)
    lines = [
        BEGIN,
        "",
        "### Headline numbers (generated)",
        "",
        "Unit: root types (R4); E0 unless stated; doublets coded by the "
        "dominant variant (R19); θ = N/ln N. Full rows: "
        "`analysis/italian/italian_tp_long.csv`; digest: "
        "`analysis/italian/italian_tp_summary.txt`.",
        "",
        "**Cross-language rule family (R17).**",
        "",
        *tp_table(rows, R17_ITEMS),
        "",
        "**Italian model fits, suffix, verb and spell-out rows.**",
        "",
        *tp_table(rows, FIT_ITEMS),
        "",
        "**Nested rule sets (Yang's order).**",
        "",
        *nested_table(rows),
        "",
        "**Simulated learners (1000 vocabularies per V, R18).**",
        "",
        *learner_table(learner_rows),
        "",
        "**Co-variation (R23): base -i plural × learned derivative, "
        "Italian-internal formations, named links.**",
        "",
        *covariation_lines(layout, records, items),
        "",
        "**No-evidence roots (R21): learned-set units by base class.**",
        "",
        *noevidence_lines(records, items),
        "",
        END,
    ]
    return "\n".join(lines) + "\n"


def update_readme(text: str, block: str) -> str:
    """*text* with the generated block replaced (or appended)."""
    start = text.find(BEGIN)
    end = text.find(END)
    if start >= 0 and end > start:
        tail = text[end + len(END) :]
        return text[:start] + block.rstrip("\n") + tail
    return text.rstrip("\n") + "\n\n" + block


def _body(ctx: RunContext, layout: Layout) -> int:
    loaded = schema.load_schema(layout.schema)
    problems = schema.header_problems(loaded, ctx.root)
    if problems:
        for problem in problems:
            print(f"  schema: {problem}")
        return 1
    ctx.write_tracked_text(
        layout.data_dict,
        schema.render_data_dict(loaded),
        sources=DICT_SOURCES,
    )
    current = (
        layout.readme.read_text(encoding="utf-8")
        if layout.readme.is_file()
        else STUB
    )
    ctx.write_tracked_text(
        layout.readme,
        update_readme(current, results_block(layout)),
        sources=RESULTS_SOURCES,
    )
    print(
        f"  data dictionary: {len(loaded)} files, "
        f"{sum(len(f.columns) for f in loaded.values())} columns"
    )
    return 0


def inputs(layout: Layout) -> list[Path]:
    return [
        layout.schema,
        layout.tp_long,
        layout.learners,
        *spellout.inputs(layout),
    ]


def outputs(layout: Layout) -> list[Path]:
    return [layout.data_dict, layout.readme]


def run(ctx: RunContext) -> int | None:
    """Run the ``docs`` stage (skipped when its manifest matches)."""
    layout = Layout(ctx)
    return run_stage(
        ctx,
        STAGE,
        inputs=inputs(layout),
        params={"version": VERSION},
        outputs=outputs(layout),
        body=lambda: _body(ctx, layout),
    )
