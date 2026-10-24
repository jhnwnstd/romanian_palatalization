"""Stage ``spellout``: spell-out TP on certain-specification roots (R20).

Gorman & Winstead's spell-out account gives a non-palatalizing
allomorph only to roots prespecified as inalterable. A root whose
masculine plural is attested with -ci/-gi is underspecified with
certainty, so the account predicts the palatalizing member of every
competing pair for it. The spell-out TP counts, per pair set:

    roots   certain roots: dorsal roots whose -i cell has the status
            ``alternates`` at τ (E0, dominant coding: an attested
            palatalized -i plural, R19 doublets resolved by tokens)
    N       certain roots with an attested derivative of the pair
            (``italian_derivations.csv``: accepted links; derivative
            attested at τ, dictionary-listed at τ = dictionary)
    e       roots that choose the non-palatalizing member only;
            roots with both members are doublets coded by tokens (R19)

Pair sets (R20: "Italian: -izzare vs -eggiare/-are/-ire; palatalized
vs velar -ista/-ismo/-ità"):

    verbalizer          -izzare vs -eggiare, -ire, -are: the
                        palatalizing member is a palatalized -izzare
                        (tipicizzare); a velar -izzare (antichizzare),
                        -eggiare, -ire or -are verb is the other member
    verbalizer_no_are   the same without -are (the default conversion
                        verbalizer of Italian, mostly Latin borrowings
                        on these stems: medicare, caricare)
    ista_ismo_ita       palatalized vs velar -ista, -ismo, -ità
    learned             palatalized vs velar, every suffix of the
                        learned set (:data:`.suffixes.LEARNED_SET`)

It is computed with :func:`lexpipe.tp.two_outcome_rows`: in the cell
named after the pair set every derivative of the pair is a ``P``
observation whose ``alternates`` means "palatalized" (``palatalize`` =
the spell-out default, ``retain`` = choosing the other member).
Subsets: ``all`` (named links, the three doublet codings),
``all+candidates`` (adds the generated candidates of R21),
``all+candidates+other_base`` (adds the candidates the palatal-base
rule rejected: :data:`OTHER_BASE_TIER`),
``base:<stratum>`` (R22, the base root's stratum),
``deriv:internal formation`` (only derivatives formed in Italian) and
``contrast:nonalt`` (roots whose -i plural is attested with -chi/-ghi:
the prespecified side, for comparison). Rows go to
``italian_tp_long.csv`` (table ``spellout``, built by the ``tp``
stage) and are summarized in ``analysis/italian/italian_spellout.txt``
(aggregates and example words; no per-word counts, R9).
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

from lexpipe import attest, csvio, tp
from lexpipe.context import RunContext, run_stage
from lexpipe.langs.italian import derivations, etymology, features, roots
from lexpipe.langs.italian import suffixes as sx
from lexpipe.langs.italian.layout import Layout
from lexpipe.langs.italian.tptables import cell, finish

STAGE: Final[str] = "spellout"
# 2: gated report writer; shared stratum classifier upstream.
# 3: the truncated derivatives and the fragility reading in the report;
#    the candidate tier with the palatal-base rejections.
# 4: the all+candidates+other_base subset.
VERSION: Final[int] = 4
SOURCE: Final[str] = "wiktionary_en"
TABLE: Final[str] = "spellout"
VERBALIZERS: Final[tuple[str, ...]] = ("-izzare", "-eggiare", "-ire", "-are")
PAIR_SETS: Final[Mapping[str, tuple[str, ...]]] = {
    "verbalizer": VERBALIZERS,
    "verbalizer_no_are": ("-izzare", "-eggiare", "-ire"),
    "ista_ismo_ita": ("-ista", "-ismo", "-ità"),
    "learned": sx.LEARNED_SET,
}
# The tier of a candidate rejected only because a frontstem lemma on its
# stem could explain its palatal spelling (calcio → calcista): a rule that
# can only remove palatal candidates, so the candidate tier is reported
# with these links (``…+other_base``) and without them (R16).
OTHER_BASE_TIER: Final[str] = "candidate:other_base"
TIER_SUBSETS: Final[Mapping[str, tuple[str, ...]]] = {
    "all": (derivations.TIER_NAMED,),
    "all+candidates": (derivations.TIER_NAMED, derivations.TIER_CANDIDATE),
    "all+candidates+other_base": (
        derivations.TIER_NAMED,
        derivations.TIER_CANDIDATE,
        OTHER_BASE_TIER,
    ),
}
INTERNAL: Final[str] = etymology.INTERNAL
EXAMPLES: Final[int] = 12


@dataclass(frozen=True, slots=True)
class Derivative:
    """An accepted derivative link, reduced to what R20–R23 need."""

    base_root_id: str
    word: str
    key: str
    suffix: str
    group: str
    palatalized: bool | None
    edge: str
    tier: str
    deriv_stratum: str
    base_stratum: str
    count_cs: int
    attested: Mapping[str, bool]
    base_lexeme_id: str
    base_word: str
    base_pos: str
    base_gender: str
    masc_counterpart: str
    ic_suffix: bool
    fact_check: str
    # The base's determiner/pronoun plural (pochi of poco): its reading
    # and the thresholds it reaches (R21).
    base_det_plural_palatal: bool | None = None
    base_det_plural_taus: frozenset[str] = frozenset()


def load_derivatives(path: Path) -> list[Derivative]:
    """Accepted links with a readable edge (palatal, velar, back), and
    the candidates rejected by the palatal-base rule, with the tier
    :data:`OTHER_BASE_TIER` (only the ``…+other_base`` sensitivity sets
    read them)."""
    out: list[Derivative] = []
    for row in csvio.read_rows(path):
        status = row["link_status"]
        if status == derivations.ACCEPTED:
            tier = row["link_tier"]
        elif status == derivations.REJECTED_PALATAL_BASE:
            tier = OTHER_BASE_TIER
        else:
            continue
        if row["edge"] not in ("palatal", "velar", "back"):
            continue
        out.append(
            Derivative(
                base_root_id=row["base_root_id"],
                word=row["derivative"],
                key=row["derivative"].lower(),
                suffix=row["suffix"],
                group=row["suffix_group"],
                palatalized=csvio.parse_bool(row["palatalized"]),
                edge=row["edge"],
                tier=tier,
                deriv_stratum=row["deriv_stratum"],
                base_stratum=row["base_stratum"],
                count_cs=int(row["count_cs"] or 0),
                attested={
                    tau: csvio.parse_bool(row[f"attested_{tau}"]) is True
                    for tau in attest.TAUS
                },
                base_lexeme_id=row["base_lexeme_id"],
                base_word=row["base_word"],
                base_pos=row["base_pos"],
                base_gender=row["base_gender"],
                masc_counterpart=row["masc_counterpart"],
                ic_suffix=csvio.parse_bool(row["ic_suffix"]) is True,
                fact_check=row["fact_check"],
                base_det_plural_palatal=csvio.parse_bool(
                    row["base_det_plural_palatal"]
                ),
                base_det_plural_taus=frozenset(
                    csvio.split_list(row["base_det_plural_taus"])
                ),
            )
        )
    return out


def by_root(items: Iterable[Derivative]) -> dict[str, list[Derivative]]:
    out: dict[str, list[Derivative]] = defaultdict(list)
    for item in items:
        out[item.base_root_id].append(item)
    return dict(out)


def chose_palatal(item: Derivative) -> bool:
    """The palatalizing member: a palatalized derivative."""
    return item.palatalized is True


@dataclass(frozen=True, slots=True)
class Filter:
    """Which derivatives a spell-out subset counts."""

    suffixes: frozenset[str]
    tiers: frozenset[str]
    internal_only: bool = False

    def keep(self, item: Derivative) -> bool:
        return (
            item.suffix in self.suffixes
            and item.tier in self.tiers
            and (not self.internal_only or item.deriv_stratum == INTERNAL)
        )


def spellout_root(
    record: roots.RootRecord,
    items: Sequence[Derivative],
    chosen: Filter,
    pair: str,
) -> tp.RootObservation:
    """The root with one cell, named after the pair set, holding its
    pair derivatives."""
    observations = tuple(
        tp.Observation(
            alternates=chose_palatal(item),
            attested_by_tau=item.attested,
            allomorph="P",
            form=item.key,
            count=item.count_cs,
        )
        for item in items
        if chosen.keep(item)
    )
    return tp.RootObservation(
        record.root_id,
        record.target_class,
        float(record.freq_paradigm),
        {pair: observations},
    )


@dataclass(frozen=True, slots=True)
class SubsetSpec:
    """A spell-out subset: roots of one -i status (and stratum), which
    derivatives count, which doublet codings are reported."""

    name: str
    status: tp.RootStatus
    tiers: tuple[str, ...]
    stratum: str = ""
    internal_only: bool = False
    codings: tuple[str, ...] = ("dominant",)


SUBSET_SPECS: Final[tuple[SubsetSpec, ...]] = (
    SubsetSpec(
        "all",
        tp.RootStatus.ALTERNATES,
        TIER_SUBSETS["all"],
        codings=tp.CODINGS,
    ),
    SubsetSpec(
        "all+candidates",
        tp.RootStatus.ALTERNATES,
        TIER_SUBSETS["all+candidates"],
        codings=tp.CODINGS,
    ),
    SubsetSpec(
        "all+candidates+other_base",
        tp.RootStatus.ALTERNATES,
        TIER_SUBSETS["all+candidates+other_base"],
    ),
    *(
        SubsetSpec(
            f"base:{name}",
            tp.RootStatus.ALTERNATES,
            TIER_SUBSETS["all"],
            stratum=name,
        )
        for name in etymology.STRATA
    ),
    SubsetSpec(
        f"deriv:{INTERNAL}",
        tp.RootStatus.ALTERNATES,
        TIER_SUBSETS["all"],
        internal_only=True,
    ),
    SubsetSpec("contrast:nonalt", tp.RootStatus.NONALT, TIER_SUBSETS["all"]),
)


def root_statuses(
    records: roots.Records,
) -> dict[tuple[str, str], tp.RootStatus]:
    """``{(root_id, tau): -i status}`` of every dorsal root (E0)."""
    out: dict[tuple[str, str], tp.RootStatus] = {}
    for record in records.roots.values():
        if not record.dorsal:
            continue
        found = roots.observation(record, "E0")
        for tau in attest.TAUS:
            out[(record.root_id, tau)] = tp.root_status(
                found,
                "-i",
                tau,
                "dominant",
                min_tokens=features.DOUBLET_MIN_TOKENS,
            )
    return out


def subset_roots(
    records: roots.Records,
    statuses: Mapping[tuple[str, str], tp.RootStatus],
    spec: SubsetSpec,
    tau: str,
) -> list[roots.RootRecord]:
    """Dorsal roots of a subset at τ."""
    return [
        record
        for record in records.roots.values()
        if record.dorsal
        and statuses[(record.root_id, tau)] == spec.status
        and (not spec.stratum or record.stratum == spec.stratum)
    ]


def spellout_rows(
    records: roots.Records, items: Sequence[Derivative]
) -> list[dict[str, Any]]:
    """Every spell-out row (module docstring)."""
    grouped = by_root(items)
    statuses = root_statuses(records)
    rows: list[dict[str, Any]] = []
    for pair, suffix_names in PAIR_SETS.items():
        for spec in SUBSET_SPECS:
            chosen = Filter(
                frozenset(suffix_names),
                frozenset(spec.tiers),
                spec.internal_only,
            )
            for tau in attest.TAUS:
                observed = [
                    spellout_root(
                        record, grouped.get(record.root_id, ()), chosen, pair
                    )
                    for record in subset_roots(records, statuses, spec, tau)
                ]
                for coding in spec.codings:
                    spec_cell = cell(
                        TABLE, "dorsal", pair, tau=tau, coding=coding
                    )
                    rows += [
                        finish(row)
                        for row in tp.two_outcome_rows(
                            observed,
                            spec_cell,
                            subset=spec.name,
                            memo=f"pair set {pair}",
                        )
                    ]
    return rows


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


def _fmt(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.1f}"
    return str(value)


def _table(rows: Sequence[Mapping[str, Any]]) -> list[str]:
    head = (
        f"  {'pair set':<18} {'subset':<26} {'tau':<10} {'coding':<18} "
        f"{'N_target':>8} {'N':>5} {'e':>5} {'theta':>6} {'tol?':>4} "
        f"{'frag':>4}"
    )
    lines = [head]
    for row in rows:
        if row["rule"] != "palatalize":
            continue
        lines.append(
            f"  {row['trigger']:<18} {row['subset']:<26} {row['tau']:<10} "
            f"{row['coding']:<18} {row['N_target']:>8} {row['N_trigger']:>5} "
            f"{row['e']:>5} {_fmt(row['theta_trigger']):>6} "
            f"{_fmt(row['tolerated_trigger']):>4} {_fmt(row['fragility']):>4}"
        )
    return lines


def _examples(
    records: roots.Records, items: Sequence[Derivative], tau: str
) -> list[str]:
    grouped = by_root(items)
    statuses = root_statuses(records)
    choice: dict[str, list[str]] = {"P": [], "N": [], "both": []}
    chosen = Filter(frozenset(VERBALIZERS), frozenset(TIER_SUBSETS["all"]))
    for record in records.roots.values():
        if (
            not record.dorsal
            or statuses[(record.root_id, tau)] != tp.RootStatus.ALTERNATES
        ):
            continue
        found = [
            d
            for d in grouped.get(record.root_id, ())
            if chosen.keep(d) and d.attested[tau]
        ]
        if not found:
            continue
        kinds = {chose_palatal(d) for d in found}
        key = "both" if len(kinds) == 2 else ("P" if True in kinds else "N")
        words = sorted({d.word for d in found})
        choice[key].append(f"{record.lemma} ({', '.join(words)})")
    out: list[str] = []
    for key, label in (
        ("P", "palatalized -izzare only"),
        ("N", "non-palatalizing member only"),
        ("both", "both"),
    ):
        listed = sorted(choice[key])
        shown = "; ".join(listed[:EXAMPLES])
        more = (
            f"; … (+{len(listed) - EXAMPLES})"
            if len(listed) > EXAMPLES
            else ""
        )
        out.append(f"  {label} [{len(listed)} roots]: {shown}{more}")
    return out


def truncated_lines(records: roots.Records, path: Path) -> list[str]:
    """Derivatives the spell-out cannot read: an accepted named link of
    a certain root (τ = dictionary) whose stem edge is truncated (the
    -ic- of tipico dropped in tipizzare): neither member of the pair
    shows, so the root is not counted for that derivative."""
    statuses = root_statuses(records)
    certain = {
        root_id
        for (root_id, tau), status in statuses.items()
        if tau == "dictionary" and status == tp.RootStatus.ALTERNATES
    }
    out: list[str] = []
    for pair, names in PAIR_SETS.items():
        found = [
            row
            for row in csvio.read_rows(path)
            if row["link_status"] == derivations.ACCEPTED
            and row["link_tier"] == derivations.TIER_NAMED
            and row["edge"] == "truncated"
            and row["suffix"] in names
            and row["base_root_id"] in certain
        ]
        if not found:
            continue
        roots_found = {row["base_root_id"] for row in found}
        examples = sorted({row["derivative"] for row in found})[:EXAMPLES]
        out.append(
            f"  {pair}: {len(found)} truncated derivatives of "
            f"{len(roots_found)} certain roots not counted "
            f"({', '.join(examples)}{', …' if len(found) > EXAMPLES else ''})"
        )
    return out


def report(
    rows: Sequence[Mapping[str, Any]],
    records: roots.Records,
    items: Sequence[Derivative],
    truncated: Sequence[str] = (),
) -> str:
    """The text report (aggregates and example words only, R9);
    *truncated* are :func:`truncated_lines`."""
    lines = [
        "Italian spell-out TP on certain-specification roots (design R20)",
        "=" * 72,
        "",
        "Certain roots: dorsal roots whose -i plural is attested with",
        "-ci/-gi at tau (E0, dominant coding). N = certain roots with an",
        "attested derivative of the pair set; e = roots choosing only the",
        "non-palatalizing member; roots with both members are doublets",
        "coded by tokens (R19). theta = N/ln N. 'all' counts named links",
        "(etymology templates, derived lists, head forms); 'all+candidates'",
        "adds generated candidates (R21). Derivatives: dictionary-listed at",
        "tau = dictionary, corpus tokens (summed over their corpus forms)",
        "at the other thresholds. A verdict with fragility 1 flips with one",
        "more (or one fewer) exception: read it as a tie.",
        "",
    ]
    for pair in PAIR_SETS:
        lines.append(f"Pair set {pair}: {', '.join(PAIR_SETS[pair])}")
        lines.append("-" * 72)
        lines += _table([r for r in rows if r["trigger"] == pair])
        lines.append("")
    lines.append("Examples (verbalizer pair, tau = dictionary; words only)")
    lines.append("-" * 72)
    lines += _examples(records, items, "dictionary")
    lines.append("")
    if truncated:
        lines.append(
            "Not counted: truncated stem edges (tau = dictionary, named links)"
        )
        lines.append("-" * 72)
        lines += list(truncated)
        lines.append("")
    lines.append(
        "Learners and fragility: fragility = fewest item changes that flip"
    )
    lines.append(
        "the verdict (lexpipe.learners.fragility); rows are in "
        "italian_tp_long.csv (table spellout)."
    )
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# The stage
# ---------------------------------------------------------------------------


def load_records(layout: Layout) -> roots.Records:
    """The root records with E0–E2 realizations (named derivative tier)."""
    return roots.load_all(
        layout.lexicon_complete,
        layout.cells,
        layout.verbs,
        layout.derivations,
    )


def _body(ctx: RunContext, layout: Layout) -> int:
    records = load_records(layout)
    items = load_derivatives(layout.derivations)
    rows = spellout_rows(records, items)
    truncated = truncated_lines(records, layout.derivations)
    ctx.write_tracked_text(
        layout.spellout,
        report(rows, records, items, truncated),
        sources=SOURCE,
    )
    print(f"  {len(rows):,} spell-out rows")
    return 0


def inputs(layout: Layout) -> list[Path]:
    return [
        layout.lexicon_complete,
        layout.cells,
        layout.verbs,
        layout.derivations,
    ]


def outputs(layout: Layout) -> list[Path]:
    return [layout.spellout]


def run(ctx: RunContext) -> int | None:
    """Run the ``spellout`` stage (skipped when its manifest matches)."""
    layout = Layout(ctx)
    return run_stage(
        ctx,
        STAGE,
        inputs=inputs(layout),
        params={"version": VERSION, "pair_sets": dict(PAIR_SETS)},
        outputs=outputs(layout),
        body=lambda: _body(ctx, layout),
    )
