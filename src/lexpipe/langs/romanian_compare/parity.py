"""Parity with the legacy Romanian TP tables (design §7, §8).

Two checks, both written to ``analysis/comparison/romanian_parity.txt``:

1. **Reproduction.** Every numeric row of
   ``analysis/romanian_tp_summary_full.csv`` and
   ``analysis/romanian_tp_summary_downsampled.csv`` (written by
   ``analysis/analyze_romanian_palatalization.R``) is recomputed from
   the lexicon's own columns with the R script's definitions
   (``tp_in_domain``, ``opportunity_tp``, cluster rows, the top-1000
   downsampling, the derivational rows and their orthographic
   detector) and passed through :func:`lexpipe.tp.tp_row`, the core
   function the Italian and Czech tables use. Every cell (N, mutated,
   non-mutated, rate, majority?, tolerated?) must match.

2. **Side by side.** Where a harmonized cell coincides with a legacy
   row (<c> + <-i> plural, <t> + <-e> plural, ...), the two counts are
   shown together, and the difference is decomposed root by root:
   each root counted by one definition and not the other gets the
   first applicable reason (the legacy exclusion it fell under, in
   the legacy code's order, or the harmonized structural NDE class).
   Legacy N + added − removed = harmonized N_trigger is checked.
"""

from __future__ import annotations

import csv
import math
import re
from collections import Counter
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

from lexpipe import tp
from lexpipe.langs.romanian_compare import lexicon as lx
from lexpipe.langs.romanian_compare import morphology as mo
from lexpipe.langs.romanian_compare.harmonized import (
    SUBSETS,
    ObservationCache,
    frame,
)
from lexpipe.langs.romanian_compare.records import RootRecord

SEGMENTS: Final[tuple[str, ...]] = ("c", "g", "t", "d", "s", "z")
NDEB: Final[tuple[str, ...]] = ("gimpe", "ochi", "paduchi")
NDEB_LABEL: Final[Mapping[str, str]] = {
    "gimpe": "gimpe type",
    "ochi": "ochi-ochi type",
    "paduchi": "paduche-paduchi type",
}
CLUSTERS: Final[tuple[str, ...]] = ("st", "sc", "ct")
FRONT_VERB_SUFFIXES: Final[tuple[str, ...]] = ("-i", "-iza", "-ifica")
Counts = tuple[int, int]  # (mutated, non-mutated)


# ---------------------------------------------------------------------------
# 1. Reproduction of the R tables
# ---------------------------------------------------------------------------


def legacy_nde(row: lx.LexiconRow) -> str:
    reason = row.exception_reason
    return reason[4:] if reason.startswith("nde:") else ""


def opportunity_tp(row: lx.LexiconRow) -> str:
    """The R script's ``opportunity_tp``."""
    nde = legacy_nde(row)
    ending = lx.plural_ending(row.plural)
    if nde in NDEB and row.opportunity == "none":
        if ending in ("i", "e"):
            return ending
        if ending == "uri":
            return "i"
    return row.opportunity


def cluster_simple(cluster: str) -> str:
    for name in CLUSTERS:
        if cluster.startswith(name):
            return name
    return ""


def _tally(
    rows: Iterable[lx.LexiconRow], key: Callable[[lx.LexiconRow], str]
) -> dict[str, Counts]:
    """mutated / non-mutated per key; NA mutation counted in neither."""
    out: dict[str, list[int]] = {}
    for row in rows:
        cell = out.setdefault(key(row), [0, 0])
        if row.mutation is True:
            cell[0] += 1
        elif row.mutation is False:
            cell[1] += 1
    return {k: (v[0], v[1]) for k, v in out.items()}


def _segment_rows(
    tp_rows: Sequence[lx.LexiconRow], suffix: str
) -> dict[str, Counts]:
    simple = [r for r in tp_rows if not r.cluster]
    out: dict[str, Counts] = {}
    by_opp = _tally(simple, lambda r: f"{r.stem_final}\t{r.opportunity}")
    for key, counts in by_opp.items():
        seg, opp = key.split("\t")
        out[f"<{seg}> + <-{opp}> plural{suffix}"] = counts
    for seg, counts in _tally(simple, lambda r: r.stem_final).items():
        out[f"<{seg}> + <-i, -e> plural{suffix}"] = counts
    return out


def detect_palatal_from_orth(lemma: str, verb: str, seg: str) -> bool | None:
    """The R script's ``detect_palatal_from_orth`` (whole list cell)."""
    text = verb.lower()
    if not text:
        return None
    tail = text[-5:]
    if seg == "c":
        return bool(re.search(r"(?<!h)c[ie]", tail))
    if seg == "g":
        return bool(re.search(r"(?<!h)g[ie]", tail))
    if seg == "t":
        return "ți" in tail or "țe" in tail
    if seg == "d":
        return text.endswith("zi") and not lemma.lower().endswith("z")
    if seg == "s":
        return "și" in tail
    if seg == "z":
        return "ji" in tail
    return None


def _distinct_lemma(
    rows: Iterable[tuple[lx.LexiconRow, bool]],
) -> list[tuple[lx.LexiconRow, bool]]:
    ordered = sorted(
        rows, key=lambda item: item[0].lemma
    )  # stable, as arrange()
    seen: set[str] = set()
    out: list[tuple[lx.LexiconRow, bool]] = []
    for row, flag in ordered:
        if row.lemma in seen:
            continue
        seen.add(row.lemma)
        out.append((row, flag))
    return out


def _base_inflect(nouns: Sequence[lx.LexiconRow]) -> list[lx.LexiconRow]:
    return [
        r for r in nouns if r.stem_final in SEGMENTS and r.mutation is not None
    ]


def _plural_class(row: lx.LexiconRow) -> str:
    if row.opportunity in ("i", "e"):
        return "front"
    if row.opportunity == "uri":
        return "back"
    return "ambig"


def _deriv_pairs(
    base: Sequence[lx.LexiconRow], column: str, *, need_ipa: bool = False
) -> list[tuple[lx.LexiconRow, bool]]:
    pairs: list[tuple[lx.LexiconRow, bool]] = []
    for row in base:
        values = row.derived_verbs if column == "verb" else row.derived_adj
        if not values or row.exception_reason.startswith("nde:"):
            continue
        if need_ipa and not row.ipa_derived_adj:
            continue
        flag = detect_palatal_from_orth(
            row.lemma, "|".join(values), row.stem_final
        )
        if flag is not None:
            pairs.append((row, flag))
    return _distinct_lemma(pairs)


def _deriv_rows(
    pairs: Sequence[tuple[lx.LexiconRow, bool]],
    summary: Sequence[tuple[lx.LexiconRow, bool]],
    label: str,
    suffix: str,
) -> dict[str, Counts]:
    out: dict[str, Counts] = {}
    for mutated, name in (
        (False, "base plural non-mut."),
        (True, "base plural mutated"),
    ):
        chosen = [flag for row, flag in summary if row.mutation is mutated]
        if chosen:
            out[f"{label} derivation, {name}{suffix}"] = (
                sum(chosen),
                len(chosen) - sum(chosen),
            )
    flags = [flag for _, flag in pairs]
    if flags:
        out[f"{label} derivation, ALL{suffix}"] = (
            sum(flags),
            len(flags) - sum(flags),
        )
    return out


def _nde_rows(opp_rows: Sequence[lx.LexiconRow]) -> dict[str, Counts]:
    tallied = _tally(
        [r for r in opp_rows if legacy_nde(r) in NDEB], legacy_nde
    )
    return {NDEB_LABEL[k]: v for k, v in tallied.items()}


def _with_ndeb(
    opp_rows: Sequence[lx.LexiconRow], label: str
) -> dict[str, Counts]:
    tallied = _tally(
        [r for r in opp_rows if not r.cluster], lambda r: r.stem_final
    )
    return {
        f"<{seg}> + <-i, -e> plural, {label}WITH NDEB": tallied[seg]
        for seg in ("c", "g")
        if seg in tallied
    }


def _cluster_rows(
    tp_rows: Sequence[lx.LexiconRow], suffix: str
) -> dict[str, Counts]:
    chosen = [r for r in tp_rows if cluster_simple(r.cluster)]
    out: dict[str, Counts] = {}
    for key, counts in _tally(
        chosen, lambda r: f"{cluster_simple(r.cluster)}\t{r.opportunity}"
    ).items():
        name, opp = key.split("\t")
        out[f"<{name}> + <-{opp}> plural{suffix}"] = counts
    for name, counts in _tally(
        chosen, lambda r: cluster_simple(r.cluster)
    ).items():
        out[f"<{name}> + <-i, -e> plural{suffix}"] = counts
    return out


def downsample(
    opp_rows: Sequence[lx.LexiconRow], n: int = 1000
) -> list[lx.LexiconRow]:
    """``build_downsampled_lexicon``: core :func:`lexpipe.tp.top_n`, then
    one row per lemma after arranging by (-freq, lemma, plural)."""
    items: list[dict[str, Any]] = [
        {"lemma": r.lemma, "freq": r.freq, "index": i}
        for i, r in enumerate(opp_rows)
    ]
    chosen = [opp_rows[item["index"]] for item in tp.top_n(items, n, "freq")]
    ordered = sorted(chosen, key=lambda r: (-r.freq, r.lemma, r.plural))
    seen: set[str] = set()
    out: list[lx.LexiconRow] = []
    for row in ordered:
        if row.lemma not in seen:
            seen.add(row.lemma)
            out.append(row)
    return out


def legacy_tables(
    rows: Sequence[lx.LexiconRow],
) -> tuple[dict[str, Counts], dict[str, Counts]]:
    """(full, downsampled) legacy rows as label → (mutated, non-mutated)."""
    nouns = [r for r in rows if r.pos == "N"]
    nouns_tp = [r for r in nouns if r.tp_in_domain]
    opp = [
        r
        for r in nouns
        if opportunity_tp(r) in ("i", "e") and r.stem_final in SEGMENTS
    ]
    base = _base_inflect(nouns)
    full: dict[str, Counts] = {}
    full.update(_nde_rows(opp))
    full.update(_segment_rows(nouns_tp, ""))
    full.update(_with_ndeb(opp, ""))
    full.update(_cluster_rows(nouns_tp, ""))
    verbs = _deriv_pairs(base, "verb")
    full.update(
        _deriv_rows(
            verbs,
            [p for p in verbs if _plural_class(p[0]) == "front"],
            "N->V",
            "",
        )
    )
    adjectives = _deriv_pairs(base, "adj")
    full.update(_deriv_rows(adjectives, adjectives, "N->Adj", ""))

    ds = downsample(opp)
    lemmas = {r.lemma for r in ds}
    ds_tp = [r for r in ds if r.tp_in_domain]
    small: dict[str, Counts] = {}
    small.update(_nde_rows(ds))
    small.update(_segment_rows(ds_tp, " (downsampled)"))
    small.update(_with_ndeb(ds, "downsampled, "))
    small.update(_cluster_rows(ds_tp, " (downsampled)"))
    base_ds = [r for r in base if r.lemma in lemmas]
    verbs_ds = _deriv_pairs(base_ds, "verb")
    small.update(_deriv_rows(verbs_ds, verbs_ds, "N->V", " (downsampled)"))
    adjectives_ds = _deriv_pairs(base_ds, "adj", need_ipa=True)
    small.update(
        _deriv_rows(adjectives_ds, adjectives_ds, "N->Adj", " (downsampled)")
    )
    return full, small


def _parse_number(text: str) -> float | None:
    return None if text in ("", "NA") else float(text)


def _parse_flag(text: str) -> bool | None:
    return None if text in ("", "NA") else text.upper() == "TRUE"


@dataclass(frozen=True, slots=True)
class ParityLine:
    """One legacy CSV row against its recomputation."""

    label: str
    legacy: tuple[Any, ...]
    computed: tuple[Any, ...]

    @property
    def ok(self) -> bool:
        a, b = self.legacy, self.computed
        if a[:3] != b[:3] or a[4:] != b[4:]:
            return False
        if a[3] is None or b[3] is None:
            return a[3] is b[3]
        return math.isclose(a[3], b[3], rel_tol=1e-12, abs_tol=1e-15)


def compare_csv(
    path: Path, computed: Mapping[str, Counts]
) -> list[ParityLine]:
    """Every numeric row of a legacy summary CSV vs :func:`tp.tp_row`."""
    out: list[ParityLine] = []
    with path.open("r", encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            n = _parse_number(row["N"])
            if n is None:
                continue
            label = row["type"]
            legacy = (
                int(n),
                int(float(row["mutated"])),
                int(float(row["non-mutated"])),
                _parse_number(row["rate"]),
                _parse_flag(row["majority?"]),
                _parse_flag(row["tolerated?"]),
            )
            if label in computed:
                mutated, non = computed[label]
                again = tp.tp_row(label, mutated + non, mutated)
                values: tuple[Any, ...] = (
                    again["N"],
                    again["mutated"],
                    again["non-mutated"],
                    again["rate"],
                    again["majority?"],
                    again["tolerated?"],
                )
            else:
                values = ("missing",) * 6
            out.append(ParityLine(label, legacy, values))
    return out


# ---------------------------------------------------------------------------
# 2. Side by side with reasons
# ---------------------------------------------------------------------------

ADD_REASONS: Final[tuple[str, ...]] = (
    "cluster (legacy reports st/sc/ct in cluster rows)",
    "legacy NDE class (assigned to non-palatalizing items only)",
    "legacy already-palatal check (non-palatalizing items only)",
    "legacy -ică suppletion / lost digraph",
    "legacy non exception (other)",
    "suffix-internal target (legacy target_is_suffix)",
    "legacy suppletion rule (non-palatalizing items only)",
    "legacy foreign-plural rule (non-palatalizing items only)",
    "other",
)


def _legacy_shared_short(record: RootRecord) -> bool:
    lemma, plural = record.lemma.lower(), record.plural.lower()
    shared = 0
    for a, b in zip(lemma, plural, strict=False):
        if a != b:
            break
        shared += 1
    return shared < max(2, len(lemma) * 2 // 5)


def add_reason(record: RootRecord) -> str:
    """Why a harmonized root is outside the legacy domain (legacy
    order: NDE, non exception, suffix target, suppletive, foreign)."""
    legacy = record.legacy
    if record.cluster:
        return ADD_REASONS[0]
    if legacy.exception_reason.startswith("nde:"):
        return ADD_REASONS[1]
    if legacy.exception_reason == "non exception":
        if "palatal_context" in record.flags:
            return ADD_REASONS[2]
        if "ica_suppletion" in record.flags or record.cluster in (
            "che",
            "ghe",
            "chi",
            "ghi",
        ):
            return ADD_REASONS[3]
        return ADD_REASONS[4]
    if legacy.target_is_suffix:
        return ADD_REASONS[5]
    if (
        legacy.mutation is False
        and record.plural
        and _legacy_shared_short(record)
    ):
        return ADD_REASONS[6]
    seg = record.seg
    if legacy.mutation is False and legacy.orth_change in (
        f"{seg}→{seg}i",
        f"{seg}o→{seg}i",
    ):
        return ADD_REASONS[7]
    return ADD_REASONS[8]


def remove_reason(record: RootRecord) -> str:
    if record.nde:
        return f"harmonized structural NDE ({record.nde})"
    if not record.is_target:
        return "not target-final in the harmonized classes"
    return "other"


@dataclass(frozen=True, slots=True)
class Reconciliation:
    """Legacy vs harmonized membership of one cell."""

    seg: str
    trigger: str
    legacy_n: int
    legacy_e: int
    harm_n: int
    harm_e: int
    added: Mapping[str, tuple[int, int]]
    removed: Mapping[str, tuple[int, int]]

    @property
    def balanced(self) -> bool:
        plus = sum(n for n, _ in self.added.values())
        minus = sum(n for n, _ in self.removed.values())
        plus_e = sum(e for _, e in self.added.values())
        minus_e = sum(e for _, e in self.removed.values())
        return (
            self.legacy_n + plus - minus == self.harm_n
            and self.legacy_e + plus_e - minus_e == self.harm_e
        )


def reconcile(
    records: Sequence[RootRecord],
    cache: ObservationCache,
    seg: str,
    opportunity: str,
    *,
    restrict: frozenset[str] | None = None,
    subset: str = "all",
) -> Reconciliation:
    """Root-by-root comparison of a legacy segment row and the
    harmonized cell (E0, dictionary, dominant coding)."""
    observed = cache.get("E0")
    trigger = f"-{opportunity}"
    member = SUBSETS[subset]

    def keep(record: RootRecord) -> bool:
        return restrict is None or record.root_id in restrict

    legacy: dict[str, bool] = {}
    for r in records:
        if not keep(r):
            continue
        lg = r.legacy
        if (
            lg.tp_in_domain
            and not r.cluster
            and r.seg == seg
            and lg.opportunity == opportunity
        ):
            if lg.mutation is not None:
                legacy[r.root_id] = lg.mutation is False
    harmonized: dict[str, bool] = {}
    for r in records:
        if not keep(r) or r.seg != seg or not r.is_target or not member(r):
            continue
        status = tp.root_status(
            observed[r.root_id], trigger, "dictionary", "dominant"
        )
        if status in (
            tp.RootStatus.ALTERNATES,
            tp.RootStatus.NONALT,
            tp.RootStatus.VARIABLE,
        ):
            harmonized[r.root_id] = status is not tp.RootStatus.ALTERNATES
    by_id = {r.root_id: r for r in records}
    added: Counter[str] = Counter()
    added_e: Counter[str] = Counter()
    for root_id in sorted(set(harmonized) - set(legacy)):
        reason = add_reason(by_id[root_id])
        added[reason] += 1
        added_e[reason] += int(harmonized[root_id])
    removed: Counter[str] = Counter()
    removed_e: Counter[str] = Counter()
    for root_id in sorted(set(legacy) - set(harmonized)):
        reason = remove_reason(by_id[root_id])
        removed[reason] += 1
        removed_e[reason] += int(legacy[root_id])
    return Reconciliation(
        seg=seg,
        trigger=trigger,
        legacy_n=len(legacy),
        legacy_e=sum(legacy.values()),
        harm_n=len(harmonized),
        harm_e=sum(harmonized.values()),
        added={k: (added[k], added_e[k]) for k in sorted(added)},
        removed={k: (removed[k], removed_e[k]) for k in sorted(removed)},
    )


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


def _fmt(values: tuple[Any, ...]) -> str:
    n, mutated, non, rate, majority, tolerated = values
    rate_text = "NA" if rate is None else f"{rate:.4f}"
    return (
        f"N={n} mut={mutated} non={non} rate={rate_text} "
        f"maj={majority} tol={tolerated}"
    )


def _harm(
    rows: Sequence[Mapping[str, Any]], **match: Any
) -> Mapping[str, Any] | None:
    found = [r for r in rows if all(r.get(k) == v for k, v in match.items())]
    return found[0] if len(found) == 1 else None


def render_report(
    *,
    lexicon_rows: Sequence[lx.LexiconRow],
    records: Sequence[RootRecord],
    cache: ObservationCache,
    tp_long: Sequence[Mapping[str, Any]],
    full_csv: Path,
    downsampled_csv: Path,
    inputs: Mapping[str, str],
) -> tuple[str, bool]:
    """The parity report and whether every legacy row was reproduced."""
    full, small = legacy_tables(lexicon_rows)
    checks = [
        ("analysis/romanian_tp_summary_full.csv", compare_csv(full_csv, full)),
        (
            "analysis/romanian_tp_summary_downsampled.csv",
            compare_csv(downsampled_csv, small),
        ),
    ]
    all_ok = all(line.ok for _, lines in checks for line in lines)
    lines: list[str] = []
    add = lines.append
    add(
        "Romanian TP parity: legacy R tables vs the shared core (design §7, §8)"
    )
    add("=" * 76)
    add("")
    add("Inputs (SHA-256):")
    for name, digest in sorted(inputs.items()):
        add(f"  {name}  {digest}")
    add("")
    add("1. Legacy rows recomputed with lexpipe.tp.tp_row")
    add("-" * 76)
    for name, results in checks:
        bad = [line for line in results if not line.ok]
        same = len(results) - len(bad)
        add(f"  {name}: {same}/{len(results)} numeric rows identical")
        for line in bad:
            add(f"    DIFF {line.label}: legacy {_fmt(line.legacy)}")
            add(f"         computed {line.computed}")
    add("")
    add(
        "  Definitions reproduced: tp_in_domain rows without a cluster grouped"
    )
    add("  by stem_final x opportunity; 'WITH NDEB' and NDE rows from the")
    add("  opportunity_tp domain; cluster rows by st/sc/ct; top-1000 by the")
    add("  lemma's maximum frequency (core top_n), one row per lemma;")
    add("  derivational rows with the R detector on the whole derived list.")
    add("")
    add("2. Harmonized cells next to the legacy rows (E0, tau = dictionary)")
    add("-" * 76)
    add("  harmonized: N = N_trigger (roots attested before the trigger),")
    add("  e = palatalize-rule exceptions. 'all' keeps every outcome-blind")
    add("  flag; 'core' removes them (suffix target, palatal context,")
    add("  suppletive stem, -ică suppletion). Structural NDE classes are")
    add("  separate rows in both.")
    add("")
    add(
        f"  {'legacy row':28s} {'legacy N/e':>12s} {'all N/e':>12s} {'core N/e':>12s}"
    )
    for seg in SEGMENTS:
        for opp, trigger in (("i", "-i"), ("e", "-e"), ("i, -e", "-i/-e")):
            label = f"<{seg}> + <-{opp}> plural"
            legacy_counts = full.get(label)
            legacy_text = (
                "-"
                if legacy_counts is None
                else f"{sum(legacy_counts)}/{legacy_counts[1]}"
            )
            cells = []
            for subset in ("all", "core"):
                row = _harm(
                    tp_long,
                    table="cells",
                    target_class=seg,
                    trigger=trigger,
                    subset=subset,
                    tau="dictionary",
                    evidence_set="E0",
                    rule="palatalize",
                )
                cells.append(
                    "-" if row is None else f"{row['N_trigger']}/{row['e']}"
                )
            add(
                f"  {label:28s} {legacy_text:>12s} {cells[0]:>12s} {cells[1]:>12s}"
            )
    add("")
    add("3. Root-by-root reconciliation (legacy row → harmonized 'all' cell)")
    add("-" * 76)
    balanced = True
    for seg in SEGMENTS:
        for opp in ("i", "e"):
            rec = reconcile(records, cache, seg, opp)
            balanced &= rec.balanced
            add(
                f"  <{seg}> + <-{opp}>: legacy N={rec.legacy_n} e={rec.legacy_e} → "
                f"harmonized N={rec.harm_n} e={rec.harm_e} "
                f"({'balances' if rec.balanced else 'DOES NOT BALANCE'})"
            )
            for reason, (n, e) in rec.added.items():
                add(f"      + {n:5d} roots ({e:5d} non-alternating)  {reason}")
            for reason, (n, e) in rec.removed.items():
                add(f"      - {n:5d} roots ({e:5d} non-alternating)  {reason}")
    add("")
    add("  Reasons marked 'non-palatalizing items only' are the outcome-")
    add("  conditioned exclusions of the legacy code (critique P0 #1): the")
    add("  harmonized cells keep those roots, as exceptions.")
    add("")
    add(
        "4. Top-1000: legacy frame vs harmonized frames (E0, dictionary, 'all')"
    )
    add("-" * 76)
    opp_rows = [
        r
        for r in lexicon_rows
        if r.pos == "N"
        and opportunity_tp(r) in ("i", "e")
        and r.stem_final in SEGMENTS
    ]
    legacy_frame = frozenset(
        lx.ROOT_PREFIX + r.lemma for r in downsample(opp_rows)
    )
    r24 = frame(records, 1000, "lemma")
    add(
        "  legacy frame: the 1,000 most frequent nouns with an -i/-e plural and"
    )
    add("  a target final (R build_downsampled_lexicon); R24 frame: the 1,000")
    add("  most frequent nouns of any final segment and plural type.")
    add(
        f"  {'cell':12s} {'legacy row':>12s} "
        f"{'harm. in legacy frame':>22s} {'harm. in R24 frame':>20s}"
    )
    for seg in SEGMENTS:
        for opp in ("i", "e"):
            legacy_counts = small.get(
                f"<{seg}> + <-{opp}> plural (downsampled)"
            )
            legacy_text = (
                "-"
                if legacy_counts is None
                else f"{sum(legacy_counts)}/{legacy_counts[1]}"
            )
            in_legacy = reconcile(
                records, cache, seg, opp, restrict=legacy_frame
            )
            in_r24 = reconcile(records, cache, seg, opp, restrict=r24)
            add(
                f"  <{seg}> + <-{opp}>  {legacy_text:>12s} "
                f"{in_legacy.harm_n:>16d}/{in_legacy.harm_e:<5d} "
                f"{in_r24.harm_n:>14d}/{in_r24.harm_e:<5d}"
            )
    add("")
    add("5. NDE rows: legacy classes vs harmonized structural classes")
    add("-" * 76)
    for label in ("gimpe type", "ochi-ochi type", "paduche-paduchi type"):
        counts = full.get(label, (0, 0))
        add(f"  legacy {label:22s} N={sum(counts)} (mutated {counts[0]})")
    structural = Counter(r.nde for r in records if r.nde)
    for name in lx.NDE_CLASSES:
        add(f"  harmonized nde:{name:18s} roots={structural[name]}")
    add("  Legacy NDE is assigned only to non-mutating rows and mixes coronal")
    add("  patterns (-iște, invariant -i loans) with the dorsal ochi/păduche")
    add("  types; the structural classes are decided by the lemma alone.")
    add("")
    add("6. Verb evidence (E1): audit verdicts vs this module's cell reading")
    add("-" * 76)
    pairs = Counter(
        (r.audit_verdict, mo.audit_like_verdict(r.seg, list(_readings(r))))
        for r in records
        if r.verb
    )
    for (audit, mine), n in sorted(pairs.items()):
        add(f"  audit {audit:15s} this module {mine:15s} {n:4d}")
    add("  Differences: (a) the audit's detector takes the z of the -ez-")
    add("  augment for the palatal reflex of d (arendezi, legendezi:")
    add("  ALTERNATES); here d + -ez- is a coronal before /e/, not a trigger;")
    add("  (b) the audit counts any cell with the reflex (amănunțești,")
    add("  îmbogățești: ț before the -esc- augment), here only cells whose")
    add(
        "  vowel triggers the class's rule count, so those are evidence of the"
    )
    add(
        "  -e cell; (c) the audit reads every augmented cell as TRIGGER-ABSENT,"
    )
    add("  here a dorsal before -ez- (blochezi) meets /e/ unpalatalized:")
    add(
        "  NO-ALTERNATION; (d) NO-READING: the verb's stem does not end in the"
    )
    add("  noun's target (picur-a, strunj-i, micșor-a), so no cell is read.")
    add("")
    add("7. Doublets (R19)")
    add("-" * 76)
    multi = sum(1 for r in records if r.duplicates > 1)
    conflicts = plural_conflicts(lexicon_rows)
    add(f"  noun lemmas with more than one lexicon row: {multi}")
    add(
        "  noun lemmas whose rows list different plurals: "
        f"{len(conflicts)}{': ' + ', '.join(conflicts) if conflicts else ''}"
    )
    add("  Each root keeps the plural of its priority row (a row with a")
    add("  plural first), so with no conflicts every root has at most one")
    add("  plural and the three doublet codings give identical E0 rows.")
    add("  Doublets arise only from E1/E2 evidence (verb cells, derivatives).")
    return "\n".join(lines) + "\n", all_ok and balanced


def plural_conflicts(rows: Sequence[lx.LexiconRow]) -> list[str]:
    """Noun lemmas whose rows list more than one distinct plural."""
    plurals: dict[str, set[str]] = {}
    for row in rows:
        if row.pos == "N" and row.plural:
            plurals.setdefault(row.lemma, set()).add(row.plural)
    return sorted(lemma for lemma, found in plurals.items() if len(found) > 1)


def _readings(record: RootRecord) -> Iterable[mo.Reading]:
    for cell in record.verb_cells:
        yield mo.Reading(cell.form, cell.vowel, cell.palatal, cell.augment)
