"""Headline tables of the Czech TP stage (``czech_tp_summary.txt``).

Reads the long rows the ``tp`` stage computed (:mod:`.tpstage`) and
prints the cells a reader needs first: (a) the rule TP and the R16
non-trigger rows, (b) the allomorph competitions against the baseline,
(c) the leave-one-cell-out spell-out TP over E0–E2, the R17 rule
family, and the simulated learners at V = 1000. Every figure is a
count of root types or a share; no per-word counts (R9). Lexeme
summaries (undergoer / resister / mixed) are recomputed from the root
records at τ = dictionary and τ ≥ 1.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from typing import Any, Final

from lexpipe.langs.czech import ortho
from lexpipe.langs.czech import suffixes as sx
from lexpipe.langs.czech import tpcells as tc
from lexpipe.langs.czech.roots import NOUN, RootRecord

RULE_FOCUS: Final[tuple[str, ...]] = ("velar", "r", "c", "coronal_inferred")
COMP_FOCUS: Final[tuple[str, ...]] = (
    "velar",
    "r",
    "coronal_inferred",
    "baseline",
)
SHOW_TAUS: Final[tuple[str, ...]] = (
    "dictionary",
    "ge1",
    "ge1_strict",
    "ge2pm",
)
RULE_WIDTH: Final[int] = 30

Key = tuple[str, str, str, str, str, str, str, str, int, str]


def _key(row: Mapping[str, Any]) -> Key:
    v = row.get("V")
    return (
        str(row["table"]),
        str(row["target_class"]),
        str(row["trigger"]),
        str(row["subset"]),
        str(row["tau"]),
        str(row["evidence_set"]),
        str(row["coding"]),
        str(row["freq_basis"]),
        -1 if v is None else int(v),
        str(row["rule"]),
    )


class Rows:
    """Long rows indexed by their labels."""

    def __init__(self, rows: Iterable[Mapping[str, Any]]) -> None:
        self._index = {_key(row): row for row in rows}

    def get(
        self,
        table: str,
        target: str,
        trigger: str,
        *,
        subset: str = "all",
        tau: str = "dictionary",
        evidence: str = "E0",
        coding: str = "dominant",
        basis: str = "lemma",
        V: int = -1,
        rule: str = "palatalize",
    ) -> Mapping[str, Any] | None:
        return self._index.get(
            (
                table,
                target,
                trigger,
                subset,
                tau,
                evidence,
                coding,
                basis,
                V,
                rule,
            )
        )


def _theta(value: Any) -> str:
    if value is None or value == "":
        return "-"
    number = float(value)
    return f"{number:.1f}" if math.isfinite(number) else "-"


def _tol(value: Any) -> str:
    if value is None:
        return "-"
    return "yes" if value is True else "no"


def _cells(row: Mapping[str, Any] | None) -> str:
    """``N  e  θ  tol  frag`` of one row ("-" when absent)."""
    if row is None:
        return f"{'-':>6} {'-':>5} {'-':>7} {'-':>4} {'-':>5}"
    frag = row.get("fragility")
    return (
        f"{row['N_trigger']:>6} {row['e']:>5} "
        f"{_theta(row['theta_trigger']):>7} "
        f"{_tol(row['tolerated_trigger']):>4} "
        f"{'-' if frag is None else frag:>5}"
    )


HEAD: Final[str] = f"{'N':>6} {'e':>5} {'theta':>7} {'tol':>4} {'frag':>5}"


def _section(title: str) -> list[str]:
    return ["", title, "-" * 72]


# ---------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------


def rule_section(rows: Rows) -> list[str]:
    lines = _section(
        "(a) Rule TP per trigger cell (palatalize; e = P form without "
        "the change)"
    )
    lines.append(
        f"  {'target':<18} {'trigger':<{RULE_WIDTH}} {'tau':<10} {HEAD}"
    )
    for trigger in tc.RULE_TRIGGERS:
        table = "rule" if trigger.pos == NOUN else "rule_adj"
        for target in RULE_FOCUS:
            for tau in SHOW_TAUS:
                row = rows.get(table, target, trigger.name, tau=tau)
                if row is None:
                    continue
                lines.append(
                    f"  {target:<18} {trigger.name:<{RULE_WIDTH}} "
                    f"{tau:<10} {_cells(row)}"
                )
    lines.append(
        "  coronal_inferred: the change of t/d/n is inferred from ě/i/í "
        "(spelling), so its e is 0 by construction."
    )
    return lines


def competition_section(rows: Rows) -> list[str]:
    lines = _section(
        "(b) Competition TP: palatalize = P-default, retain = N-default "
        "(dominant coding)"
    )
    lines.append(
        f"  {'cell':<18} {'target':<17} {'tau':<10} "
        f"{'P-default ' + HEAD:<44} {'N-default e/tol':>16}"
    )
    for comp in tc.COMPETITIONS:
        table = "competition" if comp.pos == NOUN else "competition_adj"
        for target in COMP_FOCUS:
            for tau in SHOW_TAUS:
                pal = rows.get(table, target, comp.name, tau=tau)
                ret = rows.get(
                    table, target, comp.name, tau=tau, rule="retain"
                )
                if pal is None:
                    continue
                retain = (
                    "-"
                    if ret is None
                    else f"{ret['e']}/{_tol(ret['tolerated_trigger'])}"
                )
                lines.append(
                    f"  {comp.name:<18} {target:<17} {tau:<10} "
                    f"{_cells(pal):<44} {retain:>16}"
                )
    lines += [
        "",
        "  ge1_strict: >= 1 token, counting only forms that are neither",
        "  ambiguous with another lexeme nor shared with a cell outside",
        "  their copy group (LOC.SG -u = GEN/DAT.SG; VOC -u = DAT/LOC.SG):",
        "  the ge1 counts of those cells are upper bounds.",
        "",
        "  Doublet codings (R19), velars, tau >= 1: e of the P-default rule",
    ]
    for comp in tc.COMPETITIONS:
        table = "competition" if comp.pos == NOUN else "competition_adj"
        parts = []
        for coding in ("dominant", "doublet_exception", "doublet_excluded"):
            row = rows.get(table, "velar", comp.name, tau="ge1", coding=coding)
            if row is not None:
                parts.append(f"{coding} N={row['N_trigger']} e={row['e']}")
        if parts:
            lines.append(f"    {comp.name:<18} " + "; ".join(parts))
    return lines


def nested_section(rows: Rows) -> list[str]:
    lines = _section(
        "(b) Nested classes, velars, tau >= 1 (Yang order: suffixoid "
        "-log/-urg/-graf, -sky/-cky, other suffix, monosyllable, loan, "
        "elsewhere)"
    )
    for comp in tc.COMPETITIONS:
        table = "competition" if comp.pos == NOUN else "competition_adj"
        parts = []
        for subset in (
            "nested:suffixoid",
            "nested:sky",
            "nested:suffix",
            "nested:monosyllable",
            "nested:loan",
            "nested:elsewhere",
        ):
            for rule in ("palatalize", "retain"):
                row = rows.get(
                    table,
                    "velar",
                    comp.name,
                    subset=subset,
                    tau="ge1",
                    rule=rule,
                )
                if row is not None:
                    parts.append(
                        f"{subset[7:]}:{rule[:3]} N={row['N_trigger']} "
                        f"e={row['e']} {_tol(row['tolerated_trigger'])}"
                    )
        if parts:
            lines.append(f"  {comp.name:<18} " + "; ".join(parts))
    return lines


def spellout_section(rows: Rows) -> list[str]:
    lines = _section(
        "(c) Spell-out TP, leave one cell out (velars; dominant coding)"
    )
    lines.append(
        "  underspecified = alternation evidence outside the tested cell; "
        "e = roots taking the N allomorph"
    )
    lines.append(
        "  prespecified = no such evidence; e = roots taking the P allomorph"
    )
    lines.append(
        f"  {'tested cell':<14} {'E':<3} {'tau':<10} "
        f"{'underspecified ' + HEAD:<49} {'prespecified ' + HEAD}"
    )
    for name in tc.SPELLOUT_CELLS:
        for evidence in tc.EVIDENCE_SETS:
            for tau in SHOW_TAUS:
                under = rows.get(
                    "spellout_loo",
                    "velar",
                    name,
                    subset="underspecified",
                    tau=tau,
                    evidence=evidence,
                )
                pre = rows.get(
                    "spellout_loo",
                    "velar",
                    name,
                    subset="prespecified",
                    tau=tau,
                    evidence=evidence,
                    rule="retain",
                )
                if under is None and pre is None:
                    continue
                lines.append(
                    f"  {name:<14} {evidence:<3} {tau:<10} "
                    f"{_cells(under):<49} {_cells(pre)}"
                )
    return lines


def r17_section(rows: Rows) -> list[str]:
    lines = _section(
        "R17 fixed rule family (nouns, all roots, dominant coding, "
        "palatalize row; non-trigger cells: retain row)"
    )
    lines.append(
        f"  {'table':<6} {'target':<17} {'trigger':<14} {'E':<3} "
        f"{'tau':<10} {HEAD}"
    )
    for definition in tc.R17_DEFS:
        rule = "retain" if "non-trigger" in definition.memo else "palatalize"
        for evidence in ("E0", "E2"):
            for tau in ("dictionary", "ge1"):
                row = rows.get(
                    definition.table,
                    definition.target,
                    definition.trigger,
                    tau=tau,
                    evidence=evidence,
                    rule=rule,
                )
                if row is None:
                    continue
                lines.append(
                    f"  {definition.table:<6} {definition.target:<17} "
                    f"{definition.trigger:<14} {evidence:<3} {tau:<10} "
                    f"{_cells(row)}{'  (retain)' if rule == 'retain' else ''}"
                )
    return lines


def derivation_section(rows: Rows) -> list[str]:
    lines = _section(
        "R22 derivational trigger rows (velar bases, basis wikt; "
        "palatalize row)"
    )
    lines.append(f"  {'suffix':<12} {'tau':<10} {HEAD}  memo")
    for suffix in sx.SUFFIXES:
        for tau in ("dictionary", "ge1"):
            row = rows.get(
                "derivation", "velar", suffix.label, tau=tau, evidence="E2"
            )
            if row is None:
                continue
            lines.append(
                f"  {suffix.label:<12} {tau:<10} {_cells(row)}  {row['memo']}"
            )
    return lines


def _pal_ret(row: Mapping[str, Any] | None) -> str:
    if row is None:
        return "-"
    pal = float(row["share_palatalize"])
    ret = float(row["share_retain"])
    return f"{pal:.3f}/{ret:.3f}"


def learner_section(learner_rows: Sequence[Mapping[str, Any]]) -> list[str]:
    """V = 1000 shares on the lemma basis, and the V = 5000 share of the
    palatalize (P-default) rule on each R24 basis (the sensitivity the
    audit asked for: homograph lemmas inflate the lemma basis)."""
    lines = _section(
        "Simulated learners (R18), E0, 1000 vocabularies drawn by token "
        "frequency"
    )
    lines += [
        "  V = 1000 on the lemma basis: shares of learners with palatalize,",
        "  retain, neither or a conditioned rule set productive; then at",
        "  V = 5000 the palatalize/retain shares on each R24 basis: lemma",
        "  (lemma-form tokens), paradigm (attested forms summed), unamb.",
        "  (the same without forms of other lexemes: jak 'yak' loses the",
        "  conjunction's tokens)",
        f"  {'table':<13} {'target':<17} {'trigger':<14} {'subset':<15} "
        f"{'tau':<11} {'pal':>6} {'ret':>6} {'nei':>6} {'cond':>6}  "
        f"{'5000 lemma':>11} {'paradigm':>11} {'unamb.':>11}",
    ]
    index: dict[tuple[str, ...], Mapping[str, Any]] = {}
    order: list[tuple[str, ...]] = []
    for row in learner_rows:
        key = (
            str(row["table"]),
            str(row["target_class"]),
            str(row["trigger"]),
            str(row["subset"]),
            str(row["tau"]),
        )
        if key not in index and key not in order:
            order.append(key)
        index[(*key, str(row["freq_basis"]), str(row["V"]))] = row
    for group in order:
        base = index.get((*group, "lemma", "1000"))
        if base is None:
            continue
        five = [
            index.get((*group, basis, "5000"))
            for basis in ("lemma", "paradigm", "paradigm_unambiguous")
        ]
        shares = " ".join(f"{_pal_ret(r):>11}" for r in five)
        lines.append(
            f"  {group[0]:<13} {group[1]:<17} {group[2]:<14} "
            f"{group[3]:<15} "
            f"{group[4]:<11} {float(base['share_palatalize']):>6.3f} "
            f"{float(base['share_retain']):>6.3f} "
            f"{float(base['share_neither']):>6.3f} "
            f"{float(base['share_conditioned']):>6.3f}  {shares}"
        )
    return lines


def lexeme_summaries(records: Sequence[RootRecord]) -> list[str]:
    lines = _section(
        "Root summaries (E0): undergoer / resister / mixed / no evidence"
    )
    for target in ("velar", "r"):
        for tau in ("dictionary", "ge1"):
            counts: Counter[str] = Counter()
            for record in records:
                if record.pos != NOUN or record.target_class != target:
                    continue
                counts[_summary(record, tau)] += 1
            total = sum(counts.values())
            parts = "  ".join(
                f"{name} {counts[name]:,}"
                for name in ("undergoer", "resister", "mixed", "no_evidence")
            )
            lines.append(f"  {target:<6} {tau:<10} roots {total:,}: {parts}")
    lines.append(
        "  (undergoer: an attested P form with the change; resister: only "
        "competing N allomorphs; mixed: both)"
    )
    return lines


def _summary(record: RootRecord, tau: str) -> str:
    p = n = False
    for item in record.realizations:
        if not item.attested_at(tau):
            continue
        if item.allomorph == ortho.P and item.mutation is True:
            p = True
        elif item.allomorph == ortho.N and item.n_kind == ortho.N_ALLOMORPH:
            n = True
    if p and n:
        return "mixed"
    if p:
        return "undergoer"
    if n:
        return "resister"
    return "no_evidence"


def render(
    rows: Sequence[Mapping[str, Any]],
    learner_rows: Sequence[Mapping[str, Any]],
    records: Sequence[RootRecord],
) -> str:
    """The whole summary text."""
    index = Rows(rows)
    counts = Counter(f"{r.pos}:{r.target_class}" for r in records)
    lines = [
        "Czech Tolerance Principle analyses (design §6, R15-R19, R24)",
        "=" * 72,
        "",
        "Root types (R4) by part of speech and target class; aggregates "
        "only (R9).",
        "N = N_trigger (roots attested in the cell at tau), e = exceptions "
        "of the rule,",
        "theta = N/ln N, tol = e <= theta, frag = fewest item changes that "
        "flip tol.",
        "tau: dictionary (every Wiktionary form), ge1 (>= 1 token), "
        "ge2pm (>= 2 per million)",
        "in Leipzig ces_news_2024_1M. Full grid: czech_tp_long.csv.",
        "",
        "  " + "  ".join(f"{k} {v:,}" for k, v in sorted(counts.items())),
    ]
    lines += lexeme_summaries(records)
    lines += rule_section(index)
    lines += competition_section(index)
    lines += nested_section(index)
    lines += spellout_section(index)
    lines += r17_section(index)
    lines += derivation_section(index)
    lines += learner_section(learner_rows)
    return "\n".join(lines) + "\n"
