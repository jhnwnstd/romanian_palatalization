"""Stage ``spellout``: spell-out TP on certain-specification roots (R20).

Gorman & Winstead's spell-out rules give the non-palatalizing member
of a competing pair only to roots prespecified as inalterable. A root
seen alternating is underspecified, so the account predicts the
palatalizing member for it. Czech hard feminine velars alternate in
DAT/LOC.SG, which has no non-palatalizing competitor (lize, bance,
knize, ruce), so their specification is not in doubt.

Roots
    hard feminine velars (gender f, velar stem, hard declension) with
    a DAT/LOC.SG form showing the change attested at τ; the dual
    remnant ruka is included (its singular is regular: ruce). A second
    set, ``E0-certain velars``, takes every velar noun root with an
    attested alternating inflected form in any cell.
N   roots with an attested derivative (at τ) from a competing pair
    (R20): ``-ní/-ný`` vs ``-ový/-ovní``, ``-it`` vs ``-ovat``, ``-ka``
    vs ``-yně``; ``any pair`` pools them (a root counts once per pair).
e   roots taking only the non-palatalizing member; roots with both
    are doublets coded by token share (R19: dominant ≥ 90%,
    doublet as exception, doublet excluded).

A member counts by its suffix (a -ní adjective without the change,
barokní, still chose -ní; those are listed). An -ovat verb whose stem
shows the change (zdražovat) is the imperfective of an -it verb and is
not an -ovat choice. Counting is :func:`lexpipe.tp.two_outcome_rows`
over observations whose ``alternates`` means "chose the palatalizing
member" (the encoding of the Romanian spell-out rows).

Which derivatives are choices between the members (:func:`member_side`,
audit 2026-10-23):

- ``-ka`` vs ``-yně`` are the two feminizing suffixes: only an animate
  masculine base can choose between them (Čech → Češka, žák → žačka /
  žákyně), and only a feminine agent counts (-yně, or a -ka glossed or
  linked as the base's feminine). A diminutive -ka (bleška, dečka) is
  no -yně alternative, so the pair is empty for the feminine set.
- ``-it`` vs ``-ovat``: a base derived from a verb
  (``base_deverbal_from``: záloha ← založit, sprcha ← spršet, stoka ←
  stékat) does not choose a verbalizer; its root is left out of this
  pair (the derivations stage already drops the reversed pairs).

Every row is also split by the base's R26 suffix target (subsets
``base:root-final`` and ``base:suffix``): the k of aktovka, baterka,
bankovka is the suffix -ka's, so their -ový choice says nothing about
a root's specification; the root-final rows (kniha, ruka, liga, droga)
are the test of the account.

Derivative evidence (:mod:`.derivations`) is reported per basis:
``wikt`` (main: Wiktionary links and Wiktionary headword candidates),
``linked`` (links only), ``wikt+corpus`` (adds corpus-only
candidates). Strata are Wiktionary's (R22) in the tracked report; the
section computed with czetyl and DeriNet strata and with DeriNet-only
pairs is labelled NC and holds aggregates only.

Outputs: ``analysis/czech/czech_spellout.txt`` (tracked aggregates;
example words without counts, their choice read at τ = dictionary so
that no per-word Leipzig attestation is tracked, R9) and
``data/czech/nc/czech_spellout_roots.csv`` (one row per root × pair
with counts).
"""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Final

from lexpipe import csvio, learners, tp
from lexpipe.attest import TAUS, attested_count
from lexpipe.context import RunContext, run_stage
from lexpipe.langs.czech import jsonio, leipzig, ortho
from lexpipe.langs.czech import suffixes as sx
from lexpipe.langs.czech.analysis_layout import AnalysisPaths, analysis_paths
from lexpipe.langs.czech.attestation import CORPUS_TAUS
from lexpipe.langs.czech.etymology import STRATA, UNKNOWN
from lexpipe.langs.czech.layout import CzechPaths, paths_for

STAGE: Final[str] = "spellout"
# 3: core attested_count; licence-gated report (its NC section counts
# DeriNet pairs and czetyl/DeriNet strata).
# 4: member validity (ka_yne feminine agents of m.an bases; it_ovat
# without deverbal bases); root-final vs suffix subsets; examples at
# tau = dictionary (audit 2026-10-23).
VERSION: Final[str] = "4"
BASE_SUBSETS: Final[tuple[str, ...]] = ("base:root-final", "base:suffix")
REPORT_SOURCES: Final[tuple[str, ...]] = (
    "wiktionary_en",
    "wiktionary_cs",
    "derinet",
    "czetyl",
)
LANG: Final[str] = "cs"
TABLE: Final[str] = "spellout"
ANY_PAIR: Final[str] = "any"
PAIRS: Final[tuple[str, ...]] = (*sx.PAIRS, ANY_PAIR)
REPORT_BASES: Final[tuple[str, ...]] = ("wikt", "linked", "wikt+corpus")
NC_BASIS: Final[str] = "all"
ROOT_SETS: Final[tuple[str, ...]] = ("feminine", "E0-certain")
EXAMPLES: Final[tuple[str, ...]] = (
    "liga",
    "banka",
    "droga",
    "kniha",
    "ruka",
    "řeka",
)
ROOT_FIELDS: Final[tuple[str, ...]] = (
    "root_set",
    "base_lexeme_id",
    "base_lemma",
    "base_suffix_target",
    "stratum_wikt",
    "stratum",
    "pair",
    "basis",
    "p_derivatives",
    "n_derivatives",
    "p_count",
    "n_count",
    *(f"certain_{tau}" for tau in TAUS),
    *(f"choice_{tau}" for tau in TAUS),
    "source",
    "licence_class",
)


@dataclass
class Base:
    """A candidate root of the spell-out TP (R4: root types).

    ``lexeme_id`` is the root's first lexeme (lexeme_id order), whose
    lemma and strata represent the root; ``certain`` and
    ``derivatives`` pool every lexeme of the root.
    """

    lexeme_id: str
    lemma: str
    stratum_wikt: str
    stratum: str
    certain: dict[str, bool] = field(default_factory=dict)
    derivatives: list[dict[str, str]] = field(default_factory=list)
    suffix_target: bool = False

    @property
    def base_subset(self) -> str:
        """``base:suffix`` when the root's target is an R26 suffix's
        (aktovka), else ``base:root-final`` (kniha)."""
        return BASE_SUBSETS[1] if self.suffix_target else BASE_SUBSETS[0]

    def absorb(
        self, flags: Mapping[str, bool], rows: Sequence[dict[str, str]]
    ) -> None:
        """Pool another lexeme of the same root."""
        _merge(self.certain, flags)
        seen = {(r["suffix"], r["derivative"]) for r in self.derivatives}
        for row in rows:
            if (row["suffix"], row["derivative"]) not in seen:
                seen.add((row["suffix"], row["derivative"]))
                self.derivatives.append(row)


# ---------------------------------------------------------------------------
# Certain roots
# ---------------------------------------------------------------------------


def _flags(row: Mapping[str, str]) -> dict[str, bool]:
    return {
        "dictionary": True,
        **{
            tau: csvio.parse_bool(row[f"attested_{tau}"]) is True
            for tau in CORPUS_TAUS
        },
    }


def _merge(target: dict[str, bool], flags: Mapping[str, bool]) -> None:
    for tau in TAUS:
        target[tau] = target.get(tau, False) or flags.get(tau, False)


def feminine_certain(paths: CzechPaths) -> dict[str, dict[str, bool]]:
    """``lexeme_id -> {tau: alternating DAT/LOC.SG attested}`` for the
    in-domain hard feminine velars (from ``czech_cells.csv``)."""
    out: dict[str, dict[str, bool]] = {}
    for row in jsonio.iter_csv(paths.cells):
        if (
            row["primary"] != "True"
            or row["gender_class"] != "f"
            or row["target_class"] != ortho.VELAR
            or row["cell"] not in ("DAT.SG", "LOC.SG")
        ):
            continue
        flags = out.setdefault(row["lexeme_id"], dict.fromkeys(TAUS, False))
        if row["allomorph_class"] == ortho.P and row["mutation"] == "True":
            _merge(flags, _flags(row))
    return out


def e0_certain(paths: CzechPaths) -> dict[str, dict[str, bool]]:
    """``lexeme_id -> {tau: any alternating P form attested}`` for the
    in-domain velar nouns."""
    out: dict[str, dict[str, bool]] = {}
    for row in jsonio.iter_csv(paths.cells):
        if (
            row["primary"] != "True"
            or row["pos"] != "noun"
            or row["target_class"] != ortho.VELAR
        ):
            continue
        flags = out.setdefault(row["lexeme_id"], dict.fromkeys(TAUS, False))
        if row["allomorph_class"] == ortho.P and row["mutation"] == "True":
            _merge(flags, _flags(row))
    return out


def dual_certain(
    paths: CzechPaths,
    lexemes: Mapping[str, Mapping[str, str]],
    ctx: RunContext,
) -> dict[str, dict[str, bool]]:
    """The same flags for out-of-domain dual remnants with a feminine
    velar stem (ruka): their cells are only in ``czech_cells_raw.csv``,
    so the change is checked with :func:`.ortho.align` and the forms
    are counted in the corpus here."""
    wanted = {
        lexeme_id: row
        for lexeme_id, row in lexemes.items()
        if row["ood_class"] == "dual_remnant"
        and row["target_class"] == ortho.VELAR
        and "f" in csvio.split_list(row["gender"])
    }
    if not wanted:
        return {}
    table = leipzig.attestation_table(leipzig.ensure_words_file(ctx, paths))
    out: dict[str, dict[str, bool]] = {}
    for row in jsonio.iter_csv(paths.cells_raw):
        lexeme = wanted.get(row["lexeme_id"])
        if (
            lexeme is None
            or row["primary"] != "True"
            or row["cell"] not in ("DAT.SG", "LOC.SG")
        ):
            continue
        flags = out.setdefault(row["lexeme_id"], dict.fromkeys(TAUS, False))
        alignment = ortho.align(lexeme["stem"], row["form"])
        if not (alignment.mutated and alignment.ending in ("e", "ě")):
            continue
        form = row["form"]
        count = table.count(form, case_sensitive=True) + table.count(
            form[:1].upper() + form[1:], case_sensitive=True
        )
        _merge(
            flags,
            {
                "dictionary": True,
                **{
                    tau: attested_count(count, tau, table.total_tokens)
                    for tau in CORPUS_TAUS
                },
            },
        )
    return out


# ---------------------------------------------------------------------------
# Choices
# ---------------------------------------------------------------------------


def valid_member(row: Mapping[str, str]) -> bool:
    """Is the pair row a choice between the pair's members (module
    docstring)? -ka/-yně: a feminine agent of an animate masculine
    base; -it/-ovat: a base not derived from a verb."""
    pair = row["suffix_pair"]
    if pair == sx.PAIR_KA_YNE:
        readings = csvio.split_list(row.get("base_readings", ""))
        agent = csvio.parse_bool(row.get("derivative_feminine_agent", ""))
        return "m.an" in readings and agent is True
    if pair == sx.PAIR_IT_OVAT:
        return not csvio.split_list(row.get("base_deverbal_from", ""))
    return True


def member_side(row: Mapping[str, str]) -> bool | None:
    """True for the palatalizing member of a pair, False for the other,
    None when the row is not a member choice (module docstring)."""
    if not row["suffix_pair"] or not valid_member(row):
        return None
    if row["suffix_class"] == sx.P_CLASS:
        return True
    if row["outcome"] == sx.PAL:
        return None  # -ovat imperfective of an -it verb (zdražovat)
    return False


def pair_observations(
    derivatives: Sequence[Mapping[str, str]], pair: str, basis_ok: Any
) -> tuple[tp.Observation, ...]:
    """Observations of one root in one pair (``any``: every pair)."""
    out: list[tp.Observation] = []
    for row in derivatives:
        if not basis_ok(row):
            continue
        if pair != ANY_PAIR and row["suffix_pair"] != pair:
            continue
        side = member_side(row)
        if side is None:
            continue
        out.append(
            tp.Observation(
                alternates=side,
                attested_by_tau=_flags(row),
                allomorph="P",
                form=row["derivative"],
                count=int(row["count"] or 0),
            )
        )
    return tuple(out)


def basis_test(basis: str) -> Any:
    """A row predicate for a derivative basis (nested, R14 doc)."""
    order = ("linked", "wikt", "wikt+corpus", "all")
    limit = order.index(basis)

    def test(row: Mapping[str, str]) -> bool:
        return order.index(row["basis"]) <= limit

    return test


def _spec(target: str, pair: str, tau: str, coding: str) -> tp.CellSpec:
    return tp.CellSpec(
        lang=LANG,
        table=TABLE,
        target_class=target,
        trigger=f"pair:{pair}",
        tau=tau,
        coding=coding,
    )


def rows_for(
    bases: Sequence[Base],
    root_set: str,
    basis: str,
    *,
    stratum_of: Any,
) -> list[dict[str, Any]]:
    """TP rows: pairs × τ × codings × strata (``all`` + each stratum)."""
    test = basis_test(basis)
    out: list[dict[str, Any]] = []
    for pair in PAIRS:
        for tau in TAUS:
            chosen = [b for b in bases if b.certain.get(tau, False)]
            observed = {
                b.lexeme_id: tp.RootObservation(
                    b.lexeme_id,
                    "velar",
                    0.0,
                    {
                        f"pair:{pair}": pair_observations(
                            b.derivatives, pair, test
                        )
                    },
                )
                for b in chosen
            }
            for coding in tp.CODINGS:
                spec = _spec(f"velar:{root_set}", pair, tau, coding)
                for subset in (
                    "all",
                    *BASE_SUBSETS,
                    *(f"stratum:{s}" for s in STRATA),
                ):
                    inside = [
                        observed[b.lexeme_id]
                        for b in chosen
                        if subset == "all"
                        or b.base_subset == subset
                        or f"stratum:{stratum_of(b)}" == subset
                    ]
                    for row in tp.two_outcome_rows(
                        inside, spec, subset=subset, memo=f"basis {basis}"
                    ):
                        out.append({**row, "basis": basis, "pair": pair})
    return out


# ---------------------------------------------------------------------------
# Stage
# ---------------------------------------------------------------------------


def run(ctx: RunContext) -> int | None:
    """Stage entry point (see the module docstring)."""
    paths = paths_for(ctx)
    apaths = analysis_paths(ctx)
    words = leipzig.ensure_words_file(ctx, paths)
    return run_stage(
        ctx,
        STAGE,
        inputs=[
            paths.lexicon_complete,
            paths.cells,
            paths.cells_raw,
            words,
            apaths.derivations,
            apaths.lexeme_strata,
        ],
        params={"version": VERSION},
        outputs=[apaths.spellout, apaths.spellout_roots],
        body=lambda: _spellout(ctx, paths, apaths),
    )


def load_bases(
    ctx: RunContext, paths: CzechPaths, apaths: AnalysisPaths
) -> dict[str, list[Base]]:
    """The two root sets with their derivatives and strata."""
    lexemes = {
        row["lexeme_id"]: row
        for row in jsonio.iter_csv(paths.lexicon_complete)
        if row["pos"] == "noun"
    }
    strata = {
        row["lexeme_id"]: (row["stratum_wikt"], row["stratum"])
        for row in jsonio.iter_csv(apaths.lexeme_strata)
    }
    certain = {
        "feminine": {
            **feminine_certain(paths),
            **dual_certain(paths, lexemes, ctx),
        },
        "E0-certain": e0_certain(paths),
    }
    derivatives: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in jsonio.iter_csv(apaths.derivations):
        derivatives[row["base_lexeme_id"]].append(row)
    out: dict[str, list[Base]] = {}
    for name, flags in certain.items():
        by_root: dict[str, Base] = {}
        for lexeme_id in sorted(flags):
            lexeme = lexemes.get(lexeme_id)
            if lexeme is None:
                continue
            rows = derivatives.get(lexeme_id, [])
            root = by_root.get(lexeme["root_id"])
            if root is not None:
                root.absorb(flags[lexeme_id], rows)
                continue
            wikt, combined = strata.get(lexeme_id, (UNKNOWN, UNKNOWN))
            by_root[lexeme["root_id"]] = Base(
                lexeme_id=lexeme_id,
                lemma=lexeme["lemma"],
                stratum_wikt=wikt or UNKNOWN,
                stratum=combined or UNKNOWN,
                certain=dict(flags[lexeme_id]),
                derivatives=list(rows),
                suffix_target=csvio.parse_bool(lexeme["target_is_suffix"])
                is True,
            )
        out[name] = sorted(by_root.values(), key=lambda b: b.lexeme_id)
    return out


def _spellout(
    ctx: RunContext, paths: CzechPaths, apaths: AnalysisPaths
) -> int:
    sets = load_bases(ctx, paths, apaths)
    tracked: list[dict[str, Any]] = []
    nc: list[dict[str, Any]] = []
    for name in ROOT_SETS:
        for basis in REPORT_BASES:
            tracked += rows_for(
                sets[name], name, basis, stratum_of=lambda b: b.stratum_wikt
            )
        nc += rows_for(
            sets[name], name, NC_BASIS, stratum_of=lambda b: b.stratum
        )
    csvio.write_rows(
        apaths.spellout_roots, root_rows(sets), ROOT_FIELDS, ctx=ctx
    )
    text = render_report(tracked, nc, sets)
    ctx.write_tracked_text(
        apaths.spellout, text, sources=REPORT_SOURCES, nc_aggregates=True
    )
    print(text, flush=True)
    return 0


def choice(observations: Sequence[tp.Observation], tau: str) -> str:
    """``P``, ``N``, ``both`` or "" for one root, pair and τ."""
    sides = {
        bool(obs.alternates)
        for obs in observations
        if obs.attested_by_tau.get(tau, False)
    }
    if sides == {True, False}:
        return "both"
    if sides == {True}:
        return "P"
    if sides == {False}:
        return "N"
    return ""


def root_rows(sets: Mapping[str, Sequence[Base]]) -> list[dict[str, Any]]:
    """The NC per-root file (counts included)."""
    out: list[dict[str, Any]] = []
    test = basis_test(NC_BASIS)
    for name in ROOT_SETS:
        for base in sets[name]:
            for pair in PAIRS:
                observations = pair_observations(base.derivatives, pair, test)
                if not observations:
                    continue
                rows = [
                    r
                    for r in base.derivatives
                    if member_side(r) is not None
                    and (pair == ANY_PAIR or r["suffix_pair"] == pair)
                ]
                out.append(
                    {
                        "root_set": name,
                        "base_lexeme_id": base.lexeme_id,
                        "base_lemma": base.lemma,
                        "base_suffix_target": base.suffix_target,
                        "stratum_wikt": base.stratum_wikt,
                        "stratum": base.stratum,
                        "pair": pair,
                        "basis": sorted({r["basis"] for r in rows}),
                        "p_derivatives": sorted(
                            r["derivative"] for r in rows if member_side(r)
                        ),
                        "n_derivatives": sorted(
                            r["derivative"] for r in rows if not member_side(r)
                        ),
                        "p_count": sum(
                            int(r["count"] or 0)
                            for r in rows
                            if member_side(r)
                        ),
                        "n_count": sum(
                            int(r["count"] or 0)
                            for r in rows
                            if not member_side(r)
                        ),
                        **{
                            f"certain_{tau}": base.certain.get(tau, False)
                            for tau in TAUS
                        },
                        **{
                            f"choice_{tau}": choice(observations, tau)
                            for tau in TAUS
                        },
                        "source": [
                            "wiktionary_en",
                            "wiktionary_cs",
                            "derinet",
                            "czetyl",
                            "leipzig",
                        ],
                        "licence_class": "by-nc-sa",
                    }
                )
    return out


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

SHOW_TAUS: Final[tuple[str, ...]] = ("dictionary", "ge1", "ge2pm")
PAIR_NAMES: Final[Mapping[str, str]] = {
    **sx.PAIR_LABELS,
    ANY_PAIR: "any pair",
}

Key = tuple[str, str, str, str, str, str, str]


def _index(rows: Sequence[Mapping[str, Any]]) -> dict[Key, Mapping[str, Any]]:
    return {
        (
            str(r["target_class"]),
            str(r["basis"]),
            str(r["pair"]),
            str(r["tau"]),
            str(r["coding"]),
            str(r["subset"]),
            str(r["rule"]),
        ): r
        for r in rows
    }


def _theta(value: Any) -> str:
    if value is None:
        return "-"
    number = float(value)
    return f"{number:.1f}" if math.isfinite(number) else "-"


def _tol(value: Any) -> str:
    return "-" if value is None else ("yes" if value else "no")


def _line(
    row: Mapping[str, Any] | None, retain: Mapping[str, Any] | None
) -> str:
    if row is None:
        return "   (no roots)"
    frag = row.get("fragility")
    ret = "-" if retain is None else str(retain["e"])
    return (
        f"{row['N_target']:>6} {row['N_trigger']:>5} {row['e_nonalt']:>5} "
        f"{row['e']:>5} {_theta(row['theta_trigger']):>7} "
        f"{_tol(row['tolerated_trigger']):>4} {'-' if frag is None else frag:>5}"
        f" {ret:>6}"
    )


HEAD: Final[str] = (
    f"{'roots':>6} {'N':>5} {'Nonly':>5} {'e':>5} {'theta':>7} {'tol':>4} "
    f"{'frag':>5} {'e(ret)':>6}"
)


def _with_fragility(
    rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    out = []
    for row in rows:
        item = dict(row)
        item["fragility"] = learners.fragility(
            int(row["N_trigger"]), int(row["e"])
        )
        out.append(item)
    return out


def _table(
    index: Mapping[Key, Mapping[str, Any]],
    target: str,
    basis: str,
    *,
    taus: Sequence[str] = SHOW_TAUS,
    codings: Sequence[str] = ("dominant",),
    subsets: Sequence[str] = ("all",),
) -> list[str]:
    lines = [
        f"  {'pair':<22} {'tau':<10} {'coding':<18} {'subset':<28} {HEAD}"
    ]
    for pair in PAIRS:
        for tau in taus:
            for coding in codings:
                for subset in subsets:
                    key = (target, basis, pair, tau, coding, subset)
                    row = index.get((*key, "palatalize"))
                    retain = index.get((*key, "retain"))
                    if row is None or not row["N_trigger"] and subset != "all":
                        continue
                    lines.append(
                        f"  {PAIR_NAMES[pair]:<22} {tau:<10} {coding:<18} "
                        f"{subset:<28} {_line(row, retain)}"
                    )
    return lines


def _examples(sets: Mapping[str, Sequence[Base]]) -> list[str]:
    lines = [
        "",
        "Examples (feminine set; basis wikt; tau = dictionary, so that no",
        "per-word corpus attestation is printed, R9): stratum (Wiktionary),",
        "then per pair the choice and the derivatives (P member / N member)",
    ]
    test = basis_test("wikt")
    for lemma in EXAMPLES:
        found = [b for b in sets["feminine"] if b.lemma == lemma]
        if not found:
            lines.append(f"  {lemma:<8} not among the certain roots")
            continue
        base = found[0]
        parts = []
        for pair in sx.PAIRS:
            observations = pair_observations(base.derivatives, pair, test)
            if not observations:
                continue
            rows = [
                r
                for r in base.derivatives
                if test(r)
                and r["suffix_pair"] == pair
                and member_side(r) is not None
            ]
            p_words = sorted({r["derivative"] for r in rows if member_side(r)})
            n_words = sorted(
                {r["derivative"] for r in rows if member_side(r) is False}
            )
            parts.append(
                f"{PAIR_NAMES[pair]}: "
                f"{choice(observations, 'dictionary') or '-'} "
                f"({', '.join(p_words) or '-'} / {', '.join(n_words) or '-'})"
            )
        lines.append(
            f"  {lemma:<8} {base.stratum_wikt:<20} "
            f"{base.base_subset[5:]}; "
            + ("; ".join(parts) or "no pair derivative")
        )
    return lines


def _special(sets: Mapping[str, Sequence[Base]]) -> list[str]:
    exceptions: set[str] = set()
    imperfective: set[str] = set()
    for base in sets["feminine"]:
        for row in base.derivatives:
            if not row["suffix_pair"] or row["basis"] not in (
                "linked",
                "wikt",
            ):
                continue
            if (
                row["suffix_class"] == sx.P_CLASS
                and row["outcome"] == sx.EXCEPTION
            ):
                exceptions.add(f"{base.lemma} → {row['derivative']}")
            if valid_member(row) and member_side(row) is None:
                imperfective.add(f"{base.lemma} → {row['derivative']}")
    return [
        "",
        "Palatalizing member chosen without the change (counted as P choices):",
        "  " + (", ".join(sorted(exceptions)) or "none"),
        "-ovat verbs with the changed stem (imperfectives of -it verbs; not",
        "counted as -ovat choices):",
        "  " + (", ".join(sorted(imperfective)) or "none"),
    ]


def render_report(
    tracked: Sequence[Mapping[str, Any]],
    nc: Sequence[Mapping[str, Any]],
    sets: Mapping[str, Sequence[Base]],
) -> str:
    """The report text (module docstring)."""
    index = _index(_with_fragility(tracked))
    nc_index = _index(_with_fragility(nc))
    lines = [
        "Czech spell-out TP on certain-specification roots (design R20)",
        "=" * 72,
        "",
        "Question: do roots whose specification is certain (they alternate",
        "in their own paradigm) take the palatalizing member of a competing",
        "pair, as the spell-out account predicts? Unit: root types.",
        "  roots  certain roots with any derivative of the pair at tau",
        "  N      roots with an attested derivative from the pair",
        "  Nonly  roots taking only the non-palatalizing member",
        "  e      exceptions of the palatalize (spell-out) rule: Nonly plus",
        "         doublets by coding (R19); theta = N/ln N",
        "  e(ret) exceptions of the opposite rule (non-palatalizing member)",
        "Feminine set: hard feminine velars with an alternating DAT/LOC.SG",
        "attested at tau (plus the dual remnant ruka). E0-certain set: velar",
        "nouns with any attested alternating inflected form.",
        "Derivatives (basis): wikt = Wiktionary links + Wiktionary headword",
        "candidates (main); linked = links only; wikt+corpus = adds",
        "corpus-only candidates (R21). Strata: Wiktionary (R22).",
        "",
        "Certain roots by set and tau:",
    ]
    for name in ROOT_SETS:
        counts = "  ".join(
            f"{tau} {sum(1 for b in sets[name] if b.certain.get(tau)):,}"
            for tau in TAUS
        )
        lines.append(f"  {name:<11} {counts}")
    lines += [
        "",
        "Members (module docstring): -ka vs -yne only for animate masculine",
        "bases and feminine agents (empty for the feminine set: its -ka",
        "derivatives are diminutives); -it vs -ovat without bases derived",
        "from a verb. Subsets base:root-final / base:suffix split the roots",
        "by their R26 suffix target (aktovka's k is the suffix -ka's);",
        "the root-final rows test root specification.",
        "",
        "1. Feminine set, basis wikt, all codings",
        "-" * 72,
        *_table(
            index,
            "velar:feminine",
            "wikt",
            codings=("dominant", "doublet_exception", "doublet_excluded"),
        ),
        "",
        "1b. Feminine set by R26 suffix target (basis wikt; dominant)",
        "-" * 72,
        *_table(
            index,
            "velar:feminine",
            "wikt",
            taus=("dictionary", "ge1", "ge2pm"),
            subsets=BASE_SUBSETS,
        ),
        "",
        "2. Feminine set, other bases (dominant coding)",
        "-" * 72,
    ]
    for basis in REPORT_BASES:
        if basis == "wikt":
            continue
        lines.append(f"  basis {basis}")
        lines += _table(index, "velar:feminine", basis)
    lines += [
        "",
        "3. Feminine set by base stratum (Wiktionary; basis wikt; dominant)",
        "-" * 72,
        *_table(
            index,
            "velar:feminine",
            "wikt",
            taus=("dictionary", "ge1"),
            subsets=tuple(f"stratum:{s}" for s in STRATA),
        ),
        "",
        "4. E0-certain velars (any gender), basis wikt, dominant",
        "-" * 72,
        *_table(index, "velar:E0-certain", "wikt"),
        "",
        "4b. E0-certain velars by R26 suffix target (basis wikt; dominant)",
        "-" * 72,
        *_table(
            index,
            "velar:E0-certain",
            "wikt",
            taus=("dictionary", "ge1", "ge2pm"),
            subsets=BASE_SUBSETS,
        ),
        "",
        "5. E0-certain velars by base stratum (Wiktionary), dominant",
        "-" * 72,
        *_table(
            index,
            "velar:E0-certain",
            "wikt",
            taus=("dictionary", "ge1"),
            subsets=tuple(f"stratum:{s}" for s in STRATA),
        ),
    ]
    lines += _examples(sets)
    lines += _special(sets)
    lines += [
        "",
        "6. NC sources (aggregates only): basis all (DeriNet-only pairs",
        "added), strata from Wiktionary > czetyl > DeriNet (R22); dominant",
        "-" * 72,
    ]
    for target in ("velar:feminine", "velar:E0-certain"):
        lines.append(f"  {target}")
        lines += _table(
            nc_index,
            target,
            NC_BASIS,
            taus=("dictionary", "ge1"),
            subsets=("all", *(f"stratum:{s}" for s in STRATA)),
        )
    return "\n".join(lines) + "\n"
