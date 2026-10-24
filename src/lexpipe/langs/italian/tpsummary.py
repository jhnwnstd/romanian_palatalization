"""The readable digest of the Italian TP rows (``italian_tp_summary.txt``).

Pure formatting of the rows built by :mod:`.tptables`,
:mod:`.tplearners` and :mod:`.spellout`: where the dictionary's
plurals come from (Module:it-common defaults vs typed plurals,
:func:`.roots.dictionary_provenance`), the key verdicts (pooled -i
palatalization vs retention, the non-trigger -e cell, the coronal rows,
the stress-nested model fit with R26's suffixes counted once, suffix
rows, verbs by conjugation, the pooled and category-indexed noun + verb
rows, evidence sets, doublet codings, the V sweep and the simulated
learners). Aggregates only (R9).
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import Any, Final

from lexpipe import suffixoids

Row = Mapping[str, Any]
# The R27 classes of section 9: the shared ones and Italian adjectives.
R27_CLASSES: Final[tuple[str, ...]] = (*suffixoids.CLASSES, "adjective")
_KEY_FIELDS: Final[tuple[str, ...]] = (
    "table",
    "target_class",
    "trigger",
    "subset",
    "tau",
    "evidence_set",
    "coding",
    "freq_basis",
    "V",
    "rule",
)


CONDITIONED: Final[str] = (
    "every rule tolerated: the conditioned grammar is productive"
)
NOT_CONDITIONED: Final[str] = "not every rule tolerated"


class Index:
    """Rows by their identifying fields (first match wins)."""

    def __init__(self, rows: Iterable[Row]) -> None:
        self._rows: dict[tuple[Any, ...], Row] = {}
        self._memo: dict[tuple[Any, ...], list[Row]] = {}
        for row in rows:
            key = tuple(row.get(f) for f in _KEY_FIELDS) + (row.get("memo"),)
            self._rows.setdefault(key, row)
            self._memo.setdefault(key[:-1], []).append(row)

    def get(self, memo: str | None = None, **where: Any) -> Row | None:
        """The row with these fields (V None, lemma basis, dominant
        coding and E0 unless given)."""
        defaults: dict[str, Any] = {
            "evidence_set": "E0",
            "coding": "dominant",
            "freq_basis": "lemma",
            "V": None,
        }
        defaults.update(where)
        key = tuple(defaults.get(f) for f in _KEY_FIELDS)
        found = self._memo.get(key, [])
        if memo is None:
            return found[0] if found else None
        for row in found:
            if str(row.get("memo", "")).endswith(memo):
                return row
        return None


def _num(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.1f}"
    return str(value)


def pair_line(index: Index, label: str, **where: Any) -> str:
    """palatalize and retain of one cell on one line."""
    pal = index.get(rule="palatalize", **where)
    ret = index.get(rule="retain", **where)
    if pal is None or ret is None:
        return f"  {label:<44} (no row)"
    verdict = _verdict(pal, ret)
    return (
        f"  {label:<44} N_target {pal['N_target']:>5}  N {pal['N_trigger']:>5}"
        f"  e_pal {pal['e']:>5}  e_ret {ret['e']:>5}  theta "
        f"{_num(pal['theta_trigger']):>6}  -> {verdict}"
    )


def _verdict(pal: Row, ret: Row) -> str:
    if pal["tolerated_trigger"] is None:
        return "empty cell (N <= 1)"
    good_pal = pal["tolerated_trigger"] is True
    good_ret = ret["tolerated_trigger"] is True
    if good_pal and good_ret:
        return "both tolerated"
    if good_pal:
        return f"palatalize productive (fragility {pal['fragility']})"
    if good_ret:
        return f"retain productive (fragility {ret['fragility']})"
    return "neither productive"


def nested_lines(
    rows: Sequence[Row],
    label: str,
    title: str = "",
    **where: Any,
) -> list[str]:
    """The rows of one nested set (*title* defaults to its name)."""
    chosen = [
        r
        for r in rows
        if all(r.get(k) == v for k, v in where.items())
        and str(r.get("subset", "")).startswith("nested:")
        and str(r.get("memo", "")).endswith(label)
    ]
    out = [f"  nested set {title or label}:"]
    for row in chosen:
        out.append(
            f"    {row['subset']:<26} {row['rule']:<10} N {row['N_trigger']:>5}"
            f"  e {row['e']:>5}  theta {_num(row['theta_trigger']):>6}  "
            f"tolerated {_num(row['tolerated_trigger'])}"
        )
    tolerated = chosen and all(r["tolerated_trigger"] is True for r in chosen)
    out.append(f"    -> {CONDITIONED if tolerated else NOT_CONDITIONED}")
    return out


def learner_lines(
    learners: Sequence[Row], V_values: Sequence[int]
) -> list[str]:
    """Learner shares at the given V values."""
    out = [
        f"  {'table':<10} {'trigger':<6} {'subset':<13} {'E':<3} {'tau':<10} "
        f"{'V':>5}  palatalize  retain  neither  conditioned"
    ]
    for row in learners:
        if row["V"] not in V_values:
            continue
        out.append(
            f"  {row['table']:<10} {row['trigger']:<6} {row['subset']:<13} "
            f"{row['evidence_set']:<3} {row['tau']:<10} {row['V']:>5}  "
            f"{row['share_palatalize']:>10.3f}"
            f"  {row['share_retain']:>6.3f}  {row['share_neither']:>7.3f}"
            f"  {row['share_conditioned']:>11.3f}"
        )
    return out


def _share(part: int, whole: int) -> str:
    return (
        f"{part:,} of {whole:,} ({100 * part / whole:.1f}%)" if whole else "0"
    )


def provenance_lines(counts: Mapping[str, int]) -> list[str]:
    """Section 0: where the dictionary's dorsal -i plurals come from."""
    out: list[str] = []
    for scope, label in (
        ("", "listed (tau = dictionary)"),
        ("ge1:", "attested (tau >= 1 token)"),
    ):
        forms = counts.get(f"{scope}forms", 0)
        default = counts.get(f"{scope}module_default", 0)
        out += [
            f" dorsal -i forms {label}: {forms:,}",
            f"  Module:it-common defaults: {_share(default, forms)}; "
            f"typed by an editor: {_share(counts.get(f'{scope}explicit', 0), forms)}",
            "  defaults that follow the module's rule (-ico/-aco/-ologo "
            "palatal, other -co/-go velar): "
            + _share(
                counts.get(f"{scope}module_default:follows_rule", 0), default
            ),
            "  defaults whose competing spelling (amichi for amici), not a "
            "form of another word, has a lowercase token: "
            + _share(
                counts.get(f"{scope}module_default:competitor_attested", 0),
                default,
            )
            + "; more tokens than the listed form: "
            + f"{counts.get(f'{scope}module_default:competitor_more_frequent', 0):,}",
        ]
    out += [
        " At tau = dictionary the rows mostly restate the module's spelling",
        " rule (Module:it-common-predicted); the rows at tau >= 1 rest on",
        " attestation, and tau = dictionary_explicit counts typed plurals",
        " only (a sample biased toward the module's exceptions).",
    ]
    return out


def summary(
    rows: Sequence[Row],
    learners: Sequence[Row],
    counts: Mapping[str, int],
    n_samples: int,
    seed: int,
    provenance: Mapping[str, int] | None = None,
) -> str:
    """The digest (module docstring); *provenance* is
    :func:`.roots.dictionary_provenance` (section 0, left out without
    it)."""
    index = Index(rows)
    taus = ("dictionary", "ge1", "ge0.5pm", "ge2pm")
    lines = [
        "Italian Tolerance Principle summary (design R16-R19, R24)",
        "=" * 72,
        "",
        "theta = N/ln N; a rule is productive iff e <= theta. N = roots",
        "attested before the trigger at tau; e_pal = roots that keep the",
        "velar (exceptions to palatalize), e_ret = roots that palatalize",
        "(exceptions to retain). Unit: root types (R4). E0 = nominal",
        "inflection, E1 = + verb cells, E2 = + derivatives (R15).",
        f"Roots: {counts.get('roots', 0):,}; E0 realizations "
        f"{counts.get('E0', 0):,}, E1 verb cells on noun roots "
        f"{counts.get('E1', 0):,}, E2 derivatives "
        f"{counts.get('E2', 0):,}; verb units {counts.get('verb_units', 0):,}.",
        "All rows: italian_tp_long.csv; learners: italian_learners.csv.",
        "",
    ]
    if provenance is not None:
        lines += [
            "0. Where the dictionary's dorsal -i plurals come from (R11)",
            "-" * 72,
            *provenance_lines(provenance),
            "",
        ]
    lines += [
        "1. Cross-language rule family (R17), E0, subset all, dominant coding",
        "-" * 72,
    ]
    for tau in taus:
        lines.append(f" tau = {tau}")
        lines.append(
            pair_line(
                index,
                "R17a dorsal/_[+front V] (-i/-e)",
                table="R17a",
                target_class="dorsal",
                trigger="-i/-e",
                subset="all",
                tau=tau,
            )
        )
        lines.append(
            pair_line(
                index,
                "R17a, doublet_exception coding",
                table="R17a",
                target_class="dorsal",
                trigger="-i/-e",
                subset="all",
                tau=tau,
                coding="doublet_exception",
            )
        )
        lines.append(
            pair_line(
                index,
                "R17b dorsal/_[+high,+front] (-i)",
                table="R17b",
                target_class="dorsal",
                trigger="-i",
                subset="all",
                tau=tau,
            )
        )
        lines.append(
            pair_line(
                index,
                "R17b -i, nouns only",
                table="R17b",
                target_class="dorsal",
                trigger="-i",
                subset="nouns",
                tau=tau,
            )
        )
        lines.append(
            pair_line(
                index,
                "R17b -i, adjectives only",
                table="R17b",
                target_class="dorsal",
                trigger="-i",
                subset="adjectives",
                tau=tau,
            )
        )
        lines.append(
            pair_line(
                index,
                "R17b -i, core (no -ico, no suffixoid)",
                table="R17b",
                target_class="dorsal",
                trigger="-i",
                subset="core",
                tau=tau,
            )
        )
        lines.append(
            pair_line(
                index,
                "R17b non-trigger -e (feminine plurals)",
                table="R17b",
                target_class="dorsal",
                trigger="-e",
                subset="all",
                tau=tau,
            )
        )
        lines.append(
            pair_line(
                index,
                "R17c coronal/_[+high,+front] (-i)",
                table="R17c",
                target_class="coronal",
                trigger="-i",
                subset="all",
                tau=tau,
            )
        )
        lines.append(
            pair_line(
                index,
                "R17c non-trigger -e",
                table="R17c",
                target_class="coronal",
                trigger="-e",
                subset="all",
                tau=tau,
            )
        )
    lines += [
        "",
        "2. Evidence sets and doublet codings (R17b -i, subset all)",
        "-" * 72,
    ]
    for tau in ("dictionary", "ge1"):
        for evidence in ("E0", "E1", "E2"):
            for coding in (
                "dominant",
                "doublet_exception",
                "doublet_excluded",
            ):
                lines.append(
                    pair_line(
                        index,
                        f"{tau} {evidence} {coding}",
                        table="R17b",
                        target_class="dorsal",
                        trigger="-i",
                        subset="all",
                        tau=tau,
                        evidence_set=evidence,
                        coding=coding,
                    )
                )
        lines.append(
            pair_line(
                index,
                f"{tau} E2_all dominant (+ borrowed, inherited)",
                table="R17b",
                target_class="dorsal",
                trigger="-i",
                subset="all",
                tau=tau,
                evidence_set="E2_all",
            )
        )
    lines += [
        "",
        "3. V sweep (R24 frame: V most frequent noun roots; R17b -i)",
        "-" * 72,
    ]
    for basis in ("lemma", "paradigm"):
        for V in (500, 1000, 2000, 5000):
            for tau in ("dictionary", "ge1"):
                lines.append(
                    pair_line(
                        index,
                        f"{basis} V={V} {tau} nouns",
                        table="R17b",
                        target_class="dorsal",
                        trigger="-i",
                        subset="nouns",
                        tau=tau,
                        freq_basis=basis,
                        V=V,
                    )
                )
    lines += [
        "",
        "4. Italian model fit: height-indexed -i nested by stress (it_stress)",
        "-" * 72,
    ]
    for tau in ("dictionary", "ge1"):
        lines.append(f" tau = {tau}, E0, dorsal")
        for subset in (
            "stress:post_tonic",
            "stress:not_post_tonic",
            "stress:unknown",
        ):
            lines.append(
                pair_line(
                    index,
                    subset,
                    table="it_stress",
                    target_class="dorsal",
                    trigger="-i",
                    subset=subset,
                    tau=tau,
                    memo="height-indexed /i j/ trigger",
                )
            )
        for name in (
            "giavazzi",
            "suffix_once",
            "orth_ico_once",
            "suffix_then_stress",
            "reverse",
        ):
            lines += nested_lines(
                rows,
                name,
                table="it_stress",
                target_class="dorsal",
                tau=tau,
                evidence_set="E0",
            )
        lines.append(f" tau = {tau}, E0, k and g separately")
        for target in ("c", "g"):
            lines += nested_lines(
                rows,
                "giavazzi",
                f"giavazzi, target {target}",
                table="it_stress",
                target_class=target,
                tau=tau,
                evidence_set="E0",
            )
        lines.append(f" tau = {tau}, E1 (verb cells added), dorsal")
        for name in ("giavazzi", "suffix_once"):
            lines += nested_lines(
                rows,
                name,
                table="it_stress",
                target_class="dorsal",
                tau=tau,
                evidence_set="E1",
            )
    lines.append(" tau = dictionary_explicit (typed plurals only), E0, dorsal")
    for subset in (
        "stress:post_tonic",
        "stress:not_post_tonic",
        "stress:unknown",
    ):
        lines.append(
            pair_line(
                index,
                subset,
                table="it_stress",
                target_class="dorsal",
                trigger="-i",
                subset=subset,
                tau="dictionary_explicit",
            )
        )
    for name in ("giavazzi", "suffix_once", "orth_ico_once"):
        lines += nested_lines(
            rows,
            name,
            table="it_stress",
            target_class="dorsal",
            tau="dictionary_explicit",
            evidence_set="E0",
        )
    lines += ["", "5. Suffix rows (R26; E0, -i)", "-" * 72]
    for tau in ("dictionary", "ge1"):
        for subset in (
            "suffix:-ico",
            "suffix:-iaco",
            "suffixoid:-logo",
            "suffixoid:-fago",
            "suffixoid:-urgo",
            "root",
            "orth:-ico",
            "orth:-Vico",
            "orth:-Vico post-tonic",
            "deverbal_from_are",
            "shares_1st_conj_verb",
        ):
            lines.append(
                pair_line(
                    index,
                    f"{tau} {subset}",
                    table="it_suffix",
                    target_class="dorsal",
                    trigger="-i",
                    subset=subset,
                    tau=tau,
                )
            )
    lines += [
        "",
        "6. Verbs by conjugation (verb units, E1 cells; it_verbs)",
        "-" * 72,
    ]
    for tau in ("dictionary", "ge1"):
        for subset in (
            "conj:are",
            "conj:ere_root",
            "conj:ere_stressed",
            "conj:ire",
            "conj:isc",
            "conj:dire",
            "overapplying",
            "all",
        ):
            for trigger in ("-i", "-e"):
                lines.append(
                    pair_line(
                        index,
                        f"{tau} {subset} {trigger}",
                        table="it_verbs",
                        target_class="dorsal",
                        trigger=trigger,
                        subset=subset,
                        tau=tau,
                        evidence_set="E1",
                    )
                )
    lines += [
        "",
        "7. Pooled noun + verb row and category-indexed rules (it_nounverb, -i)",
        "-" * 72,
    ]
    for tau in ("dictionary", "ge1"):
        lines.append(
            pair_line(
                index,
                f"{tau} pooled",
                table="it_nounverb",
                target_class="dorsal",
                trigger="-i",
                subset="pooled",
                tau=tau,
                evidence_set="E1",
            )
        )
        for name in ("category", "category_stress"):
            lines += nested_lines(
                rows, name, table="it_nounverb", trigger="-i", tau=tau
            )
    lines += [
        "",
        "8. Structural classes (nde; E0, dictionary and ge1)",
        "-" * 72,
    ]
    for target in ("nde:occhio", "nde:psiche", "frontstem"):
        for tau in ("dictionary", "ge1"):
            lines.append(
                pair_line(
                    index,
                    f"{target} -i {tau}",
                    table="nde",
                    target_class=target,
                    trigger="-i",
                    subset="all",
                    tau=tau,
                )
            )
    lines += ["", "9. -logo/-fago/-urgo (R27; E0, -i)", "-" * 72]
    for suffixoid in ("-logo", "-fago", "-urgo"):
        for klass in R27_CLASSES:
            for tau in ("dictionary", "ge1"):
                lines.append(
                    pair_line(
                        index,
                        f"{suffixoid} {klass} {tau}",
                        table="logtable",
                        target_class="dorsal",
                        trigger="-i",
                        subset=f"{suffixoid}:{klass}",
                        tau=tau,
                    )
                )
    lines += ["", "10. Spell-out TP (R20; italian_spellout.txt)", "-" * 72]
    for pair in (
        "verbalizer",
        "verbalizer_no_are",
        "ista_ismo_ita",
        "learned",
    ):
        for tau in ("dictionary", "ge1"):
            lines.append(
                pair_line(
                    index,
                    f"{pair} {tau} (spell-out default = palatalize)",
                    table="spellout",
                    target_class="dorsal",
                    trigger=pair,
                    subset="all",
                    tau=tau,
                )
            )
    lines += [
        "",
        f"11. Simulated learners (R18; {n_samples} vocabularies per V, "
        f"seed {seed}; lemma basis)",
        "-" * 72,
    ]
    lines += learner_lines(learners, (1000, 5000))
    return "\n".join(lines) + "\n"
