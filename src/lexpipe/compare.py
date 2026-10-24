"""Cross-language synthesis: Romanian (harmonized), Italian and Czech.

Joins the tracked outputs of the three pipelines into one table,
``analysis/comparison/crosslang_summary.csv``, and writes the generated
block of ``analysis/comparison/README.md`` (the tables of the
cross-language report; the prose around them is written by hand and
kept byte for byte, like the Romanian section that
:mod:`lexpipe.langs.romanian_compare` maintains in the same file).

Inputs
------
Only tracked files are read (:class:`Inputs`), so the summary can be
rebuilt from a checkout without the gitignored ``data/``:

  - Romanian: ``analysis/comparison/romanian_harmonized_tp_long.csv``,
    ``romanian_learners.csv``, ``romanian_spellout.txt``,
    ``romanian_noevidence.txt``, ``romanian_covariation.txt`` (the
    harmonized recount, design §7), and the legacy pipeline's
    ``analysis/romanian_tp_summary_{downsampled,full}.csv`` and
    ``analysis/palatalization_rule_validation.txt`` (read only, R1);
  - Italian: ``analysis/italian/italian_tp_long.csv``,
    ``italian_learners.csv``, ``italian_noevidence.txt``,
    ``italian_covariation.txt``, ``italian_rule_validation*.txt``;
  - Czech: ``analysis/czech/czech_tp_long.csv``, ``czech_learners.csv``,
    ``czech_spellout.txt``, ``czech_noevidence.txt``,
    ``czech_covariation.txt``, ``czech_logtable.csv``,
    ``czech_logtable_summary.txt``, ``czech_rule_validation*.txt``.

CSV rows are looked up by their columns; text reports are parsed
strictly: a section or line that is not where the report's format puts
it raises :class:`ReportFormatError` naming the file, so a change of a
report's format cannot silently drop or shift a number. Only the
Wiktionary-basis parts of the Czech reports are read: their NC sections
(DeriNet pairs, czetyl strata) are never parsed, and no row of the
summary rests on an NC source (R9).

Rows (``section`` column, in this order)
----------------------------------------
``r17``
    the cross-language rule family (R17): R17a dorsal/_[+front V] over
    -i and -e, R17b dorsal/_[+high,+front] over -i and its non-trigger
    -e cell, R17c coronal/_[+high,+front] over -i and -e for every
    coronal class a language has (Romanian t+d+s and z; Italian
    t/d/s/z; Czech t/d/n, r, s/z, l). One row per cell: both rules'
    exceptions (R18), θ, both verdicts, fragility, the shielded roots
    and the counterfactual θ(N_target).
``learners``
    simulated learners (R18): shares of 1000 learners at V = 1000 and
    V = 5000 (E0, plus Italian and Romanian E1).
``model_fit``
    the best language-specific fit: Romanian separate dorsal and
    coronal rules (the coronal rule nested: -e plural → retain,
    elsewhere → palatalize); Italian height-indexed -i nested by stress
    (giavazzi, R26 suffix_once, orth_ico_once; τ = dictionary_explicit
    too), the noun + verb category_stress set and the verbs by
    conjugation; Czech suffix-indexed rule TP (a), the nested NOM.PL
    competition and the R22 derivational suffix rows. Nested branches
    are single-rule rows (``rule``); each nested set closes with a
    ``grammar`` row whose verdict is ``conditioned`` when every branch
    is tolerated.
``spellout``
    R20 on certain-specification roots (Romanian -i vs -a/-ui;
    Italian -izzare vs -eggiare/-ire(/-are), palatal vs velar
    -ista/-ismo/-ità; Czech -ní/-ný vs -ový/-ovní, -it vs -ovat, -ka vs
    -yně, split by R26 suffix target), and the contrast sets
    (non-alternating roots).
``noevidence``
    R21: base class (alternating / non-alternating / no evidence) ×
    derivative palatalized or not, and the roots per base class.
``covariation``
    R23: the 2×2 base class × derivative choice, its sample odds ratio
    (Woolf interval) and Fisher p computed here with
    :mod:`lexpipe.stats` from the 2×2 (for a row pooled over pairs or
    suffixes, from the summed table), and as reported by each language:
    the Mantel–Haenszel OR over its strata, the stratified permutation
    p and the logistic ``base_class`` term with its status.
``logtable``
    R27: -log/-fag/-urg nouns by person/thing class.
``phonology``
    the Logical Phonology validators: rows validated and matched, per
    specification policy.
``legacy``
    the legacy Romanian TP rows the draft cites (top-1000 frame and
    full lexicon), kept apart from the harmonized ``r17`` rows.

``matched`` is True for the rows at the shared settings of the
comparison: E0, τ ∈ {dictionary, ge1}, dominant doublet coding, lemma
frequency basis, nouns, every root counted as its own root (the
default of all three languages), full lexicon (no V frame; learners:
V ∈ {1000, 5000}). Other rows are sensitivities, model fits off these
settings, or the legacy tables.

Columns (:data:`FIELDS`)
------------------------
Descriptors: ``section``, ``lang``, ``item`` (the row's label),
``target_class``, ``cell`` (the trigger, cell or pair as the language
names it), ``subset``, ``tau``, ``evidence_set``, ``coding``,
``freq_basis``, ``V``, ``rule`` (single-rule rows), ``matched``.
TP: ``N_target``, ``N`` (= N_trigger), ``n_shielded``,
``e_palatalize``, ``e_retain``, ``theta`` (N/ln N),
``tolerated_palatalize``, ``tolerated_retain``,
``fragility_palatalize``, ``fragility_retain``, ``theta_target``,
``counterfactual_palatalize`` (palatalize tolerated over N_target),
``verdict`` (two rules: palatalize / retain / both / neither /
empty; one rule: tolerated / not tolerated; nested set: conditioned /
not conditioned; learners: the modal outcome). Learners:
``share_palatalize``, ``share_retain``, ``share_neither``,
``share_conditioned``. Base class × derivative outcome (R21, R23):
``alt_pal``/``alt_plain`` (alternating, or base P), ``nonalt_pal``/
``nonalt_plain`` (non-alternating, or base N), ``noev_pal``/
``noev_plain`` (no evidence), ``n_both`` (units with both outcomes),
``roots_alternating``, ``roots_nonalternating``, ``roots_no_evidence``.
Statistics (R23): ``odds_ratio``, ``or_low``, ``or_high``,
``fisher_p``, ``mh_or``, ``mh_low``, ``mh_high``, ``mh_strata``,
``perm_p``, ``logit_base_class``, ``logit_se``, ``logit_p``,
``logit_status`` (estimated / separated / aliased / skipped / not
estimable / absent). Phonology: ``n_match``, ``n_ipa_match`` (of
``N`` rows). Provenance: ``note``, ``input`` (the tracked file the row
comes from), ``source``, ``licence_class`` (R9 gate).

Running
-------
``PYTHONPATH=src python -m lexpipe.compare`` (``--root DIR``,
``--force``; ``--check`` exits 1 when the tracked outputs are stale).
The step runs in the run context of the comparison-only registry row
``ro`` (:mod:`lexpipe.spec`), so it may write only under
``analysis/comparison/`` and ``data/romanian_compare/``; its manifest is
``data/romanian_compare/manifests/crosslang.manifest.tsv`` (R6) and it
is skipped while inputs, parameters (the version and the SHA-256 of
this module) and outputs match.
"""

from __future__ import annotations

import argparse
import csv
import io
import math
import re
import sys
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

from lexpipe import PROJECT_ROOT, csvio, stats, tp
from lexpipe.context import RunContext, run_stage, sha256_file
from lexpipe.spec import LANGUAGES

VERSION: Final[str] = "crosslang/1"
STAGE: Final[str] = "crosslang"
LANGS: Final[tuple[str, ...]] = ("ro", "it", "cs")
LANG_NAMES: Final[Mapping[str, str]] = {
    "ro": "Romanian",
    "it": "Italian",
    "cs": "Czech",
}
# The source of rows parsed from text reports (CSV rows carry their own).
DEFAULT_SOURCE: Final[Mapping[str, str]] = {
    "ro": "wiktionary_ro|wiktionary_en",
    "it": "wiktionary_en",
    "cs": "wiktionary_en|wiktionary_cs",
}
LICENCE_CLASS: Final[str] = "by-sa"
TAUS: Final[tuple[str, ...]] = ("dictionary", "ge1")
LEARNER_VS: Final[tuple[str, ...]] = ("1000", "5000")
SECTIONS: Final[tuple[str, ...]] = (
    "r17",
    "learners",
    "model_fit",
    "spellout",
    "noevidence",
    "covariation",
    "logtable",
    "phonology",
    "legacy",
)

FIELDS: Final[tuple[str, ...]] = (
    "section",
    "lang",
    "item",
    "target_class",
    "cell",
    "subset",
    "tau",
    "evidence_set",
    "coding",
    "freq_basis",
    "V",
    "rule",
    "matched",
    "N_target",
    "N",
    "n_shielded",
    "e_palatalize",
    "e_retain",
    "theta",
    "tolerated_palatalize",
    "tolerated_retain",
    "fragility_palatalize",
    "fragility_retain",
    "theta_target",
    "counterfactual_palatalize",
    "verdict",
    "share_palatalize",
    "share_retain",
    "share_neither",
    "share_conditioned",
    "alt_pal",
    "alt_plain",
    "nonalt_pal",
    "nonalt_plain",
    "noev_pal",
    "noev_plain",
    "n_both",
    "alt_avoid",
    "nonalt_avoid",
    "noev_avoid",
    "roots_alternating",
    "roots_nonalternating",
    "roots_no_evidence",
    "odds_ratio",
    "or_low",
    "or_high",
    "fisher_p",
    "mh_or",
    "mh_low",
    "mh_high",
    "mh_strata",
    "perm_p",
    "logit_base_class",
    "logit_se",
    "logit_p",
    "logit_status",
    "n_match",
    "n_ipa_match",
    "note",
    "input",
    "source",
    "licence_class",
)

Row = dict[str, Any]
CsvRow = Mapping[str, str]


class ReportFormatError(ValueError):
    """A tracked report does not have the expected format."""


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Inputs:
    """Every tracked file the summary reads, relative to a root."""

    root: Path

    def path(self, relative: str) -> Path:
        return self.root / relative

    @property
    def ro_tp(self) -> Path:
        return self.path("analysis/comparison/romanian_harmonized_tp_long.csv")

    @property
    def ro_learners(self) -> Path:
        return self.path("analysis/comparison/romanian_learners.csv")

    @property
    def ro_spellout(self) -> Path:
        return self.path("analysis/comparison/romanian_spellout.txt")

    @property
    def ro_noevidence(self) -> Path:
        return self.path("analysis/comparison/romanian_noevidence.txt")

    @property
    def ro_covariation(self) -> Path:
        return self.path("analysis/comparison/romanian_covariation.txt")

    @property
    def ro_legacy_top(self) -> Path:
        return self.path("analysis/romanian_tp_summary_downsampled.csv")

    @property
    def ro_legacy_full(self) -> Path:
        return self.path("analysis/romanian_tp_summary_full.csv")

    @property
    def ro_validation(self) -> Path:
        return self.path("analysis/palatalization_rule_validation.txt")

    @property
    def it_tp(self) -> Path:
        return self.path("analysis/italian/italian_tp_long.csv")

    @property
    def it_learners(self) -> Path:
        return self.path("analysis/italian/italian_learners.csv")

    @property
    def it_noevidence(self) -> Path:
        return self.path("analysis/italian/italian_noevidence.txt")

    @property
    def it_covariation(self) -> Path:
        return self.path("analysis/italian/italian_covariation.txt")

    def it_validation(self, variant: str) -> Path:
        suffix = "" if variant == "lexical" else f"_{variant}"
        return self.path(
            f"analysis/italian/italian_rule_validation{suffix}.txt"
        )

    @property
    def cs_tp(self) -> Path:
        return self.path("analysis/czech/czech_tp_long.csv")

    @property
    def cs_learners(self) -> Path:
        return self.path("analysis/czech/czech_learners.csv")

    @property
    def cs_spellout(self) -> Path:
        return self.path("analysis/czech/czech_spellout.txt")

    @property
    def cs_noevidence(self) -> Path:
        return self.path("analysis/czech/czech_noevidence.txt")

    @property
    def cs_covariation(self) -> Path:
        return self.path("analysis/czech/czech_covariation.txt")

    @property
    def cs_logtable(self) -> Path:
        return self.path("analysis/czech/czech_logtable.csv")

    @property
    def cs_logtable_summary(self) -> Path:
        return self.path("analysis/czech/czech_logtable_summary.txt")

    def cs_validation(self, variant: str) -> Path:
        suffix = "" if variant == "lexical" else f"_{variant}"
        return self.path(f"analysis/czech/czech_rule_validation{suffix}.txt")

    def all(self) -> list[Path]:
        """Every input, in a fixed order (the manifest's input list)."""
        fixed = [
            self.ro_tp,
            self.ro_learners,
            self.ro_spellout,
            self.ro_noevidence,
            self.ro_covariation,
            self.ro_legacy_top,
            self.ro_legacy_full,
            self.ro_validation,
            self.it_tp,
            self.it_learners,
            self.it_noevidence,
            self.it_covariation,
            self.cs_tp,
            self.cs_learners,
            self.cs_spellout,
            self.cs_noevidence,
            self.cs_covariation,
            self.cs_logtable,
            self.cs_logtable_summary,
        ]
        fixed += [self.it_validation(v) for v in IT_VARIANTS]
        fixed += [self.cs_validation(v) for v in CS_VARIANTS]
        return fixed

    def relative(self, path: Path) -> str:
        return path.relative_to(self.root).as_posix()


@dataclass(frozen=True, slots=True)
class Outputs:
    """The tracked outputs, under ``analysis/comparison/``."""

    csv: Path
    readme: Path


def outputs(root: Path) -> Outputs:
    out = root / "analysis" / "comparison"
    return Outputs(csv=out / "crosslang_summary.csv", readme=out / "README.md")


# Policies of the Logical Phonology validators (src/phonology/analyses).
IT_VARIANTS: Final[tuple[str, ...]] = (
    "lexical",
    "stress",
    "all_underspecified",
)
CS_VARIANTS: Final[tuple[str, ...]] = (
    "lexical",
    "all_underspecified",
    "lexical_loo",
    "spellout",
)


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


def _blank_row(section: str, lang: str, item: str) -> Row:
    row: Row = dict.fromkeys(FIELDS)
    row.update(section=section, lang=lang, item=item, matched=False)
    row["source"] = DEFAULT_SOURCE[lang]
    row["licence_class"] = LICENCE_CLASS
    return row


def _int(text: str) -> int:
    return int(text.replace(",", "").strip())


def _opt_int(text: str) -> int | None:
    text = text.strip()
    return None if text in ("", "-", "NA") else _int(text)


def _float(text: str) -> float:
    return float(text.strip())


def _opt_float(text: str) -> float | None:
    text = text.strip()
    if text in ("", "-", "NA", "n/a"):
        return None
    return float(text)


def _bool_cell(text: str) -> bool | None:
    """``True``/``False``/``yes``/``no``/``TRUE``; empty, ``-``, ``NA``
    and ``n/a`` are None."""
    lowered = text.strip().lower()
    if lowered in ("", "-", "na", "n/a"):
        return None
    if lowered in ("true", "yes"):
        return True
    if lowered in ("false", "no"):
        return False
    raise ValueError(f"not a boolean cell: {text!r}")


def two_rule_verdict(pal: bool | None, ret: bool | None) -> str:
    """The productive rule of a two-outcome cell (R18)."""
    if pal is None or ret is None:
        return ""
    return {
        (True, True): "both",
        (True, False): "palatalize",
        (False, True): "retain",
        (False, False): "neither",
    }[(pal, ret)]


def one_rule_verdict(tolerated: bool | None) -> str:
    if tolerated is None:
        return ""
    return "tolerated" if tolerated else "not tolerated"


def _theta_of(n: int) -> float | None:
    return tp.theta(n) if n > 1 else None


def _match(pattern: re.Pattern[str], line: str) -> re.Match[str] | None:
    return pattern.fullmatch(line.rstrip())


def _require(found: Any, what: str, label: str) -> Any:
    if found is None or found == [] or found == {}:
        raise ReportFormatError(f"{label}: {what} not found")
    return found


# ---------------------------------------------------------------------------
# TP long tables
# ---------------------------------------------------------------------------

DEFAULTS: Final[Mapping[str, str]] = {
    "evidence_set": "E0",
    "coding": "dominant",
    "freq_basis": "lemma",
    "V": "",
}


class TPTable:
    """A TP long table (:data:`lexpipe.tp.TP_LONG_TRACKED_FIELDS`)
    indexed for exact lookups; every lookup must find one row."""

    def __init__(self, rows: Sequence[CsvRow], label: str) -> None:
        self.label = label
        self._by: dict[tuple[str, str, str, str], list[CsvRow]] = {}
        for row in rows:
            key = (
                row["table"],
                row["target_class"],
                row["trigger"],
                row["tau"],
            )
            self._by.setdefault(key, []).append(row)

    @classmethod
    def read(cls, path: Path, label: str) -> TPTable:
        return cls(csvio.read_rows(path), label)

    def find_all(
        self, *, memo_contains: str | None = None, **where: str
    ) -> list[CsvRow]:
        crit = {**DEFAULTS, **where}
        key = (
            crit["table"],
            crit["target_class"],
            crit["trigger"],
            crit["tau"],
        )
        return [
            row
            for row in self._by.get(key, [])
            if all(row.get(k, "") == v for k, v in crit.items())
            and (memo_contains is None or memo_contains in row["memo"])
        ]

    def find(
        self, *, memo_contains: str | None = None, **where: str
    ) -> CsvRow:
        found = self.find_all(memo_contains=memo_contains, **where)
        if len(found) != 1:
            raise ReportFormatError(
                f"{self.label}: {len(found)} rows match "
                f"{dict(DEFAULTS, **where)}"
                + (f" with memo ~ {memo_contains!r}" if memo_contains else "")
            )
        return found[0]

    def pair(
        self, *, memo_contains: str | None = None, **where: str
    ) -> tuple[CsvRow, CsvRow]:
        """The palatalize and retain rows of one cell."""
        pal = self.find(
            rule="palatalize", memo_contains=memo_contains, **where
        )
        ret = self.find(rule="retain", memo_contains=memo_contains, **where)
        for name in ("N_target", "N_trigger", "theta_trigger"):
            if pal[name] != ret[name]:
                raise ReportFormatError(
                    f"{self.label}: palatalize and retain rows differ in "
                    f"{name} for {where}"
                )
        return pal, ret


def _descriptors(row: Row, src: CsvRow) -> None:
    row["target_class"] = src["target_class"]
    row["cell"] = src["trigger"]
    row["subset"] = src["subset"]
    row["tau"] = src["tau"]
    row["evidence_set"] = src["evidence_set"]
    row["coding"] = src["coding"]
    row["freq_basis"] = src["freq_basis"]
    row["V"] = src["V"]
    row["source"] = src["source"]
    row["licence_class"] = src["licence_class"]


def tp_pair_row(
    section: str,
    lang: str,
    item: str,
    pal: CsvRow,
    ret: CsvRow,
    *,
    input_path: str,
    matched: bool = False,
    note: str = "",
) -> Row:
    """One summary row from the palatalize and retain rows of a cell."""
    row = _blank_row(section, lang, item)
    _descriptors(row, pal)
    row.update(
        matched=matched,
        N_target=_int(pal["N_target"]),
        N=_int(pal["N_trigger"]),
        n_shielded=_int(pal["n_shielded"]),
        e_palatalize=_int(pal["e"]),
        e_retain=_int(ret["e"]),
        theta=_opt_float(pal["theta_trigger"]),
        tolerated_palatalize=_bool_cell(pal["tolerated_trigger"]),
        tolerated_retain=_bool_cell(ret["tolerated_trigger"]),
        fragility_palatalize=_opt_int(pal["fragility"]),
        fragility_retain=_opt_int(ret["fragility"]),
        theta_target=_opt_float(pal["theta_target"]),
        counterfactual_palatalize=_bool_cell(pal["tolerated_counterfactual"]),
        note=note,
        input=input_path,
    )
    row["verdict"] = two_rule_verdict(
        row["tolerated_palatalize"], row["tolerated_retain"]
    )
    return row


def tp_single_row(
    section: str,
    lang: str,
    item: str,
    src: CsvRow,
    *,
    input_path: str,
    matched: bool = False,
    note: str = "",
) -> Row:
    """One single-rule summary row (a nested branch, a rule TP row)."""
    rule = src["rule"]
    if rule not in tp.RULES:
        raise ReportFormatError(f"{input_path}: unknown rule {rule!r}")
    row = _blank_row(section, lang, item)
    _descriptors(row, src)
    tolerated = _bool_cell(src["tolerated_trigger"])
    row.update(
        rule=rule,
        matched=matched,
        N_target=_int(src["N_target"]),
        N=_int(src["N_trigger"]),
        n_shielded=_int(src["n_shielded"]),
        theta=_opt_float(src["theta_trigger"]),
        verdict=one_rule_verdict(tolerated),
        note=note,
        input=input_path,
    )
    row[f"e_{rule}"] = _int(src["e"])
    row[f"tolerated_{rule}"] = tolerated
    row[f"fragility_{rule}"] = _opt_int(src["fragility"])
    return row


def grammar_row(
    section: str,
    lang: str,
    item: str,
    branches: Sequence[Row],
    *,
    note: str = "",
) -> Row:
    """The closing row of a nested rule set: ``conditioned`` when every
    branch is tolerated (R18; a branch with N ≤ 1 is not)."""
    if not branches:
        raise ValueError(f"{item}: a nested set needs branches")
    first = branches[0]
    row = _blank_row(section, lang, item)
    for name in (
        "target_class",
        "cell",
        "tau",
        "evidence_set",
        "coding",
        "freq_basis",
        "V",
        "matched",
        "input",
        "source",
        "licence_class",
    ):
        row[name] = first[name]
    row["subset"] = "nested set"
    row["N"] = sum(int(b["N"]) for b in branches)
    ok = all(b["verdict"] == "tolerated" for b in branches)
    row["verdict"] = "conditioned" if ok else "not conditioned"
    row["note"] = note
    return row


# ---------------------------------------------------------------------------
# Section r17
# ---------------------------------------------------------------------------

R17A: Final[str] = "R17a dorsal/_[+front V]"
R17B: Final[str] = "R17b dorsal/_[+high,+front]"
R17B_NT: Final[str] = "R17b dorsal/_[+high,+front], non-trigger -e"
R17C: Final[str] = "R17c coronal/_[+high,+front]"
R17C_NT: Final[str] = "R17c coronal/_[+high,+front], non-trigger -e"


@dataclass(frozen=True, slots=True)
class R17Cell:
    item: str
    table: str
    target: str
    trigger: str
    note: str = ""


R17_CELLS: Final[Mapping[str, tuple[R17Cell, ...]]] = {
    "ro": (
        R17Cell(R17A, "R17a", "dorsal", "-i/-e"),
        R17Cell(R17B, "R17b", "dorsal", "-i"),
        R17Cell(R17B_NT, "R17b", "dorsal", "-e"),
        R17Cell(R17C, "R17c", "coronal", "-i", "coronal = t, d, s"),
        R17Cell(R17C_NT, "R17c", "coronal", "-e", "coronal = t, d, s"),
        R17Cell(R17C, "R17c", "z", "-i"),
        R17Cell(R17C_NT, "R17c", "z", "-e"),
    ),
    "it": (
        R17Cell(
            R17A,
            "R17a",
            "dorsal",
            "-i/-e",
            "a root whose -i and -e cells disagree (amici, amiche) is an "
            "exception to both rules (per-cell doublet coding; Czech and "
            "Romanian code such roots by token majority)",
        ),
        R17Cell(R17B, "R17b", "dorsal", "-i"),
        R17Cell(R17B_NT, "R17b", "dorsal", "-e"),
        R17Cell(
            R17C,
            "R17c",
            "coronal",
            "-i",
            "every aligned coronal plural is coded unchanged: N = e by "
            "construction",
        ),
        R17Cell(
            R17C_NT,
            "R17c",
            "coronal",
            "-e",
            "every aligned coronal plural is coded unchanged: N = e by "
            "construction",
        ),
    ),
    "cs": (
        R17Cell(
            R17A,
            "R17a",
            "dorsal",
            "-e/-i",
            "front vowels of every ending, triggers (-e, -i) and "
            "non-triggers (-em, -ech, -é) alike",
        ),
        R17Cell(R17B, "R17b", "dorsal", "-i"),
        R17Cell(
            R17B_NT,
            "R17b",
            "dorsal",
            "-e",
            "Czech -e both triggers (LOC/DAT -e) and does not (-em, -ech)",
        ),
        R17Cell(
            R17C,
            "R17c",
            "coronal_inferred",
            "-i",
            "t, d, n: the change is inferred from spelling, e = 0 by "
            "construction",
        ),
        R17Cell(
            R17C_NT,
            "R17c",
            "coronal_inferred",
            "-e",
            "t, d, n: the change is inferred from spelling",
        ),
        R17Cell(R17C, "R17c", "r", "-i", "r → ř"),
        R17Cell(R17C_NT, "R17c", "r", "-e", "r → ř"),
        R17Cell(R17C, "R17c", "sibilant", "-i", "s, z"),
        R17Cell(R17C_NT, "R17c", "sibilant", "-e", "s, z"),
        R17Cell(R17C, "R17c", "lateral", "-i", "l"),
        R17Cell(R17C_NT, "R17c", "lateral", "-e", "l"),
    ),
}
# The nouns-only subset of each language (matched) and its sensitivity
# subsets (Italian 'all' adds the adjective roots, most of them -ico).
R17_SUBSETS: Final[Mapping[str, tuple[tuple[str, bool], ...]]] = {
    "ro": (("all", True), ("core", False)),
    "it": (("nouns", True), ("all", False), ("core", False)),
    "cs": (("all", True), ("core", False)),
}
R17_SUBSET_NOTES: Final[Mapping[tuple[str, str], str]] = {
    ("it", "all"): "nouns and adjectives (the Italian headline subset)",
    ("ro", "core"): "without the outcome-blind flags (suffix targets, "
    "palatal context, suppletive stems, -ică)",
    ("it", "core"): "nouns and adjectives without the -ico/-iaco suffix "
    "targets and the -logo/-fago/-urgo suffixoids",
    ("cs", "core"): "without the R26 suffix targets",
}


def r17_rows(tables: Mapping[str, TPTable], inputs: Inputs) -> list[Row]:
    out: list[Row] = []
    paths = {"ro": inputs.ro_tp, "it": inputs.it_tp, "cs": inputs.cs_tp}
    for lang in LANGS:
        table = tables[lang]
        label = inputs.relative(paths[lang])
        for cell in R17_CELLS[lang]:
            for subset, matched in R17_SUBSETS[lang]:
                for V in ("", "1000"):
                    for tau in TAUS:
                        where = {
                            "table": cell.table,
                            "target_class": cell.target,
                            "trigger": cell.trigger,
                            "subset": subset,
                            "tau": tau,
                            "V": V,
                        }
                        optional = bool(V) or not matched
                        if optional and not table.find_all(
                            rule="palatalize", **where
                        ):
                            continue  # no V sweep or no such sensitivity
                        pal, ret = table.pair(**where)
                        notes = [cell.note]
                        notes.append(R17_SUBSET_NOTES.get((lang, subset), ""))
                        if V:
                            notes.append(
                                f"R24 frame: the {V} most frequent noun roots"
                            )
                        out.append(
                            tp_pair_row(
                                "r17",
                                lang,
                                cell.item,
                                pal,
                                ret,
                                input_path=label,
                                matched=matched and not V,
                                note="; ".join(n for n in notes if n),
                            )
                        )
    return out


# ---------------------------------------------------------------------------
# Section learners
# ---------------------------------------------------------------------------

# (table, target, trigger, subset) of the R17 learner cells: matched.
R17_LEARNER_CELLS: Final[
    Mapping[str, frozenset[tuple[str, str, str, str]]]
] = {
    "ro": frozenset(
        {
            ("R17a", "dorsal", "-i/-e", "with_nde"),
            ("R17b", "dorsal", "-i", "all"),
            ("R17b", "dorsal", "-e", "all"),
            ("R17c", "coronal", "-i", "all"),
            ("R17c", "coronal", "-e", "all"),
            ("R17c", "z", "-i", "all"),
            ("R17c", "z", "-e", "all"),
        }
    ),
    "it": frozenset(
        {
            ("R17a", "dorsal", "-i/-e", "with_nde"),
            ("R17b", "dorsal", "-i", "all"),
            ("R17b", "dorsal", "-e", "all"),
            ("R17c", "coronal", "-i", "all"),
        }
    ),
    "cs": frozenset(
        {
            ("R17a", "dorsal", "-e/-i", "all"),
            ("R17b", "dorsal", "-i", "all"),
            ("R17b", "dorsal", "-e", "all"),
            ("R17c", "coronal_inferred", "-i", "all"),
            ("R17c", "r", "-i", "all"),
            ("R17c", "r", "-e", "all"),
            ("R17c", "sibilant", "-i", "all"),
            ("R17c", "lateral", "-i", "all"),
        }
    ),
}
LEARNER_OUTCOMES: Final[tuple[str, ...]] = (
    "palatalize",
    "retain",
    "conditioned",
    "neither",
)


def _learner_item(row: CsvRow) -> str:
    item = f"{row['table']} {row['target_class']} {row['trigger']}"
    if row["subset"] not in ("all",):
        item += f" [{row['subset']}]"
    if row["nested"]:
        item += f" (nested: {row['nested']})"
    return item


def learner_rows(
    learner_tables: Mapping[str, Sequence[CsvRow]], inputs: Inputs
) -> list[Row]:
    out: list[Row] = []
    paths = {
        "ro": inputs.ro_learners,
        "it": inputs.it_learners,
        "cs": inputs.cs_learners,
    }
    for lang in LANGS:
        label = inputs.relative(paths[lang])
        seen: set[tuple[str, ...]] = set()
        for src in learner_tables[lang]:
            if (
                src["tau"] not in TAUS
                or src["V"] not in LEARNER_VS
                or src["coding"] != "dominant"
                or src["freq_basis"] != "lemma"
                or src["evidence_set"] not in ("E0", "E1")
            ):
                continue
            key = (
                src["table"],
                src["target_class"],
                src["trigger"],
                src["subset"],
                src["nested"],
                src["tau"],
                src["evidence_set"],
                src["V"],
            )
            if key in seen:
                raise ReportFormatError(
                    f"{label}: duplicate learner row {key}"
                )
            seen.add(key)
            cell = (
                src["table"],
                src["target_class"],
                src["trigger"],
                src["subset"],
            )
            row = _blank_row("learners", lang, _learner_item(src))
            shares = {
                name: _float(src[f"share_{name}"]) for name in LEARNER_OUTCOMES
            }
            row.update(
                target_class=src["target_class"],
                cell=src["trigger"],
                subset=src["subset"],
                tau=src["tau"],
                evidence_set=src["evidence_set"],
                coding=src["coding"],
                freq_basis=src["freq_basis"],
                V=src["V"],
                matched=(
                    cell in R17_LEARNER_CELLS[lang]
                    and src["evidence_set"] == "E0"
                ),
                verdict=max(
                    LEARNER_OUTCOMES,
                    key=lambda name: (
                        shares[name],
                        -LEARNER_OUTCOMES.index(name),
                    ),
                ),
                note=(
                    f"{src['n_samples']} vocabularies, seed {src['seed']}"
                    + (
                        f"; nested rule set {src['nested']}"
                        if src["nested"]
                        else ""
                    )
                ),
                input=label,
                source=src["source"],
                licence_class=src["licence_class"],
            )
            for name, value in shares.items():
                row[f"share_{name}"] = value
            out.append(row)
    return out


# ---------------------------------------------------------------------------
# Section model_fit
# ---------------------------------------------------------------------------


def _nested_set(
    lang: str,
    table: TPTable,
    label: str,
    title: str,
    branches: Sequence[tuple[str, str, str]],
    where: Mapping[str, str],
    *,
    memo_contains: str | None,
    matched: bool,
    note: str = "",
) -> list[Row]:
    """Branch rows (subset, rule, label) of one nested set + its grammar."""
    rows: list[Row] = []
    for subset, rule, branch in branches:
        src = table.find(
            rule=rule, subset=subset, memo_contains=memo_contains, **where
        )
        rows.append(
            tp_single_row(
                "model_fit",
                lang,
                f"{title}: {branch}",
                src,
                input_path=label,
                matched=matched,
                note=note,
            )
        )
    rows.append(
        grammar_row("model_fit", lang, f"{title}: grammar", rows, note=note)
    )
    return rows


def ro_model_rows(table: TPTable, label: str) -> list[Row]:
    out: list[Row] = []
    for tau in TAUS:
        src = table.find(
            table="R17a",
            rule="palatalize",
            target_class="dorsal",
            trigger="-i/-e",
            subset="all",
            tau=tau,
        )
        out.append(
            tp_single_row(
                "model_fit",
                "ro",
                "dorsal rule: dorsal/_[+front V] → palatalize",
                src,
                input_path=label,
                matched=True,
            )
        )
        pal, ret = table.pair(
            table="R17c",
            target_class="coronal",
            trigger="-i/-e",
            subset="all",
            tau=tau,
            memo_contains="(pooled; not in R17)",
        )
        out.append(
            tp_pair_row(
                "model_fit",
                "ro",
                "single coronal rule over -i and -e (t, d, s)",
                pal,
                ret,
                input_path=label,
                matched=True,
                note="the pooled coronal/_[+front V] rule the draft rejects",
            )
        )
        out += _nested_set(
            "ro",
            table,
            label,
            "coronal rule (t, d, s)",
            (
                ("nested:e_plural", "retain", "-e plural → retain"),
                (
                    "nested:elsewhere",
                    "palatalize",
                    "elsewhere (-i) → palatalize",
                ),
            ),
            {
                "table": "R17c",
                "target_class": "coronal",
                "trigger": "-i/-e",
                "tau": tau,
            },
            memo_contains="(pooled; not in R17); nested in all",
            matched=True,
            note="the coronal rule is restricted to -i: nested in Yang's order",
        )
        src = table.find(
            table="R17c",
            rule="retain",
            target_class="z",
            trigger="-i",
            subset="all",
            tau=tau,
        )
        out.append(
            tp_single_row(
                "model_fit",
                "ro",
                "z before -i → retain (prespecified /z/)",
                src,
                input_path=label,
                matched=True,
            )
        )
    return out


IT_STRESS_SETS: Final[tuple[tuple[str, bool, str], ...]] = (
    (
        "giavazzi",
        True,
        "every -ico adjective a root of its own (as Romanian and Czech "
        "count suffix targets)",
    ),
    (
        "suffix_once",
        False,
        "R26: each suffix (-ico, -iaco) counted once, as one unit",
    ),
    (
        "orth_ico_once",
        False,
        "every root in unstressed -ico/-ica as one unit (upper bound)",
    ),
)


def it_model_rows(table: TPTable, label: str) -> list[Row]:
    out: list[Row] = []
    for tau in (*TAUS, "dictionary_explicit"):
        for name, matched, note in IT_STRESS_SETS:
            explicit = tau == "dictionary_explicit"
            out += _nested_set(
                "it",
                table,
                label,
                f"stress-nested -i ({name})",
                (
                    (
                        "nested:post_tonic",
                        "retain",
                        "velar post-tonic → retain",
                    ),
                    (
                        "nested:elsewhere",
                        "palatalize",
                        "elsewhere → palatalize",
                    ),
                ),
                {
                    "table": "it_stress",
                    "target_class": "dorsal",
                    "trigger": "-i",
                    "tau": tau,
                },
                memo_contains=f"nested set {name}",
                matched=matched and not explicit,
                note=note
                + (
                    "; typed plurals only (no Module:it-common defaults)"
                    if explicit
                    else ""
                ),
            )
    for tau in TAUS:
        pal, ret = table.pair(
            table="it_nounverb",
            target_class="dorsal",
            trigger="-i",
            subset="pooled",
            tau=tau,
            evidence_set="E1",
        )
        out.append(
            tp_pair_row(
                "model_fit",
                "it",
                "pooled noun + verb units, -i",
                pal,
                ret,
                input_path=label,
                note="E1: nouns with their -are/-ere/-ire verb cells",
            )
        )
        out += _nested_set(
            "it",
            table,
            label,
            "category_stress (noun + verb)",
            (
                ("nested:verb_are", "retain", "-are verbs → retain"),
                ("nested:verb_isc", "retain", "-isc- verbs → retain"),
                (
                    "nested:noun_post_tonic",
                    "retain",
                    "post-tonic nouns → retain",
                ),
                ("nested:elsewhere", "palatalize", "elsewhere → palatalize"),
            ),
            {
                "table": "it_nounverb",
                "target_class": "dorsal",
                "trigger": "-i",
                "tau": tau,
                "evidence_set": "E1",
            },
            memo_contains="nested set category_stress",
            matched=False,
            note="E1 verb units",
        )
        for subset, item in (
            ("conj:are", "verbs -are (2sg, 3sg sbjv, 1pl)"),
            ("conj:ere_root", "verbs root-stressed -ere"),
            ("conj:isc", "verbs -isc-"),
            ("all", "all velar-stem verbs"),
        ):
            pal, ret = table.pair(
                table="it_verbs",
                target_class="dorsal",
                trigger="-i",
                subset=subset,
                tau=tau,
                evidence_set="E1",
            )
            out.append(
                tp_pair_row(
                    "model_fit",
                    "it",
                    item,
                    pal,
                    ret,
                    input_path=label,
                    note="verb units, -i cells",
                )
            )
    return out


CS_RULE_CELLS: Final[tuple[tuple[str, str, str], ...]] = (
    ("rule", "velar", "NOM.PL m.an -i"),
    ("rule", "velar", "DAT/LOC.SG f -e"),
    ("rule", "velar", "LOC.PL m -ích"),
    ("rule_adj", "velar", "NOM.PL.M.AN -í"),
    ("rule", "r", "NOM.PL m.an -i"),
    ("rule", "r", "DAT/LOC.SG f -e"),
)
CS_COMPETITIONS: Final[tuple[tuple[str, str], ...]] = (
    ("NOM.PL m.an", "-i vs -ové/-é"),
    ("LOC.SG m.in", "-e vs -u"),
    ("LOC.SG n", "-e vs -u"),
    ("LOC.PL m.in", "-ích vs -ách/-ech (baseline: -ech vs -ách)"),
    ("VOC.SG m.an", "-e vs -u"),
)
CS_NOMPL_NESTED: Final[tuple[tuple[str, str], ...]] = (
    ("nested:suffixoid", "-log/-urg/-graf"),
    ("nested:sky", "-ský/-cký"),
    ("nested:suffix", "other R26 suffixes (-ák, -ík, -ník, -ek)"),
    ("nested:monosyllable", "monosyllables"),
    ("nested:loan", "loans"),
    ("nested:elsewhere", "elsewhere"),
)
CS_R22_SUFFIXES: Final[tuple[str, ...]] = (
    "-it",
    "-et",
    "-ní/-ný",
    "-ský/-cký",
    "-ina",
    "-í",
    "-ek",
    "-ka",
    "-ík",
    "-an",
    "-ník",
    "-nice",
    "-ista",
    "-ismus",
)


def cs_model_rows(table: TPTable, label: str) -> list[Row]:
    out: list[Row] = []
    for tau in TAUS:
        for name, target, trigger in CS_RULE_CELLS:
            src = table.find(
                table=name,
                rule="palatalize",
                target_class=target,
                trigger=trigger,
                subset="all",
                tau=tau,
            )
            out.append(
                tp_single_row(
                    "model_fit",
                    "cs",
                    f"rule TP (a): {target} before {trigger}",
                    src,
                    input_path=label,
                    matched=True,
                    note="suffix-indexed trigger; a competing N allomorph "
                    "shields the root",
                )
            )
        for cell, pair in CS_COMPETITIONS:
            for target in ("velar", "baseline"):
                pal, ret = table.pair(
                    table="competition",
                    target_class=target,
                    trigger=cell,
                    subset="all",
                    tau=tau,
                )
                out.append(
                    tp_pair_row(
                        "model_fit",
                        "cs",
                        f"competition {cell}: {pair}",
                        pal,
                        ret,
                        input_path=label,
                        matched=True,
                        note="palatalize = the P-default rule (e = roots "
                        "choosing the competitor), retain = the N-default "
                        "rule"
                        + (
                            "; baseline = hard labial, l, s, z stems (no "
                            "alternation involved)"
                            if target == "baseline"
                            else ""
                        ),
                    )
                )
        branches: list[Row] = []
        for subset, branch in CS_NOMPL_NESTED:
            found = table.find_all(
                table="competition",
                target_class="velar",
                trigger="NOM.PL m.an",
                subset=subset,
                tau=tau,
            )
            if len(found) != 1:
                raise ReportFormatError(
                    f"{label}: {len(found)} nested rows {subset} at {tau}"
                )
            src = found[0]
            if _int(src["N_trigger"]) == 0:
                continue  # an empty class is not part of the grammar
            branches.append(
                tp_single_row(
                    "model_fit",
                    "cs",
                    f"NOM.PL m.an nested: {branch} → {src['rule']}",
                    src,
                    input_path=label,
                    matched=True,
                )
            )
        out += branches
        out.append(
            grammar_row(
                "model_fit",
                "cs",
                "NOM.PL m.an nested: grammar",
                branches,
                note="Yang's order: -log/-urg/-graf, -ský/-cký, other "
                "suffixes, monosyllables, loans, elsewhere",
            )
        )
        for suffix in CS_R22_SUFFIXES:
            src = table.find(
                table="derivation",
                rule="palatalize",
                target_class="velar",
                trigger=suffix,
                subset="all",
                tau=tau,
                evidence_set="E2",
            )
            nontrigger = suffix in ("-ista", "-ismus")
            out.append(
                tp_single_row(
                    "model_fit",
                    "cs",
                    f"R22 derivation: velar base + {suffix}",
                    src,
                    input_path=label,
                    note=(
                        "non-trigger suffix: N = e by design"
                        if nontrigger
                        else "basis wikt (Wiktionary links and headword "
                        "candidates)"
                    ),
                )
            )
    return out


# ---------------------------------------------------------------------------
# Section spellout
# ---------------------------------------------------------------------------

RO_CONTRAST_RE: Final[re.Pattern[str]] = re.compile(
    r"\s+(dorsal|coronal|z|all)\s+non-alternating-plural roots "
    r"\(-uri or resisting\): -i only (\d+), -a/-ui only (\d+), both (\d+)"
)


def parse_ro_spellout_contrast(
    text: str, label: str
) -> dict[str, tuple[int, int, int]]:
    """``{class: (-i only, -a/-ui only, both)}`` of the contrast set
    (non-alternating-plural roots, τ = dictionary)."""
    found: dict[str, tuple[int, int, int]] = {}
    for line in text.splitlines():
        m = _match(RO_CONTRAST_RE, line)
        if m:
            found[m.group(1)] = (
                _int(m.group(2)),
                _int(m.group(3)),
                _int(m.group(4)),
            )
    return dict(_require(found, "contrast lines", label))


def ro_spellout_rows(
    table: TPTable, contrast_text: str, inputs: Inputs
) -> list[Row]:
    out: list[Row] = []
    label = inputs.relative(inputs.ro_tp)
    for target in ("dorsal", "coronal"):
        for tau in TAUS:
            pal, ret = table.pair(
                table="spellout",
                target_class=target,
                trigger="verbalizer",
                subset="all",
                tau=tau,
            )
            out.append(
                tp_pair_row(
                    "spellout",
                    "ro",
                    "-i vs -a/-ui verbalizer, alternating-plural roots",
                    pal,
                    ret,
                    input_path=label,
                    matched=True,
                    note="e (palatalize) = roots taking only -a/-ui",
                )
            )
    contrast_label = inputs.relative(inputs.ro_spellout)
    contrast = parse_ro_spellout_contrast(contrast_text, contrast_label)
    for target in ("dorsal", "coronal"):
        if target not in contrast:
            raise ReportFormatError(f"{contrast_label}: no {target} contrast")
        only_i, only_a, both = contrast[target]
        row = _blank_row(
            "spellout", "ro", "contrast: non-alternating-plural roots"
        )
        row.update(
            target_class=target,
            cell="verbalizer",
            subset="contrast:nonalt",
            tau="dictionary",
            evidence_set="E0",
            coding="raw",
            matched=True,
            N=only_i + only_a + both,
            e_palatalize=only_a,
            e_retain=only_i,
            n_both=both,
            note="-uri or resisting plurals; e (palatalize) = -a/-ui only, "
            "e (retain) = -i only; both counted apart",
            input=contrast_label,
        )
        out.append(row)
    return out


IT_SPELLOUT_SETS: Final[tuple[tuple[str, str], ...]] = (
    ("verbalizer", "-izzare vs -eggiare/-ire/-are"),
    ("verbalizer_no_are", "-izzare vs -eggiare/-ire"),
    ("ista_ismo_ita", "palatal vs velar -ista/-ismo/-ità"),
)


def it_spellout_rows(table: TPTable, label: str) -> list[Row]:
    out: list[Row] = []
    for trigger, item in IT_SPELLOUT_SETS:
        for subset, prefix in (("all", ""), ("contrast:nonalt", "contrast: ")):
            for tau in TAUS:
                pal, ret = table.pair(
                    table="spellout",
                    target_class="dorsal",
                    trigger=trigger,
                    subset=subset,
                    tau=tau,
                )
                out.append(
                    tp_pair_row(
                        "spellout",
                        "it",
                        prefix + item,
                        pal,
                        ret,
                        input_path=label,
                        matched=True,
                        note=(
                            "roots with an attested -chi/-ghi plural"
                            if subset != "all"
                            else "roots with an attested -ci/-gi plural"
                        ),
                    )
                )
    return out


CS_SPELLOUT_ROW_RE: Final[re.Pattern[str]] = re.compile(
    r"  (?P<pair>\S.*?)\s+(?P<tau>dictionary|ge1|ge0\.5pm|ge2pm)\s+"
    r"(?P<coding>dominant|doublet_exception|doublet_excluded)\s+"
    r"(?P<subset>\S.*?)\s+(?P<roots>\d+)\s+(?P<N>\d+)\s+(?P<Nonly>\d+)\s+"
    r"(?P<e>\d+)\s+(?P<theta>\S+)\s+(?P<tol>\S+)\s+(?P<frag>\S+)\s+"
    r"(?P<eret>\d+)"
)
CS_SPELLOUT_SECTION_RE: Final[re.Pattern[str]] = re.compile(r"(\d+b?)\. (.*)")


@dataclass(frozen=True, slots=True)
class CsSpellout:
    section: str
    pair: str
    tau: str
    coding: str
    subset: str
    roots: int
    N: int
    n_only: int
    e: int
    theta: float | None
    tolerated: bool | None
    fragility: int | None
    e_ret: int


def parse_cs_spellout(text: str, label: str) -> list[CsSpellout]:
    """The Wiktionary-basis tables of ``czech_spellout.txt`` (sections
    1, 1b, 4, 4b); parsing stops at the NC section."""
    out: list[CsSpellout] = []
    section = ""
    for line in text.splitlines():
        m = _match(CS_SPELLOUT_SECTION_RE, line)
        if m:
            section = m.group(1)
            if m.group(2).startswith("NC sources"):
                break
            continue
        if section not in ("1", "1b", "4", "4b"):
            continue
        row = _match(CS_SPELLOUT_ROW_RE, line)
        if row is None:
            continue
        out.append(
            CsSpellout(
                section=section,
                pair=row["pair"].strip(),
                tau=row["tau"],
                coding=row["coding"],
                subset=row["subset"].strip(),
                roots=_int(row["roots"]),
                N=_int(row["N"]),
                n_only=_int(row["Nonly"]),
                e=_int(row["e"]),
                theta=_opt_float(row["theta"]),
                tolerated=_bool_cell(row["tol"]),
                fragility=_opt_int(row["frag"]),
                e_ret=_int(row["eret"]),
            )
        )
    return list(_require(out, "spell-out tables", label))


CS_SPELLOUT_ITEMS: Final[tuple[tuple[str, str, str, str], ...]] = (
    # (section, pair, subset, root set)
    ("1", "-ní/-ný vs -ový/-ovní", "all", "feminine"),
    ("1b", "-ní/-ný vs -ový/-ovní", "base:root-final", "feminine"),
    ("1b", "-ní/-ný vs -ový/-ovní", "base:suffix", "feminine"),
    ("1", "-it vs -ovat", "all", "feminine"),
    ("1b", "-it vs -ovat", "base:root-final", "feminine"),
    ("1b", "-it vs -ovat", "base:suffix", "feminine"),
    ("1", "-ka vs -yně", "all", "feminine"),
    ("4", "-ka vs -yně", "all", "E0-certain"),
    ("1", "any pair", "all", "feminine"),
)


def cs_spellout_rows(text: str, label: str) -> list[Row]:
    parsed = parse_cs_spellout(text, label)
    index = {(p.section, p.pair, p.subset, p.tau, p.coding): p for p in parsed}
    out: list[Row] = []
    for section, pair, subset, root_set in CS_SPELLOUT_ITEMS:
        for tau in TAUS:
            p = index.get((section, pair, subset, tau, "dominant"))
            if p is None:
                raise ReportFormatError(
                    f"{label}: no row {section} {pair} {subset} {tau}"
                )
            row = _blank_row("spellout", "cs", pair)
            ret_tol = tp.tolerated(p.N, p.e_ret) if p.N else None
            row.update(
                target_class="velar",
                cell=pair,
                subset=f"{root_set} roots; {subset}",
                tau=tau,
                evidence_set="E0",
                coding="dominant",
                matched=root_set == "feminine" or pair == "-ka vs -yně",
                N_target=p.roots,
                N=p.N,
                e_palatalize=p.e,
                e_retain=p.e_ret,
                theta=p.theta,
                tolerated_palatalize=p.tolerated if p.N > 1 else None,
                tolerated_retain=ret_tol,
                fragility_palatalize=p.fragility,
                fragility_retain=(
                    tp.fragility(p.N, p.e_ret) if p.N > 1 else None
                ),
                note=(
                    f"certain roots: {root_set} set; basis wikt; "
                    "tolerated_retain computed here from N and e(ret)"
                ),
                input=label,
            )
            row["verdict"] = two_rule_verdict(
                row["tolerated_palatalize"], row["tolerated_retain"]
            )
            out.append(row)
    return out


# ---------------------------------------------------------------------------
# Section noevidence
# ---------------------------------------------------------------------------


@dataclass(slots=True)
class ClassTable:
    """Base class × derivative outcome counts of one R21 table."""

    alt: tuple[int, int, int] | None = None  # (pal, plain, both)
    nonalt: tuple[int, int, int] | None = None
    noev: tuple[int, int, int] | None = None

    def complete(self) -> bool:
        return None not in (self.alt, self.nonalt, self.noev)

    def fill(self, row: Row) -> None:
        if not self.complete():
            raise ValueError("incomplete base-class table")
        assert self.alt and self.nonalt and self.noev
        row["alt_pal"], row["alt_plain"] = self.alt[0], self.alt[1]
        row["nonalt_pal"], row["nonalt_plain"] = self.nonalt[0], self.nonalt[1]
        row["noev_pal"], row["noev_plain"] = self.noev[0], self.noev[1]
        row["n_both"] = self.alt[2] + self.nonalt[2] + self.noev[2]


RO_NE_TAU_RE: Final[re.Pattern[str]] = re.compile(r"tau = (\S+)")
RO_NE_ROOTS_RE: Final[re.Pattern[str]] = re.compile(
    r"\s+target roots by base class: alternating (\d+), non-alternating "
    r"(\d+), non-alternating -uri (\d+), no evidence (\d+)"
)
RO_NE_CLASS_RE: Final[re.Pattern[str]] = re.compile(
    r"\s+class (\S+) \(base x suffix units\)"
)
RO_NE_ROW_RE: Final[re.Pattern[str]] = re.compile(
    r"\s+(alternating|non-alternating -uri|non-alternating|no evidence|"
    r"\(non-alt\. pooled\))\s+(\d+)/(\d+)\s+(\d+)/(\d+)\s+(\d+)/(\d+)\s+"
    r"(\d+)/(\d+)\s+(\d+)\s+(\d+)"
)


@dataclass(frozen=True, slots=True)
class RoNoEvidence:
    tables: Mapping[tuple[str, str], Mapping[str, tuple[int, int, int, int]]]
    roots: Mapping[str, tuple[int, int, int, int]]


def parse_ro_noevidence(text: str, label: str) -> RoNoEvidence:
    """``romanian_noevidence.txt``: per (τ, class) and base class the
    'all front' units (pal, nonpal, both, avoid); root counts per τ."""
    tables: dict[tuple[str, str], dict[str, tuple[int, int, int, int]]] = {}
    roots: dict[str, tuple[int, int, int, int]] = {}
    tau = ""
    klass = ""
    for line in text.splitlines():
        if m := _match(RO_NE_TAU_RE, line):
            tau, klass = m.group(1), ""
            continue
        if line.strip().startswith("By base stratum"):
            klass = ""
            continue
        if m := _match(RO_NE_ROOTS_RE, line):
            roots[tau] = tuple(_int(m.group(i)) for i in range(1, 5))  # type: ignore[assignment]
            continue
        if m := _match(RO_NE_CLASS_RE, line):
            klass = m.group(1)
            continue
        if klass and (m := _match(RO_NE_ROW_RE, line)):
            tables.setdefault((tau, klass), {})[m.group(1)] = (
                _int(m.group(8)),
                _int(m.group(9)),
                _int(m.group(10)),
                _int(m.group(11)),
            )
    _require(tables, "base-class tables", label)
    return RoNoEvidence(tables=tables, roots=roots)


def ro_noevidence_rows(text: str, label: str) -> list[Row]:
    parsed = parse_ro_noevidence(text, label)
    out: list[Row] = []
    for tau in TAUS:
        roots = parsed.roots.get(tau)
        if roots is None:
            raise ReportFormatError(f"{label}: no root counts at {tau}")
        row = _blank_row("noevidence", "ro", "target roots by base class")
        row.update(
            target_class="all",
            tau=tau,
            evidence_set="E0",
            matched=True,
            roots_alternating=roots[0],
            roots_nonalternating=roots[1] + roots[2],
            roots_no_evidence=roots[3],
            note="all target classes (dorsal, coronal, z); non-alternating "
            "includes the -uri plurals",
            input=label,
        )
        out.append(row)
        for klass in ("dorsal", "coronal"):
            found = parsed.tables.get((tau, klass))
            if found is None:
                raise ReportFormatError(f"{label}: no {klass} table at {tau}")
            needed = ("alternating", "(non-alt. pooled)", "no evidence")
            if any(name not in found for name in needed):
                raise ReportFormatError(f"{label}: incomplete {klass} {tau}")
            t = ClassTable(
                alt=found["alternating"][:3],
                nonalt=found["(non-alt. pooled)"][:3],
                noev=found["no evidence"][:3],
            )
            row = _blank_row(
                "noevidence",
                "ro",
                "front-vowel derivatives (verbs -i, "
                "-iza/-ifica, adjectives)",
            )
            row.update(
                target_class=klass,
                cell="front suffixes",
                subset="all",
                tau=tau,
                evidence_set="E0",
                matched=True,
                alt_avoid=found["alternating"][3],
                nonalt_avoid=found["(non-alt. pooled)"][3],
                noev_avoid=found["no evidence"][3],
                note="non-alternating = plain plural before a trigger or "
                "-uri; *_avoid = units whose derivative takes a back-vowel "
                "suffix (-a/-ui verbs, back-vowel adjective suffixes)",
                input=label,
            )
            t.fill(row)
            out.append(row)
    return out


IT_NE_ROOTS_RE: Final[re.Pattern[str]] = re.compile(
    r"\s+(dictionary|ge1|ge0\.5pm|ge2pm)\s+alternating (\d+), "
    r"non-alternating (\d+), no evidence (\d+)"
)
IT_NE_BLOCK_RE: Final[re.Pattern[str]] = re.compile(
    r"tau = (\S+); links: (\S+)"
)
IT_NE_SET_RE: Final[re.Pattern[str]] = re.compile(
    r"  (learned set|native control):"
)
IT_NE_TABLE_RE: Final[re.Pattern[str]] = re.compile(r"    (\S.*)")
IT_NE_ROW_RE: Final[re.Pattern[str]] = re.compile(
    r"      (alternating|non-alternating|no evidence)\s+(\d+)\s+(\d+)\s+(\d+)"
)
IT_NE_DESC_HEAD: Final[str] = (
    "no-evidence units by base descriptor (P/N/both):"
)
IT_NE_DESC_RE: Final[re.Pattern[str]] = re.compile(
    r"      (\S.*?)\s+(\d+)/\s*(\d+)/\s*(\d+)"
)


@dataclass(frozen=True, slots=True)
class ItNoEvidence:
    tables: Mapping[tuple[str, str, str, str], ClassTable]
    descriptors: Mapping[
        tuple[str, str, str], Mapping[str, tuple[int, int, int]]
    ]
    roots: Mapping[str, tuple[int, int, int]]


def parse_it_noevidence(text: str, label: str) -> ItNoEvidence:
    """``italian_noevidence.txt``: (τ, links, set, table) → class table;
    the no-evidence descriptors; dorsal roots per base class and τ."""
    tables: dict[tuple[str, str, str, str], ClassTable] = {}
    descriptors: dict[
        tuple[str, str, str], dict[str, tuple[int, int, int]]
    ] = {}
    roots: dict[str, tuple[int, int, int]] = {}
    tau = links = group = table = ""
    in_desc = False
    for line in text.splitlines():
        if m := _match(IT_NE_ROOTS_RE, line):
            if not tau:
                roots[m.group(1)] = (
                    _int(m.group(2)),
                    _int(m.group(3)),
                    _int(m.group(4)),
                )
            continue
        if m := _match(IT_NE_BLOCK_RE, line):
            tau, links, group, table, in_desc = (
                m.group(1),
                m.group(2),
                "",
                "",
                False,
            )
            continue
        if m := _match(IT_NE_SET_RE, line):
            group, table, in_desc = m.group(1), "", False
            continue
        if not group:
            continue
        if m := _match(IT_NE_ROW_RE, line):
            if not table:
                raise ReportFormatError(
                    f"{label}: row outside a table: {line}"
                )
            t = tables.setdefault((tau, links, group, table), ClassTable())
            counts = (_int(m.group(2)), _int(m.group(3)), _int(m.group(4)))
            if m.group(1) == "alternating":
                t.alt = counts
            elif m.group(1) == "non-alternating":
                t.nonalt = counts
            else:
                t.noev = counts
            continue
        if in_desc and (m := _match(IT_NE_DESC_RE, line)):
            descriptors.setdefault((tau, links, group), {})[m.group(1)] = (
                _int(m.group(2)),
                _int(m.group(3)),
                _int(m.group(4)),
            )
            continue
        if m := _match(IT_NE_TABLE_RE, line):
            name = m.group(1)
            in_desc = name == IT_NE_DESC_HEAD
            table = "" if in_desc or name.endswith(":") else name
            continue
    _require(tables, "base-class tables", label)
    for key, t in tables.items():
        if not t.complete():
            raise ReportFormatError(f"{label}: incomplete table {key}")
    return ItNoEvidence(tables=tables, descriptors=descriptors, roots=roots)


IT_NE_ITEMS: Final[tuple[tuple[str, str, str, bool, str], ...]] = (
    # (links, set, table, matched, item)
    ("named", "learned set", "all suffixes, all strata", True, "learned set"),
    (
        "named",
        "learned set",
        "all suffixes, internal formations",
        False,
        "learned set, internal formations",
    ),
    (
        "named+candidates",
        "learned set",
        "all suffixes, all strata",
        False,
        "learned set, with generated candidates",
    ),
    (
        "named",
        "native control",
        "all suffixes, all strata",
        False,
        "native control suffixes",
    ),
)


def it_noevidence_rows(text: str, label: str) -> list[Row]:
    parsed = parse_it_noevidence(text, label)
    out: list[Row] = []
    for tau in TAUS:
        counts = parsed.roots.get(tau)
        if counts is None:
            raise ReportFormatError(f"{label}: no root counts at {tau}")
        row = _blank_row("noevidence", "it", "dorsal roots by base class")
        row.update(
            target_class="dorsal",
            tau=tau,
            evidence_set="E0",
            matched=True,
            roots_alternating=counts[0],
            roots_nonalternating=counts[1],
            roots_no_evidence=counts[2],
            note="E0 masculine plurals of the root",
            input=label,
        )
        out.append(row)
        for links, group, table, matched, item in IT_NE_ITEMS:
            t = parsed.tables.get((tau, links, group, table))
            if t is None:
                raise ReportFormatError(
                    f"{label}: no table {table!r} ({links}, {group}) at {tau}"
                )
            row = _blank_row("noevidence", "it", item)
            row.update(
                target_class="dorsal",
                cell=group,
                subset=f"links: {links}; {table}",
                tau=tau,
                evidence_set="E0",
                matched=matched,
                note="units: base root × suffix × derivative stratum; "
                "bases include adjectives (-ico)",
                input=label,
            )
            t.fill(row)
            out.append(row)
        descriptors = parsed.descriptors.get((tau, "named", "learned set"), {})
        _require(descriptors, f"no-evidence descriptors at {tau}", label)
        for name, (pal, plain, both) in descriptors.items():
            row = _blank_row("noevidence", "it", f"no-evidence units: {name}")
            row.update(
                target_class="dorsal",
                cell="learned set",
                subset="links: named; all suffixes, all strata",
                tau=tau,
                evidence_set="E0",
                noev_pal=pal,
                noev_plain=plain,
                n_both=both,
                note="what the no-evidence bases are",
                input=label,
            )
            out.append(row)
    return out


CS_NE_TAU_RE: Final[re.Pattern[str]] = re.compile(r"tau = (\S+)")
CS_NE_ROOTS_RE: Final[re.Pattern[str]] = re.compile(
    r"  (velar|r) noun roots by base class: alternating ([\d,]+), "
    r"non-alternating ([\d,]+), no evidence ([\d,]+)"
)
CS_NE_TARGET_RE: Final[re.Pattern[str]] = re.compile(
    r"  target (velar|r) \(basis wikt\), base x suffix units"
)
CS_NE_INTERNAL_RE: Final[re.Pattern[str]] = re.compile(
    r"  target (velar|r): all P suffixes, derivatives internal formations only"
)
CS_NE_SUFFIX_RE: Final[re.Pattern[str]] = re.compile(r"    suffix (\S.*)")
CS_NE_ROW_RE: Final[re.Pattern[str]] = re.compile(
    r"\s+(alternating|non-alternating|no evidence|\(no alternation evidence\))"
    r"\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)"
)


@dataclass(frozen=True, slots=True)
class CsNoEvidence:
    tables: Mapping[tuple[str, str, str], ClassTable]
    roots: Mapping[tuple[str, str], tuple[int, int, int]]


def parse_cs_noevidence(text: str, label: str) -> CsNoEvidence:
    """``czech_noevidence.txt`` (Wiktionary basis; stops at the NC
    section): (τ, target, table) → class table; roots per base class."""
    tables: dict[tuple[str, str, str], ClassTable] = {}
    roots: dict[tuple[str, str], tuple[int, int, int]] = {}
    tau = target = table = ""
    for line in text.splitlines():
        if line.startswith("NC sources"):
            break
        if m := _match(CS_NE_TAU_RE, line):
            tau, target, table = m.group(1), "", ""
            continue
        if m := _match(CS_NE_ROOTS_RE, line):
            roots[(tau, m.group(1))] = (
                _int(m.group(2)),
                _int(m.group(3)),
                _int(m.group(4)),
            )
            continue
        if m := _match(CS_NE_TARGET_RE, line):
            target, table = m.group(1), ""
            continue
        if m := _match(CS_NE_INTERNAL_RE, line):
            target, table = m.group(1), "internal formations only"
            continue
        if line.startswith("  ") and not line.startswith("   "):
            target, table = "", ""  # another 2-space header: leave
            continue
        if target and (m := _match(CS_NE_SUFFIX_RE, line)):
            table = m.group(1)
            continue
        if target and table and (m := _match(CS_NE_ROW_RE, line)):
            name = m.group(1)
            if name.startswith("("):
                continue  # the pooled row is recomputed from its parts
            t = tables.setdefault((tau, target, table), ClassTable())
            counts = (_int(m.group(2)), _int(m.group(3)), _int(m.group(4)))
            if name == "alternating":
                t.alt = counts
            elif name == "non-alternating":
                t.nonalt = counts
            else:
                t.noev = counts
    _require(tables, "base-class tables", label)
    for key, t in tables.items():
        if not t.complete():
            raise ReportFormatError(f"{label}: incomplete table {key}")
    return CsNoEvidence(tables=tables, roots=roots)


CS_NE_ITEMS: Final[tuple[tuple[str, str, str, bool], ...]] = (
    # (target, table, item, matched)
    ("velar", "all P suffixes", "P suffixes (velar bases)", True),
    (
        "velar",
        "all P suffixes, without -ik(a) → -ický",
        "P suffixes, without learned -ik(a) → -ický",
        False,
    ),
    (
        "velar",
        "internal formations only",
        "P suffixes, derivatives internal formations only",
        False,
    ),
    ("r", "all P suffixes", "P suffixes (r bases)", True),
)
CS_NE_TAUS: Final[tuple[str, ...]] = (
    "dictionary",
    "ge1",
    "ge1_strict",
    "ge2pm",
)


def cs_noevidence_rows(text: str, label: str) -> list[Row]:
    parsed = parse_cs_noevidence(text, label)
    out: list[Row] = []
    for tau in CS_NE_TAUS:
        for target in ("velar", "r"):
            counts = parsed.roots.get((tau, target))
            if counts is None:
                raise ReportFormatError(f"{label}: no {target} roots at {tau}")
            row = _blank_row(
                "noevidence", "cs", f"{target} noun roots by base class"
            )
            row.update(
                target_class=target,
                tau=tau,
                evidence_set="E0",
                matched=tau in TAUS,
                roots_alternating=counts[0],
                roots_nonalternating=counts[1],
                roots_no_evidence=counts[2],
                note="E0 inflection, every P cell",
                input=label,
            )
            out.append(row)
        for target, table, item, matched in CS_NE_ITEMS:
            t = parsed.tables.get((tau, target, table))
            if t is None:
                raise ReportFormatError(
                    f"{label}: no {target} {table!r} at {tau}"
                )
            row = _blank_row("noevidence", "cs", item)
            row.update(
                target_class=target,
                cell="P suffixes",
                subset=f"basis wikt; {table}",
                tau=tau,
                evidence_set="E0",
                matched=matched and tau in TAUS,
                note=(
                    "units: base × P suffix; strict = forms of other "
                    "lexemes and syncretic forms not counted"
                    if tau.endswith("_strict")
                    else "units: base × P suffix"
                ),
                input=label,
            )
            t.fill(row)
            out.append(row)
    return out


# ---------------------------------------------------------------------------
# Section covariation
# ---------------------------------------------------------------------------

MH_RE: Final[re.Pattern[str]] = re.compile(
    r"(?:Mantel-Haenszel|MH) OR (\S+) \[(\S+), (\S+)\] over (\d+) strata"
)
PERM_RE: Final[re.Pattern[str]] = re.compile(
    r"stratified permutation.*?(?:p\(OR_MH >= observed\)|, p) = ([^\s;,]+)"
)
LOGIT_HEAD: Final[str] = "logistic model (R23:"
LOGIT_TERM_RE: Final[re.Pattern[str]] = re.compile(
    r"\bbase_class (\S+) \(se (\S+), p (\S+)\)"
)


@dataclass(frozen=True, slots=True)
class Pooled:
    """The stratified statistics a language reports for one test."""

    mh: tuple[float, float | None, float | None, int] | None = None
    perm_p: float | None = None


@dataclass(frozen=True, slots=True)
class Logit:
    coef: float | None
    se: float | None
    p: float | None
    status: str


def parse_mh(
    line: str,
) -> tuple[float, float | None, float | None, int] | None:
    m = MH_RE.search(line)
    if m is None:
        return None
    return (
        _float(m.group(1)),
        _nan_none(_float(m.group(2))),
        _nan_none(_float(m.group(3))),
        _int(m.group(4)),
    )


def _nan_none(value: float) -> float | None:
    return None if math.isnan(value) else value


def parse_perm(line: str) -> float | None:
    """The permutation p of a 'stratified permutation' line (None when
    it was skipped or not computed)."""
    m = PERM_RE.search(line)
    return None if m is None else _float(m.group(1))


def parse_logit(lines: Sequence[str]) -> Logit:
    """The base_class term of an R23 logistic block (``lines``: the head
    line and the lines after it)."""
    head = lines[0]
    if not head.strip().startswith(LOGIT_HEAD):
        raise ValueError(f"not a logistic model line: {head!r}")
    for marker, status in (
        ("): skipped", "skipped"),
        ("): not estimable", "not estimable"),
        ("): no units", "not estimable"),
    ):
        if marker in head:
            return Logit(None, None, None, status)
    coef = se = p = None
    status = "absent"
    for line in lines[1:]:
        text = line.strip()
        if m := LOGIT_TERM_RE.search(text):
            coef = _opt_float(m.group(1))
            se = _opt_float(m.group(2))
            p = _opt_float(m.group(3))
            status = "estimated"
            continue
        if text.startswith("aliased") and "base_class" in text:
            status = "aliased"
        if text.startswith("quasi-complete separation on"):
            terms = text.split("separation on ", 1)[1].split(":", 1)[0]
            if "base_class" in [t.strip() for t in terms.split(",")]:
                status = "separated"
    return Logit(coef, se, p, status)


def _logit_block(lines: Sequence[str], start: int) -> list[str]:
    """The logistic head line at *start* and its indented continuation."""
    indent = len(lines[start]) - len(lines[start].lstrip())
    block = [lines[start]]
    for line in lines[start + 1 :]:
        if line.strip() and len(line) - len(line.lstrip()) > indent:
            block.append(line)
        else:
            break
    return block


def _table_stats(table: tuple[int, int, int, int]) -> dict[str, float | None]:
    """Sample OR (Woolf interval) and Fisher p of a 2×2 (a, b, c, d)."""
    a, b, c, d = table
    out: dict[str, float | None] = {
        "odds_ratio": None,
        "or_low": None,
        "or_high": None,
        "fisher_p": None,
    }
    if a + b + c + d == 0:
        return out
    grid = ((a, b), (c, d))
    ratio = stats.odds_ratio(grid)
    if ratio is not None:
        out["odds_ratio"] = ratio.estimate
        out["or_low"] = _nan_none(ratio.ci_low)
        out["or_high"] = _nan_none(ratio.ci_high)
    out["fisher_p"] = stats.fisher_exact(grid)
    return out


def covariation_row(
    lang: str,
    item: str,
    table: tuple[int, int, int, int],
    *,
    both: int | None,
    pooled: Pooled,
    logit: Logit | None,
    **descr: Any,
) -> Row:
    row = _blank_row("covariation", lang, item)
    row.update(descr)
    a, b, c, d = table
    row.update(
        alt_pal=a, alt_plain=b, nonalt_pal=c, nonalt_plain=d, n_both=both
    )
    row.update(_table_stats(table))
    if pooled.mh is not None:
        row["mh_or"], row["mh_low"], row["mh_high"], row["mh_strata"] = (
            pooled.mh
        )
    row["perm_p"] = pooled.perm_p
    if logit is not None:
        row["logit_base_class"] = logit.coef
        row["logit_se"] = logit.se
        row["logit_p"] = logit.p
        row["logit_status"] = logit.status
    return row


RO_CV_BLOCK_RE: Final[re.Pattern[str]] = re.compile(
    r"(dorsal|coronal|z|all), tau = (\S+)"
)
RO_CV_VARIANT_RE: Final[re.Pattern[str]] = re.compile(
    r"    (both excluded|both as P|both as N): P:P (\d+) P:N (\d+) "
    r"N:P (\d+) N:N (\d+); .*"
)


@dataclass(frozen=True, slots=True)
class RoCovariation:
    table: tuple[int, int, int, int]
    both: int
    pooled: Pooled
    logit: Logit


def parse_ro_covariation(
    text: str, label: str
) -> dict[tuple[str, str], RoCovariation]:
    """``romanian_covariation.txt``: (class, τ) → the 'both excluded'
    2×2, its stratified statistics and the logistic term."""
    lines = text.splitlines()
    out: dict[tuple[str, str], RoCovariation] = {}
    i = 0
    while i < len(lines):
        m = _match(RO_CV_BLOCK_RE, lines[i])
        if m is None:
            i += 1
            continue
        key = (m.group(1), m.group(2))
        j = i + 1
        table: tuple[int, int, int, int] | None = None
        both = 0
        mh = None
        perm = None
        logit: Logit | None = None
        units_re = re.compile(r"\s+\d+ units: (.*)")
        while j < len(lines) and not _match(RO_CV_BLOCK_RE, lines[j]):
            line = lines[j]
            if (u := _match(units_re, line)) and "both" in u.group(1):
                both = sum(
                    _int(x)
                    for x in re.findall(r"choice both (\d+)", u.group(1))
                )
            v = _match(RO_CV_VARIANT_RE, line)
            if v and v.group(1) == "both excluded":
                table = (
                    _int(v.group(2)),
                    _int(v.group(3)),
                    _int(v.group(4)),
                    _int(v.group(5)),
                )
                k = j + 1
                while k < len(lines) and lines[k].startswith("      "):
                    if "Mantel-Haenszel" in lines[k]:
                        mh = parse_mh(lines[k])
                    if "stratified permutation" in lines[k]:
                        perm = parse_perm(lines[k])
                    k += 1
            if line.strip().startswith(LOGIT_HEAD):
                logit = parse_logit(_logit_block(lines, j))
            j += 1
        if table is None or logit is None:
            raise ReportFormatError(f"{label}: incomplete block {key}")
        out[key] = RoCovariation(table, both, Pooled(mh, perm), logit)
        i = j
    return dict(_require(out, "co-variation blocks", label))


def ro_covariation_rows(text: str, label: str) -> list[Row]:
    parsed = parse_ro_covariation(text, label)
    out: list[Row] = []
    for klass in ("dorsal", "coronal", "all"):
        for tau in TAUS:
            found = parsed.get((klass, tau))
            if found is None:
                raise ReportFormatError(f"{label}: no block {klass} {tau}")
            out.append(
                covariation_row(
                    "ro",
                    "plural class × verbalizer (-i vs -a/-ui)",
                    found.table,
                    both=found.both,
                    pooled=found.pooled,
                    logit=found.logit,
                    target_class=klass,
                    cell="plural -i (alternating) vs -uri/resisting",
                    subset="both excluded",
                    tau=tau,
                    evidence_set="E0",
                    matched=klass == "dorsal",
                    note="base P = alternating plural; derivative P = "
                    "only -i verbs; MH over class × base stratum",
                    input=label,
                )
            )
    return out


IT_CV_BLOCK_RE: Final[re.Pattern[str]] = re.compile(
    r"tau = (\S+); links: (\S+); (learned set|native control)"
)
IT_CV_TABLE_RE: Final[re.Pattern[str]] = re.compile(
    r"  (\S[^:]*?)(?: \(.*\))?:"
)
IT_CV_ROW_RE: Final[re.Pattern[str]] = re.compile(
    r"    (-\S+)\s+(\d+)/(\d+)\s+(\d+)/(\d+)\s+(\d+)\s+.*"
)
IT_CV_MH_RE: Final[re.Pattern[str]] = re.compile(r"    ([^:]+): MH OR (.*)")


@dataclass(frozen=True, slots=True)
class ItCovariation:
    tables: Mapping[str, tuple[int, int, int, int, int]]
    mh: Mapping[str, tuple[float, float | None, float | None, int] | None]
    perm_p: float | None
    logit: Logit


def parse_it_covariation(
    text: str, label: str
) -> dict[tuple[str, str, str], ItCovariation]:
    """``italian_covariation.txt``: (τ, links, set) → summed tables
    (P:P, P:N, N:P, N:N, both) per table name, the pooled MH lines, the
    permutation p and the logistic term."""
    lines = text.splitlines()
    out: dict[tuple[str, str, str], ItCovariation] = {}
    i = 0
    while i < len(lines):
        m = _match(IT_CV_BLOCK_RE, lines[i])
        if m is None:
            i += 1
            continue
        key = (m.group(1), m.group(2), m.group(3))
        sums: dict[str, list[int]] = {}
        mh: dict[str, tuple[float, float | None, float | None, int] | None] = (
            {}
        )
        perm: float | None = None
        logit: Logit | None = None
        table = ""
        in_pooled = False
        j = i + 1
        while j < len(lines) and not _match(IT_CV_BLOCK_RE, lines[j]):
            line = lines[j]
            if line.startswith("  ") and not line.startswith("   "):
                t = _match(IT_CV_TABLE_RE, line)
                in_pooled = line.strip() == "pooled:"
                table = t.group(1) if t and not in_pooled else ""
            elif table and (r := _match(IT_CV_ROW_RE, line)):
                acc = sums.setdefault(table, [0, 0, 0, 0, 0])
                for k in range(5):
                    acc[k] += _int(r.group(k + 2))
            elif in_pooled and (p := _match(IT_CV_MH_RE, line)):
                mh[p.group(1)] = parse_mh(line)
            elif in_pooled and "stratified permutation" in line:
                perm = parse_perm(line)
            elif in_pooled and line.strip().startswith(LOGIT_HEAD):
                logit = parse_logit(_logit_block(lines, j))
            j += 1
        if logit is None or not sums:
            raise ReportFormatError(f"{label}: incomplete block {key}")
        out[key] = ItCovariation(
            tables={
                k: (v[0], v[1], v[2], v[3], v[4]) for k, v in sums.items()
            },
            mh=mh,
            perm_p=perm,
            logit=logit,
        )
        i = j
    return dict(_require(out, "co-variation blocks", label))


IT_CV_ITEMS: Final[tuple[tuple[str, str, str, str, bool, bool], ...]] = (
    # (set, table, MH line, item, with perm/logit, matched)
    (
        "learned set",
        "internal formations",
        "internal formations",
        "learned suffixes, internal formations",
        True,
        True,
    ),
    (
        "learned set",
        "internal formations without -ic- bases",
        "internal, -ic- bases dropped",
        "learned suffixes, internal formations, no -ic- bases",
        False,
        False,
    ),
    (
        "native control",
        "internal formations",
        "internal formations",
        "native control suffixes, internal formations",
        True,
        False,
    ),
)


def it_covariation_rows(text: str, label: str) -> list[Row]:
    parsed = parse_it_covariation(text, label)
    out: list[Row] = []
    for group, table, mh_name, item, full, matched in IT_CV_ITEMS:
        for tau in TAUS:
            found = parsed.get((tau, "named", group))
            if found is None or table not in found.tables:
                raise ReportFormatError(
                    f"{label}: no table {table!r} ({group}) at {tau}"
                )
            if mh_name not in found.mh:
                raise ReportFormatError(f"{label}: no MH line {mh_name!r}")
            a, b, c, d, both = found.tables[table]
            out.append(
                covariation_row(
                    "it",
                    item,
                    (a, b, c, d),
                    both=both,
                    pooled=Pooled(
                        found.mh[mh_name], found.perm_p if full else None
                    ),
                    logit=found.logit if full else None,
                    target_class="dorsal",
                    cell="masculine plural -ci/-gi vs -chi/-ghi",
                    subset="links: named; summed over suffixes",
                    tau=tau,
                    evidence_set="E0",
                    matched=matched,
                    note="units: base root × suffix × derivative stratum; "
                    "OR and Fisher p of the 2×2 summed over suffixes; MH "
                    "over suffix × stratum",
                    input=label,
                )
            )
    return out


CS_CV_BLOCK_RE: Final[re.Pattern[str]] = re.compile(
    r"(LOC\.SG|NOM\.PL) class, tau = (\S+), base coding (\S+)"
)
CS_CV_GROUP_RE: Final[re.Pattern[str]] = re.compile(r"  (\S+): (\d+) units.*")
CS_CV_VARIANT_RE: Final[re.Pattern[str]] = re.compile(
    r"    (both excluded|both as P|both as N):"
)
CS_CV_PAIR_RE: Final[re.Pattern[str]] = re.compile(
    r"      (\S.*?)\s+P:P\s+(\d+) P:N\s+(\d+) N:P\s+(\d+) N:N\s+(\d+)\s+OR .*"
)


@dataclass(frozen=True, slots=True)
class CsCovariation:
    pairs: Mapping[str, tuple[int, int, int, int]]
    pooled: Pooled
    logit: Logit | None
    both: int


def parse_cs_covariation(
    text: str, label: str
) -> dict[tuple[str, str, str, str], CsCovariation]:
    """``czech_covariation.txt``: (cell, τ, coding, group) → the 'both
    excluded' pair tables, MH and permutation, the logistic term."""
    lines = text.splitlines()
    out: dict[tuple[str, str, str, str], CsCovariation] = {}
    cell = tau = coding = group = variant = ""
    pairs: dict[str, tuple[int, int, int, int]] = {}
    mh = None
    perm: float | None = None
    logit: Logit | None = None
    both = 0

    def close() -> None:
        if group:
            out[(cell, tau, coding, group)] = CsCovariation(
                pairs=dict(pairs),
                pooled=Pooled(mh, perm),
                logit=logit,
                both=both,
            )

    for j, line in enumerate(lines):
        if line.startswith("Roots with a pair derivative"):
            break
        if m := _match(CS_CV_BLOCK_RE, line):
            close()
            cell, tau, coding = m.group(1), m.group(2), m.group(3)
            group = ""
            continue
        if m := _match(CS_CV_GROUP_RE, line):
            close()
            group, variant = m.group(1), ""
            pairs, mh, perm, logit = {}, None, None, None
            b = re.search(r"both (\d+)", line)
            both = _int(b.group(1)) if b else 0
            continue
        if not group:
            continue
        if m := _match(CS_CV_VARIANT_RE, line):
            variant = m.group(1)
            continue
        if line.strip().startswith(LOGIT_HEAD):
            logit = parse_logit(_logit_block(lines, j))
            variant = ""
            continue
        if variant != "both excluded":
            continue
        if m := _match(CS_CV_PAIR_RE, line):
            pairs[m.group(1)] = (
                _int(m.group(2)),
                _int(m.group(3)),
                _int(m.group(4)),
                _int(m.group(5)),
            )
        elif "Mantel-Haenszel" in line:
            mh = parse_mh(line)
        elif "stratified permutation" in line:
            perm = parse_perm(line)
    close()
    return dict(_require(out, "co-variation blocks", label))


CS_CV_ITEMS: Final[tuple[tuple[str, str, str, bool], ...]] = (
    # (cell, coding, group, matched): the Czech headline tests (LOC.SG
    # coded 'any': under 'dominant' no velar base is P; NOM.PL coded
    # 'dominant'), the NOM.PL 'any' block for its logistic model, and
    # the baseline of hard stems without alternation.
    ("LOC.SG", "any", "velar", True),
    ("LOC.SG", "any", "baseline", False),
    ("NOM.PL", "dominant", "velar", True),
    ("NOM.PL", "any", "velar", False),
    ("NOM.PL", "dominant", "baseline", False),
)


def cs_covariation_rows(text: str, label: str) -> list[Row]:
    parsed = parse_cs_covariation(text, label)
    out: list[Row] = []
    for cell, coding, group, matched in CS_CV_ITEMS:
        for tau in TAUS:
            found = parsed.get((cell, tau, coding, group))
            if found is None:
                raise ReportFormatError(
                    f"{label}: no block {cell} {tau} {coding} {group}"
                )
            descr: dict[str, Any] = {
                "target_class": group,
                "cell": f"{cell} class (base coding {coding})",
                "tau": tau,
                "evidence_set": "E0",
                "input": label,
            }
            pairs = list(found.pairs.items())
            _require(pairs, f"pairs of {cell} {tau} {coding} {group}", label)
            for pair, table in pairs:
                out.append(
                    covariation_row(
                        "cs",
                        f"{cell} class ({coding}) × {pair}",
                        table,
                        both=None,
                        pooled=Pooled(),
                        logit=None,
                        subset="both excluded",
                        matched=False,
                        note="one competing pair",
                        **descr,
                    )
                )
            summed = tuple(sum(t[k] for _, t in pairs) for k in range(4))
            out.append(
                covariation_row(
                    "cs",
                    f"{cell} class ({coding}) × suffix choice (all pairs)",
                    (summed[0], summed[1], summed[2], summed[3]),
                    both=found.both,
                    pooled=found.pooled,
                    logit=found.logit,
                    subset="both excluded; summed over pairs",
                    matched=matched,
                    note="OR and Fisher p of the 2×2 summed over pairs; MH "
                    "over pair × base stratum"
                    + (
                        "; baseline = labial, l, s, z hard stems (no "
                        "alternation involved)"
                        if group == "baseline"
                        else ""
                    ),
                    **descr,
                )
            )
    return out


# ---------------------------------------------------------------------------
# Section logtable (R27)
# ---------------------------------------------------------------------------


def ro_logtable_rows(table: TPTable, label: str) -> list[Row]:
    out: list[Row] = []
    for klass in ("person", "thing"):
        for tau in TAUS:
            pal, ret = table.pair(
                table="logtable",
                target_class="g:-log/-fag/-urg",
                trigger="-i/-e",
                subset=klass,
                tau=tau,
            )
            out.append(
                tp_pair_row(
                    "logtable",
                    "ro",
                    f"-log/-fag/-urg {klass}",
                    pal,
                    ret,
                    input_path=label,
                    matched=True,
                    note="plural -i/-e; shielded = -uri",
                )
            )
    return out


IT_LOG_SUFFIXOIDS: Final[tuple[str, ...]] = ("-logo", "-fago", "-urgo")


def it_logtable_rows(table: TPTable, label: str) -> list[Row]:
    out: list[Row] = []
    for klass in ("person", "thing"):
        for tau in TAUS:
            parts = [
                table.pair(
                    table="logtable",
                    target_class="dorsal",
                    trigger="-i",
                    subset=f"{suffixoid}:{klass}",
                    tau=tau,
                )
                for suffixoid in IT_LOG_SUFFIXOIDS
            ]
            n_target = sum(_int(p["N_target"]) for p, _ in parts)
            n = sum(_int(p["N_trigger"]) for p, _ in parts)
            e_pal = sum(_int(p["e"]) for p, _ in parts)
            e_ret = sum(_int(r["e"]) for _, r in parts)
            row = _blank_row("logtable", "it", f"-logo/-fago/-urgo {klass}")
            row.update(
                target_class="dorsal",
                cell="-i",
                subset=klass,
                tau=tau,
                evidence_set="E0",
                coding="dominant",
                freq_basis="lemma",
                matched=True,
                N_target=n_target,
                N=n,
                n_shielded=sum(_int(p["n_shielded"]) for p, _ in parts),
                e_palatalize=e_pal,
                e_retain=e_ret,
                theta=_theta_of(n),
                tolerated_palatalize=tp.tolerated(n, e_pal),
                tolerated_retain=tp.tolerated(n, e_ret),
                fragility_palatalize=tp.fragility(n, e_pal),
                fragility_retain=tp.fragility(n, e_ret),
                note="the -logo, -fago and -urgo rows summed (disjoint "
                "roots); θ and verdicts computed here",
                input=label,
            )
            row["verdict"] = two_rule_verdict(
                row["tolerated_palatalize"], row["tolerated_retain"]
            )
            out.append(row)
            for suffixoid, (pal, ret) in zip(
                IT_LOG_SUFFIXOIDS, parts, strict=True
            ):
                if _int(pal["N_target"]) == 0:
                    continue
                out.append(
                    tp_pair_row(
                        "logtable",
                        "it",
                        f"{suffixoid} {klass}",
                        pal,
                        ret,
                        input_path=label,
                        note="one suffixoid",
                    )
                )
    return out


CS_LOG_SUFFIXOIDS: Final[frozenset[str]] = frozenset({"-log", "-urg"})
CS_LOG_GE1_HEAD: Final[str] = "tau = ge1: allomorphs of the attested cells"
CS_LOG_ALLO_RE: Final[re.Pattern[str]] = re.compile(
    r"    (-\S+)\s+(\S+)\s+(\S+)\s+(\d+)"
)
CS_LOG_PAL_RE: Final[re.Pattern[str]] = re.compile(r"  (\S+)\s+(.+)")


def _cs_log_counts(rows: Sequence[CsvRow], klass: str) -> dict[str, int]:
    """Dictionary-level NOM.PL and LOC.PL counts of the learned -log/-urg
    nouns of one class (``czech_logtable.csv``)."""
    counts = dict.fromkeys(
        (
            "nouns",
            "nom_front",
            "nom_shielded",
            "nom_pal",
            "nom_plain",
            "loc_front",
            "loc_shielded",
            "loc_pal",
            "loc_plain",
        ),
        0,
    )
    for row in rows:
        if (
            row["learned_suffixoid"] != "True"
            or row["suffixoid"] not in CS_LOG_SUFFIXOIDS
            or row["semantic_class"] != klass
        ):
            continue
        counts["nouns"] += 1
        allomorphs = csvio.split_list(row["nom_pl_allomorph"])
        forms = csvio.split_list(row["nom_pl"])
        if len(forms) != len(allomorphs):
            # The logtable lists each ending once, so two forms with the
            # same ending would shift the index below onto another form.
            raise ReportFormatError(
                f"czech_logtable.csv: {row['lemma']}: NOM.PL forms "
                f"{forms} and endings {allomorphs} are not aligned"
            )
        if "-i" in allomorphs:
            counts["nom_front"] += 1
            form = forms[allomorphs.index("-i")]
            if form == row["lemma"] + "i":
                counts["nom_plain"] += 1
            else:
                counts["nom_pal"] += 1
        elif allomorphs:
            counts["nom_shielded"] += 1
        loc = csvio.split_list(row["loc_pl_allomorph"])
        if "-ích" in loc or "-ech" in loc:
            counts["loc_front"] += 1
            if row["loc_pl_palatalized"] == "True":
                counts["loc_pal"] += 1
            else:
                counts["loc_plain"] += 1
        elif loc:
            counts["loc_shielded"] += 1
    return counts


def parse_cs_logtable_ge1(
    text: str, label: str
) -> tuple[dict[tuple[str, str, str], dict[str, int]], dict[str, int]]:
    """``czech_logtable_summary.txt`` at τ = ge1: (cell, suffixoid,
    class) → {allomorph string: nouns}, and the pooled -log/-urg LOC.PL
    palatalized counts per class."""
    lines = text.splitlines()
    allomorphs: dict[tuple[str, str, str], dict[str, int]] = {}
    palatalized: dict[str, int] = {}
    try:
        start = lines.index(CS_LOG_GE1_HEAD)
    except ValueError:
        raise ReportFormatError(f"{label}: no '{CS_LOG_GE1_HEAD}'") from None
    cell = ""
    for line in lines[start + 1 :]:
        if not line.strip() or line.startswith("tau = "):
            break
        if line.startswith("  ") and not line.startswith("   "):
            cell = line.strip()
            continue
        if m := _match(CS_LOG_ALLO_RE, line):
            key = (cell, m.group(1), m.group(2))
            allomorphs.setdefault(key, {})[m.group(3)] = _int(m.group(4))
    found = False
    for line in lines:
        if line.startswith("LOC.PL palatalized"):
            found = True
            continue
        if (
            found
            and (m := _match(CS_LOG_PAL_RE, line))
            and m.group(1) == "ge1"
        ):
            for item in m.group(2).split(","):
                klass, _, value = item.strip().rpartition(" ")
                name, _, flag = klass.partition(":")
                if flag == "True":
                    palatalized[name] = palatalized.get(name, 0) + _int(value)
            break
    if not found:
        raise ReportFormatError(f"{label}: no LOC.PL palatalized line")
    return allomorphs, palatalized


def cs_logtable_rows(
    table_rows: Sequence[CsvRow],
    summary: str,
    inputs: Inputs,
) -> list[Row]:
    out: list[Row] = []
    csv_label = inputs.relative(inputs.cs_logtable)
    text_label = inputs.relative(inputs.cs_logtable_summary)
    ge1, palatalized = parse_cs_logtable_ge1(summary, text_label)
    for klass in ("person", "thing"):
        c = _cs_log_counts(table_rows, klass)
        for cell, prefix in (("NOM.PL m.an", "nom"), ("LOC.PL", "loc")):
            row = _blank_row("logtable", "cs", f"-log/-urg {klass}")
            n = c[f"{prefix}_front"]
            row.update(
                target_class="velar",
                cell=cell,
                subset=klass,
                tau="dictionary",
                evidence_set="E0",
                matched=True,
                N_target=c["nouns"],
                N=n,
                n_shielded=c[f"{prefix}_shielded"],
                e_palatalize=c[f"{prefix}_plain"],
                e_retain=c[f"{prefix}_pal"],
                theta=_theta_of(n),
                tolerated_palatalize=tp.tolerated(n, c[f"{prefix}_plain"]),
                tolerated_retain=tp.tolerated(n, c[f"{prefix}_pal"]),
                note=(
                    "N = nouns with -i; shielded = -ové or -y only"
                    if prefix == "nom"
                    else "N = nouns with -ích or -ech"
                ),
                input=csv_label,
            )
            row["verdict"] = two_rule_verdict(
                row["tolerated_palatalize"], row["tolerated_retain"]
            )
            out.append(row)
        for cell, key, front in (
            ("NOM.PL m.an", "nom_pl", ("-i",)),
            ("LOC.PL", "loc_pl", ("-ích", "-ech")),
        ):
            n = shielded = 0
            for suffixoid in sorted(CS_LOG_SUFFIXOIDS):
                for allomorph, count in ge1.get(
                    (key, suffixoid, klass), {}
                ).items():
                    if any(a in front for a in allomorph.split("|")):
                        n += count
                    else:
                        shielded += count
            row = _blank_row("logtable", "cs", f"-log/-urg {klass}")
            row.update(
                target_class="velar",
                cell=cell,
                subset=klass,
                tau="ge1",
                evidence_set="E0",
                matched=True,
                N_target=c["nouns"],
                N=n,
                n_shielded=shielded,
                note="counts of czech_logtable_summary.txt at ≥ 1 token",
                input=text_label,
            )
            if key == "loc_pl":
                pal = palatalized.get(klass, 0)
                if pal > n:
                    raise ReportFormatError(
                        f"{text_label}: {pal} palatalized of {n} LOC.PL"
                    )
                row["e_retain"] = pal
                row["e_palatalize"] = n - pal
            elif n == 0:
                row["e_retain"] = row["e_palatalize"] = 0
            if row["e_palatalize"] is not None:
                row["theta"] = _theta_of(n)
                row["tolerated_palatalize"] = tp.tolerated(
                    n, row["e_palatalize"]
                )
                row["tolerated_retain"] = tp.tolerated(n, row["e_retain"])
                row["verdict"] = two_rule_verdict(
                    row["tolerated_palatalize"], row["tolerated_retain"]
                )
            out.append(row)
    return out


# ---------------------------------------------------------------------------
# Section phonology (LP validators)
# ---------------------------------------------------------------------------

VALIDATION_ROWS_RE: Final[re.Pattern[str]] = re.compile(
    r"Target-stem rows validated: ([\d,]+)"
)
VALIDATION_BOOL_RE: Final[re.Pattern[str]] = re.compile(
    r"BOOLEAN mutation match: ([\d,]+)/([\d,]+)\s+\(.*\)"
)
VALIDATION_IPA_RE: Final[re.Pattern[str]] = re.compile(
    r"SURFACE-IPA match: ([\d,]+)/([\d,]+)\s+\(.*\)"
)


def parse_validation(text: str, label: str) -> tuple[int, int, int]:
    """(rows, boolean matches, IPA matches) of a validator report."""
    rows = bool_match = ipa_match = None
    for line in text.splitlines()[:8]:
        if m := _match(VALIDATION_ROWS_RE, line):
            rows = _int(m.group(1))
        elif m := _match(VALIDATION_BOOL_RE, line):
            bool_match = _int(m.group(1))
        elif m := _match(VALIDATION_IPA_RE, line):
            ipa_match = _int(m.group(1))
    if rows is None or bool_match is None or ipa_match is None:
        raise ReportFormatError(f"{label}: no validator header")
    return rows, bool_match, ipa_match


PHONOLOGY_NOTES: Final[Mapping[tuple[str, str], str]] = {
    ("ro", "class"): "every dorsal, and every coronal before -i, "
    "underspecified given the attested plural suffix (the draft's rules; "
    "legacy validator)",
    ("it", "lexical"): "underspecified iff the root alternates (reads "
    "the -i outcome: circular for the -i cell)",
    ("it", "stress"): "underspecified iff the velar is not post-tonic "
    "(no lexical listing)",
    ("it", "all_underspecified"): "every dorsal underspecified",
    ("cs", "lexical"): "underspecified iff the lexeme alternates in E0 "
    "(the tested cell is its own evidence: circular)",
    ("cs", "all_underspecified"): "every target underspecified",
    ("cs", "lexical_loo"): "Inkelas learner: underspecified only with "
    "alternation evidence outside the tested cell",
    ("cs", "spellout"): "prespecified iff the lexeme takes a "
    "non-palatalizing ending somewhere (G&W's spell-out claim as a policy)",
}


def phonology_rows(inputs: Inputs) -> list[Row]:
    out: list[Row] = []
    specs: list[tuple[str, str, Path]] = [
        ("ro", "class", inputs.ro_validation)
    ]
    specs += [("it", v, inputs.it_validation(v)) for v in IT_VARIANTS]
    specs += [("cs", v, inputs.cs_validation(v)) for v in CS_VARIANTS]
    for lang, variant, path in specs:
        label = inputs.relative(path)
        n, bool_match, ipa_match = parse_validation(
            path.read_text(encoding="utf-8"), label
        )
        row = _blank_row("phonology", lang, f"LP validator: {variant}")
        row.update(
            cell="target-stem rows",
            subset=variant,
            N=n,
            n_match=bool_match,
            n_ipa_match=ipa_match,
            note=PHONOLOGY_NOTES[(lang, variant)],
            input=label,
        )
        out.append(row)
    return out


# ---------------------------------------------------------------------------
# Section legacy (Romanian pipeline, read only)
# ---------------------------------------------------------------------------

LEGACY_RE: Final[re.Pattern[str]] = re.compile(
    r"<(c|g|t|d|s|z|st|ct|sc)> \+ <-(i|e)> plural( \(downsampled\))?"
)
LEGACY_SEGMENTS: Final[tuple[str, ...]] = (
    "c",
    "g",
    "t",
    "d",
    "s",
    "z",
    "st",
    "ct",
)


def legacy_rows(inputs: Inputs) -> list[Row]:
    out: list[Row] = []
    for path, frame, top in (
        (inputs.ro_legacy_top, "legacy top-1000", True),
        (inputs.ro_legacy_full, "legacy full lexicon", False),
    ):
        label = inputs.relative(path)
        found: dict[tuple[str, str], CsvRow] = {}
        for src in csvio.read_rows(path):
            m = _match(LEGACY_RE, src["type"])
            if m is None or bool(m.group(3)) != top:
                continue
            found[(m.group(1), m.group(2))] = src
        for segment in LEGACY_SEGMENTS:
            for vowel in ("i", "e"):
                src_row = found.get((segment, vowel))
                if src_row is None:
                    continue  # the legacy table has no such cell (st + -e)
                n = _int(src_row["N"])
                mutated = _int(src_row["mutated"])
                row = _blank_row(
                    "legacy", "ro", f"<{segment}> + <-{vowel}> plural"
                )
                row.update(
                    target_class=segment,
                    cell=f"-{vowel}",
                    subset=frame,
                    tau="dictionary",
                    V="1000" if top else "",
                    freq_basis="lemma" if top else "",
                    N=n,
                    e_palatalize=n - mutated,
                    e_retain=mutated,
                    theta=_theta_of(n),
                    tolerated_palatalize=_bool_cell(src_row["tolerated?"]),
                    tolerated_retain=tp.tolerated(n, mutated),
                    note="analysis/analyze_romanian_palatalization.R "
                    "(the draft's productivity table is the top-1000 frame); "
                    "tolerated_retain computed here",
                    input=label,
                )
                row["verdict"] = two_rule_verdict(
                    row["tolerated_palatalize"], row["tolerated_retain"]
                )
                out.append(row)
    return out


# ---------------------------------------------------------------------------
# Building the summary
# ---------------------------------------------------------------------------


def build_rows(root: Path = PROJECT_ROOT) -> list[Row]:
    """Every summary row, in section order (deterministic)."""
    inputs = Inputs(root)
    missing = [inputs.relative(p) for p in inputs.all() if not p.is_file()]
    if missing:
        raise FileNotFoundError(
            "cross-language inputs missing (run the language pipelines "
            "first): " + ", ".join(missing)
        )
    tables = {
        "ro": TPTable.read(inputs.ro_tp, inputs.relative(inputs.ro_tp)),
        "it": TPTable.read(inputs.it_tp, inputs.relative(inputs.it_tp)),
        "cs": TPTable.read(inputs.cs_tp, inputs.relative(inputs.cs_tp)),
    }
    learner_tables = {
        "ro": csvio.read_rows(inputs.ro_learners),
        "it": csvio.read_rows(inputs.it_learners),
        "cs": csvio.read_rows(inputs.cs_learners),
    }

    def text(path: Path) -> str:
        return path.read_text(encoding="utf-8")

    rel = inputs.relative
    rows: list[Row] = []
    rows += r17_rows(tables, inputs)
    rows += learner_rows(learner_tables, inputs)
    rows += ro_model_rows(tables["ro"], rel(inputs.ro_tp))
    rows += it_model_rows(tables["it"], rel(inputs.it_tp))
    rows += cs_model_rows(tables["cs"], rel(inputs.cs_tp))
    rows += ro_spellout_rows(tables["ro"], text(inputs.ro_spellout), inputs)
    rows += it_spellout_rows(tables["it"], rel(inputs.it_tp))
    rows += cs_spellout_rows(text(inputs.cs_spellout), rel(inputs.cs_spellout))
    rows += ro_noevidence_rows(
        text(inputs.ro_noevidence), rel(inputs.ro_noevidence)
    )
    rows += it_noevidence_rows(
        text(inputs.it_noevidence), rel(inputs.it_noevidence)
    )
    rows += cs_noevidence_rows(
        text(inputs.cs_noevidence), rel(inputs.cs_noevidence)
    )
    rows += ro_covariation_rows(
        text(inputs.ro_covariation), rel(inputs.ro_covariation)
    )
    rows += it_covariation_rows(
        text(inputs.it_covariation), rel(inputs.it_covariation)
    )
    rows += cs_covariation_rows(
        text(inputs.cs_covariation), rel(inputs.cs_covariation)
    )
    rows += ro_logtable_rows(tables["ro"], rel(inputs.ro_tp))
    rows += it_logtable_rows(tables["it"], rel(inputs.it_tp))
    rows += cs_logtable_rows(
        csvio.read_rows(inputs.cs_logtable),
        text(inputs.cs_logtable_summary),
        inputs,
    )
    rows += phonology_rows(inputs)
    rows += legacy_rows(inputs)
    check_rows(rows)
    return rows


def check_rows(rows: Sequence[Row]) -> None:
    """Every row has a known section and language, every column, and a
    unique identity; sections appear in :data:`SECTIONS` order."""
    order = {name: i for i, name in enumerate(SECTIONS)}
    last = -1
    seen: set[tuple[Any, ...]] = set()
    identity = (
        "section",
        "lang",
        "item",
        "target_class",
        "cell",
        "subset",
        "tau",
        "evidence_set",
        "coding",
        "freq_basis",
        "V",
        "rule",
    )
    for row in rows:
        if set(row) != set(FIELDS):
            raise ValueError(f"row columns differ from FIELDS: {row['item']}")
        if row["section"] not in order or row["lang"] not in LANGS:
            raise ValueError(f"unknown section or language: {row}")
        position = order[row["section"]]
        if position < last:
            raise ValueError(f"section {row['section']} out of order")
        last = position
        key = tuple(row[name] for name in identity)
        if key in seen:
            raise ValueError(f"duplicate summary row: {key}")
        seen.add(key)


# ---------------------------------------------------------------------------
# CSV text (the bytes csvio.write_rows writes)
# ---------------------------------------------------------------------------


def csv_text(rows: Sequence[Row]) -> str:
    """The summary as :func:`lexpipe.csvio.write_rows` writes it."""
    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator=csvio.LINE_TERMINATOR)
    writer.writerow(FIELDS)
    for row in rows:
        writer.writerow([csvio.format_value(row.get(name)) for name in FIELDS])
    return buffer.getvalue()


# ---------------------------------------------------------------------------
# README block
# ---------------------------------------------------------------------------

BEGIN: Final[str] = (
    "<!-- BEGIN crosslang: generated by lexpipe.compare from "
    "crosslang_summary.csv; edit the code, not this block -->"
)
END: Final[str] = "<!-- END crosslang -->"
# The Romanian section's marker (lexpipe.langs.romanian_compare.readme):
# a new block goes in front of it.
ROMANIAN_SECTION_PREFIX: Final[str] = "<!-- BEGIN romanian_compare"
README_TITLE: Final[str] = (
    "# Cross-language comparison: Romanian, Italian, Czech\n"
)
README_SOURCES: Final[tuple[str, ...]] = (
    "wiktionary_cs",
    "wiktionary_en",
    "wiktionary_ro",
)


def merge_block(existing: str | None, body: str) -> str:
    """*existing* README text with the generated block set to *body*.

    Text outside the markers is kept byte for byte. Without markers the
    block goes in front of the Romanian section (or at the end).
    """
    block = f"{BEGIN}\n{body.rstrip()}\n{END}"
    if existing is None or not existing.strip():
        return f"{README_TITLE}\n{block}\n"
    start = existing.find(BEGIN)
    if start >= 0:
        end = existing.find(END, start)
        if end < 0:
            raise ValueError("README: crosslang BEGIN marker without END")
        return existing[:start] + block + existing[end + len(END) :]
    romanian = existing.find(ROMANIAN_SECTION_PREFIX)
    if romanian >= 0:
        return existing[:romanian] + block + "\n\n" + existing[romanian:]
    separator = (
        ""
        if existing.endswith("\n\n")
        else ("\n" if existing.endswith("\n") else "\n\n")
    )
    return existing + separator + block + "\n"


def extract_block(text: str) -> str | None:
    """The body between the markers, or None without them."""
    start = text.find(BEGIN)
    if start < 0:
        return None
    end = text.find(END, start)
    if end < 0:
        return None
    return text[start + len(BEGIN) : end].strip("\n") + "\n"


def _n(value: Any) -> str:
    return "–" if value is None or value == "" else f"{int(value):,}"


def _theta(value: Any) -> str:
    return "–" if value is None or value == "" else f"{float(value):.1f}"


def _share(value: Any) -> str:
    return f"{float(value):.2f}"


def _p(value: Any) -> str:
    if value is None or value == "":
        return "–"
    number = float(value)
    if number == 0:
        return "0"
    if number < 1e-3:
        return f"{number:.2e}"
    return f"{number:.3g}"


def _perm(value: Any) -> str:
    """Permutation p (10,000 shuffles): 0 means none reached the
    observed statistic."""
    if value is None or value == "":
        return "–"
    return "< 0.0001" if float(value) == 0 else f"{float(value):.3g}"


def _ratio(value: Any) -> str:
    if value is None or value == "":
        return "–"
    number = float(value)
    if math.isinf(number):
        return "∞"
    return f"{number:.3g}"


def _ci(low: Any, high: Any) -> str:
    if (low is None or low == "") and (high is None or high == ""):
        return ""
    return f" [{_ratio(low)}, {_ratio(high)}]"


def _pct(pal: Any, plain: Any) -> str:
    a = int(pal or 0)
    b = int(plain or 0)
    if a + b == 0:
        return "0 / 0"
    return f"{a:,} / {b:,} ({100 * a / (a + b):.0f}%)"


def _verdict(value: Any) -> str:
    text = "" if value is None else str(value)
    if text in ("palatalize", "retain", "both", "conditioned"):
        return f"**{text}**"
    return text or "–"


def _tp_cell(row: Row | None) -> str:
    """'N: e(pal) / e(ret), θ → verdict' of a two-rule row."""
    if row is None:
        return "–"
    return (
        f"{_n(row['N'])}: {_n(row['e_palatalize'])} / "
        f"{_n(row['e_retain'])}, θ {_theta(row['theta'])} → "
        f"{_verdict(row['verdict'])}"
    )


def _one_rule_cell(row: Row | None) -> str:
    """'N, e (θ) ✓/✗' of a single-rule row."""
    if row is None:
        return "–"
    rule = row["rule"]
    mark = {"tolerated": "✓", "not tolerated": "✗"}.get(row["verdict"], "–")
    return (
        f"{_n(row['N'])}, e {_n(row[f'e_{rule}'])} "
        f"(θ {_theta(row['theta'])}) {mark}"
    )


def _table(header: Sequence[str], body: Sequence[Sequence[str]]) -> list[str]:
    lines = ["| " + " | ".join(header) + " |"]
    lines.append("|" + "|".join("---" for _ in header) + "|")
    lines += ["| " + " | ".join(cells) + " |" for cells in body]
    return lines


def _select(rows: Sequence[Row], **where: Any) -> list[Row]:
    return [r for r in rows if all(r[k] == v for k, v in where.items())]


def _one(rows: Sequence[Row], **where: Any) -> Row | None:
    found = _select(rows, **where)
    if len(found) > 1:
        raise ValueError(f"{len(found)} summary rows match {where}")
    return found[0] if found else None


R17_ORDER: Final[tuple[str, ...]] = (R17A, R17B, R17B_NT, R17C, R17C_NT)
TARGET_LABELS: Final[Mapping[tuple[str, str], str]] = {
    ("ro", "coronal"): "t, d, s",
    ("ro", "z"): "z",
    ("it", "coronal"): "t, d, s, z",
    ("cs", "coronal_inferred"): "t, d, n (inferred)",
    ("cs", "r"): "r",
    ("cs", "sibilant"): "s, z",
    ("cs", "lateral"): "l",
    ("cs", "velar"): "velar",
}


def _target(lang: str, target: str) -> str:
    return TARGET_LABELS.get((lang, target), target)


def _r17_table(rows: Sequence[Row], **where: Any) -> list[str]:
    selected = _select(rows, section="r17", **where)
    groups: dict[tuple[str, str, str], dict[str, Row]] = {}
    for row in selected:
        key = (row["item"], row["lang"], row["target_class"])
        groups.setdefault(key, {})[row["tau"]] = row
    keys = sorted(
        groups,
        key=lambda k: (R17_ORDER.index(k[0]), LANGS.index(k[1])),
    )
    body = [
        [
            item,
            LANG_NAMES[lang],
            _target(lang, target),
            _tp_cell(groups[(item, lang, target)].get("dictionary")),
            _tp_cell(groups[(item, lang, target)].get("ge1")),
            " / ".join(
                _n(groups[(item, lang, target)][t]["n_shielded"])
                for t in TAUS
                if t in groups[(item, lang, target)]
            ),
        ]
        for item, lang, target in keys
    ]
    return _table(
        [
            "rule-family row",
            "language",
            "target",
            "τ = dictionary: N: e (pal) / e (ret), θ → productive",
            "τ ≥ 1 token",
            "shielded (dict. / ≥ 1)",
        ],
        body,
    )


def _shares(row: Row | None) -> str:
    if row is None:
        return "–"
    parts = [
        (float(row[f"share_{name}"]), name)
        for name in LEARNER_OUTCOMES
        if float(row[f"share_{name}"]) >= 0.005
    ]
    parts.sort(key=lambda item: (-item[0], LEARNER_OUTCOMES.index(item[1])))
    return ", ".join(f"{name} {share:.3f}" for share, name in parts)


def _learner_table(rows: Sequence[Row], evidence_set: str) -> list[str]:
    selected = _select(rows, section="learners", evidence_set=evidence_set)
    groups: dict[tuple[str, str, str], dict[str, Row]] = {}
    first: dict[tuple[str, str], int] = {}
    for row in selected:
        key = (row["lang"], row["item"], row["tau"])
        groups.setdefault(key, {})[row["V"]] = row
        first.setdefault((row["lang"], row["item"]), len(first))
    keys = sorted(
        groups,
        key=lambda k: (
            LANGS.index(k[0]),
            first[(k[0], k[1])],
            TAUS.index(k[2]),
        ),
    )
    body = [
        [
            LANG_NAMES[lang],
            item,
            tau,
            _shares(groups[(lang, item, tau)].get("1000")),
            _shares(groups[(lang, item, tau)].get("5000")),
        ]
        for lang, item, tau in keys
    ]
    return _table(
        ["language", "cell", "τ", "V = 1000", "V = 5000"],
        body,
    )


def _nested_sets(rows: Sequence[Row], lang: str) -> list[list[str]]:
    """One table row per nested set and τ: its branches and grammar."""
    out: list[list[str]] = []
    selected = _select(rows, section="model_fit", lang=lang)
    for grammar in [r for r in selected if r["item"].endswith(": grammar")]:
        title = grammar["item"][: -len(": grammar")]
        branches = [
            r
            for r in selected
            if r["item"].startswith(title + ": ")
            and r is not grammar
            and r["tau"] == grammar["tau"]
            and r["evidence_set"] == grammar["evidence_set"]
            and r["rule"]
        ]
        text = "; ".join(
            f"{b['item'][len(title) + 2 :]}: {_one_rule_cell(b)}"
            for b in branches
        )
        out.append(
            [
                title,
                grammar["tau"],
                grammar["evidence_set"],
                text,
                _verdict(grammar["verdict"]),
            ]
        )
    return out


def _model_tables(rows: Sequence[Row]) -> list[str]:
    lines: list[str] = []
    for lang in LANGS:
        sets = _nested_sets(rows, lang)
        if sets:
            lines += [
                "",
                f"Nested rule sets, {LANG_NAMES[lang]} (Yang's order; each "
                "branch: N, e (θ), ✓ = tolerated):",
                "",
            ]
            lines += _table(
                ["rule set", "τ", "evidence", "branches", "grammar"], sets
            )
    singles: list[list[str]] = []
    for row in _select(rows, section="model_fit"):
        if row["item"].endswith(": grammar") or row["subset"].startswith(
            "nested"
        ):
            continue
        if row["rule"] and str(row["note"]).startswith("non-trigger"):
            cell = _one_rule_cell(row)
            productive = "non-trigger suffix: N = e"
        elif row["rule"]:
            cell = _one_rule_cell(row)
            productive = (
                _verdict(row["rule"])
                if row["verdict"] == "tolerated"
                else "not tolerated"
            )
        else:
            cell = _tp_cell(row)
            productive = _verdict(row["verdict"])
        singles.append(
            [
                LANG_NAMES[row["lang"]],
                row["item"],
                _target(row["lang"], row["target_class"]),
                row["tau"],
                row["evidence_set"],
                cell,
                _n(row["n_shielded"]),
                productive,
            ]
        )
    lines += [
        "",
        "Single rules and pooled rows (a two-rule row reads N: e (pal) / e "
        "(ret), θ → productive). In a Czech competition row *palatalize* "
        "is the rule 'take the palatalizing allomorph' (e = roots taking "
        "the competitor) and *retain* the rule 'take the competitor'; the "
        "baseline is the hard labial, l, s and z stems, where no "
        "alternation is involved:",
        "",
    ]
    lines += _table(
        [
            "language",
            "rule",
            "target",
            "τ",
            "evidence",
            "N, e (θ)",
            "shielded",
            "verdict",
        ],
        singles,
    )
    return lines


def _spellout_table(rows: Sequence[Row]) -> list[str]:
    groups: dict[tuple[str, str, str, str], dict[str, Row]] = {}
    for row in _select(rows, section="spellout"):
        key = (row["lang"], row["item"], row["target_class"], row["subset"])
        groups.setdefault(key, {})[row["tau"]] = row
    body: list[list[str]] = []
    for (lang, item, target, subset), cells in groups.items():
        dictionary = cells.get("dictionary")
        if dictionary is not None and dictionary["coding"] == "raw":
            only_a = dictionary["e_palatalize"]
            only_i = dictionary["e_retain"]
            text = (
                f"{_n(dictionary['N'])} roots: non-palatalizing member only "
                f"{_n(only_a)}, -i only {_n(only_i)}, both "
                f"{_n(dictionary['n_both'])}"
            )
            body.append(
                [
                    LANG_NAMES[lang],
                    item,
                    _target(lang, target),
                    subset,
                    text,
                    "–",
                ]
            )
            continue
        body.append(
            [
                LANG_NAMES[lang],
                item,
                _target(lang, target),
                subset,
                _tp_cell(dictionary),
                _tp_cell(cells.get("ge1")),
            ]
        )
    return _table(
        [
            "language",
            "pair",
            "target",
            "roots",
            "τ = dictionary: N: e (pal) / e (ret), θ → productive",
            "τ ≥ 1 token",
        ],
        body,
    )


def _noevidence_tables(rows: Sequence[Row]) -> list[str]:
    roots_body: list[list[str]] = []
    units_body: list[list[str]] = []
    desc_body: list[list[str]] = []
    for row in _select(rows, section="noevidence"):
        if row["roots_alternating"] is not None:
            total = (
                row["roots_alternating"]
                + row["roots_nonalternating"]
                + row["roots_no_evidence"]
            )
            roots_body.append(
                [
                    LANG_NAMES[row["lang"]],
                    _target(row["lang"], row["target_class"]),
                    row["tau"],
                    _n(row["roots_alternating"]),
                    _n(row["roots_nonalternating"]),
                    f"{_n(row['roots_no_evidence'])} "
                    f"({100 * row['roots_no_evidence'] / total:.0f}%)",
                ]
            )
        elif row["item"].startswith("no-evidence units: "):
            desc_body.append(
                [
                    row["item"][len("no-evidence units: ") :],
                    row["tau"],
                    _pct(row["noev_pal"], row["noev_plain"]),
                ]
            )
        else:
            avoid = (
                f"{_n(row['alt_avoid'])} / {_n(row['nonalt_avoid'])} / "
                f"{_n(row['noev_avoid'])}"
                if row["alt_avoid"] is not None
                else "–"
            )
            units_body.append(
                [
                    LANG_NAMES[row["lang"]],
                    row["item"],
                    _target(row["lang"], row["target_class"]),
                    row["tau"],
                    _pct(row["alt_pal"], row["alt_plain"]),
                    _pct(row["nonalt_pal"], row["nonalt_plain"]),
                    _pct(row["noev_pal"], row["noev_plain"]),
                    _n(row["n_both"]),
                    avoid,
                ]
            )
    lines = [
        "Roots by base class (alternation evidence in inflection, E0):",
        "",
    ]
    lines += _table(
        [
            "language",
            "target",
            "τ",
            "alternating",
            "non-alternating",
            "no evidence",
        ],
        roots_body,
    )
    lines += [
        "",
        "Derivative units by base class (palatalized / not, % palatalized):",
        "",
    ]
    lines += _table(
        [
            "language",
            "units",
            "target",
            "τ",
            "alternating",
            "non-alternating",
            "no evidence",
            "both",
            "back-vowel suffix instead (alt. / non-alt. / no ev.)",
        ],
        units_body,
    )
    lines += [
        "",
        "Italian no-evidence units (learned set, named links) by what the "
        "base is:",
        "",
    ]
    lines += _table(["base", "τ", "palatalized / velar"], desc_body)
    return lines


def _logit(row: Row) -> str:
    status = row["logit_status"]
    if status is None or status == "":
        return "–"
    if status == "estimated":
        return (
            f"{_ratio(row['logit_base_class'])} (se "
            f"{_ratio(row['logit_se'])}, p {_p(row['logit_p'])})"
        )
    if row["logit_base_class"] is not None:
        return f"{_ratio(row['logit_base_class'])} ({status})"
    return str(status)


def _covariation_table(rows: Sequence[Row]) -> list[str]:
    body = [
        [
            LANG_NAMES[row["lang"]],
            row["item"],
            _target(row["lang"], row["target_class"]),
            row["tau"],
            f"{_n(row['alt_pal'])} / {_n(row['alt_plain'])}",
            f"{_n(row['nonalt_pal'])} / {_n(row['nonalt_plain'])}",
            _ratio(row["odds_ratio"]) + _ci(row["or_low"], row["or_high"]),
            _p(row["fisher_p"]),
            (
                _ratio(row["mh_or"])
                + _ci(row["mh_low"], row["mh_high"])
                + (f" ({row['mh_strata']})" if row["mh_strata"] else "")
                if row["mh_or"] is not None
                else "–"
            ),
            _perm(row["perm_p"]),
            _logit(row),
        ]
        for row in _select(rows, section="covariation")
        if "(all pairs)" in row["item"] or row["lang"] != "cs"
    ]
    return _table(
        [
            "language",
            "test",
            "base",
            "τ",
            "P base: P / N",
            "N base: P / N",
            "OR [95% CI]",
            "Fisher p",
            "MH OR [95% CI] (strata)",
            "perm. p",
            "logit base_class",
        ],
        body,
    )


def _logtable_table(rows: Sequence[Row]) -> list[str]:
    groups: dict[tuple[str, str, str], dict[str, Row]] = {}
    for row in _select(rows, section="logtable", matched=True):
        key = (row["lang"], row["item"], row["cell"])
        groups.setdefault(key, {})[row["tau"]] = row
    body: list[list[str]] = []
    for (lang, item, cell), cells in groups.items():
        first = next(iter(cells.values()))
        texts: list[str] = []
        for tau in TAUS:
            found = cells.get(tau)
            if found is None:
                texts.append("–")
                continue
            e_text = (
                f"{_n(found['e_palatalize'])} / {_n(found['e_retain'])}"
                if found["e_palatalize"] is not None
                else "?"
            )
            texts.append(
                f"{_n(found['N'])} before a front vowel ({e_text}), "
                f"{_n(found['n_shielded'])} shielded → "
                f"{_verdict(found['verdict'])}"
            )
        body.append(
            [LANG_NAMES[lang], item, cell, _n(first["N_target"]), *texts]
        )
    return _table(
        [
            "language",
            "nouns",
            "cell",
            "N",
            "τ = dictionary: front (e pal / e ret), shielded",
            "τ ≥ 1 token",
        ],
        body,
    )


def _phonology_table(rows: Sequence[Row]) -> list[str]:
    body = [
        [
            LANG_NAMES[row["lang"]],
            row["subset"],
            _n(row["N"]),
            f"{100 * row['n_match'] / row['N']:.2f}%",
            f"{100 * row['n_ipa_match'] / row['N']:.2f}%",
            row["note"],
        ]
        for row in _select(rows, section="phonology")
    ]
    return _table(
        ["language", "policy", "rows", "boolean match", "IPA match", "policy"],
        body,
    )


def _legacy_table(rows: Sequence[Row]) -> list[str]:
    body: list[list[str]] = []
    for segment in LEGACY_SEGMENTS:
        cells = []
        for vowel in ("i", "e"):
            for frame in ("legacy top-1000", "legacy full lexicon"):
                row = _one(
                    rows,
                    section="legacy",
                    target_class=segment,
                    cell=f"-{vowel}",
                    subset=frame,
                )
                cells.append(
                    "–"
                    if row is None
                    else f"{_n(row['N'])}: {_n(row['e_palatalize'])} / "
                    f"{_n(row['e_retain'])} → {_verdict(row['verdict'])}"
                )
        body.append([f"<{segment}>", *cells])
    return _table(
        [
            "stem final",
            "-i, top-1000 (draft)",
            "-i, full lexicon",
            "-e, top-1000 (draft)",
            "-e, full lexicon",
        ],
        body,
    )


def render_block(rows: Sequence[Row]) -> str:
    """The generated tables of the cross-language README."""
    lines = [
        "Unit: root types. θ = N/ln N; e (pal) = exceptions to "
        "*palatalize* (roots attested before the trigger that keep the "
        "plain consonant); e (ret) = exceptions to *retain* (roots that "
        "change); productive = the rule(s) with e ≤ θ (**both** only for "
        "very small N). τ = dictionary counts every listed form; τ ≥ 1 "
        "token counts forms with at least one corpus token. Shielded = "
        "roots attested in the cell only with a non-front or "
        "non-palatalizing allomorph. Matched settings: E0, dominant "
        "doublet coding, lemma frequency basis, nouns; every column of "
        "every row is in the CSV.",
        "",
        "### Table 1. The R17 rule family at matched settings",
        "",
        *_r17_table(rows, matched=True),
        "",
        "### Table 1b. The same in the R24 frame of the 1000 most "
        "frequent noun roots",
        "",
        *_r17_table_frame(rows),
        "",
        "### Table 1c. Italian nouns and adjectives together (the Italian "
        "headline subset)",
        "",
        *_r17_table(rows, lang="it", subset="all", V=""),
        "",
        "### Table 1d. Without suffix targets (subset core; Italian "
        "nouns and adjectives)",
        "",
        *_r17_table(rows, subset="core", V=""),
        "",
        "### Table 2. Simulated learners (1000 vocabularies by token "
        "frequency, E0, lemma basis): shares of learners by productive "
        "rule",
        "",
        *_learner_table(rows, "E0"),
        "",
        "### Table 2b. Simulated learners with verb cells added (E1)",
        "",
        *_learner_table(rows, "E1"),
        "",
        "### Table 3. Language-specific model fits",
        *_model_tables(rows),
        "",
        "### Table 4. Spell-out on certain-specification roots (R20)",
        "",
        "e (pal) = roots taking only the non-palatalizing member (doublets "
        "by the dominant variant); e (ret) = roots taking the palatalizing "
        "member. Contrast rows: roots without alternation (prespecified "
        "for the draft).",
        "",
        *_spellout_table(rows),
        "",
        "### Table 5. No alternation evidence (R21)",
        "",
        *_noevidence_tables(rows),
        "",
        "### Table 6. Co-variation of inflectional class and derivative "
        "choice (R23)",
        "",
        "P base = the base alternates (takes the palatalizing allomorph); "
        "P/N = the derivative palatalizes (takes the palatalizing member) "
        "or not; units with both outcomes left out. OR and Fisher p of the "
        "2×2 (summed over pairs or suffixes); MH = Mantel–Haenszel over "
        "each language's strata; perm. p = stratified permutation "
        "(10,000 shuffles); logit = base_class in deriv_pal ~ base_class "
        "+ log_lemma_freq + base_stratum + deriv_stratum, base-clustered "
        "errors.",
        "",
        *_covariation_table(rows),
        "",
        "### Table 7. -log/-fag/-urg nouns (R27)",
        "",
        *_logtable_table(rows),
        "",
        "### Table 8. Logical Phonology validators (share of target-stem "
        "rows whose derived surface form matches)",
        "",
        *_phonology_table(rows),
        "",
        "### Table 9. Romanian legacy tables (the draft's productivity "
        "table and the full lexicon): N: e (pal) / e (ret) → productive",
        "",
        *_legacy_table(rows),
    ]
    return "\n".join(lines) + "\n"


def _r17_table_frame(rows: Sequence[Row]) -> list[str]:
    selected = [
        r
        for r in rows
        if r["section"] == "r17"
        and r["V"] == "1000"
        and not (r["lang"] == "it" and r["subset"] == "all")
    ]
    return _r17_table(selected)


# ---------------------------------------------------------------------------
# Stage, check and command line
# ---------------------------------------------------------------------------


def write_outputs(ctx: RunContext, rows: Sequence[Row]) -> Outputs:
    """Write the CSV (licence gate, R9) and the README block."""
    out = outputs(ctx.root)
    csvio.write_tracked(out.csv, rows, FIELDS, ctx=ctx)
    existing = (
        out.readme.read_text(encoding="utf-8")
        if out.readme.is_file()
        else None
    )
    ctx.write_tracked_text(
        out.readme,
        merge_block(existing, render_block(rows)),
        sources=README_SOURCES,
    )
    return out


def context(root: Path = PROJECT_ROOT, *, force: bool = False) -> RunContext:
    """The run context: the comparison-only registry row ``ro``, so the
    step may write only under analysis/comparison/ and data/romanian_compare/.
    """
    return RunContext(
        LANGUAGES["ro"],
        root=root,
        force=frozenset({STAGE}) if force else frozenset(),
    )


def stage_params() -> dict[str, str]:
    """The manifest parameters: the version and this module's SHA-256,
    so a code change re-runs the step (it takes about a second)."""
    return {"version": VERSION, "code_sha256": sha256_file(Path(__file__))}


def run(ctx: RunContext, log: Callable[[str], None] = print) -> int | None:
    """Build and write the summary unless the manifest matches (R6)."""
    inputs = Inputs(ctx.root)
    out = outputs(ctx.root)
    missing = [inputs.relative(p) for p in inputs.all() if not p.is_file()]
    if missing:
        log(
            "  crosslang: inputs missing (run the language pipelines "
            "first): " + ", ".join(missing)
        )
        return 1

    def body() -> int:
        rows = build_rows(ctx.root)
        write_outputs(ctx, rows)
        log(
            f"  crosslang: {len(rows):,} rows → "
            f"{ctx.relative(out.csv)}; tables → {ctx.relative(out.readme)}"
        )
        return 0

    return run_stage(
        ctx,
        STAGE,
        inputs=inputs.all(),
        params=stage_params(),
        outputs=[out.csv, out.readme],
        body=body,
        log=log,
    )


def stale(root: Path = PROJECT_ROOT) -> list[str]:
    """Which tracked outputs differ from a fresh build ([] = current)."""
    rows = build_rows(root)
    out = outputs(root)
    problems: list[str] = []
    if not out.csv.is_file() or out.csv.read_text(
        encoding="utf-8"
    ) != csv_text(rows):
        problems.append(Inputs(root).relative(out.csv))
    readme = (
        out.readme.read_text(encoding="utf-8") if out.readme.is_file() else ""
    )
    if extract_block(readme) != render_block(rows):
        problems.append(
            Inputs(root).relative(out.readme) + " (generated block)"
        )
    return problems


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m lexpipe.compare",
        description="Build analysis/comparison/crosslang_summary.csv and "
        "the generated tables of analysis/comparison/README.md from the "
        "tracked Romanian, Italian and Czech outputs.",
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=PROJECT_ROOT,
        metavar="DIR",
        help="tree to work in (default: this repository)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="rebuild even if the manifest says the outputs are current",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="write nothing; exit 1 if the tracked outputs are stale",
    )
    args = parser.parse_args(argv)
    if args.check:
        problems = stale(args.root)
        for problem in problems:
            print(f"stale: {problem}")
        if problems:
            print("rebuild with: PYTHONPATH=src python -m lexpipe.compare")
        return 1 if problems else 0
    status = run(context(args.root, force=args.force))
    return 0 if status in (None, 0) else 1


if __name__ == "__main__":
    sys.exit(main())
