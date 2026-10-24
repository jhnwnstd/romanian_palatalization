"""Tests for the cross-language synthesis (lexpipe.compare).

The parsers are tested on excerpts in the exact format of each tracked
report (copied from the reports, trimmed); the end-to-end tests run the
step on a copy of the tracked inputs in a temporary tree and check the
real outputs for consistency with the per-language reports.
"""

from __future__ import annotations

import hashlib
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lexpipe import compare, csvio, licensing  # noqa: E402
from lexpipe.context import read_manifest  # noqa: E402
from lexpipe.langs.romanian_compare import readme as ro_readme  # noqa: E402

REAL = compare.Inputs(PROJECT_ROOT)
HAVE_REAL = all(path.is_file() for path in REAL.all())
needs_tracked = pytest.mark.skipif(
    not HAVE_REAL, reason="tracked language outputs are absent"
)


# ---------------------------------------------------------------------------
# Helpers and verdicts
# ---------------------------------------------------------------------------


def tp_csv_row(**over: str) -> dict[str, str]:
    row = {
        "lang": "ro",
        "table": "R17a",
        "rule": "palatalize",
        "target_class": "dorsal",
        "trigger": "-i/-e",
        "subset": "all",
        "tau": "dictionary",
        "evidence_set": "E0",
        "coding": "dominant",
        "freq_basis": "lemma",
        "V": "",
        "N_target": "100",
        "N_trigger": "50",
        "e_nonalt": "3",
        "n_shielded": "20",
        "e": "3",
        "theta_trigger": "12.78",
        "theta_target": "21.71",
        "tolerated_trigger": "True",
        "tolerated_counterfactual": "False",
        "majority": "True",
        "memo": "dorsal/_[+front V]",
        "fragility": "10",
        "source": "wiktionary_ro|wiktionary_en",
        "licence_class": "by-sa",
    }
    row.update(over)
    return row


@pytest.mark.parametrize(
    ("pal", "ret", "expected"),
    [
        (True, False, "palatalize"),
        (False, True, "retain"),
        (True, True, "both"),
        (False, False, "neither"),
        (None, None, ""),
    ],
)
def test_two_rule_verdict(
    pal: bool | None, ret: bool | None, expected: str
) -> None:
    assert compare.two_rule_verdict(pal, ret) == expected


def test_tp_pair_row_merges_both_rules() -> None:
    pal = tp_csv_row()
    ret = tp_csv_row(
        rule="retain", e="47", tolerated_trigger="False", fragility="35"
    )
    row = compare.tp_pair_row(
        "r17", "ro", "R17a", pal, ret, input_path="x.csv", matched=True
    )
    assert (row["N"], row["e_palatalize"], row["e_retain"]) == (50, 3, 47)
    assert row["verdict"] == "palatalize"
    assert (
        row["n_shielded"] == 20 and row["counterfactual_palatalize"] is False
    )
    assert (row["fragility_palatalize"], row["fragility_retain"]) == (10, 35)
    assert set(row) == set(compare.FIELDS)


def test_pair_lookup_refuses_mismatched_rules() -> None:
    table = compare.TPTable(
        [tp_csv_row(), tp_csv_row(rule="retain", N_trigger="51")], "t.csv"
    )
    with pytest.raises(compare.ReportFormatError, match="differ"):
        table.pair(
            table="R17a",
            target_class="dorsal",
            trigger="-i/-e",
            subset="all",
            tau="dictionary",
        )


def test_lookup_needs_exactly_one_row() -> None:
    table = compare.TPTable([tp_csv_row(), tp_csv_row()], "t.csv")
    with pytest.raises(compare.ReportFormatError, match="2 rows match"):
        table.find(
            table="R17a",
            rule="palatalize",
            target_class="dorsal",
            trigger="-i/-e",
            subset="all",
            tau="dictionary",
        )


def test_grammar_row_is_conditioned_only_if_every_branch_holds() -> None:
    good = compare.tp_single_row(
        "model_fit", "it", "set: a", tp_csv_row(), input_path="x"
    )
    bad = compare.tp_single_row(
        "model_fit",
        "it",
        "set: b",
        tp_csv_row(e="30", tolerated_trigger="False"),
        input_path="x",
    )
    assert good["verdict"] == "tolerated" and good["e_palatalize"] == 3
    assert bad["verdict"] == "not tolerated"
    ok = compare.grammar_row("model_fit", "it", "set: grammar", [good, good])
    assert ok["verdict"] == "conditioned" and ok["N"] == 100
    no = compare.grammar_row("model_fit", "it", "set: grammar", [good, bad])
    assert no["verdict"] == "not conditioned"


# ---------------------------------------------------------------------------
# Statistic lines shared by the reports
# ---------------------------------------------------------------------------


def test_parse_mh_and_permutation_lines() -> None:
    assert compare.parse_mh(
        "      Mantel-Haenszel OR 5.34 [2.29, 12.5] over 3 strata "
        "(class x base stratum; statsmodels)"
    ) == (5.34, 2.29, 12.5, 3)
    inf = compare.parse_mh(
        "    internal formations: MH OR inf [nan, nan] over 6 strata (lexpipe)"
    )
    assert (
        inf is not None and math.isinf(inf[0]) and inf[1:] == (None, None, 6)
    )
    assert (
        compare.parse_mh("      Mantel-Haenszel: undefined (no ...)") is None
    )
    assert (
        compare.parse_perm(
            "      stratified permutation (seed 20080715, 10,000/10,000 "
            "valid): p(OR_MH >= observed) = 0.22; null median 1.19, 97.5% 4.21"
        )
        == 0.22
    )
    assert (
        compare.parse_perm(
            "    stratified permutation (base class shuffled within suffix "
            "x stratum; 10000/10000 valid shuffles, seed 20080715): observed "
            "MH OR inf, p = 0 (< 0.0001: no shuffle reached it), null median "
            "0.964, 97.5% 2.77"
        )
        == 0.0
    )
    assert (
        compare.parse_perm(
            "    stratified permutation (...): observed MH OR inf, p = 0.0009,"
            " null median 0, 97.5% 105"
        )
        == 0.0009
    )
    assert (
        compare.parse_perm("      stratified permutation: skipped (OR_MH 0)")
        is None
    )


@pytest.mark.parametrize(
    ("lines", "expected"),
    [
        (
            [
                "    logistic model (R23: deriv_pal ~ base_class; 151 "
                "units, 151 bases, reference strata base=unknown, "
                "deriv=unknown):",
                "      Intercept -0.811 (se 0.468, p 0.0827); base_class "
                "1.57 (se 0.459, p 0.000641); log_lemma_freq -0.108 (se "
                "0.0738, p 0.144)",
                "      one level, no term: deriv_stratum (all unknown)",
            ],
            (1.57, 0.459, 0.000641, "estimated"),
        ),
        (
            [
                "    logistic model (R23: deriv_pal ~ base_class; 285 "
                "units, 244 bases, reference strata base=x, deriv=y):",
                "      Intercept -2.47 (se NA, p NA); base_class 45.7 (se "
                "NA, p NA)",
                "      quasi-complete separation on base_class, "
                "base_stratum[unknown]: the maximum-likelihood estimate "
                "does not exist",
            ],
            (45.7, None, None, "separated"),
        ),
        (
            [
                "    logistic model (R23: deriv_pal ~ base_class; 8 units, "
                "8 bases, reference strata base=x, deriv=y):",
                "      Intercept 0.315 (se 1.5, p 0.833)",
                "      aliased (constant or collinear, NA): base_class",
            ],
            (None, None, None, "aliased"),
        ),
        (
            [
                "    logistic model (R23: deriv_pal ~ base_class): skipped: "
                "no statsmodels"
            ],
            (None, None, None, "skipped"),
        ),
    ],
)
def test_parse_logit(lines: list[str], expected: tuple[Any, ...]) -> None:
    fit = compare.parse_logit(lines)
    assert (fit.coef, fit.se, fit.p, fit.status) == expected


def test_parse_logit_ignores_separation_of_other_terms() -> None:
    fit = compare.parse_logit(
        [
            "  logistic model (R23: deriv_pal ~ base_class; 282 units):",
            "    Intercept 3.03 (se 0.64, p 2.26e-06); base_class -0.0495 "
            "(se 0.604, p 0.935)",
            "    quasi-complete separation on deriv_stratum[borrowed whole]:"
            " the maximum-likelihood estimate does not exist",
        ]
    )
    assert fit.status == "estimated" and fit.coef == -0.0495


# ---------------------------------------------------------------------------
# Report parsers (excerpts in the reports' format)
# ---------------------------------------------------------------------------

RO_NOEVIDENCE = """\
tau = dictionary
------------------------------------------------------------------------
  target roots by base class: alternating 5492, non-alternating 420, \
non-alternating -uri 1716, no evidence 8106
    class dorsal (base x suffix units)
      base  verb -i  -iza/-ifica  adj. front  all front  both  avoid
                                pal/nonpal  pal/nonpal  pal/nonpal  pal/nonpal
      alternating  85/0  2/0  42/0  129/0  0  92
      non-alternating  0/0  0/0  0/0  0/0  0  2
      non-alternating -uri  8/0  0/0  2/0  10/0  0  44
      no evidence  8/0  0/0  3/0  11/0  0  19
      (non-alt. pooled)  8/0  0/0  2/0  10/0  0  46
  By base stratum (all classes):
    stratum borrowed
      base  verb -i  -iza/-ifica  adj. front  all front  both  avoid
      alternating  15/4  0/1  13/6  28/11  0  29
"""


def test_parse_ro_noevidence() -> None:
    parsed = compare.parse_ro_noevidence(RO_NOEVIDENCE, "ro.txt")
    assert parsed.roots["dictionary"] == (5492, 420, 1716, 8106)
    dorsal = parsed.tables[("dictionary", "dorsal")]
    assert dorsal["alternating"] == (129, 0, 0, 92)
    assert dorsal["(non-alt. pooled)"] == (10, 0, 0, 46)
    assert dorsal["no evidence"] == (11, 0, 0, 19)
    # The stratum tables are not class tables.
    assert set(parsed.tables) == {("dictionary", "dorsal")}


IT_NOEVIDENCE = """\
Dorsal roots by class:
  dictionary alternating 5575, non-alternating 1037, no evidence 1120
  ge1        alternating 980, non-alternating 322, no evidence 6430

tau = dictionary; links: named
  ----------------------------------------------------------------------
  learned set:
    all suffixes, all strata
      base                 P     N  both
      alternating        280     0     0
      non-alternating      9    40     0
      no evidence         10     4     0
    -issimo, all strata
      base                 P     N  both
      alternating          9     0     0
      non-alternating      1    13     0
      no evidence          0     0     0
    no-evidence units by base descriptor (P/N/both):
      -ic- suffix base                        7/   0/  0
      neither                                 3/   4/  0
    no-evidence examples (words only):
      palatalized [10]: armonica → armonicista (P)
  native control:
    all suffixes, all strata
      base                 P     N  both
      alternating          1    31     0
      non-alternating      4   120     0
      no evidence          6    77     0
"""


def test_parse_it_noevidence() -> None:
    parsed = compare.parse_it_noevidence(IT_NOEVIDENCE, "it.txt")
    assert parsed.roots == {
        "dictionary": (5575, 1037, 1120),
        "ge1": (980, 322, 6430),
    }
    learned = parsed.tables[
        ("dictionary", "named", "learned set", "all suffixes, all strata")
    ]
    assert (learned.alt, learned.nonalt, learned.noev) == (
        (280, 0, 0),
        (9, 40, 0),
        (10, 4, 0),
    )
    native = parsed.tables[
        ("dictionary", "named", "native control", "all suffixes, all strata")
    ]
    assert native.noev == (6, 77, 0)
    assert parsed.descriptors[("dictionary", "named", "learned set")] == {
        "-ic- suffix base": (7, 0, 0),
        "neither": (3, 4, 0),
    }


CS_NOEVIDENCE = """\
tau = dictionary
------------------------------------------------------------------------
  velar noun roots by base class: alternating 7,330, non-alternating 447, \
no evidence 518
  r noun roots by base class: alternating 730, non-alternating 747, no \
evidence 135
  target velar (basis wikt), base x suffix units
    suffix all P suffixes
      base class                      pal  plain   both   % pal
      alternating                    1472      8      0   99.5%
      non-alternating                  32      2      0   94.1%
      no evidence                      27      0      1   96.4%
      (no alternation evidence)        59      2      1   95.2%
  target velar: by base stratum (all P suffixes)
    base stratum inherited
      base class                      pal  plain   both   % pal
      alternating                     424      3      0   99.3%
  target velar: all P suffixes, derivatives internal formations only
    base class                      pal  plain   both   % pal
    alternating                     357      4      0   98.9%
    non-alternating                   3      1      0   75.0%
    no evidence                       5      0      0  100.0%
    (no alternation evidence)         8      1      0   88.9%
  sensitivity (velar, all P suffixes):
      alternating                     999      9      9   99.0%

NC sources (aggregates only): basis all (DeriNet-only pairs added),
  tau = dictionary, velar, all P suffixes
  target velar (basis wikt), base x suffix units
    suffix all P suffixes
      alternating                    5000      0      0  100.0%
      non-alternating                5000      0      0  100.0%
      no evidence                    5000      0      0  100.0%
"""


def test_parse_cs_noevidence_reads_wiktionary_basis_only() -> None:
    parsed = compare.parse_cs_noevidence(CS_NOEVIDENCE, "cs.txt")
    assert parsed.roots[("dictionary", "velar")] == (7330, 447, 518)
    main = parsed.tables[("dictionary", "velar", "all P suffixes")]
    assert (main.alt, main.nonalt, main.noev) == (
        (1472, 8, 0),
        (32, 2, 0),
        (27, 0, 1),
    )
    internal = parsed.tables[
        ("dictionary", "velar", "internal formations only")
    ]
    assert internal.noev == (5, 0, 0)
    # Stratum and sensitivity tables, and the NC section, are never read.
    assert set(parsed.tables) == {
        ("dictionary", "velar", "all P suffixes"),
        ("dictionary", "velar", "internal formations only"),
    }


RO_COVARIATION = """\
dorsal, tau = dictionary
------------------------------------------------------------------------
  167 units: base N choice N 28, base N choice P 9, base N choice both 2, \
base P choice N 42, base P choice P 72, base P choice both 14
    both excluded: P:P 72 P:N 42 N:P 9 N:N 28; OR 5.33 [2.3, 12.4], \
Fisher p = 5.01e-05
      Mantel-Haenszel OR 5.34 [2.29, 12.5] over 3 strata (class x base \
stratum; statsmodels)
      stratified permutation (seed 20080715, 10,000/10,000 valid): \
p(OR_MH >= observed) = 0; null median 0.98, 97.5% 2.04
    both as P: P:P 86 P:N 42 N:P 11 N:N 28; OR 5.21 [2.37, 11.5], Fisher \
p = 2.97e-05
      Mantel-Haenszel OR 5.1 [2.31, 11.2] over 3 strata (class x base \
stratum; statsmodels)
      stratified permutation (seed 20080715, 10,000/10,000 valid): \
p(OR_MH >= observed) = 0.5; null median 0.952, 97.5% 2.14
    logistic model (R23: deriv_pal ~ base_class + log_lemma_freq + \
base_stratum + deriv_stratum; 151 units, 151 bases, reference strata \
base=unknown, deriv=unknown):
      Intercept -0.811 (se 0.468, p 0.0827); base_class 1.57 (se 0.459, \
p 0.000641); log_lemma_freq -0.108 (se 0.0738, p 0.144)
      one level, no term: deriv_stratum (all unknown)
    base P choosing -a/-ui: [42] aplică → aplica
"""


def test_parse_ro_covariation() -> None:
    parsed = compare.parse_ro_covariation(RO_COVARIATION, "ro.txt")
    found = parsed[("dorsal", "dictionary")]
    assert found.table == (72, 42, 9, 28)
    assert found.both == 16
    assert found.pooled.mh == (5.34, 2.29, 12.5, 3)
    assert found.pooled.perm_p == 0.0  # the 'both excluded' variant's
    assert (found.logit.coef, found.logit.status) == (1.57, "estimated")


IT_COVARIATION = """\
tau = dictionary; links: named; learned set
  ----------------------------------------------------------------------
  units by derivative stratum: inherited 0, borrowed whole 6, internal \
formation 285, unknown 37
  internal formations (the co-variation set, §5):
    suffix      P-base:P/N  N-base:P/N  both  statistics
    -issimo            1/0         0/1     0  OR inf [nan, inf]; Fisher p = 1
    -izzare  20/0  2/3  0  OR inf [nan, inf]; Fisher p = 0.00435
    -ità  207/0  1/2  0  OR inf [nan, inf]; Fisher p = 0.000137
  internal formations without -ic- bases (251 of 285 units have an -ic- \
base, abulico → abulicità):
    suffix      P-base:P/N  N-base:P/N  both  statistics
    -izzare  0/0  2/3  0  OR undefined (empty row or column); Fisher p = 1
    -ità  5/0  1/2  0  OR inf [nan, inf]; Fisher p = 0.107
  all strata pooled (sensitivity):
    suffix      P-base:P/N  N-base:P/N  both  statistics
    -ità  215/0  3/4  0  OR inf [nan, inf]; Fisher p = 3.55e-07
  pooled:
    internal formations: MH OR inf [nan, nan] over 6 strata (lexpipe)
    internal, -ic- bases dropped: MH OR inf [nan, nan] over 6 strata (lexpipe)
    internal, feminine named bases dropped: MH OR undefined
    stratified permutation (base class shuffled within suffix x stratum; \
10000/10000 valid shuffles, seed 20080715): observed MH OR inf, p = 0 \
(< 0.0001: no shuffle reached it), null median 0.964, 97.5% 2.77
    logistic model (R23: deriv_pal ~ base_class + log_lemma_freq + \
base_stratum + deriv_stratum; 285 units, 244 bases, reference strata \
base=borrowed whole, deriv=internal formation):
      Intercept -2.47 (se NA, p NA); base_class 45.7 (se NA, p NA)
      quasi-complete separation on base_class, base_stratum[unknown]: the \
maximum-likelihood estimate does not exist
  examples (internal formations; words only):
    P base, velar derivative [0]:
"""


def test_parse_it_covariation_sums_suffix_rows() -> None:
    parsed = compare.parse_it_covariation(IT_COVARIATION, "it.txt")
    found = parsed[("dictionary", "named", "learned set")]
    assert found.tables["internal formations"] == (228, 0, 3, 6, 0)
    assert found.tables["internal formations without -ic- bases"] == (
        5,
        0,
        3,
        5,
        0,
    )
    assert found.tables["all strata pooled"] == (215, 0, 3, 4, 0)
    assert found.mh["internal, feminine named bases dropped"] is None
    mh = found.mh["internal formations"]
    assert mh is not None and math.isinf(mh[0]) and mh[3] == 6
    assert found.perm_p == 0.0
    assert (found.logit.coef, found.logit.status) == (45.7, "separated")


CS_COVARIATION = """\
LOC.SG class, tau = dictionary, base coding any
------------------------------------------------------------------------
  velar: 508 units (N 390, P 72, both 46)
    both excluded:
      -ní/-ný vs -ový/-ovní    P:P    4 P:N    5 N:P   47 N:N  344  OR \
5.86 [1.52, 22.58] Fisher p=0.0182
      -it vs -ovat             P:P    0 P:N    1 N:P   21 N:N   40  OR \
0.00 [0.00, nan] Fisher p=1
      Mantel-Haenszel OR 1.96 [0.50, 7.74] over 8 strata (pair x base \
stratum; statsmodels)
      stratified permutation (seed 20080715, 10,000/10,000 valid): \
p(OR_MH >= observed) = 0.22; null median 1.19, 97.5% 4.21
    both as P:
      -ní/-ný vs -ový/-ovní    P:P    4 P:N    5 N:P   89 N:N  344  OR \
3.09 [0.81, 11.75] Fisher p=0.0978
      Mantel-Haenszel OR 1.14 [0.29, 4.48] over 8 strata (pair x base \
stratum; statsmodels)
  logistic model (R23: deriv_pal ~ base_class + log_lemma_freq + \
base_stratum + deriv_stratum; 462 units, 443 bases, reference strata \
base=internal formation, deriv=unknown):
    Intercept -3.23 (se 0.34, p 2.28e-21); base_class -0.282 (se 1.05, p \
0.789); log_lemma_freq 0.31 (se 0.0788, p 8.54e-05)
  baseline: 632 units (N 363, P 221, both 48)
    both excluded:
      -ní/-ný vs -ový/-ovní    P:P   83 P:N   89 N:P  113 N:N  208  OR \
1.72 [1.18, 2.50] Fisher p=0.00519
      Mantel-Haenszel OR 1.63 [1.12, 2.36] over 8 strata (pair x base \
stratum; statsmodels)
      stratified permutation (seed 20080715, 10,000/10,000 valid): \
p(OR_MH >= observed) = 0.01; null median 1.01, 97.5% 1.46

Roots with a pair derivative, by target class and stratum:
  velar: inherited 254, borrowed whole 117, internal formation 430
"""


def test_parse_cs_covariation_keeps_the_both_excluded_variant() -> None:
    parsed = compare.parse_cs_covariation(CS_COVARIATION, "cs.txt")
    velar = parsed[("LOC.SG", "dictionary", "any", "velar")]
    assert velar.pairs == {
        "-ní/-ný vs -ový/-ovní": (4, 5, 47, 344),
        "-it vs -ovat": (0, 1, 21, 40),
    }
    assert velar.pooled.mh == (1.96, 0.5, 7.74, 8)
    assert velar.pooled.perm_p == 0.22
    assert velar.logit is not None and velar.logit.coef == -0.282
    assert velar.both == 46
    baseline = parsed[("LOC.SG", "dictionary", "any", "baseline")]
    assert baseline.logit is None and baseline.pooled.perm_p == 0.01


CS_SPELLOUT = """\
1. Feminine set, basis wikt, all codings
------------------------------------------------------------------------
  pair  tau  coding  subset  roots  N Nonly  e  theta  tol  frag e(ret)
  -ní/-ný vs -ový/-ovní dictionary dominant all 3786 301 271 273 52.7 no 221 30
  -ka vs -yně  ge1  dominant  all  1125  0  0  0  -  -  -  0

3. Feminine set by base stratum (Wiktionary; basis wikt; dominant)
------------------------------------------------------------------------
  -it vs -ovat ge1 dominant stratum:borrowed whole 56 3 3 3 2.7 no 1 0

4b. E0-certain velars by R26 suffix target (basis wikt; dominant)
------------------------------------------------------------------------
  -ka vs -yně  ge1  dominant  base:root-final  749  50  8  8  12.8  yes  5  42

6. NC sources (aggregates only): basis all (DeriNet-only pairs
------------------------------------------------------------------------
  -ní/-ný vs -ový/-ovní dictionary dominant all 3786 1239 1195 1204 174.0 no 1031 44
"""


def test_parse_cs_spellout_skips_strata_and_nc() -> None:
    parsed = compare.parse_cs_spellout(CS_SPELLOUT, "cs.txt")
    assert [(p.section, p.pair, p.subset) for p in parsed] == [
        ("1", "-ní/-ný vs -ový/-ovní", "all"),
        ("1", "-ka vs -yně", "all"),
        ("4b", "-ka vs -yně", "base:root-final"),
    ]
    first = parsed[0]
    assert (first.roots, first.N, first.e, first.e_ret) == (3786, 301, 273, 30)
    assert first.tolerated is False and first.theta == 52.7
    assert parsed[1].theta is None and parsed[1].tolerated is None
    assert parsed[2].tolerated is True


def test_parse_ro_spellout_contrast() -> None:
    text = (
        "3. Contrast and residue (dictionary tau)\n"
        "  dorsal   non-alternating-plural roots (-uri or resisting): -i "
        "only 9, -a/-ui only 28, both 2\n"
        "  coronal  non-alternating-plural roots (-uri or resisting): -i "
        "only 27, -a/-ui only 80, both 3\n"
    )
    assert compare.parse_ro_spellout_contrast(text, "ro.txt") == {
        "dorsal": (9, 28, 2),
        "coronal": (27, 80, 3),
    }
    with pytest.raises(compare.ReportFormatError):
        compare.parse_ro_spellout_contrast("nothing here\n", "ro.txt")


CS_LOG_SUMMARY = """\
tau = ge1: allomorphs of the attested cells
  nom_pl
    -log   person   -ové                  41
    -log   thing    -y                     5
    -urg   person   -ové                   2
  loc_pl
    -log   person   -ích                   2
    -log   thing    -ích                   4

tau = ge0.5pm: allomorphs of the attested cells
  nom_pl
    -log   person   -ové                  10

LOC.PL palatalized (-lozích, -rzích) among -log/-urg nouns:
  dictionary person:None 2, person:True 107, thing:True 9, unknown:True 1
  ge1        person:True 2, thing:True 4
"""


def test_parse_cs_logtable_ge1() -> None:
    allomorphs, palatalized = compare.parse_cs_logtable_ge1(
        CS_LOG_SUMMARY, "cs.txt"
    )
    assert allomorphs[("nom_pl", "-log", "person")] == {"-ové": 41}
    assert allomorphs[("nom_pl", "-urg", "person")] == {"-ové": 2}
    assert allomorphs[("loc_pl", "-log", "thing")] == {"-ích": 4}
    assert ("nom_pl", "-log", "person") in allomorphs
    assert palatalized == {"person": 2, "thing": 4}


def _log_row(lemma: str, forms: str, endings: str) -> dict[str, str]:
    """A czech_logtable.csv row of a learned -log person noun."""
    return {
        "lemma": lemma,
        "suffixoid": "-log",
        "learned_suffixoid": "True",
        "semantic_class": "person",
        "nom_pl": forms,
        "nom_pl_allomorph": endings,
        "loc_pl_allomorph": "-ích",
        "loc_pl_palatalized": "True",
    }


def test_cs_log_counts_read_the_ending_of_each_form() -> None:
    rows = [
        _log_row("biolog", "biologové|biolozi", "-ové|-i"),
        _log_row("geolog", "geologové", "-ové"),
        _log_row("teolog", "teologi", "-i"),
    ]
    counts = compare._cs_log_counts(rows, "person")
    assert counts["nouns"] == 3
    assert counts["nom_front"] == 2
    assert counts["nom_pal"] == 1  # biolozi
    assert counts["nom_plain"] == 1  # teologi
    assert counts["nom_shielded"] == 1


def test_cs_log_counts_refuse_forms_not_aligned_with_endings() -> None:
    # The logtable lists each ending once: two -ové forms before the -i
    # form would point "-i" at the second -ové form (counted as a
    # palatalized -i plural) unless the misalignment is refused.
    row = _log_row(
        "kardiolog", "kardiologové|kardiologovė|kardiologi", "-ové|-i"
    )
    with pytest.raises(compare.ReportFormatError, match="not aligned"):
        compare._cs_log_counts([row], "person")


def test_parse_validation() -> None:
    text = (
        "Analysis: czech-palatalization (lexical_loo)\n"
        "Target-stem rows validated: 48,143\n"
        "BOOLEAN mutation match: 25,223/48,143  (52.39%)\n"
        "SURFACE-IPA match: 25,223/48,143  (52.39%)\n"
    )
    assert compare.parse_validation(text, "v.txt") == (48143, 25223, 25223)
    with pytest.raises(compare.ReportFormatError):
        compare.parse_validation("Analysis: x\n", "v.txt")


def test_report_format_errors_name_the_file() -> None:
    with pytest.raises(compare.ReportFormatError, match="ro.txt"):
        compare.parse_ro_covariation("nothing\n", "ro.txt")
    broken = RO_COVARIATION.replace("    both excluded:", "    both gone:")
    with pytest.raises(compare.ReportFormatError, match="incomplete"):
        compare.parse_ro_covariation(broken, "ro.txt")


# ---------------------------------------------------------------------------
# The README block
# ---------------------------------------------------------------------------


def test_romanian_marker_prefix_matches_the_plugin() -> None:
    assert ro_readme.BEGIN.startswith(compare.ROMANIAN_SECTION_PREFIX)


def test_merge_block_keeps_everything_else() -> None:
    romanian = f"{ro_readme.BEGIN}\nRomanian\n{ro_readme.END}\n"
    existing = "# Title\n\nProse before.\n\n" + romanian
    merged = compare.merge_block(existing, "Table v1\n")
    assert merged.startswith("# Title\n\nProse before.\n\n")
    assert merged.endswith(romanian)
    assert compare.extract_block(merged) == "Table v1\n"
    # Replace in place: text after the block (and before) is kept.
    edited = merged.replace("Prose before.", "Prose edited.") + "\nTail.\n"
    again = compare.merge_block(edited, "Table v2\n")
    assert again == edited.replace("Table v1", "Table v2")
    assert compare.merge_block(again, "Table v2\n") == again
    # The Romanian merge leaves the cross-language block alone.
    remerged = ro_readme.merge_section(again, "Romanian v2")
    assert compare.extract_block(remerged) == "Table v2\n"
    created = compare.merge_block(None, "T\n")
    assert created.startswith(compare.README_TITLE)
    with pytest.raises(ValueError):
        compare.merge_block(compare.BEGIN + "\nno end\n", "T\n")


# ---------------------------------------------------------------------------
# End to end on a copy of the tracked inputs
# ---------------------------------------------------------------------------


def copy_inputs(root: Path) -> compare.Inputs:
    inputs = compare.Inputs(root)
    for path in REAL.all():
        target = root / REAL.relative(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
    readme = PROJECT_ROOT / "analysis" / "comparison" / "README.md"
    shutil.copyfile(readme, root / "analysis" / "comparison" / "README.md")
    return inputs


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@needs_tracked
def test_run_writes_only_its_outputs_and_skips_when_current(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    copy_inputs(tmp_path)
    before = {p: digest(p) for p in tmp_path.rglob("*") if p.is_file()}
    ctx = compare.context(tmp_path)
    assert compare.run(ctx) == 0
    out = compare.outputs(tmp_path)
    for path, sha in before.items():
        if path != out.readme:
            assert digest(path) == sha, path
    allowed = (
        tmp_path / "analysis" / "comparison",
        tmp_path / "data" / "romanian_compare",
    )
    for path in tmp_path.rglob("*"):
        if path.is_file() and path not in before:
            assert any(base in path.parents for base in allowed), path
    manifest = read_manifest(ctx.manifest_path(compare.STAGE))
    assert manifest is not None and manifest.params["version"]
    assert compare.stale(tmp_path) == []
    capsys.readouterr()
    assert compare.run(ctx) == 0
    assert "up to date" in capsys.readouterr().out
    # The hand-written text and the Romanian section survive.
    text = out.readme.read_text(encoding="utf-8")
    original = (PROJECT_ROOT / "analysis/comparison/README.md").read_text(
        encoding="utf-8"
    )
    start = original.find(compare.BEGIN)
    assert text[:start] == original[:start]
    assert text.count(ro_readme.BEGIN) == 1


@needs_tracked
def test_check_and_missing_inputs(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    inputs = copy_inputs(tmp_path)
    (tmp_path / "analysis" / "comparison" / "README.md").write_text(
        "# Comparison\n", encoding="utf-8"
    )
    assert compare.main(["--root", str(tmp_path), "--check"]) == 1
    assert "stale:" in capsys.readouterr().out
    assert compare.main(["--root", str(tmp_path)]) == 0
    assert compare.main(["--root", str(tmp_path), "--check"]) == 0
    inputs.cs_spellout.unlink()
    assert compare.main(["--root", str(tmp_path), "--force"]) == 1
    assert "inputs missing" in capsys.readouterr().out


def test_context_is_confined_to_the_comparison() -> None:
    ctx = compare.context(PROJECT_ROOT)
    with pytest.raises(PermissionError):
        ctx.check_writable(PROJECT_ROOT / "analysis" / "italian" / "x.csv")
    with pytest.raises(PermissionError):
        ctx.check_writable(PROJECT_ROOT / "analysis" / "x.csv")
    ctx.check_writable(PROJECT_ROOT / "analysis" / "comparison" / "x.csv")


# ---------------------------------------------------------------------------
# The real tracked outputs
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def rows() -> list[dict[str, Any]]:
    if not HAVE_REAL:
        pytest.skip("tracked language outputs are absent")
    return compare.build_rows(PROJECT_ROOT)


def one(rows: list[dict[str, Any]], **where: Any) -> dict[str, Any]:
    found = [r for r in rows if all(r[k] == v for k, v in where.items())]
    assert len(found) == 1, (where, len(found))
    return found[0]


@needs_tracked
def test_tracked_outputs_are_current() -> None:
    assert (
        compare.stale(PROJECT_ROOT) == []
    ), "regenerate with: PYTHONPATH=src python -m lexpipe.compare"


@needs_tracked
def test_every_language_has_the_matched_rule_family(
    rows: list[dict[str, Any]],
) -> None:
    for lang in compare.LANGS:
        for item in (compare.R17A, compare.R17B, compare.R17B_NT):
            for tau in compare.TAUS:
                row = one(
                    rows,
                    section="r17",
                    lang=lang,
                    item=item,
                    target_class="dorsal",
                    tau=tau,
                    matched=True,
                )
                assert row["evidence_set"] == "E0"
                assert row["coding"] == "dominant"
                assert row["freq_basis"] == "lemma"
                assert row["N"] > 0 and row["verdict"]
        learners = [
            r
            for r in rows
            if r["section"] == "learners"
            and r["lang"] == lang
            and r["matched"]
        ]
        assert {r["V"] for r in learners} == set(compare.LEARNER_VS)
        for r in learners:
            total = sum(
                r[f"share_{name}"] for name in compare.LEARNER_OUTCOMES
            )
            assert total == pytest.approx(1.0, abs=1e-6)


@needs_tracked
def test_rows_agree_with_the_language_reports(
    rows: list[dict[str, Any]],
) -> None:
    # Italian nouns, R17b -i at >= 1 token (italian_tp_summary.txt).
    it = one(
        rows,
        section="r17",
        lang="it",
        item=compare.R17B,
        subset="nouns",
        tau="ge1",
        V="",
    )
    tp_text = (
        PROJECT_ROOT / "analysis/italian/italian_tp_summary.txt"
    ).read_text(encoding="utf-8")
    ge1 = tp_text.split(" tau = ge1\n", 1)[1].split(" tau = ", 1)[0]
    expected = (
        rf"R17b -i, nouns only\s+N_target\s+{it['N_target']}\s+N\s+"
        rf"{it['N']}\s+e_pal\s+{it['e_palatalize']}\s+e_ret\s+"
        rf"{it['e_retain']}\s"
    )
    assert re.search(expected, ge1)
    # R23: the statistics computed here equal the printed ones.
    ro = one(
        rows,
        section="covariation",
        lang="ro",
        target_class="dorsal",
        tau="dictionary",
    )
    ro_text = REAL.ro_covariation.read_text(encoding="utf-8")
    printed = (
        f"OR {ro['odds_ratio']:.3g} [{ro['or_low']:.3g}, {ro['or_high']:.3g}],"
        f" Fisher p = {ro['fisher_p']:.3g}"
    )
    assert printed in ro_text
    cs = one(
        rows,
        section="covariation",
        lang="cs",
        item="LOC.SG class (any) × -ní/-ný vs -ový/-ovní",
        target_class="velar",
        tau="dictionary",
    )
    cs_text = REAL.cs_covariation.read_text(encoding="utf-8")
    assert (
        f"OR {cs['odds_ratio']:.2f} [{cs['or_low']:.2f}, {cs['or_high']:.2f}]"
        f" Fisher p={cs['fisher_p']:.3g}"
    ) in cs_text
    # The Italian pooled Fisher p equals the Italian README's.
    itcv = one(
        rows,
        section="covariation",
        lang="it",
        item="learned suffixes, internal formations",
        tau="dictionary",
    )
    it_readme = (PROJECT_ROOT / "docs/italian/README.md").read_text(
        encoding="utf-8"
    )
    assert (
        f"| dictionary | learned | {itcv['alt_pal']} / {itcv['alt_plain']} "
        f"| {itcv['nonalt_pal']} / {itcv['nonalt_plain']} | 0 | "
        f"{itcv['fisher_p']:.3g} |"
    ) in it_readme


@needs_tracked
def test_rows_pass_the_licence_gate_and_rest_on_wiktionary(
    rows: list[dict[str, Any]],
) -> None:
    for row in rows:
        names = licensing.split_sources(row["source"])
        licensing.check_redistributable(names)
        assert all(name.startswith("wiktionary_") for name in names)
        assert row["licence_class"] == "by-sa"
        assert row["input"] and (PROJECT_ROOT / row["input"]).is_file()
    sections = [r["section"] for r in rows]
    assert sorted(set(sections), key=compare.SECTIONS.index) == list(
        compare.SECTIONS
    )


@needs_tracked
def test_csv_round_trip(tmp_path: Path, rows: list[dict[str, Any]]) -> None:
    path = tmp_path / "analysis" / "comparison" / "s.csv"
    csvio.write_tracked(path, rows, compare.FIELDS, root=tmp_path)
    assert path.read_text(encoding="utf-8") == compare.csv_text(rows)
    back = csvio.read_rows(path)
    assert len(back) == len(rows) and list(back[0]) == list(compare.FIELDS)


@needs_tracked
@pytest.mark.slow
def test_output_does_not_depend_on_the_hash_seed() -> None:
    code = (
        "import hashlib, sys; sys.path.insert(0, 'src'); "
        "from lexpipe import compare; "
        "rows = compare.build_rows(); "
        "text = compare.csv_text(rows) + compare.render_block(rows); "
        "print(hashlib.sha256(text.encode()).hexdigest())"
    )
    digests = set()
    for seed in ("0", "1"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        result = subprocess.run(
            [sys.executable, "-c", code],
            cwd=PROJECT_ROOT,
            env=env,
            capture_output=True,
            text=True,
            check=True,
        )
        digests.add(result.stdout.strip())
    assert len(digests) == 1
