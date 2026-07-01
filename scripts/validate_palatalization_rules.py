#!/usr/bin/env python3
"""Run the palatalization validator and persist its report + CSV.

This script is a thin example of how to CALL the Python API. It is
NOT a CLI — no argparse, no flags, no ``sys.argv``. If you want to
tweak what runs (limit rows, disable normalisation, trace lemmas,
run ordering / perturbation searches), edit the constants at the top
or — better — open a REPL / notebook and call the functions in
``phonology.validation`` directly. See ``README`` and
``tests/test_validation_api.py`` for usage patterns.

Business logic lives in ``src/phonology/``:

  - :mod:`phonology.analyses.romanian_palatalization` — the analysis
    (patches, underspec, rules, oracle) and :func:`build_profile`.
  - :mod:`phonology.lexicon` — :func:`iter_lexicon_rows` +
    :class:`LexRow`.
  - :mod:`phonology.validation` — :func:`run_validation`,
    :func:`trace_lemmas`, :func:`render_report`,
    :func:`write_mismatch_csv`.
  - :mod:`phonology.diagnostics.*` — the four diagnostic modules
    (compare, distance, explain, ordering, perturb, categorise).
"""

from __future__ import annotations

import sys
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PROJECT_ROOT / "src"))

from phonology.analyses.romanian_palatalization import (  # noqa: E402
    build_profile,
)
from phonology.lexicon import iter_lexicon_rows  # noqa: E402
from phonology.validation import (  # noqa: E402
    render_report,
    run_validation,
    write_mismatch_csv,
)

CSV_PATH = _PROJECT_ROOT / "data" / "romanian_lexicon_with_freq.csv"
REPORT_PATH = _PROJECT_ROOT / "analysis" / "palatalization_rule_validation.txt"
MISMATCH_CSV = _PROJECT_ROOT / "analysis" / "palatalization_mismatches.csv"


def main() -> None:
    profile = build_profile()
    rows = iter_lexicon_rows(CSV_PATH, target_stems=profile.target_stems)
    report = run_validation(profile, rows)

    text = render_report(report, include_distance_summary=True)
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(text, encoding="utf-8")
    print(text, end="")
    print(f"\nReport persisted to {REPORT_PATH}")

    n = write_mismatch_csv(report, MISMATCH_CSV)
    print(f"Mismatch CSV: {n:,} rows → {MISMATCH_CSV}")


if __name__ == "__main__":
    main()
