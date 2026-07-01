#!/usr/bin/env python3
"""Regenerate ipa_normalized_lemma / ipa_normalized_pl in place.

Re-runs :func:`~romanian_processor_lib.ensure_ipa_fields` on every row
of the frequency-filtered lexicon so the normalized IPA columns stay
in sync with the current ``to_ipa`` implementation (notably the
consonant-geminate collapse: doubled consonants in loanword
orthography like "addenda", "kilowatt", "allemandă" now surface as
single consonants in the IPA, matching Wiktionary's independent
phonetic transcriptions).

Preserves every other column verbatim — the raw Wiktionary columns
(``ipa_raw_lemma`` / ``ipa_raw_pl``) are the ground-truth phonetic
source and stay untouched. Only the *derived* normalized columns get
refreshed. When a raw column is populated, the normalized column
picks it up unchanged; when the raw column is empty, the normalized
column comes from ``to_ipa(orthography)`` — that's where the
geminate collapse takes effect.

Idempotent: running twice is a no-op after the first run.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from romanian_processor_lib import (  # noqa: E402
    ensure_ipa_fields,
    set_ipa_normalizer,
    tweak_nominal_ipa,
)
from wiktionary_normalizer import normalize_ipa  # noqa: E402

CSV_PATH = PROJECT_ROOT / "data" / "romanian_lexicon_with_freq.csv"


def main() -> None:
    set_ipa_normalizer(normalize_ipa)

    with CSV_PATH.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)

    if fieldnames is None:
        raise RuntimeError(f"empty CSV at {CSV_PATH}")

    updated_lemma = 0
    updated_pl = 0
    for row in rows:
        old_lemma = row.get("ipa_normalized_lemma", "")
        old_pl = row.get("ipa_normalized_pl", "")
        # Clear so ensure_ipa_fields re-derives.
        row["ipa_normalized_lemma"] = ""
        row["ipa_normalized_pl"] = ""
        ensure_ipa_fields(
            row,
            orth_key="lemma",
            raw_key="ipa_raw_lemma",
            norm_key="ipa_normalized_lemma",
            tweak_fn=tweak_nominal_ipa,
        )
        ensure_ipa_fields(
            row,
            orth_key="plural",
            raw_key="ipa_raw_pl",
            norm_key="ipa_normalized_pl",
            tweak_fn=tweak_nominal_ipa,
        )
        if row["ipa_normalized_lemma"] != old_lemma:
            updated_lemma += 1
        if row["ipa_normalized_pl"] != old_pl:
            updated_pl += 1

    with CSV_PATH.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Rows: {len(rows):,}")
    print(f"ipa_normalized_lemma updated: {updated_lemma:,}")
    print(f"ipa_normalized_pl updated:    {updated_pl:,}")
    print(f"CSV written back to {CSV_PATH}")


if __name__ == "__main__":
    main()
