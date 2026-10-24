"""Hand adjudications of Czech cells (curated/czech/cell_adjudications.csv).

The automated QC (:mod:`.qc`) compares Wiktionary cells with MorfFlex
and lists the disagreements in ``analysis/czech/czech_manual_worklist.csv``
with a link to the IJP entry. A person looks each one up **by hand**
(IJP's robots.txt disallows ``/?slovo=`` and copying its content needs
ÚJČ consent, so nothing here fetches it) and records a coded decision
in the curated file, never IJP text (R9). One row per decided cell,
keyed by ``(lexeme_id, cell)`` (:func:`lexpipe.curated.load_curated`
with a composite key):

    lexeme_id         the lexeme, as in czech_lexicon_full.csv
    cell              one of :data:`CELLS`
    confirmed         the Wiktionary forms of that cell that the source
                      confirms as standard, ``|``-joined in Wiktionary
                      order; ``-`` = none of them
    citation          e.g. "IJP, prirucka.ujc.cas.cz, 2026-11-02"
    note              free text (no copied entry text)
    guard_lemma, guard_gender, guard_animacy, guard_forms
                      what the curator saw: the lexeme's lemma and
                      readings and the cell's Wiktionary forms

The guards follow :mod:`lexpipe.curated`: if a later dump changes the
lexeme (another homograph takes the id, the cell's forms change), the
features stage stops with :class:`lexpipe.curated.CuratedGuardError`
instead of applying the decision to the wrong thing. The file is
tracked and read-only to the pipeline; it starts as a header-only
template. (Until 2026-10-23 the template was wide, one row per lexeme
with ``confirmed_<CELL>`` and ``guard_forms_<CELL>`` columns, because
curated keys could not be composite.)

Applying a decision (:func:`curated_status`): forms of a decided cell
are ``confirmed`` or ``rejected``; rejected forms stay in the cells
table, flagged, and are left out of the lexeme summaries.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, Final

from lexpipe import csvio, curated

CELLS: Final[tuple[str, ...]] = (
    "NOM.PL",
    "VOC.SG",
    "LOC.SG",
    "DAT.SG",
    "LOC.PL",
    "NOM.PL.M.AN",
)
NONE_CONFIRMED: Final[str] = "-"
CONFIRMED: Final[str] = "confirmed"
REJECTED: Final[str] = "rejected"
KEY_FIELDS: Final[tuple[str, str]] = ("lexeme_id", "cell")
FIELDS: Final[tuple[str, ...]] = (
    *KEY_FIELDS,
    "confirmed",
    "citation",
    "note",
    "guard_lemma",
    "guard_gender",
    "guard_animacy",
    "guard_forms",
)

Key = tuple[str, ...]


def guard_fields(
    lemma: str,
    gender: Sequence[str],
    animacy: Sequence[str],
    forms: Sequence[str],
) -> dict[str, Any]:
    """The current values the guards of one lexeme cell are compared
    with."""
    return {
        "lemma": lemma,
        "gender": list(gender),
        "animacy": list(animacy),
        "forms": list(forms),
    }


def load(path: Path) -> dict[Key, curated.CuratedRow]:
    """Curated rows by ``(lexeme_id, cell)`` (empty when the file is
    absent or header-only).

    Raises ValueError when the header differs from :data:`FIELDS`, a
    cell is not one of :data:`CELLS`, or a key repeats.
    """
    if not path.is_file():
        return {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        header = fh.readline().rstrip("\r\n").split(",")
    if tuple(header) != FIELDS:
        raise ValueError(
            f"{path}: header {header} differs from the expected {list(FIELDS)}"
        )
    rows = curated.load_curated(path, key_fields=KEY_FIELDS)
    unknown = sorted({key[1] for key in rows if key[1] not in CELLS})
    if unknown:
        raise ValueError(f"{path}: unknown cells {unknown}; known: {CELLS}")
    return rows


def decision(row: curated.CuratedRow) -> list[str] | None:
    """The confirmed forms of one row (None: no decision; [] for ``-``)."""
    value = row.values.get("confirmed", "").strip()
    if not value:
        return None
    if value == NONE_CONFIRMED:
        return []
    return csvio.split_list(value)


def decisions(
    rows: Mapping[Key, curated.CuratedRow],
) -> dict[str, dict[str, list[str] | None]]:
    """``lexeme_id -> cell -> confirmed forms`` (None: no decision)."""
    out: dict[str, dict[str, list[str] | None]] = {}
    for (lexeme_id, cell), row in sorted(rows.items()):
        out.setdefault(lexeme_id, {})[cell] = decision(row)
    return out


def validate(row: curated.CuratedRow) -> None:
    """Confirmed forms must be among the guarded Wiktionary forms."""
    confirmed = decision(row)
    if not confirmed:
        return
    seen = csvio.split_list(row.guards.get("forms", ""))
    unknown = [form for form in confirmed if form not in seen]
    if unknown:
        lexeme_id, cell = row.key
        raise ValueError(
            f"{lexeme_id}: confirmed {cell} forms {unknown} are not among "
            f"the guarded forms {seen}"
        )


def curated_status(
    decided: Mapping[str, list[str] | None], cell: str, form: str
) -> str:
    """``confirmed``, ``rejected`` or "" for one cell form."""
    confirmed = decided.get(cell)
    if confirmed is None:
        return ""
    return CONFIRMED if form in confirmed else REJECTED
