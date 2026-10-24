"""Port of the default-inflection functions of ``Module:it-common``.

English Wiktionary generates most Italian plurals instead of storing
them: an ``{{it-noun}}`` or ``{{it-adj}}`` head without a plural
argument gets its plural from ``make_plural`` in ``Module:it-common``
when the page is rendered, and kaikki records the rendered result. Those
defaults encode the very generalisation under study (-ico → -ici,
-ologo → -ologi, other -co/-go → -chi/-ghi), so design rule R12 keeps
them apart from hand-entered plurals and requires a port of the module
plus a drift test (the ``lexicon`` stage) that checks the port against
kaikki's output for every noun without a plural argument.

Source
------
Ported from ``Module:it-common`` on English Wiktionary,
https://en.wiktionary.org/w/index.php?title=Module:it-common&oldid=79032535
— revision **79032535** (parent 78848562) of 2024-04-24T16:59:39Z,
11,812 bytes, content SHA-1 ``64837177c162e2b2506312fe1af017ea14fcd9af``
(checked against the MediaWiki API's revision record when the source
was fetched on 2026-10-23; the copy is kept, gitignored, in
``data/italian/upstream/``). The module is CC BY-SA 4.0 Wiktionary
content; this port re-implements its logic.

What is ported
--------------
:func:`make_plural`, :func:`make_feminine` and :func:`make_masculine`,
rule for rule and in the module's order:

    -io → -i;  -ologo → -ologi;  -[ia]co → -[ia]ci;  -[cg]o → -[cg]hi;
    -o → -i;   -[cg]a → -[cg]hi (m) / -[cg]he (f);  -logia → -logie (f);
    -[cg]ia → -[cg]i (m) / -[cg]e (f);  -a → -i (m) / -e (f);  -e → -i;
    anything else: no default (None).

The ``cap*``/``cap*+`` specials for capo- compounds are ported. Like the
module, the gender-specific branches raise for genders ``mf``,
``mfbysense`` and ``?`` (:class:`DefaultPluralError`), as does -logia
with a gender other than ``f``.

Not ported: multiword and hyphenated terms, which the module hands to
``Module:romance utilities`` (``handle_multiword``); :func:`make_plural`
and :func:`make_feminine` raise :class:`NotPortedError` for them. The
module's ``make_masculine`` computes that multiword result and then
discards it (it falls through to the single-word rules on the whole
term); the port reproduces that behaviour.

Headword level (reconstructed, not ported)
------------------------------------------
How ``Module:it-headword`` calls these functions was not fetched. The
cell-level defaults below are reconstructed from kaikki's rendered
heads and are validated by the same drift test:

  - a noun without plural arguments: one ``pl`` cell, ``make_plural(
    lemma, gender)``; for gender ``mf``/``mfbysense``, ``m.pl`` and
    ``f.pl`` from the masculine and feminine rules when they differ
    (*collega*: colleghi/colleghe), else one ``pl`` cell (*cantante*);
  - an adjective: ``m.pl`` = make_plural(lemma, "m") and ``f.pl`` =
    make_plural(feminine, "f") for each feminine singular (the ``f``
    argument, else make_feminine(lemma)); one ``pl`` cell when the two
    coincide (*mobile*, *verde*); no default when neither exists. A
    head with ``noforms`` shows no forms at all.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from typing import Final

MODULE_TITLE: Final[str] = "Module:it-common"
MODULE_REVISION: Final[int] = 79032535
MODULE_PARENT_REVISION: Final[int] = 78848562
MODULE_TIMESTAMP: Final[str] = "2024-04-24T16:59:39Z"
MODULE_SHA1: Final[str] = "64837177c162e2b2506312fe1af017ea14fcd9af"
MODULE_SIZE: Final[int] = 11812
MODULE_URL: Final[str] = (
    "https://en.wiktionary.org/w/index.php?title=Module:it-common"
    f"&oldid={MODULE_REVISION}"
)

CAP_SPECIALS: Final[frozenset[str]] = frozenset({"cap*", "cap*+"})
_GENDERS_WITHOUT_DEFAULT: Final[frozenset[str]] = frozenset(
    {"mf", "mfbysense", "?"}
)
_IA_CO_RE: Final[re.Pattern[str]] = re.compile(r"[ia]co\Z")
_CG_O_RE: Final[re.Pattern[str]] = re.compile(r"[cg]o\Z")
_CG_A_RE: Final[re.Pattern[str]] = re.compile(r"[cg]a\Z")
_CG_IA_RE: Final[re.Pattern[str]] = re.compile(r"[cg]ia\Z")


class DefaultPluralError(ValueError):
    """The module raises an error: no default for this term and gender."""


class NotPortedError(NotImplementedError):
    """The term needs ``Module:romance utilities`` (multiword handling)."""


def is_multiword(term: str) -> bool:
    """True if the module treats *term* as multiword (space or hyphen)."""
    return " " in term or "-" in term


def _check_no_mf(term: str, gender: str, special: str | None) -> None:
    if gender in _GENDERS_WITHOUT_DEFAULT:
        using = f" using special={special}" if special else ""
        raise DefaultPluralError(
            f"With gender={gender}, unable to pluralize term '{term}'"
            f"{using} because its plural is gender-specific"
        )


def _capo(term: str, gender: str, plspec: str) -> str:
    _check_no_mf(term, gender, plspec)
    if not term.startswith("capo"):
        raise DefaultPluralError(
            f"With special={plspec}, term '{term}' must begin with capo-"
        )
    return "capi" + term[4:] if gender == "m" else term


def make_plural(
    term: str, gender: str, special: str | None = None
) -> str | None:
    """Default plural of *term* (``ex.make_plural``); None if none.

    *special* may be ``cap*`` or ``cap*+``; any other special, and any
    multiword term, raises :class:`NotPortedError`.
    """
    plspec = special if special in CAP_SPECIALS else None
    rest = None if plspec else special
    if rest is not None or is_multiword(term):
        raise NotPortedError(
            f"multiword handling (Module:romance utilities) is not ported: "
            f"{term!r}, special={special!r}"
        )
    if plspec is not None:
        term = _capo(term, gender, plspec)
        if plspec == "cap*":
            return term
    if term.endswith("io"):
        return term[:-2] + "i"
    if term.endswith("ologo"):
        return term[:-1] + "i"
    if _IA_CO_RE.search(term):
        return term[:-1] + "i"
    if _CG_O_RE.search(term):
        return term[:-1] + "hi"
    if term.endswith("o"):
        return term[:-1] + "i"
    if _CG_A_RE.search(term):
        _check_no_mf(term, gender, plspec)
        return term[:-1] + ("hi" if gender == "m" else "he")
    if term.endswith("logia"):
        if gender != "f":
            raise DefaultPluralError(
                f"Term '{term}' ending in -logia should have gender=f if it "
                "is using the default plural"
            )
        return term[:-1] + "e"
    if _CG_IA_RE.search(term):
        _check_no_mf(term, gender, plspec)
        return term[:-2] + ("i" if gender == "m" else "e")
    if term.endswith("a"):
        _check_no_mf(term, gender, plspec)
        return term[:-1] + ("i" if gender == "m" else "e")
    if term.endswith("e"):
        return term[:-1] + "i"
    return None


def make_feminine(term: str) -> str:
    """Default feminine (``ex.make_feminine``): -o → -a, -tore →
    -trice, -one → -ona; otherwise unchanged."""
    if is_multiword(term):
        raise NotPortedError(
            f"multiword handling (Module:romance utilities) is not ported: "
            f"{term!r}"
        )
    if term.endswith("o"):
        return term[:-1] + "a"
    if term.endswith("tore"):
        return term[: -len("tore")] + "trice"
    if term.endswith("one"):
        return term[: -len("one")] + "ona"
    return term


def make_masculine(term: str) -> str:
    """Default masculine (``ex.make_masculine``): -a → -o, -trice →
    -tore; otherwise unchanged. Multiword terms get the same rules
    applied to the whole term, as in the module."""
    if term.endswith("a"):
        return term[:-1] + "o"
    if term.endswith("trice"):
        return term[: -len("trice")] + "tore"
    return term


def noun_default_cells(lemma: str, gender: str) -> list[tuple[str, str]]:
    """Default ``(cell, form)`` pairs of a noun without plural arguments.

    Reconstructed headword behaviour (module docstring). Raises like
    :func:`make_plural`; a term with no default gives ``[]``.
    """
    if gender in ("mf", "mfbysense"):
        masculine = make_plural(lemma, "m")
        feminine = make_plural(lemma, "f")
        if masculine == feminine:
            return [("pl", masculine)] if masculine else []
        cells = [("m.pl", masculine), ("f.pl", feminine)]
        return [(cell, form) for cell, form in cells if form]
    form = make_plural(lemma, gender)
    return [("pl", form)] if form else []


def adjective_default_cells(
    lemma: str, feminines: Sequence[str] = (), *, feminine_only: bool = False
) -> list[tuple[str, str]]:
    """Default ``(cell, form)`` pairs of an ``{{it-adj}}`` without
    plural arguments (reconstructed; module docstring).

    *feminines* are the feminine singulars stated by the head's ``f``
    argument (default: :func:`make_feminine` of the lemma). When the
    masculine and feminine plurals coincide the adjective has one
    ``pl`` cell (*mobile* → mobili), otherwise ``m.pl`` and ``f.pl``
    (*medico* → medici/mediche; *editore* → editori/editrici). A
    feminine-only adjective (``fonly``) has only ``f.pl``.
    """
    if is_multiword(lemma):
        raise NotPortedError(
            f"multiword handling (Module:romance utilities) is not ported: "
            f"{lemma!r}"
        )
    if feminine_only:
        only = make_plural(lemma, "f")
        return [("f.pl", only)] if only else []
    singulars = list(feminines) or [make_feminine(lemma)]
    masculine = make_plural(lemma, "m")
    feminine = [make_plural(form, "f") for form in singulars]
    if masculine is None and not any(feminine):
        return []
    if feminine == [masculine]:
        return [("pl", masculine)] if masculine else []
    cells = [("m.pl", masculine)] + [("f.pl", form) for form in feminine]
    return [(cell, form) for cell, form in cells if form]
