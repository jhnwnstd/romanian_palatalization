"""Parse the cs-ndecl declension spec, for audit and form provenance (R11, R13).

English Wiktionary declines Czech nouns with ``{{cs-ndecl|<spec>}}``,
where the spec is a dot-separated list of indicators read by
``Module:cs-noun``. The spec records only what the editor overrode;
everything else comes from the module's defaults, so it can never be
the source of the cells (R13: cells come from wiktextract's ``forms``).
It does say which cells an editor set by hand, which is what
``form_source`` needs (R11): a cell whose slot is overridden is
``explicit``, every other cell ``module_default``.

Grammar handled (indicator examples from the pinned dump):

  - slot overrides ``<slot>[pl]<endings>``, several slots joined with
    ``+``: ``gena``, ``locu:e``, ``nompli:ové``, ``nomplové``,
    ``locplách``, ``genplou:-``, ``insplama``, ``dat+locu``. The slot
    names are nom gen dat acc voc loc ins, ``pl`` marks the plural, and
    the endings are ``:``-separated alternatives (``-`` = zero ending);
  - flags: gender/animacy (m f n an in inan), stem alternations
    (``*`` ``-*`` ``#`` ``##`` ``ě`` ``#ě``), declension types (``soft``
    ``hard`` ``mixed`` ``istem`` ``tstem`` ``nstem`` ``foreign``
    ``mostlyindecl`` ``indecl`` ``+`` adjectival), number restrictions
    (``sg`` ``pl``), ``surname`` ``persname`` ``tech``, value-carrying
    flags (``decllemma:X`` ``declgender:X`` ``mixedistem:X``, kept by
    name), and ``-velar``/``c_as_k`` respelling hints;
  - ``!``: a manual table whose every cell is written out (týden), so
    all of its cells are explicit;
  - alternants ``((A,B))`` (each alternant parsed, overrides and flags
    united) and multiword specs ``Česká<+> republika<f.sg>`` (each
    ``<…>`` parsed; :attr:`DeclSpec.multiword` is set);
  - footnotes ``[rare]``, which may contain anything, are removed first.

Cells copied by the module from another slot inherit its explicit
status (:func:`explicit_cells`): VOC.PL is always NOM.PL; ACC.PL is
NOM.PL except for masculine animates; ACC.SG of a masculine animate is
GEN.SG.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Final

CASE_SLOTS: Final[dict[str, str]] = {
    "nom": "NOM",
    "gen": "GEN",
    "dat": "DAT",
    "acc": "ACC",
    "voc": "VOC",
    "loc": "LOC",
    "ins": "INS",
}
_SLOT: Final[str] = "(?:nom|gen|dat|acc|voc|loc|ins)(?:pl)?"
_OVERRIDE_RE: Final[re.Pattern[str]] = re.compile(
    rf"^(?P<slots>(?:{_SLOT}\+)*{_SLOT})(?P<endings>.+)$"
)
_VALUE_FLAGS: Final[tuple[str, ...]] = (
    "decllemma",
    "declgender",
    "mixedistem",
)
MANUAL_FLAG: Final[str] = "!"


@dataclass(frozen=True, slots=True)
class DeclSpec:
    """A parsed cs-ndecl spec (see the module docstring).

    ``overrides`` holds ``(case, number)`` pairs such as
    ``("NOM", "PL")``; ``flags`` the other indicators, sorted.
    """

    raw: str
    overrides: frozenset[tuple[str, str]] = field(default_factory=frozenset)
    flags: tuple[str, ...] = ()
    manual: bool = False
    multiword: bool = False
    alternants: int = 1

    def has_flag(self, flag: str) -> bool:
        return flag in self.flags


def strip_footnotes(spec: str) -> str:
    """Remove ``[…]`` footnotes (which may hold dots, commas, tags)."""
    out: list[str] = []
    depth = 0
    for ch in spec:
        if ch == "[":
            depth += 1
            continue
        if ch == "]" and depth:
            depth -= 1
            continue
        if depth == 0:
            out.append(ch)
    return "".join(out)


def _split_top_level(text: str) -> list[str]:
    """Split on commas outside ``<…>`` and ``((…))``."""
    parts: list[str] = []
    depth = 0
    current: list[str] = []
    for ch in text:
        if ch in "<(":
            depth += 1
        elif ch in ">)" and depth:
            depth -= 1
        if ch == "," and depth == 0:
            parts.append("".join(current))
            current = []
        else:
            current.append(ch)
    parts.append("".join(current))
    return [part for part in (p.strip() for p in parts) if part]


def _angle_specs(text: str) -> list[str]:
    """Contents of every ``<…>`` group (outermost) in *text*."""
    found: list[str] = []
    depth = 0
    start = 0
    for index, ch in enumerate(text):
        if ch == "<":
            if depth == 0:
                start = index + 1
            depth += 1
        elif ch == ">" and depth:
            depth -= 1
            if depth == 0:
                found.append(text[start:index])
    return found


def _indicators(spec: str) -> list[str]:
    return [
        token
        for token in ("".join(part.split()) for part in spec.split("."))
        if token
    ]


def _overridden_slots(indicator: str) -> list[tuple[str, str]] | None:
    match = _OVERRIDE_RE.match(indicator)
    if match is None:
        return None
    slots: list[tuple[str, str]] = []
    for slot in match.group("slots").split("+"):
        number = "PL" if slot.endswith("pl") else "SG"
        slots.append((CASE_SLOTS[slot[:3]], number))
    return slots


def _flag_name(indicator: str) -> str:
    for name in _VALUE_FLAGS:
        if indicator.startswith(name + ":"):
            return name
    return indicator


def parse_decl_spec(raw: str) -> DeclSpec:
    """Parse one cs-ndecl spec string (see the module docstring)."""
    text = strip_footnotes(raw).strip()
    alternant_text = text
    has_alternants = text.startswith("((") and text.endswith("))")
    if has_alternants:
        alternant_text = text[2:-2]
    alternants = _split_top_level(alternant_text) if has_alternants else [text]
    multiword = False
    overrides: set[tuple[str, str]] = set()
    flags: set[str] = set()
    for alternant in alternants or [""]:
        groups = _angle_specs(alternant)
        if groups:
            outside = re.sub(r"<[^<>]*>", "", alternant).strip()
            # "uši<n.pl…>" names another lemma; "a<+> b<f>" is multiword.
            multiword = multiword or len(groups) > 1 or " " in outside
            if outside and len(groups) == 1 and " " not in outside:
                flags.add("altlemma")
        for spec in groups or [alternant]:
            for indicator in _indicators(spec):
                slots = _overridden_slots(indicator)
                if slots is not None:
                    overrides.update(slots)
                else:
                    flags.add(_flag_name(indicator))
    return DeclSpec(
        raw=raw,
        overrides=frozenset(overrides),
        flags=tuple(sorted(flags)),
        manual=MANUAL_FLAG in flags,
        multiword=multiword,
        alternants=len(alternants),
    )


def explicit_cells(
    spec: DeclSpec, *, masculine_animate: bool
) -> frozenset[tuple[str, str]]:
    """Cells whose forms an editor set: overridden slots plus the cells
    the module copies from them (see the module docstring)."""
    cells = set(spec.overrides)
    if ("NOM", "PL") in cells:
        cells.add(("VOC", "PL"))
        if not masculine_animate:
            cells.add(("ACC", "PL"))
    if masculine_animate and ("GEN", "SG") in cells:
        cells.add(("ACC", "SG"))
    return frozenset(cells)
