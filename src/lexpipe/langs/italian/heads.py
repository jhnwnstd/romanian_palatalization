"""Grammar of the arguments of en.wiktionary's Italian headword templates.

``{{it-noun}}`` and ``{{it-adj}}`` (``Module:it-headword``) take their
plurals, genders and counterparts as comma-separated lists whose items
carry inline modifiers in angle brackets. kaikki keeps the arguments
verbatim in ``head_templates[*].args``; this module parses them. It is
the input side of design rule R12: plural provenance is read from these
arguments, never from kaikki's ``forms`` (which lose plurals, e.g.
*braccio*, and invent some, e.g. the "form" ``invariable`` of
*capotribù*).

Lists and items
---------------
:func:`split_items` splits a value at commas that are not inside angle
brackets, so ``braccia<g:f><q:in most meanings 1, 2, 4>,bracci<g:m>``
has two items. :func:`parse_item` separates an item's text from its
modifiers ``<key:value>``; a value may itself contain balanced brackets
(``<l:<<proscribed>><ref name=dop/>>``) and commas. An item whose
brackets do not balance, or whose modifier has no ``key:``, is kept
with ``ArgItem.error`` set (``unbalanced`` / ``bad_modifier``) rather
than guessed at: one en.wiktionary head (*mangiare*, the noun) ends in
``>>``.

Codes (``ArgItem.code``) are the item texts that are not word forms
(:data:`CODES`):

    +       the module default for the slot
    #       invariable: the plural equals the lemma (``#s`` is lemma + s)
    -       uncountable (followed by plurals: usually uncountable)
    ~       countable and uncountable
    ?       plural unknown
    !       plural unattested
    ++      plural of every part of a multiword/hyphenated lemma
    +first, +last, +each, …
            multiword specials of Module:romance utilities (which part
            to inflect); recognised, not resolved
    cap*    capo- compound: capo- → capi-, rest unchanged
    cap*+   capo- compound: capo- → capi-, rest pluralized

Modifiers recognised (others are kept as written): ``q``/``qq``
(qualifiers before/after), ``l``/``ll`` (labels), ``g`` (gender of the
form), ``ref`` and ``name`` (references, never copied into qualifiers),
``t``/``gloss``, ``pos``, ``id``, ``alt``, ``tr``, ``lit``.

Genders
-------
:func:`parse_gender` reads argument 1 of ``{{it-noun}}``: a comma list
of codes such as ``m``, ``f``, ``mf``, ``mfbysense``, ``n``, ``?``,
plurale-tantum ``m-p``/``f-p``, each with optional modifiers
(``f,m<l:colloquial,or,literary><ref:...>`` for *eco*).

Slots
-----
:func:`plural_slots` lists the plural-bearing arguments of a head in a
fixed order (R5): for ``it-noun`` ``2``, ``pl2``, ``pl3``… (cell
``pl``), then ``mpl``, ``mpl2``… (``m.pl``) and ``fpl``, ``fpl2``…
(``f.pl``); for ``it-adj`` ``mpl``…, ``fpl``…, ``pl``…. Items keep
their source order within a slot.
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Final

CODES: Final[frozenset[str]] = frozenset(
    {"+", "#", "-", "~", "?", "!", "++", "cap*", "cap*+"}
)
LEMMA_PLACEHOLDER: Final[str] = "#"
MULTIWORD_SPECIAL_RE: Final[re.Pattern[str]] = re.compile(r"\+[a-z]+")

NOUN_TEMPLATE: Final[str] = "it-noun"
ADJ_TEMPLATE: Final[str] = "it-adj"
HEAD_TEMPLATES: Final[tuple[str, str]] = (NOUN_TEMPLATE, ADJ_TEMPLATE)

# Cells of R12.
CELL_PL: Final[str] = "pl"
CELL_MPL: Final[str] = "m.pl"
CELL_FPL: Final[str] = "f.pl"
CELLS: Final[tuple[str, str, str]] = (CELL_PL, CELL_MPL, CELL_FPL)

ERROR_UNBALANCED: Final[str] = "unbalanced"
ERROR_BAD_MODIFIER: Final[str] = "bad_modifier"

QUALIFIER_KEYS: Final[tuple[str, ...]] = ("q", "qq", "l", "ll")
LABEL_KEYS: Final[tuple[str, ...]] = ("l", "ll")
REFERENCE_KEYS: Final[frozenset[str]] = frozenset({"ref", "name"})

# Words inside labels and qualifiers that are recognised as usage labels
# (lowercase). "or"/"and" are list connectors in Wiktionary labels.
KNOWN_LABELS: Final[frozenset[str]] = frozenset(
    {
        "archaic",
        "obsolete",
        "old italian",
        "literary",
        "poetic",
        "dialectal",
        "regional",
        "proscribed",
        "nonstandard",
        "colloquial",
        "dated",
        "rare",
        "uncommon",
        "less common",
        "very rare",
        "tuscan",
        "jocular",
        "informal",
        "careful style",
        "figuratively",
        "collective",
    }
)
_LABEL_CONNECTORS: Final[frozenset[str]] = frozenset({"or", "and", "_", ""})

_MODIFIER_KEY_RE: Final[re.Pattern[str]] = re.compile(
    r"[A-Za-z_][A-Za-z0-9_]*"
)
_REF_TAG_RE: Final[re.Pattern[str]] = re.compile(r"<ref\b[^<>]*/?>")
_SLOT_RE: Final[re.Pattern[str]] = re.compile(r"(2|pl|mpl|fpl)(\d*)")
_QUALIFIER_SPLIT_RE: Final[re.Pattern[str]] = re.compile(
    r"\s*,\s*|\s+or\s+|\s+and\s+"
)


@dataclass(frozen=True, slots=True)
class Modifier:
    """One inline modifier ``<key:value>`` (value as written)."""

    key: str
    value: str


@dataclass(frozen=True, slots=True)
class ArgItem:
    """One list item: its text, modifiers and parse problems.

    ``text`` is the item without modifiers, stripped; ``raw`` is the
    item exactly as it appears in the argument; ``error`` is "" or one
    of ``unbalanced`` / ``bad_modifier``.
    """

    text: str
    modifiers: tuple[Modifier, ...]
    raw: str
    error: str = ""

    @property
    def code(self) -> str:
        """The item's code (:data:`CODES` or a multiword special such as
        ``+first``), or "" for a word form."""
        if self.text in CODES or MULTIWORD_SPECIAL_RE.fullmatch(self.text):
            return self.text
        return ""

    def values(self, *keys: str) -> tuple[str, ...]:
        """Values of the modifiers with one of *keys*, in source order."""
        return tuple(m.value for m in self.modifiers if m.key in keys)

    @property
    def gender(self) -> str:
        """The first ``<g:>`` value ("" if none)."""
        found = self.values("g")
        return found[0].strip() if found else ""

    @property
    def qualifiers(self) -> tuple[str, ...]:
        """``key:value`` of every q/qq/l/ll modifier, references removed.

        ``|`` inside a value is written ``/`` so the result can be
        pipe-joined (:mod:`lexpipe.csvio`).
        """
        out: list[str] = []
        for modifier in self.modifiers:
            if modifier.key not in QUALIFIER_KEYS:
                continue
            text = clean_label_text(modifier.value)
            if text:
                out.append(f"{modifier.key}:{text.replace('|', '/')}")
        return tuple(out)

    @property
    def labels(self) -> tuple[str, ...]:
        """Normalised usage labels (lowercase, sorted, unique).

        Every label of ``<l:>``/``<ll:>`` counts; from ``<q:>``/``<qq:>``
        only the words in :data:`KNOWN_LABELS` do (qualifiers are free
        text such as "with anatomic or figurative meanings").
        """
        found: set[str] = set()
        for modifier in self.modifiers:
            if modifier.key in LABEL_KEYS:
                found.update(split_labels(modifier.value))
            elif modifier.key in ("q", "qq"):
                found.update(
                    label
                    for label in split_labels(modifier.value, qualifier=True)
                    if label in KNOWN_LABELS
                )
        return tuple(sorted(found))


@dataclass(frozen=True, slots=True)
class GenderCode:
    """One gender code of ``{{it-noun}}`` argument 1."""

    code: str
    base: str
    plural_only: bool
    labels: tuple[str, ...]
    qualifiers: tuple[str, ...]
    error: str = ""


@dataclass(frozen=True, slots=True)
class Slot:
    """A plural-bearing head argument: its key, cell and parsed items."""

    key: str
    cell: str
    items: tuple[ArgItem, ...]


def split_items(value: str) -> list[str]:
    """Split *value* at commas outside angle brackets; strip each item.

    Empty items are dropped. A stray ``>`` does not make the depth
    negative.
    """
    items: list[str] = []
    depth = 0
    current: list[str] = []
    for ch in value:
        if ch == "<":
            depth += 1
        elif ch == ">":
            depth = max(0, depth - 1)
        if ch == "," and depth == 0:
            items.append("".join(current).strip())
            current = []
            continue
        current.append(ch)
    items.append("".join(current).strip())
    return [item for item in items if item]


def _matching_close(text: str, start: int) -> int:
    """Index of the ``>`` closing the ``<`` at *start*, or -1."""
    depth = 0
    for index in range(start, len(text)):
        ch = text[index]
        if ch == "<":
            depth += 1
        elif ch == ">":
            depth -= 1
            if depth == 0:
                return index
    return -1


def parse_item(raw: str) -> ArgItem:
    """Parse one list item (see the module docstring)."""
    item = raw.strip()
    cut = item.find("<")
    if cut < 0:
        error = ERROR_UNBALANCED if ">" in item else ""
        return ArgItem(item.replace(">", "").strip(), (), raw, error)
    text = item[:cut].strip()
    modifiers: list[Modifier] = []
    error = ""
    position = cut
    while position < len(item):
        if item[position].isspace():
            position += 1
            continue
        if item[position] != "<":
            error = ERROR_UNBALANCED
            break
        close = _matching_close(item, position)
        if close < 0:
            error = ERROR_UNBALANCED
            break
        body = item[position + 1 : close]
        key, sep, value = body.partition(":")
        if not sep or not _MODIFIER_KEY_RE.fullmatch(key):
            error = ERROR_BAD_MODIFIER
        else:
            modifiers.append(Modifier(key, value))
        position = close + 1
    return ArgItem(text, tuple(modifiers), raw, error)


def parse_list(value: str) -> tuple[ArgItem, ...]:
    """:func:`parse_item` of every item of :func:`split_items`."""
    return tuple(parse_item(item) for item in split_items(value))


def clean_label_text(value: str) -> str:
    """A modifier value without reference tags and stray brackets."""
    text = _REF_TAG_RE.sub("", value)
    text = text.replace("<", "").replace(">", "")
    return " ".join(text.split())


def split_labels(value: str, *, qualifier: bool = False) -> list[str]:
    """Lowercase labels of a label (or qualifier) value, in order.

    Labels are comma-separated; Wiktionary writes list connectors as
    separate labels (``chiefly,archaic,or,idiomatic``), which are
    dropped. Qualifiers are split on commas and on " or " / " and ".
    """
    text = clean_label_text(value).lower()
    parts = (
        _QUALIFIER_SPLIT_RE.split(text)
        if qualifier
        else [part.strip() for part in text.split(",")]
    )
    return [part for part in parts if part not in _LABEL_CONNECTORS]


def parse_gender(value: str) -> tuple[GenderCode, ...]:
    """Parse ``{{it-noun}}`` argument 1 into its gender codes."""
    out: list[GenderCode] = []
    for item in parse_list(value):
        parts = item.text.split("-")
        base = parts[0]
        out.append(
            GenderCode(
                code=item.text,
                base=base,
                plural_only="p" in parts[1:],
                labels=item.labels,
                qualifiers=item.qualifiers,
                error=item.error,
            )
        )
    return tuple(out)


def head_template(entry: Mapping[str, object]) -> tuple[str, dict[str, str]]:
    """``(name, args)`` of the entry's first it-noun/it-adj head.

    ``("", {})`` when the entry has neither template.
    """
    templates = entry.get("head_templates")
    if not isinstance(templates, Sequence) or isinstance(templates, str):
        return "", {}
    for template in templates:
        if not isinstance(template, Mapping):
            continue
        name = template.get("name")
        if name in HEAD_TEMPLATES:
            args = template.get("args")
            if not isinstance(args, Mapping):
                return str(name), {}
            return str(name), {str(k): str(v) for k, v in args.items()}
    return "", {}


def _slot_order(key: str) -> tuple[int, int] | None:
    match = _SLOT_RE.fullmatch(key)
    if match is None:
        return None
    group, digits = match.groups()
    rank = {"2": 0, "pl": 0, "mpl": 1, "fpl": 2}[group]
    if group == "2":
        if digits:
            return None
        return rank, 1
    # pl2, pl3 … continue argument 2; a bare "pl" (it-adj) comes first.
    return rank, int(digits) if digits else 1


def _cell(key: str) -> str:
    if key.startswith("mpl"):
        return CELL_MPL
    if key.startswith("fpl"):
        return CELL_FPL
    return CELL_PL


def plural_slots(template: str, args: Mapping[str, str]) -> list[Slot]:
    """The plural-bearing arguments of a head, in the documented order.

    For ``it-noun``, a bare ``pl`` key is not a plural slot (argument
    2 is the first plural); for ``it-adj`` there is no argument 2. An
    empty argument (``|2=|``) counts as absent, as in the module.
    """
    keyed: list[tuple[tuple[int, int], str]] = []
    for key in args:
        order = _slot_order(key)
        if order is None or not args[key].strip():
            continue
        if template == NOUN_TEMPLATE and key == "pl":
            continue
        if template == ADJ_TEMPLATE and key == "2":
            continue
        keyed.append((order, key))
    keyed.sort()
    return [Slot(key, _cell(key), parse_list(args[key])) for _, key in keyed]


def has_plural_argument(template: str, args: Mapping[str, str]) -> bool:
    """True if the head states any plural explicitly (any slot)."""
    return bool(plural_slots(template, args))
