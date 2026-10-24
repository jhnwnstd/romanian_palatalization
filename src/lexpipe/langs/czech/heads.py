"""Gender, animacy and number of Czech headwords (design rule R13).

English Wiktionary gives a Czech noun's gender in the first argument of
its head template (``{{cs-noun|m-an}}``, ``{{cs-proper noun|f}}``),
which wiktextract keeps in ``head_templates[*].args["1"]``. The code is
a comma-separated list of readings; each reading is a gender with
optional hyphenated modifiers and an optional inline qualifier::

    m-an          masculine animate
    m-in          masculine inanimate
    f, n          feminine, neuter
    m-in-p        masculine inanimate plurale tantum (also f-p, n-p, m-an-p)
    mfbysense-an  masculine or feminine according to the referent
    m-an,m-in     both readings (animacy_ambiguous)
    f,m-in<l:rare>  a reading with a qualifier

Commas inside ``<…>`` belong to the qualifier
(``m-in,m-an<l:archaic,or,literary>``). :func:`parse_head_code` returns
every reading in source order; :class:`HeadInfo` derives the summary
fields the lexicon stores: the distinct genders and masculine
animacies (``|``-joined in source order by the writer), whether the
noun is a plurale tantum (every reading plural-only), whether gender
follows the sense (``mfbysense``) and whether the animacy or gender is
ambiguous (more than one distinct value).

The Czech Wiktionary edition has no head templates; it states gender
and animacy as entry-level tags (``["animate", "masculine"]``), read
with :func:`readings_from_tags`. Some of its entries lack animacy
(roh): the reading then has animacy ``""``.

The four MorfFlex gender letters a reading corresponds to (TR-64 tag
position 3) are given by :func:`morfflex_genders` for the R14 join on
(lemma, gender, animacy).
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any, Final

HEAD_TEMPLATES: Final[tuple[str, ...]] = (
    "cs-noun",
    "cs-proper noun",
    "cs-prop",
)
PROPER_HEAD_TEMPLATES: Final[frozenset[str]] = frozenset(
    {"cs-proper noun", "cs-prop"}
)

_GENDERS: Final[Mapping[str, tuple[str, ...]]] = {
    "m": ("m",),
    "f": ("f",),
    "n": ("n",),
    "mf": ("m", "f"),
    "mfbysense": ("m", "f"),
    "mfequiv": ("m", "f"),
}
_ANIMACY: Final[frozenset[str]] = frozenset({"an", "in", "pr"})
_PLURAL: Final[str] = "p"
UNKNOWN_GENDER: Final[str] = "?"

_TAG_GENDERS: Final[Mapping[str, str]] = {
    "masculine": "m",
    "feminine": "f",
    "neuter": "n",
}
_TAG_ANIMACY: Final[Mapping[str, str]] = {"animate": "an", "inanimate": "in"}

# MorfFlex (TR-64) gender letters: M masculine animate, I masculine
# inanimate, F, N; cover letters stand for several of them.
_MORFFLEX_COVER: Final[Mapping[str, frozenset[str]]] = {
    "M": frozenset({"M"}),
    "I": frozenset({"I"}),
    "F": frozenset({"F"}),
    "N": frozenset({"N"}),
    "Y": frozenset({"M", "I"}),
    "Z": frozenset({"M", "I", "N"}),
    "H": frozenset({"F", "N"}),
    "Q": frozenset({"F", "N"}),
    "T": frozenset({"I", "F"}),
    "X": frozenset({"M", "I", "F", "N"}),
}


@dataclass(frozen=True, slots=True)
class GenderReading:
    """One reading of a head code: ``m-in-p<l:rare>`` → (m, in, True, …)."""

    gender: str
    animacy: str = ""
    plural_only: bool = False
    qualifier: str = ""

    @property
    def label(self) -> str:
        """``m.an``, ``m.in``, ``m`` (animacy unknown), ``f`` or ``n``."""
        return f"{self.gender}.{self.animacy}" if self.animacy else self.gender


@dataclass(frozen=True, slots=True)
class HeadInfo:
    """Readings of one head code, with the derived summary fields."""

    code: str
    readings: tuple[GenderReading, ...]
    by_sense: bool = False

    @property
    def genders(self) -> tuple[str, ...]:
        return tuple(dict.fromkeys(r.gender for r in self.readings))

    @property
    def animacies(self) -> tuple[str, ...]:
        return tuple(
            dict.fromkeys(
                r.animacy
                for r in self.readings
                if r.gender == "m" and r.animacy
            )
        )

    @property
    def gender_ambiguous(self) -> bool:
        return len(self.genders) > 1

    @property
    def animacy_ambiguous(self) -> bool:
        return len(self.animacies) > 1

    @property
    def pluralia_tantum(self) -> bool:
        return bool(self.readings) and all(
            r.plural_only for r in self.readings
        )

    @property
    def labels(self) -> tuple[str, ...]:
        """Distinct reading labels (``m.an``, ``f`` …) in source order."""
        return tuple(dict.fromkeys(r.label for r in self.readings))


def _split_outside_brackets(code: str) -> list[str]:
    """Split on commas that are not inside ``<…>``."""
    parts: list[str] = []
    depth = 0
    current: list[str] = []
    for ch in code:
        if ch == "<":
            depth += 1
        elif ch == ">" and depth:
            depth -= 1
        if ch == "," and depth == 0:
            parts.append("".join(current))
            current = []
        else:
            current.append(ch)
    parts.append("".join(current))
    return [part.strip() for part in parts if part.strip()]


def _parse_reading(item: str) -> tuple[list[GenderReading], bool]:
    """Readings of one comma item and whether it was by-sense."""
    qualifier = ""
    if "<" in item:
        cut = item.index("<")
        item, qualifier = item[:cut], item[cut:]
    tokens = [token for token in item.strip().split("-") if token]
    if not tokens:
        return [GenderReading(UNKNOWN_GENDER, qualifier=qualifier)], False
    genders = _GENDERS.get(tokens[0], (UNKNOWN_GENDER,))
    animacy = next((t for t in tokens[1:] if t in _ANIMACY), "")
    plural_only = _PLURAL in tokens[1:]
    readings = [
        GenderReading(
            gender,
            animacy if gender == "m" else "",
            plural_only,
            qualifier,
        )
        for gender in genders
    ]
    return readings, tokens[0] == "mfbysense"


def parse_head_code(code: str) -> HeadInfo:
    """Parse a cs-noun gender code (see the module docstring)."""
    readings: list[GenderReading] = []
    by_sense = False
    for item in _split_outside_brackets(code):
        parsed, sense = _parse_reading(item)
        readings.extend(parsed)
        by_sense = by_sense or sense
    return HeadInfo(code=code, readings=tuple(readings), by_sense=by_sense)


def head_template(
    entry: Mapping[str, Any],
) -> tuple[str, dict[str, str]] | None:
    """Name and arguments of the entry's first Czech noun head template."""
    for template in entry.get("head_templates") or ():
        if not isinstance(template, Mapping):
            continue
        name = template.get("name")
        if name in HEAD_TEMPLATES:
            args = template.get("args")
            if not isinstance(args, Mapping):
                return str(name), {}
            return str(name), {str(k): str(v) for k, v in args.items()}
    return None


def readings_from_tags(tags: Iterable[str]) -> tuple[GenderReading, ...]:
    """Readings from gender/animacy tags (cs edition, sense tags).

    ``["animate", "masculine"]`` → (m, an); ``["feminine"]`` → (f,);
    ``["masculine"]`` alone → (m, "") (animacy unknown); both animacy
    tags give both readings.
    """
    tag_set = set(tags)
    animacies = [code for tag, code in _TAG_ANIMACY.items() if tag in tag_set]
    plural_only = "plural-only" in tag_set or "plural only" in tag_set
    readings: list[GenderReading] = []
    for tag, gender in _TAG_GENDERS.items():
        if tag not in tag_set:
            continue
        if gender == "m" and animacies:
            readings.extend(
                GenderReading(gender, animacy, plural_only)
                for animacy in animacies
            )
        else:
            readings.append(GenderReading(gender, "", plural_only))
    return tuple(readings)


def head_info(entry: Mapping[str, Any]) -> HeadInfo:
    """Gender readings of an en.wiktionary entry.

    The head template's code wins; without one, the tags of the senses
    are used (first sense with any gender tag), and failing that an
    empty :class:`HeadInfo`.
    """
    found = head_template(entry)
    code = found[1].get("1", "") if found is not None else ""
    if code:
        return parse_head_code(code)
    for sense in entry.get("senses") or ():
        if not isinstance(sense, Mapping):
            continue
        readings = readings_from_tags(
            tag for tag in sense.get("tags") or () if isinstance(tag, str)
        )
        if readings:
            return HeadInfo(code="", readings=readings)
    return HeadInfo(code="", readings=())


def morfflex_genders(reading: GenderReading) -> frozenset[str]:
    """MorfFlex gender letters (tag position 3) matching a reading.

    Masculine with unknown animacy matches both M and I.
    """
    if reading.gender == "m":
        if reading.animacy == "an":
            return frozenset({"M"})
        if reading.animacy == "in":
            return frozenset({"I"})
        return frozenset({"M", "I"})
    if reading.gender == "f":
        return frozenset({"F"})
    if reading.gender == "n":
        return frozenset({"N"})
    return frozenset()


def morfflex_gender_matches(letter: str, wanted: frozenset[str]) -> bool:
    """Does a MorfFlex gender letter (possibly a cover letter) overlap
    the wanted set of M/I/F/N?"""
    return bool(_MORFFLEX_COVER.get(letter, frozenset()) & wanted)
