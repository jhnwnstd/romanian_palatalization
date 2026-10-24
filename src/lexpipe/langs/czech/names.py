"""Proper names versus common nouns (design rule R25).

Ethnonyms and demonyms are common nouns and stay in the domain: Czech
capitalizes them (Čech, Rus, Slovák) and they carry much of the
NOM.PL -i/-ové variation (Češi/Čechové). Given names, surnames and
toponyms are out of the domain; the lexicon keeps them as rows with
``domain = out`` and records which class of name they are.

English Wiktionary
------------------
An entry is a proper name when wiktextract's ``pos`` is ``name``, when
its head template is ``cs-proper noun``/``cs-prop`` (some given names
sit under a Noun header: Inka 'a female given name' next to Inka
'Inca'), or when every gloss is a name gloss ("a male surname").
Category names are **not** used to decide this, because categories are
page-level: the common noun Čech carries "Czech proper nouns" from the
surname entry on the same page. They are used only to classify an
entry already known to be a name, preferring its own glosses, then its
sense-level categories, then the entry-level ones:

    surname      gloss "… surname …" or "Czech (male|female) surnames"
    given_name   gloss "… given name …" or "Czech … given names" /
                 "Czech diminutives of … given names"
    toponym      gloss "(a|the) country|city|town|river …" or a topical
                 place category ("cs:Places in …", "cs:Rivers …")
    other_name   anything else (deities, planets, individuals)

A **capitalized** entry under a Noun header is also a name when its
first gloss names a place or a sky object rather than a person ("Cancer
(constellation)", "Franconia (geographical region in Germany)",
"Prague (the capital city …)", "Brno Reservoir") or its sense
categories are topical place or constellation categories (Býk
'Taurus': cs:Constellations in the zodiac), unless the gloss is a
demonym gloss ("Akkadian (inhabitant of the Akkad region)", "native
or inhabitant of …"). These capitalized homographs (Rak, Býk, Liška
the constellations; Pha, Franky, Jesenicko) are not the common nouns
rak, býk, liška, whose derivatives (račí, býčí, liščí) they would
otherwise share.

Czech Wiktionary
----------------
The cs edition files proper names under ``podstatné jméno`` with the raw
tags ``vlastní jméno`` ('proper name'), but tags ethnonyms the same way
(Rus 'obyvatel Ruska'). Its entries are classified by their Czech
glosses: příjmení (surname), křestní/mužské/ženské jméno (given name),
place words (město, obec, řeka, stát, ostrov …), and ``vlastní jméno``
is ``other_name`` unless some gloss is a demonym gloss ("obyvatel…",
"příslušník…", "… pocházející z…", "… národnosti"), which makes the
entry a common noun. Many cs entries of names carry no raw tag at all
(Anuška 'ženské křestní jméno', Studénka 'město v okrese Nový Jičín',
Jindřišská 'ulice …', Azor 'jméno pro psa', Jupiter 'pátá planeta …'):
a **capitalized** cs entry whose first gloss names a person's name
type anywhere in it, a place (the place words above and ulice,
náměstí, atol, kanton, prefektura …) or a dog's name is a name too,
unless the gloss is a demonym gloss.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any, Final

from lexpipe import kaikki
from lexpipe.langs.czech.heads import PROPER_HEAD_TEMPLATES

SURNAME: Final[str] = "surname"
GIVEN_NAME: Final[str] = "given_name"
TOPONYM: Final[str] = "toponym"
OTHER_NAME: Final[str] = "other_name"
NAME_CLASSES: Final[tuple[str, ...]] = (
    SURNAME,
    GIVEN_NAME,
    TOPONYM,
    OTHER_NAME,
)

_EN_SURNAME_GLOSS: Final[re.Pattern[str]] = re.compile(
    r"\bsurnames?\b|\bpatronymic\b", re.IGNORECASE
)
_EN_GIVEN_GLOSS: Final[re.Pattern[str]] = re.compile(
    r"\bgiven names?\b|\bdiminutive of the (?:male |female )?given name",
    re.IGNORECASE,
)
_EN_PLACE_GLOSS: Final[re.Pattern[str]] = re.compile(
    r"\b(?:a|an|the)\s+(?:\w+\s+)?(?:country|city|town|village|river|"
    r"island|islands|archipelago|region|state|province|capital|mountain|"
    r"mountains|lake|sea|ocean|continent|municipality|district|peninsula|"
    r"county|territory|microstate|kingdom|republic|desert|valley|"
    r"neighbourhood|neighborhood|borough|suburb|spa town)\b",
    re.IGNORECASE,
)
_EN_NAME_GLOSS: Final[re.Pattern[str]] = re.compile(
    r"^(?:an?|the)\s+(?:\w+\s+)?(?:male|female|unisex)\s+"
    r"(?:surname|given name)|^an?\s+(?:surname|given name)|"
    r"^a diminutive of the (?:male |female )?given name",
    re.IGNORECASE,
)
_EN_SURNAME_CATEGORY: Final[re.Pattern[str]] = re.compile(
    r"^Czech (?:male |female )?surnames\b"
)
_EN_GIVEN_CATEGORY: Final[re.Pattern[str]] = re.compile(
    r"^Czech (?:male |female |unisex )?given names\b|"
    r"^Czech diminutives of (?:male |female |unisex )?given names\b"
)
_TOPICAL_PLACE: Final[re.Pattern[str]] = re.compile(
    r"^cs:(?:Places|Cities|Towns|Villages|Countries|Rivers|Islands|"
    r"Mountains|Seas|Lakes|Regions|Provinces|States|Continents|"
    r"National capitals|Geographic and cultural areas|Municipalities|"
    r"Districts|Peninsulas|Former polities|Federal|Republics|Counties|"
    r"Departments|Autonomous|Cantonal capitals|Krai capitals|"
    r"Regional capitals|Neighbourhoods|Dependent territories|"
    r"Constituent countries|Ancient settlements|Oblasts|Regions of|"
    r"Capital cities|Deserts|Valleys|Straits|Bays|Oceans|Spa towns)"
)

# Name glosses start with the name type, at most three words in:
# "české mužské příjmení", "mužské křestní jméno", "domácká podoba …".
_CS_SURNAME: Final[re.Pattern[str]] = re.compile(
    r"^(?:[\w-]+\s+){0,3}příjmení\b", re.IGNORECASE
)
_CS_GIVEN: Final[re.Pattern[str]] = re.compile(
    r"^(?:[\w-]+\s+){0,3}(?:mužské|ženské|křestní)\s+(?:křestní\s+)?"
    r"jméno\b|^domácká podoba\b",
    re.IGNORECASE,
)
_CS_PLACE: Final[re.Pattern[str]] = re.compile(
    r"\b(?:město|městys|obec|vesnice|ves|řeka|potok|stát|ostrov|souostroví|"
    r"pohoří|hora|vrch|hlavní město|území|kraj|okres|země|světadíl|moře|"
    r"oceán|jezero|poloostrov|region|provincie|oblast|čtvrť|"
    r"městská část|údolí|nížina|vrchovina|planeta|souhvězdí)\b",
    re.IGNORECASE,
)
# Demonym glosses: "obyvatel Ruska", "muž pocházející ze Slovenska",
# "člověk hlásící se k české národnosti".
_CS_DEMONYM: Final[re.Pattern[str]] = re.compile(
    r"^(?:obyvatel|obyvatelka|příslušník|příslušnice|rodák|rodačka|"
    r"občan|občanka)\b|\bpocházející\b|\bnárodnost",
    re.IGNORECASE,
)
_CS_PROPER_RAW_TAG: Final[str] = "vlastní"
# Capitalized cs entries without a raw tag (module docstring): a name
# type anywhere in the first gloss ("Anna; ženské křestní jméno",
# "mužské osobní jméno"), extra place words, an animal's name.
_CS_GIVEN_ANYWHERE: Final[re.Pattern[str]] = re.compile(
    r"\b(?:mužské|ženské|křestní|rodné|osobní)\s+(?:(?:křestní|rodné|"
    r"osobní)\s+)?jméno\b|^domácká podoba\b|\bzdrobnělina\s+jména\b",
    re.IGNORECASE,
)
_CS_SURNAME_ANYWHERE: Final[re.Pattern[str]] = re.compile(
    r"\bpříjmení\b", re.IGNORECASE
)
_CS_PLACE_EXTRA: Final[re.Pattern[str]] = re.compile(
    r"\b(?:ulice|náměstí|nábřeží|atol|kanton|prefektura|departement|"
    r"sídlo|souhvězdí|planeta|hvězda)\b",
    re.IGNORECASE,
)
_CS_ANIMAL_NAME: Final[re.Pattern[str]] = re.compile(
    r"\bjméno\s+pro\s+(?:psa|fenku|kočku|kocoura|koně|zvíře)\b",
    re.IGNORECASE,
)
# Capitalized en entries under a Noun header (module docstring).
_EN_PLACE_OR_SKY: Final[re.Pattern[str]] = re.compile(
    r"\b(?:region|reservoir|capital city|constellation|astrological sign|"
    r"zodiac(?:al)? sign|planet|country|city|town|village|river|"
    r"mountains?|lake|island|street|square|district|province)\b",
    re.IGNORECASE,
)
_EN_DEMONYM: Final[re.Pattern[str]] = re.compile(
    r"\b(?:inhabitant|inhabitants|native|natives|person|people|member|"
    r"members|citizen|descendant|someone|man|woman|female|male|speaker|"
    r"follower|adherent|tribe|resident)\b",
    re.IGNORECASE,
)
_TOPICAL_SKY: Final[re.Pattern[str]] = re.compile(
    r"^cs:(?:Constellations|Astrology|Planets|Stars)\b"
)


@dataclass(frozen=True, slots=True)
class NameVerdict:
    """Whether an entry is a proper name, which class, and why."""

    is_name: bool
    name_class: str = ""
    evidence: str = ""


COMMON: Final[NameVerdict] = NameVerdict(False)


def _glosses(entry: Mapping[str, Any]) -> list[str]:
    out: list[str] = []
    for sense in entry.get("senses") or ():
        if not isinstance(sense, Mapping):
            continue
        for gloss in sense.get("glosses") or ():
            if isinstance(gloss, str) and gloss.strip():
                out.append(gloss.strip())
                break
    return out


def _sense_categories(entry: Mapping[str, Any]) -> list[str]:
    out: list[str] = []
    for sense in entry.get("senses") or ():
        if isinstance(sense, Mapping):
            out.extend(
                c for c in sense.get("categories") or () if isinstance(c, str)
            )
    return out


def _entry_categories(entry: Mapping[str, Any]) -> list[str]:
    return [c for c in entry.get("categories") or () if isinstance(c, str)]


def _class_from_glosses(glosses: Iterable[str]) -> str:
    for gloss in glosses:
        if _EN_SURNAME_GLOSS.search(gloss):
            return SURNAME
        if _EN_GIVEN_GLOSS.search(gloss):
            return GIVEN_NAME
        if _EN_PLACE_GLOSS.search(gloss):
            return TOPONYM
    return ""


def _class_from_categories(categories: Iterable[str]) -> str:
    found = list(categories)
    if any(_EN_SURNAME_CATEGORY.match(c) for c in found):
        return SURNAME
    if any(_EN_GIVEN_CATEGORY.match(c) for c in found):
        return GIVEN_NAME
    if any(_TOPICAL_PLACE.match(c) for c in found):
        return TOPONYM
    return ""


def en_name_class(entry: Mapping[str, Any]) -> str:
    """Class of an en entry already known to be a name."""
    return (
        _class_from_glosses(_glosses(entry))
        or _class_from_categories(_sense_categories(entry))
        or _class_from_categories(_entry_categories(entry))
        or OTHER_NAME
    )


def _capitalized(entry: Mapping[str, Any]) -> bool:
    word = str(entry.get("word") or "")
    return word[:1].isupper()


def capitalized_place_en(entry: Mapping[str, Any]) -> str:
    """Name class of a capitalized en noun whose first gloss or sense
    categories name a place or a sky object (module docstring); ""
    otherwise."""
    if not _capitalized(entry):
        return ""
    glosses = _glosses(entry)
    if not glosses or _EN_DEMONYM.search(glosses[0]):
        return ""
    categories = _sense_categories(entry)
    if any(_TOPICAL_SKY.match(c) for c in categories):
        return OTHER_NAME
    if any(_TOPICAL_PLACE.match(c) for c in categories):
        return TOPONYM
    if _EN_PLACE_OR_SKY.search(glosses[0]):
        sky = re.search(
            r"constellation|astrological|zodiac|planet", glosses[0], re.I
        )
        return OTHER_NAME if sky else TOPONYM
    return ""


def classify_en(entry: Mapping[str, Any]) -> NameVerdict:
    """R25 verdict for an English Wiktionary entry."""
    evidence = ""
    if entry.get("pos") == "name":
        evidence = "pos"
    elif any(
        isinstance(t, Mapping) and t.get("name") in PROPER_HEAD_TEMPLATES
        for t in entry.get("head_templates") or ()
    ):
        evidence = "head_template"
    else:
        glosses = _glosses(entry)
        if glosses and all(_EN_NAME_GLOSS.search(g) for g in glosses):
            evidence = "gloss"
    if not evidence:
        place = capitalized_place_en(entry)
        if place:
            return NameVerdict(True, place, "capitalized_gloss")
        return COMMON
    return NameVerdict(True, en_name_class(entry), evidence)


def classify_cs(entry: Mapping[str, Any]) -> NameVerdict:
    """R25 verdict for a Czech Wiktionary entry (Czech glosses)."""
    glosses = _glosses(entry)
    for gloss in glosses:
        if _CS_SURNAME.search(gloss):
            return NameVerdict(True, SURNAME, "gloss")
        if _CS_GIVEN.search(gloss):
            return NameVerdict(True, GIVEN_NAME, "gloss")
    raw_tags = {t for t in entry.get("raw_tags") or () if isinstance(t, str)}
    if entry.get("pos") == "name":
        return NameVerdict(True, _cs_place_or_other(glosses), "pos")
    if _CS_PROPER_RAW_TAG in raw_tags:
        if any(_CS_DEMONYM.search(g) for g in glosses):
            return COMMON
        return NameVerdict(True, _cs_place_or_other(glosses), "raw_tags")
    return capitalized_name_cs(entry, glosses)


def capitalized_name_cs(
    entry: Mapping[str, Any], glosses: list[str]
) -> NameVerdict:
    """A capitalized cs entry whose first gloss names a name type, a
    place or an animal's name (module docstring); else common."""
    if not _capitalized(entry) or not glosses:
        return COMMON
    gloss = glosses[0]
    if _CS_DEMONYM.search(gloss):
        return COMMON
    if _CS_SURNAME_ANYWHERE.search(gloss):
        return NameVerdict(True, SURNAME, "capitalized_gloss")
    if _CS_GIVEN_ANYWHERE.search(gloss):
        return NameVerdict(True, GIVEN_NAME, "capitalized_gloss")
    if _CS_PLACE.search(gloss) or _CS_PLACE_EXTRA.search(gloss):
        sky = re.search(r"souhvězdí|planeta|hvězda", gloss, re.I)
        return NameVerdict(
            True, OTHER_NAME if sky else TOPONYM, "capitalized_gloss"
        )
    if _CS_ANIMAL_NAME.search(gloss):
        return NameVerdict(True, OTHER_NAME, "capitalized_gloss")
    return COMMON


def _cs_place_or_other(glosses: Iterable[str]) -> str:
    return TOPONYM if any(_CS_PLACE.search(g) for g in glosses) else OTHER_NAME


def first_gloss(entry: Mapping[str, Any]) -> str:
    """The first gloss (shared with :func:`lexpipe.ids.first_gloss`)."""
    glosses = _glosses(entry)
    return glosses[0] if glosses else ""


def czech_categories(entry: Mapping[str, Any]) -> list[str]:
    """Czech-scoped categories of an en entry (entry level, then senses)."""
    return kaikki.categories(entry, lang_name="Czech")
