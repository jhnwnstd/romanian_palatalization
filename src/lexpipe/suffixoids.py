"""The -log/-fag/-urg table's two classifiers, for every language (R27).

One learned (Greek) suffixoid, two plural behaviours in all three
languages: person nouns palatalize or take their own allomorph
(Romanian psihologi, Italian teologi, Czech biologové/biolozích),
thing nouns do not (dialoguri, cataloghi, katalogy). R27 asks for one
-log/-fag/-urg table per language "with the same person/thing
classifier", so the classifiers live here and the plugins supply only
language facts: which lemmas carry a suffixoid (their endings differ:
-log, -logo, -log), the evidence they have that it is the learned
suffixoid, extra gloss patterns for glosses that are not English, and
a hand list of decisions in ``curated/<language>/``.

Learned suffixoid or lookalike
------------------------------
Folk words and simplexes can end like the suffixoid (Romanian bârlog
'den', amurg 'dusk'; Czech profil 'profile' looks like -fil).
:func:`learned_suffixoid` decides, in this order, and returns the rule
that decided as the basis:

  1. ``manual``: a curated decision;
  2. the plugin's evidence (an etymology template or etymon naming the
     suffixoid, an Italian masculine ending — whatever the language
     plugin documents; the basis is the evidence label);
  3. ``gloss``: the English gloss carries the suffixoid's reflex
     (:data:`GLOSS_REFLEX`: log in psychologist, dialogue, catalogue;
     phag; urg in surgeon, dramaturge; graph; phil);
  4. ``gloss (profession)``: an -olog lemma glossed as a learned
     profession (africanist, political scientist; the plugin may add a
     pattern for glosses in its own language: Czech "lékař", "odborník");
  5. otherwise a lookalike (basis ``gloss``).

Lookalikes stay in the per-noun tables, flagged, and are left out of
the aggregates.

Person or thing
---------------
:func:`semantic_class` reads the first gloss (before ``;``) with these
rules, in order, and returns ``(class, basis)``:

  1. ``manual``: a curated class;
  2. ``animate`` (:data:`ANIMATE`: animal, herbivore, insect,
     scavenger, …phagous);
  3. ``thing`` (:data:`THING`: …logue (not ideologue), dialog, catalog,
     register, sarcophag, phage, cave as a word (not scavenger), …);
  4. ``person`` (:data:`PERSON`: agent nouns in -ist, -ologer, -ogist,
     -er, -or, -ian, -ant, -ent, -eon, and person words: person, man,
     expert, surgeon, playwright, cannibal, oracle, …);
  5. the plugin's extra patterns (``gloss (<language>)`` basis: Czech
     glosses of cs.wiktionary lexemes);
  6. ``unknown``.

Gender and animacy are not used: they are close to the plural
allomorph itself (Romanian masc. -i vs neut. -uri/-e; Czech -ové vs -y),
which would make the classifier circular. A plugin may add a class
outside :data:`CLASSES` before calling it (Italian adjectives in
-logo are ``adjective``).

Curated decisions
-----------------
:func:`load_decisions` reads ``curated/<language>/suffixoid_classes.csv``
with :mod:`lexpipe.curated` (guard columns, key column chosen by the
plugin: ``lexeme_id`` for Italian and Czech, ``lemma`` for Romanian,
which has no lexeme_ids). Columns: the key, ``learned_suffixoid``
(``True``, ``False`` or empty), ``semantic_class`` (one of
:data:`CLASSES` or empty), ``decided_by``, ``note`` and ``guard_*``.

History: the classifiers were written for the harmonized Romanian table
(``lexpipe.langs.romanian_compare.logtable``) with the manual lists in
code; the Czech plugin imported them, and the Italian plugin had its own
gloss rule (all glosses; classes person, thing, mixed, adjective).
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from lexpipe import csvio, curated

PERSON_CLASS: Final[str] = "person"
ANIMATE_CLASS: Final[str] = "animate"
THING_CLASS: Final[str] = "thing"
UNKNOWN_CLASS: Final[str] = "unknown"
CLASSES: Final[tuple[str, ...]] = (
    PERSON_CLASS,
    ANIMATE_CLASS,
    THING_CLASS,
    UNKNOWN_CLASS,
)

# Suffixoid (written without its inflectional vowel) → the reflex an
# English gloss of a learned formation carries.
GLOSS_REFLEX: Final[Mapping[str, str]] = {
    "log": "log",
    "fag": "phag",
    "urg": "urg",
    "graf": "graph",
    "fil": "phil",
}

PERSON: Final[re.Pattern[str]] = re.compile(
    r"\b(\w{2,}(ists?|ologer|ogist|er|or|ian|ant|ent|eon)|person|man|"
    r"woman|expert|specialist|scholar|surgeon|doctor|dentist|playwright|"
    r"dramatist|beggar|simpleton|trifler|stupid|official|thaumaturge|"
    r"demiurge|ideologue|metalworker|cadger|counterpart|cannibal|oracle|"
    r"wizard|cheat|liar|know-it-all|someone|somebody)\b"
)
ANIMATE: Final[re.Pattern[str]] = re.compile(
    r"\b(animal|herbivore|phytophage|organism|insect|scavengers?)\b|phagous"
)
THING: Final[re.Pattern[str]] = re.compile(
    r"((?<!ideo)logue|dialog\b|catalog|obituary|register|net\b|swath|"
    r"\bcaves?\b|suburb|town|dusk|sarcophag|oesophag|esophag|phage\b|"
    r"thistle|vlog)"
)
LEARNED_PROFESSION: Final[re.Pattern[str]] = re.compile(
    r"\b\w+ists?\b|scientist|scholar|expert|specialist"
)
_OLOG: Final[re.Pattern[str]] = re.compile(r"olog[oa]?\Z")

DECISIONS_FILE: Final[str] = "suffixoid_classes.csv"
LEARNED_FIELD: Final[str] = "learned_suffixoid"
CLASS_FIELD: Final[str] = "semantic_class"
DECISION_FIELDS: Final[tuple[str, ...]] = (
    LEARNED_FIELD,
    CLASS_FIELD,
    "decided_by",
    "note",
)


@dataclass(frozen=True, slots=True)
class GlossPatterns:
    """Extra gloss rules of a language whose glosses are not English
    (``label`` names them in the basis: ``gloss (cs)``)."""

    label: str
    person: re.Pattern[str] | None = None
    thing: re.Pattern[str] | None = None
    profession: re.Pattern[str] | None = None


@dataclass(frozen=True, slots=True)
class Decision:
    """A curated decision: ``learned`` is None and ``semantic_class``
    empty when the curator decided only the other one."""

    learned: bool | None = None
    semantic_class: str = ""
    note: str = ""


def first_gloss(gloss: str) -> str:
    """The first gloss (before ``;``), trimmed."""
    return gloss.split(";")[0].strip()


def learned_suffixoid(
    lemma: str,
    suffixoid: str,
    gloss: str,
    *,
    decision: Decision | None = None,
    evidence: str = "",
    extra: GlossPatterns | None = None,
) -> tuple[bool, str]:
    """``(is the learned suffixoid, basis)`` (module docstring).

    *suffixoid* is a key of :data:`GLOSS_REFLEX` (``log``, ``fag``,
    ``urg``, ``graf``, ``fil``); *gloss* the lexeme's gloss (all of it:
    the reflex may be in a later gloss); *evidence* the plugin's label
    of its own evidence ("" when it has none).
    """
    if decision is not None and decision.learned is not None:
        return decision.learned, "manual"
    if evidence:
        return True, evidence
    try:
        reflex = GLOSS_REFLEX[suffixoid]
    except KeyError:
        raise ValueError(
            f"unknown suffixoid {suffixoid!r}; known: {sorted(GLOSS_REFLEX)}"
        ) from None
    text = gloss.lower()
    if reflex in text:
        return True, "gloss"
    if _OLOG.search(lemma):
        if LEARNED_PROFESSION.search(text):
            return True, "gloss (profession)"
        if (
            extra is not None
            and extra.profession is not None
            and extra.profession.search(text)
        ):
            return True, f"gloss (profession, {extra.label})"
    return False, "gloss"


def semantic_class(
    gloss: str,
    *,
    decision: Decision | None = None,
    extra: GlossPatterns | None = None,
) -> tuple[str, str]:
    """``(class, basis)`` of a -log/-fag/-urg noun (module docstring),
    from its first gloss."""
    if decision is not None and decision.semantic_class:
        return decision.semantic_class, "manual"
    text = first_gloss(gloss).lower()
    if not text:
        return UNKNOWN_CLASS, "no gloss"
    if ANIMATE.search(text):
        return ANIMATE_CLASS, "gloss"
    if THING.search(text):
        return THING_CLASS, "gloss"
    if PERSON.search(text):
        return PERSON_CLASS, "gloss"
    if extra is not None:
        if extra.person is not None and extra.person.search(text):
            return PERSON_CLASS, f"gloss ({extra.label})"
        if extra.thing is not None and extra.thing.search(text):
            return THING_CLASS, f"gloss ({extra.label})"
    return UNKNOWN_CLASS, "gloss"


def decision_of(row: curated.CuratedRow) -> Decision:
    """The :class:`Decision` of one curated row; :class:`ValueError` on
    a value outside the allowed ones."""
    learned_text = row.values.get(LEARNED_FIELD, "")
    learned = csvio.parse_bool(learned_text) if learned_text else None
    klass = row.values.get(CLASS_FIELD, "").strip()
    if klass and klass not in CLASSES:
        raise ValueError(
            f"{'/'.join(row.key) or row.lexeme_id}: semantic_class "
            f"{klass!r} is not one of {CLASSES}"
        )
    return Decision(learned, klass, row.values.get("note", ""))


def load_decisions(
    path: Path, key_field: str = "lexeme_id"
) -> tuple[dict[str, curated.CuratedRow], dict[str, Decision]]:
    """The curated rows of *path* (keyed by *key_field*) and their
    decisions. A missing file means no decisions. The caller checks the
    guards (:func:`lexpipe.curated.check_guards`) against its own
    data."""
    rows = curated.load_curated(path, key_field)
    for key, row in rows.items():
        missing = [
            field for field in DECISION_FIELDS if field not in row.values
        ]
        if missing:
            raise ValueError(f"{path}: {key}: missing columns {missing}")
    return rows, {key: decision_of(row) for key, row in rows.items()}
