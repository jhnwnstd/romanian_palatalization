"""Stage ``verbs``: velar-stem verbs and their front-vowel cells (R15).

Design §5 ("Verbs"): -are never palatalizes; root-stressed -ere
alternates; -ére plus cuocere and nuocere over-apply; -ire is mixed;
-isc- never; dire is suppletive. The verb cells are learner evidence
E1 (R15: "Italian: 2sg, 3sg and 1pl of every conjugation") and the
units of the verb TP rows.

Inputs: the ``extract`` projections ``italian_en_verbs.jsonl.gz`` (verb
lemmas: accented canonical infinitive, categories, etymology) and the
form index (conjugation-table rows), ``italian_lexicon_complete.csv``
(dorsal noun/adjective roots, for the E1 link) and the Leipzig words
table (attestation, R11).

Velar-stem verbs
----------------
A verb lemma qualifies when its infinitive key (:mod:`.keys`; a
reflexive -arsi/-ersi/-irsi read as -are/-ere/-ire) minus the ending
ends in c or g, the orthographic velar spelling ch/gh counting as c/g
(cercare, vincere, bianchire): "infinitives with stem-final c/g". The
dire family (dire, benedire, predire …: a 1sg in -dico) is added with
the stem of its 1sg. Stems in ci/gi (lanciare) and chi/ghi (macchiare)
do not qualify. ``stem_key`` is the stem with that h removed (bianc);
``target_letters`` is sc, cc, gg, c or g.

Conjugation class (``conj_class``)
----------------------------------
are             -are (-arsi)
ere_root        -ere with a root-stressed infinitive: from the IPA of
                the infinitive (/ˈvin.t͡ʃe.re/: stress three syllables
                from the end; :func:`lexpipe.normalize.stress_from_end`),
                else the category "Italian verbs with root-stressed
                infinitive", else an accent before the ending in the
                accented canonical infinitive (vìncere, conóscere);
                ``infinitive_stress_source`` says which
ere_stressed    -ére (tacére, piacére: /taˈt͡ʃe.re/)
ere_unknown     -ere without stress information
ire             -ire without -isc- (fuggire, cucire, uscire)
isc             -ire with -isc- (a 1sg stem + isco: bianchisco)
dire            the dire family (suppletive dic- stem)

Cells
-----
From the conjugation-table rows (form-index source ``table``; the
clitic pronoun of a reflexive form is dropped, ti amìchi → amìchi;
forms labelled archaic, literary, proscribed … — the ``features``
stage's :data:`~.features.FILTERED_LABELS` — are left out):

    1sg      indicative present (back vowel -o; -isco)
    2sg      indicative present (-i)                  E1
    3sg      indicative present (-e; -a for -are)     E1
    1pl      indicative present (-iamo: /j/)          E1
    sbjv3sg  present subjunctive (-i for -are, -a otherwise)  E1
    past1sg  past historic (tacqui: velar before /w/)
    inf      infinitive

Each form is parsed at its ending (:func:`parse_form`): ``vowel`` is
the first segment of the ending (i, e, or back), ``edge`` is
``velar`` (ch/gh before a front vowel, c/g before a back vowel or qu),
``palatal`` (c/g before a front vowel, or ci/gi before a back vowel:
taccio) or ``unparsed``.

Paradigm type (verb level, over the parsed cells)
-------------------------------------------------
alternating         a velar before a back vowel and a palatal before a
                    front vowel (vinco/vinci, dico/dici)
velar_throughout    velar in every front-vowel cell (cerchi, bianchisci)
overapplying        the 1sg is palatal before -o and some cell is velar
                    (taccio/tacqui, noccio/nocqui)
palatal_throughout  no velar anywhere (cucire, agire, cuocere): no
                    evidence of a velar root; out of the TP domain
mixed               palatal and velar front-vowel cells (variants)
no_cells            no parsable cell

``in_domain`` excludes palatal_throughout and no_cells. ``noun_roots``
lists the dorsal noun/adjective roots whose stem equals the verb's
stem key (traffico ↔ trafficare, amico ↔ amicare/amicarsi, bianco ↔
bianchire): the synchronic E1 link (``shares_root_with_1st_conj_verb``
generalized to every conjugation). When several roots share the stem
(pesca 'fishing' and pesca 'peach' for pescare), only the roots with
etymological evidence of the link are kept — a lexeme of the root is
deverbal from the verb (pesca 'fishing'), or the verb's etymology
names a lemma of the root (``{{af|it|secco|-are}}``); without such
evidence for any of them, all stay linked (:func:`linked_roots`).

Outputs (gitignored: Leipzig counts)
-----------------------------------
``data/italian/italian_verbs.csv`` (:data:`VERB_FIELDS`): one row per
(verb, cell, form), sorted by verb_key, verb_id, cell order, form;
``data/italian/italian_verbs_summary.tsv``: counts.
"""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

from lexpipe import attest, csvio, normalize
from lexpipe.attest import AttestationTable
from lexpipe.context import RunContext, run_stage
from lexpipe.langs.italian import (
    etymology,
    extract,
    features,
    keys,
    leipzig,
    morphology,
    store,
)
from lexpipe.langs.italian.layout import Layout

STAGE: Final[str] = "verbs"
# 2: -ere stress from the IPA (infinitive_stress_source).
# 3: the canonical fallback reads "-ére [note]" as ending-stressed.
# 4: the shared etymology classifier (lexpipe.etymology).
# 5: names filter on the τ flags (leipzig.attested_flags); homograph
#    roots of a stem linked by etymological evidence (linked_roots).
VERSION: Final[int] = 5
SOURCE: Final[str] = "wiktionary_en|leipzig"

ARE: Final[str] = "are"
ERE_ROOT: Final[str] = "ere_root"
ERE_STRESSED: Final[str] = "ere_stressed"
ERE_UNKNOWN: Final[str] = "ere_unknown"
IRE: Final[str] = "ire"
ISC: Final[str] = "isc"
DIRE: Final[str] = "dire"
CONJ_CLASSES: Final[tuple[str, ...]] = (
    ARE,
    ERE_ROOT,
    ERE_STRESSED,
    ERE_UNKNOWN,
    IRE,
    ISC,
    DIRE,
)

ALTERNATING: Final[str] = "alternating"
VELAR_THROUGHOUT: Final[str] = "velar_throughout"
OVERAPPLYING: Final[str] = "overapplying"
PALATAL_THROUGHOUT: Final[str] = "palatal_throughout"
MIXED: Final[str] = "mixed"
NO_CELLS: Final[str] = "no_cells"
PARADIGM_TYPES: Final[tuple[str, ...]] = (
    ALTERNATING,
    VELAR_THROUGHOUT,
    OVERAPPLYING,
    PALATAL_THROUGHOUT,
    MIXED,
    NO_CELLS,
)
OUT_OF_DOMAIN: Final[frozenset[str]] = frozenset(
    {PALATAL_THROUGHOUT, NO_CELLS}
)

VELAR: Final[str] = "velar"
PALATAL: Final[str] = "palatal"
UNPARSED: Final[str] = "unparsed"

CELL_ORDER: Final[tuple[str, ...]] = (
    "inf",
    "1sg",
    "2sg",
    "3sg",
    "1pl",
    "sbjv3sg",
    "past1sg",
)
E1_CELLS: Final[frozenset[str]] = frozenset({"2sg", "3sg", "1pl", "sbjv3sg"})
# (required tags, excluded tags) of each cell in the conjugation table.
_CELL_TAGS: Final[Mapping[str, tuple[frozenset[str], frozenset[str]]]] = {
    "inf": (frozenset({"infinitive"}), frozenset({"imperative"})),
    "1sg": (
        frozenset({"indicative", "present", "first-person", "singular"}),
        frozenset({"imperative", "subjunctive"}),
    ),
    "2sg": (
        frozenset({"indicative", "present", "second-person", "singular"}),
        frozenset({"imperative", "subjunctive"}),
    ),
    "3sg": (
        frozenset({"indicative", "present", "third-person", "singular"}),
        frozenset({"imperative", "subjunctive"}),
    ),
    "1pl": (
        frozenset({"indicative", "present", "first-person", "plural"}),
        frozenset({"imperative", "subjunctive"}),
    ),
    "sbjv3sg": (
        frozenset({"subjunctive", "present", "third-person", "singular"}),
        frozenset({"imperative", "imperfect"}),
    ),
    "past1sg": (
        frozenset({"historic", "past", "first-person", "singular"}),
        frozenset(),
    ),
}
# Endings tried at the end of a form, in order, per cell.
_ENDINGS: Final[Mapping[str, tuple[str, ...]]] = {
    "inf": ("are", "ere", "ire"),
    "1sg": ("isco", "o"),
    "2sg": ("isci", "i"),
    "3sg": ("isce", "e", "a"),
    "1pl": ("iamo",),
    "sbjv3sg": ("isca", "a", "i"),
    "past1sg": ("qui", "ai", "ei", "etti", "ii", "i"),
}
_ISC_ENDINGS: Final[frozenset[str]] = frozenset(
    {"isco", "isci", "isce", "isca"}
)
_INFINITIVE_RE: Final[re.Pattern[str]] = re.compile(r"(?:ar|er|ir)(?:e|si)\Z")
_ROOT_STRESSED_CATEGORY: Final[str] = (
    "Italian verbs with root-stressed infinitive"
)
# The ending may be followed by a note: "dispiacére [auxiliary essere]".
_STRESSED_ERE_RE: Final[re.Pattern[str]] = re.compile(r"[éè]re\b")
_ACCENTED_RE: Final[re.Pattern[str]] = re.compile(r"[àèéìíòóùú]")
_CLITICS: Final[frozenset[str]] = frozenset({"mi", "ti", "si", "ci", "vi"})

VERB_FIELDS: Final[tuple[str, ...]] = (
    "verb_id",
    "verb",
    "verb_key",
    "reflexive",
    "canonical",
    "conj_class",
    "infinitive_stress",
    "infinitive_stress_source",
    "stem_key",
    "stem_final",
    "target_letters",
    "paradigm_type",
    "overapplies",
    "in_domain",
    "noun_roots",
    "stratum",
    "cell",
    "e1_cell",
    "form",
    "form_key",
    "ending",
    "vowel",
    "edge",
    "palatalized",
    "count_folded",
    "count_cs",
    *(f"attested_{tau}" for tau in attest.TAUS),
    "form_ambiguous",
    "ambiguous_with",
    "source",
)
MAX_AMBIGUOUS: Final[int] = 8


# ---------------------------------------------------------------------------
# Lemmas
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class VerbLemma:
    """A velar-stem verb lemma (module docstring)."""

    verb_id: str
    word: str
    verb_key: str
    reflexive: bool
    canonical: str
    ending: str
    stem_written: str
    stem_key: str
    root_stressed_category: bool
    stratum: str
    ipa_stress: int | None = None
    # Match keys of the Italian words its etymology names (bases).
    named: tuple[str, ...] = ()


def infinitive_parts(key: str) -> tuple[str, str, bool] | None:
    """``(ending, written stem, reflexive)`` of an infinitive key (the
    reflexive -rsi read as -re), or None."""
    if not _INFINITIVE_RE.search(key):
        return None
    reflexive = key.endswith("si")
    plain = key[:-2] + "e" if reflexive else key
    return plain[-3:], plain[:-3], reflexive


def velar_stem(stem_written: str) -> str:
    """The stem with a velar-spelling h removed (bianch → bianc), or ""
    when the stem does not end in a dorsal."""
    stem = stem_written
    if len(stem) >= 2 and stem[-1] == "h" and stem[-2] in "cg":
        stem = stem[:-1]
    if len(stem) < 2 or stem[-1] not in "cg":
        return ""
    return stem


def target_letters(stem_key: str) -> str:
    """sc, cc, gg or the final c/g of a velar stem."""
    if stem_key.endswith(("sc", "cc", "gg")):
        return stem_key[-2:]
    return stem_key[-1]


def ipa_stress(record: Mapping[str, Any]) -> int | None:
    """Stress position (from the end) of the first transcription of the
    infinitive that has one."""
    for sound in record.get("sounds") or ():
        ipa = normalize.parse_ipa(str(sound.get("ipa") or ""))
        if ipa is None:
            continue
        found = normalize.stress_from_end(ipa)
        if found is not None:
            return found
    return None


def load_verbs(path: Path) -> list[VerbLemma]:
    """Velar-stem verb lemmas (and dire-family candidates)."""
    meta, records = store.read_jsonl(path)
    if meta.get("format") != extract.EN_VERBS_FORMAT:
        raise store.StaleDerivativeError(f"{path}: unexpected {meta}")
    out: list[VerbLemma] = []
    for record in records:
        word = str(record.get("word") or "")
        key = keys.match_key(word)
        parts = infinitive_parts(key)
        if parts is None or not key.isalpha():
            continue
        ending, written, reflexive = parts
        stem = velar_stem(written)
        plain = written + ending
        if not stem and not plain.endswith("dire"):
            continue
        canonical = [str(c) for c in record.get("canonical") or ()]
        steps = etymology.steps(record.get("etymology_templates") or ())
        stratum, _ = etymology.stratum(
            steps, str(record.get("etymology_text") or ""), word=word
        )
        named = named_words(steps)
        out.append(
            VerbLemma(
                verb_id=str(record["verb_id"]),
                word=word,
                verb_key=plain,
                reflexive=reflexive,
                canonical=canonical[0] if canonical else "",
                ending=ending,
                stem_written=written,
                stem_key=stem,
                root_stressed_category=_ROOT_STRESSED_CATEGORY
                in (record.get("categories") or ()),
                stratum=stratum,
                ipa_stress=ipa_stress(record),
                named=named,
            )
        )
    out.sort(key=lambda v: v.verb_id)
    return out


_NAMING_KINDS: Final[frozenset[str]] = frozenset(
    {"af", "surf", "from", "der", "uder", "deverbal", "bf"}
)


def named_words(steps: Sequence[etymology.Step]) -> tuple[str, ...]:
    """Match keys of the Italian words an etymology names as parts or
    sources (affixes left out), in order."""
    out: dict[str, None] = {}
    for step in steps:
        if step.kind not in _NAMING_KINDS:
            continue
        for index, term in enumerate(step.terms):
            lang = step.term_lang(index) or step.lang
            if lang not in ("", etymology.PROFILE.code):
                continue
            if not term or term.startswith("-") or term.endswith("-"):
                continue
            out.setdefault(keys.match_key(term), None)
    return tuple(out)


# ---------------------------------------------------------------------------
# Cells
# ---------------------------------------------------------------------------


def cell_of(tags: Iterable[str]) -> str:
    """The cell of a conjugation-table row ("" if none of ours)."""
    found = frozenset(tags)
    if found & features.FILTERED_LABELS:
        return ""
    for cell in CELL_ORDER:
        required, excluded = _CELL_TAGS[cell]
        if required <= found and not found & excluded:
            return cell
    return ""


def strip_clitic(form: str) -> str:
    """The verb form without a leading reflexive clitic (ti amìchi)."""
    parts = form.split()
    if len(parts) == 2 and keys.match_key(parts[0]) in _CLITICS:
        return parts[1]
    return form if len(parts) == 1 else ""


@dataclass(frozen=True, slots=True)
class Parsed:
    """A verb form read at its ending (module docstring)."""

    ending: str
    vowel: str
    edge: str


def _vowel(ending: str) -> str:
    if ending == "qui":
        return "back"
    first = ending[0]
    return first if first in "ie" else "back"


def parse_form(cell: str, form_key: str) -> Parsed:
    """Read *form_key* (a cell of a velar-stem verb) at its ending."""
    for ending in _ENDINGS.get(cell, ()):
        if not form_key.endswith(ending) or len(form_key) <= len(ending):
            continue
        stem = form_key[: -len(ending)]
        vowel = _vowel(ending)
        if ending == "qui":
            edge = VELAR if stem.endswith(("c", "g")) else UNPARSED
        elif stem.endswith(("ch", "gh")):
            edge = VELAR if vowel != "back" else UNPARSED
        elif stem.endswith(("ci", "gi")) and vowel == "back":
            edge = PALATAL
        elif stem.endswith(("c", "g")):
            edge = PALATAL if vowel != "back" else VELAR
        else:
            edge = UNPARSED
        return Parsed(ending, vowel, edge)
    return Parsed("", "", UNPARSED)


@dataclass
class TableForms:
    """``{verb word: {cell: [forms as written]}}`` from the form index."""

    cells: dict[str, dict[str, list[str]]] = field(default_factory=dict)


def load_table_forms(path: Path, words: frozenset[str]) -> TableForms:
    """Conjugation-table forms of *words* (one pass over the index)."""
    meta, fields, rows = store.read_tsv(path)
    if meta.get("format") != extract.FORM_INDEX_FORMAT:
        raise store.StaleDerivativeError(f"{path}: unexpected {meta}")
    column = {name: index for index, name in enumerate(fields)}
    out: dict[str, dict[str, list[str]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for row in rows:
        if row[column["source"]] != "table" or row[column["pos"]] != "verb":
            continue
        lemma = row[column["lemma"]]
        if lemma not in words:
            continue
        cell = cell_of(row[column["tags"]].split("|"))
        form = strip_clitic(row[column["form"]])
        if cell and form:
            out[lemma][cell].append(form)
    return TableForms(
        {
            lemma: {
                cell: sorted(set(forms)) for cell, forms in by_cell.items()
            }
            for lemma, by_cell in out.items()
        }
    )


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------


def infinitive_stress(
    verb: VerbLemma, forms: Mapping[str, Sequence[str]]
) -> tuple[str, str]:
    """``(stress, source)`` of an -ere infinitive: stress ``root``,
    ``ending`` or ``unknown``; source ``ipa``, ``category``,
    ``canonical`` or "" (module docstring). ``("", "")`` for other
    endings."""
    if verb.ending != "ere":
        return "", ""
    if verb.ipa_stress is not None:
        return ("root" if verb.ipa_stress >= 3 else "ending"), "ipa"
    if verb.root_stressed_category:
        return "root", "category"
    spelled = [verb.canonical, *forms.get("inf", ())]
    for text in spelled:
        if _STRESSED_ERE_RE.search(text):
            return "ending", "canonical"
        if _ACCENTED_RE.search(text):
            return "root", "canonical"
    return "unknown", ""


def canonical_stress(verb: VerbLemma) -> str:
    """The stress the canonical accented infinitive implies ("" when it
    has no accent): the cross-check of the IPA reading."""
    if _STRESSED_ERE_RE.search(verb.canonical):
        return "ending"
    if _ACCENTED_RE.search(verb.canonical):
        return "root"
    return ""


def conj_class(
    verb: VerbLemma, forms: Mapping[str, Sequence[str]], stress: str
) -> str:
    """The conjugation class (module docstring)."""
    if verb.verb_key.endswith("dire") and any(
        keys.match_key(f).endswith("dico") for f in forms.get("1sg", ())
    ):
        return DIRE
    if verb.ending == "are":
        return ARE
    if verb.ending == "ere":
        return {"root": ERE_ROOT, "ending": ERE_STRESSED}.get(
            stress, ERE_UNKNOWN
        )
    first = [keys.match_key(f) for f in forms.get("1sg", ())]
    if any(f == verb.stem_written + "isco" for f in first):
        return ISC
    return IRE


def dire_stem(forms: Mapping[str, Sequence[str]]) -> str:
    """The dic- stem of a dire-family verb (from its 1sg)."""
    for form in forms.get("1sg", ()):
        key = keys.match_key(form)
        if key.endswith("dico"):
            return key[:-1]
    return ""


def paradigm_type(parsed: Sequence[tuple[str, Parsed]]) -> tuple[str, bool]:
    """``(paradigm type, overapplies)`` from ``(cell, Parsed)`` pairs."""
    readable = [(cell, p) for cell, p in parsed if p.edge != UNPARSED]
    if not readable:
        return NO_CELLS, False
    front = {p.edge for _, p in readable if p.vowel != "back"}
    velar_anywhere = any(p.edge == VELAR for _, p in readable)
    first = [p for cell, p in readable if cell == "1sg" and p.vowel == "back"]
    overapplies = bool(first) and all(p.edge == PALATAL for p in first)
    if not velar_anywhere:
        return PALATAL_THROUGHOUT, overapplies
    if overapplies:
        return OVERAPPLYING, True
    if front == {VELAR}:
        return VELAR_THROUGHOUT, False
    if front == {PALATAL}:
        return ALTERNATING, False
    if front:
        return MIXED, False
    return VELAR_THROUGHOUT, False


# ---------------------------------------------------------------------------
# Rows
# ---------------------------------------------------------------------------


@dataclass
class RootLink:
    """What the E1 link needs of a dorsal noun/adjective root."""

    root_id: str
    words: set[str] = field(default_factory=set)
    deverbal_from: set[str] = field(default_factory=set)


def dorsal_roots_by_stem(path: Path) -> dict[str, list[RootLink]]:
    """Dorsal noun/adjective roots by stem key (lexicon_complete), with
    their lemmas and the verbs their lexemes are deverbal from."""
    found: dict[str, dict[str, RootLink]] = defaultdict(dict)
    for row in csvio.read_rows(path):
        if row["in_domain"] != "True" or row["dorsal_target"] != "True":
            continue
        stem = morphology.analyze_lemma(keys.match_key(row["word"])).stem
        link = found[stem].setdefault(row["root_id"], RootLink(row["root_id"]))
        link.words.add(keys.match_key(row["word"]))
        if row["deverbal_verb"]:
            link.deverbal_from.add(keys.match_key(row["deverbal_verb"]))
    return {
        stem: [links[root] for root in sorted(links)]
        for stem, links in found.items()
    }


def linked_roots(verb: VerbLemma, candidates: Sequence[RootLink]) -> list[str]:
    """The roots a verb's cells are E1 evidence for (module docstring):
    every candidate, or — when several share the stem — those with
    etymological evidence, if any has some."""
    if len(candidates) < 2:
        return [link.root_id for link in candidates]
    verb_words = {keys.match_key(verb.word), verb.verb_key}
    named = set(verb.named)
    supported = [
        link.root_id
        for link in candidates
        if link.deverbal_from & verb_words or link.words & named
    ]
    return supported or [link.root_id for link in candidates]


def lemmas_of_forms(path: Path, wanted: frozenset[str]) -> dict[str, set[str]]:
    """``{form key: {"lemma:pos", …}}`` for the wanted form keys."""
    meta, fields, rows = store.read_tsv(path)
    if meta.get("format") != extract.FORM_INDEX_FORMAT:
        raise store.StaleDerivativeError(f"{path}: unexpected {meta}")
    column = {name: index for index, name in enumerate(fields)}
    out: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        key = row[column["key"]]
        if key in wanted:
            out[key].add(f"{row[column['lemma']]}:{row[column['pos']]}")
    return dict(out)


@dataclass(frozen=True, slots=True)
class VerbResult:
    """A classified verb and its parsed cells."""

    verb: VerbLemma
    stem_key: str
    conj: str
    stress: str
    stress_source: str
    kind: str
    overapplies: bool
    cells: tuple[tuple[str, str, Parsed], ...]  # (cell, form, parsed)


def classify(
    verb: VerbLemma, forms: Mapping[str, Sequence[str]]
) -> VerbResult:
    """Class, paradigm type and parsed cells of one verb."""
    stress, stress_source = infinitive_stress(verb, forms)
    conj = conj_class(verb, forms, stress)
    stem_key = verb.stem_key or (dire_stem(forms) if conj == DIRE else "")
    cells: list[tuple[str, str, Parsed]] = []
    for cell in CELL_ORDER:
        for form in forms.get(cell, ()):
            parsed = parse_form(cell, keys.match_key(form))
            if conj != ISC and parsed.ending in _ISC_ENDINGS:
                parsed = Parsed("", "", UNPARSED)
            cells.append((cell, form, parsed))
    kind, overapplies = paradigm_type([(c, p) for c, _, p in cells])
    return VerbResult(
        verb,
        stem_key,
        conj,
        stress,
        stress_source,
        kind,
        overapplies,
        tuple(cells),
    )


def verb_rows(
    result: VerbResult,
    roots_by_stem: Mapping[str, Sequence[RootLink]],
    table: AttestationTable,
    form_lemmas: Mapping[str, set[str]],
) -> Iterator[dict[str, Any]]:
    """The :data:`VERB_FIELDS` rows of one verb."""
    verb = result.verb
    own = {f"{verb.word}:verb"}
    common: dict[str, Any] = {
        "verb_id": verb.verb_id,
        "verb": verb.word,
        "verb_key": verb.verb_key,
        "reflexive": verb.reflexive,
        "canonical": verb.canonical,
        "conj_class": result.conj,
        "infinitive_stress": result.stress,
        "infinitive_stress_source": result.stress_source,
        "stem_key": result.stem_key,
        "stem_final": result.stem_key[-1:] if result.stem_key else "",
        "target_letters": (
            target_letters(result.stem_key) if result.stem_key else ""
        ),
        "paradigm_type": result.kind,
        "overapplies": result.overapplies,
        "in_domain": result.kind not in OUT_OF_DOMAIN
        and bool(result.stem_key),
        "noun_roots": linked_roots(
            verb, roots_by_stem.get(result.stem_key, ())
        ),
        "stratum": verb.stratum,
        "source": SOURCE,
    }
    for cell, form, parsed in result.cells:
        key = keys.match_key(form)
        others = sorted(form_lemmas.get(key, set()) - own)
        row = dict(common)
        row.update(
            {
                "cell": cell,
                "e1_cell": cell in E1_CELLS,
                "form": form,
                "form_key": key,
                "ending": parsed.ending,
                "vowel": parsed.vowel,
                "edge": parsed.edge,
                "palatalized": (
                    None
                    if parsed.edge == UNPARSED or parsed.vowel == "back"
                    else parsed.edge == PALATAL
                ),
                "count_folded": table.count(key),
                "count_cs": table.count(key, case_sensitive=True),
                "form_ambiguous": bool(others),
                "ambiguous_with": others[:MAX_AMBIGUOUS],
            }
        )
        for tau, flag in leipzig.attested_flags(table, key).items():
            row[f"attested_{tau}"] = flag
        yield row


def summary(results: Sequence[VerbResult]) -> str:
    """Counts by class and paradigm type."""
    counts: Counter[str] = Counter()
    for result in results:
        counts["verbs"] += 1
        counts[f"conj:{result.conj}"] += 1
        counts[f"conj_type:{result.conj}:{result.kind}"] += 1
        counts[f"paradigm_type:{result.kind}"] += 1
        if result.overapplies:
            counts[f"overapplies:{result.conj}"] += 1
        if result.stress_source:
            counts[f"stress_source:{result.stress_source}"] += 1
        check = canonical_stress(result.verb)
        if result.stress_source == "ipa" and check:
            agree = "agree" if check == result.stress else "disagree"
            counts[f"stress_ipa_vs_canonical:{agree}"] += 1
    lines = ["key\tvalue"]
    lines.extend(f"{key}\t{counts[key]}" for key in sorted(counts))
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# The stage
# ---------------------------------------------------------------------------


def _row_key(row: Mapping[str, Any]) -> tuple[str, str, int, str]:
    return (
        str(row["verb_key"]),
        str(row["verb_id"]),
        CELL_ORDER.index(str(row["cell"])),
        str(row["form"]),
    )


def _body(ctx: RunContext, layout: Layout, words: Path) -> int:
    verbs = load_verbs(layout.en_verbs)
    tables = load_table_forms(
        layout.form_index, frozenset(v.word for v in verbs)
    )
    results = [classify(v, tables.cells.get(v.word, {})) for v in verbs]
    results = [r for r in results if r.stem_key]
    wanted = frozenset(
        keys.match_key(form) for r in results for _, form, _ in r.cells
    )
    form_lemmas = lemmas_of_forms(layout.form_index, wanted)
    roots_by_stem = dorsal_roots_by_stem(layout.lexicon_complete)
    table = leipzig.attestation_table(ctx, words)
    rows = [
        row
        for result in results
        for row in verb_rows(result, roots_by_stem, table, form_lemmas)
    ]
    rows.sort(key=_row_key)
    csvio.write_rows(layout.verbs, rows, VERB_FIELDS, ctx=ctx)
    ctx.write_text(layout.verbs_summary, summary(results))
    in_domain = sum(1 for r in results if r.kind not in OUT_OF_DOMAIN)
    print(
        f"  {len(results):,} velar-stem verbs ({in_domain:,} in domain), "
        f"{len(rows):,} cell rows"
    )
    return 0


def outputs(layout: Layout) -> list[Path]:
    return [layout.verbs, layout.verbs_summary]


def run(ctx: RunContext) -> int | None:
    """Run the ``verbs`` stage (skipped when its manifest matches)."""
    layout = Layout(ctx)
    words = leipzig.words_file(ctx, layout)
    return run_stage(
        ctx,
        STAGE,
        inputs=[
            layout.en_verbs,
            layout.form_index,
            layout.lexicon_complete,
            words,
        ],
        params={
            "version": VERSION,
            "cells": list(CELL_ORDER),
            "filtered_labels": sorted(features.FILTERED_LABELS),
            "elision_clitics": sorted(ctx.language.elision_clitics),
        },
        outputs=outputs(layout),
        body=lambda: _body(ctx, layout, words),
    )
