"""Stage ``features``: cells and lexeme summary of Italian (design §5).

Inputs: ``italian_lexicon_full.csv`` and ``italian_plural_forms.csv``
(``lexicon``), ``italian_qc.csv`` (``qc``), ``italian_attestation.csv``
(``attest``), the ``extract`` projections (etymology templates and
glosses of the lexemes, the verb lemmas, the form index's lemmas) and
the curated files ``curated/italian/plural_adjudications.csv`` and
``curated/italian/suffixoid_classes.csv`` (loaded with
:mod:`lexpipe.curated`; a row whose guard columns no longer match the
data stops the stage, R4/R6).

``data/italian/italian_cells.csv``
----------------------------------
One row per (lexeme, plural cell, form), sorted by lexeme_id, cell,
item. Per form (:mod:`.morphology`): ``opportunity`` (i / e / a /
none), ``mutation`` (palatalized True, velar kept False; empty when
there is no front-vowel opportunity), ``orth_change``,
``palatal_output``, ``aligned``; per lexeme: ``stem_final`` (c/g; t/d/s/z
for the "no coronal rule" rows), ``cluster``, ``target`` (dorsal |
coronal | ""), ``frontstem``, ``nde_class`` (occhio | psiche |
gimpe_like), stress and suffix columns; ``filtered`` marks archaic, Old
Italian, literary, poetic, dialectal, nonstandard and proscribed
variants (:data:`FILTERED_LABELS`), which are reported but left out of
the cell outcome (§5); ``exception_reason``:

    undergoer       dorsal target, front-vowel opportunity, palatalized
    nde:<class>     an NDE stem (occhio, psiche) or an invariable
                    C + front vowel (gimpe_like)
    unexplained     dorsal target before a front vowel, velar kept
                    (the e cell included: amiche is a non-trigger N = e
                    row, R16), or a coronal before -i (no coronal rule)
    non exception   everything else (back-vowel or no opportunity,
                    frontstem, coronal before -e, other finals)

``cell_outcome`` (palatalized | velar | variable | unchanged | none) is
computed over the cell's unfiltered, aligned forms; for a variable cell
``doublet_class`` is :func:`lexpipe.tp.classify_doublet` of the forms'
case-sensitive corpus counts (R19: ``dominant:<form>`` at ≥ 90% of at
least :data:`DOUBLET_MIN_TOKENS` tokens, else ``true_doublet``: a
variant with 2 tokens against 0 is not dominant). Corpus counts, τ
flags and the competing spelling's counts (``competitor``) are copied
from the ``attest`` stage; ``plural_qc`` from the ``qc`` stage.

``data/italian/italian_lexicon_complete.csv``
---------------------------------------------
One row per lexeme, sorted by lexeme_id: identity, domain and target
flags, stress (``velar_post_tonic`` = the stem-final dorsal is in the
syllable right after the stressed one; ``ico_stressed`` = the i of a
final -ico/-ica carries the stress, amìco), suffix status (R26:
``lemma_suffix``, ``target_is_suffix``, ``suffix_evidence``,
``suffix_base``, ``target_is_suffix_strict`` = with an attested base),
the -logo/-fago/-urgo columns of R27 — ``suffixoid`` (the ending, with
evidence: :func:`.etymology.suffixoid_evidence`), ``learned_suffixoid``
and ``suffixoid_basis`` (the shared detector of
:mod:`lexpipe.suffixoids`; a curated ``False`` marks a lookalike such
as spurgo, whose ``lemma_suffix`` then stays empty), ``semantic_class``
and ``class_basis`` (the shared person/thing classifier on the first
gloss; adjectives are ``adjective``; a curated class wins),
``deverbal_from_are`` (etymology: "Deverbal from X-are", "short past
participle of X-are") kept apart from ``shares_root_with_1st_conj_verb``
(synchronic: an -are/-arsi verb lemma on the same stem exists), stratum
(R22), the plural summary and ``plural_qc`` (worst cell: conflict >
variant > module_default_only > confirmed), the outcomes of the i and e
opportunities, ``exception_reason`` of the i opportunity (its dominant
or first unfiltered form), ``doublet_class``, curated adjudications,
and frequencies (R24): lemma form and paradigm sum over the lexeme's
distinct attested forms, case-sensitive and folded, plus the root's
folded paradigm sum.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

from lexpipe import attest, csvio, curated, suffixoids, tp
from lexpipe.context import RunContext, run_stage
from lexpipe.langs.italian import (
    etymology,
    extract,
    keys,
    morphology,
    qc,
    store,
)
from lexpipe.langs.italian.layout import Layout

STAGE: Final[str] = "features"
# 2: etymology.stratum keeps a leading affixation (antichità).
# 3: the shared etymology classifier and -logo classifiers
#    (lexpipe.etymology, lexpipe.suffixoids); reads the full lexicon.
# 4: R26 from the stressed vowel (arcaico), suffix/confix templates and
#    verb-only back-formations; nde:psiche needs a velar singular;
#    doublets need DOUBLET_MIN_TOKENS for dominance; competitor counts.
VERSION: Final[int] = 4
# A doublet (R19) is dominant only when its variants have at least this
# many case-sensitive tokens in all (sarcofagi 2 vs sarcofaghi 0 is not).
DOUBLET_MIN_TOKENS: Final[int] = 10

FILTERED_LABELS: Final[frozenset[str]] = frozenset(
    {
        "archaic",
        "obsolete",
        "old italian",
        "literary",
        "poetic",
        "dialectal",
        "nonstandard",
        "proscribed",
    }
)
ADJECTIVE_CLASS: Final[str] = "adjective"
UNDERGOER: Final[str] = "undergoer"
UNEXPLAINED: Final[str] = "unexplained"
NON_EXCEPTION: Final[str] = "non exception"
PALATALIZED: Final[str] = "palatalized"
VELAR: Final[str] = "velar"
VARIABLE: Final[str] = "variable"
UNCHANGED: Final[str] = "unchanged"
NONE: Final[str] = "none"
CELL_ORDER: Final[Mapping[str, int]] = {"pl": 0, "m.pl": 1, "f.pl": 2}
FRONT: Final[frozenset[str]] = frozenset({"i", "e"})
SOURCE: Final[str] = "wiktionary_en"

CELL_FIELDS: Final[tuple[str, ...]] = (
    "lexeme_id",
    "root_id",
    "word",
    "pos",
    "gender",
    "cell",
    "form",
    "item",
    "provenance",
    "code",
    "labels",
    "filtered",
    "opportunity",
    "stem_final",
    "cluster",
    "target",
    "frontstem",
    "nde_class",
    "mutation",
    "orth_change",
    "palatal_output",
    "aligned",
    "exception_reason",
    "cell_outcome",
    "doublet_class",
    "adjudicated_out",
    "stress_from_end",
    "stress_class",
    "velar_post_tonic",
    "lemma_suffix",
    "target_is_suffix",
    "stratum",
    "plural_qc",
    "count_cs",
    "count_folded",
    *(f"attested_{tau}" for tau in attest.TAUS),
    "form_ambiguous",
    "form_ambiguous_cs",
    "competitor",
    "competitor_count_cs",
    "competitor_other_lemma",
    "in_domain",
    "source",
)
LEXEME_FIELDS: Final[tuple[str, ...]] = (
    "lexeme_id",
    "root_id",
    "word",
    "pos",
    "gender",
    "etymology_number",
    "first_gloss",
    "ipa_raw",
    "domain",
    "domain_flags",
    "in_domain",
    "stem_final",
    "cluster",
    "target",
    "dorsal_target",
    "frontstem",
    "nde_class",
    "ndeb_internal",
    "stress_from_end",
    "stress_class",
    "stress_source",
    "stress_missing",
    "stress_variants",
    "velar_post_tonic",
    "ico_stressed",
    "lemma_suffix",
    "suffix_kind",
    "target_is_suffix",
    "target_is_suffix_strict",
    "suffix_evidence",
    "suffix_raw",
    "suffix_base",
    "suffixoid",
    "learned_suffixoid",
    "suffixoid_basis",
    "semantic_class",
    "class_basis",
    "deverbal_from_are",
    "deverbal_verb",
    "shares_root_with_1st_conj_verb",
    "first_conj_verb",
    "stratum",
    "stratum_basis",
    "etymon_lang",
    "etymon",
    "plural_status",
    "plural_provenance",
    "plural_qc",
    "countability",
    "invariable",
    "plurale_tantum",
    "opportunities",
    "forms_i",
    "outcome_i",
    "forms_e",
    "outcome_e",
    "filtered_forms",
    "mutation_i",
    "orth_change_i",
    "palatal_output_i",
    "exception_reason",
    "doublet_class",
    "adjudicated_plurals",
    "adjudication",
    "freq_lemma_cs",
    "freq_lemma_folded",
    "freq_paradigm_cs",
    "freq_paradigm_folded",
    "freq_root_folded",
    "lemma_per_million",
    "source",
)


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Details:
    """Projection fields of one lexeme that the lexicon does not carry."""

    steps: tuple[etymology.Step, ...]
    etymology_text: str
    glosses: tuple[str, ...]


def load_details(path: Path) -> dict[str, Details]:
    """Etymology steps, etymology text and sense glosses by lexeme_id."""
    meta, records = store.read_jsonl(path)
    if meta.get("format") != extract.EN_LEXEMES_FORMAT:
        raise store.StaleDerivativeError(f"{path}: unexpected format {meta}")
    out: dict[str, Details] = {}
    for record in records:
        glosses = tuple(
            str(sense["glosses"][0])
            for sense in record.get("senses") or ()
            if sense.get("glosses")
            and not sense.get("form_of")
            and not sense.get("alt_of")
        )
        out[str(record["lexeme_id"])] = Details(
            tuple(etymology.steps(record.get("etymology_templates") or ())),
            str(record.get("etymology_text") or ""),
            glosses,
        )
    return out


def load_verb_keys(path: Path) -> frozenset[str]:
    """Match keys of the verb lemmas (``italian_en_verbs.jsonl.gz``)."""
    meta, records = store.read_jsonl(path)
    if meta.get("format") != extract.EN_VERBS_FORMAT:
        raise store.StaleDerivativeError(f"{path}: unexpected format {meta}")
    return frozenset(keys.match_key(str(r["word"])) for r in records)


def load_lemma_keys(path: Path) -> frozenset[str]:
    """Match keys of every lemma of the en split (any part of speech)."""
    meta, fields, rows = store.read_tsv(path)
    if meta.get("format") != extract.FORM_INDEX_FORMAT:
        raise store.StaleDerivativeError(f"{path}: unexpected format {meta}")
    column = {name: index for index, name in enumerate(fields)}
    key, source = column["key"], column["source"]
    return frozenset(row[key] for row in rows if row[source] == "lemma")


@dataclass(frozen=True, slots=True)
class Counts:
    """Attestation of one (lexeme, cell, form)."""

    count_cs: int
    count_folded: int
    per_million: str
    attested: Mapping[str, str]
    ambiguous: str
    ambiguous_cs: str
    competitor: str = ""
    competitor_count_cs: str = ""
    competitor_other_lemma: str = ""


def load_attestation(path: Path) -> dict[tuple[str, str, str], Counts]:
    out: dict[tuple[str, str, str], Counts] = {}
    for row in csvio.read_rows(path):
        out[(row["lexeme_id"], row["cell"], row["form"])] = Counts(
            int(row["count_cs"]),
            int(row["count_folded"]),
            row["per_million"],
            {tau: row[f"attested_{tau}"] for tau in attest.TAUS},
            row["form_ambiguous"],
            row["form_ambiguous_cs"],
            row["competitor"],
            row["competitor_count_cs"],
            row["competitor_other_lemma"],
        )
    return out


def load_plural_qc(path: Path) -> dict[tuple[str, str], str]:
    return {
        (row["lexeme_id"], row["cell"]): row["plural_qc"]
        for row in csvio.read_rows(path)
    }


def group_plurals(
    rows: Iterable[Mapping[str, str]],
) -> dict[str, list[Mapping[str, str]]]:
    out: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for row in rows:
        if row["form"]:
            out[row["lexeme_id"]].append(row)
    return dict(out)


# ---------------------------------------------------------------------------
# Curated decisions
# ---------------------------------------------------------------------------


def _current_for_adjudications(
    adjudications: Mapping[str, curated.CuratedRow],
    lexicon: Mapping[str, Mapping[str, str]],
    plurals: Mapping[str, Sequence[Mapping[str, str]]],
) -> dict[str, dict[str, Any]]:
    current: dict[str, dict[str, Any]] = {}
    for lexeme_id, decision in adjudications.items():
        row = lexicon.get(lexeme_id)
        if row is None:
            continue
        cell = decision.values.get("cell", "")
        forms = [
            p["form"]
            for p in plurals.get(lexeme_id, ())
            if not cell or p["cell"] == cell
        ]
        current[lexeme_id] = {
            "word": row["word"],
            "pos": row["pos"],
            "gender": row["gender"],
            "plurals": list(dict.fromkeys(forms)),
        }
    return current


def check_curated(
    layout: Layout,
    lexicon: Mapping[str, Mapping[str, str]],
    plurals: Mapping[str, Sequence[Mapping[str, str]]],
) -> tuple[dict[str, curated.CuratedRow], dict[str, suffixoids.Decision]]:
    """Load both curated files and stop on any guard mismatch."""
    adjudications = curated.load_curated(layout.plural_adjudications)
    curated.check_guards(
        adjudications,
        _current_for_adjudications(adjudications, lexicon, plurals),
    )
    rows, decisions = suffixoids.load_decisions(layout.suffixoid_classes)
    curated.check_guards(
        rows,
        {
            lexeme_id: {
                "word": lexicon[lexeme_id]["word"],
                "pos": lexicon[lexeme_id]["pos"],
                "first_gloss": lexicon[lexeme_id]["first_gloss"],
            }
            for lexeme_id in rows
            if lexeme_id in lexicon
        },
    )
    return adjudications, decisions


# ---------------------------------------------------------------------------
# Per-lexeme features
# ---------------------------------------------------------------------------


@dataclass
class Context:
    """Everything shared by the per-lexeme computations."""

    details: Mapping[str, Details]
    verbs: frozenset[str]
    lemmas: frozenset[str]
    attestation: Mapping[tuple[str, str, str], Counts]
    plural_qc: Mapping[tuple[str, str], str]
    adjectival_roots: frozenset[str]
    adjudications: Mapping[str, curated.CuratedRow]
    decisions: Mapping[str, suffixoids.Decision]
    root_frequency: dict[str, int] = field(default_factory=dict)


def _bool(text: str) -> bool:
    return csvio.parse_bool(text) is True


def first_conjugation_verb(key: str, verbs: frozenset[str]) -> str:
    """An -are/-arsi verb lemma on the lemma's stem ("" if none)."""
    stem = keys.strip_final_vowel(key)
    if stem == key:
        return ""
    for ending in ("are", "arsi"):
        if stem + ending in verbs:
            return stem + ending
    return ""


def suffix_columns(
    row: Mapping[str, str], details: Details | None, ctx: Context
) -> dict[str, Any]:
    """R26 suffix status and the R27 suffixoid of one lexeme."""
    key = keys.match_key(row["word"])
    steps = details.steps if details else ()
    text = details.etymology_text if details else ""
    verb = etymology.deverbal_verb(steps, text)
    evidence = etymology.ico_suffix(
        key,
        steps,
        ico_stressed=csvio.parse_bool(row["ico_stressed"]),
        adjectival=row["pos"] == "adj"
        or row["root_id"] in ctx.adjectival_roots,
        deverbal=bool(verb),
        italian_lemmas=ctx.lemmas,
    )
    candidate, ending_evidence = (
        ("", "")
        if evidence.lemma_suffix
        else etymology.suffixoid_evidence(
            key,
            steps,
            deverbal=bool(verb),
            masculine_forms=csvio.split_list(row["masculine_forms"]),
        )
    )
    learned, learned_basis = (
        suffixoids.learned_suffixoid(
            key,
            etymology.SUFFIXOID_KEYS[candidate],
            "; ".join(details.glosses if details else ()),
            decision=ctx.decisions.get(row["lexeme_id"]),
            evidence=ending_evidence,
        )
        if candidate
        else (False, "")
    )
    suffixoid = candidate if learned else ""
    lemma_suffix = evidence.lemma_suffix or suffixoid
    kind = (
        "suffix"
        if evidence.lemma_suffix
        else ("suffixoid" if suffixoid else "")
    )
    return {
        "lemma_suffix": lemma_suffix,
        "suffix_kind": kind,
        "target_is_suffix": evidence.target_is_suffix,
        "target_is_suffix_strict": evidence.target_is_suffix
        and bool(evidence.base),
        "suffix_evidence": evidence.evidence,
        "suffix_raw": evidence.raw_suffix,
        "suffix_base": evidence.base,
        "deverbal_verb": verb,
        "deverbal_from_are": bool(verb)
        and etymology.is_first_conjugation(verb),
        "suffixoid": candidate,
        "learned_suffixoid": learned if candidate else None,
        "suffixoid_basis": learned_basis,
    }


def semantic_class_columns(
    row: Mapping[str, str], candidate: str, ctx: Context
) -> dict[str, str]:
    """R27 class of a -logo/-fago/-urgo lexeme: ``adjective`` for an
    adjective, else the shared person/thing classifier on the first
    gloss (:func:`lexpipe.suffixoids.semantic_class`; a curated class
    wins). Empty for other lexemes."""
    if candidate not in etymology.SUFFIXOIDS:
        return {"semantic_class": "", "class_basis": ""}
    decision = ctx.decisions.get(row["lexeme_id"])
    if row["pos"] == "adj" and not (decision and decision.semantic_class):
        return {"semantic_class": ADJECTIVE_CLASS, "class_basis": "pos"}
    klass, basis = suffixoids.semantic_class(
        row["first_gloss"], decision=decision
    )
    return {"semantic_class": klass, "class_basis": basis}


@dataclass(frozen=True, slots=True)
class FormFeatures:
    """One plural form with its analysis."""

    plural: Mapping[str, str]
    analysis: morphology.FormAnalysis
    filtered: bool
    adjudicated_out: bool
    reason: str


def _reason(
    stem: morphology.StemAnalysis,
    analysis: morphology.FormAnalysis,
    gimpe: bool,
) -> str:
    if gimpe:
        return f"nde:{morphology.NDE_GIMPE}"
    if analysis.opportunity not in FRONT:
        return NON_EXCEPTION
    if stem.nde_class:
        return f"nde:{stem.nde_class}"
    if stem.dorsal and analysis.mutation is not None:
        return UNDERGOER if analysis.mutation else UNEXPLAINED
    if stem.target == morphology.CORONAL and analysis.opportunity == "i":
        return UNEXPLAINED if analysis.aligned else NON_EXCEPTION
    return NON_EXCEPTION


def form_features(
    row: Mapping[str, str],
    stem: morphology.StemAnalysis,
    plurals: Sequence[Mapping[str, str]],
    adjudication: curated.CuratedRow | None,
) -> list[FormFeatures]:
    lemma_key = keys.match_key(row["word"])
    allowed = (
        set(
            csvio.split_list(
                adjudication.values.get("adjudicated_plurals", "")
            )
        )
        if adjudication
        else set()
    )
    adjudicated_cell = (
        adjudication.values.get("cell", "") if adjudication else ""
    )
    out: list[FormFeatures] = []
    for plural in plurals:
        form_key = keys.match_key(plural["form"])
        analysis = morphology.analyze_form(lemma_key, stem, form_key)
        gimpe = morphology.gimpe_like(lemma_key, form_key)
        labels = set(csvio.split_list(plural["labels"]))
        out_by_curator = (
            bool(allowed)
            and (not adjudicated_cell or plural["cell"] == adjudicated_cell)
            and plural["form"] not in allowed
        )
        out.append(
            FormFeatures(
                plural,
                analysis,
                bool(labels & FILTERED_LABELS),
                out_by_curator,
                _reason(stem, analysis, gimpe),
            )
        )
    return out


def cell_outcome(
    stem: morphology.StemAnalysis, forms: Sequence[FormFeatures]
) -> str:
    """Outcome of the unfiltered, aligned front-vowel forms."""
    usable = [
        f
        for f in forms
        if not f.filtered
        and not f.adjudicated_out
        and f.analysis.aligned
        and f.analysis.opportunity in FRONT
    ]
    if not usable:
        return NONE
    if not stem.dorsal or stem.nde_class:
        return UNCHANGED
    outcomes = {f.analysis.mutation for f in usable}
    if outcomes == {True}:
        return PALATALIZED
    if outcomes == {False}:
        return VELAR
    return VARIABLE


def doublet_class(
    lexeme_id: str,
    forms: Sequence[FormFeatures],
    attestation: Mapping[tuple[str, str, str], Counts],
) -> str:
    """:func:`lexpipe.tp.classify_doublet` of a variable cell."""
    counts = {
        f.plural["form"]: attestation[
            (lexeme_id, f.plural["cell"], f.plural["form"])
        ].count_cs
        for f in forms
        if not f.filtered and f.analysis.mutation is not None
    }
    if len(counts) < 2:
        return ""
    return tp.classify_doublet(
        dict(sorted(counts.items())), min_total=DOUBLET_MIN_TOKENS
    )


@dataclass
class LexemeResult:
    """The complete row and the cell rows of one lexeme."""

    row: dict[str, Any]
    cells: list[dict[str, Any]]


# The consonant before the final vowel of a psiche-type singular: a velar.
_VELAR_IPA: Final[frozenset[str]] = frozenset({"k", "ɡ", "g"})
_IPA_MARKS: Final[str] = "ˈˌ.ː"


def velar_che(row: Mapping[str, str]) -> bool:
    """Is a singular in -che/-ghe pronounced with a velar (psiche
    /ˈpsi.ke/), as the ``nde:psiche`` class needs? From the IPA when
    there is one (brioche /briˈɔʃ/, ceviche /seˈvi.t͡ʃe/ are not); else
    from a velar plural (nuraghe → nuraghi). French loans without
    either (cloche, peluche) are left out."""
    ipa = row["ipa"].rstrip(_IPA_MARKS)
    if ipa:
        body = ipa.rstrip("eɛ")
        if body == ipa:
            return False  # the transcription does not end in -e
        consonant = body.rstrip(_IPA_MARKS)
        return bool(consonant) and consonant[-1] in _VELAR_IPA
    forms = [
        keys.match_key(form)
        for column in ("pl_forms", "mpl_forms", "fpl_forms")
        for form in csvio.split_list(row[column])
    ]
    return any(form.endswith(("chi", "ghi")) for form in forms)


def _stem(row: Mapping[str, str]) -> morphology.StemAnalysis:
    if row["domain"] != "in_domain":
        return morphology.NO_STEM
    found = morphology.analyze_lemma(keys.match_key(row["word"]))
    if found.nde_class == morphology.NDE_PSICHE and not velar_che(row):
        return morphology.NO_STEM
    return found


def _primary(
    forms: Sequence[FormFeatures], doublet: str
) -> FormFeatures | None:
    """The dominant (or first unfiltered) front-vowel form."""
    usable = [
        f
        for f in forms
        if not f.filtered
        and not f.adjudicated_out
        and f.analysis.opportunity in FRONT
    ]
    if doublet.startswith(tp.DOMINANT_PREFIX):
        dominant = doublet[len(tp.DOMINANT_PREFIX) :]
        for f in usable:
            if f.plural["form"] == dominant:
                return f
    return usable[0] if usable else None


def _frequencies(
    row: Mapping[str, str],
    attestation: Mapping[tuple[str, str, str], Counts],
    by_lexeme: Mapping[str, Sequence[tuple[str, str]]],
) -> dict[str, Any]:
    """Lemma-form counts and paradigm sums over the lexeme's distinct
    forms (R24); a form filling two cells counts once."""
    lemma = attestation.get((row["lexeme_id"], "lemma", row["word"]))
    folded: dict[str, int] = {}
    cased: dict[str, int] = {}
    for cell, form in by_lexeme.get(row["lexeme_id"], ()):
        found = attestation.get((row["lexeme_id"], cell, form))
        if found is not None:
            folded.setdefault(keys.match_key(form), found.count_folded)
            cased.setdefault(keys.case_key(form), found.count_cs)
    return {
        "freq_lemma_cs": lemma.count_cs if lemma else 0,
        "freq_lemma_folded": lemma.count_folded if lemma else 0,
        "freq_paradigm_cs": sum(cased.values()),
        "freq_paradigm_folded": sum(folded.values()),
        "lemma_per_million": lemma.per_million if lemma else "",
    }


def _rollup(lexeme_id: str, cells: Iterable[str], ctx: Context) -> str:
    found = {ctx.plural_qc.get((lexeme_id, cell), "") for cell in cells}
    for value in qc.ROLLUP_ORDER:
        if value in found:
            return value
    return ""


def lexeme_result(
    row: Mapping[str, str],
    plurals: Sequence[Mapping[str, str]],
    ctx: Context,
    forms_by_lexeme: Mapping[str, Sequence[tuple[str, str]]],
) -> LexemeResult:
    """Features of one lexeme (see the module docstring)."""
    lexeme_id = row["lexeme_id"]
    details = ctx.details.get(lexeme_id)
    stem = _stem(row)
    in_domain = row["domain"] == "in_domain" and row["pos"] in ("noun", "adj")
    stratum, basis = etymology.stratum(
        details.steps if details else (),
        details.etymology_text if details else "",
        word=row["word"],
    )
    suffix = suffix_columns(row, details, ctx)
    adjudication = ctx.adjudications.get(lexeme_id)
    forms = form_features(row, stem, plurals, adjudication)
    by_cell: dict[str, list[FormFeatures]] = defaultdict(list)
    for f in forms:
        by_cell[f.plural["cell"]].append(f)
    outcomes = {cell: cell_outcome(stem, fs) for cell, fs in by_cell.items()}
    doublets = {
        cell: (
            doublet_class(lexeme_id, fs, ctx.attestation)
            if outcomes[cell] == VARIABLE
            else ""
        )
        for cell, fs in by_cell.items()
    }
    post_tonic = (
        row["final_post_tonic"]
        if stem.dorsal and row["final_post_tonic"]
        else ""
    )
    common = {
        "stem_final": stem.stem_final,
        "cluster": stem.cluster,
        "target": stem.target,
        "frontstem": stem.frontstem,
        "nde_class": stem.nde_class,
        "stress_from_end": row["stress_from_end"],
        "stress_class": row["stress_class"],
        "velar_post_tonic": post_tonic,
        "lemma_suffix": suffix["lemma_suffix"],
        "target_is_suffix": suffix["target_is_suffix"],
        "stratum": stratum,
    }
    cells = [
        _cell_row(row, f, common, outcomes, doublets, ctx, in_domain)
        for f in forms
    ]
    complete = _complete_row(
        row,
        stem,
        forms,
        outcomes,
        doublets,
        common,
        suffix,
        details,
        ctx,
        forms_by_lexeme,
        (stratum, basis),
        in_domain,
    )
    return LexemeResult(complete, cells)


def _cell_row(
    row: Mapping[str, str],
    f: FormFeatures,
    common: Mapping[str, Any],
    outcomes: Mapping[str, str],
    doublets: Mapping[str, str],
    ctx: Context,
    in_domain: bool,
) -> dict[str, Any]:
    plural = f.plural
    cell = plural["cell"]
    counts = ctx.attestation.get((row["lexeme_id"], cell, plural["form"]))
    out: dict[str, Any] = {
        "lexeme_id": row["lexeme_id"],
        "root_id": row["root_id"],
        "word": row["word"],
        "pos": row["pos"],
        "gender": row["gender"],
        "cell": cell,
        "form": plural["form"],
        "item": int(plural["item"]),
        "provenance": plural["provenance"],
        "code": plural["code"],
        "labels": csvio.split_list(plural["labels"]),
        "filtered": f.filtered,
        "opportunity": f.analysis.opportunity,
        "mutation": f.analysis.mutation,
        "orth_change": f.analysis.orth_change,
        "palatal_output": f.analysis.palatal_output,
        "aligned": f.analysis.aligned,
        "exception_reason": f.reason,
        "cell_outcome": outcomes[cell],
        "doublet_class": doublets[cell],
        "adjudicated_out": f.adjudicated_out,
        "plural_qc": ctx.plural_qc.get((row["lexeme_id"], cell), ""),
        "count_cs": counts.count_cs if counts else None,
        "count_folded": counts.count_folded if counts else None,
        "form_ambiguous": counts.ambiguous if counts else None,
        "form_ambiguous_cs": counts.ambiguous_cs if counts else None,
        "competitor": counts.competitor if counts else "",
        "competitor_count_cs": (
            counts.competitor_count_cs if counts else None
        ),
        "competitor_other_lemma": (
            counts.competitor_other_lemma if counts else None
        ),
        "in_domain": in_domain,
        "source": SOURCE,
        **common,
    }
    for tau in attest.TAUS:
        out[f"attested_{tau}"] = counts.attested[tau] if counts else None
    return out


def _opportunity_summary(
    forms: Sequence[FormFeatures], stem: morphology.StemAnalysis, vowel: str
) -> tuple[list[str], str]:
    chosen = [f for f in forms if f.analysis.opportunity == vowel]
    names = list(dict.fromkeys(f.plural["form"] for f in chosen))
    return names, cell_outcome(stem, chosen) if chosen else ""


def _complete_row(
    row: Mapping[str, str],
    stem: morphology.StemAnalysis,
    forms: Sequence[FormFeatures],
    outcomes: Mapping[str, str],
    doublets: Mapping[str, str],
    common: Mapping[str, Any],
    suffix: Mapping[str, Any],
    details: Details | None,
    ctx: Context,
    forms_by_lexeme: Mapping[str, Sequence[tuple[str, str]]],
    strata: tuple[str, str],
    in_domain: bool,
) -> dict[str, Any]:
    lexeme_id = row["lexeme_id"]
    forms_i, outcome_i = _opportunity_summary(forms, stem, "i")
    forms_e, outcome_e = _opportunity_summary(forms, stem, "e")
    i_forms = [f for f in forms if f.analysis.opportunity == "i"]
    doublet_i = next(
        (
            doublets[f.plural["cell"]]
            for f in i_forms
            if doublets[f.plural["cell"]]
        ),
        "",
    )
    primary = _primary(i_forms, doublet_i)
    verb = first_conjugation_verb(keys.match_key(row["word"]), ctx.verbs)
    adjudication = ctx.adjudications.get(lexeme_id)
    out: dict[str, Any] = {
        key: row[key]
        for key in (
            "lexeme_id",
            "root_id",
            "word",
            "pos",
            "gender",
            "etymology_number",
            "first_gloss",
            "ipa_raw",
            "domain",
            "stress_source",
            "stress_variants",
            "etymon_lang",
            "etymon",
            "plural_status",
            "plural_provenance",
            "countability",
        )
    }
    out.update(common)
    out.update(
        {
            "domain_flags": csvio.split_list(row["domain_flags"]),
            "in_domain": in_domain,
            "dorsal_target": stem.dorsal and not stem.nde_class,
            "ndeb_internal": (
                morphology.root_internal_velar_front(
                    stem, keys.match_key(row["word"])
                )
                if row["domain"] == "in_domain"
                else None
            ),
            "stress_missing": _bool(row["stress_missing"]),
            "ico_stressed": csvio.parse_bool(row["ico_stressed"]),
            "invariable": _bool(row["invariable"]),
            "plurale_tantum": _bool(row["plurale_tantum"]),
            "plural_qc": _rollup(lexeme_id, outcomes, ctx),
            "opportunities": sorted(
                {f.analysis.opportunity for f in forms}
                - {morphology.OPPORTUNITY_NONE}
            ),
            "forms_i": forms_i,
            "outcome_i": outcome_i,
            "forms_e": forms_e,
            "outcome_e": outcome_e,
            "filtered_forms": [f.plural["form"] for f in forms if f.filtered],
            "mutation_i": (
                VARIABLE
                if outcome_i == VARIABLE
                else (primary.analysis.mutation if primary else None)
            ),
            "orth_change_i": primary.analysis.orth_change if primary else "",
            "palatal_output_i": (
                primary.analysis.palatal_output if primary else ""
            ),
            "exception_reason": primary.reason if primary else NON_EXCEPTION,
            "doublet_class": doublet_i,
            "shares_root_with_1st_conj_verb": bool(verb),
            "first_conj_verb": verb,
            "stratum": strata[0],
            "stratum_basis": strata[1],
            "adjudicated_plurals": (
                csvio.split_list(
                    adjudication.values.get("adjudicated_plurals", "")
                )
                if adjudication
                else []
            ),
            "adjudication": (
                adjudication.values.get("decision", "") if adjudication else ""
            ),
            "freq_root_folded": ctx.root_frequency.get(row["root_id"], 0),
            "source": SOURCE,
        }
    )
    out.update(suffix)
    out.update(semantic_class_columns(row, str(suffix["suffixoid"]), ctx))
    out.update(_frequencies(row, ctx.attestation, forms_by_lexeme))
    return out


def root_frequencies(
    attestation_rows: Iterable[Mapping[str, str]],
) -> dict[str, int]:
    """Folded paradigm sum per root: distinct (lowercased) forms of all
    the root's lexemes, each counted once."""
    seen: dict[str, dict[str, int]] = defaultdict(dict)
    for row in attestation_rows:
        seen[row["root_id"]].setdefault(
            keys.match_key(row["form"]), int(row["count_folded"])
        )
    return {root: sum(forms.values()) for root, forms in seen.items()}


# ---------------------------------------------------------------------------
# The stage
# ---------------------------------------------------------------------------


def _forms_by_lexeme(
    attestation_rows: Iterable[Mapping[str, str]],
) -> dict[str, list[tuple[str, str]]]:
    out: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for row in attestation_rows:
        out[row["lexeme_id"]].append((row["cell"], row["form"]))
    return dict(out)


def _body(ctx_run: RunContext, layout: Layout) -> int:
    lexicon_rows = csvio.read_rows(layout.lexicon_full)
    lexicon = {row["lexeme_id"]: row for row in lexicon_rows}
    plurals = group_plurals(csvio.read_rows(layout.plural_forms))
    adjudications, decisions = check_curated(layout, lexicon, plurals)
    attestation_rows = csvio.read_rows(layout.attestation)
    shared = Context(
        details=load_details(layout.en_lexemes),
        verbs=load_verb_keys(layout.en_verbs),
        lemmas=load_lemma_keys(layout.form_index),
        attestation=load_attestation(layout.attestation),
        plural_qc=load_plural_qc(layout.qc),
        adjectival_roots=frozenset(
            row["root_id"] for row in lexicon_rows if row["pos"] == "adj"
        ),
        adjudications=adjudications,
        decisions=decisions,
        root_frequency=root_frequencies(attestation_rows),
    )
    forms_by_lexeme = _forms_by_lexeme(attestation_rows)
    del attestation_rows
    results = [
        lexeme_result(
            row, plurals.get(row["lexeme_id"], []), shared, forms_by_lexeme
        )
        for row in lexicon_rows
    ]
    results.sort(key=lambda r: str(r.row["lexeme_id"]))
    cell_rows = [cell for result in results for cell in result.cells]
    cell_rows.sort(
        key=lambda r: (r["lexeme_id"], CELL_ORDER[r["cell"]], r["item"])
    )
    csvio.write_rows(layout.cells, cell_rows, CELL_FIELDS, ctx=ctx_run)
    csvio.write_rows(
        layout.lexicon_complete,
        (r.row for r in results),
        LEXEME_FIELDS,
        ctx=ctx_run,
    )
    print(
        f"  {len(results):,} lexemes, {len(cell_rows):,} cell rows; curated:"
        f" {len(adjudications)} adjudications, {len(decisions)} "
        "-logo/-fago/-urgo decisions (guards checked)"
    )
    return 0


def outputs(layout: Layout) -> list[Path]:
    return [layout.cells, layout.lexicon_complete]


def run(ctx: RunContext) -> int | None:
    """Run the ``features`` stage (skipped when its manifest matches)."""
    layout = Layout(ctx)
    inputs = [
        layout.lexicon_full,
        layout.plural_forms,
        layout.qc,
        layout.attestation,
        layout.en_lexemes,
        layout.en_verbs,
        layout.form_index,
    ]
    inputs += [
        path
        for path in (layout.plural_adjudications, layout.suffixoid_classes)
        if path.is_file()
    ]
    return run_stage(
        ctx,
        STAGE,
        inputs=inputs,
        params={
            "version": VERSION,
            "filtered_labels": sorted(FILTERED_LABELS),
        },
        outputs=outputs(layout),
        body=lambda: _body(ctx, layout),
    )
