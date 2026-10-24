"""Schema of every Czech CSV and the generated data dictionary (design §8).

One :class:`Table` per CSV the Czech plugin writes: its path, the stage
and module that write it, whether it is tracked, its licence class
(R9), a description, and its columns in file order. The column lists
are the writers' own field tuples (``lexicon.LEXICON_FIELDS``,
``features.CELL_FIELDS``, ``derivations.DERIVATION_FIELDS`` …), so a
writer and its schema cannot drift apart; :func:`describe` documents
every column (explicit entries first, then patterned families such as
``attested_<τ>`` or ``cell_<CELL>``) with a type, a provenance class
and, for controlled vocabularies, the allowed values.

:func:`render` produces ``docs/czech/data_dict.md``; the ``docs`` stage
writes it (:func:`run`), and ``tests/czech/test_schema.py`` checks that
every column is documented, that the files the pipeline wrote have
exactly the schema's header, that their values are in the documented
value sets (:func:`validate_file`: enums, list items with a value set,
booleans), and that the tracked dictionary is up to date. A column can
mean something narrower in one table than in another: a table's
``overrides`` replace the shared description there (the TP tables'
``target_class`` holds row groups such as ``velar:k``, ``dorsal``,
``velar+r``, ``nde:*``; their ``tau`` the strict thresholds). The
``docs`` stage validates every Czech CSV that exists when it runs and
fails on a violation.

Types: ``str``, ``int``, ``float``, ``bool`` (``True``/``False``, R10),
``list`` (``|``-joined, source or sorted order as documented), and
``enum`` (one of the listed values). Provenance: ``wiktionary``
(English or Czech Wiktionary via wiktextract), ``derived`` (computed by
the pipeline from other columns), ``corpus`` (Leipzig
ces_news_2024_1M counts), ``cnk`` (SYN2015 list), ``nc`` (MorfFlex,
DeriNet or czetyl), ``curated`` (hand decisions).
"""

from __future__ import annotations

import hashlib
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from lexpipe import csvio, tp
from lexpipe.attest import TAUS
from lexpipe.context import RunContext, run_stage
from lexpipe.langs.czech import (
    attestation,
    covariation,
    derivations,
    features,
    jsonio,
    lexicon,
    logtable,
    noevidence,
    qc,
    roots,
    spellout,
)
from lexpipe.langs.czech import tpcells as tc
from lexpipe.langs.czech import (
    tpstage,
)
from lexpipe.langs.czech.analysis_layout import analysis_paths

STAGE: Final[str] = "docs"
# 2: lean raw lexicon and full table; R27 columns; gated writer.
# 3: per-table overrides, family value sets, value validation (audit
# 2026-10-23).
VERSION: Final[str] = "3"
# What the dictionary describes (no data is copied into it).
SOURCES: Final[tuple[str, ...]] = ("wiktionary_en", "wiktionary_cs")


@dataclass(frozen=True, slots=True)
class Column:
    """One documented column."""

    name: str
    kind: str
    description: str
    provenance: str = "derived"
    values: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class Table:
    """One CSV file of the Czech plugin."""

    name: str
    path: str
    stage: str
    writer: str
    tracked: bool
    licence: str
    description: str
    fields: tuple[str, ...]
    overrides: Mapping[str, Column] | None = None

    def columns(self) -> list[Column]:
        own = self.overrides or {}
        return [own.get(name) or describe(name) for name in self.fields]


_TAU_TEXT: Final[str] = (
    "attested at τ (≥ 1 token, ≥ 0.5 or ≥ 2 per million word tokens of "
    "Leipzig ces_news_2024_1M; exact rational comparison, R11)"
)
_TARGETS: Final[tuple[str, ...]] = (
    "velar",
    "r",
    "coronal_inferred",
    "c",
    "sibilant",
    "labial",
    "lateral",
    "soft",
    "vowel",
    "other",
    "none",
)
_STRATA: Final[tuple[str, ...]] = (
    "inherited",
    "borrowed whole",
    "internal formation",
    "unknown",
)
_BASES: Final[tuple[str, ...]] = ("linked", "wikt", "wikt+corpus", "all")


def _c(
    name: str,
    kind: str,
    description: str,
    provenance: str = "derived",
    values: Sequence[str] = (),
) -> Column:
    return Column(name, kind, description, provenance, tuple(values))


# Columns documented once, used by every table that has them.
COLUMNS: Final[dict[str, Column]] = {
    c.name: c
    for c in (
        _c("lexeme_id", "str", "lang:word:pos:ety:fp identifier (R4)"),
        _c(
            "root_id",
            "str",
            "root identifier: word + etymology (+ animacy "
            "when one etymology has m.an and m.in nouns, šach) (R4)",
        ),
        _c(
            "lemma", "str", "citation form (NFC, no folding, R8)", "wiktionary"
        ),
        _c(
            "pos",
            "enum",
            "part of speech",
            "wiktionary",
            ("noun", "adj", "name"),
        ),
        _c(
            "etymology",
            "str",
            "Wiktionary etymology number (0 if none)",
            "wiktionary",
        ),
        _c(
            "edition",
            "enum",
            "edition the lexeme comes from",
            "wiktionary",
            ("en", "cs"),
        ),
        _c(
            "paradigm_source",
            "enum",
            "where the primary paradigm comes " "from (none: no table, plech)",
            "wiktionary",
            ("en", "cs", "none"),
        ),
        _c(
            "domain",
            "enum",
            "in or out of the analysis domain (R13, R25)",
            "derived",
            ("in", "out"),
        ),
        _c(
            "ood_class",
            "str",
            "why out of domain (name:*, multiword, "
            "abbreviation, hyphenated, dual_remnant, irregular_stem, tstem, "
            "indeclinable)",
        ),
        _c(
            "gender", "list", "genders of the head code readings", "wiktionary"
        ),
        _c("animacy", "list", "masculine animacies (an, in)", "wiktionary"),
        _c(
            "readings",
            "list",
            "gender/animacy readings (m.an, m.in, f, n, m)",
            "wiktionary",
        ),
        _c(
            "gender_ambiguous",
            "bool",
            "several genders in the head code",
            "wiktionary",
        ),
        _c(
            "animacy_ambiguous",
            "bool",
            "both animacies in the head code",
            "wiktionary",
        ),
        _c(
            "pluralia_tantum",
            "bool",
            "head code marks plural only",
            "wiktionary",
        ),
        _c("by_sense", "bool", "gender by sense (mfbysense)", "wiktionary"),
        _c("head_template", "str", "head template name", "wiktionary"),
        _c(
            "head_code",
            "str",
            "raw gender code of the head template",
            "wiktionary",
        ),
        _c(
            "decl_templates",
            "list",
            "declension templates (source order)",
            "wiktionary",
        ),
        _c(
            "decl_spec",
            "list",
            "raw cs-ndecl specs (audit only, R13)",
            "wiktionary",
        ),
        _c("spec_flags", "list", "flags of the specs (sorted)", "wiktionary"),
        _c(
            "explicit_cells",
            "list",
            "cells overridden in the spec",
            "wiktionary",
        ),
        _c(
            "stem_class",
            "str",
            "wiktextract stem class (irreg-stem, t-stem)",
            "wiktionary",
        ),
        _c("lemma_suffix", "str", "R26 suffix the lemma ends in"),
        _c(
            "suffix_evidence",
            "enum",
            "evidence for the R26 suffix, first that holds "
            "(lexpipe.langs.czech.morphology): etymology template, "
            "suffix category, derivational form-of sense naming a lexeme, "
            "another lexeme's headword line (f=, dim=, adj=), a -ský/-cký "
            "adjective, a base lexeme the word decomposes against (not "
            "for inherited or borrowed words), the -log/-urg/-graf "
            "spelling",
            "derived",
            (
                "etymology",
                "category",
                "formof",
                "headword",
                "adjective",
                "base",
                "orthographic",
                "",
            ),
        ),
        _c(
            "derivational_formof",
            "list",
            "relation:target of derivational form-of senses (diminutive, "
            "augmentative, female_equivalent, male_equivalent, "
            "verbal_noun; source order)",
            "wiktionary",
        ),
        _c("table_animacy", "list", "animacy of each declension table"),
        _c("capitalized", "bool", "lemma starts with a capital"),
        _c("name_evidence", "str", "why the lexeme is a proper name (R25)"),
        _c("first_gloss", "str", "first gloss of the entry", "wiktionary"),
        _c(
            "ipa",
            "list",
            "raw IPA strings (template-generated; a phrase "
            "transcription's prosodic bar is written \\|)",
            "wiktionary",
        ),
        _c(
            "head_args",
            "list",
            "the head templates as written, one item each: key=value pairs "
            "joined by ';' (backslash-escaped), prefixed with 'name:' unless "
            "the template is cs-noun (nouns) or cs-adj (adjectives); "
            "lexpipe.kaikki.head_args",
            "wiktionary",
        ),
        _c(
            "etymology_summary",
            "str",
            "the deciding etymology template, summarized: kind:lang:term "
            "for a source (inh:zlw-ocs:kniha), the kind alone for an "
            "internal formation (af); lexpipe.etymology.summary",
            "wiktionary",
        ),
        _c(
            "ipa_kind",
            "enum",
            "phonemic /…/ or phonetic […] (R8)",
            "derived",
            ("phonemic", "phonetic", ""),
        ),
        _c(
            "head_feminine",
            "list",
            "feminine equivalents on the head line",
            "wiktionary",
        ),
        _c(
            "head_masculine",
            "list",
            "masculine equivalents on the head line",
            "wiktionary",
        ),
        _c(
            "head_diminutive",
            "list",
            "diminutives on the head line",
            "wiktionary",
        ),
        _c(
            "head_adjective",
            "list",
            "relational adjectives on the head line",
            "wiktionary",
        ),
        _c("n_tables", "int", "declension tables of the entry"),
        _c("n_cells", "int", "cells with a primary form"),
        _c("n_cells_second", "int", "cells with a second-source form"),
        _c("second_source_entries", "int", "matched cs.wiktionary entries"),
        _c(
            "source",
            "list",
            "sources of the row's data (R9; tracked tables "
            "only redistributable ones)",
        ),
        _c(
            "licence_class",
            "enum",
            "licence class of the row (R9)",
            "derived",
            ("by-sa", "by", "by-nc-sa", "unverified", "citation"),
        ),
        _c(
            "cell",
            "str",
            "case.number cell (NOM.SG … INS.PL; adjectives add "
            "gender/animacy, NOM.PL.M.AN; ADV = adverb)",
        ),
        _c("form", "str", "word form as written", "wiktionary"),
        _c(
            "primary",
            "bool",
            "form of the lexeme's own edition (False: "
            "second-source cs.wiktionary form)",
        ),
        _c("table", "int", "declension table number (0-based)"),
        _c("order", "int", "position of the form in the entry"),
        _c(
            "form_source",
            "enum",
            "explicit or module-generated form (R11)",
            "derived",
            ("explicit", "module_default"),
        ),
        _c(
            "qualifiers", "list", "tags and raw tags of the form", "wiktionary"
        ),
        _c(
            "syncretic_cells",
            "list",
            "other cells of the table with the same " "form",
        ),
        _c(
            "in_other_edition", "bool", "form also in the other edition's cell"
        ),
        _c(
            "tag_error",
            "str",
            "wiktextract tag errors of the form",
            "wiktionary",
        ),
        _c("count_cs", "int", "case-sensitive tokens of the form", "corpus"),
        _c("count_cap", "int", "tokens of the capitalized form", "corpus"),
        _c(
            "cap_claimed",
            "bool",
            "the capitalized form belongs to another "
            "lexeme, so its tokens are not credited",
        ),
        _c("count_folded", "int", "case-folded tokens", "corpus"),
        _c(
            "count",
            "int",
            "tokens: case-sensitive + sentence-initial "
            "capitals (R24); derivatives: summed over their own forms",
            "corpus",
        ),
        _c("per_million", "float", "count per million word tokens", "corpus"),
        _c(
            "form_ambiguous",
            "bool",
            "the form also belongs to another lexeme " "(R11)",
        ),
        _c("form_ambiguous_folded", "bool", "the same, case-folded"),
        _c(
            "ambiguous_with", "list", "first owners the form is ambiguous with"
        ),
        _c(
            "lemma_count_cs",
            "int",
            "case-sensitive tokens of the lemma",
            "corpus",
        ),
        _c(
            "lemma_count",
            "int",
            "lemma tokens with the capital correction " "(R24 lemma basis)",
            "corpus",
        ),
        _c("lemma_count_folded", "int", "case-folded lemma tokens", "corpus"),
        _c(
            "lemma_ambiguous",
            "bool",
            "the lemma string is also a form of another lexeme (jak 'yak' "
            "vs the conjunction jak): its lemma count is not its own",
        ),
        _c(
            "paradigm_sum_unambiguous",
            "int",
            "tokens summed over distinct primary forms that are not "
            "form_ambiguous (R24 paradigm basis without homograph tokens)",
            "corpus",
        ),
        _c(
            "paradigm_sum",
            "int",
            "tokens summed over distinct primary forms "
            "(R24 paradigm basis)",
            "corpus",
        ),
        _c("paradigm_sum_folded", "int", "the same, case-folded", "corpus"),
        _c("n_forms", "int", "distinct primary forms"),
        _c("n_forms_attested", "int", "of which with a token", "corpus"),
        _c("cnk_lemma_freq", "int", "SYN2015 lemma frequency (≥ 10)", "cnk"),
        _c(
            "gender_class",
            "enum",
            "gender class of the row",
            "derived",
            ("m.an", "m.in", "m", "f", "n", "mixed", "adj"),
        ),
        _c(
            "cell_confounded",
            "bool",
            "the form also fills a cell outside its "
            "copy group (parku: GEN = DAT = LOC.SG)",
        ),
        _c(
            "stem",
            "str",
            "oblique stem (length-folded common prefix of "
            "GEN.SG, DAT.PL, INS.PL, INS.SG)",
        ),
        _c(
            "stem_method",
            "enum",
            "how the stem was found: GEN.SG minus its ending (genitive; "
            "genitive_check when the stem cells' common prefix disagrees: "
            "manual review), the common prefix without GEN.SG (lcp), one "
            "cell form stripped (strip), the lemma, an adjective lemma",
            "derived",
            (
                "genitive",
                "genitive_check",
                "lcp",
                "strip",
                "lemma",
                "adjective",
            ),
        ),
        _c(
            "declension",
            "str",
            "declension type from GEN.SG (hard, soft, "
            "i-stem, adjectival, unknown; adjectives hard/soft/possessive)",
        ),
        _c("stem_final", "str", "stem-final unit (ch is one letter)"),
        _c(
            "cluster",
            "enum",
            "stem-final cluster",
            "derived",
            ("sk", "ck", "st", "zd", "Cr", ""),
        ),
        _c(
            "target_class",
            "enum",
            "class of the stem final (§6)",
            "derived",
            _TARGETS,
        ),
        _c(
            "coronal_inferred",
            "bool",
            "t/d/n target: change inferred from " "ě/i/í (§6)",
        ),
        _c("core_target", "bool", "velar or r target"),
        _c("target_is_suffix", "bool", "target inside an R26 suffix"),
        _c("aligned", "bool", "form = (mutated) stem + ending"),
        _c("ending", "str", "ending after the stem"),
        _c("orth_change", "str", "orthographic mutation (k→c, h→z …)"),
        _c("mutation", "bool", "the form shows the stem final's reflex"),
        _c(
            "mutation_kind",
            "enum",
            "how the change is seen",
            "derived",
            ("orthographic", "inferred", "labial_glide", ""),
        ),
        _c("palatal_output", "str", "IPA of the palatalized final"),
        _c(
            "allomorph_class",
            "enum",
            "§6 class of the ending",
            "derived",
            ("P", "N", "none"),
        ),
        _c(
            "n_kind",
            "enum",
            "kind of N ending",
            "derived",
            ("allomorph", "nontrigger", "soft", ""),
        ),
        _c("front_vowel_ending", "bool", "ending starts with e é ě i í"),
        _c("opportunity", "bool", "the ending is P or N in its cell"),
        _c(
            "exception_label",
            "str",
            "undergoer, unexplained, non exception "
            "(allomorph), nde:<class>, non exception",
        ),
        _c(
            "curated_status",
            "enum",
            "hand decision on the form",
            "curated",
            ("confirmed", "rejected", ""),
        ),
        _c(
            "nde_class",
            "enum",
            "non-derived-environment class (R16)",
            "derived",
            ("tautomorphemic_dorsal", "loan_hard_coronal", ""),
        ),
        _c("ndeb_tokens", "list", "root-internal C + front vowel sequences"),
        _c("n_opportunity_rows", "int", "primary rows with a P/N choice"),
        _c("n_unexplained", "int", "P rows without the change"),
        _c("n_unaligned", "int", "rows that do not align with the stem"),
        # QC (MorfFlex, NC)
        _c("origin", "str", "wiktionary or morfflex_only row", "nc"),
        _c("status", "str", "QC status against MorfFlex", "nc"),
        _c("cell_agreement", "str", "cell agreement with MorfFlex", "nc"),
        _c("morfflex_var", "str", "MorfFlex variant digit", "nc"),
        _c("morfflex_lemmas", "list", "joined MorfFlex lemmas", "nc"),
        _c("name_labels", "list", "MorfFlex name labels", "nc"),
        _c(
            "wiktionary_forms",
            "list",
            "Wiktionary forms of the cell",
            "wiktionary",
        ),
        _c("ijp_url", "str", "percent-encoded IJP lookup URL (never fetched)"),
        # derivations
        _c("base_lexeme_id", "str", "lexeme_id of the base noun"),
        _c("base_root_id", "str", "root_id of the base noun"),
        _c("base_lemma", "str", "base lemma", "wiktionary"),
        _c(
            "base_readings",
            "list",
            "base gender/animacy readings",
            "wiktionary",
        ),
        _c("base_stem", "str", "oblique stem of the base"),
        _c(
            "base_target_class",
            "enum",
            "target class of the base",
            "derived",
            _TARGETS,
        ),
        _c("base_stem_final", "str", "stem-final unit of the base"),
        _c("base_cluster", "str", "stem-final cluster of the base"),
        _c("base_suffix_target", "bool", "base target inside an R26 suffix"),
        _c("base_domain", "str", "in, or out:dual_remnant (ruka, oko)"),
        _c(
            "base_stratum_wikt",
            "enum",
            "R22 stratum of the base from " "Wiktionary only",
            "wiktionary",
            _STRATA,
        ),
        _c(
            "base_stratum",
            "enum",
            "R22 stratum: Wiktionary > czetyl > " "DeriNet Loanword",
            "nc",
            _STRATA,
        ),
        _c(
            "base_stratum_source",
            "enum",
            "which source decided",
            "derived",
            ("wiktionary", "czetyl", "derinet", ""),
        ),
        _c("derivative", "str", "derivative lemma"),
        _c(
            "derivative_pos",
            "enum",
            "part of speech of the suffix",
            "derived",
            ("noun", "adj", "verb"),
        ),
        _c("suffix", "str", "suffix label (suffixes.SUFFIXES)"),
        _c(
            "suffix_class",
            "enum",
            "P palatalizing, N non-palatalizing, NT "
            "non-trigger (-ista, -ismus)",
            "derived",
            ("P", "N", "NT"),
        ),
        _c(
            "suffix_pair",
            "enum",
            "competing pair (R20, R23)",
            "derived",
            ("ni_ovy", "it_ovat", "ka_yne", ""),
        ),
        _c(
            "suffix_palatalizes_base",
            "bool",
            "the suffix palatalizes the " "base's target class",
        ),
        _c("prefix", "str", "prefix of a linked derivative (s-nož-it)"),
        _c(
            "outcome",
            "enum",
            "pal (changed), plain (unchanged, not "
            "expected to change), exception (P suffix, unchanged: barokní), "
            "nontarget",
            "derived",
            ("pal", "plain", "exception", "nontarget"),
        ),
        _c(
            "mutation_size",
            "enum",
            "first (big) or second (small) " "palatalization",
            "derived",
            ("big", "small", ""),
        ),
        _c(
            "evidence",
            "list",
            "wikt_ety, wikt_head, wikt_derived, "
            "wikt_related, derinet, cand_headword, cand_headword_plain "
            "(a headword candidate in the plain spelling of a palatalizing "
            "suffix), cand_corpus (sorted)",
        ),
        _c(
            "base_deverbal_from",
            "list",
            "verbs the base's own etymology derives it from (deverbal, "
            "back-formation, 'Formed from', 'From the verb'; source order)",
            "wiktionary",
        ),
        _c(
            "derivative_feminine_agent",
            "bool",
            "a feminizing derivative: -yně, or -ka linked as the base's "
            "f= or glossed as a female agent (R20, R23 ka_yne)",
        ),
        _c(
            "basis",
            "enum",
            "narrowest evidence basis (linked ⊂ wikt ⊂ " "wikt+corpus ⊂ all)",
            "derived",
            _BASES,
        ),
        _c(
            "headword_editions",
            "list",
            "Wiktionary editions with the "
            "derivative as a headword of its part of speech",
            "wiktionary",
        ),
        _c("derinet_child_id", "str", "DeriNet node id", "nc"),
        _c(
            "derinet_verbalizer",
            "str",
            "verbalizer from DeriNet's last D " "morph (R14)",
            "nc",
        ),
        _c("derinet_conj_class", "str", "DeriNet ConjugClass", "nc"),
        _c(
            "verbalizer_consistent",
            "bool",
            "lemma suffix = DeriNet " "verbalizer",
            "nc",
        ),
        _c(
            "conj_class_consistent",
            "bool",
            "ConjugClass agrees with the " "verbalizer",
            "nc",
        ),
        _c(
            "derivative_stratum_wikt",
            "enum",
            "R22 stratum of the derivative " "from Wiktionary only",
            "wiktionary",
            _STRATA,
        ),
        _c(
            "derivative_stratum",
            "enum",
            "R22 stratum: Wiktionary > czetyl",
            "nc",
            _STRATA,
        ),
        _c(
            "derivative_stratum_source",
            "enum",
            "which source decided",
            "derived",
            ("wiktionary", "czetyl", ""),
        ),
        _c(
            "n_forms_counted",
            "int",
            "derivative forms counted (forms of " "other entries left out)",
        ),
        _c(
            "stratum_wikt",
            "enum",
            "R22 stratum from Wiktionary only",
            "wiktionary",
            _STRATA,
        ),
        _c("stratum_wikt_basis", "str", "template kind that decided"),
        _c("stratum", "enum", "R22 stratum from all sources", "nc", _STRATA),
        _c(
            "stratum_source",
            "enum",
            "which source decided",
            "derived",
            ("wiktionary", "czetyl", "derinet", ""),
        ),
        _c("tau", "enum", "attestation threshold (R11)", "derived", TAUS),
        _c("n_pairs", "int", "base × derivative pairs"),
        _c("n_bases", "int", "distinct base roots"),
        # TP long
        _c("lang", "str", "language code"),
        _c(
            "rule",
            "enum",
            "palatalize (P-default) or retain (N-default)",
            "derived",
            tp.RULES,
        ),
        _c("trigger", "str", "trigger cell, competition cell or R17 class"),
        _c(
            "subset",
            "str",
            "root subset (all, core, suffix, stratum:<s>, "
            "class:<c>, nested:<c>, underspecified, prespecified)",
        ),
        _c(
            "evidence_set",
            "enum",
            "R15 evidence set",
            "derived",
            ("E0", "E1", "E2"),
        ),
        _c("coding", "enum", "doublet coding (R19)", "derived", tp.CODINGS),
        _c(
            "freq_basis",
            "enum",
            "R24 frequency basis: lemma-form tokens, the paradigm sum over "
            "attested forms, or that sum without forms of other lexemes",
            "derived",
            ("lemma", "paradigm", "paradigm_unambiguous"),
        ),
        _c("V", "int", "vocabulary size of the R24 frame (empty: all)"),
        _c("N_target", "int", "target-final roots of the row (R16)"),
        _c("N_trigger", "int", "roots attested before the trigger"),
        _c("e_nonalt", "int", "N_trigger roots that never alternate"),
        _c("n_shielded", "int", "roots attested only with N allomorphs"),
        _c("e", "int", "exceptions of the row's rule"),
        _c("theta_trigger", "float", "N_trigger / ln N_trigger"),
        _c("theta_target", "float", "N_target / ln N_target"),
        _c("tolerated_trigger", "bool", "e ≤ θ(N_trigger)"),
        _c("tolerated_counterfactual", "bool", "the counterfactual of R16"),
        _c("majority", "bool", "rule-following roots outnumber e"),
        _c("memo", "str", "(a)/(b)/(c)/R17 label of the row"),
        _c(
            "fragility",
            "int",
            "fewest item changes that flip tolerated " "(learners.fragility)",
        ),
        # learners
        _c("table", "str", "analysis table"),
        _c("n_samples", "int", "simulated vocabularies per V"),
        _c("seed", "int", "seed of the run context (R5)"),
        _c("nested", "str", "nested rules of the verdict"),
        _c(
            "share_palatalize",
            "float",
            "share of learners with palatalize " "productive",
        ),
        _c("share_retain", "float", "share with retain productive"),
        _c("share_neither", "float", "share with neither"),
        _c("share_conditioned", "float", "share with a nested rule set"),
        # per-root TP file and stage files
        _c("root", "str", "root record key (root_id; #adj for adjectives)"),
        _c("genders", "list", "gender classes of the root's rows"),
        _c("suffix_target", "bool", "target inside an R26 suffix"),
        _c("monosyllable", "bool", "one vowel nucleus"),
        _c(
            "lemma_freq",
            "float",
            "lemma tokens (max over the root's " "lexemes)",
            "corpus",
        ),
        _c(
            "paradigm_freq",
            "float",
            "paradigm tokens (max over lexemes)",
            "corpus",
        ),
        _c(
            "root_set",
            "enum",
            "spell-out root set (R20)",
            "derived",
            ("feminine", "E0-certain"),
        ),
        _c(
            "pair",
            "enum",
            "competing pair",
            "derived",
            ("ni_ovy", "it_ovat", "ka_yne", "any"),
        ),
        _c(
            "p_derivatives", "list", "derivatives with the palatalizing member"
        ),
        _c("n_derivatives", "list", "derivatives with the other member"),
        _c("p_count", "int", "their summed tokens", "corpus"),
        _c("n_count", "int", "their summed tokens", "corpus"),
        _c(
            "suffixoid",
            "enum",
            "suffixoid ending (R27)",
            "derived",
            ("-log", "-urg", "-graf", "-fil"),
        ),
        _c(
            "learned_suffixoid",
            "bool",
            "the learned suffixoid, not a lookalike (shared detector, "
            "lexpipe.suffixoids; lookalikes are left out of the aggregates)",
        ),
        _c(
            "suffixoid_basis",
            "str",
            "what decided learned_suffixoid (etymology, category, combining "
            "-o- (cs gloss), gloss, gloss (profession…), manual)",
        ),
        _c(
            "semantic_class",
            "enum",
            "person/thing classifier (R27; shared, lexpipe.suffixoids)",
            "derived",
            ("person", "animate", "thing", "unknown"),
        ),
        _c("class_basis", "str", "which classifier rule decided"),
        _c("gloss", "str", "first gloss", "wiktionary"),
        _c("nom_pl", "list", "NOM.PL forms", "wiktionary"),
        _c("nom_pl_allomorph", "list", "NOM.PL endings (-ové, -i, -y …)"),
        _c("loc_sg", "list", "LOC.SG forms", "wiktionary"),
        _c("loc_sg_allomorph", "list", "LOC.SG endings"),
        _c("loc_pl", "list", "LOC.PL forms", "wiktionary"),
        _c("loc_pl_allomorph", "list", "LOC.PL endings"),
        _c("loc_pl_palatalized", "bool", "a P LOC.PL form shows the change"),
        _c(
            "derivatives",
            "list",
            "derivative:suffix:outcome (basis linked " "or wikt)",
        ),
    )
}

# Families of patterned column names: (regex, kind, description,
# provenance).
FAMILIES: Final[tuple[tuple[str, str, str, str], ...]] = (
    (r"cell_(.+)", "list", "primary forms of cell {0}", "wiktionary"),
    (r"attested_(ge.+)", "bool", "form " + _TAU_TEXT + "; τ = {0}", "corpus"),
    (
        r"summary_E0_(.+?)(_strict)?",
        "enum",
        "E0 lexeme summary at τ = {0} "
        "(undergoer, resister, mixed, no_evidence, n/a, out){1}",
        "derived",
    ),
    (
        r"p_cells_E0_(.+)",
        "list",
        "cells with alternating P forms at τ = {0}",
        "derived",
    ),
    (
        r"n_cells_E0_(.+)",
        "list",
        "cells with competing N forms at τ = {0}",
        "derived",
    ),
    (
        r"slots_(.+)_(dictionary|ge1)",
        "str",
        "slots attested in {0} at τ = "
        "{1}: A (palatalizing / front), B (competitor), AB, empty",
        "derived",
    ),
    (
        r"evidence_(E\d)_(dictionary|ge1)",
        "bool",
        "alternation evidence in " "{0} at τ = {1} (all cells)",
        "derived",
    ),
    (
        r"certain_(.+)",
        "bool",
        "certain specification at τ = {0} (R20)",
        "derived",
    ),
    (
        r"choice_(.+)",
        "enum",
        "pair choice at τ = {0}: P, N, both, empty",
        "derived",
    ),
    (
        r"base_class_(.+)",
        "enum",
        "R21 base class at τ = {0}: alternating, "
        "non-alternating, no evidence",
        "derived",
    ),
    (
        r"(LOC\.SG|NOM\.PL)_(any|dominant)_(.+)",
        "enum",
        "{0} class (P, N, " "empty) with doublet coding {1} at τ = {2} (R23)",
        "derived",
    ),
    (
        r"(nom_pl|loc_sg|loc_pl)_count",
        "int",
        "tokens of the {0} forms",
        "corpus",
    ),
    (
        r"(nom_pl|loc_sg|loc_pl)_attested_(.+)",
        "bool",
        "a {0} form attested " "at τ = {1}",
        "corpus",
    ),
)


# Value sets of the enum families (the description lists them too).
FAMILY_VALUES: Final[Mapping[str, tuple[str, ...]]] = {
    r"summary_E0_(.+?)(_strict)?": (
        "undergoer",
        "resister",
        "mixed",
        "no_evidence",
        "n/a",
        "out",
    ),
    r"choice_(.+)": ("P", "N", "both", ""),
    r"base_class_(.+)": ("alternating", "non-alternating", "no evidence"),
    r"(LOC\.SG|NOM\.PL)_(any|dominant)_(.+)": ("P", "N", ""),
}


def describe(name: str) -> Column:
    """The documented column *name* (KeyError if undocumented)."""
    if name in COLUMNS:
        return COLUMNS[name]
    for pattern, kind, text, provenance in FAMILIES:
        match = re.fullmatch(pattern, name)
        if match is not None:
            groups = [g or "" for g in match.groups()]
            if pattern.startswith("summary"):
                groups[1] = (
                    " (strict: no ambiguous or confounded forms)"
                    if groups[1]
                    else ""
                )
            return Column(
                name,
                kind,
                text.format(*groups),
                provenance,
                FAMILY_VALUES.get(pattern, ()),
            )
    raise KeyError(f"undocumented Czech column: {name!r}")


def _tp_targets() -> tuple[str, ...]:
    """Every ``target_class`` of the TP and learner rows: the target
    classes, the velar segments and r clusters, the competition
    groups, the R17 targets, the spell-out groups and the NDE rows."""
    found: list[str] = []
    for name in (
        *tc.RULE_TARGETS,
        *tc.COMPETITION_TARGETS,
        *tc.R17_TARGETS,
        *tpstage.SPELLOUT_TARGETS,
        *tpstage.DERIVATION_TARGETS,
        *(f"nde:{name}" for name in tc.NDE_CLASSES),
    ):
        if name not in found:
            found.append(name)
    return tuple(found)


_TP_OVERRIDES: Final[dict[str, Column]] = {
    "target_class": _c(
        "target_class",
        "enum",
        "row group: a target class (§6), a velar segment or cluster "
        "(velar:k … velar:ck), r after a consonant or a vowel (r:Cr, "
        "r:Vr), a competition group (baseline and its parts), an R17 "
        "target (dorsal), the spell-out pool (velar+r) or an NDE class",
        "derived",
        _tp_targets(),
    ),
    "tau": _c(
        "tau",
        "enum",
        "attestation threshold (R11) and its strict variant (no "
        "ambiguous or cell-confounded forms)",
        "derived",
        roots.ALL_TAUS,
    ),
}
_BASIS_LIST: Final[Column] = _c(
    "basis",
    "list",
    "bases of the root's pair derivatives (sorted)",
    "derived",
    _BASES,
)
_WORKLIST_FORM_SOURCE: Final[Column] = _c(
    "form_source",
    "list",
    "form_source of the cell's forms (R11; sorted)",
    "derived",
    ("explicit", "module_default"),
)
_QC_FORM_SOURCE: Final[Column] = _c(
    "form_source",
    "enum",
    "explicit or module-generated Wiktionary form (R11); empty for a "
    "MorfFlex-only row",
    "derived",
    ("explicit", "module_default", ""),
)


def validate_rows(
    table: Table, rows: Iterable[Mapping[str, str]], limit: int = 5
) -> list[str]:
    """Values outside a column's documented set: enums, items of lists
    with a value set, and booleans (``True``/``False`` or empty, R10).
    At most *limit* examples per column, then a count."""
    columns = [c for c in table.columns() if c.values or c.kind == "bool"]
    bad: dict[str, list[str]] = {c.name: [] for c in columns}
    counts: dict[str, int] = {c.name: 0 for c in columns}
    for row in rows:
        for column in columns:
            if column.name not in row:
                continue  # a file's header is checked separately
            value = row[column.name]
            if column.kind == "bool":
                wrong = [] if value in ("True", "False", "") else [value]
            elif column.kind == "list":
                wrong = [
                    item
                    for item in csvio.split_list(value)
                    if item not in column.values
                ]
            else:
                wrong = [] if value in column.values else [value]
            for item in wrong:
                counts[column.name] += 1
                if (
                    len(bad[column.name]) < limit
                    and item not in bad[column.name]
                ):
                    bad[column.name].append(item)
    return [
        f"{table.name}.{name}: {counts[name]:,} value(s) outside the "
        f"documented set, e.g. {', '.join(repr(v) for v in bad[name])}"
        for name in bad
        if counts[name]
    ]


def validate_file(table: Table, root: Path) -> list[str]:
    """:func:`validate_rows` over the table's file under *root* ([] when
    the file is absent)."""
    path = root / table.path
    if not path.is_file():
        return []
    return validate_rows(table, jsonio.iter_csv(path))


def tables() -> tuple[Table, ...]:
    """Every Czech CSV (paths relative to the repository root)."""
    data, nc, analysis = "data/czech", "data/czech/nc", "analysis/czech"
    return (
        Table(
            "lexicon_raw",
            f"{data}/czech_lexicon_raw.csv",
            "lexicon",
            "lexicon.py",
            True,
            "by-sa",
            "the tracked raw lexicon: one row per lexeme with the columns "
            "extracted from its primary entry only (nothing derived); "
            "proper names and multiword expressions are left out (they are "
            "in lexicon_full)",
            lexicon.RAW_FIELDS,
        ),
        Table(
            "lexicon_full",
            f"{data}/czech_lexicon_full.csv",
            "lexicon",
            "lexicon.py",
            False,
            "by-sa",
            "every lexeme with every column of the lexicon stage "
            "(Wiktionary only; read by the later stages)",
            lexicon.LEXICON_FIELDS,
        ),
        Table(
            "cells_raw",
            f"{data}/czech_cells_raw.csv",
            "lexicon",
            "lexicon.py",
            False,
            "by-sa",
            "one row per lexeme × cell × form × source",
            lexicon.CELL_FIELDS,
        ),
        Table(
            "morfflex_qc",
            f"{nc}/czech_morfflex_qc.csv",
            "qc",
            "qc.py",
            False,
            "by-nc-sa",
            "per-form agreement with MorfFlex",
            qc.QC_FIELDS,
            {"form_source": _QC_FORM_SOURCE},
        ),
        Table(
            "manual_worklist",
            f"{analysis}/czech_manual_worklist.csv",
            "qc",
            "qc.py",
            True,
            "by-sa",
            "forms to check by hand in IJP",
            qc.WORKLIST_FIELDS,
            {"form_source": _WORKLIST_FORM_SOURCE},
        ),
        Table(
            "attestation",
            f"{data}/czech_attestation.csv",
            "attest",
            "attestation.py",
            False,
            "unverified",
            "corpus counts of every in-domain cell form",
            attestation.ATTESTATION_FIELDS,
        ),
        Table(
            "lexeme_frequency",
            f"{data}/czech_lexeme_frequency.csv",
            "attest",
            "attestation.py",
            False,
            "unverified",
            "lemma and paradigm frequencies (R24)",
            attestation.FREQUENCY_FIELDS,
        ),
        Table(
            "cells",
            f"{data}/czech_cells.csv",
            "features",
            "features.py",
            False,
            "unverified",
            "one row per in-domain cell form with " "its §6 classification",
            features.CELL_FIELDS,
        ),
        Table(
            "lexicon_complete",
            f"{data}/czech_lexicon_complete.csv",
            "features",
            "features.py",
            False,
            "unverified",
            "every lexeme with features, frequencies and summaries",
            features.COMPLETE_FIELDS,
        ),
        Table(
            "derivations",
            f"{nc}/czech_derivations.csv",
            "derivations",
            "derivations.py",
            False,
            "by-nc-sa",
            "base × derivative pairs (R14, R20–R22)",
            derivations.DERIVATION_FIELDS,
        ),
        Table(
            "lexeme_strata",
            f"{nc}/czech_lexeme_strata.csv",
            "derivations",
            "derivations.py",
            False,
            "by-nc-sa",
            "R22 strata of every noun and adjective lexeme",
            derivations.STRATA_FIELDS,
        ),
        Table(
            "derivations_summary",
            f"{analysis}/czech_derivations_summary.csv",
            "derivations",
            "derivations.py",
            True,
            "by-sa",
            "pair counts by basis, τ, target, suffix, outcome, stratum",
            derivations.SUMMARY_FIELDS,
        ),
        Table(
            "tp_long",
            f"{analysis}/czech_tp_long.csv",
            "tp",
            "tpstage.py",
            True,
            "by-sa",
            "every TP row of (a), (b), (c), R17 and NDE",
            tpstage.FIELDS,
            _TP_OVERRIDES,
        ),
        Table(
            "learners",
            f"{analysis}/czech_learners.csv",
            "tp",
            "tpstage.py",
            True,
            "by-sa",
            "simulated learners (R18)",
            tpstage.LEARNER_FIELDS,
            _TP_OVERRIDES,
        ),
        Table(
            "tp_roots",
            f"{data}/czech_tp_roots.csv",
            "tp",
            "tpstage.py",
            False,
            "unverified",
            "per-root slots and evidence",
            tpstage.ROOT_FIELDS,
        ),
        Table(
            "spellout_roots",
            f"{nc}/czech_spellout_roots.csv",
            "spellout",
            "spellout.py",
            False,
            "by-nc-sa",
            "spell-out roots × pairs",
            spellout.ROOT_FIELDS,
            {"basis": _BASIS_LIST},
        ),
        Table(
            "noevidence_pairs",
            f"{nc}/czech_noevidence_pairs.csv",
            "noevidence",
            "noevidence.py",
            False,
            "by-nc-sa",
            "P-suffix pairs of velar and r roots with base classes",
            noevidence.PAIR_FIELDS,
        ),
        Table(
            "covariation_units",
            f"{nc}/czech_covariation_units.csv",
            "covariation",
            "covariation.py",
            False,
            "by-nc-sa",
            "root × pair units with classes and choices",
            covariation.UNIT_FIELDS,
            {"basis": _BASIS_LIST},
        ),
        Table(
            "logtable",
            f"{analysis}/czech_logtable.csv",
            "logtable",
            "logtable.py",
            True,
            "by-sa",
            "-log/-urg/-graf/-fil nouns (lexicon facts only)",
            logtable.FIELDS,
        ),
        Table(
            "logtable_full",
            f"{data}/czech_logtable_full.csv",
            "logtable",
            "logtable.py",
            False,
            "unverified",
            "the same with corpus counts",
            logtable.FIELDS + logtable.FULL_EXTRA,
        ),
    )


def render() -> str:
    """``docs/czech/data_dict.md``."""
    lines = [
        "# Czech data dictionary",
        "",
        "Generated by `lexpipe.langs.czech.schema` (the `docs` stage); edit "
        "the schema, not this file. Types: str, int, float, bool "
        "(`True`/`False`, R10), list (`|`-joined), enum. Provenance: "
        "wiktionary, derived, corpus (Leipzig), cnk, nc (MorfFlex, DeriNet, "
        "czetyl), curated.",
        "",
        "| Table | Path | Stage | Tracked | Licence class |",
        "|---|---|---|---|---|",
    ]
    for table in tables():
        lines.append(
            f"| {table.name} | `{table.path}` | {table.stage} | "
            f"{'yes' if table.tracked else 'no'} | {table.licence} |"
        )
    for table in tables():
        lines += [
            "",
            f"## {table.name}",
            "",
            f"`{table.path}`: {table.description}. Written by "
            f"`{table.writer}` (stage `{table.stage}`); "
            f"{'tracked' if table.tracked else 'not tracked'}; licence class "
            f"{table.licence}.",
            "",
            "| Column | Type | Provenance | Description |",
            "|---|---|---|---|",
        ]
        for column in table.columns():
            text = column.description
            if column.values:
                shown = ", ".join(
                    f"`{v}`" if v else "(empty)" for v in column.values
                )
                text += f" ({shown})"
            text = text.replace("|", "\\|")
            lines.append(
                f"| `{column.name}` | {column.kind} | {column.provenance} | "
                f"{text} |"
            )
    return "\n".join(lines) + "\n"


def run(ctx: RunContext) -> int | None:
    """Stage ``docs``: write ``docs/czech/data_dict.md``.

    The dictionary depends on code only, so the manifest records the
    SHA-256 of the rendered text: a schema change re-runs the stage.
    """
    apaths = analysis_paths(ctx)
    text = render()
    digest = hashlib.sha256(text.encode("utf-8")).hexdigest()

    def body() -> int:
        problems = [
            problem
            for table in tables()
            for problem in validate_file(table, ctx.root)
        ]
        for problem in problems:
            print(f"  schema: {problem}", flush=True)
        if problems:
            return 1
        ctx.write_tracked_text(apaths.data_dict, text, sources=SOURCES)
        return 0

    return run_stage(
        ctx,
        STAGE,
        inputs=[],
        params={"version": VERSION, "content_sha256": digest},
        outputs=[apaths.data_dict],
        body=body,
    )
