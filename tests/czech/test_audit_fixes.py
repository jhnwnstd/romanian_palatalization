"""Unit tests of the fixes made after the adversarial audit (2026-10-23).

Each block names the audit finding it covers: R26 suffix targets
identified morphologically, the -ka/-yně and -it/-ovat pair members of
R20/R23, homograph lemma frequencies, strict thresholds, derivational
form-of lexemes, homograph link attachment, the direction of
derivation, adapted loans vs hybrids, stems, names, -a masculines and
unknown animacy, -icismus and -ista, the R27 classifier, plain headword
candidates and the root-internal NDEB tokens.
"""

from __future__ import annotations

from typing import Any

import pytest

from lexpipe import tp
from lexpipe.langs.czech import attestation, derivations
from lexpipe.langs.czech import etymology as ety
from lexpipe.langs.czech import extract, logtable
from lexpipe.langs.czech import morphology as morph
from lexpipe.langs.czech import names, noevidence, ortho, roots, spellout
from lexpipe.langs.czech import suffixes as sx
from lexpipe.langs.czech import tpcells as tc
from lexpipe.langs.czech import tpstage

pytestmark = pytest.mark.cs

ALL: tuple[bool, ...] = (True, True, True)


# ---------------------------------------------------------------------------
# Derivational form-of entries are lexemes (extract)
# ---------------------------------------------------------------------------


def _sense(gloss: str, target: str, tags: tuple[str, ...]) -> dict[str, Any]:
    return {
        "glosses": [gloss],
        "tags": list(tags),
        "form_of": [{"word": target}],
    }


def test_derivational_form_of_senses() -> None:
    diminutive = _sense(
        "diminutive of kluk", "kluk", ("diminutive", "form-of")
    )
    female = _sense(
        "female equivalent of Čech", "Čech", ("feminine", "form-of")
    )
    verbal = _sense("verbal noun of kvasit", "kvasit", ("form-of",))
    case = _sense("genitive plural of kniha", "kniha", ("form-of", "genitive"))
    comparative = _sense(
        "comparative degree of lehký", "lehký", ("comparative", "form-of")
    )
    assert extract.derivational_relation(diminutive) == "diminutive"
    assert extract.derivational_relation(female) == "female_equivalent"
    assert extract.derivational_relation(verbal) == "verbal_noun"
    assert extract.derivational_relation(case) == ""
    assert extract.derivational_relation(comparative) == ""
    assert extract.is_lexeme_entry({"senses": [diminutive]})
    assert not extract.is_lexeme_entry({"senses": [case]})
    assert not extract.is_lexeme_entry({"senses": [comparative]})


def test_derivational_form_of_gives_no_form_of_index_row() -> None:
    entry = {
        "word": "klouček",
        "pos": "noun",
        "senses": [
            _sense("diminutive of kluk", "kluk", ("diminutive", "form-of"))
        ],
        "forms": [],
    }
    kinds = [row[2] for row in extract.index_rows(entry, "cs")]
    assert kinds == ["lemma"]


# ---------------------------------------------------------------------------
# R26 suffix targets identified morphologically
# ---------------------------------------------------------------------------


def _context() -> morph.SuffixContext:
    index = morph.BaseIndex()
    for lemma, pos in (
        ("sedlo", "noun"),
        ("ruka", "noun"),
        ("kluk", "noun"),
        ("ban", "noun"),
        ("mat", "noun"),
        ("Čech", "noun"),
    ):
        index.add(lemma, pos)
    return morph.SuffixContext(
        lemmas=frozenset({"kluk", "Čech", "kniha"}),
        head_named=frozenset({"knížka"}),
        index=index,
    )


@pytest.mark.parametrize(
    ("lemma", "pos", "suffix", "formof", "stratum", "evidence"),
    [
        ("klouček", "noun", "ek", ("kluk",), "unknown", "formof"),
        ("knížka", "noun", "ka", (), "unknown", "headword"),
        ("typický", "adj", "cký", (), "borrowed whole", "adjective"),
        ("český", "adj", "ský", (), "inherited", "adjective"),
        ("sedlák", "noun", "ák", (), "unknown", "base"),
        ("ručník", "noun", "ník", (), "internal formation", "base"),
        ("banka", "noun", "ka", (), "borrowed whole", ""),
        ("matka", "noun", "ka", (), "inherited", ""),
        ("ruka", "noun", "", (), "inherited", ""),
    ],
)
def test_morphological_evidence(
    lemma: str,
    pos: str,
    suffix: str,
    formof: tuple[str, ...],
    stratum: str,
    evidence: str,
) -> None:
    found = morph.morphological_evidence(
        lemma,
        pos,
        suffix,
        formof_targets=formof,
        stratum=stratum,
        context=_context(),
    )
    assert found == evidence


def test_base_index_undoes_the_first_palatalization() -> None:
    context = _context()
    assert context.index.bases("ručník", "ník") == ["ruka"]
    assert context.index.bases("klouček", "ek") == ["kluk"]
    assert context.index.bases("Češka", "ka") == ["Čech"]


def test_adjective_suffix_targets_count_with_any_evidence() -> None:
    target = ortho.target_info("typick")
    assert ortho.target_is_suffix("cký", "adjective", target)
    assert ortho.lemma_suffix("typický") == "cký"


# ---------------------------------------------------------------------------
# Stems: GEN.SG first, palatal spelling, plural-only nouns, adjectival -á
# ---------------------------------------------------------------------------


def _cells(**kwargs: list[str]) -> dict[str, list[str]]:
    return {key.replace("_", "."): value for key, value in kwargs.items()}


@pytest.mark.parametrize(
    ("lemma", "cells", "stem", "method", "target"),
    [
        (
            "oheň",
            _cells(
                GEN_SG=["ohně"],
                DAT_PL=["ohňům"],
                INS_PL=["ohni"],
                INS_SG=["ohněm"],
            ),
            "ohň",
            "genitive",
            "soft",
        ),
        (
            "kůň",
            _cells(
                GEN_SG=["koně"],
                DAT_PL=["koním"],
                INS_PL=["koni"],
                INS_SG=["koněm"],
            ),
            "koň",
            "genitive",
            "soft",
        ),
        (
            "Němec",
            _cells(
                GEN_SG=["Němce"],
                DAT_PL=["Němcům"],
                INS_PL=["Němci"],
                INS_SG=["Němcem"],
            ),
            "Němc",
            "genitive",
            "c",
        ),
        (
            "míra",
            _cells(
                GEN_SG=["míry"],
                DAT_PL=["měrám", "mírám"],
                INS_PL=["měrami"],
                INS_SG=["měrou", "mírou"],
            ),
            "mír",
            "genitive_check",
            "r",
        ),
        (
            "kabátník",
            _cells(
                GEN_SG=["kabátníka"],
                DAT_PL=["kabátníkům"],
                INS_PL=["kabátníci"],
                INS_SG=["kabátníkem"],
            ),
            "kabátník",
            "genitive_check",
            "velar",
        ),
        ("fleky", _cells(DAT_PL=["flekům"]), "flek", "strip", "velar"),
        (
            "pestřec",
            _cells(
                GEN_SG=["pestřce / pestřece"],
                DAT_PL=["pestřcům"],
                INS_PL=["pestřci"],
                INS_SG=["pestřcem"],
            ),
            "pestřc",
            "genitive",
            "c",
        ),
    ],
)
def test_genitive_stems(
    lemma: str,
    cells: dict[str, list[str]],
    stem: str,
    method: str,
    target: str,
) -> None:
    found = ortho.analyse_stem("noun", lemma, cells)
    assert (found.stem, found.method, found.target.target_class) == (
        stem,
        method,
        target,
    )


def test_palatal_stems_align_with_their_spellings() -> None:
    assert ortho.align("ohň", "ohně").ending == "ě"
    assert ortho.align("ohň", "ohni").ending == "i"
    assert ortho.align("ohň", "ohňům").ending == "ům"


def test_adjectival_feminine_nouns() -> None:
    cells = _cells(
        GEN_SG=["vysoké"],
        DAT_PL=["vysokým"],
        INS_PL=["vysokými"],
        INS_SG=["vysokou"],
    )
    found = ortho.analyse_stem("noun", "vysoká", cells)
    assert (found.stem, found.declension) == ("vysok", ortho.ADJECTIVAL)


def test_ndeb_tokens_are_root_internal() -> None:
    assert ortho.ndeb_tokens("mladý", ["[ˈmladiː]"]) == ()
    assert ortho.ndeb_tokens("vysoký", ["[ˈvɪsokiː]"]) == ()
    assert ortho.ndeb_tokens("gymnázium", ["[ˈɡɪmnaːzɪjum]"]) == ("gy",)
    assert ortho.ndeb_tokens("diplomat", ["[ˈdɪplomat]"]) == ("dɪ",)


# ---------------------------------------------------------------------------
# Names (R25)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("word", "gloss", "name_class"),
    [
        ("Studénka", "město v okrese Nový Jičín", names.TOPONYM),
        ("Anuška", "Anna; ženské křestní jméno", names.GIVEN_NAME),
        ("Jindřišská", "ulice (v českém sídle), ve které …", names.TOPONYM),
        ("Azor", "jméno pro psa (samce)", names.OTHER_NAME),
        ("Liška", "jedno ze souhvězdí, Lištička", names.OTHER_NAME),
        ("Tuareg", "příslušník afrického národa", ""),
        ("kongo", "maskáčová bunda", ""),
    ],
)
def test_capitalized_cs_names(word: str, gloss: str, name_class: str) -> None:
    verdict = names.classify_cs(
        {"word": word, "pos": "noun", "senses": [{"glosses": [gloss]}]}
    )
    assert verdict.name_class == name_class
    assert verdict.is_name == bool(name_class)


def test_capitalized_en_place_and_sky_nouns() -> None:
    def entry(word: str, gloss: str, cats: tuple[str, ...] = ()) -> Any:
        return {
            "word": word,
            "pos": "noun",
            "senses": [{"glosses": [gloss], "categories": list(cats)}],
        }

    assert names.classify_en(entry("Rak", "Cancer (constellation)")).is_name
    assert names.classify_en(
        entry("Býk", "Taurus", ("cs:Constellations in the zodiac",))
    ).is_name
    assert names.classify_en(
        entry("Franky", "Franconia (geographical region in Germany)")
    ).is_name
    demonym = entry("Akkaďan", "Akkadian (inhabitant of the Akkad region)")
    assert not names.classify_en(demonym).is_name
    assert not names.classify_en(entry("rak", "crayfish")).is_name


# ---------------------------------------------------------------------------
# Frequency bases (attest)
# ---------------------------------------------------------------------------


def test_lemma_ambiguity_and_unambiguous_paradigm_sum() -> None:
    from lexpipe.attest import AttestationTable

    info = attestation.LexemeInfo(
        "cs:jak:noun:2:x", "jak", "noun", "cs:jak:r", True, "cs:jak:noun:2:x"
    )
    table = AttestationTable(
        counts_cs={"jak": 1000, "jaka": 3, "jaky": 2},
        counts_folded={"jak": 1000, "jaka": 3, "jaky": 2},
        total_tokens=1_000_000,
    )
    owners = attestation.Owners()
    owners.exact["jak"].add("w:conj:jak")
    owners.folded["jak"].add("w:conj:jak")
    compatible = {info.lexeme_id: frozenset({f"id:{info.base_id}"})}
    row = attestation._frequency_row(
        info,
        table,
        owners,
        compatible,
        {info.lexeme_id: {"jaka": 3, "jaky": 2}},
        {},
    )
    assert row["lemma_ambiguous"] is True
    assert row["lemma_count"] == 1000
    assert row["paradigm_sum"] == 1005
    assert row["paradigm_sum_unambiguous"] == 5


# ---------------------------------------------------------------------------
# TP cells: strict thresholds, animacy, -a masculines, ADV, R17
# ---------------------------------------------------------------------------


def real(
    cell: str,
    form: str,
    ending: str,
    gender: str,
    allomorph: str,
    *,
    mutation: bool | None = False,
    n_kind: str = "",
    ambiguous: bool = False,
    confounded: bool = False,
) -> roots.Realization:
    return roots.Realization(
        cell=cell,
        form=form,
        ending=ending,
        gender=gender,
        allomorph=allomorph,
        n_kind=n_kind,
        mutation=mutation,
        aligned=True,
        count=10,
        attested=ALL,
        ambiguous=ambiguous,
        confounded=confounded,
    )


def record(
    lemma: str,
    realizations: list[roots.Realization],
    *,
    pos: str = "noun",
    readings: tuple[str, ...] = ("m.an",),
    declension: str = ortho.HARD,
    target: str = "velar",
) -> roots.RootRecord:
    return roots.RootRecord(
        root_id=f"cs:{lemma}:x",
        pos=pos,
        lemma=lemma,
        lexeme_ids=[f"cs:{lemma}:{pos}:0:x"],
        target_class=target,
        stem=lemma,
        stem_final=lemma[-1],
        cluster="",
        readings=list(readings),
        suffix_target=False,
        lemma_suffix="",
        nde_class="",
        stratum="inherited",
        lemma_freq=100.0,
        paradigm_freq=200.0,
        realizations=realizations,
        declension=declension,
        unambiguous_freq=50.0,
    )


def test_strict_thresholds() -> None:
    parku = real("LOC.SG", "parku", "u", "m.in", "N", confounded=True)
    assert parku.attested_at("ge1")
    assert not parku.attested_at("ge1_strict")
    flags = parku.attested_by_tau()
    assert flags["ge2pm"] and not flags["ge2pm_strict"]
    assert set(roots.ALL_TAUS) == set(flags)
    assert roots.plain_tau("ge0.5pm_strict") == "ge0.5pm"


def test_frequency_bases() -> None:
    item = record("kluk", [])
    assert item.frequency("lemma") == 100.0
    assert item.frequency("paradigm") == 200.0
    assert item.frequency("paradigm_unambiguous") == 50.0


def test_unknown_animacy_stays_out_of_animate_cells() -> None:
    rows = [
        real("VOC.SG", "lesbizme", "e", "m", "N", n_kind="nontrigger"),
        real("LOC.PL", "lesbizmech", "ech", "m", "N", n_kind="nontrigger"),
    ]
    item = record("lesbizmus", rows, readings=("m",))
    voc = tc.COMPETITION_BY_NAME["VOC.SG m.an"]
    assert tc.competition_observations(item, voc) == ()
    pooled = next(t for t in tc.RULE_TRIGGERS if t.name == "LOC.PL m/n -ech")
    assert len(tc.rule_observations(item, pooled)) == 1


def test_a_masculines_have_no_nom_pl_choice() -> None:
    rows = [real("NOM.PL", "kolegové", "ové", "m.an", "N", n_kind="allomorph")]
    kolega = record("kolega", rows)
    assert kolega.a_masculine
    for name in ("NOM.PL m.an", "NOM.PL m.an -ové"):
        assert (
            tc.competition_observations(kolega, tc.COMPETITION_BY_NAME[name])
            == ()
        )
    assert not record("kluk", rows).a_masculine


def test_adverb_competitors_apart() -> None:
    rows = [
        real("ADV", "lehce", "e", "adj", "P", mutation=True),
        real("ADV", "lehko", "o", "adj", "N", n_kind="allomorph"),
    ]
    lehky = record("lehký", rows, pos="adj", readings=())
    pooled = tc.competition_observations(lehky, tc.COMPETITION_BY_NAME["ADV"])
    only_y = tc.competition_observations(
        lehky, tc.COMPETITION_BY_NAME["ADV -y"]
    )
    only_o = tc.competition_observations(
        lehky, tc.COMPETITION_BY_NAME["ADV -o"]
    )
    assert len(pooled) == 2 and len(only_o) == 2
    assert [o.alternates for o in only_y] == [True]


def test_adjectival_nouns_out_of_noun_r17_rows() -> None:
    rows = [real("DAT.SG", "vysoké", "é", "f", "none")]
    noun = record("vysoká", rows, readings=("f",), declension=ortho.ADJECTIVAL)
    assert tc.r17_observations(noun, frozenset({roots.FRONT_MID}), "E0") == ()


def test_sky_class_comes_before_the_etymological_classes() -> None:
    names_ = [name for name, _ in tpstage.CLASS_ORDER]
    assert names_.index("sky") < names_.index("loan")
    sky = dict(tpstage.CLASS_ORDER)["sky"]
    adjective = record("typický", [], pos="adj", readings=())
    adjective.suffix_target = True
    adjective.lemma_suffix = "cký"
    assert sky(adjective)


# ---------------------------------------------------------------------------
# Derivations: direction, homographs, feminine agents, -icismus, -ie
# ---------------------------------------------------------------------------


def _base(
    lemma: str,
    readings: tuple[str, ...] = ("m.in",),
    *,
    lexeme_id: str = "",
    root_id: str = "",
    stem: str = "",
    target: str = "velar",
    head_links: tuple[tuple[str, str], ...] = (),
) -> derivations.Base:
    return derivations.Base(
        lexeme_id=lexeme_id or f"cs:{lemma}:noun:0:aaaaaaaa",
        root_id=root_id or f"cs:{lemma}:r",
        lemma=lemma,
        readings=readings,
        stem=stem or lemma[:-1],
        target_class=target,
        stem_final="k",
        cluster="",
        suffix_target=False,
        etymology="0",
        head_links=head_links,
    )


def _pair(base: derivations.Base, word: str, *evidence: str) -> Any:
    found = sx.decompose(word, base.stem, base.lemma, allow_prefix=True)
    assert found is not None, word
    return derivations.Pair(base, word, found, set(evidence))


def test_reversed_pairs() -> None:
    utok = _base("útok", stem="útok")
    assert derivations.reversed_pair(
        _pair(utok, "útočit", "wikt_related"), ("útočit",)
    )
    leh = _base("leh", stem="leh")
    assert derivations.reversed_pair(
        _pair(leh, "ležet", "cand_headword"), ("lehat",)
    )
    poloha = _base("poloha", ("f",), stem="poloh")
    # an -ovat verb on the plain stem is formed from the noun
    assert not derivations.reversed_pair(
        _pair(poloha, "polohovat", "cand_headword"), ("položit",)
    )
    # a derivative whose own etymology names the base is never reversed
    assert not derivations.reversed_pair(
        _pair(poloha, "položit", "wikt_ety"), ("položit",)
    )


def test_deverbal_sources() -> None:
    steps = ety.steps(
        [{"name": "deverbal", "args": {"1": "cs", "2": "obložit"}}]
    )
    assert derivations.deverbal_sources(steps, "") == ("obložit",)
    assert derivations.deverbal_sources(
        [], "From the verb tušit (“to suspect”)"
    ) == ("tušit",)
    assert "páčit" in derivations.deverbal_sources(
        [], "Formed from páčil, active participle of páčit (“to lever”)."
    )


def test_homograph_dedupe_prefers_animate_bases_for_feminine_agents() -> None:
    index = derivations.WiktionaryIndex()
    rival = _base(
        "sok",
        ("m.an",),
        lexeme_id="cs:sok:noun:1:aaaaaaaa",
        root_id="cs:sok:1",
        stem="sok",
    )
    ditch = _base(
        "sok",
        ("m.in",),
        lexeme_id="cs:sok:noun:0:aaaaaaaa",
        root_id="cs:sok:0",
        stem="sok",
    )
    pairs = derivations.PairSet()
    for base in (rival, ditch):
        pair = _pair(base, "sokyně", "wikt_ety")
        pairs.pairs[pair.key] = pair
    assert derivations.dedupe_homographs(pairs, index) == 1
    (kept,) = pairs.pairs.values()
    assert kept.base.readings == ("m.an",)


def test_feminine_agents() -> None:
    index = derivations.WiktionaryIndex()
    index.feminine_agents.add("Češka")
    cech = _base("Čech", ("m.an",), stem="Čech")
    assert derivations.feminine_agent(_pair(cech, "Češka", "x"), index)
    deka = _base("deka", ("f",), stem="dek")
    assert not derivations.feminine_agent(_pair(deka, "dečka", "x"), index)


def test_icismus_is_not_k_plus_ismus() -> None:
    assert sx.decompose("katolicismus", "katolík", "katolík") is None
    plain = sx.decompose("tankista", "tank", "tank")
    assert plain is not None and plain.suffix.label == "-ista"


def test_plain_headword_candidates_are_accepted_and_flagged() -> None:
    base = _base("pantofág", ("m.an",), stem="pantofág")
    suffix = sx.BY_LABEL["-ní/-ný"]
    plain = derivations.Found(
        base, "pantofágní", suffix, "cand_headword", False, False
    )
    assert plain.plain_palatalizing
    index = derivations.WiktionaryIndex()
    assert derivations.candidate_problem(plain, index, set()) == ""
    corpus = derivations.Found(
        base, "pantofágní", suffix, "cand_corpus", False, False
    )
    assert derivations.candidate_problem(corpus, index, set()) == (
        "unexpected spelling"
    )
    pairs = derivations.PairSet()
    rejected: Any = __import__("collections").Counter()
    derivations.add_candidates(pairs, [plain], index, rejected)
    (pair,) = pairs.pairs.values()
    assert pair.evidence == {"cand_headword_plain"}
    assert pair.basis() == "wikt"


def test_ista_of_an_ie_noun_is_rejected() -> None:
    index = derivations.WiktionaryIndex()
    index.headwords["ekologie"].add(("en", "noun"))
    base = _base("ekolog", ("m.an",), stem="ekolog")
    found = derivations.Found(
        base, "ekologista", sx.BY_LABEL["-ista"], "cand_headword", True
    )
    assert derivations.candidate_problem(found, index, set()) == (
        "-ista/-ismus of an -ie noun"
    )


# ---------------------------------------------------------------------------
# Strata: adapted loans and hybrids
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("key", "stratum"),
    [
        ("-ický|tragicus", (ety.BORROWED, ety.ADAPTED_LOAN)),
        ("-ace|accumulātiō", (ety.BORROWED, ety.ADAPTED_LOAN)),
        ("-ek|fan", (ety.INTERNAL, ety.HYBRID)),
        ("-an|Welsh", (ety.INTERNAL, ety.HYBRID)),
        ("-ový|mainstream", (ety.INTERNAL, ety.HYBRID)),
        ("-ní|default", (ety.INTERNAL, ety.HYBRID)),
        ("-ýn|Baldachin", (ety.BORROWED, ety.ADAPTED_LOAN)),
        ("-ýn|blond", (ety.INTERNAL, ety.HYBRID)),
        ("-ismus|virīlis", (ety.INTERNAL, ety.HYBRID)),
    ],
)
def test_refine_stratum(key: str, stratum: tuple[str, str]) -> None:
    raw = (ety.BORROWED, ety.INTERNAL_ANALYSIS)
    assert ety.refine_stratum(raw, key, frozenset({"default"})) == stratum


def test_refine_stratum_leaves_other_verdicts() -> None:
    raw = (ety.BORROWED, "bor")
    assert ety.refine_stratum(raw, "-ek|fan", frozenset()) == raw


def test_adaptation_key() -> None:
    steps = ety.steps(
        [{"name": "af", "args": {"1": "cs", "2": "la:tragicus", "3": "-ický"}}]
    )
    assert ety.adaptation_key(steps) == "-ický|tragicus"


# ---------------------------------------------------------------------------
# R20 members
# ---------------------------------------------------------------------------


def _row(**values: str) -> dict[str, str]:
    row = {
        "suffix_pair": "",
        "suffix_class": sx.P_CLASS,
        "outcome": sx.PAL,
        "base_readings": "",
        "derivative_feminine_agent": "False",
        "base_deverbal_from": "",
    }
    row.update(values)
    return row


def test_ka_yne_needs_an_animate_masculine_and_a_feminine_agent() -> None:
    diminutive = _row(suffix_pair=sx.PAIR_KA_YNE, base_readings="f")
    assert spellout.member_side(diminutive) is None
    agent = _row(
        suffix_pair=sx.PAIR_KA_YNE,
        base_readings="m.an",
        derivative_feminine_agent="True",
    )
    assert spellout.member_side(agent) is True
    plant = _row(
        suffix_pair=sx.PAIR_KA_YNE,
        suffix_class=sx.N_CLASS,
        outcome=sx.PLAIN,
        base_readings="m.in",
        derivative_feminine_agent="True",
    )
    assert spellout.member_side(plant) is None


def test_it_ovat_leaves_out_deverbal_bases() -> None:
    row = _row(suffix_pair=sx.PAIR_IT_OVAT, base_deverbal_from="založit")
    assert spellout.member_side(row) is None
    assert spellout.member_side(_row(suffix_pair=sx.PAIR_IT_OVAT)) is True


# ---------------------------------------------------------------------------
# R21: strict thresholds and the internal-formation filter
# ---------------------------------------------------------------------------


def test_noevidence_derivative_flags_ignore_strictness() -> None:
    row = {"attested_ge1": "True", "attested_ge2pm": "False"}
    assert noevidence._attested(row, "ge1_strict")
    assert not noevidence._attested(row, "ge2pm_strict")


def test_noevidence_internal_formation_filter() -> None:
    item = record("sníh", [], readings=("m.in",))
    records = {item.key: item}
    lexeme_key = {item.lexeme_ids[0]: item.key}

    def pair(word: str, stratum: str) -> dict[str, str]:
        return {
            "base_lexeme_id": item.lexeme_ids[0],
            "basis": "linked",
            "suffix_class": sx.P_CLASS,
            "outcome": sx.PAL,
            "suffix": "-ní/-ný",
            "derivative": word,
            "base_stratum_wikt": "inherited",
            "derivative_stratum_wikt": stratum,
            "attested_ge1": "True",
            "attested_ge0.5pm": "True",
            "attested_ge2pm": "True",
        }

    rows = [pair("sněžný", "inherited"), pair("sněžnicový", "unknown")]
    every = noevidence.units(records, lexeme_key, rows, "ge1", "wikt")
    internal = noevidence.units(
        records,
        lexeme_key,
        rows,
        "ge1",
        "wikt",
        derivative_strata=frozenset({ety.INTERNAL}),
    )
    assert len(every) == 1 and internal == []


# ---------------------------------------------------------------------------
# R27 classifier
# ---------------------------------------------------------------------------


def test_r27_czech_glosses_use_the_czech_rules_only() -> None:
    gloss = "soubor deseti božích přikázání předaných Mojžíšovi na Sinaji"
    assert logtable.semantic_class(gloss, None, "cs")[0] == "thing"


def test_r27_english_rules_read_the_gloss_head() -> None:
    gloss = "section (part of a legal document, marked using section sign)"
    assert logtable.gloss_head(gloss) == "section"
    assert logtable.semantic_class(gloss, None, "en")[0] != "person"
    assert logtable.semantic_class("biologist", None, "en")[0] == "person"


def test_core_counts_strict_rows() -> None:
    rows = [
        real("LOC.SG", "parku", "u", "m.in", "N", confounded=True),
        real("LOC.SG", "parce", "e", "m.in", "P", mutation=True),
    ]
    park = record("park", rows, readings=("m.in",))
    comp = tc.COMPETITION_BY_NAME["LOC.SG m.in"]
    observed = tpstage.observe([park], tpstage.competition_triggers(comp))
    plain = tp.CellSpec("cs", "competition", "velar", comp.key, tau="ge1")
    strict = tp.CellSpec(
        "cs", "competition", "velar", comp.key, tau="ge1_strict"
    )
    assert tp.root_status(observed[0], comp.key, "ge1", "dominant") in (
        tp.RootStatus.ALTERNATES,
        tp.RootStatus.VARIABLE,
        tp.RootStatus.NONALT,
    )
    assert plain.tau != strict.tau
    assert (
        tp.root_status(observed[0], comp.key, "ge1_strict", "dominant")
        == tp.RootStatus.ALTERNATES
    )


def test_related_terms_conflicts() -> None:
    index = derivations.WiktionaryIndex()
    index.derivative_bases["vystrašit"].add("strašit")
    index.add_stratum("křičet", "verb", "0", (ety.INHERITED, "inh"))
    index.add_stratum("strašný", "adj", "0", (ety.INHERITED, "inh"))
    index.add_stratum("papoušek", "noun", "0", (ety.BORROWED, "der:ar"))
    assert index.related_conflict("vystrašit", "verb", "strach", True)
    assert index.related_conflict("křičet", "verb", "křik", True)
    assert not index.related_conflict("strašný", "adj", "strach", False)
    assert index.related_conflict("papoušek", "noun", "papouch", False)
    assert not index.related_conflict("strašnice", "noun", "strach", False)
