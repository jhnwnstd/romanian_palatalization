"""The lexicon stage's stress, domain and root-key rules (R4, R8, R26).

- stress from hyphenations: syllables split at the separators editors
  use, a hyphenation that spells another word or accents the i of an
  -ico suffix is not used (``stress_rejected``);
- ``ico_stressed``: which vowel of a penultimate-stressed -ico word
  carries the stress (amìco vs arcaico /arˈkaj.ko/);
- the ``plurale_tantum`` and ``plural_form`` domain classes;
- it.wiktionary forms: superlatives are not plurals;
- root keys: different source words stay apart (arca, arco), POS
  conversions without a stress of their own join their partner.
"""

from __future__ import annotations

from typing import Any

import pytest

from lexpipe.langs.italian import etymology, lexicon, plurals

pytestmark = pytest.mark.it


def _t(name: str, **args: str) -> dict[str, Any]:
    return {"name": name, "args": args}


# ---------------------------------------------------------------------------
# Hyphenations
# ---------------------------------------------------------------------------


def test_syllables_split_at_every_separator() -> None:
    assert lexicon.hyphenation_syllables(["post - so - vi - è - ti - co"]) == [
        "post",
        "so",
        "vi",
        "è",
        "ti",
        "co",
    ]
    assert lexicon.hyphenation_syllables(["ra", "dio", "e", "stè", "si co"])[
        -2:
    ] == ["si", "co"]
    assert lexicon.hyphenation_syllables(["bat-ti-scó-pa"]) == [
        "bat",
        "ti",
        "scó",
        "pa",
    ]
    # a footnote number is not a syllable
    assert lexicon.hyphenation_syllables(["ri", "le", "gà", "to 1"]) == [
        "ri",
        "le",
        "gà",
        "to",
    ]
    assert lexicon.syllable_stress(["post - so - vi - è - ti - co"]) == 3
    assert lexicon.syllable_stress(["ri", "le", "gà", "to 1"]) == 2


def test_a_hyphenation_must_spell_its_word() -> None:
    assert lexicon.spells("amico", ["a", "mì", "co"])
    assert lexicon.spells("amuleto", ["a", "mu", "lè", "to'"])
    # another final vowel keeps the stress position (agonistico on
    # agonistica)
    assert lexicon.spells("agonistica", ["a", "go", "nì", "sti", "co"])
    assert not lexicon.spells("biocenotico", ["bi", "o", "ce", "nò", "si"])
    assert not lexicon.spells("neonatologo", "ne o na to lo gì a".split())
    assert lexicon.hyphenation_stress("bi‧o‧ce‧nò‧si", "biocenotico") is None
    assert lexicon.hyphenation_stress("a‧mì‧co", "amico") == 2


def test_accented_vowel_and_ico_stress() -> None:
    assert lexicon.accented_vowel(["a", "mì", "co"]) == "i"
    assert lexicon.accented_vowel(["lài", "co"]) == "a"
    assert lexicon.accented_vowel(["a", "mi", "co"]) == ""
    assert lexicon.ipa_vowel_of("ˈkaj") == "a"
    assert lexicon.ipa_vowel_of("ˈmi") == "i"
    # penultimate stress on the i of -ico: amìco
    assert lexicon.ico_stressed("amico", 2, "i") is True
    # penultimate stress on the vowel before it: arcaico, laico
    assert lexicon.ico_stressed("arcaico", 2, "a") is False
    assert lexicon.ico_stressed("medico", 3, "e") is False
    assert lexicon.ico_stressed("pubblico", None, "") is None
    assert lexicon.ico_stressed("parco", 2, "a") is None  # not -ico


def test_stress_of_reads_the_stressed_vowel_of_the_ipa() -> None:
    def record(word: str, ipa: str) -> dict[str, Any]:
        return {"word": word, "sounds": [{"ipa": ipa}]}

    arcaico = lexicon.stress_of(record("arcaico", "/arˈkaj.ko/"))
    assert (arcaico.stress_from_end, arcaico.stressed_vowel) == (2, "a")
    assert lexicon.ico_stressed("arcaico", 2, arcaico.stressed_vowel) is False
    amico = lexicon.stress_of(record("amico", "/aˈmi.ko/"))
    assert (amico.stress_from_end, amico.stressed_vowel) == (2, "i")


def test_hyphenation_choice_rejections() -> None:
    template = etymology.steps(
        [_t("suffix", **{"1": "it", "2": "iperone", "3": "ico"})]
    )
    assert lexicon.hyphenation_choice(
        "iperonico", "i pe ro nì co".split(), "it_hyphenation", template
    ) == (None, "", "it_hyphenation:" + lexicon.REJECT_ICO_ACCENT)
    # an adjective with an attested base (bosone): the same rejection
    lemmas = {"bosone": False, "bosonico": True, "naso": False}
    assert lexicon.hyphenation_choice(
        "bosonico", "bo so nì co".split(), "it_hyphenation", (), lemmas
    )[2].endswith(lexicon.REJECT_ICO_ACCENT)
    # a simplex noun keeps its stressed -ìca (nasìca < Latin nāsīca)
    lemmas["nasica"] = False
    assert lexicon.hyphenation_choice(
        "nasica", "na sì ca".split(), "it_hyphenation", (), lemmas
    ) == (2, "i", "")
    assert lexicon.hyphenation_choice(
        "biocenotico", "bi o ce nò si".split(), "it_hyphenation"
    ) == (None, "", "it_hyphenation:" + lexicon.REJECT_OTHER_WORD)
    assert lexicon.hyphenation_choice(
        "arcaico", "ar cà i co".split(), "it_hyphenation"
    ) == (3, "a", "")


def test_it_fallback_records_its_rejection() -> None:
    edition = lexicon.it_edition(
        [
            {
                "word": "biocenotico",
                "hyphenations": [["bi", "o", "ce", "nò", "si"]],
            }
        ]
    )
    found = lexicon.with_it_fallback(
        lexicon.Stress(), edition, word="biocenotico"
    )
    assert found.stress_from_end is None
    assert found.rejected == ("it_hyphenation:other_word",)
    good = lexicon.it_edition(
        [{"word": "slargo", "hyphenations": [["slàr-go"]]}]
    )
    used = lexicon.with_it_fallback(lexicon.Stress(), good, word="slargo")
    assert (used.stress_from_end, used.stress_source) == (2, "it_hyphenation")


def test_it_edition_drops_superlatives() -> None:
    edition = lexicon.it_edition(
        [
            {
                "word": "timido",
                "forms": [
                    {"form": "timidi", "tags": ["masculine", "plural"]},
                    {
                        "form": "timidissimi",
                        "tags": [
                            "absolute",
                            "superlative",
                            "masculine",
                            "plural",
                        ],
                    },
                ],
            }
        ]
    )
    assert edition.forms == ("m.pl:timidi",)


# ---------------------------------------------------------------------------
# Domain
# ---------------------------------------------------------------------------


def test_plural_entries_are_out_of_domain() -> None:
    alimenti = {
        "word": "alimenti",
        "senses": [
            {
                "glosses": ["plural of alimento"],
                "form_of": ["alimento"],
                "tags": ["form-of", "plural"],
            },
            {"glosses": ["alimony"]},
        ],
    }
    assert lexicon.is_plural_form_entry(alimenti)
    crostini = {
        "word": "crostini",
        "senses": [{"glosses": ["crostini (plural of crostino)"]}],
    }
    assert lexicon.is_plural_form_entry(crostini)
    amico = {"word": "amico", "senses": [{"glosses": ["friend"]}]}
    assert not lexicon.is_plural_form_entry(amico)
    assert lexicon.domain_flags(alimenti, ()) == ["plural_form"]
    brache = plurals.resolve_plurals(
        "brache", [{"name": "it-noun", "args": {"1": "f-p"}}]
    )
    assert "plurale_tantum" in lexicon.domain_flags(
        {"word": "brache", "senses": [{"glosses": ["trousers"]}]}, (), brache
    )


# ---------------------------------------------------------------------------
# Root keys
# ---------------------------------------------------------------------------


def _lexeme(
    word: str,
    pos: str,
    *,
    etymon: str = "",
    stress: str = "",
    gender: str = "m",
    counterparts: tuple[str, ...] = (),
    ety: str = "0",
) -> lexicon.Lexeme:
    found = (
        etymology.steps([_t("inh", **{"1": "it", "2": "la", "3": etymon})])
        if etymon
        else []
    )
    return lexicon.Lexeme(
        row={"word": word, "pos": pos, "lexeme_id": f"it:{word}:{pos}:{ety}"},
        plural_rows=[],
        drift=None,
        stem_key=word[:-1],
        etymon_key=etymology.etymon_key(found),
        stress=stress,
        ety_label=ety,
        counterparts=counterparts,
        gender=gender,
        vowel_class=lexicon.final_vowel_class(word),
    )


def _roots(*lexemes: lexicon.Lexeme) -> dict[str, str]:
    lexicon.assign_root_ids("it", lexemes)
    return {
        f"{lex.row['word']}:{lex.row['pos']}": str(lex.row["root_id"])
        for lex in lexemes
    }


def test_root_keys_keep_different_words_apart() -> None:
    assert lexicon.final_vowel_class("amica") == "o/a"
    assert lexicon.final_vowel_class("dolce") == "e"
    assert lexicon.final_vowel_class("città") == ""
    roots = _roots(
        _lexeme("arca", "noun", etymon="arca", stress="2", gender="f"),
        _lexeme("arco", "noun", etymon="arcus", stress="2"),
        _lexeme("lattice", "noun", stress="3"),
        _lexeme("lattico", "adj", stress="3"),
        _lexeme(
            "amica",
            "noun",
            etymon="amīca",
            stress="2",
            gender="f",
            counterparts=("amico",),
        ),
        _lexeme("amico", "noun", etymon="amīcus", stress="2"),
        _lexeme("amico", "adj", etymon="amīcus", stress="2"),
    )
    assert roots["arca:noun"] != roots["arco:noun"]
    assert roots["lattice:noun"] != roots["lattico:adj"]
    assert roots["amica:noun"] == roots["amico:noun"] == roots["amico:adj"]


def test_a_pos_conversion_without_stress_joins_its_partner() -> None:
    adjective = _lexeme("psicotico", "adj", stress="3")
    noun = _lexeme("psicotico", "noun")
    other = _lexeme("antico", "adj", stress="2", ety="1")
    shifted = _lexeme("antico", "adj", stress="3", ety="1")
    bare = _lexeme("antico", "noun", ety="1")
    filled = lexicon.conversion_stresses(
        [adjective, noun, other, shifted, bare]
    )
    assert filled == {id(noun): "3"}  # antico: two stresses, nothing filled
    roots = _roots(adjective, noun)
    assert roots["psicotico:adj"] == roots["psicotico:noun"]
