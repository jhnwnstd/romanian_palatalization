"""Data-layer rules behind the audit fixes (features, attest, verbs,
derivations).

- ``nde:psiche`` needs a velar singular (not brioche, ceviche);
- the τ flags need a token in the form's own case (R25: Bocchi);
- the competing spelling of a dorsal -i plural (amici ↔ amichi);
- the E1 link prefers roots with etymological evidence (pescare ↔
  pesca 'fishing', not pesca 'peach').
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import pytest

from lexpipe.attest import AttestationTable
from lexpipe.langs.italian import attestation, features, leipzig, verbs

pytestmark = pytest.mark.it


def _psiche_row(word: str, ipa: str, plural: str = "") -> dict[str, str]:
    return {
        "word": word,
        "ipa": ipa,
        "pl_forms": plural,
        "mpl_forms": "",
        "fpl_forms": "",
        "domain": "in_domain",
    }


def test_psiche_needs_a_velar_singular() -> None:
    assert features.velar_che(_psiche_row("psiche", "ˈpsi.ke"))
    assert features.velar_che(_psiche_row("sineddoche", "siˈnɛd.do.ke"))
    assert not features.velar_che(_psiche_row("brioche", "briˈɔʃ"))
    assert not features.velar_che(_psiche_row("ceviche", "seˈvi.t͡ʃe"))
    # without a transcription a velar plural decides
    assert features.velar_che(_psiche_row("nuraghe", "", "nuraghi"))
    assert not features.velar_che(_psiche_row("cloche", "", ""))
    assert features._stem(_psiche_row("brioche", "briˈɔʃ")).target == ""
    assert features._stem(_psiche_row("psiche", "ˈpsi.ke")).nde_class == (
        "psiche"
    )


def _table(tmp_path: Path, counts: Mapping[str, int]) -> AttestationTable:
    """A table from a Leipzig-layout words file (made-up counts)."""
    path = tmp_path / "words.txt"
    path.write_text(
        "".join(
            f"{index}\t{word}\t{count}\n"
            for index, (word, count) in enumerate(counts.items(), 1)
        ),
        encoding="utf-8",
    )
    return AttestationTable.from_words_file(path)


def test_forms_seen_only_capitalized_are_not_attested(tmp_path: Path) -> None:
    table = _table(
        tmp_path,
        {"il": 2_000_000, "Bocchi": 78, "amici": 30, "Amici": 3, "Pochi": 2},
    )
    assert leipzig.attested_flags(table, "bocchi") == {
        "dictionary": True,
        "ge1": False,
        "ge0.5pm": False,
        "ge2pm": False,
    }
    amici = leipzig.attested_flags(table, "amici")
    assert amici["ge1"] and amici["ge2pm"]  # 33 folded tokens
    assert table.attested_by_tau("bocchi")["ge1"]  # folded alone: a name
    paradigm = leipzig.paradigm_flags(table, ["bocchi", "amici"])
    assert paradigm["ge1"]
    assert not leipzig.paradigm_flags(table, ["bocchi", "pochi"])["ge1"]


def _ref(word: str, form: str, cell: str = "pl") -> attestation.FormRef:
    return attestation.FormRef("x", "r", word, "noun", cell, form, "lemma")


def test_competing_spelling_of_dorsal_plurals() -> None:
    assert attestation.competitor(_ref("amico", "amici")) == "amichi"
    assert attestation.competitor(_ref("parco", "parchi")) == "parci"
    assert attestation.competitor(_ref("lago", "laghi")) == "lagi"
    assert attestation.competitor(_ref("bosco", "boschi")) == "bosci"
    assert attestation.competitor(_ref("amica", "amiche", "f.pl")) == ""
    assert attestation.competitor(_ref("luce", "luci")) == ""  # frontstem
    assert attestation.competitor(_ref("occhio", "occhi")) == ""  # NDE
    assert attestation.competitor(_ref("amico", "amico", "lemma")) == ""


def _verb(word: str, named: tuple[str, ...] = ()) -> verbs.VerbLemma:
    return verbs.VerbLemma(
        verb_id=f"it:{word}",
        word=word,
        verb_key=word,
        reflexive=False,
        canonical="",
        ending=word[-3:],
        stem_written=word[:-3],
        stem_key=word[:-3],
        root_stressed_category=False,
        stratum="unknown",
        named=named,
    )


def test_e1_links_prefer_roots_with_etymological_evidence() -> None:
    fishing = verbs.RootLink("it:pesc:1", {"pesca"}, {"pescare"})
    peach = verbs.RootLink("it:pesc:2", {"pesca", "pesco"}, set())
    assert verbs.linked_roots(_verb("pescare"), [fishing, peach]) == [
        "it:pesc:1"
    ]
    dry = verbs.RootLink("it:secc:1", {"secco"}, set())
    shoal = verbs.RootLink("it:secc:2", {"secca"}, set())
    assert verbs.linked_roots(_verb("seccare", ("secco",)), [dry, shoal]) == [
        "it:secc:1"
    ]
    # no evidence for any: every root stays linked; one root: linked
    assert verbs.linked_roots(_verb("giacere"), [dry, shoal]) == [
        "it:secc:1",
        "it:secc:2",
    ]
    assert verbs.linked_roots(_verb("amicare"), [dry]) == ["it:secc:1"]
