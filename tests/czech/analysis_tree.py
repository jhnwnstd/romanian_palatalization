"""A fixture tree for the Czech analysis stages.

Builds on :mod:`tests.czech.treebuilder` (the data-layer tree) and adds
what the analysis stages read:

  - the derivation fixtures (``tests/fixtures/czech/*_derivations.jsonl``,
    CC BY-SA Wiktionary lines) appended to the two kaikki inputs;
  - a **synthetic** DeriNet 2.3 TSV and czetyl TSV (made-up rows in the
    real formats; no CC BY-NC-SA data is tracked);
  - extra **synthetic** Leipzig counts for derivative forms.

Every input can be shuffled (determinism tests). The stage parameters
keep the fixture runs small: few simulated learners, one worker, V
values the tiny lexicon can supply, few permutations.
"""

from __future__ import annotations

import gzip
import io
import json
import random
import sys
import tarfile
from collections.abc import Sequence
from pathlib import Path
from typing import Final

PROJECT_ROOT: Final[Path] = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lexpipe import spec  # noqa: E402
from lexpipe.context import RunContext  # noqa: E402
from lexpipe.langs import czech  # noqa: E402
from tests.czech import treebuilder  # noqa: E402

FIXTURES: Final[Path] = treebuilder.FIXTURES
EN_DERIVATIONS: Final[Path] = FIXTURES / "en_wiktionary_cs_derivations.jsonl"
CS_DERIVATIONS: Final[Path] = FIXTURES / "cs_extract_derivations.jsonl"
DATA_STAGES: Final[tuple[str, ...]] = treebuilder.STAGES
ANALYSIS_STAGES: Final[tuple[str, ...]] = (
    "derivations",
    "tp",
    "spellout",
    "noevidence",
    "covariation",
    "logtable",
    "docs",
)
PARAMS: Final[dict[str, object]] = {
    **treebuilder.PARAMS,
    "n_samples": 10,
    "workers": 1,
    "v_values": [2, 5],
    "n_permutations": 200,
}

# Synthetic Leipzig counts of derivative forms (made up).
EXTRA_COUNTS: Final[tuple[tuple[str, int], ...]] = (
    ("knižní", 200),
    ("knižního", 18),
    ("knížka", 250),
    ("ruční", 60),
    ("ručit", 30),
    ("ligový", 15),
    ("ligové", 900),
    ("bankovní", 950),
    ("drogový", 200),
    ("vzdušný", 38),
    ("vzduchový", 50),
    ("strašný", 330),
    ("tankista", 20),
    ("šachista", 13),
    ("bioložka", 26),
    ("žačka", 5),
    ("žákyně", 5),
    ("parkovat", 90),
    ("klouček", 40),
    ("vzduchu", 700),
    ("strachu", 300),
    ("liga", 2000),
    ("banka", 1200),
    ("bance", 200),
    ("droga", 300),
    ("droze", 10),
)


def _node(
    node_id: str,
    lemma: str,
    tag: str,
    pos: str,
    feats: str,
    morphs: Sequence[tuple[str, str]],
    parent: str = "",
    relation: str = "",
    other: dict[str, object] | None = None,
) -> str:
    segmentation = (
        json.dumps(
            [{"morph": m, "type": t} for m, t in morphs], ensure_ascii=False
        )
        if morphs
        else ""
    )
    extra = json.dumps(other or {}, ensure_ascii=False)
    return "\t".join(
        [
            node_id,
            f"{lemma}#{tag}",
            lemma,
            pos,
            feats,
            segmentation,
            parent,
            relation,
            "",
            extra,
        ]
    )


DER: Final[str] = "Type=Derivation"
N_F: Final[str] = "NNF??-----A---?"
N_MI: Final[str] = "NNI??-----A---?"
N_MA: Final[str] = "NNM??-----A---?"
ADJ_T: Final[str] = "AA???----??---?"
VERB_T: Final[str] = "V????--??-??I??"

# Synthetic DeriNet trees (made up; the real file is CC BY-NC-SA).
DERINET_TREES: Final[tuple[tuple[str, ...], ...]] = (
    (
        _node("1.0", "-log-", "Affixoid", "Affixoid", "Fictitious=Yes", ()),
        _node("1.1", "kniha", N_F, "NOUN", "Gender=Fem", (), "1.0", DER),
    ),
    (
        _node(
            "2.0",
            "kniha",
            N_F,
            "NOUN",
            "Gender=Fem&Loanword=False",
            (("knih", "R"), ("a", "I")),
            other={"unmotivated": True, "techlemma": "kniha"},
        ),
        _node(
            "2.1",
            "knižní",
            ADJ_T,
            "ADJ",
            "Loanword=False",
            (("kniž", "R"), ("n", "D"), ("í", "I")),
            "2.0",
            DER,
        ),
        _node(
            "2.2",
            "knihovat",
            VERB_T,
            "VERB",
            "Aspect=Imp&ConjugClass=3&Loanword=False",
            (("knih", "R"), ("ova", "D"), ("t", "I")),
            "2.0",
            DER,
        ),
    ),
    (
        _node(
            "3.0",
            "liga",
            N_F,
            "NOUN",
            "Gender=Fem&Loanword=True",
            (("liga", "R"),),
            other={"Etymology": ["deu"], "unmotivated": True},
        ),
        _node(
            "3.1",
            "ligový",
            ADJ_T,
            "ADJ",
            "Loanword=True",
            (("lig", "R"), ("ov", "D"), ("ý", "I")),
            "3.0",
            DER,
        ),
    ),
    (
        _node(
            "4.0",
            "ruka",
            N_F,
            "NOUN",
            "Gender=Fem&Loanword=False",
            (("ruk", "R"), ("a", "I")),
        ),
        _node(
            "4.1",
            "ruční",
            ADJ_T,
            "ADJ",
            "Loanword=False",
            (("ruč", "R"), ("n", "D"), ("í", "I")),
            "4.0",
            DER,
        ),
        _node(
            "4.2",
            "ručit",
            VERB_T,
            "VERB",
            "Aspect=Imp&ConjugClass=4&Loanword=False",
            (("ruč", "R"), ("i", "D"), ("t", "I")),
            "4.0",
            DER,
        ),
    ),
    (
        _node(
            "5.0",
            "šach",
            N_MI,
            "NOUN",
            "Animacy=Inan&Gender=Masc&Loanword=True",
            (("šach", "R"),),
            other={"unmotivated": True},
        ),
        _node(
            "5.1",
            "šachista",
            N_MA,
            "NOUN",
            "Animacy=Anim&Gender=Masc&Loanword=True",
            (("šach", "R"), ("ist", "D"), ("a", "I")),
            "5.0",
            DER,
        ),
    ),
    (
        _node(
            "6.0",
            "šach",
            N_MA,
            "NOUN",
            "Animacy=Anim&Gender=Masc&Loanword=True",
            (("šach", "R"),),
        ),
    ),
    (
        _node(
            "7.0",
            "vzduch",
            N_MI,
            "NOUN",
            "Animacy=Inan&Gender=Masc&Loanword=True",
            (("vzduch", "R"),),
            other={"Etymology": ["rus"], "unmotivated": True},
        ),
        _node(
            "7.1",
            "vzdušný",
            ADJ_T,
            "ADJ",
            "Loanword=False",
            (("vzduš", "R"), ("n", "D"), ("ý", "I")),
            "7.0",
            DER,
        ),
        _node(
            "7.2",
            "vzduchový",
            ADJ_T,
            "ADJ",
            "Loanword=False",
            (("vzduch", "R"), ("ov", "D"), ("ý", "I")),
            "7.0",
            DER,
        ),
    ),
    (
        _node(
            "8.0",
            "strach",
            N_MI,
            "NOUN",
            "Animacy=Inan&Gender=Masc&Loanword=False",
            (("strach", "R"),),
        ),
        _node(
            "8.1",
            "strašit",
            VERB_T,
            "VERB",
            "Aspect=Imp&ConjugClass=4&Loanword=False",
            (("straš", "R"), ("i", "D"), ("t", "I")),
            "8.0",
            DER,
        ),
        _node(
            "8.2",
            "strašný",
            ADJ_T,
            "ADJ",
            "Loanword=False",
            (("straš", "R"), ("n", "D"), ("ý", "I")),
            "8.1",
            DER,
        ),
    ),
)

# Synthetic czetyl rows (word, donor chain, loan|native); a-2/a disagree.
CZETYL_ROWS: Final[tuple[tuple[str, str, str], ...]] = (
    ("a", "", "native"),
    ("a-2", "ell", "loan"),
    ("kniha", "", "native"),
    ("liga", "deu", "loan"),
    ("banka", "ita,gem", "loan"),
    ("vzduch", "rus", "loan"),
    ("šach", "gmh,dum", "loan"),
    ("park", "fra,lat", "loan"),
)


def _shuffled(lines: list[bytes], seed: int | None) -> list[bytes]:
    if seed is not None:
        random.Random(seed).shuffle(lines)
    return lines


def _derinet_bytes(seed: int | None) -> bytes:
    trees = list(DERINET_TREES)
    if seed is not None:
        random.Random(seed).shuffle(trees)
    blocks = ["\n".join(tree) for tree in trees]
    return ("\n\n".join(blocks) + "\n\n").encode("utf-8")


def _czetyl_bytes(seed: int | None) -> bytes:
    rows = _shuffled(
        [f"{w}\t{e}\t{loan}\n".encode() for w, e, loan in CZETYL_ROWS], seed
    )
    return b"word\tetymology\tis_loan\n" + b"".join(rows)


def _leipzig_bytes(seed: int | None) -> bytes:
    counts = [*treebuilder.LEIPZIG_COUNTS, *EXTRA_COUNTS]
    corpus = "ces_news_2024_1M"
    rows = _shuffled(
        [f"{i}\t{w}\t{f}\n".encode() for i, (w, f) in enumerate(counts, 1)],
        seed,
    )
    members = {
        f"{corpus}/{corpus}-words.txt": b"".join(rows),
        f"{corpus}/{corpus}-meta.txt": b"1\tsynthetic test corpus\n",
    }
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as tar:
        for name, data in sorted(members.items()):
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return buffer.getvalue()


def build_tree(root: Path, *, shuffle_seed: int | None = None) -> Path:
    """The data-layer tree plus the analysis inputs (module docstring)."""
    treebuilder.build_tree(root, shuffle_seed=shuffle_seed)
    sources = root / "data" / "sources"
    en_lines = _shuffled(
        [
            *treebuilder.EN_FIXTURE.read_bytes().splitlines(keepends=True),
            *EN_DERIVATIONS.read_bytes().splitlines(keepends=True),
        ],
        shuffle_seed,
    )
    (sources / "kaikki" / "en_wiktionary_cs.jsonl").write_bytes(
        b"".join(en_lines)
    )
    cs_lines = _shuffled(
        [
            *treebuilder.CS_FIXTURE.read_bytes().splitlines(keepends=True),
            *CS_DERIVATIONS.read_bytes().splitlines(keepends=True),
        ],
        shuffle_seed,
    )
    (sources / "kaikki" / "cs-extract.jsonl.gz").write_bytes(
        gzip.compress(b"".join(cs_lines), mtime=0)
    )
    (sources / "leipzig" / "ces_news_2024_1M.tar.gz").write_bytes(
        _leipzig_bytes(shuffle_seed)
    )
    derinet = sources / "derinet" / "derinet-2-3.tsv"
    derinet.parent.mkdir(parents=True, exist_ok=True)
    derinet.write_bytes(_derinet_bytes(shuffle_seed))
    czetyl = sources / "czetyl" / "czetyl.tsv"
    czetyl.parent.mkdir(parents=True, exist_ok=True)
    czetyl.write_bytes(_czetyl_bytes(shuffle_seed))
    return root


def context(root: Path) -> RunContext:
    """The RunContext of an analysis fixture run."""
    return RunContext(spec.LANGUAGES["cs"], root=root, params=PARAMS)


def run_all(root: Path) -> None:
    """Run the data-layer and analysis stages in order."""
    ctx = context(root)
    for name in (*DATA_STAGES, *ANALYSIS_STAGES):
        status = czech.STAGES[name](ctx)
        if isinstance(status, int) and status != 0:
            raise RuntimeError(f"stage {name} failed with {status}")
