"""Repository guards of the multilingual extension (R1, R2, §8).

- Import isolation (R1): importing lexpipe, every lexpipe module and
  every language plugin under lexpipe.langs leaves the Romanian
  pipeline modules out of ``sys.modules``. Runs in a fresh
  interpreter, so modules other tests imported do not count.
- Leipzig frequency directory (R2): ``data/leipzig/freq/`` holds
  exactly the Romanian table, with its pinned SHA-256;
  ``scripts/attach_frequencies.py`` attaches every ``*_freq.csv``
  there. Gitignored data: skipped when the directory is absent.
- Tracked paths: ``git check-ignore`` (patterns only, ``--no-index``)
  must not ignore the files the design tracks (docs, data licence,
  feature inventories, curated decisions, per-language analysis
  and raw lexicons), and must keep ignoring generated and
  NC-licensed data. Skipped without git.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SRC = PROJECT_ROOT / "src"
sys.path.insert(0, str(SRC))

from lexpipe import spec  # noqa: E402

FORBIDDEN_MODULES = (
    "wiktionary_normalizer",
    "romanian_processor_lib",
    "dex_utils",
    "romanian_harvester",
)
FREQ_DIR = PROJECT_ROOT / "data" / "leipzig" / "freq"
FREQ_FILE = "ron_news_2024_1M_freq.csv"
FREQ_SHA256 = (
    "790b8078a525485c0a631ae8bcc315e4a148f48ea6f0309c9dd2eaeb4a7a98e8"
)

_IMPORT_EVERYTHING = """
import importlib, json, pkgutil, sys
sys.path.insert(0, sys.argv[1])
import lexpipe
imported = ["lexpipe"]
for info in pkgutil.walk_packages(lexpipe.__path__, "lexpipe."):
    importlib.import_module(info.name)
    imported.append(info.name)
forbidden = [name for name in json.loads(sys.argv[2]) if name in sys.modules]
print(json.dumps({"imported": imported, "forbidden": forbidden}))
"""


def _import_everything(
    cwd: Path, forbidden: tuple[str, ...]
) -> dict[str, list[str]]:
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            _IMPORT_EVERYTHING,
            str(SRC),
            json.dumps(forbidden),
        ],
        cwd=cwd,
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )
    assert result.returncode == 0, result.stderr[-4000:]
    report: dict[str, list[str]] = json.loads(
        result.stdout.strip().splitlines()[-1]
    )
    return report


def test_lexpipe_never_imports_the_romanian_pipeline(tmp_path: Path) -> None:
    report = _import_everything(tmp_path, FORBIDDEN_MODULES)
    assert report["forbidden"] == []
    imported = set(report["imported"])
    core = {
        "lexpipe.attest",
        "lexpipe.context",
        "lexpipe.csvio",
        "lexpipe.frequency",
        "lexpipe.ids",
        "lexpipe.kaikki",
        "lexpipe.langs",
        "lexpipe.learners",
        "lexpipe.licensing",
        "lexpipe.normalize",
        "lexpipe.sources",
        "lexpipe.spec",
        "lexpipe.stats",
        "lexpipe.tp",
    }
    assert core <= imported
    # Every plugin package on disk is imported, whether registered or not.
    plugins = {
        f"lexpipe.langs.{path.name}"
        for path in (SRC / "lexpipe" / "langs").iterdir()
        if (path / "__init__.py").is_file()
    }
    assert plugins <= imported


def test_import_isolation_check_detects_imports(tmp_path: Path) -> None:
    # Negative control: a module lexpipe does import is reported.
    report = _import_everything(tmp_path, ("json", "no_such_module"))
    assert report["forbidden"] == ["json"]


@pytest.mark.needs_data
@pytest.mark.ro
def test_romanian_leipzig_freq_dir_is_pinned() -> None:
    if not FREQ_DIR.is_dir():
        pytest.skip(f"gitignored directory absent: {FREQ_DIR}")
    entries = sorted(path.name for path in FREQ_DIR.iterdir())
    assert entries == [FREQ_FILE], (
        f"data/leipzig/freq/ holds {entries}: attach_frequencies.py adds a "
        "column per *_freq.csv there (R2). Non-Romanian tables belong in "
        "data/<language>/leipzig/freq/."
    )
    digest = hashlib.sha256((FREQ_DIR / FREQ_FILE).read_bytes()).hexdigest()
    assert digest == FREQ_SHA256


def _language_paths() -> tuple[list[str], list[str]]:
    tracked: list[str] = []
    ignored: list[str] = []
    for info in spec.LANGUAGES.values():
        if info.comparison_only:
            # No docs, curated decisions, analysis dir, raw lexicon or
            # inventory of its own: it writes analysis/comparison/
            # (tracked, MUST_TRACK) and data/<name>/ (ignored).
            tracked.append(f"src/lexpipe/langs/{info.name}/__init__.py")
            ignored += [
                f"data/{info.name}/{info.file_name('roots', '.jsonl')}",
                f"data/{info.name}/{info.file_name('lexicon_raw')}",
                f"data/{info.name}/manifests/tp.manifest.tsv",
            ]
            continue
        tracked += [
            f"docs/{info.name}/data_dict.md",
            f"curated/{info.name}/overrides.csv",
            f"analysis/{info.name}/tp_long.csv",
            f"data/{info.name}/{info.file_name('lexicon_raw')}",
            f"src/phonology/inventories/{info.name}_features.json",
            f"src/lexpipe/langs/{info.name}/__init__.py",
        ]
        ignored += [
            f"data/{info.name}/nc/worklist.csv",
            f"data/{info.name}/manifests/tp.manifest.tsv",
            f"data/{info.name}/{info.file_name('plural_forms')}",
            f"data/{info.name}/leipzig/freq/{info.leipzig_id}_freq.csv",
        ]
    return tracked, ignored


MUST_TRACK = [
    "docs/multilingual_design.md",
    "docs/runbooks.md",
    "docs/sources_lock.tsv",
    "docs/data_dict.md",
    "DATA_LICENSE.md",
    "README.md",
    "romanian_features.json",
    "src/phonology/inventories/romanian_features.json",
    "analysis/comparison/tp_long.csv",
    "tests/fixtures/kaikki_it_sample.jsonl",
    "tests/fixtures/README.md",
]
MUST_IGNORE = [
    "data/sources/kaikki/raw-wiktextract-data.jsonl.gz",
    "data/research_v2/map.json",
    "data/romanian_lexicon_with_freq.csv",
    "docs/steriade_2008.pdf",
]


def _ignored(paths: list[str]) -> dict[str, bool]:
    """Pattern verdict of ``git check-ignore`` for each path."""
    if shutil.which("git") is None:
        pytest.skip("git is not installed")
    result = subprocess.run(
        [
            "git",
            "-C",
            str(PROJECT_ROOT),
            "check-ignore",
            "--no-index",
            "--verbose",
            "--non-matching",
            "--stdin",
        ],
        input="\n".join(paths) + "\n",
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    if result.returncode not in (0, 1):
        pytest.skip(f"git check-ignore failed: {result.stderr.strip()}")
    verdict: dict[str, bool] = {}
    for line in result.stdout.splitlines():
        source, _, path = line.partition("\t")
        pattern = source.split(":", 2)[2] if source.count(":") >= 2 else ""
        verdict[path] = bool(pattern) and not pattern.startswith("!")
    return verdict


def test_tracked_paths_are_not_ignored() -> None:
    tracked, ignored = _language_paths()
    must_track = MUST_TRACK + tracked
    must_ignore = MUST_IGNORE + ignored
    verdict = _ignored(must_track + must_ignore)
    assert set(verdict) == set(must_track + must_ignore)
    assert [p for p in must_track if verdict[p]] == []
    assert [p for p in must_ignore if not verdict[p]] == []
