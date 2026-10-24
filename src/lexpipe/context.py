"""Run context of a language pipeline: paths, seed, forcing, manifests.

Every stage is a function of one :class:`RunContext` (design rule R6):
the language's registry row, the root directory, the random seed, the
set of stages to re-run unconditionally (``--force``) and free-form
run parameters. Stages keep no module-level state; everything they
need is on the context, so ``--root DIR`` runs a whole pipeline in a
temporary tree.

Paths
-----
All directories derive from the root and the registry row::

    data_dir        data/<name>/              generated, gitignored
    nc_dir          data/<name>/nc/           NC-licensed derivatives (R9)
    manifest_dir    data/<name>/manifests/    stage manifests (below)
    analysis_dir    analysis/<name>/          tracked summaries
    comparison_dir  analysis/comparison/      cross-language tables
    curated_dir     curated/<name>/           hand decisions, read-only
    sources_dir     data/sources/             pinned downloads, read-only
    docs_dir        docs/<name>/              data dictionaries, notes
    lock_path       docs/sources_lock.tsv     sha256 lock of the sources

Writing
-------
Outputs are written atomically: to a hidden ``.<name>.<pid>.tmp`` file
in the target directory, then renamed over the target with
:func:`os.replace`, so a reader never sees half a file
(:func:`atomic_open`, :func:`atomic_write_text`,
:func:`atomic_write_bytes`). The context's own writers
(:meth:`RunContext.open_output`, :meth:`RunContext.write_text`,
:meth:`RunContext.write_bytes`) first call
:meth:`RunContext.check_writable`, which refuses the locations a
derived stage must never touch: ``data/sources/`` (managed by the
``fetch`` stage), ``curated/`` (hand decisions), ``data/leipzig/``
(Romanian; R2) and files directly in ``data/`` or ``analysis/``
(the Romanian pipeline's inputs and tracked outputs; R1).

A comparison-only language (``LanguageInfo.comparison_only``, see
:mod:`lexpipe.spec`) is further confined to ``analysis/comparison/``
and its own ``data/<name>/`` (:func:`comparison_write_dirs`), and so
is every code listed in :data:`WRITE_ONLY_UNDER`, which is derived
from the registry's comparison-only rows. Keying the second rule on
the code keeps a context confined even when it was built from a
hand-made row without the flag (:func:`write_confinements`). The
harmonized Romanian comparison (code ``ro``, design §7) reads the
Romanian pipeline's files and may write only under
``analysis/comparison/`` and ``data/romanian_compare/``.

Tracked text reports
--------------------
:meth:`RunContext.write_tracked_text` is the licence gate of R9 for
text (CSV tables use :func:`lexpipe.csvio.write_tracked`): the path
must be a tracked location (:func:`lexpipe.licensing.is_tracked_location`)
and the caller names the report's sources. By default every source
must be redistributable. With ``nc_aggregates=True`` the report may
rest on CC BY-NC-SA sources too, under the draft policy of
``DATA_LICENSE.md``: every source must have a known licence (never
``unverified``), and the standard header of
:func:`lexpipe.licensing.aggregate_header` (sources, licences,
"aggregate counts only") is written first.

Manifests
---------
After a stage succeeds, :func:`write_manifest` records
``<manifest_dir>/<stage>.manifest.tsv`` (TSV, not JSON: ``*.json`` is
gitignored and TSV diffs well)::

    kind    key             value
    stage   name            <stage>
    git     sha             <HEAD sha, empty without git>
    git     dirty           True | False (tracked files only), empty without git
    code    sha256          <code_fingerprint of the language>
    param   <name>          <canonical JSON of the value>
    input   <path>          <sha256>
    output  <path>          <sha256>

Paths are POSIX paths relative to the root (absolute when outside it);
rows of a kind are sorted by key. :func:`should_skip` is true when the
stage is not forced, the manifest exists, its code fingerprint, its
parameters and its input set equal the current ones, every input still
has its recorded SHA-256, and every recorded output exists with its
recorded SHA-256. The code fingerprint (:func:`code_fingerprint`) is a
SHA-256 over the shared core and the language's plugin source, so
editing either re-runs the stage instead of silently reusing outputs
built by older code; a manifest without it (written before it existed)
counts as stale. The git state is provenance only and never compared.
A stage that uses the seed must put it in its parameters.
:func:`run_stage` combines both steps around a stage body.
"""

from __future__ import annotations

import contextlib
import csv
import hashlib
import importlib.util
import json
import os
import subprocess
from collections.abc import Callable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import IO, Any, Final

from lexpipe import PROJECT_ROOT, licensing
from lexpipe.spec import LANGUAGES, LanguageInfo

# The seed of the Romanian permutation test (lexpipe.stats), reused as
# the default so that Romanian and new-language runs agree by default.
DEFAULT_SEED: Final[int] = 20080715

MANIFEST_SUFFIX: Final[str] = ".manifest.tsv"
MANIFEST_FIELDS: Final[tuple[str, str, str]] = ("kind", "key", "value")
_KINDS: Final[tuple[str, ...]] = (
    "stage",
    "git",
    "code",
    "param",
    "input",
    "output",
)
# Plain TSV: values (canonical JSON, paths, hashes) are written
# verbatim; one holding a tab or newline raises csv.Error.
_TSV: Final[Mapping[str, Any]] = MappingProxyType(
    {"delimiter": "\t", "quoting": csv.QUOTE_NONE, "quotechar": None}
)
_GIT_TIMEOUT_S: Final[float] = 30.0
_CHUNK: Final[int] = 1 << 20

COMPARISON_DIR: Final[str] = "analysis/comparison"


def comparison_write_dirs(info: LanguageInfo) -> tuple[str, ...]:
    """The directories (POSIX, relative to the root) a comparison-only
    language may write under: ``analysis/comparison`` and
    ``data/<name>``."""
    return (COMPARISON_DIR, f"data/{info.name}")


# Codes confined to some directories (POSIX paths relative to the
# root), derived from the registry's comparison-only rows; every other
# guard still applies. "ro" is the read-only harmonized Romanian
# comparison (lexpipe.langs.romanian_compare).
WRITE_ONLY_UNDER: Final[Mapping[str, tuple[str, ...]]] = MappingProxyType(
    {
        code: comparison_write_dirs(info)
        for code, info in LANGUAGES.items()
        if info.comparison_only
    }
)


def write_confinements(info: LanguageInfo) -> tuple[tuple[str, ...], ...]:
    """Every directory list a context of *info* is confined to.

    A path must lie under one directory of **each** list: the code's
    :data:`WRITE_ONLY_UNDER` entry, and :func:`comparison_write_dirs`
    when the row itself is comparison-only. Empty for a full language.
    """
    found: list[tuple[str, ...]] = []
    listed = WRITE_ONLY_UNDER.get(info.code)
    if listed is not None:
        found.append(listed)
    if info.comparison_only:
        own = comparison_write_dirs(info)
        if own not in found:
            found.append(own)
    return tuple(found)


# ---------------------------------------------------------------------------
# Atomic writers
# ---------------------------------------------------------------------------


def _tmp_path(path: Path) -> Path:
    return path.with_name(f".{path.name}.{os.getpid()}.tmp")


@contextlib.contextmanager
def atomic_open(
    path: Path,
    mode: str = "w",
    *,
    encoding: str = "utf-8",
    newline: str | None = None,
) -> Iterator[IO[Any]]:
    """Open a temporary file next to *path*; rename it over *path* on exit.

    *mode* is ``"w"`` (text) or ``"wb"`` (binary). If the ``with``
    block raises, the temporary file is removed and *path* is left as
    it was.
    """
    if mode not in ("w", "wb"):
        raise ValueError(f"atomic_open writes only: mode {mode!r}")
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = _tmp_path(path)
    try:
        if mode == "wb":
            with tmp.open("wb") as fh:
                yield fh
        else:
            with tmp.open("w", encoding=encoding, newline=newline) as fh:
                yield fh
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def atomic_write_text(path: Path, text: str, *, newline: str = "") -> Path:
    """Write *text* (UTF-8, line endings untranslated) atomically."""
    with atomic_open(path, "w", newline=newline) as fh:
        fh.write(text)
    return path


def atomic_write_bytes(path: Path, data: bytes) -> Path:
    """Write *data* atomically."""
    with atomic_open(path, "wb") as fh:
        fh.write(data)
    return path


def sha256_file(path: Path) -> str:
    """Hex SHA-256 of a file, read in chunks."""
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(_CHUNK):
            digest.update(chunk)
    return digest.hexdigest()


# ---------------------------------------------------------------------------
# The context
# ---------------------------------------------------------------------------


def _frozen(params: Mapping[str, Any]) -> Mapping[str, Any]:
    return MappingProxyType(dict(params))


@dataclass(frozen=True, slots=True)
class RunContext:
    """Everything a stage needs (see the module docstring).

    ``root`` is resolved to an absolute path; ``force`` and ``params``
    are frozen copies, so a context can be shared between stages.
    """

    language: LanguageInfo
    root: Path = PROJECT_ROOT
    seed: int = DEFAULT_SEED
    force: frozenset[str] = frozenset()
    params: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "root", Path(self.root).resolve())
        object.__setattr__(self, "force", frozenset(self.force))
        object.__setattr__(self, "params", _frozen(self.params))
        if isinstance(self.seed, bool) or not isinstance(self.seed, int):
            raise TypeError(f"seed must be an int, got {self.seed!r}")

    # -- paths -------------------------------------------------------------

    @property
    def data_dir(self) -> Path:
        return self.root / "data" / self.language.name

    @property
    def nc_dir(self) -> Path:
        return self.data_dir / "nc"

    @property
    def manifest_dir(self) -> Path:
        return self.data_dir / "manifests"

    @property
    def analysis_dir(self) -> Path:
        return self.root / "analysis" / self.language.name

    @property
    def comparison_dir(self) -> Path:
        return self.root / "analysis" / "comparison"

    @property
    def curated_dir(self) -> Path:
        return self.root / "curated" / self.language.name

    @property
    def sources_dir(self) -> Path:
        return self.root / "data" / "sources"

    @property
    def docs_dir(self) -> Path:
        return self.root / "docs" / self.language.name

    @property
    def lock_path(self) -> Path:
        return self.root / "docs" / "sources_lock.tsv"

    def manifest_path(self, stage: str) -> Path:
        return self.manifest_dir / f"{stage}{MANIFEST_SUFFIX}"

    def forced(self, stage: str) -> bool:
        """True if *stage* was named with ``--force``."""
        return stage in self.force

    def relative(self, path: Path) -> str:
        """POSIX path relative to the root, or absolute if outside it."""
        resolved = Path(path).resolve()
        try:
            return resolved.relative_to(self.root).as_posix()
        except ValueError:
            return resolved.as_posix()

    # -- guarded writers ---------------------------------------------------

    def check_writable(self, path: Path) -> Path:
        """Return *path* resolved, or raise :class:`PermissionError` if a
        derived stage must not write there (see the module docstring)."""
        resolved = Path(path).resolve()
        refused = (
            (self.sources_dir, "data/sources/ is managed by the fetch stage"),
            (self.root / "curated", "curated/ holds hand decisions"),
            (self.root / "data" / "leipzig", "data/leipzig/ is Romanian (R2)"),
        )
        for directory, reason in refused:
            if resolved == directory or directory in resolved.parents:
                raise PermissionError(f"refusing to write {path}: {reason}")
        if resolved.parent in (self.root / "data", self.root / "analysis"):
            raise PermissionError(
                f"refusing to write {path}: files directly in data/ and "
                "analysis/ belong to the Romanian pipeline (R1)"
            )
        for confined in write_confinements(self.language):
            if not any(
                (self.root / directory) in resolved.parents
                for directory in confined
            ):
                raise PermissionError(
                    f"refusing to write {path}: language "
                    f"{self.language.code!r} may write only under "
                    f"{', '.join(confined)}"
                )
        return resolved

    @contextlib.contextmanager
    def open_output(
        self, path: Path, mode: str = "w", *, newline: str | None = None
    ) -> Iterator[IO[Any]]:
        """:func:`atomic_open` after :meth:`check_writable`."""
        target = self.check_writable(path)
        with atomic_open(target, mode, newline=newline) as fh:
            yield fh

    def write_text(self, path: Path, text: str) -> Path:
        return atomic_write_text(self.check_writable(path), text)

    def write_bytes(self, path: Path, data: bytes) -> Path:
        return atomic_write_bytes(self.check_writable(path), data)

    def write_tracked_text(
        self,
        path: Path,
        text: str,
        *,
        sources: str | Sequence[str],
        nc_aggregates: bool = False,
    ) -> Path:
        """Write a tracked text report behind the licence gate (R9).

        *sources* are the source codes the report rests on (a
        sequence, or one ``|``-joined string). Raises
        :class:`PermissionError` where :meth:`check_writable` does, and
        :class:`lexpipe.licensing.LicenceError` when *path* is not a
        tracked location or a source is refused: without
        *nc_aggregates*, any source that is not redistributable; with
        it, any source whose terms are unverified or unknown (see the
        module docstring). With *nc_aggregates* the text is preceded
        by :func:`lexpipe.licensing.aggregate_header`. Nothing is
        written when a check fails.
        """
        target = self.check_writable(path)
        licensing.check_tracked_location(
            target, self.root, licensing.with_language(self.language)
        )
        names = licensing.split_sources(sources)
        if nc_aggregates:
            text = licensing.aggregate_header(names) + text
        else:
            licensing.check_redistributable(names)
        return atomic_write_text(target, text)


# ---------------------------------------------------------------------------
# Git provenance (read-only commands; absence of git is tolerated)
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class GitState:
    """HEAD and whether tracked files differ from it (None = unknown)."""

    sha: str | None
    dirty: bool | None


def _git(root: Path, *args: str) -> str | None:
    try:
        result = subprocess.run(
            ["git", "--no-optional-locks", "-C", str(root), *args],
            capture_output=True,
            text=True,
            timeout=_GIT_TIMEOUT_S,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return result.stdout if result.returncode == 0 else None


def git_state(root: Path) -> GitState:
    """``git rev-parse HEAD`` and ``git status --porcelain`` (tracked
    files only) of *root*; both None outside a repository."""
    sha = _git(root, "rev-parse", "--verify", "HEAD")
    if sha is None:
        return GitState(None, None)
    status = _git(root, "status", "--porcelain", "--untracked-files=no")
    dirty = None if status is None else bool(status.strip())
    return GitState(sha.strip() or None, dirty)


# ---------------------------------------------------------------------------
# Manifests
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Manifest:
    """A parsed stage manifest; mappings are keyed as in the file."""

    stage: str
    params: Mapping[str, str]
    inputs: Mapping[str, str]
    outputs: Mapping[str, str]
    git: GitState
    code: str | None = None


def _code_files(language: LanguageInfo) -> list[Path]:
    """The Python source a stage of *language* runs: core + plugin."""
    core_dir = Path(__file__).resolve().parent
    files = set(core_dir.glob("*.py"))
    try:
        spec = (
            importlib.util.find_spec(language.module)
            if language.module
            else None
        )
    except (ImportError, ValueError):
        spec = None
    if spec is not None:
        if spec.submodule_search_locations:
            for location in spec.submodule_search_locations:
                files.update(Path(location).rglob("*.py"))
        elif spec.origin and spec.origin.endswith(".py"):
            files.add(Path(spec.origin))
    return sorted(files)


def code_fingerprint(language: LanguageInfo) -> str:
    """SHA-256 over the shared core and the plugin source of *language*.

    Covers every ``src/lexpipe/*.py`` module and the language's plugin
    (its package directory, recursively, or its single module; nothing
    if the plugin cannot be found). Paths enter the digest relative to
    ``src/`` so the value does not depend on where the checkout lives.
    """
    src_dir = Path(__file__).resolve().parent.parent
    digest = hashlib.sha256()
    for path in _code_files(language):
        try:
            name = path.resolve().relative_to(src_dir).as_posix()
        except ValueError:
            name = path.name
        digest.update(name.encode("utf-8") + b"\0")
        digest.update(path.read_bytes() + b"\0")
    return digest.hexdigest()


def canonical_param(value: Any) -> str:
    """Canonical JSON of a parameter value (sorted keys, no spaces).

    Sets and other non-JSON values raise :class:`TypeError`; NaN and
    infinities raise :class:`ValueError`.
    """
    return json.dumps(
        value,
        sort_keys=True,
        ensure_ascii=False,
        separators=(",", ":"),
        allow_nan=False,
    )


def _canonical_params(params: Mapping[str, Any]) -> dict[str, str]:
    out: dict[str, str] = {}
    for name, value in params.items():
        if not name or any(ch in name for ch in "\t\r\n"):
            raise ValueError(f"invalid parameter name: {name!r}")
        out[name] = canonical_param(value)
    return out


def _digests(ctx: RunContext, paths: Sequence[Path]) -> dict[str, str]:
    out: dict[str, str] = {}
    for path in paths:
        if not Path(path).is_file():
            raise FileNotFoundError(f"stage file missing: {path}")
        out[ctx.relative(path)] = sha256_file(Path(path))
    return out


def _bool_cell(value: bool | None) -> str:
    return "" if value is None else str(value)


def _parse_bool_cell(text: str) -> bool | None:
    lowered = text.strip().lower()
    if not lowered:
        return None
    if lowered in ("true", "false"):
        return lowered == "true"
    raise ValueError(f"not a boolean: {text!r}")


def write_manifest(
    ctx: RunContext,
    stage: str,
    *,
    inputs: Sequence[Path],
    params: Mapping[str, Any],
    outputs: Sequence[Path],
    input_digests: Mapping[str, str] | None = None,
    code: str | None = None,
) -> Path:
    """Record a successful run of *stage*; return the manifest path.

    *input_digests* (as returned by an earlier hashing of *inputs*)
    saves re-hashing large inputs; *code* likewise passes a
    :func:`code_fingerprint` computed earlier in the same run.
    """
    rows: list[tuple[str, str, str]] = [("stage", "name", stage)]
    git = git_state(ctx.root)
    rows.append(("git", "dirty", _bool_cell(git.dirty)))
    rows.append(("git", "sha", git.sha or ""))
    fingerprint = code if code is not None else code_fingerprint(ctx.language)
    rows.append(("code", "sha256", fingerprint))
    canonical = _canonical_params(params)
    rows.extend(("param", k, canonical[k]) for k in sorted(canonical))
    digests = (
        dict(input_digests)
        if input_digests is not None
        else _digests(ctx, inputs)
    )
    rows.extend(("input", k, digests[k]) for k in sorted(digests))
    written = _digests(ctx, outputs)
    rows.extend(("output", k, written[k]) for k in sorted(written))
    path = ctx.manifest_path(stage)
    with atomic_open(path, "w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n", **_TSV)
        writer.writerow(MANIFEST_FIELDS)
        writer.writerows(rows)
    return path


def read_manifest(path: Path) -> Manifest | None:
    """Parse a manifest; None if it is absent or malformed."""
    if not path.is_file():
        return None
    sections: dict[str, dict[str, str]] = {kind: {} for kind in _KINDS}
    try:
        with path.open("r", encoding="utf-8", newline="") as fh:
            reader = csv.reader(fh, **_TSV)
            if tuple(next(reader, ())) != MANIFEST_FIELDS:
                return None
            for row in reader:
                if len(row) != 3 or row[0] not in sections:
                    return None
                sections[row[0]][row[1]] = row[2]
        git = GitState(
            sections["git"].get("sha") or None,
            _parse_bool_cell(sections["git"].get("dirty", "")),
        )
    except (OSError, csv.Error, ValueError):
        return None
    stage = sections["stage"].get("name")
    if not stage:
        return None
    return Manifest(
        stage=stage,
        params=MappingProxyType(sections["param"]),
        inputs=MappingProxyType(sections["input"]),
        outputs=MappingProxyType(sections["output"]),
        git=git,
        code=sections["code"].get("sha256") or None,
    )


def _outputs_intact(ctx: RunContext, outputs: Mapping[str, str]) -> bool:
    for key, digest in outputs.items():
        path = Path(key) if Path(key).is_absolute() else ctx.root / key
        if not path.is_file() or sha256_file(path) != digest:
            return False
    return True


def _up_to_date(
    ctx: RunContext,
    stage: str,
    digests: Mapping[str, str],
    params: Mapping[str, Any],
    code: str,
) -> bool:
    if ctx.forced(stage):
        return False
    manifest = read_manifest(ctx.manifest_path(stage))
    if manifest is None or manifest.stage != stage:
        return False
    if manifest.code != code:
        return False
    if dict(manifest.params) != _canonical_params(params):
        return False
    if dict(manifest.inputs) != dict(digests):
        return False
    return _outputs_intact(ctx, manifest.outputs)


def should_skip(
    ctx: RunContext,
    stage: str,
    inputs: Sequence[Path],
    params: Mapping[str, Any],
) -> bool:
    """True if *stage* need not run (see the module docstring).

    A missing input raises :class:`FileNotFoundError`: the stage could
    not run either.
    """
    return _up_to_date(
        ctx,
        stage,
        _digests(ctx, inputs),
        params,
        code_fingerprint(ctx.language),
    )


def _failed(result: object) -> bool:
    return (
        isinstance(result, int)
        and not isinstance(result, bool)
        and result != 0
    )


def run_stage(
    ctx: RunContext,
    stage: str,
    *,
    inputs: Sequence[Path],
    params: Mapping[str, Any],
    outputs: Sequence[Path],
    body: Callable[[], int | None],
    log: Callable[[str], None] = print,
) -> int | None:
    """Run *body* unless :func:`should_skip`; record the manifest after.

    Inputs are hashed once. A body that returns a non-zero int is a
    failure: its status is returned and no manifest is written. Every
    path in *outputs* must exist afterwards.
    """
    digests = _digests(ctx, inputs)
    code = code_fingerprint(ctx.language)
    if _up_to_date(ctx, stage, digests, params, code):
        log(f"  {stage}: up to date (manifest matches); skipped")
        return 0
    result = body()
    if _failed(result):
        return result
    missing = [str(p) for p in outputs if not Path(p).is_file()]
    if missing:
        raise FileNotFoundError(
            f"stage {stage!r} did not write its outputs: {missing}"
        )
    write_manifest(
        ctx,
        stage,
        inputs=inputs,
        params=params,
        outputs=outputs,
        input_digests=digests,
        code=code,
    )
    return result
