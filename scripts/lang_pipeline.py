#!/usr/bin/env python3
"""Run a language pipeline (Italian, Czech, the Romanian comparison).

Usage::

    python scripts/lang_pipeline.py --lang it              # all stages
    python scripts/lang_pipeline.py --lang cs --list       # list stages
    python scripts/lang_pipeline.py --lang it fetch verify-sources
    python scripts/lang_pipeline.py --lang it --force tp tp
    python scripts/lang_pipeline.py --lang cs --force qc-morfflex qc
    python scripts/lang_pipeline.py --lang cs --root /tmp/tree --seed 7
    python scripts/lang_pipeline.py --lang it --param learner_samples=200 tp
    python scripts/lang_pipeline.py --lang it --param 'learner_V=[500,1000]'
    python scripts/lang_pipeline.py --lang ro              # comparison only

Named stages run in pipeline order, whatever order they are given in.
Every stage is called with one :class:`lexpipe.context.RunContext`
built from the options (design rule R6):

  --root DIR      the tree to work in (default: this repository), so
                  tests can run a pipeline in a temporary directory;
  --seed N        the seed of every random step (default 20080715);
  --force STAGE   re-run STAGE even if its manifest says it is up to
                  date (repeatable). Only derived stages can be
                  forced: ``fetch`` and ``verify-sources`` manage
                  data/sources/, and no stage writes to data/sources/
                  or curated/ (RunContext.check_writable). A plugin's
                  declared sub-stage (a nested step with its own
                  manifest, e.g. Czech ``qc-morfflex`` inside ``qc``;
                  :func:`lexpipe.spec.plugin_substages`) can be
                  forced too, when its stage runs; it is never run on
                  its own;
  --param KEY=VALUE
                  a run parameter, ``RunContext.params[KEY]``
                  (repeatable, each KEY once). VALUE is typed:
                  ``true``/``false`` (any case) → bool, an integer →
                  int, a decimal number → float, a JSON list
                  (``[500,1000]``) → list; anything else stays a
                  string. Plugins document the parameters they read
                  (e.g. ``learner_samples``, ``learner_V``).

The first two stages are the same for every language and manage the
pinned external files registered in :mod:`lexpipe.sources`:

  fetch           download the language's sources that are absent;
                  a file that already has the registered size is left
                  alone, one with another size is reported, untouched.
                  A derived source (the English Wiktionary subset
                  split from the kaikki dump) is never downloaded: if
                  it is absent, fetch prints the command that makes
                  it and fails.
  verify-sources  check size and SHA-256 of every source; the first
                  verification records the hash in
                  docs/sources_lock.tsv and later ones must match it.
                  For sources served at a rolling URL (kaikki), one
                  HEAD request records the server's ETag and
                  Last-Modified, and a later change is reported.

The remaining stages come from the language plugin
(``src/lexpipe/langs/<language>/``, contract in :mod:`lexpipe.spec`).
It is imported only when one of its stages is needed, so ``fetch``
and ``verify-sources`` also work before the plugin exists.

Each stage reports its wall-clock time and a summary follows. A stage
that returns a non-zero int stops the run with exit status 1; one that
raises stops it with the traceback.

Romanian has its own pipeline (``run_pipeline.sh``); nothing here
touches it. ``--lang ro`` runs the read-only harmonized comparison
(:mod:`lexpipe.langs.romanian_compare`, a comparison-only registry
row), whose run context may write only under ``analysis/comparison/``
and ``data/romanian_compare/``.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
import time
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, Final

PROJECT_ROOT: Final[Path] = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from lexpipe import context, sources, spec  # noqa: E402
from lexpipe.context import RunContext  # noqa: E402

# Run-parameter names: they become manifest keys (no tabs or newlines).
_PARAM_KEY_RE: Final[re.Pattern[str]] = re.compile(r"[A-Za-z_][A-Za-z0-9_.-]*")
_INT_RE: Final[re.Pattern[str]] = re.compile(r"[+-]?[0-9]+")
# Decimal numbers only: "nan" and "inf" stay strings (manifests refuse
# non-finite floats).
_FLOAT_RE: Final[re.Pattern[str]] = re.compile(
    r"[+-]?(?:[0-9]+\.[0-9]*|\.[0-9]+|[0-9]+)(?:[eE][+-]?[0-9]+)?"
)


def parse_param_value(text: str) -> Any:
    """Typed value of ``--param KEY=VALUE`` (see the module docstring).

    A value starting with ``[`` must be a JSON list of scalars;
    :class:`ValueError` otherwise.
    """
    stripped = text.strip()
    lowered = stripped.lower()
    if lowered in ("true", "false"):
        return lowered == "true"
    if _INT_RE.fullmatch(stripped):
        return int(stripped)
    if _FLOAT_RE.fullmatch(stripped):
        return float(stripped)
    if stripped.startswith("["):
        try:
            value = json.loads(stripped)
        except json.JSONDecodeError as exc:
            raise ValueError(
                f"not a JSON list: {text!r} ({exc.msg})"
            ) from None
        if not isinstance(value, list) or not all(
            isinstance(item, (bool, int, float, str)) for item in value
        ):
            raise ValueError(f"not a JSON list of scalars: {text!r}")
        if any(
            isinstance(item, float) and not math.isfinite(item)
            for item in value
        ):
            raise ValueError(f"non-finite number in {text!r}")
        return value
    return text


def parse_params(items: Sequence[str]) -> dict[str, Any]:
    """``--param`` items → ``RunContext.params``; :class:`ValueError` on
    a malformed item or a repeated key."""
    params: dict[str, Any] = {}
    for item in items:
        key, sep, value = item.partition("=")
        key = key.strip()
        if not sep or not _PARAM_KEY_RE.fullmatch(key):
            raise ValueError(
                f"--param {item!r}: expected KEY=VALUE with KEY made of "
                "letters, digits, '_', '-' and '.'"
            )
        if key in params:
            raise ValueError(f"--param {key} is given twice")
        params[key] = parse_param_value(value)
    return params


def _registry(ctx: RunContext) -> dict[str, sources.Source]:
    return {s.name: s for s in sources.sources_for(ctx.language.code)}


def fetch_sources(
    ctx: RunContext, *, session: sources.HttpSession | None = None
) -> int:
    """Common stage: download the sources of the language that are absent.

    Carries on after a failure so one unreachable server does not hide
    the state of the others; returns 1 if any source failed.
    """
    registry = _registry(ctx)
    failed: list[str] = []
    for name in registry:
        try:
            sources.fetch(
                name, root=ctx.root, registry=registry, session=session
            )
        except (sources.DownloadError, OSError) as exc:
            # OSError covers requests' network errors and disk errors.
            print(f"  ERROR {name}: {exc}")
            failed.append(name)
    if failed:
        print(f"  {len(failed)} source(s) not fetched: {', '.join(failed)}")
    return 1 if failed else 0


def _remote(
    source: sources.Source, session: sources.HeadSession | None
) -> sources.RemoteMeta | None:
    """ETag/Last-Modified of a rolling source (one HEAD request)."""
    if not source.rolling:
        return None
    try:
        return sources.head(source.url, session=session)
    except (sources.DownloadError, OSError) as exc:
        print(f"  {source.name}: HEAD failed ({exc}); upstream not checked")
        return None


def verify_sources(
    ctx: RunContext, *, session: sources.HeadSession | None = None
) -> int:
    """Common stage: size and SHA-256 of each source of the language."""
    registry = _registry(ctx)
    results = [
        sources.verify(
            name,
            root=ctx.root,
            registry=registry,
            lock_path=ctx.lock_path,
            remote=_remote(source, session),
        )
        for name, source in registry.items()
    ]
    width = max((len(r.name) for r in results), default=0)
    for result in results:
        print(
            f"  {result.name:<{width}}  {result.status:<15}  {result.detail}"
        )
    failed = [r.name for r in results if not r.ok]
    if failed:
        print(f"  {len(failed)} source(s) not verified: {', '.join(failed)}")
    return 1 if failed else 0


def common_stages() -> dict[str, spec.Stage]:
    """The stages every language gets."""
    return {"fetch": fetch_sources, "verify-sources": verify_sources}


def _failed(result: object) -> bool:
    # Stages return None or 0 on success; any other value is ignored
    # unless it is a non-zero int (bools are not statuses).
    return (
        isinstance(result, int)
        and not isinstance(result, bool)
        and result != 0
    )


def run_stages(
    ctx: RunContext, names: Sequence[str], stages: Mapping[str, spec.Stage]
) -> int:
    """Run *names* in order with timing; 1 on the first failed stage."""
    code = ctx.language.code
    timings: list[tuple[str, float]] = []
    status = 0
    for name in names:
        forced = "  (forced)" if ctx.forced(name) else ""
        print(f"[{code}] {name}: start{forced}", flush=True)
        start = time.perf_counter()
        try:
            result = stages[name](ctx)
        except BaseException:
            elapsed = time.perf_counter() - start
            print(f"[{code}] {name}: raised after {elapsed:.1f}s", flush=True)
            raise
        elapsed = time.perf_counter() - start
        timings.append((name, elapsed))
        if _failed(result):
            print(
                f"[{code}] {name}: failed (status {result}) after {elapsed:.1f}s"
            )
            status = 1
            break
        print(f"[{code}] {name}: done in {elapsed:.1f}s", flush=True)
    width = max([len("total"), *(len(name) for name, _ in timings)])
    print(f"[{code}] stage timings:")
    for name, elapsed in timings:
        print(f"  {name:<{width}}  {elapsed:8.1f}s")
    print(f"  {'total':<{width}}  {sum(t for _, t in timings):8.1f}s")
    return status


def _print_stage_list(
    info: spec.LanguageInfo,
    stages: Mapping[str, spec.Stage],
    common: Mapping[str, spec.Stage],
    substages: Mapping[str, str] | None = None,
) -> None:
    print(f"Stages for {info.name} ({info.code}), in pipeline order:")
    nested = substages or {}
    for name in stages:
        print(f"  {name}{'  (common)' if name in common else ''}")
        for sub, parent in nested.items():
            if parent == name:
                print(f"    {sub}  (sub-stage of {name}; --force only)")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0] if __doc__ else None,
        epilog="Stages always run in pipeline order.",
    )
    parser.add_argument(
        "--lang",
        required=True,
        choices=sorted(spec.LANGUAGES),
        help="language code",
    )
    parser.add_argument(
        "--list", action="store_true", help="list the stages and exit"
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=PROJECT_ROOT,
        metavar="DIR",
        help="tree to work in (default: this repository)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=context.DEFAULT_SEED,
        metavar="N",
        help=f"seed of every random step (default {context.DEFAULT_SEED})",
    )
    parser.add_argument(
        "--force",
        action="append",
        default=[],
        metavar="STAGE",
        help=(
            "re-run a derived stage, or a declared sub-stage of a stage "
            "that runs, even if up to date (repeatable)"
        ),
    )
    parser.add_argument(
        "--param",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help=(
            "run parameter for RunContext.params (repeatable; VALUE "
            "typed: true/false, int, float, JSON list, else string)"
        ),
    )
    parser.add_argument(
        "stages",
        nargs="*",
        metavar="stage",
        help="stages to run (default: all)",
    )
    return parser


def _check_force(
    parser: argparse.ArgumentParser,
    forced: Sequence[str],
    selected: Sequence[str],
    stages: Mapping[str, spec.Stage],
    substages: Mapping[str, str] | None = None,
) -> None:
    nested = substages or {}
    unknown = [
        name for name in forced if name not in stages and name not in nested
    ]
    if unknown:
        parser.error(f"--force names unknown stage(s): {', '.join(unknown)}")
    idle = [name for name in forced if nested.get(name, name) not in selected]
    if idle:
        described = [
            f"{name} (sub-stage of {nested[name]})" if name in nested else name
            for name in idle
        ]
        parser.error(
            "--force names stage(s) that are not run: " + ", ".join(described)
        )


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not args.root.is_dir():
        parser.error(f"--root {args.root}: not a directory")
    try:
        params = parse_params(args.param)
    except ValueError as exc:
        parser.error(str(exc))
    info = spec.language_info(args.lang)
    common = common_stages()
    requested = list(dict.fromkeys(args.stages))
    forced = list(dict.fromkeys(args.force))
    never = [name for name in forced if name in spec.COMMON_STAGES]
    if never:
        parser.error(
            f"--force {', '.join(never)}: the common stages manage "
            "data/sources/ and are never forced"
        )
    stages: dict[str, spec.Stage] = dict(common)
    substages: Mapping[str, str] = {}
    needs_plugin = any(s not in common for s in [*requested, *forced])
    if args.list or not requested or needs_plugin:
        try:
            plugin = spec.load_plugin(info.code)
        except spec.PluginNotFoundError as exc:
            if args.list:
                _print_stage_list(info, stages, common)
                sys.stdout.flush()
            print(f"error: {exc}", file=sys.stderr)
            return 1
        stages.update(plugin.stages)
        substages = spec.plugin_substages(plugin)
    if args.list:
        _print_stage_list(info, stages, common, substages)
        return 0
    nested = [name for name in requested if name in substages]
    if nested:
        parser.error(
            "; ".join(
                f"{name} is a sub-stage of {substages[name]}: run "
                f"{substages[name]} (with --force {name} to re-run it)"
                for name in nested
            )
        )
    unknown = [name for name in requested if name not in stages]
    if unknown:
        parser.error(
            f"unknown stage(s): {', '.join(unknown)}; "
            f"available: {', '.join(stages)}"
        )
    selected = [name for name in stages if not requested or name in requested]
    _check_force(parser, forced, selected, stages, substages)
    ctx = RunContext(
        language=info,
        root=args.root,
        seed=args.seed,
        force=frozenset(forced),
        params=params,
    )
    return run_stages(ctx, selected, stages)


if __name__ == "__main__":
    raise SystemExit(main())
