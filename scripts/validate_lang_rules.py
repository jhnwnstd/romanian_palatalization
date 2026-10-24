#!/usr/bin/env python3
"""Run a non-Romanian palatalization analysis through the validator.

The Italian and Czech counterpart of
``scripts/validate_palatalization_rules.py`` (which stays the Romanian
driver and is not touched). Usage::

    python scripts/validate_lang_rules.py --lang it|cs [--variant NAME]

Each analysis module ``phonology.analyses.<language>_palatalization``
must expose:

  - ``build_profile() -> AnalysisProfile`` — the default analysis as
    data, with its own hooks (``ur_builder``, ``row_factory``,
    ``tokenizer``, ``normalisers``, ...);
  - ``iter_rows(root: Path) -> Iterable[LexRow]`` — the rows to
    validate, read from the language's generated data under
    ``root/data/<language>/``;

and may expose ``PROFILE_VARIANTS: Mapping[str, Callable[[],
AnalysisProfile]]`` for alternative analyses of the same rows (for
example lexically specified vs. uniformly underspecified stems).

Outputs go to ``analysis/<language>/<language>_rule_validation.txt``
and ``<language>_rule_mismatches.csv``; a variant adds ``_<variant>``
before the extension. Both carry only Wiktionary-derived forms and IPA
(CC BY-SA), so they may be tracked: the driver checks that licence
(design rule R9) and writes through the language's
:class:`~lexpipe.context.RunContext`, whose write guard keeps it inside
``analysis/<language>/``.
"""

from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PROJECT_ROOT / "src"))

from lexpipe import licensing  # noqa: E402
from lexpipe.context import RunContext  # noqa: E402
from lexpipe.spec import language_info  # noqa: E402
from phonology.validation import (  # noqa: E402
    render_report,
    run_validation,
    write_mismatch_csv,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--lang", required=True, help="language code (it, cs)")
    parser.add_argument(
        "--variant",
        default=None,
        help="name in the module's PROFILE_VARIANTS (default: build_profile)",
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=_PROJECT_ROOT,
        help="project root holding data/ and analysis/ (default: this repo)",
    )
    args = parser.parse_args(argv)

    info = language_info(args.lang)
    module = importlib.import_module(
        f"phonology.analyses.{info.name}_palatalization"
    )
    if args.variant is None:
        profile = module.build_profile()
        suffix = ""
    else:
        variants = getattr(module, "PROFILE_VARIANTS", {})
        if args.variant not in variants:
            known = ", ".join(sorted(variants)) or "none"
            parser.error(f"unknown variant {args.variant!r} (known: {known})")
        profile = variants[args.variant]()
        suffix = f"_{args.variant}"
    # The rows carry forms and IPA from both Wiktionary editions only.
    licensing.check_redistributable(
        ("wiktionary_en", f"wiktionary_{info.kaikki_code}")
    )
    ctx = RunContext(info, root=args.root)
    report = run_validation(profile, module.iter_rows(args.root))

    text = render_report(report, include_distance_summary=True)
    report_path = ctx.analysis_dir / f"{info.name}_rule_validation{suffix}.txt"
    ctx.write_text(report_path, text)
    print(text, end="")
    print(f"\nReport persisted to {report_path}")

    mismatch_path = ctx.check_writable(
        ctx.analysis_dir / f"{info.name}_rule_mismatches{suffix}.csv"
    )
    mismatch_path.parent.mkdir(parents=True, exist_ok=True)
    n = write_mismatch_csv(report, mismatch_path)
    print(f"Mismatch CSV: {n:,} rows → {mismatch_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
