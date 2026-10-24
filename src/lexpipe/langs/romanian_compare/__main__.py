"""Shortcut for ``python scripts/lang_pipeline.py --lang ro``.

    PYTHONPATH=src python -m lexpipe.langs.romanian_compare            # all
    PYTHONPATH=src python -m lexpipe.langs.romanian_compare --list
    PYTHONPATH=src python -m lexpipe.langs.romanian_compare tp parity
    PYTHONPATH=src python -m lexpipe.langs.romanian_compare --force tp tp

``ro`` is a comparison-only registry row, so the common pipeline runs
it (``--root``, ``--seed``, ``--force``, ``--param``, ``--list`` and
stage names mean the same there). The two options of the earlier
stand-alone command are translated: ``--learner-samples N`` is
``--param learner_samples=N`` and ``--learner-V 500,1000`` is
``--param 'learner_V=[500,1000]'``.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from collections.abc import Sequence
from types import ModuleType

from lexpipe import PROJECT_ROOT

SCRIPT = PROJECT_ROOT / "scripts" / "lang_pipeline.py"


def _pipeline() -> ModuleType:
    spec = importlib.util.spec_from_file_location("lang_pipeline", SCRIPT)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def translate(argv: Sequence[str]) -> list[str]:
    """The ``scripts/lang_pipeline.py`` arguments for *argv*."""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--learner-samples", default=None)
    parser.add_argument("--learner-V", default=None)
    known, rest = parser.parse_known_args(list(argv))
    out = ["--lang", "ro"]
    if known.learner_samples is not None:
        out += ["--param", f"learner_samples={known.learner_samples}"]
    if known.learner_V is not None:
        values = ",".join(
            str(int(v)) for v in str(known.learner_V).split(",") if v
        )
        out += ["--param", f"learner_V=[{values}]"]
    return out + rest


def main(argv: Sequence[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else list(argv)
    status = _pipeline().main(translate(args))
    return int(status)


if __name__ == "__main__":
    sys.exit(main())
