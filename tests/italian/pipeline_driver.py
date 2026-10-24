"""Run the Italian data-layer stages on a fixture tree, in a subprocess.

Usage: ``python pipeline_driver.py ROOT [SHUFFLE_SEED]``. Builds the
fixture tree under ROOT (shuffling every input when a seed is given)
and runs extract, lexicon, qc, attest and features. The determinism
tests call it with different ``PYTHONHASHSEED`` values (R5).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tests.italian.helpers import (  # noqa: E402
    build_tree,
    make_context,
    run_pipeline,
)


def main(argv: list[str]) -> int:
    root = Path(argv[1])
    seed = int(argv[2]) if len(argv) > 2 else None
    build_tree(root, seed)
    run_pipeline(make_context(root))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
