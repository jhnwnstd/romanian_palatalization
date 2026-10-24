"""Run every Italian stage on the analysis fixture tree, in a subprocess.

Usage: ``python analysis_driver.py ROOT [SHUFFLE_SEED]``. Builds the
analysis tree under ROOT (shuffling every input when a seed is given)
and runs the whole pipeline, data layer and analyses. The determinism
tests call it with different ``PYTHONHASHSEED`` values (R5).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from tests.italian import analysis_helpers, helpers  # noqa: E402


def main(argv: list[str]) -> int:
    root = Path(argv[1])
    seed = int(argv[2]) if len(argv) > 2 else None
    analysis_helpers.build_analysis_tree(root, seed)
    helpers.run_pipeline(analysis_helpers.make_context(root))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
