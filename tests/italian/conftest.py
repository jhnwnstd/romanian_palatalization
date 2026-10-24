"""pytest configuration of the Italian tests.

Puts ``src/`` on ``sys.path`` before the test modules are imported, so
they can import :mod:`lexpipe` directly (the other test packages insert
the path by hand at the top of each module).
"""

from __future__ import annotations

import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[2] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
