"""Compatibility entry point for enrich_bakta_lib.workflows.restore_translations."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

_source_root = Path(__file__).resolve().parent / "src"
if (_source_root / "enrich_bakta_lib").is_dir():
    sys.path.insert(0, str(_source_root))

_implementation = importlib.import_module(
    "enrich_bakta_lib.workflows.restore_translations"
)

if __name__ == "__main__":
    raise SystemExit(_implementation.main())
else:
    sys.modules[__name__] = _implementation
