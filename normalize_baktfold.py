"""Source-only compatibility launcher for the historical normalizer."""

from pathlib import Path
import runpy

if __name__ == "__main__":
    runpy.run_path(
        str(Path(__file__).resolve().parent / "legacy/normalize_baktfold.py"),
        run_name="__main__",
    )
