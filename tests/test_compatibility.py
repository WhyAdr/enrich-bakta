from __future__ import annotations

import importlib
import subprocess
import sys
from pathlib import Path

import pytest

PAIRS = {
    "merge_engine": "core.merge_engine",
    "graft_baktfold_additions": "sources.baktfold",
    "merge_kofamscan_bakta": "sources.kofam",
    "merge_eggnog_bakta": "sources.eggnog",
    "enrich_bakta": "workflows.enrich",
    "restore_bakta_translations": "workflows.restore_translations",
}


@pytest.mark.parametrize("legacy, canonical", PAIRS.items())
def test_alias_preserves_module_and_class_identity(legacy, canonical):
    old = importlib.import_module(legacy)
    new = importlib.import_module("enrich_bakta_lib." + canonical)
    assert old is new


@pytest.mark.parametrize("legacy", [name for name in PAIRS if name != "merge_engine"])
def test_direct_launcher_outside_checkout(tmp_path, legacy):
    root = Path(__file__).resolve().parents[1]
    result = subprocess.run(
        [sys.executable, str(root / (legacy + ".py")), "--help"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr


def test_module_entry_point(tmp_path):
    root = Path(__file__).resolve().parents[1]
    result = subprocess.run(
        [sys.executable, "-m", "enrich_bakta", "--version"],
        cwd=root,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "0.3.0" in result.stdout
