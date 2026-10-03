from __future__ import annotations

import hashlib
import json
from pathlib import Path

from tools.validate_dataset_manifest import validate_manifest


def _manifest(root: Path, *, path: str = "fixtures/sample.txt") -> Path:
    target = root / Path(*path.split("/"))
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(b"fixture\n")
    digest = hashlib.sha256(target.read_bytes()).hexdigest()
    payload = {
        "schema": "enrich-bakta.dataset-manifest.v1",
        "artifacts": [
            {
                "id": "fixture-sample",
                "path": path,
                "sample": "fixture",
                "role": "fixture",
                "size_bytes": target.stat().st_size,
                "sha256": digest,
                "lineage": [],
                "producer": {
                    "name": "synthetic",
                    "version": "1",
                    "database": "none",
                },
                "status": "validated",
                "redistribution": {"license": "MIT", "permission": "local"},
                "validation": {"code_commit": None, "checks": ["fixture"]},
            }
        ],
    }
    manifest = root / "manifest.json"
    manifest.write_text(json.dumps(payload), encoding="utf-8")
    return manifest


def test_dataset_manifest_validates_materialized_fixture(tmp_path: Path) -> None:
    assert validate_manifest(_manifest(tmp_path), tmp_path) == 0


def test_dataset_manifest_rejects_path_traversal(tmp_path: Path) -> None:
    manifest = _manifest(tmp_path, path="../outside.txt")
    assert validate_manifest(manifest, tmp_path) == 2
