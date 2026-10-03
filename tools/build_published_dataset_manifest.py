#!/usr/bin/env python3
"""Build the reviewed C14/SM/BK71A publication manifest from the local inventory."""

from __future__ import annotations

import argparse
import copy
import fnmatch
import json
from pathlib import Path
from typing import Any

SCHEMA = "enrich-bakta.dataset-manifest.v1"
CURATED_PATTERNS = (
    "data/C14/bakta/*.gbff",
    "data/C14/bakta/*.faa",
    "data/C14/baktfold/*.gbff",
    "data/C14/baktfold/*.faa",
    "data/C14/evidence/*-KofamKOALA.txt",
    "data/C14/evidence/*-emapper_annotations.xlsx",
    "data/C14/evidence/*-query.emapper.*",
    "data/C14/derived/restored/*.gbff",
    "data/C14/derived/restored/*.manifest.json",
    "data/C14/derived/enriched/*.gbff",
    "data/C14/derived/enriched/*.manifest.json",
    "data/SM/bakta/*.gbff",
    "data/SM/bakta/*.faa",
    "data/SM/baktfold/*.gbff",
    "data/SM/baktfold/*.faa",
    "data/SM/evidence/*-KofamKOALA.txt",
    "data/SM/evidence/*-emapper_annotations.xlsx",
    "data/SM/evidence/*-query.emapper.*",
    "data/SM/derived/restored/*.gbff",
    "data/SM/derived/restored/*.manifest.json",
    "data/SM/derived/enriched/*.gbff",
    "data/SM/derived/enriched/*.manifest.json",
    "data/BK71A/derived/restored/BK71A-restored.gbff",
)


def is_curated(path: str) -> bool:
    return any(fnmatch.fnmatchcase(path, pattern) for pattern in CURATED_PATTERNS)


def build(input_path: Path, output_path: Path) -> None:
    payload = json.loads(input_path.read_text(encoding="utf-8-sig"))
    if payload.get("schema") != SCHEMA:
        raise ValueError(f"unsupported dataset manifest schema: {input_path}")
    source_artifacts = payload.get("artifacts")
    if not isinstance(source_artifacts, list):
        raise ValueError("dataset manifest artifacts must be a list")

    selected = [
        copy.deepcopy(artifact)
        for artifact in source_artifacts
        if isinstance(artifact, dict) and is_curated(str(artifact.get("path", "")))
    ]
    selected.sort(key=lambda artifact: str(artifact["path"]))
    selected_ids = {str(artifact["id"]) for artifact in selected}
    for artifact in selected:
        artifact["lineage"] = [
            reference
            for reference in artifact.get("lineage", [])
            if reference in selected_ids
        ]
        artifact["status"] = "historical_unverified"
        artifact["redistribution"] = {
            "license": "CC-BY-4.0",
            "permission": "owner-authorized",
        }
        validation = artifact.setdefault("validation", {})
        checks = list(validation.get("checks", []))
        if "owner-authorized-publication" not in checks:
            checks.append("owner-authorized-publication")
        validation["checks"] = sorted(str(check) for check in checks)

    result: dict[str, Any] = {
        "schema": SCHEMA,
        "publication": {
            "license": "CC-BY-4.0",
            "storage": "git-lfs",
            "status_policy": "historical_unverified until producer metadata and lineage are independently confirmed",
            "owner_authorization": "2026-10-03",
        },
        "artifacts": selected,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(result, indent=2, ensure_ascii=False, sort_keys=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    build(args.input, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
