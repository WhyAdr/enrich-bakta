#!/usr/bin/env python3
"""Render the canonical dataset manifest as a deterministic TSV projection."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

FIELDS = (
    "id",
    "path",
    "sample",
    "role",
    "size_bytes",
    "sha256",
    "lineage",
    "producer_name",
    "producer_version",
    "producer_database",
    "status",
    "redistribution_license",
    "redistribution_permission",
    "validation_code_commit",
    "validation_checks",
)


def _cell(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, (list, dict)):
        return json.dumps(
            value, ensure_ascii=False, separators=(",", ":"), sort_keys=True
        )
    return str(value)


def render(manifest_path: Path, output_path: Path) -> None:
    payload = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    rows = []
    for artifact in payload["artifacts"]:
        producer = artifact["producer"]
        redistribution = artifact["redistribution"]
        validation = artifact["validation"]
        rows.append(
            {
                "id": artifact["id"],
                "path": artifact["path"],
                "sample": artifact["sample"],
                "role": artifact["role"],
                "size_bytes": artifact["size_bytes"],
                "sha256": artifact["sha256"],
                "lineage": artifact["lineage"],
                "producer_name": producer["name"],
                "producer_version": producer["version"],
                "producer_database": producer["database"],
                "status": artifact["status"],
                "redistribution_license": redistribution.get("license"),
                "redistribution_permission": redistribution.get("permission"),
                "validation_code_commit": validation.get("code_commit"),
                "validation_checks": validation.get("checks", []),
            }
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=FIELDS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows({key: _cell(row[key]) for key in FIELDS} for row in rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    render(args.manifest, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
