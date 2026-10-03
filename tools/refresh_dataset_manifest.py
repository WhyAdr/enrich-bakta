#!/usr/bin/env python3
"""Refresh hashes and sizes in an existing dataset manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def refresh(
    input_path: Path,
    output_path: Path,
    *,
    root: Path,
    producer_version: str | None,
    code_commit: str | None,
    excluded_paths: set[str],
) -> None:
    payload: dict[str, Any] = json.loads(input_path.read_text(encoding="utf-8-sig"))
    artifacts = payload.get("artifacts")
    if not isinstance(artifacts, list):
        raise ValueError("dataset manifest artifacts must be a list")
    for artifact in artifacts:
        if not isinstance(artifact, dict):
            raise ValueError("dataset manifest artifact must be an object")
        relative = artifact.get("path")
        if not isinstance(relative, str) or "\\" in relative:
            raise ValueError(f"invalid artifact path: {relative!r}")
        path = root.joinpath(*relative.split("/"))
        if not path.is_file():
            raise FileNotFoundError(path)
        artifact["size_bytes"] = path.stat().st_size
        artifact["sha256"] = sha256(path)
        producer = artifact.get("producer")
        if (
            producer_version
            and relative not in excluded_paths
            and isinstance(producer, dict)
        ):
            if producer.get("name") == "enrich-bakta":
                producer["version"] = producer_version
        validation = artifact.get("validation")
        if code_commit and isinstance(validation, dict):
            validation["code_commit"] = code_commit
    output_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    parser.add_argument("--producer-version")
    parser.add_argument("--code-commit")
    parser.add_argument("--exclude-path", action="append", default=[])
    args = parser.parse_args()
    refresh(
        args.input,
        args.output,
        root=args.root.resolve(),
        producer_version=args.producer_version,
        code_commit=args.code_commit,
        excluded_paths=set(args.exclude_path),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
