#!/usr/bin/env python3
"""Validate an enrich-bakta dataset manifest against materialized files."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

SCHEMA = "enrich-bakta.dataset-manifest.v1"
ROLES = {
    "upstream_input",
    "functional_evidence",
    "annotation_source",
    "derived_output",
    "fixture",
}
STATUSES = {"historical_unverified", "validated", "published", "local_only"}
HASH_LENGTH = 64


class ManifestError(ValueError):
    """Raised for a manifest or materialized-file contract failure."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _safe_relative_path(value: Any) -> str:
    if not isinstance(value, str) or not value or "\\" in value:
        raise ManifestError(f"invalid relative path: {value!r}")
    path = Path(value)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise ManifestError(f"unsafe dataset path: {value!r}")
    normalized = "/".join(path.parts)
    if normalized != value:
        raise ManifestError(f"dataset path is not normalized: {value!r}")
    return normalized


def _require_mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ManifestError(f"{label} must be an object")
    return value


def _validate_artifact(
    artifact: Any,
    *,
    root: Path,
    ids: set[str],
    paths: set[str],
) -> tuple[str, list[str]]:
    item = _require_mapping(artifact, "artifact")
    artifact_id = item.get("id")
    if not isinstance(artifact_id, str) or not artifact_id:
        raise ManifestError("artifact id must be a non-empty string")
    if artifact_id in ids:
        raise ManifestError(f"duplicate artifact id: {artifact_id}")
    ids.add(artifact_id)
    relative = _safe_relative_path(item.get("path"))
    if relative in paths:
        raise ManifestError(f"duplicate normalized artifact path: {relative}")
    paths.add(relative)
    role = item.get("role")
    status = item.get("status")
    if role not in ROLES:
        raise ManifestError(f"artifact {artifact_id}: unsupported role {role!r}")
    if status not in STATUSES:
        raise ManifestError(f"artifact {artifact_id}: unsupported status {status!r}")
    producer = _require_mapping(
        item.get("producer"), f"artifact {artifact_id} producer"
    )
    for key in ("name", "version", "database"):
        if not isinstance(producer.get(key), str) or not producer[key]:
            raise ManifestError(f"artifact {artifact_id}: producer.{key} is required")
    redistribution = _require_mapping(
        item.get("redistribution"), f"artifact {artifact_id} redistribution"
    )
    for key in ("license", "permission"):
        value = redistribution.get(key)
        if status == "published" and (
            not isinstance(value, str) or value.lower() in {"unknown", "unverified"}
        ):
            raise ManifestError(
                f"artifact {artifact_id}: published redistribution.{key} is unresolved"
            )
    lineage = item.get("lineage")
    if not isinstance(lineage, list) or any(
        not isinstance(value, str) or not value for value in lineage
    ):
        raise ManifestError(f"artifact {artifact_id}: lineage must be a list of IDs")
    size = item.get("size_bytes")
    if not isinstance(size, int) or isinstance(size, bool) or size < 0:
        raise ManifestError(f"artifact {artifact_id}: invalid size_bytes")
    digest = item.get("sha256")
    if (
        not isinstance(digest, str)
        or len(digest) != HASH_LENGTH
        or digest != digest.lower()
    ):
        raise ManifestError(f"artifact {artifact_id}: invalid sha256")
    try:
        int(digest, 16)
    except ValueError as exc:
        raise ManifestError(f"artifact {artifact_id}: invalid sha256") from exc
    candidate = root / Path(*relative.split("/"))
    resolved_root = root.resolve()
    if candidate.is_symlink():
        resolved = candidate.resolve()
        if resolved != candidate.absolute():
            raise ManifestError(f"artifact {artifact_id}: symlinks are not accepted")
    try:
        candidate.resolve(strict=True).relative_to(resolved_root)
    except (OSError, ValueError) as exc:
        raise ManifestError(
            f"artifact {artifact_id}: path escapes dataset root"
        ) from exc
    if not candidate.is_file():
        raise ManifestError(f"artifact {artifact_id}: file is missing: {relative}")
    if candidate.stat().st_size != size:
        raise ManifestError(f"artifact {artifact_id}: size does not match {relative}")
    if _sha256(candidate) != digest:
        raise ManifestError(f"artifact {artifact_id}: sha256 does not match {relative}")
    return artifact_id, lineage


def validate_manifest(manifest_path: Path, root: Path) -> int:
    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ManifestError(f"unsupported manifest schema: {manifest_path}")
        artifacts = payload.get("artifacts")
        if not isinstance(artifacts, list) or not artifacts:
            raise ManifestError("manifest artifacts must be a non-empty list")
        ids: set[str] = set()
        paths: set[str] = set()
        lineage: dict[str, list[str]] = {}
        for artifact in artifacts:
            artifact_id, references = _validate_artifact(
                artifact, root=root, ids=ids, paths=paths
            )
            lineage[artifact_id] = references
        for artifact_id, references in lineage.items():
            for reference in references:
                if reference not in ids:
                    raise ManifestError(
                        f"artifact {artifact_id}: missing lineage reference {reference!r}"
                    )
        visiting: set[str] = set()
        visited: set[str] = set()

        def visit(artifact_id: str) -> None:
            if artifact_id in visiting:
                raise ManifestError("lineage graph contains a cycle")
            if artifact_id in visited:
                return
            visiting.add(artifact_id)
            for parent in lineage[artifact_id]:
                visit(parent)
            visiting.remove(artifact_id)
            visited.add(artifact_id)

        for artifact_id in lineage:
            visit(artifact_id)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, ManifestError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(f"Validated {len(lineage)} dataset artifacts from {manifest_path}.")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="dataset root used to resolve manifest-relative artifact paths",
    )
    args = parser.parse_args()
    return validate_manifest(args.manifest, args.root.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
