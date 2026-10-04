"""Assertions shared by the durable follow-through scientific gates."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from enrich_bakta_lib.core.merge_engine import (
    TOOL_VERSION,
    load_translation_evidence,
    parse_genbank_bytes,
)

EXPECTED_POLICIES = {
    "allow_imported_translations": True,
    "baktfold_invalid_ec_policy": "skip",
    "eggnog_expected_version": "3.0.0-beta6",
    "gene_conflict_policy": "skip",
    "kofamscan_version_baktfold_kofam": None,
    "kofamscan_version_unified": "1.3.0",
    "merge_timestamp": None,
    "translation_policy": "import-faa",
}


def load_baseline(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema_version") != 1:
        raise AssertionError("scientific baseline schema_version must be 1")
    if payload.get("tool") != {"name": "enrich-bakta", "version": TOOL_VERSION}:
        raise AssertionError(
            f"scientific baseline tool identity does not match enrich-bakta {TOOL_VERSION}"
        )
    if payload.get("policies") != EXPECTED_POLICIES:
        raise AssertionError("scientific baseline policies do not match gate policies")
    return payload


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_baseline_inputs(root: Path, baseline: dict[str, Any]) -> dict[str, Any]:
    publication_path = root / "docs" / "data" / "PUBLISHED-MANIFEST.json"
    publication = json.loads(publication_path.read_text(encoding="utf-8"))
    published_by_path = {row["path"]: row for row in publication["artifacts"]}
    checked = 0
    for sample, sample_baseline in baseline["samples"].items():
        for input_name, expected in sample_baseline["inputs"].items():
            relative_path = expected["path"]
            actual_path = root / relative_path
            actual_sha256 = _sha256(actual_path)
            if actual_sha256 != expected["sha256"]:
                raise AssertionError(
                    f"{sample} {input_name} SHA-256 drift: "
                    f"expected {expected['sha256']}, got {actual_sha256}"
                )
            published = published_by_path.get(relative_path)
            if published is None:
                raise AssertionError(
                    f"baseline input is not published: {relative_path}"
                )
            for field in ("sha256", "producer"):
                if published[field] != expected[field]:
                    raise AssertionError(
                        f"{sample} {input_name} {field} disagrees with "
                        "PUBLISHED-MANIFEST.json"
                    )
            checked += 1
    return {"inputs": checked, "tool_version": TOOL_VERSION}


def assert_expected(label: str, actual: Any, expected: Any) -> None:
    if actual != expected:
        raise AssertionError(f"{label} drift: expected {expected!r}, got {actual!r}")


def assert_imported_evidence(
    manifest_path: Path, output_path: Path, expected_query_ids: list[str]
) -> dict[str, Any]:
    ledger = load_translation_evidence(
        manifest_path, parse_genbank_bytes(output_path.read_bytes())
    )
    actual_query_ids = sorted(ledger)
    expected = sorted(expected_query_ids)
    assert_expected(
        f"{manifest_path.name} imported query IDs", actual_query_ids, expected
    )
    wrong_origins = {
        query_id: evidence.origin
        for query_id, evidence in ledger.items()
        if evidence.origin != "imported_faa"
    }
    if wrong_origins:
        raise AssertionError(
            f"{manifest_path.name} contains non-imported translation origins: "
            f"{wrong_origins}"
        )
    return {"query_ids": actual_query_ids, "origin": "imported_faa"}


def assert_candidate_manifest(manifest_path: Path) -> dict[str, int]:
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    decisions = payload.get("decisions")
    entries = payload.get("entries")
    metadata = payload.get("metadata")
    if not isinstance(decisions, list) or not isinstance(entries, list):
        raise AssertionError(
            f"{manifest_path.name} lacks candidate decisions or entries"
        )
    if not isinstance(metadata, dict):
        raise AssertionError(f"{manifest_path.name} lacks metadata")

    by_id: dict[str, dict[str, Any]] = {}
    for decision in decisions:
        candidate_id = decision.get("candidate_id")
        if not isinstance(candidate_id, str) or candidate_id in by_id:
            raise AssertionError(f"{manifest_path.name} has invalid candidate IDs")
        by_id[candidate_id] = decision

    projection_count = 0
    for entry in entries:
        candidate_id = entry.get("candidate_id")
        if candidate_id is None and (
            "planned_status" in entry or "final_status" in entry
        ):
            raise AssertionError(
                f"{manifest_path.name} finalized source entry lacks candidate_id"
            )
        if candidate_id is None:
            continue
        decision = by_id.get(candidate_id)
        if decision is None:
            raise AssertionError(
                f"{manifest_path.name} entry references unknown candidate {candidate_id}"
            )
        for field in ("status", "final_status", "reason_code"):
            expected_field = "final_status" if field == "status" else field
            if entry.get(field) != decision.get(expected_field):
                raise AssertionError(
                    f"{manifest_path.name} candidate {candidate_id} projection {field} "
                    "disagrees with its final decision"
                )
        if "planned_reason_code" not in entry:
            raise AssertionError(
                f"{manifest_path.name} candidate {candidate_id} lacks planned_reason_code"
            )
        projection_count += 1

    suppressed = sum(
        isinstance(row.get("final_status"), str)
        and row["final_status"].startswith("suppressed_")
        for row in decisions
    )
    assert_expected(
        f"{manifest_path.name} suppressed_candidate_count",
        metadata.get("suppressed_candidate_count"),
        suppressed,
    )
    return {
        "candidate_decisions": len(decisions),
        "source_projections": projection_count,
        "suppressed_candidates": suppressed,
    }
