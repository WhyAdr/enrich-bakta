"""Translation evidence ledger loading, binding, and legacy COMMENT detection."""

from __future__ import annotations

import json
import re
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

from enrich_bakta_lib.core.genbank import (
    cds_by_locus,
    feature_uid,
    normalize_protein,
    protein_sha256,
)
from enrich_bakta_lib.core.types import (
    MergeError,
    RawDocument,
    RawFeature,
    TranslationEvidence,
    TranslationEvidenceLedger,
    sha256_bytes,
)

TRANSLATION_EVIDENCE_MARKER = "##enrich-bakta:translation-evidence:v1##"

LEGACY_RESTORATION_COMMENT_RE = re.compile(
    rb"(?:Post-processed by\s+)?normalize_baktfold\.py(?::|\s+to)\s+restore Bakta\s+provenance",
    re.IGNORECASE,
)


def _extract_comment_blocks(data: bytes) -> bytes:
    """Extract only GenBank COMMENT header sections, excluding other headers and features."""
    blocks: list[bytes] = []
    in_comment = False
    for line in data.splitlines(keepends=True):
        if line.startswith(b"COMMENT"):
            in_comment = True
            blocks.append(line)
        elif in_comment and line.startswith(b"            "):
            blocks.append(line)
        else:
            in_comment = False
    return b"".join(blocks)


def has_translation_evidence_marker(data: bytes) -> bool:
    return TRANSLATION_EVIDENCE_MARKER.encode("ascii") in data


def record_has_legacy_restoration_comment(header_bytes: bytes) -> bool:
    """Detect known legacy normalize_baktfold.py restoration comment in a record header."""
    return bool(
        LEGACY_RESTORATION_COMMENT_RE.search(_extract_comment_blocks(header_bytes))
    )


def has_legacy_restoration_comment(data: bytes) -> bool:
    """Detect known legacy normalize_baktfold.py restoration comment across data."""
    return bool(LEGACY_RESTORATION_COMMENT_RE.search(_extract_comment_blocks(data)))


def validate_faa_gbff(
    base: RawDocument,
    proteins: Mapping[str, str],
    query_ids: Iterable[str],
    *,
    source_name: str,
    translation_evidence: Mapping[str, TranslationEvidence] | None = None,
    allow_imported_translations: bool = False,
) -> tuple[dict[str, RawFeature], dict[str, Any]]:
    """Prove that every evidence query maps exactly to an identical GBFF CDS."""
    cds = cds_by_locus(base)
    unique_ids = set(query_ids)
    missing_faa = sorted(unique_ids - proteins.keys())
    missing_gbff = sorted(unique_ids - cds.keys())
    if missing_faa or missing_gbff:
        raise MergeError(
            f"{source_name} query IDs do not map exactly: "
            + json.dumps(
                {
                    "missing_from_faa": missing_faa[:25],
                    "missing_from_gbff": missing_gbff[:25],
                },
                indent=2,
            )
        )
    mismatches: list[str] = []
    missing_translations: list[str] = []
    for query_id in sorted(unique_ids):
        translations = cds[query_id].values("translation")
        if len(translations) != 1 or not normalize_protein(translations[0]):
            missing_translations.append(query_id)
        elif normalize_protein(translations[0]) != proteins[query_id]:
            mismatches.append(query_id)
    if missing_translations or mismatches:
        raise MergeError(
            f"{source_name} FAA/GBFF protein validation failed: "
            + json.dumps(
                {
                    "missing_or_ambiguous_translation": missing_translations[:25],
                    "sequence_mismatches": mismatches[:25],
                },
                indent=2,
            )
        )
    if translation_evidence is not None:
        base_hash = sha256_bytes(base.data)
        required_query_ids: Iterable[str] = getattr(
            translation_evidence, "required_query_ids", frozenset()
        )
        for query_id in sorted(unique_ids):
            evidence = translation_evidence.get(query_id)
            if evidence is None:
                # Restoration manifests are intentionally partial: they bind
                # imported/restored targets, while untouched translated CDSs
                # remain covered by the FAA/GBFF equality check above.
                if query_id in required_query_ids:
                    raise MergeError(
                        f"{source_name} translation evidence ledger is missing {query_id!r}"
                    )
                continue
            expected_hash = (
                evidence.bound_output_sha256 or evidence.original_base_sha256
            )
            if expected_hash != base_hash:
                raise MergeError(
                    f"{source_name} translation evidence ledger is bound to a different base"
                )
            if evidence.feature_uid != feature_uid(cds[query_id]):
                raise MergeError(
                    f"{source_name} translation evidence ledger disagrees for {query_id!r}"
                )
            if evidence.protein_sha256 != protein_sha256(proteins[query_id]):
                raise MergeError(
                    f"{source_name} translation evidence ledger disagrees for {query_id!r}"
                )
            if evidence.origin == "imported_faa" and not allow_imported_translations:
                raise MergeError(
                    f"{source_name} uses imported translation evidence for {query_id!r}; "
                    "pass --allow-imported-translations explicitly"
                )
    return cds, {
        "faa_records": len(proteins),
        "faa_gbff_sequence_matches": len(unique_ids),
        "missing_faa_query_ids": 0,
        "missing_gbff_query_ids": 0,
    }


def load_translation_evidence(
    path: Path, base: RawDocument
) -> dict[str, TranslationEvidence]:
    """Load and bind a JSON restoration manifest to the current GBFF bytes."""
    if path.suffix.lower() != ".json":
        raise MergeError(
            "translation evidence must use a JSON merge manifest; TSV manifests "
            "cannot preserve the structured evidence ledger"
        )
    try:
        snapshot = path.read_bytes()
        payload = json.loads(snapshot.decode("utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise MergeError(
            f"cannot read translation evidence manifest {path}: {exc}"
        ) from exc
    if not isinstance(payload, dict) or payload.get("schema") not in {
        "enrich-bakta.merge-manifest.v1",
        "enrich-bakta.merge-manifest.v2",
    }:
        raise MergeError(f"unsupported translation evidence manifest schema: {path}")
    metadata = payload.get("metadata")
    entries = payload.get("entries")
    if not isinstance(metadata, dict) or not isinstance(entries, list):
        raise MergeError(
            "translation evidence manifest has invalid metadata or entries"
        )
    current_hash = sha256_bytes(base.data)
    if metadata.get("output_sha256") != current_hash:
        raise MergeError(
            "translation evidence manifest does not describe the current GBFF bytes"
        )
    cds = cds_by_locus(base)
    result: dict[str, TranslationEvidence] = {}
    for entry in entries:
        if (
            not isinstance(entry, dict)
            or entry.get("entry_type") != "translation_restoration"
        ):
            continue
        query_id = entry.get("query_id")
        raw = entry.get("translation_evidence")
        if not isinstance(query_id, str) or not isinstance(raw, dict):
            raise MergeError(
                "translation restoration manifest entry lacks structured evidence"
            )
        if query_id in result:
            raise MergeError(f"translation evidence manifest repeats {query_id!r}")
        if query_id not in cds:
            raise MergeError(
                f"translation evidence manifest references unknown CDS {query_id!r}"
            )
        required = (
            "feature_uid",
            "protein_sha256",
            "origin",
            "artifact_sequence_match",
            "genomic_validation_status",
            "validation_reason",
            "original_base_sha256",
            "faa_sha256",
            "restoration_policy",
            "bound_output_sha256",
        )
        if any(key not in raw for key in required):
            raise MergeError(f"translation evidence for {query_id!r} is incomplete")
        if (
            type(raw["artifact_sequence_match"]) is not bool
            or raw["artifact_sequence_match"] is not True
        ):
            raise MergeError(
                f"translation evidence for {query_id!r} requires a true boolean artifact match"
            )
        if not isinstance(raw["origin"], str) or raw["origin"] not in {
            "imported_faa",
            "genomically_validated",
        }:
            raise MergeError(
                f"translation evidence for {query_id!r} has an unknown origin"
            )
        if not isinstance(raw["genomic_validation_status"], str) or raw[
            "genomic_validation_status"
        ] not in {
            "passed",
            "failed",
            "unresolved",
            "not_attempted",
        }:
            raise MergeError(
                f"translation evidence for {query_id!r} has an invalid genomic status"
            )
        expected_policy = (
            "import-faa" if raw["origin"] == "imported_faa" else "validated-only"
        )
        if raw["restoration_policy"] != expected_policy or (
            raw["origin"] == "genomically_validated"
            and raw["genomic_validation_status"] != "passed"
        ):
            raise MergeError(
                f"translation evidence for {query_id!r} has inconsistent origin/policy/status"
            )
        for key in (
            "protein_sha256",
            "original_base_sha256",
            "faa_sha256",
            "bound_output_sha256",
        ):
            if (
                not isinstance(raw[key], str)
                or re.fullmatch(r"[0-9a-f]{64}", raw[key]) is None
            ):
                raise MergeError(
                    f"translation evidence for {query_id!r} has an invalid {key}"
                )
        if not isinstance(raw.get("producer_lineage"), dict) or not isinstance(
            raw["validation_reason"], str
        ):
            raise MergeError(
                f"translation evidence for {query_id!r} has invalid lineage/reason fields"
            )
        try:
            evidence = TranslationEvidence(
                feature_uid=str(raw["feature_uid"]),
                protein_sha256=str(raw["protein_sha256"]),
                origin=str(raw["origin"]),
                artifact_sequence_match=bool(raw["artifact_sequence_match"]),
                genomic_validation_status=str(raw["genomic_validation_status"]),
                validation_reason=str(raw["validation_reason"]),
                original_base_sha256=str(raw["original_base_sha256"]),
                faa_sha256=str(raw["faa_sha256"]),
                restoration_policy=(
                    None
                    if raw["restoration_policy"] is None
                    else str(raw["restoration_policy"])
                ),
                producer_lineage=(
                    raw["producer_lineage"]
                    if isinstance(raw["producer_lineage"], dict)
                    else {}
                ),
                bound_output_sha256=str(raw["bound_output_sha256"]),
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise MergeError(
                f"translation evidence for {query_id!r} has invalid fields"
            ) from exc
        if evidence.bound_output_sha256 != current_hash:
            raise MergeError(
                f"translation evidence for {query_id!r} is bound to another output"
            )
        if evidence.feature_uid != feature_uid(cds[query_id]):
            raise MergeError(
                f"translation evidence for {query_id!r} has the wrong feature identity"
            )
        translations = cds[query_id].values("translation")
        if (
            len(translations) != 1
            or protein_sha256(translations[0]) != evidence.protein_sha256
        ):
            raise MergeError(
                f"translation evidence for {query_id!r} disagrees with the encoded protein"
            )
        result[query_id] = evidence
    if not result:
        raise MergeError(
            "translation evidence manifest contains no restoration entries"
        )
    metadata_evidence = metadata.get("translation_evidence")
    required_query_ids: set[str] = set()
    if metadata_evidence is not None:
        if not isinstance(metadata_evidence, dict) or any(
            not isinstance(query_id, str) or not isinstance(value, dict)
            for query_id, value in metadata_evidence.items()
        ):
            raise MergeError(
                "translation evidence manifest has invalid metadata ledger"
            )
        required_query_ids = set(metadata_evidence)
        if required_query_ids != set(result):
            raise MergeError(
                "translation evidence manifest metadata and entries disagree"
            )
        if any(
            metadata_evidence[query_id] != evidence.as_dict()
            for query_id, evidence in result.items()
        ):
            raise MergeError(
                "translation evidence manifest metadata and entries disagree in their values"
            )
    if has_translation_evidence_marker(base.data):
        marker_ids: set[str] = set()
        for record in base.records:
            prefix = base.data[record.start : record.features_offset].decode("utf-8")
            content = "\n".join(line[12:] for line in prefix.splitlines())
            for block in content.split(TRANSLATION_EVIDENCE_MARKER)[1:]:
                match = re.search(
                    r"Restored CDSs:\s*(.*?)\s*Original base SHA-256:", block, re.DOTALL
                )
                if match is None:
                    raise MergeError(
                        "translation evidence marker lacks a restored-CDS list"
                    )
                marker_ids.update(
                    token.strip()
                    for token in match.group(1).replace("\n", " ").split(",")
                )
        if not marker_ids or not marker_ids.issubset(result):
            raise MergeError(
                "translation evidence manifest omits targets identified by the GBFF marker"
            )
        required_query_ids.update(marker_ids)
    return TranslationEvidenceLedger(
        result,
        required_query_ids=required_query_ids,
        manifest_sha256=sha256_bytes(snapshot),
    )


def _bind_translation_evidence(
    rows: Iterable[dict[str, Any]], output_sha256: str
) -> list[dict[str, Any]]:
    bound_rows: list[dict[str, Any]] = []
    for row in rows:
        copied = dict(row)
        evidence = copied.get("translation_evidence")
        if isinstance(evidence, dict):
            evidence = dict(evidence)
            evidence["bound_output_sha256"] = output_sha256
            copied["translation_evidence"] = evidence
        bound_rows.append(copied)
    return bound_rows


def _bind_metadata_translation_evidence(
    metadata: dict[str, Any], output_sha256: str
) -> dict[str, Any]:
    copied = dict(metadata)
    evidence = copied.get("translation_evidence")
    if isinstance(evidence, dict):
        copied["translation_evidence"] = {
            query_id: (
                {
                    **value,
                    "bound_output_sha256": output_sha256,
                }
                if isinstance(value, dict)
                else value
            )
            for query_id, value in evidence.items()
        }
    return copied
