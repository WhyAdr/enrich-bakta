"""Export one evidence-preserving TSV row per GenBank CDS."""

from __future__ import annotations

import argparse
import csv
import io
import json
import sys
from collections import Counter, defaultdict
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from enrich_bakta_lib.core.artifacts import preflight_artifact_paths, stage_artifacts
from enrich_bakta_lib.core.decisions import ACCEPTED_STATUSES
from enrich_bakta_lib.core.genbank import feature_uid, parse_genbank_bytes
from enrich_bakta_lib.core.lineage import load_translation_evidence
from enrich_bakta_lib.core.types import (
    TOOL_VERSION,
    MergeError,
    RawFeature,
    TranslationEvidence,
    read_input_bytes,
    sha256_bytes,
)
from enrich_bakta_lib.sources.value_rules import (
    STRUCTURED_TOKEN_RE,
    parse_standalone_member_note,
    structured_note_tokens,
    validate_go_term,
    validate_interpro_accession,
)

_FIELDS = (
    "gbff_sha256",
    "record_index",
    "record_id",
    "feature_index",
    "location",
    "feature_uid",
    "locus_tag",
    "locus_tags",
    "protein_id",
    "product",
    "gene",
    "ec",
    "db_xref",
    "ko",
    "go",
    "interpro",
    "pfam",
    "tigrfam",
    "antifam_status",
    "antifam_accessions",
    "translation_present",
    "translation_origin",
    "translation_genomic_validation_status",
    "inference",
    "accepted_decisions",
)
_DECISION_FIELDS = (
    "candidate_id",
    "source_id",
    "source_sha256",
    "candidate_role",
    "qualifier",
    "normalized_value",
    "final_status",
)
_DECISION_QUALIFIERS = {"product", "gene", "EC_number", "db_xref", "note", "inference"}
_INTERPRO_CONTEXT = "enrich-bakta.interproscan-context.v1"
_EGGNOG_CONTEXT = "enrich-bakta.eggnog-context.v1"
_COMBINED_CONTEXT = "enrich-bakta.context.v1"


def _object(value: object, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise MergeError(f"{label} must be a JSON object")
    return value


def _rows(value: object, label: str) -> list[dict[str, Any]]:
    if not isinstance(value, list) or any(not isinstance(row, dict) for row in value):
        raise MergeError(f"{label} must be an array of JSON objects")
    return value


def _bound_json(data: bytes, label: str, gbff_hash: str) -> dict[str, Any]:
    try:
        payload = _object(json.loads(data), label)
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise MergeError(f"invalid {label} JSON: {exc}") from exc
    metadata = _object(payload.get("metadata"), f"{label} metadata")
    if metadata.get("output_sha256") != gbff_hash:
        raise MergeError(f"{label} output SHA-256 does not match the GBFF")
    return payload


def _accepted_decisions(
    rows: list[dict[str, Any]],
    uid_counts: Counter[str],
) -> dict[str, list[dict[str, str]]]:
    accepted: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        status, qualifier = row.get("final_status"), row.get("qualifier")
        if not isinstance(status, str) or not isinstance(qualifier, str):
            raise MergeError("manifest decision lacks a valid status or qualifier")
        if status not in ACCEPTED_STATUSES or qualifier not in _DECISION_QUALIFIERS:
            continue
        targets = row.get("target_feature_uids")
        if not isinstance(targets, list) or any(
            not isinstance(uid, str) for uid in targets
        ):
            raise MergeError("accepted manifest decision lacks valid feature UIDs")
        decision: dict[str, str] = {}
        for field in _DECISION_FIELDS:
            value = row.get(field)
            if not isinstance(value, str):
                raise MergeError(f"accepted manifest decision lacks a valid {field}")
            decision[field] = value
        for uid in dict.fromkeys(targets):
            if uid_counts[uid] > 1:
                raise MergeError(
                    f"accepted manifest decision targets an ambiguous CDS UID {uid}"
                )
            if uid in uid_counts:
                accepted[uid].append(decision)
    return accepted


def _antifam_context(context: dict[str, Any]) -> tuple[set[str], dict[str, list[str]]]:
    schema = context.get("schema")
    if schema == _EGGNOG_CONTEXT:
        return set(), {}
    if schema == _COMBINED_CONTEXT:
        context = _object(context.get("interproscan"), "combined InterProScan context")
        schema = context.get("schema")
    if schema != _INTERPRO_CONTEXT:
        raise MergeError("unsupported context report schema")

    metadata = _object(context.get("metadata"), "InterProScan context metadata")
    qc = _object(metadata.get("qc_flags", {}), "InterProScan QC flags")
    if "anti_fam" not in qc:
        return set(), {}
    queries = set()
    for row in _rows(context.get("entries"), "InterProScan context entries"):
        query = row.get("query_id")
        if not isinstance(query, str) or not query:
            raise MergeError("InterProScan context entry lacks a valid query ID")
        queries.add(query)
    hits: dict[str, list[str]] = defaultdict(list)
    for row in _rows(qc["anti_fam"], "AntiFam alerts"):
        query = row.get("query_id")
        accession = row.get("signature_accession")
        if (
            not isinstance(query, str)
            or not query
            or not isinstance(accession, str)
            or not accession
        ):
            raise MergeError(
                "AntiFam alert lacks a valid query ID or signature accession"
            )
        hits[query].append(accession)
    return queries, hits


def _annotations(cds: RawFeature) -> dict[str, list[str]]:
    fields: dict[str, list[str]] = {
        name: [] for name in ("ko", "go", "interpro", "pfam", "tigrfam")
    }
    notes = cds.values("note")
    xrefs = cds.values("db_xref")
    for token in [
        *xrefs,
        *(token for note in notes for token in structured_note_tokens(note)),
    ]:
        prefix, _, identifier = token.partition(":")
        if prefix == "KEGG" and STRUCTURED_TOKEN_RE.fullmatch(token):
            fields["ko"].append(token)
        elif validate_go_term(token) is not None:
            fields["go"].append(token)
        elif (
            prefix == "InterPro" and validate_interpro_accession(identifier) is not None
        ):
            fields["interpro"].append(token)
    for token in [*notes, *xrefs]:
        member = parse_standalone_member_note(token)
        if member is not None:
            fields[member.db.lower()].append(member.raw_note)
    return {"db_xref": xrefs, **fields}


def _cell(value: object) -> str | int:
    if isinstance(value, (str, int)):
        return value
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"), sort_keys=True)


def export_cds_tsv(
    gbff_path: Path,
    *,
    manifest_path: Path | None = None,
    context_path: Path | None = None,
) -> bytes:
    """Render CDS qualifiers and optional, byte-bound evidence as UTF-8 TSV."""
    gbff = read_input_bytes(gbff_path, "GBFF")
    gbff_hash = sha256_bytes(gbff)
    document = parse_genbank_bytes(gbff)
    cds_rows = [
        (record.index, feature_index, cds, feature_uid(cds))
        for record in document.records
        for feature_index, cds in enumerate(record.features)
        if cds.feature_type == "CDS"
    ]
    uid_counts = Counter(uid for _, _, _, uid in cds_rows)
    loci = Counter(tag for _, _, cds, _ in cds_rows for tag in set(cds.locus_tags))
    manifest_metadata: dict[str, Any] | None = None
    accepted: dict[str, list[dict[str, str]]] = {}
    translations: dict[str, TranslationEvidence] = {}
    if manifest_path is not None:
        manifest_data = read_input_bytes(manifest_path, "manifest")
        manifest = _bound_json(manifest_data, "manifest", gbff_hash)
        if manifest.get("schema") != "enrich-bakta.merge-manifest.v2":
            raise MergeError(
                "export requires an enrich-bakta.merge-manifest.v2 JSON manifest"
            )
        manifest_metadata = manifest["metadata"]
        accepted = _accepted_decisions(
            _rows(manifest.get("decisions"), "manifest decisions"), uid_counts
        )
        has_restoration = bool(manifest_metadata.get("translation_evidence")) or any(
            entry.get("entry_type") == "translation_restoration"
            for entry in _rows(manifest.get("entries"), "manifest entries")
        )
        manifest_hash = sha256_bytes(manifest_data) if has_restoration else None
        del manifest, manifest_data
        if has_restoration:
            ledger = load_translation_evidence(manifest_path, document)
            if getattr(ledger, "manifest_sha256", None) != manifest_hash:
                raise MergeError(
                    "manifest changed while reading the translation evidence"
                )
            translations = {
                evidence.feature_uid: evidence for evidence in ledger.values()
            }

    context_queries: set[str] = set()
    antifam_hits: dict[str, list[str]] = {}
    if context_path is not None:
        context_data = read_input_bytes(context_path, "context report")
        context = _bound_json(context_data, "context report", gbff_hash)
        if manifest_metadata is not None and manifest_metadata.get(
            "context_report_sha256"
        ) != sha256_bytes(context_data):
            raise MergeError("context report SHA-256 does not match the manifest link")
        context_queries, antifam_hits = _antifam_context(context)

    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(_FIELDS)
    for record_index, feature_index, cds, uid in cds_rows:
        locus = cds.locus_tag
        antifam_status = "unavailable"
        accessions: list[str] = []
        if locus is not None and (locus in context_queries or locus in antifam_hits):
            if loci[locus] != 1:
                antifam_status = "ambiguous_locus"
            else:
                accessions = antifam_hits.get(locus, [])
                antifam_status = "hit" if accessions else "not_reported"
        translation = translations.get(uid)
        origin = "unavailable" if manifest_path is None else "not_recorded"
        row = {
            "gbff_sha256": gbff_hash,
            "record_index": record_index,
            "record_id": cds.record_id,
            "feature_index": feature_index,
            "location": cds.location,
            "feature_uid": uid,
            "locus_tag": locus or "",
            "locus_tags": cds.values("locus_tag"),
            "protein_id": cds.values("protein_id"),
            "product": cds.values("product"),
            "gene": cds.values("gene"),
            "ec": cds.values("EC_number"),
            **_annotations(cds),
            "antifam_status": antifam_status,
            "antifam_accessions": accessions,
            "translation_present": "true"
            if any(cds.values("translation"))
            else "false",
            "translation_origin": translation.origin if translation else origin,
            "translation_genomic_validation_status": translation.genomic_validation_status
            if translation
            else "unavailable",
            "inference": cds.values("inference"),
            "accepted_decisions": accepted.get(uid, [])
            if manifest_path is not None
            else None,
        }
        writer.writerow(_cell(row[field]) for field in _FIELDS)
    return output.getvalue().encode("utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {TOOL_VERSION}"
    )
    parser.add_argument(
        "--gbff", type=Path, required=True, help="annotated GenBank flat file"
    )
    parser.add_argument(
        "--manifest", type=Path, help="JSON merge manifest bound to the GBFF"
    )
    parser.add_argument(
        "--context-report", type=Path, help="context JSON bound to the GBFF"
    )
    parser.add_argument(
        "--output", type=Path, help="TSV path (default: standard output)"
    )
    args = parser.parse_args(argv)
    inputs = [
        path
        for path in (args.gbff, args.manifest, args.context_report)
        if path is not None
    ]
    try:
        if args.output is not None:
            preflight_artifact_paths([args.output], inputs)
        data = export_cds_tsv(
            args.gbff, manifest_path=args.manifest, context_path=args.context_report
        )
        if args.output is None:
            sys.stdout.buffer.write(data)
        else:
            stage_artifacts([(args.output, data)], inputs)
    except (MergeError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
