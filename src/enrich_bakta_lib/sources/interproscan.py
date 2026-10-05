"""InterProScan evidence parser, layout validation, and candidate planner."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from enrich_bakta_lib.core.decisions import (
    build_candidate_ledger,
)
from enrich_bakta_lib.core.merge_engine import (
    TOOL_VERSION,
    Insertion,
    MergeError,
    RawDocument,
    TranslationEvidence,
    cds_by_locus,
    comment_insertion,
    feature_uid,
    finalize_merge,
    has_translation_evidence_marker,
    load_translation_evidence,
    normalize_protein,
    parse_faa,
    parse_genbank_bytes,
    protein_sha256,
    qualifier_insertion,
    read_input_bytes,
    reconcile_insertions,
    record_has_legacy_restoration_comment,
    sha256_bytes,
    validate_genbank_semantics,
)
from enrich_bakta_lib.sources.value_rules import (
    find_member_note_witness,
    validate_go_term,
    validate_interpro_accession,
    validate_pfam_accession,
    validate_tigrfam_accession,
)

SUPPORTED_INTERPROSCAN_VERSIONS = frozenset({"5.59-91.0"})
DEFAULT_INTERPROSCAN_VERSION = "5.59-91.0"

TSV_LAYOUT_PROFILES = ("ipr-go-pathways", "ipr-go", "ipr-pathways", "ipr-only")
DEFAULT_TSV_LAYOUT = "ipr-go-pathways"

SUPPORTED_MEMBER_DBS = frozenset({"Pfam", "TIGRFAM"})
DEFAULT_MEMBER_DBS = ("Pfam", "TIGRFAM")

INTERPROSCAN_COMMENT_MARKER = "##enrich-bakta:interproscan:v1##"

MAX_RAW_LINE_BYTES = 1024 * 1024  # 1 MiB line ceiling
MAX_AGGREGATED_QUERIES = 200_000
MAX_AGGREGATED_CANDIDATES = 500_000

MD5_HEX_RE = re.compile(r"^[0-9a-fA-F]{32}\Z")


@dataclass(frozen=True)
class InterProScanHit:
    """Parsed, validated hit from an InterProScan TSV row."""

    query_id: str
    md5: str
    length: int
    analysis: str
    signature_accession: str
    signature_description: str
    start: int
    end: int
    score: str
    status: str
    date: str
    interpro_accession: str | None
    interpro_description: str | None
    go_terms: tuple[str, ...]
    pathways: tuple[str, ...]
    row_number: int


@dataclass(frozen=True)
class InterProScanPlan:
    """Complete candidate plan produced from an InterProScan TSV input."""

    insertions: list[Insertion]
    evidence_rows: list[dict[str, Any]]
    stats: dict[str, Any]
    context_report: dict[str, Any]
    source_sha256: str
    version: str
    layout: str
    member_dbs: tuple[str, ...]


def parse_member_dbs(value: str | Sequence[str]) -> tuple[str, ...]:
    """Parse and validate member database promotion allowlist."""
    if isinstance(value, str):
        parts = [item.strip() for item in value.split(",") if item.strip()]
    else:
        parts = [item.strip() for item in value if item.strip()]
    for item in parts:
        if item not in SUPPORTED_MEMBER_DBS:
            raise MergeError(
                f"unsupported InterProScan member database {item!r}; "
                f"supported databases: {', '.join(sorted(SUPPORTED_MEMBER_DBS))}"
            )
    return tuple(dict.fromkeys(parts))


def validate_tsv_row_layout(
    fields: Sequence[str],
    layout: str,
    row_number: int,
) -> tuple[bool, str | None, str | None, tuple[str, ...], tuple[str, ...]]:
    """Validate column count and placeholder coherence for a TSV row.

    Returns (is_integrated, ipr_acc, ipr_desc, go_terms, pathways).
    """
    if layout not in TSV_LAYOUT_PROFILES:
        raise MergeError(
            f"unsupported InterProScan TSV layout {layout!r}; "
            f"supported layouts: {', '.join(TSV_LAYOUT_PROFILES)}"
        )
    num_fields = len(fields)
    if num_fields == 11:
        raise MergeError(
            f"InterProScan TSV row {row_number}: [unsupported_tsv_layout] has 11 columns (lookup disabled); "
            "this native layout is unsupported in this version"
        )
    if num_fields < 13:
        raise MergeError(
            f"InterProScan TSV row {row_number}: [malformed_tsv_row] is truncated: has {num_fields} columns, expected at least 13"
        )

    # Check InterPro accession (field 12, index 11) and description (field 13, index 12)
    raw_ipr = fields[11].strip()
    raw_desc = fields[12].strip()

    if raw_ipr == "-":
        if raw_desc != "-":
            raise MergeError(
                f"InterProScan TSV row {row_number}: [placeholder_incoherence] "
                f"InterPro accession is '-' but description is {raw_desc!r}"
            )
        is_integrated = False
        ipr_acc = None
        ipr_desc = None
    else:
        validated_ipr = validate_interpro_accession(raw_ipr)
        if validated_ipr is None:
            raise MergeError(
                f"InterProScan TSV row {row_number}: [invalid_interpro_accession] "
                f"has invalid InterPro accession {raw_ipr!r}"
            )
        if raw_desc == "-":
            raise MergeError(
                f"InterProScan TSV row {row_number}: [placeholder_incoherence] "
                f"InterPro accession is {raw_ipr!r} but description is '-'"
            )
        is_integrated = True
        ipr_acc = validated_ipr
        ipr_desc = raw_desc

    # Determine expected width for the selected layout
    if not is_integrated:
        expected = 13
    elif layout == "ipr-go-pathways":
        expected = 15
    elif layout in {"ipr-go", "ipr-pathways"}:
        expected = 14
    elif layout == "ipr-only":
        expected = 13
    else:
        expected = 13

    if num_fields != expected:
        raise MergeError(
            f"InterProScan TSV row {row_number}: [column_count_mismatch] has {num_fields} columns; "
            f"expected {expected} for layout {layout!r} ({'integrated' if is_integrated else 'unintegrated'})"
        )

    go_terms: list[str] = []
    pathways: list[str] = []

    if is_integrated:
        if layout == "ipr-go-pathways":
            raw_go = fields[13].strip()
            raw_path = fields[14].strip()
            if raw_go and raw_go != "-":
                for token in raw_go.split("|"):
                    tok = token.strip()
                    if tok:
                        if not validate_go_term(tok):
                            raise MergeError(
                                f"InterProScan TSV row {row_number}: [invalid_go_term] has invalid GO term {tok!r}"
                            )
                        go_terms.append(tok)
            if raw_path and raw_path != "-":
                for token in raw_path.split("|"):
                    tok = token.strip()
                    if tok:
                        pathways.append(tok)
        elif layout == "ipr-go":
            raw_go = fields[13].strip()
            if raw_go and raw_go != "-":
                for token in raw_go.split("|"):
                    tok = token.strip()
                    if tok:
                        if not validate_go_term(tok):
                            raise MergeError(
                                f"InterProScan TSV row {row_number}: [invalid_go_term] has invalid GO term {tok!r}"
                            )
                        go_terms.append(tok)
        elif layout == "ipr-pathways":
            raw_path = fields[13].strip()
            if raw_path and raw_path != "-":
                for token in raw_path.split("|"):
                    tok = token.strip()
                    if tok:
                        pathways.append(tok)

    return is_integrated, ipr_acc, ipr_desc, tuple(go_terms), tuple(pathways)


def stream_interproscan_tsv(
    tsv_path: Path,
    layout: str,
) -> tuple[str, list[InterProScanHit], dict[str, Any]]:
    """Stream binary InterProScan TSV in a single pass, hashing and parsing.

    Enforces 1 MiB line ceiling, strict UTF-8, placeholder coherence, and bounds.
    Returns (source_sha256, hits, diagnostics).
    """
    if not tsv_path.is_file():
        raise MergeError(f"InterProScan input file not found: {tsv_path}")

    hasher = hashlib.sha256()
    pathway_hasher = hashlib.sha256()
    hits: list[InterProScanHit] = []
    query_metadata: dict[str, tuple[str, int]] = {}
    pathway_counters: dict[str, Any] = {
        "total_occurrences": 0,
        "by_database": defaultdict(int),
    }
    anti_fam_hits: list[dict[str, Any]] = []
    integrated_count = 0
    unintegrated_count = 0
    row_number = 0

    with tsv_path.open("rb") as stream:
        while True:
            line_bytes = stream.readline(MAX_RAW_LINE_BYTES + 1)
            if not line_bytes:
                break
            if len(line_bytes) > MAX_RAW_LINE_BYTES and not (
                line_bytes.endswith(b"\n") or line_bytes.endswith(b"\r")
            ):
                raise MergeError(
                    f"InterProScan TSV line {row_number + 1} exceeds maximum allowed size of 1 MiB"
                )
            hasher.update(line_bytes)
            row_number += 1

            try:
                line_text = line_bytes.decode("utf-8")
            except UnicodeDecodeError as exc:
                raise MergeError(
                    f"InterProScan TSV line {row_number} contains invalid UTF-8 at byte {exc.start}"
                ) from exc

            if "\x00" in line_text:
                raise MergeError(
                    f"InterProScan TSV line {row_number} contains embedded NUL character"
                )

            # Strip only line endings, preserving empty fields
            if line_text.endswith("\r\n"):
                line_text = line_text[:-2]
            elif line_text.endswith("\n") or line_text.endswith("\r"):
                line_text = line_text[:-1]

            if not line_text:
                continue

            fields = line_text.split("\t")
            (
                is_integrated,
                ipr_acc,
                ipr_desc,
                go_terms,
                pathways,
            ) = validate_tsv_row_layout(fields, layout, row_number)

            if is_integrated:
                integrated_count += 1
            else:
                unintegrated_count += 1

            query_id = fields[0].strip()
            if not query_id:
                raise MergeError(
                    f"InterProScan TSV line {row_number} has an empty query ID"
                )

            raw_md5 = fields[1].strip()
            if not MD5_HEX_RE.fullmatch(raw_md5):
                raise MergeError(
                    f"InterProScan TSV line {row_number} has invalid MD5 {raw_md5!r}"
                )
            md5 = raw_md5.lower()

            try:
                length = int(fields[2].strip())
            except ValueError as exc:
                raise MergeError(
                    f"InterProScan TSV line {row_number} has non-integer length {fields[2]!r}"
                ) from exc
            if length <= 0:
                raise MergeError(
                    f"InterProScan TSV line {row_number} has non-positive length {length}"
                )

            analysis = fields[3].strip()
            if not analysis:
                raise MergeError(
                    f"InterProScan TSV line {row_number} has an empty analysis field"
                )

            sig_acc = fields[4].strip()
            if not sig_acc:
                raise MergeError(
                    f"InterProScan TSV line {row_number} has an empty signature accession"
                )
            sig_desc = fields[5].strip()

            try:
                start = int(fields[6].strip())
                end = int(fields[7].strip())
            except ValueError as exc:
                raise MergeError(
                    f"InterProScan TSV line {row_number} has non-integer coordinates"
                ) from exc

            if not (1 <= start <= end <= length):
                raise MergeError(
                    f"InterProScan TSV line {row_number} has invalid coordinates {start}..{end} "
                    f"for protein length {length}"
                )

            score = fields[8].strip()
            status = fields[9].strip()
            date = fields[10].strip()

            # Repeated-query consistency check
            if query_id in query_metadata:
                prev_md5, prev_len = query_metadata[query_id]
                if prev_md5 != md5 or prev_len != length:
                    raise MergeError(
                        f"InterProScan TSV line {row_number}: [query_conflict] query {query_id!r} "
                        f"has conflicting MD5 ({md5} vs {prev_md5}) or length ({length} vs {prev_len})"
                    )
            else:
                if len(query_metadata) >= MAX_AGGREGATED_QUERIES:
                    raise MergeError(
                        f"InterProScan TSV exceeded maximum query aggregation limit ({MAX_AGGREGATED_QUERIES})"
                    )
                query_metadata[query_id] = (md5, length)

            # Stream pathway occurrences
            for token in pathways:
                pathway_counters["total_occurrences"] += 1
                db_name = (
                    token.partition(":")[0].strip() if ":" in token else "unspecified"
                )
                pathway_counters["by_database"][db_name] += 1
                pathway_hasher.update(f"{query_id}\t{token}\n".encode("utf-8"))

            # QC flags (AntiFam)
            if analysis == "AntiFam" or sig_acc.startswith("ANF"):
                anti_fam_hits.append(
                    {
                        "query_id": query_id,
                        "signature_accession": sig_acc,
                        "signature_description": sig_desc,
                        "score": score,
                        "start": start,
                        "end": end,
                        "row_number": row_number,
                    }
                )

            if len(hits) >= MAX_AGGREGATED_CANDIDATES:
                raise MergeError(
                    f"InterProScan TSV exceeded maximum hit aggregation limit ({MAX_AGGREGATED_CANDIDATES})"
                )

            hits.append(
                InterProScanHit(
                    query_id=query_id,
                    md5=md5,
                    length=length,
                    analysis=analysis,
                    signature_accession=sig_acc,
                    signature_description=sig_desc,
                    start=start,
                    end=end,
                    score=score,
                    status=status,
                    date=date,
                    interpro_accession=ipr_acc,
                    interpro_description=ipr_desc,
                    go_terms=go_terms,
                    pathways=pathways,
                    row_number=row_number,
                )
            )

    source_sha256 = hasher.hexdigest()
    pathway_digest = pathway_hasher.hexdigest()
    pathway_counters["by_database"] = dict(pathway_counters["by_database"])

    diagnostics = {
        "row_count": row_number,
        "hits_count": len(hits),
        "integrated_count": integrated_count,
        "unintegrated_count": unintegrated_count,
        "query_metadata": query_metadata,
        "pathway_counters": pathway_counters,
        "pathway_digest": pathway_digest,
        "anti_fam_hits": anti_fam_hits,
    }
    return source_sha256, hits, diagnostics


def plan_interproscan(
    base: RawDocument,
    faa_path: Path,
    interproscan_path: Path,
    *,
    version: str,
    layout: str = DEFAULT_TSV_LAYOUT,
    member_dbs: tuple[str, ...] | Sequence[str] = DEFAULT_MEMBER_DBS,
    add_comment_note: bool = True,
    add_feature_provenance: bool = True,
    merge_timestamp: str | None = None,
    translation_evidence: Mapping[str, TranslationEvidence] | None = None,
    allow_imported_translations: bool = False,
    faa_data: bytes | None = None,
    faa_hash: str | None = None,
    proteins: Mapping[str, str] | None = None,
) -> InterProScanPlan:
    """Plan candidate insertions and evidence rows from InterProScan TSV."""
    if not version:
        raise MergeError(
            "InterProScan input requires an explicit interproscan_version assertion"
        )
    if version not in SUPPORTED_INTERPROSCAN_VERSIONS:
        raise MergeError(
            f"unsupported InterProScan version {version!r}; verified version is 5.59-91.0"
        )
    if layout not in TSV_LAYOUT_PROFILES:
        raise MergeError(
            f"unsupported InterProScan TSV layout {layout!r}; "
            f"supported layouts: {', '.join(TSV_LAYOUT_PROFILES)}"
        )
    member_dbs = tuple(member_dbs)

    # Known-legacy restoration preflight: inspect COMMENT/header per record
    for rec in base.records:
        rec_header = base.data[rec.start : rec.features_offset]
        if record_has_legacy_restoration_comment(rec_header):
            if not has_translation_evidence_marker(rec_header):
                raise MergeError(
                    "InterProScan preflight rejected base GenBank file: known legacy restoration provenance "
                    "detected (normalize_baktfold.py) without verified translation-evidence marker; "
                    "bound translation-evidence manifest is required"
                )

    if faa_data is not None:
        faa_bytes = faa_data
        actual_faa_hash = faa_hash or sha256_bytes(faa_data)
        actual_proteins = proteins if proteins is not None else parse_faa(faa_data)
        if faa_path and faa_path.is_file():
            disk_bytes = read_input_bytes(faa_path, "FAA")
            if sha256_bytes(disk_bytes) != actual_faa_hash:
                raise MergeError(f"FAA input {faa_path} was modified after capture")
    else:
        faa_bytes = read_input_bytes(faa_path, "FAA")
        actual_faa_hash = sha256_bytes(faa_bytes)
        actual_proteins = parse_faa(faa_bytes)

    source_sha256, hits, tsv_diag = stream_interproscan_tsv(interproscan_path, layout)

    cds = cds_by_locus(base)
    queried_ids = sorted(tsv_diag["query_metadata"].keys())

    missing_faa = [q for q in queried_ids if q not in actual_proteins]
    missing_gbff = [q for q in queried_ids if q not in cds]
    if missing_faa or missing_gbff:
        raise MergeError(
            f"InterProScan query IDs do not map exactly: {json.dumps({'missing_from_faa': missing_faa[:25], 'missing_from_gbff': missing_gbff[:25]}, indent=2)}"
        )

    # Validate FAA sequence against TSV MD5 and length
    mismatches: list[str] = []
    for query in queried_ids:
        seq = actual_proteins[query]
        expected_md5, expected_len = tsv_diag["query_metadata"][query]
        if (
            len(seq) != expected_len
            or hashlib.md5(seq.encode("ascii")).hexdigest() != expected_md5
        ):
            mismatches.append(query)
    if mismatches:
        raise MergeError(
            f"InterProScan FAA protein validation failed: {json.dumps({'sequence_mismatches': mismatches[:25]}, indent=2)}"
        )

    # Validate encoded GBFF /translation against FAA
    missing_translations: list[str] = []
    seq_mismatches: list[str] = []
    for query in queried_ids:
        translations = cds[query].values("translation")
        if len(translations) != 1 or not normalize_protein(translations[0]):
            missing_translations.append(query)
        elif normalize_protein(translations[0]) != actual_proteins[query]:
            seq_mismatches.append(query)
    if missing_translations or seq_mismatches:
        raise MergeError(
            f"InterProScan FAA/GBFF protein validation failed: {json.dumps({'missing_or_ambiguous_translation': missing_translations[:25], 'sequence_mismatches': seq_mismatches[:25]}, indent=2)}"
        )

    # Per-record legacy provenance check for queried CDSs
    record_by_id = {rec.record_id: rec for rec in base.records}
    for query in queried_ids:
        feature = cds[query]
        rec = record_by_id[feature.record_id]
        rec_header = base.data[rec.start : rec.features_offset]
        if record_has_legacy_restoration_comment(rec_header):
            if translation_evidence is None or query not in translation_evidence:
                raise MergeError(
                    f"InterProScan rejected unverified legacy translation for {query!r}: "
                    "genuine translation evidence is required"
                )

    # Translation evidence lifecycle checks
    if translation_evidence is not None:
        base_hash = sha256_bytes(base.data)
        required_query_ids: Iterable[str] = getattr(
            translation_evidence, "required_query_ids", frozenset()
        )
        for query_id in queried_ids:
            evidence = translation_evidence.get(query_id)
            if evidence is None:
                if query_id in required_query_ids:
                    raise MergeError(
                        f"InterProScan translation evidence ledger is missing {query_id!r}"
                    )
                continue
            expected_hash = (
                evidence.bound_output_sha256 or evidence.original_base_sha256
            )
            if expected_hash != base_hash:
                raise MergeError(
                    "InterProScan translation evidence ledger is bound to a different base"
                )
            if evidence.feature_uid != feature_uid(cds[query_id]):
                raise MergeError(
                    f"InterProScan translation evidence ledger disagrees for {query_id!r}"
                )
            if evidence.protein_sha256 != protein_sha256(actual_proteins[query_id]):
                raise MergeError(
                    f"InterProScan translation evidence ledger disagrees for {query_id!r}"
                )
            if evidence.origin == "imported_faa" and not allow_imported_translations:
                raise MergeError(
                    f"InterProScan uses imported translation evidence for {query_id!r}; "
                    "pass --allow-imported-translations explicitly"
                )

    # Group hits by query
    hits_by_query: dict[str, list[InterProScanHit]] = defaultdict(list)
    for hit in hits:
        hits_by_query[hit.query_id].append(hit)

    insertions: list[Insertion] = []
    evidence_rows: list[dict[str, Any]] = []
    order = 0
    support_hasher = hashlib.sha256()

    new_qualifier_features: set[str] = set()
    per_query_context: dict[str, dict[str, Any]] = {}

    for query in queried_ids:
        feature = cds[query]
        f_uid = feature_uid(feature)
        query_hits = hits_by_query.get(query, [])

        # Proposals: (field, qualifier, normalized_value, raw_value, role, evidence_class, first_row)
        proposals: list[tuple[str, str, str, str, str, str, int]] = []
        context_sigs: list[str] = []

        for hit in query_hits:
            # 1. InterPro accession
            if hit.interpro_accession:
                proposals.append(
                    (
                        "InterPro",
                        "db_xref",
                        f"InterPro:{hit.interpro_accession}",
                        hit.interpro_accession,
                        "functional_proposal",
                        "InterPro",
                        hit.row_number,
                    )
                )

            # 2. GO terms
            for go in hit.go_terms:
                proposals.append(
                    (
                        "GO",
                        "db_xref",
                        go,
                        go,
                        "functional_proposal",
                        "GO",
                        hit.row_number,
                    )
                )

            # 3. Member signatures (Pfam, TIGRFAM)
            if hit.analysis == "Pfam" and "Pfam" in member_dbs:
                res_pfam = validate_pfam_accession(hit.signature_accession)
                if res_pfam:
                    canonical, _ = res_pfam
                    proposals.append(
                        (
                            "Pfam",
                            "note",
                            f"PFAM:{canonical}",
                            hit.signature_accession,
                            "substantive_evidence",
                            "member:Pfam",
                            hit.row_number,
                        )
                    )
            elif hit.analysis == "TIGRFAM" and "TIGRFAM" in member_dbs:
                res_tigr = validate_tigrfam_accession(hit.signature_accession)
                if res_tigr:
                    canonical, _ = res_tigr
                    proposals.append(
                        (
                            "TIGRFAM",
                            "note",
                            f"TIGRFAM:{canonical}",
                            hit.signature_accession,
                            "substantive_evidence",
                            "member:TIGRFAM",
                            hit.row_number,
                        )
                    )
            else:
                context_sigs.append(f"{hit.analysis}:{hit.signature_accession}")

        # Unique proposals per (field, qualifier, normalized_value)
        unique_proposals: dict[
            tuple[str, str, str], tuple[str, str, str, str, str, str, int]
        ] = {}
        for prop in proposals:
            key = (prop[0], prop[1], prop[2])
            if key not in unique_proposals:
                unique_proposals[key] = prop

        for (
            field,
            qualifier,
            normalized_val,
            raw_val,
            role,
            ev_class,
            first_row,
        ) in sorted(unique_proposals.values(), key=lambda p: (p[1], p[2])):
            # Support digest accumulation
            support_hasher.update(
                f"{query}\t{field}\t{normalized_val}\n".encode("utf-8")
            )

            # Check existing base annotation
            witness: str | None = None
            if qualifier == "db_xref":
                existing = normalized_val in feature.values("db_xref")
                witness = normalized_val if existing else None
            elif qualifier == "note":
                witness = find_member_note_witness(
                    feature.values("note"), normalized_val
                )
                existing = witness is not None
            else:
                existing = False

            if existing:
                status = "supported_existing"
                reason_code = "value_already_present"
                emitted_val = ""
            else:
                status = "emitted"
                reason_code = "inserted"
                emitted_val = f'/{qualifier}="{normalized_val}"'
                insertions.append(
                    qualifier_insertion(
                        base.data,
                        feature,
                        qualifier,
                        normalized_val,
                        "InterProScan",
                        raw_val,
                        order=order,
                        candidate_role=role,
                        evidence_class=ev_class,
                    )
                )
                order += 1
                new_qualifier_features.add(f_uid)

            evidence_rows.append(
                {
                    "entry_type": "interproscan_candidate",
                    "query_id": query,
                    "record": feature.record_id,
                    "feature_type": feature.feature_type,
                    "locus_tag": feature.locus_tag,
                    "target_feature_uids": [f_uid],
                    "field": field,
                    "qualifier": qualifier,
                    "raw_value": raw_val,
                    "normalized_value": normalized_val,
                    "status": status,
                    "planned_status": status,
                    "final_status": status,
                    "candidate_role": role,
                    "evidence_class": ev_class,
                    "reason_code": reason_code,
                    "planned_reason_code": reason_code,
                    "emitted_qualifiers": emitted_val,
                    "support_witness": witness,
                    "row_number": first_row,
                    "source": "InterProScan",
                    "reason": f"{field} match {'already present' if existing else 'planned for insertion'}",
                }
            )

        if context_sigs or any(h.pathways for h in query_hits):
            per_query_context[query] = {
                "query_id": query,
                "context_signatures": sorted(set(context_sigs))[:50],
                "pathway_occurrences": sum(len(h.pathways) for h in query_hits),
            }

    # Add feature provenance inference for CDSs with new insertions
    if add_feature_provenance:
        inference_value = f"protein motif:InterProScan:{version}"
        for query in sorted(queried_ids):
            feature = cds[query]
            f_uid = feature_uid(feature)
            if f_uid not in new_qualifier_features:
                continue
            if inference_value not in feature.values("inference"):
                insertions.append(
                    qualifier_insertion(
                        base.data,
                        feature,
                        "inference",
                        inference_value,
                        "InterProScan",
                        inference_value,
                        order=order,
                        candidate_role="producer_provenance",
                        evidence_class="InterProScan:producer",
                    )
                )
                order += 1
                evidence_rows.append(
                    {
                        "entry_type": "interproscan_candidate",
                        "query_id": query,
                        "record": feature.record_id,
                        "feature_type": feature.feature_type,
                        "locus_tag": feature.locus_tag,
                        "target_feature_uids": [f_uid],
                        "field": "inference",
                        "qualifier": "inference",
                        "raw_value": inference_value,
                        "normalized_value": inference_value,
                        "status": "emitted",
                        "planned_status": "emitted",
                        "final_status": "emitted",
                        "candidate_role": "producer_provenance",
                        "evidence_class": "InterProScan:producer",
                        "reason_code": "inserted",
                        "planned_reason_code": "inserted",
                        "emitted_qualifiers": f'/inference="{inference_value}"',
                        "row_number": 0,
                        "source": "InterProScan",
                        "reason": "InterProScan feature provenance inference",
                    }
                )

    # COMMENT note insertion
    if add_comment_note:
        lines = [
            f"Source InterProScan {version}; layout={layout}; tsv sha256={source_sha256[:16]}",
            f"FAA sha256={actual_faa_hash[:16]}; member-dbs={','.join(member_dbs) or 'none'}.",
            "InterProScan calls are sequence-feature evidence, not proof of function.",
        ]
        if merge_timestamp:
            lines.append(f"Merge timestamp: {merge_timestamp}")
        for record in base.records:
            c_ins = comment_insertion(
                base,
                record,
                f"{INTERPROSCAN_COMMENT_MARKER}:{source_sha256[:12]}:{actual_faa_hash[:12]}",
                lines,
                "InterProScan provenance",
                order,
            )
            if c_ins is not None:
                insertions.append(c_ins)
                order += 1

    support_digest = support_hasher.hexdigest()

    stats = {
        "rows": tsv_diag["row_count"],
        "hits": tsv_diag["hits_count"],
        "queries": len(tsv_diag["query_metadata"]),
        "matched_queries": len(queried_ids),
        "planned_insertions": len(insertions),
        "interproscan_sha256": source_sha256,
        "faa_sha256": faa_hash,
    }

    context_report = {
        "schema": "enrich-bakta.interproscan-context.v1",
        "metadata": {
            "operation": "interproscan-context-report",
            "interproscan_version": version,
            "interproscan_tsv_layout": layout,
            "interproscan_sha256": source_sha256,
            "faa_sha256": faa_hash,
            "member_dbs": list(member_dbs),
            "rows": tsv_diag["row_count"],
            "unintegrated_rows": tsv_diag["unintegrated_count"],
            "integrated_rows": tsv_diag["integrated_count"],
            "queries": len(tsv_diag["query_metadata"]),
            "matched_queries": len(queried_ids),
            "identity_diagnostics": {
                "faa_records": len(actual_proteins),
                "gbff_records": len(base.records),
                "gbff_cds": len(cds),
                "matched_queries": len(queried_ids),
            },
            "pathway_counters": tsv_diag["pathway_counters"],
            "pathway_digest": tsv_diag["pathway_digest"],
            "support_digest": support_digest,
            "qc_flags": {
                "anti_fam": tsv_diag["anti_fam_hits"],
            },
        },
        "entries": [per_query_context[q] for q in sorted(per_query_context.keys())],
    }

    return InterProScanPlan(
        insertions=insertions,
        evidence_rows=evidence_rows,
        stats=stats,
        context_report=context_report,
        source_sha256=source_sha256,
        version=version,
        layout=layout,
        member_dbs=member_dbs,
    )


def merge(
    bakta_path: Path,
    faa_path: Path,
    interproscan_path: Path,
    output_path: Path,
    *,
    manifest_path: Path | None = None,
    interproscan_version: str | None = None,
    interproscan_tsv_layout: str = DEFAULT_TSV_LAYOUT,
    interproscan_member_dbs: str | Sequence[str] = DEFAULT_MEMBER_DBS,
    add_comment_note: bool = True,
    add_feature_provenance: bool = True,
    merge_timestamp: str | None = None,
    context_report_path: Path | None = None,
    translation_evidence_manifest: Path | None = None,
    allow_imported_translations: bool = False,
) -> dict[str, Any]:
    """Execute standalone InterProScan merge pipeline."""
    if not interproscan_version:
        raise MergeError(
            "InterProScan input requires an explicit interproscan_version assertion"
        )
    base_data = read_input_bytes(bakta_path, "Bakta")
    validate_genbank_semantics(base_data, "Bakta input")
    base = parse_genbank_bytes(base_data, "Bakta input")

    member_dbs = parse_member_dbs(interproscan_member_dbs)

    translation_evidence = (
        load_translation_evidence(translation_evidence_manifest, base)
        if translation_evidence_manifest
        else None
    )

    faa_bytes = read_input_bytes(faa_path, "FAA")
    faa_hash = sha256_bytes(faa_bytes)
    proteins = parse_faa(faa_bytes)

    plan = plan_interproscan(
        base,
        faa_path,
        interproscan_path,
        version=interproscan_version,
        layout=interproscan_tsv_layout,
        member_dbs=member_dbs,
        add_comment_note=add_comment_note,
        add_feature_provenance=add_feature_provenance,
        merge_timestamp=merge_timestamp,
        translation_evidence=translation_evidence,
        allow_imported_translations=allow_imported_translations,
        faa_data=faa_bytes,
        faa_hash=faa_hash,
        proteins=proteins,
    )

    reconciled, _rows, _stats = reconcile_insertions(plan.insertions)
    source_hashes = {"InterProScan": plan.source_sha256}

    candidate_rows, ledger_counts = build_candidate_ledger(
        base,
        plan.insertions,
        reconciled,
        plan.evidence_rows,
        source_hashes=source_hashes,
    )

    metadata: dict[str, Any] = {
        "operation": "interproscan-merge",
        "merge_timestamp": merge_timestamp or "",
        "faa_sha256": plan.stats.get("faa_sha256", faa_hash),
        "interproscan_sha256": plan.source_sha256,
        "interproscan_version": plan.version,
        "interproscan_tsv_layout": plan.layout,
        "member_dbs": list(plan.member_dbs),
        "translation_evidence_manifest": (
            str(translation_evidence_manifest) if translation_evidence_manifest else ""
        ),
        "allow_imported_translations": allow_imported_translations,
        "translation_evidence": {
            query: evidence.as_dict()
            for query, evidence in (translation_evidence or {}).items()
        },
        "translation_evidence_parent_sha256": getattr(
            translation_evidence, "manifest_sha256", ""
        ),
        "policies": {
            "add_comment_note": add_comment_note,
            "add_feature_provenance": add_feature_provenance,
        },
        "source_hashes": source_hashes,
        **ledger_counts,
    }

    other_inputs = [faa_path, interproscan_path]
    if translation_evidence_manifest:
        other_inputs.append(translation_evidence_manifest)

    return finalize_merge(
        base_path=bakta_path,
        base_data=base_data,
        output_path=output_path,
        other_inputs=other_inputs,
        insertions=reconciled,
        manifest_path=manifest_path,
        metadata=metadata,
        evidence_rows=[*candidate_rows, *plan.evidence_rows],
        sidecar_path=context_report_path,
        sidecar_payload=plan.context_report if context_report_path else None,
    )


def build_parser() -> argparse.ArgumentParser:
    """Build standalone CLI parser for InterProScan merge."""
    parser = argparse.ArgumentParser(
        prog="merge_interproscan_bakta",
        description="Enrich Bakta GenBank annotations with InterProScan evidence",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {TOOL_VERSION}"
    )
    parser.add_argument(
        "bakta_pos",
        nargs="?",
        type=Path,
        default=None,
        metavar="bakta",
        help="Path to input Bakta GenBank flatfile (.gbff)",
    )
    parser.add_argument(
        "faa_pos",
        nargs="?",
        type=Path,
        default=None,
        metavar="faa",
        help="Path to input Bakta matched protein FASTA (.faa)",
    )
    parser.add_argument(
        "interproscan_pos",
        nargs="?",
        type=Path,
        default=None,
        metavar="interproscan",
        help="Path to input InterProScan TSV output",
    )
    parser.add_argument(
        "output_pos",
        nargs="?",
        type=Path,
        default=None,
        metavar="output",
        help="Path to write enriched GenBank file",
    )
    parser.add_argument(
        "--bakta",
        type=Path,
        default=None,
        help="Path to input Bakta GenBank flatfile (.gbff)",
    )
    parser.add_argument(
        "--faa",
        type=Path,
        default=None,
        help="Path to input Bakta matched protein FASTA (.faa)",
    )
    parser.add_argument(
        "--interproscan",
        type=Path,
        default=None,
        help="Path to input InterProScan TSV output",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Path to write enriched GenBank file",
    )
    parser.add_argument(
        "--interproscan-version",
        type=str,
        default=None,
        help="Asserted InterProScan producer version (e.g. 5.59-91.0; required)",
    )
    parser.add_argument(
        "--interproscan-tsv-layout",
        type=str,
        default=DEFAULT_TSV_LAYOUT,
        choices=TSV_LAYOUT_PROFILES,
        help=f"TSV layout profile (default: {DEFAULT_TSV_LAYOUT})",
    )
    parser.add_argument(
        "--interproscan-member-dbs",
        type=str,
        default=",".join(DEFAULT_MEMBER_DBS),
        help=f"Comma-separated member databases to promote as notes (default: {','.join(DEFAULT_MEMBER_DBS)})",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="Path to write JSON merge manifest",
    )
    parser.add_argument(
        "--context-report",
        type=Path,
        default=None,
        help="Path to write JSON context sidecar",
    )
    parser.add_argument(
        "--translation-evidence-manifest",
        type=Path,
        default=None,
        help="Path to prior restoration manifest for imported translations",
    )
    parser.add_argument(
        "--allow-imported-translations",
        action="store_true",
        help="Allow functional enrichment using imported protein translations",
    )
    parser.add_argument(
        "--no-comment-note",
        action="store_true",
        help="Do not add provenance COMMENT block",
    )
    parser.add_argument(
        "--no-feature-provenance",
        action="store_true",
        help="Do not add /inference qualifiers",
    )
    parser.add_argument(
        "--merge-timestamp",
        type=str,
        default=None,
        help="Deterministic ISO timestamp for COMMENT provenance",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point for standalone InterProScan merge."""
    parser = build_parser()
    args = parser.parse_args(argv)
    bakta = args.bakta or args.bakta_pos
    faa = args.faa or args.faa_pos
    interproscan = args.interproscan or args.interproscan_pos
    output = args.output or args.output_pos
    if not (bakta and faa and interproscan and output):
        parser.error(
            "the following arguments are required: bakta, faa, interproscan, output "
            "(positional or via --bakta, --faa, --interproscan, --output)"
        )
    if not args.interproscan_version:
        parser.error("--interproscan-version is required")
    for label, path in (
        ("Bakta", bakta),
        ("FAA", faa),
        ("InterProScan", interproscan),
        ("Translation evidence", args.translation_evidence_manifest),
    ):
        if path is not None and not path.is_file():
            parser.error(f"{label} input not found: {path}")
    try:
        merge(
            bakta_path=bakta,
            faa_path=faa,
            interproscan_path=interproscan,
            output_path=output,
            manifest_path=args.manifest,
            interproscan_version=args.interproscan_version,
            interproscan_tsv_layout=args.interproscan_tsv_layout,
            interproscan_member_dbs=args.interproscan_member_dbs,
            add_comment_note=not args.no_comment_note,
            add_feature_provenance=not args.no_feature_provenance,
            merge_timestamp=args.merge_timestamp,
            context_report_path=args.context_report,
            translation_evidence_manifest=args.translation_evidence_manifest,
            allow_imported_translations=args.allow_imported_translations,
        )
    except MergeError as exc:
        print(
            f"enrich-bakta error: {exc}",
            file=sys.stderr if "sys" in globals() else None,
        )
        return 1
    return 0
