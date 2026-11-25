"""InterProScan evidence parser, layout validation, and candidate planner."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from enrich_bakta_lib.core.merge_engine import (
    TOOL_VERSION,
    Insertion,
    MergeError,
    RawDocument,
    TranslationEvidence,
    cds_by_locus,
    comment_insertion,
    feature_uid,
    has_translation_evidence_marker,
    normalize_protein,
    parse_faa,
    protein_sha256,
    qualifier_insertion,
    read_input_bytes,
    record_has_legacy_restoration_comment,
    sha256_bytes,
)
from enrich_bakta_lib.sources.base import MergeContext, SourcePlan
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
MAX_LOCATION_SAMPLES = 5
MAX_CONTEXT_SIGNATURES = 50
MAX_PATHWAY_DB_KEYS = 200

MD5_HEX_RE = re.compile(r"^[0-9a-fA-F]{32}\Z", re.ASCII)
DATE_FORMAT_RE = re.compile(
    r"^(?:\d{2}-\d{2}-\d{4}|\d{4}-\d{2}-\d{2}|-|\?)\Z", re.ASCII
)
ASCII_INT_RE = re.compile(r"^\d+\Z", re.ASCII)


def _is_valid_score(value: str) -> bool:
    """Validate native InterProScan score token ('-' or finite ASCII float)."""
    if value == "-":
        return True
    if not value or not value.isascii():
        return False
    try:
        parsed = float(value)
    except ValueError:
        return False
    return math.isfinite(parsed)


@dataclass
class InterProScanHit:
    """Parsed, validated, and aggregated hit from an InterProScan TSV row."""

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
    occurrence_count: int = 1
    location_samples: tuple[tuple[int, int], ...] = ()
    location_samples_truncated: bool = False


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
    hits_by_key: dict[tuple[str, str, str, str | None], InterProScanHit] = {}
    query_metadata: dict[str, tuple[str, int]] = {}
    pathway_counters: dict[str, Any] = {
        "total_occurrences": 0,
        "by_database": defaultdict(int),
    }
    pathway_occurrences_by_query: dict[str, int] = defaultdict(int)
    pathway_occurrences_by_query_db: dict[str, dict[str, int]] = defaultdict(
        lambda: defaultdict(int)
    )
    anti_fam_hits: list[dict[str, Any]] = []
    integrated_count = 0
    unintegrated_count = 0
    row_number = 0

    with tsv_path.open("rb") as stream:
        while True:
            line_bytes = stream.readline(MAX_RAW_LINE_BYTES + 2)
            if not line_bytes:
                break
            if line_bytes.endswith(b"\r\n"):
                raw_content_len = len(line_bytes) - 2
            elif line_bytes.endswith((b"\n", b"\r")):
                raw_content_len = len(line_bytes) - 1
            else:
                raw_content_len = len(line_bytes)
            if raw_content_len > MAX_RAW_LINE_BYTES:
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

            # Strip only line endings, preserving empty fields
            if line_text.endswith("\r\n"):
                line_text = line_text[:-2]
            elif line_text.endswith("\n") or line_text.endswith("\r"):
                line_text = line_text[:-1]

            if "\x00" in line_text:
                raise MergeError(
                    f"InterProScan TSV line {row_number} contains embedded NUL character"
                )

            if any((ord(c) < 32 and c != "\t") or ord(c) == 127 for c in line_text):
                raise MergeError(
                    f"InterProScan TSV line {row_number} contains invalid control character"
                )

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

            raw_len = fields[2].strip()
            if not ASCII_INT_RE.fullmatch(raw_len):
                raise MergeError(
                    f"InterProScan TSV line {row_number} has non-integer length {fields[2]!r}"
                )
            length = int(raw_len)
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

            raw_start = fields[6].strip()
            raw_end = fields[7].strip()
            if not ASCII_INT_RE.fullmatch(raw_start) or not ASCII_INT_RE.fullmatch(
                raw_end
            ):
                raise MergeError(
                    f"InterProScan TSV line {row_number} has non-integer coordinates"
                )
            start = int(raw_start)
            end = int(raw_end)

            if not (1 <= start <= end <= length):
                raise MergeError(
                    f"InterProScan TSV line {row_number} has invalid coordinates {start}..{end} "
                    f"for protein length {length}"
                )

            score = fields[8].strip()
            status = fields[9].strip()
            date = fields[10].strip()

            if not _is_valid_score(score):
                raise MergeError(
                    f"InterProScan TSV line {row_number} has invalid score {score!r}"
                )

            if status not in {"T", "F", "?", "-"}:
                raise MergeError(
                    f"InterProScan TSV line {row_number} has invalid status {status!r}"
                )

            if not DATE_FORMAT_RE.fullmatch(date):
                raise MergeError(
                    f"InterProScan TSV line {row_number} has invalid date {date!r}"
                )

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
                pathway_occurrences_by_query[query_id] += 1
                db_name = (
                    token.partition(":")[0].strip() if ":" in token else "unspecified"
                )
                if (
                    len(pathway_counters["by_database"]) < MAX_PATHWAY_DB_KEYS
                    or db_name in pathway_counters["by_database"]
                ):
                    pathway_counters["by_database"][db_name] += 1
                else:
                    pathway_counters["by_database"]["other"] += 1
                q_dbs = pathway_occurrences_by_query_db[query_id]
                if len(q_dbs) < MAX_PATHWAY_DB_KEYS or db_name in q_dbs:
                    q_dbs[db_name] += 1
                else:
                    q_dbs["other"] += 1
                pathway_hasher.update(f"{query_id}\t{token}\n".encode("utf-8"))

            # QC flags (AntiFam)
            if analysis == "AntiFam" or sig_acc.startswith("ANF"):
                if len(anti_fam_hits) < 1000:
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

            hit_key = (query_id, analysis, sig_acc, ipr_acc)
            if hit_key in hits_by_key:
                existing_hit = hits_by_key[hit_key]
                existing_hit.occurrence_count += 1
                if len(existing_hit.location_samples) < MAX_LOCATION_SAMPLES:
                    existing_hit.location_samples = existing_hit.location_samples + (
                        (start, end),
                    )
                else:
                    existing_hit.location_samples_truncated = True
                if go_terms:
                    existing_gos = set(existing_hit.go_terms)
                    new_gos = [g for g in go_terms if g not in existing_gos]
                    if new_gos:
                        existing_hit.go_terms = existing_hit.go_terms + tuple(new_gos)
                if pathways and len(existing_hit.pathways) < 25:
                    existing_pws = set(existing_hit.pathways)
                    new_pws = [p for p in pathways if p not in existing_pws]
                    if new_pws:
                        existing_hit.pathways = (
                            existing_hit.pathways + tuple(new_pws)
                        )[:25]
            else:
                if len(hits_by_key) >= MAX_AGGREGATED_CANDIDATES:
                    raise MergeError(
                        f"InterProScan TSV exceeded maximum unique candidate aggregation limit ({MAX_AGGREGATED_CANDIDATES})"
                    )
                hits_by_key[hit_key] = InterProScanHit(
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
                    pathways=tuple(pathways[:25]),
                    row_number=row_number,
                    occurrence_count=1,
                    location_samples=((start, end),),
                    location_samples_truncated=False,
                )

    hits = list(hits_by_key.values())
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
        "pathway_occurrences_by_query": dict(pathway_occurrences_by_query),
        "pathway_occurrences_by_query_db": {
            q: dict(sorted(dbs.items()))
            for q, dbs in pathway_occurrences_by_query_db.items()
        },
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
    starting_order: int = 0,
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
    order = starting_order
    support_hasher = hashlib.sha256()
    match_evidence_hasher = hashlib.sha256()

    new_qualifier_features: set[str] = set()
    per_query_context: dict[str, dict[str, Any]] = {}
    unique_candidate_total = 0

    for query in queried_ids:
        feature = cds[query]
        f_uid = feature_uid(feature)
        query_hits = hits_by_query.get(query, [])

        # Aggregated proposals per (field, qualifier, normalized_value)
        unique_proposals: dict[
            tuple[str, str, str],
            dict[str, Any],
        ] = {}
        context_sigs: list[str] = []

        def _record_proposal(
            field: str,
            qualifier: str,
            normalized_val: str,
            raw_val: str,
            role: str,
            ev_class: str,
            hit: InterProScanHit,
        ) -> None:
            key = (field, qualifier, normalized_val)
            entry = unique_proposals.get(key)
            if entry is None:
                unique_proposals[key] = {
                    "field": field,
                    "qualifier": qualifier,
                    "normalized_value": normalized_val,
                    "raw_value": raw_val,
                    "role": role,
                    "evidence_class": ev_class,
                    "first_row": hit.row_number,
                    "hits": [hit],
                }
            else:
                if hit.row_number < entry["first_row"]:
                    entry["first_row"] = hit.row_number
                entry["hits"].append(hit)

        for hit in query_hits:
            # 1. InterPro accession
            if hit.interpro_accession:
                _record_proposal(
                    "InterPro",
                    "db_xref",
                    f"InterPro:{hit.interpro_accession}",
                    hit.interpro_accession,
                    "functional_proposal",
                    "InterPro",
                    hit,
                )

            # 2. GO terms
            for go in hit.go_terms:
                _record_proposal(
                    "GO",
                    "db_xref",
                    go,
                    go,
                    "functional_proposal",
                    "GO",
                    hit,
                )

            # 3. Member signatures (Pfam, TIGRFAM)
            if hit.analysis == "Pfam":
                if "Pfam" in member_dbs:
                    res_pfam = validate_pfam_accession(hit.signature_accession)
                    if not res_pfam:
                        raise MergeError(
                            f"InterProScan TSV row {hit.row_number}: [invalid_pfam_accession] "
                            f"invalid Pfam signature accession {hit.signature_accession!r}"
                        )
                    canonical, _ = res_pfam
                    _record_proposal(
                        "Pfam",
                        "note",
                        f"PFAM:{canonical}",
                        hit.signature_accession,
                        "substantive_evidence",
                        "member:Pfam",
                        hit,
                    )
                else:
                    context_sigs.append(f"{hit.analysis}:{hit.signature_accession}")
            elif hit.analysis == "TIGRFAM":
                if "TIGRFAM" in member_dbs:
                    res_tigr = validate_tigrfam_accession(hit.signature_accession)
                    if not res_tigr:
                        raise MergeError(
                            f"InterProScan TSV row {hit.row_number}: [invalid_tigrfam_accession] "
                            f"invalid TIGRFAM signature accession {hit.signature_accession!r}"
                        )
                    canonical, _ = res_tigr
                    _record_proposal(
                        "TIGRFAM",
                        "note",
                        f"TIGRFAM:{canonical}",
                        hit.signature_accession,
                        "substantive_evidence",
                        "member:TIGRFAM",
                        hit,
                    )
                else:
                    context_sigs.append(f"{hit.analysis}:{hit.signature_accession}")
            else:
                context_sigs.append(f"{hit.analysis}:{hit.signature_accession}")

        unique_candidate_total += len(unique_proposals)
        if unique_candidate_total > MAX_AGGREGATED_CANDIDATES:
            raise MergeError(
                f"InterProScan planning exceeded maximum unique semantic candidate limit ({MAX_AGGREGATED_CANDIDATES})"
            )

        for prop in sorted(
            unique_proposals.values(),
            key=lambda p: (p["qualifier"], p["normalized_value"]),
        ):
            field = prop["field"]
            qualifier = prop["qualifier"]
            normalized_val = prop["normalized_value"]
            raw_val = prop["raw_value"]
            role = prop["role"]
            ev_class = prop["evidence_class"]
            first_row = prop["first_row"]
            contrib_hits: list[InterProScanHit] = prop["hits"]

            # Bounded match-evidence summary and deterministic digest
            occurrence_count = sum(h.occurrence_count for h in contrib_hits)
            location_samples: list[list[int]] = []
            location_samples_truncated = False
            cand_match_hasher = hashlib.sha256()
            for h in sorted(
                contrib_hits,
                key=lambda item: (
                    item.analysis,
                    item.signature_accession,
                    item.row_number,
                ),
            ):
                if h.location_samples_truncated:
                    location_samples_truncated = True
                for s, e in h.location_samples:
                    if len(location_samples) < MAX_LOCATION_SAMPLES:
                        location_samples.append([s, e])
                    else:
                        location_samples_truncated = True
                loc_str = ",".join(f"{s}-{e}" for s, e in h.location_samples)
                cand_match_hasher.update(
                    f"{query}\t{h.analysis}\t{h.signature_accession}\t{h.row_number}\t{h.occurrence_count}\t{loc_str}\n".encode(
                        "utf-8"
                    )
                )
            candidate_match_digest = cand_match_hasher.hexdigest()

            # Support digest (annotation set) vs match evidence digest
            support_hasher.update(
                f"{query}\t{field}\t{normalized_val}\n".encode("utf-8")
            )
            match_evidence_hasher.update(
                f"{query}\t{field}\t{normalized_val}\t{candidate_match_digest}\n".encode(
                    "utf-8"
                )
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
                    "occurrence_count": occurrence_count,
                    "location_samples": location_samples,
                    "location_samples_truncated": location_samples_truncated,
                    "match_evidence_digest": candidate_match_digest,
                    "source": "InterProScan",
                    "reason": f"{field} match {'already present' if existing else 'planned for insertion'}",
                }
            )

        pathway_occurrences = tsv_diag.get("pathway_occurrences_by_query", {}).get(
            query, 0
        )
        pathway_by_db = tsv_diag.get("pathway_occurrences_by_query_db", {}).get(
            query, {}
        )
        unique_context_sigs = sorted(set(context_sigs))
        if unique_context_sigs or pathway_occurrences > 0:
            per_query_context[query] = {
                "query_id": query,
                "context_signatures": unique_context_sigs[:MAX_CONTEXT_SIGNATURES],
                "context_signatures_total": len(unique_context_sigs),
                "context_signatures_truncated": len(unique_context_sigs)
                > MAX_CONTEXT_SIGNATURES,
                "pathway_occurrences": pathway_occurrences,
                "pathway_occurrences_by_database": pathway_by_db,
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
                        "support_witness": None,
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
    match_evidence_digest = match_evidence_hasher.hexdigest()

    stats = {
        "rows": tsv_diag["row_count"],
        "hits": tsv_diag["hits_count"],
        "queries": len(tsv_diag["query_metadata"]),
        "matched_queries": len(queried_ids),
        "planned_insertions": len(insertions),
        "interproscan_sha256": source_sha256,
        "faa_sha256": actual_faa_hash,
    }

    context_report = {
        "schema": "enrich-bakta.interproscan-context.v1",
        "metadata": {
            "operation": "interproscan-context-report",
            "interproscan_version": version,
            "interproscan_tsv_layout": layout,
            "interproscan_sha256": source_sha256,
            "faa_sha256": actual_faa_hash,
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
            "match_evidence_digest": match_evidence_digest,
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


@dataclass(frozen=True)
class InterProScanAdapter:
    """Adapter wrapping InterProScan TSV validation and candidate planning."""

    interproscan_path: Path
    version: str
    layout: str = DEFAULT_TSV_LAYOUT
    member_dbs: str | Sequence[str] = DEFAULT_MEMBER_DBS
    source_id: str = "InterProScan"

    def plan(self, ctx: MergeContext, *, starting_order: int = 0) -> SourcePlan:
        if ctx.faa_path is None:
            raise MergeError("--faa is required with --interproscan")
        # Preserve 0.4.1 insertion order=0 so multi-source GBFF splice ordering
        # remains byte-identical across baselines.
        _ = starting_order
        parsed_dbs = parse_member_dbs(self.member_dbs)
        ips_plan = plan_interproscan(
            ctx.base_doc,
            ctx.faa_path,
            self.interproscan_path,
            version=self.version,
            layout=self.layout,
            member_dbs=parsed_dbs,
            add_comment_note=ctx.add_comment_note,
            add_feature_provenance=ctx.add_feature_provenance,
            merge_timestamp=ctx.merge_timestamp,
            translation_evidence=ctx.translation_evidence,
            allow_imported_translations=ctx.allow_imported_translations,
            faa_data=ctx.faa_data,
            faa_hash=ctx.faa_sha256,
            proteins=ctx.proteins,
            starting_order=0,
        )
        return SourcePlan(
            source_id=self.source_id,
            source_sha256=ips_plan.source_sha256,
            insertions=ips_plan.insertions,
            evidence_rows=ips_plan.evidence_rows,
            stats=ips_plan.stats,
            manifest_metadata={
                "interproscan_sha256": ips_plan.source_sha256,
                "interproscan_version": ips_plan.version,
                "interproscan_tsv_layout": ips_plan.layout,
                "interproscan_member_dbs": list(ips_plan.member_dbs),
            },
            context_report=ips_plan.context_report,
            context_key="interproscan",
            input_paths=(self.interproscan_path,),
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
    from enrich_bakta_lib.workflows.pipeline import (
        capture_merge_context,
        execute_merge_pipeline,
    )

    if not interproscan_version:
        raise MergeError(
            "InterProScan input requires an explicit interproscan_version assertion"
        )
    ctx = capture_merge_context(
        bakta_path=bakta_path,
        faa_path=faa_path,
        translation_evidence_manifest=translation_evidence_manifest,
        allow_imported_translations=allow_imported_translations,
        add_comment_note=add_comment_note,
        add_feature_provenance=add_feature_provenance,
        merge_timestamp=merge_timestamp,
        require_manifest_message=(
            "restored input requires --translation-evidence-manifest and preserved translation evidence"
        ),
    )

    def _ips_metadata(
        meta: dict[str, Any],
        plans: list[SourcePlan],
        _candidate_counts: dict[str, Any],
    ) -> dict[str, Any]:
        plan = plans[0]
        return {
            **meta,
            "member_dbs": list(plan.manifest_metadata["interproscan_member_dbs"]),
            "source_hashes": {plan.source_id: plan.source_sha256},
        }

    result = execute_merge_pipeline(
        ctx,
        [
            InterProScanAdapter(
                interproscan_path=interproscan_path,
                version=interproscan_version,
                layout=interproscan_tsv_layout,
                member_dbs=interproscan_member_dbs,
            )
        ],
        operation="interproscan-merge",
        output_path=output_path,
        manifest_path=manifest_path,
        context_report_path=context_report_path,
        reconcile=True,
        metadata_override=_ips_metadata,
    )
    return result.final


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
