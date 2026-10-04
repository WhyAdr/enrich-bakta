"""InterProScan evidence parser, layout validation, and candidate planner."""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from enrich_bakta_lib.core.merge_engine import (
    Insertion,
    MergeError,
)
from enrich_bakta_lib.sources.value_rules import (
    validate_go_term,
    validate_interpro_accession,
)

SUPPORTED_INTERPROSCAN_VERSIONS = frozenset({"5.59-91.0"})
DEFAULT_INTERPROSCAN_VERSION = "5.59-91.0"

TSV_LAYOUT_PROFILES = ("ipr-go-pathways", "ipr-go", "ipr-pathways", "ipr-only")
DEFAULT_TSV_LAYOUT = "ipr-go-pathways"

SUPPORTED_MEMBER_DBS = frozenset({"Pfam", "TIGRFAM"})
DEFAULT_MEMBER_DBS = ("Pfam", "TIGRFAM")

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
