"""Check GenBank db_xref syntax and namespaces against an offline INSDC snapshot."""

from __future__ import annotations

import argparse
import csv
import io
import sys
from collections import Counter
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from enrich_bakta_lib.core.artifacts import (
    json_sidecar_bytes,
    preflight_artifact_paths,
    stage_artifacts,
)
from enrich_bakta_lib.core.genbank import feature_uid, parse_genbank_bytes
from enrich_bakta_lib.core.types import (
    TOOL_VERSION,
    MergeError,
    read_input_bytes,
    sha256_bytes,
)
from enrich_bakta_lib.workflows.insdc_db_xref import (
    DATABASE_PREFIXES,
    FEATURE_TABLE_URL,
    FEATURE_TABLE_VERSION,
    NCBI_DB_XREF_URL,
    VOCABULARY_DATE,
    VOCABULARY_ID,
    VOCABULARY_URL,
)

SCHEMA = "enrich-bakta.submission-preflight.v1"
_CANONICAL_BY_CASE = {prefix.casefold(): prefix for prefix in DATABASE_PREFIXES}
_COLUMNS = (
    "record_id",
    "record_index",
    "locus_tag",
    "feature_type",
    "feature_uid",
    "feature_location",
    "feature_byte_offset",
    "qualifier",
    "qualifier_index",
    "value",
    "database_prefix",
    "code",
    "suggested_prefix",
)


def _db_xref_finding(value: str) -> tuple[str, str, str] | None:
    prefix, separator, identifier = value.partition(":")
    if (
        not separator
        or not prefix
        or not identifier
        or prefix != prefix.strip()
        or identifier != identifier.strip()
        or not value.isprintable()
    ):
        return "malformed_db_xref", prefix if separator else "", ""
    if prefix in DATABASE_PREFIXES:
        return None
    canonical = _CANONICAL_BY_CASE.get(prefix.casefold())
    if canonical is not None:
        return "database_prefix_case", prefix, canonical
    return "unknown_database_prefix", prefix, ""


def check_submission_gbff(gbff_path: Path) -> dict[str, Any]:
    data = read_input_bytes(gbff_path, "submission GBFF")
    document = parse_genbank_bytes(data, "submission GBFF")
    findings = []
    counts: Counter[str] = Counter()
    feature_count = 0
    db_xref_count = 0
    for record in document.records:
        feature_count += len(record.features)
        for feature in record.features:
            db_xrefs = feature.values("db_xref")
            db_xref_count += len(db_xrefs)
            for index, value in enumerate(db_xrefs):
                finding = _db_xref_finding(value)
                if finding is None:
                    continue
                code, prefix, canonical = finding
                counts[code] += 1
                findings.append(
                    {
                        "record_id": record.record_id,
                        "record_index": record.index,
                        "locus_tag": feature.locus_tag or "",
                        "feature_type": feature.feature_type,
                        "feature_uid": feature_uid(feature),
                        "feature_location": feature.location,
                        "feature_byte_offset": feature.start,
                        "qualifier": "db_xref",
                        "qualifier_index": index,
                        "value": value,
                        "database_prefix": prefix,
                        "code": code,
                        "suggested_prefix": canonical,
                    }
                )
    return {
        "schema": SCHEMA,
        "scope": "db_xref value syntax and database-prefix vocabulary",
        "input_sha256": sha256_bytes(data),
        "vocabulary": {
            "id": VOCABULARY_ID,
            "retrieved_on": VOCABULARY_DATE,
            "source_url": VOCABULARY_URL,
            "feature_table_version": FEATURE_TABLE_VERSION,
            "feature_table_url": FEATURE_TABLE_URL,
            "ncbi_reference_url": NCBI_DB_XREF_URL,
            "prefix_count": len(DATABASE_PREFIXES),
        },
        "record_count": len(document.records),
        "feature_count": feature_count,
        "db_xref_count": db_xref_count,
        "finding_count": len(findings),
        "findings_by_code": dict(counts),
        "findings": findings,
    }


def _tsv_bytes(report: dict[str, Any]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow((*_COLUMNS, "input_sha256", "vocabulary_id"))
    for finding in report["findings"]:
        writer.writerow(
            (
                *(finding[column] for column in _COLUMNS),
                report["input_sha256"],
                report["vocabulary"]["id"],
            )
        )
    return stream.getvalue().encode("utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {TOOL_VERSION}"
    )
    parser.add_argument(
        "--gbff", type=Path, required=True, help="GenBank file to inspect"
    )
    parser.add_argument("--format", choices=("json", "tsv"), default="json")
    parser.add_argument("--output", type=Path, help="report path (default: stdout)")
    args = parser.parse_args(argv)
    try:
        if args.output is not None:
            preflight_artifact_paths([args.output], [args.gbff])
        report = check_submission_gbff(args.gbff)
        rendered = (
            _tsv_bytes(report) if args.format == "tsv" else json_sidecar_bytes(report)
        )
        if args.output is None:
            sys.stdout.buffer.write(rendered)
        else:
            stage_artifacts([(args.output, rendered)], [args.gbff])
    except (MergeError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return int(bool(report["finding_count"]))


if __name__ == "__main__":
    raise SystemExit(main())
