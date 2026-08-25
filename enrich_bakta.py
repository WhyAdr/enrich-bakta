#!/usr/bin/env python3
"""Canonical one-pass Baktfold, KofamScan, and eggNOG Bakta enrichment CLI."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from graft_baktfold_additions import plan_baktfold_additions
from merge_eggnog_bakta import parse_eggnog_path, plan_eggnog_additions
from merge_engine import (
    MergeError,
    finalize_merge,
    parse_faa,
    parse_genbank_bytes,
    reconcile_insertions,
    validate_genbank_semantics,
)
from merge_kofamscan_bakta import parse_kofam_table, plan_kofam_additions


def enrich(
    *,
    bakta_path: Path,
    output_path: Path,
    faa_path: Path | None = None,
    baktfold_path: Path | None = None,
    kofamscan_path: Path | None = None,
    kofamscan_version: str | None = None,
    eggnog_path: Path | None = None,
    eggnog_version: str | None = None,
    min_eggnog_confidence: str = "low",
    gene_conflict_policy: str = "skip",
    manifest_path: Path | None = None,
    add_comment_note: bool = True,
    add_feature_provenance: bool = True,
    merge_timestamp: str | None = None,
    context_report_path: Path | None = None,
) -> dict[str, Any]:
    if not any((baktfold_path, kofamscan_path, eggnog_path)):
        raise MergeError("at least one evidence source is required")
    if (kofamscan_path or eggnog_path) and faa_path is None:
        raise MergeError("--faa is required with --kofamscan or --eggnog")
    if kofamscan_version and kofamscan_path is None:
        raise MergeError("--kofamscan-version requires --kofamscan")
    if eggnog_version and eggnog_path is None:
        raise MergeError("--eggnog-version requires --eggnog")
    if context_report_path is not None and eggnog_path is None:
        raise MergeError("--context-report requires --eggnog")
    base_data = bakta_path.read_bytes()
    validate_genbank_semantics(base_data, "Bakta input")
    base = parse_genbank_bytes(base_data, "Bakta input")
    proteins = parse_faa(faa_path.read_bytes()) if faa_path else None
    insertions = []
    evidence_rows: list[dict[str, Any]] = []
    metadata: dict[str, Any] = {
        "operation": "unified-enrichment",
        "merge_timestamp": merge_timestamp or "",
    }
    other_inputs: list[Path] = []
    if faa_path:
        other_inputs.append(faa_path)
    if baktfold_path:
        data = baktfold_path.read_bytes()
        validate_genbank_semantics(data, "Baktfold input")
        source = parse_genbank_bytes(data, "Baktfold input")
        planned, stats = plan_baktfold_additions(
            base,
            source,
            baktfold_data=data,
            add_comment_note=add_comment_note,
            add_feature_provenance=add_feature_provenance,
            merge_timestamp=merge_timestamp,
            starting_order=0,
        )
        insertions.extend(planned)
        metadata["baktfold_sha256"] = stats["baktfold_sha256"]
        other_inputs.append(baktfold_path)
    if kofamscan_path:
        data = kofamscan_path.read_bytes()
        planned, rows, stats = plan_kofam_additions(
            base,
            proteins or {},
            parse_kofam_table(data),
            kofam_data=data,
            faa_data=faa_path.read_bytes() if faa_path else b"",
            kofamscan_version=kofamscan_version,
            add_comment_note=add_comment_note,
            merge_timestamp=merge_timestamp,
            starting_order=1_000_000,
        )
        insertions.extend(planned)
        evidence_rows.extend(rows)
        metadata.update(
            {
                "kofam_sha256": stats["kofam_sha256"],
                "kofamscan_version": stats["kofamscan_version"],
            }
        )
        other_inputs.append(kofamscan_path)
    if eggnog_path:
        data = eggnog_path.read_bytes()
        table = parse_eggnog_path(eggnog_path, expected_version=eggnog_version)
        planned, rows, stats = plan_eggnog_additions(
            base,
            proteins or {},
            table,
            eggnog_data=data,
            faa_data=faa_path.read_bytes() if faa_path else b"",
            min_confidence=min_eggnog_confidence,
            add_comment_note=add_comment_note,
            add_feature_provenance=add_feature_provenance,
            merge_timestamp=merge_timestamp,
            starting_order=2_000_000,
        )
        context_report = stats.pop("_context_report")
        context_report["metadata"]["operation"] = "unified-enrichment"
        insertions.extend(planned)
        evidence_rows.extend(rows)
        metadata.update(
            {
                "eggnog_sha256": stats["eggnog_sha256"],
                "eggnog_version": stats["eggnog_version"],
                "min_eggnog_confidence": min_eggnog_confidence,
            }
        )
        other_inputs.append(eggnog_path)
    insertions, reconciliation_rows, reconciliation = reconcile_insertions(
        insertions, gene_conflict_policy=gene_conflict_policy
    )
    evidence_rows.extend(reconciliation_rows)
    metadata.update(reconciliation)
    return finalize_merge(
        base_path=bakta_path,
        base_data=base_data,
        output_path=output_path,
        other_inputs=other_inputs,
        insertions=insertions,
        manifest_path=manifest_path,
        evidence_rows=evidence_rows,
        metadata=metadata,
        sidecar_path=context_report_path,
        sidecar_payload=context_report if context_report_path else None,
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Reconcile Baktfold, KofamScan, and eggNOG evidence in one byte-safe Bakta merge."
    )
    parser.add_argument("--bakta", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--faa", type=Path)
    parser.add_argument("--baktfold", type=Path)
    parser.add_argument("--kofamscan", type=Path)
    parser.add_argument("--kofamscan-version")
    parser.add_argument("--eggnog", type=Path)
    parser.add_argument("--eggnog-version")
    parser.add_argument(
        "--min-eggnog-confidence", choices=("low", "medium", "high"), default="low"
    )
    parser.add_argument(
        "--gene-conflict-policy",
        choices=("skip", "prefer-eggnog", "prefer-baktfold"),
        default="skip",
    )
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--no-comment-note", action="store_true")
    parser.add_argument("--no-feature-provenance", action="store_true")
    parser.add_argument("--merge-timestamp")
    parser.add_argument(
        "--context-report",
        type=Path,
        help="write higher-order eggNOG context as a deterministic JSON sidecar",
    )
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    for label, path in (
        ("Bakta", args.bakta),
        ("FAA", args.faa),
        ("Baktfold", args.baktfold),
        ("KofamScan", args.kofamscan),
        ("eggNOG", args.eggnog),
    ):
        if path is not None and not path.is_file():
            parser.error(f"{label} input not found: {path}")
    try:
        stats = enrich(
            bakta_path=args.bakta,
            output_path=args.output,
            faa_path=args.faa,
            baktfold_path=args.baktfold,
            kofamscan_path=args.kofamscan,
            kofamscan_version=args.kofamscan_version,
            eggnog_path=args.eggnog,
            eggnog_version=args.eggnog_version,
            min_eggnog_confidence=args.min_eggnog_confidence,
            gene_conflict_policy=args.gene_conflict_policy,
            manifest_path=args.manifest,
            add_comment_note=not args.no_comment_note,
            add_feature_provenance=not args.no_feature_provenance,
            merge_timestamp=args.merge_timestamp,
            context_report_path=args.context_report,
        )
    except MergeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(
        json.dumps(stats, indent=2, sort_keys=True)
        if args.json
        else f"Unified enrichment completed; output SHA-256: {stats['output_sha256']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
