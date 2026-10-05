#!/usr/bin/env python3
"""Canonical one-pass Baktfold, KofamScan, and eggNOG Bakta enrichment CLI."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from enrich_bakta_lib.core.decisions import (
    build_candidate_ledger,
    validate_candidate_ledger,
)
from enrich_bakta_lib.core.merge_engine import (
    TOOL_VERSION,
    MergeError,
    finalize_merge,
    has_translation_evidence_marker,
    load_translation_evidence,
    parse_faa,
    parse_genbank_bytes,
    read_input_bytes,
    reconcile_insertions,
    sha256_bytes,
    validate_genbank_semantics,
)
from enrich_bakta_lib.sources.baktfold import plan_baktfold_additions
from enrich_bakta_lib.sources.eggnog import parse_eggnog_path, plan_eggnog_additions
from enrich_bakta_lib.sources.interproscan import (
    DEFAULT_MEMBER_DBS,
    DEFAULT_TSV_LAYOUT,
    TSV_LAYOUT_PROFILES,
    parse_member_dbs,
    plan_interproscan,
)
from enrich_bakta_lib.sources.kofam import parse_kofam_table, plan_kofam_additions


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
    eggnog_schema: str | None = None,
    interproscan_path: Path | None = None,
    interproscan_version: str | None = None,
    interproscan_tsv_layout: str | None = None,
    interproscan_member_dbs: str | Sequence[str] | None = None,
    min_eggnog_confidence: str = "low",
    gene_conflict_policy: str = "skip",
    manifest_path: Path | None = None,
    add_comment_note: bool = True,
    add_feature_provenance: bool = True,
    merge_timestamp: str | None = None,
    context_report_path: Path | None = None,
    clean_gene_suffix: bool = False,
    translation_evidence_manifest: Path | None = None,
    allow_imported_translations: bool = False,
    baktfold_invalid_ec_policy: str = "reject",
) -> dict[str, Any]:
    if not any((baktfold_path, kofamscan_path, eggnog_path, interproscan_path)):
        raise MergeError("at least one evidence source is required")
    if (kofamscan_path or eggnog_path or interproscan_path) and faa_path is None:
        raise MergeError(
            "--faa is required with --kofamscan, --eggnog, or --interproscan"
        )
    if kofamscan_version and kofamscan_path is None:
        raise MergeError("--kofamscan-version requires --kofamscan")
    if eggnog_version and eggnog_path is None:
        raise MergeError("--eggnog-version requires --eggnog")
    if eggnog_schema and eggnog_path is None:
        raise MergeError("--eggnog-schema requires --eggnog")
    if interproscan_path is not None and not interproscan_version:
        raise MergeError(
            "InterProScan input requires an explicit --interproscan-version assertion"
        )
    if interproscan_version and interproscan_path is None:
        raise MergeError("--interproscan-version requires --interproscan")
    if interproscan_tsv_layout and interproscan_path is None:
        raise MergeError("--interproscan-tsv-layout requires --interproscan")
    if interproscan_member_dbs is not None and interproscan_path is None:
        raise MergeError("--interproscan-member-dbs requires --interproscan")
    if (
        context_report_path is not None
        and eggnog_path is None
        and interproscan_path is None
    ):
        raise MergeError("--context-report requires --eggnog or --interproscan")
    base_data = read_input_bytes(bakta_path, "Bakta")
    validate_genbank_semantics(base_data, "Bakta input")
    base = parse_genbank_bytes(base_data, "Bakta input")
    if (
        has_translation_evidence_marker(base_data)
        and translation_evidence_manifest is None
    ):
        raise MergeError(
            "restored input requires --translation-evidence-manifest for functional validation"
        )
    translation_evidence = (
        load_translation_evidence(translation_evidence_manifest, base)
        if translation_evidence_manifest is not None
        else None
    )
    faa_data = read_input_bytes(faa_path, "FAA") if faa_path else b""
    proteins = parse_faa(faa_data) if faa_path else None
    insertions = []
    evidence_rows: list[dict[str, Any]] = []
    metadata: dict[str, Any] = {
        "operation": "unified-enrichment",
        "merge_timestamp": merge_timestamp or "",
        "faa_sha256": sha256_bytes(faa_data) if faa_path else "",
        "translation_evidence": {
            query: evidence.as_dict()
            for query, evidence in (translation_evidence or {}).items()
        },
        "translation_evidence_parent_sha256": getattr(
            translation_evidence, "manifest_sha256", ""
        ),
        "translation_evidence_manifest": str(translation_evidence_manifest)
        if translation_evidence_manifest
        else "",
        "allow_imported_translations": allow_imported_translations,
        "policies": {
            "add_comment_note": add_comment_note,
            "add_feature_provenance": add_feature_provenance,
            "clean_gene_suffix": clean_gene_suffix,
        },
    }
    other_inputs: list[Path] = []
    if faa_path:
        other_inputs.append(faa_path)
    if translation_evidence_manifest is not None:
        other_inputs.append(translation_evidence_manifest)
    if baktfold_path:
        data = read_input_bytes(baktfold_path, "Baktfold")
        validate_genbank_semantics(data, "Baktfold input")
        source = parse_genbank_bytes(data, "Baktfold input")
        planned, baktfold_candidates, stats = plan_baktfold_additions(
            base,
            source,
            baktfold_data=data,
            add_comment_note=add_comment_note,
            add_feature_provenance=add_feature_provenance,
            merge_timestamp=merge_timestamp,
            starting_order=0,
            invalid_ec_policy=baktfold_invalid_ec_policy,
            translation_evidence=translation_evidence,
            allow_imported_translations=allow_imported_translations,
        )
        insertions.extend(planned)
        evidence_rows.extend(baktfold_candidates)
        evidence_rows.extend(stats.get("invalid_ec_values", []))
        metadata["baktfold_sha256"] = stats["baktfold_sha256"]
        metadata["baktfold_parity"] = stats["parity"]
        metadata["baktfold_version"] = stats["baktfold_version"]
        metadata["baktfold_version_detected"] = stats["baktfold_version_detected"]
        metadata["baktfold_translation_mismatch_policy"] = stats[
            "translation_mismatch_policy"
        ]
        metadata["baktfold_translation_mismatch_loci"] = stats[
            "translation_mismatch_loci"
        ]
        metadata["baktfold_translation_mismatch_suppressed_features"] = stats[
            "translation_mismatch_suppressed_features"
        ]
        metadata["baktfold_base_pair_name_conflicts"] = stats[
            "base_pair_name_conflicts"
        ]
        metadata["baktfold_source_pair_name_conflicts"] = stats[
            "source_pair_name_conflicts"
        ]
        metadata["baktfold_invalid_ec_policy"] = baktfold_invalid_ec_policy
        metadata["baktfold_invalid_ec_values"] = stats.get("invalid_ec_values", [])
        metadata["baktfold_invalid_ec_value_count"] = stats.get(
            "invalid_ec_value_count", 0
        )
        other_inputs.append(baktfold_path)
    if kofamscan_path:
        data = read_input_bytes(kofamscan_path, "KofamScan")
        planned, rows, stats = plan_kofam_additions(
            base,
            proteins or {},
            parse_kofam_table(data),
            kofam_data=data,
            faa_data=faa_data,
            kofamscan_version=kofamscan_version,
            add_comment_note=add_comment_note,
            add_feature_provenance=add_feature_provenance,
            merge_timestamp=merge_timestamp,
            starting_order=max(
                (insertion.order for insertion in insertions), default=-1
            )
            + 1,
            translation_evidence=translation_evidence,
            allow_imported_translations=allow_imported_translations,
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
    eggnog_context: dict[str, Any] | None = None
    interproscan_context: dict[str, Any] | None = None
    if eggnog_path:
        data = read_input_bytes(eggnog_path, "eggNOG")
        table = parse_eggnog_path(
            eggnog_path,
            expected_version=eggnog_version,
            data=data,
            schema_id=eggnog_schema,
        )
        plan = plan_eggnog_additions(
            base,
            proteins or {},
            table,
            eggnog_data=data,
            faa_data=faa_data,
            min_confidence=min_eggnog_confidence,
            clean_gene_suffix=clean_gene_suffix,
            add_comment_note=add_comment_note,
            add_feature_provenance=add_feature_provenance,
            merge_timestamp=merge_timestamp,
            starting_order=max(
                (insertion.order for insertion in insertions), default=-1
            )
            + 1,
            translation_evidence=translation_evidence,
            allow_imported_translations=allow_imported_translations,
        )
        planned = plan.insertions
        rows = plan.evidence_rows
        stats = plan.stats
        eggnog_context = plan.context_report
        eggnog_context["metadata"]["operation"] = "unified-enrichment"
        insertions.extend(planned)
        evidence_rows.extend(rows)
        metadata.update(
            {
                "eggnog_sha256": stats["eggnog_sha256"],
                "eggnog_version": stats["eggnog_version"],
                "eggnog_schema": stats["eggnog_schema"],
                "schema_selection_source": stats["schema_selection_source"],
                "min_eggnog_confidence": min_eggnog_confidence,
                "confidence_field_order": stats["confidence_field_order"],
                "confidence_contract_source": stats["confidence_contract_source"],
                "clean_gene_suffix": clean_gene_suffix,
                "translation_evidence_manifest": (
                    str(translation_evidence_manifest)
                    if translation_evidence_manifest
                    else ""
                ),
                "allow_imported_translations": allow_imported_translations,
            }
        )
        other_inputs.append(eggnog_path)
    if interproscan_path:
        assert faa_path is not None
        assert interproscan_version is not None
        member_dbs = parse_member_dbs(
            interproscan_member_dbs
            if interproscan_member_dbs is not None
            else DEFAULT_MEMBER_DBS
        )
        plan_ips = plan_interproscan(
            base,
            faa_path,
            interproscan_path,
            version=interproscan_version,
            layout=interproscan_tsv_layout or DEFAULT_TSV_LAYOUT,
            member_dbs=member_dbs,
            add_comment_note=add_comment_note,
            add_feature_provenance=add_feature_provenance,
            merge_timestamp=merge_timestamp,
            translation_evidence=translation_evidence,
            allow_imported_translations=allow_imported_translations,
            faa_data=faa_data,
            faa_hash=metadata["faa_sha256"],
            proteins=proteins,
        )
        insertions.extend(plan_ips.insertions)
        evidence_rows.extend(plan_ips.evidence_rows)
        interproscan_context = plan_ips.context_report
        interproscan_context["metadata"]["operation"] = "unified-enrichment"
        metadata.update(
            {
                "interproscan_sha256": plan_ips.source_sha256,
                "interproscan_version": plan_ips.version,
                "interproscan_tsv_layout": plan_ips.layout,
                "interproscan_member_dbs": list(plan_ips.member_dbs),
            }
        )
        other_inputs.append(interproscan_path)
    planned_insertions = list(insertions)
    insertions, reconciliation_rows, reconciliation = reconcile_insertions(
        insertions, gene_conflict_policy=gene_conflict_policy
    )
    evidence_rows.extend(reconciliation_rows)
    metadata.update(reconciliation)
    candidate_rows, candidate_counts = build_candidate_ledger(
        base,
        planned_insertions,
        insertions,
        evidence_rows,
        source_hashes={
            "Baktfold": str(metadata.get("baktfold_sha256", "")),
            "KofamScan": str(metadata.get("kofam_sha256", "")),
            "eggNOG": str(metadata.get("eggnog_sha256", "")),
            "InterProScan": str(metadata.get("interproscan_sha256", "")),
        },
        gene_conflict_policy=gene_conflict_policy,
    )
    evidence_rows.extend(candidate_rows)
    validate_candidate_ledger(candidate_rows, insertions, base=base)
    metadata.update(candidate_counts)
    sidecar_payload: dict[str, Any] | None = None
    if context_report_path:
        if eggnog_context is not None and interproscan_context is not None:
            sidecar_payload = {
                "schema": "enrich-bakta.context.v1",
                "metadata": {
                    "operation": "unified-enrichment",
                },
                "eggnog": eggnog_context,
                "interproscan": interproscan_context,
            }
        elif eggnog_context is not None:
            sidecar_payload = eggnog_context
        elif interproscan_context is not None:
            sidecar_payload = interproscan_context
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
        sidecar_payload=sidecar_payload,
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Reconcile Baktfold, KofamScan, eggNOG, and InterProScan evidence in one byte-safe Bakta merge."
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {TOOL_VERSION}"
    )
    parser.add_argument("--bakta", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--faa", type=Path)
    parser.add_argument("--baktfold", type=Path)
    parser.add_argument("--kofamscan", type=Path)
    parser.add_argument("--kofamscan-version")
    parser.add_argument("--eggnog", type=Path)
    parser.add_argument("--eggnog-version")
    parser.add_argument("--eggnog-schema")
    parser.add_argument(
        "--min-eggnog-confidence", choices=("low", "medium", "high"), default="low"
    )
    parser.add_argument("--clean-gene-suffix", action="store_true")
    parser.add_argument(
        "--gene-conflict-policy",
        choices=("skip", "prefer-eggnog", "prefer-baktfold"),
        default="skip",
    )
    parser.add_argument(
        "--interproscan", type=Path, help="InterProScan TSV evidence file"
    )
    parser.add_argument(
        "--interproscan-version",
        default=None,
        help="asserted InterProScan producer version (e.g. 5.59-91.0; required when --interproscan is provided)",
    )
    parser.add_argument(
        "--interproscan-tsv-layout",
        choices=TSV_LAYOUT_PROFILES,
        default=None,
        help=f"InterProScan TSV layout profile (default: {DEFAULT_TSV_LAYOUT})",
    )
    parser.add_argument(
        "--interproscan-member-dbs",
        default=None,
        help="comma-separated member databases to promote as notes (Pfam,TIGRFAM)",
    )
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--no-comment-note", action="store_true")
    parser.add_argument("--no-feature-provenance", action="store_true")
    parser.add_argument("--merge-timestamp")
    parser.add_argument("--translation-evidence-manifest", type=Path)
    parser.add_argument("--allow-imported-translations", action="store_true")
    parser.add_argument(
        "--baktfold-invalid-ec-policy",
        choices=("reject", "skip"),
        default="reject",
        help="reject non-grammar Baktfold EC tokens or record them as non-promoted",
    )
    parser.add_argument(
        "--context-report",
        type=Path,
        help="write higher-order eggNOG/InterProScan context as a deterministic JSON sidecar",
    )
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    for label, path in (
        ("Bakta", args.bakta),
        ("FAA", args.faa),
        ("Baktfold", args.baktfold),
        ("KofamScan", args.kofamscan),
        ("eggNOG", args.eggnog),
        ("InterProScan", args.interproscan),
    ):
        if path is not None and not path.is_file():
            parser.error(f"{label} input not found: {path}")
    if args.interproscan and not args.interproscan_version:
        parser.error(
            "--interproscan-version is required when --interproscan is provided"
        )
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
            eggnog_schema=args.eggnog_schema,
            interproscan_path=args.interproscan,
            interproscan_version=args.interproscan_version,
            interproscan_tsv_layout=args.interproscan_tsv_layout,
            interproscan_member_dbs=args.interproscan_member_dbs,
            min_eggnog_confidence=args.min_eggnog_confidence,
            clean_gene_suffix=args.clean_gene_suffix,
            gene_conflict_policy=args.gene_conflict_policy,
            manifest_path=args.manifest,
            add_comment_note=not args.no_comment_note,
            add_feature_provenance=not args.no_feature_provenance,
            merge_timestamp=args.merge_timestamp,
            context_report_path=args.context_report,
            translation_evidence_manifest=args.translation_evidence_manifest,
            allow_imported_translations=args.allow_imported_translations,
            baktfold_invalid_ec_policy=args.baktfold_invalid_ec_policy,
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
