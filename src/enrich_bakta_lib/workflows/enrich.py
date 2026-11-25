#!/usr/bin/env python3
"""Canonical one-pass Baktfold, KofamScan, eggNOG, and InterProScan Bakta enrichment CLI."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from enrich_bakta_lib.core.merge_engine import (
    TOOL_VERSION,
    MergeError,
)
from enrich_bakta_lib.sources.baktfold import BaktfoldAdapter
from enrich_bakta_lib.sources.base import SourceAdapter
from enrich_bakta_lib.sources.eggnog import EggnogAdapter
from enrich_bakta_lib.sources.interproscan import (
    DEFAULT_MEMBER_DBS,
    DEFAULT_TSV_LAYOUT,
    TSV_LAYOUT_PROFILES,
    InterProScanAdapter,
)
from enrich_bakta_lib.sources.kofam import KofamAdapter
from enrich_bakta_lib.workflows.pipeline import (
    capture_merge_context,
    execute_merge_pipeline,
)


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
    ctx = capture_merge_context(
        bakta_path=bakta_path,
        faa_path=faa_path,
        translation_evidence_manifest=translation_evidence_manifest,
        allow_imported_translations=allow_imported_translations,
        add_comment_note=add_comment_note,
        add_feature_provenance=add_feature_provenance,
        merge_timestamp=merge_timestamp,
        clean_gene_suffix=clean_gene_suffix,
        require_manifest_message=(
            "restored input requires --translation-evidence-manifest for functional validation"
        ),
    )
    adapters: list[SourceAdapter] = []
    if baktfold_path:
        adapters.append(
            BaktfoldAdapter(
                baktfold_path=baktfold_path,
                invalid_ec_policy=baktfold_invalid_ec_policy,
            )
        )
    if kofamscan_path:
        adapters.append(
            KofamAdapter(
                kofam_path=kofamscan_path,
                kofamscan_version=kofamscan_version,
            )
        )
    if eggnog_path:
        adapters.append(
            EggnogAdapter(
                eggnog_path=eggnog_path,
                eggnog_version=eggnog_version,
                eggnog_schema=eggnog_schema,
                min_confidence=min_eggnog_confidence,
            )
        )
    if interproscan_path:
        assert interproscan_version is not None
        adapters.append(
            InterProScanAdapter(
                interproscan_path=interproscan_path,
                version=interproscan_version,
                layout=interproscan_tsv_layout or DEFAULT_TSV_LAYOUT,
                member_dbs=(
                    interproscan_member_dbs
                    if interproscan_member_dbs is not None
                    else DEFAULT_MEMBER_DBS
                ),
            )
        )
    final, _plans = execute_merge_pipeline(
        ctx,
        adapters,
        operation="unified-enrichment",
        output_path=output_path,
        manifest_path=manifest_path,
        context_report_path=context_report_path,
        gene_conflict_policy=gene_conflict_policy,
        reconcile=True,
        include_clean_gene_suffix_policy=True,
    )
    return final


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
